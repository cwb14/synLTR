"""The forward model of `ltrtrace sim nucleus`: LTR-retrotransposon evolution through progenitor splits and
allopolyploid mergers, with exact truth.

Model (time = years before present; older = larger):
- Lineages form a tree with split times; leaves are subgenomes. Mergers join the nuclei holding the listed leaves
  at a given time (sequential mergers give nested polyploids).
- Each TE family starts with 1-3 source elements (masters) in the root lineage; masters mutate and branch into new
  masters. At a split every master and every surviving copy is inherited by each daughter at the same coordinates
  (homoeologous chromosomes are collinear) and evolves independently afterwards.
- Insertions are Poisson in time per family and lineage (lineage-specific activity, optional bursts). The source
  is a master of the inserting lineage or, with probability p_copy_source, a copy of the family younger than
  active_age anywhere in the same nucleus. Both LTRs of a new copy equal the source LTR (plus optional divergence at
  birth). After a merger a copy lands in subgenome j with probability proportional to target_weight_j x (1 if j is
  the source's subgenome else phi).
- Copies are lost (deletion, solo LTR) at a constant hazard (half_life); losses are independent after a split.
- Homoeologous exchanges after the youngest merger duplicate (non-reciprocal) or swap (reciprocal) a segment
  between homoeologous chromosomes of two subgenomes; copies inside are duplicated or moved.
- Substitutions: every object carries its own Poisson substitutions (K2P: transition share kappa/(kappa+2)) coded
  as XOR masks on 2-bit bases (A0 C1 G2 T3; transition ^2, transversions ^1 and ^3), so they commute; a
  descendant's state is its ancestor's state at its birth XOR its own events. Optional LTR gene conversion copies a
  5'-LTR tract onto the 3' LTR.
Every random draw comes from numpy's default_rng seeded with the seed (the event stream) or with [seed, object,
tag] (an object's substitutions); the order of draws is part of the model.
"""
from __future__ import annotations

import copy
import heapq
import math
import sys
from collections import OrderedDict, defaultdict, deque
from dataclasses import dataclass, field

import numpy as np

from ltrtrace.base.constants import (
    SIM_NUCLEUS,
    SIM_NUCLEUS_CACHE,
    SIM_NUCLEUS_CLADES,
    SIM_NUCLEUS_INT_LEN,
    SIM_NUCLEUS_LTR_LEN,
    SIM_NUCLEUS_SUPERFAMILIES,
)
from ltrtrace.seq.divergence import k2p_shares

MASKS = np.array([2, 1, 3], np.uint8)  # transition, transversion 1, transversion 2 (XOR on A0 C1 G2 T3)
ROOT, INSERT, DUP, BRANCH = 0, 1, 2, 3  # how an object came to be: a root master, an insertion, an inheritance, a branch


def merge_defaults(s, d=SIM_NUCLEUS):
    out = copy.deepcopy(d)
    for k, v in s.items():
        out[k] = merge_defaults(v, out[k]) if isinstance(v, dict) and isinstance(out.get(k), dict) and k != "tree" else v
    return out


# the scenario values check_scenario reads, by the domain the simulation needs
_COUNTS = ("chromosomes.n", "chromosomes.length", "families.n", "conversion_tract")  # integers >= 1
_COUNTS0 = ("exchanges.nonreciprocal", "exchanges.reciprocal")  # integers >= 0
_POSITIVE = ("root_age", "half_life", "kappa", "families.ltr_len_median", "families.int_len_median",
             "exchanges.size_min", "exchanges.size_max")  # numbers > 0
_NONNEGATIVE = ("mu", "families.branch_rate", "families.ltr_len_sigma", "families.int_len_sigma",
                "families.size_sigma", "activity.rate", "activity.lineage_sigma", "activity.burst_width",
                "activity.burst_gain", "phi", "active_age", "birth_divergence", "gene_conversion",
                "family_rate_sigma")  # numbers >= 0
_SHARES = ("families.gc", "activity.p_off", "activity.p_burst", "p_copy_source", "singleton_share")  # in [0, 1]
# the keys of what a scenario writes whole (the defaults fill in none of them), each with what it is: a split node of
# the tree, a merger, an explicit exchange
_NODE_KEYS = {"split": "when the lineage splits, years ago",
              "children": "the lineages it splits into: split nodes or leaf names"}
_MERGER_KEYS = {"time": "when the nuclei merge, years ago", "lineages": "the names of the leaves whose nuclei merge"}
_EXCHANGE_KEYS = {"type": "nonreciprocal (the donor's segment replaces the recipient's) or reciprocal (the two swap)",
                  "time": "when the segment moves, years ago", "recipient": "the receiving leaf's name",
                  "donor": "the giving leaf's name", "chrom": "its chromosome, 0 to chromosomes.n - 1",
                  "start": "the segment's start, bp", "end": "the segment's end, bp"}


def check_scenario(sc: dict) -> None:
    """ValueError naming the first value of scenario `sc` (defaults merged) that the simulation cannot use: a count
    or length that is not an integer >= 1 (an exchange count >= 0), a time, length or size that is not a number > 0,
    a rate, spread or weight below 0, a share outside [0, 1], masters that are not [low, high] with 1 <= low <= high,
    or exchanges.size_min above size_max. The tree, and the keys and names of mergers and explicit exchanges, are
    checked as the simulation builds them (Sim)."""
    def value(key):
        v = sc
        try:
            for k in key.split("."):
                v = v[k]
        except (KeyError, TypeError, IndexError):
            raise ValueError(f"{key} is missing or not in the defaults' form (ltrtrace sim nucleus --help)") from None
        return v

    def number(key, ok, want):
        v = value(key)
        if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v) or not ok(v):
            raise ValueError(f"{key} must be {want}: {v!r}")

    for key in _COUNTS:
        number(key, lambda v: v == int(v) and v >= 1, "an integer >= 1")
    for key in _COUNTS0:
        number(key, lambda v: v == int(v) and v >= 0, "an integer >= 0")
    for key in _POSITIVE:
        number(key, lambda v: v > 0, "> 0")
    for key in _NONNEGATIVE:
        number(key, lambda v: v >= 0, ">= 0")
    for key in _SHARES:
        number(key, lambda v: 0 <= v <= 1, "in [0, 1]")
    lo_hi = value("families.masters")
    if not (isinstance(lo_hi, list) and len(lo_hi) == 2 and all(isinstance(x, int) and not isinstance(x, bool)
                                                                 for x in lo_hi) and 1 <= lo_hi[0] <= lo_hi[1]):
        raise ValueError(f"families.masters must be [low, high] with 1 <= low <= high: {lo_hi!r}")
    if value("exchanges.size_min") > value("exchanges.size_max"):
        raise ValueError(f"exchanges.size_min ({value('exchanges.size_min')}) is above exchanges.size_max "
                         f"({value('exchanges.size_max')})")


def _has_keys(item, what: str, keys: dict[str, str]) -> None:
    """ValueError naming the first of `keys` that `item` (`what`: a split node, a merger, an exchange) lacks, with
    what that key is; TypeError when `item` is not a JSON object."""
    if not isinstance(item, dict):
        raise TypeError(f"{what} must be a JSON object with the keys {', '.join(keys)}: {item!r}")
    for key, meaning in keys.items():
        if key not in item:
            raise ValueError(f"{what} lacks the key {key!r} ({meaning})")


@dataclass
class Lineage:
    name: str
    parent: int
    t_from: float
    t_to: float
    children: list = field(default_factory=list)


def build_lineages(tree, root_age):
    out = []

    def rec(node, parent, t_from):
        if isinstance(node, str):
            out.append(Lineage(node, parent, t_from, 0.0))
            return len(out) - 1
        _has_keys(node, "a split node of tree", _NODE_KEYS)
        idx = len(out)
        t = float(node["split"])
        if not 0 < t < t_from:
            raise ValueError(f"split time {t} must lie between 0 and its parent's start {t_from}")
        out.append(Lineage(f"anc{idx}", parent, t_from, t))
        for ch in node["children"]:
            out[idx].children.append(rec(ch, idx, t))
        return idx
    rec(tree, -1, float(root_age))
    return out


class Sim:
    def __init__(self, scenario, seed):
        check_scenario(scenario)
        self.sc, self.seed = scenario, int(seed)
        self.rng = np.random.default_rng(self.seed)
        self.lin = build_lineages(scenario["tree"], scenario["root_age"])
        self.leaves = [i for i, l in enumerate(self.lin) if not l.children]
        self.leaf_of = {self.lin[i].name: i for i in self.leaves}
        self._check_names()
        self.n_chr, self.chr_len = int(scenario["chromosomes"]["n"]), int(scenario["chromosomes"]["length"])
        self.lam = math.log(2) / scenario["half_life"]
        self._families()
        self._activity()
        self.kind, self.fam, self.tb, self.tend, self.par, self.mode = [], [], [], [], [], []
        self.lineage, self.chrom, self.pos, self.strand, self.origin = [], [], [], [], []
        self.tins, self.exid, self.root = [], [], []
        self.masters = defaultdict(list)        # (family, lineage) -> master ids
        self.recent = defaultdict(deque)        # (family, lineage) -> copy ids in insertion order
        self.lin_copies = defaultdict(list)     # lineage -> copy ids
        self.chr_copies = defaultdict(list)     # (lineage, chromosome) -> copy ids
        self.uf = {l: l for l in self.leaves}   # nucleus union-find over leaves
        self.merge_time = {l: 0.0 for l in self.leaves}
        self.exchanges = []

    def _check_names(self):
        """Mergers and explicit exchanges hold their keys and name leaves of the tree; anything else is an error
        before simulating."""
        for mg in self.sc["mergers"]:
            _has_keys(mg, "a merger", _MERGER_KEYS)
        for e in self.sc["exchange_list"]:
            _has_keys(e, "an exchange of exchange_list", _EXCHANGE_KEYS)
        asked = [n for mg in self.sc["mergers"] for n in mg["lineages"]]
        asked += [e[k] for e in self.sc["exchange_list"] for k in ("recipient", "donor")]
        bad = sorted({n for n in asked if n not in self.leaf_of})
        if bad:
            raise ValueError(f"not leaves of the tree {sorted(self.leaf_of)}: {', '.join(map(str, bad))}")

    # ---- setup ----
    def _families(self):
        f, r = self.sc["families"], self.rng
        n = int(f["n"])
        self.ltr_len = np.clip(np.round(f["ltr_len_median"] * np.exp(r.normal(0, f["ltr_len_sigma"], n))),
                          *SIM_NUCLEUS_LTR_LEN).astype(int)
        self.int_len = np.clip(np.round(f["int_len_median"] * np.exp(r.normal(0, f["int_len_sigma"], n))),
                          *SIM_NUCLEUS_INT_LEN).astype(int)
        s = self.sc["family_rate_sigma"]
        self.mu_f = self.sc["mu"] * (np.exp(r.normal(-s * s / 2, s, n)) if s > 0 else np.ones(n))
        gc = f["gc"]
        p = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
        self.root_ltr = [r.choice(4, size=L, p=p).astype(np.uint8) for L in self.ltr_len]
        self.root_int = [r.choice(4, size=L, p=p).astype(np.uint8) for L in self.int_len]
        self.superfamily = [str(x) for x in r.choice(list(SIM_NUCLEUS_SUPERFAMILIES[0]), size=n,
                                                      p=list(SIM_NUCLEUS_SUPERFAMILIES[1]))]
        self.clade = [str(r.choice(SIM_NUCLEUS_CLADES[x])) if x in SIM_NUCLEUS_CLADES else "unknown" for x in self.superfamily]
        self.n_masters = r.integers(f["masters"][0], f["masters"][1] + 1, size=n)
        w = np.exp(r.normal(0, f["size_sigma"], n))
        self.fam_w = w / w.sum()

    def _activity(self):
        a, r = self.sc["activity"], self.rng
        nF = len(self.fam_w)
        self.act = np.zeros((nF, len(self.lin)))
        for i, l in enumerate(self.lin):                     # parents are built before children
            if l.parent < 0:
                self.act[:, i] = self.fam_w
            else:
                m = self.act[:, l.parent] * np.exp(r.normal(0, a["lineage_sigma"], nF))
                m[r.random(nF) < a["p_off"]] = 0.0
                self.act[:, i] = m
        self.act *= a["rate"] / np.maximum(self.act.sum(0, keepdims=True), 1e-300)

    def _elem_len(self, f):
        return int(2 * self.ltr_len[f] + self.int_len[f])

    # ---- object table ----
    def _new(self, kind, fam, tb, par, mode, lineage, chrom=-1, pos=-1, strand=1, origin="", tins=0.0, exid=-1,
             root=-1):
        oid = len(self.kind)
        if kind == 1:
            delay = self.rng.exponential(1.0 / self.lam)
            tend = tb - delay if tb - delay > 0 else 0.0
        else:
            tend = 0.0
        for col, v in ((self.kind, kind), (self.fam, fam), (self.tb, tb), (self.tend, tend), (self.par, par),
                       (self.mode, mode), (self.lineage, lineage), (self.chrom, chrom), (self.pos, pos),
                       (self.strand, strand), (self.origin, origin), (self.tins, tins), (self.exid, exid),
                       (self.root, root if root >= 0 else oid)):
            col.append(v)
        if kind == 1:
            self.lin_copies[lineage].append(oid)
            self.chr_copies[(lineage, chrom)].append(oid)
        return oid

    # ---- nuclei ----
    def _find(self, x):
        while self.uf[x] != x:
            self.uf[x] = self.uf[self.uf[x]]
            x = self.uf[x]
        return x

    def nucleus(self, l):
        if l not in self.uf:
            return [l]
        r = self._find(l)
        return [x for x in self.leaves if self._find(x) == r]

    # ---- events ----
    def _push(self, t, kind, *data):
        heapq.heappush(self.heap, (-float(t), self._seq, kind, data))
        self._seq += 1

    def _schedule_branch(self, m, t_now):
        b = self.sc["families"]["branch_rate"]
        if b > 0:
            t = t_now - self.rng.exponential(1.0 / b)
            if t > self.lin[self.lineage[m]].t_to:
                self._push(t, "branch", m)

    def _exchange_specs(self):
        if self.sc["exchange_list"]:
            out = []
            for e in self.sc["exchange_list"]:
                out.append(dict(e, recipient=self.leaf_of[e["recipient"]], donor=self.leaf_of[e["donor"]]))
            return out
        if not self.sc["mergers"]:
            return []
        t_young = min(m["time"] for m in self.sc["mergers"])
        uf = {l: l for l in self.leaves}

        def find(x):
            while uf[x] != x:
                x = uf[x]
            return x
        for m in self.sc["mergers"]:
            ls = [self.leaf_of[n] for n in m["lineages"]]
            for x in ls[1:]:
                uf[find(x)] = find(ls[0])
        groups = defaultdict(list)
        for l in self.leaves:
            groups[find(l)].append(l)
        nucs = [g for g in groups.values() if len(g) > 1]
        e, r, out = self.sc["exchanges"], self.rng, []
        for typ in ["nonreciprocal"] * int(e["nonreciprocal"]) + ["reciprocal"] * int(e["reciprocal"]):
            g = nucs[r.integers(len(nucs))]
            rec, don = r.choice(g, 2, replace=False)
            size = float(np.exp(r.uniform(math.log(e["size_min"]), math.log(e["size_max"]))))
            start = float(r.uniform(0, self.chr_len - size))
            out.append({"type": typ, "time": float(r.uniform(0, t_young)), "recipient": int(rec), "donor": int(don),
                            "chrom": int(r.integers(self.n_chr)), "start": start, "end": start + size})
        return out

    def run(self):
        self.heap, self._seq = [], 0
        root = next(i for i, l in enumerate(self.lin) if l.parent < 0)
        for f in range(len(self.fam_w)):
            for _ in range(self.n_masters[f]):
                m = self._new(0, f, self.lin[root].t_from, -1, ROOT, root, origin=self.lin[root].name)
                self.masters[(f, root)].append(m)
                self._schedule_branch(m, self.lin[root].t_from)
        for i, l in enumerate(self.lin):
            if l.children:
                self._push(l.t_to, "split", i)
        for j, m in enumerate(self.sc["mergers"]):
            self._push(m["time"], "merge", j)
        self.ex_specs = self._exchange_specs()
        for k, e in enumerate(self.ex_specs):
            self._push(e["time"], "exchange", k)
        a, r = self.sc["activity"], self.rng
        for i, l in enumerate(self.lin):
            for f in np.flatnonzero(self.act[:, i] > 0):
                rate, span = self.act[f, i], l.t_from - l.t_to
                for t in r.uniform(l.t_to, l.t_from, r.poisson(rate * span)):
                    self._push(t, "insert", int(f), i)
                if r.random() < a["p_burst"]:
                    c, w = r.uniform(l.t_to, l.t_from), a["burst_width"]
                    lo, hi = max(c - w, l.t_to), min(c + w, l.t_from)
                    for t in r.uniform(lo, hi, r.poisson(a["burst_gain"] * rate * (hi - lo))):
                        self._push(t, "insert", int(f), i)
        while self.heap:
            negt, _, kind, data = heapq.heappop(self.heap)
            getattr(self, "_ev_" + kind)(-negt, *data)
        return self

    def _ev_branch(self, t, m):
        new = self._new(0, self.fam[m], t, m, BRANCH, self.lineage[m], origin=self.origin[m])
        self.masters[(self.fam[m], self.lineage[m])].append(new)
        self._schedule_branch(new, t)
        self._schedule_branch(m, t)

    def _ev_split(self, t, i):
        for (f, l), ms in list(self.masters.items()):
            if l != i:
                continue
            for m in ms:
                self.tend[m] = t
                for c in self.lin[i].children:
                    d = self._new(0, f, t, m, DUP, c, origin=self.lin[c].name)
                    self.masters[(f, c)].append(d)
                    self._schedule_branch(d, t)
        for x in self.lin_copies[i]:
            if not self.tend[x] < t:
                continue
            self.tend[x] = t
            for c in self.lin[i].children:
                y = self._new(1, self.fam[x], t, x, DUP, c, self.chrom[x], self.pos[x], self.strand[x],
                              self.origin[x], self.tins[x], self.exid[x], self.root[x])
                if self.tins[x] - t < self.sc["active_age"]:
                    self.recent[(self.fam[x], c)].append(y)

    def _ev_merge(self, t, j):
        ls = [self.leaf_of[n] for n in self.sc["mergers"][j]["lineages"]]
        roots = {self._find(x) for x in ls}
        members = [x for x in self.leaves if self._find(x) in roots]
        for x in members:
            if self.merge_time[x] == 0.0:
                self.merge_time[x] = t
        r0 = self._find(ls[0])
        for x in ls[1:]:
            self.uf[self._find(x)] = r0

    def _ev_insert(self, t, f, l):
        members = self.nucleus(l)
        src = None
        if self.rng.random() < self.sc["p_copy_source"]:
            pool = []
            for x in members:
                dq = self.recent[(f, x)]
                while dq and self.tins[dq[0]] - t > self.sc["active_age"]:
                    dq.popleft()
                pool += [c for c in dq if self.tend[c] < t]
            if pool:
                src = pool[self.rng.integers(len(pool))]
        if src is None:
            ms = self.masters[(f, l)]
            if not ms:
                return
            src = ms[self.rng.integers(len(ms))]
        src_lin = l if self.kind[src] == 0 else self.lineage[src]
        if self.kind[src] == 0:
            origin = self.lin[l].name
        else:                                   # a source inherited from the ancestor transposes from its own leaf
            origin = self.origin[src] if self.origin[src] in self.leaf_of else self.lin[self.lineage[src]].name
        if len(members) > 1:
            tw = self.sc["target_weight"]
            w = np.array([tw.get(self.lin[x].name, 1.0) * (1.0 if x == src_lin else self.sc["phi"]) for x in members])
            target = members[self.rng.choice(len(members), p=w / w.sum())] if w.sum() > 0 else src_lin
        else:
            target = l
        L = self._elem_len(f)
        c = self._new(1, f, t, src, INSERT, target, int(self.rng.integers(self.n_chr)),
                      int(self.rng.integers(0, self.chr_len - L)), int(self.rng.random() < 0.5), origin, t)
        self.recent[(f, target)].append(c)

    def _ev_exchange(self, t, k):
        e = self.ex_specs[k]
        R, D, ch = e["recipient"], e["donor"], e["chrom"]

        def seg(l):
            return [x for x in self.chr_copies[(l, ch)] if self.tend[x] < t < self.tb[x]
                    and e["start"] <= self.pos[x] < e["end"]]
        sR, sD = seg(R), seg(D)
        moved = 0

        def dup(x, dest):
            self._new(1, self.fam[x], t, x, DUP, dest, ch, self.pos[x], self.strand[x], self.origin[x],
                      self.tins[x], k, self.root[x])
        if e["type"] == "nonreciprocal":
            for x in sR:
                self.tend[x] = t
            for x in sD:
                dup(x, R)
                moved += 1
        else:
            for x in sR:
                dup(x, D)
                self.tend[x] = t
            for x in sD:
                dup(x, R)
                self.tend[x] = t
            moved = len(sR) + len(sD)
        self.exchanges.append({"type": e["type"], "time": t, "recipient": self.lin[R].name, "donor": self.lin[D].name,
                                   "chrom": ch, "start": int(e["start"]), "end": int(e["end"]), "n_moved": moved, "index": k})

    def source_root(self, c):
        """The root of the copy that copy `c` was transposed from (the `root` of its original insertion), or -1
        when the source was a master."""
        src = self.par[self.root[c]]
        return self.root[src] if self.kind[src] == 1 else -1

    def survivors(self):
        leaves = set(self.leaves)
        return [c for c in range(len(self.kind)) if self.kind[c] == 1 and self.tend[c] == 0.0
                and self.lineage[c] in leaves]


# ---------------------------------------------------------------------------------------------- sequences
class Seqs:
    """Lazy region states. state(oid, t) = ancestor state at the object's birth (transformed by the birth mode)
    XOR the object's own substitutions older than t, with LTR conversions applied in time order."""

    def __init__(self, sim, cache_size=SIM_NUCLEUS_CACHE):
        self.s, self.ev, self.cache, self.cache_size = sim, {}, OrderedDict(), cache_size
        k = sim.sc["kappa"]
        self.pmask = np.array([k, 1.0, 1.0]) / (k + 2)
        sys.setrecursionlimit(max(sys.getrecursionlimit(), 100000))

    def _rng(self, oid, tag):
        return np.random.default_rng([self.s.seed, oid, tag])

    def events(self, oid):
        if oid in self.ev:
            return self.ev[oid]
        s = self.s
        f, t_hi, t_lo = s.fam[oid], s.tb[oid], s.tend[oid]
        mu, r = s.mu_f[f], self._rng(oid, 1)
        regs = ("ltr", "int") if s.kind[oid] == 0 else ("ltr5", "int", "ltr3")
        out = {}
        for reg in regs:
            L = int(s.int_len[f] if reg == "int" else s.ltr_len[f])
            n = r.poisson(mu * (t_hi - t_lo) * L)
            out[reg] = (r.uniform(t_lo, t_hi, n), r.integers(0, L, n), MASKS[r.choice(3, n, p=self.pmask)])
        conv, tract = [], 0
        g = s.sc["gene_conversion"]
        if s.kind[oid] == 1 and g > 0:
            L = int(s.ltr_len[f])
            tract = min(int(s.sc["conversion_tract"]), L)
            n = r.poisson(g * (t_hi - t_lo))
            conv = sorted(zip(r.uniform(t_lo, t_hi, n), r.integers(0, L - tract + 1, n)), reverse=True)
        self.ev[oid] = (out, conv, tract)
        return self.ev[oid]

    def state(self, oid, t):
        key = (oid, float(t))
        if key in self.cache:
            self.cache.move_to_end(key)
            return {k: v.copy() for k, v in self.cache[key].items()}
        s = self.s
        f, p = s.fam[oid], s.par[oid]
        if p < 0:
            st = {"ltr": s.root_ltr[f].copy(), "int": s.root_int[f].copy()}
        else:
            ps = self.state(p, s.tb[oid])
            if s.mode[oid] == INSERT:
                ltr = ps["ltr"] if s.kind[p] == 0 else ps["ltr5"]
                st = {"ltr5": ltr.copy(), "int": ps["int"], "ltr3": ltr.copy()}
                d0 = s.sc["birth_divergence"]
                if d0 > 0:
                    r = self._rng(oid, 2)
                    n = r.poisson(d0)
                    if n:
                        np.bitwise_xor.at(st["ltr3"], r.integers(0, ltr.size, n), MASKS[r.choice(3, n, p=self.pmask)])
            else:
                st = ps
        ev, conv, tract = self.events(oid)
        cuts = [c for c in conv if c[0] > t]
        bounds = [s.tb[oid]] + [c[0] for c in cuts] + [t]
        for i in range(len(bounds) - 1):
            hi, lo = bounds[i], bounds[i + 1]
            for reg, (tt, pp, mm) in ev.items():
                sel = (tt <= hi) & (tt > lo)
                if sel.any():
                    np.bitwise_xor.at(st[reg], pp[sel], mm[sel])
            if i < len(cuts):
                a = int(cuts[i][1])
                st["ltr3"][a:a + tract] = st["ltr5"][a:a + tract]
        if t > 0:
            self.cache[key] = {k: v.copy() for k, v in st.items()}
            if len(self.cache) > self.cache_size:
                self.cache.popitem(last=False)
        return st


def k2p(a, b):
    """Kimura 2-parameter distance between two equal-length code arrays (seq.divergence.k2p_shares): (d, se, n_ts,
    n_tv); NaN when saturated. Codes 0-3 are A, C, G, T, so a transition (A-G, C-T) has XOR 2."""
    x = a ^ b
    n = a.size
    ts, tv = int(np.count_nonzero(x == 2)), int(np.count_nonzero((x == 1) | (x == 3)))
    d, se = k2p_shares(ts / n, tv / n, n)
    return d, se, ts, tv
