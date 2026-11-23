#!/usr/bin/env python3
"""simallo.py - forward simulation of LTR-retrotransposon evolution through progenitor splits and allopolyploid
mergers, written as an LTRquest-style annotation (LTR-RT FASTA, GFF3, K2P table, .fai) with exact truth.

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
Writes OUTDIR/{ltr.fa, annotation.gff3, k2p.tsv, genome.fai} and TRUTHDIR/{copies.tsv, exchanges.tsv,
chromosomes.tsv, orthologs.tsv, paralogs.tsv, history.json}.
usage: simallo.py SCENARIO.json OUTDIR TRUTHDIR [--seed N]
"""
from __future__ import annotations

import argparse
import copy
import heapq
import json
import math
import sys
from collections import OrderedDict, defaultdict, deque
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

DEFAULTS = {
    "mu": 1.3e-8,               # substitutions/site/year
    "kappa": 2.0,               # transition/transversion rate ratio (K2P)
    "root_age": 1.0e7,          # years: start of the root lineage (copies older than ~5 half-lives are rare)
    "tree": {"split": 2.0e6, "children": ["A", "B"]},
    "mergers": [{"time": 1.0e5, "lineages": ["A", "B"]}],
    "chromosomes": {"n": 8, "length": 4.0e7},
    "families": {"n": 400, "masters": [1, 3], "branch_rate": 3e-7, "ltr_len_median": 600, "ltr_len_sigma": 0.8,
                 "int_len_median": 4500, "int_len_sigma": 0.4, "size_sigma": 2.0, "gc": 0.40},
    "activity": {"rate": 0.004, "lineage_sigma": 1.0, "p_off": 0.2, "p_burst": 0.3, "burst_width": 1.0e5,
                 "burst_gain": 10.0},
    "target_weight": {},        # leaf -> relative insertion target weight (default 1)
    "phi": 1.0,                 # cross-subgenome insertion preference after a merger (1 = by size only)
    "p_copy_source": 0.3,       # share of insertions sourced from a recent copy rather than a master
    "active_age": 1.0e5,        # copies younger than this can be sources
    "half_life": 1.5e6,         # years: loss of intact copies (deletion, solo LTR)
    "birth_divergence": 0.0,    # mean substitutions between the two LTRs at insertion
    "gene_conversion": 0.0,     # LTR-LTR conversion events per copy per year
    "conversion_tract": 200,    # bp per conversion event
    "family_rate_sigma": 0.0,   # log-sd of family substitution-rate multipliers (mean 1)
    "exchanges": {"nonreciprocal": 4, "reciprocal": 2, "size_min": 2.0e5, "size_max": 5.0e6},
    "exchange_list": [],        # explicit exchanges (override the random ones)
    "singleton_share": 0.14,    # copies written as one-copy families (LTRquest: 0.07-0.34 of copies, median 0.14
                                # over 6 libraries; older copies more often), true family kept in copies.tsv
}
CLADES = {"Gypsy": ["Tekay", "CRM", "Athila", "Retand", "Reina", "Ogre", "TatII"],
          "Copia": ["Ale", "SIRE", "Ivana", "Tork", "Angela", "Bianca", "TAR", "Ikeros"]}
MASKS = np.array([2, 1, 3], np.uint8)          # transition, transversion 1, transversion 2 (XOR on A0 C1 G2 T3)
ROOT, INSERT, DUP, BRANCH = 0, 1, 2, 3


def merge_defaults(s, d=DEFAULTS):
    out = copy.deepcopy(d)
    for k, v in s.items():
        out[k] = merge_defaults(v, out[k]) if isinstance(v, dict) and isinstance(out.get(k), dict) and k != "tree" else v
    return out


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
        self.sc, self.seed = scenario, int(seed)
        self.rng = np.random.default_rng(self.seed)
        self.lin = build_lineages(scenario["tree"], scenario["root_age"])
        self.leaves = [i for i, l in enumerate(self.lin) if not l.children]
        self.leaf_of = {self.lin[i].name: i for i in self.leaves}
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

    # ---- setup ----
    def _families(self):
        f, r = self.sc["families"], self.rng
        n = int(f["n"])
        self.ltr_len = np.clip(np.round(f["ltr_len_median"] * np.exp(r.normal(0, f["ltr_len_sigma"], n))), 100, 4000).astype(int)
        self.int_len = np.clip(np.round(f["int_len_median"] * np.exp(r.normal(0, f["int_len_sigma"], n))), 1000, 12000).astype(int)
        s = self.sc["family_rate_sigma"]
        self.mu_f = self.sc["mu"] * (np.exp(r.normal(-s * s / 2, s, n)) if s > 0 else np.ones(n))
        gc = f["gc"]
        p = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
        self.root_ltr = [r.choice(4, size=L, p=p).astype(np.uint8) for L in self.ltr_len]
        self.root_int = [r.choice(4, size=L, p=p).astype(np.uint8) for L in self.int_len]
        self.superfamily = [str(x) for x in r.choice(["Gypsy", "Copia", "unknown"], size=n, p=[0.5, 0.4, 0.1])]
        self.clade = [str(r.choice(CLADES[x])) if x in CLADES else "unknown" for x in self.superfamily]
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
            out.append(dict(type=typ, time=float(r.uniform(0, t_young)), recipient=int(rec), donor=int(don),
                            chrom=int(r.integers(self.n_chr)), start=start, end=start + size))
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
        self.exchanges.append(dict(type=e["type"], time=t, recipient=self.lin[R].name, donor=self.lin[D].name,
                                   chrom=ch, start=int(e["start"]), end=int(e["end"]), n_moved=moved, index=k))

    def survivors(self):
        leaves = set(self.leaves)
        return [c for c in range(len(self.kind)) if self.kind[c] == 1 and self.tend[c] == 0.0
                and self.lineage[c] in leaves]


# ---------------------------------------------------------------------------------------------- sequences
class Seqs:
    """Lazy region states. state(oid, t) = ancestor state at the object's birth (transformed by the birth mode)
    XOR the object's own substitutions older than t, with LTR conversions applied in time order."""

    def __init__(self, sim, cache_size=4096):
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
    """Kimura 2-parameter distance between two equal-length code arrays: (d, se, n_ts, n_tv); NaN when saturated."""
    x = a ^ b
    n = a.size
    ts, tv = int(np.count_nonzero(x == 2)), int(np.count_nonzero((x == 1) | (x == 3)))
    P, Q = ts / n, tv / n
    w1, w2 = 1 - 2 * P - Q, 1 - 2 * Q
    if w1 <= 0 or w2 <= 0:
        return float("nan"), float("nan"), ts, tv
    d = -0.5 * math.log(w1) - 0.25 * math.log(w2)
    c1, c2 = 1 / w1, 0.5 * (1 / w1 + 1 / w2)
    var = (c1 * c1 * P + c2 * c2 * Q - (c1 * P + c2 * Q) ** 2) / n
    return d, math.sqrt(max(var, 0.0)), ts, tv


# ---------------------------------------------------------------------------------------------- outputs
BASES = np.frombuffer(b"ACGT", np.uint8)
COMP = np.array([3, 2, 1, 0], np.uint8)
TSV_COLS = ["seq_id", "seq_len", "status", "ltr5_start", "ltr5_end", "ltr3_start", "ltr3_end", "ltr5_len", "ltr3_len",
            "flank5_len", "flank3_len", "aln_len", "n_sites", "n_ts", "n_tv", "n_gapcols", "identity", "p_dist", "k2p",
            "k2p_se", "bitscore", "flank_margin_bits", "cigar", "motif", "k2p_time", "orientation", "tsd",
            "tsd_offset", "tsd_input", "strand", "family", "domains", "nest_status"]


def _fmt(x):
    return "NA" if isinstance(x, float) and math.isnan(x) else (f"{x:.6g}" if isinstance(x, float) else str(x))


def _mrca_split(sim, a, b):
    def path(i):
        out = []
        while i >= 0:
            out.append(i)
            i = sim.lin[i].parent
        return out
    pa = path(a)
    for i in path(b):
        if i in pa:
            return sim.lin[i].t_to
    return float("nan")


def _unique_coordinates(sim, surv, tok):
    """Survivors with distinct output coordinates (an annotator never reports two elements at identical
    coordinates; a chance collision keeps the first) and the number dropped."""
    seen, keep = set(), []
    for c in surv:
        k = (tok[(sim.lineage[c], sim.chrom[c])], sim.pos[c], sim._elem_len(sim.fam[c]))
        if k not in seen:
            seen.add(k)
            keep.append(c)
    return keep, len(surv) - len(keep)


def input_family_labels(sim, surv):
    """Family labels as an annotator would give them: the true family, except a share of copies (more often the
    older ones: probability 2 * share * age rank / n) written as one-copy families."""
    n = len(surv)
    rank = np.argsort(np.argsort([sim.tins[c] for c in surv], kind="stable"), kind="stable")
    p = 2 * sim.sc["singleton_share"] * rank / max(n - 1, 1)
    single = np.random.default_rng([sim.seed, 7]).random(n) < p
    out, k = {}, 0
    for c, sg in zip(surv, single):
        if sg:
            k += 1
            out[c] = f"single{k:07d}"
        else:
            out[c] = f"fam{sim.fam[c]:04d}"
    return out


def write_outputs(sim, outdir, truthdir):
    outdir, truthdir = Path(outdir), Path(truthdir)
    outdir.mkdir(parents=True, exist_ok=True)
    truthdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng([sim.seed, 99])
    pairs = [(l, c) for l in sim.leaves for c in range(sim.n_chr)]
    tok = {p: f"chr{j + 1:03d}" for j, p in zip(rng.permutation(len(pairs)), pairs)}
    with open(outdir / "genome.fai", "w") as f:
        for p in sorted(pairs, key=tok.get):
            f.write(f"{tok[p]}\t{sim.chr_len}\t0\t0\t0\n")
    with open(truthdir / "chromosomes.tsv", "w") as f:
        f.write("token\tsubgenome\thomoeolog\tlength\n")
        for p in sorted(pairs, key=tok.get):
            f.write(f"{tok[p]}\t{sim.lin[p[0]].name}\t{p[1]}\t{sim.chr_len}\n")
    q = Seqs(sim)
    surv = sorted(sim.survivors(), key=lambda c: (tok[(sim.lineage[c], sim.chrom[c])], sim.pos[c]))
    surv, n_dup = _unique_coordinates(sim, surv, tok)
    label = input_family_labels(sim, surv)
    fam_size = defaultdict(int)
    for c in surv:
        fam_size[label[c]] += 1
    cid = {}
    leafname = {i: sim.lin[i].name for i in sim.leaves}
    with open(outdir / "ltr.fa", "w") as ffa, open(outdir / "k2p.tsv", "w") as ftsv, \
            open(outdir / "annotation.gff3", "w") as fgff, open(truthdir / "copies.tsv", "w") as ftr:
        ftsv.write("#" + "\t".join(TSV_COLS) + "\n")
        fgff.write("##gff-version 3\n")
        for p in sorted(pairs, key=tok.get):
            fgff.write(f"##sequence-region\t{tok[p]}\t1\t{sim.chr_len}\n")
        ftr.write("copy_id\tchrom\tsubgenome\thomoeolog\tfamily\tsuperfamily\tclade\tt_insert\tt_birth\torigin\tmode\t"
                  "post_merger\texchange_id\troot\tn_sites\tn_ts\tn_tv\tk2p\n")
        for n, c in enumerate(surv, 1):
            f = sim.fam[c]
            st = q.state(c, 0.0)
            codes = np.concatenate([st["ltr5"], st["int"], st["ltr3"]])   # element orientation, as LTRquest writes
            L, lt = codes.size, int(sim.ltr_len[f])
            chrom = tok[(sim.lineage[c], sim.chrom[c])]
            s0, e0 = sim.pos[c] + 1, sim.pos[c] + L
            name = f"{chrom}:{s0}-{e0}"
            cid[c] = name
            cls = f"LTR/{sim.superfamily[f]}/{sim.clade[f]}"
            d, se, ts, tv = k2p(st["ltr5"], st["ltr3"])
            ffa.write(f">{name}#{cls}\n{BASES[codes].tobytes().decode()}\n")
            ident = 1 - (ts + tv) / lt
            row = dict(seq_id=f"{name}#{cls}", seq_len=L, status="pass", ltr5_start=1, ltr5_end=lt,
                       ltr3_start=L - lt + 1, ltr3_end=L, ltr5_len=lt, ltr3_len=lt, flank5_len=0, flank3_len=0,
                       aln_len=lt, n_sites=lt, n_ts=ts, n_tv=tv, n_gapcols=0, identity=ident, p_dist=(ts + tv) / lt,
                       k2p=d, k2p_se=se, bitscore="NA", flank_margin_bits="NA", cigar=".", motif="tg...ca",
                       k2p_time="NA" if math.isnan(d) else round(d / (2 * 3e-8)), orientation="+" if sim.strand[c] else "-",
                       tsd=".", tsd_offset="NA", tsd_input=".", strand="+" if sim.strand[c] else "-",
                       family=label[c], domains=".", nest_status=".")
            ftsv.write("\t".join(_fmt(row[k]) for k in TSV_COLS) + "\n")
            sd = "+" if sim.strand[c] else "-"
            fgff.write(f"{chrom}\tsimallo\tLTR_retrotransposon\t{s0}\t{e0}\t.\t{sd}\t.\tID=sim{n:07d};Name={name};"
                       f"classification={cls};order=LTR;superfamily={sim.superfamily[f]};clade={sim.clade[f]};"
                       f"family={label[c]};family_size={fam_size[label[c]]};depth=0;ltr5_len={lt};ltr3_len={lt};aln_len={lt};"
                       f"identity={ident:.6g};subs={ts + tv};ti={ts};tv={tv};K2P_d={_fmt(d)}\n")
            fgff.write(f"{chrom}\tsimallo\tlong_terminal_repeat\t{s0}\t{s0 + lt - 1}\t.\t{sd}\t.\tID=sim{n:07d}.lLTR;"
                       f"Parent=sim{n:07d}\n")
            fgff.write(f"{chrom}\tsimallo\tlong_terminal_repeat\t{e0 - lt + 1}\t{e0}\t.\t{sd}\t.\tID=sim{n:07d}.rLTR;"
                       f"Parent=sim{n:07d}\n")
            leaf = leafname[sim.lineage[c]]
            mode = "own" if sim.origin[c] == leaf else "foreign" if sim.origin[c] in leafname.values() else "ancestral"
            post = sim.tins[c] < sim.merge_time[sim.lineage[c]]
            ftr.write(f"{name}\t{chrom}\t{leaf}\t{sim.chrom[c]}\tfam{f:04d}\t{sim.superfamily[f]}\t{sim.clade[f]}\t"
                      f"{sim.tins[c]:.6g}\t{sim.tb[c]:.6g}\t{sim.origin[c]}\t{mode}\t{'yes' if post else 'no'}\t"
                      f"{sim.exid[c]}\t{sim.root[c]}\t{lt}\t{ts}\t{tv}\t{_fmt(d)}\n")
    with open(truthdir / "orthologs.tsv", "w") as f:
        f.write("copy_a\tcopy_b\tt_insert\tsplit_time\n")
        by_root = defaultdict(list)
        for c in surv:
            if sim.exid[c] < 0:
                by_root[sim.root[c]].append(c)
        for v in by_root.values():
            for i in range(len(v)):
                for j in range(i + 1, len(v)):
                    a, b = v[i], v[j]
                    if sim.lineage[a] != sim.lineage[b]:
                        f.write(f"{cid[a]}\t{cid[b]}\t{sim.tins[a]:.6g}\t{_mrca_split(sim, sim.lineage[a], sim.lineage[b]):.6g}\n")
    with open(truthdir / "paralogs.tsv", "w") as f:
        f.write("copy_dup\tcopy_src\texchange_id\texchange_time\n")
        for c in surv:
            p = sim.par[c]
            if sim.exid[c] >= 0 and sim.mode[c] == DUP and p in cid:
                f.write(f"{cid[c]}\t{cid[p]}\t{sim.exid[c]}\t{sim.tb[c]:.6g}\n")
    with open(truthdir / "exchanges.tsv", "w") as f:
        f.write("index\ttype\ttime\trecipient\tdonor\tchrom\tstart\tend\tn_moved\n")
        for e in sim.exchanges:
            R = sim.leaf_of[e["recipient"]]
            rows = [(R, e["recipient"], e["donor"])]
            if e["type"] == "reciprocal":
                rows.append((sim.leaf_of[e["donor"]], e["donor"], e["recipient"]))
            for l, rec, don in rows:
                f.write(f"{e['index']}\t{e['type']}\t{e['time']:.6g}\t{rec}\t{don}\t{tok[(l, e['chrom'])]}\t"
                        f"{e['start']}\t{e['end']}\t{e['n_moved']}\n")
    true_size = defaultdict(int)
    for c in surv:
        true_size[sim.fam[c]] += 1
    hist = dict(scenario=sim.sc, seed=sim.seed,
                splits=[dict(node=l.name, time=l.t_to, children=[sim.lin[c].name for c in l.children])
                        for l in sim.lin if l.children],
                mergers=sim.sc["mergers"], merge_time={sim.lin[l].name: t for l, t in sim.merge_time.items()},
                families=[dict(name=f"fam{f:04d}", ltr_len=int(sim.ltr_len[f]), int_len=int(sim.int_len[f]),
                               mu=float(sim.mu_f[f]), superfamily=sim.superfamily[f], clade=sim.clade[f],
                               n_copies=true_size.get(f, 0)) for f in range(len(sim.fam_w))],
                n_coordinate_collisions_dropped=n_dup,
                n_copies={sim.lin[l].name: sum(sim.lineage[c] == l for c in surv) for l in sim.leaves})
    (truthdir / "history.json").write_text(json.dumps(hist, indent=1) + "\n")
    return len(surv)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("scenario", help="scenario JSON (keys override DEFAULTS)")
    ap.add_argument("outdir", help="tool inputs: ltr.fa, k2p.tsv, annotation.gff3, genome.fai")
    ap.add_argument("truthdir", help="truth tables")
    ap.add_argument("--seed", type=int, default=1, help="random seed")
    a = ap.parse_args(argv)
    sc = merge_defaults(json.loads(Path(a.scenario).read_text()))
    sim = Sim(sc, a.seed).run()
    n = write_outputs(sim, a.outdir, a.truthdir)
    print(f"{Path(a.scenario).stem} seed {a.seed}: {n} LTR-RTs, {len(sim.exchanges)} exchanges")


if __name__ == "__main__":
    main()
