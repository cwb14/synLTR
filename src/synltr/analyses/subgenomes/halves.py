"""The number of subgenomes: a divisive search on the links of all copies and of each LTR-age tercile. A split is kept
when disjoint half-libraries reproduce it and it is strong enough; a split seen in one age class only must be
age-bounded (progenitor confinement lasts a window of insertion ages)."""
from __future__ import annotations

import numpy as np

from synltr.analyses.subgenomes.links import expected
from synltr.analyses.subgenomes.partition import (
    align,
    ari,
    best_partition,
    bipartition,
    rep_stats,
    side_support,
    split_contrast,
)
from synltr.base.constants import (
    AGE_STRATA,
    CORE_FRAC,
    CORE_MIN,
    HALF_PAIRS,
    MAX_N,
    MAX_PAIRS,
    MIN_GROUP,
    MIN_REL,
    MIN_REP,
    MIN_REP_AGE,
    SE_MULT,
)


class Sources:
    """Link matrices of all phased copies ('all') and of each LTR-K2P tercile (AGE_STRATA), with complementary
    half-library pairs drawn from the same random halves for every source; more(s, n) adds n pairs to source s."""

    def __init__(self, lm, k2p, rng):
        self.lm, self.rng = lm, rng
        ph = lm.ek >= 0
        ok = ph & np.isfinite(k2p)
        self.masks = {"all": ph}
        self.edges = None
        if ok.sum() >= len(AGE_STRATA):
            q = np.quantile(k2p[ok], np.linspace(0, 1, len(AGE_STRATA) + 1)[1:-1])
            self.edges = [float(x) for x in q]
            lo = -np.inf
            for name, hi in zip(AGE_STRATA, list(q) + [np.inf]):
                self.masks[name] = ok & (k2p > lo) & (k2p <= hi)
                lo = hi
        self.names = list(self.masks)
        self.W = {s: lm.matrix(m) for s, m in self.masks.items()}
        self.pairs = {s: [] for s in self.names}
        for seed in rng.integers(0, 2 ** 62, HALF_PAIRS):
            u = np.random.default_rng(int(seed)).random(lm.n) < 0.5
            for s, m in self.masks.items():
                self.pairs[s].append(self._pair(u, m))

    def _pair(self, u, m):
        return self.lm.matrix(u & m).astype(np.float32), self.lm.matrix(~u & m).astype(np.float32)

    def more(self, s, n):
        new = [self._pair(np.random.default_rng(int(sd)).random(self.lm.n) < 0.5, self.masks[s])
               for sd in self.rng.integers(0, 2 ** 62, n)]
        self.pairs[s].extend(new)
        return new


def split_candidate(src, s, g, rng, sets=None, min_size=MIN_GROUP):
    """Best bipartition of group g on source s's links, its strength, and its replicability in s's half-library pairs
    (more pairs drawn while the mean straddles MIN_REP within 2 s.e., up to MAX_PAIRS); None without links or when a
    side would hold fewer than min_size sequences."""
    W = src.W[s]
    if W[np.ix_(g, g)].sum() <= 0:
        return None
    z = bipartition(W, g, rng, sets)
    if min(np.bincount(z, minlength=2)) < min_size:
        return None
    rep, se, hz = rep_stats(src.pairs[s], g, rng, sets)
    while len(src.pairs[s]) < MAX_PAIRS and abs(rep - MIN_REP) < SE_MULT * se:
        src.more(s, HALF_PAIRS)
        rep, se, hz = rep_stats(src.pairs[s], g, rng, sets)
    return dict(src=s, z=z, con=split_contrast(W, g, z), rep=rep, se=se, n_pairs=len(src.pairs[s]),
                borderline=abs(rep - MIN_REP) < SE_MULT * se, hz=hz)


def parent_info(src, g, z, cands):
    """What children are measured against, on each source: the strongest split of the parent group there (the
    accepted split or that source's own best split, whichever is stronger)."""
    con = {}
    for s in src.names:
        W = src.W[s]
        own = split_contrast(W, g, z) if W[np.ix_(g, g)].sum() > 0 else 0.0
        con[s] = max(own, cands[s]["con"]) if s in cands else own
    return dict(con=con)


def accept_split(cands, parent, strata_avail=None):
    """A split passes on some copies when half-libraries of those copies reproduce it (adjusted Rand >= MIN_REP), it
    has positive strength and, below the top, it is >= MIN_REL as strong as the strongest split of its parent group on
    the same copies. Pooled copies decide first. Otherwise one LTR-K2P tercile can carry it, provided the split is
    age-bounded: some tercile holds no reproducible split of the group (confinement to a progenitor has a window;
    structure present at every age is the chromosomes' own).
    The same age rule applies to the top split. Copies inserted between the progenitors' split and their merger are the
    only subgenome-specific ones; pooled with older (inherited) and younger (post-merger) copies their signal is
    diluted, so an old polyploid can show its subgenomes in one age window only. With no parent to compare strength
    with, the tercile's split must have positive strength and be the pooled copies' best split too. At every level a
    split carried by one age class must replicate at MIN_REP_AGE (three extra tests); true age-confined splits
    replicate at 0.97-1.00 in the benchmarks, within-subgenome structure of tetraploids at 0.80-0.83. strata_avail: the
    age strata the sources hold (the top has no parent to list them).
    -> (candidate, 'pooled' | 'age') or (None, reason)."""
    def strong(c):
        if c["rep"] < MIN_REP or c["con"] <= 0:
            return False
        return parent is None or (parent["con"][c["src"]] > 0 and c["con"] >= MIN_REL * parent["con"][c["src"]])
    p = cands.get("all")
    if p is not None and strong(p):
        return p, "pooled"
    if parent is None:
        strata = [s for s in AGE_STRATA if s in (strata_avail or ())]
        seen = [s for s in strata if cands.get(s) is not None and cands[s]["rep"] >= MIN_REP]
        ok = [cands[s] for s in seen if strong(cands[s]) and cands[s]["rep"] >= MIN_REP_AGE and
              (p is not None and ari(cands[s]["z"], p["z"]) >= MIN_REP)]
        if ok and len(seen) < len(strata):
            return max(ok, key=lambda c: (round(c["rep"], 2), c["con"])), "age"
        if ok:
            return None, "present at every age"
        return None, "not reproduced" if p is None or p["rep"] < MIN_REP else "no excess"
    strata = [s for s in AGE_STRATA if s in parent["con"]]
    seen = [s for s in strata if cands.get(s) is not None and cands[s]["rep"] >= MIN_REP]
    ok = [cands[s] for s in seen if strong(cands[s]) and cands[s]["rep"] >= MIN_REP_AGE]     # the strict bar
    if ok and len(seen) < len(strata):
        return max(ok, key=lambda c: (round(c["rep"], 2), c["con"] / parent["con"][c["src"]])), "age"
    if ok:
        return None, "present at every age"
    if (p is not None and p["rep"] >= MIN_REP) or seen:
        return None, "weaker than half its parent"
    return None, "not reproduced"


def core_sequences(n_el, frac=CORE_FRAC, min_core=CORE_MIN):
    """Sequences that decide the structure: >= frac x the copy-N50 LTR-RTs (the copy count such that sequences
    holding at least that many copies hold half of all phased copies). Chromosomes of a chromosome-level assembly carry
    hundreds to thousands of copies; unplaced scaffolds tens. Fewer than min_core such sequences (or a fragmented
    assembly without a break): every phased sequence is core."""
    n_el = np.asarray(n_el, float)
    if n_el.size == 0:
        return np.zeros(0, bool)
    s = np.sort(n_el)[::-1]
    c = np.cumsum(s)
    n50 = s[min(int(np.searchsorted(c, c[-1] / 2.0)), s.size - 1)]
    core = n_el >= frac * n50
    return core if core.sum() >= min_core else np.ones(n_el.size, bool)


def assign_satellites(src, lab, sat):
    """Place each non-core phased sequence in the group it shares the most excess links with (the local-search
    move on W - E over all phased sequences, groups fixed by the core); support = share of the half-library matrices
    that give the same group. Satellites never create or split a group. Returns labels for all sequences and the
    satellites' support (0 for a satellite without links to any core sequence)."""
    lab = np.asarray(lab).copy()
    sup = np.full(lab.size, np.nan)
    sat = np.asarray(sat, np.int64)
    core = np.flatnonzero(lab >= 0)
    N = int(lab[core].max()) + 1 if core.size else 0
    if not sat.size or N == 0:
        return lab, sup
    members = [core[lab[core] == g] for g in range(N)]

    def best(Wm):
        Wm = np.asarray(Wm, np.float64)
        B = Wm - expected(Wm)
        S = np.stack([B[np.ix_(sat, m)].sum(1) for m in members], 1)
        linked = Wm[np.ix_(sat, core)].sum(1) > 0
        return S.argmax(1), linked
    b, linked = best(src.W["all"])
    lab[sat] = b
    hits = [best(h)[0] == b for pair in src.pairs["all"] for h in pair]
    sup[sat] = np.mean(hits, 0) if hits else np.nan
    sup[sat[~linked]] = 0.0
    return lab, sup


def find_subgenomes(src, rng, n_force=None, min_size=MIN_GROUP, max_n=MAX_N, sets_of=None, root=None, age_ok=None):
    """Divisive search on all copies and on each LTR-K2P tercile. Inferred N: every group is tested once, splits are
    kept by accept_split, N=1 is a possible answer. Forced N (n_force): the group whose best candidate ranks first
    (accepted, then replicability, then strength relative to its parent) is split until there are n_force groups.
    sets_of(g) -> cannot-link pairs (local indices) for the bipartitions of g. Returns labels, every tested split,
    and per-sequence support (the minimum, over the splits that placed it, of its side-support).
    root = the sequences searched (the core); the others get label -1 (placed by assign_satellites).
    age_ok(g, z) -> bool: a sub-split kept by the age rule (one age class; the pooled copies do not show it) must have
    homoeologs across its parts (split_homoeologs); else it is rejected ('<class> split without homoeologs across
    it'). Oryza alta (a young 7|6 cut of one subgenome) and a Rosa gallica control took such splits. None: no check."""
    C = src.W["all"].shape[0]
    root = np.arange(C) if root is None else np.asarray(root, np.int64)
    strata_avail = [s for s in src.names if s != "all"]
    support = np.full(C, np.nan)
    tests = []

    def evaluate(g, parent, path):
        sets = sets_of(g) if sets_of else None
        cands = {}
        if len(g) >= 2 * min_size:
            for s in src.names:
                c = split_candidate(src, s, g, rng, sets, min_size)
                if c is not None:
                    cands[s] = c
        ch, how = (accept_split(cands, parent, strata_avail) if cands else
                   (None, "too small" if len(g) < 2 * min_size else "no split"))
        gated = None
        if ch is not None and how == "age" and parent is not None and age_ok is not None and not age_ok(g, ch["z"]):
            gated, ch, how = ch, None, f"{ch['src']} split without homoeologs across it"
        t = dict(members=g, path=path, cands=cands, parent=parent, chosen=ch, how=how, accepted=ch is not None,
                 forced=False, gated=gated)
        tests.append(t)
        return t

    def rank(t):
        c = t["chosen"] or max(t["cands"].values(), key=lambda c: (round(c["rep"], 2), c["con"]))
        par = t["parent"]["con"][c["src"]] if t["parent"] else 1.0
        return (t["accepted"], round(c["rep"], 2), c["con"] / max(par, 1e-12)), c

    leaves = [evaluate(root, None, "root")]
    final = []
    while leaves:
        if n_force is None:
            t = leaves.pop(0)
            if not t["accepted"] or len(final) + len(leaves) + 2 > max_n:
                final.append(t["members"])
                continue
            c = t["chosen"]
        else:
            open_ = [x for x in leaves if x["cands"]]
            if len(final) + len(leaves) >= n_force or not open_:
                final += [x["members"] for x in leaves]
                break
            t = max(open_, key=lambda x: rank(x)[0])
            leaves = [x for x in leaves if x is not t]              # identity: tests hold arrays
            if not t["accepted"]:
                t["chosen"], t["forced"] = rank(t)[1], True
            c = t["chosen"]
        g, z = t["members"], c["z"]
        sup = side_support(z, c["hz"])
        support[g] = np.fmin(support[g], sup)                      # NaN (not yet placed) gives way
        pi = parent_info(src, g, z, t["cands"])
        leaves += [evaluate(g[z == 0], pi, t["path"] + ".a"), evaluate(g[z == 1], pi, t["path"] + ".b")]
    lab = np.full(C, -1, int)
    for k, g in enumerate(sorted(final, key=lambda x: x.min())):
        lab[g] = k
    return lab, tests, support


def support_from_halves(halves, z, N, sets, rng):
    """Fraction of half-library partitions (solved from scratch with the full-data estimator) agreeing with z."""
    hits = []
    for W in halves:
        if W.sum() <= 0:
            continue
        zr = best_partition(W - expected(W), N, sets, rng)
        hits.append(align(z, zr, N) == z)
    return np.mean(hits, 0) if hits else np.full(len(z), np.nan)


def partner_pairs(W, halves, z):
    """Closest partner of each sequence within its own group (most shared lineages, observed/expected), kept when
    mutual; support = fraction of halves where it still is; excess = the pair's O/E over the sequence's median
    within its group (near-identical homologues or progenitors show large excess)."""
    C = W.shape[0]
    same = z[:, None] == z[None, :]
    np.fill_diagonal(same, False)

    def best_of(M):
        E = expected(M)
        oe = np.where(same & (E > 0), M / np.where(E > 0, E, 1.0), np.nan)
        ok = np.isfinite(oe).any(1)
        b = np.where(ok, np.nanargmax(np.where(np.isfinite(oe), oe, -np.inf), 1), -1)
        return b, oe
    b, oe = best_of(W)
    partner = np.array([j if j >= 0 and b[j] == i else -1 for i, j in enumerate(b)])
    support, excess = np.full(C, np.nan), np.full(C, np.nan)
    hb = [best_of(h)[0] for h in halves]
    for i, j in enumerate(partner):
        if j < 0:
            continue
        support[i] = float(np.mean([bb[i] == j and bb[j] == i for bb in hb])) if hb else np.nan
        rest = oe[i][np.isfinite(oe[i])]
        excess[i] = float(oe[i, j] / np.median(rest)) if rest.size > 1 and np.median(rest) > 0 else np.nan
    return partner, support, excess


def chrom_tests(W, z, N):
    """Per sequence: mean O/E to its own group and to the closest other group; one-sided Welch test with the other
    sequences as replicates."""
    from scipy.stats import ttest_ind
    C = len(z)
    E = expected(W)
    L = np.log((W + 1) / (E + 1))
    OE = W / np.where(E > 0, E, np.nan)
    own_oe, oth_oe, pval = (np.full(C, np.nan) for _ in range(3))
    for c in range(C):
        others = np.arange(C) != c
        x = L[c, others & (z == z[c])]
        best, bg = -np.inf, None
        for g in range(N):
            if g != z[c] and (z == g).any() and L[c, z == g].mean() > best:
                best, bg = L[c, z == g].mean(), g
        if bg is None:
            continue
        y = L[c, others & (z == bg)]
        own_oe[c] = np.nanmean(OE[c, others & (z == z[c])]) if x.size else np.nan
        oth_oe[c] = np.nanmean(OE[c, others & (z == bg)])
        if x.size >= 2 and y.size >= 2:
            t, p = ttest_ind(x, y, equal_var=False)
            pval[c] = p / 2 if t > 0 else 1 - p / 2
    return own_oe, oth_oe, pval
