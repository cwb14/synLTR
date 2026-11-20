"""Collinear chains of shared insertions. The pairing map finds near-identical sequences (homologues of a young
polyploid, or haplotypes), kept apart when splitting. Homoeologous chromosomes share the insertions their common
ancestor carried before the progenitors split: such an old copy and its ortholog are each other's closest relative
across the two subgenomes, at collinear positions. Most copies are younger (subgenome-specific, spread over every
chromosome of their subgenome) and say nothing about which chromosome pairs with which; the inherited ones are deleted
over time (TE turnover), so old splits leave too few."""
from __future__ import annotations

from collections import defaultdict

import numpy as np
from scipy import sparse

from ltrtrace.analyses.subgenomes import inputs
from ltrtrace.analyses.subgenomes.inputs import group_label, log, warn
from ltrtrace.analyses.subgenomes.partition import bipartition, rep_stats, side_support, split_contrast
from ltrtrace.base.constants import (
    CHAIN_MIN_CANDIDATES,
    HALF_PAIRS,
    HOMOEOLOG_ALPHA,
    MAX_PAIRS,
    MIN_ANCHORS,
    MIN_REP,
    NMAX_KNN,
    PAIR_SHARED,
    PAIR_SPAN,
    REL_BLOCK_NNZ,
    SE_MULT,
)
from ltrtrace.base.unionfind import find
from ltrtrace.seq import align, divergence
from ltrtrace.seq.chain import anchor_ranks, chain_null, chain_significant, chain_test, chain_z, longest_chain
from ltrtrace.seq.relatives import best_by_group, mutual_best, rbh_anchors


def pairing_map(lib, lm, rng, threads):
    """Sequence pairs joined by collinear chains of mutual-best relatives: insertions shared at matching positions.
    Chain length vs NULL_PERMS permutations of partner positions -> z; significance by chain_significant. Anchors
    of significant chains are aligned; a pair is 'paired' (near-identical sequences: homologues of a young
    polyploid, or haplotypes) when (1) the chain spans >= PAIR_SPAN of both sequences' copies (a local block is a
    segmental duplication or an exchange), (2) most anchors are shared insertions (between-copy K2P below both
    copies' own LTR K2P), and (3) the anchors' median K2P is below the median within-copy K2P of the two sequences'
    copies (they diverged after most of their insertions)."""
    C = lm.C
    pairs = mutual_best(lm.rel_j[:, :, 0], lm.ek)                    # each the other's closest relative there
    if len(pairs) == 0 or C < 2:
        return []
    rank = np.zeros(lib.n, np.int64)
    for s in range(C):
        idx = np.flatnonzero(lm.ek == s)
        rank[idx[np.argsort(lib.start[idx], kind="stable")]] = np.arange(idx.size)
    sa, sb = lm.ek[pairs[:, 0]], lm.ek[pairs[:, 1]]
    A = np.where(sa <= sb, pairs[:, 0], pairs[:, 1])
    Bp = np.where(sa <= sb, pairs[:, 1], pairs[:, 0])
    key = np.minimum(sa, sb) * C + np.maximum(sa, sb)
    o = np.lexsort((rank[A], key))
    key, A, Bp = key[o], A[o], Bp[o]
    starts = np.r_[0, np.flatnonzero(np.diff(key)) + 1]
    ends = np.r_[starts[1:], key.size]
    n_tests = max(C * (C - 1) // 2, 1)
    n_cp = np.bincount(lm.ek[lm.ek >= 0], minlength=C)
    rows, anchors = [], []
    for s0, e0 in zip(starts, ends):
        if e0 - s0 < CHAIN_MIN_CANDIDATES:
            continue
        x, y = rank[A[s0:e0]], rank[Bp[s0:e0]]
        ch = longest_chain(x, y)
        null = chain_null(x, y, rng)
        z = chain_z(ch.size, null)
        i, j = int(key[s0] // C), int(key[s0] % C)
        xa, yb = x[ch], y[ch]
        rows.append(dict(i=i, j=j, n_candidates=int(e0 - s0), chain=int(ch.size), null_mean=float(null.mean()),
                         z=float(z), significant=chain_significant(int(e0 - s0), int(ch.size), z, n_tests),
                         span_a=float((xa.max() - xa.min() + 1) / n_cp[i]), span_b=float((yb.max() - yb.min() + 1) / n_cp[j]),
                         shared=float("nan"), d_median=float("nan"), k_median=float("nan"), paired=False))
        anchors.append(np.stack([A[s0:e0][ch], Bp[s0:e0][ch]], 1))
    sig = [k for k, r in enumerate(rows) if r["significant"]]
    if sig:
        cnt = align.pairs(lib.store, np.concatenate([anchors[k] for k in sig]), mode="ltr", threads=threads)
        d, _ = divergence.k2p(cnt.sites, cnt.ts, cnt.tv)
        at = 0
        for k in sig:
            r, m = rows[k], len(anchors[k])
            dk = d[at:at + m]
            at += m
            within = lib.k2p[(lm.ek == r["i"]) | (lm.ek == r["j"])]
            if np.isfinite(dk).any() and np.isfinite(within).any():
                own = np.minimum(lib.k2p[anchors[k][:, 0]], lib.k2p[anchors[k][:, 1]])
                with np.errstate(invalid="ignore"):
                    r["shared"] = float(np.mean(dk < own))                     # NaN compares False
                r["d_median"], r["k_median"] = float(np.nanmedian(dk)), float(np.nanmedian(within))
                r["paired"] = bool(min(r["span_a"], r["span_b"]) >= PAIR_SPAN and r["shared"] > PAIR_SHARED
                                   and r["d_median"] < r["k_median"])
    return rows


def paired_sets(rows, C):
    """Cannot-link sets from the pairing map: connected groups of 'paired' sequences of exactly two (larger
    groups are reported only: a bipartition cannot separate three copies of a chromosome)."""
    parent = list(range(C))
    for r in rows:
        if r["paired"]:
            parent[find(parent, r["i"])] = find(parent, r["j"])
    comp = defaultdict(list)
    for x in range(C):
        comp[find(parent, x)].append(x)
    return [sorted(v) for v in comp.values() if len(v) == 2], [sorted(v) for v in comp.values() if len(v) > 2]


def mask_pairs(W, pairs):
    """Link matrix without the links between paired sequences: shared insertions join homologues whatever their
    subgenome, so those links carry no subgenome signal."""
    if not pairs:
        return W
    W = W.copy()
    for x, y in pairs:
        W[x, y] = W[y, x] = 0.0
    return W


def pair_splits(src, lab, cl_pairs, rng, support):
    """A final group holding >= 2 paired sequence pairs cannot be one subgenome (a subgenome carries one copy of each
    chromosome): split it with every such pair separated, on pooled links without the paired links, and keep the split
    when half-libraries reproduce it (adjusted Rand >= MIN_REP, more pairs when borderline). Returns the new labels
    and the tests; the support of split sequences is updated in place."""
    lab, tests, k_next = lab.copy(), [], int(lab.max()) + 1
    for k in range(int(lab.max()) + 1):
        g = np.flatnonzero(lab == k)
        inside = [(x, y) for x, y in cl_pairs if lab[x] == k and lab[y] == k]
        if len(inside) < 2:
            continue
        pos = {x: i for i, x in enumerate(g)}
        sets = [[pos[x], pos[y]] for x, y in inside]
        Wm = mask_pairs(src.W["all"], inside)
        hm = [(mask_pairs(a, inside), mask_pairs(b, inside)) for a, b in src.pairs["all"]]
        z = bipartition(Wm, g, rng, sets)
        rep, se, hz = rep_stats(hm, g, rng, sets)
        while len(hm) < MAX_PAIRS and abs(rep - MIN_REP) < SE_MULT * se:
            hm += [(mask_pairs(a, inside), mask_pairs(b, inside)) for a, b in src.more("all", HALF_PAIRS)]
            rep, se, hz = rep_stats(hm, g, rng, sets)
        ok = rep >= MIN_REP
        c = dict(src="all", z=z, con=split_contrast(Wm, g, z), rep=rep, se=se, n_pairs=len(hm),
                 borderline=abs(rep - MIN_REP) < SE_MULT * se, hz=hz)
        tests.append(dict(members=g, path=f"paired{k}", cands={"all": c}, parent=None, chosen=c if ok else None,
                          how="pairs" if ok else "pairs: not reproduced", accepted=ok, forced=False, inside=inside))
        if ok:
            lab[g[z == 1]] = k_next
            k_next += 1
            support[g] = np.fmin(support[g], side_support(z, hz))
    return lab, tests


def old_relatives(cl, old, grp, N, threads=1, block_nnz=REL_BLOCK_NNZ):
    """For every copy in old: its closest relative among the old copies of each group (-1: no shared lineage).
    Relatedness = lineages of <= NMAX_KNN carriers in the whole library shared by the two copies, each weighted by its
    specificity mult/(size-1) (top_relatives' measure), counted over old copies only; equal relatedness to the lower
    copy index (seq.relatives.best_by_group, tie "lower"). Blocks of copies run in the fork pool."""
    n = old.size
    best = np.full((n, N), -1, np.int64)
    keep = np.repeat(cl.size <= NMAX_KNN, cl.size) & old[cl.members]
    sid = np.repeat(np.arange(cl.n), cl.size)[keep]
    mem = cl.members[keep].astype(np.int64)
    ok = np.bincount(sid, minlength=cl.n)[sid] >= 2                 # lineages carried by >= 2 old copies
    sid, mem = sid[ok], mem[ok]
    if not mem.size:
        return best
    u, row = np.unique(sid, return_inverse=True)
    w = cl.mult[u] / np.maximum(cl.size[u] - 1, 1)
    A = sparse.csr_matrix((np.ones(mem.size), (row, mem)), shape=(u.size, n))
    DA = sparse.csr_matrix((w[row], (row, mem)), shape=(u.size, n))
    J, _ = best_by_group(A, DA, grp, top=1, exclude_same=False, tie="lower", n_groups=N, block_nnz=block_nnz,
                         threads=threads)
    return J[:, :, 0].astype(np.int64)


def run_homoeologs(a, lib, ph):
    """Homoeologous chromosomes from shared old insertions (module docstring). Per chromosome i and other subgenome h: the
    candidate whose anchors with i form the most collinear chain (z against NULL_PERMS permutations) is called i's
    homoeolog when it also holds the most anchors relative to expectation (degree-corrected within the subgenome pair)
    and its chain passes chain_test at Bonferroni 0.05 over the chromosomes of h; else 'unresolved'; 'none' without a
    candidate (>= MIN_ANCHORS anchors). None for one group."""
    N, C = ph["N"], ph["C"]
    edges = ph.get("strata_edges")
    if N < 2 or not edges:
        return None
    z, ek = np.asarray(ph["z"]), np.asarray(ph["ek"])
    old = (ek >= 0) & np.isfinite(lib.k2p) & (lib.k2p > edges[1])
    if inputs.verbose() and ph.get("strata_n", {}).get("old") not in (None, int(old.sum())):
        warn(f"old copies {int(old.sum())} != the phasing step's old tercile {ph['strata_n']['old']}")
    grp = np.where(old, z[np.maximum(ek, 0)], -1)
    P = rbh_anchors(old_relatives(ph["clades_all"], old, grp, N, a.threads), grp)
    si, sj = ek[P[:, 0]], ek[P[:, 1]]
    M = np.zeros((C, C))
    np.add.at(M, (si, sj), 1)
    M = M + M.T
    E = np.zeros((C, C))
    for g in range(N):
        for h in range(N):
            ia, ib = np.flatnonzero(z == g), np.flatnonzero(z == h)
            if g != h and ia.size and ib.size:
                B = M[np.ix_(ia, ib)]
                E[np.ix_(ia, ib)] = np.outer(B.sum(1), B.sum(0)) / max(B.sum(), 1e-300)
    S = np.log((M + 1) / (E + 1))
    rng = np.random.default_rng([a.seed, 9])                   # own stream: no phasing or dating draw changes
    key = np.minimum(si, sj) * C + np.maximum(si, sj)
    o = np.argsort(key, kind="stable")
    st = np.r_[0, np.flatnonzero(np.diff(key[o])) + 1] if o.size else np.zeros(0, np.int64)
    en = np.r_[st[1:], o.size] if o.size else np.zeros(0, np.int64)
    pair, in_chain = {}, np.zeros(len(P), bool)
    for s0, e0 in zip(st, en):
        idx = o[s0:e0]
        if idx.size < MIN_ANCHORS:
            continue
        i, j = int(key[idx[0]] // C), int(key[idx[0]] % C)
        ci = np.where(si[idx] == i, P[idx, 0], P[idx, 1])                 # anchor copies on i and on j
        cj = np.where(si[idx] == i, P[idx, 1], P[idx, 0])
        oi, x, y = anchor_ranks(lib.start[ci], lib.start[cj])
        idx, ci, cj = idx[oi], ci[oi], cj[oi]
        ch = longest_chain(x, y)
        null = chain_null(x, y, rng)
        zc = chain_z(ch.size, null)
        in_chain[idx[ch]] = True
        span = lambda c_, s_: float((lib.end[c_].max() - lib.start[c_].min() + 1) / max(lib.seq_len[s_], 1))
        rec = dict(n=int(idx.size), chain=int(ch.size), z=float(zc), x=x, y=y,
                   orientation="+" if ch.size < 2 or y[ch[-1]] > y[ch[0]] else "-")
        si_i, si_j = lib.seq[ci[0]], lib.seq[cj[0]]
        pair[(i, j)] = dict(rec, span=span(ci[ch], si_i), partner_span=span(cj[ch], si_j))
        pair[(j, i)] = dict(rec, span=span(cj[ch], si_j), partner_span=span(ci[ch], si_i))
    calls, tested = {}, {}
    for i in range(C):
        for h in range(N):
            if h == z[i]:
                continue
            js = np.flatnonzero(z == h)
            zs = np.array([pair[(i, j)]["z"] if (i, j) in pair else -np.inf for j in js])
            if not js.size or not np.isfinite(zs).any():
                calls[(i, h)] = dict(status="none", partner=-1, candidate=-1, p=float("nan"))
                continue
            jz, jo = int(js[np.argmax(zs)]), int(js[np.argmax(S[i, js])])
            alpha = HOMOEOLOG_ALPHA / js.size
            r = pair[(i, jz)]
            if (min(i, jz), max(i, jz), alpha) not in tested:
                tested[(min(i, jz), max(i, jz), alpha)] = chain_test(r["x"], r["y"], r["chain"], alpha, rng)
            p = tested[(min(i, jz), max(i, jz), alpha)]
            ok = jz == jo and p <= alpha
            calls[(i, h)] = dict(status="homoeolog" if ok else "unresolved", partner=jz if ok else -1, candidate=jz,
                                 p=float(p))
    for (i, h), c in calls.items():
        j = c["partner"]
        c["reciprocal"] = (None if j < 0 else
                           calls.get((j, int(z[i])), {}).get("partner") == i)
    sets = homoeolog_sets(calls, z, N, C)
    called = {(i, c["partner"]) for (i, h), c in calls.items() if c["partner"] >= 0}
    pair_called = np.array([(int(a_), int(b_)) in called or (int(b_), int(a_)) in called for a_, b_ in zip(si, sj)],
                           bool)
    out = dict(P=P, si=si, sj=sj, M=M, E=E, S=S, pair=pair, calls=calls, sets=sets, in_chain=in_chain,
               pair_called=pair_called, old_edge=float(edges[1]), n_old=int(old.sum()))
    for g in range(N):
        for h in range(g + 1, N):
            ng, nh = int((z == g).sum()), int((z == h).sum())
            cg = sum(calls[(i, h)]["status"] == "homoeolog" for i in np.flatnonzero(z == g))
            chh = sum(calls[(j, g)]["status"] == "homoeolog" for j in np.flatnonzero(z == h))
            rc = sum(bool(calls[(i, h)]["reciprocal"]) for i in np.flatnonzero(z == g))
            log(f"homoeologs {group_label(g)}-{group_label(h)}: {cg} of {ng} {group_label(g)} and {chh} of {nh} "
                f"{group_label(h)} sequences paired ({rc} reciprocal pairs)")
    log(f"homoeologs: {len(P):,} shared old insertions (reciprocal closest relatives among {int(old.sum()):,} copies "
        f"with LTR K2P > {edges[1]:.4g}); {len(sets)} homoeolog sets, {sum(s_['complete'] for s_ in sets)} complete")
    return out


def homoeolog_sets(calls, z, N, C):
    """Connected components of the homoeolog calls (either direction), ordered by their first member: members, and
    complete = one member per subgenome with every cross pair called."""
    parent = list(range(C))
    edges = {(min(i, c["partner"]), max(i, c["partner"])) for (i, h), c in calls.items() if c["partner"] >= 0}
    for i, j in edges:
        parent[find(parent, i)] = find(parent, j)
    comp = defaultdict(list)
    for i, j in edges:
        comp[find(parent, i)].append(i)
        comp[find(parent, i)].append(j)
    out = []
    for members in comp.values():
        m = sorted(set(members))
        complete = (sorted(int(z[x]) for x in m) == list(range(N)) and
                    all((min(x, y), max(x, y)) in edges for x in m for y in m if x < y))
        out.append(dict(members=m, complete=bool(complete)))
    return sorted(out, key=lambda s_: (min(int(z[x]) for x in s_["members"]), s_["members"]))
