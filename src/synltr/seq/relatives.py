"""Closest relatives of library copies.

`best_by_group`: for every copy (a column of A) and every group, its `top` closest relatives in that group. A is a
sparse lineage x copy incidence matrix (one row per lineage) and D the same with each lineage's weight; the
relatedness of copies e and f is (A^T D)[e, f], the summed weights of the lineages they share. The self-join runs in
blocks of copies of about block_nnz products each, in up to `threads` processes; results do not depend on blocks or
threads. Equal relatedness is ordered by `tie`: "lower", the lower copy index first; "product", the order in which
the sparse product (scipy's, CSR) lists them, phasing's nearest-relative table. `mutual_best` and `rbh_anchors` pair
copies that are each other's closest relatives.

`top_shared`: for focal copies, the copies in other parts sharing the most sketch hashes (seq.kmers.sketch).

`reciprocal_best`: of weighted candidate pairs between two libraries (seq.kmers.cross_pairs), those that are each
other's heaviest partner (the library screen's anchors and orthologs), ties to the lower index.
"""
from __future__ import annotations

import numpy as np
from numba import njit, prange

from synltr.base.constants import BUCKET_CAP, REL_BLOCK_NNZ, REL_MIN_SHARED, REL_TOP
from synltr.base.runtime import n_chunks, pool_map, shared

TIES = ("lower", "product")


def _best_block(e01):
    """The closest relatives of copies e0..e1-1: one block of the sparse self-join, one pool job."""
    e0, e1 = e01
    AT, D, group, G, top, exclude_same, tie = (shared[x] for x in ("rel_AT", "rel_D", "rel_group", "rel_G", "rel_top",
                                                                  "rel_exclude_same", "rel_tie"))
    J = np.full((e1 - e0, G, top), -1, np.int32)
    V = np.zeros((e1 - e0, G, top), D.dtype)
    R = (AT[e0:e1] @ D).tocoo()
    r, c, v = R.row.astype(np.int64), R.col.astype(np.int64), R.data
    g = group[c]
    k = (g >= 0) & ((g != group[e0 + r]) if exclude_same else (c != e0 + r))
    r, c, v, g = r[k], c[k], v[k], g[k]
    if r.size:
        key = r * G + g
        o = np.lexsort((-v, key)) if tie == "product" else np.lexsort((c, -v, key))
        key, c, v = key[o], c[o], v[o]
        rank = np.arange(key.size) - np.searchsorted(key, key, side="left")
        m = rank < top
        J[key[m] // G, key[m] % G, rank[m]] = c[m]
        V[key[m] // G, key[m] % G, rank[m]] = v[m]
    return e0, J, V


def best_by_group(A, D, group, *, top, exclude_same, tie, n_groups=None, block_nnz=REL_BLOCK_NNZ, threads):
    """For every copy e and group h (0..n_groups-1; default max(group) + 1): its `top` closest relatives in h, J
    (copy indices, int32; -1 = none) and their relatedness V (D's dtype; 0 = none), each (copies, groups, top), most
    related first (module docstring). A relative is a copy with a group (>= 0) and some shared lineage: with
    exclude_same, in another group than e's; else any copy but e itself. tie: "lower" or "product"."""
    if tie not in TIES:
        raise ValueError(f"tie must be one of {TIES}, not {tie!r}")
    A, D, group = A.tocsr(), D.tocsr(), np.asarray(group)
    n = A.shape[1]
    G = int(n_groups) if n_groups is not None else int(group.max()) + 1 if group.size else 0
    rows = np.diff(A.indptr)
    cost = np.bincount(A.indices, weights=np.repeat(rows, rows).astype(float), minlength=n)
    cum = np.concatenate(([0], np.cumsum(cost)))
    blocks, e0 = [], 0
    while e0 < n:
        e1 = min(max(int(np.searchsorted(cum, cum[e0] + block_nnz, side="right")), e0 + 1), n)
        blocks.append((e0, e1))
        e0 = e1
    J = np.full((n, G, top), -1, np.int32)
    V = np.zeros((n, G, top), D.dtype)
    shared.update(rel_AT=A.T.tocsr(), rel_D=D, rel_group=group, rel_G=G, rel_top=top, rel_exclude_same=exclude_same,
                  rel_tie=tie)
    try:
        for e0, bj, bv in pool_map(_best_block, blocks, threads, warm=False):
            J[e0:e0 + len(bj)] = bj
            V[e0:e0 + len(bv)] = bv
    finally:
        for x in ("rel_AT", "rel_D", "rel_group", "rel_G", "rel_top", "rel_exclude_same", "rel_tie"):
            shared.pop(x, None)
    return J, V


def mutual_best(best, group):
    """Copy pairs (a, b), a < b, each the other's closest relative in the other's group: best[a, group[b]] = b and
    best[b, group[a]] = a (best: copies x groups, -1 = none; copies with group < 0 are not paired)."""
    e = np.flatnonzero(group >= 0)
    E = np.repeat(e, best.shape[1])
    F = best[e].ravel()
    keep = F > E
    E, F = E[keep], F[keep]
    mutual = best[F, group[E]] == E
    return np.stack([E[mutual], F[mutual]], 1)


def rbh_anchors(best, grp):
    """Reciprocal closest relatives across groups -> pairs (a, b), grp[a] < grp[b]: b is a's closest relative in b's
    group and a is b's in a's group (candidate orthologous insertions)."""
    out = [np.zeros((0, 2), np.int64)]
    for g in range(best.shape[1]):
        for h in range(g + 1, best.shape[1]):
            a = np.flatnonzero((grp == g) & (best[:, h] >= 0))
            b = best[a, h]
            keep = best[b, g] == a
            out.append(np.stack([a[keep], b[keep]], 1))
    return np.concatenate(out)


def reciprocal_best(a, b, w, n_a, n_b, *, tie="lower") -> np.ndarray:
    """Weighted candidate pairs (a[i] of n_a copies, b[i] of n_b; no pair twice) -> the indices i of the pairs that are
    reciprocal best: b[i] is a[i]'s heaviest partner and a[i] is b[i]'s. tie "lower": equal weights go to the lower
    partner index, on both sides, whatever the input order. One index per a at most, in a order (int64)."""
    if tie not in ("lower",):
        raise ValueError(f"tie must be 'lower', not {tie!r}")
    a, b, w = np.asarray(a, np.int64), np.asarray(b, np.int64), np.asarray(w)
    o = np.lexsort((b, -w, a))
    first = np.ones(o.size, bool)
    first[1:] = a[o][1:] != a[o][:-1]
    best_a = o[first]
    o = np.lexsort((a, -w, b))
    first = np.ones(o.size, bool)
    first[1:] = b[o][1:] != b[o][:-1]
    best_b = o[first]
    partner = np.full(n_b, -1, np.int64)
    partner[b[best_b]] = a[best_b]
    return best_a[partner[b[best_a]] == a[best_a]]


@njit(cache=True)
def _sketch_index(S, use):
    """(hash, copy) of the sketches of copies in use, sorted by hash."""
    n, s = S.shape
    pad = np.uint64(0xFFFFFFFFFFFFFFFF)
    tot = 0
    for i in range(n):
        if use[i]:
            for j in range(s):
                if S[i, j] != pad:
                    tot += 1
    H = np.empty(tot, np.uint64)
    C = np.empty(tot, np.int64)
    k = 0
    for i in range(n):
        if use[i]:
            for j in range(s):
                if S[i, j] != pad:
                    H[k] = S[i, j]
                    C[k] = i
                    k += 1
    o = np.argsort(H)
    return H[o], C[o]


@njit(parallel=True, cache=True)
def _top_by_index(S, seq, focal, Hs, Cs, m, cap, nchunk):
    n, s = S.shape
    nf = focal.size
    J = np.full((nf, m), -1, np.int64)
    V = np.zeros((nf, m), np.int64)
    pad = np.uint64(0xFFFFFFFFFFFFFFFF)
    for ch in prange(nchunk):
        cnt = np.zeros(n, np.int32)
        touched = np.empty(n, np.int64)
        for f in range(nf * ch // nchunk, nf * (ch + 1) // nchunk):
            i = focal[f]
            nt = 0
            for j in range(s):
                h = S[i, j]
                if h == pad:
                    break
                lo = np.searchsorted(Hs, h)
                hi = np.searchsorted(Hs, h, side="right")
                if hi - lo > cap:
                    continue
                for q in range(lo, hi):
                    c = Cs[q]
                    if seq[c] == seq[i]:
                        continue
                    if cnt[c] == 0:
                        touched[nt] = c
                        nt += 1
                    cnt[c] += 1
            for t in range(nt):
                c = touched[t]
                v = cnt[c]
                cnt[c] = 0
                if v < REL_MIN_SHARED or v < V[f, m - 1] or (v == V[f, m - 1] and J[f, m - 1] >= 0 and c > J[f, m - 1]):
                    continue
                q = m - 1
                while q > 0 and (V[f, q - 1] < v or (V[f, q - 1] == v and J[f, q - 1] > c)):
                    V[f, q] = V[f, q - 1]
                    J[f, q] = J[f, q - 1]
                    q -= 1
                V[f, q] = v
                J[f, q] = c
    return J, V


def top_shared(S, part, focal, *, m=REL_TOP, cap=BUCKET_CAP):
    """For each focal copy, the m copies of other parts (part >= 0; part < 0 = not searched) sharing the most hashes of
    their sketches S (rows of seq.kmers.sketch; >= 2 shared; hashes carried by more than cap such copies ignored: they
    are family-wide), most shared first, ties to the lower copy index (results do not depend on threads): copies J and
    shared hashes V, (focal, m) int64; -1 = none."""
    part = np.asarray(part, np.int64)
    Hs, Cs = _sketch_index(S, part >= 0)
    return _top_by_index(S, part, np.asarray(focal, np.int64), Hs, Cs, m, cap, n_chunks(1, len(focal)))
