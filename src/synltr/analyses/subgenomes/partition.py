"""Partitions of sequences that maximize within-group excess links (degree-corrected modularity; exact for two
groups), adjusted Rand agreement, and the replicability of a split across half-library pairs."""
from __future__ import annotations

import math

import numpy as np
from numba import njit
from scipy.optimize import linear_sum_assignment

from synltr.analyses.subgenomes.links import expected, log_oe
from synltr.base.constants import EXACT_MAX_SETS, KMEANS_ITERS, LOCAL_ROUNDS, PARTITION_RESTARTS, SPLIT_RESTARTS


@njit(cache=True)
def gray_best(M, allsingle):
    """Maximize s^T M s over s in {-1, +1}^V with s_0 = +1, visiting configurations in Gray-code order (one sign
    flips per step; the score and the local fields update in O(V)). With all sets singletons, all-plus (one group)
    is excluded."""
    V = M.shape[0]
    s = np.ones(V)
    h = np.zeros(V)
    q = 0.0
    for i in range(V):
        for w in range(V):
            q += M[i, w]
            if w != i:
                h[i] += M[i, w]
    best = -np.inf
    bs = s.copy()
    if not allsingle:
        best = q
    nneg = 0
    total = 1 << (V - 1)
    for g in range(1, total):
        j = 1
        x = g
        while (x & 1) == 0:
            x >>= 1
            j += 1
        sj = s[j]
        q += -4.0 * sj * h[j]
        s[j] = -sj
        nneg += 1 if sj > 0 else -1
        for w in range(V):
            if w != j:
                h[w] -= 2.0 * M[w, j] * sj
        if (not allsingle or nneg > 0) and q > best:
            best = q
            bs[:] = s
    return bs


def _exhaustive2(B, sets):
    C, V = B.shape[0], len(sets)
    U = np.zeros((C, V))
    for v, st in enumerate(sets):
        U[st[0], v] = 1.0
        if len(st) == 2:
            U[st[1], v] = -1.0
    s = gray_best(U.T @ B @ U, all(len(st) == 1 for st in sets))
    return ((U @ s) < 0).astype(int)


def q_score(B, z):
    return 0.5 * float(B[z[:, None] == z[None, :]].sum())


def _local_search(B, z, sets, N, rng, max_rounds=LOCAL_ROUNDS):
    z = z.copy()
    C = B.shape[0]
    for _ in range(max_rounds):
        changed = False
        for si in rng.permutation(len(sets)):
            mem = sets[si]
            oh = np.zeros((C, N))
            oh[np.arange(C), z] = 1
            oh[mem] = 0
            A = B[mem] @ oh
            if len(mem) == 1:
                new = [int(np.argmax(A[0]))]
            else:
                r, c = linear_sum_assignment(-A)
                new = list(c[np.argsort(r)])
            for c_, g in zip(mem, new):
                if z[c_] != g:
                    z[c_] = g
                    changed = True
        if not changed:
            break
    return z


def _spectral_init(B, N, rng):
    _, V = np.linalg.eigh(B)
    X = V[:, -max(N - 1, 1):]
    if N == 2:
        return (X[:, -1] < 0).astype(int)
    cen = X[rng.choice(len(X), N, replace=False)]
    for _ in range(KMEANS_ITERS):
        z = np.argmin(((X[:, None, :] - cen[None]) ** 2).sum(-1), 1)
        cen = np.array([X[z == g].mean(0) if (z == g).any() else X[rng.integers(len(X))] for g in range(N)])
    return z


def best_partition(B, N, sets, rng, restarts=PARTITION_RESTARTS, init=None):
    """Maximum within-group excess (degree-corrected modularity) partition into N groups (exact for N=2)."""
    C = B.shape[0]
    if N == 2 and all(len(s) <= 2 for s in sets) and len(sets) <= EXACT_MAX_SETS:
        return _exhaustive2(B, sets)
    inits = [_spectral_init(B, N, rng)]
    if init is not None:
        inits.insert(0, np.asarray(init))
    inits += [rng.integers(0, N, C) for _ in range(restarts)]
    best, bz = -np.inf, None
    for z0 in inits:
        z = _local_search(B, z0, sets, N, rng)
        if len(set(z.tolist())) < N:
            continue
        q = q_score(B, z)
        if q > best + 1e-12:
            best, bz = q, z
    return bz if bz is not None else _local_search(B, inits[0], sets, N, rng)


def canon_labels(z, mass):
    """Relabel groups so that group 0 holds the most LTR-RTs (ties: lowest index)."""
    groups = sorted(set(z.tolist()), key=lambda g: (-mass[z == g].sum(), np.flatnonzero(z == g)[0]))
    m = {g: i for i, g in enumerate(groups)}
    return np.array([m[g] for g in z])


def align(ref, z, N):
    cont = np.zeros((N, N))
    np.add.at(cont, (z, ref), 1)
    r, c = linear_sum_assignment(-cont)
    m = dict(zip(r, c))
    return np.array([m[g] for g in z])


def ari(a, b):
    """Adjusted Rand index (1 = identical, ~0 = chance)."""
    from scipy.special import comb
    a, b = np.asarray(a), np.asarray(b)
    _, ai = np.unique(a, return_inverse=True)
    _, bi = np.unique(b, return_inverse=True)
    ct = np.zeros((ai.max() + 1, bi.max() + 1))
    np.add.at(ct, (ai, bi), 1)
    s = comb(ct, 2).sum()
    sa, sb = comb(ct.sum(1), 2).sum(), comb(ct.sum(0), 2).sum()
    exp = sa * sb / comb(a.size, 2)
    mx = (sa + sb) / 2
    return 1.0 if mx == exp else float((s - exp) / (mx - exp))


def bipartition(W, idx, rng, sets=None):
    """Best 2-way split of sequences idx from links among them only (null recomputed in the group); sets: cannot-
    link pairs given as local indices into idx (each pair split across the two parts)."""
    Wg = np.asarray(W[np.ix_(idx, idx)], np.float64)
    if Wg.sum() <= 0 or len(idx) < 2:
        return np.zeros(len(idx), int)
    return best_partition(Wg - expected(Wg), 2, local_sets(sets, len(idx)), rng, restarts=SPLIT_RESTARTS)


def local_sets(pairs, n):
    """Cannot-link pairs (local indices) completed with singletons into sets covering 0..n-1."""
    pairs = [list(p) for p in (pairs or [])]
    seen = {x for p in pairs for x in p}
    return pairs + [[i] for i in range(n) if i not in seen]


def split_contrast(W, idx, z):
    """Mean log(O/E) within the two parts minus between them (null recomputed inside the group)."""
    L = log_oe(W[np.ix_(idx, idx)])
    same = z[:, None] == z[None, :]
    off = ~np.eye(len(idx), dtype=bool)
    if not (same & off).any() or not (~same).any():
        return 0.0
    return float(L[same & off].mean() - L[~same].mean())


def rep_stats(pairs, g, rng, sets=None):
    """Replicability of the split of group g: mean adjusted Rand between the best splits found independently in the
    two halves of each pair, its standard error, and the halves' splits (for per-sequence support)."""
    hz = [(bipartition(Wa, g, rng, sets), bipartition(Wb, g, rng, sets)) for Wa, Wb in pairs]
    v = np.array([ari(x, y) for x, y in hz])
    return float(v.mean()), (float(v.std(ddof=1) / math.sqrt(v.size)) if v.size > 1 else 0.0), hz


def side_support(z, hz):
    """Per member of a split group: fraction of half-library splits (each side-aligned to z) that put it on its side."""
    hits = []
    for pair in hz:
        for h in pair:
            h = np.asarray(h)
            hits.append((h if (h == z).mean() >= 0.5 else 1 - h) == z)
    return np.mean(hits, 0) if hits else np.full(len(z), np.nan)
