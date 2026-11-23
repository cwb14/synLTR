"""Collinear chains of anchors in rank space.

An anchor joins a copy on one sequence to its partner on another; it is the point (x, y) of their ranks along the two
sequences (`anchor_ranks`). A chain is a run of points increasing in x and monotone in y (increasing, or decreasing for
an inversion), consecutive points at most W ranks apart on both (`chain_dp`); `longest_chain` is the longest one in
either orientation. Its significance: `chain_null`, its length when partner positions are permuted (z by `chain_z`,
decided by `chain_significant`), or `chain_test`, a p from the exact tail of the longest increasing subsequence of a
random permutation (`lis_tail`) when the window cannot bind, else from permutations. Every random draw comes from the
np.random.Generator the caller passes, one permutation of y per replicate.

LTRmix's collinear blocks (the library screen's): `chain_blocks` takes the longest chains out greedily, both
orientations, down to a minimum length set by `shuffled_chain`, the longest chain among the same anchors with their
partners' ranks shuffled (LTRmix's draws: one permutation of the anchors' indices per shuffle).
"""
from __future__ import annotations

import math
from functools import cache

import numpy as np
from numba import njit

from synltr.base.constants import CHAIN_ALPHA, CHAIN_PERMS_PER_ALPHA, CHAIN_SD_FLOOR, CHAIN_W, NULL_PERMS


@njit(cache=True)
def chain_dp(x, y, W, sgn):
    """Longest-chain DP over points sorted by x: y increasing (sgn 1) or decreasing (sgn -1), steps <= W in both.
    Returns each point's chain length and predecessor (-1: none)."""
    n = x.size
    L = np.ones(n, np.int64)
    P = np.full(n, -1, np.int64)
    for i in range(n):
        j = i - 1
        while j >= 0 and x[i] - x[j] <= W:
            dy = (y[i] - y[j]) * sgn
            if x[i] > x[j] and 0 < dy <= W and L[j] + 1 > L[i]:
                L[i] = L[j] + 1
                P[i] = j
            j -= 1
    return L, P


@njit(cache=True)
def chain_length(x, y, W):
    """Length of the longest chain (chain_dp) in either orientation; 0 without points."""
    m = 0
    for sgn in (1, -1):
        L, _ = chain_dp(x, y, W, sgn)
        if L.size and L.max() > m:
            m = L.max()
    return m


def longest_chain(x, y, W=CHAIN_W):
    """Indices (into x, y; x sorted ascending) of the longest collinear chain, in either orientation."""
    best = None
    for sgn in (1, -1):
        L, P = chain_dp(np.asarray(x, np.int64), np.asarray(y, np.int64), W, sgn)
        if L.size and (best is None or L.max() > best[0].max()):
            best = (L, P)
    if best is None:
        return np.empty(0, np.int64)
    L, P = best
    i, out = int(np.argmax(L)), []
    while i >= 0:
        out.append(i)
        i = int(P[i])
    return np.array(out[::-1], np.int64)


def chain_blocks(px, py, minlen, *, W) -> list[tuple[np.ndarray, int]]:
    """Greedy collinear blocks of anchors at positions px along one sequence and py along the other (LTRmix's): the
    anchors ranked on each (stable), then repeatedly the longest chain (chain_dp, steps <= W ranks) among the anchors
    left, in either orientation (+1 when the two tie; within one, the first point of the longest), taken out while it
    holds >= minlen anchors. Each block: (indices into px in position order, int64; orientation sign)."""
    n = len(px)
    if n < minlen:
        return []
    ox = np.argsort(px, kind="stable")
    rx = np.empty(n, np.int64)
    rx[ox] = np.arange(n)
    ry = np.empty(n, np.int64)
    ry[np.argsort(py, kind="stable")] = np.arange(n)
    alive = np.ones(n, bool)
    out = []
    while alive.sum() >= minlen:
        ids = ox[alive[ox]]
        best = None
        for sign in (1, -1):
            L, P = chain_dp(rx[ids], ry[ids], W, sign)
            i = int(np.argmax(L))
            if best is None or L[i] > best[0]:
                best = (int(L[i]), sign, P, i)
        length, sign, P, i = best
        if length < minlen:
            break
        block = []
        while i >= 0:
            block.append(ids[i])
            i = int(P[i])
        block = np.array(block[::-1], np.int64)
        alive[block] = False
        out.append((block, sign))
    return out


def shuffled_chain(pos, rng, *, n, W) -> int:
    """LTRmix's null length of a collinear chain among anchors at positions `pos` along one sequence: the longest
    chain (chain_length, steps <= W ranks) when the partners' ranks are shuffled, the maximum over n draws; each draw
    is one rng.permutation(len(pos)), assigned to the anchors in input order and read in position order (stable)."""
    o = np.argsort(pos, kind="stable")
    x = np.arange(len(pos), dtype=np.int64)
    return max((int(chain_length(x, rng.permutation(len(pos))[o], W)) for _ in range(n)), default=0)


def chain_null(x, y, rng, *, n=NULL_PERMS, W=CHAIN_W) -> np.ndarray:
    """Longest-chain lengths (floats) of n random partner orders: one rng.permutation(y) per replicate, in turn."""
    return np.array([chain_length(x, rng.permutation(y), W) for _ in range(n)], float)


def chain_z(length, null) -> float:
    """z of a chain length against its permutation null, the null's spread floored at half an anchor."""
    return (length - null.mean()) / max(null.std(), CHAIN_SD_FLOOR)


def chain_significant(n, length, z, n_tests, W=CHAIN_W):
    """Chain of `length` among n candidates significant at a Bonferroni level of 0.05 over n_tests sequence pairs.
    When n <= W the window does not constrain chains (any monotone subsequence counts) and the permutation null is too
    discrete for a normal approximation: the exact bound P(chain >= l) <= 2 C(n, l) / l! (expected monotone
    l-subsequences, two orientations) decides. Otherwise z against the permutations decides."""
    from scipy.stats import norm
    alpha = CHAIN_ALPHA / max(n_tests, 1)
    if n <= W:
        lg = math.lgamma
        return math.log(2) + lg(n + 1) - lg(length + 1) - lg(n - length + 1) - lg(length + 1) <= math.log(alpha)
    return z >= float(norm.isf(alpha))


def _partitions(n, largest=None):
    largest = n if largest is None else largest
    if n == 0:
        yield ()
        return
    for k in range(min(n, largest), 0, -1):
        for rest in _partitions(n - k, k):
            yield (k,) + rest


@cache
def lis_tail(n):
    """P(longest increasing subsequence >= l) of a uniformly random permutation of n, for l = 0..n+1, exact: the
    permutations with LIS = l number the sum, over partitions lambda of n with lambda_1 = l, of (f^lambda)^2 (RSK),
    f^lambda by the hook-length formula."""
    cnt = [0] * (n + 2)
    for lam in _partitions(n):
        conj = [sum(1 for r in lam if r > c) for c in range(lam[0])] if lam else []
        hooks = 1
        for i, r in enumerate(lam):
            for j in range(r):
                hooks *= r - j + conj[j] - i - 1
        cnt[lam[0] if lam else 0] += (math.factorial(n) // hooks) ** 2
    tail, acc = [0.0] * (n + 2), 0
    for l in range(n, -1, -1):
        acc += cnt[l]
        tail[l] = acc / math.factorial(n)
    return tuple(tail)


def chain_test(x, y, length, alpha, rng, W=CHAIN_W):
    """p of a collinear chain at least `length` long among n anchors (x = 0..n-1 along one chromosome, y their ranks
    along the other). n <= W: the window does not bind, so the chain is the longest monotone subsequence, whose tail is
    exact for one orientation and doubled for the two (a union bound). Larger n: B = ceil(2 / alpha) permutations of y,
    p = (longer-or-equal chains + 1) / (permutations + 1), stopped once p > alpha is certain (a normal approximation of
    a 50-permutation z runs 3-7x too liberal on simulated non-homoeologs with 26-60 anchors)."""
    n = x.size
    if n <= W:
        return min(1.0, 2.0 * lis_tail(n)[min(length, n + 1)])
    B, k = math.ceil(CHAIN_PERMS_PER_ALPHA / alpha), 0
    for b in range(1, B + 1):
        k += chain_length(x, rng.permutation(y), W) >= length
        if (k + 1) / (B + 1) > alpha:
            return (k + 1) / (b + 1)
    return (k + 1) / (B + 1)


def anchor_ranks(pos_a, pos_b):
    """Anchors ordered along the first sequence: the order (stable by position pos_a), x = 0..n-1 and y = each
    ordered anchor's rank by its partner's position pos_b (stable), int64."""
    o = np.argsort(pos_a, kind="stable")
    x = np.arange(o.size, dtype=np.int64)
    y = np.argsort(np.argsort(np.asarray(pos_b)[o], kind="stable"), kind="stable").astype(np.int64)
    return o, x, y
