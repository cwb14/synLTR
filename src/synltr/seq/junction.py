"""Junction reads at an insertion site: the exact placement of a window gap (`refine_gap`), and the target-site
duplication checks (`tsd_allowance`, `hamming`). Coordinates are 0-based half-open window coordinates, as in `Gap`.

The local aligner can leave a gap's edges tens of bp off. `refine_gap` re-aligns the two junction
regions base by base, so an empty insertion site (A gap = element + one TSD copy, B contiguous) can be
told from a deletion that also took flanking sequence.

Model, on windows Al (A around a0), Ar (A around a1) and Bj (B around [b0, b1)):
    left part   Al[:i] aligned to Bj[:j1], global at both starts (Al[0], Bj[0] are aligned flank)
    right part  Ar[k:] aligned to Bj[j2:], global at both ends
    the A between the parts is deleted at no cost (it is the element); Bj[j1:j2] is B-only sequence
    and costs one gap (nothing if j1 == j2)
Each part is an affine-gap (Gotoh) DP scored with constants.SCORE: match +1, mismatch -4, a gap of n
bases -(6 + 2n); N, or any non-identical pair, is a mismatch. The split (i, j1, k, j2) maximizes the
summed score. Ties: fewest gaps opened inside the parts, then smallest j2 - j1, then smallest i (the
leftmost placement), then smallest k. The first rule keeps B-only sequence at the junction in
[b0, b1): it also fits, at the same score, as a gap at the end of a part or beside A's other TSD copy.
A pure deletion is then left-normalized on the sequence (past the windows if need be) and given its
slide. No refinement runs across an N in the B window: the free deletion would swallow it.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numba import njit

from synltr.base.constants import REFINE_PAD, SCORE


@dataclass
class Gap:
    """A-unaligned [a0, a1) facing B-unaligned [b0, b1); window coordinates.

    `slide` (set by refine_gap for a pure A deletion, b0 == b1): the deletion scores
    the same at (a0 + t, a1 + t, b0 + t, b1 + t) for every t <= slide (a direct repeat at its ends,
    e.g. a TSD); the gap itself is the leftmost of these placements on the whole sequence (left-
    normalized: the base before a0 differs from the base before a1).
    """

    a0: int
    a1: int
    b0: int
    b1: int
    slide: int = 0


_N = ord("N")
_NEG = -(1 << 62)  # minus infinity for the int64 DP scores


def refine_gap(a_seq: str, b_seq: str, gap: Gap, pad: int = REFINE_PAD) -> Gap:
    """The exact gap near `gap` (A [a0, a1) facing B [b0, b1) of the windows a_seq, b_seq).

    Windows, clipped to the sequences: Bj = B[b0 - pad, b1 + pad); Al = A[a0 - pad, a0 + pad) and
    Ar = A[a1 - pad, a1 + pad), or, if a1 - a0 < 2 * pad, Al = A[a0 - pad, m) and Ar = A[m, a1 + pad)
    with m the gap midpoint, so no A base falls between them.

    A pure deletion (b0 == b1) is left-normalized: moved left one base at a time (all four
    coordinates), past the windows if need be, while the base before it equals its last base. It
    comes back at its leftmost placement with `slide` the largest k such that A[a0 + t] == A[a1 + t]
    for all t < k (up to the sequence end): each placement a0 .. a0 + slide scores the same.

    Returns `gap` itself (1) if B [b0, b1) is longer than 2 * pad: a large B insertion needs no exact
    boundary; (2) if Bj holds an N (any case): an assembly gap near the junction would be swallowed by
    the free A deletion as if it were deleted flank. Callers check for N first; this is a second guard.
    """
    if pad < 1:
        raise ValueError(f"pad must be >= 1, got {pad}")
    a0, a1, b0, b1 = gap.a0, gap.a1, gap.b0, gap.b1
    if not (0 <= a0 <= a1 <= len(a_seq) and 0 <= b0 <= b1 <= len(b_seq)):
        raise ValueError(f"gap A [{a0}, {a1}), B [{b0}, {b1}) does not fit the windows "
                         f"(A {len(a_seq)} bp, B {len(b_seq)} bp)")
    if b1 - b0 > 2 * pad:
        return gap
    s0 = max(0, b0 - pad)
    bj = b_seq[s0:min(len(b_seq), b1 + pad)]
    if "N" in bj or "n" in bj:
        return gap
    if a1 - a0 < 2 * pad:  # short gap: the A windows meet at its midpoint
        l1 = r0 = a0 + (a1 - a0) // 2
    else:
        l1, r0 = a0 + pad, a1 - pad
    l0, r1 = max(0, a0 - pad), min(len(a_seq), a1 + pad)
    i, j1, k, j2 = _junction(a_seq[l0:l1], a_seq[r0:r1], bj)
    x0, x1, y0, y1 = l0 + i, r0 + k, s0 + j1, s0 + j2
    if y0 != y1 or x1 == x0:  # B-only sequence, or no indel at all
        return Gap(x0, x1, y0, y1)
    # the DP gives the leftmost placement in the windows; a repeat can reach past them
    while x0 > 0 and y0 > 0 and a_seq[x0 - 1].upper() == a_seq[x1 - 1].upper():
        x0, x1, y0, y1 = x0 - 1, x1 - 1, y0 - 1, y1 - 1
    return Gap(x0, x1, y0, y1, _slide(a_seq, x0, x1))


def _junction(al: str, ar: str, bj: str) -> tuple[int, int, int, int]:
    """Best split (i, j1, k, j2) of the windows al, ar, bj under the module's model."""
    match, mismatch, gap_open, gap_extend = SCORE
    # Scores are scaled by w and each gap opened inside a part costs one unit more. The two parts open at
    # most w - 1 gaps (each holds a base of al, ar or bj; bj is aligned twice), so the unit only decides
    # between equal scores: the first tie rule.
    w = len(al) + len(ar) + 2 * len(bj) + 1
    sc = (match * w, mismatch * w, gap_open * w + 1, gap_extend * w)
    fbest, fi = _part(_codes(al), _codes(bj), *sc, False)
    rbest, rk = _part(_codes(ar[::-1]), _codes(bj[::-1]), *sc, True)  # the right part, reversed
    split = _split(fbest, fi, np.ascontiguousarray(rbest[::-1]), len(ar) - rk[::-1],
                   gap_open * w, gap_extend * w)
    return int(split[0]), int(split[1]), int(split[2]), int(split[3])


def _codes(s: str) -> np.ndarray:
    return np.frombuffer(s.upper().encode("ascii"), dtype=np.uint8)


@njit(cache=True)
def _part(x, y, match, mismatch, gap_open, gap_extend, prefer_last):
    """Gotoh alignment of x against y, global at both starts.

    Returns best[j] = max over i of H[i][j], the best score of x[:i] against y[:j] (any final state),
    and arg[j] = that i: the smallest on ties, the largest if prefer_last. A gap of n bases costs
    gap_open + gap_extend * n; N or a non-identical pair is a mismatch.
    """
    n, m = x.shape[0], y.shape[0]
    oe = gap_open + gap_extend
    h = np.empty(m + 1, np.int64)  # H[i-1][.] entering row i, H[i][.] after it
    d = np.full(m + 1, _NEG, np.int64)  # D[i][.]: ends with x[i-1] against a gap
    h[0] = 0
    for j in range(1, m + 1):
        h[j] = -(gap_open + gap_extend * j)
    best = h.copy()
    arg = np.zeros(m + 1, np.int64)
    for i in range(1, n + 1):
        xi = x[i - 1]
        diag = h[0]
        h[0] = -(gap_open + gap_extend * i)
        ins = _NEG  # I[i][j]: ends with y[j-1] against a gap
        for j in range(1, m + 1):
            dj = max(h[j] - oe, d[j] - gap_extend)
            ins = max(h[j - 1] - oe, ins - gap_extend)
            hj = max(diag + (match if xi == y[j - 1] and xi != _N else -mismatch), dj, ins)
            diag = h[j]
            d[j] = dj
            h[j] = hj
            if hj > best[j] or (prefer_last and hj == best[j]):
                best[j] = hj
                arg[j] = i
    return best, arg


@njit(cache=True)
def _split(fbest, fi, rbest, rk, gap_open, gap_extend):
    """(i, j1, k, j2) maximizing fbest[j1] + rbest[j2] - (0 if j1 == j2 else gap_open + gap_extend *
    (j2 - j1)) over j1 <= j2, with i = fi[j1] and k = rk[j2]; ties: smallest (j2 - j1, i, k, j1).

    For each j2 only two j1 can win: j2 itself, or the j1 < j2 with the largest
    fbest[j1] + gap_extend * j1 (the largest such j1 on ties), kept as a running maximum.
    """
    top = _NEG
    key = (0, 0, 0, 0, 0)  # (j2 - j1, i, k, j1, j2) of the best split so far
    run, rj = _NEG, -1
    for j2 in range(fbest.shape[0]):
        if j2 > 0 and fbest[j2 - 1] + gap_extend * (j2 - 1) >= run:
            run, rj = fbest[j2 - 1] + gap_extend * (j2 - 1), j2 - 1
        for j1 in (j2, rj):
            if j1 < 0:
                continue
            tot = fbest[j1] + rbest[j2] - (0 if j1 == j2 else gap_open + gap_extend * (j2 - j1))
            cand = (j2 - j1, fi[j1], rk[j2], j1, j2)
            if tot > top or (tot == top and cand < key):
                top, key = tot, cand
    return key[1], key[3], key[2], key[4]


def _slide(a: str, a0: int, a1: int) -> int:
    """Largest k with a[a0 + t] == a[a1 + t] (ignoring case) for all t < k, up to the end of a; 0 if
    a1 == a0 (no deletion)."""
    if a1 <= a0:
        return 0
    k, end, step = 0, len(a) - a1, 64
    while k < end:
        x, y = a[a0 + k:a0 + k + step].upper(), a[a1 + k:a1 + k + step].upper()
        if x[:len(y)] != y:
            return k + next(t for t, (p, q) in enumerate(zip(x, y)) if p != q)
        k += len(y)
        step *= 2
    return end


def tsd_allowance(L: int, k: float | None, ident: float) -> int:
    """TSD mismatches an empty site may show: max(1, floor(L * (K + d) + 0.5)).

    K = the element's LTR-LTR divergence (its age: TSD copies decay; 0 when unknown), d = 1 - ident (the
    A-B divergence adds its own mismatches; 0 when ident is nan).
    """
    k = 0.0 if k is None or math.isnan(k) else k
    d = 0.0 if math.isnan(ident) else 1.0 - ident
    return max(1, math.floor(L * (k + d) + 0.5))


def hamming(x: str, y: str) -> int:
    """Mismatches between equal-length strings, ignoring case; N never matches."""
    return sum(p != q or p == "N" for p, q in zip(x.upper(), y.upper()))
