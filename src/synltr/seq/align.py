"""Global alignment of two LTR-RT sequences (uint8 codes: A C G T = 0-3, any other base 4), and the alignment counts
of library copy pairs.

`global_counts`: an edlib unit-cost path fixes a band, an affine-gap realignment inside it (match, mismatch, gap open,
gap extend = `scores`) counts aligned sites (columns with both bases A, C, G or T), transitions, transversions, gap
events (runs of gap columns in one sequence, leading and trailing gaps included) and gap bases; unit costs undercount
substitutions. Positions are 0-based.

`unit_cost_path`, `position_map`: the edlib unit-cost global path as runs, and the position each base of the first
sequence takes in the second (the library screen's aligner, which counts the columns of this path; unit costs
undercount substitutions).

`pairs`: those counts for copy pairs of a library (`codes`, `offs`, `orient`, `l5` and `l3` of io.library's Library:
LTR record offsets 0-based half-open, 5' LTR first, -1 = none), by mode: "ltr", the 5'-5' plus the 3'-3' LTR
alignments (whole records when either copy has no LTR bounds); "whole", the whole records; "self", each copy's own 5'
against its 3' LTR. Stranded records share element orientation; when either strand is unknown (orient 0), the second
copy's other strand is used when it is closer (`orient`). One fork pool over chunks of ALIGN_CHUNK[mode] items (at
most ALIGN_JOBS per thread), every chunk in a worker; the kernels compile before the fork, on ALIGN_WARM bases of the
library as a slice and as a reverse-complement view (the two array layouts the workers pass), so no worker compiles.
"""
from __future__ import annotations

import re
from typing import NamedTuple

import edlib
import numpy as np
from numba import njit

from synltr.base.constants import ALIGN_BAND, ALIGN_CHUNK, ALIGN_JOBS, ALIGN_WARM, SCORES
from synltr.base.runtime import pool_map, shared

ALPH = np.frombuffer(b"ACGTN", np.uint8)
RCODE = np.array([3, 2, 1, 0, 4], np.uint8)  # complement of each code
_CIGAR = re.compile(rb"(\d+)([=XID])")
_OPCODE = {ord("="): 0, ord("X"): 0, ord("I"): 1, ord("D"): 2}
MODES = ("ltr", "whole", "self")


class Counts(NamedTuple):
    """Counts of one global alignment (ints), or of many (int64 arrays, one row per alignment): aligned sites with both
    bases A, C, G or T, transitions, transversions, gap events (runs of gap columns in one sequence) and gap bases."""

    sites: int | np.ndarray
    ts: int | np.ndarray
    tv: int | np.ndarray
    gap_events: int | np.ndarray
    gap_bases: int | np.ndarray


@njit(cache=True)
def _band(n, m, lens, ops, B):
    """Per query row i (0..n): target window [lo, hi] spanning every diagonal (j - i) from the lowest to the highest
    among 0, m - n and those the edlib path visits, +- B. A long indel's single-gap path runs between diagonals 0 and
    m - n, which edlib's unit-cost path (many short gaps along a drifting diagonal) can leave by far more than B."""
    i = 0
    j = 0
    dmin = min(0, m - n)
    dmax = max(0, m - n)
    for k in range(lens.size):
        for _ in range(lens[k]):
            if ops[k] == 0:
                i += 1
                j += 1
            elif ops[k] == 1:
                i += 1
            else:
                j += 1
            dmin = min(dmin, j - i)
            dmax = max(dmax, j - i)
    lo = np.empty(n + 1, np.int64)
    hi = np.empty(n + 1, np.int64)
    for r in range(n + 1):
        lo[r] = max(0, r + dmin - B)
        hi[r] = min(m, r + dmax + B)
    lo[0] = 0
    hi[n] = m
    return lo, hi


@njit(cache=True)
def _gotoh_band(q, t, lo, hi, sm, sx, go, ge):
    """Global affine-gap alignment restricted to the band (lo, hi nondecreasing) -> (sites with both bases ACGT,
    transitions, transversions, touched: the best path runs along an inner band edge, so a wider band may score better;
    gap events: runs of gap columns in one sequence; gap bases). Ties prefer the diagonal (substitution), then a gap in
    the query, then in the target. Each row runs in three vectorizable passes: gap-in-target and diagonal scores from
    the previous row, the gap-in-query chain as a running maximum (E[j] = max over k < j of Hd[k] - go - ge (j - k), Hd
    = best without E), then the row's scores. int32 scores (no overflow below ~10^8 bp); two rolling rows; only the
    traceback (1 byte per band cell) is kept."""
    n = q.size
    m = t.size
    NEG = np.int32(-(1 << 29))
    goe = np.int32(go + ge)
    ge = np.int32(ge)
    prof = np.zeros((5, m + 1), np.int32)             # prof[x, j]: score of query base x against t[j - 1]
    for x in range(4):
        for j in range(1, m + 1):
            y = t[j - 1]
            prof[x, j] = 0 if y > 3 else (sm if x == y else sx)
    base = np.zeros(n + 2, np.int64)
    W = 0
    for i in range(n + 1):
        base[i + 1] = base[i] + hi[i] - lo[i] + 1
        W = max(W, hi[i] - lo[i] + 1)
    P = np.empty(base[n + 1], np.uint8)              # bits 0-1: H source (0 diag, 1 E, 2 F); bit 2: E extends; 3: F extends
    Hp = np.full(W + 2, NEG, np.int32)               # rows indexed 1 + (j - lo[i]); padded with NEG on both sides
    Fp = np.full(W + 2, NEG, np.int32)
    Hc = np.full(W + 2, NEG, np.int32)
    Fc = np.full(W + 2, NEG, np.int32)
    D = np.empty(W, np.int32)
    Hd = np.empty(W, np.int32)
    E = np.empty(W, np.int32)
    w0 = hi[0] - lo[0] + 1
    for c in range(w0):                               # row 0: leading gap in the query
        j = lo[0] + c
        if j == 0:
            Hc[1 + c] = 0
            P[base[0] + c] = 0
        else:
            Hc[1 + c] = -(go + ge * j)
            P[base[0] + c] = 1 | (4 if j > 1 else 0)
    wold_c = w0                                       # widths last written into the two row buffers
    wold_p = 0
    for i in range(1, n + 1):
        Hp, Hc = Hc, Hp
        Fp, Fc = Fc, Fp
        wold_c, wold_p = wold_p, wold_c
        l0 = lo[i]
        w = hi[i] - l0 + 1
        sh = l0 - lo[i - 1]
        pr = prof[q[i - 1], l0:l0 + w]
        HpA = Hp[1 + sh:1 + sh + w]                   # previous row at column l0 + c
        HpB = Hp[sh:sh + w]                           # previous row at column l0 + c - 1
        FpA = Fp[1 + sh:1 + sh + w]
        FcV = Fc[1:1 + w]
        HcV = Hc[1:1 + w]
        PV = P[base[i]:base[i] + w]
        Dv = D[:w]
        Hdv = Hd[:w]
        Ev = E[:w]
        for c in range(w):
            fo = HpA[c] - goe
            fe = FpA[c] - ge
            f = max(fe, fo)
            FcV[c] = f
            d = HpB[c] + pr[c]
            Dv[c] = d
            Hdv[c] = max(d, f)
            PV[c] = 8 if fe >= fo else 0
        if l0 == 0:                                   # column 0: no diagonal
            Dv[0] = NEG
            Hdv[0] = FcV[0]
        run = NEG
        Ev[0] = NEG
        for c in range(1, w):
            run = max(run - ge, Hdv[c - 1] - goe)
            Ev[c] = run
        for c in range(w):
            e = Ev[c]
            d = Dv[c]
            f = FcV[c]
            ef = max(e, f)
            HcV[c] = max(d, ef)
            PV[c] = PV[c] | ((d < ef) * (2 - (e >= f)))
        for c in range(1, w):
            PV[c] = PV[c] | (4 * (Ev[c - 1] - ge >= Hdv[c - 1] - goe))
        for c in range(w + 1, wold_c + 1):            # clear what an older, wider row left beyond this one
            Hc[c] = NEG
            Fc[c] = NEG
        wold_c = w
    i = n
    j = m
    state = 0
    sites = 0
    ts = 0
    tv = 0
    touched = False
    gap_events = 0
    gap_bases = 0
    while i > 0 or j > 0:
        if (j == lo[i] and lo[i] > 0) or (j == hi[i] and hi[i] < m):
            touched = True
        pb = P[base[i] + j - lo[i]]
        if state == 0:
            arg = pb & 3
            if arg == 0:
                x = q[i - 1]
                y = t[j - 1]
                if x < 4 and y < 4:
                    sites += 1
                    if x != y:
                        if (x ^ y) == 2:
                            ts += 1
                        else:
                            tv += 1
                i -= 1
                j -= 1
            else:
                state = arg
                gap_events += 1
        elif state == 1:
            j -= 1
            gap_bases += 1
            if (pb & 4) == 0:
                state = 0
        else:
            i -= 1
            gap_bases += 1
            if (pb & 8) == 0:
                state = 0
    return sites, ts, tv, touched, gap_events, gap_bases


def unit_cost_path(q, t) -> np.ndarray:
    """edlib's unit-cost global (NW) path of two code arrays as runs, int64 (2, n): row 0 the run lengths, row 1 the
    operation (0 an aligned column, match or mismatch; 1 a base of q against a gap; 2 a base of t against a gap). An
    empty sequence aligns as gap runs."""
    if q.size == 0 or t.size == 0:
        runs = [(n, op) for n, op in ((q.size, 1), (t.size, 2)) if n]
        return np.array(runs, np.int64).reshape(-1, 2).T.copy()
    r = edlib.align(ALPH[q].tobytes(), ALPH[t].tobytes(), mode="NW", task="path")
    m = _CIGAR.findall(r["cigar"].encode())
    return np.array([[int(x) for x, _ in m], [_OPCODE[o[0]] for _, o in m]], np.int64).reshape(2, -1)


def position_map(path) -> np.ndarray:
    """For each base of q in a unit_cost_path, the 0-based position of the base of t aligned with it (a match or a
    mismatch), -1 where it faces a gap: int64."""
    lens, ops = path
    out = np.full(int(lens[ops != 2].sum()), -1, np.int64)
    i = j = 0
    for n, op in zip(lens.tolist(), ops.tolist()):
        if op == 0:
            out[i:i + n] = np.arange(j, j + n)
            i += n
            j += n
        elif op == 1:
            i += n
        else:
            j += n
    return out


def global_counts(q, t, *, band=ALIGN_BAND, scores=SCORES) -> Counts:
    """Global alignment counts of two code arrays (module docstring): the edlib path fixes a band of +-band columns,
    widened 4-fold while the best affine path rides its edge. An empty sequence aligns as one gap."""
    if q.size == 0 or t.size == 0:
        g = int(q.size + t.size)
        return Counts(0, 0, 0, int(g > 0), g)
    lens, ops = unit_cost_path(q, t)
    B = band
    while True:                                 # widen 4-fold while the best path rides the band edge
        lo, hi = _band(q.size, t.size, lens, ops, B)
        sites, ts, tv, touched, gap_events, gap_bases = _gotoh_band(q, t, lo, hi, *scores)
        if not touched or (lo[-1] == 0 and hi[0] == t.size):
            return Counts(sites, ts, tv, gap_events, gap_bases)
        B *= 4


def _rc(codes):
    """The reverse complement of a code array (a reversed view)."""
    return RCODE[codes][::-1]


def orient(x, y, rev=None, *, mode="NW", head=None, slack=0) -> int:
    """+1 when y aligns to x as given, -1 when y's other strand does: `rev` (default y's reverse complement) at a
    strictly smaller edlib edit distance; a tie keeps y as given. By default (phasing's rule) whole sequences, global
    (NW) distances. LTRmix's rule: mode "HW" (infix: y's ends free), head 3000, slack 200: x's first k bases, k =
    min(len x, len y, head), against the first k + slack of y and of the other strand."""
    other = _rc(y) if rev is None else rev
    if head is not None:
        k = min(x.size, y.size, head)
        x, y, other = x[:k], y[:k + slack], other[:k + slack]
    closer = edlib.align(ALPH[x].tobytes(), ALPH[other].tobytes(), mode=mode)["editDistance"] < \
        edlib.align(ALPH[x].tobytes(), ALPH[y].tobytes(), mode=mode)["editDistance"]
    return -1 if closer else 1


def _record(lib, i):
    return lib.codes[lib.offs[i]:lib.offs[i + 1]]


def _counts(lib, x, mode) -> Counts:
    """The counts of one item: copy x's own LTRs (mode "self") or the copy pair x = (a, b) (module docstring)."""
    if mode == "self":
        r = _record(lib, x)
        return global_counts(r[lib.l5[x, 0]:lib.l5[x, 1]], r[lib.l3[x, 0]:lib.l3[x, 1]])
    a, b = x
    ra, rb = _record(lib, a), _record(lib, b)
    unknown = lib.orient[a] == 0 or lib.orient[b] == 0
    if mode == "ltr" and lib.l5[a, 0] >= 0 and lib.l5[b, 0] >= 0:
        a5, a3 = ra[lib.l5[a, 0]:lib.l5[a, 1]], ra[lib.l3[a, 0]:lib.l3[a, 1]]
        b5, b3 = rb[lib.l5[b, 0]:lib.l5[b, 1]], rb[lib.l3[b, 0]:lib.l3[b, 1]]
        if unknown:
            rev = _rc(b3)                             # b's other strand starts with its 3' LTR, reverse-complemented
            if orient(a5, b5, rev) < 0:
                b5, b3 = rev, _rc(b5)
        s5, s3 = global_counts(a5, b5), global_counts(a3, b3)
        return Counts(*(u + v for u, v in zip(s5, s3)))
    if unknown and orient(ra, rb) < 0:
        rb = _rc(rb)
    return global_counts(ra, rb)


def _job(chunk):
    lib, mode = shared["align_lib"], shared["align_mode"]
    out = np.zeros((len(chunk), 5), np.int64)
    for n, x in enumerate(chunk):
        out[n] = _counts(lib, x, mode)
    return out


def _compile(lib) -> None:
    """Compile the kernels in this process for both array layouts the workers pass them: a record slice and a
    reverse-complement view (module docstring)."""
    x = lib.codes[:ALIGN_WARM]
    global_counts(x, x)
    global_counts(x, _rc(x))


def pairs(lib, pairs, *, mode: str, threads: int) -> Counts:
    """The counts (int64 arrays, one row per item, in input order) of copy pairs of a library (rows a, b of `pairs`),
    or, for mode "self", of copies' own LTRs (`pairs` lists the copies): module docstring. Up to `threads` workers."""
    if mode not in MODES:
        raise ValueError(f"mode must be one of {MODES}, not {mode!r}")
    items = np.asarray(pairs, np.int64)
    items = items.reshape(-1) if mode == "self" else items.reshape(-1, 2)
    if len(items) == 0:
        return Counts(*np.zeros((5, 0), np.int64))
    chunks = np.array_split(items, max(1, min(len(items) // ALIGN_CHUNK[mode], ALIGN_JOBS * max(threads, 1))))
    shared.update(align_lib=lib, align_mode=mode)
    try:
        _compile(lib)
        out = np.concatenate(pool_map(_job, chunks, threads, warm=False))
    finally:
        shared.pop("align_lib", None)
        shared.pop("align_mode", None)
    return Counts(*out.T.copy())
