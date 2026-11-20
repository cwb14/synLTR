"""CIGAR strings and their walks.

A CIGAR is held as two arrays, the op codes and their lengths. A code indexes OPS ("MIDNSHP=X", SAM's codes 0-8, as
pysam and mappy report them). `parse` reads two dialects: "sam" (SAM and PAF cg:Z: CIGARs; every op has a count)
and "lastz" (LASTZ's cigarx, which writes a count of 1 without digits).

Walks are 0-based: `walk` gives the query and reference positions before each op and after the last, `blocks` the
aligned (M, = or X) ops as (query start, reference start, length), each the half-open interval [start, start +
length) on its sequence. I and S ops advance the query, D and N the reference, H and P neither.

Read alignments (pysam's `cigartuples`: (op, length) pairs, the alignment starting at 0-based reference position
r0) have three queries: `aligned_in` (aligned bases inside a reference interval), `qpos_at` (the read position at a
reference position, hard-clipped bases counted, as in the read before clipping) and `deletion_at` (the deletion
holding a reference position). They walk the ops themselves rather than `blocks`: blocks drop hard clips from
query positions and the order of the insertions and deletions between two aligned ops, which `qpos_at` (inside a
deletion) and `deletion_at` (each deletion op on its own) keep.
"""
from __future__ import annotations

import re

import numpy as np

OPS = "MIDNSHP=X"
ALIGNMENT_OPS = (0, 1, 2, 3, 7, 8)  # M, I, D, N, =, X: the ops of an alignment of two sequences, without clips
_QUERY = np.array([1, 1, 0, 0, 1, 0, 0, 1, 1], dtype=np.int64)  # M I S = X advance the query
_REFERENCE = np.array([1, 0, 1, 1, 0, 0, 0, 1, 1], dtype=np.int64)  # M D N = X advance the reference
_DIALECTS = {  # dialect -> (whole-string pattern, one-op pattern, what an error calls it)
    "sam": (re.compile(r"(?:\d+[MIDNSHP=X])+"), re.compile(r"(\d+)([MIDNSHP=X])"), "CIGAR"),
    "lastz": (re.compile(r"(?:\d*[MIDX=])*"), re.compile(r"(\d*)([MIDX=])"), "LASTZ cigar"),
}


def parse(text: str, dialect: str) -> tuple[np.ndarray, np.ndarray]:
    """The op codes (uint8, indexing OPS) and lengths (int64) of CIGAR `text` in `dialect` ("sam" or "lastz").

    A string that is not a CIGAR of the dialect is a ValueError ("malformed CIGAR ..." or "malformed LASTZ cigar
    ..."); an empty string is one only in "sam".
    """
    try:
        whole, one, what = _DIALECTS[dialect]
    except KeyError:
        raise ValueError(f"unknown CIGAR dialect {dialect!r}; known: {', '.join(_DIALECTS)}") from None
    if not whole.fullmatch(text):
        raise ValueError(f"malformed {what} {text[:40]!r}")
    pairs = one.findall(text)
    ops = np.fromiter((OPS.index(c) for _, c in pairs), dtype=np.uint8, count=len(pairs))
    lens = np.fromiter((int(n) if n else 1 for n, _ in pairs), dtype=np.int64, count=len(pairs))
    return ops, lens


def walk(ops, lens, q0: int = 0, r0: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """The query and reference positions before each op and after the last (len(ops) + 1 each), from (q0, r0)."""
    ops = np.asarray(ops, dtype=np.intp)
    lens = np.asarray(lens, dtype=np.int64)
    if len(ops) and not 0 <= ops.min() <= ops.max() < len(OPS):
        raise ValueError(f"CIGAR op code outside 0-{len(OPS) - 1}")
    q = np.concatenate(([q0], q0 + np.cumsum(lens * _QUERY[ops])))
    r = np.concatenate(([r0], r0 + np.cumsum(lens * _REFERENCE[ops])))
    return q, r


def blocks(ops, lens, q0: int, r0: int) -> list[tuple[int, int, int]]:
    """(query start, reference start, length) of each aligned (M, = or X) op, in op order, walking from (q0, r0).

    `ops` and `lens` are arrays or sequences of ints; an op code outside OPS is a ValueError.
    """
    out = []
    q, r = q0, r0
    for op, n in zip(_ints(ops), _ints(lens), strict=True):
        if op == 0 or op == 7 or op == 8:
            out.append((q, r, n))
            q += n
            r += n
        elif op == 1 or op == 4:
            q += n
        elif op == 2 or op == 3:
            r += n
        elif op != 5 and op != 6:
            raise ValueError(f"CIGAR op code {op} is not one of {OPS}")
    return out


def aligned_in(cigar, r0: int, lo: int, hi: int) -> int:
    """Reference bases in [lo, hi) covered by aligned (M, = or X) ops of the alignment `cigar` ((op, length)
    pairs) that starts at reference position `r0`."""
    r, n = r0, 0
    for op, ln in cigar:
        if op == 0 or op == 7 or op == 8:
            n += max(0, min(hi, r + ln) - max(lo, r))
            r += ln
        elif op == 2 or op == 3:
            r += ln
        if r >= hi:
            break
    return n


def qpos_at(cigar, r0: int, refpos: int) -> int | None:
    """The query position (soft and hard clips counted: a position in the read as sequenced, in the stored
    orientation) aligned to reference position `refpos`, for the alignment `cigar` ((op, length) pairs) starting at
    reference position `r0`; inside a deletion or skip, the query position at that op. None outside the alignment."""
    if refpos < r0:
        return None
    q, r = 0, r0
    for op, ln in cigar:
        if op == 0 or op == 7 or op == 8:
            if r + ln > refpos:
                return q + (refpos - r)
            q += ln
            r += ln
        elif op == 1 or op == 4 or op == 5:
            q += ln
        elif op == 2 or op == 3:
            if r + ln > refpos:
                return q
            r += ln
    return None


def deletion_at(cigar, r0: int, refpos: int, min_len: int) -> tuple[int, int] | None:
    """(start, end), 0-based half-open on the reference, of the deletion or skip op of at least `min_len` bp that
    holds reference position `refpos` in the alignment `cigar` ((op, length) pairs) starting at `r0`; None when
    there is none (each op counts alone: two deletions split by an insertion are two)."""
    r = r0
    for op, ln in cigar:
        if op == 0 or op == 7 or op == 8:
            r += ln
        elif op == 2 or op == 3:
            if r <= refpos < r + ln and ln >= min_len:
                return (r, r + ln)
            r += ln
        if r > refpos:
            return None
    return None


def _ints(values) -> list[int]:
    """`values` as a list of Python ints (numpy arrays through tolist, which is fast)."""
    return values.tolist() if isinstance(values, np.ndarray) else [int(v) for v in values]
