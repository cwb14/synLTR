"""Window alignment: A (the window holding annotated elements) against its orthologous B window, with LASTZ.

`align_window` aligns A (query) to B (reference, already on A's strand) and returns an `AlnMap`: the collinear,
gap-free aligned segments, plus any inverted hits. All coordinates are 0-based half-open window coordinates.
LASTZ scores every window as minimap2 asm20 (SCORE), whatever the pair (lastz_args says why).
"""
from __future__ import annotations

import logging
import math
import os
import re
import shutil
import subprocess
import tempfile
from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path
from typing import NamedTuple

import numpy as np

from synltr.base.checkpoint import tool_version
from synltr.base.constants import LASTZ_SEED, MIN_SEG, SCORE
from synltr.seq.cigar import ALIGNMENT_OPS, blocks, parse, walk
from synltr.seq.junction import Gap

log = logging.getLogger(__name__)

_N = ord("N")


@dataclass(eq=False)
class AlnMap:
    """Collinear gap-free segments sorted by A, never overlapping on A or B, plus inverted hits.

    Segment i aligns A[seg_a[i], seg_a[i] + seg_len[i]) to B[seg_b[i], seg_b[i] + seg_len[i]) with
    seg_mm[i] mismatches (a base facing N counts as a mismatch). `inv` holds minus-strand hits as
    (a0, a1, b0, b1), B on its forward strand.
    """

    seg_a: np.ndarray
    seg_b: np.ndarray
    seg_len: np.ndarray
    seg_mm: np.ndarray
    inv: list[tuple[int, int, int, int]] = field(default_factory=list)

    def __post_init__(self) -> None:
        cols = [np.asarray(x, dtype=np.int64).reshape(-1)
                for x in (self.seg_a, self.seg_b, self.seg_len, self.seg_mm)]
        if len({len(c) for c in cols}) != 1:
            raise ValueError("AlnMap segment arrays differ in length")
        self.seg_a, self.seg_b, self.seg_len, self.seg_mm = cols
        self._a_end = self.seg_a + self.seg_len
        self._b_end = self.seg_b + self.seg_len
        if (self.seg_len <= 0).any() or (self.seg_mm < 0).any() or (self.seg_mm > self.seg_len).any():
            raise ValueError("AlnMap segments need length > 0 and 0 <= mismatches <= length")
        if (self.seg_a < 0).any() or (self.seg_b < 0).any():
            raise ValueError("AlnMap segments need non-negative coordinates")
        if (self.seg_a[1:] < self._a_end[:-1]).any() or (self.seg_b[1:] < self._b_end[:-1]).any():
            raise ValueError("AlnMap segments must be sorted and must not overlap on A or B")
        self.inv = [tuple(int(v) for v in x) for x in self.inv]
        self._runs_cache: dict[int, tuple[np.ndarray, ...]] = {}

    @classmethod
    def empty(cls) -> AlnMap:
        z = np.zeros(0, dtype=np.int64)
        return cls(z, z, z, z)

    def _overlapping(self, a0: int, a1: int) -> tuple[int, int]:
        """Index range of the segments overlapping A [a0, a1)."""
        if a1 <= a0:
            return 0, 0
        return (int(np.searchsorted(self._a_end, a0, side="right")),
                int(np.searchsorted(self.seg_a, a1, side="left")))

    def _overlap(self, i0: int, i1: int, a0: int, a1: int) -> np.ndarray:
        return np.minimum(self._a_end[i0:i1], a1) - np.maximum(self.seg_a[i0:i1], a0)

    def aligned_bp(self, a0: int, a1: int) -> int:
        """Number of aligned A bases in [a0, a1)."""
        i0, i1 = self._overlapping(a0, a1)
        return int(self._overlap(i0, i1, a0, a1).sum()) if i0 < i1 else 0

    def identity(self, a0: int, a1: int) -> float:
        """1 - mismatches/aligned over the aligned A bases in [a0, a1); nan if none is aligned.

        A segment's mismatches are prorated by the fraction of it inside [a0, a1).
        """
        i0, i1 = self._overlapping(a0, a1)
        if i0 >= i1:
            return math.nan
        ov = self._overlap(i0, i1, a0, a1)
        mm = (self.seg_mm[i0:i1] * ov / self.seg_len[i0:i1]).sum()
        return float(1.0 - mm / ov.sum())

    def _runs(self, min_seg: int) -> tuple[np.ndarray, ...]:
        """(a0, a1, b0, b1) of the runs holding >= min_seg aligned bases.

        A run is a maximal stretch of segments in which both the A gap and the B gap between
        consecutive segments are < min_seg, so small indels never split it.
        """
        if min_seg in self._runs_cache:
            return self._runs_cache[min_seg]
        if len(self.seg_a) == 0:
            runs = (self.seg_a,) * 4
        else:
            brk = (((self.seg_a[1:] - self._a_end[:-1]) >= min_seg)
                   | ((self.seg_b[1:] - self._b_end[:-1]) >= min_seg))
            first = np.concatenate(([0], np.flatnonzero(brk) + 1))
            last = np.append(first[1:] - 1, len(self.seg_a) - 1)
            keep = np.add.reduceat(self.seg_len, first) >= min_seg
            first, last = first[keep], last[keep]
            runs = (self.seg_a[first], self._a_end[last], self.seg_b[first], self._b_end[last])
        self._runs_cache[min_seg] = runs
        return runs

    def _b_at(self, a: int) -> int:
        """B position of `a` inside a run, or the B end of the segment before an A-only gap."""
        i = int(np.searchsorted(self.seg_a, a, side="right")) - 1
        return int(self.seg_b[i] + min(a - int(self.seg_a[i]), int(self.seg_len[i])))

    def gap_around(self, a: int, min_seg: int = MIN_SEG) -> Gap | None:
        """The unaligned gap holding A position `a`, between flanking runs of >= min_seg aligned bases.

        Returns the Gap from the end of the last such run ending at or before `a` to the start of the
        first such run starting at or after `a`; None if either side has no such run; a zero-length
        Gap at `a` (B side at its aligned position) if `a` lies inside such a run.
        """
        ra0, ra1, rb0, rb1 = self._runs(min_seg)
        k = int(np.searchsorted(ra0, a, side="right")) - 1  # last run starting at or before a
        if k >= 0 and a < ra1[k]:
            b = self._b_at(a)
            return Gap(a, a, b, b)
        left = int(np.searchsorted(ra1, a, side="right")) - 1  # last run ending at or before a
        right = int(np.searchsorted(ra0, a, side="left"))  # first run starting at or after a
        if left < 0 or right >= len(ra0):
            return None
        return Gap(int(ra1[left]), int(ra0[right]), int(rb1[left]), int(rb0[right]))

    def b_gaps_within(self, a0: int, a1: int) -> list[Gap]:
        """Gaps with B-only sequence (B gap > 0) between consecutive segments, starting in A [a0, a1).

        A gap starts where its left segment ends on A; the caller filters the gaps by size.
        """
        ea, eb = self._a_end[:-1], self._b_end[:-1]
        na, nb = self.seg_a[1:], self.seg_b[1:]
        sel = np.flatnonzero((ea >= a0) & (ea < a1) & (nb > eb))
        return [Gap(int(ea[i]), int(na[i]), int(eb[i]), int(nb[i])) for i in sel]


def _chain(hits: list) -> list:
    """Collinear chain of plus-strand hits maximizing total mlen.

    Hit j may precede hit i only if q_st, q_en, r_st and r_en of i all exceed those of j.
    Ties go to the hit that sorts first by (q_st, q_en, r_st, r_en).
    """
    hs = sorted(hits, key=lambda h: (h.q_st, h.q_en, h.r_st, h.r_en))
    if not hs:
        return []
    score = [h.mlen for h in hs]
    prev = [-1] * len(hs)
    for i, hi in enumerate(hs):
        for j in range(i):
            hj = hs[j]
            if (hi.q_st > hj.q_st and hi.q_en > hj.q_en and hi.r_st > hj.r_st and hi.r_en > hj.r_en
                    and score[j] + hi.mlen > score[i]):
                score[i], prev[i] = score[j] + hi.mlen, j
    i = max(range(len(hs)), key=lambda k: (score[k], -k))
    out = []
    while i >= 0:
        out.append(hs[i])
        i = prev[i]
    return out[::-1]


def _segments(chain: list) -> list[list[int]]:
    """Gap-free [a, b, len] segments of the chained hits' CIGARs (their aligned blocks, seq.cigar).

    Each match run is clipped on the diagonal so it starts at or after the end of the previous
    segment on both axes (parts overlapping an earlier hit are dropped); abutting runs merge.
    """
    segs: list[list[int]] = []
    last_a = last_b = 0
    for h in chain:
        bad = [int(op) for op in h.ops if op not in ALIGNMENT_OPS]
        if bad:
            raise ValueError(f"unsupported CIGAR op {bad[0]} in the hit at A {h.q_st}-{h.q_en}")
        for qa, rb, n in blocks(h.ops, h.lens, h.q_st, h.r_st):
            d = max(0, last_a - qa, last_b - rb)
            if d < n:
                sa, sb, ln = qa + d, rb + d, n - d
                if segs and segs[-1][0] + segs[-1][2] == sa and segs[-1][1] + segs[-1][2] == sb:
                    segs[-1][2] += ln
                else:
                    segs.append([sa, sb, ln])
                last_a, last_b = sa + ln, sb + ln
        q, r = walk(h.ops, h.lens, h.q_st, h.r_st)
        qa, rb = int(q[-1]), int(r[-1])
        if (qa, rb) != (h.q_en, h.r_en):
            raise ValueError(f"CIGAR of the hit at A {h.q_st}-{h.q_en}, B {h.r_st}-{h.r_en} "
                             f"ends at A {qa}, B {rb}")
    return segs


def alnmap_from_hits(hits: Iterable, a_seq: str, b_seq: str) -> AlnMap:
    """AlnMap from local alignment hits of A (query) on B (reference).

    Hits need the fields q_st, q_en, r_st, r_en, strand, mlen and the CIGAR as op codes and lengths (ops, lens;
    seq.cigar). Minus-strand hits go to `inv`, distinct and sorted. Plus-strand hits are chained (`_chain`) and
    cut into collinear segments (`_segments`) whose mismatches are counted base by base.
    """
    plus, inv = [], []
    for h in hits:
        if h.strand > 0:
            plus.append(h)
        else:
            inv.append((h.q_st, h.q_en, h.r_st, h.r_en))
    segs = _segments(_chain(plus))
    a = np.frombuffer(a_seq.upper().encode("ascii"), dtype=np.uint8)
    b = np.frombuffer(b_seq.upper().encode("ascii"), dtype=np.uint8)
    if segs and (segs[-1][0] + segs[-1][2] > len(a) or segs[-1][1] + segs[-1][2] > len(b)):
        raise ValueError("alignment hits extend beyond the window sequences")
    mm = []
    for sa, sb, ln in segs:
        x, y = a[sa:sa + ln], b[sb:sb + ln]
        mm.append(np.count_nonzero((x != y) | (x == _N) | (y == _N)))
    cols = np.array(segs, dtype=np.int64).reshape(-1, 3)
    return AlnMap(cols[:, 0], cols[:, 1], cols[:, 2], np.array(mm, dtype=np.int64), sorted(set(inv)))


# ---------------------------------------------------------------- LASTZ

_LASTZ_FIELDS = "zstart1,end1,zstart2+,end2+,strand2,nmatch,cigarx"  # B, then A (forward) coordinates
_SHM = "/dev/shm" if os.access("/dev/shm", os.W_OK) else None  # RAM-backed: two small files per window
_TMP_PREFIX = "synltr_lastz."  # + "<pid>.": sweep_lastz_tmp finds a dead process's temp folders by it
_LIVE: dict = {}  # this process's running lastz child ("proc") and temp folder ("tmp"), for kill_lastz
_NO_OPS = np.zeros(0, dtype=np.uint8)
_NO_LENS = np.zeros(0, dtype=np.int64)


class _Hit(NamedTuple):
    """A LASTZ hit as the fields alnmap_from_hits reads (A = query, B = reference); its CIGAR as op codes and
    lengths (seq.cigar), empty on the minus strand."""

    q_st: int
    q_en: int
    r_st: int
    r_en: int
    strand: int
    mlen: int
    ops: np.ndarray
    lens: np.ndarray


def lastz_args() -> list[str]:
    """LASTZ options for a window, shared by both strands (no strand, files) and by every genome pair.

    Scoring: minimap2 asm20's (SCORE: match 1, mismatch 4, gap open 6, extend 2): an alignment scores above
    zero only above 80% identity (Wicker's rule), so paralogous TE copies 70-80% apart do not align as the
    element. asm10's and asm5's stricter scoring fragments the alignment of a locally divergent window (an
    introgressed tract in a near-identical pair) into wrong D, O, T and S calls, so every pair is scored as
    asm20. --ydrop: the cost of a MIN_SEG-bp gap, so indels shorter than MIN_SEG are bridged
    inside one alignment and an element-sized indel ends it (the flanks come out as separate hits).
    --gappedthresh MIN_SEG: less aligned sequence counts as none. --step 1: every B position seeds (windows
    are small). --chain: only the best collinear chain of seed hits (HSPs) is extended, with no diagonal
    penalties (they would fight the element-sized shifts a window exists to read); unchained, every hit
    between two satellite arrays is extended (a centromeric window then takes minutes, not seconds).
    --ambiguous=iupac: IUPAC codes (N included) are ambiguous bases, never a rejected file.
    """
    match, mismatch, gap_open, gap_extend = SCORE
    return [f"--match={match},{mismatch}", f"--gap={gap_open},{gap_extend}",
            f"--ydrop={gap_open + gap_extend * MIN_SEG}", f"--gappedthresh={MIN_SEG}", "--step=1",
            f"--seed={LASTZ_SEED}", "--transition", "--chain", "--gapped", "--ambiguous=iupac",
            f"--format=general-:{_LASTZ_FIELDS}"]


def lastz_version() -> str | None:
    """The lastz binary's version (checkpoint.tool_version); the genotype stage keys hold it, so an aligner
    upgrade recomputes the calls."""
    return tool_version("lastz")


def _lastz_hits(exe: str, target: str, query: str, args: list[str], strand: str) -> list[_Hit]:
    """LASTZ of `query` (A) on `target` (B) on one strand (chained, lastz_args): plus-strand hits with =/X
    CIGARs, minus-strand hits with coordinates only (they go to `inv`)."""
    cmd = [exe, target, query, f"--strand={strand}", *args]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.DEVNULL,
                            text=True)
    _LIVE["proc"] = proc
    try:
        out, err = proc.communicate()
    except BaseException:
        proc.kill()
        proc.wait()
        raise
    finally:
        _LIVE.pop("proc", None)
    if proc.returncode != 0:
        raise RuntimeError(f"lastz exit {proc.returncode}: {err.strip()[-300:]}")
    hits = []
    for line in out.splitlines():
        if not line or line.startswith("#"):
            continue
        b0, b1, a0, a1, a_strand, nmatch, cigar = line.split("\t")
        plus = a_strand == "+"
        ops, lens = parse(cigar, "lastz") if plus else (_NO_OPS, _NO_LENS)
        hits.append(_Hit(int(a0), int(a1), int(b0), int(b1), 1 if plus else -1, int(nmatch), ops, lens))
    return hits


def align_window(a_seq: str, b_seq: str) -> AlnMap:
    """LASTZ on the window: B (target) and A (query) as temp FASTA files, run once per strand (lastz_args).

    Each strand is chained (--chain), so a repeat-dense window hands alnmap_from_hits one collinear chain of
    hits per strand, and satellite arrays do not stall a window. The minus strand runs on its own, so its
    chain's hits (the inverted copies) reach `inv`. A window LASTZ fails on (an error exit, a malformed CIGAR)
    is left unaligned (an empty AlnMap: its elements are U) with one warning: one pathological window does not
    stop a genome pair.
    """
    args = lastz_args()
    if not a_seq or not b_seq:
        return AlnMap.empty()
    exe = shutil.which("lastz")
    if exe is None:
        raise RuntimeError("the window aligner needs lastz on PATH (not found)")
    with tempfile.TemporaryDirectory(prefix=f"{_TMP_PREFIX}{os.getpid()}.", dir=_SHM) as tmp:
        _LIVE["tmp"] = tmp
        try:
            target, query = os.path.join(tmp, "b.fa"), os.path.join(tmp, "a.fa")
            with open(target, "w") as fh:
                fh.write(f">b\n{b_seq.upper()}\n")
            with open(query, "w") as fh:
                fh.write(f">a\n{a_seq.upper()}\n")
            hits = (_lastz_hits(exe, target, query, args, "plus")
                    + _lastz_hits(exe, target, query, args, "minus"))
        except (RuntimeError, ValueError) as exc:
            log.warning("lastz: window A %d bp, B %d bp left unaligned: %s", len(a_seq), len(b_seq), exc)
            return AlnMap.empty()
        finally:
            _LIVE.pop("tmp", None)
    return alnmap_from_hits(hits, a_seq, b_seq)


def kill_lastz() -> None:
    """Kill this process's running lastz child, if any, and remove its temp folder: the cleanup of a dying
    genotype worker (base.runtime.pool_map's `on_exit`). A worker killed outright (SIGKILL) leaves its folder to
    sweep_lastz_tmp; its lastz child then dies at its first write to the closed pipe, or when it finishes."""
    proc = _LIVE.get("proc")
    if proc is not None and proc.poll() is None:
        proc.kill()
    tmp = _LIVE.get("tmp")
    if tmp is not None:
        shutil.rmtree(tmp, ignore_errors=True)


def sweep_lastz_tmp() -> int:
    """Remove the LASTZ temp folders of this user's dead processes (workers killed outright): the folders
    <temp dir>/synltr_lastz.<pid>.* whose pid no longer runs. Returns how many it removed."""
    root, uid, n = Path(_SHM or tempfile.gettempdir()), os.getuid(), 0
    for d in root.glob(f"{_TMP_PREFIX}*"):
        m = _TMP_NAME.fullmatch(d.name)
        try:
            if m is None or not d.is_dir() or d.stat().st_uid != uid or pid_alive(int(m[1])):
                continue
        except OSError:
            continue
        shutil.rmtree(d, ignore_errors=True)
        n += 1
    return n


_TMP_NAME = re.compile(re.escape(_TMP_PREFIX) + r"(\d+)\..+")


def pid_alive(pid: int) -> bool:
    """True if a process `pid` runs (another user's included)."""
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True
