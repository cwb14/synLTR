"""Chaining and the syntenic 1:1 net over alignment blocks (0-based half-open).

`chain_blocks` links raw blocks from any aligner (paralogous hits included) into collinear chains with a
minimap2-style dynamic program. `net_blocks` then keeps a syntenic 1:1 set, best chain first: every A and
every B position belongs to at most one chain, every kept piece is in place in both genomes (top level; its span
in the other genome inside its parent chain's matching gap; or in a chain that is top level there), and each
kept block carries its margin over the best competing chain.
"""

from __future__ import annotations

import logging
import math
from bisect import bisect_right, insort
from collections import Counter, defaultdict
from operator import itemgetter

import numpy as np
from numba import njit

from synltr.base.constants import (
    ANCHOR_MIN,
    CHAIN_INDEL_COST,
    CHAIN_LOOKBACK,
    CHAIN_MAX_GAP,
    CHAIN_OVERLAP,
    CLAIM_BUCKET_BITS,
    MIN_CHAIN_SCORE,
)
from synltr.routes.assembly.blocks import BLOCK_DTYPE, a_to_b

log = logging.getLogger(__name__)

_END = itemgetter(1)
_A_SPAN, _B_SPAN = (1, 2), (4, 5)  # fields of a piece row's (start, end) on A and on B
# a piece's place in one genome's frame (`net_blocks`): no parent; in the parent's gap on its strand or the
# other; out of place but in a chain with a top-level piece; out of place
_TOP, _SYN, _INV, _CHAIN, _OUT = "top", "syn", "inv", "chain", "out"


# --- chaining ---------------------------------------------------------------------------------------------


def chain_blocks(blocks: np.ndarray) -> np.ndarray:
    """Chain blocks into collinear chains; every block joins exactly one chain.

    Blocks are grouped by (a_chrom, b_chrom, strand) and sorted by (a_start, a_end) = (x0, x1). On B,
    minus-strand blocks use (y0, y1) = (-b_end, -b_start), so collinear means increasing on x and on y.
    DP over each group, f[i] = best score of a chain ending at block i:
        f[i] = matches[i]; for j = i-1 down to i-CHAIN_LOOKBACK, stopping once a_start[i] - a_start[j]
        exceeds CHAIN_MAX_GAP + the group's longest block, j may precede i only if both ends advance by
        >= ANCHOR_MIN on both genomes (x0, x1, y0 and y1; minimap2's rule for anchors) and the gaps
        da = x0[i] - x1[j], db = y0[i] - y1[j] are <= CHAIN_MAX_GAP; then
        f[i] = max(f[i], f[j] + matches[i] - max(0, -min(da, db))
                         - CHAIN_INDEL_COST * |da - db| - log2(1 + max(0, min(da, db)))).
    Overlapping blocks thus chain without counting the overlapped bases twice: the two flank alignments of
    a full-length element in A vs a solo LTR in B (they overlap by the LTR on B), or the chunk-edge pieces
    of one alignment. A copy whose start or end does not advance (a tandem paralog) is not linked.

    Extraction: visit blocks by f descending; backtrack through the predecessors until a block already used
    (or the start). A chain's score is its own, f[end] - f[first] + matches[first]: its blocks and the links
    between them, not the link into a block a better chain already used, so a tandem paralog branching off
    an orthologous chain competes at full strength in the net. At a chain start this equals f[end].

    No chain is dropped here: `net_blocks` applies MIN_CHAIN_SCORE after weak (paralogous) chains have
    counted as competitors. Chain ids run 0..n-1 by decreasing score, ties by the input index of the chain's
    first block in A order. Returns a new array of the blocks with `chain` and `chain_score` set, sorted by
    (chain, a_start, a_end); `margin` is untouched. Blocks with an empty A or B span cannot be mapped and
    are skipped with a warning. Exact duplicates (same coordinates and strand: one alignment reported twice,
    never independent evidence of paralogy) are first collapsed to the copy with the most matches; then a
    block lying inside another block on A and on B, on the same diagonal, is dropped (`_uncontained`).
    """
    ok = (blocks["a_end"] > blocks["a_start"]) & (blocks["b_end"] > blocks["b_start"])
    if not ok.all():
        log.warning("skipped %d block(s) with an empty A or B span", int((~ok).sum()))
    valid = np.flatnonzero(ok)
    if len(valid) == 0:
        return np.zeros(0, dtype=BLOCK_DTYPE)
    keep = _distinct(blocks, valid)
    if len(keep) < len(valid):
        log.debug("dropped %d exact duplicate block(s)", len(valid) - len(keep))
    n_distinct = len(keep)
    keep = _uncontained(blocks, keep)
    if len(keep) < n_distinct:
        log.debug("dropped %d block(s) contained in a same-diagonal block", n_distinct - len(keep))
    sub = blocks[keep]
    order = np.lexsort((sub["a_end"], sub["a_start"], sub["strand"], sub["b_chrom"], sub["a_chrom"]))
    srt = sub[order]
    src = keep[order]  # input index of each sorted block

    a0, a1 = srt["a_start"].astype(np.int64), srt["a_end"].astype(np.int64)
    plus = srt["strand"] > 0
    y0 = np.where(plus, srt["b_start"], -srt["b_end"]).astype(np.int64)
    y1 = np.where(plus, srt["b_end"], -srt["b_start"]).astype(np.int64)
    new_group = ((np.diff(srt["a_chrom"]) != 0) | (np.diff(srt["b_chrom"]) != 0)
                 | (np.diff(srt["strand"]) != 0))
    starts = np.concatenate(([0], np.flatnonzero(new_group) + 1))
    sizes = np.diff(np.append(starts, len(srt)))
    gstart = np.repeat(starts, sizes)
    maxlen = np.repeat(np.maximum.reduceat(a1 - a0, starts), sizes)

    matches = srt["matches"].astype(np.int64)
    f, p = _chain_dp(a0, a1, y0, y1, matches, gstart, maxlen,
                     CHAIN_MAX_GAP, CHAIN_INDEL_COST, ANCHOR_MIN, CHAIN_LOOKBACK)
    label, score, first = _extract_chains(f, p, matches, np.argsort(-f, kind="stable"))

    cid = np.empty(len(score), np.int64)
    cid[np.lexsort((src[first], -score))] = np.arange(len(score))
    srt["chain"] = cid[label]
    srt["chain_score"] = score[label]
    log.debug("chained %d blocks into %d chains", len(srt), len(score))
    return srt[np.argsort(srt["chain"], kind="stable")]


_DUP_KEY = ("a_chrom", "a_start", "a_end", "b_chrom", "b_start", "b_end", "strand")


def _distinct(blocks: np.ndarray, idx: np.ndarray) -> np.ndarray:
    """Indices `idx` (non-empty, input order) without exact duplicates, still in input order.

    Rows equal on _DUP_KEY are one alignment; the one with the most matches stays, the first on ties.
    """
    sub = blocks[idx]
    keys = tuple(sub[k] for k in reversed(_DUP_KEY))  # lexsort: last key is primary
    order = np.lexsort((idx, -sub["matches"].astype(np.int64)) + keys)
    srt = sub[order]
    first = np.ones(len(srt), bool)
    first[1:] = np.logical_or.reduce([srt[k][1:] != srt[k][:-1] for k in _DUP_KEY])
    return np.sort(idx[order[first]])


def _uncontained(blocks: np.ndarray, idx: np.ndarray) -> np.ndarray:
    """Indices `idx` (distinct blocks, input order) minus same-diagonal contained blocks, in input order.

    Within an (a_chrom, b_chrom, strand) group, block Y is dropped when another block X spans it on A and on B
    and their diagonals agree within CHAIN_OVERLAP: b_start - a_start on '+', b_end + a_start on '-'. Such a Y
    is a partial copy of X's alignment (chunked mapping reports a chunk-edge alignment again, left-truncated,
    from the next chunk); it cannot chain with X (its A end does not advance) and would compete with X as a
    paralog. A copy on another diagonal (a tandem repeat) stays. Sweep over A start, testing each block
    against the blocks still overlapping it: O(n log n + n * depth).
    """
    sub = blocks[idx]
    a0, a1 = sub["a_start"].astype(np.int64), sub["a_end"].astype(np.int64)
    b0, b1 = sub["b_start"].astype(np.int64), sub["b_end"].astype(np.int64)
    diag = np.where(sub["strand"] > 0, b0 - a0, b1 + a0)
    # a container sorts before what it contains: A start ascending, then A end and B length descending
    order = np.lexsort((b0 - b1, -a1, a0, sub["strand"], sub["b_chrom"], sub["a_chrom"]))
    srt = sub[order]
    new_group = np.ones(len(srt), np.bool_)
    new_group[1:] = ((np.diff(srt["a_chrom"]) != 0) | (np.diff(srt["b_chrom"]) != 0)
                     | (np.diff(srt["strand"]) != 0))
    drop = _contained(a0[order], a1[order], b0[order], b1[order], diag[order], new_group, CHAIN_OVERLAP)
    return np.sort(idx[order[~drop]]) if drop.any() else idx


@njit(cache=True)
def _contained(a0, a1, b0, b1, diag, new_group, tol):
    """drop[i]: block i lies inside an earlier block of its group on A and on B, diagonals within `tol`.

    Blocks are sorted by group, A start ascending, then A end and B length descending, so any container of
    block i comes before it. The active list holds the group's blocks still overlapping the current A start;
    a block ending at or before it can contain no later block and leaves the list for good. Dropped blocks
    stay active: the rule tests each block against every other.
    """
    n = a0.shape[0]
    drop = np.zeros(n, np.bool_)
    active = np.empty(n, np.int64)
    na = 0
    for i in range(n):
        if new_group[i]:
            na = 0
        k = 0
        for t in range(na):
            j = active[t]
            if a1[j] <= a0[i]:
                continue
            active[k] = j
            k += 1
            if (not drop[i] and a1[i] <= a1[j] and b0[j] <= b0[i] and b1[i] <= b1[j]
                    and abs(diag[i] - diag[j]) <= tol):
                drop[i] = True
        active[k] = i
        na = k + 1
    return drop


@njit(cache=True)
def _chain_dp(a_start, a_end, y0, y1, matches, gstart, maxlen, max_gap, indel_cost, min_step, lookback):
    """Best chain score f[i] ending at each block and its predecessor p[i] (-1: chain start).

    Blocks are sorted by group, then (a_start, a_end); gstart[i] is the index of i's group start, maxlen[i]
    its group's longest A span and min_step the advance required of both ends on both genomes. The
    recurrence is in `chain_blocks`.
    """
    n = a_start.shape[0]
    f = np.empty(n, np.float64)
    p = np.full(n, -1, np.int64)
    for i in range(n):
        f[i] = matches[i]
        for j in range(i - 1, max(gstart[i], i - lookback) - 1, -1):
            if a_start[i] - a_start[j] > max_gap + maxlen[i]:
                break
            if (a_start[i] - a_start[j] < min_step or a_end[i] - a_end[j] < min_step
                    or y0[i] - y0[j] < min_step or y1[i] - y1[j] < min_step):
                continue
            da = a_start[i] - a_end[j]
            db = y0[i] - y1[j]
            if da > max_gap or db > max_gap:
                continue
            m = min(da, db)
            s = f[j] + matches[i] - max(0, -m) - (indel_cost * abs(da - db) + math.log2(1.0 + max(0, m)))
            if s > f[i]:
                f[i] = s
                p[i] = j
    return f, p


@njit(cache=True)
def _extract_chains(f, p, matches, visit):
    """Split blocks into chains: from each unused block in `visit` order, backtrack until a used block.

    Returns the chain label of each block, each chain's own score f[end] - f[first] + matches[first] (see
    `chain_blocks`) and each chain's first block (lowest in A order, the last one reached).
    """
    n = f.shape[0]
    label = np.full(n, -1, np.int64)
    score = np.empty(n, np.float64)
    first = np.empty(n, np.int64)
    nc = 0
    for t in range(n):
        end = visit[t]
        if label[end] >= 0:
            continue
        k = end
        last = end
        while k >= 0 and label[k] < 0:
            label[k] = nc
            last = k
            k = p[k]
        score[nc] = f[end] - f[last] + matches[last]
        first[nc] = last
        nc += 1
    return label, score[:nc], first[:nc]


# --- net --------------------------------------------------------------------------------------------------


def net_blocks(blocks: np.ndarray, label: str = "") -> np.ndarray:
    """Keep a syntenic 1:1 set of chained blocks, best chain first, with each kept block's margin.

    Chains are visited by decreasing chain_score, the blocks of a chain in A order. For each block:
      1. Competition: every piece claimed by an earlier, kept chain that the block's A span or B span
         overlaps by more than CHAIN_OVERLAP bp (alignment-edge slop, not competition) records this
         chain's score.
      2. Trimming: the block's free A parts (its A span minus claimed pieces, including the chain's own
         earlier pieces) are mapped to B with `a_to_b`, claimed B is removed, and the rest is mapped back
         to A by the inverse linear map (interval boundaries, clipped to the free A part). Parts of
         >= ANCHOR_MIN bp on A become pieces: new rows with the block's chain, strand and chain_score, and
         matches/alen scaled by the kept A fraction. A block shorter than ANCHOR_MIN yields no piece; an
         unchallenged block is a single part, kept whole (no mapping round trip).
    Then, before the chain claims anything, each piece passes the synteny test below or is dropped. A dropped
    piece claims nothing, so later chains may still fill its positions; its competition credit stands (it is a
    rival alignment). A chain whose kept pieces carry >= MIN_CHAIN_SCORE matches claims them on A and B;
    otherwise it is dropped, and the competition it recorded stands (a chain whose pieces carry fewer matches
    before the test is dropped untested: same outcome). margin = 1 - (best competitor's chain_score) /
    chain_score, 1.0 if never challenged. Input: `chain_blocks` output (each chain on one A sequence, one B
    sequence and one strand). Returns the pieces sorted by (a_chrom, a_start). Logs, at DEBUG, one summary
    line (the kept pieces by place (below), and the dropped pieces with their A bp summed; a position a
    drop frees may be tried, and dropped, again by a later chain) and the per-frame details.

    Synteny test (UCSC's chain/net hierarchy with netSyntenic and `netFilter -syn`, Kent et al. 2003, PNAS
    100:11484, without its size escapes). A kept chain's A extent runs from its first kept piece's a_start to
    its last kept piece's a_end; its B extent likewise on B. A piece's A-parent is the kept (so better) chain
    with the smallest A extent containing the piece's A span on its A sequence. In A's frame the piece is
    supported if
      (a) it has no A-parent (top level),
      (b) its B span lies on the A-parent's B sequence inside the B gap between the A-parent's two kept pieces
          flanking it on A (syn on the parent's strand, an inversion on the other), or
      (c) its chain has a trimmed piece with no A-parent: top-level support counts at the chain level.
    B's frame is the mirror (B-parent, the A span inside the B-parent's A gap). A piece is kept only with
    support in both frames. A paralogous fill (an A-only TE copy paired with a B-only copy elsewhere) sits in
    a better chain's A gap, lands outside that chain's B gap, and its chain is top level in neither genome: it
    is dropped. B's frame alone keeps no piece: top level in B sequence that no chain covers means no context,
    not synteny. An inversion is kept as an in-gap inverted fill, as its own top-level chain, or, where a
    better chain's extension covers part of it, through its chain's top-level pieces. Known costs:
    real transpositions shorter than CHAIN_MAX_GAP are dropped like paralogs, and so is an ortholog on an
    unplaced scaffold of either genome whose other side falls in a chromosome chain's gap; their elements come
    out missing (D/O/U), never P/E. There is no parameter, threshold, slop or size escape: a trimmed piece overlaps no
    kept piece, so it lies wholly inside or wholly outside any parent's gap and extent.
    Cost: per frame and sequence, the kept extents are filed in 64-kb buckets (`_Extents`); an extent that
    contains a piece contains the piece's start, so finding the parent scans that one bucket, and bisecting
    the parent's pieces finds the two flanking it: O(b + log k) per piece, b being the kept extents touching
    the bucket (the nesting depth plus the chains with pieces there: a few) and k the parent's pieces. Filing
    a kept chain costs O(k log k + extent / 64 kb).
    """
    if len(blocks) == 0:
        return np.zeros(0, dtype=BLOCK_DTYPE)
    if (blocks["chain"] < 0).any():
        raise ValueError("net_blocks needs chained blocks: run chain_blocks first")
    srt = blocks[np.lexsort((blocks["a_end"], blocks["a_start"], blocks["chain"], -blocks["chain_score"]))]
    same = np.diff(srt["chain"]) == 0
    if (same & ((np.diff(srt["a_chrom"]) != 0) | (np.diff(srt["b_chrom"]) != 0)
                | (np.diff(srt["strand"]) != 0))).any():
        raise ValueError("net_blocks needs chain_blocks output: a chain on one A sequence, one B sequence "
                         "and one strand")
    claimed_a: defaultdict[int, _Claims] = defaultdict(_Claims)
    claimed_b: defaultdict[int, _Claims] = defaultdict(_Claims)
    kept = _KeptChains()
    tally: Counter = Counter()
    rows: list[tuple] = []
    comp: list[float] = []  # best competitor score per kept row
    for chain in np.split(srt, np.flatnonzero(~same) + 1):
        _net_chain(chain, claimed_a, claimed_b, rows, comp, kept, tally)
    _log_net(tally, len(rows), label)
    if not rows:
        return np.zeros(0, dtype=BLOCK_DTYPE)

    out = np.array(rows, dtype=BLOCK_DTYPE)
    best = np.array(comp)
    margin = np.ones(len(out))
    hit = best > 0  # a positive competitor came later, so chain_score >= best > 0
    margin[hit] = 1 - best[hit] / out["chain_score"][hit]
    out["margin"] = margin
    return out[np.lexsort((out["a_start"], out["a_chrom"]))]


def _net_chain(chain: np.ndarray, claimed_a: dict, claimed_b: dict, rows: list, comp: list,
               kept: _KeptChains, tally: Counter) -> None:
    """Net one chain (blocks in A order): credit and trim each block, then test the pieces for synteny
    (`net_blocks`); append and claim the kept pieces if they keep >= MIN_CHAIN_SCORE matches."""
    own_a: defaultdict[int, _Claims] = defaultdict(_Claims)
    own_b: defaultdict[int, _Claims] = defaultdict(_Claims)
    pieces = []
    for blk in chain:
        ac, a0, a1, bc, b0, b1, strand, matches, alen, cid, score, _ = blk.item()
        hits_a, hits_b = claimed_a[ac].overlapping(a0, a1), claimed_b[bc].overlapping(b0, b1)
        _credit(hits_a, a0, a1, score, comp)
        _credit(hits_b, b0, b1, score, comp)
        for p0, p1, q0, q1 in _trim(blk, hits_a, hits_b, own_a[ac], own_b[bc]):
            frac = (p1 - p0) / (a1 - a0)
            pieces.append((ac, p0, p1, bc, q0, q1, strand, round(matches * frac), round(alen * frac),
                           cid, score, np.nan))
    if sum(pc[7] for pc in pieces) < MIN_CHAIN_SCORE:
        return
    places = [kept.places(pc) for pc in pieces]
    top_a, top_b = any(a == _TOP for a, _ in places), any(b == _TOP for _, b in places)
    keep = []
    for pc, (a, b) in zip(pieces, places):
        a = _CHAIN if a == _OUT and top_a else a
        b = _CHAIN if b == _OUT and top_b else b
        tally[a, b] += 1
        if a == _OUT or b == _OUT:
            tally["dropped"] += 1
            tally["dropped_bp"] += pc[2] - pc[1]
        else:
            keep.append((pc, _kind(a, b)))
    tally["tested"] += 1
    tally["cut"] += len(keep) < len(pieces)
    if sum(pc[7] for pc, _ in keep) < MIN_CHAIN_SCORE:  # it had >= MIN_CHAIN_SCORE before the drops
        tally["lost"] += 1
        return
    for pc, kind in keep:
        claimed_a[pc[0]].add(pc[1], pc[2], len(rows))
        claimed_b[pc[3]].add(pc[4], pc[5], len(rows))
        rows.append(pc)
        comp.append(0.0)
        tally[kind] += 1
    tally["kept"] += 1
    kept.add([pc for pc, _ in keep])


def _kind(a: str, b: str) -> str:
    """A kept piece's place in the net summary, from its places in A's and B's frames: chain-supported if
    either frame needs its chain, else an inversion or syn if either frame says so, else top level in both."""
    for kind in (_CHAIN, _INV, _SYN):
        if kind in (a, b):
            return kind
    return _TOP


def _log_net(tally: Counter, n_rows: int, label: str = "") -> None:
    """The net summary line (kept pieces by place, dropped pieces), prefixed by `label` (the pair), then the
    per-frame details; both logged at DEBUG."""
    prefix = f"{label}: " if label else ""
    log.debug("%ssyntenic net: kept %d pieces (%d top in both genomes, %d syn, %d inv, %d chain-supported); "
             "dropped %d out of place, %d bp of A in total", prefix, n_rows, tally[_TOP], tally[_SYN],
             tally[_INV], tally[_CHAIN], tally["dropped"], tally["dropped_bp"])
    if log.isEnabledFor(logging.DEBUG):
        frames = ", ".join(f"{a}/{b} {n}" for (a, b), n in sorted(
            ((k, n) for k, n in tally.items() if isinstance(k, tuple)), key=lambda kv: -kv[1]))
        log.debug("%ssyntenic net: %d chains tested, %d kept, %d chains lost pieces, %d of these fell below "
                  "MIN_CHAIN_SCORE; pieces by place in A/B: %s", prefix, tally["tested"], tally["kept"],
                  tally["cut"], tally["lost"], frames or "none")


def _credit(hits: list, start: int, end: int, score: float, comp: list) -> None:
    """Record `score` as a competitor of each claimed piece that [start, end) overlaps by > CHAIN_OVERLAP bp."""
    for cs, ce, owner in hits:
        if min(ce, end) - max(cs, start) > CHAIN_OVERLAP and score > comp[owner]:
            comp[owner] = score


def _trim(blk: np.void, hits_a: list, hits_b: list, own_a: _Claims, own_b: _Claims) -> list:
    """Unclaimed (a0, a1, b0, b1) pieces of one block, each also added to the chain's own claims.

    hits_a / hits_b: pieces of other chains claimed on the block's A / B span. Own claims are re-queried
    per free A part, so pieces of this block never overlap each other either. Every piece is >= ANCHOR_MIN
    on A, so a shorter block yields none.
    """
    a0, a1, b0, b1 = int(blk["a_start"]), int(blk["a_end"]), int(blk["b_start"]), int(blk["b_end"])
    if a1 - a0 < ANCHOR_MIN:
        return []
    taken_a = sorted(hits_a + own_a.overlapping(a0, a1))
    if not taken_a and not hits_b and not own_b.overlapping(b0, b1):
        own_a.add(a0, a1, -1)
        own_b.add(b0, b1, -1)
        return [(a0, a1, b0, b1)]
    pieces = []
    for x0, x1 in _subtract(a0, a1, taken_a):
        y0, y1 = _a_span_to_b(blk, x0, x1)
        for z0, z1 in _subtract(y0, y1, sorted(hits_b + own_b.overlapping(y0, y1))):
            p0, p1 = _b_span_to_a(blk, z0, z1)
            p0, p1 = max(p0, x0), min(p1, x1)
            if p1 - p0 >= ANCHOR_MIN:
                pieces.append((p0, p1, z0, z1))
                own_a.add(p0, p1, -1)
                own_b.add(z0, z1, -1)
    return pieces


def _subtract(start: int, end: int, taken: list) -> list[tuple[int, int]]:
    """Parts of [start, end) not covered by `taken` ((start, end, ...) intervals sorted by start)."""
    parts = []
    cur = start
    for ts, te, _ in taken:
        if te <= cur:
            continue
        if ts >= end:
            break
        if ts > cur:
            parts.append((cur, ts))
        cur = te
        if cur >= end:
            break
    if cur < end:
        parts.append((cur, end))
    return parts


def _a_span_to_b(blk: np.void, x0: int, x1: int) -> tuple[int, int]:
    """B image of A interval [x0, x1) under the block's linear map (`a_to_b` at the end positions)."""
    if blk["strand"] > 0:
        return a_to_b(blk, x0), a_to_b(blk, x1 - 1) + 1
    return a_to_b(blk, x1 - 1), a_to_b(blk, x0) + 1


def _b_span_to_a(blk: np.void, z0: int, z1: int) -> tuple[int, int]:
    """A interval of B interval [z0, z1): the inverse linear map applied to its boundaries.

    Boundaries (not end positions) keep the images of disjoint B intervals disjoint; callers clip to the source
    A part.
    """
    a0, a1, b0, b1 = int(blk["a_start"]), int(blk["a_end"]), int(blk["b_start"]), int(blk["b_end"])
    r = (a1 - a0) / (b1 - b0)
    if blk["strand"] > 0:
        return a0 + round((z0 - b0) * r), a0 + round((z1 - b0) * r)
    return a0 + round((b1 - z1) * r), a0 + round((b1 - z0) * r)


class _Claims:
    """Disjoint claimed intervals on one sequence, each with its owner (output row; -1 = own chain, pending).

    Each interval is filed in every fixed-width bucket it touches, a bucket being a list sorted by start; queries
    and inserts touch only the buckets an interval spans. (One sorted list per chromosome would make every
    insert O(claims on that chromosome), quadratic at maize scale.)
    """

    __slots__ = ("buckets",)

    def __init__(self) -> None:
        self.buckets: dict[int, list[tuple[int, int, int]]] = {}

    def add(self, start: int, end: int, owner: int) -> None:
        item = (start, end, owner)
        for k in range(start >> CLAIM_BUCKET_BITS, ((end - 1) >> CLAIM_BUCKET_BITS) + 1):
            bucket = self.buckets.get(k)
            if bucket is None:
                self.buckets[k] = [item]
            else:
                insort(bucket, item)

    def overlapping(self, start: int, end: int) -> list[tuple[int, int, int]]:
        """Claimed (start, end, owner) intervals overlapping [start, end), sorted by start, each once."""
        hits: list[tuple[int, int, int]] = []
        if not self.buckets:
            return hits
        last = -1  # start of the last reported interval; a long interval recurs in later buckets
        for k in range(start >> CLAIM_BUCKET_BITS, ((end - 1) >> CLAIM_BUCKET_BITS) + 1):
            bucket = self.buckets.get(k)
            if bucket is None:
                continue
            for i in range(bisect_right(bucket, start, key=_END), len(bucket)):
                item = bucket[i]
                if item[0] >= end:
                    break
                if item[0] > last:
                    hits.append(item)
                    last = item[0]
        return hits


class _KeptChains:
    """The kept chains as the synteny test reads them (`net_blocks`): per A sequence their A extents, per B
    sequence their B extents (`_Extents`), each with the chain's sequence in the other genome, its strand and
    its pieces in that frame's order (with their starts, for bisection)."""

    __slots__ = ("a", "b")

    def __init__(self) -> None:
        self.a: defaultdict[int, _Extents] = defaultdict(_Extents)
        self.b: defaultdict[int, _Extents] = defaultdict(_Extents)

    def add(self, pieces: list[tuple]) -> None:
        """File one kept chain: its piece rows (one A sequence, one B sequence, one strand)."""
        by_a, by_b = sorted(pieces, key=itemgetter(1)), sorted(pieces, key=itemgetter(4))
        ac, bc, strand = by_a[0][0], by_a[0][3], by_a[0][6]
        self.a[ac].add(by_a[0][1], by_a[-1][2], (bc, strand, [p[1] for p in by_a], by_a))
        self.b[bc].add(by_b[0][4], by_b[-1][5], (ac, strand, [p[4] for p in by_b], by_b))

    def places(self, pc: tuple) -> tuple[str, str]:
        """The piece row's place in A's frame and in B's: _TOP, _SYN, _INV or _OUT (`_place`)."""
        return (_place(self.a.get(pc[0]), pc, _A_SPAN, _B_SPAN, pc[3]),
                _place(self.b.get(pc[3]), pc, _B_SPAN, _A_SPAN, pc[0]))


def _place(extents: _Extents | None, pc: tuple, own: tuple[int, int], other: tuple[int, int],
           other_seq: int) -> str:
    """The piece row's place in one genome's frame: _TOP without a parent (the innermost kept extent on its
    sequence containing its span there, fields `own`); _SYN or _INV if its span in the other genome (fields
    `other`, on sequence `other_seq`) lies inside the parent's gap there, on the parent's strand or the other;
    else _OUT.

    The parent's two pieces flanking the piece in this frame bound that gap. They exist: the piece lies inside
    the parent's extent and overlaps none of its (claimed) pieces.
    """
    start, end = pc[own[0]], pc[own[1]]
    parent = extents.innermost(start, end) if extents is not None else None
    if parent is None:
        return _TOP
    seq, strand, starts, pieces = parent
    if seq != other_seq:
        return _OUT
    i = bisect_right(starts, start)
    left, right = pieces[i - 1], pieces[i]
    lo, hi = other
    g0, g1 = (left[hi], right[lo]) if strand > 0 else (right[hi], left[lo])
    if g0 <= pc[lo] and pc[hi] <= g1:
        return _SYN if pc[6] == strand else _INV
    return _OUT


class _Extents:
    """Extents [start, end) of kept chains on one sequence, each with its chain's record (`_KeptChains`).

    As in `_Claims`, each extent is filed in every 64-kb bucket it touches, in filing (chain) order; an extent
    containing [start, end) contains `start`, so `innermost` scans that one bucket.
    """

    __slots__ = ("buckets",)

    def __init__(self) -> None:
        self.buckets: dict[int, list[tuple[int, int, tuple]]] = {}

    def add(self, start: int, end: int, rec: tuple) -> None:
        item = (start, end, rec)
        for k in range(start >> CLAIM_BUCKET_BITS, ((end - 1) >> CLAIM_BUCKET_BITS) + 1):
            bucket = self.buckets.get(k)
            if bucket is None:
                self.buckets[k] = [item]
            else:
                bucket.append(item)

    def innermost(self, start: int, end: int) -> tuple | None:
        """The record of the shortest extent containing [start, end) (the first filed on ties), or None."""
        best, size = None, 0
        for s, e, rec in self.buckets.get(start >> CLAIM_BUCKET_BITS, ()):
            if s <= start and end <= e and (best is None or e - s < size):
                best, size = rec, e - s
        return best
