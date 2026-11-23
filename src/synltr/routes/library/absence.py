"""Each record's syntenic window in the other genome of a pair: where its ortholog would sit (LTRmix's).

A record of genome x lies in a collinear block when its midpoint is on the block's x chromosome, between the block's
first and last anchor; the first such block (blocks in their order) is its window's. The window: the other genome's
positions (anchor midpoints, 1-based) between the anchors flanking the record, one more anchor of slack on each side
(lo, hi); room: the span between the two immediately flanking anchors' partners, how much sequence the other genome
holds where the insertion would be. A record in no block has no window (block -1, the rest -1). Inside a window,
no ortholog is an inferred absence (A).
"""
from __future__ import annotations

from typing import NamedTuple

import numpy as np


class Windows(NamedTuple):
    """Per record of one genome: block index (-1 none), window lo and hi on the other genome, room (-1 none)."""

    bid: np.ndarray
    lo: np.ndarray
    hi: np.ndarray
    room: np.ndarray


def windows(lx, blocks, side: str) -> Windows:
    """The windows of the records of lx (Records), genome a (`side` "a") or b ("b") of the pair's blocks."""
    n = lx.n
    bid, lo, hi, room = (np.full(n, -1, np.int64) for _ in range(4))
    for k, blk in enumerate(blocks):
        ax, ay, cx = (blk.ax, blk.bx, blk.chrom_a) if side == "a" else (blk.bx, blk.ax, blk.chrom_b)
        o = np.argsort(ax)
        ax, ay = ax[o], ay[o]
        sel = np.where((lx.chrom == cx) & (lx.mid >= ax[0]) & (lx.mid <= ax[-1]))[0]
        if len(sel) == 0:
            continue
        j = np.searchsorted(ax, lx.mid[sel], side="right")
        j0 = np.clip(j - 2, 0, len(ax) - 1)   # one anchor of slack on each side
        j1 = np.clip(j + 1, 0, len(ax) - 1)
        y0, y1 = np.minimum(ay[j0], ay[j1]), np.maximum(ay[j0], ay[j1])
        t0 = np.clip(j - 1, 0, len(ax) - 1)   # the immediate flanking anchors
        t1 = np.clip(j, 0, len(ax) - 1)
        tr = np.abs(ay[t1] - ay[t0])
        new = bid[sel] < 0
        s2 = sel[new]
        bid[s2], lo[s2], hi[s2], room[s2] = k, y0[new], y1[new], tr[new]
    return Windows(bid, lo, hi, room)
