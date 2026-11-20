"""Anchor-to-window construction: A and B sequence windows around each top-level element.

Coordinates: blocks and windows are 0-based half-open; `Element` is 1-based inclusive.
Minus-strand B windows are returned in genome coordinates (the caller reverse-complements).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, replace

import numpy as np

from ltrtrace.base.constants import ANCHOR_MIN, NCAND, WINDOW_CAP, WINDOW_PAD
from ltrtrace.io.elements import Element
from ltrtrace.routes.assembly.blocks import a_to_b

# the longest window make_windows builds: an anchor span aR - aL of WINDOW_CAP, padded to [aL - WINDOW_PAD + 1,
# aR + WINDOW_PAD); a merged window may be as long (_try_merge)
PADDED_CAP = WINDOW_CAP + 2 * WINDOW_PAD - 1


@dataclass
class Window:
    wid: int
    a_seqid: str
    a_start: int
    a_end: int  # 0-based half-open
    b_seqid: str | None
    b_start: int
    b_end: int
    strand: int  # +1/-1; B window is reverse-complemented by the caller when -1
    eids: list[str]  # top-level elements (descendants implied)
    status: str  # ok | no_anchor | no_left_anchor | no_right_anchor | breakpoint | span_cap
    left_anchor: int  # A position (0-based) of the left anchor's last anchored base (ok windows)
    right_anchor: int  # A position of the right anchor's first anchored base (ok windows)
    margin: float  # min anchor-block margin

    def anchor_distances(self, start: int, end: int) -> tuple[int, int]:
        """bp from an element of the window (1-based inclusive `start`, `end`) to the window's left and right
        anchors: each element's own, in a merged window or nested in a host too; (0, 0) without anchors (a
        window not ok)."""
        if self.status != "ok":
            return 0, 0
        return start - 2 - self.left_anchor, self.right_anchor - end


def _left_cands(starts, ends, s0):
    """Nearest-first blocks ending (clipped) at <= s0 with clipped length >= ANCHOR_MIN.

    Returns (row index, last anchored A base). Blocks are non-overlapping and sorted,
    so ends decrease going left and the WINDOW_CAP search can stop early.
    """
    out = []
    i = int(np.searchsorted(starts, s0, side="left")) - 1
    while i >= 0 and len(out) < NCAND:
        ce = min(int(ends[i]), s0)
        if ce < s0 - WINDOW_CAP:
            break
        if ce - int(starts[i]) >= ANCHOR_MIN:
            out.append((i, ce - 1))
        i -= 1
    return out


def _right_cands(starts, ends, t0):
    """Mirror of `_left_cands`: (row index, first anchored A base)."""
    out = []
    i = int(np.searchsorted(ends, t0, side="right"))
    n = len(starts)
    while i < n and len(out) < NCAND:
        cs = max(int(starts[i]), t0)
        if cs > t0 + WINDOW_CAP:
            break
        if int(ends[i]) - cs >= ANCHOR_MIN:
            out.append((i, cs))
        i += 1
    return out


def _unresolved(wid, seqid, s0, t0, eid, status):
    return Window(wid, seqid, s0, t0, None, 0, 0, 0, [eid], status, 0, 0, math.nan)


def _try_merge(cur: Window, new: Window) -> bool:
    """Merge `new` into `cur` (in place) if both are ok, overlap on A, share B's seqid and strand, touch or
    overlap on B, and the merged window spans at most PADDED_CAP on A and on B: the cap make_windows puts on
    one window's anchor span (WINDOW_CAP), on the padded spans windows carry. The merged window's anchors are
    the outermost ones (each element's distances to them: Window.anchor_distances). True if merged."""
    if (cur.status != "ok" or new.status != "ok" or new.a_start >= cur.a_end or new.a_end <= cur.a_start
            or cur.b_seqid != new.b_seqid or cur.strand != new.strand
            or new.b_start > cur.b_end or cur.b_start > new.b_end):
        return False
    a0, a1 = min(cur.a_start, new.a_start), max(cur.a_end, new.a_end)
    b0, b1 = min(cur.b_start, new.b_start), max(cur.b_end, new.b_end)
    if a1 - a0 > PADDED_CAP or b1 - b0 > PADDED_CAP:
        return False
    cur.left_anchor = min(cur.left_anchor, new.left_anchor)
    cur.right_anchor = max(cur.right_anchor, new.right_anchor)
    cur.a_start, cur.a_end, cur.b_start, cur.b_end = a0, a1, b0, b1
    cur.eids += new.eids
    cur.margin = min(cur.margin, new.margin)
    return True


def make_windows(elements: list[Element], blocks: np.ndarray, a_names: list[str],
                 b_names: list[str], b_lengths: dict[str, int],
                 a_lengths: dict[str, int]) -> list[Window]:
    """One window per top-level element, merged where overlapping windows share a B locus.

    `blocks` must be 1:1 and sorted by (a_chrom, a_start) (the `net_blocks` output).
    Unresolved elements get one never-merged window with a reason in `status`
    (A interval = the element itself, b_seqid None, strand 0).
    """
    a_idx = {n: i for i, n in enumerate(a_names)}
    ach = blocks["a_chrom"]
    by_seq: dict[str, list[Element]] = {}
    for e in elements:
        if e.level == 0:
            if e.seqid not in a_idx or e.seqid not in a_lengths:
                raise ValueError(f"element {e.eid}: seqid {e.seqid!r} not in a_names and a_lengths")
            by_seq.setdefault(e.seqid, []).append(e)

    out: list[Window] = []
    for seqid in sorted(by_seq, key=a_idx.__getitem__):
        els = sorted(by_seq[seqid], key=lambda e: (e.start, e.end))
        lo, hi = np.searchsorted(ach, [a_idx[seqid], a_idx[seqid] + 1])
        blk = blocks[lo:hi]
        starts, ends = blk["a_start"], blk["a_end"]
        a_len = a_lengths[seqid]
        chrom_wins: list[Window] = []
        last_ok: Window | None = None
        for e in els:
            s0, t0 = e.start - 1, e.end
            if len(blk) == 0:
                chrom_wins.append(_unresolved(0, seqid, s0, t0, e.eid, "no_anchor"))
                continue
            lc = _left_cands(starts, ends, s0)
            rc = _right_cands(starts, ends, t0)
            if not lc or not rc:
                chrom_wins.append(_unresolved(0, seqid, s0, t0, e.eid,
                                              "no_left_anchor" if not lc else "no_right_anchor"))
                continue
            pairs = sorted(((i, j) for i in range(len(lc)) for j in range(len(rc))),
                           key=lambda p: (p[0] + p[1], p[0]))
            win = None
            reason = None
            for n, (i, j) in enumerate(pairs):
                bl, br = blk[lc[i][0]], blk[rc[j][0]]
                aL, aR = lc[i][1], rc[j][1]
                same = (bl["b_chrom"] == br["b_chrom"] and bl["strand"] == br["strand"])
                bL, bR = a_to_b(bl, aL), a_to_b(br, aR)
                ordered = bL < bR if bl["strand"] > 0 else bL > bR
                small = aR - aL <= WINDOW_CAP and abs(bR - bL) <= WINDOW_CAP
                if same and ordered and small:
                    bname = b_names[int(bl["b_chrom"])]
                    b_len = b_lengths[bname]
                    if bl["strand"] > 0:
                        b0, b1 = bL - WINDOW_PAD + 1, bR + WINDOW_PAD
                    else:
                        b0, b1 = bR - WINDOW_PAD + 1, bL + WINDOW_PAD
                    win = Window(
                        0, seqid, max(0, aL - WINDOW_PAD + 1), min(a_len, aR + WINDOW_PAD),
                        bname, max(0, b0), min(b_len, b1), int(bl["strand"]), [e.eid], "ok",
                        aL, aR, float(min(bl["margin"], br["margin"])))
                    break
                if n == 0:
                    reason = "breakpoint" if not (same and ordered) else "span_cap"
            if win is None:
                chrom_wins.append(_unresolved(0, seqid, s0, t0, e.eid, reason))
                continue
            if last_ok is not None and _try_merge(last_ok, win):
                continue
            chrom_wins.append(win)
            last_ok = win
        chrom_wins.sort(key=lambda w: w.a_start)  # stable: element order breaks ties
        out.extend(chrom_wins)
    return [replace(w, wid=i) for i, w in enumerate(out)]
