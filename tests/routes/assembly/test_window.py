"""Windows around the top-level elements (routes.assembly.window)."""
import math

import numpy as np

from ltrtrace.base.constants import WINDOW_CAP, WINDOW_PAD
from ltrtrace.io.elements import Element
from ltrtrace.routes.assembly.blocks import BLOCK_DTYPE
from ltrtrace.routes.assembly.window import PADDED_CAP, Window, _try_merge, make_windows

A_NAMES = ["a0", "a1"]
B_NAMES = ["b0", "b1"]
A_LEN = {"a0": 2_000_000, "a1": 100_000}
B_LEN = {"b0": 2_000_000, "b1": 2_000_000}


def blk(*rows):
    """rows: (a_chrom, a0, a1, b_chrom, b0, b1, strand[, margin])."""
    out = np.zeros(len(rows), dtype=BLOCK_DTYPE)
    for i, r in enumerate(rows):
        m = r[7] if len(r) > 7 else 1.0
        out[i] = (r[0], r[1], r[2], r[3], r[4], r[5], r[6], 0, r[2] - r[1], 0, 0.0, m)
    return out[np.lexsort((out["a_start"], out["a_chrom"]))]


def el(start, end, seqid="a0", level=0):
    return Element("A", f"A|{seqid}:{start}-{end}", "n", seqid, start, end, "+", (start, start),
                   (end, end), None, None, None, None, None, None, None if level == 0 else "p",
                   level)


def run(els, b):
    return make_windows(els, b, A_NAMES, B_NAMES, B_LEN, A_LEN)


def test_window_basic():
    b = blk((0, 10_000, 20_000, 0, 110_000, 120_000, 1, 0.9),
            (0, 30_000, 40_000, 0, 130_000, 140_000, 1, 0.7))
    [w] = run([el(25_001, 26_000)], b)
    assert w.status == "ok" and w.strand == 1 and w.b_seqid == "b0" and w.wid == 0
    assert (w.a_start, w.a_end) == (19_999 - WINDOW_PAD + 1, 30_000 + WINDOW_PAD)
    assert (w.b_start, w.b_end) == (119_999 - WINDOW_PAD + 1, 130_000 + WINDOW_PAD)
    assert (w.left_anchor, w.right_anchor) == (19_999, 30_000)
    assert w.anchor_distances(25_001, 26_000) == (5_000, 4_000)
    assert w.margin == np.float32(0.7) and w.eids == ["A|a0:25001-26000"]


def test_window_breakpoint():
    b = blk((0, 10_000, 20_000, 0, 110_000, 120_000, 1), (0, 30_000, 40_000, 1, 130_000, 140_000, 1))
    [w] = run([el(25_001, 26_000)], b)
    assert w.status == "breakpoint" and w.b_seqid is None and w.strand == 0
    assert math.isnan(w.margin) and w.anchor_distances(25_001, 26_000) == (0, 0)


def test_window_missing_right():
    [w] = run([el(25_001, 26_000)], blk((0, 10_000, 20_000, 0, 110_000, 120_000, 1)))
    assert w.status == "no_right_anchor"


def test_window_span_cap():
    b = blk((0, 10_000, 20_000, 0, 100_000, 110_000, 1),
            (0, 30_000, 40_000, 0, 710_000, 720_000, 1))
    [w] = run([el(25_001, 26_000)], b)
    assert w.status == "span_cap"


def test_windows_merge_overlap():
    b = blk((0, 10_000, 20_000, 0, 110_000, 120_000, 1, 0.9),
            (0, 30_000, 40_000, 0, 130_000, 140_000, 1, 0.6),
            (0, 50_000, 60_000, 0, 150_000, 160_000, 1, 0.8))
    es = [el(25_001, 26_000), el(27_001, 28_000)]
    ws = run(es, b)
    assert len(ws) == 1 and ws[0].eids == [es[0].eid, es[1].eid] and ws[0].status == "ok"
    assert (ws[0].left_anchor, ws[0].right_anchor) == (19_999, 30_000)
    # each element's own distances to the merged window's anchors, not one pair for both
    assert [ws[0].anchor_distances(e.start, e.end) for e in es] == [(5_000, 4_000), (7_000, 2_000)]
    assert ws[0].margin == np.float32(0.6)


def test_merge_refused_cases():
    def w(a0, a1, bs="b0", b0=0, b1=100, strand=1, status="ok"):
        return Window(0, "a0", a0, a1, bs, b0, b1, strand, ["x"], status, a0, a1 - 1, 0.5)

    assert _try_merge(w(0, 100), w(50, 150, b0=100, b1=200))          # B touching
    assert not _try_merge(w(0, 100), w(50, 150, b0=101, b1=200))      # B gap
    assert not _try_merge(w(0, 100), w(50, 150, bs="b1"))
    assert not _try_merge(w(0, 100), w(50, 150, strand=-1))
    assert not _try_merge(w(0, 100), w(100, 150))                     # A only touching
    assert not _try_merge(w(100, 200), w(50, 100))                    # ... on the other side
    assert not _try_merge(w(0, 100), w(50, 150, status="breakpoint"))
    assert _try_merge(w(0, 100), w(50, PADDED_CAP))                   # as long as the longest single window
    assert not _try_merge(w(0, 100), w(50, PADDED_CAP + 1))           # span cap on A
    assert not _try_merge(w(0, 100), w(50, 150, b0=100, b1=PADDED_CAP + 1))  # and on B


def test_padded_cap_matches_make_windows():
    """make_windows caps the anchor span aR - aL at WINDOW_CAP and pads it to WINDOW_CAP + 2 * WINDOW_PAD - 1
    bases: the longest window it builds is PADDED_CAP long."""
    def window(span):  # the element's window between anchors aL = 19_999 and aR = aL + span (A and B)
        b = blk((0, 10_000, 20_000, 0, 110_000, 120_000, 1),
                (0, 19_999 + span, 29_999 + span, 0, 119_999 + span, 129_999 + span, 1))
        [w] = run([el(25_001, 26_000)], b)
        return w

    w = window(WINDOW_CAP)
    assert w.status == "ok" and w.a_end - w.a_start == w.b_end - w.b_start == PADDED_CAP
    assert window(WINDOW_CAP + 1).status == "span_cap"


def test_identical_windows_under_the_cap_merge():
    """Two elements between the same anchors, an anchor span of WINDOW_CAP - 100: their windows are identical and
    longer than WINDOW_CAP once padded, and merge into one (unmerged, each would be aligned once per
    element)."""
    a_r = 19_999 + WINDOW_CAP - 100  # the first anchored base on the right; aL = 19_999
    b = blk((0, 10_000, 20_000, 0, 110_000, 120_000, 1), (0, a_r, a_r + 10_000, 0, a_r + 100_000, a_r + 110_000, 1))
    es = [el(25_001, 26_000), el(300_001, 301_000)]
    [w] = run(es, b)
    assert w.status == "ok" and w.eids == [es[0].eid, es[1].eid]
    assert (w.a_start, w.a_end) == (20_000 - WINDOW_PAD, a_r + WINDOW_PAD)
    assert w.a_end - w.a_start == PADDED_CAP - 100 > WINDOW_CAP


def test_window_minus_strand():
    b = blk((0, 10_000, 20_000, 0, 180_000, 190_000, -1), (0, 30_000, 40_000, 0, 160_000, 170_000, -1))
    [w] = run([el(25_001, 26_000)], b)
    assert w.status == "ok" and w.strand == -1 and w.b_start < w.b_end
    # aL=19999 -> bL=180000; aR=30000 -> bR=169999
    assert (w.b_start, w.b_end) == (169_999 - WINDOW_PAD + 1, 180_000 + WINDOW_PAD)
    assert w.b_start < 175_000 < w.b_end


def test_window_at_sequence_start():
    [w] = run([el(4, 1000)], blk((0, 5_000, 6_000, 0, 5_000, 6_000, 1)))
    assert w.status == "no_left_anchor"


def test_no_blocks_on_seqid():
    [w] = run([el(1000, 2000, "a1")], blk((0, 10_000, 20_000, 0, 110_000, 120_000, 1)))
    assert w.status == "no_anchor"


def test_single_block_spans_element():
    b = blk((0, 10_000, 40_000, 0, 110_000, 140_000, 1))
    [w] = run([el(25_001, 26_000)], b)
    assert w.status == "ok" and w.anchor_distances(25_001, 26_000) == (0, 0)
    assert (w.a_start, w.a_end) == (25_000 - WINDOW_PAD, 26_000 + WINDOW_PAD)
    assert (w.b_start, w.b_end) == (125_000 - WINDOW_PAD, 126_000 + WINDOW_PAD)


def test_try_next_nearest_rescue():
    b = blk((0, 10_000, 20_000, 0, 110_000, 120_000, 1),
            (0, 21_000, 24_000, 1, 900_000, 903_000, 1),   # nearest left: wrong chromosome
            (0, 30_000, 40_000, 0, 130_000, 140_000, 1))
    [w] = run([el(25_001, 26_000)], b)
    assert w.status == "ok" and w.b_seqid == "b0" and w.anchor_distances(25_001, 26_000)[0] == 5_000


def test_short_block_not_an_anchor():
    b = blk((0, 10_000, 20_000, 0, 110_000, 120_000, 1), (0, 24_950, 25_000, 0, 124_950, 125_000, 1),
            (0, 30_000, 40_000, 0, 130_000, 140_000, 1))
    [w] = run([el(25_001, 26_000)], b)
    assert w.anchor_distances(25_001, 26_000)[0] == 5_000


def test_clip_at_sequence_ends():
    b = blk((1, 0, 300, 0, 100_000, 100_300, 1), (1, 1000, 1300, 0, 100_800, 101_100, 1))
    [w] = run([el(401, 900, "a1")], b)  # a1 length 100_000; anchors close to start
    assert w.status == "ok" and w.a_start == 0
    b = blk((1, 98_000, 98_300, 1, 1_998_000, 1_998_300, 1), (1, 99_700, 100_000, 1, 1_999_700, 2_000_000, 1))
    [w] = run([el(99_001, 99_500, "a1")], b)
    assert w.a_end == 100_000 and w.b_end == 2_000_000 and w.status == "ok"


def test_nested_ignored_and_ordering():
    b = blk((0, 10_000, 20_000, 0, 110_000, 120_000, 1), (0, 30_000, 40_000, 0, 130_000, 140_000, 1))
    ws = run([el(25_501, 25_800, level=1), el(25_001, 26_000), el(1000, 2000, "a1")], b)
    assert [w.a_seqid for w in ws] == ["a0", "a1"] and [w.wid for w in ws] == [0, 1]
    assert len(ws[0].eids) == 1


def test_merge_keeps_the_outermost_anchors():
    def w(a0, a1):  # anchors WINDOW_PAD inside the ends, as make_windows pads them
        return Window(0, "a0", a0, a1, "b0", 0, 100, 1, ["x"], "ok", a0 + WINDOW_PAD - 1, a1 - WINDOW_PAD, 0.5)

    cur = w(5_000, 10_000)
    assert _try_merge(cur, w(6_000, 9_000))  # inside: keeps both anchors
    assert (cur.left_anchor, cur.right_anchor) == (5_499, 9_500)
    assert _try_merge(cur, w(9_000, 15_000))  # extends right only
    assert (cur.left_anchor, cur.right_anchor) == (5_499, 14_500)
    assert _try_merge(cur, w(1_000, 6_000))  # extends left only
    assert (cur.left_anchor, cur.right_anchor) == (1_499, 14_500)
    assert not _try_merge(cur, w(0, 1_000))  # ends where cur starts: no overlap on A


def test_anchor_search_stops_at_window_cap():
    """Anchor candidates are searched only within WINDOW_CAP of the element on A."""
    right = (0, 600_000, 610_000, 0, 1_600_000, 1_610_000, 1)
    for gap, status in ((WINDOW_CAP, "span_cap"), (WINDOW_CAP + 1, "no_left_anchor")):
        start = 20_000 + gap + 1  # 1-based: s0 - (the left block's end) = gap
        [w] = run([el(start, start + 999)], blk((0, 10_000, 20_000, 0, 1_010_000, 1_020_000, 1), right))
        assert w.status == status, gap  # within reach the pair is tried (and exceeds the span cap)
        assert (w.a_start, w.a_end) == (start - 1, start + 999)  # an unresolved window is the element itself


def test_three_windows_merge():
    b = blk((0, 10_000, 20_000, 0, 110_000, 120_000, 1), (0, 30_000, 40_000, 0, 130_000, 140_000, 1))
    es = [el(21_001, 22_000), el(24_001, 25_000), el(27_001, 28_000)]
    [w] = run(es, b)
    assert w.eids == [e.eid for e in es] and (w.left_anchor, w.right_anchor) == (19_999, 30_000)
    assert [w.anchor_distances(e.start, e.end) for e in es] == [(1_000, 8_000), (4_000, 5_000), (7_000, 2_000)]


def test_candidate_pairs_tried_by_i_plus_j_then_i():
    """When the nearest pair fails (B images out of order), (0, 1) is tried before (1, 0):
    the nearest left anchor with the second right one, though the second left with the nearest right also
    works."""
    b = blk((0, 10_000, 15_000, 0, 50_000, 55_000, 1),  # L1: works with R0
            (0, 18_000, 20_000, 0, 200_000, 202_000, 1),  # L0, the nearest left: works with R1 only
            (0, 30_000, 32_000, 0, 100_000, 102_000, 1),  # R0, the nearest right
            (0, 35_000, 40_000, 0, 300_000, 305_000, 1))  # R1
    [w] = run([el(25_001, 26_000)], b)
    assert w.status == "ok" and (w.left_anchor, w.right_anchor) == (19_999, 35_000)
