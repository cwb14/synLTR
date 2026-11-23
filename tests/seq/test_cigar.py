"""CIGAR parsing and walks (seq.cigar): op codes, lengths, positions and aligned blocks, 0-based."""
import numpy as np
import pysam
import pytest

from synltr.seq.cigar import OPS, aligned_in, blocks, deletion_at, parse, qpos_at, walk


def _pairs(ops, lens):
    return list(zip(lens.tolist(), ops.tolist()))


def test_parse_lastz_cigarx():
    """LASTZ writes a count of 1 without digits ("X=" is 1X1=)."""
    assert _pairs(*parse("1001=12I3X=", "lastz")) == [(1001, 7), (12, 1), (3, 8), (1, 7)]
    assert _pairs(*parse("X=", "lastz")) == [(1, 8), (1, 7)] and _pairs(*parse("5M2D", "lastz")) == [(5, 0), (2, 2)]
    ops, lens = parse("", "lastz")
    assert ops.dtype == np.uint8 and lens.dtype == np.int64 and len(ops) == len(lens) == 0
    for bad in ("12", "3Q", "=4", "5N"):
        with pytest.raises(ValueError, match="malformed LASTZ cigar"):
            parse(bad, "lastz")


def test_parse_sam_cigar():
    """SAM and PAF cg:Z: CIGARs: every op has a count; all nine ops."""
    assert _pairs(*parse("10M2I3D4N5S6H7P8=9X", "sam")) == [(10, 0), (2, 1), (3, 2), (4, 3), (5, 4), (6, 5), (7, 6),
                                                           (8, 7), (9, 8)]
    assert OPS == "MIDNSHP=X"
    for bad in ("", "=4", "4", "4Q", "X="):
        with pytest.raises(ValueError, match="malformed CIGAR"):
            parse(bad, "sam")
    with pytest.raises(ValueError, match="dialect"):
        parse("4M", "bam")


def test_walk_positions_before_each_op():
    ops, lens = parse("3=2I4D1X2S5N", "sam")
    q, r = walk(ops, lens, 10, 100)
    assert q.tolist() == [10, 13, 15, 15, 16, 18, 18] and r.tolist() == [100, 103, 103, 107, 108, 108, 113]
    q, r = walk(*parse("", "lastz"))
    assert q.tolist() == [0] and r.tolist() == [0]


def test_blocks_are_the_aligned_ops():
    """(query start, reference start, length) per M, = or X op: I and S move the query, D and N the reference,
    H and P neither."""
    ops, lens = parse("2H3S10=2I5X3D4M1N6=1P", "sam")
    assert blocks(ops, lens, 0, 1000) == [(3, 1000, 10), (15, 1010, 5), (20, 1018, 4), (24, 1023, 6)]
    assert blocks(ops.tolist(), lens.tolist(), 0, 1000) == blocks(ops, lens, 0, 1000)  # lists work alike
    assert blocks(*parse("", "lastz"), 5, 5) == []
    with pytest.raises(ValueError, match="op code 9"):
        blocks([9], [1], 0, 0)


def _aln(start: int, cigar: str, length: int = 100_000) -> pysam.AlignedSegment:
    """A read alignment at 0-based `start` of a sequence c, query bases A (hard clips hold none)."""
    a = pysam.AlignedSegment(pysam.AlignmentHeader.from_dict({"SQ": [{"SN": "c", "LN": length}]}))
    a.reference_id, a.reference_start = 0, start
    a.cigarstring = cigar
    a.query_sequence = "A" * a.query_length
    return a


def test_qpos_at():
    a = _aln(100, "10S20M5D20M10I", 1000)  # query: 10 clip, 20 match, (5 del), 20 match, 10 ins (no trailing ref)
    ct, r0 = a.cigartuples, a.reference_start
    assert qpos_at(ct, r0, 100) == 10 and qpos_at(ct, r0, 119) == 29
    assert qpos_at(ct, r0, 122) == 30  # inside the deletion -> query position after it
    assert qpos_at(ct, r0, 125) == 30 and qpos_at(ct, r0, 144) == 49 and qpos_at(ct, r0, 145) is None
    assert qpos_at(ct, r0, 99) is None


def test_qpos_at_counts_hard_clips_and_keeps_the_op_order():
    """A supplementary alignment's hard clip counts (the position in the read as sequenced); inside a deletion the
    position is the one at that op, so an insertion before the deletion counts and one after it does not."""
    ct = _aln(100, "50H10M3I5D2I10M").cigartuples
    assert qpos_at(ct, 100, 100) == 50 and qpos_at(ct, 100, 109) == 59
    assert qpos_at(ct, 100, 112) == 63  # in the deletion: after the 3-bp insertion, before the 2-bp one
    assert qpos_at(ct, 100, 115) == 65  # first base after the deletion, past both insertions
    assert blocks(*zip(*ct), 0, 100) == [(0, 100, 10), (15, 115, 10)]  # blocks count neither the clip nor the op order


def test_anchor_inside_deletion_is_invalid():
    """A read whose anchor falls inside a large alignment deletion lacks that sequence (benchmark bug 2)."""
    a = _aln(1000, "100M2000D100M10D10M")
    ct, r0 = a.cigartuples, a.reference_start
    assert deletion_at(ct, r0, 1500, 50) and deletion_at(ct, r0, 3099, 50)
    assert not deletion_at(ct, r0, 1050, 50) and not deletion_at(ct, r0, 3150, 50)
    assert not deletion_at(ct, r0, 3205, 50)  # 10-bp deletion ignored


def test_deletion_at_interval():
    a = _aln(1000, "100M2000D100M")
    ct, r0 = a.cigartuples, a.reference_start
    assert deletion_at(ct, r0, 1500, 50) == (1100, 3100)
    assert deletion_at(ct, r0, 1050, 50) is None and deletion_at(ct, r0, 3150, 50) is None


def test_deletion_ops_count_one_by_one():
    """Two deletions split by an insertion are two deletions: neither reaches 50 bp although together they do."""
    ct = _aln(0, "100M30D2I30D100M").cigartuples
    assert deletion_at(ct, 0, 110, 50) is None and deletion_at(ct, 0, 140, 50) is None
    assert deletion_at(ct, 0, 110, 30) == (100, 130) and deletion_at(ct, 0, 140, 30) == (130, 160)
    assert deletion_at(_aln(0, "10M60N10M").cigartuples, 0, 40, 50) == (10, 70)  # a skip (N) counts as well


def test_aligned_in():
    """Aligned (M, = and X) reference bases inside [lo, hi): deletions, skips and insertions add none."""
    ct = _aln(1000, "5S100M20D10I50=30X").cigartuples
    assert aligned_in(ct, 1000, 1000, 1100) == 100 and aligned_in(ct, 1000, 1050, 1150) == 50 + 30
    assert aligned_in(ct, 1000, 1100, 1120) == 0 and aligned_in(ct, 1000, 1110, 1200) == 80
    assert aligned_in(ct, 1000, 0, 10_000) == 180 and aligned_in(ct, 1000, 2000, 3000) == 0
