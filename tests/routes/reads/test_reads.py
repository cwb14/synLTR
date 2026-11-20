"""Per-read allele evidence at one element (routes.reads.reads): allele lengths and junction reads."""
import pysam

from ltrtrace.routes.reads.reads import Site, Target, _one_sided


def test_solo_with_unequal_ltrs():
    """A solo LTR keeps one of two LTRs of different length: either length change is solo (benchmark bug 1)."""
    el = Target("x", "c", 1000, 14000, "+", (1000, 2429), (13844, 14000))  # LTRs 1,429 and 156 bp
    st = Site(el, 5)
    L = 13000
    assert st.classify_delta(-(L - 1429)) == "S" and st.classify_delta(-(L - 156)) == "S"
    assert st.classify_delta(-(L - 800)) == "S"
    assert st.classify_delta(0) == "F" and st.classify_delta(-(L + 5)) == "E"
    assert st.classify_delta(-3000) == "O"


def test_boundary_tolerance_without_tsd():
    """Empty and solo windows widen only by the requested extra (annotation error learned from the sample)."""
    el = Target("x", "c", 10000, 20000, "+", (10000, 10500), (19500, 20000))
    st = Site(el, 5)
    d = -(10000 + 5) + 80
    assert st.classify_delta(d) == "O"
    assert st.classify_delta(d, extra=60) == "E"
    assert st.classify_delta(0, extra=200) == "F" and st.classify_delta(400, extra=200) == "O"  # full unaffected


def _site():
    # element 10,000-20,000; LTRs 500 bp: left LTR ends at a = 10,500, right LTR starts at b = 19,500
    el = Target("x", "c", 10000, 20000, "+", (10000, 10500), (19500, 20000))
    return Site(el, 5)


def _aln(start, cigar):
    hdr = pysam.AlignmentHeader.from_dict({"SQ": [{"SN": "c", "LN": 100000}]})
    a = pysam.AlignedSegment(hdr)
    a.reference_id, a.reference_start, a.mapping_quality = 0, start, 60
    a.cigarstring = cigar
    a.query_sequence = "A" * a.query_length
    return a


def test_junction_read_with_internal_deletion_is_full():
    """A read of a full allele with a >= 50-bp deletion just past the LTR still crosses the junction (v0.1 known
    issue: such reads were dropped because the base 100 bp into the internal region lay inside the deletion)."""
    st = _site()
    # 2 kb of left flank, the whole left LTR, 20 bp internal, 80-bp deletion, 400 bp internal
    x = _aln(8000, "2000M500M20M80D400M")
    assert _one_sided(st, None, "c", [x], "L") == "F_left"


def test_solo_like_read_is_not_full():
    """A read that jumps from the LTR over the internal region (deletion) is not full-junction evidence."""
    st = _site()
    x = _aln(8000, "2000M500M8990D50M")  # LTR, then deleted internal region and right LTR
    assert _one_sided(st, None, "c", [x], "L") != "F_left"
    y = _aln(8000, "1950M3000D300M")  # deletion over the junction (another structure)
    assert _one_sided(st, None, "c", [y], "L") is None or _one_sided(st, None, "c", [y], "L") != "F_left"


def test_target_geometry():
    """Target holds ltrkin's element, 0-based half-open: length, mean LTR length, internal region and the TSD
    length assumed without a TSD (5 bp)."""
    el = Target("x", "c", 100, 1100, "-", (100, 300), (1000, 1100), tsd="")
    assert (el.length, el.ltr_mean, el.internal, el.tsd_len) == (1000, 150.0, (300, 1000), 5)
    assert Target("y", "c", 0, 50, "+", (0, 10), (40, 50), tsd="ACG").tsd_len == 3
