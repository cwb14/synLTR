"""State classifier (routes.assembly.classify): the toy genotype cases (tests/builders.py), aligned with LASTZ, and
the rules one by one."""
import functools
import math
import random
from dataclasses import replace

import builders as toy
import numpy as np
import pytest
from markers import needs_lastz

from ltrtrace.base.constants import MIN_SEG, WINDOW_PAD
from ltrtrace.calls.records import CALL_COLUMNS
from ltrtrace.calls.states import PAIR_STATES
from ltrtrace.io.annotation import assign_nesting
from ltrtrace.io.elements import ElementTable
from ltrtrace.io.seqio import revcomp
from ltrtrace.routes.assembly.classify import Call, call_row, classify_window
from ltrtrace.routes.assembly.lastz import AlnMap, align_window
from ltrtrace.routes.assembly.window import Window

# 7 is the base seed; the others exercise chance microhomology at the junctions (measured with
# refine_gap): 4 and 9 shift E's leftmost gap 1 bp left of its TSD copy, 6 shifts adjacent_E's gap 2 bp
# left, 11 shifts nested_E's 2 bp, 14 gives E an 8-bp slide (3 bp of chance beyond the TSD). Seed 7
# shifts adjacent_E's gap 1 bp left.
SEEDS = (7, 0, 4, 6, 9, 11, 14)
TSD = toy.T


def _window(c, **kw):
    """A window over a whole toy case (wid 0, A = c.a on seqid a, B = c.b on seqid b, strand 1),
    with fields `kw` replaced."""
    win = Window(wid=0, a_seqid="a", a_start=0, a_end=len(c.a), b_seqid="b", b_start=0, b_end=len(c.b),
                 strand=1, eids=[e.eid for e in c.elements if e.parent is None], status="ok", left_anchor=0,
                 right_anchor=len(c.a) - 1, margin=1.0)
    return replace(win, **kw)


def _children(elements):
    out = {}
    for e in elements:
        if e.parent is not None:
            out.setdefault(e.parent, []).append(e.eid)
    return out


def run_case(c, win=None, aln="align", b=None):
    """classify_window over the case's window; `aln` defaults to align_window(c.a, b) (LASTZ)."""
    b = c.b if b is None else b
    win = _window(c, b_end=len(b)) if win is None else win
    if isinstance(aln, str):
        aln = align_window(c.a, b)
    return classify_window(win, aln, c.a, b, {e.eid: e for e in c.elements}, _children(c.elements))


@functools.cache
def _calls(name, seed):
    c = toy.case(name, random.Random(seed))
    return c, {x.eid: x for x in run_case(c)}


def _by_name(c, calls):
    return {e.name: calls[e.eid] for e in c.elements}


def amap(*segs, inv=()):
    """AlnMap from (a, b, len[, mm]) tuples."""
    rows = [s if len(s) == 4 else (*s, 0) for s in segs]
    return AlnMap(*(np.array(col, dtype=np.int64) for col in zip(*rows)), inv=list(inv))


# ---------------------------------------------------------------- states, flags and tiers of the toy cases


@pytest.mark.parametrize("seed", SEEDS)
@pytest.mark.parametrize("name", toy.CASES)
@needs_lastz
def test_case_state(name, seed):
    c, calls = _calls(name, seed)
    assert set(calls) == {e.eid for e in c.elements}  # one call per element, descendants included
    for eid, state in c.expect.items():
        assert calls[eid].state == state, (name, seed, eid, calls[eid])
        assert calls[eid].state in PAIR_STATES


@pytest.mark.parametrize("name", toy.CASES)
@needs_lastz
def test_flags_and_tiers(name):
    for seed in SEEDS:
        c, calls = _calls(name, seed)
        for eid, flags in c.flags.items():
            assert flags <= set(calls[eid].flags), (name, seed, calls[eid])
        for eid, flags in c.noflags.items():
            assert not flags & set(calls[eid].flags), (name, seed, calls[eid])
        for eid, tier in c.tiers.items():
            assert calls[eid].tier == tier, (name, seed, calls[eid])
        for eid, reason in c.reasons.items():
            assert calls[eid].ev["reason"] == reason, (name, seed, calls[eid])
        for x in calls.values():  # G and U get tier "."; every other state a level; only U has a reason
            if x.state in ("G", "U"):
                assert x.tier == "."
            else:
                assert x.tier in ("high", "medium", "low") and x.ev["reason"] is None


@pytest.mark.parametrize("name", ["host_E", "nested_E", "adjacent_E", "P"])
@needs_lastz
def test_window_status_propagates(name):
    c = toy.case(name, random.Random(7))
    # an unresolved window as make_windows builds it: no B locus, strand 0, margin nan
    win = _window(c, b_seqid=None, b_start=0, b_end=0, strand=0, status="breakpoint", margin=math.nan)
    for aln in (None, "align"):
        calls = run_case(c, win=win, aln=aln)
        assert len(calls) == len(c.elements)
        for x in calls:
            e = next(e for e in c.elements if e.eid == x.eid)
            assert (x.state, x.tier, x.ev["reason"]) == ("U", ".", "breakpoint")
            assert (x.b_seqid, x.b_start, x.b_end, x.b_strand) == (None, None, None, None)
            assert x.flags == (["inherit"] if e.parent else [])


@needs_lastz
def test_minus_strand_window_flags_inv():
    """B's genome carries the case on its minus strand; the caller hands classify_window B reverse-
    complemented (A's orientation), so only the B coordinates and the inv flag change."""
    c = toy.case("E", random.Random(7))
    rng = random.Random(70)
    pre, post = toy.rand_seq(1000, rng), toy.rand_seq(700, rng)
    genome = pre + revcomp(c.b) + post
    b0, b1 = len(pre), len(pre) + len(c.b)
    assert revcomp(genome[b0:b1]) == c.b
    plus = genome[:b0] + c.b + genome[b1:]  # the same window on the plus strand
    for strand, g, want in ((1, plus, TSD), (-1, genome, revcomp(TSD))):
        [x] = run_case(c, win=_window(c, b_start=b0, b_end=b1, strand=strand))
        assert x.state == "E" and "tsd_match" in x.flags and ("inv" in x.flags) == (strand == -1)
        assert x.b_seqid == "b" and x.b_strand == ("+" if strand == 1 else "-")
        assert x.b_end - x.b_start + 1 == len(TSD) and g[x.b_start - 1:x.b_end] == want
        if strand == 1:
            assert (x.b_start, x.b_end) == (b0 + 3001, b0 + 3005)  # the B copy at window [3000, 3005)
        else:
            assert (x.b_start, x.b_end) == (b1 - 3004, b1 - 3000)  # genome = b_end - 1 - local, flipped


# ---------------------------------------------------------------- B coordinates and evidence


@needs_lastz
def test_b_coordinates_per_state():
    """Window B starts at 500 of its sequence: P/S/T the aligned element image, E the TSD copy, D the two
    flank bases around the junction, O the B insertion, G none. 1-based inclusive."""
    off = 500
    seen = {}
    for name in ("P", "S", "T", "E", "D_left", "O", "G"):
        c, _ = _calls(name, 7)
        [x] = run_case(c, win=_window(c, b_start=off, b_end=off + len(c.b)))
        seen[name] = (x.state, x.b_start, x.b_end)
        if x.state == "G":
            assert (x.b_seqid, x.b_start, x.b_end, x.b_strand) == (None, None, None, None)
        else:
            assert (x.b_seqid, x.b_strand) == ("b", "+")
    assert seen["P"] == ("P", off + 3006, off + 7805)  # the element itself (equal lengths, no indel)
    assert seen["S"] == ("S", off + 3006, off + 3405)  # the solo LTR
    assert seen["T"] == ("T", off + 3006, off + 3005 + 400 + 1600 + 400)
    assert seen["E"] == ("E", off + 3001, off + 3005)  # B's TSD copy
    assert seen["D_left"] == ("D", off + 2500, off + 2501)  # B = L[:-500] + R: L[2499] and R[0]
    c, _ = _calls("O", 7)  # the 2-kb insertion between B's TSD copies, as the local aligner bounds it
    g = align_window(c.a, c.b).gap_around(5405)
    assert abs(g.b0 - 3005) <= 5 and abs(g.b1 - 5005) <= 5
    assert seen["O"] == ("O", off + g.b0 + 1, off + g.b1)


@needs_lastz
def test_evidence_and_row():
    c, calls = _calls("E", 7)
    e = c.elements[0]
    x = calls[e.eid]
    assert (x.ev["dl"], x.ev["dr"]) == (2, 3)  # balanced in the 5-bp slide: the TSD split 2 / 3
    assert (x.ev["b_ins"], x.ev["tsd_a"], x.ev["tsd_b"], x.ev["tsd_mm"]) == (0, TSD, TSD, 0)
    assert x.ev["cov_left"] == x.ev["cov_right"] == x.ev["cov_int"] == x.ev["cov_elem"] == 0.0
    assert x.ev["ident"] == 1.0 and x.ev["margin"] == 1.0 and x.ev["reason"] is None
    win = _window(c, margin=np.float32(0.75), left_anchor=np.int64(12))
    [y] = run_case(c, win=win)
    row = call_row(y, e, win)
    assert (row["anchor_left"], row["anchor_right"]) == (e.start - 2 - 12, len(c.a) - 1 - e.end)
    assert tuple(row) == CALL_COLUMNS
    assert row["eid"] == e.eid and (row["a_seqid"], row["a_start"], row["a_end"]) == ("a", e.start, e.end)
    assert (row["parent"], row["level"], row["state"], row["tier"]) == (None, 0, "E", "high")
    assert row["flags"] == "tsd_match" and row["reason"] is None and row["wid"] == 0
    assert row["margin"] == 0.75 and type(row["margin"]) is float and type(row["anchor_left"]) is int
    for v in row.values():
        assert v is None or type(v) in (str, int, float)
    with pytest.raises(ValueError, match="eid"):
        call_row(y, replace(e, eid="other"), win)


@needs_lastz
def test_row_flags_and_absent_values():
    c, calls = _calls("host_G", 7)
    win = _window(c)
    by = _by_name(c, calls)
    h, ch = sorted(c.elements, key=lambda e: e.level)
    row = call_row(by["c"], ch, win)
    assert (row["state"], row["tier"], row["flags"], row["parent"], row["level"]) == ("G", ".", "inherit",
                                                                                      h.eid, 1)
    assert row["b_seqid"] is None and row["b_start"] is None and row["tsd_b"] is None and row["dl"] is None
    assert call_row(by["h"], h, win)["flags"] is None  # no flags: None, not an empty string


# ---------------------------------------------------------------- rules


def test_constants():
    assert CALL_COLUMNS[:10] == ("eid", "a_seqid", "a_start", "a_end", "parent", "level", "state", "tier",
                                 "flags", "reason")
    assert len(CALL_COLUMNS) == len(set(CALL_COLUMNS)) == 29
    assert set(Call.__dataclass_fields__) == {"eid", "state", "tier", "flags", "b_seqid", "b_start", "b_end",
                                              "b_strand", "ev"}


@needs_lastz
def test_inheritance_copies_host_site():
    """An inner element of an E or G host takes the host's state, tier and B site (the host's site is
    its site); its own evidence is its own."""
    for name, flag in (("host_E", "host_empty"), ("host_G", "inherit")):
        c, calls = _calls(name, 7)
        by = _by_name(c, calls)
        h, ch = by["h"], by["c"]
        assert (ch.state, ch.tier, ch.flags) == (h.state, h.tier, [flag])
        assert (ch.b_seqid, ch.b_start, ch.b_end, ch.b_strand) == (h.b_seqid, h.b_start, h.b_end, h.b_strand)
        assert ch.ev["cov_elem"] == 0.0 and ch.ev["tsd_a"] == toy.T2 and ch.ev["dl"] is None


def test_no_alignment():
    c = toy.case("host_E", random.Random(7))
    calls = run_case(c, aln=None)
    assert [(x.state, x.ev["reason"]) for x in calls] == [("U", "no_alignment")] * 2
    assert calls[1].flags == ["inherit"]


@needs_lastz
def test_tsd_allowance_k_rescues_e_mm2():
    """An old element (K = 0.3) allows 2 TSD mismatches: E_mm2's B copy (2 off) then matches."""
    c = toy.case("E_mm2", random.Random(7))
    old = replace(c, elements=[replace(c.elements[0], k=0.3)])
    [x] = run_case(old)
    assert x.state == "E" and "tsd_match" in x.flags and x.tier == "high" and x.ev["tsd_mm"] == 2


@needs_lastz
def test_small_ins_and_filler():
    """B-only sequence of w < n < MIN_SEG bp: O + small_ins at an otherwise empty site, D + filler where
    flank is lost as well."""
    rng = random.Random(5)
    c = toy.case("E", rng)
    L, R = c.a[:3000], c.a[-3000:]
    filler = "A" * 30  # cannot align to random flank
    for b, state, flag in ((L + TSD + filler + R, "O", "small_ins"), (L[:-500] + filler + R, "D", "filler")):
        [x] = run_case(c, b=b)
        assert x.state == state and flag in x.flags and x.ev["b_ins"] > 10, x
        assert x.ev["b_ins"] < MIN_SEG


@needs_lastz
def test_balanced_placement_ignores_reported_placement():
    """An aligner can report the empty-site deletion at any of its equal placements; the gap is refined
    to the leftmost and E/D is read on the balanced one, so every report gives the same call."""
    c, calls = _calls("E", 7)
    [want] = calls.values()
    n = len(c.a) - len(c.b)  # 4805: the element and one TSD copy
    for a0 in range(2995, 3011):  # off the equal placements (3000..3005) too: refinement fixes those
        m = amap((0, 0, a0), (a0 + n, a0, len(c.a) - a0 - n))
        [x] = run_case(c, aln=m)
        assert (x.state, x.flags, x.ev["dl"], x.ev["dr"], x.ev["tsd_b"]) == ("E", ["tsd_match"], 2, 3, TSD)
        assert (x.b_start, x.b_end) == (want.b_start, want.b_end)


@needs_lastz
def test_flank_loss_beyond_slop_is_d():
    """B = L[:-n] + T + R loses n bp of flank next to the TSD, so the gap overhangs the element by n + 5 on
    the left (the lost flank and A's left TSD copy) and by 0 on the right: E while n + 5 <= w = 10 (the
    balanced placement only lowers the left side), D well beyond it."""
    c = toy.case("E", random.Random(3))
    L, R = c.a[:3000], c.a[-3000:]
    for n, state in ((2, "E"), (5, "E"), (6, "D"), (20, "D"), (300, "D")):
        [x] = run_case(c, b=L[:-n] + TSD + R)
        assert x.state == state, (n, x)
        if n in (5, 6):  # the boundary itself: dl = w is E, w + 1 is D
            assert x.ev["dl"] == n + 5, (n, x)


@needs_lastz
def test_adjacent_requires_neighbor_inside_gap():
    """The overhang is explained only by another element deleted with e; if B keeps the neighbor (a
    deletion of e1 that also took 300 bp of e2) or the neighbor is e's host, e is D."""
    c = toy.case("adjacent_E", random.Random(7))
    e1, e2 = c.elements
    L = c.a[:3000]
    s2 = e2.start - 1
    partial = L + TSD + c.a[s2 + 300:]  # B keeps all of e2 but its first 300 bp
    calls = {x.eid: x for x in run_case(c, b=partial)}
    assert calls[e2.eid].state == "T"  # B lost e2's first 300 bp, most of its left LTR
    assert calls[e1.eid].state == "D" and "adjacent" not in calls[e1.eid].flags
    # the host is not "another element in the gap": host_S's child is D (its host is S, not P)
    hc, hcalls = _calls("host_S", 7)
    assert "adjacent" not in _by_name(hc, hcalls)["c"].flags


def test_zero_length_gap_is_ambiguous():
    """A host whose own sequence is gone while its inner element aligns: the host's midpoint lies inside
    the child's aligned run, so there is no gap to read: U ambiguous."""
    c = toy.case("host_E", random.Random(7))
    h, ch = sorted(c.elements, key=lambda e: e.level)
    s_c, t_c = ch.start - 1, ch.end
    n = len(c.a)
    # B = L T c R: left flank, the child, right flank aligned; the host's own parts unaligned
    m = amap((0, 0, 3005), (s_c, 3005, t_c - s_c), (n - 3000, 3005 + t_c - s_c, 3000))
    b = c.a[:3005] + c.a[s_c:t_c] + c.a[-3000:]
    calls = {x.eid: x for x in run_case(c, aln=m, b=b)}
    assert (calls[h.eid].state, calls[h.eid].ev["reason"]) == ("U", "ambiguous")
    assert calls[ch.eid].state == "U" and calls[ch.eid].flags == ["inherit"]


def test_flank_rules():
    """Rule 2: < MIN_SEG aligned bp on a side -> U flank_unaligned; gap_around without a flanking run
    (aligned bases only in pieces too short to make one) -> U flank_unaligned as well."""
    c = toy.case("E", random.Random(7))
    n = len(c.a)
    m = amap((0, 0, 3000), (n - 70, 3005 + 2930, 70))  # 70 bp on the right
    [x] = run_case(c, aln=m)
    assert (x.state, x.ev["reason"]) == ("U", "flank_unaligned")
    # 3 x 40 bp on the right separated by B gaps >= MIN_SEG: >= MIN_SEG aligned, but no run
    m = amap((0, 0, 3000), (n - 2000, 3100, 40), (n - 1500, 3300, 40), (n - 1000, 3500, 40))
    [x] = run_case(c, aln=m)
    assert (x.state, x.ev["reason"]) == ("U", "flank_unaligned")


@needs_lastz
def test_inverted_hit_rules():
    """inv_site needs an inverted hit of >= MIN_SEG A bases overlapping the element whose B interval meets
    the site [lo, hi] (the B images of the last aligned A base before the element and the first after it);
    any other inverted hit is ignored (the site is then read as usual)."""
    c, _ = _calls("E", 7)
    n = len(c.a)
    base = [(0, 0, 3000), (7805, 3000, n - 7805)]  # site: lo = 2999 (A 2999), hi = 3000 (A 7805)
    for b0, b1 in ((2990, 2990 + MIN_SEG), (3000, 3000 + MIN_SEG), (3000 - MIN_SEG, 3000)):  # meet [2999, 3000]
        [x] = run_case(c, aln=amap(*base, inv=[(5000, 5000 + MIN_SEG, b0, b1)]))
        assert (x.state, x.ev["reason"]) == ("U", "inv_site"), (b0, b1)
    for b0, b1 in ((10, 10 + MIN_SEG), (3001, 3001 + MIN_SEG), (2999 - MIN_SEG, 2999)):  # miss the site
        [x] = run_case(c, aln=amap(*base, inv=[(5000, 5000 + MIN_SEG, b0, b1)]))
        assert x.state == "E", (b0, b1)
    [x] = run_case(c, aln=amap(*base, inv=[(5000, 5000 + MIN_SEG - 1, 2990, 2989 + MIN_SEG)]))  # too short
    assert x.state == "E"
    [x] = run_case(c, aln=amap(*base, inv=[(100, 2000, 2000, 3900)]))  # outside the element
    assert x.state == "E"


@pytest.mark.parametrize("seed", SEEDS)
@needs_lastz
def test_inverted_copy_off_the_site_is_not_inv_site(seed):
    """B also carries an unrelated reverse-complemented copy of the element's LTR 2 kb past the site (say,
    another element's fragment). Both of A's LTRs align to it inverted, but away from the site: E stays E."""
    c = toy.case("E", random.Random(seed))
    [e] = c.elements
    ltr = c.a[e.ltr_left[0] - 1:e.ltr_left[1]]
    site = toy.FLANK + len(toy.T)  # B = L T R: R starts here
    b = c.b[:site + 2000] + revcomp(ltr) + c.b[site + 2000:]
    aln = align_window(c.a, b)
    s, t = e.start - 1, e.end
    assert any(a1 - a0 >= MIN_SEG and a0 < t and a1 > s and b0 >= site + 2000 - 10 and b1 <= site + 2410
               for a0, a1, b0, b1 in aln.inv), aln.inv  # the inverted copy is seen, off the site
    [x] = run_case(c, aln=aln, b=b)
    assert x.state == "E" and "tsd_match" in x.flags, (x.state, x.ev["reason"], x.flags)


@pytest.mark.parametrize("pad", [None, WINDOW_PAD])
@needs_lastz
def test_inv_site_at_every_seed(pad):
    """The toy inverted site stays U inv_site at every seed: its inverted hit lies at the site. Whole case,
    and cut to WINDOW_PAD bp of flank on each side as make_windows cuts a window."""
    for seed in SEEDS:
        c = toy.case("inv_site", random.Random(seed))
        cut = 0 if pad is None else toy.FLANK - pad
        a, b = c.a[cut:len(c.a) - cut], c.b[cut:len(c.b) - cut]
        win = _window(c, a_start=cut, a_end=len(c.a) - cut, b_start=cut, b_end=len(c.b) - cut)
        [x] = classify_window(win, align_window(a, b), a, b, {e.eid: e for e in c.elements}, {})
        assert (x.state, x.ev["reason"]) == ("U", "inv_site"), (seed, x.state, x.ev["reason"])


@needs_lastz
def test_low_identity_and_margin_lower_the_tier():
    c = toy.case("E", random.Random(7))
    [x] = run_case(c, win=_window(c, margin=0.4))  # placement not trusted, evidence holds
    assert x.state == "E" and x.tier == "medium"
    c2 = toy.case("E_mm2", random.Random(7))
    [x] = run_case(c2, win=_window(c2, margin=0.4))  # neither
    assert x.state == "E" and x.tier == "low"


@needs_lastz
def test_window_offset():
    """Elements in genome coordinates, the window starting at 10 kb of the chromosome: same calls."""
    off = 10_000
    for name in ("E", "nested_E", "S", "adjacent_E"):
        c, calls = _calls(name, 7)
        els = list(assign_nesting(ElementTable.from_rows([
            replace(e, start=e.start + off, end=e.end + off, ltr_left=(e.ltr_left[0] + off, e.ltr_left[1] + off),
                    ltr_right=(e.ltr_right[0] + off, e.ltr_right[1] + off), eid=f"{e.eid}+", parent=None, level=0)
            for e in c.elements])))
        moved = replace(c, elements=els)
        win = _window(moved, a_start=off, a_end=off + len(c.a))
        got = {x.eid[:-1]: x for x in run_case(moved, win=win)}
        for eid, x in calls.items():
            assert (got[eid].state, got[eid].flags, got[eid].tier, got[eid].b_start) == (
                x.state, x.flags, x.tier, x.b_start)


@needs_lastz
def test_elements_outside_window_ignored():
    """`elements` may hold the whole genome's records: only the window's elements are called."""
    c = toy.case("E", random.Random(7))
    other = replace(c.elements[0], eid="elsewhere", start=50_000, end=55_000)
    win = _window(c)
    calls = classify_window(win, align_window(c.a, c.b), c.a, c.b,
                            {e.eid: e for e in [*c.elements, other]}, {})
    assert [x.eid for x in calls] == [c.elements[0].eid]


@needs_lastz
def test_element_outside_window_rejected():
    c = toy.case("E", random.Random(7))
    win = _window(c, a_end=5000)
    with pytest.raises(ValueError, match="window"):
        classify_window(win, align_window(c.a[:5000], c.b), c.a[:5000], c.b,
                        {e.eid: e for e in c.elements}, {})


@pytest.mark.parametrize(("name", "seed"), [("S", 15), ("S", 30), ("S", 115), ("S_notsd", 75),
                                            ("S_notsd", 125), ("S_notsd", 126)])
@needs_lastz
def test_solo_tsd_read_from_b(name, seed):
    """solo_tsd is read from B beside the solo LTR, not from how the TSD copies align. Seeds where chance
    microhomology moves the reported deletion 2-3 bp into S's left TSD copy (15, 30, 115) keep the flag;
    seeds where S_notsd's LTR ends in TCTC (125, 126) or C (75), so the leftmost placement aligns 4 or 1
    bases of the TSD copy B lacks instead of the LTR's end, do not get it."""
    c, calls = _calls(name, seed)
    x = calls[c.elements[0].eid]
    assert x.state == "S" and ("solo_tsd" in x.flags) == (name == "S")
    assert x.tier == ("high" if name == "S" else "medium")


@needs_lastz
def test_solo_tsd_reads_beside_the_annotated_ltr_ends():
    """An LTR end base left unaligned (chance microhomology can give it to the TSD) does not shift the read:
    the solo LTR's B ends are extrapolated from the nearest aligned LTR base to the annotated LTR ends, so the
    TSD reads sit right beside the LTR. S (B = L T LTR T R) with the left LTR's last, or first, base
    unaligned in a hand-built map."""
    c, _ = _calls("S", 7)
    n, ltr = len(c.a), (3005, 3405)  # the left LTR in A, and in B
    assert c.b[ltr[1]:ltr[1] + 5] == TSD == c.b[ltr[0] - 5:ltr[0]]
    # reads starting one base off, at the last aligned base, would miss the TSD by >= 4 mismatches
    assert min(sum(p != q for p, q in zip(c.b[3404:3409], TSD)),
               sum(p != q for p, q in zip(c.b[3001:3006], TSD))) >= 4
    last_off = amap((0, 0, 3404), (7805, 3405, n - 7805))  # A, B 3404 (the LTR's last base) unaligned
    first_off = amap((0, 0, 3005), (3006, 3006, 399), (7805, 3405, n - 7805))  # A, B 3005 unaligned
    for m in (last_off, first_off):
        [x] = run_case(c, aln=m)
        assert x.state == "S" and "solo_tsd" in x.flags and x.tier == "high", x


@needs_lastz
def test_element_remnant_is_d_not_e():
    """B keeps 50-70 bp of the element between its TSD copies (below MIN_SEG, so the absent branch): the gap
    ends that far inside the element (dl' or dr' < -w), so B carries part of the insertion: D, never E."""
    for seed in (7, 1):
        c = toy.case("E", random.Random(seed))
        L, R, e = c.a[:3000], c.a[-3000:], c.a[3005:-3005]
        for piece in (e[:50], e[-50:], e[:70], e[-70:]):
            [x] = run_case(c, b=L + TSD + piece + TSD + R)
            assert x.state == "D" and x.tier == "high" and min(x.ev["dl"], x.ev["dr"]) < -10, (seed, x)
            assert "tsd_match" not in x.flags and x.ev["b_ins"] == 0


@pytest.mark.parametrize("seed", [6, 7])
@needs_lastz
def test_tsd_read_over_slide_range(seed):
    """adjacent_E: the refined leftmost gap starts 2 (seed 6) or 1 (seed 7) bp left of e1, by chance
    microhomology, and ends past e1, so no TSD copy of e1 lies outside it; at placement x = s - g.a0 (<= the
    slide) the gap starts at e1 and e1's left copy is read: B[3000, 3005) = T. e2's right copy is read at the
    leftmost placement: B[3005, 3010) = T2."""
    c, calls = _calls("adjacent_E", seed)
    by = _by_name(c, calls)
    for name, tsd, (b0, b1) in (("e1", TSD, (3001, 3005)), ("e2", toy.T2, (3006, 3010))):
        x = by[name]
        assert x.state == "E" and "tsd_match" in x.flags and x.tier == "high", x
        assert (x.ev["tsd_b"], x.ev["tsd_mm"], x.b_start, x.b_end) == (tsd, 0, b0, b1)


# ---------------------------------------------------------------- rule 5a: present calls are anchored at their site


def _paralog_case(seed, tsd):
    """A = L T e T R. B = L X1 [T] e [T] X2 R: the element's bases sit in B between random 500-bp pieces, away from
    the flank images (rule 5a's trigger: a same-family copy elsewhere in the window); the TSD copies T are
    in B beside it iff `tsd`. Returns (case, hand-built AlnMap, B): L, e and R aligned, the gaps >= MIN_SEG."""
    rng = random.Random(seed)
    c = toy.case("P", rng)
    [e] = c.elements
    s, t, n = e.start - 1, e.end, len(c.a)
    x1, x2 = toy.rand_seq(500, rng), toy.rand_seq(500, rng)
    mid = (TSD if tsd else "") + c.a[s:t] + (TSD if tsd else "")
    b = c.a[:3000] + x1 + mid + x2 + c.a[-3000:]
    b_e = 3000 + 500 + (len(TSD) if tsd else 0)
    m = amap((0, 0, 3000), (s, b_e, t - s), (n - 3000, len(b) - 3000, 3000))
    return c, m, b


def test_unanchored_present_is_u():
    """The element's bases align to a copy between B's flank images, both flanks at their own images and no TSD
    beside the copy: not the element at its site -> U unanchored (rule 5 alone would call it P)."""
    for seed in SEEDS:
        c, m, b = _paralog_case(seed, tsd=False)
        [x] = run_case(c, aln=m, b=b)
        assert (x.state, x.tier, x.ev["reason"]) == ("U", ".", "unanchored"), (seed, x)
        assert x.ev["cov_elem"] == 1.0 and x.b_start is None  # the evidence stays; no B site


def test_unanchored_with_tsd_stands():
    """The same geometry with B's bases on both sides of the copy matching A's TSD (the solo_tsd read): the
    element sits in a site of its own, so the call stands."""
    for seed in SEEDS:
        c, m, b = _paralog_case(seed, tsd=True)
        [x] = run_case(c, aln=m, b=b)
        assert x.state == "P" and x.tier == "high" and x.ev["reason"] is None, (seed, x)


def _one_sided_case(name, side):
    """A = L T e T R with B's element junction broken on one side (an insertion, 500 bp, or a truncation), by a
    hand-built AlnMap. side "left": the element joins the left flank only (P: insertion after it; T: a right
    truncation keeping A[s, s + 2000)); "right": the right flank only (P: insertion before it; T: a left
    truncation keeping A[s + 2800, t)). Returns (case, map, B)."""
    c = toy.case("P", random.Random(7))
    [e] = c.elements
    s, t, n = e.start - 1, e.end, len(c.a)
    ins = toy.rand_seq(500, random.Random(8))
    if name == "P" and side == "left":
        b = c.a[:t] + ins + c.a[t:]
        m = amap((0, 0, t), (t, t + 500, n - t))
    elif name == "P":
        b = c.a[:3000] + ins + c.a[3000:]
        m = amap((0, 0, 3000), (3000, 3500, n - 3000))
    elif side == "left":  # keep A[s, s + 2000): the left LTR and 1.6 kb of the internal region
        b = c.a[:s + 2000] + ins + c.a[t:]
        m = amap((0, 0, s + 2000), (t, s + 2500, n - t))
    else:  # keep A[s + 2800, t)
        b = c.a[:3000] + ins + c.a[s + 2800:]
        m = amap((0, 0, 3000), (s + 2800, 3500 + 5, n - s - 2800))
    return c, m, b


@pytest.mark.parametrize("side", ["left", "right"])
@pytest.mark.parametrize("name", ["P", "T"])
def test_anchored_on_one_side_keeps_present(name, side):
    """The element's bases join the flank on one side only (the other junction broken by an insertion or a
    truncation): the call stands, P for the whole element, T for a truncated copy."""
    c, m, b = _one_sided_case(name, side)
    [x] = run_case(c, aln=m, b=b)
    assert x.state == name and x.ev["reason"] is None and x.b_start is not None, (name, side, x)


def test_unanchored_host_is_inherited_by_nested():
    """An element nested in an unanchored host is U with the host's reason and the inherit flag (rule 0)."""
    rng = random.Random(7)
    c = toy.case("nested_E", rng)
    h, ch = sorted(c.elements, key=lambda e: e.level)
    hs, ht, cs, ce, n = h.start - 1, h.end, ch.start - 1, ch.end, len(c.a)
    x1, x2 = toy.rand_seq(500, rng), toy.rand_seq(500, rng)
    host_b = c.a[hs:cs] + c.a[ce + len(toy.T2):ht]  # the host without its child and one T2 copy
    b = c.a[:3000] + x1 + host_b + x2 + c.a[-3000:]
    b0 = 3500
    m = amap((0, 0, 3000), (hs, b0, cs - hs), (ce + len(toy.T2), b0 + cs - hs, ht - ce - len(toy.T2)),
             (n - 3000, len(b) - 3000, 3000))
    calls = {x.eid: x for x in run_case(c, aln=m, b=b)}
    ph, pc = calls[h.eid], calls[ch.eid]
    assert (ph.state, ph.ev["reason"], ph.flags) == ("U", "unanchored", [])
    assert (pc.state, pc.ev["reason"], pc.flags) == ("U", "unanchored", ["inherit"])


@needs_lastz
def test_present_nested_in_anchored_present_host():
    """A present element nested in a present, anchored host is called on its own (rule 0
    passes it on, rule 5a finds its own runs): both P, no inherit flag."""
    for seed in SEEDS:
        c = toy.case("nested_E", random.Random(seed))
        calls = run_case(c, b=c.a)  # B holds the host and its child
        assert [(x.state, x.flags, x.ev["reason"]) for x in calls] == [("P", [], None)] * 2, seed


def _thin_case(k, side):
    """A = L T e T R; B holds the element at its site, its junction on `side` aligned by a run that reaches k bases
    into the element (then a 100-bp B insertion breaks the run), its other junction broken by a 500-bp B insertion
    whose bases beside the element are not a TSD (no TSD exception). Returns (case, map, B)."""
    c = toy.case("P", random.Random(7))
    [e] = c.elements
    s, t, n, a = e.start - 1, e.end, len(c.a), c.a
    rng = random.Random(8)
    x, y = toy.rand_seq(100, rng), "AAAAA" + toy.rand_seq(495, rng)  # "AAAAA": 5 mismatches to the TSD
    if side == "left":
        b = a[:s + k] + x + a[s + k:t] + y + a[t:]
        m = amap((0, 0, s + k), (s + k, s + k + 100, t - s - k), (t, t + 600, n - t))
    else:
        b = a[:s] + y[::-1] + a[s:t - k] + x + a[t - k:]
        m = amap((0, 0, s), (s, s + 500, t - k - s), (t - k, t - k + 600, n - t + k))
    return c, m, b


@pytest.mark.parametrize("side", ["left", "right"])
def test_thin_crossing_is_unanchored(side):
    """A run crossing the junction with its flank bases but fewer than MIN_SEG of the element's own does not join
    the element to its flank (rule 5a): U unanchored, though rule 5 calls P."""
    c, m, b = _thin_case(MIN_SEG - 1, side)
    [x] = run_case(c, aln=m, b=b)
    assert (x.state, x.ev["reason"]) == ("U", "unanchored") and x.ev["cov_elem"] == 1.0, x


@pytest.mark.parametrize("side", ["left", "right"])
def test_crossing_by_min_seg_passes(side):
    c, m, b = _thin_case(MIN_SEG, side)
    [x] = run_case(c, aln=m, b=b)
    assert x.state == "P" and x.ev["reason"] is None and x.b_start is not None, x


def _tsd_exception_case(outside):
    """A = L T e T R; B holds e without its first 100 bases, unaligned on both sides (no run crosses a junction;
    its far side broken by a 100-bp insertion after its TSD copy), with a TSD copy beside the image of its LTR
    ends on both sides (the TSD exception's read). `outside`: B lacks T and e[:100] and holds a TSD-like 5 bp in
    L, so that image starts 100 bp inside the left flank's image (B = L' e[100:] T X R); else B = L T e T X R and
    the image lies between the flanks' images. Returns (case, map, B)."""
    c = toy.case("P", random.Random(7))
    [e] = c.elements
    s, t, n, a = e.start - 1, e.end, len(c.a), c.a
    x = toy.rand_seq(100, random.Random(8))
    if outside:
        b = a[:2895] + TSD + a[2900:3000] + a[s + 100:t] + TSD + x + a[t + 5:]
        b_e = 3000  # B position of e[100:]
    else:
        b = a[:3000] + TSD + a[s:t] + TSD + x + a[t + 5:]
        b_e = 3105
    m = amap((0, 0, 3000), (s + 100, b_e, t - s - 100), (t + 5, b_e + t - s - 100 + 105, n - t - 5))
    return c, m, b


def test_tsd_exception_outside_the_flank_images_is_u():
    """The TSD exception's image (the B images of the LTR ends, extrapolated) must lie between the B images of
    the element's flanks (rule 3's lo and hi): here it starts 100 bp inside the left flank's image, so the TSD
    read there is flank sequence, not a site of the element's own -> U unanchored (rule 5 calls T)."""
    c, m, b = _tsd_exception_case(outside=True)
    [x] = run_case(c, aln=m, b=b)
    assert (x.state, x.ev["reason"]) == ("U", "unanchored"), x
    assert b[2895:2900] == b[7700:7705] == TSD  # the reads beside the image would match


def test_tsd_exception_inside_the_flank_images_passes():
    c, m, b = _tsd_exception_case(outside=False)
    [x] = run_case(c, aln=m, b=b)
    assert x.state == "T" and x.ev["reason"] is None and x.b_start is not None, x


@needs_lastz
def test_anchor_distances_are_each_elements_own():
    """Every element's anchor distances are measured from itself to its window's anchors, a nested element's
    from its own ends (not its host's) and so are an interior element's in a merged window; a window that is not
    ok has none (0, 0)."""
    c = toy.case("nested_E", random.Random(7))
    h, ch = sorted(c.elements, key=lambda e: e.level)
    win = _window(c, left_anchor=100, right_anchor=len(c.a) - 200)
    calls = {x.eid: x for x in run_case(c, win=win)}
    for e in (h, ch):
        assert (calls[e.eid].ev["anchor_left"], calls[e.eid].ev["anchor_right"]) == (
            e.start - 2 - 100, len(c.a) - 200 - e.end)
    assert calls[ch.eid].ev["anchor_left"] > calls[h.eid].ev["anchor_left"]
    bad = replace(win, status="breakpoint")
    assert [(x.ev["anchor_left"], x.ev["anchor_right"]) for x in classify_window(
        bad, None, "", "", {e.eid: e for e in c.elements}, _children(c.elements))] == [(0, 0), (0, 0)]


def test_low_flank_identity_lowers_the_tier():
    """Identity below IDENT_MIN over the element +- IDENT_FLANK means placement is not
    trusted: a P with full coverage drops from high to medium (evidence alone)."""
    c = toy.case("P", random.Random(7))
    n = len(c.a)
    m = amap((0, 0, n, n // 4))  # 75% identity everywhere
    [x] = run_case(c, aln=m, b=c.a)
    assert (x.state, x.tier) == ("P", "medium") and x.ev["ident"] == pytest.approx(0.75, abs=0.01)


@needs_lastz
def test_window_members_are_visited_once():
    """A window listing an element twice, or a children map repeating a child, still calls
    each element once."""
    c = toy.case("nested_E", random.Random(7))
    h, ch = sorted(c.elements, key=lambda e: e.level)
    win = _window(c, eids=[h.eid, h.eid])
    calls = classify_window(win, align_window(c.a, c.b), c.a, c.b, {e.eid: e for e in c.elements},
                            {h.eid: [ch.eid, ch.eid]})
    assert [x.eid for x in calls] == [h.eid, ch.eid]
