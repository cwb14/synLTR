"""LASTZ window alignment (routes.assembly.lastz): the alignment map (AlnMap) and its queries, hits to the map,
the toy genotype cases (tests/builders.py) aligned by LASTZ and called by the classifier over the classifier's
seeds, LASTZ's options and failure handling, and the cleanup of a dying worker's LASTZ child."""
from __future__ import annotations

import gzip
import logging
import math
import os
import random
import signal
import subprocess
import sys
import time
from collections import namedtuple
from dataclasses import replace
from pathlib import Path

import builders
import mappy
import numpy as np
import pytest
from markers import needs_lastz

from synltr.base.constants import MIN_SEG, MM_F_EQX, SCORE, WINDOW_PAD
from synltr.io.elements import Element, make_eid
from synltr.io.samples import Sample
from synltr.io.seqio import Fasta, revcomp
from synltr.routes.assembly import lastz
from synltr.routes.assembly.anchor import PafRoute, PairContext, anchor_pair
from synltr.routes.assembly.classify import classify_window
from synltr.routes.assembly.lastz import AlnMap, align_window, alnmap_from_hits, lastz_args
from synltr.routes.assembly.pair import genotype_pair
from synltr.routes.assembly.window import Window
from synltr.seq.junction import Gap

LOG = "synltr.routes.assembly.lastz"
SEEDS = (7, 0, 4, 6, 9, 11, 14)  # the classifier's seeds (tests/routes/assembly/test_classify.py): chance microhomology
CUT = ("P", "E", "T", "O", "inv_site")  # also called on window-sized flanks (WINDOW_PAD bp: make_windows's cut)
_Hit = namedtuple("_Hit", "q_st q_en r_st r_en strand mlen ops lens")  # the hit fields alnmap_from_hits reads


def Hit(q_st, q_en, r_st, r_en, strand, mlen, cigar):
    """A hit of A (query) on B (reference) with its CIGAR given as (length, op) pairs."""
    return _Hit(q_st, q_en, r_st, r_en, strand, mlen, np.array([op for _, op in cigar], dtype=np.uint8),
                np.array([n for n, _ in cigar], dtype=np.int64))


def amap(*segs, inv=()):
    """AlnMap from (a, b, len[, mm]) tuples."""
    rows = [s if len(s) == 4 else (*s, 0) for s in segs]
    cols = list(zip(*rows)) if rows else [(), (), (), ()]
    return AlnMap(*(np.array(c, dtype=np.int64) for c in cols), inv=list(inv))


def _hamming(x, y):
    return sum(p != q for p, q in zip(x, y))


# ---------------------------------------------------------------- simple pairs


@needs_lastz
def test_identical_one_segment():
    s = builders.rand_seq(20_000, random.Random(1))
    m = align_window(s, s)
    assert m.aligned_bp(0, 20_000) >= 19_950
    assert len(m.seg_a) == 1 and m.identity(0, 20_000) == 1.0 and m.inv == []


@needs_lastz
def test_insertion_makes_gap():
    rng = random.Random(2)
    b = builders.rand_seq(20_000, rng)
    el, _, _ = builders.make_element(rng, 500, 5000)
    a = b[:10_000] + el + b[10_000:]
    g = align_window(a, b).gap_around(13_000)
    assert abs(g.a0 - 10_000) <= 50 and abs(g.a1 - 16_000) <= 50
    assert abs(g.b0 - 10_000) <= 50 and abs(g.b1 - 10_000) <= 50


@needs_lastz
def test_divergence_identity():
    rng = random.Random(3)
    a = builders.rand_seq(20_000, rng)
    b = builders.mutate(a, 0.02, rng)
    assert abs(align_window(a, b).identity(0, len(a)) - 0.98) <= 0.005


@needs_lastz
def test_inversion_recorded():
    rng = random.Random(4)
    a = builders.rand_seq(23_000, rng)
    b = a[:10_000] + revcomp(a[10_000:13_000]) + a[13_000:]
    m = align_window(a, b)
    assert any(a0 <= 10_050 and a1 >= 12_950 and b0 <= 10_050 and b1 >= 12_950 for a0, a1, b0, b1 in m.inv)
    assert m.aligned_bp(10_050, 12_950) == 0  # the inverted piece is not in the collinear segments
    assert m.aligned_bp(0, 10_000) >= 9_950 and m.aligned_bp(13_000, 23_000) >= 9_950


# ---------------------------------------------------------------- AlnMap queries (hand-built maps)


def test_gap_around_runs():
    # run 1 = two segments split by a 1-bp A insertion (a flank indel); 3-kb A insertion; run 2
    m = amap((0, 0, 1000), (1001, 1000, 999), (5000, 1999, 1000))
    assert m.gap_around(3000) == Gap(2000, 5000, 1999, 1999)
    assert m.gap_around(2000) == Gap(2000, 5000, 1999, 1999)  # run 1 ends at 2000 (half-open)
    # inside a run: zero-length gap at a
    assert m.gap_around(500) == Gap(500, 500, 500, 500)
    assert m.gap_around(1000) == Gap(1000, 1000, 1000, 1000)  # on the 1-bp indel: B insertion point
    assert m.gap_around(5000) == Gap(5000, 5000, 1999, 1999)
    assert m.gap_around(5500) == Gap(5500, 5500, 2499, 2499)
    # no run on one side
    assert m.gap_around(6000) is None and m.gap_around(-5) is None
    assert amap().gap_around(10) is None


def test_gap_around_run_joining():
    flank, far = (0, 0, 1000), (3000, 1119, 1000)
    # B gap < MIN_SEG: the 40-bp segment joins the flank run, which then ends at its end
    m = amap(flank, (1000, 1000 + MIN_SEG - 1, 40), far)
    assert m.gap_around(2000) == Gap(1040, 3000, 1119, 1119)
    # B gap == MIN_SEG: the 40-bp segment is its own run, too short to be a flank
    m = amap(flank, (1000, 1000 + MIN_SEG, 40), (3000, 1200, 1000))
    assert m.gap_around(2000) == Gap(1000, 3000, 1000, 1200)
    assert m.gap_around(1010) == Gap(1000, 3000, 1000, 1200)  # inside the short run: not a flank
    # A gap < MIN_SEG joins as well
    m = amap(flank, (1000 + MIN_SEG - 1, 1000, 40), far)
    assert m.gap_around(2000) == Gap(1119, 3000, 1040, 1119)
    # short segments that join into a run of >= MIN_SEG aligned bases make a flank
    m = amap((0, 0, 50), (51, 50, 50), (3000, 100, 1000))
    assert m.gap_around(1000) == Gap(101, 3000, 100, 100)
    # a side holding only a short run has no flank
    m = amap((0, 0, MIN_SEG - 1), (3000, 100, 1000))
    assert m.gap_around(1000) is None


def test_b_gaps_within():
    m = amap((0, 0, 1000), (1000, 1100, 1000), (2005, 2100, 500), (2600, 2700, 400))
    assert m.b_gaps_within(500, 3000) == [Gap(1000, 1000, 1000, 1100), Gap(2505, 2600, 2600, 2700)]
    assert m.b_gaps_within(1000, 1001) == [Gap(1000, 1000, 1000, 1100)]
    assert m.b_gaps_within(0, 1000) == []  # the gap starts at the interval end
    assert m.b_gaps_within(1500, 2500) == []  # only an A insertion (B gap 0) starts inside
    assert amap().b_gaps_within(0, 100) == []


def test_aligned_bp():
    m = amap((0, 0, 1000), (1001, 1000, 999), (5000, 1999, 1000))
    assert m.aligned_bp(0, 6000) == 2999
    assert m.aligned_bp(500, 1500) == 999  # A base 1000 is unaligned
    assert m.aligned_bp(2000, 5000) == 0 and m.aligned_bp(10, 10) == 0 and m.aligned_bp(20, 10) == 0
    assert m.aligned_bp(5999, 7000) == 1 and m.aligned_bp(-100, 1) == 1


def test_identity_prorated():
    m = amap((0, 0, 1000, 10), (1000, 1100, 1000, 0))
    assert m.identity(0, 1000) == pytest.approx(0.99)
    assert m.identity(500, 1500) == pytest.approx(0.995)  # 5 prorated mismatches over 1000 bp
    assert m.identity(0, 5000) == pytest.approx(0.995)
    assert math.isnan(m.identity(3000, 4000)) and math.isnan(amap().identity(0, 10))


def test_alnmap_rejects_bad_segments():
    for segs in ([(0, 0, 1000), (900, 1000, 100)],  # overlap on A
                 [(0, 0, 1000), (1000, 900, 100)],  # overlap on B
                 [(0, 0, 0)],  # empty segment
                 [(0, 0, 10, 11)]):  # more mismatches than bases
        with pytest.raises(ValueError):
            amap(*segs)


# ---------------------------------------------------------------- hits -> map: chain DP and clipping


def test_chain_picks_best_collinear_hits_and_clips():
    rng = random.Random(9)
    a, b = builders.rand_seq(3000, rng), builders.rand_seq(3000, rng)
    h1 = Hit(0, 1000, 0, 1000, 1, 1000, [(1000, 0)])
    hits = [
        Hit(2100, 2600, 900, 1400, 1, 500, [(500, 0)]),  # after h1 only: B goes back after h2
        Hit(500, 1500, 2500, 3000, 1, 400, [(100, 0), (500, 1), (400, 0)]),  # after h1 only; weaker
        Hit(970, 2000, 990, 2040, 1, 1030, [(30, 0), (20, 2), (1000, 0)]),  # 1st M inside h1 on both axes
        Hit(2200, 2400, 100, 300, -1, 200, [(200, 0)]),
        Hit(2200, 2400, 100, 300, -1, 200, [(200, 0)]),  # mappy can repeat an inverted hit as secondary
        h1,
    ]
    m = alnmap_from_hits(hits, a, b)
    assert list(zip(m.seg_a.tolist(), m.seg_b.tolist(), m.seg_len.tolist())) == [(0, 0, 1000), (1000, 1040, 1000)]
    assert m.seg_mm.tolist() == [_hamming(a[:1000], b[:1000]), _hamming(a[1000:2000], b[1040:2040])]
    assert m.inv == [(2200, 2400, 100, 300)]
    # partial overlap: the second hit is clipped on the diagonal, not dropped
    m = alnmap_from_hits([h1, Hit(950, 2000, 1050, 2100, 1, 1000, [(1050, 0)])], a, b)
    assert list(zip(m.seg_a.tolist(), m.seg_b.tolist(), m.seg_len.tolist())) == [(0, 0, 1000), (1000, 1100, 1000)]
    # hits that do not advance on every axis cannot chain: the higher-mlen one wins
    m = alnmap_from_hits([h1, Hit(0, 1500, 1200, 2700, 1, 1400, [(1500, 0)])], a, b)
    assert list(zip(m.seg_a.tolist(), m.seg_b.tolist(), m.seg_len.tolist())) == [(0, 1200, 1500)]


def test_hits_cigar_ops():
    s = "ACGTACGTAC" * 3
    t = s[:10] + "T" + s[11:]  # one substitution at 10
    m = alnmap_from_hits([Hit(0, 30, 0, 30, 1, 29, [(10, 7), (1, 8), (19, 7)])], s, t)
    assert m.seg_len.tolist() == [30] and m.seg_mm.tolist() == [1]  # =/X runs merge into one segment
    m = alnmap_from_hits([Hit(0, 30, 0, 29, 1, 29, [(10, 0), (1, 1), (19, 0)])], s, s[:10] + s[11:])
    assert list(zip(m.seg_a.tolist(), m.seg_b.tolist(), m.seg_len.tolist())) == [(0, 0, 10), (11, 10, 19)]
    with pytest.raises(ValueError, match="unsupported CIGAR op 4 in the hit at A 0-10"):  # S: soft clip
        alnmap_from_hits([Hit(0, 10, 0, 10, 1, 10, [(10, 4)])], s, s)
    with pytest.raises(ValueError, match="ends at A 10, B 10"):  # the CIGAR spans less than the hit
        alnmap_from_hits([Hit(0, 12, 0, 12, 1, 10, [(10, 7)])], s, s)
    m = alnmap_from_hits([Hit(0, 20, 0, 25, 1, 20, [(10, 7), (5, 3), (10, 7)])], s, s[:10] + "AAAAA" + s[10:])
    assert list(zip(m.seg_a.tolist(), m.seg_b.tolist())) == [(0, 0), (10, 15)]  # N: B only, like D
    assert len(alnmap_from_hits([], s, s).seg_a) == 0


def test_n_counts_as_mismatch():
    s = "ACGTACGTAC" * 3
    t = s[:5] + "NN" + s[7:]
    m = alnmap_from_hits([Hit(0, 30, 0, 30, 1, 28, [(30, 0)])], s.lower(), t)
    assert m.seg_mm.tolist() == [2]


# ---------------------------------------------------------------- the toy cases, called on LASTZ alignments


def calls(c: builders.ToyCase) -> dict:
    """eid -> Call for toy case `c` aligned by LASTZ (asm20) over one window spanning the whole case."""
    win = Window(wid=0, a_seqid="a", a_start=0, a_end=len(c.a), b_seqid="b", b_start=0, b_end=len(c.b),
                 strand=1, eids=[e.eid for e in c.elements if e.parent is None], status="ok", left_anchor=0,
                 right_anchor=len(c.a) - 1, margin=1.0)
    children: dict[str, list[str]] = {}
    for e in c.elements:
        if e.parent is not None:
            children.setdefault(e.parent, []).append(e.eid)
    aln = align_window(c.a, c.b)
    return {x.eid: x for x in classify_window(win, aln, c.a, c.b, {e.eid: e for e in c.elements}, children)}


def cut(c: builders.ToyCase) -> builders.ToyCase:
    """The case with WINDOW_PAD bp of flank kept on each side of A and B (elements shifted; eids kept)."""
    k = builders.FLANK - WINDOW_PAD
    els = [replace(e, start=e.start - k, end=e.end - k, ltr_left=(e.ltr_left[0] - k, e.ltr_left[1] - k),
                   ltr_right=(e.ltr_right[0] - k, e.ltr_right[1] - k)) for e in c.elements]
    return replace(c, a=c.a[k:len(c.a) - k], b=c.b[k:len(c.b) - k], elements=els)


def check(c: builders.ToyCase, got: dict, where: tuple) -> None:
    """Every element's expected state, required and forbidden flags, tier and U reason."""
    for eid, state in c.expect.items():
        x = got[eid]
        assert x.state == state, (where, x)
        flags = set(x.flags)
        assert c.flags.get(eid, set()) <= flags and not c.noflags.get(eid, set()) & flags, (where, x)
        assert c.tiers.get(eid, x.tier) == x.tier, (where, x)
        assert c.reasons.get(eid, x.ev["reason"]) == x.ev["reason"], (where, x)


@needs_lastz
@pytest.mark.parametrize("name", builders.CASES)
def test_toy_case(name):
    """Each of the 21 toy cases, on every seed, over a window spanning it (3-kb flanks): the expected calls,
    P, E, D_left and D_both among them."""
    for seed in SEEDS:
        c = builders.case(name, random.Random(seed))
        check(c, calls(c), (name, seed))


@needs_lastz
@pytest.mark.parametrize("name", CUT)
def test_toy_case_on_window_sized_flanks(name):
    """P, E, T, O and inv_site with only WINDOW_PAD bp of flank on each side, as make_windows cuts a window: both
    flanks still seed and align around the element-sized indel, and the chained minus strand still reports the
    inverted site."""
    for seed in SEEDS:
        c = cut(builders.case(name, random.Random(seed)))
        check(c, calls(c), (name, seed))


# ---------------------------------------------------------------- close genome pairs


def close_pair(rng: random.Random, sites=("E", "P", "E", "P"), spacer=300_000, patch=1500, rate=0.05):
    """A close pair: per site L T e T R (3-kb
    flanks, a 5-kb element, a 5-bp TSD) between identical random spacers; B is A with each site's element and
    `patch` bp on each side (TSDs included) substituted at `rate`, and an E site's element and right TSD copy
    deleted. Returns (A, B, A's elements, eid -> true state). The spacers dominate the anchored divergence d
    (about 0.0007), though every site's window is 5% divergent: an introgressed tract inside a near-identical
    pair."""
    def seq(n):
        return "".join(rng.choices("ACGT", k=n))

    a_parts, b_parts, elements, truth = [seq(spacer)], [], [], {}
    b_parts.append(a_parts[0])
    for state in sites:
        left, right, tsd = seq(builders.FLANK), seq(builders.FLANK), seq(5)
        e, ltr, _ = builders.make_element(rng, 400, 4200)
        start = sum(map(len, a_parts)) + len(left) + len(tsd) + 1  # 1-based
        el = Element("A", make_eid("A", "chr1", start, start + len(e) - 1), f"{state}{len(elements)}", "chr1",
                     start, start + len(e) - 1, "+", (start, start + ltr - 1),
                     (start + len(e) - ltr, start + len(e) - 1), tsd, 0.0, None, None, None, None)
        elements.append(el)
        truth[el.eid] = state
        a_parts.append(left + tsd + e + tsd + right)
        kept = tsd + e + tsd if state == "P" else tsd  # B = L T R at an empty site
        b_parts.append(left[:-patch] + builders.mutate(left[-patch:] + kept + right[:patch], rate, rng) + right[patch:])
        a_parts.append(seq(spacer))
        b_parts.append(a_parts[-1])
    return "".join(a_parts), "".join(b_parts), elements, truth


def _whole_genome_paf(a_path: Path, b_path: Path, out: Path) -> None:
    """A's whole-genome alignment to B as out/A__B.paf: every mappy asm20 hit, with its =/X CIGAR."""
    aligner = mappy.Aligner(str(b_path), preset="asm20", extra_flags=MM_F_EQX)
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "A__B.paf", "w") as fh:
        for name, seq, _ in mappy.fastx_read(str(a_path)):
            for h in aligner.map(seq):
                fields = (name, len(seq), h.q_st, h.q_en, "+" if h.strand > 0 else "-", h.ctg, h.ctg_len, h.r_st,
                          h.r_en, h.mlen, h.blen, h.mapq)
                fh.write("\t".join(map(str, fields)) + f"\tcg:Z:{h.cigar_str}\n")


@needs_lastz
def test_close_pair_with_divergent_patches(tmp_path):
    """A near-identical pair (d <= 0.001 over its whole-genome anchors) whose sites sit in 5%-divergent patches:
    each element is called right, E or P (under asm5 scoring LASTZ fragments these windows' alignments into D, O,
    T and S calls, some at high tier). LASTZ scores every window with asm20's scoring, whatever the pair's
    divergence."""
    a, b, elements, truth = close_pair(random.Random(1))
    d = tmp_path
    fa, fb = d / "A.fa", d / "B.fa"
    for path, s in ((fa, a), (fb, b)):
        builders.write_fasta(path, {"chr1": s})
    _whole_genome_paf(fa, fb, d / "paf")
    route = PafRoute(d / "paf")
    with Fasta(fa) as A, Fasta(fb) as B:
        ctx = PairContext(Sample("A", fa, None, None, None), Sample("B", fb, None, None, None), A, B, elements, None,
                          d / "out", 1)
        div = anchor_pair(ctx, route)[3]
        out = genotype_pair(ctx, route)
    assert div <= 0.001, div
    with gzip.open(out, "rt") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        got = {r["eid"]: r for r in (dict(zip(head, line.rstrip("\n").split("\t"))) for line in fh)}
    assert {eid: got[eid]["state"] for eid in truth} == truth, [(got[e]["state"], got[e]["tier"]) for e in truth]


# ---------------------------------------------------------------- options, defaults and failures


def test_lastz_args_are_asm20_for_every_pair():
    """asm20's scoring (the engine's SCORE), ydrop = the cost of a MIN_SEG-bp gap (166), gapped threshold MIN_SEG,
    every position seeded with the 12of19 transition seed, chained, IUPAC codes accepted: one option set for
    every pair."""
    match, mismatch, gap_open, gap_extend = SCORE
    assert SCORE == (1, 4, 6, 2)  # minimap2 asm20
    args = lastz_args()
    assert args[:5] == [f"--match={match},{mismatch}", f"--gap={gap_open},{gap_extend}", "--ydrop=166",
                        f"--gappedthresh={MIN_SEG}", "--step=1"], args
    assert gap_open + gap_extend * MIN_SEG == 166
    assert {"--seed=12of19", "--transition", "--chain", "--gapped", "--ambiguous=iupac"} <= set(args)
    assert not any(x.startswith("--strand") for x in args)  # set per run


def test_both_strands_run_with_lastz_args(tmp_path, monkeypatch):
    """Both strands run with lastz_args(), asm20's scoring, whatever the pair."""
    log = tmp_path / "argv.txt"
    monkeypatch.setenv("PATH", stub(tmp_path, f'echo "$@" >> {log}; exit 0'))
    align_window("ACGT" * 50, "ACGT" * 50)
    runs = [line.split()[3:] for line in log.read_text().splitlines()]  # after target, query and strand
    assert runs == [lastz_args(), lastz_args()]


@needs_lastz
def test_ambiguous_bases_align_through():
    """IUPAC codes and an N run in B are ambiguous bases (--ambiguous=iupac), never a failed window: the
    sequence aligns through them and the N bases count as mismatches."""
    rng = random.Random(5)
    a = builders.rand_seq(6000, rng)
    b = a[:2000] + "R" + a[2001:3000] + "N" * 20 + a[3020:4000] + "y" + a[4001:]
    m = align_window(a, b)
    assert m.aligned_bp(0, 6000) >= 5900 and m.seg_mm.sum() >= 20


def stub(tmp_path, script: str | None):
    """PATH holding only a `lastz` running `script` (none if None)."""
    d = tmp_path / "bin"
    d.mkdir()
    if script is not None:
        exe = d / "lastz"
        exe.write_text(f"#!/bin/sh\n{script}\n")
        exe.chmod(0o755)
    return str(d)


def test_each_strand_runs_chained_on_its_own(tmp_path, monkeypatch):
    """Two LASTZ runs per window, plus then minus, both chained: the minus strand's chain is its own, so its
    hits (the inverted copies) reach `inv` (unchained, the minus strand of a
    satellite window runs for minutes)."""
    log = tmp_path / "argv.txt"
    monkeypatch.setenv("PATH", stub(tmp_path, f'echo "$@" >> {log}; '
                                              r"[ $3 = --strand=minus ] && printf '0\t50\t0\t50\t-\t50\t50=\n'"
                                              "; exit 0"))
    m = align_window("ACGT" * 50, "ACGT" * 50)
    runs = [line.split() for line in log.read_text().splitlines()]
    assert [r[2] for r in runs] == ["--strand=plus", "--strand=minus"]
    assert all("--chain" in r for r in runs) and m.inv == [(0, 50, 0, 50)] and len(m.seg_a) == 0


def test_lastz_error_leaves_the_window_unaligned(tmp_path, monkeypatch, caplog):
    """A window LASTZ fails on gives an empty AlnMap (its elements are U) and one warning, not a failed
    pair."""
    monkeypatch.setenv("PATH", stub(tmp_path, "echo 'bad sequence' >&2; exit 3"))
    with caplog.at_level(logging.WARNING, logger=LOG):
        m = align_window("ACGT" * 50, "ACGT" * 50)
    assert len(m.seg_a) == 0 and m.inv == []
    [rec] = caplog.records
    assert "left unaligned" in rec.getMessage() and "lastz exit 3: bad sequence" in rec.getMessage()


def test_lastz_malformed_cigar_leaves_the_window_unaligned(tmp_path, monkeypatch, caplog):
    monkeypatch.setenv("PATH", stub(tmp_path, r"printf '0\t10\t0\t10\t+\t10\t10Q\n'"))
    with caplog.at_level(logging.WARNING, logger=LOG):
        m = align_window("ACGT" * 50, "ACGT" * 50)
    assert len(m.seg_a) == 0 and "malformed LASTZ cigar" in caplog.text


def test_lastz_missing_fails(tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", stub(tmp_path, None))
    with pytest.raises(RuntimeError, match="needs lastz on PATH"):
        align_window("ACGT", "ACGT")
    assert len(align_window("", "ACGT").seg_a) == 0  # nothing to align: no binary needed


# ---------------------------------------------------------------- temp folders and children of dying workers


def test_sweep_removes_only_dead_processes_folders(tmp_path, monkeypatch):
    """A worker killed outright leaves its synltr_lastz.<pid>.* folder: the sweep removes those of dead pids
    only, never a live process's, an unrelated name's (such as synltr_lastz_<random>, which a
    running job may still use)."""
    monkeypatch.setattr(lastz, "_SHM", str(tmp_path))
    gone = subprocess.run(["true"], check=True)  # a pid that no longer runs
    dead = subprocess.Popen(["true"])
    dead.wait()
    names = [f"synltr_lastz.{dead.pid}.x1", f"synltr_lastz.{os.getpid()}.x2", "synltr_lastz_12345_ab",
             "synltr_lastz.notapid.x3", "other.1.x"]
    for name in names:
        (tmp_path / name).mkdir()
        (tmp_path / name / "a.fa").write_text(">a\nACGT\n")
    assert gone.returncode == 0 and lastz.sweep_lastz_tmp() == 1
    assert sorted(p.name for p in tmp_path.iterdir()) == sorted(names[1:])


_DYING_WORKER = r"""
import os, sys
from synltr.base import runtime
from synltr.routes.assembly import lastz
lastz._SHM = sys.argv[1]
runtime._init_worker(os.getppid(), lastz.kill_lastz)  # a pool worker's setup: SIGTERM runs kill_lastz, then exits
lastz.align_window("ACGT" * 50, "ACGT" * 50)
"""


@pytest.mark.skipif(not sys.platform.startswith("linux"), reason="lists processes through /proc")
def test_sigterm_kills_the_lastz_child_and_removes_its_folder(tmp_path):
    """A worker killed by SIGTERM (its parent died, or the pool terminated it) kills its
    running lastz child and removes its temp folder before it exits."""
    pidfile, shm = tmp_path / "lastz.pid", tmp_path / "shm"
    shm.mkdir()
    bin_dir = stub(tmp_path, f"echo $$ > {pidfile}; exec sleep 60")
    env = {**os.environ, "PATH": f"{bin_dir}{os.pathsep}{os.environ['PATH']}"}
    proc = subprocess.Popen([sys.executable, "-c", _DYING_WORKER, str(shm)], env=env)
    try:
        t0 = time.monotonic()
        while not (pidfile.exists() and pidfile.read_text().strip()):
            assert proc.poll() is None and time.monotonic() - t0 < 30
            time.sleep(0.05)
        child = int(pidfile.read_text())
        [folder] = shm.iterdir()
        assert folder.name.startswith(f"synltr_lastz.{proc.pid}.")
        proc.terminate()
        assert proc.wait(timeout=10) == 128 + signal.SIGTERM
        t0 = time.monotonic()
        while Path(f"/proc/{child}").exists() and time.monotonic() - t0 < 5:
            time.sleep(0.05)
        assert not Path(f"/proc/{child}").exists() and not folder.exists()
    finally:
        if proc.poll() is None:
            proc.kill()
