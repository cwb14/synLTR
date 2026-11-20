"""One-pass genome builder (ltrtrace.sim.pop.build): sample sequences, per-element copy records, lift maps and assembly
gaps, on exact toy genealogies (poptoy.toy_ts) with hand-placed events and exact expected sequences."""
from __future__ import annotations

import pytest

pytest.importorskip("msprime")

import math
import random
from dataclasses import replace

import numpy as np
from sim import poptoy
from sim.poptoy import DONORS, EXACT, build, events, ins_at, loss_at, other, toy7

from ltrtrace.base.constants import SIM_GAP_LEN, SIM_GAP_NEAR, SLOP_TSD_MULT, TSD_DEFAULT
from ltrtrace.sim.pop.build import Copy, PopGenome, _Chrom, add_gaps, build_chrom, lift_pos
from ltrtrace.sim.pop.donors import Donor, jitter_bounds
from ltrtrace.sim.pop.events import Ins, donor_body, draw_events
from ltrtrace.sim.pop.genome import revcomp_bytes
from ltrtrace.sim.pop.model import Rates, ancestry

L = poptoy.L
W = SLOP_TSD_MULT * TSD_DEFAULT
T5 = TSD_DEFAULT
def rng(seed: int) -> np.random.Generator:
    return np.random.default_rng(seed)


def expected_body(x: Ins, kids: dict, s: int, donors=DONORS) -> bytes:
    """x's body in sample s with everything nested in it there, in x's own (donor) order, without substitutions:
    each nested payload s carries (its body in x's order, then its second TSD copy, x's bases before it) spliced in
    at its offset, the youngest nearest its target at one offset."""
    body = donor_body(donors[x.donor])
    out = bytearray(body)
    for y in sorted((y for y in kids.get(x.iid, []) if s in y.carriers), key=lambda y: (-y.pos, -y.time)):
        inner = expected_body(y, kids, s, donors)
        out[y.pos:y.pos] = (inner if y.strand == "+" else revcomp_bytes(inner)) + body[y.pos - T5:y.pos]
    return bytes(out)


def genome_image(x: Ins, ins: dict, kids: dict, s: int, donors=DONORS) -> bytes:
    """x's expected body in sample s, in genome orientation (the strands of x and its host chain)."""
    flips, y = 0, x
    while y is not None:
        flips += y.strand == "-"
        y = ins[y.host] if y.host is not None else None
    b = expected_body(x, kids, s, donors)
    return b if flips % 2 == 0 else revcomp_bytes(b)


def kids_of(xs) -> dict:
    out: dict = {}
    for x in xs:
        if x.host is not None:
            out.setdefault(x.host, []).append(x)
    return out


# ------------------------------------------------------------------------------------------------ sequence and TSDs

def test_tsd_both_sides():
    ts, tree, n = toy7()
    bb, _ = poptoy.toy_backbone(L, seed=1)
    tsd = bb[99_995:100_000]
    for strand in "+-":
        x = ins_at(tree, 1, None, 100_000, n["01"], 4_000.0, strand=strand)
        body = donor_body(DONORS[0])
        body = body if strand == "+" else revcomp_bytes(body)
        got = build(ts, bb, [], events(bb, 7, [x]))
        for s in (0, 1):
            seq, copies, _ = got[s]
            c = copies[x.iid]
            assert c.fate == "intact" and c.inherited and c.site == -1 and c.tsd == tsd.decode()
            assert c.strand == strand and c.target == (0, None, 100_000) and c.host is None
            assert (c.start, c.end) == (100_000, 100_000 + len(body))
            assert seq[c.start - T5:c.start] == tsd and seq[c.start:c.end] == body and seq[c.end:c.end + T5] == tsd
            assert seq == bb[:100_000] + body + tsd + bb[100_000:]
        for s in range(2, 7):
            seq, copies, _ = got[s]
            c = copies[x.iid]
            assert c.fate == "absent" and not c.inherited and (c.start, c.end) == (-1, -1) and c.site == 99_995
            assert seq == bb and seq[c.site:c.site + T5] == tsd


def test_tsd_at_insertion_time():
    # the second TSD copy holds the target bases at the insertion time: a SNP above the insertion (its derived samples
    # include the carriers) is in both copies; one below it, or one the carriers do not share, only in the first
    ts, tree, n = toy7()
    bb, _ = poptoy.toy_backbone(L, seed=19)
    p = 100_000
    x = ins_at(tree, 1, None, p, n["01"], 4_000.0)  # carriers 0, 1; branch from 12,000 down to 3,000
    snps = {p - 5: {s: other(bb[p - 5]) for s in (0, 1, 2, 3)}, p - 4: {0: other(bb[p - 4])},
            p - 2: {1: other(bb[p - 2]), 4: other(bb[p - 2])}}
    got = build(ts, bb, [], events(bb, 7, [x], snps=snps))
    at = bytearray(bb[p - 5:p])
    at[0] = ord(other(bb[p - 5]))
    for s in (0, 1):
        seq, copies, _ = got[s]
        c = copies[x.iid]
        left = bytearray(at)
        if s == 0:
            left[1] = ord(other(bb[p - 4]))
        else:
            left[3] = ord(other(bb[p - 2]))
        assert c.tsd == at.decode() and seq[c.end:c.end + T5] == at and seq[c.start - T5:c.start] == left
    # a SNP derived in exactly the carriers lies on the insertion's own branch: above it with probability
    # (12,000 - 4,000) / (12,000 - 3,000) = 8/9
    xs = [ins_at(tree, k + 1, None, 1_000 + 600 * k, n["01"], 4_000.0) for k in range(300)]
    snps = {x.pos - 5: {0: other(bb[x.pos - 5]), 1: other(bb[x.pos - 5])} for x in xs}
    copies = build(ts, bb, [], events(bb, 7, xs, snps=snps), seed=3)[0][1]
    above = sum(copies[x.iid].tsd[0] == other(bb[x.pos - 5]) for x in xs)
    assert sum(copies[x.iid].tsd[0] in (other(bb[x.pos - 5]), chr(bb[x.pos - 5])) for x in xs) == 300
    assert abs(above / 300 - 8 / 9) <= 3 * math.sqrt(8 / 9 * 1 / 9 / 300), above


def test_identical_ltrs_and_k2p():
    ts, tree, n = toy7()
    bb, _ = poptoy.toy_backbone(L, seed=2)
    xs = [ins_at(tree, k + 1, None, 20_000 + 15_000 * k, n["0123"], 15_000.0, strand="+-"[k % 2], donor=k % 3)
          for k in range(8)]
    got = build(ts, bb, [], events(bb, 7, xs), Rates(mu=0.0), seed=3)  # copy noise 0.005, before the LTR copy
    for s in range(4):
        seq, copies, _ = got[s]
        for x in xs:
            c = copies[x.iid]
            assert c.k2p == 0.0 and c.ltr_own == c.ltr_kept
            assert seq[slice(*c.ltr_left)] == seq[slice(*c.ltr_right)]
    # mu = 1e-6 per generation, insertions 10,000 generations old on the two leaf branches of a 2-leaf tree: the two
    # LTRs diverge independently, so the expected K2P distance is 2 mu age = 0.02
    ts2 = poptoy.toy_ts((20_000.0, 0, 1), L)
    tree2 = ts2.first()
    long = [Donor(poptoy.rand_bases(3_000, random.Random(k)).decode(), 1_000, 1_000) for k in range(2)]
    xs = [ins_at(tree2, k + 1, None, 1_000 + 2_400 * k, k % 2, 10_000.0, donor=k % 2) for k in range(80)]
    got = build(ts2, bb, [], events(bb, 2, xs), Rates(mu=1e-6), seed=4, donors=long)
    ks = np.array([got[k % 2][1][x.iid].k2p for k, x in enumerate(xs)])
    se = ks.std(ddof=1) / math.sqrt(len(ks))
    assert abs(ks.mean() - 0.02) <= 3 * se, (ks.mean(), se)


def test_orthologous_copies_share_mutations():
    # an insertion 15,000 generations old above the 5,000-generation MRCA of samples 0 and 1: each copy is 15,000
    # generations from the inserted body, the two copies 2 x 5,000 from each other (shared branch, shared mutations)
    ts = poptoy.toy_ts((20_000.0, (5_000.0, 0, 1), 2), L)
    tree = ts.first()
    bb, _ = poptoy.toy_backbone(L, seed=3)
    xs = [ins_at(tree, k + 1, None, 5_000 + 9_000 * k, tree.parent(0), 15_000.0, donor=k % 3) for k in range(20)]
    mu = 1e-6
    got = build(ts, bb, [], events(bb, 3, xs), Rates(mu=mu, copy_noise=0.0), seed=5)
    to_ref = between = total = 0
    for x in xs:
        a, b = (got[s][0][got[s][1][x.iid].start:got[s][1][x.iid].end] for s in (0, 1))
        ref = donor_body(DONORS[x.donor])
        to_ref += sum(u != v for u, v in zip(a, ref))
        between += sum(u != v for u, v in zip(a, b))
        total += len(ref)

    def jc(t):
        return 0.75 * -math.expm1(-4 / 3 * mu * t)

    for obs, t in ((to_ref, 15_000), (between, 10_000)):
        p = jc(t)
        assert abs(obs / total - p) <= 3 * math.sqrt(p * (1 - p) / total), (obs / total, p)


def test_copy_noise():
    ts, tree, n = toy7()
    bb, _ = poptoy.toy_backbone(L, seed=4)
    xs = [ins_at(tree, k + 1, None, 10_000 + 9_000 * k, n["01"], 4_000.0, donor=k % 3) for k in range(20)]
    got = build(ts, bb, [], events(bb, 7, xs), Rates(mu=0.0, copy_noise=0.005), seed=6)
    diffs = total = 0
    for x in xs:
        (s0, c0), (s1, c1) = ((got[s][0], got[s][1][x.iid]) for s in (0, 1))
        body = s0[c0.start:c0.end]
        assert body == s1[c1.start:c1.end]  # one draw per insertion, shared by its carriers
        ltr = DONORS[x.donor].ltr_left
        ref = donor_body(DONORS[x.donor])
        assert body[:ltr] == body[-ltr:] and c0.k2p == 0.0  # the left LTR copied over the right
        m = len(ref) - ltr  # the left LTR and the internal region carry the draws
        diffs += sum(u != v for u, v in zip(body[:m], ref[:m]))
        total += m
    se = math.sqrt(0.005 * 0.995 / total)
    assert abs(diffs / total - 0.005) <= 3 * se, diffs / total


def test_snp_bases_applied():
    ts, _, _ = toy7()
    bb, _ = poptoy.toy_backbone(L, seed=6)
    ev = draw_events(ts, bb, [], [], Rates(lam=0, delta=0, mu=2e-7, indel=0), rng(6))
    assert len(ev.snp_pos) > 1_000
    ref = np.frombuffer(bb, dtype=np.uint8)
    for s, (seq, copies, _) in enumerate(build(ts, bb, [], ev)):
        a = np.frombuffer(seq, dtype=np.uint8)
        assert len(a) == L and not copies
        assert np.array_equal(a[ev.snp_pos], ev.snp_base[:, s])
        other_sites = np.ones(L, dtype=bool)
        other_sites[ev.snp_pos] = False
        assert np.array_equal(a[other_sites], ref[other_sites])
    derived = ev.snp_base != ref[ev.snp_pos][:, None]
    assert derived.any() and not derived.all()


def test_copy_coordinates_hold_payload():
    ts, tree, n = toy7()
    bb, _ = poptoy.toy_backbone(L, seed=5)
    x = ins_at(tree, 1, None, 60_000, n["0123"], 35_000.0, strand="+", donor=0)
    y = ins_at(tree, 2, x.iid, 700, n["0123"], 30_000.0, depth=1, strand="-", donor=1)
    z = ins_at(tree, 3, x.iid, 2_000, n["01"], 6_000.0, depth=1, strand="+", donor=2)
    w = ins_at(tree, 4, y.iid, 300, n["0123"], 20_000.0, depth=2, strand="+", donor=0)
    x2 = ins_at(tree, 5, None, 120_000, n["456"], 25_000.0, strand="-", donor=1)
    old = ins_at(tree, 8, None, 150_000, 6, 15_000.0, strand="+", donor=2)  # sample 6: two at one position
    young = ins_at(tree, 9, None, 150_000, 6, 5_000.0, strand="-", donor=1)
    v = ins_at(tree, 6, x2.iid, 900, n["45"], 10_000.0, depth=1, strand="+", donor=2)
    u = ins_at(tree, 7, x2.iid, 900, 4, 1_000.0, depth=1, strand="-", donor=0)  # v's offset, younger: first
    xs = [x, y, z, w, x2, v, u, old, young]
    indels = [(59_000, 4, False, {0, 1, 2, 3}), (59_500, 3, True, {0}), (119_990, 2, False, {4})]
    got = build(ts, bb, [], events(bb, 7, xs, indels=indels))
    ins, kids = {i.iid: i for i in xs}, kids_of(xs)
    for s, (seq, copies, _) in enumerate(got):
        for i in xs:
            c = copies[i.iid]
            if c.fate == "intact":
                assert seq[c.start:c.end] == genome_image(i, ins, kids, s), (s, i.iid)
    seq, copies, _ = got[0]
    assert copies[x.iid].start == 60_000 - 4 + 3
    assert [copies[i.iid].fate for i in xs] == ["intact"] * 4 + ["absent"] * 5
    assert [copies[i.iid].strand for i in (x, y, z, w)] == ["+", "-", "+", "-"]
    assert seq[:copies[x.iid].start - T5] == bb[:59_000] + bb[59_004:59_500] + seq[59_496:59_499] + bb[59_500:59_995]
    _, copies4, _ = got[4]
    assert copies4[x2.iid].start == 120_000 - 2 and copies4[x2.iid].strand == "-"
    assert [copies4[i.iid].strand for i in (v, u)] == ["-", "+"]
    # absent nested copies: the site lies in their host's image (host present, window intact)
    seq2, copies2, _ = got[2]
    z2 = copies2[z.iid]
    assert z2.fate == "absent" and z2.site >= 0 and seq2[z2.site:z2.site + T5] == donor_body(DONORS[0])[1_995:2_000]
    seq6, copies6, _ = got[6]
    pay_old = donor_body(DONORS[2]) + bb[149_995:150_000]
    pay_young = revcomp_bytes(donor_body(DONORS[1])) + bb[149_995:150_000]
    pay_x2 = revcomp_bytes(donor_body(DONORS[1])) + bb[119_995:120_000]  # sample 6 carries x2 alone
    assert seq6 == bb[:120_000] + pay_x2 + bb[120_000:150_000] + pay_young + pay_old + bb[150_000:]
    assert copies6[young.iid].start == 150_000 + len(pay_x2)
    assert copies6[old.iid].start == copies6[young.iid].start + len(pay_young)
    seq5, copies5, _ = got[5]
    u5 = copies5[u.iid]
    assert u5.fate == "absent" and seq5[u5.site:u5.site + T5] == revcomp_bytes(donor_body(DONORS[1])[895:900])
    assert u5.tsd == revcomp_bytes(donor_body(DONORS[1])[895:900]).decode()  # in genome orientation
    # every intact copy of a drawn landscape with nesting (no losses), with small indels in the backbone
    ev = draw_events(ts, bb, [], DONORS, Rates(lam=5e-9, delta=0, mu=0, indel=5e-9, copy_noise=0), rng(7))
    assert sum(i.depth > 0 for i in ev.ins) > 10 and len(ev.indel_pos) > 20
    ins, kids = {i.iid: i for i in ev.ins}, kids_of(ev.ins)
    checked = 0
    for s, (seq, copies, _) in enumerate(build(ts, bb, [], ev, EXACT)):
        for i in ev.ins:
            c = copies[i.iid]
            assert c.fate in ("intact", "absent")
            if c.fate == "intact":
                assert seq[c.start:c.end] == genome_image(i, ins, kids, s)
                checked += 1
    assert checked > 100


def test_lift_pos():
    ts, tree, n = toy7()
    bb, els = poptoy.toy_backbone(L, [(50_000, 300, 3_000), (150_000, 300, 4_000)], seed=12)
    e2, e = els  # [50,000, 53,600) and [150,000, 154,600)
    x = ins_at(tree, 1, None, 100_000, n["0123"], 30_000.0, donor=0)  # payload 2,600 + 5
    losses = [loss_at(tree, e2.eid, n["01"], 10_000.0, "replace", 500),
              loss_at(tree, e.eid, n["01"], 10_000.0, "del_flank", 200, 300)]  # removes [149,800, 154,900)
    indels = [(70_000, 6, False, {0, 1}), (80_000, 4, True, {0})]
    got = build(ts, bb, els, events(bb, 7, [x], losses, indels=indels))
    seq, _, lift = got[0]
    rep, pay = 500 - 3_600, 2_600 + T5
    expect = {0: 0, 49_999: 49_999, 50_000: 50_500, 53_599: 50_500, 53_600: 50_500, 69_999: 69_999 + rep,
              70_000: 70_000 + rep, 70_005: 70_000 + rep, 70_006: 70_000 + rep, 79_999: 79_993 + rep,
              80_000: 79_998 + rep, 99_999: 99_997 + rep, 100_000: 99_998 + rep + pay,
              149_799: 149_797 + rep + pay, 149_800: 149_798 + rep + pay, 154_899: 149_798 + rep + pay,
              154_900: 149_798 + rep + pay, L: len(seq)}
    assert {p: lift_pos(lift, p) for p in expect} == expect
    assert len(seq) == L + rep - 6 + 4 + pay - 5_100
    kept = [p for p in range(0, L, 997) if not (50_000 <= p < 53_600 or 70_000 <= p < 70_006
                                                  or 149_800 <= p < 154_900)]
    assert all(seq[lift_pos(lift, p)] == bb[p] for p in kept)
    seq4, _, lift4 = got[4]
    assert seq4 == bb and all(lift_pos(lift4, p) == p for p in range(0, L + 1, 1_009))
    with pytest.raises(ValueError):
        lift_pos(lift, L + 1)


def test_site_window_ignores_other_insertion():
    # x's empty-site window is [p - 15, p + 10): the target and w = 10 bases each side. Only removals of its bases
    # break it; SNPs, small indels and insertions at other sites inside it do not
    ts, tree, n = toy7()
    bb, els = poptoy.toy_backbone(L, [(97_000, 300, 2_000), (100_020, 300, 2_000)], seed=13)
    e0, e = els  # [97,000, 99,600) and [100,020, 102,620)
    p = 100_000
    x = ins_at(tree, 1, None, p, 0, 1_000.0)  # sample 0 only
    w = ins_at(tree, 2, None, p + 4, 4, 1_000.0)  # sample 4: 4 bp past the target
    w2 = ins_at(tree, 3, None, p - 4, 5, 1_000.0)  # sample 5: inside the target bases
    losses = [loss_at(tree, e.eid, n["23"], 6_000.0, "del_flank", 11, 0),  # [100,009, ...): the window's last base
              loss_at(tree, e.eid, 6, 2_000.0, "del_flank", 10, 0),  # [100,010, ...): just past the window
              loss_at(tree, e0.eid, 1, 2_000.0, "del_flank", 0, 386),  # [..., 99,986): the window's first base
              loss_at(tree, e0.eid, 5, 500.0, "del_flank", 0, 385)]  # [..., 99,985): just before the window
    snps = {p - 3: {4: other(bb[p - 3])}}
    indels = [(p + 2, 2, False, {4}), (p - 8, 1, True, {4})]
    got = build(ts, bb, els, events(bb, 7, [x, w, w2], losses, snps=snps, indels=indels))
    sites = [got[s][1][x.iid].site for s in range(7)]
    assert sites == [-1, -1, -1, -1, p - 4, p - 5 - 2_985, p - 5]
    target = bytearray(bb[p - 5:p])
    target[2] = ord(other(bb[p - 3]))
    assert got[4][0][p - 4:p + 1] == target
    assert all(got[s][1][x.iid].fate == "absent" for s in range(1, 7))


# ------------------------------------------------------------------------------------------------ losses

def test_solo_structure():
    ts, tree, n = toy7()
    bb, els = poptoy.toy_backbone(L, [(10_000, 300, 4_000)], seed=7)
    e = els[0]  # [10,000, 14,600), LTRs of 300
    d = DONORS[1]
    body, ltr, tsd = donor_body(d), d.ltr_left, bb[99_995:100_000]
    x = ins_at(tree, 1, None, 100_000, n["0123"], 30_000.0, strand="-", donor=1)
    losses = [loss_at(tree, x.iid, n["01"], 10_000.0, "solo", 100),
              loss_at(tree, e.eid, n["01"], 10_000.0, "solo", 120)]  # removes [10,120, 14,420)
    got = build(ts, bb, els, events(bb, 7, [x], losses))
    solo_x = revcomp_bytes(body[:100] + body[len(body) - ltr + 100:])
    for s in (0, 1):
        seq, copies, _ = got[s]
        ce, cx = copies[e.eid], copies[x.iid]
        assert seq == bb[:10_120] + bb[14_420:100_000] + solo_x + tsd + bb[100_000:]
        assert ce.fate == "solo" and (ce.start, ce.end) == (10_000, 10_300) and ce.ltr_left == ce.ltr_right == (10_000,
                                                                                                             10_300)
        assert own_kept(ce) == 300 and ce.ltr_kept == (120, 180) and ce.k2p is None and ce.events == ["solo"]
        assert seq[ce.start - T5:ce.start] == seq[ce.end:ce.end + T5] == bb[9_995:10_000]
        assert cx.fate == "solo" and cx.start == 100_000 - 4_300 and cx.end - cx.start == ltr
        assert seq[cx.start - T5:cx.start] == tsd == seq[cx.end:cx.end + T5]
        assert cx.ltr_left == cx.ltr_right == (cx.start, cx.end) and own_kept(cx) == ltr
        assert cx.ltr_kept == (ltr - 100, 100)  # genome orientation: the body's right LTR comes first on strand -
        # kept own runs in each element's own system (left LTR, own, right LTR; LTRs in genome order)
        assert ce.kept_runs == (((10_000, 10_120),), ((10_000, 10_120), (14_420, 14_600)), ((14_420, 14_600),))
        assert cx.kept_runs == (((1_850, 2_000),), ((0, 100), (1_850, 2_000)), ((0, 100),))
    for s in (2, 3):
        seq, copies, _ = got[s]
        assert seq == bb[:100_000] + revcomp_bytes(body) + tsd + bb[100_000:]
        assert copies[e.eid].fate == copies[x.iid].fate == "intact"


def test_solo_then_covered():
    # e's solo leaves an LTR; a younger del_flank of its neighbor takes it: deleted, after the solo
    ts, tree, n = toy7()
    bb, els = poptoy.toy_backbone(L, [(10_000, 300, 4_000), (15_000, 300, 1_000)], seed=22)
    e, f = els  # [10,000, 14,600) and [15,000, 16,600)
    losses = [loss_at(tree, e.eid, n["01"], 10_000.0, "solo", 120),
              loss_at(tree, f.eid, 0, 1_000.0, "del_flank", 5_000, 0)]  # [10,000, 16,600)
    got = build(ts, bb, els, events(bb, 7, [], losses))
    c0, c1 = got[0][1][e.eid], got[1][1][e.eid]
    assert got[0][0] == bb[:10_000] + bb[16_600:]
    assert c0.fate == "deleted" and c0.events == ["solo", "deleted"] and (c0.start, c0.end) == (-1, -1)
    assert c1.fate == "solo" and c1.events == ["solo"]


def test_flank_deletion():
    ts, tree, n = toy7()
    bb, els = poptoy.toy_backbone(L, [(10_000, 300, 4_000)], seed=8)
    e = els[0]
    x = ins_at(tree, 1, None, 100_000, n["0123"], 30_000.0)
    tsd, body = bb[99_995:100_000], donor_body(DONORS[0])
    losses = [loss_at(tree, e.eid, n["01"], 10_000.0, "del_flank", 300, 0),
              loss_at(tree, x.iid, n["01"], 10_000.0, "del_flank", 300, 0),
              loss_at(tree, e.eid, n["23"], 7_000.0, "del_flank", 0, 400)]
    got = build(ts, bb, els, events(bb, 7, [x], losses))
    for s in (0, 1):  # the element and 300 bp before it gone; x's start-side flank takes its target bases
        seq, copies, _ = got[s]
        assert seq == bb[:9_700] + bb[14_600:99_700] + tsd + bb[100_000:]
        for c in (copies[e.eid], copies[x.iid]):
            assert c.fate == "deleted" and (c.start, c.end) == (-1, -1) and c.events == ["del_flank"]
            assert own_kept(c) == 0
    for s in (2, 3):
        seq, copies, _ = got[s]
        assert seq == bb[:10_000] + bb[15_000:100_000] + body + tsd + bb[100_000:]
        assert copies[e.eid].fate == "deleted" and copies[x.iid].fate == "intact"
    seq, copies, _ = got[4]
    assert copies[x.iid].fate == "absent" and copies[x.iid].site == 99_995 and copies[e.eid].fate == "intact"


def test_partial_deletion():
    ts, tree, n = toy7()
    bb, els = poptoy.toy_backbone(L, [(10_000, 300, 4_000)], seed=9)
    e = els[0]
    d = DONORS[2]
    body, ltr, tsd = donor_body(d), d.ltr_left, bb[99_995:100_000]
    x = ins_at(tree, 1, None, 100_000, n["0123"], 30_000.0, donor=2)
    losses = [loss_at(tree, e.eid, n["0123"], 20_000.0, "del_partial", 500, 1_500),  # [10,800, 12,300)
              loss_at(tree, e.eid, n["01"], 5_000.0, "del_partial", 1_200, 1_600),  # [11,500, 13,100): union 2,300
              loss_at(tree, x.iid, n["23"], 6_000.0, "del_partial", 100, 800)]  # body [ltr + 100, ltr + 900)
    indels = [(12_000, 3, False, {4})]  # a small deletion inside e: the kept own bases drop, the fate does not change
    got = build(ts, bb, els, events(bb, 7, [x], losses, indels=indels))
    seq, copies, _ = got[0]
    ce = copies[e.eid]
    assert own_kept(ce) == 4_600 - 2_300 and ce.fate == "partial"
    assert ce.ltr_kept == (300, 300) and ce.events == ["del_partial", "del_partial"]
    assert seq == bb[:10_800] + bb[13_100:100_000] + body + tsd + bb[100_000:]
    seq, copies, _ = got[2]
    ce, cx = copies[e.eid], copies[x.iid]
    assert own_kept(ce) == 4_600 - 1_500 and ce.fate == "partial"
    assert cx.fate == "partial" and own_kept(cx) == len(body) - 800 and cx.k2p == 0.0
    assert seq[cx.start:cx.end] == body[:ltr + 100] + body[ltr + 900:] and cx.ltr_kept == (ltr, ltr)
    assert got[4][1][e.eid].fate == "intact" and own_kept(got[4][1][e.eid]) == 4_600 - 3
    assert got[0][1][e.eid].kept_runs == (((10_000, 10_300),), ((10_000, 10_800), (13_100, 14_600)),
                                          ((14_300, 14_600),))
    assert cx.kept_runs == (((0, 400),), ((0, 500), (1_300, 3_300)), ((2_900, 3_300),))
    assert got[4][1][e.eid].kept_runs[1] == ((10_000, 12_000), (12_003, 14_600))
    assert got[4][1][x.iid].kept_runs == ((), (), ())  # absent
    assert got[5][1][e.eid].fate == "intact" and own_kept(got[5][1][e.eid]) == 4_600


def test_replace():
    ts, tree, n = toy7()
    bb, els = poptoy.toy_backbone(L, [(10_000, 300, 4_000), (30_000, 300, 10_000), (32_000, 300, 1_000)], seed=10)
    e, e1, e2 = els  # e2 [32,000, 33,600) nests in e1's internal region [30,300, 40,300)
    x = ins_at(tree, 1, None, 100_000, n["0123"], 30_000.0, strand="-")
    losses = [loss_at(tree, e.eid, n["01"], 10_000.0, "replace", 500),
              loss_at(tree, x.iid, n["01"], 10_000.0, "replace", 700),
              loss_at(tree, e2.eid, n["0123"], 20_000.0, "replace", 300),  # older than e1's solo
              loss_at(tree, e1.eid, n["01"], 10_000.0, "solo", 100)]  # takes e2's replacing sequence on 0, 1
    got = build(ts, bb, els, events(bb, 7, [x], losses))
    seq, copies, _ = got[0]
    tsd = bb[99_995:100_000]
    assert seq[:10_000] == bb[:10_000] and set(seq[10_000:10_500]) <= set(b"ACGT")
    assert seq[10_500:30_100 - 4_100] == bb[14_600:30_100] and seq[26_000:26_000 + 59_600] == bb[40_400:100_000]
    p = 100_000 - 4_100 - 10_300
    assert seq[p + 700:] == tsd + bb[100_000:] and len(seq) == L - 4_100 - 10_300 + 700 + T5
    ce, cx = copies[e.eid], copies[x.iid]
    assert ce.fate == cx.fate == "replaced" and (ce.start, ce.end) == (10_000, 10_500)
    assert (cx.start, cx.end) == (p, p + 700) and own_kept(ce) == own_kept(cx) == 0 and ce.ltr_left == (-1, -1)
    assert ce.events == ["replace"] and ce.k2p is None
    assert copies[e2.eid].fate == "deleted" and copies[e2.eid].events == ["replace", "deleted"]
    assert copies[e1.eid].fate == "solo"
    assert got[1][0] == seq  # one replacing sequence per loss, shared by its carriers
    seq2, copies2, _ = got[2]
    assert copies2[e2.eid].fate == "replaced" and copies2[x.iid].fate == "intact"
    assert seq2[:32_000] == bb[:32_000] and seq2[32_300:32_300 + 66_400] == bb[33_600:100_000]


def test_deletion_takes_nested():
    ts, tree, n = toy7()
    bb, els = poptoy.toy_backbone(L, [(10_000, 300, 4_000)], seed=11)
    e = els[0]
    x = ins_at(tree, 1, None, 100_000, n["0123"], 35_000.0)
    y = ins_at(tree, 2, x.iid, 1_000, n["0123"], 30_000.0, depth=1, strand="-", donor=1)
    z = ins_at(tree, 3, None, 12_000, n["0123"], 30_000.0, donor=2)  # inside e's internal region
    losses = [loss_at(tree, x.iid, n["01"], 10_000.0, "del_flank", 300, 200),
              loss_at(tree, e.eid, n["01"], 10_000.0, "solo", 100)]
    got = build(ts, bb, els, events(bb, 7, [x, y, z], losses))
    for s in (0, 1):
        seq, copies, _ = got[s]
        assert seq == bb[:10_100] + bb[14_400:99_700] + bb[100_195:]
        assert [copies[i.iid].fate for i in (x, y, z)] == ["deleted"] * 3
        assert copies[x.iid].events == ["del_flank"] and copies[y.iid].events == copies[z.iid].events == ["deleted"]
    for s in (2, 3):
        assert [got[s][1][i.iid].fate for i in (x, y, z)] == ["intact"] * 3


def test_nested_in_solo_ltr_stays():
    # x's solo (crossover 150) leaves y, nested at 50 in x's left LTR, with its own younger loss
    ts, tree, n = toy7()
    bb, _ = poptoy.toy_backbone(L, seed=20)
    x = ins_at(tree, 1, None, 100_000, n["0123"], 35_000.0, donor=0)
    y = ins_at(tree, 2, x.iid, 50, n["0123"], 30_000.0, depth=1, donor=1)
    losses = [loss_at(tree, x.iid, n["01"], 10_000.0, "solo", 150), loss_at(tree, y.iid, 0, 1_000.0, "replace", 200)]
    got = build(ts, bb, [], events(bb, 7, [x, y], losses))
    bx, by = donor_body(DONORS[0]), donor_body(DONORS[1])
    seq, copies, _ = got[1]
    assert copies[x.iid].fate == "solo" and copies[y.iid].fate == "intact"
    assert seq == bb[:100_000] + bx[:50] + by + bx[45:50] + bx[50:150] + bx[len(bx) - 300 + 150:] + bb[99_995:]
    seq, copies, _ = got[0]
    assert copies[y.iid].fate == "replaced" and copies[x.iid].fate == "solo"


def test_indel_against_removal():
    ts, tree, n = toy7()
    bb, els = poptoy.toy_backbone(L, [(10_000, 300, 4_000)], seed=18)
    e = els[0]
    solo = loss_at(tree, e.eid, n["01"], 10_000.0, "solo", 100)  # removes [10,100, 14,400) on samples 0 and 1
    edge = ins_at(tree, 1, None, 10_100, n["0123"], 30_000.0)  # anchored before the solo's first removed base
    indels = [(10_098, 5, False, {0}),  # carriers within the solo's: below it, overlaps it: does not happen
              (14_398, 5, False, {0, 1, 2, 3}),  # carriers include the solo's: older, goes on outside it
              (12_000, 3, True, {0, 2}),  # inside the removed stretch: gone with it
              (14_400, 2, True, {1}),  # at its edge (14,399 removed, 14,400 kept): stays
              (10_100, 1, True, {1})]  # at its other edge (10,099 kept, 10,100 removed): stays
    got = build(ts, bb, els, events(bb, 7, [], [solo], indels=indels))
    assert got[0][0] == bb[:10_100] + bb[14_403:]
    s1 = got[1][0]
    assert s1[:10_100] == bb[:10_100] and s1[10_103:] == bb[14_403:] and set(s1[10_100:10_103]) <= set(b"ACGT")
    s2 = got[2][0]
    assert s2[:12_000] == bb[:12_000] and s2[12_003:] == bb[12_000:14_398] + bb[14_403:]
    assert got[3][0] == bb[:14_398] + bb[14_403:] and got[4][0] == bb
    # an LTR insertion anchored before the solo's first removed base stays too
    got = build(ts, bb, els, events(bb, 7, [edge], [solo]))
    pay = donor_body(DONORS[0]) + bb[10_095:10_100]
    assert got[0][0] == bb[:10_100] + pay + bb[14_400:] and got[0][1][edge.iid].fate == "intact"


def test_held_indels():
    """build_chrom's `held`: the indels at least one sample holds. A deletion the solo above it
    pre-empts and inserted bases every carrier loses with the solo are in no genome; inserted bases one carrier
    keeps, and a deletion that happens, are held. The sequences are those of the same events without `held`."""
    ts, tree, n = toy7()
    bb, els = poptoy.toy_backbone(L, [(10_000, 300, 4_000)], seed=18)
    solo = loss_at(tree, els[0].eid, n["01"], 10_000.0, "solo", 100)  # removes [10,100, 14,400) on samples 0, 1
    indels = [(10_098, 5, False, {0}),       # below the solo, overlapping it: does not happen
              (14_398, 5, False, {0, 1, 2, 3}),  # above it: happens
              (12_000, 3, True, {0, 1}),     # inside the removed stretch on both carriers: gone
              (12_500, 2, True, {1, 2}),     # gone on 1, kept on 2
              (14_400, 2, True, {1})]        # at the stretch's edge: stays
    ev = events(bb, 7, [], [solo], indels=indels)
    held: set[int] = set()
    got = build_chrom(0, bb, list(els), ts, ev, DONORS, EXACT, rng(0), held=held)
    assert held == {1, 3, 4}
    assert [s for s, _, _ in got] == [s for s, _, _ in build(ts, bb, els, ev)]


def test_clipped_event():
    # "clipped" follows a copy's own loss kinds when another element's removal took own bases of a copy still there;
    # its own removals never add it, on the backbone or in an insertion's system
    ts, tree, n = toy7()
    bb, els = poptoy.toy_backbone(L, [(10_000, 300, 4_000), (15_000, 300, 1_000)], seed=23)
    e, f = els  # [10,000, 14,600) and [15,000, 16,600)
    x = ins_at(tree, 1, None, 100_000, n["0123"], 30_000.0, donor=0)  # body 2,600
    y = ins_at(tree, 2, x.iid, 1_000, n["0123"], 25_000.0, depth=1, donor=1)
    losses = [loss_at(tree, e.eid, n["01"], 10_000.0, "del_partial", 500, 1_000),  # e's own: [10,800, 11,800)
              loss_at(tree, f.eid, 0, 1_000.0, "del_flank", 600, 0),  # [14,400, 16,600): 200 bp of e's right LTR
              loss_at(tree, f.eid, n["23"], 6_000.0, "del_flank", 600, 0),
              loss_at(tree, y.iid, 2, 1_000.0, "del_flank", 100, 0)]  # y's body and x's bases [900, 1,000)
    got = build(ts, bb, els, events(bb, 7, [x, y], losses))
    ce = [got[s][1][e.eid] for s in range(7)]
    assert [c.events for c in ce] == [["del_partial", "clipped"], ["del_partial"], ["clipped"], ["clipped"], [], [],
                                      []]
    assert [c.fate for c in ce] == ["partial"] * 4 + ["intact"] * 3
    assert (own_kept(ce[0]), ce[0].ltr_kept) == (4_600 - 1_000 - 200, (300, 100))
    assert (own_kept(ce[2]), ce[2].ltr_kept) == (4_600 - 200, (300, 100))
    cx, cy = got[2][1][x.iid], got[2][1][y.iid]
    assert cx.events == ["clipped"] and cx.fate == "partial" and own_kept(cx) == 2_600 - 100
    assert cy.events == ["del_flank"] and cy.fate == "deleted"
    assert got[3][1][x.iid].events == got[3][1][y.iid].events == [] and got[3][1][x.iid].fate == "intact"
    assert got[2][1][f.eid].events == ["del_flank"] and got[4][1][x.iid].events == []


def test_replace_anchored_at_end():
    # a backbone element's replacing sequence goes after what is anchored at the element's start and before what is
    # anchored at its end
    ts, tree, n = toy7()
    bb, els = poptoy.toy_backbone(L, [(10_000, 300, 4_000)], seed=25)
    e = els[0]  # [10,000, 14,600)
    at_start = ins_at(tree, 1, None, 10_000, n["0123"], 30_000.0, donor=0)
    at_end = ins_at(tree, 2, None, 14_600, n["0123"], 30_000.0, strand="-", donor=1)  # target: e's last 5 bases
    losses = [loss_at(tree, e.eid, n["01"], 10_000.0, "replace", 500)]
    got = build(ts, bb, els, events(bb, 7, [at_start, at_end], losses))
    pay_s = donor_body(DONORS[0]) + bb[9_995:10_000]
    pay_e = revcomp_bytes(donor_body(DONORS[1])) + bb[14_595:14_600]
    seq, copies, _ = got[0]
    a = 10_000 + len(pay_s)
    rep = seq[a:a + 500]
    assert set(rep) <= set(b"ACGT") and seq == bb[:10_000] + pay_s + rep + pay_e + bb[14_600:]
    ce = copies[e.eid]
    assert ce.fate == "replaced" and (ce.start, ce.end) == (a, a + 500)
    assert copies[at_start.iid].end + T5 == a and copies[at_end.iid].start == a + 500
    assert got[2][0] == bb[:10_000] + pay_s + bb[10_000:14_600] + pay_e + bb[14_600:]


def test_indel_bases_after_ltr_payloads():
    # at one anchor, inserted indel bases go after the LTR payloads, so the payload's TSD copies stay beside its body
    ts, tree, n = toy7()
    bb, _ = poptoy.toy_backbone(L, seed=26)
    p = 100_000
    x = ins_at(tree, 1, None, p, n["01"], 4_000.0, donor=0)
    got = build(ts, bb, [], events(bb, 7, [x], indels=[(p, 3, True, {0, 4})]))
    pay = donor_body(DONORS[0]) + bb[p - T5:p]
    seq, copies, _ = got[0]
    ib = seq[p + len(pay):p + len(pay) + 3]
    assert set(ib) <= set(b"ACGT") and seq == bb[:p] + pay + ib + bb[p:]
    c = copies[x.iid]
    assert (c.start, c.end) == (p, p + len(pay) - T5) and seq[c.start - T5:c.start] == seq[c.end:c.end + T5]
    assert got[1][0] == bb[:p] + pay + bb[p:] and got[4][0] == bb[:p] + ib + bb[p:]  # bases drawn once per event


# ------------------------------------------------------------------------------------------------ existence, K2P, gaps

def test_site_lost():
    # (a) through draw_events: insertions into sequence an older del_flank removed are discarded upstream,
    # so the builder never sees them and no sample has a copy of them
    ts, _, _ = toy7()
    bb, els = poptoy.toy_backbone(L, [(5_000 + 20_000 * i, 500, 12_000) for i in range(9)], seed=4)
    rates = Rates(lam=1e-8, delta=2e-5, loss_split=(0.0, 1.0, 0.0, 0.0), mu=0.0, indel=0.0, copy_noise=0.0)
    ev = draw_events(ts, bb, els, DONORS, rates, rng(4))
    assert ev.dropped["ins_after_loss"] >= 1
    drawn = {f"i00{k:06d}" for k in range(1, max(int(i.iid[3:]) for i in ev.ins) + 1)}
    gone = drawn - {i.iid for i in ev.ins}
    assert gone
    for _, copies, _ in build(ts, bb, els, ev, rates):
        assert set(copies) == {e.eid for e in els} | {i.iid for i in ev.ins} and not gone & set(copies)
        assert all(c.fate != "site_lost" for c in copies.values())
    # (b) hand-built: an insertion below an older del_flank that covers its target fails fast
    _, tree, n = toy7()
    flank = loss_at(tree, els[0].eid, n["0123"], 20_000.0, "del_flank", 500, 0)  # removes [4,500, 17,000)
    x = ins_at(tree, 1, None, 4_800, n["01"], 5_000.0)
    with pytest.raises(ValueError, match="existence"):
        build(ts, bb, els, events(bb, 7, [x], [flank]))


def test_backbone_k2p_adds_snps():
    ts, _, _ = toy7()
    bb, els = poptoy.toy_backbone(L, [(10_000, 300, 4_000)], seed=14)
    e = replace(els[0], k=0.02)  # identical LTRs in the backbone: the transition share is unknown, so JC
    snps = {q: {0: other(bb[q])} for q in (10_050, 10_100, 10_200)}  # 3 derived SNPs in the left LTR, sample 0
    got = build(ts, bb, [e], events(bb, 7, snps=snps))
    c0, c1 = got[0][1][e.eid], got[1][1][e.eid]
    p = 0.75 * -math.expm1(-4 / 3 * 0.02) + 3 / 300
    assert c0.k2p > 0.02 and c0.k2p == pytest.approx(-0.75 * math.log(1 - 4 / 3 * p), rel=1e-9)
    assert c1.k2p == pytest.approx(0.02, rel=1e-9)
    # a known share: two transitions and a transversion between the source's LTRs; the round trip keeps K
    bb2 = bytearray(bb)
    for q, swap in ((14_310, {"A": "G", "G": "A", "C": "T", "T": "C"}), (14_320, {"A": "G", "G": "A", "C": "T",
                                                                                 "T": "C"}),
                    (14_330, {"A": "C", "C": "A", "G": "T", "T": "G"})):
        bb2[q] = ord(swap[chr(bb2[q])])
    got2 = build(ts, bytes(bb2), [e], events(bytes(bb2), 7, snps=snps))
    assert got2[1][1][e.eid].k2p == pytest.approx(0.02, rel=1e-9)
    assert got2[0][1][e.eid].k2p > 0.02 and got2[0][1][e.eid].k2p != pytest.approx(c0.k2p, rel=1e-6)
    no_k = build(ts, bb, els, events(bb, 7, snps=snps))[0][1][e.eid]
    assert no_k.k2p is None  # no source K, no distance


def test_backbone_k2p_skips_removed_snps():
    # derived SNPs count only on own LTR bases the sample still has: not on bases a removal or a small indel took
    ts, tree, _ = toy7()
    bb, els = poptoy.toy_backbone(L, [(10_000, 300, 4_000), (15_000, 300, 1_000)], seed=24)
    e, f = replace(els[0], k=0.02), els[1]  # e [10,000, 14,600), identical LTRs: JC's share
    snps = {q: {0: other(bb[q]), 1: other(bb[q])} for q in (10_050, 10_101, 14_450, 14_500)}
    losses = [loss_at(tree, f.eid, 0, 1_000.0, "del_flank", 600, 0)]  # [14,400, 16,600): takes 14,450 and 14,500
    got = build(ts, bb, [e, f], events(bb, 7, [], losses, snps=snps, indels=[(10_100, 3, False, {0})]))
    c0, c1 = got[0][1][e.eid], got[1][1][e.eid]
    assert c0.fate == "partial" and c0.ltr_kept == (297, 100) and c1.fate == "intact"

    def jc(p):
        return -0.75 * math.log(1 - 4 / 3 * p)

    p0 = 0.75 * -math.expm1(-4 / 3 * 0.02)
    assert c0.k2p == pytest.approx(jc(p0 + 1 / 300), rel=1e-9)  # 10,050 only; the mean own LTR length is 300
    assert c1.k2p == pytest.approx(jc(p0 + 4 / 300), rel=1e-9)


def own_kept(c: Copy) -> int:
    """The own bases a copy still has."""
    return sum(hi - lo for lo, hi in c.kept_runs[1])


def mk_copy(eid: str, fate: str, start: int = -1, end: int = -1, site: int = -1) -> Copy:
    return Copy(eid, 0, fate != "absent", fate, start, end, (-1, -1), (-1, -1), (0, 0), (0, 0), site,
                (0, None, -1), None, "", None, [])


def test_gap_rate_and_near_share():
    size = 1_000_000
    seq = poptoy.rand_bases(size, random.Random(15))
    copies = {f"bb{i}": mk_copy(f"bb{i}", "intact", 50_000 + 90_000 * i, 55_000 + 90_000 * i) for i in range(10)}
    copies |= {f"i{i}": mk_copy(f"i{i}", "absent", site=95_000 + 180_000 * i) for i in range(5)}
    copies["gone"] = mk_copy("gone", "deleted")
    lift = (np.array([0, size]), np.array([0, size]))
    g = PopGenome("A1", "A", ["chr1"], [seq], copies, [], [lift])
    twin = PopGenome("A2", "A", ["chr1"], [seq], dict(copies), [], [lift])
    add_gaps([g, twin], 50, [rng(15), rng(15)])
    assert len(g.gaps) == 50 and g.gaps == twin.gaps and g.seqs == twin.seqs  # one stream per genome
    mask = np.zeros(size, dtype=bool)
    for ch, a, b in g.gaps:
        assert ch == 0 and b - a == SIM_GAP_LEN and g.seqs[0][a:b] == b"N" * SIM_GAP_LEN
        mask[a:b] = True
    out, ref = np.frombuffer(g.seqs[0], dtype=np.uint8), np.frombuffer(seq, dtype=np.uint8)
    assert np.array_equal(out[~mask], ref[~mask])
    targets = [x for c in copies.values() if c.fate == "intact" for x in (c.start, c.end)]
    targets += [c.site for c in copies.values() if c.fate == "absent"]

    def near_share(gaps) -> float:
        near = sum(any(abs(a - t) <= SIM_GAP_NEAR or abs(b - t) <= SIM_GAP_NEAR for t in targets) for _, a, b in gaps)
        return near / len(gaps)

    assert abs(near_share(g.gaps) - 1 / 3) <= 3 * math.sqrt(1 / 3 * 2 / 3 / 50)
    many = PopGenome("A3", "A", ["chr1"], [seq], dict(copies), [], [lift])  # 300 gaps: a tighter share
    assert add_gaps([many], 300, [rng(16)]) == [0]  # no near gap falls off the chromosome
    assert len(many.gaps) == 300 and abs(near_share(many.gaps) - 1 / 3) <= 3 * math.sqrt(1 / 3 * 2 / 3 / 300)


def test_gap_off_chromosome_counted():
    # the only near targets are a chromosome's two ends (an intact copy over all of it): a near gap past the start's
    # left side or the end's right side falls wholly off it, half of them; each is counted, none written
    size = 20_000
    seq = poptoy.rand_bases(size, random.Random(17))
    lift = (np.array([0, size]), np.array([0, size]))
    g = PopGenome("A1", "A", ["chr1"], [seq], {"bb1": mk_copy("bb1", "intact", 0, size)}, [], [lift])
    n = 900  # 45,000 per Mb on 20 kb; near: 300
    (off,) = add_gaps([g], n / size * 1e6, [rng(17)])
    assert off + len(g.gaps) == n and abs(off - 150) <= 3 * math.sqrt(300 / 4), off
    assert all(0 <= a < b <= size and b - a <= SIM_GAP_LEN for _, a, b in g.gaps)


def same(a, b) -> bool:
    return all(sa == sb and ca == cb and all(np.array_equal(u, v) for u, v in zip(la, lb))
               for (sa, ca, la), (sb, cb, lb) in zip(a, b)) and len(a) == len(b)


def test_deterministic_and_order_free():
    ts = ancestry(poptoy.small_demog(), Rates(), L, seed=2)
    bb, els = poptoy.toy_backbone(L, [(10_000 + 30_000 * i, 300, 6_000) for i in range(6)], seed=16)
    rates = Rates(lam=3e-9, delta=1e-5, mu=1e-7, indel=1e-8)
    ev = draw_events(ts, bb, els, DONORS, rates, rng(16))
    assert ev.ins and ev.losses and any(i and len(c) > 1 for i, c in zip(ev.indel_ins, ev.indel_carriers))
    a = build_chrom(0, bb, els, ts, ev, DONORS, rates, rng(17))
    assert same(a, build_chrom(0, bb, els, ts, ev, DONORS, rates, rng(17)))
    for _, copies, _ in a:  # the kept own runs of the LTRs hold exactly the counted bases
        for c in copies.values():
            runs = [sum(hi - lo for lo, hi in r) for r in c.kept_runs]
            assert (runs[0], runs[2]) == c.ltr_kept, c.eid
            assert all(lo < hi for r in c.kept_runs for lo, hi in r)
    chrom = _Chrom(0, bb, els, ts, ev, DONORS, rates, rng(17))
    assert same(a, [chrom.sample(s) for s in reversed(range(ts.num_samples))][::-1])
    assert [s for s, _, _ in build_chrom(0, bb, els, ts, ev, DONORS, rates, rng(18))] != [s for s, _, _ in a]


def test_bad_input_fails():
    ts, tree, _ = toy7()
    bb, els = poptoy.toy_backbone(L, [(10_000, 300, 4_000)], seed=21)
    with pytest.raises(ValueError, match="length"):
        build(ts, bb[:-1], els, events(bb[:-1]))
    with pytest.raises(ValueError, match="unknown element"):
        build(ts, bb, els, events(bb, 7, [], [loss_at(tree, "bb0099", 0, 10.0, "solo", 10)]))
    with pytest.raises(ValueError, match="host"):
        build(ts, bb, els, events(bb, 7, [ins_at(tree, 1, "i00000099", 100, 0, 10.0, depth=1)]))


def test_positions_and_times_checked():
    # an insertion's target bases and base pos lie in its host's system (the backbone, or the host's donor body), and
    # every event lies on its node's branch (node time <= time < parent time) on the tree at its position
    ts, tree, n = toy7()
    bb, els = poptoy.toy_backbone(L, [(10_000, 300, 4_000)], seed=21)
    e = els[0]
    host = ins_at(tree, 1, None, 100_000, n["0123"], 30_000.0, donor=0)  # body 2,600

    def fails(match, ins=(), losses=()):
        with pytest.raises(ValueError, match=match):
            build(ts, bb, els, events(bb, 7, ins, losses))

    for pos in (T5 - 1, L):
        fails("pos", [ins_at(tree, 2, None, pos, n["01"], 4_000.0)])
    for pos in (T5 - 1, 2_600):
        fails("pos", [host, ins_at(tree, 2, host.iid, pos, n["0123"], 20_000.0, depth=1)])
    for t in (2_999.0, 12_000.0):  # n["01"] at 3,000, its parent at 12,000
        fails("branch", [ins_at(tree, 2, None, 50_000, n["01"], t)])
        fails("branch", [host, ins_at(tree, 2, host.iid, 1_000, n["01"], t, depth=1)])
        fails("branch", losses=[loss_at(tree, e.eid, n["01"], t, "solo", 100)])
        fails("branch", [host], [loss_at(tree, host.iid, n["01"], t, "solo", 100)])
    fails("branch", [ins_at(tree, 2, None, 50_000, tree.root, 50_000.0)])  # above the root: no branch
    fails("not in the tree", [replace(ins_at(tree, 2, None, 50_000, n["01"], 4_000.0), node=99)])
    # the bounds themselves hold
    ok = [ins_at(tree, 2, None, T5, n["01"], 3_000.0), ins_at(tree, 3, None, L - 1, n["01"], 11_999.0), host,
          ins_at(tree, 4, host.iid, T5, n["0123"], 12_000.0, depth=1),
          ins_at(tree, 5, host.iid, 2_599, n["0123"], 20_000.0, depth=1)]
    got = build(ts, bb, els, events(bb, 7, ok, [loss_at(tree, e.eid, n["01"], 3_000.0, "solo", 100)]))
    seq, copies, _ = got[0]
    assert copies[ok[0].iid].start == T5 and seq[-1:] == bb[-1:] and seq[-1 - T5 - 2_600:-1] == (
        donor_body(DONORS[0]) + bb[L - 1 - T5:L - 1])
    assert [copies[i.iid].fate for i in ok] == ["intact"] * 5


def test_jitter_bounds_clamps():
    """jitter_bounds keeps both LTRs non-empty (inner ends fixed: start moves at most to the left LTR's last base,
    end to the right LTR's first) and the element on the sequence; without jitter nothing moves."""
    jb = jitter_bounds
    assert jb(100, 200, (100, 102), (199, 200), 1_000, True, 3, -3) == (101, 200, 101, 200)  # 2-bp and 1-bp LTRs
    assert jb(100, 200, (100, 103), (197, 200), 1_000, True, 3, -3) == (102, 198, 102, 198)  # 3-bp LTRs: 2 each
    assert jb(100, 200, (100, 101), (199, 200), 1_000, True, 3, -3) == (100, 200, 100, 200)  # 1-bp LTRs
    assert jb(100, 200, (100, 102), (199, 200), 1_000, True, -3, 3) == (97, 203, 97, 203)  # outward: free
    assert jb(0, 50, (0, 10), (40, 50), 50, True, -3, 3) == (0, 50, 0, 50)  # at both sequence ends
    assert jb(1, 49, (1, 10), (40, 49), 50, True, -3, 3) == (0, 50, 0, 50)
    assert jb(100, 200, (100, 102), (199, 200), 1_000, False, 3, -3) == (100, 200, 100, 200)
