"""Population truth, annotations, lifted tracts and the insertion table (synltr.sim.pop.truth), on genomes build
builds from exact toy genealogies (poptoy.toy_ts) with hand-placed events: one test per truth rule, exact states."""
from __future__ import annotations

import pytest

pytest.importorskip("msprime")

import csv
import math
import random
from dataclasses import replace
from pathlib import Path

import numpy as np
from sim import poptoy
from sim.poptoy import DONORS, L, build, events, ins_at, loss_at, other, sim_main, toy7

from synltr.base.constants import COV_PRESENT, MIN_SEG, REFINE_PAD, TSD_DEFAULT
from synltr.calls.records import CALL_COLUMNS
from synltr.io.annotation import load_annotation
from synltr.io.annotation.gff3 import write_gff3
from synltr.io.elements import make_eid, write_elements
from synltr.io.seqio import Fasta
from synltr.sim.pop.build import Copy, PopGenome
from synltr.sim.pop.events import donor_body
from synltr.sim.pop.model import Rates
from synltr.sim.pop.truth import INSERTION_COLUMNS, alive, annotate, insertion_rows, pair_truth, tract_bed
from synltr.sim.score import TRUTH_COLUMNS

T5 = TSD_DEFAULT


def pop(got, seqid: str = "chr1") -> list[PopGenome]:
    """PopGenomes S0..S6 (group A) of one built chromosome."""
    return [PopGenome(f"S{s}", "A", [seqid], [seq], copies, [], [lift]) for s, (seq, copies, lift) in enumerate(got)]


def states(x: PopGenome, y: PopGenome) -> dict[str, str]:
    """eid -> state (state:flags when flagged) of every alive element of x, read in y."""
    return {r["eid"]: r["state"] + (f":{r['flags']}" if r["flags"] else "") for r in pair_truth(x, y)}


def n_at(g: PopGenome, q: int) -> PopGenome:
    """g with an N at position q of its first chromosome."""
    s = g.seqs[0]
    return replace(g, seqs=[s[:q] + b"N" + s[q + 1:]])


def pay(x) -> int:
    """Length of x's payload without anything nested: donor body and second TSD copy."""
    return len(donor_body(DONORS[x.donor])) + T5


def own_kept(c: Copy) -> int:
    """The own bases a copy still has."""
    return sum(hi - lo for lo, hi in c.kept_runs[1])


def mk_copy(eid: str = "e", fate: str = "intact", ltr_kept=(300, 300), start: int = 0) -> Copy:
    """A 1,000-bp backbone copy at start with 300-bp LTRs (own bases all kept unless ltr_kept says otherwise)."""
    return Copy(eid, 0, fate != "absent", fate, start, start + 1_000, (start, start + 300),
                (start + 700, start + 1_000), (300, 300), ltr_kept, -1, (0, None, -1), None, "",
                0.01, [])


# ------------------------------------------------------------------------------------------------ alive

def test_alive():
    # alive: intact or partial, and each LTR keeps >= COV_PRESENT of its own bases
    assert COV_PRESENT == 0.8
    assert alive(mk_copy()) and alive(mk_copy(fate="partial", ltr_kept=(240, 300)))
    assert not alive(mk_copy(fate="partial", ltr_kept=(239, 300))) and not alive(mk_copy(fate="partial",
                                                                                         ltr_kept=(300, 239)))
    assert not any(alive(mk_copy(fate=f)) for f in ("solo", "deleted", "replaced", "absent"))


# ------------------------------------------------------------------------------------------------ truth rules

def test_truth_P_inherited_intact():
    ts, tree, n = toy7()
    bb, els = poptoy.toy_backbone(L, [(10_000, 300, 4_000)], seed=31)
    e = replace(els[0], k=0.02)
    x = ins_at(tree, 1, None, 100_000, n["0123"], 30_000.0, donor=0)  # body 2,600
    gs = pop(build(ts, bb, [e], events(bb, 7, [x])))
    rows = pair_truth(gs[0], gs[1])
    common = {"genome": "S0", "seqid": "chr1", "state": "P", "flags": "", "events_x": "", "events_y": ""}
    assert rows == [dict(common, eid=e.eid, name=e.eid, start=10_001, end=14_600, origin="backbone"),
                    dict(common, eid=x.iid, name=x.iid, start=100_001, end=102_600, origin="sim")]
    assert all(tuple(r) == TRUTH_COLUMNS for r in rows)
    assert all(states(gs[0], gs[y]) == {e.eid: "P", x.iid: "P"} for y in (1, 2, 3))


def test_truth_P_shared_loss():
    # Y is read against X's own bases ("Y holds >= COV_PRESENT of e's own ids"): a loss X and Y share leaves
    # e P in Y however much it took (the kept share of its own bases is 0.57 and 0.55 here), while X = an
    # intact copy reads it T
    ts, tree, n = toy7()
    bb, els = poptoy.toy_backbone(L, [(10_000, 300, 4_000)], seed=32)
    e = els[0]
    x = ins_at(tree, 1, None, 100_000, n["0123"], 30_000.0, donor=2)  # body 3,300, LTRs 400
    losses = [loss_at(tree, e.eid, n["0123"], 20_000.0, "del_partial", 500, 2_000),  # [10,800, 12,800)
              loss_at(tree, x.iid, n["0123"], 20_000.0, "del_partial", 100, 1_500)]
    gs = pop(build(ts, bb, [e], events(bb, 7, [x], losses)))
    c0 = gs[0].copies
    assert (own_kept(c0[e.eid]), own_kept(c0[x.iid])) == (2_600, 1_800) and alive(c0[e.eid]) and alive(c0[x.iid])
    assert states(gs[0], gs[1]) == {e.eid: "P", x.iid: "P"}
    assert states(gs[0], gs[4]) == {e.eid: "P", x.iid: "E"}  # S4 holds all of S0's own bases, and more
    assert states(gs[4], gs[0]) == {e.eid: "T"}  # S0 holds 2,600 of S4's 4,600
    assert pair_truth(gs[0], gs[1])[0]["events_x"] == "del_partial"


def test_truth_E_never_inserted():
    ts, tree, n = toy7()
    bb, els = poptoy.toy_backbone(L, [(10_000, 300, 4_000)], seed=35)
    e = els[0]
    x = ins_at(tree, 1, None, 100_000, n["01"], 4_000.0, donor=0)
    y = ins_at(tree, 2, None, 12_000, n["01"], 4_000.0, donor=1)  # inside e: its empty site lies in e
    gs = pop(build(ts, bb, [e], events(bb, 7, [x, y])))
    rows = {r["eid"]: r for r in pair_truth(gs[0], gs[4])}
    assert {k: r["state"] for k, r in rows.items()} == {e.eid: "P", x.iid: "E", y.iid: "E"}
    assert (rows[x.iid]["origin"], rows[x.iid]["flags"], rows[x.iid]["events_y"]) == ("sim", "", "")
    assert gs[4].copies[x.iid].site == 99_995 and gs[4].copies[y.iid].site == 11_995


def test_truth_S_solo_in_y():
    ts, tree, n = toy7()
    bb, els = poptoy.toy_backbone(L, [(10_000, 300, 4_000), (15_000, 300, 1_000)], seed=34)
    e, f = els  # [10,000, 14,600) and [15,000, 16,600)
    x = ins_at(tree, 1, None, 100_000, n["0123"], 30_000.0, strand="-", donor=1)  # LTRs 250
    losses = [loss_at(tree, e.eid, n["23"], 6_000.0, "solo", 100),  # keeps [10,000, 10,100) + [14,400, 14,600)
              loss_at(tree, x.iid, n["23"], 6_000.0, "solo", 100),
              loss_at(tree, f.eid, 3, 1_000.0, "del_flank", 750, 0)]  # [14,250, 16,600): S3's solo keeps 100 bp
    gs = pop(build(ts, bb, els, events(bb, 7, [x], losses)))
    assert states(gs[0], gs[2]) == {e.eid: "S", f.eid: "P", x.iid: "S"}
    assert pair_truth(gs[0], gs[2])[0]["events_y"] == "solo"
    # S needs >= COV_PRESENT of an LTR's worth of X's own bases: S3's clipped solo holds 100 of 300, T
    assert gs[3].copies[e.eid].fate == "solo" and states(gs[0], gs[3]) == {e.eid: "T", f.eid: "D", x.iid: "S"}


def test_truth_T_partial():
    ts, tree, _ = toy7()
    bb, els = poptoy.toy_backbone(L, [(10_000, 300, 4_000), (15_000, 300, 1_000), (20_000, 100, 200),
                                      (20_500, 100, 300)], seed=33)
    e, h, g, f = els  # [10,000, 14,600), [15,000, 16,600), [20,000, 20,400) (own 400) and [20,500, 21,000)
    losses = [loss_at(tree, e.eid, 5, 1_000.0, "del_partial", 0, 1_500),  # S5 holds 3,100 of 4,600
              loss_at(tree, e.eid, 0, 1_000.0, "del_partial", 0, 1_000),  # S0: [10,300, 11,300)
              loss_at(tree, e.eid, 2, 1_000.0, "del_partial", 2_000, 1_000),  # S2: [12,300, 13,300)
              loss_at(tree, h.eid, 1, 1_000.0, "del_flank", 461, 0),  # [14,539, 16,600): e's right LTR keeps 239
              loss_at(tree, h.eid, 3, 1_000.0, "del_flank", 460, 0),  # keeps 240 of 300
              loss_at(tree, f.eid, 3, 1_000.0, "del_flank", 420, 0),  # [20,080, 21,000): g keeps 80 bp
              loss_at(tree, f.eid, 6, 1_000.0, "del_flank", 421, 0)]  # g keeps 79
    gs = pop(build(ts, bb, els, events(bb, 7, [], losses)))
    assert states(gs[4], gs[5])[e.eid] == "T"
    # P needs each LTR too: S1 keeps 4,539 of 4,600 own bases but 239 of the right LTR's 300
    assert states(gs[4], gs[1])[e.eid] == "T" and states(gs[4], gs[3])[e.eid] == "P"
    # each lost its own 1,000 bp: S2 holds 2,600 of S0's 3,600 own bases (0.72), though both keep 3,600
    assert own_kept(gs[0].copies[e.eid]) == own_kept(gs[2].copies[e.eid]) == 3_600
    assert states(gs[0], gs[2])[e.eid] == states(gs[2], gs[0])[e.eid] == "T"
    assert MIN_SEG == 80 and states(gs[4], gs[3])[g.eid] == "T" and states(gs[4], gs[6])[g.eid] == "D"


def test_truth_D_deleted_in_y():
    ts, tree, n = toy7()
    bb, els = poptoy.toy_backbone(L, [(10_000, 300, 4_000)], seed=36)
    e = els[0]
    x = ins_at(tree, 1, None, 100_000, n["0123"], 30_000.0, donor=0)
    y = ins_at(tree, 2, x.iid, 1_000, n["0123"], 20_000.0, depth=1, donor=1)
    losses = [loss_at(tree, e.eid, 4, 1_000.0, "del_flank", 100, 100),
              loss_at(tree, x.iid, 2, 1_000.0, "del_flank", 100, 100)]  # y goes with x's body
    gs = pop(build(ts, bb, [e], events(bb, 7, [x, y], losses)))
    rows = {r["eid"]: r for r in pair_truth(gs[0], gs[2])}
    assert {k: r["state"] for k, r in rows.items()} == {e.eid: "P", x.iid: "D", y.iid: "D"}
    assert (rows[x.iid]["events_y"], rows[y.iid]["events_y"]) == ("del_flank", "deleted")
    assert states(gs[0], gs[4]) == {e.eid: "D", x.iid: "E", y.iid: "E:host_empty"}


def test_truth_D_site_deleted():
    # no site_lost state: e not in Y, and Y's empty-site window [p - 15, p + 10) broken by a removal: D
    ts, tree, n = toy7()
    bb, els = poptoy.toy_backbone(L, [(97_000, 300, 2_000)], seed=37)
    e0 = els[0]  # [97,000, 99,600)
    p = 100_000
    x = ins_at(tree, 1, None, p, n["01"], 4_000.0, donor=0)
    h = ins_at(tree, 2, None, 150_000, n["0123"], 30_000.0, donor=0)
    z = ins_at(tree, 3, h.iid, 1_000, 0, 1_000.0, depth=1, donor=1)  # S0 only
    losses = [loss_at(tree, e0.eid, 4, 1_000.0, "del_flank", 0, 386),  # [97,000, 99,986): the window's first base
              loss_at(tree, e0.eid, 5, 1_000.0, "del_flank", 0, 385),  # [97,000, 99,985): just before it
              loss_at(tree, h.iid, 2, 1_000.0, "del_flank", 100, 100)]  # S2: z's site goes with h
    gs = pop(build(ts, bb, els, events(bb, 7, [x, h, z], losses)))
    assert gs[4].copies[x.iid].site == -1 and states(gs[0], gs[4])[x.iid] == "D"
    assert states(gs[0], gs[5])[x.iid] == "E"
    assert gs[2].copies[z.iid].site == -1 and states(gs[0], gs[2])[z.iid] == "D"
    assert states(gs[0], gs[3])[z.iid] == "E"


def test_truth_O_replaced():
    ts, tree, n = toy7()
    bb, els = poptoy.toy_backbone(L, [(10_000, 300, 4_000)], seed=38)
    e = els[0]
    x = ins_at(tree, 1, None, 100_000, n["0123"], 30_000.0, donor=0)
    losses = [loss_at(tree, e.eid, 4, 1_000.0, "replace", 500), loss_at(tree, x.iid, 2, 1_000.0, "replace", 700)]
    gs = pop(build(ts, bb, [e], events(bb, 7, [x], losses)))
    assert states(gs[0], gs[4]) == {e.eid: "O", x.iid: "E"} and states(gs[0], gs[2]) == {e.eid: "P", x.iid: "O"}
    assert gs[2].copies[x.iid].fate == "replaced" and own_kept(gs[2].copies[x.iid]) == 0


def test_truth_O_same_target():
    # e not in Y, and Y holds another insertion whose target bases ([pos - 5, pos) of the same host's system) overlap
    # e's: O. Not when the targets only touch, when the host differs, or when Y's copy of the other is gone
    ts, tree, n = toy7()
    bb, _ = poptoy.toy_backbone(L, seed=39)
    p = 100_000
    x = ins_at(tree, 1, None, p, n["01"], 4_000.0, donor=0)
    x2 = ins_at(tree, 2, None, 1_000, n["01"], 4_000.0, donor=0)
    near = [ins_at(tree, 3, None, p + 4, 2, 1_000.0, donor=1),  # S2: the targets share base p - 1
            ins_at(tree, 4, None, p + 5, 3, 1_000.0, donor=1),  # S3: they touch
            ins_at(tree, 5, None, p - 4, 4, 1_000.0, donor=2),  # S4
            ins_at(tree, 6, None, p + 2, 5, 1_000.0, donor=1)]  # S5: its body deleted below
    h = ins_at(tree, 7, None, 150_000, 6, 2_000.0, donor=0)
    y = ins_at(tree, 8, h.iid, 1_000, 6, 1_000.0, depth=1, donor=1)  # S6: x2's pos, in h's system
    losses = [loss_at(tree, near[3].iid, 5, 500.0, "del_flank", 0, 0)]  # its body only: S5's window stays intact
    gs = pop(build(ts, bb, [], events(bb, 7, [x, x2, *near, h, y], losses)))
    got = [states(gs[0], gs[s]) for s in range(2, 7)]
    assert [g[x.iid] for g in got] == ["O", "E", "O", "E", "E"]
    assert gs[5].copies[near[3].iid].fate == "deleted" and got[4][x2.iid] == "E"


def test_truth_host_empty():
    # e nested in a simulated host Y never got: the outermost such host's truth, plus host_empty when that is E
    ts, tree, n = toy7()
    bb, els = poptoy.toy_backbone(L, [(97_000, 300, 2_000)], seed=40)
    e0 = els[0]  # [97,000, 99,600)
    h = ins_at(tree, 1, None, 100_000, n["01"], 4_000.0, donor=2)
    e = ins_at(tree, 2, h.iid, 1_000, 0, 2_000.0, depth=1, donor=0)  # S0 only
    e2 = ins_at(tree, 3, e.iid, 500, 0, 1_000.0, depth=2, donor=1)
    losses = [loss_at(tree, e0.eid, 5, 1_000.0, "del_flank", 0, 386)]  # breaks h's window in S5
    gs = pop(build(ts, bb, els, events(bb, 7, [h, e, e2], losses)))
    assert states(gs[0], gs[1]) == {e0.eid: "P", h.iid: "P", e.iid: "E", e2.iid: "E:host_empty"}
    assert states(gs[0], gs[4]) == {e0.eid: "P", h.iid: "E", e.iid: "E:host_empty", e2.iid: "E:host_empty"}
    assert states(gs[0], gs[5]) == {e0.eid: "D", h.iid: "D", e.iid: "D", e2.iid: "D"}
    gapped = n_at(gs[6], gs[6].copies[h.iid].site)
    assert states(gs[0], gapped) == {e0.eid: "P", h.iid: "G", e.iid: "G", e2.iid: "G"}


def test_truth_G_gap_near_site():
    # an N in Y within REFINE_PAD of e's site (absent: its target bases) or of its image there
    ts, tree, n = toy7()
    bb, els = poptoy.toy_backbone(L, [(10_000, 300, 4_000)], seed=41)
    e = els[0]
    x = ins_at(tree, 1, None, 100_000, n["01"], 4_000.0, donor=0)
    gs = pop(build(ts, bb, [e], events(bb, 7, [x])))
    s = gs[4].copies[x.iid].site
    assert s == 99_995
    for q, want in ((s - REFINE_PAD, "G"), (s - REFINE_PAD - 1, "E"), (s + T5 - 1 + REFINE_PAD, "G"),
                    (s + T5 + REFINE_PAD, "E")):
        assert states(gs[0], n_at(gs[4], q))[x.iid] == want, q
    c = gs[1].copies[e.eid]
    assert (c.start, c.end) == (10_000, 14_600)
    for q, want in ((c.start - REFINE_PAD, "G"), (c.start - REFINE_PAD - 1, "P"), (12_000, "G"),
                    (c.end - 1 + REFINE_PAD, "G"), (c.end + REFINE_PAD, "P")):
        assert states(gs[0], n_at(gs[1], q))[e.eid] == want, q


def test_truth_E_with_snp_and_indel_beside_site():
    # a SNP in the TSD and 2-bp indels 5 bp from the site leave it E (point mutations are noise)
    ts, tree, n = toy7()
    bb, _ = poptoy.toy_backbone(L, seed=42)
    p = 100_000
    x = ins_at(tree, 1, None, p, n["01"], 4_000.0, donor=0)
    snps = {p - 3: {4: other(bb[p - 3])}}  # in x's target bases, S4 only
    indels = [(p + 5, 2, False, {4}), (p - 10, 2, True, {5})]  # 5 bp past the target, and 5 before it
    gs = pop(build(ts, bb, [], events(bb, 7, [x], snps=snps, indels=indels)))
    assert gs[4].seqs[0][p - 3] != bb[p - 3] and len(gs[4].seqs[0]) == L - 2 and len(gs[5].seqs[0]) == L + 2
    assert states(gs[0], gs[4]) == {x.iid: "E"} and states(gs[0], gs[5]) == {x.iid: "E"}


def test_truth_E_with_other_insertion_in_window():
    # Y-only insertions whose targets lie 4 bp from e's, either side, inside e's window [p - 15, p + 10): still E
    ts, tree, n = toy7()
    bb, _ = poptoy.toy_backbone(L, seed=43)
    p = 100_000
    x = ins_at(tree, 1, None, p, n["01"], 4_000.0, donor=0)
    w = ins_at(tree, 2, None, p + 9, 4, 1_000.0, donor=1)  # target [p + 4, p + 9)
    w2 = ins_at(tree, 3, None, p - 9, 5, 1_000.0, donor=1)  # target [p - 14, p - 9)
    gs = pop(build(ts, bb, [], events(bb, 7, [x, w, w2])))
    assert gs[4].copies[w.iid].fate == gs[5].copies[w2.iid].fate == "intact"
    assert states(gs[0], gs[4]) == {x.iid: "E"} and states(gs[0], gs[5]) == {x.iid: "E"}


def test_truth_input_checks():
    ts, tree, n = toy7()
    bb, _ = poptoy.toy_backbone(L, seed=53)
    gs = pop(build(ts, bb, [], events(bb, 7, [ins_at(tree, 1, None, 100_000, n["01"], 4_000.0)])))
    with pytest.raises(ValueError, match="label"):
        pair_truth(gs[0], replace(gs[1], label="S0"))
    with pytest.raises(ValueError, match="chromosomes"):
        pair_truth(gs[0], replace(gs[1], seqids=["chr2"]))
    with pytest.raises(ValueError, match="copy"):
        pair_truth(gs[0], replace(gs[1], copies={}))


# ------------------------------------------------------------------------------------------------ annotation

def test_annotation_dropout_bounds():
    n = 2_000
    seq = poptoy.rand_bases(1_000 * n + 1_000, random.Random(44))
    copies = {f"e{i}": mk_copy(f"e{i}", start=500 + 1_000 * i) for i in range(n)}
    dead = {f"d{i}": mk_copy(f"d{i}", fate=f, start=500 + 1_000 * i)
            for i, f in enumerate(("solo", "deleted", "replaced", "absent"))}
    dead["d4"] = mk_copy("d4", fate="partial", ltr_kept=(239, 300), start=4_500)
    lift = (np.array([0, len(seq)]), np.array([0, len(seq)]))
    g = PopGenome("A1", "A", ["c01"], [seq], copies | dead, [], [lift])
    els = annotate(g, 0.1, 0.0, np.random.default_rng(45))
    assert abs(len(els) - 0.9 * n) <= 3 * math.sqrt(n * 0.1 * 0.9), len(els)
    assert {x.name for x in els} <= set(copies)
    assert all((x.start, x.end) == (copies[x.name].start + 1, copies[x.name].end) for x in els)
    jit = annotate(g, 0.0, 1.0, np.random.default_rng(46))  # every one jittered: ends move by U{-3..3}
    assert len(jit) == n
    moved = 0
    for x in jit:
        c = copies[x.name]
        ds, dt = x.start - 1 - c.start, x.end - c.end
        assert abs(ds) <= 3 and abs(dt) <= 3
        assert x.ltr_left == (c.ltr_left[0] + 1 + ds, c.ltr_left[1]) and x.ltr_right == (c.ltr_right[0] + 1,
                                                                                         c.ltr_right[1] + dt)
        moved += ds != 0 or dt != 0
    assert moved > n * 0.9  # 1 - (1/7)^2 of them


def write_fasta(path: Path, g: PopGenome) -> None:
    with open(path, "wb") as fh:
        fh.writelines(b">" + sid.encode() + b"\n" + b"".join(s[i:i + 80] + b"\n" for i in range(0, len(s), 80))
                      for sid, s in zip(g.seqids, g.seqs))


def test_annotation_loads(tmp_path):
    ts, tree, n = toy7()
    bb, els = poptoy.toy_backbone(L, [(10_000, 300, 4_000), (30_000, 300, 6_000), (32_000, 300, 1_000),
                                      (60_000, 300, 2_000)], seed=47)
    els = [replace(e, k=0.01 * (i + 1), family=f"bb{i}", superfamily="Copia", strand="+-"[i % 2])
           for i, e in enumerate(els)]  # els[2] nests in els[1]
    xs = [ins_at(tree, 1, None, 100_000, n["0123"], 30_000.0, donor=0),
          ins_at(tree, 2, None, 120_000, n["0123"], 30_000.0, strand="-", donor=1),
          ins_at(tree, 3, None, 12_000, n["01"], 4_000.0, donor=2)]  # inside els[0]
    nested = ins_at(tree, 4, xs[0].iid, 1_000, n["01"], 5_000.0, depth=1, strand="-", donor=2)
    losses = [loss_at(tree, els[3].eid, n["01"], 10_000.0, "solo", 100),  # not alive: no record
              loss_at(tree, xs[1].iid, 0, 1_000.0, "del_partial", 100, 500)]  # partial, alive
    g = pop(build(ts, bb, els, events(bb, 7, [*xs, nested], losses), rates=Rates(mu=0.0, copy_noise=0.0)))[0]
    live = sorted((c for c in g.copies.values() if alive(c)), key=lambda c: (c.chrom, c.start, -c.end, c.eid))
    ann = annotate(g, 0.0, 0.0, np.random.default_rng(48))
    assert len(live) == 7 and [x.name for x in ann] == [c.eid for c in live]
    write_fasta(tmp_path / "S0.fa", g)
    write_gff3(ann, tmp_path / "S0.gff3")
    with Fasta(tmp_path / "S0.fa") as fa:
        loaded = load_annotation(str(tmp_path / "S0.gff3"), "S0", fa)
    by = {x.name: x for x in loaded}
    assert len(loaded) == 7 and set(by) == {c.eid for c in live}
    for c in live:
        x = by[c.eid]  # Name = the simulator's eid
        assert x.eid == make_eid("S0", "chr1", c.start + 1, c.end) and (x.seqid, x.strand) == ("chr1", c.strand)
        assert x.ltr_left == (c.ltr_left[0] + 1, c.ltr_left[1]) and x.ltr_right == (c.ltr_right[0] + 1, c.ltr_right[1])
        assert x.tsd == c.tsd and len(x.tsd) == T5 and c.k2p is not None and x.k == c.k2p  # K2P_d
        assert (x.family, x.superfamily) == (c.family, c.superfamily)
    assert by[nested.iid].parent == by[xs[0].iid].eid and by[els[2].eid].parent == by[els[1].eid].eid


# ------------------------------------------------------------------------------------------------ tracts

def test_tracts_lifted():
    ts, tree, n = toy7()
    bb, els = poptoy.toy_backbone(L, [(20_000, 300, 4_000), (45_000, 300, 2_000), (120_000, 300, 3_000)], seed=49)
    _, b, c = els  # [20,000, 24,600), [45,000, 47,600), [120,000, 123,600)
    xs = [ins_at(tree, 1, None, 30_000, n["0123"], 30_000.0, donor=0),  # inside tract 1
          ins_at(tree, 2, None, 50_000, n["0123"], 30_000.0, donor=1),  # anchored at tract 1's end: in it
          ins_at(tree, 3, None, 10_000, n["0123"], 30_000.0, donor=2),  # anchored at tract 1's start: not
          ins_at(tree, 4, None, 80_000, n["0123"], 30_000.0, donor=0)]  # between tracts
    nested = ins_at(tree, 5, xs[0].iid, 1_000, n["01"], 5_000.0, depth=1, donor=1)
    losses = [loss_at(tree, b.eid, n["01"], 10_000.0, "del_flank", 0, 1_400),  # [45,000, 49,000)
              loss_at(tree, c.eid, n["01"], 10_000.0, "del_flank", 0, 0)]  # [120,000, 123,600): tract 3 goes
    indels = [(15_000, 4, False, {0}), (35_000, 3, True, {0})]
    tracts = [(10_000, 50_000), (100_000, 110_000), (121_000, 122_000)]
    gs = pop(build(ts, bb, els, events(bb, 7, [*xs, nested], losses, indels=indels)), seqid="c07")
    p0, p1, p2, p3, pn = (pay(x) for x in (*xs, nested))
    before = p2 - 4 + (p0 + pn) + 3 - 4_000 + p1  # S0's shift of backbone positions 50,000-79,999
    assert tract_bed(gs[0], {0: tracts}) == [("c07", 10_000 + p2, 50_000 + before),
                                            ("c07", 100_000 + before + p3, 110_000 + before + p3)]
    q = p2 + p0 + p1 + p3  # S2: the four primaries; no nested insertion, loss or indel
    assert tract_bed(gs[2], {0: tracts}) == [("c07", 10_000 + p2, 50_000 + p2 + p0 + p1),
                                            ("c07", 100_000 + q, 110_000 + q), ("c07", 121_000 + q, 122_000 + q)]
    assert tract_bed(gs[4], {0: tracts}) == [("c07", *t) for t in tracts] and tract_bed(gs[4], {}) == []
    # every copy whose backbone point lies in a tract lies in its lifted tract; one wholly outside every tract
    # lies outside every lifted one (an insertion's point: its anchor, after base pos - 1; nested: its primary's)
    span = {e.eid: (e.start, e.end) for e in els}
    anchor = {x.iid: x.pos for x in xs} | {nested.iid: xs[0].pos}
    checked = 0
    for g in gs:
        bed = tract_bed(g, {0: tracts})
        for cp in g.copies.values():
            if cp.start < 0 or cp.fate == "replaced":
                continue
            lo, hi = span[cp.eid] if cp.eid in span else (anchor[cp.eid] - 1, anchor[cp.eid])
            if any(u <= lo and hi <= v for u, v in tracts):
                assert any(s <= cp.start and cp.end <= t for _, s, t in bed), (g.label, cp.eid)
                checked += 1
            elif not any(u < hi and lo < v for u, v in tracts):
                assert not any(s < cp.end and cp.start < t for _, s, t in bed), (g.label, cp.eid)
                checked += 1
    assert checked == 6 * 4 + 2 * 3  # S0-S3: six copies each (one overlaps a tract's end); S4-S6: a and b
    with pytest.raises(ValueError, match="chromosome"):
        tract_bed(gs[0], {1: tracts})


def test_tracts_deleted_end():
    # the named limit: a tract end inside a stretch the sample deleted lifts to the next kept base, so a payload
    # anchored at that base falls on the wrong side: outside the record at a left end, inside it at a right end
    ts, tree, n = toy7()
    bb, els = poptoy.toy_backbone(L, [(20_000, 300, 1_000), (60_000, 300, 1_000)], seed=54)
    e1, e2 = els  # [20,000, 21,600) and [60,000, 61,600)
    left = ins_at(tree, 1, None, 22_000, n["0123"], 30_000.0, donor=0)  # first kept base after [20,000, 22,000)
    right = ins_at(tree, 2, None, 61_600, n["0123"], 30_000.0, donor=1)  # first kept base after [59,600, 61,600)
    losses = [loss_at(tree, e1.eid, n["01"], 10_000.0, "del_flank", 0, 400),  # [20,000, 22,000)
              loss_at(tree, e2.eid, n["01"], 10_000.0, "del_flank", 400, 0)]  # [59,600, 61,600)
    gs = pop(build(ts, bb, els, events(bb, 7, [left, right], losses)))
    tract = {0: [(21_000, 60_500)]}  # both ends deleted in S0, kept in S2
    pl, pr = pay(left), pay(right)
    c0 = gs[0].copies
    assert (c0[left.iid].start, c0[right.iid].start) == (20_000, 57_600 + pl)
    assert tract_bed(gs[0], tract) == [("chr1", 20_000 + pl, 57_600 + pl + pr)]  # left out; right in
    c2 = gs[2].copies
    assert (c2[left.iid].start, c2[right.iid].start) == (22_000, 61_600 + pl)
    assert tract_bed(gs[2], tract) == [("chr1", 21_000, 60_500 + pl)]  # kept ends: the anchor convention


# ------------------------------------------------------------------------------------------------ insertion table

def test_insertion_rows():
    ts, tree, n = toy7()
    bb, els = poptoy.toy_backbone(L, [(10_000, 300, 4_000)], seed=52)
    x = ins_at(tree, 1, None, 100_000, n["0123"], 30_000.0, donor=0)
    y = ins_at(tree, 2, x.iid, 1_000, n["01"], 5_000.0, depth=1, strand="-", donor=1)
    losses = [loss_at(tree, x.iid, 2, 1_000.0, "solo", 100), loss_at(tree, x.iid, 3, 1_000.0, "del_flank", 50, 50),
              loss_at(tree, y.iid, 1, 1_000.0, "del_partial", 100, 1_000)]  # S1's y: partial, LTRs whole: alive
    gs = pop(build(ts, bb, els, events(bb, 7, [x, y], losses)))
    rows = insertion_rows(gs, {0: [x, y]})
    fates = tuple(f"fate_S{s}" for s in range(7))
    assert [tuple(r) for r in rows] == [INSERTION_COLUMNS + fates] * 2
    assert rows[0] == dict(zip(INSERTION_COLUMNS + fates, (
        x.iid, "chr1", ".", 100_000, 30_000.0, n["0123"], 0, "fam0", "S0,S1,S2,S3", "S0,S1",
        "intact", "intact", "solo", "deleted", "absent", "absent", "absent")))
    assert rows[1] == dict(zip(INSERTION_COLUMNS + fates, (
        y.iid, "chr1", x.iid, 1_000, 5_000.0, n["01"], 1, "fam1", "S0,S1", "S0,S1",
        "intact", "partial", "absent", "absent", "absent", "absent", "absent")))
    with pytest.raises(ValueError, match="chromosome"):
        insertion_rows(gs, {1: [x]})
    with pytest.raises(ValueError, match="carrier"):
        insertion_rows(gs[:3], {0: [x]})


# ------------------------------------------------------------------------------------------------ score join

def write_tsv(path: Path, columns, rows) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(columns)
        w.writerows([r.get(c, "") for c in columns] for r in rows)


def test_score_joins(tmp_path):
    # truth rows of X in Y, X's annotation and calls equal to the truth: `synltr sim score` joins them (call eid ->
    # element Name -> truth name), every state right
    ts, tree, n = toy7()
    bb, els = poptoy.toy_backbone(L, [(10_000, 300, 4_000), (30_000, 300, 4_000)], seed=50)
    xs = [ins_at(tree, 1, None, 100_000, n["01"], 4_000.0, donor=0),
          ins_at(tree, 2, None, 120_000, n["0123"], 30_000.0, donor=1)]
    losses = [loss_at(tree, els[1].eid, n["45"], 10_000.0, "solo", 100)]
    gs = pop(build(ts, bb, els, events(bb, 7, xs, losses)))
    truth = pair_truth(gs[0], gs[4]) + pair_truth(gs[4], gs[0])  # both directions in one file
    assert sorted(r["state"] for r in truth if r["genome"] == "S0") == ["E", "E", "P", "S"]
    write_tsv(tmp_path / "truth.tsv", TRUTH_COLUMNS, truth)
    ann = annotate(gs[0], 0.0, 0.0, np.random.default_rng(51))
    write_elements(ann, tmp_path / "S0.tsv")
    state = {r["name"]: r["state"] for r in truth if r["genome"] == "S0"}
    write_tsv(tmp_path / "calls.tsv", CALL_COLUMNS,
              [{"eid": x.eid, "state": state[x.name], "tier": "high"} for x in ann])
    assert sim_main(["score", "--calls", str(tmp_path / "calls.tsv"), "--elements", str(tmp_path / "S0.tsv"),
                     "--truth", str(tmp_path / "truth.tsv"), "--out", str(tmp_path / "s")]) == 0
    with open(tmp_path / "s.score.tsv") as fh:
        rows = {(r["state"], r["tier_set"]): r for r in csv.DictReader(fh, delimiter="\t")}
    for st in ("P", "E", "S"):
        assert (rows[(st, "all")]["precision"], rows[(st, "all")]["recall"]) == ("1.0000", "1.0000"), st
    assert rows[("summary", "all")]["n_scored"] == "4" and rows[("summary", "all")]["n_missing_call"] == "0"
