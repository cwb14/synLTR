"""Events on one chromosome's genealogy (synltr.sim.pop.events): primary and nested LTR-RT insertions, structural
losses, SNPs and small indels, checked against msprime's own genotypes and exact toy genealogies (poptoy.toy_ts)."""
from __future__ import annotations

import pytest

pytest.importorskip("msprime")

import math
import random

import msprime
import numpy as np
from sim import poptoy

from synltr.base.constants import MIN_SEG, SIM_FLANK_DEL, SIM_PARTIAL, SIM_PARTIAL_MIN, SLOP_TSD_MULT, TSD_DEFAULT
from synltr.sim.pop.donors import Donor
from synltr.sim.pop.events import ChromEvents, Ins, Loss, donor_body, draw_events, removed, thin_events
from synltr.sim.pop.model import Rates, ancestry

W = SLOP_TSD_MULT * TSD_DEFAULT
END = W + TSD_DEFAULT
ENDING = {"solo", "del_flank", "replace"}
KINDS = ("solo", "del_flank", "del_partial", "replace")
QUIET = {"mu": 0.0, "indel": 0.0}  # rates for tests of insertions and losses alone

# exact genealogies: 7 samples (total branch length 128,000) and 16 (218,000)
TREE7 = (40_000.0, (12_000.0, (3_000.0, 0, 1), (5_000.0, 2, 3)), (20_000.0, (8_000.0, 4, 5), 6))
TREE16 = (40_000.0, *[(25_000.0, *[(12_000.0, *[(5_000.0, 8 * a + 4 * b + 2 * c, 8 * a + 4 * b + 2 * c + 1)
                                                for c in range(2)]) for b in range(2)]) for a in range(2)])


def balanced(times: list[float], first: int = 0):
    """A balanced binary toy tree with 2 ** len(times) leaves from `first`; times[0] is the root's."""
    if not times:
        return first
    half = 2 ** (len(times) - 1)
    return (times[0], balanced(times[1:], first), balanced(times[1:], first + half))


TREE64 = balanced([45_000.0, 40_000.0, 33_000.0, 25_000.0, 12_000.0, 5_000.0])


def rng(seed: int) -> np.random.Generator:
    return np.random.default_rng(seed)


@pytest.fixture
def captured(monkeypatch):
    """The tree sequences msprime.sim_mutations returns during a test, keyed by model class name, plus
    "_discrete" for a discrete-genome call; "seeds": their random seeds, in call order."""
    out = {"seeds": []}
    real = msprime.sim_mutations

    def spy(ts, *args, **kw):
        m = real(ts, *args, **kw)
        key = type(kw["model"]).__name__ + ("_discrete" if kw.get("discrete_genome") else "")
        assert key not in out, f"two sim_mutations calls with {key}"
        assert ts.num_mutations == 0 and kw.get("keep") is False  # each call starts from no mutations
        out[key] = m
        out["seeds"].append(kw["random_seed"])
        return m

    monkeypatch.setattr(msprime, "sim_mutations", spy)
    return out


def same(a: ChromEvents, b: ChromEvents) -> bool:
    return (a.ins == b.ins and a.losses == b.losses and a.dropped == b.dropped
            and a.indel_carriers == b.indel_carriers
            and all(np.array_equal(getattr(a, k), getattr(b, k))
                    for k in ("snp_pos", "snp_base", "indel_pos", "indel_len", "indel_ins")))


def below(tree, node: int, time: float, top_node: int, top_time: float) -> bool:
    """Is the point (node, time) strictly below the point (top_node, top_time) on `tree` (each point on the
    branch above its node)?"""
    if node == top_node:
        return time < top_time
    return tree.is_descendant(node, top_node)


def samples_below(tree, node: int) -> frozenset[int]:
    index = {int(s): i for i, s in enumerate(tree.tree_sequence.samples())}
    return frozenset(index[int(s)] for s in tree.samples(node))


def iupac_block(n: int) -> bytes:
    return (b"RY" * n)[:n]


def test_carriers_match_genotypes(captured):
    ts = ancestry(poptoy.small_demog(), Rates(), poptoy.L, seed=1)
    bb, els = poptoy.toy_backbone(poptoy.L, [(20_000, 300, 5_000), (90_000, 500, 30_000)], seed=1)
    ev = draw_events(ts, bb, els, poptoy.toy_donors(), Rates(lam=4e-9, delta=0, **QUIET), rng(1))
    mts = captured["BinaryMutationModel"]
    geno: dict[int, list[frozenset[int]]] = {}
    for v in mts.variants():
        geno.setdefault(math.floor(v.site.position), []).append(frozenset(np.flatnonzero(v.genotypes == 1).tolist()))
    primary = [i for i in ev.ins if i.depth == 0]
    assert len(primary) > 100 and all(i.host is None for i in primary)
    for i in primary:
        assert i.carriers in geno[i.pos], i
    # every mutation is an insertion, or dropped near an end (none here) or on a non-ACGT target (none here)
    assert len(primary) == mts.num_mutations and ev.dropped["ins_end"] == ev.dropped["ins_n"] == 0
    # inside a backbone element by position only: host None, depth 0
    assert any(e.start < i.pos < e.end for i in primary for e in els)
    assert len(captured["seeds"]) == len(set(captured["seeds"])) == 3  # insertions, SNPs, indels: own seeds


def test_target_boundaries(monkeypatch):
    # msprime's insertion mutations replaced by one at every half base (x.5) of a 200-bp toy: an insertion is
    # kept exactly when 15 <= pos <= 185 (W + TSD_DEFAULT from each end) and [pos - 5, pos) holds only ACGT
    length = 200
    ts = poptoy.toy_ts(TREE7, length)
    bb = bytearray(poptoy.toy_backbone(length, seed=16)[0])
    bb[60] = ord("N")
    bb[100:103] = b"RYN"
    bb[150] = ord("k")  # lowercase IUPAC
    tables = ts.dump_tables()
    for p in range(length):
        tables.mutations.add_row(site=tables.sites.add_row(p + 0.5, "0"), node=p % 7, derived_state="1", time=0.0)
    planted = tables.tree_sequence()
    real = msprime.sim_mutations

    def fake(ts0, *args, **kw):
        return planted if kw.get("discrete_genome") is False else real(ts0, *args, **kw)

    monkeypatch.setattr(msprime, "sim_mutations", fake)
    ev = draw_events(ts, bytes(bb), [], poptoy.toy_donors(), Rates(lam=1e-12, **QUIET), rng(16))
    clean = {p for p in range(END, length - END + 1) if bytes(bb[p - TSD_DEFAULT:p]).upper().translate(None, b"ACGT")
             == b""}
    assert {i.pos for i in ev.ins} == clean and END == W + TSD_DEFAULT == 15
    assert ev.dropped["ins_end"] == 2 * END - 1 and ev.dropped["ins_n"] == length - 2 * END + 1 - len(clean)
    assert all(i.carriers == frozenset({i.pos % 7}) and i.time == 0.0 for i in ev.ins)


def test_no_insertion_target_on_n(captured):
    length, n0, n1 = 300_000, 100_000, 150_000
    bb, els = poptoy.toy_backbone(length, [(30_000, 300, 4_000)], seed=2)
    bb = bb[:n0] + b"N" * (n1 - n0) + bb[n1:]
    ts = poptoy.toy_ts(TREE7, length)
    ev = draw_events(ts, bb, els, poptoy.toy_donors(ltr=(100, 200), internal=(500, 1_000)),
                     Rates(lam=1.5e-8, delta=0, **QUIET), rng(2))
    primary = [i for i in ev.ins if i.depth == 0]
    assert primary and all(b"N" not in bb[i.pos - TSD_DEFAULT:i.pos] for i in primary)
    assert all(END <= i.pos <= length - END for i in primary)
    # dropped: exactly the mutations whose 5 target bases [pos - 5, pos) touch the N block (once clear of the ends)
    pos = np.floor(captured["BinaryMutationModel"].sites_position).astype(int)
    end = (pos < END) | (pos > length - END)
    on_n = int(np.sum(~end & (pos - TSD_DEFAULT < n1) & (pos > n0)))
    assert ev.dropped["ins_n"] == on_n > 0 and ev.dropped["ins_end"] == int(end.sum())
    assert len(primary) == len(pos) - on_n - int(end.sum())


def test_backbone_case_and_iupac(captured):
    length, r0, r1 = 300_000, 120_000, 140_000
    bb, els = poptoy.toy_backbone(length, [(30_000, 300, 4_000), (200_000, 400, 20_000)], seed=3)
    bb = bb[:r0] + iupac_block(r1 - r0) + bb[r1:]
    ts = poptoy.toy_ts(TREE7, length)
    donors = poptoy.toy_donors(ltr=(100, 200), internal=(500, 1_000))
    rates = Rates(lam=1.5e-8, delta=2e-5, mu=1e-6, indel=1e-7)
    upper = draw_events(ts, bb, els, donors, rates, rng(3))
    sites = np.floor(captured.pop("JC69_discrete").sites_position).astype(int)
    captured.clear()
    captured["seeds"] = []
    lower = draw_events(ts, bb.lower(), els, donors, rates, rng(3))
    assert same(upper, lower)  # soft-masked lowercase reads as the same bases
    ev = upper
    # no SNP on an IUPAC base: those sites are dropped, every other site kept
    assert not np.any((ev.snp_pos >= r0) & (ev.snp_pos < r1)) and len(ev.snp_pos)
    assert ev.dropped["snp_n"] == int(np.sum((sites >= r0) & (sites < r1))) > 0
    assert len(ev.snp_pos) == len(sites) - ev.dropped["snp_n"]
    assert set(np.unique(ev.snp_base).tolist()) <= set(b"ACGT")
    # no target read across the block
    primary = [i for i in ev.ins if i.depth == 0]
    assert primary and not any(i.pos - TSD_DEFAULT < r1 and i.pos > r0 for i in primary)
    assert ev.dropped["ins_n"] > 0


def test_snp_ancestral_is_backbone(captured):
    ts = ancestry(poptoy.small_demog(), Rates(), poptoy.L, seed=2)
    bb, els = poptoy.toy_backbone(poptoy.L, [(50_000, 300, 6_000)], seed=4)
    ev = draw_events(ts, bb, els, [], Rates(lam=0, delta=0, indel=0), rng(4))
    jts = captured["JC69_discrete"]
    assert len(ev.snp_pos) == jts.num_sites > 100 and ev.dropped["snp_n"] == 0
    assert ev.snp_base.shape == (jts.num_sites, ts.num_samples) and ev.snp_base.dtype == np.uint8
    order = "ACGT"
    tree = jts.first()
    for k, var in enumerate(jts.variants()):
        site = var.site
        pos = int(site.position)
        assert ev.snp_pos[k] == pos
        tree.seek(pos)
        hits = [0] * ts.num_samples
        for m in site.mutations:
            for s in samples_below(tree, m.node):
                hits[s] += 1
        base = chr(bb[pos])
        for s in range(ts.num_samples):
            got = chr(ev.snp_base[k, s])
            if hits[s] == 0:
                assert got == base  # on no mutated lineage: the backbone base
            elif hits[s] == 1:
                assert got != base  # JC69 always changes the base
        # the relabeling: the permutation of ACGT sending msprime's ancestral state to the backbone base and
        # keeping the order of the other three
        anc = site.ancestral_state
        perm = {anc: base, **dict(zip([b for b in order if b != anc], [b for b in order if b != base]))}
        assert [chr(x) for x in ev.snp_base[k]] == [perm[var.alleles[g]] for g in var.genotypes]


def test_loss_split():
    # 80 elements 600 bp apart: a del_flank (up to 2,000 bp) reaches into a neighbor but never over its whole 5.6-kb
    # span, so by the whole-span rule every element exists wherever it was drawn and every drawn loss is kept
    ts = poptoy.toy_ts(TREE64, 500_000)
    bb, _ = poptoy.toy_backbone(500_000, seed=5)
    els = [poptoy.toy_element(f"bb{i:04d}", 1_000 + 6_200 * i, 300, 5_000) for i in range(80)]
    assert all(e.ltr_right[0] - e.ltr_left[1] >= SIM_PARTIAL_MIN for e in els)
    split = Rates().loss_split
    ev = draw_events(ts, bb, els, [], Rates(lam=0, delta=1e-4, **QUIET), rng(5))
    n = len(ev.losses)
    assert n > 1000 and ev.dropped["loss_after_loss"] == 0
    for kind, p in zip(KINDS, split):
        got = sum(x.kind == kind for x in ev.losses) / n
        assert abs(got - p) <= 3 * math.sqrt(p * (1 - p) / n), (kind, got, p)


def test_loss_ends_lineage():
    length = 300_000
    ts = poptoy.toy_ts(TREE16, length)
    tree = ts.first()
    bb, els = poptoy.toy_backbone(length, [(10_000 + 20_000 * i, 300, 6_000) for i in range(12)], seed=6)
    ev = draw_events(ts, bb, els, poptoy.toy_donors(), Rates(lam=3e-9, delta=1e-4, **QUIET), rng(6))
    by: dict[str, list] = {}
    for x in ev.losses:
        by.setdefault(x.element, []).append(x)
    assert set(by) & {e.eid for e in els} and set(by) & {i.iid for i in ev.ins}  # both kinds of element lose
    after_partial = after_end_elsewhere = 0
    for xs in by.values():
        for end in (x for x in xs if x.kind in ENDING):
            for x in xs:
                assert x is end or not below(tree, x.node, x.time, end.node, end.time), (end, x)
                after_end_elsewhere += x is not end
        for p in (x for x in xs if x.kind == "del_partial"):
            after_partial += any(below(tree, x.node, x.time, p.node, p.time) for x in xs)
    assert after_partial > 0 and after_end_elsewhere > 0  # del_partial does not end the element; others go on


def test_replace_length():
    ts = poptoy.toy_ts(TREE16, 400_000)
    bb, _ = poptoy.toy_backbone(400_000, seed=7)
    els = [poptoy.toy_element(f"bb{i:04d}", 1_000 + 900 * i, 300, 1_000) for i in range(400)]
    rates = Rates(lam=0, delta=1e-4, loss_split=(0, 0, 0, 1), **QUIET)
    ev = draw_events(ts, bb, els, [], rates, rng(7))
    lengths = np.array([x.a for x in ev.losses])
    assert len(lengths) > 800 and all(x.kind == "replace" and x.b == 0 for x in ev.losses)
    assert lengths.min() >= MIN_SEG + 1
    se = lengths.std(ddof=1) / math.sqrt(len(lengths))
    assert abs(lengths.mean() - rates.replace_mean) <= 3 * se, (lengths.mean(), se)
    ev = draw_events(ts, bb, els, [], Rates(lam=0, delta=1e-4, **QUIET), rng(8))  # the default mix too
    assert all(x.a >= MIN_SEG + 1 for x in ev.losses if x.kind == "replace")


def test_loss_parameters():
    length = 300_000
    ts = poptoy.toy_ts(TREE16, length)
    tree = ts.first()
    bb, els = poptoy.toy_backbone(length, [(5_000 + 25_000 * i, 200 + 50 * i, 3_000 + 1_000 * i)
                                           for i in range(11)] + [(285_000, 300, 800)], seed=9)
    donors = poptoy.toy_donors()
    ev = draw_events(ts, bb, els, donors, Rates(lam=3e-9, delta=1e-4, **QUIET), rng(9))
    assert any(x.element == els[-1].eid for x in ev.losses)  # its 800-bp internal region takes no del_partial
    shape = {e.eid: (min(e.ltr_left[1] - e.ltr_left[0], e.ltr_right[1] - e.ltr_right[0]),
                     e.ltr_right[0] - e.ltr_left[1]) for e in els}
    for i in ev.ins:
        d = donors[i.donor]
        shape[i.iid] = (d.ltr_left, len(d.seq) - d.ltr_left - d.ltr_right)
    origin = {i.iid: (i.node, i.time) for i in ev.ins}
    flank, partial_after = set(), 0
    for x in ev.losses:
        ltr, internal = shape[x.element]
        assert x.carriers == samples_below(tree, x.node)
        assert tree.time(x.node) <= x.time <= tree.time(tree.parent(x.node))
        if x.element in origin:
            assert below(tree, x.node, x.time, *origin[x.element])
        if x.kind == "solo":
            assert 1 <= x.a <= ltr - 1 and x.b == 0
        elif x.kind == "del_flank":
            assert (x.a, x.b) != (0, 0) and all(v == 0 or SIM_FLANK_DEL[0] <= v <= SIM_FLANK_DEL[1] for v in (x.a, x.b))
            flank.add((x.a > 0, x.b > 0))
        elif x.kind == "del_partial":  # [a, a + b) of the original internal region, 20-60% of it
            assert internal >= SIM_PARTIAL_MIN and 0 <= x.a and x.a + x.b <= internal
            assert SIM_PARTIAL[0] * internal - 1 <= x.b <= SIM_PARTIAL[1] * internal + 1
            partial_after += any(y.element == x.element and y.kind == "del_partial" and y is not x
                                 and below(tree, x.node, x.time, y.node, y.time) for y in ev.losses)
        else:
            assert x.kind == "replace" and x.a >= MIN_SEG + 1 and x.b == 0
    assert flank == {(True, False), (False, True), (True, True)}  # left, right, both
    assert {x.kind for x in ev.losses} == set(KINDS) and partial_after > 0


def test_nested_carriers_and_offsets():
    length = 300_000
    ts = poptoy.toy_ts(TREE16, length)
    tree = ts.first()
    bb, els = poptoy.toy_backbone(length, [(50_000, 300, 6_000)], seed=10)
    donors = poptoy.toy_donors()
    for seed, rates in ((10, Rates(lam=6e-9, delta=2e-5, **QUIET)),
                        (11, Rates(lam=6e-9, delta=1e-4, loss_split=(0.1, 0.1, 0.7, 0.1), **QUIET))):
        ev = draw_events(ts, bb, els, donors, rates, rng(seed))
        host = {i.iid: i for i in ev.ins}
        losses: dict[str, list] = {}
        for x in ev.losses:
            losses.setdefault(x.element, []).append(x)
        nested = [i for i in ev.ins if i.host is not None]
        assert {i.depth for i in nested} == {1, 2}
        assert all(type(i.time) is float for i in ev.ins) and all(type(x.time) is float for x in ev.losses)
        for i in nested:
            h = host[i.host]
            assert i.depth == h.depth + 1 and i.iid != h.iid
            assert TSD_DEFAULT <= i.pos < len(donor_body(donors[h.donor]))  # its target in the host's body
            assert i.carriers <= h.carriers and i.carriers == samples_below(tree, i.node)
            assert below(tree, i.node, i.time, h.node, h.time)  # below the host's own event
            # never below an ending loss of its direct host (B_X stops there; a host's host's ending loss may
            # leave the host's position, and the host with it)
            assert not any(x.kind in ENDING and below(tree, i.node, i.time, x.node, x.time)
                           for x in losses.get(h.iid, [])), (i, h)
            # no target base in what an older kept loss removed from the host's body
            for x in ev.losses:
                if below(tree, i.node, i.time, x.node, x.time):
                    for system, lo, hi in removed(x, {e.eid: e for e in els}, host, donors, length):
                        assert not (system == h.iid and lo < i.pos and i.pos - TSD_DEFAULT < hi), (i, x)
        if seed == 11:
            assert ev.dropped["nested_after_loss"] > 0


def test_nested_target_not_on_n():
    # donors made of a 7-bp period with one non-ACGT base (N, or a lowercase IUPAC r): a nested insertion's target
    # [pos - 5, pos) of its host's body is clean, so offsets = 5 or 6 (mod 7) are kept and the rest dropped; at 6 the
    # next base is the N, so the check is one-sided whatever the seed
    length = 300_000
    ts = poptoy.toy_ts(TREE16, length)
    bb, els = poptoy.toy_backbone(length, seed=17)
    donors = [Donor(unit * 900, 700, 700) for unit in ("ACGTACN", "ACGTACr")]  # body 6,300 = the seq
    assert all(len(donor_body(d)) == 6_300 for d in donors)
    ev = draw_events(ts, bb, els, donors, Rates(lam=6e-9, delta=0, **QUIET), rng(17))
    host = {i.iid: i for i in ev.ins}
    nested = [i for i in ev.ins if i.host is not None]
    assert len(nested) > 50
    for i in nested:
        body = donor_body(donors[host[i.host].donor])
        assert not body[i.pos - TSD_DEFAULT:i.pos].translate(None, b"ACGT"), i
        assert i.pos % 7 in (5, 6)
    assert any(i.pos % 7 == 6 for i in nested)  # its next base is non-ACGT, which the target does not read
    assert ev.dropped["ins_n"] > 0  # all from nested ones: the backbone is ACGT


def test_nested_rate():
    # delta = 0, so nothing is thinned: the depth-1 count against sum over hosts of lam * (len(body) - 5) * B_X, B_X
    # the branch length below the host's event; donors with LTRs far apart in length, so len(body) differs from
    # len(seq) by half
    length, lam = 300_000, 6e-9
    ts = poptoy.toy_ts(TREE16, length)
    tree = ts.first()
    bb, els = poptoy.toy_backbone(length, seed=19)
    donors = [Donor(poptoy.rand_bases(5_100, random.Random(k)).decode(), 3_000, 100) for k in range(3)]
    assert len(donor_body(donors[0])) == 8_000

    def below_point(node: int, t: float) -> float:
        return (t - tree.time(node)) + sum(tree.branch_length(u) for u in tree.nodes(node) if u != node)

    obs, exp = {1: 0, 2: 0}, {1: 0.0, 2: 0.0}
    for seed in range(19, 24):
        ev = draw_events(ts, bb, els, donors, Rates(lam=lam, delta=0, **QUIET), rng(seed))
        assert ev.dropped["nested_after_loss"] == ev.dropped["ins_n"] == 0
        for x in ev.ins:
            if x.depth < 2:
                exp[x.depth + 1] += lam * (len(donor_body(donors[x.donor])) - TSD_DEFAULT) * below_point(x.node, x.time)
            if x.depth:
                obs[x.depth] += 1
    for d in (1, 2):
        assert exp[d] > 200 and abs(obs[d] - exp[d]) <= 3 * math.sqrt(exp[d]), (d, obs[d], exp[d])


def test_insertion_count_poisson():
    length, lam = 1_000_000, 3e-9
    ts = poptoy.toy_ts(TREE7, length)
    expected = lam * length * poptoy.tree_length(TREE7)
    assert expected == pytest.approx(384)
    bb, els = poptoy.toy_backbone(length, seed=12)
    ev = draw_events(ts, bb, els, poptoy.toy_donors(ltr=(100, 200), internal=(500, 1_000)),
                     Rates(lam=lam, delta=0, **QUIET), rng(12))
    n = sum(i.depth == 0 for i in ev.ins) + ev.dropped["ins_end"]
    assert abs(n - expected) <= 3 * math.sqrt(expected), (n, expected)


def test_indels(captured):
    length = 300_000
    ts = poptoy.toy_ts(TREE7, length)
    bb, els = poptoy.toy_backbone(length, seed=13)
    ev = draw_events(ts, bb, els, [], Rates(lam=0, delta=0, mu=0, indel=1e-7), rng(13))
    its = captured["BinaryMutationModel_discrete"]
    m = its.num_mutations
    assert m > 1000 and len(ev.indel_pos) == len(ev.indel_len) == len(ev.indel_ins) == len(ev.indel_carriers) == m
    tree = ts.first()
    assert ev.indel_carriers == [samples_below(tree, int(u)) for u in its.mutations_node]
    assert ev.indel_pos.tolist() == its.sites_position[its.mutations_site].astype(int).tolist()
    assert ev.indel_len.min() >= 1 and np.all(ev.indel_pos + ev.indel_len * ~ev.indel_ins <= length)
    share = ev.indel_ins.mean()
    assert abs(share - 0.5) <= 3 * math.sqrt(0.25 / m)
    mean, var = 3.0, 6.0  # geometric, p = 1/3
    assert abs(ev.indel_len.mean() - mean) <= 3 * math.sqrt(var / m)


def test_seed_determinism():
    length = 300_000
    ts = ancestry(poptoy.small_demog(), Rates(), length, seed=3)
    bb, els = poptoy.toy_backbone(length, [(10_000 + 30_000 * i, 300, 6_000) for i in range(8)], seed=14)
    donors = poptoy.toy_donors()
    rates = Rates(lam=3e-9, delta=1e-5, mu=1e-7, indel=1e-8)
    a, b = (draw_events(ts, bb, els, donors, rates, rng(14)) for _ in range(2))
    assert a.ins and a.losses and len(a.snp_pos) and len(a.indel_pos)
    assert same(a, b)
    assert not same(a, draw_events(ts, bb, els, donors, rates, rng(15)))


def test_bad_input_fails():
    ts = poptoy.toy_ts(TREE7, 1_000)
    bb, _ = poptoy.toy_backbone(1_000, seed=0)
    with pytest.raises(ValueError, match="length"):
        draw_events(ts, bb[:999], [], [], Rates(), rng(0))
    with pytest.raises(ValueError, match="donor"):
        draw_events(ts, bb, [], [], Rates(lam=1e-6), rng(0))
    with pytest.raises(ValueError, match="bb0001"):
        draw_events(ts, bb, [poptoy.toy_element("bb0001", 900, 100, 100)], [], Rates(lam=0), rng(0))


# --------------------------------------------------------------- existence (thin_events, removed)

TOY_L = 300_000
TOY_DONOR = Donor("".join("ACGT"[(7 * i) % 4] for i in range(5_000)), 500, 500)  # body 5,000 bp, LTRs 500
E1 = poptoy.toy_element("bb0001", 10_000, 300, 4_000)  # span [10,000, 14,600), internal [10,300, 14,300)


def toy7():
    """TREE7 on TOY_L, its tree, and its internal nodes by the samples below them."""
    ts = poptoy.toy_ts(TREE7, TOY_L)
    tree = ts.first()
    n = {"01": tree.parent(0), "23": tree.parent(2), "45": tree.parent(4)}
    n["0123"], n["456"] = tree.parent(n["01"]), tree.parent(n["45"])
    assert [tree.time(n[k]) for k in ("01", "23", "45", "0123", "456")] == [3_000, 5_000, 8_000, 12_000, 20_000]
    return tree, n


def ins_at(tree, k: int, host: str | None, pos: int, node: int, time: float, depth: int = 0,
           strand: str = "+") -> Ins:
    return Ins(f"i00{k:06d}", host, pos, node, float(time), 0, strand, depth, samples_below(tree, node))


def loss_at(tree, element: str, node: int, time: float, kind: str, a: int, b: int = 0) -> Loss:
    return Loss(element, node, float(time), kind, a, b, samples_below(tree, node))


def thin(ins, losses, elements=(E1,)):
    return thin_events(list(ins), list(losses), list(elements), [TOY_DONOR], TOY_L)


def counts(ins_after=0, nested_after=0, loss_after=0) -> dict:
    return {"ins_after_loss": ins_after, "nested_after_loss": nested_after, "loss_after_loss": loss_after}


def test_thin_primary_below_solo():
    tree, n = toy7()
    solo = loss_at(tree, "bb0001", n["0123"], 20_000, "solo", 100)  # on the branch above samples 0-3
    assert removed(solo, {"bb0001": E1}, {}, [TOY_DONOR], TOY_L) == [(None, 10_100, 14_400)]
    below = ins_at(tree, 1, None, 12_000, n["01"], 5_000)  # target in E1's internal region, below the solo
    sister = ins_at(tree, 2, None, 12_000, n["45"], 5_000)  # the same insertion on the sister clade
    older = ins_at(tree, 3, None, 12_000, n["0123"], 30_000)  # above the solo on its branch
    right = ins_at(tree, 4, None, 14_405, n["01"], 5_000)  # target [14,400, 14,405): just past what it removed
    edge = ins_at(tree, 5, None, 14_404, n["01"], 5_000)  # target [14,399, 14,404): one base removed
    ins, losses, dropped = thin([below, sister, older, right, edge], [solo])
    assert ins == [sister, older, right] and losses == [solo] and dropped == counts(ins_after=2)


def test_thin_older_insertion_younger_solo():
    tree, n = toy7()
    x = ins_at(tree, 1, None, 12_000, n["0123"], 30_000)  # inside E1's internal region, older than its solo
    solo = loss_at(tree, "bb0001", n["0123"], 20_000, "solo", 100)
    x_below = loss_at(tree, x.iid, n["01"], 2_000, "replace", 500)  # below the solo: x is gone there
    x_above = loss_at(tree, x.iid, n["0123"], 25_000, "del_partial", 2_000, 1_500)  # between x and the solo
    x_sister = loss_at(tree, x.iid, 2, 1_000, "solo", 50)  # sample 2 is below the solo too
    # at the removed stretch's edges [10,100, 14,400) an insertion sits at a junction and stays
    at_lo = ins_at(tree, 2, None, 10_100, n["0123"], 30_000)
    at_hi = ins_at(tree, 3, None, 14_400, n["0123"], 30_000)
    edge_losses = [loss_at(tree, i.iid, 0, 1_000, "replace", 500) for i in (at_lo, at_hi)]
    # y, nested in x before the solo, goes with x: its own anchor in x's body is intact, but its host is gone
    y = ins_at(tree, 4, x.iid, 1_000, n["0123"], 28_000, depth=1)
    y_above = loss_at(tree, y.iid, n["0123"], 22_000, "del_partial", 100, 500)  # above the solo
    y_below = loss_at(tree, y.iid, 1, 1_000, "replace", 500)
    ins, losses, dropped = thin([x, at_lo, at_hi, y], [solo, x_below, x_above, x_sister, *edge_losses, y_above,
                                                       y_below])
    assert ins == [x, at_lo, at_hi, y] and losses == [solo, x_above, *edge_losses, y_above]
    assert dropped == counts(loss_after=3)
    assert [(s, lo, hi) for s, lo, hi in removed(solo, {"bb0001": E1}, {}, [TOY_DONOR], TOY_L)
            if s is None and lo < x.pos < hi] == [(None, 10_100, 14_400)]  # it covers x's position


def depth2_case(tree, n):
    """x primary, its del_flank (ending) on the branch above samples 0-1, y nested in x above that loss."""
    x = ins_at(tree, 1, None, 50_000, n["0123"], 35_000)
    flank = loss_at(tree, x.iid, n["01"], 10_000, "del_flank", 300)
    y = ins_at(tree, 2, x.iid, 1_000, n["0123"], 30_000, depth=1)
    return x, flank, y


def test_thin_depth2_below_host_of_host():
    tree, n = toy7()
    x, flank, y = depth2_case(tree, n)
    z_below = ins_at(tree, 3, y.iid, 500, n["01"], 5_000, depth=2)  # below x's del_flank: y went with x
    z_sister = ins_at(tree, 4, y.iid, 500, n["23"], 4_000, depth=2)
    ins, losses, dropped = thin([x, y, z_below, z_sister], [flank])
    assert ins == [x, y, z_sister] and losses == [flank] and dropped == counts(nested_after=1)


def test_thin_nested_loss_below_host_loss():
    tree, n = toy7()
    x, flank, y = depth2_case(tree, n)
    y_below = loss_at(tree, y.iid, 0, 1_000, "solo", 50)  # sample 0 is below x's del_flank
    y_sister = loss_at(tree, y.iid, 2, 1_000, "solo", 50)
    ins, losses, dropped = thin([x, y], [flank, y_below, y_sister])
    assert ins == [x, y] and losses == [flank, y_sister] and dropped == counts(loss_after=1)


def test_thin_nested_in_solo_ltr():
    # x's solo (crossover 300) removes [300, 4,800) of its body. y sits in x's left LTR at 50, which the solo
    # leaves: y stays present below it and keeps its own younger losses and nested insertions there (a host's
    # ending loss ends only what it removes). y2 sits at 1,000, in the removed stretch: gone below the solo. A new
    # insertion into x itself below its solo is discarded: x no longer exists as an element there.
    tree, n = toy7()
    x = ins_at(tree, 1, None, 50_000, n["0123"], 35_000)
    y = ins_at(tree, 2, x.iid, 50, n["0123"], 30_000, depth=1)
    y2 = ins_at(tree, 3, x.iid, 1_000, n["0123"], 30_000, depth=1)
    solo = loss_at(tree, x.iid, n["01"], 10_000, "solo", 300)
    assert removed(solo, {}, {x.iid: x}, [TOY_DONOR], TOY_L) == [(x.iid, 300, 4_800)]  # 50 is not in it
    y_below, y_sister = (loss_at(tree, y.iid, s, 1_000, "replace", 500) for s in (0, 2))
    y2_below, y2_sister = (loss_at(tree, y2.iid, s, 1_000, "replace", 500) for s in (0, 2))
    z_below = ins_at(tree, 4, y.iid, 400, 1, 1_000, depth=2)  # into y, below x's solo
    x_new = ins_at(tree, 5, x.iid, 60, 1, 1_000, depth=1)  # into x's solo LTR: x has ended there
    ins, losses, dropped = thin([x, y, y2, z_below, x_new], [solo, y_below, y_sister, y2_below, y2_sister])
    assert ins == [x, y, y2, z_below] and losses == [solo, y_below, y_sister, y2_sister]
    assert dropped == counts(nested_after=1, loss_after=1)


def test_thin_same_position_del_flank():
    # the junction rule: w, anchored before backbone base 100,000 like x and younger than x, goes with x's del_flank
    # only when it takes both sides of that position; a one-sided del_flank or x's replace leaves it
    tree, n = toy7()
    x = ins_at(tree, 1, None, 100_000, n["0123"], 30_000)
    w = ins_at(tree, 2, None, 100_000, n["0123"], 25_000)
    w_loss = loss_at(tree, w.iid, 0, 1_000, "replace", 500)
    cases = {
        "two-sided": (loss_at(tree, x.iid, n["01"], 10_000, "del_flank", 1_000, 500), [], 1),
        "end side only": (loss_at(tree, x.iid, n["01"], 10_000, "del_flank", 0, 500), [w_loss], 0),
        "start side only": (loss_at(tree, x.iid, n["01"], 10_000, "del_flank", 1_000, 0), [w_loss], 0),
        "replace": (loss_at(tree, x.iid, n["01"], 10_000, "replace", 800), [w_loss], 0),
    }
    for name, (xl, kept_w, gone) in cases.items():
        ins, losses, dropped = thin([x, w], [xl, w_loss], ())
        assert ins == [x, w] and losses == [xl, *kept_w] and dropped == counts(loss_after=gone), name


def test_thin_del_flank_flanks():
    tree, n = toy7()
    x = ins_at(tree, 1, None, 100_000, n["0123"], 30_000)
    flank = loss_at(tree, x.iid, n["0123"], 20_000, "del_flank", 1_000, 500)
    # x's body and its second TSD copy; 1,000 bp before it (its target bases first); 500 after it, the TSD copy's
    # 5 included
    assert removed(flank, {}, {x.iid: x}, [TOY_DONOR], TOY_L) == \
        [(x.iid, 0, 5_005), (None, 99_000, 100_000), (None, 100_000, 100_495)]
    left_below = ins_at(tree, 2, None, 99_500, n["01"], 5_000)
    left_above = ins_at(tree, 3, None, 99_500, n["0123"], 25_000)  # older than the del_flank
    right_below = ins_at(tree, 4, None, 100_300, n["23"], 4_000)
    past_right = ins_at(tree, 5, None, 100_500, n["01"], 5_000)  # target [100,495, 100,500)
    before_left = ins_at(tree, 6, None, 99_000, n["01"], 5_000)  # target [98,995, 99,000)
    ins, losses, dropped = thin([x, left_below, left_above, right_below, past_right, before_left], [flank], ())
    assert ins == [x, left_above, past_right, before_left] and losses == [flank] and dropped == counts(ins_after=2)


def test_removed_flank_runs_into_host():
    # a del_flank of a nested insertion near its host's body end goes on through the host's second TSD copy into
    # the host's own host; on strand - the host's body end lies on its host's start side
    tree, n = toy7()
    # 500 bp past y's end: y's TSD copy (5), x's last 100 body bases, then on +: x's TSD copy (5) and 390 bp of
    # backbone after x; on -: 395 bp of backbone before x
    for strand, out in (("+", [(None, 100_000, 100_390)]), ("-", [(None, 99_605, 100_000)])):
        x = ins_at(tree, 1, None, 100_000, n["0123"], 30_000, strand=strand)
        y = ins_at(tree, 2, x.iid, 4_900, n["0123"], 25_000, depth=1)  # 100 bp from x's body end
        flank = loss_at(tree, y.iid, n["01"], 10_000, "del_flank", 0, 500)
        got = removed(flank, {}, {x.iid: x, y.iid: y}, [TOY_DONOR], TOY_L)
        tsd = [(x.iid, 5_000, 5_005)] if strand == "+" else []
        assert got == [(y.iid, 0, 5_005), (x.iid, 4_900, 5_000)] + tsd + out, strand
    # a backbone element's del_flank is clipped at the backbone's ends; solo, del_partial and replace of an insertion
    # are in its body
    e = poptoy.toy_element("bb0002", 500, 300, 2_000)
    assert removed(loss_at(tree, "bb0002", 0, 10, "del_flank", 1_000, 2_000), {"bb0002": e}, {}, [TOY_DONOR],
                   TOY_L) == [(None, 0, 5_100)]
    x = ins_at(tree, 1, None, 100_000, n["0123"], 30_000)
    m = {x.iid: x}
    assert removed(loss_at(tree, x.iid, 0, 10, "solo", 40), {}, m, [TOY_DONOR], TOY_L) == [(x.iid, 40, 4_540)]
    assert removed(loss_at(tree, x.iid, 0, 10, "del_partial", 100, 800), {}, m, [TOY_DONOR], TOY_L) == \
        [(x.iid, 600, 1_400)]
    assert removed(loss_at(tree, x.iid, 0, 10, "replace", 9_000), {}, m, [TOY_DONOR], TOY_L) == [(x.iid, 0, 5_000)]


def test_thin_overlapping_del_partials():
    tree, n = toy7()
    p1 = loss_at(tree, "bb0001", n["0123"], 20_000, "del_partial", 500, 1_500)
    p2 = loss_at(tree, "bb0001", n["01"], 5_000, "del_partial", 1_200, 1_600)
    got = [removed(p, {"bb0001": E1}, {}, [TOY_DONOR], TOY_L) for p in (p1, p2)]
    assert got == [[(None, 10_800, 12_300)], [(None, 11_500, 13_100)]]
    union = np.zeros(TOY_L, dtype=bool)  # on the lineage of samples 0 and 1, below both
    for (_, lo, hi), in got:
        union[lo:hi] = True
    assert np.flatnonzero(np.diff(union.astype(int))).tolist() == [10_800 - 1, 13_100 - 1]  # [10,800, 13,100)
    only_p1 = ins_at(tree, 1, None, 12_800, n["23"], 4_000)  # target [12,795, 12,800): p2 alone removed it
    both = ins_at(tree, 2, None, 12_800, 0, 1_000)
    tail = ins_at(tree, 3, None, 13_102, 0, 1_000)  # target [13,097, 13,102) reaches into p2's
    past = ins_at(tree, 4, None, 13_105, 0, 1_000)
    ins, losses, dropped = thin([only_p1, both, tail, past], [p1, p2])
    assert ins == [only_p1, past] and losses == [p1, p2] and dropped == counts(ins_after=2)


def test_thin_backbone_element_inside_removed_one():
    # a backbone element nested in another's internal region is gone below the host element's solo, not below a
    # del_partial that takes only part of it
    tree, n = toy7()
    outer = poptoy.toy_element("bb0001", 10_000, 300, 20_000)  # internal [10,300, 30,300)
    inner = poptoy.toy_element("bb0002", 15_000, 200, 1_000)  # [15,000, 16,400)
    solo = loss_at(tree, "bb0001", n["0123"], 20_000, "solo", 100)
    part = loss_at(tree, "bb0001", n["456"], 25_000, "del_partial", 5_000, 1_000)  # [15,300, 16,300) of 'inner'
    inner_below_solo = loss_at(tree, "bb0002", n["01"], 2_000, "replace", 300)
    inner_below_part = loss_at(tree, "bb0002", n["45"], 2_000, "replace", 300)
    ins, losses, dropped = thin([], [solo, part, inner_below_solo, inner_below_part], (outer, inner))
    assert ins == [] and losses == [solo, part, inner_below_part] and dropped == counts(loss_after=1)


def test_thin_density():
    # insertions drawn without losses, then fixed replacements: E1 (40% of the backbone) above samples 0-3 from
    # 39,000 generations, E2 above sample 4 from 6,000. Kept primaries ~ Poisson(lam * sum over branches of
    # present positions x branch length); a position is lost when its 5 target bases touch the removed span
    length, lam = 1_000_000, 1e-8
    ts = poptoy.toy_ts(TREE7, length)
    tree = ts.first()
    big = poptoy.toy_element("bb0001", 200_000, 1_000, 398_000)  # [200,000, 600,000)
    small = poptoy.toy_element("bb0002", 700_000, 1_000, 98_000)  # [700,000, 800,000)
    bb, _ = poptoy.toy_backbone(length, seed=18)
    donors = poptoy.toy_donors(ltr=(100, 200), internal=(500, 1_000))
    ev = draw_events(ts, bb, [big, small], donors, Rates(lam=lam, delta=0, **QUIET), rng(18))
    n01 = tree.parent(0)
    n0123 = tree.parent(n01)
    fixed = [loss_at(tree, "bb0001", n0123, 39_000, "replace", 13_538),
             loss_at(tree, "bb0002", 4, 6_000, "replace", 13_538)]
    ins, losses, dropped = thin_events(ev.ins, fixed, [big, small], donors, length)
    below_big = (39_000 - 12_000) + sum(poptoy.tree_length(k) + (12_000 - (0 if isinstance(k, int) else k[0]))
                                        for k in TREE7[1][1:])  # 27,000 + 32,000
    assert below_big == 59_000
    lost = 400_004 * below_big + 100_004 * 6_000  # positions p in (start, end + 5)
    expected = lam * (length * poptoy.tree_length(TREE7) - lost)
    kept = sum(i.depth == 0 for i in ins)
    drawn = sum(i.depth == 0 for i in ev.ins)
    assert abs(kept - expected) <= 3 * math.sqrt(expected), (kept, expected)
    assert dropped["ins_after_loss"] == drawn - kept > 0 and losses == fixed
    # without the thinning the count would sit far above: lam * L * T
    assert lam * length * poptoy.tree_length(TREE7) - expected > 6 * math.sqrt(expected)
