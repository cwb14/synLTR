"""Population model of `synltr sim pop` (synltr.sim.pop.model): the msprime demography (three taxa A, B and C,
an outgroup, optional lineage populations, a pulse, a bottleneck), the ancestry of one haploid genome per
sample, truth tracts from the pulse's migration records, the recombination map and the --demography JSON."""
from __future__ import annotations

import pytest

pytest.importorskip("msprime")

import json
import os
import subprocess
import sys
from dataclasses import asdict, dataclass, fields
from pathlib import Path

import numpy as np
from sim import poptoy

from synltr.sim.pop.model import (
    MSPRIME_HELP,
    Demog,
    Rates,
    ancestry,
    demog_from_json,
    demography,
    flat_fields,
    recombination_map,
    sample_labels,
    tracts,
)

L = poptoy.L
PULSES = ("C>A", "A>C", "C>A1", "C1>A1", "A1>C1")
LINEAGES = {"sub_a": (2000, 1000), "sub_c": (30000, 1000)}  # A1p splits from A at 2,000, C1p from C at 30,000


def pop_names(ts) -> list[str]:
    """The population of each sample node, by name."""
    return [ts.population(ts.node(s).population).metadata["name"] for s in ts.samples()]


def merge(intervals) -> list[tuple[int, int]]:
    """Sorted, disjoint, non-adjacent int intervals: the union of `intervals`."""
    out: list[tuple[int, int]] = []
    for lo, hi in sorted((int(a), int(b)) for a, b in intervals):
        if out and lo <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], hi))
        else:
            out.append((lo, hi))
    return out


def pulse_records(ts, d: Demog) -> np.ndarray:
    """Indices of the migration records of the pulse: time ta, recipient to donor (backward in time)."""
    donor, recipient = {"C>A": ("C", "A"), "A>C": ("A", "C"), "C>A1": ("C", "A1p"), "C1>A1": ("C1p", "A1p"),
                        "A1>C1": ("A1p", "C1p")}[d.pulse]
    pid = {ts.population(i).metadata["name"]: i for i in range(ts.num_populations)}
    return np.flatnonzero((ts.migrations_source == pid[recipient]) & (ts.migrations_dest == pid[donor])
                          & (ts.migrations_time == d.ta))


def test_sample_times_zero():
    d = poptoy.small_demog()
    assert sample_labels(d) == [("A1", "A"), ("A2", "A"), ("B1", "B"), ("B2", "B"), ("C1", "C"), ("C2", "C"),
                                ("O1", "O")]
    ts = ancestry(d, Rates(), L, seed=1)
    assert ts.num_samples == 7 and ts.sequence_length == L
    assert all(ts.node(s).time == 0 for s in ts.samples())
    assert pop_names(ts) == ["A", "A", "B", "B", "C", "C", "O"]
    # A1 and C1 are sampled from their own lineage populations when those exist
    d = poptoy.small_demog(**LINEAGES)
    ts = ancestry(d, Rates(), L, seed=1)
    assert pop_names(ts) == ["A1p", "A", "B", "B", "C1p", "C", "O"]
    assert all(ts.node(s).time == 0 for s in ts.samples())
    # per_taxon 3, and per_taxon 1 with both lineage populations (A and C keep no sample of their own)
    assert [lab for lab, _ in sample_labels(Demog(per_taxon=3))] == \
        ["A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3", "O1"]
    d = poptoy.small_demog(per_taxon=1, **LINEAGES)
    assert sample_labels(d) == [("A1", "A"), ("B1", "B"), ("C1", "C"), ("O1", "O")]
    assert pop_names(ancestry(d, Rates(), L, seed=1)) == ["A1p", "B", "C1p", "O"]


def test_no_pulse_no_tracts():
    for code in PULSES:
        extra = LINEAGES if "1" in code else {}
        d = poptoy.small_demog(f=0.0, ta=500, pulse=code, **extra)
        assert tracts(ancestry(d, Rates(), L, seed=2), d) == {}, code
    # the fail-fast stays strict at f = 0: a lineage pulse still needs its lineage population
    for code in ("C>A1", "C1>A1", "A1>C1"):
        with pytest.raises(ValueError, match="sub_"):
            poptoy.small_demog(f=0.0, ta=500, pulse=code)


def test_full_pulse_covers_recipient():
    d = poptoy.small_demog(f=1.0, ta=100, pulse="C>A")
    ts = ancestry(d, Rates(), L, seed=4)
    assert tracts(ts, d) == {0: [(0, L)], 1: [(0, L)]}  # A1, A2 whole; B, C and O none


def test_reverse_direction():
    d = poptoy.small_demog(f=1.0, ta=100, pulse="A>C")
    ts = ancestry(d, Rates(), L, seed=5)
    assert tracts(ts, d) == {4: [(0, L)], 5: [(0, L)]}  # C1, C2 whole; A none


def test_lineage_pulse():
    d = poptoy.small_demog(sub_a=(2000, 1000), sub_c=(30000, 1000), pulse="C1>A1", f=1.0, ta=500)
    ts = ancestry(d, Rates(), L, seed=6)
    assert tracts(ts, d) == {0: [(0, L)]}  # A1 (in A1p) whole; A2 (in A) none


def test_tract_through_ancestor():
    # A1 and A2 coalesce below ta in every tree, so at each position the one A lineage the pulse may move (a coin
    # flip at f = 0.5) is their ancestor, not a sample, and the tract of both is the union of the pulse records'
    # intervals. Recombination spreads that ancestry over a few lineages (ne_a 100 on 1 Mb), so a pulse can move
    # part of it: the case is the first seed from 3 whose pulse moves part of the sequence, not all of it (seed 3
    # with msprime 1.4.2); every seed scanned is checked the same way
    d = poptoy.small_demog(ne_a=100, ta=1000, f=0.5)
    length = 1_000_000
    for seed in range(3, 23):
        ts = ancestry(d, Rates(), length, seed)
        if not all(tree.tmrca(0, 1) < d.ta for tree in ts.trees()):
            continue
        rec = pulse_records(ts, d)
        union = merge(zip(ts.migrations_left[rec], ts.migrations_right[rec]))
        t = tracts(ts, d)
        assert t == ({0: union, 1: union} if union else {}), seed
        if union and union != [(0, length)]:
            break
    else:
        pytest.fail("no seed from 3 to 22 moves part of A's ancestry")
    assert all(ts.node(u).time > 0 for u in ts.migrations_node[rec])  # every move is of an ancestor
    for l, r in t[0]:
        assert 0 <= l < r <= length and isinstance(l, int) and isinstance(r, int)
    assert all(r1 < l2 for (_, r1), (l2, _) in zip(t[0], t[0][1:]))  # merged: sorted, disjoint, not adjacent


def migrant_truth(ts, d: Demog) -> dict[int, list[tuple[int, int]]]:
    """Sample index -> the merged intervals whose path to the root passes a MIGRANT node of the pulse (time ta,
    the donor population): msprime's own nodes for the move, read tree by tree without migration records."""
    import msprime

    donor = {"C>A": "C", "A>C": "A"}[d.pulse]
    pid = {ts.population(i).metadata["name"]: i for i in range(ts.num_populations)}[donor]
    mark = ((ts.nodes_flags & msprime.NODE_IS_MIG_EVENT) > 0) & (ts.nodes_time == d.ta) \
        & (ts.nodes_population == pid)
    time = ts.nodes_time
    pieces: dict[int, list[tuple[int, int]]] = {}
    for tree in ts.trees():
        parent = tree.parent_array
        for i, s in enumerate(ts.samples()):
            u = s
            while u != -1 and time[u] <= d.ta and not mark[u]:
                u = parent[u]
            if u != -1 and mark[u]:
                pieces.setdefault(i, []).append(tree.interval)
    return {i: merge(v) for i, v in sorted(pieces.items())}


def test_tracts_three_per_taxon():
    # three samples per taxon: the samples below a pulse record's node change along it, from tree to tree, so a
    # tracer must read every tree the record spans (one reading only its first tree is wrong here)
    import msprime

    d = poptoy.small_demog(f=0.3, ta=2000, per_taxon=3)
    length = 2_000_000
    samples = [msprime.SampleSet(n, population=p, ploidy=1) for p, n in (("A", 3), ("B", 3), ("C", 3), ("O", 1))]
    assert [lab for lab, _ in sample_labels(d)] == ["A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3", "O1"]
    moved = set()
    for seed in range(1, 11):
        ts = msprime.sim_ancestry(samples=samples, demography=demography(d), ploidy=2,
                                  recombination_rate=recombination_map(Rates(), length), record_migrations=True,
                                  random_seed=seed, additional_nodes=msprime.NodeType.MIGRANT,
                                  coalescing_segments_only=False)
        truth = migrant_truth(ts, d)
        assert tracts(ts, d) == truth, seed
        moved |= set(truth)
    assert moved == {0, 1, 2}  # each A sample has a tract in some seed; no other sample ever does


def test_fraction_matches_f():
    length, f = 2_000_000, 0.18
    d = poptoy.small_demog(f=f, ta=500)
    frac = [sum(r - l for l, r in tracts(ancestry(d, Rates(), length, seed), d).get(0, [])) / length
            for seed in range(1, 21)]
    se = np.std(frac, ddof=1) / np.sqrt(len(frac))
    assert se > 0
    assert abs(np.mean(frac) - f) <= 3 * se, (np.mean(frac), se)


def test_bottleneck_sizes():
    dem = demography(Demog(bottleneck=(5000, 2000)))  # A at Ne 5,000 for the 2,000 generations after t_ab 9,000
    a = dem["A"].id
    assert dem.debug().population_size_trajectory([1000, 8000])[:, a].tolist() == [30000, 5000]
    assert demography(Demog()).debug().population_size_trajectory([1000, 8000])[:, a].tolist() == [30000, 30000]


def test_recombination_map():
    rm = recombination_map(Rates(), L)  # r_low on [0.3 L, 0.7 L), r elsewhere
    assert rm.sequence_length == L
    assert rm.get_rate([0.5 * L, 0.1 * L]).tolist() == [1e-9, 1e-8]
    assert rm.get_rate([0.3 * L - 1, 0.3 * L, 0.7 * L - 1, 0.7 * L]).tolist() == [1e-8, 1e-9, 1e-9, 1e-8]
    flat = recombination_map(Rates(low=(0.0, 0.0)), L)  # no low-recombination interval
    assert flat.get_rate([0, L - 1]).tolist() == [1e-8, 1e-8]
    edge = recombination_map(Rates(low=(0.0, 0.5)), L)
    assert edge.get_rate([0, L / 2 - 1, L / 2]).tolist() == [1e-9, 1e-9, 1e-8]


def test_demog_from_json(tmp_path):
    d = Demog(ne_a=31_000, ne_b=101_000, ne_c=71_000, ne_bc=151_000, ne_o=152_000, ne_root=153_000, t_ab=8_000,
              t_bc=61_000, t_o=1_200_000, f=0.25, ta=400.0, pulse="C1>A1", sub_a=(1_500.0, 900.0),
              sub_c=(20_000.0, 800.0), bottleneck=(4_000.0, 1_000.0), per_taxon=3)
    r = Rates(mu=2e-8, r=2e-8, r_low=2e-9, low=(0.25, 0.75), lam=1e-10, delta=1 / 500_000, indel=2e-9,
              copy_noise=0.01, loss_split=(0.4, 0.3, 0.2, 0.1), replace_mean=10_000.0)
    assert all(getattr(d, f.name) != f.default for f in fields(Demog))  # every field leaves its default
    assert all(getattr(r, f.name) != f.default for f in fields(Rates))
    p = tmp_path / "demography.json"
    p.write_text(json.dumps({**asdict(d), **asdict(r)}, indent=1))
    assert demog_from_json(p) == (d, r)
    p.write_text(json.dumps({"f": 0.1, "low": [0.2, 0.6]}))  # fields left out keep their defaults
    assert demog_from_json(p) == (Demog(f=0.1), Rates(low=(0.2, 0.6)))
    p.write_text(json.dumps({"f": 0.1, "Ne_a": 5000}))
    with pytest.raises(ValueError, match="Ne_a"):
        demog_from_json(p)
    p.write_text(json.dumps([0.1]))
    with pytest.raises(ValueError, match="object"):
        demog_from_json(p)
    p.write_text(json.dumps({"f": 2.0}))  # the fields are checked as Demog checks them
    with pytest.raises(ValueError, match="f must be"):
        demog_from_json(p)


def test_flat_fields_disjoint():
    # demog_from_json reads Demog and Rates from one flat object: a field name both had would feed both
    assert flat_fields(Demog, Rates) == [{f.name for f in fields(Demog)}, {f.name for f in fields(Rates)}]

    @dataclass(frozen=True)
    class Other:
        mu: float = 0.0
        depth: int = 0

    with pytest.raises(ValueError, match="Rates and Other share field.* mu"):
        flat_fields(Demog, Rates, Other)


def test_bad_demog_fails():
    bad = [
        {"ta": 9000, "f": 0.1},                                         # ta not below t_ab with a pulse
        {"pulse": "A-C"},                                               # unknown pulse code
        {"f": 1.5}, {"f": -0.1},                                        # f outside [0, 1]
        {"pulse": "C1>A1", "sub_a": (2000, 1000), "ta": 500, "f": 0.5},  # C1>A1 without sub_c
        {"pulse": "C>A1", "ta": 500, "f": 0.5},                         # C>A1 without sub_a
        {"pulse": "C>A1", "sub_a": (2000, 1000), "ta": 2000, "f": 0.5},  # ta not below A1's lineage split
        {"pulse": "A1>C1", "sub_a": (2000, 1000), "sub_c": (400, 1000), "ta": 500, "f": 0.5},  # nor C1's
        {"per_taxon": 0},
        {"bottleneck": (5000, 9001)},                                   # longer than t_ab
    ]
    for kw in bad:
        with pytest.raises(ValueError):
            Demog(**kw)
    Demog(ta=9000, f=0.0)  # no pulse: ta is not read, against t_ab or a lineage split (sub_* is still required)
    Demog(pulse="C>A1", sub_a=(2000, 1000), f=0.0)  # the default ta, 5,000, above A1p's split
    Demog(pulse="A1>C1", sub_a=(2000, 1000), sub_c=(400, 1000), ta=500, f=0.0)  # ta above C1p's split
    Demog(bottleneck=(5000, 9000))  # as long as t_ab


def test_import_without_msprime(monkeypatch):
    # a fresh interpreter in which msprime and tskit cannot be imported: the package, the simulator and the
    # command lines still import; the population model says how to install the extra
    src = Path(sys.modules["synltr"].__file__).parents[1]
    code = """
import sys
sys.modules["msprime"] = None
sys.modules["tskit"] = None
import synltr, synltr.sim, synltr.cli, synltr.sim.cli, synltr.sim.pop.events, synltr.sim.pop.build
import synltr.sim.pop.truth
from synltr.sim.pop import model as popmodel
assert synltr.__file__.startswith(sys.argv[1]), synltr.__file__
try:
    popmodel.ancestry(popmodel.Demog(), popmodel.Rates(), 1000, 1)
except RuntimeError as e:
    print(e)
"""
    env = {**os.environ, "PYTHONPATH": os.pathsep.join(p for p in (str(src), os.environ.get("PYTHONPATH")) if p)}
    out = subprocess.run([sys.executable, "-c", code, str(src)], env=env, capture_output=True, text=True,
                         check=False)
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip() == MSPRIME_HELP
    # in this interpreter too, once msprime cannot be imported
    monkeypatch.setitem(sys.modules, "msprime", None)
    with pytest.raises(RuntimeError) as e:
        ancestry(poptoy.small_demog(), Rates(), L, seed=1)
    assert str(e.value) == MSPRIME_HELP
    with pytest.raises(RuntimeError, match=r"pip install -e '\.\[sim\]'"):
        demography(Demog())
    with pytest.raises(RuntimeError, match=r"pip install -e '\.\[sim\]'"):
        recombination_map(Rates(), L)
