"""Trio roles and tests (analyses.geneflow.trio): D_pav on the loci the donor discovers, read from its own
calls, with its strata (age bins, high tier, nested clusters collapsed, no outgroup filter), the age
gradient, the outgroup filter, resolution, roles, the f4-ratio and the block-size scan of the headline D, on toy
run directories (catalog_toy.py)."""
from __future__ import annotations

import math

import numpy as np
import pytest
from analyses.geneflow.catalog_toy import GENOMES, TRACTS, block_size, dmat, loc, simulate_run, write_run

from ltrtrace.analyses.geneflow.catalog import AGE_LABELS, TIER_CODE, AnalysisError
from ltrtrace.analyses.geneflow.relations import discover, pav_distances
from ltrtrace.analyses.geneflow.trio import NOT_TESTABLE, Trio, _f4, d_pav, second_donors, trio_tests, trios
from ltrtrace.stats.report import Findings
from ltrtrace.stats.resample import block_ids


def f4_ratio(cat, trio, d2, size, dist):
    """The f4_ratio rows of the second donor d2 ([] when no genome may discover)."""
    c = discover(cat, trio.R, trio.S, trio.D, dist, "E")
    return [] if c is None else _f4(cat, trio, d2, c, size, Findings())


def _trio_tests(cat, trio):
    return trio_tests(cat, trio, block_size(cat), pav_distances(cat)[0])


def _row(rows, test, stratum, statistic):
    return next(r for r in rows if (r.test, r.stratum, r.statistic) == (test, stratum, statistic))


def test_counted_rule_and_d(tmp_path, toy_fastas):
    loci = ([loc(10_000 * i, "R:P S:E D:P* O:E") for i in range(1, 5)]
            + [loc(10_000 * i, "R:E S:P D:P* O:E") for i in range(5, 7)]
            + [loc(70_000, "R:P S:E D:P* O:P"), loc(80_000, "R:P S:E D:P* O:U"),
               loc(90_000, "R:P* S:E D:P O:E"), loc(100_000, "R:P S:U D:P* O:E")])
    cat = write_run(tmp_path, loci, fasta_dir=toy_fastas).load()
    c = discover(cat, "R", "S", "D", pav_distances(cat)[0], "E")
    jk, n_rd, n_sd = d_pav(cat, "R", "S", "D", c, block_size(cat))
    assert (n_rd, n_sd) == (4, 2) and jk.est == pytest.approx(2 / 6) and c.cand.sum() == 7
    assert c.loci == "disc=D;present=D;res=R,S;outgroup=E"


def test_counted_reads_discoverer_calls(tmp_path, toy_fastas):
    loci = ([loc(10_000 * i, "R:P S:E D:P* O:E") for i in range(1, 5)]                               # counted
            + [loc(10_000 * i, "R:P* S:E D:P* O:E", calls={("D", "R"): "U"}) for i in range(5, 8)]  # R P/annot merged
            + [loc(10_000 * i, "R:P* S:E D:P* O:E", calls={("D", "S"): "U"}) for i in range(8, 10)] # S E merged (R's call)
            + [loc(10_000 * i, "R:E S:P D:P* O:E") for i in range(10, 14)])                         # counted
    cat = write_run(tmp_path, loci, fasta_dir=toy_fastas).load()
    c = discover(cat, "R", "S", "D", pav_distances(cat)[0], "E")
    _, n_rd, n_sd = d_pav(cat, "R", "S", "D", c, block_size(cat))
    assert (n_rd, n_sd) == (4, 4)
    assert (cat.ptier_by(c.gstar, "R")[c.rows] == TIER_CODE["high"]).all()        # call tiers; never "annot"


def test_null_with_unequal_annotation(tmp_path, toy_fastas):
    # a depth-leak case: no tracts, R annotated more completely than S
    for seed in (1, 2, 3):
        cat = simulate_run(tmp_path / str(seed), seed=seed, n_per_seq=3000, p_member={"R": 0.95, "S": 0.5},
                           p_missing=0.3, fasta_dir=toy_fastas).load()
        c = discover(cat, "R", "S", "D", pav_distances(cat)[0], "E")
        jk, _, _ = d_pav(cat, "R", "S", "D", c, block_size(cat))
        assert abs(jk.z) < 3
        if seed == 1:                                # the same rule on merged loci.tsv states: the leak it avoids
            m = (cat.member[:, cat.g("D")] & (cat.cls("D") == 1) & (cat.cls("O") == 0)
                 & (cat.cls("R") >= 0) & (cat.cls("S") >= 0))
            a = (m & (cat.cls("R") == 1) & (cat.cls("S") == 0)).sum()
            b = (m & (cat.cls("S") == 1) & (cat.cls("R") == 0)).sum()
            assert (a - b) / (a + b) > 0.05          # expected about +0.13


def test_roles_from_sign(tmp_path, toy_fastas):
    loci = ([loc(1000 * i, "R:P* S:P D:E O:E") for i in range(1, 13)]           # sisters R, S
            + [loc(20_000 + 1000 * i, "R:E S:P D:P* O:E") for i in range(4)]   # S+D x4
            + [loc(30_000, "R:P S:E D:P* O:E")])                               # R+D x1: D(R,S;D) = -0.6
    cat = write_run(tmp_path, loci, fasta_dir=toy_fastas).load()
    found, rows, _ = trios(cat, None, block_size(cat), pav_distances(cat)[0])
    assert found == [Trio("S", "R", "D", False)] and _row(rows, "roles", "all", "source").value == "data"


def test_fixed_trio(tmp_path, toy_fastas):
    cat = write_run(tmp_path, [loc(1000, "R:P* S:E D:P* O:E")], fasta_dir=toy_fastas).load()
    size, dist = block_size(cat), pav_distances(cat)[0]
    found, rows, _ = trios(cat, ("S", "R", "D"), size, dist)
    assert found == [Trio("S", "R", "D", True)] and _row(rows, "roles", "all", "source").value == "--trio"
    for bad in (("R", "R", "D"), ("R", "S", "O")):
        with pytest.raises(AnalysisError):
            trios(cat, bad, size, dist)


def test_null_and_admixed(tmp_path, toy_fastas):
    for seed in range(1, 6):
        cat = simulate_run(tmp_path / f"n{seed}", seed=seed, fasta_dir=toy_fastas).load()
        dist = pav_distances(cat)[0]
        (trio,), _, _ = trios(cat, None, block_size(cat), dist)
        jk, _, _ = d_pav(cat, trio.R, trio.S, trio.D, discover(cat, trio.R, trio.S, trio.D, dist, "E"), block_size(cat))
        assert {trio.R, trio.S} == {"R", "S"} and trio.D == "D" and abs(jk.z) < 3
    cat = simulate_run(tmp_path / "adm", seed=1, tracts=TRACTS, n_per_seq=3000, fasta_dir=toy_fastas).load()
    dist = pav_distances(cat)[0]
    (trio,), _, _ = trios(cat, None, block_size(cat), dist)
    rows, f = trio_tests(cat, trio, block_size(cat), dist)
    assert trio == Trio("R", "S", "D", False)
    assert _row(rows, "D_pav", "all", "D").value > 0 and _row(rows, "D_pav", "all", "D").z > 3
    assert _row(rows, "D_pav", "K=0", "D").value > _row(rows, "D_pav", "K>=0.01", "D").value
    # K=0 is 92/0, so the jackknife sees none of that bin's variance, and a dD.z > 0 would rest on it alone. dD
    # stands (+0.871), its se and z are nan, and the line says why
    g = _row(rows, "D_pav_gradient", "K=0 minus K>=0.01", "dD")
    assert g.value == pytest.approx(0.8712, abs=1e-4) and math.isnan(g.se) and math.isnan(g.z)
    line = next(t for s, t in f.items() if s == "Dating")
    assert NOT_TESTABLE in line and "K=0: no S+D locus" in line and "gene flow" not in line


def test_fallback_discovery(tmp_path, toy_fastas):
    genomes = [("R", "taxR", None, True), ("S", "taxS", None, True), ("D", "taxD", None, False),
               ("X", "taxD", None, True), ("O", "taxO", "outgroup", False)]
    cat = write_run(tmp_path, [loc(1000, "R:P S:E D:P X:P* O:E")], genomes=genomes, fasta_dir=toy_fastas).load()
    near_d = dmat({("R", "S"): 0.2, ("D", "X"): 0.1}, 0.5, cat.labels)
    near_r = dmat({("R", "S"): 0.2, ("R", "X"): 0.1, ("S", "X"): 0.25}, 0.5, cat.labels)
    assert discover(cat, "R", "S", "D", near_d, "E").disc == ("X",)
    assert discover(cat, "R", "S", "D", near_r, "E") is None
    _, f = trio_tests(cat, Trio("R", "S", "D", True), block_size(cat), near_r)
    assert any(t.startswith("WARNING:") and "skipped" in t and "X" in t for _, t in f.items())


def test_resolution_warning(tmp_path, toy_fastas):
    cat = simulate_run(tmp_path, seed=1, p_missing={"S": 0.6}, fasta_dir=toy_fastas).load()
    rows, f = _trio_tests(cat, Trio("R", "S", "D", True))
    assert _row(rows, "resolution", "all", "resolved_R").value == pytest.approx(0.85, abs=0.03)
    assert _row(rows, "resolution", "all", "resolved_S").value == pytest.approx(0.40, abs=0.03)
    assert any(t.startswith("WARNING:") and "resolved" in t for _, t in f.items())


def test_nested_collapse(tmp_path, toy_fastas):
    loci = []
    for h in range(10):                              # 10 RD hosts in D, each with 3 nested RD children
        host = len(loci)
        loci += [loc(50_000 * (h + 1), "R:P S:E D:P* O:E", length=6000)]
        loci += [loc(50_000 * (h + 1) + 1000 * (j + 1), "R:P S:E D:P* O:E", parent=host) for j in range(3)]
    host = len(loci)                                 # a host present in R, S and D with 3 RD children
    loci += [loc(100_000, "R:P S:P D:P* O:E", seqid="chr2", length=6000)]
    loci += [loc(101_000 + 1000 * j, "R:P S:E D:P* O:E", seqid="chr2", parent=host) for j in range(3)]
    cat = write_run(tmp_path, loci, fasta_dir=toy_fastas).load()
    rows, _ = _trio_tests(cat, Trio("R", "S", "D", True))
    assert _row(rows, "D_pav", "all", "n_RD").value == 43
    # 10 + 1; outermost rows alone would give 10
    assert _row(rows, "D_pav", "nested_collapsed", "n_RD").value == 11


def test_age_strata_without_k(tmp_path, toy_fastas):
    cat = simulate_run(tmp_path, seed=1, k_known=False, fasta_dir=toy_fastas).load()
    rows, f = _trio_tests(cat, Trio("R", "S", "D", True))
    assert not [r for r in rows if r.stratum in AGE_LABELS or r.test == "D_pav_gradient"]
    assert sum("age strata skipped" in t for _, t in f.items()) == 1


def test_outgroup_mostly_missing(tmp_path, toy_fastas):
    cat = simulate_run(tmp_path, seed=1, p_missing={"O": 0.8}, fasta_dir=toy_fastas).load()
    rows, f = _trio_tests(cat, Trio("R", "S", "D", True))
    removed = _row(rows, "outgroup_filter", "all", "removed_frac").value
    assert removed == pytest.approx(0.85, abs=0.03)
    assert any(r.test == "D_pav" and r.stratum == "no_outgroup_filter" for r in rows)
    assert any(f"{removed:.2f}" in t for _, t in f.items())


def test_second_donors_and_f4(tmp_path, toy_fastas):
    genomes = [GENOMES[0], GENOMES[1], GENOMES[2], ("D2", "taxD", None, True), GENOMES[3]]
    loci = ([loc(1000 * i, "R:P S:E D:P* D2:P O:E") for i in range(1, 4)]
            + [loc(5000, "R:E S:P D:P* D2:E O:E")]
            + [loc(6000 + 1000 * i, "R:E S:E D:P* D2:P O:E") for i in range(2)])
    cat = write_run(tmp_path, loci, genomes=genomes, fasta_dir=toy_fastas).load()
    trio, dist = Trio("R", "S", "D", True), dmat({("D", "D2"): 0.1}, 0.3, cat.labels)
    assert second_donors(cat, trio, dist) == ["D2"]
    assert second_donors(cat, trio, dmat({("D", "D2"): 0.1, ("D2", "S"): 0.05}, 0.3, cat.labels)) == []
    got = {r.statistic: r.value for r in f4_ratio(cat, trio, "D2", block_size(cat), dist)}
    assert (got["num"], got["den"], got["f"], got["resolved_D2"]) == (2, 4, 0.5, 1.0)


# ---------------------------------------------------------------- edge cases


def _blocks_jackknife(stat, blocks):
    """Delete-one-block jackknife by hand: (est, se) over the distinct blocks."""
    est = stat(np.ones(len(blocks), dtype=bool))
    p = np.array([stat(blocks != b) for b in np.unique(blocks)])
    return est, math.sqrt((len(p) - 1) / len(p) * ((p - p.mean()) ** 2).sum())


def test_gradient_is_one_joint_jackknife(tmp_path, toy_fastas):
    # 40-kb blocks (4 Mb genomes / 100); K=0 and K>=0.01 loci share blocks 0-3, so the two D's move together
    young = [(10_000, "RD"), (20_000, "SD"), (50_000, "RD"), (90_000, "RD"), (130_000, "RD")]
    old = [(30_000, "SD"), (60_000, "RD"), (70_000, "SD"), (100_000, "RD"), (140_000, "SD"), (150_000, "SD")]
    spec = {"RD": "R:P S:E D:P* O:E", "SD": "R:E S:P D:P* O:E"}
    loci = ([loc(p, spec[x], k=0.0) for p, x in young] + [loc(p, spec[x], k=0.02) for p, x in old]
            + [loc(170_000, "R:P S:E D:P* O:E", k=0.004)])                   # 0<K<0.01: in no gradient bin
    cat = write_run(tmp_path, loci, fasta_dir=toy_fastas).load()
    rows, f = _trio_tests(cat, Trio("R", "S", "D", True))
    g = _row(rows, "D_pav_gradient", "K=0 minus K>=0.01", "dD")
    a = np.array([x == "RD" for _, x in young + old])
    b = ~a
    y = np.array([True] * len(young) + [False] * len(old))
    blocks = block_ids(np.zeros(len(a), dtype=np.int64), np.array([p + 74.5 for p, _ in young + old]), 40_000)

    def d(keep, sel):
        sa, sb = (a & sel)[keep].sum(), (b & sel)[keep].sum()
        return (sa - sb) / (sa + sb)
    est, se = _blocks_jackknife(lambda keep: d(keep, y) - d(keep, ~y), blocks)
    assert g.value == pytest.approx(est) and est == pytest.approx(0.6 - (-1 / 3))
    assert g.se == pytest.approx(se) and g.z == pytest.approx(est / se) and g.n == len(a)
    # both bins hold both classes over 4 blocks: testable; a young-biased excess, Z +2.00 (se sqrt(0.2175))
    assert g.z > 0 and g.z == pytest.approx(2.0013, abs=1e-4) and se == pytest.approx(math.sqrt(0.2175))
    assert _row(rows, "D_pav", "0<K<0.01", "n_RD").value == 1
    assert any(s == "Dating" and "K=0" in t and "not significant" in t for s, t in f.items())


def test_gradient_sparse_young_bin_not_testable(tmp_path, toy_fastas):
    # one young R+D locus beside 60 old loci split 30/30 over chr2's blocks: the old bin's se alone would give
    # dD +1.000 with Z 13.81, a gene-flow reading
    rd, sd = "R:P S:E D:P* O:E", "R:E S:P D:P* O:E"
    loci = ([loc(20_000, rd, k=0.0)]
            + [loc(5_000 + 13_000 * j, rd if j % 2 == 0 else sd, k=0.02, seqid="chr2") for j in range(60)])
    cat = write_run(tmp_path, loci, fasta_dir=toy_fastas).load()
    rows, f = _trio_tests(cat, Trio("R", "S", "D", True))
    g = _row(rows, "D_pav_gradient", "K=0 minus K>=0.01", "dD")
    assert g.value == pytest.approx(1.0) and math.isnan(g.se) and math.isnan(g.z)
    young, old = _row(rows, "D_pav", "K=0", "D"), _row(rows, "D_pav", "K>=0.01", "D")
    assert young.value == 1.0 and math.isnan(young.se) and math.isnan(young.z)       # no S+D locus at K=0
    assert old.value == 0.0 and old.se > 0 and old.z == 0.0          # 30/30 over 20 blocks: testable
    line = next(t for s, t in f.items() if s == "Dating")
    assert NOT_TESTABLE in line and "K=0: no S+D locus" in line and "gene flow" not in line
    assert any(f"K=0 +1.000 ({NOT_TESTABLE}: no S+D locus; 1/0)" in t for _, t in f.items())


def test_gradient_balanced_bins(tmp_path, toy_fastas):
    # no excess at either age (3/3 in each bin over the same 3 blocks): dD 0, testable, Z 0
    young = [(5_000, "RD"), (15_000, "RD"), (45_000, "SD"), (55_000, "SD"), (85_000, "RD"), (95_000, "SD")]
    old = [(25_000, "SD"), (35_000, "SD"), (65_000, "RD"), (75_000, "RD"), (105_000, "RD"), (115_000, "SD")]
    spec = {"RD": "R:P S:E D:P* O:E", "SD": "R:E S:P D:P* O:E"}
    loci = [loc(p, spec[x], k=0.0) for p, x in young] + [loc(p, spec[x], k=0.02) for p, x in old]
    cat = write_run(tmp_path, loci, fasta_dir=toy_fastas).load()
    rows, f = _trio_tests(cat, Trio("R", "S", "D", True))
    g = _row(rows, "D_pav_gradient", "K=0 minus K>=0.01", "dD")
    # pseudovalues -1, +1, 0 (deleting blocks 0, 1, 2): se = sqrt(2/3 * 2)
    assert g.value == 0.0 and g.se == pytest.approx(math.sqrt(4 / 3)) and g.z == 0.0
    assert any(s == "Dating" and "not significant" in t and "gene flow" not in t for s, t in f.items())


def test_strata_tiers_and_loci(tmp_path, toy_fastas):
    # D's call into R is T at one locus (tier medium): counted as present, but not in the high_tier stratum
    loci = ([loc(10_000 * i, "R:P S:E D:P* O:E") for i in range(1, 4)]
            + [loc(40_000, "R:T S:E D:P* O:E"), loc(50_000, "R:E S:P D:P* O:E"),
               loc(60_000, "R:E S:P D:P* O:U")])                                # outgroup missing: filtered
    cat = write_run(tmp_path, loci, fasta_dir=toy_fastas).load()
    rows, f = _trio_tests(cat, Trio("R", "S", "D", True))
    strata = ("all", "high_tier", "no_outgroup_filter")
    assert [_row(rows, "D_pav", s, "n_RD").value for s in strata] == [4, 3, 4]
    assert _row(rows, "D_pav", "no_outgroup_filter", "n_SD").value == 2
    assert _row(rows, "D_pav", "all", "D").n == 5 and _row(rows, "D_pav", "all", "n_RD").n == 5
    assert _row(rows, "outgroup_filter", "all", "removed_frac").value == pytest.approx(1 / 6)
    assert _row(rows, "D_pav", "no_outgroup_filter", "D").loci == "disc=D;present=D;res=R,S;outgroup=none"
    assert all(r.loci not in (".", "") for r in rows)
    assert all(r.layer == ("trio" if r.test == "block_scan" else "pav") for r in rows)
    assert any("0.17" in t and "missing" in t for _, t in f.items())


def test_trios_skip_and_undefined_sign(tmp_path, toy_fastas):
    two = write_run(tmp_path / "a", [loc(1000, "R:P* S:E O:E")], genomes=[GENOMES[0], GENOMES[1], GENOMES[3]],
                    fasta_dir=toy_fastas).load()
    found, rows, f = trios(two, None, block_size(two), pav_distances(two)[0])
    assert (found, rows) == ([], [])
    assert ("Trio tests", "Trio tests skipped: fewer than three ingroup genomes") in f.items()
    # sisters R, S (merged sharing), but no locus D discovers is informative: D_pav undefined, R stays first
    cat = write_run(tmp_path / "b", [loc(1000 * i, "R:P* S:P D:E O:E") for i in range(1, 4)]
                    + [loc(9000, "R:P S:P D:P* O:E")], fasta_dir=toy_fastas).load()
    found, rows, f = trios(cat, None, block_size(cat), pav_distances(cat)[0])
    assert found == [Trio("R", "S", "D", False)]
    assert any(t.startswith("WARNING:") and "recipient is not identified" in t for _, t in f.items())
    assert _row(rows, "roles", "all", "source").value == "sheet order"              # not "data"
    # a tie (one R+D and one S+D locus D discovers): D_pav 0, sheet order decides too
    tie = write_run(tmp_path / "c", [loc(1000 * i, "R:P* S:P D:E O:E") for i in range(1, 4)]
                    + [loc(9000, "R:P S:E D:P* O:E"), loc(10_000, "R:E S:P D:P* O:E")],
                    fasta_dir=toy_fastas).load()
    found, rows, f = trios(tie, None, block_size(tie), pav_distances(tie)[0])
    assert found == [Trio("R", "S", "D", False)]
    assert _row(rows, "roles", "all", "source").value == "sheet order"
    assert any(t.startswith("WARNING:") and "a tie" in t for _, t in f.items())


def test_fallback_counts_from_the_fallback_calls(tmp_path, toy_fastas):
    # D has no annotation; X splits with D from R and S: X's own calls are read, its copy's K bins the ages
    genomes = [("R", "taxR", None, True), ("S", "taxS", None, True), ("D", "taxD", None, False),
               ("X", "taxD", None, True), ("O", "taxO", "outgroup", False)]
    loci = ([loc(1000 * i, "R:P S:E D:P X:P* O:E", k=0.0) for i in range(1, 4)]
            + [loc(5000, "R:E S:P D:P X:P* O:E", k=0.02),
               loc(6000, "R:P S:E D:U X:P* O:E", k=0.02)])                       # D missing in X's call
    cat = write_run(tmp_path, loci, genomes=genomes, fasta_dir=toy_fastas).load()
    near_d = dmat({("R", "S"): 0.2, ("D", "X"): 0.1}, 0.5, cat.labels)
    rows, f = trio_tests(cat, Trio("R", "S", "D", True), block_size(cat), near_d)
    assert (_row(rows, "D_pav", "all", "n_RD").value, _row(rows, "D_pav", "all", "n_SD").value) == (3, 1)
    assert _row(rows, "D_pav", "all", "D").loci == "disc=X;present=D;res=R,S;outgroup=E"
    assert _row(rows, "D_pav", "K=0", "n_RD").value == 3 and _row(rows, "D_pav", "K>=0.01", "n_SD").value == 1
    assert any("four-point" in t for _, t in f.items())


def test_f4_skips_and_warns(tmp_path, toy_fastas):
    # D2 qualifies by distance but shares fewer insertions with D than S does: den < 0 is not a fraction
    genomes = [GENOMES[0], GENOMES[1], GENOMES[2], ("D2", "taxD", None, True), GENOMES[3]]
    loci = ([loc(1000 * i, "R:P S:P D:P* D2:E O:E") for i in range(1, 4)]
            + [loc(5000, "R:P S:E D:P* D2:E O:E")])
    cat = write_run(tmp_path, loci, genomes=genomes, fasta_dir=toy_fastas).load()
    dist = dmat({("D", "D2"): 0.1}, 0.3, cat.labels)
    rows, f = trio_tests(cat, Trio("R", "S", "D", True), block_size(cat), dist)
    assert (_row(rows, "f4_ratio", "all", "den").value, _row(rows, "f4_ratio", "all", "num").value) == (-3, 1)
    assert any(t.startswith("WARNING:") and "f4-ratio" in t and "D2" in t for _, t in f.items())
    f4 = _row(rows, "f4_ratio", "all", "f")             # two of its four classes are empty: no se, no z
    assert f4.value == pytest.approx(-1 / 3) and math.isnan(f4.se) and math.isnan(f4.z)
    why = f"{NOT_TESTABLE}: no [S present, R E] and no [D2 present, S E] locus"
    assert any(why in t for _, t in f.items())
    rows, f = trio_tests(cat, Trio("R", "S", "D", True), block_size(cat), dmat({}, 0.3, cat.labels))
    assert not [r for r in rows if r.test == "f4_ratio"]
    assert any("f4-ratio skipped" in t for _, t in f.items())
    dist[cat.g("D2"), cat.g("S")] = dist[cat.g("S"), cat.g("D2")] = np.nan       # unresolved: never qualifies
    assert second_donors(cat, Trio("R", "S", "D", True), dist) == []


def test_d_pav_keep_must_be_a_mask(tmp_path, toy_fastas):
    cat = write_run(tmp_path, [loc(1000, "R:P S:E D:P* O:E"), loc(2000, "R:E S:P D:P* O:E")],
                    fasta_dir=toy_fastas).load()
    c = discover(cat, "R", "S", "D", pav_distances(cat)[0], "E")
    jk, n_rd, n_sd = d_pav(cat, "R", "S", "D", c, block_size(cat), keep=np.array([True, False]))
    assert (n_rd, n_sd, jk.est) == (1, 0, 1.0)
    with pytest.raises(ValueError, match="mask"):
        d_pav(cat, "R", "S", "D", c, block_size(cat), keep=np.array([0, 1]))


def test_f4_second_donor_that_discovers(tmp_path, toy_fastas):
    # the same true insertions (R S D X), with D annotated or not. With D
    # unannotated, X discovers by fallback and is also the second donor, present at every counted locus by
    # construction: f would be 2/7 instead of 2/5, so the f4-ratio is skipped with a line
    truth = ([("P", "E", "P", "P")] * 3 + [("P", "E", "P", "E")] * 2 + [("E", "P", "P", "E")] * 2
             + [("E", "P", "P", "P")] + [("E", "E", "P", "P")] * 4)
    got = {}
    for name, annotated in (("D_annotated", True), ("D_fallback", False)):
        genomes = [("R", "taxR", None, False), ("S", "taxS", None, False), ("D", "taxD", None, annotated),
                   ("X", "taxD", None, True), ("O", "taxO", "outgroup", False)]
        star = {True: "*", False: ""}
        loci = [loc(10_000 * (j + 1), f"R:{r} S:{s} D:{d}{star[annotated]} X:{x}{star[x == 'P']} O:E")
                for j, (r, s, d, x) in enumerate(truth) if annotated or x == "P"]   # X finds what it holds
        cat = write_run(tmp_path / name, loci, genomes=genomes, fasta_dir=toy_fastas).load()
        dist = dmat({("R", "S"): 0.2, ("D", "X"): 0.1}, 0.5, cat.labels)
        assert second_donors(cat, Trio("R", "S", "D", True), dist) == ["X"]
        rows, f = trio_tests(cat, Trio("R", "S", "D", True), block_size(cat), dist)
        got[name] = ({r.statistic: r.value for r in rows if r.test == "f4_ratio"}, f.items(),
                     discover(cat, "R", "S", "D", dist, "E").disc)
    rows, _, disc = got["D_annotated"]
    assert disc == ("D",) and (rows["num"], rows["den"], rows["f"]) == (2, 5, pytest.approx(0.4))
    rows, lines, disc = got["D_fallback"]
    assert disc == ("X",) and rows == {}
    assert any("second donor X skipped: X discovers the counted loci" in t for _, t in lines)
    assert f4_ratio(cat, Trio("R", "S", "D", True), "X", block_size(cat), dist) == []


def test_blocks_on_d_sites(tmp_path, toy_fastas):
    # D's chr1 sites mirrored over [1, 420_000] (x -> 420_001 - x), so its 40-kb blocks group the
    # loci unlike R's, S's or X's (blocks 10, 9, 9, 8, 8, 7 against 0, 0, 1, 1, 2, 2). The jackknife must use
    # D's site as the discoverer reads it: D's own element, or X's call into D under fallback
    pats = [(5_000, "RD"), (25_000, "RD"), (45_000, "SD"), (65_000, "RD"), (85_000, "SD"), (105_000, "RD")]
    a = np.array([x == "RD" for _, x in pats])

    def d(keep):
        return (a[keep].sum() - (~a)[keep].sum()) / keep.sum()
    hand = {name: _blocks_jackknife(d, mids // 40_000)
            for name, mids in (("D", np.array([420_001 - p - 74.5 for p, _ in pats])),     # D's mirrored mids
                               ("R", np.array([p + 74.5 for p, _ in pats])))}              # anyone else's
    assert hand["D"][1] == pytest.approx(math.sqrt(0.0675)) and hand["R"][1] == pytest.approx(1 / 3)
    spec = {"RD": "R:P S:E D:P* O:E", "SD": "R:E S:P D:P* O:E"}
    cat = write_run(tmp_path / "own", [loc(p, spec[x]) for p, x in pats], inversion=("D", "chr1", 1, 420_000),
                    fasta_dir=toy_fastas).load()
    c = discover(cat, "R", "S", "D", pav_distances(cat)[0], "E")
    jk, n_rd, n_sd = d_pav(cat, "R", "S", "D", c, block_size(cat))
    assert (n_rd, n_sd, jk.nblocks) == (4, 2, 4) and jk.est == pytest.approx(1 / 3)
    assert jk.se == pytest.approx(hand["D"][1])
    genomes = [("R", "taxR", None, True), ("S", "taxS", None, True), ("D", "taxD", None, False),
               ("X", "taxD", None, True), ("O", "taxO", "outgroup", False)]
    spec = {"RD": "R:P S:E D:P X:P* O:E", "SD": "R:E S:P D:P X:P* O:E"}
    cat = write_run(tmp_path / "fallback", [loc(p, spec[x]) for p, x in pats], genomes=genomes,
                    inversion=("D", "chr1", 1, 420_000), fasta_dir=toy_fastas).load()
    c = discover(cat, "R", "S", "D", dmat({("R", "S"): 0.2, ("D", "X"): 0.1}, 0.5, cat.labels), "E")
    jk, n_rd, n_sd = d_pav(cat, "R", "S", "D", c, block_size(cat))
    assert c.disc == ("X",) and (n_rd, n_sd, jk.nblocks) == (4, 2, 4) and jk.se == pytest.approx(hand["D"][1])


def _uniform(n_blocks, n_rd, n_sd, k=0.02, offset=0):
    """n_blocks 40-kb blocks on chr1, each holding n_rd R+D then n_sd S+D loci, 2 kb apart: every pseudovalue
    is the same ratio."""
    rd, sd = "R:P S:E D:P* O:E", "R:E S:P D:P* O:E"
    return [loc(40_000 * blk + 1_000 + offset + 2_000 * j, rd if j < n_rd else sd, k=k)
            for blk in range(n_blocks) for j in range(n_rd + n_sd)]


def test_identical_pseudovalues_are_not_testable(tmp_path, toy_fastas):
    # 3 blocks x (3 R+D, 2 S+D) give the pseudovalue 0.2 three times, whose float mean is
    # 0.20000000000000004: a se of 3.9e-17 would give Z 5.1e15, read as significant
    cat = write_run(tmp_path / "d", _uniform(3, 3, 2), fasta_dir=toy_fastas).load()
    rows, f = _trio_tests(cat, Trio("R", "S", "D", True))
    d = _row(rows, "D_pav", "all", "D")
    assert d.value == pytest.approx(0.2) and math.isnan(d.se) and math.isnan(d.z)
    line = next(t for _, t in f.items() if t.startswith("D_pav(R,S;D) ="))
    assert f"{NOT_TESTABLE}: the block pseudovalues do not vary" in line and "significant" not in line
    # the same uniform blocks in both age bins (3/2 young, 2/3 old): dD 0.4 in every pseudovalue
    cat = write_run(tmp_path / "g", _uniform(3, 3, 2, k=0.0) + _uniform(3, 2, 3, offset=20_000),
                    fasta_dir=toy_fastas).load()
    rows, f = _trio_tests(cat, Trio("R", "S", "D", True))
    g = _row(rows, "D_pav_gradient", "K=0 minus K>=0.01", "dD")
    assert g.value == pytest.approx(0.4) and math.isnan(g.se) and math.isnan(g.z)
    line = next(t for s, t in f.items() if s == "Dating")
    assert f"{NOT_TESTABLE}: the block pseudovalues do not vary" in line and "gene flow" not in line


def gradient_toy_loci():
    """A gradient toy: on each of 4 sequences, 10 repeats of 4 young S+D, 1 young
    R+D, 4 old R+D, 4 old S+D and 3 R+S loci (2.5 kb apart). D_pav(R,S;D) = (200 - 320) / 520 < 0: S carries the
    excess, all of it young (K=0: 40 R+D, 160 S+D; K>=0.01: 160 and 160)."""
    loci = []
    for seqid in ("chr1", "chr2", "chr3", "chr4"):
        pos = 10_000
        for _ in range(10):
            for n, spec, k in ((4, "R:E S:P D:P* O:E", 0.0), (1, "R:P S:E D:P* O:E", 0.0),
                               (4, "R:P S:E D:P* O:E", 0.02), (4, "R:E S:P D:P* O:E", 0.02),
                               (3, "R:P* S:P D:E O:E", 0.02)):
                for _ in range(n):
                    loci.append(loc(pos, spec, seqid=seqid, k=k))
                    pos += 2_500
    return loci


@pytest.mark.parametrize("roles, excess", [(("R", "S", "D"), -1), (("S", "R", "D"), 1)])
def test_gradient_reads_the_excess_direction(tmp_path, toy_fastas, roles, excess):
    """Under --trio roles where S carries the donor-shared excess (D_pav < 0), that excess is
    young, so dD < 0; the Dating verdict reads dD in the excess's direction and names its owner, the genome
    labeled S in both role orders, never "older insertions"."""
    cat = write_run(tmp_path, gradient_toy_loci(), fasta_dir=toy_fastas).load()
    rows, f = _trio_tests(cat, Trio(*roles, True))
    d, g = _row(rows, "D_pav", "all", "D"), _row(rows, "D_pav_gradient", "K=0 minus K>=0.01", "dD")
    assert d.value == pytest.approx(excess * 120 / 520) and abs(d.z) >= 3
    assert g.value == pytest.approx(excess * 0.6) and g.z * excess >= 3            # testable, toward the excess
    line = next(t for s, t in f.items() if s == "Dating" and t.startswith("Age gradient"))
    assert ("read in the direction of the excess, S's (S shares more insertions with D than R does: D_pav "
            f"{excess * 120 / 520:+.3f}): that excess is concentrated in young insertions, consistent with gene flow "
            "from D-like donors into S rather than lineage sorting alone") in line
    assert "older insertions" not in line


# ---------------------------------------------------------------- the block-size scan of the headline D


def test_block_scan_rows(tmp_path, toy_fastas):
    """The headline D_pav jackknifed again with 1, 2, 5, 10 and 20 Mb blocks: one row per size (layer trio, test
    block_scan, stratum "<size> Mb", statistic D), the estimate, n and loci the headline's, the se, z and block
    count each size's own. Here every 1-Mb toy sequence is one block at every size: 4 blocks, too few to warn."""
    cat = simulate_run(tmp_path, seed=1, tracts=TRACTS, n_per_seq=3000, fasta_dir=toy_fastas).load()
    rows, f = _trio_tests(cat, Trio("R", "S", "D", True))
    d = _row(rows, "D_pav", "all", "D")
    scan = [r for r in rows if r.test == "block_scan"]
    assert [r.stratum for r in scan] == ["1 Mb", "2 Mb", "5 Mb", "10 Mb", "20 Mb"]
    assert all((r.layer, r.scope, r.statistic, r.value, r.n, r.loci) == ("trio", "R,S;D", "D", d.value, d.n, d.loci)
               for r in scan)
    assert [r.note.split(" on ")[0] for r in scan] == [f"4 jackknife blocks of {s:,} bp" for s in
                                                       (1_000_000, 2_000_000, 5_000_000, 10_000_000, 20_000_000)]
    assert all(r.note.endswith("D's sites, 4 holding informative loci") for r in scan)
    assert all(math.isfinite(r.se) and r.z == pytest.approx(r.value / r.se) for r in scan)
    assert d.note.endswith("100 jackknife blocks of 40,000 bp on D's sites") and abs(d.z) >= 3
    assert not any("block-size scan" in t for _, t in f.items())


def _scan_toy(n_inf: int, n_uninf: int):
    """n_inf + n_uninf 30-kb sequences (1.05 Mb in all with 35: the rule's blocks are 10.5 kb). The first n_inf are
    split by the rule's blocks into three that hold the same 6 informative loci: (R+D, S+D) per block (5, 1),
    (1, 5), (4, 2), (4, 2), (4, 2) on consecutive sequences. The donor-shared excess varies between sequences and not
    inside them, so the rule's blocks see a third of its variance at a time (Z >= 3), while 1-20 Mb blocks, one per
    sequence, see all of it (Z < 3). The other n_uninf hold 3 counted loci each with R and S both present: a
    jackknife block each, none informative."""
    cycle = ((5, 1), (1, 5), (4, 2), (4, 2), (4, 2))
    rd, sd, both = "R:P S:E D:P* O:E", "R:E S:P D:P* O:E", "R:P S:P D:P* O:E"
    seqids = tuple(f"s{i:02d}" for i in range(n_inf + n_uninf))
    loci = []
    for i, seqid in enumerate(seqids[:n_inf]):
        a, b = cycle[i % len(cycle)]
        for base in (1_000, 11_000, 22_000):
            loci += [loc(base + 1_000 * j, rd if j < a else sd, seqid=seqid) for j in range(a + b)]
    for seqid in seqids[n_inf:]:
        loci += [loc(1_000 * j, both, seqid=seqid) for j in (1, 2, 3)]
    return loci, seqids


@pytest.mark.parametrize("n_inf, n_uninf, warned", [(35, 0, True), (20, 15, False)])
def test_block_scan_warning(tmp_path, toy_fastas, n_inf, n_uninf, warned):
    """The verdict at |Z| >= 3 flips at 20 Mb (and every scan size), where both toys have 35 jackknife blocks. With
    informative loci in 35 of them, findings.txt warns once in Trio tests, naming the block-size scan; with
    informative loci in 20, fewer than BLOCK_SCAN_MIN_BLOCKS (the other 15 hold uninformative loci only), it does
    not."""
    loci, seqids = _scan_toy(n_inf, n_uninf)
    cat = write_run(tmp_path, loci, seqids=seqids, seq_len=30_000, fasta_dir=toy_fastas).load()
    rows, f = _trio_tests(cat, Trio("R", "S", "D", True))
    d = _row(rows, "D_pav", "all", "D")
    scan = {r.stratum: r for r in rows if r.test == "block_scan"}
    assert block_size(cat) == 10_500
    assert d.note.endswith(f"{3 * n_inf + n_uninf} jackknife blocks of 10,500 bp on D's sites")
    assert abs(d.z) >= 3 and all(abs(r.z) < 3 for r in scan.values())
    assert scan["20 Mb"].note == f"35 jackknife blocks of 20,000,000 bp on D's sites, {n_inf} holding informative loci"
    lines = [(s, t) for s, t in f.items() if "block-size scan" in t]
    assert bool(lines) == warned and all(s == "Trio tests" and t.startswith("WARNING:") for s, t in lines)
    if warned:
        assert len(lines) == 1 and "20 Mb" in lines[0][1] and "informative loci in 35 blocks" in lines[0][1]
