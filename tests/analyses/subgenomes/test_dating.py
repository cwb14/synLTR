"""Dating: split estimators on copy pairs, the inherited-insertion split, the reported merger and dating.tsv."""
import math
from types import SimpleNamespace

import numpy as np
from analyses.subgenomes import sgtools as sg

from ltrtrace.analyses.subgenomes import dating, inputs


def test_merger_report_rule():
    dated = dict(datable=True, bound=False, estimate=1e-2, lo=5e-3, hi=2e-2, landing_bound=None)
    young = dict(datable=True, bound=False, estimate=1e-4, lo=5e-5, hi=3e-4, landing_bound=None)
    bound = dict(datable=True, bound=True, estimate=1e-4, lo=0.0, hi=5e-4, landing_bound=None)
    nd = dict(datable=False, bound=False, estimate=1e-3, lo=0.0, hi=1.0, landing_bound=3e-3)
    t = 1e-3                                                  # one substitution per LTR
    assert dating.merger_report(nd, dated, t)["source"] == "genealogy_clock"
    r = dating.merger_report(dated, bound, t)                     # genealogy only bounded: the LTR clock's date
    assert r["source"] == "ltr_clock" and not r["bound"]
    r = dating.merger_report(nd, young, t)                        # lower limit below the LTR resolution: a bound
    assert r["bound"] and r["lo"] == 0.0 and r["hi"] == 3e-4
    r = dating.merger_report(nd, bound, t)
    assert r["source"] == "genealogy_clock" and r["bound"]
    loose = dict(bound, hi=5e-3)
    r = dating.merger_report(loose, bound, t)                     # both bounds: the tighter one
    assert r["source"] == "genealogy_clock" and r["hi"] == 5e-4
    r = dating.merger_report(bound, dict(bound, hi=9e-4), t)
    assert r["source"] == "ltr_clock" and r["hi"] == 5e-4
    r = dating.merger_report(nd, nd, t)
    assert r["source"] == "ltr_clock" and not r["datable"] and r["landing_bound"] == 3e-3
    assert dating.merger_report(None, None, t) is None


def test_reported_row_matches_reported_merger(tmp_path):
    """dating.tsv's reported merger row carries the reported (possibly demoted) interval, as summary.json does."""
    raw = dict(estimate=5e-4, lo=2e-4, hi=1.2e-3, datable=True, bound=False, n=900, n_cross=60, landing_bound=None)
    rep = dating.merger_report(None, raw, 1e-3)                    # lower limit below 1/LTR: reported as a bound
    assert rep["bound"]
    dt = dict(mu=1.3e-8, mu_given=True, n_copy_pairs=1000, pairs=[dict(
        groups=(0, 1), split={}, split_reported=None, merger=rep, merger_ltr=None, merger_genealogy=raw, n_pairs=900)])
    dating.write_dating(SimpleNamespace(outdir=str(tmp_path)), dt)
    R = [l.split("\t") for l in (tmp_path / "dating.tsv").read_text().splitlines() if not l.startswith("#")]
    rows = [dict(zip(R[0], r)) for r in R[1:]]
    m = [r for r in rows if r["event"] == "merger" and r["reported"] == "yes"][0]
    s = dating.dating_summary(dt)["pairs"][0]["merger"]
    assert m["status"] == s["status"] == "bound" and m["k2p_low"] == "0" and float(m["k2p_high"]) == 1.2e-3


def test_half_rise_within_a_bin_of_the_rise():
    """The lineage-floor estimator (not interpolated: interpolation keeps a ~15% young bias and narrows intervals):
    a bin centre within one 0.1-log10 bin of a sharp rise at 0.03."""
    rng = np.random.default_rng(7)
    d = 10 ** rng.uniform(-3.5, -0.3, 40000)
    cross = (rng.random(d.size) < np.where(d < 0.03, 0.03, 0.5)).astype(float)
    assert abs(math.log10(dating.half_rise(d, cross) / 0.03)) <= 0.1 + 1e-9


def _planted_profile(rng, Dm=1e-3, Ds=0.04, n=40000):
    """Pairs over log-uniform divergence: cross fraction 0.5 below Dm, 0.03 on the confined plateau, 0.45 above Ds."""
    d = 10 ** rng.uniform(-4, -0.5, n)
    p = np.where(d < Dm, 0.5, np.where(d < Ds, 0.03, 0.45))
    return d, (rng.random(n) < p).astype(float)


def test_profile_tools_on_planted_profile():
    rng = np.random.default_rng(6)
    d, cross = _planted_profile(rng)
    assert 0.5e-3 <= dating.young_cut(d, cross) <= 2e-3
    assert 0.025 <= dating.plateau_end(d, cross) <= 0.06
    assert 0.03 <= dating.half_rise(d, cross) <= 0.06
    assert np.isnan(dating.half_rise(d, np.where(d > 0.04, 0.03, cross)))       # no rise: no lineage-floor split


def test_duplicated_insertions_and_ortholog_split():
    """Orthologous pairs (divergence 0.05, far below both copies' own 0.15) give the split; pairs of independent
    insertions (divergence above the older copy's own) are not duplicated insertions."""
    rng = np.random.default_rng(7)
    lib = inputs.Library()
    lib.n = 200
    lib.k2p = np.full(200, 0.15)
    lib.k2p[100:] = 0.01
    lib.k2p_se = np.full(200, np.nan)
    lib.n_sites = np.full(200, 1000.0)
    lib.ltr5_len = np.full(200, 1000.0)
    P = np.stack([np.arange(0, 100, 2), np.arange(1, 100, 2)], 1)              # old copies: orthologues
    Q = np.stack([np.arange(100, 200, 2), np.arange(101, 200, 2)], 1)          # young copies, older MRCA
    PP = np.concatenate([P, Q])
    d = np.concatenate([rng.normal(0.05, 0.003, 50), np.full(50, 0.05)])
    se = np.full(100, 0.003)
    dup = dating.duplicated_insertions(lib, PP, d, se)
    assert dup[:50].all() and not dup[50:].any()
    r = dating.split_orthologs(d, dup, np.arange(100), rng)
    assert abs(r["estimate"] - 0.05) < 0.003 and r["lo"] < 0.05 < r["hi"] and r["n"] == 50
    assert dating.split_orthologs(d[:10], dup[:10], np.arange(10), rng) is None


def test_split_confinement_planted():
    """Share of lineage evidence pointing to the other subgenome: low for copies younger than 0.04 (confined
    lineages), high for older ones (lineages present in both progenitors): the step is found at ~0.04."""
    rng = np.random.default_rng(11)
    k = 10 ** rng.uniform(-3, -0.7, 5000)
    share = np.clip(np.where(k < 0.04, 0.1, 0.5) + rng.normal(0, 0.1, k.size), 0, 1)
    r = dating.split_confinement(k, np.full(k.size, 0.003), share, np.ones(k.size), rng)
    assert r is not None and 0.03 <= r["estimate"] <= 0.05 and r["lo"] <= 0.04 <= r["hi"]
    assert dating.split_confinement(k, np.full(k.size, 0.003), np.full(k.size, 0.1), np.ones(k.size), rng) is None


def test_sample_pairs_and_date_pair_degenerate():
    """Too few phased copies: no pair table; a subgenome pair with too few copy pairs: no estimates, no crash."""
    rng = np.random.default_rng(10)
    lib = sg.lib_from_seqs(sg.random_family(rng, 10, 600), [0] * 5 + [1] * 5)
    lib.k2p = np.full(10, 0.01)
    assert dating.sample_pairs(lib, lib.seq.astype(np.int64), rng, 1) is None
    pt = dict(P=np.array([[0, 5]]), d=np.array([0.01]), se=np.array([0.001]), stratum=np.zeros(1, np.int8))
    r = dating.date_pair(lib, pt, np.array([0] * 5 + [1] * 5), 0, 1, None, rng)
    assert r["split"] == {} and r["merger"] is None and r["n_pairs"] == 1


def test_dating_texts_say_what_was_measured():
    """Log wording: a bound is a bound on the onset of cross-subgenome transposition (not on the merger); 'not
    datable' says why (too few sibling pairs, or no excess of cross-subgenome siblings); a single split estimator is
    not reported as a disagreement."""
    mu = 1e-8
    bound = dict(estimate=1e-4, lo=0.0, hi=2e-4, datable=True, bound=True, n=500, n_cross=40)
    t = dating.merger_text(bound, 12, mu)
    assert "transposition began <" in t and "not on the merger" in t and "10.0 kyr" in t
    assert "too few sibling copy pairs (12)" in dating.merger_text(None, 12, mu)
    assert "no detectable excess" in dating.merger_text(dict(bound, datable=False), 500, mu)
    one = dict(groups=(0, 1), split={"orthologs": None, "confinement": dict(estimate=0.04, lo=0.03, hi=0.05, n=9),
                                     "profile": None}, split_reported="confinement", split_agreement="single")
    t = dating.split_text(one, mu)
    assert "one estimator" in t and "disagree" not in t


def test_write_dating_bound_row_has_no_point_estimate(tmp_path):
    a = sg.parse_args(["--ltr_fasta", "x", "--annotation", "y", "--outdir", str(tmp_path)])
    mg = dict(estimate=1e-4, lo=0.0, hi=2e-4, datable=True, bound=True, n=500, n_cross=40)
    dt = dict(mu=1e-8, mu_given=True, n_copy_pairs=100, pairs=[dict(
        groups=(0, 1), split={"orthologs": None, "confinement": None, "profile": None}, split_reported=None,
        split_agreement="none", n_sibling_pairs=500, merger_ltr=None, merger_genealogy=mg,
        merger=dict(mg, source="genealogy_clock"))])
    dating.write_dating(a, dt)
    rows = [r for r in sg.dating_rows(tmp_path) if r["event"] == "merger"]
    assert [r["estimator"] for r in rows] == ["ltr_clock", "genealogy_clock"]
    row = [r for r in rows if r["reported"] == "yes"][0]
    assert row["estimator"] == "genealogy_clock"
    assert row["status"] == "bound" and row["k2p"] == "NA" and row["years"] == "NA"
    assert row["years_low"] == "0" and float(row["years_high"]) == 1e4


def test_landing_bound_reported_beside_the_merger(tmp_path):
    """dating.tsv carries the random-landing bound as its own row; the merger rows keep their status; texts say what
    the bound assumes."""
    mg = dict(estimate=0.0, lo=0.0, hi=0.0, phi=0.03, eps=0.03, datable=False, bound=False, n=500, n_cross=15,
              landing_bound=0.0013)
    t = dating.merger_text(mg, 500, 1.3e-8)
    assert "not datable" in t and "at random" in t and "50.0 kyr" in t and "own subgenome" in t
    dt = dict(mu=1.3e-8, mu_given=True, n_copy_pairs=1000,
              pairs=[dict(groups=(0, 1), split={}, merger=dict(mg, source="ltr_clock"), merger_ltr=mg,
                          merger_genealogy=None, n_pairs=900, split_reported=None)])
    dating.write_dating(SimpleNamespace(outdir=str(tmp_path)), dt)
    R = [l.split("\t") for l in (tmp_path / "dating.tsv").read_text().splitlines() if not l.startswith("#")]
    rows = [dict(zip(R[0], r)) for r in R[1:]]
    m = [r for r in rows if r["event"] == "merger" and r["reported"] == "yes"][0]
    lb = [r for r in rows if r["event"] == "merger_if_random"][0]
    assert m["estimator"] == "ltr_clock" and m["status"] == "not_datable"
    assert lb["status"] == "bound" and abs(float(lb["years_high"]) - 50000) < 1
    s = dating.dating_summary(dt)
    assert s["pairs"][0]["merger"]["status"] == "not_datable"
    assert abs(s["pairs"][0]["merger"]["if_random_landing_below_years"] - 50000) < 1


def test_inherited_estimate_recovers_median_and_drops_exchanges():
    rng = np.random.default_rng(3)
    n = 200
    k2p = rng.normal(0.08, 0.008, n)
    k2p[:12] = 0.002                                       # homogenized (exchanged) pairs
    sites = np.full(n, 800)
    blocks = np.repeat(np.arange(10), 20)
    r = dating.inherited_estimate(k2p, sites, blocks, np.random.default_rng(4))
    assert r is not None
    assert abs(r["estimate"] - 0.08) < 0.004
    assert r["n_flagged"] >= 10
    assert r["lo"] < r["estimate"] < r["hi"]


def test_inherited_estimate_needs_enough_anchors():
    k2p = np.full(10, 0.05)
    assert dating.inherited_estimate(k2p, np.full(10, 500), np.zeros(10), np.random.default_rng(0)) is None


def test_inherited_estimate_ignores_short_alignments():
    k2p = np.r_[np.full(30, 0.05), np.full(30, 0.5)]
    sites = np.r_[np.full(30, 600), np.full(30, 150)]       # the 0.5 values come from < 200 aligned sites
    r = dating.inherited_estimate(k2p, sites, np.repeat([0, 1, 2], 20), np.random.default_rng(0))
    assert r is not None and abs(r["estimate"] - 0.05) < 1e-9


def test_inherited_estimate_drops_post_merger_anchors():
    """Old split: true inherited insertions gone, exchange duplicates (young, collinear) remain -> no estimate."""
    rng = np.random.default_rng(5)
    k2p = np.r_[rng.normal(0.004, 0.0005, 100), rng.normal(0.30, 0.03, 8)]    # exchanges + a few survivors
    sites = np.full(k2p.size, 700)
    blocks = np.r_[np.repeat([0, 1], 50), np.arange(8) + 2]
    assert dating.inherited_estimate(k2p, sites, blocks, np.random.default_rng(1), floor=0.008) is None
    r = dating.inherited_estimate(k2p, sites, blocks, np.random.default_rng(1), floor=0.0)
    assert r is not None and r["estimate"] < 0.01                      # without the floor it would date the exchanges


def test_inherited_estimate_keeps_duplicated_insertions_only():
    """Orthologs: pair divergence below both copies' ages; paralogs: above -> excluded."""
    rng = np.random.default_rng(6)
    n_o, n_p = 40, 40
    d = np.r_[rng.normal(0.10, 0.008, n_o), rng.normal(0.18, 0.02, n_p)]
    older = np.r_[rng.uniform(0.14, 0.30, n_o), rng.uniform(0.11, 0.16, n_p)]   # paralogs: copies younger than d
    se = np.full(d.size, 0.005)
    blocks = np.repeat(np.arange(8), 10)
    r = dating.inherited_estimate(d, np.full(d.size, 800), blocks, np.random.default_rng(2), se=se, older=older)
    assert r is not None and abs(r["estimate"] - 0.10) < 0.01 and r["n"] <= n_o
