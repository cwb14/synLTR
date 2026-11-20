"""Painting (analyses.geneflow.paint): the classes from the recipient's own calls, the forward-backward against
brute-force enumeration, the tract rule, target-decoy q-values, tract recovery on a simulated admixed run
against a null run, the output files, the skips and the youngest insertion per tract, on toy run directories
(catalog_toy.py)."""
from __future__ import annotations

import itertools
import json
import math

import numpy as np
import pytest
from analyses.geneflow.catalog_toy import GENOMES, TRACTS, loc, read_tsv, simulate_run, write_run

from ltrtrace.analyses.geneflow.inversions import inverted_blocks
from ltrtrace.analyses.geneflow.paint import (
    PAINT_CLASSES,
    POSTERIOR_COLUMNS,
    TRACT_COLUMNS,
    HmmFit,
    _floor_line,
    _n_loci,
    call_tracts,
    classes,
    fit_hmm,
    forward_backward,
    paint,
    paint_rows,
    shortest_callable,
    target_decoy,
    write_paintings,
)
from ltrtrace.analyses.geneflow.route import ROUTE_PROPS
from ltrtrace.analyses.geneflow.trio import Trio
from ltrtrace.seq.divergence import garwood_u95 as u95

PROPS = ROUTE_PROPS["assembly"]

TOTAL = sum(b - a + 1 for _, a, b in TRACTS)

def enumerate_paths(pos, emit, pi, lam, e0, e1):
    """Brute force over the 2**n state paths of one chain: posterior of state 1 and log-likelihood."""
    n, tot, p1 = len(pos), 0.0, np.zeros(len(pos))
    for path in itertools.product((0, 1), repeat=n):
        w = pi if path[0] else 1 - pi
        for t in range(1, n):
            q = math.exp(-lam * (pos[t] - pos[t - 1])); p01, p10 = pi * (1 - q), (1 - pi) * (1 - q)
            w *= (p01 if path[t] else 1 - p01) if path[t - 1] == 0 else (1 - p10 if path[t] else p10)
        for t, k in enumerate(emit):
            if k >= 0:
                w *= (e1 if path[t] else e0)[k]
        tot += w; p1 += w * np.array(path)
    return p1 / tot, math.log(tot)

def covered_bp(tracts, cat, truth):
    """bp of the truth intervals (R seqid name, a, b) covered by called tracts on the same R sequence."""
    return sum(max(0, min(b, t.end) - max(a, t.start) + 1) for name, a, b in truth for t in tracts
               if cat.seqids["R"][t.seq] == name)

def test_forward_backward_matches_enumeration():
    pos = np.array([0.0, 5e4, 6e4, 3e5]); emit = np.array([0, 1, -1, 0])
    e0, e1 = np.array([0.3, 0.7]), np.array([0.8, 0.2])
    post, ll, _, _ = forward_backward(pos, emit, np.array([0, 4]), 0.2, 1 / 1e5, e0, e1)
    exp_post, exp_ll = enumerate_paths(pos, emit, 0.2, 1 / 1e5, e0, e1)
    assert np.allclose(post, exp_post, atol=1e-9) and ll == pytest.approx(exp_ll, abs=1e-9)

def test_tract_rule():
    tr = call_tracts(np.zeros(6, dtype=int), np.arange(6) * 1000 + 1, np.arange(6) * 1000 + 500,
                     np.array([0.9, 0.95, 0.4, 0.8, 0.99, 0.6]), pi=0.5)     # prior odds 1: post >= 0.75
    assert [(t.first, t.last) for t in tr] == [(0, 1), (3, 4)]

def test_target_decoy_plus_one():
    q, p = target_decoy(np.arange(30, 0, -1.0), np.array([5.5]), 1.0)
    assert np.allclose(q[:25], 1 / 25) and np.allclose(q[25:], 2 / 30)
    assert p[0] == pytest.approx(0.5) and p[-1] == pytest.approx(1.0)
    assert target_decoy(np.array([10.0]), np.array([1.0, 2.0]), 1.0)[0][0] == 1.0   # a lone target is not significant

def test_classes(tmp_path, toy_fastas):
    loci = [loc(1000, "R:P* S:E D:P O:E"), loc(2000, "R:P* S:P D:E O:E"), loc(3000, "R:P* S:P D:P O:E"),
            loc(4000, "R:P* S:E D:E O:E"), loc(5000, "R:P* S:U D:P O:E"), loc(6000, "R:P* S:E D:P O:P"),
            loc(7000, "R:P* S:E D:P* O:E", calls={("R", "S"): "U"})]          # merged S E (D's call); R's own call U
    cat = write_run(tmp_path, loci, fasta_dir=toy_fastas).load()
    idx, cls = classes(cat, "R", "S", "D", alphabet="insertion")
    assert idx.tolist() == list(range(7))
    assert [PAINT_CLASSES[c] for c in cls] == ["RD", "RS", "both", "neither", "missing", "old", "missing"]

def test_paint_finds_tracts(tmp_path, toy_fastas):
    cat = simulate_run(tmp_path, seed=1, tracts=TRACTS, fasta_dir=toy_fastas).load()
    p = paint(cat, Trio("R", "S", "D", False), PROPS)
    assert p.donor_state and covered_bp(p.tracts, cat, TRACTS) >= 0.8 * TOTAL and len(p.decoy_tracts) < len(p.tracts)

def test_paint_null_and_positive_control(tmp_path, toy_fastas):
    null = paint(simulate_run(tmp_path / "n", seed=1, fasta_dir=toy_fastas).load(), Trio("R", "S", "D", False), PROPS)
    assert sum(q is not None and q < 0.05 for q in null.q) == 0
    assert null.post.mean() < null.decoy_post.mean() + 0.05        # the fitted proportion stays at the decoy's level
    cat = simulate_run(tmp_path / "a", seed=1, tracts=TRACTS, fasta_dir=toy_fastas).load()
    adm = paint(cat, Trio("R", "S", "D", False), PROPS)                # positive control: known tracts recovered
    assert covered_bp(adm.tracts, cat, TRACTS) >= 0.8 * TOTAL and adm.post.mean() > adm.decoy_post.mean() + 0.05

def test_outputs(tmp_path, toy_fastas):
    run = simulate_run(tmp_path, seed=1, tracts=TRACTS, inversion=("S", "chr2", 300_001, 700_000),
                       fasta_dir=toy_fastas)
    cat = run.load(); blocks = inverted_blocks(cat)
    p = paint(cat, Trio("R", "S", "D", False), PROPS, blocks=blocks)
    write_paintings([p], cat, run.outdir)
    post = read_tsv(run.outdir / "posterior.R.tsv")
    assert len(post) == cat.member[:, cat.g("R")].sum() and {r["class"] for r in post} <= set(PAINT_CLASSES)
    beds = read_tsv(run.outdir / "tracts.R.bed")
    t = p.tracts[0]
    assert int(beds[0]["start0"]) == t.start - 1 and int(beds[0]["end"]) == t.end and beds[0]["youngest_u95"] == "."
    rd = [p.idx[i] for i in range(t.first, t.last + 1) if PAINT_CLASSES[p.cls[i]] == "RD"]
    assert float(beds[0]["youngest_k"]) == pytest.approx(min(cat.k[i, cat.g("R")] for i in rd))
    (blk,) = [b for b in blocks if {b.x, b.y} == {"R", "S"}]
    flagged = {cat.lid[i] for i in p.idx if cat.seq[i, cat.g("R")] == blk.seq_x
               and blk.start_x <= cat.start[i, cat.g("R")] <= blk.end_x}
    assert {r["lid"] for r in post if r["inv"] == "1"} == flagged and flagged
    assert (run.outdir / "decoy_tracts.S.bed").is_file() and (run.outdir / "decoy_posterior.S.tsv").is_file()

def test_no_decoy_without_sister_annotation(tmp_path, toy_fastas):
    genomes = [GENOMES[0], ("S", "taxS", None, False), GENOMES[2], GENOMES[3]]
    loci = ([loc(1000 * i, "R:P* S:E D:P* O:E") for i in range(1, 40)]
            + [loc(50_000 + 1000 * i, "R:P* S:P D:E O:E") for i in range(40)])
    run = write_run(tmp_path, loci, genomes=genomes, fasta_dir=toy_fastas)
    cat = run.load(); p = paint(cat, Trio("R", "S", "D", True), PROPS)
    write_paintings([p], cat, run.outdir); _, f = paint_rows(p, cat)
    beds = read_tsv(run.outdir / "tracts.R.bed")
    assert len(beds) >= 1 and all(b["q"] == "none" and b["p_decoy"] == "none" for b in beds)
    assert any(t.startswith("WARNING:") and "decoy" in t for _, t in f.items())

def test_tracts_without_k(tmp_path, toy_fastas):
    run = simulate_run(tmp_path, seed=1, tracts=TRACTS, k_known=False, fasta_dir=toy_fastas)
    cat = run.load(); p = paint(cat, Trio("R", "S", "D", False), PROPS)
    write_paintings([p], cat, run.outdir); _, f = paint_rows(p, cat)
    assert p.tracts and {b["youngest_k"] for b in read_tsv(run.outdir / "tracts.R.bed")} == {"."}
    assert sum("youngest insertion per tract skipped" in t for s, t in f.items() if s == "Dating") == 1

def test_paint_many_scaffolds(tmp_path, toy_fastas):
    seqids = ("chr1",) + tuple(f"scaf{i:03d}" for i in range(300))
    loci = ([loc(1000 * i, "R:P* S:E D:P O:E") for i in range(1, 60)]
            + [loc(70_000 + 1000 * i, "R:P* S:P D:E O:E") for i in range(60)]
            + [loc(1001, "R:P* S:E D:P O:E" if i % 2 else "R:P* S:P D:E O:E", seqid=f"scaf{i:03d}") for i in range(300)])
    cat = write_run(tmp_path, loci, seqids=seqids, seq_len=5_000, seq_lens={"chr1": 200_000},
                    fasta_dir=toy_fastas).load()
    p = paint(cat, Trio("R", "S", "D", True), PROPS)
    assert p.skipped is None and len(p.idx) == 419 and np.isfinite(p.post).all()
    assert all(t.start <= t.end and t.first <= t.last for t in p.tracts)  # one-locus chains give well-formed tracts

def test_shortest_callable():
    fit = HmmFit(pi=0.1, L=1e5, e0=np.array([0.1, 0.9]), e1=np.array([0.9, 0.1]), post=np.zeros(0), loglik=0.0)
    assert shortest_callable(fit, 5000.0) == ((2, 10_000.0), (4, 20_000.0))   # worked out by hand


# ---------------------------------------------------------------- the HMM, tracts, q-values, skips and files


def _row(rows, test, statistic, stratum="all"):
    return next(r for r in rows if (r.test, r.stratum, r.statistic) == (test, stratum, statistic))


def test_forward_backward_chains_restart_and_counts():
    """Two chains in one call equal two separate calls (each restarts at pi; log-likelihoods add), and the
    expected emission counts are the posterior summed per class over the emitting loci."""
    pos = np.array([0.0, 2e4, 9e4, 1e3, 5e4]); emit = np.array([0, 0, 1, 1, -1])
    e0, e1 = np.array([0.25, 0.75]), np.array([0.7, 0.3])
    post, ll, c1, c0 = forward_backward(pos, emit, np.array([0, 3, 5]), 0.3, 1 / 5e4, e0, e1)
    a = forward_backward(pos[:3], emit[:3], np.array([0, 3]), 0.3, 1 / 5e4, e0, e1)
    b = forward_backward(pos[3:], emit[3:], np.array([0, 2]), 0.3, 1 / 5e4, e0, e1)
    assert np.allclose(post, np.r_[a[0], b[0]]) and ll == pytest.approx(a[1] + b[1])
    assert np.allclose(c1, [post[0] + post[1], post[2] + post[3]])
    assert np.allclose(c0, [2 - post[0] - post[1], 2 - post[2] - post[3]])
    one = forward_backward(np.array([7.0]), np.array([0]), np.array([0, 1]), 0.3, 1 / 5e4, e0, e1)[0]
    assert one[0] == pytest.approx(0.3 * 0.7 / (0.3 * 0.7 + 0.7 * 0.25))        # a one-locus chain: Bayes


def test_fit_hmm_two_blocks():
    """A donor-like block in a sister-like background: EM keeps state 1 donor-like, pi within its bounds, and
    the block's loci called as one tract."""
    rd = np.zeros(30, dtype=np.int64); rs = np.ones(120, dtype=np.int64)
    emit = np.r_[rs[:60], rd, rs[60:]]
    pos = np.arange(len(emit)) * 2000.0
    e0 = np.array([0.3, 0.7]); e1 = e0 * (5, 0.2); e1 /= e1.sum()
    fit = fit_hmm(pos, emit, np.array([0, len(emit)]), e0, e1)
    assert fit.e1[0] > fit.e0[0] and 1e-4 <= fit.pi <= 0.5 and fit.L in np.geomspace(1e4, 5e7, 41)
    tr = call_tracts(np.zeros(len(emit), dtype=int), pos.astype(int) + 1, pos.astype(int) + 150, fit.post, fit.pi)
    assert [(t.first, t.last) for t in tr] == [(60, 89)]
    with pytest.raises(ValueError, match="no emitting locus"):
        fit_hmm(pos[:3], np.array([-1, -1, -1]), np.array([0, 3]), e0, e1)


def test_fit_hmm_reference_values():
    """test_fit_hmm_two_blocks' input with every third locus missing: fixed reference values that pin the
    emission update and pi averaged over the emitting loci only."""
    emit = np.r_[np.ones(60), np.zeros(30), np.ones(60)].astype(np.int64)
    emit[::3] = -1
    pos = np.arange(len(emit)) * 2000.0
    e0 = np.array([0.3, 0.7]); e1 = e0 * (5, 0.2); e1 /= e1.sum()
    fit = fit_hmm(pos, emit, np.array([0, len(emit)]), e0, e1)
    assert fit.pi == pytest.approx(0.20021184711333548, rel=1e-6)
    assert fit.L == pytest.approx(44393.5372535804, rel=1e-6)
    assert fit.e0 == pytest.approx([0.0062588872383283626, 0.9937411127616715], rel=1e-6)
    assert fit.e1 == pytest.approx([0.9748814354581047, 0.02511856454189529], rel=1e-6)
    assert fit.loglik == pytest.approx(-9.952418349028346, rel=1e-6)


def test_tracts_break_at_sequences():
    """A run of called loci never crosses a sequence boundary."""
    tr = call_tracts(np.array([0, 0, 1, 1]), np.array([1, 1001, 1, 1001]), np.array([500, 1500, 500, 1500]),
                     np.full(4, 0.9), pi=0.5)
    assert [(t.seq, t.first, t.last) for t in tr] == [(0, 0, 1), (1, 2, 3)]


def test_shortest_callable_middle_locus():
    """An asymmetric fit: the rule reads the run's middle locus (k // 2);
    its first locus would need 9 loci between the flanks."""
    fit = HmmFit(pi=0.2, L=2e5, e0=np.array([0.3, 0.7]), e1=np.array([0.6, 0.4]), post=np.zeros(0), loglik=0.0)
    assert shortest_callable(fit, 20_000.0) == ((3, 60_000.0), (7, 140_000.0))


def test_tied_tract_sums_compare_equal():
    """Tract scores are exactly rounded sums: the same posteriors score the same after a long uncalled prefix
    and in reverse order (cumulative-sum differences and np.sum both break this tie), so a tied decoy tract
    counts as >= the target: p_decoy 1, not 1/2."""
    tail = np.array([0.7, 0.9, 0.8, 0.6, 0.95, 0.85])

    def tracts(post):
        n = len(post)
        return call_tracts(np.zeros(n, dtype=int), np.arange(n) * 1000 + 1, np.arange(n) * 1000 + 500, post, 0.2)
    (t,), (d,) = tracts(np.r_[np.full(1000, 0.3), tail]), tracts(tail[::-1].copy())
    assert t.sum_post == d.sum_post == math.fsum(tail.tolist())
    q, p = target_decoy(np.array([t.sum_post]), np.array([d.sum_post]), 1.0)
    assert p[0] == 1.0 and q[0] == 1.0


def test_floor_line_names_each_reason():
    """The +1 corrections' floors, each with its own reason; and "1 locus" in the findings' counts."""
    both = _floor_line("R", 3, 0)
    assert "q >= 1/3" in both and "too few tracts" in both and "p_decoy >= 1/1" in both
    assert "too few decoy tracts" in both
    only_p = _floor_line("R", 100, 5)          # q can reach 0.01, p_decoy only 1/6
    assert "p_decoy >= 1/6" in only_p and "too few decoy tracts" in only_p
    assert "q >=" not in only_p and "too few tracts" not in only_p
    only_q = _floor_line("R", 5, 100)
    assert "q >= 1/5" in only_q and "p_decoy" not in only_q
    assert _floor_line("R", 100, 100) is None and _floor_line("R", 0, 0) is None
    assert (_n_loci(1), _n_loci(7)) == ("1 locus", "7 loci")


def test_collapsed_fit_reported_as_such(tmp_path, toy_fastas):
    """Seed 1's null: EM drives pi to its floor and calls no tract. donor_state, L, the D-like run and the
    shortest callable tract are NA (pseudocount artifacts otherwise), and one line says no D-like ancestry
    is detectable."""
    cat = simulate_run(tmp_path, seed=1, fasta_dir=toy_fastas).load()
    p = paint(cat, Trio("R", "S", "D", False), PROPS)
    assert p.fit.pi == 1e-4 and not p.tracts and p.collapsed
    rows, f = paint_rows(p, cat)
    for s in ("donor_state", "mean_tract_bp", "mean_donor_run_bp", "shortest_callable_loci_lo",
              "shortest_callable_bp_lo", "shortest_callable_loci_hi", "shortest_callable_bp_hi"):
        assert math.isnan(_row(rows, "paint", s).value), s
    assert _row(rows, "paint", "pi").value == 1e-4 and _row(rows, "paint", "n_tracts").value == 0
    lines = [t for s, t in f.items() if s == "Painting"]
    assert sum("the fit collapsed" in t and "no D-like ancestry is detectable" in t for t in lines) == 1
    assert not any("50,000.0 kb" in t or "loci (about" in t or "no D-like state" in t for t in lines)


def test_q_path_scores_sums_at_the_informative_scale(tmp_path, toy_fastas, monkeypatch):
    """paint hands target_decoy each tract's sum_post, the decoy tracts', and n informative R / n
    informative S; q and p are what it returns."""
    seen = []

    def spy(target, decoy, scale):
        seen.append((np.array(target), np.array(decoy), scale))
        return target_decoy(target, decoy, scale)
    monkeypatch.setattr("ltrtrace.analyses.geneflow.paint.target_decoy", spy)
    cat = simulate_run(tmp_path, seed=3, tracts=TRACTS, fasta_dir=toy_fastas).load()
    p = paint(cat, Trio("R", "S", "D", False), PROPS)
    ((target, decoy, scale),) = seen
    assert p.tracts and p.decoy_tracts                         # seed 3: 4 tracts, 1 decoy tract
    assert target.tolist() == [t.sum_post for t in p.tracts]
    assert decoy.tolist() == [t.sum_post for t in p.decoy_tracts]
    n_r, n_s = (int(np.count_nonzero(c <= 1)) for c in (p.cls, p.decoy_cls))
    assert scale == n_r / n_s and scale != 1.0
    q, pv = target_decoy(target, decoy, scale)
    assert p.q == q.tolist() and p.p == pv.tolist()


def test_target_decoy_edges():
    q, p = target_decoy(np.zeros(0), np.array([1.0]), 2.0)
    assert len(q) == len(p) == 0
    q, p = target_decoy(np.array([3.0, 1.0]), np.zeros(0), 1.0)      # no decoy tract: raw 1 / #target >= x
    assert np.allclose(q, [0.5, 0.5]) and np.allclose(p, [1.0, 1.0])
    target = np.array([5.0, 4.0, 3.0, 2.0])
    q, p = target_decoy(target, np.array([3.5]), 2.0)               # the scale multiplies the decoy counts
    assert np.allclose(q, [0.5, 0.5, 0.75, 0.75]) and np.allclose(p, [0.5, 0.5, 1.0, 1.0])
    assert np.allclose(target_decoy(target, np.array([3.5]), 1.0)[0], 0.5)
    _, p = target_decoy(np.array([4.0, 3.0, 2.0, 1.0]), np.array([3.0]), 1.0)   # a decoy tied with a target
    assert np.allclose(p, [0.5, 1.0, 1.0, 1.0])                                 # counts as >= it


def test_classes_recipient_discovery(tmp_path, toy_fastas):
    """Painting reads the recipient's members, each classed by the recipient's own calls."""
    loci = [loc(1000, "R:P* S:E D:P O:E"), loc(2000, "R:P S:P* D:E O:E"),
            loc(3000, "R:P* S:E D:P* O:E", calls={("R", "S"): "U"}), loc(4000, "R:E S:P* D:P O:E")]
    cat = write_run(tmp_path, loci, fasta_dir=toy_fastas).load()
    idx, cls = classes(cat, "R", "S", "D", alphabet="insertion")
    assert idx.tolist() == [0, 2] and [PAINT_CLASSES[c] for c in cls] == ["RD", "missing"]


def test_skips(tmp_path, toy_fastas):
    """No RS locus, and a recipient without annotation: nothing painted, the reason kept, no file written."""
    loci = [loc(1000 * i, "R:P* S:E D:P O:E") for i in range(1, 30)]
    run = write_run(tmp_path / "a", loci, fasta_dir=toy_fastas)
    cat = run.load()
    p = paint(cat, Trio("R", "S", "D", True), PROPS)
    assert p.skipped and "RS" in p.skipped and p.fit is None and not p.tracts
    rows, f = paint_rows(p, cat)
    assert _row(rows, "paint", "n_RD").value == 29 and _row(rows, "paint", "n_RS").value == 0
    assert any(t.startswith("WARNING:") and "skipped" in t for _, t in f.items())
    assert write_paintings([p], cat, run.outdir) == [] and not list(run.outdir.glob("posterior.*"))
    genomes = [("R", "taxR", None, False), GENOMES[1], GENOMES[2], GENOMES[3]]
    loci = [loc(1000 * i, "R:P S:E D:P* O:E") for i in range(1, 10)]
    cat = write_run(tmp_path / "b", loci, genomes=genomes, fasta_dir=toy_fastas).load()
    p = paint(cat, Trio("R", "S", "D", True), PROPS)
    assert p.skipped and "annotation" in p.skipped and len(p.idx) == 0


def test_youngest_u95_from_self_rows(tmp_path, toy_fastas):
    """youngest_u95: the smallest u95 over a tract's RD loci, from R's own copy, else from D's copy; loci
    without a self row left out; its own Dating line."""
    loci = ([loc(1000 * i, "R:P* S:E D:P* O:E", k=0.004) for i in range(1, 40)]
            + [loc(50_000 + 1000 * i, "R:P* S:P D:E O:E") for i in range(40)]
            + [loc(20_500, "R:P* S:P D:P* O:E", k=0.0)])  # inside the tract, not RD: dates nothing
    cat = write_run(tmp_path, loci, fasta_dir=toy_fastas).load()
    lid = cat.lid.tolist()
    counts = {(lid[0], "R"): (300, 2), (lid[0], "D"): (300, 0),     # R's copy first, though D's is younger
              (lid[1], "D"): (400, 1),                              # no R self row: D's copy dates it
              (lid[45], "R"): (500, 0)}                             # an RS locus: never counted
    p = paint(cat, Trio("R", "S", "D", True), PROPS, self_counts=counts)
    assert len(p.tracts) == 1 and p.tracts[0].first == 0
    t = p.tracts[0]
    assert t.last >= 1
    _, f = paint_rows(p, cat)
    out = tmp_path / "run"
    write_paintings([p], cat, out)
    (bed,) = read_tsv(out / "tracts.R.bed")
    assert float(bed["youngest_u95"]) == pytest.approx(min(u95(2, 300), u95(1, 400)), rel=1e-5)
    assert float(bed["youngest_k"]) == pytest.approx(0.004) and bed["n_identical"] == "0"
    assert (bed["n_loci"], bed["n_RD"], bed["n_RS"]) == ("40", "39", "0")
    assert sum("youngest_u95" in t for s, t in f.items() if s == "Dating") == 1


def test_rows_and_files(tmp_path, toy_fastas):
    """paint_rows' statistics, its per-sequence rows and loci declarations; the files' columns; stale files
    removed; the decoy painted with the target's parameters, its scale the ratio of informative loci."""
    run = simulate_run(tmp_path, seed=1, tracts=TRACTS, fasta_dir=toy_fastas)
    cat = run.load()
    p = paint(cat, Trio("R", "S", "D", False), PROPS)
    rows, f = paint_rows(p, cat)
    didx, dcls = classes(cat, "S", "R", "D", alphabet="insertion")  # the decoy: S's members and S's own calls
    assert np.array_equal(p.decoy_idx, didx) and np.array_equal(p.decoy_cls, dcls) and not p.collapsed
    n_inf = int(np.count_nonzero(np.isin(p.cls, (0, 1))))
    assert _row(rows, "paint", "n_informative").value == n_inf
    assert _row(rows, "paint", "mean_tract_bp").value == p.fit.L
    assert _row(rows, "paint", "mean_donor_run_bp").value == pytest.approx(p.fit.L / (1 - p.fit.pi))
    assert _row(rows, "paint", "excess").value == pytest.approx(p.post.mean() - p.decoy_post.mean())
    assert _row(rows, "paint", "n_tracts").value == len(p.tracts)
    assert _row(rows, "paint", "n_tracts_q05").value == sum(q <= 0.05 for q in p.q)
    assert p.scale == pytest.approx(n_inf / np.count_nonzero(np.isin(p.decoy_cls, (0, 1))))
    assert all(r.loci != "." and r.layer == "pav" and r.scope == "R,S;D" for r in rows)
    assert _row(rows, "paint", "n_RD").loci == "disc=R;present=R;res=S,D;outgroup=not_present"
    assert {r.statistic for r in rows if r.test == "paint_resolution"} == {"resolved_S", "resolved_D"}
    seq_rows = [r for r in rows if r.test == "paint_seqid" and r.statistic == "class_share"]
    assert [r.stratum for r in seq_rows] == ["chr1", "chr2", "chr3", "chr4"]
    sp, lo, hi = (_row(rows, "paint", s).value for s in ("spacing_median", "shortest_callable_loci_lo",
                                                           "shortest_callable_loci_hi"))
    assert ((lo, lo * sp), (hi, hi * sp)) == shortest_callable(p.fit, sp)
    pos, starts = _chain(cat, "S", p.decoy_idx)
    emit = np.where(p.decoy_cls <= 1, p.decoy_cls, -1).astype(np.int64)
    decoy = forward_backward(pos, emit, starts, p.fit.pi, 1 / p.fit.L, p.fit.e0, p.fit.e1)[0]
    assert np.allclose(decoy, p.decoy_post)                      # the fitted parameters, unchanged
    stale, other = run.outdir / "tracts.Old.bed", run.outdir / "tracts.mine.bed"   # listed by the last marker; not
    stale.write_text("x\n")
    other.write_text("x\n")
    paths = write_paintings([p], cat, run.outdir, previous=[stale])
    assert not stale.exists() and other.exists() and sorted(x.name for x in paths) == [
        "decoy_posterior.S.tsv", "decoy_tracts.S.bed", "posterior.R.tsv", "tracts.R.bed"]
    with open(run.outdir / "tracts.R.bed") as fh:
        assert tuple(fh.readline().rstrip("\n").split("\t")) == TRACT_COLUMNS
    with open(run.outdir / "posterior.R.tsv") as fh:
        assert tuple(fh.readline().rstrip("\n").split("\t")) == POSTERIOR_COLUMNS
    beds = read_tsv(run.outdir / "tracts.R.bed")
    assert len(beds) == len(p.tracts) and all(b["name"] == "D-like" and b["strand"] == "." for b in beds)
    assert [int(b["score"]) for b in beds] == [min(1000, int(1000 * t.mean_post)) for t in p.tracts]
    assert len(read_tsv(run.outdir / "decoy_posterior.S.tsv")) == len(p.decoy_idx)
    assert any(s == "Painting" and "D-like posterior mean" in t for s, t in f.items())
    assert _row(rows, "paint", "n_RD").loci == "disc=R;present=R;res=S,D;outgroup=not_present"


def _chain(cat, x, idx):
    """(pos, starts) of a painting's loci: x's site starts, one chain per sequence (paint's own layout)."""
    seq, pos = cat.seq[idx, cat.g(x)], cat.start[idx, cat.g(x)].astype(np.float64)
    return pos, np.r_[0, np.flatnonzero(np.diff(seq)) + 1, len(idx)].astype(np.int64)


def test_write_paintings_removes_only_painting_files(tmp_path, toy_fastas):
    """The analyses remove only the painting files named for the run's labels and those the last run recorded
    writing, never another file: a genome labeled `tracts` keeps the route's tracts.pav.bed, a user's
    tracts.mine.bed stays, and a recorded painting of a genome no longer in the sheet goes."""
    from ltrtrace.analyses.geneflow import run as geneflow
    from ltrtrace.analyses.geneflow.paint import PAINTING_FILES

    genomes = [("R", "tR", None, True), ("S", "tS", None, True), ("D", "tD", None, True),
               ("tracts", "tO", "outgroup", False)]
    rd, sd, rs = "R:P* S:E D:P tracts:E", "R:E S:P* D:P tracts:E", "R:P* S:P D:E tracts:E"
    loci = [loc(5_000 + 6_000 * j, (rd, sd, rs, rs)[j % 4], seqid=s) for s in ("chr1", "chr2") for j in range(60)]
    run = write_run(tmp_path, loci, genomes=genomes, fasta_dir=toy_fastas)
    out = run.outdir
    keep = [out / "tracts.pav.bed", out / "tracts.mine.bed"]            # the engine's, the user's
    gone = out / "tracts.X.bed"                                         # a painting the last run listed
    for p in (*keep, gone):
        p.write_text("x\n")
    record = out / geneflow.PAINTINGS
    record.parent.mkdir(parents=True)
    record.write_text(json.dumps(["tracts.X.bed"]))
    geneflow.analyses(run.sheet, out, 1, trio=("R", "S", "D"))
    assert all(p.read_text() == "x\n" for p in keep) and not gone.exists()
    assert (out / "tracts.R.bed").is_file() and (out / "decoy_tracts.S.bed").is_file()
    assert json.loads(record.read_text()) == ["tracts.R.bed", "posterior.R.tsv", "decoy_tracts.S.bed",
                                              "decoy_posterior.S.tsv"]
    assert PAINTING_FILES[0].format("tracts") == "tracts.tracts.bed"


# ---------------------------------------------------------------- the decoy-anchored emission model


def _emission_rows(rows):
    """{statistic: row} of a painting's emission row and its four e* rows."""
    return {r.statistic: r for r in rows if r.test == "paint" and r.stratum == "all"
            and r.statistic in ("emission", "e0_RD", "e0_RS", "e1_RD", "e1_RS")}


def _segregated(tmp_path, toy_fastas):
    """Sparse, perfectly class-segregated loci on each of the four toy sequences (4.1 Mb): 10 RD runs and 10 RS runs
    of 3 R loci, alternating, runs 200 kb apart (loci 10 kb apart inside a run); one S locus midway between
    consecutive runs, alternately D-like and R-like: a decoy of 40 D-like and 40 R-like loci, no class bias. The
    pattern is repeated on four sequences because the EM update's pseudocounts cap an emission at (n + 0.5) / (n + 1):
    30 loci of a class reach only 0.984, 120 reach 0.996."""
    rd, rs = "R:P* S:E D:P O:E", "R:P* S:P D:E O:E"
    sd, sr = "R:E S:P* D:P O:E", "R:P S:P* D:E O:E"
    loci = []
    for seqid in ("chr1", "chr2", "chr3", "chr4"):
        for r in range(20):
            start = 100_000 + 200_000 * r
            loci += [loc(start + 10_000 * j, rd if r % 2 == 0 else rs, seqid=seqid) for j in range(3)]
            loci.append(loc(start + 100_000, sd if r % 2 == 0 else sr, seqid=seqid))
    return write_run(tmp_path, loci, seq_len=4_100_000, fasta_dir=toy_fastas)


def test_em_degenerates_on_segregated_loci(tmp_path, toy_fastas):
    """On sparse, class-segregated loci the EM fit degenerates: two pure states (e1[RD] > 0.99, e0[RS] > 0.99), so
    the posteriors copy the classes."""
    cat = _segregated(tmp_path, toy_fastas).load()
    em = paint(cat, Trio("R", "S", "D", True), PROPS)
    assert em.fit.e1[0] > 0.99 and em.fit.e0[1] > 0.99
    inf = em.cls <= 1
    assert not np.any((em.post[inf] > 0.05) & (em.post[inf] < 0.95))       # 0/1 copies of the classes


def test_emission_rows(tmp_path, toy_fastas):
    """Every painted trio has the emission row and the e0_RD, e0_RS, e1_RD, e1_RS rows, beside donor_state: the fitted
    emissions, each pair summing to 1. NA for a collapsed fit (seed 1's null) and for a skipped painting."""
    cat = simulate_run(tmp_path / "a", seed=1, tracts=TRACTS, fasta_dir=toy_fastas).load()
    p = paint(cat, Trio("R", "S", "D", False), PROPS)
    rows, _ = paint_rows(p, cat)
    e = _emission_rows(rows)
    assert not p.collapsed and e["emission"].value == "em"
    assert [e[s].value for s in ("e0_RD", "e0_RS", "e1_RD", "e1_RS")] == [*p.fit.e0.tolist(), *p.fit.e1.tolist()]
    assert e["e0_RD"].value + e["e0_RS"].value == pytest.approx(1.0)
    assert e["e1_RD"].value + e["e1_RS"].value == pytest.approx(1.0)
    names = [r.statistic for r in rows if r.test == "paint"]
    i = names.index("donor_state")
    assert names[i + 1:i + 6] == ["emission", "e0_RD", "e0_RS", "e1_RD", "e1_RS"]
    null = simulate_run(tmp_path / "n", seed=1, fasta_dir=toy_fastas).load()
    p = paint(null, Trio("R", "S", "D", False), PROPS)
    assert p.collapsed
    e = _emission_rows(paint_rows(p, null)[0])
    assert e["emission"].value == "em" and all(math.isnan(e[s].value) for s in ("e0_RD", "e0_RS", "e1_RD", "e1_RS"))
    loci = [loc(1000 * i, "R:P* S:E D:P O:E") for i in range(1, 30)]       # no RS locus: skipped
    cat = write_run(tmp_path / "s", loci, fasta_dir=toy_fastas).load()
    p = paint(cat, Trio("R", "S", "D", True), PROPS)
    assert p.skipped
    e = _emission_rows(paint_rows(p, cat)[0])
    assert all(math.isnan(e[s].value) for s in ("e0_RD", "e0_RS", "e1_RD", "e1_RS"))


def test_no_decoy_warns_once(tmp_path, toy_fastas):
    """Without a decoy the painting still fits, writes no q-values and gives one WARNING naming why: a sister
    without annotation, and a sister whose loci hold no informative locus."""
    genomes = [GENOMES[0], ("S", "taxS", None, False), GENOMES[2], GENOMES[3]]
    loci = ([loc(1000 * i, "R:P* S:E D:P* O:E") for i in range(1, 40)]
            + [loc(50_000 + 1000 * i, "R:P* S:P D:E O:E") for i in range(40)])
    cat = write_run(tmp_path / "a", loci, genomes=genomes, fasta_dir=toy_fastas).load()
    loci_b = loci + [loc(200_000 + 1000 * i, "R:E S:P* D:E O:E") for i in range(10)]   # S's own loci: private
    cat_b = write_run(tmp_path / "b", loci_b, fasta_dir=toy_fastas).load()
    for c, why in ((cat, "S has no annotation"), (cat_b, "no informative")):
        p = paint(c, Trio("R", "S", "D", True), PROPS)
        assert p.scale is None and p.skipped is None and all(q is None for q in p.q)
        f = paint_rows(p, c)[1]
        warned = [t for _, t in f.items() if t.startswith("WARNING:") and "no decoy" in t]
        assert len(warned) == 1 and why in warned[0], warned
        assert [t for _, t in f.items() if t.startswith("WARNING:") and why in t] == warned   # one warning, not two
