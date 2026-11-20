"""Sequence-identity layer (analyses.geneflow.identity): exact LTR counts through the route's aligner (LASTZ),
the rate-free split (isolation fit, mixture guard and CDF check), the layer end to end on toy runs with real
sequence (catalog_toy.py), its alignment checkpoint and the pair-divergence layer it writes."""
from __future__ import annotations

import gzip
import math
import os
import random
import signal

import builders
import numpy as np
import pytest
from analyses.geneflow.catalog_toy import GENOMES, block_size, loc, read_tsv, write_run
from markers import needs_lastz
from scipy.special import gammaincc

import ltrtrace.analyses.geneflow.identity as identity_module
from ltrtrace.analyses.geneflow.identity import (
    IDENTITY_COLUMNS,
    _elements,
    _job,
    _mix_params,
    fit_mix,
    fit_split,
    identity,
    ltr_counts,
    pred_cdf,
    self_counts,
    split_nll,
    write_identity,
    write_pair_divergences,
)
from ltrtrace.analyses.geneflow.relations import pav_distances
from ltrtrace.analyses.geneflow.trio import Trio
from ltrtrace.base.constants import IDENT_MIN_SITES, MIX_MIN_PAIRS, WINDOW_PAD
from ltrtrace.calls.records import CALL_COLUMNS, read_calls_meta, read_divergences, write_divergences
from ltrtrace.io.seqio import revcomp
from ltrtrace.seq.divergence import garwood_u95 as u95


def substitute(seq: str, positions: list[int]) -> str:
    """Each listed base replaced by the next one in ACGT order (always a change)."""
    s = list(seq)
    for p in positions:
        s[p] = "ACGT"[("ACGT".index(s[p]) + 1) % 4]
    return "".join(s)

INT_POSITIONS_12 = [500 + 250 * i for i in range(12)]           # inside the internal region [400, 3400)

@needs_lastz
def test_ltr_counts_exact():
    rng = random.Random(5); e = builders.make_element(rng, 400, 3000)[0]          # LTRs [0, 400) and [3400, 3800)
    y = builders.rand_seq(800, rng) + substitute(e, [10, 50, 90, 120, 3405, 3700] + INT_POSITIONS_12) + builders.rand_seq(800, rng)
    assert ltr_counts(e, y, [(0, 400), (3400, 3800)]) == (800, 6)

@needs_lastz
def test_self_counts():
    rng = random.Random(6); ltr = builders.rand_seq(400, rng)
    assert self_counts(ltr, substitute(ltr, [20, 200, 380])) == (400, 3)

def test_fit_split_recovers():
    rng = np.random.default_rng(1); n = np.full(50_000, 800.0)
    k = rng.poisson(n * (0.002 + rng.exponential(0.01, n.size)))
    s, m = fit_split(k, n)
    assert s == pytest.approx(0.002, rel=0.15) and m == pytest.approx(0.01, rel=0.15)

@pytest.mark.slow          # about 5 min (18 Nelder-Mead starts on 50,000 pairs; 10,000 is too few): CI deselects it
def test_fit_mix_guards_against_gene_flow():
    rng = np.random.default_rng(2); n = np.full(50_000, 800.0)
    flow = rng.random(n.size) < 0.04                                          # 4% post-split gene flow
    k = rng.poisson(n * np.where(flow, rng.exponential(3e-4, n.size), 0.002 + rng.exponential(0.01, n.size)))
    s_iso, m_iso = fit_split(k, n)
    w, _m_flow, s_mix, _m_mix, _m_tail, _nll = fit_mix(k, n, s_iso, m_iso)
    assert abs(s_mix - 0.002) < abs(s_iso - 0.002) and 0.01 < w[1] < 0.1    # fit: w_flow 0.042, s 0.00210 vs 0.00115

def test_pred_cdf_at_zero():
    n = np.full(1000, 800.0)                                                  # P(k = 0) = e^(-ns) / (1 + nm)
    assert pred_cdf(0.0, n, 0.002, 0.01) == pytest.approx(math.exp(-1.6) / 9, abs=0.01)

@needs_lastz
def test_identity_end_to_end(tmp_path, toy_fastas):
    rng = random.Random(9)
    els = [builders.make_element(rng, 400, 2000)[0] for _ in range(7)]          # 2,800 bp; LTRs [0, 400) and [2400, 2800)
    flanks = [builders.rand_seq(2000, rng) for _ in range(8)]
    copies = {"D": els,
              "R": [substitute(e, list(range(10, 10 * i + 1, 10))) for i, e in enumerate(els[:6])] + [els[6]],
              "S": [substitute(e, [2500]) for e in els[:6]] + [builders.rand_seq(2800, rng)]}   # S's 7th site: unrelated
    fastas = {}
    for g, cs in copies.items():
        fastas[g] = tmp_path / f"{g}.fa"
        fastas[g].write_text(">chr1\n" + "".join(f + c for f, c in zip(flanks, cs)) + flanks[7] + "\n")
    pos = [2001 + 4800 * i for i in range(7)]
    loci = ([loc(p, "R:P* S:P* D:P* O:E", length=2800, ltr=400) for p in pos[:6]]
            + [loc(pos[6], "R:P* S:P D:P* O:E", length=2800, ltr=400, strand={"D": "."})])
    cat = write_run(tmp_path / "run", loci, seqids=("chr1",), fastas=fastas, fasta_dir=toy_fastas).load()
    rows, f, out = identity(cat, Trio("R", "S", "D", True), block_size(cat), 2, pav_distances(cat)[0])
    pairs = sorted((r[1], r[4], r[5], r[6]) for r in out if r[0] == "pair")       # (lid, partner, sites, diffs)
    assert [(s, d) for _, y, s, d in pairs if y == "R"] == [(800, i) for i in range(6)] + [(800, 0)]
    # the 7th S pair, D's strand ".", is aligned on its call's b_strand
    # ("+"), not skipped; its site is unrelated, so it holds fewer than IDENT_MIN_SITES LTR sites (div ".")
    s_pairs = [(s, d) for _, y, s, d in pairs if y == "S"]
    assert s_pairs[:6] == [(800, 1)] * 6 and len(s_pairs) == 7 and s_pairs[6][0] < IDENT_MIN_SITES
    assert [r[7] for r in out if r[0] == "pair" and r[4] == "S"][6] is None
    selfs = {(r[3], r[1]): (r[5], r[6]) for r in out if r[0] == "self"}          # (genome, lid) -> (sites, diffs)
    assert len(selfs) == 20 and {s for s, _ in selfs.values()} == {400}
    assert [selfs[("R", lid)][1] for lid in cat.lid] == [0, 1, 2, 3, 4, 5, 0]
    assert {selfs[("S", lid)][1] for lid in cat.lid[:6]} == {1} and {selfs[("D", lid)][1] for lid in cat.lid} == {0}
    assert rows and all(r.layer == "identity" for r in rows)
    assert not any("pair(s) skipped" in t for _, t in f.items())     # every P call carries its b_strand


# ---------------------------------------------------------------- counts, fits, strands and failures


@needs_lastz
def test_ltr_counts_skips_ambiguous_columns():
    """A column with N, or another ambiguity code, on either side is neither a site nor a substitution."""
    rng = random.Random(5)
    e = builders.make_element(rng, 400, 3000)[0]
    y = list(substitute(e, [10, 3405]))
    for p, b in ((10, "N"), (30, "N"), (3500, "N"), (3600, "R")):   # the N at 10 hides a substitution
        y[p] = b
    y = builders.rand_seq(800, rng) + "".join(y) + builders.rand_seq(800, rng)
    assert ltr_counts(e, y, [(0, 400), (3400, 3800)]) == (796, 1)
    assert ltr_counts(e, y, [(0, 400), (300, 450)]) == (448, 0)     # overlapping intervals count once


def _ref_nll(par, k, n):
    """Reference negative log-likelihood of the rate-free split, written apart from split_nll, with Q
    clamped at 1e-300."""
    s, m = np.exp(par)
    c = 1.0 + 1.0 / (n * m)
    q = gammaincc(k + 1, c * n * s)
    lp = s / m - np.log(m * n) - (k + 1) * np.log(c) + np.log(np.maximum(q, 1e-300))
    return -lp.sum()


def test_split_nll_matches_reference_and_is_bounded():
    """The reference's value wherever Q does not underflow; the Poisson limit, not -inf, as m -> 0 (the
    reference's clamp of Q at 1e-300 lets s/m grow without bound there), so a fit to data without spread stays
    finite."""
    rng = np.random.default_rng(3)
    n = rng.integers(200, 3000, 2000).astype(float)
    k = rng.poisson(n * (0.001 + rng.exponential(0.005, n.size))).astype(float)
    for par in ([math.log(1e-3), math.log(1e-2)], [math.log(5e-4), math.log(3e-3)], [math.log(0.02), math.log(1)]):
        assert split_nll(np.array(par), k, n) == _ref_nll(np.array(par), k, n)
    k1, n1 = np.ones(6), np.full(6, 800.0)                     # no spread: the likelihood peaks as m -> 0
    poisson = -6 * math.log(math.exp(-1.0))                   # k = 1 at mean n s = 1, six times: 6.0
    for log_m in (-14.0, -20.0, -40.0):
        assert split_nll(np.array([math.log(1 / 800), log_m]), k1, n1) == pytest.approx(poisson, abs=1e-4)
    assert _ref_nll(np.array([math.log(1 / 800), -20.0]), k1, n1) < -1e6   # the clamp's runaway
    s, m = fit_split(k1, n1)
    assert s == pytest.approx(1 / 800, rel=1e-3) and 0 < m < 1e-6


def test_identity_skipped_without_donor_annotation(tmp_path, toy_fastas):
    genomes = [g if g[0] != "D" else ("D", "taxD", None, False) for g in GENOMES]
    cat = write_run(tmp_path, [loc(10_001, "R:P* S:E D:P O:E")], genomes=genomes, fasta_dir=toy_fastas).load()
    rows, f, out = identity(cat, Trio("R", "S", "D", True), block_size(cat), 1, pav_distances(cat)[0])
    assert rows == [] and out == []
    assert len(f.items()) == 1 and "skipped: D has no annotation" in f.items()[0][1]


def test_write_identity(tmp_path):
    rows = [("self", "L0000001", "R,S;D", "D", ".", 400, 0, 0.0, 0.00748933, 1),
            ("pair", "L0000001", "R,S;D", "D", "S", 150, 2, None, None, None)]
    write_identity(rows, tmp_path / "identity.tsv")
    back = read_tsv(tmp_path / "identity.tsv")
    assert tuple(back[0]) == IDENTITY_COLUMNS
    assert [list(r.values()) for r in back] == [
        ["self", "L0000001", "R,S;D", "D", ".", "400", "0", "0", "0.00748933", "1"],
        ["pair", "L0000001", "R,S;D", "D", "S", "150", "2", ".", ".", "."]]
    with pytest.raises(ValueError, match="IDENTITY_COLUMNS"):
        write_identity([("self", "L1")], tmp_path / "bad.tsv")


@needs_lastz
def test_identity_young_and_below_split(tmp_path, toy_fastas):
    """S-D pairs 8 substitutions in 800 LTR sites apart fix s = 0.01 (no spread: m -> 0). R+D loci whose D copy
    has identical 400-bp LTRs (u95 0.0075 < s) are young; S+D loci whose D copy has 3 LTR-LTR differences (u95
    0.019) are not. One S copy sits on the reverse strand of its window (b_strand '-'), as the engine reports
    it: the pair is read from the reverse complement."""
    rng = random.Random(11)
    els = [builders.make_element(rng, 400, 2000)[0] for _ in range(10)]
    flanks = [builders.rand_seq(2000, rng) for _ in range(11)]
    sub8 = [5, 60, 130, 250, 2450, 2600, 2700, 2790]             # 4 + 4 LTR positions
    d_copy = els[:7] + [substitute(e, [2410, 2500, 2610]) for e in els[7:]]   # 3 LTR-LTR differences at S+D loci
    r_copy = list(d_copy[:7]) + [builders.rand_seq(2800, rng) for _ in els[7:]]
    s_copy = [substitute(e, sub8) for e in d_copy[:4]] + [builders.rand_seq(2800, rng) for _ in els[4:7]] + [
        substitute(e, sub8) for e in d_copy[7:]]
    s_copy[1] = revcomp(s_copy[1])                               # S's second copy: reverse strand
    fastas = {}
    for g, cs in (("D", d_copy), ("R", r_copy), ("S", s_copy)):
        fastas[g] = tmp_path / f"{g}.fa"
        fastas[g].write_text(">chr1\n" + "".join(f + c for f, c in zip(flanks, cs)) + flanks[10] + "\n")
    pos = [2001 + 4800 * i for i in range(10)]
    loci = ([loc(p, "R:P S:P D:P* O:E", length=2800, ltr=400, strand={"S": "-"} if j == 1 else {})
             for j, p in enumerate(pos[:4])]
            + [loc(p, "R:P S:E D:P* O:E", length=2800, ltr=400) for p in pos[4:7]]     # R+D
            + [loc(p, "R:E S:P D:P* O:E", length=2800, ltr=400) for p in pos[7:]])     # S+D
    cat = write_run(tmp_path / "run", loci, seqids=("chr1",), fastas=fastas, fasta_dir=toy_fastas).load()
    rows, f, out = identity(cat, Trio("R", "S", "D", True), block_size(cat), 1, pav_distances(cat)[0])
    val = {r.statistic: r for r in rows}
    sd = sorted((r[1], r[5], r[6]) for r in out if r[0] == "pair" and r[4] == "S")
    assert [(s, d) for _, s, d in sd] == [(800, 8)] * 7                    # 4 shared + 3 S+D loci
    assert val["n_pairs"].value == 7 and val["s"].value == pytest.approx(0.01, rel=1e-3)
    assert val["s"].scope == "S;D" and val["n_RD_young"].scope == "R,S;D"
    assert (val["n_RD_young"].value, val["n_RD_young"].n) == (3, 3)
    assert (val["n_SD_young"].value, val["n_SD_young"].n) == (0, 3)
    assert val["D_young"].value == 1.0 and math.isnan(val["D_young"].z)      # an empty class: not testable
    assert val["sites_needed"].value == pytest.approx(u95(np.array([0]), np.array([1.0]))[0] / val["s"].value)
    assert val["frac_RD_below"].value == 1.0 and val["frac_RD_below"].n == 7  # R copies identical to D's
    young = {r[1]: r[9] for r in out if r[0] == "self" and r[3] == "D"}
    assert [young[lid] for lid in cat.lid] == [None] * 4 + [1] * 3 + [0] * 3
    d_self = {r[1]: (r[5], r[6], r[8]) for r in out if r[0] == "self" and r[3] == "D"}
    assert d_self[cat.lid[9]][:2] == (400, 3) and d_self[cat.lid[9]][2] == pytest.approx(u95(np.array([3]),
                                                                                            np.array([400]))[0])
    assert all(r.layer == "identity" for r in rows) and IDENT_MIN_SITES == 200
    # with no spread, the gene-flow and tail weights are driven to their bound: their means are NA
    assert val["w_flow"].value * 7 < MIX_MIN_PAIRS and val["w_tail"].value * 7 < MIX_MIN_PAIRS
    assert math.isnan(val["m_flow"].value) and math.isnan(val["m_tail"].value)
    assert "component is empty" in val["m_tail"].note and val["s_mix"].value == pytest.approx(0.01, rel=1e-3)
    assert any("empty, its parameters NA" in t and "gene flow, tail" in t for _, t in f.items())
    assert not any("orientation not inferred" in t for _, t in f.items())


@needs_lastz
def test_identity_missing_fasta_is_one_runtime_error(tmp_path, toy_fastas):
    """A FASTA removed after the catalog was loaded ends the layer with one RuntimeError naming it (it fails while
    the alignment is keyed, before any worker starts), which the CLI reports as one error line."""
    rng = random.Random(13)
    e = builders.make_element(rng, 400, 2000)[0]
    flank = builders.rand_seq(2000, rng)
    fastas = {}
    for g in ("R", "S", "D"):
        fastas[g] = tmp_path / f"{g}.fa"
        fastas[g].write_text(">chr1\n" + flank + e + flank + "\n")
    cat = write_run(tmp_path / "run", [loc(2001, "R:P S:P D:P* O:E", length=2800, ltr=400)], seqids=("chr1",),
                    fastas=fastas, fasta_dir=toy_fastas).load()
    fastas["S"].unlink()
    with pytest.raises(RuntimeError, match="identity layer: an alignment failed: FileNotFoundError"):
        identity(cat, Trio("R", "S", "D", True), block_size(cat), 2, pav_distances(cat)[0])


def _three_loci(tmp_path, toy_fastas):
    """Three loci D discovers, each with D's P calls into R and S: three alignment jobs, so with two threads three
    chunks, each run in a forked worker."""
    return write_run(tmp_path, [loc(p, "R:P S:P D:P* O:E", length=2800, ltr=400) for p in (2001, 10_001, 20_001)],
                     fasta_dir=toy_fastas).load()


def _plant(monkeypatch, fail):
    """identity._locus replaced: in this process it gives zero counts, in a forked worker it calls `fail`."""
    parent = os.getpid()

    def locus(job):
        if os.getpid() != parent:
            fail()
        selfs, pair = job
        return tuple((0, 0) for _ in selfs), tuple((0, 0, 1) for _ in (pair[5] if pair else ()))
    monkeypatch.setattr(identity_module, "_locus", locus)


def test_identity_worker_error_is_one_runtime_error(tmp_path, toy_fastas, monkeypatch):
    """A failure inside a forked worker (here a planted error in its alignment) ends the layer with one
    RuntimeError naming it, which the CLI reports as one error line."""
    cat = _three_loci(tmp_path, toy_fastas)

    def fail():
        raise ValueError("a planted alignment error")
    _plant(monkeypatch, fail)
    with pytest.raises(RuntimeError, match="identity layer: an alignment failed: ValueError: a planted alignment error"):
        identity(cat, Trio("R", "S", "D", True), block_size(cat), 2, pav_distances(cat)[0])


def test_identity_dead_worker_is_one_runtime_error(tmp_path, toy_fastas, monkeypatch):
    """A worker that dies (here it kills itself mid-alignment) ends the layer with one RuntimeError naming the layer
    and the dead worker's chunk (base.runtime.PoolError), never a hang."""
    cat = _three_loci(tmp_path, toy_fastas)
    _plant(monkeypatch, lambda: os.kill(os.getpid(), signal.SIGKILL))
    with pytest.raises(RuntimeError, match=r"identity layer: a worker process died or could not start .*; item [012] "):
        identity(cat, Trio("R", "S", "D", True), block_size(cat), 2, pav_distances(cat)[0])
    assert not (cat.outdir / "work" / "geneflow").exists() or not list((cat.outdir / "work" / "geneflow").glob("*.done.json"))


def test_identity_aligns_every_chunk_in_a_worker(tmp_path, toy_fastas, monkeypatch):
    """Several chunks and threads: this process aligns no locus, every chunk runs in a forked worker (as soon as one
    is free, none waiting for this process)."""
    cat = _three_loci(tmp_path / "run", toy_fastas)
    _plant(monkeypatch, lambda: None)  # zero counts, in a worker as in this process
    ran, real, parent = tmp_path / "ran", identity_module._locus, os.getpid()

    def locus(job):
        with open(ran, "a") as fh:  # the process that aligns each locus
            fh.write(f"{os.getpid()}\n")
        return real(job)

    monkeypatch.setattr(identity_module, "_locus", locus)
    identity(cat, Trio("R", "S", "D", True), block_size(cat), 2, pav_distances(cat)[0])
    pids = [int(x) for x in ran.read_text().split()]
    assert len(pids) == 3 and parent not in pids


def _genome(flanks: list[str], copies: list[str]) -> list[str]:
    """flank 0, copy 0, flank 1, copy 1, ..., the last flank: one sequence as a list of bases."""
    return list("".join(f + c for f, c in zip(flanks, copies)) + flanks[len(copies)])


def _reverse_window(seq: list[str], start: int, length: int) -> None:
    """Reverse-complement the element at 1-based `start` (`length` bp) with WINDOW_PAD bp of flank on each side,
    in place: a partner window the engine reports on the reverse strand, flanks included."""
    a, b = start - 1 - WINDOW_PAD, start - 1 + length + WINDOW_PAD
    seq[a:b] = list(revcomp("".join(seq[a:b])))


@needs_lastz
def test_identity_pair_strands(tmp_path, toy_fastas):
    """The orientation of each pair is its call's b_strand, never the product of the two element strands. D on
    '-' with the partner on '-' (same orientation, b_strand '+'); D on '-' with the partner on '+' (b_strand
    '-'); D '.' with the partner's window reversed (b_strand '-'). Each S copy differs from D's at 2, 3 and 4
    LTR positions; reading the strand product as the relative strand would flip the wrong windows (0 sites)."""
    rng = random.Random(17)
    els = [builders.make_element(rng, 400, 2000)[0] for _ in range(3)]
    flanks = [builders.rand_seq(2000, rng) for _ in range(4)]
    subs = [[20, 2420], [20, 300, 2700], [5, 150, 2405, 2795]]
    pos = [2001 + 4800 * i for i in range(3)]
    d_seq = _genome(flanks, els)
    s_seq = _genome(flanks, [substitute(e, p) for e, p in zip(els, subs)])
    _reverse_window(s_seq, pos[1], 2800)
    _reverse_window(s_seq, pos[2], 2800)
    fastas = {"D": tmp_path / "D.fa", "S": tmp_path / "S.fa", "R": tmp_path / "R.fa"}
    fastas["D"].write_text(">chr1\n" + "".join(d_seq) + "\n")
    fastas["S"].write_text(">chr1\n" + "".join(s_seq) + "\n")
    fastas["R"].write_text(">chr1\n" + "".join(d_seq) + "\n")
    loci = [loc(pos[0], "R:E S:P D:P* O:E", length=2800, ltr=400, strand={"D": "-", "S": "-"}),
            loc(pos[1], "R:E S:P D:P* O:E", length=2800, ltr=400, strand={"D": "-"}),
            loc(pos[2], "R:E S:P D:P* O:E", length=2800, ltr=400, strand={"D": "."})]
    inv = ("S", "chr1", pos[2] - WINDOW_PAD, pos[2] + 2799 + WINDOW_PAD)   # symmetric: the site keeps its place
    cat = write_run(tmp_path / "run", loci, seqids=("chr1",), fastas=fastas, fasta_dir=toy_fastas,
                    inversion=inv).load()
    gs, gd = cat.g("S"), cat.g("D")
    assert cat.pbstrand["D"][:, gs].tolist() == [1, -1, -1]
    assert cat.strand[:, gd].tolist() == [-1, -1, 0]
    _rows, f, out = identity(cat, Trio("R", "S", "D", True), block_size(cat), 1, pav_distances(cat)[0])
    assert [(r[5], r[6]) for r in out if r[0] == "pair"] == [(800, 2), (800, 3), (800, 4)]
    assert not any("pair(s) skipped" in t for _, t in f.items())


@needs_lastz
def test_identity_counts_a_call_without_b_strand(tmp_path, toy_fastas):
    """A P call whose b_strand is unusable cannot be oriented: the pair is skipped, counted and reported
    with a WARNING, never aligned on a guessed strand."""
    rng = random.Random(19)
    e = builders.make_element(rng, 400, 2000)[0]
    flank = builders.rand_seq(2000, rng)
    fastas = {}
    for g in ("R", "S", "D"):
        fastas[g] = tmp_path / f"{g}.fa"
        fastas[g].write_text(">chr1\n" + flank + e + flank + "\n")
    run = write_run(tmp_path / "run", [loc(2001, "R:P S:P D:P* O:E", length=2800, ltr=400, strand={"D": "."})],
                    seqids=("chr1",), fastas=fastas, fasta_dir=toy_fastas)
    pair = run.outdir / "pairs" / "D__S.calls.tsv.gz"
    head, row = gzip.decompress(pair.read_bytes()).decode().splitlines()
    fields = row.split("\t")
    fields[CALL_COLUMNS.index("b_strand")] = "."
    pair.write_bytes(gzip.compress((head + "\n" + "\t".join(fields) + "\n").encode()))
    cat = run.load()
    assert cat.pbstrand["D"][0, cat.g("S")] == 0 and cat.pbstrand["D"][0, cat.g("R")] == 1
    _rows, f, out = identity(cat, Trio("R", "S", "D", True), block_size(cat), 1, pav_distances(cat)[0])
    assert [(r[4], r[5], r[6]) for r in out if r[0] == "pair"] == [("R", 800, 0)]
    assert any(t.startswith("WARNING: 1 pair(s) skipped: no usable b_strand") for _, t in f.items())


def test_job_pads_the_donor_copy(tmp_path, toy_fastas):
    """x is D's element +- WINDOW_PAD, its LTR intervals shifted by the left pad actually fetched (clipped at
    a sequence start); the partner's site keeps its +- WINDOW_PAD."""
    cat = write_run(tmp_path, [loc(10_001, "R:P S:E D:P* O:E", length=2800, ltr=400),
                               loc(201, "R:P S:E D:P* O:E", length=2800, ltr=400, seqid="chr2")],
                    fasta_dir=toy_fastas).load()
    gd = cat.g("D")
    loci = np.arange(2)
    elements = {"D": _elements(cat, "D", loci)}
    pmask = {"R": cat.member[:, gd].copy()}
    _, (d, seqid, x0, x1, ltrs, ys) = _job(cat, 0, ["D"], elements, pmask, "D")
    assert (d, seqid, x0, x1) == ("D", "chr1", 10_000 - WINDOW_PAD, 12_800 + WINDOW_PAD)
    assert ltrs == ((WINDOW_PAD, WINDOW_PAD + 400), (WINDOW_PAD + 2400, WINDOW_PAD + 2800))
    assert ys == (("R", "chr1", 10_000 - WINDOW_PAD, 12_800 + WINDOW_PAD, 1),)
    _, (_, _, x0, _, ltrs, _) = _job(cat, 1, ["D"], elements, pmask, "D")
    assert x0 == 0 and ltrs == ((200, 600), (2600, 3000))                # only 200 bp of left flank exist


def test_mix_params_empty_components():
    """A component whose weight is worth fewer than MIX_MIN_PAIRS of the pairs is at its weight's bound:
    its parameters are NA, the weights themselves kept."""
    w = np.array([0.966089, 0.033911, 0.0])                                  # a 6,000-pair fit
    par, empty = _mix_params(w, 0.0004, 0.002, 0.01, 6.4e284, 6000)
    assert empty == ["tail"] and math.isnan(par["m_tail"])
    assert (par["s_mix"], par["m_mix"], par["m_flow"]) == (0.002, 0.01, 0.0004)
    par, empty = _mix_params(np.array([0.534, 5.3e-15, 0.466]), 0.00049, 0.0025, 0.0148, 0.0149, 23_739)
    assert empty == ["gene flow"] and math.isnan(par["m_flow"]) and par["m_tail"] == 0.0149
    par, empty = _mix_params(np.array([1e-9, 0.5, 0.5 - 1e-9]), 0.0005, 0.002, 0.01, 0.02, 10_000)
    assert empty == ["isolation"] and math.isnan(par["s_mix"]) and math.isnan(par["m_mix"])
    assert _mix_params(np.array([0.98, 0.01, 0.01]), 0.0005, 0.002, 0.01, 0.02, 100)[1] == []   # one pair each
    assert _mix_params(np.array([0.99, 0.0099, 0.0001]), 0.0005, 0.002, 0.01, 0.02, 100)[1] == ["gene flow",
                                                                                                "tail"]


def test_pair_divergence_layer(tmp_path, toy_fastas):
    """The identity layer's pair rows as the catalog's pair-divergence layer: pairs/<D>__<partner>.div.tsv.gz
    (DIV_COLUMNS: D's element, the partner's member element or none, the LTR sites compared and differing, d their
    p-distance, none without a site), each pair once whatever the trios that aligned it, meta engine assembly and
    d_kind p. An earlier assembly layer is replaced; another route's divergence file stays."""
    cat = write_run(tmp_path, [loc(10_001, "R:P* S:P D:P* O:E"), loc(20_001, "R:P S:P* D:P* O:E")],
                    fasta_dir=toy_fastas).load()
    lid = cat.lid.tolist()
    rows = [("self", lid[0], "R,S;D", "D", ".", 40, 0, 0.0, 0.07, None),
            ("pair", lid[0], "R,S;D", "D", "R", 800, 4, 0.005, None, None),
            ("pair", lid[0], "R,S;D", "D", "S", 150, 2, None, None, None),
            ("pair", lid[1], "R,S;D", "D", "S", 0, 0, None, None, None),
            ("pair", lid[0], "S,R;D", "D", "R", 800, 4, 0.005, None, None)]       # another trio, the same pair
    stale = cat.outdir / "pairs" / "D__O.div.tsv.gz"
    other = cat.outdir / "pairs" / "R__S.div.tsv.gz"
    write_divergences(stale, [], {"engine": "assembly", "d_kind": "p"})
    write_divergences(other, [], {"engine": "library", "d_kind": "k2p"})
    paths = write_pair_divergences(cat, rows)
    assert [p.name for p in paths] == ["D__R.div.tsv.gz", "D__S.div.tsv.gz"]
    assert not stale.exists() and not (stale.parent / (stale.name + ".done.json")).exists() and other.exists()
    d, r, s = (cat.g(x) for x in "DRS")
    assert list(read_divergences(paths[0])) == [
        {"eid": cat.eid[0, d], "b_eid": cat.eid[0, r], "n_sites": 800, "n_diff": 4, "d": 0.005}]
    assert list(read_divergences(paths[1])) == [
        {"eid": cat.eid[0, d], "b_eid": None, "n_sites": 150, "n_diff": 2, "d": pytest.approx(2 / 150, rel=1e-5)},
        {"eid": cat.eid[1, d], "b_eid": cat.eid[1, s], "n_sites": 0, "n_diff": 0, "d": None}]
    assert read_calls_meta(paths[0]) == {"engine": "assembly", "d_kind": "p"}
    assert write_pair_divergences(cat, []) == [] and not list((cat.outdir / "pairs").glob("D__*.div.tsv.gz"))
