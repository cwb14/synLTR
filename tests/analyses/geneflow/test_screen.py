"""The gene-flow engine on the library route (`synltr geneflow --screen`): LTRmix's painting and inversion tests on
its synthetic trios, the roles, classes, initial emissions, fit, posteriors, tracts and decoys against LTRmix's own
(stored outputs, tests/golden/ltrmix_paint.npz), the route's rules on constructed inputs (the donor-state check, tract
ends, ties in the painting order, the 1 Mb count of chain blocks), the screen's files and rows, the engine's honest
outputs (a skipped painting, a collapsed fit), the sheet routing, the trio checkpoint's key and the resume after a kill
inside the trio loop, and the overview figure on a screen run."""
from __future__ import annotations

import itertools
import math
import os
import re
import shutil
import signal
import subprocess
import sys
import time
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from synltr import cli
from synltr.analyses.geneflow import inversions as inv
from synltr.analyses.geneflow import library_stats as ls
from synltr.analyses.geneflow import paint as pt
from synltr.analyses.geneflow import run as gf
from synltr.analyses.geneflow.inversions import INV_COLUMNS, chain_inversion_rows, chain_inversions
from synltr.analyses.geneflow.library_stats import LibraryRun, pair_calls, pick_roles, roles
from synltr.analyses.geneflow.route import ROUTE_PROPS
from synltr.analyses.geneflow.trio import Trio
from synltr.base.checkpoint import code_files
from synltr.base.constants import PI_MIN
from synltr.io.samples import read_samples
from synltr.routes.library.run import BLOCK_COLUMNS
from synltr.stats.report import read_stats

GOLDEN = Path(__file__).resolve().parents[2] / "golden" / "ltrmix_paint.npz"
PROPS = ROUTE_PROPS["library"]
MODES = ("ancestry", "recent")
CATS = ("near_donor", "near_sister", "near_both", "donor_only", "sister_only", "shared", "private")  # LTRmix's
CASES = {"admixed": "admixed", "null": "null_screen", "four": "four_screen", "unplaced": "unplaced_screen"}
HMM_FIELDS = ("sister", "donor", "pi", "mean_tract_bp", "loglik", "n_tracts", "n_decoy_tracts", "n_tracts_q05",
              "posterior_mean", "decoy_posterior_mean", "excess_over_decoy", "near_threshold")


@pytest.fixture(scope="module")
def stored():
    with np.load(GOLDEN) as z:
        return {k: z[k] for k in z.files}


def _tsv(path: Path) -> list[dict[str, str]]:
    """A text TSV's rows by its header (a leading '#' dropped)."""
    lines = path.read_text().splitlines()
    head = lines[0].lstrip("#").split("\t")
    return [dict(zip(head, line.split("\t"), strict=True)) for line in lines[1:]]


def _geneflow(sc, src: Path, out: Path, *extra: str) -> Path:
    """`synltr geneflow --screen` on scenario sc's sheet, into a copy of its finished screen `src` (the pairs are
    done: only the analyses run)."""
    shutil.copytree(src, out)
    assert cli.main(["geneflow", "--screen", str(sc.sheet), "-o", str(out), "-t", "2", *extra]) == 0
    return out


@pytest.fixture(scope="module")
def screened(admixed, tmp_path_factory):
    sc, run = admixed
    return sc, _geneflow(sc, run, tmp_path_factory.mktemp("screened") / "run")


@pytest.fixture(scope="module")
def screened_null(null_screen, tmp_path_factory):
    sc, run = null_screen
    return sc, _geneflow(sc, run, tmp_path_factory.mktemp("screened_null") / "run")


def _lib(request, case: str) -> LibraryRun:
    sc, run = request.getfixturevalue(CASES[case])
    return LibraryRun(run, tuple(read_samples(sc.sheet)), 2, {})


def _trios(stored, case: str) -> list[tuple[str, str, str]]:
    return [tuple(stored[f"{case}_t{k}_roles"].tolist()) for k in range(int(stored[f"{case}_ntrios"]))]


def _same(x, y) -> bool:
    return (math.isnan(x) and math.isnan(y)) if (math.isnan(x) or math.isnan(y)) else x == y


# ---------------------------------------------------------------- LTRmix's tests, on the screen


def test_sister_pair_inferred(screened):
    """LTRmix's admixed trio: R and S share the most single-copy k-mers, D is the third genome; the roles row (the
    engine's form, source data) and findings say so."""
    _, out = screened
    rows = [r for r in read_stats(out / "stats.tsv") if r.test == "roles"]
    assert [(r.layer, r.scope, r.statistic, r.value) for r in rows] == [("library", "R,S;D", "source", "data")]
    text = (out / "findings.txt").read_text()
    assert "sisters R and S (the pair sharing the most single-copy k-mers the third genome lacks), donor D" in text


def test_tract_recovered(screened):
    """The recent painting's tracts of R cover more than 60% of each planted tract, and no tract outside them
    reaches q <= 0.05."""
    sc, out = screened
    b = [r for r in _tsv(out / "tracts.R.bed") if r["alphabet"] == "recent"]
    inside = np.zeros(len(b), bool)
    for c, s0, e0 in sc.tracts:
        ov = np.array([r["chrom"] == f"R_chr{c}" and int(r["end"]) > s0 and int(r["start0"]) < e0 for r in b])
        inside |= ov
        covered = sum(min(int(r["end"]), e0) - max(int(r["start0"]), s0) for r, o in zip(b, ov, strict=True) if o)
        assert covered / (e0 - s0) > 0.6, (c, s0, e0, covered)
    assert sum(float(r["q"]) <= 0.05 for r, i in zip(b, inside, strict=True) if not i) == 0


def test_no_confident_tracts_when_not_admixed(screened_null):
    """Without gene flow no tract of either painting reaches q <= 0.05."""
    _, out = screened_null
    assert not [r for r in _tsv(out / "tracts.R.bed") if float(r["q"]) <= 0.05]


def test_inversion_found(screened):
    """The planted inversion (D holds a segment of R's chromosome 2 reversed) is an inverted chain block of at least
    10 anchors."""
    sc, out = screened
    v = _tsv(out / "inversions.tsv")
    assert [r for r in v if r["seqid_x"].endswith(f"chr{sc.inversion[0]}") and int(r["n_loci"]) >= 10]


def test_overview_figure(screened):
    """geneflow --screen draws the overview figure next to its tables (cli, after the command returns): the painting
    panels from the library paintings, and every panel whose data a screen run lacks says why."""
    _, out = screened
    for f in ("overview.png", "overview.pdf", "overview_numbers.tsv"):
        assert (out / f).stat().st_size > 0, f
    text = (out / "overview_numbers.tsv").read_text()
    assert "needs assemblies: run geneflow without --screen" in text


def test_overview_wording(screened, monkeypatch):
    """The screen overview's key and legend say what is true on both rows of a sequence pair, the recipient R's and
    the decoy S's: the lower ticks mark a near-identical copy in the other sister only (S on R's rows, R on S's).
    Every sentence of its legends starts with a capital letter, also when the lead painting was skipped."""
    from synltr.plot import geneflow as gx
    from synltr.plot import style

    _, out = screened

    def legends():
        nums = style.Numbers()
        m = gx._screen_model(gx._Run(out, gx._ScreenCat(out)), nums, None)
        return m, [p for page in gx._screen_legends(m, nums) for p in page]
    m, painted = legends()
    assert [text for kind, _, text in gx._key_items(m) if kind == "tick"] == [
        "near-identical copy in D only", "near-identical copy in the other sister only"]
    ticks = next(t for p in painted for t, _ in p if "Ticks:" in t)
    assert "in D only (above) or in the other sister only (below: S on R's rows, R on S's)" in ticks
    monkeypatch.setattr(gx, "_screen_painting", lambda *a, **k: None)
    _, skipped = legends()
    for para in painted + skipped:
        text = "".join(t for t, _ in para)
        assert all(s[0].isupper() or s[0].isdigit() for s in re.split(r"(?<=\.) +", text) if s), text


# ---------------------------------------------------------------- LTRmix's numbers


@pytest.mark.parametrize("case", list(CASES))
def test_roles_match_ltrmix(case, request, stored):
    """Each trio of the sheet (its combinations in sheet order): the roles LTRmix's Trio inferred (sisters by the most
    single-copy k-mers shared, the recipient by the sign of D_kmer) and its D_kmer, exactly."""
    run = _lib(request, case)
    trios = list(itertools.combinations(run.labels, 3))
    assert len(trios) == int(stored[f"{case}_ntrios"])
    for k, genomes in enumerate(trios):
        r = roles(run, genomes, PROPS)
        assert (r.R, r.S, r.D) == tuple(stored[f"{case}_t{k}_roles"].tolist()), (case, k)
        want = stored[f"{case}_t{k}_dkmer"].tolist()
        got = (r.d.est, r.d.se, r.d.z, r.n_a, r.n_b, r.d.nblocks)
        assert all(_same(float(a), float(b)) for a, b in zip(got, want, strict=True)), (case, k, got, want)
        assert not r.tie and r.source == "data"


def test_roles_on_ties_and_nan(stored):
    """Constructed single-copy tallies through LTRmix's Trio: a tie between the top sharing pairs keeps the first in
    the order {a,b}, {a,c}, {b,c}; a NaN D_kmer gives the second sister, D_kmer 0 the first, a negative one the
    second (D_kmer recomputed from its side). pick_roles gives the same roles and D_kmer, flags the tie and sets the
    source to sheet order where a tie or D_kmer (nan or 0) decided."""
    order = ("a", "b", "c")
    for k in range(int(stored["roles_n"])):
        name = str(stored[f"roles{k}_name"])
        tallies = {}
        for x in order:
            o1, o2 = [y for y in order if y != x]
            t = stored[f"roles{k}_{x}_tally"]
            tallies[x] = {(): t[:, 0], (o1,): t[:, 1], (o2,): t[:, 2], (o1, o2): t[:, 3]}
        bp = int(stored[f"roles{k}_block_bp"])

        def blocks(g, k=k, bp=bp):
            names = [f"{c}:{m // bp}" for c, m in zip(stored[f"roles{k}_{g}_chrom"].tolist(),
                                                     stored[f"roles{k}_{g}_mid"].tolist(), strict=True)]
            return np.unique(np.array(names, object), return_inverse=True)[1].astype(np.int64)
        r = pick_roles(order, tallies, blocks)
        assert (r.R, r.S, r.D) == tuple(stored[f"roles{k}_roles"].tolist()), name
        assert list(r.share.values()) == stored[f"roles{k}_share"].tolist(), name
        got = (r.d.est, r.d.se, r.d.z, r.n_a, r.n_b, r.d.nblocks)
        assert all(_same(float(a), float(b)) for a, b in zip(got, stored[f"roles{k}_dkmer"].tolist(), strict=True)), \
            (name, got)
        assert r.tie == name.startswith("tie_"), name
        assert r.source == ("sheet order" if name.startswith("tie_") or name in ("nan_dkmer", "zero_dkmer")
                            else "data"), name


@pytest.mark.parametrize("case", list(CASES))
def test_classes_match_ltrmix(case, request, stored):
    """Each trio's near threshold (the 5th percentile of S's ortholog divergences to D) and the classes of R's and S's
    records in LTRmix's painting order (sequence name in text order, then midpoint; records whose names hold no
    coordinates on their own sequences), exactly; both alphabets read the same classes."""
    run = _lib(request, case)
    for k, (R, S, D) in enumerate(_trios(stored, case)):
        p = f"{case}_t{k}"
        near = pt.near_threshold(run, S, D)
        assert near == float(stored[f"{p}_thr"]), p
        for X, Y, o, cat in ((R, S, "o", "catX"), (S, R, "oy", "catY")):
            want = [CATS[c] if c >= 0 else "unresolved" for c in stored[f"{p}_{cat}"][stored[f"{p}_{o}"]]]
            for a in MODES:
                idx, cls = pt.classes(run, X, Y, D, alphabet=a, near=near)
                assert idx.tolist() == stored[f"{p}_{o}"].tolist(), (p, X)
                assert [pt.LIBRARY_CLASSES[c] for c in cls] == want, (p, X, a)


@pytest.mark.parametrize("case", list(CASES))
def test_initial_emissions_match_ltrmix(case, request, stored):
    """EM's initial emissions from the decoy's classes (each emitting class + 1; e1 = e0 x 5 for donor-like classes
    and x 0.2 for the others; both normalized), bit for bit."""
    run = _lib(request, case)
    for k, (R, S, D) in enumerate(_trios(stored, case)):
        p = f"{case}_t{k}"
        _, dcls = pt.classes(run, S, R, D, alphabet="ancestry", near=pt.near_threshold(run, S, D))
        for a in MODES:
            e0, e1 = pt.initial_emissions(pt.emissions(dcls, a), a)
            assert e0.tolist() == stored[f"{p}_{a}_e0i"].tolist() and e1.tolist() == stored[f"{p}_{a}_e1i"].tolist()


def _f(spec: str, values) -> list[str]:
    return [format(float(v), spec) for v in values]


def _close(x, y) -> bool:
    x, y = np.asarray(x, float), np.asarray(y, float)
    return x.shape == y.shape and bool(np.allclose(x, y, rtol=1e-9, atol=1e-12, equal_nan=True))


def _check_tracts(tracts, run, X, stored, prefix, q=None, p=None) -> None:
    names = run.seqids[X]
    assert [names[t.seq] for t in tracts] == stored[f"{prefix}_chrom"].tolist(), prefix
    assert [t.start for t in tracts] == stored[f"{prefix}_start"].tolist(), prefix
    assert [t.end for t in tracts] == stored[f"{prefix}_end"].tolist(), prefix
    assert [t.last - t.first + 1 for t in tracts] == stored[f"{prefix}_n"].tolist(), prefix
    mean, total = [t.mean_post for t in tracts], [t.sum_post for t in tracts]
    assert _f(".3f", mean) == _f(".3f", stored[f"{prefix}_mean"]) and _close(mean, stored[f"{prefix}_mean"])
    assert _f(".2f", total) == _f(".2f", stored[f"{prefix}_sum"]) and _close(total, stored[f"{prefix}_sum"])
    assert [min(1000, int(1000 * m)) for m in mean] == [min(1000, int(1000 * m)) for m in stored[f"{prefix}_mean"]]


def _check_skipped(p, stored, case: str, k: int, a: str) -> None:
    """A painting the engine skips (no donor-like or no sister-like element) is one LTRmix fitted from its
    pseudocounts alone: R holds no element of that kind in LTRmix's classes, and its fit collapsed (pi at its floor,
    no tract), so the skip removes no tract LTRmix reported."""
    assert (case, k) == ("four", 3), (case, k, p.skipped)
    alphabet = pt.ALPHABETS[a]
    cat = stored[f"{case}_t{k}_catX"]          # LTRmix's codes are LIBRARY_CLASSES' (its -1: unresolved)
    codes = alphabet.codes[np.where(cat >= 0, cat, pt.UNRESOLVED)]
    donor = [alphabet.emits.index(c) for c in alphabet.emits if c in alphabet.donor]
    assert not np.isin(codes, donor).any() and "no donor-like element" in p.skipped
    pre = f"{case}_t{k}_{a}"
    assert float(stored[f"{pre}_pi"]) == PI_MIN and len(stored[f"{pre}_tr_start"]) == 0


@pytest.mark.parametrize("case", list(CASES))
def test_painting_matches_ltrmix(case, request, stored):
    """Each trio's two paintings against LTRmix's own (its HMM fit, its tract calls and the target-decoy q and p,
    stored in tests/golden/ltrmix_paint.npz) at LTRmix's printed precision: pi, L and the log-likelihood (6
    significant digits), the fitted emissions (3), every element's posterior (4 decimals), each tract's extent and
    element count exactly, its mean (3 decimals; its BED score), sum (2), q and p (3 significant digits), the decoy's
    posteriors and tracts alike, the length-weighted posterior means (6 significant digits), the tracts at q <= 0.05
    and the donor-state check. The engine's shared arithmetic (its geomspace L grid, 3 pi / (1 - pi), exactly rounded
    sums) leaves at most last bits apart: every value also agrees to 1e-9 relative."""
    run = _lib(request, case)
    for k, (R, S, D) in enumerate(_trios(stored, case)):
        for a in MODES:
            pre = f"{case}_t{k}_{a}"
            p = pt.paint(run, Trio(R, S, D, False), PROPS, alphabet=a)
            if p.skipped is not None:      # the engine's honest output (the only case: four, trio D,E;S)
                _check_skipped(p, stored, case, k, a)
                continue
            assert p.idx.tolist() == stored[f"{case}_t{k}_o"].tolist(), pre
            fit = p.fit
            for got, key in ((fit.pi, "pi"), (fit.L, "L"), (fit.loglik, "ll")):
                assert _f(".6g", [got]) == _f(".6g", [stored[f"{pre}_{key}"]]) and _close(got, stored[f"{pre}_{key}"])
            for e in ("e0", "e1"):
                assert _f(".3g", getattr(fit, e)) == _f(".3g", stored[f"{pre}_{e}"]), (pre, e)
                assert _close(getattr(fit, e), stored[f"{pre}_{e}"]), (pre, e)
            assert _f(".4f", p.post) == _f(".4f", stored[f"{pre}_post"]) and _close(p.post, stored[f"{pre}_post"])
            assert p.donor_state == bool(stored[f"{pre}_donor_state"]), pre
            _check_tracts(p.tracts, run, R, stored, f"{pre}_tr")
            assert _f(".3g", p.q) == _f(".3g", stored[f"{pre}_q"]) and _close(p.q, stored[f"{pre}_q"]), pre
            assert _f(".3g", p.p) == _f(".3g", stored[f"{pre}_p"]) and _close(p.p, stored[f"{pre}_p"]), pre
            assert sum(q <= 0.05 for q in p.q) == int(stored[f"{pre}_nq"]), pre
            assert p.decoy_idx.tolist() == stored[f"{case}_t{k}_oy"].tolist(), pre
            assert _f(".4f", p.decoy_post) == _f(".4f", stored[f"{pre}_py"]) and _close(p.decoy_post,
                                                                                        stored[f"{pre}_py"])
            _check_tracts(p.decoy_tracts, run, S, stored, f"{pre}_dec")
            assert p.scale == int(stored[f"{case}_t{k}_nR"]) / int(stored[f"{case}_t{k}_nS"]), pre
            for got, key in ((pt.posterior_mean(p, run), "frac_post"),
                             (pt.posterior_mean(p, run, decoy=True), "frac_dec")):
                assert _f(".6g", [got]) == _f(".6g", [stored[f"{pre}_{key}"]]) and _close(got, stored[f"{pre}_{key}"])


def test_inversions_match_ltrmix(stored, tmp_path):
    """LTRmix's Trio.inversions on constructed chain blocks (both signs; a chromosome pair whose two signs hold
    equally many anchors: + is main; a pair whose main orientation is -; pairs under 50 anchors skipped, exactly 50
    kept; equal n_anchor across pairs, in pair order): the same blocks in the same order, in the engine's
    inversions.tsv columns."""
    rows = stored["inv_blocks"].tolist()
    path = tmp_path / "blocks.tsv"
    path.write_text("\t".join(BLOCK_COLUMNS) + "\n" + "".join("\t".join(r) + "\n" for r in rows))
    seqids = {g: tuple(sorted({r[1] for r in rows if r[0] == g} | {r[5] for r in rows if r[4] == g})) for g in "RSD"}
    blocks = chain_inversions(path, seqids)
    want = [line.split("\t") for line in str(stored["inv_tsv"]).splitlines()[1:]]
    got = [[b.x, seqids[b.x][b.seq_x], str(b.start_x), str(b.end_x), b.y, seqids[b.y][b.seq_y], str(b.start_y),
            str(b.end_y), str(b.n_loci), str(b.end_x - b.start_x)] for b in blocks]
    assert got == want and len(got) == 8
    assert [b.dominant for b in blocks if (b.x, b.y) == ("R", "D") and b.seq_x == seqids["R"].index("chr1")] == [-1]


def test_chain_block_1mb_by_anchor_distance(tmp_path):
    """A chain block's extent runs from its first to its last anchor's midpoint, so it counts in n_blocks_1Mb when
    that distance, end - start, reaches 1 Mb, as LTRmix counted: 999,999 bp does not, 1,000,000 bp does.
    bp_inverted_x sums end - start + 1, and the assembly route counts a block of loci by that closed extent."""
    blocks_tsv = [("R", "chr1", 1, 5_000_000, "S", "chr1", 1, 5_000_000, 1, 100),
                  ("R", "chr1", 6_000_001, 7_000_000, "S", "chr1", 6_000_001, 7_000_000, -1, 10),
                  ("R", "chr1", 8_000_001, 9_000_001, "S", "chr1", 8_000_001, 9_000_001, -1, 10)]
    path = tmp_path / "blocks.tsv"
    path.write_text("\t".join(BLOCK_COLUMNS) + "\n" + "".join("\t".join(map(str, r)) + "\n" for r in blocks_tsv))
    seqids = {"R": ("chr1",), "S": ("chr1",)}
    blocks = chain_inversions(path, seqids)
    assert sorted(b.end_x - b.start_x for b in blocks) == [999_999, 1_000_000]
    rows, _ = chain_inversion_rows(blocks, path, ["R", "S"], seqids)
    assert {r.statistic: r.value for r in rows if r.scope == "R~S"} == {
        "n_blocks": 2, "n_blocks_1Mb": 1, "bp_inverted_x": 1_000_000 + 1_000_001}
    loci = [replace(blocks[0], start_x=1, end_x=end) for end in (999_999, 1_000_000)]
    assert {r.statistic: r.value for r in inv._count_rows("pav", "R~S", loci, 0, "x")}["n_blocks_1Mb"] == 1


# ---------------------------------------------------------------- the screen's files and rows


def test_screen_files(screened):
    """The screen's paintings in the engine's file names: posterior.R.tsv (every element of R per trio: eid, seqid,
    start, end, sister, donor, the class, one posterior per alphabet, the inverted-block flag), tracts.R.bed and
    decoy_tracts.S.bed holding both alphabets (the alphabet in the BED name and in its own column),
    decoy_posterior.S.tsv; stats rows hmm_ancestry and hmm_recent with LTRmix's fields; inversions.tsv in the
    engine's columns, "." where a chain block has no such value."""
    sc, out = screened
    post = _tsv(out / "posterior.R.tsv")
    assert tuple(post[0]) == pt.library_posterior_columns(MODES) == (
        "eid", "seqid", "start", "end", "sister", "donor", "class", "post_ancestry", "post_recent", "inv")
    run = LibraryRun(out, tuple(read_samples(sc.sheet)), 2, {})
    eids = [e for e in run.rec("R").lib.elements.column("eid").tolist() if e is not None]
    assert sorted(r["eid"] for r in post) == sorted(eids) and {(r["sister"], r["donor"]) for r in post} == {("S", "D")}
    assert {r["class"] for r in post} <= set(pt.LIBRARY_CLASSES)
    for name, sister in (("tracts.R.bed", "S"), ("decoy_tracts.S.bed", "R")):
        rows = _tsv(out / name)
        assert tuple(rows[0]) == tuple(c.lstrip("#") for c in pt.LIBRARY_TRACT_COLUMNS), name
        assert {r["alphabet"] for r in rows} == set(MODES), name
        assert all(r["name"] == f"D-{r['alphabet']}" and r["sister"] == sister for r in rows), name
    assert {r["q"] for r in _tsv(out / "decoy_tracts.S.bed")} == {"."}
    assert (out / "decoy_posterior.S.tsv").is_file()
    stats = read_stats(out / "stats.tsv")
    for a, names in (("ancestry", ("near_donor", "donor_only", "near_sister", "sister_only")),
                     ("recent", ("near_donor", "near_sister"))):
        mine = [r for r in stats if r.test == f"hmm_{a}"]
        assert {(r.layer, r.scope, r.stratum) for r in mine} == {("library", "R,S;D", "R")}
        assert [r.statistic for r in mine] == [*HMM_FIELDS, *(f"e1_{c}" for c in names),
                                               *(f"e0_{c}" for c in names), "donor_state"]
    inv = _tsv(out / "inversions.tsv")
    assert tuple(inv[0]) == INV_COLUMNS
    assert {(r["bp_left_x"], r["bp_right_y"], r["n_strand"], r["n_discordant"]) for r in inv} == {(".",) * 4}
    assert {r.scope for r in stats if r.test == "inversions"} == {"R~S", "R~D", "S~D"}


def test_posteriors_read_back_exactly(screened):
    """posterior.R.tsv and the tract files write every float exactly: read back, each equals the painting's own
    value, so LTRmix's 4-decimal elements table and its tract BEDs are rebuilt exactly."""
    sc, out = screened
    run = LibraryRun(out, tuple(read_samples(sc.sheet)), 2, {})
    eid = run.rec("R").lib.elements.column("eid")
    rows = {r["eid"]: r for r in _tsv(out / "posterior.R.tsv")}
    tracts = _tsv(out / "tracts.R.bed")
    for a in MODES:
        p = pt.paint(run, Trio("R", "S", "D", False), PROPS, alphabet=a)
        got = [float(rows[eid[i]][f"post_{a}"]) for i in p.idx.tolist() if eid[i] is not None]
        assert got == [v for i, v in zip(p.idx.tolist(), p.post.tolist(), strict=True) if eid[i] is not None]
        mine = [r for r in tracts if r["alphabet"] == a]
        assert [(float(r["mean_post"]), float(r["sum_post"]), float(r["q"]), float(r["p_decoy"])) for r in mine] == [
            (t.mean_post, t.sum_post, q, pv) for t, q, pv in zip(p.tracts, p.q, p.p, strict=True)]


def test_assembly_analyses_named_as_skipped(screened):
    """One findings line names the analyses a screen does not run (they need observed empty sites from assemblies)
    and points to geneflow without --screen; the run records its route."""
    _, out = screened
    text = (out / "findings.txt").read_text()
    assert text.count(gf.SCREEN_SKIPS) == 1
    assert "`synltr geneflow` without --screen" in gf.SCREEN_SKIPS and "(E)" in gf.SCREEN_SKIPS
    assert not [r for r in read_stats(out / "stats.tsv") if r.test in ("D_pav", "relations", "mechanisms")]


def test_collapsed_fit_prints_na(screened_null):
    """LTRmix's null trio: both fits collapse (pi at its floor, no tract): L (mean_tract_bp), the emissions and the
    donor-state check are NA (state 1 holds no data), pi, the log-likelihood and the posterior means stay, and
    findings say the fit collapsed."""
    _, out = screened_null
    stats = read_stats(out / "stats.tsv")
    for a in MODES:
        row = {r.statistic: r for r in stats if r.test == f"hmm_{a}"}
        assert row["pi"].value == PI_MIN and row["n_tracts"].value == 0
        assert all(math.isnan(r.value) and r.note.startswith("NA: the fit collapsed") for s, r in row.items()
                   if s == "mean_tract_bp" or s.startswith(("e0_", "e1_")))
        assert math.isnan(row["donor_state"].value)
        assert math.isfinite(row["loglik"].value) and math.isfinite(row["posterior_mean"].value)
    assert "the fit collapsed (pi at its floor" in (out / "findings.txt").read_text()


def test_painting_skipped_without_both_kinds(admixed, monkeypatch):
    """With no element near either genome (a near threshold of 0) the recent painting has no donor-like and no
    sister-like element: it is skipped with its reason (rows NA, one WARNING), while the ancestry painting runs."""
    sc, run_dir = admixed
    run = LibraryRun(run_dir, tuple(read_samples(sc.sheet)), 2, {})
    monkeypatch.setattr(pt, "near_threshold", lambda *a: 0.0)
    trio = Trio("R", "S", "D", False)
    recent = pt.paint(run, trio, PROPS, alphabet="recent")
    assert recent.skipped is not None and recent.fit is None
    assert "no donor-like element (near_donor) and no sister-like element (near_sister)" in recent.skipped
    rows, f = pt.library_rows(recent, run)
    by = {r.statistic: r for r in rows}
    assert by["sister"].value == "S" and by["near_threshold"].value == 0.0
    assert all(math.isnan(r.value) for s, r in by.items() if s not in ("sister", "donor", "near_threshold"))
    assert [t for _, t in f.items() if t.startswith("WARNING: Recent shared ancestry of R")]
    assert pt.paint(run, trio, PROPS, alphabet="ancestry").skipped is None


def test_donor_state_needs_near_sister_down():
    """The library route's donor-state check has two halves: the fitted state 1 emits near_donor above state 0 and
    near_sister below it. On constructed ancestry fits (near_donor, donor_only, near_sister, sister_only) a fit where
    near_donor goes up but near_sister does not go down fails the check, as one where near_donor does not go up."""
    def ok(e1) -> bool:
        fit = pt.HmmFit(pi=0.1, L=1e5, e0=np.full(4, 0.25), e1=np.array(e1), post=np.zeros(0), loglik=0.0)
        return pt._donor_state(PROPS.donor_state, "ancestry", fit)
    assert ok([0.4, 0.2, 0.1, 0.3])
    assert not ok([0.4, 0.05, 0.3, 0.25]) and not ok([0.4, 0.1, 0.25, 0.25])
    assert not ok([0.25, 0.35, 0.1, 0.3]) and not ok([0.2, 0.45, 0.1, 0.25])


def test_tract_end_on_nested_ends():
    """A tract whose last element ends before an earlier one (a short element inside a long one): on the library
    route it ends at its last element's end (tract_end "last_end", LTRmix's tracts), on the assembly route at its
    elements' largest end ("largest_end": the interval holds every locus)."""
    seq = np.zeros(4, np.int64)
    start, end = np.array([1_000, 2_000, 5_000, 10_000]), np.array([1_500, 9_000, 6_200, 10_500])
    post = np.array([0.99, 0.99, 0.99, 0.1])
    for props, stop in ((PROPS, 6_200), (ROUTE_PROPS["assembly"], 9_000)):
        [t] = pt.call_tracts(seq, start, end, post, 0.1, tract_end=props.tract_end)
        assert (t.start, t.end, t.first, t.last) == (1_000, stop, 0, 2), props.tract_end


def test_painting_order_keeps_library_order_on_a_midpoint_tie(admixed):
    """Records on one sequence with the same midpoint keep their library order in the painting (a stable sort, as
    LTRmix's lexsort), whatever their starts and ends: 40 of R's records given one midpoint, their starts and ends
    in another order than the library's."""
    sc, run_dir = admixed
    sheet = tuple(read_samples(sc.sheet))
    rec = LibraryRun(run_dir, sheet, 2, {}).rec("R")
    chrom = rec.chrom.astype(str)
    tie = np.flatnonzero(chrom == chrom[0])[:40]
    half = np.random.default_rng(1).permutation(len(tie)) * 10 + 100
    mid, start, end = rec.mid.copy(), rec.start.copy(), rec.end.copy()
    m = int(np.median(mid[tie]))
    mid[tie], start[tie], end[tie] = m, m - half, m + half
    run = LibraryRun(run_dir, sheet, 2, {"R": replace(rec, mid=mid, start=start, end=end)})
    idx, _ = pt.classes(run, "R", "S", "D", alphabet="recent", near=pt.near_threshold(run, "S", "D"))
    assert idx.tolist() == sorted(range(rec.n), key=lambda i: (chrom[i], int(mid[i])))  # Python's sort is stable
    at = [idx.tolist().index(i) for i in tie.tolist()]
    assert at == list(range(at[0], at[0] + len(tie)))


# ---------------------------------------------------------------- routing


def _copy(admixed, tmp_path: Path, sheet_text: str) -> tuple[Path, Path]:
    """A copy of the admixed trio's finished screen and a sheet beside it."""
    _, run = admixed
    out = tmp_path / "run"
    shutil.copytree(run, out)
    sheet = tmp_path / "sheet.tsv"
    sheet.write_text(sheet_text)
    return sheet, out


def _rows(sc, fasta: dict[str, str] | None = None, library: bool = True) -> str:
    cols = ["label"] + (["fasta"] if fasta is not None else []) + (["library"] if library else [])
    lines = ["\t".join(cols)]
    for g in "RSD":
        lines.append("\t".join([g] + ([fasta.get(g, ".")] if fasta is not None else [])
                               + ([sc.libraries[g]] if library else [])))
    return "\n".join(lines) + "\n"


def test_sheet_without_fastas_runs_the_screen(admixed, tmp_path, capsys, monkeypatch):
    """A sheet without genome FASTAs runs the screen, detected before calls.genotype.read_sheet, and says so: one
    INFO line and one findings line."""
    sc, _ = admixed
    sheet, out = _copy(admixed, tmp_path, _rows(sc))

    def refuse(*args, **kwargs):
        raise AssertionError("read_sheet reached")
    monkeypatch.setattr(gf.genotype, "read_sheet", refuse)
    assert cli.main(["geneflow", str(sheet), "-o", str(out), "-t", "2"]) == 0
    err = capsys.readouterr().err
    assert sum("holds no genome FASTA, so geneflow runs the library screen" in ln for ln in err.splitlines()) == 1
    assert (out / "findings.txt").read_text().count(gf.SCREEN_ROUTED) == 1
    assert (out / "tracts.R.bed").is_file()


def test_sheet_with_some_fastas_is_a_usage_error(admixed, tmp_path, capsys):
    """FASTAs on some rows only: exit 2, one error line naming the rows without one and both remedies, nothing
    written."""
    sc, _ = admixed
    fa = sc.libraries["R"].split(",")[0]
    sheet, out = _copy(admixed, tmp_path, _rows(sc, {"R": fa}))
    before = sorted(p.name for p in out.iterdir())
    assert cli.main(["geneflow", str(sheet), "-o", str(out)]) == 2
    err = capsys.readouterr().err.strip().splitlines()
    assert len(err) == 1 and err[0].startswith("error: ") and "row(s) S, D have no genome FASTA" in err[0]
    assert "give every row a fasta" in err[0] and "--screen" in err[0]
    assert sorted(p.name for p in out.iterdir()) == before


def test_screen_refuses_identity_and_paf(admixed, tmp_path, capsys):
    """--screen with --identity or --paf (each realigns assemblies), and a sheet without FASTAs with either: exit 2,
    one error line."""
    sc, _ = admixed
    sheet, out = _copy(admixed, tmp_path, _rows(sc))
    for args in (["--screen", "--identity"], ["--screen", "--paf", str(tmp_path)], ["--identity"]):
        assert cli.main(["geneflow", str(sheet), "-o", str(out), *args]) == 2, args
        err = capsys.readouterr().err.strip().splitlines()
        assert len(err) == 1 and "realigns genome assemblies" in err[0], (args, err)


def test_screen_on_a_sheet_with_fastas_reads_the_libraries(admixed, tmp_path, capsys):
    """--screen on a sheet whose every row has a FASTA reads the libraries only (the FASTAs are never opened: here
    they are not even FASTA files); a row without a library is the screen's usage error."""
    sc, _ = admixed
    junk = tmp_path / "not_a_genome.fa"
    junk.write_text("this is not a FASTA\n")
    sheet, out = _copy(admixed, tmp_path, _rows(sc, dict.fromkeys("RSD", str(junk))))
    assert cli.main(["geneflow", "--screen", str(sheet), "-o", str(out), "-t", "2"]) == 0
    assert (out / "stats.tsv").is_file()
    capsys.readouterr()
    bare = tmp_path / "bare.tsv"
    bare.write_text("label\tfasta\tlibrary\n" + "".join(f"{g}\t{junk}\t{sc.libraries[g] if g != 'D' else '.'}\n"
                                                        for g in "RSD"))
    assert cli.main(["geneflow", "--screen", str(bare), "-o", str(tmp_path / "x")]) == 2
    assert "D has no library" in capsys.readouterr().err


def test_failed_screen_rerun_leaves_no_record_of_the_earlier_run(admixed, tmp_path, monkeypatch, capsys):
    """In a finished screen's folder, a rerun stopped by a usage error (a row without a library) leaves its
    manifest.json and overview as they were; a rerun whose analyses fail (exit 3) leaves neither, nor stats.tsv and
    findings.txt, so nothing there describes results the folder no longer holds."""
    sc, _ = admixed
    sheet, out = _copy(admixed, tmp_path, _rows(sc))
    argv = ["geneflow", "--screen", str(sheet), "-o", str(out), "-t", "2"]
    assert cli.main(argv) == 0
    record = ("manifest.json", "overview.png", "overview.pdf", "overview_numbers.tsv")
    before = {name: (out / name).read_bytes() for name in record}
    junk = tmp_path / "not_a_genome.fa"  # --screen never opens a FASTA
    junk.write_text("this is not a FASTA\n")
    bare = tmp_path / "bare.tsv"
    bare.write_text("label\tfasta\tlibrary\n" + "".join(f"{g}\t{junk}\t{sc.libraries[g] if g != 'D' else '.'}\n"
                                                        for g in "RSD"))
    assert cli.main(["geneflow", "--screen", str(bare), "-o", str(out)]) == 2
    assert "D has no library" in capsys.readouterr().err
    assert {name: (out / name).read_bytes() for name in record} == before

    def boom(*args, **kwargs):
        raise ValueError("a planted analysis error")

    monkeypatch.setattr(gf, "compute", boom)
    assert cli.main(argv) == 3 and (out / "blocks.tsv").is_file()
    assert not [name for name in (*record, "stats.tsv", "findings.txt") if (out / name).exists()]


def test_screen_trio_option(screened, tmp_path):
    """--trio on the screen: the roles as given (source --trio), checked against the sheet like the assembly route
    (check_trio_option)."""
    sc, src = screened
    out = _geneflow(sc, src, tmp_path / "run", "--trio", "R,S,D")
    rows = [r for r in read_stats(out / "stats.tsv") if r.test == "roles"]
    assert [(r.scope, r.value) for r in rows] == [("R,S;D", "--trio")]
    assert cli.main(["geneflow", "--screen", str(sc.sheet), "-o", str(out), "--trio", "R,S,Z"]) == 2


# ---------------------------------------------------------------- resume after a kill inside the trio loop


def _outputs(run: Path) -> dict[str, bytes]:
    """The analyses' outputs (findings.txt without its header line, which names the folder)."""
    out = {}
    for name in ("stats.tsv", "inversions.tsv", "findings.txt", *(p.name for p in run.glob("*tracts.*.bed")),
                 *(p.name for p in run.glob("*posterior.*.tsv"))):
        data = (run / name).read_bytes()
        out[name] = data.split(b"\n", 1)[1] if name == "findings.txt" else data
    return out


def _trip(run: LibraryRun, R: str, S: str, D: str, n: int = 8) -> np.ndarray:
    """The trio's first n ortholog triples as the statistics take them (rows: D's, R's and S's record index)."""
    rd, rs, sd = (pair_calls(run, x, y).orth for x, y in ((R, D), (R, S), (S, D)))
    tri = np.flatnonzero((rd >= 0) & (rs >= 0))
    tri = tri[sd[rs[tri]] == rd[tri]][:n]
    return np.stack([rd[tri], tri, rs[tri]], 1)


@pytest.mark.parametrize("module", ["io/library.py", "routes/library/__init__.py", "routes/library/clock.py"])
def test_trio_checkpoint_keyed_on_the_code_of_its_counts(four_screen, tmp_path, monkeypatch, module):
    """Each trio's checkpointed alignments are keyed on the code their counts depend on, as the screen's pair files
    are: the record reader (io/library.py: the records' sequences), the records' LTR bounds
    (routes/library/__init__.py) and the alignments (routes/library/clock.py). With the code unchanged a rerun
    reuses every trio; after an edit to any of these modules (its source digest changed) every trio is aligned
    again."""
    sc, src = four_screen
    shutil.copytree(src, tmp_path / "run")
    shutil.rmtree(tmp_path / "run" / "work" / "geneflow", ignore_errors=True)
    sheet = tuple(read_samples(sc.sheet))
    trios = list(itertools.combinations([s.label for s in sheet], 3))

    def reused() -> list[bool]:
        run = LibraryRun(tmp_path / "run", sheet, 2, {})
        return [ls.triples(run, *t, _trip(run, *t))[2] for t in trios]
    assert len(trios) == 4 and reused() == [False] * 4
    assert reused() == [True] * 4
    real = ls.code_fingerprint
    monkeypatch.setattr(ls, "code_fingerprint", lambda unit: real(unit) + ("+" if module in code_files(unit) else ""))
    assert reused() == [False] * 4


@pytest.mark.slow
def test_resume_after_a_kill_in_the_trio_loop(four_screen, tmp_path):
    """A screen of four genomes (four trios) killed (SIGKILL) inside the trio loop, once the first trio's alignments
    were checkpointed: the rerun reuses every finished trio's alignments (their files untouched) and aligns only the
    rest, and its outputs are byte-identical to an uninterrupted run's."""
    sc, src = four_screen
    whole = _geneflow(sc, src, tmp_path / "whole")
    killed = tmp_path / "killed"
    shutil.copytree(src, killed)
    env = {**os.environ, "PYTHONPATH": os.pathsep.join(sys.path)}
    cmd = [sys.executable, "-c", "import sys; from synltr.cli import main; sys.exit(main(sys.argv[1:]))",
           "geneflow", "--screen", str(sc.sheet), "-o", str(killed), "-t", "2"]
    proc = subprocess.Popen(cmd, env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    triples = killed / "work" / "geneflow" / "triples"
    deadline = time.monotonic() + 600
    while not list(triples.glob("*.done.json")) and proc.poll() is None and time.monotonic() < deadline:
        time.sleep(0.05)
    proc.send_signal(signal.SIGKILL)
    proc.wait()
    done = {p.name: p.stat().st_mtime_ns for p in triples.glob("*.npz")
            if (triples / f"{p.name}.done.json").is_file()}
    assert 1 <= len(done) < 4 and not (killed / "stats.tsv").exists(), done
    assert cli.main(["geneflow", "--screen", str(sc.sheet), "-o", str(killed), "-t", "2"]) == 0
    after = {p.name: p.stat().st_mtime_ns for p in triples.glob("*.npz")}
    assert len(after) == 4 and all(after[n] == t for n, t in done.items())
    assert _outputs(killed) == _outputs(whole)
