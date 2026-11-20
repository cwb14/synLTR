"""The geneflow command and its analyses (analyses.geneflow.run): every pair genotyped first, then the analyses of
the finished run, rerun every time; stats.tsv and findings.txt written last; the identity layer, its alignment
checkpoint and the pair-divergence layer; usage and analysis errors. Toy runs from catalog_toy.py; worlds with real
sequence from builders.py."""
from __future__ import annotations

import gzip
import hashlib
import json
import logging
import random
import re
import subprocess
import sys
from pathlib import Path

import builders
import pytest
from analyses.geneflow.catalog_toy import GENOMES, TRACTS, loc, simulate_run, write_run
from markers import needs_lastz

from ltrtrace import __version__, cli
from ltrtrace.analyses.geneflow import identity, run
from ltrtrace.analyses.geneflow.catalog import AnalysisError
from ltrtrace.analyses.geneflow.identity import IDENTITY_COLUMNS
from ltrtrace.calls.records import write_divergences
from ltrtrace.io.elements import Element, make_eid
from ltrtrace.stats.report import SECTIONS, STATS_COLUMNS, Findings, StatRow, read_stats


def _world(d: Path) -> Path:
    """Genomes A (annotated: the E, P, S and host_E toy cases), B (the cases' B) and C (a copy of A), and their
    samples sheet."""
    pair = builders.build_pair(random.Random(5), ["E", "P", "S", "host_E"])
    builders.write_fasta(d / "A.fa", {"chr1": pair.a})
    builders.write_fasta(d / "B.fa", {"chr1": pair.b})
    builders.write_fasta(d / "C.fa", {"chr1": pair.a})
    builders.write_gff3(d / "A.gff3", pair.elements)
    (d / "samples.tsv").write_text("label\tfasta\tannotation\nA\tA.fa\tA.gff3\nB\tB.fa\t.\nC\tC.fa\t.\n")
    return d / "samples.tsv"


def _mini(d: Path) -> Path:
    """A (annotated, one element) and B, unrelated 6-kb genomes, and their samples sheet."""
    rng = random.Random(1)
    builders.write_fasta(d / "A.fa", {"chr1": builders.rand_seq(6000, rng)})
    builders.write_fasta(d / "B.fa", {"chr1": builders.rand_seq(6000, rng)})
    e = Element(genome="A", eid=make_eid("A", "chr1", 1001, 2000), name="e", seqid="chr1", start=1001, end=2000,
                strand="+", ltr_left=(1001, 1100), ltr_right=(1901, 2000), tsd=builders.T, k=0.0, family=None,
                superfamily=None, clade=None, confidence=None)
    builders.write_gff3(d / "A.gff3", [e])
    sheet = d / "samples.tsv"
    sheet.write_text("label\tfasta\tannotation\nA\tA.fa\tA.gff3\nB\tB.fa\t.\n")
    return sheet


def _analyses(run_, **kw):
    return run.analyses(run_.sheet, run_.outdir, 1, **kw)


@needs_lastz
def test_geneflow_genotypes_then_analyses(tmp_path, capsys):
    """`ltrtrace geneflow` genotypes every pair, then runs the analyses; a rerun skips the finished pairs and runs
    the analyses again (they are never checkpointed)."""
    sheet, out = _world(tmp_path), tmp_path / "out"
    argv = ["geneflow", str(sheet), "-o", str(out), "-t", "2"]
    assert cli.main(argv) == 0
    assert all((out / f).is_file() for f in ("loci.tsv", "stats.tsv", "findings.txt", "inversions.tsv",
                                             "mechanisms.tsv"))
    assert "Trio tests" in (out / "findings.txt").read_text().splitlines()
    first = (out / "stats.tsv").read_bytes()
    capsys.readouterr()
    assert cli.main(argv) == 0
    err = capsys.readouterr().err
    assert "A vs B: calls done, skip" in err and "A vs C: calls done, skip" in err
    assert "geneflow: analyses of" in err and (out / "stats.tsv").read_bytes() == first


def test_analyses_on_a_toy_run(tmp_path, toy_fastas):
    run_ = simulate_run(tmp_path, seed=1, tracts=TRACTS, fasta_dir=toy_fastas)
    stats, findings = _analyses(run_)
    assert (stats, findings) == (run_.outdir / "stats.tsv", run_.outdir / "findings.txt")
    assert {"sharing", "D_pav", "paint"} <= {r.test for r in read_stats(stats)}
    assert stats.read_text().splitlines()[0].split("\t") == list(STATS_COLUMNS)
    lines = findings.read_text().splitlines()
    assert lines[0] == f"ltrtrace {__version__} findings for the run in {run_.outdir}"
    heads = [ln for ln in lines if ln in SECTIONS]
    assert heads == sorted(heads, key=SECTIONS.index) and (run_.outdir / "tracts.R.bed").is_file()
    assert "run" not in {r.layer for r in read_stats(stats)}               # no sentinel rows


def test_trio_validation(tmp_path, toy_fastas, capsys):
    """A --trio that names an unknown genome, one genome twice or an outgroup is a usage error found before any
    genotyping: one error line, exit 2."""
    run_ = simulate_run(tmp_path, seed=1, fasta_dir=toy_fastas)
    out = tmp_path / "out"
    for bad in ("R,R,D", "R,S,X", "R,S,O"):
        assert cli.main(["geneflow", str(run_.sheet), "-o", str(out), "--trio", bad]) == 2
        err = capsys.readouterr().err
        assert err.strip().splitlines()[-1].startswith("error: --trio") and err.count("\n") == 1
        assert "Traceback" not in err and not out.exists()


def test_missing_loci(tmp_path, toy_fastas):
    run_ = simulate_run(tmp_path, seed=1, fasta_dir=toy_fastas)
    (run_.outdir / "loci.tsv").unlink()
    with pytest.raises(AnalysisError, match="ltrtrace genotype"):
        _analyses(run_)


def test_help_is_short(capsys):
    assert cli.main(["geneflow", "--help"]) == 0
    text = capsys.readouterr().out
    assert len(text.splitlines()) < 40 and all(opt in text for opt in ("--trio", "--identity", "--paf"))


def _copies_world(d: Path) -> tuple[Path, Path]:
    """Genome A (annotated: the E, P, S and host_E toy cases) and B and C, exact copies of it; their samples sheet,
    and a folder of whole-genome PAFs, A's sequence against each copy end to end."""
    pair = builders.build_pair(random.Random(5), ["E", "P", "S", "host_E"])
    for g in "ABC":
        builders.write_fasta(d / f"{g}.fa", {"chr1": pair.a})
    builders.write_gff3(d / "A.gff3", pair.elements)
    (d / "samples.tsv").write_text("label\tfasta\tannotation\nA\tA.fa\tA.gff3\nB\tB.fa\t.\nC\tC.fa\t.\n")
    paf = d / "paf"
    paf.mkdir()
    n = len(pair.a)
    for b in "BC":
        (paf / f"A__{b}.paf").write_text(f"chr1\t{n}\t0\t{n}\t+\tchr1\t{n}\t0\t{n}\t{n}\t{n}\t60\tcg:Z:{n}=\n")
    return d / "samples.tsv", paf


@needs_lastz
def test_geneflow_paf_reuses_the_pairs_genotype_made(tmp_path, capsys):
    """geneflow --paf DIR anchors every pair with the given whole-genome alignments, as genotype --paf does, so it
    reuses the pairs a genotype --paf run made (a pair's key holds its route), then runs the analyses."""
    sheet, paf = _copies_world(tmp_path)
    out = tmp_path / "out"
    assert cli.main(["genotype", str(sheet), "-o", str(out), "--paf", str(paf)]) == 0
    pairs = {p.name: p.read_bytes() for p in sorted((out / "pairs").glob("*.calls.tsv.gz"))}
    assert sorted(pairs) == ["A__B.calls.tsv.gz", "A__C.calls.tsv.gz"]
    alignments = {name: hashlib.sha256((paf / name).read_bytes()).hexdigest() for name in ("A__B.paf", "A__C.paf")}
    assert json.loads((out / "manifest.json").read_text())["params"]["paf_checksums"] == alignments
    capsys.readouterr()
    assert cli.main(["geneflow", str(sheet), "-o", str(out), "--paf", str(paf)]) == 0
    err = capsys.readouterr().err
    assert "A vs B: calls done, skip" in err and "A vs C: calls done, skip" in err
    assert {p.name: p.read_bytes() for p in sorted((out / "pairs").glob("*.calls.tsv.gz"))} == pairs
    assert (out / "stats.tsv").is_file() and "Trio tests" in (out / "findings.txt").read_text().splitlines()
    params = json.loads((out / "manifest.json").read_text())["params"]
    assert (params["command"], params["paf"], params["paf_checksums"]) == ("geneflow", str(paf), alignments)


_NO_SCIPY = "import sys; sys.modules['scipy'] = None; sys.modules['matplotlib'] = None; "


def test_cli_imports_without_scipy():
    """The command line and the geneflow command module import with scipy and matplotlib blocked: the heavy
    imports wait for a command to run."""
    for mod in ("ltrtrace.cli", "ltrtrace.analyses.geneflow.run"):
        r = subprocess.run([sys.executable, "-c", _NO_SCIPY + f"import {mod}"], capture_output=True, text=True,
                           check=False)
        assert r.returncode == 0, r.stderr


def test_shared_seqids(tmp_path, toy_fastas):
    a = simulate_run(tmp_path / "a", seed=1, tracts=TRACTS, fasta_dir=toy_fastas)
    b = simulate_run(tmp_path / "b", seed=1, tracts=TRACTS, prefix=True, fasta_dir=toy_fastas)
    for run_ in (a, b):
        _analyses(run_)

    def norm(run_):
        return [re.sub(r"\b[RSDO]_(chr\d)\b", r"\1", ln) for ln in (run_.outdir / "stats.tsv").read_text().splitlines()]
    assert norm(a) == norm(b)


def test_two_genomes(tmp_path, toy_fastas):
    run_ = write_run(tmp_path, [loc(1000 * i, "R:P* S:E O:E") for i in range(1, 20)],
                     genomes=[GENOMES[0], GENOMES[1], GENOMES[3]], fasta_dir=toy_fastas)
    _analyses(run_)
    assert "Trio tests skipped: fewer than three ingroup genomes" in (run_.outdir / "findings.txt").read_text()
    assert (run_.outdir / "mechanisms.tsv").is_file()


def test_analyses_write_their_files(tmp_path, toy_fastas):
    """The analyses write stats.tsv, findings.txt, inversions.tsv, mechanisms.tsv and the paintings beside the
    genotyping outputs, and under work/ only the record of the painting files they wrote."""
    run_ = simulate_run(tmp_path, seed=1, tracts=TRACTS, fasta_dir=toy_fastas)
    _analyses(run_)
    painted = ["tracts.R.bed", "posterior.R.tsv", "decoy_tracts.S.bed", "decoy_posterior.S.tsv"]
    assert {p.name for p in run_.outdir.iterdir() if p.is_file()} == (
        {"loci.tsv", "qc.tsv", "stats.tsv", "findings.txt", "inversions.tsv", "mechanisms.tsv"} | set(painted))
    assert json.loads((run_.outdir / run.PAINTINGS).read_text()) == painted
    assert sorted(p.name for p in (run_.outdir / "work").iterdir()) == ["elements", "geneflow"]
    assert sorted(p.name for p in (run_.outdir / "work" / "geneflow").iterdir()) == ["paintings.json"]


def test_without_identity_the_identity_layer_goes(tmp_path, toy_fastas, caplog):
    """Analyses without --identity delete an earlier run's identity.tsv and its assembly pair-divergence files,
    with one WARNING that --identity keeps them, and say in Dating what needs the identity layer; a divergence
    file of another route stays."""
    run_ = simulate_run(tmp_path, seed=1, fasta_dir=toy_fastas)
    ident = run_.outdir / "identity.tsv"
    ident.write_text("\t".join(IDENTITY_COLUMNS) + "\n")      # as a --identity run with no rows leaves it
    mine, other = run_.outdir / "pairs" / "D__R.div.tsv.gz", run_.outdir / "pairs" / "R__S.div.tsv.gz"
    write_divergences(mine, [], {"engine": "assembly", "d_kind": "p"})
    write_divergences(other, [], {"engine": "library", "d_kind": "k2p"})
    with caplog.at_level(logging.INFO, logger="ltrtrace"):
        _analyses(run_)
    warn = [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING and "identity.tsv" in r.getMessage()]
    assert len(warn) == 1 and warn[0].startswith("geneflow: deleting identity.tsv") and "with --identity to keep" in warn[0]
    assert not ident.exists() and not mine.exists() and other.exists()
    text = (run_.outdir / "findings.txt").read_text()
    assert "Counts of shared insertions younger than the split need the sequence-identity layer" in text
    assert "identity" not in {r.layer for r in read_stats(run_.outdir / "stats.tsv")}


def _sequence_run(tmp_path: Path, toy_fastas: Path, n: int = 4):
    """R, S and D holding the same four 2,800-bp elements (S's with 0.2% substitutions) in real sequence, the
    outgroup O empty: a run the identity layer can align. n: the first n of the four loci are written."""
    rng = random.Random(9)
    els = [builders.make_element(rng, 400, 2000)[0] for _ in range(4)]
    flanks = [builders.rand_seq(2000, rng) for _ in range(5)]
    copies = {"R": els, "S": [builders.mutate(e, 0.002, rng) for e in els], "D": els}
    fastas = {}
    for g, cs in copies.items():
        fastas[g] = tmp_path / f"{g}.fa"
        fastas[g].write_text(">chr1\n" + "".join(f + c for f, c in zip(flanks, cs)) + flanks[4] + "\n")
    loci = [loc(2001 + 4800 * i, "R:P* S:P* D:P* O:E", length=2800, ltr=400) for i in range(n)]
    return write_run(tmp_path / "w", loci, seqids=("chr1",), fastas=fastas, fasta_dir=toy_fastas)


def _div(path: Path) -> list[str]:
    return gzip.decompress(path.read_bytes()).decode().splitlines()


@needs_lastz
def test_identity_through_the_analyses(tmp_path, toy_fastas):
    """With --identity: identity.tsv, the identity rows in stats.tsv, and the pair-divergence layer of each donor
    and partner (pairs/<D>__<Y>.div.tsv.gz, one row per aligned pair); without it, all three go."""
    run_ = _sequence_run(tmp_path, toy_fastas)
    _analyses(run_, identity=True)
    rows = (run_.outdir / "identity.tsv").read_text().splitlines()
    assert rows[0].split("\t") == list(IDENTITY_COLUMNS) and sum(r.startswith("self\t") for r in rows) == 12
    assert any(r.layer == "identity" for r in read_stats(run_.outdir / "stats.tsv"))
    pairs = [r.split("\t") for r in rows if r.startswith("pair\t")]
    div = {y: _div(run_.outdir / "pairs" / f"D__{y}.div.tsv.gz") for y in "RS"}
    for y, lines in div.items():
        assert lines[0].split("\t") == ["eid", "b_eid", "n_sites", "n_diff", "d"] and len(lines) == 5
        got = [line.split("\t")[2:4] for line in lines[1:]]
        assert got == [[p[5], p[6]] for p in pairs if p[4] == y]                  # the identity layer's counts
    assert div["R"][1].split("\t")[:2] == ["D|chr1:2001-4800", "R|chr1:2001-4800"]
    _analyses(run_)
    assert not (run_.outdir / "identity.tsv").exists() and not list((run_.outdir / "pairs").glob("*.div.tsv.gz"))


@needs_lastz
def test_analyses_rerun_identity_cached(tmp_path, toy_fastas, monkeypatch, caplog):
    """The identity alignment is the one checkpointed analysis: a second run reuses it (its outputs byte-identical),
    while every other analysis reruns; another LASTZ version realigns."""
    run_ = _sequence_run(tmp_path, toy_fastas)
    _analyses(run_, identity=True)
    files = ("identity.tsv", "stats.tsv", "findings.txt", "pairs/D__R.div.tsv.gz", "pairs/D__S.div.tsv.gz")
    before = {f: (run_.outdir / f).read_bytes() for f in files}

    def refuse(*args, **kwargs):
        raise AssertionError("the identity alignment ran again")
    real = identity._align
    monkeypatch.setattr(identity, "_align", refuse)
    with caplog.at_level(logging.INFO, logger="ltrtrace"):
        _analyses(run_, identity=True)
    assert "identity R,S;D: alignment done, skip" in caplog.text
    assert {f: (run_.outdir / f).read_bytes() for f in files} == before
    monkeypatch.setattr(identity, "_align", real)
    monkeypatch.setattr(identity, "lastz_version", lambda: "lastz 0.0 (another build)")
    caplog.clear()
    with caplog.at_level(logging.INFO, logger="ltrtrace"):
        _analyses(run_, identity=True)
    assert "alignment done, skip" not in caplog.text and "identity R,S;D: aligning" in caplog.text
    assert (run_.outdir / "identity.tsv").read_bytes() == before["identity.tsv"]


@needs_lastz
@pytest.mark.parametrize("change", ["fasta", "jobs"])
def test_identity_realigns_when_its_inputs_change(tmp_path, toy_fastas, caplog, change):
    """The alignment checkpoint's key holds the FASTAs read and the alignment jobs: S's FASTA changed at the same
    coordinates (one base of its first flank), or the run re-merged with a locus fewer (another job list), each
    realign."""
    run_ = _sequence_run(tmp_path, toy_fastas)
    _analyses(run_, identity=True)
    if change == "fasta":
        fa = tmp_path / "S.fa"
        text = fa.read_text()
        i = text.index("\n") + 10
        fa.write_text(text[:i] + ("C" if text[i] == "A" else "A") + text[i + 1:])
    else:
        _sequence_run(tmp_path, toy_fastas, n=3)
    caplog.clear()
    with caplog.at_level(logging.INFO, logger="ltrtrace"):
        _analyses(run_, identity=True)
    assert "alignment done, skip" not in caplog.text and "identity R,S;D: aligning" in caplog.text
    pairs = [line for line in (run_.outdir / "identity.tsv").read_text().splitlines() if line.startswith("pair\t")]
    assert len(pairs) == 2 * (4 if change == "fasta" else 3)


def test_reads_pair_files_without_done_markers(tmp_path, toy_fastas):
    """The analyses read a finished run as data: the toy's pair files carry no done markers."""
    run_ = write_run(tmp_path, [loc(10_001, "R:P* S:E D:P* O:E"), loc(20_001, "R:E S:P* D:P* O:E")],
                     fasta_dir=toy_fastas)
    assert not list((run_.outdir / "pairs").glob("*.done.json"))
    _analyses(run_)
    assert (run_.outdir / "stats.tsv").is_file()


def test_write_results(tmp_path, monkeypatch):
    """stats.tsv: the header row, then the rows. findings.txt: the version and run directory, each section with
    findings under its name in SECTIONS order (a line's tabs kept, its newlines spaces), then Warnings repeating
    every WARNING after its section's name ("- none" without one)."""
    out = tmp_path / "run"
    out.mkdir()
    f = Findings()
    f.add("Trio tests", "t1")
    f.warn("Run", "w1")
    f.add("Run", "r1")
    f.add("Painting", "p1\twith a tab\nand a newline")
    rows = [StatRow("pav", "sharing", "R,S,D", "all", "n_R+S", 3, n=7), StatRow("pav", "paint", "R,S;D", "all", "pi", 0.25)]
    stats, findings = run.write_results(out, rows, f)
    assert (stats, findings) == (out / "stats.tsv", out / "findings.txt")
    assert stats.read_text().splitlines() == ["\t".join(STATS_COLUMNS), "pav\tsharing\tR,S,D\tall\tn_R+S\t3\t.\t.\t7\t.\t.",
                                             "pav\tpaint\tR,S;D\tall\tpi\t0.25\t.\t.\t.\t.\t."]
    assert findings.read_text().splitlines() == [
        f"ltrtrace {__version__} findings for the run in {out}", "",
        "Run", "- WARNING: w1", "- r1", "",
        "Trio tests", "- t1", "",
        "Painting", "- p1\twith a tab and a newline", "",
        "Warnings", "- Run: WARNING: w1"]
    monkeypatch.chdir(tmp_path)                                  # a relative run directory: named absolute
    g = Findings()
    g.add("Run", "r1")
    lines = run.write_results(Path("run"), [], g)[1].read_text().splitlines()
    assert lines[0].endswith(f" {out}") and lines[-2:] == ["Warnings", "- none"]
    with pytest.raises(ValueError, match="not a findings section"):
        h = Findings()
        h.add("Nowhere", "x")
        run.write_results(out, [], h)


def test_geneflow_trio_fails_fast(tmp_path, capsys):
    """geneflow checks --trio against the sheet before any genotyping: one error line, nothing written; a --trio
    that is not three labels is a usage error too."""
    sheet = _mini(tmp_path)
    base = ["geneflow", str(sheet), "-o", str(tmp_path / "out")]
    for bad, message in (("A,B,Z", "'Z' is not a label"), ("A,A,B", "three different genomes")):
        assert cli.main([*base, "--trio", bad]) == 2
        err = capsys.readouterr().err
        assert err.startswith(f"error: --trio {bad}: ") and message in err and err.count("\n") == 1, err
    assert cli.main([*base, "--trio", "A,B"]) == 2
    assert "expected three labels" in capsys.readouterr().err and not (tmp_path / "out").exists()


@needs_lastz
@pytest.mark.parametrize("planted", ["run_summary", "write_paintings"])
def test_analysis_error_exits_3(tmp_path, capsys, monkeypatch, planted):
    """The genotyping outputs are complete when the analyses start, and they are kept. An analysis failure logs
    one line naming the kept outputs, ends stderr with exactly one error line (no traceback without -v) and exits 3;
    stats.tsv and findings.txt, written last, are not there."""
    def boom(*args, **kwargs):
        raise ValueError("a planted analysis error")
    monkeypatch.setattr(run, planted, boom)
    out = tmp_path / "out"
    assert cli.main(["geneflow", str(_mini(tmp_path)), "-o", str(out), "-t", "1"]) == 3
    err = capsys.readouterr().err
    lines = err.strip().splitlines()
    assert lines[-1] == "error: analyses failed: ValueError: a planted analysis error"
    assert sum(ln.startswith("error:") for ln in lines) == 1 and "Traceback" not in err
    assert sum("the genotyping outputs in" in ln and "are complete" in ln for ln in lines) == 1
    assert (out / "loci.tsv").is_file() and not (out / "stats.tsv").exists() and not (out / "findings.txt").exists()


def test_geneflow_has_no_variant_options(tmp_path, toy_fastas, capsys):
    """One painting model and one block rule: geneflow takes neither --discovery, --emission nor --block-size (a
    usage error, exit 2), and every painting's emission row is "em"."""
    assert cli.main(["geneflow", "--help"]) == 0
    helptext = capsys.readouterr().out
    run_ = simulate_run(tmp_path, seed=1, tracts=TRACTS, fasta_dir=toy_fastas)
    for option in (["--discovery", "any"], ["--emission", "decoy"], ["--block-size", "1000000"]):
        assert option[0] not in helptext
        assert cli.main(["geneflow", str(run_.sheet), "-o", str(run_.outdir), *option]) == 2
    _analyses(run_)
    assert {r.value for r in read_stats(run_.outdir / "stats.tsv")
            if r.test == "paint" and r.statistic == "emission"} == {"em"}


@needs_lastz
def test_one_block_size_reaches_every_jackknife(tmp_path, toy_fastas, monkeypatch):
    """The analyses hand one block length, the route's rule (here 1/100 of the 4-Mb toy genome, 40 kb), to the
    trio tests and to the identity layer; stats.tsv says it in its block_bp row."""
    run_ = simulate_run(tmp_path, seed=1, tracts=TRACTS, fasta_dir=toy_fastas)
    seen: dict[str, set] = {}
    for name in ("trio_tests", "identity_layer"):
        real = getattr(run, name)

        def spy(*args, _real=real, _name=name, **kwargs):
            seen.setdefault(_name, set()).add(args[2])
            return ([], Findings(), []) if _name == "identity_layer" else _real(*args, **kwargs)   # no alignment
        monkeypatch.setattr(run, name, spy)
    _analyses(run_, identity=True)
    assert seen == {"trio_tests": {40_000}, "identity_layer": {40_000}}
    (bp,) = [r for r in read_stats(run_.outdir / "stats.tsv") if (r.test, r.statistic) == ("jackknife", "block_bp")]
    assert (bp.layer, bp.scope, bp.stratum, bp.value, bp.loci, bp.note) == ("pav", "all", "all", 40_000, "none", "rule")
