"""The genotype command (calls.genotype): every pair of a sheet of assemblies genotyped by the assembly route, then
the loci, tables, VCF and tracks; --pair and --paf; a run killed between pairs resumes."""
from __future__ import annotations

import gzip
import os
import random
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

import builders
import pytest
from markers import needs_lastz

from ltrtrace import cli

OUTPUTS = ("loci.tsv", "pav_matrix.tsv", "qc.tsv", "pav.vcf.gz", "pav.vcf.gz.tbi", "pav.unplaced.tsv",
           *(f"{g}.pav.{x}" for g in "ABC" for x in ("gff3", "bed")))
PAIRS = ("A__B.calls.tsv.gz", "A__C.calls.tsv.gz")


@pytest.fixture(scope="module")
def sheet(tmp_path_factory):
    """Genomes A (annotated: the E, P, S and host_E toy cases), B (the cases' B) and C (a copy of A), and their
    samples sheet."""
    d = tmp_path_factory.mktemp("sheet")
    pair = builders.build_pair(random.Random(5), ["E", "P", "S", "host_E"])
    builders.write_fasta(d / "A.fa", {"chr1": pair.a})
    builders.write_fasta(d / "B.fa", {"chr1": pair.b})
    builders.write_fasta(d / "C.fa", {"chr1": pair.a})
    builders.write_gff3(d / "A.gff3", pair.elements)
    (d / "samples.tsv").write_text("label\tfasta\tannotation\nA\tA.fa\tA.gff3\nB\tB.fa\t.\nC\tC.fa\t.\n")
    return d / "samples.tsv", pair


def _files(out: Path) -> dict[str, bytes]:
    """The deliverables of a genotype run: the pair call files and the outputs, by name."""
    return {**{f"pairs/{p}": (out / "pairs" / p).read_bytes() for p in PAIRS},
            **{f: (out / f).read_bytes() for f in OUTPUTS}}


@pytest.fixture(scope="module")
def ref_run(sheet, tmp_path_factory):
    """An uninterrupted run with 2 threads."""
    if shutil.which("lastz") is None:
        pytest.skip("lastz (the window aligner) not on PATH")
    out = tmp_path_factory.mktemp("ref") / "out"
    assert cli.main(["genotype", str(sheet[0]), "-o", str(out), "-t", "2"]) == 0
    return out


@needs_lastz
def test_genotype_writes_the_calls_and_the_loci(sheet, ref_run):
    _, pair = sheet
    assert sorted(p.name for p in (ref_run / "pairs").glob("*.gz")) == list(PAIRS)
    for name in OUTPUTS:
        assert (ref_run / name).is_file(), name
    names = {e.eid: e.name for e in pair.elements}
    calls = builders.read_calls(ref_run / "pairs" / "A__B.calls.tsv.gz")
    want = {names[eid]: state for eid, state in pair.expect.items()}
    got = {names[r["eid"].replace("A|", "toy|")]: r["state"] for r in calls.values()}
    assert got == want
    assert {r["state"] for r in builders.read_calls(ref_run / "pairs" / "A__C.calls.tsv.gz").values()} == {"P"}
    loci = (ref_run / "loci.tsv").read_text().splitlines()
    assert len(loci) == 1 + len(pair.elements)  # one locus per element: B and C hold no other annotation
    with gzip.open(ref_run / "pav.vcf.gz", "rt") as fh:
        assert [line.split("\t")[9:] for line in fh if not line.startswith("#")]  # genotypes on A's coordinates


@needs_lastz
def test_threads_do_not_change_the_outputs(sheet, ref_run, tmp_path):
    out = tmp_path / "out"
    assert cli.main(["genotype", str(sheet[0]), "-o", str(out), "-t", "1"]) == 0
    assert _files(out) == _files(ref_run)


@needs_lastz
def test_pair_option_genotypes_one_pair_and_stops(sheet, ref_run, tmp_path):
    out = tmp_path / "out"
    assert cli.main(["genotype", str(sheet[0]), "-o", str(out), "--pair", "A,C"]) == 0
    assert sorted(p.name for p in (out / "pairs").glob("*.gz")) == ["A__C.calls.tsv.gz"]
    assert not (out / "loci.tsv").exists()
    assert (out / "pairs" / "A__C.calls.tsv.gz").read_bytes() == (ref_run / "pairs" / "A__C.calls.tsv.gz").read_bytes()
    assert cli.main(["genotype", str(sheet[0]), "-o", str(out)]) == 0  # the rest of the run, then the loci
    assert _files(out) == _files(ref_run)


@needs_lastz
def test_paf_option_anchors_with_the_given_alignments(sheet, ref_run, tmp_path):
    """Whole-genome PAFs (here A's sequence aligned to C, an exact copy) anchor the pairs instead of flank pieces;
    a pair without a PAF fails naming the files it looked for."""
    paf = tmp_path / "paf"
    paf.mkdir()
    n = len(sheet[1].a)
    (paf / "A__C.paf").write_text(f"chr1\t{n}\t0\t{n}\t+\tchr1\t{n}\t0\t{n}\t{n}\t{n}\t60\tcg:Z:{n}=\n")
    out = tmp_path / "out"
    assert cli.main(["genotype", str(sheet[0]), "-o", str(out), "--pair", "A,C", "--paf", str(paf)]) == 0
    assert {r["state"] for r in builders.read_calls(out / "pairs" / "A__C.calls.tsv.gz").values()} == {"P"}
    assert (out / "work" / "anchors" / "A__C.paf.npz").is_file()
    assert cli.main(["genotype", str(sheet[0]), "-o", str(out), "--paf", str(paf)]) == 1  # no PAF for A vs B


@pytest.mark.parametrize("argv, msg", [
    (["--pair", "A,Z"], "'Z' is not a label"),
    (["--pair", "B,A"], "B has no annotation"),
    (["--pair", "A"], "expected two different labels"),
    (["--pair", "A,A"], "expected two different labels"),
    (["--paf", "/no/such/dir"], "is not a directory"),
])
def test_bad_options_are_usage_errors(sheet, tmp_path, capsys, argv, msg):
    assert cli.main(["genotype", str(sheet[0]), "-o", str(tmp_path / "out"), *argv]) == 2
    assert msg in capsys.readouterr().err


def test_missing_sheet_is_a_usage_error(tmp_path, capsys):
    assert cli.main(["genotype", str(tmp_path / "none.tsv"), "-o", str(tmp_path)]) == 2
    assert "samples sheet not found" in capsys.readouterr().err


def _reads_sheet(d: Path, extra: str = "") -> Path:
    """A sheet of reference R (one element) and read set S, whose BAM is an empty placeholder file."""
    builders.write_fasta(d / "R.fa", {"chr1": "ACGT" * 100})
    (d / "r.bam").write_text("")
    (d / "R.gff3").write_text("##gff-version 3\nchr1\tt\tLTR_retrotransposon\t101\t300\t.\t+\t.\tID=e1\n"
                              "chr1\tt\tlong_terminal_repeat\t101\t150\t.\t+\t.\tParent=e1\n"
                              "chr1\tt\tlong_terminal_repeat\t251\t300\t.\t+\t.\tParent=e1\n")
    (d / "s.tsv").write_text("label\tfasta\tannotation\trole\treads\nR\tR.fa\tR.gff3\treference\t.\n"
                             "S\t.\t.\t.\tr.bam\n" + extra)
    return d / "s.tsv"


@pytest.mark.parametrize("argv, msg", [
    (["--pair", "R,S"], "--pair and --paf genotype assemblies"),
    (["--paf", "."], "--pair and --paf genotype assemblies"),
    (["--ploidy", "0"], "--ploidy: must be at least 1"),
    (["--min-mapq", "-1"], "--min-mapq: must be at least 0"),
    (["--ploidy", "two"], "--ploidy: invalid integer"),
])
def test_bad_options_on_read_sets_are_usage_errors(tmp_path, capsys, argv, msg):
    assert cli.main(["genotype", str(_reads_sheet(tmp_path)), "-o", str(tmp_path / "out"), *argv]) == 2
    assert msg in capsys.readouterr().err
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize("argv", [["--ploidy", "4"], ["--min-mapq", "10"]])
def test_read_options_on_assemblies_are_usage_errors(sheet, tmp_path, capsys, argv):
    assert cli.main(["genotype", str(sheet[0]), "-o", str(tmp_path / "out"), *argv]) == 2
    assert "--ploidy and --min-mapq apply to read sets" in capsys.readouterr().err
    assert not (tmp_path / "out").exists()


def test_reads_sheet_row_without_reads_is_a_usage_error(tmp_path, capsys):
    """A sheet of read sets holds one reference and read sets: another assembly there would be silently unused."""
    sheet = _reads_sheet(tmp_path, "B\tR.fa\t.\t.\t.\n")
    assert cli.main(["genotype", str(sheet), "-o", str(tmp_path / "out")]) == 2
    err = capsys.readouterr().err.splitlines()
    assert len(err) == 1 and "row(s) B have no reads" in err[0]


def test_reads_sheet_goes_to_the_reads_route(tmp_path, capsys):
    """A sheet with read sets is genotyped from the reads: the reference and its annotation load, then the
    placeholder BAM fails as a BAM (one error line naming it), before any output."""
    assert cli.main(["genotype", str(_reads_sheet(tmp_path)), "-o", str(tmp_path / "out")]) == 2
    err = capsys.readouterr().err.splitlines()
    assert len(err) == 1 and "r.bam" in err[0] and "assemblies" not in err[0]
    assert not (tmp_path / "out").exists()


def test_rows_without_fasta_are_a_usage_error(tmp_path, capsys):
    """A sheet with rows that carry a library but no FASTA (a geneflow --screen sheet) cannot be genotyped: one
    error line naming those rows, exit 2, before any work."""
    d = tmp_path
    builders.write_fasta(d / "A.fa", {"chr1": "ACGT" * 100})
    (d / "lib.fa").write_text(">x\nACGT\n")
    (d / "A.gff3").write_text("##gff-version 3\n")
    (d / "s.tsv").write_text("label\tfasta\tannotation\tlibrary\nA\tA.fa\tA.gff3\tlib.fa\nB\t.\t.\tlib.fa\n"
                             "C\t.\t.\tlib.fa\n")
    assert cli.main(["genotype", str(d / "s.tsv"), "-o", str(d / "out")]) == 2
    err = capsys.readouterr().err.splitlines()
    assert len(err) == 1 and "every row needs a fasta; none for B, C" in err[0]
    assert not (d / "out").exists()


# a genotype run whose second pair (A vs C) stalls in its anchoring: argv = the command line
_SLOW_SECOND_PAIR = r"""
import sys, time
from ltrtrace import cli
from ltrtrace.routes.assembly import anchor

real = anchor.FlankMm2Route.raw_blocks


def raw_blocks(self, ctx):
    if ctx.b.label == "C":
        time.sleep(120)
    return real(self, ctx)


anchor.FlankMm2Route.raw_blocks = raw_blocks
sys.exit(cli.main(sys.argv[1:]))
"""


@needs_lastz
def test_resume_after_kill(sheet, ref_run, tmp_path, capsys):
    """A run killed (SIGKILL) after the first pair's done marker, while the second pair runs: the rerun skips the
    first pair, genotypes the second, and writes the same files, byte for byte, as an uninterrupted run."""
    out = tmp_path / "out"
    first = out / "pairs" / "A__B.calls.tsv.gz"
    done = first.with_name(first.name + ".done.json")
    script = tmp_path / "slow.py"
    script.write_text(_SLOW_SECOND_PAIR)
    proc = subprocess.Popen([sys.executable, str(script), "genotype", str(sheet[0]), "-o", str(out), "-t", "2"],
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    try:
        t0 = time.monotonic()
        while not done.exists():
            assert proc.poll() is None and time.monotonic() - t0 < 120
            time.sleep(0.05)
        os.kill(proc.pid, signal.SIGKILL)
        proc.wait()
    finally:
        if proc.poll() is None:
            proc.kill()
            proc.wait()
    assert not (out / "pairs" / "A__C.calls.tsv.gz").exists() and not (out / "loci.tsv").exists()
    before = (first.stat().st_mtime_ns, done.read_bytes())
    capsys.readouterr()
    assert cli.main(["genotype", str(sheet[0]), "-o", str(out), "-t", "2"]) == 0
    err = capsys.readouterr().err
    assert "A vs B: calls done, skip" in err and "A vs C: calls done" not in err
    assert err.count("A: elements done, skip") == 1  # read once per run: the outputs reuse what the route read
    assert (first.stat().st_mtime_ns, done.read_bytes()) == before  # skipped, not rewritten
    assert _files(out) == _files(ref_run)
