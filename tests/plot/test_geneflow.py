"""The gene-flow overview figure (plot.geneflow) and `ltrtrace plot`: files, legend pages, numbers, fallbacks (two
genomes, many scaffolds, no painting), the example's structural checks on real sequence, and the command line
(`plot`, the figure drawn after `geneflow`, a figure failure only a warning)."""
from __future__ import annotations

import datetime
import fnmatch
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import PIL.Image
import pytest
from analyses.geneflow.catalog_toy import GENOMES, TRACTS, loc, read_tsv, simulate_run, write_run
from builders import write_fasta
from markers import needs_lastz

from ltrtrace import cli
from ltrtrace.analyses.geneflow import run as geneflow_run
from ltrtrace.analyses.geneflow.catalog import load_catalog
from ltrtrace.analyses.geneflow.paint import PAINTING_FILES
from ltrtrace.base.constants import MAX_ROWS, NUMBERS_COLUMNS
from ltrtrace.base.log import UsageError
from ltrtrace.calls.manifest import write_manifest
from ltrtrace.io.samples import read_samples
from ltrtrace.io.seqio import Fasta, revcomp
from ltrtrace.plot import geneflow as ov
from ltrtrace.plot import style
from ltrtrace.plot.geneflow import _ends, _ident, _site, overview
from ltrtrace.stats.report import read_stats

FIGURE = ("overview.png", "overview.pdf", "overview_numbers.tsv")


def _record(outdir, sheet, route="assembly"):
    """manifest.json as `ltrtrace geneflow` writes it when it finishes: plot reads the samples sheet and the route
    from it."""
    now = datetime.datetime.now(datetime.UTC)
    write_manifest(Path(outdir) / "manifest.json", ["ltrtrace", "geneflow", str(sheet), "-o", str(outdir)],
                   {"command": "geneflow", "sheet": os.path.abspath(sheet), "outdir": os.path.abspath(outdir)},
                   read_samples(sheet), {}, now, now, route=route)


def _analyses(run, trio=None):
    geneflow_run.analyses(run.sheet, run.outdir, 1, trio=trio)
    _record(run.outdir, run.sheet)


def _analyzed(path, toy_fastas, **kw):
    run = simulate_run(path, seed=1, fasta_dir=toy_fastas, **kw)
    _analyses(run)
    return run


def _numbers(path):
    return {r["key"]: r for r in read_tsv(path)}


_PAGES: list[str] = []      # the legend pages the generator drew (style.text_page's paragraphs), the latest last


@pytest.fixture(autouse=True)
def _record_legend_pages(monkeypatch):
    """Record each legend page's text as the generator writes it, so the legend tests read the
    generator's text where pdftotext is missing (and the PDF's, by pdftotext, where it is installed)."""
    real = style.text_page
    _PAGES.clear()

    def spy(paras, family=None):
        _PAGES.append(" ".join(text for para in paras for text, _ in para))
        return real(paras, family)
    monkeypatch.setattr(style, "text_page", spy)


def _pdf_pages(pdf: Path, first: int, last: int) -> str:
    return subprocess.run(["pdftotext", "-f", str(first), "-l", str(last), str(pdf), "-"], capture_output=True,
                          text=True, check=True).stdout


def _legend_text(pdf: Path) -> str:
    """The short and long legend pages of the figure just drawn: from the PDF by pdftotext when installed, else the
    generator's text."""
    return _pdf_pages(pdf, 2, 3) if shutil.which("pdftotext") else " ".join(_PAGES[-2:])


def test_overview_files(tmp_path, toy_fastas):
    run = _analyzed(tmp_path, toy_fastas, tracts=TRACTS, inversion=("S", "chr2", 300_001, 700_000))
    png, pdf, nums = overview(run.outdir)
    with PIL.Image.open(png) as im:                                           # closed: -W error
        assert min(im.info["dpi"]) >= 299.5                                    # Pillow reports 299.9994 for 300
    assert len(re.findall(rb"/Type\s*/Page[^s]", pdf.read_bytes())) == 3
    rows = read_tsv(nums)
    assert list(rows[0]) == list(NUMBERS_COLUMNS) and len(rows) >= 10
    kinds = ("genome index:", "genome sequence:", "display:")   # sources that are no run file
    assert all(r["source"].startswith(kinds) or (run.outdir / r["source"].split(":")[0]).is_file() for r in rows)
    assert "D_pav" in _numbers(nums)


def test_overview_legend_text(tmp_path, toy_fastas):
    run = _analyzed(tmp_path, toy_fastas, tracts=TRACTS)
    _, pdf, nums = overview(run.outdir)
    text = _legend_text(pdf)
    assert _numbers(nums)["D_pav"]["value"] in text and "Analysis:" in text


def test_every_number_is_in_the_legends(tmp_path, toy_fastas):
    """Each numbers row's value is written, as given, on the legend pages (pages 2-3; whitespace as one space)."""
    run = simulate_run(tmp_path, seed=1, fasta_dir=toy_fastas, tracts=TRACTS)
    _analyses(run)
    _check_numbers_in_legends(*overview(run.outdir)[1:])


def _check_numbers_in_legends(pdf: Path, nums: Path) -> None:
    text = " ".join(_legend_text(pdf).split())
    missing = [(r["key"], r["value"]) for r in read_tsv(nums) if " ".join(r["value"].split()) not in text]
    assert not missing


def test_overview_two_genomes(tmp_path, toy_fastas):
    run = write_run(tmp_path, [loc(1000 * i, "R:P* S:E O:E") for i in range(1, 20)],
                    genomes=[GENOMES[0], GENOMES[1], GENOMES[3]], fasta_dir=toy_fastas)
    _analyses(run)
    png, pdf, nums = overview(run.outdir)
    assert png.is_file() and _numbers(nums)["no_trio"]["value"] == "fewer than three ingroup genomes"
    _check_numbers_in_legends(pdf, nums)


def test_overview_many_scaffolds(tmp_path, toy_fastas):
    seqids = ("chr1",) + tuple(f"scaf{i:03d}" for i in range(300))
    loci = ([loc(1000 * i, "R:P* S:E D:P O:E") for i in range(1, 60)]
            + [loc(70_000 + 1000 * i, "R:P* S:P D:E O:E") for i in range(60)]
            + [loc(1001, "R:P* S:E D:P O:E", seqid=f"scaf{i:03d}") for i in range(300)])
    run = write_run(tmp_path, loci, seqids=seqids, seq_len=5_000, seq_lens={"chr1": 200_000}, fasta_dir=toy_fastas)
    _analyses(run)
    num = _numbers(overview(run.outdir)[2])
    assert int(num["sequences_drawn"]["value"]) <= MAX_ROWS and int(num["sequences_not_drawn"]["value"]) >= 277


def test_overview_deterministic(tmp_path, toy_fastas):
    run = _analyzed(tmp_path, toy_fastas, tracts=TRACTS)
    first = overview(run.outdir)[2].read_text()
    assert overview(run.outdir)[2].read_text() == first


def test_overview_numbers_layout(tmp_path, toy_fastas):
    """overview_numbers.tsv: the columns key, value and source, one row per key, each source a run file (or one of the
    figure's own kinds), written atomically beside the figure."""
    run = _analyzed(tmp_path, toy_fastas, tracts=TRACTS)
    overview(run.outdir)
    lines = (run.outdir / "overview_numbers.tsv").read_text().splitlines()
    assert lines[0].split("\t") == ["key", "value", "source"]
    rows = [x.split("\t") for x in lines[1:]]
    assert all(len(r) == 3 and all(r) for r in rows) and len({r[0] for r in rows}) == len(rows)
    assert not list(run.outdir.glob("*.tmp*")) and not list(run.outdir.glob(".overview*"))


def test_plot_help_is_short(capsys):
    assert cli.main(["plot", "--help"]) == 0
    out = capsys.readouterr().out
    assert len(out.splitlines()) < 40 and "RUN" in out


# ---------------------------------------------------------------- the plot stage: options, skips and errors


def test_overview_trio_option(tmp_path, toy_fastas):
    run = _analyzed(tmp_path, toy_fastas, tracts=TRACTS)
    with pytest.raises(ValueError, match="not a trio of this run"):
        overview(run.outdir, trio=("S", "R", "D"))
    nums = overview(run.outdir, trio=("R", "S", "D"))[2]
    assert _numbers(nums)["D_pav"]["source"] == "stats.tsv:D_pav/R,S;D/all/D"


def test_plot_trio_option_names_the_lead_trio(tmp_path, toy_fastas):
    """`ltrtrace plot RUN --trio R,S,D` redraws a run of several trios around the trio named: its numbers come from
    that trio's rows and the headline says the user named it; without --trio the strongest trio leads."""
    run = _four_ingroup_run(tmp_path, toy_fastas)
    _analyses(run)
    trios = [r.scope for r in read_stats(run.outdir / "stats.tsv") if r.test == "roles" and r.statistic == "source"]
    nums, pdf = run.outdir / "overview_numbers.tsv", run.outdir / "overview.pdf"
    assert len(trios) == 4 and cli.main(["plot", str(run.outdir)]) == 0
    lead = _numbers(nums)["D_pav"]["source"].split("/")[1]
    assert lead in trios and "This is the strongest of the run's 4 trios" in _short_legend(pdf)
    other = next(t for t in trios if t != lead)
    assert cli.main(["plot", str(run.outdir), "--trio", other.replace(";", ",")]) == 0
    assert _numbers(nums)["D_pav"]["source"] == f"stats.tsv:D_pav/{other}/all/D"
    assert "This trio was named by the user" in _short_legend(pdf)


@pytest.mark.parametrize("trio, message", [
    ("R,S,Q", "'Q' is not a label of {sheet} (R, S, D, O)"),
    ("R,S,O", "O has role outgroup in {sheet}"),
    ("R,R,S", "the recipient, sister and donor must be three different genomes"),
    ("S,R,D", "not a trio this run tested (its trios: R,S,D)"),
], ids=["unknown-label", "outgroup", "repeated", "untested"])
def test_plot_trio_is_checked_as_geneflow_checks_it(tmp_path, toy_fastas, capsys, trio, message):
    """plot --trio names three ingroup labels of the run's samples sheet, as geneflow's --trio does, and a trio the
    run tested: anything else is a usage error (exit 2, one `error:` line naming the sheet or the run's trios) and
    leaves the figure as it was (here: none)."""
    run = _analyzed(tmp_path, toy_fastas, tracts=TRACTS)  # geneflow with roles from the data: one trio, R,S;D
    assert cli.main(["plot", str(run.outdir), "--trio", trio]) == 2
    err = capsys.readouterr().err.strip().splitlines()
    assert len(err) == 1 and err[0].startswith(f"error: --trio {trio}: "), err
    assert message.format(sheet=run.sheet) in err[0], err
    assert not [name for name in FIGURE if (run.outdir / name).exists()]


def test_plot_trio_on_a_subgenomes_run_is_a_usage_error(tmp_path, capsys):
    """--trio leads a gene-flow run's figure: on a subgenomes run it is a usage error naming the run."""
    now = datetime.datetime.now(datetime.UTC)
    write_manifest(tmp_path / "manifest.json", ["ltrtrace", "subgenomes", "s.tsv"], {"command": "subgenomes"}, [], {},
                   now, now)
    assert cli.main(["plot", str(tmp_path), "--trio", "R,S,D"]) == 2
    err = capsys.readouterr().err.strip().splitlines()
    assert len(err) == 1 and err[0].startswith("error: --trio R,S,D: ") and "subgenomes" in err[0], err


def test_plot_command_redraws_the_figure(tmp_path, toy_fastas, capsys):
    """`ltrtrace plot RUN` draws the three files of a finished run from the sheet the run recorded, and draws the same
    figure again."""
    run = _analyzed(tmp_path, toy_fastas, tracts=TRACTS)
    files = [run.outdir / f"overview{x}" for x in (".png", ".pdf", "_numbers.tsv")]
    assert cli.main(["plot", str(run.outdir)]) == 0 and all(f.is_file() for f in files)
    first = files[2].read_text()
    files[0].unlink()
    assert cli.main(["plot", str(run.outdir)]) == 0 and files[0].is_file() and files[2].read_text() == first
    assert "Analysis:" in _legend_text(files[1])


def test_plot_reads_the_run_from_its_manifest(tmp_path, toy_fastas, capsys):
    """`ltrtrace plot` finds the samples sheet and the route of a finished run in its manifest.json, and leaves the
    manifest as it was; a run without one is a usage error (exit 2): one line naming the folder and what to run."""
    run = _analyzed(tmp_path, toy_fastas, tracts=TRACTS)
    record = run.outdir / "manifest.json"
    before = record.read_bytes()
    assert cli.main(["plot", str(run.outdir)]) == 0 and record.read_bytes() == before
    assert not (run.outdir / "work" / "geneflow" / "run.json").exists()
    record.unlink()
    capsys.readouterr()
    assert cli.main(["plot", str(run.outdir)]) == 2
    err = capsys.readouterr().err.strip().splitlines()
    assert len(err) == 1 and err[0].startswith("error: ") and str(run.outdir) in err[0], err
    assert "manifest.json" in err[0] and "ltrtrace geneflow" in err[0], err
    with pytest.raises(UsageError, match="manifest.json"):
        overview(run.outdir)


def test_plot_without_stats(tmp_path, toy_fastas, capsys):
    run = simulate_run(tmp_path, seed=1, fasta_dir=toy_fastas)
    assert cli.main(["plot", str(run.outdir)]) == 2
    assert capsys.readouterr().err.strip().splitlines()[-1].startswith("error: ") and not list(
        run.outdir.glob("overview*"))
    assert cli.main(["plot", str(tmp_path / "missing")]) == 2


def test_cli_draws_the_figure_after_geneflow(tmp_path, capsys):
    """`ltrtrace geneflow` draws the overview after the analyses return 0."""
    from analyses.geneflow.test_run import _world

    out = tmp_path / "out"
    assert cli.main(["geneflow", str(_world(tmp_path)), "-o", str(out), "-t", "1"]) == 0
    assert all((out / f"overview{x}").is_file() for x in (".png", ".pdf", "_numbers.tsv"))


def _painting_files(out: Path) -> list[str]:
    """The painting files in a run folder (analyses.geneflow.paint.PAINTING_FILES), by name, sorted."""
    return sorted(p.name for p in out.iterdir() if any(fnmatch.fnmatchcase(p.name, f.format("*"))
                                                       for f in PAINTING_FILES))


@pytest.mark.parametrize("route", ["assembly", "library"])
def test_figure_error_in_run_is_a_warning(tmp_path, monkeypatch, capsys, request, route):
    """After `ltrtrace geneflow` a figure that fails halfway through saving (the PNG written, the PDF not) is one
    WARNING and the exit code stays 0; the analyses' outputs and the run's manifest stay complete, the painting files
    are exactly those work/geneflow/paintings.json records (none on the assembly route's toy world, the trio's on the
    screen of the admixed trio), and no figure or temporary file is left."""
    out = tmp_path / "out"
    if route == "assembly":
        if shutil.which("lastz") is None:
            pytest.skip("lastz (the window aligner) not on PATH")
        from analyses.geneflow.test_run import _world

        argv = ["geneflow", str(_world(tmp_path)), "-o", str(out), "-t", "1"]
        outputs = ("loci.tsv", "stats.tsv", "findings.txt", "inversions.tsv", "mechanisms.tsv", "manifest.json")
    else:
        sc, screened = request.getfixturevalue("admixed")
        shutil.copytree(screened, out)  # its pairs are done: the analyses run
        argv = ["geneflow", "--screen", str(sc.sheet), "-o", str(out), "-t", "2"]
        outputs = ("blocks.tsv", "stats.tsv", "findings.txt", "inversions.tsv", "manifest.json")

    def broken(*args, **kwargs):
        raise KeyError("a planted figure error")

    monkeypatch.setattr(style, "legend_page", broken)
    assert cli.main(argv) == 0
    err = capsys.readouterr().err
    warn = [x for x in err.splitlines() if "figure" in x and "failed" in x]
    assert len(warn) == 1 and "a planted figure error" in warn[0]
    painted = json.loads((out / geneflow_run.PAINTINGS).read_text())
    assert sorted(painted) == _painting_files(out) and bool(painted) == (route == "library"), painted
    assert [name for name in outputs if not (out / name).is_file()] == []
    assert not [name for name in FIGURE if (out / name).exists()] and not list(out.glob(".*.tmp"))


@needs_lastz
def test_failed_geneflow_run_draws_no_figure(tmp_path, monkeypatch, capsys):
    """A geneflow run whose analyses fail (exit 3) draws no figure and records no finished run (no manifest.json), so
    `ltrtrace plot` finds nothing to redraw there (exit 2)."""
    from analyses.geneflow.test_run import _world

    def boom(*args, **kwargs):
        raise ValueError("a planted analysis error")

    monkeypatch.setattr(geneflow_run, "run_summary", boom)
    out = tmp_path / "out"
    assert cli.main(["geneflow", str(_world(tmp_path)), "-o", str(out), "-t", "1"]) == 3
    assert (out / "loci.tsv").is_file() and not (out / "stats.tsv").exists()
    assert not [name for name in (*FIGURE, "manifest.json") if (out / name).exists()]
    capsys.readouterr()
    assert cli.main(["plot", str(out)]) == 2
    assert not [name for name in FIGURE if (out / name).exists()]


@needs_lastz
def test_failed_rerun_leaves_no_record_of_the_earlier_run(tmp_path, monkeypatch, capsys):
    """A geneflow rerun in a finished run's folder whose analyses fail (exit 3) leaves neither the earlier run's
    manifest.json nor its overview files (nor stats.tsv and findings.txt), so nothing there describes results the
    folder no longer holds, and `ltrtrace plot` finds no finished run to redraw (exit 2)."""
    from analyses.geneflow.test_run import _world

    out = tmp_path / "out"
    argv = ["geneflow", str(_world(tmp_path)), "-o", str(out), "-t", "1"]
    assert cli.main(argv) == 0 and all((out / name).is_file() for name in (*FIGURE, "manifest.json"))

    def boom(*args, **kwargs):
        raise ValueError("a planted analysis error")

    monkeypatch.setattr(geneflow_run, "run_summary", boom)
    capsys.readouterr()
    assert cli.main(argv) == 3 and (out / "loci.tsv").is_file()
    assert not [name for name in (*FIGURE, "manifest.json", "stats.tsv", "findings.txt") if (out / name).exists()]
    assert cli.main(["plot", str(out)]) == 2
    assert not [name for name in FIGURE if (out / name).exists()]


# ---------------------------------------------------------------- the example on real sequence


def _bases(rng, n):
    return "".join(np.array(list("ACGT"))[rng.integers(0, 4, n)])


def _sequence_run(tmp_path):
    """A run whose FASTAs hold real sequence: every R+D insertion (R, D members; S empty) is a 150-bp element
    TG...CA between two copies of the toy's TSD ACGTA in R and D, and S holds the empty site with one copy; R+S
    insertions elsewhere. The R+D cluster on chr1 40-80 kb paints a tract."""
    rng = np.random.default_rng(7)
    length, n = 150, 200_000
    rd = [40_001 + 1000 * i for i in range(40)]
    rs = [100_001 + 1000 * i for i in range(60)] + [5_001 + 1000 * i for i in range(30)]
    r = list(_bases(rng, n))
    for p in rd:
        r[p - 6:p - 1] = "ACGTA"
        r[p - 1:p - 1 + length] = "TG" + _bases(rng, length - 4) + "CA"
        r[p - 1 + length:p + 4 + length] = "ACGTA"
    s = list(r)
    for p in rd:   # S: the single TSD copy at p-5..p-1, then R's right flank (no element)
        s[p - 1:p - 1 + 30] = r[p + 4 + length:p + 4 + length + 30]
    fa = tmp_path / "fa"
    fa.mkdir()
    fastas = {"R": write_fasta(fa / "R.fa", {"chr1": "".join(r)}), "D": write_fasta(fa / "D.fa", {"chr1": "".join(r)}),
              "S": write_fasta(fa / "S.fa", {"chr1": "".join(s)}),
              "O": write_fasta(fa / "O.fa", {"chr1": "".join(r)})}
    loci = ([loc(p, "R:P* S:E D:P* O:E") for p in rd] + [loc(p, "R:P* S:P* D:E O:E") for p in rs]
            + [loc(160_001 + 1000 * i, "R:P* S:E D:P* O:E") for i in range(3)])
    return write_run(tmp_path, sorted(loci, key=lambda x: x.pos), seqids=("chr1",), seq_len=n, fastas=fastas)


def test_example_from_sequence(tmp_path):
    run = _sequence_run(tmp_path)
    _analyses(run)
    _, pdf, nums = overview(run.outdir)
    num = _numbers(nums)
    assert "no_example" not in num, num.get("no_example")
    assert num["C_TSD"]["value"] == "ACGTA" and num["C_flank_identity"]["value"] == "20/20 and 20/20"
    start = int(num["C_site"]["value"].split(":")[1].split("-")[0].replace(",", ""))
    assert 40_001 <= start <= 79_001                     # an R+D insertion inside the tract
    assert num["C_lid"]["value"] in _legend_text(pdf)


def test_ends_and_site(tmp_path):
    """The example's sequence reads: an element's flanks, TSD copies and ends, in R's orientation (reversed
    for o = -1), and an empty site's single TSD copy with its flanks (reversed on a '-' window)."""
    rng = np.random.default_rng(3)
    left, right, body = _bases(rng, 40), _bases(rng, 40), "TG" + _bases(rng, 96) + "CA"
    seq = left + "ACGTA" + body + "ACGTA" + right
    fa = Fasta(write_fasta(tmp_path / "x.fa", {"s": seq}))
    start, end = 46, 45 + len(body)                    # 1-based element
    e = _ends(fa, "s", start, end, 5, 1)
    assert (e["ltsd"], e["rtsd"], e["body5"], e["body3"]) == ("ACGTA", "ACGTA", "TG", "CA")
    assert e["lflank"] == left[-20:] and e["rflank"] == right[:20] and e["inner"] == len(body) - 20
    r = _ends(fa, "s", start, end, 5, -1)
    assert r["ltsd"] == revcomp("ACGTA") and r["lflank"] == revcomp(right[:20]) and r["body5"] == "TG"
    empty = write_fasta(tmp_path / "e.fa", {"s": left + "ACGTA" + right})
    site = _site(Fasta(empty), {"b_seqid": "s", "b_start": "41", "b_end": "45", "b_strand": "+"}, 5)
    assert site["tsd"] == "ACGTA" and _ident(site["lflank"], e["lflank"]) == 1.0
    assert _ident(site["rflank"], e["rflank"]) == 1.0
    flip = _site(Fasta(empty), {"b_seqid": "s", "b_start": "41", "b_end": "45", "b_strand": "-"}, 5)
    assert flip["tsd"] == revcomp("ACGTA") and flip["lflank"] == revcomp(right[:20])
    assert _ident("ACGT", "ACGA") == 0.75 and _ident("", "") == 0.0 and _ident("NNNN", "NNNN") == 0.0


# ---------------------------------------------------------------- legend wording: warnings, open points, every trio


def _lead_scope(outdir: Path) -> str:
    return next(r.scope for r in read_stats(outdir / "stats.tsv") if r.test == "roles")


def test_warnings_reach_the_legend_in_plain_words(tmp_path, toy_fastas):
    """A warning on a drawn comparison becomes one plain sentence; every other one is only counted; no stage,
    metric or file name reaches the legend."""
    run = _analyzed(tmp_path, toy_fastas, tracts=TRACTS)
    scope = _lead_scope(run.outdir)
    planted = [
        ("- Mechanisms: WARNING: branches_resolution all: resolved fractions differ by more than 0.1 (R 0.60, S 0.58, "
         "D 0.63, O 0.25); the comparison favors the better-resolved genome (heuristic threshold)"),
        (f"- Painting: WARNING: paint_resolution {scope}: resolved fractions differ by more than 0.1 (S 0.80, D 0.60); "
         "the comparison favors the better-resolved genome (heuristic threshold)"),
        ("- Run: WARNING: the paint part is left out of this file and stats.tsv: its inputs changed (loci.tsv); rerun "
         "ltrtrace geneflow")]
    path = run.outdir / "findings.txt"
    lines = path.read_text().splitlines()
    k = lines.index("Warnings")
    path.write_text("\n".join(lines[:k + 1] + planted) + "\n")
    _, pdf, nums = overview(run.outdir)
    text = " ".join(_legend_text(pdf).split())
    for word in ("branches_resolution", "paint_resolution", "WARNING", "ltrtrace geneflow", "findings.txt", "heuristic",
                 "stats.tsv"):
        assert word not in text
    assert "resolved unequally among the genomes (S 80%, D 60%)" in text
    assert "2 further warnings are listed with the run's findings" in text
    assert _numbers(nums)["warnings_other"]["value"] == "2"


def _open_points(outdir: Path) -> bool:
    rows = read_stats(outdir / "stats.tsv")
    return any(r.test == "D_pav" and r.statistic == "D" and np.isfinite(r.value) and not np.isfinite(r.z)
               for r in rows)


@pytest.mark.parametrize("k_known, expect", [(True, True), (False, False)])
def test_open_points_sentence_only_when_drawn(tmp_path, toy_fastas, k_known, expect):
    """The legend explains open (untestable) points only when panel F draws one."""
    run = _analyzed(tmp_path, toy_fastas, tracts=TRACTS, k_known=k_known)
    assert _open_points(run.outdir) is expect                  # the toy holds the case the test is about
    text = " ".join(_legend_text(overview(run.outdir)[1]).split())
    assert ("Open points: not testable by the block jackknife." in text) is expect


def _four_ingroup_run(tmp_path, toy_fastas):
    """Four annotated ingroup genomes and an outgroup: four trios."""
    genomes = [("R", "tR", None, True), ("S", "tS", None, True), ("D", "tD", None, True), ("X", "tX", None, True),
               ("O", "tO", "outgroup", False)]
    rng = np.random.default_rng(5)
    pats = ["R:P* S:P* D:E X:E O:E", "R:P* S:E D:P* X:E O:E", "R:E S:P* D:P* X:E O:E", "R:P* S:E D:E X:P* O:E",
            "R:E S:E D:P* X:P* O:E", "R:P* S:P* D:P* X:P* O:P", "R:E S:P* D:E X:P* O:E"]
    weights = [0.25, 0.2, 0.1, 0.1, 0.15, 0.1, 0.1]
    loci = [loc(2000 * i, pats[rng.choice(len(pats), p=weights)], seqid=sid)
            for sid in ("chr1", "chr2") for i in range(1, 400)]
    return write_run(tmp_path, loci, genomes=genomes, seqids=("chr1", "chr2"), fasta_dir=toy_fastas)


def test_overview_every_trio(tmp_path, toy_fastas):
    """A run with several trios draws panel H, and the numbers table holds one D_pav row per trio."""
    run = _four_ingroup_run(tmp_path, toy_fastas)
    _analyses(run)
    trios = {r.scope for r in read_stats(run.outdir / "stats.tsv") if r.test == "roles"}
    png, _, nums = overview(run.outdir)
    rows = [r for r in read_tsv(nums) if r["key"].startswith("H_")]
    assert png.is_file() and len(trios) == 4
    assert sorted(r["source"].split("/")[1] for r in rows) == sorted(trios)
    assert all(r["source"].endswith("/all/D") for r in rows)


# ---------------------------------------------------------------- the probe paths


def _probe_run(base: Path, toy_fastas, case: str):
    """The probe runs: normal, swapped (--trio S,R,D), stale_paint (paint rows
    dropped from stats.tsv, files kept), no_outgroup, d_unannotated (fallback discovery) and scaffolds.
    Returns (run, trio option)."""
    trio = None
    if case in ("normal", "swapped", "stale_paint"):
        run = simulate_run(base, seed=1, fasta_dir=toy_fastas, tracts=TRACTS)
        if case == "swapped":
            trio = ("S", "R", "D")
    elif case in ("no_outgroup", "d_unannotated"):
        genomes = [("R", "tR", None, True), ("S", "tS", None, True), ("D", "tD", None, case != "d_unannotated")]
        pats = ["R:P* S:P* D:E", "R:P* S:E D:P", "R:E S:P* D:P", "R:P* S:E D:E", "R:E S:P* D:E"]
        w = [0.35, 0.25, 0.1, 0.15, 0.15]
        if case == "d_unannotated":
            genomes.append(("X", "tX", None, True))
            pats = [x + " X:E" for x in pats] + ["R:P* S:P* D:P X:P*"]
            w.append(0.05)
        rng = np.random.default_rng(5)
        prob = np.array(w) / sum(w)
        loci = [loc(2000 * i, pats[rng.choice(len(pats), p=prob)], seqid=sid)
                for sid in ("chr1", "chr2") for i in range(1, 400)]
        run = write_run(base, loci, genomes=genomes, seqids=("chr1", "chr2"), fasta_dir=toy_fastas)
    else:
        seqids = ("chr1",) + tuple(f"scaf{i:03d}" for i in range(300))
        loci = ([loc(1000 * i, "R:P* S:E D:P O:E") for i in range(1, 60)]
                + [loc(70_000 + 1000 * i, "R:P* S:P D:E O:E") for i in range(60)]
                + [loc(1001, "R:P* S:E D:P O:E", seqid=f"scaf{i:03d}") for i in range(300)])
        run = write_run(base, loci, seqids=seqids, seq_len=5_000, seq_lens={"chr1": 200_000}, fasta_dir=toy_fastas)
    _analyses(run, trio)
    if case == "stale_paint":   # what assemble does with a stale paint part: its rows leave stats.tsv, files stay
        st = run.outdir / "stats.tsv"
        st.write_text("".join(line for line in st.read_text().splitlines(True)
                              if line.split("\t")[1] not in ("paint", "paint_seqid", "paint_resolution")))
    return run, trio


_BAD = re.compile(r"\b(None|nan|NaN|inf|Inf)\b")


@pytest.mark.parametrize("case", ["normal", "swapped", "stale_paint", "no_outgroup", "d_unannotated", "scaffolds"])
def test_no_missing_value_reaches_a_legend(tmp_path, toy_fastas, case):
    """On every probe path the legend pages hold no None, nan or inf."""
    run, trio = _probe_run(tmp_path, toy_fastas, case)
    _, pdf, nums = overview(run.outdir, trio=trio)
    text = " ".join(_legend_text(pdf).split())
    assert not _BAD.findall(text), _BAD.findall(text)
    assert not any(_BAD.search(r["value"]) for r in read_tsv(nums))


def test_verdicts_follow_their_statistics(tmp_path, toy_fastas):
    """A significantly negative D (roles set by the user against the data) is said so; an untestable age
    gradient is said untestable; the sharing counts' order and the roles' source are stated with the counts."""
    run, trio = _probe_run(tmp_path / "s", toy_fastas, "swapped")
    text = " ".join(_legend_text(overview(run.outdir, trio=trio)[1]).split())
    assert "R, not S, carries more of the insertions D shares with one sister (268 against 155): the test points " \
           "against the roles, which were set by the user." in text
    assert "about equally" not in text and "gene flow from" not in text
    assert "R, not S, shares more with D (434 against 248)." in text            # no "makes ... the recipient"
    run, _ = _probe_run(tmp_path / "n", toy_fastas, "normal")
    text = " ".join(_legend_text(overview(run.outdir)[1]).split())
    assert "R carries 1.7 times as many as S: more than lineage sorting alone is expected to give; the age " \
           "gradient could not be tested." in text
    assert "not concentrated in young" not in text and "holds in all 3 checks" in text
    assert "which makes R the recipient" not in text                          # E's counts don't make R the recipient
    f_letter = re.search(r"\b([A-H]) D_pav = \(n_R\+D", text).group(1)     # the letter panel F has here
    assert f"R, not S, shares more with D (434 against 248); the donor-sharing test (panel {f_letter}) makes R " \
           "the recipient." in text
    assert ("roles from the data: the sisters are the pair sharing the most insertions (panel "
            f"{chr(ord(f_letter) - 1)}), the recipient the sister with more donor-shared insertions in the "
            f"donor-sharing test (panel {f_letter})") in text


def test_stale_painting_is_not_drawn(tmp_path, toy_fastas):
    """Paint rows gone from stats.tsv (a stale part) leave the painting files undrawn, said in plain words."""
    run, _ = _probe_run(tmp_path, toy_fastas, "stale_paint")
    _, pdf, nums = overview(run.outdir)
    num = _numbers(nums)
    assert num["no_painting"]["value"] == "no current painting for this trio"
    assert "posterior_mean" not in num and "n_tracts" not in num and "B_high_target" not in num
    assert "There is no current painting for this trio, so no posterior or tract is drawn." in " ".join(
        _legend_text(pdf).split())


def test_kksc_counts_map_by_label(tmp_path, toy_fastas):
    """The donor-read counts of panel E are mapped by genome label: with the recipient after its sister in the
    sheet (--trio S,R,D), R+D is the {S-label, D} count of the note."""
    run, trio = _probe_run(tmp_path, toy_fastas, "swapped")
    note = next(r.note for r in read_stats(run.outdir / "stats.tsv")
                if r.test == "sharing" and r.statistic == "p_kksc")
    counts = {tuple(sorted(m.group(1).split("+"))): int(m.group(2)) for m in re.finditer(r"n_(\S+?) (\d+)", note)}
    num = _numbers(overview(run.outdir, trio=trio)[2])
    assert num["E_kksc_RD"]["value"] == f"{counts[('D', 'S')]:,}" and num["E_kksc_SD"]["value"] == f"{counts[('D', 'R')]:,}"
    assert counts[("D", "S")] != counts[("D", "R")]


@pytest.mark.parametrize("case", ["d_unannotated", "no_outgroup"])
def test_discovery_and_outgroup_wording(tmp_path, toy_fastas, case):
    """Panel F names the discovering genome from the D_pav row (fallback discovery by X), and the outgroup
    clauses appear only with an outgroup."""
    run, _ = _probe_run(tmp_path, toy_fastas, case)
    text = " ".join(_legend_text(overview(run.outdir)[1]).split())
    assert "outgroup an empty site" not in text and "not present in the outgroup" not in text
    if case == "d_unannotated":
        assert "loci X discovered (fallback discovery: D has no annotation)" in text and "X's own calls" in text
        assert "loci D discovered" not in text
    else:
        assert "has no result: its counted set is empty" in text


def test_numbers_equal_their_tables(tmp_path, toy_fastas):
    """Numbers rows hold the values of their sources, compared as values (not as substrings)."""
    run, _ = _probe_run(tmp_path, toy_fastas, "normal")
    num = _numbers(overview(run.outdir)[2])
    st = {(r.test, r.stratum, r.statistic): r for r in read_stats(run.outdir / "stats.tsv")
          if r.scope in ("R,S;D", "R,S,D", "R,D,S", "S,R,D") or r.test == "sharing"}
    d = st[("D_pav", "all", "D")]
    assert float(num["D_pav"]["value"]) == pytest.approx(d.value, abs=0.005)
    assert float(num["D_pav_z"]["value"]) == pytest.approx(d.z, abs=0.05)
    assert int(num["D_pav_n_RD"]["value"].replace(",", "")) == st[("D_pav", "all", "n_RD")].value
    assert int(num["D_pav_n_SD"]["value"].replace(",", "")) == st[("D_pav", "all", "n_SD")].value
    for key, stat in (("E_R+S", "n_R+S"), ("E_R+D", "n_R+D"), ("E_S+D", "n_S+D")):
        assert int(num[key]["value"].replace(",", "")) == st[("sharing", "all", stat)].value
    assert float(num["posterior_mean"]["value"]) == pytest.approx(st[("paint", "all", "posterior_mean")].value,
                                                                  abs=0.0005)
    assert float(num["decoy_posterior_mean"]["value"]) == pytest.approx(
        st[("paint", "all", "decoy_posterior_mean")].value, abs=0.0005)
    for key, name in (("B_high_target", "posterior.R.tsv"), ("B_high_decoy", "decoy_posterior.S.tsv")):
        post = np.array([float(r["post"]) for r in read_tsv(run.outdir / name)])
        assert float(num[key]["value"].rstrip("%")) == pytest.approx(100 * np.mean(post >= 0.9), abs=0.05)


@pytest.mark.skipif(shutil.which("pdftotext") is None, reason="poppler missing")
def test_long_paths_wrap(tmp_path):
    """A token wider than the page (a long analysis path) breaks after "/" or "_" and stays recoverable."""
    from matplotlib.backends.backend_pdf import PdfPages

    path = "/" + "/".join(f"dir_{i:03d}_{'x' * 9}" for i in range(10)) + "/run"
    assert len(path) > 150
    page = style.text_page([[(f"Analysis: host:{path} · Script: ltrtrace.plot.geneflow", "")]], "DejaVu Sans")
    out = tmp_path / "p.pdf"
    with PdfPages(out) as pp:
        pp.savefig(page)
    text = subprocess.run(["pdftotext", str(out), "-"], capture_output=True, text=True, check=True).stdout
    assert path in "".join(text.split())


@pytest.mark.skipif(os.geteuid() == 0, reason="root reads a file whatever its mode")
def test_unreadable_fasta_omits_the_example(tmp_path, toy_fastas):
    """A FASTA that cannot be read leaves panel C out, said in plain words."""
    run = _sequence_run(tmp_path)
    _analyses(run)
    cat = load_catalog(run.outdir, run.samples, run.sheet)
    fasta = next(s.fasta for s in run.samples if s.label == "D")
    fasta.chmod(0)
    try:
        num = _numbers(overview(run.outdir, cat=cat)[2])
    finally:
        fasta.chmod(0o644)
    assert num["no_example"]["value"] == "no base-pair example: D's genome sequence cannot be read"
    assert "C_lid" not in num


# ---------------------------------------------------------------- significance wording, numbers keys, stale stats.tsv


def test_not_significant_is_not_equal(tmp_path, toy_fastas):
    """At a Z just below Z_SIG with unequal counts, the headline and panel F give both counts and say the
    difference is not significant, never that the counts are about equal."""
    from ltrtrace.base.constants import Z_SIG

    run, _ = _probe_run(tmp_path, toy_fastas, "normal")
    scope = _lead_scope(run.outdir)
    st = run.outdir / "stats.tsv"
    rows = [line.rstrip("\n").split("\t") for line in st.read_text().splitlines(True)]
    for r in rows:
        if r[1:5] == ["D_pav", scope, "all", "D"]:
            r[7] = repr(Z_SIG - 0.1)                                          # z: just below the bar
    st.write_text("".join("\t".join(r) + "\n" for r in rows))
    _, pdf, nums = overview(run.outdir)
    num = _numbers(nums)
    a, b = num["D_pav_n_RD"]["value"], num["D_pav_n_SD"]["value"]
    assert a != b and num["D_pav_z"]["value"] == f"{Z_SIG - 0.1:.1f}"
    text = " ".join(_legend_text(pdf).split())
    verdict = f"R carries {a}, S {b}; the difference is not significant: no evidence of gene flow."   # no Z
    assert text.count(verdict) == 2                                           # the headline and panel F
    assert f"Of the insertions D shares with only one of the sisters, {verdict}" in text
    assert "about equally" not in text and "times as many" not in text


def test_a_numbers_key_keeps_one_text():
    """A legend key written again with the same text returns it; with a different text it is a ValueError."""
    leg = object.__new__(ov._Legend)                     # num and text use only the numbers registry
    leg.n = style.Numbers()
    assert leg.num("n", 1317, "int", "stats.tsv:n_RD") == "1,317" == leg.num("n", 1317.0, "int", "stats.tsv:n_RD")
    assert leg.num("n", None, "int", "stats.tsv:n_RD") is None                # missing: no text, nothing checked
    with pytest.raises(ValueError, match="'n' written twice with different texts: '1,317' and '1,318'"):
        leg.num("n", 1318, "int", "stats.tsv:n_RD")
    assert leg.text("site", "chr1:1-5", "loci.tsv:R_start") == "chr1:1-5" == leg.text("site", "chr1:1-5", "x")
    with pytest.raises(ValueError, match="'site'"):
        leg.text("site", "chr1:1-6", "loci.tsv:R_start")
    assert [k for k, _, _ in leg.n.rows] == ["n", "site"]


# ---------------------------------------------------------------- the gradient read in the excess's direction


def test_gradient_text_follows_the_excess(tmp_path, toy_fastas):
    """With --trio R,S,D and S carrying the excess (D_pav < 0), its young concentration gives dD < 0; panel F
    reads the gradient in the excess's direction and names S, and the headline says S's excess is young."""
    from analyses.geneflow.test_trio import gradient_toy_loci

    run = write_run(tmp_path, gradient_toy_loci(), fasta_dir=toy_fastas)
    _analyses(run, ("R", "S", "D"))
    _, pdf, nums = overview(run.outdir, trio=("R", "S", "D"))
    text = " ".join(_legend_text(pdf).split())
    dd, dz = _numbers(nums)["D_pav_gradient"]["value"], _numbers(nums)["D_pav_gradient_z"]["value"]
    assert dd == "-0.60" and float(dz) <= -3
    assert (f"Youngest minus oldest: {dd} (Z {dz}), read in the direction of the excess, S's: concentrated in young "
            "insertions.") in text
    assert ("S's excess is concentrated in young insertions, consistent with gene flow from D-like donors into S "
            "rather than incomplete lineage sorting alone.") in text
    assert "concentrated in old insertions" not in text


# ---------------------------------------------------------------- stats.tsv's own sentinels


def test_figure_reads_no_findings_txt(tmp_path, toy_fastas):
    """The figure is drawn from the tables: a current run draws with findings.txt deleted."""
    run = _analyzed(tmp_path, toy_fastas, tracts=TRACTS)
    (run.outdir / "findings.txt").unlink()
    assert cli.main(["plot", str(run.outdir)]) == 0 and (run.outdir / "overview.png").is_file()


# ---------------------------------------------------------------- panel H scales with the trios


def _many_genomes_run(tmp_path, toy_fastas, n_in: int):
    """n_in annotated ingroup genomes G1.. and an outgroup, 600 loci with random presence (seed 1): C(n_in, 3)
    trios."""
    labels = [f"G{i}" for i in range(1, n_in + 1)]
    genomes = [(x, f"t{x}", None, True) for x in labels] + [("O", "tO", "outgroup", False)]
    rng = np.random.default_rng(1)
    loci = []
    for seqid in ("chr1", "chr2", "chr3", "chr4"):
        for j in range(150):
            pres = rng.random(n_in) < 0.5
            if not pres.any():
                pres[rng.integers(n_in)] = True
            first = int(np.flatnonzero(pres)[0])
            spec = " ".join(f"{x}:P{'*' if g == first else ''}" if pres[g] else f"{x}:E" for g, x in enumerate(labels))
            loci.append(loc(5_000 + 6_000 * j, spec + " O:E", seqid=seqid, k=float(rng.choice([0.0, 0.005, 0.02]))))
    return write_run(tmp_path, loci, genomes=genomes, fasta_dir=toy_fastas)


def _overlaps(fig, axes) -> list[tuple[str, str]]:
    """Pairs of drawn texts of `axes` (tick labels in view, key and row texts) whose boxes overlap, and texts that
    leave the figure (paired with "figure")."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    rend = canvas.get_renderer()
    texts = []
    for ax in axes:
        lo, hi = sorted(ax.get_xlim())
        texts += [t for t in ax.texts if t.get_text().strip()]
        texts += [t for t in ax.get_yticklabels() if t.get_text().strip()]
        texts += [t for t, x in zip(ax.get_xticklabels(), ax.get_xticks(), strict=True)
                  if t.get_text().strip() and lo <= x <= hi]
    boxes = [(t.get_text(), t.get_window_extent(rend)) for t in texts]
    W, H = fig.bbox.width, fig.bbox.height
    out = [(s, "figure") for s, b in boxes if b.x0 < 0 or b.y0 < 0 or b.x1 > W or b.y1 > H]
    for i, (a, ba) in enumerate(boxes):
        for b, bb in boxes[i + 1:]:
            if min(ba.x1, bb.x1) - max(ba.x0, bb.x0) > 0.5 and min(ba.y1, bb.y1) - max(ba.y0, bb.y0) > 0.5:
                out.append((a, b))
    return out


def test_panel_h_twenty_trios(tmp_path, toy_fastas, monkeypatch):
    """Six ingroup genomes (20 trios): every trio is a labeled row of panel
    H, and no text of the panel collides with another or leaves the figure."""
    run = _many_genomes_run(tmp_path, toy_fastas, 6)
    _analyses(run)
    figs, overlapped = [], []
    real = ov._draw

    def draw(m, mono):
        """Draw, then look for collisions while the figure's font settings are still in force (a family name
        resolves when the text is drawn; a redraw after overview() returns would use matplotlib's default
        fonts)."""
        fig = real(m, mono)
        figs.append(fig)
        overlapped.append(_overlaps(fig, [ax for ax in fig.axes if ax.get_gid() == "H"]))
        return fig

    monkeypatch.setattr(ov, "_draw", draw)
    nums = _numbers(overview(run.outdir)[2])
    axes = [ax for ax in figs[0].axes if ax.get_gid() == "H"]
    labels = [t.get_text() for ax in axes for t in ax.get_yticklabels() if t.get_text()]
    with_d = {r.scope for r in read_stats(run.outdir / "stats.tsv")          # trios with a D_pav value
              if (r.test, r.stratum, r.statistic) == ("D_pav", "all", "D") and np.isfinite(r.value)}
    assert len(labels) == 20 == len(set(labels)) == int(nums["trios_n"]["value"])
    assert sorted(r["source"].split("/")[1] for k, r in nums.items() if k.startswith("H_")) == sorted(with_d)
    assert overlapped == [[]]


def test_panel_h_84_trios(tmp_path, toy_fastas, capsys):
    """Nine ingroup genomes (84 trios): `ltrtrace plot` draws, with panel H growing to fit the trios."""
    run = _many_genomes_run(tmp_path, toy_fastas, 9)
    _analyses(run)
    assert cli.main(["plot", str(run.outdir)]) == 0
    nums = _numbers(run.outdir / "overview_numbers.tsv")
    assert int(nums["trios_n"]["value"]) == 84 and 0 < sum(k.startswith("H_") for k in nums) <= 84


# ---------------------------------------------------------------- the headline, the short legend


def _short_legend(pdf: Path) -> str:
    return " ".join((_pdf_pages(pdf, 2, 2) if shutil.which("pdftotext") else _PAGES[-2]).split())


_STATISTICS = re.compile(r"q ≤|\|Z\||\bZ\b|K = 0|K ≥|\bD [+−-]?\d|[+−]\d\.\d")   # thresholds, Z, D values


def _unannotated_run(base: Path, toy_fastas, case: str):
    """S genotype-only (no decoy), or R genotype-only (nothing painted)."""
    annotated, specs = {
        "S_unannotated": ({"R", "D"}, [("R:P S:E D:P* O:E", 0.0, 3), ("R:E S:P D:P* O:E", 0.02, 1),
                                       ("R:P* S:P D:E O:E", 0.02, 4), ("R:P* S:E D:P O:E", 0.0, 2),
                                       ("R:P* S:P D:E O:E", 0.01, 2)]),
        "R_unannotated": ({"S", "D"}, [("R:P S:E D:P* O:E", 0.0, 3), ("R:E S:P D:P* O:E", 0.02, 1),
                                       ("R:P S:P* D:E O:E", 0.02, 4), ("R:E S:P* D:P O:E", 0.0, 1),
                                       ("R:P S:P* D:E O:E", 0.01, 2)])}[case]
    genomes = [("R", "taxR", None, "R" in annotated), ("S", "taxS", None, "S" in annotated),
               ("D", "taxD", None, "D" in annotated), ("O", "taxO", "outgroup", False)]
    loci = []
    for seqid in ("chr1", "chr2", "chr3", "chr4"):
        pos = 5_000
        for _ in range(6):
            for spec, k, n in specs:
                for _ in range(n):
                    loci.append(loc(pos, spec, seqid=seqid, k=k))
                    pos += 6_000
    run = write_run(base, loci, genomes=genomes, fasta_dir=toy_fastas)
    _analyses(run)
    return run


@pytest.mark.parametrize("case", ["normal", "swapped", "no_outgroup"])
def test_short_legend_gives_no_statistic(tmp_path, toy_fastas, case):
    """The short legend states results in plain biology: no threshold, q-value, Z or D value (the numbers stay
    on the long legend)."""
    run, trio = _probe_run(tmp_path, toy_fastas, case)
    short = _short_legend(overview(run.outdir, trio=trio)[1])
    assert not _STATISTICS.findall(short), _STATISTICS.findall(short)


def test_headline_says_consistent_with(tmp_path, toy_fastas):
    """A significant excess concentrated in young insertions reads "consistent with gene flow ... rather than
    incomplete lineage sorting alone" (D_pav's null is not calibrated)."""
    from analyses.geneflow.test_trio import gradient_toy_loci

    run = write_run(tmp_path, gradient_toy_loci(), fasta_dir=toy_fastas)
    _analyses(run, ("S", "R", "D"))
    short = _short_legend(overview(run.outdir, trio=("S", "R", "D"))[1])
    assert ("S carries 1.6 times as many as R, most of all among young insertions: consistent with gene flow from "
            "D-like donors into S rather than incomplete lineage sorting alone.") in short
    assert "not incomplete lineage sorting" not in short


def test_headline_names_the_lead_among_several_trios(tmp_path, toy_fastas):
    """With several trios the headline says the lead trio is the strongest of n, or that the user named it."""
    run = _many_genomes_run(tmp_path, toy_fastas, 6)
    _analyses(run)
    short = _short_legend(overview(run.outdir)[1])
    assert re.search(r"This is the strongest of the run's 20 trios; panel [A-H] shows them all\.", short)
    short = _short_legend(overview(run.outdir, trio=("G6", "G1", "G2"))[1])
    assert re.search(r"This trio was named by the user; panel [A-H] shows all 20 of the run's trios\.", short)


def test_headline_states_an_uncertain_sister_pair(tmp_path, toy_fastas):
    """A sister pair the sharing counts do not settle (p_top >= P_WARN) is stated in the headline, with the two
    counts, not only counted among the warnings."""
    from ltrtrace.base.constants import P_WARN
    run = _unannotated_run(tmp_path, toy_fastas, "S_unannotated")
    p_top = next(r.value for r in read_stats(run.outdir / "stats.tsv") if r.statistic == "p_top")
    assert p_top >= P_WARN                                                     # the toy holds the case
    short = _short_legend(overview(run.outdir)[1])
    assert ("Which two of R, S and D are sisters is uncertain: R and S share 144 insertions, R and D 120, a difference "
            "within chance, so the roles may be swapped.") in short


# ---------------------------------------------------------------- generality paths


def test_genotype_only_sister_text(tmp_path, toy_fastas):
    """S genotype-only, so there is no decoy painting: panel A never says the model was applied to S or that
    tracts are tested against a decoy, and the collapsed fit is said in plain words."""
    run = _unannotated_run(tmp_path, toy_fastas, "S_unannotated")
    ds = next(r for r in read_stats(run.outdir / "stats.tsv") if r.test == "paint" and r.statistic == "donor_state")
    assert not np.isfinite(ds.value)                                          # the toy's fit collapsed
    _, pdf, _ = overview(run.outdir)
    text = " ".join(_legend_text(pdf).split())
    assert "applied unchanged to S" not in text and "against the decoy" not in text
    assert "with no decoy painting (S discovers no insertions to paint), no tract has a q-value." in text
    assert "No D-like ancestry is detectable along R's genome." in _short_legend(pdf)


def test_genotype_only_recipient_text(tmp_path, toy_fastas):
    """R genotype-only, so nothing is painted: panel A never says S is painted, says why R holds no insertions
    of its own, and panel E names only the panels that exist."""
    run = _unannotated_run(tmp_path, toy_fastas, "R_unannotated")
    _, pdf, nums = overview(run.outdir)
    text = " ".join(_legend_text(pdf).split())
    assert "painted the same way" not in text and "same painting of its sister" not in text
    assert "R has no annotation, so it holds no insertions of its own to class or paint" in text
    assert _numbers(nums)["no_painting"]["value"] == "no painting: R has no annotation, so it discovers no insertions to paint"
    assert not re.search(r"panels ([A-H]) and \1\b", text) and re.search(r"colored as in panel [A-H] \(blue", text)


def test_sequences_or_chromosomes(tmp_path, toy_fastas):
    """The legends call the drawn sequences chromosomes only when they are chromosome-scale by the drawing rule
    (never by name): 300 scaffolds of which the 24 largest are drawn are sequences; four 1-Mb ones are chromosomes."""
    run, _ = _probe_run(tmp_path / "scaf", toy_fastas, "scaffolds")
    short = _short_legend(overview(run.outdir)[1])
    assert "R's sequences" in short and "chromosome" not in short
    run, _ = _probe_run(tmp_path / "chr", toy_fastas, "normal")
    short = _short_legend(overview(run.outdir)[1])
    assert "R's chromosomes" in short and "chromosome pairs" in short


# ---------------------------------------------------------------- a figure on disk is whole


def test_failed_redraw_keeps_the_figure_that_was_there(tmp_path, toy_fastas, monkeypatch, capsys):
    """A redraw that fails, while drawing or halfway through saving (the PNG written, the PDF not), exits 1 and
    leaves the figure files as they were, whole: the earlier figure, or none; no temporary file stays. A run folder
    without stats.tsv gets no figure."""
    run = _analyzed(tmp_path, toy_fastas, tracts=TRACTS)
    figure = [run.outdir / name for name in FIGURE]
    assert cli.main(["plot", str(run.outdir)]) == 0
    before = {p.name: p.read_bytes() for p in figure}

    def fail(*args, **kwargs):
        raise RuntimeError("a planted drawing error")

    for earlier in (True, False):
        if not earlier:
            for p in figure:
                p.unlink()
        for module, name in ((ov, "_legends"), (style, "legend_page")):  # while drawing; inside the PDF's writing
            with monkeypatch.context() as m:
                m.setattr(module, name, fail)
                capsys.readouterr()
                assert cli.main(["plot", str(run.outdir)]) == 1
            assert capsys.readouterr().err.strip().splitlines()[-1] == "error: RuntimeError: a planted drawing error"
            assert {p.name: p.read_bytes() for p in figure if p.exists()} == (before if earlier else {}), name
            assert not [p.name for p in run.outdir.iterdir() if p.name.endswith(".tmp")], name
    (run.outdir / "stats.tsv").unlink()
    assert cli.main(["plot", str(run.outdir)]) == 2 and not any(p.exists() for p in figure)


# ---------------------------------------------------------------- the script's source path and the fonts


def test_last_line_gives_the_script_path(tmp_path, toy_fastas):
    """The long legend's last line names the script by its source path, not by its module name."""
    run = _analyzed(tmp_path, toy_fastas, tracts=TRACTS)
    script = os.path.abspath(ov.__file__)
    assert f"· Script: {script} ·" in ov._analysis_line(run.outdir)
    text = "".join(_legend_text(overview(run.outdir)[1]).split())
    assert f"Script:{script}·" in text and "Script:ltrtrace.plot.geneflow" not in text


def test_png_is_byte_deterministic(tmp_path, toy_fastas):
    """Two renders in fresh interpreters whose hash seeds order the system's font list
    differently (seeds 0 and 1 can order duplicate font files, such as two LiberationSans-Regular.ttf
    versions, differently) give byte-identical PNGs."""

    run = _analyzed(tmp_path, toy_fastas, tracts=TRACTS)
    code = ("import sys; from pathlib import Path; from ltrtrace.plot.geneflow import overview; "
            "overview(Path(sys.argv[2]))")
    pngs = []
    for seed in ("0", "1"):
        subprocess.run([sys.executable, "-c", code, str(run.sheet), str(run.outdir)], check=True, capture_output=True,
                       env={**os.environ, "PYTHONHASHSEED": seed})
        pngs.append((run.outdir / "overview.png").read_bytes())
    assert pngs[0] == pngs[1]


def test_font_fallback_is_one_warning():
    """Without Liberation or Arial fonts (none in the system's folders, none in matplotlib's font list, which
    matplotlib 3.10 fills from them when it first runs) the figure falls back to DejaVu with one WARNING
    naming both fallbacks."""

    code = ("import logging, matplotlib.font_manager as fm; fm.findSystemFonts = lambda *a, **k: []; "
            "fm.fontManager.ttflist[:] = [f for f in fm.fontManager.ttflist if 'Liberation' not in f.name "
            "and 'Arial' not in f.name]; "
            "logging.basicConfig(level=logging.WARNING, format='%(levelname)s %(message)s'); "
            "from ltrtrace.plot import style; print((style.setup_fonts(), style.rc_params()['mathtext.tt']))")
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    assert r.stdout.strip() == "('DejaVu Sans', 'DejaVu Sans Mono')"
    warn = [x for x in r.stderr.splitlines() if x.startswith("WARNING")]
    assert len(warn) == 1 and "DejaVu Sans for Liberation Sans or Arial and DejaVu Sans Mono for Liberation Mono" in warn[0]

