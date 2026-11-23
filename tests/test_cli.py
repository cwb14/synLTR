"""The command line (synltr.cli): the thread environment set before numpy or numba load, exit codes and error lines,
the commands and their options, help that fits one screen, routing and input checks, the manifest.json a finished run
writes, and `synltr genotype` end to end: three toy genomes (A and B annotated, C genotype-only), reruns that skip
every stage, and the fail-fast errors."""
from __future__ import annotations

import argparse
import contextlib
import gzip
import hashlib
import importlib.util
import json
import logging
import os
import random
import re
import shutil
import signal
import subprocess
import sys
import time
import types
from dataclasses import dataclass, replace
from pathlib import Path

import builders
import pysam
import pytest
from markers import needs_lastz

from synltr import __version__, cli
from synltr.analyses.geneflow.run import SCREEN_ROUTED
from synltr.base import checkpoint
from synltr.base.constants import TOOL_VERSION_ARGS
from synltr.base.log import AnalysisError, SynltrError
from synltr.calls import manifest
from synltr.io.elements import Element, make_eid, read_elements
from synltr.io.samples import Sample

ENV = {**os.environ, "NUMPY_MADVISE_HUGEPAGE": "0", "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"}
REPO = Path(__file__).resolve().parents[1]
COMMANDS = ("subgenomes", "genotype", "geneflow", "source", "plot", "sim")
FIGURE = ("overview.png", "overview.pdf", "overview_numbers.tsv")


# ---------------------------------------------------------------- the thread environment


def test_thread_env_sets_numba_from_t():
    assert cli.thread_env(["subgenomes", "-t", "3"])["NUMBA_NUM_THREADS"] == "3"
    assert cli.thread_env(["subgenomes", "--threads=7"])["NUMBA_NUM_THREADS"] == "7"
    assert cli.thread_env(["subgenomes", "-t5"])["NUMBA_NUM_THREADS"] == "5"
    assert cli.thread_env(["subgenomes"])["NUMBA_NUM_THREADS"] == "4"


def test_thread_env_single_threaded_libraries():
    """BLAS and OpenMP run one thread and numba's pool is -t, whatever the environment says: a user's
    OMP_NUM_THREADS would break the -t budget."""
    user = {"OMP_NUM_THREADS": "8", "OPENBLAS_NUM_THREADS": "8", "MKL_NUM_THREADS": "8", "NUMBA_NUM_THREADS": "64"}
    env = cli.thread_env([], user)
    assert {k: env[k] for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")} == {
        "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    assert env["NUMBA_NUM_THREADS"] == "4"
    assert env["OMP_WAIT_POLICY"] == "PASSIVE"
    assert env["KMP_BLOCKTIME"] == "0"


def test_thread_env_honors_a_user_hugepage_setting():
    """NUMPY_MADVISE_HUGEPAGE is 0 unless the user set it; a value the user set stands."""
    assert cli.thread_env([], {})["NUMPY_MADVISE_HUGEPAGE"] == "0"
    assert cli.thread_env([], {"NUMPY_MADVISE_HUGEPAGE": "1"})["NUMPY_MADVISE_HUGEPAGE"] == "1"
    assert cli.thread_env(["-t", "2"], {"NUMPY_MADVISE_HUGEPAGE": "0"})["NUMPY_MADVISE_HUGEPAGE"] == "0"


def test_import_cli_loads_no_numpy():
    out = subprocess.run([sys.executable, "-c", "import synltr.cli, sys; print('numpy' in sys.modules)"],
                         capture_output=True, text=True, check=True, env=ENV).stdout
    assert out.strip() == "False"


def test_building_parser_and_sim_help_load_no_heavy_modules():
    """Every subparser is built, and `sim pop --help` printed with the defaults of the simulator, without numpy,
    msprime or the simulator's own modules."""
    code = ("import sys\nfrom synltr import cli\ncli._parser()\n"
            "try:\n    cli.main(['sim', 'pop', '--help'])\nexcept SystemExit:\n    pass\n"
            "print([m for m in ('numpy', 'msprime', 'numba', 'synltr.sim.pop.model', "
            "'synltr.sim.nucleus.model') if m in sys.modules])")
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True, env=ENV)
    assert r.stdout.strip().splitlines()[-1] == "[]"
    assert "ne_a 30000" in " ".join(r.stdout.split())


def test_numba_threads_capped():
    env = {**ENV, **cli.thread_env(["x", "-t", "1"])}
    code = "import numba; print(numba.get_num_threads())"
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True, env=env).stdout
    assert out.strip() == "1"


# ---------------------------------------------------------------- exit codes and error lines


def _command(monkeypatch, exc, result=0):
    """Register a command named `boom` that raises `exc`, or returns `result` when `exc` is None."""
    mod = types.ModuleType("fake_cmd")

    def run(args):
        if exc is not None:
            raise exc
        return result(args) if callable(result) else result

    mod.run = run
    monkeypatch.setitem(sys.modules, "fake_cmd", mod)
    monkeypatch.setitem(cli.COMMANDS, "boom", "fake_cmd:run")


def test_exit_codes(monkeypatch, capsys):
    _command(monkeypatch, None)
    assert cli.main(["boom"]) == 0
    assert cli.main(["boom", "--no-such-option"]) == 2
    capsys.readouterr()
    _command(monkeypatch, SynltrError("bad input"))
    assert cli.main(["boom"]) == 1
    assert capsys.readouterr().err.splitlines() == ["error: bad input"]
    _command(monkeypatch, AnalysisError("analyses failed"))
    assert cli.main(["boom"]) == 3
    capsys.readouterr()
    _command(monkeypatch, KeyboardInterrupt())
    assert cli.main(["boom"]) == 130
    assert capsys.readouterr().err.splitlines() == ["error: interrupted"]


def test_warnings_are_marked(monkeypatch, capsys):
    """A command's log lines on stderr: a milestone as it is, a warning after "warning: " (as an error follows
    "error: "), per-step progress only with -v."""
    def run(args):
        log = logging.getLogger("synltr.boom")
        log.info("a milestone")
        log.warning("a warning")
        log.debug("a step")
        return 0

    _command(monkeypatch, None, run)
    assert cli.main(["boom"]) == 0
    assert capsys.readouterr().err.splitlines() == ["a milestone", "warning: a warning"]
    assert cli.main(["boom", "-v"]) == 0
    assert capsys.readouterr().err.splitlines() == ["a milestone", "warning: a warning", "a step"]


def test_verbose_error_adds_traceback(monkeypatch, capsys):
    _command(monkeypatch, SynltrError("bad input"))
    assert cli.main(["boom", "-v"]) == 1
    err = capsys.readouterr().err
    assert err.startswith("error: bad input") or "error: bad input" in err.splitlines()
    assert "Traceback" in err


def test_no_command_is_usage_error(capsys):
    assert cli.main([]) == 2
    assert capsys.readouterr().err.startswith("error:")


def test_version():
    out = subprocess.run([sys.executable, "-m", "synltr", "--version"], capture_output=True, text=True,
                         check=True, env=ENV).stdout
    assert out.strip() == "synltr 0.2.0"


@pytest.mark.parametrize("flag", ["-t", "--threads", "-o", "--outdir", "-v", "--verbose"])
def test_common_options_on_every_command(monkeypatch, flag):
    _command(monkeypatch, None)
    value = [] if flag in ("-v", "--verbose") else ["2" if "t" in flag[:3] else "out"]
    assert cli.main(["boom", flag, *value]) == 0


def test_unexpected_exception_is_one_error_line(monkeypatch, capsys):
    _command(monkeypatch, KeyError("missing"))
    assert cli.main(["boom"]) == 1
    assert capsys.readouterr().err.splitlines() == ["error: KeyError: 'missing'"]
    assert cli.main(["boom", "-v"]) == 1
    err = capsys.readouterr().err
    assert err.startswith("error: KeyError: 'missing'") and "Traceback" in err
    _command(monkeypatch, OSError("disk full"))
    assert cli.main(["boom"]) == 1
    assert capsys.readouterr().err.splitlines() == ["error: OSError: disk full"]


def test_abbreviated_option_is_rejected(monkeypatch, capsys):
    _command(monkeypatch, None)
    assert cli.main(["boom", "--thr", "3"]) == 2
    assert cli.main(["boom", "--verb"]) == 2


@pytest.mark.parametrize("argv", [["-vt3"], ["-vt", "3"], ["-t3", "-v"], ["-vvt3"]])
def test_clustered_short_options_agree_with_argparse(monkeypatch, argv):
    _command(monkeypatch, None)
    assert cli.thread_env(["boom", *argv])["NUMBA_NUM_THREADS"] == "3"
    assert cli.main(["boom", *argv]) == 0
    seen = {}
    mod = sys.modules["fake_cmd"]
    mod.run = lambda args: seen.setdefault("t", args.threads) and 0
    cli.main(["boom", *argv])
    assert seen["t"] == 3


@pytest.mark.parametrize("value", ["0", "-1"])
def test_threads_below_one_is_usage_error(monkeypatch, capsys, value):
    _command(monkeypatch, None)
    assert cli.main(["boom", "-t", value]) == 2
    assert cli.main(["boom", f"--threads={value}"]) == 2


# ---------------------------------------------------------------- the commands, their options and their help


def _subparsers(parser: argparse.ArgumentParser) -> dict[str, argparse.ArgumentParser]:
    return next(a for a in parser._actions if isinstance(a, argparse._SubParsersAction)).choices


def test_commands_and_options_as_designed():
    """The six commands, and each one's options: -t/--threads, -v/--verbose, -o/--outdir (not plot: its run folder
    is both its input and its output) and --version, then its own."""
    own = {"subgenomes": {"-n", "--assign", "--mu", "--seed"},
           "genotype": {"--pair", "--paf", "--ploidy", "--min-mapq"},
           "geneflow": {"--screen", "--trio", "--identity", "--paf"}, "source": set(), "plot": {"--trio"}}
    commands = _subparsers(cli._parser())
    assert set(commands) == set(COMMANDS)
    for name, extra in own.items():
        options = {s for a in commands[name]._actions for s in a.option_strings}
        common = {"-h", "--help", "-t", "--threads", "-v", "--verbose", "--version"}
        assert options == common | ({"-o", "--outdir"} if name != "plot" else set()) | extra, name
    sim = _subparsers(commands["sim"])
    assert set(sim) == {"pop", "nucleus", "score"}
    for name, sub in sim.items():
        options = {s for a in sub._actions for s in a.option_strings}
        assert {"-t", "--threads", "-v", "--verbose", "-o", "--version"} <= options, name


@pytest.mark.parametrize("width", ["80", "100"])
@pytest.mark.parametrize("argv", [[], *([c] for c in COMMANDS)], ids=["synltr", *COMMANDS])
def test_help_fits_one_screen(monkeypatch, capsys, argv, width):
    """`synltr --help` and each command's --help: at most 40 lines of at most 100 characters, on a terminal 80 or
    100 columns wide."""
    monkeypatch.setenv("COLUMNS", width)
    assert cli.main([*argv, "--help"]) == 0
    lines = capsys.readouterr().out.splitlines()
    assert 5 < len(lines) <= 40 and max(map(len, lines)) <= 100, "\n".join(lines)


DROPPED = (["--anchor", "mm2"], ["--vcf-ref", "A"], ["--block-size", "1000000"], ["--no-figure"], ["--label", "x"],
           ["--dev_scale", "1"])


@pytest.mark.parametrize("command", ["subgenomes", "genotype", "geneflow", "source", "plot"])
def test_dropped_options_absent(tmp_path, capsys, src_root, command):
    """The options the four tools had and synltr drops are usage errors on every command, and src/ and the README
    name none of them."""
    for option in DROPPED:
        assert cli.main([command, str(tmp_path / "sheet.tsv"), *option]) == 2, option
        assert f"unrecognized arguments: {' '.join(option)}" in capsys.readouterr().err
    assert not list(tmp_path.iterdir())
    text = "\n".join(p.read_text() for p in [*src_root.rglob("*.py"), REPO / "README.md"])
    assert [o[0] for o in DROPPED if o[0] in text] == []


def test_user_facing_text_names_synltr_commands(src_root):
    """No message, log or findings line of the package names a command of the tools it replaces (ltrpav-stats,
    ltrpav-merge and the other ltrpav-* commands, ltrmix.py, ltrharvest_phasing10.py, ltrkin's subcommands): the line
    asking for the sequence-identity layer names `synltr geneflow` and --identity."""
    from synltr.analyses.geneflow import run as geneflow

    old = re.compile(r"ltrpav-\w|\bltrpav (run|stats|merge|anchor)\b|ltrmix\.py|ltrharvest_phasing(10|\d+\.py)|"
                     r"\bltrkin (genotype|source)\b")
    hits = [f"{p.relative_to(src_root)}:{n}: {line.strip()}" for p in sorted(src_root.rglob("*.py"))
            for n, line in enumerate(p.read_text().splitlines(), 1) if old.search(line)]
    assert hits == []
    assert "rerun synltr geneflow with --identity" in geneflow._IDENTITY_LINE


def test_one_console_script_and_every_command_has_version(capsys):
    """One console script, `synltr` (pyproject's [project.scripts]); --version prints `synltr <version>` on its
    own and after every command and simulation."""
    scripts = re.search(r"\[project\.scripts\]\n(.*?)\n\n", (REPO / "pyproject.toml").read_text(), re.DOTALL).group(1)
    assert [tuple(re.findall(r'[^\s="]+', ln)) for ln in scripts.splitlines()] == [("synltr", "synltr.cli:entry")]
    for argv in ([], *([c] for c in COMMANDS), ["sim", "pop"], ["sim", "nucleus"], ["sim", "score"]):
        assert cli.main([*argv, "--version"]) == 0, argv
        assert capsys.readouterr().out == f"synltr {__version__}\n", argv


def test_program_label_is_built_in_one_place(src_root):
    """Every output that names the program and its version (VCF ##source, findings headers, figure footers, GFF3
    source column, --version) takes it from base.constants (PROGRAM, GFF3_SOURCE): no other module builds
    "synltr <version>" itself."""
    built = re.compile(r"synltr[ -]\{")
    hits = [f"{p.relative_to(src_root)}:{n}" for p in sorted(src_root.rglob("*.py")) if p.name != "constants.py"
            for n, line in enumerate(p.read_text().splitlines(), 1) if built.search(line)]
    assert hits == []


def test_pyproject_version_is_the_package_version():
    """The version is written twice, in pyproject.toml and as synltr.__version__: the two are the same."""
    import tomllib

    assert tomllib.loads((REPO / "pyproject.toml").read_text())["project"]["version"] == __version__


# ---------------------------------------------------------------- routing and input checks


def test_sheet_without_fasta_runs_screen(tmp_path):
    """`synltr geneflow` (the console entry, in a fresh interpreter) on a sheet without genome FASTAs (label,
    library, role) runs the library screen and says so: one log line, one findings line, and route library in its
    manifest.json, with the checksums of the library files the screen hashed."""
    sc = builders.trio(tmp_path / "in", seed=3, admix=True)
    out = tmp_path / "run"
    r = subprocess.run([sys.executable, "-m", "synltr", "geneflow", str(sc.sheet), "-o", str(out), "-t", "2"],
                       capture_output=True, text=True, env=ENV, check=False)
    assert r.returncode == 0, r.stderr
    assert sum("holds no genome FASTA, so geneflow runs the library screen" in ln
               for ln in r.stderr.splitlines()) == 1, r.stderr
    assert (out / "findings.txt").read_text().count(SCREEN_ROUTED) == 1
    doc = json.loads((out / "manifest.json").read_text())
    assert doc["route"] == "library" and doc["params"]["command"] == "geneflow" and doc["params"]["screen"] is False
    for s in doc["samples"]:
        files = [Path(p) for p in s["library"].split(",")]
        assert s["checksums"] == {f"library:{f.name}": _sha256(f) for f in files}, s["label"]
    assert all((out / name).is_file() for name in FIGURE)


def _reads_sheets(d: Path) -> tuple[Path, Path]:
    """A sheet of one reference and a read set, and the same with an assembly row (B: a FASTA, no reads)."""
    builders.write_fasta(d / "ref.fa", {"chr1": "ACGT" * 100})
    builders.write_fasta(d / "asm.fa", {"chr1": "ACGT" * 100})
    (d / "ref.gff3").write_text("##gff-version 3\n")
    (d / "s1.bam").write_bytes(b"")
    head = "label\tfasta\tannotation\treads\trole\nref\tref.fa\tref.gff3\t.\treference\ns1\t.\t.\ts1.bam\t.\n"
    (d / "reads.tsv").write_text(head)
    (d / "mixed.tsv").write_text(head + "B\tasm.fa\t.\t.\t.\n")
    return d / "reads.tsv", d / "mixed.tsv"


def test_mixed_reads_and_assemblies_rejected(tmp_path, capsys):
    """A sheet holding read sets and an assembly (a row with a FASTA but no reads, besides the reference): `genotype`
    and `geneflow` exit 2 with one error line naming the row; `geneflow` on a sheet of read sets is a usage error too
    (it genotypes assemblies or compares libraries). Nothing is written."""
    reads, mixed = _reads_sheets(tmp_path)
    for command in ("genotype", "geneflow"):
        assert cli.main([command, str(mixed), "-o", str(tmp_path / "out")]) == 2
        err = capsys.readouterr().err.strip().splitlines()
        assert len(err) == 1 and err[0].startswith("error: ") and "row(s) B have no reads" in err[0], err
    assert cli.main(["geneflow", str(reads), "-o", str(tmp_path / "out")]) == 2
    err = capsys.readouterr().err.strip().splitlines()
    assert len(err) == 1 and "read sets" in err[0] and "synltr genotype" in err[0], err
    assert not (tmp_path / "out").exists()


SUBGENOMES_TOY = {"root_age": 6e6, "chromosomes": {"n": 2, "length": 4e6}, "families": {"n": 60},
                  "activity": {"rate": 0.000225, "p_burst": 0.0}, "exchanges": {"nonreciprocal": 0, "reciprocal": 0},
                  "tree": {"split": 2e6, "children": ["A", "B"]}, "mergers": [{"time": 1e5, "lineages": ["A", "B"]}]}


@pytest.fixture(scope="module")
def subgenomes_run(tmp_path_factory):
    """`synltr subgenomes` on a small simulated allotetraploid with an empty PATH (no lastz, no external program
    at all): (sheet, run folder, exit code)."""
    d = tmp_path_factory.mktemp("subgenomes")
    lib, _ = builders.simulate(SUBGENOMES_TOY, d / "in", d / "truth", seed=1)
    sheet = d / "sheet.tsv"
    sheet.write_text("label\tlibrary\tannotation\tfasta\n"
                     f"toy\t{lib / 'ltr.fa'}\t{lib / 'annotation.gff3'}\t{lib / 'genome.fai'}\n")
    (d / "bin").mkdir()
    with pytest.MonkeyPatch.context() as mp:
        mp.setenv("PATH", str(d / "bin"))
        code = cli.main(["subgenomes", str(sheet), "-o", str(d / "run"), "-t", "2"])
    logging.getLogger("synltr").handlers.clear()
    return sheet, d / "run", code


def test_external_programs_checked_when_needed(subgenomes_run, tmp_path, capsys, monkeypatch):
    """`subgenomes` runs to the end with no lastz (no program at all) on PATH; `genotype` on assemblies, which aligns
    with lastz, exits 1 with one error line naming it, before anything is written."""
    _, run, code = subgenomes_run
    assert code == 0 and (run / "findings.txt").is_file() and (run / "report.pdf").is_file()
    sheet = _mini_sheet(tmp_path)
    (tmp_path / "bin").mkdir()
    monkeypatch.setenv("PATH", str(tmp_path / "bin"))
    assert cli.main(["genotype", str(sheet), "-o", str(tmp_path / "out")]) == 1
    err = capsys.readouterr().err.strip().splitlines()
    assert len(err) == 1 and err[0].startswith("error: ") and "lastz" in err[0], err
    assert not (tmp_path / "out").exists()


# ---------------------------------------------------------------- manifest.json


def test_a_finished_run_records_itself(monkeypatch, tmp_path, capsys):
    """A command that returns its finished run's record gets manifest.json in that run's folder: the command line as
    given, the parsed options (paths absolute), the run's own values over them, its samples and route, the code
    fingerprint of the command's module, the tools and the times. A command that returns an exit code, or fails,
    gets none."""
    sheet = _mini_sheet(tmp_path)
    out = tmp_path / "out"
    samples = (Sample("A", tmp_path / "A.fa", str(tmp_path / "A.gff3"), None, None),
               Sample("B", tmp_path / "B.fa", None, None, None))
    _command(monkeypatch, None, lambda args: manifest.Finished(Path(args.outdir), samples, route="assembly",
                                                                 params={"ploidy": 2}))
    monkeypatch.chdir(tmp_path)
    assert cli.main(["boom", "-o", "out", "-t", "2"]) == 0
    doc = json.loads((out / "manifest.json").read_text())
    assert doc["command"] == ["synltr", "boom", "-o", "out", "-t", "2"]
    assert doc["params"] == {"command": "boom", "outdir": str(out), "threads": 2, "verbose": False, "ploidy": 2}
    assert doc["route"] == "assembly" and doc["code"] == {}  # a module outside the package has no fingerprint
    assert [s["label"] for s in doc["samples"]] == ["A", "B"] and doc["samples"][1]["annotation"] is None
    assert doc["samples"][0]["fasta"] == str(sheet.parent / "A.fa") and doc["samples"][1]["reads"] is None
    assert doc["tool"] == "synltr" and doc["version"] == __version__ and doc["started"] <= doc["finished"]
    assert set(TOOL_VERSION_ARGS) <= set(doc["dependencies"])
    (out / "manifest.json").unlink()
    for exc, result in ((None, 0), (SynltrError("bad input"), 0)):
        _command(monkeypatch, exc, result)
        cli.main(["boom", "-o", "out"])
        assert not (out / "manifest.json").exists()


def test_subgenomes_run_writes_its_manifest_and_plot_reads_it(subgenomes_run, capsys):
    """A finished `subgenomes` run records itself (its sheet's row, its options, no route); `synltr plot` reads
    which run it is from the manifest, redraws report.pdf and leaves the manifest as it was."""
    sheet, run, code = subgenomes_run
    assert code == 0
    doc = json.loads((run / "manifest.json").read_text())
    assert doc["params"]["command"] == "subgenomes" and doc["params"]["sheet"] == str(sheet)
    assert doc["params"]["seed"] == 1 and doc["params"]["n"] is None and doc["route"] is None
    assert [(s["label"], s["library"], s["checksums"]) for s in doc["samples"]] == [
        ("toy", str(sheet.parent / "in" / "ltr.fa"), {})]
    assert list(doc["code"]) == ["analyses.subgenomes.run"]
    before = (run / "manifest.json").read_bytes()
    (run / "report.pdf").unlink()
    assert cli.main(["plot", str(run)]) == 0 and (run / "report.pdf").is_file()
    assert (run / "manifest.json").read_bytes() == before


def test_simulations_write_no_manifest(tmp_path):
    """The simulators record their runs in their own files (sim nucleus: truth/history.json; sim pop: params.json),
    so `sim` writes no manifest.json."""
    sc = tmp_path / "scenario.json"
    sc.write_text(json.dumps({**SUBGENOMES_TOY, "families": {"n": 10}}))
    assert cli.main(["sim", "nucleus", str(sc), "-o", str(tmp_path / "lib"), "--truthdir", str(tmp_path / "truth")]) == 0
    assert (tmp_path / "truth" / "history.json").is_file()
    assert not list(tmp_path.rglob("manifest.json"))


@needs_lastz
def test_checks_and_pair_jobs_keep_the_finished_run_record(world, first, tmp_path, monkeypatch):
    """manifest.json describes the run that finished in its folder until another run there passes its checks: a
    --pair job (it finishes no run), a --paf that is no folder (exit 2) and a missing lastz (exit 1, genotype or
    geneflow) leave it as it was."""
    out = tmp_path / "out"
    shutil.copytree(first[0], out)
    before = (out / "manifest.json").read_bytes()
    assert cli.main(["genotype", str(world.sheet), "-o", str(out), "--pair", "A,B", "-t", "2"]) == 0
    assert (out / "manifest.json").read_bytes() == before
    assert cli.main(["genotype", str(world.sheet), "-o", str(out), "--paf", str(tmp_path / "none")]) == 2
    assert (out / "manifest.json").read_bytes() == before
    real = shutil.which
    monkeypatch.setattr(shutil, "which", lambda exe, *a, **k: None if exe == "lastz" else real(exe, *a, **k))
    for command in ("genotype", "geneflow"):
        assert cli.main([command, str(world.sheet), "-o", str(out)]) == 1
        assert (out / "manifest.json").read_bytes() == before, command


def test_failed_subgenomes_rerun_leaves_no_record_of_the_earlier_run(subgenomes_run, tmp_path, monkeypatch):
    """A subgenomes rerun in a finished run's folder that fails once its checks have passed leaves neither the
    earlier run's manifest.json nor its figure (report_data.json, report.png, report.pdf, report_numbers.tsv)."""
    from synltr.analyses.subgenomes import run as subgenomes

    sheet, done, _ = subgenomes_run
    out = tmp_path / "run"
    shutil.copytree(done, out)
    record = ("manifest.json", "report_data.json", "report.png", "report.pdf", "report_numbers.tsv")
    assert all((out / name).is_file() for name in record)

    def boom(*args, **kwargs):
        raise RuntimeError("a planted failure")

    monkeypatch.setattr(subgenomes, "load_or_build_clades", boom)
    assert cli.main(["subgenomes", str(sheet), "-o", str(out), "-t", "2"]) == 1
    assert not [name for name in record if (out / name).exists()]


@needs_lastz
def test_geneflow_records_its_route_and_inputs(world, first, tmp_path):
    """`geneflow` on the toy world (its pairs done by `genotype`): manifest.json names route assembly, the options
    and every input checksum the run computed, as `genotype`'s does."""
    out = tmp_path / "out"
    shutil.copytree(first[0], out)
    assert cli.main(["geneflow", str(world.sheet), "-o", str(out), "-t", "2"]) == 0
    doc = json.loads((out / "manifest.json").read_text())
    assert doc["route"] == "assembly" and doc["params"]["command"] == "geneflow"
    assert doc["params"]["sheet"] == str(world.sheet) and doc["params"]["identity"] is False
    assert list(doc["code"]) == ["analyses.geneflow.run"]
    genotyped = json.loads((first[0] / "manifest.json").read_text())
    assert [s["checksums"] for s in doc["samples"]] == [s["checksums"] for s in genotyped["samples"]]
    assert all((out / name).is_file() for name in FIGURE)


# ---------------------------------------------------------------- `synltr genotype` end to end

SEED = 7
CASES = ["P", "E", "S", "T", "D_left", "G", "nested_E"]
LABELS = ("A", "B", "C")
PAIRS = ("A__B", "A__C", "B__A", "B__C")  # (annotated A, any other B), sheet order
OUTPUTS = ("loci.tsv", "pav_matrix.tsv", "qc.tsv", "pav.vcf.gz", "pav.vcf.gz.tbi", "pav.unplaced.tsv",
           *(f"{x}.pav.{ext}" for x in LABELS for ext in ("gff3", "bed")), "manifest.json")


@dataclass
class World:
    """Three toy genomes in `dir` (C's FASTA in the read-only directory C/) and their samples sheet; the
    annotated elements as the run names them; the expected states in (A, B, C) of every locus, by its first
    member's eid; B's element that A lacks, and the 1-based A position of the target site (TSD) it lacks it at."""

    dir: Path
    sheet: Path
    elements: dict[str, list[Element]]
    expect: dict[str, tuple[str, str, str]]
    b_only: Element
    a_site: int


def element(label: str, name: str, seqid: str, start: int, length: int, ltr: int) -> Element:
    """An element of genome `label`, 1-based closed from `start`, with `ltr`-bp LTRs and TSD builders.T."""
    end = start + length - 1
    return Element(genome=label, eid=make_eid(label, seqid, start, end), name=name, seqid=seqid, start=start,
                   end=end, strand="+", ltr_left=(start, start + ltr - 1), ltr_right=(end - ltr + 1, end),
                   tsd=builders.T, k=0.0, family=None, superfamily=None, clade=None, confidence=None)


def make_world(d: Path, seqids: tuple[str, str, str] = ("seqA", "seqB", "seqC")) -> World:
    """A and B: builders.build_pair over CASES, then a segment where only B holds an element (A: L T R; B: L T e T
    R). A is annotated with build_pair's elements; B with its copies of the P case's element and of nested_E's host
    (without the child), and its own element. C (genotype-only) is B with 0.5% substitutions."""
    sid = dict(zip(LABELS, seqids, strict=True))
    rng = random.Random(SEED)
    pair = builders.build_pair(rng, CASES)
    draw = random.Random(SEED)
    built = [builders.case(name, draw) for name in CASES]  # build_pair draws its cases first, in order: these
    a_els = [replace(e, genome="A", eid=make_eid("A", sid["A"], e.start, e.end), seqid=sid["A"])
             for e in pair.elements]
    expect = {make_eid("A", sid["A"], e.start, e.end): ("P", pair.expect[e.eid], pair.expect[e.eid])
              for e in pair.elements}
    b_els = []
    off_a = off_b = builders.SPACER  # each case's 0-based start in A and in B
    for name, c in zip(CASES, built, strict=True):
        if name in ("P", "nested_E"):  # B holds the element (nested_E: the host, without its child)
            h = next(e for e in c.elements if e.level == 0)
            ltr = h.ltr_left[1] - h.ltr_left[0] + 1
            length = h.end - h.start + 1 - sum(e.end - e.start + 1 + len(builders.T2) for e in c.elements if e.level)
            be = element("B", f"{name}.b", sid["B"], off_b + h.start, length, ltr)
            b_els.append(be)
            a_seq, b_seq = pair.a[off_a + h.start - 1:off_a + h.end], pair.b[be.start - 1:be.end]
            if name == "P":  # B's copy carries 1% substitutions
                assert len(a_seq) == len(b_seq) and sum(map(str.__eq__, a_seq, b_seq)) > 0.95 * len(a_seq)
            else:  # B's host: A's LTRs at both ends
                assert b_seq[:ltr] == b_seq[-ltr:] == a_seq[:ltr]
        off_a += len(c.a) + builders.SPACER
        off_b += len(c.b) + builders.SPACER
    assert (off_a, off_b) == (len(pair.a), len(pair.b))  # build_pair's layout
    flank_l, flank_r = builders.rand_seq(builders.FLANK, rng), builders.rand_seq(builders.FLANK, rng)
    e, ltr, _ = builders.make_element(rng, 400, 4000)
    tail = builders.rand_seq(builders.SPACER, rng)
    a = pair.a + flank_l + builders.T + flank_r + tail
    b = pair.b + flank_l + builders.T + e + builders.T + flank_r + tail
    b_only = element("B", "b_only", sid["B"], len(pair.b) + builders.FLANK + len(builders.T) + 1, len(e), ltr)
    assert b[b_only.start - 1:b_only.end] == e
    b_els.append(b_only)
    expect[b_only.eid] = ("E", "P", "P")
    c = builders.mutate(b, 0.005, rng)

    builders.write_fasta(d / "A.fa", {sid["A"]: a})
    builders.write_fasta(d / "B.fa", {sid["B"]: b})
    builders.write_fasta(d / "C" / "C.fa", {sid["C"]: c})
    builders.write_gff3(d / "A.gff3", a_els)
    builders.write_gff3(d / "B.gff3", b_els)
    sheet = d / "samples.tsv"
    sheet.write_text("label\tfasta\tannotation\nA\tA.fa\tA.gff3\nB\tB.fa\tB.gff3\nC\tC/C.fa\t.\n")
    (d / "C").chmod(0o555)  # read-only: C's FASTA index must go to the run's work/fai
    return World(d, sheet, {"A": a_els, "B": b_els}, expect, b_only, len(pair.a) + builders.FLANK + 1)


def read_tsv(path: Path) -> list[dict[str, str]]:
    lines = path.read_text().splitlines()
    head = lines[0].split("\t")
    return [dict(zip(head, ln.split("\t"), strict=True)) for ln in lines[1:]]


def loci_by_first(out: Path) -> dict[str, dict[str, str]]:
    """loci.tsv rows by the locus's first member eid (sheet order)."""
    return {next(r[f"{x}_eid"] for x in LABELS if r[f"{x}_eid"] != "."): r for r in read_tsv(out / "loci.tsv")}


def qc_values(out: Path) -> dict[tuple[str, str], str]:
    return {(r["scope"], r["metric"]): r["value"] for r in read_tsv(out / "qc.tsv")}


def outputs(out: Path) -> dict[Path, bytes]:
    """Every file under `out` but manifest.json and the work directory: {relative path: bytes}."""
    return {p.relative_to(out): p.read_bytes() for p in sorted(out.rglob("*"))
            if p.is_file() and p.name != "manifest.json" and p.relative_to(out).parts[0] != "work"}


def vcf_text(out: Path) -> str:
    with gzip.open(out / "pav.vcf.gz", "rt") as fh:
        return fh.read()


def _sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    d = tmp_path_factory.mktemp("world")
    yield make_world(d)
    (d / "C").chmod(0o755)


@pytest.fixture(scope="module")
def first(world, tmp_path_factory):
    """`synltr genotype` on the world: (outdir, argv)."""
    if shutil.which("lastz") is None:
        pytest.skip("lastz (the window aligner) not on PATH")
    out = tmp_path_factory.mktemp("run") / "out"
    argv = ["genotype", str(world.sheet), "--outdir", str(out), "--threads", "2"]
    assert cli.main(argv) == 0
    logging.getLogger("synltr").handlers.clear()
    return out, argv


def _mini_sheet(d: Path) -> Path:
    """A (annotated: one element) and B, unrelated random 6-kb genomes, and their samples sheet."""
    rng = random.Random(1)
    builders.write_fasta(d / "A.fa", {"chr1": builders.rand_seq(6000, rng)})
    builders.write_fasta(d / "B.fa", {"chr1": builders.rand_seq(6000, rng)})
    builders.write_gff3(d / "A.gff3", [element("A", "e", "chr1", 1001, 1000, 100)])
    sheet = d / "samples.tsv"
    sheet.write_text("label\tfasta\tannotation\nA\tA.fa\tA.gff3\nB\tB.fa\t.\n")
    return sheet


@pytest.fixture
def mini(tmp_path):
    """A (annotated: one element) and B, unrelated random 6-kb genomes: enough for the checks made before any
    alignment, and for a run whose one pair aligns nowhere. Returns (directory, samples sheet)."""
    return tmp_path, _mini_sheet(tmp_path)


@pytest.fixture
def logged(monkeypatch, caplog):
    """cli.main's log records in caplog: log.setup only sets the level (DEBUG with -v, else INFO), and the records
    reach caplog by propagation, as many as the console would print."""
    def setup(verbose):
        logging.getLogger("synltr").setLevel(logging.DEBUG if verbose else logging.INFO)

    monkeypatch.setattr(cli.log, "setup", setup)
    caplog.set_level(logging.DEBUG)
    return caplog


@needs_lastz
def test_run_toy_three_genomes(world, first):
    out, argv = first
    for name in OUTPUTS:
        assert (out / name).is_file(), name
    for pair in PAIRS:
        assert (out / "pairs" / f"{pair}.calls.tsv.gz").is_file(), pair
    assert sorted(p.name for p in (out / "work" / "elements").glob("*.tsv")) == ["A.tsv", "B.tsv"]
    # C's directory is read-only: its index was built in the run's cache, not beside it
    assert len(list((out / "work" / "fai").glob("*.C.fa.fai"))) == 1
    assert not (world.dir / "C" / "C.fa.fai").exists()

    m = json.loads((out / "manifest.json").read_text())
    assert (m["version"], m["command"]) == (__version__, ["synltr", *argv])
    assert m["code"] == {"calls.genotype": checkpoint.code_fingerprint("calls.genotype")}
    assert m["dependencies"]["lastz"].startswith("lastz")
    assert all(m["dependencies"][x] for x in ("scipy", "matplotlib"))
    assert ("msprime" in m["dependencies"]) == (importlib.util.find_spec("msprime") is not None)
    assert set(TOOL_VERSION_ARGS) <= set(m["dependencies"])  # lastz, minimap2 and bcftools included
    assert m["route"] == "assembly" and m["params"]["command"] == "genotype" and m["params"]["pair"] is None
    by = {s["label"]: s for s in m["samples"]}
    for x in LABELS:
        assert by[x]["checksums"]["fasta"] == _sha256(by[x]["fasta"])
    for x in ("A", "B"):
        assert by[x]["checksums"]["annotation"] == _sha256(world.dir / f"{x}.gff3")
    assert set(by["C"]["checksums"]) == {"fasta"}
    assert m["started"] <= m["finished"]

    assert {eid: tuple(r[f"{x}_state"] for x in LABELS) for eid, r in loci_by_first(out).items()} == world.expect

    qc = qc_values(out)
    assert qc[("A~B", "reciprocity")] == "1"
    assert qc[("A>B", "annot_agree")] == "1" and qc[("B>A", "annot_agree")] == "1"
    assert qc[("A>B", "e_in_annot")] == "0"  # nested_E's child: E inside its host's ortholog, not counted
    assert qc[("B", "gap_rate")] == "0.125"  # the G case: 1 of A's 8 calls in B
    assert ("A>C", "annot_agree") not in qc and ("A~C", "reciprocity") not in qc  # C has no annotation


@needs_lastz
def test_run_shared_seqids(first, tmp_path):
    """All three genomes name their sequence chr1. The outputs are the distinct-name run's
    with the names replaced, and each genome's coordinates are its own."""
    d = tmp_path / "world"
    d.mkdir()
    world = make_world(d, seqids=("chr1", "chr1", "chr1"))
    out = tmp_path / "out"
    try:
        assert cli.main(["genotype", str(world.sheet), "--outdir", str(out), "--threads", "2"]) == 0
    finally:
        (d / "C").chmod(0o755)

    def rename(text: str) -> str:
        for name in ("seqA", "seqB", "seqC"):
            text = text.replace(name, "chr1")
        return text

    ref = first[0]
    for name in OUTPUTS:
        if name.endswith((".tsv", ".gff3", ".bed")):
            assert (out / name).read_text() == rename((ref / name).read_text()), name
    assert vcf_text(out) == rename(vcf_text(ref))

    loci = loci_by_first(out)
    a_host = next(e for e in world.elements["A"] if e.name == "nested_E.h")
    b_host = next(e for e in world.elements["B"] if e.name == "nested_E.b")
    assert a_host.start != b_host.start  # the two genomes' coordinates differ here
    host = loci[a_host.eid]
    assert [host[k] for k in ("A_start", "A_end", "B_eid", "B_start", "B_end")] == [
        str(a_host.start), str(a_host.end), b_host.eid, str(b_host.start), str(b_host.end)]
    empty = loci[world.b_only.eid]
    assert [empty[k] for k in ("A_state", "A_seqid", "A_start", "A_end", "B_start", "B_end")] == [
        "E", "chr1", str(world.a_site), str(world.a_site + 4), str(world.b_only.start), str(world.b_only.end)]
    with pysam.VariantFile(str(out / "pav.vcf.gz")) as vcf:
        recs = {r.id: (r.contig, r.pos, [r.samples[x]["GT"] for x in LABELS]) for r in vcf}
    # the VCF is on A's coordinates: the host at its start - 1; the empty site at its TSD copy's last base
    assert recs[host["lid"]] == ("chr1", a_host.start - 1, [(1,), (1,), (1,)])
    assert recs[empty["lid"]] == ("chr1", world.a_site + 4, [(0,), (1,), (1,)])


@needs_lastz
def test_rerun_skips(first, logged):
    out, argv = first
    before = outputs(out)
    assert cli.main(argv) == 0
    text = logged.text
    assert text.count("elements done, skip") == 2  # A and B
    assert text.count("anchors done, skip") == len(PAIRS)
    assert text.count("calls done, skip") == len(PAIRS)
    assert "indexing" not in text
    assert outputs(out) == before  # every output but manifest.json, byte for byte


@needs_lastz
def test_run_warns_on_unrelated_pair(mini, logged):
    d, sheet = mini
    assert cli.main(["genotype", str(sheet), "--outdir", str(d / "out"), "--threads", "1"]) == 0
    warn = [r.getMessage() for r in logged.records
            if r.name == "synltr.calls.genotype" and r.levelno == logging.WARNING]
    assert len(warn) == 1 and "no orthologous alignment between A and B?" in warn[0]
    assert (d / "out" / "manifest.json").is_file()


@needs_lastz
def test_pair_jobs_then_a_run_in_a_moved_folder(world, first, tmp_path, logged):
    """A run folder stays valid when moved. There, `genotype --pair B,A` finds its pair done (anchors and calls) and
    leaves manifest.json as it was (a pair job finishes no run); the run without --pair skips every pair, writes the
    same outputs and records itself."""
    out = tmp_path / "out"
    shutil.copytree(first[0], out)
    before = outputs(out)
    recorded = (out / "manifest.json").read_bytes()
    base = [str(world.sheet), "--outdir", str(out), "--threads", "2"]
    assert cli.main(["genotype", *base, "--pair", "B,A"]) == 0
    assert (out / "manifest.json").read_bytes() == recorded
    assert cli.main(["genotype", *base]) == 0
    text = logged.text
    assert text.count("anchors done, skip") == 1 + len(PAIRS)  # B,A; then the run checks every pair's
    assert text.count("calls done, skip") == 1 + len(PAIRS)
    assert outputs(out) == before
    assert json.loads((out / "manifest.json").read_text())["command"] == ["synltr", "genotype", *base]


@needs_lastz
def test_outdir_of_older_code_is_recomputed_stage_by_stage(mini, logged, monkeypatch, tmp_path):
    """Each stage key holds the code fingerprint of the stage's module and its imports, here read from an edited
    copy of the package.
    - The same code: every stage is skipped.
    - A comment in calls/vcf.py (the VCF writer, outside every stage): still skipped.
    - An edit of routes/assembly/classify.py: a rerun recomputes the calls only.
    - An edit of base/constants.py: every stage reruns.
    - Another annotation of A: its elements, the anchors (mapped from its elements' flanks) and the calls are
      recomputed."""
    d, sheet = mini
    argv = ["genotype", str(sheet), "--outdir", str(d / "out"), "--threads", "1"]
    skips = ("elements done, skip", "anchors done, skip", "calls done, skip")
    pair = d / "out" / "pairs" / "A__B.calls.tsv.gz"
    assert cli.main(argv) == 0
    code = checkpoint.done_meta(checkpoint.marker(pair))["code"]
    assert code == checkpoint.code_fingerprint("routes.assembly.pair")
    copy = tmp_path / "synltr"
    shutil.copytree(Path(checkpoint.__file__).parents[1], copy, ignore=shutil.ignore_patterns("__pycache__"))
    monkeypatch.setattr(checkpoint, "_ROOT", copy)
    checkpoint._fingerprint.cache_clear()

    def edit(module):
        (copy / module).write_text((copy / module).read_text() + "\n# an edit\n")
        checkpoint._fingerprint.cache_clear()

    def rerun():
        logged.clear()
        assert cli.main(argv) == 0
        return [s for s in skips if s in logged.text]

    try:
        assert rerun() == list(skips)  # the copy's bytes: the same keys
        edit("calls/vcf.py")
        assert rerun() == list(skips)
        edit("routes/assembly/classify.py")
        assert rerun() == ["elements done, skip", "anchors done, skip"]
        edit("base/constants.py")
        assert rerun() == []  # every stage recomputed
        builders.write_gff3(d / "A.gff3", [element("A", "e", "chr1", 2001, 1000, 100)])  # another annotation of A
        assert rerun() == []
    finally:
        checkpoint._fingerprint.cache_clear()


@needs_lastz
def test_elements_checkpoint_follows_the_annotation(mini, logged):
    d, sheet = mini
    argv = ["genotype", str(sheet), "--outdir", str(d / "out"), "--pair", "A,B", "--threads", "1"]
    for _ in range(2):
        assert cli.main(argv) == 0
    builders.write_gff3(d / "A.gff3", [element("A", "e", "chr1", 2001, 1000, 100)])  # the annotation changed
    logged.clear()  # only the third run's records
    assert cli.main(argv) == 0
    assert "A: 1 elements" in logged.text and "elements done, skip" not in logged.text
    assert [e.start for e in read_elements(d / "out" / "work" / "elements" / "A.tsv")] == [2001]


def test_checksums_name_each_input(tmp_path):
    """The manifest names each input of a sample by what it is: "fasta"; "annotation" for a GFF3 or BED file, or
    "annotation:<name>" for each table of an LTRquest prefix; "library:<name>" for each library file. It records the
    checksums the run's stages keyed on (checkpoint.keyed_checksums), of the kinds of input its route reads: an input
    no stage keyed on is left out, and nothing is hashed for the manifest."""
    fa = builders.write_fasta(tmp_path / "g.fa", {"chr1": "ACGT" * 10})
    tables = [tmp_path / f"g_depth{n}_clean_ltr.tsv" for n in (0, 1)]
    for n, t in enumerate(tables):
        t.write_text(f"{n}\n")
    gff = tmp_path / "g.gff3"
    gff.write_text("##gff-version 3\n")
    lib = tmp_path / "lib.fa"
    lib.write_text(">x\nACGT\n")
    with checkpoint.keyed_checksums() as keyed:
        for f in (fa, *tables, gff):
            checkpoint.file_checksum(f, tmp_path / "work" / "checksums.json")
    samples = [Sample("G", fa, str(tmp_path / "g"), None, None), Sample("H", fa, str(gff), None, None),
               Sample("I", fa, None, None, None, library=str(lib))]
    got = manifest.sample_checksums(samples, keyed, ("fasta", "annotation", "library"))
    assert got == {"G": {"fasta": _sha256(fa), **{f"annotation:{t.name}": _sha256(t) for t in tables}},
                   "H": {"fasta": _sha256(fa), "annotation": _sha256(gff)},
                   "I": {"fasta": _sha256(fa)}}  # no stage keyed on lib.fa
    assert manifest.sample_checksums(samples, keyed, ("annotation",)) == {
        "G": {f"annotation:{t.name}": _sha256(t) for t in tables}, "H": {"annotation": _sha256(gff)}, "I": {}}
    keyed[str(lib.resolve())] = "the checksum a stage keyed on"  # read, never computed
    assert manifest.sample_checksums(samples[2:], keyed, ("library",)) == {
        "I": {"library:lib.fa": "the checksum a stage keyed on"}}


@needs_lastz
def test_manifest_records_the_checksums_the_stages_keyed_on(world, first, tmp_path, monkeypatch):
    """manifest.json records each input's checksum as the run's stages keyed on it, never what the run's checksum
    cache holds at the end: C's FASTA rewritten after its last stage is recorded as the stages read it (not left
    out), and B's annotation rewritten and hashed into the run's cache by another process (a --pair job of the same
    folder) is recorded as the stages read it (not as the other writer's checksum)."""
    from synltr.calls import genotype as command

    inputs, out = tmp_path / "in", tmp_path / "out"
    shutil.copytree(world.dir, inputs)
    (inputs / "C").chmod(0o755)
    shutil.copytree(first[0], out)  # its pairs are done: the rerun keys every stage and skips it
    c_fa, b_gff = inputs / "C" / "C.fa", inputs / "B.gff3"
    read = {"C": _sha256(c_fa), "B": _sha256(b_gff)}
    real = command.write_outputs

    def change_inputs_then_write(samples, outdir, genomes):  # every stage has keyed on its inputs by now
        text = c_fa.read_text()
        at = text.index("\n") + 1
        c_fa.write_text(text[:at] + ("C" if text[at] == "A" else "A") + text[at + 1:])
        b_gff.write_text(b_gff.read_text() + "# edited\n")
        code = "import sys; from synltr.base.checkpoint import file_checksum; file_checksum(sys.argv[1], sys.argv[2])"
        subprocess.run([sys.executable, "-c", code, str(b_gff), str(Path(outdir) / "work" / "checksums.json")],
                       check=True, env=ENV)
        return real(samples, outdir, genomes)

    monkeypatch.setattr(command, "write_outputs", change_inputs_then_write)
    assert cli.main(["genotype", str(inputs / world.sheet.name), "-o", str(out), "-t", "2"]) == 0
    cache = json.loads((out / "work" / "checksums.json").read_text())
    assert cache[str(b_gff.resolve())]["sha256"] == _sha256(b_gff) != read["B"]  # the other writer's entry
    by = {s["label"]: s["checksums"] for s in json.loads((out / "manifest.json").read_text())["samples"]}
    assert (by["C"], by["B"].get("annotation")) == ({"fasta": read["C"]}, read["B"]) and read["C"] != _sha256(c_fa)


def _changed_during_the_run(caplog, path: Path) -> list[str]:
    """The warnings that `path` changed between two stages of the run (base.checkpoint.keyed_checksums)."""
    return [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING
            and r.getMessage().startswith(f"{path.resolve()} changed during the run")]


@needs_lastz
def test_geneflow_records_its_inputs_as_its_analyses_keyed_them(tmp_path, monkeypatch, logged):
    """geneflow's manifest.json records the checksums its stages keyed on in both of its phases, the analyses'
    included: C's FASTA, edited as the analyses start (the genotyping has keyed it), is recorded as the identity layer
    keyed it when it realigned the trio's copies, with the warning that it changed during the run; the PAF checksums
    the anchors keyed on stay beside --paf."""
    from analyses.geneflow.test_run import _copies_world

    from synltr.analyses.geneflow import run as gf

    sheet, paf = _copies_world(tmp_path)
    c_fa = tmp_path / "C.fa"
    genotyped = _sha256(c_fa)
    real = gf.analyses

    def edit_then_analyse(*args, **kwargs):
        text = c_fa.read_text()
        at = text.index("\n") + 1
        c_fa.write_text(text[:at] + ("C" if text[at] == "A" else "A") + text[at + 1:])
        return real(*args, **kwargs)

    monkeypatch.setattr(gf, "analyses", edit_then_analyse)
    out = tmp_path / "out"
    argv = ["geneflow", str(sheet), "-o", str(out), "--paf", str(paf), "--identity", "--trio", "B,C,A", "-t", "1"]
    assert cli.main(argv) == 0
    assert (out / "identity.tsv").is_file()
    doc = json.loads((out / "manifest.json").read_text())
    by = {s["label"]: s["checksums"] for s in doc["samples"]}
    assert by["C"]["fasta"] == _sha256(c_fa) != genotyped
    assert len(_changed_during_the_run(logged, c_fa)) == 1
    assert doc["params"]["paf_checksums"] == {n: _sha256(paf / n) for n in ("A__B.paf", "A__C.paf")}


def test_screen_records_its_libraries_as_its_trio_stage_keyed_them(tmp_path, monkeypatch, logged):
    """The screen's manifest.json records the checksums its stages keyed on in both of its phases, the trio stage's
    included: R's library FASTA, edited as the analyses start (the pairs have keyed it), is recorded as the trio stage
    keyed it, with the warning that it changed during the run; the other files of the libraries as both phases keyed
    them."""
    from synltr.analyses.geneflow import run as gf

    sc = builders.trio(tmp_path / "in", seed=3, admix=True)
    r_fa = Path(sc.libraries["R"].split(",")[0])
    screened = _sha256(r_fa)
    real = gf.screen_analyses

    def edit_then_analyse(*args, **kwargs):
        with open(r_fa, "a") as fh:
            fh.write("\n")  # a blank line: the same records, another checksum
        return real(*args, **kwargs)

    monkeypatch.setattr(gf, "screen_analyses", edit_then_analyse)
    out = tmp_path / "out"
    assert cli.main(["geneflow", "--screen", str(sc.sheet), "-o", str(out), "-t", "2"]) == 0
    by = {s["label"]: s["checksums"] for s in json.loads((out / "manifest.json").read_text())["samples"]}
    assert by["R"][f"library:{r_fa.name}"] == _sha256(r_fa) != screened
    assert len(_changed_during_the_run(logged, r_fa)) == 1
    files = {g: [Path(f) for f in sc.libraries[g].split(",")] for g in "RSD"}
    assert {g: by[g] for g in "SD"} == {g: {f"library:{f.name}": _sha256(f) for f in files[g]} for g in "SD"}


@needs_lastz
def test_paf_run_records_its_alignments(tmp_path):
    """A --paf run records, beside the PAF folder among its options, the checksum of each whole-genome alignment its
    anchors keyed on, by its name in the folder (a symlinked PAF too; another file of the folder is not one)."""
    from analyses.geneflow.test_run import _copies_world

    sheet, paf = _copies_world(tmp_path)
    kept = tmp_path / "kept_elsewhere.paf"
    (paf / "A__C.paf").rename(kept)
    (paf / "A__C.paf").symlink_to(kept)
    (paf / "notes.txt").write_text("not an alignment\n")
    out = tmp_path / "out"
    assert cli.main(["genotype", str(sheet), "-o", str(out), "--paf", str(paf), "-t", "1"]) == 0
    params = json.loads((out / "manifest.json").read_text())["params"]
    assert params["paf"] == str(paf)
    assert params["paf_checksums"] == {"A__B.paf": _sha256(paf / "A__B.paf"), "A__C.paf": _sha256(kept)}


@needs_lastz
def test_unfinished_pair_file_is_recomputed(world, first, tmp_path, logged):
    """A pair file without its done marker (left by a killed run, or from other inputs) is never built into the loci:
    the run genotypes that pair again, and the outputs are those of the uninterrupted run."""
    out = tmp_path / "out"
    shutil.copytree(first[0], out)

    def written():  # the outputs and pair files; the done markers record times
        return {p: data for p, data in outputs(out).items() if not p.name.endswith(".done.json")}

    before = written()
    stale = out / "pairs" / "A__C.calls.tsv.gz"
    checkpoint.marker(stale).unlink()
    stale.write_bytes(gzip.compress(b"stale\n"))
    assert cli.main(["genotype", str(world.sheet), "--outdir", str(out), "--threads", "2"]) == 0
    assert logged.text.count("calls done, skip") == len(PAIRS) - 1
    assert written() == before and checkpoint.marker(stale).is_file()


# `synltr genotype` with one window per chunk, so each pair's windows make a pool of several chunks; with "slow",
# every chunk a worker takes announces itself (busy.<pid> in ROOT) and waits. argv: slow|fast ROOT ARGV...
_ONE_WINDOW_CHUNKS = r"""
import os, sys, time
from pathlib import Path

import synltr.routes.assembly.pair as genotyping
from synltr import cli

mode, root, argv = sys.argv[1], Path(sys.argv[2]), sys.argv[3:]
genotyping.CHUNK_MIN = 1
if mode == "slow":
    real, parent = genotyping._Genotyper.rows, os.getpid()

    def rows(self, win):
        if os.getpid() != parent:
            (root / f"busy.{os.getpid()}").touch()
            time.sleep(120)
        return real(self, win)

    genotyping._Genotyper.rows = rows
sys.exit(cli.main(argv))
"""


@needs_lastz
@pytest.mark.skipif(not sys.platform.startswith("linux"), reason="Ctrl-C as a terminal sends it: a process group")
def test_ctrl_c_is_one_line_and_a_rerun_resumes(world, first, tmp_path):
    """Ctrl-C (SIGINT to the process group, as a terminal sends it) while a genotype run's pool genotypes chunks: the
    workers end without a word, the run adds one line to its log, `error: interrupted`, exits 130 and leaves no
    process behind; the same command again resumes the run (its finished stages skipped) and writes what an
    uninterrupted run writes."""
    script = tmp_path / "run.py"
    script.write_text(_ONE_WINDOW_CHUNKS)
    argv = ["genotype", str(world.sheet), "-o", str(tmp_path / "out"), "-t", "2"]
    proc = subprocess.Popen([sys.executable, str(script), "slow", str(tmp_path), *argv], stdout=subprocess.DEVNULL,
                            stderr=subprocess.PIPE, text=True, env=ENV, start_new_session=True)
    try:
        t0 = time.monotonic()
        while len(list(tmp_path.glob("busy.*"))) < 2:  # both workers inside a chunk
            assert proc.poll() is None and time.monotonic() - t0 < 120, proc.communicate()[1]
            time.sleep(0.05)
        os.killpg(proc.pid, signal.SIGINT)
        err = proc.communicate(timeout=60)[1]
        t0 = time.monotonic()
        while time.monotonic() - t0 < 10:  # every process of the group gone (a worker could still be reaped)
            try:
                os.killpg(proc.pid, 0)
            except ProcessLookupError:
                break
            time.sleep(0.05)
        else:
            raise AssertionError("a process of the interrupted run is still alive")
    finally:
        with contextlib.suppress(ProcessLookupError):
            os.killpg(proc.pid, signal.SIGKILL)
        proc.wait()
    log_lines = ("genotype: ", "A: ", "B: ", "C: ", "pair ", "A vs ", "B vs ")  # the run's own milestones
    assert proc.returncode == 130 and err.splitlines()[-1] == "error: interrupted", err
    assert [ln for ln in err.splitlines() if not ln.startswith(log_lines)] == ["error: interrupted"], err
    again = subprocess.run([sys.executable, str(script), "fast", str(tmp_path), *argv], capture_output=True,
                           text=True, env=ENV, check=False)
    assert again.returncode == 0 and "A: elements done, skip" in again.stderr, again.stderr
    written = [{p: data for p, data in outputs(out).items() if not p.name.endswith(".done.json")}  # markers hold times
               for out in (tmp_path / "out", first[0])]
    assert written[0] == written[1]


# ---------------------------------------------------------------- errors and help


@pytest.mark.parametrize("argv, message", [
    (["--pair", "A,Z"], "'Z' is not a label"),
    (["--pair", "B,A"], "B has no annotation"),
])
def test_fail_fast_on_labels(mini, capsys, argv, message):
    d, sheet = mini
    assert cli.main(["genotype", str(sheet), "--outdir", str(d / "out"), *argv]) == 2
    err = capsys.readouterr().err
    assert err.startswith("error: ") and message in err and err.count("\n") == 1, err  # one line
    assert not (d / "out").exists()  # nothing written


@pytest.mark.parametrize("argv", [[], ["--pair", "A,B"]], ids=["run", "pair"])
def test_fail_fast_on_missing_binary(mini, tmp_path, capsys, monkeypatch, argv):
    """PATH holds stubs of every executable the engine runs but lastz (the test sets it: the caller's PATH does not
    matter): one error line naming only lastz, before anything is written."""
    d, sheet = mini
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    for exe in {*TOOL_VERSION_ARGS, "samtools"} - {"lastz"}:
        (bin_dir / exe).write_text("#!/bin/sh\nexit 1\n")
        (bin_dir / exe).chmod(0o755)
    monkeypatch.setenv("PATH", str(bin_dir))
    assert cli.main(["genotype", str(sheet), "--outdir", str(d / "out"), *argv]) == 1
    err = capsys.readouterr().err
    assert err.startswith("error: ") and err.count("\n") == 1, err
    assert "needs lastz on PATH" in err and "minimap2" not in err, err
    assert not (d / "out").exists()  # nothing written


@pytest.fixture
def lastz_stub(tmp_path, monkeypatch):
    """Put a stub `lastz` (it prints a version line) first on PATH, so the tests that check inputs reach those
    checks on a host without lastz. Returns the real lastz's path, or None."""
    real = shutil.which("lastz")
    bin_dir = tmp_path / "lastz_stub"
    bin_dir.mkdir()
    (bin_dir / "lastz").write_text("#!/bin/sh\necho 'LASTZ version 1.04.00 (stub)'\n")
    (bin_dir / "lastz").chmod(0o755)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ.get('PATH', '')}")
    return real


def test_fail_fast_on_bad_inputs(mini, capsys, lastz_stub):
    d, sheet = mini
    bad = d / "bad.tsv"
    bad.write_text(sheet.read_text().replace("B.fa", "missing.fa"))
    assert cli.main(["genotype", str(bad), "--outdir", str(d / "out")]) == 2  # the sheet: a usage error
    assert "FASTA file(s) not found" in capsys.readouterr().err
    gff = d / "A.gff3"
    gff.write_text(gff.read_text().replace("chr1\t", "chrX\t"))
    assert cli.main(["genotype", str(sheet), "--outdir", str(d / "out")]) == 2  # the annotation: a usage error
    err = capsys.readouterr().err.splitlines()
    assert err[-1].startswith("error: ") and "absent from the FASTA" in err[-1]
    assert not any("Traceback" in ln for ln in err)


@needs_lastz
def test_malformed_annotation_records(mini, capsys, logged):
    """A GFF3 whose gene line has a bad coordinate runs (the line is never read), a malformed LTR-RT line is warned
    with its file and line and counted in qc.tsv and manifest.json; an LTRquest table without ltr3_end is a usage
    error: one error line naming the file and the column, no traceback."""
    d, sheet = mini
    gff = d / "A.gff3"
    text = gff.read_text()
    gff.write_text(text + "chr1\tx\tgene\t2000\t2O100\t.\t+\t.\tID=g1\n"
                   "chr1\tx\tLTR_retrotransposon\t3000\t35OO\t.\t+\t.\tID=bad\n")
    assert cli.main(["genotype", str(sheet), "--outdir", str(d / "out"), "--threads", "1"]) == 0
    warned = [r.getMessage() for r in logged.records if r.name.startswith("synltr.io.annotation")
              and r.levelno == logging.WARNING]
    assert warned[0].startswith(f"{gff}:{len(text.splitlines()) + 2}: record skipped, malformed_coordinates")
    assert qc_values(d / "out")[("A", "skipped_malformed_coordinates")] == "1"
    by = {s["label"]: s for s in json.loads((d / "out" / "manifest.json").read_text())["samples"]}
    assert by["A"]["skipped"] == {"malformed_coordinates": 1} and by["B"]["skipped"] == {}

    header = "seq_id\tseq_len\tstatus\tltr5_start\tltr5_end\tltr3_start\ttsd\tk2p\tstrand"
    (d / "lq_depth0_clean_ltr.tsv").write_text(f"#{header}\nchr1:1001-2000#LTR/Copia/SIRE\t1000\tpass\t1\t100\t901"
                                               "\tCACTC\t0\t+\n")
    sheet.write_text(sheet.read_text().replace("A.gff3", "lq"))
    capsys.readouterr()
    assert cli.main(["genotype", str(sheet), "--outdir", str(d / "out2")]) == 2
    err = capsys.readouterr().err.splitlines()
    assert err[-1].startswith("error: ") and "lacks required column(s) ltr3_end" in err[-1]
    assert "lq_depth0_clean_ltr.tsv" in err[-1] and not any("Traceback" in ln for ln in err)


def test_verbose_error_has_traceback(mini, capsys):
    d, sheet = mini
    assert cli.main(["genotype", str(sheet), "--outdir", str(d / "out"), "--pair", "A,Z", "-v"]) == 2
    err = capsys.readouterr().err
    assert "Traceback" in err and [ln for ln in err.splitlines() if ln.startswith("error: ")] == [
        "error: --pair A,Z: 'Z' is not a label of the samples sheet (A, B)"]


@pytest.mark.parametrize("argv", [
    [], ["genotype", "SHEET", "-o", "OUT", "--threads", "0"], ["genotype", "SHEET", "-o", "OUT", "--anchor", "bogus"],
    ["genotype", "SHEET", "-o", "OUT", "--pair", "A"], ["genotype", "SHEET", "-o", "OUT", "--pair", "A,A"],
    ["plot", "OUT", "--outdir", "OUT"],
], ids=["no-command", "threads-0", "anchor", "pair-one-label", "pair-twice", "plot-outdir"])
def test_usage_errors_exit_2(mini, argv):
    """Usage errors exit 2 before anything is written: no command, -t 0, a dropped option, a --pair that is not two
    different labels, an option the command does not take (plot has no -o: its run folder is its output)."""
    d, sheet = mini
    assert cli.main([{"SHEET": str(sheet), "OUT": str(d / "out")}.get(a, a) for a in argv]) == 2
    assert not (d / "out").exists()


@pytest.mark.parametrize("command, program", [("genotype", "lastz"), ("geneflow", "lastz"), ("source", "minimap2")])
def test_missing_program_is_an_error_not_a_usage_error(mini, tmp_path, capsys, monkeypatch, command, program):
    """A command whose external program is not on PATH (lastz for genotype and geneflow on assemblies, minimap2 for
    source) stops before any work with exit 1, one `error:` line naming it and the install route; the output
    folder's earlier manifest.json stays as it was."""
    d, sheet = mini
    out = tmp_path / "out"
    out.mkdir()
    (out / "manifest.json").write_text('{"params": {"command": "' + command + '"}}\n')  # an earlier run's record
    (tmp_path / "bin").mkdir()
    monkeypatch.setenv("PATH", str(tmp_path / "bin"))
    argv = ["source", str(d), "-o", str(out)] if command == "source" else [command, str(sheet), "-o", str(out)]
    assert cli.main(argv) == 1
    err = capsys.readouterr().err.strip().splitlines()
    assert len(err) == 1 and err[0].startswith("error: ") and program in err[0] and "environment.yml" in err[0], err
    assert sorted(p.name for p in out.iterdir()) == ["manifest.json"]
    assert (out / "manifest.json").read_text() == '{"params": {"command": "' + command + '"}}\n'


def test_missing_lastz_names_install_route(mini, capsys, monkeypatch):
    """With lastz hidden from the lookup the command line uses, `synltr genotype` stops before any work with one
    error line that names lastz and the conda environment file."""
    d, sheet = mini
    real = shutil.which
    monkeypatch.setattr(shutil, "which", lambda exe, *a, **k: None if exe == "lastz" else real(exe, *a, **k))
    assert cli.main(["genotype", str(sheet), "--outdir", str(d / "out")]) == 1
    err = capsys.readouterr().err
    assert err.startswith("error: ") and err.count("\n") == 1, err
    assert "lastz" in err and "environment.yml" in err, err
    assert not (d / "out").exists()


# per-pair anchoring detail
DETAIL = ("anchoring", "raw blocks", "pieces", "syntenic net", "reusing index", "indexed in")


@needs_lastz
def test_default_log_is_minimal(world, tmp_path, logged):
    """At default verbosity the log has the pair counter and "done" but no per-pair anchoring detail; -v adds
    it."""
    texts = {}
    for name, extra in (("quiet", []), ("verbose", ["-v"])):
        logged.clear()
        argv = ["genotype", str(world.sheet), "--outdir", str(tmp_path / name), "--threads", "2", *extra]
        assert cli.main(argv) == 0
        texts[name] = [r.getMessage() for r in logged.records if r.name.startswith("synltr")]
    quiet = "\n".join(texts["quiet"])
    assert "pair 1/" in quiet and "done in" in quiet
    assert not [d for d in DETAIL if d in quiet], quiet
    verbose = "\n".join(texts["verbose"])
    assert all(d in verbose for d in ("raw blocks", "syntenic net: kept")), verbose


@pytest.mark.parametrize("command, blocked", [("subgenomes", "scipy"), ("geneflow", "scipy")])
def test_missing_python_module_is_one_error_line(mini, command, blocked):
    """A Python module that cannot be imported gives one `error:` line naming it and the install route, exit 1, no
    traceback; `synltr geneflow` finds out before it starts genotyping."""
    d, sheet = mini
    args = [command, str(sheet), "--outdir", str(d / "out")]
    code = f"import sys; sys.modules[{blocked!r}] = None; from synltr.cli import main; sys.exit(main({args!r}))"
    res = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, cwd=d, check=False, env=ENV)
    assert res.returncode == 1, res.stderr
    lines = [ln for ln in res.stderr.splitlines() if ln.startswith("error:")]
    assert len(lines) == 1 and blocked in lines[0] and "environment.yml" in lines[0], res.stderr
    assert "Traceback" not in res.stderr and "pair(s)" not in res.stderr, res.stderr  # no "N genome(s), N pair(s)"
    assert not (d / "out").exists()  # nothing was started
