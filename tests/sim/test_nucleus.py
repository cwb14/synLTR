"""`ltrtrace sim nucleus` (ltrtrace.sim.nucleus and its subcommand): a tiny seeded scenario (three leaves, two
nested mergers, exchanges, gene conversion, birth divergence) whose outputs are pinned byte for byte in
tests/golden/nucleus_tiny/ (annotation.gff3's source column without the version it names), and the `source` column of
copies.tsv."""
from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest
from builders import unversioned_gff3
from sim.poptoy import sim_main

from ltrtrace import __version__

GOLDEN = Path(__file__).resolve().parent.parent / "golden" / "nucleus_tiny"
SEED = 7
LIB = ("ltr.fa", "k2p.tsv", "annotation.gff3", "genome.fai")
TRUTH = ("copies.tsv", "exchanges.tsv", "chromosomes.tsv", "orthologs.tsv", "paralogs.tsv", "history.json")


def scenario() -> dict:
    return json.loads((GOLDEN / "scenario.json").read_text())


def run_cli(tmp_path: Path, sc: dict, seed: int = SEED) -> tuple[Path, Path]:
    f = tmp_path / "scenario.json"
    f.write_text(json.dumps(sc))
    lib, truth = tmp_path / "lib", tmp_path / "truth"
    assert sim_main(["nucleus", str(f), "-o", str(lib), "--truthdir", str(truth), "--seed", str(seed)]) == 0
    return lib, truth


def rows(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@pytest.fixture(scope="module")
def tiny(tmp_path_factory):
    return run_cli(tmp_path_factory.mktemp("tiny"), scenario())


def without_source(path: Path) -> bytes:
    """copies.tsv as a table without its last column (`source`), tab- and line-exact."""
    lines = path.read_bytes().split(b"\n")
    assert lines[0].endswith(b"\tk2p\tsource")
    return b"\n".join(ln.rsplit(b"\t", 1)[0] for ln in lines[:-1]) + b"\n"


@pytest.mark.parametrize("name", LIB)
def test_library_is_pinned(tiny, name):
    pin = unversioned_gff3 if name.endswith(".gff3") else bytes
    assert pin((tiny[0] / name).read_bytes()) == pin((GOLDEN / "lib" / name).read_bytes())


def test_annotation_names_ltrtrace(tiny):
    """annotation.gff3's source column names the package and its version on every feature line."""
    lines = [ln for ln in (tiny[0] / "annotation.gff3").read_text().splitlines() if ln[:1] != "#"]
    assert lines and {ln.split("\t")[1] for ln in lines} == {f"ltrtrace-{__version__}"}


@pytest.mark.parametrize("name", [t for t in TRUTH if t != "copies.tsv"])
def test_truth_is_pinned(tiny, name):
    assert (tiny[1] / name).read_bytes() == (GOLDEN / "truth" / name).read_bytes()


def test_copies_are_pinned_without_source(tiny):
    assert without_source(tiny[1] / "copies.tsv") == (GOLDEN / "truth" / "copies.tsv").read_bytes()


def test_golden_is_not_trivial():
    """The pinned run exercises exchanges, both source kinds, post-merger copies and every leaf."""
    copies = rows(GOLDEN / "truth" / "copies.tsv")
    assert len(copies) > 100
    assert {r["subgenome"] for r in copies} == {"A", "C", "D"}
    assert {r["mode"] for r in copies} == {"own", "foreign", "ancestral"}
    assert any(r["exchange_id"] != "-1" for r in copies)
    assert len(rows(GOLDEN / "truth" / "exchanges.tsv")) >= 5


def test_source_column_follows_the_insertion(tiny):
    """A copy's source is the root of the copy it was transposed from: older than the copy's own root, of the
    same family when that copy is a row, and shared by every copy of one insertion; "." is a master."""
    copies = rows(tiny[1] / "copies.tsv")
    assert list(copies[0])[-1] == "source"
    by_root, src_of = {}, {}
    for r in copies:
        by_root.setdefault(r["root"], r)
        assert src_of.setdefault(r["root"], r["source"]) == r["source"]
    from_copy = [r for r in copies if r["source"] != "."]
    assert from_copy and len(from_copy) < len(copies)
    for r in from_copy:
        assert int(r["source"]) < int(r["root"])
        if r["source"] in by_root:
            assert by_root[r["source"]]["family"] == r["family"]


@pytest.mark.parametrize(("share", "expect"), [(0.0, 0), (1.0, 1)])
def test_source_follows_p_copy_source(tmp_path, share, expect):
    """No copy-sourced insertion at p_copy_source 0; at 1 nearly every insertion is copy-sourced (the first of a
    family in a lineage has no copy to draw from)."""
    _, truth = run_cli(tmp_path, {**scenario(), "p_copy_source": share})
    copies = rows(truth / "copies.tsv")
    frac = sum(r["source"] != "." for r in copies) / len(copies)
    assert (frac == 0) if expect == 0 else frac > 0.5


def test_table_ltr_offsets_run_along_the_record(tiny):
    """Both LTRs of a copy have its family's length, the left one first, on either strand."""
    with open(tiny[0] / "k2p.tsv") as fh:
        head = fh.readline().lstrip("#").rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(head, line.rstrip("\n").split("\t")))
            n, lt = int(r["seq_len"]), int(r["ltr5_len"])
            assert (int(r["ltr5_start"]), int(r["ltr5_end"])) == (1, lt)
            assert (int(r["ltr3_start"]), int(r["ltr3_end"])) == (n - lt + 1, n)
            assert r["ltr3_len"] == r["ltr5_len"]


@pytest.mark.parametrize("text, names", [
    ("{", "not valid JSON"),
    ("[1]", "the scenario must be a JSON object"),
    (json.dumps({"mergers": [{"time": 1e5, "lineages": ["A", "Z"]}]}), "not leaves of the tree"),
    (json.dumps({"chromosomes": {"n": -3}}), "chromosomes.n must be an integer >= 1: -3"),
    (json.dumps({"chromosomes": {"length": 0}}), "chromosomes.length must be an integer >= 1: 0"),
    (json.dumps({"families": {"n": 0}}), "families.n must be an integer >= 1: 0"),
    (json.dumps({"families": {"masters": [3, 1]}}), "families.masters must be [low, high]"),
    (json.dumps({"half_life": -1}), "half_life must be > 0: -1"),
    (json.dumps({"activity": {"p_burst": 1.5}}), "activity.p_burst must be in [0, 1]: 1.5"),
    (json.dumps({"tree": {"split": 3e7, "children": ["A", "B"]}}), "split time"),
    (json.dumps({"tree": {"split": 2e6}}),
     "a split node of tree lacks the key 'children' (the lineages it splits into"),
    (json.dumps({"mergers": [{"lineages": ["A", "B"]}]}),
     "a merger lacks the key 'time' (when the nuclei merge, years ago)"),
    (json.dumps({"exchange_list": [{"type": "reciprocal", "time": 5e4, "recipient": "B", "donor": "A", "start": 0,
                                    "end": 1e6}]}),
     "an exchange of exchange_list lacks the key 'chrom' (its chromosome"),
], ids=["json", "not-object", "unknown-leaf", "chromosomes-negative", "length-zero", "families-zero", "masters",
        "half-life", "share", "split-before-root", "tree-key", "merger-key", "exchange-key"])
def test_bad_scenario_is_a_usage_error(tmp_path, capsys, text, names):
    """A scenario file that is not a JSON object, or holds a value the simulation cannot use, is a usage error:
    exit 2 and one `error:` line naming the file and what is wrong, before anything is simulated or written."""
    f = tmp_path / "s.json"
    f.write_text(text)
    assert sim_main(["nucleus", str(f), "-o", str(tmp_path / "l"), "--truthdir", str(tmp_path / "t")]) == 2
    err = capsys.readouterr().err.splitlines()
    assert len(err) == 1 and err[0].startswith(f"error: {f}: ") and names in err[0], err
    assert not (tmp_path / "l").exists() and not (tmp_path / "t").exists()


def test_seed_changes_the_run(tmp_path):
    lib, _ = run_cli(tmp_path, scenario(), seed=8)
    assert (lib / "ltr.fa").read_bytes() != (GOLDEN / "lib" / "ltr.fa").read_bytes()


def test_help_prints_defaults_without_the_simulator():
    code = ("import sys\nfrom ltrtrace import cli\n"
            "try:\n    cli.main(['sim', 'nucleus', '--help'])\nexcept SystemExit:\n    pass\n"
            "print([m for m in ('numpy', 'ltrtrace.sim.nucleus.model', 'ltrtrace.sim.nucleus.write') "
            "if m in sys.modules])")
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    assert r.stdout.strip().splitlines()[-1] == "[]"
    text = " ".join(r.stdout.split())
    assert "half_life 1500000.0" in text and "singleton_share 0.14" in text
