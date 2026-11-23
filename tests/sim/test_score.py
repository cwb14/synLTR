"""synltr.sim.score and `synltr sim score` on a tiny hand-made table (genome A scored; B rows are noise)."""
from __future__ import annotations

import csv
import gzip
from pathlib import Path

import pytest
from sim.poptoy import sim_main

from synltr.calls.records import CALL_COLUMNS
from synltr.io.elements import Element, make_eid, write_elements
from synltr.sim.score import SCORE_COLUMNS, TRUTH_COLUMNS, confusion, score

# name: (truth state, truth flags, call state, tier, call flags); call None = no call row
TABLE = {
    "n1": ("P", "inv", "P", "high", "inv"),
    "n2": ("P", "", "P", "high", "inv"),
    "n3": ("E", "", "P", "high", ""),
    "n4": ("E", "host_empty", "E", "high", ""),
    "n5": ("D", "", "D", "medium", ""),
    "n6": ("S", "", "U", ".", ""),
    "n7": ("T", "inv", "T", "medium", "inv"),
    "n8": ("P", "", None, None, None),
}
DROPPED = "n9"


def eid(name: str) -> str:
    return make_eid("A", "chr", int(name[1:]) * 100, int(name[1:]) * 100 + 50)


def write_tsv(path: Path, columns, rows) -> None:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(columns)
        w.writerows([r.get(c, "") for c in columns] for r in rows)


@pytest.fixture
def inputs(tmp_path):
    els = [Element("A", eid(n), n, "chr", int(n[1:]) * 100, int(n[1:]) * 100 + 50, "+", (1, 2), (3, 4), None,
                   None, None, None, None, None) for n in TABLE]
    write_elements(els, tmp_path / "A.tsv")
    truth = [{"genome": "A", "eid": n, "name": n, "state": t[0], "flags": t[1]} for n, t in TABLE.items()]
    truth += [{"genome": "A", "eid": DROPPED, "name": DROPPED, "state": "E"},
              {"genome": "B", "eid": "n1", "name": "n1", "state": "D"}]
    write_tsv(tmp_path / "truth.tsv", TRUTH_COLUMNS, truth)
    calls = [{"eid": eid(n), "state": t[2], "tier": t[3], "flags": t[4]} for n, t in TABLE.items() if t[2]]
    write_tsv(tmp_path / "calls.tsv.gz", CALL_COLUMNS, calls)
    return tmp_path / "calls.tsv.gz", tmp_path / "A.tsv", tmp_path / "truth.tsv"


def pick(rows, state, tier_set):
    (row,) = [r for r in rows if r["state"] == state and r["tier_set"] == tier_set]
    return row


def test_state_rows_all(inputs):
    rows = score(*inputs)
    assert all(tuple(r) == SCORE_COLUMNS for r in rows)
    p = pick(rows, "P", "all")
    assert (p["tp"], p["fp"], p["fn"]) == (2, 1, 1)
    assert p["precision"] == pytest.approx(2 / 3) and p["recall"] == pytest.approx(2 / 3)
    e = pick(rows, "E", "all")
    assert (e["tp"], e["fp"], e["fn"]) == (1, 0, 1) and e["recall"] == 0.5
    t = pick(rows, "T", "all")
    assert (t["tp"], t["fp"], t["fn"]) == (1, 0, 0)
    o = pick(rows, "O", "all")  # no truth, no calls
    assert (o["tp"], o["precision"], o["recall"]) == (0, None, None)


def test_high_tier_restriction(inputs):
    rows = score(*inputs)
    t = pick(rows, "T", "high")  # the T call is medium: a miss, not a call
    assert (t["tp"], t["fp"], t["fn"], t["recall"]) == (0, 0, 1, 0.0) and t["precision"] is None
    d = pick(rows, "D", "high")
    assert (d["tp"], d["fn"]) == (0, 1)
    p = pick(rows, "P", "high")
    assert (p["tp"], p["fp"], p["fn"]) == (2, 1, 1)


def test_class_rows(inputs):
    rows = score(*inputs)
    present, absent = pick(rows, "present", "all"), pick(rows, "absent", "all")
    assert (present["tp"], present["fp"], present["fn"]) == (3, 1, 2)
    assert (absent["tp"], absent["fp"], absent["fn"]) == (1, 0, 1)
    hp = pick(rows, "present", "high")
    assert (hp["tp"], hp["fp"], hp["fn"]) == (2, 1, 3)


def test_flag_rows(inputs):
    rows = score(*inputs)
    inv = pick(rows, "flag:inv", "all")  # over n1, n2, n4, n5, n7 (call == truth)
    assert (inv["tp"], inv["fp"], inv["fn"]) == (2, 1, 0)
    he = pick(rows, "flag:host_empty", "all")
    assert (he["tp"], he["fp"], he["fn"]) == (0, 0, 1)
    hinv = pick(rows, "flag:inv", "high")  # n1, n2, n4 only
    assert (hinv["tp"], hinv["fp"], hinv["fn"]) == (1, 1, 0)


def test_summary(inputs):
    rows = score(*inputs)
    a, h = pick(rows, "summary", "all"), pick(rows, "summary", "high")
    assert a["n_scored"] == 8 and a["n_unannotated"] == 1 and a["n_missing_call"] == 1
    assert a["correct_resolved_frac"] == pytest.approx(4 / 7)  # n1 n2 n4 n7 of 7 resolved truths
    assert h["correct_resolved_frac"] == pytest.approx(3 / 7)  # n7 is medium
    assert a["unresolved_frac"] == 0.25 and h["unresolved_frac"] == 0.25  # n6 (U) and n8 (no call)
    assert a["tp"] is None and a["precision"] is None


def test_confusion(inputs):
    got = {(r["tier_set"], r["truth"], r["call"]): r["n"] for r in confusion(*inputs)}
    assert got[("all", "P", "P")] == 2 and got[("all", "E", "P")] == 1 and got[("all", "P", "U")] == 1
    assert got[("all", "S", "U")] == 1 and ("all", "D", "not_high") not in got
    assert got[("high", "D", "not_high")] == 1 and got[("high", "T", "not_high")] == 1
    assert got[("high", "P", "not_high")] == 1 and got[("high", "P", "P")] == 2
    assert sum(v for k, v in got.items() if k[0] == "high") == 8
    order = [(r["tier_set"], r["truth"], r["call"]) for r in confusion(*inputs)]
    assert order[0] == ("all", "P", "P") and order[1] == ("all", "P", "U")


def test_plain_calls(inputs, tmp_path):
    calls, els, truth = inputs
    plain = tmp_path / "plain.tsv"
    plain.write_bytes(gzip.decompress(calls.read_bytes()))
    assert score(plain, els, truth) == score(calls, els, truth)


def test_call_eid_not_in_elements(inputs, tmp_path):
    _, els, truth = inputs
    write_tsv(tmp_path / "bad.tsv", CALL_COLUMNS, [{"eid": "A|chr:1-2", "state": "P", "tier": "high"}])
    with pytest.raises(ValueError, match="A.chr:1-2"):
        score(tmp_path / "bad.tsv", els, truth)


def test_element_without_truth(inputs, tmp_path):
    calls, els, _ = inputs
    write_tsv(tmp_path / "t2.tsv", TRUTH_COLUMNS, [{"genome": "A", "name": "n1", "state": "P"}])
    with pytest.raises(ValueError, match=r"7 annotated.*n2"):
        score(calls, els, tmp_path / "t2.tsv")


def test_missing_column(inputs, tmp_path):
    calls, els, _ = inputs
    write_tsv(tmp_path / "t3.tsv", ("genome", "name"), [])
    with pytest.raises(ValueError, match="t3.tsv.*state"):
        score(calls, els, tmp_path / "t3.tsv")


def test_cli(inputs, tmp_path, capsys):
    calls, els, truth = inputs
    out = tmp_path / "out" / "ab"
    assert sim_main(["score", "--calls", str(calls), "--elements", str(els), "--truth", str(truth),
                     "--out", str(out)]) == 0
    with open(f"{out}.score.tsv") as fh:
        lines = [ln.rstrip("\n").split("\t") for ln in fh]
    assert tuple(lines[0]) == SCORE_COLUMNS
    p = dict(zip(lines[0], next(ln for ln in lines if ln[:2] == ["P", "all"])))
    assert p["precision"] == "0.6667" and p["tp"] == "2"
    o = dict(zip(lines[0], next(ln for ln in lines if ln[:2] == ["O", "all"])))
    assert o["precision"] == "NA"
    with open(f"{out}.confusion.tsv") as fh:
        assert fh.readline().split() == ["tier_set", "truth", "call", "n"]
    assert sim_main(["score", "--calls", str(tmp_path / "nope"), "--elements", str(els),
                     "--truth", str(truth), "--out", str(out)]) == 1
    assert capsys.readouterr().err.count("error:") == 1


def _other_genome_truth(d, calls, els, truth):
    write_tsv(d / "other.tsv", TRUTH_COLUMNS, [{"genome": "A", "name": "n1", "state": "P"}])
    return calls, els, d / "other.tsv", d / "other.tsv"


def _calls_of_other_elements(d, calls, els, truth):
    write_tsv(d / "c.tsv", CALL_COLUMNS, [{"eid": "A|chr:1-2", "state": "P", "tier": "high"}])
    return d / "c.tsv", els, truth, d / "c.tsv"


def _elements_of_two_genomes(d, calls, els, truth):
    write_elements([Element(g, f"{g}|c:1-9", g.lower(), "c", 1, 9, "+", (1, 2), (8, 9), None, None, None, None, None,
                            None) for g in "AB"], d / "two.tsv")
    return calls, d / "two.tsv", truth, d / "two.tsv"


def _not_a_call_file(d, calls, els, truth):
    write_tsv(d / "x.tsv", ("eid", "state"), [{"eid": eid("n1"), "state": "P"}])
    return d / "x.tsv", els, truth, d / "x.tsv"


@pytest.mark.parametrize("make", [_other_genome_truth, _calls_of_other_elements, _elements_of_two_genomes,
                                  _not_a_call_file], ids=lambda f: f.__name__.strip("_"))
def test_cli_files_that_do_not_match_are_usage_errors(inputs, tmp_path, capsys, make):
    """Calls, elements and truth that do not belong together, or a calls file that is not a pair call file: a usage
    error, exit 2, one `error:` line naming the file, nothing written."""
    calls, els, truth, named = make(tmp_path, *inputs)
    out = tmp_path / "out" / "ab"
    assert sim_main(["score", "--calls", str(calls), "--elements", str(els), "--truth", str(truth),
                     "--out", str(out)]) == 2
    err = [ln for ln in capsys.readouterr().err.splitlines() if "error:" in ln]
    assert len(err) == 1 and err[0].startswith("error: ") and str(named) in err[0], err
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize("which, rows, message", [
    ("truth", [{"genome": "A", "name": "n1", "state": "P"}, {"genome": "A", "name": "n1", "state": "E"}],
     "name 'n1' appears more than once"),
    ("truth", [{"genome": "A", "name": "n1", "state": "Q"}], "name 'n1' has state 'Q'"),
    ("calls", [{"eid": "A|chr:100-150", "state": "P", "tier": "high"}] * 2, "appears more than once"),
    ("calls", [{"eid": "A|chr:100-150", "state": "", "tier": "high"}], "has state ''"),
])
def test_duplicate_or_unknown_rows_fail(inputs, tmp_path, which, rows, message):
    """A duplicate truth or call row, or a state outside the table, fails instead of being
    overwritten or ignored."""
    calls, els, truth = inputs
    if which == "truth":
        write_tsv(tmp_path / "t.tsv", TRUTH_COLUMNS, rows)
        truth = tmp_path / "t.tsv"
    else:
        write_tsv(tmp_path / "c.tsv", CALL_COLUMNS, rows)
        calls = tmp_path / "c.tsv"
    with pytest.raises(ValueError, match=message):
        score(calls, els, truth)


def test_gzip_told_by_content(inputs, tmp_path):
    """A gzip calls file without .gz in its name is read as gzip."""
    calls, els, truth = inputs
    renamed = tmp_path / "calls.tsv"
    renamed.write_bytes(calls.read_bytes())
    assert score(renamed, els, truth) == score(calls, els, truth)


def test_elements_of_one_genome(inputs, tmp_path):
    """The elements file must hold exactly one genome's elements."""
    calls, _, truth = inputs
    two = tmp_path / "two.tsv"
    write_elements([Element("A", "A|c:1-9", "a", "c", 1, 9, "+", (1, 2), (8, 9), None, None, None, None, None,
                            None), Element("B", "B|c:1-9", "b", "c", 1, 9, "+", (1, 2), (8, 9), None, None, None,
                                           None, None, None)], two)
    with pytest.raises(ValueError, match=r"exactly one genome, found \['A', 'B'\]"):
        score(calls, two, truth)
    empty = tmp_path / "empty.tsv"
    write_elements([], empty)
    with pytest.raises(ValueError, match="exactly one genome, found none"):
        score(calls, empty, truth)
