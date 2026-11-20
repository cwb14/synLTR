import gzip
import json
import math

import numpy as np
import pytest

from ltrtrace.base.checkpoint import is_done, marker
from ltrtrace.calls.records import (
    CALL_COLUMNS,
    DIV_COLUMNS,
    call_line,
    read_calls,
    read_calls_meta,
    read_divergences,
    stream_calls,
    tsv_field,
    write_call_lines,
    write_calls,
    write_divergences,
)


def header(path):
    with gzip.open(path, "rt") as fh:
        return fh.readline().rstrip("\n").split("\t")


def row(eid="A|c1:1-9", **kw):
    base = {"eid": eid, "a_seqid": "c1", "a_start": 1, "a_end": 9, "parent": None, "level": 0, "state": "P",
            "tier": "high", "flags": "nestB", "reason": None, "b_seqid": "c1", "b_start": 11, "b_end": 19,
            "b_strand": "+", "wid": 3, "cov_left": 1.0, "cov_right": 0.5, "cov_int": 0.25, "cov_elem": 0.0,
            "dl": None, "dr": 2, "b_ins": 0, "tsd_a": "ACGTA", "tsd_b": None, "tsd_mm": 1, "ident": 0.978916,
            "anchor_left": 0, "anchor_right": 5, "margin": 1.0}
    return {**base, **kw}


def test_calls_roundtrip(tmp_path):
    rows = [row(), row("A|c1:20-30", state="E", flags=None, b_strand=None, ident=None, parent="A|c1:1-9", level=1)]
    path = tmp_path / "pairs" / "A__B.calls.tsv.gz"
    meta = {"engine": "library", "n": 2}
    write_calls(path, rows, meta)
    assert list(read_calls(path)) == rows
    assert read_calls_meta(path) == meta
    assert is_done(marker(path), json.loads(marker(path).read_text())["key"])
    assert header(path) == list(CALL_COLUMNS)


def test_calls_file_is_deterministic_and_key_is_honored(tmp_path):
    write_calls(tmp_path / "x.gz", [row()], {"engine": "assembly", "key": "k1"})
    write_calls(tmp_path / "y.gz", [row()], {"engine": "assembly"})
    assert (tmp_path / "x.gz").read_bytes() == (tmp_path / "y.gz").read_bytes()
    assert json.loads(marker(tmp_path / "x.gz").read_text())["key"] == "k1"
    assert read_calls_meta(tmp_path / "x.gz") == {"engine": "assembly"}  # the key is not meta


def test_call_lines_make_the_same_file(tmp_path):
    """Rows formatted by call_line, or lines copied from another call file, write the same bytes and marker."""
    rows = [row(), row("A|c1:20-30", state="E", flags=None, ident=None)]
    write_calls(tmp_path / "rows.gz", rows, {"engine": "assembly", "key": "k"})
    write_call_lines(tmp_path / "lines.gz", (call_line(r[c] for c in CALL_COLUMNS) for r in rows),
                     {"engine": "assembly", "key": "k"})
    with gzip.open(tmp_path / "rows.gz", "rt", newline="") as fh:
        fh.readline()
        write_call_lines(tmp_path / "copied.gz", fh, {"engine": "assembly", "key": "k"})
    assert (tmp_path / "rows.gz").read_bytes() == (tmp_path / "lines.gz").read_bytes()
    assert (tmp_path / "rows.gz").read_bytes() == (tmp_path / "copied.gz").read_bytes()
    assert read_calls_meta(tmp_path / "lines.gz") == {"engine": "assembly"}


def test_engine_is_required(tmp_path):
    for meta in ({}, {"engine": "reads"}):
        with pytest.raises(ValueError, match="engine"):
            write_calls(tmp_path / "a.gz", [row()], meta)
        assert not (tmp_path / "a.gz").exists()


def test_meta_needs_a_finished_file(tmp_path):
    with pytest.raises(ValueError, match="no readable done marker"):
        read_calls_meta(tmp_path / "none.gz")


def test_stream_calls_returns_the_asked_fields_and_names_bad_files(tmp_path):
    path = tmp_path / "c.gz"
    write_calls(path, [row(), row("A|c1:20-30")], {"engine": "assembly"})
    assert list(stream_calls(path, ("eid", "state"))) == [(2, ("A|c1:1-9", "P")), (3, ("A|c1:20-30", "P"))]
    path.write_bytes(path.read_bytes()[:-10])
    with pytest.raises(ValueError, match="corrupt pair file"):
        list(stream_calls(path))
    with gzip.open(path, "wt") as fh:
        fh.write("\t".join(CALL_COLUMNS) + "\nx\ty\n")
    with pytest.raises(ValueError, match=r"line 2 has 2 column\(s\); a pair file has 29"):
        list(read_calls(path))


def test_divergences_roundtrip(tmp_path):
    rows = [{"eid": "A|c1:1-9", "b_eid": "B|c1:11-19", "n_sites": 300, "n_diff": 6, "d": 0.02},
            {"eid": "A|c1:20-30", "b_eid": "B|c2:1-11", "n_sites": 280, "n_diff": 0, "d": 0.0}]
    path = tmp_path / "pairs" / "A__B.div.tsv.gz"
    meta = {"engine": "library", "d_kind": "K2P"}
    write_divergences(path, rows, meta)
    assert list(read_divergences(path)) == rows
    assert read_calls_meta(path) == meta
    assert header(path) == list(DIV_COLUMNS)
    with pytest.raises(ValueError, match="d_kind"):
        write_divergences(tmp_path / "bad.gz", rows, {"engine": "library"})


@pytest.mark.parametrize("value, text", [(None, "."), (float("nan"), "."), (True, "1"), (0.1234567, "0.123457"),
                                         (7, "7"), ("x", "x")])
def test_tsv_field(value, text):
    assert tsv_field(value) == text


def test_tsv_field_values():
    vals = [None, math.nan, np.float32("nan"), True, False, np.bool_(True), 3, np.int64(-4), 0.1234567, 1.0,
            2.5e-7, np.float32(0.7), "E"]
    assert [tsv_field(v) for v in vals] == [".", ".", ".", "1", "0", "1", "3", "-4", "0.123457", "1",
                                            "2.5e-07", "0.7", "E"]
