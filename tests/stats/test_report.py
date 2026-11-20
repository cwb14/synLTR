import math

import numpy as np
import pytest

from ltrtrace.stats.report import (
    SECTIONS,
    STATS_COLUMNS,
    Findings,
    StatRow,
    fmt,
    guard,
    read_stats,
    resolution,
    why,
    write_findings,
    write_stats,
)
from ltrtrace.stats.resample import Jk


def test_resolution_rows_and_warning():
    f = Findings()
    rows = resolution("resolution", "R,S;D", {"R": 0.85, "S": 0.40}, 1000, "disc=D", f, "Trio tests")
    assert [(r.statistic, r.value, r.n, r.loci) for r in rows] == [("resolved_R", 0.85, 1000, "disc=D"),
                                                                ("resolved_S", 0.40, 1000, "disc=D")]
    assert len(f.items()) == 1 and f.items()[0][1].startswith("WARNING:") and "heuristic" in f.items()[0][1]
    g = Findings()
    resolution("resolution", "R,S;D", {"R": 0.85, "S": 0.80}, 1000, "disc=D", g, "Trio tests")
    assert g.items() == []


def test_stats_roundtrip(tmp_path):
    rows = [StatRow("pav", "D_pav", "R,S;D", "all", "D", 0.123456789, 0.01, 12.3456789, 400, "disc=D", "."),
            StatRow("pav", "tree", "all", "all", "newick", "(A,B);")]
    write_stats(rows, tmp_path / "s.tsv")
    lines = (tmp_path / "s.tsv").read_text().splitlines()
    assert lines[0].split("\t") == list(STATS_COLUMNS) and len(lines) == 3
    assert lines[1].split("\t")[5:9] == ["0.123457", "0.01", "12.3457", "400"]
    assert lines[2].split("\t")[5:9] == ["(A,B);", ".", ".", "."]
    back = read_stats(tmp_path / "s.tsv")
    assert back[0].value == pytest.approx(0.123457) and math.isnan(back[1].se) and back[1].value == "(A,B);"


def test_findings_order():
    f = Findings()
    f.add("Painting", "p1"); f.warn("Run", "w1"); f.add("Run", "r1")
    assert f.items() == [("Painting", "p1"), ("Run", "WARNING: w1"), ("Run", "r1")]


def test_text_that_looks_like_a_number(tmp_path):
    """read_stats reads a value as an int, else a float, else text: a text value that looks like a number comes
    back as one (no row writes such a text value)."""
    rows = [StatRow("pav", "x", "all", "all", "key", "1234e5"), StatRow("pav", "x", "all", "all", "n", "12"),
            StatRow("pav", "x", "all", "all", "name", "R+S")]
    write_stats(rows, tmp_path / "s.tsv")
    assert [r.value for r in read_stats(tmp_path / "s.tsv")] == [1.234e8, 12, "R+S"]


def test_why_and_guard():
    blocks = np.array([0, 0, 1, 1, 2])
    assert why({"R": 3, "S": 0}, blocks) == "no S locus"
    assert why({"R": 0, "S": 0}, blocks) == "no R and no S locus"
    assert why({"R": 3, "S": 2}, np.array([0, 0, 1])) == "its informative loci lie in 2 blocks, fewer than 3"
    assert why({"R": 3}, np.array([4])) == "its informative loci lie in 1 block, fewer than 3"
    assert why({"R": 3, "S": 2}, blocks) is None
    ok = Jk(0.5, 0.1, 5.0, 4)
    assert guard(ok, None) == (ok, None)
    assert guard(ok, "r") == (Jk(0.5, math.nan, math.nan, 4), "r")
    assert guard(Jk(math.nan, math.nan, math.nan, 4), None)[1] == "the estimate is undefined"
    assert guard(Jk(0.5, 0.0, math.nan, 3), None)[1] == "the block pseudovalues do not vary"
    assert guard(Jk(0.5, math.nan, math.nan, 2), None)[1] == "the jackknife gives no se"


def test_fmt():
    assert fmt(math.nan, ".3f") == "nan" and fmt(0.12345, ".3f") == "0.123"


def test_write_findings_sections_layout(tmp_path):
    """Gene flow's layout: the header, a blank line, each section with findings in `sections` order under its name
    (one "- " line each, a newline inside one written as a space, a blank line after), then Warnings repeating every
    WARNING after its section's name; returns the lines written; a section outside `sections` writes nothing."""
    f = Findings()
    f.add("Painting", "p1\nsecond line")
    f.warn("Run", "w1")
    f.add("Run", "r1")
    path = tmp_path / "findings.txt"
    lines = write_findings(path, f, header="ltrtrace x findings for the run in /r", sections=SECTIONS)
    want = ["ltrtrace x findings for the run in /r", "", "Run", "- WARNING: w1", "- r1", "", "Painting",
            "- p1 second line", "", "Warnings", "- Run: WARNING: w1"]
    assert lines == want and path.read_bytes() == ("\n".join(want) + "\n").encode()
    g = Findings()
    g.add("Run", "r1")
    assert write_findings(path, g, header="h", sections=SECTIONS)[-2:] == ["Warnings", "- none"]
    h = Findings()
    h.add("Nowhere", "x")
    with pytest.raises(ValueError, match="not a findings section"):
        write_findings(tmp_path / "other.txt", h, header="h", sections=SECTIONS)
    assert sorted(p.name for p in tmp_path.iterdir()) == ["findings.txt"]       # nothing else, no temporary file


def test_write_findings_lines_layout(tmp_path):
    """The subgenomes layout (no sections): the header, a blank line, then one statement per line in the order
    added, under any section name, WARNINGs in place; UTF-8."""
    f = Findings()
    f.add("Homoeologs", "Homoeologs SG1-SG2: 8 of 8 SG1 and 8 of 8 SG2 sequences paired.")
    f.add("Other", "Orientation: every called pair runs the same way \u2265 1 Mb.")
    f.warn("Homoeologs", "w")
    path = tmp_path / "findings.txt"
    lines = write_findings(path, f, header="What the pairs say.")
    assert lines == ["What the pairs say.", "", "Homoeologs SG1-SG2: 8 of 8 SG1 and 8 of 8 SG2 sequences paired.",
                     "Orientation: every called pair runs the same way \u2265 1 Mb.", "WARNING: w"]
    assert path.read_bytes() == ("\n".join(lines) + "\n").encode("utf-8")
