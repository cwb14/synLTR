"""stats.tsv rows and findings.txt: the row format, the testability guard of the block jackknife (`why`, `guard`),
the resolution-symmetry warning, the findings text helpers and the one findings.txt writer (`write_findings`). Each
route keeps its findings.txt layout as the writer's named parameters: gene flow passes its header line and SECTIONS
(section blocks and a Warnings block); subgenomes passes its fixed header line and no sections (one statement per
line)."""
from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..base.checkpoint import atomic_write
from ..base.constants import JK_MIN_BLOCKS, RESOLVE_WARN
from ..calls.records import tsv_field
from .resample import Jk

STATS_COLUMNS = ("layer", "test", "scope", "stratum", "statistic", "value", "se", "z", "n", "loci", "note")


@dataclass
class StatRow:
    """One line of stats.tsv; "." is a missing text field, nan a missing number."""
    layer: str
    test: str
    scope: str
    stratum: str
    statistic: str
    value: float | int | str
    se: float = math.nan
    z: float = math.nan
    n: int | None = None
    loci: str = "."
    note: str = "."


def resolution(test: str, scope: str, fracs: dict[str, float], n: int, loci: str, findings: Findings,
               section: str) -> list[StatRow]:
    """One row resolved_<genome> per genome (layer pav, stratum all). Warns in `findings` when the largest and
    smallest resolved fractions differ by more than RESOLVE_WARN (a heuristic): a test
    that compares genomes of unequal resolution is biased toward the better-resolved one."""
    rows = [StatRow("pav", test, scope, "all", f"resolved_{g}", f, n=n, loci=loci) for g, f in fracs.items()]
    if fracs and max(fracs.values()) - min(fracs.values()) > RESOLVE_WARN:
        spread = ", ".join(f"{g} {f:.2f}" for g, f in fracs.items())
        findings.warn(section, f"{test} {scope}: resolved fractions differ by more than {RESOLVE_WARN} "
                               f"({spread}); the comparison favors the better-resolved genome "
                               f"(heuristic threshold)")
    return rows


def format_row(r: StatRow) -> str:
    """One stats.tsv line (`tsv_field` values), with its newline."""
    return "\t".join(tsv_field(getattr(r, c)) for c in STATS_COLUMNS) + "\n"


def write_stats(rows: list[StatRow], path: Path) -> None:
    """stats.tsv: a header row of STATS_COLUMNS only, then one line per row (`format_row`)."""
    with open(path, "w", newline="") as fh:
        fh.write("\t".join(STATS_COLUMNS) + "\n")
        fh.writelines(format_row(r) for r in rows)


def num(v) -> float:
    """A number that may be missing: nan for None or a stats.tsv "." field, else float(v)."""
    return math.nan if v is None or v == "." else float(v)


def read_stats(path: Path) -> list[StatRow]:
    """Inverse of `write_stats` (values to 6 significant digits; a non-numeric value stays text). A value is read as
    an int, else a float ("." as nan), else text, so a text value that looks like a number comes back as one ("12"
    -> 12) and an integral float as an int when written so; no row writes such a text value."""
    rows = []
    with open(path, newline="") as fh:
        rd = csv.reader(fh, delimiter="\t")
        if tuple(next(rd)) != STATS_COLUMNS:
            raise ValueError(f"{path}: not a stats.tsv (header is not {STATS_COLUMNS})")
        for f in rd:
            if len(f) != len(STATS_COLUMNS):
                raise ValueError(f"{path}: expected {len(STATS_COLUMNS)} fields, got {len(f)}: {f}")
            value: float | int | str
            try:
                value = int(f[5])
            except ValueError:
                try:
                    value = num(f[5])
                except ValueError:
                    value = f[5]
            rows.append(StatRow(*f[:5], value, num(f[6]), num(f[7]), None if f[8] == "." else int(f[8]),
                                f[9], f[10]))
    return rows


SECTIONS = ("Run", "Relationships", "Trio tests", "Painting", "Dating", "Inversions", "Mechanisms",
            "Sequence-identity layer (outside the PAV-only core)")


class Findings:
    """Statements for findings.txt, each under a section of its route's layout (gene flow: SECTIONS), kept in
    insertion order."""

    def __init__(self) -> None:
        self._items: list[tuple[str, str]] = []

    def add(self, section: str, text: str) -> None:
        self._items.append((section, text))

    def warn(self, section: str, text: str) -> None:
        self._items.append((section, "WARNING: " + text))

    def extend(self, other: Findings) -> None:
        self._items.extend(other._items)

    def items(self) -> list[tuple[str, str]]:
        return list(self._items)


def write_findings(path: Path, findings: Findings, *, header: str, sections: tuple[str, ...] | None = None
                   ) -> list[str]:
    """findings.txt at `path` (UTF-8, written atomically): `header`, a blank line, then the statements, each on one
    line (a newline inside one written as a space). Returns the lines written.

    With `sections` (gene flow's layout): each section that has statements, in `sections` order, under its name, one
    "- " line per statement and a blank line after; then "Warnings", repeating every WARNING statement after its
    section's name ("- none" without one). A statement under a section outside `sections` is a ValueError, raised
    before anything is written. Without (the subgenomes layout): the statements in the order added, no section
    names."""
    found = []
    for section, text in findings.items():
        if sections is not None and section not in sections:
            raise ValueError(f"{section!r} is not a findings section ({', '.join(sections)})")
        found.append((section, " ".join(text.splitlines())))
    lines = [header, ""]
    if sections is None:
        lines += [text for _, text in found]
    else:
        for section in sections:
            mine = [text for sec, text in found if sec == section]
            if mine:
                lines += [section, *(f"- {text}" for text in mine), ""]
        warnings = [f"- {sec}: {text}" for sec, text in found if text.startswith("WARNING: ")]
        lines += ["Warnings", *(warnings or ["- none"])]
    atomic_write(Path(path), "\n".join(lines) + "\n")
    return lines


def fmt(v: float, spec: str) -> str:
    """A number for findings.txt; nan as "nan"."""
    return "nan" if math.isnan(v) else format(v, spec)


def guard(jk: Jk, why: str | None) -> tuple[Jk, str | None]:
    """jk with se and z nan, and the reason, when `why` says the jackknife cannot estimate the variance, or
    the estimate is undefined, or the pseudovalues do not vary; the estimate stands. (jk, None) otherwise."""
    if why is None and not math.isfinite(jk.est):
        why = "the estimate is undefined"
    if why is None and not jk.se > 0:
        why = "the block pseudovalues do not vary" if jk.se == 0 else "the jackknife gives no se"
    return (jk, None) if why is None else (Jk(jk.est, math.nan, math.nan, jk.nblocks), why)


def why(classes: dict[str, int], blocks: np.ndarray) -> str | None:
    """Why a delete-one block jackknife cannot estimate the variance of a statistic comparing `classes`
    (name -> informative loci), whose informative loci lie in `blocks` (their block ids): an empty class,
    whose sampling variance never shows among the blocks, or fewer than JK_MIN_BLOCKS distinct blocks. None
    when neither."""
    empty = [name for name, n in classes.items() if not n]
    if empty:
        return "no " + " and no ".join(empty) + " locus"
    nb = len(np.unique(blocks))
    if nb < JK_MIN_BLOCKS:
        return f"its informative loci lie in {nb} block{'' if nb == 1 else 's'}, fewer than {JK_MIN_BLOCKS}"
    return None
