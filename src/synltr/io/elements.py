"""The element model: one record per LTR-retrotransposon copy (`Element`), a columnar table of them
(`ElementTable`) and the elements.tsv file.

Coordinates: `start`, `end` and the LTR extents are genome coordinates, 1-based closed (GFF3, VCF); `zstart` is
the 0-based start of the same interval, for code that slices sequences. An element without a known LTR pair (an
annotation read under ltr_rule "optional" or table_rule "phasing", or a library record without annotation) has both
LTR extents NO_LTR.
"""
from __future__ import annotations

import hashlib
import math
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, fields
from operator import attrgetter
from pathlib import Path

import numpy as np

NO_LTR = (0, 0)  # the LTR extent of an element whose annotation gives no LTR pair (no genome position is 0)


@dataclass(frozen=True, slots=True)
class Element:
    genome: str
    eid: str  # `genome|seqid:start-end`
    name: str  # the annotation's own ID
    seqid: str
    start: int  # 1-based closed
    end: int
    strand: str  # '+', '-' or '.'
    ltr_left: tuple[int, int]  # extent in genome orientation, 1-based closed, nests included; NO_LTR if unknown
    ltr_right: tuple[int, int]
    tsd: str | None  # genome orientation
    k: float | None  # divergence between the two LTRs, labelled by `k_kind`
    family: str | None
    superfamily: str | None
    clade: str | None
    confidence: str | None
    parent: str | None = None  # eid of the innermost containing element (assign_nesting); None = top-level
    level: int = 0
    ltr_left_len: int | None = None  # the LTRs' own lengths; None when the annotation does not give them
    ltr_right_len: int | None = None
    tsd_source: str | None = None  # "annot" or "inferred"; None when unknown (e.g. read back from elements.tsv)
    k_kind: str | None = "K2P"  # "K2P" or "JC69"; None when unknown (read back from elements.tsv)
    k_se: float | None = None
    n_sites: float | None = None
    n_ts: float | None = None
    n_tv: float | None = None

    @property
    def zstart(self) -> int:
        """0-based start of the element (`start` - 1), the first base of a half-open slice."""
        return self.start - 1


ALL_FIELDS: tuple[str, ...] = tuple(f.name for f in fields(Element))
# the columns of elements.tsv, in order: the fields up to `level`
_FIELDS_OF = attrgetter(*ALL_FIELDS)
ELEMENT_COLUMNS: tuple[str, ...] = ALL_FIELDS[:ALL_FIELDS.index("level") + 1]
_TUPLES = {"ltr_left", "ltr_right"}
_INTS = {"start", "end", "level"}
_OPT_INTS = {"ltr_left_len", "ltr_right_len"}  # int64 values in a float column, NaN for None
_FLOATS = {"k", "k_se", "n_sites", "n_ts", "n_tv"}


def make_eid(genome: str, seqid: str, start: int, end: int) -> str:
    return f"{genome}|{seqid}:{start}-{end}"


class ElementTable:
    """Elements as numpy columns, one per `Element` field: strings and None in object arrays, `start`, `end`,
    `level` as int64, the two LTR extents as (n, 2) int64, and the floats and the LTRs' own lengths as float64 with
    NaN for None."""

    __slots__ = ("_cols", "_n")

    def __init__(self, cols: dict[str, np.ndarray]):
        if set(cols) != set(ALL_FIELDS):
            raise ValueError(f"element columns {sorted(set(cols) ^ set(ALL_FIELDS))} missing or unknown")
        self._cols = cols
        self._n = len(cols["eid"])

    @classmethod
    def from_rows(cls, rows: Iterable[Element]) -> ElementTable:
        rows = list(rows)
        by_field = list(zip(*map(_FIELDS_OF, rows))) if rows else [()] * len(ALL_FIELDS)
        cols: dict[str, np.ndarray] = {}
        for name, vals in zip(ALL_FIELDS, by_field):
            if name in _TUPLES:
                cols[name] = np.array(vals, dtype=np.int64).reshape(len(rows), 2)
            elif name in _INTS:
                cols[name] = np.array(vals, dtype=np.int64)
            elif name in _FLOATS | _OPT_INTS:
                cols[name] = np.array([math.nan if v is None else v for v in vals], dtype=np.float64)
            else:
                col = np.empty(len(rows), dtype=object)
                col[:] = vals
                cols[name] = col
        return cls(cols)

    @classmethod
    def concat(cls, tables: Iterable[ElementTable]) -> ElementTable:
        """The rows of one or more tables, in order, as one table."""
        tables = list(tables)
        return cls({name: np.concatenate([t.column(name) for t in tables]) for name in ALL_FIELDS})

    def __len__(self) -> int:
        return self._n

    def __iter__(self) -> Iterator[Element]:
        return (self.row(i) for i in range(self._n))

    def column(self, name: str) -> np.ndarray:
        """The column of field `name`."""
        return self._cols[name]

    def row(self, i: int) -> Element:
        values = {}
        for name in ALL_FIELDS:
            v = self._cols[name][i]
            if name in _TUPLES:
                v = (int(v[0]), int(v[1]))
            elif name in _INTS:
                v = int(v)
            elif name in _FLOATS | _OPT_INTS:
                v = None if math.isnan(v) else (int(v) if name in _OPT_INTS else float(v))
            values[name] = v
        return Element(**values)

    def select(self, mask) -> ElementTable:
        """The rows picked by a boolean mask or an integer index array, as a new table."""
        return ElementTable({name: col[mask] for name, col in self._cols.items()})

    def replace(self, **columns: np.ndarray) -> ElementTable:
        """A table with the named columns replaced (same length)."""
        for name, col in columns.items():
            if name not in self._cols or len(col) != self._n:
                raise ValueError(f"cannot replace column {name!r}: unknown, or not {self._n} rows")
        return ElementTable({**self._cols, **columns})


def elements_digest(elements: Iterable[Element] | None) -> str:
    """sha256 of the element records (the fields of elements.tsv, in order); "none" without an annotation."""
    if elements is None:
        return "none"
    h = hashlib.sha256()
    for e in elements:
        h.update(f"Element({', '.join(f'{c}={getattr(e, c)!r}' for c in ELEMENT_COLUMNS)})".encode())
        h.update(b"\n")
    return h.hexdigest()


def _fmt(name: str, v) -> str:
    if v is None:
        return "."
    if name in _TUPLES:
        return f"{v[0]}-{v[1]}"
    if isinstance(v, float):
        return repr(v)
    return str(v)


def _parse(name: str, s: str):
    if name in _TUPLES:
        a, b = s.split("-")
        return (int(a), int(b))
    if name in _INTS:
        return int(s)
    if s == "." and name != "strand":
        return None
    if name == "k":
        return float(s)
    return s


def write_elements(elements: Iterable[Element], path: Path) -> None:
    """Write elements.tsv: the columns of ELEMENT_COLUMNS only (the other fields are not written)."""
    with open(path, "w") as fh:
        fh.write("\t".join(ELEMENT_COLUMNS) + "\n")
        fh.writelines("\t".join(_fmt(c, getattr(e, c)) for c in ELEMENT_COLUMNS) + "\n" for e in elements)


def read_elements(path: Path) -> ElementTable:
    """Read elements.tsv. The fields it does not hold are None; `k_kind` and `tsd_source` are unknown (None), not
    assumed."""
    out = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if tuple(header) != ELEMENT_COLUMNS:
            raise ValueError(f"{path}: unexpected header {header}")
        for line in fh:
            if not line.strip():
                continue
            vals = line.rstrip("\n").split("\t")
            out.append(Element(**{c: _parse(c, v) for c, v in zip(ELEMENT_COLUMNS, vals)}, k_kind=None))
    return ElementTable.from_rows(out)
