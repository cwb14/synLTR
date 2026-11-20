"""Pair call files and the pair-divergence layer: gzip TSVs, each with a done marker holding its meta.

A pair call file `<A>__<B>.calls.tsv.gz` holds every element of A once with its call in B (CALL_COLUMNS; "." for a
missing value; A and B coordinates 1-based closed). `<A>__<B>.div.tsv.gz` is the sparse pair-divergence layer
(DIV_COLUMNS): for each pair of orthologous copies, the sites compared and the differences among them. The meta of
both records `engine` ("assembly" or "library"); the divergence meta also records `d_kind`, what `d` is.

What the columns mean on each engine:
- assembly (routes.assembly): states P S T E D O G U; b_seqid, b_start, b_end the call's site in B; b_strand the
  window's orientation; wid the alignment window; b_ins the B bases inserted at the site; tier, flags and the evidence
  columns as routes.assembly.classify defines them; tsd_a, tsd_b in genome orientation.
- library (routes.library): states P (1:1 ortholog), A (inferred absence, never an observed empty site) and U (reason
  no_block); a P call's b_seqid, b_start and b_end are the ortholog's coordinates ("." when the ortholog is a record
  whose name holds no coordinates, as is then the b_eid of its divergence row), an A call's the window on B searched
  for the element (a search window between two anchors' midpoints, not a site); b_strand the window block's sign;
  wid the block's row among the pair's rows of blocks.tsv (0-based); b_ins the room, B's span between the midpoints
  of the two anchors flanking the element (not inserted bases); tsd_a, tsd_b in canonical orientation (the lesser of
  a TSD and its reverse complement); tier "." and the other evidence columns ".".
- divergence layer: d_kind "p" (assembly: the LASTZ p-distance over the aligned LTR sites) or "p_whole_unit"
  (library: mismatches over the aligned columns of the two whole N-stripped records, edlib unit cost); n_sites and
  n_diff the sites compared and differing; d = n_diff / n_sites is written to 6 significant digits, so a reader that
  needs the exact divergence takes n_diff / n_sites.

Files are gzip with mtime 0 and no stored name (the same rows give the same bytes), written to a temporary file
and renamed into place, then marked done (`checkpoint.marker`).

TSV values (`tsv_field`): None and nan -> "."; bool -> 1/0; float -> 6 significant digits; int, str as is.
"""
from __future__ import annotations

import contextlib
import gzip
import io
import math
import zlib
from collections.abc import Iterable, Iterator
from pathlib import Path

import numpy as np

from ltrtrace.base.checkpoint import atomic_path, clear_done, done_meta, mark_done, marker, sha256_file
from ltrtrace.base.constants import GZ_LEVEL
from ltrtrace.io.files import open_input

CALL_COLUMNS = ("eid", "a_seqid", "a_start", "a_end", "parent", "level", "state", "tier", "flags", "reason",
                "b_seqid", "b_start", "b_end", "b_strand", "wid", "cov_left", "cov_right", "cov_int",
                "cov_elem", "dl", "dr", "b_ins", "tsd_a", "tsd_b", "tsd_mm", "ident", "anchor_left",
                "anchor_right", "margin")
DIV_COLUMNS = ("eid", "b_eid", "n_sites", "n_diff", "d")
ENGINES = ("assembly", "library")
_INT_COLUMNS = frozenset(("a_start", "a_end", "level", "b_start", "b_end", "wid", "dl", "dr", "b_ins", "tsd_mm",
                          "anchor_left", "anchor_right", "n_sites", "n_diff"))
_FLOAT_COLUMNS = frozenset(("cov_left", "cov_right", "cov_int", "cov_elem", "ident", "margin", "d"))


def tsv_field(v) -> str:
    """A TSV value: None and nan -> "."; bool -> "1"/"0"; float -> 6 significant digits; int, str as is."""
    if v is None:
        return "."
    if isinstance(v, (bool, np.bool_)):
        return "1" if v else "0"
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if isinstance(v, (float, np.floating)):
        return "." if math.isnan(v) else format(float(v), ".6g")
    return str(v)


def write_calls(path: Path, rows: Iterable[dict], meta: dict) -> None:
    """Write the pair call file `path` (rows: dicts keyed by CALL_COLUMNS, in the order given) and its done
    marker. `meta` must hold `engine` ("assembly" or "library"); it is stored beside the file's key, `meta["key"]`
    when given, else the sha256 of the file."""
    write_call_lines(path, _lines(rows, CALL_COLUMNS), meta)


def write_call_lines(path: Path, lines: Iterable[str], meta: dict) -> None:
    """write_calls from rows already formatted: each line the tsv_field values of CALL_COLUMNS joined by tabs, with
    its newline (as call_line makes it, or as read from another call file)."""
    _write_table(Path(path), CALL_COLUMNS, lines, meta)


def call_line(values: Iterable) -> str:
    """One TSV line: the tsv_field values joined by tabs, with its newline."""
    return "\t".join(map(tsv_field, values)) + "\n"


def read_calls(path: Path) -> Iterator[dict]:
    """The rows of a pair call file as dicts: "." is None; a_start, a_end, level, b_start, b_end, wid and the
    count columns are int, the coverages, ident and margin float, the rest str."""
    return (dict(zip(CALL_COLUMNS, map(_parse, CALL_COLUMNS, f), strict=True)) for _, f in stream_calls(path))


def read_calls_meta(path: Path) -> dict:
    """The meta `write_calls` stored for `path`; a ValueError when the file has no readable done marker."""
    return _read_meta(Path(path))


def write_divergences(path: Path, rows: Iterable[dict], meta: dict) -> None:
    """Write the pair-divergence file `path` (rows keyed by DIV_COLUMNS: eid, b_eid, the sites compared and
    differing, d; module docstring) and its done marker. `meta` must hold `engine` and `d_kind`, what `d` means."""
    missing = [k for k in ("engine", "d_kind") if not meta.get(k)]
    if missing:
        raise ValueError(f"divergence meta lacks {', '.join(missing)}")
    _write_table(Path(path), DIV_COLUMNS, _lines(rows, DIV_COLUMNS), meta)


def read_divergences(path: Path) -> Iterator[dict]:
    """The rows of a pair-divergence file as dicts: n_sites and n_diff int, d float, the eids str."""
    return _read_divergences(Path(path))


def stream_calls(path: Path, fields: tuple[str, ...] = CALL_COLUMNS) -> Iterator[tuple[int, tuple[str, ...]]]:
    """(line number, the raw text of `fields`) of each row of a pair call file (`call_lines`), streamed; a row is
    split only as far as the last field asked for. Fails on what call_lines fails on, and on a short row."""
    path = Path(path)
    cols = [CALL_COLUMNS.index(f) for f in fields]
    nsplit = max(cols) + 1
    whole = nsplit >= len(CALL_COLUMNS)
    for n, line in enumerate(call_lines(path), 2):
        f = line.rstrip("\r\n").split("\t") if whole else line.split("\t", nsplit)
        short = len(f) != len(CALL_COLUMNS) if whole else len(f) <= nsplit
        if short:
            raise ValueError(f"{path}: line {n} has {len(f)} column(s); a pair file has {len(CALL_COLUMNS)}")
        yield n, tuple(f[c] for c in cols)


def call_lines(path: Path) -> Iterator[str]:
    """The data lines of a pair call file, or of a chunk of one, as written (line ends kept), streamed: gzip, or plain
    TSV, told by its first bytes. Fails on a header that is not CALL_COLUMNS or a corrupt file (named: an EOFError
    alone would not say which)."""
    try:
        with open_input(path, encoding="utf-8", newline="") as fh:
            if fh.readline().rstrip("\r\n").split("\t") != list(CALL_COLUMNS):
                raise ValueError(f"{path}: the header is not CALL_COLUMNS: not a pair call file")
            yield from fh
    except (EOFError, gzip.BadGzipFile, zlib.error, UnicodeDecodeError) as exc:
        raise ValueError(f"{path}: corrupt pair file ({type(exc).__name__}: {exc})") from exc


@contextlib.contextmanager
def gz_text(path: Path) -> Iterator[io.TextIOWrapper]:
    """A text handle on a new gzip file `path` (its folder made) with mtime 0 and no stored name, so the same text
    gives the same bytes, written to a temporary file renamed into place on success: a reader never sees a partial
    file."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with atomic_path(path) as tmp, open(tmp, "wb") as raw, \
            gzip.GzipFile(filename="", mode="wb", compresslevel=GZ_LEVEL, fileobj=raw, mtime=0) as gz, \
            io.TextIOWrapper(gz, encoding="utf-8", newline="") as fh:
        yield fh


def _lines(rows: Iterable[dict], columns: tuple[str, ...]) -> Iterator[str]:
    """The rows as TSV lines of `columns` (call_line; a missing key is None)."""
    return (call_line(r.get(c) for c in columns) for r in rows)


def _write_table(path: Path, columns: tuple[str, ...], lines: Iterable[str], meta: dict) -> None:
    if meta.get("engine") not in ENGINES:
        raise ValueError(f"meta['engine'] must be one of {ENGINES}, not {meta.get('engine')!r}")
    clear_done(marker(path))
    with gz_text(path) as fh:
        fh.write("\t".join(columns) + "\n")
        fh.writelines(lines)
    key = meta.get("key") or sha256_file(path)
    mark_done(marker(path), key, [path], {k: v for k, v in meta.items() if k != "key"})


def _read_divergences(path: Path) -> Iterator[dict]:
    try:
        with gzip.open(path, "rt", encoding="utf-8", newline="") as fh:
            if fh.readline().rstrip("\r\n").split("\t") != list(DIV_COLUMNS):
                raise ValueError(f"{path}: the header is not DIV_COLUMNS")
            for n, line in enumerate(fh, 2):
                f = line.rstrip("\r\n").split("\t")
                if len(f) != len(DIV_COLUMNS):
                    raise ValueError(f"{path}: line {n} has {len(f)} column(s); a divergence file has "
                                     f"{len(DIV_COLUMNS)}")
                yield dict(zip(DIV_COLUMNS, map(_parse, DIV_COLUMNS, f), strict=True))
    except (EOFError, gzip.BadGzipFile, zlib.error, UnicodeDecodeError) as exc:
        raise ValueError(f"{path}: corrupt divergence file ({type(exc).__name__}: {exc})") from exc


def _parse(column: str, text: str):
    if text == ".":
        return None
    if column in _INT_COLUMNS:
        return int(text)
    if column in _FLOAT_COLUMNS:
        return float(text)
    return text


def _read_meta(path: Path) -> dict:
    meta = done_meta(marker(path))
    if meta is None:
        raise ValueError(f"{path}: no readable done marker; the file is unfinished or not from a finished run")
    return meta
