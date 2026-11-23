"""LTR-RT libraries: the FASTA records of one genome's LTR-RT copies joined to its annotation, as one columnar store
(`Library`), and the reading of record names.

Coordinates (each conversion below is made here and nowhere else):
- A record name `seqid:start-end` or `seqid:start..end` (LTRquest, EDTA; the seqid may be PanSN, HG002#1#chr1; an
  optional (+), (-) or (.) and a #class may follow) is 1-based closed, like the elements. A bedtools name
  `name::seqid:start-end` has a 0-based start and is read as start + 1 (`parse_record_name`).
- Positions in a record (`offs`, `l5`, `l3`) are 0-based half-open.
- `record_offsets` places genome intervals (1-based closed) of an element in its record: a minus-strand record is the
  reverse complement of the element's genome span, so a genome position p sits at end - p; any other record is in
  genome orientation, p at p - start. LTR offsets an annotation counts along the span (LTRquest's table, on either
  strand) reach this module as genome intervals (io.annotation.ltrquest).
- `record_ends` places the two LTRs at the ends of a record shorter than its genome span (nests cut out), from the
  LTRs' own lengths. Under the screen's rule it places them on every record, at the ends of the N-stripped record, and
  l5 and l3 are offsets into that stripped record.
- A record whose name gives no interval (kept, unplaced, under the screen's rule only) has no coordinates at all: no
  eid, seqid None, start and end 0; nothing writes it to a file of elements, calls or BED rows.

`library_files`, `read_screen_tables` and `read_library` read a samples sheet's library spec (an LTRquest prefix, a
FASTA, or FASTA and tables) for the library screen, its tables under the screen's table rule (LTRmix's).
"""
from __future__ import annotations

import glob
import logging
import re
import zlib
from dataclasses import dataclass
from pathlib import Path
from typing import NamedTuple

import numpy as np
from numba import njit, prange

from synltr.base.constants import CLEAN_CAL, MIN_MATCH, UNPLACED_MAX
from synltr.base.log import UsageError
from synltr.io.elements import NO_LTR, Element, ElementTable, make_eid
from synltr.io.files import gzip_kind, open_input
from synltr.io.seqio import encode, read_fai, stream_fasta
from synltr.seq.divergence import clean_tolerance, self_ltr_z

log = logging.getLogger(__name__)

RULES = ("phasing", "screen")
_NAME = re.compile(r"^(?P<bed>.*::)?(?P<seqid>.+):(?P<start>[0-9]+)(?:-|\.\.)(?P<end>[0-9]+)(?:\([-+.]\))?"
                   r"(?:#(?P<cls>.*))?$")


class LibraryError(UsageError, ValueError):
    """Unusable library: an unreadable FASTA, no record with an interval in its name, too many records repeating an
    interval, or too few records annotated."""


class RecordName(NamedTuple):
    """The genome interval a record name gives (1-based closed) and its class (the text after '#'; None without)."""

    seqid: str
    start: int
    end: int
    cls: str | None


def parse_record_name(token: str) -> RecordName | None:
    """The interval of a record name: `seqid:start-end` or `seqid:start..end`, 1-based closed (LTRquest, EDTA), the
    seqid possibly PanSN (`HG002#1#chr1`), then an optional `(+)`, `(-)` or `(.)` and an optional `#class`. A bedtools
    name `name::seqid:start-end` has a 0-based start: start + 1 is returned. None when the name holds no interval."""
    m = _NAME.match(token)
    if m is None:
        return None
    return RecordName(m["seqid"], int(m["start"]) + (1 if m["bed"] else 0), int(m["end"]), m["cls"])


def _in_record(a: int, b: int, start: int, end: int, strand: str) -> tuple[int, int]:
    """The record offsets (0-based half-open) of the genome interval a..b (1-based closed) of the element at
    start..end: counted back from end on a minus-strand record (the reverse complement of the span), else from
    start."""
    return (end - b, end - a + 1) if strand == "-" else (a - start, b - start + 1)


def record_offsets(ivs, start: int, end: int, strand: str) -> tuple[tuple[int, int], tuple[int, int]]:
    """Genome intervals (1-based closed) of the two LTRs of the element at start..end -> their offsets in its record
    (0-based half-open), 5' LTR first. Minus-strand records are reverse complements (LTRquest writes element
    orientation); unknown-strand records are in genome orientation."""
    out = sorted(_in_record(a, b, start, end, strand) for a, b in ivs)
    return out[0], out[1]


def record_ends(left_len, right_len, strand: str, length: int) -> tuple[tuple[int, int], tuple[int, int]] | None:
    """The two LTRs of a record shorter than its genome span (nests cut out, so genome LTR intervals no longer map into
    it), at its ends: `left_len` and `right_len` are the own lengths of the genome-left and genome-right LTRs
    (LTRquest's ltr5_len and ltr3_len on either strand, nests removed); a minus-strand record (element orientation)
    starts with the genome-right LTR. Offsets in the record of `length` bases (0-based half-open, 5' LTR first); None
    when a length is unknown (None or NaN) or not positive, or the two do not fit."""
    if left_len is None or right_len is None or not (np.isfinite(left_len) and np.isfinite(right_len)) \
            or left_len <= 0 or right_len <= 0:
        return None
    first, last = (int(right_len), int(left_len)) if strand == "-" else (int(left_len), int(right_len))
    if first + last > length:
        return None
    return (0, first), (length - last, length)


@dataclass(frozen=True, eq=False)
class Library:
    """One genome's LTR-RT records in FASTA order, columnar. Record offsets are 0-based half-open."""

    names: np.ndarray  # object: each record's name, the first word of its header line
    cls: np.ndarray  # object: the class its name gives (after '#'); None without one
    codes: np.ndarray  # uint8: the records, concatenated; A C G T (either case) = 0-3, any other letter 4
    offs: np.ndarray  # int64, n + 1 boundaries: record i is codes[offs[i]:offs[i + 1]]
    orient: np.ndarray  # int8: +1 plus strand; -1 minus strand (record = reverse complement of the span); 0 unknown
    #                     strand (record in genome orientation)
    l5: np.ndarray  # int64 (n, 2): offsets of the 5' LTR in the record (element orientation; the screen's rule: in
    #                 the N-stripped record); -1 = none
    l3: np.ndarray  # int64 (n, 2): the 3' LTR, likewise
    bounds: np.ndarray  # str: how l5 and l3 were found: "mapped" (genome LTR intervals), "ends" (record ends), "none"
    clean: np.ndarray  # str: "cut" (shorter than its genome span: nests cut out), "masked" (nests N-masked), "none"
    placed: np.ndarray  # bool: the name places the record in the genome; an unplaced record keeps its sequence but
    #                     has no eid and no BED row
    elements: ElementTable  # one row per record; a record its annotation lacks has a row from its name (no K2P,
    #                         strand '.', LTRs NO_LTR)

    @property
    def n(self) -> int:
        return len(self.offs) - 1

    def seq(self, i: int, *, stripped: bool = False) -> np.ndarray:
        """The codes of record i; `stripped` drops every base that is not A, C, G or T (code 4)."""
        r = self.codes[self.offs[i]:self.offs[i + 1]]
        return r[r < 4] if stripped else r


class _Codes:
    """Record codes, concatenated: written into an array of `size` codes allocated up front (`_capacity`), or gathered
    and joined at the end when there is no size or the records outgrow it."""

    def __init__(self, size: int | None):
        self._buf = np.empty(size or 0, np.uint8)
        self._used = 0
        self._rest: list[np.ndarray] = []

    def add(self, codes: np.ndarray) -> None:
        end = self._used + codes.size
        if not self._rest and end <= self._buf.size:
            self._buf[self._used:end] = codes
            self._used = end
        else:
            self._rest.append(codes)

    def codes(self) -> np.ndarray:
        head = self._buf[:self._used]
        return np.concatenate([head, *self._rest]) if self._rest else head


def _fai_total(fastas: list[Path]) -> int | None:
    """The bases the FASTAs' .fai indexes list in all; None unless every FASTA has a readable one beside it (the total
    only sizes the code array: without it the array is sized from the files)."""
    total = 0
    for p in fastas:
        fai = Path(f"{p}.fai")
        if not fai.is_file():
            return None
        try:
            total += sum(read_fai(fai).values())
        except ValueError:
            return None
    return total


def _capacity(fastas: list[Path]) -> int | None:
    """Room for the FASTAs' codes: the .fai total when every FASTA has one, else the files' sizes when none is gzip
    (an upper bound: headers and line ends hold no bases); None for compressed FASTAs without .fai."""
    total = _fai_total(fastas)
    if total is not None:
        return total
    size = 0
    for p in fastas:
        if gzip_kind(p):
            return None
        size += p.stat().st_size
    return size


def _read_records(fastas: list[Path], keep_unnamed: bool = False):
    """The records whose name gives an interval, the first of each interval, and with `keep_unnamed` those whose name
    gives none (interval None): (names, classes, intervals, lengths, codes, records without an interval, records
    repeating an interval)."""
    names, cls, keys, lens = [], [], [], []
    sink = _Codes(_capacity(fastas))
    seen = set()
    unnamed = repeats = 0
    for path in fastas:
        try:
            for name, seq in stream_fasta(path):
                x = parse_record_name(name)
                if x is None:
                    unnamed += 1
                    if not keep_unnamed:
                        continue
                    key = None
                else:
                    key = (x.seqid, x.start, x.end)
                    if key in seen:
                        repeats += 1
                        continue
                    seen.add(key)
                names.append(name)
                cls.append(None if x is None else x.cls)
                keys.append(key)
                lens.append(len(seq))
                sink.add(encode(seq))
        except (ValueError, OSError, EOFError, zlib.error) as exc:  # UnicodeDecodeError is a ValueError
            raise LibraryError(f"cannot read library FASTA {path}: {exc}") from exc
    return names, cls, keys, np.array(lens, np.int64), sink.codes(), unnamed, repeats


def _objects(values: list) -> np.ndarray:
    out = np.empty(len(values), dtype=object)
    out[:] = values
    return out


def _join(elements: ElementTable, names: list[str], keys: list[tuple | None], min_match: float,
          genome: str | None = None) -> ElementTable:
    """One annotation row per record, matched by seqid:start-end (the first row of an interval); a record without one
    gets a row from its name (genome `genome`, default the annotation's), and a record without an interval (key None)
    a row without eid or coordinates (seqid None, start and end 0). Fewer than min_match of the placed records
    annotated is a LibraryError."""
    row: dict[tuple, int] = {}
    seqid, start, end = elements.column("seqid"), elements.column("start").tolist(), elements.column("end").tolist()
    for j, key in enumerate(zip(seqid, start, end)):
        row.setdefault(key, j)
    idx = [-1 if k is None else row.get(k, -1) for k in keys]
    hit, placed = sum(j >= 0 for j in idx), sum(k is not None for k in keys)
    if hit < min_match * placed:
        miss = [f"{k[0]}:{k[1]}-{k[2]}" for k, j in zip(keys, idx) if j < 0 and k is not None][:5]
        raise LibraryError(f"only {hit} of {placed} FASTA copies are annotated (need {min_match:.0%}); unmatched "
                           f"e.g. {', '.join(miss)} (annotation and FASTA must share seqid:start-end coordinates)")
    if hit == len(keys):
        return elements.select(np.array(idx, np.int64))
    if genome is None:
        genome = str(elements.column("genome")[0]) if len(elements) else ""
    extra = [Element(genome=genome, eid=None if k is None else make_eid(genome, *k), name=names[i],
                     seqid=None if k is None else k[0], start=0 if k is None else k[1], end=0 if k is None else k[2],
                     strand=".", ltr_left=NO_LTR, ltr_right=NO_LTR, tsd=None, k=None, family=None,
                     superfamily=None, clade=None, confidence=None)
             for i, (k, j) in enumerate(zip(keys, idx)) if j < 0]
    table = ElementTable.concat([elements, ElementTable.from_rows(extra)])
    pick = np.array(idx, np.int64)
    pick[pick < 0] = np.arange(len(elements), len(elements) + len(extra))
    return table.select(pick)


def _clean_kinds(table: ElementTable, elements: ElementTable, codes: np.ndarray, offs: np.ndarray,
                 placed: np.ndarray) -> np.ndarray:
    """Each record's clean kind: "cut" when it is shorter than its genome span; "masked" when it is as long and every
    nest of its element (an annotation row whose parent the element is) reads N (code 4) in it; else "none", and
    "none" for a record without an interval."""
    lens = np.diff(offs)
    start, end, eid, strand = (table.column(c) for c in ("start", "end", "eid", "strand"))
    span = end - start + 1
    out = np.full(len(table), "none", dtype="<U6")
    out[(lens < span) & placed] = "cut"
    nests: dict[str, list[tuple[int, int]]] = {}
    for a, b, parent in zip(elements.column("start").tolist(), elements.column("end").tolist(),
                            elements.column("parent")):
        if parent is not None:
            nests.setdefault(parent, []).append((a, b))
    for i in np.flatnonzero((lens == span) & placed):
        inner = nests.get(eid[i])
        if inner:
            r = codes[offs[i]:offs[i + 1]]
            holes = (_in_record(a, b, int(start[i]), int(end[i]), strand[i]) for a, b in inner)
            if all((r[x:y] == 4).all() for x, y in holes):
                out[i] = "masked"
    return out


def _phasing_bounds(lib: Library, threads: int) -> None:
    """Phasing's LTR bounds, set in place: a copy's genome LTR intervals mapped into its record when they fit there
    ("mapped"); else, for a record shorter than its genome span, its two record ends (`record_ends`), kept ("ends")
    when their own K2P reproduces the annotated K2P within the tolerance `clean_tolerance` calibrates on up to
    CLEAN_CAL copies with mapped LTRs (drawn with seed 0). LTRquest's LTR lengths can still include a nest inside an
    LTR; there the record-end windows reach into the internal region and their divergence disagrees (4-6% of clean
    records). A copy left without bounds is used whole."""
    e = lib.elements
    left, right = e.column("ltr_left").tolist(), e.column("ltr_right").tolist()
    start, end, strand = e.column("start").tolist(), e.column("end").tolist(), e.column("strand")
    n5, n3 = e.column("ltr_left_len"), e.column("ltr_right_len")
    lens = np.diff(lib.offs).tolist()
    n_bad, at_ends = 0, []
    for i in range(lib.n):
        if tuple(left[i]) == NO_LTR:
            continue
        o5, o3 = record_offsets((left[i], right[i]), start[i], end[i], strand[i])
        L = lens[i]
        if 0 <= o5[0] < o5[1] <= o3[0] < o3[1] <= L:
            lib.l5[i], lib.l3[i] = o5, o3
            lib.bounds[i] = "mapped"
        elif L < end[i] - start[i] + 1 and (ends := record_ends(n5[i], n3[i], strand[i], L)) is not None:
            lib.l5[i], lib.l3[i] = ends
            at_ends.append(i)
        else:
            n_bad += 1
    if at_ends:
        at_ends = np.array(at_ends, np.int64)
        k = e.column("k")
        mapped = np.setdiff1d(np.flatnonzero((lib.l5[:, 0] >= 0) & np.isfinite(k)), at_ends)
        if mapped.size > CLEAN_CAL:
            mapped = np.sort(np.random.default_rng(0).choice(mapped, CLEAN_CAL, replace=False))
        tol = clean_tolerance(self_ltr_z(lib, mapped, k, threads))
        ok = self_ltr_z(lib, at_ends, k, threads) <= tol
        lib.l5[at_ends[~ok]] = lib.l3[at_ends[~ok]] = -1
        lib.bounds[at_ends[ok]] = "ends"
        n_bad += int((~ok).sum())
        log.info("%d copies are shorter than their genome span (nested insertions cut out): LTRs at the record ends "
                 "(%d more whose record-end LTRs do not reproduce the annotated LTR divergence; tolerance z <= %.2f "
                 "from %d copies with mapped LTRs)", int(ok.sum()), int((~ok).sum()), tol, mapped.size)
    if n_bad:
        log.warning("%d copies have LTR bounds outside their record; whole elements used for them", n_bad)


def _screen_elements(names: list[str], keys: list[tuple | None], rows: dict[str, dict[str, str]],
                     genome: str) -> ElementTable:
    """One element per record, in record order, from its name (interval, eid) and the row of LTRmix's table whose
    seq_id equals the name exactly (none: unannotated): k2p from the column k2p, else K2P_d (a column K2P is not read),
    a finite number or unknown; strand + or -, else '.'; the TSD as written (upper case) when it is A, C, G and T only;
    ltr5_len and ltr3_len truncated to integers, unknown when not finite numbers. Elements are nested over the
    library's placed records (assign_nesting); a record without an interval has no eid, no coordinates, no parent."""
    from synltr.io.annotation.nesting import assign_nesting
    from synltr.io.annotation.records import number
    from synltr.io.annotation.tsd import canon_tsd

    def length(v):
        x = number(v)
        return None if x is None else int(x)

    out = []
    for name, key in zip(names, keys):
        r = rows.get(name) or {}
        st, raw = r.get("strand", "."), r.get("tsd", "")
        tsd = raw.upper() if canon_tsd(raw) else None
        out.append(Element(
            genome=genome, eid=None if key is None else make_eid(genome, *key), name=name,
            seqid=None if key is None else key[0], start=0 if key is None else key[1], end=0 if key is None else key[2],
            strand=st if st in ("+", "-") else ".", ltr_left=NO_LTR, ltr_right=NO_LTR, tsd=tsd,
            tsd_source="annot" if tsd else None, k=number(r.get("k2p", r.get("K2P_d"))), family=None,
            superfamily=None, clade=None, confidence=None, ltr_left_len=length(r.get("ltr5_len")),
            ltr_right_len=length(r.get("ltr3_len"))))
    table = ElementTable.from_rows(out)
    placed = np.flatnonzero([k is not None for k in keys])
    if placed.size:
        nested = assign_nesting(table.select(placed))
        at = {e: i for i, e in enumerate(nested.column("eid"))}
        order = np.array([at[e] for e in table.column("eid")[placed]], np.int64)
        parent, level = table.column("parent").copy(), table.column("level").copy()
        parent[placed], level[placed] = nested.column("parent")[order], nested.column("level")[order]
        table = table.replace(parent=parent, level=level)
    return table


@njit(parallel=True, cache=True)
def _acgt_lengths(codes, offs):
    """Each record's count of A, C, G and T codes (its N-stripped length)."""
    out = np.zeros(offs.size - 1, np.int64)
    for i in prange(offs.size - 1):
        c = 0
        for p in range(offs[i], offs[i + 1]):
            if codes[p] < 4:
                c += 1
        out[i] = c
    return out


def _screen_bounds(lib: Library) -> None:
    """The screen's LTR bounds, set in place: on every record, its two LTRs at the ends of the N-stripped record (every
    base but A, C, G and T removed: the two clean conventions give one sequence), from the LTRs' own lengths
    (`record_ends`), offsets into `seq(i, stripped=True)`, bounds "ends"; none when a length is unknown or the two do
    not fit."""
    e = lib.elements
    n5, n3, strand = e.column("ltr_left_len"), e.column("ltr_right_len"), e.column("strand")
    stripped = _acgt_lengths(lib.codes, lib.offs).tolist()
    for i in range(lib.n):
        ends = record_ends(n5[i], n3[i], strand[i], stripped[i])
        if ends is not None:
            lib.l5[i], lib.l3[i] = ends
            lib.bounds[i] = "ends"


def load_library(fastas, elements, *, rule: str = "phasing", min_match: float = MIN_MATCH, threads: int = 1,
                 genome: str | None = None) -> Library:
    """The records of the library FASTAs (plain, gzip or bgzip, streamed; the code array is allocated from their .fai
    when each has one) joined to the genome's annotation `elements`.

    Both rules drop a record repeating an earlier record's interval (the first is kept; more than 1 - min_match of the
    records repeating fails). rule "phasing": `elements` is an ElementTable (io.annotation), joined by seqid:start-end;
    a record the annotation lacks gets a row from its name, under `genome` (default the annotation's); a record whose
    name gives no interval is dropped (one warning); fewer than min_match of the records annotated fails, naming some;
    LTR bounds are phasing's (`_phasing_bounds`, alignments in up to `threads` processes). Every record is placed.
    rule "screen" (LTRmix's reading): `elements` is the rows of the library's tables (`read_screen_tables`), joined by
    the exact record name (`_screen_elements`; no annotated share is required; LTRmix's warning when rows exist and no
    record gets a K2P); a record whose name gives no interval is kept unplaced (no eid, no coordinates; one warning),
    more than UNPLACED_MAX of the records so fails; LTR bounds at the N-stripped record's ends on every record
    (`_screen_bounds`: l5 and l3 are then offsets into seq(i, stripped=True)). A failure is one LibraryError line; an
    unknown rule raises ValueError.
    """
    if rule not in RULES:
        raise ValueError(f"rule must be one of {RULES}, not {rule!r}")
    fastas = [Path(p) for p in fastas]
    screen = rule == "screen"
    names, cls, keys, lens, codes, unnamed, repeats = _read_records(fastas, keep_unnamed=screen)
    if unnamed and not screen:
        log.warning("skipped %d FASTA records whose name has no seqid:start-end", unnamed)
    if not names:
        raise LibraryError("no usable LTR-RT records in " + ", ".join(map(str, fastas)))
    if unnamed and screen:
        if unnamed > UNPLACED_MAX * len(names):
            raise LibraryError(f"{unnamed} of {len(names)} FASTA records have no seqid:start-end in their name (a "
                               "library needs LTR_peel-style names, chrom:start-end)")
        log.warning("%d of %d FASTA records have no seqid:start-end; kept unplaced (no element, no call)", unnamed,
                    len(names))
    total = len(names) + repeats
    if repeats > (1 - min_match) * total:
        raise LibraryError(f"{repeats} of {total} FASTA records repeat the coordinates of an earlier record")
    if repeats:
        log.warning("%d FASTA records repeat the coordinates of an earlier record; the first is kept", repeats)
    if screen:
        table = _screen_elements(names, keys, elements, genome or "")
        if elements and not np.isfinite(table.column("k")).any():
            log.warning("%s: no K2P values: age strata disabled for this genome", genome)
    else:
        table = _join(elements, names, keys, min_match, genome)
    n = len(names)
    offs = np.concatenate(([0], np.cumsum(lens))).astype(np.int64)
    orient = np.array([1 if s == "+" else -1 if s == "-" else 0 for s in table.column("strand")], np.int8)
    placed = np.array([k is not None for k in keys], bool)
    lib = Library(names=_objects(names), cls=_objects(cls), codes=codes, offs=offs, orient=orient,
                  l5=np.full((n, 2), -1, np.int64), l3=np.full((n, 2), -1, np.int64),
                  bounds=np.full(n, "none", dtype="<U6"),
                  clean=_clean_kinds(table, table if screen else elements, codes, offs, placed), placed=placed,
                  elements=table)
    if screen:
        _screen_bounds(lib)
    else:
        _phasing_bounds(lib, threads)
    return lib


_LIBRARY_FASTA = (".fa", ".fasta", ".fa.gz", ".fasta.gz", ".fna")
_DEPTH_FA = re.compile(r"_depth(\d+)_clean_ltr\.fa$")


def library_files(spec: str) -> tuple[list[Path], list[Path]]:
    """The FASTAs and tables of a samples sheet's library: "FASTA,TABLE[,TABLE ...]", a FASTA alone when it ends in
    .fa, .fasta, .fa.gz, .fasta.gz or .fna, else an LTRquest prefix: its <prefix>_depth<N>_clean_ltr.fa in depth order,
    each with the .tsv beside it. A missing file, or a prefix with no library, is a LibraryError naming it."""
    if "," in spec or spec.endswith(_LIBRARY_FASTA):
        fa, *names = spec.split(",")
        fastas, tables = [Path(fa)], [Path(t) for t in names]
    else:
        found = [(int(m[1]), p) for p in glob.glob(glob.escape(spec) + "_depth*_clean_ltr.fa")
                 if (m := _DEPTH_FA.search(p))]
        if not found:
            raise LibraryError(f"no LTR_peel libraries {spec}_depth*_clean_ltr.fa")
        fastas = [Path(p) for _, p in sorted(found)]
        tables = [p.with_suffix(".tsv") for p in fastas]
    missing = [str(p) for p in (*fastas, *tables) if not p.is_file()]
    if missing:
        raise LibraryError(f"missing library file(s): {', '.join(missing)}")
    return fastas, tables


def read_screen_tables(paths) -> dict[str, dict[str, str]]:
    """The screen's table reading (LTRmix's): each table (plain or gzip, told by its first bytes; read as text, so
    CRLF line ends read as LF) has its header on the first line, a leading '#' stripped, names as written (case
    kept); every later row becomes a dict of those names keyed by its first field (seq_id); a seq_id repeated, in one
    table or across tables, keeps its last row."""
    rows: dict[str, dict[str, str]] = {}
    for path in paths:
        try:
            with open_input(path) as fh:
                header = fh.readline().rstrip("\n").lstrip("#").split("\t")
                for line in fh:
                    fields = line.rstrip("\n").split("\t")
                    rows[fields[0]] = dict(zip(header, fields))
        except (OSError, EOFError, UnicodeDecodeError, zlib.error) as exc:
            raise LibraryError(f"cannot read library table {path}: {exc}") from exc
    return rows


def read_library(spec: str, genome: str, *, threads: int = 1) -> Library:
    """Genome `genome`'s library for the screen from a samples sheet's spec (`library_files`): its FASTAs joined, under
    rule "screen", to the rows of its tables (`read_screen_tables`; a FASTA alone: every record unannotated)."""
    fastas, tables = library_files(spec)
    return load_library(fastas, read_screen_tables(tables), rule="screen", threads=threads, genome=genome)
