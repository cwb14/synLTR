"""Load LTRquest `{prefix}_depth*_clean_ltr.tsv` tables, a single such table, or phasing's K2P table into Elements.

Conversions: the `seq_id` header `seqid:start-end#class` is 1-based closed, like the elements (no change). The
columns ltr5_start..ltr3_end are 1-based offsets along the element's genome span in genome orientation, counted from
the span's left end (its seq_id start) on either strand and in every clean mode, whatever the FASTA record's
orientation or length (a minus-strand record is the reverse complement of the span, a cut record is shorter than
it); `strand` only says which strand the element is on, and ltr5 is the genome-left LTR. So a genome coordinate is
`span start + offset - 1` (`_to_genome`) for both strands; io.library places the LTRs in the records. Under
table_rule "phasing" the seq_id is read as a library record name (io.library.parse_record_name)."""
from __future__ import annotations

import glob
import math
import re
from collections import Counter

from ltrtrace.io.annotation.records import (
    CLASS_RULES,
    TSD_RULES,
    RecordSkips,
    check_rule,
    class_fields,
    number,
    own_length,
)
from ltrtrace.io.annotation.tsd import annotated_tsd
from ltrtrace.io.elements import NO_LTR, Element, make_eid
from ltrtrace.io.files import open_input
from ltrtrace.io.library import parse_record_name

_SEQ_ID = re.compile(r"^(?P<seqid>.+):(?P<start>\d+)-(?P<end>\d+)#(?P<cls>.+)$")
_DEPTH = re.compile(r"_depth(\d+)_clean_ltr\.tsv$")
_NULL = {"", "NA", "."}
_CLASS_UNKNOWN = {"unknown"}  # the table maps only this class name to none
_COUNTS = ("n_sites", "n_ts", "n_tv")  # optional numeric columns, read as given
REQUIRED = ("seq_id", "status", "ltr5_start", "ltr5_end", "ltr3_start", "ltr3_end")  # the rest is optional
_BOUNDS = REQUIRED[2:]
_OPTIONAL = ("k2p", "strand", "tsd", "family", "confidence", "ltr5_len", "ltr3_len", "k2p_se") + _COUNTS
TABLE_RULES = ("ltrquest", "phasing")


def _opt(v: str | None) -> str | None:
    return None if v is None or v in _NULL else v


def _to_genome(span_start: int, offset: int) -> int:
    """Genome coordinate (1-based) of a 1-based offset along the genome span that starts at `span_start` (genome
    orientation)."""
    return span_start + offset - 1


def load_ltrquest(prefix: str, genome: str, class_rule: str = "split", tsd_rule: str = "keep",
                  table_rule: str = "ltrquest") -> tuple[list[Element], Counter]:
    """Return (elements, skip counts) of the tables {prefix}_depth<N>_clean_ltr.tsv, in depth order, each read as
    `load_table` reads one."""
    for name, value, allowed in (("class_rule", class_rule, CLASS_RULES), ("tsd_rule", tsd_rule, TSD_RULES),
                                 ("table_rule", table_rule, TABLE_RULES)):
        check_rule(name, value, allowed)
    pattern = f"{glob.escape(prefix)}_depth*_clean_ltr.tsv"
    files = sorted((f for f in glob.glob(pattern) if _DEPTH.search(f)), key=lambda f: int(_DEPTH.search(f).group(1)))
    if not files:
        raise FileNotFoundError(f"no LTRquest files match {prefix}_depth*_clean_ltr.tsv")
    return _load(files, genome, class_rule, tsd_rule, table_rule)


def load_table(path, genome: str, class_rule: str = "split", tsd_rule: str = "keep",
               table_rule: str = "ltrquest") -> tuple[list[Element], Counter]:
    """Return (elements, skip counts) of one table file, plain or gzip (its first bytes tell). `tsd_rule` is read by
    `annotated_tsd`.

    `table_rule` "ltrquest" (LTRpav): the header must hold the REQUIRED columns (ValueError naming the file and the
    missing ones); `class_rule` is read by `class_fields` (only "unknown" is unknown here). A row that cannot be read
    is skipped with a warning naming the file and line: fewer fields than its required columns need
    (malformed_columns), a seq_id not <seqid>:<start>-<end>#<class> (bad_seq_id) or with start > end
    (malformed_coordinates), LTR bounds that are not integers (malformed_ltr_bounds), a k2p that is not a number or
    NA (malformed_k2p). Rows with status other than pass are counted (status_not_pass). A k2p that is not finite is
    unknown.

    "phasing" (phasing's K2P table): the first nonblank line is the header (column names stripped of
    a leading '#' and of blanks, in lower case), which must name a k2p column (ValueError naming the file); later
    lines starting with '#' are comments. The seq_id column, else the first, is a library record name; a row whose
    name gives no interval (bad_seq_id) or start > end (malformed_coordinates) is skipped as above. Every row is kept
    whatever its status. k2p, k2p_se, n_ts and n_tv that are not numbers are unknown; n_sites is the n_sites column,
    or ltr5_len when the table has none; superfamily and clade are not read; the four LTR columns, when all are whole
    numbers, give the LTR extents (as under "ltrquest"), else both extents are NO_LTR.
    """
    for name, value, allowed in (("class_rule", class_rule, CLASS_RULES), ("tsd_rule", tsd_rule, TSD_RULES),
                                 ("table_rule", table_rule, TABLE_RULES)):
        check_rule(name, value, allowed)
    return _load([path], genome, class_rule, tsd_rule, table_rule)


def _load(paths: list, genome: str, class_rule: str, tsd_rule: str, table_rule: str) -> tuple[list[Element], Counter]:
    elements: list[Element] = []
    skips = RecordSkips()
    tsds: dict[str | None, str | None] = {}  # the few thousand distinct TSDs, read once
    if table_rule == "phasing":
        for path in paths:
            _read_phasing(path, genome, tsd_rule, elements, skips, tsds)
        return elements, skips.counts
    classes: dict[str, tuple[str | None, str | None]] = {}  # the few distinct classifications, read once
    for path in paths:
        _read_ltrquest(path, genome, class_rule, tsd_rule, elements, skips, classes, tsds)
    return elements, skips.counts


def _read_phasing(path, genome: str, tsd_rule: str, elements: list[Element], skips: RecordSkips,
                  tsds: dict) -> None:
    with open_input(path) as fh:
        header = None
        for n, line in enumerate(fh, 1):
            if not line.strip():
                continue
            p = line.rstrip("\r\n").split("\t")
            if header is None:
                header = [h.lstrip("#").strip().lower() for h in p]
                if "k2p" not in header:
                    raise ValueError(f"{path}: no 'k2p' column in the header (columns found by name)")
                continue
            if line.startswith("#"):
                continue
            r = dict(zip(header, p))
            seq_name = r.get("seq_id", p[0])
            x = parse_record_name(seq_name)
            if x is None:
                skips.malformed(path, n, "bad_seq_id", f"seq_id {seq_name[:80]!r} is not <seqid>:<start>-<end>")
                continue
            if not 1 <= x.start <= x.end:
                skips.malformed(path, n, "malformed_coordinates", f"seq_id start {x.start} > end {x.end} or < 1")
                continue
            bounds = [r.get(c, "") for c in _BOUNDS]
            if all(b.isdigit() for b in bounds):
                l5s, l5e, l3s, l3e = (int(b) for b in bounds)
                left = (_to_genome(x.start, l5s), _to_genome(x.start, l5e))
                right = (_to_genome(x.start, l3s), _to_genome(x.start, l3e))
            else:
                left = right = NO_LTR
            raw = _opt(r.get("tsd"))
            if raw not in tsds:
                tsds[raw] = annotated_tsd(raw, tsd_rule, _NULL)
            tsd, strand = tsds[raw], r.get("strand")
            elements.append(Element(
                genome=genome, eid=make_eid(genome, x.seqid, x.start, x.end), name=seq_name,
                seqid=x.seqid, start=x.start, end=x.end, strand=strand if strand in ("+", "-") else ".",
                ltr_left=left, ltr_right=right, tsd=tsd, tsd_source="annot" if tsd else None,
                k=number(r.get("k2p")), family=_opt(r.get("family")), superfamily=None, clade=None,
                confidence=_opt(r.get("confidence")),
                ltr_left_len=own_length(r.get("ltr5_len")), ltr_right_len=own_length(r.get("ltr3_len")),
                k_se=number(r.get("k2p_se")), n_sites=number(r["n_sites"] if "n_sites" in r else r.get("ltr5_len")),
                n_ts=number(r.get("n_ts")), n_tv=number(r.get("n_tv")),
            ))


def _read_ltrquest(path, genome: str, class_rule: str, tsd_rule: str, elements: list[Element], skips: RecordSkips,
                   classes: dict, tsds: dict) -> None:
    with open_input(path) as fh:
        header = fh.readline().rstrip("\r\n").lstrip("#").split("\t")
        col = {name: i for i, name in enumerate(header)}
        missing = [c for c in REQUIRED if c not in col]
        if missing:
            raise ValueError(f"{path}: LTRquest table header lacks required column(s) {', '.join(missing)}")
        need = max(col[c] for c in REQUIRED) + 1
        width = len(header)
        ix = {name: col.get(name, -1) for name in _OPTIONAL}  # -1 reads the None appended to each row
        (i_k2p, i_strand, i_tsd, i_family, i_conf, i_l5len, i_l3len, i_se, *i_counts) = (ix[c] for c in _OPTIONAL)
        i_seq, i_status = col["seq_id"], col["status"]
        i_bounds = [col[c] for c in _BOUNDS]
        for n, line in enumerate(fh, 2):
            if not line.strip():
                continue
            r = line.rstrip("\r\n").split("\t")
            if len(r) < need:
                skips.malformed(path, n, "malformed_columns", f"{len(r)} fields, the required columns need "
                                f"{need}")
                continue
            if len(r) < width:
                r += [None] * (width - len(r))
            r.append(None)
            seq_name = r[i_seq]
            m = _SEQ_ID.match(seq_name)
            if not m:
                skips.malformed(path, n, "bad_seq_id", f"seq_id {seq_name[:80]!r} is not "
                                "<seqid>:<start>-<end>#<class>")
                continue
            if r[i_status] != "pass":
                skips.filtered("status_not_pass")
                continue
            seqid, start, end = m["seqid"], int(m["start"]), int(m["end"])
            if not 1 <= start <= end:
                skips.malformed(path, n, "malformed_coordinates", f"seq_id start {start} > end {end} or < 1")
                continue
            try:
                l5s, l5e, l3s, l3e = (int(r[i]) for i in i_bounds)
            except ValueError:
                skips.malformed(path, n, "malformed_ltr_bounds", "LTR bounds " + ", ".join(
                    f"{c} {r[col[c]]!r}" for c in _BOUNDS) + " are not all integers")
                continue
            k = _opt(r[i_k2p])
            try:
                k = float(k) if k is not None else None
            except ValueError:
                skips.malformed(path, n, "malformed_k2p", f"k2p {k[:40]!r} is not a number or NA")
                continue
            cls = m["cls"]
            fields = classes.get(cls)
            if fields is None:
                fields = classes[cls] = class_fields(cls, class_rule, unknown=_CLASS_UNKNOWN)
            sf, clade = fields
            strand = r[i_strand]
            raw = _opt(r[i_tsd])
            if raw in tsds:
                tsd = tsds[raw]
            else:
                tsd = tsds[raw] = annotated_tsd(raw, tsd_rule, _NULL)
            elements.append(Element(
                genome=genome,
                eid=make_eid(genome, seqid, start, end),
                name=seq_name,
                seqid=seqid, start=start, end=end,
                strand=strand if strand in ("+", "-") else ".",
                ltr_left=(_to_genome(start, l5s), _to_genome(start, l5e)),
                ltr_right=(_to_genome(start, l3s), _to_genome(start, l3e)),
                tsd=tsd, tsd_source="annot" if tsd else None,
                k=k if k is not None and math.isfinite(k) else None,
                family=_opt(r[i_family]), superfamily=sf, clade=clade,
                confidence=_opt(r[i_conf]),
                ltr_left_len=own_length(r[i_l5len]), ltr_right_len=own_length(r[i_l3len]),
                k_se=number(r[i_se]), **{c: number(r[i]) for c, i in zip(_COUNTS, i_counts)},
            ))
