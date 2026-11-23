"""Annotation dispatcher: load any supported annotation (GFF3, BED, an LTRquest prefix, one table file) into a
validated, de-duplicated, nested and TSD-annotated ElementTable."""
from __future__ import annotations

import glob
import logging
import zlib
from collections import Counter
from pathlib import Path

import mappy
import numpy as np

from synltr.base.constants import TSD_ORDER
from synltr.base.log import UsageError
from synltr.io.annotation.bed import load_bed
from synltr.io.annotation.gff3 import HOST_RULES, LTR_RULES, NAME_RULES, load_gff3
from synltr.io.annotation.ltrquest import TABLE_RULES, load_ltrquest, load_table
from synltr.io.annotation.nesting import assign_nesting
from synltr.io.annotation.records import CLASS_RULES, TSD_RULES, RecordSkips, check_rule
from synltr.io.annotation.tsd import infer_tsd_at
from synltr.io.elements import Element, ElementTable
from synltr.io.files import open_input
from synltr.io.seqio import Fasta

__all__ = ["AnnotationError", "RecordSkips", "annotation_files", "assign_nesting", "load_annotation",
           "load_annotations", "loader_versions"]

log = logging.getLogger(__name__)

_GFF_EXT = (".gff3", ".gff", ".gff3.gz", ".gff.gz")
_BED_EXT = (".bed", ".bed.gz")
MISSING_RULES = ("error", "skip")  # an annotation seqid absent from the FASTA: LTRpav fails, ltrkin skips


class AnnotationError(UsageError, ValueError):
    """Unusable annotation: unknown format, unreadable file, a missing required column, no usable record, or
    inconsistent with the FASTA."""


def loader_versions(spec: str) -> dict[str, str]:
    """The versions of the libraries the loader of `spec` aligns with, for a stage key: mappy for a BED (its
    LTRs are found by self-alignment, bed.find_ltrs); none for GFF3 and LTRquest tables."""
    return {"mappy": mappy.__version__} if spec.lower().endswith(_BED_EXT) else {}


def annotation_files(spec: str) -> list[Path]:
    """The files `spec` names, as load_annotation reads it: the GFF3/BED file, the LTRquest prefix's
    {prefix}_depth*_clean_ltr.tsv tables (sorted), else the one file `spec` is; [] if there are none."""
    if spec.lower().endswith(_GFF_EXT + _BED_EXT):
        return [Path(spec)] if Path(spec).is_file() else []
    tables = sorted(Path(p) for p in glob.glob(glob.escape(spec) + "_depth*_clean_ltr.tsv"))
    return tables or ([Path(spec)] if Path(spec).is_file() else [])


def _is_gff3(path: Path) -> bool:
    """A file's first nonblank line starts ##gff, or has GFF3's nine tab-separated columns and is no seq_id header."""
    with open_input(path) as fh:
        first = next((line for line in fh if line.strip()), "")
    return first.startswith("##gff") or (len(first.split("\t")) == 9
                                         and not first.lower().lstrip("#").startswith("seq_id"))


def _dispatch(spec: str, genome: str, fasta: Fasta | None, ltr_rule: str, class_rule: str, tsd_rule: str,
              name_family: bool, k_keys: tuple[str, ...], name_rule: str, host_rule: str,
              table_rule: str) -> tuple[list[Element], Counter]:
    low = spec.lower()
    gff3 = (ltr_rule, class_rule, tsd_rule, name_family, k_keys, name_rule, host_rule)
    if low.endswith(_GFF_EXT):
        if not Path(spec).is_file():
            raise AnnotationError(f"annotation file not found: {spec}")
        return load_gff3(Path(spec), genome, *gff3)
    if low.endswith(_BED_EXT):
        if not Path(spec).is_file():
            raise AnnotationError(f"annotation file not found: {spec}")
        if fasta is None:
            raise AnnotationError(f"{spec}: a BED annotation needs the genome FASTA (its LTRs are found in the "
                                  "sequence)")
        return load_bed(Path(spec), genome, fasta)
    if glob.glob(glob.escape(spec) + "_depth*_clean_ltr.tsv"):
        return load_ltrquest(spec, genome, class_rule=class_rule, tsd_rule=tsd_rule, table_rule=table_rule)
    if Path(spec).is_file():
        if _is_gff3(Path(spec)):
            return load_gff3(Path(spec), genome, *gff3)
        return load_table(Path(spec), genome, class_rule=class_rule, tsd_rule=tsd_rule, table_rule=table_rule)
    raise AnnotationError(
        f"cannot detect annotation format of {spec!r}: expected .gff3/.gff/.bed (optionally .gz), "
        "an LTRquest prefix with {prefix}_depth*_clean_ltr.tsv files, or a table file"
    )


def _validate(table: ElementTable, lengths: dict[str, int], genome: str, spec: str, missing_rule: str,
              fasta: Fasta) -> ElementTable:
    """`table` checked against the FASTA's sequence `lengths`: annotation seqids absent from the FASTA fail
    (missing_rule "error") or drop their elements with ltrkin's warning ("skip"); an element ending beyond its
    sequence fails. Returns the table kept."""
    where = f"{genome}: annotation {spec}: "
    missing = sorted({s for s in table.column("seqid") if s not in lengths})
    if missing and missing_rule == "error":
        raise AnnotationError(
            f"{where}{len(missing)} annotation seqid(s) absent from the FASTA: {', '.join(missing[:5])}"
            + (" ..." if len(missing) > 5 else "")
        )
    if missing:
        log.warning("%d annotation sequence(s) not in the reference, skipped: %s", len(missing), missing[:5])
        table = table.select(np.fromiter((s in lengths for s in table.column("seqid")), dtype=bool,
                                         count=len(table)))
        if not len(table):
            raise AnnotationError(f"{where}no annotated element is on a sequence of {fasta.path}")
    seqid, end = table.column("seqid"), table.column("end")
    limit = np.fromiter((lengths[s] for s in seqid), dtype=np.int64, count=len(seqid))
    over = np.flatnonzero(end > limit)
    if len(over):
        e = table.row(int(over[0]))
        raise AnnotationError(
            f"{where}element {e.name} ({e.seqid}:{e.start}-{e.end}) ends beyond its sequence length "
            f"{lengths[e.seqid]}"
        )
    return table


def load_annotation(spec: str, genome: str, fasta: Fasta | None = None, *, ltr_rule: str = "two",
                    class_rule: str = "split", skipped: RecordSkips | None = None,
                    tsd_order: tuple[int, ...] = TSD_ORDER, tsd_rule: str = "keep", name_family: bool = False,
                    k_keys: tuple[str, ...] = ("K2P_d", "ltr_identity"), name_rule: str = "name",
                    missing_rule: str = "error", host_rule: str = "type",
                    table_rule: str = "ltrquest") -> ElementTable:
    """Load, validate, de-duplicate, nest and TSD-annotate one genome's LTR-RT annotation.

    `spec` is a GFF3/BED path, an LTRquest prefix, or another file: GFF3 when its first line says so (or has GFF3's
    nine columns), else a table. `ltr_rule` ("two": exactly two LTRs; "outermost": the outermost of two or more,
    span extended to the LTR ends; "optional": phasing, exactly two LTRs or none known) and `host_rule` ("type":
    the LTR-RT types; "k2p": phasing, LTR_retrotransposon and repeat_region lines with a K2P key) are read by the
    GFF3 loader, `table_rule` ("ltrquest": LTRpav, LTRquest's columns and pass rows; "phasing": phasing's K2P table)
    by the table loader; `class_rule` ("split": superfamily is field 2, clade field 3; "last": clade is the last
    field) and `tsd_rule` ("keep": LTRpav, any annotated TSD; "clean": ltrkin, ACGT 3-8 bp; "acgtn": phasing) by
    the GFF3 and table loaders; `name_family` (ltrkin's family from Name), `k_keys` (the attributes that give k, in
    order) and `name_rule` ("name": LTRpav, the Name attribute first; "id": ltrkin, the ID) by the GFF3 loader.
    `missing_rule` says what an annotation seqid absent from the FASTA does: "error" (LTRpav) fails the load, "skip"
    (ltrkin) drops its elements with one warning line, ltrkin's ("N annotation sequence(s) not in the reference,
    skipped: [first five]"). An unknown rule raises ValueError. A TSD the annotation lacks is inferred from the
    flanks, trying `tsd_order` lengths in turn.
    A malformed record is skipped by its loader with a warning naming its file and line; skip reasons (the
    loader's plus `duplicate`) are logged as one warning line and added to `skipped.counts` if given. Without a
    `fasta` the records are not checked against it and no TSD is inferred; the rows come back sorted by
    (seqid, start, -end) instead of in FASTA order. Structural problems fail with one AnnotationError line: an
    unknown format, an unreadable or corrupt file, a missing required column, records outside the FASTA (on a
    sequence it lacks, under missing_rule "error", or past a sequence's end), or no usable record at all.
    """
    check_rule("ltr_rule", ltr_rule, LTR_RULES)
    check_rule("class_rule", class_rule, CLASS_RULES)
    check_rule("tsd_rule", tsd_rule, TSD_RULES)
    check_rule("name_rule", name_rule, NAME_RULES)
    check_rule("missing_rule", missing_rule, MISSING_RULES)
    check_rule("host_rule", host_rule, HOST_RULES)
    check_rule("table_rule", table_rule, TABLE_RULES)
    try:
        elements, reasons = _dispatch(spec, genome, fasta, ltr_rule, class_rule, tsd_rule, name_family, k_keys,
                                      name_rule, host_rule, table_rule)
    except AnnotationError:
        raise
    except (ValueError, OSError, EOFError, zlib.error) as exc:  # UnicodeDecodeError is a ValueError
        detail = str(exc) if isinstance(exc, (ValueError, FileNotFoundError)) else f"{type(exc).__name__}: {exc}"
        raise AnnotationError(f"{genome}: cannot read annotation {spec}: {detail}") from exc
    reasons = Counter(reasons)
    if not elements:
        why = ", ".join(f"{k}={reasons[k]}" for k in sorted(reasons)) or "no LTR-RT record in it"
        raise AnnotationError(f"{genome}: annotation {spec} holds no usable LTR-RT record ({why})")
    lengths = fasta.lengths() if fasta is not None else None
    table = ElementTable.from_rows(elements)
    if lengths is not None:
        table = _validate(table, lengths, genome, spec, missing_rule, fasta)

    seen: set[tuple[str, int, int]] = set()
    keep = []
    for i, key in enumerate(zip(table.column("seqid"), table.column("start").tolist(), table.column("end").tolist())):
        if key in seen:
            reasons["duplicate"] += 1
            continue
        seen.add(key)
        keep.append(i)
    if len(keep) < len(table):
        table = table.select(keep)

    if reasons:
        detail = ", ".join(f"{k}={reasons[k]}" for k in sorted(reasons))
        log.warning("%s: skipped %d annotation records (%s)", genome, sum(reasons.values()), detail)
        if skipped is not None:
            skipped.counts.update(reasons)

    if fasta is not None:
        table = _fill_tsds(table, fasta, tsd_order)
    nested = assign_nesting(table)
    if lengths is None:
        return nested
    order = {s: i for i, s in enumerate(lengths)}
    seqid, start, end = nested.column("seqid"), nested.column("start"), nested.column("end")
    rows = sorted(range(len(nested)), key=lambda i: (order[seqid[i]], start[i], -end[i]))
    return nested if rows == list(range(len(rows))) else nested.select(rows)


def load_annotations(specs, genome: str, **rules) -> ElementTable:
    """Several annotation files of one genome (each read by load_annotation with `rules`, no FASTA) as one table,
    nested as a whole: a host and its nest annotated in different files (LTRquest's depth tables) keep their parent
    link. One file comes back as load_annotation gives it; none gives an empty table."""
    tables = [load_annotation(str(p), genome, **rules) for p in specs]
    if len(tables) == 1:
        return tables[0]
    return assign_nesting(ElementTable.concat(tables) if tables else ElementTable.from_rows([]))


def _fill_tsds(table: ElementTable, fasta: Fasta, order: tuple[int, ...]) -> ElementTable:
    """Infer the TSD of each element that has none; `tsd_source` records which TSDs were inferred."""
    tsd, source = table.column("tsd").copy(), table.column("tsd_source").copy()
    seqid, start, end = table.column("seqid"), table.column("start").tolist(), table.column("end").tolist()
    for i in (i for i, v in enumerate(tsd) if v is None):
        found = infer_tsd_at(fasta, seqid[i], start[i], end[i], order)
        if found is not None:
            tsd[i], source[i] = found, "inferred"
    return table.replace(tsd=tsd, tsd_source=source)
