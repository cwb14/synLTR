"""Read element bounds from BED and recover each element's LTRs by self-alignment.

BED is 0-based half-open; Element coordinates are 1-based closed. The two conversions, `from_bed` (on read) and
`to_bed` (on write), are the only place the shift happens.
"""
from __future__ import annotations

from collections import Counter
from pathlib import Path

import mappy

from synltr.base.constants import BED_LTR_EDGE, BED_LTR_MAX, COV_PRESENT, MIN_SEG
from synltr.io.annotation.records import RecordSkips
from synltr.io.elements import Element, ElementTable, make_eid
from synltr.io.files import open_input
from synltr.io.seqio import Fasta

_SKIP_PREFIX = ("#", "track", "browser")


def from_bed(start: int, end: int) -> tuple[int, int]:
    """A BED interval (0-based half-open) as an element's (start, end), 1-based closed."""
    return start + 1, end


def to_bed(start: int, end: int) -> tuple[int, int]:
    """An element's (start, end), 1-based closed, as a BED interval (0-based half-open)."""
    return start - 1, end


def write_bed(elements: ElementTable, path: Path) -> None:
    """Write BED6 (chrom, start, end, name, score 0, strand), one line per element."""
    with open(path, "w") as fh:
        for e in elements:
            start, end = to_bed(e.start, e.end)
            fh.write(f"{e.seqid}\t{start}\t{end}\t{e.name}\t0\t{e.strand}\n")


def find_ltrs(seq: str) -> tuple[int, int] | None:
    """Return (left, right) LTR lengths by aligning the element start to its end, or None."""
    n = min(len(seq) // 2, BED_LTR_MAX)
    if n < MIN_SEG:
        return None
    aligner = mappy.Aligner(seq=seq[-n:], preset="asm20")
    if not aligner.seq_names:  # index failed to build
        return None
    best = None
    for h in aligner.map(seq[:n]):
        if (h.strand == 1 and h.q_st <= BED_LTR_EDGE and h.r_en >= n - BED_LTR_EDGE
                and h.blen > 0 and h.mlen / h.blen >= COV_PRESENT
                and h.q_en >= MIN_SEG and n - h.r_st >= MIN_SEG
                and (best is None or h.mlen > best.mlen)):
            best = h
    return None if best is None else (best.q_en, n - best.r_st)


def _read(path: Path, skips: RecordSkips) -> list[tuple[int, str, int, int, str, str]]:
    """(line, chrom, start, end, name, strand) of each record; a line with fewer than 3 columns
    (malformed_columns) or a start and end that are not integers with 0 <= start < end
    (malformed_coordinates) is skipped with a warning naming the file and line."""
    rows = []
    with open_input(path) as fh:
        for n, line in enumerate(fh, 1):
            if not line.strip() or line.startswith(_SKIP_PREFIX):
                continue
            f = line.rstrip("\r\n").split("\t")
            if len(f) < 3:
                skips.malformed(path, n, "malformed_columns", f"{len(f)} tab-separated column(s), BED needs 3")
                continue
            try:
                start, end = int(f[1]), int(f[2])
            except ValueError:
                start = end = -1
            if not 0 <= start < end:
                skips.malformed(path, n, "malformed_coordinates", f"start {f[1]!r} and end {f[2]!r} are not "
                                "integers with 0 <= start < end")
                continue
            chrom = f[0]
            name = f[3] if len(f) > 3 and f[3] else f"{chrom}:{from_bed(start, end)[0]}-{end}"
            strand = f[5] if len(f) > 5 and f[5] in ("+", "-") else "."
            rows.append((n, chrom, start, end, name, strand))
    return rows


def load_bed(path: Path, genome: str, fasta: Fasta) -> tuple[list[Element], Counter]:
    """Return (elements, skip counts). Skip reasons: no_ltrs, and the malformed records of `_read`. Raises
    ValueError on records outside the FASTA (an unknown seqid, or an end past the sequence: another genome's
    annotation), naming the first ones' lines."""
    path = Path(path)
    skips = RecordSkips()
    rows = _read(path, skips)
    lengths = fasta.lengths()
    bad = [r for r in rows if r[1] not in lengths or r[3] > lengths[r[1]]]
    if bad:
        shown = "; ".join(f"line {n} {c}:{s}-{e}" for n, c, s, e, _, _ in bad[:5])
        raise ValueError(f"{path}: {len(bad)} BED record(s) not in the FASTA (unknown seqid or out of range); "
                         f"first: {shown}")
    elements: list[Element] = []
    for _, chrom, start, end, name, strand in rows:
        ltrs = find_ltrs(fasta.fetch(chrom, start, end))
        if ltrs is None:
            skips.filtered("no_ltrs")
            continue
        l5, l3 = ltrs
        s, _ = from_bed(start, end)
        elements.append(Element(
            genome=genome, eid=make_eid(genome, chrom, s, end), name=name,
            seqid=chrom, start=s, end=end, strand=strand,
            ltr_left=(s, start + l5), ltr_right=(end - l3 + 1, end),
            tsd=None, k=None, family=None, superfamily=None, clade=None, confidence=None,
        ))
    return elements, skips.counts
