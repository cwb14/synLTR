"""loci.tsv, pav_matrix.tsv and qc.tsv. Deterministic text: 1-based inclusive coordinates, "." for none, floats
to 6 significant digits, genomes in sheet order, rows in the order given (build_loci's lid order)."""
from __future__ import annotations

from pathlib import Path

from ltrtrace.io.samples import Sample

from .loci import Locus, QcRow
from .records import tsv_field

_ATTRS = ("lid", "family", "superfamily", "clade", "length", "ltr_len", "tsd")
_PER_GENOME = ("eid", "state", "tier", "seqid", "start", "end", "strand", "k", "flags")


def _join(items) -> str:
    return ",".join(items) if items else "."


def _write(path: Path, rows) -> None:
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.writelines("\t".join(row) + "\n" for row in rows)


def write_loci_tsv(loci: list[Locus], samples: list[Sample], path: Path) -> None:
    """One row per locus: its attributes, then for each genome its member eid and state fields."""
    labels = [s.label for s in samples]
    header = [*_ATTRS, "discovery", "flags", *(f"{g}_{c}" for g in labels for c in _PER_GENOME)]

    def rows():
        yield header
        for loc in loci:
            row = [tsv_field(getattr(loc, a)) for a in _ATTRS]
            row += [_join(loc.discovery), _join(loc.flags)]
            for g in labels:
                st = loc.states[g]
                row += [tsv_field(loc.members.get(g)), st.state, tsv_field(st.tier), tsv_field(st.seqid),
                        tsv_field(st.start), tsv_field(st.end), tsv_field(st.strand), tsv_field(st.k),
                        _join(st.flags)]
            yield row

    _write(Path(path), rows())


def write_pav_matrix(loci: list[Locus], samples: list[Sample], path: Path) -> None:
    """lid x genome matrix of state letters."""
    labels = [s.label for s in samples]
    _write(Path(path), [["lid", *labels], *([loc.lid, *(loc.states[g].state for g in labels)] for loc in loci)])


def write_qc(rows: list[QcRow], path: Path) -> None:
    """QC rows as scope, metric, value."""
    _write(Path(path), [["scope", "metric", "value"], *([s, m, tsv_field(v)] for s, m, v in rows)])
