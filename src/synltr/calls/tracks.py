"""Per-genome tracks: <G>.pav.gff3 and <G>.pav.bed. A genome's rows are the loci where it has coordinates and
state P, S, T, E or D, in its FASTA index order. GFF3 is 1-based inclusive, its source column GFF3_SOURCE (synltr
and its version), seqids and attribute values escaped as GFF3 requires (io.annotation.gff3: escape_seqid,
escape_value); BED is 0-based half-open (BED6 + state, tier, family; a tab, line break or % in a family name is
%XX)."""
from __future__ import annotations

import re
from pathlib import Path

from synltr.base.constants import GFF3_SOURCE
from synltr.io.annotation.gff3 import escape_seqid, escape_value
from synltr.io.samples import Sample

from .loci import Locus
from .states import OBSERVED_ABSENT, PRESENT

_TYPE = {**dict.fromkeys(PRESENT, "LTR_retrotransposon"), **dict.fromkeys(OBSERVED_ABSENT, "insertion_site"),
         "D": "deletion_junction"}
_BED_BAD = re.compile(r"[%\t\n\r\x00-\x1f\x7f]")  # BED columns: tabs and line breaks (and %) as %XX


def _bed(s: str) -> str:
    return _BED_BAD.sub(lambda m: f"%{ord(m.group()):02X}", s)


def write_tracks(loci: list[Locus], samples: list[Sample], outdir: Path, seq_order: dict[str, dict[str, int]],
                 seq_lengths: dict[str, dict[str, int]]) -> None:
    """Write <G>.pav.gff3 and <G>.pav.bed into `outdir` for every genome G. seq_order: genome -> seqid ->
    fai index; seq_lengths: genome -> seqid -> length (for the ##sequence-region lines)."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    labels = [s.label for s in samples]
    for g in labels:
        order = seq_order[g]
        rows = []
        for loc in loci:
            st = loc.states[g]
            if st.state in _TYPE and st.seqid is not None and st.start is not None and st.end is not None:
                if st.seqid not in order or st.seqid not in seq_lengths[g]:
                    raise ValueError(f"{loc.lid}: {g}'s seqid {st.seqid!r} is not in {g}'s FASTA index")
                rows.append((order[st.seqid], st.start, loc.lid, loc))
        rows.sort(key=lambda t: t[:3])
        used = sorted({loc.states[g].seqid for *_, loc in rows}, key=order.__getitem__)
        with open(outdir / f"{g}.pav.gff3", "w", encoding="utf-8", newline="") as gff, \
                open(outdir / f"{g}.pav.bed", "w", encoding="utf-8", newline="") as bed:
            gff.write("##gff-version 3\n")
            gff.writelines(f"##sequence-region {escape_seqid(sid)} 1 {seq_lengths[g][sid]}\n" for sid in used)
            for _, _, _, loc in rows:
                st = loc.states[g]
                strand = st.strand or "."
                attrs = [f"ID={escape_value(loc.lid)}", f"state={st.state}", f"tier={escape_value(st.tier)}"]
                if loc.family:
                    attrs.append(f"family={escape_value(loc.family)}")
                attrs.append("states=" + ",".join(f"{escape_value(h)}:{loc.states[h].state}" for h in labels))
                gff.write("\t".join([escape_seqid(st.seqid), GFF3_SOURCE, _TYPE[st.state], str(st.start),
                                     str(st.end), ".", strand, ".", ";".join(attrs)]) + "\n")
                bed.write("\t".join([st.seqid, str(st.start - 1), str(st.end), loc.lid, "0", strand, st.state,
                                     _bed(st.tier), _bed(loc.family) if loc.family else "."]) + "\n")
