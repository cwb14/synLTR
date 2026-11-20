"""VCFs, bgzipped and tabix-indexed, of two record types. The two share their contig lines, the order of contigs and
records, the escaping of values, the unplaced table and the bgzip and tabix code; the rest of each header is the
record type's own.

Catalog VCF (`write_vcf`, LTRpav): the locus catalog with haploid genotypes, on the reference genome's coordinates.
An insertion is anchored VCF-style on the base before it: a locus present in the reference genome (P, S, T) has
POS = start - 1; one empty there (E) has POS = the last base of the target-site copy, or, for an E site without a
readable TSD copy (its site is then the flank junction), the junction's last base: the insertion point to within
1 bp (the header's LTRpav_POS line says so). REF is the reference base at POS, ALT <INS:ME:LTR>, END = POS, SVLEN
the locus length. GT: 1 present, 0 absent (E), missing for D, O, G, U. Loci the VCF cannot hold go to a side file
(lid, ref_state, reason).

Reference VCF (`write_reference_vcf`, ltrkin): the reads route's calls of the reference's elements; GT 0 is the
reference's element (full), the ALT alleles are <DEL:ME:LTR> (the empty site) and <DEL:ME:SOLO> (a solo LTR).

Both: IDs and INFO strings are percent-encoded (`escape`); contig lines and records follow the reference FASTA's .fai
order; a record is anchored on the base before its event, so an event starting at base 1 (no such base) is not
written: it goes to the VCF's unplaced table (<stem>.unplaced.tsv: lid, ref_state, reason), reason starts_at_base_1
(`anchored`; a POS 0 record would be dropped by htslib's indexed readers).
"""
from __future__ import annotations

import os
import re
import statistics
from collections.abc import Sequence
from pathlib import Path
from typing import TYPE_CHECKING

import pysam

from ltrtrace.base.constants import PROGRAM
from ltrtrace.io.samples import Sample
from ltrtrace.io.seqio import Fasta

from .loci import Locus
from .states import OBSERVED_ABSENT, call_class

if TYPE_CHECKING:
    from ltrtrace.routes.reads.genotype import SampleCalls
    from ltrtrace.routes.reads.reads import Target

# the header's description of POS (##LTRpav_POS): where each site is anchored, the E site's too
POS_NOTE = ("POS is the base before the insertion: start - 1 where the reference genome holds the element (ST P, "
            "S or T); at its empty site (ST E), the last base of the target-site (TSD) copy, or, without a "
            "readable TSD copy, the flank junction, to within 1 bp (within the few extra bases, at most two TSD "
            "lengths, the reference may hold there)")
_ESC = re.compile(r"[%;=,\s]")  # characters with a meaning in VCF ID and INFO values (VCF 4.3, 1.2), and whitespace

_GT = {"present": (1,), "absent": (0,), "missing": (None,)}
UNPLACED_COLUMNS = ("lid", "ref_state", "reason")  # the unplaced table of either VCF
BASE_1 = "starts_at_base_1"  # unplaced reason: the event starts at base 1, so no base before it can pad its record
REF_STATE = "P"  # the reference's state at a reads-route element: it holds the element, an intact copy

# the reference VCF's ALT, INFO and FORMAT lines
REF_ALTS = "<DEL:ME:LTR>,<DEL:ME:SOLO>"
_REF_META = [
    '##ALT=<ID=DEL:ME:LTR,Description="Empty site: the element and one TSD copy absent (pre-insertion allele)">',
    '##ALT=<ID=DEL:ME:SOLO,Description="Solo LTR: internal region and one LTR absent">',
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Structural variant type">',
    ('##INFO=<ID=END,Number=1,Type=Integer,Description="Last base of the reference sequence the empty site lacks: the '
     'element and one TSD copy (POS plus their length)">'),
    '##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="Length change of each ALT allele">',
    '##INFO=<ID=ELEMENT,Number=1,Type=String,Description="Element ID from the annotation">',
    '##INFO=<ID=FAMILY,Number=1,Type=String,Description="Family">',
    '##INFO=<ID=CLADE,Number=1,Type=String,Description="Clade">',
    ('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype: 0 full, 1 empty site, 2 solo LTR, . null (haplotype '
     'without the site)">'),
    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality (phred)">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Informative reads per allele (full, empty, solo)">',
    '##FORMAT=<ID=FT,Number=1,Type=String,Description="Per-sample flags">',
]
_REF_COLUMNS = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]


def escape(value) -> str:
    """`value` as text with these characters percent-encoded: % ; = , (they have a meaning in VCF 4.3, 1.2) and
    every whitespace character (Python's \\s): tab, newline, carriage return, form feed, vertical tab, space,
    \\x1c-\\x1f, \\x85, \\xa0 and the other Unicode spaces and line separators. Each becomes % and its code point in
    upper-case hex: two digits up to U+00FF (%3B, %0C, %A0), four above it (%2028, %3000)."""
    return _ESC.sub(lambda m: f"%{ord(m.group()):02X}", str(value))


def contig_order(fasta: Fasta) -> dict[str, int]:
    """Each sequence of the reference FASTA by its index in the .fai: the order of contigs and records."""
    return {sid: i for i, sid in enumerate(fasta.lengths())}


def _contigs(fasta: Fasta, held: set[str] | None = None) -> list[tuple[str, int]]:
    """(name, length) of the reference's sequences in .fai order: all of them, or those in `held`."""
    return [(sid, n) for sid, n in fasta.lengths().items() if held is None or sid in held]


def _bgzip_index(plain: Path) -> Path:
    """Compress the text VCF `plain` with bgzip into <plain>.gz, index that (.tbi), remove `plain`; returns the
    .gz path."""
    return Path(pysam.tabix_index(str(plain), preset="vcf", force=True))


def unplaced_path(vcf: Path) -> Path:
    """The unplaced table beside the VCF `vcf`: <stem>.unplaced.tsv (pav.vcf.gz -> pav.unplaced.tsv)."""
    vcf = Path(vcf)
    return vcf.with_name(vcf.name.removesuffix(".vcf.gz").removesuffix(".gz") + ".unplaced.tsv")


def anchored(pos: int) -> bool:
    """True if `pos`, the 1-based position of the base a record is anchored on (the base before its event), is on
    the sequence. An event starting at base 1 has no base before it: both record types leave it out of the VCF and
    list it in the unplaced table with reason BASE_1 (VCF 4.3's POS 0, a telomere, is dropped without a message by
    htslib's indexed readers: bcftools merge and view -r)."""
    return pos >= 1


def _write_unplaced(path: Path, rows: list[tuple[str, str, str]]) -> None:
    """The unplaced table at `path`: a header of UNPLACED_COLUMNS, then (lid, ref_state, reason) per row."""
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(UNPLACED_COLUMNS) + "\n")
        fh.writelines(f"{lid}\t{state}\t{reason}\n" for lid, state, reason in rows)


def _header(samples: list[Sample], fasta: Fasta) -> pysam.VariantHeader:
    h = pysam.VariantHeader()
    h.add_meta("fileformat", "VCFv4.2")
    h.add_meta("source", PROGRAM)
    h.add_meta("LTRpav_POS", POS_NOTE)
    for seqid, n in _contigs(fasta):
        h.contigs.add(seqid, length=n)
    h.add_meta("ALT", items=[("ID", "INS:ME:LTR"), ("Description", "LTR retrotransposon insertion")])
    for name, num, typ, desc in (
            ("SVTYPE", 1, "String", "Structural variant type"),
            ("SVLEN", 1, "Integer", "Length of the LTR retrotransposon"),
            ("END", 1, "Integer", "End position (equal to POS)"),
            ("FAM", 1, "String", "Family"),
            ("SF", 1, "String", "Superfamily"),
            ("K", 1, "Float", "Median K2P divergence of the member elements' LTRs"),
            ("DISC", ".", "String", "Genomes whose annotation holds the insertion"),
            ("NE", 1, "Integer", "Genomes with the empty site (state E)"),
            ("AA", 1, "String", "Ancestral allele: the reference base, given when an empty site is observed")):
        h.info.add(name, num, typ, desc)
    h.formats.add("GT", 1, "String", "Genotype: 1 present, 0 absent, . missing")
    h.formats.add("ST", 1, "String", "Call state (P S T E D O G U)")
    h.formats.add("TR", 1, "String", "Call tier")
    for s in samples:
        h.add_sample(s.label)
    return h


def write_vcf(loci: list[Locus], samples: list[Sample], ref_label: str, fasta: Fasta, path: Path, *,
              absent: frozenset[str] | set[str] = OBSERVED_ABSENT) -> tuple[Path, Path]:
    """Write `path` (.vcf.gz) and its .tbi, from the loci in the coordinates of genome `ref_label` (`fasta`).
    Returns (vcf.gz path, unplaced TSV path: <stem>.unplaced.tsv next to it). A genotype is 1 for a present state,
    0 for a state in `absent` (call_class; E alone by default), missing otherwise; a reference state of the
    missing class is unplaced."""
    path = Path(path)
    if path.suffix != ".gz":
        raise ValueError(f"VCF path must end in .gz: {path}")
    plain = path.with_suffix("")
    unplaced = unplaced_path(path)
    order = contig_order(fasta)
    labels = [s.label for s in samples]
    header = _header(samples, fasta)

    keep, skipped = [], []
    for loc in loci:
        st = loc.states[ref_label]
        klass = call_class(st.state, absent=absent)
        if klass == "missing":
            reason = "ref_state_not_placeable"
        elif st.seqid is None or st.start is None or st.end is None:
            reason = "no_coordinates"
        elif st.seqid not in order:
            reason = "seqid_not_in_reference"
        else:
            pos = st.start - 1 if klass == "present" else st.end
            if not anchored(pos):
                skipped.append((loc.lid, st.state, BASE_1))
            else:
                keep.append((order[st.seqid], pos, loc.lid, st.seqid, loc))
            continue
        skipped.append((loc.lid, st.state, reason))
    keep.sort(key=lambda t: t[:3])

    with pysam.VariantFile(str(plain), "w", header=header) as out:
        for _, pos, _, seqid, loc in keep:
            ref = fasta.fetch(seqid, pos - 1, pos)
            ne = sum(loc.states[g].state in absent for g in labels)
            ks = [loc.states[g].k for g in loc.members if g in loc.states and loc.states[g].k is not None]
            rec = header.new_record(contig=seqid, start=pos - 1, stop=pos, alleles=(ref, "<INS:ME:LTR>"),
                                    id=loc.lid)
            rec.info["SVTYPE"] = "INS"
            rec.info["SVLEN"] = int(loc.length)
            if loc.family is not None:
                rec.info["FAM"] = escape(loc.family)
            if loc.superfamily is not None:
                rec.info["SF"] = escape(loc.superfamily)
            if ks:
                rec.info["K"] = float(statistics.median(ks))
            if loc.discovery:
                rec.info["DISC"] = tuple(escape(g) for g in loc.discovery)
            rec.info["NE"] = ne
            if ne:
                rec.info["AA"] = ref
            for g in labels:
                st = loc.states[g]
                call = rec.samples[g]
                call["GT"] = _GT[call_class(st.state, absent=absent)]
                call.phased = False
                call["ST"] = st.state
                call["TR"] = st.tier
            out.write(rec)
    _bgzip_index(plain)
    _write_unplaced(unplaced, skipped)
    return path, unplaced


def write_reference_vcf(path: Path, targets: Sequence[Target], samples: Sequence[SampleCalls], fasta: Fasta,
                        reference: Path) -> tuple[Path, Path]:
    """Write the reads route's VCF `path` (.vcf.gz, with its .tbi): one record per element of `targets` (routes.reads
    Target, 0-based half-open) that at least one of `samples` holds, one column per sample. Returns (vcf.gz path,
    unplaced table path: <stem>.unplaced.tsv next to it).

    POS is the base before the element: its 0-based start read as a 1-based position. An element starting at base 1
    has none: it is left out and listed in the unplaced table (its ID, ref_state P, reason starts_at_base_1). END =
    POS + the element's length + its TSD's length, the reference sequence the empty site lacks, so END - POS =
    -SVLEN of <DEL:ME:LTR>; SVLEN of <DEL:ME:SOLO> is minus the element's length less its mean LTR length. Contig
    lines (sequences holding an element) and records follow the .fai order of `fasta`, equal POS keeping the order of
    `targets`. A sample without data at an element gets a no-call of its ploidy flagged no_data. `reference` is the
    FASTA's path for the header."""
    path = Path(path)
    if path.suffix != ".gz":
        raise ValueError(f"VCF path must end in .gz: {path}")
    plain = path.with_suffix("")
    order = contig_order(fasta)
    hdr = ["##fileformat=VCFv4.3", f"##source={PROGRAM}", f"##reference={os.path.abspath(reference)}"]
    hdr += [f"##contig=<ID={c},length={n}>" for c, n in _contigs(fasta, {t.chrom for t in targets})]
    hdr += _REF_META
    hdr.append("#" + "\t".join(_REF_COLUMNS + [s.label for s in samples]))
    held = [i for i in range(len(targets)) if any(i in s.calls for s in samples)]
    held.sort(key=lambda i: (order[targets[i].chrom], targets[i].start))
    unplaced = []
    with open(plain, "w") as fh:
        fh.write("\n".join(hdr) + "\n")
        for i in held:
            el = targets[i]
            pos = el.start
            if not anchored(pos):
                unplaced.append((el.id, REF_STATE, BASE_1))
                continue
            ref = fasta.fetch(el.chrom, pos - 1, pos)
            solo = -int(el.length - el.ltr_mean)
            info = (f"SVTYPE=DEL;END={pos + el.length + el.tsd_len};SVLEN={-(el.length + el.tsd_len)},{solo};"
                    f"ELEMENT={escape(el.id)};FAMILY={escape(el.family)};CLADE={escape(el.clade)}")
            cells = []
            for s in samples:
                c = s.calls.get(i)
                if c is None:
                    cells.append("/".join("." * s.ploidy) + ":0:0,0,0:no_data")
                else:
                    cells.append(f"{c.gt}:{c.gq}:{c.ad[0]},{c.ad[1]},{c.ad[2]}:{';'.join(c.flags) or 'PASS'}")
            fh.write("\t".join([el.chrom, str(pos), escape(el.id), ref, REF_ALTS, ".", ".", info, "GT:GQ:AD:FT"]
                               + cells) + "\n")
    _bgzip_index(plain)
    _write_unplaced(unplaced_path(path), unplaced)
    return path, unplaced_path(path)
