#!/usr/bin/env python3
"""ltr_tsv_to_gff3.py

Turn the annotated depth-bucketed LTR-RT tables into GFF3.

Writes:
  {prefix}_all_depth_LTR_cleaned.gff3
      Every element from every depth, pooled into one coordinate-sorted GFF3.
  {prefix}_all_depth_protein_LTR_cleaned.gff3
      The same, with the round-1 miniprot protein alignments merged in.
      Only written when a genic GFF exists, i.e. when the run had --proteins.

Per element:

  LTR_retrotransposon                       the element, carrying all attributes
    long_terminal_repeat  (.lLTR, .rLTR)    genomic left/right, not 5'/3' --
                                            on a '-' element the left LTR is the
                                            biological 3' one
    protein_match         (.GENE.n)         one per TEsorter domain

Nested elements stay flat top-level features carrying depth= and nest_status=;
GFF3 Parent means part-of, and a nested LTR-RT is not a component of its host.

Coordinates for the two LTRs come from ltr5_end/ltr3_start, never LTR_len: those
two are mutually consistent on every row, while LTR_len disagrees with both on a
minority of elements.

Usage:
  python ltr_tsv_to_gff3.py --prefix Athal_tair10_chr2_LTRs [--indir .]
                            [--genome genome.fa[.gz]] [--miniprot-gff FILE] [-v]

Input tables must already carry the strand and family columns added by
ltr_annotate.py; missing columns degrade to '.' rather than aborting.
"""
from __future__ import annotations

import argparse
import glob
import gzip
import math
import os
import sys
from typing import Dict, List, NamedTuple, Optional, Sequence, Tuple

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from ltr_annotate import (  # noqa: E402
    UNKNOWN,
    ELEMENT_RE,
    collect_elements,
    element_key,
    header_names,
    load_families,
    load_pass2_links,
    load_tesorter_strands,
    load_unannotated,
    lookup_family,
    parse_domains_field,
    read_table,
    resolve_strands,
    select_annotation_set,
)

SOURCE = "synLTR"

# Attribute values needing escaping in GFF3 column 9. '%' is in the table so it
# is escaped too, which keeps the encoding reversible.
_ESCAPES = {ch: "%%%02X" % ord(ch) for ch in "\t\n\r;=%&,"}

LTR_TYPE = "LTR_retrotransposon"
LTR_REPEAT_TYPE = "long_terminal_repeat"
DOMAIN_TYPE = "protein_match"


def warn(msg: str) -> None:
    print(f"[ltr_tsv_to_gff3] WARNING: {msg}", file=sys.stderr)


def log(msg: str) -> None:
    print(f"[ltr_tsv_to_gff3] {msg}")


class Raw(str):
    """An attribute value that is already column-9 safe and must not be escaped.

    Used for values whose ',' is meant as GFF3's own multi-value separator and
    whose '%' is meant to read as a percent sign. '%' followed by a non-hex pair
    is not a valid percent-escape, so a spec-compliant reader passes it through
    unchanged rather than mis-decoding it.
    """


def escape(value: object) -> str:
    if isinstance(value, Raw):
        return str(value)
    return "".join(_ESCAPES.get(ch, ch) for ch in str(value))


def render_attributes(pairs: Sequence[Tuple[str, object]]) -> str:
    """Join attribute pairs, dropping ones with no information ('.' or empty).

    Omitting absent attributes rather than writing 'tsd=.' is the GFF3 norm and
    keeps lines short; absence is the encoding for 'not known'.
    """
    out = []
    for key, value in pairs:
        text = "" if value is None else str(value)
        if text in ("", UNKNOWN):
            continue
        out.append(f"{escape(key)}={escape(value)}")
    return ";".join(out)


def format_clade_composition(clades: Sequence[Tuple[str, int]], total: int) -> str:
    """Family clade makeup as '50%_SIRE,30%_unknown,15%_CRM,5%_Tekay'.

    Percentages are of family size, rounded half-up to whole numbers; a clade
    that rounds to 0% is dropped rather than written as '0%_'. Because of that
    (and of rounding generally) the printed percentages need not total 100.
    """
    if total <= 0:
        return ""
    parts = []
    for clade, count in clades:
        percent = int(math.floor(100.0 * count / total + 0.5))
        if percent > 0:
            parts.append((percent, clade))
    parts.sort(key=lambda item: (-item[0], item[1]))
    return ",".join(f"{percent}%_{clade}" for percent, clade in parts)


def gff_line(seqid: str, ftype: str, start: int, end: int, strand: str,
             attributes: str, score: str = ".", phase: str = ".") -> str:
    return "\t".join((escape(seqid), SOURCE, ftype, str(start), str(end),
                      score, strand, phase, attributes))


# -----------------------------
# seqid ordering
# -----------------------------
class SeqidRanker:
    """Assigns a sort rank per sequence: genome order when known, else the
    order the sequence is first seen."""

    def __init__(self, known: Sequence[str] = ()):
        self._rank: Dict[str, int] = {name: i for i, name in enumerate(known)}
        self._next = len(self._rank)

    def rank(self, seqid: str) -> int:
        if seqid not in self._rank:
            self._rank[seqid] = self._next
            self._next += 1
        return self._rank[seqid]


def is_gzip(path: str) -> bool:
    """Sniff the magic bytes. The worker's genome is always named '.fa' whether
    it is a gzipped original or a plain FP-masked FASTA, so the extension lies."""
    try:
        with open(path, "rb") as fh:
            return fh.read(2) == b"\x1f\x8b"
    except OSError:
        return False


def read_seq_lengths(path: str) -> "Dict[str, int]":
    """Sequence name -> length, by streaming. Lengths only; never resident."""
    opener = gzip.open if is_gzip(path) else open
    lengths: Dict[str, int] = {}
    name: Optional[str] = None
    total = 0
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = total
                name = line[1:].split()[0].rstrip() if len(line) > 1 else ""
                total = 0
            else:
                total += len(line.strip())
    if name is not None:
        lengths[name] = total
    return lengths


# -----------------------------
# Element blocks
# -----------------------------
class Block(NamedTuple):
    """A parent feature and its children, kept together and sorted as a unit.

    `payload` is either the rendered GFF3 lines (element blocks, built here) or
    the byte offsets of the lines in the source file (miniprot blocks, which are
    re-read and passed through verbatim at emit time).
    """
    seqid: str
    rank: int
    start: int
    end: int
    order: int
    payload: Tuple


def _field(row: Sequence[str], names: Sequence[str], column: str,
           default: str = UNKNOWN) -> str:
    """Value of a named column, or `default` when the column or field is absent."""
    try:
        index = names.index(column)
    except ValueError:
        return default
    if index >= len(row):
        return default
    return row[index] or default


def _as_int(text: str) -> Optional[int]:
    try:
        return int(float(text))
    except (TypeError, ValueError):
        return None


def _line_sort_key(line: str) -> Tuple[int, int]:
    """(start, -end) of a rendered GFF3 line, for ordering within a block."""
    cols = line.split("\t")
    start = _as_int(cols[3]) if len(cols) > 4 else None
    end = _as_int(cols[4]) if len(cols) > 4 else None
    return (start or 0, -(end or 0))


def order_block_lines(lines: Sequence[str]) -> List[str]:
    """Parent line first, then the rest by position.

    Neither source is positionally ordered on its own: the TSV `domains` column
    lists domains in discovery order, and miniprot lists CDS in transcript order,
    which runs high-to-low for a '-' strand model. Sorting here keeps each block
    internally consistent with coordinate order.
    """
    if len(lines) < 2:
        return list(lines)
    return [lines[0]] + sorted(lines[1:], key=_line_sort_key)


def build_element_blocks(prefix: str, tables, ranker: SeqidRanker,
                         provenance: Dict[str, Tuple[str, str]],
                         families, verbose: bool = False
                         ) -> Tuple[List[Block], int]:
    """One Block per element across all depth tables. Returns (blocks, skipped).

    Tables are read as written by ltr_annotate.py, i.e. with the strand and
    family columns present; a table lacking them degrades to '.' per column
    rather than failing.
    """
    family_by_name, family_by_coord = families
    blocks: List[Block] = []
    skipped = 0
    serial = 0
    width = 5

    for table in tables:
        header, rows = read_table(table.path)
        names = header_names(header)

        for row in rows:
            name = row[0] if row else ""
            key = element_key(name)
            match = ELEMENT_RE.match(key) if key else None
            if match is None:
                skipped += 1
                continue
            seqid, start, end = match.group(1), int(match.group(2)), int(match.group(3))
            if end < start:
                skipped += 1
                continue

            serial += 1
            elem_id = f"{prefix}_LTRRT_{serial:0{width}d}"
            strand = _field(row, names, "strand")
            if strand not in ("+", "-"):
                strand = "."

            cascade_strand, tier = provenance.get(key, (UNKNOWN, "none"))
            if strand == ".":
                source_label = "none"
            elif strand == cascade_strand:
                source_label = tier
            else:
                source_label = "table"

            family = lookup_family(name, key, family_by_name, family_by_coord)
            # The table's own label is authoritative -- it is what the user
            # reads -- but fall back to the cluster lookup so a table that was
            # never annotated still gets a consistent set of family attributes.
            family_label = _field(row, names, "family")
            if family_label == UNKNOWN and family is not None:
                family_label = family.family_id

            classification = name.split("#", 1)[1] if "#" in name else ""
            bits = classification.split("/") if classification else []
            order_name = bits[0] if len(bits) > 0 else ""
            superfam = bits[1] if len(bits) > 1 else ""
            clade = bits[2] if len(bits) > 2 else ""

            attributes = render_attributes([
                ("ID", elem_id),
                ("Name", key),
                ("classification", classification),
                ("order", order_name),
                ("superfamily", superfam),
                ("clade", clade),
                ("family", family_label),
                ("family_rep", family.representative if family else ""),
                ("family_size", family.size if family else ""),
                ("family_clades",
                 Raw(format_clade_composition(family.clades, family.size))
                 if family else ""),
                ("depth", table.depth),
                ("ltr_len", _field(row, names, "LTR_len")),
                ("aln_len", _field(row, names, "aln_len")),
                ("subs", _field(row, names, "subs")),
                ("ti", _field(row, names, "ti")),
                ("tv", _field(row, names, "tv")),
                ("K2P_d", _field(row, names, "K2P_d")),
                ("K2P_T", _field(row, names, "K2P_T")),
                ("tsd", _field(row, names, "tsd")),
                ("strand_source", source_label),
                ("nest_status", _field(row, names, "nest_status")),
            ])
            lines = [gff_line(seqid, LTR_TYPE, start, end, strand, attributes)]

            # Two LTRs, from ltr5_end / ltr3_start (element-relative, 1-based).
            ltr5_end = _as_int(_field(row, names, "ltr5_end"))
            ltr3_start = _as_int(_field(row, names, "ltr3_start"))
            if ltr5_end and ltr5_end > 0:
                left_end = min(end, start + ltr5_end - 1)
                if left_end >= start:
                    lines.append(gff_line(
                        seqid, LTR_REPEAT_TYPE, start, left_end, strand,
                        render_attributes([("ID", f"{elem_id}.lLTR"),
                                           ("Parent", elem_id)])))
            if ltr3_start and ltr3_start > 0:
                right_start = max(start, start + ltr3_start - 1)
                if right_start <= end:
                    lines.append(gff_line(
                        seqid, LTR_REPEAT_TYPE, right_start, end, strand,
                        render_attributes([("ID", f"{elem_id}.rLTR"),
                                           ("Parent", elem_id)])))

            # Protein domains, already in genomic coordinates.
            seen_genes: Dict[str, int] = {}
            for gene, dom_start, dom_end in parse_domains_field(
                    _field(row, names, "domains", "")):
                if dom_end < dom_start:
                    dom_start, dom_end = dom_end, dom_start
                seen_genes[gene] = seen_genes.get(gene, 0) + 1
                lines.append(gff_line(
                    seqid, DOMAIN_TYPE, dom_start, dom_end, strand,
                    render_attributes([
                        ("ID", f"{elem_id}.{gene}.{seen_genes[gene]}"),
                        ("Parent", elem_id),
                        ("Name", gene),
                        ("clade", clade),
                    ])))

            blocks.append(Block(seqid, ranker.rank(seqid), start, end, serial,
                                tuple(order_block_lines(lines))))

    if verbose:
        log(f"elements: {len(blocks)} block(s) from {len(tables)} table(s)")
    return blocks, skipped


# -----------------------------
# miniprot passthrough
# -----------------------------
def find_miniprot_gff(prefix: str, indir: str) -> Optional[str]:
    """The round-1 genic GFF, i.e. miniprot's protein alignments.

    Round 1 specifically: later rounds run on LTR-masked genomes, so their
    alignments are incomplete over exactly the regions of interest. Same glob
    ltrharvest_plots.sh uses for TEGV's gene track.
    """
    hits = sorted(glob.glob(os.path.join(indir, f"{prefix}_r1.work", "*.genic.gff")))
    hits = [p for p in hits if os.path.getsize(p) > 0]
    if len(hits) == 1:
        return hits[0]
    if len(hits) > 1:
        warn(f"{len(hits)} round-1 genic GFFs found; skipping the protein GFF3")
    return None


def _attr_value(attributes: str, key: str) -> Optional[str]:
    for field in attributes.split(";"):
        field = field.strip()
        if field.startswith(key + "="):
            return field[len(key) + 1:]
    return None


def index_miniprot_blocks(path: str, ranker: SeqidRanker,
                          verbose: bool = False) -> List[Block]:
    """Index the genic GFF by block without holding its lines in memory.

    A line with ID= and no Parent= opens a block; Parent= lines join the block
    their parent opened. Only (rank, span, byte offset, line count) is retained,
    so peak memory scales with the number of gene models rather than the size of
    the file -- these can run to millions of lines on a large genome. Lines are
    re-read from disk at emit time and passed through verbatim, preserving
    miniprot's own ID/Parent links and its unescaped Target= values.
    """
    index: List[Block] = []
    open_blocks: Dict[str, int] = {}
    # Parallel to `index`: byte offsets of each line belonging to that block.
    offsets: List[List[int]] = []

    with open(path, "rb") as fh:
        position = 0
        for raw in fh:
            start_offset = position
            position += len(raw)
            if raw.startswith(b"#"):
                continue
            line = raw.decode("utf-8", "replace").rstrip("\n")
            cols = line.split("\t")
            if len(cols) < 9:
                continue
            seqid = cols[0]
            begin = _as_int(cols[3])
            finish = _as_int(cols[4])
            if begin is None or finish is None:
                continue

            parent = _attr_value(cols[8], "Parent")
            slot = open_blocks.get(parent) if parent else None
            if slot is None:
                slot = len(index)
                index.append(Block(seqid, ranker.rank(seqid), begin, finish,
                                   slot, ()))
                offsets.append([start_offset])
                own_id = _attr_value(cols[8], "ID")
                if own_id:
                    open_blocks[own_id] = slot
            else:
                offsets[slot].append(start_offset)
                block = index[slot]
                index[slot] = block._replace(start=min(block.start, begin),
                                             end=max(block.end, finish))

    packed = [block._replace(payload=tuple(offsets[i]))
              for i, block in enumerate(index)]
    if verbose:
        log(f"miniprot: {len(packed)} block(s) indexed from "
            f"{os.path.basename(path)}")
    return packed


def read_lines_at(handle, offsets: Sequence[int]) -> List[str]:
    """Re-read a miniprot block's lines verbatim, in coordinate order.

    Verbatim matters: miniprot's own ID/Parent links and its unescaped
    `Target=NP_x 1 302` values (spaces are legal in GFF3 attribute values) must
    survive the round trip untouched.
    """
    out = []
    for offset in offsets:
        handle.seek(offset)
        out.append(handle.readline().decode("utf-8", "replace").rstrip("\n"))
    return order_block_lines(out)


# -----------------------------
# Writing
# -----------------------------
def write_gff3(path: str,
               element_blocks: Sequence[Block],
               seq_lengths: Optional[Dict[str, int]],
               ranker: SeqidRanker,
               miniprot_path: Optional[str] = None,
               miniprot_blocks: Sequence[Block] = (),
               verbose: bool = False) -> int:
    """Emit one GFF3. Returns the number of feature lines written."""
    tagged = [(b, False) for b in element_blocks] + \
             [(b, True) for b in miniprot_blocks]
    tagged.sort(key=lambda item: (item[0].rank, item[0].start, -item[0].end,
                                  item[1], item[0].order))

    seqids_used = {block.seqid for block, _ in tagged}

    n_lines = 0
    handle = open(miniprot_path, "rb") if miniprot_path else None
    try:
        with open(path, "w") as out:
            out.write("##gff-version 3\n")
            if seq_lengths:
                for name, length in seq_lengths.items():
                    if name in seqids_used:
                        out.write(f"##sequence-region\t{name}\t1\t{length}\n")
            for block, is_miniprot in tagged:
                lines = (read_lines_at(handle, block.payload) if is_miniprot
                         else list(block.payload))
                for line in lines:
                    out.write(line + "\n")
                n_lines += len(lines)
                out.write("###\n")
    finally:
        if handle is not None:
            handle.close()

    if verbose:
        log(f"wrote {os.path.basename(path)}: {len(tagged)} block(s), "
            f"{n_lines} feature line(s)")
    return n_lines


# -----------------------------
# Strand provenance
# -----------------------------
def strand_provenance(prefix: str, indir: str, tables, families,
                      verbose: bool = False) -> Dict[str, Tuple[str, str]]:
    """Recompute the strand cascade purely to label each call's origin.

    Returns element key -> (strand the cascade would give, tier name). The
    strand *value* written to the GFF3 always comes from the table, which is
    what the user sees in the TSV; this only supplies the strand_source
    attribute. A disagreement -- a hand-edited table, say -- is reported as
    'table' rather than a tier name.

    """
    elements = collect_elements(load_unannotated(t.path) for t in tables)

    # warn_missing is off: ltr_annotate already reported any missing input when
    # it wrote these tables, and this pass reads exactly the same files.
    tesorter = load_tesorter_strands(prefix, indir, verbose=False,
                                     warn_missing=False)
    target_of, orientation = load_pass2_links(prefix, indir, verbose=False,
                                              warn_missing=False)
    strand, source = resolve_strands(elements, tesorter, target_of, orientation,
                                     verbose=False)
    return {key: (strand[key], source[key]) for key in elements
            if key in strand and key in source}


# -----------------------------
# Driver
# -----------------------------
def convert(prefix: str, indir: str = ".", genome: Optional[str] = None,
            miniprot_gff: Optional[str] = None, verbose: bool = False) -> int:
    variant, tables = select_annotation_set(prefix, indir)
    if not tables:
        print(f"[ltr_tsv_to_gff3] ERROR: no {prefix}_depth<N>[_clean]_ltr.tsv "
              f"found in {indir}", file=sys.stderr)
        return 1
    if variant == "raw":
        warn("no *_clean_ltr.tsv found; building the GFF3 from the raw depth "
             "tables (false-positive families not purged)")
    if verbose:
        log(f"annotation set ({variant}): "
            + ", ".join(os.path.basename(t.path) for t in tables))

    seq_lengths: Optional[Dict[str, int]] = None
    if genome:
        if os.path.isfile(genome):
            seq_lengths = read_seq_lengths(genome)
            if verbose:
                log(f"genome: {len(seq_lengths)} sequence(s) measured")
        else:
            warn(f"genome not found: {genome}; omitting ##sequence-region")

    ranker = SeqidRanker(list(seq_lengths) if seq_lengths else ())
    families = load_families(prefix, indir, verbose=False, warn_missing=False)
    provenance = strand_provenance(prefix, indir, tables, families, verbose)
    element_blocks, skipped = build_element_blocks(prefix, tables, ranker,
                                                   provenance, families, verbose)
    if skipped:
        warn(f"{skipped} row(s) had an unparseable element name and were skipped")
    if not element_blocks:
        # Every element purged as a false positive, or every name unparseable.
        # Detection itself succeeded, so emit a valid empty GFF3 rather than
        # failing the run at the last step.
        warn("no usable elements in the depth tables; writing an empty GFF3")

    ltr_out = os.path.join(indir, f"{prefix}_all_depth_LTR_cleaned.gff3")
    write_gff3(ltr_out, element_blocks, seq_lengths, ranker, verbose=verbose)
    log(f"LTR-RT GFF3: {os.path.basename(ltr_out)} ({len(element_blocks)} elements)")

    genic = miniprot_gff or find_miniprot_gff(prefix, indir)
    if genic and not os.path.isfile(genic):
        warn(f"miniprot GFF not found: {genic}")
        genic = None
    if not genic:
        if verbose:
            log("no round-1 genic GFF (run with --proteins to produce one); "
                "skipping the protein GFF3")
        return 0

    miniprot_blocks = index_miniprot_blocks(genic, ranker, verbose)
    protein_out = os.path.join(indir,
                               f"{prefix}_all_depth_protein_LTR_cleaned.gff3")
    write_gff3(protein_out, element_blocks, seq_lengths, ranker,
               miniprot_path=genic, miniprot_blocks=miniprot_blocks,
               verbose=verbose)
    log(f"LTR-RT + miniprot GFF3: {os.path.basename(protein_out)} "
        f"({len(element_blocks)} elements, {len(miniprot_blocks)} gene models)")
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Write GFF3 annotation from the depth-bucketed LTR-RT "
                    "tables produced by ltrharvest_wrapper2.sh.")
    parser.add_argument("--prefix", required=True,
                        help="Wrapper output prefix, e.g. Athal_tair10_chr2_LTRs")
    parser.add_argument("--indir", default=".",
                        help="Directory holding the wrapper outputs (default: .)")
    parser.add_argument("--genome", default=None,
                        help="Genome FASTA (plain or gzipped). Streamed once for "
                             "##sequence-region lines; omitted if not given.")
    parser.add_argument("--miniprot-gff", default=None,
                        help="miniprot genic GFF to merge (default: autodetect "
                             "<prefix>_r1.work/*.genic.gff)")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Per-step progress and per-file counts")
    args = parser.parse_args(argv)

    if not os.path.isdir(args.indir):
        print(f"[ltr_tsv_to_gff3] ERROR: not a directory: {args.indir}",
              file=sys.stderr)
        return 1
    return convert(args.prefix, args.indir, args.genome, args.miniprot_gff,
                   args.verbose)


if __name__ == "__main__":
    sys.exit(main())
