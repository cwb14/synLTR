"""The library loader (io.library): record names, the columnar store, the LTR record offsets of both LTRquest clean
conventions read from a table and from GFF3, the join of the records to their annotation, the screen's rule (LTRmix's
reading) and the library specs of a samples sheet."""
import random
import shutil
from pathlib import Path

import numpy as np
import pysam
import pytest
from builders import ElementSpec, build_scene, rand_seq, write_fasta, write_gff3, write_ltrquest

from synltr.io.annotation import load_annotation, load_annotations
from synltr.io.elements import NO_LTR, Element, ElementTable, make_eid
from synltr.io.library import (
    LibraryError,
    RecordName,
    library_files,
    load_library,
    parse_record_name,
    read_library,
    record_ends,
    record_offsets,
)
from synltr.io.seqio import encode, revcomp

# phasing's annotation readers (analyses.subgenomes.inputs reads with these)
PHASING = {"ltr_rule": "optional", "host_rule": "k2p", "table_rule": "phasing", "tsd_rule": "acgtn",
           "k_keys": ("K2P_d", "k2p", "K2P")}


@pytest.mark.parametrize("name,want", [
    ("chr1:10-20#LTR/Gypsy", RecordName("chr1", 10, 20, "LTR/Gypsy")),
    ("HG002#1#chr1:10-20#LTR/Copia", RecordName("HG002#1#chr1", 10, 20, "LTR/Copia")),   # PanSN seqid
    ("chr1:10..20(-)", RecordName("chr1", 10, 20, None)),
    ("chr1:10-20(.)#LTR", RecordName("chr1", 10, 20, "LTR")),
    ("elem::chr1:9-20", RecordName("chr1", 10, 20, None)),                                 # bedtools: 0-based start
    ("x::HG002#1#chr1:9-20(+)#LTR/Gypsy", RecordName("HG002#1#chr1", 10, 20, "LTR/Gypsy")),
    ("chr1:10-20#", RecordName("chr1", 10, 20, "")),
])
def test_parse_record_name(name, want):
    assert parse_record_name(name) == want


def test_parse_record_name_variants():
    """The record names of LTRquest, EDTA and bedtools getfasta -name, with an optional strand suffix."""
    assert parse_record_name("chr1:10-500#LTR/Copia/Ale") == ("chr1", 10, 500, "LTR/Copia/Ale")
    assert parse_record_name("chr1:10..500") == ("chr1", 10, 500, None)
    assert parse_record_name("elt7::Chr_2:5-90(+)") == ("Chr_2", 6, 90, None)
    assert parse_record_name("scaf:1:2:10-20") == ("scaf:1:2", 10, 20, None)
    assert parse_record_name("nonsense") is None
    assert parse_record_name("chr1:10-x") is None


def test_record_offsets_strands():
    """Genome LTR intervals (1-based closed) of the element at 100..1099 -> 0-based half-open record offsets, 5' LTR
    first: a minus-strand record is the reverse complement of the span; any other record is in genome orientation."""
    assert record_offsets([(100, 149), (1000, 1099)], 100, 1099, "+") == ((0, 50), (900, 1000))
    assert record_offsets([(100, 149), (1000, 1099)], 100, 1099, "-") == ((0, 100), (950, 1000))
    assert record_offsets([(100, 149), (1000, 1099)], 100, 1099, ".") == ((0, 50), (900, 1000))


def test_record_ends():
    """The LTRs of a record shorter than its span, at its ends (0-based half-open, 5' first): the genome-left LTR
    first on a plus or unknown strand, the genome-right LTR first on a minus-strand record (reverse complement)."""
    assert record_ends(12, 10, "+", 60) == ((0, 12), (50, 60))
    assert record_ends(12, 10, ".", 60) == ((0, 12), (50, 60))
    assert record_ends(12, 10, "-", 60) == ((0, 10), (48, 60))
    for left, right, length in ((None, 10, 60), (float("nan"), 10, 60), (0, 10, 60), (12, -1, 60), (40, 30, 60)):
        assert record_ends(left, right, "+", length) is None


def _minus_host(seed):
    """One minus-strand host (genome-left LTR 250 bp, genome-right 320 bp, annotated K2P 0) holding a plus-strand
    child in its internal region: (genome records, elements)."""
    child = ElementSpec(left=150, right=150, internal=400, strand="+", tsd="GATC", k=0.0)
    host = ElementSpec(left=250, right=320, internal=1500, strand="-", tsd="ACGTA", k=0.0, nests=[(600, child)])
    return build_scene(random.Random(seed), [host])


def _library(tmp_path, elements, records, clean, source, rules):
    prefix = str(tmp_path / "lib")
    write_ltrquest(prefix, elements, records, clean=clean)
    spec = prefix if source == "table" else str(write_gff3(tmp_path / "a.gff3", elements, lengths=True))
    table = load_annotation(spec, "G", **(PHASING if rules == "phasing" else {}))
    return load_library(sorted(tmp_path.glob("lib_depth*_clean_ltr.fa")), table)


def _row(lib, e):
    return int(np.flatnonzero(lib.elements.column("eid") == e.eid)[0])


@pytest.mark.parametrize("rules", ["shared", "phasing"])
@pytest.mark.parametrize("source", ["table", "gff3"])
@pytest.mark.parametrize("clean", ["cut", "masked"])
def test_minus_strand_unequal_ltrs(tmp_path, clean, source, rules):
    """A minus-strand record whose two LTRs differ in length (250 and 320 bp), its nest cut out or N-masked, read
    from the LTRquest table (LTR offsets along the genome span) or from GFF3: the record is the reverse complement
    of the span, so l5 covers the genome-right LTR at its start and l3 the genome-left LTR at its end, exactly."""
    records, elements = _minus_host(5)
    lib = _library(tmp_path, elements, records, clean, source, rules)
    host = next(e for e in elements if e.level == 0)
    i = _row(lib, host)
    chrom = records[host.seqid]
    left, right = (chrom[a - 1:b] for a, b in (host.ltr_left, host.ltr_right))
    assert (len(left), len(right)) == (250, 320)
    rec = lib.seq(i)
    (a5, b5), (a3, b3) = lib.l5[i], lib.l3[i]
    assert np.array_equal(rec[a5:b5], encode(revcomp(right).encode()))
    assert np.array_equal(rec[a3:b3], encode(revcomp(left).encode()))
    assert (lib.orient[i], lib.clean[i], lib.bounds[i]) == (-1, clean, "ends" if clean == "cut" else "mapped")
    child = next(e for e in elements if e.level == 1)
    j = _row(lib, child)
    assert (lib.orient[j], lib.clean[j], lib.bounds[j]) == (1, "none", "mapped")


def test_record_left_whole_is_clean_none(tmp_path):
    """A host whose record keeps its nest is clean "none"; its LTRs map from the genome."""
    records, elements = _minus_host(6)
    lib = _library(tmp_path, elements, records, "none", "table", "shared")
    i = _row(lib, next(e for e in elements if e.level == 0))
    assert (lib.clean[i], lib.bounds[i]) == ("none", "mapped")
    assert tuple(lib.l5[i]) == (0, 320) and tuple(lib.l3[i]) == (lib.seq(i).size - 250, lib.seq(i).size)


def test_equal_length_ltrs_record_relative_table(tmp_path):
    """A minus-strand element whose two LTRs have one length, from its table: the strand-aware offsets equal the
    table's own LTR columns read as record offsets (both LTRs sit at the record's ends)."""
    records, elements = build_scene(random.Random(7), [ElementSpec(left=300, right=300, strand="-", k=0.0)])
    prefix = str(tmp_path / "lib")
    write_ltrquest(prefix, elements, records)
    lib = load_library([Path(f"{prefix}_depth0_clean_ltr.fa")], load_annotation(prefix, "G"))
    head, row = (line.lstrip("#").split("\t") for line in
                 Path(f"{prefix}_depth0_clean_ltr.tsv").read_text().splitlines())
    c = {h: int(row[head.index(h)]) for h in ("ltr5_start", "ltr5_end", "ltr3_start", "ltr3_end")}
    assert tuple(lib.l5[0]) == (c["ltr5_start"] - 1, c["ltr5_end"])
    assert tuple(lib.l3[0]) == (c["ltr3_start"] - 1, c["ltr3_end"])
    assert (lib.orient[0], lib.bounds[0]) == (-1, "mapped")


def _el(seqid, start, end, strand="+", k=0.01, ltrs=None):
    left, right = ltrs or (NO_LTR, NO_LTR)
    return Element(genome="G", eid=make_eid("G", seqid, start, end), name=f"{seqid}:{start}-{end}", seqid=seqid,
                   start=start, end=end, strand=strand, ltr_left=left, ltr_right=right, tsd=None, k=k, family=None,
                   superfamily=None, clade=None, confidence=None)


def _flat(tmp_path, n=30, seed=1):
    """n 60-bp records c{i % 3}:{100 i + 1}-{100 i + 60} (name -> sequence) and their elements, annotated K2P 0.01."""
    rng = random.Random(seed)
    recs = {f"c{i % 3}:{100 * i + 1}-{100 * i + 60}": rand_seq(60, rng) for i in range(n)}
    return recs, ElementTable.from_rows(_el(f"c{i % 3}", 100 * i + 1, 100 * i + 60) for i in range(n))


def test_duplicate_coordinates_keep_first(tmp_path, caplog):
    """A record repeating an earlier record's coordinates (here written as start..end) is dropped with one warning,
    the first kept; more than 1 - min_match of such records fails."""
    recs, table = _flat(tmp_path)
    first = recs["c0:1-60"]
    fa = write_fasta(tmp_path / "l.fa", {**recs, "c0:1..60": rand_seq(60, random.Random(9))})
    lib = load_library([fa], table)
    assert lib.n == 30 and np.array_equal(lib.seq(0), encode(first.encode()))
    assert "1 FASTA records repeat the coordinates of an earlier record" in caplog.text
    dups = {f"c{i % 3}:{100 * i + 1}..{100 * i + 60}": "ACGT" for i in range(5)}
    fa = write_fasta(tmp_path / "l.fa", {**recs, **dups})
    with pytest.raises(LibraryError, match="5 of 35 FASTA records repeat the coordinates"):
        load_library([fa], table)


def test_join_below_min_match_fails(tmp_path):
    """Fewer than min_match of the records annotated (by seqid:start-end) fails, naming unmatched records; at a lower
    min_match an unannotated record keeps a row named from its record: no K2P, strand '.', no LTRs."""
    recs, table = _flat(tmp_path, n=3)
    fa = write_fasta(tmp_path / "l.fa", recs)
    some = table.select([0])
    with pytest.raises(LibraryError, match=r"only 1 of 3 FASTA copies are annotated .*c1:101-160"):
        load_library([fa], some)
    lib = load_library([fa], some, min_match=0.3)
    e = lib.elements
    assert list(e.column("eid")) == [make_eid("G", c, s, t) for c, s, t in (("c0", 1, 60), ("c1", 101, 160),
                                                                             ("c2", 201, 260))]
    assert e.column("k")[0] == 0.01 and np.isnan(e.column("k")[1:]).all()
    assert list(e.column("strand")) == ["+", ".", "."] and list(lib.orient) == [1, 0, 0]
    assert e.column("name")[1] == "c1:101-160" and tuple(e.column("ltr_left")[1]) == NO_LTR
    assert list(lib.bounds) == ["none"] * 3 and lib.placed.all()


def test_records_without_coordinates_skipped(tmp_path, caplog):
    """A record whose name holds no seqid:start-end is dropped with one warning; none left fails."""
    recs, table = _flat(tmp_path, n=3)
    fa = write_fasta(tmp_path / "l.fa", {"contig_7": "ACGT", **recs})
    lib = load_library([fa], table)
    assert lib.n == 3 and list(lib.names) == list(recs)
    assert "skipped 1 FASTA records whose name has no seqid:start-end" in caplog.text
    with pytest.raises(LibraryError, match="no usable LTR-RT records"):
        load_library([write_fasta(tmp_path / "x.fa", {"contig_7": "ACGT"})], table)


def test_store_columns(tmp_path):
    """Records of several FASTAs in order: codes (ACGT either case 0-3, anything else 4), record offsets, the class
    of each name, orientation from the element's strand, stripped records without their non-ACGT codes."""
    table = ElementTable.from_rows([_el("c1", 1, 12, "+"), _el("c1", 21, 24, "-"), _el("c2", 1, 4, ".")])
    a = write_fasta(tmp_path / "a.fa", {"c1:1-12#LTR/Gypsy": "ACGTNNacgtRY", "c1:21-24": "GGCC"})
    b = write_fasta(tmp_path / "b.fa", {"c2:1-4#LTR": "TTAA"})
    lib = load_library([a, b], table)
    assert lib.n == 3 and list(lib.offs) == [0, 12, 16, 20]
    assert lib.seq(0).tolist() == [0, 1, 2, 3, 4, 4, 0, 1, 2, 3, 4, 4]
    assert lib.seq(0, stripped=True).tolist() == [0, 1, 2, 3, 0, 1, 2, 3]
    assert list(lib.cls) == ["LTR/Gypsy", None, "LTR"] and list(lib.orient) == [1, -1, 0]
    assert lib.codes.dtype == np.uint8 and lib.offs.dtype == np.int64 and lib.orient.dtype == np.int8


def test_codes_preallocated_from_fai(tmp_path):
    """With a .fai beside the FASTA the codes go into one array of its total length; a stale .fai (fewer or more
    bases than the FASTA holds) changes nothing in the result."""
    recs, table = _flat(tmp_path)
    fa = write_fasta(tmp_path / "l.fa", recs)
    plain = load_library([fa], table)
    pysam.faidx(str(fa))
    indexed = load_library([fa], table)
    assert np.array_equal(indexed.codes, plain.codes) and np.array_equal(indexed.offs, plain.offs)
    fai = Path(f"{fa}.fai")
    lines = fai.read_text().splitlines(keepends=True)
    for stale in (lines[:3], [*lines, "extra\t5000\t0\t60\t61\n"]):
        fai.write_text("".join(stale))
        lib = load_library([fa], table)
        assert np.array_equal(lib.codes, plain.codes) and np.array_equal(lib.offs, plain.offs)


def test_unknown_rule_fails(tmp_path):
    recs, table = _flat(tmp_path, n=3)
    with pytest.raises(ValueError, match="rule"):
        load_library([write_fasta(tmp_path / "l.fa", recs)], table, rule="nonsense")


# ---- the screen's rule: LTRmix's library reading ----
def _screen_table(path, rows):
    """A screen table (LTRmix's columns) of rows (seq_id, k2p, tsd, strand, ltr5_len, ltr3_len)."""
    path.write_text("seq_id\tk2p\ttsd\tstrand\tltr5_len\tltr3_len\n" + "".join("\t".join(map(str, r)) + "\n"
                                                                                for r in rows))
    return path


def test_screen_rule_keeps_unplaced_records(tmp_path, caplog):
    """Under the screen rule a record whose name holds no interval is kept, unplaced: its sequence stays (its k-mers
    count), it has no eid and no coordinates (seqid none, start and end 0), one warning says so; more than half of the
    records unplaced fails, naming the share."""
    recs, _ = _flat(tmp_path, n=3)
    rows = {name: {"seq_id": name, "k2p": "0.01"} for name in recs}
    fa = write_fasta(tmp_path / "l.fa", {"contig_7": "ACGTNacgt", **recs})
    lib = load_library([fa], rows, rule="screen", genome="G")
    e = lib.elements
    assert lib.n == 4 and lib.placed.tolist() == [False, True, True, True] and lib.names[0] == "contig_7"
    assert (e.column("eid")[0], e.column("seqid")[0], e.column("start")[0], e.column("end")[0]) == (None, None, 0, 0)
    assert lib.seq(0).tolist() == [0, 1, 2, 3, 4, 0, 1, 2, 3] and lib.clean[0] == "none" and lib.bounds[0] == "none"
    assert list(e.column("eid")[1:]) == [make_eid("G", f"c{i}", 100 * i + 1, 100 * i + 60) for i in range(3)]
    assert "1 of 4 FASTA records have no seqid:start-end; kept unplaced" in caplog.text
    assert e.column("k").tolist()[1:] == [0.01] * 3 and np.isnan(e.column("k")[0])
    two = write_fasta(tmp_path / "x.fa", {"a": "ACGT", "b": "ACGT", "c0:1-60": recs["c0:1-60"]})
    with pytest.raises(LibraryError, match="2 of 3 FASTA records have no seqid:start-end"):
        load_library([two], rows, rule="screen", genome="G")
    lib = load_library([write_fasta(tmp_path / "y.fa", {"a": "ACGT", "c0:1-60": recs["c0:1-60"]})], rows,
                       rule="screen", genome="G")
    assert lib.placed.tolist() == [False, True]


def test_screen_rule_ltrs_at_the_stripped_record_ends(tmp_path):
    """Under the screen rule every record's LTRs sit at the ends of its N-stripped record, from the LTRs' own lengths
    (ltr5_len for the genome-left LTR; a minus-strand record starts with the genome-right one): offsets in
    seq(i, stripped=True), 0-based half-open, bounds "ends"; an unknown length, or two that do not fit the stripped
    record, leave none. A record the tables lack, or no table at all, is kept (no annotated share is required) under
    the genome label given."""
    fa = write_fasta(tmp_path / "l.fa", {"c1:1-40#LTR": "AAAAANNNNNCCCCCCCCCCGGGGGGGGGGTTTTTTTTTT",
                                         "c1:101-140": "AAAAACCCCCNNNNNGGGGGTTTTTTTTTTACGTACGTAC",
                                         "c1:201-210": "ACGTACGTAC", "c1:301-310": "ACGTACGTAC",
                                         "c1:401-410": "ACGTACGTAC"})
    tb = _screen_table(tmp_path / "l.tsv", [("c1:1-40#LTR", 0.01, "ACGTA", "+", 5, 10),
                                            ("c1:101-140", 0.02, ".", "-", 5, 10),
                                            ("c1:201-210", 0.03, "NA", ".", 6, 5),
                                            ("c1:301-310", "NA", "acgtt", "+", 0, 4)])
    lib = read_library(f"{fa},{tb}", "G")
    assert lib.bounds.tolist() == ["ends", "ends", "none", "none", "none"]
    assert [tuple(map(tuple, x)) for x in zip(lib.l5[:2].tolist(), lib.l3[:2].tolist())] == \
        [((0, 5), (25, 35)), ((0, 10), (30, 35))]
    assert lib.elements.column("eid")[4] == make_eid("G", "c1", 401, 410) and np.isnan(lib.elements.column("k")[4])
    assert lib.placed.all() and list(lib.elements.column("tsd")[:4]) == ["ACGTA", None, None, "ACGTT"]
    bare = read_library(str(fa), "X")
    assert bare.elements.column("eid")[0] == make_eid("X", "c1", 1, 40) and bare.bounds.tolist() == ["none"] * 5


def test_library_files(tmp_path):
    """A samples sheet's library: FASTA,TABLE[,TABLE ...]; a FASTA alone when it ends in .fa, .fasta, .fa.gz,
    .fasta.gz or .fna; else an LTRquest prefix: its <prefix>_depth<N>_clean_ltr.fa in depth order (10 after 9), each
    with the .tsv beside it. A missing file, or a prefix with no library, is a LibraryError naming it."""
    pre = tmp_path / "G_LTRs"
    for d in (10, 0, 9):
        for ext in ("fa", "tsv"):
            (tmp_path / f"G_LTRs_depth{d}_clean_ltr.{ext}").write_text("")
    fas, tables = library_files(str(pre))
    assert [p.name for p in fas] == [f"G_LTRs_depth{d}_clean_ltr.fa" for d in (0, 9, 10)]
    assert tables == [p.with_suffix(".tsv") for p in fas]
    fa, tb = tmp_path / "x.fa", tmp_path / "x.tsv"
    fa.write_text(">c:1-4\nACGT\n")
    tb.write_text("seq_id\tk2p\n")
    assert library_files(f"{fa},{tb}") == ([fa], [tb])
    assert library_files(str(fa)) == ([fa], [])
    with pytest.raises(LibraryError, match="y.tsv"):
        library_files(f"{fa},{tmp_path / 'y.tsv'}")
    with pytest.raises(LibraryError, match="no LTR_peel libraries"):
        library_files(str(tmp_path / "H_LTRs"))
    (tmp_path / "G_LTRs_depth9_clean_ltr.tsv").unlink()
    with pytest.raises(LibraryError, match="G_LTRs_depth9_clean_ltr.tsv"):
        library_files(str(pre))


def test_nesting_across_annotation_files(tmp_path):
    """A host and its nest annotated in two files (LTRquest's depth tables read one by one, as the subgenomes analysis
    reads its annotation list): load_annotations nests them together, so the host's N-masked record is clean
    "masked"; each file nested alone left the host without its nest ("none")."""
    from synltr.analyses.subgenomes.inputs import read_annotation

    records, elements = _minus_host(5)
    prefix = str(tmp_path / "lib")
    write_ltrquest(prefix, elements, records, clean="masked")
    tables = sorted(tmp_path.glob("lib_depth*_clean_ltr.tsv"))
    fastas = sorted(tmp_path.glob("lib_depth*_clean_ltr.fa"))
    assert len(tables) == 2
    host = next(e for e in elements if e.level == 0)
    for table in (load_annotations(tables, "G", **PHASING), read_annotation(tables, "G")):
        child = next(e for e in table if e.level == 1)
        assert child.parent == host.eid
        lib = load_library(fastas, table)
        assert lib.clean[_row(lib, host)] == "masked"
    alone = load_library(fastas, ElementTable.concat([load_annotation(str(t), "G", **PHASING) for t in tables]))
    assert alone.clean[_row(alone, host)] == "none"


# ---- the screen's table rule: LTRmix's load_library outputs on the same files (stored in the golden npz) ----
@pytest.fixture(scope="module")
def table_variants(tmp_path_factory):
    from builders import screen_table_variants

    return screen_table_variants(tmp_path_factory.mktemp("variants"))


@pytest.mark.parametrize("variant", ["as_is", "hash_header", "gz", "crlf_table", "no_k2p_column", "only_K2P_d",
                                     "upper_K2P", "header_only", "repeated_seq_id", "non_integer_ltr_len",
                                     "unplaced_with_row", "seq_id_without_class", "odd_values", "fasta_alone",
                                     "k2p_na_beside_K2P_d"])
def test_screen_reads_tables_as_ltrmix(table_variants, variant, caplog):
    """The screen's reading of a library's tables (io.library.read_library, then the screen's view of each record,
    routes.library.Records) gives, record by record, what LTRmix's load_library gave on the same files: names,
    chromosome, start, end, strand, K2P (k2p, else K2P_d; a column spelled K2P is not read), canonical TSD, the LTR
    lengths as integers (truncated; 0 unknown), the record-end LTR bounds and the stripped length. K2P_d is read only
    when the table has no k2p column: a row whose k2p is not a number is unknown even beside a K2P_d value. Rows join
    records by their exact name (an unplaced record by its row; a row naming the same interval under another name joins
    nothing); a repeated seq_id keeps its last row; a header-only table annotates nothing; LTRmix's warning when rows
    hold no K2P."""
    from synltr.routes.library import Records

    with np.load(Path(__file__).resolve().parent.parent / "golden" / "ltrmix_route.npz") as z:
        want = {f: z[f"rd_{variant}_{f}"] for f in ("names", "chrom", "start", "end", "strand", "k2p", "tsd", "l5",
                                                       "l3", "bounds", "stripped")}
    with caplog.at_level("WARNING", logger="synltr"):
        rec = Records.of("R", read_library(table_variants[variant], "R"))
    assert rec.names.tolist() == want["names"].tolist()
    assert [str(c) for c in rec.chrom] == want["chrom"].tolist()
    assert rec.start.tolist() == want["start"].tolist() and rec.end.tolist() == want["end"].tolist()
    assert rec.strand.tolist() == want["strand"].tolist() and rec.tsd.tolist() == want["tsd"].tolist()
    assert np.array_equal(rec.k2p, want["k2p"], equal_nan=True)
    assert rec.ltr_len.dtype == np.int64 and rec.ltr_len.tolist() == np.stack([want["l5"], want["l3"]], 1).tolist()
    got = [(-1, -1, -1, -1) if b is None else (*b[0], *b[1]) for b in map(rec.ltrs, range(rec.n))]
    assert got == [tuple(b) for b in want["bounds"].tolist()]
    assert [rec.seq(i).size for i in range(rec.n)] == want["stripped"].tolist()
    warned = "R: no K2P values: age strata disabled for this genome" in caplog.text
    assert warned == (variant in ("no_k2p_column", "upper_K2P", "seq_id_without_class"))


def test_screen_header_only_tables_in_prefix_form(table_variants, tmp_path, caplog):
    """A prefix whose tables hold only their header reads as the same files given as FASTA,TABLE (LTRmix's reading,
    above): every record unannotated (no K2P, TSD, strand or LTR lengths) and no K2P warning."""
    from synltr.routes.library import Records

    fa, tb = table_variants["header_only"].split(",")
    pre = tmp_path / "R_LTRs"
    shutil.copy(fa, f"{pre}_depth0_clean_ltr.fa")
    shutil.copy(tb, f"{pre}_depth0_clean_ltr.tsv")
    with caplog.at_level("WARNING", logger="synltr"):
        got = Records.of("R", read_library(str(pre), "R"))
    want = Records.of("R", read_library(table_variants["header_only"], "R"))
    for f in ("names", "chrom", "start", "end", "strand", "k2p", "tsd", "ltr_len"):
        assert np.array_equal(getattr(got, f), getattr(want, f), equal_nan=f == "k2p"), f
    assert np.isnan(got.k2p).all() and (got.tsd == "").all() and (got.strand == ".").all()
    assert (got.ltr_len == 0).all() and got.n > 10 and "no K2P values" not in caplog.text


def test_screen_bounds_at_the_stripped_ends_with_empty_and_masked_records(tmp_path):
    """Records empty, all N, with an N run and soft-masked: each record's N-stripped length decides its record-end
    LTRs (2 + 2 bp here): none where they do not fit, else at the stripped ends."""
    names = ["c1:1-4", "c1:11-14", "c1:21-30", "c1:31-38"]
    fa = write_fasta(tmp_path / "l.fa", dict(zip(names, ["", "NNNN", "ACGTNNACGT", "acgtACGT"])))
    tb = _screen_table(tmp_path / "l.tsv", [(n, 0.01, ".", "+", 2, 2) for n in names])
    lib = read_library(f"{fa},{tb}", "G")
    assert [lib.seq(i, stripped=True).size for i in range(4)] == [0, 0, 8, 8]
    assert lib.bounds.tolist() == ["none", "none", "ends", "ends"]
    assert lib.l5[2:].tolist() == [[0, 2], [0, 2]] and lib.l3[2:].tolist() == [[6, 8], [6, 8]]
