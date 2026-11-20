import gzip
import logging
from dataclasses import replace

import pytest

from ltrtrace import __version__
from ltrtrace.base.constants import RECORD_WARNINGS
from ltrtrace.io.annotation import AnnotationError, load_annotation
from ltrtrace.io.annotation.gff3 import load_gff3, write_gff3
from ltrtrace.io.elements import NO_LTR, Element, make_eid

PHASING_K = ("K2P_d", "k2p", "K2P")

LTRQUEST = (
    "##gff-version 3\n"
    "B73v5_chr1\tLTRquest\tLTR_retrotransposon\t51712\t70182\t.\t-\t.\t"
    "ID=B73v5_LTRs_LTRRT_95050;Name=B73v5_chr1:51712-70182;classification=LTR/Copia/SIRE;order=LTR;"
    "superfamily=Copia;clade=SIRE;family=merged_fam00001;depth=1;identity=0.99661;K2P_d=0.00340137;"
    "motif=tg...ca;tsd=CACCT;confidence=high;strand_source=tesorter\n"
    "B73v5_chr1\tLTRquest\tlong_terminal_repeat\t51712\t53005\t.\t-\t.\t"
    "ID=B73v5_LTRs_LTRRT_95050.lLTR;Parent=B73v5_LTRs_LTRRT_95050\n"
    "B73v5_chr1\tLTRquest\tlong_terminal_repeat\t59462\t70182\t.\t-\t.\t"
    "ID=B73v5_LTRs_LTRRT_95050.rLTR;Parent=B73v5_LTRs_LTRRT_95050\n"
)

EDTA = (
    "##gff-version 3\n"
    "chr1\tEDTA\trepeat_region\t60659\t70135\t.\t-\t.\tID=repeat_region_1\n"
    "chr1\tEDTA\ttarget_site_duplication\t60659\t60663\t.\t-\t.\tID=lTSD_1;Parent=repeat_region_1\n"
    "chr1\tEDTA\tlong_terminal_repeat\t60664\t61985\t.\t-\t.\tID=lLTR_1;Parent=repeat_region_1;Name=x;"
    "Classification=LTR/Copia;ltr_identity=1.0000;Method=structural;motif=TGCA;tsd=GAGTG\n"
    "chr1\tEDTA\tCopia_LTR_retrotransposon\t60664\t70130\t.\t-\t.\tID=Copia_1;Parent=repeat_region_1;"
    "Name=x;Classification=LTR/Copia;Sequence_ontology=SO:0002264;ltr_identity=1.0000;"
    "Method=structural;motif=TGCA;tsd=gagtg\n"
    "chr1\tEDTA\tlong_terminal_repeat\t68826\t70130\t.\t-\t.\tID=rLTR_1;Parent=repeat_region_1;"
    "Classification=LTR/Copia;ltr_identity=1.0000;Method=structural\n"
    "chr1\tEDTA\ttarget_site_duplication\t70131\t70135\t.\t-\t.\tID=rTSD_1;Parent=repeat_region_1\n"
    "chr2\tEDTA\tLTR_retrotransposon\t100\t900\t.\t+\t.\tID=frag1;Name=frag1;Classification=LTR/Gypsy;"
    "Method=homology\n"
)


def test_ltrquest_convention(tmp_path):
    p = tmp_path / "a.gff3"
    p.write_text(LTRQUEST)
    (e,), skipped = load_gff3(p, "B73")
    assert (e.start, e.end, e.ltr_left, e.ltr_right, e.tsd) == (51712, 70182, (51712, 53005), (59462, 70182), "CACCT")
    assert abs(e.k - 0.00340137) < 1e-9 and (e.superfamily, e.clade) == ("Copia", "SIRE")
    assert (e.family, e.confidence, e.strand, e.seqid) == ("merged_fam00001", "high", "-", "B73v5_chr1")
    assert e.name == "B73v5_chr1:51712-70182" and e.eid == make_eid("B73", "B73v5_chr1", 51712, 70182)
    assert not skipped


def test_edta_convention_and_fragment(tmp_path):
    p = tmp_path / "e.gff3"
    p.write_text(EDTA)
    els, skipped = load_gff3(p, "B73")
    assert len(els) == 1 and els[0].ltr_right == (68826, 70130) and els[0].tsd == "GAGTG"
    assert els[0].ltr_left == (60664, 61985)
    assert els[0].k == 0.0 and els[0].superfamily == "Copia" and els[0].family == "x"
    assert els[0].clade is None and els[0].confidence is None
    assert skipped["no_ltr_pair"] == 1


def test_jc69_and_escapes_and_fasta_stop(tmp_path):
    p = tmp_path / "j.gff3"
    p.write_text(
        "##gff-version 3\n# comment\n"
        "c\tX\tLTR_retrotransposon\t1\t100\t.\t.\t.\tID=a%3Bb;ltr_identity=0.9\n"
        "c\tX\tlong_terminal_repeat\t1\t10\t.\t.\t.\tID=a1;Parent=a%3Bb\n"
        "c\tX\tlong_terminal_repeat\t91\t100\t.\t.\t.\tID=a2;Parent=a%3Bb\n"
        "##FASTA\n"
        "c\tX\tLTR_retrotransposon\t1\t50\t.\t.\t.\tID=after\n"
        "c\tX\tlong_terminal_repeat\t1\t5\t.\t.\t.\tID=z1;Parent=after\n"
        "c\tX\tlong_terminal_repeat\t45\t50\t.\t.\t.\tID=z2;Parent=after\n"
    )
    els, _ = load_gff3(p, "G")
    assert len(els) == 1  # the feature after ##FASTA is not read
    (e,) = els
    assert e.name == "a;b" and e.strand == "."
    assert abs(e.k - 0.107326) < 1e-5


def _els():
    a = Element("G", make_eid("G", "chr1", 10, 200), "n;1=x", "chr1", 10, 200, "+", (10, 30), (180, 200),
                "ACGTA", 0.1 + 0.2, "fam,1", "Gypsy", "Tekay", "high", tsd_source="annot")
    b = Element("G", make_eid("G", "chr2", 5, 99), "chr2:5-99", "chr2", 5, 99, "-", (5, 14), (90, 99),
                None, None, None, None, None, None)
    return [a, b]


def test_write_read_round_trip(tmp_path):
    els = _els()
    for name in ("o.gff3", "o.gff3.gz"):
        p = tmp_path / name
        write_gff3(els, p)
        got, skipped = load_gff3(p, "G")
        assert got == els and not skipped
    text = (tmp_path / "o.gff3").read_text().splitlines()
    assert text[0] == "##gff-version 3" and len(text) == 7
    assert "\tLTR_retrotransposon\t" in text[1] and "ID=G|chr1:10-200.lLTR;Parent=G|chr1:10-200" in text[2]
    assert {line.split("\t")[1] for line in text[1:]} == {f"ltrtrace-{__version__}"}  # the source column
    with gzip.open(tmp_path / "o.gff3.gz", "rt") as fh:
        assert fh.read() == (tmp_path / "o.gff3").read_text()


def test_gz_deterministic_and_parent_level_ignored(tmp_path):
    els = [replace(e, parent="x", level=2) for e in _els()]
    write_gff3(els, tmp_path / "1.gz")
    write_gff3(els, tmp_path / "2.gz")
    assert (tmp_path / "1.gz").read_bytes() == (tmp_path / "2.gz").read_bytes()
    got, _ = load_gff3(tmp_path / "1.gz", "G")
    assert got == _els()


def test_multi_parent_and_escaped_comma(tmp_path):
    p = tmp_path / "m.gff3"
    p.write_text(
        "c\tX\tLTR_retrotransposon\t1\t100\t.\t+\t.\tID=h1\n"
        "c\tX\tlong_terminal_repeat\t1\t10\t.\t+\t.\tID=l1;Parent=h1,other\n"
        "c\tX\tlong_terminal_repeat\t91\t100\t.\t+\t.\tID=l2;Parent=other,h1\n"
    )
    (e,), skipped = load_gff3(p, "G")
    assert (e.ltr_left, e.ltr_right) == ((1, 10), (91, 100)) and not skipped
    a = Element("G,x", make_eid("G,x", "chr,1", 10, 200), "n", "chr,1", 10, 200, "+", (10, 30), (180, 200),
                None, None, None, None, None, None)
    write_gff3([a], tmp_path / "c.gff3")
    got, skipped = load_gff3(tmp_path / "c.gff3", "G,x")
    assert got == [a] and not skipped


# ---------------------------------------------------------------- malformed records

GOOD = ("c\tX\tLTR_retrotransposon\t1001\t2000\t.\t+\t.\tID=ok\n"
        "c\tX\tlong_terminal_repeat\t1001\t1100\t.\t+\t.\tParent=ok\n"
        "c\tX\tlong_terminal_repeat\t1901\t2000\t.\t+\t.\tParent=ok\n")


def test_malformed_records_are_skipped_with_file_and_line(tmp_path, caplog):
    """A malformed LTR-RT or LTR line is skipped with a warning naming the file and line; a malformed line of a
    feature type the loader never uses (a gene) is not even read: no warning, no count, no abort."""
    p = tmp_path / "bad.gff3"
    p.write_text(
        "##gff-version 3\n"
        "c\tX\tgene\t20000\t2O100\t.\t+\t.\tID=g1\n"  # line 2: unused type, bad coordinate
        "c\tX\tmRNA\t5\n"  # line 3: unused type, too few columns
        + GOOD  # lines 4-6
        + "c\tX\tLTR_retrotransposon\t3001\t4OOO\t.\t+\t.\tID=b1\n"  # line 7: bad coordinate
        "c\tX\tLTR_retrotransposon\t6000\t5001\t.\t+\t.\tID=b2\n"  # line 8: start > end
        "c\tX\tLTR_retrotransposon\t7001\t8000\t.\t+\t.\tID=b3\n"  # line 9: its LTRs are malformed
        "c\tX\tlong_terminal_repeat\t7001\t7100\t.\t+\n"  # line 10: 7 columns
        "c\tX\tlong_terminal_repeat\t0\t8000\t.\t+\t.\tParent=b3\n"  # line 11: start < 1
    )
    with caplog.at_level(logging.WARNING, logger="ltrtrace.io.annotation"):
        els, skipped = load_gff3(p, "G")
    assert [e.name for e in els] == ["ok"]
    assert skipped == {"malformed_coordinates": 3, "malformed_columns": 1, "no_ltr_pair": 1}
    warned = [r.getMessage() for r in caplog.records]
    assert [w.split(": ")[0] for w in warned] == [f"{p}:{n}" for n in (7, 8, 10, 11)]
    assert "malformed_coordinates" in warned[0] and "'4OOO'" in warned[0] and "malformed_columns" in warned[2]


def test_malformed_record_warnings_are_capped(tmp_path, caplog):
    p = tmp_path / "many.gff3"
    p.write_text(GOOD + "".join(f"c\tX\tLTR_retrotransposon\t{i}x\t9\t.\t+\t.\tID=b{i}\n" for i in range(12)))
    with caplog.at_level(logging.WARNING, logger="ltrtrace.io.annotation"):
        els, skipped = load_gff3(p, "G")
    assert len(els) == 1 and skipped == {"malformed_coordinates": 12}
    assert len(caplog.records) == RECORD_WARNINGS and "only counted" in caplog.records[-1].getMessage()


def test_ltr_pair_negatives_crlf_nan_k2p_and_escaped_seqid(tmp_path):
    """Three LTRs or one LTR is no_ltr_pair; CRLF line ends are read; a K2P of nan is unknown; a seqid is
    percent-decoded (GFF3 escapes column 1)."""
    p = tmp_path / "n.gff3"
    p.write_text(
        "c\tX\tLTR_retrotransposon\t1\t100\t.\t+\t.\tID=three\n"
        + "".join(f"c\tX\tlong_terminal_repeat\t{s}\t{s + 9}\t.\t+\t.\tParent=three\n" for s in (1, 45, 91))
        + "c\tX\tLTR_retrotransposon\t201\t300\t.\t+\t.\tID=one\n"
        "c\tX\tlong_terminal_repeat\t201\t210\t.\t+\t.\tParent=one\n"
        "chr%2C9\tX\tLTR_retrotransposon\t401\t500\t.\t-\t.\tID=crlf;K2P_d=nan\r\n"
        "chr%2C9\tX\tlong_terminal_repeat\t401\t410\t.\t-\t.\tParent=crlf\r\n"
        "chr%2C9\tX\tlong_terminal_repeat\t491\t500\t.\t-\t.\tParent=crlf\r\n")
    (e,), skipped = load_gff3(p, "G")
    assert skipped == {"no_ltr_pair": 2}
    assert (e.seqid, e.strand, e.k, e.ltr_right) == ("chr,9", "-", None, (491, 500))


def test_edta_family_from_name(tmp_path):
    """EDTA keeps the family in Name and a unique ID: a record with Classification and no `family` takes Name
    as its family and ID as its name; superfamily comes from Classification."""
    p = tmp_path / "edta.gff3"
    p.write_text(
        "##gff-version 3\n"
        "chr1\tEDTA\trepeat_region\t100\t900\t.\t+\t.\tID=repeat_region_1\n"
        "chr1\tEDTA\tlong_terminal_repeat\t100\t199\t.\t+\t.\tID=lLTR_1;Parent=repeat_region_1;"
        "Name=TE_00000123;Classification=LTR/Gypsy\n"
        "chr1\tEDTA\tLTR_retrotransposon\t100\t900\t.\t+\t.\tID=Gypsy_1;Parent=repeat_region_1;"
        "Name=TE_00000123;Classification=LTR/Gypsy;Method=structural\n"
        "chr1\tEDTA\tlong_terminal_repeat\t801\t900\t.\t+\t.\tID=rLTR_1;Parent=repeat_region_1;"
        "Name=TE_00000123;Classification=LTR/Gypsy\n")
    (e,), skipped = load_gff3(p, "G")
    assert (e.family, e.superfamily, e.name) == ("TE_00000123", "Gypsy", "Gypsy_1")
    assert not skipped


# ---------------------------------------------------------------- types, rules and attribute keys

def _record(host_type="LTR_retrotransposon", attrs="ID=h", ltrs=((1, 10), (91, 100)), span=(1, 100), parent="h"):
    rows = [f"c\tX\t{host_type}\t{span[0]}\t{span[1]}\t.\t+\t.\t{attrs}\n"]
    rows += [f"c\tX\tlong_terminal_repeat\t{s}\t{e}\t.\t+\t.\tParent={parent}\n" for s, e in ltrs]
    return "".join(rows)


def test_dante_and_ltrharvest_types(tmp_path):
    """DANTE_LTR's transposable_element (read only with its Final_Classification) and LTRharvest's records (an
    LTR_retrotransposon under a repeat_region, TSDs as sibling features) are elements."""
    p = tmp_path / "d.gff3"
    p.write_text(
        "##gff-version 3\n"
        + _record("transposable_element", "ID=TE_1;Final_Classification=Class_I|LTR|Ty1/copia|Ale", parent="TE_1")
        + "c\tX\ttransposable_element\t300\t400\t.\t+\t.\tID=other\n"  # not an LTR-RT: no Final_Classification
        "c\tLTRharvest\trepeat_region\t1001\t2000\t.\t+\t.\tID=repeat_region1\n"
        "c\tLTRharvest\ttarget_site_duplication\t1001\t1005\t.\t+\t.\tParent=repeat_region1\n"
        "c\tLTRharvest\tLTR_retrotransposon\t1006\t1995\t.\t+\t.\tID=LTR_retrotransposon1;Parent=repeat_region1;"
        "ltr_similarity=98.5;seq_number=0\n"
        "c\tLTRharvest\tlong_terminal_repeat\t1006\t1105\t.\t+\t.\tParent=LTR_retrotransposon1\n"
        "c\tLTRharvest\tlong_terminal_repeat\t1896\t1995\t.\t+\t.\tParent=LTR_retrotransposon1\n"
        "c\tLTRharvest\ttarget_site_duplication\t1996\t2000\t.\t+\t.\tParent=repeat_region1\n")
    els, skipped = load_gff3(p, "G")
    assert [e.name for e in els] == ["TE_1", "LTR_retrotransposon1"] and not skipped
    d, h = els
    assert (d.superfamily, d.clade) == ("Ty1/copia", "Ale") and (d.start, d.end) == (1, 100)
    assert (h.ltr_left, h.ltr_right) == ((1006, 1105), (1896, 1995))


def _three_ltrs():
    return ("c\tX\tLTR_retrotransposon\t11\t90\t.\t+\t.\tID=h\n"
            + "".join(f"c\tX\tlong_terminal_repeat\t{s}\t{e}\t.\t+\t.\tParent=h\n"
                      for s, e in ((5, 20), (40, 50), (80, 95))))


def test_outermost_rule_matches_ltrkin(tmp_path):
    """Three LTRs: rule "two" skips the record; "outermost" takes the outermost two and extends the span to their
    ends, as ltrkin does."""
    p = tmp_path / "t.gff3"
    p.write_text(_three_ltrs())
    els, skipped = load_gff3(p, "G")
    assert not els and skipped["no_ltr_pair"] == 1
    (e,), skipped = load_gff3(p, "G", ltr_rule="outermost")
    assert not skipped
    assert (e.ltr_left, e.ltr_right) == ((5, 20), (80, 95)) and (e.start, e.end) == (5, 95)
    assert e.eid == make_eid("G", "c", 5, 95)


def test_outermost_rule_takes_siblings_when_one_ltr_is_a_child(tmp_path):
    p = tmp_path / "s.gff3"
    p.write_text(
        "c\tX\trepeat_region\t1\t100\t.\t+\t.\tID=rr\n"
        "c\tX\tLTR_retrotransposon\t1\t100\t.\t+\t.\tID=h;Parent=rr\n"
        "c\tX\tlong_terminal_repeat\t1\t10\t.\t+\t.\tParent=h\n"
        "c\tX\tlong_terminal_repeat\t91\t100\t.\t+\t.\tParent=rr\n")
    (e,), _ = load_gff3(p, "G", ltr_rule="outermost")
    assert (e.ltr_left, e.ltr_right) == ((1, 10), (91, 100))


def test_duplicate_element_ids_are_rejected(tmp_path):
    """Under the outermost rule a repeated element ID (concatenated GFF3 files; another chromosome's element)
    would merge two elements' LTRs, so it is an error."""
    p = tmp_path / "dup.gff3"
    p.write_text(_record().replace("c\t", "chr1\t") + _record().replace("c\t", "chr2\t"))
    with pytest.raises(ValueError, match="more than once"):
        load_gff3(p, "G", ltr_rule="outermost")
    with pytest.raises(AnnotationError, match=r"G: cannot read annotation .*more than once"):
        load_annotation(str(p), "G", ltr_rule="outermost")
    assert len(load_annotation(str(p), "G")) == 2  # rule "two" links LTRs per seqid: both are fine


def test_class_rules(tmp_path):
    p = tmp_path / "c.gff3"
    p.write_text(_record(attrs="ID=h;Classification=LTR/Gypsy"))
    (e,), _ = load_gff3(p, "G")
    assert (e.superfamily, e.clade) == ("Gypsy", None)
    (e,), _ = load_gff3(p, "G", class_rule="last")
    assert (e.superfamily, e.clade) == (None, "Gypsy")
    p.write_text(_record(attrs="ID=h;Classification=LTR/Copia/Ale"))
    (e,), _ = load_gff3(p, "G", class_rule="last")
    assert e.clade == "Ale"
    p.write_text(_record(attrs="ID=h;Classification=LTR/Copia/Ale;clade=SIRE;superfamily=Copia"))
    (e,), _ = load_gff3(p, "G", class_rule="last")
    assert (e.superfamily, e.clade) == ("Copia", "SIRE")  # explicit attributes win under both rules


def test_phasing_attribute_keys(tmp_path):
    """k2p_se, n_sites (aln_len when n_sites is absent), ti and tv (as n_ts, n_tv), k2p and K2P as k, and the LTRs'
    own lengths."""
    p = tmp_path / "k.gff3"
    p.write_text(_record(attrs="ID=h;k2p=0.05;k2p_se=0.01;aln_len=300;ti=7;tv=3;ltr5_len=8;ltr3_len=9"))
    (e,), _ = load_gff3(p, "G")
    assert e.k is None  # LTRpav reads K2P_d and ltr_identity only
    (e,), _ = load_gff3(p, "G", k_keys=PHASING_K)
    assert (e.k, e.k_kind, e.k_se, e.n_sites, e.n_ts, e.n_tv) == (0.05, "K2P", 0.01, 300.0, 7.0, 3.0)
    assert (e.ltr_left_len, e.ltr_right_len) == (8, 9)
    p.write_text(_record(attrs="ID=h;K2P=0.07;n_sites=250;aln_len=300"))
    (e,), _ = load_gff3(p, "G", k_keys=PHASING_K)
    assert (e.k, e.n_sites, e.k_se, e.n_ts) == (0.07, 250.0, None, None)
    assert (e.ltr_left_len, e.ltr_right_len) == (None, None)  # not given: unknown, not the extent lengths


def test_edta_identity_is_jc69(tmp_path):
    p = tmp_path / "j.gff3"
    p.write_text(_record(attrs="ID=h;ltr_identity=0.9"))
    (e,), _ = load_gff3(p, "G")
    assert e.k_kind == "JC69" and abs(e.k - 0.107326) < 1e-5
    p.write_text(_record(attrs="ID=h;K2P_d=0.1"))
    (e,), _ = load_gff3(p, "G")
    assert e.k_kind == "K2P"


def test_tsd_rules(tmp_path):
    """LTRpav keeps any annotated TSD (upper-cased); ltrkin's rule drops one that is not 3-8 bp of ACGT; phasing's
    accepts ACGTN of any length."""
    p = tmp_path / "t.gff3"
    cases = {"acgta": ("ACGTA", "ACGTA", "ACGTA"), "ACGTACGTA": ("ACGTACGTA", None, "ACGTACGTA"),
             "ACNTA": ("ACNTA", None, "ACNTA"), "AC": ("AC", None, "AC"), "AXGTA": ("AXGTA", None, None),
             "NA": (None, None, None), ".": (None, None, None)}
    for tsd, want in cases.items():
        p.write_text(_record(attrs=f"ID=h;tsd={tsd}"))
        for rule, w in zip(("keep", "clean", "acgtn"), want):
            (e,), _ = load_gff3(p, "G", tsd_rule=rule)
            assert e.tsd == w and e.tsd_source == ("annot" if w else None), (tsd, rule)


def test_tsd_key_upper_case_only_outside_ltrpavs_rule(tmp_path):
    p = tmp_path / "t.gff3"
    p.write_text(_record(attrs="ID=h;TSD=ACGTA"))
    assert load_gff3(p, "G")[0][0].tsd is None
    assert load_gff3(p, "G", tsd_rule="clean")[0][0].tsd == "ACGTA"


def test_name_family_matches_ltrkin(tmp_path):
    """ltrkin: family is the family or Family attribute, else the Name unless it is only coordinates; off (LTRpav),
    a Name is not a family."""
    p = tmp_path / "n.gff3"
    for attrs, off, on in (("ID=h;Name=Tekay_fam3", None, "Tekay_fam3"), ("ID=h;Name=chr1:1-100", None, None),
                           ("ID=h;Name=chr1:1..100", None, None), ("ID=h;Family=famB;Name=x", None, "famB"),
                           ("ID=h;family=famA;Name=x", "famA", "famA"), ("ID=h", None, None)):
        p.write_text(_record(attrs=attrs))
        assert load_gff3(p, "G")[0][0].family == off, attrs
        assert load_gff3(p, "G", name_family=True)[0][0].family == on, attrs


# ---------------------------------------------------------------- each tool's rule gives its original's result

def test_last_rule_keeps_literal_unknown_like_ltrkin(tmp_path):
    """ltrkin treats only a missing, empty or "." clade as absent: "unknown" is kept, and the last classification
    field is kept as written. LTRpav's rule maps "unknown" to none."""
    p = tmp_path / "u.gff3"
    p.write_text(_record(attrs="ID=h;clade=unknown;Classification=LTR/Copia/Ale"))
    assert load_gff3(p, "G", class_rule="last")[0][0].clade == "unknown"
    assert load_gff3(p, "G")[0][0].clade == "Ale"  # LTRpav: the attribute is unknown, the classification fills it
    assert load_gff3(p, "G")[0][0].superfamily == "Copia"
    p.write_text(_record(attrs="ID=h;Classification=LTR/unknown/unknown"))
    (e,), _ = load_gff3(p, "G", class_rule="last")
    assert e.clade == "unknown"
    (e,), _ = load_gff3(p, "G")
    assert (e.superfamily, e.clade) == (None, None)
    p.write_text(_record(attrs="ID=h;clade=.;Classification=LTR/Copia/Ale"))
    assert load_gff3(p, "G", class_rule="last")[0][0].clade == "Ale"


def test_name_family_is_ltrkins_truthiness_test(tmp_path):
    p = tmp_path / "f.gff3"
    for attrs, want in (("ID=h;family=unknown;Name=x", "unknown"), ("ID=h;family=.;Name=x", None),
                        ("ID=h;family=NA", "NA"), ("ID=h;Family=famB", "famB")):
        p.write_text(_record(attrs=attrs))
        assert load_gff3(p, "G", name_family=True)[0][0].family == want, attrs


def test_dante_classification_clade_is_the_last_field(tmp_path):
    p = tmp_path / "d.gff3"
    for cls, sf, clade in (("Class_I|LTR|Ty3/gypsy|chromovirus|Tekay", "Ty3/gypsy", "Tekay"),
                           ("Class_I|LTR|Ty3/gypsy|non-chromovirus|OTA|Athila", "Ty3/gypsy", "Athila"),
                           ("Class_I|LTR|Ty1/copia|Ale", "Ty1/copia", "Ale")):
        p.write_text(_record("transposable_element", f"ID=TE_1;Final_Classification={cls}", parent="TE_1"))
        (e,), _ = load_gff3(p, "G")
        assert (e.superfamily, e.clade) == (sf, clade), cls
        (e,), _ = load_gff3(p, "G", class_rule="last")
        assert e.clade == clade, cls


def test_final_classification_is_read_on_ltr_retrotransposons_only_under_last(tmp_path):
    p = tmp_path / "x.gff3"
    p.write_text(_record(attrs="ID=h;Final_Classification=LTR/Copia/Ale"))
    assert load_gff3(p, "G")[0][0].clade is None
    assert load_gff3(p, "G", class_rule="last")[0][0].clade == "Ale"


def test_default_types_are_case_sensitive_and_ltrkins_are_not(tmp_path):
    p = tmp_path / "c.gff3"
    p.write_text(_record("ltr_retrotransposon").replace("long_terminal_repeat", "Long_Terminal_Repeat"))
    assert load_gff3(p, "G")[0] == []
    (e,), _ = load_gff3(p, "G", ltr_rule="outermost")
    assert (e.ltr_left, e.ltr_right) == ((1, 10), (91, 100))
    p.write_text(_record("transposable_element"))  # no Final_Classification
    assert load_gff3(p, "G")[0] == [] and len(load_gff3(p, "G", ltr_rule="outermost")[0]) == 1


def test_k_key_precedence_and_nan_identity(tmp_path):
    p = tmp_path / "k.gff3"
    p.write_text(_record(attrs="ID=h;k2p=0.05;ltr_identity=0.9"))
    (e,), _ = load_gff3(p, "G")  # LTRpav: K2P_d, then ltr_identity; k2p is not read
    assert e.k_kind == "JC69" and abs(e.k - 0.107326) < 1e-5
    (e,), _ = load_gff3(p, "G", k_keys=PHASING_K)
    assert (e.k, e.k_kind) == (0.05, "K2P")
    p.write_text(_record(attrs="ID=h;K2P_d=0.2;ltr_identity=0.9"))
    assert load_gff3(p, "G")[0][0].k == 0.2
    p.write_text(_record(attrs="ID=h;ltr_identity=nan"))
    (e,), _ = load_gff3(p, "G")
    assert e.k != e.k and e.k_kind == "JC69"  # nan, as LTRpav's code returned it
    p.write_text(_record(attrs="ID=h;K2P_d=nan"))
    assert load_gff3(p, "G")[0][0].k is None


def test_keep_rule_is_ltrpavs_exactly(tmp_path):
    """No stripping and only the GFF3 unknown spellings: a blank-padded TSD is kept; "unknown" is none."""
    p = tmp_path / "t.gff3"
    p.write_text(_record(attrs="ID=h;tsd= acgta"))
    assert load_gff3(p, "G")[0][0].tsd == " ACGTA"
    p.write_text(_record(attrs="ID=h;tsd=unknown"))
    assert load_gff3(p, "G")[0][0].tsd is None
    assert load_gff3(p, "G", tsd_rule="acgtn")[0][0].tsd is None


def test_acgtn_rule_is_phasings(tmp_path):
    """Phasing reads only `tsd`, strips blanks, and drops "." and NA in any case."""
    p = tmp_path / "t.gff3"
    for attrs, want in (("tsd= acgnt ", "ACGNT"), ("tsd=na", None), ("tsd=Na", None), ("TSD=ACGTA", None)):
        p.write_text(_record(attrs=f"ID=h;{attrs}"))
        assert load_gff3(p, "G", tsd_rule="acgtn")[0][0].tsd == want, attrs


def test_clean_rule_is_ltrkins(tmp_path):
    p = tmp_path / "t.gff3"
    p.write_text(_record(attrs="ID=h;TSD=acgta"))
    assert load_gff3(p, "G", tsd_rule="clean")[0][0].tsd == "ACGTA"
    p.write_text(_record(attrs="ID=h;tsd= ACGTA"))  # a blank is not ACGT
    assert load_gff3(p, "G", tsd_rule="clean")[0][0].tsd is None


def test_repeated_coordinate_id_is_an_error_under_the_outermost_rule(tmp_path):
    """Elements without an ID are keyed by their coordinates, as in ltrkin: two at the same place are a repeat."""
    p = tmp_path / "d.gff3"
    p.write_text(_record(attrs="Name=x", parent="c:1-100") * 2)
    with pytest.raises(ValueError, match=r"c:1-100 occurs more than once"):
        load_gff3(p, "G", ltr_rule="outermost")


def test_unknown_rules_raise(tmp_path):
    p = tmp_path / "r.gff3"
    p.write_text(_record())
    for kw in ({"ltr_rule": "x"}, {"class_rule": "x"}, {"tsd_rule": "x"}, {"k_keys": ("nope",)},
               {"host_rule": "x"}):
        with pytest.raises(ValueError, match="must be"):
            load_gff3(p, "G", **kw)


PHASING_HOSTS = (
    "##gff-version 3\n"
    "c1\tx\tLTR_retrotransposon\t100\t1099\t.\t+\t.\tID=a;K2P_d=0.01\n"
    "c1\tx\tlong_terminal_repeat\t100\t149\t.\t+\t.\tID=a.l;Parent=a\n"
    "c1\tx\tlong_terminal_repeat\t1000\t1099\t.\t+\t.\tID=a.r;Parent=a\n"
    "c1\tx\tLTR_retrotransposon\t2000\t2999\t.\t+\t.\tID=b\n"
    "c1\tx\tlong_terminal_repeat\t2000\t2049\t.\t+\t.\tID=b.l;Parent=b\n"
    "c1\tx\tlong_terminal_repeat\t2950\t2999\t.\t+\t.\tID=b.r;Parent=b\n"
    "c1\tx\trepeat_region\t4000\t4999\t.\t-\t.\tID=c;k2p=0.02\n"
    "c1\tx\tLTR_retrotransposon\t6000\t6999\t.\t+\t.\tID=d;K2P=0.03\n"
    "c1\tx\tlong_terminal_repeat\t6000\t6049\t.\t+\t.\tID=d.l;Parent=d\n"
    "c1\tx\tCopia_LTR_retrotransposon\t8000\t8999\t.\t+\t.\tID=e;K2P_d=0.04\n"
    "c1\tx\tlong_terminal_repeat\t8000\t8049\t.\t+\t.\tID=e.l;Parent=e\n"
    "c1\tx\tlong_terminal_repeat\t8950\t8999\t.\t+\t.\tID=e.r;Parent=e\n"
)


def test_phasing_host_rule_and_optional_ltrs(tmp_path):
    """Phasing's GFF3 reading: host_rule "k2p" takes LTR_retrotransposon and repeat_region lines (those types
    exactly) that carry K2P_d, k2p or K2P; ltr_rule "optional" keeps a host without exactly two LTRs, with no LTR
    extents (NO_LTR)."""
    p = tmp_path / "p.gff3"
    p.write_text(PHASING_HOSTS)
    els, skipped = load_gff3(p, "G", ltr_rule="optional", host_rule="k2p", k_keys=PHASING_K)
    assert {e.name: (e.start, e.strand, e.k, e.ltr_left, e.ltr_right) for e in els} == {
        "a": (100, "+", 0.01, (100, 149), (1000, 1099)), "c": (4000, "-", 0.02, NO_LTR, NO_LTR),
        "d": (6000, "+", 0.03, NO_LTR, NO_LTR)}
    assert not skipped
    els, skipped = load_gff3(p, "G")                     # LTRpav's: types ending in LTR_retrotransposon, two LTRs
    assert [(e.name, e.k) for e in els] == [("a", 0.01), ("b", None), ("e", 0.04)]
    assert skipped == {"no_ltr_pair": 1}


def test_jc69_is_written_as_ltr_identity_and_read_back(tmp_path):
    """A JC69 distance is not written as K2P_d: the file holds the identity it came from, and reading it back keeps
    the label and the value."""
    a = Element("G", make_eid("G", "c", 1, 100), "a", "c", 1, 100, "+", (1, 10), (91, 100), None, 0.107326, None,
                None, None, None, k_kind="JC69")
    b = replace(a, eid=make_eid("G", "c", 201, 300), start=201, end=300, ltr_left=(201, 210), ltr_right=(291, 300),
                name="b", k_kind="K2P", k=0.05)
    p = tmp_path / "j.gff3"
    write_gff3([a, b], p)
    text = p.read_text()
    assert "ltr_identity=" in text and text.count("K2P_d=") == 1 and "K2P_d=0.05" in text
    got, _ = load_gff3(p, "G")
    assert (got[0].k_kind, got[1].k_kind) == ("JC69", "K2P")
    assert abs(got[0].k - 0.107326) < 1e-12 and got[1].k == 0.05


def test_id_name_rule_is_ltrkins(tmp_path):
    """name_rule "id" (ltrkin): an element's name is its ID attribute, else its line's coordinates seqid:start-end
    (before the outermost rule extends the span to its LTRs); "name" (the default, LTRpav) prefers the Name."""
    p = tmp_path / "a.gff3"
    p.write_text("##gff-version 3\n"
                 "chr1\tt\tLTR_retrotransposon\t100\t900\t.\t+\t.\tID=E1;Name=chr1:100-900\n"
                 "chr1\tt\tlong_terminal_repeat\t100\t200\t.\t+\t.\tParent=E1\n"
                 "chr1\tt\tlong_terminal_repeat\t800\t900\t.\t+\t.\tParent=E1\n"
                 "chr1\tt\tLTR_retrotransposon\t2001\t2900\t.\t-\t.\tfamily=f\n"
                 "chr1\tt\tlong_terminal_repeat\t1990\t2100\t.\t-\t.\tParent=chr1:2001-2900\n"
                 "chr1\tt\tlong_terminal_repeat\t2800\t2900\t.\t-\t.\tParent=chr1:2001-2900\n")
    els, _ = load_gff3(p, "G", ltr_rule="outermost", name_rule="id")
    assert [(e.name, e.start, e.end) for e in els] == [("E1", 100, 900), ("chr1:2001-2900", 1990, 2900)]
    assert [e.name for e in load_gff3(p, "G", ltr_rule="outermost")[0]] == ["chr1:100-900", "G|chr1:1990-2900"]
    assert [e.name for e in load_annotation(str(p), "G", ltr_rule="outermost", name_rule="id")] == [
        "E1", "chr1:2001-2900"]
    with pytest.raises(ValueError, match="name_rule"):
        load_gff3(p, "G", name_rule="Name")


PANSN = "HG002#1#chr1"  # a PanSN sequence name: '#' lies outside GFF3's seqid characters


def _write_io(d):
    """io.annotation.gff3.write_gff3 (sim pop's annotations): one element on PANSN."""
    e = Element(genome="G", eid=make_eid("G", PANSN, 101, 900), name="e;1", seqid=PANSN, start=101, end=900,
                strand="+", ltr_left=(101, 200), ltr_right=(801, 900), tsd="ACGTA", k=0.01, family="f=1",
                superfamily=None, clade=None, confidence=None)
    write_gff3([e], d / "io.gff3")
    return d / "io.gff3", True


def _write_tracks(d):
    """calls.tracks.write_tracks (<G>.pav.gff3): A's P site and B's E site on PANSN."""
    from builders import abc_samples, make_locus
    from builders import genome_state as gs

    from ltrtrace.calls.tracks import write_tracks

    loci = [make_locus("L0000001", gs("P", "annot", PANSN, 101, 900, "+"), gs("E", "high", PANSN, 50, 54, "+"),
                       gs("U", "."))]
    order, lengths = {g: {PANSN: 0} for g in "ABC"}, {g: {PANSN: 1000} for g in "ABC"}
    write_tracks(loci, abc_samples(), d, order, lengths)
    return d / "A.pav.gff3", True


def _write_exchanges(d):
    """analyses.subgenomes.painting.write_painting (exchanges.gff3): one called exchange on PANSN."""
    from types import SimpleNamespace

    from ltrtrace.analyses.subgenomes.painting import write_painting

    seg = {"seq": PANSN, "start": 101, "end": 900, "own": 0, "origin": 1, "left_flank": 50, "right_flank": 950,
           "n_ltr": 5, "n_support": 5, "n_against": 0, "n_lineages": 3, "llr": 20.0, "post": 0.9, "called": True}
    write_painting(SimpleNamespace(outdir=str(d)), {"segs": [seg], "res": [], "N": 2, "post": [], "n_needed": 10})
    return d / "exchanges.gff3", False                                         # regions: no LTR-RT to load


def _write_nucleus(d, monkeypatch):
    """sim.nucleus.write (annotation.gff3) of a tiny scenario whose sequences are named like PANSN."""
    from ltrtrace.sim.nucleus import write

    monkeypatch.setattr(write, "SIM_NUCLEUS_SEQID", "HG002#1#chr{:01d}", raising=False)
    sc = {"root_age": 2e6, "tree": {"split": 8e5, "children": ["A", "B"]}, "mergers": [{"time": 1e5,
          "lineages": ["A", "B"]}], "chromosomes": {"n": 1, "length": 200000}, "families": {"n": 4, "masters": [1, 1],
          "ltr_len_median": 120, "int_len_median": 1000}, "activity": {"rate": 4e-5}, "exchanges":
          {"nonreciprocal": 0, "reciprocal": 0}}
    assert write.simulate(sc, d / "lib", d / "truth", 7) > 0
    return d / "lib" / "annotation.gff3", True


@pytest.mark.parametrize("writer", ["io", "tracks", "exchanges", "nucleus"])
def test_every_gff3_writer_escapes_the_seqid(tmp_path, monkeypatch, writer):
    """Every GFF3 the package writes escapes a seqid as GFF3 requires (characters outside [a-zA-Z0-9.:^*$@!+_?|-] as
    %XX, in column 1 and in ##sequence-region), with the one helper the reader inverts: a PanSN name such as
    HG002#1#chr1 comes back whole from the package's GFF3 reader."""
    from ltrtrace.io.annotation import gff3

    make = {"io": _write_io, "tracks": _write_tracks, "exchanges": _write_exchanges,
            "nucleus": lambda d: _write_nucleus(d, monkeypatch)}[writer]
    path, loadable = make(tmp_path)
    lines = path.read_text().splitlines()
    named = [ln.split("\t")[0] for ln in lines if ln and not ln.startswith("#")]
    named += [ln.replace("\t", " ").split(" ")[1] for ln in lines if ln.startswith("##sequence-region")]
    assert named and all(n.startswith("HG002%231%23chr") for n in named), named
    assert {gff3.unescape(n)[:len(PANSN) - 1] for n in named} == {PANSN[:-1]}
    if loadable:
        seqids = {e.seqid for e in load_gff3(path, "G", ltr_rule="optional")[0]}
        assert seqids and all(s.startswith(PANSN[:-1]) for s in seqids), seqids


def test_escape_helpers_round_trip():
    """escape_seqid and escape_value give what GFF3 requires and unescape inverts both."""
    from ltrtrace.io.annotation.gff3 import escape_seqid, escape_value, unescape

    for text in (PANSN, "chr1", "seq..1(.)", "a b\té%"):
        assert unescape(escape_seqid(text)) == text and unescape(escape_value(text)) == text
    assert escape_seqid(PANSN) == "HG002%231%23chr1" and escape_seqid("chr1|x:y") == "chr1|x:y"
    assert escape_value("a;b=c&d,e%f\tg\rh") == "a%3Bb%3Dc%26d%2Ce%25f%09g%0Dh" and escape_value("x/y #1") == "x/y #1"
