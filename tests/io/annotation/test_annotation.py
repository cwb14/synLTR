import gzip
import logging
import random
from collections import Counter
from dataclasses import replace

import pytest

import ltrtrace.io.annotation as ann
from ltrtrace.io.annotation import AnnotationError, RecordSkips, assign_nesting, load_annotation
from ltrtrace.io.annotation.gff3 import write_gff3
from ltrtrace.io.annotation.nesting import containment
from ltrtrace.io.annotation.tsd import canon_tsd, infer_tsd
from ltrtrace.io.elements import Element, ElementTable, make_eid
from ltrtrace.io.seqio import Fasta


def _el(name, start, end, seqid="chr1", tsd=None, genome="G"):
    return Element(genome=genome, eid=make_eid(genome, seqid, start, end), name=name, seqid=seqid,
                   start=start, end=end, strand="+", ltr_left=(start, start + 9),
                   ltr_right=(end - 9, end), tsd=tsd, k=None, family=None, superfamily=None,
                   clade=None, confidence=None)


def _fasta(tmp_path, seqs):
    fa = tmp_path / "g.fa"
    fa.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
    return Fasta(fa)


def _nest(els):
    return {e.name: e for e in assign_nesting(ElementTable.from_rows(els))}


def test_nesting_and_levels():
    e1, e2, e3 = _el("e1", 100, 10000), _el("e2", 2000, 5000), _el("e3", 2500, 3000)
    out = _nest([e3, e1, e2])
    assert out["e1"].parent is None
    assert out["e2"].parent == out["e1"].eid and out["e3"].parent == out["e2"].eid
    assert [out[n].level for n in ("e1", "e2", "e3")] == [0, 1, 2]


def test_nesting_partial_overlap():
    out = _nest([_el("e1", 100, 10000), _el("e4", 9000, 12000)])
    assert out["e4"].parent is None and out["e4"].level == 0
    # e5 lies inside both e1 and e4; the topmost containing stack entry (e4) is its parent
    out = _nest([_el("e1", 100, 10000), _el("e4", 9000, 12000), _el("e5", 9100, 9500)])
    assert out["e5"].parent == out["e4"].eid and out["e5"].level == 1


def test_containment_counts_every_container():
    """ltrkin's counts: an element inside two partly overlapping ones has two containers (assign_nesting gives it
    level 1), and each container lists the spans of everything inside it, 1-based closed, in sweep order."""
    els = [_el("e5", 9100, 9500), _el("e1", 100, 10000), _el("x", 100, 1000, seqid="chr2"), _el("e4", 9000, 12000),
           _el("e3", 2500, 3000), _el("e2", 2000, 5000), _el("e6", 100, 900)]
    table = ElementTable.from_rows(els)
    n_inner, n_outer, inner = containment(table)
    got = {e.name: (n_inner[i], n_outer[i], inner[i]) for i, e in enumerate(table)}
    assert got["e1"] == (4, 0, [(100, 900), (2000, 5000), (2500, 3000), (9100, 9500)])
    assert got["e4"] == (1, 0, [(9100, 9500)]) and got["e2"] == (1, 1, [(2500, 3000)])
    assert got["e5"] == (0, 2, []) and got["e3"] == (0, 2, []) and got["e6"] == (0, 1, [])  # same start: e1 holds e6
    assert got["x"] == (0, 0, [])
    assert _nest(els)["e5"].level == 1


def test_containment_matches_its_definition():
    """On random intervals (no two alike): n_outer counts the others that hold the element (start <= its start, end
    >= its end), n_inner the others it holds."""
    rng = random.Random(3)
    spans = set()
    while len(spans) < 300:
        s = rng.randint(1, 50_000)
        spans.add((rng.choice(("c1", "c2")), s, s + rng.randint(10, 8000)))
    table = ElementTable.from_rows(_el(f"e{i}", s, e, seqid=c) for i, (c, s, e) in enumerate(sorted(spans)))
    n_inner, n_outer, inner = containment(table)
    rows = list(table)
    for i, e in enumerate(rows):
        holds = [o for o in rows if o is not e and o.seqid == e.seqid and o.start <= e.start and e.end <= o.end]
        held = [o for o in rows if o is not e and o.seqid == e.seqid and e.start <= o.start and o.end <= e.end]
        assert (n_outer[i], n_inner[i]) == (len(holds), len(held))
        assert sorted(inner[i]) == sorted((o.start, o.end) for o in held)


def test_nesting_separate_seqids_and_siblings():
    out = _nest([_el("a", 100, 1000), _el("b", 2000, 3000), _el("c", 100, 1000, seqid="chr2")])
    assert all(e.parent is None and e.level == 0 for e in out.values())


def test_infer_tsd(tmp_path):
    fa = _fasta(tmp_path, {"chr1": "GG" + "ACGTA" + "TTTTTTTTTT" + "ACGTA" + "CC",
                           "chr2": "GGCATG" + "TTTTTTTTTT" + "AAAAAC"})
    assert infer_tsd(fa, _el("x", 8, 17), (5, 4, 6)) == "ACGTA"
    assert infer_tsd(fa, _el("y", 7, 16, seqid="chr2"), (5, 4, 6)) is None


def test_infer_tsd_tries_5_then_4_then_6(tmp_path):
    fa = _fasta(tmp_path, {"chr1": "AAAAAA" + "TTTTTTTTTT" + "AAAAAA" + "G"})
    assert infer_tsd(fa, _el("x", 7, 16), (5, 4, 6)) == "AAAAA"  # 5 and 6 both match; 5 wins
    fa = _fasta(tmp_path, {"chr1": "CTGACT" + "TTTTTTTTTT" + "CTGACT" + "G"})
    assert infer_tsd(fa, _el("x", 7, 16), (5, 4, 6)) == "CTGACT"  # only 6 matches
    fa2 = _fasta(tmp_path, {"chr1": "AAAC" + "TTTTTTTTTT" + "AAACG"})
    assert infer_tsd(fa2, _el("x", 5, 14), (5, 4, 6)) == "AAAC"  # L=5 runs off the start; L=4 matches


def test_infer_tsd_at_sequence_start(tmp_path):
    fa = _fasta(tmp_path, {"chr1": "A" + "TTTTTTTTTT" + "ACGTA"})
    assert infer_tsd(fa, _el("x", 2, 11), (5, 4, 6)) is None
    fa = _fasta(tmp_path, {"chr1": "ACGTA" + "TTTTTTTTTT" + "ACG"})  # right flank too short
    assert infer_tsd(fa, _el("x", 6, 15), (5, 4, 6)) is None


def test_infer_tsd_rejects_n(tmp_path):
    fa = _fasta(tmp_path, {"chr1": "NNNNN" + "TTTTTTTTTT" + "NNNNN"})
    assert infer_tsd(fa, _el("x", 6, 15), (5, 4, 6)) is None


def _gff(tmp_path, elements, name="a.gff3"):
    p = tmp_path / name
    write_gff3(elements, p)
    return str(p)


def test_unknown_seqid_fails_fast(tmp_path):
    fa = _fasta(tmp_path, {"chr1": "A" * 100})
    spec = _gff(tmp_path, [_el("x", 10, 50, seqid="chrZ"), _el("y", 10, 50, seqid="chrY")])
    with pytest.raises(AnnotationError, match=r"chrZ.*2|2.*chrZ"):
        load_annotation(spec, "G", fa)


def test_missing_seqid_rule_skip_is_ltrkins(tmp_path, caplog):
    """missing_rule "skip" (ltrkin): the elements on sequences the FASTA lacks are dropped with ltrkin's one
    warning line (the first five names, sorted); the rest load and are still checked against their sequence ends;
    nothing left is an error. The default, "error", fails (test_unknown_seqid_fails_fast)."""
    fa = _fasta(tmp_path, {"chr1": "A" * 100})
    els = [_el(f"z{i}", 10, 50, seqid=f"chrU{i}") for i in range(6, -1, -1)] + [_el("x", 10, 50)]
    with caplog.at_level(logging.WARNING, logger="ltrtrace"):
        table = load_annotation(_gff(tmp_path, els), "G", fa, missing_rule="skip")
    assert [e.name for e in table] == ["x"]
    assert ("7 annotation sequence(s) not in the reference, skipped: ['chrU0', 'chrU1', 'chrU2', 'chrU3', 'chrU4']"
            in caplog.messages)
    with pytest.raises(AnnotationError, match=r"element y \(chr1:10-150\) ends beyond"):
        load_annotation(_gff(tmp_path, [_el("y", 10, 150), _el("z", 10, 50, seqid="chrZ")], "b.gff3"), "G", fa,
                        missing_rule="skip")
    with pytest.raises(AnnotationError, match=f"no annotated element is on a sequence of {fa.path}"):
        load_annotation(_gff(tmp_path, [_el("z", 10, 50, seqid="chrZ")], "c.gff3"), "G", fa, missing_rule="skip")
    with pytest.raises(ValueError, match="missing_rule"):
        load_annotation(_gff(tmp_path, [_el("x", 10, 50)], "d.gff3"), "G", fa, missing_rule="warn")


def test_end_beyond_sequence_fails(tmp_path):
    fa = _fasta(tmp_path, {"chr1": "A" * 100})
    spec = _gff(tmp_path, [_el("x", 10, 150)])
    with pytest.raises(AnnotationError, match=r"element x \(chr1:10-150\) ends beyond its sequence length 100"):
        load_annotation(spec, "G", fa)


def test_many_unknown_seqids_are_truncated(tmp_path):
    """The error names 5 unknown seqids and says there are more."""
    fa = _fasta(tmp_path, {"chr1": "A" * 100})
    spec = _gff(tmp_path, [_el(f"x{i}", 10, 50, seqid=f"chrU{i}") for i in range(7)])
    with pytest.raises(AnnotationError, match=r"7 annotation seqid\(s\) absent from the FASTA: chrU0, chrU1, "
                                              r"chrU2, chrU3, chrU4 \.\.\."):
        load_annotation(spec, "G", fa)


def test_nesting_past_a_buried_stale_entry():
    """An element that ended (e2) stays on the stack under a long-lived partial overlap (e3)
    until popped; a later element (e4) still finds its true container (e1), never the stale entry."""
    els = [_el("e1", 100, 10_000), _el("e2", 200, 300), _el("e3", 250, 9_000), _el("e4", 5_000, 5_100)]
    out = _nest(els)
    assert out["e2"].parent == out["e1"].eid and out["e3"].parent == out["e1"].eid
    assert out["e4"].parent == out["e3"].eid and out["e4"].level == 2


def test_skip_counts_without_a_counter(tmp_path, caplog):
    """Skip reasons are logged even when the caller passes no counter."""
    fa = _fasta(tmp_path, {"chr1": "A" * 200})
    spec = _gff(tmp_path, [_el("x", 10, 50), replace(_el("x2", 10, 50), eid="G|dup")])
    with caplog.at_level(logging.WARNING):
        assert len(load_annotation(spec, "G", fa)) == 1
    assert "G: skipped 1 annotation records (duplicate=1)" in caplog.text


def test_duplicates_dropped(tmp_path, caplog):
    fa = _fasta(tmp_path, {"chr1": "A" * 200})
    spec = _gff(tmp_path, [_el("x", 10, 50), replace(_el("x2", 10, 50), eid="G|dup"), _el("z", 60, 90)])
    skipped = RecordSkips()
    with caplog.at_level(logging.WARNING):
        out = load_annotation(spec, "G", fa, skipped=skipped)
    assert [(e.start, e.end) for e in out] == [(10, 50), (60, 90)]
    assert skipped.counts["duplicate"] == 1
    assert any("G: skipped 1 annotation records (duplicate=1)" in r.getMessage() for r in caplog.records)


def test_sorted_in_fai_order_with_nesting(tmp_path):
    fa = _fasta(tmp_path, {"chr2": "A" * 300, "chr1": "A" * 300})
    spec = _gff(tmp_path, [_el("a", 20, 100, seqid="chr1"), _el("b", 20, 200, seqid="chr1"),
                           _el("c", 50, 60, seqid="chr2")])
    out = load_annotation(spec, "G", fa)
    assert [e.name for e in out] == ["c", "b", "a"]
    assert out.row(1).level == 0 and out.row(2).parent == out.row(1).eid and out.row(2).level == 1


def test_tsd_filled_only_when_missing(tmp_path):
    fa = _fasta(tmp_path, {"chr1": "GG" + "ACGTA" + "T" * 10 + "ACGTA" + "CC"})
    spec = _gff(tmp_path, [_el("x", 8, 17)])
    assert load_annotation(spec, "G", fa).row(0).tsd == "ACGTA"
    spec = _gff(tmp_path, [_el("y", 8, 17, tsd="GGGGG")], "b.gff3")
    assert load_annotation(spec, "G", fa).row(0).tsd == "GGGGG"


def test_dispatch_by_format(tmp_path, monkeypatch):
    fa = _fasta(tmp_path, {"chr1": "A" * 200})
    e = _el("x", 10, 50)
    gz = _gff(tmp_path, [e], "a.gff3.gz")
    assert len(load_annotation(gz, "G", fa)) == 1
    assert len(load_annotation(_gff(tmp_path, [e], "a.gff"), "G", fa)) == 1
    bed = tmp_path / "a.bed"
    bed.write_text("chr1\t9\t50\n")
    seen = {}

    def fake_bed(path, genome, fasta):
        seen["bed"] = str(path)
        return [e], Counter({"no_ltr_pair": 2})

    def fake_lq(prefix, genome, **rules):
        seen["lq"] = prefix
        return [e], Counter()

    monkeypatch.setattr(ann, "load_bed", fake_bed)
    monkeypatch.setattr(ann, "load_ltrquest", fake_lq)
    skipped = RecordSkips()
    load_annotation(str(bed), "G", fa, skipped=skipped)
    assert seen["bed"] == str(bed) and skipped.counts["no_ltr_pair"] == 2
    prefix = tmp_path / "run"
    (tmp_path / "run_depth1_clean_ltr.tsv").write_text("x\n")
    load_annotation(str(prefix), "G", fa)
    assert seen["lq"] == str(prefix)
    with pytest.raises(AnnotationError, match="cannot detect annotation format"):
        load_annotation(str(tmp_path / "nothing"), "G", fa)
    with pytest.raises(AnnotationError, match="missing.gff3"):
        load_annotation(str(tmp_path / "missing.gff3"), "G", fa)
    with pytest.raises(AnnotationError, match="missing.bed"):
        load_annotation(str(tmp_path / "missing.bed"), "G", fa)


def test_annotation_files(tmp_path):
    """The files a spec names, as load_annotation reads them (they key the CLI's elements checkpoint)."""
    gff = tmp_path / "a.gff3.gz"
    gff.write_bytes(b"")
    assert ann.annotation_files(str(gff)) == [gff]
    assert ann.annotation_files(str(tmp_path / "missing.bed")) == []
    for depth in (1, 0):
        (tmp_path / f"run_depth{depth}_clean_ltr.tsv").write_text("x\n")
    assert ann.annotation_files(str(tmp_path / "run")) == [tmp_path / "run_depth0_clean_ltr.tsv",
                                                           tmp_path / "run_depth1_clean_ltr.tsv"]
    assert ann.annotation_files(str(tmp_path / "nothing")) == []
    table = tmp_path / "k.tsv"
    table.write_text("seq_id\tk2p\n")
    assert ann.annotation_files(str(table)) == [table]


def test_annotation_file_sniffed_when_not_gff3_or_bed(tmp_path):
    """A file named neither .gff3/.gff nor .bed: GFF3 when its first line starts ##gff or has GFF3's nine columns
    (and is no seq_id header), else a table read under table_rule."""
    els = [_el("x", 10, 50), _el("y", 100, 150)]
    gff = tmp_path / "ann.txt"
    write_gff3(els, gff)
    assert [e.name for e in load_annotation(str(gff), "G")] == ["x", "y"]
    gff.write_text("".join(l for l in gff.read_text().splitlines(keepends=True) if not l.startswith("#")))
    assert [e.name for e in load_annotation(str(gff), "G")] == ["x", "y"]
    t = tmp_path / "k.tsv"
    t.write_text("#seq_id\tk2p\tltr5_start\tltr5_end\tltr3_start\tltr3_end\tstrand\n"
                 "chr1:10-50\t0.01\t1\t10\t32\t41\t-\n")
    (e,) = load_annotation(str(t), "G", table_rule="phasing")
    assert (e.start, e.end, e.k, e.strand, e.ltr_left, e.ltr_right) == (10, 50, 0.01, "-", (10, 19), (41, 50))
    with pytest.raises(AnnotationError, match="status"):
        load_annotation(str(t), "G")                              # LTRquest's columns are required by default


def test_loader_errors_wrapped(tmp_path):
    fa = _fasta(tmp_path, {"chr1": "A" * 200})
    bed = tmp_path / "bad.bed"
    bed.write_text("chrQ\t1\t50\n")
    with pytest.raises(AnnotationError, match="chrQ") as ei:
        load_annotation(str(bed), "G", fa)
    assert isinstance(ei.value.__cause__, ValueError)


def test_no_usable_record_fails(tmp_path):
    """An annotation without one usable LTR-RT record (a wrong file, every record malformed or
    without its LTRs) fails with one line naming the genome, the file and the skip counts."""
    fa = _fasta(tmp_path, {"chr1": "A" * 200})
    gff = tmp_path / "genes.gff3"
    gff.write_text("##gff-version 3\nchr1\tX\tgene\t1\t100\t.\t+\t.\tID=g\n"
                   "chr1\tX\tLTR_retrotransposon\t1\t1OO\t.\t+\t.\tID=e\n")
    with pytest.raises(AnnotationError, match=r"G: annotation .*genes\.gff3 holds no usable LTR-RT record "
                                              r"\(malformed_coordinates=1\)"):
        load_annotation(str(gff), "G", fa)
    gff.write_text("##gff-version 3\n")
    with pytest.raises(AnnotationError, match="no LTR-RT record in it"):
        load_annotation(str(gff), "G", fa)


def test_unreadable_annotation_fails_with_one_line(tmp_path):
    """A truncated gzip GFF3 (EOFError) is an AnnotationError naming the genome and the file, never a traceback; a
    file named .gz that holds plain text is read as text (its first bytes decide), here finding no LTR-RT record."""
    fa = _fasta(tmp_path, {"chr1": "A" * 200})
    good = tmp_path / "a.gff3.gz"
    write_gff3([_el("x", 10, 50)], good)
    cut = tmp_path / "cut.gff3.gz"
    cut.write_bytes(good.read_bytes()[:-12])
    with pytest.raises(AnnotationError, match=r"G: cannot read annotation .*cut\.gff3\.gz: EOFError"):
        load_annotation(str(cut), "G", fa)
    fake = tmp_path / "fake.gff3.gz"
    fake.write_text("not gzip\n")
    with pytest.raises(AnnotationError, match=r"G: annotation .*fake\.gff3\.gz holds no usable LTR-RT record"):
        load_annotation(str(fake), "G", fa)


def test_annotation_compression_is_told_by_content(tmp_path):
    """Plain and gzip annotations read the same whatever their names, with or without .gz: the opener sniffs the
    first bytes."""
    fa = _fasta(tmp_path, {"chr1": "A" * 200})
    plain = tmp_path / "a.gff3"
    write_gff3([_el("x", 10, 50), _el("y", 100, 160)], plain)
    gz = tmp_path / "b.gff3.gz"
    gz.write_bytes(gzip.compress(plain.read_bytes()))
    gz_unnamed = tmp_path / "c.gff3"
    gz_unnamed.write_bytes(gz.read_bytes())
    plain_named_gz = tmp_path / "d.gff3.gz"
    plain_named_gz.write_bytes(plain.read_bytes())
    want = [(e.start, e.end, e.ltr_left, e.ltr_right) for e in load_annotation(str(plain), "G", fa)]
    assert len(want) == 2
    for f in (gz, gz_unnamed, plain_named_gz):
        assert [(e.start, e.end, e.ltr_left, e.ltr_right) for e in load_annotation(str(f), "G", fa)] == want, f


def test_seqid_mismatch_error_names_genome_and_file(tmp_path):
    fa = _fasta(tmp_path, {"chr1": "A" * 100})
    spec = _gff(tmp_path, [_el("x", 10, 50, seqid="scaf1")])
    with pytest.raises(AnnotationError) as exc:
        load_annotation(spec, "B73", fa)
    want = f"B73: annotation {spec}: 1 annotation seqid(s) absent from the FASTA: scaf1"
    assert str(exc.value).startswith(want)


def test_infer_tsd_orders(tmp_path):
    """A site whose 5-bp and 6-bp flanks both repeat: the order (5, 4, 6) returns 5 bp, (6, 5, 4) returns 6."""
    fa = _fasta(tmp_path, {"chr1": "AAAAAA" + "TTTTTTTTTT" + "AAAAAA" + "G"})
    e = _el("x", 7, 16)
    assert infer_tsd(fa, e, (5, 4, 6)) == "AAAAA"
    assert infer_tsd(fa, e, (6, 5, 4)) == "AAAAAA"
    fa = _fasta(tmp_path, {"chr1": "AAAC" + "TTTTTTTTTT" + "AAACG"})
    assert infer_tsd(fa, _el("x", 5, 14), (6, 5, 4)) == "AAAC"


def test_load_annotation_tsd_order(tmp_path):
    fa = _fasta(tmp_path, {"chr1": "AAAAAA" + "TTTTTTTTTT" + "AAAAAA" + "G"})
    spec = _gff(tmp_path, [_el("x", 7, 16)])
    t = load_annotation(spec, "G", fa)
    assert (t.row(0).tsd, t.row(0).tsd_source) == ("AAAAA", "inferred")
    assert load_annotation(spec, "G", fa, tsd_order=(6, 5, 4)).row(0).tsd == "AAAAAA"


def test_annotated_tsd_rules():
    from ltrtrace.io.annotation.tsd import annotated_tsd
    assert annotated_tsd("acgtacgtac") == "ACGTACGTAC" and annotated_tsd("acgtacgtac", "clean") is None
    assert annotated_tsd("AXGTA") == "AXGTA" and annotated_tsd("AXGTA", "acgtn") is None
    assert annotated_tsd("ACNTA", "acgtn") == "ACNTA" and annotated_tsd(None) is None


def test_load_annotation_tsd_rule(tmp_path):
    fa = _fasta(tmp_path, {"chr1": "ACGGTCATTGC" * 30})  # period 11: the flanks of 100-200 never repeat
    spec = _gff(tmp_path, [_el("x", 100, 200, tsd="ACGTACGTAC")])
    assert load_annotation(spec, "G", fa).row(0).tsd == "ACGTACGTAC"
    assert load_annotation(spec, "G", fa, tsd_rule="clean").row(0).tsd is None
    with pytest.raises(ValueError, match="tsd_rule"):
        load_annotation(spec, "G", fa, tsd_rule="x")


def test_unknown_rules_raise_at_every_entry(tmp_path):
    from ltrtrace.io.annotation.records import class_fields
    from ltrtrace.io.annotation.tsd import annotated_tsd
    with pytest.raises(ValueError, match="class_rule"):
        class_fields("LTR/Copia", "x")
    with pytest.raises(ValueError, match="tsd_rule"):
        annotated_tsd("ACGTA", "x")
    fa = _fasta(tmp_path, {"chr1": "A" * 200})
    spec = _gff(tmp_path, [_el("x", 10, 50)])
    for kw in ({"ltr_rule": "x"}, {"class_rule": "x"}, {"tsd_rule": "x"}, {"host_rule": "x"}, {"table_rule": "x"}):
        with pytest.raises(ValueError, match="must be"):
            load_annotation(spec, "G", fa, **kw)


def test_clean_tsd():
    from ltrtrace.io.annotation.tsd import clean_tsd
    assert clean_tsd("acgta") == "ACGTA" and clean_tsd("ACG") == "ACG" and clean_tsd("ACGTACGT") == "ACGTACGT"
    for bad in (None, "", ".", "NA", "AC", "ACGTACGTA", "ACNTA", "ACG TA"):
        assert clean_tsd(bad) is None


def test_ltr_and_class_rules_reach_the_loader(tmp_path):
    fa = _fasta(tmp_path, {"c": "A" * 200})
    gff = tmp_path / "r.gff3"
    gff.write_text("c\tX\tLTR_retrotransposon\t11\t90\t.\t+\t.\tID=h;Classification=LTR/Gypsy\n"
                   + "".join(f"c\tX\tlong_terminal_repeat\t{s}\t{e}\t.\t+\t.\tParent=h\n"
                             for s, e in ((5, 20), (40, 50), (80, 95))))
    with pytest.raises(AnnotationError, match="no usable"):
        load_annotation(str(gff), "G", fa)
    t = load_annotation(str(gff), "G", fa, ltr_rule="outermost", class_rule="last")
    assert (t.row(0).start, t.row(0).end, t.row(0).clade) == (5, 95, "Gypsy")


@pytest.mark.parametrize(("tsd", "want"), [
    ("ACGTA", "ACGTA"), ("TACGT", "ACGTA"),   # a TSD and its reverse complement are one site
    ("ttaca", "TGTAA"), ("TGTAA", "TGTAA"),   # case ignored
    ("ACGT", "ACGT"),                         # its own reverse complement
    ("", ""), (None, ""), (".", ""), ("NA", ""), ("ACGNA", ""), ("AC GT", ""),   # not a TSD: empty
])
def test_canon_tsd_is_the_smaller_of_the_two_orientations(tsd, want):
    assert canon_tsd(tsd) == want
