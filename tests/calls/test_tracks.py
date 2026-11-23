import pytest
from builders import ABC_LENGTHS as LENGTHS
from builders import ABC_ORDER as ORDER
from builders import abc_samples, make_locus
from builders import genome_state as gs

from synltr import __version__
from synltr.calls.tracks import write_tracks

SAMPLES = abc_samples()


def test_tracks(tmp_path, loci):
    write_tracks(loci, SAMPLES, tmp_path, ORDER, LENGTHS)
    lines = (tmp_path / "A.pav.gff3").read_text().splitlines()
    assert lines[0] == "##gff-version 3"
    assert lines[1] == "##sequence-region chr1 1 1000" and not lines[2].startswith("##")  # only seqids used
    feats = [ln.split("\t") for ln in lines if not ln.startswith("#")]
    assert [f[2] for f in feats] == ["LTR_retrotransposon", "LTR_retrotransposon", "insertion_site"]
    assert feats[0][:8] == ["chr1", f"synltr-{__version__}", "LTR_retrotransposon", "1", "400", ".", "-",
                            "."]  # sorted by start
    assert {f[1] for f in feats} == {f"synltr-{__version__}"}  # the source column: the package and its version
    byid = {dict(kv.split("=") for kv in f[8].split(";"))["ID"]: f for f in feats}
    attrs = dict(kv.split("=", 1) for kv in byid["L0000001"][8].split(";"))
    assert attrs["states"] == "A:P,B:E,C:D" and attrs["state"] == "P" and attrs["tier"] == "annot"
    assert byid["L0000002"][2] == "insertion_site"
    assert "family=A%3BB%2CC" in byid["L0000002"][8]  # escaped ; and ,
    gb = [ln.split("\t") for ln in (tmp_path / "C.pav.gff3").read_text().splitlines() if not ln.startswith("#")]
    assert [f[2] for f in gb] == ["LTR_retrotransposon"]  # D without coordinates and U are omitted
    bed = [ln.split("\t") for ln in (tmp_path / "A.pav.bed").read_text().splitlines()]
    assert bed[0] == ["chr1", "0", "400", "L0000004", "0", "-", "P", "annot", "RLG_1"]
    b = [ln.split("\t") for ln in (tmp_path / "B.pav.bed").read_text().splitlines()]
    assert ["chr1", "49", "54", "L0000001", "0", "+", "E", "high", "RLG_1"] in b
    assert any(r[3] == "L0000003" and r[8] == "." for r in b)
    assert (tmp_path / "B.pav.gff3").read_text().count("deletion_junction") == 0


def test_tracks_deletion_junction(tmp_path):
    ls = [make_locus("L0000001", gs("P", "annot", "chr1", 5, 50, "+"), gs("D", "high", "chr2", 7, 8, "+"),
                     gs("U", "."))]
    write_tracks(ls, SAMPLES, tmp_path, ORDER, LENGTHS)
    assert "deletion_junction" in (tmp_path / "B.pav.gff3").read_text()


def test_tracks_follow_fai_order_and_escape_bed_columns(tmp_path):
    """Rows in the genome's FASTA index order (chr10 before chr2 here, against lexicographic
    order); a family with a tab or newline cannot break a BED line; an unknown seqid is a clear error."""
    order = {g: {"chr10": 0, "chr2": 1} for g in "ABC"}
    lengths = {g: {"chr10": 9000, "chr2": 9000} for g in "ABC"}
    ls = [make_locus("L0000001", gs("P", "annot", "chr2", 5, 50, "+"), gs("U", "."), gs("U", "."), family="a\tb\nc"),
          make_locus("L0000002", gs("P", "annot", "chr10", 900, 950, "+"), gs("U", "."), gs("U", "."))]
    write_tracks(ls, SAMPLES, tmp_path, order, lengths)
    bed = [ln.split("\t") for ln in (tmp_path / "A.pav.bed").read_text().splitlines()]
    assert [r[0] for r in bed] == ["chr10", "chr2"] and len(bed) == 2
    assert bed[1][8] == "a%09b%0Ac" and all(len(r) == 9 for r in bed)
    gff = (tmp_path / "A.pav.gff3").read_text().splitlines()
    assert gff[1:3] == ["##sequence-region chr10 1 9000", "##sequence-region chr2 1 9000"]
    bad = [make_locus("L0000003", gs("P", "annot", "chrX", 5, 50, "+"), gs("U", "."), gs("U", "."))]
    with pytest.raises(ValueError, match="L0000003: A's seqid 'chrX' is not in A's FASTA index"):
        write_tracks(bad, SAMPLES, tmp_path, order, lengths)
