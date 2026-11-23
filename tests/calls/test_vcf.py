import gzip
import shutil
import subprocess
from pathlib import Path

import pysam
import pytest
from builders import abc_samples, make_locus, write_fasta
from builders import genome_state as gs

from synltr import __version__
from synltr.calls.states import INFERRED_ABSENT, OBSERVED_ABSENT
from synltr.calls.vcf import escape, write_reference_vcf, write_vcf
from synltr.io.seqio import Fasta
from synltr.routes.reads.genotype import Call, SampleCalls
from synltr.routes.reads.reads import Target

SAMPLES = abc_samples()


def test_vcf_encoding(tmp_path, loci, fa):
    vcf, _ = write_vcf(loci, SAMPLES, "A", fa, tmp_path / "pav.vcf.gz")
    v = pysam.VariantFile(str(vcf))
    recs = list(v)
    assert [r.id for r in recs] == ["L0000001", "L0000002"]
    rec = recs[0]
    assert rec.alts == ("<INS:ME:LTR>",) and rec.info["SVTYPE"] == "INS"
    assert [rec.samples[s]["GT"] for s in ("A", "B", "C")] == [(1,), (0,), (None,)]
    assert [rec.samples[s]["ST"] for s in ("A", "B", "C")] == ["P", "E", "D"]
    assert [rec.samples[s]["TR"] for s in ("A", "B", "C")] == ["annot", "high", "low"]
    assert rec.info["AA"] == rec.ref
    assert rec.info["NE"] == 1 and rec.info["SVLEN"] == 5000 and rec.info["FAM"] == "RLG_1"
    assert rec.info["K"] == pytest.approx(0.02) and rec.info["DISC"] == ("A",)
    assert Path(str(vcf) + ".tbi").is_file()
    assert list(v.header.contigs) == ["chr1", "chr2"] and v.header.contigs["chr2"].length == 500
    [note] = [r.value for r in v.header.records if r.key == "LTRpav_POS"]  # the E site's POS
    assert [r.value for r in v.header.records if r.key == "source"] == [f"synltr {__version__}"]
    assert "without a readable TSD copy, the flank junction, to within 1 bp" in note


def test_vcf_pos_conventions(tmp_path, loci, fa):
    vcf, _ = write_vcf(loci, SAMPLES, "A", fa, tmp_path / "pav.vcf.gz")
    r1, r2 = list(pysam.VariantFile(str(vcf)))
    assert r1.pos == 100 and r1.ref == "ACGT"[99 % 4]  # present: POS = start - 1
    assert r2.pos == 204 and r2.ref == "ACGT"[203 % 4]  # E: last base of the site copy
    assert r2.info["K"] == pytest.approx(0.04) and "AA" in r2.info and r2.info["DISC"] == ("B", "C")
    assert r2.info["FAM"] == "A%3BB%2CC"
    assert r1.stop == r1.pos
    assert [r2.samples[s]["GT"] for s in "ABC"] == [(0,), (1,), (1,)]


def test_unplaced_side_file(tmp_path, loci, fa):
    vcf, unp = write_vcf(loci, SAMPLES, "A", fa, tmp_path / "pav.vcf.gz")
    rows = [ln.split("\t") for ln in Path(unp).read_text().splitlines()]
    assert rows[0] == ["lid", "ref_state", "reason"]
    assert [r[:2] for r in rows[1:]] == [["L0000003", "D"], ["L0000004", "P"]]
    assert "L0000003" not in [r.id for r in pysam.VariantFile(str(vcf))]


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools missing")
def test_bcftools_reads_vcf(tmp_path, loci, fa):
    vcf, _ = write_vcf(loci, SAMPLES, "A", fa, tmp_path / "pav.vcf.gz")
    out = subprocess.run(["bcftools", "view", "-H", str(vcf)], capture_output=True, text=True, check=True)
    assert len(out.stdout.splitlines()) == 2
    q = subprocess.run(["bcftools", "query", "-f", "%ID\t[%GT ]\n", str(vcf)], capture_output=True, text=True,
                       check=True)
    assert q.stdout.splitlines()[0].split("\t")[1].split() == ["1", "0", "."]


def test_inferred_absence_is_genotyped_only_when_accepted(tmp_path, fa):
    """A is missing by default (as for the statistics validated on E); with the inferred set it is absent."""
    ls = [make_locus("L0000001", gs("P", "annot", "chr1", 101, 500, "+"), gs("A", "high"), gs("E", "high", "chr1", 9, 13))]
    for absent, want in ((OBSERVED_ABSENT, [(1,), (None,), (0,)]), (OBSERVED_ABSENT | INFERRED_ABSENT, [(1,), (0,), (0,)])):
        vcf, _ = write_vcf(ls, SAMPLES, "A", fa, tmp_path / f"{len(absent)}.vcf.gz", absent=absent)
        [rec] = list(pysam.VariantFile(str(vcf)))
        assert [rec.samples[s]["GT"] for s in "ABC"] == want
        assert [rec.samples[s]["ST"] for s in "ABC"] == ["P", "A", "E"]


# ---- the reference VCF of the reads route (ltrkin's record type)

def _target(eid, chrom, start, length, ltr=100, tsd="ACGTA", family="fam1", clade="Ale"):
    """A reads-route element, 0-based half-open: `start` is its first base, its LTRs `ltr` bp each."""
    end = start + length
    return Target(eid, chrom, start, end, "+", (start, start + ltr), (end - ltr, end), family, clade, tsd)


def _records(path: Path) -> list[list[str]]:
    with gzip.open(path, "rt") as fh:
        return [line.rstrip("\n").split("\t") for line in fh if not line.startswith("#")]


def _info(rec: list[str]) -> dict[str, str]:
    return dict(kv.split("=", 1) for kv in rec[7].split(";"))


def _one_sample(targets, gt="0/1"):
    return [SampleCalls("S1", 2, {i: Call(gt, 40, (10, 5, 0), ()) for i in range(len(targets))})]


def test_vcf_end_matches_svlen(tmp_path, fa):
    """END is POS plus the length of the reference sequence the empty allele lacks (the element and one TSD
    copy), the first ALT's |SVLEN|."""
    targets = [_target("e1", "chr1", 99, 200), _target("e2", "chr1", 499, 300, tsd="")]  # e2: no TSD, 5 bp assumed
    vcf, _ = write_reference_vcf(tmp_path / "S1.vcf.gz", targets, _one_sample(targets), fa, tmp_path / "A.fa")
    r1, r2 = _records(vcf)
    assert (r1[1], _info(r1)["END"], _info(r1)["SVLEN"]) == ("99", "304", "-205,-100")
    assert (r2[1], _info(r2)["END"], _info(r2)["SVLEN"]) == ("499", "804", "-305,-200")
    for rec in (r1, r2):
        assert int(_info(rec)["END"]) - int(rec[1]) == -int(_info(rec)["SVLEN"].split(",")[0])
    assert [r.stop for r in pysam.VariantFile(str(vcf))] == [304, 804]


def test_vcf_pos_at_base_one(tmp_path, fa):
    """An element starting at base 1 has no base before it to anchor its record on: it is left out of the VCF and
    listed in the unplaced table, reason starts_at_base_1, as the catalog VCF does. A later element keeps the base
    before it."""
    targets = [_target("first", "chr1", 0, 200), _target("next", "chr1", 300, 200)]
    vcf, unplaced = write_reference_vcf(tmp_path / "S1.vcf.gz", targets, _one_sample(targets), fa, tmp_path / "A.fa")
    (nxt,) = _records(vcf)
    assert nxt[1:4] == ["300", "next", "ACGT"[299 % 4]]
    assert [r.id for r in pysam.VariantFile(str(vcf))] == ["next"]
    assert unplaced == tmp_path / "S1.unplaced.tsv"
    assert unplaced.read_text() == "lid\tref_state\treason\nfirst\tP\tstarts_at_base_1\n"
    with gzip.open(vcf, "rt") as fh:
        assert "##contig=<ID=chr1,length=1000>\n" in fh.read()  # its sequence still holds an element


def test_vcf_fai_order(tmp_path):
    """Contig lines and records follow the reference's .fai order, not the text order of the names; only contigs
    holding elements get a line; equal POS keep the element order."""
    fa = Fasta(write_fasta(tmp_path / "R.fa", {"chrB": "ACGT" * 300, "chrX": "ACGT" * 10, "chrA": "TTGCA" * 200}))
    targets = [_target("a2", "chrA", 500, 200), _target("b1", "chrB", 100, 300), _target("a1", "chrA", 50, 150),
               _target("b_long", "chrB", 600, 400), _target("b_short", "chrB", 600, 250)]
    vcf, _ = write_reference_vcf(tmp_path / "S1.vcf.gz", targets, _one_sample(targets), fa, tmp_path / "R.fa")
    with gzip.open(vcf, "rt") as fh:
        contigs = [line.strip() for line in fh if line.startswith("##contig")]
    assert contigs == ["##contig=<ID=chrB,length=1200>", "##contig=<ID=chrA,length=1000>"]
    assert [r[2] for r in _records(vcf)] == ["b1", "b_long", "b_short", "a1", "a2"]
    assert [r.id for r in pysam.VariantFile(str(vcf)).fetch("chrA")] == ["a1", "a2"]


def test_reference_vcf_records(tmp_path, fa):
    """ltrkin's record: REF the base before the element, ALT the empty site and the solo LTR, FORMAT GT:GQ:AD:FT
    (FT flags joined by ';', PASS without any); a sample without data at an element gets a no-data call of its
    ploidy; an element no sample holds is left out; ID, family and clade are percent-encoded."""
    targets = [_target("e;1", "chr1", 99, 200, ltr=50, family="fam=1", clade="A le"), _target("e2", "chr2", 9, 100),
               _target("none", "chr2", 300, 100)]
    samples = [SampleCalls("S1", 2, {0: Call("0/1", 40, (10, 5, 0), ("low_depth", "low_mapq")),
                                     1: Call("./.", 0, (0, 0, 0), ("no_reads",))}),
               SampleCalls("S2", 4, {0: Call("0/0/0/2", 99, (30, 0, 7), ())})]
    vcf, unplaced = write_reference_vcf(tmp_path / "S.vcf.gz", targets, samples, fa, tmp_path / "A.fa")
    assert unplaced.read_text() == "lid\tref_state\treason\n"
    with gzip.open(vcf, "rt") as fh:
        head = [line.rstrip("\n") for line in fh if line.startswith("#")]
    assert head[0] == "##fileformat=VCFv4.3" and head[-1].split("\t")[9:] == ["S1", "S2"]
    assert f"##reference={tmp_path / 'A.fa'}" in head and f"##source=synltr {__version__}" in head
    assert any(line.startswith("##ALT=<ID=DEL:ME:SOLO,") for line in head)
    r1, r2 = _records(vcf)
    assert r1[2:9] == ["e%3B1", "ACGT"[98 % 4], "<DEL:ME:LTR>,<DEL:ME:SOLO>", ".", ".",
                       "SVTYPE=DEL;END=304;SVLEN=-205,-150;ELEMENT=e%3B1;FAMILY=fam%3D1;CLADE=A%20le", "GT:GQ:AD:FT"]
    assert r1[9:] == ["0/1:40:10,5,0:low_depth;low_mapq", "0/0/0/2:99:30,0,7:PASS"]
    assert r2[9:] == ["./.:0:0,0,0:no_reads", "./././.:0:0,0,0:no_data"]
    assert [r.id for r in pysam.VariantFile(str(vcf))] == ["e%3B1", "e2"]


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools missing")
def test_bcftools_index_merge_and_view_lose_nothing(tmp_path, fa):
    """Two reads-route VCFs (one with an element at base 1, which goes to the unplaced table): bcftools index, merge
    and indexed region views keep every record either VCF holds, and print no warning."""
    targets = [_target("first", "chr1", 0, 200), _target("e1", "chr1", 299, 200), _target("e2", "chr2", 9, 100),
               _target("e3", "chr2", 199, 150)]
    s1 = [SampleCalls("S1", 2, {i: Call("0/1", 40, (10, 5, 0), ()) for i in range(4)})]
    s2 = [SampleCalls("S2", 2, {i: Call("0/0", 30, (20, 0, 0), ()) for i in (0, 2, 3)})]
    vcfs = [write_reference_vcf(tmp_path / f"{s[0].label}.vcf.gz", targets, s, fa, tmp_path / "A.fa")[0]
            for s in (s1, s2)]

    def bcftools(*args):
        r = subprocess.run(["bcftools", *map(str, args)], capture_output=True, text=True, check=True)
        assert r.stderr == "", r.stderr
        return r.stdout

    def ids(text):
        return [line.split("\t")[2] for line in text.splitlines() if not line.startswith("#")]

    for v in vcfs:
        bcftools("index", "-f", v)
    assert ids(bcftools("view", "-H", vcfs[0])) == ["e1", "e2", "e3"] and ids(bcftools("view", "-H", vcfs[1])) == [
        "e2", "e3"]
    merged = bcftools("merge", "-Ov", *vcfs)
    assert ids(merged) == ["e1", "e2", "e3"]
    assert next(line.split("\t")[9:] for line in merged.splitlines() if not line.startswith("#")) == [
        "0/1:40:10,5,0:PASS", "./.:.:.:."]
    for contig, want in (("chr1", ["e1"]), ("chr2", ["e2", "e3"])):
        assert ids(bcftools("view", "-H", "-r", contig, vcfs[0])) == want
    assert ids(bcftools("view", "-H", "-r", "chr1:1-50", vcfs[0])) == []  # nothing is anchored before base 1


def test_escape_is_shared_and_encodes_every_whitespace():
    """One escaping for both record types: % ; = , and whitespace as %XX (VCF 4.3, 1.2)."""
    assert escape("a%b;c=d,e f\tg\nh\ri") == "a%25b%3Bc%3Dd%2Ce%20f%09g%0Ah%0Di"
    assert escape("x\x0cy\x0bz") == "x%0Cy%0Bz" and escape(5) == "5" and escape("Ale") == "Ale"
    assert escape("a\u2028b\u3000c\xa0d\x1ce") == "a%2028b%3000c%A0d%1Ce"  # four hex digits above U+00FF
