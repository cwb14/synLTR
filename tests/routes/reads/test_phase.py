"""Read phasing: distinct haplotypes carrying each allele (synthetic diploid with 1% heterozygosity)."""
import builders
import pysam
import pytest

from synltr.io.seqio import Fasta
from synltr.routes.reads.phase import haplotypes
from synltr.routes.reads.reads import evidence
from synltr.routes.reads.run import load_targets


@pytest.fixture(scope="module")
def het(tmp_path_factory):
    return builders.read_set(tmp_path_factory.mktemp("het"), 2, seed=4, het=0.01)


def _reads(d, eid):
    """(FASTA, element, labeled reads) of element `eid` in read set `d`."""
    fa = Fasta(d.ref)
    el = next(t for t in load_targets(str(d.gff), fa, "ref")[0] if t.id == eid)
    rd = []
    with pysam.AlignmentFile(str(d.bam)) as bam:
        evidence(bam, fa, el, el.tsd_len, 20, reads_out=rd)
    return fa, el, rd


def _h(d, eid, ploidy=2):
    fa, el, rd = _reads(d, eid)
    return haplotypes(rd, fa, el.chrom, el.start, el.end, ploidy)


def test_two_haplotypes_carry_full(het):
    h, n = _h(het, "E3_hom")
    assert h["F"] == 2 and n >= 2


def test_hemizygous_has_one(het):
    h, _ = _h(het, "E4_hemi")
    assert h["F"] <= 1


def test_heterozygous_empty(het):
    h, _ = _h(het, "E1_het_empty")
    assert h["F"] == 1 and h["E"] == 1


def test_identical_haplotypes_uninformative(data):
    h, _ = _h(data, "E3_hom")  # no heterozygosity: never more than one distinct haplotype
    assert h["F"] <= 1


def test_haploid_off(het):
    h, _ = _h(het, "E3_hom", ploidy=1)
    assert h == {"F": 0, "S": 0, "E": 0}


def test_low_quality_bases_are_ignored(het):
    """Bases below the quality threshold are neither reference nor alternative: with a threshold above every base
    quality there is nothing to phase."""
    fa, el, rd = _reads(het, "E3_hom")
    assert haplotypes(rd, fa, el.chrom, el.start, el.end, 2)[0]["F"] == 2
    assert haplotypes(rd, fa, el.chrom, el.start, el.end, 2, min_q=60) == ({"F": 0, "S": 0, "E": 0}, 0)
