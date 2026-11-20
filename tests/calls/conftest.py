import pytest
from builders import genome_state as gs
from builders import make_locus

from ltrtrace.io.seqio import Fasta


@pytest.fixture
def fa(tmp_path):
    """Genome A's FASTA: chr1 (1000 bp) and chr2 (500 bp)."""
    p = tmp_path / "A.fa"
    p.write_text(">chr1\n" + "ACGT" * 250 + "\n>chr2\n" + "TTGCA" * 100 + "\n")
    return Fasta(p)


@pytest.fixture
def loci():
    """Four loci: present in the reference with B empty and C missing; empty in the reference with B present; a
    reference state D and one present at base 1 (both unplaced)."""
    return [
        make_locus("L0000001", gs("P", "annot", "chr1", 101, 5100, "+", 0.02), gs("E", "high", "chr1", 50, 54, "+"),
                   gs("D", "low"), members={"A": "A:e1"}),
        make_locus("L0000002", gs("E", "high", "chr1", 200, 204, "+"), gs("S", "medium", "chr2", 10, 5010, "-"),
                   gs("P", "annot", "chr1", 1, 50, "+", 0.04), members={"B": "B:e2", "C": "C:e3"},
                   discovery=["B", "C"], family="A;B,C"),
        make_locus("L0000003", gs("D", "low"), gs("P", "annot", "chr2", 5, 400, "+", 0.01),
                   gs("U", ".", flags=["no_call"]), members={"B": "B:e3"}, discovery=["B"], family=None),
        make_locus("L0000004", gs("P", "annot", "chr1", 1, 400, "-", 0.03), gs("U", "."), gs("U", "."),
                   members={"A": "A:e4"}),
    ]
