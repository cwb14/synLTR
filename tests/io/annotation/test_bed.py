import logging
import random
from pathlib import Path

import builders as toy
import pytest

from ltrtrace.io.annotation.bed import find_ltrs, load_bed
from ltrtrace.io.seqio import Fasta


def test_find_ltrs_planted():
    rng = random.Random(1)
    ltr = toy.rand_seq(400, rng)
    el = ltr + toy.rand_seq(3000, rng) + toy.mutate(ltr, 0.005, rng)
    l5, l3 = find_ltrs(el)
    assert abs(l5 - 400) <= 5 and abs(l3 - 400) <= 5


def test_find_ltrs_random_is_none():
    assert find_ltrs(toy.rand_seq(4000, random.Random(2))) is None


def test_find_ltrs_too_short_is_none():
    assert find_ltrs("ACGT" * 20) is None


def _genome(tmp_path):
    rng = random.Random(3)
    ltr = toy.rand_seq(400, rng)
    el = ltr + toy.rand_seq(3000, rng) + toy.mutate(ltr, 0.005, rng)
    chr1 = toy.rand_seq(1000, rng) + el + toy.rand_seq(1000, rng)
    fa = tmp_path / "g.fa"
    fa.write_text(f">chr1\n{chr1}\n>chr2\n{toy.rand_seq(5000, rng)}\n")
    return fa, len(el)


def test_load_bed(tmp_path):
    fa, n = _genome(tmp_path)
    bed = tmp_path / "a.bed"
    bed.write_text(
        "# comment\ntrack name=x\n\n"
        f"chr1\t1000\t{1000 + n}\tTE1\t0\t-\n"
        f"chr1\t1000\t{1000 + n}\n"
        "chr2\t500\t4500\trand\n"
    )
    with Fasta(fa) as f:
        els, skipped = load_bed(bed, "G", f)
    assert skipped == {"no_ltrs": 1}
    assert len(els) == 2
    e, e2 = els
    assert (e.start, e.end, e.strand, e.name) == (1001, 1000 + n, "-", "TE1")
    assert abs(e.ltr_left[1] - 1400) <= 5 and e.ltr_left[0] == 1001
    assert abs(e.ltr_right[0] - (1000 + n - 399)) <= 5 and e.ltr_right[1] == e.end
    assert e.tsd is None and e.k is None and e.family is None and e.confidence is None
    assert e2.strand == "." and e2.name == f"chr1:1001-{1000 + n}"


def test_load_bed_bad_records(tmp_path):
    fa, _ = _genome(tmp_path)
    bed = tmp_path / "b.bed"
    bed.write_text("nochr\t0\t100\nchr2\t0\t999999\nchr2\t0\t100\n")
    with Fasta(fa) as f, pytest.raises(ValueError, match=r"2 BED record\(s\) not in the FASTA .*first: line 1 "
                                                         r"nochr:0-100; line 2 chr2:0-999999"):
        load_bed(Path(bed), "G", f)


def test_load_bed_gz_browser_line_and_reverse_element(tmp_path):
    """A gzipped BED with a browser line; an element stored reverse-complemented in the FASTA
    (its LTRs found on its own strand: the self-alignment is strand-blind)."""
    import gzip

    from ltrtrace.io.seqio import revcomp
    rng = random.Random(4)
    ltr = toy.rand_seq(400, rng)
    el = revcomp(ltr + toy.rand_seq(3000, rng) + toy.mutate(ltr, 0.005, rng))
    fa = tmp_path / "r.fa"
    fa.write_text(f">chr1\n{toy.rand_seq(1000, rng)}{el}{toy.rand_seq(1000, rng)}\n")
    bed = tmp_path / "r.bed.gz"
    with gzip.open(bed, "wt") as fh:
        fh.write(f"browser position chr1:1-100\nchr1\t1000\t{1000 + len(el)}\trev\t0\t-\n")
    with Fasta(fa) as f:
        [e], skipped = load_bed(bed, "G", f)
    assert not skipped and (e.name, e.strand) == ("rev", "-")
    assert abs(e.ltr_left[1] - e.ltr_left[0] + 1 - 400) <= 5 and abs(e.ltr_right[1] - e.ltr_right[0] + 1 - 400) <= 5


def test_load_bed_malformed_records_warn_and_skip(tmp_path, caplog):
    """A line with fewer than 3 columns or with coordinates that are not integers with
    0 <= start < end is skipped with a warning naming the file and line; CRLF line ends are read (strand
    "-\\r" is "-")."""
    fa, n = _genome(tmp_path)
    bed = tmp_path / "m.bed"
    bed.write_text(f"chr1\t1000\t{1000 + n}\tTE1\t0\t-\r\n"  # line 1: good, CRLF
                   "chr1\t9999\t15OOO\te1\n"  # line 2: not an integer
                   "chr1 1000 2000\n"  # line 3: spaces, one column
                   "chr2\t500\t500\tempty\n"  # line 4: start == end
                   "chr2\t-5\t100\tneg\n")  # line 5: start < 0
    with Fasta(fa) as f, caplog.at_level(logging.WARNING, logger="ltrtrace.io.annotation"):
        els, skipped = load_bed(bed, "G", f)
    assert [(e.name, e.strand) for e in els] == [("TE1", "-")]
    assert skipped == {"malformed_coordinates": 3, "malformed_columns": 1}
    assert [r.getMessage().split(": ")[0] for r in caplog.records] == [f"{bed}:{i}" for i in (2, 3, 4, 5)]
    assert "'15OOO'" in caplog.records[0].getMessage()


def test_load_bed_records_outside_the_fasta_name_their_lines(tmp_path):
    fa, _ = _genome(tmp_path)
    bed = tmp_path / "o.bed"
    bed.write_text("#c\nchr2\t0\t100\nnochr\t0\t100\n")
    with Fasta(fa) as f, pytest.raises(ValueError, match=r"o\.bed: 1 BED record.*first: line 3 nochr:0-100"):
        load_bed(bed, "G", f)


def test_gff3_one_based_bed_zero_based(tmp_path):
    """One element, GFF3 (1-based closed) to table to BED (0-based half-open): the BED start is the GFF3 start
    minus one and the end is unchanged; reading that BED line back gives the GFF3 interval."""
    from ltrtrace.io.annotation import load_annotation
    from ltrtrace.io.annotation.bed import from_bed, to_bed, write_bed
    from ltrtrace.io.annotation.gff3 import write_gff3
    from ltrtrace.io.elements import Element, ElementTable, make_eid
    e = Element("G", make_eid("G", "chr1", 11, 20), "x", "chr1", 11, 20, "-", (11, 12), (19, 20), None, None, None,
                None, None, None)
    gff = tmp_path / "a.gff3"
    write_gff3([e], gff)
    table = load_annotation(str(gff), "G")
    assert isinstance(table, ElementTable) and (table.row(0).start, table.row(0).end) == (11, 20)
    out = tmp_path / "a.bed"
    write_bed(table, out)
    assert out.read_text() == "chr1\t10\t20\tx\t0\t-\n"
    assert to_bed(11, 20) == (10, 20) and from_bed(10, 20) == (11, 20)
    assert table.row(0).zstart == 10
