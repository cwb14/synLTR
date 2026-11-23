from pathlib import Path

import pytest

from synltr.io.samples import Sample, SampleSheetError, read_samples

HEADER = "label\tfasta\tannotation\tgroup\trole"


def write_sheet(tmp_path, rows, header=HEADER, name="samples.tsv"):
    (tmp_path / "inputs").mkdir(exist_ok=True)
    for fa in ("B73.fa", "Zd.fa", "Mo17.fa"):
        (tmp_path / "inputs" / fa).write_text(">c\nACGT\n")
    p = tmp_path / name
    p.write_text("\n".join([header, *rows]) + "\n")
    return p


GOOD = [
    "B73\tinputs/B73.fa\trun/B73_LTRs\tmaize\t.",
    "Mo17\tinputs/Mo17.fa\trun/Mo17_LTRs\tmaize\t.",
    "Zd\tinputs/Zd.fa\t.\tdiploperennis\toutgroup",
]


def test_read_valid_sheet(tmp_path):
    s = read_samples(write_sheet(tmp_path, ["# comment", "", *GOOD]))
    assert [x.label for x in s] == ["B73", "Mo17", "Zd"]
    assert s[0] == Sample("B73", tmp_path / "inputs/B73.fa", str(tmp_path / "run/B73_LTRs"), "maize", None)
    assert s[2].annotation is None and s[2].role == "outgroup"
    assert all(x.fasta.is_absolute() for x in s)


def test_optional_columns_and_column_order(tmp_path):
    p = write_sheet(tmp_path, ["inputs/B73.fa\tB73\trun/a", "inputs/Zd.fa\tZd\t."],
                    header="fasta\tlabel\tannotation")
    s = read_samples(p)
    assert s[0].group is None and s[0].role is None and s[1].annotation is None


def err(tmp_path, rows, match, header=HEADER):
    with pytest.raises(SampleSheetError, match=match):
        read_samples(write_sheet(tmp_path, rows, header=header))


def test_errors(tmp_path):
    err(tmp_path, [GOOD[0], GOOD[0]], "duplicate label.*B73")
    err(tmp_path, ["B 73\tinputs/B73.fa\trun/a\t.\t.", GOOD[2]], "B 73")
    err(tmp_path, ["a|b\tinputs/B73.fa\trun/a\t.\t.", GOOD[2]], r"a\|b")
    err(tmp_path, [GOOD[0], "X\tinputs/nope.fa\t.\t.\t.", "Y\tinputs/nope2.fa\t.\t.\t."], "nope.fa.*nope2.fa|nope2.fa.*nope.fa")
    err(tmp_path, [GOOD[0]], "at least 2")
    err(tmp_path, [GOOD[2], "Q\tinputs/B73.fa\t.\t.\t."], "annotation")
    err(tmp_path, [GOOD[0], "Zd\tinputs/Zd.fa\t.\t.\tdonor"], "donor")
    err(tmp_path, [GOOD[0], "Zd\tinputs/Zd.fa\t."], "line 3")
    err(tmp_path, [GOOD[0], GOOD[2]], "colour", header=HEADER + "\tcolour")
    err(tmp_path, [GOOD[0], GOOD[2]], "label", header="fasta\tannotation")


def test_label_pair_separator(tmp_path):
    # "__" joins two labels in cache file names (<A>__<B>.<route>.npz); a label may not contain it, nor
    # start or end with "_", which would make the join ambiguous ("a_" + "b" and "a" + "_b" both give "a___b")
    for bad in ("B73__x", "B73_", "_B73"):
        err(tmp_path, [f"{bad}\tinputs/B73.fa\trun/a\t.\t.", GOOD[2]], f"'{bad}' \\(line 2\\).*'__'")
    s = read_samples(write_sheet(tmp_path, ["B_7.3-x\tinputs/B73.fa\trun/a\t.\t.", GOOD[2]]))
    assert [x.label for x in s] == ["B_7.3-x", "Zd"]


def test_missing_sheet(tmp_path):
    with pytest.raises(SampleSheetError, match="not found"):
        read_samples(Path(tmp_path / "none.tsv"))


def test_sheet_details(tmp_path):
    """A byte-order mark and CRLF line ends are read; "." in group is none; a
    relative path resolves against the sheet's directory without following symlinks (a FASTA's index may sit
    beside its link); an empty FASTA field, duplicate header columns and the reserved label "all" fail."""
    p = write_sheet(tmp_path, GOOD)
    p.write_bytes(b"\xef\xbb\xbf" + p.read_bytes().replace(b"\n", b"\r\n"))
    s = read_samples(p)
    assert [x.label for x in s] == ["B73", "Mo17", "Zd"] and s[0].group == "maize"
    (tmp_path / "links").mkdir()
    (tmp_path / "links" / "B73.fa").symlink_to(tmp_path / "inputs" / "B73.fa")
    s = read_samples(write_sheet(tmp_path, ["B73\tlinks/../links/B73.fa\trun/a\t.\t.", GOOD[2]]))
    assert s[0].fasta == tmp_path / "links" / "B73.fa" and s[0].group is None
    err(tmp_path, ["B73\t.\trun/a\t.\t.", GOOD[2]], r"line 2 \(B73\) has no FASTA")
    err(tmp_path, [GOOD[0], GOOD[2]], "duplicate column", header=HEADER + "\tgroup")
    err(tmp_path, ["all\tinputs/B73.fa\trun/a\t.\t.", GOOD[2]], "'all' .* reserved")


def test_sheet_without_header_names_the_header(tmp_path):
    """A first line that is a data row is not a header: the error names the header the sheet needs."""
    (tmp_path / "samples.tsv").write_text("A\ta.fa\ta.gff3\nB\tb.fa\t.\n")
    with pytest.raises(SampleSheetError, match=r"first line.*header.*label.*fasta.*annotation"):
        read_samples(tmp_path / "samples.tsv")


H2 = "label\tfasta\tannotation\tgroup\trole\tlibrary\treads"


def test_library_and_reads_columns(tmp_path):
    (tmp_path / "x.bam").write_text("")
    rows = ["B73\tinputs/B73.fa\trun/B73_LTRs\tmaize\treference\t.\t.",
            "Mo17\tinputs/Mo17.fa\trun/Mo17_LTRs\tmaize\t.\tlib/fam.fa,lib/fam.tsv\t.",
            "R1\t.\t.\t.\t.\trun/lib\tx.bam"]
    s = read_samples(write_sheet(tmp_path, rows, header=H2))
    assert s[0].library is None and s[0].reads is None and s[0].role == "reference"
    assert s[1].library == f"{tmp_path / 'lib/fam.fa'},{tmp_path / 'lib/fam.tsv'}"
    assert s[2].fasta is None and s[2].reads == tmp_path / "x.bam" and s[2].library == str(tmp_path / "run/lib")
    err(tmp_path, [rows[0], "R1\t.\t.\t.\t.\t.\tnope.bam"], "reads file.*nope.bam", header=H2)


def test_row_without_fasta_ok_for_sheet(tmp_path):
    rows = ["L1\t.\t.\t.\t.\tlib1.fa\t.", "L2\t.\t.\t.\t.\tlib2.fa,t.tsv\t."]
    s = read_samples(write_sheet(tmp_path, rows, header=H2))
    assert [x.fasta for x in s] == [None, None]
    err(tmp_path, [rows[0], "L2\t.\t.\t.\t.\t.\t."], r"line 3 \(L2\) has no FASTA", header=H2)
    (tmp_path / "inputs" / "B73.fa.fai").write_text("c\t4\t3\t4\t5\n")
    s = read_samples(write_sheet(tmp_path, ["B73\tinputs/B73.fa.fai\trun/a\t.\t.\t.\t.", GOOD[2] + "\t.\t."],
                                 header=H2))
    assert s[0].fasta.name == "B73.fa.fai"


def test_reads_row_needs_reference(tmp_path):
    (tmp_path / "x.bam").write_text("")
    r1 = "R1\t.\t.\t.\t.\t.\tx.bam"
    ref = "B73\tinputs/B73.fa\trun/a\t.\treference\t.\t."
    other = "Mo17\tinputs/Mo17.fa\trun/b\t.\t.\t.\t."
    err(tmp_path, [r1, other], "reference", header=H2)
    err(tmp_path, [r1, ref.replace("run/a", ".")], "reference", header=H2)
    err(tmp_path, [r1, ref, ref.replace("B73", "Mo17")], "one .*reference", header=H2)
    assert len(read_samples(write_sheet(tmp_path, [r1, ref], header=H2))) == 2


def test_min_rows_keyword(tmp_path):
    p = write_sheet(tmp_path, [GOOD[0]])
    assert [x.label for x in read_samples(p, min_rows=1)] == ["B73"]
    with pytest.raises(SampleSheetError, match="at least 2"):
        read_samples(p)
