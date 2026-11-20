"""ltrtrace.io.files: an input's compression told by its first bytes, never by its name."""
import gzip

import pysam
import pytest

from ltrtrace.io.files import gzip_kind, open_input


def test_compression_is_told_by_content(tmp_path):
    """Plain, gzip and bgzip files read the same text whatever their names, with or without .gz."""
    text = "a\tb\nc\td\n"
    plain, plain_gz_name = tmp_path / "x.tsv", tmp_path / "y.tsv.gz"
    plain.write_text(text)
    plain_gz_name.write_text(text)
    gz, gz_plain_name = tmp_path / "x.tsv.gz", tmp_path / "z.tsv"
    gz.write_bytes(gzip.compress(text.encode()))
    gz_plain_name.write_bytes(gz.read_bytes())
    bgz = tmp_path / "w.tsv.bgz"
    pysam.tabix_compress(str(plain), str(bgz))
    assert [gzip_kind(p) for p in (plain, plain_gz_name, gz, gz_plain_name, bgz)] == [None, None, "gzip", "gzip",
                                                                                       "bgzip"]
    for p in (plain, plain_gz_name, gz, gz_plain_name, bgz):
        with open_input(p) as fh:
            assert fh.read() == text, p
        with open_input(p, "rb") as fh:
            assert fh.read() == text.encode(), p
    with open_input(plain, newline="") as fh:
        assert fh.readline() == "a\tb\n"


def test_open_input_reads_only(tmp_path):
    (tmp_path / "x").write_text("")
    with pytest.raises(ValueError, match="mode"):
        open_input(tmp_path / "x", "wt")
    with pytest.raises(FileNotFoundError):
        open_input(tmp_path / "missing")
