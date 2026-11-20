import gzip
import os

import numpy as np
import pysam
import pytest

from ltrtrace.base.log import UsageError
from ltrtrace.io.seqio import Fasta, encode, n_runs, open_fasta, read_fai, revcomp, revcomp_codes, stream_fasta


def test_fetch_plain_and_bgzip(tmp_path):
    fa = tmp_path / "g.fa"
    fa.write_text(">s1\nACGTNNNNacgt\n>s2\nGGGG\n")
    gz = tmp_path / "g2.fa.gz"
    pysam.tabix_compress(str(fa), str(gz))
    for p in (fa, gz):
        with Fasta(p, cache_dir=tmp_path / "cache") as f:
            assert f.lengths() == {"s1": 12, "s2": 4} and (f.length("s1"), f.length("s2")) == (12, 4)
            assert f.fetch("s1", 0, 12) == "ACGTNNNNACGT"
            assert f.fetch("s1", 10, 99) == "GT"


def test_fetch_uppercases(tmp_path):
    fa = tmp_path / "g.fa"
    fa.write_text(">s1\nacgtACGT\n")
    with Fasta(fa) as f:
        assert f.fetch("s1", 0, 4) == "ACGT"


def test_fetch_clipping_and_unknown(tmp_path):
    fa = tmp_path / "g.fa"
    fa.write_text(">s1\nACGT\n")
    with Fasta(fa) as f:
        assert f.fetch("s1", -5, 2) == "AC"
        assert f.fetch("s1", 3, 3) == ""
        assert f.fetch("s1", 9, 12) == ""
        with pytest.raises(KeyError):
            f.fetch("nope", 0, 1)
        with pytest.raises(KeyError):
            f.length("nope")


@pytest.mark.skipif(os.geteuid() == 0, reason="root ignores directory permissions")
def test_index_into_cache_when_dir_readonly(tmp_path):
    d = tmp_path / "ro"
    d.mkdir()
    fa = d / "g.fa"
    fa.write_text(">s1\nACGT\n")
    d.chmod(0o555)
    try:
        with Fasta(fa, cache_dir=tmp_path / "cache") as f:
            assert f.fetch("s1", 1, 3) == "CG"
        assert any((tmp_path / "cache").glob("*.fai"))
        with Fasta(fa, cache_dir=tmp_path / "cache") as f:  # reuse
            assert f.fetch("s1", 0, 1) == "A"
    finally:
        d.chmod(0o755)


@pytest.mark.skipif(os.geteuid() == 0, reason="root ignores directory permissions")
def test_unwritable_cache_is_named(tmp_path):
    """An index cache that cannot be created is named as the cause, not the FASTA."""
    d, out = tmp_path / "ro", tmp_path / "out"
    d.mkdir()
    out.mkdir()
    (d / "g.fa").write_text(">s1\nACGT\n")
    d.chmod(0o555)
    out.chmod(0o555)
    try:
        with pytest.raises(OSError, match="cannot create the FASTA index cache .*out/work/fai"):
            Fasta(d / "g.fa", cache_dir=out / "work" / "fai")
    finally:
        d.chmod(0o755)
        out.chmod(0o755)


@pytest.mark.skipif(os.geteuid() == 0, reason="root ignores directory permissions")
def test_readonly_dir_without_cache_raises(tmp_path):
    d = tmp_path / "ro"
    d.mkdir()
    fa = d / "g.fa"
    fa.write_text(">s1\nACGT\n")
    d.chmod(0o555)
    try:
        with pytest.raises(ValueError, match="cache"):
            Fasta(fa)
    finally:
        d.chmod(0o755)


def test_plain_gzip_rejected(tmp_path):
    gz = tmp_path / "g.fa.gz"
    gz.write_bytes(gzip.compress(b">s1\nACGT\n"))
    with pytest.raises(ValueError, match="recompress with bgzip"):
        Fasta(gz, cache_dir=tmp_path / "cache")


def test_n_runs_and_revcomp():
    assert n_runs("ACNNNGTnA") == [(2, 5), (7, 8)]
    assert n_runs("ACGT") == []
    assert revcomp("ACGTN") == "NACGT"
    assert revcomp("acgRYk") == "mRYcgt"
    assert revcomp("ACGTRYKMBVDHNacgtrykmbvdhn") == "ndhbvkmryacgtNDHBVKMRYACGT"  # the full IUPAC table


@pytest.mark.skipif(os.geteuid() == 0, reason="root ignores directory permissions")
def test_cached_bgzip_index_follows_the_fasta(tmp_path):
    """A bgzip FASTA in a read-only directory is indexed (.fai and .gzi) in the cache, under a
    name keyed by the FASTA's size and mtime, so a rewritten FASTA is never read through its old index."""
    d = tmp_path / "ro"
    d.mkdir()
    gz, cache = d / "g.fa.gz", tmp_path / "cache"
    pysam.tabix_compress(str(_plain(tmp_path, ">s1\nACGTACGT\n")), str(gz))
    d.chmod(0o555)
    try:
        with Fasta(gz, cache_dir=cache) as f:
            assert f.lengths() == {"s1": 8} and f.fetch("s1", 2, 5) == "GTA"
        assert len(list(cache.glob("*.fai"))) == len(list(cache.glob("*.gzi"))) == 1
        d.chmod(0o755)
        gz.unlink()
        pysam.tabix_compress(str(_plain(tmp_path, ">s1\nAC\n>s2\nGGG\n")), str(gz))
        d.chmod(0o555)
        with Fasta(gz, cache_dir=cache) as f:
            assert f.lengths() == {"s1": 2, "s2": 3}
        assert len(list(cache.glob("*.fai"))) == 2
    finally:
        d.chmod(0o755)


def _plain(tmp_path, text):
    p = tmp_path / "plain.fa"
    p.write_text(text)
    return p


def test_index_build_is_atomic_and_errors_name_the_file(tmp_path):
    """The index is built under a temporary name and renamed (nothing else is left); a FASTA htslib cannot index
    is a ValueError naming it (not pysam's SamtoolsError, which workers would not catch); close() closes."""
    fa = tmp_path / "g.fa"
    fa.write_text(">s1\nACGT\n")
    f = Fasta(fa)
    assert sorted(p.name for p in tmp_path.iterdir()) == ["g.fa", "g.fa.fai"]
    f.close()
    with pytest.raises(ValueError, match="closed"):
        f.fetch("s1", 0, 2)
    bad = tmp_path / "bad.fa"
    bad.write_text(">s1\nACGT\nACGTACGT\nAC\n")  # lines of unequal length
    with pytest.raises(ValueError, match=r"bad\.fa: cannot index the FASTA"):
        Fasta(bad)
    assert not list(tmp_path.glob("*.tmp"))


def _write_three(tmp_path):
    text = ">s1 desc\nACGTNNNNacgt\nAC\n>s2\nGGGG\n"
    plain = tmp_path / "g.fa"
    plain.write_text(text)
    gz = tmp_path / "g.fa.gz"
    gz.write_bytes(gzip.compress(text.encode()))
    bgz = tmp_path / "g.fa.bgz"
    pysam.tabix_compress(str(plain), str(bgz))
    return plain, gz, bgz


def test_stream_fasta_plain_gzip_bgzip(tmp_path):
    """The same records from plain, gzip and bgzip files (sniffed from magic bytes, not the name); the name is the
    first word of the header; a record without sequence is empty."""
    want = [("s1", b"ACGTNNNNacgtAC"), ("s2", b"GGGG")]
    plain, gz, bgz = _write_three(tmp_path)
    renamed = tmp_path / "noext"
    renamed.write_bytes(gz.read_bytes())
    for p in (plain, gz, bgz, renamed):
        assert list(stream_fasta(p)) == want
    empty = tmp_path / "e.fa"
    empty.write_text(">a\n>b\nAC\n>c\n")
    assert list(stream_fasta(empty)) == [("a", b""), ("b", b"AC"), ("c", b"")]


def test_fasta_plain_gzip_names_bgzip(tmp_path):
    """A library FASTA in plain gzip streams, while indexed access asks for bgzip by name and names the file."""
    _, gz, bgz = _write_three(tmp_path)
    assert len(list(stream_fasta(gz))) == 2
    with pytest.raises(ValueError, match=r"g\.fa\.gz is gzip-compressed but not bgzip; recompress with bgzip"):
        Fasta(gz, cache_dir=tmp_path / "cache")
    with Fasta(bgz, cache_dir=tmp_path / "cache") as f:
        assert f.lengths() == {"s1": 14, "s2": 4}


def test_encode_and_revcomp_codes():
    codes = encode(b"ACGTacgtNnRy")
    assert codes.dtype == np.uint8 and codes.tolist() == [0, 1, 2, 3, 0, 1, 2, 3, 4, 4, 4, 4]
    assert revcomp_codes(encode(b"AACGTN")).tolist() == [4, 0, 1, 2, 3, 3]
    assert revcomp_codes(encode(b"")).tolist() == []
    assert revcomp_codes(codes).tolist() == encode(revcomp("ACGTacgtNnRy").encode()).tolist()


def test_stream_fasta_rejects_non_fasta(tmp_path):
    """Text before the first header is an error naming the file, not an empty stream."""
    bad = tmp_path / "notfasta.txt"
    bad.write_text("chr1\t100\t200\n")
    with pytest.raises(ValueError, match=r"notfasta\.txt: not a FASTA file"):
        list(stream_fasta(bad))
    gz = tmp_path / "bad.gz"
    gz.write_bytes(gzip.compress(b"ACGT\n>s\nAC\n"))
    with pytest.raises(ValueError, match=r"bad\.gz: not a FASTA file"):
        list(stream_fasta(gz))


def test_open_fasta_names_the_genome(tmp_path):
    """A FASTA pysam cannot open or index is one usage error naming the genome and the file; a good one opens."""
    bad = tmp_path / "bad.fa.gz"
    with gzip.open(bad, "wt") as fh:
        fh.write(">s1\nACGT\n")
    with pytest.raises(UsageError, match=r"^G1: cannot open FASTA .*bad\.fa\.gz: .*bgzip"):
        open_fasta(bad, "G1")
    good = tmp_path / "good.fa"
    good.write_text(">s1\nACGT\n")
    with open_fasta(good, "G2", tmp_path / "cache") as f:
        assert f.lengths() == {"s1": 4}


def test_encode_every_byte():
    """Every byte value: A C G T in either case are 0-3, anything else 4; the input is not changed."""
    raw = bytes(range(256))
    want = [{65: 0, 67: 1, 71: 2, 84: 3, 97: 0, 99: 1, 103: 2, 116: 3}.get(b, 4) for b in raw]
    got = encode(raw)
    assert got.dtype == np.uint8 and got.tolist() == want and raw == bytes(range(256))
    assert encode(b"").size == 0


def _line_reader(path):
    """The line-by-line reading stream_fasta must equal: each header line starts a record named by its first word;
    each other line, stripped of surrounding whitespace, is appended; non-blank text before the first header is an
    error."""
    opener = gzip.open if path.read_bytes()[:2] == b"\x1f\x8b" else open
    name, parts, out = None, [], []
    with opener(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if name is not None:
                    out.append((name, b"".join(parts)))
                fields = line[1:].split()
                name, parts = (fields[0].decode() if fields else ""), []
            elif name is not None:
                parts.append(line.strip())
            elif line.strip():
                raise ValueError("text before the first header")
    if name is not None:
        out.append((name, b"".join(parts)))
    return out


def test_stream_fasta_equals_line_by_line_reading(tmp_path, monkeypatch):
    """Records equal the line-by-line reading on random FASTAs: wrapped and unwrapped lines, CRLF ends, blank lines,
    spaces and tabs around and inside lines, empty records, headers without a name, no final newline, leading blank
    lines, gzip; read in blocks smaller than one record too."""
    from ltrtrace.io import seqio

    rng = np.random.default_rng(41)
    for case in range(40):
        recs = []
        for r in range(int(rng.integers(0, 12))):
            seq = "".join(rng.choice(list("ACGTNacgtn"), int(rng.integers(0, 300))))
            width = int(rng.integers(1, 90))
            lines = [seq[i:i + width] for i in range(0, len(seq), width)]
            if r % 4 == 1:
                lines = [f" {x}\t" for x in lines]
            if r % 5 == 2:
                lines.insert(len(lines) // 2, "")
            if r % 7 == 3 and lines:
                lines[0] = lines[0][:3] + " " + lines[0][3:]
            head = ">" if r % 9 == 4 else f">rec{r} some description"
            recs.append("\n".join([head, *lines]))
        text = ("\n\n" if case % 3 == 0 else "") + "\n".join(recs) + ("" if case % 2 else "\n")
        if case % 4 == 1:
            text = text.replace("\n", "\r\n")
        path = tmp_path / f"c{case}.fa"
        if case % 5 == 0:
            path = tmp_path / f"c{case}.fa.gz"
            path.write_bytes(gzip.compress(text.encode()))
        else:
            path.write_bytes(text.encode())
        assert list(seqio.stream_fasta(path)) == _line_reader(path), case
        with monkeypatch.context() as m:
            m.setattr(seqio, "FASTA_BLOCK", 7)
            assert list(seqio.stream_fasta(path)) == _line_reader(path), case


def test_read_fai(tmp_path):
    """A .fai as samtools writes it, or with just name and length, gives the lengths in file order; blank lines are
    skipped; any other line without a whole-number length fails, naming the file and the line."""
    fa = tmp_path / "g.fa"
    fa.write_text(">c2\nACGTA\n>c1\nAC\n")
    pysam.faidx(str(fa))
    assert list(read_fai(f"{fa}.fai").items()) == [("c2", 5), ("c1", 2)]
    two = tmp_path / "two.fai"
    two.write_text("s1\t100\n\ns2\t7\n")
    assert read_fai(two) == {"s1": 100, "s2": 7}
    for bad in ("s1\n", "s1\tlong\n", "s1\t-5\n"):
        (tmp_path / "bad.fai").write_text("s0\t3\n" + bad)
        with pytest.raises(ValueError, match=r"bad\.fai: line 2 is not a FASTA index line"):
            read_fai(tmp_path / "bad.fai")
