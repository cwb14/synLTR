"""Checkpoint helpers: cached file checksums, stage keys and done markers."""
import hashlib
import json
import logging
import os
import shutil
from pathlib import Path

import pysam
import pytest

from ltrtrace.base import checkpoint, constants
from ltrtrace.base.checkpoint import (
    bam_key,
    clear_done,
    code_fingerprint,
    done_outputs,
    file_checksum,
    is_done,
    mark_done,
    stage_key,
)


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def test_checkpoint_roundtrip(tmp_path, monkeypatch):
    monkeypatch.setattr(constants, "CHECKSUM_BLOCK", 3)  # stream the file in several blocks
    fa = tmp_path / "a.fa"
    fa.write_bytes(b">a\nACGTACGT\n")
    sha = file_checksum(fa, tmp_path / "sums.json")
    assert sha == _sha(fa.read_bytes())

    params = {"chunk": 100_000, "preset": "asm20", "nested": {"x": 1, "y": [1, 2]}}
    fp = "code"
    key = stage_key("anchor", {"A": sha, "B": "b" * 64}, params, fp)
    assert key == stage_key("anchor", {"B": "b" * 64, "A": sha}, dict(reversed(params.items())), fp)  # canonical
    assert key != stage_key("anchor", {"A": "0" * 64, "B": "b" * 64}, params, fp)  # an input checksum changed
    assert key != stage_key("anchor", {"A": sha, "B": "b" * 64}, {**params, "chunk": 50_000}, fp)  # a param
    assert key != stage_key("raw", {"A": sha, "B": "b" * 64}, params, fp)  # another stage, the same code
    assert key != stage_key("anchor", {"A": sha, "B": "b" * 64}, params, "other code")
    monkeypatch.setattr(checkpoint, "__version__", "99.0")
    assert key != stage_key("anchor", {"A": sha, "B": "b" * 64}, params, fp)  # another ltrtrace version
    monkeypatch.undo()

    out = tmp_path / "x.npz"
    marker = tmp_path / "x.npz.done.json"
    assert not is_done(marker, key)
    out.write_bytes(b"12345")
    mark_done(marker, key, [out])
    assert is_done(marker, key)
    assert not is_done(marker, "another key")
    assert sorted(p.name for p in tmp_path.iterdir()) == ["a.fa", "sums.json", "x.npz", "x.npz.done.json"]
    out.write_bytes(b"123456")  # an output rewritten after its marker
    assert not is_done(marker, key)
    out.write_bytes(b"12345")
    assert is_done(marker, key)
    clear_done(marker)
    assert not marker.exists() and not is_done(marker, key)
    clear_done(marker)  # no marker: nothing to do
    mark_done(marker, key, [out])
    out.unlink()
    assert not is_done(marker, key)


def test_done_outputs_lists_a_markers_files_whatever_its_key(tmp_path):
    """done_outputs: the files a marker lists, resolved against its folder, while each keeps its recorded size; None
    without a marker, with a malformed one, or once a listed file is gone or rewritten to another size."""
    a, b, mark = tmp_path / "a.npz", tmp_path / "sub" / "b.npz", tmp_path / "stage.done.json"
    assert done_outputs(mark) is None
    b.parent.mkdir()
    a.write_bytes(b"123")
    b.write_bytes(b"4567")
    mark_done(mark, "any key", [a, b])
    assert done_outputs(mark) == [a, b]
    b.write_bytes(b"45678")
    assert done_outputs(mark) is None and not is_done(mark, "any key")
    b.write_bytes(b"4567")
    a.unlink()
    assert done_outputs(mark) is None
    mark.write_text("[1, 2]")
    assert done_outputs(mark) is None and not is_done(mark, "any key")


def test_stage_key_needs_a_fingerprint():
    with pytest.raises(TypeError):
        stage_key("anchor", {}, {})


def test_mark_done_needs_its_outputs(tmp_path):
    marker = tmp_path / "x.done.json"
    with pytest.raises(FileNotFoundError, match="missing.npz"):
        mark_done(marker, "k", [tmp_path / "missing.npz"])
    assert list(tmp_path.iterdir()) == []


def test_marker_survives_moving_the_directory(tmp_path):
    work = tmp_path / "work"
    work.mkdir()
    (work / "x.npz").write_bytes(b"abc")
    mark_done(work / "x.npz.done.json", "k", [work / "x.npz"])
    moved = tmp_path / "moved"
    work.rename(moved)
    assert is_done(moved / "x.npz.done.json", "k")


def test_file_checksum_cache_follows_size_and_mtime(tmp_path, caplog):
    f, cache = tmp_path / "a.fa", tmp_path / "sums.json"
    f.write_bytes(b"AAAA")
    st = f.stat()
    first = file_checksum(f, cache)
    assert list(json.loads(cache.read_text())) == [str(f.resolve())]  # keyed by absolute path
    f.write_bytes(b"CCCC")  # same size, mtime restored: the cache cannot tell, the old sum comes back
    os.utime(f, ns=(st.st_atime_ns, st.st_mtime_ns))
    assert file_checksum(f, cache) == first
    os.utime(f, ns=(st.st_atime_ns, st.st_mtime_ns + 1_000_000))  # a new mtime: the file is read again
    assert file_checksum(f, cache) == _sha(b"CCCC")
    f.write_bytes(b"GGGGG")  # a new size, same mtime as the cached entry
    os.utime(f, ns=(st.st_atime_ns, st.st_mtime_ns + 1_000_000))
    assert file_checksum(f, cache) == _sha(b"GGGGG")
    assert file_checksum(tmp_path / "." / "a.fa", cache) == _sha(b"GGGGG")  # same file, another spelling
    assert len(json.loads(cache.read_text())) == 1

    cache.write_text("{not json")  # a damaged cache is rebuilt, with a warning
    with caplog.at_level(logging.WARNING, logger="ltrtrace.checkpoint"):
        assert file_checksum(f, cache) == _sha(b"GGGGG")
    assert "sums.json" in caplog.text
    assert json.loads(cache.read_text())[str(f.resolve())]["sha256"] == _sha(b"GGGGG")


def test_keyed_checksums_collect_what_the_stages_keyed_on(tmp_path, caplog):
    """keyed_checksums collects every checksum file_checksum returns while it is open, by absolute path, the cache's
    hits included, and nothing before or after; a file whose checksum changes inside it (it changed between two
    stages) keeps the latest, with one warning naming it."""
    f, g, h, cache = tmp_path / "a.fa", tmp_path / "b.fa", tmp_path / "c.fa", tmp_path / "sums.json"
    for path, data in ((f, b"AAAA"), (g, b"CCCC"), (h, b"TTTT")):
        path.write_bytes(data)
    file_checksum(g, cache)  # before: cached, not collected
    with checkpoint.keyed_checksums() as keyed:
        file_checksum(f, cache)
        file_checksum(tmp_path / "." / "b.fa", cache)  # a cache hit, by another spelling
        assert keyed == {str(f.resolve()): _sha(b"AAAA"), str(g.resolve()): _sha(b"CCCC")}
        f.write_bytes(b"GGGGG")
        with caplog.at_level(logging.WARNING, logger="ltrtrace.checkpoint"):
            for _ in range(2):
                file_checksum(f, cache)
        assert keyed[str(f.resolve())] == _sha(b"GGGGG")
    warned = [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]
    assert len(warned) == 1 and str(f.resolve()) in warned[0] and "changed during the run" in warned[0]
    file_checksum(h, cache)  # after: not collected
    assert set(keyed) == {str(f.resolve()), str(g.resolve())}


# ---------------------------------------------------------------- code fingerprints: a unit's import closure

_TREE = {  # a small package tree: relative path -> source
    "__init__.py": "VERSION = 1\n",
    "cli.py": "from ltrtrace.routes.assembly import pair\n",
    "base/__init__.py": "",
    "base/constants.py": "X = 1\n",
    "base/checkpoint.py": "from ltrtrace.base import constants\n",
    "seq/__init__.py": "",
    "seq/cigar.py": "import re\n",
    "seq/junction.py": "from ltrtrace.base.constants import X\n",
    "io/__init__.py": "",
    "io/seqio.py": "from ltrtrace.base.checkpoint import atomic_path\n",
    "routes/__init__.py": "",
    "routes/assembly/__init__.py": "",
    "routes/assembly/pair.py": ("from typing import TYPE_CHECKING\n\nfrom ltrtrace.io.seqio import Fasta\n\n"
                                "from .classify import call_row\n\nif TYPE_CHECKING:\n    from ltrtrace import cli\n"),
    "routes/assembly/classify.py": "from ltrtrace.seq import junction\n\n\ndef f():\n    import ltrtrace.seq.cigar\n",
    "routes/assembly/index.py": "from ltrtrace.base import checkpoint\n",
    "routes/library.py": "from ltrtrace.seq.cigar import parse\n",
    "calls/__init__.py": "",
    "calls/records.py": "Y = 2\n",
}
_CLOSURES = {  # unit -> the files of its import closure
    "routes.assembly.pair": {"routes/assembly/pair.py", "routes/assembly/classify.py", "io/seqio.py",
                             "base/checkpoint.py", "base/constants.py", "seq/junction.py", "seq/cigar.py"},
    "routes.assembly.index": {"routes/assembly/index.py", "base/checkpoint.py", "base/constants.py"},
    "routes.library": {"routes/library.py", "seq/cigar.py"},
    "calls": {"calls/__init__.py", "calls/records.py"},
}
_CLOSURES["routes.assembly"] = {"routes/assembly/__init__.py", *_CLOSURES["routes.assembly.pair"],
                                *_CLOSURES["routes.assembly.index"]}


@pytest.fixture
def pkg(tmp_path, monkeypatch):
    """The small package tree, read by code_fingerprint instead of the installed one."""
    root = tmp_path / "ltrtrace"
    for rel, text in _TREE.items():
        (root / rel).parent.mkdir(parents=True, exist_ok=True)
        (root / rel).write_text(text)
    monkeypatch.setattr(checkpoint, "_ROOT", root)
    checkpoint._fingerprint.cache_clear()
    return root


def _edit(path: Path, text: bytes = b"\n# an edit\n") -> None:
    """Append to a source file of the tree; fingerprints are cached per process (sources never change in a
    run), so the cache is cleared."""
    path.write_bytes(path.read_bytes() + text)
    checkpoint._fingerprint.cache_clear()


def _fps():
    return {u: code_fingerprint(u) for u in _CLOSURES}


def test_code_files_are_the_import_closure(pkg):
    """A unit (a package or a single module) and every ltrtrace module it imports, transitively: absolute,
    relative and in-function imports count, `from pkg import name` counts the module pkg.name, an import under
    `if TYPE_CHECKING:` does not, and neither does a package's __init__ that nothing imports from."""
    for unit, files in _CLOSURES.items():
        assert checkpoint.code_files(unit) == sorted(files), unit


@pytest.mark.parametrize("edited, changed", [
    ("seq/cigar.py", {"routes.assembly.pair", "routes.assembly", "routes.library"}),  # imported in a function
    ("base/checkpoint.py", {"routes.assembly.pair", "routes.assembly.index", "routes.assembly"}),  # via io.seqio
    ("base/constants.py", {"routes.assembly.pair", "routes.assembly.index", "routes.assembly"}),
    ("routes/assembly/index.py", {"routes.assembly.index", "routes.assembly"}),
    ("calls/records.py", {"calls"}),
    ("cli.py", set()),  # imported only under TYPE_CHECKING
    ("__init__.py", set()),  # nothing imports from the package's own __init__
])
def test_code_fingerprint_changes_with_its_import_closure(pkg, edited, changed):
    """An edit to a module a unit imports, directly or not, changes the unit's fingerprint; an edit to an
    unrelated module does not."""
    before = _fps()
    _edit(pkg / edited)
    after = _fps()
    assert {u for u in _CLOSURES if after[u] != before[u]} == changed


def test_code_fingerprint_names_an_unknown_unit(pkg):
    with pytest.raises(ValueError, match="no ltrtrace unit 'routes.nope'"):
        code_fingerprint("routes.nope")


def test_code_fingerprint_is_unchanged_by_a_copy(pkg, tmp_path, monkeypatch):
    """The fingerprint depends on the sources, not on where the package sits."""
    before = _fps()
    copy = tmp_path / "elsewhere" / "ltrtrace"
    shutil.copytree(pkg, copy)
    monkeypatch.setattr(checkpoint, "_ROOT", copy)
    checkpoint._fingerprint.cache_clear()
    assert _fps() == before


def test_code_fingerprint_of_the_installed_package():
    """The real sources are readable; the pair stage's closure holds the junction refinement and the call file
    writer, and not the command line or the elements stage."""
    assert len(code_fingerprint("base")) == 64
    files = checkpoint.code_files("routes.assembly.pair")
    assert {"routes/assembly/pair.py", "seq/junction.py", "seq/cigar.py", "calls/records.py",
            "io/seqio.py"} <= set(files)
    assert "cli.py" not in files and "routes/assembly/run.py" not in files


def test_bam_key_ignores_content_but_not_header(tmp_path):
    def write(path, name, seq):
        header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": name, "LN": 100}]}
        with pysam.AlignmentFile(path, "wb", header=header) as out:
            a = pysam.AlignedSegment(out.header)
            a.query_name, a.query_sequence, a.flag = "r", seq, 4
            out.write(a)

    bam = tmp_path / "a.bam"
    write(bam, "chr1", "ACGT")
    st = bam.stat()
    first = bam_key(bam)
    assert first == bam_key(bam)
    other = tmp_path / "b.bam"
    write(other, "chr1", "TTTT")  # other content, same size, same header
    os.utime(other, ns=(st.st_atime_ns, st.st_mtime_ns))
    assert other.stat().st_size == st.st_size and bam_key(other) == first
    os.utime(other, ns=(st.st_atime_ns, st.st_mtime_ns + 1))  # a new mtime changes the key
    assert bam_key(other) != first
    renamed = tmp_path / "c.bam"
    write(renamed, "chr2", "ACGT")  # another reference name, same size
    os.utime(renamed, ns=(st.st_atime_ns, st.st_mtime_ns))
    assert renamed.stat().st_size == st.st_size and bam_key(renamed) != first


def test_marker_records_meta(tmp_path):
    out, marker = tmp_path / "x.tsv", tmp_path / "x.tsv.done.json"
    out.write_text("x\n")
    mark_done(marker, "k", [out], meta={"elements": "e1", "code": "c1"})
    assert checkpoint.done_meta(marker) == {"elements": "e1", "code": "c1"}
    mark_done(marker, "k", [out])
    assert checkpoint.done_meta(marker) == {}
    assert checkpoint.done_meta(tmp_path / "none.done.json") is None
    marker.write_text("[1]")
    assert checkpoint.done_meta(marker) is None


def _exe(d, name, body):
    p = d / name
    p.write_text(f"#!/bin/sh\n{body}\n")
    p.chmod(0o755)
    return p


def test_tool_version(tmp_path, monkeypatch):
    """The first line a tool prints for its version flag, asked again when the binary changes (the cache is keyed
    by its size and modification time); an absent tool is None."""
    d = tmp_path / "bin"
    d.mkdir()
    lastz = _exe(d, "lastz", 'echo "lastz (version 9.9)"; echo "built with gcc"')
    monkeypatch.setenv("PATH", str(d))
    assert checkpoint.tool_version("lastz") == "lastz (version 9.9)"
    lastz.write_text('#!/bin/sh\necho "lastz (version 10.0.1)"\n')
    assert checkpoint.tool_version("lastz") == "lastz (version 10.0.1)"
    assert checkpoint.tool_version("minimap2") is None
    assert checkpoint.tool_version("tool-with-no-entry") is None  # not on PATH either


def test_checksum_cache_survives_concurrent_threads(tmp_path):
    """Threads of one process share a cold checksum cache; each writes through its own temp
    file (pid and thread id), so none fails on another's rename."""
    from concurrent.futures import ThreadPoolExecutor
    files = []
    for i in range(32):
        f = tmp_path / f"f{i}.fa"
        f.write_bytes(b"ACGT" * (i + 1))
        files.append(f)
    with ThreadPoolExecutor(8) as pool:
        sums = list(pool.map(lambda f: file_checksum(f, tmp_path / "sums.json"), files))
    assert sums == [_sha(f.read_bytes()) for f in files] and not list(tmp_path.glob(".*.tmp"))


def test_sha256_file_streams_the_file(tmp_path):
    """sha256_file is hashlib's sha256 of the whole file, read in CHECKSUM_BLOCK pieces (here more than one)."""
    data = os.urandom(constants.CHECKSUM_BLOCK + 12345)
    (tmp_path / "f").write_bytes(data)
    assert checkpoint.sha256_file(tmp_path / "f") == hashlib.sha256(data).hexdigest()
    (tmp_path / "e").write_bytes(b"")
    assert checkpoint.sha256_file(tmp_path / "e") == hashlib.sha256(b"").hexdigest()
