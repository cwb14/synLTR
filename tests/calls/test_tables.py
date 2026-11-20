import json
import subprocess
import sys
from datetime import UTC, datetime

from builders import ABC_LENGTHS as LENGTHS
from builders import ABC_ORDER as ORDER
from builders import abc_samples

from ltrtrace.calls.manifest import write_manifest
from ltrtrace.calls.tables import write_loci_tsv, write_pav_matrix, write_qc
from ltrtrace.calls.tracks import write_tracks
from ltrtrace.calls.vcf import write_vcf

SAMPLES = abc_samples()


def test_matrix_and_loci_tsv(tmp_path, loci):
    write_pav_matrix(loci, SAMPLES, tmp_path / "m.tsv")
    m = [ln.split("\t") for ln in (tmp_path / "m.tsv").read_text().splitlines()]
    assert m[0] == ["lid", "A", "B", "C"] and m[1] == ["L0000001", "P", "E", "D"]
    write_loci_tsv(loci, SAMPLES, tmp_path / "l.tsv")
    rows = [ln.split("\t") for ln in (tmp_path / "l.tsv").read_text().splitlines()]
    head = rows[0]
    assert head[:9] == ["lid", "family", "superfamily", "clade", "length", "ltr_len", "tsd", "discovery", "flags"]
    assert head[9:18] == ["A_eid", "A_state", "A_tier", "A_seqid", "A_start", "A_end", "A_strand", "A_k", "A_flags"]
    r = dict(zip(head, rows[1], strict=True))
    assert r["A_eid"] == "A:e1" and r["B_eid"] == "." and r["A_k"] == "0.02" and r["B_k"] == "."
    assert r["flags"] == "." and r["C_state"] == "D" and r["B_start"] == "50"
    r3 = dict(zip(head, rows[3], strict=True))
    assert r3["family"] == "." and r3["C_flags"] == "no_call" and r3["discovery"] == "B"
    r2 = dict(zip(head, rows[2], strict=True))
    assert r2["discovery"] == "B,C"


def test_qc(tmp_path):
    write_qc([("all", "n_loci", 4), ("all", "frac", 0.123456789)], tmp_path / "qc.tsv")
    assert (tmp_path / "qc.tsv").read_text().splitlines() == ["scope\tmetric\tvalue", "all\tn_loci\t4",
                                                              "all\tfrac\t0.123457"]


def test_manifest(tmp_path):
    t0 = datetime(2026, 1, 1, tzinfo=UTC)
    write_manifest(tmp_path / "manifest.json", ["ltrtrace", "run"], {"threads": 4}, SAMPLES,
                   {"A": {"fasta": "abc"}}, t0, t0, units=["calls"], route="assembly")
    m = json.loads((tmp_path / "manifest.json").read_text())
    assert m["tool"] == "ltrtrace" and m["command"] == ["ltrtrace", "run"] and m["params"] == {"threads": 4}
    assert m["route"] == "assembly"
    assert m["samples"][0]["label"] == "A" and m["samples"][0]["checksums"] == {"fasta": "abc"}
    assert m["samples"][1]["checksums"] == {} and m["samples"][1]["library"] is None
    assert m["started"].startswith("2026-01-01") and "pysam" in m["dependencies"] and "host" in m
    assert list(m["code"]) == ["calls"] and len(m["code"]["calls"]) == 64


_FAILING_WRITE = """
import resource, signal, sys
from datetime import UTC, datetime
from pathlib import Path

from ltrtrace.calls.manifest import write_manifest
from ltrtrace.io.samples import Sample

t0 = datetime(2026, 1, 1, tzinfo=UTC)
signal.signal(signal.SIGXFSZ, signal.SIG_IGN)  # a write past the limit fails with EFBIG instead of a signal
resource.setrlimit(resource.RLIMIT_FSIZE, (64, 64))
try:
    write_manifest(Path(sys.argv[1]), ["ltrtrace", "run"], {}, [Sample("A", None, None, None, None)], {}, t0, t0)
except OSError as exc:
    sys.exit(f"the write failed: {exc}")
"""


def test_manifest_write_that_fails_midway_keeps_the_earlier_manifest(tmp_path):
    """A manifest write that fails midway (a child process whose file size limit cuts every file at 64 bytes) leaves
    the earlier manifest.json as it was and no temporary file: the manifest is written under a temporary name and
    renamed into place."""
    path = tmp_path / "manifest.json"
    path.write_text('{"the earlier": "run"}\n')
    res = subprocess.run([sys.executable, "-c", _FAILING_WRITE, str(path)], capture_output=True, text=True,
                         check=False)
    assert res.returncode == 1 and "the write failed" in res.stderr, res.stderr
    assert path.read_text() == '{"the earlier": "run"}\n'
    assert sorted(p.name for p in tmp_path.iterdir()) == ["manifest.json"]


def test_deterministic(tmp_path, loci, fa):
    outs = []
    for d in ("x", "y"):
        o = tmp_path / d
        o.mkdir()
        write_loci_tsv(loci, SAMPLES, o / "loci.tsv")
        write_pav_matrix(loci, SAMPLES, o / "m.tsv")
        write_qc([("all", "n", 1)], o / "qc.tsv")
        write_vcf(loci, SAMPLES, "A", fa, o / "pav.vcf.gz")
        write_tracks(loci, SAMPLES, o, ORDER, LENGTHS)
        outs.append(o)
    names = sorted(p.name for p in outs[0].iterdir())
    assert "pav.vcf.gz" in names and "A.pav.gff3" in names
    for n in names:
        assert (outs[0] / n).read_bytes() == (outs[1] / n).read_bytes(), n
