"""Catalog loader (analyses.geneflow.catalog): toy run directories written through the route's writers and merge
(catalog_toy.py) read back as arrays, merged states beside each discoverer's own calls; fail-fast errors on a
missing or mismatched run; age bins; the toy writing only what the route writes; the simulated toy's determinism;
and an analysis edit that leaves every genotyping stage key unchanged."""
from __future__ import annotations

import gzip
import logging
import os
import shutil
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
from analyses.geneflow.catalog_toy import loc, simulate_run, write_run

import synltr
from synltr.analyses.geneflow.catalog import (
    FLAG_BIT,
    STATE_CODE,
    TIER_CODE,
    AnalysisError,
    age_bin,
    cluster_keys,
    load_catalog,
)
from synltr.base import checkpoint


def test_load_roundtrip(tmp_path, toy_fastas):
    run = write_run(tmp_path, [
        loc(10_001, "R:P* S:E D:P* O:E", k=0.0, calls={("D", "S"): "U"}),
        loc(20_001, "R:U S:P* D:P* O:P", k=0.02),
        loc(5_001, "R:P* S:D D:E O:E", seqid="chr2", k=None, length=6000, flags={"S": ["inherit"]}),
        loc(9_001, "R:P* S:P D:G O:U", seqid="chr2", parent=2)], fasta_dir=toy_fastas)
    cat = run.load()
    assert cat.labels == ("R", "S", "D", "O") and cat.outgroups == ("O",) and cat.ingroup == ("R", "S", "D")
    assert cat.sheet == run.sheet and cat.tsd[0] == "ACGTA"
    assert cat.cls("S").tolist() == [0, 1, -1, 1] and cat.cls("D").tolist() == [1, 1, 0, -1]
    assert cat.member[:, cat.g("D")].tolist() == [True, True, False, False]
    assert cat.k[0, cat.g("D")] == 0.0 and np.isnan(cat.k[0, cat.g("S")])         # K only for members
    assert cat.flags[2, cat.g("S")] & FLAG_BIT["inherit"]
    assert cat.top["R"][cat.eid[3, cat.g("R")]] == cat.eid[2, cat.g("R")]         # nested in locus 2's element
    assert cat.lengths["R"].tolist() == [1_000_000] * 4
    # the discoverer's own calls: D's call into S is U although the merged S is E (R's call)
    assert cat.pcls("D", "S")[0] == -1 and cat.pcls("R", "S")[0] == 0 and cat.pcls("D", "D")[0] == 1
    assert cat.ptier["D"][0, cat.g("R")] == TIER_CODE["high"]                     # a member's call tier, never "annot"
    assert cat.pcls("S", "R")[0] == -1                                            # S has no member at locus 0: no call
    gstar = np.array([cat.g("D"), cat.g("D"), -1, -1], dtype=np.int16)
    assert cat.pcls_by(gstar, "S").tolist() == [-1, 1, -1, -1]
    seq, mid = cat.proj_site(gstar, "D")
    assert seq.tolist()[:2] == [0, 0] and mid[0] == pytest.approx(10_001 + 149 / 2) and seq[2] == -1


def test_missing_or_mismatched_run(tmp_path, toy_fastas):
    run = write_run(tmp_path, [loc(10_001, "R:P* S:E D:P* O:E")], fasta_dir=toy_fastas)
    loci = run.outdir / "loci.tsv"
    text = loci.read_text()
    loci.rename(run.outdir / "loci.bak")
    with pytest.raises(AnalysisError, match="synltr genotype"):
        run.load()
    loci.write_text(text)
    with pytest.raises(AnalysisError, match="X"):
        load_catalog(run.outdir, [replace(s, label="X") if s.label == "O" else s for s in run.samples], run.sheet)
    loci.write_text(text.replace("\tE\t", "\tQ\t", 1))
    with pytest.raises(AnalysisError, match="line 2"):
        run.load()
    loci.write_text(text)
    pair = run.outdir / "pairs" / "D__S.calls.tsv.gz"
    with gzip.open(pair, "wt") as fh:
        fh.write("eid\tstate\n")
    with pytest.raises(AnalysisError, match="D__S.*CALL_COLUMNS"):
        run.load()
    pair.unlink()
    with pytest.raises(AnalysisError, match="D__S"):
        run.load()


def test_own_calls_strands_clusters_and_mismatches(tmp_path, toy_fastas):
    run = write_run(tmp_path, [
        loc(10_001, "R:P* S:P D:E O:E"),                              # S inverted here: R's window into S is "-"
        loc(20_001, "R:P* S:P* D:E O:E", strand={"S": "."}),          # S's own strand unknown (a BED annotation)
        loc(60_001, "R:P* S:E D:P* O:E", length=6000),                # a host in R ...
        loc(62_001, "R:P* S:E D:E O:E", parent=2)],                   # ... and an element nested in it
        inversion=("S", "chr1", 1, 50_000), fasta_dir=toy_fastas)
    cat = run.load()
    r, s, d = cat.g("R"), cat.g("S"), cat.g("D")
    assert cat.strand[0, s] == -1
    assert cat.flags[0, s] & FLAG_BIT["inv"]                                    # decided by a reversed window
    assert (cat.pstart["R"][0, s], cat.pend["R"][0, s]) == (1 + 50_000 - 10_150, 1 + 50_000 - 10_001)
    assert cat.pbstrand["S"][1, r] == -1 and cat.pbstrand["R"][0, s] == -1      # S window known, own strand not
    assert cat.pbstrand["R"][0, r] == 1 and cat.pbstrand["R"][0, d] == 1        # own column; an E call is sited
    assert cat.pbstrand["D"][0, r] == 0                                         # D has no member here: no call
    keys = cluster_keys(cat, np.full(4, r, dtype=np.int16), np.array([False, False, True, True]))
    assert keys.tolist() == [("R", cat.eid[2, r])] * 2                          # host and nested: one cluster
    with pytest.raises(ValueError, match="no member element"):
        cluster_keys(cat, np.full(4, d, dtype=np.int16), np.array([True, False, False, False]))
    with pytest.raises(ValueError, match="boolean mask"):                       # an index array, not a mask
        cluster_keys(cat, np.full(4, r, dtype=np.int16), np.array([2, 3]))
    pair = run.outdir / "pairs" / "R__D.calls.tsv.gz"
    lines = gzip.decompress(pair.read_bytes()).decode().splitlines(keepends=True)
    extra = lines[1].replace(cat.eid[0, r], "R|chr1:1-150", 1)    # no locus's first member element
    pair.write_bytes(gzip.compress("".join([*lines, extra]).encode()))
    assert run.load().proj["R"][0, d] == STATE_CODE["E"]                        # skipped, not an error
    pair.write_bytes(gzip.compress("".join(lines[:1] + lines[2:]).encode()))    # a member element's call lost
    with pytest.raises(AnalysisError, match="R__D.*1 member element.*no call"):
        run.load()
    pair.write_bytes(gzip.compress("".join(lines).encode())[:-20])
    with pytest.raises(AnalysisError, match="R__D.*corrupt"):
        run.load()


def _pair(run, a, b):
    """A's pair file into B as {eid: row}."""
    path = run.outdir / "pairs" / f"{a}__{b}.calls.tsv.gz"
    lines = gzip.decompress(path.read_bytes()).decode().splitlines()
    head = lines[0].split("\t")
    return {row["eid"]: row for row in (dict(zip(head, line.split("\t"), strict=True)) for line in lines[1:])}


def test_toy_writes_what_the_engine_writes(tmp_path, toy_fastas):
    run = write_run(tmp_path, [
        loc(10_001, "R:P* S:O D:D O:T", evidence={"S": (2, -1, 300)}),
        loc(20_001, "R:P* S:S D:E O:G"),
        loc(30_001, "R:P* S:P D:U O:E")], inversion=("S", "chr1", 25_001, 40_000), fasta_dir=toy_fastas)
    cat = run.load()
    r, s, d, o = (cat.g(x) for x in "RSDO")
    calls = {b: _pair(run, "R", b) for b in "SDO"}
    first, second, third = (cat.eid[i, r] for i in range(3))
    assert [calls["S"][first][c] for c in ("state", "tier", "b_start", "b_end", "b_strand")] == \
        ["O", "medium", "10001", "10300", "+"]                                 # O: B's inserted sequence
    assert [calls["D"][first][c] for c in ("state", "tier", "b_start", "b_end")] == \
        ["D", "medium", "10000", "10001"]                                      # D: the junction
    assert calls["O"][first]["tier"] == "medium"                               # T and O are never high
    assert [calls["S"][second][c] for c in ("tier", "flags", "b_start", "b_end")] == \
        ["high", "solo_tsd", "20001", "20040"]                                 # S: its solo LTR
    assert [calls["D"][second][c] for c in ("tier", "flags", "b_start", "b_end")] == \
        ["high", "tsd_match", "19996", "20000"]                                # E: its TSD copy
    assert [calls["O"][second][c] for c in ("tier", "b_seqid", "wid")] == [".", ".", "1"]
    assert [calls["S"][third][c] for c in ("b_start", "b_end", "b_strand", "flags")] == \
        [str(65_001 - 30_150), str(65_001 - 30_001), "-", "inv"]               # mirrored: a reversed window
    assert [calls["D"][third][c] for c in ("tier", "reason", "wid", "flags", "b_strand")] == \
        [".", "no_anchor", "2", ".", "."]                                      # U: a window id, no inv
    every = [row for rows in calls.values() for row in rows.values()]
    assert all((row["b_strand"] in ("+", "-")) == (row["b_seqid"] != ".") for row in every)
    assert all(("inv" in row["flags"].split(",")) == (row["b_strand"] == "-") for row in every
               if row["b_seqid"] != ".")                                       # a sited call's window
    # loci.tsv is the route's merge of these calls; R is the only member, so each state is R's call
    for x in (s, d, o):
        for merged, own in ((cat.state, cat.proj), (cat.tier, cat.ptier), (cat.start, cat.pstart)):
            assert (merged[:, x] == own["R"][:, x]).all()
    assert cat.flags[2, s] & FLAG_BIT["inv"] and cat.strand[2, s] == -1
    for bad, match in (([loc(10_001, "R:P* S:P D:E O:E", strand={"S": "."})], "only for a member"),
                       ([loc(49_990, "R:P* S:P D:E O:E")], "leaves the sequence"),   # mirrored below 1
                       ([loc(10_001, "R:P* S:E D:P* O:E", calls={("R", "D"): "U", ("D", "R"): "U"})],
                        "splits"),
                       ([loc(10_001, "R:P* S:E D:P* O:E", calls={("R", "S"): "P"})], "merge to"),
                       ([loc(10_001, "R:P* S:E D:E O:E", flags={"S": ["inv"]})], "set by the engine")):
        with pytest.raises(ValueError, match=match):
            write_run(tmp_path / "bad", bad, inversion=("S", "chr1", 1, 50_000), fasta_dir=toy_fastas)


@pytest.mark.skipif(os.name != "posix" or os.geteuid() == 0, reason="directory permissions do not bind root")
def test_read_only_run(tmp_path, toy_fastas, caplog):
    """A run directory that cannot be written is read all the same: the loader writes nothing into it."""
    run = write_run(tmp_path, [loc(10_001, "R:P* S:E D:P* O:E")], fasta_dir=toy_fastas)
    run.load()  # the FASTA indexes exist from here on
    for d in (run.outdir / "work", run.outdir):
        d.chmod(0o555)
    try:
        with caplog.at_level(logging.WARNING, logger="synltr"):
            cat = run.load()
    finally:
        for d in (run.outdir, run.outdir / "work"):
            d.chmod(0o755)
    assert len(cat.lid) == 1 and cat.sheet == run.sheet and not caplog.records


def test_age_bin():
    assert age_bin(np.array([0.0, 0.004, 0.01, 0.05, np.nan])).tolist() == [0, 1, 2, 2, -1]


def test_simulate_run_deterministic(tmp_path, toy_fastas):
    a = simulate_run(tmp_path / "a", seed=3, fasta_dir=toy_fastas)
    b = simulate_run(tmp_path / "b", seed=3, fasta_dir=toy_fastas)
    assert (a.outdir / "loci.tsv").read_bytes() == (b.outdir / "loci.tsv").read_bytes()
    cat = a.load()
    assert len(cat.lid) > 4000 and cat.member.any(axis=1).all()     # every locus has a discoverer
    with pytest.raises(ValueError, match="spacing"):                # 998_000 // 4000 = 249 bp < 2 x 150 bp
        simulate_run(tmp_path / "c", seed=3, n_per_seq=4000, fasta_dir=toy_fastas)


ENGINE_UNITS = ("io.annotation.stage", "routes.assembly.index", "routes.assembly.anchor", "routes.assembly.pair")


def test_analysis_edit_keeps_engine_keys(tmp_path, monkeypatch):
    """An edit to an analysis module changes the identity alignment's code fingerprint and no genotyping stage's:
    the routes never import the analyses (tests/test_layers.py)."""
    copy = tmp_path / "synltr"
    shutil.copytree(Path(synltr.__file__).parent, copy, ignore=shutil.ignore_patterns("__pycache__"))
    monkeypatch.setattr(checkpoint, "_ROOT", copy)
    checkpoint._fingerprint.cache_clear()
    try:
        def prints():
            return ({u: checkpoint.code_fingerprint(u) for u in ENGINE_UNITS},
                    checkpoint.code_fingerprint("analyses.geneflow.identity"))
        e0, a0 = prints()
        with open(copy / "analyses" / "geneflow" / "catalog.py", "a") as fh:
            fh.write("\n# an edit to an analysis module\n")
        checkpoint._fingerprint.cache_clear()
        e1, a1 = prints()
    finally:
        checkpoint._fingerprint.cache_clear()
    assert a1 != a0 and e1 == e0                                    # every genotyping stage key unchanged
