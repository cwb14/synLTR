"""ltrtrace source on a planted genealogy with known answers."""
import csv
import json
import re
import shutil

import numpy as np
import pytest

from ltrtrace import cli
from ltrtrace.analyses import source
from ltrtrace.base.log import UsageError
from ltrtrace.io.annotation import load_annotation
from ltrtrace.io.seqio import Fasta

pytestmark = pytest.mark.skipif(shutil.which("minimap2") is None, reason="minimap2 (the source screen) not on PATH")

BASES = np.frombuffer(b"ACGT", dtype=np.uint8)
RC = bytes.maketrans(b"ACGT", b"TGCA")


def _mut(rng, a, n):
    a = a.copy()
    pos = rng.choice(len(a), n, replace=False)
    a[pos] = (a[pos] + rng.integers(1, 4, n)) % 4
    return a, set(pos.tolist())


@pytest.fixture(scope="module")
def files(tmp_path_factory):
    d = tmp_path_factory.mktemp("src")
    rng = np.random.default_rng(3)
    ltr0, int0 = rng.integers(0, 4, 500).astype(np.uint8), rng.integers(0, 4, 5000).astype(np.uint8)
    cps = {}
    for i in range(6):                                   # background fixed copies (consensus = ancestor)
        cps[f"bg{i}"] = (ltr0, ltr0, _mut(rng, int0, 12)[0], "F/F")
    M_int, _ = _mut(rng, int0, 5)                  # master: 5 private mutations
    M_l3, _ = _mut(rng, ltr0, 3)                        # master LTRs differ at 3 sites (older)
    cps["master"] = (ltr0, M_l3, M_int, "F/F")
    cps["child"] = (M_l3, M_l3, _mut(rng, M_int, 2)[0], "F/E")      # copied from master after its mutations
    P_int, _ = _mut(rng, int0, 3)                        # lost parent with 3 private mutations
    cps["sibA"] = (ltr0, ltr0, _mut(rng, P_int, 4)[0], "F/F")
    cps["sibB"] = (ltr0, ltr0, _mut(rng, P_int, 4)[0], "F/E")
    T_int, _ = _mut(rng, int0, 6)
    cps["twin1"], cps["twin2"] = (ltr0, ltr0, T_int, "F/F"), (ltr0, ltr0, T_int.copy(), "F/F")
    cps["twin_child"] = (ltr0, ltr0, _mut(rng, T_int, 1)[0], "F/E")
    R1_int, _ = _mut(rng, int0, 15)                      # two fixed parents and a template-switch child
    R2_int, _ = _mut(rng, int0, 15)
    cps["rp1"], cps["rp2"] = (ltr0, ltr0, R1_int, "F/F"), (ltr0, ltr0, R2_int, "F/F")
    cps["rchild"] = (ltr0, ltr0, np.concatenate([R1_int[:2500], R2_int[2500:]]), "F/E")
    H_int, _ = _mut(rng, int0, 7)                        # fixed host with a new copy of master nested inside
    inner_el = np.concatenate([M_l3, _mut(rng, M_int, 1)[0], M_l3])
    t_in = H_int[2000:2005]
    cps["host"] = (ltr0, ltr0, np.concatenate([H_int[:2005], inner_el, t_in, H_int[2005:]]), "F/F")
    D_int, _ = _mut(rng, int0, 8)
    cps["dup_src"] = (ltr0, ltr0, D_int, "F/F")
    cps["dup_copy"] = (ltr0, ltr0, D_int.copy(), "F/E")
    flank_dup = rng.integers(0, 4, 800).astype(np.uint8)
    inv = np.random.default_rng(99)                      # the inverted duplicate has its own random stream
    V_int, _ = _mut(inv, int0, 9)
    cps["inv_src"] = (ltr0, ltr0, V_int, "F/F")
    cps["inv_copy"] = (ltr0, ltr0, V_int.copy(), "F/E")
    tsd_inv = np.array([0, 1, 1, 2, 3], np.uint8)        # ACCGT: not its own reverse complement
    seq, gff, gts = [], [], []
    pos = 0
    tsd_dup = rng.integers(0, 4, 5).astype(np.uint8)
    for name, (l5, l3, internal, gt) in cps.items():
        dup = name.startswith("dup")
        r = inv if name.startswith("inv") else rng
        sp = r.integers(0, 4, 4000).astype(np.uint8)
        fl = flank_dup if dup else r.integers(0, 4, 800).astype(np.uint8)
        tsd = tsd_dup if dup else r.integers(0, 4, 5).astype(np.uint8)
        if name.startswith("inv"):                       # same site, opposite orientation, unrelated flanks
            tsd = tsd_inv if name == "inv_src" else (3 - tsd_inv[::-1]).astype(np.uint8)
        el = np.concatenate([l5, internal, l3])
        s = pos + len(sp) + len(fl) + 5
        e = s + len(el)
        seq += [sp, fl, tsd, el, tsd, fl[::-1].copy()]
        pos = e + 5 + len(fl)
        fam = "famV" if name.startswith("inv") else "famX"      # the inverted pair is its own family
        gff.append(f"chr1\tt\tLTR_retrotransposon\t{s + 1}\t{e}\t.\t+\t.\tID={name};family={fam};clade=cl")
        gff.append(f"chr1\tt\tlong_terminal_repeat\t{s + 1}\t{s + 500}\t.\t+\t.\tID={name}.l;Parent={name}")
        gff.append(f"chr1\tt\tlong_terminal_repeat\t{e - 499}\t{e}\t.\t+\t.\tID={name}.r;Parent={name}")
        gts.append((name, "s1", gt, "99", "PASS"))
        if name == "host":                               # annotate the nested inner element (new)
            i0 = s + 500 + 2005
            i1 = i0 + len(inner_el)
            gff.append(f"chr1\tt\tLTR_retrotransposon\t{i0 + 1}\t{i1}\t.\t+\t.\tID=inner;family=famX;clade=cl")
            gff.append(f"chr1\tt\tlong_terminal_repeat\t{i0 + 1}\t{i0 + 500}\t.\t+\t.\tID=inner.l;Parent=inner")
            gff.append(f"chr1\tt\tlong_terminal_repeat\t{i1 - 499}\t{i1}\t.\t+\t.\tID=inner.r;Parent=inner")
            gts.append(("inner", "s1", "F/E", "99", "PASS"))
    seq.append(rng.integers(0, 4, 3000).astype(np.uint8))
    s = BASES[np.concatenate(seq)].tobytes().decode()
    (d / "g.fa").write_text(">chr1\n" + "\n".join(s[i:i + 80] for i in range(0, len(s), 80)) + "\n")
    (d / "g.gff3").write_text("##gff-version 3\n" + "\n".join(gff) + "\n")
    return {"dir": d, "fa": d / "g.fa", "gff": d / "g.gff3", "gt": gts}


HEADER = ("element_id", "sample", "genotype", "gq", "flags")


def _run_dir(files, tmp, gff=None, gt=None, columns=HEADER, name="run"):
    """A finished reads run: its manifest (the reference row) and reads.genotypes.tsv."""
    run = tmp / name
    run.mkdir(parents=True, exist_ok=True)
    gt = files["gt"] if gt is None else gt
    manifest = {"samples": [{"label": "ref", "fasta": str(files["fa"]), "annotation": str(gff or files["gff"]),
                             "role": "reference"}, {"label": "s1", "fasta": None, "annotation": None,
                                                    "role": None}]}
    (run / "manifest.json").write_text(json.dumps(manifest))
    (run / "reads.genotypes.tsv").write_text(
        "\n".join("\t".join(map(str, r)) for r in [columns, *gt]) + "\n")
    return run


def _source(files, out, gff=None, gt=None, columns=HEADER, threads=2):
    run = _run_dir(files, out.parent, gff, gt, columns, name=out.name + "_run")
    assert source.run(run, out, threads) == 0
    with open(out / "reads.sources.tsv") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@pytest.fixture(scope="module")
def genome(files, tmp_path_factory):
    return _source(files, tmp_path_factory.mktemp("out") / "out")


def _by(rows, new):
    return [r for r in rows if r["new_id"] == new]


def test_parent_found_with_shared_mutations(genome):
    (r,) = _by(genome, "child")
    assert r["source_id"] == "master"
    assert int(float(r["shared_derived"])) == 5 and int(float(r["source_private"])) == 0
    assert r["source_older"] == "True" and r["flags"] == "PASS"


def test_sibling_has_source_private_mutations(genome):
    (r,) = _by(genome, "sibB")
    assert r["source_id"] == "sibA"
    assert int(float(r["shared_derived"])) == 3 and int(float(r["source_private"])) == 4


def test_ties_share_weight(genome):
    rs = _by(genome, "twin_child")
    assert {r["source_id"] for r in rs} == {"twin1", "twin2"}
    assert all(float(r["weight"]) == 0.5 and "tied" in r["flags"] for r in rs)


def test_segmental_duplication_flagged(genome):
    (r,) = _by(genome, "dup_copy")
    assert r["source_id"] == "dup_src" and "segdup" in r["flags"] and r["same_tsd"] == "True"


def test_inverted_segdup_flagged(genome):
    """The same TSD read from the other strand (reverse complement) is the same site: the copy is a duplicate, though
    its flanks are unrelated and its TSD is a different string."""
    (r,) = _by(genome, "inv_copy")
    assert r["source_id"] == "inv_src"
    assert r["same_tsd"] == "True" and "segdup" in r["flags"]
    assert float(r["flank_identity"]) < 0.9                # flagged by the TSD alone


def test_two_parent_mosaic(genome):
    rs = _by(genome, "rchild")
    assert rs and all("recombinant" in r["flags"] for r in rs)
    r = rs[0]
    assert {r["mosaic_left"], r["mosaic_right"]} == {"rp1", "rp2"}
    assert 1500 <= float(r["mosaic_breakpoint"]) <= 3500


def test_single_parent_not_mosaic(genome):
    for new in ("child", "sibB", "twin_child"):
        assert all("recombinant" not in r["flags"] for r in _by(genome, new)), new


def test_nested_host_is_never_the_source(genome):
    rs = _by(genome, "inner")
    assert rs and all(r["source_id"] != "host" for r in rs)
    assert any(r["source_id"] == "master" for r in rs)


def test_tsd_is_compared_by_orientation():
    assert source.same_site("ACCGT", "ACGGT") is True
    assert source.same_site("ACCGT", "ACCGT") is True
    assert source.same_site("ACCGT", "ACCGA") is False
    assert source.same_site("", "") is False and source.same_site("ACCGT", "") is False


# ---- inputs and caching ------------------------------------------------------------------------------------------
def test_no_new_copies_gives_empty_tables(files, tmp_path):
    """No new copy is a result, not an error: header-only tables."""
    gt = [(e, s, "F/F", q, f) for e, s, _, q, f in files["gt"]]
    out = tmp_path / "none"
    assert _source(files, out, gt=gt) == []
    for suf in (".sources.tsv", ".copies.tsv", ".locality.tsv"):
        assert len((out / f"reads{suf}").read_text().splitlines()[0].split("\t")) > 3, suf


def test_genotypes_of_another_annotation_are_rejected(files, tmp_path):
    gt = [("x_" + e, *rest) for e, *rest in files["gt"]]
    with pytest.raises(UsageError, match="annotation used for genotyping"):
        _source(files, tmp_path / "o", gt=gt)


def test_moved_coordinates_are_rejected(files, tmp_path):
    table = load_annotation(str(files["gff"]), "g", Fasta(files["fa"], cache_dir=tmp_path))
    els = {e.name: e for e in table}
    gt = [(e, els[e].seqid, els[e].start + (100 if e == "child" else 0), els[e].end + (100 if e == "child" else 0),
           *rest) for e, *rest in files["gt"]]
    with pytest.raises(UsageError, match="other coordinates"):
        _source(files, tmp_path / "m", gt=gt, columns=("element_id", "chrom", "start", "end", *HEADER[1:]))


def test_coordinates_of_the_run_are_accepted(files, tmp_path):
    table = load_annotation(str(files["gff"]), "g", Fasta(files["fa"], cache_dir=tmp_path))
    els = {e.name: e for e in table}
    gt = [(e, els[e].seqid, els[e].start, els[e].end, *rest) for e, *rest in files["gt"]]
    rows = _source(files, tmp_path / "ok", gt=gt, columns=("element_id", "chrom", "start", "end", *HEADER[1:]))
    assert any(r["new_id"] == "child" for r in rows)


def test_no_fixed_relative_in_any_family(files, tmp_path):
    """Every new copy alone in its family: nothing to assign, and no crash."""
    new = {r[0] for r in files["gt"] if r[2] == "F/E"}
    lines = []
    for ln in files["gff"].read_text().splitlines(keepends=True):
        eid = ln.split("ID=")[1].split(";")[0] if "\tLTR_retrotransposon\t" in ln else None
        lines.append(re.sub(r"family=\w+", f"family=fam_{eid}", ln) if eid in new else ln)
    g = tmp_path / "alone.gff3"
    g.write_text("".join(lines))
    assert _source(files, tmp_path / "a", gff=g) == []


def test_work_cache_follows_family_labels(files, tmp_path):
    """A rerun in the same folder after relabeling a family must not reuse the earlier alignments."""
    out = tmp_path / "c"
    assert any(r["new_id"] == "child" for r in _source(files, out))
    g = tmp_path / "relabel.gff3"
    g.write_text(files["gff"].read_text().replace("ID=child;family=famX", "ID=child;family=famY"))
    rerun = _source(files, out, gff=g)
    assert rerun == _source(files, tmp_path / "fresh", gff=g)
    assert not any(r["new_id"] == "child" for r in rerun)


def test_source_cache_keys_code(files, tmp_path, monkeypatch):
    """The cached alignments are reused by the same code and parameters, and recomputed when the code that makes
    them changes or when MIN_COV or TOP_K do."""
    calls = []
    real = source.screen
    monkeypatch.setattr(source, "screen", lambda *a, **k: calls.append(1) or real(*a, **k))
    out = tmp_path / "k"

    def again():
        calls.clear()
        rows = _source(files, out, threads=1)
        return len(calls), rows

    n0, rows0 = again()
    assert n0 > 0
    assert again() == (0, rows0)                                       # nothing changed: reused
    monkeypatch.setattr(source, "code_fingerprint", lambda unit: "edited")
    n1, rows1 = again()
    assert n1 == n0 and rows1 == rows0                                 # other code: recomputed, same result here
    assert again()[0] == 0
    for name, value in (("TOP_K", source.TOP_K - 1), ("MIN_COV", source.MIN_COV - 0.1)):
        monkeypatch.setattr(source, name, value)
        assert again()[0] == n0, name
        assert again()[0] == 0, name


def test_source_needs_a_reads_run(files, tmp_path):
    """source reads RUN's manifest.json through the package's one manifest reader (calls.manifest.read_manifest),
    with its own messages for a folder without one and a manifest without a reference row; a manifest that cannot
    be read, or is no manifest, is the reader's usage error naming the file."""
    with pytest.raises(UsageError, match=re.escape(f"{tmp_path} is not a finished read genotyping run: no "
                                                   "manifest.json (ltrtrace genotype on read sets writes it)")):
        source.run(tmp_path, tmp_path / "o", 1)
    run = _run_dir(files, tmp_path)
    (run / "reads.genotypes.tsv").unlink()
    with pytest.raises(UsageError, match="reads.genotypes.tsv"):
        source.run(run, tmp_path / "o", 1)
    run = _run_dir(files, tmp_path, name="two")
    (run / "manifest.json").write_text(json.dumps({"samples": [{"label": "a", "role": None}]}))
    with pytest.raises(UsageError, match=re.escape(f"{run / 'manifest.json'}: no sample of role reference with a "
                                                   "fasta and an annotation")):
        source.run(run, tmp_path / "o", 1)
    for text, message in (("{not json", "is unreadable"), ("[]", "is not an ltrtrace manifest")):
        (run / "manifest.json").write_text(text)
        with pytest.raises(UsageError, match=re.escape(f"{run / 'manifest.json'} {message}")):
            source.run(run, tmp_path / "o", 1)


def test_command_line_writes_the_run_folder(files, tmp_path):
    """`ltrtrace source RUN` writes into RUN, whose manifest.json (the read genotyping run's) stays as it was; or into
    -o, which gets its own manifest.json: the source run, on RUN's samples."""
    run = _run_dir(files, tmp_path)
    before = (run / "manifest.json").read_bytes()
    assert cli.main(["source", str(run), "-t", "2"]) == 0
    assert (run / "reads.sources.tsv").is_file() and (run / "reads.copies.tsv").is_file()
    assert (run / "manifest.json").read_bytes() == before
    other = tmp_path / "elsewhere"
    assert cli.main(["source", str(run), "-o", str(other), "-t", "2"]) == 0
    assert (other / "reads.locality.tsv").is_file()
    assert (other / "reads.sources.tsv").read_text() == (run / "reads.sources.tsv").read_text()
    doc = json.loads((other / "manifest.json").read_text())
    assert doc["params"]["command"] == "source" and doc["params"]["run"] == str(run) and doc["route"] is None
    assert [(s["label"], s["role"]) for s in doc["samples"]] == [("ref", "reference"), ("s1", None)]
    assert doc["samples"][0]["fasta"] == str(files["fa"]) and list(doc["code"]) == ["analyses.source"]
    assert (run / "manifest.json").read_bytes() == before


def test_source_keeps_the_read_genotyping_runs_record(files, tmp_path, monkeypatch, capsys):
    """`source` never removes RUN's manifest.json, the read genotyping run's record, even when it fails; into another
    folder it removes that folder's earlier manifest.json once its checks have passed, so a run that fails there leaves
    no record of another run."""
    run = _run_dir(files, tmp_path)
    before = (run / "manifest.json").read_bytes()
    other = tmp_path / "other"
    other.mkdir()
    (other / "manifest.json").write_text('{"params": {"command": "source"}}\n')  # an earlier source run's record

    def boom(*args, **kwargs):
        raise RuntimeError("a planted failure")

    monkeypatch.setattr(source, "load_targets", boom)  # after the checks: the reference's elements
    assert cli.main(["source", str(run), "-t", "1"]) == 1
    assert (run / "manifest.json").read_bytes() == before
    assert cli.main(["source", str(run), "-o", str(other), "-t", "1"]) == 1
    assert not (other / "manifest.json").exists() and (run / "manifest.json").read_bytes() == before


def test_table_cells_print_as_ltrkin_printed_them(genome, files, tmp_path):
    """Booleans as True/False, integer columns without a decimal point, a column holding a NaN as floats with the
    NaN empty."""
    out = tmp_path / "fmt"
    _source(files, out)
    text = (out / "reads.sources.tsv").read_text().splitlines()
    head = text[0].split("\t")
    row = dict(zip(head, text[1 + [r["new_id"] for r in genome].index("child")].split("\t"), strict=True))
    assert row["same_chrom"] == "True" and row["n_tied"] == "1" and row["aln_len"].isdigit()
    assert row["distance"] not in ("", "nan") and row["mosaic_breakpoint"] == ""
    assert row["shared_derived"].endswith(".0")
