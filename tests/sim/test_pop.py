"""`ltrtrace sim pop` (ltrtrace.sim.pop.pop and its subcommand) on a toy real-like input: a FASTA of three 300-kb sequences
with 40 planted LTR-RTs (one nested) in a GFF3. One module run (seven genomes, three chromosomes, a C>A pulse)
serves the output tests; the others run small regions or fail before simulating."""
from __future__ import annotations

import pytest

pytest.importorskip("msprime")

import csv
import dataclasses
import gzip
import hashlib
import json
import os
import random
import subprocess
import sys
from pathlib import Path

import numpy as np
from builders import write_fasta
from sim import poptoy
from sim.poptoy import sim_main

from ltrtrace.base.constants import WINDOW_PAD
from ltrtrace.io.annotation import load_annotation
from ltrtrace.io.annotation.gff3 import write_gff3
from ltrtrace.io.elements import Element, make_eid
from ltrtrace.io.samples import read_samples
from ltrtrace.io.seqio import Fasta
from ltrtrace.sim.pop.build import PopGenome, add_gaps, build_chrom
from ltrtrace.sim.pop.donors import backbone_elements, pick_donors
from ltrtrace.sim.pop.events import draw_events
from ltrtrace.sim.pop.model import MSPRIME_HELP, Demog, Rates, ancestry, tracts
from ltrtrace.sim.pop.pop import MAX_SEED, SIMPLIFICATIONS, PopParams, check_outdir, simulate_pop, write_pop
from ltrtrace.sim.pop.truth import INSERTION_COLUMNS, annotate, tract_bed
from ltrtrace.sim.score import TRUTH_COLUMNS

SEQ_LEN = 300_000
STARTS = [10_000 + 22_000 * k for k in range(13)]  # 13 elements per sequence; s1 also holds a nested one: 40
LABELS = ["A1", "A2", "B1", "B2", "C1", "C2", "O1"]
GROUPS = ["A", "A", "B", "B", "C", "C", "O"]
ANNOTATED = LABELS[:6]
SEED = 3
NAMED = ("per-chromosome rule", "exponential", "receive no insertions", "any-carrier thinning", "whole span",
         "take no substitutions", "one mutation per SNP site", "indel times", "K at the root", "never G",
         "does not ask whether X carries", "tract end")  # each in one params.json simplification


def _toy(root: Path) -> tuple[Path, Path, list[Element]]:
    """s1, s2, s3 (300 kb each) with LTR-RTs planted at STARTS (LTRs 200-600 bp, internal 2-6 kb, every third on
    the minus strand), and one more nested in s1's first; written as root/toy.fa and root/toy.gff3."""
    rng = random.Random(7)
    seqs, recs = {}, []
    for j, seqid in enumerate(("s1", "s2", "s3")):
        spans = [(s, rng.randint(200, 600), rng.randint(2_000, 6_000)) for s in STARTS]
        if seqid == "s1":
            spans.append((STARTS[0] + spans[0][1] + 500, 150, 800))  # inside the first element's internal region
        seq, els = poptoy.toy_backbone(SEQ_LEN, spans, seed=j)
        seqs[seqid] = seq.decode()
        for i, e in enumerate(els):
            s, t = e.start + 1, e.end
            recs.append(Element("ref", make_eid("ref", seqid, s, t), f"{seqid}e{i}", seqid, s, t,
                                "-" if i % 3 == 2 else "+", (e.ltr_left[0] + 1, e.ltr_left[1]),
                                (e.ltr_right[0] + 1, e.ltr_right[1]), e.tsd, 0.01, f"fam{i % 4}", "Gypsy", None,
                                "high"))
    root.mkdir(parents=True, exist_ok=True)
    fa, gff = root / "toy.fa", root / "toy.gff3"
    write_fasta(fa, seqs)
    write_gff3(recs, gff)
    return fa, gff, recs


def _crossing(recs: list[Element]) -> tuple[Element, int]:
    """s3's seventh element and a region end 1,000 bp into it."""
    e = [r for r in recs if r.seqid == "s3"][6]
    return e, e.start + 1_000


@pytest.fixture(scope="module")
def toy(tmp_path_factory):
    return _toy(tmp_path_factory.mktemp("toy"))


def _pop(fa: Path, gff: Path, out: Path, regions: list[str], *extra: str) -> int:
    argv = ["pop", "--fasta", str(fa), "--annotation", str(gff), "--outdir", str(out)]
    for r in regions:
        argv += ["--region", r]
    return sim_main([*argv, *extra])


def _main_regions(recs: list[Element]) -> list[str]:
    _, cut = _crossing(recs)
    return [f"s1:1-{SEQ_LEN}", f"s2:1-{SEQ_LEN}", f"s3:1-{cut}"]


MAIN_ARGS = ("--seed", str(SEED), "--f", "0.2", "--ta", "2000")


@pytest.fixture(scope="module")
def run(toy, tmp_path_factory):
    """The module run: seven genomes on s1, s2 and s3 up to a cut through s3's seventh element; threads 1."""
    fa, gff, recs = toy
    out = tmp_path_factory.mktemp("pop") / "out"
    assert _pop(fa, gff, out, _main_regions(recs), *MAIN_ARGS, "-t", "1") == 0
    return out


def _files(root: Path) -> list[str]:
    """The files under root, but FASTA indexes (reading a FASTA writes one)."""
    return sorted(str(p.relative_to(root)) for p in root.rglob("*") if p.is_file() and p.suffix != ".fai")


def _tsv(path: Path) -> list[dict[str, str]]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _headers(fa: Path) -> list[str]:
    return [ln[1:].strip() for ln in fa.read_text().splitlines() if ln.startswith(">")]


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _error_lines(err: str) -> list[str]:
    return [ln for ln in err.splitlines() if "error:" in ln]


# ---------------------------------------------------------------- outputs


def test_pop_writes_outputs(run, toy):
    fa, gff, recs = toy
    truth = [f"truth/{x}__{y}.tsv" for x in ANNOTATED for y in LABELS if y != x]
    assert len(truth) == 36
    want = sorted([f"{lab}.fa" for lab in LABELS] + [f"{lab}.gff3" for lab in ANNOTATED] + truth
                  + ["samples.tsv", "params.json", "truth/insertions.tsv", "truth/fractions.tsv",
                     "truth/tracts.A1.bed", "truth/tracts.A2.bed"])
    assert _files(run) == want
    # the samples sheet: read_samples loads it; paths relative to it; O1 is the unannotated outgroup
    assert (run / "samples.tsv").read_text().splitlines()[0] == "label\tfasta\tannotation\tgroup\trole"
    sheet = read_samples(run / "samples.tsv")
    assert [s.label for s in sheet] == LABELS and [s.group for s in sheet] == GROUPS
    assert [s.role for s in sheet] == [None] * 6 + ["outgroup"]
    assert [s.annotation for s in sheet] == [str(run / f"{lab}.gff3") for lab in ANNOTATED] + [None]
    assert [s.fasta for s in sheet] == [run / f"{lab}.fa" for lab in LABELS]
    assert all(_headers(run / f"{lab}.fa") == ["c01", "c02", "c03"] for lab in LABELS)
    # annotations load, names unique; every element of a truth file is X's
    for lab in ANNOTATED:
        with Fasta(run / f"{lab}.fa") as f:
            els = load_annotation(str(run / f"{lab}.gff3"), lab, f)
        assert els and len({e.name for e in els}) == len(els)
    for path in truth:
        rows = _tsv(run / path)
        x = Path(path).stem.split("__")[0]
        assert list(rows[0]) == list(TRUTH_COLUMNS) and {r["genome"] for r in rows} == {x}
        assert len({r["name"] for r in rows}) == len(rows)
    ins = _tsv(run / "truth/insertions.tsv")
    assert list(ins[0]) == [*INSERTION_COLUMNS, *(f"fate_{lab}" for lab in LABELS)]
    assert {r["chrom"] for r in ins} == {"c01", "c02", "c03"}
    # params.json: inputs with their digests, versions, trees, event counts, truth counts, simplifications
    p = json.loads((run / "params.json").read_text())
    assert p["seed"] == SEED
    assert p["regions"] == [["s1", 1, SEQ_LEN], ["s2", 1, SEQ_LEN], ["s3", 1, _crossing(recs)[1]]]
    assert p["fasta"] == {"path": str(fa.resolve()), "sha256": _sha(fa)}
    assert p["annotation"] == {"path": str(gff), "files": [{"path": str(gff.resolve()), "sha256": _sha(gff)}]}
    assert p["params"]["demog"]["f"] == 0.2 and p["params"]["demog"]["ta"] == 2000
    assert p["params"]["dropout"] == 0.1 and p["params"]["jitter"] == 0.05 and p["params"]["gap_per_mb"] == 2.5
    assert set(p["versions"]) == {"ltrtrace", "msprime", "tskit"}
    assert set(p["trees"]) == {"c01", "c02", "c03"} and all(n >= 1 for n in p["trees"].values())
    simple = p["simplifications"]  # every approximation the sim modules' docstrings name, one line each
    assert simple == list(SIMPLIFICATIONS) and len(simple) == len(NAMED)
    assert all(sum(key in s for s in simple) == 1 for key in NAMED), simple
    ev = p["events"]
    assert set(ev["dropped"]["total"]) == {"ins_n", "ins_end", "snp_n", "ins_after_loss", "nested_after_loss",
                                           "loss_after_loss", "indel_in_no_genome", "gap_off_chromosome"}
    for kind in ("requested", "realized", "dropped"):
        assert set(ev[kind]["per_chromosome"]) == {"c01", "c02", "c03"}
        total = ev[kind]["total"]
        for k in total:  # totals are the sums of the per-chromosome and per-genome details
            assert total[k] == sum(d.get(k, 0) for d in ev[kind]["per_chromosome"].values()) \
                + sum(d.get(k, 0) for d in ev[kind].get("per_genome", {}).values())
    t, r, d = ev["requested"]["total"], ev["realized"]["total"], ev["dropped"]["total"]
    assert t["insertions"] == r["insertions"] + d["ins_n"] + d["ins_end"] + d["ins_after_loss"] \
        + d["nested_after_loss"]
    assert t["losses"] == r["losses"] + d["loss_after_loss"] and t["snps"] == r["snps"] + d["snp_n"]
    assert t["gaps"] == r["gaps"] + d["gap_off_chromosome"]
    # realized indels are those some genome holds; the rest are dropped as in no genome
    assert t["indels"] == r["indels"] + d["indel_in_no_genome"] and r["indels"] > 0
    assert r["indels"] == r["indel_insertions"] + r["indel_deletions"]
    # the backbone, known by construction: 14 elements on s1 (one nested), 13 on s2, 6 before the cut on s3
    cut = _crossing(recs)[1]
    assert p["backbone"]["annotation_elements"] == 40
    assert p["backbone"]["chromosomes"] == {
        "c01": {"region": f"s1:1-{SEQ_LEN}", "length": SEQ_LEN, "elements": 14},
        "c02": {"region": f"s2:1-{SEQ_LEN}", "length": SEQ_LEN, "elements": 13},
        "c03": {"region": f"s3:1-{cut}", "length": cut, "elements": 6}}
    assert r["insertions"] == len(ins) > 0
    for path in truth:  # truth state counts equal the files'
        x, y = Path(path).stem.split("__")
        counts: dict[str, int] = {}
        for row in _tsv(run / path):
            counts[row["state"]] = counts.get(row["state"], 0) + 1
        assert p["truth"][x][y] == counts
    # tracts: BED records inside their chromosomes; fractions.tsv's lifted bp are the BED's
    fr = {row["label"]: row for row in _tsv(run / "truth/fractions.tsv")}
    assert list(fr) == LABELS
    for lab in ("A1", "A2"):
        bed = [ln.split("\t") for ln in (run / f"truth/tracts.{lab}.bed").read_text().splitlines()]
        assert int(fr[lab]["lifted_bp"]) == sum(int(b) - int(a) for _, a, b in bed)
    for lab in LABELS:  # fraction: of the backbone (the regions' bp); lifted_fraction: of the genome's own length
        row = fr[lab]
        assert int(row["backbone_bp"]) == 2 * SEQ_LEN + cut
        assert int(row["genome_bp"]) == sum(len(ln) for ln in (run / f"{lab}.fa").read_text().splitlines()
                                            if not ln.startswith(">"))
        assert row["fraction"] == f"{int(row['tract_bp']) / int(row['backbone_bp']):.6f}"
        assert row["lifted_fraction"] == f"{int(row['lifted_bp']) / int(row['genome_bp']):.6f}"
    assert all(fr[lab]["tract_bp"] == "0" for lab in LABELS[2:])
    assert sum(int(fr[lab]["tract_bp"]) for lab in ("A1", "A2")) > 0  # f 0.2 at 2,000 generations on 900 kb


def test_pop_help_epilog(capsys, monkeypatch):
    """`pop --help` ends with every Demog and Rates field (the --demography JSON's keys) at its default, where the
    defaults come from, what the lineage pulse codes need, the convention of insertions.tsv's pos, and where the
    simplifications are listed."""
    monkeypatch.setenv("COLUMNS", "80")
    assert sim_main(["pop", "--help"]) == 0
    out = capsys.readouterr().out
    epilog = " ".join(out[out.rindex("defaults:"):].split())  # after the last option, unwrapped
    for obj in (Demog(), Rates()):
        for f in dataclasses.fields(obj):
            v = getattr(obj, f.name)
            shown = ("none" if v is None else v if isinstance(v, str)
                     else f"[{', '.join(f'{x:g}' for x in v)}]" if isinstance(v, tuple) else f"{v:g}")
            assert f"{f.name} {shown}" in epilog, f.name
    assert "defaults: the demography and rates below" in epilog and "measured on real genomes" in epilog
    assert "C>A1 needs sub_a, C1>A1 and A1>C1 need sub_a and sub_c" in epilog
    assert "insertions.tsv: pos is 0-based" in epilog and "params.json lists the simplifications" in epilog
    assert max(len(ln) for ln in out.splitlines()) <= 80


def test_pop_deterministic(run, toy, tmp_path):
    fa, gff, recs = toy
    out = tmp_path / "out"
    assert _pop(fa, gff, out, _main_regions(recs), *MAIN_ARGS, "-t", "2") == 0
    assert _files(out) == _files(run)
    for f in _files(run):  # every file, FASTA, GFF3, truth and params.json included, byte for byte
        assert (out / f).read_bytes() == (run / f).read_bytes(), f


def test_pop_demography_json(toy, tmp_path):
    fa, gff, _ = toy
    dem = tmp_path / "demography.json"
    dem.write_text(json.dumps({"t_o": 300_000, "f": 0.05, "ta": 1_000, "per_taxon": 1, "mu": 2e-8}))
    out = tmp_path / "out"
    assert _pop(fa, gff, out, ["s1:1-30000"], "--demography", str(dem), "--f", "0.2", "--dropout-genome",
                "A1=0") == 0
    p = json.loads((out / "params.json").read_text())
    assert p["params"]["demog"]["t_o"] == 300_000 and p["params"]["demog"]["ta"] == 1_000
    assert p["params"]["demog"]["f"] == 0.2  # the command line overrides the JSON
    assert p["params"]["demog"]["per_taxon"] == 1 and p["params"]["rates"]["mu"] == 2e-8
    assert [s.label for s in read_samples(out / "samples.tsv")] == ["A1", "B1", "C1", "O1"]
    # --dropout-genome A1=0: A1's annotation holds every alive element; B1 and C1 keep --dropout
    assert p["params"]["dropout_genome"] == [["A1", 0.0]]
    a1 = p["genomes"]["A1"]
    assert a1["dropout"] == 0.0 and a1["annotated"] == a1["alive"] > 0
    assert p["genomes"]["B1"]["dropout"] == p["genomes"]["C1"]["dropout"] == 0.1


def test_pop_deterministic_across_processes(toy, tmp_path):
    """Two interpreters with different string-hash seeds (and thread counts)
    write the same bytes."""
    fa, gff, _ = toy
    src = Path(sys.modules["ltrtrace"].__file__).parents[1]
    code = "import sys\nfrom ltrtrace.cli import main\nsys.exit(main(sys.argv[1:]))"
    outs = []
    for hashseed, threads in (("0", "1"), ("12345", "2")):
        out = tmp_path / f"out{hashseed}"
        env = {**os.environ, "PYTHONHASHSEED": hashseed,
               "PYTHONPATH": os.pathsep.join(p for p in (str(src), os.environ.get("PYTHONPATH")) if p)}
        res = subprocess.run([sys.executable, "-c", code, "sim", "pop", "--fasta", str(fa), "--annotation", str(gff),
                              "--region", "s1:1-60000", "--region", "s2:1-60000", "--per-taxon", "1", "--seed", "4",
                              "--f", "0.3", "--ta", "2000", "-t", threads, "--outdir", str(out)],
                             env=env, capture_output=True, text=True, check=False)
        assert res.returncode == 0, res.stderr
        outs.append(out)
    a, b = outs
    assert _files(a) == _files(b) and "truth/tracts.A1.bed" in _files(a)
    for f in _files(a):
        assert (a / f).read_bytes() == (b / f).read_bytes(), f


def test_pop_streams(toy):
    """The streams of pop.py's docstring, rebuilt step by step: chromosome c from [seed, c] (ancestry's seed its
    first draw, then draw_events and build_chrom), gaps of genome g from [seed, 100, g], its annotation from
    [seed, 200, g], donors from [seed, 0]. Dense gaps and regions that end and start exactly at an element (a copy
    end at a chromosome end) make near gaps fall off: the gap counts of params.json."""
    fa, gff, recs = toy
    first, last = [r for r in recs if r.seqid == "s2"][2], [r for r in recs if r.seqid == "s1"][10]
    seed, regions = 5, [("s2", 1, first.end), ("s1", last.start, last.start + 59_999)]
    params = PopParams(Demog(per_taxon=1, f=1.0, ta=3_000), Rates(), gap_per_mb=1_000)  # f = 1: A1 all tract
    res = simulate_pop(regions, fa, str(gff), params, seed, 1)
    built, chrom_tracts = [], []
    with Fasta(fa) as f:
        els = load_annotation(str(gff), "ref", f)
        donors = pick_donors(els, f, seed)
        for c, (seqid, a, b) in enumerate(regions, 1):
            bb, be = f.fetch(seqid, a - 1, b).encode(), backbone_elements(els, seqid, a, b)
            rng = np.random.default_rng([seed, c])
            ts = ancestry(params.demog, params.rates, len(bb), int(rng.integers(1, 2**32 - 1)))
            ev = draw_events(ts, bb, be, donors, params.rates, rng, chrom=c)
            built.append(build_chrom(c - 1, bb, be, ts, ev, donors, params.rates, rng))
            chrom_tracts.append(tracts(ts, params.demog))
    assert [g.label for g in res.genomes] == ["A1", "B1", "C1", "O1"]
    gaps = res.meta["events"]
    off_total = 0
    for i, g in enumerate(res.genomes):
        mine = PopGenome(g.label, g.group, ["c01", "c02"], [x[i][0] for x in built],
                         {e: cp for x in built for e, cp in x[i][1].items()}, [], [x[i][2] for x in built])
        n = int(1_000 * sum(len(x) for x in mine.seqs) / 1e6 + 0.5)
        (off,) = add_gaps([mine], params.gap_per_mb, [np.random.default_rng([seed, 100, i])])
        assert mine.seqs == g.seqs and mine.gaps == g.gaps and mine.copies == g.copies
        assert gaps["requested"]["per_genome"][g.label] == {"gaps": n}
        assert gaps["realized"]["per_genome"][g.label] == {"gaps": n - off} and len(g.gaps) == n - off
        assert gaps["dropped"]["per_genome"][g.label] == {"gap_off_chromosome": off}
        off_total += off
        if g.group != "O":
            assert annotate(mine, 0.1, 0.05, np.random.default_rng([seed, 200, i])) == res.annotations[g.label]
    assert off_total > 0 and gaps["dropped"]["total"]["gap_off_chromosome"] == off_total
    # the C>A pulse's tracts of A1 (sample 0): whole chromosomes at f = 1, lifted; no other sample has any
    a1 = res.genomes[0]
    assert list(res.tracts) == ["A1"]
    assert res.tracts["A1"] == [("c01", 0, len(a1.seqs[0])), ("c02", 0, len(a1.seqs[1]))]
    assert res.tracts["A1"] == tract_bed(a1, {k: t.get(0, []) for k, t in enumerate(chrom_tracts)})
    assert [set(t) for t in chrom_tracts] == [{0}, {0}]


def test_pop_flags_override_json(toy, tmp_path):
    """--pulse, --bottleneck and --per-taxon override the JSON's values, as --f does."""
    fa, gff, _ = toy
    dem = tmp_path / "demography.json"
    dem.write_text(json.dumps({"pulse": "C>A", "f": 0.1, "ta": 3_000, "per_taxon": 2, "bottleneck": None}))
    out = tmp_path / "out"
    assert _pop(fa, gff, out, ["s1:1-30000"], "--demography", str(dem), "--pulse", "A>C", "--bottleneck",
                "2000,500", "--per-taxon", "1", "--f", "0.5", "--ta", "1000") == 0
    d = json.loads((out / "params.json").read_text())["params"]["demog"]
    assert (d["pulse"], d["bottleneck"], d["per_taxon"], d["f"], d["ta"]) == ("A>C", [2000.0, 500.0], 1, 0.5, 1000)
    assert [s.label for s in read_samples(out / "samples.tsv")] == ["A1", "B1", "C1", "O1"]
    assert sorted(x.name for x in (out / "truth").glob("tracts.*.bed")) == ["tracts.C1.bed"]  # A>C: C receives


def test_pop_lineage_pulse(toy, tmp_path):
    """A pulse into A1's own lineage (C>A1, sub_a) at f = 1 makes all of A1 a tract and none of
    A2, whose tract file is written empty."""
    fa, gff, _ = toy
    dem = tmp_path / "demography.json"
    dem.write_text(json.dumps({"sub_a": [4_000, 10_000]}))
    out = tmp_path / "out"
    assert _pop(fa, gff, out, ["s1:1-30000"], "--demography", str(dem), "--pulse", "C>A1", "--f", "1",
                "--ta", "1000") == 0
    assert (out / "truth/tracts.A2.bed").read_text() == ""
    a1 = (out / "truth/tracts.A1.bed").read_text().splitlines()
    fr = {row["label"]: row for row in _tsv(out / "truth/fractions.tsv")}
    assert a1 and fr["A1"]["fraction"] == "1.000000" and fr["A2"]["tract_bp"] == "0"
    assert sum(int(b) - int(a) for _, a, b in (ln.split("\t") for ln in a1)) == int(fr["A1"]["genome_bp"])


def test_pop_two_regions_one_sequence(toy, tmp_path):
    fa, gff, recs = toy
    out = tmp_path / "out"
    regions = ["s1:1-150000", "s1:150001-300000"]
    assert _pop(fa, gff, out, regions, "--per-taxon", "1") == 0
    assert _headers(out / "A1.fa") == ["c01", "c02"]
    regions_json = json.loads((out / "params.json").read_text())["regions"]
    assert regions_json == [["s1", 1, 150_000], ["s1", 150_001, 300_000]]
    # backbone rows: each element named by its source span sits on its own window's chromosome
    inside = {f"s1:{r.start}-{r.end}": "c01" if r.end <= 150_000 else "c02" for r in recs if r.seqid == "s1"}
    rows = [r for r in _tsv(out / "truth/A1__B1.tsv") if r["origin"] == "backbone"]
    assert rows and all(inside[r["name"]] == r["seqid"] for r in rows)


def test_pop_edge_element_excluded(toy, run):
    fa, gff, recs = toy
    edge, cut = _crossing(recs)
    with Fasta(fa) as f:
        els = load_annotation(str(gff), "ref", f)
    got = backbone_elements(els, "s3", 1, cut)
    want = [r for r in recs if r.seqid == "s3" and r.end <= cut]
    assert len(want) == 6 and edge.start <= cut < edge.end
    assert [(e.eid, e.start, e.end, e.ltr_left, e.ltr_right, e.strand, e.origin, e.parent) for e in got] == [
        (f"s3:{r.start}-{r.end}", r.start - 1, r.end, (r.ltr_left[0] - 1, r.ltr_left[1]),
         (r.ltr_right[0] - 1, r.ltr_right[1]), r.strand, "backbone", None) for r in want]
    # shifted: a window starting at base 20,001 moves every coordinate 20,000 down; the nested one keeps its host
    s1 = backbone_elements(els, "s1", 1, 40_000)
    assert [(e.start, e.parent) for e in s1][:2] == [(STARTS[0], None), (s1[1].start, s1[0].eid)]
    assert [e.start for e in backbone_elements(els, "s1", 20_001, 40_000)] == [STARTS[1] - 20_000]
    # and in the run: the crossing element is nobody's backbone element
    name = f"s3:{edge.start}-{edge.end}"
    names = {r["name"] for f in (run / "truth").glob("*__*.tsv") for r in _tsv(f)}
    assert name not in names and any(n.startswith("s3:") for n in names)


# ---------------------------------------------------------------- failing fast


def test_pop_regions_overlap_fails(toy, tmp_path, capsys):
    fa, gff, _ = toy
    out = tmp_path / "out"
    assert _pop(fa, gff, out, ["s1:1-100000", "s2:1-50000", "s1:100000-200000"]) == 2
    err = capsys.readouterr().err.splitlines()
    assert len(err) == 1 and err[0].startswith("error: --region ") and "overlap" in err[0]
    assert not out.exists()


@pytest.mark.parametrize("region, names", [
    ("s9:1-1000", "sequence 's9' is not in {fa}"),
    ("s1:1000", "s1:1000: expected SEQ:START-END"),
    (f"s1:1-{SEQ_LEN + 1}", f"is outside s1 (1-{SEQ_LEN})"),
], ids=["unknown-sequence", "malformed", "off-the-sequence"])
def test_pop_bad_region_is_a_usage_error(toy, tmp_path, capsys, region, names):
    """A --region that names no sequence of the FASTA, is malformed or runs off its sequence is a usage error: exit
    2 and one `error:` line naming --region (and the FASTA), before anything is simulated or written."""
    fa, gff, _ = toy
    out = tmp_path / "out"
    assert _pop(fa, gff, out, [region], "--per-taxon", "1") == 2
    err = capsys.readouterr().err.splitlines()
    assert len(err) == 1 and err[0].startswith("error: --region ") and names.format(fa=fa) in err[0], err
    assert not (out / "samples.tsv").exists()


def test_pop_bad_demography_is_a_usage_error(toy, tmp_path, capsys):
    """A --demography file that is not valid JSON is a usage error naming --demography and the file."""
    fa, gff, _ = toy
    (tmp_path / "d.json").write_text('{"f": 0.2,')
    assert _pop(fa, gff, tmp_path / "out", ["s1:1-30000"], "--demography", str(tmp_path / "d.json")) == 2
    err = capsys.readouterr().err.splitlines()
    assert len(err) == 1 and err[0].startswith(f"error: --demography {tmp_path / 'd.json'}: not valid JSON"), err


def test_pop_dropout_guards(toy, tmp_path, capsys):
    fa, gff, _ = toy
    out = tmp_path / "out"
    for spec in ("A1=1", "A1=1.5", "A1=x"):  # a value of 1 or more, or not a number: argparse's usage error
        code = _pop(fa, gff, out, [f"s1:1-{SEQ_LEN}"], "--dropout-genome", "A2=0.3", spec)
        lines = _error_lines(capsys.readouterr().err)
        assert code == 2 and len(lines) == 1 and "--dropout-genome" in lines[0]
    for spec in ("O1=0.2", "Z9=0.2"):  # the outgroup or an unknown label: checked after parsing, before simulating
        assert _pop(fa, gff, out, [f"s1:1-{SEQ_LEN}"], "--dropout-genome", "A2=0.3", spec) == 2
        err = capsys.readouterr().err.splitlines()
        assert len(err) == 1 and err[0].startswith("error: --dropout-genome") and spec.split("=")[0] in err[0], err
    assert not out.exists()


@pytest.mark.parametrize("extra, named, json_text", [
    (["--f", "0.2", "--ta", "20000"], ["--f", "--ta"], None),
    (["--f", "0.2", "--pulse", "C>A1"], ["--f", "--pulse"], None),
    (["--bottleneck", "1000,20000"], ["--bottleneck"], None),
    (["--ta", "20000"], ["--ta", "--demography"], '{"f": 0.2}'),
], ids=["ta-above-the-split", "lineage-pulse-without-its-lineage", "bottleneck-longer-than-the-split",
        "option-against-the-file"])
def test_pop_refused_demography_is_a_usage_error(toy, tmp_path, capsys, monkeypatch, extra, named, json_text):
    """A --f, --ta, --pulse or --bottleneck value the demography refuses (alone, or with a --demography file) is a
    usage error found before anything is simulated: exit 2 and one `error:` line naming the options."""
    fa, gff, _ = toy
    out = tmp_path / "out"
    if json_text is not None:
        (tmp_path / "d.json").write_text(json_text)
        extra = [*extra, "--demography", str(tmp_path / "d.json")]

    def never(*a, **k):
        raise AssertionError("simulated despite a refused demography")
    monkeypatch.setattr("ltrtrace.sim.pop.pop.simulate_pop", never)
    assert _pop(fa, gff, out, [f"s1:1-{SEQ_LEN}"], "--per-taxon", "1", *extra) == 2
    err = capsys.readouterr().err.splitlines()
    assert len(err) == 1 and err[0].startswith("error: --") and all(o in err[0] for o in named), err
    assert not out.exists()


def test_pop_region_too_short(toy, tmp_path, capsys):
    fa, gff, _ = toy
    with Fasta(fa) as f:
        longest = max(len(d.seq) for d in pick_donors(load_annotation(str(gff), "ref", f), f, 1))
    need = 2 * WINDOW_PAD + longest
    out = tmp_path / "out"
    assert _pop(fa, gff, out, [f"s2:1-{SEQ_LEN}", f"s1:1-{need - 1}"], "--per-taxon", "1") == 2
    err = capsys.readouterr().err.splitlines()
    assert len(err) == 1 and err[0].startswith(f"error: --region s1:1-{need - 1} ") and str(need) in err[0]
    assert not out.exists()
    assert _pop(fa, gff, out, [f"s1:1-{need}"], "--per-taxon", "1") == 0  # exactly long enough


def test_pop_seed_range(toy, tmp_path, capsys):
    """A seed is one 32-bit word, so streams of different seeds never alias: outside [0, 2^32 - 1]
    --seed is a usage error (exit 2) and simulate_pop a ValueError; 2^32 - 1 itself runs."""
    fa, gff, _ = toy
    out = tmp_path / "out"
    for bad in (str(MAX_SEED + 1), "-1"):
        argv = ["pop", "--fasta", str(fa), "--annotation", str(gff), "--outdir", str(out), "--seed", bad,
                "--region", "s1:1-30000"]
        code = sim_main(argv)
        lines = _error_lines(capsys.readouterr().err)
        assert code == 2 and len(lines) == 1 and "--seed" in lines[0] and str(MAX_SEED) in lines[0]
    with pytest.raises(ValueError, match="seed"):
        simulate_pop([("s1", 1, 30_000)], fa, str(gff), PopParams(Demog(per_taxon=1), Rates()), MAX_SEED + 1, 1)
    assert not out.exists()
    assert MAX_SEED == 2**32 - 1
    assert _pop(fa, gff, out, ["s1:1-30000"], "--per-taxon", "1", "--seed", str(MAX_SEED)) == 0
    assert json.loads((out / "params.json").read_text())["seed"] == MAX_SEED


def test_pop_refuses_stale_outdir(toy, tmp_path, capsys, monkeypatch):
    """An outdir holding anything but its own FASTA index cache (fai/) is a usage error (exit 2, one error line
    naming --outdir), found before anything is simulated, and write_pop refuses it too; fai/ alone, or no outdir, is
    fine."""
    fa, gff, _ = toy
    out = tmp_path / "out"
    (out / "truth").mkdir(parents=True)
    stale = out / "truth" / "tracts.A2.bed"
    stale.write_text("c01\t0\t10\n")

    def never(*a, **k):
        raise AssertionError("simulated despite a stale outdir")
    with monkeypatch.context() as m:
        m.setattr("ltrtrace.sim.pop.pop.simulate_pop", never)
        assert _pop(fa, gff, out, ["s1:1-30000"], "--per-taxon", "1") == 2
    err = capsys.readouterr().err.splitlines()
    assert len(err) == 1 and err[0].startswith(f"error: --outdir {out}") and "not empty" in err[0], err
    assert sorted(x.name for x in out.rglob("*")) == ["tracts.A2.bed", "truth"]  # untouched
    result = simulate_pop([("s1", 1, 30_000)], fa, str(gff), PopParams(Demog(per_taxon=1), Rates()), 1, 1)
    with pytest.raises(ValueError, match="not empty"):
        write_pop(out, result)
    assert stale.read_text() == "c01\t0\t10\n"
    cached = tmp_path / "cached"
    (cached / "fai").mkdir(parents=True)
    (cached / "fai" / "x.toy.fa.fai").write_text("")
    assert _pop(fa, gff, cached, ["s1:1-30000"], "--per-taxon", "1") == 0
    assert (cached / "samples.tsv").is_file() and (cached / "fai" / "x.toy.fa.fai").is_file()
    check_outdir(tmp_path / "new")  # missing: fine
    odd = tmp_path / "odd"
    odd.mkdir()
    (odd / "fai").write_text("")  # a file named fai is not the cache
    with pytest.raises(ValueError, match="not empty"):
        check_outdir(odd)
    with pytest.raises(ValueError, match="not a directory"):
        check_outdir(odd / "fai")


def test_pop_builds_each_annotation_row_once(toy, monkeypatch):
    """The donor and backbone scans pass over the annotation's elements several times (backbone: once per chromosome):
    the loader's rows are built into elements once, not once per pass over its columnar table."""
    from ltrtrace.sim.pop import pop as pop_module

    fa, gff, _ = toy
    passes, real = [], pop_module.load_annotation

    class Counted:  # the loader's table, counting the passes over its rows
        def __init__(self, table):
            self.table = table

        def __len__(self):
            return len(self.table)

        def __iter__(self):
            passes.append(1)
            return iter(self.table)

    monkeypatch.setattr(pop_module, "load_annotation", lambda *args, **kwargs: Counted(real(*args, **kwargs)))
    simulate_pop([("s1", 1, 30_000), ("s2", 1, 30_000)], fa, str(gff), PopParams(Demog(per_taxon=1), Rates()), 1, 1)
    assert len(passes) == 1


def test_pop_without_msprime(toy, tmp_path):
    # a fresh interpreter that cannot import msprime: one error line naming the extra, nothing written
    fa, gff, _ = toy
    src = Path(sys.modules["ltrtrace"].__file__).parents[1]
    code = ("import sys; sys.modules['msprime'] = None; sys.modules['tskit'] = None\n"
            "from ltrtrace.cli import main\n"
            "sys.exit(main(sys.argv[1:]))")
    out = tmp_path / "out"
    env = {**os.environ, "PYTHONPATH": os.pathsep.join(p for p in (str(src), os.environ.get("PYTHONPATH")) if p)}
    res = subprocess.run([sys.executable, "-c", code, "sim", "pop", "--fasta", str(fa), "--annotation", str(gff),
                          "--region", "s1:1-30000", "--outdir", str(out)], env=env, capture_output=True, text=True,
                         check=False)
    assert res.returncode == 1 and res.stderr.splitlines() == [f"error: {MSPRIME_HELP}"]
    assert not out.exists()


# ---------------------------------------------------------------- end to end


@pytest.mark.slow
def test_pop_end_to_end(toy, tmp_path):
    from ltrtrace.cli import main as cli_main

    fa, gff, _ = toy
    sim = tmp_path / "sim"
    assert _pop(fa, gff, sim, [f"s1:1-{SEQ_LEN}", f"s2:1-{SEQ_LEN}", f"s3:1-{SEQ_LEN}"], "--per-taxon", "1",
                "--seed", "2", "--f", "0.2", "--ta", "2000") == 0
    run = tmp_path / "run"
    assert cli_main(["geneflow", str(sim / "samples.tsv"), "-o", str(run), "-t", "8", "--trio", "A1,B1,C1"]) == 0
    stats = _tsv(run / "stats.tsv")
    assert [r for r in stats if r["test"] == "D_pav" and r["statistic"] == "D"]
    pair = run / "pairs" / "C1__A1.calls.tsv.gz"
    with gzip.open(pair, "rt") as fh:
        assert fh.readline().startswith("eid")
    prefix = tmp_path / "score" / "C1__A1"
    assert sim_main(["score", "--calls", str(pair), "--elements", str(run / "work" / "elements" / "C1.tsv"),
                     "--truth", str(sim / "truth" / "C1__A1.tsv"), "--out", str(prefix)]) == 0
    rows = {(r["state"], r["tier_set"]): r for r in _tsv(Path(f"{prefix}.score.tsv"))}
    for state in ("P", "E"):
        r = rows[(state, "all")]
        assert int(r["tp"]) + int(r["fn"]) > 0 and int(r["tp"]) > 0, r


def test_sim_cli_has_only_pop_and_score(capsys):
    """`ltrtrace sim` offers exactly the subcommands pop and score; any other (pair) is a usage error."""
    assert sim_main(["--help"]) == 0
    out = capsys.readouterr().out
    assert all(f"\n    {c} " in out for c in ("pop", "score"))
    assert "\n    pair " not in out
    assert sim_main(["pair"]) == 2
