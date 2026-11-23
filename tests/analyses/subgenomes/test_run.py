"""The subgenomes command end to end on small simulated polyploids, and its command line."""
import csv
import json
import os
import re
import shutil
import subprocess
import sys

import builders
import numpy as np
import pytest
from analyses.subgenomes import sgtools as sg

from synltr import __version__, cli
from synltr.analyses.subgenomes import biology, dating, inputs, lineages, pipeline, run
from synltr.base import constants as cst
from synltr.sim.nucleus.model import merge_defaults

ENV = {**os.environ, "NUMPY_MADVISE_HUGEPAGE": "0", "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"}


def test_help_is_short_and_hides_dev_flags(capsys):
    assert cli.main(["subgenomes", "--help"]) == 0
    out = capsys.readouterr().out
    assert "SHEET" in out and "--assign" in out and "--mu" in out and "--seed" in out
    assert "--dev_" not in out and "--label" not in out and len(out.splitlines()) < 45


def test_missing_input_fails_fast(tmp_path):
    """A library or annotation the sheet names that does not exist is a usage error: exit 2, one `error:` line."""
    sh = sg.sheet(tmp_path / "s.tsv", [tmp_path / "no.fa"], [tmp_path / "no.gff3"])
    r = subprocess.run([sys.executable, "-m", "synltr", "subgenomes", str(sh), "-o", str(tmp_path / "o")],
                       capture_output=True, text=True, env=ENV)
    assert r.returncode == 2 and r.stderr.strip().splitlines() == [f"error: input not found: {tmp_path / 'no.fa'}"]
    assert not (tmp_path / "o").exists()


@pytest.mark.parametrize("extra, message", [
    (["-n", "1"], f"-n must be between 2 and {cst.MAX_N} (omit it to infer the structure)"),
    (["--mu", "-1"], "--mu must be > 0 (substitutions per site per year)"),
    (["--assign", "missing.tsv"], "input not found: {tmp}/missing.tsv"),
    (["--assign", "part.tsv", "-n", "2"], "--assign fixes the partition: -n does not apply"),
], ids=["n-1", "mu-negative", "assign-missing", "assign-with-n"])
def test_argument_checks_are_usage_errors(tmp_path, capsys, extra, message):
    """subgenomes' own argument checks stop it before anything is written: exit 2 and one `error:` line."""
    (tmp_path / "x.fa").write_text(">c1:1-10\nACGTACGTAC\n")
    (tmp_path / "part.tsv").write_text("c1\tA\n")
    sh = sg.sheet(tmp_path / "s.tsv", [tmp_path / "x.fa"], [tmp_path / "x.fa"])
    extra = [str(tmp_path / a) if a.endswith(".tsv") else a for a in extra]
    assert cli.main(["subgenomes", str(sh), "-o", str(tmp_path / "o"), *extra]) == 2
    assert capsys.readouterr().err.strip().splitlines() == [f"error: {message.format(tmp=tmp_path)}"]
    assert not (tmp_path / "o").exists()


def _truncated_library(d):
    """The small library gzipped and cut short (no end-of-stream marker)."""
    import gzip

    fa, ann, fai = sg.small_library(d)
    data = gzip.compress(fa.read_bytes())
    (d / "ltr.fa.gz").write_bytes(data[:len(data) // 2])
    return d / "ltr.fa.gz", ann, fai


def _empty_annotation(d):
    fa, _, fai = sg.small_library(d)
    (d / "empty.gff3").write_text("##gff-version 3\n")
    return fa, d / "empty.gff3", fai


# input fault: (inputs from a folder, extra options, the file or option its error line names: "lib", "ann", "fai",
# "assign" or an option)
INPUT_FAULTS = {
    "library-truncated": (_truncated_library, [], "lib"),
    "annotation-empty": (_empty_annotation, [], "ann"),
    "k2p-like-identity": (lambda d: sg.small_library(d, k2p=lambda s, j: 0.98 - 0.001 * j), [], "ann"),
    "k2p-missing": (lambda d: sg.small_library(d, k2p=lambda s, j: "NA" if j % 2 else 0.01), [], "ann"),
    "fai-lacks-a-sequence": (lambda d: sg.small_library(d, fai_seqs=["c1"]), [], "fai"),
    "one-sequence-phasable": (lambda d: sg.small_library(d, counts=(12, 5)), [], "lib"),
    "n-above-the-sequences": (sg.small_library, ["-n", "3"], "-n 3"),
    "assign-two-labels": (sg.small_library, ["--assign", "c1 A\nc2 B\nc1 B\n"], "assign"),
    "assign-no-lines": (sg.small_library, ["--assign", "# nothing\n"], "assign"),
    "assign-no-library-sequence": (sg.small_library, ["--assign", "x1 A\nx2 B\n"], "assign"),
    "assign-too-many-groups": (lambda d: sg.small_library(d, fai_seqs=[f"c{i}" for i in range(1, 15)]),
                               ["--assign", "".join(f"c{i} G{i}\n" for i in range(1, 15))], "assign"),
    "assign-one-sequence-phasable": (sg.small_library, ["--assign", "c1 A\nnosuch B\n"], "assign"),
    "no-shared-kmers": (lambda d: sg.small_library(d, shared=False), [], "lib"),
}


@pytest.mark.parametrize("fault", list(INPUT_FAULTS))
def test_input_faults_are_one_usage_error(tmp_path, capsys, fault):
    """A fault in subgenomes' inputs or options found while the run reads them (the library and its annotation, the
    .fai, --assign, -n against the library, k-mers that link no two sequences) is a usage error, whichever step finds
    it: exit 2 and exactly one `error:` line, naming the file or option, no traceback, no run record."""
    make, extra, names = INPUT_FAULTS[fault]
    fa, ann, fai = make(tmp_path / "in")
    if extra[:1] == ["--assign"]:
        (tmp_path / "assign.txt").write_text(extra[1])
        extra = ["--assign", str(tmp_path / "assign.txt")]
    sh = sg.sheet(tmp_path / "s.tsv", [fa], [ann], [fai], label="toy")
    assert cli.main(["subgenomes", str(sh), "-o", str(tmp_path / "o"), "-t", "2", *extra]) == 2
    err = capsys.readouterr().err
    errors = [ln for ln in err.splitlines() if "error:" in ln.lower()]
    named = {"lib": str(fa), "ann": str(ann), "fai": str(fai), "assign": str(tmp_path / "assign.txt")}.get(names, names)
    assert len(errors) == 1 and errors[0].startswith("error: ") and named in errors[0], err
    assert "Traceback" not in err and not (tmp_path / "o" / "manifest.json").exists()


@pytest.mark.parametrize("fault", [f for f in INPUT_FAULTS if f.startswith("assign-")])
def test_assign_faults_fail_before_the_kmer_stage(tmp_path, capsys, monkeypatch, fault):
    """--assign is read and checked against the library right after the library is read, before the k-mer stage: a
    fault in it fails in seconds, as the same one usage error."""
    make, extra, _ = INPUT_FAULTS[fault]
    fa, ann, fai = make(tmp_path / "in")
    (tmp_path / "assign.txt").write_text(extra[1])

    def kmers(*a, **k):
        raise AssertionError("the k-mer stage ran before --assign was checked")
    monkeypatch.setattr(run, "load_or_build_clades", kmers)
    sh = sg.sheet(tmp_path / "s.tsv", [fa], [ann], [fai], label="toy")
    assert cli.main(["subgenomes", str(sh), "-o", str(tmp_path / "o"), "-t", "2", "--assign",
                     str(tmp_path / "assign.txt")]) == 2
    errors = [ln for ln in capsys.readouterr().err.splitlines() if "error:" in ln.lower()]
    assert len(errors) == 1 and errors[0].startswith(f"error: --assign {tmp_path / 'assign.txt'}"), errors


def test_negative_seed_is_a_usage_error(tmp_path, capsys):
    """--seed below 0 is refused with the other option checks, before any input is read: exit 2, one error line
    naming --seed, nothing written."""
    fa, ann, fai = sg.small_library(tmp_path / "in")
    sh = sg.sheet(tmp_path / "s.tsv", [fa], [ann], [fai], label="toy")
    assert cli.main(["subgenomes", str(sh), "-o", str(tmp_path / "o"), "--seed", "-5"]) == 2
    errors = [ln for ln in capsys.readouterr().err.splitlines() if "error:" in ln]
    assert len(errors) == 1 and "--seed" in errors[0], errors
    assert not (tmp_path / "o").exists()


@pytest.mark.parametrize("argv,want", [(["-t", "8"], 8), (["--threads=3"], 3), (["--x", "-t"], 4), ([], 4)])
def test_threads_read_before_numba(argv, want):
    assert cli.thread_env(argv)["NUMBA_NUM_THREADS"] == str(want)


def test_threads_not_capped_by_environment(tmp_path):
    """NUMBA_NUM_THREADS and OMP_NUM_THREADS in the environment cap nothing: numba's pool is sized before numba loads,
    from -t in argv (default 4), BLAS runs single-threaded; setup() caps a -t above the pool, with a warning."""
    def run_(argv):
        code = ("import os, sys; os.environ['NUMBA_NUM_THREADS'] = '2'; os.environ['OMP_NUM_THREADS'] = '8'; "
                f"from synltr import cli; os.environ.update(cli.thread_env({argv!r})); import numba; "
                "from synltr.analyses.subgenomes import run; "
                f"a = run.options([{str(tmp_path / 'x.fa')!r}], [{str(tmp_path / 'x.fa')!r}], {str(tmp_path / 'o')!r}, "
                "threads=3); run.setup(a); print(numba.get_num_threads(), os.environ['OMP_NUM_THREADS'])")
        r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
        assert r.returncode == 0, r.stderr[-500:]
        return r
    (tmp_path / "x.fa").write_text(">c1:1-10\nACGTACGTAC\n")
    assert run_(["subgenomes"]).stdout.split() == ["3", "1"]
    r = run_(["subgenomes", "-t", "2"])
    assert r.stdout.split() == ["2", "1"] and "exceeds numba's pool of 2 threads" in r.stderr


def test_n_validation_and_removed_options(tmp_path, capsys):
    """-n outside 2..MAX_N stops with a clear message; options of earlier versions (--config, --features, --no_plots,
    --dev_*), the input flags and --label (taken from the sheet) are rejected."""
    (tmp_path / "x.fa").write_text(">c1:1-10\nACGTACGTAC\n")
    x = str(tmp_path / "x.fa")
    for n in ("1", "13"):
        with pytest.raises(SystemExit):
            sg.main(["--ltr_fasta", x, "--annotation", x, "--outdir", str(tmp_path / "o"), "-n", n])
        assert f"-n must be between 2 and {cst.MAX_N}" in capsys.readouterr().err
    sh = str(sg.sheet(tmp_path / "s.tsv", [x], [x]))
    for opt in (["--config", "c.txt"], ["--features", "f.bed"], ["--no_plots"], ["--dev_rule", "R0"],
                ["--dev_split_rule", "orthologs"], ["--dev_skip_dating"], ["--dev_scale", "2"], ["--label", "L"],
                ["--ltr_fasta", x]):
        assert cli.main(["subgenomes", sh, "-o", str(tmp_path / "o")] + opt) == 2


@pytest.mark.slow
def _sample(tmp_path, **kw):
    from synltr.io.samples import Sample
    return Sample(**dict(dict(label="g", fasta=None, annotation=str(tmp_path / "a.gff3"), group=None, role=None,
                              library=str(tmp_path / "l.fa")), **kw))


def test_sheet_row_gives_library_annotation_and_fai(tmp_path):
    """The sheet's library is the LTR-RT FASTA (several, comma-separated), its annotation the annotation; the fasta
    column gives the .fai: the genome FASTA's index beside it, or the .fai itself."""
    (tmp_path / "g.fa").write_text(">c1\nACGT\n")
    (tmp_path / "g.fa.fai").write_text("c1\t4\t4\t4\t5\n")
    lib, ann, fai = run.sample_inputs(_sample(tmp_path, fasta=tmp_path / "g.fa"))
    assert lib == [str(tmp_path / "l.fa")] and ann == [str(tmp_path / "a.gff3")] and fai == [str(tmp_path / "g.fa.fai")]
    assert run.sample_inputs(_sample(tmp_path, fasta=tmp_path / "g.fa.fai"))[2] == [str(tmp_path / "g.fa.fai")]
    assert run.sample_inputs(_sample(tmp_path))[2] is None
    two = run.sample_inputs(_sample(tmp_path, library=f"{tmp_path / 'l.fa'},{tmp_path / 'm.fa'}"))[0]
    assert two == [str(tmp_path / "l.fa"), str(tmp_path / "m.fa")]


def test_sheet_problems_are_usage_errors(tmp_path):
    from synltr.base.log import UsageError
    (tmp_path / "g.fa").write_text(">c1\nACGT\n")
    with pytest.raises(UsageError, match="fai"):
        run.sample_inputs(_sample(tmp_path, fasta=tmp_path / "g.fa"))                # no index beside the FASTA
    with pytest.raises(UsageError, match="library"):
        run.sample_inputs(_sample(tmp_path, library=None))
    with pytest.raises(UsageError, match="annotation"):
        run.sample_inputs(_sample(tmp_path, annotation=None))
    (tmp_path / "two.tsv").write_text("label\tlibrary\tannotation\na\tl.fa\ta.gff3\nb\tl.fa\ta.gff3\n")
    assert cli.main(["subgenomes", str(tmp_path / "two.tsv"), "-o", str(tmp_path / "o")]) == 2


def test_mu_must_be_positive(tmp_path):
    d, truth = sg.simulate(tmp_path, tree={"split": 2e6, "children": ["A", "B"]},
                         mergers=[{"time": 1e5, "lineages": ["A", "B"]}])
    for mu in ("0", "-1e-8"):
        with pytest.raises(SystemExit):
            sg.main(["--ltr_fasta", str(d / "ltr.fa"), "--annotation", str(d / "annotation.gff3"), "--outdir",
                     str(tmp_path / "o"), "--mu", mu])


@pytest.mark.slow
def test_assign_conflicts_fail_fast(tmp_path):
    """--assign with -n, or a file naming no sequence, stops with a clear message."""
    d, truth = sg.simulate(tmp_path, tree={"split": 2e6, "children": ["A", "B"]},
                         mergers=[{"time": 1e5, "lineages": ["A", "B"]}])
    (tmp_path / "part.tsv").write_text(f"{sorted(truth)[0]}\tX\n")
    base = ["--ltr_fasta", str(d / "ltr.fa"), "--annotation", str(d / "annotation.gff3"), "--outdir",
            str(tmp_path / "o"), "--assign", str(tmp_path / "part.tsv")]
    with pytest.raises(SystemExit):
        sg.main(base + ["-n", "2"])
    (tmp_path / "empty.tsv").write_text("# nothing\n")
    with pytest.raises(SystemExit):
        sg.main(base[:-1] + [str(tmp_path / "empty.tsv")])


@pytest.mark.slow
def test_e2e_allotetraploid(tmp_path):
    d, truth = sg.simulate(tmp_path, tree={"split": 2e6, "children": ["A", "B"]},
                         mergers=[{"time": 1e5, "lineages": ["A", "B"]}])
    rows, s = sg.run_tool(d, tmp_path / "out")
    assert s["n_subgenomes"] == 2 and sg.wrong(rows, truth) == 0
    assert all(r["status"] == "assigned" for r in rows)
    assert json.loads((tmp_path / "out" / "summary.json").read_text())["version"] == __version__


@pytest.mark.slow
def test_e2e_diploid_n1(tmp_path):
    d, truth = sg.simulate(tmp_path, tree="A", mergers=[])
    rows, s = sg.run_tool(d, tmp_path / "out")
    assert s["n_subgenomes"] == 1 and all(r["status"] == "no_structure" for r in rows)


@pytest.mark.slow
def test_e2e_unphased_scaffolds(tmp_path):
    d, truth = sg.simulate(tmp_path, tree={"split": 2e6, "children": ["A", "B"]},
                         mergers=[{"time": 1e5, "lineages": ["A", "B"]}])
    fa = (d / "ltr.fa").read_text()
    seq = "ACGT" * 2000
    (d / "ltr.fa").write_text(fa + f">scaf9:1-8000#LTR/Gypsy/Tekay\n{seq}\n")
    (d / "annotation.gff3").write_text((d / "annotation.gff3").read_text() +
                                       "scaf9\ts\tLTR_retrotransposon\t1\t8000\t.\t+\t.\tID=x;Name=scaf9:1-8000;K2P_d=0.01\n")
    (d / "genome.fai").write_text((d / "genome.fai").read_text() + "scaf9\t9000\t0\t0\t0\n")
    rows, s = sg.run_tool(d, tmp_path / "out")
    sc = [r for r in rows if r["sequence"] == "scaf9"][0]
    assert sc["status"] == "unphased" and sc["subgenome"] == "NA" and s["n_subgenomes"] == 2


@pytest.mark.slow
def test_forced_n_reports_the_subgenomes_found(tmp_path):
    d, truth = sg.simulate(tmp_path, tree={"split": 2e6, "children": ["A", "B"]},
                         mergers=[{"time": 1e5, "lineages": ["A", "B"]}])
    rows, s = sg.run_tool(d, tmp_path / "out", "-n", "8")
    assert all(v > 0 for v in s["subgenome_sizes"].values())
    assert s["n_subgenomes"] == len(s["subgenome_sizes"]) and s["n_subgenomes_requested"] == 8


@pytest.mark.slow
def test_links_oe_is_observed_over_expected(tmp_path):
    d, truth = sg.simulate(tmp_path, tree={"split": 2e6, "children": ["A", "B"]},
                         mergers=[{"time": 1e5, "lineages": ["A", "B"]}])
    sg.run_tool(d, tmp_path / "out")
    rows = [l.split("\t") for l in (tmp_path / "out" / "links_oe.tsv").read_text().splitlines()[1:]]
    assert all(r[1 + i] == "NA" for i, r in enumerate(rows)) and float(rows[0][2]) > 0


def _sim_pairs(tmp_path, split):
    d, truth = sg.simulate(tmp_path, tree={"split": split, "children": ["A", "B"]},
                         mergers=[{"time": min(1e4, split / 5), "lineages": ["A", "B"]}])
    rows = [l.split("\t") for l in (tmp_path / "truth" / "chromosomes.tsv").read_text().splitlines()[1:]]
    hom = {r[0]: r[2] for r in rows}
    return d, truth, hom


@pytest.mark.slow
def test_e2e_verdicts_and_outputs(tmp_path):
    """Verdict, pairing map, and the structure records: structure.tsv holds summary.json's tree rows (the reported
    candidate and every link source's numbers); summary.json gives the LTR-K2P tercile edges and copies, no 'rule'."""
    d, truth = sg.simulate(tmp_path, tree={"split": 2e6, "children": ["A", "B"]},
                         mergers=[{"time": 1e5, "lineages": ["A", "B"]}])
    rows, s = sg.run_tool(d, tmp_path / "out")
    assert s["verdict"] == "subgenomes found" and "n_informative" in rows[0]
    assert (tmp_path / "out" / "pairing_map.tsv").read_text().startswith("sequence_a\tsequence_b\t")
    st = (tmp_path / "out" / "structure.tsv").read_text().splitlines()[0].split("\t")
    src = ("all",) + cst.AGE_STRATA
    per_source = {f"{k}_{x}" for x in src for k in ("replicability", "strength", "agrees", "parent_strength")}
    assert {"copies", "replicability_se", "n_pairs", "borderline", "forced"} | per_source <= set(st)
    assert set(st) == set(s["structure_tree"][0]) | {"n_a", "n_b"}
    assert "rule" not in s and len(s["strata_k2p_edges"]) == 2 and set(s["strata_copies"]) == set(src)
    assert sum(s["strata_copies"][x] for x in cst.AGE_STRATA) == s["strata_copies"]["all"]


@pytest.mark.slow
def test_e2e_pairs_split_close_subgenomes(tmp_path):
    """Subgenomes split 300 kyr ago: lineage links alone give N=1 (the root lineage split is rejected); pairing the
    homologous chromosomes (shared insertions at matching positions) and separating each pair, on by default,
    recovers both subgenomes."""
    d, truth, hom = _sim_pairs(tmp_path, 3e5)
    rows, s = sg.run_tool(d, tmp_path / "out")
    assert not any(t["accepted"] for t in s["structure_tree"] if t["rule"] != "pairs")
    assert any(t["accepted"] for t in s["structure_tree"] if t["rule"] == "pairs")
    assert s["n_subgenomes"] == 2 and sg.wrong(rows, truth) == 0


@pytest.mark.slow
def test_links_oe_unmasked_after_pair_split(tmp_path):
    """links_oe.tsv reports the observed links: homologues stay the most linked pairs after the pair split."""
    d, truth, hom = _sim_pairs(tmp_path, 3e5)
    sg.run_tool(d, tmp_path / "out")
    lines = (tmp_path / "out" / "links_oe.tsv").read_text().splitlines()
    names = lines[0].split("\t")[1:]
    M = {r.split("\t")[0]: dict(zip(names, r.split("\t")[1:])) for r in lines[1:]}
    vals = [float(M[x][y]) for x in names for y in names if x != y and hom[x] == hom[y]]
    assert vals and min(vals) > 1


@pytest.mark.slow
def test_forced_n_applies_pairing(tmp_path):
    """-n 2: paired homologues are kept apart as when N is inferred, so the verdict (computed on the final groups)
    reports no paired sequences sharing a group."""
    d, truth, hom = _sim_pairs(tmp_path, 3e5)
    rows, s = sg.run_tool(d, tmp_path / "out", "-n", "2")
    assert s["n_subgenomes"] == 2 and sg.wrong(rows, truth) == 0 and "share a group" not in s["verdict"]


@pytest.mark.slow
def test_e2e_dating_allotetraploid(tmp_path):
    """Split 2 Myr, merger 100 kyr (mu 1.3e-8: K2P 0.052 and 0.0026): the reported split is within 10% and the merger
    interval (or bound) contains the truth."""
    d, truth = sg.simulate(tmp_path, tree={"split": 2e6, "children": ["A", "B"]},
                         mergers=[{"time": 1e5, "lineages": ["A", "B"]}])
    rows, s = sg.run_tool(d, tmp_path / "out", "--mu", "1.3e-8", "-t", "8")
    D = sg.dating_rows(tmp_path / "out")
    sp = [r for r in D if r["event"] == "split" and r["reported"] == "yes"]
    assert sp and abs(float(sp[0]["years"]) / 2e6 - 1) < 0.1
    mg = [r for r in D if r["event"] == "merger" and r["reported"] == "yes"][0]
    assert mg["status"] in ("dated", "bound")
    assert float(mg["years_low"]) <= 1e5 <= float(mg["years_high"])
    sd = s["dating"]["pairs"][0]
    assert sd["merger"]["status"] == mg["status"] and abs(sd["split"]["years"] / 2e6 - 1) < 0.1


@pytest.mark.slow
def test_e2e_no_dating_for_one_subgenome(tmp_path):
    d, truth = sg.simulate(tmp_path, tree="A", mergers=[])
    rows, s = sg.run_tool(d, tmp_path / "out")
    assert s["n_subgenomes"] == 1 and not (tmp_path / "out" / "dating.tsv").exists()


@pytest.mark.slow
def test_e2e_biology_outputs(tmp_path):
    """Allotetraploid (split 2 Myr, merger 100 kyr): clade table with copies after the merger in both subgenomes when
    the merger is dated (a bound is reported for this small LTR-poor genome: its genealogy-clock interval reaches below
    one substitution per LTR; then no copy can be called post-merger), insertion history per subgenome (observed copies
    sum to the subgenome's copies), nesting check."""
    d, truth = sg.simulate(tmp_path, tree={"split": 2e6, "children": ["A", "B"]},
                         mergers=[{"time": 1e5, "lineages": ["A", "B"]}])
    rows, s = sg.run_tool(d, tmp_path / "out", "--mu", "1.3e-8", "-t", "8")
    rd = lambda f: [l.split("\t") for l in (tmp_path / "out" / f).read_text().splitlines() if not l.startswith("#")]
    C = rd("clades.tsv")
    after = {r[1]: 0.0 for r in C[1:]}
    for r in C[1:]:
        after[r[1]] += float(r[8]) if r[8] != "NA" else 0.0
    sgs = {r["sequence"]: r["subgenome"] for r in rows}
    T = [l.split("\t") for l in (tmp_path / "truth" / "copies.tsv").read_text().splitlines()]
    h = T[0]
    true_after = {g: sum(1 for r in T[1:] if sgs.get(r[h.index("chrom")]) == g and r[h.index("post_merger")] == "yes")
                  for g in ("SG1", "SG2")}
    assert set(after) == {"SG1", "SG2"}
    if s["dating"]["pairs"][0]["merger"]["status"] == "dated":
        for g in ("SG1", "SG2"):                            # expected count from each copy's LTR counts
            assert 0.75 * true_after[g] <= after[g] <= 1.33 * true_after[g], (g, after[g], true_after[g])
    else:
        assert all(v == 0 for v in after.values())
    H = rd("insertion_history.tsv")
    for g in ("SG1", "SG2"):
        obs = sum(int(r[6]) for r in H[1:] if r[0] == g and r[1] == "all")
        assert obs == sum(int(r[2]) for r in C[1:] if r[1] == g)            # every simulated copy has a K2P
    assert s["biology"]["nesting"]["n_nested"] == len(rd("nesting.tsv")) - 1


@pytest.mark.slow
def test_e2e_biology_for_one_subgenome(tmp_path):
    """One subgenome: the clade table, insertion history and nesting check are written; no merger, so no copies
    'after the merger'."""
    d, truth = sg.simulate(tmp_path, tree="A", mergers=[])
    rows, s = sg.run_tool(d, tmp_path / "out")
    assert s["n_subgenomes"] == 1 and s["biology"]["nesting"]["n_nested"] >= 0
    assert len((tmp_path / "out" / "insertion_history.tsv").read_text().splitlines()) > 10
    C = [l.split("\t") for l in (tmp_path / "out" / "clades.tsv").read_text().splitlines()]
    assert len(C) > 2 and all(r[8] == "NA" for r in C[1:])


@pytest.mark.slow
def test_e2e_painting_finds_a_simulated_exchange(tmp_path):
    """A 1-Mb non-reciprocal exchange (A -> B, 50 kyr after a 100-kyr merger; a quarter of a 4-Mb chromosome) is
    called on the right chromosome with the right origin; the painting outputs are written."""
    ex = [dict(type="nonreciprocal", time=5e4, recipient="B", donor="A", chrom=1, start=1.0e6, end=2.0e6)]
    d, truth = sg.simulate(tmp_path, tree={"split": 2e6, "children": ["A", "B"]},
                         mergers=[{"time": 1e5, "lineages": ["A", "B"]}], exchange_list=ex)
    rows, s = sg.run_tool(d, tmp_path / "out", "-t", "4")
    T = [l.split("\t") for l in (tmp_path / "truth" / "exchanges.tsv").read_text().splitlines()]
    h = T[0]
    te = [dict(zip(h, r)) for r in T[1:]]
    g2t = {}
    for r in rows:
        if r["subgenome"] != "NA":
            g2t.setdefault(r["subgenome"], truth[r["sequence"]])
    S = [l.split("\t") for l in (tmp_path / "out" / "segments.tsv").read_text().splitlines()]
    called = [dict(zip(S[0], r)) for r in S[1:] if r[-1] == "yes"]
    hit = [c for c in called for e in te if c["sequence"] == e["chrom"] and g2t[c["origin"]] == e["donor"]
           and int(c["start"]) < int(float(e["end"])) and int(c["end"]) > int(float(e["start"]))]
    assert hit and s["painting"]["n_called"] >= 1
    for f in ("exchanges.bed", "exchanges.gff3", "resolution.bedGraph", "ancestry_SG1.bedGraph"):
        assert (tmp_path / "out" / f).exists()


@pytest.mark.slow
def test_e2e_elements_table(tmp_path):
    """elements.tsv: one row per LTR-RT with its sequence's subgenome, ancestry posteriors and a call; among copies
    with decisive evidence, those inside a planted 1-Mb exchange are called foreign, those elsewhere rarely (copies
    that moved after the merger)."""
    ex = [dict(type="nonreciprocal", time=5e4, recipient="B", donor="A", chrom=1, start=1.0e6, end=2.0e6)]
    d, truth = sg.simulate(tmp_path, tree={"split": 2e6, "children": ["A", "B"]},
                         mergers=[{"time": 1e5, "lineages": ["A", "B"]}], exchange_list=ex)
    rows, s = sg.run_tool(d, tmp_path / "out", "-t", "8", "--mu", "1.3e-8")
    E = [l.split("\t") for l in (tmp_path / "out" / "elements.tsv").read_text().splitlines()]
    h, R = E[0], [dict(zip(E[0], r)) for r in E[1:]]
    n_fa = sum(1 for l in (d / "ltr.fa").read_text().splitlines() if l.startswith(">"))
    assert len(R) == n_fa
    assert {"element", "sequence", "start", "end", "subgenome", "clade", "k2p", "age_years", "informative", "P_SG1",
            "P_SG2", "call", "favored", "p_after_merger"} <= set(h)
    T = [l.split("\t") for l in (tmp_path / "truth" / "exchanges.tsv").read_text().splitlines()]
    e = dict(zip(T[0], T[1]))
    inf = [r for r in R if r["informative"] == "yes" and r["subgenome"] != "NA"]
    inside = [r for r in inf if r["sequence"] == e["chrom"] and float(e["start"]) <= int(r["start"])
              and int(r["end"]) <= float(e["end"])]
    outside = [r for r in inf if r["sequence"] != e["chrom"]]
    share = lambda rs: sum(r["call"] == "foreign" for r in rs) / max(1, sum(r["call"] in ("own", "foreign") for r in rs))
    assert share(inside) >= 0.8 and share(outside) <= share(inside) / 4    # among decisive calls ('weak' = undecided)
    assert all(r["favored"] != r["subgenome"] for r in inf if r["call"] == "foreign")


@pytest.mark.slow
def test_assign_fixes_the_partition(tmp_path, capsys):
    """--assign: the given partition is used as is (labels become SG1, SG2 in the file's order, the mapping kept in
    summary.json), the structure search is skipped, sequences left out of the file are unphased, outputs say 'user'."""
    d, truth = sg.simulate(tmp_path, tree={"split": 2e6, "children": ["A", "B"]},
                         mergers=[{"time": 1e5, "lineages": ["A", "B"]}])
    seqs = sorted(truth, key=inputs.natural_key)
    keep = seqs[1:]                                                  # the first sequence left out: unphased
    (tmp_path / "part.tsv").write_text("# user partition\n" + "".join(  # file order: last sequence first
        f"{s}\t{'X' if truth[s] == truth[keep[0]] else 'Y'}\n" for s in keep[::-1]) + "nosuchseq\tX\n")
    rows, s = sg.run_tool(d, tmp_path / "out", "--assign", str(tmp_path / "part.tsv"))
    lab = {r["sequence"]: r for r in rows}
    assert s["n_subgenomes_source"] == "user" and s["n_subgenomes"] == 2 and s["structure_tree"] == []
    assert lab[seqs[0]]["status"] == "unphased" and all(lab[x]["status"] == "user" for x in keep)
    assert all((lab[x]["subgenome"] == lab[y]["subgenome"]) == (truth[x] == truth[y]) for x in keep for y in keep)
    first = "X" if truth[keep[-1]] == truth[keep[0]] else "Y"            # the file's first label is SG1
    assert lab[keep[-1]]["subgenome"] == "SG1" and s["user_labels"] == {"SG1": first, "SG2": "XY".replace(first, "")}
    assert "nosuchseq" in capsys.readouterr().err


@pytest.mark.slow
def test_unphased_sequence_ancestry(tmp_path):
    """A scaffold too poor in LTR-RTs to phase but carrying copies of one subgenome's lineages is called that
    subgenome; a scaffold of unrelated sequence gets no call."""
    d, truth = sg.simulate(tmp_path, tree={"split": 2e6, "children": ["A", "B"]},
                         mergers=[{"time": 1e5, "lineages": ["A", "B"]}])
    recs, name, seq = [], None, []
    for line in (d / "ltr.fa").read_text().splitlines():
        if line.startswith(">"):
            if name:
                recs.append((name, "".join(seq)))
            name, seq = line[1:], []
        else:
            seq.append(line)
    recs.append((name, "".join(seq)))
    src = sorted(truth, key=inputs.natural_key)[0]
    picks = [r for r in recs if r[0].split(":")[0] == src][::7][:6]
    rng = np.random.default_rng(3)
    fa, gff, pos = [], [], 1
    for k, (nm, sq) in enumerate(picks):
        fa.append(f">scaf9:{pos}-{pos + len(sq) - 1}#LTR/Gypsy/Tekay\n{sq}\n")
        gff.append(f"scaf9\ts\tLTR_retrotransposon\t{pos}\t{pos + len(sq) - 1}\t.\t+\t.\tID=u{k};K2P_d=0.01\n")
        pos += len(sq) + 1000
    for k in range(3):
        sq = "".join(rng.choice(list("ACGT"), 5000))
        fa.append(f">scaf8:{1 + k * 6000}-{5000 + k * 6000}#LTR/Gypsy/Tekay\n{sq}\n")
        gff.append(f"scaf8\ts\tLTR_retrotransposon\t{1 + k * 6000}\t{5000 + k * 6000}\t.\t+\t.\tID=r{k};K2P_d=0.01\n")
    (d / "ltr.fa").write_text((d / "ltr.fa").read_text() + "".join(fa))
    (d / "annotation.gff3").write_text((d / "annotation.gff3").read_text() + "".join(gff))
    (d / "genome.fai").write_text((d / "genome.fai").read_text() + f"scaf9\t{pos + 1000}\t0\t0\t0\nscaf8\t20000\t0\t0\t0\n")
    rows, s = sg.run_tool(d, tmp_path / "out")
    sgs = {r["sequence"]: r["subgenome"] for r in rows}
    U = [l.split("\t") for l in (tmp_path / "out" / "unphased_sequences.tsv").read_text().splitlines()]
    u = {r[0]: dict(zip(U[0], r)) for r in U[1:]}
    assert u["scaf9"]["ancestry"] == sgs[src] and u["scaf9"]["confidence"] in ("decisive", "substantial")
    assert u["scaf8"]["confidence"] == "none" and u["scaf8"]["ancestry"] == "NA"


@pytest.mark.slow
def test_biology_draws_history_before_after_merger(tmp_path):
    """Reproducibility across versions: the insertion history is drawn first from the biology generator, as before
    elements.tsv existed (whether or not the merger is dated)."""
    d, truth = sg.simulate(tmp_path, tree={"split": 2e6, "children": ["A", "B"]},
                         mergers=[{"time": 1e5, "lineages": ["A", "B"]}])
    a = sg.parse_args(["--ltr_fasta", str(d / "ltr.fa"), "--annotation", str(d / "annotation.gff3"), "--outdir",
                       str(tmp_path / "o"), "-t", "8", "--mu", "1.3e-8"])
    run.setup(a)
    lib = inputs.load_library(a)
    ph = pipeline.run_phasing(a, lib, lineages.load_or_build_clades(a, lib))
    dt = dating.run_dating(a, lib, ph)
    bio = biology.run_biology(a, lib, ph, dt)
    loc = np.where(ph["ek"] >= 0, ph["z"][np.maximum(ph["ek"], 0)], -1)
    ref = biology.insertion_history(lib, loc, ph["N"], np.random.default_rng([a.seed, 11]))
    assert [r["deconvolved"] for r in bio["history"]] == [r["deconvolved"] for r in ref]


HEX = {"root_age": 8e6, "chromosomes": {"n": 5, "length": 6e6}, "families": {"n": 80},
       "activity": {"rate": 0.001, "p_burst": 0.0}, "exchanges": {"nonreciprocal": 0, "reciprocal": 0},
       "tree": {"split": 4.0e6, "children": [{"split": 1.0e6, "children": ["A", "B"]}, "C"]},
       "mergers": [{"time": 2.0e4, "lineages": ["A", "B", "C"]}]}


@pytest.mark.slow
def test_e2e_recent_split_found_in_young_copies(tmp_path):
    """Hexaploid: A and B split 1 Myr ago, C 4 Myr ago, all merged 20 kyr ago. On pooled copies A|B is reproducible
    but under half as strong as C|AB (the pooled rule stops at 2 groups); in young copies it is as strong as its
    parent, and old copies (older than the A/B split) show no reproducible split: the three subgenomes are found, every
    sequence right."""
    builders.simulate(HEX, tmp_path / "in", tmp_path / "truth", seed=1)
    d, out = tmp_path / "in", tmp_path / "out"
    sg.main(["--ltr_fasta", str(d / "ltr.fa"), "--annotation", str(d / "annotation.gff3"), "--fai",
             str(d / "genome.fai"), "--outdir", str(out), "-t", "4"])
    s = json.loads((out / "summary.json").read_text())
    assert s["n_subgenomes"] == 3
    sub = [t for t in s["structure_tree"] if t["level"] == "root.a"][0]
    assert sub["accepted"] and sub["rule"] == "age" and sub["copies"] == "young"
    assert sub["strength_all"] < cst.MIN_REL * sub["parent_strength_all"]               # the pooled rule rejects it
    truth = {l.split("\t")[0]: l.split("\t")[1] for l in (tmp_path / "truth" / "chromosomes.tsv").read_text().splitlines()[1:]}
    rows = [l.split("\t") for l in (out / "chromosomes.tsv").read_text().splitlines()[1:]]
    pairs = {(r[4], truth[r[0]]) for r in rows}
    assert len(pairs) == 3 and len({p[0] for p in pairs}) == 3                           # one label per subgenome
    assert (out / "report.pdf").exists()


def _tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


NEW_FILES = {"homoeologs.tsv", "homoeolog_sets.tsv", "homoeolog_anchors.tsv", "homoeolog_exchanges.tsv",
             "homoeolog_halflife.tsv", "homoeolog_rearrangements.tsv", "findings.txt"}


def _same_tables(a, b):
    """Every table of two output dirs equal: summary.json without version, timings and homoeologs; exchanges.gff3
    without its source column; the homoeolog files, report.pdf and report.png aside (the page names its run folder),
    and manifest.json (it names the run's folder, options and times)."""
    skip = {"report.pdf", "report.png", "stderr.txt", "manifest.json"} | NEW_FILES
    fa = sorted(p.name for p in a.iterdir() if p.is_file() and p.name not in skip)
    fb = sorted(p.name for p in b.iterdir() if p.is_file() and p.name not in skip)
    assert fa == fb
    for f in fa:
        if f == "summary.json":
            ja, jb = json.loads((a / f).read_text()), json.loads((b / f).read_text())
            for k in ("version", "timings", "homoeologs"):
                ja.pop(k, None), jb.pop(k, None)
            assert ja == jb, f
        elif f == "exchanges.gff3":
            cut = lambda p: [l.split("\t")[:1] + l.split("\t")[2:] for l in p.read_text().splitlines()]
            assert cut(a / f) == cut(b / f), f
        else:
            assert (a / f).read_bytes() == (b / f).read_bytes(), f


@pytest.mark.slow
def test_e2e_pairs_right_and_tables_unchanged_by_verbose(tmp_path, capsys):
    """Allotetraploid, 2-Myr split: every chromosome paired with its true homoeolog, reciprocally, in four complete
    sets; every other table as a run without -v writes it."""
    d, truth = sg.simulate(tmp_path, tree={"split": 2e6, "children": ["A", "B"]},
                           mergers=[{"time": 1e5, "lineages": ["A", "B"]}])
    sg.run_tool(d, tmp_path / "plain", cache=False)                   # same inputs, same paths in summary.json
    capsys.readouterr()
    sg.run_tool(d, tmp_path / "verbose", "-v")
    err = capsys.readouterr().err
    assert "old copies" not in err                                    # old mask = the phasing step's old tercile
    _same_tables(tmp_path / "plain", tmp_path / "verbose")
    hom = {r["token"]: r["homoeolog"] for r in _tsv(tmp_path / "truth" / "chromosomes.tsv")}
    rows = _tsv(tmp_path / "verbose" / "homoeologs.tsv")
    assert len(rows) == 8
    for r in rows:
        assert r["status"] == "homoeolog" and r["reciprocal"] == "yes"
        assert hom[r["homoeolog"]] == hom[r["sequence"]]
    sets = _tsv(tmp_path / "verbose" / "homoeolog_sets.tsv")
    assert len(sets) == 4 and all(s["complete"] == "yes" for s in sets)
    anchors = _tsv(tmp_path / "verbose" / "homoeolog_anchors.tsv")
    assert anchors and all(a["subgenome_a"] != a["subgenome_b"] for a in anchors)
    s = json.loads((tmp_path / "verbose" / "summary.json").read_text())
    assert s["version"] == __version__
    assert s["homoeologs"]["pairs"][0]["reciprocal"] == 4
    assert not _tsv(tmp_path / "verbose" / "homoeolog_exchanges.tsv")       # no exchange simulated: no stretch
    assert (tmp_path / "verbose" / "findings.txt").read_text().startswith("What the homoeologous chromosome pairs say")
    assert {r["lost_in"] for r in _tsv(tmp_path / "verbose" / "homoeolog_halflife.tsv")} == {"SG1", "SG2"}
    pdf = tmp_path / "verbose" / "report.pdf"                               # one page; B is the homoeolog dot plot
    if shutil.which("pdftotext") and shutil.which("pdfinfo"):
        info = subprocess.run(["pdfinfo", str(pdf)], capture_output=True, text=True).stdout
        assert re.search(r"^Pages:\s+1$", info, re.M)
        txt = " ".join(subprocess.run(["pdftotext", str(pdf), "-"], capture_output=True, text=True).stdout.split())
        assert "Old insertions shared by" in txt and "homoeolog set" in txt and "pair called" in txt
        words = re.findall(r">([^<]*)</word>", subprocess.run(["pdftotext", "-bbox", str(pdf), "-"],
                                                              capture_output=True, text=True).stdout)
        assert set("ABCDEFGH") <= set(words)                          # every panel letter, as its own word


@pytest.mark.slow
def test_e2e_planted_replacement_is_found_where_it_is(tmp_path):
    """One replacing exchange (1 Mb of chromosome 1, B <- A, 50 kyr ago) -> one near-identical stretch on that
    homoeolog pair, inside the exchanged interval; nothing elsewhere."""
    d, truth = sg.simulate(tmp_path, seed=3, tree={"split": 2e6, "children": ["A", "B"]},
                           mergers=[{"time": 1e5, "lineages": ["A", "B"]}],
                           activity=dict(merge_defaults(sg.SIM_TINY)["activity"], rate=0.002),
                           exchange_list=[dict(type="nonreciprocal", time=5e4, recipient="B", donor="A", chrom=1,
                                               start=1.5e6, end=2.5e6)])
    sg.run_tool(d, tmp_path / "out")
    hom = {r["token"]: r["homoeolog"] for r in _tsv(tmp_path / "truth" / "chromosomes.tsv")}
    st = _tsv(tmp_path / "out" / "homoeolog_exchanges.tsv")
    assert len(st) == 1
    x = st[0]
    assert hom[x["sequence_a"]] == hom[x["sequence_b"]] == "1"
    for k in ("a", "b"):                                               # a copy moves when its start is inside
        assert 1.5e6 <= int(x[f"start_{k}"]) and int(x[f"end_{k}"]) <= 2.5e6 + 2e4
    assert "Exchanges: 1 stretch" in (tmp_path / "out" / "findings.txt").read_text()


@pytest.mark.slow
def test_e2e_old_split_makes_no_wrong_call(tmp_path):
    """15-Myr split (half-life 1.5 Myr): inherited insertions are gone; whatever is called must be right."""
    d, truth = sg.simulate(tmp_path, tree={"split": 1.5e7, "children": ["A", "B"]}, root_age=3e7,
                           mergers=[{"time": 1e5, "lineages": ["A", "B"]}])
    rows = _tsv(tmp_path / "truth" / "chromosomes.tsv")
    (tmp_path / "assign.txt").write_text("".join(f"{r['token']} {r['subgenome']}\n" for r in rows))
    sg.run_tool(d, tmp_path / "out", "--assign", str(tmp_path / "assign.txt"))
    hom = {r["token"]: r["homoeolog"] for r in rows}
    out = _tsv(tmp_path / "out" / "homoeologs.tsv")
    assert len(out) == 8
    assert all(hom[r["homoeolog"]] == hom[r["sequence"]] for r in out if r["status"] == "homoeolog")


@pytest.mark.slow
def test_e2e_one_group_writes_no_homoeologs(tmp_path):
    d, _ = sg.simulate(tmp_path, tree="A", mergers=[])
    sg.run_tool(d, tmp_path / "out")
    assert not any((tmp_path / "out" / f).exists() for f in NEW_FILES)
    assert json.loads((tmp_path / "out" / "summary.json").read_text())["homoeologs"] is None


@pytest.mark.slow
def test_e2e_failed_homoeologs_cost_nothing_else(tmp_path, monkeypatch, capsys):
    d, _ = sg.simulate(tmp_path, tree={"split": 2e6, "children": ["A", "B"]},
                       mergers=[{"time": 1e5, "lineages": ["A", "B"]}])

    def boom(*a, **k):
        raise RuntimeError("planted failure")
    monkeypatch.setattr(run, "run_homoeologs", boom)
    sg.main(["--ltr_fasta", str(d / "ltr.fa"), "--annotation", str(d / "annotation.gff3"), "--fai",
             str(d / "genome.fai"), "--outdir", str(tmp_path / "out"), "-t", "4"])
    assert "homoeolog pairing failed (RuntimeError: planted failure)" in capsys.readouterr().err
    out = tmp_path / "out"
    assert (out / "elements.tsv").exists() and (out / "chromosomes.tsv").exists() and (out / "report.pdf").exists()
    assert not any((out / f).exists() for f in NEW_FILES)


SMALL = {"root_age": 8e6, "chromosomes": {"n": 5, "length": 6e6}, "families": {"n": 80},
         "activity": {"rate": 0.002, "p_burst": 0.0}, "tree": {"split": 2.0e6, "children": ["A", "B"]},
         "mergers": [{"time": 1.0e5, "lineages": ["A", "B"]}], "exchanges": {"nonreciprocal": 0, "reciprocal": 0}}


def _run_k2p(d, out):
    """A run on the K2P-table annotation (k2p.tsv): ({sequence: chromosomes.tsv row}, summary)."""
    rows, s = sg.run_tool(d, out, "--seed", "1", annotation="k2p.tsv")
    return {x["sequence"]: x for x in rows}, s


@pytest.fixture(scope="module")
def tetraploid(tmp_path_factory):
    base = tmp_path_factory.mktemp("sim4x")
    out, truth = builders.simulate(SMALL, base / "sim", base / "truth", seed=3)
    planted, parent, decoys = builders.plant_scaffolds(out, base / "planted", seed=3)
    return out, truth, planted, parent, decoys


@pytest.mark.slow
def test_scaffolds_and_decoy_compartment(tetraploid, tmp_path):
    """An allotetraploid whose chromosomes lost their ends to 30 unplaced scaffolds, plus a decoy compartment of 8
    small sequences sharing a private subfamily: two subgenomes, every chromosome right, decoys never a group."""
    out, truth, planted, parent, decoys = tetraploid
    rows, js = _run_k2p(planted, tmp_path / "out")
    tsg = builders.truth_subgenomes(truth)
    assert js["n_subgenomes"] == 2
    chroms = [c for c in tsg]
    pred = {c: rows[c]["subgenome"] for c in chroms}
    m = {}
    for c in chroms:                                      # predicted label -> true label, must be one-to-one
        m.setdefault(pred[c], set()).add(tsg[c])
    assert all(len(v) == 1 for v in m.values()) and len(m) == 2
    assert all(rows[c]["core"] == "yes" for c in chroms)
    assert all(rows[d]["core"] == "no" for d in decoys if rows[d]["status"] != "unphased")
    placed = [s for s in parent if rows.get(s, {}).get("status") == "assigned"]
    for s in placed:                                      # confidently placed scaffolds follow their chromosome
        assert rows[s]["subgenome"] == rows[parent[s]]["subgenome"], s
    assert len(placed) >= len(parent) // 3


@pytest.mark.slow
def test_diploid_has_no_structure(tmp_path):
    """A diploid (one lineage, no merger): no structure."""
    sc = dict(SMALL, tree="A", mergers=[])
    out, truth = builders.simulate(sc, tmp_path / "sim", tmp_path / "truth", seed=5)
    rows, js = _run_k2p(out, tmp_path / "out")
    assert js["n_subgenomes"] == 1


@pytest.mark.slow
def test_full_run_reports_inherited_split(tetraploid, tmp_path):
    """Full dating on an allotetraploid: the inherited-insertion split is reported and close to the truth."""
    out, truth, _, _, _ = tetraploid
    rows, js = _run_k2p(out, tmp_path / "full")
    dt = sg.dating_rows(tmp_path / "full")
    rep = [r for r in dt if r["event"] == "split" and r["reported"] == "yes"]
    assert rep
    inh = [r for r in dt if r["event"] == "split" and r["estimator"] == "inherited"]
    assert inh and inh[0]["status"] == "dated"
    truth_k2p = 2 * 1.3e-8 * 2.0e6
    est = float(inh[0]["k2p"])
    assert abs(est / truth_k2p - 1) < 0.15, (est, truth_k2p)
    assert rep[0]["estimator"] == "inherited"


@pytest.mark.slow
@pytest.mark.skipif(not os.environ.get("SYNLTR_LONG_TESTS"), reason="set SYNLTR_LONG_TESTS=1 for this simulation")
def test_young_merger_end_to_end(tmp_path):
    """Two subgenomes split 2 Myr, merged 10 kyr ago (40-Mb chromosomes x 4): the reported merger interval contains
    the truth (2 mu T = 2.6e-4) or is a bound above it; the LTR-only clock alone cannot date it."""
    sc = {"mergers": [{"time": 1.0e4, "lineages": ["A", "B"]}], "chromosomes": {"n": 4, "length": 4.0e7}}
    d, _ = builders.simulate(sc, tmp_path / "data", tmp_path / "truth", seed=11)
    out = tmp_path / "out"
    sg.main(["--ltr_fasta", str(d / "ltr.fa"), "--annotation", str(d / "annotation.gff3"), "--fai",
             str(d / "genome.fai"), "--outdir", str(out), "-t", "4"])
    rows = sg.dating_rows(out)
    rep = [r for r in rows if r["event"] == "merger" and r["reported"] == "yes"]
    assert len(rep) == 1
    lo = float(rep[0]["k2p_low"]) if rep[0]["k2p_low"] != "NA" else None
    hi = float(rep[0]["k2p_high"]) if rep[0]["k2p_high"] != "NA" else None
    assert rep[0]["status"] in ("dated", "bound") and hi is not None and (lo or 0) <= 2.6e-4 <= hi
