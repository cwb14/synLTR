"""The phasing pipeline: verdicts, the --assign file, homoeologs across a candidate split."""
import builders
import numpy as np
import pytest
from analyses.subgenomes import sgtools as sg

from ltrtrace.analyses.subgenomes import inputs, lineages, pipeline, run
from ltrtrace.base import constants as cst
from ltrtrace.base.log import UsageError


def test_verdict_names_near_identical_sequences():
    """N=1 with paired sequences or larger groups of near-identical sequences: the verdict says what they may be
    (homologues of a young polyploid or both haplotypes of a diploid assembly), never 'no structure detected'; splits
    forced by -n are counted."""
    t = [dict(cands={"all": dict(borderline=False)}, chosen=None, accepted=False, forced=False)]
    for inside, big in ((3, []), (0, [[0, 1, 2]])):
        v = pipeline.verdict_text(1, inside, big, t, "inferred")
        assert v.startswith("cannot tell") and "near-identical" in v and "haplotypes" in v
    assert pipeline.verdict_text(1, 0, [], t, "inferred") == "no structure detected"
    assert pipeline.verdict_text(2, 0, [], t, "inferred") == "subgenomes found"
    assert "1 of the requested splits fail the split rule" in pipeline.verdict_text(2, 0, [], [dict(t[0], forced=True)], "given")


def test_read_assign_validates_the_file(tmp_path, caplog):
    """Conflicting labels for one sequence, no sequence of the library, or more groups than MAX_N are usage errors
    naming --assign and the file; one-token lines are counted and skipped; a single group is warned about."""
    import logging

    names = [f"chr{i}" for i in range(1, 21)]
    f = tmp_path / "a.tsv"
    f.write_text("chr1 A\nchr2 B\nchr1 B\n")
    with pytest.raises(UsageError, match=f"--assign {f}: .*chr1"):
        pipeline.read_assign(f, names)
    f.write_text("scaffold_9 A\nscaffold_8 B\n")
    with pytest.raises(UsageError, match=f"--assign {f}: .*scaffold_9"):
        pipeline.read_assign(f, names)
    f.write_text("".join(f"chr{i} G{i}\n" for i in range(1, 15)))
    with pytest.raises(UsageError, match=f"--assign {f}: 14 groups; at most {cst.MAX_N}"):
        pipeline.read_assign(f, names)
    caplog.set_level(logging.WARNING, logger="ltrtrace")
    f.write_text("chr1 A\nchr2\nchr3 B\nchr4\n")
    assert pipeline.read_assign(f, names) == {"chr1": "A", "chr3": "B"}
    assert [r.getMessage() for r in caplog.records] == ["--assign: 2 line(s) with a single field skipped"]
    caplog.clear()
    f.write_text("chr1 0 100\nchr2 0 200\n")                         # a BED file by mistake: one label
    assert pipeline.read_assign(f, names) == {"chr1": "0", "chr2": "0"}
    assert [r.getMessage() for r in caplog.records] == ["--assign: one group only; nothing to date or paint"]


@pytest.mark.slow
def test_split_homoeologs_true_split_vs_cut_inside_one_subgenome(tmp_path):
    sc = {"root_age": 8e6, "chromosomes": {"n": 5, "length": 6e6}, "families": {"n": 80},
          "activity": {"rate": 0.002, "p_burst": 0.0}, "tree": {"split": 2.0e6, "children": ["A", "B"]},
          "mergers": [{"time": 1.0e5, "lineages": ["A", "B"]}], "exchanges": {"nonreciprocal": 0, "reciprocal": 0}}
    out, truth = builders.simulate(sc, tmp_path / "sim", tmp_path / "truth", seed=3)
    a = sg.parse_args(["--ltr_fasta", str(out / "ltr.fa"), "--annotation", str(out / "k2p.tsv"),
                        "--fai", str(out / "genome.fai"), "--outdir", str(tmp_path / "o"), "-t", "4", "--seed", "1"])
    run.setup(a)
    lib = inputs.load_library(a)
    cl = lineages.load_or_build_clades(a, lib)
    n_el = np.bincount(lib.seq, minlength=len(lib.seqnames))
    phased = [i for i in range(len(lib.seqnames)) if n_el[i] >= cst.MIN_ELEMENTS]
    pos = {s: k for k, s in enumerate(phased)}
    ek = np.array([pos.get(int(s), -1) for s in lib.seq], np.int64)
    ok = (ek >= 0) & np.isfinite(lib.k2p)
    edges = [float(x) for x in np.quantile(lib.k2p[ok], [1 / 3, 2 / 3])]
    tsg = builders.truth_subgenomes(truth)
    names = [lib.seqnames[i] for i in phased]
    true = np.array([0 if tsg[n] == sorted(set(tsg.values()))[0] else 1 for n in names])
    g = np.arange(len(names))
    assert pipeline.split_homoeologs(a, lib, cl, ek, edges, g, true) >= 1
    ga = g[true == 0]
    cut = (np.arange(ga.size) >= ga.size // 2).astype(int)                # one subgenome cut in two
    assert pipeline.split_homoeologs(a, lib, cl, ek, edges, ga, cut) == 0
