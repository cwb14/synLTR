"""k-mer lineages of a run: the clades of the run's library, a checkpoint stage in its output folder."""
from pathlib import Path

import numpy as np
import pytest

from ltrtrace.analyses.subgenomes import inputs, lineages, run
from ltrtrace.seq import kmers

TINY = Path(__file__).resolve().parents[2] / "golden" / "nucleus_tiny" / "lib"


def test_clades_are_a_checkpoint_of_the_run(tmp_path, monkeypatch):
    """The run's clades are written once under <outdir>/work/subgenomes with their done marker (no other cache
    folder); a rerun on the same library loads them, nothing rebuilt; the FracMinHash scale is recorded on them
    (exact for a small library)."""
    a = run.options([TINY / "ltr.fa"], [TINY / "annotation.gff3"], tmp_path / "o", fai=[TINY / "genome.fai"], threads=2)
    run.setup(a)
    lib = inputs.load_library(a)
    cl = lineages.load_or_build_clades(a, lib)
    work = tmp_path / "o" / "work" / "subgenomes"
    assert (work / "clades.npz").is_file() and (work / "clades.npz.done.json").is_file() and cl.scale == 1
    assert sorted(p.name for p in (tmp_path / "o").iterdir()) == ["work"]
    monkeypatch.setattr(kmers, "_build", lambda *x, **kw: pytest.fail("rebuilt"))
    again = lineages.load_or_build_clades(a, lib)
    assert all(np.array_equal(getattr(cl, f), getattr(again, f)) for f in ("members", "size", "mult", "nseq"))
    assert (again.n_occ, again.n_multi, again.scale) == (cl.n_occ, cl.n_multi, 1)
