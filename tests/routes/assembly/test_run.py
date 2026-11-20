"""Running the assembly route over a sheet (routes.assembly.run): the elements of the annotated genomes, the ordered
pairs, the pair option, and what a run returns."""
from __future__ import annotations

import random
from collections import Counter
from contextlib import ExitStack
from dataclasses import replace
from pathlib import Path

import builders
import pytest
from markers import needs_lastz

from ltrtrace.base.log import LtrtraceError, UsageError
from ltrtrace.io.samples import Sample
from ltrtrace.routes.assembly import run


def _genome(tmp_path, label, seq, gff=None):
    fa = builders.write_fasta(tmp_path / f"{label}.fa", {"chr1": seq})
    return Sample(label, fa, None if gff is None else str(gff), None, None)


@pytest.fixture
def small(tmp_path):
    """A (annotated: the E and P toy cases), B (A's target) and C (unannotated, a copy of A)."""
    pair = builders.build_pair(random.Random(3), ["E", "P"])
    gff = builders.write_gff3(tmp_path / "A.gff3", pair.elements)
    return pair, [_genome(tmp_path, "A", pair.a, gff), _genome(tmp_path, "B", pair.b), _genome(tmp_path, "C", pair.a)]


def test_ordered_pairs_follow_the_sheet(small):
    _, samples = small
    assert run.ordered_pairs(samples) == [("A", "B"), ("A", "C")]
    b = replace(samples[1], annotation=samples[0].annotation)
    assert run.ordered_pairs([samples[0], b, samples[2]]) == [("A", "B"), ("A", "C"), ("B", "A"), ("B", "C")]


def test_load_elements_reads_the_annotated_genomes(small, tmp_path):
    """The elements stage (io.annotation.stage) for each annotated genome of the sheet, none for the others."""
    pair, samples = small
    with ExitStack() as stack:
        fastas = run.open_fastas(samples, tmp_path / "out", stack)
        elements, skipped = run.load_elements(samples, fastas, tmp_path / "out")
    assert set(elements) == {"A"} and skipped == {"A": Counter()}
    assert [(e.start, e.end) for e in elements["A"]] == [(e.start, e.end) for e in pair.elements]
    assert (tmp_path / "out" / "work" / "elements" / "A.tsv").is_file()


def test_pair_option_is_checked(small, tmp_path):
    _, samples = small
    for pair, msg in ((("A", "Z"), "'Z' is not a label"), (("B", "A"), "B has no annotation")):
        with pytest.raises(UsageError, match=msg):
            run.genotype(samples, tmp_path / "out", 1, pair=pair)
    with pytest.raises(UsageError, match="not a directory"):
        run.genotype(samples, tmp_path / "out", 1, paf=tmp_path / "none")


def test_missing_lastz_fails_first(small, tmp_path, monkeypatch):
    _, samples = small
    monkeypatch.setenv("PATH", str(tmp_path))
    with pytest.raises(LtrtraceError, match="needs lastz"):
        run.genotype(samples, tmp_path / "out", 1)
    assert not (tmp_path / "out").exists()


@needs_lastz
def test_genotype_one_pair_then_all(small, tmp_path, caplog):
    """--pair genotypes that pair only; a run over every pair then skips it."""
    _, samples = small
    out = tmp_path / "out"
    run.genotype(samples, out, 2, pair=("A", "C"))
    assert sorted(p.name for p in (out / "pairs").glob("*.gz")) == ["A__C.calls.tsv.gz"]
    calls = builders.read_calls(out / "pairs" / "A__C.calls.tsv.gz")
    assert {r["state"] for r in calls.values()} == {"P"}  # C is a copy of A
    with caplog.at_level("INFO", logger="ltrtrace"):
        genomes = run.genotype(samples, out, 2)
    assert sorted(p.name for p in (out / "pairs").glob("*.gz")) == ["A__B.calls.tsv.gz", "A__C.calls.tsv.gz"]
    assert "A vs C: calls done, skip" in caplog.text
    assert set(genomes.lengths) == {"A", "B", "C"} and genomes.lengths["B"] == {"chr1": len(samples_b_seq(samples))}
    assert set(genomes.elements) == {"A"} and genomes.skipped == {"A": Counter()}


def samples_b_seq(samples):
    """B's sequence, read back from its FASTA."""
    return "".join(line.strip() for line in Path(samples[1].fasta).read_text().splitlines()[1:])
