"""The elements stage (io.annotation.stage): one genome's annotation loaded once, then read back from its table,
each element's k_kind restored from the done marker; and the checksums that key it."""
from __future__ import annotations

import json
import math
import random
from collections import Counter

import builders
import pytest

from synltr.base.checkpoint import code_files
from synltr.io.annotation.stage import elements_stage, input_checksums
from synltr.io.elements import read_elements
from synltr.io.samples import Sample
from synltr.io.seqio import Fasta


@pytest.fixture
def genome(tmp_path):
    """Genome A (the E and P toy cases) and its GFF3: (toy pair, sample)."""
    pair = builders.build_pair(random.Random(3), ["E", "P"])
    fa = builders.write_fasta(tmp_path / "A.fa", {"chr1": pair.a})
    gff = builders.write_gff3(tmp_path / "A.gff3", pair.elements)
    return pair, Sample("A", fa, str(gff), None, None)


def _stage(sample, work):
    with Fasta(sample.fasta) as fa:
        return elements_stage(sample.label, sample.annotation, fa, work, input_checksums(sample, work / "sums.json"))


def test_elements_stage_is_cached_and_keeps_k_kind(genome, tmp_path, caplog):
    """The elements are read back from work/elements/<label>.tsv on every run; k_kind, which the file does not
    hold, comes back from the stage's done marker. A second run reuses the file."""
    pair, sample = genome
    work = tmp_path / "work"
    elements, skipped = _stage(sample, work)
    tsv = work / "elements" / "A.tsv"
    assert skipped == Counter()
    assert [(e.start, e.end) for e in elements] == [(e.start, e.end) for e in pair.elements]
    assert all(e.k_kind == "K2P" for e in elements)
    assert [e.k_kind for e in read_elements(tsv)] == [None] * len(pair.elements)  # the file has no kind
    assert json.loads((tsv.parent / "A.tsv.done.json").read_text())["meta"] == {"k_kind": "K2P"}
    mtime = tsv.stat().st_mtime_ns
    with caplog.at_level("INFO", logger="synltr"):
        again, _ = _stage(sample, work)
    assert again == elements and tsv.stat().st_mtime_ns == mtime and "elements done, skip" in caplog.text


def test_elements_stage_restores_mixed_kinds(genome, tmp_path):
    """EDTA's ltr_identity gives a JC69 distance, labeled JC69 and never K2P; a genome whose records mix the two
    kinds keeps each element's own label."""
    _, sample = genome
    gff = tmp_path / "A.gff3"
    text = gff.read_text()
    assert "K2P_d=0.0\n" in text
    gff.write_text(text.replace("K2P_d=0.0\n", "ltr_identity=0.95\n", 1))
    for _ in range(2):  # computed, then read back
        elements, _ = _stage(sample, tmp_path / "work")
        assert [e.k_kind for e in elements] == ["JC69", "K2P"]
        assert elements[0].k == pytest.approx(-0.75 * math.log(1 - 4 * 0.05 / 3))


def test_elements_stage_follows_the_annotation(genome, tmp_path):
    pair, sample = genome
    first, _ = _stage(sample, tmp_path / "work")
    builders.write_gff3(tmp_path / "A.gff3", pair.elements[:1])
    second, _ = _stage(sample, tmp_path / "work")
    assert len(first) == 2 and len(second) == 1


def test_elements_stage_key_holds_its_own_code_and_the_loaders():
    files = code_files("io.annotation.stage")
    assert {"io/annotation/stage.py", "io/annotation/gff3.py", "io/annotation/ltrquest.py", "io/annotation/bed.py",
            "io/elements.py"} <= set(files)
    assert not any(f.startswith("routes/") for f in files)


def test_input_checksums_name_each_input(genome, tmp_path):
    _, sample = genome
    sums = tmp_path / "sums.json"
    got = input_checksums(sample, sums)
    assert set(got) == {"fasta", "annotation"} and all(len(v) == 64 for v in got.values())
    assert set(input_checksums(Sample("B", sample.fasta, None, None, None), sums)) == {"fasta"}
