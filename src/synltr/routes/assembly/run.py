"""The assembly route over a samples sheet: each annotated genome's elements, then each ordered pair's call file.

Stages, each checkpointed (a rerun skips what is done and resumes what is not):
  elements  each annotated genome's annotation, loaded and normalized: <outdir>/work/elements/<label>.tsv
            (+ <label>.skipped.json), the cached elements stage (io.annotation.stage)
  pairs     each ordered pair (annotated A, any other B; sheet order), or the one pair asked for:
            <outdir>/pairs/<A>__<B>.calls.tsv.gz (pair.genotype_pair, which anchors the pair first)
Every FASTA's .fai index is written next to the FASTA when its folder is writable, else to <outdir>/work/fai,
so a FASTA in a read-only folder works in place.
"""
from __future__ import annotations

import logging
import shutil
import time
from collections import Counter
from contextlib import ExitStack
from dataclasses import dataclass
from pathlib import Path

from synltr.base.constants import CHECKSUM_CACHE
from synltr.base.log import SynltrError, UsageError
from synltr.io.annotation.stage import elements_stage, input_checksums
from synltr.io.elements import Element
from synltr.io.samples import PAIR_SEP, Sample
from synltr.io.seqio import Fasta, open_fasta
from synltr.routes.assembly.anchor import FlankMm2Route, PafRoute, PairContext
from synltr.routes.assembly.pair import genotype_pair

log = logging.getLogger(__name__)


@dataclass
class Genomes:
    """What a run read about the genomes of its pairs: each one's sequence lengths (fai order), and each annotated
    one's elements and its loader's skip counts."""

    lengths: dict[str, dict[str, int]]
    elements: dict[str, list[Element]]
    skipped: dict[str, Counter]


def genotype(samples: list[Sample], outdir: Path, threads: int, pair: tuple[str, str] | None = None,
             paf: Path | None = None) -> Genomes:
    """Write the call file of every ordered pair of `samples` (annotated A, any other B; sheet order), or of
    `pair` only, under `outdir`, with `threads` in all; anchored by flank pieces (FlankMm2Route), or by the
    whole-genome alignments in `paf` (PafRoute). Returns what it read of the genomes in those pairs.

    Fails before any work (`check_run`) on a `pair` label not in the sheet or whose first genome has no annotation
    (UsageError), a `paf` that is not a directory (UsageError) and a missing lastz (SynltrError).
    """
    t0 = time.perf_counter()
    by = {s.label: s for s in samples}
    check_run(samples, pair=pair, paf=paf)
    route = FlankMm2Route() if paf is None else PafRoute(paf)
    pairs = [pair] if pair is not None else ordered_pairs(samples)
    used = [s for s in samples if any(s.label in p for p in pairs)]
    log.info("genotype: %d genome(s), %d pair(s), route %s, %d thread(s) -> %s", len(used), len(pairs), route.name,
             threads, outdir)
    with ExitStack() as stack:
        fastas = open_fastas(used, outdir, stack)
        elements, skipped = load_elements(used, fastas, outdir)
        for k, (a, b) in enumerate(pairs, 1):
            log.info("pair %d/%d: %s vs %s", k, len(pairs), a, b)
            ctx = PairContext(a=by[a], b=by[b], a_fa=fastas[a], b_fa=fastas[b], a_elements=elements[a],
                              b_elements=elements.get(b), outdir=Path(outdir), threads=threads)
            genotype_pair(ctx, route)
        lengths = {x: fa.lengths() for x, fa in fastas.items()}
    log.info("genotype: %d pair(s) done in %.1f s", len(pairs), time.perf_counter() - t0)
    return Genomes(lengths, elements, skipped)


def check_run(samples: list[Sample], *, pair: tuple[str, str] | None = None, paf: Path | None = None) -> None:
    """The route's checks, made before any work: a `pair` label not in the sheet or whose first genome has no
    annotation (UsageError), a `paf` that is not a directory (UsageError) and a missing lastz (SynltrError)."""
    if pair is not None:
        _check_pair(pair, {s.label: s for s in samples})
    if paf is not None and not Path(paf).is_dir():
        raise UsageError(f"--paf {paf} is not a directory")
    if shutil.which("lastz") is None:
        raise SynltrError("genotyping assemblies needs lastz on PATH (not found); install it with the conda "
                            "environment in environment.yml (bioconda's lastz)")


def ordered_pairs(samples: list[Sample]) -> list[tuple[str, str]]:
    """Every ordered pair (annotated A, any other B), in sheet order."""
    return [(a.label, b.label) for a in samples if a.annotation is not None for b in samples if b is not a]


def pair_path(outdir: Path, a: str, b: str) -> Path:
    """The call file of pair (a, b): <outdir>/pairs/<a>__<b>.calls.tsv.gz."""
    return Path(outdir) / "pairs" / f"{a}{PAIR_SEP}{b}.calls.tsv.gz"


def _check_pair(pair: tuple[str, str], by: dict[str, Sample]) -> None:
    """Fails (UsageError) on a pair label not in the sheet or a first genome without annotation."""
    labels = ", ".join(by)
    for x in pair:
        if x not in by:
            raise UsageError(f"--pair {','.join(pair)}: {x!r} is not a label of the samples sheet ({labels})")
    if by[pair[0]].annotation is None:
        raise UsageError(f"--pair {','.join(pair)}: {pair[0]} has no annotation; the first genome of a pair must "
                         "be annotated (its elements are called)")


def open_fastas(samples: list[Sample], outdir: Path, stack: ExitStack) -> dict[str, Fasta]:
    """Each genome's FASTA, its index cached in <outdir>/work/fai when its own directory is read-only."""
    cache = Path(outdir) / "work" / "fai"
    out = {}
    for s in samples:
        out[s.label] = stack.enter_context(open_fasta(s.fasta, s.label, cache))
    return out


def load_elements(samples: list[Sample], fastas: dict[str, Fasta],
                  outdir: Path) -> tuple[dict[str, list[Element]], dict[str, Counter]]:
    """The elements and skip counts of each annotated genome among `samples` (the elements stage)."""
    work = Path(outdir) / "work"
    elements, skipped = {}, {}
    for s in samples:
        if s.annotation is not None:
            inputs = input_checksums(s, work / CHECKSUM_CACHE)
            elements[s.label], skipped[s.label] = elements_stage(s.label, s.annotation, fastas[s.label], work,
                                                                 inputs)
    return elements, skipped
