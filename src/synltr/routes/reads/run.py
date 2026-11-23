"""The reads route over a samples sheet: one reference (its FASTA and annotation) and read sets aligned to it.

`genotype` loads the reference row's elements with ltrkin's reader rules (`load_targets`), checks every read set's
BAM against the reference before any work (`check_bams`), then genotypes the read sets in sheet order
(genotype.genotype_sample, each with its evidence checkpoint <outdir>/reads.<S>.evidence.json.gz, S the row's label)
and writes ltrkin's one table for the run, <outdir>/reads.genotypes.tsv, every read set a sample. A rerun reuses each
read set's finished evidence, so a run killed between read sets resumes where it stopped. It returns the calls;
calls.genotype writes the run's VCF, <outdir>/reads.vcf.gz, from them.

Coordinates: the io Elements (1-based closed) become route Targets (0-based half-open) here, in `load_targets`.
"""
from __future__ import annotations

import logging
import time
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import pysam

from synltr.base.constants import READS_STEM, READS_TSD_ORDER
from synltr.base.log import UsageError
from synltr.io.annotation import RecordSkips, load_annotation
from synltr.io.annotation.nesting import containment
from synltr.io.elements import Element
from synltr.io.samples import Sample
from synltr.io.seqio import Fasta, open_fasta
from synltr.routes.reads.genotype import SampleCalls, genotype_sample, write_genotypes
from synltr.routes.reads.reads import Target

log = logging.getLogger(__name__)


@dataclass(frozen=True)
class ReadsRun:
    """What a reads run genotyped: the reference FASTA, its elements (route order) and each read set's calls (sheet
    order)."""

    reference: Path
    targets: list[Target]
    samples: list[SampleCalls]


def genotype(samples: list[Sample], outdir: Path, threads: int, *, ploidy: int, min_mapq: int) -> ReadsRun:
    """Genotype every read set of `samples` (rows with reads) at the elements of the reference row, with `threads`
    in all; ltrkin's table <outdir>/reads.genotypes.tsv (read sets in sheet order) and each read set's evidence
    checkpoint. Fails before any work (UsageError) on a reference or annotation that cannot be read or does not
    match, and on a BAM without an index or aligned to another reference."""
    t0 = time.perf_counter()
    ref = next(s for s in samples if s.role == "reference")
    sets = [s for s in samples if s.reads is not None]
    outdir = Path(outdir)
    with open_fasta(ref.fasta, ref.label, outdir / "work" / "fai") as fasta:
        targets, skipped = load_targets(ref.annotation, fasta, ref.label)
        check_bams(sets, ref.fasta, fasta, sorted({t.chrom for t in targets}))
    log.info("genotype: %d read set(s) on %s (%d elements; %d skipped: fewer than two LTRs), ploidy %d, %d "
             "thread(s) -> %s", len(sets), ref.label, len(targets), skipped["no_ltr_pair"], ploidy, threads, outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    done = [genotype_sample(s.label, s.reads, ref.fasta, targets, threads, ploidy=ploidy, min_mapq=min_mapq,
                            outdir=outdir) for s in sets]
    write_genotypes(outdir / f"{READS_STEM}.genotypes.tsv", [lines for _, lines in done])
    log.info("genotype: %d read set(s) done in %.1f s", len(sets), time.perf_counter() - t0)
    return ReadsRun(ref.fasta, targets, [calls for calls, _ in done])


def load_targets(annotation: str, fasta: Fasta, genome: str) -> tuple[list[Target], Counter]:
    """The elements of `annotation` on `fasta` read with ltrkin's rules (the outermost two of two or more LTRs
    bound an element, the clade is the last field of a classification, a TSD is the annotated one if it is 3-8 bp
    of ACGT, else the longest repeated flank of 6, 5 or 4 bp; the family is the family attribute, else a Name that
    is not coordinates; the element's ID names it; elements on sequences the FASTA lacks are skipped with ltrkin's
    warning), with ltrkin's nesting counts, as Targets in ltrkin's order (sequence name as text, start, longest
    first). Returns them and the loader's skip counts."""
    skipped = RecordSkips()
    table = load_annotation(annotation, genome, fasta, skipped=skipped, ltr_rule="outermost", class_rule="last",
                            tsd_rule="clean", tsd_order=READS_TSD_ORDER, name_family=True, name_rule="id",
                            missing_rule="skip")
    n_inner, n_outer, inner = containment(table)
    rows = list(table)
    order = sorted(range(len(rows)), key=lambda i: (rows[i].seqid, rows[i].start, -rows[i].end))
    return [_target(rows[i], n_inner[i], n_outer[i], inner[i]) for i in order], skipped.counts


def _target(e: Element, n_inner: int, n_outer: int, inner: list[tuple[int, int]]) -> Target:
    """The route's element (0-based half-open) of io Element `e` (1-based closed) and its nesting counts."""
    return Target(e.name, e.seqid, e.zstart, e.end, e.strand, (e.ltr_left[0] - 1, e.ltr_left[1]),
                  (e.ltr_right[0] - 1, e.ltr_right[1]), e.family or ".", e.clade or ".", e.tsd or "", n_inner,
                  n_outer, tuple((s - 1, t) for s, t in inner))


def check_bams(sets: list[Sample], reference: Path, fasta: Fasta, chroms: list[str]) -> None:
    """Fail before any work (UsageError): every read set's BAM readable, indexed and aligned to this reference
    (the same length for each element sequence it holds). Warn about a BAM lacking element sequences (their
    elements are not genotyped in that read set) or holding several samples' read groups."""
    for s in sets:
        try:
            bam = pysam.AlignmentFile(str(s.reads), reference_filename=str(reference))
        except (OSError, ValueError) as exc:
            raise UsageError(f"{s.label}: cannot read {s.reads}: {exc}") from None
        with bam:
            if not bam.has_index():
                raise UsageError(f"{s.label}: {s.reads} has no index (samtools index it)")
            refs = set(bam.references)
            held = [c for c in chroms if c in refs]
            if not held:
                raise UsageError(f"{s.label}: no element sequence is in {s.reads} (sequence names differ from the "
                                 "reference?)")
            diff = [c for c in held if bam.get_reference_length(c) != fasta.length(c)]
            if diff:
                c = diff[0]
                raise UsageError(f"{s.label}: {s.reads} was aligned to another reference: {len(diff)} sequence(s) "
                                 f"differ in length, e.g. {c} {bam.get_reference_length(c)} bp in the BAM, "
                                 f"{fasta.length(c)} in the FASTA")
            if len(held) < len(chroms):
                log.warning("%s: %s lacks %d element sequence(s); their elements are not genotyped in this read set",
                            s.label, s.reads, len(chroms) - len(held))
            sm = {rg.get("SM") for rg in bam.header.to_dict().get("RG", []) if rg.get("SM")}
            if len(sm) > 1:
                log.warning("%s: %s has read groups of %d samples %s: genotyped as one sample (give one BAM per "
                            "sample)", s.label, s.reads, len(sm), sorted(sm)[:3])
