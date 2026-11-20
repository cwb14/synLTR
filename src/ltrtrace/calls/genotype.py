"""The genotype command: each annotated element's presence or absence in every other genome (assemblies), or its
genotype in each read set (reads), then the outputs.

On a sheet of assemblies it runs the assembly route (routes.assembly.run.genotype) over every ordered pair
(annotated A, any other B; sheet order), or with --pair A,B over that pair only, to spread a run over jobs, and
stops there. Without --pair it then builds the loci over every pair's call file (calls.loci.build_loci) and writes,
in the outdir: loci.tsv, pav_matrix.tsv and qc.tsv (calls.tables), pav.vcf.gz with its .tbi and pav.unplaced.tsv
on the first sheet row's coordinates (calls.vcf), and <label>.pav.gff3 and <label>.pav.bed per genome
(calls.tracks). Finished pairs are skipped on a rerun, so a run killed between pairs resumes where it stopped. The
run's record (calls.manifest.Finished, written as manifest.json by the command line) names the route assembly, the
checksums of the FASTAs and annotations the route's stages keyed on (base.checkpoint.keyed_checksums), with --paf
each PAF file's, and each annotated genome's skip counts. A run removes the folder's earlier manifest.json once the
route's checks have passed (routes.assembly.run.check_run: the --pair labels, the --paf folder, lastz); a --pair job
finishes no run, returns its exit code and leaves the folder's manifest.json as it was.

On a sheet of one reference (role reference, with a FASTA and an annotation) and read sets (the reads column: a
BAM aligned to that reference) it runs the reads route (routes.reads.run.genotype, --ploidy and --min-mapq), which
writes ltrkin's reads.genotypes.tsv, then the run's reads.vcf.gz with its .tbi and reads.unplaced.tsv (elements the
VCF cannot anchor: those starting at base 1; calls.vcf.write_reference_vcf), every read set a sample, in sheet order.
Its record (route reads, --ploidy and --min-mapq as used, the checksum of the reference FASTA, which every read set's
evidence is keyed on) becomes the run's manifest.json, from which `source` finds the reference and its annotation; the
folder's earlier manifest.json is removed once the options are checked. Finished read sets' evidence is reused on a
rerun.
"""
from __future__ import annotations

import argparse
import logging
import time
from contextlib import ExitStack
from pathlib import Path

from ltrtrace.base.checkpoint import keyed_checksums
from ltrtrace.base.constants import DEFAULT_MIN_MAPQ, DEFAULT_PLOIDY, READS_STEM, RESOLVED_WARN
from ltrtrace.base.log import UsageError
from ltrtrace.calls.loci import build_loci
from ltrtrace.calls.manifest import Finished, paf_checksums, remove_manifest, sample_checksums
from ltrtrace.calls.tables import write_loci_tsv, write_pav_matrix, write_qc
from ltrtrace.calls.tracks import write_tracks
from ltrtrace.calls.vcf import write_reference_vcf, write_vcf
from ltrtrace.io.samples import Sample, SampleSheetError, read_samples
from ltrtrace.io.seqio import Fasta
from ltrtrace.routes.assembly import run as assembly
from ltrtrace.routes.reads import run as reads

log = logging.getLogger(__name__)


def run(args: argparse.Namespace) -> int | Finished:
    """`ltrtrace genotype SHEET`: args.sheet, args.outdir, args.threads, args.pair (A, B) or None, args.paf (a
    directory) or None (assemblies), args.ploidy and args.min_mapq or None (read sets). Returns the run's record, or
    the exit code of a --pair job."""
    t0 = time.perf_counter()
    samples = read_sheet(args.sheet)
    outdir = Path(args.outdir)
    if any(s.reads is not None for s in samples):
        if args.pair is not None or args.paf is not None:
            raise UsageError("--pair and --paf genotype assemblies; this sheet holds read sets")
        ploidy = DEFAULT_PLOIDY if args.ploidy is None else args.ploidy
        min_mapq = DEFAULT_MIN_MAPQ if args.min_mapq is None else args.min_mapq
        remove_manifest(outdir)
        with keyed_checksums() as keyed:
            n = genotype_reads(samples, outdir, args.threads, ploidy=ploidy, min_mapq=min_mapq)
        log.info("done in %.1f s: %d read set(s) -> %s", time.perf_counter() - t0, n, outdir)
        return Finished(outdir, tuple(samples), route="reads", checksums=sample_checksums(samples, keyed, ("fasta",)),
                        params={"ploidy": ploidy, "min_mapq": min_mapq})
    if args.ploidy is not None or args.min_mapq is not None:
        raise UsageError("--ploidy and --min-mapq apply to read sets (the reads column); this sheet holds assemblies")
    done = genotype_samples(samples, outdir, args.threads, pair=args.pair, paf=args.paf)
    if done is None:
        log.info("done in %.1f s: pair %s; run without --pair to build the loci", time.perf_counter() - t0,
                 ",".join(args.pair))
        return 0
    log.info("done in %.1f s -> %s", time.perf_counter() - t0, outdir)
    return done


def read_sheet(path: Path) -> list[Sample]:
    """The samples sheet at `path`, checked for genotyping: a UsageError for a sheet that cannot be read, a sheet
    of assemblies with a row without a FASTA, or a sheet of read sets with a row (other than the reference) without
    reads."""
    try:
        samples = read_samples(path)
    except SampleSheetError as exc:
        raise UsageError(str(exc)) from exc
    if any(s.reads is not None for s in samples):
        idle = [s.label for s in samples if s.reads is None and s.role != "reference"]
        if idle:
            raise UsageError(f"{path}: a sheet of read sets holds one reference and the read sets aligned to it; "
                             f"row(s) {', '.join(idle)} have no reads")
        return samples
    no_fasta = [s.label for s in samples if s.fasta is None]
    if no_fasta:
        raise UsageError(f"{path}: genotype aligns assemblies, so every row needs a fasta; none for "
                         f"{', '.join(no_fasta)}")
    return samples


def genotype_samples(samples: list[Sample], outdir: Path, threads: int, *, pair: tuple[str, str] | None = None,
                     paf: Path | None = None) -> Finished | None:
    """The assembly route over every ordered pair of `samples` (or `pair` only; routes.assembly.run.genotype),
    finished pairs skipped, then without `pair` the loci and every output (`write_outputs`). Without `pair`, the
    folder's earlier manifest.json is removed once the route's checks have passed (routes.assembly.run.check_run).
    Returns the run's record (route assembly, the FASTA and annotation checksums the route's stages keyed on, with
    `paf` each PAF file's beside it, each annotated genome's skip counts); None with `pair`, which builds no loci."""
    outdir = Path(outdir)
    assembly.check_run(samples, pair=pair, paf=paf)
    if pair is None:
        remove_manifest(outdir)
    with keyed_checksums() as keyed:
        genomes = assembly.genotype(samples, outdir, threads, pair=pair, paf=paf)
    if pair is not None:
        return None
    write_outputs(samples, outdir, genomes)
    return assembly_record(samples, outdir, keyed, skipped=genomes.skipped, paf=paf)


def assembly_record(samples: list[Sample], outdir: Path, keyed: dict[str, str], *,
                    skipped: dict[str, dict[str, int]], paf: Path | None = None) -> Finished:
    """The record of a finished run on assemblies in `outdir`: route assembly, the checksums of the FASTAs and
    annotations its stages keyed on (`keyed`, a base.checkpoint.keyed_checksums collection), with `paf` each PAF
    file's beside it, and each annotated genome's skip counts (`skipped`)."""
    return Finished(Path(outdir), tuple(samples), route="assembly",
                    checksums=sample_checksums(samples, keyed, ("fasta", "annotation")), skipped=skipped,
                    params={} if paf is None else {"paf_checksums": paf_checksums(paf, keyed)})


def write_outputs(samples: list[Sample], outdir: Path, genomes: assembly.Genomes) -> int:
    """build_loci over the call file of every ordered pair, and every output, from what the route read of the
    genomes (`genomes`: every genome of the sheet); returns the number of loci."""
    outdir = Path(outdir)
    calls = {(a, b): assembly.pair_path(outdir, a, b) for a, b in assembly.ordered_pairs(samples)}
    order = {x: {sid: i for i, sid in enumerate(lengths)} for x, lengths in genomes.lengths.items()}
    loci, qc = build_loci(samples, genomes.elements, calls, order, skipped=genomes.skipped)
    for scope, metric, value in qc:
        if metric == "resolved_frac" and value < RESOLVED_WARN:  # nan (no calls) never warns
            a, b = scope.split(">")
            log.warning("%s vs %s: %.2f%% of the calls resolved (P/S/T/E): no orthologous alignment between "
                        "%s and %s?", a, b, 100 * value, a, b)
    write_loci_tsv(loci, samples, outdir / "loci.tsv")
    write_pav_matrix(loci, samples, outdir / "pav_matrix.tsv")
    write_qc(qc, outdir / "qc.tsv")
    ref = samples[0]
    with ExitStack() as stack:  # the VCF's coordinates and REF bases
        fasta = assembly.open_fastas([ref], outdir, stack)[ref.label]
        _, unplaced = write_vcf(loci, samples, ref.label, fasta, outdir / "pav.vcf.gz")
    write_tracks(loci, samples, outdir, order, genomes.lengths)
    with open(unplaced) as fh:
        n_unplaced = sum(1 for _ in fh) - 1
    log.info("%d loci -> loci.tsv, pav_matrix.tsv, qc.tsv, pav.vcf.gz (%d on %s's coordinates, %d unplaced), "
             "<label>.pav.gff3/.bed", len(loci), len(loci) - n_unplaced, ref.label, n_unplaced)
    return len(loci)


def genotype_reads(samples: list[Sample], outdir: Path, threads: int, *, ploidy: int, min_mapq: int) -> int:
    """The reads route over every read set of `samples` (routes.reads.run.genotype: reads.genotypes.tsv), then the
    run's reads.vcf.gz on the reference's coordinates, every read set a sample. Returns the number of read sets."""
    outdir = Path(outdir)
    done = reads.genotype(samples, outdir, threads, ploidy=ploidy, min_mapq=min_mapq)
    with Fasta(done.reference, cache_dir=outdir / "work" / "fai") as fasta:
        _, unplaced = write_reference_vcf(outdir / f"{READS_STEM}.vcf.gz", done.targets, done.samples, fasta,
                                          done.reference)
    with open(unplaced) as fh:
        n_unplaced = sum(1 for _ in fh) - 1
    log.info("%d read set(s) -> %s.genotypes.tsv, %s.vcf.gz (%d element(s) unplaced: %s)", len(done.samples),
             READS_STEM, READS_STEM, n_unplaced, unplaced.name)
    return len(done.samples)
