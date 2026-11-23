"""The geneflow command: every pair of a sheet of assemblies genotyped first (the genotype driver,
calls.genotype; finished pairs are skipped by their checkpoints; --paf anchors them with the user's whole-genome
alignments, as genotype --paf does), then the gene-flow analyses of the finished run (`analyses`). With --screen, or
on a sheet without genome FASTAs (said in the log and in findings.txt), the library screen instead (`screen`).

Routes (`route`, read before genotype.read_sheet, which refuses a sheet without FASTAs): --screen runs the screen on
the sheet's libraries (every row needs one; FASTAs are not read); a sheet without genome FASTAs runs it too; a sheet
with FASTAs on some rows only, and --screen with --identity or --paf, are usage errors. A sheet of read sets is a
usage error without --screen (`synltr genotype` genotypes read sets), after genotype.read_sheet's own checks.
Before any of this, the modules the analyses import where they use them (ANALYSIS_IMPORTS) are imported, so a broken
environment fails before any genotyping.

A finished run returns its record (calls.manifest.Finished: route assembly or library, the checksums of the inputs
its stages keyed on in both of its phases, the analyses' included: the identity layer keys the FASTAs again, the
screen's trio stage the libraries), which the command line writes as manifest.json; `synltr plot` reads the samples
sheet and the route from it. A run removes the folder's earlier manifest.json once its checks have passed (the
genotype driver; the screen before its pairs), and the analyses remove the earlier overview figure with stats.tsv and
findings.txt, so a run that fails leaves no record of other results.

The analyses, on the run's catalog (catalog.load_catalog, read once), under the assembly route's properties
(route.ROUTE_PROPS), in this order:
  the run summary (relations.run_summary) and the jackknife block length (row jackknife block_bp, the route's
    rule: route.block_bp), which every jackknife here takes
  relationships (relations.relations): distances, the tree, shared derived insertions, sister pairs, KKSC
  the trios (trio.trios: --trio, else roles from the data), then each trio's tests (trio.trio_tests)
  inverted blocks (inversions.tsv) and their rows
  mechanisms (mechanisms.tsv)
  with --identity, the sequence-identity layer of each trio (identity.tsv) and the catalog's pair-divergence layer
    (pairs/<D>__<Y>.div.tsv.gz); without it an earlier run's identity.tsv and assembly divergence files are
    deleted, with one WARNING, and one Dating line says what needs the layer
  each trio's painting (tracts.<R>.bed, posterior.<R>.tsv, decoy_tracts.<S>.bed, decoy_posterior.<S>.tsv), with
    the identity layer's self rows when it ran; the painting files written are recorded in
    work/geneflow/paintings.json, so the next run removes those of a genome no longer in the sheet
then findings.txt and stats.tsv (`write_results`), each written atomically, last: they are removed when the
analyses start, with the overview figure's files (overview.png, overview.pdf, overview_numbers.tsv, which the command
line draws once the run is recorded), so a failed run leaves none beside results they would not describe.

The library screen (`screen`): every pair's evidence from the libraries (routes.library.run.screen_pairs, finished
pairs skipped), then its analyses (`screen_analyses`) under the library route's properties (route.ROUTE_PROPS), in
this order:
  the jackknife block length (row jackknife block_bp, layer library: library_stats.block_bp)
  the trios: --trio, else every three ingroup genomes of the sheet (sheet order), roles by the single-copy k-mers
    (library_stats.roles, role_rows); each trio's library statistics with the gate (library_stats.compute, its
    alignments checkpointed per trio under work/geneflow/triples: library_stats.triples) and their findings
  inverted blocks from the screen's chain blocks (inversions.chain_inversions: inversions.tsv) and their rows
  each trio's paintings, one per alphabet of the route (paint.paint on the screen's records; files as on the
    assembly route, recorded in work/geneflow/paintings.json)
then findings.txt and stats.tsv (`write_results`). Relationships, D_pav and its strata, the assembly route's
f4-ratio, the identity layer and mechanisms need observed empty sites (E) from assemblies: they are not run, and
findings.txt says so (SCREEN_SKIPS).

The analyses rerun every time (they take seconds); only the identity alignment and the screen's trio alignments are
checkpointed (identity.py, library_stats.triples).
Exit status: 2 a usage error found before any genotyping (the samples sheet, --trio, the route's options); 1 a
genotyping or screen error; 3 an analysis error after the genotyping (or the screen's pairs) finished (one line
"error: analyses failed: ..."; the genotyping or screen outputs are kept).
"""
from __future__ import annotations

import argparse
import importlib
import itertools
import json
import logging
import os
import shutil
import time
from pathlib import Path

from synltr.analyses.geneflow.catalog import AnalysisError, Catalog, load_catalog
from synltr.analyses.geneflow.identity import SECTION as IDENTITY_SECTION
from synltr.analyses.geneflow.identity import assembly_divergences, write_identity, write_pair_divergences
from synltr.analyses.geneflow.identity import identity as identity_layer
from synltr.analyses.geneflow.inversions import (
    InvBlock,
    chain_inversion_rows,
    chain_inversions,
    inversion_rows,
    inverted_blocks,
    write_inversions,
)
from synltr.analyses.geneflow.library_stats import LibraryRun, compute, findings, role_rows, roles
from synltr.analyses.geneflow.library_stats import block_bp as library_block_bp
from synltr.analyses.geneflow.mechanisms import mechanisms, write_mechanisms
from synltr.analyses.geneflow.paint import Painting, library_rows, paint, paint_rows, write_paintings
from synltr.analyses.geneflow.relations import relations, run_summary
from synltr.analyses.geneflow.route import ROUTE_PROPS, block_bp, require
from synltr.analyses.geneflow.trio import Trio, trio_tests, trios
from synltr.base import constants
from synltr.base.checkpoint import atomic_path, atomic_write, keyed_checksums
from synltr.base.log import AnalysisError as AnalysesFailed
from synltr.base.log import UsageError
from synltr.calls import genotype
from synltr.calls.manifest import Finished, remove_manifest, sample_checksums
from synltr.io.samples import Sample, SampleSheetError, read_samples
from synltr.routes.library.run import check_screen, screen_pairs
from synltr.stats.report import SECTIONS, Findings, StatRow, write_findings, write_stats

log = logging.getLogger(__name__)

PROPS = ROUTE_PROPS["assembly"]  # the route whose pair files geneflow reads
SCREEN_PROPS = ROUTE_PROPS["library"]  # the screen's
PAINTINGS = Path("work") / "geneflow" / "paintings.json"  # under the run directory: the painting files written
_IDENTITY_LINE = ("Counts of shared insertions younger than the split need the sequence-identity "
                  "layer: rerun synltr geneflow with --identity")
SCREEN_SKIPS = ("Relationships (distances, the tree, sister pairs by shared insertions), D_pav with its strata and age "
                "gradient, the assembly route's f4-ratio, the sequence-identity layer and mechanisms are not run on the "
                "library screen: they need observed empty sites (E) from genome assemblies, which inferred absences "
                "(A) never stand for; run `synltr geneflow` without --screen on the genome assemblies for them")
SCREEN_ROUTED = ("The samples sheet holds no genome FASTA, so geneflow ran the library screen on its LTR-RT libraries "
                 "(as with --screen)")


def main(args: argparse.Namespace) -> Finished:
    """`synltr geneflow SHEET`: args.sheet, args.outdir, args.threads, args.paf (a directory) or None, args.trio
    (R, S, D) or None, args.identity and args.screen. Returns the run's record."""
    return geneflow(args.sheet, args.outdir, args.threads, trio=args.trio, identity=args.identity, paf=args.paf,
                    screen=args.screen)


def geneflow(samples: Path, outdir: Path, threads: int, *, trio: tuple[str, str, str] | None = None,
             identity: bool = False, paf: Path | None = None, screen: bool = False) -> Finished:
    """Genotype every pair of the samples sheet at `samples` into `outdir` with `threads` in all (the genotype
    driver: finished pairs are skipped), each anchored by the whole-genome alignments in `paf` when given (as
    `synltr genotype --paf`, whose pairs are then reused), then run the analyses (`analyses`) with `trio` (R, S,
    D; None: roles from the data) and the identity layer when `identity`. With `screen`, or on a sheet without
    genome FASTAs, the library screen instead (`run_screen`; module docstring, Routes). The analyses' imports, the
    sheet, the route's options and the trio are checked before any genotyping (an ImportError, a UsageError); an
    analysis failure is an AnalysisError (exit status 3), logged with the outputs kept. Returns the run's record."""
    t0 = time.perf_counter()
    for name in constants.ANALYSIS_IMPORTS:
        importlib.import_module(name)
    which = route(Path(samples), screen=screen, identity=identity, paf=paf)
    if which != "assembly":
        return run_screen(samples, outdir, threads, trio=trio, routed=which == "routed")
    sheet = genotype.read_sheet(samples)
    if any(s.reads is not None for s in sheet):
        raise UsageError(f"{samples} holds read sets (column reads): geneflow genotypes genome assemblies, or with "
                         "--screen compares LTR-RT libraries; synltr genotype genotypes read sets")
    check_trio_option(sheet, trio, Path(samples))
    with keyed_checksums() as keyed:  # both phases: the identity layer keys the FASTAs again
        done = genotype.genotype_samples(sheet, Path(outdir), threads, paf=paf)
        try:
            analyses(samples, outdir, threads, trio=trio, identity=identity)
        except Exception as exc:  # the genotyping outputs are complete: exit 3, with one line naming the failure
            what = " ".join(f"{type(exc).__name__}: {exc}".split())
            log.error("geneflow: analyses failed (%s); the genotyping outputs in %s are complete: fix the cause and "
                      "rerun synltr geneflow (finished pairs are skipped)", what, outdir)
            raise AnalysesFailed(f"analyses failed: {what}") from exc
    log.info("geneflow: done in %.1f s -> %s", time.perf_counter() - t0, outdir)
    return genotype.assembly_record(sheet, Path(outdir), keyed, skipped=done.skipped, paf=paf)


def route(sheet: Path, *, screen: bool, identity: bool, paf: Path | None) -> str:
    """The route of a geneflow run (module docstring, Routes): "screen" (--screen), "routed" (a sheet without genome
    FASTAs: the screen) or "assembly". UsageErrors: --screen, or a sheet without FASTAs, with --identity or --paf
    (both realign assemblies); a sheet with FASTAs on some rows only; an unreadable sheet."""
    try:
        samples = read_samples(sheet)
    except SampleSheetError as exc:
        raise UsageError(str(exc)) from exc
    no_fasta = [s.label for s in samples if s.fasta is None]
    if not screen and (any(s.reads is not None for s in samples) or not no_fasta):
        return "assembly"
    if not screen and len(no_fasta) < len(samples):
        raise UsageError(f"{sheet}: row(s) {', '.join(no_fasta)} have no genome FASTA while the others have one: give "
                         "every row a fasta to genotype the assemblies, or run synltr geneflow --screen to compare "
                         "the sheet's LTR-RT libraries (column library)")
    why = "--screen compares LTR-RT libraries" if screen else f"{sheet} holds no genome FASTA (the library screen)"
    for flag, given in (("--identity", identity), ("--paf", paf is not None)):
        if given:
            raise UsageError(f"{why}: {flag} realigns genome assemblies; run synltr geneflow without --screen on a "
                             "sheet with a fasta on every row for it")
    return "screen" if screen else "routed"


def run_screen(samples: Path, outdir: Path, threads: int, *, trio: tuple[str, str, str] | None = None,
               routed: bool = False) -> Finished:
    """The library screen of the samples sheet at `samples` into `outdir` with `threads` in all: every pair's
    evidence from the sheet's libraries (routes.library.run.screen_pairs; finished pairs are skipped), then its
    analyses (`screen_analyses`) with `trio` (R, S, D; None: roles from the data). routed: the sheet holds no genome
    FASTA (said in the log and in findings.txt). The sheet, the trio and every row's library are checked first (a
    UsageError), then the folder's earlier manifest.json is removed; an analysis failure is an AnalysisError (exit
    status 3), logged with the screen's outputs kept. Returns the run's record (route library, the checksums of the
    library files the stages of the pairs and of the analyses keyed on)."""
    t0 = time.perf_counter()
    try:
        sheet = read_samples(samples)
    except SampleSheetError as exc:
        raise UsageError(str(exc)) from exc
    check_trio_option(sheet, trio, Path(samples))
    check_screen(sheet)
    if routed:
        log.info("geneflow: %s holds no genome FASTA, so geneflow runs the library screen on its LTR-RT libraries (as "
                 "with --screen)", samples)
    remove_manifest(Path(outdir))
    records: dict = {}
    with keyed_checksums() as keyed:  # both phases: the trio stage keys the libraries again
        screen_pairs(sheet, Path(outdir), threads, records=records)
        try:
            screen_analyses(samples, outdir, threads, trio=trio, records=records, routed=routed)
        except Exception as exc:  # the screen's outputs are complete: exit 3, with one line naming the failure
            what = " ".join(f"{type(exc).__name__}: {exc}".split())
            log.error("geneflow: analyses failed (%s); the screen's outputs in %s are complete: fix the cause and "
                      "rerun synltr geneflow --screen (finished pairs and trios are skipped)", what, outdir)
            raise AnalysesFailed(f"analyses failed: {what}") from exc
    log.info("geneflow: done in %.1f s -> %s", time.perf_counter() - t0, outdir)
    return Finished(Path(outdir), tuple(sheet), route="library",
                    checksums=sample_checksums(sheet, keyed, ("library",)))


def check_trio_option(samples: list[Sample], trio: tuple[str, str, str] | None, sheet: Path) -> None:
    """UsageError unless `trio` (--trio) names three distinct labels of the sheet, none an outgroup; None: nothing
    to check."""
    if trio is None:
        return
    labels = [s.label for s in samples]
    spec = ",".join(trio)
    unknown = [x for x in trio if x not in labels]
    if unknown:
        raise UsageError(f"--trio {spec}: {unknown[0]!r} is not a label of {sheet} ({', '.join(labels)})")
    if len(set(trio)) != 3:
        raise UsageError(f"--trio {spec}: the recipient, sister and donor must be three different genomes")
    out = [s.label for s in samples if s.label in trio and s.role == "outgroup"]
    if out:
        raise UsageError(f"--trio {spec}: {', '.join(out)} has role outgroup in {sheet}; a trio is three ingroup "
                         "genomes")


def analyses(samples: Path, outdir: Path, threads: int, *, trio: tuple[str, str, str] | None = None,
             identity: bool = False) -> tuple[Path, Path]:
    """The analyses (module docstring) of the finished genotyping run in `outdir`, made with the samples sheet at
    `samples`: `trio` (R, S, D) or roles from the data, the identity layer when `identity` (with `threads`).
    Returns the paths of stats.tsv and findings.txt."""
    t0 = time.perf_counter()
    outdir = Path(outdir)
    sheet = read_samples(samples)
    if identity and shutil.which("lastz") is None:
        raise AnalysisError("the identity layer realigns copies with lastz, which is not on PATH; install it with "
                            "the conda environment in environment.yml (bioconda's lastz)")
    _remove_results(outdir)
    log.info("geneflow: analyses of %s: %d genomes, %s%s", outdir.absolute(), len(sheet),
             "roles from the data" if trio is None else f"--trio {','.join(trio)}",
             ", with the identity layer" if identity else "")
    cat = load_catalog(outdir, sheet, Path(samples))
    log.info("geneflow: catalog of %d loci, %d genomes (%.1f s)", len(cat.lid), len(cat.labels),
             time.perf_counter() - t0)
    require("roles", PROPS.roles, ("shared_derived",))
    require("inversions", PROPS.inversions, ("reversed_runs",))
    size = block_bp(cat, PROPS)
    rows, f = run_summary(cat)
    rows.append(StatRow("pav", "jackknife", "all", "all", "block_bp", size, loci="none", note="rule"))
    r, g, tree, dist = relations(cat)
    rows += r
    f.extend(g)
    found, r, g = trios(cat, trio, size, dist)
    rows += r
    f.extend(g)
    for t in found:
        r, g = trio_tests(cat, t, size, dist)
        rows += r
        f.extend(g)
    log.debug("geneflow: relationships and %d trio(s) tested", len(found))
    blocks = inverted_blocks(cat)
    write_inversions(blocks, cat, outdir / "inversions.tsv")
    r, g = inversion_rows(blocks, cat)
    rows += r
    f.extend(g)
    table, r, g = mechanisms(cat, tree)
    write_mechanisms(table, outdir / "mechanisms.tsv")
    rows += r
    f.extend(g)
    log.debug("geneflow: %d inverted block(s), mechanisms", len(blocks))
    self_counts = None
    if identity:
        self_counts = _identity(cat, found, size, threads, dist, rows, f)
    else:
        _no_identity(cat, f)
    r, g = _paintings(cat, found, blocks, self_counts)
    rows += r
    f.extend(g)
    stats, findings = write_results(outdir, rows, f)
    log.info("geneflow: analyses done in %.1f s -> %s, %s", time.perf_counter() - t0, stats.name, findings.name)
    return stats, findings


def screen_analyses(samples: Path, outdir: Path, threads: int, *, trio: tuple[str, str, str] | None = None,
                    records: dict | None = None, routed: bool = False) -> tuple[Path, Path]:
    """The analyses (module docstring, the library screen) of the finished screen in `outdir`, made with the samples
    sheet at `samples`: `trio` (R, S, D) or roles from the data; `records` (label -> routes.library.Records) lends the
    libraries the screen read. Returns the paths of stats.tsv and findings.txt."""
    t0 = time.perf_counter()
    outdir = Path(outdir)
    sheet = read_samples(samples)
    _remove_results(outdir)
    props = SCREEN_PROPS
    require("inversions", props.inversions, ("chain_blocks",))
    run = LibraryRun(outdir, tuple(sheet), threads, {} if records is None else records)
    f = Findings()
    f.add("Run", f"Library screen of {len(run.labels)} genomes ({', '.join(run.labels)}) from their LTR-RT libraries "
                 "alone (the library route's evidence: pairs/, blocks.tsv)")
    if routed:
        f.add("Run", SCREEN_ROUTED)
    f.add("Run", SCREEN_SKIPS)
    size = library_block_bp(run, props)
    rows = [StatRow("library", "jackknife", "all", "all", "block_bp", size, loci="none", note="rule")]
    found, r, g = _screen_trios(run, sheet, trio, props)
    rows += r
    f.extend(g)
    for t in found:
        r = compute(run, (t.R, t.S, t.D), props)
        rows += r
        f.extend(findings(r))
    log.debug("geneflow: %d trio(s) tested", len(found))
    blocks = chain_inversions(outdir / "blocks.tsv", run.seqids)
    write_inversions(blocks, run, outdir / "inversions.tsv")
    r, g = chain_inversion_rows(blocks, outdir / "blocks.tsv", run.labels, run.seqids)
    rows += r
    f.extend(g)
    r, g = _screen_paintings(run, found, blocks, props)
    rows += r
    f.extend(g)
    stats, findings_path = write_results(outdir, rows, f)
    log.info("geneflow: analyses done in %.1f s -> %s, %s", time.perf_counter() - t0, stats.name, findings_path.name)
    return stats, findings_path


def _remove_results(outdir: Path) -> None:
    """Remove an earlier run's stats.tsv, findings.txt and overview figure (overview.png, overview.pdf,
    overview_numbers.tsv) as the analyses start: a run that fails leaves none beside results they would not
    describe."""
    figure = [f"{constants.OVERVIEW_STEM}{end}" for end in constants.FIGURE_SUFFIXES]
    for name in ("stats.tsv", "findings.txt", *figure):
        (Path(outdir) / name).unlink(missing_ok=True)


def _screen_trios(run: LibraryRun, sheet: list[Sample], fixed: tuple[str, str, str] | None, props
                  ) -> tuple[list[Trio], list[StatRow], Findings]:
    """The screen's trios and their roles rows: --trio as given, else every three ingroup genomes of the sheet (sheet
    order) with roles by the route's rule (library_stats.roles)."""
    f = Findings()
    if fixed is not None:
        t = Trio(*fixed, fixed=True)
        f.add("Trio tests", f"Trio {t.scope}: roles fixed by --trio (recipient {t.R}, sister {t.S}, donor {t.D}); "
                            "nothing inferred")
        return [t], [StatRow("library", "roles", t.scope, "all", "source", "--trio", loci="none",
                             note="roles fixed by --trio, as R,S;D")], f
    ingroup = [s.label for s in sheet if s.role != "outgroup"]
    if len(ingroup) < 3:
        f.add("Trio tests", "Trio tests skipped: fewer than three ingroup genomes")
        return [], [], f
    found, rows = [], []
    for genomes in itertools.combinations(ingroup, 3):
        r = roles(run, genomes, props)
        rr, g = role_rows(r, genomes)
        rows += rr
        f.extend(g)
        found.append(Trio(r.R, r.S, r.D, False))
    return found, rows, f


def _screen_paintings(run: LibraryRun, found: list[Trio], blocks: list[InvBlock], props
                      ) -> tuple[list[StatRow], Findings]:
    """Each trio's paintings, one per alphabet of the route, their rows and findings; the painting files written (those
    of the run before removed) and recorded in PAINTINGS."""
    f = Findings()
    if not found:
        f.add("Painting", "Painting skipped: no trio (fewer than three ingroup genomes)")
    paintings: list[Painting] = [paint(run, t, props, alphabet=a, blocks=blocks) for t in found
                                 for a in props.alphabet]
    rows: list[StatRow] = []
    for p in paintings:
        r, g = library_rows(p, run)
        rows += r
        f.extend(g)
    written = write_paintings(paintings, run, run.outdir, _recorded(run.outdir))
    record = run.outdir / PAINTINGS
    record.parent.mkdir(parents=True, exist_ok=True)
    atomic_write(record, json.dumps([p.name for p in written]) + "\n")
    return rows, f


def _identity(cat: Catalog, found: list[Trio], size: int, threads: int, dist, rows: list[StatRow],
              f: Findings) -> dict[tuple[str, str], tuple[int, int]]:
    """The identity layer of each trio: its rows and findings added to `rows` and `f`, identity.tsv and the
    pair-divergence layer written. Returns its self rows as paint's self_counts: (lid, genome) -> (sites, diffs),
    a copy's first row kept (the same in every trio that reads it)."""
    id_rows: list[tuple] = []
    if not found:
        f.add(IDENTITY_SECTION, "Sequence-identity layer skipped: no trio (fewer than three ingroup genomes)")
    for t in found:
        r, g, mine = identity_layer(cat, t, size, threads, dist)
        rows += r
        f.extend(g)
        id_rows += mine
    write_identity(id_rows, cat.outdir / "identity.tsv")
    write_pair_divergences(cat, id_rows)
    counts: dict[tuple[str, str], tuple[int, int]] = {}
    for kind, lid, _, genome, _, sites, diffs, *_ in id_rows:
        if kind == "self":
            counts.setdefault((lid, genome), (sites, diffs))
    return counts


def _no_identity(cat: Catalog, f: Findings) -> None:
    """Without the identity layer: an earlier run's identity.tsv and assembly divergence files deleted (one
    WARNING), and the Dating line on what needs the layer."""
    ident = cat.outdir / "identity.tsv"
    old = (["identity.tsv"] if ident.is_file() else []) + [f"pairs/{p.name}" for p in assembly_divergences(cat.outdir)]
    if old:  # say what a run without --identity throws away
        log.warning("geneflow: deleting %s, the sequence-identity layer of an earlier run with --identity: this run is "
                    "without it (rerun synltr geneflow with --identity to keep it)", ", ".join(old))
        ident.unlink(missing_ok=True)
        write_pair_divergences(cat, [])
    f.add("Dating", _IDENTITY_LINE)


def _paintings(cat: Catalog, found: list[Trio], blocks: list[InvBlock],
               self_counts: dict[tuple[str, str], tuple[int, int]] | None) -> tuple[list[StatRow], Findings]:
    """Each trio's painting, its rows and findings; the painting files written (those of the run before removed,
    PAINTINGS listing them) and recorded in PAINTINGS."""
    f = Findings()
    if not found:
        f.add("Painting", "Painting skipped: no trio (fewer than three ingroup genomes)")
    paintings = [paint(cat, t, PROPS, blocks=blocks, self_counts=self_counts) for t in found]
    rows: list[StatRow] = []
    for p in paintings:
        r, g = paint_rows(p, cat)
        rows += r
        f.extend(g)
    record = cat.outdir / PAINTINGS
    written = write_paintings(paintings, cat, cat.outdir, _recorded(cat.outdir))
    record.parent.mkdir(parents=True, exist_ok=True)
    atomic_write(record, json.dumps([p.name for p in written]) + "\n")
    return rows, f


def _recorded(outdir: Path) -> list[Path]:
    """The painting files in `outdir` the run before wrote (its PAINTINGS record); [] without a readable record."""
    try:
        names = json.loads((outdir / PAINTINGS).read_text())
        return [outdir / n for n in names if isinstance(n, str)]
    except (OSError, ValueError, TypeError):
        return []


def write_results(outdir: Path, rows: list[StatRow], findings: Findings) -> tuple[Path, Path]:
    """findings.txt, then stats.tsv, in `outdir`, each written atomically; returns (stats.tsv, findings.txt).

    findings.txt: stats.report.write_findings in gene flow's layout, its header "synltr <version> findings for the
    run in <outdir, absolute>" and its sections SECTIONS (a section outside SECTIONS is a ValueError, and nothing is
    written). stats.tsv: the header row STATS_COLUMNS, then the rows (stats.report.write_stats)."""
    outdir = Path(outdir)
    stats_path, findings_path = outdir / "stats.tsv", outdir / "findings.txt"
    header = f"{constants.PROGRAM} findings for the run in {os.path.abspath(outdir)}"
    write_findings(findings_path, findings, header=header, sections=SECTIONS)
    with atomic_path(stats_path) as tmp:
        write_stats(rows, tmp)
    return stats_path, findings_path
