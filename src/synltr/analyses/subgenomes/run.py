"""`synltr subgenomes`: one genome's LTR-RT library and annotation through phasing, homoeologs, dating, derived
biology, painting and what the homoeolog pairs say, into the run's tables and report_data.json (the command line
records the finished run in manifest.json, then draws report.pdf from report_data.json). The settings and inputs are
checked first (usage errors); then the folder's earlier manifest.json, report_data.json and report figure are removed,
so a run that fails leaves no record of other results."""
from __future__ import annotations

import argparse
import time
from pathlib import Path

from synltr.analyses.subgenomes.biology import run_biology, write_biology
from synltr.analyses.subgenomes.dating import run_dating, write_dating
from synltr.analyses.subgenomes.exchanges import run_homoeolog_biology
from synltr.analyses.subgenomes.homoeologs import run_homoeologs
from synltr.analyses.subgenomes.inputs import load_library, log, warn
from synltr.analyses.subgenomes.lineages import load_or_build_clades
from synltr.analyses.subgenomes.painting import run_painting, write_painting
from synltr.analyses.subgenomes.pipeline import phasable, run_phasing
from synltr.analyses.subgenomes.report import (
    REPORT_DATA,
    write_elements,
    write_homoeologs,
    write_phasing,
    write_report_data,
)
from synltr.base.constants import DEFAULT_THREADS, FIGURE_SUFFIXES, MAX_N, REPORT_STEM
from synltr.base.log import UsageError
from synltr.base.runtime import set_numba_threads
from synltr.calls.manifest import Finished, remove_manifest


def options(ltr_fasta, annotation, outdir, *, fai=None, n=None, assign=None, mu=None, label=None,
            threads=DEFAULT_THREADS, seed=1) -> argparse.Namespace:
    """The run's settings as every step reads them: ltr_fasta and annotation (lists of paths), fai (list or None),
    outdir, n_subgenomes, assign, mu (None: the default rate, flagged), label, threads, seed."""
    return argparse.Namespace(ltr_fasta=[str(p) for p in ltr_fasta], annotation=[str(p) for p in annotation],
                              fai=[str(p) for p in fai] if fai else None, outdir=str(outdir), n_subgenomes=n,
                              assign=str(assign) if assign else None, mu=mu, label=label, threads=threads, seed=seed)


def sample_inputs(sample) -> tuple[list[str], list[str], list[str] | None]:
    """A sheet row's library (the LTR-RT FASTA; several comma-separated), annotation and .fai: the fasta column
    names the .fai itself, or the genome FASTA whose .fai sits beside it (no fasta: lengths from the library)."""
    if not sample.library:
        raise UsageError(f"sample {sample.label}: subgenomes needs a library (the LTR-RT FASTA)")
    if not sample.annotation:
        raise UsageError(f"sample {sample.label}: subgenomes needs an annotation (LTRquest GFF3 or a K2P table)")
    fai = None
    if sample.fasta is not None:
        fai = Path(sample.fasta) if str(sample.fasta).endswith(".fai") else Path(f"{sample.fasta}.fai")
        if not fai.is_file():
            raise UsageError(f"sample {sample.label}: no .fai beside {sample.fasta} (samtools faidx makes it; or "
                             f"give the .fai in the fasta column)")
    return sample.library.split(","), [sample.annotation], [str(fai)] if fai else None


def setup(a):
    """Check the settings and inputs before any work (each a UsageError: exit 2, one line), size numba's pool to -t,
    make the output folder."""
    for p in list(a.ltr_fasta) + list(a.annotation) + list(a.fai or []):
        if not Path(p).exists():
            raise UsageError(f"input not found: {p}")
    if a.mu is not None and a.mu <= 0:
        raise UsageError("--mu must be > 0 (substitutions per site per year)")
    if a.threads < 1:
        raise UsageError("--threads must be >= 1")
    if a.n_subgenomes is not None and not 2 <= a.n_subgenomes <= MAX_N:
        raise UsageError(f"-n must be between 2 and {MAX_N} (omit it to infer the structure)")
    if a.assign:
        if not Path(a.assign).exists():
            raise UsageError(f"input not found: {a.assign}")
        if a.n_subgenomes:
            raise UsageError("--assign fixes the partition: -n does not apply")
    import numba
    if a.threads > numba.config.NUMBA_NUM_THREADS:
        warn(f"-t {a.threads} exceeds numba's pool of {numba.config.NUMBA_NUM_THREADS} threads (sized at "
             f"start-up); numba uses {numba.config.NUMBA_NUM_THREADS}, process pools use {a.threads}")
    set_numba_threads(a.threads)
    Path(a.outdir).mkdir(parents=True, exist_ok=True)


def subgenomes(sample, outdir: Path, threads: int, *, n: int | None = None, assign: Path | None = None,
               mu: float | None = None, seed: int = 1) -> int:
    """Run the analysis of one genome (a samples-sheet row: library, annotation, optional .fai; its label names the
    figure) into outdir with threads in all; n fixes the number of subgenomes, assign the partition, mu the
    substitution rate (None: the default, flagged). Returns 0; a setting or input that fails its check is a
    UsageError, whether setup finds it or a later step reading the inputs does (the library and its annotation, the
    .fai, --assign, -n against the library); a runtime failure an exception (exit 1). Once the checks pass, the
    folder's earlier manifest.json, report_data.json and report figure are removed."""
    ltr_fasta, annotation, fai = sample_inputs(sample)
    a = options(ltr_fasta, annotation, outdir, fai=fai, n=n, assign=assign, mu=mu, label=sample.label,
                threads=threads, seed=seed)
    setup(a)
    remove_manifest(Path(a.outdir))                               # never the record or the figure of an earlier run
    for name in (REPORT_DATA, *(f"{REPORT_STEM}{end}" for end in FIGURE_SUFFIXES)):
        (Path(a.outdir) / name).unlink(missing_ok=True)
    t = {}
    t0 = time.time()
    lib = load_library(a)
    sel = phasable(a, lib)                                        # --assign and -n checked before the k-mer stage
    t["inputs_s"] = round(time.time() - t0, 1)
    t0 = time.time()
    cl = load_or_build_clades(a, lib)
    t["kmers_s"] = round(time.time() - t0, 1)
    t0 = time.time()
    ph = run_phasing(a, lib, cl, sel)
    t["phasing_s"] = round(time.time() - t0, 1)
    t0 = time.time()
    try:
        ph["homoeologs"] = run_homoeologs(a, lib, ph)                 # nothing else reads or changes it
    except Exception as e:                                               # noqa: BLE001 (never costs the other tables)
        ph["homoeologs"] = None
        warn(f"homoeolog pairing failed ({type(e).__name__}: {e}); every other output is complete")
    t["homoeologs_s"] = round(time.time() - t0, 1)
    t0 = time.time()
    dt = run_dating(a, lib, ph)
    t["dating_s"] = round(time.time() - t0, 1)
    t0 = time.time()
    bio = run_biology(a, lib, ph, dt)
    t["biology_s"] = round(time.time() - t0, 1)
    t0 = time.time()
    pa = dt["painting"] if dt and dt.get("painting") else run_painting(a, lib, ph, dt)
    t["painting_s"] = round(time.time() - t0, 1)
    t0 = time.time()
    if ph.get("homoeologs"):
        try:
            ph["homoeologs"]["bio"] = run_homoeolog_biology(a, lib, ph, dt, pa)     # reads, changes nothing
        except Exception as e:                                           # noqa: BLE001 (never costs the other tables)
            ph["homoeologs"]["bio"] = None
            warn(f"homoeolog biology failed ({type(e).__name__}: {e}); every other output is complete")
    t["homoeolog_biology_s"] = round(time.time() - t0, 1)
    write_phasing(a, lib, cl, ph, t, dt, bio, pa)
    if dt:
        write_dating(a, dt)
    write_biology(a, lib, bio)
    if pa:
        write_painting(a, pa)
    write_elements(a, lib, ph, dt, bio, pa)
    write_homoeologs(a, lib, ph)
    try:
        write_report_data(a, lib, ph, dt, bio, pa)
    except Exception as e:                                               # noqa: BLE001 (a figure never costs a run)
        warn(f"{REPORT_DATA} failed ({type(e).__name__}: {e}); the tables are complete, report.pdf is not "
             f"drawn")
    log(f"done: {ph['N']} subgenome(s); outputs in {a.outdir}")
    return 0


def main(args: argparse.Namespace) -> Finished:
    """`synltr subgenomes SHEET`: the sheet's one row (label, library, annotation, optional fasta) into args.outdir.
    Returns the run's record (no route; the run computes no input checksum)."""
    from synltr.io.samples import SampleSheetError, read_samples

    try:
        rows = read_samples(args.sheet, min_rows=1)
    except SampleSheetError as e:
        raise UsageError(str(e)) from None
    if len(rows) != 1:
        raise UsageError(f"{args.sheet}: subgenomes takes a sheet with one row (one genome); found {len(rows)}")
    subgenomes(rows[0], Path(args.outdir), args.threads, n=args.n, assign=args.assign, mu=args.mu, seed=args.seed)
    return Finished(Path(args.outdir), (rows[0],))
