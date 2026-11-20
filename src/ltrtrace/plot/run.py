"""The plot command: redraw the figure of a finished run from the files the run wrote, its manifest.json saying which
run it was (calls.manifest). A subgenomes run gets report.pdf, report.png and report_numbers.tsv from its
report_data.json (plot.subgenomes); a gene-flow run gets the overview figure from its tables, with the samples sheet
and the route its manifest records (plot.geneflow), led by the trio --trio names (one the run tested) or else by the
overview's default. A drawing that fails leaves the figure's files as they were
(plot.style.save writes all three under temporary names and renames them into place together at the end), so a
figure on disk is always whole: the earlier one stays until a redraw succeeds."""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

from ltrtrace.base.constants import REPORT_DATA
from ltrtrace.base.log import UsageError
from ltrtrace.calls.manifest import read_manifest

log = logging.getLogger(__name__)

REMEDY = "run ltrtrace geneflow or ltrtrace subgenomes there first"  # what a folder without a finished run needs


def main(args: argparse.Namespace) -> int:
    """`ltrtrace plot RUN [--trio R,S,D]`: args.run, the folder of a finished run; args.trio, the lead trio of a
    gene-flow run's figure (None: the overview's default). Returns the exit code."""
    return plot(Path(args.run), trio=args.trio)


def plot(outdir: Path, *, trio: tuple[str, str, str] | None = None) -> int:
    """Draw the figure of the finished run in `outdir`: report.pdf, report.png and report_numbers.tsv for a
    subgenomes run; overview.png, overview.pdf and overview_numbers.tsv for a gene-flow run, led by `trio` (R, S, D;
    None: the trio with the largest D_pav |Z|, on a screen the first of stats.tsv's). Returns 0; a folder that holds
    no finished run (no manifest.json), or one of another command or without its figure's data, is a UsageError, and
    so is a `trio` the run did not test (plot.geneflow.check_trio) or a `trio` for a subgenomes run; a failure while
    drawing an exception, the figure's files left as they were."""
    outdir = Path(outdir)
    if not outdir.is_dir():
        raise UsageError(f"{outdir} is not a folder")
    doc = read_manifest(outdir, REMEDY)
    command = doc["params"].get("command")
    if command == "subgenomes":
        if trio is not None:
            raise UsageError(f"--trio {','.join(trio)}: {outdir} holds a subgenomes run; --trio names the lead trio "
                             "of a geneflow run's figure")
        if not (outdir / REPORT_DATA).is_file():
            raise UsageError(f"{outdir}: the subgenomes run wrote no {REPORT_DATA} (its log says why); rerun ltrtrace "
                             "subgenomes there")
        from ltrtrace.plot.subgenomes import report

        report(outdir)
        return 0
    if command == "geneflow":
        if not (outdir / "stats.tsv").is_file():
            raise UsageError(f"{outdir} holds no finished gene-flow analyses (no stats.tsv); rerun ltrtrace geneflow "
                             "there")
        from ltrtrace.plot.geneflow import check_trio, overview

        if trio is not None:
            check_trio(outdir, trio)
        overview(outdir, trio=trio)
        return 0
    raise UsageError(f"{outdir} holds a finished `{' '.join(doc.get('command', ['ltrtrace', str(command)])[:2])}` "
                     "run, which draws no figure; ltrtrace plot redraws the figure of a geneflow or subgenomes run")
