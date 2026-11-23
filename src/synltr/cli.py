"""The synltr command line. Imports only the standard library at module level; commands load on demand.

A command's function takes the parsed arguments and returns its exit code, or, when it finished a run in its output
folder, the run's record (calls.manifest.Finished): `main` writes it there as manifest.json, then draws the run's
figure (FIGURES). A run that stops on an error writes none. A command that will write one removes the folder's
earlier manifest.json itself once its own checks have passed (calls.manifest.remove_manifest), never here, so a
usage error leaves a finished run's record as it was."""

import argparse
import importlib
import os
import sys
import traceback
from collections.abc import Mapping
from datetime import UTC, datetime
from pathlib import Path

from synltr.base import constants, log

# subcommand -> "module:function"; the function takes the parsed args and returns the exit code or the run's record
COMMANDS: dict[str, str] = {
    "subgenomes": "synltr.analyses.subgenomes.run:main",
    "genotype": "synltr.calls.genotype:run",
    "geneflow": "synltr.analyses.geneflow.run:main",
    "source": "synltr.analyses.source:main",
    "plot": "synltr.plot.run:main",
    "sim": "synltr.sim.cli:main",
}

# commands that draw their figure in the output folder once their run is recorded (a figure failure only warns)
FIGURES = ("geneflow", "subgenomes")
_INSTALL = ("the conda environment in environment.yml (mamba env create -f environment.yml), or pip install -e . from "
            "the repository")


def _short_cluster_threads(arg: str, following: str | None) -> str | None:
    """The -t value in a short-option cluster such as `-vt3` or `-vt 3`, as argparse reads it; None if absent."""
    for j, c in enumerate(arg[1:], 1):
        if c == "t":
            return arg[j + 1:] or following
        if c != "v":  # any other option takes the rest of the cluster as its value
            return None
    return None


def _at_least(low: int):
    """argparse type: an integer of at least `low`."""

    def parse(text: str) -> int:
        try:
            n = int(text)
        except ValueError:
            raise argparse.ArgumentTypeError(f"invalid integer: {text!r}") from None
        if n < low:
            raise argparse.ArgumentTypeError(f"must be at least {low}")
        return n

    return parse


_threads = _at_least(1)  # argparse type for -t


def _pair(text: str) -> tuple[str, str]:
    """argparse type for --pair: two different labels A,B."""
    a, sep, b = text.partition(",")
    if not (sep and a and b) or "," in b or a == b:
        raise argparse.ArgumentTypeError(f"expected two different labels A,B: {text!r}")
    return a, b


def _trio(text: str) -> tuple[str, str, str]:
    """argparse type for --trio: three labels R,S,D (checked against the sheet by the command)."""
    parts = tuple(text.split(","))
    if len(parts) != 3 or not all(parts):
        raise argparse.ArgumentTypeError(f"expected three labels R,S,D (recipient, sister, donor): {text!r}")
    return parts


def _sheet_option(p: argparse.ArgumentParser, library: bool = False) -> None:
    p.add_argument("sheet", type=Path, metavar="SHEET",
                   help="samples sheet: TSV with the header line 'label fasta annotation' and optional 'group' and "
                        + ("'role', or 'label library' (an LTR_peel prefix, FASTA or FASTA,TABLE) for the screen"
                           if library else
                           "'role'; for read sets 'reads' (a BAM of long reads aligned to the row of role reference)")
                        + "; '.' or empty means none; relative paths resolve against the sheet's folder")


def _paf_option(p: argparse.ArgumentParser) -> None:
    p.add_argument("--paf", type=Path, metavar="DIR",
                   help="anchor each pair with your whole-genome alignments, DIR/<A>__<B>.paf[.gz], instead of "
                        "mapping the flanks of the elements")


def _genotype_options(p: argparse.ArgumentParser) -> None:
    _sheet_option(p)
    p.add_argument("--pair", type=_pair, metavar="A,B",
                   help="only this ordered pair (A annotated), to spread a run over jobs; a run without --pair "
                        "skips the finished pairs and builds the loci")
    _paf_option(p)
    p.add_argument("--ploidy", type=_at_least(1), metavar="N",
                   help=f"read sets: copies of each chromosome per sample (default: {constants.DEFAULT_PLOIDY})")
    p.add_argument("--min-mapq", type=_at_least(0), metavar="Q",
                   help="read sets: lowest mapping quality of a read that anchors at an element (default: "
                        f"{constants.DEFAULT_MIN_MAPQ})")


def _geneflow_options(p: argparse.ArgumentParser) -> None:
    _sheet_option(p, library=True)
    _paf_option(p)
    p.add_argument("--trio", type=_trio, metavar="R,S,D",
                   help="the trio's roles: recipient, sister, donor (default: every trio of ingroup genomes, roles "
                        "from the data)")
    p.add_argument("--identity", action="store_true",
                   help="add the sequence-identity layer: LTR divergence of shared insertions, realigned with "
                        "LASTZ (slow; its alignment is kept for reruns)")
    p.add_argument("--screen", action="store_true",
                   help="compare the sheet's LTR-RT libraries (column library) instead of genotyping assemblies; a "
                        "sheet without genome FASTAs runs it too")


def _subgenomes_options(p: argparse.ArgumentParser) -> None:
    p.add_argument("sheet", type=Path, metavar="SHEET",
                   help="samples sheet with one row: label, library (the LTR-RT FASTA), annotation (LTR_peel GFF3 or a "
                        "K2P table) and optionally fasta (the genome FASTA with its .fai beside it, or the .fai); "
                        "relative paths resolve against the sheet's folder")
    p.add_argument("-n", type=int, metavar="N",
                   help=f"number of subgenomes, 2-{constants.MAX_N} (default: inferred)")
    p.add_argument("--assign", type=Path, metavar="FILE",
                   help="fixed partition instead of phasing: 'sequence label' per line (labels become SG1, SG2... in "
                        "file order)")
    p.add_argument("--mu", type=float, metavar="RATE",
                   help=f"substitutions per site per year for ages in years (default: {constants.MU_DEFAULT:g}, "
                        "flagged in the outputs)")
    p.add_argument("--seed", type=_at_least(0), default=1, metavar="N",
                   help="random seed, 0 or more (default: %(default)s)")


def _plot_options(p: argparse.ArgumentParser) -> None:
    p.add_argument("run", type=Path, metavar="RUN",
                   help="folder of a finished geneflow or subgenomes run (its manifest.json says which), whose figure "
                        "is redrawn there")
    p.add_argument("--trio", type=_trio, metavar="R,S,D",
                   help="the lead trio of a geneflow run's figure (recipient, sister, donor), one the run tested "
                        "(default: the strongest, or on a screen the first)")


def _source_options(p: argparse.ArgumentParser) -> None:
    p.add_argument("run", type=Path, metavar="RUN",
                   help="folder of a finished read genotyping run (synltr genotype on read sets): its "
                        "reads.genotypes.tsv and manifest.json; the tables reads.{sources,copies,locality}.tsv "
                        "are written there")
    p.add_argument("-o", "--outdir", type=Path, metavar="DIR",
                   help="write the tables to DIR instead, with the source run's manifest.json (default: RUN)")


# subcommand -> (one-line help, function adding its own options)
_SUBCOMMANDS = {
    "subgenomes": (("subgenomes of one polyploid genome from its LTR-RT library: homoeologs, split and merger ages, "
                    "exchanges"), _subgenomes_options),
    "genotype": (("presence/absence calls of every annotated element, then loci, VCF and tracks; or, on read sets, "
                  "genotypes of the reference's elements"), _genotype_options),
    "geneflow": (("genotype every pair (or compare libraries: --screen), then relationships, D tests, painting, "
                  "inversions and mechanisms"), _geneflow_options),
    "source": ("likely source copy of each new insertion of a read genotyping run, with the evidence to judge it",
               _source_options),
    "plot": ("redraw the figure of a finished run", _plot_options),
    "sim": ("simulate genomes with known truth (pop, nucleus), or score calls against it (score)", None),
}


def thread_env(argv: list[str], environ: Mapping[str, str] | None = None) -> dict[str, str]:
    """The environment variables that size every thread pool from `-t/--threads` (default 4): numba's pool from -t
    and BLAS and OpenMP single-threaded, whatever `environ` (default os.environ) says, and NUMPY_MADVISE_HUGEPAGE 0
    unless `environ` sets it (a value the user set stands).

    Called before numpy or numba load, so it reads argv itself.
    """
    threads = str(constants.DEFAULT_THREADS)
    for i, arg in enumerate(argv):
        value = None
        if arg in ("-t", "--threads") and i + 1 < len(argv):
            value = argv[i + 1]
        elif arg.startswith("--threads="):
            value = arg.split("=", 1)[1]
        elif arg.startswith("-") and not arg.startswith("--"):
            value = _short_cluster_threads(arg, argv[i + 1] if i + 1 < len(argv) else None)
        if value is not None and value.isdigit() and int(value) >= 1:
            threads = value
    environ = os.environ if environ is None else environ
    return {
        "NUMBA_NUM_THREADS": threads,
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "OMP_WAIT_POLICY": "PASSIVE",
        "KMP_BLOCKTIME": "0",
        "NUMPY_MADVISE_HUGEPAGE": environ.get("NUMPY_MADVISE_HUGEPAGE", "0"),
    }


def _parser() -> argparse.ArgumentParser:
    """The top-level parser with one subparser per command, each carrying the common options (plot: no -o, its run
    folder is both its input and its output; sim: --version, its simulations carry the rest)."""
    version = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    version.add_argument("--version", action="version", version=constants.PROGRAM)
    parser = argparse.ArgumentParser(prog="synltr", description="LTR retrotransposon tracing.", parents=[version],
                                     allow_abbrev=False)
    common = argparse.ArgumentParser(add_help=False, parents=[version], allow_abbrev=False)
    common.add_argument("-t", "--threads", type=_threads, default=constants.DEFAULT_THREADS,
                        help="total CPU budget across all processes (default: %(default)s)")
    common.add_argument("-v", "--verbose", action="store_true", help="log progress and sanity checks")
    outdir = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    outdir.add_argument("-o", "--outdir", type=Path, default=Path("."), help="output folder (default: current folder)")
    sub = parser.add_subparsers(dest="command", metavar="command")
    for name in COMMANDS:
        what, options = _SUBCOMMANDS.get(name, (f"{name} command", None))
        parents = [common] if name in ("plot", "source") else [version] if name == "sim" else [common, outdir]
        p = sub.add_parser(name, parents=parents, allow_abbrev=False, help=what, description=what)
        if options is not None:
            options(p)
        if name == "sim":  # its subcommands (pop, score) carry the common options; their own options live in sim.cli
            from synltr.sim.cli import add_options

            add_options(p, common)
    return parser


def _figure(outdir: Path) -> None:
    """Draw the figure of the run in `outdir`; a failure is one warning, never a change of the exit code."""
    from synltr.plot.style import safe_figure

    @safe_figure
    def draw() -> None:
        from synltr.plot.run import plot

        plot(outdir)

    draw()


def _options(args: argparse.Namespace) -> dict:
    """The parsed options as JSON values: paths absolute, pairs and trios as lists; functions (sim's) left out."""
    out = {}
    for key, value in vars(args).items():
        if not callable(value):
            out[key] = (os.path.abspath(value) if isinstance(value, Path) else list(value) if isinstance(value, tuple)
                        else value)
    return out


def _record(args: argparse.Namespace, argv: list[str], run, started: datetime) -> None:
    """Write manifest.json of the finished run `run` (calls.manifest.record_run): the command line, the parsed
    options and the code fingerprint of the command's module with its imports."""
    from synltr.calls.manifest import record_run

    module = COMMANDS[args.command].split(":")[0]
    unit = module.removeprefix("synltr.") if module.startswith("synltr.") else None
    record_run(run, ["synltr", *argv], _options(args), started, unit)


def main(argv: list[str] | None = None) -> int:
    """Run the command line and return the exit code: 0 done, 1 error, 2 usage, 3 analyses failed, 130 interrupted.
    A command that returns its finished run's record gets manifest.json in that run's folder, then (FIGURES) the
    run's figure."""
    argv = sys.argv[1:] if argv is None else list(argv)
    verbose = False
    try:
        args = _parser().parse_args(argv)
        verbose = getattr(args, "verbose", False)
        log.setup(verbose)
        if args.command is None:
            raise log.UsageError("no command given (see synltr --help)")
        module, func = COMMANDS[args.command].split(":")
        started = datetime.now(UTC).replace(microsecond=0)
        result = getattr(importlib.import_module(module), func)(args)
        if result is None or isinstance(result, int):
            return int(result or 0)
        _record(args, argv, result, started)
        if args.command in FIGURES:
            _figure(Path(result.outdir))
        return 0
    except SystemExit as exc:  # argparse: usage error (2) or --help/--version (0)
        return exc.code if isinstance(exc.code, int) else 0
    except KeyboardInterrupt:  # Ctrl-C: one line (the pool's workers stop without a word)
        print("error: interrupted", file=sys.stderr)
        return 130
    except ImportError as exc:  # a broken environment: name the module and the way to install it
        print(f"error: cannot import the Python module {exc.name or exc}: install synltr's dependencies with "
              f"{_INSTALL}", file=sys.stderr)
        if verbose:
            traceback.print_exc()
        return 1
    except log.SynltrError as exc:
        print(f"error: {exc}", file=sys.stderr)
        if verbose:
            traceback.print_exc()
        return exc.exit_code
    except Exception as exc:  # noqa: BLE001  a command's own failure: one line, traceback with -v
        print(f"error: {type(exc).__name__}: {exc}", file=sys.stderr)
        if verbose:
            traceback.print_exc()
        return 1


def entry() -> None:
    """Console script: size the thread pools from argv before numpy or numba load, then run."""
    os.environ.update(thread_env(sys.argv[1:]))
    sys.exit(main())
