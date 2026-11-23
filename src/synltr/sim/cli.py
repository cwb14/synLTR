"""`synltr sim pop`: a population of genomes evolved on real backbone regions, with the truth of every annotated
genome's elements in every other genome (synltr.sim.pop). `synltr sim nucleus`: LTR-RT evolution through
progenitor splits and polyploid mergers, written as an LTRquest-style library with truth (synltr.sim.nucleus).
`synltr sim score`: one ordered pair's calls against truth. `add_options` adds the subcommands' options to the `sim` parser (synltr.cli builds it); `main` runs one."""
from __future__ import annotations

import argparse
import csv
import logging
import math
import re
import textwrap
import time
from dataclasses import replace
from pathlib import Path

from synltr.base.constants import (
    SIM_DEMOG,
    SIM_DROPOUT,
    SIM_FAI_CACHE,
    SIM_JITTER,
    SIM_JITTER_BP,
    SIM_MAX_CHROMS,
    SIM_MAX_DONORS,
    SIM_MAX_SEED,
    SIM_NUCLEUS,
    SIM_PULSES,
    SIM_RATES,
    WINDOW_PAD,
)
from synltr.base.log import SynltrError, UsageError

log = logging.getLogger(__name__)

_REGION = re.compile(r"^(?P<seqid>.+):(?P<start>[0-9][0-9,]*)-(?P<end>[0-9][0-9,]*)$")


class _Lines(argparse.HelpFormatter):
    """Keeps a description's or epilog's line breaks and fills each line to the width on its own; an indented line
    ("  label  text") continues under its text."""

    def _fill_text(self, text: str, width: int, indent: str) -> str:
        out = []
        for line in text.splitlines():
            label = re.match(r" +\S+ +", line)
            out.append(textwrap.fill(line, width, initial_indent=indent, break_on_hyphens=False,
                                     subsequent_indent=indent + " " * (label.end() if label else 0)))
        return "\n".join(out)


def _show(*names: str, of=SIM_DEMOG) -> str:
    """Each field with its default, as the --demography JSON writes it: a number, [a list] or none (null)."""
    def value(v) -> str:
        if v is None:
            return "none"
        if isinstance(v, tuple):
            return f"[{', '.join(f'{x:g}' for x in v)}]"
        return v if isinstance(v, str) else f"{v:g}"
    return ", ".join(f"{n} {value(of.get(n))}" for n in names)


_POP_EPILOG = "\n".join([
    ("defaults: the demography and rates below, small indels at mu/10, and three values measured on real "
     "genomes: t_o from their outgroup conflict, loss_split and replace_mean from their element "
     "losses. A --demography JSON object sets any of these fields by name (a list as a JSON array, none as null)."),
    (f"  Ne          diploid effective sizes: {_show('ne_a', 'ne_b', 'ne_c')}, {_show('ne_bc')} (the ancestor of B "
     f"and C), {_show('ne_o')}, {_show('ne_root')} (the ancestor of BC and O)"),
    (f"  splits      generations ago: {_show('t_ab')} (A from B), {_show('t_bc')} (B from C), {_show('t_o')} "
     "(O from BC)"),
    (f"  pulse       {_show('f')} (the share of the recipient's lineages it moves), {_show('ta')} generations ago, "
     f"{_show('pulse')} (donor>recipient); {_show('per_taxon')} genomes per taxon"),
    (f"  lineages    {_show('sub_a', 'sub_c')}, else [split generations, Ne] of A1's or C1's own population: pulse "
     "C>A1 needs sub_a, C1>A1 and A1>C1 need sub_a and sub_c, and ta must lie below their splits when f > 0"),
    f"  bottleneck  {_show('bottleneck')}, else [Ne, generations] of A right after its split from B",
    (f"  rates       per generation: {_show('mu', of=SIM_RATES)} substitutions per bp; {_show('r', of=SIM_RATES)} and "
     f"{_show('r_low', of=SIM_RATES)} crossovers per bp, r_low between the fractions {_show('low', of=SIM_RATES)} of "
     f"each region; {_show('lam', of=SIM_RATES)} insertions per bp; {_show('delta', of=SIM_RATES)} (1/{1 / SIM_RATES["delta"]:g}) "
     f"losses per element; {_show('indel', of=SIM_RATES)} small indels per bp; {_show('copy_noise', of=SIM_RATES)} "
     "substitutions per base of a new copy"),
    (f"  losses      {_show('loss_split', of=SIM_RATES)}, the shares of solo, del_flank, del_partial and replace; "
     f"{_show('replace_mean', of=SIM_RATES)}, the mean replacement length in bp"),
    "",
    ("truth/insertions.tsv: pos is 0-based, as in BED: the copy goes before base pos of its chromosome's backbone "
     "(0 at its region's START) when host is \".\", else of its host insertion's donor body, so pos is a coordinate "
     "of no output genome; time is in generations ago."),
    "params.json lists the simplifications the simulation makes.",
])


def _flat(d: dict) -> str:
    """A scenario default as the JSON a scenario file writes it, compact."""
    import json

    return json.dumps(d, separators=(",", ":"))


_NUCLEUS_EPILOG = "\n".join([
    ("scenario: a JSON object; a key left out keeps the default below (nested objects are merged key by key, "
     "except tree, which is replaced whole). Times are years before present; the leaves of tree are the "
     "subgenomes."),
    *(f"  {k:<17} {_flat(v)}" for k, v in SIM_NUCLEUS.items()),
    "",
    ("truth/copies.tsv: root identifies the original insertion a copy descends from (copies inherited across a "
     "split or moved by an exchange share it); source is the root of the copy it was transposed from, or \".\" when "
     "the source was a master. LTR offsets in k2p.tsv are 1-based positions along the copy's genome span, counted "
     "from its start (as LTR_peel writes them), the genome-left LTR first, on either strand."),
])


def add_options(parser: argparse.ArgumentParser, common: argparse.ArgumentParser) -> None:
    """Add the `sim` subcommands `pop`, `nucleus` and `score` to `parser`; `common` holds -t, -v and --version."""
    sub = parser.add_subparsers(dest="sim_command", required=True, metavar="SIMULATION")
    p = sub.add_parser(
        "pop", parents=[common], allow_abbrev=False, help="a population of genomes evolved on real backbone regions, with truth",
        formatter_class=_Lines, epilog=_POP_EPILOG,
        description="Evolve haploid genomes A1..An, B1..Bn, C1..Cn (n = --per-taxon) and the outgroup O1 on real "
                    "backbone regions, each region one chromosome (c01, c02, ...) with its annotated elements: an "
                    "msprime genealogy per chromosome (taxa A, B and C, outgroup O, a pulse of gene flow), LTR-RT "
                    "insertions of annotated donors (nested up to two deep), solo-LTRs, flank and partial "
                    "deletions, replacements, SNPs, small indels and assembly gaps. Writes <label>.fa, "
                    "<label>.gff3 for A, B and C (dropout, jitter), samples.tsv (ready for synltr), truth/ "
                    "(<X>__<Y>.tsv for every annotated X and other Y, insertions.tsv, tracts.<label>.bed of the "
                    "recipient taxon, fractions.tsv) and params.json.")
    p.add_argument("--fasta", required=True, type=Path, metavar="FA", help="FASTA holding the backbone regions")
    p.add_argument("--region", required=True, action="append", metavar="SEQ:START-END",
                   help="a backbone region, 1-based inclusive; repeat it for more chromosomes (disjoint, at most "
                        f"{SIM_MAX_CHROMS}, each at least 2 x {WINDOW_PAD} bp plus the longest donor)")
    p.add_argument("--annotation", required=True, metavar="ANNOT",
                   help="the FASTA's LTR-RT annotation (GFF3, BED or LTR_peel prefix): the backbone elements (those "
                        f"wholly inside a region) and up to {SIM_MAX_DONORS} donors (top-level, childless, "
                        "high-confidence where the annotation rates confidence)")
    p.add_argument("-o", "--outdir", required=True, type=Path, metavar="DIR",
                   help=f"output directory: new, or empty but for its own FASTA index cache {SIM_FAI_CACHE}/")
    p.add_argument("--demography", type=Path, metavar="FILE.json",
                   help="a JSON object of population (Demog) and rate (Rates) fields; a field left out keeps its "
                        "default, and the options below override the file")
    p.add_argument("--seed", type=_seed, default=1, metavar="N", help=f"random seed, 0 to {SIM_MAX_SEED} (default: 1)")
    p.add_argument("--f", type=_real, metavar="F",
                   help=f"fraction of the recipient's lineages the pulse moves (default: {SIM_DEMOG["f"]:g}, "
                        "no gene flow)")
    p.add_argument("--ta", type=_real, metavar="GEN", help=f"pulse time, generations ago (default: {SIM_DEMOG["ta"]:g})")
    p.add_argument("--pulse", choices=list(SIM_PULSES), metavar="CODE",
                   help=f"donor>recipient: {', '.join(SIM_PULSES)} (A1, C1: their own lineage populations; default: "
                        f"{SIM_DEMOG["pulse"]})")
    p.add_argument("--per-taxon", type=_positive, metavar="N",
                   help=f"genomes of each of A, B and C (default: {SIM_DEMOG["per_taxon"]})")
    p.add_argument("--bottleneck", type=_ne_gen, metavar="NE,GEN",
                   help="A at effective size NE for GEN generations right after its split from B (default: none)")
    p.add_argument("--dropout", type=_below_one, default=SIM_DROPOUT, metavar="X",
                   help=f"chance an alive element is missing from its genome's annotation, in [0, 1) (default: "
                        f"{SIM_DROPOUT:g})")
    p.add_argument("--dropout-genome", type=_genome_dropout, action="extend", nargs="+", default=[],
                   metavar="LABEL=X", help="an annotated genome's own dropout, in place of --dropout")
    p.add_argument("--jitter", type=_unit, default=SIM_JITTER, metavar="X",
                   help=f"chance an annotated element's ends move by up to {SIM_JITTER_BP} bp, in [0, 1] (default: "
                        f"{SIM_JITTER:g})")
    p.set_defaults(func=_pop)
    n = sub.add_parser(
        "nucleus", parents=[common], allow_abbrev=False,
        help="LTR-RTs through progenitor splits and polyploid mergers, as an LTR_peel-style library with truth",
        formatter_class=_Lines, epilog=_NUCLEUS_EPILOG,
        description="Forward-simulate LTR-retrotransposon families through a tree of progenitor lineages (leaves = "
                    "subgenomes) and the mergers that join their nuclei, with insertions sourced from masters or "
                    "recent copies, copy loss, homoeologous exchanges and per-object substitutions. Writes "
                    "ltr.fa, k2p.tsv, annotation.gff3 and genome.fai in OUTDIR, and copies.tsv (with each copy's "
                    "source), exchanges.tsv, chromosomes.tsv, orthologs.tsv, paralogs.tsv and history.json in "
                    "--truthdir.")
    n.add_argument("scenario", type=Path, metavar="SCENARIO.json", help="scenario file; its keys override the defaults")
    n.add_argument("-o", "--outdir", required=True, type=Path, metavar="DIR", help="library folder (made if missing)")
    n.add_argument("--truthdir", required=True, type=Path, metavar="DIR", help="truth folder (made if missing)")
    n.add_argument("--seed", type=_seed, default=1, metavar="N", help=f"random seed, 0 to {SIM_MAX_SEED} (default: 1)")
    n.set_defaults(func=_nucleus)
    s = sub.add_parser(
        "score", parents=[common], allow_abbrev=False, help="score one ordered pair's calls against its truth file",
        description="Score `synltr genotype`'s calls for X in Y against the simulator's truth: per-state and "
                    "present/absent precision and recall (all calls and high tier), flag agreement, summary "
                    "fractions, and a truth-by-call confusion table. Writes PREFIX.score.tsv and "
                    "PREFIX.confusion.tsv.")
    s.add_argument("--calls", required=True, type=Path, metavar="FILE",
                   help="<outdir>/pairs/<X>__<Y>.calls.tsv.gz (or plain TSV)")
    s.add_argument("--elements", required=True, type=Path, metavar="FILE",
                   help="X's elements, <outdir>/work/elements/<X>.tsv")
    s.add_argument("--truth", required=True, type=Path, metavar="FILE",
                   help="the simulation's truth/<X>__<Y>.tsv, for the same ordered pair X, Y")
    s.add_argument("-o", "--out", required=True, type=Path, metavar="PREFIX", help="output prefix")
    s.set_defaults(func=_score)


def main(args: argparse.Namespace) -> int:
    """Run the chosen `sim` subcommand. Every fault found before simulating or scoring is a usage error, exit 2, one
    `error:` line naming the option or file: a bad --region, --demography file or scenario file (not JSON, or a value
    the simulation cannot use), a --f, --ta, --pulse or --bottleneck value the demography refuses, an unknown
    --dropout-genome label, an --outdir that is not empty, score inputs that do not belong together, and an input
    file the package's readers refuse; other bad input and failed workers are one `error:` line, exit 1."""
    try:
        args.func(args)
    except SynltrError:
        raise
    except (ValueError, RuntimeError, OSError) as exc:
        raise SynltrError(str(exc)) from exc
    return 0


def _positive(text: str) -> int:
    try:
        n = int(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"not an integer: {text!r}") from None
    if n < 1:
        raise argparse.ArgumentTypeError(f"must be at least 1: {n}")
    return n


def _seed(text: str) -> int:
    try:
        n = int(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"not an integer: {text!r}") from None
    if not 0 <= n <= SIM_MAX_SEED:
        raise argparse.ArgumentTypeError(f"must be in [0, {SIM_MAX_SEED}] (one 32-bit word): {n}")
    return n


def _real(text: str) -> float:
    try:
        v = float(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"not a number: {text!r}") from None
    if not math.isfinite(v):
        raise argparse.ArgumentTypeError(f"not a finite number: {text!r}")
    return v


def _ne_gen(text: str) -> tuple[float, float]:
    parts = text.split(",")
    if len(parts) != 2:
        raise argparse.ArgumentTypeError(f"expected NE,GEN: {text!r}")
    return _real(parts[0]), _real(parts[1])


def _below_one(text: str) -> float:
    v = _real(text)
    if not 0 <= v < 1:
        raise argparse.ArgumentTypeError(f"must be in [0, 1): {text!r}")
    return v


def _unit(text: str) -> float:
    v = _real(text)
    if not 0 <= v <= 1:
        raise argparse.ArgumentTypeError(f"must be in [0, 1]: {text!r}")
    return v


def _genome_dropout(text: str) -> tuple[str, float]:
    label, eq, value = text.partition("=")
    if not eq or not label:
        raise argparse.ArgumentTypeError(f"expected LABEL=X: {text!r}")
    try:
        return label, _below_one(value)
    except argparse.ArgumentTypeError as e:
        raise argparse.ArgumentTypeError(f"{label}: {e}") from None


def _parse_region(text: str) -> tuple[str, int, int]:
    """(seqid, start, end) of SEQ:START-END, 1-based inclusive (not yet checked against a FASTA); a UsageError
    naming --region otherwise."""
    m = _REGION.match(text)
    if not m:
        raise UsageError(f"--region {text}: expected SEQ:START-END (1-based inclusive)")
    return m["seqid"], int(m["start"].replace(",", "")), int(m["end"].replace(",", ""))


def _pop(args: argparse.Namespace) -> None:
    from synltr.sim.pop.model import Demog, Rates, demog_from_json
    from synltr.sim.pop.pop import PopParams, check_outdir, check_regions, simulate_pop, write_pop

    t0 = time.perf_counter()
    try:
        demog, rates = demog_from_json(args.demography) if args.demography is not None else (Demog(), Rates())
    except ValueError as exc:  # its message names the file
        raise UsageError(f"--demography {exc}") from None
    given = {k: v for k, v in {"f": args.f, "ta": args.ta, "pulse": args.pulse, "per_taxon": args.per_taxon,
                               "bottleneck": args.bottleneck}.items() if v is not None}
    try:
        demog = replace(demog, **given)  # Demog checks the combination
    except ValueError as exc:
        named = [f"--{k.replace('_', '-')}" for k in given]
        if args.demography is not None:
            named.append(f"--demography {args.demography}")
        raise UsageError(f"{', '.join(named)}: {exc}") from None
    try:
        params = PopParams(demog, rates, dropout=args.dropout, jitter=args.jitter,
                           dropout_genome=tuple(args.dropout_genome))
    except ValueError as exc:  # argparse checked the ranges: what is left is a --dropout-genome label
        raise UsageError(f"--dropout-genome: {str(exc).removeprefix('dropout_genome: ')}") from None
    regions = [_parse_region(r) for r in args.region]
    check_regions(regions)
    try:
        check_outdir(args.outdir)  # before anything is simulated; write_pop checks again
    except ValueError as exc:
        raise UsageError(f"--outdir {exc}") from None
    result = simulate_pop(regions, args.fasta, args.annotation, params, args.seed, args.threads,
                          cache_dir=args.outdir / SIM_FAI_CACHE)
    write_pop(args.outdir, result)
    log.info("done in %.1f s -> %s", time.perf_counter() - t0, args.outdir)


def _nucleus(args: argparse.Namespace) -> None:
    """The scenario file read and every value checked (sim.nucleus.model.Sim) before anything is simulated: a fault
    is a UsageError naming the file."""
    import json

    from synltr.sim.nucleus.model import Sim, merge_defaults
    from synltr.sim.nucleus.write import write_outputs

    t0 = time.perf_counter()
    try:
        scenario = json.loads(args.scenario.read_text())
    except json.JSONDecodeError as exc:
        raise UsageError(f"{args.scenario}: not valid JSON ({exc})") from None
    if not isinstance(scenario, dict):
        raise UsageError(f"{args.scenario}: the scenario must be a JSON object")
    try:
        sim = Sim(merge_defaults(scenario), args.seed)
    except (ValueError, KeyError, TypeError) as exc:
        raise UsageError(f"{args.scenario}: {exc}") from None
    n = write_outputs(sim.run(), args.outdir, args.truthdir)
    log.info("%s seed %d: %d LTR-RTs in %.1f s -> %s, %s", args.scenario.stem, args.seed, n,
             time.perf_counter() - t0, args.outdir, args.truthdir)


def _cell(v) -> str:
    return "NA" if v is None else f"{v:.4f}" if isinstance(v, float) else str(v)


def _write(path: Path, columns: tuple[str, ...], rows: list[dict]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(columns)
        w.writerows([_cell(r[c]) for c in columns] for r in rows)


def _score(args: argparse.Namespace) -> None:
    from synltr.sim.score import CONFUSION_COLUMNS, SCORE_COLUMNS, confusion, score

    log.info("sim score: %s against %s", args.calls, args.truth)
    paths = (args.calls, args.elements, args.truth)
    try:
        rows, conf = score(*paths), confusion(*paths)
    except ValueError as exc:  # files that do not belong together, or one that is not what its option says: named
        raise UsageError(str(exc)) from None
    args.out.parent.mkdir(parents=True, exist_ok=True)
    _write(Path(f"{args.out}.score.tsv"), SCORE_COLUMNS, rows)
    _write(Path(f"{args.out}.confusion.tsv"), CONFUSION_COLUMNS, conf)
    log.info("done -> %s.score.tsv, %s.confusion.tsv", args.out, args.out)
