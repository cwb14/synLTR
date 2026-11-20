"""The overview figure of a finished run: overview.png (300 dpi), overview.pdf (the figure, a short and a
long legend, text kept as text) and overview_numbers.tsv (every number of the legends and the run file it
comes from).

What it shows: insertions the recipient R shares with the donor D, at sites where the sister S has an empty
target site, outnumber the mirror count (S sharing with D at R's empty sites), most of all among young
insertions, and cluster along R's chromosomes beyond what the same painting of S shows. Lineage sorting
gives D = 0 at every age and the same clustering in the sister's mirror. The claim is tested per run, one
row per trio in panel H, never across runs. Everything is drawn from the run's files (stats.tsv, the tracts
and posterior tables, the locus catalog, mechanisms.tsv, inversions.tsv, the pair calls and the FASTAs,
fetched by `io.seqio.Fasta.fetch`, 0-based half-open); no genome is loaded whole.

Panels (letters run in reading order and skip absent panels):
  A: R's sequences to scale, each above its counterpart in S, painted the same way (the decoy). Every
     informative locus is a tick (RD vermillion above the bar, RS blue below); the bar is colored by the
     posterior of D-like ancestry, hatched where a display bin holds no informative locus; tracts are a dark
     line on its top edge (black at q <= Q_TRACT); inverted blocks are gray bands. The posterior is the
     bar's color, not a separate line: a line above the ticks reads as noise at whole-chromosome scale.
  B: the example tract +- one tract length, one row per genome with its state from R's own call.
  C: the example RD locus at base-pair scale: R's and D's copies against S's empty site.
  D: what an absence is: the absence classes of mechanisms.tsv (analyses/geneflow/mechanisms.py) for R's insertions
     in S, by age (log share); the other ingroup pairs as context.
  E: the three sharing counts and the lineage-sorting expectation (equal minority counts).
  F: D_pav overall, by age bin and in the check strata, with 95% intervals; D = 0 drawn.
  G: the mean posterior per drawn sequence pair, decoy to recipient, and the genome-wide means.
  H (more than one trio): every trio as one row, sorted by D_pav, with D among the youngest and oldest
     insertions.
The example is typical, never the most dramatic: tracts (those at q <= Q_TRACT when any, else all) are
ranked by their summed mid-rank percentile distance from the medians of what panel B draws (loci, length,
mean posterior) and taken in that order until one holds a clean candidate. Candidates are the tract's RD
loci with D's call at tier high or D annotated, S's E carrying tsd_match, and R's own TSD known; among the
clean ones the one nearest the median K of all RD loci is drawn. Clean, checked on the FASTAs: R's two TSD
copies identical to each other and to R's annotated TSD; R's element TG...CA; both LTRs whole (the shorter
>= LTR_RATIO of the longer); length inside its family's interquartile range (the superfamily's when the
family has < FAMILY_MIN elements); D's copy with the same two TSD copies and TG...CA; S's single TSD copy
equal to R's TSD and its flanks >= FLANK_IDENT identical to R's over FLANK_CHECK bp on each side (heuristic
thresholds).

Fallbacks, each also a numbers row: no trio (fewer than three ingroup genomes, `no_trio`); no painting for
the lead trio (`no_painting`); no tract (`no_tract`); no FASTA or no clean candidate, so panel C is omitted
with a legend note (`no_example`).

A library screen (`geneflow --screen`; its manifest.json says route library) draws the same layout from the screen's
files (`_ScreenCat` stands for the catalog): the lead trio is --trio, else the first of stats.tsv's roles rows; panels A,
B and G from its recent shared ancestry painting (SCREEN_ALPHABET: posterior.<R>.tsv, the tract files; ticks are the
elements near-identical to D's copy only, above, and to the other sister's only, below (S's on R's rows, R's on the
decoy S's); sequences in natural name order, their lengths the last element's end; each R sequence over the sister's
sequence holding most of its 1:1 orthologs; B's states from R's pair calls, an inferred absence drawn open and keyed
as such, never as an empty site; the typical tract as on the assembly route, with no base-pair example); every other
panel (C, D, E, F; H with several trios) needs observed empty sites or assemblies and is a placeholder saying so
(NEEDS), also a numbers row.

Style: matplotlib (object API, no pyplot state); Liberation Sans or Arial (DejaVu Sans with a warning); PDF
fonttype 42; 180 mm wide; no titles; bold panel letters at each panel's top-left corner; 6-7 pt text.
"""
from __future__ import annotations

import csv
import datetime
import inspect
import logging
import math
import os
import re
import socket
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ltrtrace.analyses.geneflow.catalog import (
    AGE_LABELS,
    STATE_CODE,
    TIER_CODE,
    AnalysisError,
    Catalog,
    load_catalog,
)
from ltrtrace.analyses.geneflow.inversions import syntenic_pairs
from ltrtrace.analyses.geneflow.library_stats import chrom_key
from ltrtrace.analyses.geneflow.paint import ALPHABETS, PAINT_CLASSES
from ltrtrace.analyses.geneflow.paint import classes as paint_classes
from ltrtrace.analyses.geneflow.route import ROUTE_PROPS
from ltrtrace.base.constants import (
    ACC,
    BAND,
    BF_TRACT,
    BINS,
    BLUE,
    DARK,
    FAMILY_MIN,
    FLANK_CHECK,
    FLANK_IDENT,
    FLANK_SHOW,
    FONT_SIZE,
    FONT_SIZE_SMALL,
    GRAY,
    H_FOREST_MIN,
    H_PITCH,
    H_ROWS,
    INNER,
    KEY_LINE_MM,
    LIGHT,
    LTR_RATIO,
    LTR_SHOW,
    LTRC,
    MAX_ROWS,
    MIN_FRAC,
    OVERVIEW_STEM,
    P_WARN,
    Q_TRACT,
    RIGHT_MARGIN_MM,
    SCREEN_ALPHABET,
    Z_SIG,
)
from ltrtrace.base.log import UsageError
from ltrtrace.calls.manifest import read_manifest
from ltrtrace.calls.records import CALL_COLUMNS, stream_calls
from ltrtrace.calls.states import ALLELES, INFERRED_ABSENT, PAIR_STATES, PRESENT, call_class
from ltrtrace.io.elements import Element, read_elements
from ltrtrace.io.samples import PAIR_SEP, read_samples
from ltrtrace.io.seqio import Fasta, revcomp
from ltrtrace.plot import style
from ltrtrace.stats.report import read_stats

log = logging.getLogger(__name__)

FS, FS_SMALL = FONT_SIZE, FONT_SIZE_SMALL
STEM = OVERVIEW_STEM
WIDTH_MM = style.WIDTH_MM
RIGHT = WIDTH_MM - RIGHT_MARGIN_MM  # the right edge (mm) shared by the data areas of A, B and D: one grid
RD, RS, _, _, MISSING, _ = range(len(PAINT_CLASSES))
REMEDY = "run ltrtrace geneflow there first"  # what a folder without a finished gene-flow run needs
_PRESENT = tuple(STATE_CODE[x] for x in sorted(PRESENT, key=ALLELES.index))
_ABSENCE = (("frac_E", "empty site"), ("frac_D_inside+frac_O_inside", "cut inside, remnant shown"),
            ("frac_D_inside_unshown+frac_O_inside_unshown", "cut inside, remnant not shown"),
            ("frac_D_beyond", "deletion past the ends"), ("frac_O_site", "other sequence at site"),
            ("frac_O_complex", "complex other sequence"))
_AGE_GRAYS = ("#1A1A1A", "#7A7A7A", "#C4C4C4")


# ---------------------------------------------------------------- the entry point


def overview(outdir: Path, sheet: Path | None = None, trio: tuple[str, str, str] | None = None, *,
             cat: Catalog | None = None) -> list[Path]:
    """Draw the overview of the finished gene-flow run in `outdir`: [overview.png, overview.pdf,
    overview_numbers.tsv]. sheet: the samples sheet the run was made with (default: the one its manifest.json
    records; a run without one is a UsageError).
    trio: the lead trio (R, S, D), else the trio with the largest D_pav |Z| in stats.tsv. cat: the run's catalog
    when already loaded (else it is loaded here)."""
    import matplotlib as mpl

    outdir = Path(outdir)
    if not (outdir / "stats.tsv").is_file():
        raise AnalysisError(f"no stats.tsv in {outdir}: run ltrtrace geneflow first")
    if cat is None and _screen_run(outdir):
        return _screen_overview(outdir, trio)
    if cat is None:
        sheet = _sheet(outdir) if sheet is None else Path(sheet)
        cat = load_catalog(outdir, read_samples(sheet), sheet)
    run = _Run(outdir, cat)
    nums = style.Numbers()
    style.setup_fonts()
    mono = style.rc_params()["mathtext.tt"]
    model = _model(run, nums, trio)
    with mpl.rc_context(style.rc_params()):
        fig = _draw(model, mono)
        short, long_ = _legends(model, nums)
        style.save(fig, outdir / STEM, nums, legend=[short, long_])
    log.info("overview: %s.png, %s.pdf, %s_numbers.tsv", STEM, STEM, STEM)
    return style.figure_files(outdir / STEM)


def _sheet(outdir: Path) -> Path:
    """The samples sheet the run in `outdir` was made with, as its manifest.json records it; a UsageError when the
    manifest names none or the sheet is gone."""
    sheet = read_manifest(outdir, REMEDY)["params"].get("sheet")
    if not isinstance(sheet, str):
        raise UsageError(f"the manifest.json in {outdir} names no samples sheet; {REMEDY}")
    if not Path(sheet).is_file():
        raise UsageError(f"the samples sheet of the run in {outdir}, {sheet}, is gone; put it back, or {REMEDY}")
    return Path(sheet)


# ---------------------------------------------------------------- the run's files


class _Run:
    """The run's result files, read once: stats.tsv indexed by (test, scope, stratum, statistic)."""

    def __init__(self, outdir: Path, cat: Catalog):
        self.outdir, self.cat = outdir, cat
        self._tables: dict[str, list[dict[str, str]] | None] = {}
        self.stats = read_stats(outdir / "stats.tsv")
        self.by: dict[tuple[str, str, str, str], object] = {}
        for r in self.stats:
            self.by.setdefault((r.test, r.scope, r.stratum, r.statistic), r)

    def row(self, test: str, scope: str, stratum: str, statistic: str):
        return self.by.get((test, scope, stratum, statistic))

    def value(self, test: str, scope: str, stratum: str, statistic: str, default=math.nan):
        r = self.row(test, scope, stratum, statistic)
        return default if r is None else r.value

    def table(self, name: str) -> list[dict[str, str]] | None:
        """A TSV result file of the run as rows (header keys; a leading # dropped); None if absent. Parsed once."""
        if name not in self._tables:
            path = self.outdir / name
            rows = None
            if path.is_file():
                with open(path, encoding="utf-8", newline="") as fh:
                    rd = csv.reader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)
                    header = [h.lstrip("#") for h in next(rd, [])]
                    rows = [dict(zip(header, row, strict=False)) for row in rd if row]
            self._tables[name] = rows
        return self._tables[name]


def _src(test: str, scope: str, stratum: str, statistic: str) -> str:
    return f"stats.tsv:{test}/{scope}/{stratum}/{statistic}"


@dataclass(frozen=True)
class _Trio:
    R: str
    S: str
    D: str

    @property
    def scope(self) -> str:
        return f"{self.R},{self.S};{self.D}"


def _trios(run: _Run) -> list[_Trio]:
    """The run's trios, from stats.tsv's roles rows (stats order)."""
    return _tested(run.stats)


def _tested(stats) -> list[_Trio]:
    """The trios of stats.tsv's rows `stats`: its roles rows' scopes R,S;D (stats order)."""
    out = []
    for r in stats:
        if r.test == "roles" and r.statistic == "source":
            m = re.fullmatch(r"([^,;]+),([^,;]+);([^,;]+)", r.scope)
            if m and _Trio(*m.groups()) not in out:
                out.append(_Trio(*m.groups()))
    return out


def check_trio(outdir: Path, trio: tuple[str, str, str]) -> None:
    """UsageError unless `trio` (plot --trio R,S,D) is a trio the finished gene-flow run in `outdir` tested: three
    ingroup labels of its samples sheet, as geneflow's --trio must be (analyses.geneflow.run.check_trio_option), in
    the roles of one of its stats.tsv's roles rows. The error line names the sheet or the run's trios."""
    from ltrtrace.analyses.geneflow.run import check_trio_option

    sheet = _sheet(outdir)
    check_trio_option(read_samples(sheet), trio, sheet)
    tested = _tested(read_stats(outdir / "stats.tsv"))
    if _Trio(*trio) not in tested:
        have = "; ".join(f"{t.R},{t.S},{t.D}" for t in tested) or "none"
        raise UsageError(f"--trio {','.join(trio)}: not a trio this run tested (its trios: {have}); rerun ltrtrace "
                         f"geneflow with --trio {','.join(trio)} to test it")


def _lead(run: _Run, trios: list[_Trio], trio: tuple[str, str, str] | None) -> _Trio | None:
    """--trio as given (it must be one of the run's trios), else the largest D_pav |Z| (nan last, stats
    order on ties); None without a trio."""
    if trio is not None:
        want = _Trio(*trio)
        if want not in trios:
            have = ", ".join(t.scope for t in trios) or "none"
            raise AnalysisError(f"--trio {','.join(trio)}: not a trio of this run's stats.tsv (trios: {have}); "
                                "rerun ltrtrace geneflow with --trio")
        return want
    if not trios:
        return None
    z = [abs(_num(run.row("D_pav", t.scope, "all", "D"), "z")) for t in trios]
    best = max(range(len(trios)), key=lambda j: (z[j] if math.isfinite(z[j]) else -1.0, -j))
    return trios[best]


def _num(row, attr: str = "value") -> float:
    return math.nan if row is None else _float(getattr(row, attr))


@dataclass
class _Painting:
    """One painted genome X of the lead trio (R, or S as the decoy): its painted loci in (fai index, start)
    order and its tracts. file: the posterior file's name (None when the painting was not run: classes from
    paint.classes, no posterior)."""

    X: str
    lid: np.ndarray
    idx: np.ndarray            # catalog indexes
    seq: np.ndarray            # fai index of X
    start: np.ndarray
    end: np.ndarray
    cls: np.ndarray            # PAINT_CLASSES index
    post: np.ndarray           # nan without a painting
    tracts: list[dict]
    file: str | None
    tract_file: str | None


def _painting(run: _Run, trio: _Trio, decoy: bool) -> _Painting | None:
    """The lead trio's painting of R (decoy: of S) from its posterior and tract files; None if absent."""
    cat = run.cat
    X, Y = (trio.S, trio.R) if decoy else (trio.R, trio.S)
    stem = "decoy_" if decoy else ""
    name = f"{stem}posterior.{X}.tsv"
    rows = run.table(name)
    if rows is None:
        return None
    rows = [r for r in rows if r.get("sister") == Y and r.get("donor") == trio.D]
    if not rows:
        return None
    order = {n: j for j, n in enumerate(cat.seqids[X])}
    lid_index = _lid_index(cat)
    lid = np.array([r["lid"] for r in rows])
    idx = np.array([lid_index[x] for x in lid], dtype=np.int64)
    seq = np.array([order[r["seqid"]] for r in rows], dtype=np.int64)
    start = np.array([int(r["start"]) for r in rows], dtype=np.int64)
    end = np.array([int(r["end"]) for r in rows], dtype=np.int64)
    code = {c: j for j, c in enumerate(PAINT_CLASSES)}
    cls = np.array([code[r["class"]] for r in rows], dtype=np.int8)
    post = np.array([_float(r["post"]) for r in rows])
    o = np.lexsort((start, seq))
    tname = f"{stem}tracts.{X}.bed"
    tracts = []
    for r in run.table(tname) or []:
        if r.get("sister") != Y or r.get("donor") != trio.D:
            continue
        tracts.append({"seq": order[r["chrom"]], "seqid": r["chrom"], "start": int(r["start0"]) + 1,
                       "end": int(r["end"]), "n_loci": int(r["n_loci"]), "n_RD": int(r["n_RD"]),
                       "n_RS": int(r["n_RS"]), "mean_post": _float(r["mean_post"]), "q": _float(r.get("q", ".")),
                       "youngest_k": _float(r.get("youngest_k", "."))})
    return _Painting(X, lid[o], idx[o], seq[o], start[o], end[o], cls[o], post[o], tracts, name, tname)


def _classes_only(cat: Catalog, trio: _Trio, decoy: bool) -> _Painting:
    """Without a painting: X's loci and classes from paint.classes (the painting's own class rule)."""
    X, Y = (trio.S, trio.R) if decoy else (trio.R, trio.S)
    idx, cls = paint_classes(cat, X, Y, trio.D, alphabet=ROUTE_PROPS["assembly"].alphabet[0])
    g = cat.g(X)
    return _Painting(X, cat.lid[idx], idx, cat.seq[idx, g].astype(np.int64), cat.start[idx, g], cat.end[idx, g],
                     cls, np.full(len(idx), np.nan), [], None, None)


def _float(text: str | None) -> float:
    try:
        return float(text)
    except (TypeError, ValueError):
        return math.nan


def _lid_index(cat: Catalog) -> dict[str, int]:
    """lid -> catalog index, cached on the catalog object (never by id(): a freed catalog's id is reused)."""
    idx = getattr(cat, "_overview_lid_index", None)
    if idx is None:
        idx = {x: i for i, x in enumerate(cat.lid.tolist())}
        cat._overview_lid_index = idx
    return idx


# ---------------------------------------------------------------- what is drawn


@dataclass
class _Row:
    """One drawn sequence of a painted genome in panel A."""

    genome: str
    seq: int
    name: str
    label: str
    length: int
    rd: np.ndarray             # positions of RD loci (1-based starts)
    rs: np.ndarray
    edges: np.ndarray          # display bin edges (bp)
    post: np.ndarray           # mean posterior per bin (nan: no painted locus)
    gaps: list[tuple[int, int]]     # stretches of bins without an informative locus
    tracts: list[tuple[int, int, bool]]   # (start, end, q <= Q_TRACT)
    bands: list[tuple[int, int]]   # inverted blocks
    mean_post: float = math.nan


@dataclass
class _Model:
    """Everything the figure and its legends show."""

    run: _Run
    trio: _Trio | None
    trios: list[_Trio]
    groups: dict[str, str | None]
    rows: list[tuple[_Row, _Row | None]] = field(default_factory=list)
    n_drawn: int = 0
    n_not: int = 0
    bin_bp: int = 0
    xmax: int = 0
    target: _Painting | None = None
    decoy: _Painting | None = None
    painted: bool = False
    window: dict | None = None
    example: dict | None = None
    example_note: str = ""
    panels: list[str] = field(default_factory=list)
    letters: dict[str, str] = field(default_factory=dict)
    no_trio: dict | None = None
    text: dict = field(default_factory=dict)
    named: bool = False             # the lead trio was named by the user (--trio), not chosen by its |Z|
    screen: bool = False            # a library screen's run: its painting rows and placeholders
    placeholders: set = field(default_factory=set)   # panels drawn as placeholders (a screen run lacks their data)


def _model(run: _Run, nums: style.Numbers, trio: tuple[str, str, str] | None) -> _Model:
    cat = run.cat
    trios = _trios(run)
    lead = _lead(run, trios, trio)
    m = _Model(run=run, trio=lead, trios=trios, groups={s.label: s.group for s in cat.samples},
               named=trio is not None)
    if lead is None:
        _no_trio_model(m, nums)
    else:
        _trio_model(m, nums)
    m.letters = {p: "ABCDEFGHIJ"[j] for j, p in enumerate(m.panels)}
    return m


def _drawn(cat: Catalog, X: str) -> tuple[list[int], int]:
    """X's sequences drawn in panel A (fai indexes in fai order: the MAX_ROWS longest of those >= MIN_FRAC of
    the longest, ties by fai order) and how many are not drawn."""
    lengths = cat.lengths[X]
    if not len(lengths):
        return [], 0
    ok = np.flatnonzero(lengths >= MIN_FRAC * lengths.max())
    ok = ok[np.lexsort((ok, -lengths[ok]))][:MAX_ROWS]
    return sorted(int(s) for s in ok), int(len(lengths) - len(ok))


def _partners(cat: Catalog, X: str, Y: str, drawn: list[int]) -> dict[int, int | None]:
    """Each drawn X sequence's counterpart in Y: of its syntenic partners (`inversions.syntenic_pairs`, the
    inversion caller's own definition of homologous sequences), the one holding most of the loci both place
    (merged sites, any discovery), each Y sequence used once (X order; ties by Y's fai order); None when it has
    no partner left."""
    by: dict[int, list[tuple[int, int]]] = {}
    for (s, t), c in syntenic_pairs(cat, X, Y).items():
        by.setdefault(s, []).append((-c, t))
    used: set[int] = set()
    out: dict[int, int | None] = {}
    for s in drawn:
        out[s] = next((y for _, y in sorted(by.get(s, [])) if y not in used), None)
        if out[s] is not None:
            used.add(out[s])
    return out


def _chromosome_scale(cat: Catalog, drawn: dict[str, list[int]]) -> bool:
    """Whether every drawn sequence is chromosome-scale by the drawing rule, so that a legend may
    call them chromosomes (a name is never read): for each genome drawn, its sequences at least MIN_FRAC of its
    longest number at most MAX_ROWS (the rule could draw them all: they are the genome's large sequences, not the
    largest of many scaffolds), and every sequence drawn of it is one of them."""
    for X, seqs in drawn.items():
        lengths = cat.lengths[X]
        if not len(lengths):
            return False
        big = set(np.flatnonzero(lengths >= MIN_FRAC * lengths.max()).tolist())
        if len(big) > MAX_ROWS or not set(seqs) <= big:
            return False
    return True


def _nice_floor(x: float) -> int:
    """The largest 1, 2 or 5 x 10^k <= x (at least 1)."""
    if x < 1:
        return 1
    e = 10 ** math.floor(math.log10(x))
    return int(max(m * e for m in (1, 2, 5) if m * e <= x))


def _short_names(names: list[str]) -> list[str]:
    """Labels: names without the longest common prefix that ends in a delimiter (_ . | - :), when every name
    keeps some text."""
    if len(names) < 2:
        cut = names[0].rfind("_") + 1 if names and "_" in names[0] else 0
        return [n[cut:] or n for n in names]
    pre = os.path.commonprefix(names)
    cut = max((pre.rfind(c) for c in "_.|-:"), default=-1) + 1
    return [n[cut:] if cut and all(len(x) > cut for x in names) else n for n in names]


def _row(p: _Painting, cat: Catalog, s: int, label: str, bin_bp: int, inv: dict[int, list[tuple[int, int]]],
         weighted: bool = False) -> _Row:
    """Panel A's drawing data of sequence s of painting p; its mean posterior weighted by element length when
    `weighted` (the library route's posterior mean)."""
    X = p.X
    length = int(cat.lengths[X][s])
    on = p.seq == s
    start, cls, post = p.start[on], p.cls[on], p.post[on]
    nb = max(1, -(-length // bin_bp))
    b = np.clip((start - 1) // bin_bp, 0, nb - 1)
    finite = np.isfinite(post)
    tot = np.bincount(b[finite], weights=post[finite], minlength=nb)
    cnt = np.bincount(b[finite], minlength=nb)
    mean = np.where(cnt > 0, tot / np.maximum(cnt, 1), np.nan)
    informative = np.bincount(b[cls <= RS], minlength=nb)
    gaps, j = [], 0
    while j < nb:
        if informative[j] == 0:
            k = j
            while k + 1 < nb and informative[k + 1] == 0:
                k += 1
            gaps.append((j * bin_bp, min((k + 1) * bin_bp, length)))
            j = k + 1
        else:
            j += 1
    edges = np.minimum(np.arange(nb + 1) * bin_bp, length)
    tracts = [(t["start"], t["end"], math.isfinite(t["q"]) and t["q"] <= Q_TRACT) for t in p.tracts
              if t["seq"] == s]
    if not finite.any():
        mean_post = math.nan
    elif weighted:
        mean_post = float(np.average(post[finite], weights=(p.end[on] - start + 1)[finite]))
    else:
        mean_post = float(np.mean(post[finite]))
    return _Row(X, s, cat.seqids[X][s], label, length, start[cls == RD], start[cls == RS], edges, mean, gaps,
                tracts, inv.get(s, []), mean_post)


def _inverted(run: _Run, X: str, partners: set[str], bin_bp: int, tally: dict
               ) -> dict[int, list[tuple[int, int]]]:
    """X's inverted blocks against `partners` (inversions.tsv), by X's fai index: those whose loci average at
    least one per display bin over X's span of the block (a display-resolution rule: a sparser block is mostly
    unsupported area at this scale, typically a few scattered loci chained over tens of Mb). tally counts the
    blocks drawn and not drawn, and the largest not drawn (span bp, loci)."""
    cat = run.cat
    order = {n: j for j, n in enumerate(cat.seqids[X])}
    out: dict[int, list[tuple[int, int]]] = {}
    for r in run.table("inversions.tsv") or []:
        for me, other in (("x", "y"), ("y", "x")):
            if r.get(f"genome_{me}") == X and r.get(f"genome_{other}") in partners:
                s = order.get(r.get(f"seqid_{me}", ""))
                a, b, n = _float(r.get(f"start_{me}")), _float(r.get(f"end_{me}")), _float(r.get("n_loci"))
                if s is None or not (math.isfinite(a) and math.isfinite(b) and math.isfinite(n)):
                    continue
                if n * bin_bp >= b - a + 1:
                    out.setdefault(s, []).append((int(a), int(b)))
                    tally["drawn"] = tally.get("drawn", 0) + 1
                else:
                    tally["sparse"] = tally.get("sparse", 0) + 1
                    if b - a + 1 > tally.get("largest", (0, 0))[0]:
                        tally["largest"] = (int(b - a + 1), int(n))
    return out


# ---------------------------------------------------------------- the trio model


def _trio_model(m: _Model, nums: style.Numbers) -> None:
    run, cat, t = m.run, m.run.cat, m.trio
    # a painting is drawn only when stats.tsv holds its current rows: assemble leaves a stale paint part out
    current = run.row("paint", t.scope, "all", "posterior_mean") is not None
    m.target = _painting(run, t, decoy=False) if current else None
    m.decoy = _painting(run, t, decoy=True) if current and m.target is not None else None
    m.painted = m.target is not None
    if m.target is None:
        m.target = _classes_only(cat, t, decoy=False)
        why = ("no current painting for this trio" if cat.annotated.get(t.R) else
               f"no painting: {t.R} has no annotation, so it discovers no insertions to paint")
        nums.add("no_painting", why, _src("paint", t.scope, "all", "posterior_mean"))
    if m.decoy is None and cat.annotated.get(t.S) and not m.painted:
        m.decoy = _classes_only(cat, t, decoy=True)
    drawn, m.n_not = _drawn(cat, t.R)
    m.n_drawn = len(drawn)
    part = _partners(cat, t.R, t.S, drawn) if m.decoy is not None else {s: None for s in drawn}
    lens = [int(cat.lengths[t.R][s]) for s in drawn] + [int(cat.lengths[t.S][y]) for y in part.values()
                                                        if y is not None]
    m.xmax = max(lens) if lens else 1
    m.bin_bp = _nice_floor(m.xmax / BINS)
    m.text["inv"] = tally = {}
    inv_r, inv_s = _inverted(run, t.R, {t.S, t.D}, m.bin_bp, tally), _inverted(run, t.S, {t.R, t.D}, m.bin_bp, tally)
    rl = _short_names([cat.seqids[t.R][s] for s in drawn])
    pl = _short_names([cat.seqids[t.S][y] for y in part.values() if y is not None])
    pl_it = iter(pl)
    for s, lab in zip(drawn, rl, strict=True):
        r = _row(m.target, cat, s, f"{t.R} {lab}", m.bin_bp, inv_r)
        y = part[s]
        d = _row(m.decoy, cat, y, f"{t.S} {next(pl_it)}", m.bin_bp, inv_s) if y is not None else None
        m.rows.append((r, d))
    shown = {t.R: drawn} | ({t.S: [y for y in part.values() if y is not None]} if m.decoy is not None else {})
    m.text["unit"] = "chromosome" if _chromosome_scale(cat, shown) else "sequence"
    nums.add("sequences_drawn", m.n_drawn, f"genome index:{t.R} sequence lengths")
    nums.add("sequences_not_drawn", m.n_not, f"genome index:{t.R} sequence lengths")
    _choose_example(m, nums)
    m.panels = ["A"] + (["B"] if m.window else []) + (["C"] if m.example else [])
    if _absence_rows(run, t):
        m.panels.append("D")
    m.panels += ["E", "F"] + (["G"] if m.painted and m.decoy is not None else [])
    if len(m.trios) > 1:
        m.panels.append("H")


# ---------------------------------------------------------------- the example tract and locus


def _midpct(x: np.ndarray, among: np.ndarray | None = None) -> np.ndarray:
    """Mid-rank percentile of each value of x among `among` (default x itself): (#less + #equal / 2) / n."""
    s = np.sort(x if among is None else among)
    lo, hi = np.searchsorted(s, x, "left"), np.searchsorted(s, x, "right")
    return (lo + (hi - lo) / 2) / max(len(s), 1)


def _pct_of(v: float, values: np.ndarray) -> float:
    """The mid-rank percentile (0-1) of v among `values`; nan without values."""
    values = values[np.isfinite(values)]
    if not len(values) or not math.isfinite(v):
        return math.nan
    return float(_midpct(np.array([v]), values)[0])


_TRACT_KEYS = ("n_loci", "length", "mean_post")   # what panel B shows of a tract: its loci, its span, its plateau


def _tract_value(x: dict, key: str) -> float:
    return float(x["end"] - x["start"] + 1) if key == "length" else float(x[key])


def _choose_example(m: _Model, nums: style.Numbers) -> None:
    """The example tract and locus (module docstring), with their structural checks; sets m.window, m.example
    and m.example_note. Tracts are ranked by their summed mid-rank distance from the medians of the three things
    panel B shows (loci, length, mean posterior; Q_TRACT tracts first) and taken in that order until one holds a
    clean R+D insertion."""
    cat, t, p = m.run.cat, m.trio, m.target
    gR = cat.g(t.R)
    rd = p.cls == RD
    k_all = cat.k[p.idx[rd], gR]
    k_ok = k_all[np.isfinite(k_all)]
    med_k = float(np.median(k_ok)) if len(k_ok) else math.nan
    tracts = p.tracts
    order = _typical(m, tracts)
    checker = _Checker(m)
    inside = {j: np.flatnonzero((p.seq == tracts[j]["seq"]) & (p.start >= tracts[j]["start"])
                                & (p.start <= tracts[j]["end"]) & rd) for j in order}
    checker.prefetch([i for js in inside.values() for i in js.tolist()])   # one pass over the pair file
    seen = screened = passed = 0
    chosen = None
    for j in order:
        seen += 1
        cands = checker.clean(inside[j], med_k)
        screened += len(inside[j])
        if cands:
            chosen, passed = (tracts[j], cands[0]), len(cands)
            break
    if chosen is None and not tracts:
        allrd = np.flatnonzero(rd)
        if len(allrd):
            checker.prefetch(allrd.tolist())
            cands = checker.clean(allrd, med_k)
            screened = len(allrd)
            if cands:
                chosen, passed = (None, cands[0]), len(cands)
    m.text.update(med_k=med_k, n_rd=int(rd.sum()), tracts_seen=seen, screened=screened, passed=passed,
                  checker=checker)
    if chosen is not None:
        x, ex = chosen
        m.example = ex
        m.window = _window(m, x if x is not None else _locus_window(m, ex), x)
    else:
        m.example_note = checker.why_none(bool(tracts))
        if tracts and order:
            m.window = _window(m, tracts[order[0]], tracts[order[0]])
        elif rd.any():
            j = int(np.flatnonzero(rd)[np.argmin(np.abs(np.nan_to_num(k_all, nan=np.inf) - med_k))]
                    if len(k_ok) else np.flatnonzero(rd)[0])
            m.window = _window(m, _locus_window(m, {"j": j}), None)
    if m.window is not None and m.window.get("tract") is not None:
        x = m.window["tract"]
        m.window["pct"] = {k: _pct_of(_tract_value(x, k), np.array([_tract_value(y, k) for y in tracts]))
                           for k in _TRACT_KEYS}
        dec = m.decoy.tracts if m.decoy is not None else []
        m.window["pct"]["decoy_mean_post"] = _pct_of(x["mean_post"], np.array([y["mean_post"] for y in dec]))
    if not tracts:
        nums.add("no_tract", "no tract called", f"{p.tract_file or 'stats.tsv'}:"
                 + ("n_loci" if p.tract_file else f"paint/{t.scope}/all/n_tracts"))
    if m.example is None:
        nums.add("no_example", m.example_note, f"{p.file or 'loci.tsv'}:class")


def _typical(m: _Model, tracts: list[dict]) -> list[int]:
    """The tracts in the example's order (module docstring): those at q <= Q_TRACT when any, else all, by their
    summed mid-rank distance from the medians of what panel B shows; m.text records the pool and its medians."""
    q05 = [j for j, x in enumerate(tracts) if math.isfinite(x["q"]) and x["q"] <= Q_TRACT]
    pool = q05 or list(range(len(tracts)))
    m.text["pool_q05"] = bool(q05)
    if not pool:
        return []
    cols = {k: np.array([_tract_value(tracts[j], k) for j in pool]) for k in _TRACT_KEYS}
    dist = sum(np.abs(_midpct(np.nan_to_num(v, nan=np.nanmedian(v) if np.isfinite(v).any() else 0.0)) - 0.5)
               for v in cols.values())
    m.text["tract_medians"] = ({k: float(np.nanmedian(v)) for k, v in cols.items()}, len(pool))
    return [pool[j] for j in np.lexsort((np.arange(len(pool)), dist))]


def _locus_window(m: _Model, ex: dict) -> dict:
    """A window around one locus when there is no tract: +- 5 median informative spacings (at least 20 kb)."""
    p = m.target
    j = ex["j"]
    sp = _num(m.run.row("paint", m.trio.scope, "all", "spacing_median"))
    half = int(max(5 * sp, 20_000)) if math.isfinite(sp) else 100_000
    return {"seq": int(p.seq[j]), "start": int(p.start[j]) - half, "end": int(p.end[j]) + half, "locus": True}


def _window(m: _Model, w: dict, tract: dict | None) -> dict:
    """Panel B's window: the tract (or a locus window) +- one tract length, clipped to the sequence."""
    cat, t = m.run.cat, m.trio
    seq = int(w["seq"])
    length = int(cat.lengths[t.R][seq])
    if tract is not None:
        pad = tract["end"] - tract["start"] + 1
        a, b = tract["start"] - pad, tract["end"] + pad
    else:
        a, b = w["start"], w["end"]
    return {"seq": seq, "start": max(1, int(a)), "end": min(length, int(b)), "tract": tract,
            "name": cat.seqids[t.R][seq]}


class _Checker:
    """The example's candidates and structural checks (module docstring). Each failed check is counted under a
    plain label (the long legend lists the commonest)."""

    def __init__(self, m: _Model):
        self.m = m
        self.run, self.cat, self.t = m.run, m.run.cat, m.trio
        self.gR, self.gS, self.gD = self.cat.g(self.t.R), self.cat.g(self.t.S), self.cat.g(self.t.D)
        self.fails: dict[str, int] = {}
        self._els: dict[str, dict[str, Element]] = {}
        self._fam: dict[str, tuple[float, float, int, str]] | None = None
        self._calls: dict[str, dict[str, str]] = {}
        self._read: set[str] = set()
        self._fa: dict[str, Fasta | None] = {}
        self.no_fasta: list[str] = []

    def _fail(self, label: str) -> None:
        """Count one candidate failing the check `label`; None (no example from it)."""
        self.fails[label] = self.fails.get(label, 0) + 1

    def elements(self, X: str) -> dict[str, Element]:
        if X not in self._els:
            path = self.run.outdir / "work" / "elements" / f"{X}.tsv"
            self._els[X] = {e.eid: e for e in read_elements(path)} if path.is_file() else {}
        return self._els[X]

    def family(self, e: Element) -> tuple[float, float, int, str]:
        """(q25, q75, n, what) of the element lengths of e's family in R's annotation (its superfamily's
        when the family has < FAMILY_MIN elements)."""
        if self._fam is None:
            fam: dict[str, list[int]] = {}
            sup: dict[str, list[int]] = {}
            for x in self.elements(self.t.R).values():
                fam.setdefault(x.family or ".", []).append(x.end - x.start + 1)
                sup.setdefault(x.superfamily or ".", []).append(x.end - x.start + 1)
            self._fam = {}
            for name, lens in fam.items():
                self._fam["f:" + name] = (*np.percentile(lens, [25, 75]).tolist(), len(lens), "family")
            for name, lens in sup.items():
                self._fam["s:" + name] = (*np.percentile(lens, [25, 75]).tolist(), len(lens), "superfamily")
        f = self._fam.get("f:" + (e.family or "."))
        if f is None or f[2] < FAMILY_MIN:
            f = self._fam.get("s:" + (e.superfamily or "."), (math.nan, math.nan, 0, "superfamily"))
        return f

    def prefetch(self, js: list[int]) -> None:
        """Read R's own calls into S for the elements of painting rows `js` in one pass over the pair file."""
        p, cat = self.m.target, self.cat
        want = {cat.eid[int(p.idx[j]), self.gR] for j in js} - self._read
        want.discard(None)
        if not want:
            return
        path = self.run.outdir / "pairs" / f"{self.t.R}{PAIR_SEP}{self.t.S}.calls.tsv.gz"
        for _, f in stream_calls(path):
            if f[0] in want:
                self._calls[f[0]] = dict(zip(CALL_COLUMNS, f, strict=True))
        self._read |= want

    def fasta(self, X: str) -> Fasta | None:
        if X not in self._fa:
            s = next(x for x in self.cat.samples if x.label == X)
            try:
                self._fa[X] = Fasta(s.fasta, cache_dir=self.run.outdir / "work" / "fai")
            except (OSError, RuntimeError, ValueError) as exc:  # a missing or unreadable FASTA: no base-pair panel
                log.warning("overview: cannot read %s's FASTA %s (%s): no base-pair example", X, s.fasta, exc)
                self._fa[X] = None
                self.no_fasta.append(X)
        return self._fa[X]

    def clean(self, js: np.ndarray, med_k: float) -> list[dict]:
        """The clean candidates among painting rows `js` (R+D loci), nearest the median K first (position when
        no K is known)."""
        cat, t, p = self.cat, self.t, self.m.target
        out = []
        for j in js.tolist():
            i = int(p.idx[j])
            if not (cat.ptier[t.R][i, self.gD] == TIER_CODE["high"] or cat.member[i, self.gD]):
                self._fail("donor's call not high-confidence")
                continue
            e = self.elements(t.R).get(cat.eid[i, self.gR])
            if e is None:
                self._fail("element missing from the annotation")
                continue
            if not e.tsd:
                self._fail("no annotated TSD")
                continue
            if cat.proj[t.R][i, self.gS] != STATE_CODE["E"]:
                self._fail("sister's site not an empty site")
                continue
            c = self._calls.get(e.eid)
            if c is None or "tsd_match" not in (c.get("flags") or "").split(","):
                self._fail("sister's empty site without a matching TSD")
                continue
            ex = self.verify(j, e, c)
            if ex is not None:
                out.append(ex)
        key = (lambda e: (abs(e["k"] - med_k) if math.isfinite(e["k"]) and math.isfinite(med_k) else math.inf,
                          e["start"]))
        return sorted(out, key=key)

    def verify(self, j: int, e: Element, call: dict[str, str]) -> dict | None:
        """The FASTA checks of one candidate (module docstring) against R's own annotated TSD: its drawing data,
        or None (the failed check counted)."""
        cat, t, p = self.cat, self.t, self.m.target
        i = int(p.idx[j])
        tsd = e.tsd.upper()
        L = len(tsd)
        fr, fs, fd = self.fasta(t.R), self.fasta(t.S), self.fasta(t.D)
        if fr is None or fs is None or fd is None:
            return self._fail("a genome sequence unreadable")
        r = _ends(fr, e.seqid, e.start, e.end, L, 1)
        if r is None:
            return self._fail("sequence ends")
        if r["ltsd"] != r["rtsd"]:
            return self._fail("recipient's TSD copies differ")
        if r["ltsd"] != tsd:
            return self._fail("recipient's TSD copies differ from the annotated TSD")
        if not (r["body5"].startswith("TG") and r["body3"].endswith("CA")):
            return self._fail("recipient's ends not TG...CA")
        l1, l2 = e.ltr_left[1] - e.ltr_left[0] + 1, e.ltr_right[1] - e.ltr_right[0] + 1
        if min(l1, l2) <= 0 or min(l1, l2) < LTR_RATIO * max(l1, l2):
            return self._fail("an LTR truncated")
        q25, q75, nfam, what = self.family(e)
        length = e.end - e.start + 1
        if not (q25 <= length <= q75):
            return self._fail("length atypical of its family")
        d = self.donor(i)
        if d is None:
            return self._fail("donor's copy without a site")
        dfa = _ends(fd, d["seqid"], d["start"], d["end"], L, d["o"])
        if dfa is None or not (dfa["ltsd"] == dfa["rtsd"] == tsd):
            return self._fail("donor's TSD copies not identical to the recipient's")
        if not (dfa["body5"].startswith("TG") and dfa["body3"].endswith("CA")):
            return self._fail("donor's ends not TG...CA")
        s = _site(fs, call, L)
        if s is None or s["tsd"] != tsd:
            return self._fail("sister's TSD copy differs")
        il = _ident(r["lflank"][-FLANK_CHECK:], s["lflank"][-FLANK_CHECK:])
        ir = _ident(r["rflank"][:FLANK_CHECK], s["rflank"][:FLANK_CHECK])
        if min(il, ir) < FLANK_IDENT:
            return self._fail("sister's flanks differ")
        return {"j": j, "i": i, "lid": str(cat.lid[i]), "el": e, "tsd": tsd, "k": float(cat.k[i, self.gR]),
                "start": e.start, "end": e.end, "seqid": e.seqid, "length": length, "ltr": (l1, l2),
                "fam": (q25, q75, nfam, what), "r": r, "d": dfa, "dsite": d, "s": s, "s_ident": (il, ir),
                "family": e.family or e.superfamily or "unknown"}

    def donor(self, i: int) -> dict | None:
        """D's copy: its annotated element when D annotated the locus, else R's call into D (its site)."""
        cat, t = self.cat, self.t
        o = int(cat.pbstrand[t.R][i, self.gD]) or 1
        if cat.member[i, self.gD]:
            de = self.elements(t.D).get(cat.eid[i, self.gD])
            if de is not None:
                return {"seqid": de.seqid, "start": de.start, "end": de.end, "o": o, "annotated": True,
                        "k": de.k}
        sq = int(cat.pseq[t.R][i, self.gD])
        if sq < 0:
            return None
        return {"seqid": cat.seqids[t.D][sq], "start": int(cat.pstart[t.R][i, self.gD]),
                "end": int(cat.pend[t.R][i, self.gD]), "o": o, "annotated": False, "k": None}

    def why_none(self, any_tract: bool) -> str:
        """Why no base-pair example is drawn, in plain words."""
        if self.no_fasta:
            return f"no base-pair example: {', '.join(sorted(set(self.no_fasta)))}'s genome sequence cannot be read"
        if not self.fails:
            return "no base-pair example: no R+D insertion" + (" inside a tract" if any_tract else "")
        top = max(self.fails.items(), key=lambda kv: (kv[1], kv[0]))
        return f"no base-pair example: no R+D insertion passed the structural checks (most often: {top[0]})"


def _ends(fa: Fasta, seqid: str, start: int, end: int, L: int, o: int) -> dict | None:
    """An element [start, end] (1-based) with TSD length L, in R's orientation (o = -1: reverse-complemented):
    flanks (max(FLANK_CHECK, FLANK_SHOW) bp), the TSD copies, the outer LTR_SHOW bp of each end and the body's
    first and last two bases (body5, body3); None where the sequence ends too soon."""
    F = max(FLANK_CHECK, FLANK_SHOW)
    if start - 1 - L - F < 0:
        return None
    left = fa.fetch(seqid, start - 1 - L - F, start - 1)
    right = fa.fetch(seqid, end, end + L + F)
    head = fa.fetch(seqid, start - 1, start - 1 + LTR_SHOW)
    tail = fa.fetch(seqid, end - LTR_SHOW, end)
    if len(left) != L + F or len(right) != L + F or len(head) != LTR_SHOW or len(tail) != LTR_SHOW:
        return None
    if o < 0:
        left, right, head, tail = revcomp(right), revcomp(left), revcomp(tail), revcomp(head)
    return {"lflank": left[:F], "ltsd": left[F:], "head": head, "tail": tail, "rtsd": right[:L],
            "rflank": right[L:], "body5": head[:2], "body3": tail[-2:], "inner": end - start + 1 - 2 * LTR_SHOW}


def _site(fa: Fasta, call: dict[str, str], L: int) -> dict | None:
    """S's empty site from R's E call (its single TSD copy b_start..b_end), in R's orientation."""
    try:
        a, b = int(call["b_start"]), int(call["b_end"])
    except (KeyError, ValueError):
        return None
    if b - a + 1 != L:
        return None
    F = max(FLANK_CHECK, FLANK_SHOW)
    sid = call["b_seqid"]
    if a - 1 - F < 0:
        return None
    left, tsd, right = fa.fetch(sid, a - 1 - F, a - 1), fa.fetch(sid, a - 1, b), fa.fetch(sid, b, b + F)
    if len(left) != F or len(right) != F:
        return None
    if call.get("b_strand") == "-":
        left, tsd, right = revcomp(right), revcomp(tsd), revcomp(left)
    return {"lflank": left, "tsd": tsd, "rflank": right, "seqid": sid, "start": a, "end": b,
            "strand": call.get("b_strand")}


def _ident(a: str, b: str) -> float:
    if not a or len(a) != len(b):
        return 0.0
    return sum(x == y and x in "ACGT" for x, y in zip(a, b, strict=True)) / len(a)


def _absence_table(run: _Run) -> Iterator[dict[str, str]]:
    """mechanisms.tsv's absence rows over every superfamily and family ("all")."""
    for r in run.table("mechanisms.tsv") or []:
        if r.get("table") == "absence" and r.get("superfamily") == "all" and r.get("family") == "all":
            yield r


def _absence_rows(run: _Run, t: _Trio) -> dict[str, dict[str, float]]:
    """The absence rows of R's insertions in S: age bin -> metric -> value."""
    out: dict[str, dict[str, float]] = {}
    for r in _absence_table(run):
        if r.get("scope") == f"{t.R}>{t.S}":
            out.setdefault(r["age_bin"], {})[r["metric"]] = _float(r["value"])
    return out


def _absence_others(run: _Run, t: _Trio) -> dict[str, dict[str, float]]:
    """The same, all ages, for every other ordered pair of ingroup genomes: scope -> metric -> value."""
    ing = set(run.cat.ingroup)
    out: dict[str, dict[str, float]] = {}
    for r in _absence_table(run):
        sc = r.get("scope", "")
        a, _, b = sc.partition(">")
        if a in ing and b in ing and sc != f"{t.R}>{t.S}" and r.get("age_bin") == "all":
            out.setdefault(sc, {})[r["metric"]] = _float(r["value"])
    return out


def _absence_value(row: dict[str, float], metric: str) -> float:
    return sum(row.get(x, math.nan) for x in metric.split("+"))


# ---------------------------------------------------------------- the no-trio model


def _no_trio_model(m: _Model, nums: style.Numbers) -> None:
    """Fewer than three ingroup genomes: the first annotated ingroup genome X's loci by their state in the other
    ingroup genome Y (panel A) and the state composition of every ordered pair (panel B)."""
    cat = m.run.cat
    nums.add("no_trio", "fewer than three ingroup genomes", "stats.tsv:catalog/all/all/n_loci")
    ann = [x for x in cat.ingroup if cat.annotated[x]]
    pairs = []
    for a in [x for x in cat.labels if cat.annotated[x]]:
        for b in cat.labels:
            if b != a:
                pairs.append((a, b, _composition(cat, a, b)))
    m.no_trio = {"pairs": pairs, "X": None, "Y": None}
    if ann and len(cat.ingroup) == 2:
        X = ann[0]
        Y = next(y for y in cat.ingroup if y != X)
        drawn, m.n_not = _drawn(cat, X)
        m.n_drawn = len(drawn)
        m.xmax = max([int(cat.lengths[X][s]) for s in drawn] or [1])
        m.bin_bp = _nice_floor(m.xmax / BINS)
        gx = cat.g(X)
        rows = np.flatnonzero(cat.member[:, gx] & (cat.seq[:, gx] >= 0))
        st = cat.proj[X][rows, cat.g(Y)]
        cls = np.where(np.isin(st, _PRESENT), RS, np.where(st == STATE_CODE["E"], RD, MISSING)).astype(np.int8)
        p = _Painting(X, cat.lid[rows], rows, cat.seq[rows, gx].astype(np.int64), cat.start[rows, gx],
                      cat.end[rows, gx], cls, np.full(len(rows), np.nan), [], None, None)
        labels = _short_names([cat.seqids[X][s] for s in drawn])
        for s, lab in zip(drawn, labels, strict=True):
            m.rows.append((_row(p, cat, s, f"{X} {lab}", m.bin_bp, {}), None))
        m.no_trio.update(X=X, Y=Y, n_present=int((cls == RS).sum()), n_empty=int((cls == RD).sum()),
                         n=len(cls))
        m.text["unit"] = "chromosome" if _chromosome_scale(cat, {X: drawn}) else "sequence"
        nums.add("sequences_drawn", m.n_drawn, f"loci.tsv:{X}_seqid")
        nums.add("sequences_not_drawn", m.n_not, f"loci.tsv:{X}_seqid")
        m.panels = ["A", "B"]
    else:
        nums.add("sequences_drawn", 0, "loci.tsv:lid")
        nums.add("sequences_not_drawn", sum(len(v) for v in cat.lengths.values()), "loci.tsv:lid")
        m.panels = ["B"]


def _composition(cat: Catalog, a: str, b: str) -> dict[str, int]:
    """Counts of a's own calls into b at a's member loci, by state."""
    st = cat.proj[a][cat.member[:, cat.g(a)], cat.g(b)]
    inv = {v: k for k, v in STATE_CODE.items()}
    vals, cnt = np.unique(st, return_counts=True)
    return {inv.get(int(v), "none"): int(c) for v, c in zip(vals.tolist(), cnt.tolist(), strict=True)}


# ---------------------------------------------------------------- formatting


def _int(n) -> str:
    return f"{int(n):,}"


def _pct(x: float) -> str:
    """A finite share as a percentage, its precision by its size (callers pass finite values: `_fmt`)."""
    if x == 0:
        return "0%"
    v = 100 * x
    return f"{v:.0f}%" if v >= 9.5 else f"{v:.1f}%" if v >= 0.95 else f"{v:.2f}%"


def _bp(x: float) -> str:
    if x >= 1e7:
        return f"{x / 1e6:.0f} Mb"
    if x >= 1e6:
        return f"{x / 1e6:.2f} Mb".replace(".00 Mb", " Mb")
    if x >= 1e5:
        return f"{x / 1e3:.0f} kb"
    if x >= 1e3:
        return f"{x / 1e3:.1f} kb".replace(".0 kb", " kb")
    return f"{x:.0f} bp"


def _role(m: _Model, x: str) -> str:
    g = m.groups.get(x)
    return f"{x} ({g})" if g and g != x else x


# ---------------------------------------------------------------- the page


class _Page:
    """A figure laid out in mm from its top-left corner."""

    def __init__(self, height_mm: float):
        self.W, self.H = WIDTH_MM, height_mm
        self.fig = style.page(height_mm)

    def ax(self, x: float, y: float, w: float, h: float, **kw):
        return style.axes_mm(self.fig, x, y, w, h, **kw)

    def letter(self, x: float, y: float, s: str) -> None:
        style.panel_letter(self.fig, s, x_mm=x, y_mm=y)

    def text(self, x: float, y: float, s: str, **kw) -> None:
        kw.setdefault("fontsize", FS)
        self.fig.text(x / self.W, 1 - y / self.H, s, **kw)


_bare = style.bare
_text_mm = style.text_mm


def _spines(ax, keep=("left", "bottom")) -> None:
    for name, s in ax.spines.items():
        s.set_visible(name in keep)


def _bp_axis(ax, lo: float, hi: float, label: str | None = None) -> None:
    """x ticks in Mb, kb or bp by the span."""
    from matplotlib.ticker import FuncFormatter, MaxNLocator

    span = hi - lo
    unit, name = (1e6, "Mb") if span >= 2e6 or hi >= 5e6 else (1e3, "kb") if hi >= 5e3 else (1, "bp")
    ax.xaxis.set_major_locator(MaxNLocator(nbins=6, steps=[1, 2, 2.5, 5, 10]))
    digits = 0 if span / unit >= 5 else 1 if span / unit >= 0.5 else 2
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v / unit:,.{digits}f}"))
    if label is not None:
        ax.set_xlabel(f"{label} ({name})")


def _draw(m: _Model, mono: str):
    """The figure (page 1 of the PDF and the PNG)."""
    if m.trio is None:
        return _draw_no_trio(m) if not m.screen else _draw_text(m, "no trio: fewer than three ingroup genomes")
    n_pairs = len(m.rows)
    ncol = 2 if n_pairs > 3 else 1
    per = -(-n_pairs // ncol) if n_pairs else 1
    pair_h = min(14.0, 68.0 / per)
    lay = style.key_layout(_key_items(m), WIDTH_MM - 6.0 - 1)
    key_h = (max(ln for *_, ln in lay) + 1) * KEY_LINE_MM + 1.2
    a_h = 14.0 if "A" in m.placeholders else key_h + per * pair_h + 6.0
    sec = [("A", a_h)]
    if "B" in m.panels:
        sec.append(("B", 40.0 + 1.6 * max(0, _n_sub(m) - 2) * (3 + len(m.run.cat.outgroups))))
    if m.screen:
        sec.append(("screen", 42.0 if "G" in m.panels else 14.0 * -(-len(m.placeholders & set(m.panels)) // 3)))
    if ("C" in m.panels or "D" in m.panels) and not m.screen:
        sec.append(("CD", 30.0))
    if not m.screen:
        sec.append(("EFG", 42.0))
    if "H" in m.panels and not m.screen:
        sec.append(("H", _h_layout(m)["height"]))
    gap = 5.0
    H = sum(h for _, h in sec) + gap * (len(sec) - 1) + 3.0
    page = _Page(H)
    y = 1.0
    for name, h in sec:
        if name == "A" and "A" in m.placeholders:
            _placeholder(page, m, "A", 0, y, WIDTH_MM, h)
        elif name == "A":
            _panel_a(page, m, y, ncol, per, pair_h, key_h)
        elif name == "B":
            _panel_b(page, m, y, h)
        elif name == "screen":
            _screen_row(page, m, y, h)
        elif name == "CD":
            has_c, has_d = "C" in m.panels, "D" in m.panels
            if has_c and has_d:
                _panel_c(page, m, 0, y, 108.0, h, mono)
                _panel_d(page, m, 109.0, y, h)
            elif has_c:
                _panel_c(page, m, 0, y, WIDTH_MM, h, mono)
            else:
                _panel_d(page, m, 0, y, h)
        elif name == "EFG":
            _panels_efg(page, m, y, h)
        else:
            _panel_h(page, m, y)
        y += h + gap
    return page.fig


# ---------------------------------------------------------------- panel A


def _panel_a(page: _Page, m: _Model, y0: float, ncol: int, per: int, pair_h: float, key_h: float
             ) -> None:
    page.letter(0.5, y0, m.letters["A"])
    _key_a(page, m, 6.0, y0 + 0.2)
    label_w = 15.0
    col_gap = 4.0
    col_w = (RIGHT - 1.0 - ncol * label_w - (ncol - 1) * col_gap) / ncol
    top = y0 + key_h
    row_h = (pair_h - 2.2) / 2
    for k, (r, d) in enumerate(m.rows):
        c, j = divmod(k, per)
        x = 1.0 + label_w + c * (col_w + label_w + col_gap)
        yy = top + j * pair_h
        ax_r = page.ax(x, yy, col_w, row_h)
        _draw_row(ax_r, r, m, m.painted)
        if d is not None:
            ax_d = page.ax(x, yy + row_h + 0.4, col_w, row_h)
            _draw_row(ax_d, d, m, m.painted)
        else:
            ax_d = None
        win = m.window
        if win is not None and r.seq == win["seq"]:
            _zoom_box(ax_r, win, m, col_w, top=j == 0)
        last = j == per - 1 or k == len(m.rows) - 1
        target = ax_d if ax_d is not None else ax_r
        if last:
            target.spines["bottom"].set_visible(True)
            target.spines["bottom"].set_position(("outward", 1.0))
            _bp_axis(target, 0, m.xmax, "position")
            target.tick_params(axis="x", length=2.0)
            target.xaxis.label.set_size(FS)


def _key_items(m: _Model) -> list[tuple[str, str | None, str]]:
    t = m.trio
    roles = [("role", "R", f"recipient {_role(m, t.R)}"),
             ("role", "S", f"sister {_role(m, t.S)}, the decoy"),
             ("role", "D", f"donor {_role(m, t.D)}")]
    if m.run.cat.outgroups:
        roles.append(("role", "O", "outgroup " + ", ".join(_role(m, o) for o in m.run.cat.outgroups)))
    ticks = ((("tick", ACC, "near-identical copy in D only"),
              ("tick", BLUE, "near-identical copy in the other sister only"))
             if m.screen else (("tick", ACC, "insertion shared with D, the other sister's site empty"),
                               ("tick", BLUE, "shared with the other sister, D's site empty")))
    items = roles + [*ticks,
                     ("cmap", None, "posterior of D-like ancestry"), ("tract", None, "tract"),
                     ("strong", None, f"tract, q \u2264 {Q_TRACT:g} ({_n_text(m)})"), ("band", BAND,
                                                                                    "inverted block"),
                     ("hatch", GRAY, "no informative element" if m.screen else "no informative insertion")]
    if not m.painted:
        items = [x for x in items if x[0] not in ("cmap", "tract", "strong")]
    return items


def _key_a(page: _Page, m: _Model, x0: float, y0: float) -> float:
    """Panel A's key (the roles, then the symbols, the posterior's color scale among them); its height (mm)."""
    width = WIDTH_MM - x0 - 1
    items = _key_items(m)
    n = max(ln for *_, ln in style.key_layout(items, width)) + 1
    h = n * KEY_LINE_MM + 0.4
    style.key_row(page.ax(x0, y0, width, h), items, cmap=_post_cmap())
    return h


def _n_text(m: _Model) -> str:
    """The n of the q <= Q_TRACT tract class for its key entry: "none" or "n = k"."""
    v = _num(_paint_row(m, "n_tracts_q05"))
    return "n unknown" if not math.isfinite(v) else "none" if v == 0 else f"n = {int(v)}"


_CMAP = None


def _post_cmap():
    """The posterior's color scale: one hue, light (0) to dark (1), colorblind-safe (ColorBrewer Oranges)."""
    global _CMAP
    if _CMAP is None:
        from matplotlib.colors import LinearSegmentedColormap

        _CMAP = LinearSegmentedColormap.from_list(
            "ltrtrace_post", ["#FFF5EB", "#FDD0A2", "#FD8D3C", "#D94801", "#7F2704"])
    return _CMAP


def _draw_row(ax, r: _Row, m: _Model, painted: bool) -> None:
    """One painted sequence (panel A), y 0-10: RS ticks 0.3-3.0, the bar 3.5-6.1 (its color the mean posterior
    per display bin, hatched where no informative insertion), tracts as a dark line on its top edge (black with
    a mark above at q <= Q_TRACT), RD ticks 6.9-9.7; inverted blocks as gray bands behind."""
    from matplotlib.collections import LineCollection, PatchCollection
    from matplotlib.patches import Rectangle

    ax.set_xlim(0, m.xmax)
    ax.set_ylim(0, 10)
    _bare(ax)
    for a, b in r.bands:
        ax.add_patch(Rectangle((a, 0), max(b - a, 1), 10, fc=BAND, ec="none", zorder=0))
    lo, hi = 3.5, 6.1
    if painted and np.isfinite(r.post).any():
        cm = _post_cmap()
        ok = np.isfinite(r.post)
        cells = [Rectangle((r.edges[k], lo), r.edges[k + 1] - r.edges[k], hi - lo) for k in np.flatnonzero(ok)]
        ax.add_collection(PatchCollection(cells, fc=cm(r.post[ok]), ec="none", lw=0, zorder=1,
                                          antialiased=False))
        if (~ok).any():
            ax.add_collection(PatchCollection([Rectangle((r.edges[k], lo), r.edges[k + 1] - r.edges[k], hi - lo)
                                               for k in np.flatnonzero(~ok)], fc=LIGHT, ec="none", zorder=1))
    else:
        ax.add_patch(Rectangle((0, lo), r.length, hi - lo, fc=LIGHT, ec="none", zorder=1))
    if r.gaps:
        ax.add_collection(PatchCollection([Rectangle((a, lo), b - a, hi - lo) for a, b in r.gaps], fc="white",
                                          ec=GRAY, lw=0.0, hatch="//////", zorder=2))
    ax.add_patch(Rectangle((0, lo), r.length, hi - lo, fill=False, ec=GRAY, lw=0.25, zorder=4))
    for a, b, strong in r.tracts:
        if strong:
            continue
        ax.plot([a - 1, b], [hi + 0.25, hi + 0.25], color=DARK, lw=0.7, solid_capstyle="butt", zorder=5)
    for a, b, strong in r.tracts:
        if strong:
            ax.plot([a - 1, b], [hi + 0.3, hi + 0.3], color="black", lw=1.4, solid_capstyle="butt", zorder=6)
            ax.scatter([(a + b) / 2], [9.9], marker="v", s=6, color="black", lw=0, zorder=6, clip_on=False)
    for xs, y0, y1, color in ((r.rd, 6.9, 9.7, ACC), (r.rs, 0.3, 3.0, BLUE)):
        if len(xs):
            seg = np.empty((len(xs), 2, 2))
            seg[:, :, 0] = xs[:, None]
            seg[:, 0, 1], seg[:, 1, 1] = y0, y1
            ax.add_collection(LineCollection(seg, colors=color, linewidths=0.3, zorder=3))
    ax.text(-0.006 * m.xmax, (lo + hi) / 2, r.label, ha="right", va="center", fontsize=FS_SMALL, clip_on=False)


def _zoom_box(ax, win: dict, m: _Model, col_w: float, top: bool = False) -> None:
    """The window of panel B as a box on its recipient row, its letter clear of the data: above the box in the
    gap between rows, or (on a column's first row, under the key) right of it on a white backing."""
    from matplotlib.patches import Rectangle

    a, b = win["start"], win["end"]
    minw = m.xmax * 1.4 / col_w
    if b - a < minw:
        c = (a + b) / 2
        a, b = c - minw / 2, c + minw / 2
    ax.add_patch(Rectangle((a, 0.1), b - a, 9.8, fill=False, ec=DARK, lw=0.6, zorder=6, clip_on=False))
    xy, va = ((b + m.xmax * 0.004, 9.6), "top") if top else ((b, 10.15), "bottom")
    ax.text(*xy, m.letters.get("B", "B"), fontsize=FS, fontweight="bold", va=va, ha="left", zorder=7,
            clip_on=False, bbox={"fc": "white", "ec": "none", "pad": 0.2})


# ---------------------------------------------------------------- panel B


def _window_loci(m: _Model) -> dict:
    """Panel B's loci (R's painted loci starting inside the window) and their packing into rows."""
    w, p = m.window, m.target
    on = np.flatnonzero((p.seq == w["seq"]) & (p.end >= w["start"]) & (p.start <= w["end"]))
    s, e = p.start[on], p.end[on]
    gap = max(1, int((w["end"] - w["start"]) * 0.003))
    order = np.lexsort((-e, s))
    sub = np.zeros(len(on), dtype=int)
    ends: list[int] = []
    for j in order.tolist():
        for k, last in enumerate(ends):
            if s[j] > last + gap:
                sub[j], ends[k] = k, int(e[j])
                break
        else:
            sub[j] = len(ends)
            ends.append(int(e[j]))
    return {"rows": on, "sub": sub, "n_sub": max(1, len(ends))}


def _n_sub(m: _Model) -> int:
    if m.window is None:
        return 1
    if "loci" not in m.window:
        m.window["loci"] = _window_loci(m)
    return m.window["loci"]["n_sub"]


def _panel_b(page: _Page, m: _Model, y0: float, h: float) -> None:
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Rectangle

    page.letter(0.5, y0, m.letters["B"])
    cat, t, w, p = m.run.cat, m.trio, m.window, m.target
    loci = w.setdefault("loci", _window_loci(m))
    rows, sub, n_sub = loci["rows"], loci["sub"], loci["n_sub"]
    side_w, side_gap = 19.0, 4.0           # the right column: the control's histogram, the state key
    x0 = 16.0
    wd = RIGHT - x0 - side_w - side_gap
    xs_side = x0 + wd + side_gap
    genomes = [t.R, t.S, t.D, *cat.outgroups]
    post_h = 11.0 if m.painted else 0.0
    mat_top = y0 + 1.0 + post_h + (2.0 if m.painted else 0.0)
    mat_h = h - (mat_top - y0) - 7.0
    if m.painted:
        ax = page.ax(x0, y0 + 1.0, wd, post_h)
        ax.set_xlim(w["start"], w["end"])
        ax.set_ylim(0, 1.25)
        _spines(ax, ("left",))
        ax.set_yticks([0, 0.5, 1])
        ax.set_yticklabels(["0", "0.5", "1"])
        ax.set_xticks([])
        o = np.argsort(p.start[rows], kind="stable")
        xs, ys = p.start[rows][o], p.post[rows][o]
        ax.fill_between(xs, 0, ys, color=ACC, alpha=0.18, lw=0)
        ax.plot(xs, ys, color=ACC, lw=0.6)
        ax.scatter(xs, ys, s=1.5, color=ACC, lw=0, zorder=3)
        ax.text(-0.035, 0.5, "posterior", transform=ax.transAxes, ha="right", va="center", fontsize=FS_SMALL)
        tr = w.get("tract")
        if tr is not None:
            strong = math.isfinite(tr["q"]) and tr["q"] <= Q_TRACT
            ax.add_patch(Rectangle((tr["start"], 1.1), tr["end"] - tr["start"] + 1, 0.13, fc=ACC,
                                   ec=DARK if strong else "none", lw=0.5, alpha=1.0 if strong else 0.6,
                                   clip_on=False))
            ax.text(tr["end"] + (w["end"] - w["start"]) * 0.006, 1.165, "tract", fontsize=FS_SMALL, va="center",
                    ha="left")
        _posterior_histograms(page, m, xs_side, y0 + 1.0, side_w, post_h)
    ax = page.ax(x0, mat_top, wd, mat_h)
    ax.set_xlim(w["start"], w["end"])
    unit = 1.0
    block = n_sub * unit + 0.9
    total = block * len(genomes)
    ax.set_ylim(total, -0.2)
    _spines(ax, ("bottom",))
    ax.set_yticks([])
    label = next((r.label for r, _ in m.rows if r.seq == w["seq"]), f"{t.R} {w['name']}")
    _bp_axis(ax, w["start"], w["end"], label)
    color = {RD: ACC, RS: BLUE}
    proj = cat.proj[t.R]
    pres, empty, miss = [], [], []
    for gi, g in enumerate(genomes):
        gg = cat.g(g)
        ybase = gi * block
        ax.text(-0.012, ybase + n_sub * unit / 2, g, transform=ax.get_yaxis_transform(), ha="right",
                va="center", fontsize=FS)
        for j, rw in enumerate(rows.tolist()):
            i = int(p.idx[rw])
            c = color.get(int(p.cls[rw]), GRAY)
            st = STATE_CODE["P"] if g == t.R else int(proj[i, gg])
            rect = (max(p.start[rw], w["start"]) - 0.5, ybase + sub[j] * unit + 0.08,
                    max(min(p.end[rw], w["end"]) - max(p.start[rw], w["start"]) + 1, 1), unit * 0.78)
            if st in _PRESENT:
                pres.append((rect, c))
            elif st == STATE_CODE["E"] or st == _INFERRED:
                empty.append((rect, c))
            else:
                miss.append(rect)
    if pres:
        ax.add_collection(PatchCollection([Rectangle(r[:2], *r[2:]) for r, _ in pres],
                                          fc=[c for _, c in pres], ec=[c for _, c in pres], lw=0.3, zorder=3))
    if empty:
        ax.add_collection(PatchCollection([Rectangle(r[:2], *r[2:]) for r, _ in empty], fc="white",
                                          ec=[c for _, c in empty], lw=0.6, zorder=3))
    if miss:
        ax.add_collection(PatchCollection([Rectangle(r[:2], *r[2:]) for r in miss], fc="white", ec=GRAY,
                                          lw=0.0, hatch="//////", zorder=2))
    _side_key(page, [("fill", "present"), ("open", "absent, inferred" if m.screen else "empty site"),
                     ("hatch", "missing")], xs_side, mat_top + (6.0 if m.painted else 1.5))
    ex = m.example
    if ex is not None:
        a, b = ex["start"], ex["end"]
        span = w["end"] - w["start"]
        pad = span * 0.004
        ax.add_patch(Rectangle((a - pad, -0.15), b - a + 2 * pad, total + 0.1, fill=False, ec=DARK, lw=0.6,
                               zorder=6))
        ax.text(b + pad * 2, -0.1, m.letters.get("C", "C"), fontsize=FS, fontweight="bold", va="top",
                ha="left", zorder=6)


def _posterior_histograms(page: _Page, m: _Model, x: float, y: float, w: float, h: float) -> None:
    """Beside B's posterior track, on its posterior axis: the posterior of every painted insertion genome-wide in
    bins of 0.1, the recipient (vermillion) against the sister's painting with the roles swapped (gray): how
    often the mirror control, in which lineage sorting acts as in R, reaches a given posterior. Not a per-locus
    test (the tract's q is the test)."""
    ax = page.ax(x, y, w, h)
    ax.set_ylim(0, 1.25)
    edges = np.linspace(0, 1, 11)
    top = 0.0
    for paint, color in ((m.decoy, GRAY), (m.target, ACC)):
        if paint is None or not np.isfinite(paint.post).any():
            continue
        v = paint.post[np.isfinite(paint.post)]
        share = np.histogram(np.clip(v, 0, 1), edges)[0] / len(v)
        top = max(top, float(share.max()))
        ax.stairs(share, edges, orientation="horizontal", color=color, lw=0.9, baseline=None)
    ax.set_xlim(0, max(0.05, top * 1.15))
    _spines(ax, ("left", "bottom"))
    ax.set_yticks([0, 0.5, 1])
    ax.set_yticklabels([])
    from matplotlib.ticker import FuncFormatter, MaxNLocator

    ax.xaxis.set_major_locator(MaxNLocator(nbins=2, steps=[1, 2, 5, 10]))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{100 * v:.0f}"))
    ax.set_xlabel("% of elements" if m.screen else "% of insertions", fontsize=FS_SMALL, labelpad=1.0)
    on_y = ax.get_yaxis_transform()  # x in axes fraction, y on the posterior axis: the empty upper right
    ax.text(0.98, 0.95, m.trio.R, color=ACC, ha="right", va="center", fontsize=FS_SMALL, transform=on_y)
    if m.decoy is not None:
        ax.text(0.98, 0.72, "decoy", color=GRAY, ha="right", va="center", fontsize=FS_SMALL, transform=on_y)


def _side_key(page: _Page, items: list[tuple[str, str]], x: float, y: float) -> None:
    """A vertical state key from (x, y) mm: present filled, empty site open, missing hatched (gray: the colors
    are the class's)."""
    from matplotlib.patches import Rectangle

    for k, (kind, text) in enumerate(items):
        yy = y + k * 3.0
        kw = {"fill": {"fc": GRAY, "ec": GRAY, "lw": 0.3}, "open": {"fc": "white", "ec": GRAY, "lw": 0.6},
              "hatch": {"fc": "white", "ec": GRAY, "lw": 0.0, "hatch": "//////"}}[kind]
        page.fig.add_artist(Rectangle((x / page.W, 1 - (yy + 0.6) / page.H), 2.4 / page.W, 1.2 / page.H,
                                      transform=page.fig.transFigure, **kw))
        page.text(x + 3.2, yy, text, va="center", ha="left", fontsize=FS_SMALL)


# ---------------------------------------------------------------- panel C


def _panel_c(page: _Page, m: _Model, x0: float, y0: float, w: float, h: float, mono: str) -> None:
    from matplotlib.patches import Rectangle

    page.letter(x0 + 0.5, y0, m.letters["C"])
    ex, t = m.example, m.trio
    r, d, s = ex["r"], ex["d"], ex["s"]
    F, L, M = FLANK_SHOW, len(ex["tsd"]), LTR_SHOW
    G = 17
    N = 2 * F + 2 * L + 2 * M + G
    lab_w = 11.0
    ax = page.ax(x0 + lab_w, y0 + 1.0, w - lab_w - 2.0, h - 1.5)
    ax.set_xlim(-0.5, N - 0.5)
    ax.set_ylim(10.6, 0)
    _bare(ax)
    # the element to scale (cartoon)
    el = ex["el"]
    length = el.end - el.start + 1
    ea, eb = 0.2 * N, 0.8 * N
    sc = (eb - ea) / length
    yc = 1.25
    ax.plot([-0.5, N - 0.5], [yc, yc], color=DARK, lw=0.5)
    l1 = (el.ltr_left[1] - el.ltr_left[0] + 1) * sc
    l2 = (el.ltr_right[1] - el.ltr_right[0] + 1) * sc
    ax.add_patch(Rectangle((ea, yc - 0.45), eb - ea, 0.9, fc=INNER, ec="none"))
    ax.add_patch(Rectangle((ea, yc - 0.45), l1, 0.9, fc=LTRC, ec="none"))
    ax.add_patch(Rectangle((eb - l2, yc - 0.45), l2, 0.9, fc=LTRC, ec="none"))
    ax.text(ea + l1 / 2, yc - 0.7, "LTR", ha="center", va="bottom", fontsize=FS_SMALL)
    ax.text(eb - l2 / 2, yc - 0.7, "LTR", ha="center", va="bottom", fontsize=FS_SMALL)
    ax.text((ea + eb) / 2, yc - 0.7, f"{_bp(length)} element", ha="center", va="bottom", fontsize=FS_SMALL)
    kb = 1000 * sc if length >= 2000 else 100 * sc
    ax.plot([N - 1.0 - kb, N - 1.0], [yc + 0.55, yc + 0.55], color=DARK, lw=0.6, solid_capstyle="butt")
    ax.text(N - 1.0, yc + 0.7, "1 kb" if length >= 2000 else "100 bp", ha="right", va="top", fontsize=FS_SMALL)
    c_head, c_tail = F + L, F + L + 2 * M + G
    for cx, ex_x in ((c_head - 0.5, ea), (c_tail - 0.5, eb)):
        ax.plot([ex_x, cx], [yc + 0.45, 3.35], color=GRAY, lw=0.4)
    # the base-pair rows: R, D (insertion), S (empty site)
    rows = [(t.R, r, "R"), (t.D, d, "D"), (t.S, s, "S")]
    ys = (4.3, 6.55, 8.8)

    def columns(seq: dict, kind: str) -> list:
        """The row's bases by column (None: not drawn): flank, TSD, LTR heads, gap, LTR tails, TSD, flank; an
        empty site has its one TSD copy and no element."""
        if kind == "S":
            return list(seq["lflank"][-F:]) + list(seq["tsd"]) + [None] * (2 * M + G + L) + list(seq["rflank"][:F])
        return (list(seq["lflank"][-F:]) + list(seq["ltsd"]) + list(seq["head"]) + [None] * G
                + list(seq["tail"]) + list(seq["rtsd"]) + list(seq["rflank"][:F]))

    for (g, seq, kind), yy in zip(rows, ys, strict=True):
        page.text(x0 + lab_w - 1.0, y0 + 1.0 + (h - 1.5) * yy / 10.6, g, ha="right", va="center",
                  fontsize=FS)
        for c, base in enumerate(columns(seq, kind)):
            if base is None:
                continue
            ltr = kind != "S" and (F + L <= c < F + L + M or F + L + M + G <= c < F + L + 2 * M + G)
            if ltr:
                ax.add_patch(Rectangle((c - 0.5, yy - 0.62), 1.0, 1.24, fc=LTRC, ec="none"))
            ax.text(c, yy, base, ha="center", va="center", fontsize=FS, family=mono,
                    color="white" if ltr else DARK)
        tsds = [F] if kind == "S" else [F, F + L + 2 * M + G]
        for c in tsds:
            ax.add_patch(Rectangle((c - 0.5, yy - 0.7), L, 1.4, fill=False, ec=DARK, lw=0.6))
        if kind == "S":
            a, b = F + L, F + 2 * L + 2 * M + G
            ax.plot([a - 0.3, b - 0.7], [yy, yy], color=GRAY, lw=0.5, ls=(0, (2, 2)))
            ax.text((a + b) / 2, yy, " empty site ", ha="center", va="center", fontsize=FS_SMALL,
                    bbox={"fc": "white", "ec": "none", "pad": 0.4})
            ax.text(F + L / 2 - 0.5, yy + 0.95, "single TSD copy", ha="center", va="top", fontsize=FS_SMALL)
        else:
            inner = seq["inner"]
            ax.text(F + L + M + G / 2 - 0.5, yy, f"... {_int(inner)} bp ...", ha="center", va="center",
                    fontsize=FS_SMALL)
    ax.text(F + L / 2 - 0.5, ys[0] - 0.95, "TSD", ha="center", va="bottom", fontsize=FS_SMALL)
    ax.text(F + L + 2 * M + G + L / 2 - 0.5, ys[0] - 0.95, "TSD", ha="center", va="bottom", fontsize=FS_SMALL)
    # match ticks between rows
    rc, dc, sc_ = columns(r, "R"), columns(d, "D"), columns(s, "S")
    for c in range(N):
        if rc[c] is not None and rc[c] == dc[c]:
            ax.plot([c, c], [ys[0] + 0.7, ys[1] - 0.7], color=GRAY, lw=0.35)
        if dc[c] is not None and sc_[c] is not None and dc[c] == sc_[c]:
            ax.plot([c, c], [ys[1] + 0.7, ys[2] - 0.7], color=GRAY, lw=0.35)


# ---------------------------------------------------------------- panel D: absences


def _panel_d(page: _Page, m: _Model, x0: float, y0: float, h: float) -> None:
    page.letter(x0 + 0.5, y0, m.letters["D"])
    t = m.trio
    rows = _absence_rows(m.run, t)
    others = _absence_others(m.run, t)
    lab_w = 31.0
    ax = page.ax(x0 + lab_w, y0 + 4.5, RIGHT - x0 - lab_w, h - 13.5)
    n = len(_ABSENCE)
    ax.set_ylim(n - 0.4, -0.6)
    ax.set_xscale("log")
    values = [_absence_value(r, metric) for r in (*rows.values(), *others.values()) for metric, _ in _ABSENCE]
    pos = [v for v in values if math.isfinite(v) and v > 0]
    dec = 10.0 ** min(-4, math.floor(math.log10(min(pos)))) if pos else 1e-4
    zero = dec / 10 ** 0.6                  # a zero share has no place on a log axis: its own "0" tick
    has_zero = any(math.isfinite(v) and v <= 0 for v in values)
    ax.set_xlim(zero / 10 ** 0.3 if has_zero else dec / 10 ** 0.25, 1.4)
    _spines(ax, ("bottom",))
    ax.set_yticks(range(n))
    ax.set_yticklabels([lab for _, lab in _ABSENCE])
    ax.tick_params(axis="y", length=0)
    decades = [10.0 ** e for e in range(round(math.log10(dec)), 1)]
    ax.set_xticks([zero, *decades] if has_zero else decades)
    ax.set_xticklabels((["0"] if has_zero else []) + [f"{100 * d:g}%" for d in decades])
    ax.xaxis.set_minor_locator(_minor_log())
    ax.set_xlabel(f"share of absences\n({t.R} insertions in {t.S})")

    def at(v: float) -> float:
        return zero if v <= 0 else v

    for k in range(n):
        ax.axhline(k, color=LIGHT, lw=0.3, zorder=0)
    for k, (metric, _) in enumerate(_ABSENCE):
        xs = [at(v) for v in (_absence_value(o, metric) for o in others.values()) if math.isfinite(v)]
        if xs:
            ax.scatter(xs, [k + 0.28] * len(xs), marker="|", s=10, color=GRAY, lw=0.5, zorder=2, clip_on=False)
        vals = [_absence_value(rows.get(a, {}), metric) for a in AGE_LABELS]
        pts = [(at(v), c) for v, c in zip(vals, _AGE_GRAYS, strict=True) if math.isfinite(v)]
        if len(pts) > 1:
            ax.plot([v for v, _ in pts], [k] * len(pts), color=LIGHT, lw=0.8, zorder=1)
        for v, c in pts:
            ax.scatter([v], [k], s=11, color=c, ec=DARK, lw=0.3, zorder=3, clip_on=False)
    aw = RIGHT - x0 - lab_w
    xk = -(lab_w - 6.0) / aw                # the key starts over the class labels' column, past the letter
    for c, lab in zip((*_AGE_GRAYS, None), ("K = 0", "< 0.01", "\u2265 0.01", "other pairs"), strict=True):
        if c is None:
            ax.scatter([xk + 0.9 / aw], [1.07], marker="|", s=10, color=GRAY, lw=0.5, transform=ax.transAxes,
                       clip_on=False)
        else:
            ax.scatter([xk + 0.9 / aw], [1.07], s=11, color=c, ec=DARK, lw=0.3, transform=ax.transAxes,
                       clip_on=False)
        ax.text(xk + 2.0 / aw, 1.07, lab, transform=ax.transAxes, va="center", ha="left", fontsize=FS_SMALL)
        xk += (2.0 + _text_mm(lab, FS_SMALL) + 2.4) / aw


def _minor_log():
    """Minor ticks at 2-9 x each decade."""
    from matplotlib.ticker import LogLocator

    return LogLocator(base=10, subs=tuple(np.arange(2, 10) * 1.0), numticks=50)


# ---------------------------------------------------------------- panels E, F, G


def _panels_efg(page: _Page, m: _Model, y0: float, h: float) -> None:
    have = [p for p in ("E", "F", "G") if p in m.panels]
    widths = {"E": 50.0, "F": 72.0, "G": 52.0}
    total = sum(widths[p] for p in have)
    scale = (WIDTH_MM - 3.0 * (len(have) - 1)) / total
    x = 0.0
    for p in have:
        w = widths[p] * scale
        {"E": _panel_e, "F": _panel_f, "G": _panel_g}[p](page, m, x, y0, w, h)
        x += w + 3.0


def _sharing(m: _Model) -> dict | None:
    """The sharing counts of the lead trio's three genomes (stats.tsv rows sharing)."""
    t = m.trio
    want = {t.R, t.S, t.D}
    for r in m.run.stats:
        if r.test == "sharing" and set(r.scope.split(",")) == want and r.statistic == "sister_pair":
            scope = r.scope
            out = {"scope": scope, "pairs": {}}
            for a, b in ((t.R, t.S), (t.R, t.D), (t.S, t.D)):
                for name in (f"n_{a}+{b}", f"n_{b}+{a}"):
                    row = m.run.row("sharing", scope, "all", name)
                    if row is not None:
                        out["pairs"][(a, b)] = (int(row.value), name, row.n)
                        break
            return out if len(out["pairs"]) == 3 else None
    return None


def _panel_e(page: _Page, m: _Model, x0: float, y0: float, w: float, h: float) -> None:
    page.letter(x0 + 0.5, y0, m.letters["E"])
    t = m.trio
    sh = _sharing(m)
    lab_w = 17.0
    ax = page.ax(x0 + lab_w, y0 + 5.0, w - lab_w - 3.0, h - 13.0)
    _spines(ax, ("bottom",))
    if sh is None:
        ax.text(0.5, 0.5, "no sharing counts", transform=ax.transAxes, ha="center", va="center")
        return
    pairs = [(t.R, t.S), (t.R, t.D), (t.S, t.D)]
    vals = [sh["pairs"][p][0] for p in pairs]
    colors = [BLUE, ACC, "open"]  # R+S, R+D, S+D (the decoy's style in A and G)
    names = [f"{a} + {b}" for a, b in pairs]
    ax.set_ylim(2.6, -0.6)
    ax.set_yticks(range(3))
    ax.set_yticklabels(names)
    ax.tick_params(axis="y", length=0)
    hi = max(max(vals) * 1.12, 1)
    ax.set_xlim(-0.03 * hi, hi)
    for k, (v, c) in enumerate(zip(vals, colors, strict=True)):
        ax.plot([0, v], [k, k], color=LIGHT, lw=0.5, zorder=1)
        if c == "open":
            ax.scatter([v], [k], s=16, color="white", ec=ACC, lw=0.9, zorder=3, clip_on=False)
        else:
            ax.scatter([v], [k], s=16, color=c, ec=DARK, lw=0.3, zorder=3, clip_on=False)
    mean = (vals[1] + vals[2]) / 2
    ax.plot([mean, mean], [0.75, 2.35], color=DARK, lw=0.6, ls=(0, (2, 1.5)), zorder=2)
    ax.text(mean, 0.7, "equal under\nlineage sorting", ha="center", va="bottom", fontsize=FS_SMALL,
            linespacing=1.0)
    ax.text(vals[0], -0.42, "sisters", ha="center", va="bottom", fontsize=FS_SMALL)
    ax.set_xlabel("insertions shared by two genomes,\nthe third an empty site")
    from matplotlib.ticker import FuncFormatter, MaxNLocator

    ax.xaxis.set_major_locator(MaxNLocator(nbins=4, steps=[1, 2, 2.5, 5, 10]))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))


_STRATA = (("all", "all"), ("K=0", "K = 0"), ("0<K<0.01", "0 < K < 0.01"), ("K>=0.01", "K \u2265 0.01"),
           ("high_tier", "high tier"), ("nested_collapsed", "nested collapsed"),
           ("no_outgroup_filter", "no outgroup filter"))


def _dpav(m: _Model, t: _Trio) -> list[tuple[str, str, float, float, float, int, int]]:
    """(stratum, label, D, se, z, n_RD, n_SD) of the D_pav rows of trio t that exist."""
    out = []
    for st, lab in _STRATA:
        d = m.run.row("D_pav", t.scope, st, "D")
        if d is None:
            continue
        a, b = m.run.row("D_pav", t.scope, st, "n_RD"), m.run.row("D_pav", t.scope, st, "n_SD")
        out.append((st, lab, _num(d), _num(d, "se"), _num(d, "z"), int(_num(a)) if a else 0,
                    int(_num(b)) if b else 0))
    return out


def _forest_range(rows) -> tuple[float, float]:
    """The D axis range of forest rows: their points and intervals, 0 included, in quarter steps within [-1, 1],
    at least -0.25..0.25."""
    lo = min([0.0] + [r[2] - 1.96 * r[3] for r in rows if math.isfinite(r[2]) and math.isfinite(r[3])]
             + [r[2] for r in rows if math.isfinite(r[2])])
    hi = max([0.0] + [r[2] + 1.96 * r[3] for r in rows if math.isfinite(r[2]) and math.isfinite(r[3])]
             + [r[2] for r in rows if math.isfinite(r[2])])
    lo, hi = max(-1.0, math.floor((lo - 0.05) * 4) / 4), min(1.0, math.ceil((hi + 0.05) * 4) / 4)
    return min(lo, -0.25), max(hi, 0.25)


def _forest(ax, rows, with_n: bool = True, labels: bool = True, span: tuple[float, float] | None = None) -> None:
    """D_pav points with +-1.96 se, D = 0 dashed; y = row order with a gap before the check strata. span: a
    shared D range, else the rows' own."""
    ys, y = [], 0.0
    for k, row in enumerate(rows):
        if k and row[0] in ("high_tier",) and rows[k - 1][0] != "high_tier":
            y += 0.5
        ys.append(y)
        y += 1.0
    lo, hi = span if span is not None else _forest_range(rows)
    pad = 0.04 * (hi - lo)
    ax.set_xlim(lo - pad, hi + pad)
    ax.set_ylim((ys[-1] if ys else 0) + 0.6, -0.7)
    ax.axvline(0, color=DARK, lw=0.6, ls=(0, (2, 1.5)), zorder=1)
    for (st, lab, d, se, z, nrd, nsd), yy in zip(rows, ys, strict=True):
        age = st in ("K=0", "0<K<0.01", "K>=0.01")
        c = _AGE_GRAYS[("K=0", "0<K<0.01", "K>=0.01").index(st)] if age else DARK
        if math.isfinite(d) and math.isfinite(se):
            ax.plot([max(lo, d - 1.96 * se), min(hi, d + 1.96 * se)], [yy, yy], color=DARK, lw=0.7, zorder=2)
        if math.isfinite(d):
            filled = math.isfinite(z)
            ax.scatter([min(max(d, lo), hi)], [yy], s=22 if st == "all" else 14, color=c if filled else "white",
                       ec=DARK, lw=0.5, zorder=3, marker="D" if st == "all" else "o", clip_on=False)
        if with_n:
            ax.text(1.02, yy, f"{_int(nrd)} / {_int(nsd)}", transform=ax.get_yaxis_transform(), ha="left",
                    va="center", fontsize=FS_SMALL)
    ax.set_yticks(ys)
    ax.set_yticklabels([r[1] for r in rows] if labels else [""] * len(rows))
    ax.tick_params(axis="y", length=0)
    _spines(ax, ("bottom",))


def _panel_f(page: _Page, m: _Model, x0: float, y0: float, w: float, h: float) -> None:
    page.letter(x0 + 0.5, y0, m.letters["F"])
    t = m.trio
    rows = _dpav(m, t)
    lab_w, n_w = 19.0, 13.0
    ax = page.ax(x0 + lab_w, y0 + 5.0, w - lab_w - n_w, h - 13.0)
    if not rows:
        _bare(ax)
        ax.text(0.5, 0.5, "D_pav not computed", transform=ax.transAxes, ha="center", va="center")
        return
    _forest(ax, rows)
    ax.set_xlabel(r"$D_\mathrm{pav}$ = (R+D $-$ S+D) / (R+D + S+D)")
    ax.text(1.02, 1.02, "R+D / S+D", transform=ax.transAxes, ha="left", va="bottom", fontsize=FS_SMALL)
    ax.annotate("lineage sorting", xy=(0, 1.0), xycoords=("data", "axes fraction"), xytext=(0, 1.03),
                textcoords=("data", "axes fraction"), ha="center", va="bottom", fontsize=FS_SMALL)


def _panel_g(page: _Page, m: _Model, x0: float, y0: float, w: float, h: float) -> None:
    page.letter(x0 + 0.5, y0, m.letters["G"])
    t = m.trio
    lab_w = 11.0
    ax = page.ax(x0 + lab_w, y0 + 5.0, w - lab_w - 13.0, h - 13.0)
    pairs = [(r.mean_post, d.mean_post) for r, d in m.rows if d is not None]
    pairs = [(a, b) for a, b in pairs if math.isfinite(a) and math.isfinite(b)]
    g_r = _num(_paint_row(m, "posterior_mean"))
    g_d = _num(_paint_row(m, "decoy_posterior_mean"))
    for a, b in pairs:
        ax.plot([0, 1], [b, a], color=GRAY, lw=0.5, zorder=1)
        ax.scatter([0], [b], s=9, color="white", ec=ACC, lw=0.8, zorder=3)
        ax.scatter([1], [a], s=9, color=ACC, ec=DARK, lw=0.3, zorder=3)
    for xx, v in ((0, g_d), (1, g_r)):
        if math.isfinite(v):
            ax.plot([xx - 0.22, xx + 0.22], [v, v], color=DARK, lw=1.3, zorder=4)
    top = max([0.5] + [v for ab in pairs for v in ab] + [v for v in (g_r, g_d) if math.isfinite(v)])
    ax.set_ylim(0, math.ceil((top + 0.05) * 10) / 10)
    ax.set_xlim(-0.45, 1.45)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"{t.S}\n(decoy)", t.R])
    ax.set_ylabel("mean posterior of\nD-like ancestry")
    _spines(ax)
    if math.isfinite(g_r):
        ax.text(1.27, g_r, "genome-\nwide", ha="left", va="center", fontsize=FS_SMALL, linespacing=1.0)




def _h_rows(m: _Model) -> list[tuple[_Trio, list]]:
    """Panel H's trios with their D_pav rows (all, K=0, K>=0.01), sorted by D overall, largest first (nan last;
    stats order on ties)."""
    keep = ("all", "K=0", "K>=0.01")
    out = [(tr, [r for r in _dpav(m, tr) if r[0] in keep]) for tr in m.trios]
    d = [next((r[2] for r in rows if r[0] == "all"), math.nan) for _, rows in out]
    order = sorted(range(len(out)), key=lambda j: (not math.isfinite(d[j]), -d[j] if math.isfinite(d[j]) else 0, j))
    return [out[j] for j in order]


def _h_layout(m: _Model) -> dict:
    """Panel H's grid: one row per trio (H_PITCH mm), in columns of at most H_ROWS rows while each column keeps
    H_FOREST_MIN mm of D axis beside its widest label (beyond that, longer columns); its height in mm."""
    n = max(len(m.trios), 1)
    label_w = max(_text_mm(t.scope, FS_SMALL) for t in m.trios) + 1.5 if m.trios else 10.0
    gap = 3.0
    fit = max(1, int((WIDTH_MM - 1.0 + gap) // (label_w + H_FOREST_MIN + gap)))
    ncol = max(1, min(fit, -(-n // H_ROWS)))
    per = -(-n // ncol)
    col_w = (WIDTH_MM - 1.0 - (ncol - 1) * gap) / ncol
    key_h, axis_h = 4.5, 7.5
    return {"ncol": ncol, "per": per, "label_w": label_w, "col_w": col_w, "gap": gap, "key_h": key_h,
            "height": key_h + per * H_PITCH + axis_h}


def _panel_h(page: _Page, m: _Model, y0: float) -> None:
    """Every trio of the run, one row each, sorted by D_pav: its 95% interval and point (filled when testable),
    and D among the youngest (K = 0) and oldest (K >= 0.01) insertions; one D scale for every column."""
    from matplotlib.lines import Line2D

    page.letter(0.5, y0, m.letters["H"])
    lay = _h_layout(m)
    rows = _h_rows(m)
    span = _forest_range([r for _, rs in rows for r in rs])
    lo, hi = span
    pad = 0.04 * (hi - lo)
    young, old = _AGE_GRAYS[0], _AGE_GRAYS[2]
    # the key, on the panel's first line
    kx = page.ax(6.0, y0 + 0.2, WIDTH_MM - 7.0, lay["key_h"] - 0.6)
    kx.set_gid("H")
    kx.set_xlim(0, WIDTH_MM - 7.0)
    kx.set_ylim(0, 1)
    _bare(kx)
    x = 0.0
    for marker, fc, ec, text in (("D", DARK, DARK, "all insertions, 95% interval (open: not testable)"),
                                 ("o", young, DARK, "K = 0"), ("o", old, DARK, "K \u2265 0.01")):
        kx.add_line(Line2D([x + 1.0], [0.5], marker=marker, ms=2.6 if marker == "D" else 2.2, mfc=fc, mec=ec,
                           mew=0.4, ls="none"))
        kx.text(x + 2.4, 0.5, text, va="center", ha="left", fontsize=FS_SMALL)
        x += 2.4 + _text_mm(text, FS_SMALL) + 3.5
    kx.text(x, 0.5, "bold: the lead trio", va="center", ha="left", fontsize=FS_SMALL)
    top = y0 + lay["key_h"]
    for c in range(lay["ncol"]):
        mine = rows[c * lay["per"]:(c + 1) * lay["per"]]
        if not mine:
            continue
        x0 = 1.0 + c * (lay["col_w"] + lay["gap"]) + lay["label_w"]
        ax = page.ax(x0, top, lay["col_w"] - lay["label_w"], lay["per"] * H_PITCH)
        ax.set_gid("H")
        ax.set_xlim(lo - pad, hi + pad)
        ax.set_ylim(lay["per"] - 0.5, -0.5)
        ax.axvline(0, color=DARK, lw=0.6, ls=(0, (2, 1.5)), zorder=1)
        for k, (tr, rs) in enumerate(mine):
            ax.axhline(k, color=LIGHT, lw=0.3, zorder=0)
            by = {r[0]: r for r in rs}
            for st, color in (("K>=0.01", old), ("K=0", young)):
                r = by.get(st)
                if r is not None and math.isfinite(r[2]):
                    ax.scatter([min(max(r[2], lo), hi)], [k], s=4.5, color=color, ec=DARK, lw=0.3, zorder=3,
                               clip_on=False)
            r = by.get("all")
            if r is None or not math.isfinite(r[2]):
                ax.text(lo, k, "no result", ha="left", va="center", fontsize=FS_SMALL, color=GRAY, zorder=5,
                        bbox={"fc": "white", "ec": "none", "pad": 0.3})
                continue
            d, se, z = r[2], r[3], r[4]
            if math.isfinite(se):
                ax.plot([max(lo, d - 1.96 * se), min(hi, d + 1.96 * se)], [k, k], color=DARK, lw=0.7, zorder=2)
            ax.scatter([min(max(d, lo), hi)], [k], s=9, marker="D", color=DARK if math.isfinite(z) else "white",
                       ec=DARK, lw=0.5, zorder=4, clip_on=False)
        ax.set_yticks(range(len(mine)))
        ax.set_yticklabels([tr.scope for tr, _ in mine], fontsize=FS_SMALL)
        for lab, (tr, _) in zip(ax.get_yticklabels(), mine, strict=True):
            if tr == m.trio:
                lab.set_fontweight("bold")
        ax.tick_params(axis="y", length=0, pad=1.0)
        from matplotlib.ticker import MaxNLocator

        ax.xaxis.set_major_locator(MaxNLocator(nbins=4 if lay["col_w"] - lay["label_w"] > 30 else 2,
                                               steps=[1, 2, 2.5, 5, 10]))
        _spines(ax, ("bottom",))
    page.text(WIDTH_MM / 2, top + lay["per"] * H_PITCH + 3.8, r"$D_\mathrm{pav}$", ha="center", va="top", fontsize=FS)


# ---------------------------------------------------------------- the no-trio figure


def _draw_no_trio(m: _Model):
    n_pairs = len(m.rows)
    ncol = 2 if n_pairs > 3 else 1
    per = -(-n_pairs // ncol) if n_pairs else 0
    pair_h = min(9.0, 68.0 / per) if per else 0.0
    a_h = 8.0 + per * pair_h + 6.0 if per else 0.0
    b_h = 12.0 + 4.5 * len(m.no_trio["pairs"])
    H = (a_h + 5.0 if a_h else 0.0) + b_h + 3.0
    page = _Page(H)
    y = 1.0
    if a_h:
        page.letter(0.5, y, m.letters["A"])
        nt = m.no_trio
        ax = page.ax(6.0, y + 3.6, 100, 3.0)
        ax.set_xlim(0, 100)
        ax.set_ylim(0, 3)
        _bare(ax)
        ax.plot([1, 1], [0.6, 2.4], color=BLUE, lw=0.8)
        ax.text(3, 1.5, f"present in {nt['Y']}", va="center", fontsize=FS)
        ax.plot([30, 30], [0.6, 2.4], color=ACC, lw=0.8)
        ax.text(32, 1.5, f"empty site in {nt['Y']}", va="center", fontsize=FS)
        top = y + 8.0
        label_w, col_gap = 15.0, 4.0
        col_w = (RIGHT - 1.0 - ncol * label_w - (ncol - 1) * col_gap) / ncol
        for k, (r, _) in enumerate(m.rows):
            c, j = divmod(k, per)
            ax = page.ax(1.0 + label_w + c * (col_w + label_w + col_gap), top + j * pair_h, col_w, pair_h - 2.0)
            _draw_row(ax, r, m, False)
            if j == per - 1 or k == len(m.rows) - 1:
                ax.spines["bottom"].set_visible(True)
                _bp_axis(ax, 0, m.xmax, "position")
        y += a_h + 5.0
    page.letter(0.5, y, m.letters["B"])
    ax = page.ax(30.0, y + 5.0, WIDTH_MM - 60.0, b_h - 12.0)
    pairs = m.no_trio["pairs"]
    states = [("present", tuple(sorted(PRESENT, key=ALLELES.index)), DARK, "o"), ("empty site", ("E",), ACC, "o"),
              ("deletion", ("D",), GRAY, "s"), ("other sequence", ("O",), GRAY, "^"),
              ("gap or unresolved", ("G", "U", "none"), LIGHT, "X")]
    ax.set_ylim(len(pairs) - 0.4, -0.6)
    ax.set_xlim(0, 1)
    _spines(ax, ("bottom",))
    ax.set_yticks(range(len(pairs)))
    ax.set_yticklabels([f"{a}'s insertions in {b}" for a, b, _ in pairs])
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("share of calls")
    for k, (_, _, comp) in enumerate(pairs):
        tot = sum(comp.values())
        ax.axhline(k, color=LIGHT, lw=0.3, zorder=0)
        if not tot:
            ax.text(0.01, k, "no insertions", va="center", ha="left", fontsize=FS_SMALL, color=GRAY)
            continue
        for name, codes, c, mk in states:
            v = sum(comp.get(x, 0) for x in codes) / tot
            if v > 0:  # a state with no call is not drawn
                ax.scatter([v], [k], s=14, color=c, ec=DARK, lw=0.3, marker=mk, zorder=3, clip_on=False)
    for j, (name, _, c, mk) in enumerate(states):
        ax.scatter([0.02 + 0.2 * j], [1.08], s=12, color=c, ec=DARK, lw=0.3, marker=mk, transform=ax.transAxes,
                   clip_on=False)
        ax.text(0.04 + 0.2 * j, 1.08, name, transform=ax.transAxes, va="center", fontsize=FS_SMALL)
    return page.fig


# ---------------------------------------------------------------- legends


def _legends(m: _Model, nums: style.Numbers) -> tuple[list, list]:
    """(short, long): the legend pages' paragraphs, each a list of (text, style) runs; style "", "b" or "i"."""
    if m.trio is None:
        return _legends_no_trio(m, nums)
    return _Legend(m, nums).build()


def _fmt(x, style: str) -> str | None:
    """The legends' one number formatter: the text of a finite number in `style`, else None (the clause holding
    it is then dropped, or says "not available"). Styles: int (1,317), sign2 (+0.29), sign3 (+0.053), f1 (7.4),
    f2 (0.92), f3 (0.370), f4 (0.0202), pct (36%, 5.0%), pct1 (10.2%), bp (302 kb), pctl (a 0-1 fraction as a
    percentile, 52)."""
    if x is None or isinstance(x, (str, bytes)):
        return None
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(v):
        return None
    if style == "int":
        return f"{round(v):,}"
    if style in ("sign2", "sign3"):
        return f"{v:+.{style[-1]}f}"
    if style in ("f1", "f2", "f3", "f4"):
        return f"{v:.{style[1]}f}"
    if style == "pct":
        return _pct(v)
    if style == "pct1":
        return f"{100 * v:.1f}%"
    if style == "bp":
        return _bp(v)
    if style == "pctl":
        return f"{100 * v:.0f}"
    raise ValueError(f"no number style {style!r}")


def _direction(est, z, sig: float = Z_SIG) -> str:
    """The verdict on an estimate and its Z by the analyses' bar Z_SIG: "up" (estimate > 0, Z >= sig), "down"
    (estimate < 0, Z <= -sig), "none" (|Z| < sig: not shown) or "untestable" (estimate or Z missing)."""
    try:
        e, zz = float(est), float(z)
    except (TypeError, ValueError):
        return "untestable"
    if not (math.isfinite(e) and math.isfinite(zz)):
        return "untestable"
    if zz >= sig and e > 0:
        return "up"
    if zz <= -sig and e < 0:
        return "down"
    return "none"


def _compare(a, b) -> str | None:
    """"more", "less" or "equal" for two finite numbers; None when either is missing."""
    try:
        x, y = float(a), float(b)
    except (TypeError, ValueError):
        return None
    if not (math.isfinite(x) and math.isfinite(y)):
        return None
    return "more" if x > y else "less" if x < y else "equal"


def _disc(row) -> list[str]:
    """The discovering genome(s) a stats row declares in its loci field ("disc=X" or "disc=X,Y")."""
    text = getattr(row, "loci", "") or ""
    mm = re.search(r"(?:^|;)disc=([^;]+)", text)
    return [] if mm is None or mm.group(1) in ("any", "none") else mm.group(1).split(",")


def _pair_counts(note: str, labels: tuple[str, ...]) -> dict[frozenset, int]:
    """Counts written "n_A+B k" in a note, keyed by their pair of genome labels (labels may hold "+")."""
    out: dict[frozenset, int] = {}
    for mm in re.finditer(r"n_(\S+) (\d+)", note):
        tok, k = mm.group(1).rstrip(","), int(mm.group(2))
        for i, ch in enumerate(tok):
            if ch == "+" and tok[:i] in labels and tok[i + 1:] in labels:
                out[frozenset((tok[:i], tok[i + 1:]))] = k
                break
    return out


def _plain_reason(why: str) -> str:
    """trio.py's "not testable" reason in plain words."""
    why = why.strip().rstrip(".")
    mm = re.match(r"its informative loci lie in (\d+) blocks?, fewer than (\d+)", why)
    if mm:
        return (f"its informative insertions lie in {mm.group(1)} genome block(s), fewer than the "
                f"{mm.group(2)} a block jackknife needs")
    if why.startswith("no ") and why.endswith("locus"):
        return why.replace("n_RD", "R+D").replace("n_SD", "S+D").replace("locus", "insertion")
    return {"the block pseudovalues do not vary": "its estimate does not vary between genome blocks",
            "the jackknife gives no se": "no standard error could be estimated",
            "the estimate is undefined": "the estimate is undefined"}.get(why, "too few informative insertions")


# how the roles were set (stats.tsv roles row "source"): (short, long); {E} and {F} are those panels' letters. From
# the data, the sister pair is the pair sharing the most insertions (E's merged counts) and the recipient is the
# sister the donor-sharing test (D_pav's sign, F) favors (trio.trios), never E's merged counts
_HOW = {"data": ("from the data", ("from the data: the sisters are the pair sharing the most insertions (panel {E}), "
                                   "the recipient the sister with more donor-shared insertions in the donor-sharing "
                                   "test (panel {F})")),
        "--trio": ("by the user", "as given by the user"),
        "sheet order": ("in sheet order", "in sheet order, the data not identifying the recipient")}
_CHECKS = (("high_tier", "high-confidence calls only"), ("nested_collapsed", "one insertion per nested cluster"),
           ("no_outgroup_filter", "without the outgroup rule"))


class _Legend:
    """The legend texts of a trio figure. Every number passes `num` (registered as written; a missing or
    non-finite value returns None and its clause is dropped), and every verdict word comes from `_direction` or
    `_compare` on the numbers it rests on."""

    def __init__(self, m: _Model, nums: style.Numbers):
        self.m, self.n, self.run, self.t = m, nums, m.run, m.trio
        self.L = m.letters
        src = self.run.row("roles", self.t.scope, "all", "source")
        self.how_short, how_long = _HOW.get(str(src.value) if src else "", ("as set", "as set"))
        self.how_long = how_long.format(E=self.L.get("E", "E"), F=self.L.get("F", "F"))
        self.has_out = bool(self.run.cat.outgroups)

    # ------------------------------------------------ numbers

    def num(self, key: str, x, style: str, source: str) -> str | None:
        """Register x as written in `style` under key (`text`'s rule for a key already taken); None if missing."""
        text = _fmt(x, style)
        return None if text is None else self.text(key, text, source)

    def text(self, key: str, value: str, source: str) -> str:
        """Register a text value (a number as written, a coordinate, a sequence). A key already taken returns its
        text when the new one is the same (one number written in two sentences) and is a ValueError otherwise:
        one key never stands for two values."""
        old = self.n.get(key)
        if old is None:
            return self.n.add(key, value, source)
        if old != str(value):
            raise ValueError(f"numbers key {key!r} written twice with different texts: {old!r} and {value!r}")
        return old

    def stat(self, key: str, test: str, stratum: str, statistic: str, style: str, attr: str = "value",
             scope: str | None = None) -> str | None:
        scope = self.t.scope if scope is None else scope
        row = self.run.row(test, scope, stratum, statistic)
        if row is None:
            return None
        src = _src(test, scope, stratum, statistic) + ("" if attr == "value" else f":{attr}")
        return self.num(key, getattr(row, attr), style, src)

    # ------------------------------------------------ pages

    def build(self) -> tuple[list, list]:
        m, L = self.m, self.L
        f = self.facts()
        short = [[(self.headline(f), "")]]
        long_ = [[(self.definitions(), "")]]
        for p in m.panels:
            s, lg = getattr(self, f"panel_{p.lower()}")(f)
            short.append([(L[p], "b"), (s, "")])
            long_.append([(L[p], "b"), (lg, "")])
        long_.append([(self.caveats(), "")])
        long_.append([(_analysis_line(self.run.outdir), "")])
        return short, long_

    def facts(self) -> dict:
        """The statistics the verdicts rest on, read once."""
        t, run = self.t, self.run
        f: dict = {}
        d = run.row("D_pav", t.scope, "all", "D")
        f["d_row"] = d
        f["D"], f["z"] = _num(d), _num(d, "z")
        f["n_rd"] = _num(run.row("D_pav", t.scope, "all", "n_RD"))
        f["n_sd"] = _num(run.row("D_pav", t.scope, "all", "n_SD"))
        f["dir"] = "missing" if d is None else _direction(f["D"], f["z"])
        disc = _disc(d) or _disc(run.row("D_pav", t.scope, "all", "n_RD"))
        f["disc"] = disc or [t.D]
        f["fallback"] = bool(disc) and disc != [t.D]
        g = run.row("D_pav_gradient", t.scope, "K=0 minus K>=0.01", "dD")
        has_age = any(run.row("D_pav", t.scope, a, "D") is not None for a in AGE_LABELS)
        # the gradient is read in the direction of the overall excess: sign(D) x dD > 0 says the
        # excess of the sister that carries it ("owner": R when D > 0, S when D < 0) is concentrated in young insertions
        sign = 1 if f["D"] > 0 else -1 if f["D"] < 0 else 0      # nan: 0
        f["owner"] = t.R if sign > 0 else t.S if sign < 0 else None
        if g is None:
            f["grad"] = "untestable" if has_age else "no K"
        elif _direction(_num(g), _num(g, "z")) == "untestable":
            f["grad"] = "untestable"
        else:
            f["grad"] = _direction(sign * _num(g), sign * _num(g, "z")) if sign else "no excess"
        f["checks"] = [(lab, _direction(_num(run.row("D_pav", t.scope, st, "D")),
                                        _num(run.row("D_pav", t.scope, st, "D"), "z")))
                       for st, lab in _CHECKS if run.row("D_pav", t.scope, st, "D") is not None]
        f["n_q05"] = _num(run.row("paint", t.scope, "all", "n_tracts_q05")) if self.m.painted else math.nan
        ds = run.row("paint", t.scope, "all", "donor_state")   # nan (NA) exactly when the fit collapsed
        f["collapsed"] = self.m.painted and ds is not None and not math.isfinite(_num(ds))
        return f

    def counts(self, f: dict) -> tuple[str | None, str | None]:
        return (self.stat("D_pav_n_RD", "D_pav", "all", "n_RD", "int"),
                self.stat("D_pav_n_SD", "D_pav", "all", "n_SD", "int"))

    def not_significant(self, a: str | None, b: str | None) -> str:
        """The verdict at |Z| < Z_SIG: both counts and "not significant", never "equal" (a difference the test does
        not show is not shown to be absent)."""
        t = self.t
        lead = (f"Of the insertions {t.D} shares with only one of the sisters, {t.R} carries {a}, {t.S} {b}; the "
                "difference is" if a and b else
                f"The difference between {t.R} and {t.S} in the insertions {t.D} shares with one sister is")
        return f"{lead} not significant: no evidence of gene flow."

    def untestable_reason(self, f: dict) -> str:
        note = getattr(f["d_row"], "note", "") or ""
        mm = re.search(r"not testable by the block jackknife: ([^;]+)", note)
        return _plain_reason(mm.group(1)) if mm else "too few informative insertions"

    def skip_reason(self) -> str:
        """Why the donor-sharing test has no rows, in plain words."""
        t = self.t
        for section, text in _findings_lines(self.run.outdir):
            if section == "Trio tests" and f"D_pav({t.scope}) skipped" in text and "has no annotation" in text:
                return f"{t.D} has no annotation and no other annotated genome could discover its insertions"
        return (f"its counted set is empty: no insertion its discovering genome holds is present in {t.D} with "
                f"{t.R} and {t.S} both resolved" + (" and every outgroup an empty site" if self.has_out else ""))

    def headline(self, f: dict) -> str:
        """The headline sentence (`verdict`), then the caveats its reading depends on: which
        of several trios it is, and a sister pair the data do not settle."""
        s = self.verdict(f)
        n = len(self.m.trios)
        if n > 1:
            nt = self.num("trios_n", n, "int", "stats.tsv:roles")
            h = self.L.get("H", "H")
            s += (f" This trio was named by the user; panel {h} shows all {nt} of the run's trios." if self.m.named
                  else f" This is the strongest of the run's {nt} trios; panel {h} shows them all.")
        amb = self.ambiguity()
        return s + (f" {amb}" if amb else "")

    def ambiguity(self) -> str | None:
        """When the sister pair is uncertain (the sharing test's p_top >= P_WARN, or nothing shared): a sentence with
        the two largest sharing counts, saying the roles may be swapped; else None."""
        sh = _sharing(self.m)
        if sh is None:
            return None
        p = self.run.row("sharing", sh["scope"], "all", "p_top")
        if p is None or _num(p) < P_WARN:
            return None
        ranked = sorted(sh["pairs"].items(), key=lambda kv: -kv[1][0])
        (a1, a2), (n1, name1, _) = ranked[0]
        (b1, b2), (n2, name2, _) = ranked[1]
        k1 = self.num("sister_top", n1, "int", _src("sharing", sh["scope"], "all", name1))
        k2 = self.num("sister_second", n2, "int", _src("sharing", sh["scope"], "all", name2))
        t = self.t
        return (f"Which two of {t.R}, {t.S} and {t.D} are sisters is uncertain: {a1} and {a2} share {k1} insertions, "
                f"{b1} and {b2} {k2}, a difference within chance, so the roles may be swapped.")

    def verdict(self, f: dict) -> str:
        """The headline sentence: the donor-sharing test's verdict, its age profile and the painting."""
        m, t = self.m, self.t
        R, S, D = t.R, t.S, t.D
        a, b = self.counts(f)
        if f["dir"] == "missing":
            return (f"The donor-sharing test of {R}, {S} and {D} has no result: {self.skip_reason()}. The figure "
                    "shows the genotypes that exist.")
        if f["dir"] == "untestable":
            if f["n_rd"] == 0 and f["n_sd"] == 0:
                return (f"No insertion counted by the donor-sharing test is shared by {D} with exactly one of {R} "
                        f"and {S}, so the test has nothing to compare.")
            return f"The donor-sharing test of {R}, {S} and {D} could not be run: {self.untestable_reason(f)}."
        if f["dir"] == "down":   # the excess is S's: its age profile reads in that direction
            return (f"{S}, not {R}, carries more of the insertions {D} shares with one sister "
                    + (f"({b} against {a})" if a and b else "") + f": the test points against the roles, which "
                    f"were set {self.how_short}."
                    + (f" {S}'s excess is concentrated in young insertions, consistent with gene flow from {D}-like "
                       f"donors into {S} rather than incomplete lineage sorting alone." if f["grad"] == "up" else ""))
        if f["dir"] == "none":
            return self.not_significant(a, b)
        ratio = (self.num("D_pav_ratio", f["n_rd"] / f["n_sd"], "f1", _src("D_pav", t.scope, "all", "n_RD")
                          + " / n_SD") if f["n_sd"] else None)
        s = (f"Of the insertions {D} shares with only one of the sisters, {R} carries "
             + (f"{ratio} times as many as {S}" if ratio else f"more than {S}"))
        # D_pav's null Z is not calibrated, so a gene-flow verdict says "consistent with"
        s += {"up": f", most of all among young insertions: consistent with gene flow from {D}-like donors into {R} "
                    "rather than incomplete lineage sorting alone.",
              "down": ", concentrated in older insertions: more than lineage sorting alone is expected to give, but "
                      "not the young excess expected from recent gene flow.",
              "none": ": more than lineage sorting alone is expected to give; young and old insertions do not differ "
                      "clearly, so the test does not date it.",
              "untestable": ": more than lineage sorting alone is expected to give; the age gradient could not be "
                            "tested."
              }.get(f["grad"], ": more than lineage sorting alone is expected to give (insertion ages unknown).")
        if m.painted and math.isfinite(f["n_q05"]):
            nq = self.stat("n_tracts_q05", "paint", "all", "n_tracts_q05", "int")
            pm = self.stat("posterior_mean", "paint", "all", "posterior_mean", "f3")
            dm = self.stat("decoy_posterior_mean", "paint", "all", "decoy_posterior_mean", "f3")
            if f["n_q05"] > 0 and nq:
                s += f" {nq} stretch(es) of {R}'s genome stand out as {D}-like beyond the sister's control painting."
            elif pm and dm:
                s += (f" No stretch of {R}'s genome stands out from the sister's control painting; their mean "
                      f"posteriors are {pm} and {dm}.")
        return s

    def definitions(self) -> str:
        m, t, run = self.m, self.t, self.run
        cat = run.cat
        out = (f"R, S and D are the recipient {_role(m, t.R)}, its sister {_role(m, t.S)} and the donor "
               f"{_role(m, t.D)}, roles {self.how_long}")
        if cat.outgroups:
            out += f"; O is the outgroup {', '.join(_role(m, o) for o in cat.outgroups)}"
        nl = self.stat("n_loci", "catalog", "all", "n_loci", "int", scope="all")
        out += (". Each of the " + (f"{nl} " if nl else "") + "loci is one LTR-retrotransposon insertion, genotyped "
                "in every genome as present, an empty site (the intact target site with one copy of the target-site "
                "duplication, TSD) or missing (deleted, other sequence, gap or unresolved; never counted as absent). "
                "R+D: present in R and D, an empty site in S; S+D and R+S likewise. K: the LTR-LTR divergence (K2P) of "
                "an insertion's copy, a relative age. Intervals: 95%, ±1.96 block-jackknife standard errors")
        d = run.row("D_pav", t.scope, "all", "D")
        mm = re.search(r"of ([\d,]+) bp", d.note) if d is not None and isinstance(d.note, str) else None
        if mm:
            bs = self.num("block_bp", float(mm.group(1).replace(",", "")), "bp", _src("D_pav", t.scope, "all", "D")
                          + ":note")
            out += f" over {bs} blocks"
        return out + "."

    # ------------------------------------------------ panels

    def panel_a(self, f: dict) -> tuple[str, str]:
        m, t, n = self.m, self.t, self.n
        R, S, D = t.R, t.S, t.D
        p, q = m.target, m.decoy
        nd, nnd = n.get("sequences_drawn"), n.get("sequences_not_drawn")
        src_r = f"{p.file or 'loci.tsv'}:class"
        rd_r, rs_r = int((p.cls == RD).sum()), int((p.cls == RS).sum())
        a_rd, a_rs = self.num("A_RD_target", rd_r, "int", src_r), self.num("A_RS_target", rs_r, "int", src_r)
        sh_r = rd_r / (rd_r + rs_r) if rd_r + rs_r else math.nan
        a_sh = self.num("A_RD_share_target", sh_r, "pct", src_r)
        lg = (f"{R}'s {nd} longest sequences (each ≥ {MIN_FRAC:.0%} of the longest, at most {MAX_ROWS}; "
              f"{nnd} shorter ones not drawn) to scale")
        short = f"{R}'s {m.text.get('unit', 'sequence')}s"
        sh_d = d_sh = None
        decoy = m.painted and q is not None and q.file is not None   # a decoy painting, not classes only
        if q is not None:
            src_d = f"{q.file or 'loci.tsv'}:class"
            rd_d, rs_d = int((q.cls == RD).sum()), int((q.cls == RS).sum())
            d_rd, d_rs = self.num("A_RD_decoy", rd_d, "int", src_d), self.num("A_RS_decoy", rs_d, "int", src_d)
            sh_d = rd_d / (rd_d + rs_d) if rd_d + rs_d else math.nan
            d_sh = self.num("A_RD_share_decoy", sh_d, "pct", src_d)
            lg += (f", each above its {S} counterpart (its syntenic partner, the sequence sharing most loci with it) "
                   + (f"painted the same way as a control (the decoy: {R} and {S} swapped)" if decoy else
                      f"with {S}'s insertions classed the same way as a control ({R} and {S} swapped; nothing is "
                      "painted)")
                   + f". Ticks: every informative insertion, R+D above in vermillion (S+D in the control), R+S below "
                   f"in blue: {R} {a_rd} and {a_rs}, {S} {d_rd} and {d_rs}.")
            if a_sh and d_sh:
                lg += f" {D}-shared insertions make up {a_sh} of {R}'s informative insertions and {d_sh} of {S}'s."
            short += (f", each above the same painting of its sister {S} as a control" if decoy else
                      f", each above its sister {S}'s, whose insertions are classed the same way as a control")
        else:
            lg += f". Ticks: every informative insertion, R+D above in vermillion, R+S below in blue ({a_rd} and {a_rs})."
        if not self.run.cat.annotated.get(R):
            lg += (f" {R} has no annotation, so it holds no insertions of its own to class or paint (a painting reads "
                   "the insertions the recipient discovers).")
            short += f"; {R}, without an annotation, has no insertions of its own to show"
        word = _compare(sh_r, sh_d)
        if word and a_sh and d_sh:
            short += (f": insertions shared with {D} (vermillion) make up "
                      + {"more": f"a larger share of {R}'s informative insertions ({a_sh}) than of the control's",
                         "less": f"a smaller share of {R}'s informative insertions ({a_sh}) than of the control's",
                         "equal": f"the same share of {R}'s informative insertions ({a_sh}) as of the control's"}[word]
                      + f" ({d_sh})")
        elif a_sh:
            short += f": insertions shared with {D} (vermillion) make up {a_sh} of {R}'s informative insertions"
        k = sum(len(r.rd) > 0 for r, _ in m.rows)
        if m.rows and a_sh:
            if k == len(m.rows):
                short += " and lie on every drawn sequence"
            else:
                kk = self.num("A_rows_with_RD", k, "int", src_r)
                short += f" and lie on {kk} of the {nd} drawn sequences"
        if m.painted:
            ell = self.stat("mean_tract_bp", "paint", "all", "mean_tract_bp", "bp")
            nt = self.stat("n_tracts", "paint", "all", "n_tracts", "int")
            ndt = self.stat("n_decoy_tracts", "paint", "all", "n_decoy_tracts", "int")
            nq = self.stat("n_tracts_q05", "paint", "all", "n_tracts_q05", "int")
            bins = self.text("A_bin", _bp(m.bin_bp), "display:bin, 1/250 of the longest drawn sequence")
            lg += (f" Bar color: the posterior of {D}-like ancestry from the painting (a two-state hidden Markov "
                   "model over the informative insertions" + (f", switch scale {ell}" if ell else "")
                   + (f", fitted on {R} and applied unchanged to {S}" if decoy else f", fitted on {R}")
                   + f"), averaged per {bins} display bin; hatched: bins without an informative insertion; dark line "
                   f"on the bar: a tract (posterior > 0.5 and posterior odds ≥ {BF_TRACT:g} times the prior odds"
                   + (f"; {nt} in {R}" if nt else "") + (f", {ndt} in {S}" if ndt else "") + ")"
                   + ((f", black when q ≤ {Q_TRACT:g} against the decoy" + (f" ({nq})" if nq else "") + ".") if decoy
                      else f"; with no decoy painting ({S} discovers no insertions to paint), no tract has a q-value."))
            if f["collapsed"]:   # a fit at its floor holds no data in its donor-like state (paint.Painting.collapsed)
                lg += (f" The fit collapsed: no {D}-like ancestry is detectable along {R}'s genome, so the posterior "
                       "is flat and no tract is called.")
                short += f". No {D}-like ancestry is detectable along {R}'s genome."
            elif nt and ndt and nq:   # the short legend gives no threshold
                short += (f". {R} has {nt} donor-like tracts and the control {ndt}; "
                          + (f"none of {R}'s stands" if f["n_q05"] == 0 else f"{nq} of {R}'s stand")
                          + " out from the control's tracts.")
            else:
                short += "."
        else:
            msg = n.get("no_painting") or "no current painting for this trio"
            lg += (f" There is {msg}, so no posterior or tract is drawn." if self.run.cat.annotated.get(R)
                   else f" No posterior or tract is drawn ({msg}).")
            short += "."
        inv = m.text.get("inv", {})
        drawn, sparse = inv.get("drawn", 0), inv.get("sparse", 0)
        if drawn or sparse:
            k1 = self.num("A_inv_drawn", drawn, "int", "inversions.tsv:n_loci")
            lg += (" Gray bands: blocks inverted against the trio's other genomes, drawn where their loci average at "
                   f"least one per display bin ({k1} intervals")
            if sparse:
                k2 = self.num("A_inv_sparse", sparse, "int", "inversions.tsv:n_loci")
                span, nl = inv["largest"]
                k3 = self.text("A_inv_sparse_largest", f"{_bp(span)} on {_int(nl)} loci",
                               "inversions.tsv:start_x,end_x,n_loci")
                lg += f"; {k2} sparser blocks not drawn, the largest {k3}"
            lg += "); the painting flags them but does not use them."
        return short, lg

    def panel_b(self, f: dict) -> tuple[str, str]:
        m, t, n = self.m, self.t, self.n
        R, S = t.R, t.S
        w = m.window
        tr = w.get("tract")
        src_t = f"{m.target.tract_file}:"
        ex = m.example is not None
        if tr is not None:
            meds, npool = m.text["tract_medians"]
            mn = self.num("B_median_n_loci", meds["n_loci"], "int", src_t + "n_loci")
            ml = self.num("B_median_length", meds["length"], "bp", src_t + "end - start0")
            mp = self.num("B_median_mean_post", meds["mean_post"], "f2", src_t + "mean_post")
            npool_t = self.num("B_tracts_pool", npool, "int", src_t + ("q" if m.text["pool_q05"] else "name"))
            pool = f"tracts at q ≤ {Q_TRACT:g}" if m.text["pool_q05"] else "tracts"
            pos = self.text("B_tract", f"{w['name']}:{_int(tr['start'])}-{_int(tr['end'])}",
                            src_t + "chrom, start0 + 1, end")
            tl = self.num("B_tract_length", tr["end"] - tr["start"] + 1, "bp", src_t + "end - start0")
            tn = self.num("B_tract_n_loci", tr["n_loci"], "int", src_t + "n_loci")
            trd = self.num("B_tract_n_RD", tr["n_RD"], "int", src_t + "n_RD")
            trs = self.num("B_tract_n_RS", tr["n_RS"], "int", src_t + "n_RS")
            tq = self.num("B_tract_q", tr["q"], "f2", src_t + "q")
            tmp = self.num("B_tract_mean_post", tr["mean_post"], "f2", src_t + "mean_post")
            seen = self.num("B_tracts_examined", m.text["tracts_seen"], "int", src_t + "name")
            meds_t = _and([x for x in (f"loci ({mn})" if mn else "", f"length ({ml})" if ml else "",
                                       f"mean posterior ({mp})" if mp else "") if x])
            lg = f"The {pool} of {R} ({npool_t}) were ranked by their distance from the medians of {meds_t}"
            if ex:
                lg += (f" and taken in that order until one held a verified R+D insertion ({seen} examined"
                       + (f"; panel {self.L['C']}" if "C" in self.L else "") + ")")
            else:
                scr = self.num("C_screened", m.text["screened"], "int", f"{m.target.file}:class")
                lg += (f"; no tract held a structurally verified R+D insertion ({scr} screened in {seen} tracts), so "
                       "the window is the first tract in that order and no base-pair example is drawn")
            lg += (f". Shown with one tract length on each side: {pos} (" + ", ".join(
                x for x in (tl, f"{tn} loci" if tn else "", f"{trd} R+D" if trd else "", f"{trs} R+S" if trs else "",
                            f"mean posterior {tmp}" if tmp else "", f"q = {tq}" if tq else "") if x) + ").")
            pc = w.get("pct", {})
            parts = [(k, lab) for k, lab in (("n_loci", "loci"), ("length", "length"), ("mean_post", "mean posterior"))]
            col = {"length": "end - start0"}   # the BED columns each percentile ranks
            pts = [f"{lab} {self.num(f'B_pct_{k}', pc.get(k), 'pctl', src_t + col.get(k, k))}" for k, lab in parts
                   if _fmt(pc.get(k), "pctl")]
            if pts:
                lg += f" Its percentiles among {R}'s tracts: " + _and(pts)
                pdm = self.num("B_pct_decoy_mean_post", pc.get("decoy_mean_post"), "pctl",
                               f"{m.decoy.tract_file}:mean_post" if m.decoy is not None and m.decoy.tract_file
                               else src_t + "mean_post")
                lg += (f"; its mean posterior's among the decoy's tracts: {pdm}." if pdm else ".")
            typical = all(0.25 <= pc.get(k, math.nan) <= 0.75 for k in ("n_loci", "length", "mean_post"))
            short = ("One typical donor-like stretch" if typical
                     else "One donor-like stretch, the nearest to typical that holds a verified insertion" if ex
                     else "One donor-like stretch, the nearest to typical")
        else:
            nt = n.get("no_tract") or "no tract called"
            lg = (f"With {nt}, the window is centered on "
                  + ("the example insertion." if ex else "the R+D insertion nearest the median K."))
            short = "One window around " + ("the example insertion" if ex else "an R+D insertion")
        short += ", insertion by insertion: which genomes carry each insertion and which have its empty site"
        lg += (f" Each {R} insertion is a box at its position (nested ones stacked); each row gives its state in one "
               f"genome from {R}'s own calls: filled present, open empty site, hatched missing; vermillion R+D, blue "
               "R+S, gray other (both, neither, missing" + (", or carried by the outgroup)." if self.has_out else ")."))
        if m.painted:
            lg += " Above: the posterior per insertion" + (" and the tract." if tr is not None else ".")
            hi_r = self.text("B_high_target", _share_high(m.target.post), f"{m.target.file}:post")
            if m.decoy is not None and m.decoy.file and np.isfinite(m.decoy.post).any():
                share_d = float(np.mean(m.decoy.post[np.isfinite(m.decoy.post)] >= 0.9))
                share_r = float(np.mean(m.target.post[np.isfinite(m.target.post)] >= 0.9))
                hi_d = self.text("B_high_decoy", _share_high(m.decoy.post), f"{m.decoy.file}:post")
                lg += (f" Right, on the same posterior axis: every painted insertion genome-wide in bins of 0.1, {R} "
                       f"(vermillion) against {S}'s painting with the roles swapped (gray). Posteriors of 0.9 or more "
                       f"are {hi_r} of {R}'s insertions and ")
                if share_d > 0 and share_r > 0:
                    lg += (f"{hi_d} of the decoy's: plateaus this high also occur in the sister's control painting "
                           f"(roles swapped), {_ratio_words(share_d, share_r)}"
                           + ("; the tract's q weighs this one against the decoy's tracts." if tr is not None else "."))
                    short += "; posteriors this high also occur in the control painting"
                elif share_d > 0:
                    lg += (f"{hi_d} of the decoy's: the sister's control painting (roles swapped) reaches 0.9, while "
                           f"{R}'s painting never does.")
                    short += f"; the control painting reaches posteriors {R}'s never does"
                else:
                    lg += f"{hi_d} of the decoy's: the sister's control painting (roles swapped) never reaches 0.9."
                    short += "; the control painting never reaches such posteriors"
            else:
                lg += (f" Right: {R}'s posteriors in bins of 0.1 ({hi_r} reach 0.9 or more); with no decoy painting, "
                       "no control is drawn.")
        return short + ".", lg

    def panel_c(self, f: dict) -> tuple[str, str]:
        m, t = self.m, self.t
        R, S, D = t.R, t.S, t.D
        ex = m.example
        lid = self.text("C_lid", ex["lid"], "loci.tsv:lid")
        pos = self.text("C_site", f"{ex['seqid']}:{_int(ex['start'])}-{_int(ex['end'])}",
                        f"work/elements/{R}.tsv:seqid,start,end")
        ln = self.num("C_length", ex["length"], "int", f"work/elements/{R}.tsv:start,end")
        k = self.num("C_K", ex["k"], "f4", f"loci.tsv:{R}_k")
        mk = self.num("C_median_K_RD", m.text["med_k"], "f4", f"loci.tsv:{R}_k")
        nrd = self.num("C_n_RD", m.text["n_rd"], "int", f"{m.target.file or 'loci.tsv'}:class")
        tsd = self.text("C_TSD", ex["tsd"], f"work/elements/{R}.tsv:tsd")
        q25, q75, _, what = ex["fam"]
        fr = self.text("C_family_IQR", f"{_int(q25)}-{_int(q75)}", f"work/elements/{R}.tsv:start,end")
        l1, l2 = ex["ltr"]
        ltr = self.text("C_LTRs", f"{_int(l1)} and {_int(l2)}", f"work/elements/{R}.tsv:ltr_left,ltr_right")
        ds = ex["dsite"]
        dpos = self.text("C_donor_site", f"{ds['seqid']}:{_int(ds['start'])}-{_int(ds['end'])}",
                         (f"work/elements/{D}.tsv:seqid,start,end" if ds["annotated"] else
                          f"pairs/{R}{PAIR_SEP}{D}.calls.tsv.gz:b_seqid,b_start,b_end"))
        s = ex["s"]
        spos = self.text("C_sister_site", f"{s['seqid']}:{_int(s['start'])}-{_int(s['end'])}",
                         f"pairs/{R}{PAIR_SEP}{S}.calls.tsv.gz:b_seqid,b_start,b_end")
        il, ir = ex["s_ident"]
        idt = self.text("C_flank_identity", f"{round(il * FLANK_CHECK)}/{FLANK_CHECK} and "
                        f"{round(ir * FLANK_CHECK)}/{FLANK_CHECK}", f"genome sequence:{R},{S}")
        src_chk = f"genome sequence:{R},{S},{D} (structural checks)"
        scr = self.num("C_screened", m.text["screened"], "int", f"{m.target.file or 'loci.tsv'}:class")
        npass = self.num("C_passed", m.text["passed"], "int", src_chk)
        fails = sorted(m.text["checker"].fails.items(), key=lambda kv: (-kv[1], kv[0]))[:3]
        why = ", ".join(f"{label} ({self.num(f'C_fail_{i}', c, 'int', src_chk)})" for i, (label, c) in
                        enumerate(fails, 1))
        where = (f"in the tract of panel {self.L['B']}" if "B" in self.L and m.window and m.window.get("tract")
                 else "among all R+D insertions")
        rule = (f"the one with K nearest the median K of all {nrd} R+D insertions ({mk})" if mk
                else "the first by position (no K is known)")
        lg = (f"Insertion {lid} ({pos}; {ln} bp, family {ex['family']}" + (f", K = {k}" if k else "")
              + f"): of the R+D insertions {where} that passed every check, {rule}. In {R}'s orientation: "
              f"{FLANK_SHOW} bp of flank, the TSD (box) and {LTR_SHOW} bp of each LTR end (dark) in {R} and in "
              f"{D}'s copy ({dpos}), and {S}'s empty site with one TSD copy ({spos}); ticks join identical bases; "
              f"above, the element to scale. Checks: TSD copies identical to {R}'s annotated "
              f"TSD ({tsd}) in {R} and {D}; ends TG...CA in both; both LTRs whole ({ltr} bp); length within the "
              f"{what}'s interquartile range ({fr} bp); {S}'s TSD copy identical and its flanks "
              f"≥ {FLANK_IDENT:.0%} identical to {R}'s over {FLANK_CHECK} bp ({idt}). Of the {scr} R+D "
              f"insertions screened in that order, {npass} passed every check, so the example is drawn from that "
              "clean subset" + (f"; commonest failures: {why}." if why else "."))
        short = (f"One insertion shared by {R} and {D}, base by base: the element between two identical target-site "
                 f"copies in both, and the undisturbed empty site with a single copy in {S}.")
        return short, lg

    def panel_d(self, f: dict) -> tuple[str, str]:
        m, t = self.m, self.t
        R, S = t.R, t.S
        rows = _absence_rows(m.run, t)
        src = f"mechanisms.tsv:absence/{R}>{S}/all/all/" + "{age}/{metric}"
        total = rows.get("all", {})
        na = self.num("D_n_absent", total.get("n_absent"), "int", src.format(age="all", metric="n_absent"))
        ages = []
        for a, lab in zip(AGE_LABELS, ("K = 0", "0 < K < 0.01", "K ≥ 0.01"), strict=True):
            v = self.num(f"D_n_absent_{a}", rows.get(a, {}).get("n_absent"), "int", src.format(age=a, metric="n_absent"))
            if v:
                ages.append(f"{lab}: {v}")
        e_all = total.get("frac_E", math.nan)
        ea = self.num("D_E_all", e_all, "pct", src.format(age="all", metric="frac_E"))
        e_y, e_o = rows.get("K=0", {}).get("frac_E", math.nan), rows.get("K>=0.01", {}).get("frac_E", math.nan)
        ey = self.num("D_E_young", e_y, "pct", src.format(age="K=0", metric="frac_E"))
        eo = self.num("D_E_old", e_o, "pct", src.format(age="K>=0.01", metric="frac_E"))
        lo = self.num("D_lower_bound", _absence_value(total, "frac_D_inside+frac_O_inside"), "pct",
                      src.format(age="all", metric="frac_D_inside+frac_O_inside"))
        up = self.num("D_upper_bound", 1 - e_all, "pct", src.format(age="all", metric="1-frac_E"))
        es = [v for v in (o.get("frac_E", math.nan) for o in _absence_others(m.run, t).values()) if math.isfinite(v)]
        rng = (self.text("D_E_other_pairs", _pct(min(es)) if _pct(min(es)) == _pct(max(es))
                         else f"{_pct(min(es))}-{_pct(max(es))}", "mechanisms.tsv:absence/<ingroup pairs>/all/all/all/frac_E")
               if es else None)
        used = ", ".join(self.L[x] for x in ("A", "B", "E", "F") if x in self.L)
        lg = ((f"The {na} absences" if na else "The absences") + f" of {R} insertions from {S}"
              + (f" ({'; '.join(ages)})" if ages else "")
              + f", by class: empty site ({S} never carried it); a deletion or other sequence ending inside the "
              "element beside a remnant of it (shown derived), or with the remnant not shown; a deletion past the "
              "element's ends; other sequence at the site; complex other sequence. Gray ticks: the other ordered "
              "ingroup pairs, all ages.")
        if ea:
            lg += f" Empty sites are {ea} of absences" + (f" (other pairs {rng})" if rng else "")
            if ey and eo:
                lg += f", {ey} of the youngest and {eo} of the oldest"
            lg += "."
        if lo and up:
            lg += f" The derived share lies between {lo} (shown) and {up} (not shown ancestral)."
        lg += f" Only empty sites count as absent in {used}."
        short = (f"{ea} of absences are undisturbed empty sites" if ea else "Absences by class")
        word = _compare(e_o, e_y)
        if ey and eo and word:
            short += {"less": f", fewer among old insertions ({eo}) than among young ones ({ey})",
                      "more": f", more among old insertions ({eo}) than among young ones ({ey})",
                      "equal": f", as many among old insertions as among young ones ({eo})"}[word]
        short += "; the rest are deletions or other sequence, never counted as absences."
        return short, lg

    def panel_e(self, f: dict) -> tuple[str, str]:
        m, t = self.m, self.t
        R, S, D = t.R, t.S, t.D
        sh = _sharing(m)
        if sh is None:
            return "No sharing counts.", "The sharing counts were not computed."
        val = {pair: v for pair, (v, _, _) in sh["pairs"].items()}
        out = {pair: self.num(f"E_{pair[0]}+{pair[1]}", v, "int", _src("sharing", sh["scope"], "all", name))
               for pair, (v, name, _) in sh["pairs"].items()}
        nl = sh["pairs"][(R, S)][2]
        nl_t = self.num("E_n_loci", nl, "int", _src("sharing", sh["scope"], "all", sh["pairs"][(R, S)][1]) + ":n")
        where = ((f", on the {nl_t} loci resolved in all three" if nl_t else "")
                 + (" and not present in the outgroup" if self.has_out and nl_t else ""))
        mt = self.num("E_expected", (val[(R, D)] + val[(S, D)]) / 2, "int",
                      _src("sharing", sh["scope"], "all", "n_*") + " minority mean")
        top = val[(R, S)] >= max(val[(R, D)], val[(S, D)])
        lg = (f"Insertions in exactly two of {R}, {S} and {D}, the third an empty site{where} (merged calls: a "
              f"description; the test is panel {self.L['F']}): {R} + {S} {out[(R, S)]} ("
              + ("the most: the sisters" if top else f"the sisters by roles set {self.how_short}, not the most")
              + f"), {R} + {D} {out[(R, D)]}, {S} + {D} {out[(S, D)]}; colored as in "
              + _panels([self.L[x] for x in ("A", "B") if x in self.L])
              + f" (blue R+S, vermillion R+D, open vermillion S+D, the decoy's donor-shared class). Lineage sorting "
              f"alone makes R+D and S+D equal (dashed: their mean, {mt}).")
        kk = self.run.row("sharing", sh["scope"], "all", "p_kksc")
        pc = _pair_counts(kk.note, self.run.cat.labels) if kk is not None and isinstance(kk.note, str) else {}
        rd_k, sd_k = pc.get(frozenset((R, D))), pc.get(frozenset((S, D)))
        if rd_k is not None and sd_k is not None:
            src = _src("sharing", sh["scope"], "all", "p_kksc") + ":note"
            a = self.num("E_kksc_RD", rd_k, "int", src)
            b = self.num("E_kksc_SD", sd_k, "int", src)
            lg += f" From {D}'s own calls they are {a} (R+D) and {b} (S+D)."
        word = _compare(val[(R, D)], val[(S, D)])
        short = (f"{R} and {S} share the most insertions ({out[(R, S)]}: the sisters); " if top else
                 f"{R} and {S}, the sisters by roles set {self.how_short}, share {out[(R, S)]}, not the most; ")
        short += {"more": f"{R}, not {S}, shares more with {D} ({out[(R, D)]} against {out[(S, D)]})",
                  "less": f"{S}, not {R}, shares more with {D} ({out[(S, D)]} against {out[(R, D)]})",
                  "equal": f"{R} and {S} share equally with {D} ({out[(R, D)]})"}[word]
        # from the data, the recipient is set by D_pav's sign on the donor's own calls, not by these merged counts
        short += (f"; the donor-sharing test (panel {self.L['F']}) makes {R} the recipient."
                  if self.how_short == "from the data" else ".")
        return short, lg

    def panel_f(self, f: dict) -> tuple[str, str]:
        m, t = self.m, self.t
        R, S, D = t.R, t.S, t.D
        rows = _dpav(m, t)
        if not rows:
            return (f"The donor-sharing test has no result: {self.skip_reason()}.",
                    f"D_pav has no result for this trio: {self.skip_reason()}.")
        disc = ", ".join(f["disc"])
        cnt = m.run.row("D_pav", t.scope, "all", "n_RD")
        nl = self.num("F_n_loci", cnt.n if cnt is not None else None, "int",
                      _src("D_pav", t.scope, "all", "n_RD") + ":n")
        a = self.stat("F_resolved_R", "resolution", "all", f"resolved_{R}", "pct")
        b = self.stat("F_resolved_S", "resolution", "all", f"resolved_{S}", "pct")
        o = self.stat("F_outgroup_removed", "outgroup_filter", "all", "removed_frac", "pct") if self.has_out else None
        rd, sd = self.counts(f)
        lg = ("D_pav = (n_R+D − n_S+D) / (n_R+D + n_S+D) on the " + (f"{nl} " if nl else "")
              + f"loci {disc} discovered"
              + (f" (fallback discovery: {D} has no annotation)" if f["fallback"] else "")
              + f", with {D} present, {R} and {S} resolved" + (" and every outgroup an empty site" if self.has_out
                                                                else "")
              + f", every state from {disc}'s own calls")
        extra = [x for x in ((f"{R} resolved at {a}, {S} at {b}" if a and b else ""),
                             (f"the outgroup rule removed {o} of candidates" if o else "")) if x]
        lg += (f" ({'; '.join(extra)})" if extra else "") + ": "
        dv = self.stat("D_pav", "D_pav", "all", "D", "sign2")
        zv = self.stat("D_pav_z", "D_pav", "all", "D", "f1", attr="z")
        if dv:
            ci = None
            if math.isfinite(f["D"]) and math.isfinite(_num(f["d_row"], "se")):
                se = _num(f["d_row"], "se")
                ci = self.text("D_pav_CI", f"{f['D'] - 1.96 * se:+.2f} to {f['D'] + 1.96 * se:+.2f}",
                               _src("D_pav", t.scope, "all", "D") + ":value+-1.96se")
            lg += dv + (f" (95% interval {ci}; Z {zv}" if ci and zv else " (not testable by the block jackknife")
        else:
            lg += "not available (no informative insertion"
        lg += (f"; {rd} R+D, {sd} S+D)." if rd and sd else ").")
        parts = {}
        for st, lab, dv2, _, _, nrd, nsd in rows:
            if st == "all":
                continue
            key = "F_" + re.sub(r"[^A-Za-z0-9]+", "_", st).strip("_")
            d_txt = self.num(key, dv2, "sign2", _src("D_pav", t.scope, st, "D"))
            nn = self.text(key + "_n", f"{_int(nrd)}/{_int(nsd)}", _src("D_pav", t.scope, st, "n_RD") + ",n_SD")
            parts[st] = f"{dict(_STRATA)[st]} " + (f"{d_txt} ({nn})" if d_txt else f"not available ({nn})")
        age = [parts[x] for x in AGE_LABELS if x in parts]
        if age:
            lg += f" By the K of {disc}'s copy: " + ", ".join(age) + "."
            if f["grad"] in ("up", "down", "none", "no excess"):   # a finite Z: numbers written with their verdict
                dd = self.stat("D_pav_gradient", "D_pav_gradient", "K=0 minus K>=0.01", "dD", "sign2")
                dz = self.stat("D_pav_gradient_z", "D_pav_gradient", "K=0 minus K>=0.01", "dD", "f1", attr="z")
                lg += f" Youngest minus oldest: {dd} (Z {dz})"
                if f["grad"] == "no excess":
                    lg += "; with no excess overall, there is no excess whose age profile it could show."
                else:   # read in the direction of the excess
                    lg += (f", read in the direction of the excess, {f['owner']}'s: " + {
                        "up": "concentrated in young insertions.", "down": "concentrated in old insertions.",
                        "none": "no clear difference between young and old."}[f["grad"]])
            else:
                lg += " The age gradient could not be tested."
        lg += " Lineage sorting gives D = 0 at every age (dashed)."
        names = dict(_CHECKS)
        chk = [parts[st].replace(dict(_STRATA)[st], names[st], 1) for st in names if st in parts]
        if chk:
            lg += " Checks: " + ", ".join(chk) + "."
        if any(math.isfinite(r[2]) and not math.isfinite(r[4]) for r in rows):
            lg += " Open points: not testable by the block jackknife."
        # the short sentence: the verdict, its age profile and its checks, each from its own statistic
        if f["dir"] == "up":
            short = f"{R} carries more of the insertions {D} shares with one sister ({rd} against {sd})"
            if f["grad"] == "up":   # its D values stay on the long legend
                short += ", most of all young ones"
            elif f["grad"] == "down":
                short += ", most of all old ones"
            elif f["grad"] == "none":
                short += ", with no clear difference between young and old"
            same = [lab for lab, d in f["checks"] if d == "up"]
            if f["checks"]:
                if len(same) == len(f["checks"]):
                    short += f", and the excess holds in all {len(f['checks'])} checks."
                else:
                    miss = [lab for lab, d in f["checks"] if d != "up"]
                    short += (f", and the excess holds in {len(same)} of {len(f['checks'])} checks (not "
                              + ", ".join(miss) + ").")
            else:
                short += "."
        elif f["dir"] == "down":
            short = (f"{S}, not {R}, carries more of the insertions {D} shares with one sister ({sd} against {rd})"
                     + {"up": ", most of all young ones", "down": ", most of all old ones"}.get(f["grad"], "")
                     + f": the test points against the roles, which were set {self.how_short}.")
        elif f["dir"] == "none":
            short = self.not_significant(rd, sd)
        else:
            short = f"The donor-sharing test could not be run: {self.untestable_reason(f)}."
        return short, lg

    def panel_g(self, f: dict) -> tuple[str, str]:
        m, t = self.m, self.t
        R, S, D = t.R, t.S, t.D
        pairs = [(r, d) for r, d in m.rows if d is not None and math.isfinite(r.mean_post)
                 and math.isfinite(d.mean_post)]
        up = sum(r.mean_post > d.mean_post for r, d in pairs)
        nup = self.num("G_pairs_higher", up, "int", f"{m.target.file}:post")
        npair = self.num("G_pairs", len(pairs), "int", f"{m.decoy.file}:seqid")
        pm = self.stat("posterior_mean", "paint", "all", "posterior_mean", "f3")
        dm = self.stat("decoy_posterior_mean", "paint", "all", "decoy_posterior_mean", "f3")
        ex = self.stat("excess", "paint", "all", "excess", "sign3")
        nq = self.stat("n_tracts_q05", "paint", "all", "n_tracts_q05", "int")
        nt = self.stat("n_tracts", "paint", "all", "n_tracts", "int")
        q = [x["q"] for x in m.target.tracts if math.isfinite(x["q"])]
        qmin = self.num("G_min_q", min(q), "f2", f"{m.target.tract_file}:q") if q else None
        lg = (f"Mean posterior of {D}-like ancestry over each drawn sequence's painted insertions, {S} (decoy, "
              f"open) and {R} (filled), joined per pair of panel {self.L['A']}")
        if pm and dm:
            lg += f"; bars: genome-wide means ({dm} and {pm}" + (f"; excess {ex})" if ex else ")")
        lg += "."
        if nup and npair:
            lg += f" {R} is higher on {nup} of {npair} pairs."
        if nq and nt:
            lg += (f" {nq} of its {nt} tracts reach q ≤ {Q_TRACT:g} against the decoy"
                   + (f" (smallest q {qmin})." if qmin else "."))
        exv = _num(m.run.row("paint", t.scope, "all", "excess"))
        if pm and dm and math.isfinite(exv):
            if exv > 0:
                short = (f"{R}'s mean posterior is above its sister's control painting ({pm} against {dm})"
                         + (f", higher on {nup} of {npair} {m.text.get('unit', 'sequence')} pairs." if nup and npair
                            else "."))
            else:
                short = f"{R}'s painting is not more donor-like than its sister's control painting ({pm} against {dm})."
        else:
            short = "No genome-wide comparison with the control painting is available."
        return short, lg

    def panel_h(self, f: dict) -> tuple[str, str]:
        m = self.m
        parts, up, young = [], 0, 0
        for tr, rs in _h_rows(m):
            d = m.run.row("D_pav", tr.scope, "all", "D")
            key = "H_" + re.sub(r"[^A-Za-z0-9]+", "_", tr.scope).strip("_")
            v = self.num(key, _num(d), "sign2", _src("D_pav", tr.scope, "all", "D")) if d is not None else None
            parts.append(f"{tr.scope} {v}" if v else f"{tr.scope} not available")
            if d is not None and _direction(_num(d), _num(d, "z")) == "up":
                up += 1
                g = m.run.row("D_pav_gradient", tr.scope, "K=0 minus K>=0.01", "dD")
                young += g is not None and _direction(_num(g), _num(g, "z")) == "up"
        n = len(m.trios)
        nt = self.num("trios_n", n, "int", "stats.tsv:roles")
        nu = self.num("trios_excess", up, "int", "stats.tsv:D_pav/<trio>/all/D")
        ny = self.num("trios_young", young, "int", "stats.tsv:D_pav_gradient/<trio>/K=0 minus K>=0.01/dD")
        lg = (f"Every trio of the run (R,S;D), one row each, sorted by D_pav: its 95% interval and point (open: not "
              "testable by the block jackknife), and D among the youngest (K = 0, dark) and oldest (K ≥ 0.01, light) "
              f"insertions, as in panel {self.L['F']}, on one shared axis; the lead trio's label is bold, the "
              f"trio with the largest |Z| unless the user named it. In {nu} of the {nt} trios R carries significantly "
              f"more donor-shared insertions than S (Z ≥ {Z_SIG:g}); in {ny} of those the excess is concentrated in "
              "young insertions. By trio: " + "; ".join(parts) + ".")
        short = (f"Every trio of the run, the same test: in {nu} of the {nt} the recipient carries significantly more "
                 "insertions shared with the donor than its sister does"
                 + (f", most of all young ones in {ny} of them." if up else "."))
        return short, lg

    def caveats(self) -> str:
        m = self.m
        out = []
        if m.painted:
            sc = self.stat("shortest_callable", "paint", "all", "shortest_callable_loci_lo", "int")
            sp = self.stat("spacing_median", "paint", "all", "spacing_median", "bp")
            if sc and sp:
                one = "a single insertion" if sc == "1" else f"{sc} insertions"
                out.append(f"the painting calls a tract from as few as {one} (at the median spacing of informative "
                           f"insertions, {sp}), so its q-value, not the number of tracts, carries the evidence")
            ninv = self.stat("n_inv_loci", "paint", "all", "n_inv_loci", "int")
            if ninv:
                out.append(f"{ninv} painted insertions lie inside inverted blocks and stay in the model")
        out.append("K = 0 is not age 0, and gene conversion between the two LTRs makes some copies look young")
        if m.example is None and m.example_note:
            out.append(m.example_note)
        out += self.warnings()
        return "Caveats: " + "; ".join(out) + "."

    def warnings(self) -> list[str]:
        """The run's findings warnings in plain words: a resolution-symmetry warning on a drawn comparison as one
        sentence naming its panel; every other warning only counted (no stage, file or metric names)."""
        t, L = self.t, self.L
        drawn = {("Trio tests", "resolution", t.scope): ("F", "the donor-sharing test"),
                 ("Painting", "paint_resolution", t.scope): ("A", "the painting"),
                 ("Relationships", "kksc_resolution", None): ("E", "the sharing counts read from the donor's calls")}
        out, other = [], 0
        for section, body in _findings_warnings(self.run.outdir):
            mm = re.match(r"(\S+) (\S+): resolved fractions differ by more than [\d.]+ \(([^)]*)\)", body)
            hit = None
            if mm:
                test, scope, _ = mm.groups()
                trio_set = set(scope.split(",")) == {t.R, t.S, t.D}
                hit = drawn.get((section, test, scope)) or (drawn.get((section, test, None)) if trio_set else None)
            fr = _fractions(mm.group(3)) if mm else None
            if hit is None or hit[0] not in L or fr is None:
                other += 1
                continue
            text = self.text(f"warning_{len(out) + 1}", fr, "findings.txt:Warnings")
            out.append(f"calls are resolved unequally among the genomes ({text}), which favors the better-resolved "
                       f"genome in {hit[1]} (panel {L[hit[0]]})")
        if other:
            k = self.num("warnings_other", other, "int", "findings.txt:Warnings")
            out.append(f"{k} further warning{'s are' if other > 1 else ' is'} listed with the run's findings")
        return out


def _panels(letters: list[str]) -> str:
    """"panel A" or "panels A and B"."""
    return ("panel " if len(letters) == 1 else "panels ") + _and(letters)


def _and(items: list[str]) -> str:
    """"a", "a and b", "a, b and c"."""
    return items[0] if len(items) == 1 else ", ".join(items[:-1]) + " and " + items[-1] if items else ""


def _ratio_words(share_d: float, share_r: float) -> str:
    """How often the control reaches a posterior compared with the recipient, from the two printed shares."""
    if not share_r:
        return "though the recipient does not"
    r = share_d / share_r
    if r < 0.15:
        return "far less often"
    if r < 0.4:
        return "less than half as often"
    if r < 0.6:
        return "about half as often"
    if r < 0.85:
        return "somewhat less often"
    if r <= 1.15:
        return "about as often"
    return "more often"


def _share_high(post: np.ndarray) -> str:
    """The share of painted insertions with a posterior of 0.9 or more, as "10.2%"."""
    v = post[np.isfinite(post)]
    return f"{100 * float(np.mean(v >= 0.9)):.1f}%" if len(v) else "none"


def _analysis_line(outdir: Path) -> str:
    """The long legend's last line: where the analysis ran, the script's source path, the date."""
    script = os.path.abspath(inspect.getsourcefile(_analysis_line) or __file__)
    return (f"Analysis: {socket.gethostname()}:{os.path.abspath(outdir)} · Script: {script} "
            f"· {datetime.date.today().isoformat()}")


def _findings_file(outdir: Path) -> list[str]:
    """The lines of findings.txt; none without the file."""
    path = outdir / "findings.txt"
    return path.read_text(encoding="utf-8").splitlines() if path.is_file() else []


def _findings_lines(outdir: Path) -> list[tuple[str, str]]:
    """findings.txt's lines as (section, text), the Warnings section excluded."""
    out, section = [], ""
    for line in _findings_file(outdir):
        if line == "Warnings":
            break
        if line and not line.startswith("- "):
            section = line
        elif line.startswith("- "):
            out.append((section, line[2:]))
    return out


def _findings_warnings(outdir: Path) -> list[tuple[str, str]]:
    """findings.txt's Warnings section as (section, text) pairs ("- <section>: WARNING: <text>")."""
    lines = _findings_file(outdir)
    if "Warnings" not in lines:
        return []
    out = []
    for line in lines[lines.index("Warnings") + 1:]:
        if line.startswith("- ") and line != "- none":
            section, sep, text = line[2:].partition(": WARNING: ")
            out.append((section, text) if sep else ("", line[2:]))
    return out


def _fractions(spread: str) -> str | None:
    """"A 0.60, B 0.58" -> "A 60%, B 58%" (None if it does not parse)."""
    parts = []
    for item in spread.split(","):
        name, _, v = item.strip().rpartition(" ")
        try:
            parts.append(f"{name} {100 * float(v):.0f}%")
        except ValueError:
            return None
    return ", ".join(parts) if parts else None


def _legends_no_trio(m: _Model, nums: style.Numbers) -> tuple[list, list]:
    cat = m.run.cat
    nt = m.no_trio
    L = m.letters
    head = (f"With {nums.get('no_trio')} ({', '.join(cat.ingroup)}) there is no trio, so no donor-sharing test or "
            "painting: the figure shows what each genome holds at the other's insertions.")
    short = [[(head, "")]]
    intro = (f"The run has {len(cat.labels)} genomes ({', '.join(_role(m, x) for x in cat.labels)}), "
             f"{len(cat.ingroup)} of them ingroup: fewer than the three a trio needs. States as genotyped: "
             "present (intact, solo LTR or truncated), empty site (the target site with one TSD copy), "
             "deletion, other sequence, assembly gap or unresolved.")
    long_ = [[(intro, "")]]
    if "A" in L:
        X, Y = nt["X"], nt["Y"]
        a = nums.add("A_present", _int(nt["n_present"]), f"pairs/{X}{PAIR_SEP}{Y}.calls.tsv.gz:state")
        b = nums.add("A_empty", _int(nt["n_empty"]), f"pairs/{X}{PAIR_SEP}{Y}.calls.tsv.gz:state")
        nd, nnd = nums.get("sequences_drawn"), nums.get("sequences_not_drawn")
        short.append([(L["A"], "b"), (f" {X}'s insertions along its {m.text['unit']}s, by what {Y} holds there.",
                                       "")])
        text_a = (f" {X}'s insertions by their state in {Y}: {X}'s {nd} longest sequences ({nnd} shorter ones not "
                  f"drawn); every insertion {X} annotated is a tick: blue below where {Y} carries it ({a}), vermillion "
                  f"above where {Y} has the empty site ({b}); hatched: display bins with neither.")
        long_.append([(L["A"], "b"), (text_a, "")])
    parts = []
    for a_, b_, comp in nt["pairs"]:
        tot = sum(comp.values())
        e = comp.get("E", 0)
        share = (nums.add(f"B_{a_}_{b_}_E", _pct(e / tot), f"pairs/{a_}{PAIR_SEP}{b_}.calls.tsv.gz:state")
                 if tot else None)
        parts.append(f"{a_} in {b_}: {share} empty sites of {_int(tot)}" if share else f"{a_} in {b_}: no insertions")
    short.append([(L["B"], "b"), (" What each genome holds at the other genomes' insertions.", "")])
    text_b = (" Share of each annotated genome's insertions by their state in each other genome (that genome's "
              "own calls; a state with no call is not drawn): " + "; ".join(parts) + ".")
    long_.append([(L["B"], "b"), (text_b, "")])
    long_.append([(_analysis_line(m.run.outdir), "")])
    return short, long_


# ---------------------------------------------------------------- the library screen

NEEDS = "needs assemblies: run geneflow without --screen"  # a placeholder panel's text (module docstring)
# what each placeholder stands for: the assembly route's panel (its box's first line, and the legends)
NEEDED = {"C": "one insertion at base-pair scale", "D": "absence classes by age", "E": "shared-insertion counts",
          "F": "D_pav with its strata", "H": "every trio's D_pav"}
_INFERRED = 99  # panel B's code of an inferred absence (A): drawn open, never keyed as an empty site
# panel B's code of a screen's call by its class (calls.states.call_class with the screen's inferred absences)
_SCREEN_DRAW = {"present": STATE_CODE["P"], "absent": _INFERRED, "missing": STATE_CODE["U"]}


def _paint_row(m: _Model, statistic: str):
    """The lead trio's painting row `statistic`: test paint (assembly), or the screen's hmm_<SCREEN_ALPHABET>."""
    if m.screen:
        return m.run.row(f"hmm_{SCREEN_ALPHABET}", m.trio.scope, m.trio.R, statistic)
    return m.run.row("paint", m.trio.scope, "all", statistic)


def _screen_run(outdir: Path) -> bool:
    """Whether the run in `outdir` is a library screen: its manifest.json says route library (a run without one is a
    UsageError)."""
    return read_manifest(outdir, REMEDY).get("route") == "library"


class _ScreenCat:
    """What the overview's shared drawing code reads of a run's catalog, from a library screen's files: the sheet's
    samples, groups and roles; each genome's sequence names in natural order (library_stats.chrom_key) and their
    lengths, the last end of the genome's elements on each (work/elements); and `proj`, filled by
    `_screen_painting`: the recipient's state in each genome at its painted elements."""

    def __init__(self, outdir: Path):
        self.outdir = outdir
        self.samples = read_samples(_sheet(outdir))
        self.labels = tuple(s.label for s in self.samples)
        self.ingroup = tuple(s.label for s in self.samples if s.role != "outgroup")
        self.outgroups = tuple(s.label for s in self.samples if s.role == "outgroup")
        self.annotated = dict.fromkeys(self.labels, True)
        self.seqids: dict[str, tuple[str, ...]] = {}
        self.lengths: dict[str, np.ndarray] = {}
        for x in self.labels:
            ends: dict[str, int] = {}
            for e in read_elements(outdir / "work" / "elements" / f"{x}.tsv"):
                ends[e.seqid] = max(ends.get(e.seqid, 0), e.end)
            names = sorted(ends, key=chrom_key)
            self.seqids[x] = tuple(names)
            self.lengths[x] = np.array([ends[n] for n in names], dtype=np.int64)
        self.proj: dict[str, np.ndarray] = {}

    def g(self, x: str) -> int:
        return self.labels.index(x)


def _screen_overview(outdir: Path, trio: tuple[str, str, str] | None) -> list[Path]:
    """`overview` of a library screen (module docstring)."""
    import matplotlib as mpl

    run = _Run(outdir, _ScreenCat(outdir))
    nums = style.Numbers()
    style.setup_fonts()
    mono = style.rc_params()["mathtext.tt"]
    m = _screen_model(run, nums, trio)
    with mpl.rc_context(style.rc_params()):
        fig = _draw(m, mono)
        short, long_ = _screen_legends(m, nums)
        style.save(fig, outdir / STEM, nums, legend=[short, long_])
    log.info("overview: %s.png, %s.pdf, %s_numbers.tsv", STEM, STEM, STEM)
    return style.figure_files(outdir / STEM)


def _screen_painting(run: _Run, trio: _Trio, decoy: bool) -> _Painting | None:
    """The lead trio's SCREEN_ALPHABET painting of R (decoy: of S) from the screen's posterior and tract files: an
    element of a donor-like class is RD, of a sister-like class RS, any other MISSING; None without a row. For R it
    fills the catalog stand-in's `proj` from R's pair calls, each drawn by its class (calls.states.call_class with
    the screen's inferred absences: present, absent inferred, missing; R itself present)."""
    cat = run.cat
    X, Y = (trio.S, trio.R) if decoy else (trio.R, trio.S)
    stem = "decoy_" if decoy else ""
    name, tname = f"{stem}posterior.{X}.tsv", f"{stem}tracts.{X}.bed"
    rows = [r for r in run.table(name) or [] if r.get("sister") == Y and r.get("donor") == trio.D]
    if not rows:
        return None
    a = ALPHABETS[SCREEN_ALPHABET]
    code = {c: RD if c in a.donor else RS for c in a.emits}
    order = {n: j for j, n in enumerate(cat.seqids[X])}
    cls = np.array([code.get(r["class"], MISSING) for r in rows], dtype=np.int8)
    post = np.array([_float(r[f"post_{SCREEN_ALPHABET}"]) for r in rows])
    seq = np.array([order[r["seqid"]] for r in rows], dtype=np.int64)
    start = np.array([int(r["start"]) for r in rows], dtype=np.int64)
    end = np.array([int(r["end"]) for r in rows], dtype=np.int64)
    eid = np.array([r["eid"] for r in rows])
    if not decoy:
        at = {e: i for i, e in enumerate(eid.tolist())}
        proj = np.zeros((len(rows), len(cat.labels)), dtype=np.uint8)
        proj[:, cat.g(X)] = STATE_CODE["P"]
        state = {s: _SCREEN_DRAW[call_class((s,), absent=INFERRED_ABSENT)] for s in PAIR_STATES}
        for g in cat.labels:
            path = run.outdir / "pairs" / f"{X}{PAIR_SEP}{g}.calls.tsv.gz"
            if g != X and path.is_file():
                for _, (e, st) in stream_calls(path, ("eid", "state")):
                    if e in at:
                        proj[at[e], cat.g(g)] = state.get(st, 0)
        cat.proj[X] = proj
    o = np.lexsort((start, seq))
    tracts = [{"seq": order[r["chrom"]], "seqid": r["chrom"], "start": int(r["start0"]) + 1, "end": int(r["end"]),
               "n_loci": int(r["n_loci"]), "n_RD": int(r["n_RD"]), "n_RS": int(r["n_RS"]),
               "mean_post": _float(r["mean_post"]), "q": _float(r.get("q", ".")),
               "youngest_k": _float(r.get("youngest_k", "."))}
              for r in run.table(tname) or []
              if r.get("sister") == Y and r.get("donor") == trio.D and r.get("alphabet") == SCREEN_ALPHABET]
    return _Painting(X, eid[o], np.arange(len(rows))[o], seq[o], start[o], end[o], cls[o], post[o], tracts, name,
                     tname)


def _screen_partners(run: _Run, X: str, Y: str, drawn: list[int]) -> dict[int, int | None]:
    """Each drawn X sequence's counterpart in Y: the Y sequence holding most of its 1:1 orthologs (X's P calls in
    Y), each Y sequence used once (X order; ties by Y's order); None when none is left."""
    cat = run.cat
    ox, oy = ({n: j for j, n in enumerate(cat.seqids[g])} for g in (X, Y))
    counts: dict[tuple[int, int], int] = {}
    path = run.outdir / "pairs" / f"{X}{PAIR_SEP}{Y}.calls.tsv.gz"
    for _, (st, a, b) in stream_calls(path, ("state", "a_seqid", "b_seqid")):
        if st == "P" and a in ox and b in oy:
            counts[ox[a], oy[b]] = counts.get((ox[a], oy[b]), 0) + 1
    used: set[int] = set()
    out: dict[int, int | None] = {}
    for s in drawn:
        best = sorted((-c, t) for (s2, t), c in counts.items() if s2 == s and t not in used)
        out[s] = best[0][1] if best else None
        if out[s] is not None:
            used.add(out[s])
    return out


def _screen_model(run: _Run, nums: style.Numbers, trio: tuple[str, str, str] | None) -> _Model:
    """The model of a screen's overview (module docstring)."""
    cat = run.cat
    trios = _trios(run)
    if trio is not None and _Trio(*trio) not in trios:
        have = ", ".join(t.scope for t in trios) or "none"
        raise AnalysisError(f"--trio {','.join(trio)}: not a trio of this run's stats.tsv (trios: {have})")
    lead = _Trio(*trio) if trio is not None else (trios[0] if trios else None)
    m = _Model(run=run, trio=lead, trios=trios, groups={s.label: s.group for s in cat.samples},
               named=trio is not None, screen=True)
    m.placeholders = {"C", "D", "E", "F", "H"}
    if lead is None:
        nums.add("no_trio", "fewer than three ingroup genomes", "stats.tsv:roles")
        m.panels, m.letters = [], {}
        return m
    t = lead
    m.target = _screen_painting(run, t, decoy=False)
    m.decoy = _screen_painting(run, t, decoy=True) if m.target is not None else None
    m.painted = m.target is not None
    for p in ("C", "D", "E", "F"):
        nums.add(f"panel_{p}", NEEDS, "manifest.json:route library")
    if len(trios) > 1:
        nums.add("panel_H", NEEDS, "manifest.json:route library")
    if not m.painted:
        nums.add("no_painting", f"no {SCREEN_ALPHABET} painting of {t.R} (skipped: findings.txt, Painting)",
                 _src(f"hmm_{SCREEN_ALPHABET}", t.scope, t.R, "pi"))
        m.panels = ["A", "C", "D", "E", "F"] + (["H"] if len(trios) > 1 else [])
        m.placeholders |= {"A"}
        m.letters = {p: "ABCDEFGHIJ"[j] for j, p in enumerate(m.panels)}
        return m
    drawn, m.n_not = _drawn(cat, t.R)
    m.n_drawn = len(drawn)
    part = _screen_partners(run, t.R, t.S, drawn) if m.decoy is not None else {s: None for s in drawn}
    lens = [int(cat.lengths[t.R][s]) for s in drawn] + [int(cat.lengths[t.S][y]) for y in part.values()
                                                        if y is not None]
    m.xmax = max(lens) if lens else 1
    m.bin_bp = _nice_floor(m.xmax / BINS)
    m.text["inv"] = tally = {}
    inv_r, inv_s = _inverted(run, t.R, {t.S, t.D}, m.bin_bp, tally), _inverted(run, t.S, {t.R, t.D}, m.bin_bp, tally)
    rl = _short_names([cat.seqids[t.R][s] for s in drawn])
    pl = iter(_short_names([cat.seqids[t.S][y] for y in part.values() if y is not None]))
    for s, lab in zip(drawn, rl, strict=True):
        r = _row(m.target, cat, s, f"{t.R} {lab}", m.bin_bp, inv_r, weighted=True)
        y = part[s]
        d = _row(m.decoy, cat, y, f"{t.S} {next(pl)}", m.bin_bp, inv_s, weighted=True) if y is not None else None
        m.rows.append((r, d))
    shown = {t.R: drawn} | ({t.S: [y for y in part.values() if y is not None]} if m.decoy is not None else {})
    m.text["unit"] = "chromosome" if _chromosome_scale(cat, shown) else "sequence"
    nums.add("sequences_drawn", m.n_drawn, f"work/elements/{t.R}.tsv:seqid")
    nums.add("sequences_not_drawn", m.n_not, f"work/elements/{t.R}.tsv:seqid")
    order = _typical(m, m.target.tracts)
    if order:
        x = m.target.tracts[order[0]]
        m.window = _window(m, x, x)
    else:
        nums.add("no_tract", "no tract called", f"{m.target.tract_file}:n_loci")
    m.panels = (["A"] + (["B"] if m.window else []) + (["G"] if m.decoy is not None else [])
                + ["C", "D", "E", "F"] + (["H"] if len(trios) > 1 else []))
    m.letters = {p: "ABCDEFGHIJ"[j] for j, p in enumerate(m.panels)}
    return m


def _screen_row(page: _Page, m: _Model, y0: float, h: float) -> None:
    """A screen's last row: panel G (the posterior per sequence pair) at the left, then the placeholders in a grid of
    up to three columns, in reading order."""
    boxes = [p for p in m.panels if p in m.placeholders and p != "A"]
    x = 0.0
    if "G" in m.panels:
        _panel_g(page, m, 0.0, y0, 56.0, h)
        x = 59.0
    ncol = min(3, max(1, len(boxes))) if "G" not in m.panels else 2 if len(boxes) <= 4 else 3
    nrow = max(1, -(-len(boxes) // ncol))
    bw, bh = (WIDTH_MM - x - 2.0 * (ncol - 1)) / ncol, (h - 2.0 * (nrow - 1)) / nrow
    for k, p in enumerate(boxes):
        r, c = divmod(k, ncol)
        _placeholder(page, m, p, x + c * (bw + 2.0), y0 + r * (bh + 2.0), bw, bh)


def _placeholder(page: _Page, m: _Model, p: str, x0: float, y0: float, w: float, h: float) -> None:
    """A panel whose data the run lacks: its letter, a light box and why (NEEDS; with no painting, why not)."""
    from matplotlib.patches import Rectangle

    page.letter(x0 + 0.5, y0, m.letters[p])
    ax = page.ax(x0 + 6.0, y0 + 1.0, max(w - 7.0, 1.0), max(h - 3.0, 1.0))
    _bare(ax)
    ax.add_patch(Rectangle((0, 0), 1, 1, transform=ax.transAxes, fc="white", ec=LIGHT, lw=0.6))
    if p == "A":
        text = f"no {SCREEN_ALPHABET} painting of {m.trio.R}: skipped (findings.txt)"
    else:
        text = f"{NEEDED[p]}\n{NEEDS}"
    ax.text(0.5, 0.5, text, transform=ax.transAxes, ha="center", va="center", fontsize=FS_SMALL, color=GRAY,
            linespacing=1.5)


def _screen_legends(m: _Model, nums: style.Numbers) -> tuple[list, list]:
    """(short, long) legend paragraphs of a screen's overview (runs of (text, style)), panels in their order."""
    L, t = m.letters, m.trio
    run = m.run
    end = [[(_analysis_line(run.outdir), "")]]
    if t is None:
        text = (f"With {nums.get('no_trio')} there is no trio, so no painting: the library screen's statistics are "
                "in stats.tsv and findings.txt.")
        return [[(text, "")]], [[(text, "")], *end]
    src = run.row("roles", t.scope, "all", "source")
    how = {"data": "roles from the single-copy k-mers", "--trio": "roles as set", "sheet order":
           "roles left to the sheet's order (a tie)"}.get(str(src.value) if src else "", "roles as set")
    lead = (f"Library screen of {t.R} (recipient), {t.S} (sister) and {t.D} (donor), {how}, from their LTR-RT "
            "libraries alone: the recent shared ancestry painting.")
    boxes = [p for p in m.panels if p in m.placeholders and p != "A"]
    needs = (f"On genome assemblies {_panels([L[p] for p in boxes])} show "
             + _and([NEEDED[p] for p in boxes])
             + "; they need observed empty sites or the genome sequences, which a library screen does not have (run "
               "geneflow without --screen on the assemblies).")
    if not m.painted:
        head = f"{lead} There is {nums.get('no_painting')}, so no posterior or tract is drawn."
        return [[(head, "")]], [[(head, "")], [(needs, "")], *end]

    def stat(key: str, statistic: str, style_: str) -> str | None:
        row = _paint_row(m, statistic)
        return None if row is None else nums.add(key, _fmt(row.value, style_) or "nan",
                                                 _src(f"hmm_{SCREEN_ALPHABET}", t.scope, t.R, statistic))
    pm, dm = stat("posterior_mean", "posterior_mean", "f3"), stat("decoy_posterior_mean", "decoy_posterior_mean", "f3")
    nt, nq = stat("n_tracts", "n_tracts", "int"), stat("n_tracts_q05", "n_tracts_q05", "int")
    head = (f"{lead} Mean posterior of {t.D}-like ancestry {pm} in {t.R} against {dm} in the decoy {t.S} painted "
            f"with {t.R}'s model; {nt} tracts, {nq} at q \u2264 {Q_TRACT:g} against the decoy's.")
    nd, nn = nums.get("sequences_drawn"), nums.get("sequences_not_drawn")
    texts: dict[str, tuple[str, str]] = {
        "A": (f" {t.R}'s {m.text.get('unit', 'sequence')}s painted, each over its sister counterpart.",
              (f" {t.R}'s {nd} longest sequences ({nn} shorter not drawn), each above the sister's sequence holding "
               f"most of its 1:1 orthologs, painted the same way. Ticks: elements whose ortholog in {t.D} only (above) "
               f"or in the other sister only (below: {t.S} on {t.R}'s rows, {t.R} on {t.S}'s) is near-identical "
               f"(below the 5th percentile of {t.S}-{t.D} ortholog divergence); the bar's color is the posterior of "
               "D-like ancestry, hatched where a bin holds no such element; dark lines are tracts."))}
    if "B" in L:
        w = m.window
        coord = nums.add("B_window", f"{w['name']}:{w['start']:,}-{w['end']:,}", f"{m.target.tract_file}:chrom")
        texts["B"] = (" One typical tract and the genomes' states along it.",
                      (f" The typical tract ({coord}: the tract and one tract length each side), one row per genome "
                       f"with its state from {t.R}'s pair calls: present, absent (inferred: inside a collinear window "
                       "without an ortholog, never an observed empty site) or unresolved; the posterior above, beside "
                       "every element's posterior genome-wide against the decoy's."))
    if "G" in L:
        texts["G"] = (" Mean posterior per sequence pair, decoy to recipient.",
                      (f" The mean posterior per drawn sequence pair, weighted by element length, the decoy {t.S} to "
                       f"the recipient {t.R}; bars are the genome-wide means."))
    short, long_ = [[(head, "")]], [[(head, "")]]
    for p in m.panels:
        if p in texts:
            short.append([(L[p], "b"), (texts[p][0], "")])
            long_.append([(L[p], "b"), (texts[p][1], "")])
        elif p in boxes:
            short.append([(L[p], "b"), (f" {NEEDED[p].capitalize()}: needs assemblies.", "")])
    long_.append([(needs, "")])
    long_.append([(("The screen's own statistics, their false-signal rates and the gate that holds back those uneven "
                    "annotation can fake are in findings.txt and stats.tsv."), "")])
    return short, [*long_, *end]


def _draw_text(m: _Model, text: str):
    """A one-line page: what the run has no figure for."""
    page = _Page(20.0)
    page.text(WIDTH_MM / 2, 10.0, text, ha="center", va="center", color=GRAY)
    return page.fig
