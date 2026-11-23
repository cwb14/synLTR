"""Painting of a recipient's genome and the youngest insertion per tract, on either evidence route.

Inputs. The assembly route paints the loci of a run's catalog (catalog.Catalog); the library route paints the records
of a finished library screen (library_stats.LibraryRun: its pair calls, every record of the recipient in library
order). Both go through one HMM core (`forward_backward`, `fit_hmm`, `call_tracts`, `target_decoy`).

Route properties (route.RouteProps): `paint` honors each painting property, and refuses a value this engine does not
implement:
  alphabet        the painting's classes (ALPHABETS: each alphabet's classes, those that emit and the donor-like
                  ones); one painting per alphabet of the route
  decoy_scale     "informative": scale = the target's emitting units over the decoy's, no decoy without one;
                  "all": every unit of each, the decoy always painted
  donor_state     the classes whose fitted state-1 emission must exceed state 0's (and fall below it): DONOR_STATE
  tract_position  "site_start" (a unit sits at its start) or "midpoint" ((start + end) // 2): the HMM's position
                  and the inverted-block flags'
  tract_end       "largest_end" (a tract ends at its units' largest end) or "last_end" (its last unit's end)
  posterior_mean  "unweighted" (each unit once) or "length_weighted" (by end - start + 1): `posterior_mean`
The assembly route's are described first below: the alphabet "insertion", the decoy scale "informative", the donor
state "rd_up", loci at their site start, tracts ending at their loci's largest end, unweighted posterior means. Then
the library route's (section "The library route").

Classes (`classes`, PAINT_CLASSES; alphabet "insertion"). Along the recipient X's sequences, X's own loci, each
classed from X's own calls (`Catalog.pcls`, never loci.tsv's merged states), checked in this order: `old` (some
outgroup present: older than the ingroup split, set aside; a missing outgroup keeps the locus), `missing` (the
sister Y or the donor D unresolved), then `RD` (D present, Y E: donor-like), `RS` (Y present, D E:
sister-like), `both`, `neither`. Only RD and RS emit (codes 0 and 1); the rest are missing data, since emitted
private classes become fake states.

HMM (`forward_backward`, `fit_hmm`, `call_tracts`, `target_decoy`). Two states, 0 sister-like and 1
donor-like; one chain per X sequence (fai index), each starting at the stationary distribution (1 - pi, pi);
between loci d bp apart the state is redrawn from it with probability 1 - q, q = exp(-d / L), so p01 = pi (1 -
q) and p10 = (1 - pi)(1 - q); a missing emission is 1, so the chain is stationary across missing loci.
Positions are X's site starts (1-based; the insertion point), in idx order, so they never decrease along a
chain. `fit_hmm`: start at pi 0.1, L 1e5; each EM
round takes the L of the HMM_L_GRID that maximizes the log-likelihood, then e = (expected counts + 0.5),
normalized, and pi = the mean posterior of the emitting loci clipped to [PI_MIN, PI_MAX]; it stops when the
log-likelihood moves by < HMM_TOL of itself, or after HMM_ITERS rounds. Initial emissions: e0 = the
decoy's class frequencies + 1 (the target's own without a decoy), e1 = e0 * HMM_INIT (5, 0.2) so that state 1
starts RD-enriched and RS-depleted, both normalized. The fitted
state 1 is donor-like when e1[RD] > e0[RD] (`donor_state`); without it no tract is called.

Emissions: EM fits e0 and e1 (`fit_hmm`), as above. On sparse, class-segregated loci EM makes state 0 pure RS
and state 1 pure RD; the target-decoy q-values below are what separate that from real ancestry.

Lengths: L, the switch scale, is reported as `mean_tract_bp`;
the mean donor-like run under the model, L / (1 - pi) (the rate out of state 1 is (1 - pi) / L), as
`mean_donor_run_bp`.

Tracts (`call_tracts`): maximal runs of loci on one sequence with posterior > 0.5 and posterior odds >=
BF_TRACT x prior odds pi / (1 - pi). A tract runs from the start of its first locus to the largest end of
its loci (the last one's, unless a nested host reaches further), so the interval holds every locus. Its
score, sum_post, is the exactly rounded sum of its posteriors (math.fsum): the same posteriors give the same
score wherever they lie and in any order, so a true tie between a target and a decoy tract compares equal.

Decoy and q: the sister Y's loci (those Y discovers), classed by Y's own
calls with X and Y swapped, are painted with the fitted parameters unchanged, and their tracts called by
the same rule. Scores are each tract's sum of posteriors. `target_decoy` (q with the +1
correction): raw = ((#decoy >= x) * scale + 1) / (#target >= x), q = the smallest raw at or below the
tract's score, capped at 1; p = (1 + #decoy >= x) / (1 + #decoy). scale = n informative(X) / n informative(Y)
(informative loci, not all elements: decoy scale "informative"). There is no decoy when Y has no annotation
or no informative locus: q and p are written "none", with a WARNING.

Resolution: the spacing of consecutive informative loci on one sequence, and `shortest_callable` at its
median: the fewest consecutive RD loci whose middle locus passes the tract rule, alone (missing flanks change
nothing: the chain is stationary across them) and between CALLABLE_FLANK RS loci on each side (a heuristic
flank).

Per tract: youngest_k = the smallest K2P of X's copy (the annotation's) over its RD loci: the admixture is
younger than any insertion it transferred. youngest_u95 needs identity.tsv's self rows (`self_counts`): the
smallest u95(diffs, sites) (`seq.divergence.garwood_u95`) over its RD loci, from X's copy, else from D's copy of the same
insertion; "." without them (no K-based stand-in). n_identical = RD loci with K = 0; inv_frac = the share of
its loci inside an inverted block of X against Y or D (`inversions.inside`, keyed on X's fai indexes; the
posterior files flag each locus). Inverted blocks are flags only: the fit ignores them.

Files (`write_paintings`): tracts.<X>.bed and posterior.<X>.tsv per recipient, decoy_tracts.<Y>.bed and
decoy_posterior.<Y>.tsv per sister with a decoy, every trio in one file with its sister and donor columns
(in a decoy file, the sister is X), rows by trio scope then position. Sequence names come from
`Catalog.seqids` (output only). A posterior row's `tract` is its tract as seqid:start-end (1-based), "."
outside one; a decoy tract has q and p_decoy ".". A skipped painting writes nothing, so a recipient whose
every trio is skipped gets no file.

The library route (alphabets "ancestry" and "recent", LTRmix's painting).
Classes (`classes`, LIBRARY_CLASSES), of every record of the painted genome X with sister Y and donor D, from X's pair
calls into D and into Y (library_stats.pair_calls; statuses by the absences ABSENCES["painting"] accepts: P present,
A absent, inferred): unresolved unless resolved (P or A) in both; then near_donor (X's ortholog in D diverges below
the near threshold, its ortholog in Y does not), near_sister (the mirror), near_both, and among the rest donor_only (P
in D, not in Y), sister_only, shared (P in both) and private (in neither). A record whose name holds no coordinates
has no call: unresolved. The near threshold (`near_threshold`) is the NEAR_QUANTILE quantile of the trio sister S's
ortholog divergences to D (n_diff / n_sites over S's records; 0.0 when none is finite), one threshold for the target
and the decoy, so both genomes' classes mean the same. "ancestry" emits near_donor, donor_only, near_sister and
sister_only (codes 0-3; the first two donor-like); "recent" emits near_donor and near_sister. Every other class is
missing data: shared, private and near_both are exchangeable between sisters without gene flow.
Painted set and positions: every record of X, records whose names hold no coordinates included (each its own
sequence, from 0 to its length: a one-record chain), in chain order: its sequence name in text order, then its
midpoint (stable); chains break where the name changes. A record sits at its midpoint (start + end) // 2.
Initial emissions: e0 = the decoy's emitting classes over all its records + 1, e1 = e0 x HMM_INIT (5 for the
donor-like classes, 0.2 for the others), both normalized. Donor state: every near_donor emission of state 1 above
state 0's and every near_sister one below. Decoy: every record of the trio's sister Y, classed by its own calls with X
and Y swapped (the same threshold), painted with the fitted parameters; scale = records of X / records of Y. A tract
runs from its first element's start to its last element's end. Posterior means weight each record by its length
(end - start + 1), over all records, the decoy's alike.
Skips and NA as on the assembly route: a painting is skipped, with its reason, when X has no element of a donor-like
class or none of a sister-like one; a collapsed fit (pi at its floor, no tract) prints NA for L, the donor state and
the emissions.
Rows (`library_rows`; layer library, test hmm_<alphabet>, scope the trio, stratum X): sister, donor, pi,
mean_tract_bp (L), loglik, n_tracts, n_decoy_tracts, n_tracts_q05, posterior_mean, decoy_posterior_mean,
excess_over_decoy, near_threshold, e1_<class> and e0_<class> per emitted class in the alphabet's order, donor_state.
Files (`write_paintings`): one set per genome holding every trio and both alphabets: posterior.<X>.tsv
(LIBRARY_POSTERIOR_COLUMNS: eid, seqid, start, end, sister, donor, the class, one posterior per alphabet and the
inverted-block flag) per element of X, records whose names hold no coordinates left out (no element: no row);
tracts.<X>.bed (LIBRARY_TRACT_COLUMNS: TRACT_COLUMNS with the alphabet last; name <D>-<alphabet>; n_RD and n_RS the
tract's elements of a donor-like and of a sister-like class; youngest_k and n_identical from X's annotated K2P;
youngest_u95 "."); decoy_tracts.<Y>.bed and decoy_posterior.<Y>.tsv alike. Rows by trio scope, then alphabet, then
position. Floats are written exactly (the shortest text that reads back as the same number), so every number LTRmix
printed is recoverable at its precision.
"""
from __future__ import annotations

import fnmatch
import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from numba import njit

from synltr.analyses.geneflow.catalog import AnalysisError, Catalog, age_bin
from synltr.analyses.geneflow.inversions import InvBlock, inside
from synltr.analyses.geneflow.library_stats import ABSENCES, LibraryRun, pair_calls
from synltr.analyses.geneflow.route import RouteProps, require
from synltr.analyses.geneflow.trio import Trio
from synltr.base.constants import (
    BF_TRACT,
    CALLABLE_FLANK,
    CALLABLE_MAX,
    EMIT_PSEUDO,
    HMM_INIT,
    HMM_ITERS,
    HMM_L0,
    HMM_L_GRID,
    HMM_PI0,
    HMM_TOL,
    NEAR_QUANTILE,
    PI_MAX,
    PI_MIN,
    Q_TRACT,
    TRACT_POST,
)
from synltr.calls.records import tsv_field
from synltr.seq.divergence import garwood_u95
from synltr.stats.report import Findings, StatRow, num, resolution


@dataclass(frozen=True)
class Alphabet:
    """A painting alphabet: the classes of a painted unit (its class code is the index here), those that emit (their
    emission codes in this order; every other class is missing data) and, of those, the donor-like ones (state 1
    starts enriched for them, HMM_INIT)."""

    classes: tuple[str, ...]
    emits: tuple[str, ...]
    donor: frozenset[str]

    @property
    def codes(self) -> np.ndarray:
        """int64 emission code per class code (-1: missing data)."""
        return np.array([self.emits.index(c) if c in self.emits else -1 for c in self.classes], dtype=np.int64)


# the assembly route's classes: emissions RD 0, RS 1, the rest -1 (missing data)
PAINT_CLASSES = ("RD", "RS", "both", "neither", "missing", "old")
# the library route's classes (LTRmix's categories, then a record unresolved in the sister or the donor)
LIBRARY_CLASSES = ("near_donor", "near_sister", "near_both", "donor_only", "sister_only", "shared", "private",
                   "unresolved")
ALPHABETS = {  # the painting alphabets this engine paints (route.RouteProps.alphabet)
    "insertion": Alphabet(PAINT_CLASSES, ("RD", "RS"), frozenset({"RD"})),
    "ancestry": Alphabet(LIBRARY_CLASSES, ("near_donor", "donor_only", "near_sister", "sister_only"),
                         frozenset({"near_donor", "donor_only"})),
    "recent": Alphabet(LIBRARY_CLASSES, ("near_donor", "near_sister"), frozenset({"near_donor"})),
}
CATALOG_ALPHABETS = ("insertion",)  # the alphabets of a run's catalog (assembly loci)
RECORD_ALPHABETS = ("ancestry", "recent")  # the alphabets of a library screen's records
# donor-state check: (classes whose state-1 emission must exceed state 0's, classes whose must fall below it)
DONOR_STATE = {"rd_up": (("RD",), ()), "near_donor_up_near_sister_down": (("near_donor",), ("near_sister",))}
PAINT_LABELS = {"ancestry": "Ancestry painting", "recent": "Recent shared ancestry"}  # library findings
POSTERIOR_COLUMNS = ("lid", "seqid", "start", "end", "sister", "donor", "class", "post", "tract", "inv")
# a genome's painting files, by its label: as the recipient (tracts, posterior) and as a decoy
PAINTING_FILES = ("tracts.{}.bed", "posterior.{}.tsv", "decoy_tracts.{}.bed", "decoy_posterior.{}.tsv")
TRACT_COLUMNS = ("#chrom", "start0", "end", "name", "score", "strand", "sister", "donor", "n_loci", "n_RD",
                 "n_RS", "mean_post", "sum_post", "q", "p_decoy", "youngest_k", "youngest_u95", "n_identical",
                 "inv_frac")
LIBRARY_TRACT_COLUMNS = (*TRACT_COLUMNS, "alphabet")
LIBRARY_LAYER = "library"  # the library route's painting rows
RD, RS, BOTH, NEITHER, MISSING, OLD = range(len(PAINT_CLASSES))
UNRESOLVED = LIBRARY_CLASSES.index("unresolved")
PAINTING, DATING = "Painting", "Dating"  # findings sections (stats.report.SECTIONS)


def library_posterior_columns(alphabets: tuple[str, ...]) -> tuple[str, ...]:
    """The library route's posterior file columns for its alphabets (module docstring)."""
    return ("eid", "seqid", "start", "end", "sister", "donor", "class", *(f"post_{a}" for a in alphabets), "inv")


@dataclass(eq=False)
class HmmFit:
    """A fitted painting HMM: pi (the donor-like fraction), L (the switch scale, bp), emissions e0 (state 0)
    and e1 (state 1) over the alphabet's emitting classes, the posterior of state 1 per locus and the
    log-likelihood."""

    pi: float
    L: float
    e0: np.ndarray
    e1: np.ndarray
    post: np.ndarray
    loglik: float


@dataclass
class Tract:
    """A run of loci called donor-like: seq is X's sequence code (a fai index; on the library route the rank of
    the sequence name in text order), start/end 1-based on X; first/last are positions in the painting's idx
    (inclusive)."""

    seq: int
    start: int
    end: int
    first: int
    last: int
    mean_post: float
    sum_post: float


@dataclass(eq=False)
class Painting:
    """One trio's painting of its recipient R (trio.R) and of the decoy, its sister S (trio.S), under the route
    properties `props`, in one alphabet. idx: R's painted units in chain order (catalog indexes in R's (fai index,
    start) order; on the library route record indexes, module docstring), with their class (the alphabet's class
    index), posterior and inverted-block flag; q and p per tract (None without a decoy). skipped: why nothing was
    painted (fit None, post nan). near: the library route's near threshold (None on the assembly route)."""

    trio: Trio
    idx: np.ndarray
    cls: np.ndarray
    post: np.ndarray
    inv: np.ndarray
    fit: HmmFit | None
    donor_state: bool
    tracts: list[Tract]
    q: list[float | None]
    p: list[float | None]
    decoy_idx: np.ndarray
    decoy_cls: np.ndarray
    decoy_post: np.ndarray
    decoy_tracts: list[Tract]
    skipped: str | None
    # what the files and rows need besides the catalog
    props: RouteProps | None = None
    decoy_inv: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=bool))
    scale: float | None = None                       # target-decoy scale; None: no decoy
    self_counts: dict[tuple[str, str], tuple[int, int]] | None = None
    alphabet: str = "insertion"
    near: float | None = None

    @property
    def collapsed(self) -> bool:
        """The fit collapsed: pi at its floor and no tract. State 1 then holds no data, so its emissions, L
        and the donor-state rule come from the pseudocounts: no donor-like ancestry is detectable."""
        return self.fit is not None and self.fit.pi <= PI_MIN and not self.tracts


# ---------------------------------------------------------------- classes


def classes(table: Catalog | LibraryRun, X: str, Y: str, D: str, *, alphabet: str, near: float | None = None
            ) -> tuple[np.ndarray, np.ndarray]:
    """(idx, cls): X's painted units in chain order and their int8 class (the index in the alphabet's classes;
    module docstring). A catalog (alphabet "insertion"): X's loci as int64 catalog indexes in X's (fai index, start)
    order, ties in catalog order: X's members, every state from X's own calls. A library screen (alphabets
    "ancestry" and "recent", which share LIBRARY_CLASSES): every record of X as int64 record indexes in (sequence
    name in text order, midpoint) order, stable, classed with the near threshold `near` (`near_threshold`). A
    ValueError for an alphabet outside ALPHABETS or another input's."""
    _input_alphabet(table, alphabet)
    if not isinstance(table, Catalog):
        if near is None:
            raise ValueError("the library route's painting classes need the near threshold (near_threshold)")
        return _record_classes(table, X, Y, D, near)
    cat = table
    _check(cat, X, Y, D)
    gx = cat.g(X)
    rows = cat.member[:, gx].copy()
    y, d = cat.pcls(X, Y), cat.pcls(X, D)
    out = [cat.pcls(X, o) for o in cat.outgroups]
    rows &= cat.seq[:, gx] >= 0  # a present site has coordinates
    idx = np.flatnonzero(rows)
    idx = idx[np.lexsort((cat.start[idx, gx], cat.seq[idx, gx]))].astype(np.int64)
    y, d = y[idx], d[idx]
    old = np.zeros(len(idx), dtype=bool)
    for o in out:
        old |= o[idx] == 1
    cls = np.select([old, (y < 0) | (d < 0), (d == 1) & (y == 0), (y == 1) & (d == 0), (y == 1) & (d == 1)],
                    [OLD, MISSING, RD, RS, BOTH], default=NEITHER).astype(np.int8)
    return idx, cls


def _input_alphabet(table: Catalog | LibraryRun, alphabet: str) -> None:
    """A ValueError unless `alphabet` is one this engine paints on `table`'s kind of input."""
    require("painting alphabet", alphabet, tuple(ALPHABETS))
    mine, kind = ((CATALOG_ALPHABETS, "a run's catalog of assembly loci") if isinstance(table, Catalog)
                  else (RECORD_ALPHABETS, "a library screen's records"))
    if alphabet not in mine:
        raise ValueError(f"painting alphabet {alphabet!r} does not paint {kind} (its alphabets: {', '.join(mine)})")


def _record_classes(run: LibraryRun, X: str, Y: str, D: str, near: float) -> tuple[np.ndarray, np.ndarray]:
    """`classes` of a library screen's records (module docstring, the library route)."""
    _check(run, X, Y, D)
    xd, xy = pair_calls(run, X, D), pair_calls(run, X, Y)
    st_d, st_y = xd.status(ABSENCES["painting"]), xy.status(ABSENCES["painting"])
    near_d = np.isfinite(xd.div) & (xd.div < near)
    near_y = np.isfinite(xy.div) & (xy.div < near)
    in_d, in_y = st_d == 2, st_y == 2
    resolved = (st_d > 0) & (st_y > 0)
    rest = resolved & ~near_d & ~near_y
    code = {c: np.int8(i) for i, c in enumerate(LIBRARY_CLASSES)}
    cls = np.select([resolved & near_d & ~near_y, resolved & near_y & ~near_d, resolved & near_y & near_d,
                     rest & in_d & ~in_y, rest & in_y & ~in_d, rest & in_y & in_d, rest],
                    [code[c] for c in LIBRARY_CLASSES[:7]], default=code["unresolved"]).astype(np.int8)
    rec = run.rec(X)
    order = np.lexsort((rec.mid, run.sequences(X)[1])).astype(np.int64)
    return order, cls[order]


def near_threshold(run: LibraryRun, S: str, D: str) -> float:
    """The library route's near threshold of a trio (module docstring): the NEAR_QUANTILE quantile of the finite
    divergences of S's records to their orthologs in D (n_diff / n_sites), 0.0 when none is finite."""
    d = pair_calls(run, S, D).div
    return float(np.nanquantile(d, NEAR_QUANTILE)) if np.isfinite(d).any() else 0.0


def _check(table: Catalog | LibraryRun, X: str, Y: str, D: str) -> None:
    """Three distinct genomes of the run (an AnalysisError names an unknown one)."""
    for g in (X, Y, D):
        if isinstance(table, Catalog):
            table.g(g)
        elif g not in table.labels:
            raise AnalysisError(f"{g} is not a genome of the screen ({', '.join(table.labels)})")
    if len({X, Y, D}) != 3:
        raise ValueError(f"painting needs three distinct genomes, not {X}, {Y}, {D}")


def emissions(cls: np.ndarray, alphabet: str) -> np.ndarray:
    """int64 emission per class code under `alphabet` (Alphabet.codes; -1 missing data). The assembly route's: RD 0,
    RS 1, the rest -1."""
    return ALPHABETS[alphabet].codes[np.asarray(cls, dtype=np.int64)]


def initial_emissions(emit: np.ndarray, alphabet: str) -> tuple[np.ndarray, np.ndarray]:
    """EM's initial (e0, e1) from the emissions `emit` of the decoy (of the target without a decoy): e0 = each
    emitting class's count + 1, e1 = e0 x HMM_INIT (its first factor for a donor-like class, its second for the
    others), each normalized."""
    a = ALPHABETS[alphabet]
    emit = np.asarray(emit, dtype=np.int64)
    counts = np.bincount(emit[emit >= 0], minlength=len(a.emits)).astype(np.float64)
    e0 = counts + 1.0                         # initial e0: the class frequencies + 1
    e0 /= e0.sum()
    e1 = e0 * np.asarray([HMM_INIT[0] if c in a.donor else HMM_INIT[1] for c in a.emits])
    e1 /= e1.sum()
    return e0, e1


def _sites(cat: Catalog, X: str, idx: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(fai index, start, end) of X's site at each locus (1-based)."""
    g = cat.g(X)
    return cat.seq[idx, g].astype(np.int64), cat.start[idx, g], cat.end[idx, g]


def sites(table: Catalog | LibraryRun, X: str, idx: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(sequence code, start, end) of X's painted units `idx`, 1-based: a catalog's sites of X (fai index); a library
    screen's records of X (the rank of the sequence name in text order; a record whose name holds no coordinates is
    its own sequence, from 0 to its length)."""
    if isinstance(table, Catalog):
        return _sites(table, X, idx)
    rec = table.rec(X)
    return table.sequences(X)[1][idx], rec.start[idx], rec.end[idx]


def _position(props: RouteProps, start: np.ndarray, end: np.ndarray) -> np.ndarray:
    """Each unit's position on its sequence under the route's tract position: its start, or its midpoint."""
    return start if props.tract_position == "site_start" else (np.asarray(start) + np.asarray(end)) // 2


def _starts(seq: np.ndarray) -> np.ndarray:
    """Each chain's (sequence's) first index, then a final len(seq) sentinel; seq sorted."""
    return np.r_[0, np.flatnonzero(np.diff(seq)) + 1, len(seq)].astype(np.int64)


# ---------------------------------------------------------------- the HMM


@njit(cache=True)
def forward_backward(pos, emit, starts, pi, lam, e0, e1):
    """Scaled forward-backward over chains [starts[c], starts[c + 1]): categorical
    emissions (emit < 0: missing, emission 1), the switch probability 1 - exp(-lam * d) over d bp, the chain
    starting at (1 - pi, pi). Returns the posterior of state 1, the log-likelihood and the expected emission
    counts under state 1 and under state 0."""
    n = len(pos)
    K = len(e0)
    post = np.zeros(n)
    c1 = np.zeros(K)
    c0 = np.zeros(K)
    ll = 0.0
    for c in range(len(starts) - 1):
        s0 = starts[c]
        m = starts[c + 1] - s0
        if m <= 0:
            continue
        fa = np.zeros((m, 2))
        for t in range(m):
            k = emit[s0 + t]
            if k >= K:
                raise ValueError("forward_backward: an emission code is out of range")
            em0 = e0[k] if k >= 0 else 1.0
            em1 = e1[k] if k >= 0 else 1.0
            if t == 0:
                a0 = 1.0 - pi
                a1 = pi
            else:
                q = math.exp(-lam * max(pos[s0 + t] - pos[s0 + t - 1], 0.0))
                p01 = pi * (1.0 - q)
                p10 = (1.0 - pi) * (1.0 - q)
                a0 = fa[t - 1, 0] * (1.0 - p01) + fa[t - 1, 1] * p10
                a1 = fa[t - 1, 0] * p01 + fa[t - 1, 1] * (1.0 - p10)
            f0 = a0 * em0
            f1 = a1 * em1
            sc = f0 + f1
            fa[t, 0] = f0 / sc
            fa[t, 1] = f1 / sc
            ll += math.log(sc)
        b0 = 1.0
        b1 = 1.0
        for t in range(m - 1, -1, -1):
            if t < m - 1:
                q = math.exp(-lam * max(pos[s0 + t + 1] - pos[s0 + t], 0.0))
                p01 = pi * (1.0 - q)
                p10 = (1.0 - pi) * (1.0 - q)
                k = emit[s0 + t + 1]
                em0 = e0[k] if k >= 0 else 1.0
                em1 = e1[k] if k >= 0 else 1.0
                nb0 = (1.0 - p01) * em0 * b0 + p01 * em1 * b1
                nb1 = p10 * em0 * b0 + (1.0 - p10) * em1 * b1
                ss = nb0 + nb1
                b0 = nb0 / ss
                b1 = nb1 / ss
            g0 = fa[t, 0] * b0
            g1 = fa[t, 1] * b1
            g = g1 / (g0 + g1)
            post[s0 + t] = g
            k = emit[s0 + t]
            if k >= 0:
                c1[k] += g
                c0[k] += 1.0 - g
    return post, ll, c1, c0


def fit_hmm(pos: np.ndarray, emit: np.ndarray, starts: np.ndarray, e0: np.ndarray, e1: np.ndarray) -> HmmFit:
    """EM for the painting HMM (module docstring), from the initial emissions e0, e1.
    A ValueError when no locus emits."""
    pos = np.ascontiguousarray(pos, dtype=np.float64)
    emit = np.ascontiguousarray(emit, dtype=np.int64)
    starts = np.ascontiguousarray(starts, dtype=np.int64)
    e0, e1 = np.asarray(e0, dtype=np.float64), np.asarray(e1, dtype=np.float64)
    informative = emit >= 0
    if not informative.any():
        raise ValueError("fit_hmm: no emitting locus (every emission is missing data)")
    grid = np.geomspace(*HMM_L_GRID)
    pi, L, ll_old = HMM_PI0, HMM_L0, -math.inf
    for _ in range(HMM_ITERS):
        best = None
        for lc in grid:
            post, ll, c1, c0 = forward_backward(pos, emit, starts, pi, 1.0 / lc, e0, e1)
            if best is None or ll > best[1]:
                best = (float(lc), ll, post, c1, c0)
        L, ll, post, c1, c0 = best
        e1 = (c1 + EMIT_PSEUDO) / (c1 + EMIT_PSEUDO).sum()
        e0 = (c0 + EMIT_PSEUDO) / (c0 + EMIT_PSEUDO).sum()
        pi = float(np.clip(post[informative].mean(), PI_MIN, PI_MAX))
        if abs(ll - ll_old) < HMM_TOL * abs(ll):
            break
        ll_old = ll
    post, ll, _, _ = forward_backward(pos, emit, starts, pi, 1.0 / L, e0, e1)
    return HmmFit(pi=pi, L=L, e0=e0, e1=e1, post=post, loglik=float(ll))


def _called(post: np.ndarray, pi: float) -> np.ndarray:
    """The tract rule per locus: posterior > TRACT_POST and posterior odds >= BF_TRACT x prior odds."""
    post = np.asarray(post, dtype=np.float64)
    return (post > TRACT_POST) & (post / np.maximum(1.0 - post, 1e-12) >= BF_TRACT * pi / (1.0 - pi))


def call_tracts(seq: np.ndarray, start: np.ndarray, end: np.ndarray, post: np.ndarray, pi: float, *,
                tract_end: str = "largest_end") -> list[Tract]:
    """Maximal runs, on one sequence, of loci passing the tract rule: from the
    first locus's start to the largest end of the run's loci (tract_end "largest_end"), or to its last locus's end
    ("last_end"); sum_post exactly rounded (math.fsum). Arrays in position order."""
    require("tract end", tract_end, ("largest_end", "last_end"))
    seq, start, end = np.asarray(seq), np.asarray(start), np.asarray(end)
    post = np.asarray(post, dtype=np.float64)
    on = _called(post, pi)
    if not on.any():
        return []
    same = seq[1:] == seq[:-1]
    head = on.copy()
    head[1:] &= ~(on[:-1] & same)
    tail = on.copy()
    tail[:-1] &= ~(on[1:] & same)
    out = []
    for a, b in zip(np.flatnonzero(head).tolist(), np.flatnonzero(tail).tolist(), strict=True):
        s = math.fsum(post[a:b + 1].tolist())  # exactly rounded: ties stay ties (module docstring)
        stop = int(end[a:b + 1].max()) if tract_end == "largest_end" else int(end[b])
        out.append(Tract(seq=int(seq[a]), start=int(start[a]), end=stop, first=a, last=b,
                         mean_post=s / (b - a + 1), sum_post=s))
    return out


def target_decoy(target: np.ndarray, decoy: np.ndarray, scale: float) -> tuple[np.ndarray, np.ndarray]:
    """(q, p) per target score (module docstring): q from the decoy counts times `scale`
    with the +1 correction, so a lone target above every decoy is not significant; p the empirical p-value
    against the decoy scores."""
    target = np.asarray(target, dtype=np.float64)
    null = np.sort(np.asarray(decoy, dtype=np.float64))
    obs = np.sort(target)
    n_null = len(null) - np.searchsorted(null, target, side="left")   # decoy scores >= x
    n_obs = len(obs) - np.searchsorted(obs, target, side="left")      # target scores >= x
    raw = (n_null * scale + 1.0) / np.maximum(n_obs, 1)
    order = np.argsort(-target, kind="stable")
    q = np.empty(len(target))
    q[order] = np.minimum(1.0, np.minimum.accumulate(raw[order][::-1])[::-1])
    return q, (1.0 + n_null) / (1.0 + len(null))


def shortest_callable(fit: HmmFit, spacing: float
                      ) -> tuple[tuple[int | None, float | None], tuple[int | None, float | None]]:
    """((loci, bp) without flanks, (loci, bp) between CALLABLE_FLANK RS loci on each side): the fewest
    consecutive RD loci, `spacing` bp apart, whose middle locus (index k // 2 of the run) passes the tract
    rule under `fit`; bp = loci x spacing. (None, None) when no k up to CALLABLE_MAX does, or without a
    finite spacing."""
    if not (math.isfinite(spacing) and spacing >= 0):
        return (None, None), (None, None)
    e0, e1 = np.asarray(fit.e0, dtype=np.float64), np.asarray(fit.e1, dtype=np.float64)
    out = []
    for flank in (0, CALLABLE_FLANK):
        found: tuple[int | None, float | None] = (None, None)
        for k in range(1, CALLABLE_MAX + 1):
            emit = np.r_[np.full(flank, RS), np.full(k, RD), np.full(flank, RS)].astype(np.int64)
            pos = np.arange(len(emit), dtype=np.float64) * spacing
            post = forward_backward(pos, emit, np.array([0, len(emit)], dtype=np.int64), fit.pi, 1.0 / fit.L,
                                    e0, e1)[0]
            if _called(post[flank + k // 2:flank + k // 2 + 1], fit.pi)[0]:
                found = (k, k * float(spacing))
                break
        out.append(found)
    return out[0], out[1]


# ---------------------------------------------------------------- painting


def paint(table: Catalog | LibraryRun, trio: Trio, props: RouteProps, *, alphabet: str | None = None,
          blocks: list[InvBlock] = (), self_counts: dict[tuple[str, str], tuple[int, int]] | None = None
          ) -> Painting:
    """Paint trio.R (module docstring) in `alphabet` (one of the route's; None: the route's only one) under the
    route properties `props`: classes, the HMM fit, tracts, the decoy painting of trio.S with the fitted
    parameters, q and p, and the inverted-block flags (blocks: `inverted_blocks` or `chain_inversions`).
    self_counts: (lid, genome) -> (sites, diffs) from the identity layer's self rows, when it ran. table: a run's
    catalog or a library screen (module docstring). Skipped (Painting.skipped) when R has no annotation, or holds
    no unit of a donor-like class or none of a sister-like one (on the assembly route: no RD or no RS locus). A
    ValueError, before anything is painted, for a route property this engine does not implement."""
    if alphabet is None:
        if len(props.alphabet) != 1:
            raise ValueError(f"the route paints the alphabets {', '.join(props.alphabet)}: name one")
        alphabet = props.alphabet[0]
    elif alphabet not in props.alphabet:
        raise ValueError(f"painting alphabet {alphabet!r} is not one of the route's ({', '.join(props.alphabet)})")
    _input_alphabet(table, alphabet)
    _require(props, alphabet)
    X, Y, D = trio.R, trio.S, trio.D
    _check(table, X, Y, D)
    catalog = isinstance(table, Catalog)
    near = None if catalog else near_threshold(table, Y, D)
    none = np.zeros(0, dtype=np.int64)
    blank = {"decoy_idx": none, "decoy_cls": none.astype(np.int8), "decoy_post": np.zeros(0),
             "decoy_tracts": [], "fit": None, "donor_state": False, "tracts": [], "q": [], "p": [],
             "self_counts": self_counts, "props": props, "alphabet": alphabet, "near": near}
    if catalog and not table.annotated[X]:
        why = (f"{X} has no annotation, so it discovers no loci to paint (painting reads the "
               "loci the recipient discovers)")
        return Painting(trio=trio, idx=none, cls=none.astype(np.int8), post=np.zeros(0),
                        inv=np.zeros(0, dtype=bool), skipped=why, **blank)
    idx, cls = classes(table, X, Y, D, alphabet=alphabet, near=near)
    seq, start, end = sites(table, X, idx)
    pos = _position(props, start, end)
    inv = inside(blocks, X, {Y, D}, seq, pos)
    emit = emissions(cls, alphabet)
    why = _lacking(table, alphabet, X, Y, D, cls, emit)
    if why is not None:
        return Painting(trio=trio, idx=idx, cls=cls, post=np.full(len(idx), np.nan), inv=inv, skipped=why, **blank)
    decoy = None
    if not catalog or table.annotated[Y]:
        didx, dcls = classes(table, Y, X, D, alphabet=alphabet, near=near)
        demit = emissions(dcls, alphabet)
        if props.decoy_scale == "all" or np.count_nonzero(demit >= 0):
            decoy = (didx, dcls, demit)
    e0, e1 = initial_emissions(decoy[2] if decoy is not None else emit, alphabet)
    fit = fit_hmm(pos.astype(np.float64), emit, _starts(seq), e0, e1)
    donor_state = _donor_state(props.donor_state, alphabet, fit)
    tracts = call_tracts(seq, start, end, fit.post, fit.pi, tract_end=props.tract_end) if donor_state else []
    if decoy is None:
        return Painting(trio=trio, idx=idx, cls=cls, post=fit.post, inv=inv, fit=fit, donor_state=donor_state,
                        tracts=tracts, q=[None] * len(tracts), p=[None] * len(tracts), decoy_idx=none,
                        decoy_cls=none.astype(np.int8), decoy_post=np.zeros(0), decoy_tracts=[], skipped=None,
                        self_counts=self_counts, props=props, alphabet=alphabet, near=near)
    didx, dcls, demit = decoy
    dseq, dstart, dend = sites(table, Y, didx)
    dpos = _position(props, dstart, dend)
    dpost = forward_backward(dpos.astype(np.float64), demit, _starts(dseq), fit.pi, 1.0 / fit.L, fit.e0,
                             fit.e1)[0]
    dtracts = (call_tracts(dseq, dstart, dend, dpost, fit.pi, tract_end=props.tract_end) if donor_state
               else [])
    if props.decoy_scale == "informative":
        scale = np.count_nonzero(emit >= 0) / np.count_nonzero(demit >= 0)
    else:
        scale = len(emit) / max(len(demit), 1)
    q, p = target_decoy(np.array([t.sum_post for t in tracts]), np.array([t.sum_post for t in dtracts]),
                        scale)
    return Painting(trio=trio, idx=idx, cls=cls, post=fit.post, inv=inv, fit=fit, donor_state=donor_state,
                    tracts=tracts, q=q.tolist(), p=p.tolist(), decoy_idx=didx, decoy_cls=dcls,
                    decoy_post=dpost, decoy_tracts=dtracts, skipped=None,
                    decoy_inv=inside(blocks, Y, {X, D}, dseq, dpos), scale=float(scale), self_counts=self_counts,
                    props=props, alphabet=alphabet, near=near)


def _lacking(table: Catalog | LibraryRun, alphabet: str, X: str, Y: str, D: str, cls: np.ndarray,
             emit: np.ndarray) -> str | None:
    """Why X cannot be painted: no unit of a donor-like class or none of a sister-like one (the HMM needs both);
    None when it can."""
    a = ALPHABETS[alphabet]
    if isinstance(table, Catalog):
        meaning = {RD: f"{D} present, {Y} E", RS: f"{Y} present, {D} E"}
        lacking = [c for c in (RD, RS) if not np.count_nonzero(cls == c)]
        if not lacking:
            return None
        why = " and no ".join(f"{PAINT_CLASSES[c]} locus ({meaning[c]})" for c in lacking)
        return f"no {why} in {X}'s own calls at the {len(cls)} loci of {X}: the HMM needs both classes"
    kinds = (("donor-like", [c for c in a.emits if c in a.donor]), ("sister-like", [c for c in a.emits
                                                                                 if c not in a.donor]))
    lacking = [f"{kind} element ({' or '.join(names)})" for kind, names in kinds
               if not np.isin(emit, [a.emits.index(c) for c in names]).any()]
    if not lacking:
        return None
    return (f"no {' and no '.join(lacking)} among the {len(cls)} records of {X} (sister {Y}, donor {D}): the HMM "
            "needs both kinds")


def _donor_state(rule: str, alphabet: str, fit: HmmFit) -> bool:
    """The route's donor-state check (DONOR_STATE): state 1's emission above state 0's for each `up` class and below
    it for each `down` class."""
    up, down = DONOR_STATE[rule]
    k = ALPHABETS[alphabet].emits.index
    return bool(all(fit.e1[k(c)] > fit.e0[k(c)] for c in up) and all(fit.e1[k(c)] < fit.e0[k(c)] for c in down))


def _require(props: RouteProps, alphabet: str) -> None:
    """A ValueError for a painting property of the route this engine does not implement (module docstring), or a
    donor-state check that reads a class the alphabet does not emit."""
    require("painting alphabet", alphabet, tuple(ALPHABETS))
    require("decoy scale", props.decoy_scale, ("informative", "all"))
    require("donor state", props.donor_state, tuple(DONOR_STATE))
    require("tract position", props.tract_position, ("site_start", "midpoint"))
    require("tract end", props.tract_end, ("largest_end", "last_end"))
    require("posterior mean", props.posterior_mean, ("unweighted", "length_weighted"))
    missing = [c for c in sum(DONOR_STATE[props.donor_state], ()) if c not in ALPHABETS[alphabet].emits]
    if missing:
        raise ValueError(f"donor state {props.donor_state!r} reads the class(es) {', '.join(missing)}, which the "
                         f"painting alphabet {alphabet!r} does not emit")


def posterior_mean(p: Painting, table: Catalog | LibraryRun, *, decoy: bool = False) -> float:
    """The mean posterior over the painted units of R (decoy: of S) under the route's posterior mean: each unit
    once ("unweighted"), or weighted by its length end - start + 1 ("length_weighted"); nan without a unit."""
    idx, post, X = (p.decoy_idx, p.decoy_post, p.trio.S) if decoy else (p.idx, p.post, p.trio.R)
    if not len(post):
        return math.nan
    if p.props is None or p.props.posterior_mean == "unweighted":
        return float(post.mean())
    _, start, end = sites(table, X, idx)
    return float(np.average(post, weights=end - start + 1))


# ---------------------------------------------------------------- per tract


def _tract_table(p: Painting, cat: Catalog, decoy: bool = False) -> list[dict]:
    """Per tract of the painting (or of its decoy): n_loci, n_RD, n_RS, youngest_k, youngest_u95, n_identical,
    inv_frac (module docstring; None where a value does not exist)."""
    X, D = (p.trio.S if decoy else p.trio.R), p.trio.D
    idx, cls, inv, tracts = ((p.decoy_idx, p.decoy_cls, p.decoy_inv, p.decoy_tracts) if decoy
                             else (p.idx, p.cls, p.inv, p.tracts))
    gx = cat.g(X)
    out = []
    for t in tracts:
        sl = slice(t.first, t.last + 1)
        c = cls[sl]
        rd = idx[sl][c == RD]
        k = cat.k[rd, gx]
        known = k[np.isfinite(k)]
        out.append({"n_loci": t.last - t.first + 1, "n_RD": len(rd), "n_RS": int(np.count_nonzero(c == RS)),
                    "youngest_k": float(known.min()) if known.size else None,
                    "youngest_u95": _youngest_u95(cat, p.self_counts, rd, X, D),
                    "n_identical": int(np.count_nonzero(age_bin(k) == 0)),
                    "inv_frac": float(np.count_nonzero(inv[sl])) / (t.last - t.first + 1)})
    return out


def _youngest_u95(cat: Catalog, counts: dict[tuple[str, str], tuple[int, int]] | None, rd: np.ndarray, X: str,
                  D: str) -> float | None:
    """The smallest u95(diffs, sites) over the RD loci `rd`, from X's self row, else D's; loci without a
    self row (with sites) are left out. None without self_counts or without any such row."""
    if counts is None:
        return None
    sites, diffs = [], []
    for lid in cat.lid[rd].tolist():
        for g in (X, D):
            c = counts.get((lid, g))
            if c is not None and c[0] > 0:
                sites.append(c[0])
                diffs.append(c[1])
                break
    return float(np.min(garwood_u95(np.array(diffs), np.array(sites)))) if sites else None


def _record_tract_table(p: Painting, run: LibraryRun, decoy: bool = False) -> list[dict]:
    """`_tract_table` on the library route: n_RD and n_RS the tract's elements of a donor-like and of a sister-like
    class, youngest_k and n_identical from the painted genome's annotated K2P over its donor-like elements,
    youngest_u95 None (no identity layer)."""
    X = p.trio.S if decoy else p.trio.R
    idx, cls, inv, tracts = ((p.decoy_idx, p.decoy_cls, p.decoy_inv, p.decoy_tracts) if decoy
                             else (p.idx, p.cls, p.inv, p.tracts))
    a = ALPHABETS[p.alphabet]
    emit = emissions(cls, p.alphabet)
    donor = np.isin(emit, [a.emits.index(c) for c in a.emits if c in a.donor])
    k2p = run.rec(X).k2p
    out = []
    for t in tracts:
        sl = slice(t.first, t.last + 1)
        k = k2p[idx[sl][donor[sl]]]
        known = k[np.isfinite(k)]
        out.append({"n_loci": t.last - t.first + 1, "n_RD": int(np.count_nonzero(donor[sl])),
                    "n_RS": int(np.count_nonzero((emit[sl] >= 0) & ~donor[sl])),
                    "youngest_k": float(known.min()) if known.size else None, "youngest_u95": None,
                    "n_identical": int(np.count_nonzero(age_bin(k) == 0)),
                    "inv_frac": float(np.count_nonzero(inv[sl])) / (t.last - t.first + 1)})
    return out


# ---------------------------------------------------------------- files


def write_paintings(paintings: list[Painting], cat: Catalog | LibraryRun, outdir: Path, previous=()) -> list[Path]:
    """Remove the old painting files, then write tracts.<R>.bed and posterior.<R>.tsv per recipient and
    decoy_tracts.<S>.bed and decoy_posterior.<S>.tsv per sister with a decoy (module docstring); the paths
    written. Skipped paintings write nothing. Removed: the painting files named for the run's labels, and those
    in `outdir` that `previous` lists (the painting files the run before recorded writing, run.PAINTINGS:
    work/geneflow/paintings.json; a painting of a genome no longer in the sheet), never another file, so a genome
    labeled "tracts" keeps the engine's tracts.pav.bed and a file of the user's own stays. On a library screen
    (`cat` a LibraryRun) the files hold every alphabet (module docstring, the library route)."""
    outdir = Path(outdir)
    old = {outdir / f.format(x) for x in cat.labels for f in PAINTING_FILES}
    old |= {Path(q) for q in previous if Path(q).resolve().parent == outdir.resolve()
            and any(fnmatch.fnmatchcase(Path(q).name, f.format("*")) for f in PAINTING_FILES)}
    for q in old:
        q.unlink(missing_ok=True)
    if not isinstance(cat, Catalog):
        return _write_library(paintings, cat, outdir)
    target: dict[str, list[Painting]] = {}
    mirror: dict[str, list[Painting]] = {}
    for p in sorted((p for p in paintings if p.skipped is None), key=lambda p: p.trio.scope):
        target.setdefault(p.trio.R, []).append(p)
        if p.scale is not None:
            mirror.setdefault(p.trio.S, []).append(p)
    paths = []
    for files, decoy in ((target, False), (mirror, True)):
        for x, ps in files.items():
            stem = "decoy_" if decoy else ""
            paths.append(_write_tracts(outdir / f"{stem}tracts.{x}.bed", ps, cat, decoy))
            paths.append(_write_posterior(outdir / f"{stem}posterior.{x}.tsv", ps, cat, decoy))
    return paths


def _sides(p: Painting, decoy: bool) -> tuple[str, str]:
    """(the painted genome, its sister): (R, S), or (S, R) for the decoy."""
    return (p.trio.S, p.trio.R) if decoy else (p.trio.R, p.trio.S)


def _write_tracts(path: Path, ps: list[Painting], cat: Catalog, decoy: bool) -> Path:
    with open(path, "w") as fh:
        fh.write("\t".join(TRACT_COLUMNS) + "\n")
        for p in ps:
            x, y = _sides(p, decoy)
            tracts = p.decoy_tracts if decoy else p.tracts
            for j, (t, v) in enumerate(zip(tracts, _tract_table(p, cat, decoy), strict=True)):
                if decoy:
                    q = pv = None                   # a decoy tract is the null, not a test
                elif p.scale is None:
                    q = pv = "none"                 # no decoy (WARNING in findings)
                else:
                    q, pv = p.q[j], p.p[j]
                row = (cat.seqids[x][t.seq], t.start - 1, t.end, f"{p.trio.D}-like",
                       min(1000, int(1000 * t.mean_post)), ".", y, p.trio.D, v["n_loci"], v["n_RD"],
                       v["n_RS"], t.mean_post, t.sum_post, q, pv, v["youngest_k"], v["youngest_u95"],
                       v["n_identical"], v["inv_frac"])
                fh.write("\t".join(tsv_field(c) for c in row) + "\n")
    return path


def _write_posterior(path: Path, ps: list[Painting], cat: Catalog, decoy: bool) -> Path:
    with open(path, "w") as fh:
        fh.write("\t".join(POSTERIOR_COLUMNS) + "\n")
        for p in ps:
            x, y = _sides(p, decoy)
            idx, cls, post, inv, tracts = (
                (p.decoy_idx, p.decoy_cls, p.decoy_post, p.decoy_inv, p.decoy_tracts) if decoy
                else (p.idx, p.cls, p.post, p.inv, p.tracts))
            names = cat.seqids[x]
            seq, start, end = _sites(cat, x, idx)
            label = np.full(len(idx), ".", dtype=object)
            for t in tracts:
                label[t.first:t.last + 1] = f"{names[t.seq]}:{t.start}-{t.end}"
            mid = f"\t{y}\t{p.trio.D}\t"
            cols = zip(cat.lid[idx].tolist(), seq.tolist(), start.tolist(), end.tolist(), cls.tolist(),
                       post.tolist(), label.tolist(), inv.tolist(), strict=True)
            fh.writelines(f"{lid}\t{names[s]}\t{a}\t{b}{mid}{PAINT_CLASSES[c]}\t{tsv_field(v)}\t{t}\t{int(i)}\n"
                          for lid, s, a, b, c, v, t, i in cols)
    return path


def exact_field(v) -> str:
    """A TSV value written exactly: None and nan -> "."; bool -> "1"/"0"; int as is; float as the shortest text that
    reads back as the same number (repr); str as is."""
    if isinstance(v, (float, np.floating)) and not isinstance(v, bool):
        return "." if math.isnan(v) else repr(float(v))
    return tsv_field(v)


def _write_library(paintings: list[Painting], run: LibraryRun, outdir: Path) -> list[Path]:
    """The library route's painting files (module docstring): per genome one file of each kind, holding every trio
    (scope order) and every alphabet; a trio whose paintings are all skipped writes no row."""
    trios: dict[str, list[Painting]] = {}
    for p in sorted(paintings, key=lambda p: p.trio.scope):
        trios.setdefault(p.trio.scope, []).append(p)
    alphabets = tuple(dict.fromkeys(p.alphabet for p in paintings))
    target: dict[str, list[list[Painting]]] = {}
    mirror: dict[str, list[list[Painting]]] = {}
    for ps in trios.values():
        if any(p.skipped is None for p in ps):
            target.setdefault(ps[0].trio.R, []).append(ps)
        if any(p.scale is not None for p in ps):
            mirror.setdefault(ps[0].trio.S, []).append(ps)
    paths = []
    for files, decoy in ((target, False), (mirror, True)):
        for x, groups in files.items():
            stem = "decoy_" if decoy else ""
            paths.append(_write_record_tracts(outdir / f"{stem}tracts.{x}.bed", groups, run, decoy))
            paths.append(_write_record_posterior(outdir / f"{stem}posterior.{x}.tsv", groups, run, decoy, alphabets))
    return paths


def _write_record_tracts(path: Path, groups: list[list[Painting]], run: LibraryRun, decoy: bool) -> Path:
    with open(path, "w") as fh:
        fh.write("\t".join(LIBRARY_TRACT_COLUMNS) + "\n")
        for ps in groups:
            for p in ps:
                x, y = _sides(p, decoy)
                tracts = p.decoy_tracts if decoy else p.tracts
                names = run.seqids[x]
                for j, (t, v) in enumerate(zip(tracts, _record_tract_table(p, run, decoy), strict=True)):
                    q, pv = (None, None) if decoy else (p.q[j], p.p[j])
                    row = (names[t.seq], t.start - 1, t.end, f"{p.trio.D}-{p.alphabet}",
                           min(1000, int(1000 * t.mean_post)), ".", y, p.trio.D, v["n_loci"], v["n_RD"], v["n_RS"],
                           t.mean_post, t.sum_post, q, pv, v["youngest_k"], v["youngest_u95"], v["n_identical"],
                           v["inv_frac"], p.alphabet)
                    fh.write("\t".join(exact_field(c) for c in row) + "\n")
    return path


def _write_record_posterior(path: Path, groups: list[list[Painting]], run: LibraryRun, decoy: bool,
                            alphabets: tuple[str, ...]) -> Path:
    with open(path, "w") as fh:
        fh.write("\t".join(library_posterior_columns(alphabets)) + "\n")
        for ps in groups:
            painted = [p for p in ps if (p.scale is not None if decoy else p.skipped is None)]
            ref = painted[0]
            x, y = _sides(ref, decoy)
            idx, cls, inv = (ref.decoy_idx, ref.decoy_cls, ref.decoy_inv) if decoy else (ref.idx, ref.cls, ref.inv)
            post = {p.alphabet: p.decoy_post if decoy else p.post for p in painted}
            for p in painted:
                if not np.array_equal(p.decoy_idx if decoy else p.idx, idx):
                    raise ValueError(f"{p.trio.scope}: the alphabets painted different records")
            cols = [post.get(a, np.full(len(idx), np.nan)).tolist() for a in alphabets]
            rec = run.rec(x)
            eid = rec.lib.elements.column("eid")
            mid = f"\t{y}\t{ref.trio.D}\t"
            for j, i in enumerate(idx.tolist()):
                if eid[i] is None:            # a record whose name holds no coordinates: no element, no row
                    continue
                values = "\t".join(exact_field(c[j]) for c in cols)
                fh.write(f"{eid[i]}\t{rec.chrom[i]}\t{rec.start[i]}\t{rec.end[i]}{mid}{LIBRARY_CLASSES[cls[j]]}\t"
                         f"{values}\t{int(inv[j])}\n")
    return path


# ---------------------------------------------------------------- rows and findings


def _loci(p: Painting, x: str, rest: str) -> str:
    """A row's loci declaration (relations.py's grammar) for loci painted on x."""
    return f"disc={x};present={x};{rest}"


def _spacing(cat: Catalog, p: Painting) -> np.ndarray:
    """bp between consecutive informative (RD/RS) loci on one sequence of R (site starts)."""
    keep = p.cls <= RS
    seq, start, _ = _sites(cat, p.trio.R, p.idx[keep])
    gap = np.diff(start)
    return gap[seq[1:] == seq[:-1]].astype(np.float64)


def paint_rows(p: Painting, cat: Catalog) -> tuple[list[StatRow], Findings]:
    """Rows (layer pav; test paint, paint_resolution and paint_seqid; scope the trio) and findings (sections
    Painting and Dating) of one painting. A skipped painting gives its class counts and one WARNING."""
    f = Findings()
    X, Y, D, scope = p.trio.R, p.trio.S, p.trio.D, p.trio.scope
    inf_loci = _loci(p, X, f"res={Y},{D};outgroup=not_present")
    all_loci = _loci(p, X, "outgroup=none")
    n_all = len(p.idx)
    counts = {c: int(np.count_nonzero(p.cls == c)) for c in range(len(PAINT_CLASSES))}
    n_rd, n_rs = counts[RD], counts[RS]
    n_inf = n_rd + n_rs
    rows: list[StatRow] = []
    if p.skipped is not None:
        f.warn(PAINTING, f"Painting of {X} (sister {Y}, donor {D}) skipped: {p.skipped}")
        if n_all:
            rows += [StatRow("pav", "paint", scope, "all", s, v, n=n_all, loci=inf_loci, note=note)
                     for s, v, note in (("n_informative", n_inf, "RD and RS loci"),
                                        ("n_RD", n_rd, f"{D} present, {Y} E"),
                                        ("n_RS", n_rs, f"{Y} present, {D} E"))]
            rows += [StatRow("pav", "paint", scope, "all", s, v, n=n, loci=loci, note=note)
                     for s, v, n, loci, note in _emission_stats(p, True, n_inf, inf_loci,
                                                                "NA: the painting was skipped")]
            rows += _resolution(cat, p, scope, all_loci, f)
        return rows, f
    fit = p.fit
    decoy = p.scale is not None
    post_mean = posterior_mean(p, cat)
    dec_mean = posterior_mean(p, cat, decoy=True) if decoy else math.nan
    excess = post_mean - dec_mean if decoy else math.nan
    n_q05 = sum(q is not None and q <= Q_TRACT for q in p.q) if decoy else math.nan
    gaps = _spacing(cat, p)
    sp_med, sp_mean = (float(np.median(gaps)), float(gaps.mean())) if gaps.size else (math.nan, math.nan)
    na = p.collapsed  # NA, never a pseudocount artifact (Painting.collapsed)
    (lo_k, lo_bp), (hi_k, hi_bp) = (shortest_callable(fit, sp_med) if p.donor_state and not na
                                    else ((None, None), (None, None)))
    n_inv = int(np.count_nonzero(p.inv))
    run_bp = math.nan if na else fit.L / (1.0 - fit.pi)
    em = (f"e1 (RD, RS) {fit.e1[RD]:.4g}, {fit.e1[RS]:.4g}; e0 {fit.e0[RD]:.4g}, {fit.e0[RS]:.4g}; state 1 "
          f"is {D}-like when e1[RD] > e0[RD]")
    na_note = (f"NA: the fit collapsed (pi at its floor {PI_MIN:g}, no tract), so state 1 holds no data and "
               "this value would come from the pseudocounts")
    q_note = (f"tracts at q <= {Q_TRACT:g} (target-decoy against {Y}'s painting; scale = n informative {X} / n "
              f"informative {Y} = {p.scale:.4g})" if decoy else f"no decoy: {_no_decoy(cat, p)}")
    callable_note = (f"the fewest consecutive RD loci at the median spacing whose middle locus passes the "
                     f"tract rule; up to {CALLABLE_MAX} loci")
    spacing_note = "bp between consecutive informative loci on one sequence"
    stats = [
        ("pi", fit.pi, n_inf, inf_loci,
         f"the fitted {D}-like fraction (EM; clipped to [{PI_MIN:g}, {PI_MAX:g}])"
         + ("; at its floor: the fit collapsed" if na else "")),
        ("mean_tract_bp", math.nan if na else fit.L, n_inf, inf_loci,
         na_note if na else
         ("L, the switch scale: between loci d bp apart the state is redrawn from (1 - pi, pi) with "
          "probability 1 - exp(-d / L)")),
        ("mean_donor_run_bp", run_bp, n_inf, inf_loci,
         na_note if na else f"L / (1 - pi): the mean {D}-like run under the model"),
        ("donor_state", math.nan if na else int(p.donor_state), n_inf, inf_loci,
         f"{na_note}; {em}" if na else em),
        *_emission_stats(p, na, n_inf, inf_loci, na_note),
        ("loglik", fit.loglik, n_inf, inf_loci, "log-likelihood of the fitted HMM"),
        ("n_informative", n_inf, n_all, inf_loci, "RD and RS loci: the emitting loci"),
        ("n_RD", n_rd, n_all, inf_loci, f"{D} present, {Y} E"),
        ("n_RS", n_rs, n_all, inf_loci, f"{Y} present, {D} E"),
        ("posterior_mean", post_mean, n_all, all_loci,
         f"mean posterior of the {D}-like state over {X}'s painted loci, unweighted"),
        ("decoy_posterior_mean", dec_mean, len(p.decoy_idx), _loci(p, Y, "outgroup=none"),
         f"the same over the decoy {Y}, painted with {X}'s fitted parameters" if decoy
         else f"no decoy: {_no_decoy(cat, p)}"),
        ("excess", excess, n_all, all_loci, "posterior_mean - decoy_posterior_mean"),
        ("n_tracts", len(p.tracts), n_inf, inf_loci,
         f"runs with posterior > 0.5 and posterior odds >= {BF_TRACT:g} x prior odds"),
        ("n_decoy_tracts", len(p.decoy_tracts) if decoy else math.nan, len(p.decoy_idx),
         _loci(p, Y, f"res={X},{D};outgroup=not_present"), "the same rule on the decoy"),
        ("n_tracts_q05", n_q05, n_inf, inf_loci, q_note),
        ("spacing_median", sp_med, n_inf, inf_loci, spacing_note),
        ("spacing_mean", sp_mean, n_inf, inf_loci, spacing_note),
        ("shortest_callable_loci_lo", num(lo_k), n_inf, inf_loci,
         na_note if na else f"{callable_note}, without flanks"),
        ("shortest_callable_bp_lo", num(lo_bp), n_inf, inf_loci,
         na_note if na else "shortest_callable_loci_lo x spacing_median"),
        ("shortest_callable_loci_hi", num(hi_k), n_inf, inf_loci,
         na_note if na else
         f"{callable_note}, between {CALLABLE_FLANK} RS loci on each side (a heuristic flank)"),
        ("shortest_callable_bp_hi", num(hi_bp), n_inf, inf_loci,
         na_note if na else "shortest_callable_loci_hi x spacing_median"),
        ("n_inv_loci", n_inv, n_all, all_loci,
         (f"painted loci inside an inverted block of {X} against {Y} or {D} (flags only; the fit ignores "
          "them)")),
    ]
    rows += [StatRow("pav", "paint", scope, "all", s, v, n=n, loci=loci, note=note)
             for s, v, n, loci, note in stats]
    rows += _resolution(cat, p, scope, all_loci, f)
    rows += _seqid_rows(cat, p, scope, inf_loci, all_loci)
    _painting_lines(cat, p, f, post_mean, dec_mean, excess, n_q05, counts, (sp_med, sp_mean),
                    ((lo_k, lo_bp), (hi_k, hi_bp)), n_inv, run_bp)
    _dating_lines(cat, p, f)
    return rows, f


def library_rows(p: Painting, run: LibraryRun) -> tuple[list[StatRow], Findings]:
    """Rows (layer library, test hmm_<alphabet>, scope the trio, stratum the recipient; module docstring, the library
    route) and findings (section Painting) of one library painting. A collapsed fit gives NA for L, the emissions
    and the donor state; a skipped painting NA for every fitted value, and one WARNING."""
    f = Findings()
    X, Y, D, scope, a = p.trio.R, p.trio.S, p.trio.D, p.trio.scope, p.alphabet
    label, names = PAINT_LABELS[a], ALPHABETS[a].emits

    def row(stat: str, value, note: str = ".") -> StatRow:
        return StatRow(LIBRARY_LAYER, f"hmm_{a}", scope, X, stat, value, note=note)
    near = (f"divergence below it is near: the {NEAR_QUANTILE:g} quantile of {Y}'s ortholog divergences to {D} "
            f"(n_diff / n_sites; 0 when none), one threshold for {X} and the decoy {Y}")
    head = [row("sister", Y, f"the decoy, painted with {X}'s fitted model"), row("donor", D)]
    emission = (f"state 0's ({Y}-like) emission probability of", f"state 1's ({D}-like) emission probability of")
    if p.skipped is not None:
        f.warn(PAINTING, f"{label} of {X} (sister {Y}, donor {D}) skipped: {p.skipped}")
        na = f"NA: the painting was skipped ({p.skipped})"
        fitted = ("pi", "mean_tract_bp", "loglik", "n_tracts", "n_decoy_tracts", "n_tracts_q05", "posterior_mean",
                  "decoy_posterior_mean", "excess_over_decoy")
        return [*head, *(row(s, math.nan, na) for s in fitted), row("near_threshold", p.near, near),
                *(row(f"e{s}_{c}", math.nan, na) for s in (1, 0) for c in names),
                row("donor_state", math.nan, na)], f
    fit, na = p.fit, p.collapsed
    post_mean, dec_mean = posterior_mean(p, run), posterior_mean(p, run, decoy=True)
    excess = post_mean - dec_mean
    n_q05 = sum(q <= Q_TRACT for q in p.q)
    na_note = (f"NA: the fit collapsed (pi at its floor {PI_MIN:g}, no tract), so state 1 holds no data and this "
               "value would come from the pseudocounts")
    up, down = DONOR_STATE[p.props.donor_state]
    state = (f"state 1 is {D}-like when its {' and '.join(up)} emission exceeds state 0's"
             + (f" and its {' and '.join(down)} emission falls below it" if down else ""))
    rows = [*head,
            row("pi", fit.pi, f"the fitted {D}-like fraction (EM; clipped to [{PI_MIN:g}, {PI_MAX:g}])"
                + ("; at its floor: the fit collapsed" if na else "")),
            row("mean_tract_bp", math.nan if na else fit.L,
                na_note if na else ("L, the switch scale: between elements d bp apart the state is redrawn from "
                                    "(1 - pi, pi) with probability 1 - exp(-d / L)")),
            row("loglik", fit.loglik, "log-likelihood of the fitted HMM"),
            row("n_tracts", len(p.tracts), f"runs with posterior > 0.5 and posterior odds >= {BF_TRACT:g} x prior "
                                           "odds"),
            row("n_decoy_tracts", len(p.decoy_tracts), "the same rule on the decoy"),
            row("n_tracts_q05", n_q05, f"tracts at q <= {Q_TRACT:g} (target-decoy against {Y}'s painting; scale = "
                                       f"records of {X} / records of {Y} = {p.scale:.4g})"),
            row("posterior_mean", post_mean, f"mean posterior of the {D}-like state over {X}'s records, weighted by "
                                             "length" if p.props.posterior_mean == "length_weighted" else
                f"mean posterior of the {D}-like state over {X}'s records"),
            row("decoy_posterior_mean", dec_mean, f"the same over the decoy {Y}, painted with {X}'s fitted model"),
            row("excess_over_decoy", excess, "posterior_mean - decoy_posterior_mean"),
            row("near_threshold", p.near, near),
            *(row(f"e{s}_{c}", math.nan if na else float((fit.e1 if s == 1 else fit.e0)[k]),
                  na_note if na else f"{emission[s]} {c}, fitted by EM")
              for s in (1, 0) for k, c in enumerate(names)),
            row("donor_state", math.nan if na else int(p.donor_state), f"{na_note}; {state}" if na else state)]
    if not p.donor_state and not na:
        f.warn(PAINTING, f"{label} of {X}: no {D}-like state (the fitted second state is not enriched for "
                         f"near-identical {D} orthologs): no tracts")
    f.add(PAINTING, f"{label} of {X} (sister {Y}, donor {D}): {D}-like posterior mean {post_mean:.3f} vs "
                    f"{dec_mean:.3f} in the decoy {Y} (excess {excess:+.3f}); "
                    + (f"pi at its floor {PI_MIN:g}" if na else
                       f"switch scale L (mean_tract_bp) {_kb(fit.L)}, pi {fit.pi:.3f}")
                    + f"; {len(p.tracts)} tract(s) ({len(p.decoy_tracts)} decoy tracts), {n_q05} at q <= "
                      f"{Q_TRACT:g}")
    if na:
        f.add(PAINTING, f"{label} of {X}: the fit collapsed (pi at its floor {PI_MIN:g}, no tract): no {D}-like "
                        "ancestry is detectable in this trio. L (mean_tract_bp), the emissions and donor_state are "
                        "NA: state 1 holds no data, so they would come from the pseudocounts")
    if fit.pi >= PI_MAX:
        f.add(PAINTING, f"{label} of {X}: the fitted {D}-like fraction pi reached its cap PI_MAX = {PI_MAX:g} (a "
                        f"recipient more {D}-like than {Y}-like would be {D}'s sister: check the roles)")
    return rows, f


def _no_decoy(cat: Catalog, p: Painting) -> str:
    """Why a painting has no decoy."""
    Y = p.trio.S
    if not cat.annotated[Y]:
        return f"{Y} has no annotation, so it discovers no loci to paint"
    return f"{Y}'s loci hold no informative (RD or RS) locus"


def _emission_stats(p: Painting, na: bool, n_inf: int, inf_loci: str, na_note: str) -> list[tuple]:
    """The emission row and the e0_RD, e0_RS, e1_RD, e1_RS rows (statistic, value, n, loci, note) of a painting:
    the fitted emissions; NA (with na_note) for a skipped or collapsed painting."""
    Y, D = p.trio.S, p.trio.D
    if p.fit is None:
        model = "emission model: nothing fitted (the painting was skipped)"
    else:
        model = "emission model: both states' emissions fitted by EM"
    rows = [("emission", "em", n_inf, inf_loci, model)]
    for state, k in (("e0", RD), ("e0", RS), ("e1", RD), ("e1", RS)):
        name = f"{state}_{PAINT_CLASSES[k]}"
        what = (f"state 0's ({Y}-like) emission probability of {PAINT_CLASSES[k]}" if state == "e0"
                else f"state 1's ({D}-like) emission probability of {PAINT_CLASSES[k]}")
        if na or p.fit is None:
            rows.append((name, math.nan, n_inf, inf_loci, na_note))
        else:
            rows.append((name, float(getattr(p.fit, state)[k]), n_inf, inf_loci, f"{what}, fitted by EM"))
    return rows


def _resolution(cat: Catalog, p: Painting, scope: str, loci: str, f: Findings) -> list[StatRow]:
    """paint_resolution rows: the share of R's painted loci where S and D are resolved, from R's own calls."""
    X, Y, D = p.trio.R, p.trio.S, p.trio.D
    if not len(p.idx):
        return []
    fracs = {}
    for g in (Y, D):
        fracs[g] = int(np.count_nonzero(cat.pcls(X, g)[p.idx] >= 0)) / len(p.idx)
    return resolution("paint_resolution", scope, fracs, len(p.idx), loci, f, PAINTING)


def _seqid_rows(cat: Catalog, p: Painting, scope: str, inf_loci: str, all_loci: str) -> list[StatRow]:
    """paint_seqid rows, one set per R sequence holding painted loci (fai order; stratum = its name)."""
    X = p.trio.R
    seq, _, _ = _sites(cat, X, p.idx)
    rows = []
    for s in np.unique(seq).tolist():
        m = seq == s
        c = p.cls[m]
        n_rd, n_rs = int(np.count_nonzero(c == RD)), int(np.count_nonzero(c == RS))
        name, n = cat.seqids[X][s], int(np.count_nonzero(m))
        share = n_rd / (n_rd + n_rs) if n_rd + n_rs else math.nan
        rows += [StatRow("pav", "paint_seqid", scope, name, "posterior_mean", float(p.post[m].mean()), n=n,
                         loci=all_loci, note="mean posterior over the sequence's painted loci"),
                 StatRow("pav", "paint_seqid", scope, name, "n_RD", n_rd, n=n, loci=inf_loci),
                 StatRow("pav", "paint_seqid", scope, name, "n_RS", n_rs, n=n, loci=inf_loci),
                 StatRow("pav", "paint_seqid", scope, name, "class_share", share, n=n_rd + n_rs,
                         loci=inf_loci, note="RD / (RD + RS)")]
    return rows


def _kb(v: float | None) -> str:
    return "nan" if v is None or not math.isfinite(v) else f"{v / 1e3:,.1f} kb"


def _n_loci(k: int) -> str:
    return f"{k} locus" if k == 1 else f"{k} loci"


def _floor_line(X: str, n_tracts: int, n_decoy: int) -> str | None:
    """The findings line on the +1 corrections' floors, naming each floor above Q_TRACT with its own reason:
    q >= 1/n_tracts (too few tracts), p_decoy >= 1/(1 + n_decoy) (too few decoy tracts). None when neither
    floor is above Q_TRACT, or without a tract."""
    if not n_tracts:
        return None
    parts = []
    if 1 / n_tracts > Q_TRACT:
        parts.append(f"q >= 1/{n_tracts} = {1 / n_tracts:.3g}: too few tracts for a target-decoy q to reach "
                     f"{Q_TRACT:g}")
    if 1 / (1 + n_decoy) > Q_TRACT:
        parts.append(f"p_decoy >= 1/{1 + n_decoy} = {1 / (1 + n_decoy):.3g}: too few decoy tracts for a "
                     f"per-tract p-value to reach {Q_TRACT:g}")
    if not parts:
        return None
    return (f"Painting of {X}: {n_tracts} tract(s) and {n_decoy} decoy tract(s), and the +1 corrections keep "
            + "; ".join(parts) + ", whatever the evidence")


def _painting_lines(cat: Catalog, p: Painting, f: Findings, post_mean: float, dec_mean: float, excess: float,
                    n_q05: float, counts: dict[int, int], spacing: tuple[float, float],
                    callable_: tuple[tuple[int | None, float | None], tuple[int | None, float | None]],
                    n_inv: int, run_bp: float) -> None:
    X, Y, D = p.trio.R, p.trio.S, p.trio.D
    fit, decoy = p.fit, p.scale is not None
    who = f"the {len(p.idx)} loci {X} discovers, each classed by {X}'s own calls"
    line = (f"Painting of {X} (sister {Y}, donor {D}; {who}; {counts[RD] + counts[RS]} informative: "
            f"{counts[RD]} RD with {D} present and {Y} E, {counts[RS]} RS with {Y} present and {D} E): "
            f"{D}-like posterior mean {post_mean:.3f}")
    if decoy:
        line += (f" against {dec_mean:.3f} in the decoy {Y} (excess {excess:+.3f}); {len(p.tracts)} "
                 f"tract(s), {len(p.decoy_tracts)} in the decoy, {n_q05} at q <= {Q_TRACT:g}")
    else:
        line += f"; {len(p.tracts)} tract(s), no q-values"
    if p.collapsed:
        line += f"; pi at its floor {PI_MIN:g}"
    else:
        line += (f"; switch scale L (mean_tract_bp) {_kb(fit.L)}, mean {D}-like run L / (1 - pi) "
                 f"(mean_donor_run_bp) {_kb(run_bp)}, pi {fit.pi:.3f}")
    f.add(PAINTING, line)
    if p.collapsed:
        f.add(PAINTING, f"Painting of {X}: the fit collapsed (pi at its floor {PI_MIN:g}, no tract): no "
                        f"{D}-like ancestry is detectable in this trio. donor_state, L (mean_tract_bp), "
                        "mean_donor_run_bp and the shortest callable tract are NA: state 1 holds no data, so "
                        "they would come from the pseudocounts")
    floor = _floor_line(X, len(p.tracts), len(p.decoy_tracts)) if decoy else None
    if floor:
        f.add(PAINTING, floor)
    if not p.donor_state and not p.collapsed:
        f.warn(PAINTING, f"Painting of {X}: no {D}-like state (the fitted second state is not enriched for "
                         f"RD loci: e1[RD] {fit.e1[RD]:.3g} <= e0[RD] {fit.e0[RD]:.3g}), so no tract is "
                         f"called and the posterior is not a {D}-like probability")
    if not decoy:
        f.warn(PAINTING, f"Painting of {X}: no decoy ({_no_decoy(cat, p)}), so tract q-values and decoy "
                         "p-values are written none and the excess over the decoy is not estimated")
    if fit.pi >= PI_MAX:
        f.add(PAINTING, f"Painting of {X}: the fitted {D}-like fraction pi reached its cap PI_MAX = "
                        f"{PI_MAX:g} (a recipient more {D}-like than {Y}-like would be {D}'s sister: check "
                        "the roles)")
    f.add(PAINTING, f"Painting of {X}: {counts[OLD]} loci an outgroup carries in {X}'s own call (class old: "
                    f"older than the ingroup split), {counts[MISSING]} with {Y} or {D} unresolved, "
                    f"{counts[BOTH]} both and {counts[NEITHER]} neither enter the HMM as missing data")
    (lo_k, lo_bp), (hi_k, hi_bp) = callable_

    def span(k: int | None, bp: float | None) -> str:
        return f"not callable within {CALLABLE_MAX} loci" if k is None else f"{_n_loci(k)} (about {_kb(bp)})"
    head = (f"Painting resolution of {X}: informative loci every {_kb(spacing[0])} along a sequence (median; "
            f"mean {_kb(spacing[1])}); ")
    if p.collapsed:
        tail = "shortest callable tract NA (the fit collapsed)"
    elif p.donor_state:
        tail = (f"shortest callable tract {span(lo_k, lo_bp)} without flanks, {span(hi_k, hi_bp)} between "
                f"{CALLABLE_FLANK} RS loci on each side ({CALLABLE_FLANK} is a heuristic flank)")
    else:
        tail = f"no shortest callable tract without a {D}-like state"
    f.add(PAINTING, head + tail)
    f.add(PAINTING, f"Painting of {X}: {n_inv} of the {len(p.idx)} painted loci lie inside inverted blocks "
                    f"of {X} against {Y} or {D} (posterior column inv; tract column inv_frac); the fit "
                    "ignores them (flags only)")


def _dating_lines(cat: Catalog, p: Painting, f: Findings) -> None:
    """The Dating lines: the median youngest K over tracts and the tracts holding an identical-LTR RD
    insertion, with their caveats; the youngest_u95 line only with identity.tsv's self rows."""
    X, D = p.trio.R, p.trio.D
    gx = cat.g(X)
    table = _tract_table(p, cat)
    if not np.isfinite(cat.k[p.idx[p.cls == RD], gx]).any():
        f.add(DATING, f"Painting of {X}: youngest insertion per tract skipped: no K2P for {X}'s copies at "
                      "its RD loci (an annotation without LTR identity)")
    elif not p.tracts:
        f.add(DATING, f"Painting of {X}: youngest insertion per tract: no tract called")
    else:
        ks = [v["youngest_k"] for v in table if v["youngest_k"] is not None]
        n_ident = sum(v["n_identical"] > 0 for v in table)
        med = f"median {np.median(ks):.4g} over the {len(ks)} of {len(table)} tracts with a K" if ks else (
            f"none of the {len(table)} tracts has a K at its RD loci")
        f.add(DATING, f"Youngest insertion per tract of {X} (donor {D}; the smallest K2P of {X}'s copy over "
                      f"each tract's RD loci, from the annotation): {med}; {n_ident} tract(s) hold an "
                      "identical-LTR (K = 0) RD insertion. The admixture is younger than any insertion it "
                      "transferred. Caveats: K = 0 is not age 0 (a few hundred LTR sites without a "
                      "difference bound the age only from above), and LTR-LTR gene conversion, frequent in "
                      "large plant genomes such as maize (Cossu et al. 2017, doi:10.1093/gbe/evx260), "
                      "makes K too young, so the bound can be too recent")
    if p.self_counts is not None:
        us = [v["youngest_u95"] for v in table if v["youngest_u95"] is not None]
        text = (f"median {np.median(us):.4g} over the {len(us)} of {len(table)} tracts with a self row"
                if us else f"no self row at the RD loci of the {len(table)} tract(s)")
        f.add(DATING, f"youngest_u95 per tract of {X} (identity.tsv's self rows, exact counts: the smallest "
                      "one-sided 95% Poisson upper bound of LTR-LTR divergence over each tract's RD loci, "
                      f"from {X}'s copy, else {D}'s): {text}")
