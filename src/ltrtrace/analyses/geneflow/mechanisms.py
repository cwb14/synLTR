"""Mechanisms and missingness: what genome B holds at genome A's elements, by superfamily, family
and age; how many absences are shown ancestral or derived; absence events; the reasons calls are missing;
and the parsimony (Dollo) placement of insertions and losses on the rooted PAV tree.

Sources (catalog.py):
  pair files  each ordered pair file pairs/<A>__<B>.calls.tsv.gz (A annotated) is read once, streamed
              (`pair_calls`, through calls.records.stream_calls): one call per element of A's annotation, every
              element (not only a locus's first member), with the evidence columns the catalog never parses. The evidence exists only per pair,
              so the tables states, absence, events and missingness have scope "A>B". A call's superfamily,
              family and age bin come from A's element (work/elements/<A>.tsv): AGE_LABELS by its LTR-LTR K2P,
              an element without K counted in the age "all" only; "." an unknown superfamily or family.
  merged      loci.tsv's states (`Catalog.cls`, any discovery), for the Dollo placement only.

O and D split by their junctions. w = SLOP_TSD_MULT x L, two TSD lengths (L from the call's tsd_a, TSD_DEFAULT
without one): the classifier's own E slop. dl and dr are the balanced overhangs of the gap past the element's
left and right ends, negative inside it (classify.py, rule 4). A gap end lies inside when it is more than w
inside the element (dl < -w or dr < -w). The stretch it cuts off (the element's first -dl or last -dr bases) is
a shown remnant of the element's own sequence when both hold:
  - it is shorter than MIN_SEG: rule 4 is reached only when fewer than MIN_SEG own bases align, so a longer
    stretch cannot be aligned own sequence (the gap is read around the element's midpoint, and an isolated hit
    can bound it deep inside);
  - no base of it belongs to an annotated child (work/elements/<A>.tsv, the `parent` links): what B keeps
    there may be the child's sequence (the child, or another copy of it), which does not show the element's.
  D_inside           a gap end inside with a shown remnant: B kept part of the element, so it carried it: a
                     derived absence, shown.
  D_inside_unshown   a gap end inside, none with a shown remnant: undetermined.
  D_beyond           any other D: undetermined (a deletion in B, of the element or of an empty site, or an
                     A-only insertion of a larger unit the annotation lacks).
  O_inside           as D_inside, the rest replaced by B-only sequence: a derived absence, shown.
  O_inside_unshown   as D_inside_unshown, with B-only sequence: undetermined.
  O_site             no gap end inside, |dl| <= w and |dr| <= w: B-only sequence in the element's place
                     (small_ins included): an independent insertion at the site, or a replacement.
  O_complex          any other O (a gap end more than w past the element, none inside it): a complex event.
  "Shown" assumes exact annotated ends, as the classifier's E/D split does: an end misplaced more than w into
  the flank gives the same signal. A D or O passed from its host (rule 0, flag inherit) has no junction of its
  own (dl, dr "."): the host's gap covers the nested element and runs past both its ends, and the host's history
  is not the nested element's. It counts as D_beyond or O_complex, never as shown derived or at the site.
Route constants read (base.constants): SLOP_TSD_MULT and TSD_DEFAULT (w), and MIN_SEG (the longest remnant that
can be aligned own sequence). They are read when the analysis runs, so they must equal those of the route that made
the run.

Tables of mechanisms.tsv (MECH_COLUMNS, "all" a rolled-up column; written sorted by every column but value):
  states       per (superfamily, family, age bin), with the rollups (superfamily, all) and (all, all), each
               per age bin and "all": n_P n_S n_T n_E n_D_inside n_D_inside_unshown n_D_beyond n_O_site
               n_O_inside n_O_inside_unshown n_O_complex n_G n_U, for the strata holding a call. A named family
               lists its nonzero counts only (a metric it lacks is 0): most families hold one or two elements,
               so with every metric the table is mostly zeros. Every
               rollup (family all) lists every metric, and (all, all, all) is always there.
  absence      per age bin holding a call, and "all" (superfamily, family all): n_absent = E + D + O, and
               frac_E, frac_D_inside, frac_O_inside, frac_D_inside_unshown, frac_O_inside_unshown,
               frac_D_beyond, frac_O_site, frac_O_complex of it. E is ancestral, shown (an empty site:
               retrotransposons do not excise precisely); D_inside and O_inside are derived, shown; the other
               five are undetermined (D and O say nothing about whether the insertion happened).
               So 1 - frac_E bounds the derived share of absences from above and frac_D_inside + frac_O_inside
               from below: the bounds bracket it, never measure it.
  events       (all, all, all): n_E_single n_E_cluster n_E_cluster_elements n_D_single n_D_multi
               n_D_multi_elements. An E event is an E call without host_empty with its descendants (by the
               `parent` links) that are E with host_empty; a D event a D call without inherit with its
               descendants that are D with inherit; E calls flagged adjacent in one window (`wid`) form one
               event, and D calls likewise. An event of more than one element, or an adjacent one, is
               multi-element: for E a nested cluster absent (inserted) as a unit, for D a multi-element
               deletion whose own origin stays open as for any D. The pair files keep no gap ids, so two
               adjacent clusters in one window count as one event.
  missingness  (all, all, all): U_<reason> for each U reason found (the window's status, rules 1-4 and
               unanchored from rule 5a; "unknown" for none) and G, as fractions of the pair's calls.
  branches     the Dollo placement, below: per branch (scope), per superfamily (loci.tsv's) and "all":
               n_insertion n_solo n_truncation n_deletion; scope "incongruent" (superfamily all): n_loci.

Dollo placement (needs an outgroup and the tree rooted on it by relations; else one "skipped" line), on merged
states, each distinct state pattern placed once: present = the genomes of class 1 (P/S/T); the clade is their
MRCA. An E genome in the clade makes the locus incongruent (E is a site that never held the insertion, and
retrotransposons do not excise precisely): counted, not placed. Otherwise the insertion goes on the branch
above the clade (before_root when the clade is the whole tree), and for each of S, T and D one event on each
maximal subclade of the clade whose genomes are all that state or missing, with at least one of the state.
"Missing" is D/O/G/U for S and T events and O/G/U for D events; a genome outside the clade is never placed
(it is missing). A branch is named by its leaves joined by "+" in sheet order.

Stats rows (section "Mechanisms") mirror the tables: tests mechanisms (absence), absence_events, missingness
(plus the pair's resolved fraction of B) and branches (plus branches_resolution, each genome's resolved
fraction over all loci, through stats.report.resolution). Findings: per pair the bounds by age, the multi-element
split and the main U reasons; the insertions and losses per branch.
"""
from __future__ import annotations

import math
from collections import Counter
from collections.abc import Iterator
from pathlib import Path

import numpy as np

from ltrtrace.analyses.geneflow.catalog import AGE_LABELS, STATE_CODE, AnalysisError, Catalog, age_bin
from ltrtrace.analyses.geneflow.relations import Node
from ltrtrace.base.constants import MIN_SEG, SLOP_TSD_MULT, TSD_DEFAULT
from ltrtrace.calls.records import CALL_COLUMNS, stream_calls, tsv_field
from ltrtrace.calls.states import PRESENT
from ltrtrace.io.elements import read_elements
from ltrtrace.io.samples import PAIR_SEP
from ltrtrace.stats.report import Findings, StatRow, fmt, resolution

SECTION = "Mechanisms"  # findings section (stats.report.SECTIONS)
MECH_COLUMNS = ("table", "scope", "superfamily", "family", "age_bin", "metric", "value")
STATE_METRICS = ("n_P", "n_S", "n_T", "n_E", "n_D_inside", "n_D_inside_unshown", "n_D_beyond", "n_O_site",
                 "n_O_inside", "n_O_inside_unshown", "n_O_complex", "n_G", "n_U")
# n_absent's parts, in this order: ancestral shown, derived shown, undetermined
ABSENT = ("n_E", "n_D_inside", "n_O_inside", "n_D_inside_unshown", "n_O_inside_unshown", "n_D_beyond", "n_O_site",
          "n_O_complex")
ABSENCE_METRICS = ("frac_E", "frac_D_inside", "frac_O_inside", "frac_D_inside_unshown", "frac_O_inside_unshown",
                   "frac_D_beyond", "frac_O_site", "frac_O_complex")
EVENT_METRICS = ("n_E_single", "n_E_cluster", "n_E_cluster_elements", "n_D_single", "n_D_multi",
                 "n_D_multi_elements")
BRANCH_METRICS = ("n_insertion", "n_solo", "n_truncation", "n_deletion")
ALL = "all"

ABSENCE_NOTES = {
    "n_absent": "E + D + O calls",
    "frac_E": "ancestral, shown: an empty site (retrotransposons do not excise precisely)",
    "frac_D_inside": f"derived, shown: a gap end more than two TSD lengths inside the element cuts off a remnant "
                     f"of its own sequence (under {MIN_SEG} bp, no annotated child in it), so B carried the "
                     "element; assumes exact annotated ends, as the E/D split does",
    "frac_O_inside": f"derived, shown: as D_inside (a remnant of the element's own sequence, under {MIN_SEG} bp, "
                     "no annotated child in it), the rest replaced by B-only sequence",
    "frac_D_inside_unshown": f"not determined: a gap end more than two TSD lengths inside the element, but the "
                             f"stretch it cuts off is {MIN_SEG} bp or longer (more own sequence than rule 4 lets "
                             "align) or holds an annotated child, so it does not show the element's own sequence",
    "frac_O_inside_unshown": "not determined: as D_inside_unshown, with B-only sequence in the gap",
    "frac_D_beyond": "not determined: a deletion in B (of the element or of an empty site) or an A-only "
                     "insertion of a larger unit the annotation lacks; D inherited from a host included",
    "frac_O_site": "not determined: B-only sequence in the element's place (an independent insertion at the "
                   "site, or a replacement)",
    "frac_O_complex": "not determined: B-only sequence with a gap end more than two TSD lengths past the "
                      "element's ends and none inside it, or O inherited from a host: a complex event",
}
EVENT_NOTES = {
    "n_E_single": "an E call alone (no host_empty link, not adjacent)",
    "n_E_cluster": "a nested cluster absent (inserted) as a unit: an E host with its host_empty descendants, "
                   "or E calls flagged adjacent in one window",
    "n_E_cluster_elements": "the elements of those clusters",
    "n_D_single": "a D call alone (no inherit link, not adjacent)",
    "n_D_multi": "a multi-element deletion: a D host with its inherit descendants, or D calls flagged adjacent "
                 "in one window; its own origin stays open, as for any D",
    "n_D_multi_elements": "the elements of those deletions",
}
DOLLO_LOCI = "disc=any;merged states;every locus"
DOLLO_NOTE = ("Dollo parsimony on the rooted PAV tree: the insertion on the branch above the MRCA of the "
              "present genomes, S, T and D losses on maximal subclades; a missing state never constrains")
INCONGRUENT_NOTE = "an E genome inside the clade of the present genomes: counted, not placed"

_SIMPLE = {x: f"n_{x}" for x in STATE_CODE if x not in ("D", "O")}  # states not split by their junction
_MI = {m: i for i, m in enumerate(STATE_METRICS)}
_LINKED = {"E": "host_empty", "D": "inherit"}  # the flag that ties a call to its host's event
_PRESENT = frozenset(STATE_CODE[x] for x in PRESENT)
_EMPTY = STATE_CODE["E"]
_MISS_ST = frozenset((0, *(STATE_CODE[x] for x in "DOGU")))  # missing for S and T events (0: no call)
_MISS_D = frozenset((0, *(STATE_CODE[x] for x in "OGU")))    # missing for D events
_LOSSES = (("n_solo", STATE_CODE["S"], _MISS_ST), ("n_truncation", STATE_CODE["T"], _MISS_ST),
           ("n_deletion", STATE_CODE["D"], _MISS_D))
_MISMATCH = "the run's element table and pair files do not match (rerun ltrtrace genotype)"


# ---------------------------------------------------------------- calls


def pair_calls(path: Path) -> Iterator[dict[str, str]]:
    """The calls of one pair file (pairs/<A>__<B>.calls.tsv.gz), streamed (calls.records.stream_calls): one dict
    per line, CALL_COLUMNS -> the field as written ("." for none). An AnalysisError naming the file on a header
    that is not CALL_COLUMNS, a line of another width, or a file that cannot be read (missing, corrupt or
    truncated)."""
    path = Path(path)
    try:
        for _, f in stream_calls(path):
            yield dict(zip(CALL_COLUMNS, f))
    except (OSError, ValueError) as exc:
        raise AnalysisError(f"{path}: cannot read the pair file ({type(exc).__name__}: {exc})") from exc


def _w(tsd: str | None) -> int:
    """The classifier's E slop: SLOP_TSD_MULT x the TSD length (TSD_DEFAULT without a TSD)."""
    return SLOP_TSD_MULT * (len(tsd) if tsd else TSD_DEFAULT)


_NO_CHILD = (math.inf, math.inf)  # `own` of an element without annotated children


def _inside(dl: int, dr: int, w: int, own: tuple[float, float]) -> tuple[bool, bool]:
    """(a gap end lies more than w inside the element, one such end cuts off a shown remnant): fewer than
    MIN_SEG bases, all the element's own (no deeper than own[0] at the left end, own[1] at the right)."""
    left, right = dl < -w, dr < -w
    shown = (left and -dl < MIN_SEG and -dl <= own[0]) or (right and -dr < MIN_SEG and -dr <= own[1])
    return left or right, shown


def o_kind(dl: int, dr: int, tsd: str | None, own: tuple[float, float] = _NO_CHILD) -> str:
    """With a gap end more than w inside the element (dl < -w or dr < -w), whatever the other end does:
    "O_inside" if one such end cuts off a shown remnant (under MIN_SEG bp of the element's own sequence: B
    kept part of the element, so it carried it), else "O_inside_unshown". Otherwise "O_site" if both gap ends
    lie within w of the element's ends (|dl| <= w and |dr| <= w), else "O_complex". w = SLOP_TSD_MULT x L
    (L = len(tsd), or TSD_DEFAULT). own: the element's own bases at its left and right ends before its first
    annotated child there (_NO_CHILD without children)."""
    w = _w(tsd)
    inside, shown = _inside(dl, dr, w, own)
    if inside:
        return "O_inside" if shown else "O_inside_unshown"
    return "O_site" if abs(dl) <= w and abs(dr) <= w else "O_complex"


def d_kind(dl: int, dr: int, tsd: str | None, own: tuple[float, float] = _NO_CHILD) -> str:
    """"D_inside" if a gap end more than w inside the element (dl < -w or dr < -w) cuts off a shown remnant
    (B kept part of the element, so it carried it: a derived absence), "D_inside_unshown" if gap ends lie
    inside but none does, else "D_beyond". w and own as in o_kind."""
    inside, shown = _inside(dl, dr, _w(tsd), own)
    return "D_inside" if shown else "D_inside_unshown" if inside else "D_beyond"


def _metric(c: dict[str, str], path: Path, line: int, own: tuple[float, float]) -> int:
    """The STATE_METRICS index of call c: P S T E G U as they are; D by d_kind and O by o_kind on its own
    junction (own: its element's, as in o_kind), D_beyond or O_complex without one (inherited from a host: dl,
    dr ".")."""
    st = c["state"]
    m = _SIMPLE.get(st)
    if m is not None:
        return _MI[m]
    if st not in ("D", "O"):
        raise AnalysisError(f"{path}: line {line}: unknown state {st!r}")
    if c["dl"] == "." or c["dr"] == ".":
        return _MI["n_D_beyond" if st == "D" else "n_O_complex"]
    try:
        dl, dr = int(c["dl"]), int(c["dr"])
    except ValueError:
        raise AnalysisError(f"{path}: line {line}: dl, dr ({c['dl']!r}, {c['dr']!r}) are not integers") from None
    tsd = None if c["tsd_a"] == "." else c["tsd_a"]
    return _MI["n_" + (d_kind(dl, dr, tsd, own) if st == "D" else o_kind(dl, dr, tsd, own))]


# ---------------------------------------------------------------- per pair


_Info = tuple[dict[str, tuple[str, str, int]], dict[str, tuple[int, int]]]  # _element_info's (keys, own)


def _element_info(cat: Catalog, a: str) -> _Info:
    """Genome a's element table (work/elements/<a>.tsv) as (keys, own). keys: eid -> (superfamily, family, age
    bin code) of every element; "." an unknown superfamily or family, -1 an element without K. own: eid ->
    (its own bases at its left end before its first annotated child, at its right end after its last one), for
    the elements with children (direct children, by the `parent` links, as the classifier's holes)."""
    if not cat.annotated.get(a):
        raise AnalysisError(f"{a!r} is not an annotated genome of the run, so it has no pair files")
    path = cat.outdir / "work" / "elements" / f"{a}.tsv"
    try:
        els = read_elements(path)
    except (OSError, ValueError, TypeError) as exc:
        raise AnalysisError(f"{path}: cannot read the element table of {a} ({exc})") from exc
    ages = age_bin(np.array([math.nan if e.k is None else e.k for e in els], dtype=np.float64)).tolist()
    keys = {e.eid: (e.superfamily or ".", e.family or ".", age) for e, age in zip(els, ages, strict=True)}
    span = {e.eid: (e.start, e.end) for e in els}
    own: dict[str, tuple[int, int]] = {}
    for e in els:
        host = span.get(e.parent) if e.parent is not None else None
        if host is not None:
            left, right = own.get(e.parent, (host[1] - host[0] + 1,) * 2)
            own[e.parent] = (min(left, e.start - host[0]), min(right, host[1] - e.end))
    return keys, own


def _pair(cat: Catalog, a: str, b: str, info: _Info, f: Findings) -> tuple[list[tuple], list[StatRow]]:
    """The states, absence, events and missingness tables of the ordered pair a>b (module docstring), read from
    pairs/<a>__<b>.calls.tsv.gz, with a's element table read (`info`), and their stats rows: (rows of
    MECH_COLUMNS, StatRows). Its findings lines are added to f."""
    keys, own = info
    cat.g(b)
    if b == a:
        raise AnalysisError(f"a pair needs two genomes, not {a} twice")
    scope = f"{a}>{b}"
    path = cat.outdir / "pairs" / f"{a}{PAIR_SEP}{b}.calls.tsv.gz"
    counts: dict[tuple[str, str, int], list[int]] = {}
    reasons: Counter[str] = Counter()
    events = _Events()
    seen: set[str] = set()
    for line, c in enumerate(pair_calls(path), 2):
        eid = c["eid"]
        key = keys.get(eid)
        if key is None:
            raise AnalysisError(f"{path}: line {line}: {eid} is not an element of {a} (work/elements/{a}.tsv): "
                                f"{_MISMATCH}")
        if eid in seen:
            raise AnalysisError(f"{path}: line {line}: {eid} is called more than once")
        seen.add(eid)
        acc = counts.get(key)
        if acc is None:
            acc = counts[key] = [0] * len(STATE_METRICS)
        acc[_metric(c, path, line, own.get(eid, _NO_CHILD))] += 1
        st = c["state"]
        if st == "U":
            reasons["unknown" if c["reason"] == "." else c["reason"]] += 1
        elif st in _LINKED:
            events.add(st, eid, c["parent"], c["flags"], c["wid"])
    if len(seen) != len(keys):
        raise AnalysisError(f"{path}: {len(keys) - len(seen)} of the {len(keys)} elements of {a} have no call: "
                            f"{_MISMATCH}")
    n = len(seen)
    agg = _rollup(counts)
    rows = [("states", scope, sf, fam, age, m, x) for (sf, fam, age), v in agg.items()
            for m, x in zip(STATE_METRICS, v) if x or fam == ALL]   # a named family: nonzero counts only
    out_rows, stat, absence_line = _absence(agg, a, b, scope)
    rows += out_rows
    ev = events.counts()
    loci = f"disc={a};{b} E or D in {a}'s calls (every element of {a})"
    for m in EVENT_METRICS:
        rows.append(("events", scope, ALL, ALL, ALL, m, ev[m]))
        stat.append(StatRow("pav", "absence_events", scope, ALL, m, ev[m], n=ev["calls_" + m[2]], loci=loci,
                            note=EVENT_NOTES[m]))
    total = agg[(ALL, ALL, ALL)]
    miss = {f"U_{r}": k / n if n else math.nan for r, k in sorted(reasons.items())}
    miss["G"] = total[_MI["n_G"]] / n if n else math.nan
    loci = f"disc={a};every element of {a}, its call into {b}"
    for m, x in miss.items():
        rows.append(("missingness", scope, ALL, ALL, ALL, m, x))
        stat.append(StatRow("pav", "missingness", scope, ALL, m, x, n=n, loci=loci,
                            note="fraction of the pair's calls"))
    if not n:
        f.add(SECTION, f"{scope}: no call (the element table of {a} is empty)")
        return rows, stat
    resolved = sum(total[_MI[m]] for m in ("n_P", "n_S", "n_T", "n_E"))
    stat += resolution("missingness", scope, {b: resolved / n}, n, loci, f, SECTION)
    f.add(SECTION, absence_line)
    f.add(SECTION, f"{scope} absence events (by parent links and the adjacent flag): E {ev['events_E']:,} = "
                   f"{ev['n_E_single']:,} single + {ev['n_E_cluster']:,} nested clusters absent as a unit "
                   f"({ev['n_E_cluster_elements']:,} elements); D {ev['events_D']:,} = {ev['n_D_single']:,} single "
                   f"+ {ev['n_D_multi']:,} multi-element deletions ({ev['n_D_multi_elements']:,} elements; their "
                   "origin stays open, as for any D)")
    top = ", ".join(f"{r} {k / n:.3f}" for r, k in sorted(reasons.items(), key=lambda kv: (-kv[1], kv[0]))[:3])
    f.add(SECTION, f"{scope} missing: U {total[_MI['n_U']] / n:.3f} of the {n:,} calls"
          + (f" ({top})" if top else "") + f", G {miss['G']:.4f}")
    return rows, stat


def _rollup(counts: dict[tuple[str, str, int], list[int]]) -> dict[tuple[str, str, str], list[int]]:
    """(superfamily, family, age label) -> STATE_METRICS counts, for every stratum holding a call and its
    rollups (superfamily, all) and (all, all), per age bin and "all"; (all, all, all) always."""
    agg: dict[tuple[str, str, str], list[int]] = {(ALL, ALL, ALL): [0] * len(STATE_METRICS)}
    for (sf, fam, age), v in counts.items():
        ages = (ALL,) if age < 0 else (AGE_LABELS[age], ALL)
        for s, fm in dict.fromkeys(((sf, fam), (sf, ALL), (ALL, ALL))):  # a family named "all": counted once
            for ag in ages:
                acc = agg.setdefault((s, fm, ag), [0] * len(STATE_METRICS))
                for i, x in enumerate(v):
                    acc[i] += x
    return agg


def _absence(agg: dict[tuple[str, str, str], list[int]], a: str, b: str, scope: str
             ) -> tuple[list[tuple], list[StatRow], str]:
    """The absence table per age bin holding a call and "all", its stats rows (test mechanisms) and the
    findings line with the two bounds."""
    rows: list[tuple] = []
    stat: list[StatRow] = []
    loci = f"disc={a};{b} E, D or O in {a}'s calls (every element of {a})"
    frac: dict[str, dict[str, float]] = {}
    n_abs: dict[str, int] = {}
    for age in (*AGE_LABELS, ALL):
        if (ALL, ALL, age) not in agg:
            continue
        v = agg[(ALL, ALL, age)]
        parts = [v[_MI[m]] for m in ABSENT]
        n_abs[age] = k = sum(parts)
        frac[age] = {m: p / k if k else math.nan for m, p in zip(ABSENCE_METRICS, parts)}
        rows.append(("absence", scope, ALL, ALL, age, "n_absent", k))
        stat.append(StatRow("pav", "mechanisms", scope, age, "n_absent", k, n=sum(v), loci=loci,
                            note=ABSENCE_NOTES["n_absent"]))
        for m, x in frac[age].items():
            rows.append(("absence", scope, ALL, ALL, age, m, x))
            stat.append(StatRow("pav", "mechanisms", scope, age, m, x, n=k, loci=loci, note=ABSENCE_NOTES[m]))
    bins = [age for age in AGE_LABELS if age in frac]
    if bins:
        by_age = ", ".join(f"{age} {n_abs[age]:,} / {fmt(1 - frac[age]['frac_E'], '.3g')} / "
                           f"{fmt(_lower(frac[age]), '.3g')}" for age in bins)
        rest = n_abs[ALL] - sum(n_abs[age] for age in bins)
        if rest:
            by_age += f" ({rest:,} absences without K2P count in all only)"
    else:
        by_age = f"skipped (no K2P for the elements of {a})"
    fa = frac[ALL]
    calls = f"the {sum(agg[(ALL, ALL, ALL)]):,} calls of {a}'s elements in {b}"
    if not n_abs[ALL]:
        return rows, stat, f"{scope}: no absence (E, D or O) among {calls}"
    g = {m: fmt(x, ".3g") for m, x in fa.items()}  # 3 significant digits: a nonzero fraction never prints as 0
    return rows, stat, (
        f"{scope}: {n_abs[ALL]:,} absences (E, D or O) among {calls}: absences not shown ancestral = 1 - frac_E "
        f"= {fmt(1 - fa['frac_E'], '.3g')} (an upper bound on derived absences); shown derived = frac_D_inside + "
        f"frac_O_inside = {g['frac_D_inside']} + {g['frac_O_inside']} = {fmt(_lower(fa), '.3g')} (a lower "
        f"bound: a gap end inside the element cuts off a remnant of its own sequence, under {MIN_SEG} bp with no "
        f"annotated child, so {b} carried it); by age, absences / upper / lower: {by_age}. A gap end inside "
        f"without such a remnant (D {g['frac_D_inside_unshown']}, O {g['frac_O_inside_unshown']}), "
        f"D beyond the element ({g['frac_D_beyond']}), O at the site ({g['frac_O_site']}) and complex "
        f"O ({g['frac_O_complex']}) do not say whether {b} carried the insertion, so the "
        "bounds bracket the derived share; they do not measure it")


def _lower(fr: dict[str, float]) -> float:
    """The lower bound on the derived share of absences: those shown derived, D_inside + O_inside."""
    return fr["frac_D_inside"] + fr["frac_O_inside"]


class _Events:
    """Absence events of one pair (module docstring): union-find over its E and D calls, joined by the
    host_empty / inherit links to their hosts and by the adjacent flag within a window."""

    def __init__(self) -> None:
        self.up: list[int] = []                          # union-find parent of each call
        self.state: list[str] = []                       # "E" or "D"
        self.adjacent: list[bool] = []
        self.node: dict[tuple[str, str], int] = {}       # (state, eid) -> call
        self.links: list[tuple[int, str]] = []           # (call, its host's eid): joined once all are read
        self.first: dict[tuple[str, str], int] = {}      # (state, wid) -> its first adjacent call

    def add(self, state: str, eid: str, parent: str, flags: str, wid: str) -> None:
        i = len(self.up)
        words = flags.split(",")
        self.up.append(i)
        self.state.append(state)
        self.adjacent.append("adjacent" in words)
        self.node[(state, eid)] = i
        if _LINKED[state] in words and parent != ".":
            self.links.append((i, parent))
        if self.adjacent[i]:
            self._union(i, self.first.setdefault((state, wid), i))

    def counts(self) -> dict[str, int]:
        """EVENT_METRICS, plus events_E, events_D (events) and calls_E, calls_D (calls)."""
        for i, host in self.links:
            j = self.node.get((self.state[i], host))  # a host of another state (never in the engine): no link
            if j is not None:
                self._union(i, j)
        size: Counter[int] = Counter()
        joined: set[int] = set()
        for i in range(len(self.up)):
            r = self._find(i)
            size[r] += 1
            if self.adjacent[i]:
                joined.add(r)
        out = dict.fromkeys((*EVENT_METRICS, "events_E", "events_D", "calls_E", "calls_D"), 0)
        multi = {"E": ("n_E_cluster", "n_E_cluster_elements"), "D": ("n_D_multi", "n_D_multi_elements")}
        for r, k in size.items():
            st = self.state[r]
            out[f"events_{st}"] += 1
            out[f"calls_{st}"] += k
            if k > 1 or r in joined:
                out[multi[st][0]] += 1
                out[multi[st][1]] += k
            else:
                out[f"n_{st}_single"] += 1
        return out

    def _find(self, i: int) -> int:
        r = i
        while self.up[r] != r:
            r = self.up[r]
        while self.up[i] != r:
            self.up[i], i = r, self.up[i]
        return r

    def _union(self, i: int, j: int) -> None:
        ri, rj = self._find(i), self._find(j)
        if ri != rj:
            self.up[max(ri, rj)] = min(ri, rj)


# ---------------------------------------------------------------- Dollo placement


def _dollo(cat: Catalog, tree: Node) -> tuple[list[tuple], list[StatRow], dict]:
    """The branches table of the Dollo placement on `tree` (module docstring), merged states, and its stats
    rows, plus a summary for the findings: branches (pre-order) -> counts over every superfamily, and the
    placed, incongruent and total loci. A ValueError without an outgroup, or unless `tree` holds every genome
    once and is rooted between the outgroups and the ingroup (two children at its top, one of outgroups only), as
    relations roots it."""
    if not cat.outgroups:
        raise ValueError("dollo: needs an outgroup: without one the tree is unrooted, so no branch lies above "
                         "a clade")
    masks, kids = _clades(tree, cat)
    names = ["before_root", *(_label(cat, m) for m in masks[1:])]
    n = len(cat.lid)
    pats, inv = np.unique(cat.state, axis=0, return_inverse=True)
    sf_names, sf_code = np.unique(np.array(["." if s is None else s for s in cat.superfamily.tolist()],
                                           dtype=object).astype(str), return_inverse=True)
    width = max(len(sf_names), 1)
    keys, sizes = np.unique(np.asarray(inv).reshape(-1).astype(np.int64) * width
                            + np.asarray(sf_code).reshape(-1), return_counts=True)
    places = [_place(p, masks, kids) for p in pats.tolist()]
    table: dict[tuple[str, str], list[int]] = {}
    placed: Counter[str] = Counter()
    incongruent = 0
    for key, k in zip(keys.tolist(), sizes.tolist()):
        p, s = divmod(key, width)
        kind, events = places[p]
        if kind == "incongruent":
            incongruent += k
        elif kind == "placed":
            cols = tuple(dict.fromkeys((str(sf_names[s]), ALL)))  # a superfamily named "all": counted once
            for col in cols:
                placed[col] += k
            for node, metric in events:
                for col in cols:
                    table.setdefault((names[node], col), [0] * len(BRANCH_METRICS))[
                        BRANCH_METRICS.index(metric)] += k
    order = {}
    for name in names:
        order.setdefault(name, len(order))
    rows = [("branches", br, sf, ALL, ALL, m, table[(br, sf)][j])
            for br, sf in sorted(table, key=lambda t: (order[t[0]], t[1] != ALL, t[1]))
            for j, m in enumerate(BRANCH_METRICS)]
    stat = [StatRow("pav", "branches", br, sf, m, v, n=placed[sf], loci=DOLLO_LOCI, note=DOLLO_NOTE)
            for _, br, sf, _, _, m, v in rows]
    rows.append(("branches", "incongruent", ALL, ALL, ALL, "n_loci", incongruent))
    stat.append(StatRow("pav", "branches", "incongruent", ALL, "n_loci", incongruent, n=n, loci=DOLLO_LOCI,
                        note=INCONGRUENT_NOTE))
    summary = {"branches": [(br, table[(br, ALL)]) for br in order if (br, ALL) in table],
               "placed": placed[ALL], "incongruent": incongruent, "n": n}
    return rows, stat, summary


def _clades(tree: Node, cat: Catalog) -> tuple[list[int], list[list[int]]]:
    """The leaf masks (bit g: sheet genome g) and child indexes of `tree`'s nodes in pre-order (the root
    first); a ValueError unless its leaves are the run's genomes, each once, and it is rooted between the
    outgroups and the ingroup."""
    gi = {x: g for g, x in enumerate(cat.labels)}
    masks: list[int] = []
    kids: list[list[int]] = []

    def walk(x: Node) -> int:
        i = len(masks)
        masks.append(0)
        kids.append([])
        if not x.children:
            g = gi.get(x.name)
            if g is None:
                raise ValueError(f"dollo: the tree's leaf {x.name!r} is not a genome of the run "
                                 f"({', '.join(cat.labels)})")
            masks[i] = 1 << g
        for child in x.children:
            j = walk(child)
            if masks[i] & masks[j]:
                raise ValueError("dollo: a genome is a leaf of the tree more than once")
            masks[i] |= masks[j]
            kids[i].append(j)
        return i

    walk(tree)
    lacking = [x for g, x in enumerate(cat.labels) if not (masks[0] >> g) & 1]
    if lacking:
        raise ValueError(f"dollo: the tree lacks genome(s) {', '.join(lacking)}")
    out = sum(1 << gi[o] for o in cat.outgroups)
    if len(kids[0]) != 2 or not any(not (masks[j] & ~out) for j in kids[0]):
        raise ValueError("dollo: the tree must be rooted between the outgroups and the ingroup (two children at "
                         "its top, one of outgroups only), as relations roots it")
    return masks, kids


def _label(cat: Catalog, mask: int) -> str:
    """A branch's name: the leaves below it, joined by "+" in sheet order."""
    return "+".join(x for g, x in enumerate(cat.labels) if (mask >> g) & 1)


def _place(pat: list[int], masks: list[int], kids: list[list[int]]) -> tuple[str, list[tuple[int, str]]]:
    """One merged-state pattern (state codes in sheet order) on the tree: ("none", []) without a present
    genome, ("incongruent", []) with an E genome in the clade of the present ones, else ("placed", events):
    (node, BRANCH_METRICS name) for the insertion and each loss (module docstring)."""
    present = empty = 0
    for g, s in enumerate(pat):
        if s in _PRESENT:
            present |= 1 << g
        elif s == _EMPTY:
            empty |= 1 << g
    if not present:
        return "none", []
    c = min((i for i, m in enumerate(masks) if (m & present) == present), key=lambda i: masks[i].bit_count())
    clade = masks[c]
    if clade & empty:
        return "incongruent", []
    events = [(c, "n_insertion")]
    for metric, code, missing in _LOSSES:
        state = ok = 0
        for g, s in enumerate(pat):
            if s == code:
                state |= 1 << g
            elif s in missing:
                ok |= 1 << g
        if not (clade & state):
            continue
        ok |= state
        stack = [c]
        while stack:
            i = stack.pop()
            if masks[i] & ~ok:          # a genome here is neither the state nor missing: look below
                stack.extend(kids[i])
            elif masks[i] & state:      # all the state or missing, with the state: a maximal subclade
                events.append((i, metric))
    return "placed", events


# ---------------------------------------------------------------- the section


def mechanisms(cat: Catalog, tree: Node | None) -> tuple[list[tuple], list[StatRow], Findings]:
    """The mechanisms section (module docstring): the tables of every ordered pair from an
    annotated genome (sheet order), each pair file read once and each element table once, then the Dollo
    placement on `tree` (relations' tree: rooted on the outgroups; None when it was skipped), or one
    "skipped" line without an outgroup or a tree. Returns (rows of MECH_COLUMNS, StatRows, findings)."""
    f = Findings()
    rows: list[tuple] = []
    stat: list[StatRow] = []
    for a in cat.labels:
        if not cat.annotated[a]:
            continue
        info = _element_info(cat, a)
        for b in cat.labels:
            if b != a:
                r, s = _pair(cat, a, b, info, f)
                rows += r
                stat += s
    if not cat.outgroups:
        f.add(SECTION, "Dollo placement skipped: no outgroup (without one the tree is unrooted, so no branch "
                       "lies above a clade)")
    elif tree is None:
        f.add(SECTION, "Dollo placement skipped: no tree (Relationships says why)")
    else:
        r, s, summary = _dollo(cat, tree)
        rows += r
        stat += s
        per = "; ".join(f"{br} {'/'.join(str(x) for x in v)}" for br, v in summary["branches"]) or "none"
        f.add(SECTION, f"Dollo placement on the rooted PAV tree (merged states, any discovery; a missing state "
                       f"never constrains): {summary['placed']:,} of {summary['n']:,} loci placed, "
                       f"{summary['incongruent']:,} incongruent ({INCONGRUENT_NOTE}); insertions / solo-LTR "
                       f"formations / truncations / deletions per branch: {per}")
        if summary["n"]:
            fracs = {g: int(np.count_nonzero(cat.cls(g) >= 0)) / summary["n"] for g in cat.labels}
            stat += resolution("branches_resolution", ALL, fracs, summary["n"], DOLLO_LOCI, f, SECTION)
    return rows, stat, f


def write_mechanisms(rows: list[tuple], path: Path) -> None:
    """mechanisms.tsv: a header row of MECH_COLUMNS, then the rows sorted by (table, scope, superfamily,
    family, age_bin, metric); values as `tsv_field` writes them (nan ".")."""
    with open(path, "w", newline="") as fh:
        fh.write("\t".join(MECH_COLUMNS) + "\n")
        fh.writelines("\t".join(tsv_field(v) for v in r) + "\n" for r in sorted(rows, key=lambda r: r[:6]))
