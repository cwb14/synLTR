"""Trio roles and tests: which sister is the recipient, D_pav with its check strata,
the age gradient, the outgroup filter, the resolution of the sisters and the f4-ratio.

Counted loci: `relations.discover(cat, R, S, D, dist, "E")`, the one
discovery function. They are the loci the donor D discovers (with D unannotated, the eligible fallback
genomes), with D present, every outgroup E and R and S resolved. Every state is read from
the discovering genome's own calls at any tier (`Catalog.pcls_by`), never from loci.tsv's merged states,
where a member is P whatever the discoverer calls it and a member's own calls give the other genomes extra
chances to be resolved. gstar is each row's discovering genome.

D_pav = (n_RD - n_SD) / (n_RD + n_SD). n_RD counts the loci with R present (P/S/T) and S E; n_SD is the
mirror. It is positive when R shares more insertions with D than S does; lineage sorting alone gives 0 at
every age. Delete-one block jackknife (`stats.jackknife`) over the counted rows, on D's site as each row's
discoverer reads it (D's own element when D discovers), in blocks of `size` bp that never span two sequences.

Strata (rows D_pav, column stratum), each a subset of the primary count unless said otherwise:
  all                 the primary count
  K=0, 0<K<0.01, K>=0.01
                      by the LTR-LTR K2P of the discovering genome's copy (`age_bin`); no age rows when no
                      counted locus has a K ("age strata skipped: no K2P")
  high_tier           the discoverer's calls into R and S both at tier high. These are call tiers: a member's
                      loci.tsv tier "annot" never enters.
  nested_collapsed    the informative counted loci (R+D or S+D), one per nested cluster: the first of each
                      cluster (`cluster_keys`) in the discovering genome's (seq, start, -end) order
                      (informative loci are collapsed, so a host that is not informative hides no child)
  no_outgroup_filter  the same count without the outgroup rule (`discover(..., "none")`)
A stratum with no counted locus writes no rows; one whose counted loci are all uninformative writes D nan
(n 0) beside n_RD 0 and n_SD 0.

Age gradient (row D_pav_gradient): dD = D(K=0) - D(K>=0.01) as one statistic, jackknifed by
deleting each block of either bin from both at once. Outgroup filter (row outgroup_filter): the share of the
candidates without the outgroup rule (D discovers, D present) that "every outgroup E" removes, split into an
outgroup present and an outgroup missing.

Testability: a delete-one block jackknife sees only the variation
among blocks, and a variance it cannot estimate is not zero. Every jackknifed ratio here (D_pav and its
strata, dD, f) gets se and z nan, with its estimate kept and a findings line "not testable by the block
jackknife: <reason>", when a class it compares is empty (D_pav: n_RD or n_SD 0, where D sits at +-1 in every
pseudovalue; dD: either bin's; f: any of its four classes), when its informative loci lie in fewer than
JK_MIN_BLOCKS distinct blocks, or when its pseudovalues do not vary. A nan z is never read as gene flow.

Roles (`trios`): with --trio, as given (source "--trio"). Otherwise, for each three ingroup genomes,
`relations.sister_pair` gives the sister pair (x, y) and the donor; the recipient R is x when D_pav(x, y;
donor) >= 0, else y (source "data"). A nan D_pav (no informative locus, or no discovering genome) or a tie
at 0 keeps x, the first in sheet order (source "sheet order"), with a WARNING.

f4-ratio: for each second donor D2 (`second_donors`: D and D2 are each
other's closest genomes among R, S, D and D2 by PAV distance), on the D_pav rows with D2
resolved, from the same calls: f = sum(num) / sum(den), num = [R present, S E] - [S present, R E], den =
[D2 present, S E] - [S present, D2 E], jackknifed on D's blocks. A D2 that discovers the counted loci (a
fallback genome, D unannotated) is skipped with a findings line: it is present at every one by
construction, so den would be another statistic.

Block-size scan: the headline D_pav (stratum all) is jackknifed again with blocks of each BLOCK_SCAN_BP length
(1, 2, 5, 10 and 20 Mb), rows layer "trio", test "block_scan", stratum "<size> Mb", whose estimate is the
headline's and whose se, z and block counts are each length's own. When the verdict at |Z| >= Z_SIG differs from
the headline's at a length whose informative (R+D, S+D) loci lie in at least BLOCK_SCAN_MIN_BLOCKS blocks (the
jackknife se moves only through the blocks holding them, as stats.report.why counts testability; fewer make the se
itself unstable: a heuristic), one WARNING in "Trio tests" names the block-size scan and each such length.

Every other row is layer "pav" and declares its counted loci in `loci` (relations.py's grammar). Findings go to
the sections "Trio tests" and "Dating".
"""
from __future__ import annotations

import itertools
import math
from dataclasses import dataclass

import numpy as np

from synltr.analyses.geneflow.catalog import (
    AGE_LABELS,
    TIER_CODE,
    AnalysisError,
    Catalog,
    age_bin,
    check_trio,
    cluster_keys,
    patterns,
    site_blocks,
)
from synltr.analyses.geneflow.relations import Counted, discover, eligible_fallback, fallback_reasons, sister_pair
from synltr.base.constants import BLOCK_SCAN_BP, BLOCK_SCAN_MIN_BLOCKS, Z_SIG
from synltr.stats import report
from synltr.stats.report import Findings, StatRow, resolution
from synltr.stats.resample import Jk, d_ratio, jackknife

TRIO, DATING = "Trio tests", "Dating"  # findings sections (stats.report.SECTIONS)
GRADIENT = f"{AGE_LABELS[0]} minus {AGE_LABELS[2]}"  # the stratum of row D_pav_gradient
NOT_TESTABLE = "not testable by the block jackknife"  # findings and notes of a statistic with se and z nan
SCAN = "block_scan"  # the test of the block-size scan's rows
_YOUNG, _OLD = 0, 2  # the age_bin codes of the gradient's two bins
_HIGH = TIER_CODE["high"]


@dataclass(frozen=True)
class Trio:
    """A trio's roles: recipient R, its sister S and the donor D; `fixed` when given by --trio."""

    R: str
    S: str
    D: str
    fixed: bool

    @property
    def scope(self) -> str:
        """The trio's scope in stats.tsv: "R,S;D"."""
        return f"{self.R},{self.S};{self.D}"


# ---------------------------------------------------------------- D_pav


def d_pav(cat: Catalog, R: str, S: str, D: str, c: Counted, size: int, keep: np.ndarray | None = None
          ) -> tuple[Jk, int, int]:
    """D_pav(R, S; D) on the rows `c.rows & keep` (keep: a boolean mask over the n loci; None keeps all),
    every state from each row's discovering genome's own calls (c.gstar): (jackknife, n_RD, n_SD). n_RD
    counts R present and S E, n_SD the mirror. Blocks of `size` bp on D's site in each row's discoverer call.
    se and z are nan when the jackknife cannot estimate the variance (module docstring, Testability)."""
    jk, n_rd, n_sd, _, _ = _d_pav(cat, R, S, D, c, size, keep)
    return jk, n_rd, n_sd


def _d_pav(cat: Catalog, R: str, S: str, D: str, c: Counted, size: int, keep: np.ndarray | None
           ) -> tuple[Jk, int, int, str | None, int]:
    """`d_pav`, why it is not testable by the block jackknife (None when it is), and the number of blocks holding
    its informative loci."""
    rows = c.rows if keep is None else c.rows & _mask(cat, keep)
    a, b = patterns(cat, c.gstar, R, S)
    a, b = a[rows], b[rows]
    blocks = site_blocks(cat, c.gstar, D, rows, size)
    n_rd, n_sd = int(np.count_nonzero(a)), int(np.count_nonzero(b))
    informative = blocks[a | b]
    why = report.why({f"{R}+{D}": n_rd, f"{S}+{D}": n_sd}, informative)
    jk, why = report.guard(jackknife(d_ratio(a, b), blocks), why)
    return jk, n_rd, n_sd, why, len(np.unique(informative))


def _mask(cat: Catalog, keep: np.ndarray) -> np.ndarray:
    keep = np.asarray(keep)
    if keep.dtype != bool or keep.shape != cat.lid.shape:
        raise ValueError(f"keep must be a boolean mask over the {len(cat.lid)} loci, not {keep.dtype} of "
                         f"shape {keep.shape} (an index array?)")
    return keep


# ---------------------------------------------------------------- roles


def trios(cat: Catalog, fixed: tuple[str, str, str] | None, size: int, dist: np.ndarray
          ) -> tuple[list[Trio], list[StatRow], Findings]:
    """The trios to test and their `roles` rows (statistic source: "data", "--trio" or "sheet order").

    fixed (--trio): three distinct ingroup genomes R, S, D, else an AnalysisError; taken as given.
    Otherwise each three ingroup genomes in sheet order: `sister_pair` gives the pair (x, y) and the donor;
    R is x when D_pav(x, y; donor) >= 0, else y; a nan D_pav or a tie at 0 keeps x (source "sheet order")
    with a WARNING. Fewer than three ingroup genomes: no trio, and one findings line."""
    f = Findings()
    if fixed is not None:
        fixed = tuple(fixed)
        if len(fixed) != 3:
            raise AnalysisError(f"--trio: three genomes R,S,D are needed, not {', '.join(fixed) or 'none'}")
        check_trio(cat, fixed, "--trio")
        t = Trio(*fixed, fixed=True)
        f.add(TRIO, f"Trio {t.scope}: roles fixed by --trio (recipient {t.R}, sister {t.S}, donor {t.D}); "
                    "nothing inferred")
        return [t], [StatRow("pav", "roles", t.scope, "all", "source", "--trio", loci="none",
                             note="roles fixed by --trio, as R,S;D")], f
    if len(cat.ingroup) < 3:
        f.add(TRIO, "Trio tests skipped: fewer than three ingroup genomes")
        return [], [], f
    found, rows = [], []
    for genomes in itertools.combinations(cat.ingroup, 3):
        sc = sister_pair(cat, genomes, dist)
        (x, y), d = sc.pair, sc.donor
        c = discover(cat, x, y, d, dist, "E")
        if c is None:
            est, n, loci = math.nan, 0, "none"
            why = f"{d} has no annotation and no eligible fallback discovery genome"
        else:
            jk, n_xd, n_yd = d_pav(cat, x, y, d, c, size)
            est, n, loci = jk.est, n_xd + n_yd, c.loci
            why = (f"a tie, {n_xd} {x}+{d} and {n_yd} {y}+{d} loci" if est == 0 else
                   f"none of the {np.count_nonzero(c.rows)} counted loci is {x}+{d} or {y}+{d}")
        t = Trio(y, x, d, False) if est < 0 else Trio(x, y, d, False)  # nan < 0 is False: x, sheet order
        found.append(t)
        sign = f"D_pav({x},{y};{d}) = {report.fmt(est, '+.3f')}"
        pair = f"{x}+{y}" + ("; ambiguous, see Relationships" if sc.ambiguous else "")
        if math.isnan(est) or est == 0:  # nothing in the data tells the recipient: sheet order decides
            source = "sheet order"
            note = (f"sister pair {pair} (rows sharing: merged states, any discovery); recipient not "
                    f"identified ({sign}: {why}), {x}, the first in sheet order, taken as R")
            f.warn(TRIO, f"trio {','.join(genomes)}: the recipient is not identified ({sign}: {why}); {x}, "
                         "the first of the sister pair in sheet order, is taken as R")
        else:
            source = "data"
            note = (f"sister pair {pair} (rows sharing: merged states, any discovery); recipient {t.R} by "
                    f"the sign of {sign} on these loci")
            f.add(TRIO, f"Trio {','.join(genomes)}: roles from the data: sisters {x} and {y} (the pair "
                        f"sharing the most insertions{', an ambiguous pair' if sc.ambiguous else ''}), donor "
                        f"{d}; recipient {t.R}, the sister sharing more insertions with {d} ({sign})")
        rows.append(StatRow("pav", "roles", t.scope, "all", "source", source, n=n, loci=loci, note=note))
    return found, rows, f


# ---------------------------------------------------------------- the tests of one trio


@dataclass
class _Stratum:
    name: str
    jk: Jk
    n_rd: int
    n_sd: int
    n: int              # the stratum's loci: counted loci (nested_collapsed: the informative loci kept)
    loci: str
    note: str
    why: str | None     # why it is not testable by the block jackknife; None when it is


def trio_tests(cat: Catalog, trio: Trio, size: int, dist: np.ndarray) -> tuple[list[StatRow], Findings]:
    """The tests of one trio (module docstring): rows resolution, outgroup_filter, D_pav per stratum, the
    block-size scan of the primary D_pav, D_pav_gradient and f4_ratio per second donor, with their findings.
    Skipped with a WARNING when no genome may discover (D unannotated, no eligible fallback genome)."""
    f = Findings()
    R, S, D, scope = trio.R, trio.S, trio.D, trio.scope
    c = discover(cat, R, S, D, dist, "E")
    if c is None:
        why = fallback_reasons(eligible_fallback(cat, R, S, D, dist)[1])
        f.warn(TRIO, f"D_pav({scope}) skipped, with its strata, age gradient and f4-ratio: {D} has no "
                     f"annotation and no annotated ingroup genome outside the trio splits with it from {R} "
                     f"and {S} by the four-point condition ({why})")
        return [], f
    loose = discover(cat, R, S, D, dist, "none")
    disc = ",".join(c.disc)
    strata, ages_line = _strata(cat, trio, c, loose, size)
    rows: list[StatRow] = []
    if strata and strata[0].name == "all":
        f.add(TRIO, _primary_line(trio, c, strata[0], size))
    else:
        f.add(TRIO, f"D_pav({scope}): no counted locus ({np.count_nonzero(c.cand)} loci {disc} discovers "
                    f"with {D} present and every outgroup E; none has {R} and {S} both resolved)")
    n_cand = int(np.count_nonzero(c.cand))
    if n_cand:
        rows += resolution("resolution", scope, _resolved(cat, c, (R, S)), n_cand,
                           f"disc={disc};present={D};outgroup=E", f, TRIO)
    rows.append(_outgroup_filter(cat, trio, c, loose, f))
    for s in strata:
        rows += _stratum_rows(scope, s, D, size)
    if strata and strata[0].name == "all":
        rows += _scan_rows(cat, trio, c, strata[0], size, f)
    checks = [s for s in strata if s.name != "all"]
    if checks:
        f.add(TRIO, f"D_pav({scope}) check strata, D (Z; n_RD/n_SD): "
                    + "; ".join(f"{s.name} {report.fmt(s.jk.est, '+.3f')} ({_z_text(s)}; {s.n_rd}/{s.n_sd})"
                                for s in checks)
                    + ". high_tier: the discoverer's calls into the sisters at tier high; nested_collapsed: "
                      "one informative locus per nested cluster; no_outgroup_filter: without the outgroup "
                      "rule")
    if ages_line:
        f.add(TRIO, ages_line)
    if c.rows.any():
        d_all = strata[0].jk.est if strata and strata[0].name == "all" else math.nan
        rows += _gradient(cat, trio, c, size, f, d_all)
    donors = second_donors(cat, trio, dist)
    if not donors:
        f.add(TRIO, f"f4-ratio skipped: no second donor (an ingroup genome outside the trio that is closer "
                    f"to {D}, and {D} to it, than either is to {R} or {S}, by PAV distance)")
    for d2 in donors:
        rows += _f4(cat, trio, d2, c, size, f)
    return rows, f


def _strata(cat: Catalog, trio: Trio, c: Counted, loose: Counted, size: int
            ) -> tuple[list[_Stratum], str | None]:
    """D_pav per stratum (module docstring), in row order, a stratum with no counted locus left out; and the
    findings line on the age strata when some or all counted loci have no K2P (else None)."""
    R, S, D, scope = trio.R, trio.S, trio.D, trio.scope
    disc = ",".join(c.disc)
    a, b = patterns(cat, c.gstar, R, S)
    ages = _ages(cat, c.gstar)
    n_rows = int(np.count_nonzero(c.rows))
    n_aged = int(np.count_nonzero(c.rows & (ages >= 0)))
    plan: list[tuple[str, np.ndarray | None, Counted, str]] = [("all", None, c, "the primary count")]
    line = None
    if n_rows and not n_aged:
        line = (f"D_pav({scope}) age strata skipped: no K2P for {disc}'s copies at the {n_rows} counted loci "
                "(an annotation without LTR identity)")
    else:
        if n_aged < n_rows:
            line = (f"D_pav({scope}) age strata: {n_rows - n_aged} of the {n_rows} counted loci have no K2P "
                    f"for {disc}'s copy and enter no age bin")
        plan += [(label, ages == code, c, f"K2P of {disc}'s copy: {label}")
                 for code, label in enumerate(AGE_LABELS)]
    high = (cat.ptier_by(c.gstar, R) == _HIGH) & (cat.ptier_by(c.gstar, S) == _HIGH)
    plan += [("high_tier", high, c, f"{disc}'s calls into {R} and {S} at tier high"),
             ("nested_collapsed", _collapsed(cat, c, a | b), c,
              (f"the informative loci, one per nested cluster of {disc}'s elements (the first in position "
               "order)")),
             ("no_outgroup_filter", None, loose, "without the outgroup rule")]
    out = []
    for name, keep, cc, note in plan:
        n = int(np.count_nonzero(cc.rows if keep is None else cc.rows & keep))
        if n:
            jk, n_rd, n_sd, why, _ = _d_pav(cat, R, S, D, cc, size, keep)
            out.append(_Stratum(name, jk, n_rd, n_sd, n, cc.loci, note, why))
    return out, line


def _ages(cat: Catalog, gstar: np.ndarray) -> np.ndarray:
    """(n,) int8 age bin (`age_bin`) of each row's discovering copy; -1 without a K or a discoverer."""
    k = np.full(len(cat.lid), np.nan)
    i = np.flatnonzero(gstar >= 0)
    k[i] = cat.k[i, gstar[i]]
    return age_bin(k)


def _collapsed(cat: Catalog, c: Counted, informative: np.ndarray) -> np.ndarray:
    """(n,) mask: of the informative counted rows, the first of each nested cluster (`cluster_keys`: the
    discovering genome and its outermost annotated element) in that genome's (seq, start, -end) order."""
    sel = c.rows & informative
    idx = np.flatnonzero(sel)
    keep = np.zeros(len(cat.lid), dtype=bool)
    if not idx.size:
        return keep
    g = c.gstar[idx].astype(np.intp)
    order = np.lexsort((-cat.end[idx, g], cat.start[idx, g], cat.seq[idx, g]))  # stable: ties keep row order
    keys = cluster_keys(cat, c.gstar, sel)  # aligned with idx (catalog order)
    seen = set()
    for j in order.tolist():
        if keys[j] not in seen:
            seen.add(keys[j])
            keep[idx[j]] = True
    return keep


def _stratum_rows(scope: str, s: _Stratum, D: str, size: int) -> list[StatRow]:
    blocks = f"{s.jk.nblocks} jackknife blocks of {size:,} bp on {D}'s sites"
    if s.why:
        blocks += f"; {NOT_TESTABLE}: {s.why}"
    n_is = ("the informative loci kept, one per nested cluster" if s.name == "nested_collapsed"
            else "the counted loci")
    return [StatRow("pav", "D_pav", scope, s.name, "D", s.jk.est, s.jk.se, s.jk.z, n=s.n_rd + s.n_sd,
                    loci=s.loci, note=f"{s.note}; {blocks}"),
            StatRow("pav", "D_pav", scope, s.name, "n_RD", s.n_rd, n=s.n, loci=s.loci,
                    note=f"R present, S E (n: {n_is})"),
            StatRow("pav", "D_pav", scope, s.name, "n_SD", s.n_sd, n=s.n, loci=s.loci,
                    note=f"S present, R E (n: {n_is})")]


def block_scan(cat: Catalog, trio: Trio, c: Counted) -> list[tuple[int, Jk, str | None, int]]:
    """The primary D_pav of `c`'s counted loci jackknifed with blocks of each BLOCK_SCAN_BP length: (length,
    jackknife, why it is not testable by the block jackknife or None, the blocks holding informative loci) per
    length."""
    out = []
    for bp in BLOCK_SCAN_BP:
        jk, _, _, why, informative = _d_pav(cat, trio.R, trio.S, trio.D, c, bp, None)
        out.append((bp, jk, why, informative))
    return out


def significant(jk: Jk) -> bool:
    """The verdict at |Z| >= Z_SIG; a nan Z (not testable) is never significant."""
    return bool(abs(jk.z) >= Z_SIG)


def _scan_rows(cat: Catalog, trio: Trio, c: Counted, head: _Stratum, size: int, f: Findings) -> list[StatRow]:
    """The block-size scan's rows of the primary D_pav `head` (jackknifed with the rule's `size`-bp blocks), and
    the WARNING when its verdict differs at a length whose informative loci lie in at least BLOCK_SCAN_MIN_BLOCKS
    blocks (module docstring)."""
    rows, flips = [], []
    sig = significant(head.jk)
    for bp, jk, why, informative in block_scan(cat, trio, c):
        mb = f"{bp // 1_000_000} Mb"
        note = (f"{jk.nblocks} jackknife blocks of {bp:,} bp on {trio.D}'s sites, {informative} holding informative "
                "loci" + (f"; {NOT_TESTABLE}: {why}" if why else ""))
        rows.append(StatRow("trio", SCAN, trio.scope, mb, "D", jk.est, jk.se, jk.z, n=head.n_rd + head.n_sd,
                            loci=head.loci, note=note))
        if informative >= BLOCK_SCAN_MIN_BLOCKS and significant(jk) != sig:
            flips.append(f"{mb} (Z {report.fmt(jk.z, '.2f')}, informative loci in {informative} blocks)")
    if flips:
        verdict = "significant" if sig else "not significant"
        f.warn(TRIO, f"D_pav({trio.scope}) block-size scan: {verdict} at |Z| >= {Z_SIG:g} with the rule's "
                     f"{size:,}-bp blocks (Z {report.fmt(head.jk.z, '.2f')}), {'not ' if sig else ''}significant "
                     f"with blocks of {'; '.join(flips)}: the jackknife se depends on the block length here, so "
                     f"read the verdict with care (lengths whose informative loci lie in fewer than "
                     f"{BLOCK_SCAN_MIN_BLOCKS} blocks are not compared, a heuristic)")
    return rows


def _z_text(s: _Stratum) -> str:
    """A stratum's Z for a findings line, or why it has none."""
    return f"{NOT_TESTABLE}: {s.why}" if s.why else f"Z {report.fmt(s.jk.z, '.2f')}"


def _primary_line(trio: Trio, c: Counted, s: _Stratum, size: int) -> str:
    R, S, D = trio.R, trio.S, trio.D
    disc = ",".join(c.disc)
    if s.why:
        value = f"{report.fmt(s.jk.est, '+.4f')} ({NOT_TESTABLE}: {s.why})"
    else:
        sig = (f"significant, |Z| >= {Z_SIG:g}" if abs(s.jk.z) >= Z_SIG
               else f"not significant, |Z| < {Z_SIG:g}")
        value = f"{report.fmt(s.jk.est, '+.4f')} ± {report.fmt(s.jk.se, '.4f')} (Z {report.fmt(s.jk.z, '.2f')}; {sig})"
    fallback = "" if c.disc == (D,) else (
        f" ({D} has no annotation; {disc} split with it from {R} and {S} by the four-point condition"
        + (f"; not eligible: {fallback_reasons(c.excluded)}" if c.excluded else "") + ")")
    chosen = "" if trio.fixed or not s.jk.est > 0 else (
        f"; the roles come from the data, {R} being the sister with the larger count, so D_pav >= 0 by "
        "construction and |Z| is the test")
    return (f"D_pav({trio.scope}) = {value}: n_RD {s.n_rd}, n_SD {s.n_sd}, on {s.n} loci {disc} "
            f"discovers{fallback} with {D} present, {R} and {S} resolved and every outgroup E, each state "
            f"from {disc}'s own calls at any tier; jackknife over {s.jk.nblocks} blocks of {size:,} bp on "
            f"{D}'s sites. Positive: {R} shares more insertions with {D} than {S} does; lineage sorting "
            f"alone gives 0 in expectation{chosen}")


def _resolved(cat: Catalog, c: Counted, genomes: tuple[str, ...]) -> dict[str, float]:
    """The fraction of the candidates (c.cand) where each genome is resolved in gstar's own call."""
    n = int(np.count_nonzero(c.cand))
    return {g: int(np.count_nonzero(c.cand & (cat.pcls_by(c.gstar, g) >= 0))) / n for g in genomes}


def _outgroup_filter(cat: Catalog, trio: Trio, c: Counted, loose: Counted, f: Findings) -> StatRow:
    """Row outgroup_filter (removed_frac) and its findings line: the share of the candidates without the
    outgroup rule that "every outgroup E" removes, split into an outgroup present and an outgroup missing."""
    D, scope = trio.D, trio.scope
    disc = ",".join(c.disc)
    n = int(np.count_nonzero(loose.cand))
    dropped = loose.cand & ~c.cand
    present = np.zeros(len(cat.lid), dtype=bool)
    for o in cat.outgroups:
        present |= cat.pcls_by(c.gstar, o) == 1
    n_present = int(np.count_nonzero(dropped & present))
    n_missing = int(np.count_nonzero(dropped)) - n_present
    frac = (n_present + n_missing) / n if n else math.nan
    split = (f"{n_present} with an outgroup present (older than the ingroup split) and {n_missing} with an "
             "outgroup missing (D/O/G/U or no call)")
    if cat.outgroups:
        f.add(TRIO, f"Outgroup filter of D_pav({scope}): {frac:.2f} of the {n} candidate loci ({disc} "
                    f"discovers, {D} present) removed by the rule that every outgroup "
                    f"({', '.join(cat.outgroups)}) be E in {disc}'s own call: {split}; the "
                    "no_outgroup_filter stratum keeps them")
    else:
        f.add(TRIO, f"Outgroup filter of D_pav({scope}): {frac:.2f} of the {n} candidate loci removed: the "
                    "samples sheet names no outgroup, so old insertions are not set aside")
    return StatRow("pav", "outgroup_filter", scope, "all", "removed_frac", frac, n=n,
                   loci=f"disc={disc};present={D};outgroup=none",
                   note=f"the share of these candidates the rule every outgroup E removes: {split}")


def _gradient(cat: Catalog, trio: Trio, c: Counted, size: int, f: Findings, d_all: float) -> list[StatRow]:
    """Row D_pav_gradient: dD = D(K=0) - D(K>=0.01), one statistic, jackknifed by deleting each block of
    either bin from both bins at once; findings in "Dating". Skipped with one line without a K or an
    informative locus in either bin (dD undefined). se and z are nan, dD kept, when either bin is not
    testable by the block jackknife (module docstring, Testability): its variance would drop out of the se.

    The verdict reads dD in the direction of the overall excess: d_all, D_pav over all
    counted loci, names the sister with the excess (R when > 0, S when < 0), and sign(d_all) x dD > 0 says that
    sister's excess is concentrated in young insertions. With --trio roles S may carry the excess (D_pav < 0):
    a young-concentrated excess of S's then gives dD < 0. No direction is read when d_all is 0 or nan."""
    R, S, D, scope = trio.R, trio.S, trio.D, trio.scope
    disc = ",".join(c.disc)
    ages = _ages(cat, c.gstar)
    young, old = c.rows & (ages == _YOUNG), c.rows & (ages == _OLD)
    union = young | old
    a, b = patterns(cat, c.gstar, R, S)
    bins = [(label, m, int(np.count_nonzero(a & m)), int(np.count_nonzero(b & m)))
            for label, m in ((AGE_LABELS[_YOUNG], young), (AGE_LABELS[_OLD], old))]
    if not np.count_nonzero(c.rows & (ages >= 0)):
        f.add(DATING, f"Age gradient of D_pav({scope}) skipped: no K2P for {disc}'s copies")
        return []
    empty = [label for label, _, n_rd, n_sd in bins if not n_rd + n_sd]
    if empty:
        f.add(DATING, f"Age gradient of D_pav({scope}) skipped: no informative locus at "
                      f"{' or '.join(empty)}")
        return []
    blocks = site_blocks(cat, c.gstar, D, union, size)
    d_young = d_ratio((a & young)[union], (b & young)[union])
    d_old = d_ratio((a & old)[union], (b & old)[union])
    reasons = []
    for label, m, n_rd, n_sd in bins:  # the rule applies to either bin
        why = report.why({f"{R}+{D}": n_rd, f"{S}+{D}": n_sd}, blocks[((a | b) & m)[union]])
        if why:
            reasons.append(f"{label}: {why} (n_RD {n_rd}, n_SD {n_sd})")
    jk = jackknife(lambda keep: d_young(keep) - d_old(keep), blocks)
    jk, why = report.guard(jk, "; ".join(reasons) or None)
    every = np.ones(np.count_nonzero(union), dtype=bool)
    d0, d2 = d_young(every), d_old(every)
    n_inf = sum(n_rd + n_sd for _, _, n_rd, n_sd in bins)
    parts = [f"D({label}) {report.fmt(d, '+.3f')} (n_RD {n_rd}, n_SD {n_sd})"
             for (label, _, n_rd, n_sd), d in zip(bins, (d0, d2), strict=True)]
    head = (f"Age gradient of D_pav({scope}), by the LTR-LTR K2P of {disc}'s copy: {parts[0]} minus "
            f"{parts[1]} = dD {report.fmt(jk.est, '+.3f')}")
    null = "Lineage sorting alone gives D = 0 at every age"
    sign = 1 if d_all > 0 else -1 if d_all < 0 else 0   # nan: 0
    owner, other = (R, S) if sign > 0 else (S, R)
    excess = (f"the excess, {owner}'s ({owner} shares more insertions with {D} than {other} does: D_pav "
              f"{report.fmt(d_all, '+.3f')})")
    if why:  # never read as gene flow: nothing here measures dD's sampling variance
        f.add(DATING, f"{head}, {NOT_TESTABLE}: {why}. {null}")
    else:
        zz = sign * jk.z    # in the direction of the excess
        if not sign:
            sig = (f"significant, |Z| >= {Z_SIG:g}" if abs(jk.z) >= Z_SIG
                   else f"not significant, |Z| < {Z_SIG:g}")
            verdict = (f"{sig}, but D_pav over all counted loci is {report.fmt(d_all, '+.3f')}: no excess, so no age "
                       "profile of one to read")
        elif zz >= Z_SIG:
            verdict = (f"significant at |Z| >= {Z_SIG:g}, read in the direction of {excess}: that excess is "
                       f"concentrated in young insertions, consistent with gene flow from {D}-like donors into "
                       f"{owner} rather than lineage sorting alone")
        elif zz <= -Z_SIG:
            verdict = (f"significant at |Z| >= {Z_SIG:g}, read in the direction of {excess}: that excess is "
                       "concentrated in older insertions, not the age signature of gene flow")
        else:
            verdict = f"not significant, |Z| < {Z_SIG:g}; {excess} shows no clear age profile"
        f.add(DATING, f"{head} ± {report.fmt(jk.se, '.3f')} (Z {report.fmt(jk.z, '.2f')}; one jackknife deleting each "
                      f"of {jk.nblocks} blocks from both bins): {verdict}. {null}")
    note = (f"D({AGE_LABELS[_YOUNG]}) {report.fmt(d0, '+.4f')} minus D({AGE_LABELS[_OLD]}) {report.fmt(d2, '+.4f')}; n: "
            f"the informative loci of both bins; one delete-one block jackknife over the {jk.nblocks} blocks "
            f"of both bins ({size:,} bp, {D}'s sites)" + (f"; {NOT_TESTABLE}: {why}" if why else ""))
    return [StatRow("pav", "D_pav_gradient", scope, GRADIENT, "dD", jk.est, jk.se, jk.z, n=n_inf, loci=c.loci,
                    note=note)]


# ---------------------------------------------------------------- f4-ratio


def second_donors(cat: Catalog, trio: Trio, dist: np.ndarray) -> list[str]:
    """The ingroup genomes g outside the trio that qualify as a second donor (on PAV
    distances): dist[D, g] < min(dist[x, y] for x in (D, g) for y in (R, S)), so D and g are each other's
    closest genomes among R, S, D and g. Sheet order; a nan distance never qualifies."""
    gd = cat.g(trio.D)
    out = []
    for g in cat.ingroup:
        if g in (trio.R, trio.S, trio.D):
            continue
        gg = cat.g(g)
        others = [float(dist[x, cat.g(y)]) for x in (gd, gg) for y in (trio.R, trio.S)]
        near = float(dist[gd, gg])
        # every distance must be finite: min() with a nan is order-dependent
        if all(map(math.isfinite, (near, *others))) and near < min(others):
            out.append(g)
    return out


def _f4(cat: Catalog, trio: Trio, d2: str, c: Counted, size: int, f: Findings) -> list[StatRow]:
    """Rows f4_ratio (scope R,S;D;D2) for the second donor d2 on the D_pav count `c`: f (jackknifed), num, den,
    and the resolution of S and d2 among the D_pav candidates (module docstring); [] when d2 discovers the
    counted loci. Its findings: the f4 line (or why d2 is skipped), the resolution WARNING and a den <= 0
    WARNING."""
    R, S, D = trio.R, trio.S, trio.D
    if d2 not in cat.ingroup or d2 in (R, S, D):
        raise AnalysisError(f"f4-ratio: the second donor {d2!r} must be an ingroup genome outside "
                            f"{trio.scope}")
    if d2 in c.disc:  # its own column is P at each of its members: [S present, d2 E] could never be counted
        f.add(TRIO, f"f4-ratio of {trio.scope} with second donor {d2} skipped: {d2} discovers the counted "
                    f"loci ({D} has no annotation), so {d2} is present at every one by construction")
        return []
    scope, disc = f"{trio.scope};{d2}", ",".join(c.disc)
    pr, ps, p2 = (cat.pcls_by(c.gstar, x) for x in (R, S, d2))
    rows = c.rows & (p2 >= 0)
    classes = {f"[{R} present, {S} E]": (pr == 1) & (ps == 0),
               f"[{S} present, {R} E]": (ps == 1) & (pr == 0),
               f"[{d2} present, {S} E]": (p2 == 1) & (ps == 0),
               f"[{S} present, {d2} E]": (ps == 1) & (p2 == 0)}
    rd, sr, xd, sx = (m[rows] for m in classes.values())
    num = rd.astype(np.int64) - sr
    den = xd.astype(np.int64) - sx

    def ratio(keep: np.ndarray) -> float:
        d = den[keep].sum()
        return float(num[keep].sum() / d) if d else math.nan
    blocks = site_blocks(cat, c.gstar, D, rows, size)
    counts = {name: int(np.count_nonzero(m[rows])) for name, m in classes.items()}
    jk, why = report.guard(jackknife(ratio, blocks), report.why(counts, blocks[rd | sr | xd | sx]))
    n, s_num, s_den = int(np.count_nonzero(rows)), int(num.sum()), int(den.sum())
    loci = f"disc={disc};present={D};res={R},{S},{d2};outgroup=E"
    out = [StatRow("pav", "f4_ratio", scope, "all", "f", jk.est, jk.se, jk.z, n=n, loci=loci,
                   note=f"f = num / den, the admixture fraction of {R} from {D}'s lineage; {jk.nblocks} "
                        f"jackknife blocks of {size:,} bp on {D}'s sites"
                        + (f"; {NOT_TESTABLE}: {why}" if why else "")),
           StatRow("pav", "f4_ratio", scope, "all", "num", s_num, n=n, loci=loci,
                   note=f"[{R} present, {S} E] - [{S} present, {R} E]"),
           StatRow("pav", "f4_ratio", scope, "all", "den", s_den, n=n, loci=loci,
                   note=f"[{d2} present, {S} E] - [{S} present, {d2} E]")]
    value = f"f = {s_num} / {s_den} = {report.fmt(jk.est, '.3f')}" + (
        f" ({NOT_TESTABLE}: {why})" if why else
        f" ± {report.fmt(jk.se, '.3f')} (Z {report.fmt(jk.z, '.2f')}; {jk.nblocks} jackknife blocks)")
    f.add(TRIO, f"f4-ratio of {trio.scope} with second donor {d2} (the admixture fraction of {R} from {D}'s "
                f"lineage): {value}, on {n} loci {disc} discovers with {D} present, every outgroup E and "
                f"{R}, {S}, {d2} resolved ({disc}'s own calls); classes "
                + ", ".join(f"{k} {v}" for k, v in counts.items()))
    n_cand = int(np.count_nonzero(c.cand))
    if n_cand:  # den compares d2 with S: their resolution among the D_pav candidates (WARNING after the line)
        out += resolution("f4_ratio", scope, _resolved(cat, c, (S, d2)), n_cand,
                          f"disc={disc};present={D};outgroup=E", f, TRIO)
    if s_den <= 0:
        f.warn(TRIO, f"f4-ratio of {trio.scope} with second donor {d2}: den {s_den} <= 0, so {d2} shares no "
                     f"more insertions with {D} than {S} does on these loci and f is not an admixture "
                     "fraction")
    return out


