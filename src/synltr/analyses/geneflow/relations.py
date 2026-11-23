"""Relationships between the genomes of a run: the run summary, the one discovery function
(`discover`, with its four-point fallback), PAV distances and a neighbor-joining tree, shared derived
insertions, the sharing counts that pick each trio's sister pair, the lineage-sorting (KKSC) symmetry test and
outgroup conflicts.

Which states each part reads (catalog.py):
  merged  loci.tsv's states (`Catalog.cls`), for the any-discovery descriptions only: the run summary,
          distances and the tree, shared derived insertions, and the sharing counts that pick the sister
          pair. Their rows declare `disc=any`.
  own     each row's discovering genome's own calls (`Catalog.pcls_by`), for every statistic with a discovery
          rule. `discover` is the one discovery function: D_pav, the f4-ratio and the identity layer call it
          with the outgroup rule "E" (an outgroup must be E), the KKSC symmetry test with
          "not_present" (loci an outgroup carries are set aside). The outgroup-conflict bound
          reads each genome's state from other ingroup genomes' discoveries (`_conflicts`).

A missing state (D/O/G/U, or no call) is never absence and is never imputed: each count takes the loci where
the genomes it names are resolved (complete cases). Each StatRow names its counted loci in `loci`: `disc=` the
discovering genome(s), in the order each row takes its first member (`any`: merged states; `any_but_<g>`: any
genome but g), `present=` / `absent=` / `res=` the genomes required present (P/S/T), absent (E) or resolved
(P/S/T or E), and `outgroup=` the outgroup rule: `E` every outgroup E, `not_present` no outgroup P/S/T,
`present` some outgroup P/S/T. `none`: the row counts no loci (a property of a genome, a count of annotation
records, or a test skipped for want of a discovering genome).

Tree: neighbor joining (Saitou & Nei 1987, in the Q-matrix form of Studier & Keppler 1988) on every genome's
PAV distances, ties to the first pair in sheet order, negative branch lengths set to 0; rooted at the midpoint
of the branch that splits the outgroups from the ingroup (the root's place on that branch is not identifiable,
and the midpoint keeps every leaf-to-leaf path length): one outgroup's own branch, or several outgroups' clade
(when they do not form one, the first outgroup's branch, with a WARNING); unrooted without an outgroup.
"""
from __future__ import annotations

import itertools
import math
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np

from synltr.analyses.geneflow.catalog import AnalysisError, Catalog, check_trio, first_member
from synltr.base.constants import P_WARN
from synltr.stats.report import Findings, StatRow, fmt, resolution

RUN, SECTION = "Run", "Relationships"  # findings sections (stats.report.SECTIONS)
OUTGROUP_RULES = ("E", "not_present", "none")  # discover's outgroup rules (module docstring)
CONFLICT_NOTE = ("an upper bound on the combined rate of ingroup false E and outgroup false presence; deep "
                 "lineage sorting and parallel insertions also give it")
KKSC_NOTE = ("it ignores block dependence (nested clusters, linkage), so its p is anti-conservative; the "
             "jackknifed D_pav is the test of record")
_TIE = 1e-12  # NJ: Q values this close (relative) to the minimum tie: floating-point rounding, no threshold
_QC_HEADER = ["scope", "metric", "value"]


@dataclass(frozen=True)
class Node:
    """A tree node: a leaf (a genome label, no children) or an inner node (name None); `length` is the branch
    above it (0 at the top)."""

    name: str | None
    children: tuple[Node, ...] = ()
    length: float = 0.0

    def leaves(self) -> frozenset[str]:
        """The labels of the leaves at or below this node."""
        if not self.children:
            return frozenset() if self.name is None else frozenset((self.name,))
        return frozenset().union(*(c.leaves() for c in self.children))

    def newick(self) -> str:
        """Newick text: leaf labels, branch lengths to 6 significant digits (none above the top), ";" at the
        end. Labels need no quoting: the samples sheet allows only A-Z a-z 0-9 _ . -."""
        return self._newick() + ";"

    def _newick(self) -> str:
        inner = ("(" + ",".join(f"{c._newick()}:{c.length:.6g}" for c in self.children) + ")"
                 if self.children else "")
        return inner + (self.name or "")


@dataclass
class Counted:
    """The loci a statistic with a discovery rule counts (`discover`): masks over the n loci (catalog order)
    and each row's discovering genome, whose own calls every state is read from."""

    rows: np.ndarray          # (n,) bool: candidates with x and y resolved in gstar's calls
    cand: np.ndarray          # (n,) bool: gstar exists, z present, outgroup rule (gstar's calls)
    gstar: np.ndarray         # (n,) int16: each row's discovering genome (z, or a fallback genome); -1 none
    disc: tuple[str, ...]     # z, or the eligible fallback genomes
    loci: str                 # e.g. "disc=D;present=D;res=R,S;outgroup=E"
    excluded: dict[str, str]  # fallback genomes left out -> why


@dataclass(frozen=True)
class SisterCall:
    """A trio's sister pair, from the sharing counts (merged states), and the KKSC symmetry test, from the
    third genome's own calls (`sister_pair`)."""

    pair: tuple[str, str]     # sheet order
    donor: str                # the third genome
    counts: dict[str, int]    # "x+y" -> loci with x and y present and the third E (`trio_sharing`)
    p_top: float              # two-sided binomial p, the top count against the second; nan if both are 0
    p_kksc: float             # two-sided binomial p, n_x+z against n_y+z; nan without them or a discoverer
    kksc_loci: str            # the KKSC test's loci (`discover`'s), or "none"
    ambiguous: bool           # p_top >= P_WARN, or nan


# ---------------------------------------------------------------- run summary


def run_summary(cat: Catalog) -> tuple[list[StatRow], Findings]:
    """Rows (test "catalog") and findings (section "Run"): the number of loci; per genome whether it is
    annotated and its resolved fraction (merged state P/S/T or E, over all loci); and qc.tsv's
    annot_sensitivity and every skipped_<reason> count of each genome, quoted as written. A
    missing or unreadable qc.tsv is one findings line."""
    f = Findings()
    n = len(cat.lid)
    rows = [StatRow("pav", "catalog", "all", "all", "n_loci", n, n=n, loci="disc=any")]
    quoted, problem = _read_qc(cat.outdir / "qc.tsv", cat.labels)
    f.add(RUN, f"{n} loci; genomes {', '.join(cat.labels)} (sheet order); outgroup(s): "
               f"{', '.join(cat.outgroups) or 'none'}")
    for x in cat.labels:
        frac = np.count_nonzero(cat.cls(x) >= 0) / n if n else math.nan
        rows.append(StatRow("pav", "catalog", x, "all", "annotated", int(cat.annotated[x]), loci="none"))
        rows.append(StatRow("pav", "catalog", x, "all", "resolved_frac", frac, n=n, loci="disc=any"))
        for metric, text in quoted.get(x, ()):
            sensitivity = metric == "annot_sensitivity"
            rows.append(StatRow(
                "pav", "catalog", x, "all", metric, _qc_value(text),
                loci=f"disc=any_but_{x};present={x}" if sensitivity else "none",
                note=(f"qc.tsv: among loci with a high-tier P call into {x} from another genome's member, "
                      f"the fraction {x}'s annotation holds" if sensitivity
                      else "qc.tsv: annotation records the loader skipped")))
        line = (f"{x}: {'annotated' if cat.annotated[x] else 'genotype-only'}"
                f"{', outgroup' if x in cat.outgroups else ''}; resolved (P/S/T or E) at {fmt(frac, '.3g')} of the "
                "loci (merged states)")
        if quoted.get(x):
            line += "; qc.tsv: " + ", ".join(f"{m} {t}" for m, t in quoted[x])
        f.add(RUN, line)
    if problem:
        f.add(RUN, problem)
    return rows, f


def _read_qc(path: Path, labels: tuple[str, ...]) -> tuple[dict[str, list[tuple[str, str]]], str | None]:
    """qc.tsv's annot_sensitivity and skipped_<reason> rows per genome, (metric, value as written) in file
    order; and a findings line when qc.tsv is missing, unreadable or not a qc.tsv (else None)."""
    if not path.is_file():
        return {}, (f"qc.tsv not found in {path.parent}: the annotation sensitivity and the annotation "
                    "loaders' skip counts are not quoted")
    out: dict[str, list[tuple[str, str]]] = {}
    try:
        with open(path, encoding="utf-8", newline="") as fh:
            if fh.readline().rstrip("\r\n").split("\t") != _QC_HEADER:
                return {}, (f"{path}: not a qc.tsv of synltr genotype (header not scope, metric, value); "
                            "nothing quoted")
            for line in fh:
                fields = line.rstrip("\r\n").split("\t")
                if (len(fields) == 3 and fields[0] in labels
                        and (fields[1] == "annot_sensitivity" or fields[1].startswith("skipped_"))):
                    out.setdefault(fields[0], []).append((fields[1], fields[2]))
    except (OSError, UnicodeDecodeError) as exc:
        return {}, f"cannot read {path} ({exc}); nothing quoted"
    return out, None


def _qc_value(text: str) -> float | int | str:
    """A qc.tsv value: "." is nan; an integer, a float, else the text."""
    if text == ".":
        return math.nan
    for kind in (int, float):
        try:
            return kind(text)
        except ValueError:
            pass
    return text


# ---------------------------------------------------------------- distances, shared derived insertions, tree


def pav_distances(cat: Catalog) -> tuple[np.ndarray, np.ndarray]:
    """(G, G) PAV distances and (G, G) loci counted, merged states (any discovery, every locus, outgroups
    included): d(x, y) = the loci where x and y are both resolved and differ (present against E), over the
    loci where both are resolved; nan where none is. The diagonal is 0, its n each genome's resolved loci."""
    G = len(cat.labels)
    cls = [cat.cls(x) for x in cat.labels]
    res = [c >= 0 for c in cls]
    d = np.zeros((G, G))
    n = np.zeros((G, G), dtype=np.int64)
    for a in range(G):
        n[a, a] = np.count_nonzero(res[a])
        for b in range(a + 1, G):
            both = res[a] & res[b]
            m = int(np.count_nonzero(both))
            n[a, b] = n[b, a] = m
            d[a, b] = d[b, a] = np.count_nonzero(both & (cls[a] != cls[b])) / m if m else math.nan
    return d, n


def shared_derived(cat: Catalog) -> np.ndarray:
    """(G, G) int64 counts of the derived insertions two ingroup genomes share,
    merged states: present in both, no outgroup present, and E in at least one other ingroup genome (so not
    fixed in the ingroup; with two ingroup genomes there is none). Outgroup rows and columns and the diagonal
    are 0."""
    G = len(cat.labels)
    out = np.zeros((G, G), dtype=np.int64)
    cls = {x: cat.cls(x) for x in cat.ingroup}
    free = _no_outgroup_present(cat)
    for x, y in itertools.combinations(cat.ingroup, 2):
        absent = np.zeros(len(cat.lid), dtype=bool)
        for w in cat.ingroup:
            if w not in (x, y):
                absent |= cls[w] == 0
        both = free & absent & (cls[x] == 1) & (cls[y] == 1)
        out[cat.g(x), cat.g(y)] = out[cat.g(y), cat.g(x)] = np.count_nonzero(both)
    return out


def _no_outgroup_present(cat: Catalog) -> np.ndarray:
    """(n,) bool, merged states: no outgroup is P/S/T (each E or missing; all True without an outgroup)."""
    keep = np.ones(len(cat.lid), dtype=bool)
    for o in cat.outgroups:
        keep &= cat.cls(o) != 1
    return keep


def neighbor_joining(d: np.ndarray, labels: list[str]) -> Node:
    """Neighbor joining of `labels` on the symmetric (G, G) distances `d` (zero diagonal, finite): Saitou &
    Nei's Q-matrix, Q(i, j) = (r - 2) d(i, j) - sum_k d(i, k) - sum_k d(j, k) over the r active nodes; the
    minimum pair joins, ties to the first pair in sheet order (a joined node takes the place of its first
    member); negative branch lengths are set to 0. Unrooted: the top node has 3 children (2 for 2 leaves);
    children are in sheet order. A ValueError on fewer than two labels or a matrix that does not fit them."""
    d = np.array(d, dtype=np.float64)  # a copy: the joins rewrite it
    labels = list(labels)
    m = len(labels)
    if m < 2 or len(set(labels)) != m or d.shape != (m, m):
        raise ValueError(f"neighbor_joining: needs at least two distinct labels and their square matrix, not "
                         f"{m} label(s) and shape {d.shape}")
    if not np.isfinite(d).all():
        raise ValueError("neighbor_joining: every distance must be finite (nan: no locus resolved in both)")
    if not (np.allclose(d, d.T) and np.allclose(np.diag(d), 0.0)):
        raise ValueError("neighbor_joining: the distances must be symmetric with a zero diagonal")
    d = (d + d.T) / 2  # exact symmetry and diagonal, past rounding
    np.fill_diagonal(d, 0.0)
    nodes = [Node(x) for x in labels]
    while len(nodes) > 3:
        r = len(nodes)
        s = d.sum(axis=1)
        i, j = _first_min_pair((r - 2) * d - s[:, None] - s[None, :])
        li = d[i, j] / 2 + (s[i] - s[j]) / (2 * (r - 2))
        joined = Node(None, (replace(nodes[i], length=_nonneg(li)),
                             replace(nodes[j], length=_nonneg(d[i, j] - li))))
        du = (d[i] + d[j] - d[i, j]) / 2
        d[i, :] = du
        d[:, i] = du
        d[i, i] = 0.0
        keep = np.arange(r) != j
        d = d[np.ix_(keep, keep)]
        nodes[i] = joined
        del nodes[j]
    if len(nodes) == 3:
        ab, ac, bc = d[0, 1], d[0, 2], d[1, 2]
        lengths = ((ab + ac - bc) / 2, (ab + bc - ac) / 2, (ac + bc - ab) / 2)
    else:
        lengths = (d[0, 1] / 2, d[0, 1] / 2)
    return Node(None, tuple(replace(x, length=_nonneg(v)) for x, v in zip(nodes, lengths, strict=True)))


def _first_min_pair(q: np.ndarray) -> tuple[int, int]:
    """The first (i, j), i < j in row-major (sheet) order, whose q is the minimum, to within rounding."""
    iu, ju = np.triu_indices(len(q), 1)
    v = q[iu, ju]
    best = float(v.min())
    k = int(np.flatnonzero(v <= best + _TIE * max(1.0, abs(best)))[0])
    return int(iu[k]), int(ju[k])


def _nonneg(v: float) -> float:
    return float(v) if v > 0 else 0.0


def root_on(tree: Node, outgroup: str) -> Node:
    """`tree` rooted on the outgroup's branch: (outgroup, rest), the root at the branch's midpoint, so every
    leaf-to-leaf path keeps its length. A node of degree 2 (an earlier root) is merged into its branch;
    children are ordered by their first leaf in `tree`'s leaf order. A ValueError unless `outgroup` is a leaf
    of a tree with at least two distinct leaves and at least two children at its top."""
    rooted = _root_on_split(tree, frozenset((outgroup,)))
    if rooted is None:
        raise ValueError(f"root_on: {outgroup!r} is not a leaf of the tree "
                         f"({', '.join(sorted(tree.leaves()))})")
    return rooted


def _root_on_split(tree: Node, side: frozenset[str]) -> Node | None:
    """`tree` rooted at the midpoint of the branch that splits the leaves `side` from the others: (side,
    rest), as `root_on` roots on one leaf's branch; None when no branch does (`side` is not a clade of the
    unrooted tree). A ValueError on a tree `root_on` refuses."""
    nodes = list(_walk(tree))
    names = [x.name for x in nodes if not x.children]
    if None in names or len(set(names)) != len(names) or len(names) < 2 or len(tree.children) < 2:
        raise ValueError(f"root_on: a tree with at least two distinct leaf labels and two children at its "
                         f"top is needed, not leaves {names}")
    every = frozenset(names)
    adj: dict[int, list[tuple[Node, float]]] = {id(x): [] for x in nodes}
    branch = None  # (the node on side's end of the branch, the node on the other end, the branch's length)
    for x in nodes:
        for c in x.children:
            adj[id(x)].append((c, c.length))
            adj[id(c)].append((x, c.length))
            if branch is None:
                below = c.leaves()
                if below == side:
                    branch = (c, x, c.length)
                elif every - below == side:
                    branch = (x, c, c.length)
    if branch is None:
        return None
    near, far, length = branch
    near_end, near_from, length = _extend(near, far, length, adj)
    far_end, far_from, length = _extend(far, near, length, adj)
    rank = {name: k for k, name in enumerate(names)}
    return Node(None, (_hang(near_end, near_from, length / 2, adj, rank),
                       _hang(far_end, far_from, length / 2, adj, rank)))


def _extend(x: Node, toward: Node, length: float, adj: dict[int, list[tuple[Node, float]]]
            ) -> tuple[Node, Node, float]:
    """From x, at one end of a branch whose other end is `toward`, walk away through nodes of degree 2 (an
    earlier root): the first node of another degree, its neighbor along the branch, the branch's length."""
    while len(adj[id(x)]) == 2:
        y, w = next(e for e in adj[id(x)] if e[0] is not toward)
        toward, x, length = x, y, length + w
    return x, toward, length


def _hang(x: Node, parent: Node, length: float, adj: dict[int, list[tuple[Node, float]]],
          rank: dict[str, int]) -> Node:
    """x hung below `parent` on a branch of `length`: its other neighbors become its children (ordered by
    their first leaf); a node left with one child is merged into that child's branch."""
    kids = [_hang(y, x, w, adj, rank) for y, w in adj[id(x)] if y is not parent]
    if not kids:
        return Node(x.name, (), length)
    if len(kids) == 1:
        return replace(kids[0], length=kids[0].length + length)
    kids.sort(key=lambda k: min(rank[name] for name in k.leaves()))
    return Node(x.name, tuple(kids), length)


def _walk(t: Node):
    """Pre-order: t, then each child's subtree in order."""
    yield t
    for c in t.children:
        yield from _walk(c)


# ---------------------------------------------------------------- discovery


def fallback_reasons(excluded: dict[str, str]) -> str:
    """The fallback genomes left out and why (`eligible_fallback`), as one text."""
    return ("; ".join(f"{w}: {why}" for w, why in excluded.items())
            or "no annotated ingroup genome outside the trio")


def four_point_ok(dist: np.ndarray, cat: Catalog, r: str, s: str, x: str, d: str) -> bool:
    """The four-point condition d(r,s) + d(x,d) < min(d(r,x) + d(s,d), d(r,d) + d(s,x)): x splits with d
    from r and s. `cat` maps a label to its index (`cat.g`); a nan distance gives False."""
    rs_xd, rx_sd, rd_sx = _four_point_sums(dist, cat, r, s, x, d)
    return bool(rs_xd < rx_sd and rs_xd < rd_sx)


def _four_point_sums(dist: np.ndarray, cat: Catalog, r: str, s: str, x: str, d: str
                     ) -> tuple[float, float, float]:
    """The three four-point sums of (r, s, x, d): rs|xd, rx|sd, rd|sx."""
    gr, gs, gx, gd = (cat.g(k) for k in (r, s, x, d))
    return (float(dist[gr, gs] + dist[gx, gd]), float(dist[gr, gx] + dist[gs, gd]),
            float(dist[gr, gd] + dist[gs, gx]))


def eligible_fallback(cat: Catalog, x: str, y: str, z: str, dist: np.ndarray
                      ) -> tuple[tuple[str, ...], dict[str, str]]:
    """The fallback discovery genomes for (x, y; z), used when z has no annotation: the annotated ingroup
    genomes w outside {x, y, z} that split with z from x and y, `four_point_ok(dist, cat, x, y, w, z)`, in
    sheet order. Terminal-branch effects of annotation depth cancel in each four-point sum. The other
    candidates -> why: "groups with <x>" or "groups with <y>" (the smaller of the other two sums pairs w with
    that genome), or "distance unresolved" (a nan distance, or a tie)."""
    dist = np.asarray(dist, dtype=np.float64)
    G = len(cat.labels)
    if dist.shape != (G, G):
        raise ValueError(f"eligible_fallback: dist has shape {dist.shape}; the run has {G} genomes")
    ok, excluded = [], {}
    for w in cat.ingroup:
        if w in (x, y, z) or not cat.annotated[w]:
            continue
        if four_point_ok(dist, cat, x, y, w, z):
            ok.append(w)
            continue
        sums = _four_point_sums(dist, cat, x, y, w, z)
        if not all(map(math.isfinite, sums)) or sums[1] == sums[2]:
            excluded[w] = "distance unresolved"
        else:
            excluded[w] = f"groups with {x}" if sums[1] < sums[2] else f"groups with {y}"
    return tuple(ok), excluded


def discover(cat: Catalog, x: str, y: str, z: str, dist: np.ndarray, outgroup: str) -> Counted | None:
    """The loci a statistic with a discovery rule counts for x, y and the discovering side z,
    each state read from the row's discovering genome's own calls (`pcls_by`); the one discovery function.

    disc = (z,) if z is annotated, else `eligible_fallback`. Per row, gstar = z where z is a member; under
    fallback, the first eligible genome (sheet order) that is a member there. cand = gstar exists, z present
    in gstar's call, and the outgroup rule on gstar's call into every outgroup o: "E" o is E (D_pav),
    "not_present" o is not P/S/T (E or missing: the KKSC test's set-aside), "none" no rule. rows = cand with x
    and y resolved (P/S/T or E). None when no genome may discover (z unannotated, no eligible fallback): the
    caller writes "skipped" with a WARNING naming the excluded genomes and why (`eligible_fallback`). x, y and
    z are distinct ingroup genomes."""
    if outgroup not in OUTGROUP_RULES:
        raise ValueError(f"discover: outgroup rule {outgroup!r} is not one of {', '.join(OUTGROUP_RULES)}")
    check_trio(cat, (x, y, z), "discover")
    if cat.annotated[z]:
        disc, excluded = (z,), {}
    else:
        disc, excluded = eligible_fallback(cat, x, y, z, dist)
        if not disc:
            return None
    gstar = first_member(cat, disc)
    cand = (gstar >= 0) & (cat.pcls_by(gstar, z) == 1)
    if outgroup != "none":
        for o in cat.outgroups:
            po = cat.pcls_by(gstar, o)
            cand &= (po == 0) if outgroup == "E" else (po != 1)
    rows = cand & (cat.pcls_by(gstar, x) >= 0) & (cat.pcls_by(gstar, y) >= 0)
    loci = f"disc={','.join(disc)};present={z};res={x},{y};outgroup={outgroup}"
    return Counted(rows=rows, cand=cand, gstar=gstar, disc=disc, loci=loci, excluded=excluded)


# ---------------------------------------------------------------- sharing, sister pairs, KKSC


def trio_sharing(cat: Catalog, genomes: tuple[str, str, str]) -> tuple[dict[str, int], int]:
    """A trio's sharing counts, which only pick its sister pair, on merged states (any
    discovery): over the loci resolved in all three with no outgroup present, "x+y" -> the loci where x and y
    are present and the third E, for the three pairs in sheet order. Returns (counts, loci counted)."""
    a, b, c = _sheet_trio(cat, genomes, "trio_sharing")
    cls = {x: cat.cls(x) for x in (a, b, c)}
    keep = _no_outgroup_present(cat) & (cls[a] >= 0) & (cls[b] >= 0) & (cls[c] >= 0)
    counts = {f"{x}+{y}": int(np.count_nonzero(keep & (cls[x] == 1) & (cls[y] == 1) & (cls[t] == 0)))
              for x, y, t in ((a, b, c), (a, c, b), (b, c, a))}
    return counts, int(np.count_nonzero(keep))


def sister_pair(cat: Catalog, genomes: tuple[str, str, str], dist: np.ndarray) -> SisterCall:
    """The sister pair of three ingroup genomes (any order; sheet order inside) and the KKSC symmetry test.

    pair = the largest sharing count (`trio_sharing`; ties to the first pair in sheet order), donor = the
    third genome z. p_top = two-sided binomtest(n_top, n_top + n_second); ambiguous when p_top >= P_WARN, or
    when no pair shares an insertion (p_top nan). p_kksc = two-sided binomtest(n_x+z, n_x+z + n_y+z) on
    `discover(cat, x, y, z, dist, "not_present")` (z's own calls; the fallback when z has no annotation):
    n_x+z its rows with x present and y E, n_y+z the mirror; kksc_loci its `loci`, or "none" (p_kksc nan)
    when `discover` returns None. `relations` writes the WARNING an ambiguous pair needs."""
    return _sister(cat, genomes, dist).call


@dataclass
class _Sharing:
    """A trio's sister call with what its rows and findings need beyond SisterCall."""

    call: SisterCall
    scope: str                  # "a,b,c", sheet order
    n: int                      # loci the sharing counts take
    second: str                 # the second pair, "x+y"
    kksc: Counted | None
    n_xz: int
    n_yz: int
    resolved: dict[str, float]  # x, y resolved (P/S/T or E) among the KKSC candidates, from gstar's calls
    excluded: dict[str, str]    # fallback genomes left out


def _sister(cat: Catalog, genomes: tuple[str, str, str], dist: np.ndarray) -> _Sharing:
    trio = _sheet_trio(cat, genomes, "sister_pair")
    counts, n = trio_sharing(cat, trio)
    pairs = list(itertools.combinations(trio, 2))  # (a, b), (a, c), (b, c): trio_sharing's order
    ranked = sorted(pairs, key=lambda p: -counts[f"{p[0]}+{p[1]}"])  # stable: ties keep sheet order
    (x, y), second = ranked[0], f"{ranked[1][0]}+{ranked[1][1]}"
    z = next(g for g in trio if g not in (x, y))
    n_top = counts[f"{x}+{y}"]
    p_top = _binom_p(n_top, n_top + counts[second])
    kk = discover(cat, x, y, z, dist, "not_present")
    n_xz = n_yz = 0
    resolved: dict[str, float] = {}
    if kk is None:
        p_kksc, kksc_loci, excluded = math.nan, "none", eligible_fallback(cat, x, y, z, dist)[1]
    else:
        px, py = cat.pcls_by(kk.gstar, x), cat.pcls_by(kk.gstar, y)
        n_xz = int(np.count_nonzero(kk.rows & (px == 1) & (py == 0)))
        n_yz = int(np.count_nonzero(kk.rows & (py == 1) & (px == 0)))
        p_kksc, kksc_loci, excluded = _binom_p(n_xz, n_xz + n_yz), kk.loci, kk.excluded
        n_cand = int(np.count_nonzero(kk.cand))
        if n_cand:
            resolved = {g: np.count_nonzero(p[kk.cand] >= 0) / n_cand for g, p in ((x, px), (y, py))}
    call = SisterCall(pair=(x, y), donor=z, counts=counts, p_top=p_top, p_kksc=p_kksc, kksc_loci=kksc_loci,
                      ambiguous=not p_top < P_WARN)  # nan (nothing shared) is ambiguous too
    return _Sharing(call, ",".join(trio), n, second, kk, n_xz, n_yz, resolved, excluded)


def _sheet_trio(cat: Catalog, genomes: tuple[str, str, str], who: str) -> tuple[str, ...]:
    """Three distinct ingroup genomes, in sheet order."""
    genomes = tuple(genomes)
    if len(genomes) != 3:
        raise AnalysisError(f"{who}: a trio is three genomes, not {', '.join(genomes)}")
    check_trio(cat, genomes, who)
    return tuple(sorted(genomes, key=cat.g))


def _binom_p(k: int, n: int) -> float:
    """Two-sided exact binomial p of k in n trials at 1/2; nan when n is 0."""
    if n == 0:
        return math.nan
    from scipy.stats import binomtest
    return float(binomtest(k, n).pvalue)


# ---------------------------------------------------------------- the section


def relations(cat: Catalog) -> tuple[list[StatRow], Findings, Node | None, np.ndarray]:
    """The relationships of a run: rows `distance` and `shared_derived` per pair, `tree`, per
    trio of ingroup genomes `sharing` (n_x+y for the three pairs, sister_pair, p_top, p_kksc) and
    `kksc_resolution`, and `outgroup_conflict` per ingroup genome, with their findings (section
    "Relationships"). Returns (rows, findings, tree, distances): the tree is rooted between the outgroups and
    the ingroup (module docstring), unrooted without an outgroup, None when skipped (fewer than three genomes,
    or a nan distance)."""
    f = Findings()
    dist, n_res = pav_distances(cat)
    rows = [StatRow("pav", "distance", f"{x}~{y}", "all", "d", float(dist[a, b]), n=int(n_res[a, b]),
                    loci=f"disc=any;res={x},{y}")
            for (a, x), (b, y) in itertools.combinations(enumerate(cat.labels), 2)]
    sd = shared_derived(cat)
    rows += [StatRow("pav", "shared_derived", f"{x}~{y}", "all", "n", int(sd[cat.g(x), cat.g(y)]),
                     loci=f"disc=any;present={x},{y};absent=another ingroup genome;outgroup=not_present")
             for x, y in itertools.combinations(cat.ingroup, 2)]
    tree, tree_rows = _tree(cat, dist, f)
    rows += tree_rows
    if len(cat.ingroup) < 3:
        f.add(SECTION, "Sharing counts, sister pairs and the KKSC symmetry test skipped: fewer than three "
                       "ingroup genomes")
    for trio in itertools.combinations(cat.ingroup, 3):
        rows += _trio_rows(_sister(cat, trio, dist), f)
    rows += _conflicts(cat, f)
    return rows, f, tree, dist


def _tree(cat: Catalog, dist: np.ndarray, f: Findings) -> tuple[Node | None, list[StatRow]]:
    """The NJ tree on every genome, rooted between the outgroups and the ingroup (module docstring), with its
    row and findings line; or None, no row and the line saying why it is skipped."""
    G = len(cat.labels)
    if G < 3:
        f.add(SECTION, f"Tree skipped: fewer than three genomes ({G})")
        return None, []
    bad = [(cat.labels[a], cat.labels[b]) for a, b in itertools.combinations(range(G), 2)
           if not math.isfinite(dist[a, b])]
    if bad:
        more = f" (and {len(bad) - 1} more pair(s))" if len(bad) > 1 else ""
        f.add(SECTION, f"Tree skipped: no locus is resolved in both {bad[0][0]} and {bad[0][1]}, so their "
                       f"PAV distance is undefined{more}")
        return None, []
    tree = neighbor_joining(dist, list(cat.labels))
    outs = cat.outgroups
    rooted = _root_on_split(tree, frozenset(outs)) if len(outs) > 1 else None
    if not outs:
        how = "unrooted (no outgroup)"
    elif len(outs) == 1:
        tree, how = root_on(tree, outs[0]), f"rooted on the outgroup {outs[0]}"
    elif rooted is not None:
        tree, how = rooted, f"rooted between the outgroups ({', '.join(outs)}) and the ingroup"
    else:
        tree, how = root_on(tree, outs[0]), f"rooted on the first outgroup, {outs[0]}"
        f.warn(SECTION, f"the outgroups {', '.join(outs)} do not form a clade in the PAV tree, so it is "
                        f"rooted on the first, {outs[0]}: check their roles in the samples sheet")
    text = tree.newick()
    f.add(SECTION, f"Neighbor-joining tree on PAV distances (merged states, any discovery), {how}: {text}")
    return tree, [StatRow("pav", "tree", "all", "all", "newick", text, loci="disc=any;res=each pair",
                          note=f"neighbor joining on the PAV distances (rows distance); {how}")]


def _trio_rows(sh: _Sharing, f: Findings) -> list[StatRow]:
    """One trio's `sharing` and `kksc_resolution` rows, its findings lines and the WARNINGs: an ambiguous
    sister pair, a KKSC test without a discovery genome."""
    sc, scope = sh.call, sh.scope
    (x, y), z = sc.pair, sc.donor
    top, n_top, n_second = f"{x}+{y}", sc.counts[f"{x}+{y}"], sc.counts[sh.second]
    loci = f"disc=any;res={scope};outgroup=not_present"
    rows = [StatRow("pav", "sharing", scope, "all", f"n_{k}", v, n=sh.n, loci=loci)
            for k, v in sc.counts.items()]
    rows += [StatRow("pav", "sharing", scope, "all", "sister_pair", top, n=sh.n, loci=loci,
                     note=f"the pair present together most often with the third genome E; third genome {z}"
                          + ("; ambiguous" if sc.ambiguous else "")),
             StatRow("pav", "sharing", scope, "all", "p_top", sc.p_top, n=n_top + n_second, loci=loci,
                     note=f"two-sided binomial test, {top} {n_top} against {sh.second} {n_second}")]
    shared = ", ".join(f"{k} {v}" for k, v in sc.counts.items())
    f.add(SECTION, f"Trio {scope}: sister pair {top}, third genome {z} (insertions two share with the third "
                   f"E, merged states, on {sh.n} loci resolved in all three with no outgroup present: "
                   f"{shared}; p_top {fmt(sc.p_top, '.3g')})")
    if sc.ambiguous:
        why = (f"no two of the three share an insertion the third lacks, so {top} is only the first pair in "
               "sheet order" if math.isnan(sc.p_top) else
               f"{top} {n_top} against {sh.second} {n_second}: p_top {fmt(sc.p_top, '.3g')} >= {P_WARN}, two-sided "
               "binomial")
        f.warn(SECTION, f"trio {scope}: the sister pair {top} is ambiguous ({why}); annotation depth can "
                        "tilt a near-tie, so the roles of this trio are uncertain")
    if sh.kksc is None:
        rows.append(StatRow("pav", "sharing", scope, "all", "p_kksc", math.nan, n=0, loci="none",
                            note=f"skipped: {z} has no annotation and no eligible fallback discovery genome "
                                 f"({fallback_reasons(sh.excluded)})"))
        f.warn(SECTION, f"KKSC symmetry test of {x},{y};{z} skipped: {z} has no annotation and no annotated "
                        f"ingroup genome outside the trio splits with it from {x} and {y} by the four-point "
                        f"condition ({fallback_reasons(sh.excluded)})")
        return rows
    disc = ",".join(sh.kksc.disc)
    rows.append(StatRow("pav", "sharing", scope, "all", "p_kksc", sc.p_kksc, n=sh.n_xz + sh.n_yz,
                        loci=sc.kksc_loci,
                        note=f"n_{x}+{z} {sh.n_xz}, n_{y}+{z} {sh.n_yz} from {disc}'s own calls; "
                             f"{KKSC_NOTE}"))
    unequal = f"; the two counts differ at p < {P_WARN}" if sc.p_kksc < P_WARN else ""
    fallback = "" if sh.kksc.disc == (z,) else (
        f" ({z} has no annotation; {disc} split with it from {x} and {y} by the four-point condition"
        + (f"; not eligible: {fallback_reasons(sh.excluded)}" if sh.excluded else "") + ")")
    f.add(SECTION, f"KKSC symmetry of {x},{y};{z} (equal minority counts under lineage sorting), on loci "
                   f"{disc} discovers{fallback} with {z} present, no outgroup present and {x}, {y} resolved, "
                   f"from {disc}'s own calls: n_{x}+{z} {sh.n_xz}, n_{y}+{z} {sh.n_yz}, p_kksc "
                   f"{fmt(sc.p_kksc, '.3g')}{unequal}; {KKSC_NOTE}")
    if sh.resolved:  # resolution of x and y among the KKSC candidates; its WARNING follows the line above
        rows += resolution("kksc_resolution", scope, sh.resolved, int(np.count_nonzero(sh.kksc.cand)),
                           f"disc={disc};present={z};outgroup=not_present", f, SECTION)
    return rows


def _conflicts(cat: Catalog, f: Findings) -> list[StatRow]:
    """Rows `outgroup_conflict` (statistic frac_E) per ingroup genome g: an error-rate estimate,
    reported as an upper bound, from the discoverers' own calls.

    Over the loci discovered by an annotated ingroup genome other than g (gstar: the first member among them,
    sheet order) where gstar's own call puts some outgroup present (P/S/T) and g resolved, the fraction where
    gstar calls g E. g never discovers its own count: a locus only g discovered (P there by construction)
    never enters, and a locus another genome also discovered enters with g's state from that genome's call,
    so g's annotation depth cannot dilute the fraction. Annotated outgroups never discover here: the locus
    must be one an ingroup genome carries too, since an insertion private to an outgroup leaves every ingroup
    genome E without any error."""
    if not cat.outgroups:
        f.add(SECTION, "No outgroup: old insertions are not set aside from the sharing and shared-derived "
                       "counts, and the outgroup-conflict bound is not computed")
        return []
    rows, parts = [], []
    for g in cat.ingroup:
        disc = [a for a in cat.ingroup if a != g and cat.annotated[a]]
        gstar = first_member(cat, disc)
        some = np.zeros(len(cat.lid), dtype=bool)
        for o in cat.outgroups:
            some |= cat.pcls_by(gstar, o) == 1
        c = cat.pcls_by(gstar, g)
        m = (gstar >= 0) & some & (c >= 0)
        n = int(np.count_nonzero(m))
        e = int(np.count_nonzero(m & (c == 0)))
        rows.append(StatRow("pav", "outgroup_conflict", g, "all", "frac_E", e / n if n else math.nan, n=n,
                            loci=f"disc={','.join(disc)};res={g};outgroup=present" if disc else "none",
                            note=CONFLICT_NOTE))
        parts.append(f"{g} {e} of {n}" if disc else f"{g}: no other annotated ingroup genome discovers")
    f.add(SECTION, f"Ingroup E at loci an outgroup carries, each genome read from the other ingroup genomes' "
                   f"discoveries and own calls: {', '.join(parts)}; {CONFLICT_NOTE}")
    return rows
