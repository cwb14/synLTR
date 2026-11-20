"""Inverted blocks between genome pairs: stretches of resolved loci whose order on one genome is
reversed on the other, found from the loci's merged coordinates (any discovery: only positions are used), and
only between homologous sequences.

Homologous sequences (`syntenic_pairs`, the one definition: the overview draws each recipient sequence above
its counterpart by the same rule). For genomes x and y, the loci both place (merged sites, any discovery: a
state with a site in each) are counted per (x sequence, y sequence), keyed on fai indexes, never on names, so
two genomes that both call a sequence "chr2" never collide. A y sequence t is x sequence s's dominant partner
when it holds most of the loci both place on s (ties to the lower fai index), and alike from y's side. (s, t)
are syntenic partners when t is s's dominant partner or s is t's (reciprocal dominance): every sequence is
paired with its dominant partner, and a translocated piece keeps its pair when it is the larger share of
either of its two sequences. A handful of loci placed on another sequence (a paralog anchored elsewhere, a
local misassembly) never makes a pair, so they are never searched: an inversion is a reversal between
homologous sequences, and their reversed order is noise. The cost: a piece that is the smaller share of both
its sequences (say, two short arms swapped by a reciprocal translocation) is not searched, so an inversion
inside it is missed.

For each pair (x, y), x before y in sheet order (outgroups included: they polarize arrangements), the loci
resolved (P/S/T or E) in both with coordinates in both are grouped by (x sequence, y sequence), syntenic
partner pairs only. In a group, sorted by x, a run is a maximal stretch of rows whose y steps go against the
group's dominant direction (equal y joins without breaking it); runs of >= 2 rows are candidates;
consecutive candidates merge across a few discordant rows (`_find`: u + v + the rows between < min(|A|, |B|)); a merged block needs
INV_MIN_LOCI rows of distinct y outside the discordant ones. INV_MIN_LOCI is a heuristic (named in
findings.txt). A whole reversed sequence (every step against the sign of the sum) is not an inversion: the
dominant direction is the reference. A locus placed singly among a group's collinear loci makes at most a
two-row candidate, far below INV_MIN_LOCI, so loci misplaced one by one within a homologous pair form no block.

Contiguity on both genomes: an inverted segment is one stretch of each genome, homologous as a whole to one
stretch of the other, so inside a block's extent on x and on y the pair's loci placed elsewhere (resolved in
both, any discovery) must be fewer than the block's loci outside its discordant ones (INV_FOREIGN_MAX times
them: the merge rule's own tolerance, an interruption shorter than what it interrupts). Else the run is
dropped, and counted in findings.txt. Never counted: loci of the block's own sequence pair (the same homologous
region in another order: an inversion covering most of a pair makes its collinear remainder the block, whose y
extent spans the inverted rows), and, when the pair's sequences are each other's dominant partners, loci on
another syntenic partner of either (an alternate haplotig or a duplicated region holds homologous sequence, not
misplaced loci). Reciprocal dominance alone keeps one sparse artifact, which this removes: a
scaffold whose few loci lie spread along the chromosome that dominates it, a pair kept from one side only, so
every other sequence counts (its tests). The rule's cost: a real inversion whose region's second
copy lies on a non-partner sequence, or one on a pair kept from one side only that holds more loci of other
sequences in its extent than its own, is dropped.

Breakpoint intervals are bounded by the loci just outside the block: on x by the rows before and after it;
on y by the rows bounding the block's y extent from below and above (the neighbors in x order are not those
when the block is reversed). None where the block meets a group's end.

The library route (route property inversions "chain_blocks"; `chain_inversions`). The screen's collinear chain blocks
(blocks.tsv: every genome pair's blocks, pair order then block order, each with its sign within its chromosome pair,
its anchors and its extents, anchor midpoints, 1-based) are grouped per genome pair and chromosome pair, in order of
first appearance. A chromosome pair's main orientation is the sign holding more anchors (a tie: +); a pair with fewer
than CHAIN_INV_MIN anchors in all is not searched; every block of the other sign is an inverted block. The blocks of
every pair are sorted by anchors, most first (stable). They are written in INV_COLUMNS: x and y the pair's genomes,
the block's extents on each, n_loci its anchors; breakpoints, n_strand and n_discordant "." (a chain block has
neither). A block's extent runs from its first to its last anchor's midpoint, so it counts in n_blocks_1Mb when that
distance, end - start, is >= 1 Mb (the assembly route: a block's span end - start + 1). Painted elements inside one
are flagged as on the assembly route (`inside`, by the sequence codes of LibraryRun.sequences).
"""
from __future__ import annotations

import csv
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ltrtrace.analyses.geneflow.catalog import AnalysisError, Catalog
from ltrtrace.base.constants import CHAIN_INV_MIN, CHAIN_MIN, INV_FOREIGN_MAX, INV_MIN_LOCI
from ltrtrace.stats.report import Findings, StatRow

INV_COLUMNS = ("genome_x", "seqid_x", "start_x", "end_x", "bp_left_x", "bp_right_x",
               "genome_y", "seqid_y", "start_y", "end_y", "bp_left_y", "bp_right_y", "n_loci", "n_strand",
               "n_discordant")
SECTION = "Inversions"
MB = 1_000_000


@dataclass(frozen=True)
class InvBlock:
    x: str
    y: str
    seq_x: int
    seq_y: int  # genomes and sequence codes (fai indexes; on the library route LibraryRun.sequences): the keys
    start_x: int
    end_x: int
    start_y: int
    end_y: int  # extent of the block's loci, 1-based
    bp_left_x: tuple[int, int] | None
    bp_right_x: tuple[int, int] | None  # breakpoint intervals; None at a sequence end
    bp_left_y: tuple[int, int] | None
    bp_right_y: tuple[int, int] | None
    n_loci: int
    n_strand: int | None
    n_discordant: int | None  # block loci; of them, those whose strands also say inverted; rows merged over
    dominant: int = 1  # the group's dominant y direction (+1, or -1: the block is the collinear remainder)


def _find(y: np.ndarray) -> list[tuple[np.ndarray, int, np.ndarray, int]]:
    """One (seq_x, seq_y) group's mid_y in x order -> (row indexes, n_discordant, rows outside the discordant
    ones, dominant) of each kept block.

    step_i = sign(y[i+1] - y[i]); the dominant direction is +1 if the steps sum to >= 0, else -1. A run is a
    maximal stretch of rows joined by steps equal to -dominant (0-steps join without breaking it) with at
    least one such step. Consecutive candidates A, B (x order; a merged block is the next A) merge when the
    fewest rows u trimmed from A's end and v from B's start leave A's last and B's first remaining rows in
    reversed order, with u + v + (rows between) < min(|A|, |B|); those rows are the block's discordant rows.
    A merged block is kept with >= INV_MIN_LOCI rows of distinct y outside the discordant rows."""
    y = np.asarray(y, dtype=float)
    if len(y) < 2:
        return []
    step = np.sign(np.diff(y)).astype(int)
    dom = 1 if step.sum() >= 0 else -1
    join = (step == -dom) | (step == 0)
    cands: list[list[int]] = []
    i = 0
    while i < len(step):
        if not join[i]:
            i += 1
            continue
        j = i
        while j + 1 < len(step) and join[j + 1]:
            j += 1
        if (step[i:j + 1] == -dom).any():
            cands.append([i, j + 1])
        i = j + 1
    bad = np.zeros(len(y), dtype=bool)  # discordant rows
    merged: list[list[int]] = []
    for b in cands:
        if merged and _merge(y, dom, merged[-1], b, bad):
            merged[-1][1] = b[1]
        else:
            merged.append(list(b))
    out = []
    for lo, hi in merged:
        rows = np.arange(lo, hi + 1)
        keep = rows[~bad[rows]]
        if len(np.unique(y[keep])) >= INV_MIN_LOCI:
            out.append((rows, int(bad[rows].sum()), keep, dom))
    return out


def _merge(y: np.ndarray, dom: int, a: list[int], b: list[int], bad: np.ndarray) -> bool:
    """Merge candidate b into a (the merge rule in `_find`'s docstring); marks the discordant rows in `bad`."""
    na, nb, gap = a[1] - a[0] + 1, b[1] - b[0] + 1, b[0] - a[1] - 1
    for t in range(max(min(na, nb) - gap, 0)):  # u + v = t < min(|A|, |B|) - gap
        for u in range(t + 1):
            v = t - u
            if u >= na or v >= nb:
                continue
            ya, yb = y[a[1] - u], y[b[0] + v]
            if (yb < ya) if dom == 1 else (yb > ya):
                bad[a[1] - u + 1:b[0] + v] = True  # A's trimmed end, the rows between, B's trimmed start
                return True
    return False


def syntenic_pairs(cat: Catalog, x: str, y: str) -> dict[tuple[int, int], int]:
    """The syntenic partner pairs of genome x's and genome y's sequences (module docstring: t is s's dominant
    partner, or s is t's) -> the loci both genomes place on the pair (merged sites, any discovery); keyed
    (x fai index, y fai index), in that order. {} when no locus is placed in both."""
    return _dominance(cat, x, y)[0]


def _dominance(cat: Catalog, x: str, y: str) -> tuple[dict[tuple[int, int], int], dict[int, int], dict[int, int]]:
    """(`syntenic_pairs`, each x sequence's dominant y partner, each y sequence's dominant x partner)."""
    gx, gy = cat.g(x), cat.g(y)
    both = (cat.seq[:, gx] >= 0) & (cat.seq[:, gy] >= 0)
    sx, sy = cat.seq[both, gx].astype(np.int64), cat.seq[both, gy].astype(np.int64)
    ny = max(len(cat.lengths[y]), 1)
    keys, counts = np.unique(sx * ny + sy, return_counts=True)   # sorted by s, then t
    pairs = [(int(k) // ny, int(k) % ny, int(c)) for k, c in zip(keys.tolist(), counts.tolist(), strict=True)]
    dom_x: dict[int, tuple[int, int]] = {}  # s -> (count, t): the first maximum is the lowest t
    dom_y: dict[int, tuple[int, int]] = {}  # t -> (count, s): s comes in ascending order, so likewise
    for s, t, c in pairs:
        if s not in dom_x or c > dom_x[s][0]:
            dom_x[s] = (c, t)
        if t not in dom_y or c > dom_y[t][0]:
            dom_y[t] = (c, s)
    partners = {(s, t): c for s, t, c in pairs if dom_x[s][1] == t or dom_y[t][1] == s}
    return partners, {s: t for s, (_, t) in dom_x.items()}, {t: s for t, (_, s) in dom_y.items()}


def _pair_loci(cat: Catalog, x: str, y: str) -> tuple[np.ndarray, np.ndarray]:
    """The loci a pair's search reads (resolved in both, coordinates in both; catalog order) and, of them, those
    on a syntenic partner pair (bool, aligned)."""
    gx, gy = cat.g(x), cat.g(y)
    idx = np.flatnonzero((cat.cls(x) >= 0) & (cat.cls(y) >= 0) & (cat.seq[:, gx] >= 0) & (cat.seq[:, gy] >= 0))
    ny = max(len(cat.lengths[y]), 1)
    keys = cat.seq[idx, gx].astype(np.int64) * ny + cat.seq[idx, gy].astype(np.int64)
    partners = np.array([s * ny + t for s, t in syntenic_pairs(cat, x, y)], dtype=np.int64)
    return idx, np.isin(keys, partners)


def inverted_blocks(cat: Catalog) -> list[InvBlock]:
    """Every inverted block of every pair (x before y in sheet order), by pair then x position; groups of syntenic
    partner sequences only, blocks contiguous on both genomes (module docstring)."""
    return _search(cat)[0]


def _search(cat: Catalog) -> tuple[list[InvBlock], dict[tuple[str, str], int]]:
    """`inverted_blocks`, and per pair the reversed runs dropped as not contiguous on both genomes."""
    out: list[InvBlock] = []
    dropped: dict[tuple[str, str], int] = {}
    labels = cat.labels
    for ix, x in enumerate(labels):
        gx = cat.g(x)
        mid_x = cat.mid(x)
        for y in labels[ix + 1:]:
            gy = cat.g(y)
            mid_y = cat.mid(y)
            every, on = _pair_loci(cat, x, y)
            idx = every[on]
            if len(idx) < INV_MIN_LOCI:
                continue
            partners, dom_x, dom_y = _dominance(cat, x, y)
            near = (_by_sequence(cat.seq[every, gx], mid_x[every], cat.seq[every, gy]),
                    _by_sequence(cat.seq[every, gy], mid_y[every], cat.seq[every, gx]))
            sx, sy = cat.seq[idx, gx], cat.seq[idx, gy]
            mx, my = mid_x[idx], mid_y[idx]
            order = np.lexsort((my, mx, sy, sx))
            idx, sx, sy, my = idx[order], sx[order], sy[order], my[order]
            cut = np.flatnonzero((np.diff(sx) != 0) | (np.diff(sy) != 0)) + 1
            for lo, hi in zip(np.r_[0, cut], np.r_[cut, len(idx)], strict=True):
                if hi - lo < INV_MIN_LOCI:
                    continue
                grp = idx[lo:hi]
                s, t = int(sx[lo]), int(sy[lo])
                # loci that do not count against the block: its own pair's, and, when s and t are each other's
                # dominant partners, those on another syntenic partner of s (or of t): an alternate haplotig or a
                # duplicated region holds homologous sequence, not misplaced loci
                mutual = dom_x.get(s) == t and dom_y.get(t) == s
                ok_y = {t2 for (s2, t2) in partners if s2 == s} if mutual else {t}
                ok_x = {s2 for (s2, t2) in partners if t2 == t} if mutual else {s}
                for rows, nd, good, dom in _find(my[lo:hi]):
                    blk = _block(cat, x, y, gx, gy, grp, rows, nd, s, t, good, dom)
                    foreign = (_count(near[0].get(s), blk.start_x, blk.end_x, ok_y),
                               _count(near[1].get(t), blk.start_y, blk.end_y, ok_x))
                    if max(foreign) < INV_FOREIGN_MAX * len(good):
                        out.append(blk)
                    else:
                        dropped[(x, y)] = dropped.get((x, y), 0) + 1
    return out, dropped


def _by_sequence(seq: np.ndarray, mid: np.ndarray, other: np.ndarray) -> dict[int, tuple[np.ndarray, np.ndarray]]:
    """Fai index -> (the sorted mids of the loci on that sequence, the other genome's sequence of each)."""
    out = {}
    for s in np.unique(seq):
        on = np.flatnonzero(seq == s)
        o = np.argsort(mid[on], kind="stable")
        out[int(s)] = (mid[on][o], other[on][o].astype(np.int64))
    return out


def _count(near: tuple[np.ndarray, np.ndarray] | None, lo: int, hi: int, ok: set[int]) -> int:
    """How many loci of `near` (sorted mids, the other genome's sequence of each) lie in [lo, hi] placed on an
    other-genome sequence outside `ok` (0 for None)."""
    if near is None:
        return 0
    mids, other = near
    a, b = np.searchsorted(mids, lo, "left"), np.searchsorted(mids, hi, "right")
    return int(np.count_nonzero(~np.isin(other[a:b], np.fromiter(ok, dtype=np.int64))))


def _block(cat: Catalog, x: str, y: str, gx: int, gy: int, grp: np.ndarray, rows: np.ndarray, nd: int,
           seq_x: int, seq_y: int, good: np.ndarray, dom: int) -> InvBlock:
    a, b = int(rows[0]), int(rows[-1])
    last = len(grp) - 1
    sub = grp[rows]
    ok = grp[good]  # the extents leave out displaced (discordant) rows, which may lie Mb away on y
    sxs, exs, sys_, eys = cat.start[:, gx], cat.end[:, gx], cat.start[:, gy], cat.end[:, gy]

    def span(xs, es, lo_row, hi_row):  # interval between the end of lo_row and the start of hi_row (group rows)
        if lo_row < 0 or hi_row < 0 or lo_row > last or hi_row > last:
            return None
        e, s = int(es[grp[lo_row]]), int(xs[grp[hi_row]])
        return (min(e, s), max(e, s))

    if dom == 1:  # the block's y extent runs from row b (low y) to row a (high y)
        left_y, right_y = span(sys_, eys, a - 1, b), span(sys_, eys, a, b + 1)
    else:         # from row a (low y) to row b (high y); the group's descending neighbors bound it
        left_y, right_y = span(sys_, eys, b + 1, a), span(sys_, eys, b, a - 1)
    strand = cat.strand[sub, gx].astype(int) * cat.strand[sub, gy].astype(int)
    return InvBlock(
        x, y, seq_x, seq_y,
        int(sxs[ok].min()), int(exs[ok].max()), int(sys_[ok].min()), int(eys[ok].max()),
        span(sxs, exs, a - 1, a), span(sxs, exs, b, b + 1), left_y, right_y,
        len(rows), int((strand == -dom).sum()), nd, dom)


def _iv(v: tuple[int, int] | None) -> str:
    return "." if v is None else f"{v[0]}-{v[1]}"


def write_inversions(blocks: list[InvBlock], cat: Catalog, path: Path) -> None:
    """inversions.tsv: INV_COLUMNS, one row per block; names from cat.seqids (a catalog's, or a library screen's
    LibraryRun.seqids); "." for None."""
    with open(path, "w") as fh:
        fh.write("\t".join(INV_COLUMNS) + "\n")
        for k in blocks:
            row = [k.x, cat.seqids[k.x][k.seq_x], k.start_x, k.end_x, _iv(k.bp_left_x), _iv(k.bp_right_x),
                   k.y, cat.seqids[k.y][k.seq_y], k.start_y, k.end_y, _iv(k.bp_left_y), _iv(k.bp_right_y),
                   k.n_loci, k.n_strand, k.n_discordant]
            fh.write("\t".join("." if v is None else str(v) for v in row) + "\n")


def inside(blocks: list[InvBlock], label: str, partners: set[str], seq: np.ndarray, pos: np.ndarray
           ) -> np.ndarray:
    """Bool per position: inside a block of `label` against a partner (seq: label's fai indexes, pos: 1-based
    positions)."""
    seq, pos = np.asarray(seq), np.asarray(pos)
    out = np.zeros(len(pos), dtype=bool)
    for k in blocks:
        if k.x == label and k.y in partners:
            s, lo, hi = k.seq_x, k.start_x, k.end_x
        elif k.y == label and k.x in partners:
            s, lo, hi = k.seq_y, k.start_y, k.end_y
        else:
            continue
        out |= (seq == s) & (pos >= lo) & (pos <= hi)
    return out


def inversion_rows(blocks: list[InvBlock], cat: Catalog) -> tuple[list[StatRow], Findings]:
    """Rows (test "inversions", scope x~y: n_blocks, n_blocks_1Mb, bp_inverted_x; n: the loci searched, those on
    syntenic partner pairs) and findings (section "Inversions"): each pair's three largest blocks (by x span) and
    the loci left unsearched off its partner pairs, the partner rule, the INV_MIN_LOCI heuristic, and what the
    method misses."""
    rows: list[StatRow] = []
    f = Findings()
    labels = cat.labels
    dropped = _search(cat)[1]
    for ix, x in enumerate(labels):
        for y in labels[ix + 1:]:
            mine = [k for k in blocks if k.x == x and k.y == y]
            idx, on = _pair_loci(cat, x, y)
            n, n_off = int(np.count_nonzero(on)), int(np.count_nonzero(~on))
            scope = f"{x}~{y}"
            rows += _count_rows("pav", scope, mine, n, f"disc=any;res={x},{y};syntenic partners")
            top = sorted(mine, key=lambda k: -(k.end_x - k.start_x))[:3]
            desc = "; ".join(f"{x} {cat.seqids[k.x][k.seq_x]}:{k.start_x}-{k.end_x} "
                             f"({k.n_loci} loci, {k.n_discordant} discordant)" for k in top)
            off = ""
            if n_off:
                gx, gy = cat.g(x), cat.g(y)
                n_pairs = len(np.unique(cat.seq[idx[~on], gx].astype(np.int64) * max(len(cat.lengths[y]), 1)
                                        + cat.seq[idx[~on], gy]))
                off = (f"; {n_off} of the {n + n_off} loci resolved in both lie on {n_pairs} sequence pair(s) that "
                       "are not syntenic partners and are not searched")
            if dropped.get((x, y)):
                off += (f"; {dropped[(x, y)]} reversed run(s) dropped as not contiguous on both genomes (at least as many "
                        "loci placed elsewhere inside their extent as their own)")
            f.add(SECTION, f"{scope}: {len(mine)} inverted block(s)" + (f"; largest: {desc}" if top else "") + off)
    flipped = sorted({f"{k.x}~{k.y}" for k in blocks if k.dominant < 0})
    if flipped:
        f.add(SECTION, f"{', '.join(flipped)}: in some sequence pair the dominant y direction is reversed, so "
                       "the blocks reported there are the collinear stretches (an inversion covering more than "
                       "half of a group's loci, or a reversed sequence, reports the remainder): read them as "
                       "relative orientation")
    f.add(SECTION, "Blocks are searched only between syntenic partner sequences: a sequence and the sequence of "
                   "the other genome that holds most of the loci both genomes place on it, from either side. "
                   "Loci placed on another sequence are not searched, as an inversion is a reversal between "
                   "homologous sequences; a translocated piece that is the smaller share of both its sequences is "
                   "not searched either, so an inversion inside it is missed. A block must be contiguous on both "
                   "genomes: inside its extent on each, fewer loci placed elsewhere than its own (INV_FOREIGN_MAX = "
                   f"{INV_FOREIGN_MAX:g} times them, the merge rule's tolerance for an interruption); on a pair of "
                   "mutually dominant sequences, loci on another syntenic partner of either (an alternate haplotig, a "
                   "duplicated region) do not count. A real inversion whose region's second copy lies on a "
                   "non-partner sequence is still dropped")
    f.add(SECTION, f"A block needs >= {INV_MIN_LOCI} loci of distinct position outside the discordant ones; "
                   f"INV_MIN_LOCI = {INV_MIN_LOCI} is a heuristic, not a tested threshold")
    f.add(SECTION, "Inversions shorter than that, and breakpoints inside unresolved stretches, are not "
                   "blocks; their loci show up as U inv_site (counted in missingness)")
    return rows, f


def _count_rows(layer: str, scope: str, mine: list[InvBlock], n: int, loci: str, *, chain: bool = False
                ) -> list[StatRow]:
    """A pair's inversion rows: n_blocks, n_blocks_1Mb and bp_inverted_x (the x spans end - start + 1, summed). A block
    counts in n_blocks_1Mb when its x span end - start + 1 is >= 1 Mb; a chain block (`chain`), whose extent runs from
    its first to its last anchor's midpoint, when the distance between them, end - start, is."""
    span = [k.end_x - k.start_x + 1 for k in mine]
    long_ = [k.end_x - k.start_x if chain else s for k, s in zip(mine, span, strict=True)]
    return [StatRow(layer, "inversions", scope, "all", stat, v, n=n, loci=loci)
            for stat, v in (("n_blocks", len(mine)), ("n_blocks_1Mb", sum(s >= MB for s in long_)),
                            ("bp_inverted_x", sum(span)))]


# ---------------------------------------------------------------- the library route: chain blocks


def _chain_groups(path: Path) -> dict[tuple[str, str, str, str], list[dict[str, str]]]:
    """blocks.tsv's blocks by (genome_a, genome_b, chrom_a, chrom_b), in order of first appearance."""
    groups: dict[tuple[str, str, str, str], list[dict[str, str]]] = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            groups.setdefault((r["genome_a"], r["genome_b"], r["chrom_a"], r["chrom_b"]), []).append(r)
    return groups


def _main_sign(blocks: list[dict[str, str]]) -> tuple[int, int]:
    """(main orientation, anchors) of one chromosome pair's blocks: the sign holding more anchors (a tie: +1)."""
    w = {1: 0, -1: 0}
    for b in blocks:
        w[int(b["sign"])] += int(b["n_anchor"])
    return (1 if w[1] >= w[-1] else -1), w[1] + w[-1]


def chain_inversions(path: Path, seqids: Mapping[str, Sequence[str]]) -> list[InvBlock]:
    """The library route's inverted blocks (module docstring) from the screen's blocks.tsv at `path`; seqids: label
    -> the genome's sequence names, a name's index its sequence code (LibraryRun.seqids). An AnalysisError names a
    sequence of blocks.tsv that a genome's records lack."""
    codes: dict[str, dict[str, int]] = {}
    out = []
    for (ga, gb, ca, cb), blocks in _chain_groups(path).items():
        dom, total = _main_sign(blocks)
        if total < CHAIN_INV_MIN:
            continue
        for g, c in ((ga, ca), (gb, cb)):
            if g not in codes:
                codes[g] = {name: i for i, name in enumerate(seqids[g])}
            if c not in codes[g]:
                raise AnalysisError(f"{path}: sequence {c} of {g} is no sequence of {g}'s library records: blocks.tsv "
                                    "does not match the libraries (rerun the screen)")
        out += [InvBlock(ga, gb, codes[ga][ca], codes[gb][cb], int(b["start_a"]), int(b["end_a"]), int(b["start_b"]),
                         int(b["end_b"]), None, None, None, None, int(b["n_anchor"]), None, None, dom)
                for b in blocks if int(b["sign"]) != dom]
    return sorted(out, key=lambda k: -k.n_loci)


def chain_inversion_rows(blocks: list[InvBlock], path: Path, labels: Sequence[str],
                         seqids: Mapping[str, Sequence[str]]) -> tuple[list[StatRow], Findings]:
    """Rows (layer library, test "inversions", scope x~y for every pair of the sheet, x before y: n_blocks,
    n_blocks_1Mb, the blocks whose first and last anchor midpoints on x lie >= 1 Mb apart (end - start), and
    bp_inverted_x; n the anchors of the chromosome pairs searched) and findings (section "Inversions":
    each pair's blocks and its three largest, then the rule) of the library route's inverted blocks `blocks`
    (`chain_inversions` of blocks.tsv at `path`)."""
    searched: dict[tuple[str, str], int] = {}
    for (ga, gb, _, _), group in _chain_groups(path).items():
        total = _main_sign(group)[1]
        searched[ga, gb] = searched.get((ga, gb), 0) + (total if total >= CHAIN_INV_MIN else 0)
    rows: list[StatRow] = []
    f = Findings()
    for ix, x in enumerate(labels):
        for y in labels[ix + 1:]:
            mine = [k for k in blocks if k.x == x and k.y == y]
            scope = f"{x}~{y}"
            rows += _count_rows("library", scope, mine, searched.get((x, y), 0), "chain blocks", chain=True)
            top = sorted(mine, key=lambda k: -(k.end_x - k.start_x))[:3]
            desc = "; ".join(f"{x} {seqids[x][k.seq_x]}:{k.start_x}-{k.end_x} ({k.n_loci} anchors)" for k in top)
            f.add(SECTION, f"{scope}: {len(mine)} inverted block(s)" + (f"; largest: {desc}" if top else ""))
    flipped = sorted({f"{k.x}~{k.y}" for k in blocks if k.dominant < 0})
    if flipped:
        f.add(SECTION, f"{', '.join(flipped)}: in some chromosome pair the main orientation is reversed, so the "
                       "blocks reported there are its + blocks: read them as relative orientation")
    f.add(SECTION, "Inverted blocks are the collinear chain blocks of the library screen whose orientation is "
                   "opposite to their chromosome pair's main orientation (the sign holding more anchors; a tie "
                   f"counts as +). A chromosome pair with fewer than CHAIN_INV_MIN = {CHAIN_INV_MIN} anchors is not "
                   f"searched (too few to tell its main orientation); a block needs at least CHAIN_MIN = {CHAIN_MIN} "
                   "anchors, and more than the longest chain among shuffled anchors")
    return rows, f
