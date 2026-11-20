"""Events on one chromosome's genealogy for `ltrtrace sim pop`: LTR-RT insertions, nested up to
depth 2 in simulated insertions; structural losses of elements; SNPs; small indels. Each event sits at a point of
the tree at its position (on the branch above `node`, at `time`, in generations); its carriers are the sample
indices (ts.samples() order) below that point, before any loss.

Rules, in order (W = SLOP_TSD_MULT * TSD_DEFAULT; the backbone is upper-cased, and a base other than A, C, G or T
is never a SNP site, nor one of an insertion's 5 target bases):
  1 primary insertions  msprime binary mutations at lam on a continuous genome, pos = floor(position); one within
                        W + TSD_DEFAULT of an end is dropped (ins_end), one whose target holds a non-ACGT base too
                        (ins_n); donor uniform, strand + or - (1/2 each). One inside a backbone element nests in it
                        by position (host None, depth 0).
  2 losses              backbone elements on the tree at their start, every branch; insertions on the tree at
                        their position (nested ones on their host's), every branch below their own event. A
                        Poisson process at delta per branch, mechanisms by loss_split. solo, del_flank and replace
                        end the element on that lineage (nothing more is drawn below them); del_partial does not.
                        the parameters: solo crossover U{1..min(LTR) - 1}; del_flank SIM_FLANK_DEL bp beyond the
                        start side, the end side or both (1/3 each); del_partial a contiguous SIM_PARTIAL fraction of
                        the element's original internal region at a uniform start, for an element whose original
                        internal region is at least SIM_PARTIAL_MIN (two that overlap on a lineage remove their union);
                        replace length exponential with mean replace_mean, at least MIN_SEG + 1 (the exponential
                        shape is a named heuristic; the mean is measured). A mechanism an element cannot take (solo
                        with an LTR under 2 bp; del_partial under SIM_PARTIAL_MIN) is not drawn there.
  3 nested insertions   for each insertion X after its losses: Poisson(lam * (len(X) - TSD_DEFAULT) * B_X), B_X the
                        branch length of X's subtree below its event without the parts below an ending loss of X;
                        each takes a point uniform by length on that subtree and an offset uniform in
                        [TSD_DEFAULT, len(X)) of X's donor body, so that its target bases lie in X's own sequence
                        (as the end rule keeps primary targets on the backbone); one with a non-ACGT base among its
                        target bases is dropped (ins_n); then its own losses. Depth 2 repeats this once.
  4 existence           an event exists only on sequence present on its lineage at its time (thin_events): from
                        the oldest event down, an insertion is discarded when its host does not exist there or a
                        target base lies in sequence an older kept loss removed (primary: ins_after_loss; nested,
                        any depth and cause: nested_after_loss), a loss when its element does not exist there
                        (loss_after_loss). A Poisson process thinned by a set the older events fix is a Poisson
                        process on what remains, so insertions arrive at lam per present bp per generation, except
                        on two kinds of present sequence that receive none (named simplifications): the replacing
                        sequence of a replace (its random bases join no coordinate system) and the solo LTR a solo
                        leaves of a simulated insertion (B_X stops at X's ending losses).
  5 SNPs                msprime JC69 at mu on a discrete genome; at each site the alleles are relabeled by the
                        permutation of ACGT that sends msprime's ancestral state to the backbone base and keeps the
                        order of the other three (JC69 is symmetric, so the law is unchanged); non-ACGT sites are
                        dropped (snp_n). Not thinned: the builder writes only the sequence present on a lineage.
  6 small indels        msprime binary mutations at indel on a discrete genome, each one event: an insertion of
                        random bases (drawn by the builder) or a deletion, 1/2 each, length geometric with mean 3;
                        a deletion is clipped at the chromosome end. Backbone sequence only; not thinned.
Each msprime call takes the mutation-free tree sequence (keep=False) and a seed drawn from `rng` first.

Coordinates, for the builder. A coordinate system is the backbone (None; [0, L), genome orientation) or a simulated
insertion (its iid): its donor body [0, len) in donor orientation (donor_body), and [len, len + TSD_DEFAULT) its
second TSD copy, counted from the body outward.
  Ins.pos      the copy goes before base pos of its host's system (backbone, or the host's donor body). Its 5
               target bases are [pos - TSD_DEFAULT, pos) of that system; their second copy follows the body in the
               host's order, and the body runs in the host's order forward on strand +, reverse-complemented on -.
  Loss (a, b)  solo: the crossover, an offset into the element's left LTR in its own system, 0. del_flank: bp
               beyond the element's start side, beyond its end side, in its host's order (the backbone's for a
               backbone element). del_partial: the stretch [a, a + b) of the element's original internal region
               (offsets from the end of its left LTR). replace: the length of the replacing sequence, 0.
  removed()    what each kept loss removes, as intervals of these systems; on a lineage the removed bases of a
               system are the union of its kept removals there. The junction rule: an insertion anchored before
               base s of its host's system (pos = s) is removed where bases s - 1 and s are both removed, with its
               body, what is nested in it and its second TSD copy; its host's removal (the host's own anchor removed
               by this rule, or its host's) removes it too, while a host's solo removes only the bases it takes.
Existence at a point q of the tree, read on the carriers of an event (on one tree, exactly the path from q up):
  - a backbone element exists unless older kept losses removed the whole of its span (a named heuristic);
  - a simulated insertion is present below its own event while its anchor survives the older kept removals of its
    host's system (the junction rule) and its host is present, up the host chain; it exists there while, in
    addition, no older kept ending loss of its own lies above q. A host's ending loss does not end what sits in
    sequence the loss leaves (a nested insertion in a host's solo LTR keeps its own losses).
Reading the carriers is the approximation that reading each element on one tree forces: a loss's carriers are fixed
on its element's tree and the builder applies it to them, so across the trees of different positions an event is
discarded when the lineage of any one of its carriers lacks its sequence. That has two effects: insertions whose
sequence is present on some carriers' lineages are discarded with the rest, and a discarded ending loss still
truncated the draw below it, so an element present on some of its carriers' lineages there gets no later losses or
nested insertions on them.
"""
from __future__ import annotations

import bisect
import numbers
from collections import Counter, defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from ltrtrace.base.constants import (
    MIN_SEG,
    SIM_FLANK_DEL,
    SIM_INDEL_P,
    SIM_MAX_CHROMS,
    SIM_MAX_DEPTH,
    SIM_PARTIAL,
    SIM_PARTIAL_MIN,
    SLOP_TSD_MULT,
    TSD_DEFAULT,
)
from ltrtrace.sim.pop.donors import Donor
from ltrtrace.sim.pop.genome import _BASES, _CODE, SimElement
from ltrtrace.sim.pop.model import Rates, _msprime

if TYPE_CHECKING:
    import tskit

KINDS = ("solo", "del_flank", "del_partial", "replace")  # loss mechanisms, in loss_split's order
_ENDING = frozenset({"solo", "del_flank", "replace"})  # these end the element on their lineage
_END = SLOP_TSD_MULT * TSD_DEFAULT + TSD_DEFAULT  # an insertion closer than this to a chromosome end is dropped

# _PERM[a, b, x]: the base code allele x takes at a site whose msprime ancestral state is a and backbone base b
_PERM = np.empty((4, 4, 4), dtype=np.uint8)
for _a in range(4):
    for _b in range(4):
        _PERM[_a, _b, _a] = _b
        _PERM[_a, _b, [x for x in range(4) if x != _a]] = [y for y in range(4) if y != _b]


@dataclass(frozen=True)
class Ins:
    """One LTR-RT insertion event."""

    iid: str  # "i<chrom:02d><n:06d>"
    host: str | None  # None: before backbone base `pos`; else the iid of the simulated insertion it landed in
    pos: int  # backbone position, or offset in the host's donor body
    node: int
    time: float
    donor: int
    strand: str
    depth: int  # 0 primary, 1-2 nested in a simulated insertion
    carriers: frozenset[int]  # sample indices below the event, before any loss


@dataclass(frozen=True)
class Loss:
    """One structural loss of an element on a lineage (the module docstring gives the coordinates of a and b)."""

    element: str  # iid or backbone eid
    node: int
    time: float
    kind: str  # "solo" | "del_flank" | "del_partial" | "replace"
    a: int  # solo: crossover offset, 0; del_flank: bp beyond left end, beyond right end; del_partial: start
    b: int  # offset in the internal region, length; replace: length of the replacing sequence, 0
    carriers: frozenset[int]


@dataclass
class ChromEvents:
    """Every event of one chromosome. snp_base holds ASCII base codes (b"A"[0] ...), one row per SNP site and one
    column per sample; indel_ins is True for an insertion, False for a deletion; dropped counts ins_n, ins_end,
    snp_n, ins_after_loss, nested_after_loss and loss_after_loss. The iids of discarded insertions are not
    reused."""

    ins: list[Ins]
    losses: list[Loss]
    snp_pos: np.ndarray  # (m,) int64
    snp_base: np.ndarray  # (m, n_samples) uint8
    indel_pos: np.ndarray  # int64
    indel_len: np.ndarray  # int64
    indel_ins: np.ndarray  # bool
    indel_carriers: list[frozenset[int]]
    dropped: Counter


def donor_body(d: Donor) -> bytes:
    """What a copy of `d` inserts, before any mutation: the left LTR, the internal region and the left LTR again
    (LTRs identical at insertion, a young copy), upper-cased, in the donor's own orientation."""
    s = d.seq.upper().encode()
    ltr = s[:d.ltr_left]
    return ltr + s[d.ltr_left:len(s) - d.ltr_right] + ltr


class _Tree:
    """One marginal tree, small enough for Python: parent (-1 above a root), time, children and the sample
    indices below each node."""

    __slots__ = ("below", "children", "parent", "roots", "time")

    def __init__(self, tree: tskit.Tree, index: dict[int, int]):
        self.parent: dict[int, int] = {}
        self.time: dict[int, float] = {}
        self.children: dict[int, tuple[int, ...]] = {}
        self.below: dict[int, frozenset[int]] = {}
        for u in tree.postorder().tolist():
            kids = tuple(int(c) for c in tree.children(u))
            self.parent[u] = int(tree.parent(u))
            self.time[u] = float(tree.time(u))
            self.children[u] = kids
            s = {index[u]} if u in index else set()
            for c in kids:
                s |= self.below[c]
            self.below[u] = frozenset(s)
        self.roots = tuple(int(r) for r in tree.roots)


def _trees_at(ts: tskit.TreeSequence, xs: list[float], index: dict[int, int]) -> list[_Tree]:
    """The tree at each position of xs, in one pass over the tree sequence."""
    out: list = [None] * len(xs)
    order = sorted(range(len(xs)), key=xs.__getitem__)
    k = 0
    if not order:
        return out
    for tree in ts.trees():
        right = tree.interval.right
        if xs[order[k]] < right:
            t = _Tree(tree, index)
            while k < len(order) and xs[order[k]] < right:
                out[order[k]] = t
                k += 1
            if k == len(order):
                break
    return out


def _draw_losses(eid: str, ltr_min: int, internal: int, tree: _Tree, origin: tuple[int, float] | None,
                 rates: Rates, cum: np.ndarray, rng: np.random.Generator,
                 out: list[Loss]) -> list[tuple[int, float, float]]:
    """Draw the losses of one element on `tree` below `origin` ((node, time): from that point down; None: every
    branch), appending them to `out`. Returns the segments (node, top, bottom) on which the element still exists:
    its subtree without the parts below an ending loss."""
    alive: list[tuple[int, float, float]] = []
    if origin is None:
        stack = [(c, tree.time[r]) for r in reversed(tree.roots) for c in reversed(tree.children[r])]
    else:
        stack = [origin]
    while stack:
        node, top = stack.pop()
        bottom = tree.time[node]
        n = int(rng.poisson(rates.delta * (top - bottom))) if top > bottom else 0
        end = None
        if n:
            times = np.sort(rng.uniform(bottom, top, n))[::-1]  # oldest first
            kinds = np.searchsorted(cum, rng.random(n), side="right")
            for s, k in zip(times.tolist(), kinds.tolist()):
                kind = KINDS[k]
                if kind == "solo":
                    if ltr_min < 2:
                        continue
                    a, b = int(rng.integers(1, ltr_min)), 0
                elif kind == "del_flank":
                    side = int(rng.integers(3))  # 0: left flank, 1: right flank, 2: both
                    left, right = (int(v) for v in rng.integers(SIM_FLANK_DEL[0], SIM_FLANK_DEL[1] + 1, 2))
                    a, b = (left if side != 1 else 0), (right if side != 0 else 0)
                elif kind == "del_partial":
                    if internal < SIM_PARTIAL_MIN:
                        continue
                    b = int(rng.uniform(*SIM_PARTIAL) * internal + 0.5)
                    a = int(rng.integers(0, internal - b + 1))
                else:
                    a, b = max(MIN_SEG + 1, int(rng.exponential(rates.replace_mean) + 0.5)), 0
                out.append(Loss(eid, node, s, kind, a, b, tree.below[node]))
                if kind in _ENDING:
                    end = s
                    break
        if end is not None:
            alive.append((node, top, end))
        else:
            alive.append((node, top, bottom))
            stack.extend((c, bottom) for c in reversed(tree.children[node]))
    return alive


def _body_len(d: Donor) -> int:
    """len(donor_body(d))."""
    return len(d.seq) + d.ltr_left - d.ltr_right


def _spill(system: str | None, lo: int, hi: int, ins: Mapping[str, Ins], donors: Sequence[Donor], length: int,
           out: list[tuple[str | None, int, int]]) -> None:
    """[lo, hi) of `system` into `out`; a part past a donor body's start or end runs on into its host's system, past
    the end side in the host's order through the body's second TSD copy first; clipped to the backbone."""
    if hi <= lo:
        return
    if system is None:
        lo, hi = max(lo, 0), min(hi, length)
        if hi > lo:
            out.append((None, lo, hi))
        return
    h = ins[system]
    n = _body_len(donors[h.donor])
    if min(hi, n) > max(lo, 0):
        out.append((system, max(lo, 0), min(hi, n)))
    before, after = max(0, -lo), max(0, hi - n)  # bp past the body's start (offset 0) and past its end
    left, right = (before, after) if h.strand == "+" else (after, before)  # in the host's order
    if left:
        _spill(h.host, h.pos - left, h.pos, ins, donors, length, out)
    if right:
        out.append((system, n, n + min(right, TSD_DEFAULT)))
        _spill(h.host, h.pos, h.pos + right - TSD_DEFAULT, ins, donors, length, out)


def removed(loss: Loss, elements: Mapping[str, SimElement], ins: Mapping[str, Ins], donors: Sequence[Donor],
            length: int) -> list[tuple[str | None, int, int]]:
    """What `loss` removes: half-open intervals (system, start, end) of the coordinate systems of the module
    docstring (None: the backbone; an iid: that insertion's body, then its second TSD copy). The one definition of
    each kind, for an element whose own system spans [start, end) with LTRs [ltr_left0, ltr_left1) and
    [ltr_right0, ltr_right1) (an insertion: [0, len), LTRs [0, ltr) and [len - ltr, len)):
      solo         [ltr_left0 + a, ltr_right0 + a)
      del_flank    a backbone element: [start - a, end + b), clipped to the backbone (its
                   TSDs are backbone bases). An insertion at pos of its host's system: its body [0, len) and, when
                   b > 0, its second TSD copy; then [pos - a, pos) of the host (the start-side flank, from its
                   target bases outward) and, after the TSD copy, [pos, pos + b - TSD_DEFAULT) (the end-side
                   flank), in the host's order. A flank running past a donor body's start or end goes on into that
                   body's host, through the body's second TSD copy on its end side.
      del_partial  [ltr_left1 + a, ltr_left1 + a + b): two on one lineage that overlap remove their union
      replace      [start, end): a backbone element's span, or an insertion's body (its TSDs stay; no backbone
                   sequence goes)
    Every interval is measured in its system's original coordinates. Bases an older removal already took count in
    it, as for the union of two del_partials: a flank of an insertion that sits in its host's solo LTR runs on
    through the host's removed stretch and stops where its length ends. It stops short of the LTR's far part when
    the flank (a, or b - TSD_DEFAULT, in the host) is no longer than the host bases between its start and that
    part, as always when the host's removed stretch (ltr + internal bp) is at least the flank's length; a longer
    flank takes bases of the far part too.
    What a lineage loses is the union of its kept removals per system, read with the junction rule: an insertion
    anchored before base s of its host's system (pos = s) is removed where bases s - 1 and s are both removed, and
    with its host. The host counts as removed when its own anchor is removed by the same rule, or its own host is,
    up the chain. A host's own loss removes a nested insertion only through the bases it takes: a solo leaves one
    in the LTR it keeps, while a del_flank or a replace takes the whole body and every insertion anchored in it.
    So an insertion W at the same position as X goes with X's del_flank when it takes both sides (a > 0 and
    b > TSD_DEFAULT: [pos - a, pos) and [pos, ...)), and stays under a one-sided one or X's replace.
    Nested offsets are drawn in [TSD_DEFAULT, len), so a nested insertion's anchor and target lie in its host's
    body. `elements` maps backbone eids to SimElements, `ins` iids to Ins (the host chain of an insertion
    included), `donors` are draw_events' donors and `length` the backbone's."""
    if loss.element in elements:
        e = elements[loss.element]
        system, start, end, l0, l1, r0 = None, e.start, e.end, e.ltr_left[0], e.ltr_left[1], e.ltr_right[0]
    else:
        x = ins[loss.element]
        d = donors[x.donor]
        system, start, end = x.iid, 0, _body_len(d)
        l0, l1, r0 = 0, d.ltr_left, end - d.ltr_left
    a, b = loss.a, loss.b
    if loss.kind == "solo":
        return [(system, l0 + a, r0 + a)]
    if loss.kind == "del_partial":
        return [(system, l1 + a, l1 + a + b)]
    if loss.kind == "replace":
        return [(system, start, end)]
    if loss.kind != "del_flank":
        raise ValueError(f"unknown loss kind {loss.kind!r}")
    if system is None:
        return [(None, max(0, start - a), min(length, end + b))]
    out: list[tuple[str | None, int, int]] = [(system, 0, end + min(b, TSD_DEFAULT))]
    _spill(x.host, x.pos - a, x.pos, ins, donors, length, out)
    _spill(x.host, x.pos, x.pos + b - TSD_DEFAULT, ins, donors, length, out)
    return out


class _Removals:
    """The kept removals of one coordinate system, sorted by start: (start, end, time, carriers)."""

    __slots__ = ("items", "maxlen", "starts")

    def __init__(self):
        self.starts: list[int] = []
        self.items: list[tuple[int, int, float, frozenset[int]]] = []
        self.maxlen = 0

    def add(self, lo: int, hi: int, time: float, carriers: frozenset[int]) -> None:
        k = bisect.bisect_right(self.starts, lo)
        self.starts.insert(k, lo)
        self.items.insert(k, (lo, hi, time, carriers))
        self.maxlen = max(self.maxlen, hi - lo)

    def hits(self, lo: int, hi: int, time: float, carriers: frozenset[int]) -> list[tuple[int, int, frozenset[int]]]:
        """Removals overlapping [lo, hi), older than `time`, on the lineage of one of `carriers`."""
        k0, k1 = bisect.bisect_left(self.starts, lo - self.maxlen), bisect.bisect_left(self.starts, hi)
        return [(a, b, c) for a, b, t, c in self.items[k0:k1] if b > lo and t > time and c & carriers]


def _covers(spans: list[tuple[int, int]], lo: int, hi: int) -> bool:
    """Does the union of `spans` cover [lo, hi)?"""
    for a, b in sorted(spans):
        if a > lo:
            return False
        lo = max(lo, b)
        if lo >= hi:
            return True
    return lo >= hi


def thin_events(ins: list[Ins], losses: list[Loss], elements: list[SimElement], donors: list[Donor],
                length: int) -> tuple[list[Ins], list[Loss], Counter]:
    """The events that exist (module docstring, rule 4), in their input order, and the counts of those discarded:
    ins_after_loss (primary insertions), nested_after_loss (nested ones) and loss_after_loss. Events are decided
    from the oldest down, each against the older kept ones; a discarded event removes nothing."""
    els = {e.eid: e for e in elements}
    by_iid = {i.iid: i for i in ins}
    rem: dict[str | None, _Removals] = defaultdict(_Removals)
    ending: dict[str, list[tuple[float, frozenset[int]]]] = defaultdict(list)
    kept_ins: set[str] = set()
    dropped = Counter({"ins_after_loss": 0, "nested_after_loss": 0, "loss_after_loss": 0})

    def present(x: Ins, carriers: frozenset[int], time: float) -> bool:
        """Is x's sequence there: below its event, its anchor kept (the junction rule) and its host's too?"""
        if x.iid not in kept_ins or not time < x.time or not carriers <= x.carriers:
            return False
        s = x.pos  # anchored before base s: gone where bases s - 1 and s are both removed
        hits = rem[x.host].hits(s - 1, s + 1, time, carriers)
        if hits and any(any(a <= s - 1 < b and c in cs for a, b, cs in hits)
                        and any(a <= s < b and c in cs for a, b, cs in hits) for c in carriers):
            return False
        return x.host is None or present(by_iid[x.host], carriers, time)

    def exists(eid: str, carriers: frozenset[int], time: float) -> bool:
        if eid in els:  # a backbone element: gone on a lineage whose older removals cover its whole span
            e = els[eid]
            hits = rem[None].hits(e.start, e.end, time, carriers)
            return not hits or not any(_covers([(a, b) for a, b, c in hits if s in c], e.start, e.end)
                                       for s in carriers)
        return present(by_iid[eid], carriers, time) and not any(t > time and c & carriers for t, c in ending[eid])

    order = sorted([(-i.time, k) for k, i in enumerate(ins)] + [(-y.time, len(ins) + k) for k, y in enumerate(losses)])
    keep = [False] * (len(ins) + len(losses))
    for _, k in order:
        if k < len(ins):
            x = ins[k]
            ok = x.host is None or exists(x.host, x.carriers, x.time)
            lo, hi = max(0, x.pos - TSD_DEFAULT), x.pos
            if ok and hi > lo and rem[x.host].hits(lo, hi, x.time, x.carriers):
                ok = False
            if ok:
                keep[k] = True
                kept_ins.add(x.iid)
            else:
                dropped["ins_after_loss" if x.host is None else "nested_after_loss"] += 1
        else:
            y = losses[k - len(ins)]
            if not exists(y.element, y.carriers, y.time):
                dropped["loss_after_loss"] += 1
                continue
            keep[k] = True
            if y.kind in _ENDING:
                ending[y.element].append((y.time, y.carriers))
            for system, lo, hi in removed(y, els, by_iid, donors, length):
                if hi > lo:
                    rem[system].add(lo, hi, y.time, y.carriers)
    return ([i for k, i in enumerate(ins) if keep[k]], [y for k, y in enumerate(losses) if keep[len(ins) + k]],
            dropped)


def _check(ts: tskit.TreeSequence, backbone: bytes, elements: list[SimElement], donors: list[Donor],
           rates: Rates, chrom: int) -> None:
    if ts.sequence_length != len(backbone):
        raise ValueError(f"backbone length {len(backbone)} != the tree sequence's length {ts.sequence_length:g}")
    if isinstance(chrom, bool) or not isinstance(chrom, numbers.Integral) or not 0 <= chrom <= SIM_MAX_CHROMS:
        raise ValueError(f"chrom must be an integer 0-{SIM_MAX_CHROMS}: {chrom!r}")
    seen: set[str] = set()
    for e in elements:
        (a, b), (c, d) = e.ltr_left, e.ltr_right
        if e.eid in seen:
            raise ValueError(f"duplicate element id {e.eid}")
        seen.add(e.eid)
        if not 0 <= e.start <= a < b <= c < d <= e.end <= len(backbone):
            raise ValueError(f"element {e.eid}: LTRs {e.ltr_left} and {e.ltr_right} do not lie in order inside "
                             f"[{e.start}, {e.end}) on the {len(backbone)}-bp backbone")
    for i, d in enumerate(donors):
        if not (0 < d.ltr_left and 0 < d.ltr_right and d.ltr_left + d.ltr_right < len(d.seq)):
            raise ValueError(f"donor {i}: LTRs of {d.ltr_left} and {d.ltr_right} bp leave no internal region in "
                             f"{len(d.seq)} bp")
    if rates.lam > 0 and not donors:
        raise ValueError("no donor elements to insert (lam > 0)")


def draw_events(ts: tskit.TreeSequence, backbone: bytes, elements: list[SimElement], donors: list[Donor],
                rates: Rates, rng: np.random.Generator, *, chrom: int = 0) -> ChromEvents:
    """Every event of one chromosome on its genealogy `ts` (mutation-free; its length the backbone's), by the
    rules of the module docstring. `elements` are the backbone's (origin "backbone"; 0-based half-open
    coordinates), `donors` the templates of insertions (Ins.donor indexes them), `chrom` the number in the
    insertion ids. Deterministic for a given `rng` state. Bad input raises ValueError."""
    _check(ts, backbone, elements, donors, rates, chrom)
    msp = _msprime()
    seeds = [int(s) for s in rng.integers(1, 2**32 - 1, size=3)]  # insertions, SNPs, indels
    bb = bytes(backbone).upper()
    length = len(bb)
    bb_code = _CODE[np.frombuffer(bb, dtype=np.uint8)]
    bad = np.zeros(length + 1, dtype=np.int64)
    np.cumsum(bb_code > 3, out=bad[1:])  # bad[j] - bad[i]: non-ACGT bases in [i, j)
    index = {int(u): i for i, u in enumerate(ts.samples())}
    dropped = Counter({"ins_n": 0, "ins_end": 0, "snp_n": 0})
    split = np.asarray(rates.loss_split, dtype=float)
    cum = np.cumsum(split / split.sum())
    cum[-1] = 1.0
    ins: list[Ins] = []
    losses: list[Loss] = []
    count = 0

    def iid() -> str:
        nonlocal count
        count += 1
        return f"i{chrom:02d}{count:06d}"

    # 1. primary insertions
    mts = msp.sim_mutations(ts, rate=rates.lam, model=msp.BinaryMutationModel(), discrete_genome=False,
                            random_seed=seeds[0], keep=False)
    xs = mts.sites_position[mts.mutations_site].tolist()  # one mutation per site on a continuous genome
    keep = []
    for k, x in enumerate(xs):
        pos = int(x)
        if pos < _END or pos > length - _END:
            dropped["ins_end"] += 1
        elif bad[pos] - bad[pos - TSD_DEFAULT]:
            dropped["ins_n"] += 1
        else:
            keep.append(k)
    els = sorted(elements, key=lambda e: (e.start, -e.end, e.eid))
    trees = _trees_at(ts, [xs[k] for k in keep] + [e.start for e in els], index)
    ins_tree = trees[:len(keep)]
    node, when = mts.mutations_node.tolist(), mts.mutations_time.tolist()
    for k, tree in zip(keep, ins_tree):
        d = int(rng.integers(len(donors)))
        strand = "+" if rng.random() < 0.5 else "-"
        ins.append(Ins(iid(), None, int(xs[k]), node[k], when[k], d, strand, 0, tree.below[node[k]]))

    # 2. losses of backbone elements
    for e, tree in zip(els, trees[len(keep):]):
        _draw_losses(e.eid, min(e.ltr_left[1] - e.ltr_left[0], e.ltr_right[1] - e.ltr_right[0]),
                     e.ltr_right[0] - e.ltr_left[1], tree, None, rates, cum, rng, losses)

    # 2-3. losses of each insertion, then the insertions nested in it, theirs, and so on to SIM_MAX_DEPTH
    bodies: dict[int, bytes] = {}

    def body(d: int) -> bytes:
        if d not in bodies:
            bodies[d] = donor_body(donors[d])
        return bodies[d]

    def settle(x: Ins, tree: _Tree) -> None:
        ltr = donors[x.donor].ltr_left
        size = len(body(x.donor))
        alive = _draw_losses(x.iid, ltr, size - 2 * ltr, tree, (x.node, x.time), rates, cum, rng, losses)
        if x.depth >= SIM_MAX_DEPTH:
            return
        seg = np.array([top - bottom for _, top, bottom in alive])
        total = float(seg.sum())
        n = int(rng.poisson(rates.lam * (size - TSD_DEFAULT) * total)) if total > 0 and size > TSD_DEFAULT else 0
        ends = np.cumsum(seg)
        kids = []
        for _ in range(n):
            u = float(rng.uniform(0, total))
            j = min(int(np.searchsorted(ends, u, side="right")), len(alive) - 1)
            v, top, _ = alive[j]
            t = float(top - (u - (ends[j] - seg[j])))
            off = int(rng.integers(TSD_DEFAULT, size))
            if body(x.donor)[off - TSD_DEFAULT:off].translate(None, b"ACGT"):
                dropped["ins_n"] += 1
                continue
            d = int(rng.integers(len(donors)))
            strand = "+" if rng.random() < 0.5 else "-"
            kids.append(Ins(iid(), x.iid, off, v, t, d, strand, x.depth + 1, tree.below[v]))
        ins.extend(kids)
        for y in kids:
            settle(y, tree)

    for x, tree in zip(list(ins), ins_tree):
        settle(x, tree)

    # 4. existence
    ins, losses, gone = thin_events(ins, losses, elements, donors, length)
    dropped.update(gone)

    # 5. SNPs
    sts = msp.sim_mutations(ts, rate=rates.mu, model=msp.JC69(), discrete_genome=True, random_seed=seeds[1],
                            keep=False)
    snp_pos = sts.sites_position.astype(np.int64)
    snp_base = np.empty((0, ts.num_samples), dtype=np.uint8)
    if len(snp_pos):
        sites = sts.tables.sites
        if not np.all(np.diff(sites.ancestral_state_offset) == 1):
            raise RuntimeError("JC69 ancestral states are not single bases")
        anc = _CODE[sites.ancestral_state.view(np.uint8)]
        ref = bb_code[snp_pos]
        ok = ref < 4
        dropped["snp_n"] = int((~ok).sum())
        geno = sts.genotype_matrix(alleles=("A", "C", "G", "T"), isolated_as_missing=False)[ok]
        snp_pos = snp_pos[ok]
        snp_base = _BASES[_PERM[anc[ok, None], ref[ok, None], geno]]

    # 6. small indels
    its = msp.sim_mutations(ts, rate=rates.indel, model=msp.BinaryMutationModel(), discrete_genome=True,
                            random_seed=seeds[2], keep=False)
    indel_pos = its.sites_position[its.mutations_site].astype(np.int64)
    carriers: list[frozenset[int]] = []
    for tree in its.trees():
        if not tree.num_sites:
            continue
        seen: dict[int, frozenset[int]] = {}
        for site in tree.sites():
            for m in site.mutations:
                if m.node not in seen:
                    seen[m.node] = frozenset(index[int(s)] for s in tree.samples(m.node))
                carriers.append(seen[m.node])
    m = len(indel_pos)
    indel_ins = rng.random(m) < 0.5
    indel_len = rng.geometric(SIM_INDEL_P, m).astype(np.int64)
    indel_len[~indel_ins] = np.minimum(indel_len[~indel_ins], length - indel_pos[~indel_ins])
    return ChromEvents(ins, losses, snp_pos, snp_base, indel_pos, indel_len, indel_ins, carriers, dropped)
