"""One-pass genome builder of `synltr sim pop`: each sample's sequence of one chromosome and the record
of every element as it stands there (Copy), from that chromosome's events (popevents.ChromEvents); then assembly gaps.

Coordinates and existence are popevents' (its module docstring, "Coordinates, for the builder" and "Existence at a
point q"); removed() is the one definition of what each loss removes. Every event handed to build_chrom must exist
on each of its carriers (popevents rule 4, existence): the builder checks that once per chromosome with
thin_events and raises ValueError otherwise. A sample holds an event when it is one of its carriers.

Per sample, on each coordinate system (the backbone; each simulated insertion's donor body and second TSD copy):
  removed      the union of the removed() intervals of the kept losses the sample carries
  insertion    present while its host is present and its anchor survives (the junction rule: not both of bases
               pos - 1 and pos removed); placed before base pos of its host's system, the youngest first at one
               position, as its body (forward on +, reverse-complemented on -) then its second TSD copy, both in
               the host's order
  replace      `a` random bases (drawn once per loss), in place of the element's span or body, after whatever is
               anchored at its start and before whatever is anchored at its end; gone with its insertion, and for a
               backbone element where the sample's other removals cover the whole span (they took its place)
  SNP          the sample's base (ChromEvents.snp_base) at every backbone base still present
  small indel  backbone sequence only. Inserted bases are drawn once per event and go before base pos, after the
               insertions anchored there; gone where bases pos - 1 and pos are both removed. A deletion takes
               [pos, pos + len) of the backbone, never a payload, unless a removal that overlaps it lies above it
               (the removal's carriers include the indel's: thin_events' reading of "above" on one tree): then it
               does not happen, on any carrier. Indel event times are not in ChromEvents, so an overlapping removal
               with the same carriers counts as above (a named approximation; a few bp per event).
The sequence is emitted once, as pieces in order (backbone slices with the sample's SNP bases, payloads, replacing
and indel bases), and joined once; element coordinates are cumulative piece lengths. No per-event copy is made.

Payloads. An insertion's body is its donor body (donor_body) with copy_noise substitutions, one draw per insertion,
then its left LTR copied over its right. Its second TSD copy is its target bases at the insertion time. Body and TSD
copy then evolve down the event's subtree, one JC69 step per branch segment (genome.substitute with
p = 3/4 (1 - exp(-4/3 mu t)), the chance that a base differs after t generations at mu), so orthologous copies share
the substitutions of their shared branches; each sample keeps its differences from the inserted content. A nested
insertion's target bases are its host's content at the nested event's point. A primary insertion's are the
backbone's with, at each SNP site, the derived base when every carrier holds it and its samples include the carriers
(the mutation lies above the insertion; when they equal the carriers it lies on the insertion's own branch, above it
with probability (top - time) / (top - bottom)); SNP times are not in ChromEvents, so this reading assumes one
mutation per site (a named approximation). The first TSD copy is the host's sequence in place. Replacing sequences
and inserted indel bases take no substitutions (a named simplification): every carrier holds the bases drawn for
the event.

Copy records (one per backbone element and per insertion, in every sample):
  fate      absent (not inherited); else, in this order: replaced (its own replace, the replacing bases still there;
            deleted when another removal took them), deleted (its own del_flank; an insertion's anchor or host
            removed; a backbone element's whole span or an insertion's whole body removed), solo (its own solo),
            partial (a removal took own bases), intact. Small indel deletions change the kept own bases, never the fate.
  image     start, end (half-open): genome coordinates of the first and last kept base of a backbone element's span
            or an insertion's body, what is nested between included; a replaced copy's are its replacing bases'.
            LTR images, ltr_own and ltr_kept are in genome order (left: the lower coordinates); a solo's LTRs are its
            whole image. No image: (-1, -1).
  site      for an absent copy whose empty-site window is intact: the genome coordinate of its first target base
            (in genome order; on the backbone a base an indel deleted maps to the next kept one, as lift_pos). The
            window is the target bases and w = SLOP_TSD_MULT * TSD_DEFAULT bases each side in the host's system
            (clipped to a host's body); only removals of its bases break it, never SNPs, small indels or insertions
            at other sites. An absent host leaves no window. Else -1.
  k2p       intact and partial copies: an insertion's two LTRs compared position by position where both bases are
            present (K2P); a backbone element's p-distance from its source K, plus the sample's derived SNPs in its
            kept own LTR bases per mean own LTR length, back to K2P with the source's transition share (read from
            its two backbone LTRs when they have equal length and differ; JC otherwise). None without a source K.
            Source K is read as K at the root (the backbone is the root genome), so each sample adds the derived
            SNPs of its whole root-to-sample path: about 2 mu T_root above source K (a named approximation).
  tsd       an insertion's target bases at the insertion time in genome order; a backbone element's source TSD.
  strand, family, superfamily: the element's genome orientation, and its donor's or source's family and
            superfamily, for the annotation.
  kept_runs: the own bases the sample still has, as sorted half-open runs in the element's own system
            (backbone coordinates; an insertion's donor-body offsets), for its left LTR, its whole own sequence and
            its right LTR (LTRs in genome order, as ltr_kept): the own sequence is a backbone element's span
            without the backbone elements nested wholly in it, or an insertion's donor body. Truth reads one
            sample's copy against another's own bases from them ("Y holds >= COV_PRESENT of e's own ids": a
            loss two samples share takes the same runs from both).
  events    the kinds of its own losses the sample carries, oldest first; then "deleted" when it is gone other than
            by its own del_flank, or "clipped" when another element's removal took own bases of a copy still there.
Lift: the backbone tiled into pieces, each piece's start in backbone and in genome coordinates, then (L, genome
length): a kept piece maps base by base, a deleted piece (the same genome start as the next entry) to the next kept
base. lift_pos reads it.

Gaps (add_gaps): the rule per genome, at per_mb 100-bp N runs per Mb: a third start or end within SIM_GAP_NEAR of
the end of an intact or partial copy or of an absent copy's site, the rest at uniform starts; recorded in gaps. A
near run that falls wholly off a chromosome is not written; add_gaps returns how many did, per genome.

Input checks (ValueError, before anything is drawn): lengths and shapes; backbone elements on the backbone; ids,
hosts, depths, donors, strands and carriers; an insertion's base pos and target bases in its host's system
([TSD_DEFAULT, L) of the backbone, [TSD_DEFAULT, len) of a host's donor body); every event on its node's branch
(node time <= time < parent time) on the tree at its position; then existence (thin_events).
"""
from __future__ import annotations

import bisect
import math
import numbers
from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from synltr.base.constants import SIM_GAP_LEN, SIM_GAP_NEAR, SIM_GAP_NEAR_DIV, SLOP_TSD_MULT, TSD_DEFAULT
from synltr.seq.divergence import k2p_shares
from synltr.seq.intervals import union
from synltr.sim.pop.donors import Donor
from synltr.sim.pop.events import KINDS, ChromEvents, Ins, _body_len, _trees_at, donor_body, removed, thin_events
from synltr.sim.pop.genome import SimElement, revcomp_bytes, substitute
from synltr.sim.pop.model import Rates

if TYPE_CHECKING:
    import tskit

W = SLOP_TSD_MULT * TSD_DEFAULT  # empty-site window flank, each side of the target bases
_ENDING = frozenset({"solo", "del_flank", "replace"})
_ACGT = np.frombuffer(b"ACGT", dtype=np.uint8)
_IS_ACGT = np.zeros(256, dtype=bool)
_IS_ACGT[_ACGT] = True
_PURINE = np.zeros(256, dtype=bool)
_PURINE[[ord("A"), ord("G")]] = True
_NONE = (-1, -1)


@dataclass
class Copy:
    """One element (backbone eid or iid) as it stands in one built genome (module docstring)."""

    eid: str
    chrom: int
    inherited: bool  # the sample lies below the element's origin (a backbone element: always)
    fate: str  # absent, replaced, deleted, solo, partial or intact
    start: int  # image in the genome (0-based half-open, chromosome coordinates), -1 if none
    end: int
    ltr_left: tuple[int, int]
    ltr_right: tuple[int, int]
    ltr_own: tuple[int, int]
    ltr_kept: tuple[int, int]
    site: int  # the target site of an absent copy whose empty-site window is intact, else -1
    target: tuple[int, str | None, int]  # (chrom, host, pos) of an insertion's target; backbone (chrom, None, -1)
    host: str | None
    tsd: str
    k2p: float | None
    events: list[str]
    strand: str = "+"  # genome orientation
    family: str | None = None
    superfamily: str | None = None
    kept_runs: tuple[tuple[tuple[int, int], ...], ...] = ((), (), ())  # own bases kept: left LTR, own, right LTR


Lift = tuple[np.ndarray, np.ndarray]  # the backbone's pieces: start in backbone and in genome coordinates (int64,
#                                       ascending), then (L, genome length); a deleted piece shares the next start


@dataclass
class PopGenome:
    """One sample's genome: chromosome i is seqids[i], seqs[i] and lift[i]; copies of every chromosome by id."""

    label: str
    group: str
    seqids: list[str]
    seqs: list[bytes]
    copies: dict[str, Copy]
    gaps: list[tuple[int, int, int]]  # (chrom, start, end)
    lift: list[Lift]  # per chromosome


# ------------------------------------------------------------------------------------------------ intervals

class _Ivs:
    """A sorted union of half-open intervals (touching ones merged)."""

    __slots__ = ("his", "los")

    def __init__(self, ivs=()):
        self.los: list[int] = []
        self.his: list[int] = []
        for lo, hi in union(ivs):
            self.los.append(lo)
            self.his.append(hi)

    def has(self, p: int) -> bool:
        k = bisect.bisect_right(self.los, p) - 1
        return k >= 0 and p < self.his[k]

    def junction(self, p: int) -> bool:
        """Are bases p - 1 and p both in it (the junction rule)?"""
        k = bisect.bisect_right(self.los, p - 1) - 1
        return k >= 0 and p < self.his[k]

    def hits(self, lo: int, hi: int) -> bool:
        k = bisect.bisect_right(self.los, lo) - 1
        if k >= 0 and self.his[k] > lo:
            return True
        return k + 1 < len(self.los) and self.los[k + 1] < hi

    def covers(self, lo: int, hi: int) -> bool:
        k = bisect.bisect_right(self.los, lo) - 1
        return k >= 0 and self.his[k] >= hi

    def kept(self, lo: int, hi: int) -> list[tuple[int, int]]:
        """[lo, hi) without the union, as sorted runs."""
        out = []
        k = max(bisect.bisect_right(self.los, lo) - 1, 0)
        cur = lo
        while k < len(self.los) and self.los[k] < hi:
            if self.his[k] > cur:
                if self.los[k] > cur:
                    out.append((cur, self.los[k]))
                cur = self.his[k]
            k += 1
        if cur < hi:
            out.append((cur, hi))
        return out


_EMPTY = _Ivs()


class _Rows:
    """Removal intervals of one system, each with the index of its loss, sorted by start: overlap queries."""

    __slots__ = ("maxlen", "rows", "starts")

    def __init__(self, rows):
        self.rows = sorted(rows)
        self.starts = [r[0] for r in self.rows]
        self.maxlen = max((hi - lo for lo, hi, _ in self.rows), default=0)

    def over(self, lo: int, hi: int) -> list[tuple[int, int, int]]:
        i0 = bisect.bisect_left(self.starts, lo - self.maxlen)
        i1 = bisect.bisect_left(self.starts, hi)
        return [r for r in self.rows[i0:i1] if r[1] > lo]


def _sweep(runs: list[tuple[int, int]], anchors: list[tuple]):
    """Kept runs and anchored items of one system in its order, an anchor at p before base p: yields (lo, hi, None)
    for runs and (None, None, anchor) for anchors (sorted by their first fields)."""
    i, start = 0, 0
    for a in anchors:
        p = a[0]
        while i < len(runs):
            lo, hi = runs[i]
            lo = max(lo, start)
            if lo >= p:
                break
            if hi <= p:
                yield lo, hi, None
                i += 1
            else:
                yield lo, p, None
                start = p
                break
        yield None, None, a
    while i < len(runs):
        lo, hi = runs[i]
        lo = max(lo, start)
        if hi > lo:
            yield lo, hi, None
        i += 1


def _index(marks: list[tuple[int, int, int, int]]) -> tuple[list[int], list[int], list[int], list[int]]:
    """(lo, hi, genome start, orientation) runs of one system, sorted, as four lists."""
    marks = sorted(marks)
    return ([m[0] for m in marks], [m[1] for m in marks], [m[2] for m in marks], [m[3] for m in marks])


def _image(ix, a: int, b: int) -> tuple[int, int, int]:
    """(first, last + 1, count) of the genome positions of the kept bases [a, b) of a system; (-1, -1, 0) if none."""
    if ix is None:
        return -1, -1, 0
    los, his, offs, sgs = ix
    k = max(bisect.bisect_right(los, a) - 1, 0)
    g0, g1, n = -1, -1, 0
    while k < len(los) and los[k] < b:
        lo, hi = los[k], his[k]
        u, v = max(a, lo), min(b, hi)
        if v > u:
            if sgs[k] > 0:
                x0, x1 = offs[k] + u - lo, offs[k] + v - lo
            else:
                x0, x1 = offs[k] + hi - v, offs[k] + hi - u
            g0 = x0 if g0 < 0 else min(g0, x0)
            g1 = max(g1, x1)
            n += v - u
        k += 1
    return g0, g1, n


# ------------------------------------------------------------------------------------------------ distances

def _k2p(p_ts: float, q_tv: float) -> float | None:
    """Kimura 2-parameter distance from the transition and transversion shares (seq.divergence.k2p_shares); None when
    saturated, never -0.0."""
    d = k2p_shares(p_ts, q_tv)[0]
    return None if math.isnan(d) else d + 0.0


def _p_of_k2p(k: float, f: float) -> float:
    """The p-distance whose K2P distance, at transition share f, is k (bisection; K2P grows with p)."""
    hi = min(1 / (1 + f), 1 / (2 * (1 - f)) if f < 1 else math.inf)  # where K2P diverges
    lo = 0.0
    for _ in range(100):
        mid = (lo + hi) / 2
        d = _k2p(mid * f, mid * (1 - f))
        if d is None or d > k:
            hi = mid
        else:
            lo = mid
    return lo


def _counts(a: np.ndarray, b: np.ndarray) -> tuple[int, int, int]:
    """(compared, transitions, transversions) between two equal-length base arrays, over A/C/G/T pairs."""
    ok = _IS_ACGT[a] & _IS_ACGT[b]
    a, b = a[ok], b[ok]
    diff = a != b
    same_class = _PURINE[a] == _PURINE[b]
    return int(ok.sum()), int((diff & same_class).sum()), int((diff & ~same_class).sum())


def _evolve(c: bytes, dt: float, mu: float, rng: np.random.Generator) -> bytes:
    """c after dt generations of JC69 substitution at mu per base."""
    if mu <= 0 or dt <= 0:
        return c
    return substitute(c, 0.75 * -math.expm1(-4.0 / 3.0 * mu * dt), rng)


# ------------------------------------------------------------------------------------------------ one chromosome

def _own_intervals(elements: Sequence[SimElement]) -> dict[str, list[tuple[int, int]]]:
    """Each backbone element's span without the spans of the elements nested wholly in it (in
    (start, -end, eid) order, a later element ending inside an earlier one is nested in it)."""
    els = sorted(elements, key=lambda e: (e.start, -e.end, e.eid))
    out = {}
    for i, e in enumerate(els):
        inner = []
        for f in els[i + 1:]:
            if f.start >= e.end:
                break
            if f.end <= e.end:
                inner.append((f.start, f.end))
        out[e.eid] = _Ivs(inner).kept(e.start, e.end)
    return out


def _clip(ivs: list[tuple[int, int]], lo: int, hi: int) -> list[tuple[int, int]]:
    return [(max(a, lo), min(b, hi)) for a, b in ivs if min(b, hi) > max(a, lo)]


class _Chrom:
    """One chromosome: checks, every random draw (in a fixed order, before any sample is built), and sample(s), the
    built sample s. Samples may be built in any order with identical results."""

    def __init__(self, chrom: int, backbone: bytes, elements: list[SimElement], ts: tskit.TreeSequence,
                 ev: ChromEvents, donors: list[Donor], rates: Rates, rng: np.random.Generator):
        self.chrom = chrom
        self.bb = bytes(backbone).upper()
        self.L = len(self.bb)
        self.n = ts.num_samples
        self.elements = list(elements)
        self.els = {e.eid: e for e in self.elements}
        self.ins = list(ev.ins)
        self.by_iid = {x.iid: x for x in self.ins}
        self.losses = list(ev.losses)
        self.donors = list(donors)
        self.ev = ev
        self._check(ts)
        kept_ins, kept_losses, _ = thin_events(self.ins, self.losses, self.elements, self.donors, self.L)
        if len(kept_ins) != len(self.ins) or len(kept_losses) != len(self.losses):
            raise ValueError(
                f"chromosome {chrom}: {len(self.ins) - len(kept_ins)} insertion(s) and "
                f"{len(self.losses) - len(kept_losses)} loss(es) do not exist on their lineage (popevents rule 4, "
                "existence: thin_events would discard them); build from draw_events' output, or thin hand-built "
                "events first")
        self.rem = [[(s, lo, hi) for s, lo, hi in removed(y, self.els, self.by_iid, self.donors, self.L) if hi > lo]
                    for y in self.losses]
        self.order = sorted(self.ins, key=lambda x: x.depth)  # hosts before what is nested in them
        self.kids: dict[str, list[Ins]] = defaultdict(list)
        for x in self.ins:
            if x.host is not None:
                self.kids[x.host].append(x)
        self.orient: dict[str | None, int] = {None: 1}
        for x in self.order:
            self.orient[x.iid] = self.orient[x.host] * (1 if x.strand == "+" else -1)
        bodies = {i: donor_body(d) for i, d in enumerate(self.donors)}
        self.body = {x.iid: bodies[x.donor] for x in self.ins}
        self.own_losses: dict[str, list[int]] = defaultdict(list)  # oldest first
        for k in sorted(range(len(self.losses)), key=lambda k: (-self.losses[k].time, k)):
            self.own_losses[self.losses[k].element].append(k)
        self.own_iv = _own_intervals(self.elements)
        self.ltr_iv = {e.eid: (_clip(self.own_iv[e.eid], *e.ltr_left), _clip(self.own_iv[e.eid], *e.ltr_right))
                       for e in self.elements}
        self.ts_share = {e.eid: self._ts_share(e) for e in self.elements}
        self.p0 = {e.eid: _p_of_k2p(float(e.k), self.ts_share[e.eid]) for e in self.elements if e.k is not None}
        self.snp_pos = np.asarray(ev.snp_pos, dtype=np.int64)
        self.snp_base = np.asarray(ev.snp_base, dtype=np.uint8)
        self.snp_ref = np.frombuffer(self.bb, dtype=np.uint8)[self.snp_pos]
        # draws: indel bases, replacing sequences, then payloads, each in event order
        self.indel_bases = {k: _ACGT[rng.integers(0, 4, int(ev.indel_len[k]))].tobytes()
                            for k in range(len(ev.indel_pos)) if ev.indel_ins[k]}
        self.rep_seq = {k: _ACGT[rng.integers(0, 4, y.a)].tobytes() for k, y in enumerate(self.losses)
                        if y.kind == "replace"}
        self._payloads(ts, rates, rng)
        rows = _Rows((lo, hi, k) for k, r in enumerate(self.rem) for s, lo, hi in r if s is None)
        self.indel_skip = {k for k in range(len(ev.indel_pos)) if not ev.indel_ins[k] and any(
            self.losses[j].carriers >= ev.indel_carriers[k]
            for _, _, j in rows.over(int(ev.indel_pos[k]), int(ev.indel_pos[k] + ev.indel_len[k])))}
        self.held: set[int] = set()  # indels some built sample holds (build_chrom's `held`)

    # ---------------------------------------------------------------------------------------- setup

    def _check(self, ts) -> None:
        ev, n = self.ev, self.n
        if ts.sequence_length != self.L:
            raise ValueError(f"backbone length {self.L} != the tree sequence's length {ts.sequence_length:g}")
        if np.shape(ev.snp_base) != (len(ev.snp_pos), n):
            raise ValueError(f"snp_base has shape {np.shape(ev.snp_base)}, expected ({len(ev.snp_pos)}, {n})")
        if len(ev.snp_pos) and (np.any(np.diff(ev.snp_pos) <= 0) or ev.snp_pos[0] < 0 or ev.snp_pos[-1] >= self.L):
            raise ValueError("SNP positions must be distinct, ascending and on the backbone")
        m = len(ev.indel_pos)
        if not len(ev.indel_len) == len(ev.indel_ins) == len(ev.indel_carriers) == m:
            raise ValueError("indel arrays of different lengths")
        if m and (np.any(ev.indel_pos < 0) or np.any(ev.indel_pos > self.L) or np.any(ev.indel_len < 1)):
            raise ValueError("an indel lies off the backbone or has no length")
        if len(self.els) != len(self.elements):
            raise ValueError("duplicate backbone element ids")
        if len(self.by_iid) != len(self.ins) or set(self.by_iid) & set(self.els):
            raise ValueError("duplicate insertion ids, or an insertion id equal to a backbone element's")
        for e in self.elements:
            if not 0 <= e.start < e.end <= self.L:
                raise ValueError(f"backbone element {e.eid} [{e.start}, {e.end}) is not on the backbone [0, {self.L})")
        everyone = frozenset(range(n))
        for x in self.ins:
            if x.host is not None and x.host not in self.by_iid:
                raise ValueError(f"insertion {x.iid}: unknown host {x.host}")
            if x.depth != (0 if x.host is None else self.by_iid[x.host].depth + 1):
                raise ValueError(f"insertion {x.iid}: depth {x.depth} does not follow its host chain")
            if not 0 <= x.donor < len(self.donors) or x.strand not in "+-" or len(x.strand) != 1:
                raise ValueError(f"insertion {x.iid}: bad donor {x.donor} or strand {x.strand!r}")
            if not x.carriers <= everyone:
                raise ValueError(f"insertion {x.iid}: carriers outside the {n} samples")
        for x in self.ins:  # base pos and the target bases [pos - TSD_DEFAULT, pos) lie in the host's system
            hi = self.L if x.host is None else _body_len(self.donors[self.by_iid[x.host].donor])
            if not TSD_DEFAULT <= x.pos < hi:
                where = "the backbone" if x.host is None else f"its host {x.host}'s donor body"
                raise ValueError(f"insertion {x.iid}: pos {x.pos} outside [{TSD_DEFAULT}, {hi}) of {where} (its "
                                 "target bases and base pos must lie in it)")
        for y in self.losses:
            if y.element not in self.els and y.element not in self.by_iid:
                raise ValueError(f"loss of an unknown element {y.element}")
            if y.kind not in KINDS:
                raise ValueError(f"loss of {y.element}: unknown kind {y.kind!r}")
            if not y.carriers <= everyone:
                raise ValueError(f"loss of {y.element}: carriers outside the {n} samples")
        self._check_times(ts)

    def _check_times(self, ts) -> None:
        """Every event lies on its node's branch (node time <= time < parent time) on the tree at its position: an
        insertion's (a nested one's: its primary's), a loss's (a backbone element's start, else its insertion's).
        Keeps the primaries' trees (self.trees) for _payloads."""
        index = {int(u): i for i, u in enumerate(ts.samples())}
        root: dict[str, str] = {}
        for x in sorted(self.ins, key=lambda x: x.depth):
            root[x.iid] = x.iid if x.host is None else root[x.host]
        prim = [x.iid for x in self.ins if x.host is None]
        lost = sorted({y.element for y in self.losses if y.element in self.els})
        xs = [float(self.by_iid[i].pos) for i in prim] + [float(self.els[e].start) for e in lost]
        trees = dict(zip(prim + lost, _trees_at(ts, xs, index)))

        def on_branch(what: str, key: str, node: int, time: float) -> None:
            tree = trees[key]
            if node not in tree.time:
                raise ValueError(f"{what}: node {node} is not in the tree at its position")
            bottom, top = tree.time[node], tree.time.get(tree.parent.get(node, -1))
            if top is None or not bottom <= time < top:
                why = "the node is a root" if top is None else f"[{bottom:g}, {top:g})"
                raise ValueError(f"{what}: time {time:g} is not on the branch above node {node} ({why}, on the tree "
                                 "at its position)")

        for x in self.ins:
            on_branch(f"insertion {x.iid}", root[x.iid], x.node, x.time)
        for y in self.losses:
            on_branch(f"loss of {y.element}", y.element if y.element in self.els else root[y.element], y.node, y.time)
        self.trees = {i: trees[i] for i in prim}

    def _ts_share(self, e: SimElement) -> float:
        """Transition share of e's source divergence: read position by position from its two backbone LTRs when
        they have equal length and differ; else JC's 1/3."""
        (a, b), (c, d) = e.ltr_left, e.ltr_right
        if b - a != d - c:
            return 1 / 3
        x = np.frombuffer(self.bb, dtype=np.uint8)
        _, ts_n, tv_n = _counts(x[a:b], x[c:d])
        return ts_n / (ts_n + tv_n) if ts_n + tv_n else 1 / 3

    def _payloads(self, ts, rates: Rates, rng: np.random.Generator) -> None:
        """Each insertion's inserted content (c0: body then second TSD copy) and, per carrier, its differences
        after evolving down the event's subtree. Hosts first, so a nested insertion reads its target bases from
        its host's content at its point."""
        index = {int(u): i for i, u in enumerate(ts.samples())}
        root: dict[str, str] = {}
        for x in self.order:
            root[x.iid] = x.iid if x.host is None else root[x.host]
        self.c0: dict[str, bytes] = {}
        self.diff: dict[str, dict[int, tuple[np.ndarray, np.ndarray]]] = {}
        self.tsd0: dict[str, bytes] = {}
        target: dict[str, bytes] = {}
        for x in self.order:
            d = self.donors[x.donor]
            body = self.body[x.iid]
            if rates.copy_noise > 0:
                body = substitute(body, rates.copy_noise, rng)
            ltr = d.ltr_left
            body = body[:len(body) - ltr] + body[:ltr]
            tree = self.trees[root[x.iid]]  # _check_times: x.node lies in it, x.time on its branch
            if x.host is None:
                tsd = self._target_at(x, tree, rng)
            elif x.iid in target:
                tsd = target.pop(x.iid)
            else:
                raise ValueError(f"insertion {x.iid} does not lie on the subtree of its host {x.host}")
            self.tsd0[x.iid] = tsd
            self.c0[x.iid] = body + tsd
            self.diff[x.iid] = self._walk(x, body + tsd, tree, index, rates.mu, rng, target)
            if set(self.diff[x.iid]) != x.carriers:
                raise ValueError(f"insertion {x.iid}: its carriers are not the samples below its point")

    def _target_at(self, x: Ins, tree, rng: np.random.Generator) -> bytes:
        """A primary insertion's target bases at its time (module docstring)."""
        lo, hi = x.pos - TSD_DEFAULT, x.pos
        t = bytearray(self.bb[lo:hi])
        j0, j1 = (int(j) for j in np.searchsorted(self.snp_pos, [lo, hi]))
        if j1 == j0:
            return bytes(t)
        car = np.array(sorted(x.carriers), dtype=np.int64)
        for j in range(j0, j1):
            q, col, b0 = int(self.snp_pos[j]), self.snp_base[j], int(self.snp_ref[j])
            held = np.unique(col[car])
            if len(held) != 1 or int(held[0]) == b0:
                continue  # a carrier keeps the backbone base: the mutation lies below the insertion
            if frozenset(np.flatnonzero(col != b0).tolist()) == x.carriers:  # on the insertion's own branch
                up = tree.parent.get(x.node, -1)
                top, bottom = tree.time.get(up), tree.time[x.node]
                if top is None or top <= bottom or rng.random() * (top - bottom) >= top - x.time:
                    continue
            t[q - lo] = int(held[0])
        return bytes(t)

    def _walk(self, x: Ins, c0: bytes, tree, index: dict[int, int], mu: float, rng: np.random.Generator,
              target: dict[str, bytes]) -> dict[int, tuple[np.ndarray, np.ndarray]]:
        """Evolve x's content down its subtree (depth first, children in tree order); record the target bases of
        each insertion nested in x at its point, and each sample's differences from c0."""
        points: dict[int, list[Ins]] = defaultdict(list)
        for y in self.kids[x.iid]:
            points[y.node].append(y)
        ref = np.frombuffer(c0, dtype=np.uint8)
        out = {}
        stack = [(x.node, x.time, c0)]
        while stack:
            v, top, c = stack.pop()
            bottom = tree.time[v]
            for y in sorted(points.pop(v, ()), key=lambda y: -y.time):
                c = _evolve(c, top - y.time, mu, rng)
                top = y.time
                target[y.iid] = c[y.pos - TSD_DEFAULT:y.pos]
            c = _evolve(c, top - bottom, mu, rng)
            if v in index:
                a = np.frombuffer(c, dtype=np.uint8)
                where = np.flatnonzero(a != ref)
                out[index[v]] = (where, a[where].copy())
            for ch in reversed(tree.children[v]):
                stack.append((ch, bottom, c))
        return out

    # ---------------------------------------------------------------------------------------- one sample

    def sample(self, s: int) -> tuple[bytes, dict[str, Copy], Lift]:
        """Sample s's sequence, its copy of every element and its lift map."""
        ev, losses = self.ev, self.losses
        mine = [k for k, y in enumerate(losses) if s in y.carriers]
        raw: dict[str | None, list[tuple[int, int, int]]] = defaultdict(list)
        for k in mine:
            for sys_, lo, hi in self.rem[k]:
                raw[sys_].append((lo, hi, k))
        rs = {sys_: _Ivs((lo, hi) for lo, hi, _ in v) for sys_, v in raw.items()}
        rows = {sys_: _Rows(v) for sys_, v in raw.items()}
        ending = {losses[k].element: k for k in mine if losses[k].kind in _ENDING}
        present: dict[str, bool] = {}
        for x in self.order:
            present[x.iid] = (s in x.carriers and (x.host is None or present[x.host])
                              and not rs.get(x.host, _EMPTY).junction(x.pos))
        rep = set()  # elements whose replacing sequence the sample holds
        for el, k in ending.items():
            if losses[k].kind != "replace":
                continue
            if el in self.els:
                e = self.els[el]
                others = _Ivs((lo, hi) for lo, hi, j in rows[None].over(e.start, e.end) if j != k)
                if not others.covers(e.start, e.end):
                    rep.add(el)
            elif present[el]:
                rep.add(el)
        rb = rs.get(None, _EMPTY)
        mine_dels = [k for k in range(len(ev.indel_pos))
                     if not ev.indel_ins[k] and s in ev.indel_carriers[k] and k not in self.indel_skip]
        dels = [(int(ev.indel_pos[k]), int(ev.indel_pos[k] + ev.indel_len[k])) for k in mine_dels]
        self.held.update(mine_dels)
        gone = _Ivs(list(zip(rb.los, rb.his)) + dels)  # backbone bases the sample lacks
        arr = np.frombuffer(self.bb, dtype=np.uint8).copy()
        if len(self.snp_pos):
            arr[self.snp_pos] = self.snp_base[:, s]
        view = memoryview(arr)
        anchors: list[tuple] = []
        for el in rep:
            if el in self.els:
                anchors.append((self.els[el].end, 0, 0.0, ending[el], "rep", el))
        for k, x in enumerate(self.ins):
            if x.host is None and present[x.iid]:
                anchors.append((x.pos, 1, x.time, k, "ins", x))
        for k in range(len(ev.indel_pos)):
            p = int(ev.indel_pos[k])
            if ev.indel_ins[k] and s in ev.indel_carriers[k] and not rb.junction(p):
                anchors.append((p, 2, 0.0, k, "indel", k))
                self.held.add(k)
        anchors.sort(key=lambda a: a[:4])
        ctx = _Sample(s, rs, present, rep, ending, self)
        pieces: list = []
        marks: list[tuple] = []
        rmarks: list[tuple[str, int, int]] = []
        bmarks: list[tuple[int, int, int]] = []
        off = 0
        for lo, hi, a in _sweep(gone.kept(0, self.L), anchors):
            if a is None:
                pieces.append(view[lo:hi])
                bmarks.append((lo, hi, off))
                off += hi - lo
            elif a[4] == "ins":
                b, m, r = self._payload(a[5], ctx)
                pieces.append(b)
                marks += [(sy, u, v, o + off, sg) for sy, u, v, o, sg in m]
                rmarks += [(e, o + off, ln) for e, o, ln in r]
                off += len(b)
            else:
                b = self.rep_seq[a[3]] if a[4] == "rep" else self.indel_bases[a[3]]
                if a[4] == "rep":
                    rmarks.append((a[5], off, len(b)))
                pieces.append(b)
                off += len(b)
        seq = b"".join(pieces)
        del pieces, view, arr
        groups: dict[str | None, list[tuple[int, int, int, int]]] = defaultdict(list)
        for sy, u, v, o, sg in marks:
            groups[sy].append((u, v, o, sg))
        ix: dict[str | None, tuple] = {sy: _index(g) for sy, g in groups.items()}
        ix[None] = ([m[0] for m in bmarks], [m[1] for m in bmarks], [m[2] for m in bmarks], [1] * len(bmarks))
        lift = _lift(bmarks, self.L, len(seq))
        reps = {e: (o, o + ln) for e, o, ln in rmarks}
        copies: dict[str, Copy] = {}
        for e in self.elements:
            copies[e.eid] = self._bb_copy(e, ctx, rows, ix, reps, gone)
        for x in self.ins:
            copies[x.iid] = self._ins_copy(x, ctx, rows, ix, reps, lift)
        return seq, copies, lift

    def _payload(self, x: Ins, ctx: _Sample) -> tuple[bytes, list[tuple], list[tuple]]:
        """x's payload in its host's order (body oriented by strand, then its second TSD copy), the runs of each
        system in it (system, lo, hi, offset, orientation) and its replacing sequences (element, offset, length)."""
        n = len(self.body[x.iid])
        c = ctx.content(x.iid)
        rx = ctx.rs.get(x.iid, _EMPTY)
        anchors: list[tuple] = [(y.pos, 1, y.time, k, "ins", y) for k, y in enumerate(self.kids[x.iid])
                                if ctx.present[y.iid]]
        if x.iid in ctx.rep:
            anchors.append((n, 0, 0.0, 0, "rep", x.iid))
        anchors.sort(key=lambda a: a[:4])
        pieces, marks, rmarks, off = [], [], [], 0
        for lo, hi, a in _sweep(rx.kept(0, n), anchors):
            if a is None:
                pieces.append(c[lo:hi])
                marks.append((x.iid, lo, hi, off, 1))
                off += hi - lo
            elif a[4] == "ins":
                b, m, r = self._payload(a[5], ctx)
                pieces.append(b)
                marks += [(sy, u, v, o + off, sg) for sy, u, v, o, sg in m]
                rmarks += [(e, o + off, ln) for e, o, ln in r]
                off += len(b)
            else:
                b = self.rep_seq[ctx.ending[x.iid]]
                pieces.append(b)
                rmarks.append((x.iid, off, len(b)))
                off += len(b)
        body = b"".join(pieces)
        if x.strand == "-":
            nb = len(body)
            body = revcomp_bytes(body)
            marks = [(sy, u, v, nb - o - (v - u), -sg) for sy, u, v, o, sg in marks]
            rmarks = [(e, nb - o - ln, ln) for e, o, ln in rmarks]
        tail, off = [], len(body)
        for lo, hi in rx.kept(n, n + TSD_DEFAULT):
            tail.append(c[lo:hi])
            marks.append((x.iid, lo, hi, off, 1))
            off += hi - lo
        return body + b"".join(tail), marks, rmarks

    def _events(self, eid: str, s: int, rows: dict, own: list[tuple[int, int]], sys_: str | None,
                fate: str, ending: str | None) -> list[str]:
        out = [self.losses[k].kind for k in self.own_losses.get(eid, ()) if s in self.losses[k].carriers]
        if fate == "deleted":
            if ending != "del_flank":
                out.append("deleted")
        elif fate != "absent" and sys_ in rows and any(
                self.losses[j].element != eid for lo, hi in own for _, _, j in rows[sys_].over(lo, hi)):
            out.append("clipped")
        return out

    def _bb_copy(self, e: SimElement, ctx: _Sample, rows: dict, ix: dict, reps: dict, gone: _Ivs) -> Copy:
        s, rb = ctx.s, ctx.rs.get(None, _EMPTY)
        own_iv = self.own_iv[e.eid]
        k = ctx.ending.get(e.eid)
        kind = self.losses[k].kind if k is not None else None
        if kind == "replace":
            fate = "replaced" if e.eid in ctx.rep else "deleted"
        elif kind == "del_flank" or rb.covers(e.start, e.end):
            fate = "deleted"
        elif kind == "solo":
            fate = "solo"
        elif any(rb.hits(lo, hi) for lo, hi in own_iv):
            fate = "partial"
        else:
            fate = "intact"
        bix = ix[None]
        ltr_own = tuple(sum(hi - lo for lo, hi in iv) for iv in self.ltr_iv[e.eid])
        ltr_kept = tuple(sum(_image(bix, lo, hi)[2] for lo, hi in iv) for iv in self.ltr_iv[e.eid])
        start, end = _NONE
        left = right = _NONE
        k2p = None
        if fate in ("intact", "partial", "solo"):
            start, end, _ = _image(bix, e.start, e.end)
            if fate == "solo":
                left = right = (start, end)
            else:
                left, right = (_image(bix, *e.ltr_left)[:2], _image(bix, *e.ltr_right)[:2])
                k2p = self._bb_k2p(e, s, ltr_own, gone)
        elif fate == "replaced":
            start, end = reps.get(e.eid, _NONE)
        runs = tuple(tuple(r for lo, hi in iv for r in gone.kept(lo, hi))
                     for iv in (self.ltr_iv[e.eid][0], own_iv, self.ltr_iv[e.eid][1]))
        return Copy(e.eid, self.chrom, True, fate, start, end, left, right, ltr_own, ltr_kept, -1,
                    (self.chrom, None, -1), None, e.tsd or "", k2p,
                    self._events(e.eid, s, rows, own_iv, None, fate, kind), e.strand, e.family, e.superfamily, runs)

    def _bb_k2p(self, e: SimElement, s: int, ltr_own: tuple[int, int], gone: _Ivs) -> float | None:
        if e.k is None or not sum(ltr_own):
            return None
        f = self.ts_share[e.eid]
        derived = 0
        for iv in self.ltr_iv[e.eid]:
            for lo, hi in iv:
                j0, j1 = (int(j) for j in np.searchsorted(self.snp_pos, [lo, hi]))
                for j in range(j0, j1):
                    if self.snp_base[j, s] != self.snp_ref[j] and not gone.has(int(self.snp_pos[j])):
                        derived += 1
        p = self.p0[e.eid] + derived / (sum(ltr_own) / 2)
        return _k2p(p * f, p * (1 - f))

    def _ins_copy(self, x: Ins, ctx: _Sample, rows: dict, ix: dict, reps: dict, lift: Lift) -> Copy:
        s = ctx.s
        d = self.donors[x.donor]
        n, ltr = len(self.body[x.iid]), d.ltr_left
        o = self.orient[x.iid]
        tsd = self.tsd0[x.iid] if self.orient[x.host] > 0 else revcomp_bytes(self.tsd0[x.iid])
        flip = (lambda a, b: (a, b)) if o > 0 else (lambda a, b: (b, a))
        common = {"target": (self.chrom, x.host, x.pos), "host": x.host, "tsd": tsd.decode(),
                  "strand": "+" if o > 0 else "-", "family": d.family, "superfamily": d.superfamily}
        if s not in x.carriers:
            return Copy(x.iid, self.chrom, False, "absent", -1, -1, _NONE, _NONE, (ltr, ltr), (0, 0),
                        self._site(x, ctx, ix, lift), k2p=None, events=[], **common)
        rx = ctx.rs.get(x.iid, _EMPTY)
        k = ctx.ending.get(x.iid)
        kind = self.losses[k].kind if k is not None else None
        if not ctx.present[x.iid]:
            fate = "deleted"
        elif kind == "replace":
            fate = "replaced"
        elif kind == "del_flank" or rx.covers(0, n):
            fate = "deleted"
        elif kind == "solo":
            fate = "solo"
        elif rx.hits(0, n):
            fate = "partial"
        else:
            fate = "intact"
        xi = ix.get(x.iid) if ctx.present[x.iid] else None
        sl, sr = (0, ltr), (n - ltr, n)
        ltr_kept = flip(_image(xi, *sl)[2], _image(xi, *sr)[2])
        start, end = _NONE
        left = right = _NONE
        k2p = None
        if fate in ("intact", "partial", "solo"):
            start, end, _ = _image(xi, 0, n)
            if fate == "solo":
                left = right = (start, end)
            else:
                left, right = flip(_image(xi, *sl)[:2], _image(xi, *sr)[:2])
                k2p = self._ins_k2p(x, ctx, rx, n, ltr)
        elif fate == "replaced":
            start, end = reps.get(x.iid, _NONE)
        runs = ((), (), ())
        if ctx.present[x.iid]:
            lo_ltr, hi_ltr = flip(tuple(rx.kept(*sl)), tuple(rx.kept(*sr)))
            runs = (lo_ltr, tuple(rx.kept(0, n)), hi_ltr)
        return Copy(x.iid, self.chrom, True, fate, start, end, left, right, (ltr, ltr), ltr_kept, -1,
                    k2p=k2p, events=self._events(x.iid, s, rows, [(0, n)], x.iid, fate, kind), kept_runs=runs,
                    **common)

    def _ins_k2p(self, x: Ins, ctx: _Sample, rx: _Ivs, n: int, ltr: int) -> float | None:
        a = np.frombuffer(ctx.content(x.iid), dtype=np.uint8)
        keep = np.zeros(n, dtype=bool)
        for lo, hi in rx.kept(0, n):
            keep[lo:hi] = True
        both = keep[:ltr] & keep[n - ltr:n]
        m, ts_n, tv_n = _counts(a[:ltr][both], a[n - ltr:n][both])
        return _k2p(ts_n / m, tv_n / m) if m else None

    def _site(self, x: Ins, ctx: _Sample, ix: dict, lift: Lift) -> int:
        """The site of x in a sample that did not inherit it (module docstring), or -1."""
        h = x.host
        if h is not None and not ctx.present[h]:
            return -1
        lo, hi = x.pos - TSD_DEFAULT - W, x.pos + W
        lo, hi = max(lo, 0), min(hi, self.L if h is None else len(self.body[h]))
        if ctx.rs.get(h, _EMPTY).hits(lo, hi):
            return -1
        if h is None:
            return lift_pos(lift, x.pos - TSD_DEFAULT)
        return _image(ix.get(h), x.pos - TSD_DEFAULT, x.pos)[0]


class _Sample:
    """The state of one sample while it is built: its removals per system, present insertions, replacing
    sequences it holds, ending losses, and its payload contents (made on demand)."""

    __slots__ = ("chrom", "ending", "memo", "present", "rep", "rs", "s")

    def __init__(self, s: int, rs: dict, present: dict, rep: set, ending: dict, chrom: _Chrom):
        self.s, self.rs, self.present, self.rep, self.ending, self.chrom = s, rs, present, rep, ending, chrom
        self.memo: dict[str, bytes] = {}

    def content(self, iid: str) -> bytes:
        """The insertion's content (body then second TSD copy, its own order) as it stands in this sample."""
        if iid not in self.memo:
            a = np.frombuffer(self.chrom.c0[iid], dtype=np.uint8).copy()
            where, vals = self.chrom.diff[iid][self.s]
            a[where] = vals
            self.memo[iid] = a.tobytes()
        return self.memo[iid]


def _lift(bmarks: list[tuple[int, int, int]], length: int, glen: int) -> Lift:
    bb, g = [], []
    cur = 0
    for lo, hi, off in bmarks:
        if lo > cur:
            bb.append(cur)
            g.append(off)  # a deleted piece: its bases map to the next kept one
        bb.append(lo)
        g.append(off)
        cur = hi
    if cur < length:
        bb.append(cur)
        g.append(glen)
    bb.append(length)
    g.append(glen)
    return np.array(bb, dtype=np.int64), np.array(g, dtype=np.int64)


# ------------------------------------------------------------------------------------------------ public

def build_chrom(chrom: int, backbone: bytes, elements: list[SimElement], ts: tskit.TreeSequence, ev: ChromEvents,
                donors: list[Donor], rates: Rates,
                rng: np.random.Generator, *, held: set[int] | None = None
                ) -> list[tuple[bytes, dict[str, Copy], Lift]]:
    """Every sample's (sequence, copies, lift) of one chromosome, in ts.samples() order, by the rules of the module
    docstring. `elements`, `donors` and `ev` are those of draw_events on `ts` and `backbone`; rates.mu and
    rates.copy_noise are read. Deterministic for a given `rng` state, whatever order the samples are built in.
    Events that do not exist on their lineage, or malformed input, raise ValueError. held: when given, the small
    indels (indexes into ev's indel arrays) that at least one sample holds are added to it: a deletion the rule of
    the module docstring lets happen (it then happens on every carrier), inserted bases a carrier keeps (not gone
    with both junction bases); an indel outside it is in no genome."""
    c = _Chrom(chrom, backbone, elements, ts, ev, donors, rates, rng)
    out = [c.sample(s) for s in range(c.n)]
    if held is not None:
        held.update(c.held)
    return out


def lift_pos(lift: Lift, pos: int) -> int:
    """Genome coordinate of backbone position pos (0 <= pos <= L): a deleted position maps to the next kept base,
    and L to the genome's length."""
    bb, g = lift
    if not bb[0] <= pos <= bb[-1]:
        raise ValueError(f"position {pos} outside the backbone [0, {int(bb[-1])}]")
    k = int(np.searchsorted(bb, pos, side="right")) - 1
    if k == len(bb) - 1:
        return int(g[-1])
    if g[k] == g[k + 1]:
        return int(g[k])
    return int(g[k] + pos - bb[k])


def add_gaps(genomes: list[PopGenome], per_mb: float, rngs: list[np.random.Generator]) -> list[int]:
    """Assembly gaps (module docstring), each genome from its own stream: round(per_mb x its length / 1 Mb) runs of
    SIM_GAP_LEN Ns, a third (rounded; SIM_GAP_NEAR_DIV) near the end of an intact or partial copy or an absent copy's
    site (the start SIM_GAP_NEAR or fewer bp past it, or the end SIM_GAP_NEAR or fewer bp before it, 1/2 each), the
    rest at uniform starts that fit the chromosome. A copy's chrom indexes the genome's seqs. Writes N runs into seqs and appends
    (chrom, start, end) to gaps; a near run is clipped at a chromosome end, and one that falls wholly off it is not
    written. Returns, per genome, how many near runs fell off (so a genome holds that many fewer runs)."""
    if len(rngs) != len(genomes):
        raise ValueError(f"{len(genomes)} genomes but {len(rngs)} random streams")
    if isinstance(per_mb, bool) or not isinstance(per_mb, numbers.Real) or not per_mb >= 0:
        raise ValueError(f"per_mb must be a number >= 0: {per_mb!r}")
    dropped = []
    for g, rng in zip(genomes, rngs):
        sizes = [len(x) for x in g.seqs]
        n = int(per_mb * sum(sizes) / 1e6 + 0.5)
        dropped.append(0)
        if not n:
            continue
        near = int(n / SIM_GAP_NEAR_DIV + 0.5)
        ends = sorted({(c.chrom, x) for c in g.copies.values() if c.fate in ("intact", "partial") and c.start >= 0
                       for x in (c.start, c.end)} | {(c.chrom, c.site) for c in g.copies.values()
                                                       if c.fate == "absent" and c.site >= 0})
        room = np.array([max(1, z - SIM_GAP_LEN + 1) for z in sizes], dtype=np.int64)
        cum = np.cumsum(room)
        bufs = [bytearray(x) for x in g.seqs]
        for i in range(n):
            if i < near and ends:
                ch, t = ends[int(rng.integers(len(ends)))]
                u = int(rng.integers(0, SIM_GAP_NEAR + 1))
                st = t + u if rng.random() < 0.5 else t - u - SIM_GAP_LEN
            else:
                x = int(rng.integers(0, int(cum[-1])))
                ch = int(np.searchsorted(cum, x, side="right"))
                st = x - int(cum[ch] - room[ch])
            a, b = max(st, 0), min(st + SIM_GAP_LEN, sizes[ch])
            if b > a:
                bufs[ch][a:b] = b"N" * (b - a)
                g.gaps.append((ch, a, b))
            else:
                dropped[-1] += 1
        g.seqs = [bytes(x) for x in bufs]
    return dropped
