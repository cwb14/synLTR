"""State classifier: what genome B holds at each annotated element of genome A in an aligned window.

The window's elements are called outermost-first (a host before what is nested in it). For each, the rules
run in order and the first that applies gives the call (thresholds from base.constants):

    0 inheritance   host E -> E + host_empty; host D/O/G/U -> the host's state + inherit; host P/S/T -> go on
    1 window        window status not ok -> U (reason = status); no alignment -> U no_alignment
    2 flanks        < MIN_SEG aligned bp on either side of the element -> U flank_unaligned
    3 assembly gap  an N in B within REFINE_PAD of the site's image -> G
    4 absent        < MIN_SEG of the element's own bases aligned: an inverted hit over the element whose B
                    interval meets the site (between the flanks' B images) is U inv_site; else the gap
                    around the element is refined, balanced in its slide range and read: a B insertion
                    >= MIN_SEG is O; both gap ends within -w..w of the element ends (an overhang past the
                    element filled by other deleted elements counts as within: adjacent) is E (O small_ins if
                    B adds > w bp); anything else is D (filler if B adds > w bp)
    5 present       both LTRs and the element >= COV_PRESENT aligned -> P (nestB: a B insertion inside);
                    one LTR's worth aligned and next to no internal region -> S (solo_tsd: B's bases on both
                    sides of the solo LTR match A's TSD); else T
    5a anchored     a P, S or T call needs the element's own aligned bases joined to its flank on one side:
                    an alignment run (segments joined across A and B gaps < MIN_SEG) crosses the element's
                    start or end holding >= MIN_SEG aligned A bases on each side of it, in the flank and
                    among the element's own bases (a run that barely enters the element does not join it
                    to its flank); or (solo_tsd) the element's B image (the images of its LTR ends) lies
                    between rule 3's lo and hi and B's bases beside it match an A TSD copy on both sides;
                    else U unanchored (its nested elements inherit it, rule 0). Bases that align to a
                    paralogous copy elsewhere in the window are not the element at its site.

Notation: s, t = the element's window interval (0-based half-open); L = its TSD length (TSD_DEFAULT when
not annotated); w = SLOP_TSD_MULT * L. The spans of an element's annotated children are left out of every
aligned-bp count and part length of it. ident = local A-B identity over [s - IDENT_FLANK, t + IDENT_FLANK).
Tier: placement (window margin >= MARGIN_HIGH and ident >= IDENT_MIN) plus the state's own evidence (P:
always; S: solo_tsd; E: tsd_match; D: an exactly refined junction; T, O: none) -> high if both, medium if
one, low if neither; G and U get ".".

`b_seq` is always in A's orientation (the caller reverse-complements a strand -1 window); `win.strand`
only sets the B coordinates (genome, 1-based inclusive) and the inv flag. inv is orientation, not biology by
itself: B's window lies on A's reverse strand, through an inversion or only a contig assembled the other way
round (unoriented scaffolds); inversions come from the net's inverted fills (`chain.net_blocks`).
"""
from __future__ import annotations

import math
from bisect import bisect_left, bisect_right
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from synltr.base.constants import (
    COV_PRESENT,
    IDENT_FLANK,
    IDENT_MIN,
    MARGIN_HIGH,
    MIN_SEG,
    REFINE_PAD,
    SLOP_TSD_MULT,
    SOLO_LEN,
    TSD_DEFAULT,
)
from synltr.calls.records import CALL_COLUMNS
from synltr.calls.states import PRESENT
from synltr.io.elements import Element
from synltr.routes.assembly.lastz import AlnMap
from synltr.routes.assembly.window import Window
from synltr.seq.intervals import union
from synltr.seq.junction import Gap, hamming, refine_gap, tsd_allowance

_EV_KEYS = ("reason", "cov_left", "cov_right", "cov_int", "cov_elem", "dl", "dr", "b_ins", "tsd_a", "tsd_b",
            "tsd_mm", "ident", "anchor_left", "anchor_right", "margin")


@dataclass
class Call:
    """An element's state in B.

    tier: high | medium | low, or "." for G and U. B coordinates: B's genome, 1-based inclusive; None
    without a site. ev: the evidence columns of CALL_COLUMNS (reason included), None where absent.
    """

    eid: str
    state: str
    tier: str
    flags: list[str]
    b_seqid: str | None
    b_start: int | None
    b_end: int | None
    b_strand: str | None
    ev: dict[str, float | int | str | None]


def classify_window(win: Window, aln: AlnMap | None, a_seq: str, b_seq: str,
                    elements: dict[str, Element], children: dict[str, list[str]]) -> list[Call]:
    """Calls for the window's elements (win.eids and their annotated descendants), parents first.

    `elements` maps eid -> Element (it may hold more than the window); `children` maps eid -> the eids of
    its direct children. a_seq, b_seq: the window sequences, B in A's orientation.
    """
    ctx = _Ctx(win, aln, a_seq, b_seq, elements, children)
    calls: dict[str, Call] = {}
    for e in ctx.members:
        x = _Site(ctx, e)
        calls[e.eid] = (_inherit(x, calls.get(e.parent)) or _window_rule(x) or _flank_rule(x)
                        or _gap_rule(x) or (_absent(x) if x.own_aligned < MIN_SEG else _present(x)))
    return [calls[e.eid] for e in ctx.members]


def call_row(call: Call, e: Element, win: Window) -> dict:
    """The CALL_COLUMNS row of `call` (element `e`, window `win`): raw Python values, None where absent,
    flags comma-joined (None if there are none). A coordinates are e's (1-based inclusive)."""
    if call.eid != e.eid:
        raise ValueError(f"call eid {call.eid} does not match element eid {e.eid}")
    row = dict(call.ev, eid=e.eid, a_seqid=e.seqid, a_start=e.start, a_end=e.end, parent=e.parent,
               level=e.level, state=call.state, tier=call.tier, flags=",".join(call.flags) or None,
               b_seqid=call.b_seqid, b_start=call.b_start, b_end=call.b_end, b_strand=call.b_strand,
               wid=int(win.wid))
    return {k: row[k] for k in CALL_COLUMNS}


# ---------------------------------------------------------------- rules (first match wins)


def _inherit(x: _Site, parent: Call | None) -> Call | None:
    """Rule 0, inheritance: a host's E makes its inner elements E + host_empty (B never
    carried the host, so never what landed in it); a host's D/O/G/U passes its state + inherit. Tier, B
    site and U reason are the host's (its site is theirs). A host P/S/T leaves the element to the rules."""
    if parent is None or parent.state in PRESENT:
        return None
    x.ev["reason"] = parent.ev["reason"]
    flag = "host_empty" if parent.state == "E" else "inherit"
    return Call(x.e.eid, parent.state, parent.tier, [flag, *x.inv], parent.b_seqid, parent.b_start,
                parent.b_end, parent.b_strand, dict(x.ev))


def _window_rule(x: _Site) -> Call | None:
    """Rule 1, no usable alignment: window status not ok -> U (reason = the status); no alignment ->
    U no_alignment."""
    if x.ctx.win.status != "ok":
        return x.call("U", reason=x.ctx.win.status)
    if x.ctx.aln is None:
        return x.call("U", reason="no_alignment")
    return None


def _flank_rule(x: _Site) -> Call | None:
    """Rule 2, flanks: fewer than MIN_SEG aligned A bases on either side of the element -> U
    flank_unaligned (nothing places the site)."""
    aln = x.ctx.aln
    if aln.aligned_bp(0, x.s) < MIN_SEG or aln.aligned_bp(x.t, len(x.ctx.a)) < MIN_SEG:
        return x.call("U", reason="flank_unaligned")
    return None


def _gap_rule(x: _Site) -> Call | None:
    """Rule 3, assembly gap: an N in B between the images of the last aligned A base before the
    element and the first after it, widened by REFINE_PAD -> G. Runs before any refinement: LASTZ aligns N
    runs through as mismatches (--ambiguous=iupac), and refinement would read one as deleted flank."""
    aln, b = x.ctx.aln, x.ctx.b
    lo, hi = _b_last_before(aln, x.s), _b_first_from(aln, x.t)
    region = b[max(0, min(lo, hi) - REFINE_PAD):min(len(b), max(lo, hi) + REFINE_PAD + 1)]
    return x.call("G") if "N" in region or "n" in region else None


def _absent(x: _Site) -> Call:
    """Rule 4, absent branch (< MIN_SEG own bases aligned): U inv_site / flank_unaligned / ambiguous, or
    the gap around the element read as O, E or D (rules 4a-4h). E needs both balanced overhangs within
    -w..w: a gap ending more than w inside the element means B keeps part of it, which is D, never E."""
    aln, s, t, w = x.ctx.aln, x.s, x.t, x.w
    if _inverted_at_site(aln, s, t):
        return x.call("U", reason="inv_site")  # 4a: an inverted site looks like an O gap to gap_around
    g0 = aln.gap_around((s + t) // 2)  # 4b
    if g0 is None:
        return x.call("U", reason="flank_unaligned")
    if g0.a1 <= g0.a0:
        return x.call("U", reason="ambiguous")
    g = refine_gap(x.ctx.a, x.ctx.b, g0)  # 4c: leftmost exact placement + slide (input back if not refined)
    refined = g is not g0
    b_ins = g.b1 - g.b0
    dl, dr, sh = _balanced(s - g.a0, g.a1 - t, g.slide)  # 4d, 4e
    # 4f: within -w..w; below -w, B keeps part of the element; above w, only other elements explain it
    left = -w <= dl <= w or (dl > w and _explained(x, g, g.a0 + sh, s))
    right = -w <= dr <= w or (dr > w and _explained(x, g, t, g.a1 + sh))
    flags = ["adjacent"] if (dl > w and left) or (dr > w and right) else []
    ev = {"dl": dl, "dr": dr, "b_ins": b_ins}
    if b_ins >= MIN_SEG:  # 4g
        return x.call("O", flags, b=(g.b0, g.b1), **ev)
    if left and right:
        if b_ins > w:
            return x.call("O", [*flags, "small_ins"], b=(g.b0, g.b1), **ev)
        return _empty_site(x, g, flags, ev)
    return x.call("D", flags + (["filler"] if b_ins > w else []), evidence=refined, b=_junction(x, g), **ev)


def _inverted_at_site(aln: AlnMap, s: int, t: int) -> bool:
    """Rule 4a, inverted site: an inverted hit of >= MIN_SEG A bases overlapping the element [s, t) whose B
    interval meets the site [min(lo, hi), max(lo, hi)]; lo, hi = the B images of the last aligned A base
    before the element and the first from its end (rule 3's, not widened). An inverted copy elsewhere in the
    B window (another element's fragment, a paralogous LTR) says nothing about the site."""
    lo, hi = _b_last_before(aln, s), _b_first_from(aln, t)
    lo, hi = min(lo, hi), max(lo, hi)
    return any(a1 - a0 >= MIN_SEG and a0 < t and a1 > s and b0 <= hi and b1 > lo for a0, a1, b0, b1 in aln.inv)


def _balanced(dl: int, dr: int, slide: int) -> tuple[int, int, int]:
    """Rule 4e, balanced placement: every placement a0 .. a0 + slide of a deletion scores the same, so E/D
    must not depend on the one an aligner reports. Shift x = clamp(floor((dl - dr) / 2 + 0.5), 0, slide)
    evens the overhangs: returns (dl - x, dr + x, x)."""
    sh = min(max(math.floor((dl - dr) / 2 + 0.5), 0), slide)
    return dl - sh, dr + sh, sh


def _explained(x: _Site, g: Gap, lo: int, hi: int) -> bool:
    """Rule 4f, overhang explanation: overhang [lo, hi) is explained when other annotated elements deleted
    with e cover all but <= w bp of it, each widened by its own TSD length on both sides (the TSD between
    neighbors). Deleted with e: inside the gap at some equal placement, [g.a0, g.a1 + g.slide); never e
    or one of its ancestors (a host trivially covers its children's overhangs)."""
    spans = []
    for o, os_, ot in x.ctx.inside(g.a0, g.a1 + g.slide):
        if o.eid != x.e.eid and o.eid not in x.ancestors:
            n = len(o.tsd) if o.tsd else TSD_DEFAULT
            if max(lo, os_ - n) < min(hi, ot + n):
                spans.append((max(lo, os_ - n), min(hi, ot + n)))
    return hi - lo - sum(b - a for a, b in union(spans)) <= x.w


def _empty_site(x: _Site, g: Gap, flags: list[str], ev: dict) -> Call:
    """Rule 4h, E and its TSD read: B's copy is the exact image of an A TSD copy wholly outside the gap at
    a placement in its slide range: T2 = A[t, t + L) at the leftmost placement if g.a1 <= t; else T1 =
    A[s - L, s) at placement x = s - g.a0 (the gap starting at the element) if x <= g.slide; else unknown.
    tsd_mm = its fewest mismatches to either A copy; tsd_match if tsd_mm <= tsd_allowance. B site: the
    copy, or without one the two flank bases around the junction."""
    b, s, t, L = x.ctx.b, x.s, x.t, x.L
    if g.a1 <= t:
        lo = g.b1 + (t - g.a1)
    elif s - g.a0 <= g.slide:
        lo = g.b0 + (s - g.a0) - L
    else:
        lo = None
    copy = None if lo is None else _slice(b, lo, lo + L)
    mm = x.tsd_mm(copy)
    match = mm is not None and mm <= x.allowance
    site = (lo, lo + L) if copy is not None else _junction(x, g)
    return x.call("E", flags + (["tsd_match"] if match else []), evidence=match, b=site, tsd_b=copy,
                  tsd_mm=mm, **ev)


def _junction(x: _Site, g: Gap) -> tuple[int, int]:
    """B site of a junction: the last aligned B base before gap g and the first after it, [g.b0 - 1,
    g.b1 + 1) (2 bases for a pure deletion), kept inside the B window."""
    return max(0, g.b0 - 1), min(len(x.ctx.b), g.b1 + 1)


def _present(x: _Site) -> Call:
    """Rule 5, then rule 5a: the present call stands only if it is anchored at the site, else U unanchored."""
    call = _present_call(x)
    return call if _anchored(x) else x.call("U", reason="unanchored")


def _anchored(x: _Site) -> bool:
    """Rule 5a: some alignment run (AlnMap._runs, the runs gap_around reads) crosses the element's start s
    or its end t holding >= MIN_SEG aligned A bases on each side of it: in the flank (before s, after t) and
    among the element's own bases (after s, before t; children's spans left out). Failing that, the TSD
    exception: the element's B image (_ltr_image) lies strictly between lo and hi, the B images of the last
    aligned A base before the element and of the first from its end (rule 3's), and B's bases beside the
    image match an A TSD copy on both sides (_tsd_beside). The window's alignment is collinear, so only the
    extrapolated ends of that image can cross lo or hi. Residual: a paralogous copy whose own TSD bases align
    as short tails beside it passes, because lo and hi then sit at the copy itself."""
    aln, s, t = x.ctx.aln, x.s, x.t
    ra0, ra1 = aln._runs(MIN_SEG)[:2]
    for i in np.flatnonzero(((ra0 < s) & (ra1 > s)) | ((ra0 < t) & (ra1 > t))):
        a0, a1 = int(ra0[i]), int(ra1[i])
        if a0 < s < a1 and min(aln.aligned_bp(a0, s), _own(aln, s, min(a1, t), x.holes)[0]) >= MIN_SEG:
            return True
        if a0 < t < a1 and min(aln.aligned_bp(t, a1), _own(aln, max(a0, s), t, x.holes)[0]) >= MIN_SEG:
            return True
    image = _ltr_image(x)
    if image is None:
        return False
    lo, hi = _b_last_before(aln, s), _b_first_from(aln, t)
    return lo < image[0] and image[1] <= hi and _tsd_beside(x, *image)


def _present_call(x: _Site) -> Call:
    """Rule 5, present branch (>= MIN_SEG own bases aligned): P if both LTRs and the element are each
    >= COV_PRESENT aligned (nestB if B holds >= MIN_SEG bp of its own inside the element); S if the
    internal region is next to unaligned (< max(MIN_SEG, (1 - COV_PRESENT) * its length)) and the LTRs
    hold SOLO_LEN x one mean LTR length (solo_tsd: _solo_tsd); else T."""
    aln, s, t = x.ctx.aln, x.s, x.t
    site = _b_image(aln, s, t)
    cov = (x.ev["cov_left"], x.ev["cov_right"], x.ev["cov_elem"])
    if all(c is not None and c >= COV_PRESENT for c in cov):
        nest = any(g.b1 - g.b0 >= MIN_SEG for g in aln.b_gaps_within(s, t))
        return x.call("P", ["nestB"] if nest else [], evidence=True, b=site)
    (int_al, int_len), (ll_al, ll_len), (rl_al, rl_len) = x.parts["int"], x.parts["ltr_l"], x.parts["ltr_r"]
    mean_ltr = (ll_len + rl_len) / 2
    if (int_al < max(MIN_SEG, (1 - COV_PRESENT) * int_len)
            and SOLO_LEN[0] * mean_ltr <= ll_al + rl_al <= SOLO_LEN[1] * mean_ltr):
        tsd = _solo_tsd(x)
        return x.call("S", ["solo_tsd"] if tsd else [], evidence=tsd, b=site)
    return x.call("T", b=site)


def _solo_tsd(x: _Site) -> bool:
    """Rule 5, solo_tsd, read from B: the L bases on each side of the solo LTR in B both match an A TSD
    copy within tsd_allowance (_tsd_beside the LTR's image, _ltr_image). Read from B, not from how the TSD
    copies align: the reported deletion can reach into a copy B keeps, or align a copy B lacks."""
    image = _ltr_image(x)
    return image is not None and _tsd_beside(x, *image)


def _ltr_image(x: _Site) -> tuple[int, int] | None:
    """[bmin, bmax + 1): the B image of the element's LTR ends, the B positions of the annotated LTR ends
    extrapolated from the nearest aligned base inside the LTR (the first aligned LTR base back to its LTR's
    first base, the last forward to its LTR's last base), so a read beside it sits right beside the LTR even
    when chance microhomology leaves an end base unaligned; None if no LTR base is aligned. Computed once per
    site (rule 5's solo_tsd and rule 5a's TSD exception both read it)."""
    if x.image is not _UNSET:
        return x.image
    x.image = _image(x)
    return x.image


_UNSET = object()  # _Site.image not computed yet


def _image(x: _Site) -> tuple[int, int] | None:
    """_ltr_image, computed."""
    aln = x.ctx.aln
    first = last = None  # (A position, B position, the LTR's start / end)
    for lo, hi in x.ltrs:
        for p0, p1 in _pieces(lo, hi, x.holes):
            ends = _aligned_ends(aln, p0, p1)
            if ends is None:
                continue
            (af, bf), (al, bl) = ends
            if first is None or af < first[0]:
                first = (af, bf, lo)
            if last is None or al > last[0]:
                last = (al, bl, hi)
    if first is None:
        return None
    return first[1] - (first[0] - first[2]), last[1] + (last[2] - last[0])


def _tsd_beside(x: _Site, bmin: int, bend: int) -> bool:
    """The L bases of B on each side of B interval [bmin, bend) both match an A TSD copy within
    tsd_allowance."""
    b, L = x.ctx.b, x.L
    mms = [x.tsd_mm(side) for side in (_slice(b, bmin - L, bmin), _slice(b, bend, bend + L))]
    return all(mm is not None and mm <= x.allowance for mm in mms)


# ---------------------------------------------------------------- per-window and per-element context


class _Ctx:
    """The window: its sequences, alignment and elements (outermost-first, level by level)."""

    def __init__(self, win: Window, aln: AlnMap | None, a: str, b: str, elements: dict[str, Element],
                 children: dict[str, list[str]]):
        self.win, self.aln, self.a, self.b = win, aln, a, b
        self.elements, self.children = elements, children
        if win.status == "ok":
            if len(a) != win.a_end - win.a_start or len(b) != win.b_end - win.b_start:
                raise ValueError(f"window {win.wid}: sequence lengths A {len(a)}, B {len(b)} do not match "
                                 f"the window spans A {win.a_end - win.a_start}, B {win.b_end - win.b_start}")
            if win.strand not in (1, -1):
                raise ValueError(f"window {win.wid}: strand {win.strand} must be 1 or -1")
        self.members: list[Element] = []
        seen: set[str] = set()  # each element once, whatever duplicates or cycles `children` holds
        level = [self._get(eid) for eid in dict.fromkeys(win.eids)]
        while level:
            seen.update(e.eid for e in level)
            self.members.extend(level)
            level = [self._get(c) for c in dict.fromkeys(c for e in level for c in children.get(e.eid, ()))
                     if c not in seen]
        # window intervals of the members, by start, for rule 4f
        a0 = win.a_start
        spans = sorted((e.start - 1 - a0, e.end - a0, i) for i, e in enumerate(self.members))
        self._starts = [sp[0] for sp in spans]
        self._spans = [(self.members[i], s, t) for s, t, i in spans]

    def _get(self, eid: str) -> Element:
        if eid not in self.elements:
            raise ValueError(f"window {self.win.wid}: element {eid} is not in `elements`")
        return self.elements[eid]

    def inside(self, lo: int, hi: int) -> list[tuple[Element, int, int]]:
        """(element, s, t) of the members whose window interval lies in [lo, hi)."""
        i0, i1 = bisect_left(self._starts, lo), bisect_right(self._starts, hi)
        return [sp for sp in self._spans[i0:i1] if sp[2] <= hi]


class _Site:
    """One element in its window: [s, t), L, w, ancestors, own (children-excluded) parts and evidence."""

    def __init__(self, ctx: _Ctx, e: Element):
        win = ctx.win
        self.ctx, self.e = ctx, e
        self.s, self.t = e.start - 1 - win.a_start, e.end - win.a_start
        if not 0 <= self.s < self.t <= win.a_end - win.a_start:
            raise ValueError(f"element {e.eid} ({e.seqid}:{e.start}-{e.end}) lies outside window {win.wid} "
                             f"({win.a_seqid}:{win.a_start + 1}-{win.a_end})")
        self.L = len(e.tsd) if e.tsd else TSD_DEFAULT
        self.w = SLOP_TSD_MULT * self.L
        self.inv = ["inv"] if win.strand == -1 else []
        self.ancestors: set[str] = set()
        p = e.parent
        while p is not None and p in ctx.elements and p not in self.ancestors:
            self.ancestors.add(p)
            p = ctx.elements[p].parent
        self.ev: dict = dict.fromkeys(_EV_KEYS)
        left, right = win.anchor_distances(e.start, e.end)
        self.ev.update(tsd_a=e.tsd, anchor_left=int(left), anchor_right=int(right), margin=_num(win.margin))
        self.ident = math.nan
        self.own_aligned = 0
        self.parts: dict[str, tuple[int, int]] = {}
        self.ltrs: list[tuple[int, int]] = []  # the LTRs' window intervals, clipped to the element
        self.holes: list[tuple[int, int]] = []  # the children's spans, merged
        self.image: object = _UNSET  # _ltr_image, once computed
        if win.status == "ok" and ctx.aln is not None:
            self._measure(ctx.aln)
        self.placement = bool(_num(win.margin, 0.0) >= MARGIN_HIGH and self.ident >= IDENT_MIN)

    def _measure(self, aln: AlnMap) -> None:
        """Own aligned bp and length of the LTRs, internal region and whole element; coverages; ident."""
        a0 = self.ctx.win.a_start
        kids = [self.ctx.elements[c] for c in self.ctx.children.get(self.e.eid, ())]  # checked by _Ctx
        holes = self.holes = union([(c.start - 1 - a0, c.end - a0) for c in kids])
        ll = (self.e.ltr_left[0] - 1 - a0, self.e.ltr_left[1] - a0)
        rl = (self.e.ltr_right[0] - 1 - a0, self.e.ltr_right[1] - a0)
        parts = {"ltr_l": ll, "ltr_r": rl, "int": (ll[1], rl[0]), "elem": (self.s, self.t)}
        for name, (lo, hi) in parts.items():
            self.parts[name] = _own(aln, max(lo, self.s), min(hi, self.t), holes)
        self.ltrs = [(max(lo, self.s), min(hi, self.t)) for lo, hi in (ll, rl)]
        for key, name in zip(("cov_left", "cov_right", "cov_int", "cov_elem"), parts):
            al, n = self.parts[name]
            self.ev[key] = al / n if n > 0 else None
        self.own_aligned = self.parts["elem"][0]
        self.ident = aln.identity(self.s - IDENT_FLANK, self.t + IDENT_FLANK)
        self.ev["ident"] = _num(self.ident)

    @property
    def allowance(self) -> int:
        """TSD mismatches allowed at this element (tsd_allowance with its L, K and local identity)."""
        return tsd_allowance(self.L, self.e.k, self.ident)

    def tsd_mm(self, copy: str | None) -> int | None:
        """Fewest mismatches of B sequence `copy` to either A TSD copy, A[s - L, s) or A[t, t + L); None
        without a copy or with neither A copy inside the window."""
        a, L = self.ctx.a, self.L
        own = [c for c in (_slice(a, self.s - L, self.s), _slice(a, self.t, self.t + L)) if c is not None]
        return min(hamming(copy, c) for c in own) if copy is not None and own else None

    def call(self, state: str, flags: Sequence[str] = (), evidence: bool = False,
             b: tuple[int, int] | None = None, **ev) -> Call:
        """The call: tier from placement + evidence ("." for G, U), B site from window interval b."""
        self.ev.update(ev)
        tier = "." if state in ("G", "U") else ("high", "medium", "low")[2 - self.placement - bool(evidence)]
        site = _b_site(self.ctx.win, *b) if b is not None else (None, None, None, None)
        return Call(self.e.eid, state, tier, [*flags, *self.inv], *site, dict(self.ev))


# ---------------------------------------------------------------- helpers


def _pieces(lo: int, hi: int, holes: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """The parts of window interval [lo, hi) outside the holes (merged, sorted children spans)."""
    out = []
    for h0, h1 in holes:
        if h0 > lo:
            out.append((lo, min(hi, h0)))
        lo = max(lo, h1)
    out.append((lo, hi))
    return [(a, b) for a, b in out if a < b]


def _own(aln: AlnMap, lo: int, hi: int, holes: list[tuple[int, int]]) -> tuple[int, int]:
    """(aligned bp, length) of window interval [lo, hi) outside the holes (merged children spans)."""
    pcs = _pieces(lo, hi, holes)
    return sum(aln.aligned_bp(a, b) for a, b in pcs), sum(b - a for a, b in pcs)


def _b_last_before(aln: AlnMap, a: int) -> int:
    """B position of the last aligned A position < a (one exists: rule 2)."""
    i = int(np.searchsorted(aln.seg_a, a, side="left")) - 1
    last = min(int(aln.seg_a[i] + aln.seg_len[i]), a) - 1
    return int(aln.seg_b[i]) + last - int(aln.seg_a[i])


def _b_first_from(aln: AlnMap, a: int) -> int:
    """B position of the first aligned A position >= a (one exists: rule 2)."""
    i = int(np.searchsorted(aln.seg_a + aln.seg_len, a, side="right"))
    return int(aln.seg_b[i]) + max(int(aln.seg_a[i]), a) - int(aln.seg_a[i])


def _aligned_ends(aln: AlnMap, s: int, t: int) -> tuple[tuple[int, int], tuple[int, int]] | None:
    """((A, B) of the first aligned A position in [s, t), (A, B) of the last), or None if none is."""
    i0 = int(np.searchsorted(aln.seg_a + aln.seg_len, s, side="right"))  # first segment ending after s
    i1 = int(np.searchsorted(aln.seg_a, t, side="left")) - 1  # last segment starting before t
    if i0 > i1:
        return None
    af = max(int(aln.seg_a[i0]), s)
    al = min(int(aln.seg_a[i1] + aln.seg_len[i1]), t) - 1
    return ((af, int(aln.seg_b[i0]) + af - int(aln.seg_a[i0])),
            (al, int(aln.seg_b[i1]) + al - int(aln.seg_a[i1])))


def _b_image(aln: AlnMap, s: int, t: int) -> tuple[int, int] | None:
    """[min, max + 1) of the B positions aligned to A [s, t) (collinear: those of the first and last
    aligned A positions), or None if no base of it is aligned."""
    ends = _aligned_ends(aln, s, t)
    return None if ends is None else (ends[0][1], ends[1][1] + 1)


def _b_site(win: Window, lo: int, hi: int) -> tuple[str | None, int, int, str]:
    """(seqid, start, end, strand) in B's genome, 1-based inclusive, of window interval [lo, hi); on strand
    -1 window position j is genome position win.b_end - 1 - j, so the interval flips."""
    if win.strand == -1:
        return win.b_seqid, win.b_end - hi + 1, win.b_end - lo, "-"
    return win.b_seqid, win.b_start + lo + 1, win.b_start + hi, "+"


def _slice(seq: str, lo: int, hi: int) -> str | None:
    return seq[lo:hi] if 0 <= lo <= hi <= len(seq) else None


def _num(v, nan=None):
    """float(v), or `nan` for a missing value (None or NaN)."""
    return nan if v is None or math.isnan(v) else float(v)
