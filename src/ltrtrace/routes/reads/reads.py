"""Per-read allele evidence at one element in one read set (BAM).

Labels, one per read: F full, S solo LTR, E empty site, O other allele length, U unexplained junction clip, P present
(crosses a junction into an LTR only: full or solo).

Coordinates: 0-based half-open throughout, as pysam's. `Target` is the route's view of an element, converted from
the 1-based closed io Element once, at the route boundary (run.load_targets).
"""
from __future__ import annotations

from dataclasses import dataclass

import parasail

from ltrtrace.base.constants import ANCHOR as A
from ltrtrace.base.constants import (
    CLIP_MIN,
    CLIP_PROBE,
    JUNC_WIN,
    SLOP,
    SW_MIN_ALN,
    SW_MIN_ID,
    SW_SCORE,
    TAU_ABS,
    TAU_REL,
    TSD_DEFAULT,
)
from ltrtrace.seq.cigar import aligned_in, deletion_at, qpos_at

_SC = parasail.matrix_create("ACGT", SW_SCORE[0], SW_SCORE[1])


@dataclass(frozen=True, slots=True)
class Target:
    """One annotated element as the reads route genotypes it, 0-based half-open: `start` is the first base of the
    5'-most LTR, `end` the base after the 3'-most one; `lltr` and `rltr` are the left and right LTR (by position).
    `inner` holds the spans of the elements nested inside it (n_inner of them, at any depth); n_outer counts those
    that hold it."""

    id: str  # the annotation's element ID
    chrom: str
    start: int
    end: int
    strand: str
    lltr: tuple[int, int]
    rltr: tuple[int, int]
    family: str = "."
    clade: str = "."
    tsd: str = ""  # annotated or inferred from the flanks; "" when neither
    n_inner: int = 0
    n_outer: int = 0
    inner: tuple[tuple[int, int], ...] = ()

    @property
    def length(self) -> int:
        return self.end - self.start

    @property
    def ltr_mean(self) -> float:
        return ((self.lltr[1] - self.lltr[0]) + (self.rltr[1] - self.rltr[0])) / 2

    @property
    def internal(self) -> tuple[int, int]:
        return (self.lltr[1], self.rltr[0])

    @property
    def tsd_len(self) -> int:
        """The TSD's length, TSD_DEFAULT without one."""
        return len(self.tsd) if self.tsd else TSD_DEFAULT


def sw_identity(x: str, y: str) -> float:
    """Smith-Waterman identity of x within y over >= SW_MIN_ALN aligned columns (0 if shorter or empty)."""
    if len(x) < SW_MIN_ALN or len(y) < SW_MIN_ALN:
        return 0.0
    tb = parasail.sw_trace_scan_16(x, y, SW_SCORE[2], SW_SCORE[3], _SC).traceback
    aln = sum(c != " " for c in tb.comp)
    return sum(c == "|" for c in tb.comp) / aln if aln >= SW_MIN_ALN else 0.0


def _read_lacks(x, refpos: int) -> bool:
    """True if refpos falls inside a deletion or skip of >= TAU_ABS bp in alignment x (the read lacks that
    sequence, so it cannot anchor there)."""
    return deletion_at(x.cigartuples, x.reference_start, refpos, TAU_ABS) is not None


def _qpos(x, refpos: int) -> int | None:
    """The read position of alignment x at reference position refpos (seq.cigar.qpos_at)."""
    return qpos_at(x.cigartuples, x.reference_start, refpos)


class Site:
    """Expected allele geometry of one element (reference = full allele)."""

    def __init__(self, el: Target, tsd_len: int):
        self.el = el
        self.s, self.e = el.start, el.end
        self.a, self.b = el.lltr[1], el.rltr[0]  # end of left LTR, start of right LTR
        self.L = self.e - self.s
        self.t = tsd_len
        l5, l3 = el.lltr[1] - el.lltr[0], el.rltr[1] - el.rltr[0]
        span_ref = self.L + 2 * A
        self.span_ref = span_ref
        # expected length change (query - reference span between the two anchors) and tolerance per allele.
        # A solo LTR keeps one of the two LTRs (they can differ in length): any change between the two is solo.
        self.delta = {"F": 0.0, "E": -(self.L + self.t)}
        self.solo = (-(self.L - min(l5, l3)), -(self.L - max(l5, l3)))  # (most negative, least negative)
        spans = {"F": span_ref, "S": max(l5, l3) + 2 * A, "E": 2 * A - self.t}
        self.tau = {k: max(TAU_ABS, TAU_REL * v) for k, v in spans.items()}
        self.slop = SLOP
        self.tsd_found = bool(el.tsd)  # ends confirmed by a TSD; others may be off by the annotation error

    def classify_delta(self, d: float, extra: float = 0.0) -> str:
        """Nearest allele for a length change d; 'extra' widens the empty and solo windows (annotation boundary
        error: it shifts their expected lengths, not the full one, whose reads match the reference)."""
        lo, hi = self.solo
        tau = {"F": self.tau["F"], "E": self.tau["E"] + extra, "S": self.tau["S"] + extra}
        z = {k: abs(d - self.delta[k]) / tau[k] for k in "FE"}
        z["S"] = 0.0 if lo <= d <= hi else min(abs(d - lo), abs(d - hi)) / tau["S"]
        k = min(z, key=z.get)
        return k if z[k] <= 1.0 else "O"


def carries_left(x, st: Site) -> bool:
    """Alignment x carries the left element end: the junction base and the end of the left LTR aligned (not inside
    a deletion of >= TAU_ABS) and >= ANCHOR aligned bp in the JUNC_WIN bp of internal region after the LTR."""
    s, a, b = st.s, st.a, st.b
    return (aligned_in(x.cigartuples, x.reference_start, a, min(b, a + JUNC_WIN)) >= min(A, b - a)
            and not _read_lacks(x, s) and not _read_lacks(x, a - 1))


def carries_right(x, st: Site) -> bool:
    """Alignment x carries the right element end (mirror of carries_left)."""
    e, a, b = st.e, st.a, st.b
    return (aligned_in(x.cigartuples, x.reference_start, max(a, b - JUNC_WIN), b) >= min(A, b - a)
            and not _read_lacks(x, e - 1) and not _read_lacks(x, b))


def _probe(fasta, chrom: str, clip: str, far: tuple[int, int], near: tuple[int, int]) -> str | None:
    """Label a junction clip: 'far' if it continues into the far side, 'near' if it matches the element, else None."""
    i_far = sw_identity(clip, fasta.fetch(chrom, max(0, far[0]), far[1]).upper())
    i_near = sw_identity(clip, fasta.fetch(chrom, max(0, near[0]), near[1]).upper())
    if i_far >= SW_MIN_ID and i_far >= i_near:
        return "far"
    if i_near >= SW_MIN_ID:
        return "near"
    return None


def evidence(bam, fasta, el: Target, tsd_len: int, min_mapq: int, reads_out: list | None = None) -> dict:
    """Classify every read at element el. Returns a dict of counts and the other-allele length changes.
    If reads_out is a list, (label, alignments) of every labeled read is appended to it (for read phasing)."""
    st = Site(el, tsd_len)
    s, e, chrom = st.s, st.e, el.chrom
    pieces, lowq = {}, set()
    sl = st.slop
    for lo, hi in ((s - sl, s + sl), (e - sl, e + sl)):
        for x in bam.fetch(chrom, max(0, lo), max(1, hi)):
            if x.is_unmapped or x.is_secondary or x.is_qcfail or x.is_duplicate:
                continue
            if x.mapping_quality < min_mapq:
                lowq.add(x.query_name)  # placed here, but ambiguously (repeat): not evidence either way
                continue
            key = (x.reference_start, x.flag)
            lst = pieces.setdefault(x.query_name, {})
            lst[key] = x
    c = {"F": 0, "F_span": 0, "F_other": 0, "F_left": 0, "F_right": 0, "S": 0, "S_split": 0, "E": 0, "E_split": 0,
         "O": 0, "U": 0, "P": 0}
    other_d, es_d = [], []
    for d in pieces.values():
        alns = list(d.values())
        # an anchor inside a large alignment deletion is not covered by the read: that read comes from a haplotype
        # with another structure there (e.g. a larger deletion) and is counted as 'other', never as F/S/E
        # ... but only if that deletion reaches into the element (> SLOP bp): a deletion that merely abuts it (a
        # polymorphic neighbor) says nothing about this element, and the read is used through its valid anchor
        cover = [x for x in alns if x.reference_start <= s - A and x.reference_end >= e + A]
        dels = [d for x in cover for d in (deletion_at(x.cigartuples, x.reference_start, s - A, TAU_ABS),
                                           deletion_at(x.cigartuples, x.reference_start, e + A - 1, TAU_ABS)) if d]
        if any(d[0] < e - sl and d[1] > s + sl for d in dels):
            c["O"] += 1
            if reads_out is not None:
                reads_out.append(("O", alns))
            continue
        left = [x for x in alns if x.reference_start <= s - A < x.reference_end and not _read_lacks(x, s - A)]
        right = [x for x in alns if x.reference_start < e + A <= x.reference_end and not _read_lacks(x, e + A - 1)]
        both = [x for x in left if any(x is y for y in right)]
        lab, dlt, kind = None, None, None
        if both:  # one alignment spans the site
            x = both[0]
            dlt = (_qpos(x, e + A - 1) + 1) - _qpos(x, s - A) - st.span_ref
            lab, kind = st.classify_delta(dlt), "span"
            if lab == "O" and carries_left(x, st) and carries_right(x, st):
                lab, kind = "F", "span_other"  # both element ends present: full allele, other internal length
        elif left and right:  # two pieces of one read, joined across the site
            x1, x2 = left[0], right[0]
            if x1.is_reverse == x2.is_reverse:
                q1, q2 = _qpos(x1, s - A), _qpos(x2, e + A - 1) + 1
                if q2 > q1:
                    dlt = (q2 - q1) - st.span_ref
                    lab, kind = st.classify_delta(dlt), "split"
                    if lab == "O" and carries_left(x1, st) and carries_right(x2, st):
                        lab, kind = "F", "span_other"
            if lab is None:
                lab, kind = "U", "split"
        elif left or right:  # anchored on one side only
            lab, kind = _one_sided(st, fasta, chrom, left or right, "L" if left else "R"), "one"
        if lab is None:
            continue
        if reads_out is not None:
            reads_out.append(("F" if lab in ("F_left", "F_right") else lab, alns))
        if lab in ("F_left", "F_right"):
            c[lab] += 1
            c["F"] += 1
            continue
        c[lab] += 1
        if lab == "F":
            c["F_span"] += 1
            if kind == "span_other":
                c["F_other"] += 1
        elif lab in ("E", "S") and kind != "span":
            c[lab + "_split"] += 1
        elif lab == "O":
            other_d.append(int(dlt))
        if lab in ("E", "S") and dlt is not None:
            es_d.append((lab, int(dlt)))
    c["other_delta"] = other_d
    c["es_delta"] = es_d  # length changes of spanning empty/solo reads (boundary precision)
    c["lowq"] = len(lowq - set(pieces))  # reads at the junctions with only low-MAPQ alignments
    c["highq"] = len(pieces)
    return c


def _one_sided(st: Site, fasta, chrom: str, xs: list, side: str) -> str | None:
    """A read anchored on one side only: junction crossing (F/P) or one piece of a split read (E/S/U)."""
    s, e, a, b = st.s, st.e, st.a, st.b
    SL = st.slop
    best = None
    for x in xs:
        ct = x.cigartuples
        if side == "L":
            p = x.reference_end
            clip_op, clip_len = ct[-1] if ct[-1][0] in (4, 5) else (None, 0)
            at = "s" if abs(p - s) <= SL else ("a" if abs(p - a) <= SL else None)
            if at and abs(p - a) < abs(p - s):
                at = "a"
            if at and clip_len >= CLIP_MIN:
                if clip_op == 5 or x.query_sequence is None:
                    return "U"
                q = x.query_sequence[x.query_alignment_end:x.query_alignment_end + CLIP_PROBE].upper()
                near = (s, s + CLIP_PROBE + SL) if at == "s" else (a, a + CLIP_PROBE + SL)
                r = _probe(fasta, chrom, q, (e - SL, e + CLIP_PROBE + SL), near)
                return {"far": "E" if at == "s" else "S"}.get(r, "U")
            # full: the junction base and the LTR end aligned, and >= A aligned bases of internal region next to
            # the LTR (a read crossing a deletion of the element, or a solo allele, has none there)
            if carries_left(x, st):
                best = "F"
            elif p >= s + A and not _read_lacks(x, s + A - 1) and best is None:
                best = "P"
        else:
            p = x.reference_start
            clip_op, clip_len = ct[0] if ct[0][0] in (4, 5) else (None, 0)
            at = "e" if abs(p - e) <= SL else ("b" if abs(p - b) <= SL else None)
            if at and abs(p - b) < abs(p - e):
                at = "b"
            if at and clip_len >= CLIP_MIN:
                if clip_op == 5 or x.query_sequence is None:
                    return "U"
                qs = x.query_alignment_start
                q = x.query_sequence[max(0, qs - CLIP_PROBE):qs].upper()
                near = (e - CLIP_PROBE - SL, e) if at == "e" else (b - CLIP_PROBE - SL, b)
                r = _probe(fasta, chrom, q, (s - CLIP_PROBE - SL, s + SL), near)
                return {"far": "E" if at == "e" else "S"}.get(r, "U")
            if carries_right(x, st):
                best = "F"
            elif p <= e - A and not _read_lacks(x, e - A) and best is None:
                best = "P"
    if best == "F":
        return "F_left" if side == "L" else "F_right"
    return best
