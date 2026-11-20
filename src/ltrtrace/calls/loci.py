"""Locus catalog: one locus per insertion event, across every genome pair.

Input: the pair call files (calls.records) that the routes write (genotype_pair in routes.assembly.pair),
call_paths[(A, B)] = <A>__<B>.calls.tsv.gz for every ordered pair with an annotated A (each element of A once,
with its state in B). They are read line by line; each call
keeps only its state, tier, flags, reason and B site, in compact per-pair columns (memory proportional to
the calls).

Linking: union-find over the annotated elements, keyed (genome, eid). A P, S or T call of A's element links
it to every annotated element of B (any level) on the call's B seqid whose span has reciprocal overlap
>= COV_PRESENT with the call's B interval (1-based inclusive). A locus is a connected set of elements, its
member elements.

Members and attributes: members[g] = the eid of genome g's first member element (fai seqid order, then
start); g's other member elements go to extra[g] and flag the locus multi_member. discovery: the genomes
with a member, in sheet order. family, superfamily, clade: each the majority over the member elements that
have one (a tie goes to the first in sheet, then position order). length, ltr_len: the medians, rounded
half up, of the member elements' spans and mean LTR lengths. tsd: the first member element's (sheet, then
position order) that has one.

States, one per genome, in sheet order:
  member genome  P, tier "annot", at its first member element (coordinates, strand, K2P); flag
                 annot_conflict if a member element of another genome calls this genome E, D or O.
  other genome   the candidates are the calls of the member elements into it; those of the highest tier
                 (high > medium > low > ".") decide. If they hold more than one genotype class (present
                 P/S/T, absent E, missing A/D/O/G/U: calls.states.call_class; build_loci's `absent` can move A
                 into absent) -> U at that tier, flag conflict, no site;
                 else their most frequent state (a tie: the order of calls.states.ALLELES: P > S > T,
                 E > A, D > O > G > U) with the site, flags and
                 reason of its first candidate (sheet, then position order of its element). No candidate
                 (the pair was not run) -> U, tier ".", flag no_call.
  The strand is the element's own in that genome: a member's annotated strand; for a call, its element's
  strand, reversed when the B site is (b_strand "-", an inverted window).

lids: L0000001, ... in the order of (sheet index of the first genome with coordinates, that genome's fai
seqid index, start, end, then the first member element: unique). Input order (element lists, call_paths,
rows) never changes them.

QC rows (scope, metric, value), scope by scope in this order, each scope's metrics in registry order
(_QC_METRICS: one line per metric function):
  all    n_loci; n_conflict and n_annot_conflict, the loci with a genome so flagged; n_multi_member.
  A>B    each ordered pair read (sheet order), over A's calls in B: frac_<state> for the 8 states,
         resolved_frac (P/S/T/E), high_frac (tier high); if B is annotated, annot_agree and e_in_annot
         (how far the calls agree with B's own annotation). A call's B site is [b_start, b_end] for
         P/S/T, E (its TSD copy) and D (the junction); O's B coordinates are the sequence B holds in the
         element's place, not a site, and G and U have none. A B element counts for a call of
         element a unless it is in the locus of an ancestor of a (a nested element's site lies inside its
         host's ortholog). annot_agree: among the calls with such a B element at the site, the fraction
         P/S/T (nan if none). A B element is at the site when it overlaps it (>= 1 bp) and does not span it;
         it spans the site when it reaches more than w bp past both of its ends (w = SLOP_TSD_MULT x a's TSD
         length, the classifier's tolerance for "at the element's ends"): such an element holds the site
         inside it, mostly a host annotated only in B, not an insertion at the site. e_in_annot: the E calls
         whose site overlaps (>= 1 bp) such a B element that overlaps (>= 1 bp) the B site of no ancestor of
         a called P, S or T in B (the host's present call places its ortholog there, even where B's
         annotation gives it other boundaries and the two never link into one locus), over all E calls (0 if
         none); spanning elements included, so on real data an upper bound on E error.
  A~B    each pair of annotated genomes read both ways. A link (a, b) is seen in A->B when a P/S/T call of a
         links it to b (the linking rule above), and in B->A likewise. reciprocity (A->B and
         B->A agree among loci resolved both ways): among the links seen either way whose two calls (a's in
         B, b's in A) are both resolved (P/S/T/E), the fraction where both are present (nan if none).
         link_symmetry: |links seen both ways| / |links seen either way| (nan without links), whether each
         projection finds the other's element: lower than reciprocity by the links one of whose calls is
         unresolved (D/O/G/U) and by boundary differences that fail the reciprocal overlap one way.
  genome each genome (sheet order): annot_sensitivity (annotated genomes): among the loci with a P call of
         tier high into the genome from some member element, the fraction where the genome has a member (nan
         if none); gap_rate: the fraction G among the calls into the genome (nan if none); skipped_<reason>:
         its annotation loader's skip counts (build_loci's `skipped`), nonzero ones, by reason.
"""
from __future__ import annotations

import logging
import math
import statistics
import time
from array import array
from bisect import bisect_left, bisect_right
from collections import Counter
from collections.abc import Callable, Iterable, Iterator
from dataclasses import dataclass, field
from itertools import accumulate
from operator import itemgetter
from pathlib import Path

import numpy as np

from ltrtrace.base.constants import COV_PRESENT, SLOP_TSD_MULT, TSD_DEFAULT
from ltrtrace.io.elements import Element
from ltrtrace.io.samples import PAIR_SEP, Sample

from .records import stream_calls
from .states import ALLELES, OBSERVED_ABSENT, PAIR_STATES, PRESENT, call_class

log = logging.getLogger("ltrtrace.loci")

QcRow = tuple[str, str, float]  # (scope, metric, value)

_TIERS = (".", "low", "medium", "high")  # call tiers, weakest first: the highest present wins at a locus
_TIER = {t: r for r, t in enumerate(_TIERS)}
_STATE = {s: c for c, s in enumerate(PAIR_STATES, 1)}  # state codes; 0 = no call
_LINKS = frozenset(_STATE[s] for s in PRESENT)  # B holds the element: its site may be a B element
_DENIES = frozenset(_STATE[s] for s in OBSERVED_ABSENT | {"D", "O"})  # B lacks it at the site: contradicts B
_SITED = _LINKS | {_STATE[s] for s in OBSERVED_ABSENT | {"D"}}  # states whose B coordinates are a site (annot_agree)
_RESOLVED = frozenset(_STATE[s] for s in PRESENT | OBSERVED_ABSENT)  # resolved calls (resolved_frac, reciprocity)
_STRAND = {".": 0, "+": 1, "-": 2}
_STRANDS = (None, "+", "-")
_PARSED = ("eid", "state", "tier", "flags", "reason", "b_seqid", "b_start", "b_end", "b_strand")
_LEFT = (1 - COV_PRESENT) / COV_PRESENT  # how far, per query bp, a partner may start before the query
_RIGHT = 1 - COV_PRESENT  # ... and after it


@dataclass(slots=True)
class GenomeState:
    """A genome's state at a locus. Coordinates: that genome's, 1-based inclusive; None without a site.

    tier: "annot" for a member genome (its annotation holds the insertion), else the winning call's tier
    (high | medium | low | "."). strand: the element's strand in this genome (module docstring). k: a member's
    K2P, None elsewhere. flags: the winning call's flags, or conflict | no_call | annot_conflict. reason: the
    winning call's reason (why a U is unresolved), None for a member, a conflict or no_call.
    """

    state: str
    tier: str
    seqid: str | None
    start: int | None
    end: int | None
    strand: str | None
    k: float | None
    flags: list[str]
    reason: str | None = None


@dataclass(slots=True)
class Locus:
    """One insertion event: its member elements (genome -> eid of the first; extra: the others), attributes,
    and every genome's state (sheet order). flags: multi_member when a genome holds more than one member."""

    lid: str
    members: dict[str, str]
    discovery: list[str]
    family: str | None
    superfamily: str | None
    clade: str | None
    length: int
    ltr_len: int
    tsd: str | None
    states: dict[str, GenomeState]
    extra: dict[str, list[str]] = field(default_factory=dict)
    flags: list[str] = field(default_factory=list)


def build_loci(samples: list[Sample], elements: dict[str, list[Element]],
               call_paths: dict[tuple[str, str], Path], seq_order: dict[str, dict[str, int]],
               skipped: dict[str, Counter] | None = None,
               absent: frozenset[str] | set[str] = OBSERVED_ABSENT) -> tuple[list[Locus], list[QcRow]]:
    """The loci of all annotated elements, in lid order, and their QC rows (module docstring).

    samples: the sheet (its order is the sheet order). elements: genome -> its annotated elements (an
    unannotated genome is absent). call_paths: (A, B) -> the pair file of A's calls in B. seq_order: genome
    -> seqid -> fai index, for every genome. skipped: genome -> its annotation loader's skip counts (QC
    skipped_<reason>). absent: the alleles that make a genotype class absent (call_class; default E alone, so a
    top-tier A and E conflict, as present and absent do; pass OBSERVED_ABSENT | INFERRED_ABSENT to treat A as
    an absence too). Fails on inputs that disagree: a pair file whose header is not CALL_COLUMNS, that
    calls an element A does not have, misses one or calls one twice, or names a B seqid not in seq_order.
    """
    t0 = time.perf_counter()
    cat = _Catalog(samples, elements, seq_order, absent)
    cat.check_pairs(call_paths)
    log.info("loci: %d element(s) of %d annotated genome(s), %d pair file(s)", len(cat.elems),
             sum(cat.annotated), len(call_paths))
    for (a, b), path in sorted(call_paths.items(), key=lambda kv: (cat.gi[kv[0][0]], cat.gi[kv[0][1]])):
        cat.read(a, b, Path(path))
    log.info("loci: %d call(s) read, %d link(s), %.1f s", cat.n_calls, cat.n_links, time.perf_counter() - t0)
    loci = cat.loci()
    qc = _qc_rows(_Qc(cat, loci, skipped or {}))
    log.info("loci: %s in %.1f s", ", ".join(f"{k} {v}" for s, k, v in qc if s == "all"),
             time.perf_counter() - t0)
    return loci, qc


# ---------------------------------------------------------------- QC: one registry line per metric function


class _Qc:
    """What the QC metrics read: the catalog (its calls, spans and finished union-find), the loci, the skip
    counts; the pairs read (sheet order) and the pairs of annotated genomes read both ways (A before B)."""

    def __init__(self, cat: _Catalog, loci: list[Locus], skipped: dict[str, Counter]):
        self.cat, self.loci, self.skipped = cat, loci, skipped
        n = len(cat.labels)
        self.pairs = [(a, b) for a in range(n) for b in range(n) if cat.calls[a][b] is not None]
        self.both = [(a, b) for a, b in self.pairs if a < b and cat.calls[b][a] is not None]
        self._roots: np.ndarray | None = None

    def roots(self) -> np.ndarray:
        """Each node's union-find root: its locus."""
        if self._roots is None:
            find, n = self.cat.uf.find, len(self.cat.elems)
            self._roots = np.fromiter((find(x) for x in range(n)), dtype=np.int64, count=n)
        return self._roots


def _ratio(k: int, n: int) -> float:
    return k / n if n else math.nan


def _locus_counts(q: _Qc) -> Iterator[tuple[str, float]]:
    """all: n_loci; n_conflict, the loci where the top-tier calls into some genome disagree on the genotype
    class; n_multi_member; n_annot_conflict, the loci where a member genome is called E, D or O by another
    genome's member element."""
    loci = q.loci
    yield "n_loci", len(loci)
    yield "n_conflict", sum(any("conflict" in s.flags for s in x.states.values()) for x in loci)
    yield "n_multi_member", sum("multi_member" in x.flags for x in loci)
    yield "n_annot_conflict", sum(any("annot_conflict" in s.flags for s in x.states.values()) for x in loci)


def _state_fractions(q: _Qc, a: int, b: int) -> Iterator[tuple[str, float]]:
    """A>B: frac_<state> (each state), resolved_frac (P/S/T/E) and high_frac (tier high) of A's calls in B."""
    pc = q.cat.calls[a][b]
    n = len(pc.state)
    counts = {s: pc.state.count(_STATE[s]) for s in PAIR_STATES}
    for s in PAIR_STATES:
        if s != "A" or counts[s]:  # A calls only come from the library route
            yield f"frac_{s}", _ratio(counts[s], n)
    yield "resolved_frac", _ratio(sum(counts[s] for s in PAIR_STATES if _STATE[s] in _RESOLVED), n)
    yield "high_frac", _ratio(pc.tier.count(_TIER["high"]), n)


def _annotation_agreement(q: _Qc, a: int, b: int) -> Iterator[tuple[str, float]]:
    """A>B, B annotated: annot_agree and e_in_annot (module docstring)."""
    cat = q.cat
    spans = cat.spans[b]
    if spans is None:
        return  # nothing to agree with
    pc, names, off, find = cat.calls[a][b], cat.names[b], cat.offset[a], cat.uf.find
    met = agree = n_e = in_annot = 0
    for i, st in enumerate(pc.state):
        if st == _STATE["E"]:
            n_e += 1
        if st not in _SITED or pc.seqid[i] < 0:
            continue
        s, e = pc.start[i], pc.end[i]
        hits = spans.overlapping(names[pc.seqid[i]], s, e)
        ancestors = cat.ancestors(off + i) if hits else []
        if ancestors:
            avoid = {find(x) for x in ancestors}
            hits = [h for h in hits if find(h[0]) not in avoid]
        if not hits:
            continue
        if st == _STATE["E"]:
            # the B sites of the ancestors called present on this B seqid: their orthologs' places in B
            images = [(pc.start[j], pc.end[j]) for j in (x - off for x in ancestors)
                      if pc.state[j] in _LINKS and pc.seqid[j] == pc.seqid[i]]
            if any(not any(lo <= i1 and i0 <= hi for i0, i1 in images) for _, lo, hi in hits):
                in_annot += 1
        tsd = cat.elems[off + i].tsd
        w = SLOP_TSD_MULT * (len(tsd) if tsd else TSD_DEFAULT)
        if all(lo < s - w and hi > e + w for _, lo, hi in hits):
            continue  # every B element here spans the site: none is at it
        met += 1
        if st in _LINKS:
            agree += 1
    yield "annot_agree", _ratio(agree, met)
    yield "e_in_annot", in_annot / n_e if n_e else 0.0


def _reciprocity(q: _Qc, a: int, b: int) -> Iterator[tuple[str, float]]:
    """A~B: reciprocity and link_symmetry (module docstring)."""
    cat = q.cat
    ab = cat.links(a, b)
    ba = {(x, y) for y, x in cat.links(b, a)}
    either = ab | ba
    sa, sb, oa, ob = cat.calls[a][b].state, cat.calls[b][a].state, cat.offset[a], cat.offset[b]
    both = [(sa[x - oa], sb[y - ob]) for x, y in either if sa[x - oa] in _RESOLVED and sb[y - ob] in _RESOLVED]
    yield "reciprocity", _ratio(sum(x in _LINKS and y in _LINKS for x, y in both), len(both))
    yield "link_symmetry", _ratio(len(ab & ba), len(either))


def _annot_sensitivity(q: _Qc, g: int) -> Iterator[tuple[str, float]]:
    """genome, annotated: among the loci with a P call of tier high into it, the fraction with a member in it."""
    cat = q.cat
    if not cat.annotated[g]:
        return
    roots = q.roots()
    member = np.zeros(len(roots), dtype=bool)
    member[roots[cat.offset[g]:cat.offset[g] + len(cat.index[g])]] = True
    projected = np.zeros(len(roots), dtype=bool)
    for a, by in enumerate(cat.calls):
        pc = by[g]
        if pc is None:
            continue
        sel = ((np.frombuffer(pc.state, dtype=np.uint8) == _STATE["P"])
               & (np.frombuffer(pc.tier, dtype=np.uint8) == _TIER["high"]))
        projected[roots[cat.offset[a]:cat.offset[a] + len(sel)][sel]] = True
    yield "annot_sensitivity", _ratio(int((projected & member).sum()), int(projected.sum()))


def _gap_rate(q: _Qc, g: int) -> Iterator[tuple[str, float]]:
    """genome: the fraction G among the calls into it."""
    pcs = [by[g] for by in q.cat.calls if by[g] is not None]
    yield "gap_rate", _ratio(sum(pc.state.count(_STATE["G"]) for pc in pcs), sum(len(pc.state) for pc in pcs))


def _skipped(q: _Qc, g: int) -> Iterator[tuple[str, float]]:
    """genome: skipped_<reason>, the annotation loader's nonzero skip counts, by reason."""
    counts = q.skipped.get(q.cat.labels[g]) or {}
    for reason in sorted(counts):
        if counts[reason]:
            yield f"skipped_{reason}", counts[reason]


# (scope kind, metric function): "all" takes the _Qc, "pair" and "both" also a pair (A, B) of sheet indexes,
# "genome" a sheet index; each yields (metric, value)
_QC_METRICS: tuple[tuple[str, Callable[..., Iterable[tuple[str, float]]]], ...] = (
    ("all", _locus_counts),
    ("pair", _state_fractions),
    ("pair", _annotation_agreement),
    ("both", _reciprocity),
    ("genome", _annot_sensitivity),
    ("genome", _gap_rate),
    ("genome", _skipped),
)


def _qc_rows(q: _Qc) -> list[QcRow]:
    """Every metric's rows, scope by scope: all, A>B (pairs read), A~B (pairs read both ways), each genome."""
    labels = q.cat.labels
    scopes = (("all", "all", ()),
              *(("pair", f"{labels[a]}>{labels[b]}", (a, b)) for a, b in q.pairs),
              *(("both", f"{labels[a]}~{labels[b]}", (a, b)) for a, b in q.both),
              *(("genome", label, (g,)) for g, label in enumerate(labels)))
    return [(scope, name, value) for kind, scope, args in scopes for k, metric in _QC_METRICS if k == kind
            for name, value in metric(q, *args)]


# ---------------------------------------------------------------- the catalog


class _Catalog:
    """The annotated elements as union-find nodes (genome by genome in sheet order), each genome's spans, and
    the calls read so far."""

    def __init__(self, samples: list[Sample], elements: dict[str, list[Element]],
                 seq_order: dict[str, dict[str, int]], absent: frozenset[str] | set[str] = OBSERVED_ABSENT):
        # genotype classes (call_class) by state code, and each state's rank among its class's (the earlier
        # letter of ALLELES wins a tie: P > S > T, E > A, D > O > G > U)
        self.klass = {c: call_class(s, absent=absent) for s, c in _STATE.items()}
        self.rank = {c: ALLELES.index(s) for s, c in _STATE.items()}
        self.labels = [s.label for s in samples]
        self.gi = {x: g for g, x in enumerate(self.labels)}  # label -> sheet index
        if len(self.gi) != len(self.labels):
            raise ValueError(f"duplicate sample label in {self.labels}")
        unknown = sorted(set(elements) - set(self.gi))
        if unknown:
            raise ValueError(f"elements of genome(s) not in the samples: {', '.join(unknown)}")
        missing = [x for x in self.labels if x not in seq_order]
        if missing:
            raise ValueError(f"seq_order lacks genome(s) {', '.join(missing)}")
        self.order = [seq_order[x] for x in self.labels]  # sheet index -> seqid -> fai index
        self.names = [{i: s for s, i in o.items()} for o in self.order]  # sheet index -> fai index -> seqid
        self.annotated = [elements.get(x) is not None for x in self.labels]
        self.elems: list[Element] = []  # node -> element
        self.genome: list[int] = []  # node -> sheet index
        self.offset: list[int] = []  # sheet index -> its first node
        self.index: list[dict[str, int]] = []  # sheet index -> eid -> node - offset
        self.spans: list[_Spans | None] = []  # sheet index -> its elements' spans (None: unannotated)
        for g, label in enumerate(self.labels):
            els = elements.get(label) or []
            index: dict[str, int] = {}
            for e in els:
                if e.seqid not in self.order[g]:
                    raise ValueError(f"element {e.eid}: seqid {e.seqid!r} is not in {label}'s seq_order "
                                     "(fai)")
                if e.eid in index:
                    raise ValueError(f"{label}: duplicate element id {e.eid}")
                index[e.eid] = len(index)
            first = len(self.elems)
            self.offset.append(first)
            self.index.append(index)
            self.spans.append(_Spans(els, first) if self.annotated[g] else None)
            self.elems.extend(els)
            self.genome.extend([g] * len(els))
        self.uf = _UnionFind(len(self.elems))
        # [A][B] (sheet indexes): A's calls in B, None if that pair was not read
        self.calls: list[list[_PairCalls | None]] = [[None] * len(self.labels) for _ in self.labels]
        self.flag_words, self.reason_words = _Words(), _Words()
        self.n_calls = self.n_links = 0

    def check_pairs(self, call_paths: dict[tuple[str, str], Path]) -> None:
        """Fail on a pair of unknown genomes, of one genome, or from an unannotated one; warn on missing
        ones."""
        for a, b in call_paths:
            if a not in self.gi or b not in self.gi:
                raise ValueError(f"call_paths pair {a}{PAIR_SEP}{b}: genome not in the samples")
            if a == b:
                raise ValueError(f"call_paths pair {a}{PAIR_SEP}{b}: a genome against itself")
            if not self.annotated[self.gi[a]]:
                raise ValueError(f"call_paths pair {a}{PAIR_SEP}{b}: {a} has no annotation, so no calls")
        missing = [f"{a}{PAIR_SEP}{b}" for g, a in enumerate(self.labels) if self.annotated[g]
                   for b in self.labels if b != a and (a, b) not in call_paths]
        if missing:
            log.warning("loci: no calls for %d pair(s) (%s): a genome that no member element of a locus was "
                        "called in is U no_call there", len(missing), ", ".join(missing))

    def read(self, a: str, b: str, path: Path) -> None:
        """Keep pair file (a, b)'s calls; link each P, S or T call to the B elements its site overlaps."""
        ga, gb = self.gi[a], self.gi[b]
        index, off, b_order, spans = self.index[ga], self.offset[ga], self.order[gb], self.spans[gb]
        pc = _PairCalls(len(index))
        flag_code, reason_code, union = self.flag_words.code, self.reason_words.code, self.uf.union
        links = 0
        for n, (eid, state, tier, flags, reason, seqid, start, end, strand) in stream_calls(path, _PARSED):
            i = index.get(eid)
            if i is None:
                raise ValueError(f"{path}: line {n}: {eid} is not an element of {a}: the calls do not match "
                                 "the annotation (rerun the pair)")
            if pc.state[i]:
                raise ValueError(f"{path}: line {n}: {eid} is called more than once")
            st, tr, sd = _STATE.get(state), _TIER.get(tier), _STRAND.get(strand)
            if st is None or tr is None or sd is None:
                raise ValueError(f"{path}: line {n}: unknown state, tier or b_strand ({state!r}, {tier!r}, "
                                 f"{strand!r})")
            if seqid == ".":
                if st in _LINKS:
                    raise ValueError(f"{path}: line {n}: a {state} call without a B site: not a pair call file "
                                     "(rerun the pair)")
                sq, s, e = -1, 0, 0
            else:
                sq = b_order.get(seqid)
                if sq is None:
                    raise ValueError(f"{path}: line {n}: B seqid {seqid!r} is not in {b}'s seq_order (fai): "
                                     f"the calls do not match {b}'s FASTA")
                try:
                    s, e = int(start), int(end)
                except ValueError:
                    s = e = 0
                if not 1 <= s <= e:
                    raise ValueError(f"{path}: line {n}: B site {seqid}:{start}-{end} is not two integers with "
                                     "1 <= start <= end")
            pc.set(i, st, tr, flag_code(flags), reason_code(reason), sq, s, e, sd)
            if st in _LINKS and sq >= 0 and spans is not None:
                for node in spans.partners(seqid, s, e):
                    union(off + i, node)
                    links += 1
        absent = pc.state.count(0)
        if absent:
            some = [eid for eid, i in index.items() if not pc.state[i]][:5]
            raise ValueError(f"{path}: {absent} element(s) of {a} have no call ({', '.join(some)}"
                             f"{' ...' if absent > 5 else ''}): the calls do not match the annotation (rerun "
                             "the pair)")
        self.calls[ga][gb] = pc
        self.n_calls += len(index)
        self.n_links += links
        log.debug("%s%s%s: %d call(s), %d link(s)", a, PAIR_SEP, b, len(index), links)

    def links(self, ga: int, gb: int) -> set[tuple[int, int]]:
        """The links (node in A, node in B) that pair (A, B)'s P, S and T calls make: read's rule, replayed
        from the kept calls (empty if the pair was not read or B is unannotated)."""
        pc, spans = self.calls[ga][gb], self.spans[gb]
        out: set[tuple[int, int]] = set()
        if pc is None or spans is None:
            return out
        names, off = self.names[gb], self.offset[ga]
        for i, st in enumerate(pc.state):
            if st in _LINKS and pc.seqid[i] >= 0:
                out.update((off + i, node) for node in spans.partners(names[pc.seqid[i]], pc.start[i], pc.end[i]))
        return out

    def ancestors(self, node: int) -> list[int]:
        """The nodes of the element's ancestors in its own genome: parent, grandparent, ... (as far as they are
        annotated)."""
        g = self.genome[node]
        index, off = self.index[g], self.offset[g]
        out: list[int] = []
        seen = {node}
        parent = self.elems[node].parent
        while parent is not None:
            i = index.get(parent)
            if i is None or off + i in seen:  # not annotated here, or a cycle
                break
            seen.add(off + i)
            out.append(off + i)
            parent = self.elems[off + i].parent
        return out

    def loci(self) -> list[Locus]:
        """Every connected set of elements as a locus, numbered and returned in lid order."""
        keyed = []
        for group in self.uf.groups():
            group.sort(key=self._key)
            locus = self._locus(group)
            keyed.append((self._lid_key(locus, self._key(group[0])), locus))
        keyed.sort(key=itemgetter(0))  # keys are unique: the first member element ends each
        for n, (_, locus) in enumerate(keyed, 1):
            locus.lid = f"L{n:07d}"
        return [locus for _, locus in keyed]

    def _key(self, node: int) -> tuple[int, int, int, int, str]:
        """Sheet, then position order of an element: (sheet index, fai seqid index, start, end, eid)."""
        e, g = self.elems[node], self.genome[node]
        return g, self.order[g][e.seqid], e.start, e.end, e.eid

    def _lid_key(self, locus: Locus, first: tuple) -> tuple:
        """(sheet index of the first genome with coordinates, its fai seqid index, start, end, the first
        member element's _key)."""
        for g, st in enumerate(locus.states.values()):  # sheet order
            if st.seqid is not None:
                return g, self.order[g][st.seqid], st.start, st.end, first
        raise AssertionError(f"locus of {first[-1]} has no coordinates")  # a member always has them

    def _locus(self, group: list[int]) -> Locus:
        """The locus of a connected set of nodes sorted by _key (lid unset)."""
        els = [self.elems[x] for x in group]
        src = [(self.genome[x], x - self.offset[self.genome[x]]) for x in group]  # (genome, index in it)
        members: dict[str, str] = {}
        extra: dict[str, list[str]] = {}
        first: dict[int, Element] = {}  # sheet index -> its first member element
        for (g, _), e in zip(src, els, strict=True):
            if g in first:
                extra.setdefault(self.labels[g], []).append(e.eid)
            else:
                first[g] = e
                members[self.labels[g]] = e.eid
        states = {label: self._member_state(g, first[g], src) if g in first else self._state(g, els, src)
                  for g, label in enumerate(self.labels)}
        family, superfamily, clade = (_majority(getattr(e, f) for e in els)
                                      for f in ("family", "superfamily", "clade"))
        return Locus(lid="", members=members, discovery=list(members), family=family, superfamily=superfamily,
                     clade=clade, length=_median(e.end - e.start + 1 for e in els),
                     ltr_len=_median(map(_mean_ltr, els)), tsd=next((e.tsd for e in els if e.tsd), None),
                     states=states, extra=extra, flags=["multi_member"] if extra else [])

    def _member_state(self, g: int, e: Element, src: list[tuple[int, int]]) -> GenomeState:
        """P (tier annot) at member genome g's first member element e; annot_conflict if a member element of
        another genome (src: each member element's (genome, index in it)) calls g E, D or O."""
        deny = False
        for a, i in src:
            pc = self.calls[a][g] if a != g else None
            if pc is not None and pc.state[i] in _DENIES:
                deny = True
                break
        return GenomeState("P", "annot", e.seqid, e.start, e.end, e.strand, e.k,
                           ["annot_conflict"] if deny else [])

    def _state(self, g: int, els: list[Element], src: list[tuple[int, int]]) -> GenomeState:
        """Genome g's state from the member elements' calls into it (module docstring): the calls of the
        highest tier decide; more than one genotype class among them is a conflict, else their most frequent
        state (a tie: self.rank) at its first candidate (els order). src: each member element's (genome,
        index in it)."""
        rank, top = -1, []  # the highest tier so far and its candidates: (position in els, calls, index)
        for k, (a, i) in enumerate(src):
            pc = self.calls[a][g]
            if pc is None:
                continue
            if pc.tier[i] > rank:
                rank, top = pc.tier[i], [(k, pc, i)]
            elif pc.tier[i] == rank:
                top.append((k, pc, i))
        if not top:
            return GenomeState("U", ".", None, None, None, None, None, ["no_call"])
        states = [pc.state[i] for _, pc, i in top]
        state = states[0]
        if any(s != state for s in states):
            if len({self.klass[s] for s in states}) > 1:
                return GenomeState("U", _TIERS[rank], None, None, None, None, None, ["conflict"])
            state = min(set(states), key=lambda s: (-states.count(s), self.rank[s]))
        best, pc, i = top[states.index(state)]
        words, reason = self.flag_words.words[pc.flags[i]], self.reason_words.words[pc.reason[i]]
        flags = words.split(",") if words else []
        if pc.seqid[i] < 0:
            return GenomeState(PAIR_STATES[state - 1], _TIERS[rank], None, None, None, None, None, flags, reason)
        return GenomeState(PAIR_STATES[state - 1], _TIERS[rank], self.names[g][pc.seqid[i]], pc.start[i],
                           pc.end[i], _turn(els[best].strand, _STRANDS[pc.strand[i]]), None, flags, reason)


class _PairCalls:
    """One pair file's calls as compact columns indexed by A's element (its index in A's annotation): state
    code (0 = no call), tier rank, flags and reason codes (_Words), B seqid (fai index; -1 = no site), start,
    end and strand code. About 30 bytes per call."""

    __slots__ = ("end", "flags", "reason", "seqid", "start", "state", "strand", "tier")

    def __init__(self, n: int):
        self.state, self.tier, self.strand = bytearray(n), bytearray(n), bytearray(n)
        self.flags, self.reason, self.seqid = array("i", [0]) * n, array("i", [0]) * n, array("i", [0]) * n
        self.start, self.end = array("q", [0]) * n, array("q", [0]) * n

    def set(self, i: int, state: int, tier: int, flags: int, reason: int, seqid: int, start: int, end: int,
            strand: int) -> None:
        self.state[i], self.tier[i], self.flags[i], self.reason[i] = state, tier, flags, reason
        self.seqid[i], self.start[i], self.end[i], self.strand[i] = seqid, start, end, strand


class _Words:
    """The few distinct strings of a column (flags, reason) as small ints; "." is 0 (None)."""

    def __init__(self):
        self.codes: dict[str, int] = {".": 0}
        self.words: list[str | None] = [None]

    def code(self, word: str) -> int:
        c = self.codes.get(word)
        if c is None:
            c = self.codes[word] = len(self.words)
            self.words.append(word)
        return c


class _Spans:
    """One genome's elements per seqid, sorted by start: the reciprocal-overlap partners of a B interval, and
    the elements overlapping one."""

    def __init__(self, elements: list[Element], first: int):
        """`elements` are nodes first, first + 1, ..."""
        rows: dict[str, list[tuple[int, int, int]]] = {}
        for node, e in enumerate(elements, first):
            rows.setdefault(e.seqid, []).append((e.start, e.end, node))
        # per seqid: starts, ends, nodes, and the running maximum of the ends (non-decreasing)
        self._by: dict[str, tuple[list[int], list[int], list[int], list[int]]] = {}
        for seqid, r in rows.items():
            r.sort()
            ends = [x[1] for x in r]
            self._by[seqid] = ([x[0] for x in r], ends, [x[2] for x in r], list(accumulate(ends, max)))

    def overlapping(self, seqid: str, start: int, end: int) -> list[tuple[int, int, int]]:
        """(node, start, end) of the elements overlapping [start, end] (1-based inclusive) by >= 1 bp.

        Those before the first whose running maximum end reaches `start` end before it (however long a
        nested host is), those after the last starting by `end` start after it: bisects bound the scan."""
        hit = self._by.get(seqid)
        if hit is None or end < start:
            return []
        starts, ends, nodes, reach = hit
        lo, hi = bisect_left(reach, start), bisect_right(starts, end)
        return [(nodes[k], starts[k], ends[k]) for k in range(lo, hi) if ends[k] >= start]

    def partners(self, seqid: str, start: int, end: int) -> list[int]:
        """The nodes whose span has reciprocal overlap >= COV_PRESENT with [start, end] (1-based inclusive).

        A partner starts within [start - (1 - f) / f * n, start + (1 - f) * n] (n = end - start + 1, f =
        COV_PRESENT): its bases outside the overlap are at most (1 - f) of its length, itself at most n / f;
        the query's, at most (1 - f) * n. A bisect on the starts bounds the scan (widened 1 bp for rounding).
        """
        hit = self._by.get(seqid)
        if hit is None or end < start:
            return []
        starts, ends, nodes, _ = hit
        n = end - start + 1
        lo = bisect_left(starts, start - math.ceil(_LEFT * n) - 1)
        hi = bisect_right(starts, start + math.ceil(_RIGHT * n) + 1)
        return [nodes[k] for k in range(lo, hi) if _reciprocal(start, end, starts[k], ends[k]) >= COV_PRESENT]


class _UnionFind:
    """Disjoint sets over 0 .. n - 1: union by size, path compression."""

    def __init__(self, n: int):
        self.parent = list(range(n))
        self.size = [1] * n

    def find(self, x: int) -> int:
        parent = self.parent
        root = x
        while parent[root] != root:
            root = parent[root]
        while x != root:  # compress: point the whole path at the root
            nxt = parent[x]
            parent[x] = root
            x = nxt
        return root

    def union(self, x: int, y: int) -> None:
        x, y = self.find(x), self.find(y)
        if x != y:
            if self.size[x] < self.size[y]:
                x, y = y, x
            self.parent[y] = x
            self.size[x] += self.size[y]

    def groups(self) -> list[list[int]]:
        """The sets, each in increasing node order."""
        out: dict[int, list[int]] = {}
        for x in range(len(self.parent)):
            out.setdefault(self.find(x), []).append(x)
        return list(out.values())


# ---------------------------------------------------------------- helpers


def _reciprocal(s1: int, e1: int, s2: int, e2: int) -> float:
    """The smaller of the overlap's fractions of [s1, e1] and of [s2, e2] (1-based inclusive); 0 if
    disjoint."""
    ov = min(e1, e2) - max(s1, s2) + 1
    return min(ov / (e1 - s1 + 1), ov / (e2 - s2 + 1)) if ov > 0 else 0.0


def _turn(strand: str, site: str | None) -> str | None:
    """An element's strand seen at a B site of orientation `site` ("-": reversed relative to A); None
    without a site."""
    if site is None:
        return None
    if site == "+" or strand not in ("+", "-"):
        return strand
    return "-" if strand == "+" else "+"


def _majority(values: Iterable[str | None]) -> str | None:
    """The most frequent value other than None, a tie going to the first seen; None if there is none."""
    counts: dict[str, int] = {}  # in first-seen order
    for v in values:
        if v is not None:
            counts[v] = counts.get(v, 0) + 1
    return max(counts, key=counts.__getitem__) if counts else None  # max keeps the first of equal counts


def _median(values: Iterable[float]) -> int:
    """The median, rounded half up."""
    return math.floor(statistics.median(values) + 0.5)


def _mean_ltr(e: Element) -> float:
    return (e.ltr_left[1] - e.ltr_left[0] + e.ltr_right[1] - e.ltr_right[0] + 2) / 2
