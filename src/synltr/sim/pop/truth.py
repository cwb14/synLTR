"""Truth, annotations, introgression tracts and the insertion table of `synltr sim pop`, read from the
genomes popbuild built (PopGenome and its Copy records); nothing the builder decided is derived again here.

alive(c): an alive element: fate intact or partial, and each LTR keeps >= COV_PRESENT of its own bases.

Truth (pair_truth): each alive element e of X, read in Y; an element counts as present in Y when Y inherited it (Y
lies below e's origin: Copy.inherited). cx and cy are e's copies in X and Y.
  G   an N in Y within REFINE_PAD of e's image in Y (cy.start, cy.end: a replaced copy's is its replacing bases) or
      of its site there (cy.site, its target bases [site, site + TSD_DEFAULT))
  e not inherited by Y:
      with a simulated host Y did not inherit (Copy.host, up the chain in X): the outermost such host's truth,
      plus host_empty when that is E (a dead host's site test runs as an alive one's)
      O   Y holds another insertion (a copy with an image) whose target bases overlap e's: same chrom and host
          system, positions less than TSD_DEFAULT apart (target bases [pos - TSD_DEFAULT, pos), popevents)
      E   cy.site >= 0: Y's empty-site window is intact (popbuild: only removals of its bases break it; SNPs,
          small indels and insertions at other sites in it never do)
      D   otherwise (there is no site_lost state: a broken window, or a host Y holds but lost, reads D)
  e inherited by Y, against X's own bases ("Y holds >= COV_PRESENT of e's own ids", read through
  Copy.kept_runs: held = own bases both copies keep; a loss X and Y share takes the same bases from both, so it
  leaves e P, where Y's kept share of its own bases alone would read T):
      S   cy.fate solo, and held >= COV_PRESENT x the shorter of X's two kept own LTRs
      O   cy.fate replaced: unrelated sequence of at least MIN_SEG at the intact site (classifier rule 4)
      P   Y keeps >= COV_PRESENT of X's kept own bases, and of each of X's kept own LTR bases
      T   held >= MIN_SEG
      D   otherwise
  Inversions are not simulated, so no truth carries inv. Columns are score.TRUTH_COLUMNS, which `synltr sim
  score` reads: name = the eid (the annotation's Name), start/end 1-based inclusive in X, origin backbone or sim,
  events_x and events_y the copies' events. Named limits: (1) a copy of e that Y holds without an image (deleted),
  or an absent one whose window is broken, has no site in Y, so it is never G (its flanks are not read there); (2) O
  does not ask whether X carries the other insertion too: of two insertions with overlapping targets that X both
  carries, the one a Y lacks reads O there when Y holds the other.

Annotation (annotate): the alive copies in (chrom, start, -end, eid) order, through donors.jitter_bounds and
donors.read_tsd: per copy, dropped with probability dropout, and with probability jitter its
start and end moved by U{-3..3} (two uniforms, then two integers); the TSD read beside the
reported ends; K2P_d = Copy.k2p; family and superfamily the donor's or the source's; Name = the eid.

Tracts (tract_bed): a sample's backbone tracts [l, r) lifted to its genome as [lift_pos(l), lift_pos(r)). Every
backbone base of a tract lies in its record and every other backbone base outside it; a tract whose bases are all
removed is left out. Payloads (sequence anchored before a base pos: insertions, replacing and indel bases) between
kept ends lie in the record. At a kept end a payload is read by the anchor convention, base pos - 1:
one anchored exactly at r lies in the tract and one anchored exactly at l does not, the reverse of its genealogy
(popevents draws an insertion's tree at base pos). Named limit: a tract end that falls in a stretch the sample
deleted lifts to the next kept base, so a payload anchored at that base falls on the wrong side: inside the record
at a right end, outside it at a left end.

Insertion table (insertion_rows): one row per insertion in INSERTION_COLUMNS, then fate_<label> per genome (its
Copy.fate there): chrom is the seqid; carriers are the labels of the samples below its event, survivors those whose
copy is alive (comma-separated, in sample order); host, family and survivors are "." when there are none. pos is
Ins.pos, 0-based as in BED: the copy goes before base pos of the chromosome's backbone (0 at its region's start) when
host is ".", else before offset pos of its host's donor body (in donor orientation), so pos is a coordinate of no
output genome; time is the event's, in generations ago, and node the tree node below it.
"""
from __future__ import annotations

import bisect
from collections import defaultdict
from typing import TYPE_CHECKING

import numpy as np

from synltr.base.constants import COV_PRESENT, MIN_SEG, REFINE_PAD, SIM_JITTER_BP, TSD_DEFAULT
from synltr.io.elements import Element, make_eid
from synltr.sim.pop.build import Copy, PopGenome, lift_pos
from synltr.sim.pop.donors import jitter_bounds, read_tsd
from synltr.sim.score import TRUTH_COLUMNS

if TYPE_CHECKING:
    from synltr.sim.pop.events import Ins

INSERTION_COLUMNS = ("iid", "chrom", "host", "pos", "time", "node", "depth", "family", "carriers", "survivors")
_LIVE = ("intact", "partial")


def alive(c: Copy) -> bool:
    """An alive element: fate intact or partial, and each LTR keeps >= COV_PRESENT of its own bases."""
    return c.fate in _LIVE and all(k >= COV_PRESENT * o for k, o in zip(c.ltr_kept, c.ltr_own))


def _order(c: Copy) -> tuple:
    return c.chrom, c.start, -c.end, c.eid


def _overlap(a: tuple[tuple[int, int], ...], b: tuple[tuple[int, int], ...]) -> int:
    """Bases two sorted lists of disjoint half-open runs share."""
    i = j = n = 0
    while i < len(a) and j < len(b):
        lo, hi = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if hi > lo:
            n += hi - lo
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return n


def _size(runs: tuple[tuple[int, int], ...]) -> int:
    return sum(hi - lo for lo, hi in runs)


class _Truth:
    """pair_truth's rules for X's elements in Y (module docstring); state(eid) is memoized."""

    def __init__(self, x: PopGenome, y: PopGenome):
        self.x, self.y = x, y
        self.held: dict[tuple[int, str | None], list[tuple[int, str]]] = defaultdict(list)
        for c in y.copies.values():  # Y's simulated insertions with an image, by (chrom, host system)
            chrom, host, pos = c.target
            if pos >= 0 and c.start >= 0:
                self.held[(chrom, host)].append((pos, c.eid))
        for v in self.held.values():
            v.sort()
        self.memo: dict[str, tuple[str, tuple[str, ...]]] = {}

    def gapped(self, c: Copy) -> bool:
        if c.start >= 0:
            a, b = c.start, c.end
        elif c.site >= 0:
            a, b = c.site, c.site + TSD_DEFAULT
        else:
            return False
        seq = self.y.seqs[c.chrom]
        return seq.find(b"N", max(a - REFINE_PAD, 0), b + REFINE_PAD) >= 0

    def outer_host(self, cx: Copy) -> str | None:
        """The outermost simulated host of cx (X's host chain) that Y did not inherit, or None."""
        out, h = None, cx.host
        while h is not None:
            if not self.y.copies[h].inherited:
                out = h
            h = self.x.copies[h].host
        return out

    def overlap(self, cx: Copy) -> bool:
        chrom, host, pos = cx.target
        v = self.held.get((chrom, host), ())
        k = bisect.bisect_left(v, (pos - TSD_DEFAULT + 1, ""))
        return any(eid != cx.eid for _, eid in v[k:bisect.bisect_left(v, (pos + TSD_DEFAULT, ""))])

    def state(self, eid: str) -> tuple[str, tuple[str, ...]]:
        if eid not in self.memo:
            self.memo[eid] = self._state(eid)
        return self.memo[eid]

    def _state(self, eid: str) -> tuple[str, tuple[str, ...]]:
        cx, cy = self.x.copies[eid], self.y.copies[eid]
        if self.gapped(cy):
            return "G", ()
        if not cy.inherited:
            host = self.outer_host(cx)
            if host is not None:
                state, flags = self.state(host)
                return state, (*flags, "host_empty") if state == "E" else flags
            if self.overlap(cx):
                return "O", ()
            return ("E", ()) if cy.site >= 0 else ("D", ())
        (xl, xo, xr), (yl, yo, yr) = cx.kept_runs, cy.kept_runs
        held = _overlap(xo, yo)
        ltr = min(_size(xl), _size(xr))
        if cy.fate == "solo" and held >= COV_PRESENT * ltr:
            return "S", ()
        if cy.fate == "replaced":
            return "O", ()
        frac = [_overlap(a, b) / _size(a) if _size(a) else 0.0 for a, b in ((xo, yo), (xl, yl), (xr, yr))]
        if min(frac) >= COV_PRESENT:
            return "P", ()
        if held >= MIN_SEG:
            return "T", ()
        return "D", ()


def _check_pair(x: PopGenome, y: PopGenome) -> None:
    if x.label == y.label:
        raise ValueError(f"the two genomes need different labels, got {x.label!r} twice")
    if x.seqids != y.seqids:
        raise ValueError(f"genomes {x.label} and {y.label} hold different chromosomes ({x.seqids} and {y.seqids})")
    missing = x.copies.keys() - y.copies.keys()
    if missing:
        raise ValueError(f"genome {y.label} has no copy of {len(missing)} element(s) of {x.label} (first: "
                         f"{min(missing)}); were they built from one set of events?")


def pair_truth(x: PopGenome, y: PopGenome) -> list[dict]:
    """TRUTH_COLUMNS rows, one per alive element of x in (chrom, start, -end, eid) order: its state in y by the rules
    of the module docstring. Genomes with one label, different chromosomes or missing copies raise ValueError."""
    _check_pair(x, y)
    t = _Truth(x, y)
    rows = []
    for c in sorted((c for c in x.copies.values() if alive(c)), key=_order):
        state, flags = t.state(c.eid)
        rows.append(dict(zip(TRUTH_COLUMNS, (
            x.label, c.eid, c.eid, x.seqids[c.chrom], c.start + 1, c.end, "backbone" if c.target[2] < 0 else "sim",
            state, ",".join(flags), ",".join(c.events), ",".join(y.copies[c.eid].events)))))
    return rows


def annotate(g: PopGenome, dropout: float, jitter: float, rng: np.random.Generator) -> list[Element]:
    """g's alive copies as an annotator would report them (module docstring), genome g.label, eids by make_eid."""
    out = []
    for c in sorted((c for c in g.copies.values() if alive(c)), key=_order):
        drop, jit = rng.random(2) < (dropout, jitter)
        ds, dt = (int(v) for v in rng.integers(-SIM_JITTER_BP, SIM_JITTER_BP + 1, 2))
        if drop:
            continue
        seq, seqid = g.seqs[c.chrom], g.seqids[c.chrom]
        s, t, a, d = jitter_bounds(c.start, c.end, c.ltr_left, c.ltr_right, len(seq), jit, ds, dt)
        out.append(Element(genome=g.label, eid=make_eid(g.label, seqid, s + 1, t), name=c.eid, seqid=seqid,
                           start=s + 1, end=t, strand=c.strand, ltr_left=(a + 1, c.ltr_left[1]),
                           ltr_right=(c.ltr_right[0] + 1, d), tsd=read_tsd(seq, s, t), k=c.k2p, family=c.family,
                           superfamily=c.superfamily, clade=None, confidence=None))
    return out


def tract_bed(g: PopGenome, chrom_tracts: dict[int, list[tuple[int, int]]]) -> list[tuple[str, int, int]]:
    """(seqid, start, end) records, 0-based half-open, of g's introgressed tracts in its genome (module docstring),
    by chromosome index then tract: chrom_tracts maps a chromosome index of g to this sample's backbone tracts
    (popmodel.tracts of that chromosome, at the sample's index). Backbone bases land exactly; a payload anchored at
    a kept end goes by the anchor convention, and one at the next kept base after a deleted end lands on the wrong
    side of it (the named limit)."""
    out = []
    for ch in sorted(chrom_tracts):
        if not 0 <= ch < len(g.seqids):
            raise ValueError(f"genome {g.label} has no chromosome index {ch} (it has {len(g.seqids)})")
        for lo, hi in chrom_tracts[ch]:
            a, b = lift_pos(g.lift[ch], lo), lift_pos(g.lift[ch], hi)
            if b > a:
                out.append((g.seqids[ch], a, b))
    return out


def insertion_rows(genomes: list[PopGenome], ins: dict[int, list[Ins]]) -> list[dict]:
    """One row per insertion (module docstring), by chromosome index then the order given: `genomes` are every
    sample in sample order (Ins.carriers indexes them), with the same chromosomes; `ins` maps a chromosome index
    to the insertions build_chrom built there. A missing copy, carrier or chromosome raises ValueError."""
    if not genomes or any(g.seqids != genomes[0].seqids for g in genomes):
        raise ValueError("insertion_rows needs genomes with the same chromosomes")
    seqids, labels = genomes[0].seqids, [g.label for g in genomes]
    rows = []
    for ch in sorted(ins):
        if not 0 <= ch < len(seqids):
            raise ValueError(f"no chromosome index {ch}: the genomes have {len(seqids)}")
        for x in ins[ch]:
            if x.carriers and max(x.carriers) >= len(genomes):
                raise ValueError(f"insertion {x.iid}: carrier {max(x.carriers)} is not one of the {len(genomes)} "
                                 "genomes")
            copies = [g.copies.get(x.iid) for g in genomes]
            if any(c is None or c.chrom != ch for c in copies):
                raise ValueError(f"insertion {x.iid}: not built on chromosome index {ch} of every genome")
            row = dict(zip(INSERTION_COLUMNS, (
                x.iid, seqids[ch], x.host or ".", x.pos, x.time, x.node, x.depth, copies[0].family or ".",
                ",".join(labels[s] for s in sorted(x.carriers)),
                ",".join(lab for lab, c in zip(labels, copies) if alive(c)) or ".")))
            row.update({f"fate_{lab}": c.fate for lab, c in zip(labels, copies)})
            rows.append(row)
    return rows
