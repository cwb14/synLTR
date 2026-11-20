"""1:1 orthologs of one genome pair from its libraries alone, and each record's status in the other genome
(LTRmix's orthology).

Steps, for records of a and b (Records):
  candidates  record pairs sharing k-mers with at most CMAX_PAIR copies in each library, weighted
              (seq.kmers.cross_pairs)
  anchors     reciprocal best candidates by weight (seq.relatives.reciprocal_best, ties to the lower index) that pass
              the LTR clock (clock.clock_filter) or carry the same annotated TSD (canonical orientation): two
              independent certificates of one insertion event
  blocks      chromosome pairs enriched for anchors and their collinear blocks (correspondence)
  windows     each record's syntenic window in the other genome (absence)
  orthologs   candidates inside the window of a's record (the other record on the block's chromosome, its midpoint
              between lo and hi), or anchors inside a block, that also outweigh every candidate either record has
              outside its window (a same-family copy that merely sits in the window when the ortholog is absent does
              not beat the element's out-of-window paralogs, anchors in blocks exempt); of those the reciprocal best;
              each then passes the clock (tested here unless it was a reciprocal best pair) or shares its TSD
  status      2 present (an ortholog), 1 absent (inside a window, no ortholog), 0 unresolved (no window)
The divergence of an ortholog pair is its whole-record mismatches / aligned columns from the clock's alignment.
"""
from __future__ import annotations

import logging
import time
from dataclasses import dataclass

import numpy as np

from ltrtrace.base.constants import SCREEN_SEED
from ltrtrace.routes.library.absence import Windows, windows
from ltrtrace.routes.library.clock import clock_filter
from ltrtrace.routes.library.correspondence import Block, collinear_blocks
from ltrtrace.seq.kmers import PairWeights, cross_pairs
from ltrtrace.seq.relatives import reciprocal_best

log = logging.getLogger(__name__)


@dataclass(eq=False)
class PairOrtho:
    """The orthology of genomes a and b: per record of a, its ortholog in b (-1 none), its status in b (2 present, 1
    absent, 0 unresolved), its window, and with an ortholog the whole-record divergence (NaN otherwise, or when the
    pair did not align) and the aligned columns and mismatches behind it; the same per record of b; the blocks."""

    a: str
    b: str
    orth_b: np.ndarray
    orth_a: np.ndarray
    status_a: np.ndarray
    status_b: np.ndarray
    win_a: Windows
    win_b: Windows
    div_a: np.ndarray
    counts_a: np.ndarray  # int64 (n_a, 2): aligned columns and mismatches with the ortholog; 0 without one
    blocks: list[Block]

    @property
    def room_a(self) -> np.ndarray:
        return self.win_a.room

    @property
    def room_b(self) -> np.ndarray:
        return self.win_b.room


def _same_tsd(la, lb, ia, ib) -> np.ndarray:
    """True where both records have an annotated TSD and they are equal (canonical orientation)."""
    return (la.tsd[ia] != "") & (la.tsd[ia] == lb.tsd[ib])


def orthology(la, lb, *, threads: int, candidates: PairWeights | None = None) -> PairOrtho:
    """The 1:1 orthologs of records of la and lb (Records) and each record's status (module docstring), from the
    pair's candidate record pairs (`candidates`, as seq.kmers.cross_pairs gives them; computed here when None), with
    up to `threads` threads for the k-mer scan and the alignments."""
    t0 = time.perf_counter()
    a, b, w = cross_pairs(la.lib, lb.lib, threads=threads) if candidates is None else candidates
    log.debug("%s-%s: %d candidate record pairs (%.1f s)", la.label, lb.label, a.size, time.perf_counter() - t0)
    recip = reciprocal_best(a, b, w, la.n, lb.n)
    kclock, dr, cr = clock_filter(la, lb, a[recip], b[recip], threads)
    anc = recip[kclock | _same_tsd(la, lb, a[recip], b[recip])]
    clock_ok = dict(zip(recip.tolist(), kclock.tolist()))
    div_pair = dict(zip(recip.tolist(), dr.tolist()))
    counts_pair = dict(zip(recip.tolist(), cr[:, :2].tolist()))
    blocks, anchor_in_block, kept = collinear_blocks(la, lb, a, b, anc, np.random.default_rng(SCREEN_SEED))
    log.debug("%s-%s: %d reciprocal best pairs, %d pass the clock or share their TSD; %d enriched chromosome pairs; "
              "%d anchors in %d collinear blocks", la.label, lb.label, recip.size, anc.size, kept,
              int(anchor_in_block.sum()), len(blocks))
    win_a, win_b = windows(la, blocks, "a"), windows(lb, blocks, "b")
    chrom_b = np.array([blk.chrom_b for blk in blocks] + [None], object)
    ok = (win_a.bid[a] >= 0) & (lb.chrom[b] == chrom_b[win_a.bid[a]]) & (lb.mid[b] >= win_a.lo[a]) \
        & (lb.mid[b] <= win_a.hi[a])
    ok |= anchor_in_block
    off = ~ok
    w_off_a, w_off_b = np.zeros(la.n), np.zeros(lb.n)
    np.maximum.at(w_off_a, a[off], w[off])
    np.maximum.at(w_off_b, b[off], w[off])
    ok &= ((w > w_off_a[a]) & (w > w_off_b[b])) | anchor_in_block
    ci = np.flatnonzero(ok)
    mut = ci[reciprocal_best(a[ci], b[ci], w[ci], la.n, lb.n)]
    new = np.array([k for k in mut.tolist() if k not in clock_ok], np.int64)
    if new.size:
        kp, dn, cn = clock_filter(la, lb, a[new], b[new], threads)
        clock_ok.update(zip(new.tolist(), kp.tolist()))
        div_pair.update(zip(new.tolist(), dn.tolist()))
        counts_pair.update(zip(new.tolist(), cn[:, :2].tolist()))
    clock_pass = np.array([clock_ok[k] for k in mut.tolist()], bool)
    tsd_eq = _same_tsd(la, lb, a[mut], b[mut])
    tested = clock_pass | tsd_eq
    rescued, rejected = int((~clock_pass & tsd_eq).sum()), int((~tested).sum())
    mut = mut[tested]
    orth_b, orth_a = np.full(la.n, -1, np.int64), np.full(lb.n, -1, np.int64)
    orth_b[a[mut]] = b[mut]
    orth_a[b[mut]] = a[mut]
    div_a = np.full(la.n, np.nan)
    div_a[a[mut]] = [div_pair[k] for k in mut.tolist()]
    counts_a = np.zeros((la.n, 2), np.int64)
    counts_a[a[mut]] = np.array([counts_pair[k] for k in mut.tolist()], np.int64).reshape(-1, 2)
    status_a = np.where(orth_b >= 0, 2, np.where(win_a.bid >= 0, 1, 0)).astype(np.int8)
    status_b = np.where(orth_a >= 0, 2, np.where(win_b.bid >= 0, 1, 0)).astype(np.int8)
    log.info("%s-%s: %d 1:1 orthologs (%d rescued by an identical TSD after failing the clock), %d rejected; %s "
             "present/absent/unresolved %s (%.1f s)", la.label, lb.label, mut.size, rescued, rejected, la.label,
             np.bincount(status_a, minlength=3)[::-1].tolist(), time.perf_counter() - t0)
    return PairOrtho(la.label, lb.label, orth_b, orth_a, status_a, status_b, win_a, win_b, div_a, counts_a, blocks)
