"""The library route's alignment jobs: whole-record and LTR counts of record pairs and the LTR clock filter on them,
ortholog triples projected onto the donor's record, and each record's own-LTR check.

Records are N-stripped (Records.seq); the aligner is edlib's unit-cost global path (seq.align.unit_cost_path and
position_map). Two records of known strands are both in element orientation; when either strand is unknown the
second takes its closer strand by LTRmix's rule (seq.align.orient with mode "HW", head ORIENT_HEAD, slack
ORIENT_SLACK). Records whose lengths differ more than twofold (MIN_ALN_FRAC) are not aligned. Positions are 0-based
in the stripped records.

Scheduling: each job list in chunks of PAIR_CHUNK items (TRIPLE_CHUNK triples), at most POOL_JOBS per
thread, through base.runtime.pool_map's default: the first chunk in this process, the others in forked workers;
results in input order, whatever the threads.
"""
from __future__ import annotations

import logging
import time

import numpy as np

from synltr.base.constants import (
    CLOCK_LTR_FLOOR,
    CLOCK_LTR_SITES,
    MIN_ALN_FRAC,
    ORIENT_HEAD,
    ORIENT_SLACK,
    PAIR_CHUNK,
    POOL_JOBS,
    SELF_Z,
    TRIPLE_CHUNK,
)
from synltr.base.runtime import pool_map, shared
from synltr.io.seqio import revcomp_codes
from synltr.seq.align import orient, position_map, unit_cost_path
from synltr.seq.divergence import k2p

log = logging.getLogger(__name__)


def _oriented(x, y, sx, sy):
    """y in x's orientation: as given when both strands are known, else the closer strand (LTRmix's rule)."""
    if sx != "." and sy != ".":
        return y
    return y if orient(x, y, mode="HW", head=ORIENT_HEAD, slack=ORIENT_SLACK) > 0 else revcomp_codes(y)


def _short(*lengths) -> bool:
    """True when a sequence is empty or the lengths differ more than 1 / MIN_ALN_FRAC fold."""
    return min(lengths) == 0 or min(lengths) / max(lengths) < MIN_ALN_FRAC


def _chunks(items, size, threads):
    return np.array_split(items, max(1, min(len(items) // size, POOL_JOBS * max(threads, 1))))


def _pooled(job, chunks, threads, width, what, **state):
    """pool_map(job, chunks) with `state` in base.runtime.shared for the workers; the rows stacked (`width` columns);
    the pool's wall time logged."""
    t0 = time.perf_counter()
    shared.update(state)
    try:
        out = np.concatenate(pool_map(job, chunks, threads)).reshape(-1, width)
    finally:
        for k in state:
            shared.pop(k, None)
    log.debug("%s: %d item(s) in %d chunk(s) aligned in %.1f s", what, len(out), len(chunks),
              time.perf_counter() - t0)
    return out


# ---------------------------------------------------------------- record pairs and the clock


def _pair_job(chunk):
    la, lb = shared["screen_a"], shared["screen_b"]
    out = np.zeros((len(chunk), 4), np.int64)
    for n, (i, j) in enumerate(chunk):
        x, y = la.seq(i), lb.seq(j)
        if _short(x.size, y.size):
            continue
        y = _oriented(x, y, la.strand[i], lb.strand[j])
        m = position_map(unit_cost_path(x, y))
        ok = m >= 0
        diff = np.zeros(x.size, bool)
        diff[ok] = x[ok] != y[m[ok]]
        out[n, 0], out[n, 1] = int(ok.sum()), int(diff.sum())
        ltrs = la.ltrs(i)
        if ltrs is not None:
            lt = np.zeros(x.size, bool)
            for a0, a1 in ltrs:
                lt[a0:a1] = True
            out[n, 2], out[n, 3] = int((ok & lt).sum()), int((diff & lt).sum())
    return out


def pair_counts(la, lb, ia, ib, threads: int) -> np.ndarray:
    """For record pairs (ia[k] of la, ib[k] of lb; Records): aligned columns and mismatches over the two whole records,
    then the same over la's record's two LTRs (0 without its bounds), int64 (n, 4); a pair not aligned (empty or
    lengths more than twofold apart) is all 0. Up to `threads` processes."""
    if len(ia) == 0:
        return np.zeros((0, 4), np.int64)
    pairs = np.stack([np.asarray(ia, np.int64), np.asarray(ib, np.int64)], 1)
    return _pooled(_pair_job, _chunks(pairs, PAIR_CHUNK, threads), threads, 4, f"{la.label}-{lb.label} pairs",
                   screen_a=la, screen_b=lb)


def clock_filter(la, lb, ai, bi, threads: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The clock test of record pairs (ai[k] of la, bi[k] of lb): orthologous copies coalesced after the insertion, so
    their divergence cannot exceed the element's own LTR-LTR divergence beyond sampling error; paralogs have no such
    bound. Divergence: LTR with LTR when >= CLOCK_LTR_SITES LTR columns align, else whole records (mismatches /
    columns); bound: the larger annotated K2P of the two, with se = sqrt(d / columns + K2P / max(shorter LTR,
    CLOCK_LTR_FLOOR)). A pair passes when it aligned and its divergence is <= K2P + SELF_Z se, or at most one
    mismatch, or no K2P is known. Returns (passes, whole-record divergence (NaN unaligned), pair_counts' rows)."""
    c = pair_counts(la, lb, ai, bi, threads)
    sites, mism = c[:, 0].astype(float), c[:, 1].astype(float)
    lsites, lmism = c[:, 2].astype(float), c[:, 3].astype(float)
    d = np.where(sites > 0, mism / np.maximum(sites, 1), np.nan)
    use_ltr = lsites >= CLOCK_LTR_SITES
    dt = np.where(use_ltr, lmism / np.maximum(lsites, 1), d)
    nt = np.where(use_ltr, lsites, sites)
    mt = np.where(use_ltr, lmism, mism)
    kk = np.fmax(la.k2p[ai], lb.k2p[bi])
    lt = np.fmax(np.minimum(la.ltr_len[ai, 0], la.ltr_len[ai, 1]), CLOCK_LTR_FLOOR).astype(float)
    se = np.sqrt(np.nan_to_num(dt) / np.maximum(nt, 1) + np.nan_to_num(kk) / lt)
    ok = (sites > 0) & (~np.isfinite(kk) | (dt - kk <= SELF_Z * se) | (mt <= 1))
    return ok, d, c


# ---------------------------------------------------------------- ortholog triples and the own-LTR check


def _triple_job(chunk):
    ld, lr, ls, self_ok = (shared[k] for k in ("screen_d", "screen_r", "screen_s", "screen_self_ok"))
    out = np.zeros((len(chunk), 7), np.int64)
    for n, (k, i, j) in enumerate(chunk):
        x, y, z = ld.seq(k), lr.seq(i), ls.seq(j)
        if _short(x.size, y.size, z.size):
            continue
        y = _oriented(x, y, ld.strand[k], lr.strand[i])
        z = _oriented(x, z, ld.strand[k], ls.strand[j])
        my = position_map(unit_cost_path(x, y))
        mz = position_map(unit_cost_path(x, z))
        ok = (my >= 0) & (mz >= 0)
        if not ok.any():
            continue
        xs, ys, zs = x[ok], y[my[ok]], z[mz[ok]]
        out[n, :4] = (int(ok.sum()), int((ys != xs).sum()), int((zs != xs).sum()), int((ys != zs).sum()))
        ltrs = ld.ltrs(k)
        if not self_ok[k] or ltrs is None:
            continue
        (a0, a1), (b0, b1) = ltrs
        fwd = position_map(unit_cost_path(x[a0:a1], x[b0:b1]))
        rev = position_map(unit_cost_path(x[b0:b1], x[a0:a1]))
        for p0, p1, q0, mp in ((a0, a1, b0, fwd), (b0, b1, a0, rev)):
            pos = np.arange(p0, p1)
            q = mp[pos - p0]
            g = q >= 0
            pos, q = pos[g], q[g] + q0
            if pos.size == 0:
                continue
            v = (my[pos] >= 0) & (mz[pos] >= 0) & (my[q] >= 0) & (mz[q] >= 0)
            pos, q = pos[v], q[v]
            d_, r_, s_ = x[pos], y[my[pos]], z[mz[pos]]
            o1, o2, o3 = x[q], y[my[q]], z[mz[q]]
            og = (o1 == o2) & (o2 == o3)
            d_, r_, s_, o = d_[og], r_[og], s_[og], o1[og]
            out[n, 4] += int(((r_ == d_) & (s_ != d_) & (s_ == o)).sum())
            out[n, 5] += int(((s_ == d_) & (r_ != d_) & (r_ == o)).sum())
            out[n, 6] += int(((r_ == s_) & (d_ != r_) & (d_ == o)).sum())
    return out


def triple_counts(ld, lr, ls, trip, self_ok, threads: int) -> np.ndarray:
    """For ortholog triples (rows donor d of ld, recipient r of lr, sister s of ls; Records), both sisters projected
    onto the donor's record (alignment noise then touches r and s alike), int64 (n, 7): columns aligned in all three,
    r-d, s-d and r-s mismatches; then, where the donor's own LTRs pass the self check (self_ok, by donor record),
    ABBA, BABA and BBAA at LTR sites polarized by the donor element's other LTR (the outgroup base: the other-LTR base,
    equal in all three). Up to `threads` processes."""
    trip = np.asarray(trip, np.int64).reshape(-1, 3)
    if len(trip) == 0:
        return np.zeros((0, 7), np.int64)
    return _pooled(_triple_job, _chunks(trip, TRIPLE_CHUNK, threads), threads, 7,
                   f"{lr.label}-{ls.label}-{ld.label} triples", screen_d=ld, screen_r=lr, screen_s=ls,
                   screen_self_ok=np.asarray(self_ok))


def _self_job(chunk):
    lib = shared["screen_self"]
    out = np.zeros(len(chunk), np.int8)
    for n, i in enumerate(chunk):
        ltrs = lib.ltrs(i)
        if ltrs is None or not np.isfinite(lib.k2p[i]):
            continue
        x = lib.seq(i)
        (a0, a1), (b0, b1) = ltrs
        mp = position_map(unit_cost_path(x[a0:a1], x[b0:b1]))
        ok = mp >= 0
        p, q = x[a0:a1][ok], x[b0:b1][mp[ok]]
        diff = p != q
        ts = int(((p ^ q) == 2)[diff].sum())  # A<->G and C<->T: codes 0-2 and 1-3
        d, se = (float(v) for v in k2p(int(ok.sum()), ts, int(diff.sum()) - ts))
        if np.isfinite(d) and abs(d - lib.k2p[i]) <= max(SELF_Z * se, 1.0 / max(int(ok.sum()), 1)):
            out[n] = 1
    return out


def self_ltr_ok(lib, idx, threads: int) -> np.ndarray:
    """1 for each record of `idx` (Records) whose own two LTRs (at its record ends) give a K2P within SELF_Z standard
    errors, or one substitution, of its annotated K2P: usable as its own outgroup; 0 elsewhere, int8 per record of
    lib. Up to `threads` processes."""
    ok = np.zeros(lib.n, np.int8)
    idx = np.asarray(idx, np.int64)
    if idx.size:
        ok[idx] = _pooled(_self_job, _chunks(idx, PAIR_CHUNK, threads), threads, 1, f"{lib.label} own LTRs",
                          screen_self=lib).ravel()
    return ok
