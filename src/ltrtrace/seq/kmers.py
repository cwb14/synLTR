"""k-mer lineages of an LTR-RT library and record sketches.

Records are uint8 codes (A C G T = 0-3, any other base 4), concatenated, with record boundaries `offs` (the `codes` and
`offs` of io.library's Library). A k-mer is canonical as the smaller of its forward and reverse-complement 2-bit values;
a k-mer holding another base is skipped. `mix64` (MurmurHash3's fmix64) is the package's one 64-bit mixer: of canonical
k-mers, of carrier-set signatures and of sketch values.

`clades`: every canonical k-mer of the library hashed exactly (or FracMinHash-sampled: hashes <= 2^64 / scale), in
bounded passes over ranges of the hashes' top 16 bits; k-mers carried by 2..cap distinct copies on >= min_seq
sequences define lineages, kept as their distinct carrier sets (`Clades`). Results depend on the emission order, which
is fixed: passes in ascending hash ranges; within a pass, k-mers in hash order (buckets of the top 16 bits); a carrier
set is kept at its first occurrence, its multiplicity counting the k-mers that repeat it. Carrier sets are signed by
the xor of per-copy keys (drawn once from CLADE_KEYS_SEED) and mix64(size). With `cache`, the clades are a checkpoint
stage (base.checkpoint): <cache>/clades.npz and its done marker, keyed on a digest of what they are built from (record
boundaries, each copy's sequence, the codes), the parameters and this module's code.

`cross_pairs` and `single_copy_tallies`: canonical k-mers shared across libraries, counted exactly through the same
bounded passes run over two or more libraries at once (copies numbered through the libraries in turn): the record pairs
of two libraries sharing low-copy k-mers, with weights (the library screen's orthology candidates), and per record of
one library its single-copy k-mers by the other libraries that hold them. `shared_kmers` gives both for every pair and
every trio of a set of libraries from one scan whose passes each hash the libraries once (sized by one counting scan;
buckets sorted in parallel; the tallies counted by groups of trios times chunks of the buckets, so every thread works
at any number of libraries). mix64 is a bijection, so a run of one hash is one k-mer; the results do not depend on the
passes, the scan order or the threads.

`sketch`: bottom-s MinHash sketch of each record's canonical k-mers.
"""
from __future__ import annotations

import hashlib
import logging
import math
import zipfile
from itertools import combinations
from pathlib import Path
from typing import NamedTuple

import numba
import numpy as np
from numba import njit, prange

from ltrtrace.base.checkpoint import atomic_path, clear_done, code_fingerprint, is_done, mark_done, marker, stage_key
from ltrtrace.base.constants import (
    CAP,
    CLADE_KEYS_SEED,
    CMAX_PAIR,
    KMER_BUDGET,
    KMER_BYTES,
    MEM_GB,
    PAIR_CHUNKS,
    RUN_LIM,
    SKETCH_K,
    SKETCH_S,
    TALLY_CHUNK_BYTES,
    K,
)
from ltrtrace.base.runtime import set_numba_threads

log = logging.getLogger(__name__)


@njit(cache=True, inline="always")
def mix64(x):
    """MurmurHash3's 64-bit finalizer of a uint64, or elementwise of a uint64 array."""
    x = x ^ (x >> np.uint64(33))
    x = x * np.uint64(0xFF51AFD7ED558CCD)
    x = x ^ (x >> np.uint64(33))
    x = x * np.uint64(0xC4CEB9FE1A85EC53)
    return x ^ (x >> np.uint64(33))


@njit(parallel=True, cache=True)
def _kmer_pass(codes, offs, k, blo, bhi, cut, fill, cnt_or_off, H, E):
    """Canonical k-mer hashes of every copy whose top 16 bits lie in [blo, bhi) and that are <= cut (FracMinHash;
    cut = 2^64-1 keeps all). fill=False: count per copy into cnt_or_off; fill=True: write (hash, copy) from
    cnt_or_off[copy] on."""
    n = offs.size - 1
    mask = (np.uint64(1) << np.uint64(2 * k)) - np.uint64(1)
    sh = np.uint64(2 * (k - 1))
    for i in prange(n):
        f = np.uint64(0)
        r = np.uint64(0)
        run = 0
        o = cnt_or_off[i] if fill else 0
        cnt = 0
        for p in range(offs[i], offs[i + 1]):
            c = codes[p]
            if c > 3:
                run = 0
                f = np.uint64(0)
                r = np.uint64(0)
                continue
            cu = np.uint64(c)
            f = ((f << np.uint64(2)) | cu) & mask
            r = (r >> np.uint64(2)) | ((np.uint64(3) - cu) << sh)
            run += 1
            if run >= k:
                h = mix64(min(f, r))
                b = np.int64(h >> np.uint64(48))
                if h <= cut and b >= blo and b < bhi:
                    if fill:
                        H[o] = h
                        E[o] = i
                        o += 1
                    else:
                        cnt += 1
        if not fill:
            cnt_or_off[i] = cnt


@njit(cache=True)
def _bucket_sort(H, E, blo, nb):
    """(hash, copy) records moved into top-16-bit bucket order (counting sort), and the bucket bounds: each bucket
    is then a contiguous slice that sorts in cache."""
    bounds = np.zeros(nb + 1, np.int64)
    for j in range(H.size):
        bounds[np.int64(H[j] >> np.uint64(48)) - blo + 1] += 1
    for b in range(nb):
        bounds[b + 1] += bounds[b]
    pos = bounds[:-1].copy()
    Hs = np.empty_like(H)
    Es = np.empty_like(E)
    for j in range(H.size):
        b = np.int64(H[j] >> np.uint64(48)) - blo
        Hs[pos[b]] = H[j]
        Es[pos[b]] = E[j]
        pos[b] += 1
    return Hs, Es, bounds


@njit(cache=True)
def _distinct(vals, m, out):
    """Sorted distinct values of vals[:m] written to out; returns their number. In place, no allocation: runs are
    short, so insertion sort (array sort above 32 values)."""
    for i in range(m):
        out[i] = vals[i]
    if m > 32:
        out[:m].sort()
    else:
        for i in range(1, m):
            x = out[i]
            j = i - 1
            while j >= 0 and out[j] > x:
                out[j + 1] = out[j]
                j -= 1
            out[j + 1] = x
    nd = 0
    for t in range(m):
        if t == 0 or out[t] != out[t - 1]:
            out[nd] = out[t]
            nd += 1
    return nd


@njit(parallel=True, cache=True)
def _scan_buckets(H, E, bounds, seqof, rnd, cap, min_seq, run_lim, write, nq, nm, qoff, moff, sig, size, nseq, mem,
                  nchunk):
    """Per bucket: sort by hash; runs of one k-mer (2..run_lim occurrences) carried by 2..cap distinct copies on
    >= min_seq sequences qualify.
    write=False: count qualifying runs (nq) and members (nm) per bucket; write=True: write signature (xor of per-copy
    random keys ^ mix64(size)), size, n sequences and members at the bucket offsets. Buckets are processed in nchunk
    parallel chunks, each with its own scratch buffers (allocated once per chunk, not per bucket)."""
    nb = bounds.size - 1
    for ch in prange(nchunk):
        cps = np.empty(run_lim, np.int64)
        buf = np.empty(run_lim, np.int64)
        sq = np.empty(run_lim, np.int64)
        sbuf = np.empty(run_lim, np.int64)
        for b in range(nb * ch // nchunk, nb * (ch + 1) // nchunk):
            s, e = bounds[b], bounds[b + 1]
            q = 0
            mm = 0
            if e - s >= 2:
                seg = np.argsort(H[s:e]) + s          # bucket is contiguous (bucket order)
                qo = qoff[b] if write else 0
                mo = moff[b] if write else 0
                L = e - s
                j = 0
                while j < L:
                    t = j + 1
                    hj = H[seg[j]]
                    while t < L and H[seg[t]] == hj:
                        t += 1
                    m = t - j
                    if 2 <= m <= run_lim:
                        for u in range(m):
                            cps[u] = E[seg[j + u]]
                        nd = _distinct(cps, m, buf)
                        if 2 <= nd <= cap:
                            for u in range(nd):
                                sq[u] = seqof[buf[u]]
                            ns = _distinct(sq, nd, sbuf)
                            if ns >= min_seq:
                                if write:
                                    sg = mix64(np.uint64(nd))
                                    for u in range(nd):
                                        sg = sg ^ rnd[buf[u]]
                                    sig[qo] = sg
                                    size[qo] = nd
                                    nseq[qo] = ns
                                    for u in range(nd):
                                        mem[mo + u] = buf[u]
                                    qo += 1
                                    mo += nd
                                else:
                                    q += 1
                                    mm += nd
                    j = t
            if not write:
                nq[b] = q
                nm[b] = mm


class Clades:
    """Distinct carrier sets of shared k-mers (2..cap copies on >= min_seq sequences). CSR: members of clade i are
    members[indptr[i]:indptr[i+1]]; mult = k-mers with exactly that carrier set; nseq = sequences spanned; n_occ =
    k-mer occurrences hashed, n_multi = qualifying k-mers; scale = the FracMinHash scale (1 = exact)."""

    def __init__(self, members, size, mult, nseq, n_occ=0, n_multi=0, scale=1):
        self.members = members.astype(np.int32)
        self.size = size.astype(np.int32)
        self.mult = mult.astype(np.int64)
        self.nseq = nseq.astype(np.int32)
        self.indptr = np.concatenate(([0], np.cumsum(self.size, dtype=np.int64)))
        self.n = int(self.size.size)
        self.n_occ, self.n_multi = int(n_occ), int(n_multi)
        self.scale = scale

    def sid(self):
        return np.repeat(np.arange(self.n), self.size)

    def restrict(self, copy_mask, seq=None):
        """Members inside copy_mask; clades keep >= 2 members (on >= 2 sequences when seq is given). Clades whose
        kept carrier sets coincide merge (mult summed): a lineage reaching removed copies is not counted twice."""
        keep_m = copy_mask[self.members]
        sid = self.sid()[keep_m]
        mem = self.members[keep_m]
        sz = np.bincount(sid, minlength=self.n)
        ok = sz >= 2
        if seq is not None and mem.size:
            u = np.unique(sid.astype(np.int64) * (int(seq.max()) + 1) + seq[mem])
            ns = np.bincount(u // (int(seq.max()) + 1), minlength=self.n)
            ok &= ns >= 2
        else:
            ns = self.nseq.copy()
        km = ok[sid]
        mem, size, mult, nseq = mem[km], sz[ok], self.mult[ok], (ns if seq is not None else self.nseq)[ok]
        if size.size and not copy_mask.all():                    # sets differing only by removed copies coincide now
            _, size, nseq, mult, mem = _dedupe(_set_signatures(mem, size), size, nseq, mult, mem)
        return Clades(mem, size, mult, nseq, self.n_occ, self.n_multi, self.scale)


def _set_signatures(mem, size):
    """Order-free signature of each carrier set (xor of per-copy keys ^ mix64(size)), clades concatenated in mem."""
    keys = mix64(mem.astype(np.uint64) * np.uint64(0x9E3779B97F4A7C15) + np.uint64(1))
    starts = np.concatenate(([0], np.cumsum(size)[:-1])).astype(np.int64)
    return np.bitwise_xor.reduceat(keys, starts) ^ mix64(np.asarray(size, np.uint64))


def _dedupe(sig, size, nseq, mult, mem):
    """Merge entries with the same carrier-set signature (mult summed). Entries keep the order of their first
    occurrence; mem (members concatenated in entry order) keeps the members of each first occurrence."""
    u, first, inv = np.unique(sig, return_index=True, return_inverse=True)
    mult_u = np.bincount(inv, weights=mult, minlength=u.size).astype(np.int64)
    o = np.argsort(first)
    keep = np.zeros(sig.size, bool)
    keep[first] = True
    return u[o], size[first[o]], nseq[first[o]], mult_u[o], mem[np.repeat(keep, size)]


def auto_scale(n_bp):
    """FracMinHash scale keeping at most KMER_BUDGET k-mers (1 = exact)."""
    return max(1, math.ceil(n_bp / KMER_BUDGET))


def _cut(scale):
    """The largest hash kept at FracMinHash `scale` (1 = every hash)."""
    return np.uint64(2 ** 64 - 1) if scale <= 1 else np.uint64((1 << 64) // scale - 1)


def _passes(libs, k, cut, mem_gb):
    """(total k-mers, hash-range edges): the k-mers of all `libs` counted once, then split into ranges of the hashes'
    top 16 bits holding at most mem_gb of (hash, copy) records each."""
    eH, eE = np.empty(0, np.uint64), np.empty(0, np.int32)
    total = 0
    for lib in libs:
        counts = np.zeros(len(lib.offs) - 1, np.int64)
        _kmer_pass(lib.codes, lib.offs, k, 0, 65536, cut, False, counts, eH, eE)
        total += int(counts.sum())
    per_pass = max(1, int(mem_gb * 1e9 / KMER_BYTES))
    n_pass = min(65536, max(1, math.ceil(total / per_pass)))
    return total, np.linspace(0, 65536, n_pass + 1).astype(np.int64)


def _pass_records(libs, k, blo, bhi, cut):
    """The (hash, copy) records of the libraries' k-mers whose hash's top 16 bits lie in [blo, bhi) and that are <=
    cut, in bucket order (`_bucket_sort`), with the bucket bounds: (H, E, bounds); None when there are none. Copies
    are numbered through the libraries in turn (the second library's first copy follows the first's last)."""
    eH, eE = np.empty(0, np.uint64), np.empty(0, np.int32)
    cnts = []
    for lib in libs:
        cnt = np.zeros(len(lib.offs) - 1, np.int64)
        _kmer_pass(lib.codes, lib.offs, k, blo, bhi, cut, False, cnt, eH, eE)
        cnts.append(cnt)
    tot = sum(int(c.sum()) for c in cnts)
    if tot == 0:
        return None
    H, E = np.empty(tot, np.uint64), np.empty(tot, np.int32)
    start = first = 0
    for lib, cnt in zip(libs, cnts):
        off = np.zeros(cnt.size, np.int64)
        off[1:] = np.cumsum(cnt)[:-1]
        off += start
        _kmer_pass(lib.codes, lib.offs, k, blo, bhi, cut, True, off, H, E)
        n = int(cnt.sum())
        if first:
            E[start:start + n] += np.int32(first)
        start += n
        first += cnt.size
    return _bucket_sort(H, E, blo, bhi - blo)


def _build(lib, seqof, *, k, cap, min_seq, scale, mem_gb):
    """The clades of the library (module docstring), in passes holding at most mem_gb of (hash, copy) records."""
    n = len(lib.offs) - 1
    rng = np.random.default_rng(CLADE_KEYS_SEED)
    rnd = (rng.integers(1, 2 ** 62, n, dtype=np.int64).astype(np.uint64) << np.uint64(1)) | np.uint64(1)
    cut = _cut(scale)
    total, edges = _passes([lib], k, cut, mem_gb)
    P_sig, P_size, P_nseq, P_mult, P_mem = [], [], [], [], []
    n_multi = 0
    for p in range(edges.size - 1):
        blo, bhi = int(edges[p]), int(edges[p + 1])
        if bhi <= blo:
            continue
        rec = _pass_records([lib], k, blo, bhi, cut)
        if rec is None:
            continue
        H, E, bounds = rec
        nb = bhi - blo
        nq, nm = np.zeros(nb, np.int64), np.zeros(nb, np.int64)
        z64, zu, z32 = np.zeros(1, np.int64), np.zeros(1, np.uint64), np.zeros(1, np.int32)
        nchunk = max(1, min(nb, 8 * numba.get_num_threads()))
        _scan_buckets(H, E, bounds, seqof, rnd, cap, min_seq, RUN_LIM, False, nq, nm, z64, z64, zu, z32, z32, z64,
                      nchunk)
        qoff = np.concatenate(([0], np.cumsum(nq)[:-1]))
        moff = np.concatenate(([0], np.cumsum(nm)[:-1]))
        Q, M = int(nq.sum()), int(nm.sum())
        n_multi += Q
        if Q == 0:
            continue
        sig, size, nseq, mem = np.empty(Q, np.uint64), np.empty(Q, np.int32), np.empty(Q, np.int32), np.empty(M, np.int64)
        _scan_buckets(H, E, bounds, seqof, rnd, cap, min_seq, RUN_LIM, True, nq, nm, qoff, moff, sig, size, nseq, mem,
                      nchunk)
        del H, E
        sig, size, nseq, mult, mem = _dedupe(sig, size, nseq, np.ones(Q, np.int64), mem)
        P_sig.append(sig)
        P_size.append(size)
        P_nseq.append(nseq)
        P_mult.append(mult)
        P_mem.append(mem)
    if not P_sig:
        return Clades(np.empty(0, np.int32), np.empty(0, np.int32), np.empty(0, np.int64), np.empty(0, np.int32), total, 0)
    _, size, nseq, mult, members = _dedupe(*[np.concatenate(x) for x in (P_sig, P_size, P_nseq, P_mult, P_mem)])
    return Clades(members, size, mult, nseq, total, n_multi)


# ---- k-mers shared across libraries (exact counts): candidate record pairs and single-copy tallies ----


class PairWeights(NamedTuple):
    """Record pairs of two libraries sharing low-copy k-mers, sorted by (a, b): record indices a and b (int64) and
    w, the summed weight of their shared k-mers (float64)."""

    a: np.ndarray
    b: np.ndarray
    w: np.ndarray


@njit(parallel=True, cache=True)
def _sort_buckets(H, E, bounds, nchunk):
    """Each bucket's (hash, copy) records sorted by hash, in place (copies moved alike)."""
    nb = bounds.size - 1
    for ch in prange(nchunk):
        for b in range(nb * ch // nchunk, nb * (ch + 1) // nchunk):
            s, e = bounds[b], bounds[b + 1]
            if e - s >= 2:
                o = np.argsort(H[s:e])
                h, c = H[s:e][o], E[s:e][o]
                H[s:e] = h
                E[s:e] = c


@njit(parallel=True, cache=True)
def _pair_runs(H, E, bounds, n_a, cmax, write, cnt, B, W, nchunk):
    """Over buckets sorted by hash (_sort_buckets), runs of one hash (one k-mer): with ca occurrences in library a
    (copies < n_a) and cb in library b, both 1..cmax, each of the ca x cb occurrence pairs gets the weight 1/(ca cb)
    rounded to float32. write=False: count them per (chunk, a copy) into cnt; write=True: write b's copy and the
    weight at cnt[chunk, a copy], the next free slot of a's block (a-major, chunks in turn), advancing it."""
    nb = bounds.size - 1
    for ch in prange(nchunk):
        for b in range(nb * ch // nchunk, nb * (ch + 1) // nchunk):
            s, e = bounds[b], bounds[b + 1]
            j = s
            while j < e:
                t = j + 1
                while t < e and H[t] == H[j]:
                    t += 1
                ca = 0
                for u in range(j, t):
                    if E[u] < n_a:
                        ca += 1
                cb = t - j - ca
                if 1 <= ca <= cmax and 1 <= cb <= cmax:
                    w = np.float32(1.0 / (ca * cb))
                    for u in range(j, t):
                        x = E[u]
                        if x < n_a:
                            if write:
                                for v in range(j, t):
                                    if E[v] >= n_a:
                                        o = cnt[ch, x]
                                        B[o] = E[v] - n_a
                                        W[o] = w
                                        cnt[ch, x] = o + 1
                            else:
                                cnt[ch, x] += cb
                j = t


@njit(parallel=True, cache=True)
def _pair_sums(B, W, start, write, nd, off, UB, UW):
    """Per a copy (its block B[start[a]:start[a + 1]], W alike): the distinct b copies, ascending, and the float64 sums
    of their weights. write=False: sort the block in place and count them into nd; write=True: write them from
    off[a]."""
    for a in prange(start.size - 1):
        s, e = start[a], start[a + 1]
        if e == s:
            continue
        if not write:
            o = np.argsort(B[s:e])
            bb, ww = B[s:e][o], W[s:e][o]
            B[s:e] = bb
            W[s:e] = ww
            m = 1
            for i in range(s + 1, e):
                if B[i] != B[i - 1]:
                    m += 1
            nd[a] = m
        else:
            q = off[a] - 1
            for i in range(s, e):
                if i == s or B[i] != B[i - 1]:
                    q += 1
                    UB[q] = B[i]
                    UW[q] = 0.0
                UW[q] += W[i]


@njit(cache=True)
def _merge_sum(u1, s1, u2, s2):
    """Two sorted distinct-key arrays with their sums, merged: keys in both get the sum of both."""
    u = np.empty(u1.size + u2.size, np.uint64)
    sw = np.empty(u1.size + u2.size, np.float64)
    i = j = m = 0
    while i < u1.size or j < u2.size:
        if j == u2.size or (i < u1.size and u1[i] < u2[j]):
            u[m], sw[m] = u1[i], s1[i]
            i += 1
        elif i == u1.size or u2[j] < u1[i]:
            u[m], sw[m] = u2[j], s2[j]
            j += 1
        else:
            u[m], sw[m] = u1[i], s1[i] + s2[j]
            i += 1
            j += 1
        m += 1
    return u[:m].copy(), sw[:m].copy()


def _exact(scale):
    if scale != 1:
        raise ValueError(f"cross-library k-mer counts are exact only (scale 1), not scale {scale}")


def _cross_pairs(lib_a, lib_b, *, k, cmax, mem_gb):
    n_a, n_b = len(lib_a.offs) - 1, len(lib_b.offs) - 1
    if n_a + n_b >= 2 ** 31:
        raise ValueError("more than 2^31 - 1 records in the two libraries")
    cut = _cut(1)
    total, edges = _passes([lib_a, lib_b], k, cut, mem_gb)
    log.debug("cross pairs: %d k-mers in %d pass(es)", total, edges.size - 1)
    U, SW = np.empty(0, np.uint64), np.empty(0, np.float64)
    nchunk = numba.get_num_threads()
    for p in range(edges.size - 1):
        blo, bhi = int(edges[p]), int(edges[p + 1])
        rec = None if bhi <= blo else _pass_records([lib_a, lib_b], k, blo, bhi, cut)
        if rec is None:
            continue
        H, E, bounds = rec
        _sort_buckets(H, E, bounds, nchunk)
        cnt = np.zeros((nchunk, n_a), np.int64)
        e32, ef = np.zeros(1, np.int32), np.zeros(1, np.float32)
        _pair_runs(H, E, bounds, n_a, cmax, False, cnt, e32, ef, nchunk)
        per_a = cnt.sum(0)
        start = np.concatenate(([0], np.cumsum(per_a))).astype(np.int64)
        cnt = start[:-1] + np.cumsum(cnt, 0) - cnt                    # each chunk's first slot in each a block
        B, W = np.empty(int(start[-1]), np.int32), np.empty(int(start[-1]), np.float32)
        _pair_runs(H, E, bounds, n_a, cmax, True, cnt, B, W, nchunk)
        del H, E, cnt
        nd = np.zeros(n_a, np.int64)
        _pair_sums(B, W, start, False, nd, nd, B, np.zeros(1, np.float64))
        off = np.concatenate(([0], np.cumsum(nd)[:-1])).astype(np.int64)
        UB, UW = np.empty(int(nd.sum()), np.int32), np.empty(int(nd.sum()), np.float64)
        _pair_sums(B, W, start, True, nd, off, UB, UW)
        del B, W
        pk = (np.repeat(np.arange(n_a, dtype=np.uint64), nd) << np.uint64(32)) | UB.astype(np.uint64)
        U, SW = _merge_sum(U, SW, pk, UW)
    return PairWeights((U >> np.uint64(32)).astype(np.int64), (U & np.uint64(0xFFFFFFFF)).astype(np.int64), SW)


def cross_pairs(lib_a, lib_b, *, k=K, cmax=CMAX_PAIR, scale=1, mem_gb=MEM_GB, threads=1) -> PairWeights:
    """Record pairs (a of lib_a, b of lib_b) sharing canonical k-mers that occur 1..cmax times in each library, with
    weight w = the sum over those k-mers of 1/(ca cb) per occurrence pair (ca, cb: the k-mer's occurrences in each
    library, a record's repeats included; each term rounded to float32, the sum in float64, which is exact for these
    terms, so the order of the scan does not matter). Libraries are codes and offs; bounded passes of at most mem_gb
    of (hash, copy) records, numba loops on `threads` threads; exact counts only (`scale` 1, else ValueError)."""
    _exact(scale)
    return _with_threads(threads, _cross_pairs, lib_a, lib_b, k=k, cmax=cmax, mem_gb=mem_gb)


@njit(parallel=True, cache=True)
def _tally_runs(H, E, bounds, lo, focal, out, nchunk):
    """Over buckets sorted by hash (_sort_buckets), runs of one hash: a k-mer occurring once in the focal library and at
    most once in each other library (library l holds copies lo[l]..lo[l + 1] - 1) adds 1 to out[chunk, focal copy,
    pattern], pattern = the bits (in library order, the focal one skipped) of the other libraries holding it."""
    nb = bounds.size - 1
    nlib = lo.size - 1
    c = np.zeros((nchunk, nlib), np.int64)
    for ch in prange(nchunk):
        for b in range(nb * ch // nchunk, nb * (ch + 1) // nchunk):
            s, e = bounds[b], bounds[b + 1]
            j = s
            while j < e:
                t = j + 1
                while t < e and H[t] == H[j]:
                    t += 1
                for q in range(nlib):
                    c[ch, q] = 0
                x = -1
                for u in range(j, t):
                    y = E[u]
                    q = 0
                    while y >= lo[q + 1]:
                        q += 1
                    c[ch, q] += 1
                    if q == focal:
                        x = y - lo[focal]
                ok = c[ch, focal] == 1
                pat = 0
                bit = 0
                for q in range(nlib):
                    if q == focal:
                        continue
                    if c[ch, q] > 1:
                        ok = False
                    elif c[ch, q] == 1:
                        pat |= 1 << bit
                    bit += 1
                if ok:
                    out[ch, x, pat] += 1
                j = t


def _tallies(libs, focal, k, mem_gb):
    lo = np.concatenate(([0], np.cumsum([len(lib.offs) - 1 for lib in libs]))).astype(np.int64)
    if lo[-1] >= 2 ** 31:
        raise ValueError("more than 2^31 - 1 records in the libraries")
    cut = _cut(1)
    total, edges = _passes(libs, k, cut, mem_gb)
    log.debug("single-copy tallies: %d k-mers in %d pass(es)", total, edges.size - 1)
    n = int(lo[focal + 1] - lo[focal])
    total = np.zeros((n, 1 << (len(libs) - 1)), np.int64)
    for p in range(edges.size - 1):
        blo, bhi = int(edges[p]), int(edges[p + 1])
        rec = None if bhi <= blo else _pass_records(libs, k, blo, bhi, cut)
        if rec is None:
            continue
        H, E, bounds = rec
        nchunk = max(1, min(bhi - blo, numba.get_num_threads()))
        _sort_buckets(H, E, bounds, nchunk)
        out = np.zeros((nchunk, n, total.shape[1]), np.int32)
        _tally_runs(H, E, bounds, lo, focal, out, nchunk)
        total += out.sum(0, dtype=np.int64)
    return total


def single_copy_tallies(libs, *, focal, k=K, scale=1, mem_gb=MEM_GB, threads=1) -> dict[tuple[str, ...], np.ndarray]:
    """Per record of the focal library (a key of `libs`, label -> library of codes and offs): its canonical k-mers
    that occur once in it and at most once in each other library, counted by the set of other libraries holding them.
    Returns {tuple of those libraries' labels, in `libs` order: int64 counts per focal record}, every subset of the
    others from none up, the libraries' first label as the low bit ((), (o1,), (o2,), (o1, o2) for two others).
    Bounded passes as cross_pairs; exact counts only (`scale` 1, else ValueError)."""
    _exact(scale)
    labels = list(libs)
    f = labels.index(focal)
    others = [x for x in labels if x != focal]
    counts = _with_threads(threads, _tallies, [libs[x] for x in labels], f, k, mem_gb)
    return {tuple(x for j, x in enumerate(others) if pat >> j & 1): counts[:, pat] for pat in range(counts.shape[1])}


# ---- one scan over every library: each pair's candidates and each trio's single-copy tallies ----


class SharedKmers(NamedTuple):
    """Every library pair's candidate record pairs (keyed (x, y), x before y in the libraries' order; as cross_pairs
    gives them), every trio's single-copy tallies of each member (keyed trio labels -> focal -> as
    single_copy_tallies on the trio's three libraries gives them), and the bounded passes the scan made."""

    pairs: dict[tuple[str, str], PairWeights]
    tallies: dict[tuple[str, str, str], dict[str, dict[tuple[str, ...], np.ndarray]]]
    passes: int


@njit(parallel=True, cache=True)
def _valid_kmers(codes, offs, k):
    """Each record's count of k-mers made of A, C, G and T only (no hashing)."""
    n = offs.size - 1
    out = np.zeros(n, np.int64)
    for i in prange(n):
        run = 0
        c = 0
        for p in range(offs[i], offs[i + 1]):
            if codes[p] > 3:
                run = 0
            else:
                run += 1
                if run >= k:
                    c += 1
        out[i] = c
    return out


@njit(parallel=True, cache=True)
def _pass_counts(codes, offs, k, pass_of, out):
    """out[i, p] += record i's canonical k-mers whose hash's top 16 bits fall in pass p (pass_of: top 16 bits ->
    pass), all passes in one scan."""
    n = offs.size - 1
    mask = (np.uint64(1) << np.uint64(2 * k)) - np.uint64(1)
    sh = np.uint64(2 * (k - 1))
    for i in prange(n):
        f = np.uint64(0)
        r = np.uint64(0)
        run = 0
        for p in range(offs[i], offs[i + 1]):
            c = codes[p]
            if c > 3:
                run = 0
                f = np.uint64(0)
                r = np.uint64(0)
                continue
            cu = np.uint64(c)
            f = ((f << np.uint64(2)) | cu) & mask
            r = (r >> np.uint64(2)) | ((np.uint64(3) - cu) << sh)
            run += 1
            if run >= k:
                h = mix64(min(f, r))
                out[i, pass_of[np.int64(h >> np.uint64(48))]] += 1


@njit(parallel=True, cache=True)
def _bucket_sort_par(H, E, blo, nb, nchunk):
    """_bucket_sort in parallel: each chunk of the input counts its buckets, then writes its records at its own
    offsets in each bucket (records keep their input order within a bucket). Returns the moved records and the bucket
    bounds."""
    n = H.size
    hist = np.zeros((nchunk, nb), np.int64)
    for ch in prange(nchunk):
        for j in range(n * ch // nchunk, n * (ch + 1) // nchunk):
            hist[ch, np.int64(H[j] >> np.uint64(48)) - blo] += 1
    bounds = np.zeros(nb + 1, np.int64)
    pos = np.empty((nchunk, nb), np.int64)
    for b in range(nb):
        q = bounds[b]
        for ch in range(nchunk):
            pos[ch, b] = q
            q += hist[ch, b]
        bounds[b + 1] = q
    Hs = np.empty_like(H)
    Es = np.empty_like(E)
    for ch in prange(nchunk):
        for j in range(n * ch // nchunk, n * (ch + 1) // nchunk):
            b = np.int64(H[j] >> np.uint64(48)) - blo
            q = pos[ch, b]
            Hs[q] = H[j]
            Es[q] = E[j]
            pos[ch, b] = q + 1
    return Hs, Es, bounds


@njit(parallel=True, cache=True)
def _joint_pairs(H, E, bounds, lo, cmax, write, cnt, B, W, nchunk):
    """Over buckets sorted by hash, runs of one hash (one k-mer) with c[l] occurrences in library l (copies lo[l] ..
    lo[l + 1] - 1, n in all): for every two libraries l < m whose counts are both 1..cmax, each of their c[l] x c[m]
    occurrence pairs gets the weight 1/(c[l] c[m]) rounded to float32. Slot m * n + x holds the pairs of l's copy x
    (a global copy) in library m. write=False: count per (chunk, slot) into cnt; write=True: write m's own copy index
    and the weight at cnt[chunk, slot], advancing it."""
    nb = bounds.size - 1
    nlib = lo.size - 1
    n = lo[nlib]
    for ch in prange(nchunk):
        c = np.zeros(nlib, np.int64)
        occ = np.empty((nlib, cmax), np.int64)
        for b in range(nb * ch // nchunk, nb * (ch + 1) // nchunk):
            s, e = bounds[b], bounds[b + 1]
            j = s
            while j < e:
                t = j + 1
                while t < e and H[t] == H[j]:
                    t += 1
                if t - j >= 2:
                    for q in range(nlib):
                        c[q] = 0
                    for u in range(j, t):
                        y = np.int64(E[u])
                        q = 0
                        while y >= lo[q + 1]:
                            q += 1
                        if c[q] < cmax:
                            occ[q, c[q]] = y
                        c[q] += 1
                    for l1 in range(nlib):
                        if c[l1] < 1 or c[l1] > cmax:
                            continue
                        for l2 in range(l1 + 1, nlib):
                            if c[l2] < 1 or c[l2] > cmax:
                                continue
                            w = np.float32(1.0 / (c[l1] * c[l2]))
                            for u in range(c[l1]):
                                slot = l2 * n + occ[l1, u]
                                if write:
                                    for v in range(c[l2]):
                                        o = cnt[ch, slot]
                                        B[o] = occ[l2, v] - lo[l2]
                                        W[o] = w
                                        cnt[ch, slot] = o + 1
                                else:
                                    cnt[ch, slot] += c[l2]
                j = t


@njit(parallel=True, cache=True)
def _joint_tallies(H, E, bounds, lo, trios, base, gstart, out, nchunk):
    """Over buckets sorted by hash, runs of one hash: for each trio of libraries (a row of trios) and each member f, a
    k-mer occurring once in f and at most once in each of the other two adds 1 to out[chunk, base[trio, member] +
    f's own copy index, pattern], pattern bit 0 for the trio's first other library holding it, bit 1 for the
    second. Work items: the trios in groups (group g: rows gstart[g] .. gstart[g + 1] - 1 of trios, whose rows of out
    are their own) by the buckets in nchunk chunks; item (g, chunk) counts group g's trios over the chunk's buckets,
    so no two items write one row of one chunk's counters."""
    nb = bounds.size - 1
    nlib = lo.size - 1
    ngroup = gstart.size - 1
    for item in prange(ngroup * nchunk):
        g = item // nchunk
        ch = item - g * nchunk
        t0, t1 = gstart[g], gstart[g + 1]
        c = np.zeros(nlib, np.int64)
        first = np.zeros(nlib, np.int64)
        for b in range(nb * ch // nchunk, nb * (ch + 1) // nchunk):
            s, e = bounds[b], bounds[b + 1]
            j = s
            while j < e:
                t = j + 1
                while t < e and H[t] == H[j]:
                    t += 1
                for q in range(nlib):
                    c[q] = 0
                for u in range(j, t):
                    y = np.int64(E[u])
                    q = 0
                    while y >= lo[q + 1]:
                        q += 1
                    if c[q] == 0:
                        first[q] = y
                    c[q] += 1
                for tt in range(t0, t1):
                    for r in range(3):
                        f = trios[tt, r]
                        if c[f] != 1:
                            continue
                        o1 = trios[tt, 1] if r == 0 else trios[tt, 0]
                        o2 = trios[tt, 1] if r == 2 else trios[tt, 2]
                        if c[o1] > 1 or c[o2] > 1:
                            continue
                        out[ch, base[tt, r] + first[f] - lo[f], c[o1] + 2 * c[o2]] += 1
                j = t


def _tally_split(n_trios, width, threads):
    """(groups, chunks) of the tallies' work items: chunks of the buckets, each with an int32 copy of every trio's
    rows (16 bytes a row) while those copies fit TALLY_CHUNK_BYTES, by groups of trios, so that groups x chunks comes
    as close to `threads` as it can (without passing it), with the fewest groups (each group rescans its chunk's
    k-mers). One trio: one group, as many chunks as fit."""
    fit = TALLY_CHUNK_BYTES // (16 * width)
    best = (1, 1)
    for chunks in range(1, max(1, min(threads, fit)) + 1):
        groups = max(1, min(n_trios, threads // chunks))
        if groups * chunks >= best[0] * best[1]:
            best = (groups, chunks)
    return best


def _shared(libs, labels, k, cmax, mem_gb):
    m = len(libs)
    sizes = [len(lib.offs) - 1 for lib in libs]
    lo = np.concatenate(([0], np.cumsum(sizes))).astype(np.int64)
    n = int(lo[-1])
    if n >= 2 ** 31:
        raise ValueError("more than 2^31 - 1 records in the libraries")
    valid = [_valid_kmers(lib.codes, lib.offs, k) for lib in libs]
    total = int(sum(int(v.sum()) for v in valid))
    n_pass = min(65536, max(1, math.ceil(total / max(1, int(mem_gb * 1e9 / KMER_BYTES)))))
    edges = np.linspace(0, 65536, n_pass + 1).astype(np.int64)
    log.debug("shared k-mers of %d libraries: %d k-mers in %d pass(es)", m, total, n_pass)
    if n_pass == 1:
        counts = [v[:, None] for v in valid]
    else:
        pass_of = (np.searchsorted(edges, np.arange(65536), side="right") - 1).astype(np.int64)
        counts = []
        for lib, size in zip(libs, sizes):
            counts.append(np.zeros((size, n_pass), np.int64))
            _pass_counts(lib.codes, lib.offs, k, pass_of, counts[-1])
    trios = np.array(list(combinations(range(m), 3)), np.int64).reshape(-1, 3)
    widths = np.array([[sizes[f] for f in t] for t in trios.tolist()], np.int64).reshape(-1, 3)
    base = (np.cumsum(widths.ravel()) - widths.ravel()).reshape(-1, 3).astype(np.int64)
    width = int(widths.sum())
    tally = np.zeros((width, 4), np.int64)
    pairs = [(i, j) for i in range(m) for j in range(i + 1, m)]
    found = {p: (np.empty(0, np.uint64), np.empty(0, np.float64)) for p in pairs}
    threads, cut = numba.get_num_threads(), _cut(1)
    split = _tally_split(len(trios), width, threads) if width else (1, 1)
    gstart = (np.arange(split[0] + 1) * len(trios) // split[0]).astype(np.int64)
    if width:
        log.debug("single-copy tallies of %d trio(s): %d group(s) x %d chunk(s) of the buckets", len(trios), *split)
    for p in range(n_pass):
        blo, bhi = int(edges[p]), int(edges[p + 1])
        tot = sum(int(c[:, p].sum()) for c in counts)
        if bhi <= blo or tot == 0:
            continue
        H, E = np.empty(tot, np.uint64), np.empty(tot, np.int32)
        start = 0
        for lib, c, first in zip(libs, counts, lo[:-1].tolist()):
            off = np.zeros(c.shape[0], np.int64)
            off[1:] = np.cumsum(c[:, p])[:-1]
            off += start
            _kmer_pass(lib.codes, lib.offs, k, blo, bhi, cut, True, off, H, E)
            used = int(c[:, p].sum())
            if first:
                E[start:start + used] += np.int32(first)
            start += used
        H, E, bounds = _bucket_sort_par(H, E, blo, bhi - blo, threads)
        _sort_buckets(H, E, bounds, threads)
        if width:
            nct = split[1]
            if nct == 1:                                # each group adds into its own rows of the tallies
                _joint_tallies(H, E, bounds, lo, trios, base, gstart, tally.reshape(1, width, 4), 1)
            else:
                out = np.zeros((nct, width, 4), np.int32)
                _joint_tallies(H, E, bounds, lo, trios, base, gstart, out, nct)
                tally += out.sum(0, dtype=np.int64)
                del out
        if not pairs:
            continue
        ncp = min(threads, PAIR_CHUNKS)
        cnt = np.zeros((ncp, m * n), np.int64)
        e32, ef = np.zeros(1, np.int32), np.zeros(1, np.float32)
        _joint_pairs(H, E, bounds, lo, cmax, False, cnt, e32, ef, ncp)
        slots = np.concatenate(([0], np.cumsum(cnt.sum(0)))).astype(np.int64)
        cnt = slots[:-1] + np.cumsum(cnt, 0) - cnt
        B, W = np.empty(int(slots[-1]), np.int32), np.empty(int(slots[-1]), np.float32)
        _joint_pairs(H, E, bounds, lo, cmax, True, cnt, B, W, ncp)
        del H, E, cnt
        nd = np.zeros(m * n, np.int64)
        _pair_sums(B, W, slots, False, nd, nd, B, np.zeros(1, np.float64))
        off = np.concatenate(([0], np.cumsum(nd))).astype(np.int64)
        UB, UW = np.empty(int(off[-1]), np.int32), np.empty(int(off[-1]), np.float64)
        _pair_sums(B, W, slots, True, nd, off, UB, UW)
        del B, W
        for i, j in pairs:
            s0, s1 = j * n + lo[i], j * n + lo[i + 1]
            a = np.repeat(np.arange(sizes[i], dtype=np.uint64), nd[s0:s1])
            pk = (a << np.uint64(32)) | UB[off[s0]:off[s1]].astype(np.uint64)
            found[i, j] = _merge_sum(*found[i, j], pk, UW[off[s0]:off[s1]])
    out_pairs = {(labels[i], labels[j]): PairWeights((u >> np.uint64(32)).astype(np.int64),
                                                     (u & np.uint64(0xFFFFFFFF)).astype(np.int64), w)
                 for (i, j), (u, w) in found.items()}
    out_tallies = {}
    for t, trio in enumerate(trios.tolist()):
        names = tuple(labels[f] for f in trio)
        out_tallies[names] = {}
        for r, f in enumerate(trio):
            block = tally[base[t, r]:base[t, r] + sizes[f]]
            others = [names[q] for q in range(3) if q != r]
            out_tallies[names][labels[f]] = {tuple(x for b, x in enumerate(others) if pat >> b & 1): block[:, pat].copy()
                                             for pat in range(4)}
    return SharedKmers(out_pairs, out_tallies, n_pass)


def shared_kmers(libs, *, k=K, cmax=CMAX_PAIR, scale=1, mem_gb=MEM_GB, threads=1) -> SharedKmers:
    """One bounded scan over every library of `libs` (label -> library of codes and offs, in order) that gives every
    pair's candidate record pairs (equal to cross_pairs on that pair) and, for every trio, each member's single-copy
    tallies (equal to single_copy_tallies on the trio's three libraries alone: "at most once in each other library"
    looks only at the trio's). Passes of at most mem_gb of (hash, copy) records, planned from a count of the k-mers that
    hashes nothing and sized by one counting scan for all passes, so each pass hashes the libraries once; numba loops on
    `threads` threads; exact counts only (`scale` 1, else ValueError)."""
    _exact(scale)
    return _with_threads(threads, _shared, list(libs.values()), list(libs), k, cmax, mem_gb)


def _with_threads(threads, fn, *args, **kw):
    """fn(*args, **kw) with numba's parallel loops on `threads` threads (capped at numba's pool), restored after."""
    prev = numba.get_num_threads()
    set_numba_threads(threads)
    try:
        return fn(*args, **kw)
    finally:
        numba.set_num_threads(prev)


def _digest(lib, seqof):
    """sha256 of what the clades are built from: record boundaries, each copy's sequence index and the codes."""
    h = hashlib.sha256()
    for a in (np.asarray(lib.offs, np.int64), seqof, np.asarray(lib.codes, np.uint8)):
        a = np.ascontiguousarray(a)
        h.update(a.size.to_bytes(8, "little"))
        h.update(a)
    return h.hexdigest()


def _load(path):
    """The clades saved at path (np.savez of their arrays), or None when the file cannot be read."""
    try:
        with np.load(path) as z:
            return Clades(z["members"], z["size"], z["mult"], z["nseq"], int(z["n_occ"]), int(z["n_multi"]))
    except (OSError, ValueError, KeyError, EOFError, zipfile.BadZipFile) as e:
        log.warning("unreadable clade checkpoint %s (%s); rebuilding it", path, type(e).__name__)
        return None


def clades(lib, seq_of, *, k=K, cap=CAP, min_seq=2, scale="auto", mem_gb=MEM_GB, threads, cache=None) -> Clades:
    """The clades of a library (codes and offs) whose copy i lies on sequence seq_of[i] (module docstring): k-mer
    length k, carrier sets of 2..cap copies on >= min_seq sequences, FracMinHash scale ("auto": auto_scale of the
    library's bases; 1 = exact), passes of at most mem_gb of (hash, copy) records, numba loops on `threads` threads.
    With `cache` (a folder), a checkpoint stage there; the result does not depend on mem_gb or threads."""
    seqof = np.asarray(seq_of, np.int64)
    sc = auto_scale(int(lib.codes.size)) if scale == "auto" else int(scale)
    params = {"k": int(k), "cap": int(cap), "min_seq": int(min_seq), "scale": sc}
    if cache is None:
        cl = _with_threads(threads, _build, lib, seqof, mem_gb=mem_gb, **params)
    else:
        out = Path(cache) / "clades.npz"
        done = marker(out)
        key = stage_key("clades", {"library": _digest(lib, seqof)}, params, code_fingerprint("seq.kmers"))
        cl = _load(out) if is_done(done, key) else None
        if cl is not None:
            log.debug("clades: checkpoint %s holds; loaded", out)
        else:
            clear_done(done)
            cl = _with_threads(threads, _build, lib, seqof, mem_gb=mem_gb, **params)
            out.parent.mkdir(parents=True, exist_ok=True)
            with atomic_path(out) as tmp, open(tmp, "wb") as fh:
                np.savez(fh, members=cl.members, size=cl.size, mult=cl.mult, nseq=cl.nseq, n_occ=cl.n_occ,
                         n_multi=cl.n_multi)
            mark_done(done, key, [out])
    cl.scale = sc
    return cl


@njit(parallel=True, cache=True)
def _sketches(codes, offs, k, s):
    """Bottom-s MinHash of each record's canonical k-mers (rows padded with the maximum uint64)."""
    n = offs.size - 1
    pad = np.uint64(0xFFFFFFFFFFFFFFFF)
    out = np.full((n, s), pad, np.uint64)
    mask = np.uint64((1 << (2 * k)) - 1)
    for i in prange(n):
        a0 = offs[i]
        a1 = offs[i + 1]
        m = a1 - a0 - k + 1
        if m <= 0:
            continue
        h = np.empty(m, np.uint64)
        nh = 0
        fw = np.uint64(0)
        rv = np.uint64(0)
        good = 0
        for p in range(a0, a1):
            c = codes[p]
            if c > 3:
                good = 0
                fw = np.uint64(0)
                rv = np.uint64(0)
                continue
            fw = ((fw << np.uint64(2)) | np.uint64(c)) & mask
            rv = (rv >> np.uint64(2)) | (np.uint64(3 - c) << np.uint64(2 * (k - 1)))
            good += 1
            if good >= k:
                h[nh] = mix64(min(fw, rv))
                nh += 1
        if nh == 0:
            continue
        h = np.sort(h[:nh])
        j = 0
        last = pad
        for q in range(nh):
            if h[q] != last:
                out[i, j] = h[q]
                last = h[q]
                j += 1
                if j == s:
                    break
    return out


def sketch(lib, *, k=SKETCH_K, s=SKETCH_S):
    """Bottom-s MinHash sketch of each record of a library (codes and offs): the s smallest distinct mix64 values of
    its canonical k-mers, ascending, rows padded with the largest uint64."""
    return _sketches(lib.codes, lib.offs, k, s)
