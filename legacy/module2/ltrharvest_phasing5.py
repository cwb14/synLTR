#!/usr/bin/env python3
"""
v5.py - subgenome phasing and ancestry painting of an allopolyploid from its LTR-RT library.

Idea
----
While the progenitors were apart, each amplified its own LTR-RT lineages. Copies of such a lineage share
rare sequence variants (k-mers carried by few copies = clades) and are each other's closest relatives,
and they sit on chromosomes of one subgenome only. v5 reads that signal with no homoeolog configuration,
no number of subgenomes and no progenitors.

1. Structure. Chromosomes are linked by (a) every clade of <= w copies spanning them (1/(m-1) per pair,
   each clade counted once; w chosen per genome by split-half reproducibility) and (b) nearest-relative
   votes (each copy votes for the 3 other chromosomes holding its closest relatives), added with equal
   weight; copies on the same chromosome never link. Links are compared with a degree-corrected null and
   chromosomes are split by maximum within-group excess (exact for N=2). The number of subgenomes is
   found by divisive splitting: a split is kept when independent halves of the library reproduce it
   (adjusted Rand >= 0.8) and, below the top, when it is >= 0.5x as strong as its parent. Every tested
   split is reported (structure_tree), so nested or partly autopolyploid structure is visible rather than
   forced. Support = fraction of 100 half-libraries reproducing each chromosome's assignment. Each
   chromosome's mutually closest partner within its subgenome is listed with the fraction of
   half-libraries that reproduce the pair (partner_support; pairs >= 0.9 are reproducible and count as one
   chromosome in all ancestry scores); a large excess marks near-identical chromosomes (an autopolyploid
   component, or two progenitors too close for LTR-RTs to separate).
2. Copy ancestry. Each LTR-RT's k-mers are scored two ways against the copies on the OTHER chromosomes:
   abundance (how many copies of each subgenome carry them) and breadth (how many chromosomes of each
   subgenome carry them, ignoring k-mers a homoeolog could explain). Each score is calibrated against the
   chromosomes' own assignments (temperature, an evidence exponent for correlated k-mers, intercepts); the
   two are pooled log-linearly and recalibrated. Calls at posterior >= 0.95: own / foreign / unresolved;
   no_relatives (no informative k-mer); assigned_unphased (copies on unphased sequences).
3. Ancestry along chromosomes. A hidden Markov model over the copies gives LTR-only ancestry segments;
   a candidate exchange needs a decisive likelihood ratio (>= 100) counting each independent lineage
   once, from >= 2 lineages. With --genome, every 50-kb window is painted from genome-wide k-mers confined
   to one subgenome among the other chromosomes (all repeat classes; homoeologous single-copy sequence is
   neutral by construction); emissions are the empirical frequencies of window compositions on each
   subgenome's chromosomes (cross-fitted, so a chromosome's own exchanges cannot look typical); a candidate
   exchange needs a likelihood ratio >= 100 and mean posterior >= 0.99 (segments with the ratio alone are
   listed as 'weak'). Unplaced sequences >= 50 kb are painted too; from their LTR-RTs alone, unphased
   sequences with informative copies get an ancestry call (NA when the evidence is weak). Exchanges and
   assembly switch errors cannot be told apart from sequence.
4. Ages (with --k2p). The age window of own-ancestry copies per subgenome (the SubPhaser divergence-
   hybridization window) and the cross-subgenome transposition clock (step in the probability that a
   copy carries another subgenome's ancestry, fitted at the resolution of each copy's LTR length): a
   minimum merger age (reported only when its interval excludes zero), and, when supported, a tentative
   minimum age of progenitor divergence (where lineages become shared again).
5. TE biology from the same evidence: cross-subgenome transposition by direction (copies called foreign
   outside candidate exchange regions, with those younger than a dated merger counted separately; a lower
   bound, since calls need posterior >= 0.95), each subgenome's insertion history (observed LTR divergence of
   the copies on its chromosomes; young ages blurred by LTR length)
   and a per-family table (copies, subgenome enrichment, calls, ages). Not estimated: insertion vs removal
   rates or half-lives - intact copies alone cannot separate insertion history from survival.

Chromosome status: assigned (support >= 0.9), ambiguous, mosaic (candidate exchanges), conflict
(painted majority differs from the assignment). Nothing assumes equal subgenome sizes, homoeolog pairs
or one ancestry per chromosome.

Input FASTA headers: chrom:start-end[#Class/Superfamily/Family] (LTRquest, LTR_retriever, EDTA style).

Example
-------
  python v5.py --ltr_fasta ltr.fa --outdir out -t 8                           # structure + ancestry
  python v5.py --ltr_fasta ltr.fa --outdir out --k2p ltr.tsv --genome asm.fa  # + ages, genome painting
  python v5.py --ltr_fasta subA_ltr.fa subB_ltr.fa --outdir out -n 2          # several files, N given

Outputs: chromosomes.tsv, elements.tsv, segments.tsv, links_oe.tsv, summary.json, families.tsv,
cross_subgenome.tsv, insertion_history.tsv (with --k2p), unphased_sequences.tsv (when some sequences are
too small to phase), fig_structure, fig_painting, fig_ages (.pdf/.png) with
<figure>_legend.md; with --genome also segments_genome.tsv, genome_sequences.tsv, genome_windows.tsv.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
import re
import sys
import time

import numpy as np
from scipy import sparse
from scipy.optimize import linear_sum_assignment

__version__ = "5.0"

# ---------------------------------------------------------------- constants (one place, one reason each)
K_DEFAULT = 21          # k-mer length: unique in plant genomes yet short enough to span one SNP per clade variant
WINDOWS = (10, 25, 50, 100, 250)   # candidate carrier-count caps for phasing links (chosen per genome by split-half r)
MIN_ELEMENTS = 10       # sequences with fewer LTR-RTs are not phased (cannot be assessed)
MIN_REP = 0.8           # auto-N: split-half adjusted Rand index needed to accept a split (heuristic, from v4)
MIN_REL = 0.5           # auto-N: a sub-split must be >= this fraction as strong as its parent (heuristic, from v4)
BF_CALL = 3.0           # Bayes factor for a per-copy call ('positive evidence', Kass & Raftery 1995)
BF_DECISIVE = 100.0     # likelihood ratio for a segment call ('decisive', Jeffreys 1961 App. B)
MU_DEFAULT = 1.3e-8     # LTR substitutions/site/year if --mu is not given (rice LTRs, Ma & Bennetzen 2004 PNAS)
CLOCK_BOOT = 200        # bootstrap replicates of the binned clock (split-step interval and support; as v4)
MERGER_BOOT = 50        # bootstrap replicates of the pooled merger fit (interval stable to +-1 grid step; each re-fits the age NPMLE)

# ---------------------------------------------------------------- logging
_T0 = time.time()
VERBOSE = False


def log(msg):
    print(f"[v5 {time.time() - _T0:7.1f}s] {msg}", file=sys.stderr, flush=True)


def vlog(msg):
    if VERBOSE:
        log(msg)


def die(msg):
    print(f"ERROR: {msg}", file=sys.stderr, flush=True)
    sys.exit(1)


def natural_key(s):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", s)]


# ---------------------------------------------------------------- input
def read_fasta(path):
    """Yield (id, sequence bytes); id is the first header token. Plain or gzip/bgzip."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rb") as f:
        h, buf = None, []
        for line in f:
            if line[:1] == b">":
                if h is not None:
                    yield h, b"".join(buf)
                tok = line[1:].split()
                h, buf = (tok[0].decode() if tok else ""), []
            else:
                buf.append(line.strip())
        if h is not None:
            yield h, b"".join(buf)


def parse_id(h):
    """chrom:start-end[#cls] or chrom:start..end[#cls] -> (chrom, start, end, cls) or None."""
    core, _, cls = h.partition("#")
    chrom, sep, coords = core.rpartition(":")
    s, dash, e = coords.partition("..") if ".." in coords else coords.partition("-")
    if not sep or not dash or not chrom or not s.isdigit() or not e.isdigit():
        return None
    return chrom, int(s), int(e), cls or "NA"


def read_config(path):
    sets = []
    with open(path) as f:
        for raw in f:
            toks = raw.split("#", 1)[0].split()
            if toks:
                sets.append(toks)
    if not sets:
        die(f"no chromosome sets in config {path}")
    return sets


def _isnum(x):
    try:
        float(x)
        return True
    except ValueError:
        return False


def read_k2p(paths, k2p_col=None):
    """LTR id -> (K2P, s.e., aligned LTR sites, substitutions). Tables with a header naming seq_id/k2p
    (LTRquest, Kmer2LTR; k2p_se, n_sites|ltr5_len, n_ts+n_tv used when present); repeated headers
    (concatenated tables) are fine. Header-less tables: id in column 1, K2P in the last column, or in
    k2p_col (1-based; overrides any header). The K2P column is checked: values must lie in [0, 1]."""
    out, bad = {}, 0
    nan = float("nan")
    for path in paths:
        opener = gzip.open if path.endswith(".gz") else open
        ix, ki, ci = None, None, 0
        first = True
        with opener(path, "rt") as f:
            for line in f:
                if not line.strip():
                    continue
                is_first, first = first, False
                p = line.rstrip("\n").split("\t")
                head = [x.lstrip("#") for x in p]
                if k2p_col is None and "k2p" in [h.lower() for h in head]:   # (re)defines columns
                    head = [h.lower() for h in head]
                    ix = {c: head.index(c) for c in ("k2p_se", "n_sites", "ltr5_len", "n_ts", "n_tv") if c in head}
                    ki = head.index("k2p")
                    ci = head.index("seq_id") if "seq_id" in head else 0
                    continue
                if line.startswith("#"):
                    continue
                if k2p_col is not None:                     # user-given column (ids from column 1)
                    p = line.split()
                    kk, cc, jx = k2p_col - 1, 0, {}
                    if is_first and len(p) > kk and not _isnum(p[kk]) and p[kk].upper() not in ("NA", "NAN", ""):
                        continue                            # a header line (first line only)
                elif ki is None:                            # header-less: id + K2P last column
                    p = line.split()
                    kk, cc, jx = len(p) - 1, 0, {}
                else:
                    kk, cc, jx = ki, ci, ix
                if len(p) <= max(kk, cc):
                    bad += 1
                    continue
                if p[kk].upper() in ("NA", "NAN", "."):     # missing value: not an error
                    continue
                try:
                    k2p = float(p[kk])
                except ValueError:
                    bad += 1
                    continue

                def num(c):
                    try:
                        return float(p[jx[c]]) if c in jx else nan
                    except (ValueError, IndexError):
                        return nan
                if not 0.0 <= k2p <= 1.0:                   # not a divergence (wrong column?)
                    bad += 1
                    continue
                L = num("n_sites")
                if not L == L:
                    L = num("ltr5_len")
                out[p[cc].strip().lstrip(">")] = (k2p, num("k2p_se"), L, num("n_ts") + num("n_tv"))
    if not out:
        die("no K2P values read from " + ", ".join(paths)
            + (" (check --k2p_col)" if k2p_col is not None else " (no 'k2p' header and no numeric last column?)"))
    if bad > 0.05 * max(len(out) + bad, 1):
        die(f"{bad} of {len(out) + bad} K2P rows are unreadable or outside [0, 1]: the K2P column is probably wrong "
            + ("(check --k2p_col)" if k2p_col is not None else
               "(header-less tables are read from the last column; give --k2p_col)"))
    return out, bad


# ---------------------------------------------------------------- k-mer sketch
LUT = np.full(256, 4, np.uint8)
for _i, _b in enumerate(b"ACGT"):
    LUT[_b] = _i
    LUT[_b + 32] = _i
_M1 = np.uint64(0xFF51AFD7ED558CCD)
_M2 = np.uint64(0xC4CEB9FE1A85EC53)
_S33 = np.uint64(33)
_SH = np.uint64(24)            # low 24 bits of a key hold the element id
_EMASK = np.uint64((1 << 24) - 1)


def mix64(x):
    """murmur3 finalizer: a bijective 64-bit mixer."""
    x = x ^ (x >> _S33)
    x = x * _M1
    x = x ^ (x >> _S33)
    x = x * _M2
    return x ^ (x >> _S33)


def kmer_codes(c, k):
    """Forward and reverse-complement 2-bit codes of every k-window of c (uint64 0..3), built by
    doubling (log2 k vector passes instead of k)."""
    n = c.size - k + 1
    f, r = {1: c}, {1: np.uint64(3) - c}
    p = 1
    while p * 2 <= k:
        a, b = f[p], r[p]
        f[2 * p] = (a[:a.size - p] << np.uint64(2 * p)) | a[p:]
        r[2 * p] = b[:b.size - p] | (b[p:] << np.uint64(2 * p))
        p *= 2
    F = R = None
    off = 0
    for q in sorted((x for x in f if k & x), reverse=True):
        fq, rq = f[q][off:off + n], r[q][off:off + n]
        if F is None:
            F, R = fq.copy(), rq.copy()
        else:
            F = (F << np.uint64(2 * q)) | fq
            R = R | (rq << np.uint64(2 * off))
        off += q
    return F, R


def sampled_hashes(codes, k, scale):
    """Positions and hashes of FracMinHash-sampled canonical k-mers (hash < 2^64/scale) of a code
    array (4 = non-ACGT; k-mers spanning it are skipped)."""
    n = codes.size - k + 1
    if n <= 0:
        return np.empty(0, np.int64), np.empty(0, np.uint64)
    bad = codes == 4
    cs = np.zeros(codes.size + 1, np.int32)
    np.cumsum(bad, out=cs[1:])
    ok = (cs[k:k + n] - cs[:n]) == 0
    c = codes.astype(np.uint64)
    c[bad] = 0
    F, R = kmer_codes(c, k)
    h = mix64(np.minimum(F, R))
    keep = ok if scale <= 1 else ok & (h < np.uint64((1 << 64) // scale - 1))
    pos = np.flatnonzero(keep)
    return pos, h[pos]


def _sketch_chunk(args):
    codes, eids, k, scale = args
    pos, h = sampled_hashes(codes, k, scale)
    b = np.uint64(max(int(math.floor(math.log2(scale))), 0)) if scale > 1 else np.uint64(0)
    h = h << b                                   # sampled hashes are < 2^64/scale: keep their informative bits
    return np.unique(((h >> _SH) << _SH) | eids[pos].astype(np.uint64))


def pool_map(fn, jobs, threads):
    """Ordered map over jobs with a fork process pool (serial if threads <= 1 or fork unavailable)."""
    if threads <= 1:
        for j in jobs:
            yield fn(j)
        return
    import multiprocessing as mp
    try:
        ctx = mp.get_context("fork")
    except ValueError:
        for j in jobs:
            yield fn(j)
        return
    with ctx.Pool(threads) as pool:
        for r in pool.imap(fn, jobs, chunksize=1):
            yield r


class Library:
    """LTR-RT library: element table + sorted (k-mer hash, element) keys."""

    def __init__(self, paths, k, scale, threads=1, chunk=1 << 24):
        self.ids, self.chrom, self.start, self.end, self.cls, self.length = [], [], [], [], [], []
        bad_hdr = short = 0

        def batches():
            nonlocal bad_hdr, short
            buf, ebuf, nb = [], [], 0
            for path in paths:
                for h, s in read_fasta(path):
                    p = parse_id(h)
                    if p is None:
                        bad_hdr += 1
                        continue
                    if len(s) < k:
                        short += 1
                        continue
                    i = len(self.ids)
                    if i >= (1 << 24):
                        die("more than 16.7 million elements are not supported")
                    self.ids.append(h)
                    self.chrom.append(p[0])
                    self.start.append(p[1])
                    self.end.append(p[2])
                    self.cls.append(p[3])
                    self.length.append(len(s))
                    codes = LUT[np.frombuffer(s, np.uint8)]
                    buf += [codes, np.full(1, 4, np.uint8)]
                    ebuf.append(np.full(codes.size + 1, i, np.uint32))
                    nb += codes.size + 1
                    if nb >= chunk:
                        yield np.concatenate(buf), np.concatenate(ebuf), k, scale
                        buf, ebuf, nb = [], [], 0
            if buf:
                yield np.concatenate(buf), np.concatenate(ebuf), k, scale

        keys = list(pool_map(_sketch_chunk, batches(), threads))
        if bad_hdr:
            log(f"WARNING: skipped {bad_hdr} records whose header is not chrom:start-end[#class]")
        if short:
            vlog(f"skipped {short} records shorter than k={k}")
        if not self.ids:
            die("no usable LTR-RT records in " + ", ".join(paths))
        keys = np.concatenate(keys) if keys else np.empty(0, np.uint64)
        keys.sort()
        self.kmer = keys >> _SH
        self.eid = (keys & _EMASK).astype(np.int32)
        self.n = len(self.ids)
        self.start = np.array(self.start, np.int64)
        self.end = np.array(self.end, np.int64)
        self.length = np.array(self.length, np.int64)
        self.chroms = sorted(set(self.chrom), key=natural_key)


# ---------------------------------------------------------------- clades ("splits")
class Splits:
    """Distinct carrier sets of shared k-mers that span >= 2 chromosomes (clades of the copy tree).
    CSR: members[indptr[s]:indptr[s+1]] = elements carrying split s; mult[s] = sampled k-mers with
    exactly that carrier set; nch[s] = chromosomes spanned."""

    def __init__(self, kmer, eid, echrom, n_chrom, nmax, emask=None):
        if emask is not None:
            keep = emask[eid]
            kmer, eid = kmer[keep], eid[keep]
        self.n_shared_kmers = 0
        if kmer.size == 0:
            self._empty()
            return
        brk = np.flatnonzero(kmer[1:] != kmer[:-1]) + 1
        starts = np.concatenate(([0], brk))
        sizes = np.diff(np.concatenate((starts, [kmer.size])))
        E = int(eid.max()) + 1
        erand = mix64((np.arange(E, dtype=np.uint64) + np.uint64(1)) * np.uint64(0x9E3779B97F4A7C15))
        sh = np.bitwise_xor.reduceat(erand[eid], starts) ^ mix64(sizes.astype(np.uint64))
        gi = np.flatnonzero((sizes >= 2) & (sizes <= nmax))
        _, first, mult = np.unique(sh[gi], return_index=True, return_counts=True)
        o = np.argsort(first)
        sel, mult = gi[first[o]], mult[o]
        gmask = np.zeros(starts.size, bool)
        gmask[sel] = True
        members = eid[np.repeat(gmask, sizes)]
        ssz = sizes[sel]
        sid = np.repeat(np.arange(sel.size), ssz)
        ch = echrom[members]
        uk = np.unique(sid.astype(np.int64) * n_chrom + ch)
        nch = np.bincount(uk // n_chrom, minlength=sel.size)
        keep = nch >= 2
        kp = np.repeat(keep, ssz)
        self.members = members[kp].astype(np.int32)
        self.size = ssz[keep]
        self.indptr = np.concatenate(([0], np.cumsum(self.size)))
        self.mult = mult[keep]
        self.nch = nch[keep]
        self.n = int(keep.sum())
        self.n_shared_kmers = int((sizes >= 2).sum())
        self.sid = np.repeat(np.arange(self.n), self.size)

    def _empty(self):
        self.members = np.empty(0, np.int32)
        self.size = np.empty(0, np.int64)
        self.indptr = np.zeros(1, np.int64)
        self.mult = np.empty(0, np.int64)
        self.nch = np.empty(0, np.int64)
        self.n = 0
        self.sid = np.empty(0, np.int64)

    def view(self, mask):
        v = Splits.__new__(Splits)
        kp = np.repeat(mask, self.size)
        v.members, v.size, v.mult, v.nch = self.members[kp], self.size[mask], self.mult[mask], self.nch[mask]
        v.n = int(mask.sum())
        v.indptr = np.concatenate(([0], np.cumsum(v.size)))
        v.sid = np.repeat(np.arange(v.n), v.size)
        v.n_shared_kmers = self.n_shared_kmers
        return v


def affinity(sp, echrom, C, eweight=None):
    """C x C link matrix: each clade adds 1/(m-1) to every pair of the m chromosomes it spans, so a
    clade counts once however many copies or k-mers it has; same-chromosome copies never link."""
    data = np.ones(sp.members.size) if eweight is None else eweight[sp.members].astype(float)
    X = sparse.csr_matrix((data, (sp.sid, echrom[sp.members])), shape=(sp.n, C))
    X.data = (X.data > 0).astype(float)
    X.eliminate_zeros()
    m = np.asarray(X.sum(1)).ravel()
    w = np.where(m >= 2, 1.0 / np.maximum(m - 1, 1), 0.0)
    W = (X.T @ sparse.diags(w) @ X).toarray()
    np.fill_diagonal(W, 0)
    return W


def expected(W, iters=200):
    """Degree-corrected null for a zero-diagonal matrix: E_ij = b_i b_j with sum_{j!=i} E_ij = row sum."""
    d = W.sum(1)
    b = d / math.sqrt(max(d.sum(), 1e-300))
    for _ in range(iters):
        nb = np.sqrt(b * d / np.maximum(b.sum() - b, 1e-300))
        if np.allclose(nb, b, rtol=1e-12, atol=0):
            b = nb
            break
        b = nb
    E = np.outer(b, b)
    np.fill_diagonal(E, 0)
    return E


def log_oe(W):
    return np.log((W + 1) / (expected(W) + 1))


# ---------------------------------------------------------------- partition
def make_sets(chroms, config_sets, N):
    """Cannot-link sets (lists of chromosome indices); every chromosome in exactly one set."""
    ix = {c: i for i, c in enumerate(chroms)}
    sets, seen = [], set()
    for s in config_sets or []:
        mem = [ix[c] for c in s if c in ix and ix[c] not in seen]
        if len(mem) > N:
            die(f"config set {' '.join(s)} has more chromosomes than subgenomes (N={N})")
        if mem:
            sets.append(mem)
            seen.update(mem)
    sets += [[i] for i in range(len(chroms)) if i not in seen]
    return sets


def q_score(B, z):
    return 0.5 * float(B[z[:, None] == z[None, :]].sum())


def _local_search(B, z, sets, N, rng, max_rounds=200):
    z = z.copy()
    C = B.shape[0]
    for _ in range(max_rounds):
        changed = False
        for si in rng.permutation(len(sets)):
            mem = sets[si]
            oh = np.zeros((C, N))
            oh[np.arange(C), z] = 1
            oh[mem] = 0
            A = B[mem] @ oh
            if len(mem) == 1:
                new = [int(np.argmax(A[0]))]
            else:
                r, c = linear_sum_assignment(-A)
                new = list(c[np.argsort(r)])
            for c_, g in zip(mem, new):
                if z[c_] != g:
                    z[c_] = g
                    changed = True
        if not changed:
            break
    return z


def _exhaustive2(B, sets):
    """Exact optimum for N=2 by enumerating set orientations (<= 2^22)."""
    V = len(sets)
    C = B.shape[0]
    best, bz = -np.inf, None
    total = 1 << (V - 1)
    base = np.zeros(C)
    for c in sets[0]:
        base[c] = 1.0 if c == sets[0][0] else -1.0
    for lo in range(0, total, 1 << 16):
        codes = np.arange(lo, min(total, lo + (1 << 16)), dtype=np.int64)
        S = np.tile(base, (codes.size, 1))
        for v in range(1, V):
            bit = ((codes >> (v - 1)) & 1) * 2 - 1
            mem = sets[v]
            S[:, mem[0]] = bit
            if len(mem) == 2:
                S[:, mem[1]] = -bit
        q = np.einsum("ij,ij->i", S @ B, S)
        q[np.abs(S.sum(1)) == C] = -np.inf
        j = int(np.argmax(q))
        if q[j] > best:
            best, bz = q[j], S[j].copy()
    return (bz < 0).astype(int)


def _spectral_init(B, N, rng):
    _, V = np.linalg.eigh(B)
    X = V[:, -max(N - 1, 1):]
    if N == 2:
        return (X[:, -1] < 0).astype(int)
    cen = X[rng.choice(len(X), N, replace=False)]
    for _ in range(50):
        z = np.argmin(((X[:, None, :] - cen[None]) ** 2).sum(-1), 1)
        cen = np.array([X[z == g].mean(0) if (z == g).any() else X[rng.integers(len(X))] for g in range(N)])
    return z


def best_partition(B, N, sets, rng, restarts=30, init=None, exact=True):
    """Maximum within-group excess (degree-corrected modularity) partition into N groups."""
    C = B.shape[0]
    if exact and N == 2 and all(len(s) <= 2 for s in sets) and len(sets) <= 23:
        return _exhaustive2(B, sets)
    inits = [_spectral_init(B, N, rng)]
    if init is not None:
        inits.insert(0, np.asarray(init))
    inits += [rng.integers(0, N, C) for _ in range(restarts)]
    best, bz = -np.inf, None
    for z0 in inits:
        z = _local_search(B, z0, sets, N, rng)
        if len(set(z.tolist())) < N:
            continue
        q = q_score(B, z)
        if q > best + 1e-12:
            best, bz = q, z
    return bz if bz is not None else _local_search(B, inits[0], sets, N, rng)


def canon_labels(z, mass):
    """Relabel groups so that subgenome 1 holds the most LTR-RTs (ties: lowest index)."""
    groups = sorted(set(z.tolist()), key=lambda g: (-mass[z == g].sum(), np.flatnonzero(z == g)[0]))
    m = {g: i for i, g in enumerate(groups)}
    return np.array([m[g] for g in z])


def align(ref, z, N):
    cont = np.zeros((N, N))
    np.add.at(cont, (z, ref), 1)
    r, c = linear_sum_assignment(-cont)
    m = dict(zip(r, c))
    return np.array([m[g] for g in z])


def ari(a, b):
    """Adjusted Rand index (1 = identical, ~0 = chance)."""
    from scipy.special import comb
    a, b = np.asarray(a), np.asarray(b)
    _, ai = np.unique(a, return_inverse=True)
    _, bi = np.unique(b, return_inverse=True)
    ct = np.zeros((ai.max() + 1, bi.max() + 1))
    np.add.at(ct, (ai, bi), 1)
    s = comb(ct, 2).sum()
    sa, sb = comb(ct.sum(1), 2).sum(), comb(ct.sum(0), 2).sum()
    exp = sa * sb / comb(a.size, 2)
    mx = (sa + sb) / 2
    return 1.0 if mx == exp else float((s - exp) / (mx - exp))


# ---------------------------------------------------------------- parallel helpers (fork shares these)
_G = {}


def partner_pairs(W, halves, z):
    """Closest partner of each chromosome within its own subgenome: the chromosome it shares the most
    lineages with (observed/expected), kept when the two are each other's closest partner in the full
    library. Support = fraction of half-libraries in which they still are. Chromosome pairs that share
    young lineages far beyond the rest of their subgenome are near-identical copies (homologues of an
    autopolyploid component, or of two progenitors too close for LTR-RTs to separate); excess = the pair's
    observed/expected over the chromosome's median within its subgenome. Returns partner (-1 = none),
    support, excess."""
    C = W.shape[0]
    same = z[:, None] == z[None, :]
    np.fill_diagonal(same, False)

    def best_of(M):
        E = expected(M)
        oe = np.where(same & (E > 0), M / np.where(E > 0, E, 1.0), np.nan)
        ok = np.isfinite(oe).any(1)
        b = np.where(ok, np.nanargmax(np.where(np.isfinite(oe), oe, -np.inf), 1), -1)
        return b, oe
    b, oe = best_of(W)
    partner = np.array([j if j >= 0 and b[j] == i else -1 for i, j in enumerate(b)])
    support, excess = np.full(C, np.nan), np.full(C, np.nan)
    hb = [best_of(h)[0] for h in halves]
    for i, j in enumerate(partner):
        if j < 0:
            continue
        support[i] = float(np.mean([bb[i] == j and bb[j] == i for bb in hb])) if hb else np.nan
        rest = oe[i][np.isfinite(oe[i])]
        excess[i] = float(oe[i, j] / np.median(rest)) if rest.size > 1 and np.median(rest) > 0 else np.nan
    return partner, support, excess


def chrom_tests(W, z, N):
    """Per chromosome: mean O/E to own subgenome and to the closest other one; one-sided Welch test
    with the other chromosomes as replicates."""
    from scipy.stats import ttest_ind
    C = len(z)
    E = expected(W)
    L = np.log((W + 1) / (E + 1))
    OE = W / np.where(E > 0, E, np.nan)
    own_oe, oth_oe, pval = (np.full(C, np.nan) for _ in range(3))
    for c in range(C):
        others = np.arange(C) != c
        x = L[c, others & (z == z[c])]
        best, bg = -np.inf, None
        for g in range(N):
            if g != z[c] and (z == g).any() and L[c, z == g].mean() > best:
                best, bg = L[c, z == g].mean(), g
        if bg is None:
            continue
        y = L[c, others & (z == bg)]
        own_oe[c] = np.nanmean(OE[c, others & (z == z[c])]) if x.size else np.nan
        oth_oe[c] = np.nanmean(OE[c, others & (z == bg)])
        if x.size >= 2 and y.size >= 2:
            t, p = ttest_ind(x, y, equal_var=False)
            pval[c] = p / 2 if t > 0 else 1 - p / 2
    return own_oe, oth_oe, pval


# ---------------------------------------------------------------- dating (from v4)
def _clock_bins(a, se, y, m, edges):
    """Per K2P bin: median K2P, weighted mean share y, total weight, median s.e."""
    bi = np.clip(np.searchsorted(edges, a, side="right") - 1, 0, edges.size - 2)
    nb = edges.size - 1
    W = np.bincount(bi, weights=m, minlength=nb)
    Y = np.bincount(bi, weights=m * y, minlength=nb) / np.maximum(W, 1e-300)
    A = np.array([np.median(a[bi == k]) if (bi == k).any() else np.nan for k in range(nb)])
    S = np.array([np.median(se[bi == k]) if (bi == k).any() else np.nan for k in range(nb)])
    ok = W > 0
    return A[ok], Y[ok], W[ok], S[ok]


def _one_step(A, Y, W, S, grid, young_high):
    """Weighted least-squares single step between two levels, blurred by K2P
    s.e.: young_high -> share = hi*Phi((t-a)/s) + lo*(1-Phi), else mirrored.
    Returns (t, level_young, level_old) with the high side >= low side."""
    from scipy.stats import norm
    s = np.maximum(S, 1e-5)
    w = norm.cdf((grid[:, None] - A[None, :]) / s[None, :])          # P(younger than t)
    X1, X2 = w, 1 - w
    a11 = (W * X1 * X1).sum(1)
    a12 = (W * X1 * X2).sum(1)
    a22 = (W * X2 * X2).sum(1)
    b1 = (W * X1 * Y).sum(1)
    b2 = (W * X2 * Y).sum(1)
    det = a11 * a22 - a12 ** 2
    det = np.where(np.abs(det) < 1e-12, np.nan, det)
    ly = np.clip((b1 * a22 - b2 * a12) / det, 0, 1)
    lo = np.clip((b2 * a11 - b1 * a12) / det, 0, 1)
    sse = (W[None, :] * (Y[None, :] - X1 * ly[:, None] - X2 * lo[:, None]) ** 2).sum(1)
    bad = ~np.isfinite(sse) | ((ly < lo) if young_high else (ly > lo))
    sse[bad] = np.inf
    if not np.isfinite(sse).any():
        return None
    j = int(np.argmin(sse))
    return float(grid[j]), float(ly[j]), float(lo[j])


def _clock_from_bins(A, Y, W, S, grid):
    """Locate the confined trough (minimum of the 3-bin smoothed share), then fit
    the merger step on younger bins and the divergence step on older bins."""
    k = np.convolve(W * Y, np.ones(3), "same") / np.maximum(np.convolve(W, np.ones(3), "same"), 1e-300)
    t = int(np.argmin(k))
    young, old = slice(0, t + 1), slice(t, A.size)
    gy = grid[grid <= (A[min(t + 1, A.size - 1)] if t + 1 < A.size else A[-1])]
    go = grid[grid >= A[t]]
    m = _one_step(A[young], Y[young], W[young], S[young], gy, True) if t >= 1 and gy.size else None
    d = _one_step(A[old], Y[old], W[old], S[old], go, False) if A.size - t >= 2 and go.size else None
    return m, d, float(A[t])


def fit_clock(k2p, se, share, m, rng, boot=200, n_bins=50):
    """Cross-subgenome transposition clock (population level).
    share = fraction of an LTR-RT's closest relatives on other-subgenome
    chromosomes; m = evidence weight. Copies inserted after the merger come from
    lineages free to land in any subgenome (high share); copies inserted between
    progenitor divergence and the merger belong to lineages confined to one
    progenitor (low share, the trough); copies older than the divergence belong
    to lineages present in both progenitors (share rises again, tentative: older
    copies also have noisier relatives). Each step is a weighted least-squares
    fit on K2P bins blurred by K2P s.e.; 95% CIs by bootstrap over LTR-RTs."""
    ok = np.isfinite(k2p) & np.isfinite(se) & (m > 0)
    if ok.sum() < 200:
        return None
    a, s, y, w = k2p[ok], np.maximum(se[ok], 1e-5), share[ok], m[ok]
    hi = float(np.quantile(a, 0.995))
    q = np.unique(np.quantile(a[a > 0], np.linspace(0, 1, n_bins))) if (a > 0).any() else np.array([])
    edges = np.unique(np.concatenate([[0.0, 1e-12], q[(q > 1e-12) & (q < hi)], [hi, np.inf]]))
    grid = np.unique(np.concatenate([[0.0], np.linspace(0, hi, 400), edges[1:-1]]))

    def one(idx):
        return _clock_from_bins(*_clock_bins(a[idx], s[idx], y[idx], w[idx], edges), grid)

    mg, dv, trough = one(np.arange(a.size))
    R = []
    for _ in range(boot):
        mm, dd, _t = one(rng.integers(0, a.size, a.size))
        R.append((mm[0] if mm else np.nan, (mm[1] - mm[2]) if mm else np.nan,
                  dd[0] if dd else np.nan, (dd[2] - dd[1]) if dd else np.nan))
    R = np.array(R, float).reshape(-1, 4)
    ci = lambda c: [float(np.nanquantile(R[:, c], 0.025)), float(np.nanquantile(R[:, c], 0.975))] \
        if np.isfinite(R[:, c]).any() else [None, None]
    A, Y, Wt, S = _clock_bins(a, s, y, w, edges)
    out = dict(n=int(ok.sum()), k2p_trough=trough, k2p_max=hi,
               bins=dict(k2p=A.tolist(), share=Y.tolist(), weight=Wt.tolist()))
    if mg:
        out.update(tau_merger=mg[0], tau_merger_ci=ci(0), share_post=mg[1], share_confined=mg[2],
                   step_merger_ci=ci(1), merger_supported=bool(ci(1)[0] is not None and ci(1)[0] > 0))
    if dv:
        out.update(tau_divergence=dv[0], tau_divergence_ci=ci(2), share_ancestral=dv[2],
                   step_divergence_ci=ci(3), divergence_supported=bool(ci(3)[0] is not None and ci(3)[0] > 0))
    return out


def _age_lik(k, L, grid):
    """Poisson likelihood of k substitutions over L LTR sites at each divergence on the grid (rows scaled)."""
    from scipy.special import gammaln
    lam = L[:, None] * grid[None, :]
    logP = k[:, None] * np.log(np.maximum(lam, 1e-300)) - lam - gammaln(k + 1)[:, None]
    if grid[0] == 0:
        logP[:, 0] = np.where(k == 0, 0.0, -np.inf)
    return np.exp(logP - logP.max(1, keepdims=True))


def _npmle_ip(P, w, mu_end=1e-12):
    """Grid NPMLE by a primal log-barrier interior-point method (the approach of Koenker & Mizera 2014, REBayes):
    min -w.log(Px) + 1'x - mu sum(log x) over x > 0 (its optimum has 1'x = 1 + m mu), Newton steps in scaled
    variables x = X u (well conditioned however peaked the likelihoods), fraction-to-boundary 0.99 and Armijo
    backtracking; mu falls 10-fold whenever the Newton decrement is small, down to m mu <= mu_end (the duality
    gap). ~100-150 Newton steps, each one m x m solve."""
    m = P.shape[1]
    x = np.full(m, 1.0 / m)
    mu = 0.1 / m
    sw = np.sqrt(w)

    def f(z, mu_):
        return -float(w @ np.log(np.maximum(P @ z, 1e-300))) + z.sum() - mu_ * float(np.log(z).sum())
    while True:
        for _ in range(100):
            d = np.maximum(P @ x, 1e-300)
            g = 1.0 - P.T @ (w / d) - mu / x
            A = P * (sw / d)[:, None] * x[None, :]
            M = A.T @ A
            M[np.diag_indices(m)] += mu
            try:
                du = -np.linalg.solve(M, x * g)
            except np.linalg.LinAlgError:
                du = -np.linalg.lstsq(M, x * g, rcond=None)[0]
            dx = x * du
            lam2 = -float(g @ dx)                  # squared Newton decrement
            if not lam2 > 1e-14:
                break
            neg = dx < 0
            t = min(1.0, 0.99 * float(np.min(-x[neg] / dx[neg]))) if neg.any() else 1.0
            f0 = f(x, mu)
            while f(x + t * dx, mu) > f0 - 0.25 * t * lam2 and t > 1e-12:
                t *= 0.5
            x = x + t * dx
            if lam2 < 1e-10 * mu:
                break
        if m * mu <= mu_end:
            break
        mu /= 10
    return x / x.sum()


def _age_em(P, cnt, pi=None, sqp=True):
    """Nonparametric age distribution (Kiefer-Wolfowitz NPMLE on the grid) from likelihood rows P of groups of
    identical copies (multiplicity cnt). Returns grid weights, per-group posteriors, the optimality gap
    max_j sum_i w_i P_ij / (P pi)_i - 1 (0 at the optimum; it bounds the per-copy distance to the maximum
    log-likelihood) and whether mix-SQP certified the solution.
    Solver: mix-SQP (Kim et al. 2020; fast) accepted only when it certifies the optimum (gap <= 1e-6); otherwise
    the interior-point method (_npmle_ip), which converges whatever the conditioning. On very peaked likelihoods
    (long LTRs, ~1e5 copies; Aegilops) mix-SQP stalls far from the optimum and EM/SQUAREM converges sublinearly
    (Ae. ventricosa: gap 4e-5 after 1,000 SQUAREM iterations, hours per fit), while the interior point reaches
    gap < 1e-12 in ~10 s. Grid points whose directional derivative is clearly negative (< -1e-6; KKT: no mass
    at the optimum) are set to zero. sqp=False skips mix-SQP (bootstrap re-fits of a matrix it failed on)."""
    from scipy.optimize import nnls
    m = P.shape[1]
    w = cnt / cnt.sum()

    def deriv(z):
        return P.T @ (w / np.maximum(P @ z, 1e-300)) - 1
    gp = np.inf
    if sqp:                                    # 1. mix-SQP
        x = np.full(m, 1.0 / m) if pi is None else pi.copy()
        eps = 1e-8
        phi = lambda z: -float(w @ np.log(P @ z + eps)) + z.sum()
        for _ in range(50):
            d = P @ x + eps
            g = 1.0 - P.T @ (w / d)
            if -g.min() < 1e-8:
                break
            A = P * (np.sqrt(w) / d)[:, None]
            H = A.T @ A
            c = g - H @ x
            try:
                R = np.linalg.cholesky(H + 1e-10 * np.trace(H) / m * np.eye(m)).T
                y = nnls(R, -np.linalg.solve(R.T, c), maxiter=50 * m)[0]
            except (np.linalg.LinAlgError, RuntimeError):
                break
            st, f0, a = y - x, phi(x), 1.0
            while phi(x + a * st) > f0 + 1e-2 * a * (g @ st) and a > 1e-10:
                a *= 0.5
            x = x + a * st
        x = np.maximum(x, 0)
        x = x / x.sum() if x.sum() > 0 else np.full(m, 1.0 / m)
        gp = float(deriv(x).max())
    ok = gp <= 1e-6
    if not ok:                                 # 2. interior point
        x = _npmle_ip(P, w)
        x = np.where(deriv(x) < -1e-6, 0.0, x)
        x = x / x.sum()
        gp = float(deriv(x).max())
    R_ = P * x
    return x, R_ / np.maximum(R_.sum(1, keepdims=True), 1e-300), gp, ok


def _merger_sse_grouped(post, grid, Sw, Swy, Swyy, s_post):
    """For every candidate tau: share = s_post * P(d < tau) + sum_{d >= tau} P(d) (c0 + c1 d) fitted by
    weighted least squares, using per-group sums of weights, weighted shares and squared shares
    (copies with the same substitution count and LTR length share one age posterior)."""
    cp = np.cumsum(post, 1)
    cpd = np.cumsum(post * grid[None, :], 1)
    A = np.concatenate([np.zeros((post.shape[0], 1)), cp[:, :-1]], 1)
    B = 1 - A
    Cd = cpd[:, -1:] - np.concatenate([np.zeros((post.shape[0], 1)), cpd[:, :-1]], 1)
    Tm = Swy[:, None] - s_post * A * Sw[:, None]
    a11, a12, a22 = (Sw[:, None] * B * B).sum(0), (Sw[:, None] * B * Cd).sum(0), (Sw[:, None] * Cd * Cd).sum(0)
    b1, b2 = (B * Tm).sum(0), (Cd * Tm).sum(0)
    syy = (Swyy[:, None] - 2 * s_post * A * Swy[:, None] + s_post ** 2 * A * A * Sw[:, None]).sum(0)
    det = a11 * a22 - a12 ** 2
    det = np.where(np.abs(det) < 1e-12, np.nan, det)
    c0 = (b1 * a22 - b2 * a12) / det
    c1 = (b2 * a11 - b1 * a12) / det
    sse = syy - 2 * c0 * b1 - 2 * c1 * b2 + c0 * c0 * a11 + 2 * c0 * c1 * a12 + c1 * c1 * a22
    sse[~np.isfinite(sse) | (c0 > s_post)] = np.inf
    return sse, c0, c1


def _tau_hat(sse, grid):
    """Step time from an SSE profile. The age distribution has few support points, so the fit is flat between
    them: the estimate is the flat minimum [lo, hi]; point = its geometric middle (its low end is biased young),
    or 0 when it reaches the bottom of the grid (below resolution). Returns (point, lo, hi, argmin)."""
    j = int(np.argmin(sse))
    flat = sse <= sse[j] + 1e-9 * abs(sse[j]) + 1e-12
    lo = hi = j
    while lo > 0 and flat[lo - 1]:
        lo -= 1
    while hi < sse.size - 1 and flat[hi + 1]:
        hi += 1
    mid = 0.0 if lo <= 1 else float(np.sqrt(grid[lo] * grid[hi]))
    return mid, float(grid[lo]), float(grid[hi]), j


_MB = {}                 # merger-bootstrap state shared with forked workers


def _merger_boot(nu):
    """One Poisson-bootstrap replicate of the pooled merger fit: (optimality gap, (flat-minimum low, high))."""
    b = _MB
    inv, G, w, y = b["inv"], b["G"], b["w"], b["y"]
    cb, swb, swyb, swyyb = (np.bincount(inv, nu, G), np.bincount(inv, nu * w, G), np.bincount(inv, nu * w * y, G),
                            np.bincount(inv, nu * w * y * y, G))
    _, pb, gpb, _ok = _age_em(b["P"], np.maximum(cb, 1e-12), pi=b["pi"], sqp=b["sqp"])
    return gpb, _tau_hat(_merger_sse_grouped(pb, b["grid"], swb, swyb, swyyb, b["s_null"])[0], b["grid"])[1:3]


def fit_merger_pooled(k2p, ksub, L, share, m, s_null, trough, rng, boot=50, threads=1):
    """Merger time from all copies younger than the confined trough at once, at the resolution of their
    pooled LTR sites. Each copy's age posterior comes from its substitution count and LTR length under
    the nonparametric age distribution of all such copies; copies inserted after the merger carry the
    no-confinement share (s_null), older copies a lower share that may drift with age. tau by weighted
    least squares; 95% CI and one-sided 95% bound ('upper', useful below resolution) by Poisson bootstrap
    over copies with the age distribution re-fitted each time (the SSE-profile bound undercovers: 72-76% in
    simulation). Copies of one burst are correlated; lineages cannot be resampled as clusters (one giant
    connected lineage), so intervals may be optimistic. Computed on groups of copies with identical
    (substitution count, LTR length), which is exact and fast; bootstrap re-fits run in parallel (weights drawn
    in order beforehand, so results do not depend on the number of threads)."""
    k = np.where(np.isfinite(ksub), ksub, np.round(np.nan_to_num(k2p) * np.nan_to_num(L, nan=500.0)))
    Lf = np.where(np.isfinite(L) & (L > 0), L, 500.0)
    ok = np.isfinite(k2p) & (k2p <= trough) & (Lf >= 50) & (m > 0) & np.isfinite(share) & (k >= 0)
    if ok.sum() < 100:
        return None
    k, Lf, y, w, kobs = k[ok], np.round(Lf[ok]), share[ok], m[ok], k2p[ok]
    top = max(float(np.quantile(k / Lf, 0.999)), 1e-3) * 1.5
    grid = np.concatenate([[0.0], np.geomspace(1e-6, top, 140)])
    key = np.unique(np.column_stack([k, Lf]), axis=0, return_inverse=True)
    ukl, inv = key[0], np.asarray(key[1]).ravel()
    P = _age_lik(ukl[:, 0], ukl[:, 1], grid)
    G = ukl.shape[0]

    def sums(nu):
        return (np.bincount(inv, nu, G), np.bincount(inv, nu * w, G), np.bincount(inv, nu * w * y, G),
                np.bincount(inv, nu * w * y * y, G))
    cnt, Sw, Swy, Swyy = sums(np.ones(k.size))
    step = grid[2] / grid[1]                   # one grid step (ratio)

    pi, post, gp0, sqp_ok = _age_em(P, cnt)
    gaps = [gp0]
    sse, c0, c1 = _merger_sse_grouped(post, grid, Sw, Swy, Swyy, s_null)
    tau, flat_lo, flat_hi, j = _tau_hat(sse, grid)
    tau_raw = tau
    nus = [rng.poisson(1.0, k.size).astype(float) for _ in range(boot)]
    _MB.update(P=P, inv=inv, G=G, w=w, y=y, grid=grid, s_null=s_null, pi=pi, sqp=sqp_ok)
    reps = []
    for gpb, rb in pool_map(_merger_boot, nus, threads if boot > 1 else 1):
        gaps.append(gpb)
        reps.append(rb)
    _MB.clear()
    if reps:                                    # percentile CI and one-sided bound from the ends of each
        reps = np.array(reps)                   # replicate's flat minimum (conservative), widened one grid step
        lo = np.quantile(reps[:, 0], 0.025) / step
        hi, upper = np.quantile(reps[:, 1], 0.975) * step, np.quantile(reps[:, 1], 0.95) * step
    else:                                       # no bootstrap: profile bound (undercovers; Hinkley 1970)
        sig2 = sse[j] / max(w.sum() - 3, 1) * w.mean()
        lo = hi = np.nan
        upper = float(grid[sse <= sse[j] + 3.84 * sig2].max())
    if tau > 0 and not lo > 0:                  # a date is reported only when its interval excludes zero
        tau = 0.0
    # fitted model in data space: each copy's expected share under its own age posterior (for plotting)
    Aj = post[:, :j].sum(1)
    pred = s_null * Aj + c0[j] * (1 - Aj) + c1[j] * (post[:, j:] * grid[None, j:]).sum(1)
    old_ = kobs >= tau_raw                         # observed share of copies older than the step (data, not model)
    sh_old = float(np.sum(w[old_] * y[old_]) / np.sum(w[old_])) if old_.any() and flat_lo > grid[1] else float("nan")
    return dict(tau=tau, tau_raw=tau_raw, flat_minimum=[flat_lo, flat_hi], ci=[float(lo), float(hi)], upper=float(upper),
                share_post=float(s_null), share_old_observed=sh_old,
                share_confined=float(c0[j]), share_slope=float(c1[j]), n=int(k.size),
                resolution=float(1.0 / np.median(Lf)), pred_k2p=kobs, pred=pred[inv], pred_w=w,
                max_optimality_gap=float(np.max(gaps)))


def k2p_se_guess(k2p, ltr_len=None):
    """K2P standard error when the table lacks one: binomial on the LTR length
    (default 500 bp), floored at one substitution."""
    L = np.where(np.isfinite(ltr_len), ltr_len, 500.0) if ltr_len is not None else 500.0
    p = np.maximum(np.nan_to_num(k2p, nan=0.0), 1.0 / L)
    return np.sqrt(p * (1 - p) / L)


# ---------------------------------------------------------------- genome scans (optional --genome)
def _iter_sequences(paths, want=None):
    for p in paths:
        for name, s in read_fasta(p):
            if want is None or name in want:
                yield name, s


A_PSEUDO = 0.5           # Jeffreys pseudo-count for empirical pattern counts


def confined_to(m):
    """For k-mers with per-subgenome counts m of OTHER chromosomes carrying them (the scored sequence's own
    chromosome left out): the subgenome they are confined to (present on >= 2 of its chromosomes and on
    no chromosome of any other subgenome), and a mask of confined k-mers. A k-mer of a lineage that
    amplified in one progenitor is confined to that subgenome wherever this copy sits; homoeologous
    single-copy k-mers (one other chromosome) and family-wide repeats are never confined."""
    g = m.argmax(1)
    tot = m.sum(1)
    mx = m[np.arange(m.shape[0]), g]
    return g, (mx >= 2) & (tot == mx)


PURITY_EDGES = np.array([0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0 + 1e-9])   # dominant-subgenome share of confined k-mers


NB_COUNT = 8            # log2 bins of confined k-mers per window (1, 2-3, ..., >= 128)


def _cells(X):
    """Composition cell of each window: (count bin, dominant subgenome, purity bin); -1 if empty."""
    n = X.sum(1)
    d = X.argmax(1)
    f = X[np.arange(X.shape[0]), d] / np.maximum(n, 1)
    fb = np.clip(np.searchsorted(PURITY_EDGES, f, side="right") - 1, 0, PURITY_EDGES.size - 2)
    nb = np.clip(np.floor(np.log2(np.maximum(n, 1))).astype(int), 0, NB_COUNT - 1)
    cell = (nb * X.shape[1] + d) * (PURITY_EDGES.size - 1) + fb
    return np.where(n > 0, cell, -1)


def _cond_logp(H, N):
    """log P(dominant, purity | count bin, origin) from cell counts H (origins x cells), Jeffreys-smoothed.
    Conditioning on the count bin: how many confined k-mers a window holds says nothing about its origin."""
    per = N * (PURITY_EDGES.size - 1)
    Hs = (H + 0.5).reshape(H.shape[0], NB_COUNT, per)
    return np.log(Hs / Hs.sum(2, keepdims=True)).reshape(H.shape[0], -1)


def paint_hmm(X_list, own, N, unit=None):
    """Ancestry along sequences from per-window confined compositions. Emission for origin g = how often
    windows of that composition (given their k-mer count) occur on g's chromosomes; for a chromosome's
    own subgenome the frequencies come from the OTHER chromosomes (cross-fitting), so its own
    exchanged windows cannot make themselves look typical (leave-one-chromosome-out, as in LOCO mixed-model
    GWAS, Listgarten et al. 2012; a near-identical partner is left out with it, unit[i]). One switch
    probability per window (Baum-Welch). own[i] = subgenome of sequence i (-1 unphased). Returns
    emissions, posteriors, Viterbi paths, switch rate."""
    own = np.asarray(own)
    ncell = NB_COUNT * N * (PURITY_EDGES.size - 1)
    cells = [_cells(X) for X in X_list]
    Hx = [np.bincount(c[c >= 0], minlength=ncell) for c in cells]
    H = np.zeros((N, ncell))
    for c, h, o in zip(cells, Hx, own):
        if o >= 0:
            H[o] += h
    unit = list(range(len(X_list))) if unit is None else list(unit)
    Hu = {}
    for h, o, u in zip(Hx, own, unit):
        if o >= 0:
            Hu[u] = Hu.get(u, 0) + h
    E = []
    for c, h, o, u in zip(cells, Hx, own, unit):
        Hc = H.copy()
        if o >= 0:
            Hc[o] -= Hu[u]                               # leave this chromosome (and its partner) out of its own origin
        lp = _cond_logp(Hc, N)
        E.append(np.where(c[:, None] >= 0, lp[:, np.maximum(c, 0)].T, 0.0))
    l0 = np.array([[math.log(0.99) if g == o else math.log(0.01 / max(N - 1, 1)) for g in range(N)]
                   if o >= 0 else [-math.log(N)] * N for o in own])
    ph = [i for i in range(len(X_list)) if own[i] >= 0]
    sw = estimate_switch([E[i] for i in ph], l0[ph])
    post, *_ = forward_backward(E, l0, sw)
    return E, post, viterbi(E, l0, sw), sw


def _scan_counts(args):
    """distinct sampled hashes of one sequence with their occurrence counts"""
    name, s, k, scale = args
    codes = LUT[np.frombuffer(s, np.uint8)]
    hs = []
    step = 1 << 24
    for st in range(0, codes.size, step):
        hs.append(sampled_hashes(codes[st:st + step + k - 1], k, scale)[1])
    h = np.concatenate(hs) if hs else np.empty(0, np.uint64)
    u, c = np.unique(h, return_counts=True)
    return name, u, c.astype(np.int32), int(codes.size)


def genome_counts(paths, sg_of, N, k, scale, threads, unit_of=None):
    """Pass 1: for every sampled k-mer carried by >= 2 phased chromosomes, the number of chromosomes of
    each subgenome carrying it (near-identical partners, unit_of, count once). Returns sorted hashes U,
    counts CNT (|U| x N), chromosomes per subgenome, sequence lengths."""
    unit_of = unit_of or {}
    H, G, Q, Cn, lengths = [], [], [], [], {}
    jobs = ((n, s, k, scale) for n, s in _iter_sequences(paths, set(sg_of)))
    for name, u, c, n in pool_map(_scan_counts, jobs, threads):
        H.append(u)
        Cn.append(c)
        G.append(np.full(u.size, sg_of[name], np.int16))
        Q.append(np.full(u.size, unit_of.get(name, len(lengths) + 1_000_000), np.int64))
        lengths[name] = n
        vlog(f"genome pass 1: {name} {n / 1e6:.1f} Mb")
    if not H:
        return np.empty(0, np.uint64), np.zeros((0, N), np.int32), np.zeros(N), lengths
    H = np.concatenate(H)
    G = np.concatenate(G)
    Q = np.concatenate(Q)
    Cn = np.concatenate(Cn)
    T = np.bincount(np.array([sg_of[nm] for nm in lengths]), minlength=N).astype(float)   # chromosomes per subgenome
    o = np.lexsort((Q, H))
    H, G, Q = H[o], G[o], Q[o]
    first = np.concatenate(([True], (H[1:] != H[:-1]) | (Q[1:] != Q[:-1])))   # one count per unit
    H, G = H[first], G[first]
    brk = np.flatnonzero(H[1:] != H[:-1]) + 1
    st = np.concatenate(([0], brk))
    sz = np.diff(np.concatenate((st, [H.size])))
    U = H[st]
    CNT = np.zeros((U.size, N), np.int32)
    np.add.at(CNT, (np.repeat(np.arange(U.size), sz), G), 1)          # chromosomes carrying the k-mer
    rep = CNT.sum(1) >= 2
    return U[rep], CNT[rep], T, lengths


def _scan_windows(args):
    """Pass 2: per window, the number of k-mers confined to each subgenome (confined_to) and their total."""
    name, s, k, scale, win, own = args
    U, CNT = _G["U"], _G["CNT"]
    N = CNT.shape[1]
    codes = LUT[np.frombuffer(s, np.uint8)]
    nw = codes.size // win + 1
    P, H = [], []
    step = 1 << 24
    for st in range(0, codes.size, step):
        p, h = sampled_hashes(codes[st:st + step + k - 1], k, scale)
        P.append(p + st)
        H.append(h)
    S = np.zeros((nw, N))
    n = np.zeros(nw, np.int64)
    if not P or U.size == 0:
        return name, S, n, int(codes.size)
    p, h = np.concatenate(P), np.concatenate(H)
    i = np.minimum(np.searchsorted(U, h), U.size - 1)
    hit = U[i] == h
    p, h, i = p[hit], h[hit], i[hit]
    m = CNT[i].astype(np.int64)
    if own >= 0:
        m[:, own] -= 1                              # this sequence's own presence
    g, ok = confined_to(m)
    w = p[ok] // win
    np.add.at(S, (w, g[ok]), 1.0)
    n = np.bincount(w, minlength=nw)
    return name, S, n, int(codes.size)


def genome_window_scores(paths, sg_of, U, CNT, k, scale, win, threads, min_len):
    """Window scores for every phased sequence and every other sequence >= min_len."""
    _G.update(U=U, CNT=CNT)
    out = {}
    jobs = ((nm, s, k, scale, win, sg_of.get(nm, -1)) for nm, s in _iter_sequences(paths)
            if nm in sg_of or len(s) >= min_len)
    for name, S, n, L in pool_map(_scan_windows, jobs, threads):
        out[name] = (S, n, L)
        vlog(f"genome pass 2: {name} {L / 1e6:.1f} Mb, {int(n.sum())} informative k-mers")
    return out


# ---------------------------------------------------------------- combined link matrix (v5)
NMAX_KNN = 50           # relatives for nearest-relative votes come from clades of <= 50 copies (as accurate as 250, 3.5x cheaper)
K_NEIGH = 3             # each copy votes for the 3 other chromosomes holding its closest relatives


def top_relatives(sp, ek, C, n_el, T, block_nnz=2e7):
    """For every copy and every other chromosome: its T closest relatives there (copy index, relatedness),
    relatedness = shared clades weighted by specificity (sum of mult/(size-1)); -1 = none.
    Computed once; half-libraries then use only relatives inside the same half."""
    S = sp.n
    A = sparse.csr_matrix((np.ones(sp.members.size, np.float32), sp.members, sp.indptr), shape=(S, n_el))
    d = (sp.mult / np.maximum(sp.size - 1, 1)).astype(np.float32)
    DA = sparse.csr_matrix((np.repeat(d, sp.size), sp.members, sp.indptr), shape=(S, n_el))
    AT = A.T.tocsr()
    cost = np.bincount(sp.members, weights=np.repeat(sp.size, sp.size).astype(float), minlength=n_el)
    cum = np.concatenate(([0], np.cumsum(cost)))
    rel_j = np.full((n_el, C, T), -1, np.int32)
    rel_r = np.zeros((n_el, C, T), np.float32)
    e0 = 0
    while e0 < n_el:
        e1 = int(np.searchsorted(cum, cum[e0] + block_nnz, side="right"))
        e1 = min(max(e1, e0 + 1), n_el)
        R = (AT[e0:e1] @ DA).tocoo()
        rows, cols, vals = R.row, R.col, R.data
        ch = ek[cols]
        k = ch != ek[e0 + rows]
        rows, cols, ch, vals = rows[k], cols[k], ch[k], vals[k]
        if rows.size:
            key = rows.astype(np.int64) * C + ch
            o = np.lexsort((-vals, key))
            key, cols, vals = key[o], cols[o], vals[o]
            st = np.searchsorted(key, key, side="left")
            rank = np.arange(key.size) - st
            m = rank < T
            r_, c_ = key[m] // C, key[m] % C
            rel_j[e0 + r_, c_, rank[m]] = cols[m]
            rel_r[e0 + r_, c_, rank[m]] = vals[m]
        e0 = e1
    return rel_j, rel_r


def knn_votes(rel_j, rel_r, ek, C, mask):
    """C x C symmetric nearest-relative votes: every copy in mask adds 1 to its chromosome x each of the
    K_NEIGH other chromosomes holding its closest relatives within mask. Each copy votes once, so
    large families cannot dominate."""
    valid = (rel_j >= 0) & mask[np.maximum(rel_j, 0)]
    best = np.where(valid, rel_r, 0).max(2)
    best[~mask] = 0
    k = min(K_NEIGH, C - 1)
    top = np.argpartition(-best, k - 1, axis=1)[:, :k] if k < C else np.argsort(-best, 1)[:, :k]
    tv = np.take_along_axis(best, top, 1)
    W = np.zeros((C, C))
    ok = tv > 0
    rows = np.repeat(np.arange(best.shape[0])[:, None], k, 1)
    np.add.at(W, (ek[rows[ok]], top[ok]), 1.0)
    W = W + W.T
    np.fill_diagonal(W, 0)
    return W


class LinkModel:
    """Everything needed to build chromosome link matrices for the whole library or any half of it:
    clades (carrier sets <= max(WINDOWS)) binned by size, and each copy's closest relatives."""

    def __init__(self, lib, ek, keep, C):
        self.ek, self.C, self.n = ek, C, lib.n
        self.sp = Splits(lib.kmer, lib.eid, ek, C, max(WINDOWS), keep)
        lo = 1
        self.bins = []
        for w in WINDOWS:
            self.bins.append((w, self.sp.view((self.sp.size > lo) & (self.sp.size <= w))))
            lo = w
        T = int(max(1, min(4, 5e7 // max(lib.n * C, 1))))
        self.rel_j, self.rel_r = top_relatives(self.sp.view(self.sp.size <= NMAX_KNN), ek, C, lib.n, T)

    def links_by_window(self, mask):
        """{window: clade-link matrix with clades of <= window carriers}, carriers restricted to mask."""
        out, acc = {}, np.zeros((self.C, self.C))
        wt = mask.astype(float)
        for w, v in self.bins:
            if v.n:
                acc = acc + affinity(v, self.ek, self.C, wt)
            out[w] = acc.copy()
        return out

    def matrix(self, mask, window, L=None):
        """Clade links (<= window) + nearest-relative votes rescaled to the same total (equal weight, no
        tuning). Benchmark (60 subsamples of 4 weak genomes): 138 -> 86 wrong chromosomes vs links alone."""
        L = self.links_by_window(mask)[window] if L is None else L
        K = knn_votes(self.rel_j, self.rel_r, self.ek, self.C, mask)
        if L.sum() > 0 and K.sum() > 0:
            return L + K * (L.sum() / K.sum())
        return L if L.sum() > 0 else K


def _half_links(seed):
    g = _G
    lm = g["lm"]
    u = np.random.default_rng(seed).random(lm.n) < 0.5
    out = []
    for m in (u & g["keep"], ~u & g["keep"]):
        Ls = lm.links_by_window(m)
        out.append((Ls, knn_votes(lm.rel_j, lm.rel_r, lm.ek, lm.C, m)))
    return out


def half_pairs(lm, keep, pairs, rng, threads):
    """For `pairs` complementary random halves of the library: link matrices per window and votes."""
    _G.update(lm=lm, keep=keep)
    return list(pool_map(_half_links, list(rng.integers(0, 2 ** 62, pairs)), threads))


def choose_window_halves(hp, C):
    """Carrier-count cap whose clade-link pattern (log O/E) is most reproducible between disjoint
    halves (mean Pearson r over the half pairs)."""
    iu = np.triu_indices(C, 1)
    scores = {}
    for w in WINDOWS:
        r = []
        for (La, _), (Lb, _) in hp:
            a, b = log_oe(La[w])[iu], log_oe(Lb[w])[iu]
            if a.std() > 0 and b.std() > 0:
                r.append(np.corrcoef(a, b)[0, 1])
        scores[w] = float(np.mean(r)) if r else -1.0
    return max(scores, key=lambda x: (round(scores[x], 3), x)), scores


def combine(L, K):
    if L.sum() > 0 and K.sum() > 0:
        return L + K * (L.sum() / K.sum())
    return L if L.sum() > 0 else K


def bipartition(W, idx, rng):
    """Best 2-way split of chromosomes idx from links among them only (null recomputed in the group)."""
    Wg = W[np.ix_(idx, idx)]
    if Wg.sum() <= 0 or len(idx) < 2:
        return np.zeros(len(idx), int)
    return best_partition(Wg - expected(Wg), 2, [[i] for i in range(len(idx))], rng, restarts=10)


def split_contrast(W, idx, z):
    """Mean log(O/E) within the two parts minus between them (null recomputed inside the group)."""
    Wg = W[np.ix_(idx, idx)]
    L = log_oe(Wg)
    same = z[:, None] == z[None, :]
    off = ~np.eye(len(idx), dtype=bool)
    if not (same & off).any() or not (~same).any():
        return 0.0
    return float(L[same & off].mean() - L[~same].mean())


def find_subgenomes(W, pairs, rng, min_size=2, max_n=12):
    """Divisive search for subgenome structure. A group of chromosomes is split in two when the
    best splits found independently in disjoint halves of the library agree (mean adjusted Rand
    index >= MIN_REP) and, below the top level, the split is at least MIN_REL as strong as the
    split that created the group (lineage structure inside a subgenome is weaker than the
    subgenome split itself). N=1 is a possible answer. Returns labels and every tested split
    (the structure tree: accepted and rejected levels are both reported)."""
    C = W.shape[0]
    groups = [(np.arange(C), None, "root")]
    final, tests = [], []
    while groups:
        g, parent, path = groups.pop(0)
        if len(g) < 2 * min_size or len(final) + len(groups) + 1 >= max_n:
            final.append(g)
            continue
        z = bipartition(W, g, rng)
        if min(np.bincount(z, minlength=2)) < min_size:
            tests.append(dict(members=g, split=z, rep=0.0, contrast=0.0, parent_contrast=parent, accepted=False, path=path))
            final.append(g)
            continue
        con = split_contrast(W, g, z)
        reps = [ari(bipartition(Wa, g, rng), bipartition(Wb, g, rng)) for Wa, Wb in pairs]
        rep = float(np.mean(reps))
        ok = rep >= MIN_REP and con > 0 and (parent is None or con >= MIN_REL * parent)
        tests.append(dict(members=g, split=z, rep=rep, contrast=con, parent_contrast=parent, accepted=ok, path=path))
        if ok:
            groups += [(g[z == 0], con, path + ".a"), (g[z == 1], con, path + ".b")]
        else:
            final.append(g)
    lab = np.zeros(C, int)
    for k, g in enumerate(sorted(final, key=lambda x: x.min())):
        lab[g] = k
    return lab, tests


def support_from_halves(halves, z, N, sets, rng):
    """Fraction of half-library partitions (each solved from scratch with the full-data estimator) agreeing
    with z, per chromosome."""
    hits = []
    for W in halves:
        if W.sum() <= 0:
            continue
        zr = best_partition(W - expected(W), N, sets, rng)       # same estimator as the full data, no warm start
        hits.append(align(z, zr, N) == z)
    return np.mean(hits, 0) if hits else np.full(len(z), np.nan)


# ---------------------------------------------------------------- ancestry HMM (shared by LTR-only and genome painting)
def _pad(E_list):
    n = max(e.shape[0] for e in E_list)
    N = E_list[0].shape[1]
    X = np.zeros((len(E_list), n, N))
    M = np.zeros((len(E_list), n), bool)
    for i, e in enumerate(E_list):
        X[i, :e.shape[0]] = e
        M[i, :e.shape[0]] = True
    return X, M


def _lT(N, sw):
    lT = np.full((N, N), math.log(sw / max(N - 1, 1)) if N > 1 else 0.0)
    np.fill_diagonal(lT, math.log1p(-sw) if N > 1 else 0.0)
    return lT


def forward_backward(E_list, l0, sw):
    """Posterior state probabilities for many sequences at once (padded; log space).
    E_list: per sequence (n_i x N) log emissions; l0: (S x N) log initial; sw: switch prob/step.
    Returns posteriors (list), expected switches, expected stays, total log-likelihood."""
    from scipy.special import logsumexp
    X, M = _pad(E_list)
    S, n, N = X.shape
    lT = _lT(N, sw)
    fw = np.empty((S, n, N))
    fw[:, 0] = l0 + X[:, 0]
    for t in range(1, n):
        fw[:, t] = logsumexp(fw[:, t - 1][:, :, None] + lT[None], 1) + X[:, t]
    bw = np.zeros((S, n, N))
    for t in range(n - 2, -1, -1):
        bw[:, t] = logsumexp(lT[None] + (X[:, t + 1] + bw[:, t + 1])[:, None, :], 2)
    ll = logsumexp(fw[np.arange(S), M.sum(1) - 1], 1)
    post = np.exp(fw + bw - ll[:, None, None])
    sw_e = st_e = 0.0
    for t in range(n - 1):
        m = M[:, t + 1]
        if not m.any():
            break
        x = fw[m, t][:, :, None] + lT[None] + (X[m, t + 1] + bw[m, t + 1])[:, None, :] - ll[m, None, None]
        x = np.exp(x).sum(0)
        st_e += np.trace(x)
        sw_e += x.sum() - np.trace(x)
    return [post[i, :M[i].sum()] for i in range(S)], sw_e, st_e, float(ll.sum())


def viterbi(E_list, l0, sw):
    X, M = _pad(E_list)
    S, n, N = X.shape
    lT = _lT(N, sw)
    vt = l0 + X[:, 0]
    bp = np.zeros((S, n, N), np.int8 if N < 128 else np.int16)
    for t in range(1, n):
        sc = vt[:, :, None] + lT[None]
        bp[:, t] = sc.argmax(1)
        vt = sc.max(1) + X[:, t]
    out = []
    for i in range(S):
        L = int(M[i].sum())
        path = np.empty(L, int)
        # padded steps carry zero emissions: the best end state at L-1 is the argmax of the running score
        # recomputed by back-tracking from the global end through padding
        st = int(vt[i].argmax())
        for t in range(n - 1, L - 1, -1):
            st = int(bp[i, t, st]) if t > 0 else st
        path[L - 1] = st
        for t in range(L - 1, 0, -1):
            path[t - 1] = bp[i, t, path[t]]
        out.append(path)
    return out


def estimate_switch(E_list, l0, sw=1e-3, iters=30):
    """One genome-wide switch probability per step by Baum-Welch (emissions fixed)."""
    for _ in range(iters):
        _, swe, ste, _ = forward_backward(E_list, l0, sw)
        new = float(np.clip(swe / max(swe + ste, 1e-300), 1e-7, 0.2))
        done = abs(new - sw) <= 0.01 * sw
        sw = new
        if done:
            break
    return sw


def calibrate(S, n, y):
    """Calibrated ancestry: P(o | unit) = softmax(lam * S / n^beta + b). S = summed k-mer scores,
    n = informative k-mers in the unit; beta in [0, 1] absorbs the correlation of k-mers within a
    unit (0: independent evidence, 1: one effective observation), lam is the temperature (Guo et
    al. 2017), b the subgenome biases; all by maximum likelihood against the units' chromosome
    labels y (exchanged or foreign units are rare label noise and only make the fit conservative)."""
    from scipy.optimize import minimize
    from scipy.special import logsumexp, expit
    N = S.shape[1]
    ln = np.log(np.maximum(n, 1.0))
    scale = float(np.median(np.abs(S).sum(1))) + 1e-12          # puts log-temperature near 0 at the optimum

    def unpack(par):
        return math.exp(min(par[0], 50.0)) / scale, float(expit(par[1])), np.concatenate(([0.0], par[2:]))

    def nll(par):
        lam, beta, b = unpack(par)
        Z = lam * S * np.exp(-beta * ln)[:, None] + b[None, :]
        return -(Z[np.arange(len(y)), y] - logsumexp(Z, 1)).sum()
    best = None
    for b0 in (-2.0, 0.0, 2.0):                     # convex in (temperature, biases) for fixed beta: start beta only
        r = minimize(nll, np.concatenate(([0.0, b0], np.zeros(N - 1))), method="L-BFGS-B")
        if best is None or r.fun < best.fun:
            best = r
    lam, beta, b = unpack(best.x)
    return lam, beta, b


def apply_calibration(S, n, cal):
    """log P(o | unit) under a calibration (lam, beta, b)."""
    from scipy.special import logsumexp
    lam, beta, b = cal
    Z = lam * S * np.exp(-beta * np.log(np.maximum(n, 1.0)))[:, None] + b[None, :]
    return Z - logsumexp(Z, 1, keepdims=True)


def _kmer_groups(lib, keep, ek, zc):
    """k-mer groups of the phased copies: sorted unique k-mers U, per-occurrence group index kid, and the
    occurrences' element ids / chromosomes / subgenomes."""
    sel = keep[lib.eid]
    km, ei = lib.kmer[sel], lib.eid[sel]
    cc = ek[ei]
    brk = np.flatnonzero(km[1:] != km[:-1]) + 1
    st = np.concatenate(([0], brk))
    sz = np.diff(np.concatenate((st, [km.size])))
    return km[st], np.repeat(np.arange(st.size), sz), ei, cc, zc[cc]


def _query(lib, keep, U):
    """Occurrences of unphased copies' k-mers among the phased copies' k-mer groups."""
    q = ~keep[lib.eid]
    qk, qe = lib.kmer[q], lib.eid[q]
    if not qk.size or not U.size:
        return np.empty(0, np.int64), np.empty(0, np.int64)
    j = np.minimum(np.searchsorted(U, qk), U.size - 1)
    hit = U[j] == qk
    return j[hit], qe[hit]


def copy_scores_abundance(lib, ek, keep, zc, N, unit=None):
    """Per copy: for each sampled k-mer, how many copies carry it in each subgenome on the OTHER
    chromosomes (own chromosome and its near-identical partner left out, so the copy and its local or
    homologous duplicates cannot vote), as a log share of that subgenome's copies (Jeffreys); summed over
    the copy's k-mers, centered over subgenomes. Captures lineage abundance. Returns S (E x N), n (E)."""
    E = lib.n
    S, n = np.zeros((E, N)), np.zeros(E)
    U, kid, ei, cc, gg = _kmer_groups(lib, keep, ek, zc)
    if not U.size:
        return S, n
    C = int(zc.size)
    unit = np.arange(C) if unit is None else unit
    Mg = np.zeros((U.size, N), np.int32)
    np.add.at(Mg, (kid, gg), 1)
    _, inv, kc = np.unique(kid.astype(np.int64) * C + unit[cc], return_inverse=True, return_counts=True)
    T = np.bincount(zc[ek[keep]], minlength=N).astype(float)
    for lo in range(0, kid.size, 1 << 24):
        hi = min(kid.size, lo + (1 << 24))
        m = Mg[kid[lo:hi]].astype(float)
        m[np.arange(hi - lo), gg[lo:hi]] -= kc[inv[lo:hi]]
        inf = m.sum(1) > 0
        sc = np.log((m[inf] + A_PSEUDO) / T[None, :])
        sc -= sc.mean(1, keepdims=True)
        np.add.at(S, ei[lo:hi][inf], sc)
        n += np.bincount(ei[lo:hi][inf], minlength=E)
    j, qe = _query(lib, keep, U)                          # unphased copies: query only
    if j.size:
        m = Mg[j].astype(float)
        inf = m.sum(1) > 0
        sc = np.log((m[inf] + A_PSEUDO) / T[None, :])
        sc -= sc.mean(1, keepdims=True)
        np.add.at(S, qe[inf], sc)
        n += np.bincount(qe[inf], minlength=E)
    return S, n


def copy_scores_breadth(lib, ek, keep, zc, N, unit=None):
    """Per copy: for each sampled k-mer, the rate at which each subgenome's OTHER chromosomes carry it
    (Jeffreys); if the copy came from origin o, its chromosome is one more draw from o, so log rate =
    log-likelihood of o. k-mers whose other carriers could all be homoeologs (<= 1 other chromosome per
    subgenome) say nothing about which subgenome this copy came from and are skipped, as are k-mers
    whose own evidence is below BF_CALL. Captures lineage breadth across chromosomes. Near-identical
    partner chromosomes (unit) count as one chromosome. Returns S, n."""
    E = lib.n
    S, n = np.zeros((E, N)), np.zeros(E)
    U, kid, ei, cc, gg = _kmer_groups(lib, keep, ek, zc)
    if not U.size:
        return S, n
    C = int(zc.size)
    unit = np.arange(C) if unit is None else unit
    uk = np.unique(kid.astype(np.int64) * C + unit[cc])
    Mg = np.zeros((U.size, N), np.int32)
    np.add.at(Mg, (uk // C, zc[uk % C]), 1)                         # chromosomes (units) per subgenome carrying k
    Cn = np.bincount(zc[np.unique(unit)], minlength=N).astype(float)

    def score(m, avail):
        sc = np.log((m + A_PSEUDO) / (avail + 2 * A_PSEUDO))
        ok = (m.max(1) >= 2) & ((sc.max(1) - sc.min(1)) >= math.log(BF_CALL))
        return sc - sc.mean(1, keepdims=True), ok
    for lo in range(0, kid.size, 1 << 24):
        hi = min(kid.size, lo + (1 << 24))
        m = Mg[kid[lo:hi]].astype(float)
        g = gg[lo:hi]
        m[np.arange(hi - lo), g] -= 1                                # leave the copy's own chromosome out
        sc, ok = score(m, Cn[None, :] - np.eye(N)[g])
        np.add.at(S, ei[lo:hi][ok], sc[ok])
        n += np.bincount(ei[lo:hi][ok], minlength=E)
    j, qe = _query(lib, keep, U)
    if j.size:
        sc, ok = score(Mg[j].astype(float), Cn[None, :])
        np.add.at(S, qe[ok], sc[ok])
        n += np.bincount(qe[ok], minlength=E)
    return S, n


def copy_ancestry(lib, ek, keep, zc, N, unit=None):
    """Calibrated per-copy ancestry log-posteriors from two complementary evidence scores (lineage
    abundance, lineage breadth), each temperature-calibrated against the copies' chromosome labels, then
    averaged and recalibrated (a product of experts; benchmark: the only per-copy model better than v4 in
    every test genome). Returns log P (E x N), number of informative k-mers (E), calibrations."""
    own = zc[ek]
    parts, ns, cals = [], [], []
    for fn in (copy_scores_abundance, copy_scores_breadth):
        S, n = fn(lib, ek, keep, zc, N, unit)
        fit = keep & (n > 0)
        cal = calibrate(S[fit], n[fit], own[fit])
        lp = apply_calibration(S, n, cal)
        lp[n == 0] = -math.log(N)
        parts.append(lp)
        ns.append(n)
        cals.append(cal)
    Z = np.mean(parts, 0)
    Z -= Z.mean(1, keepdims=True)
    ninf = np.maximum(ns[0], ns[1])
    fit = keep & (ninf > 0)
    cal = calibrate(Z[fit], np.ones(fit.sum()), own[fit])
    lP = apply_calibration(Z, np.ones(lib.n), cal)
    # posterior of a copy whose k-mers are spread over the subgenomes in proportion to their size (both
    # scores 0 by construction): what an unconfined (post-merger) copy is expected to carry
    z0 = np.mean([apply_calibration(np.zeros((1, N)), np.ones(1), c)[0] for c in cals], 0)
    lP0 = apply_calibration((z0 - z0.mean())[None, :], np.ones(1), cal)[0]
    return lP, ninf, dict(abundance=cals[0], breadth=cals[1], combined=cal, unconfined=lP0)


# ---------------------------------------------------------------- driver
HALF_PAIRS = 50          # complementary half-library pairs (100 halves): auto-N agreement, support, clade cap


def run_phasing(a, lib):
    """Chromosome-level structure: returns dict with chroms, ek, keep, z, N, support, tree, W."""
    rng = np.random.default_rng(a.seed)
    cnt = {}
    for c in lib.chrom:
        cnt[c] = cnt.get(c, 0) + 1
    me = getattr(a, "min_elements", MIN_ELEMENTS)
    chroms = [c for c in lib.chroms if cnt[c] >= me]
    cfg = read_config(a.config) if a.config else None
    N = a.n_subgenomes or (max(len(s) for s in cfg) if cfg else None)
    if cfg and not a.n_subgenomes and N < 2:
        die("config sets have one chromosome each: they define no cannot-link constraint (N would be 1)")
    if len(chroms) < max(N or 2, 2):
        die(f"only {len(chroms)} sequences carry >= {me} LTR-RTs; cannot form {max(N or 2, 2)} groups "
            f"(lower --min_elements for small libraries)")
    cix = {c: i for i, c in enumerate(chroms)}
    C = len(chroms)
    if cfg:
        miss = sorted({c for st_ in cfg for c in st_} - set(chroms))
        if miss:
            log(f"WARNING: {len(miss)} config chromosome(s) absent from the library or with < {me} LTR-RTs (ignored): "
                + ", ".join(miss[:10]) + (" ..." if len(miss) > 10 else ""))
    echrom = np.array([cix.get(c, -1) for c in lib.chrom], np.int32)
    keep = echrom >= 0
    ek = np.where(keep, echrom, 0)
    dropped = len(lib.chroms) - C
    log(f"{lib.n} LTR-RTs on {len(lib.chroms)} sequences; phasing {C} with >= {me} LTR-RTs"
        + (f" ({dropped} smaller sequences: ancestry from their LTR-RTs; painted too if >= {WIN // 1000} kb and a genome "
           f"is given)" if dropped else ""))
    lm = LinkModel(lib, ek, keep, C)
    if lm.sp.n == 0:
        die("no shared k-mers link two chromosomes; is this an LTR-RT library?")
    vlog(f"{lm.sp.n} clades; relatives computed")
    hp = half_pairs(lm, keep, HALF_PAIRS, rng, a.threads)
    window, wscores = choose_window_halves(hp, C)
    log("clade-size cap (split-half r): " + ", ".join(f"<={k}: {v:.3f}" for k, v in wscores.items()) + f" -> {window}")
    W = lm.matrix(keep, window)
    pairs = [(combine(La[window], Ka), combine(Lb[window], Kb)) for (La, Ka), (Lb, Kb) in hp]
    sp = lm.sp
    tree, init = [], None
    if N is None:
        lab, tests = find_subgenomes(W, pairs, rng)
        for t in tests:
            pa = [chroms[i] for i in t["members"][t["split"] == 0]]
            pb = [chroms[i] for i in t["members"][t["split"] == 1]]
            tree.append(dict(level=t["path"], part_a=pa, part_b=pb, replicability=t["rep"], strength=t["contrast"],
                             parent_strength=t.get("parent_contrast"), accepted=bool(t["accepted"])))
            log(f"split {len(pa)}|{len(pb)}: replicability {t['rep']:.2f}, strength {t['contrast']:.2f}"
                + (f" (parent {t['parent_contrast']:.2f})" if t.get("parent_contrast") is not None else "")
                + (" -> accepted" if t["accepted"] else " -> rejected"))
        N = int(lab.max()) + 1
        init = lab
        log(f"inferred number of subgenomes: {N}")
    nel = np.bincount(ek[keep], minlength=C)
    halves = [w for p in pairs for w in p]
    if N == 1:
        z0 = np.zeros(C, int)
        return dict(chroms=chroms, C=C, ek=ek, keep=keep, echrom=echrom, z=z0, N=1, support=np.full(C, np.nan),
                    tree=tree, W=W, window=window, wscores=wscores, n_el=nel, sp=sp,
                    own_oe=np.full(C, np.nan), oth_oe=np.full(C, np.nan), pval=np.full(C, np.nan),
                    partners=(pt0 := _partners_logged(W, halves, z0, chroms)), unit=_units(pt0), min_elements=me)
    sets = make_sets(chroms, cfg, N)
    z = best_partition(W - expected(W), N, sets, rng, init=init)
    z = canon_labels(z, nel)
    support = support_from_halves(halves, z, N, sets, rng)
    for g in range(N):
        log(f"SG{g + 1}: {int((z == g).sum())} sequences, {int(nel[z == g].sum()):,} LTR-RTs")
    weak = [chroms[i] for i in np.flatnonzero(support < SUPPORT_OK)]
    if weak:
        log(f"support < {SUPPORT_OK:g} for {len(weak)} sequence(s): " + ", ".join(weak[:12]) + (" ..." if len(weak) > 12 else ""))
    own_oe, oth_oe, pval = chrom_tests(W, z, N)
    return dict(chroms=chroms, C=C, ek=ek, keep=keep, echrom=echrom, z=z, N=N, support=support, tree=tree, W=W, window=window,
                wscores=wscores, n_el=nel, sp=sp, own_oe=own_oe, oth_oe=oth_oe, pval=pval,
                partners=(pt := _partners_logged(W, halves, z, chroms)), unit=_units(pt), min_elements=me)


def _units(pt):
    """Counting unit of each chromosome: itself, or the lower index of a reproducible near-identical
    partner pair (both then count as one chromosome wherever chromosomes are counted as evidence)."""
    u = np.arange(pt["partner"].size)
    for i, j in enumerate(pt["partner"]):
        if j >= 0 and pt["support"][i] >= SUPPORT_OK:
            u[i] = min(i, j)
    return u


def _partners_logged(W, halves, z, chroms):
    partner, psup, pex = partner_pairs(W, halves, z)
    rep = [(i, j) for i, j in enumerate(partner) if j > i and psup[i] >= SUPPORT_OK]
    if rep:
        ex = np.array([min(pex[i], pex[j]) for i, j in rep])
        log(f"{len(rep)} reproducible within-subgenome partner pairs (support >= {SUPPORT_OK:g}; lineage excess over the "
            f"subgenome median: median {np.median(ex):.1f}x, range {ex.min():.1f}-{ex.max():.1f}x)"
            + (": " + ", ".join(f"{chroms[i]}/{chroms[j]}" for i, j in rep[:6]) + (" ..." if len(rep) > 6 else "")))
    return dict(partner=partner, support=psup, excess=pex)


P_CALL = 0.95            # posterior needed to call a copy's ancestry (calibrated posteriors; 1 error in 20)


def lineage_labels(sp, members):
    """Lineage label of each copy in `members`: connected components of copies sharing any clade (k-mer
    carrier set of <= NMAX_KNN copies). Copies of one lineage are one piece of evidence."""
    members = np.asarray(members)
    lab = np.arange(members.size)
    if members.size <= 1 or sp.members.size == 0:
        return lab
    pos = -np.ones(int(max(sp.members.max(), members.max())) + 1, np.int64)
    pos[members] = np.arange(members.size)
    hit = (pos[sp.members] >= 0) & np.repeat(sp.size <= NMAX_KNN, sp.size)
    sids, loc = sp.sid[hit], pos[sp.members[hit]]
    parent = np.arange(members.size)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    o = np.argsort(sids, kind="stable")
    sids, loc = sids[o], loc[o]
    brk = np.flatnonzero(np.diff(sids)) + 1
    for grp in np.split(loc, brk):
        r0 = find(grp[0])
        for e in grp[1:]:
            r = find(e)
            if r != r0:
                parent[r] = r0
    return np.array([find(i) for i in range(members.size)])


def lineage_groups(sp, members):
    """Number of independent lineages among `members`."""
    return int(np.unique(lineage_labels(sp, members)).size) if len(members) else 0


def scaffold_ancestry(lib, ph, LR, n, keep):
    """Ancestry of each sequence that was not phased (too few LTR-RTs, e.g. unplaced scaffolds), from its
    copies' calibrated evidence with each independent lineage counted once (mean over its copies).
    Returns {name: (n_ltr, n_lineages, best subgenome, log-odds best vs next)}."""
    out = {}
    names = np.array(lib.chrom)
    for nm in sorted(set(names[~keep]), key=natural_key):
        idx = np.flatnonzero((names == nm) & (n > 0))
        if idx.size == 0:
            continue
        labs = lineage_labels(ph["sp"], idx)
        _, inv = np.unique(labs, return_inverse=True)
        tot = np.zeros(LR.shape[1])
        for g in range(LR.shape[1]):
            tot[g] = (np.bincount(inv, weights=LR[idx, g]) / np.bincount(inv)).sum()
        o = np.argsort(-tot)
        out[nm] = (int((names == nm).sum()), int(inv.max() + 1), int(o[0]), float(tot[o[0]] - tot[o[1]]) if tot.size > 1 else 0.0)
    return out


def run_copies(a, lib, ph):
    """Per-copy ancestry posteriors, calls and LTR-only ancestry segments."""
    N, zc, ek, keep = ph["N"], ph["z"], ph["ek"], ph["keep"]
    lP, n, cals = copy_ancestry(lib, ek, keep, zc, N, ph.get("unit"))
    own = zc[ek]
    P = np.exp(lP)
    cal = cals["combined"]
    # class frequencies of the calibration population: posterior / prior = scaled likelihood (HMM evidence),
    # and the posterior of a copy without ancestry signal (clock null). The calibration intercepts are not
    # priors: they also absorb systematic score shifts (e.g. a subgenome holding duplicated homologues).
    fitted = keep & (n > 0)
    prior = (np.bincount(own[fitted], minlength=N) + 0.5) / (fitted.sum() + 0.5 * N)
    best = P.argmax(1)
    bp = P[np.arange(lib.n), best]
    call = np.full(lib.n, "unresolved", dtype=object)
    call[n == 0] = "no_relatives"
    sure = (n > 0) & (bp >= P_CALL)
    call[keep & sure & (best == own)] = "own"
    call[keep & sure & (best != own)] = "foreign"
    call[~keep & sure] = "assigned_unphased"
    lin = np.where(sure, best, -1)
    scaf = scaffold_ancestry(lib, ph, lP - np.log(prior)[None, :], n, keep)
    log(f"copy ancestry (abundance x breadth): calls "
        + ", ".join(f"{v}={c}" for v, c in zip(*np.unique(call, return_counts=True))))
    return dict(P=P, n_inf=n, cal=cal, cals=cals, prior=prior, call=call, lin=lin, bp=bp, scaffolds=scaf)


def age_informativeness(k2p, fit):
    """Weight in [0, 1] of each copy as evidence about the origin of its neighbourhood, read off the
    fitted clock curve: 1 at ages whose copies are most confined to one subgenome, 0 at ages whose
    copies carry other subgenomes' ancestry as often as unconfined copies would (post-merger
    insertions, lineages older than the progenitor split). Copies without K2P get 1."""
    b = fit.get("bins") if fit else None
    s_null = fit.get("s_null", np.nan) if fit else np.nan
    if not b or len(b["k2p"]) < 3 or not np.isfinite(s_null):
        return None
    A, Y = np.array(b["k2p"]), np.array(b["share"])
    o = np.argsort(A)
    A, Y = A[o], Y[o]
    if s_null - Y.min() < 0.05:                       # no age-dependent confinement to exploit
        return None
    w = np.clip((s_null - np.interp(np.nan_to_num(k2p, nan=-1.0), A, Y)) / (s_null - Y.min()), 0.0, 1.0)
    w[~np.isfinite(k2p)] = 1.0
    return w


def run_segments(a, lib, ph, cp, ck):
    """LTR-only ancestry along chromosomes: HMM over copies in positional order; with K2P, each copy's
    evidence is tempered by the subgenome confinement of its age class. A candidate exchange needs a
    likelihood ratio >= BF_DECISIVE with each independent lineage counted once, from >= 2 lineages each
    favouring the origin >= BF_CALL-fold."""
    N, zc = ph["N"], ph["z"]
    lP, prior = np.log(np.maximum(cp["P"], 1e-300)), cp["prior"]
    w = np.ones(lib.n)
    aw = age_informativeness(ck["k2p"], ck["fit"]) if ck else None
    if aw is not None:
        w = w * aw
        log(f"LTR-only HMM: copies weighted by the subgenome confinement of their age class (median {np.median(aw):.2f})")
    order = {c: np.flatnonzero((ph["echrom"] == c))[np.argsort(lib.start[ph["echrom"] == c])] for c in range(ph["C"])}
    cs = [c for c in range(ph["C"]) if order[c].size]
    E = [np.where(cp["n_inf"][order[c]][:, None] > 0, (lP[order[c]] - np.log(prior)[None, :]) * w[order[c]][:, None], 0.0)
         for c in cs]
    l0 = np.array([[math.log(0.99) if g == zc[c] else math.log(0.01 / max(N - 1, 1)) for g in range(N)] for c in cs])
    sw = estimate_switch(E, l0)
    paths = viterbi(E, l0, sw)
    posts, *_ = forward_backward(E, l0, sw)
    state = np.full(lib.n, -1)
    spost = np.full(lib.n, np.nan)
    segs = []
    for c, path, e, pst in zip(cs, paths, E, posts):
        idx = order[c]
        state[idx] = path
        spost[idx] = pst[np.arange(idx.size), path]
        brk = np.flatnonzero(np.diff(path)) + 1
        for s0, e0 in zip(np.concatenate(([0], brk)), np.concatenate((brk, [idx.size]))):
            g = int(path[s0])
            ii = idx[s0:e0]
            d = e[s0:e0, g] - e[s0:e0, zc[c]]
            llr = float(d.sum())
            sup = ii[d >= math.log(BF_CALL)]                  # copies whose own evidence favours this origin >= 3-fold
            n_against = int((d <= -math.log(BF_CALL)).sum())     # ... and copies favouring the chromosome's own subgenome
            nlin = lineage_groups(ph["sp"], sup) if g != zc[c] else 0
            if g != zc[c]:                                   # each lineage counts once: mean evidence per lineage
                _, invl = np.unique(lineage_labels(ph["sp"], ii), return_inverse=True)
                llr_lin = float((np.bincount(invl, weights=d) / np.bincount(invl)).sum())
            else:
                llr_lin = 0.0
            segs.append(dict(chrom=ph["chroms"][c], start=int(lib.start[ii].min()), end=int(lib.end[ii].max()),
                             origin=g, chrom_sg=int(zc[c]), n_ltr=int(ii.size), n_support=int(sup.size), n_against=n_against,
                             n_lineages=nlin, llr=llr, llr_lineage=llr_lin, mean_post=float(np.nanmean(pst[s0:e0, g])),
                             flag=bool(g != zc[c] and llr_lin >= math.log(BF_DECISIVE) and nlin >= 2)))
    log(f"LTR-only ancestry HMM: switch rate {sw:.2e} per copy; {sum(x['flag'] for x in segs)} candidate exchange segments "
        f"(likelihood ratio >= {BF_DECISIVE:g}, >= 2 independent lineages)")
    cp.update(segs=segs, state=state, spost=spost, switch=sw, order=order, posts=dict(zip(cs, posts)))
    return cp


WIN = 50_000            # genome painting window (bp): >= hundreds of sampled k-mers even in repeat-poor genomes
MIN_SEG_WIN = 2         # a genome segment must span >= 2 windows (> 50 kb, longer than any LTR-RT: not one insertion)


def genome_scale(paths):
    """Genome sampling: ~200 M k-mers at most (memory), at least 1/16."""
    size = 0
    for p in paths:
        fai = p + ".fai"
        if os.path.exists(fai):
            with open(fai) as f:
                size += sum(int(l.split("\t")[1]) for l in f)
        else:
            size += os.path.getsize(p) * (4 if p.endswith(".gz") else 1)
    return max(16, int(math.ceil(size / 2e8)))


def run_genome(a, ph):
    """Ancestry painting of every genome window from genome-wide k-mers whose subgenome distribution
    (leave-own-sequence-out) is learned from the phased chromosomes; calibrated against chromosome
    labels; HMM along each sequence; unphased sequences painted too."""
    N, zc = ph["N"], ph["z"]
    sg_of = {ph["chroms"][i]: int(zc[i]) for i in range(ph["C"])}
    unit_of = {ph["chroms"][i]: int(u) for i, u in enumerate(ph.get("unit", np.arange(ph["C"])))}
    gscale = genome_scale(a.genome)
    import pickle
    cf = os.path.join(a.outdir, ".v5cache", f"genome_{cache_key(a.genome, a.k, gscale, WIN, WIN, sorted(sg_of.items()), sorted(unit_of.items()))}.pkl")
    if os.path.exists(cf):
        WS = pickle.load(open(cf, "rb"))
        log(f"reloaded genome windows from {cf}")
    else:
        U, CNT, T, lengths = genome_counts(a.genome, sg_of, N, a.k, gscale, a.threads, unit_of)
        found = [c for c in ph["chroms"] if c in lengths]
        if not found:
            log("WARNING: no genome sequence names match the LTR-RT chromosome names; genome painting skipped")
            return None
        log(f"genome: {len(found)}/{ph['C']} phased chromosomes found, 1/{gscale} of k-mers, {U.size} repeated k-mers")
        WS = genome_window_scores(a.genome, sg_of, U, CNT, a.k, gscale, WIN, a.threads, min_len=WIN)
        del U, CNT
        pickle.dump(WS, open(cf + ".tmp", "wb"))
        os.replace(cf + ".tmp", cf)
    names = list(WS)
    own = [sg_of.get(nm, -1) for nm in names]
    E, posts, paths, sw = paint_hmm([WS[nm][0] for nm in names], own, N,
                                        [unit_of.get(nm, -1 - i) for i, nm in enumerate(names)])
    alpha = np.array([np.concatenate([WS[c][0] for c in WS if sg_of.get(c, -1) == g]).sum(0) + 0.5 for g in range(N)])
    segs, seq = [], {}
    for nm, path, e, pst in zip(names, paths, E, posts):
        L = WS[nm][2]
        own = sg_of.get(nm, -1)
        brk = np.flatnonzero(np.diff(path)) + 1
        bp_by = np.zeros(N)
        for s0, e0 in zip(np.concatenate(([0], brk)), np.concatenate((brk, [path.size]))):
            g = int(path[s0])
            st_, en_ = s0 * WIN, min(e0 * WIN, L)
            bp_by[g] += en_ - st_
            if own >= 0 and g != own:
                llr = float((e[s0:e0, g] - e[s0:e0, own]).sum())
                mp = float(pst[s0:e0, g].mean())
                # decisive = the data favour the foreign origin 100:1 (likelihood ratio) AND so does the posterior,
                # which also weighs the genome's own switch rate; 'weak' = likelihood ratio alone
                ev = ("decisive" if mp >= 1 - 1 / BF_DECISIVE else "weak") \
                    if llr >= math.log(BF_DECISIVE) and e0 - s0 >= MIN_SEG_WIN else "none"
                segs.append(dict(chrom=nm, start=st_, end=en_, origin=g, chrom_sg=own, n_windows=int(e0 - s0),
                                 n_kmers=int(WS[nm][1][s0:e0].sum()),
                                 llr=llr, mean_post=mp, evidence=ev, flag=ev == "decisive"))
        tot = e.sum(0)
        seq[nm] = dict(length=L, phased=own >= 0, own=own, bp=bp_by, majority=int(bp_by.argmax()),
                       logodds=tot, path=path, post=pst)
    nflag = sum(x["flag"] for x in segs)
    nweak = sum(x["evidence"] == "weak" for x in segs)
    log("genome painting: confined k-mer composition per window (empirical emissions, Baum-Welch); mean composition "
        + "; ".join(f"SG{g + 1} chromosomes " + "/".join(f"{x:.2f}" for x in alpha[g] / alpha[g].sum()) for g in range(N))
        + f"; switch {sw:.2e}/window; "
        f"{nflag} candidate exchange segments (likelihood ratio >= {BF_DECISIVE:g} and mean posterior >= 0.99, "
        f">= {MIN_SEG_WIN * WIN // 1000} kb); {nweak} weak (likelihood ratio only)")
    return dict(segs=segs, seq=seq, cal=alpha, switch=sw, scale=gscale, WS=WS)


def run_clock(a, lib, ph, cp, rng):
    """Ages and the cross-subgenome transposition clock from per-copy ancestry posteriors."""
    km, bad = read_k2p(a.k2p, getattr(a, "k2p_col", None))
    k2p = np.full(lib.n, np.nan)
    k2se = np.full(lib.n, np.nan)
    ltrl = np.full(lib.n, np.nan)
    ksub = np.full(lib.n, np.nan)
    for i, id_ in enumerate(lib.ids):
        v = km.get(id_)
        if v is not None:
            k2p[i], k2se[i], ltrl[i], ksub[i] = v
    nk = int(np.isfinite(k2p).sum())
    log(f"K2P for {nk}/{lib.n} LTR-RTs" + (f" ({bad} unreadable rows)" if bad else ""))
    out = dict(k2p=k2p, k2se=k2se, ltrl=ltrl, ksub=ksub, fit=None, window=None)
    if nk == 0:
        log("WARNING: no K2P ids match FASTA ids; ages skipped")
        return out
    keep = ph["keep"]
    own = ph["z"][ph["ek"]]
    use = keep & (cp["n_inf"] > 0)
    share = np.where(use, 1 - cp["P"][np.arange(lib.n), own], np.nan)
    m = use.astype(float)
    se = np.where(np.isfinite(k2se) & (k2se > 0), k2se, k2p_se_guess(k2p, ltrl))
    fit = fit_clock(k2p, se, np.nan_to_num(share), m, rng, boot=CLOCK_BOOT)
    p0 = np.exp(cp["cals"]["unconfined"])
    s_null = float(np.median(1 - p0[own[use]])) if use.any() else np.nan
    if fit:
        fit["pooled"] = fit_merger_pooled(k2p, ksub, ltrl, share, m, s_null, fit.get("k2p_trough", np.inf), rng,
                                          boot=MERGER_BOOT, threads=a.threads)
        fit["s_null"] = s_null
    out["fit"] = fit
    out["share"], out["w"] = share, m                # per-copy data behind the clock (for plotting)
    # SubPhaser-style divergence-hybridization window: central 95% of own-ancestry copies' K2P, per subgenome
    win = {}
    for g in range(ph["N"]):
        v = k2p[(cp["call"] == "own") & (own == g) & np.isfinite(k2p)]
        if v.size >= 20:
            win[g] = (float(np.quantile(v, 0.025)), float(np.median(v)), float(np.quantile(v, 0.975)), int(v.size))
    out["window"] = win
    if fit and fit.get("divergence_supported"):
        # confined lineages amplify after the progenitors split, so most own-ancestry copies must be younger than
        # the split: a rise in sharing younger than the median own-ancestry copy of any subgenome is not the split
        td = fit["tau_divergence"]
        fit["split_consistent"] = bool(win and td > fit.get("k2p_trough", 0) and all(td >= v[1] for v in win.values()))
    if not fit or not fit.get("pooled"):
        log(f"too few dated LTR-RTs with ancestry evidence to fit the cross-subgenome transposition clock "
            f"({int((use & np.isfinite(k2p)).sum())} dated copies with informative k-mers; >= 200 needed, and >= 100 "
            f"younger than the confinement trough)")
    if fit and fit.get("pooled") and fit["pooled"]["max_optimality_gap"] > 1e-4:
        log(f"WARNING: the age-distribution fit was not certified optimal in every bootstrap replicate (largest "
            f"optimality gap {fit['pooled']['max_optimality_gap']:.1e}); treat the merger interval as approximate")
    if fit and fit.get("pooled"):
        pm = fit["pooled"]
        ky = lambda x: f"{x / (2 * a.mu) / 1e3:,.1f}"
        if pm["tau"] > 0:
            log(f"merger (pooled clock, {pm['n']} young copies): K2P {pm['tau']:.3g} = {ky(pm['tau'])} kyr "
                f"(95% CI {ky(pm['ci'][0])}-{ky(pm['ci'][1])} kyr, mu={a.mu:g}{'' if getattr(a, 'mu_given', True) else ' rice default'}); "
                f"a minimum age: LTR gene conversion, copy loss and a lag before the first cross-subgenome "
                f"transposition all make the clock young")
        else:
            log(f"cross-subgenome transposition began < {ky(pm['upper'])} kyr ago (95% bound; step below resolution, "
                f"mu={a.mu:g}{'' if getattr(a, 'mu_given', True) else ' rice default'}); the merger itself may be older "
                f"(the clock dates the onset of cross-subgenome transposition, a minimum merger age)")
    if fit and fit.get("divergence_supported"):
        my = lambda x: f"{x / (2 * a.mu) / 1e6:,.2f}"
        dci = fit["tau_divergence_ci"]
        if fit.get("split_consistent"):
            log(f"progenitor divergence (tentative, a minimum age): K2P {fit['tau_divergence']:.3g} = "
                f"{my(fit['tau_divergence'])} Myr (95% CI {my(dci[0])}-{my(dci[1])} Myr): beyond this age lineages are "
                f"shared between subgenomes again - the progenitor split, or the age limit of surviving LTR-RTs if the "
                f"split is older")
        else:
            log(f"sharing rises again beyond K2P {fit['tau_divergence']:.3g} ({my(fit['tau_divergence'])} Myr), but most "
                f"own-ancestry copies of a subgenome are older: not interpretable as the progenitor split (not reported)")
    return out


def cache_key(paths, *params):
    """Checksum of input identity (path, size, mtime) and parameters, for resumable caches."""
    h = hashlib.sha1()
    for p in paths:
        st = os.stat(p)
        h.update(f"{os.path.abspath(p)}|{st.st_size}|{int(st.st_mtime)}".encode())
    h.update(("|".join(map(str, params)) + "|" + __version__).encode())
    return h.hexdigest()[:16]


def load_library(a):
    """Sketch the LTR-RT library, or reload it from outdir/.v5cache when the inputs are unchanged."""
    cdir = os.path.join(a.outdir, ".v5cache")
    os.makedirs(cdir, exist_ok=True)
    f = os.path.join(cdir, f"library_{cache_key(a.ltr_fasta, a.k, a.scale)}.npz")
    if os.path.exists(f):
        z = np.load(f, allow_pickle=False)
        lib = Library.__new__(Library)
        lib.kmer, lib.eid = z["kmer"], z["eid"]
        lib.ids, lib.chrom, lib.cls = z["ids"].tolist(), z["chrom"].tolist(), z["cls"].tolist()
        lib.start, lib.end, lib.length = z["start"], z["end"], z["length"]
        lib.n = len(lib.ids)
        lib.chroms = sorted(set(lib.chrom), key=natural_key)
        log(f"reloaded sketch from {f}")
        return lib
    lib = Library(a.ltr_fasta, a.k, a.scale, a.threads)
    np.savez(f + ".tmp.npz", kmer=lib.kmer, eid=lib.eid, ids=np.array(lib.ids), chrom=np.array(lib.chrom),
             cls=np.array(lib.cls), start=lib.start, end=lib.end, length=lib.length)
    os.replace(f + ".tmp.npz", f)
    return lib


def auto_scale(paths, k):
    """FracMinHash scale: keep ~75 Mbp worth of k-mers (memory-bound), at least 1/4."""
    size = sum(os.path.getsize(p) * (4 if p.endswith(".gz") else 1) for p in paths)
    return max(4, int(math.ceil(size / 75e6)))


def main(argv=None):
    global VERBOSE
    ap = argparse.ArgumentParser(description="Subgenome phasing of an allopolyploid from its LTR-RT library (v5).",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--ltr_fasta", nargs="+", required=True, help="LTR-RT FASTA(.gz) file(s); headers chrom:start-end[#Class/Superfamily/Family]")
    ap.add_argument("--outdir", required=True, help="output directory")
    ap.add_argument("--k2p", nargs="+", default=None, help="optional LTR divergence table(s) (LTRquest/Kmer2LTR TSV with seq_id,k2p[,k2p_se,n_sites|ltr5_len,n_ts,n_tv])")
    ap.add_argument("--genome", nargs="+", default=None, help="optional genome FASTA(.gz) of the same assembly: paints ancestry along every sequence")
    ap.add_argument("-n", "--n_subgenomes", type=int, default=None, help="number of subgenomes (default: inferred)")
    ap.add_argument("--config", default=None, help="optional: sets of chromosomes known to be in different subgenomes (one set per line)")
    ap.add_argument("-t", "--threads", type=int, default=4, help="worker processes")
    ap.add_argument("--genome_fai", nargs="+", default=None,
                    help="FASTA index(es) of the assembly: chromosome lengths for tables and figures without --genome")
    ap.add_argument("--min_elements", type=int, default=MIN_ELEMENTS,
                    help=f"sequences with fewer LTR-RTs are not phased (default {MIN_ELEMENTS}; lower for small libraries)")
    ap.add_argument("--k2p_col", type=int, default=None,
                    help="1-based K2P column, ids in column 1 (overrides any header; default: the 'k2p' header column, "
                         "else the last column)")
    ap.add_argument("--mu", type=float, default=None,
                    help="LTR substitutions/site/year for K2P -> years (default 1.3e-8, a rice LTR rate; set it for your taxon)")
    ap.add_argument("--seed", type=int, default=1, help="random seed")
    ap.add_argument("--no_plots", action="store_true", help="skip figures")
    ap.add_argument("-v", "--verbose", action="store_true", help="per-step progress and sanity checks")
    a = ap.parse_args(argv)
    VERBOSE = a.verbose
    a.mu_given = a.mu is not None
    a.mu = a.mu if a.mu_given else MU_DEFAULT
    for p in a.ltr_fasta + (a.k2p or []) + (a.genome or []) + (a.genome_fai or []):
        if not os.path.exists(p):
            die(f"not found: {p}")
    a.lengths = {}                                   # sequence lengths: --genome_fai, else the genome's .fai if present
    for fai in (a.genome_fai or []) + [g + ".fai" for g in (a.genome or []) if os.path.exists(g + ".fai")]:
        with open(fai) as f:
            for l in f:
                q = l.split("\t")
                if len(q) > 1:
                    a.lengths[q[0]] = int(q[1])
    if a.n_subgenomes is not None and a.n_subgenomes < 2:
        die("-n must be >= 2 (omit it to infer the structure)")
    os.makedirs(a.outdir, exist_ok=True)
    a.k = K_DEFAULT
    a.scale = auto_scale(a.ltr_fasta, a.k)
    log(f"v{__version__}: {' '.join(a.ltr_fasta)} -> {a.outdir} (N={a.n_subgenomes or ('config' if a.config else 'auto')}, "
        f"k={a.k}, 1/{a.scale} of k-mers, {a.threads} threads)")
    lib = load_library(a)
    log(f"sketched {lib.n} LTR-RTs ({lib.kmer.size} sampled k-mer occurrences)")
    ph = run_phasing(a, lib)
    rng = np.random.default_rng(a.seed + 1)
    if ph["N"] == 1:
        write_null(a, lib, ph)
        log("done")
        return
    cp = run_copies(a, lib, ph)
    ck = run_clock(a, lib, ph, cp, rng) if a.k2p else None
    run_segments(a, lib, ph, cp, ck)
    gm = run_genome(a, ph) if a.genome else None
    if gm:
        cross_check(ph, cp, gm)
    cp["bio"] = run_biology(a, lib, ph, cp, ck, gm)
    write_outputs(a, lib, ph, cp, ck, gm)
    if not a.no_plots:
        try:
            make_figures(a, lib, ph, cp, ck, gm)
        except Exception as ex:          # figures never block the tables
            log(f"WARNING: plotting failed ({ex!r}); tables are complete")
    log("done")


# ---------------------------------------------------------------- outputs
SUPPORT_OK = 0.9         # support at which a chromosome counts as assigned (v4 benchmark: >= 97% correct)


def chrom_status(ph, gm, segs_ltr):
    """Per chromosome: status and painted ancestry fractions (genome painting if available)."""
    rows = []
    for i, c in enumerate(ph["chroms"]):
        g0 = int(ph["z"][i])
        frac = None
        nseg = 0
        flagged = [x for x in segs_ltr if x["chrom"] == c and x["flag"]]
        foreign = None
        if gm and c in gm["seq"]:
            bp = gm["seq"][c]["bp"]
            frac = bp / max(bp.sum(), 1)
            fg = [x for x in gm["segs"] if x["chrom"] == c and x["flag"]]
            foreign = sum(x["end"] - x["start"] for x in fg) / max(gm["seq"][c]["length"], 1)
            flagged += fg
        ends = []                                # candidate regions: LTR-only and genome calls, overlaps merged
        for x in sorted(flagged, key=lambda x: x["start"]):
            if ends and x["start"] < ends[-1]:
                ends[-1] = max(ends[-1], x["end"])
            else:
                ends.append(x["end"])
        nseg = len(ends)
        sup = ph["support"][i]
        if frac is not None and int(np.argmax(frac)) != g0:
            status = "conflict"
        elif nseg:
            status = "mosaic"
        elif np.isfinite(sup) and sup < SUPPORT_OK:
            status = "ambiguous"
        else:
            status = "assigned"
        rows.append(dict(status=status, frac=frac, n_seg=nseg, foreign=foreign))
    return rows


def run_biology(a, lib, ph, cp, ck, gm):
    """Per-genome TE biology from the same evidence (no new assumptions): cross-subgenome transposition by
    direction (copies called foreign, outside candidate exchange regions, split by age relative to the merger),
    each subgenome's insertion history (observed LTR divergence of the copies on its chromosomes), and a
    per-family table. Writes cross_subgenome.tsv, insertion_history.tsv, families.tsv;
    returns a summary dict."""
    od, N, keep = a.outdir, ph["N"], ph["keep"]
    sgn = [f"SG{g + 1}" for g in range(N)]
    host = np.where(keep, ph["z"][ph["ek"]], -1)
    call, lin = cp["call"], cp["lin"]
    # candidate exchange regions (LTR-only flags, decisive genome segments): foreign copies there moved with the
    # exchanged DNA, not by transposition
    regions = [(x["chrom"], x["start"], x["end"]) for x in cp["segs"] if x["flag"]]
    if gm:
        regions += [(x["chrom"], x["start"], x["end"]) for x in gm["segs"] if x["flag"]]
    in_ex = np.zeros(lib.n, bool)
    chrom = np.array(lib.chrom)
    for c_, s_, e_ in regions:
        in_ex |= (chrom == c_) & (lib.end > s_) & (lib.start < e_)
    k2p = ck["k2p"] if ck else np.full(lib.n, np.nan)
    pm = (ck["fit"] or {}).get("pooled") if ck else None
    young_cut = pm["tau"] if (pm and pm["tau"] > 0) else None    # only a dated merger (an onset bound is not one)
    own_n = np.array([int(((call == "own") & (host == g)).sum()) for g in range(N)])
    rows, summ = [], {}
    with open(os.path.join(od, "cross_subgenome.tsv"), "w") as f:
        f.write("donor_subgenome\thost_subgenome\tn_foreign\tn_in_exchange_regions\tn_transposed\t"
                "n_transposed_younger_than_merger\ttransposed_per_1000_donor_own_copies\n")
        for d in range(N):
            for h in range(N):
                if d == h:
                    continue
                fo = (call == "foreign") & (lin == d) & (host == h)
                tr = fo & ~in_ex
                yg = int((tr & (k2p <= young_cut)).sum()) if young_cut is not None else None
                rate = 1000 * tr.sum() / max(own_n[d], 1)
                f.write(f"{sgn[d]}\t{sgn[h]}\t{int(fo.sum())}\t{int((fo & in_ex).sum())}\t{int(tr.sum())}\t"
                        f"{'NA' if yg is None else yg}\t{rate:.2f}\n")
                rows.append(dict(donor=sgn[d], host=sgn[h], transposed=int(tr.sum()), younger_than_merger=yg,
                                 per_1000_donor_own=round(float(rate), 3)))
    summ["cross_subgenome_transposition"] = rows
    if rows:
        log("cross-subgenome transposition (copies called foreign outside exchange regions; a lower bound, calls need "
            "posterior >= 0.95): " + "; ".join(f"{r['donor']}->{r['host']} {r['transposed']}"
                                              + (f" ({r['younger_than_merger']} younger than the merger)"
                                                 if r["younger_than_merger"] is not None else "") for r in rows))
    # insertion history per subgenome: observed LTR divergence of the copies on its chromosomes. Not deconvolved:
    # deconvolving the whole age distribution is ill-posed (its shape depended on optimisation effort, total
    # variation up to 0.14 between runs); young ages are blurred by LTR length (about 1/L per substitution)
    if ck is not None and np.isfinite(k2p).any():
        ok = np.isfinite(k2p)
        top = max(float(np.nanquantile(k2p[ok], 0.995)), 1e-3) * 1.2
        edges = np.concatenate([[0.0], np.geomspace(2e-5, top, 40)])
        hist = {}
        with open(os.path.join(od, "insertion_history.tsv"), "w") as f:
            f.write("subgenome\tk2p_lo\tk2p_hi\tyears_lo\tyears_hi\tfraction_of_copies\tcopies\n")
            for g in range(N):
                v = k2p[ok & (host == g)]
                if v.size < 20:
                    continue
                cnt = np.histogram(np.clip(v, 0, edges[-1]), edges)[0]
                for lo, hi, c_ in zip(edges[:-1], edges[1:], cnt):
                    f.write(f"{sgn[g]}\t{lo:.3g}\t{hi:.3g}\t{lo / (2 * a.mu):.0f}\t{hi / (2 * a.mu):.0f}\t"
                            f"{c_ / v.size:.4f}\t{int(c_)}\n")
                q = {f"q{int(100 * x)}_years": float(np.quantile(v, x) / (2 * a.mu)) for x in (0.1, 0.5, 0.9)}
                hist[sgn[g]] = dict(n_copies=int(v.size), **q,
                                    fraction_younger_than_merger=(float((v <= young_cut).mean())
                                                                  if young_cut is not None else None))
        summ["insertion_history"] = hist
        if hist:
            log("insertion history (observed LTR divergence, copies on each subgenome's chromosomes): "
                + "; ".join(f"{g_}: median {v['q50_years'] / 1e6:,.2f} Myr (10-90%: {v['q10_years'] / 1e6:,.2f}-"
                            f"{v['q90_years'] / 1e6:,.2f})" for g_, v in hist.items()))
    # families
    fam = np.array([("/".join(c_.split("/")[1:3]) if c_.count("/") >= 2 else c_) for c_ in lib.cls])
    tot = np.array([int((host == g).sum()) for g in range(N)], float)
    with open(os.path.join(od, "families.tsv"), "w") as f:
        f.write("family\tn_copies\t" + "\t".join(f"n_{x}" for x in sgn) + "\t" + "\t".join(f"enrichment_{x}" for x in sgn)
                + "\tn_own\tn_foreign\tn_unresolved\tmedian_k2p\n")
        for fm in sorted(set(fam[keep]), key=lambda x: (-int(((fam == x) & keep).sum()), x)):   # deterministic
            m = (fam == fm) & keep
            ng = np.array([int((m & (host == g)).sum()) for g in range(N)], float)
            enr = (ng / max(ng.sum(), 1)) / np.maximum(tot / max(tot.sum(), 1), 1e-12)
            kk = k2p[m & np.isfinite(k2p)]
            f.write(f"{fm}\t{int(m.sum())}\t" + "\t".join(str(int(x)) for x in ng) + "\t"
                    + "\t".join(f"{x:.2f}" for x in enr) + f"\t{int((m & (call == 'own')).sum())}\t"
                    f"{int((m & (call == 'foreign')).sum())}\t{int((m & (call == 'unresolved')).sum())}\t"
                    f"{(f'{np.median(kk):.4g}' if kk.size else 'NA')}\n")
    return summ


def cross_check(ph, cp, gm):
    """Agreement between the two independent evidence tracks (LTR-RT lineages; genome-wide k-mers)."""
    both = [i for i, c in enumerate(ph["chroms"]) if c in gm["seq"] and gm["seq"][c]["bp"].sum() > 0]
    agree = sum(int(gm["seq"][ph["chroms"][i]]["majority"] == ph["z"][i]) for i in both)
    log(f"genome painting majority agrees with the LTR-RT assignment for {agree}/{len(both)} phased chromosomes")
    extra = [nm for nm, d in gm["seq"].items() if d["own"] < 0 and (gm["WS"][nm][1] > 0).any()]
    if extra:
        log(f"{len(extra)} unphased sequence(s) >= {WIN // 1000} kb painted from the genome (genome_sequences.tsv)")
    fl = [x for x in cp["segs"] if x["flag"]]
    if fl:
        sup = [genome_support(x, gm) for x in fl]
        log(f"LTR-only candidate exchanges overlapping a genome segment of the same origin: "
            f"{sup.count('decisive')} decisive, {sup.count('weak')} weak, {sup.count('none')} none (of {len(fl)})")


def genome_support(x, gm):
    """Independent check of an LTR-only segment by genome painting: the strongest evidence tier ('decisive',
    'weak') of an overlapping genome segment of the same origin, 'none' if there is none, NA without --genome."""
    if not gm or x["chrom"] not in gm["seq"] or x["origin"] == x["chrom_sg"]:
        return "NA"
    tiers = [y["evidence"] for y in gm["segs"] if y["chrom"] == x["chrom"] and y["origin"] == x["origin"]
             and y["start"] < x["end"] and y["end"] > x["start"] and y["evidence"] != "none"]
    return "decisive" if "decisive" in tiers else "weak" if tiers else "none"


def _f(x, fmt="{:.3f}"):
    return "NA" if x is None or (isinstance(x, float) and not np.isfinite(x)) else fmt.format(x)


def write_outputs(a, lib, ph, cp, ck, gm):
    od = a.outdir
    N = ph["N"]
    sgn = [f"SG{g + 1}" for g in range(N)]
    st = chrom_status(ph, gm, cp["segs"])
    ph["status"] = st
    with open(os.path.join(od, "chromosomes.tsv"), "w") as f:
        f.write("chrom\tlength\tsubgenome\tstatus\tsupport\tpvalue\toe_own\toe_other\tn_ltr\tn_own\tn_foreign\t"
                "n_unresolved\tn_no_relatives\tn_exchange_segments\tgenome_exchange_fraction\tpartner\tpartner_support\t"
                "partner_excess\t" + "\t".join(f"painted_{x}" for x in sgn) + "\n")
        for i, c in enumerate(ph["chroms"]):
            mi = ph["echrom"] == i
            cc = cp["call"][mi]
            fr = st[i]["frac"]
            ln = gm["seq"][c]["length"] if gm and c in gm["seq"] else a.lengths.get(c) if getattr(a, "lengths", None) else None
            f.write(f"{c}\t{ln if ln else 'NA'}\t{sgn[ph['z'][i]]}\t{st[i]['status']}\t{_f(ph['support'][i])}\t"
                    f"{_f(ph['pval'][i], '{:.3g}')}\t{_f(ph['own_oe'][i])}\t{_f(ph['oth_oe'][i])}\t{int(mi.sum())}\t"
                    f"{int((cc == 'own').sum())}\t{int((cc == 'foreign').sum())}\t{int((cc == 'unresolved').sum())}\t"
                    f"{int((cc == 'no_relatives').sum())}\t{st[i]['n_seg']}\t{_f(st[i]['foreign'])}\t{_partner_cols(ph, i)}\t"
                    + "\t".join(_f(x) if fr is not None else "NA" for x in (fr if fr is not None else [None] * N)) + "\n")
    k2p = ck["k2p"] if ck else np.full(lib.n, np.nan)
    cix = {c: i for i, c in enumerate(ph["chroms"])}
    with open(os.path.join(od, "elements.tsv"), "w") as f:
        f.write("id\tchrom\tstart\tend\tclass\tchrom_subgenome\tcall\tancestry_subgenome\tposterior\t"
                + "\t".join(f"P_{x}" for x in sgn) + "\tn_informative_kmers\tsegment_origin\tsegment_posterior\tk2p\tage_years\n")
        for i in range(lib.n):
            ci = cix.get(lib.chrom[i])
            csg = sgn[ph["z"][ci]] if ci is not None else "NA"
            anc = sgn[cp["lin"][i]] if cp["lin"][i] >= 0 else "NA"
            P = cp["P"][i]
            so = sgn[cp["state"][i]] if cp["state"][i] >= 0 else "NA"
            kv = f"{k2p[i]:.5g}" if np.isfinite(k2p[i]) else "NA"
            age = f"{k2p[i] / (2 * a.mu):.0f}" if np.isfinite(k2p[i]) else "NA"
            noinf = cp["n_inf"][i] == 0                  # no informative k-mers: no posterior (not evidence)
            f.write(f"{lib.ids[i]}\t{lib.chrom[i]}\t{lib.start[i]}\t{lib.end[i]}\t{lib.cls[i]}\t{csg}\t{cp['call'][i]}\t{anc}\t"
                    f"{'NA' if noinf else _f(cp['bp'][i])}\t" + "\t".join('NA' if noinf else _f(x) for x in P)
                    + f"\t{int(cp['n_inf'][i])}\t{so}\t{_f(cp['spost'][i])}\t{kv}\t{age}\n")
    with open(os.path.join(od, "segments.tsv"), "w") as f:
        f.write("chrom\tstart\tend\torigin_subgenome\tchrom_subgenome\tn_ltr\tn_support\tn_against\tn_lineages\tllr\t"
                "llr_per_lineage\tmean_posterior\tcandidate_exchange\tgenome_support\n")
        for x in cp["segs"]:
            f.write(f"{x['chrom']}\t{x['start']}\t{x['end']}\t{sgn[x['origin']]}\t{sgn[x['chrom_sg']]}\t{x['n_ltr']}\t"
                    f"{x['n_support']}\t{x['n_against']}\t{x['n_lineages']}\t{x['llr']:.1f}\t{x['llr_lineage']:.1f}\t"
                    f"{x['mean_post']:.3f}\t{'yes' if x['flag'] else 'no'}\t{genome_support(x, gm)}\n")
    if cp["scaffolds"]:
        with open(os.path.join(od, "unphased_sequences.tsv"), "w") as f:
            f.write("sequence\tn_ltr\tn_lineages\tancestry\tlog_odds_vs_next\tconfidence\n")
            for nm, (nl, nlin, g, lo) in cp["scaffolds"].items():
                conf = "decisive" if lo >= math.log(BF_DECISIVE) and nlin >= 2 else ("positive" if lo >= math.log(BF_CALL) else "none")
                f.write(f"{nm}\t{nl}\t{nlin}\t{sgn[g] if conf != 'none' else 'NA'}\t{lo:.2f}\t{conf}\n")
    write_links_oe(od, ph)
    if gm:
        with open(os.path.join(od, "segments_genome.tsv"), "w") as f:
            f.write("chrom\tstart\tend\torigin_subgenome\tchrom_subgenome\tn_windows\tn_kmers\tllr\tmean_posterior\t"
                    "evidence\tcandidate_exchange\n")
            for x in gm["segs"]:
                f.write(f"{x['chrom']}\t{x['start']}\t{x['end']}\t{sgn[x['origin']]}\t{sgn[x['chrom_sg']]}\t{x['n_windows']}\t"
                        f"{x['n_kmers']}\t{x['llr']:.1f}\t{x['mean_post']:.3f}\t{x['evidence']}\t{'yes' if x['flag'] else 'no'}\n")
        with open(os.path.join(od, "genome_sequences.tsv"), "w") as f:
            f.write("sequence\tlength\tltr_subgenome\tn_windows\tn_informative_windows\tpainted_majority\t"
                    + "\t".join(f"bp_{x}" for x in sgn) + "\tlog_odds_majority_vs_next\n")
            for nm, d in gm["seq"].items():
                lo = np.sort(d["logodds"])[::-1]
                ninf = int((gm["WS"][nm][1] > 0).sum())
                none_ = ninf == 0                          # no confined k-mer anywhere: no call
                f.write(f"{nm}\t{d['length']}\t{sgn[d['own']] if d['own'] >= 0 else 'NA'}\t{d['path'].size}\t{ninf}\t"
                        f"{'NA' if none_ else sgn[d['majority']]}\t"
                        + "\t".join("NA" if none_ else str(int(x)) for x in d["bp"])
                        + f"\t{'NA' if none_ else format((lo[0] - lo[1]) if N > 1 else 0, '.1f')}\n")
        with open(os.path.join(od, "genome_windows.tsv"), "w") as f:
            f.write("sequence\tstart\tend\tn_kmers\t" + "\t".join(f"kmers_{x}" for x in sgn) + "\tstate\t"
                    + "\t".join(f"P_{x}" for x in sgn) + "\n")
            for nm, d in gm["seq"].items():
                S_, n = gm["WS"][nm][0], gm["WS"][nm][1]
                none_ = not (n > 0).any()                    # no information on the whole sequence: no state
                for w in range(d["path"].size):
                    f.write(f"{nm}\t{w * WIN}\t{min((w + 1) * WIN, d['length'])}\t{int(n[w])}\t"
                            + "\t".join(str(int(x)) for x in S_[w]) + f"\t{'NA' if none_ else sgn[d['path'][w]]}\t"
                            + "\t".join("NA" if none_ else f"{x:.3f}" for x in d["post"][w]) + "\n")
    summ = dict(version=__version__, ltr_fasta=a.ltr_fasta, k=a.k, scale=a.scale, n_ltr=lib.n, n_phased_sequences=ph["C"],
                n_subgenomes=N, n_subgenomes_source=("given" if a.n_subgenomes else "config" if a.config else "inferred"),
                structure_tree=ph["tree"], clade_cap=ph["window"], clade_cap_reproducibility=ph["wscores"],
                chromosome_status={k_: int(v) for k_, v in zip(*np.unique([x["status"] for x in st], return_counts=True))},
                min_support=float(np.nanmin(ph["support"])),
                copy_calibration={k_: dict(temperature=v[0], evidence_exponent=v[1], intercepts=list(map(float, v[2])))
                                  for k_, v in cp["cals"].items() if k_ != "unconfined"},
                n_clades_total=int(ph["sp"].n), n_clades_within_cap=int((ph["sp"].size <= ph["window"]).sum()),
                subgenome_sizes={sgn[g]: dict(n_sequences=int((ph["z"] == g).sum()),
                                                                        n_ltr=int(ph["n_el"][ph["z"] == g].sum()))
                                                           for g in range(N)},
                mean_oe_within=float(np.nanmean(ph["own_oe"])), mean_oe_between=float(np.nanmean(ph["oth_oe"])),
                unconfined_copy_posterior=np.exp(cp["cals"]["unconfined"]).tolist(), copy_class_frequencies=cp["prior"].tolist(),
                n_unphased_assigned=sum(1 for v in cp["scaffolds"].values() if v[3] >= math.log(BF_CALL)),
                copy_calls={str(v): int(c) for v, c in zip(*np.unique(cp["call"], return_counts=True))},
                ltr_switch_rate=cp["switch"], n_candidate_exchanges_ltr=int(sum(x["flag"] for x in cp["segs"])), mu=a.mu,
                mu_source="user" if getattr(a, "mu_given", True) else "default (rice LTR rate, Ma & Bennetzen 2004)",
                partner_pairs=_partner_pairs_json(ph), **cp.get("bio", {}))
    if gm:
        summ.update(genome_scale=gm["scale"], genome_switch_rate=gm["switch"],
                    genome_composition_alpha=gm["cal"].tolist(),
                    n_candidate_exchanges_genome=int(sum(x["flag"] for x in gm["segs"])))
    if ck and ck["fit"]:
        fit = ck["fit"]
        yrs = lambda x: None if x is None else x / (2 * a.mu)
        pm = fit.get("pooled")
        if pm:
            summ["merger_years"] = dict(estimate=yrs(pm["tau"]) if pm["tau"] > 0 else None, ci95=[yrs(x) for x in pm["ci"]],
                                        upper95=yrs(pm["upper"]), resolved=bool(pm["tau"] > 0), n_copies=pm["n"],
                                        method="pooled cross-subgenome transposition clock (minimum merger age)",
                                        note=("dated: the onset of cross-subgenome transposition, a minimum merger age"
                                              if pm["tau"] > 0 else "not dated: upper95 bounds the onset of "
                                              "cross-subgenome transposition, not the merger (which can be older)"))
        summ["transposition_clock"] = dict(
            n_dated_copies=fit.get("n"), k2p_confinement_trough=fit.get("k2p_trough"), k2p_max=fit.get("k2p_max"),
            share_unconfined=fit.get("s_null"),
            pooled_merger=None if not pm else dict(
                note="the reported merger estimate; tau is 0 (unresolved) unless its interval excludes 0",
                tau_k2p=pm["tau"], tau_k2p_before_resolution_rule=pm["tau_raw"], flat_minimum_k2p=pm["flat_minimum"],
                ci_k2p=pm["ci"], upper_k2p=pm["upper"], n_young_copies=pm["n"], resolution_k2p=pm["resolution"],
                share_old_copies_observed=pm["share_old_observed"],
                max_optimality_gap=pm["max_optimality_gap"]),
            binned_v4_method=dict(
                note="coarse binned step fit inherited from v4 (merger shown for comparison; split step used for the "
                     "tentative progenitor split)",
                tau_merger_k2p=fit.get("tau_merger"), tau_merger_ci=fit.get("tau_merger_ci"),
                merger_step_ci=fit.get("step_merger_ci"), share_post=fit.get("share_post"),
                share_confined=fit.get("share_confined"), tau_split_k2p=fit.get("tau_divergence"),
                tau_split_ci=fit.get("tau_divergence_ci"), split_step_ci=fit.get("step_divergence_ci"),
                split_step_supported=fit.get("divergence_supported"), split_consistent=fit.get("split_consistent"),
                share_ancestral=fit.get("share_ancestral")))
        if fit.get("divergence_supported"):
            summ["progenitor_split_years_tentative"] = dict(
                estimate=yrs(fit["tau_divergence"]), ci95=[yrs(x) for x in fit["tau_divergence_ci"]],
                k2p=fit["tau_divergence"], minimum_age=bool(fit.get("split_consistent")),
                consistent=bool(fit.get("split_consistent")),
                consistency_rule="older than the confinement trough and than the median own-ancestry copy of every subgenome; "
                                 "only consistent estimates are logged and drawn",
                note="age beyond which LTR-RT lineages are shared between subgenomes again: the progenitor split, or "
                     "the age limit of surviving LTR-RTs (copy loss, saturation, homoplasy) if the split is older")
        summ["own_ancestry_k2p_window"] = {sgn[g]: dict(p2_5=v[0], median=v[1], p97_5=v[2], n=v[3],
                                                        years=[yrs(v[0]), yrs(v[2])]) for g, v in ck["window"].items()}
    with open(os.path.join(od, "summary.json"), "w") as f:
        json.dump(summ, f, indent=2, default=float)


def _partner_pairs_json(ph):
    """reproducible within-subgenome partner pairs (support >= SUPPORT_OK)"""
    pt = ph.get("partners")
    if not pt:
        return []
    return [dict(a=ph["chroms"][i], b=ph["chroms"][j], support=float(pt["support"][i]),
                 excess=float(min(pt["excess"][i], pt["excess"][j])))
            for i, j in enumerate(pt["partner"]) if j > i and pt["support"][i] >= SUPPORT_OK]


def _partner_cols(ph, i):
    pt = ph.get("partners")
    if not pt or pt["partner"][i] < 0:
        return "NA\tNA\tNA"
    return f"{ph['chroms'][pt['partner'][i]]}\t{_f(pt['support'][i], '{:.2f}')}\t{_f(pt['excess'][i], '{:.2f}')}"


def write_links_oe(od, ph):
    """Chromosome x chromosome shared lineages, observed / expected under the degree-corrected null."""
    W = ph["W"]
    E = expected(W)
    OE = W / np.where(E > 0, E, np.nan)
    with open(os.path.join(od, "links_oe.tsv"), "w") as f:
        f.write("chrom\t" + "\t".join(ph["chroms"]) + "\n")
        for i, c in enumerate(ph["chroms"]):
            f.write(c + "\t" + "\t".join("NA" if not np.isfinite(x) else f"{x:.4f}" for x in OE[i]) + "\n")


def write_null(a, lib, ph):
    """No reproducible subgenome structure: every chromosome in SG1, plus the best (rejected) two-way
    split as a tentative grouping with its replicability, so a user can judge it or force -n."""
    log("no reproducible subgenome structure (diploid-like, autopolyploid, or too few informative LTR-RTs); "
        "the best split is reported as tentative; use -n to force a partition")
    root = next((t for t in ph["tree"] if t["level"] == "root"), None)
    tent = {}
    if root:
        tent.update({c: "a" for c in root["part_a"]})
        tent.update({c: "b" for c in root["part_b"]})
    with open(os.path.join(a.outdir, "chromosomes.tsv"), "w") as f:
        f.write("chrom\tlength\tsubgenome\tstatus\ttentative_group\tn_ltr\tpartner\tpartner_support\tpartner_excess\n")
        for i, c in enumerate(ph["chroms"]):
            ln = getattr(a, "lengths", {}).get(c)
            f.write(f"{c}\t{ln if ln else 'NA'}\tSG1\tno_structure\t{tent.get(c, 'NA')}\t{int((ph['echrom'] == i).sum())}\t"
                    f"{_partner_cols(ph, i)}\n")
    write_links_oe(a.outdir, ph)
    with open(os.path.join(a.outdir, "summary.json"), "w") as f:
        json.dump(dict(version=__version__, ltr_fasta=a.ltr_fasta, n_subgenomes=1,
                       n_subgenomes_source=("given" if a.n_subgenomes else "config" if a.config else "inferred"),
                       structure_tree=ph["tree"], n_ltr=lib.n, clade_cap=ph["window"],
                       n_phased_sequences=ph["C"],
                       tentative_split=(None if not root else dict(replicability=root["replicability"],
                                                                   strength=root["strength"])),
                       partner_pairs=_partner_pairs_json(ph)),
                  f, indent=2, default=float)


# ---------------------------------------------------------------- figures
PALETTE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000", "#F0E442"]  # Okabe-Ito, CVD-validated order
GREY = "#B0B0B0"
LABEL_MAX = 12            # at most this many point labels per scatter (legibility)
INK2 = "#444444"


def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none", "font.size": 7,
        "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5, "axes.spines.top": False,
        "axes.spines.right": False, "legend.frameon": False})
    return plt


def _panel(ax, letter, x=-0.02):
    ax.text(x, 1.02, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", ha="right", va="bottom")


def _save(fig, out):
    fig.savefig(out + ".pdf", dpi=300)
    fig.savefig(out + ".png", dpi=300)


def fig_structure(out, ph, sgn):
    """A: chromosome link matrix (O/E) ordered by subgenome; B: support and own-vs-other preference."""
    plt = _plt()
    from matplotlib.colors import LinearSegmentedColormap
    C, N, z = ph["C"], ph["N"], ph["z"]
    chroms = ph["chroms"]
    order = sorted(range(C), key=lambda i: (z[i], natural_key(chroms[i])))
    W = ph["W"]
    L = np.log2((W + 1) / (expected(W) + 1))[np.ix_(order, order)]
    np.fill_diagonal(L, np.nan)
    fs = 6 if C <= 30 else (4 if C <= 60 else 2.5)
    fig = plt.figure(figsize=(7.1, 3.5))
    rend = fig.canvas.get_renderer()
    t_ = fig.text(0, 0, "", fontsize=fs)
    nm = 0.0                                             # widest chromosome name, inches
    for n_ in set(chroms):
        t_.set_text(n_)
        nm = max(nm, t_.get_window_extent(rend).width / fig.dpi)
    t_.remove()
    H = 3.5 + max(0.0, nm + 0.1 - 0.455)                 # room below the matrix for rotated names
    fig.set_size_inches(7.1, H)
    b_ = (H - 3.5 + 0.455) / H
    l_ = max(0.08, (nm + 0.1) / 7.1)
    ax = fig.add_axes([l_, b_, min(0.42 * 3.5 / 7.1 * 2.0, 0.5 - l_), 0.75 * 3.5 / H])
    lim = max(0.5, float(np.nanpercentile(np.abs(L), 98)))
    cmap = LinearSegmentedColormap.from_list("div", ["#2166AC", "#F2F2F2", "#B2182B"])
    im = ax.imshow(L, cmap=cmap, vmin=-lim, vmax=lim, interpolation="nearest")
    ax.set_xticks(range(C))
    ax.set_xticklabels([chroms[i] for i in order], rotation=90, fontsize=fs)
    ax.set_yticks(range(C))
    ax.set_yticklabels([chroms[i] for i in order], fontsize=fs)
    for sp_ in ax.spines.values():
        sp_.set_visible(False)
    ax.tick_params(length=0)
    for g in range(N):
        idx = [k for k, i in enumerate(order) if z[i] == g]
        if idx:
            ax.add_patch(plt.Rectangle((min(idx) - 0.5, -0.04 * C - 1.4), len(idx), 0.03 * C + 0.6, color=PALETTE[g], clip_on=False))
            ax.text((min(idx) + max(idx)) / 2, -0.05 * C - 1.6, f"{sgn[g]} (n={len(idx)})", ha="center", va="bottom",
                    fontsize=7, color="black")
    ax.set_xlim(-0.5, C - 0.5)
    ax.set_ylim(C - 0.5, -0.05 * C - 2.6)
    cax = fig.add_axes([ax.get_position().x1 + 0.01, b_, 0.008, 0.3 * 3.5 / H])
    cb = fig.colorbar(im, cax=cax)
    oe = (0.25, 0.5, 1, 2, 4) if lim >= 1 else (0.5, 0.7, 1, 1.4, 2)     # O/E ticks inside the colour range
    ticks = [np.log2(t) for t in oe if abs(np.log2(t)) <= lim]
    cb.set_ticks(ticks)
    cb.set_ticklabels([f"{2.0 ** t:.2g}" for t in ticks])
    cb.set_label("shared lineages,\nobserved / expected", fontsize=6)
    cb.ax.tick_params(labelsize=6)
    _panel(ax, "A")
    ax2 = fig.add_axes([0.66, b_ + 0.07 / H, 0.31, 0.75 * 3.5 / H])
    own, oth, sup = ph["own_oe"], ph["oth_oe"], ph["support"]
    lo = np.nanmin(np.concatenate([own, oth]))
    hi = np.nanmax(np.concatenate([own, oth]))
    pad = 0.05 * (hi - lo + 1e-9)
    ax2.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color=GREY, lw=0.6, ls="--", zorder=0)
    st = [x["status"] for x in ph.get("status", [{"status": "assigned"}] * C)]
    for i in range(C):
        if not np.isfinite(own[i]):
            continue
        col = PALETTE[z[i]]
        solid = st[i] == "assigned"
        ax2.scatter(oth[i], own[i], s=14, facecolor=col if solid else "white", edgecolor=col, lw=0.8, zorder=3)
    ax2.set_xlabel("observed/expected shared lineages\nwith the closest other subgenome")
    ax2.set_ylabel("observed/expected shared lineages with own subgenome")
    ax2.set_xlim(lo - pad, hi + pad)
    ax2.set_ylim(lo - pad, hi + pad)
    # label open points, weakest support first, wherever a label fits without touching another
    cand = sorted((i for i in range(C) if np.isfinite(own[i]) and st[i] != "assigned"),
                  key=lambda i: (st[i] != "conflict", np.nan_to_num(sup[i], nan=0.0)))
    from matplotlib.transforms import Bbox
    from matplotlib.text import Text
    fin = np.isfinite(own) & np.isfinite(oth)
    pts = ax2.transData.transform(np.column_stack([oth[fin], own[fin]]))
    rpx = 2.2 * fig.dpi / 72                                        # marker radius in pixels
    obst = [Bbox([[x - rpx, y - rpx], [x + rpx, y + rpx]]) for x, y in pts]
    boxes, nlab = [], 0
    for i in cand[:LABEL_MAX]:
        for dx, dy, ha in ((4, 0, "left"), (-4, 0, "right"), (8, 8, "left"), (8, -8, "left"), (-8, 8, "right"),
                           (-8, -8, "right"), (12, 16, "left"), (12, -16, "left"), (-12, 16, "right"),
                           (-12, -16, "right"), (16, 26, "left"), (16, -26, "left")):
            t = ax2.annotate(f"{chroms[i]} ({st[i]}, {sup[i]:.2f})", (oth[i], own[i]), xytext=(dx, dy),
                             textcoords="offset points", fontsize=5, color=INK2, ha=ha, va="center", zorder=5,
                             bbox=dict(fc="white", ec="none", pad=0.3), arrowprops=None if dy == 0 else dict(arrowstyle="-", lw=0.3, color=GREY, shrinkA=2.5,
                                                                  shrinkB=1.5))
            t.update_positions(rend)
            bb = Text.get_window_extent(t, rend).expanded(1.05, 1.15)      # the label alone, not its leader line
            inside = ax2.get_window_extent(rend).contains(bb.x0, bb.y0) and ax2.get_window_extent(rend).contains(bb.x1, bb.y1)
            if inside and not any(bb.overlaps(b) for b in boxes) and not any(bb.overlaps(b) for b in obst):
                boxes.append(bb)
                nlab += 1
                break
            t.remove()
    nopen = len(cand)
    ax2.text(0.97, 0.03, f"filled: assigned (support ≥ {SUPPORT_OK:g})\nopen: ambiguous, mosaic or conflict"
             + (f"\n({nlab} of {nopen} open points labelled, weakest first)" if nlab < nopen else ""),
             transform=ax2.transAxes, ha="right", va="bottom", fontsize=5.5, color=INK2)
    for g in range(N):
        ax2.text(0.02 + 0.16 * g, 1.02, f"● {sgn[g]}", transform=ax2.transAxes, color=PALETTE[g], fontsize=7,
                 fontweight="bold", va="bottom")
    _panel(ax2, "B")
    _save(fig, out)
    plt.close(fig)


def _segments_bar(ax, y, h, path_runs, colors, alpha=1.0):
    for x0, x1, g in path_runs:
        ax.add_patch(__import__("matplotlib").patches.Rectangle((x0, y), max(x1 - x0, 1e-4), h, color=colors[g], lw=0, alpha=alpha))


def fig_painting(out, lib, ph, cp, gm, sgn):
    """Chromosomes to scale; every LTR-RT drawn as a tick at its position, coloured by its ancestry
    call; ancestry bars (LTR-only HMM; genome painting); candidate exchanges boxed; plus a zoom."""
    plt = _plt()
    C, N, z = ph["C"], ph["N"], ph["z"]
    chroms = ph["chroms"]
    rows = sorted(range(C), key=lambda c: (z[c], natural_key(chroms[c])))
    clen = np.zeros(C)
    np.maximum.at(clen, ph["ek"][ph["keep"]], lib.end[ph["keep"]])
    for c in range(C):                                   # assembly lengths when known (--genome_fai / genome .fai)
        clen[c] = max(clen[c], ph.get("lengths", {}).get(chroms[c], 0))
    if gm:
        for c in range(C):
            if chroms[c] in gm["seq"]:
                clen[c] = gm["seq"][chroms[c]]["length"]
    flagged_ltr = [x for x in cp["segs"] if x["flag"]]
    flagged_g = [x for x in gm["segs"] if x["flag"]] if gm else []
    zoom = None
    pool = flagged_g or flagged_ltr
    zoom_src = "genome" if flagged_g else "ltr"
    if pool:                                    # a typical (median-evidence) candidate, not the strongest
        srt = sorted(pool, key=lambda x: x["llr"])
        zoom = srt[len(srt) // 2]
    R = len(rows)
    rh = 0.22 if R <= 30 else max(0.07, 6.5 / R)            # inches per chromosome row
    key_h, xl_h, gap, zh, zb = 0.55, 0.42, 0.3, 1.0 + (0.45 if gm else 0.0), 0.62
    H = key_h + rh * R + xl_h + ((gap + zh + zb) if zoom else 0.05)
    fig = plt.figure(figsize=(7.1, H))
    y0 = H - key_h - rh * R
    ax = fig.add_axes([0.12, y0 / H, 0.86, rh * R / H])
    fsz = 5.5 if R <= 30 else max(3.0, 5.5 * rh / 0.22)
    maxlen = clen.max() / 1e6
    for k, c in enumerate(rows):
        y = -k
        ax.plot([0, clen[c] / 1e6], [y, y], color=GREY, lw=0.5, zorder=1)
        idx = cp["order"][c]
        if idx.size:
            x = (lib.start[idx] + lib.end[idx]) / 2e6
            call = cp["call"][idx]
            lin = cp["lin"][idx]
            grey = ~np.isin(call, ["own", "foreign"])
            ax.vlines(x[grey], y + 0.04, y + 0.22, colors="#D0D0D0", lw=0.25, zorder=2, rasterized=True)
            for g in range(N):
                sel = (~grey) & (lin == g)
                if sel.any():
                    ax.vlines(x[sel], y + 0.04, y + 0.22, colors=PALETTE[g], lw=0.3, zorder=3, rasterized=True)
            st = cp["state"][idx]
            brk = np.flatnonzero(np.diff(st)) + 1
            runs = []
            for s0, e0 in zip(np.concatenate(([0], brk)), np.concatenate((brk, [idx.size]))):
                x0 = lib.start[idx[s0]] / 1e6 if s0 > 0 else 0.0
                x1 = lib.end[idx[e0 - 1]] / 1e6 if e0 < idx.size else clen[c] / 1e6
                runs.append((x0, x1, int(st[s0])))
            _segments_bar(ax, y - 0.13, 0.1, runs, PALETTE)
        if gm and chroms[c] in gm["seq"]:
            d = gm["seq"][chroms[c]]
            pth = d["path"]
            brk = np.flatnonzero(np.diff(pth)) + 1
            runs = [(s0 * WIN / 1e6, min(e0 * WIN, d["length"]) / 1e6, int(pth[s0]))
                    for s0, e0 in zip(np.concatenate(([0], brk)), np.concatenate((brk, [pth.size])))]
            _segments_bar(ax, y - 0.27, 0.1, runs, PALETTE)
        for x in flagged_ltr:
            if x["chrom"] == chroms[c]:
                ax.add_patch(plt.Rectangle((x["start"] / 1e6, y - 0.15), (x["end"] - x["start"]) / 1e6, 0.41, fill=False,
                                           edgecolor="black", lw=0.6, zorder=4))
        for x in flagged_g:
            if x["chrom"] == chroms[c]:
                ax.add_patch(plt.Rectangle((x["start"] / 1e6, y - 0.3), (x["end"] - x["start"]) / 1e6, 0.56, fill=False,
                                           edgecolor="black", lw=0.5, ls=(0, (2, 1.2)), zorder=4))
        sup = ph["support"][c]
        ax.text(-0.01 * maxlen, y, f"{chroms[c]}  {sup:.2f}" if np.isfinite(sup) else chroms[c], ha="right",
                va="center", fontsize=fsz, color=PALETTE[z[c]])
    ax.set_xlim(-0.005 * maxlen, maxlen * 1.01)
    ax.set_ylim(-len(rows) + 0.5, 0.5)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("position (Mb)")
    keytxt = ("ticks: every LTR-RT, coloured by called ancestry (grey: unresolved)   upper bar: ancestry along the "
              "chromosome from LTR-RTs" + ("   lower bar: genome painting" if gm else "") + "\nboxes: candidate exchanges"
              + (" (solid: LTR-RTs; dashed: genome)" if gm else "") + "   number after name: support")
    for g in range(N):
        fig.text(0.12 + 0.07 * g, (H - 0.05) / H, f"\u25a0 {sgn[g]}", color=PALETTE[g], fontsize=6.5, va="top",
                 ha="left", fontweight="bold")
    fig.text(0.12, (H - 0.2) / H, keytxt, fontsize=5.5, color=INK2, va="top", ha="left")
    _panel(ax, "A")
    if zoom:
        axz = fig.add_axes([0.12, zb / H, 0.86, zh / H])
        c = chroms.index(zoom["chrom"])
        span = zoom["end"] - zoom["start"]
        x0, x1 = max(0, zoom["start"] - span), min(zoom["end"] + span, clen[c])
        yr_ = -rows.index(c)                         # zoomed window outlined in A and joined to B
        ax.add_patch(plt.Rectangle((x0 / 1e6, yr_ - 0.36), (x1 - x0) / 1e6, 0.64, fill=False, edgecolor=GREY, lw=0.6,
                                   zorder=5))
        from matplotlib.patches import ConnectionPatch
        idx = cp["order"][c]
        idx = idx[(lib.end[idx] >= x0) & (lib.start[idx] <= x1)]
        # every element as its own box, stacked into rows when they overlap (IGV expanded)
        lanes = []
        for i in idx[np.argsort(lib.start[idx])]:
            for li, last in enumerate(lanes):
                if lib.start[i] > last + span * 0.002:
                    lanes[li] = lib.end[i]
                    lane = li
                    break
            else:
                lanes.append(lib.end[i])
                lane = len(lanes) - 1
            g = cp["lin"][i]
            col = PALETTE[g] if g >= 0 else "#D0D0D0"
            axz.add_patch(plt.Rectangle((lib.start[i] / 1e6, 0.1 + lane * 0.25), (lib.end[i] - lib.start[i]) / 1e6, 0.18,
                                        color=col, lw=0))
        lab_tr = axz.get_yaxis_transform()
        axz.text(-0.01, 0.1 + 0.125 * max(len(lanes), 1), "LTR-RTs", transform=lab_tr, ha="right", va="center", fontsize=5.5,
                 color=INK2)
        tracks = []                                  # posterior ancestry tracks, each 0..1 tall
        if cp.get("posts") is not None and c in cp["posts"]:
            o_ = cp["order"][c]
            m_ = (lib.end[o_] >= x0) & (lib.start[o_] <= x1)
            xm_ = (lib.start[o_][m_] + lib.end[o_][m_]) / 2e6
            tracks.append(("LTR-RT HMM\nancestry", xm_, cp["posts"][c][m_]))
        if gm and zoom["chrom"] in gm["seq"]:
            d = gm["seq"][zoom["chrom"]]
            w0, w1 = int(x0 // WIN), int(x1 // WIN) + 1
            xs = (np.arange(w0, min(w1, d["post"].shape[0])) + 0.5) * WIN / 1e6
            tracks.append(("genome window\nancestry", xs, d["post"][w0:w0 + xs.size]))
        for ti, (tl_, xs, Pt) in enumerate(tracks):
            base = -1.2 - 1.25 * ti
            for g in range(N):
                axz.plot(xs, base + Pt[:, g], color=PALETTE[g], lw=0.9, drawstyle="steps-mid" if ti == 0 else "default")
            for yv, tl in ((base, "0"), (base + 1, "1")):
                axz.plot([0, 0.008], [yv, yv], transform=lab_tr, color=INK2, lw=0.6, clip_on=False)
                axz.text(-0.003, yv, tl, transform=lab_tr, ha="right", va="center", fontsize=5, color=INK2)
            axz.text(-0.02, base + 0.5, tl_ + "\n(posterior)", transform=lab_tr, ha="right", va="center",
                     fontsize=5.5, color=INK2)
        ybot = -1.3 - 1.25 * max(len(tracks) - 1, 0)
        axz.add_patch(plt.Rectangle((zoom["start"] / 1e6, ybot), span / 1e6, -ybot + 0.1 + 0.25 * max(len(lanes), 1),
                                    fill=False, edgecolor="black", lw=0.6, ls=(0, (2, 1.2)) if zoom_src == "genome" else "-"))
        axz.set_xlim(x0 / 1e6, x1 / 1e6)
        axz.set_ylim(ybot - 0.05, 0.35 + 0.25 * max(len(lanes), 1))
        from matplotlib.transforms import blended_transform_factory
        yb_ = (ax.get_tightbbox(fig.canvas.get_renderer()).y0 - 2) / fig.bbox.height   # below A's axis label
        xf_ = blended_transform_factory(ax.transData, fig.transFigure)
        for xa in (x0, x1):                          # zoom lines from the outlined window in A to B
            fig.add_artist(ConnectionPatch((xa / 1e6, yb_), (xa / 1e6, axz.get_ylim()[1]), coordsA=xf_,
                                           coordsB="data", axesB=axz, color="#CFCFCF", lw=0.4, zorder=0))
        axz.set_yticks([])
        axz.spines["left"].set_visible(False)
        axz.set_xlabel(f"{zoom['chrom']} position (Mb)\n{'dashed' if zoom_src == 'genome' else 'solid'} box: candidate "
                       f"exchange of {sgn[zoom['origin']]} ancestry from {'the genome' if zoom_src == 'genome' else 'LTR-RTs'} "
                       f"(median evidence); small boxes: LTR-RTs coloured by called ancestry")
        _panel(axz, "B")
    _save(fig, out)
    plt.close(fig)
    return zoom


def fig_ages(out, a, ph, cp, ck, sgn):
    """One divergence axis (linear below K2P 0.002, log above), two stacked panels. A: when each subgenome's
    own lineages were inserted (copies of own-subgenome ancestry, per subgenome; bars: SubPhaser-style 95%
    window). B: cross-subgenome transposition clock - the share of other-subgenome ancestry per age, the
    fitted step and the merger (a date only when its interval excludes zero, else a 'younger than' band)."""
    plt = _plt()
    N = ph["N"]
    k2p = ck["k2p"]
    own = ph["z"][ph["ek"]]
    ok = np.isfinite(k2p)
    hi = max(float(np.nanquantile(k2p[ok], 0.995)), 1e-2)
    lin = 0.002
    edges = np.concatenate([np.linspace(0, lin, 5), np.geomspace(lin, hi, 33)[1:]])
    fit = ck["fit"]
    pm = fit.get("pooled") if fit else None
    yr = lambda x: x / (2 * a.mu)
    fig, axs = plt.subplots(2, 1, figsize=(7.1, 4.6), sharex=True, gridspec_kw=dict(height_ratios=[1, 1.15]))
    fig.subplots_adjust(left=0.1, right=0.97, bottom=0.11, top=0.88, hspace=0.12)
    for ax in axs:
        ax.set_xscale("symlog", linthresh=lin, linscale=0.6)
        ax.set_xlim(0, hi)
    # merger band drawn behind both panels
    if pm:
        if pm["tau"] > 0:
            for ax in axs:
                ax.axvspan(pm["ci"][0], pm["ci"][1], color="black", alpha=0.08, lw=0, zorder=0)
                ax.axvline(pm["tau"], color="black", lw=0.9, zorder=1)
        else:
            for ax in axs:
                ax.axvspan(0, pm["upper"], color="black", alpha=0.08, lw=0, zorder=0, hatch="////", ec=GREY)
    if fit and fit.get("split_consistent"):
        td, (dl, dh) = fit["tau_divergence"], fit["tau_divergence_ci"]
        for ax in axs:
            ax.axvspan(dl, dh, color=GREY, alpha=0.12, lw=0, zorder=0)
            ax.axvline(td, color="black", lw=0.8, ls=(0, (3, 2)), zorder=1)
        fx = axs[1].transAxes.inverted().transform(axs[1].transData.transform((td, 0)))[0]   # keep the label inside
        axs[1].text(td, 0.83, f" progenitor split ≥ {yr(td) / 1e6:,.1f} Myr \n (tentative; dashed) ", fontsize=6,
                    ha="left" if fx < 0.75 else "right", va="top", transform=axs[1].get_xaxis_transform(), zorder=6,
                    bbox=dict(fc="white", ec="none", pad=0.6))
    ax = axs[0]
    ymax = 0
    labs = []
    for g in range(N):
        v = k2p[(cp["call"] == "own") & (own == g) & ok]
        if v.size < 5:
            continue
        hc, _ = np.histogram(np.clip(v, 0, hi), edges)
        f = hc / v.size
        ax.step(edges[:-1], f, where="post", color=PALETTE[g], lw=1.0)
        ymax = max(ymax, f.max())
        labs.append((g, v.size, edges[int(np.argmax(f))]))
    vf = k2p[(cp["call"] == "foreign") & ok]                # copies carrying another subgenome's ancestry
    if vf.size >= 5:
        hc, _ = np.histogram(np.clip(vf, 0, hi), edges)
        ax.step(edges[:-1], hc / vf.size, where="post", color=INK2, lw=0.9, ls=(0, (3, 1.5)))
        ymax = max(ymax, (hc / vf.size).max())
    for i, (g, n, xpk) in enumerate(labs):          # 95% windows as bars above the curves, labelled directly
        yb = ymax * (1.12 + 0.18 * i)
        if g in ck["window"]:
            w = ck["window"][g]
            ax.plot([w[0], w[2]], [yb, yb], color=PALETTE[g], lw=1.6, solid_capstyle="butt")
            ax.plot([w[1]], [yb], "|", color=PALETTE[g], ms=5, mew=1.2)
        ax.text(hi, yb + 0.04 * ymax, f"  {sgn[g]}: own-ancestry copies in {sgn[g]} (n={n:,})", color=PALETTE[g],
                fontsize=6, ha="right", va="bottom", zorder=6, bbox=dict(fc="white", ec="none", pad=0.6))
    if vf.size >= 5:
        yb = ymax * (1.12 + 0.18 * len(labs))
        ax.text(hi, yb, f"  dashed: foreign-ancestry copies, all subgenomes (n={vf.size:,})", color=INK2, fontsize=6,
                ha="right", va="bottom", zorder=6, bbox=dict(fc="white", ec="none", pad=0.6))
    ax.set_ylim(0, ymax * (1.3 + 0.18 * (len(labs) + (vf.size >= 5))))
    ax.set_ylabel("fraction of copies\nper divergence bin")
    top = ax.twiny()                                  # insertion-time ticks on the same (symlog) scale
    top.set_xscale("symlog", linthresh=lin, linscale=0.6)
    top.set_xlim(ax.get_xlim())
    tt = [(t, lab_) for t, lab_ in ((1e3, "1 kyr"), (1e4, "10 kyr"), (1e5, "100 kyr"), (1e6, "1 Myr"), (1e7, "10 Myr"),
                                    (1e8, "100 Myr")) if lin / 2 <= 2 * a.mu * t <= hi]
    top.set_xticks([0] + [2 * a.mu * t for t, _ in tt])
    top.set_xticklabels(["0"] + [l_ for _, l_ in tt], fontsize=6)
    top.minorticks_off()
    top.set_xlabel(f"insertion time (μ = {a.mu:g} per site per year" + ("" if getattr(a, "mu_given", True) else
                   ", rice default") + ")", fontsize=6)
    _panel(ax, "A")
    ax = axs[1]
    xc = np.where(edges[:-1] < lin, (edges[:-1] + edges[1:]) / 2, np.sqrt(np.maximum(edges[:-1], 1e-12) * edges[1:]))

    def binned(xk, val, wt):
        bi = np.clip(np.searchsorted(edges, xk, side="right") - 1, 0, edges.size - 2)
        sw_ = np.bincount(bi, wt, edges.size - 1)
        return sw_, np.bincount(bi, wt * val, edges.size - 1) / np.where(sw_ > 0, sw_, np.nan)
    if fit and ck.get("share") is not None:
        use = np.isfinite(ck["k2p"]) & np.isfinite(ck["share"]) & (ck["w"] > 0)
        Wt, Y = binned(np.clip(ck["k2p"][use], 0, hi), ck["share"][use], ck["w"][use])
        keep_ = Wt > 0
        ax.scatter(xc[keep_], Y[keep_], s=4 + 40 * Wt[keep_] / Wt.max(), color=INK2, zorder=3, lw=0)
        sn = fit.get("s_null", np.nan)
        if np.isfinite(sn):
            ax.axhline(sn, color=GREY, ls="--", lw=0.6)
            ax.text(hi, sn, "no subgenome confinement ", fontsize=5.5, color=INK2, ha="right", va="bottom", zorder=6,
                    bbox=dict(fc="white", ec="none", pad=0.6))
        if pm:
            tau = pm["tau"]
            sw_, mp = binned(pm["pred_k2p"], pm["pred"], pm["pred_w"])
            ax.plot(xc[sw_ > 0], mp[sw_ > 0], color="black", lw=0.9, zorder=4, marker="_", ms=4)
            if tau > 0:
                lab = (f"merger {yr(tau) / 1e3:,.1f} kyr (95% CI {yr(pm['ci'][0]) / 1e3:,.1f}–{yr(pm['ci'][1]) / 1e3:,.1f}); "
                       f"black line: fitted step model, averaged like the points")
            else:
                lab = (f"step not resolved: cross-subgenome transposition began < {yr(pm['upper']) / 1e3:,.1f} kyr ago "
                       f"(hatched; 95% bound); black line: fitted model, averaged like the points")
            ax.text(0.01, 0.97, lab, transform=ax.transAxes, fontsize=6, ha="left", va="top", color="black", zorder=6,
                    bbox=dict(facecolor="white", edgecolor="none", pad=1.0, alpha=0.9))
        ax.set_ylim(0, max(0.75, float(np.nanmax(Y[keep_])) * 1.15))
    else:
        ax.text(0.5, 0.5, "too few dated LTR-RTs for the clock", transform=ax.transAxes, ha="center", va="center",
                fontsize=6.5, color=INK2)
    ax.set_xlabel("LTR divergence (K2P; linear below 0.002, log above)")
    ax.set_ylabel("mean probability that a copy\nhas another subgenome's ancestry")
    _panel(ax, "B")
    _save(fig, out)
    plt.close(fig)


def make_figures(a, lib, ph, cp, ck, gm):
    sgn = [f"SG{g + 1}" for g in range(ph["N"])]
    od = a.outdir
    ph["lengths"] = getattr(a, "lengths", {})
    fig_structure(os.path.join(od, "fig_structure"), ph, sgn)
    zoom = fig_painting(os.path.join(od, "fig_painting"), lib, ph, cp, gm, sgn)
    if ck and np.isfinite(ck["k2p"]).any():
        fig_ages(os.path.join(od, "fig_ages"), a, ph, cp, ck, sgn)
    write_legends(a, lib, ph, cp, ck, gm, sgn, zoom)


def write_legends(a, lib, ph, cp, ck, gm, sgn, zoom):
    od = a.outdir
    C, N = ph["C"], ph["N"]
    st = [x["status"] for x in ph["status"]]
    nass = sum(x == "assigned" for x in st)
    calls = {v: int(c) for v, c in zip(*np.unique(cp["call"][ph["keep"]], return_counts=True))}   # drawn copies only
    shared = ("Terms. *LTR-RT*: an intact long terminal repeat retrotransposon copy. *Shared lineage (link)*: two copies on "
              f"different chromosomes that carry the same rare {a.k}-mer (carried by at most {ph['window']} copies: a "
              "clade-defining sequence variant), or that are each "
              "other's closest relatives; copies on the same chromosome never count. *Support*: fraction of 100 random "
              "half-libraries (50% of copies; each copy's closest relatives taken from the full library and restricted to the "
              "half) that place a chromosome in the same "
              "subgenome. *Ancestry* of a copy or genome window: which subgenome's lineages its k-mers belong to, judged "
              "only from how those k-mers are spread over the other chromosomes and calibrated against the chromosomes' "
              "own assignments; k-mers that could be homoeologous copies are ignored.")
    with open(os.path.join(od, "fig_structure_legend.md"), "w") as f:
        f.write(f"# fig_structure\n\n{shared}\n\n"
                f"**A.** Every pair of the {C} chromosomes with at least {ph.get('min_elements', MIN_ELEMENTS)} LTR-RTs, coloured by shared lineages "
                f"observed/expected under a null in which each chromosome keeps its total (red: more than expected; blue: "
                f"fewer). Chromosomes are ordered by inferred subgenome (coloured bars). Each progenitor amplified its own "
                f"lineages before the merger, so chromosomes of one subgenome share more lineages "
                f"(mean observed/expected {np.nanmean(ph['own_oe']):.2f}) than chromosomes of different subgenomes "
                f"({np.nanmean(ph['oth_oe']):.2f})."
                + (f" {len(_partner_pairs_json(ph))} chromosome pairs are each other's closest partner within their "
                   f"subgenome in ≥ {SUPPORT_OK:.0%} of half-libraries (bright cells off the diagonal; median "
                   f"{np.median([x['excess'] for x in _partner_pairs_json(ph)]):.1f}× the subgenome's typical sharing); "
                   f"strong pairs are near-identical chromosomes (homologues of an autopolyploid component, or of two "
                   f"progenitors too close for LTR-RTs to separate) and count as one chromosome in all ancestry scores."
                   if _partner_pairs_json(ph) else "")
                + "\n\n"
                f"**B.** Each chromosome's mean observed/expected shared lineages with its own subgenome against the "
                f"closest other subgenome; points above the dashed diagonal prefer their own subgenome. Filled: assigned "
                f"(support ≥ {SUPPORT_OK:g} and no conflict with ancestry painting; {nass}/{C}); open: ambiguous, mosaic or "
                f"conflicting, labelled with status and support (weakest first, up to {LABEL_MAX}).\n")
    with open(os.path.join(od, "fig_painting_legend.md"), "w") as f:
        nf = sum(x["flag"] for x in cp["segs"])
        ng = sum(x["flag"] for x in gm["segs"]) if gm else 0
        f.write(f"# fig_painting\n\n{shared}\n\n"
                f"**A.** Chromosomes drawn to scale (Mb), grouped by subgenome. Ticks above each chromosome: every LTR-RT at "
                f"its position, coloured by its called ancestry (posterior ≥ {P_CALL:g}; {calls.get('own', 0):,} copies of their "
                f"chromosome's subgenome, {calls.get('foreign', 0):,} of another subgenome's ancestry) or grey when unresolved "
                f"({calls.get('unresolved', 0) + calls.get('no_relatives', 0):,}). Bar below: ancestry along the chromosome "
                f"from a hidden Markov model over the copies"
                + (" ; lower bar: ancestry painted from genome-wide k-mers in 50-kb windows" if gm else "")
                + f". Solid boxes: candidate exchanges from LTR-RTs alone (likelihood ratio ≥ {BF_DECISIVE:g} over the "
                f"chromosome's subgenome and copies from ≥ 2 independent lineages; n={nf})"
                + (f"; dashed boxes: candidate exchanges from the genome (likelihood ratio ≥ {BF_DECISIVE:g} and mean posterior ≥ 0.99, "
                   f"≥ {MIN_SEG_WIN * WIN // 1000} kb; n={ng}; weaker segments are listed in "
                   f"segments_genome.tsv but not drawn)" if gm else "")
                + ". Candidates are homoeologous exchanges or assembly switch errors; the data cannot tell them apart. "
                f"Numbers after chromosome names: support.\n\n"
                + (f"**B.** Zoom on one candidate exchange ({zoom['chrom']}:{zoom['start']:,}-{zoom['end']:,}, "
                   f"{sgn[zoom['origin']]} ancestry on a {sgn[zoom['chrom_sg']]} chromosome), chosen as the candidate with the "
                   f"median evidence, not the strongest. Each LTR-RT is drawn as its own box (overlapping copies stacked), "
                   f"coloured by called ancestry"
                   + ("; lines: posterior ancestry of each 50-kb genome window" if gm else "") + ".\n" if zoom else ""))
    if ck and np.isfinite(ck["k2p"]).any():
        fit = ck["fit"]
        pm = fit.get("pooled") if fit else None
        ky = lambda x: f"{x / (2 * a.mu) / 1e3:,.1f}"
        t = ("too few dated copies to fit the clock" if not pm else
             (f"cross-subgenome transposition began {ky(pm['tau'])} kyr ago (95% bootstrap CI {ky(pm['ci'][0])}-"
              f"{ky(pm['ci'][1])} kyr; {pm['n']:,} young copies), a minimum age for the merger" if pm["tau"] > 0 else
              f"cross-subgenome transposition began less than {ky(pm['upper'])} kyr ago (one-sided 95% bound; below the "
              f"resolution of LTR divergence); the merger itself can be older if transposition was quiet after it"))
        win = "; ".join(f"{sgn[g]}: {v[0] / (2 * a.mu) / 1e6:.2f}-{v[2] / (2 * a.mu) / 1e6:.2f} Myr (n={v[3]:,})"
                        for g, v in ck["window"].items())
        with open(os.path.join(od, "fig_ages_legend.md"), "w") as f:
            f.write(f"# fig_ages\n\n{shared} Ages assume {a.mu:g} substitutions per site per year (insertion time = K2P / 2μ"
                    + ("" if getattr(a, "mu_given", True) else "; the default rice LTR rate - set --mu for other taxa")
                    + "); lineage rates differ, so ages are relative. LTR gene conversion, loss of old copies and a lag "
                    "before the first cross-subgenome transposition all make the clock young: it gives minimum ages.\n\n"
                    f"**A.** When each subgenome's own lineages were inserted: LTR divergence of copies with their own subgenome's "
                    f"ancestry, per subgenome (fraction of copies per bin; one axis for both panels, linear below K2P 0.002 and "
                    f"log above, so young bursts stay visible). Bars above the curves: central 95% of these ages with the "
                    f"median ticked (the SubPhaser divergence-hybridization window; {win}). Dashed grey: copies carrying "
                    f"another subgenome's ancestry (foreign calls), expected mostly after the merger or inside exchanges. "
                    f"The window is rough: lineages may "
                    f"burst at any time, old copies decay and lineages confined to a subgenome can keep transposing after "
                    f"the merger.\n\n"
                    f"**B.** Cross-subgenome transposition clock. Points: copies binned by divergence (area = number of copies); y: "
                    f"mean posterior that a copy belongs to another subgenome's lineages. Before the merger a lineage could "
                    f"only spread within its own progenitor, so copies of that age rarely carry another subgenome's ancestry; "
                    f"after it, copies can land anywhere (dashed: expectation without confinement). Black line: the step "
                    f"fitted to every young copy at the resolution of its own LTR length (no confinement before the step, "
                    f"a confined share that may drift with age after it); grey band: 95% bootstrap interval, hatched when the "
                    f"step is not resolved (interval reaches zero) and only an upper bound is given. "
                    + (f"Dashed line and light band (both panels): beyond K2P {fit['tau_divergence']:.3g} "
                       f"({fit['tau_divergence'] / (2 * a.mu) / 1e6:,.1f} Myr; 95% CI "
                       f"{fit['tau_divergence_ci'][0] / (2 * a.mu) / 1e6:,.1f}-{fit['tau_divergence_ci'][1] / (2 * a.mu) / 1e6:,.1f}) "
                       f"lineages are shared between the subgenomes again - a tentative minimum age of progenitor divergence "
                       f"(the split itself, or the age limit of surviving LTR-RTs when the split is older: old copies are lost "
                       f"or saturated); shown only when the rise is supported and older than the median own-ancestry copy of "
                       f"every subgenome. " if fit and fit.get("split_consistent") else "")
                    + f"Here {t}.\n")




if __name__ == "__main__":
    main()
