#!/usr/bin/env python3
"""
config-free subgenome phasing of an allopolyploid from its LTR-RT library.

Idea
----
Between progenitor divergence and allopolyploidization, each progenitor grew
its own LTR-RT lineages. Copies of such a lineage share derived mutations, so
they share k-mers that no other copy carries. v4 treats every LTR-RT as one
unit of evidence and every k-mer carried by a small set of copies (a "split",
i.e. a clade of the family tree) as a link between the chromosomes those
copies sit on. Copies on the same chromosome are ignored, so local tandem
duplication, centromeric arrays and other chromosome-specific amplification
cannot create a signal. The chromosome x chromosome link matrix is compared
with a degree-corrected null (observed / expected), and chromosomes are split
into N subgenomes by maximizing within-subgenome excess links (modularity of a
degree-corrected block model; exact for N=2 with <=23 free units). Support is
the fraction of random half-libraries (50% of LTR-RTs) giving the same
assignment; in benchmarks it is calibrated (support >= 0.9 -> >= 97% correct).

No homoeolog config, number of subgenomes, genome or jellyfish is needed.
The number of subgenomes is inferred by divisive splitting: a group of
chromosomes is split only if two disjoint halves of the library reproduce the
split (adjusted Rand index >= 0.8) and, below the top level, the split is at
least half as strong as its parent split. N=1 ("no subgenome structure", e.g.
a diploid or autopolyploid) is a possible answer; uneven numbers of
chromosomes per subgenome and missing homoeologs are fine. -n forces N. A
config, if given, only adds "these chromosomes are in different subgenomes"
constraints.

The copy-number window of informative k-mers (rare = specific, common =
numerous) is chosen per genome as the one whose chromosome link pattern is
most reproducible between disjoint halves of the library; each distinct
carrier set (clade) counts once.

Downstream:
  * each LTR-RT is classified by where its closest relatives on other
    chromosomes sit: own-subgenome lineage, other-subgenome ('foreign')
    lineage, or shared lineage (empirical-Bayes mixture; a call also needs a
    Bayes factor >= 3 from the copy's own data);
  * with LTR K2P divergences (--k2p), the cross-subgenome transposition clock
    dates the merger: copies inserted after it have relatives in any
    subgenome, copies inserted while the progenitors were apart do not. All
    young copies are fitted jointly at the resolution of their pooled LTR
    sites (nonparametric age distribution from substitution counts), so
    mergers younger than one substitution per LTR can be dated or bounded;
  * an HMM along each chromosome flags runs of copies whose relatives point to
    another subgenome (candidate homoeologous exchanges or switch errors;
    likelihood ratio >= 100); correlated local copies and copies of age
    classes without subgenome confinement are down-weighted;
  * with --genome, every 50-kb window is painted with k-mers of lineages
    confined to one subgenome (fragmented and solo LTRs included): a second,
    denser exchange track, and a subgenome for every sequence (e.g. unplaced
    scaffolds without intact LTR-RTs).

Remaining heuristics (explicit options): auto-N accepts a split at split-half
replicability >= 0.8 (--min_rep) and, below the top, only if it is >= 0.5x as
strong as its parent split (subgenomes are assumed to be the dominant
structure).

Input FASTA headers: chrom:start-end#Class/Superfamily/Family (LTRquest,
LTR_retriever, EDTA-style); '#...' is optional.

Example
-------
  python v4.py --ltr_fasta ltr.fa --outdir out --k2p ltr.tsv      # N inferred, no config
  python v4.py --ltr_fasta ltr.fa --outdir out -n 3                # force N=3
  python v4.py --ltr_fasta ltr.fa --outdir out --k2p ltr.tsv --genome asm.fa   # + genome painting

Outputs: chromosomes.tsv, elements.tsv, segments.tsv, links_oe.tsv,
summary.json, fig_phasing/fig_painting/fig_ages (.pdf/.png), FIGURE_LEGEND.md;
with --genome also segments_genome.tsv, genome_sequences.tsv, genome_windows.tsv.
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
import os
import sys
import time

import numpy as np
from scipy import sparse
from scipy.optimize import linear_sum_assignment

__version__ = "4.2"

# ---------------------------------------------------------------- logging
_T0 = time.time()
VERBOSE = False


def log(msg):
    print(f"[v4 {time.time() - _T0:7.1f}s] {msg}", file=sys.stderr, flush=True)


def vlog(msg):
    if VERBOSE:
        log(msg)


def die(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


# ---------------------------------------------------------------- input
def read_fasta(path):
    """Yield (id, sequence bytes); id is the first header token."""
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
    """chrom:start-end[#cls/sf/fam] or chrom:start..end[#...] (LTR_retriever)
    -> (chrom, start, end, cls) or None."""
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


def read_k2p(path, k2p_col=None):
    """LTR id -> (K2P, K2P s.e., aligned LTR sites, substitutions). Accepts an
    LTRquest/Kmer2LTR table (header with seq_id and k2p; k2p_se, n_sites (else
    ltr5_len), n_ts + n_tv used if present) or any whitespace table with the id
    in column 1 and K2P in --k2p_col (1-based)."""
    out, bad = {}, 0
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        first = f.readline()
        hdr = first.lstrip("#").rstrip("\n").split("\t")
        ix = {}
        if k2p_col is None and "k2p" in hdr:
            ix = {c: hdr.index(c) for c in ("k2p_se", "n_sites", "ltr5_len", "n_ts", "n_tv") if c in hdr}
            ci = hdr.index("seq_id") if "seq_id" in hdr else 0
            ki = hdr.index("k2p")
            lines = f
            sep = "\t"
        else:
            ci, ki = 0, (k2p_col or 11) - 1
            lines = [first] + list(f)
            sep = None
        nan = float("nan")
        for line in lines:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split(sep)
            if len(p) <= max(ci, ki):
                bad += 1
                continue
            try:
                k2p = float(p[ki])
            except ValueError:
                bad += 1
                continue

            def num(c):
                try:
                    return float(p[ix[c]]) if c in ix else nan
                except (ValueError, IndexError):
                    return nan
            L = num("n_sites")
            if not L == L:
                L = num("ltr5_len")
            ksub = num("n_ts") + num("n_tv")
            out[p[ci].strip().lstrip(">")] = (k2p, num("k2p_se"), L, ksub)
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


def _chunk_keys(codes, eids, k, scale):
    """Canonical k-mers of concatenated, N-separated elements -> unique
    (hash40 << 24 | element) keys, FracMinHash-subsampled 1/scale."""
    n = codes.size - k + 1
    if n <= 0:
        return np.empty(0, np.uint64)
    bad = codes == 4
    cs = np.zeros(codes.size + 1, np.int32)
    np.cumsum(bad, out=cs[1:])
    ok = (cs[k:k + n] - cs[:n]) == 0
    c = codes.astype(np.uint64)
    c[bad] = 0
    f = np.zeros(n, np.uint64)
    r = np.zeros(n, np.uint64)
    two, three = np.uint64(2), np.uint64(3)
    for j in range(k):
        cj = c[j:j + n]
        f <<= two
        f |= cj
        r |= (three - cj) << np.uint64(2 * j)
    h = mix64(np.minimum(f, r)[ok])
    e = eids[:n][ok]
    if scale > 1:
        m = (h % np.uint64(scale)) == 0
        h, e = h[m], e[m]
    return np.unique(((h >> _SH) << _SH) | e.astype(np.uint64))


class Library:
    """Parsed LTR-RT library: element table + sorted k-mer/element keys."""

    def __init__(self, path, k, scale, keep_chroms=None, chunk=1 << 22, min_len=None):
        self.ids, self.chrom, self.start, self.end, self.cls, self.length = [], [], [], [], [], []
        min_len = min_len or k
        bad_hdr = short = 0
        keys, buf, ebuf, nb = [], [], [], 0
        for h, s in read_fasta(path):
            p = parse_id(h)
            if p is None:
                bad_hdr += 1
                continue
            if keep_chroms is not None and p[0] not in keep_chroms:
                continue
            if len(s) < min_len:
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
            buf.append(codes)
            buf.append(np.full(1, 4, np.uint8))
            ebuf.append(np.full(codes.size + 1, i, np.uint32))
            nb += codes.size + 1
            if nb >= chunk:
                keys.append(_chunk_keys(np.concatenate(buf), np.concatenate(ebuf), k, scale))
                buf, ebuf, nb = [], [], 0
        if buf:
            keys.append(_chunk_keys(np.concatenate(buf), np.concatenate(ebuf), k, scale))
        if bad_hdr:
            log(f"WARNING: skipped {bad_hdr} records whose header is not chrom:start-end[#class]")
        if short:
            vlog(f"skipped {short} records shorter than {min_len} bp")
        if not self.ids:
            die(f"no usable LTR-RT records in {path}")
        keys = np.concatenate(keys) if keys else np.empty(0, np.uint64)
        keys.sort()
        self.kmer = keys >> _SH
        self.eid = (keys & _EMASK).astype(np.int32)
        self.n = len(self.ids)
        self.start = np.array(self.start, np.int64)
        self.end = np.array(self.end, np.int64)
        self.length = np.array(self.length, np.int64)
        self.chroms = sorted(set(self.chrom), key=natural_key)


def natural_key(s):
    import re
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", s)]


# ---------------------------------------------------------------- splits
class Splits:
    """Distinct carrier sets of shared k-mers that span >= 2 chromosomes.
    CSR: members[indptr[s]:indptr[s+1]] are the elements carrying split s;
    mult[s] = number of sampled k-mers with exactly that carrier set."""

    def __init__(self, kmer, eid, echrom, n_chrom, nmax, emask=None):
        if emask is not None:
            keep = emask[eid]
            kmer, eid = kmer[keep], eid[keep]
        if kmer.size == 0:
            die("no k-mers left")
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
        self.n_shared_kmers = int(sizes[sizes >= 2].size)
        self.sid = np.repeat(np.arange(self.n), self.size)


def split_view(sp, mask):
    """Subset of a Splits object (e.g. carrier sets with <= n copies)."""
    v = Splits.__new__(Splits)
    kp = np.repeat(mask, sp.size)
    v.members, v.size, v.mult, v.nch = sp.members[kp], sp.size[mask], sp.mult[mask], sp.nch[mask]
    v.n = int(mask.sum())
    v.indptr = np.concatenate(([0], np.cumsum(v.size)))
    v.sid = np.repeat(np.arange(v.n), v.size)
    v.n_shared_kmers = sp.n_shared_kmers
    return v


def choose_window(sp, echrom, C, n_el, keep, cands, rng, reps=4):
    """Copy-number window for phasing chosen from the data: the largest carrier
    count whose chromosome link pattern (log observed/expected) is most
    reproducible between disjoint halves of the library (mean Pearson r).
    Rare k-mers mark the most specific lineages, common ones add counts;
    the best balance differs between genomes."""
    iu = np.triu_indices(C, 1)
    scores = {}
    halves = [rng.random(n_el) < 0.5 for _ in range(reps)]
    for nm in cands:
        v = split_view(sp, sp.size <= nm)
        if v.n == 0:
            continue
        r = []
        for u in halves:
            Ls = []
            for part in (u & keep, ~u & keep):
                W = affinity(v, echrom, C, part.astype(float), "chrom", 0.0)
                Ls.append(np.log((W + 1) / (expected(W) + 1))[iu])
            if Ls[0].std() > 0 and Ls[1].std() > 0:
                r.append(np.corrcoef(Ls[0], Ls[1])[0, 1])
        scores[nm] = float(np.mean(r)) if r else -1.0
    best = max(scores, key=lambda x: (round(scores[x], 3), x))
    return best, scores


# ---------------------------------------------------------------- affinity
def chrom_counts(sp, echrom, C, eweight=None):
    """S x C sparse carrier counts per split and chromosome."""
    data = np.ones(sp.members.size) if eweight is None else eweight[sp.members].astype(float)
    return sparse.csr_matrix((data, (sp.sid, echrom[sp.members])), shape=(sp.n, C))


def affinity(sp, echrom, C, eweight=None, wmode="chrom", mult_pow=0.0):
    """C x C link matrix W: each split adds w to every pair of chromosomes it
    spans; w = 1/(m-1) for m spanned chromosomes ('chrom'), 1/(n-1) for n
    carriers ('elem') or 1 ('one'); times mult**mult_pow. Diagonal is zero."""
    X = chrom_counts(sp, echrom, C, eweight)
    P = X.copy()
    P.data = (P.data > 0).astype(float)
    P.eliminate_zeros()
    m = np.asarray(P.sum(1)).ravel()
    if wmode == "chrom":
        w = 1.0 / np.maximum(m - 1, 1)
    elif wmode == "elem":
        n = np.asarray(X.sum(1)).ravel()
        w = 1.0 / np.maximum(n - 1, 1)
    else:
        w = np.ones(sp.n)
    if mult_pow:
        w = w * sp.mult ** mult_pow
    w[m < 2] = 0
    W = (P.T @ sparse.diags(w) @ P).toarray()
    np.fill_diagonal(W, 0)
    return W


def expected(W, iters=200):
    """Degree-corrected null for a zero-diagonal matrix: E_ij = b_i b_j with
    sum_{j!=i} E_ij = row sum of W."""
    d = W.sum(1)
    b = d / math.sqrt(max(d.sum(), 1e-300))
    for _ in range(iters):
        nb = np.sqrt(b * d / np.maximum(b.sum() - b, 1e-300))   # damped fixed point
        if np.allclose(nb, b, rtol=1e-12, atol=0):
            b = nb
            break
        b = nb
    E = np.outer(b, b)
    np.fill_diagonal(E, 0)
    return E


# ---------------------------------------------------------------- partition
def make_sets(chroms, config_sets, N):
    """Cannot-link sets as lists of chromosome indices; every chromosome in one set."""
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
    same = z[:, None] == z[None, :]
    return 0.5 * float(B[same].sum())


def _local_search(B, z, sets, N, rng):
    z = z.copy()
    for _ in range(200):
        changed = False
        for si in rng.permutation(len(sets)):
            mem = sets[si]
            oh = np.zeros((B.shape[0], N))
            oh[np.arange(B.shape[0]), z] = 1
            oh[mem] = 0
            A = B[mem] @ oh                      # |set| x N affinity to groups
            if len(mem) == 1:
                g = int(np.argmax(A[0]))
                new = [g]
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
    nvar = V - 1
    total = 1 << nvar
    step = 1 << 16
    base = np.zeros(C)
    for c in sets[0]:
        base[c] = 1.0 if c == sets[0][0] else -1.0
    for lo in range(0, total, step):
        codes = np.arange(lo, min(total, lo + step), dtype=np.int64)
        S = np.tile(base, (codes.size, 1))
        for v in range(1, V):
            bit = ((codes >> (v - 1)) & 1) * 2 - 1
            mem = sets[v]
            S[:, mem[0]] = bit
            if len(mem) == 2:
                S[:, mem[1]] = -bit
        q = np.einsum("ij,ij->i", S @ B, S)
        q[np.abs(S.sum(1)) == C] = -np.inf          # both subgenomes must be non-empty
        j = int(np.argmax(q))
        if q[j] > best:
            best, bz = q[j], S[j].copy()
    return (bz < 0).astype(int)


def _spectral_init(B, N, rng):
    ev, V = np.linalg.eigh(B)
    X = V[:, -max(N - 1, 1):]
    if N == 2:
        return (X[:, -1] < 0).astype(int)
    cen = X[rng.choice(len(X), N, replace=False)]
    for _ in range(50):
        z = np.argmin(((X[:, None, :] - cen[None]) ** 2).sum(-1), 1)
        cen = np.array([X[z == g].mean(0) if (z == g).any() else X[rng.integers(len(X))] for g in range(N)])
    return z


def best_partition(B, N, sets, rng, restarts=30, init=None):
    C = B.shape[0]
    if N == 2 and all(len(s) <= 2 for s in sets) and len(sets) <= 23:
        return _exhaustive2(B, sets)
    inits = [_spectral_init(B, N, rng)]
    if init is not None:
        inits.append(np.asarray(init))
    inits += [rng.integers(0, N, C) for _ in range(restarts)]
    best, bz = -np.inf, None
    for z0 in inits:
        z = _local_search(B, z0, sets, N, rng)
        if len(set(z.tolist())) < N:
            continue
        q = q_score(B, z)
        if q > best + 1e-12:
            best, bz = q, z
    if bz is None:
        bz = _local_search(B, inits[0], sets, N, rng)
    return bz


def canon_labels(z, sizes_key):
    """Relabel groups so that subgenome 1 is the one holding the most chromosome
    mass (ties: lowest index); deterministic output."""
    groups = sorted(set(z.tolist()), key=lambda g: (-sizes_key[z == g].sum(), np.flatnonzero(z == g)[0]))
    m = {g: i for i, g in enumerate(groups)}
    return np.array([m[g] for g in z])


def align(ref, z, N):
    cont = np.zeros((N, N))
    np.add.at(cont, (z, ref), 1)
    r, c = linear_sum_assignment(-cont)
    m = dict(zip(r, c))
    return np.array([m[g] for g in z])


# ---------------------------------------------------------------- number of subgenomes
def ari(a, b):
    """Adjusted Rand index of two labelings (1 = identical, ~0 = chance)."""
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


def bipartition(W, idx, rng):
    """Best 2-way split of chromosomes idx using links among them only (null
    recomputed within the group). Returns 0/1 labels."""
    Wg = W[np.ix_(idx, idx)]
    if Wg.sum() <= 0:
        return np.zeros(len(idx), int)
    B = Wg - expected(Wg)
    return best_partition(B, 2, [[i] for i in range(len(idx))], rng, restarts=10)


def half_pairs(sp, echrom, C, n_el, keep, R, rng, wmode="chrom", mult_pow=0.0):
    """R pairs of link matrices from complementary random halves of the LTR-RTs."""
    out = []
    for _ in range(R):
        u = rng.random(n_el) < 0.5
        out.append((affinity(sp, echrom, C, (u & keep).astype(float), wmode, mult_pow),
                    affinity(sp, echrom, C, (~u & keep).astype(float), wmode, mult_pow)))
    return out


def split_contrast(W, idx, z):
    """Mean log(observed/expected) of links within the two parts minus between
    them, with the null recomputed inside the group."""
    Wg = W[np.ix_(idx, idx)]
    L = np.log((Wg + 1) / (expected(Wg) + 1))
    same = z[:, None] == z[None, :]
    off = ~np.eye(len(idx), dtype=bool)
    return float(L[same & off].mean() - L[~same].mean())


def find_subgenomes(W, pairs, rng, min_size=2, min_rep=0.8, min_rel=0.5, max_n=8):
    """Divisive search for the number of subgenomes. A group is split in two
    when (i) the split is reproducible: the best splits found independently in
    two disjoint halves of the library agree (mean adjusted Rand index over the
    half pairs >= min_rep), and (ii) for groups below the top, the split is at
    least min_rel as strong (within-minus-between log O/E) as the split that
    created the group, so lineage structure inside a subgenome (e.g. a few
    chromosomes sharing centromere-targeting lineages) is not mistaken for
    another subgenome. N=1 (no reproducible structure) is a possible answer.
    Returns labels and the tested splits (for the report)."""
    C = W.shape[0]
    groups = [(np.arange(C), None)]
    final, tests = [], []
    while groups:
        g, parent = groups.pop(0)
        if len(g) < 2 * min_size or len(final) + len(groups) + 1 >= max_n:
            final.append(g)
            continue
        z = bipartition(W, g, rng)
        if min(np.bincount(z, minlength=2)) < min_size:
            tests.append(dict(members=g, split=z, rep=0.0, contrast=0.0, accepted=False))
            final.append(g)
            continue
        con = split_contrast(W, g, z)
        reps = [ari(bipartition(Wa, g, rng), bipartition(Wb, g, rng)) for Wa, Wb in pairs]
        rep = float(np.mean(reps))
        ok = rep >= min_rep and con > 0 and (parent is None or con >= min_rel * parent)
        tests.append(dict(members=g, split=z, rep=rep, contrast=con, parent_contrast=parent, accepted=ok))
        if ok:
            groups += [(g[z == 0], con), (g[z == 1], con)]
        else:
            final.append(g)
    lab = np.zeros(C, int)
    for k, g in enumerate(sorted(final, key=lambda x: x.min())):
        lab[g] = k
    return lab, tests


# ---------------------------------------------------------------- element lineages
def nearest_chromosomes(sp, echrom, C, n_el, top_m, block_nnz=2e7, want_local=False):
    """For every element: the top_m other chromosomes holding its closest
    relatives. Relatedness r_ij = sum over shared splits of mult/(size-1), so
    k-mers shared with few copies weigh most; per chromosome the best relative
    counts once (local copy clusters cannot inflate evidence).
    Returns (E x top_m) chromosome indices (-1 = none) and relatedness; with
    want_local also the number of same-chromosome copies at least as close as
    the best relative elsewhere (local siblings: one cluster, one witness)."""
    S = sp.n
    A = sparse.csr_matrix((np.ones(sp.members.size), sp.members, sp.indptr), shape=(S, n_el))
    d = sp.mult / np.maximum(sp.size - 1, 1)
    DA = sparse.csr_matrix((np.repeat(d, sp.size), sp.members, sp.indptr), shape=(S, n_el))
    AT = A.T.tocsr()                                   # E x S
    cost = np.bincount(sp.members, weights=np.repeat(sp.size, sp.size).astype(float), minlength=n_el)
    top_c = np.full((n_el, top_m), -1, np.int32)
    top_r = np.zeros((n_el, top_m))
    n_local = np.zeros(n_el, np.int64)
    e0 = 0
    cum = np.concatenate(([0], np.cumsum(cost)))
    while e0 < n_el:
        e1 = int(np.searchsorted(cum, cum[e0] + block_nnz, side="right"))
        e1 = min(max(e1, e0 + 1), n_el)
        R = (AT[e0:e1] @ DA).tocoo()
        rows, cols, vals = R.row, R.col, R.data
        ch = echrom[cols]
        own = echrom[e0 + rows]
        k = ch != own
        if want_local:
            bo = np.zeros(e1 - e0)
            np.maximum.at(bo, rows[k], vals[k])
            loc = ~k & (cols != e0 + rows) & (bo[rows] > 0) & (vals >= bo[rows])
            n_local[e0:e1] = np.bincount(rows[loc], minlength=e1 - e0)
        rows, ch, vals = rows[k], ch[k], vals[k]
        if rows.size:
            key = rows.astype(np.int64) * C + ch
            o = np.lexsort((-vals, key))
            key, vals = key[o], vals[o]
            first = np.concatenate(([True], key[1:] != key[:-1]))
            key, vals = key[first], vals[first]          # best relative per (element, chrom)
            r_ = key // C
            c_ = (key % C).astype(np.int32)
            o = np.lexsort((-vals, r_))
            r_, c_, vals = r_[o], c_[o], vals[o]
            start = np.searchsorted(r_, r_, side="left")
            rank = np.arange(r_.size) - start
            k = rank < top_m
            top_c[e0 + r_[k], rank[k]] = c_[k]
            top_r[e0 + r_[k], rank[k]] = vals[k]
        e0 = e1
    return (top_c, top_r, n_local) if want_local else (top_c, top_r)


def relative_informativeness(top_c, top_r, z, own, p0, nbins=20):
    """Self-calibrated weight of each relative, by its closeness relative to the
    copy's closest relative (r/r1): how much more often relatives of that
    closeness sit in the copy's own subgenome than random copies would
    ((q - p0) / (1 - p0), clipped to [0, 1]), measured over the whole library.
    Lineage-level relatives get ~1, family-level relatives ~0, without a fixed
    cut-off and without external truth."""
    valid = top_c >= 0
    ratio = np.where(valid, top_r / np.maximum(top_r[:, :1], 1e-300), np.nan)
    same = valid & (z[np.maximum(top_c, 0)] == own[:, None])
    v = ratio[valid]
    edges = np.unique(np.quantile(v, np.linspace(0, 1, nbins + 1))) if v.size else np.array([0.0, 1.0])
    if edges.size < 2:
        edges = np.array([0.0, 1.0])
    b = np.clip(np.searchsorted(edges, np.nan_to_num(ratio), side="right") - 1, 0, edges.size - 2)
    P0 = np.broadcast_to(p0[:, None], ratio.shape)
    inf = np.zeros(edges.size - 1)
    mid = np.zeros(edges.size - 1)
    for j in range(edges.size - 1):
        sel = valid & (b == j)
        if sel.sum() < 20:
            continue
        q, pn = same[sel].mean(), P0[sel].mean()
        inf[j] = np.clip((q - pn) / max(1 - pn, 1e-9), 0, 1)
        mid[j] = np.median(ratio[sel])
    w = np.where(valid, inf[b], 0.0)
    return w, dict(ratio=mid.tolist(), informativeness=inf.tolist())


def lineage_mixture(top_c, echrom, z, N, nel_chrom, top_r=None, min_rel=0.0, iters=200, eps0=0.05):
    """Empirical-Bayes mixture over elements. Data: the subgenomes of the
    chromosomes holding each element's closest relatives (counts h over N).
    Classes: 'own' (lineage confined to the element's own subgenome),
    'foreign' (confined to another subgenome; one class per other subgenome),
    'shared' (relatives distributed like the elements themselves).
    Returns posteriors (E x 2+N-1 ordered own, shared, foreign_g...), lineage
    subgenome, h counts, and fitted parameters."""
    E = top_c.shape[0]
    C = len(z)
    valid = top_c >= 0
    sgm = np.where(valid, z[np.maximum(top_c, 0)], -1)
    own = z[echrom]
    # null: relatives' subgenome ~ element mass on other chromosomes
    mass = np.array([nel_chrom[z == g].sum() for g in range(N)], float)
    p_null = mass[None, :] - np.eye(N)[own] * nel_chrom[echrom][:, None]
    p_null = p_null / p_null.sum(1, keepdims=True)
    curve = None
    if top_r is None:
        wk = valid.astype(float)
    else:
        wk, curve = relative_informativeness(top_c, top_r, z, own, p_null[np.arange(E), own])
    h = np.stack([((sgm == g) * wk).sum(1) for g in range(N)], 1)
    m = h.sum(1)
    hown = h[np.arange(E), own]
    from scipy.special import gammaln
    lcomb = gammaln(m + 1) - gammaln(h + 1).sum(1)
    pi = np.array([0.3, 0.6, 0.1])  # own, shared, foreign(total)
    eps = eps0
    use = m > 0
    for _ in range(iters):
        ll_own = lcomb + hown * np.log(1 - eps) + (m - hown) * np.log(eps / max(N - 1, 1))
        ll_sh = lcomb + (h * np.log(np.maximum(p_null, 1e-12))).sum(1)
        ll_for = []
        for g in range(N):
            hg = h[:, g]
            ll = lcomb + hg * np.log(1 - eps) + (m - hg) * np.log(eps / max(N - 1, 1))
            ll_for.append(np.where(own == g, -np.inf, ll))
        ll_for = np.stack(ll_for, 1)
        Lik = np.column_stack([ll_own, ll_sh, ll_for])          # data only, no prior
        L = np.column_stack([ll_own + np.log(pi[0]), ll_sh + np.log(pi[1]),
                             ll_for + np.log(pi[2] / max(N - 1, 1))])
        mx = L.max(1, keepdims=True)
        P = np.exp(L - mx)
        P /= P.sum(1, keepdims=True)
        Pu = P[use]
        new_pi = np.array([Pu[:, 0].mean(), Pu[:, 1].mean(), Pu[:, 2:].sum(1).mean()])
        new_pi = np.clip(new_pi, 1e-4, None)
        new_pi /= new_pi.sum()
        # eps: share of relatives outside the confining subgenome, own+foreign classes
        conf = (Pu[:, 0] * (m[use] - hown[use])).sum()
        tot = (Pu[:, 0] * m[use]).sum()
        for g in range(N):
            conf += (Pu[:, 2 + g] * (m[use] - h[use, g])).sum()
            tot += (Pu[:, 2 + g] * m[use]).sum()
        new_eps = float(np.clip(conf / max(tot, 1e-12), 1e-3, 0.3))
        done = np.abs(new_pi - pi).max() < 1e-6 and abs(new_eps - eps) < 1e-6
        pi, eps = new_pi, new_eps
        if done:
            break
    P[~use] = np.nan
    return P, h, m, dict(pi_own=pi[0], pi_shared=pi[1], pi_foreign=pi[2], eps=eps, p_null=p_null, loglik=Lik,
                         relative_curve=curve, rel_w=wk)


def call_lineages(P, own, N, min_post, loglik=None, min_bf=3.0):
    """own | shared | foreign | unresolved, plus lineage subgenome index. A call
    needs posterior >= min_post and, when loglik is given, the copy's own data
    must favour the called class >= min_bf-fold over every alternative (Bayes
    factor 3 = 'positive evidence', Kass & Raftery 1995), so a strong population
    prior alone cannot label an uninformative copy."""
    E = P.shape[0]
    call = np.full(E, "unresolved", dtype=object)
    lin = np.full(E, -1)
    ok = np.isfinite(P[:, 0])
    best = np.where(ok[:, None], P, -1).argmax(1)
    bp = np.where(ok, P[np.arange(E), best], 0)
    sure = ok & (bp >= min_post)
    if loglik is not None:
        Lb = loglik[np.arange(E), best]
        alt = loglik.copy()
        alt[np.arange(E), best] = -np.inf
        sure &= (Lb - alt.max(1)) >= math.log(min_bf)
    call[sure & (best == 0)] = "own"
    lin[sure & (best == 0)] = own[sure & (best == 0)]
    call[sure & (best == 1)] = "shared"
    fo = sure & (best >= 2)
    call[fo] = "foreign"
    lin[fo] = best[fo] - 2
    call[~ok] = "no_relatives"
    return call, lin, bp


# ---------------------------------------------------------------- phasing driver
def chrom_tests(W, E, z, N, pseudo=1.0):
    """Per chromosome: mean O/E to its own subgenome vs the best other one, and
    a one-sided Welch t-test treating the other chromosomes as replicates."""
    from scipy.stats import ttest_ind
    C = len(z)
    L = np.log((W + pseudo) / (E + pseudo))
    OE = W / np.where(E > 0, E, np.nan)
    own_oe, oth_oe, pval, tstat = (np.full(C, np.nan) for _ in range(4))
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
            tstat[c] = t
            pval[c] = p / 2 if t > 0 else 1 - p / 2
    return own_oe, oth_oe, tstat, pval


def phase(lib, echrom, C, N, sets, nmax, wmode, mult_pow, boot, rng, emask=None,
          boot_mode="half", sp=None, init=None):
    """Partition chromosomes; support = fraction of resamples keeping each
    assignment. boot_mode 'half': random 50% of LTR-RTs per replicate (as
    conservative as a half-size library), 'poisson': classic bootstrap.
    sp: precomputed Splits; init: starting labels (e.g. from find_subgenomes)."""
    if sp is None:
        sp = Splits(lib.kmer, lib.eid, echrom, C, nmax, emask)
    if sp.n == 0:
        die("no shared k-mers link two chromosomes; is this an LTR-RT library?")
    eweight = None if emask is None else emask.astype(float)
    W = affinity(sp, echrom, C, eweight, wmode, mult_pow)
    E = expected(W)
    B = W - E
    z = best_partition(B, N, sets, rng, init=init)
    nel = np.bincount(echrom if emask is None else echrom[emask], minlength=C)
    z = canon_labels(z, nel)
    hits = np.zeros((boot, C), bool)
    for r in range(boot):
        if boot_mode == "half":
            wts = (rng.random(lib.n) < 0.5).astype(float)
        else:
            wts = rng.poisson(1.0, lib.n).astype(float)
        if emask is not None:
            wts *= emask
        Wr = affinity(sp, echrom, C, wts, wmode, mult_pow)
        Br = Wr - expected(Wr)
        zr = best_partition(Br, N, sets, rng, restarts=8, init=z)
        hits[r] = align(z, zr, N) == z
    support = hits.mean(0) if boot else np.full(C, np.nan)
    own_oe, oth_oe, tstat, pval = chrom_tests(W, E, z, N)
    return dict(splits=sp, W=W, E=E, B=B, z=z, support=support, boot_hits=hits,
                own_oe=own_oe, oth_oe=oth_oe, tstat=tstat, pval=pval, n_el=nel)


# ---------------------------------------------------------------- dating
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
    R = np.array(R)
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


def age_posterior(k, L, grid, iters=500, tol=1e-7, pi=None):
    """Nonparametric age distribution of copies (Kiefer-Wolfowitz NPMLE on a
    K2P grid) from their LTR substitution counts k ~ Poisson(L d); returns the
    grid weights and each copy's posterior over the grid. No parametric age
    model: bursts, declines and detection bias are absorbed by the weights."""
    from scipy.special import gammaln
    lam = L[:, None] * grid[None, :]
    logP = k[:, None] * np.log(np.maximum(lam, 1e-300)) - lam - gammaln(k + 1)[:, None]
    if grid[0] == 0:
        logP[:, 0] = np.where(k == 0, 0.0, -np.inf)
    P = np.exp(logP - logP.max(1, keepdims=True))
    pi = np.full(grid.size, 1.0 / grid.size) if pi is None else pi.copy()
    for _ in range(iters):
        R = P * pi
        R /= np.maximum(R.sum(1, keepdims=True), 1e-300)
        new = R.mean(0)
        done = np.abs(new - pi).max() < tol
        pi = new
        if done:
            break
    R = P * pi
    R /= np.maximum(R.sum(1, keepdims=True), 1e-300)
    return pi, R


def _merger_sse(post, grid, y, w, s_post):
    """For every candidate tau on the grid: share_i = s_post * P(d_i < tau) +
    sum_{d >= tau} P(d_i = d) (c0 + c1 d); (c0, c1) by weighted least squares,
    so the pre-merger share may drift with age."""
    cp = np.cumsum(post, 1)
    cpd = np.cumsum(post * grid[None, :], 1)
    A = np.concatenate([np.zeros((post.shape[0], 1)), cp[:, :-1]], 1)            # P(d < tau_t)
    B = 1 - A
    Cd = cpd[:, -1:] - np.concatenate([np.zeros((post.shape[0], 1)), cpd[:, :-1]], 1)
    yp = y[:, None] - s_post * A
    W = w[:, None]
    a11, a12, a22 = (W * B * B).sum(0), (W * B * Cd).sum(0), (W * Cd * Cd).sum(0)
    b1, b2 = (W * B * yp).sum(0), (W * Cd * yp).sum(0)
    det = a11 * a22 - a12 ** 2
    det = np.where(np.abs(det) < 1e-12, np.nan, det)
    c0 = (b1 * a22 - b2 * a12) / det
    c1 = (b2 * a11 - b1 * a12) / det
    sse = (W * (yp - B * c0 - Cd * c1) ** 2).sum(0)
    sse[~np.isfinite(sse) | (c0 > s_post)] = np.inf
    return sse, c0, c1


def fit_merger_pooled(k2p, ksub, L, share, m, s_null, trough, rng, boot=50):
    """Merger time from all copies younger than the confined trough at once, at
    the resolution of their pooled LTR sites. Each copy's age posterior comes
    from its own substitution count and LTR length under the nonparametric age
    distribution of all such copies; copies inserted after the merger carry the
    no-confinement share of other-subgenome relatives (s_null, observed wherever
    the post-merger plateau is resolvable); older copies a lower share that may
    drift with age. tau by weighted least squares; 95% CI by bootstrap over
    copies; 'upper' = one-sided 95% bound (useful below resolution)."""
    k = np.where(np.isfinite(ksub), ksub, np.round(np.nan_to_num(k2p) * np.nan_to_num(L, nan=500.0)))
    Lf = np.where(np.isfinite(L) & (L > 0), L, 500.0)
    ok = np.isfinite(k2p) & (k2p <= trough) & (Lf >= 50) & (m > 0) & np.isfinite(share) & (k >= 0)
    if ok.sum() < 100:
        return None
    k, Lf, y, w = k[ok], Lf[ok], share[ok], m[ok]
    top = max(float(np.quantile(k / Lf, 0.999)), 1e-3) * 1.5
    grid = np.concatenate([[0.0], np.geomspace(1e-6, top, 140)])
    pi, post = age_posterior(k, Lf, grid)
    sse, c0, c1 = _merger_sse(post, grid, y, w, s_null)
    j = int(np.argmin(sse))
    sig2 = sse[j] / max(w.sum() - 3, 1) * w.mean()
    upper = float(grid[sse <= sse[j] + 3.84 * sig2].max())
    reps = []
    for _ in range(boot):                       # full re-fit of the age distribution per replicate
        i = rng.integers(0, k.size, k.size)
        _, pb = age_posterior(k[i], Lf[i], grid, iters=500, pi=pi)
        reps.append(int(np.argmin(_merger_sse(pb, grid, y[i], w[i], s_null)[0])))
    if reps:                                    # percentile CI, widened by one grid step (discretization)
        qlo, qhi = np.quantile(reps, [0.025, 0.975])
        lo, hi = grid[max(int(np.floor(qlo)) - 1, 0)], grid[min(int(np.ceil(qhi)) + 1, grid.size - 1)]
    else:
        lo = hi = np.nan
    return dict(tau=float(grid[j]), ci=[float(lo), float(hi)], upper=upper, share_post=float(s_null),
                share_confined=float(c0[j]), share_slope=float(c1[j]), n=int(k.size),
                resolution=float(1.0 / np.median(Lf)))


def age_informativeness(k2p, fit, s_null):
    """Weight in [0,1] of each copy as evidence about which progenitor a segment
    came from, read off the observed share-vs-age curve: 1 at ages where
    lineages are most confined to one subgenome, 0 where copies' relatives are
    spread like random copies (post-merger insertions, ancestral lineages).
    Copies without K2P get weight 1."""
    b = fit.get("bins") if fit else None
    if not b or len(b["k2p"]) < 3:
        return None
    A, Y = np.array(b["k2p"]), np.array(b["share"])
    o = np.argsort(A)
    A, Y = A[o], Y[o]
    ymin = Y.min()
    if s_null - ymin < 0.05:                      # no age-dependent confinement to exploit
        return None
    sh = np.interp(np.nan_to_num(k2p, nan=-1.0), A, Y)
    w = np.clip((s_null - sh) / (s_null - ymin), 0.0, 1.0)
    w[~np.isfinite(k2p)] = 1.0
    return w


def k2p_se_guess(k2p, ltr_len=None):
    """K2P standard error when the table lacks one: binomial on the LTR length
    (default 500 bp), floored at one substitution."""
    L = np.where(np.isfinite(ltr_len), ltr_len, 500.0) if ltr_len is not None else 500.0
    p = np.maximum(np.nan_to_num(k2p, nan=0.0), 1.0 / L)
    return np.sqrt(p * (1 - p) / L)


# ---------------------------------------------------------------- exchange HMM
def _isotonic(y, w):
    """Weighted pool-adjacent-violators: the non-decreasing sequence closest to y."""
    vals, wts, cnt = [], [], []
    for yi, wi in zip(y, np.maximum(w, 1e-12)):
        vals.append(yi)
        wts.append(wi)
        cnt.append(1)
        while len(vals) > 1 and vals[-2] > vals[-1]:
            w2 = wts[-2] + wts[-1]
            v2 = (vals[-2] * wts[-2] + vals[-1] * wts[-1]) / w2
            c2 = cnt[-2] + cnt[-1]
            vals[-2:], wts[-2:], cnt[-2:] = [v2], [w2], [c2]
    return np.repeat(vals, cnt)


def lineage_model(top_c, top_r, echrom, z, N, nel_chrom, nbins=12, iters=500, tol=1e-8):
    """Two-level mixture over copies, fitted by EM (no thresholds, no truth).
    Data per copy: for every other chromosome, the closeness x of its best
    relative there (x = r / r_closest; 0 if none). A copy's lineage is either
    confined to one subgenome g - then closeness on g's chromosomes follows
    f_in(x) and on the others f_out(x) - or shared, with f_sh(x) everywhere.
    The densities (histograms over x) and class priors are learned; the
    likelihood ratio f_in/f_out is constrained to be non-decreasing in x (a
    closer relative never argues less for the same lineage), which removes the
    mirror solution. Modelling closeness rather than counts per subgenome makes
    the model indifferent to how many chromosomes or LTR-RTs each subgenome has.
    Returns per-copy class log-likelihoods (columns: confined to g, shared),
    priors, the learned curves and per-relative lineage weights."""
    E = top_c.shape[0]
    C = len(z)
    own = z[echrom]
    X = np.zeros((E, C))
    rows = np.repeat(np.arange(E), top_c.shape[1])
    tc, tr = top_c.ravel(), top_r.ravel()
    ok = tc >= 0
    r1 = np.maximum(top_r[:, 0], 1e-300)
    X[rows[ok], tc[ok]] = tr[ok] / r1[rows[ok]]
    mask = np.ones((E, C), bool)
    mask[np.arange(E), echrom] = False               # own chromosome excluded
    xv = X[mask & (X > 0)]
    q = np.unique(np.quantile(xv, np.linspace(0, 1, nbins + 1))) if xv.size else np.array([0.0, 1.0])
    edges = np.concatenate([[0.0], q[1:-1], [np.inf]]) if q.size > 2 else np.array([0.0, np.inf])
    Bn = edges.size                                  # bin 0 = no relative; 1.. = closeness bins
    bidx = np.where(X > 0, np.clip(np.searchsorted(edges, X, side="right"), 1, Bn - 1), 0)
    ing = np.stack([z[None, :] == g for g in range(N)], 0)          # N x 1 x C (broadcast)
    cnt_all = np.zeros((E, Bn))
    cnt_in = np.zeros((N, E, Bn))
    for bb in range(Bn):
        hit = (bidx == bb) & mask
        cnt_all[:, bb] = hit.sum(1)
        for g in range(N):
            cnt_in[g, :, bb] = (hit & ing[g]).sum(1)
    use = (mask & (X > 0)).any(1)
    tot = cnt_all.sum(0) + 1.0
    f_sh = tot / tot.sum()
    f_in = f_sh * np.linspace(0.5, 1.5, Bn)
    f_in /= f_in.sum()
    f_out = f_sh * np.linspace(1.5, 0.5, Bn)
    f_out /= f_out.sum()
    pi_own, pi_sh, pi_for = 0.3, 0.6, 0.1
    for _ in range(iters):
        lin, lout, lsh = np.log(f_in), np.log(f_out), np.log(f_sh)
        L = np.empty((E, N + 1))
        for g in range(N):
            L[:, g] = cnt_in[g] @ lin + (cnt_all - cnt_in[g]) @ lout
        L[:, N] = cnt_all @ lsh
        prior = np.empty((E, N + 1))
        for g in range(N):
            prior[:, g] = np.where(own == g, pi_own, pi_for / max(N - 1, 1))
        prior[:, N] = pi_sh
        lp = L + np.log(np.maximum(prior, 1e-300))
        lp -= lp.max(1, keepdims=True)
        R = np.exp(lp)
        R /= R.sum(1, keepdims=True)
        R[~use] = 0
        new_own = R[np.arange(E), own][use].mean()
        new_sh = R[use, N].mean()
        new_for = max(1.0 - new_own - new_sh, 1e-9)
        n_in = sum(R[:, g] @ cnt_in[g] for g in range(N)) + 0.5
        n_out = sum(R[:, g] @ (cnt_all - cnt_in[g]) for g in range(N)) + 0.5
        n_sh = R[:, N] @ cnt_all + 0.5
        g_out = n_out / n_out.sum()
        lr = np.log((n_in / n_in.sum()) / g_out)
        lr = _isotonic(lr, n_in + n_out)                 # monotone likelihood ratio in closeness
        g_in = g_out * np.exp(lr)
        g_in /= g_in.sum()
        g_sh = n_sh / n_sh.sum()
        done = (np.abs(g_in - f_in).max() < tol and np.abs(g_out - f_out).max() < tol
                and abs(new_own - pi_own) < tol and abs(new_sh - pi_sh) < tol)
        f_in, f_out, f_sh = g_in, g_out, g_sh
        pi_own, pi_sh, pi_for = new_own, new_sh, new_for
        if done:
            break
    # per-relative lineage weight: posterior that a relative at this closeness is a lineage relative
    wbin = np.clip(1.0 - f_out / f_in, 0.0, 1.0)
    wbin[0] = 0.0
    W = np.where(mask, wbin[bidx], 0.0)
    mids = [0.0] + [float(np.median(X[mask & (bidx == j)])) if (mask & (bidx == j)).any() else float("nan")
                    for j in range(1, Bn)]
    return dict(L=L, use=use, pi_own=float(pi_own), pi_shared=float(pi_sh), pi_foreign=float(pi_for),
                a=wbin.tolist(), a_closeness=mids, W=W, f_in=f_in.tolist(), f_out=f_out.tolist(),
                f_shared=f_sh.tolist())


def lineage_calls(model, echrom, z, N, min_post, min_bf=3.0):
    """own | shared | foreign | unresolved | no_relatives from the lineage model;
    posterior >= min_post and Bayes factor >= min_bf (copy's own data) vs every
    alternative class."""
    E = model["L"].shape[0]
    own = z[echrom]
    L = model["L"]
    # class order for output: own, shared, foreign_g
    LL = np.column_stack([L[np.arange(E), own], L[:, N]] + [np.where(own == g, -np.inf, L[:, g]) for g in range(N)])
    pri = np.array([model["pi_own"], model["pi_shared"]] + [model["pi_foreign"] / max(N - 1, 1)] * N)
    lp = LL + np.log(np.maximum(pri, 1e-300))[None, :]
    lp -= lp.max(1, keepdims=True)
    P = np.exp(lp)
    P /= P.sum(1, keepdims=True)
    P[~model["use"]] = np.nan
    call, lin, bp = call_lineages(P, own, N, min_post, LL, min_bf)
    return P, call, lin, bp


def model_emissions(model, N):
    """HMM emission for origin g: the copy's lineage is confined to g, shared,
    or confined to another subgenome (a transposed copy), weighted by the
    fitted class priors. The shared alternative caps any single copy's weight."""
    from scipy.special import logsumexp
    L = model["L"]
    E = L.shape[0]
    po, ps, pf = model["pi_own"], model["pi_shared"], model["pi_foreign"]
    em = np.empty((E, N))
    for g in range(N):
        others = [L[:, x] for x in range(N) if x != g]
        lo = logsumexp(np.stack(others, 1), 1) - math.log(max(N - 1, 1)) if others else np.full(E, -np.inf)
        em[:, g] = logsumexp(np.stack([math.log(max(po, 1e-6)) + L[:, g], math.log(max(ps, 1e-6)) + L[:, N],
                                       math.log(max(pf, 1e-6)) + lo], 1), 1)
    em[~model["use"]] = 0.0
    return em


def relative_emissions(top_c, rel_w, z, p_null, N):
    """Per copy, log-likelihood of each origin subgenome g from its relatives:
    a relative of measured informativeness I (relative_informativeness) sits in
    the copy's origin subgenome with probability I and is otherwise a random
    copy (null share), so P(relative in s | origin g) = I [s = g] + (1 - I) p0(s).
    Uninformative relatives cost nothing; informative ones count fully."""
    E = top_c.shape[0]
    valid = top_c >= 0
    s = np.where(valid, z[np.maximum(top_c, 0)], 0)
    p0 = p_null[np.arange(E)[:, None], s]                   # E x K
    em = np.zeros((E, N))
    for g in range(N):
        lik = rel_w * (s == g) + (1 - rel_w) * p0
        em[:, g] = np.where(valid, np.log(np.maximum(lik, 1e-12)), 0.0).sum(1)
    return em


def segment_hmm(order_by_chrom, em, z, N, switch=2e-3, weight=None):
    """Per chromosome, a hidden 'origin subgenome' runs along the elements in
    positional order; em = per-copy log-likelihood of each origin (from
    relative_emissions). Returns per element the Viterbi state, its posterior and
    the (tempered) emissions."""
    from scipy.special import logsumexp
    E = em.shape[0]
    em = em - em.max(1, keepdims=True)
    if weight is not None:          # temper correlated copies (local clusters, uninformative age classes)
        em = em * weight[:, None]
    lT = np.full((N, N), np.log(switch / max(N - 1, 1)))
    np.fill_diagonal(lT, np.log(1 - switch))
    state = np.full(E, -1)
    post = np.full(E, np.nan)
    for c, idx in order_by_chrom.items():
        if idx.size == 0:
            continue
        e = em[idx]
        n = idx.size
        l0 = np.full(N, np.log(0.1 / max(N - 1, 1)))
        l0[z[c]] = np.log(0.9) if N > 1 else 0.0
        fw = np.empty((n, N))
        fw[0] = l0 + e[0]
        for t in range(1, n):
            fw[t] = logsumexp(fw[t - 1][:, None] + lT, 0) + e[t]
        bw = np.zeros((n, N))
        for t in range(n - 2, -1, -1):
            bw[t] = logsumexp(lT + (e[t + 1] + bw[t + 1])[None, :], 1)
        pp = fw + bw
        pp = np.exp(pp - logsumexp(pp, 1, keepdims=True))
        # Viterbi
        vt = np.empty((n, N))
        bp = np.zeros((n, N), int)
        vt[0] = l0 + e[0]
        for t in range(1, n):
            sc = vt[t - 1][:, None] + lT
            bp[t] = sc.argmax(0)
            vt[t] = sc.max(0) + e[t]
        path = np.empty(n, int)
        path[-1] = int(vt[-1].argmax())
        for t in range(n - 1, 0, -1):
            path[t - 1] = bp[t, path[t]]
        state[idx] = path
        post[idx] = pp[np.arange(n), path]
    return state, post, em


def segments_from_states(order_by_chrom, state, post, em, call, lin, start, end, z, chroms):
    """Runs of one HMM state per chromosome. llr = summed per-copy log-likelihood
    of the run's origin subgenome over the chromosome's own subgenome (natural
    log; 0 for runs matching the chromosome)."""
    rows = []
    confined = (call == "own") | (call == "foreign")
    for c, idx in order_by_chrom.items():
        if idx.size == 0:
            continue
        st = state[idx]
        brk = np.flatnonzero(np.diff(st)) + 1
        for s, e in zip(np.concatenate(([0], brk)), np.concatenate((brk, [idx.size]))):
            ii = idx[s:e]
            g = int(st[s])
            rows.append(dict(chrom=chroms[c], start=int(start[ii].min()), end=int(end[ii].max()),
                             origin=g, chrom_sg=int(z[c]), n_elements=int(ii.size),
                             n_support=int((confined[ii] & (lin[ii] == g)).sum()),
                             n_conflict=int((confined[ii] & (lin[ii] != g) & (lin[ii] >= 0)).sum()),
                             llr=float((em[ii, g] - em[ii, z[c]]).sum()),
                             mean_post=float(np.nanmean(post[ii]))))
    return rows


# ---------------------------------------------------------------- genome painting (optional)
def subgenome_markers(lib, echrom, keep, z, N, nmax):
    """k-mers (sampled hashes) carried by 2..nmax LTR-RTs that sit on >= 2
    chromosomes, all of one subgenome: markers of lineages confined to it."""
    sel = keep[lib.eid]
    km, ei = lib.kmer[sel], lib.eid[sel]
    brk = np.flatnonzero(km[1:] != km[:-1]) + 1
    st = np.concatenate(([0], brk))
    size = np.diff(np.concatenate((st, [km.size])))
    ch = echrom[ei]
    sg = z[ch]
    ok = (size >= 2) & (size <= nmax)
    one_sg = np.minimum.reduceat(sg, st) == np.maximum.reduceat(sg, st)
    multi = np.minimum.reduceat(ch, st) != np.maximum.reduceat(ch, st)
    g = sg[st]
    m = ok & one_sg & multi
    return [np.sort(km[st][m & (g == j)]) for j in range(N)]


def _sampled(codes, k, scale):
    n = codes.size - k + 1
    if n <= 0:
        return np.empty(0, np.int64), np.empty(0, np.uint64)
    bad = codes == 4
    cs = np.zeros(codes.size + 1, np.int32)
    np.cumsum(bad, out=cs[1:])
    ok = (cs[k:k + n] - cs[:n]) == 0
    c = codes.astype(np.uint64)
    c[bad] = 0
    f = np.zeros(n, np.uint64)
    r = np.zeros(n, np.uint64)
    for j in range(k):
        cj = c[j:j + n]
        f <<= np.uint64(2)
        f |= cj
        r |= (np.uint64(3) - cj) << np.uint64(2 * j)
    h = mix64(np.minimum(f, r))
    sel = ok & ((h % np.uint64(scale)) == 0)
    return np.flatnonzero(sel), h[sel] >> _SH


def genome_paint(paths, chroms, markers, k, scale, win, chunk=1 << 23, min_len=None):
    """Counts of each subgenome's marker k-mers per window along the genome
    (fragmented and solo LTRs included). Paints every sequence >= min_len (all
    of chroms regardless). Returns {name: (starts, counts NxW)}."""
    want = set(chroms)
    min_len = win if min_len is None else min_len
    allm = np.concatenate(markers)
    lab = np.concatenate([np.full(m.size, j) for j, m in enumerate(markers)])
    o = np.argsort(allm)
    allm, lab = allm[o], lab[o]
    out = {}
    for path in paths:
        for name, sq in read_fasta(path):
            if name not in want and len(sq) < min_len:
                continue
            codes = LUT[np.frombuffer(sq, np.uint8)]
            nw = codes.size // win + 1
            cnt = np.zeros((nw, len(markers)), np.int64)
            for s0 in range(0, codes.size, chunk):
                pos, h = _sampled(codes[s0:s0 + chunk + k - 1], k, scale)
                i = np.searchsorted(allm, h)
                i = np.minimum(i, allm.size - 1)
                hit = allm[i] == h
                np.add.at(cnt, ((pos[hit] + s0) // win, lab[i[hit]]), 1)
            out[name] = (np.arange(nw) * win, cnt)
            vlog(f"genome: {name} {codes.size / 1e6:.1f} Mb, {int(cnt.sum())} marker hits")
    return out


def window_votes(cnt, min_markers=3):
    """Per window: the subgenome holding most marker k-mers (-1 if < min_markers or tied)."""
    tot = cnt.sum(1)
    tie = (cnt == cnt.max(1, keepdims=True)).sum(1) > 1
    return np.where((tot >= min_markers) & ~tie, cnt.argmax(1), -1)


def assign_sequences(paint, N, min_markers=3):
    """Subgenome of every painted sequence from its window votes: majority
    subgenome, its share of informative windows, and a binomial p-value against
    an even split (useful for scaffolds without intact LTR-RTs)."""
    from scipy.stats import binomtest
    out = {}
    for name, (st, cnt) in paint.items():
        v = window_votes(cnt, min_markers)
        n = np.bincount(v[v >= 0], minlength=N)
        tot = int(n.sum())
        if tot == 0:
            out[name] = (-1, n, np.nan, np.nan)
            continue
        g = int(n.argmax())
        p = binomtest(int(n[g]), tot, 1.0 / N, alternative="greater").pvalue
        out[name] = (g, n, n[g] / tot, p)
    return out


def genome_segments(paint, chroms, z, N, min_win=1, switch=1e-3, min_markers=3, min_llr=math.log(100)):
    """HMM along windows. Marker k-mers in one window are not independent (one
    copy or fragment carries many), so each window casts one vote: the subgenome
    with most markers (>= min_markers; ties abstain). Vote error rates are
    estimated from the windows of chromosomes assigned to each subgenome; a run
    of the other origin is reported when its likelihood ratio is decisive
    (>= 100, Kass & Raftery 1995)."""
    cix = {c: i for i, c in enumerate(chroms)}
    votes = {}
    conf = np.ones((N, N))                       # pseudocount
    for c, (_, cnt) in paint.items():
        votes[c] = window_votes(cnt, min_markers)
        if c not in cix:
            continue
        g0 = z[cix[c]]
        for j in range(N):
            conf[g0, j] += (votes[c] == j).sum()
    lE = np.log(conf / conf.sum(1, keepdims=True))   # log P(vote j | origin o) = lE[o, j]
    lT = np.full((N, N), math.log(switch / max(N - 1, 1)))
    np.fill_diagonal(lT, math.log(1 - switch))
    rows, states = [], {}
    for c, (st, cnt) in paint.items():
        if c not in cix:
            continue
        g0 = z[cix[c]]
        v = votes[c]
        em = np.where(v[:, None] >= 0, lE[:, np.maximum(v, 0)].T, 0.0)   # W x N
        n = em.shape[0]
        vt = np.empty((n, N))
        bp = np.zeros((n, N), int)
        l0 = np.full(N, math.log(0.01 / max(N - 1, 1)))
        l0[g0] = math.log(0.99)
        vt[0] = l0 + em[0]
        for t in range(1, n):
            sc = vt[t - 1][:, None] + lT
            bp[t] = sc.argmax(0)
            vt[t] = sc.max(0) + em[t]
        path = np.empty(n, int)
        path[-1] = int(vt[-1].argmax())
        for t in range(n - 1, 0, -1):
            path[t - 1] = bp[t, path[t]]
        states[c] = path
        step = st[1] - st[0] if n > 1 else 0
        brk = np.flatnonzero(np.diff(path)) + 1
        for s0, e0 in zip(np.concatenate(([0], brk)), np.concatenate((brk, [n]))):
            g = int(path[s0])
            if g == g0:
                continue
            agree = int((v[s0:e0] == g).sum())
            llr = float((em[s0:e0, g] - em[s0:e0, g0]).sum())
            rows.append(dict(chrom=c, start=int(st[s0]), end=int(st[e0 - 1] + step), origin=g, chrom_sg=int(g0),
                             n_windows=int(e0 - s0), markers=int(cnt[s0:e0].sum()), agree=agree, llr=llr,
                             flag=agree >= min_win and llr >= min_llr))
    return rows, states


# ---------------------------------------------------------------- figures
PALETTE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000"]  # Okabe-Ito
GREY = "#9E9E9E"


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


def _panel(ax, letter):
    ax.text(-0.02, 1.02, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", ha="right", va="bottom")


def _save(fig, out):
    fig.savefig(out + ".pdf")
    fig.savefig(out + ".png", dpi=300)


def plot_phasing(out, chroms, z, support, W, E, own_oe, oth_oe, N, sgn):
    plt = _plt()
    from matplotlib.colors import ListedColormap
    order = sorted(range(len(chroms)), key=lambda i: (z[i], natural_key(chroms[i])))
    C = len(order)
    L = np.log2((W + 1) / (E + 1))[np.ix_(order, order)]
    np.fill_diagonal(L, np.nan)
    fig = plt.figure(figsize=(7.1, 3.4))
    ax = fig.add_axes([0.07, 0.12, 0.38, 0.38 * 7.1 / 3.4 * 0.95])
    lim = max(0.5, float(np.nanpercentile(np.abs(L), 98)))
    im = ax.imshow(L, cmap="RdBu_r", vmin=-lim, vmax=lim, interpolation="nearest")
    fs = 6 if C <= 30 else 4
    ax.set_xticks(range(C))
    ax.set_xticklabels([chroms[i] for i in order], rotation=90, fontsize=fs)
    ax.set_yticks(range(C))
    ax.set_yticklabels([chroms[i] for i in order], fontsize=fs)
    for sp_ in ax.spines.values():
        sp_.set_visible(False)
    ax.tick_params(length=0)
    cmap = ListedColormap(PALETTE[:N])
    for g in range(N):
        idx = [k for k, i in enumerate(order) if z[i] == g]
        if idx:
            ax.add_patch(plt.Rectangle((min(idx) - 0.5, -2.2), len(idx), 1.2, color=PALETTE[g], clip_on=False))
            ax.text((min(idx) + max(idx)) / 2, -2.6, sgn[g], ha="center", va="bottom", fontsize=7, color=PALETTE[g], fontweight="bold")
    ax.set_xlim(-0.5, C - 0.5)
    ax.set_ylim(C - 0.5, -2.4)
    cax = fig.add_axes([0.465, 0.12, 0.01, 0.3])
    cb = fig.colorbar(im, cax=cax)
    ticks = [t for t in (-2, -1, 0, 1, 2) if abs(t) <= lim]
    cb.set_ticks(ticks)
    cb.set_ticklabels([f"{2.0 ** t:g}" for t in ticks])
    cb.set_label("links between chromosomes,\nobserved / expected", fontsize=6)
    cb.ax.tick_params(labelsize=6)
    _panel(ax, "A")
    ax2 = fig.add_axes([0.62, 0.14, 0.35, 0.78])
    lo = np.nanmin(np.concatenate([own_oe, oth_oe]))
    hi = np.nanmax(np.concatenate([own_oe, oth_oe]))
    pad = 0.05 * (hi - lo + 1e-9)
    ax2.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color=GREY, lw=0.6, ls="--", zorder=0)
    for i in range(len(chroms)):
        if not np.isfinite(own_oe[i]):
            continue
        col = PALETTE[z[i]]
        solid = support[i] >= 0.95 if np.isfinite(support[i]) else True
        ax2.scatter(oth_oe[i], own_oe[i], s=16, facecolor=col if solid else "white", edgecolor=col, lw=0.8, zorder=3)
        if not solid:
            ax2.annotate(chroms[i], (oth_oe[i], own_oe[i]), xytext=(3, -1), textcoords="offset points",
                         fontsize=5.5, color="#444444")
    ax2.set_xlabel("mean observed/expected links to the\nclosest other subgenome")
    ax2.set_ylabel("mean observed/expected links to own subgenome")
    ax2.set_xlim(lo - pad, hi + pad)
    ax2.set_ylim(lo - pad, hi + pad)
    ax2.text(0.97, 0.03, "filled: support ≥ 0.95\nopen (named): support < 0.95", transform=ax2.transAxes,
             ha="right", va="bottom", fontsize=5.5, color="#444444")
    for g in range(N):
        mm = (z == g) & np.isfinite(own_oe)
        if mm.any():
            ax2.annotate(f"{sgn[g]} (n={int(mm.sum())})", (np.nanmax(oth_oe[mm]), np.nanmedian(own_oe[mm])),
                         xytext=(10, 0), textcoords="offset points", color=PALETTE[g], fontsize=7,
                         fontweight="bold", va="center")
    _panel(ax2, "B")
    _save(fig, out)
    plt.close(fig)


def plot_painting(out, chroms, z, support, clen, order_by_chrom, h, state, start, end, segs, N, sgn, purity=0.75,
                  gstates=None, gsegs=None, win=None):
    """Every LTR-RT with relatives on other chromosomes, drawn at its position and
    coloured by the subgenome holding >= purity of its closest relatives (grey:
    mixed); bar under each chromosome: HMM origin; black boxes: candidate exchanges."""
    plt = _plt()
    rows = sorted([c for c in order_by_chrom if order_by_chrom[c].size], key=lambda c: (z[c], natural_key(chroms[c])))
    H = 0.24 * len(rows) + 0.75
    fig, ax = plt.subplots(figsize=(7.1, H))
    fig.subplots_adjust(left=0.13, right=0.98, top=1 - 0.5 / H, bottom=0.35 / H)
    maxlen = max(clen[c] for c in rows) / 1e6
    m = h.sum(1)
    frac = h / np.maximum(m, 1e-12)[:, None]
    point = np.where((m > 0) & (frac.max(1) >= purity), frac.argmax(1), -1)
    cols = np.array(PALETTE[:N] + ["#C8C8C8"])
    for k, c in enumerate(rows):
        y = -k
        idx = order_by_chrom[c]
        ax.plot([0, clen[c] / 1e6], [y, y], color="#9E9E9E", lw=0.5, solid_capstyle="butt", zorder=1)
        st = state[idx]
        brk = np.flatnonzero(np.diff(st)) + 1
        for s0, e0 in zip(np.concatenate(([0], brk)), np.concatenate((brk, [idx.size]))):
            x0 = start[idx[s0]] / 1e6 if s0 > 0 else 0.0
            x1 = end[idx[e0 - 1]] / 1e6 if e0 < idx.size else clen[c] / 1e6
            ax.add_patch(plt.Rectangle((x0, y - 0.2), max(x1 - x0, 1e-3), 0.09, color=PALETTE[st[s0]], lw=0, zorder=1))
        for sg in segs:
            if sg["chrom"] == chroms[c] and sg["flag"]:
                ax.add_patch(plt.Rectangle((sg["start"] / 1e6, y - 0.27), (sg["end"] - sg["start"]) / 1e6, 0.62,
                                           fill=False, edgecolor="black", lw=0.7, zorder=4))
        if gstates and chroms[c] in gstates:
            gp = gstates[chroms[c]]
            gb = np.flatnonzero(np.diff(gp)) + 1
            for s0, e0 in zip(np.concatenate(([0], gb)), np.concatenate((gb, [gp.size]))):
                ax.add_patch(plt.Rectangle((s0 * win / 1e6, y - 0.33), (e0 - s0) * win / 1e6, 0.09,
                                           color=PALETTE[gp[s0]], lw=0, zorder=1, alpha=0.75))
            for g in gsegs or []:
                if g["chrom"] == chroms[c] and g["flag"]:
                    ax.add_patch(plt.Rectangle((g["start"] / 1e6, y - 0.37), (g["end"] - g["start"]) / 1e6, 0.72,
                                               fill=False, edgecolor="black", lw=0.6, ls=(0, (2, 1.2)), zorder=4))
        has = m[idx] > 0
        ii = idx[has]
        x = (start[ii] + end[ii]) / 2e6
        pt = point[ii]
        grey = pt < 0
        ax.vlines(x[grey], y + 0.03, y + 0.27, colors=cols[-1], lw=0.3, zorder=2)
        ax.vlines(x[~grey], y + 0.03, y + 0.27, colors=cols[pt[~grey]], lw=0.35, zorder=3)
        sup = f"{support[c]:.2f}" if np.isfinite(support[c]) else ""
        ax.text(-0.01 * maxlen, y, f"{chroms[c]}  {sup}", ha="right", va="center", fontsize=5.5, color=PALETTE[z[c]])
    ax.set_xlim(-0.005 * maxlen, maxlen * 1.01)
    ax.set_ylim(-len(rows) + 0.5, 1.05)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("position (Mb)")
    hx = 0.0
    for col, lab in [(PALETTE[g], f"relatives in {sgn[g]}") for g in range(N)] + [("#C8C8C8", "mixed")]:
        ax.text(hx * maxlen, 0.78, "|", color=col, fontsize=8, fontweight="bold", va="center")
        ax.text(hx * maxlen + 0.012 * maxlen, 0.78, lab, fontsize=6, va="center")
        hx += 0.16
    ax.text(hx * maxlen, 0.78, f"ticks: LTR-RTs (\u2265{purity:.0%} of closest relatives in one subgenome)\n"
            + ("bars: origin (upper LTR-RTs, lower genome); boxes: exchanges (dashed: genome)"
               if gstates else "bar: inferred origin (HMM); box: candidate exchange; number: support"),
            fontsize=5.5, va="center", color="#444444")
    _save(fig, out)
    plt.close(fig)


def plot_ages(out, k2p, call, lin, mu, fit, share, null_share):
    plt = _plt()
    ok = np.isfinite(k2p)
    hi = float(np.nanquantile(k2p[ok], 0.99)) if ok.any() else 0.1
    hi = max(hi, 1e-3)
    lin_t = 0.002
    fig, axs = plt.subplots(1, 2, figsize=(7.1, 2.7))
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.19, top=0.8, wspace=0.32)
    ax = axs[0]
    bins = np.linspace(0, hi, 41)
    ny = 0
    for kind, col in (("own", "#009E73"), ("shared", GREY), ("foreign", "#CC79A7")):
        v = k2p[ok & (call == kind)]
        if v.size < 5:
            continue
        hcount, _ = np.histogram(np.clip(v, 0, hi), bins)
        ax.step(bins[:-1], hcount / v.size, where="post", color=col, lw=1.0)
        ax.text(0.97, 0.95 - 0.09 * ny, f"{kind} lineage (n={v.size:,})", transform=ax.transAxes, color=col,
                fontsize=6, ha="right", va="top")
        ny += 1
    ax.set_xlim(0, hi)
    ax.set_xlabel("LTR divergence (K2P)")
    ax.set_ylabel("fraction of copies per bin")
    sec = ax.secondary_xaxis("top", functions=(lambda x: x / (2 * mu) / 1e6, lambda t: t * 2 * mu * 1e6))
    sec.set_xlabel(f"insertion time (Myr, \u03bc = {mu:g})", fontsize=6)
    _panel(ax, "A")
    ax = axs[1]
    if fit and fit.get("bins"):
        A, Y, Wt = (np.array(fit["bins"][k]) for k in ("k2p", "share", "weight"))
        ax.scatter(A, Y, s=6 + 30 * Wt / Wt.max(), color="#555555", zorder=3, lw=0)
        ax.axhline(null_share, color=GREY, ls="--", lw=0.6)
        ax.text(hi, null_share, "no subgenome confinement", fontsize=5.5, color=GREY, ha="right", va="bottom")
        pm = fit.get("pooled")
        if pm:
            xm = pm["tau"]
            lo_, hi_ = pm["ci"]
            if xm > 0:
                ax.axvspan(lo_, max(hi_, lo_ + 1e-7), color="black", alpha=0.12, lw=0)
                ax.axvline(xm, color="black", lw=1.0)
                lab = (f"merger {xm / (2 * mu) / 1e3:,.1f} kyr\n(95% CI {lo_ / (2 * mu) / 1e3:,.1f}\u2013"
                       f"{hi_ / (2 * mu) / 1e3:,.1f})")
            else:
                ax.axvspan(0, pm["upper"], color="black", alpha=0.12, lw=0)
                lab = f"merger < {pm['upper'] / (2 * mu) / 1e3:,.1f} kyr\n(95% upper bound)"
            ax.text(max(xm, 0), 1.02, lab, color="black", fontsize=6, ha="left", va="bottom",
                    transform=ax.get_xaxis_transform())
        ax.set_xscale("symlog", linthresh=lin_t, linscale=0.6)
        ax.set_xlim(0, hi)
        ax.set_ylim(0, max(0.7, float(np.nanmax(Y)) * 1.1))
    else:
        ax.text(0.5, 0.5, "too few dated LTR-RTs\nfor the transposition clock", transform=ax.transAxes,
                ha="center", va="center", fontsize=6.5, color="#444444")
    ax.set_xlabel(f"LTR divergence (K2P; linear below {lin_t:g}, log above)")
    ax.set_ylabel("share of closest relatives\non other-subgenome chromosomes")
    _panel(ax, "B")
    _save(fig, out)
    plt.close(fig)


# ---------------------------------------------------------------- main
def main(argv=None):
    global VERBOSE
    ap = argparse.ArgumentParser(
        description="Config-free subgenome phasing of an allopolyploid from its LTR-RT library (v4).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--ltr_fasta", required=True, help="LTR-RT FASTA(.gz); headers chrom:start-end[#Class/Superfamily/Family]")
    ap.add_argument("--outdir", required=True, help="output directory")
    ap.add_argument("-n", "--n_subgenomes", type=int, default=None, help="number of subgenomes (default: inferred from the data by reproducible divisive splitting; with --config: config width)")
    ap.add_argument("--min_rep", type=float, default=0.8, help="auto N: split-half replicability (adjusted Rand index) needed to accept a split")
    ap.add_argument("--auto_reps", type=int, default=20, help="auto N: number of disjoint half-library pairs")
    ap.add_argument("--config", default=None, help="optional: whitespace-separated sets of chromosomes known to be in different subgenomes (one set per line)")
    ap.add_argument("--k2p", default=None, help="optional: LTR divergence table (LTRquest/Kmer2LTR TSV with seq_id,k2p[,k2p_se,ltr5_len], or id + K2P in --k2p_col)")
    ap.add_argument("--k2p_col", type=int, default=None, help="1-based K2P column for header-less tables (e.g. 11)")
    ap.add_argument("--mu", type=float, default=1.3e-8, help="substitutions/site/year for K2P -> years")
    ap.add_argument("--genome_fai", default=None, help="optional genome .fai for chromosome lengths in plots")
    ap.add_argument("--genome", nargs="+", default=None, help="optional polyploid genome FASTA(s): paint every window with subgenome-marker k-mers (fragmented/solo LTRs included) and call exchanges genome-wide")
    ap.add_argument("--win", type=int, default=50_000, help="window size (bp) for --genome painting")
    ap.add_argument("-k", type=int, default=21, help="k-mer size (<=31)")
    ap.add_argument("--scale", type=int, default=0, help="keep 1/scale of k-mers (FracMinHash); 0 = auto (4, more for very large libraries)")
    ap.add_argument("--nmax", type=int, default=0, help="phasing: ignore k-mers carried by more than this many LTR-RTs (0 = choose from 10,25,50,100,250 by split-half reproducibility)")
    ap.add_argument("--nmax_rel", type=int, default=250, help="per-copy relatives: ignore k-mers carried by more than this many LTR-RTs")
    ap.add_argument("--min_elements", type=int, default=10, help="phase only sequences carrying at least this many LTR-RTs (skips small unplaced scaffolds)")
    ap.add_argument("--boot", type=int, default=100, help="half-sampling replicates for chromosome support (0 = off)")
    ap.add_argument("--min_post", type=float, default=0.9, help="posterior needed to call an element's lineage")
    ap.add_argument("--clock_boot", type=int, default=200, help="bootstrap replicates for the transposition-clock dates")
    ap.add_argument("--age_weight", type=int, default=1, choices=(0, 1), help="with --k2p, weight each copy's exchange evidence by how confined to one subgenome lineages of its age are (1) or count all copies equally (0)")
    ap.add_argument("--min_bf_exchange", type=float, default=100.0, help="likelihood ratio of a run's origin over the chromosome's own subgenome needed to report a candidate exchange (100 = 'decisive')")
    ap.add_argument("--seed", type=int, default=1, help="random seed")
    ap.add_argument("--no_plots", action="store_true", help="skip figures")
    ap.add_argument("-v", "--verbose", action="store_true", help="per-step progress and sanity checks")
    a = ap.parse_args(argv)
    VERBOSE = a.verbose
    if not 5 <= a.k <= 31:
        die("-k must be between 5 and 31")
    if not os.path.exists(a.ltr_fasta):
        die(f"not found: {a.ltr_fasta}")
    os.makedirs(a.outdir, exist_ok=True)
    rng = np.random.default_rng(a.seed)
    cfg = read_config(a.config) if a.config else None
    N = a.n_subgenomes or (max(len(s) for s in cfg) if cfg else None)
    if N is not None and N < 2:
        die("-n must be >= 2 (omit it to infer the number of subgenomes)")
    if a.scale <= 0:
        size = os.path.getsize(a.ltr_fasta) * (4 if a.ltr_fasta.endswith(".gz") else 1)
        a.scale = max(4, int(math.ceil(size / 75e6)))
    log(f"v{__version__}: {a.ltr_fasta} -> {a.outdir} (N={N or 'auto'}, k={a.k}, 1/{a.scale} of k-mers)")

    lib = Library(a.ltr_fasta, a.k, a.scale)
    cnt = {}
    for c in lib.chrom:
        cnt[c] = cnt.get(c, 0) + 1
    chroms = [c for c in lib.chroms if cnt[c] >= a.min_elements]
    if cfg:
        missing = sorted({c for s in cfg for c in s} - set(chroms))
        if missing:
            log(f"WARNING: {len(missing)} config chromosomes absent or below --min_elements: {' '.join(missing[:10])}")
    dropped = len(lib.chroms) - len(chroms)
    if len(chroms) < max(N or 2, 2):
        die(f"only {len(chroms)} sequences carry >= {a.min_elements} LTR-RTs; cannot form {max(N or 2, 2)} subgenomes")
    cix = {c: i for i, c in enumerate(chroms)}
    C = len(chroms)
    echrom = np.array([cix.get(c, -1) for c in lib.chrom], np.int32)
    keep = echrom >= 0
    log(f"read {lib.n} LTR-RTs on {len(lib.chroms)} sequences; phasing {C} with >= {a.min_elements} LTR-RTs"
        + (f" ({dropped} smaller sequences ignored)" if dropped else ""))
    # 1. number of subgenomes (unless given) and chromosome phasing
    ek0 = np.where(keep, echrom, 0)
    cands = [a.nmax] if a.nmax > 0 else [10, 25, 50, 100, 250]
    sp_all = Splits(lib.kmer, lib.eid, ek0, C, max(max(cands), a.nmax_rel), keep)
    if sp_all.n == 0:
        die("no shared k-mers link two chromosomes; is this an LTR-RT library?")
    if len(cands) > 1:
        a.nmax, wscores = choose_window(sp_all, ek0, C, lib.n, keep, cands, rng)
        log("copy-number window (split-half reproducibility r): "
            + ", ".join(f"<={k}: {v:.3f}" for k, v in wscores.items()) + f" -> {a.nmax}")
    else:
        wscores = {}
    a.window_scores = wscores
    sp = split_view(sp_all, sp_all.size <= a.nmax)
    sp_rel = split_view(sp_all, sp_all.size <= max(a.nmax_rel, a.nmax))
    init, tree = None, []
    if N is None:
        W0 = affinity(sp, ek0, C, keep.astype(float), "chrom", 0.0)
        pairs = half_pairs(sp, ek0, C, lib.n, keep, a.auto_reps, rng)
        lab, tests = find_subgenomes(W0, pairs, rng, min_rep=a.min_rep)
        for t in tests:
            pa = [chroms[i] for i in t["members"][t["split"] == 0]]
            pb = [chroms[i] for i in t["members"][t["split"] == 1]]
            tree.append(dict(part_a=pa, part_b=pb, replicability=t["rep"], contrast=t["contrast"],
                             parent_contrast=t.get("parent_contrast"), accepted=bool(t["accepted"])))
            log(f"split {len(pa)}|{len(pb)} chromosomes: replicability {t['rep']:.2f}, strength {t['contrast']:.2f}"
                + (f" (parent {t['parent_contrast']:.2f})" if t.get("parent_contrast") is not None else "")
                + (" -> accepted" if t["accepted"] else " -> rejected"))
        N = int(lab.max()) + 1
        log(f"inferred number of subgenomes: {N}")
        if N == 1:
            write_null(a, chroms, W0, lib, echrom, tree)
            return
        init = lab
    sets = make_sets(chroms, cfg, N)
    res = phase(lib, ek0, C, N, sets, a.nmax, "chrom", 0.0, a.boot, rng,
                emask=keep, boot_mode="half", sp=sp, init=init)
    sp, z, support = res["splits"], res["z"], res["support"]
    vlog(f"{sp.n_shared_kmers} shared k-mers -> {sp.n} distinct carrier sets spanning >= 2 sequences")
    sgn = [f"SG{g + 1}" for g in range(N)]
    for g in range(N):
        mem = [chroms[i] for i in range(C) if z[i] == g]
        log(f"{sgn[g]}: {len(mem)} sequences, {int(res['n_el'][z == g].sum())} LTR-RTs")
    weak = [chroms[i] for i in range(C) if np.isfinite(support[i]) and support[i] < 0.95]
    if weak:
        log(f"support < 0.95 for {len(weak)} sequences: {' '.join(weak[:12])}")

    # 2. element lineages
    echk = np.where(keep, echrom, 0)
    top_c, top_r, n_local = nearest_chromosomes(sp_rel, echk, C, lib.n, max(C - 1, 1), want_local=True)
    top_c[~keep] = -1
    lm = lineage_model(top_c, top_r, echk, z, N, res["n_el"])
    P, call, lin, bp = lineage_calls(lm, echk, z, N, a.min_post)
    h = np.stack([lm["W"][:, z == g].sum(1) for g in range(N)], 1)      # lineage-weighted relatives per SG
    m = h.sum(1)
    # no-confinement expectation for the share of other-subgenome lineage relatives: other chromosomes per SG
    nch = np.array([(z == g).sum() for g in range(N)], float)
    p_null = nch[None, :] - np.eye(N)[z[echk]]
    p_null = p_null / p_null.sum(1, keepdims=True)
    par = dict(pi_own=lm["pi_own"], pi_shared=lm["pi_shared"], pi_foreign=lm["pi_foreign"], p_null=p_null,
               relative_curve=dict(closeness=lm["a_closeness"], lineage_probability=lm["a"]))
    call[~keep] = "not_phased"
    lin[~keep] = -1
    own_sg = z[echk]
    share = np.where(m > 0, 1 - h[np.arange(lib.n), own_sg] / np.maximum(m, 1e-12), np.nan)
    share[~keep] = np.nan
    vals, nn = np.unique(call, return_counts=True)
    log("LTR-RT lineage calls: " + ", ".join(f"{v}={n}" for v, n in zip(vals, nn)))

    # 3. ages and dating
    k2p = np.full(lib.n, np.nan)
    k2se = np.full(lib.n, np.nan)
    ltrl = np.full(lib.n, np.nan)
    ksub = np.full(lib.n, np.nan)
    fit = None
    se = None
    s_null = np.nan
    if a.k2p:
        km, bad = read_k2p(a.k2p, a.k2p_col)
        for i, id_ in enumerate(lib.ids):
            v = km.get(id_)
            if v is not None:
                k2p[i], k2se[i], ltrl[i], ksub[i] = v
        nk = int(np.isfinite(k2p).sum())
        log(f"K2P for {nk}/{lib.n} LTR-RTs" + (f" ({bad} unreadable rows)" if bad else ""))
        if nk == 0:
            log("WARNING: no K2P ids match FASTA ids; ages skipped")
        se = np.where(np.isfinite(k2se) & (k2se > 0), k2se, k2p_se_guess(k2p, ltrl))
        fit = fit_clock(k2p, se, np.nan_to_num(share), np.where(keep, m, 0), rng, boot=a.clock_boot)
        s_null = float(np.nanmedian(1 - par["p_null"][np.arange(lib.n), own_sg][keep]))
        if fit:
            pooled = fit_merger_pooled(k2p, ksub, ltrl, share, np.where(keep, m, 0), s_null,
                                       fit.get("k2p_trough", np.inf), rng, boot=min(a.clock_boot, 50))
            fit["pooled"] = pooled
        if fit and fit.get("pooled"):
            pm = fit["pooled"]
            ky = lambda x: f"{x / (2 * a.mu) / 1e3:,.1f}"
            if pm["tau"] > 0:
                log(f"merger (pooled transposition clock, {pm['n']} young copies): K2P {pm['tau']:.3g} = {ky(pm['tau'])} kyr "
                    f"(95% CI {ky(pm['ci'][0])}-{ky(pm['ci'][1])} kyr, mu={a.mu:g})")
            else:
                log(f"merger below resolution: younger than {ky(pm['upper'])} kyr (95% upper bound, mu={a.mu:g})")
        else:
            log("too few dated LTR-RTs to fit the transposition clock")

    # 4. exchange segmentation
    order = {c: np.flatnonzero(echrom == c)[np.argsort(lib.start[echrom == c])] for c in range(C)}
    wseg = 1.0 / (1.0 + n_local)
    agew = age_informativeness(k2p, fit, s_null) if (a.age_weight and a.k2p) else None
    if agew is not None:
        wseg = wseg * agew
        log(f"exchange HMM: copies weighted by the subgenome confinement of their age class (median weight {np.nanmedian(agew):.2f})")
    state, post, em = segment_hmm(order, model_emissions(lm, N), z, N, weight=wseg)
    segs = segments_from_states(order, state, post, em, call, lin, lib.start, lib.end, z, chroms)
    for s in segs:
        s["flag"] = s["origin"] != s["chrom_sg"] and s["llr"] >= math.log(a.min_bf_exchange)
    nflag = sum(s["flag"] for s in segs)
    log(f"{nflag} candidate exchange segments (origin subgenome differs from the chromosome's)")
    gsegs, gstates, paint = [], {}, None
    if a.genome:
        markers = subgenome_markers(lib, echk, keep, z, N, max(a.nmax_rel, a.nmax))
        log("genome painting with " + ", ".join(f"{m.size} {sgn[j]} markers" for j, m in enumerate(markers)))
        paint = genome_paint(a.genome, chroms, markers, a.k, a.scale, a.win)
        if not paint:
            log("WARNING: no genome sequence names match the LTR-RT chromosome names; genome painting skipped")
        else:
            gsegs, gstates = genome_segments(paint, chroms, z, N, min_llr=math.log(a.min_bf_exchange))
            a.seq_assign = assign_sequences(paint, N)
            agree = [a.seq_assign[c][0] == z[i] for i, c in enumerate(chroms) if c in a.seq_assign]
            extra = [n for n in a.seq_assign if n not in set(chroms) and a.seq_assign[n][0] >= 0]
            log(f"genome markers agree with the LTR-RT assignment for {sum(agree)}/{len(agree)} phased sequences; "
                f"{len(extra)} further sequences assigned from genome markers (genome_sequences.tsv)")
            for sg_ in segs:
                sg_["genome_support"] = (None if sg_["chrom"] not in paint else
                                         any(g["flag"] and g["chrom"] == sg_["chrom"] and g["origin"] == sg_["origin"]
                                             and min(g["end"], sg_["end"]) > max(g["start"], sg_["start"]) for g in gsegs))
            log(f"{sum(g['flag'] for g in gsegs)} genome-based candidate exchanges; "
                f"{sum(1 for x in segs if x['flag'] and x.get('genome_support'))}/{nflag} LTR-based candidates supported")

    # 5. outputs
    clen = np.zeros(C)
    np.maximum.at(clen, echk[keep], lib.end[keep])
    if a.genome_fai:
        with open(a.genome_fai) as f:
            for line in f:
                p = line.split("\t")
                if p[0] in cix:
                    clen[cix[p[0]]] = float(p[1])
    a.split_tree = tree
    a.gsegs, a.paint = gsegs, paint
    write_outputs(a, lib, chroms, z, res, call, lin, bp, h, k2p, segs, fit, par, sgn, N, clen, keep, echrom, share)
    if not a.no_plots:
        try:
            plot_phasing(os.path.join(a.outdir, "fig_phasing"), chroms, z, support, res["W"], res["E"],
                         res["own_oe"], res["oth_oe"], N, sgn)
            plot_painting(os.path.join(a.outdir, "fig_painting"), chroms, z, support, clen, order, h, state,
                          lib.start, lib.end, segs, N, sgn, gstates=gstates, gsegs=gsegs, win=a.win)
            if a.k2p and np.isfinite(k2p).any():
                plot_ages(os.path.join(a.outdir, "fig_ages"), k2p, call, lin, a.mu, fit, share,
                          float(np.nanmedian(1 - par["p_null"][np.arange(lib.n), own_sg][keep])))
            write_legend(a, chroms, z, support, res, call, k2p, fit, segs, sgn, N, par)
        except Exception as ex:  # figures never block the tables
            log(f"WARNING: plotting failed ({ex}); tables are complete")
    log("done")


def write_null(a, chroms, W, lib, echrom, tree):
    """Outputs when no reproducible subgenome structure is found (N=1)."""
    od = a.outdir
    log("no reproducible subgenome structure: the library behaves like a diploid or an autopolyploid "
        "(use -n to force a partition); writing chromosomes.tsv, links_oe.tsv, summary.json")
    with open(os.path.join(od, "chromosomes.tsv"), "w") as f:
        f.write("chrom\tsubgenome\tsupport\tn_ltr\n")
        for i, c in enumerate(chroms):
            f.write(f"{c}\tSG1\tNA\t{int((echrom == i).sum())}\n")
    E = expected(W)
    OE = W / np.where(E > 0, E, np.nan)
    with open(os.path.join(od, "links_oe.tsv"), "w") as f:
        f.write("chrom\t" + "\t".join(chroms) + "\n")
        for i, c in enumerate(chroms):
            f.write(c + "\t" + "\t".join("NA" if not np.isfinite(x) else f"{x:.4f}" for x in OE[i]) + "\n")
    with open(os.path.join(od, "summary.json"), "w") as f:
        json.dump(dict(version=__version__, ltr_fasta=a.ltr_fasta, n_subgenomes=1, n_subgenomes_source="inferred",
                       split_tree=tree, n_ltr=lib.n, n_phased_sequences=len(chroms)), f, indent=2, default=float)
    log("done")


def write_outputs(a, lib, chroms, z, res, call, lin, bp, h, k2p, segs, fit, par, sgn, N, clen, keep, echrom, share):
    od = a.outdir
    C = len(chroms)
    with open(os.path.join(od, "chromosomes.tsv"), "w") as f:
        f.write("chrom\tsubgenome\tsupport\tpvalue\toe_own\toe_other\tn_ltr\tlength\tn_own\tn_foreign\tn_shared\tn_unresolved\tn_exchange_segments\n")
        for i, c in enumerate(chroms):
            mi = echrom == i
            cc = call[mi]
            nseg = sum(1 for s in segs if s["chrom"] == c and s["flag"])
            f.write(f"{c}\t{sgn[z[i]]}\t{res['support'][i]:.3f}\t{res['pval'][i]:.3g}\t{res['own_oe'][i]:.3f}\t"
                    f"{res['oth_oe'][i]:.3f}\t{int(mi.sum())}\t{int(clen[i])}\t{int((cc == 'own').sum())}\t"
                    f"{int((cc == 'foreign').sum())}\t{int((cc == 'shared').sum())}\t"
                    f"{int(np.isin(cc, ['unresolved', 'no_relatives']).sum())}\t{nseg}\n")
    cix = {c: i for i, c in enumerate(chroms)}
    with open(os.path.join(od, "elements.tsv"), "w") as f:
        f.write("id\tchrom\tstart\tend\tclass\tchrom_subgenome\tlineage_call\tlineage_subgenome\tposterior\t"
                + "\t".join(f"relatives_{s}" for s in sgn) + "\tother_sg_share\tk2p\tage_years\n")
        for i in range(lib.n):
            ci = cix.get(lib.chrom[i])
            csg = sgn[z[ci]] if ci is not None else "NA"
            ls = sgn[lin[i]] if lin[i] >= 0 else "NA"
            post = f"{bp[i]:.3f}" if call[i] in ("own", "foreign", "shared", "unresolved") else "NA"
            rel = "\t".join(f"{x:.2f}" for x in h[i]) if keep[i] else "\t".join("NA" for _ in sgn)
            kv = f"{k2p[i]:.5g}" if np.isfinite(k2p[i]) else "NA"
            age = f"{k2p[i] / (2 * a.mu):.0f}" if np.isfinite(k2p[i]) else "NA"
            sh = f"{share[i]:.3f}" if np.isfinite(share[i]) else "NA"
            f.write(f"{lib.ids[i]}\t{lib.chrom[i]}\t{lib.start[i]}\t{lib.end[i]}\t{lib.cls[i]}\t{csg}\t{call[i]}\t"
                    f"{ls}\t{post}\t{rel}\t{sh}\t{kv}\t{age}\n")
    with open(os.path.join(od, "links_oe.tsv"), "w") as f:
        OE = res["W"] / np.where(res["E"] > 0, res["E"], np.nan)
        f.write("chrom\t" + "\t".join(chroms) + "\n")
        for i, c in enumerate(chroms):
            f.write(c + "\t" + "\t".join("NA" if not np.isfinite(x) else f"{x:.4f}" for x in OE[i]) + "\n")
    with open(os.path.join(od, "segments.tsv"), "w") as f:
        f.write("chrom\tstart\tend\torigin_subgenome\tchrom_subgenome\tn_ltr\tn_confined_support\tn_confined_conflict\tllr\tmean_posterior\tcandidate_exchange\tgenome_support\n")
        for s in segs:
            gs = s.get("genome_support")
            f.write(f"{s['chrom']}\t{s['start']}\t{s['end']}\t{sgn[s['origin']]}\t{sgn[s['chrom_sg']]}\t{s['n_elements']}\t"
                    f"{s['n_support']}\t{s['n_conflict']}\t{s['llr']:.1f}\t{s['mean_post']:.3f}\t{'yes' if s['flag'] else 'no'}\t"
                    f"{'NA' if gs is None else ('yes' if gs else 'no')}\n")
    if getattr(a, "gsegs", None) is not None and a.genome:
        with open(os.path.join(od, "segments_genome.tsv"), "w") as f:
            f.write("chrom\tstart\tend\torigin_subgenome\tchrom_subgenome\tn_windows\tmarkers\tllr\tcandidate_exchange\n")
            for g in a.gsegs:
                f.write(f"{g['chrom']}\t{g['start']}\t{g['end']}\t{sgn[g['origin']]}\t{sgn[g['chrom_sg']]}\t{g['n_windows']}\t"
                        f"{g['markers']}\t{g['llr']:.1f}\t{'yes' if g['flag'] else 'no'}\n")
        with open(os.path.join(od, "genome_sequences.tsv"), "w") as f:
            f.write("sequence\tlength\tltr_subgenome\tgenome_subgenome\t" + "\t".join(f"windows_{x}" for x in sgn)
                    + "\tmajority_fraction\tpvalue\n")
            cix_ = {c: i for i, c in enumerate(chroms)}
            for name, (g, n, fr, p) in a.seq_assign.items():
                ln = int(a.paint[name][0][-1] + a.win)
                ls = sgn[z[cix_[name]]] if name in cix_ else "NA"
                f.write(f"{name}\t{ln}\t{ls}\t{sgn[g] if g >= 0 else 'NA'}\t" + "\t".join(str(int(x)) for x in n)
                        + f"\t{fr:.3f}\t{p:.3g}\n")
        with open(os.path.join(od, "genome_windows.tsv"), "w") as f:
            f.write("chrom\tstart\tend\t" + "\t".join(f"markers_{x}" for x in sgn) + "\n")
            for c, (st, cnt) in a.paint.items():
                for i in range(len(st)):
                    f.write(f"{c}\t{st[i]}\t{st[i] + a.win}\t" + "\t".join(str(int(x)) for x in cnt[i]) + "\n")
    summ = dict(version=__version__, ltr_fasta=a.ltr_fasta, n_subgenomes=N, nmax_phasing=a.nmax,
                nmax_relatives=max(a.nmax_rel, a.nmax), window_reproducibility=getattr(a, "window_scores", {}),
                n_subgenomes_source=("given" if a.n_subgenomes else "config" if a.config else "inferred"),
                split_tree=getattr(a, "split_tree", []), k=a.k, scale=a.scale, nmax=a.nmax,
                n_ltr=lib.n, n_phased_sequences=C, n_splits=int(res["splits"].n),
                mean_oe_within=float(np.nanmean(res["own_oe"])), mean_oe_between=float(np.nanmean(res["oth_oe"])),
                min_support=float(np.nanmin(res["support"])) if a.boot else None,
                lineage_mixture={k: float(v) for k, v in par.items() if k not in ("p_null", "loglik", "relative_curve", "rel_w")},
                relative_informativeness=par.get("relative_curve"),
                lineage_calls={str(v): int(n) for v, n in zip(*np.unique(call, return_counts=True))},
                n_candidate_exchanges=int(sum(s["flag"] for s in segs)), mu=a.mu,
                n_genome_candidate_exchanges=(int(sum(g["flag"] for g in a.gsegs)) if a.genome and a.gsegs is not None else None),
                transposition_clock=(None if not fit else {k: v for k, v in fit.items() if k != "bins"}))
    if fit:
        yrs = lambda x: None if x is None else x / (2 * a.mu)
        pm = fit.get("pooled")
        if pm:
            summ["merger_years"] = dict(estimate=yrs(pm["tau"]), ci95=[yrs(x) for x in pm["ci"]],
                                        upper95=yrs(pm["upper"]), method="pooled per-copy transposition clock")
        if "tau_merger" in fit:
            summ["merger_binned_years"] = dict(estimate=yrs(fit["tau_merger"]), ci95=[yrs(x) for x in fit["tau_merger_ci"]])
        if "tau_divergence" in fit:
            summ["progenitor_split_years_tentative"] = dict(estimate=yrs(fit["tau_divergence"]),
                                                            ci95=[yrs(x) for x in fit["tau_divergence_ci"]],
                                                            supported=fit.get("divergence_supported"))
    conf = (call == "own") | (call == "foreign")
    if np.isfinite(k2p).any():
        own_age = k2p[(call == "own") & np.isfinite(k2p)]
        if own_age.size >= 20:
            summ["subphaser_style_window_k2p"] = dict(own_p2_5=float(np.quantile(own_age, 0.025)),
                                                     own_median=float(np.median(own_age)),
                                                     own_p97_5=float(np.quantile(own_age, 0.975)))
    with open(os.path.join(od, "summary.json"), "w") as f:
        json.dump(summ, f, indent=2, default=float)


def write_legend(a, chroms, z, support, res, call, k2p, fit, segs, sgn, N, par):
    C = len(chroms)
    n_hi = int(np.sum(np.asarray(support) >= 0.95))
    calls = {v: int(n) for v, n in zip(*np.unique(call, return_counts=True))}
    lines = [
        "# Figure legends (v4 LTR-RT subgenome phasing)", "",
        "Shared terms. *Link*: two LTR retrotransposon (LTR-RT) copies on different chromosomes that share a "
        f"{a.k}-mer carried by at most {a.nmax} copies in the library (a rare, clade-defining sequence variant); copies "
        "on the same chromosome are never counted. *Observed/expected*: links between two chromosomes relative to a "
        "null in which each chromosome keeps its total number of links. *Support*: fraction of "
        f"{a.boot} random half-libraries (50% of LTR-RTs) that place a chromosome in the same subgenome. "
        "*Lineage call*: whether a copy's closest relatives on other chromosomes lie in its own subgenome (own), "
        "in another subgenome (foreign), or in both (shared); posterior >= "
        f"{a.min_post}.", "",
        f"**fig_phasing.** (A) Links between all {C} chromosomes, observed/expected (red: more links than expected, "
        f"blue: fewer), ordered by inferred subgenome (coloured bars). Chromosomes of one subgenome share more "
        f"rare LTR-RT variants with each other (mean observed/expected {np.nanmean(res['own_oe']):.2f}) than with the "
        f"other subgenome ({np.nanmean(res['oth_oe']):.2f}), because each progenitor amplified its own lineages before "
        f"the merger. (B) Each chromosome's mean observed/expected links to its own subgenome versus the closest other "
        f"subgenome; points above the dashed diagonal prefer their own subgenome. Filled points: support >= 0.95 "
        f"({n_hi}/{C} chromosomes).", "",
        f"**fig_painting.** Chromosomes to scale (Mb). Ticks: every LTR-RT with relatives on other chromosomes, at its "
        f"position, coloured by the subgenome that holds >= 75% of its closest relatives (grey: mixed). Coloured bar under "
        f"each chromosome: the origin inferred along it by a hidden Markov model. Black boxes: "
        f"runs of copies whose relatives point to another subgenome (hidden Markov model along the "
        f"chromosome; likelihood ratio >= {a.min_bf_exchange:g} over the chromosome's own subgenome; "
        f"n={sum(s['flag'] for s in segs)}), candidate homoeologous exchanges or "
        f"assembly switch errors. Numbers after chromosome names: support."
        + ("" if not a.genome else
           f" Lower bar (with --genome): origin along the genome from {a.win // 1000}-kb windows, each voting for the "
           f"subgenome whose marker k-mers (k-mers of lineages confined to one subgenome, also found in fragmented and "
           f"solo LTRs) it holds most; dashed boxes: runs with likelihood ratio >= {a.min_bf_exchange:g} "
           f"(n={sum(g['flag'] for g in (a.gsegs or []))})."), ""]
    if a.k2p and np.isfinite(k2p).any():
        ky = lambda x: f"{x / (2 * a.mu) / 1e3:,.0f}"
        if fit and fit.get("pooled"):
            pm = fit["pooled"]
            c = pm["ci"]
            t = (f"Fitting every young copy at the resolution of its own LTR length ({pm['n']:,} copies; post-merger "
                 f"copies expected at the no-confinement share {pm['share_post']:.2f}, older ones at "
                 f"{pm['share_confined']:.2f}) places the merger at K2P {pm['tau']:.3g} ({ky(pm['tau'])} kyr, 95% bootstrap "
                 f"CI {ky(c[0])}-{ky(c[1])} kyr; one-sided 95% upper bound {ky(pm['upper'])} kyr)")

        else:
            t = "too few dated copies to fit the clock"
        lines += [f"**fig_ages.** (A) LTR divergence (K2P between the two LTRs of a copy; insertion time on the top "
                  f"axis assumes {a.mu:g} substitutions/site/year) of copies in own-, shared- and foreign-lineage "
                  f"classes. (B) Cross-subgenome transposition clock: for copies binned by K2P (points; area = evidence), "
                  f"the share of each copy's closest relatives that sit on chromosomes of another subgenome (dashed: "
                  f"share expected without subgenome confinement). Before the merger a lineage that grew in one "
                  f"progenitor could not reach the other subgenome, so copies inserted between the progenitor split and "
                  f"the merger have low shares; after the merger, copies land in either subgenome. {t}. Shaded: 95% "
                  f"CIs (bootstrap over copies).", ""]
    with open(os.path.join(a.outdir, "FIGURE_LEGEND.md"), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
