"""Links between sequences: each clade adds 1/(m-1) to every pair of the m sequences it spans (a clade counts once;
same-sequence copies never link), plus nearest-relative votes; compared with a degree-corrected null."""
from __future__ import annotations

import math

import numpy as np
from numba import njit, prange
from scipy import sparse

from ltrtrace.analyses.subgenomes.inputs import warn
from ltrtrace.base.constants import EXPECTED_ITERS, K_NEIGH, NMAX_KNN, REL_BLOCK_NNZ, REL_KEEP_CELLS, REL_KEEP_MAX
from ltrtrace.base.runtime import n_chunks
from ltrtrace.seq.kmers import Clades
from ltrtrace.seq.relatives import best_by_group


@njit(parallel=True, cache=True)
def _link_kernel(members, indptr, code, C, nchunk):
    """Per chunk of clades (summed by the caller): each clade adds 1/(m-1) to every pair of the m distinct
    sequences spanned by its members with code >= 0 (code = sequence index of copies in the mask, else -1)."""
    n = indptr.size - 1
    out = np.zeros((nchunk, C, C))
    for ch in prange(nchunk):
        seen = np.full(C, -1, np.int64)
        buf = np.empty(C, np.int64)
        for k in range(n * ch // nchunk, n * (ch + 1) // nchunk):
            m = 0
            for p in range(indptr[k], indptr[k + 1]):
                sq = code[members[p]]
                if sq >= 0 and seen[sq] != k:
                    seen[sq] = k
                    buf[m] = sq
                    m += 1
            if m >= 2:
                w = 1.0 / (m - 1)
                for i in range(m):
                    for j in range(i + 1, m):
                        out[ch, buf[i], buf[j]] += w
                        out[ch, buf[j], buf[i]] += w
    return out


def link_matrix(cl, ek, C, mask):
    """C x C link matrix from the clades' members in mask: each clade adds 1/(m-1) to every pair of the m sequences
    it spans (a clade counts once however many copies or k-mers it has); same-sequence copies never link; copies
    with ek < 0 are ignored."""
    if cl.n == 0:
        return np.zeros((C, C))
    code = np.where(np.asarray(mask, bool), np.asarray(ek, np.int64), -1)
    nchunk = n_chunks(C, cl.n)
    return _link_kernel(cl.members, cl.indptr, code, C, nchunk).sum(0)


_EXPECTED_WARNED = False


def expected(W, iters=EXPECTED_ITERS):
    """Degree-corrected null for a zero-diagonal matrix: E_ij = b_i b_j with sum_{j!=i} E_ij = row sum. Warns
    (once per run) when the fixed-point iteration has not converged."""
    global _EXPECTED_WARNED
    d = W.sum(1)
    b = d / math.sqrt(max(d.sum(), 1e-300))
    converged = False
    for _ in range(iters):
        nb = np.sqrt(b * d / np.maximum(b.sum() - b, 1e-300))
        if np.allclose(nb, b, rtol=1e-12, atol=0):
            b, converged = nb, True
            break
        b = nb
    E = np.outer(b, b)
    np.fill_diagonal(E, 0)
    if not converged and not _EXPECTED_WARNED:
        _EXPECTED_WARNED = True
        err = float(np.max(np.abs(E.sum(1) - d)) / max(float(d.max()), 1e-300))
        warn(f"the degree-corrected null did not converge in {iters} iterations (largest row-sum error "
             f"{err:.2g} of the largest row); the last iterate is used")
    return E


def log_oe(W):
    return np.log((W + 1) / (expected(W) + 1))


def top_relatives(cl, ek, C, n_el, T, block_nnz=REL_BLOCK_NNZ, threads=1):
    """For every copy and every other sequence: its T closest relatives there (copy, relatedness; float32),
    relatedness = shared clades weighted by specificity (sum of mult/(size-1)); -1 = none; equal relatedness in the
    order of the sparse product (seq.relatives.best_by_group, tie "product"). Blocks of copies run in a fork pool."""
    A = sparse.csr_matrix((np.ones(cl.members.size, np.float32), cl.members, cl.indptr), shape=(cl.n, n_el))
    d = (cl.mult / np.maximum(cl.size - 1, 1)).astype(np.float32)
    DA = sparse.csr_matrix((np.repeat(d, cl.size), cl.members, cl.indptr), shape=(cl.n, n_el))
    return best_by_group(A, DA, ek, top=T, exclude_same=True, tie="product", n_groups=C, block_nnz=block_nnz,
                         threads=threads)


@njit(parallel=True, cache=True)
def _knn_kernel(rel_j, rel_r, ek, mask, C, k, nchunk):
    """Per chunk of copies: each copy in mask (phased) adds 1 to (its sequence, each of the k sequences holding its
    closest relatives within mask); ties go to the lower sequence index."""
    n, _, T = rel_j.shape
    out = np.zeros((nchunk, C, C))
    for ch in prange(nchunk):
        best = np.zeros(C)
        for e in range(n * ch // nchunk, n * (ch + 1) // nchunk):
            if not mask[e] or ek[e] < 0:
                continue
            for c in range(C):
                b = 0.0
                for t in range(T):
                    f = rel_j[e, c, t]
                    if f >= 0 and mask[f] and rel_r[e, c, t] > b:
                        b = rel_r[e, c, t]
                best[c] = b
            for _ in range(k):
                m, arg = 0.0, -1
                for c in range(C):
                    if best[c] > m:
                        m, arg = best[c], c
                if arg < 0:
                    break
                out[ch, ek[e], arg] += 1.0
                best[arg] = 0.0
    return out


def knn_votes(rel_j, rel_r, ek, C, mask):
    """C x C symmetric nearest-relative votes: every copy in mask adds 1 to its sequence x each of the K_NEIGH
    other sequences holding its closest relatives within mask (each copy votes once)."""
    k = min(K_NEIGH, C - 1)
    nchunk = n_chunks(C, rel_j.shape[0])
    W = _knn_kernel(rel_j, rel_r, np.asarray(ek, np.int64), np.asarray(mask, np.bool_), C, k, nchunk).sum(0)
    W = W + W.T
    np.fill_diagonal(W, 0)
    return W


def combine(L, Kv):
    if L.sum() > 0 and Kv.sum() > 0:
        return L + Kv * (L.sum() / Kv.sum())
    return L if L.sum() > 0 else Kv


class LinkModel:
    """Clades restricted to phased copies, and each copy's closest relatives (clades <= NMAX_KNN)."""

    def __init__(self, cl, ek, C, threads=1):
        self.ek, self.C = ek, C
        self.n = ek.size
        phased = ek >= 0
        self.cl = cl.restrict(phased, ek.clip(min=0))
        small = Clades(*_subset(self.cl, self.cl.size <= NMAX_KNN))
        T = int(max(1, min(REL_KEEP_MAX, REL_KEEP_CELLS // max(self.n * C, 1))))
        self.rel_j, self.rel_r = top_relatives(small, ek, C, self.n, T, threads=threads)

    def links(self, mask):
        return link_matrix(self.cl, self.ek, self.C, mask)

    def matrix(self, mask):
        return combine(self.links(mask), knn_votes(self.rel_j, self.rel_r, self.ek, self.C, mask))


def _subset(cl, keep):
    km = np.repeat(keep, cl.size)
    return cl.members[km], cl.size[keep], cl.mult[keep], cl.nseq[keep], cl.n_occ, cl.n_multi
