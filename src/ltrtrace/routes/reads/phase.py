"""Read phasing at one element: how many distinct haplotypes carry each allele.

Reads anchored at an element are compared over its unique flanks. Only bases with quality >= min_q count (as the
reference or an alternative base; others are missing). Heterozygous sites are substitutions seen in >= 2 reads and
in a minor-allele fraction of the covering reads of at least min_af, kept only if linked to another such site by the
four-gamete test (real haplotype differences are linked; independent sequencing errors are not). Reads that agree
at the sites they share are joined into groups; two groups are distinct when some pair of their reads disagrees at
>= min_dis shared sites and at >= 70% of them. For each allele class (F, S, E) the result is the largest number of
pairwise distinct groups that each hold >= 2 reads of that class: a lower bound on the haplotypes carrying that
allele. Identical haplotypes (no heterozygous sites) give at most 1; the bound is uninformative in homozygous or
hemizygous regions. Coordinates are 0-based half-open.
"""
from __future__ import annotations

import zlib

import numpy as np

from ltrtrace.base.constants import (
    PH_DISTINCT_DIS,
    PH_FLANK,
    PH_GAMETE_FRAC,
    PH_MAX_READS,
    PH_MAX_SITES,
    PH_MIN_AF,
    PH_MIN_ALLELE_READS,
    PH_MIN_DIS,
    PH_MIN_LINK_READS,
    PH_MIN_Q,
    PH_MIN_READS,
    PH_SAME_DIS,
    PH_SAME_SITES,
)
from ltrtrace.base.unionfind import find
from ltrtrace.seq.cigar import blocks

_ACGT = np.zeros(256, dtype=bool)
for _c in b"ACGT":
    _ACGT[_c] = True


def _aligned_blocks(x) -> list[tuple[int, int, int]]:
    """(query start, reference start, length) of every aligned (M, =, X) block of one alignment (seq.cigar.blocks:
    query positions in the stored sequence, hard clips not counted)."""
    ct = x.cigartuples
    if not ct:
        return []
    ops, lens = zip(*ct)
    return blocks(ops, lens, 0, x.reference_start)


def _read_bases(alns, wins, width: int, min_q: int) -> np.ndarray:
    """One read over the concatenated flank windows: array of the read's base per window position (0 = not covered
    by an aligned base of quality >= min_q)."""
    row = np.zeros(width, dtype=np.uint8)
    for x in alns:
        seq = x.query_sequence
        if seq is None:
            continue
        qb = np.frombuffer(seq.encode(), dtype=np.uint8)
        qual = x.query_qualities
        qa = None
        if min_q > 0 and qual is not None and len(qual) == len(seq):
            qa = np.asarray(qual, dtype=np.int16)
        for q0, r0, ln in _aligned_blocks(x):
            for w0, w1, off in wins:
                lo, hi = max(r0, w0), min(r0 + ln, w1)
                if hi <= lo:
                    continue
                b = qb[q0 + lo - r0:q0 + hi - r0].copy()
                b[~_ACGT[b]] = 0
                if qa is not None:
                    b[qa[q0 + lo - r0:q0 + hi - r0] < min_q] = 0
                row[off + lo - w0:off + hi - w0] = b
    return row


def _compatible(col_i: np.ndarray, col_j: np.ndarray) -> bool:
    """Four-gamete test on two site columns (+1 alt, -1 ref, 0 missing): True if the sites are linked."""
    m = (col_i != 0) & (col_j != 0)
    if m.sum() < PH_MIN_LINK_READS:
        return False
    a, b = col_i[m] > 0, col_j[m] > 0
    cells = np.array([(a & b).sum(), (a & ~b).sum(), (~a & b).sum(), (~a & ~b).sum()])
    if min(a.sum(), (~a).sum(), b.sum(), (~b).sum()) < PH_MIN_ALLELE_READS:
        return False
    return cells.min() <= PH_GAMETE_FRAC * m.sum()


def haplotypes(reads, fasta, chrom: str, s: int, e: int, ploidy: int, min_q: int = PH_MIN_Q,
               min_af: float | None = None, min_dis: int = PH_MIN_DIS) -> tuple[dict[str, int], int]:
    """reads: list of (label, alignments). Returns ({F, S, E: distinct haplotypes, lower bound}, n_sites).
    min_af defaults to max(PH_MIN_AF, 0.5 / ploidy): half the fraction of one haplotype in `ploidy`."""
    out = {"F": 0, "S": 0, "E": 0}
    if ploidy < 2 or len(reads) < PH_MIN_READS:
        return out, 0
    if len(reads) > PH_MAX_READS:  # deterministic subset (read-name hash): bounds time and memory
        reads = sorted(reads, key=lambda r: zlib.crc32(r[1][0].query_name.encode()))[:PH_MAX_READS]
    if min_af is None:
        min_af = max(PH_MIN_AF, 0.5 / ploidy)
    L = fasta.length(chrom)
    wins, refs, width = [], [], 0
    for w0, w1 in ((max(0, s - PH_FLANK), s), (e, min(L, e + PH_FLANK))):
        if w1 > w0:
            wins.append((w0, w1, width))
            refs.append(np.frombuffer(fasta.fetch(chrom, w0, w1).upper().encode(), dtype=np.uint8))
            width += w1 - w0
    if not width:
        return out, 0
    ref = np.concatenate(refs)
    B = np.stack([_read_bases(alns, wins, width, min_q) for _, alns in reads])  # reads x positions
    labels = np.array([lab for lab, _ in reads])
    cov = B != 0
    ncov = cov.sum(0)
    alt = cov & (B != ref[None, :]) & _ACGT[ref][None, :]
    nalt = alt.sum(0)
    cand = np.nonzero((nalt >= PH_MIN_ALLELE_READS) & (nalt >= min_af * ncov) & (nalt <= (1 - min_af) * ncov))[0]
    if len(cand) < 2:
        return out, 0
    # the most frequent alternative base per candidate; reads with another base there are missing
    sites = []
    for p in cand:
        bases = B[alt[:, p], p]
        vals, cnt = np.unique(bases, return_counts=True)
        a = vals[np.argmax(cnt)]
        if cnt.max() >= PH_MIN_ALLELE_READS and cnt.max() >= min_af * ncov[p]:
            sites.append((p, a))
    if len(sites) < 2:
        return out, 0
    edge = s - max(0, s - PH_FLANK)  # offset of the element in the concatenated windows
    pos = np.array([p for p, _ in sites])
    order = np.argsort(np.abs(pos - edge), kind="stable")[:PH_MAX_SITES]  # sites nearest the element
    sites = [sites[k] for k in order]
    M = np.zeros((len(reads), len(sites)), dtype=np.int8)
    for j, (p, a) in enumerate(sites):
        col = B[:, p]
        M[col == a, j] = 1
        M[(col == ref[p]) & (col != 0), j] = -1
    keep = [j for j in range(len(sites)) if any(_compatible(M[:, j], M[:, k]) for k in range(len(sites)) if k != j)]
    if len(keep) < 2:
        return out, len(keep)
    M = M[:, keep].astype(np.int32)
    A, R = (M > 0).astype(np.int32), (M < 0).astype(np.int32)
    agree = A @ A.T + R @ R.T
    dis = A @ R.T + R @ A.T
    shared = agree + dis
    same = (shared >= PH_SAME_SITES) & (dis <= PH_SAME_DIS * shared)
    distinct = (dis >= min_dis) & (dis >= PH_DISTINCT_DIS * shared)
    n = len(reads)
    parent = list(range(n))
    for i, j in zip(*np.nonzero(np.triu(same, 1))):
        parent[find(parent, i)] = find(parent, j)
    comp = np.array([find(parent, i) for i in range(n)])
    for a in out:
        groups = [g for g in set(comp) if int(((comp == g) & (labels == a)).sum()) >= PH_MIN_ALLELE_READS]
        groups.sort(key=lambda g: -(comp == g).sum())
        chosen = []
        for g in groups:
            ig = np.nonzero(comp == g)[0]
            if all(distinct[np.ix_(ig, np.nonzero(comp == h)[0])].any() for h in chosen):
                chosen.append(g)
            if len(chosen) >= ploidy:
                break
        out[a] = len(chosen)
    return out, len(keep)
