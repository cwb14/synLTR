"""Resampling: jackknife blocks over sequence coordinates and the delete-one block jackknife; the lineage-cluster
bootstrap and the widening of a bootstrap interval by a model error (the subgenome dating intervals).

The jackknife is the unweighted delete-one form.
"""
from __future__ import annotations

import math
from collections.abc import Callable, Iterator, Sequence
from dataclasses import dataclass

import numpy as np

from ..base.constants import BLOCK_BP, BLOCK_FRAC, BLOCK_MIN, JK_MIN_BLOCKS, Z95


@dataclass(frozen=True)
class Jk:
    """A delete-one block jackknife: estimate, standard error, z = est / se, and the number of blocks."""
    est: float
    se: float
    z: float
    nblocks: int


def block_size(lengths: Sequence[int]) -> int:
    """Jackknife block length from the genome lengths (bp, one total per genome): the smaller of BLOCK_BP and
    1/BLOCK_FRAC of the smallest genome (at least BLOCK_MIN). The rule on no genomes raises a bare ValueError (min
    of nothing)."""
    smallest = min(int(n) for n in lengths)
    return min(BLOCK_BP, max(smallest // BLOCK_FRAC, BLOCK_MIN))


def block_ids(seq: np.ndarray, mid: np.ndarray, size: int) -> np.ndarray:
    """int64 block id (seq << 32) | (mid // size): a block never spans two sequences. seq >= 0. The ids order the
    blocks, and the jackknife sums its pseudovalues in that order, so the last bits of a standard error follow it: the
    library screen keeps its own ids (analyses.geneflow.library_stats._blocks: the same blocks, ranked by their
    names' text order, as LTRmix ranked them) so that its standard errors stay bit-identical to LTRmix's."""
    return (np.asarray(seq).astype(np.int64) << 32) | (np.asarray(mid) // size).astype(np.int64)


def jackknife(stat: Callable[[np.ndarray], float], blocks: np.ndarray) -> Jk:
    """Delete-one block jackknife of `stat`, a function of a boolean keep-mask over the rows.

    est = stat(all rows); a pseudovalue stat(blocks != b) per distinct block, non-finite ones dropped;
    se = sqrt((m - 1) / m * sum((p - mean)^2)) over the m finite pseudovalues, and exactly 0 when they are all
    equal (their float mean can be off by an ulp, which would give se ~1e-17); se and z are nan with fewer
    than JK_MIN_BLOCKS (3) distinct blocks (counted before dropping), with fewer than 2 finite
    pseudovalues, and at once, with no pseudovalue computed, when the estimate is nan; z =
    est / se when se > 0.
    """
    blocks = np.asarray(blocks)
    est = float(stat(np.ones(len(blocks), dtype=bool)))
    uniq = np.unique(blocks)
    nb = len(uniq)
    if math.isnan(est) or nb < JK_MIN_BLOCKS:
        return Jk(est, math.nan, math.nan, nb)
    p = np.array([stat(blocks != b) for b in uniq], dtype=float)
    p = p[np.isfinite(p)]
    m = len(p)
    if m < 2:
        return Jk(est, math.nan, math.nan, nb)
    # identical pseudovalues: no variation among the blocks (equal ratios of counts are the same double, so
    # the test is exact and needs no tolerance)
    se = 0.0 if p.min() == p.max() else math.sqrt((m - 1) / m * float(((p - p.mean()) ** 2).sum()))
    z = est / se if se > 0 else math.nan
    return Jk(est, se, z, nb)


def d_ratio(a: np.ndarray, b: np.ndarray) -> Callable[[np.ndarray], float]:
    """keep -> (sum a - sum b) / (sum a + sum b) over the kept rows; nan when the denominator is 0."""
    a = np.asarray(a)
    b = np.asarray(b)

    def stat(keep: np.ndarray) -> float:
        sa, sb = float(a[keep].sum()), float(b[keep].sum())
        return (sa - sb) / (sa + sb) if sa + sb else math.nan
    return stat


def cluster_bootstrap(clusters: np.ndarray, n_boot: int, rng: np.random.Generator) -> Iterator[np.ndarray]:
    """n_boot cluster-bootstrap replicates of the rows labeled by `clusters`: each draws as many clusters as there are
    distinct labels, uniformly with replacement (one rng.integers call, made when the replicate is taken), and yields
    the drawn clusters' row indices concatenated in draw order (each cluster's rows ascending; a cluster drawn twice
    appears twice). The caller computes its statistic on each index array."""
    u, inv = np.unique(clusters, return_inverse=True)
    idx = [np.flatnonzero(inv == k) for k in range(u.size)]
    for _ in range(n_boot):
        yield np.concatenate([idx[k] for k in rng.integers(0, u.size, u.size)])


def widen(ci: tuple[float, float, float], model_err: float) -> tuple[float, float]:
    """The bootstrap interval ci = (est, lo, hi) widened by a model error, a log-scale sd: each side's log-ratio to
    the estimate, log(est / lo) and log(hi / est), combined in quadrature with 1.96 x model_err -> (lo, hi). A lower
    bound lo <= 0 gives 0.0, and an estimate of 0 gives (0.0, 0.0)."""
    est, lo, hi = ci
    me = Z95 * model_err
    zl = np.log(est / lo) if lo > 0 else np.inf
    zu = np.log(hi / est) if est > 0 else 0.0
    return (float(est * np.exp(-np.hypot(zl, me))) if np.isfinite(zl) else 0.0), float(est * np.exp(np.hypot(zu, me)))
