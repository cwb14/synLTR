"""Partitions maximizing within-group excess links, and the half-library side support."""
import itertools

import numpy as np
import pytest

from ltrtrace.analyses.subgenomes import links, partition


def _brute_best(B, sets):
    C, best, bz = B.shape[0], -np.inf, None
    for bits in itertools.product([1, -1], repeat=len(sets) - 1):
        s = np.zeros(C)
        for v, sg in enumerate((1,) + bits):
            s[sets[v][0]] = sg
            if len(sets[v]) == 2:
                s[sets[v][1]] = -sg
        if abs(s.sum()) == C:
            continue
        q = s @ B @ s
        if q > best + 1e-9:
            best, bz = q, s.copy()
    return best


@pytest.mark.parametrize("seed", range(5))
def test_gray_code_equals_brute_force(seed):
    rng = np.random.default_rng(seed)
    C = 9
    A = rng.normal(size=(C, C))
    B = A + A.T
    np.fill_diagonal(B, 0)
    sets = [[0], [1, 2], [3], [4], [5, 6], [7], [8]] if seed % 2 else [[i] for i in range(C)]
    z = partition.best_partition(B, 2, sets, rng)
    s = np.where(z == 0, 1.0, -1.0)
    assert abs(s @ B @ s - _brute_best(B, sets)) < 1e-6
    for st in sets:
        if len(st) == 2:
            assert z[st[0]] != z[st[1]]


def test_planted_partition_recovered():
    rng = np.random.default_rng(1)
    C, z0 = 12, np.repeat([0, 1, 2], 4)
    W = rng.random((C, C)) * 0.2 + (z0[:, None] == z0[None, :]) * 2.0
    W = W + W.T
    np.fill_diagonal(W, 0)
    z = partition.best_partition(W - links.expected(W), 3, [[i] for i in range(C)], rng)
    assert partition.ari(z, z0) == 1.0


def test_side_support_aligns_sides():
    z = np.array([0, 0, 1, 1])
    hz = [(np.array([1, 1, 0, 0]), np.array([0, 0, 1, 1])), (np.array([0, 1, 1, 1]), np.array([0, 0, 1, 1]))]
    assert np.allclose(partition.side_support(z, hz), [1.0, 0.75, 1.0, 1.0])
