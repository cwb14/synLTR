"""Collinear chains of anchors in rank space (seq.chain): the chain DP, the permutation and exact-LIS nulls, the chain
test and anchor ranks; LTRmix's greedy blocks and shuffled null."""
import collections
import itertools
from bisect import bisect_left
from pathlib import Path

import numpy as np
import pytest

from synltr.base import constants as cst
from synltr.seq import chain

GOLDEN = Path(__file__).resolve().parent.parent / "golden"


def test_longest_chain_planted():
    rng = np.random.default_rng(3)
    x = np.arange(200)
    y = rng.permutation(200)
    y[50:90] = np.arange(500, 540)                                  # 40 collinear points
    o = np.argsort(x)
    ch = chain.longest_chain(x[o], y[o])
    assert len(ch) >= 40 and set(range(50, 90)) <= set(ch.tolist())
    ch_rev = chain.longest_chain(x[o], -y[o])                          # inverted segment
    assert len(ch_rev) >= 40


def test_chain_significance_exact_for_few_candidates():
    """Few candidates: the normal approximation is too lenient (a chain of 8 among 10 candidates has p ~ 2e-3);
    the exact bound 2 C(n,l)/l! decides. Many candidates: the z score against permutations decides."""
    assert not chain.chain_significant(10, 8, 5.0, 1225)
    assert chain.chain_significant(10, 10, 9.0, 1225)
    assert chain.chain_significant(200, 40, 6.0, 1225)
    assert not chain.chain_significant(200, 40, 3.0, 1225)


def test_chain_null_seeded():
    """The permutation null draws only from the generator it is given, in phasing's order (one permutation of y per
    replicate): the same seed gives the same null, another seed another one, and the global numpy state is untouched."""
    rng = np.random.default_rng(1)
    x = np.arange(80, dtype=np.int64)
    y = rng.permutation(80).astype(np.int64)
    np.random.seed(5)
    state = np.random.get_state()[1].copy()
    a = chain.chain_null(x, y, np.random.default_rng(7))
    b = chain.chain_null(x, y, np.random.default_rng(7))
    c = chain.chain_null(x, y, np.random.default_rng(8))
    assert np.array_equal(np.random.get_state()[1], state)
    assert a.dtype == float and a.size == cst.NULL_PERMS and np.array_equal(a, b) and not np.array_equal(a, c)
    g = np.random.default_rng(7)
    assert a.tolist() == [chain.chain_length(x, g.permutation(y), cst.CHAIN_W) for _ in range(cst.NULL_PERMS)]
    assert chain.chain_null(x, y, np.random.default_rng(7), n=5, W=3).size == 5
    z = chain.chain_z(30, a)
    assert z == pytest.approx((30 - a.mean()) / max(a.std(), 0.5))
    assert chain.chain_z(3, np.full(10, 3.0)) == 0.0                             # a flat null: spread floored
    p1 = chain.chain_test(x, y, 12, 0.01, np.random.default_rng(3))
    assert p1 == chain.chain_test(x, y, 12, 0.01, np.random.default_rng(3))


def _lis(p):
    t = []
    for x in p:
        i = bisect_left(t, x)
        t[i:i + 1] = [x]
    return len(t)


def test_lis_tail_is_exact():
    for n in range(1, 8):
        lens = [_lis(p) for p in itertools.permutations(range(n))]
        for l in range(n + 2):
            assert chain.lis_tail(n)[l] == pytest.approx(np.mean([x >= l for x in lens]), abs=1e-12)


def test_chain_test_holds_its_level_on_random_anchors():
    """Non-homoeologous pairs = anchors in random order: p <= alpha in at most ~alpha of them (exact bound for
    n <= CHAIN_W, permutations above), while a collinear set passes."""
    rng = np.random.default_rng(3)
    for n, alpha in ((12, 0.05), (60, 0.05)):
        hits = 0
        for _ in range(150):
            x = np.arange(n, dtype=np.int64)
            y = rng.permutation(n).astype(np.int64)
            L = chain.chain_length(x, y, cst.CHAIN_W)
            hits += chain.chain_test(x, y, L, alpha, rng) <= alpha
        assert hits / 150 <= alpha + 0.04
        x = np.arange(n, dtype=np.int64)
        y = x.copy()
        y[: n // 4] = rng.permutation(y[: n // 4])                    # mostly collinear
        assert chain.chain_test(x, y, chain.chain_length(x, y, cst.CHAIN_W), alpha, rng) <= alpha


def test_anchor_ranks():
    """Anchors ordered by position on the first sequence (stable on ties): x = 0..n-1, y = each anchor's rank along
    the second (stable)."""
    pa = np.array([500, 100, 300, 100, 900])
    pb = np.array([7, 9, 2, 9, 1])
    o, x, y = chain.anchor_ranks(pa, pb)
    assert o.tolist() == [1, 3, 2, 0, 4] and x.tolist() == [0, 1, 2, 3, 4]
    assert y.tolist() == [3, 4, 1, 2, 0] and x.dtype == y.dtype == np.int64
    o, x, y = chain.anchor_ranks(np.zeros(0, np.int64), np.zeros(0, np.int64))
    assert o.size == x.size == y.size == 0


def _frozen_anchors(z):
    """The frozen anchor sets: (positions along one sequence, partner positions), in input order."""
    ends = np.cumsum(z["ch_sizes"])
    return [(px.astype(np.int64), py.astype(np.int64))
            for px, py in zip(np.split(z["ch_px"], ends[:-1]), np.split(z["ch_py"], ends[:-1]))]


def test_chain_blocks_match_ltrmix():
    """Greedy collinear blocks: the longest chain (chain_dp, steps <= W ranks) among the anchors left, either
    orientation (+1 on a tie), taken out while it holds >= minlen anchors; each block = (indices in position order,
    sign). LTRmix's blocks exactly on 30 frozen anchor sets (ties in both positions, planted runs both ways) at minlen
    5 and 8, W = CHAIN_GAP."""
    with np.load(GOLDEN / "ltrmix_route.npz") as z:
        cases, rows = _frozen_anchors(z), z["ch_blocks"].tolist()
    want = collections.defaultdict(list)
    for c, minlen, b, sign, i in rows:
        want[c, minlen, b, sign].append(i)
    got = {}
    for c, (px, py) in enumerate(cases):
        for minlen in (5, 8):
            for b, (ch, sign) in enumerate(chain.chain_blocks(px, py, minlen, W=cst.CHAIN_GAP)):
                assert ch.dtype == np.int64
                got[c, minlen, b, sign] = ch.tolist()
    assert got == dict(want) and any(k[3] == -1 for k in got)
    assert chain.chain_blocks(np.arange(3), np.arange(3), 5, W=20) == []


def test_shuffled_chain_matches_ltrmix():
    """The longest chain among anchors whose partner ranks are shuffled, the maximum over n draws (each
    rng.permutation(n) assigned in input order, read in position order): LTRmix's null length for each frozen
    anchor set with default_rng(case)."""
    with np.load(GOLDEN / "ltrmix_route.npz") as z:
        cases, want = _frozen_anchors(z), z["ch_null"].tolist()
    got = [chain.shuffled_chain(px, np.random.default_rng(c), n=cst.NULL_SHUFFLES, W=cst.CHAIN_GAP)
           for c, (px, _) in enumerate(cases)]
    assert got == want and all(type(v) is int for v in got)
