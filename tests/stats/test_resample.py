import math

import numpy as np
import pytest

from synltr.stats.resample import block_ids, block_size, cluster_bootstrap, d_ratio, jackknife, widen


def test_jackknife_matches_leave_one_out():
    rng = np.random.default_rng(0)
    a = rng.integers(0, 2, 200); b = 1 - a; blocks = rng.integers(0, 20, 200)
    jk = jackknife(d_ratio(a, b), blocks)
    loo = np.array([(a[blocks != x].sum() - b[blocks != x].sum()) / (blocks != x).sum() for x in range(20)])
    assert jk.est == pytest.approx((a.sum() - b.sum()) / 200) and jk.nblocks == 20
    assert jk.se == pytest.approx(math.sqrt(19 / 20 * ((loo - loo.mean()) ** 2).sum()))


def test_jackknife_two_blocks_has_no_se():
    jk = jackknife(d_ratio(np.array([1, 1, 0]), np.array([0, 0, 1])), np.array([0, 0, 1]))
    assert jk.est == pytest.approx(1 / 3) and jk.nblocks == 2 and math.isnan(jk.se) and math.isnan(jk.z)


def test_block_size():
    assert block_size([1_000_000, 2_000_000]) == 10_000
    assert block_size([2_100_000_000]) == 5_000_000 and block_size([300_000_000]) == 3_000_000
    assert block_size([500_000]) == 10_000


def test_block_ids_split_by_seqid_and_size():
    ids = block_ids(np.array([0, 1, 0, 0, 0, 0]), np.array([5e3, 5e3, 19_999.5, 20_000.5, 21e3, 29e3]), 10_000)
    assert ids[0] != ids[1] and ids[2] != ids[3] and ids[3] == ids[4] == ids[5]


def test_jackknife_identical_pseudovalues_have_no_se():
    # three blocks of 3 a and 2 b: every pseudovalue is 0.2, whose float mean is 0.20000000000000004
    a = np.array([1, 1, 1, 0, 0] * 3)
    jk = jackknife(d_ratio(a, 1 - a), np.repeat([0, 1, 2], 5))
    assert jk.est == pytest.approx(0.2) and jk.se == 0.0 and math.isnan(jk.z) and jk.nblocks == 3


def test_jackknife_nan_estimate_short_circuits():
    """A nan estimate returns nan se and z at once, with no pseudovalue computed."""
    calls = []

    def stat(keep):
        calls.append(int(keep.sum()))
        return math.nan
    jk = jackknife(stat, np.repeat([0, 1, 2, 3], 2))
    assert math.isnan(jk.est) and math.isnan(jk.se) and math.isnan(jk.z) and jk.nblocks == 4 and calls == [8]


def test_jackknife_fewer_than_two_finite_pseudovalues():
    """With three blocks but one finite pseudovalue, se and z are nan; the estimate stands."""
    def stat(keep):     # finite with every block, or with block 0 deleted; nan otherwise
        return 3.0 if keep.all() else 1.0 if not keep[0] else math.nan
    jk = jackknife(stat, np.array([0, 1, 2]))
    assert jk.est == 3.0 and math.isnan(jk.se) and math.isnan(jk.z) and jk.nblocks == 3


def test_jackknife_drops_non_finite_pseudovalues():
    """A non-finite pseudovalue is dropped, the se taken over the m finite ones with (m - 1) / m."""
    v = np.array([1.0, 2.0, 4.0, 8.0])

    def stat(keep):     # the sum over the kept blocks; inf when block 3 is deleted
        return float(v[keep].sum()) if keep[3] else math.inf
    jk = jackknife(stat, np.arange(4))
    p = np.array([14.0, 13.0, 11.0])                          # the three finite pseudovalues
    se = math.sqrt(2 / 3 * ((p - p.mean()) ** 2).sum())
    assert jk.est == 15.0 and jk.nblocks == 4 and jk.se == pytest.approx(se) and jk.z == pytest.approx(15.0 / se)


def test_cluster_bootstrap_seeded():
    """Each replicate draws as many clusters as there are labels, with replacement, and yields their rows in draw
    order (one integers call per replicate): on this fixture and seed, these four replicates."""
    clusters = np.array([7, 3, 7, 9, 3, 3, 5])                  # clusters 3: rows 1, 4, 5; 5: 6; 7: 0, 2; 9: 3
    rng = np.random.default_rng(11)
    gen = cluster_bootstrap(clusters, 4, rng)
    assert iter(gen) is gen                                       # lazy: a replicate is drawn when it is taken
    assert [d.tolist() for d in gen] == [[1, 4, 5, 1, 4, 5, 3, 6], [0, 2, 0, 2, 0, 2, 1, 4, 5], [6, 1, 4, 5, 6, 3],
                                        [0, 2, 1, 4, 5, 0, 2, 1, 4, 5]]
    assert rng.integers(0, 1000) == 754                           # the stream continues right after the last draw


def test_widen_known_values():
    """1.96 x model_err in quadrature with each side's log-ratio to the estimate; a lower bound <= 0 gives 0.0."""
    # lo: 0.1 exp(-hypot(ln 1.25, 0.196)) = 0.1 exp(-0.2970001) = 0.07430440
    # hi: 0.1 exp(+hypot(ln 2, 0.196)) = 0.1 exp(0.7203256) = 0.2055102
    lo, hi = widen((0.1, 0.08, 0.2), 0.1)
    assert type(lo) is float and type(hi) is float
    assert lo == pytest.approx(0.07430440, rel=1e-6) and hi == pytest.approx(0.2055102, rel=1e-6)
    assert widen((0.1, 0.05, 0.2), 0.0) == pytest.approx((0.05, 0.2), rel=1e-12)    # no model error: unchanged
    assert widen((0.1, 0.1, 0.1), 0.1) == pytest.approx((0.1 / 1.2165269, 0.12165269), rel=1e-6)   # exp(0.196)
    for lo_ in (0.0, -0.01):                                      # lower bound at or below 0: 0.0, upper as before
        assert widen((0.1, lo_, 0.2), 0.1) == (0.0, pytest.approx(0.2055102, rel=1e-6))
    assert widen((0.0, 0.0, 0.0), 0.24) == (0.0, 0.0)             # estimate 0
    with np.errstate(divide="ignore"):                            # log(0 / lo) is -inf: not finite, so 0.0
        assert widen((0.0, 0.01, 0.02), 0.24) == (0.0, 0.0)
