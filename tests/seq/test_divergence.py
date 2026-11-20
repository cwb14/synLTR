"""Divergence measures (seq.divergence): the Garwood upper bound, the JC69 distance of an identity, K2P, the standard
error of a copy's LTR K2P and the self-LTR check."""
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.stats import chi2

from ltrtrace.base import constants as cst
from ltrtrace.seq.divergence import (
    clean_tolerance,
    garwood_u95,
    jc69_from_identity,
    k2p,
    k2p_se_copy,
    k2p_shares,
    self_ltr_z,
)


def test_u95():
    assert garwood_u95(np.array([0]), np.array([1000]))[0] == pytest.approx(0.0029957, rel=1e-4)
    assert garwood_u95(np.array([2]), np.array([500]))[0] == pytest.approx(chi2.ppf(0.95, 6) / 2 / 500)


def test_u95_elementwise():
    """garwood_u95 with an array of trials is elementwise."""
    got = garwood_u95(np.array([0, 2, 5]), np.array([1000, 500, 250]))
    want = [chi2.ppf(0.95, 2 * (d + 1)) / 2 / n for d, n in ((0, 1000), (2, 500), (5, 250))]
    assert got.shape == (3,) and np.allclose(got, want)


def test_u95_of_scalars():
    assert float(garwood_u95(0, 1000)) == pytest.approx(0.0029957, rel=1e-4)


def test_jc69_from_identity():
    """-3/4 ln(1 - 4p/3) with p = 1 - identity; None where undefined (identity above 1, p >= 3/4); nan stays nan."""
    assert jc69_from_identity(0.9) == pytest.approx(0.107326, abs=1e-6)
    assert jc69_from_identity(1.0) == 0.0 and math.copysign(1.0, jc69_from_identity(1.0)) == 1.0  # never -0.0
    assert jc69_from_identity(1.01) is None and jc69_from_identity(0.25) is None
    assert jc69_from_identity(0.2501) == pytest.approx(-0.75 * math.log(1 - 4 * 0.7499 / 3))
    assert math.isnan(jc69_from_identity(math.nan))


def test_k2p_distance_and_se():
    """K2P = -1/2 ln(1 - 2P - Q) - 1/4 ln(1 - 2Q) for transition and transversion fractions P and Q, with the
    delta-method standard error."""
    d, se = k2p([100], [5], [2])
    P, Q = 0.05, 0.02
    assert d[0] == pytest.approx(-0.5 * math.log(1 - 2 * P - Q) - 0.25 * math.log(1 - 2 * Q), rel=1e-12)
    c1, c2 = 1 / (1 - 2 * P - Q), 0.5 * (1 / (1 - 2 * P - Q) + 1 / (1 - 2 * Q))
    var = (c1 * c1 * P + c2 * c2 * Q - (c1 * P + c2 * Q) ** 2) / 100
    assert se[0] == pytest.approx(math.sqrt(var), rel=1e-12)


def test_k2p_elementwise_and_undefined_is_nan():
    """Saturated (1 - 2P - Q <= 0 or 1 - 2Q <= 0) and empty (no sites) comparisons are NaN, d and se alike."""
    d, se = k2p([100, 100, 100, 0, 50], [5, 40, 0, 0, 0], [2, 20, 60, 0, 3])
    assert d.shape == se.shape == (5,)
    assert d[0] > 0 and se[0] > 0
    assert np.isnan(d[1:4]).all() and np.isnan(se[1:4]).all()
    assert d[4] == pytest.approx(-0.5 * math.log(1 - 0.06) - 0.25 * math.log(1 - 0.12))


def test_k2p_no_differences_is_zero():
    d, se = k2p([500], [0], [0])
    assert d[0] == 0.0 and se[0] == 0.0


def test_k2p_product_form_is_never_negative_zero():
    """The product form -1/2 ln((1 - 2P - Q) sqrt(1 - 2Q)) + 0.0 prints 0.0, never -0.0, for identical sequences;
    away from zero it equals the sum form to rounding."""
    sites, ts, tv = np.array([500, 400, 1000]), np.array([0, 3, 12]), np.array([0, 1, 5])
    d, _ = k2p(sites, ts, tv, form="product")
    assert d[0] == 0.0 and math.copysign(1.0, d[0]) == 1.0
    assert np.allclose(d, k2p(sites, ts, tv)[0], rtol=1e-12, atol=0)
    assert math.copysign(1.0, k2p([500], [0], [0])[0][0]) == -1.0   # the sum form of no difference is -0.0
    p = 3 / 400
    q = 1 / 400
    assert d[1] == -0.5 * math.log((1 - 2 * p - q) * math.sqrt(1 - 2 * q))


def test_k2p_values():
    """Known counts -> known d and se (50-digit decimal arithmetic of the same formulas); phasing's sum form prints
    identical sequences as -0.0; saturated or empty comparisons are NaN."""
    known = {(1000, 20, 10): (3.069732402315512882e-2, 5.649541502794334039e-3),
             (300, 6, 0): (2.041099726012756478e-2, 8.419691425682042399e-3),
             (2500, 150, 60): (9.004001246789040986e-2, 6.389521586585884052e-3),
             (80, 0, 3): (3.860099177752684617e-2, 2.251570257689479375e-2)}
    sites, ts, tv = (np.array(x) for x in zip(*known))
    d, se = k2p(sites, ts, tv)
    for (dk, sk), x, y in zip(known.values(), d, se):
        assert x == pytest.approx(dk, rel=1e-13) and y == pytest.approx(sk, rel=1e-12)
    d0, se0 = k2p([700], [0], [0])
    assert d0[0] == 0.0 and math.copysign(1.0, d0[0]) == -1.0 and se0[0] == 0.0
    d, se = k2p([100, 0], [50, 0], [1, 0])
    assert np.isnan(d).all() and np.isnan(se).all()


def test_k2p_counts_matches_formula():
    d, se = k2p(np.array([1000, 1000, 0]), np.array([20, 0, 0]), np.array([10, 0, 0]))
    P, Q = 0.02, 0.01
    assert abs(d[0] - (-0.5 * np.log(1 - 2 * P - Q) - 0.25 * np.log(1 - 2 * Q))) < 1e-12 and se[0] > 0
    assert d[1] == 0 and np.isnan(d[2])


def test_k2p_se_copy():
    """The annotation's standard error when it gives one (finite, > 0); else Poisson on the aligned LTR sites
    (sqrt(max(K2P, 1e-6) / sites)), the 5' LTR length standing in for missing sites, at least one site."""
    k = np.array([0.04, 0.04, 0.0, 0.09, np.nan, 0.01])
    k_se = np.array([0.003, np.nan, np.nan, 0.0, np.nan, np.nan])
    n_sites = np.array([400.0, 400.0, 250.0, np.nan, 300.0, np.nan])
    ltr_len = np.array([500.0, 500.0, 300.0, 900.0, 300.0, np.nan])
    se = k2p_se_copy(k, k_se, n_sites, ltr_len)
    want = [0.003, math.sqrt(0.04 / 400), math.sqrt(1e-6 / 250), math.sqrt(0.09 / 900), np.nan, math.sqrt(0.01)]
    assert se[:4] == pytest.approx(want[:4], rel=1e-15) and np.isnan(se[4]) and se[5] == pytest.approx(want[5])


def test_k2p_rejects_an_unknown_form():
    with pytest.raises(ValueError, match="form"):
        k2p([10], [1], [0], form="other")


def test_clean_tolerance_from_mapped_copies():
    assert clean_tolerance(np.zeros(50)) == cst.CLEAN_Z                 # too few mapped copies: the floor
    assert clean_tolerance(np.zeros(500)) == cst.CLEAN_Z                # aligners agree: sampling error only
    z = np.r_[np.zeros(400), np.full(100, 5.0)]                           # 20% of mapped copies differ by 5 SE
    assert clean_tolerance(z) == pytest.approx(5.0)


def test_self_ltr_z():
    """|K2P of a copy's own two LTRs - its annotated K2P| / SE: 0 when they agree (within one substitution), large when
    the annotation disagrees, inf without an annotated K2P; LTRs are the record slices l5 and l3 (0-based half-open)."""
    rng = np.random.default_rng(1)
    ltr, inner = rng.integers(0, 4, 300).astype(np.uint8), rng.integers(0, 4, 500).astype(np.uint8)
    other = ltr.copy()
    other[[10, 60, 110, 160, 210, 260]] = (other[[10, 60, 110, 160, 210, 260]] + 2) % 4    # six transitions
    recs = [np.concatenate((ltr, inner, ltr)), *[np.concatenate((ltr, inner, other))] * 3]
    lib = SimpleNamespace(codes=np.concatenate(recs), offs=np.arange(5, dtype=np.int64) * 1100,
                          l5=np.tile([0, 300], (4, 1)), l3=np.tile([800, 1100], (4, 1)))
    d6 = float(k2p([300], [6], [0])[0][0])
    z = self_ltr_z(lib, [0, 1, 2, 3], np.array([0.0, d6, 0.5, np.nan]))
    assert z[0] == 0 and z[1] == 0 and z[2] > 50 and z[3] == np.inf
    assert self_ltr_z(lib, [], np.zeros(4)).size == 0


def test_k2p_matches_ltrmix():
    """K2P and its standard error on LTRmix's inputs (whole counts; none, saturated and empty comparisons), one
    comparison at a time as the library route calls it: LTRmix's scalar math form bit for bit, d and se. On arrays d
    still agrees; se may differ in its last bit (numpy squares an array where Python's ** 2 calls pow; 1 in 3,000
    here)."""
    with np.load(Path(__file__).resolve().parent.parent / "golden" / "ltrmix_route.npz") as z:
        sites, ts, tv, want_d, want_se = (z[f"kp_{c}"] for c in ("sites", "ts", "tv", "d", "se"))
    one = [k2p(int(s), int(a), int(b)) for s, a, b in zip(sites, ts, tv)]
    d, se = (np.array([float(v[j]) for v in one]) for j in (0, 1))
    assert np.array_equal(d, want_d, equal_nan=True) and np.array_equal(se, want_se, equal_nan=True)
    assert np.array_equal(np.signbit(d), np.signbit(want_d)) and np.isnan(want_d).sum() > 100
    d, se = k2p(sites.astype(np.int64), ts.astype(np.int64), tv.astype(np.int64))
    assert np.array_equal(d, want_d, equal_nan=True) and np.allclose(se, want_se, rtol=1e-15, atol=0, equal_nan=True)


def test_k2p_shares_is_k2p_on_one_comparison():
    """k2p_shares (one comparison's shares, math's log) gives k2p's sum-form values bit for bit, and (nan, nan) where
    the comparison is saturated."""
    rng = np.random.default_rng(4)
    n = rng.integers(1, 3000, 2000)
    ts, tv = np.floor(rng.random(2000) * 0.4 * n), np.floor(rng.random(2000) * 0.3 * n)
    d, se = k2p(n, ts, tv)
    for i in range(n.size):
        got = k2p_shares(float(ts[i]) / int(n[i]), float(tv[i]) / int(n[i]), int(n[i]))
        if np.isnan(d[i]):
            assert all(math.isnan(x) for x in got)
        else:
            assert got == (float(d[i]), float(se[i])) or abs(got[1] - se[i]) <= 2 * np.spacing(se[i]), i
            assert got[0] == float(d[i]), i
    assert all(math.isnan(x) for x in k2p_shares(0.5, 0.0))  # 1 - 2P - Q = 0
    assert all(math.isnan(x) for x in k2p_shares(0.0, 0.5))  # 1 - 2Q = 0
    assert k2p_shares(0.0, 0.0) == (0.0, 0.0)
