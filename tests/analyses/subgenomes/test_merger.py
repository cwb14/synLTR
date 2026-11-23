"""The merger clocks: age posteriors, the sibling genealogy, the two-rate fit and its bounds."""
import math

import numpy as np
import pytest

from synltr.analyses.subgenomes import inputs, merger
from synltr.base import constants as cst


def _counts(age, n, rng):
    """LTR substitution counts of copies of the given ages (K2P units) in n sites each."""
    p = 0.75 * (1 - np.exp(-4 * np.asarray(age) / 3))
    return rng.binomial(n, p).astype(float), np.full(np.size(age), float(n))


def test_younger_cdf_ltr_only_is_one_minus_both_survivals():
    rng = np.random.default_rng(1)
    k, n = _counts(rng.uniform(0, 0.02, 40), 600, rng)
    g = np.full(merger.AGE_GRID.size, 1.0 / merger.AGE_GRID.size)
    B = merger.younger_cdf(k[:20], n[:20], k[20:], n[20:], g)
    post = lambda kk, nn: merger.age_posteriors(kk, nn, g)
    Su = np.cumsum(post(k[:20], n[:20])[:, ::-1], 1)[:, ::-1]
    Sv = np.cumsum(post(k[20:], n[20:])[:, ::-1], 1)[:, ::-1]
    assert np.allclose(B, 1 - Su * Sv)
    assert np.all(np.diff(B, axis=1) >= -1e-12)          # P(younger < T) never decreases with T


def test_pair_divergence_only_makes_younger_copies_younger():
    """D bounds the younger copy's age from above: B with D >= B without, and near 1 above a tiny divergence."""
    rng = np.random.default_rng(2)
    k, n = _counts(np.full(30, 0.001), 500, rng)
    g = np.full(merger.AGE_GRID.size, 1.0 / merger.AGE_GRID.size)
    B0 = merger.younger_cdf(k[:15], n[:15], k[15:], n[15:], g)
    B1 = merger.younger_cdf(k[:15], n[:15], k[15:], n[15:], g, np.zeros(15), np.full(15, 8000.0))
    assert np.all(B1 >= B0 - 1e-9)
    j = np.searchsorted(merger.AGE_GRID, 2e-3)
    assert B1[:, j].min() > 0.95                         # 0 of 8,000 sites differ: younger copy < ~1e-3


def _sibling_pairs(r_draw, m, rng, n_ltr=800):
    delta = rng.uniform(0.006, 0.05, m)
    a_young = r_draw(m) * delta
    ku, nu = _counts(a_young, n_ltr, rng)
    kv, nv = _counts(delta, n_ltr, rng)                    # the older copy carries the pair's whole history
    return ku, nu, kv, nv, delta


@pytest.mark.parametrize("kind", ["copy_of_copy", "one_source"])
def test_sibling_genealogy_recovers_r(kind):
    rng = np.random.default_rng(3)
    draw = (lambda m: np.ones(m)) if kind == "copy_of_copy" else (lambda m: rng.uniform(0, 1, m))
    f = merger.sibling_genealogy(*_sibling_pairs(draw, 3000, rng))
    mean_r = float(np.dot(merger.R_GRID, f))
    assert abs(f.sum() - 1) < 1e-9
    if kind == "copy_of_copy":
        assert mean_r > 0.85
    else:
        assert 0.35 < mean_r < 0.65


def test_genealogy_cdf_monotone_and_identical_pairs_young():
    rng = np.random.default_rng(4)
    f = np.full(merger.R_GRID.size, 1.0 / merger.R_GRID.size)
    k, n = _counts(np.zeros(10), 600, rng)
    B = merger.younger_cdf_genealogy(k, n, k.copy(), n.copy(), np.zeros(10), np.full(10, 6000.0), 1.0, f)
    assert np.all(np.diff(B, axis=1) >= -1e-9)
    assert np.all(B[:, -1] <= 1 + 1e-9)
    j = np.searchsorted(merger.AGE_GRID, 1e-3)
    assert B[:, j].min() > 0.9                           # identical over 6 kb: the younger copy is < 1e-3


def _clock_data(rng, m=3000, T=2e-3, phi=0.5, eps=0.03, step=True):
    a = rng.uniform(0, 0.02, m)                          # younger copy ages, known exactly (B is a step at a)
    B = (merger.AGE_GRID[None, :] > a[:, None]).astype(float)
    p = np.where(a < T, phi, eps) if step else np.full(m, eps)
    y = rng.random(m) < p
    clusters = rng.integers(0, 300, m)
    return B, y, clusters


def test_clock_dates_a_real_step():
    rng = np.random.default_rng(5)
    B, y, cl = _clock_data(rng)
    mg = merger.merger_clock(np.ascontiguousarray(B), y, cl, rng, 0.05, 1e-5)
    assert mg["datable"] and not mg["bound"]
    assert mg["lo"] <= 2e-3 <= mg["hi"]
    assert mg["log_lr"] > math.log(cst.BF_DECISIVE)


def test_clock_needs_decisive_evidence():
    """No step (cross pairs at the background rate at every age): never datable, whatever the bootstrap does."""
    rng = np.random.default_rng(6)
    B, y, cl = _clock_data(rng, step=False, eps=0.01)
    mg = merger.merger_clock(np.ascontiguousarray(B), y, cl, rng, 0.05, 1e-5)
    assert not mg["datable"]


def test_ltr_counts_fallback_from_k2p():
    """Annotations without substitution counts (K2P and LTR length only): counts come from K2P and length."""
    lib = inputs.Library()
    lib.k2p = np.array([0.01, 0.02])
    lib.n_ts, lib.n_tv = np.array([np.nan, 3.0]), np.array([np.nan, 1.0])
    lib.n_sites = np.array([np.nan, 400.0])
    lib.ltr5_len = np.array([500.0, 400.0])
    k, n = merger.ltr_counts(lib)
    assert n.tolist() == [500.0, 400.0] and k[1] == 4.0 and k[0] == round(0.75 * (1 - np.exp(-4 * 0.01 / 3)) * 500)


def test_ltr_counts_missing_k2p_is_missing():
    """A copy without K2P and without counts has no age: NaN, so it never enters the age prior as age 0."""
    lib = inputs.Library()
    lib.k2p = np.array([np.nan, 0.02])
    lib.n_ts, lib.n_tv = np.array([np.nan, np.nan]), np.array([np.nan, np.nan])
    lib.n_sites = np.array([np.nan, np.nan])
    lib.ltr5_len = np.array([450.0, 450.0])
    k, n = merger.ltr_counts(lib)
    assert np.isnan(k[0]) and k[1] == round(0.75 * (1 - np.exp(-4 * 0.02 / 3)) * 450)


def _planted_clock(rng, T, phi=0.5, eps=0.03, n=6000, n_sites=1000):
    """Sibling pairs: ages (LTR divergence units) log-uniform 1e-4..0.05; cross with probability phi when the
    younger copy is younger than T, else eps; binomial LTR substitution counts."""
    a_u = 10 ** rng.uniform(-4, np.log10(0.05), n)
    a_v = 10 ** rng.uniform(-4, np.log10(0.05), n)
    young = np.minimum(a_u, a_v) < T
    cross = rng.random(n) < np.where(young, phi, eps)
    p = lambda a: 0.75 * (1 - np.exp(-4 * a / 3))
    k_u, k_v = rng.binomial(n_sites, p(a_u)), rng.binomial(n_sites, p(a_v))
    lib_k = rng.binomial(n_sites, p(10 ** rng.uniform(-4, np.log10(0.05), 4000)))
    return k_u, np.full(n, float(n_sites)), k_v, np.full(n, float(n_sites)), cross, lib_k


def _ltr_clock(k_u, n_u, k_v, n_v, cross, cl, g, rng, t_max, phi_land=None):
    """The LTR clock: the younger copy's age matrix built first (younger_cdf)."""
    t_res = 1.0 / np.median(np.concatenate([n_u, n_v]))
    return merger.merger_clock(merger.younger_cdf(k_u, n_u, k_v, n_v, g), cross, cl, rng, t_max, t_res, phi_land=phi_land)


def test_merger_clock_planted_step():
    rng = np.random.default_rng(8)
    T = 0.004
    k_u, n_u, k_v, n_v, cross, lib_k = _planted_clock(rng, T)
    g = merger.age_prior(lib_k, np.full(lib_k.size, 1000.0))
    r = _ltr_clock(k_u, n_u, k_v, n_v, cross, np.arange(k_u.size) // 20, g, rng, 0.05)
    assert r["datable"] and not r["bound"] and r["log_lr"] >= math.log(cst.BF_DECISIVE)
    assert 0.75 * T <= r["estimate"] <= 1.33 * T and r["lo"] <= T <= r["hi"]
    assert 0.4 <= r["phi"] <= 0.6 and r["eps"] < 0.06


def test_merger_clock_no_step_and_too_young():
    rng = np.random.default_rng(9)
    k_u, n_u, k_v, n_v, cross, lib_k = _planted_clock(rng, 0.004, phi=0.03, n=3000)  # no enrichment: not datable
    g = merger.age_prior(lib_k, np.full(lib_k.size, 1000.0))
    r = _ltr_clock(k_u, n_u, k_v, n_v, cross, np.arange(k_u.size) // 20, g, rng, 0.05)
    assert not r["datable"]
    k_u, n_u, k_v, n_v, cross, lib_k = _planted_clock(rng, 1e-6, n=3000)            # merger younger than resolution
    cross = np.where(rng.random(cross.size) < 0.01, ~cross, cross)                 # a few post-merger crossings
    k_u[:30], k_v[:30], cross[:30] = 0, 0, True
    r = _ltr_clock(k_u, n_u, k_v, n_v, cross, np.arange(k_u.size) // 20, g, rng, 0.05)
    assert not r["datable"] or (r["bound"] and r["lo"] == 0.0)


def _ll_two_rates(phi, eps, B, y):
    p = np.clip(phi * B + eps * (1 - B), 1e-12, 1 - 1e-12)
    return float((y * np.log(p) + (1 - y) * np.log1p(-p)).sum())


def test_phi_eps_fit_reaches_the_constrained_maximum():
    """The two-rate fit (cross probability phi for post-merger pairs, eps otherwise, phi >= eps) reaches the maximum
    likelihood, including at the bounds (eps = 0: no exchange noise; phi = 1), and says it converged."""
    from scipy.optimize import minimize
    rng = np.random.default_rng(17)
    for it in range(60):
        n = int(rng.choice([60, 300, 3000]))
        B = (rng.random(n), rng.beta(0.2, 5, n), rng.beta(5, 0.2, n), (rng.random(n) < 0.3).astype(float))[it % 4]
        phi_t, eps_t = rng.choice([0.0, 0.01, 0.1, 0.5, 0.95]), rng.choice([0.0, 0.001, 0.02, 0.1])
        y = (rng.random(n) < phi_t * B + eps_t * (1 - B)).astype(float)
        ll, phi, eps, ok = merger._phi_eps_fit(B, y)
        g = np.concatenate([np.geomspace(1e-6, 0.5, 60), 1 - np.geomspace(1e-6, 0.5, 60)[::-1]])
        best, x0 = max((_ll_two_rates(a, b, B, y), (a, b)) for a in g for b in g[g <= a][::3])
        r = minimize(lambda x: -_ll_two_rates(x[0], x[1], B, y), x0, method="L-BFGS-B", bounds=[(1e-9, 1 - 1e-9)] * 2)
        if r.x[0] >= r.x[1]:
            best = max(best, -r.fun)
        best = max(best, _ll_two_rates(y.mean(), y.mean(), B, y))
        assert ok and phi >= eps and ll >= best - 1e-3, (it, ll, best, phi, eps)


def test_merger_clock_landing_bound():
    """No excess of moves (merger younger than every sibling pair): not datable, but if new copies land at random
    the merger is younger than an informative bound that contains it; a datable merger carries no such bound; no
    random-landing shares given, no bound."""
    rng = np.random.default_rng(12)
    k_u, n_u, k_v, n_v, cross, lib_k = _planted_clock(rng, 5e-5, n=6000)
    g = merger.age_prior(lib_k, np.full(lib_k.size, 1000.0))
    cl = np.arange(k_u.size) // 20
    r = _ltr_clock(k_u, n_u, k_v, n_v, cross, cl, g, rng, 0.05, phi_land=np.full(k_u.size, 0.5))
    assert not r["datable"] and r["landing_bound"] is not None and 5e-5 <= r["landing_bound"] <= 2e-3
    assert _ltr_clock(k_u, n_u, k_v, n_v, cross, cl, g, rng, 0.05)["landing_bound"] is None
    k_u, n_u, k_v, n_v, cross, lib_k = _planted_clock(rng, 0.004)
    r = _ltr_clock(k_u, n_u, k_v, n_v, cross, cl, g, rng, 0.05, phi_land=np.full(k_u.size, 0.5))
    assert r["datable"] and r["landing_bound"] is None


def test_landing_bound_needs_room_and_a_signal():
    """The random-landing bound is dropped when it would not lie below the split (t_max: nothing learned) and when
    random landing predicts fewer than twice the background cross pairs (a bound resting on a deficit)."""
    rng = np.random.default_rng(12)
    k_u, n_u, k_v, n_v, cross, lib_k = _planted_clock(rng, 5e-5, n=6000)
    g = merger.age_prior(lib_k, np.full(lib_k.size, 1000.0))
    cl = np.arange(k_u.size) // 20
    r = _ltr_clock(k_u, n_u, k_v, n_v, cross, cl, g, rng, 5e-4, phi_land=np.full(k_u.size, 0.5))
    assert r is not None and not r["datable"] and r["landing_bound"] is None      # floor 1e-3 >= t_max 5e-4
    r = _ltr_clock(k_u, n_u, k_v, n_v, cross, cl, g, rng, 0.05, phi_land=np.full(k_u.size, 0.05))
    assert not r["datable"] and r["landing_bound"] is None                         # 0.05 < 2 x eps (0.03)
