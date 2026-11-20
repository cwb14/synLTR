"""The merger clocks: a sibling pair of copies sits in two subgenomes only if its younger copy was inserted after the
merger. The LTR clock times the younger copy by its own LTRs; the genealogy clock by the pair's whole-record divergence
read through the genome's sibling genealogy. Profile likelihood over AGE_GRID (LTR divergence), lineage-cluster
bootstrap (stats.resample.cluster_bootstrap) widened by the model error (stats.resample.widen), as every dating
interval."""
from __future__ import annotations

import math

import numpy as np
from numba import njit, prange

from ltrtrace.analyses.subgenomes.inputs import warn
from ltrtrace.base.constants import (
    BF_DECISIVE,
    BOOT,
    CI95,
    CLOCK_MIN_AGES,
    EXCESS_FOLD,
    FIT_GTOL,
    FIT_ITERS,
    G_MODEL_DF,
    G_MODEL_OFFSET,
    G_MODEL_RIDGE,
    GENEALOGY_CHUNK,
    GENEALOGY_D_GRID,
    GENEALOGY_G_GRID,
    GENEALOGY_ITERS,
    GENEALOGY_R_GRID,
    LR95,
    MERGER_AGE_GRID,
    MODEL_ERR,
    SIB_MIN,
    Z95,
)
from ltrtrace.stats.resample import cluster_bootstrap, widen

AGE_GRID = np.concatenate(([0.0], np.geomspace(*MERGER_AGE_GRID)))   # LTR divergence grid of the merger clock


def ltr_counts(lib):
    """Substitutions and aligned sites between each copy's two LTRs (annotation counts; else from K2P and length)."""
    n = np.where(np.isfinite(lib.n_sites) & (lib.n_sites > 0), lib.n_sites, np.nan_to_num(lib.ltr5_len, nan=0.0))
    k = lib.n_ts + lib.n_tv
    with np.errstate(invalid="ignore"):
        p = _p_of(lib.k2p)                                       # NaN K2P stays NaN: no age, not age 0
    k = np.where(np.isfinite(k), k, np.round(p * n))
    return k, n


def _age_lik(k, n):
    """Binomial likelihood of k substitutions in n LTR sites at each grid divergence (rows: copies): _div_lik of the
    LTRs themselves."""
    return _div_lik(k, n, 1.0)


def age_prior(k, n, df=G_MODEL_DF, c=G_MODEL_RIDGE):
    """Smooth distribution of LTR divergence (insertion age) over AGE_GRID, by Efron g-modeling: g proportional to
    exp(Q beta), Q cubic B-splines on log(age + 1e-4), ridge penalty c |beta|^2, from the copies' counts."""
    from scipy.interpolate import BSpline
    from scipy.optimize import minimize
    L = _age_lik(k, n)
    x = np.log(AGE_GRID + G_MODEL_OFFSET)
    t = np.linspace(x.min(), x.max(), df - 2)
    Q = BSpline.design_matrix(np.clip(x, t[0], t[-1]), np.concatenate(([t[0]] * 3, t, [t[-1]] * 3)), 3).toarray()
    Q -= Q.mean(0)

    def f(b):
        e = Q @ b
        g = np.exp(e - e.max())
        g /= g.sum()
        m = L @ g
        w = g * (L.T @ (1.0 / m))
        return -(np.log(m).sum() - c * b @ b), -(Q.T @ w - L.shape[0] * (Q.T @ g) - 2 * c * b)
    r = minimize(f, np.zeros(Q.shape[1]), jac=True, method="L-BFGS-B")
    e = Q @ r.x
    g = np.exp(e - e.max())
    return g / g.sum()


@njit(cache=True)
def _two_rates_ll(B, y, phi, eps):
    """Log-likelihood, gradient and Hessian of P(cross) = phi B + eps (1 - B) (concave in phi, eps)."""
    f = 0.0
    g1 = 0.0
    g2 = 0.0
    h11 = 0.0
    h12 = 0.0
    h22 = 0.0
    for i in range(B.size):
        b = B[i]
        p = min(max(phi * b + eps * (1 - b), 1e-12), 1 - 1e-12)
        f += y[i] * np.log(p) + (1 - y[i]) * np.log(1 - p)
        r = y[i] / p - (1 - y[i]) / (1 - p)
        c = y[i] / (p * p) + (1 - y[i]) / ((1 - p) * (1 - p))
        g1 += r * b
        g2 += r * (1 - b)
        h11 -= c * b * b
        h12 -= c * b * (1 - b)
        h22 -= c * (1 - b) * (1 - b)
    return f, g1, g2, h11, h12, h22


@njit(cache=True)
def _phi_eps_fit(B, y):
    """Maximum likelihood (phi, eps) in [1e-9, 1 - 1e-9]^2 with phi >= eps: projected Newton with an active set (a
    parameter at a bound whose gradient points outward stays there) and backtracking (every step raises the
    likelihood). If the unconstrained maximum has phi < eps, the constrained one has phi = eps = mean(y).
    Returns (log-likelihood, phi, eps, converged)."""
    lo = 1e-9
    hi = 1 - 1e-9
    m = 0.0
    for i in range(y.size):
        m += y[i]
    m = min(max(m / max(y.size, 1), lo), hi)
    phi = min(max(2 * m, 0.05), 0.95)
    eps = min(max(m / 2, lo), hi)
    f, g1, g2, h11, h12, h22 = _two_rates_ll(B, y, phi, eps)
    converged = False
    for _ in range(FIT_ITERS):
        free1 = not ((phi <= lo and g1 < 0) or (phi >= hi and g1 > 0))
        free2 = not ((eps <= lo and g2 < 0) or (eps >= hi and g2 > 0))
        scale = 1.0 + abs(f)
        if (not free1 or abs(g1) < FIT_GTOL * scale) and (not free2 or abs(g2) < FIT_GTOL * scale):
            converged = True
            break
        d1 = 0.0
        d2 = 0.0
        det = h11 * h22 - h12 * h12
        if free1 and free2 and det > 1e-12 * (abs(h11 * h22) + 1e-300):
            d1 = -(h22 * g1 - h12 * g2) / det
            d2 = -(h11 * g2 - h12 * g1) / det
        elif free1 and h11 < 0:
            d1 = -g1 / h11
        elif free2 and h22 < 0:
            d2 = -g2 / h22
        else:                                                    # flat direction: gradient step
            d1 = g1 if free1 else 0.0
            d2 = g2 if free2 else 0.0
        step = 1.0
        moved = False
        while step > 1e-12:
            n1 = min(max(phi + step * d1, lo), hi)
            n2 = min(max(eps + step * d2, lo), hi)
            fn, a1, a2, b11, b12, b22 = _two_rates_ll(B, y, n1, n2)
            if fn >= f - 1e-12 * scale:
                moved = abs(n1 - phi) + abs(n2 - eps) > 1e-15
                phi, eps, f, g1, g2, h11, h12, h22 = n1, n2, fn, a1, a2, b11, b12, b22
                break
            step *= 0.5
        if not moved:
            converged = True
            break
    if phi < eps:                                                # one-sided: no young enrichment
        phi = eps = m
        f = _two_rates_ll(B, y, phi, eps)[0]
    return f, phi, eps, converged


@njit(parallel=True, cache=True)
def _clock_profile(B, y, cols):
    """Log-likelihood, phi, eps, converged at each candidate merger AGE_GRID[cols]; B[i, j] = P(pair i's younger copy
    was inserted less than AGE_GRID[j] ago)."""
    out = np.empty((cols.size, 4))
    for c in prange(cols.size):
        ll, phi, eps, ok = _phi_eps_fit(np.ascontiguousarray(B[:, cols[c]]), y)
        out[c, 0] = ll
        out[c, 1] = phi
        out[c, 2] = eps
        out[c, 3] = 1.0 if ok else 0.0
    return out


def _landing_bound(B, y, cols, phi, clusters):
    """Upper bound on the merger (LTR divergence units) when the clock sees no excess of moves, IF new copies land in
    either subgenome in proportion to its size: phi_i = the other subgenome's genome share seen from pair i's older
    copy. A merger at T predicts P(split) = phi_i P(younger < T) + eps (1 - P(younger < T)), eps = the observed split
    share; per lineage-cluster bootstrap draw (own generator: the dating's draws stay untouched) the oldest T within
    1.92 log-likelihood units of the best (T = 0 included); bound = 97.5% quantile widened by the merger model error
    (merger_clock floors it at the clock's resolution); None when even the oldest candidate survives in > 2.5% of
    draws. Fails when transposition stays within each subgenome after the merger (stated wherever it is reported)."""
    eps = min(max(float(y.mean()), 1e-6), 1 - 1e-6)
    Bm = B[:, cols]
    P = np.clip(phi[:, None] * Bm + eps * (1 - Bm), 1e-12, 1 - 1e-12)
    LL = y[:, None] * np.log(P) + (1 - y[:, None]) * np.log1p(-P)
    l0 = y * math.log(eps) + (1 - y) * math.log1p(-eps)
    rng = np.random.default_rng([int(y.size), int(y.sum())])
    W = np.array([np.bincount(sel, minlength=y.size) for sel in cluster_bootstrap(clusters, BOOT, rng)], float)
    llT, ll0 = W @ LL, W @ l0
    ok = llT >= (np.maximum(llT.max(1), ll0) - LR95)[:, None]
    last = np.array([cols[np.flatnonzero(o)[-1]] if o.any() else cols[0] for o in ok])
    if (last == cols[-1]).mean() > CI95[0]:
        return None
    return float(np.quantile(AGE_GRID[last], CI95[1]) * math.exp(Z95 * MODEL_ERR["merger"]))


def age_posteriors(k, n, g):
    """Each copy's insertion-age posterior over AGE_GRID (rows) from its LTR counts and the g-model prior g."""
    post = _age_lik(k, n) * g[None]
    return post / np.maximum(post.sum(1, keepdims=True), 1e-300)


def _div_lik(k, n, scale):
    """Binomial likelihood of k substitutions in n sites of a region evolving at 1/scale of the LTR rate, at each LTR
    divergence of AGE_GRID (rows max-normalised; n = 0 rows are flat)."""
    ll = _binom_ll(np.asarray(k, float)[:, None], np.asarray(n, float)[:, None], _p_of(AGE_GRID, scale)[None])
    return np.exp(ll - ll.max(1, keepdims=True))


def younger_cdf(k_u, n_u, k_v, n_v, g, kD=None, nD=None, scale=1.0):
    """B[i, j] = P(the younger copy of pair i was inserted less than AGE_GRID[j] ago).
    LTR ages only: 1 - S_u S_v. With the pair's own substitutions (kD of nD sites): the younger copy cannot
    be older than the pair's divergence (both copies' changes since the younger insertion are in it), so the joint
    age posterior is weighted by max over delta >= min(a_u, a_v) of Binomial(kD | nD, delta / scale), scale = LTR
    rate / rate of the aligned region (measured; LTRs diverge faster than whole elements): no assumption about how
    long before the younger copy the two lineages separated (pairs without counts, nD = 0, are unconstrained)."""
    Fu, Fv = age_posteriors(k_u, n_u, g), age_posteriors(k_v, n_v, g)
    Su = np.cumsum(Fu[:, ::-1], 1)[:, ::-1]                      # P(age >= AGE_GRID[j])
    Sv = np.cumsum(Fv[:, ::-1], 1)[:, ::-1]
    if kD is None:
        return np.ascontiguousarray(1.0 - Su * Sv)
    h = (Fu * Sv + Fv * Su - Fu * Fv) * np.maximum.accumulate(_div_lik(kD, nD, scale)[:, ::-1], 1)[:, ::-1]
    h /= np.maximum(h.sum(1, keepdims=True), 1e-300)             # P(min age = AGE_GRID[j] | LTRs, pair)
    return np.ascontiguousarray(np.clip(np.cumsum(h, 1) - h, 0.0, 1.0))


R_GRID = np.linspace(*GENEALOGY_R_GRID)     # the younger copy's age as a fraction of its pair's divergence
G_GRID = np.concatenate(([0.0], np.geomspace(*GENEALOGY_G_GRID)))   # pair divergence grid (LTR units)


def _binom_ll(k, n, p):
    """Binomial log-likelihood of k successes in n trials at probability p (broadcast)."""
    from scipy.special import gammaln
    return gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1) + k * np.log(p) + (n - k) * np.log1p(-p)


def _p_of(a, scale=1.0):
    """Expected substitution share at divergence a of a region evolving at 1/scale of the LTR rate (Jukes-Cantor),
    clipped away from 0 and 1."""
    return np.clip(0.75 * (1 - np.exp(-4 * np.asarray(a, float) / (3 * scale))), 1e-12, 1 - 1e-12)


def sibling_genealogy(k_u, n_u, k_v, n_v, delta, iters=GENEALOGY_ITERS):
    """f(r) over R_GRID: how far back, as a fraction r of the pair's divergence delta (LTR units), the younger copy of a
    sibling pair was inserted. Learned by EM from sibling pairs whose LTRs carry several substitutions; each pair's
    likelihood is L(r) = 1/2 sum over which copy is the younger of Binom(k_young | r delta) x max over d >= delta of
    Binom(k_older | d). r = 1: copied from its partner; r uniform: both drawn from one source lineage at random times.
    The genome's own genealogy, not an assumed one, then maps a young pair's divergence to its younger copy's age."""
    dgrid = np.geomspace(*GENEALOGY_D_GRID)
    lik = lambda k, n, a: np.exp(_binom_ll(k[:, None], n[:, None], _p_of(a)))
    Ly_u, Ly_v = lik(k_u, n_u, R_GRID[None] * delta[:, None]), lik(k_v, n_v, R_GRID[None] * delta[:, None])
    above = dgrid[None, :] >= delta[:, None]
    Lo_u = np.where(above, lik(k_u, n_u, dgrid[None]), 0.0).max(1)
    Lo_v = np.where(above, lik(k_v, n_v, dgrid[None]), 0.0).max(1)
    L = 0.5 * (Ly_u * Lo_v[:, None] + Ly_v * Lo_u[:, None]) + 1e-300
    f = np.full(R_GRID.size, 1.0 / R_GRID.size)
    for _ in range(iters):
        post = L * f[None]
        post /= post.sum(1, keepdims=True)
        f = post.mean(0)
    return f


def younger_cdf_genealogy(k_u, n_u, k_v, n_v, kD, nD, scale, f, chunk=GENEALOGY_CHUNK):
    """B[i, j] = P(the younger copy of pair i was inserted less than AGE_GRID[j] ago): the pair's divergence delta (LTR
    units; its kD of nD whole-record sites at delta / scale, uniform prior over G_GRID), the younger copy's age r delta
    with r from the learned genealogy f updated by both copies' LTR counts (either may be the younger)."""
    n = k_u.size
    B = np.zeros((n, AGE_GRID.size))
    w0 = np.diff(np.concatenate(([0.0], G_GRID)))
    A = R_GRID[None, :] * G_GRID[:, None]                         # younger ages (G, R)
    col = np.searchsorted(AGE_GRID, A, side="right").ravel()       # first grid age above each younger age
    for s in range(0, n, chunk):
        e = min(n, s + chunk)
        m = e - s
        Pd = np.exp(_binom_ll(kD[s:e, None], nD[s:e, None], _p_of(G_GRID[None] / scale))) * w0[None]
        Pd /= np.maximum(Pd.sum(1, keepdims=True), 1e-300)
        lk = lambda k, nn: np.exp(_binom_ll(k[:, None, None], nn[:, None, None], _p_of(A[None])))
        Ou = np.exp(_binom_ll(k_u[s:e, None], n_u[s:e, None], _p_of(G_GRID[None])))[:, :, None]
        Ov = np.exp(_binom_ll(k_v[s:e, None], n_v[s:e, None], _p_of(G_GRID[None])))[:, :, None]
        Pr = f[None, None, :] * 0.5 * (lk(k_u[s:e], n_u[s:e]) * Ov + lk(k_v[s:e], n_v[s:e]) * Ou)
        Pr /= np.maximum(Pr.sum(2, keepdims=True), 1e-300)          # P(r | delta, LTRs)
        Wt = (Pd[:, :, None] * Pr).reshape(m, -1)
        H = np.zeros((m, AGE_GRID.size + 1))
        np.add.at(H, (np.repeat(np.arange(m), col.size), np.tile(col, m)), Wt.ravel())
        B[s:e] = np.cumsum(H[:, :-1], 1)
    return np.ascontiguousarray(np.clip(B, 0.0, 1.0))


def merger_clock(B, cross, clusters, rng, t_max, t_res, phi_land=None):
    """Pairwise sibling clock: a sibling pair sits in two subgenomes only when its younger copy was inserted after the
    merger (eps: exchange or misassignment noise, age-independent). T (LTR divergence units) by profile likelihood
    over AGE_GRID below t_max, given B = P(younger copy inserted < T) per pair (younger_cdf); lineage-cluster bootstrap
    widened by MODEL_ERR["merger"]. Not datable unless cross pairs are >= 2-fold enriched among post-merger pairs
    (phi >= 2 eps) in >= 97.5% of resamples; an upper bound only ("younger than") when the interval reaches the
    clock's resolution t_res (one substitution per LTR for LTR ages, per aligned pair with the pair divergence)."""
    y = cross.astype(np.float64)
    cols = np.flatnonzero((AGE_GRID > 0) & (AGE_GRID < t_max)).astype(np.int64)
    if cols.size < CLOCK_MIN_AGES or y.size < SIB_MIN:
        return None
    prof = _clock_profile(B, y, cols)
    i = int(np.argmax(prof[:, 0]))
    T, phi, eps = AGE_GRID[cols[i]], prof[i, 1], prof[i, 2]
    Tb, gap, bad = [], [], int((prof[:, 3] == 0).sum())
    for sel in cluster_bootstrap(clusters, BOOT, rng):
        pb = _clock_profile(np.ascontiguousarray(B[sel]), y[sel], cols)
        bad += int((pb[:, 3] == 0).sum())
        j = int(np.argmax(pb[:, 0]))
        Tb.append(AGE_GRID[cols[j]])
        gap.append(pb[j, 1] - EXCESS_FOLD * pb[j, 2])
    lo, hi = widen((T, *np.quantile(Tb, CI95)), MODEL_ERR["merger"])
    bound = bool(np.quantile(Tb, CI95[0]) <= t_res)              # resolution judged on the bootstrap alone
    m_ = min(max(float(y.mean()), 1e-12), 1 - 1e-12)
    log_lr = float(prof[i, 0] - (y * math.log(m_) + (1 - y) * math.log1p(-m_)).sum())   # step vs no step
    datable = bool(phi >= EXCESS_FOLD * eps and np.quantile(gap, CI95[0]) > 0 and log_lr >= math.log(BF_DECISIVE))
    if bad:
        warn(f"the merger clock's rate fit did not converge in {bad} of {cols.size * (BOOT + 1)} fits")
    lb = _landing_bound(B, y, cols, np.asarray(phi_land, float), clusters) if (
        phi_land is not None and not datable and float(np.mean(phi_land)) >= EXCESS_FOLD * float(y.mean())) else None
    if lb is not None:                                     # (random landing must predict >= 2x the background,
        lb = max(lb, t_res)                                #  the clock's rule; resolution floor: one substitution
        lb = lb if lb < t_max else None                    #  per LTR; a bound at or above the split says nothing)
    return dict(landing_bound=lb, estimate=float(T), lo=0.0 if bound else float(lo), hi=float(hi), phi=float(phi),
                eps=float(eps), datable=datable, bound=bound, n=int(y.size), n_cross=int(y.sum()), t_res=float(t_res),
                log_lr=log_lr)
