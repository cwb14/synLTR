"""Divergence measures: the exact (Garwood) upper bound of a Poisson rate, the JC69 distance of an identity, the
Kimura two-parameter (K2P) distance of aligned sites with its standard error (`k2p` on arrays, `k2p_shares` on one
comparison's shares: one formula), the standard error of a copy's annotated LTR K2P (`k2p_se_copy`), and the check of
a copy's own two LTRs against its annotated K2P (`self_ltr_z`, `clean_tolerance`).

K2P has two algebraically equal forms that differ in the last bits of the result: "sum", -1/2 ln w1 - 1/4 ln w2, and
"product", -1/2 ln(w1 sqrt(w2)) + 0.0, with w1 = 1 - 2P - Q and w2 = 1 - 2Q. The sum form gives -0.0 for identical
sequences; the product form gives 0.0. A tool whose printed values were validated in one form keeps it through `form`.
"""
from __future__ import annotations

import math

import numpy as np

from ltrtrace.base.constants import CLEAN_CAL_MIN, CLEAN_PASS, CLEAN_Z, K2P_SE_FLOOR


def garwood_u95(k, n):
    """One-sided 95% upper bound of a Poisson rate from k events in n trials: chi2.ppf(0.95, 2(k + 1)) / 2 / n
    (exact, Garwood); k = 0 gives about 3 / n. Elementwise on arrays."""
    from scipy.stats import chi2

    k = np.asarray(k)
    return chi2.ppf(0.95, 2 * (k + 1)) / 2 / np.asarray(n)


def jc69_from_identity(identity: float) -> float | None:
    """The JC69 distance -3/4 ln(1 - 4p/3) of an identity, p = 1 - identity (EDTA's ltr_identity); None where it is
    undefined (p < 0 or 4p/3 >= 1). A nan identity gives nan. Never -0.0."""
    p = 1.0 - identity
    if p < 0 or 4 * p / 3 >= 1:
        return None
    return -0.75 * math.log(1 - 4 * p / 3) + 0.0


def k2p(sites, ts, tv, *, form: str = "sum") -> tuple[np.ndarray, np.ndarray]:
    """K2P distance and its standard error (delta method) from aligned sites, transitions and transversions,
    elementwise on arrays: P = ts / sites, Q = tv / sites, d = -1/2 ln(1 - 2P - Q) - 1/4 ln(1 - 2Q).

    `form` is "sum" or "product" (the module docstring). d and the standard error are NaN where the comparison is
    saturated (1 - 2P - Q <= 0 or 1 - 2Q <= 0) or has no sites.
    """
    if form not in ("sum", "product"):
        raise ValueError(f"form must be 'sum' or 'product', not {form!r}")
    sites = np.asarray(sites, float)
    n = np.maximum(sites, 1.0)
    P, Q = np.asarray(ts) / n, np.asarray(tv) / n
    with np.errstate(all="ignore"):
        d, v = _k2p_terms(P, Q, np.log, np.sqrt, form)
        var = v / n
    bad = (1 - 2 * P - Q <= 0) | (1 - 2 * Q <= 0) | (sites == 0)
    return np.where(bad, np.nan, d), np.where(bad, np.nan, np.sqrt(np.maximum(var, 0.0)))


def k2p_shares(P: float, Q: float, sites: float = 1.0) -> tuple[float, float]:
    """K2P distance (the "sum" form) and its standard error over `sites` sites from one comparison's transition and
    transversion shares P and Q, as Python floats with math's log (quick in a loop over copies, where k2p's arrays
    would cost ~40x more per call); (nan, nan) where saturated (1 - 2P - Q <= 0 or 1 - 2Q <= 0)."""
    if 1 - 2 * P - Q <= 0 or 1 - 2 * Q <= 0:
        return math.nan, math.nan
    d, v = _k2p_terms(P, Q, math.log, math.sqrt, "sum")
    return d, math.sqrt(max(v / sites, 0.0))


def _k2p_terms(P, Q, log, sqrt, form):
    """(K2P distance, sites x its delta-method variance) from transition and transversion shares: the formula once,
    for arrays (numpy's log and sqrt) and for numbers (math's)."""
    w1, w2 = 1 - 2 * P - Q, 1 - 2 * Q
    d = -0.5 * log(w1) - 0.25 * log(w2) if form == "sum" else -0.5 * log(w1 * sqrt(w2)) + 0.0
    c1, c2 = 1 / w1, 0.5 * (1 / w1 + 1 / w2)
    return d, c1 * c1 * P + c2 * c2 * Q - (c1 * P + c2 * Q) ** 2


def k2p_se_copy(k, k_se, n_sites, ltr_len):
    """Standard error of each copy's LTR K2P `k`: the annotation's own (`k_se`) where finite and > 0, else Poisson on
    the copy's aligned LTR sites, sqrt(max(k, 1e-6) / max(n, 1)), n = `n_sites` where finite and > 0, else the 5' LTR
    length (`ltr_len`; unknown = 0). Elementwise on arrays."""
    k, k_se, n_sites, ltr_len = (np.asarray(x, float) for x in (k, k_se, n_sites, ltr_len))
    n = np.where(np.isfinite(n_sites) & (n_sites > 0), n_sites, np.nan_to_num(ltr_len, nan=0.0))
    se = np.sqrt(np.maximum(k, K2P_SE_FLOOR) / np.maximum(n, 1.0))
    return np.where(np.isfinite(k_se) & (k_se > 0), k_se, se)


def self_ltr_z(lib, idx, k, threads=1):
    """For the copies `idx` of a library (codes, offs and the LTR record offsets l5, l3, 0-based half-open): |K2P
    between a copy's own two LTRs (at its current bounds) - its annotated LTR K2P (`k`, indexed by copy)| / its standard
    error; 0 within one substitution, inf when either is missing or saturated. The alignments run in up to `threads`
    processes; the kernels compile here first (a run's first numba compile happens before any fork: workers compiling at
    once read package metadata through shared file handles, 'BadZipFile: Overlapped entries')."""
    from ltrtrace.seq.align import pairs  # loads numba and edlib: only where the check runs

    idx = np.asarray(idx, np.int64)
    if idx.size == 0:
        return np.zeros(0)
    c = pairs(lib, idx, mode="self", threads=threads)
    d, se = k2p(c.sites, c.ts, c.tv)
    dk = np.abs(d - np.asarray(k, float)[idx])
    with np.errstate(invalid="ignore", divide="ignore"):
        z = np.where(dk <= 1.0 / np.maximum(c.sites, 1.0), 0.0, dk / se)
    return np.where(np.isfinite(dk), np.nan_to_num(z, nan=np.inf), np.inf)


def clean_tolerance(z_mapped):
    """The z a record-end LTR pair may reach: the 95th percentile over copies whose LTRs map from the genome (their
    bounds are right, so their z is aligner difference plus sampling error), never below CLEAN_Z; CLEAN_Z from fewer
    than 100 copies."""
    z = np.asarray(z_mapped, float)
    return CLEAN_Z if z.size < CLEAN_CAL_MIN else max(CLEAN_Z, float(np.quantile(z, CLEAN_PASS)))
