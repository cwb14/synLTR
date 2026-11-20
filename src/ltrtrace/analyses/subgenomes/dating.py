"""Dating of each pair of subgenomes from copy pairs (closest relatives across phased sequences, whole records
aligned): split estimators (inherited insertions, lineage profile, confinement, duplicated cross pairs; reported:
inherited, else their log-scale median) and the merger (sibling pairs on two clocks), all on the LTR clock (K2P = 2 mu
years). dating.tsv and summary.json carry every estimate."""
from __future__ import annotations

from pathlib import Path

import numpy as np

from ltrtrace.analyses.subgenomes.exchanges import flag_homogenized
from ltrtrace.analyses.subgenomes.inputs import _f, group_label, log, vlog
from ltrtrace.analyses.subgenomes.merger import (
    R_GRID,
    age_prior,
    ltr_counts,
    merger_clock,
    sibling_genealogy,
    younger_cdf,
    younger_cdf_genealogy,
)
from ltrtrace.analyses.subgenomes.painting import copy_ancestry, run_painting
from ltrtrace.base.constants import (
    BOOT,
    BOOT_MIN,
    CI95,
    CLUSTER_D,
    CONF_BINS,
    CONF_GRID_POINTS,
    CONF_MIN_COPIES,
    CONF_SE_FLOOR,
    CONF_TOP_QUANTILE,
    DATE_MIN_COPIES,
    DATE_MIN_PAIRS,
    FOCAL,
    GEN_MIN_PAIRS,
    GEN_MIN_SUBS,
    HALF_COUNT,
    JC_P_CAP,
    MODEL_ERR,
    MU_DEFAULT,
    N_ORTH,
    PAIR_MIN_SITES,
    PLATEAU_ALPHA,
    PRIOR_POOL_MAX,
    PROF_LOG10_EDGES,
    PROF_MIN_PAIRS,
    RATE_RATIO_MAX,
    REL_TOP,
    RISE_MIN,
    SCALE_NEAR_FACTOR,
    SCALE_NEAR_MAX,
    SCALE_NEAR_MIN,
    SE_MULT,
    SG_LTR_MIN_SITES,
    SIB_MIN,
    SKETCH_K,
    SKETCH_S,
    SUPPORT_OK,
    TROUGH_POOL,
    YOUNG_CUT_TOL,
    YOUNG_MAX,
    YOUNG_QUANTILE,
)
from ltrtrace.seq import align, divergence
from ltrtrace.seq.divergence import k2p_se_copy
from ltrtrace.seq.kmers import sketch
from ltrtrace.seq.relatives import top_shared
from ltrtrace.stats.resample import cluster_bootstrap, widen


def _counts(lib, P, mode, threads):
    """(aligned sites, transitions, transversions) of the copy pairs P (rows), seq.align.pairs' mode."""
    return np.stack(align.pairs(lib.store, P, mode=mode, threads=threads)[:3], 1)


def _se(lib):
    """Standard error of each copy's LTR K2P (seq.divergence.k2p_se_copy of the annotation's columns)."""
    return k2p_se_copy(lib.k2p, lib.k2p_se, lib.n_sites, lib.ltr5_len)


def sample_pairs(lib, ek, rng, threads):
    """Copy pairs for dating: FOCAL random phased copies plus FOCAL among the youngest 20% (LTR K2P), each with its
    REL_TOP closest relatives on other phased sequences (shared 15-mer sketch hashes), aligned as whole records.
    Returns dict(P=pairs (focal first), d=pair K2P, se, stratum (0 random, 1 young))."""
    phased = (ek >= 0) & np.isfinite(lib.k2p)
    cand = np.flatnonzero(phased)
    if cand.size < DATE_MIN_COPIES:
        return None
    S = sketch(lib.store, k=SKETCH_K, s=SKETCH_S)
    young = cand[lib.k2p[cand] <= np.quantile(lib.k2p[cand], YOUNG_QUANTILE)]
    fa = rng.choice(cand, min(FOCAL, cand.size), replace=False)
    fy = rng.choice(young, min(FOCAL, young.size), replace=False)
    focal = np.concatenate([fa, fy])
    stratum = np.concatenate([np.zeros(fa.size, np.int8), np.ones(fy.size, np.int8)])
    J, _ = top_shared(S, np.where(phased, ek, -1).astype(np.int64), focal, m=REL_TOP)
    r, q = np.nonzero(J >= 0)
    P = np.stack([focal[r], J[r, q]], 1)
    key = np.minimum(P[:, 0], P[:, 1]) * lib.n + np.maximum(P[:, 0], P[:, 1])
    _, first = np.unique(key, return_index=True)
    first.sort()
    P, st = P[first], stratum[r[first]]
    cnt = _counts(lib, P, "whole", threads)
    d, se = divergence.k2p(cnt[:, 0], cnt[:, 1], cnt[:, 2])
    ok = np.isfinite(d) & (cnt[:, 0] >= PAIR_MIN_SITES)     # whole elements
    return dict(P=P[ok], d=d[ok], se=se[ok], stratum=st[ok], cnt=cnt[ok], sketch=S, young=young, focal_young=fy)


def sample_young_more(lib, ek, pt, rng, threads):
    """More young focal copies for the merger clock. The youngest 20% of phased copies beyond pt's FOCAL young ones,
    up to YOUNG_MAX in all (drawn on their own stream: pt's draws and table are untouched), each with its REL_TOP
    closest relatives on other sequences, aligned like pt's pairs; pairs already in pt are not repeated.
    LTR-RT-rich libraries thus contribute more sibling pairs (narrower intervals, younger mergers resolved); the
    budget bounds run time. Returns pt extended (stratum 2 = extra young focal copies)."""
    phased = (ek >= 0) & np.isfinite(lib.k2p)
    rest = np.setdiff1d(pt["young"], pt["focal_young"])
    n_extra = int(min(max(YOUNG_MAX - pt["focal_young"].size, 0), rest.size))
    out = dict(P=pt["P"], d=pt["d"], se=pt["se"], stratum=pt["stratum"], cnt=pt["cnt"])
    if n_extra == 0:
        return out
    fy2 = np.sort(rng.choice(rest, n_extra, replace=False))
    J, _ = top_shared(pt["sketch"], np.where(phased, ek, -1).astype(np.int64), fy2, m=REL_TOP)
    r, q = np.nonzero(J >= 0)
    P = np.stack([fy2[r], J[r, q]], 1)
    key = np.minimum(P[:, 0], P[:, 1]) * lib.n + np.maximum(P[:, 0], P[:, 1])
    k7 = np.minimum(pt["P"][:, 0], pt["P"][:, 1]) * lib.n + np.maximum(pt["P"][:, 0], pt["P"][:, 1])
    _, first = np.unique(key, return_index=True)
    first = np.sort(first[~np.isin(key[first], k7)])
    P = P[first]
    if not len(P):
        return out
    cnt = _counts(lib, P, "whole", threads)
    d, se = divergence.k2p(cnt[:, 0], cnt[:, 1], cnt[:, 2])
    ok = np.isfinite(d) & (cnt[:, 0] >= PAIR_MIN_SITES)
    return dict(P=np.concatenate([pt["P"], P[ok]]), d=np.concatenate([pt["d"], d[ok]]),
                se=np.concatenate([pt["se"], se[ok]]), stratum=np.concatenate([pt["stratum"], np.full(ok.sum(), 2, np.int8)]),
                cnt=np.concatenate([pt["cnt"], cnt[ok]]))


def ltr_scale_near(lib, P, d, cnt, est, rng, threads, n_max=SCALE_NEAR_MAX):
    """LTR / whole-record K2P of the pairs whose whole-record divergence is within a factor of 2 of est (pooled
    counts of their 5'-5' + 3'-3' LTR alignments over those of their whole records); None with < 30 such pairs."""
    sel = np.flatnonzero((d > est / SCALE_NEAR_FACTOR) & (d < est * SCALE_NEAR_FACTOR))
    if sel.size < SCALE_NEAR_MIN:
        return None
    sel = rng.choice(sel, min(sel.size, n_max), replace=False)
    lc = _counts(lib, P[sel], "ltr", threads)
    ok = lc[:, 0] >= SG_LTR_MIN_SITES
    if ok.sum() < SCALE_NEAR_MIN:
        return None
    dl = divergence.k2p(lc[ok, 0].sum(), lc[ok, 1].sum(), lc[ok, 2].sum())[0]
    dw = divergence.k2p(cnt[sel][ok, 0].sum(), cnt[sel][ok, 1].sum(), cnt[sel][ok, 2].sum())[0]
    return float(dl / dw) if np.isfinite(dl) and np.isfinite(dw) and dw > 0 else None


def ltr_rate_ratio(lib, P, cnt, rng, threads, n_max=RATE_RATIO_MAX):
    """LTR rate / whole-record rate of sibling copy pairs: pooled substitutions per site of their 5'-5' + 3'-3' LTR
    alignments over those of their whole-record alignments (up to n_max pairs, both with >= 200 LTR sites and some
    substitutions); 1 when undeterminable. Puts pair divergences on the LTR-age scale of the clock."""
    if not len(P):
        return 1.0
    i = rng.choice(len(P), min(len(P), n_max), replace=False)
    lc = _counts(lib, P[i], "ltr", threads)
    ok = lc[:, 0] >= SG_LTR_MIN_SITES
    pl = (lc[ok, 1] + lc[ok, 2]).sum() / max(lc[ok, 0].sum(), 1)
    pw = (cnt[i][ok, 1] + cnt[i][ok, 2]).sum() / max(cnt[i][ok, 0].sum(), 1)
    return float(pl / pw) if pl > 0 and pw > 0 else 1.0


def pair_units(pt, n):
    """Bootstrap units from the pair table pt: components of copies joined by pairs closer than CLUSTER_D (pure
    families in simulations), and which copies appear in it. Denser sampling must not chain families into one unit, so
    pairs added for the merger clock take the unit of their focal copy if it is in pt, else of their partner."""
    comp = lineage_clusters(pt["P"], pt["d"], n, CLUSTER_D)
    in7 = np.zeros(n, bool)
    in7[pt["P"].ravel()] = True
    return comp, in7


PROF_EDGES = np.linspace(*PROF_LOG10_EDGES)  # log10 pair-divergence bins of the cross-subgenome profile


def cross_profile(d, cross):
    """Pairs and cross-subgenome pairs per log10-divergence bin (pairs at d = 0 go to the first bin)."""
    j = np.clip(np.searchsorted(PROF_EDGES, np.log10(np.maximum(d, 1e-12))) - 1, 0, PROF_EDGES.size - 2)
    return np.bincount(j, minlength=PROF_EDGES.size - 1), np.bincount(j, weights=cross, minlength=PROF_EDGES.size - 1)


def _profile_trough(d, cross):
    """Bins (n, cross fraction f, bins with >= 20 pairs) and the trough: the bin of lowest 3-bin pooled cross
    fraction (the confined plateau between the merger and the split), sought over every bin: below the median pair
    divergence, which young-rich libraries (and the merger table's extra young focal copies) pull below the plateau,
    the trough would fall on noise among post-merger pairs (no profile split, a near-zero merger cap)."""
    n, k = cross_profile(d, cross)
    f = np.where(n >= PROF_MIN_PAIRS, k / np.maximum(n, 1), np.nan)
    ok = np.flatnonzero(np.isfinite(f))
    if not ok.size:
        return n, k, f, ok, None
    w = np.ones(TROUGH_POOL)
    fs = np.convolve(k, w, "same") / np.maximum(np.convolve(n, w, "same"), 1)
    return n, k, f, ok, int(ok[np.argmin(fs[ok])])


def young_cut(d, cross):
    """Upper edge of the young (post-merger) regime: the first bin from below whose cross fraction is within 0.02 of
    the trough."""
    n, k, f, ok, b0 = _profile_trough(d, cross)
    if b0 is None:
        return 0.0
    b = next(b for b in ok if f[b] <= f[b0] + YOUNG_CUT_TOL)
    return 10 ** PROF_EDGES[b]


def plateau_end(d, cross):
    """End of the confined plateau: scanning up from the trough, bins join the plateau while their cross fraction is
    not significantly above the plateau's pooled fraction (one-sided binomial, Bonferroni over the bins scanned).
    Below it no pair needs to predate the split."""
    from scipy.stats import norm
    n, k, f, ok, b0 = _profile_trough(d, cross)
    if b0 is None:
        return float(np.median(d)) if d.size else 0.0
    up = ok[ok > b0]
    zc = float(norm.isf(PLATEAU_ALPHA / max(up.size, 1)))
    kk, nn, last = k[b0], n[b0], b0
    for b in up:
        p0 = min(max(kk / nn, HALF_COUNT / nn), 1 - HALF_COUNT / nn)
        if (f[b] - p0) / np.sqrt(p0 * (1 - p0) / n[b]) >= zc:
            break
        kk, nn, last = kk + k[b], nn + n[b], b
    return 10 ** PROF_EDGES[last + 1]


def _iso_increasing(y, w):
    """Pool-adjacent-violators (non-decreasing fit)."""
    blocks = []
    for v, ww in zip(map(float, y), map(float, w)):
        blocks.append([v, ww, 1])
        while len(blocks) > 1 and blocks[-2][0] > blocks[-1][0]:
            a, b = blocks.pop(), blocks.pop()
            W = a[1] + b[1]
            blocks.append([(a[0] * a[1] + b[0] * b[1]) / max(W, 1e-12), W, a[2] + b[2]])
    return np.concatenate([[v] * c for v, _, c in blocks]) if blocks else np.array([])


def half_rise(d, cross):
    """Lineage-floor split: the divergence where the (monotone) cross fraction above the young cut climbs halfway
    from the confined plateau to the ancestral level; NaN when it does not climb by >= 0.05."""
    n, k, f, ok, b0 = _profile_trough(d, cross)
    if b0 is None:
        return float("nan")
    above = ok[ok >= b0]
    iso = _iso_increasing(f[above], n[above])
    if iso.size < 2 or iso[-1] - iso[0] < RISE_MIN:
        return float("nan")
    b = above[np.argmax(iso >= iso[0] + 0.5 * (iso[-1] - iso[0]))]
    return float(10 ** (0.5 * (PROF_EDGES[b] + PROF_EDGES[b + 1])))


def lineage_clusters(P, d, n, max_d):
    """Bootstrap units: connected components of copies joined by pairs closer than max_d (pure families in sims)."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    m = d < max_d
    A = coo_matrix((np.ones(m.sum()), (P[m, 0], P[m, 1])), shape=(n, n))
    return connected_components(A, directed=False)[1]


def duplicated_insertions(lib, P, d, se):
    """Pairs whose divergence is significantly below both copies' own LTR divergence: one insertion duplicated
    later (inherited by two homoeologues: orthologues; or copied by an exchange or segmental duplication)."""
    ks = _se(lib)
    a, b = P[:, 0], P[:, 1]
    with np.errstate(invalid="ignore"):
        return (d + SE_MULT * se < lib.k2p[a] - SE_MULT * ks[a]) & (d + SE_MULT * se < lib.k2p[b] - SE_MULT * ks[b])


def split_orthologs(d, orth, clusters, rng):
    """Orthologous-insertion split: median divergence of orthologous pairs, lineage-cluster bootstrap interval."""
    v = d[orth]
    if v.size < N_ORTH:
        return None
    bs = [np.median(v[sel]) for sel in cluster_bootstrap(clusters[orth], BOOT, rng)]
    est = float(np.median(v))
    lo, hi = widen((est, *np.quantile(bs, CI95)), MODEL_ERR["orthologs"])
    return dict(estimate=est, lo=lo, hi=hi, n=int(v.size))


def split_profile(d, cross, clusters, rng):
    """Lineage-floor split (profile half-rise), lineage-cluster bootstrap interval."""
    est = half_rise(d, cross)
    if not np.isfinite(est):
        return None
    bs = [half_rise(d[sel], cross[sel]) for sel in cluster_bootstrap(clusters, BOOT, rng)]
    bs = np.array([x for x in bs if np.isfinite(x)])
    if bs.size < BOOT_MIN:
        return None
    lo, hi = widen((est, *np.quantile(bs, CI95)), MODEL_ERR["profile"])
    return dict(estimate=est, lo=lo, hi=hi, n=int(d.size))


def _clock_bins(a, se, y, m, edges):
    """Per K2P bin: median K2P, weighted mean share y, total weight, median s.e."""
    bi = np.clip(np.searchsorted(edges, a, side="right") - 1, 0, edges.size - 2)
    nb = edges.size - 1
    W = np.bincount(bi, weights=m, minlength=nb)
    Y = np.bincount(bi, weights=m * y, minlength=nb) / np.maximum(W, 1e-300)
    A = np.array([np.median(a[bi == k]) if (bi == k).any() else np.nan for k in range(nb)])
    S = np.array([np.median(se[bi == k]) if (bi == k).any() else np.nan for k in range(nb)])
    ok = W > 0
    return A[ok], Y[ok], W[ok], S[ok]


def _one_step(A, Y, W, S, grid, young_high):
    """Weighted least-squares step between two levels blurred by K2P s.e."""
    from scipy.stats import norm
    w = norm.cdf((grid[:, None] - A[None, :]) / np.maximum(S, CONF_SE_FLOOR)[None, :])
    X1, X2 = w, 1 - w
    a11, a12, a22 = (W * X1 * X1).sum(1), (W * X1 * X2).sum(1), (W * X2 * X2).sum(1)
    b1, b2 = (W * X1 * Y).sum(1), (W * X2 * Y).sum(1)
    det = a11 * a22 - a12 ** 2
    det = np.where(np.abs(det) < 1e-12, np.nan, det)
    ly = np.clip((b1 * a22 - b2 * a12) / det, 0, 1)
    lo = np.clip((b2 * a11 - b1 * a12) / det, 0, 1)
    sse = (W[None, :] * (Y[None, :] - X1 * ly[:, None] - X2 * lo[:, None]) ** 2).sum(1)
    sse[~np.isfinite(sse) | ((ly < lo) if young_high else (ly > lo))] = np.inf
    if not np.isfinite(sse).any():
        return None
    j = int(np.argmin(sse))
    return float(grid[j]), float(ly[j]), float(lo[j])


def _divergence_step(A, Y, W, S, grid):
    k = np.convolve(W * Y, np.ones(TROUGH_POOL), "same") / np.maximum(np.convolve(W, np.ones(TROUGH_POOL), "same"),
                                                                   1e-300)
    t = int(np.argmin(k))
    go = grid[grid >= A[t]]
    return _one_step(A[t:], Y[t:], W[t:], S[t:], go, False) if A.size - t >= 2 and go.size else None


def split_confinement(k2p, se, share, m, rng, n_bins=CONF_BINS):
    """Confinement split: beyond the confined trough, the share of each copy's lineage evidence pointing to the other
    subgenome rises again (lineages present in both progenitors). Step fit on LTR K2P bins; bootstrap over copies.
    Reported only when the step is supported (rise > 0 in >= 97.5% of resamples)."""
    ok = np.isfinite(k2p) & np.isfinite(se) & (m > 0)
    if ok.sum() < CONF_MIN_COPIES:
        return None
    a, s, y, w = k2p[ok], np.maximum(se[ok], CONF_SE_FLOOR), share[ok], m[ok]
    hi = float(np.quantile(a, CONF_TOP_QUANTILE))
    q = np.unique(np.quantile(a[a > 0], np.linspace(0, 1, n_bins))) if (a > 0).any() else np.array([])
    edges = np.unique(np.concatenate([[0.0, 1e-12], q[(q > 1e-12) & (q < hi)], [hi, np.inf]]))
    grid = np.unique(np.concatenate([[0.0], np.linspace(0, hi, CONF_GRID_POINTS), edges[1:-1]]))
    fit = lambda idx: _divergence_step(*_clock_bins(a[idx], s[idx], y[idx], w[idx], edges), grid)
    dv = fit(np.arange(a.size))
    if dv is None:
        return None
    R = [fit(rng.integers(0, a.size, a.size)) for _ in range(BOOT)]
    t = np.array([r[0] for r in R if r is not None])
    rise = np.array([r[2] - r[1] for r in R if r is not None])
    if t.size < BOOT_MIN or np.quantile(rise, CI95[0]) <= 1e-9:             # rise > 0 beyond rounding
        return None
    lo, hi_ = widen((dv[0], *np.quantile(t, CI95)), MODEL_ERR["confinement"])
    return dict(estimate=dv[0], lo=lo, hi=hi_, n=int(ok.sum()))


def split_inherited(lib, ho, z, x, y, rng, threads, floor=0.0):
    """Split from inherited insertions. Anchors of the homoeolog step (reciprocal closest relatives among old
    copies; run_homoeologs, read only) that lie in a collinear chain on a called homoeolog pair of subgenomes x and y
    are insertions the progenitors' common ancestor carried: since the split each pair's two copies diverged
    independently, so their LTR-LTR divergence (5'-5' + 3'-3' alignments, >= 200 sites) is the split on the LTR clock
    of the confinement estimator and the merger clocks. An inherited insertion predates the split, and the split the
    merger: anchors no more diverged than the post-merger regime (floor = the cross-subgenome profile's young cut on the
    LTR clock; homoeologous exchanges and post-merger relatives) are dropped, then pairs homogenized after the split
    (flag_homogenized, no merger limit). Estimate = median; interval = bootstrap over homoeolog chromosome pairs
    (regions differ in rate), widened by MODEL_ERR['inherited']. None with fewer than N_ORTH usable anchors (old
    splits, where true inherited insertions are gone and exchange duplicates would dominate)."""
    if not ho or not len(ho.get("P", [])):
        return None
    P, si, sj = ho["P"], np.asarray(ho["si"]), np.asarray(ho["sj"])
    zi, zj = z[si], z[sj]
    hot = ho["in_chain"] & ho["pair_called"] & (((zi == x) & (zj == y)) | ((zi == y) & (zj == x)))
    idx = np.flatnonzero(hot)
    if idx.size < N_ORTH:
        return None
    cnt = _counts(lib, P[idx], "ltr", threads)
    k2p, se = divergence.k2p(cnt[:, 0], cnt[:, 1], cnt[:, 2])
    ks = _se(lib)
    a_, b_ = P[idx, 0], P[idx, 1]
    with np.errstate(invalid="ignore"):                           # each copy's own age, its lower limit (2 s.e.)
        older = np.fmin(lib.k2p[a_] - SE_MULT * ks[a_], lib.k2p[b_] - SE_MULT * ks[b_])
    C = int(max(si.max(), sj.max())) + 1
    blk = (np.minimum(si, sj) * C + np.maximum(si, sj))[idx]
    return inherited_estimate(k2p, cnt[:, 0], blk, rng, floor, se=se, older=older)


def inherited_estimate(k2p, sites, blocks, rng, floor=0.0, se=None, older=None):
    """The inherited-insertion split from anchor divergences (K2P), aligned sites and block labels (homoeolog
    chromosome pair): >= 200 sites, above the post-merger floor, one insertion duplicated by the split (pair divergence
    + 2 s.e. below both copies' own LTR divergence - 2 s.e., 'older'; independent insertions of one family and copies
    inserted after the split fail it), homogenized pairs dropped (flag_homogenized, no merger limit), median, block
    bootstrap widened by MODEL_ERR['inherited']; None with fewer than N_ORTH usable anchors."""
    k2p = np.asarray(k2p, float)
    ok = (np.asarray(sites) >= SG_LTR_MIN_SITES) & (k2p > floor)
    if older is not None:
        with np.errstate(invalid="ignore"):
            ok &= (k2p + SE_MULT * np.nan_to_num(np.asarray(se, float), nan=0.0)) < np.asarray(older, float)
    k2p = np.where(ok, k2p, np.nan)
    if np.isfinite(k2p).sum() < N_ORTH:
        return None
    flags, _ = flag_homogenized(k2p, sites, 0.0)
    use = np.isfinite(k2p) & ~flags
    if use.sum() < N_ORTH:
        return None
    v, b = k2p[use], np.asarray(blocks)[use]
    nb = int(np.unique(b).size)
    est = float(np.median(v))
    if nb >= 2:
        bs = [np.median(v[sel]) for sel in cluster_bootstrap(b, BOOT, rng)]
        lo, hi = widen((est, *np.quantile(bs, CI95)), MODEL_ERR["inherited"])
    else:
        lo, hi = widen((est, est, est), MODEL_ERR["inherited"])
    return dict(estimate=est, lo=lo, hi=hi, n=int(use.sum()), n_blocks=nb, n_flagged=int(flags.sum()))


def date_pair(lib, pt, loc, x, y, ancestry, rng, sizes=None, ptm=None, units=None, rng_m=None, threads=1, ho=None,
              z=None, rng_i=None):
    """Split and merger of subgenomes x and y (group indices) from the copy-pair table pt and copy locations loc.
    K2P units (pair divergence = 2 mu T; LTR divergence of a copy = 2 mu age). The merger uses the extended table ptm
    (pt's pairs + extra young focal copies; see merger_pair); the split estimators use pt.
    Every split estimate is on the LTR clock (orthologs/profile rescaled by the measured LTR/whole-record ratio, the
    genome-wide ratio when too few pairs lie near the estimate). Reported: the inherited-insertion estimator (ho = the
    homoeolog step's output, z = sequence labels; own stream rng_i) when enough anchors survive, else the log-scale
    median of the indirect estimators (consensus of profile, confinement and orthologs); the merger clocks are capped by
    the reported split on the same clock."""
    P, d, se = pt["P"], pt["d"], pt["se"]
    la, lb = loc[P[:, 0]], loc[P[:, 1]]
    m0 = np.isin(la, [x, y]) & np.isin(lb, [x, y])
    out = dict(groups=(x, y), split={}, merger=None, n_pairs=int(m0.sum()))
    if m0.sum() < DATE_MIN_PAIRS:
        return out
    Pm, dm, sm = P[m0], d[m0], se[m0]
    cross = (la[m0] != lb[m0]).astype(np.float64)
    clusters = lineage_clusters(Pm, dm, lib.n, CLUSTER_D)[Pm[:, 0]]
    dup = duplicated_insertions(lib, Pm, dm, sm)
    orth = (cross > 0) & dup & (dm > young_cut(dm, cross))
    out["split"]["orthologs"] = split_orthologs(dm, orth, clusters, rng)
    out["split"]["profile"] = split_profile(dm, cross, clusters, rng)
    lP, ninf = ancestry
    inside = np.isin(loc, [x, y]) & (ninf > 0)
    px, py = np.exp(lP[:, x]), np.exp(lP[:, y])
    share = np.where(loc == x, py, px) / np.maximum(px + py, 1e-300)
    out["split"]["confinement"] = split_confinement(lib.k2p, _se(lib), share, inside.astype(float), rng)
    cr = cross > 0                               # orthologs/profile: whole-record pair K2P -> LTR K2P (confinement's)
    g_sc = None
    for k in ("orthologs", "profile"):
        v = out["split"][k]
        if not v:
            continue
        sc = ltr_scale_near(lib, Pm[cr], dm[cr], pt["cnt"][m0][cr], v["estimate"], rng_m or rng, threads)
        if not sc:                                # too few pairs near the estimate -> the genome-wide ratio
            if g_sc is None:
                g_sc = ltr_rate_ratio(lib, Pm, pt["cnt"][m0], rng_m or rng, threads)
            sc = g_sc
        out["split"][k] = dict(v, estimate=v["estimate"] * sc, lo=v["lo"] * sc, hi=v["hi"] * sc, ltr_scale=sc)
    if ho is not None and z is not None:                          # anchors must predate the post-merger regime
        yc = young_cut(dm, cross)
        sc_y = ltr_scale_near(lib, Pm, dm, pt["cnt"][m0], yc, rng_m or rng, threads) if yc > 0 else None
        if not sc_y:
            sc_y = ltr_rate_ratio(lib, Pm, pt["cnt"][m0], rng_m or rng, threads)
        out["split"]["inherited"] = split_inherited(lib, ho, np.asarray(z), x, y, rng_i or rng, threads,
                                                    floor=yc * sc_y)
        out["inherited_floor"] = float(yc * sc_y)
    else:
        out["split"]["inherited"] = None
    comp = [k for k in ("profile", "confinement", "orthologs") if out["split"].get(k)]   # the indirect estimators
    if comp:
        lg = lambda key: float(np.exp(np.median([np.log(max(out["split"][k][key], 1e-12)) for k in comp])))
        out["split"]["consensus"] = dict(estimate=lg("estimate"), lo=lg("lo"), hi=lg("hi"), components=comp,
                                         n=int(sum(out["split"][k]["n"] for k in comp)))
    else:
        out["split"]["consensus"] = None
    rep = next((k for k in ("inherited", "consensus") if out["split"].get(k)), None)
    out["split_reported"] = rep
    ests = [v for k, v in out["split"].items() if v and k != "consensus" and np.isfinite(v["lo"])]
    out["split_agreement"] = ("none" if not ests else "single" if len(ests) == 1 else
                              "agree" if max(v["lo"] for v in ests) <= min(v["hi"] for v in ests) else "disagree")
    cap = out["split"][rep]["estimate"] if rep else None              # the reported split, LTR clock
    out.update(merger_pair(lib, ptm or pt, loc, x, y, rng if rng_m is None else rng_m, cap, sizes, units, threads,
                           plateau=plateau_end(dm, cross)))
    return out


def merger_pair(lib, pt, loc, x, y, rng, split, sizes, units, threads=1, plateau=None):
    """Merger of subgenomes x and y from sibling copy pairs (MRCA = the older copy's insertion; no duplicated
    insertions) younger than the confined plateau's end and the reported split. Two clocks on the same
    pairs and bootstrap units: the LTR clock, the younger copy's age from its LTR counts; the genealogy clock, the
    pair's divergence over the whole records (one substitution per aligned pair instead of per LTR) read through the
    genome's sibling genealogy f(r) learned from pairs with informative LTRs (sibling_genealogy). Reported:
    merger_report."""
    P, d, se = pt["P"], pt["d"], pt["se"]
    la, lb = loc[P[:, 0]], loc[P[:, 1]]
    m0 = np.isin(la, [x, y]) & np.isin(lb, [x, y])
    out = dict(merger=None, merger_ltr=None, merger_genealogy=None, n_sibling_pairs=0)
    if m0.sum() < DATE_MIN_PAIRS:
        return out
    Pm, dm, sm, cm = P[m0], d[m0], se[m0], pt["cnt"][m0]
    cross = (la[m0] != lb[m0]).astype(np.float64)
    comp, in7 = units
    clusters = np.where(in7[Pm[:, 0]], comp[Pm[:, 0]], comp[Pm[:, 1]])
    dup = duplicated_insertions(lib, Pm, dm, sm)
    ks = _se(lib)
    ka, kb = lib.k2p[Pm[:, 0]], lib.k2p[Pm[:, 1]]
    old = np.where(ka >= kb, Pm[:, 0], Pm[:, 1])
    with np.errstate(invalid="ignore"):
        sib = dm <= np.fmax(ka, kb) + SE_MULT * np.sqrt(ks[old] ** 2 + sm ** 2)
    kc, nc = ltr_counts(lib)
    pre = ~dup & sib & np.isfinite(kc[Pm[:, 0]]) & np.isfinite(kc[Pm[:, 1]]) & (nc[Pm[:, 0]] > 0) & (nc[Pm[:, 1]] > 0)
    if pre.sum() < SIB_MIN:
        return out
    scale = ltr_rate_ratio(lib, Pm[pre], cm[pre], rng, threads)     # measured first, so both caps share a clock
    t_pair = plateau_end(dm, cross) if plateau is None else plateau   # whole-record pair divergence
    if split:
        t_pair = min(t_pair, split / scale)                        # the split arrives on the LTR clock
    t_max = t_pair * scale                                         # the clocks' cap, LTR divergence
    good = pre & (dm < t_pair)
    out["n_sibling_pairs"] = int(good.sum())
    pool = np.flatnonzero(np.isin(loc, [x, y]) & np.isfinite(kc) & (nc > 0))
    if good.sum() < SIB_MIN or pool.size < SIB_MIN:
        return out
    pool = rng.choice(pool, min(pool.size, PRIOR_POOL_MAX), replace=False)
    g = age_prior(kc[pool], nc[pool])
    u, v = Pm[good, 0], Pm[good, 1]
    phi_land = None
    if sizes is not None:                              # random landing: the other subgenome's genome share
        og = loc[old[good]]
        phi_land = np.where(og == x, sizes[y], sizes[x]) / (sizes[x] + sizes[y])
    y_ = cross[good] > 0
    mg = merger_clock(younger_cdf(kc[u], nc[u], kc[v], nc[v], g), y_, clusters[good], rng, t_max,
                      1.0 / np.median(np.concatenate([nc[u], nc[v]])), phi_land=phi_land)
    kD, nD = (cm[good, 1] + cm[good, 2]).astype(float), cm[good, 0].astype(float)
    p_ = np.clip(kD / np.maximum(nD, 1.0), 0, JC_P_CAP)
    # JC69 on arrays (log1p) beside seq.divergence.jc69_from_identity's scalar form (ln): the two differ in the last
    # bits, so each caller keeps its own form and its outputs
    dlt = -0.75 * np.log1p(-4 * p_ / 3) * scale                  # pair divergence, LTR units
    teach = (dlt >= GEN_MIN_SUBS / np.median(np.minimum(nc[u], nc[v]))) & (dlt < t_max)
    f = (sibling_genealogy(kc[u][teach], nc[u][teach], kc[v][teach], nc[v][teach], dlt[teach])
         if teach.sum() >= GEN_MIN_PAIRS else np.full(R_GRID.size, 1.0 / R_GRID.size))
    mp = merger_clock(younger_cdf_genealogy(kc[u], nc[u], kc[v], nc[v], kD, nD, scale, f), y_, clusters[good], rng,
                      t_max, scale / np.median(nD))
    out.update(ltr_rate_ratio=scale, genealogy=f, n_teach=int(teach.sum()))
    out["merger_ltr"], out["merger_genealogy"] = mg, mp
    out["merger"] = merger_report(mg, mp, 1.0 / float(np.median(nc[np.isfinite(nc) & (nc > 0)])))
    return out


def merger_report(mg, mp, t_ltr):
    """The reported merger of a subgenome pair, by the order fixed on simulations: the genealogy clock's
    two-sided interval; else the LTR clock's; else the tighter of their bounds (the genealogy clock's upper limit is
    usually far tighter: one substitution per aligned pair instead of per LTR); else 'not datable' (the LTR clock's
    record, with its supplementary random-landing bound). A two-sided interval whose lower limit lies below t_ltr (one
    substitution per LTR) is reported as a bound [0, upper limit]: below the LTRs' own resolution a lower limit rests on
    the learned genealogy alone, and simulations put those too old (TE-poor genomes, hexaploids, lineages shared before
    the merger). Each dict gains 'source' (ltr_clock / genealogy_clock)."""
    if mp and mp["datable"] and not mp["bound"]:
        rep = dict(mp, source="genealogy_clock", landing_bound=None)
    elif mg and mg["datable"] and not mg["bound"]:
        rep = dict(mg, source="ltr_clock")
    else:
        bounds = [dict(m, source=k) for k, m in (("genealogy_clock", mp), ("ltr_clock", mg)) if m and m["datable"]]
        if not bounds:
            return dict(mg, source="ltr_clock") if mg else None
        rep = min(bounds, key=lambda m: m["hi"])
        if rep["source"] == "genealogy_clock":
            rep["landing_bound"] = None
    if not rep["bound"] and rep["lo"] < t_ltr:
        rep = dict(rep, bound=True, lo=0.0)
    return rep


def split_text(r, mu):
    """Log text of a subgenome pair's split (years)."""
    sp = r["split"].get(r.get("split_reported")) if r.get("split_reported") else None
    if not sp:
        return "split not datable"
    yr = lambda k: k / (2 * mu) / 1e6
    note = {"single": "one estimator", "agree": "estimators agree", "disagree": "estimators disagree"}
    return (f"split {yr(sp['estimate']):.2f} Myr ({yr(sp['lo']):.2f}-{yr(sp['hi']):.2f}; {r['split_reported']}, "
            f"{note.get(r.get('split_agreement'), '')})")


def merger_text(mg, n_sib, mu):
    """Log text of the merger clock: a date, a bound on the onset of cross-subgenome transposition, or why not."""
    kyr = lambda k: f"{k / (2 * mu) / 1e3:,.1f} kyr"
    if mg is None:
        return f"merger not datable: too few sibling copy pairs ({n_sib})"
    if not mg["datable"]:
        return ("merger not datable: no detectable excess of cross-subgenome sibling pairs (no copies moved between "
                "these subgenomes since the merger, or too few to tell)"
                + (f"; if new copies land in either subgenome at random, the merger is younger than "
                   f"{kyr(mg['landing_bound'])} (wrong where new copies stay in their own subgenome)"
                   if mg.get("landing_bound") else ""))
    src = {"ltr_clock": "LTR clock", "genealogy_clock": "genealogy clock"}.get(mg.get("source"), "clock")
    if mg["bound"]:
        return (f"cross-subgenome transposition began < {kyr(mg['hi'])} ago ({src}; a bound on its onset, not on the "
                f"merger)")
    return (f"merger {kyr(mg['estimate'])} ({kyr(mg['lo'])}-{kyr(mg['hi'])}; {src}, onset of cross-subgenome "
            f"transposition)")


def run_dating(a, lib, ph):
    """Dating of every pair of subgenomes (N >= 2): copy pairs, three split estimators, the merger clock."""
    if ph["N"] < 2:
        return None
    rng = np.random.default_rng([a.seed, 7])
    ek = ph["ek"]
    si = {s_: i for i, s_ in enumerate(lib.seqnames)}
    sizes = np.bincount(ph["z"], weights=[float(lib.seq_len[si[s_]]) for s_ in ph["names"]], minlength=ph["N"])
    amb = ~np.asarray(ph.get("core", np.ones(ph["C"], bool))) & ~(np.asarray(ph["support"]) >= SUPPORT_OK)
    loc = np.where((ek >= 0) & ~amb[np.maximum(ek, 0)], ph["z"][np.maximum(ek, 0)], -1)   # no unsure satellites
    pt = sample_pairs(lib, ek, rng, a.threads)
    if pt is None or len(pt["P"]) < DATE_MIN_PAIRS:
        log("dating skipped: too few copy pairs")
        return None
    vlog(f"dating: {len(pt['P']):,} copy pairs aligned")
    ancestry = copy_ancestry(ph["clades"], ek, ph["z"], ph["N"], rng)
    rng_m = np.random.default_rng([a.seed, 80])               # the merger table's own stream: pt's draws unchanged
    ptm = sample_young_more(lib, ek, pt, rng_m, a.threads)
    pt.pop("sketch", None)
    vlog(f"dating: {len(ptm['P']) - len(pt['P']):,} more copy pairs for the merger clock (young focal copies)")
    units = pair_units(pt, lib.n)
    pa = run_painting(a, lib, ph, dict(ancestry=ancestry))     # painted first, on the same ancestry
    res = []
    for x in range(ph["N"]):
        for y in range(x + 1, ph["N"]):
            res.append(date_pair(lib, pt, loc, x, y, ancestry, rng, sizes=sizes, ptm=ptm, units=units, rng_m=rng_m,
                                 threads=a.threads, ho=ph.get("homoeologs"), z=ph["z"],
                                 rng_i=np.random.default_rng([a.seed, 10, x, y])))
    mu = a.mu or MU_DEFAULT
    for r in res:
        x, y = (group_label(g) for g in r["groups"])
        log(f"{x}-{y}: {split_text(r, mu)}; {merger_text(r['merger'], r.get('n_sibling_pairs', 0), mu)} "
            f"[mu={mu:g}{'' if a.mu else ', default'}]")
    return dict(mu=mu, mu_given=bool(a.mu), pairs=res, n_copy_pairs=int(len(pt["P"])), ancestry=ancestry, pt=pt,
                n_copy_pairs_merger=int(len(ptm["P"])), painting=pa)


def dating_summary(dt):
    """Compact dating record for summary.json (years; K2P in dating.tsv)."""
    yr = lambda k: None if k is None else k / (2 * dt["mu"])
    pairs = []
    for r in dt["pairs"]:
        sp = r["split"].get(r.get("split_reported")) if r.get("split_reported") else None
        mg = r["merger"]
        st = "none" if mg is None else ("not_datable" if not mg["datable"] else ("bound" if mg["bound"] else "dated"))
        pairs.append(dict(subgenomes=[group_label(g) for g in r["groups"]],
                          split=dict(estimator=r.get("split_reported"), years=yr(sp["estimate"]),
                                     low=yr(sp["lo"]), high=yr(sp["hi"]), agreement=r.get("split_agreement"))
                          if sp else None,
                          merger=dict(status=st, years=yr(mg["estimate"]) if st == "dated" else None,
                                      low=yr(mg["lo"]) if st == "dated" else (0.0 if st == "bound" else None),
                                      high=yr(mg["hi"]) if st in ("dated", "bound") else None,
                                      clock=mg.get("source") if mg else None,
                                      if_random_landing_below_years=yr(mg["landing_bound"]) if mg and
                                      mg.get("landing_bound") else None),
                          split_estimators={k: dict(k2p=v["estimate"], low=v["lo"], high=v["hi"], n=v.get("n"),
                                                    ltr_scale=v.get("ltr_scale"), components=v.get("components"))
                                            for k, v in r["split"].items() if v},
                          inherited_floor=r.get("inherited_floor"), ltr_rate_ratio=r.get("ltr_rate_ratio"),
                          genealogy_mean_r=float(np.dot(R_GRID, r["genealogy"])) if r.get("genealogy") is not None
                          else None, n_sibling_pairs=r.get("n_sibling_pairs")))
    return dict(mu=dt["mu"], mu_given=dt["mu_given"], n_copy_pairs=dt["n_copy_pairs"],
                n_copy_pairs_merger=dt.get("n_copy_pairs_merger"), pairs=pairs)


def write_dating(a, dt):
    od = Path(a.outdir)
    with open(od / "dating.tsv", "w") as f:
        f.write(f"# mu={dt['mu']:g} substitutions/site/year{'' if dt['mu_given'] else ' (default)'}; K2P = 2 mu years\n")
        f.write("subgenome_a\tsubgenome_b\tevent\testimator\treported\tstatus\tk2p\tk2p_low\tk2p_high\tyears\t"
                "years_low\tyears_high\tn\n")
        for r in dt["pairs"]:
            x, y = (group_label(g) for g in r["groups"])
            for k in ("inherited", "consensus", "orthologs", "confinement", "profile"):
                v = r["split"].get(k)
                st = "dated" if v else "none"
                row = [x, y, "split", k, "yes" if r.get("split_reported") == k else "no", st]
                vals = [v["estimate"], v["lo"], v["hi"]] if v else [None] * 3
                f.write("\t".join(row + [_f(z_) for z_ in vals] + [_f(z_ / (2 * dt["mu"]) if z_ is not None else None)
                                                                     for z_ in vals] + [str(v["n"]) if v else "0"]) + "\n")
            mg = r["merger"]
            src = mg.get("source") if mg else None
            for k, m_ in (("ltr_clock", r.get("merger_ltr")), ("genealogy_clock", r.get("merger_genealogy"))):
                if k == src:
                    m_ = mg                      # the reported row carries what is reported (a demoted interval: bound)
                st = "none" if m_ is None else ("not_datable" if not m_["datable"] else
                                                ("bound" if m_["bound"] else "dated"))
                vals = ([None if m_["bound"] else m_["estimate"], m_["lo"], m_["hi"]] if m_ and m_["datable"]
                        else [None] * 3)                   # a bound has no point estimate
                f.write("\t".join([x, y, "merger", k, "yes" if src == k else "no", st] + [_f(z_) for z_ in vals] +
                                   [_f(z_ / (2 * dt["mu"]) if z_ is not None else None) for z_ in vals] +
                                   [str(m_["n"]) if m_ else "0"]) + "\n")
            mg = r.get("merger_ltr")
            lb = mg.get("landing_bound") if mg else None             # supplementary: assumes random landing
            vals = [None, 0.0, lb] if lb else [None] * 3
            f.write("\t".join([x, y, "merger_if_random", "landing", "no", "bound" if lb else "none"]
                               + [_f(z_) for z_ in vals] + [_f(z_ / (2 * dt["mu"]) if z_ is not None else None)
                                                            for z_ in vals] + [str(mg["n"]) if mg else "0"]) + "\n")
