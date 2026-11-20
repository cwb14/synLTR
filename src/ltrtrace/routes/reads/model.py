"""The genotype model of one read set.

Alleles: F full, S solo LTR, E empty site, N null (the haplotype gives no informative reads: the site and its
flanks are missing or rearranged on it). A haplotype carrying allele a yields informative reads at rate d_hap * r_a,
with r_a from the sample's read-span distribution (r_N = 0). Likelihood of genotype g:
  NB(total informative reads | d_hap * sum_a n_g(a) r_a, k)  x  multinomial(F/S/E counts | allele mixture, error eps)
Flat prior over genotypes unless the sample's own genotype frequencies can be learned (genotype_prior);
GQ = phred(1 - posterior of the best genotype).
"""
from __future__ import annotations

import math
import random
from itertools import combinations_with_replacement

import numpy as np

from ltrtrace.base.constants import ANCHOR as A
from ltrtrace.base.constants import (
    CAL_MIN_EMPTY,
    CAL_MIN_SITES,
    CONF_GQ,
    DEPTH_COLLAPSE,
    DEPTH_FINE,
    DEPTH_GRID,
    DEPTH_MIN_POINTS,
    DEPTH_MIN_READS,
    EB_ITERS,
    EB_MIN_PER_CLASS,
    EPS_CAP,
    EPS_FLOOR,
    EPS_MIN_READS,
    GQ_CAP,
    MIXTURE_ITERS,
    MU_MIN,
    N_DEPTH_POINTS,
    NB_SIZE_FINE,
    NB_SIZE_GRID,
    NB_SIZE_MAX,
    NB_SIZE_MIN,
    PAV_P,
    PH_ERR,
    PROFILE_READS,
    SEED,
)

CLASSES = "FSE"


def fit_copy_mixture(depths, ploidy: int) -> tuple[float, float, np.ndarray]:
    """Depth at random positions as a mixture of 1..ploidy copies of one haplotype depth d: copy number c has
    NB(c d, size c k) (a sum of c haplotypes). Returns (d, k, weights). Grid search over d and k (vectorized over d),
    EM for the weights, a mild Dirichlet prior (alpha 2 on the full-ploidy component) to break the tie a single peak
    leaves between "ploidy copies of d" and "fewer copies of a larger d"."""
    x = np.asarray(depths, float)
    med = float(np.median(x))
    x = np.round(x[x <= DEPTH_COLLAPSE * med])  # collapsed repeats are not copies of the genome
    xi = x.astype(np.int64)
    logfact = np.concatenate(([0.0], np.cumsum(np.log(np.arange(1, int(x.max()) + 1)))))
    cs = np.arange(1, ploidy + 1)
    d_grid = med / ploidy * np.exp(np.linspace(np.log(DEPTH_GRID[0]), np.log(DEPTH_GRID[1]), DEPTH_GRID[2]))
    k_grid = np.exp(np.linspace(np.log(NB_SIZE_GRID[0]), np.log(NB_SIZE_GRID[1]), NB_SIZE_GRID[2]))
    prior = np.where(cs == ploidy, 1.0, 0.0)
    best = _mixture_search(x, xi, logfact, cs, d_grid, k_grid, prior, ploidy)
    # refine around the best grid point (+-1 step in d and k)
    d0, k0 = best[1], best[2]
    d_fine = d0 * np.exp(np.linspace(*DEPTH_FINE))
    k_fine = k0 * np.exp(np.linspace(*NB_SIZE_FINE))
    fine = _mixture_search(x, xi, logfact, cs, d_fine, k_fine, prior, ploidy)
    best = fine if fine[0] > best[0] else best
    return best[1], best[2], best[3]


def _mixture_search(x, xi, logfact, cs, d_grid, k_grid, prior, ploidy):
    best = None
    for k in k_grid:
        L = np.empty((len(x), len(d_grid), ploidy))
        for j, c in enumerate(cs):
            kc = c * k
            lg = np.concatenate(([0.0], np.cumsum(np.log(kc + np.arange(int(x.max()))))))[xi]
            p = kc / (kc + c * d_grid)  # per d
            L[:, :, j] = (lg - logfact[xi])[:, None] + kc * np.log(p)[None, :] + x[:, None] * np.log1p(-p)[None, :]
        m = L.max(2, keepdims=True)
        P = np.exp(L - m)  # points x d x copies
        w = np.full((len(d_grid), ploidy), 1.0 / ploidy)
        for _ in range(MIXTURE_ITERS):
            R = P * w[None]
            R /= R.sum(2, keepdims=True)
            w = R.sum(0) + prior[None]
            w /= w.sum(1, keepdims=True)
        ll = (m[:, :, 0] + np.log((P * w[None]).sum(2))).sum(0) + np.log(w[:, -1])
        j = int(np.argmax(ll))
        if best is None or ll[j] > best[0]:
            best = (float(ll[j]), float(d_grid[j]), float(k), w[j].copy())
    return best


class ReadModel:
    """Per-sample depth and read-span model, estimated from random positions."""

    def __init__(self, depths, spans, ploidy: int):
        depths = np.asarray(depths, float)
        self.spans = np.sort(np.asarray(spans, float))
        self.ploidy = ploidy
        pos = depths[depths > 0]
        if len(pos) < DEPTH_MIN_POINTS or len(self.spans) < DEPTH_MIN_READS:
            raise ValueError(f"too little coverage to model depth ({len(pos)} covered points, {len(self.spans)} reads)")
        # the genome may be partly hemizygous (outcrossers: up to ~40% of sequence on one haplotype): fit depth as a
        # mixture of 1..ploidy copies of one haplotype depth instead of taking a single mode
        self.d_hap, k_hap, self.copy_weights = fit_copy_mixture(pos, ploidy)
        self.depth = self.d_hap * ploidy
        self.k = float(min(NB_SIZE_MAX, max(NB_SIZE_MIN, k_hap * ploidy)))  # NB size of all-haplotype depth
        self.mean_span = float(self.spans.mean())

    def r(self, x: float) -> float:
        """Expected reads covering an interval of length x per unit depth: E[max(0, l - x) / l] over a
        coverage-weighted (length-biased) sample of read spans."""
        return float(np.mean(np.clip(self.spans - x, 0, None) / self.spans))

    def rates(self, site) -> dict[str, float]:
        """Informative-read rates per allele for one element (reads per unit per-haplotype depth)."""
        el = site.el
        l5, l3 = el.lltr[1] - el.lltr[0], el.rltr[1] - el.rltr[0]
        rF = self.r(l5 + 2 * A) + self.r(l3 + 2 * A) - self.r(site.L + 2 * A)
        rS = self.r(el.ltr_mean + 2 * A + site.t)
        rE = self.r(2 * A + site.t)
        return {"F": max(rF, 1e-9), "S": max(rS, 1e-9), "E": max(rE, 1e-9)}


def sample_points(chrom_lengths: dict[str, int], n: int = N_DEPTH_POINTS, seed: int = SEED) -> list[tuple[str, int]]:
    """Random positions, uniform over the given sequences (dict name -> length)."""
    rng = random.Random(seed)
    names = list(chrom_lengths)
    tot = sum(chrom_lengths.values())
    pts = []
    for _ in range(n):
        u = rng.randrange(tot)
        for c in names:
            if u < chrom_lengths[c]:
                pts.append((c, u))
                break
            u -= chrom_lengths[c]
    return pts


def depth_and_spans(bam, points, min_mapq: int) -> tuple[list[int], list[int]]:
    """Depth (primary alignments, MAPQ >= min_mapq) at each point, and the reference spans of those alignments
    (a coverage-weighted sample of read spans)."""
    depths, spans = [], []
    for c, p in points:
        n = 0
        for x in bam.fetch(c, p, p + 1):
            if (x.is_unmapped or x.is_secondary or x.is_supplementary or x.is_duplicate or x.is_qcfail
                    or x.mapping_quality < min_mapq):  # the reads evidence() can use
                continue
            n += 1
            spans.append(x.reference_end - x.reference_start)
        depths.append(n)
    return depths, spans


def read_profile(bam, points, min_mapq: int, n_max: int = PROFILE_READS) -> tuple[float, float]:
    """Share of reads with base qualities and median read divergence (NM / aligned bases) over primary alignments at
    the sampled points: decides whether read phasing is safe for the sample (see phase.py)."""
    seen, qual, div = set(), 0, []
    for c, p in points:
        for x in bam.fetch(c, p, p + 1):
            if (x.is_unmapped or x.is_secondary or x.is_supplementary or x.is_duplicate or x.is_qcfail
                    or x.mapping_quality < min_mapq
                    or x.query_name in seen):
                continue
            seen.add(x.query_name)
            qual += x.query_qualities is not None
            if x.has_tag("NM") and x.query_alignment_length:
                div.append(x.get_tag("NM") / x.query_alignment_length)
            if len(seen) >= n_max:
                break
        if len(seen) >= n_max:
            break
    if not seen:
        return 0.0, float("nan")
    return qual / len(seen), (float(np.median(div)) if div else float("nan"))


def genotypes(ploidy: int) -> list[str]:
    return ["/".join(g) for g in combinations_with_replacement("FSEN", ploidy)]


def _nb_logpmf(x: float, mu: float, k: float) -> float:
    p = k / (k + mu)
    return (math.lgamma(x + k) - math.lgamma(k) - math.lgamma(x + 1) + k * math.log(p) + x * math.log1p(-p))


def loglik(counts: dict, gt: str, rates: dict[str, float], d_hap: float, eps: float, k: float) -> float:
    alle = gt.split("/")
    lam = {a: alle.count(a) * rates[a] for a in CLASSES}
    tot = sum(lam.values())
    T = sum(counts[c] for c in CLASSES)
    ll = _nb_logpmf(T, max(d_hap * tot, MU_MIN), k)
    if T:
        for c in CLASSES:
            if counts[c]:
                p = sum((lam[a] / tot) * ((1 - eps) if c == a else eps / 2) for a in CLASSES) if tot > 0 else 1 / 3
                ll += counts[c] * math.log(p)
    # read phasing: h distinct haplotypes seen carrying allele a need >= h copies of a (each missing copy = a
    # spurious haplotype group, probability PH_ERR)
    for a in CLASSES:
        short = counts.get("hap_" + a, 0) - alle.count(a)
        if short > 0:
            ll += short * math.log(PH_ERR)
    return ll


def allele_posteriors(lls: dict[str, float]) -> dict[str, float]:
    """Posterior probability that at least one haplotype carries F, E or S (sum over genotypes; lls are the
    log-posteriors returned by call(), so the sample's genotype prior is included)."""
    mx = max(lls.values())
    post = {g: math.exp(v - mx) for g, v in lls.items()}
    z = sum(post.values())
    return {a: sum(p for g, p in post.items() if a in g.split("/")) / z for a in CLASSES}


def pav_call(pa: dict[str, float]) -> str:
    """Presence/absence from allele posteriors: present (an F allele), solo (a solo LTR, no F), absent (an empty
    site, no F and no solo LTR), else uncertain. A solo LTR beside an empty site is solo: one chromosome carries the
    insertion."""
    if pa["F"] >= PAV_P:
        return "present"
    if pa["F"] <= 1 - PAV_P and pa["S"] >= PAV_P:
        return "solo"
    if pa["F"] <= 1 - PAV_P and pa["E"] >= PAV_P:
        return "absent"
    return "uncertain"


def genotype_prior(lls_list: list[dict[str, float]], gts: list[str], n_iter: int = EB_ITERS) -> dict[str, float] | None:
    """Empirical-Bayes genotype frequencies of one sample: EM over the per-site genotype log-likelihoods, with
    add-one (Dirichlet 1) smoothing. Returns {genotype: log prior}, or None (flat prior) with fewer than
    EB_MIN_PER_CLASS sites per genotype class, where the frequencies cannot be estimated. Near-ties then follow the
    genotypes common in this sample (e.g. F/N in an outcrosser, F/F in an inbred line) instead of a flat prior."""
    L = np.array([[ll[g] for g in gts] for ll in lls_list], float)
    if len(L) < EB_MIN_PER_CLASS * len(gts):
        return None
    L -= L.max(1, keepdims=True)
    P0 = np.exp(L)
    pi = np.full(len(gts), 1.0 / len(gts))
    for _ in range(n_iter):
        P = P0 * pi
        P /= P.sum(1, keepdims=True)
        pi = (P.sum(0) + 1.0) / (len(P) + len(gts))
    return {g: float(np.log(p)) for g, p in zip(gts, pi)}


def call(counts: dict, rates: dict[str, float], model, eps: float,
         logprior: dict[str, float] | None = None) -> tuple[str, int, dict[str, float], float]:
    """Return (best genotype, GQ, dict of log-posteriors up to a constant (log-likelihood + log prior; flat prior if
    logprior is None), expected informative reads under the best genotype).
    No informative read but reads placed here with low MAPQ: no call ('./.', GQ 0, flat likelihoods); the
    absence of evidence there reflects mappability, not the alleles."""
    gts = genotypes(model.ploidy)
    if sum(counts[c] for c in CLASSES) == 0 and counts.get("lowq", 0) > 0:
        return "/".join("." * model.ploidy), 0, {g: 0.0 for g in gts}, 0.0
    lls = {g: loglik(counts, g, rates, model.d_hap, eps, model.k) + (logprior[g] if logprior else 0.0)
           for g in gts}
    mx = max(lls.values())
    post = {g: math.exp(v - mx) for g, v in lls.items()}
    z = sum(post.values())
    best = max(post, key=post.get)
    pb = post[best] / z
    gq = GQ_CAP if pb >= 1 - 10 ** (-GQ_CAP / 10) else round(-10 * math.log10(1 - pb))
    alle = best.split("/")
    mu = model.d_hap * sum(alle.count(a) * rates[a] for a in CLASSES)
    return best, min(gq, GQ_CAP), lls, mu


def site_calibration(records, model) -> tuple[float, float, int] | None:
    """Per-haplotype depth and overdispersion of informative reads, learned from the sample's own sites whose full-
    allele dosage the reads fix without depth: sites called with an empty or solo allele beside >= 1 full allele (GQ
    >= CONF_GQ, >= 2 empty/solo reads), where the full reads come from the called number of full haplotypes.
    Genome-wide depth misleads where much of the genome is on one haplotype only (outcrossers) or depth is very
    overdispersed (ONT). records: (counts, rates, genotype, gq). Returns (d_hap, k, n_sites) or None (< CAL_MIN_SITES
    sites)."""
    obs, base = [], []
    for counts, rates, gt, gq in records:
        al = gt.split("/")
        nf = al.count("F")
        if gq < CONF_GQ or nf == 0 or not ({"E", "S"} & set(al)) or counts["E"] + counts["S"] < CAL_MIN_EMPTY:
            continue
        obs.append(counts["F"])
        base.append(nf * rates["F"])
    if len(obs) < CAL_MIN_SITES:
        return None
    obs, base = np.asarray(obs, float), np.asarray(base, float)
    d_hap = float(np.median(obs / base))
    mu = d_hap * base
    inv_k = float(np.mean(((obs - mu) ** 2 - mu) / mu ** 2))
    k = NB_SIZE_MAX if inv_k <= 1 / NB_SIZE_MAX else float(min(NB_SIZE_MAX, max(NB_SIZE_MIN, 1 / inv_k)))
    return d_hap, k, len(obs)


def estimate_eps(records) -> float | None:
    """Misclassification rate from confident homozygous calls: discordant / total informative reads.
    records: iterable of (counts, genotype, gq)."""
    disc = tot = 0
    for counts, gt, gq in records:
        alle = set(gt.split("/")) - {"N"}
        if "." in alle:
            continue
        T = sum(counts[c] for c in CLASSES)
        if len(alle) != 1 or T < EPS_MIN_READS or gq < CONF_GQ:
            continue
        (a,) = alle
        disc += T - counts[a]
        tot += T
    if tot == 0:
        return None
    return min(EPS_CAP, max(EPS_FLOOR, disc / tot))
