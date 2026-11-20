"""The genotype model (routes.reads.model): read model, likelihoods, priors, calibration and calls."""
import math

import numpy as np
import pytest

from ltrtrace.routes.reads.model import (
    ReadModel,
    allele_posteriors,
    call,
    genotype_prior,
    genotypes,
    pav_call,
    site_calibration,
)


class _M:  # minimal read model for likelihood tests
    def __init__(self, d_hap=15.0, k=50.0, ploidy=2):
        self.d_hap, self.k, self.ploidy = d_hap, k, ploidy


RATES = {"F": 1.6, "S": 1.0, "E": 1.0}


@pytest.mark.parametrize("counts,expect", [
    ({"F": 48, "S": 0, "E": 0}, "F/F"),
    ({"F": 24, "S": 0, "E": 15}, "F/E"),
    ({"F": 24, "S": 15, "E": 0}, "F/S"),
    ({"F": 0, "S": 0, "E": 30}, "E/E"),
    ({"F": 24, "S": 0, "E": 0}, "F/N"),  # half the reads a two-copy site would give: one haplotype lacks it
    ({"F": 0, "S": 0, "E": 0}, "N/N"),
    ({"F": 47, "S": 0, "E": 1}, "F/F"),  # one stray read does not make a heterozygote
])
def test_model_calls(counts, expect):
    gt, gq, _lls, _mu = call(counts, RATES, _M(), 0.02)
    assert gt == expect
    assert 0 <= gq <= 99


def test_genotype_space():
    assert len(genotypes(2)) == 10 and "F/N" in genotypes(2)


def test_read_model_rates():
    m = ReadModel([30] * 100, [12000] * 200, 2)
    assert abs(m.depth - 30) < 2 and abs(m.d_hap - 15) < 1
    assert m.r(0) == pytest.approx(1.0) and m.r(12000) == 0.0 and m.r(6000) == pytest.approx(0.5)


def test_low_mapq_only_is_no_call():
    """No informative read but reads placed with low MAPQ: a no-call, not a confident N/N (mappability)."""
    gt, gq, _lls, _mu = call({"F": 0, "S": 0, "E": 0, "lowq": 6}, RATES, _M(), 0.02)
    assert gt == "./." and gq == 0
    gt, gq, _lls, _mu = call({"F": 0, "S": 0, "E": 0, "lowq": 0}, RATES, _M(), 0.02)
    assert gt == "N/N"  # no reads at all: the site is missing from both haplotypes


def test_phase_lower_bound_in_likelihood():
    """Two distinct haplotypes carrying the full element tip a borderline depth call to F/F (moderate weight:
    PH_ERR 0.05, so clear depth evidence still wins)."""
    c = {"F": 30, "S": 0, "E": 0}  # F/N expects 24, F/F 48 reads
    assert call(c, RATES, _M(), 0.02)[0] == "F/N"
    assert call(dict(c, hap_F=2), RATES, _M(), 0.02)[0] == "F/F"
    gq0 = call({"F": 20, "S": 0, "E": 0}, RATES, _M(), 0.02)[1]
    gq2 = call({"F": 20, "S": 0, "E": 0, "hap_F": 2}, RATES, _M(), 0.02)[1]
    assert gq2 < gq0  # conflicting phase evidence lowers confidence


def test_site_calibration_recovers_haplotype_depth():
    """Full-read counts at sites with an empty/solo allele fix the per-haplotype depth at element sites, whatever
    the genome-wide estimate (which hemizygous sequence and overdispersion bias)."""
    rng = np.random.default_rng(3)
    rates = {"F": 1.5, "S": 1.0, "E": 1.0}
    recs = []
    for _ in range(60):  # F/E sites: one full haplotype at d = 13
        recs.append(({"F": int(rng.poisson(13 * 1.5)), "S": 0, "E": int(rng.poisson(13))}, rates, "F/E", 60))
    for _ in range(200):  # F/F sites carry no calibration information
        recs.append(({"F": int(rng.poisson(26 * 1.5)), "S": 0, "E": 0}, rates, "F/F", 60))
    d, k, n = site_calibration(recs, _M(d_hap=9.0))
    assert n == 60 and abs(d - 13) < 1.0 and k >= 100  # Poisson counts: little overdispersion
    assert site_calibration(recs[:20], _M()) is None  # too few sites: keep the genome-wide depth


def test_genotype_prior_learns_frequencies():
    """EM recovers genotype frequencies from per-site likelihoods; ties then follow the common genotype."""
    gts = genotypes(2)
    common = {g: (0.0 if g == "F/F" else -8.0) for g in gts}  # 80 sites clearly F/F
    rare = {g: (0.0 if g == "F/N" else -8.0) for g in gts}  # 20 sites clearly F/N
    tie = {g: (0.0 if g in ("F/F", "F/N") else -8.0) for g in gts}  # a site that cannot tell them apart
    lp = genotype_prior([common] * 80 + [rare] * 20, gts)
    assert abs(math.exp(lp["F/F"]) - 81 / 110) < 0.02 and abs(math.exp(lp["F/N"]) - 21 / 110) < 0.02
    call({"F": 30, "S": 0, "E": 0}, RATES, _M(), 0.02, lp)
    assert max(tie, key=lambda g: tie[g] + lp[g]) == "F/F"
    assert genotype_prior([common] * 40 + [rare] * 9, gts) is None  # < 5 sites per class (50): flat prior


def test_solo_empty_is_solo():
    """One haplotype with a solo LTR and one with the empty site: one chromosome carries the insertion, so the
    presence/absence call is solo, not absent; empty sites alone stay absent."""
    gt, _gq, lls, _mu = call({"F": 0, "S": 15, "E": 15}, RATES, _M(), 0.02)
    pa = allele_posteriors(lls)
    assert gt == "S/E" and pa["S"] >= 0.99 and pa["E"] >= 0.99 and pa["F"] <= 0.01
    assert pav_call(pa) == "solo"
    assert pav_call(allele_posteriors(call({"F": 0, "S": 0, "E": 30}, RATES, _M(), 0.02)[2])) == "absent"
    assert pav_call(allele_posteriors(call({"F": 0, "S": 30, "E": 0}, RATES, _M(), 0.02)[2])) == "solo"
    assert pav_call({"F": 0.0, "E": 1.0, "S": 1.0}) == "solo" and pav_call({"F": 0.0, "E": 1.0, "S": 0.0}) == "absent"
    assert pav_call({"F": 0.995, "E": 1.0, "S": 1.0}) == "present" and pav_call({"F": 0.5, "E": 1.0, "S": 1.0}) == (
        "uncertain")
