"""The number of subgenomes: split acceptance on pooled copies and age classes, the divisive search, core
sequences and satellites."""
from types import SimpleNamespace

import numpy as np
import pytest

from synltr.analyses.subgenomes import halves, partition
from synltr.base import constants as cst


def _c(src, rep, con):
    return dict(src=src, rep=rep, con=con, z=np.array([0, 0, 1, 1]), se=0.0, n_pairs=50, borderline=False, hz=[])


def _parent(con):
    return dict(con=dict(zip(("all",) + cst.AGE_STRATA, con)))


def test_pooled_split_decides_first():
    c, how = halves.accept_split({"all": _c("all", 1.0, 0.6), "young": _c("young", 1.0, 0.1)}, _parent((1.0, 1, 1, 1)))
    assert how == "pooled" and c["src"] == "all"


def test_age_bounded_split_is_kept_on_its_tercile():
    """Weak on pooled copies (< half its parent), strong in young copies, absent in old ones: kept on young."""
    cands = {"all": _c("all", 1.0, 0.3), "young": _c("young", 1.0, 0.9), "middle": _c("middle", 1.0, 0.1),
             "old": _c("old", 0.2, 0.1)}
    c, how = halves.accept_split(cands, _parent((1.6, 1.0, 2.0, 1.7)))
    assert how == "age" and c["src"] == "young"


def test_split_present_at_every_age_is_rejected():
    """Reproducible in every tercile but weaker than half its parent pooled: the chromosomes' own structure."""
    cands = {"all": _c("all", 1.0, 0.56), "young": _c("young", 0.9, 0.3), "middle": _c("middle", 0.95, 0.66),
             "old": _c("old", 0.91, 0.93)}
    c, how = halves.accept_split(cands, _parent((1.19, 1.17, 1.47, 1.10)))
    assert c is None and how == "present at every age"


def test_no_rescue_at_the_top_or_against_an_absent_parent():
    cands = {"all": _c("all", 0.5, 0.05), "young": _c("young", 1.0, 0.5), "old": _c("old", 0.1, 0.0)}
    assert halves.accept_split(cands, None) == (None, "not reproduced")
    c, how = halves.accept_split(cands, _parent((1.0, 0.0, 1.0, 1.0)))          # parent invisible in young copies
    assert c is None


class _PlantedSources:
    """Sources stand-in: one planted matrix for every source, noisy half-library pairs, more() adding pairs."""

    def __init__(self, W, rng):
        self.names, self.rng, self.W0 = ["all"] + list(cst.AGE_STRATA), rng, W
        self.W = {s: W for s in self.names}
        self.pairs = {s: [self._pair() for _ in range(20)] for s in self.names}

    def _noisy(self):
        M = self.W0 + self.rng.random(self.W0.shape) * 0.3
        return (M + M.T) * (1 - np.eye(self.W0.shape[0]))

    def _pair(self):
        return self._noisy(), self._noisy()

    def more(self, s, n):
        new = [self._pair() for _ in range(n)]
        self.pairs[s].extend(new)
        return new


@pytest.mark.parametrize("n_force,groups,forced", [(None, 3, 0), (3, 3, 0), (5, 5, 2)])
def test_find_subgenomes_inferred_and_forced(n_force, groups, forced):
    """Planted ((B,C),A): inferred and -n 3 give the three blocks; -n 5 forces two splits the rule rejects."""
    rng = np.random.default_rng(0)
    z0 = np.repeat([0, 1, 2], 4)
    W = rng.random((12, 12)) * 0.2 + (z0[:, None] == z0[None, :]) * 2.0 + ((z0[:, None] > 0) & (z0[None, :] > 0)) * 1.5
    W = W + W.T
    np.fill_diagonal(W, 0)
    lab, tests, sup = halves.find_subgenomes(_PlantedSources(W, rng), rng, n_force=n_force)
    assert int(lab.max()) + 1 == groups and sum(t["forced"] for t in tests) == forced
    if groups == 3:
        assert partition.ari(lab, z0) == 1.0 and np.nanmin(sup) >= 0.9


def test_core_drops_scaffolds():
    n = np.array([1500, 1200, 1800, 900, 1100, 1300, 25, 12, 40, 15])
    core = halves.core_sequences(n)
    assert core.tolist() == [True] * 6 + [False] * 4


def test_core_all_when_no_break():
    n = np.array([300, 250, 280, 310, 260, 290])
    assert halves.core_sequences(n).all()


def test_core_falls_back_with_too_few():
    n = np.array([5000, 4000, 3000, 10, 12, 11])          # only 3 sequences pass -> everything is core
    assert halves.core_sequences(n).all()


def test_core_keeps_order_of_magnitude():
    n = np.array([1000] * 8 + [150, 99])                  # copy-N50 1000: 150 (>= 100) in, 99 out
    core = halves.core_sequences(n)
    assert core[8] and not core[9]


class _Src:
    def __init__(self, W, halves):
        self.W = {"all": W}
        self.pairs = {"all": halves}


def _block_matrix(rng, lab, strength=5.0, noise=1.0):
    C = lab.size
    W = rng.random((C, C)) * noise
    W = W + W.T
    same = lab[:, None] == lab[None, :]
    W[same] += strength
    np.fill_diagonal(W, 0)
    return W


def test_satellites_join_their_group():
    rng = np.random.default_rng(1)
    truth = np.array([0, 0, 0, 1, 1, 1, 0, 1])            # last two are satellites
    W = _block_matrix(rng, truth)
    W[6:, :] *= 0.05                                       # satellites: few links, same pattern
    W[:, 6:] *= 0.05
    halfs = [(W * rng.uniform(0.4, 0.6, W.shape), W * rng.uniform(0.4, 0.6, W.shape)) for _ in range(5)]
    halfs = [(0.5 * (a + a.T), 0.5 * (b + b.T)) for a, b in halfs]
    lab = np.array([0, 0, 0, 1, 1, 1, -1, -1])
    out, sup = halves.assign_satellites(_Src(W, halfs), lab, np.array([6, 7]))
    assert out.tolist() == truth.tolist()
    assert np.all(sup[6:] >= 0.9)
    assert np.all(np.isnan(sup[:6]))


def test_unlinked_satellite_has_zero_support():
    rng = np.random.default_rng(2)
    truth = np.array([0, 0, 1, 1, 0])
    W = _block_matrix(rng, truth)
    W[4, :] = 0
    W[:, 4] = 0
    lab = np.array([0, 0, 1, 1, -1])
    out, sup = halves.assign_satellites(_Src(W, [(W, W)]), lab, np.array([4]))
    assert sup[4] == 0.0
    assert out[4] in (0, 1)


def _cand(src, rep, con, z):
    return dict(src=src, rep=rep, con=con, z=np.asarray(z), se=0.01, hz=[], n_pairs=50, borderline=False)


Z = [0, 0, 0, 1, 1, 1]


ZX = [0, 1, 0, 1, 0, 1]


STRATA = ["young", "middle", "old"]


def test_top_pooled_pass():
    c = {"all": _cand("all", 0.95, 0.5, Z)}
    ch, how = halves.accept_split(c, None, STRATA)
    assert how == "pooled" and ch is c["all"]


def test_top_rescued_by_one_age_class():
    c = {"all": _cand("all", 0.6, 0.1, Z), "young": _cand("young", 0.1, 0.05, ZX),
         "middle": _cand("middle", 0.95, 0.4, Z), "old": _cand("old", 0.2, 0.05, ZX)}
    ch, how = halves.accept_split(c, None, STRATA)
    assert how == "age" and ch["src"] == "middle"


def test_age_rescue_needs_strict_replication():
    c = {"all": _cand("all", 0.6, 0.1, Z), "young": _cand("young", 0.1, 0.05, ZX),
         "middle": _cand("middle", 0.85, 0.4, Z), "old": _cand("old", 0.2, 0.05, ZX)}
    assert halves.accept_split(c, None, STRATA)[0] is None                 # top: 0.85 < MIN_REP_AGE
    parent = {"con": {"all": 1.0, "young": 0.3, "middle": 0.5, "old": 0.4}}
    assert halves.accept_split(c, parent)[0] is None                      # below the top too
    c["middle"] = _cand("middle", 0.97, 0.4, Z)
    ch, how = halves.accept_split(c, parent)
    assert how == "age" and ch["src"] == "middle"


def test_top_rescue_needs_agreement_with_pooled():
    c = {"all": _cand("all", 0.6, 0.1, ZX), "middle": _cand("middle", 0.95, 0.4, Z),
         "young": _cand("young", 0.1, 0.05, ZX), "old": _cand("old", 0.1, 0.05, ZX)}
    ch, how = halves.accept_split(c, None, STRATA)
    assert ch is None


def test_top_not_rescued_when_present_at_every_age():
    c = {"all": _cand("all", 0.6, 0.1, Z)}
    for s in STRATA:
        c[s] = _cand(s, 0.95, 0.3, Z)
    ch, how = halves.accept_split(c, None, STRATA)
    assert ch is None and how == "present at every age"


def test_top_rescue_needs_positive_strength():
    c = {"all": _cand("all", 0.6, 0.1, Z), "middle": _cand("middle", 0.95, -0.1, Z)}
    ch, how = halves.accept_split(c, None, STRATA)
    assert ch is None


def test_below_top_split_needs_half_its_parent():
    parent = {"con": {"all": 1.0, "young": 1.0, "middle": 1.0, "old": 1.0}}
    c = {"all": _cand("all", 0.95, 0.3, Z)}                 # 0.3 < 0.5 x 1.0: weaker than half its parent
    ch, how = halves.accept_split(c, parent)
    assert ch is None and how == "weaker than half its parent"


class _AgeSources:
    """Sources stand-in: 'all' carries A (0-3) vs the rest (4-11); inside the rest, 4-7 vs 8-11 shows only in one age
    stratum's copies (the other strata flat): a sub-split the age rule keeps on that stratum."""

    def __init__(self, rng, stratum="young"):
        self.rng = rng
        self.names = ["all"] + list(cst.AGE_STRATA)
        z = np.repeat([0, 1, 2], 4)
        top = (z[:, None] == 0) == (z[None, :] == 0)
        inner = (z[:, None] == z[None, :]) & (z[:, None] > 0)
        self.W0 = {s: top * 2.0 for s in cst.AGE_STRATA}
        self.W0[stratum] = top * 1.0 + inner * 3.0
        self.W0["all"] = top * 2.0 + inner * 0.3
        self.W = {s: self._sym(self.W0[s], 0.0) for s in self.names}
        self.pairs = {s: [self._pair(s) for _ in range(30)] for s in self.names}

    def _sym(self, M, noise):
        M = M + self.rng.random(M.shape) * noise
        return (M + M.T) * (1 - np.eye(M.shape[0]))

    def _pair(self, s):
        return self._sym(self.W0[s], 0.4), self._sym(self.W0[s], 0.4)

    def more(self, s, n):
        new = [self._pair(s) for _ in range(n)]
        self.pairs[s].extend(new)
        return new


@pytest.mark.parametrize("stratum", list(cst.AGE_STRATA))
def test_age_subsplit_kept_with_homoeologs_across(stratum):
    lab, tests, _ = halves.find_subgenomes(_AgeSources(np.random.default_rng(3), stratum), np.random.default_rng(3),
                                        age_ok=lambda g, z: True)
    acc = [t for t in tests if t["how"] == "age" and t["accepted"]]
    assert int(lab.max()) + 1 == 3 and acc and acc[0]["chosen"]["src"] == stratum


@pytest.mark.parametrize("stratum", list(cst.AGE_STRATA))
def test_age_subsplit_rejected_without_homoeologs_across(stratum):
    lab, tests, _ = halves.find_subgenomes(_AgeSources(np.random.default_rng(3), stratum), np.random.default_rng(3),
                                        age_ok=lambda g, z: False)
    assert int(lab.max()) + 1 == 2
    assert any(t["how"] == f"{stratum} split without homoeologs across it" and not t["accepted"] for t in tests)


def test_no_check_given_keeps_age_splits():
    lab, _, _ = halves.find_subgenomes(_AgeSources(np.random.default_rng(3)), np.random.default_rng(3))
    assert int(lab.max()) + 1 == 3


def test_top_split_is_never_gated(monkeypatch):
    """An age-rescued top split (parent None) is not checked: it must already be the pooled copies' best split."""
    real = halves.accept_split

    def top_by_age(cands, parent, strata_avail=None):
        if parent is None:
            return cands["all"], "age"
        return real(cands, parent, strata_avail)
    monkeypatch.setattr(halves, "accept_split", top_by_age)
    seen = []

    def never(g, z):
        seen.append(len(g))
        return False
    lab, tests, _ = halves.find_subgenomes(_AgeSources(np.random.default_rng(3)), np.random.default_rng(3), age_ok=never)
    assert tests[0]["accepted"] and tests[0]["how"] == "age"
    assert 12 not in seen and int(lab.max()) + 1 == 2


@pytest.mark.parametrize("noise,drawn", [(2.5, False), (0.6, True)])
def test_split_candidate_draws_more_pairs_when_borderline(noise, drawn):
    """Replicability within 2 s.e. of MIN_REP: the source draws more half-library pairs (up to MAX_PAIRS) until it is
    decisive; a decisive start draws none. The candidate records the pairs used, the s.e. and the borderline flag."""
    rng = np.random.default_rng(0)
    C = 8
    W = rng.random((C, C)) + np.kron(np.eye(2), np.ones((4, 4))) * 0.3
    W = W + W.T
    np.fill_diagonal(W, 0)
    noisy = lambda: (lambda M: (M + M.T) * (1 - np.eye(C)))(W + rng.random((C, C)) * noise)
    calls = []
    src = SimpleNamespace(W={"all": W}, pairs={"all": [(noisy(), noisy()) for _ in range(10)]})
    src.more = lambda s, n: calls.append(n) or src.pairs[s].extend((noisy(), noisy()) for _ in range(n))
    rep0, se0, _ = partition.rep_stats(list(src.pairs["all"]), np.arange(C), rng)          # the first 10 pairs
    c = halves.split_candidate(src, "all", np.arange(C), rng)
    assert c["n_pairs"] == len(src.pairs["all"]) >= 10 and "se" in c and "borderline" in c
    assert bool(calls) == drawn == (abs(rep0 - cst.MIN_REP) < 2 * se0)
    assert not c["borderline"] or c["n_pairs"] >= cst.MAX_PAIRS
