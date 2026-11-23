"""Painting: the ancestry HMM over copies, lineages as evidence, resolution, unphased sequences, copy ancestry."""
import itertools
import math
from types import SimpleNamespace

import numpy as np
import pytest
from analyses.subgenomes import sgtools as sg

from synltr import __version__
from synltr.analyses.subgenomes import inputs, lineages, painting, pipeline, run
from synltr.base import constants as cst
from synltr.seq import kmers


def test_hmm_matches_brute_force():
    """Posteriors and the Viterbi path equal brute-force enumeration over all state paths (3 states, 6 copies)."""
    rng = np.random.default_rng(41)
    n, N, sw = 6, 3, 0.15
    E = rng.normal(0, 1.5, (n, N))
    l0 = np.log(np.array([0.9, 0.05, 0.05]))
    lt = lambda i, j: math.log(1 - sw) if i == j else math.log(sw / (N - 1))
    paths = list(itertools.product(range(N), repeat=n))
    lp = np.array([l0[p[0]] + E[0, p[0]] + sum(lt(p[t - 1], p[t]) + E[t, p[t]] for t in range(1, n)) for p in paths])
    w = np.exp(lp - lp.max())
    w /= w.sum()
    marg = np.zeros((n, N))
    for p, x in zip(paths, w):
        marg[np.arange(n), list(p)] += x
    post, ll, _ = painting.hmm_posterior(E, l0, sw)
    assert np.allclose(post, marg, atol=1e-10) and abs(ll - (lp.max() + math.log(np.exp(lp - lp.max()).sum()))) < 1e-9
    assert tuple(painting.hmm_viterbi(E, l0, sw)) == paths[int(np.argmax(lp))]


def test_switch_rate_recovers_planted():
    """Baum-Welch recovers a planted switch rate (0.002 per copy) from noisy emissions (90% informative copies)."""
    rng = np.random.default_rng(42)
    seqs = []
    for _ in range(10):
        n = 3000
        st = np.cumsum(rng.random(n) < 0.002) % 2
        good = rng.random(n) < 0.9
        obs = np.where(good, st, 1 - st)
        E = np.where(np.arange(2)[None] == obs[:, None], math.log(0.9), math.log(0.1))
        seqs.append((np.ascontiguousarray(E), np.log(np.array([0.5, 0.5]))))
    sw = painting.switch_rate(seqs, 2)
    assert 0.0014 < sw < 0.0028, sw


def test_lineages_are_components_of_small_clades():
    """Copies sharing a small clade (<= NMAX_KNN carriers), directly or through others, are one lineage; a large clade
    (family-wide k-mers) does not join lineages."""
    big = list(range(0, 60))
    members = [[0, 1], [1, 2], [5, 6], big]
    cl = kmers.Clades(np.concatenate([np.array(m) for m in members]), np.array([len(m) for m in members]),
                   np.ones(4), np.full(4, 2))
    lab = painting.Lineages(cl, 100).labels([0, 1, 2, 5, 6, 7])
    assert lab[0] == lab[1] == lab[2] and lab[3] == lab[4] and len({lab[0], lab[3], lab[5]}) == 3


def test_resolution_track_follows_copy_density():
    """Copies every 10 kb on the first half of a sequence and every 40 kb on the second: the smallest detectable
    exchange (span of 5 consecutive copies) is 40 kb, then 160 kb."""
    mid = np.concatenate([np.arange(5_000, 2_000_000, 10_000), np.arange(2_020_000, 4_000_000, 40_000)]).astype(float)
    st, en, bp = painting.resolution_track(mid, 5, 4_000_000)
    assert bp[st < 1_800_000].max() == 40_000 and bp[(st > 2_300_000) & (st < 3_700_000)].min() == 160_000


def test_resolution_track_needs_enough_copies():
    """A sequence with fewer informative copies than needed has no resolution (no exchange detectable)."""
    st, en, bp = painting.resolution_track(np.array([10.0, 20.0, 30.0]), 5, 1_000_000)
    assert st.size == en.size == bp.size == 0


def test_resolution_track_sorts_nested_copies():
    """Copies given in start order with a nested one (midpoints not sorted): every value is a positive span."""
    mid = np.array([50_000.0, 20_000, 150_000, 250_000, 350_000])
    st, en, bp = painting.resolution_track(mid, 2, 400_000)
    assert bp.size and (bp > 0).all()


def test_resolution_track_spans_gaps():
    """Copies every 10 kb on 0-1 Mb and 4-5 Mb, none between: around 2.5 Mb the smallest detectable exchange spans
    the whole gap (>= 3 Mb), not the density of the copies before it."""
    mid = np.concatenate([np.arange(5_000, 1_000_000, 10_000), np.arange(4_005_000, 5_000_000, 10_000)]).astype(float)
    st, en, bp = painting.resolution_track(mid, 3, 5_000_000)
    w = (st <= 2_500_000) & (en > 2_500_000)
    assert w.sum() == 1 and bp[w][0] >= 3_000_000


def test_resolution_track_stops_at_the_outermost_copies():
    """Copies only in the first 100 kb of a 1-Mb sequence: no resolution rows beyond them (no promise without copies)."""
    st, en, bp = painting.resolution_track(np.arange(5_000, 100_000, 10_000).astype(float), 3, 1_000_000)
    assert st.size and (st + cst.PAINT_RES_WIN / 2 <= 95_000).all()


def test_posterior_tiles_sorted_and_unknown_in_gaps():
    """Ancestry tiles: in midpoint order, non-overlapping (a nested copy given out of order); a stretch without copies
    longer than the resolution R stays unknown (no tile), and tiles stop R/2 from the outermost copies."""
    mid = np.array([3_000.0, 1_000, 5_000, 5_000_000, 5_002_000])
    post = np.array([[0.9, 0.1], [0.8, 0.2], [0.7, 0.3], [0.1, 0.9], [0.2, 0.8]])
    st, en, val = painting.posterior_tiles(mid, post, 10_000_000, 100_000)
    assert (np.diff(st) > 0).all() and (st[1:] >= en[:-1]).all()
    assert not ((st < 4_000_000) & (en > 1_000_000)).any()
    assert st[0] == 0 and en[-1] == 5_002_000 + 50_000 and val.shape == (5, 2) and val[0, 0] == 0.8


def test_age_weights_discount_unconfined_ages():
    """Copies younger than 0.01 carry other-subgenome ancestry as often as uninformative copies (post-merger
    insertions): weight ~0; copies 0.03-0.05 are confined: weight 1; copies without K2P: weight 1."""
    rng = np.random.default_rng(43)
    n = 4000
    lib = inputs.Library()
    lib.n = n
    lib.k2p = rng.uniform(0, 0.05, n)
    lib.k2p[:10] = np.nan
    own = np.arange(n) % 2
    p_own = np.where(lib.k2p < 0.01, 0.5, 0.95)
    p_own[:10] = 0.95
    lP = np.log(np.where(np.arange(2)[None] == own[:, None], p_own[:, None], 1 - p_own[:, None]))
    w = painting.age_weights(lib, own, lP, np.ones(n), 2)
    assert np.median(w[lib.k2p < 0.008]) < 0.1 and np.median(w[lib.k2p > 0.03]) > 0.95 and (w[:10] == 1).all()


def test_age_weights_with_identical_ages():
    """Every copy with the same LTR K2P (e.g. a very young library): weights 1, no crash."""
    lib = inputs.Library()
    lib.n = 300
    lib.k2p = np.zeros(300)
    own = np.arange(300) % 2
    lP = np.log(np.where(np.arange(2)[None] == own[:, None], 0.9, 0.1))
    assert (painting.age_weights(lib, own, lP, np.ones(300), 2) == 1).all()


def _planted_painting(burst=30):
    """Sequence s1 (subgenome 0) and s2 (subgenome 1), 400 copies each, every copy 4-fold for its home; on s1, copies
    200-239 favour subgenome 1 and belong to 5 lineages; copies 300-(300 + burst) favour subgenome 1 but are one
    lineage (one clade holds them all)."""
    n = 800
    lib = inputs.Library()
    lib.n = n
    lib.seq = np.repeat([0, 1], 400)
    lib.start = np.tile(np.arange(400) * 10_000 + 1, 2)
    lib.end = lib.start + 5_000
    lib.seqnames, lib.seq_len = ["s1", "s2"], np.array([4_010_000, 4_010_000])
    lib.k2p = np.full(n, np.nan)
    home = lib.seq.copy()
    fav = home.copy()
    fav[200:240] = 1
    fav[300:300 + burst] = 1
    lP = np.where(np.arange(2)[None] == fav[:, None], math.log(0.8), math.log(0.2))
    groups = [list(range(200 + 8 * k, 208 + 8 * k)) for k in range(5)] + [list(range(300, 300 + burst))]
    groups += [[i, 400 + i] for i in range(800) if i < 400 and not (200 <= i < 240 or 300 <= i < 300 + burst)]
    cl = kmers.Clades(np.concatenate([np.array(g) for g in groups]), np.array([len(g) for g in groups]),
                   np.ones(len(groups)), np.full(len(groups), 2))
    ph = dict(N=2, names=["s1", "s2"], z=np.array([0, 1]), ek=lib.seq.copy(), clades=cl)
    return SimpleNamespace(seed=1), lib, ph, dict(ancestry=(lP, np.ones(n)))


def test_paint_counts_lineages():
    """The 5-lineage block is called an exchange from subgenome 1 at its copies' span; the one-lineage block is a
    segment of the best path but not called (one lineage is one piece of evidence)."""
    a, lib, ph, dt = _planted_painting()
    pa = painting.run_painting(a, lib, ph, dt)
    segs = [x for x in pa["segs"] if x["seq"] == "s1"]
    called = [x for x in segs if x["called"]]
    assert len(called) == 1 and called[0]["origin"] == 1 and called[0]["n_lineages"] == 5
    assert called[0]["start"] == lib.start[200] and called[0]["end"] == lib.end[239]
    one = [x for x in segs if x["start"] == lib.start[300]]
    assert one and not one[0]["called"] and one[0]["n_lineages"] == 1


def test_painting_none_for_one_subgenome():
    """One subgenome: nothing to paint."""
    assert painting.run_painting(SimpleNamespace(seed=1), inputs.Library(), dict(N=1), None) is None


def test_paint_counts_a_large_local_burst_once():
    """An 80-copy burst of one lineage (one clade larger than NMAX_KNN, most of its carriers inside the segment) is
    one piece of evidence: not called."""
    a, lib, ph, dt = _planted_painting(burst=80)
    pa = painting.run_painting(a, lib, ph, dt)
    one = [x for x in pa["segs"] if x["seq"] == "s1" and x["start"] == lib.start[300]]
    assert one and one[0]["n_lineages"] == 1 and not one[0]["called"]


def test_paint_segments_carry_breakpoint_flanks():
    """A called exchange's breakpoints lie between its outermost copies and the neighbouring informative copies: the
    flanks are the previous copy's end and the next copy's start; summary.json names the needed copies as copies."""
    a, lib, ph, dt = _planted_painting()
    pa = painting.run_painting(a, lib, ph, dt)
    c = [x for x in pa["segs"] if x["called"]][0]
    assert c["left_flank"] == lib.end[199] and c["right_flank"] == lib.start[240]
    s = painting.painting_summary(pa)
    assert s["copies_needed"] == pa["n_needed"] and "lineages_needed" not in s


def test_unphased_rows_without_evidence_are_na(tmp_path):
    """unphased_sequences.tsv: a sequence without informative copies has no odds (NA, not 0.00); an uncalled
    sequence with evidence names the subgenome its odds favor."""
    rows = [dict(seq="u1", n_ltr=3, n_inf=0, n_lin=0, g=-1, fav=-1, lo=0.0, conf="none"),
            dict(seq="u2", n_ltr=9, n_inf=4, n_lin=1, g=-1, fav=1, lo=0.4, conf="none"),
            dict(seq="u3", n_ltr=30, n_inf=20, n_lin=3, g=0, fav=0, lo=7.0, conf="decisive")]
    painting.write_unphased(tmp_path, rows)
    R = [l.split("\t") for l in (tmp_path / "unphased_sequences.tsv").read_text().splitlines()]
    d = {r[0]: dict(zip(R[0], r)) for r in R[1:]}
    assert d["u1"]["log10_odds_vs_next"] == "NA" and d["u1"]["favored"] == "NA"
    assert d["u2"]["ancestry"] == "NA" and d["u2"]["favored"] == "SG2" and d["u2"]["log10_odds_vs_next"] == "0.17"
    assert d["u3"]["ancestry"] == "SG1" and d["u3"]["confidence"] == "decisive"


@pytest.mark.slow
def test_copy_ancestry_matches_location(tmp_path):
    """Calibrated per-copy ancestry on a simulated allotetraploid: copies of lineages that arose after the split
    are assigned to their true subgenome of origin (>= 85%); copies older than the split are ambiguous."""
    import csv
    d, truth = sg.simulate(tmp_path, tree={"split": 2e6, "children": ["A", "B"]},
                         mergers=[{"time": 1e5, "lineages": ["A", "B"]}])
    a = sg.parse_args(["--ltr_fasta", str(d / "ltr.fa"), "--annotation", str(d / "annotation.gff3"), "--fai",
                       str(d / "genome.fai"), "--outdir", str(tmp_path / "o"), "-t", "4"])
    run.setup(a)
    lib = inputs.load_library(a)
    ph = pipeline.run_phasing(a, lib, lineages.load_or_build_clades(a, lib))
    lP, ninf = painting.copy_ancestry(ph["clades"], ph["ek"], ph["z"], ph["N"], np.random.default_rng(1))
    with open(tmp_path / "truth" / "copies.tsv") as fh:
        origin = {r["copy_id"]: r["origin"] for r in csv.DictReader(fh, delimiter="\t")}
    o = np.array([origin[i] for i in lib.ids])
    loc = np.array([truth[c.split(":")[0]] for c in lib.ids])
    own = ph["z"][np.maximum(ph["ek"], 0)]
    g2t = {g: max(set(loc[own == g]), key=list(loc[own == g]).count) for g in range(ph["N"])}
    call = np.array([g2t[g] for g in lP.argmax(1)])
    post = np.isin(o, ["A", "B"])
    assert post.mean() > 0.3 and (call[post] == o[post]).mean() >= 0.85


def test_exchanges_gff3_names_synltr(tmp_path):
    """exchanges.gff3 holds one feature per called exchange (1-based closed), its source column naming the package
    and its version; exchanges.bed the same exchange 0-based half-open."""
    seg = dict(seq="chr1", start=1001, end=5000, left_flank=0, right_flank=0, own=1, origin=0, n_ltr=4, n_support=3,
               n_against=0, n_lineages=2, llr=10.0, post=0.99, called=True)
    pa = dict(segs=[seg, dict(seg, start=8001, end=9000, called=False)], res=[], N=2, post=[], n_needed=5)
    painting.write_painting(SimpleNamespace(outdir=str(tmp_path)), pa)
    rows = [ln.split("\t") for ln in (tmp_path / "exchanges.gff3").read_text().splitlines() if ln[:1] != "#"]
    assert [r[1:5] for r in rows] == [[f"synltr-{__version__}", "region", "1001", "5000"]]
    assert (tmp_path / "exchanges.bed").read_text().split("\t")[:3] == ["chr1", "1000", "5000"]
