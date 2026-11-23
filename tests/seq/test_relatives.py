"""Closest relatives (seq.relatives): per group from shared lineages (sparse self-join, two tie rules), mutual and
reciprocal best pairs, per focal copy from shared sketch hashes, and LTRmix's reciprocal best candidate pairs."""
from pathlib import Path

import numpy as np
import pytest
from analyses.subgenomes import sgtools as sg
from scipy import sparse

from synltr.base import constants as cst
from synltr.seq import kmers, relatives


def _incidence(clades, weights, n, dtype):
    """Clade x copy incidence (ones) and weighted incidence (the clade's weight on each member), CSR."""
    rows = np.repeat(np.arange(len(clades)), [len(c) for c in clades])
    cols = np.concatenate([np.array(c) for c in clades])
    A = sparse.csr_matrix((np.ones(cols.size, dtype), (rows, cols)), shape=(len(clades), n))
    D = sparse.csr_matrix((np.asarray(weights, dtype)[rows], (rows, cols)), shape=(len(clades), n))
    return A, D


def test_best_by_group_ties():
    """Copy 0 shares one lineage of weight 1 with each of copies 4, 2 and 3 (group 1): a three-way tie. tie="lower"
    ranks tied relatives by copy index; tie="product" keeps the order the sparse product lists them in (phasing's
    nearest-relative table). exclude_same drops relatives in the copy's own group, else only the copy itself; copies
    without a group (-1) are never relatives (a focal copy without one has relatives in every group); missing ranks
    are -1 with relatedness 0."""
    group = np.array([0, 0, 1, 1, 1, -1])
    clades = [[0, 4], [0, 2], [0, 3], [1, 3], [0, 1], [0, 5]]
    A, D = _incidence(clades, [1, 1, 1, 2, 3, 5], 6, np.float32)
    J, V = relatives.best_by_group(A, D, group, top=3, exclude_same=True, tie="lower", threads=1)
    assert J.shape == V.shape == (6, 2, 3) and J.dtype == np.int32 and V.dtype == np.float32
    assert J[0, 1].tolist() == [2, 3, 4] and V[0, 1].tolist() == [1, 1, 1] and J[0, 0].tolist() == [-1, -1, -1]
    assert J[1, 1].tolist() == [3, -1, -1] and V[1, 1].tolist() == [2, 0, 0]
    assert J[5].tolist() == [[0, -1, -1], [-1, -1, -1]]                          # no group: every group is other
    R = (A.T.tocsr() @ D).tocoo()                                                 # the product's own listing
    row0 = [(c, v) for r, c, v in zip(R.row, R.col, R.data) if r == 0 and group[c] == 1]
    want = [c for c, v in sorted(row0, key=lambda cv: -cv[1])]                   # stable: ties keep that order
    P, _ = relatives.best_by_group(A, D, group, top=3, exclude_same=True, tie="product", threads=1)
    assert P[0, 1].tolist() == want and sorted(want) == [2, 3, 4]
    S, W = relatives.best_by_group(A, D, group, top=2, exclude_same=False, tie="lower", threads=1)
    assert S[0, 0].tolist() == [1, -1] and W[0, 0].tolist() == [3, 0] and S[0, 1].tolist() == [2, 3]
    assert S[1, 0].tolist() == [0, -1] and S[3, 1].tolist() == [-1, -1] and S[3, 0].tolist() == [1, 0]
    A64, D64 = _incidence(clades, [1, 1, 1, 2, 3, 5], 6, np.float64)
    J64, V64 = relatives.best_by_group(A64, D64, group, top=1, exclude_same=False, tie="lower", n_groups=3,
                                       threads=1)
    assert J64.shape == (6, 3, 1) and V64.dtype == np.float64 and J64[0, :, 0].tolist() == [1, 2, -1]
    with pytest.raises(ValueError, match="tie"):
        relatives.best_by_group(A, D, group, top=1, exclude_same=True, tie="first", threads=1)


@pytest.mark.parametrize("tie", ["lower", "product"])
def test_best_by_group_blocks_and_threads(tie):
    """Blocks of the self-join and pool workers change nothing; every relative is in another group, ranked by
    relatedness (brute force)."""
    rng = np.random.default_rng(3)
    n, G = 80, 4
    group = rng.integers(-1, G, n)
    clades = [sorted(rng.choice(n, int(rng.integers(2, 7)), replace=False).tolist()) for _ in range(150)]
    A, D = _incidence(clades, rng.integers(1, 5, len(clades)), n, np.float32)
    one = relatives.best_by_group(A, D, group, top=3, exclude_same=True, tie=tie, n_groups=G, threads=1)
    many = relatives.best_by_group(A, D, group, top=3, exclude_same=True, tie=tie, n_groups=G, block_nnz=40,
                                   threads=3)
    assert np.array_equal(one[0], many[0]) and np.array_equal(one[1], many[1])
    R = (A.T @ D).toarray()
    J, V = one
    for e in range(n):
        for g in range(G):
            cand = sorted((-R[e, f], f) for f in range(n) if group[f] == g and g != group[e] and R[e, f] > 0)
            vals = [-r for r, _ in cand[:3]]
            assert V[e, g, :len(vals)].tolist() == pytest.approx(vals) and (J[e, g, len(vals):] == -1).all()
            if tie == "lower":
                assert J[e, g, :len(vals)].tolist() == [f for _, f in cand[:3]]


def test_mutual_best_and_rbh_anchors():
    """mutual_best: copy pairs a < b each the other's closest relative in the other's group (any two groups, a
    copy's own included). rbh_anchors: pairs (a, b) across groups g < h, b a's closest relative in h and a b's in g."""
    best = np.array([[1, 2, 4], [0, 3, -1], [0, 3, 4], [0, 2, -1], [1, 2, -1]])
    group = np.array([0, 0, 1, 1, 2])
    assert relatives.mutual_best(best, group).tolist() == [[0, 1], [0, 2], [2, 3], [2, 4]]
    assert relatives.rbh_anchors(best, group).tolist() == [[0, 2], [2, 4]]
    none = relatives.rbh_anchors(np.full((3, 2), -1), np.array([0, 1, -1]))
    assert none.shape == (0, 2)


def test_sketch_relatives_find_family_at_old_divergence():
    """Families of 2-kb copies, 12% divergent between two sequences (the depth of an old progenitor split): each
    copy's top sketch relatives on the other sequence belong to its own family; a copy is never its own relative
    and relatives never sit on the focal copy's sequence."""
    rng = np.random.default_rng(4)
    seqs, seqof, fam = [], [], []
    for f in range(6):
        base = rng.choice(list("ACGT"), 2000)
        for side in (0, 1):
            anc = base.copy()
            hit = rng.random(anc.size) < 0.06                     # each side 6% from the ancestor: 12% apart
            anc[hit] = rng.choice(list("ACGT"), hit.sum())
            for _ in range(5):
                s = anc.copy()
                h2 = rng.random(s.size) < 0.01
                s[h2] = rng.choice(list("ACGT"), h2.sum())
                seqs.append("".join(s))
                seqof.append(side)
                fam.append(f)
    lib = sg.lib_from_seqs(seqs, seqof)
    S = kmers.sketch(lib, k=cst.SKETCH_K, s=cst.SKETCH_S)
    focal = np.arange(lib.n, dtype=np.int64)
    J, _ = relatives.top_shared(S, lib.seq.astype(np.int64), focal, m=3)
    fam = np.array(fam)
    assert (J >= 0).all()
    assert (fam[J] == fam[:, None]).all() and (lib.seq[J] != lib.seq[:, None]).all()


def test_top_shared_matches_brute_force():
    """Per focal copy: copies in another part (part >= 0) sharing >= 2 sketch hashes, hashes carried by more than cap
    copies ignored, most shared first, ties to the lower copy index; -1 = none."""
    rng = np.random.default_rng(11)
    n, s, pad = 60, 12, np.iinfo(np.uint64).max
    S = np.full((n, s), pad, np.uint64)
    for i in range(n):
        v = np.unique(rng.integers(1, 40, int(rng.integers(0, s + 1))).astype(np.uint64))
        S[i, :v.size] = v
    part = rng.integers(-1, 4, n)
    focal = np.flatnonzero(part >= 0)[::2]
    for cap in (2000, 6):
        J, V = relatives.top_shared(S, part, focal, m=4, cap=cap)
        carriers = {}
        for i in np.flatnonzero(part >= 0):
            for h in S[i][S[i] != pad].tolist():
                carriers[h] = carriers.get(h, 0) + 1
        for r, i in enumerate(focal):
            own = {h for h in S[i][S[i] != pad].tolist() if carriers.get(h, 0) <= cap}
            cnt = {c: len(own & set(S[c].tolist())) for c in range(n) if part[c] >= 0 and part[c] != part[i]}
            want = sorted(((-v, c) for c, v in cnt.items() if v >= 2))[:4]
            assert J[r, :len(want)].tolist() == [c for _, c in want] and (J[r, len(want):] == -1).all(), (cap, i)
            assert V[r, :len(want)].tolist() == [-v for v, _ in want]


def test_reciprocal_best_matches_both_ltrmix_copies():
    """Candidate pairs (a, b) with weights: the pairs that are each other's heaviest partner, equal weights to the
    lower index, one per a in a order. LTRmix's two copies exactly on 30 frozen candidate sets (sorted by (a, b), many
    equal weights): over all candidates (anchors), and over a subset's indices (pairs inside syntenic windows)."""
    with np.load(Path(__file__).resolve().parent.parent / "golden" / "ltrmix_route.npz") as z:
        ends = np.cumsum(z["rb_sizes"])[:-1]
        a, b, w, ok = (np.split(z[f"rb_{c}"], ends) for c in ("a", "b", "w", "ok"))
        sizes, na, nb = z["rb_sizes"], z["rb_na"].tolist(), z["rb_nb"].tolist()
        recip, mutual = z["rb_recip"], z["rb_mutual"]
    for c in range(sizes.size):
        ac, bc = a[c].astype(np.int64), b[c].astype(np.int64)
        got = relatives.reciprocal_best(ac, bc, w[c], na[c], nb[c])
        assert got.dtype == np.int64 and got.tolist() == recip[recip[:, 0] == c, 1].tolist()
        ci = np.flatnonzero(ok[c])
        sub = ci[relatives.reciprocal_best(ac[ci], bc[ci], w[c][ci], na[c], nb[c])]
        assert sub.tolist() == mutual[mutual[:, 0] == c, 1].tolist()


def test_reciprocal_best_ties_and_rule():
    """Equal weights go to the lower partner index on both sides, whatever the input order; tie names the rule."""
    a, b, w = np.array([1, 0, 0, 1]), np.array([0, 1, 0, 1]), np.array([2.0, 2.0, 2.0, 1.0])
    assert relatives.reciprocal_best(a, b, w, 2, 2).tolist() == [2]
    assert relatives.reciprocal_best(np.array([], np.int64), np.array([], np.int64), np.array([]), 1, 1).size == 0
    with pytest.raises(ValueError, match="tie"):
        relatives.reciprocal_best(a, b, w, 2, 2, tie="product")
