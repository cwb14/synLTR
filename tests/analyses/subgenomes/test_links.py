"""Links between sequences: clades shared, the degree-corrected null, nearest-relative votes."""
import numpy as np
import pytest
from analyses.subgenomes import sgtools as sg

from ltrtrace.analyses.subgenomes import links
from ltrtrace.base import constants as cst
from ltrtrace.base import runtime
from ltrtrace.seq import kmers


def test_link_matrix_counts_each_clade_once_and_ignores_same_sequence():
    cl = kmers.Clades(np.array([0, 1, 2, 0, 3]), np.array([3, 2]), np.array([50, 1]), np.array([2, 2]))
    ek = np.array([0, 0, 1, 2])          # copies 0 and 1 share sequence 0
    W = links.link_matrix(cl, ek, 3, np.ones(len(ek), bool))       # every copy: the all-copies links
    # clade 0 spans sequences {0, 1}: 1/(2-1) = 1 for the pair; clade 1 spans {0, 2}: 1
    assert np.allclose(W, [[0, 1, 1], [1, 0, 0], [1, 0, 0]])


def test_expected_matches_row_sums():
    rng = np.random.default_rng(0)
    W = rng.random((6, 6))
    W = W + W.T
    np.fill_diagonal(W, 0)
    E = links.expected(W)
    assert np.allclose(E.sum(1), W.sum(1), rtol=1e-6) and np.allclose(np.diag(E), 0)


def test_knn_votes_symmetric():
    rel_j = np.array([[[-1], [1], [2]], [[0], [-1], [2]], [[0], [1], [-1]]], np.int32)
    rel_r = np.array([[[0], [3], [1]], [[3], [0], [1]], [[1], [1], [0]]], np.float32)
    W = links.knn_votes(rel_j, rel_r, np.array([0, 1, 2]), 3, np.ones(3, bool))
    assert np.allclose(W, W.T) and W[0, 1] > 0 and np.allclose(np.diag(W), 0)


def _link_matrix_ref(cl, ek, C, mask):
    """Reference link matrix from scipy sparse products: clades restricted to the members in mask, each clade adding
    1/(m-1) to every pair of the m sequences it spans."""
    from scipy import sparse
    sub = cl.restrict(mask)
    sid = sub.sid()
    e = ek[sub.members]
    ok = e >= 0
    X = sparse.csr_matrix((np.ones(int(ok.sum())), (sid[ok], e[ok])), shape=(sub.n, C))
    X.data = (X.data > 0).astype(float)
    X.eliminate_zeros()
    m = np.asarray(X.sum(1)).ravel()
    w = np.where(m >= 2, 1.0 / np.maximum(m - 1, 1), 0.0)
    W = (X.T @ sparse.diags(w) @ X).toarray()
    np.fill_diagonal(W, 0)
    return W


@pytest.mark.parametrize("seed", range(3))
def test_link_matrix_matches_reference(seed):
    rng = np.random.default_rng(seed)
    n_el, C = 400, 7
    sizes = rng.integers(2, 30, 300)
    members = np.concatenate([np.sort(rng.choice(n_el, s, replace=False)) for s in sizes])
    cl = kmers.Clades(members, sizes, np.ones(sizes.size, np.int64), np.full(sizes.size, 2))
    ek = rng.integers(-1, C, n_el)
    for mask in (np.ones(n_el, bool), rng.random(n_el) < 0.5):
        assert np.allclose(links.link_matrix(cl, ek, C, mask), _link_matrix_ref(cl, ek, C, mask))


def test_link_chunks_bounded_and_thread_independent():
    import numba
    assert runtime.n_chunks(3000, 10 ** 6) * 3000 ** 2 * 8 <= 1e9 and runtime.n_chunks(10, 10 ** 6) == 64
    rng = np.random.default_rng(7)
    sizes = rng.integers(2, 20, 2000)
    members = np.concatenate([np.sort(rng.choice(500, s, replace=False)) for s in sizes])
    cl = kmers.Clades(members, sizes, np.ones(sizes.size, np.int64), np.full(sizes.size, 2))
    ek = rng.integers(0, 9, 500)
    mask = rng.random(500) < 0.6
    out = []
    for t in (1, min(4, numba.config.NUMBA_NUM_THREADS)):
        numba.set_num_threads(t)
        out.append(links.link_matrix(cl, ek, 9, mask))
    assert np.array_equal(out[0], out[1])


def test_expected_warns_when_not_converged(caplog):
    W = np.ones((6, 6))
    W[0, 1:] = W[1:, 0] = 1e6                                     # one hub holds nearly all links
    np.fill_diagonal(W, 0)
    links._EXPECTED_WARNED = False
    with caplog.at_level("WARNING", logger="ltrtrace"):
        links.expected(W, iters=3)
    assert [r.levelname for r in caplog.records if "did not converge" in r.getMessage()] == ["WARNING"]


def _knn_ref(rel_j, rel_r, ek, C, mask):
    """A numpy version (reference)."""
    valid = (rel_j >= 0) & mask[np.maximum(rel_j, 0)]
    best = np.where(valid, rel_r, 0).max(2)
    best[~mask] = 0
    k = min(cst.K_NEIGH, C - 1)
    top = np.argpartition(-best, k - 1, axis=1)[:, :k]
    tv = np.take_along_axis(best, top, 1)
    ok = tv > 0
    rows = np.repeat(np.arange(best.shape[0])[:, None], k, 1)
    W = np.zeros((C, C))
    np.add.at(W, (ek[rows[ok]], top[ok]), 1.0)
    W = W + W.T
    np.fill_diagonal(W, 0)
    return W


@pytest.mark.parametrize("seed", range(3))
def test_knn_kernel_matches_reference(seed):
    rng = np.random.default_rng(seed)
    n, C, T = 300, 6, 3
    ek = rng.integers(-1, C, n)
    rel_j = rng.integers(-1, n, (n, C, T)).astype(np.int32)
    for e in range(n):
        if ek[e] >= 0:
            rel_j[e, ek[e]] = -1                                  # never a relative on its own sequence
    rel_r = rng.random((n, C, T)).astype(np.float32)              # continuous: no ties
    for mask in (ek >= 0, (ek >= 0) & (rng.random(n) < 0.5)):
        assert np.allclose(links.knn_votes(rel_j, rel_r, ek, C, mask), _knn_ref(rel_j, rel_r, ek, C, mask))


def test_top_relatives_parallel_equals_serial():
    rng = np.random.default_rng(4)
    seqs = sg.random_family(rng, 60)
    lib = sg.lib_from_seqs(seqs, (list(range(4)) * 15)[:60])
    cl = kmers.clades(lib, lib.seq, k=15, cap=40, threads=1)
    ek = lib.seq.astype(np.int64)
    a = links.top_relatives(cl, ek, 4, lib.n, 3, block_nnz=50, threads=1)
    b = links.top_relatives(cl, ek, 4, lib.n, 3, block_nnz=50, threads=4)
    assert np.array_equal(a[0], b[0]) and np.array_equal(a[1], b[1]) and (a[0] >= 0).sum() > 100
