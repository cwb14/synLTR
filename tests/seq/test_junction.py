"""Exact junction refinement and the TSD checks (seq.junction)."""
import itertools
import math
import random

import builders
import pytest

from ltrtrace.base.constants import REFINE_PAD, SCORE
from ltrtrace.seq.junction import Gap, _junction, _slide, hamming, refine_gap, tsd_allowance

T = "CTCTC"  # target-site duplication
FLANK = 2000  # length of L and R
GRID = (-40, 0, 40)  # edge errors the local aligner can leave (tens of bp)


def _seq(n, rng, alphabet="ACGT"):
    return "".join(rng.choice(alphabet) for _ in range(n))


def _element(rng, avoid=""):
    """5-kb element: identical 500-bp LTRs (start TG, end CA) around a 4-kb INT.

    `avoid`: letters kept out of the 40 bp inside each terminal TG / CA, so a filler made of them
    cannot align there.
    """
    near = "".join(c for c in "ACGT" if c not in avoid)
    ltr = "TG" + _seq(40, rng, near) + _seq(416, rng) + _seq(40, rng, near) + "CA"
    return ltr + _seq(4000, rng) + ltr


def _site(seed, avoid=""):
    """L, e, R and A = L + T + e + T + R."""
    rng = random.Random(seed)
    L, R = _seq(FLANK, rng), _seq(FLANK, rng)
    e = _element(rng, avoid)
    return L, e, R, L + T + e + T + R


def _jitter(g, rng, d=40):
    """`g` with each edge moved by up to +-d bp (B edges kept in order), as the local aligner leaves it."""
    b0 = g.b0 + rng.randint(-d, d)
    return Gap(g.a0 + rng.randint(-d, d), g.a1 + rng.randint(-d, d), b0, max(b0, g.b1 + rng.randint(-d, d)))


def _grid(g):
    """`g` with every combination of GRID errors on its four edges (B edges kept in order)."""
    for da0, da1, db0, db1 in itertools.product(GRID, repeat=4):
        if g.b0 + db0 <= g.b1 + db1:
            yield Gap(g.a0 + da0, g.a1 + da1, g.b0 + db0, g.b1 + db1)


def _check_deletion(a, b, g):
    """`g` is a pure A deletion at its leftmost placement; every placement in its slide range turns A
    into B, and the range is maximal. (A and B share the left flank coordinates.)"""
    assert g.b0 == g.b1 == g.a0
    assert all(a[:g.a0 + t] + a[g.a1 + t:] == b for t in range(g.slide + 1))
    assert a[g.a0 - 1] != a[g.a1 - 1]  # no equal placement further left
    assert a[g.a0 + g.slide] != a[g.a1 + g.slide]  # none further right


# ---------------------------------------------------------------- empty sites, flank deletions and fillers


def test_empty_site():
    """A = L+T+e+T+R, B = L+T+R: the gap is e plus one TSD copy, B contiguous; the slide range holds
    both TSD-copy placements (and any chance microhomology beyond them)."""
    chance = 0
    for seed in range(20):
        L, e, R, a = _site(seed)
        b = L + T + R
        g = refine_gap(a, b, _jitter(Gap(len(L), len(L) + len(e) + 5, len(L), len(L)), random.Random(seed)))
        assert g.a1 - g.a0 == len(e) + 5 and g.b1 == g.b0
        assert g.a0 <= len(L) and g.a0 + g.slide >= len(L) + 5
        _check_deletion(a, b, g)
        chance += g.a0 < len(L) or g.slide > 5
    assert chance >= 5  # chance microhomology was exercised


def test_flank_deletion():
    """B = L[:-300] + R: the deletion also took 300 bp of left flank."""
    for seed in range(5):
        L, e, R, a = _site(seed)
        b = L[:-300] + R
        x = len(L) - 300
        g = refine_gap(a, b, _jitter(Gap(x, len(L) + len(e) + 10, x, x), random.Random(seed)))
        assert g.a1 - g.a0 == 300 + 5 + len(e) + 5 and g.b1 == g.b0
        assert g.a0 <= x <= g.a0 + g.slide
        _check_deletion(a, b, g)


def test_filler():
    """B = L + 7-bp filler + R: B-only sequence at the junction is reported exactly, slide 0.

    The filler is GGGGGGG beside element ends with no G within reach, so no filler base can align
    profitably and the insertion is the unique optimum (a GATTACA filler aligns better than it
    inserts: test_filler_that_aligns_better).
    """
    for seed in range(5):
        L, e, R, a = _site(seed, avoid="G")
        b = L + "G" * 7 + R
        true = Gap(len(L), len(L) + len(e) + 10, len(L), len(L) + 7)
        g = refine_gap(a, b, _jitter(true, random.Random(seed)))
        assert g.b1 - g.b0 == 7 and g.slide == 0
        assert g == true


def test_filler_that_aligns_better():
    """GATTACA facing the 7 A bases before R ('CA' + T = CACTCTC) has 2 matches, 5 mismatches:
    2 - 5*4 = -18, better than inserting it, -(6 + 2*7) = -20. The exact optimum is then a pure
    deletion 7 bp shorter than the generating event, with the filler aligned in the right flank."""
    for seed in range(5):
        L, e, R, a = _site(seed)
        b = L + "GATTACA" + R
        event = Gap(len(L), len(L) + len(e) + 10, len(L), len(L) + 7)  # the generating event
        g = refine_gap(a, b, _jitter(event, random.Random(seed)))
        assert g.b1 == g.b0 == g.a0 and g.a1 - g.a0 == len(e) + 3
        assert g.a0 <= len(L) <= g.a0 + g.slide


def test_large_b_insertion_unchanged():
    L, e, R, a = _site(3)
    b = L + T + builders.rand_seq(2000, random.Random(4)) + T + R
    gap = Gap(len(L) + 5, len(L) + 5 + len(e), len(L) + 5, len(L) + 2005)
    assert refine_gap(a, b, gap) is gap
    # the limit: more than 2 * pad bases of B-only sequence
    edge = Gap(gap.a0, gap.a1, gap.b0, gap.b0 + 2 * 50)
    assert refine_gap(a, b, edge, pad=50) is not edge
    over = Gap(gap.a0, gap.a1, gap.b0, gap.b0 + 2 * 50 + 1)
    assert refine_gap(a, b, over, pad=50) is over


# ---------------------------------------------------------------- short gaps, microhomology and perturbed inputs


def test_short_gap_midpoint_split():
    """120-bp deletion, so a1 - a0 < 2 * pad: the two A windows meet at the gap midpoint.

    With the edges pulled in to (x + 30, x + 90), the unsplit windows A [a0 - pad, a1) and
    [a0 + pad, a1 + pad) would skip A [x + 90, x + 180), the deletion's end included.
    """
    rng = random.Random(7)
    L, D, R = builders.rand_seq(FLANK, rng), builders.rand_seq(120, rng), builders.rand_seq(FLANK, rng)
    a, b, x = L + D + R, L + R, len(L)
    edges = ((x + 30, x + 90), (x + 40, x + 80), (x - 40, x + 80), (x + 40, x + 160), (x - 40, x + 160))
    for a0, a1 in edges:
        for b0, b1 in ((x, x), (x - 40, x), (x, x + 40), (x - 40, x + 40)):
            g = refine_gap(a, b, Gap(a0, a1, b0, b1))
            assert g.a1 - g.a0 == 120 and g.a0 <= x <= g.a0 + g.slide
            _check_deletion(a, b, g)


def test_slide_microhomology():
    """A = X + M + D + M + Y, B = X + M + Y, M a 4-bp direct repeat: the deletion (D and one M) has
    4 + chance further placements that give B; the gap is the leftmost, its slide reaches the rightmost."""
    M = "ACGT"
    for seed in range(10):
        rng = random.Random(seed)
        X, D, Y = builders.rand_seq(FLANK, rng), builders.rand_seq(1000, rng), builders.rand_seq(FLANK, rng)
        if seed == 0:  # no chance extension: X[-1] != D[-1] and D[0] != Y[0]
            X, D, Y = X[:-1] + "A", "C" + D[1:-1] + "G", "T" + Y[1:]
        a, b, n = X + M + D + M + Y, X + M + Y, len(D) + 4
        ok = [p for p in range(len(X) - 50, len(X) + 55) if a[:p] + a[p + n:] == b]  # every placement
        assert ok == list(range(ok[0], ok[-1] + 1)) and ok[0] <= len(X) and ok[-1] >= len(X) + 4
        g = refine_gap(a, b, _jitter(Gap(len(X), len(X) + n, len(X), len(X)), rng))
        assert (g.a0, g.a1, g.b0, g.b1, g.slide) == (ok[0], ok[0] + n, ok[0], ok[0], ok[-1] - ok[0])
        _check_deletion(a, b, g)
        if seed == 0:
            assert g.a0 == len(X) and g.slide == 4


def test_perturbed_inputs():
    """Edges off by up to 40 bp either way on A and B (every corner of the grid) give one exact gap."""
    L, e, R, a = _site(11, avoid="G")
    x, y = len(L), len(L) + len(e) + 10  # A: end of L, start of R
    filler = L + "G" * 7 + R
    for b, true in ((L + T + R, Gap(x, y - 5, x, x)),  # empty site
                    (L[:-300] + R, Gap(x - 300, y, x - 300, x - 300)),  # flank deletion
                    (filler, Gap(x, y, x, x + 7))):
        want = refine_gap(a, b, true)
        if b is filler:
            assert want == true
        else:
            assert want.a1 - want.a0 == true.a1 - true.a0
            _check_deletion(a, b, want)
        for g in _grid(true):
            assert refine_gap(a, b, g) == want, g


# ---------------------------------------------------------------- ties, N, edge cases


def test_b_only_sequence_beside_tsd():
    """B-only sequence next to B's TSD copy belongs to the junction ([b0, b1)), not to a flank indel
    beside A's other TSD copy, which scores the same. The fillers cannot align: AAAAAAA after the TSD
    (it faces 'TG' + no A, or CTCTC), GGGGGGG before it (it faces CTCTC + 'TG', or 'CA' + no G)."""
    for seed in range(3):
        L, e, R, a = _site(seed, avoid="A")
        x, y = len(L), len(L) + len(e) + 10
        after = Gap(x + 5, y, x + 5, x + 12)  # B's T faces A's first copy: A loses e + the second copy
        assert refine_gap(a, L + T + "A" * 7 + R, after) == after
        L, e, R, a = _site(seed, avoid="G")
        before = Gap(x, y - 5, x, x + 7)  # B's T faces A's second copy: A loses the first copy + e
        assert refine_gap(a, L + "G" * 7 + T + R, before) == before


def test_n_in_b_window_returns_input():
    """An N (any case) anywhere in the B window Bj = B[b0 - pad, b1 + pad) returns the input gap: the
    free A deletion would swallow an assembly gap near the junction as deleted flank."""
    L, e, R, a = _site(0)
    x, n = len(L), len(e) + 5
    gap = Gap(x, x + n, x, x)

    def one_n(s, i):
        return s[:i] + "N" + s[i + 1:]

    for b in (L[:-70] + "N" * 50 + L[-20:] + T + R,  # G_near-like: N run 20-70 bp before the TSD
              L[:-70] + "n" * 50 + L[-20:] + T + R,  # soft-masked
              L + T + "N" * 100 + R,  # N run at the junction
              one_n(L, x - 150) + T + R, L + T + one_n(R, 144)):  # one N at either window end
        assert refine_gap(a, b, gap) is gap
    for b in (one_n(L, x - 151) + T + R, L + T + one_n(R, 145),  # one N just outside the window
              L[:-550] + "N" * 50 + L[-500:] + T + R):  # N run 500 bp away
        g = refine_gap(a, b, gap)
        assert g is not gap and g.a1 - g.a0 == n and g.b0 == g.b1


def test_long_repeat_left_normalized():
    """Equal placements spanning more than the window: from the rightmost placement (where the local
    aligner can leave a deletion) or inside the range, the gap comes back at the leftmost placement
    with the full slide. A solo LTR (400-bp LTRs > pad) and a tandem array missing one 7-bp unit."""
    rng = random.Random(100)
    for case in range(4):
        if case < 3:  # solo LTR
            L, R = _seq(FLANK, rng), _seq(FLANK, rng)
            ltr = "TG" + _seq(396, rng) + "CA"
            a, b = L + T + ltr + _seq(4000, rng) + ltr + T + R, L + T + ltr + T + R
        else:  # tandem array
            X, Y = _seq(FLANK, rng), _seq(FLANK, rng)
            a, b = X + "ACGTTGC" * 100 + Y, X + "ACGTTGC" * 99 + Y
        n = len(a) - len(b)
        ok = [p for p in range(FLANK - 100, FLANK + 800) if a[:p] + a[p + n:] == b]  # every placement
        assert ok == list(range(ok[0], ok[-1] + 1)) and ok[-1] - ok[0] > 2 * REFINE_PAD
        want = Gap(ok[0], ok[0] + n, ok[0], ok[0], ok[-1] - ok[0])
        r, mid = ok[-1], (ok[0] + ok[-1]) // 2
        # off-diagonal edges (a0 and b0 moved apart) only outside the array: inside it the left window
        # start realigns by whole periods, which the free deletion then absorbs; the tandem input is
        # moved along the diagonals instead, as a local aligner's run ends are
        off, on = Gap(r - 40, r + n + 30, r - 20, r - 20), Gap(r - 40, r + n + 40, r - 40, r + 40)
        for g in (Gap(r, r + n, r, r), off if case < 3 else on, Gap(mid, mid + n, mid, mid)):
            assert refine_gap(a, b, g) == want, g
        _check_deletion(a, b, want)


def test_no_indel():
    """An empty gap inside collinear sequence stays empty, slide 0."""
    s = builders.rand_seq(3000, random.Random(9))
    assert refine_gap(s, s, Gap(1500, 1500, 1500, 1500)) == Gap(1500, 1500, 1500, 1500)


def test_soft_masked_input():
    """Lowercase (soft-masked) bases count as uppercase, in the DP and in the slide."""
    crossed = 0
    for seed in range(8):
        L, e, R, a = _site(seed)
        b = L + T + R
        g = Gap(len(L), len(L) + len(e) + 5, len(L), len(L))
        want = refine_gap(a, b, g)
        masked = L + T + e.lower() + T + R  # only the element masked
        assert refine_gap(masked, b, g) == refine_gap(a.lower(), b, g) == refine_gap(a, b.lower(), g) == want
        crossed += want.a0 < len(L) or want.slide > 5  # the slide range crosses the mask edge
    assert crossed


def test_rejects_bad_input():
    s = "ACGT" * 100
    for g in (Gap(10, 5, 0, 0), Gap(0, 10, 5, 4), Gap(-1, 10, 0, 0), Gap(0, 401, 0, 0), Gap(0, 10, 0, 401)):
        with pytest.raises(ValueError, match="gap"):
            refine_gap(s, s, g)
    with pytest.raises(ValueError, match="pad"):
        refine_gap(s, s, Gap(10, 20, 10, 10), pad=0)


# ---------------------------------------------------------------- the DP against brute force


def _ref_matrix(x, y):
    """H[i][j] = (score, -gaps opened) of the best global alignment of x[:i] and y[:j], from full Gotoh
    matrices with tuple scores compared lexicographically. A gap of n bases costs open + extend * n;
    a pair is a match only if identical and not N."""
    ma, mm, go, ge = SCORE
    neg = (float("-inf"), 0)
    n, m = len(x), len(y)
    M = [[neg] * (m + 1) for _ in range(n + 1)]
    D = [[neg] * (m + 1) for _ in range(n + 1)]  # x[i-1] against a gap
    In = [[neg] * (m + 1) for _ in range(n + 1)]  # y[j-1] against a gap
    H = [[neg] * (m + 1) for _ in range(n + 1)]
    H[0][0] = (0, 0)
    for i in range(n + 1):
        for j in range(m + 1):
            if i and j:
                s = ma if x[i - 1] == y[j - 1] != "N" else -mm
                M[i][j] = (H[i - 1][j - 1][0] + s, H[i - 1][j - 1][1])
            if i:  # open a gap after any state, or extend one
                h, g = H[i - 1][j], D[i - 1][j]
                D[i][j] = max((h[0] - go - ge, h[1] - 1), (g[0] - ge, g[1]))
            if j:
                h, g = H[i][j - 1], In[i][j - 1]
                In[i][j] = max((h[0] - go - ge, h[1] - 1), (g[0] - ge, g[1]))
            if i or j:
                H[i][j] = max(M[i][j], D[i][j], In[i][j])
    return H


def _ref_junction(al, ar, bj):
    """(i, j1, k, j2) by enumerating every split: al[:i] ~ bj[:j1] and ar[k:] ~ bj[j2:] (each a global
    alignment, the right one computed forward on the suffixes), bj[j1:j2] one B gap. Best score, then
    fewest gaps opened inside the parts, then smallest j2 - j1, i, k, j1."""
    go, ge = SCORE[2], SCORE[3]
    fwd = _ref_matrix(al, bj)
    best = None
    for k in range(len(ar) + 1):
        for j2 in range(len(bj) + 1):
            r = _ref_matrix(ar[k:], bj[j2:])[-1][-1]
            for i in range(len(al) + 1):
                for j1 in range(j2 + 1):
                    f = fwd[i][j1]
                    pen = 0 if j1 == j2 else go + ge * (j2 - j1)
                    key = (-(f[0] + r[0] - pen), -(f[1] + r[1]), j2 - j1, i, k, j1)
                    if best is None or key < best[0]:
                        best = (key, (i, j1, k, j2))
    return best[1]


def test_dp_matches_brute_force():
    """Small random windows over few letters (many ties) and N: the DP picks the brute-force split."""
    rng = random.Random(13)
    for _ in range(250):
        abc = rng.choice(("AC", "ACG", "ACN", "ACGT"))
        al, ar, bj = (_seq(rng.randint(0, hi), rng, abc) for hi in (6, 6, 8))
        assert _junction(al, ar, bj) == _ref_junction(al, ar, bj), (al, ar, bj)
    # structured: a deletion with a direct repeat, a B insertion,
    # a mismatch run; then cases where several i (or k) tie at one column and the smallest must win for the
    # split to be right
    for al, ar, bj in (("ACGTTG", "CATGCA", "ACGTGCA"), ("AACC", "GGTT", "AACCNNGGTT"),
                       ("ACGT", "TTGCA", "ACGAAAGCA"),
                       ("AACCAAA", "A", "CAAAA"), ("GACACGG", "G", "AGGGCCCGA"),  # tied i
                       ("", "NAAANAA", "ACANCACN"), ("C", "TGAAACC", "TCACATAT")):  # tied k
        assert _junction(al, ar, bj) == _ref_junction(al, ar, bj), (al, ar, bj)


@pytest.mark.parametrize("k", [0, 1, 63, 64, 65, 191, 192, 193, 1000])
def test_slide_doubling_steps(k):
    """_slide compares in doubling steps (64, then 128, ...); a direct repeat of k bases
    after the deleted segment is measured exactly across the step edges (64, 192), ignoring case, and up to
    the sequence end."""
    rng = random.Random(k)
    rep, filler, after = builders.rand_seq(k, rng), "C" + builders.rand_seq(299, rng), "G" + builders.rand_seq(49, rng)
    a0, a1 = 10, 10 + k + 300  # deletes rep + filler; the next bases are rep again, then "G" against "C"
    assert _slide("T" * 10 + rep + filler + rep.lower() + after, a0, a1) == k
    assert _slide("T" * 10 + rep + filler + rep, a0, a1) == k  # the repeat ends the sequence
    assert _slide("T" * 10 + rep + filler, a1, a1) == 0  # no deletion


# ---------------------------------------------------------------- TSD checks


def test_tsd_allowance():
    """max(1, floor(L * (K + d) + 0.5)), d = 1 - identity (0 when nan)."""
    assert tsd_allowance(5, 0.0, 1.0) == 1
    assert tsd_allowance(5, None, math.nan) == 1 and tsd_allowance(5, math.nan, 1.0) == 1  # K unknown
    assert tsd_allowance(5, 0.3, 1.0) == 2  # 1.5 + 0.5
    assert tsd_allowance(5, 0.2, 0.9) == 2  # d = 0.1: 5 * 0.3 + 0.5 = 2.0
    assert tsd_allowance(5, 0.2, 1.0) == 1  # 1.0 + 0.5 = 1.5
    assert tsd_allowance(6, 0.5, 1.0) == 3  # 3.0 + 0.5


def test_hamming_ignores_case_and_never_matches_n():
    assert hamming("ACGTA", "ACGTA") == 0 and hamming("acgta", "ACGTT") == 1
    assert hamming("ACNTA", "ACNTA") == 1 and hamming("NNNNN", "nnnnn") == 5
