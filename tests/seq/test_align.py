"""Phasing's global aligner (seq.align): an edlib path fixes a band, an affine-gap realignment inside it counts aligned
sites, transitions, transversions, gap events and gap bases; the orientation rule and the pooled copy-pair counts;
LTRmix's orientation rule, unit-cost path and position map."""
import os
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from analyses.subgenomes import sgtools as sg

from ltrtrace.base import constants as cst
from ltrtrace.io.seqio import encode
from ltrtrace.seq import align

GOLDEN = Path(__file__).resolve().parent.parent / "golden"

def codes(x):
    return encode(x.encode())


def test_global_counts_known_alignment():
    """A hand-checked pair: q (40 bp, an N at 20) and t = q with G>A at 4 and T>C at 10 (transitions), G>C at 31 (a
    transversion) and AGT at 25-27 deleted. Aligned columns with both bases ACGT: 40 - 3 deleted - 1 N = 36; one gap
    of three bases. Any band and the default scores give it; a trailing or a whole-sequence gap is one gap event."""
    q = "ATGCGTACCTTAGGCATCGANTGCCAGTTACGGATCCTAG"
    t = "ATGCATACCTCAGGCATCGANTGCCTACCGATCCTAG"
    c = align.global_counts(codes(q), codes(t))
    assert c == (36, 2, 1, 1, 3) and (c.sites, c.ts, c.tv, c.gap_events, c.gap_bases) == (36, 2, 1, 1, 3)
    assert align.global_counts(codes(q), codes(t), band=1, scores=cst.SCORES) == c
    assert align.global_counts(codes(t), codes(q)) == (36, 2, 1, 1, 3)
    assert align.global_counts(codes(q), codes(q[:-3])) == (36, 0, 0, 1, 3)        # trailing gap
    assert align.global_counts(codes(""), codes(t)) == (0, 0, 0, 1, 37)
    assert align.global_counts(codes(""), codes("")) == (0, 0, 0, 0, 0)


def test_align_counts_planted():
    rng = np.random.default_rng(1)
    q = "".join(rng.choice(list("ACGT"), 300))
    t = sg.mutate(rng, q, 5, 3)
    t = t[:150] + "GG" + t[150:]                                    # a 2-base insertion
    assert align.global_counts(codes(q), codes(t))[:3] == (300, 5, 3)


def test_align_counts_unbiased_at_split_divergence():
    """13.5% substitutions, no indels: unit-cost alignment swaps mismatch pairs for indel pairs and undercounts;
    the affine realignment recovers every planted substitution."""
    rng = np.random.default_rng(3)
    q = "".join(rng.choice(list("ACGT"), 3000))
    t = sg.mutate(rng, q, 270, 135)
    assert align.global_counts(codes(q), codes(t))[:3] == (3000, 270, 135)
    t2 = t[:1500] + "ACGTACGTAC" + t[1500:]                            # plus a 10-base insertion
    assert align.global_counts(codes(q), codes(t2))[:3] == (3000, 270, 135)


def test_align_counts_long_indel():
    """A 500-bp insertion (a nested element or internal deletion) in a 5-kb copy with 3% substitutions: the counts
    equal the planted ones (the band must hold the single-gap path, not only edlib's drifting path)."""
    rng = np.random.default_rng(12)
    q = "".join(rng.choice(list("ACGT"), 5000))
    t = sg.mutate(rng, q, 100, 50)
    t = t[:2500] + "".join(rng.choice(list("ACGT"), 500)) + t[2500:]
    assert align.global_counts(codes(q), codes(t)) == (5000, 100, 50, 1, 500)


def _unbanded(a, b):
    return align._gotoh_band(a, b, np.zeros(a.size + 1, np.int64), np.full(a.size + 1, b.size, np.int64), *cst.SCORES)


def test_align_counts_equals_unbanded_with_several_indels():
    """Random pairs with 1-4 indels of 1-600 bp (insertions and deletions, possibly cancelling) and 2-12%
    substitutions: the banded counts equal the unbanded affine alignment's."""
    rng = np.random.default_rng(13)
    for _ in range(25):
        q = "".join(rng.choice(list("ACGT"), int(rng.integers(1500, 4000))))
        t = list(q)
        hit = rng.random(len(t)) < rng.uniform(0.02, 0.12)
        t = np.array(t)
        t[hit] = rng.choice(list("ACGT"), hit.sum())
        t = "".join(t)
        for _ in range(int(rng.integers(1, 5))):
            p, L = int(rng.integers(0, len(t))), int(rng.integers(1, 600))
            t = t[:p] + "".join(rng.choice(list("ACGT"), L)) + t[p:] if rng.random() < 0.5 else t[:p] + t[p + L:]
        a, b = codes(q), codes(t)
        full = _unbanded(a, b)
        assert align.global_counts(a, b) == tuple(full[:3]) + tuple(full[4:])


def test_align_counts_equals_unbanded_for_fragments():
    """A 1.5-3-kb fragment against a 6-8-kb element (truncated copies are common in real libraries), 8-20%
    substitutions and a few indels: banded counts equal the unbanded affine alignment's (the band widens whenever
    the best path touches its edge)."""
    rng = np.random.default_rng(14)
    for _ in range(30):
        long_ = "".join(rng.choice(list("ACGT"), int(rng.integers(6000, 8000))))
        s0 = int(rng.integers(0, len(long_) - 3000))
        frag = np.array(list(long_[s0:s0 + int(rng.integers(1500, 3000))]))
        hit = rng.random(frag.size) < rng.uniform(0.08, 0.2)
        frag[hit] = rng.choice(list("ACGT"), hit.sum())
        frag = "".join(frag)
        for _ in range(int(rng.integers(0, 3))):
            p, L = int(rng.integers(0, len(frag))), int(rng.integers(1, 200))
            frag = frag[:p] + frag[p + L:]
        a, b = codes(frag), codes(long_)
        full = _unbanded(a, b)
        assert align.global_counts(a, b) == tuple(full[:3]) + tuple(full[4:])


def _gotoh_reference(q, t, lo, hi, sm, sx, go, ge):
    """Plain full-matrix affine alignment inside a band (cells outside [lo[i], hi[i]] are unreachable), ties: diagonal,
    then gap in the query (E), then gap in the target (F); a gap extends rather than reopens on ties. Returns sites,
    transitions, transversions, gap events (runs of gap columns) and gap bases of its traceback."""
    n, m, NEG = len(q), len(t), -10 ** 12
    H = [[NEG] * (m + 1) for _ in range(n + 1)]
    E = [[NEG] * (m + 1) for _ in range(n + 1)]
    F = [[NEG] * (m + 1) for _ in range(n + 1)]
    src, eext, fext = {}, {}, {}
    H[0][0] = 0
    for i in range(n + 1):
        for j in range(max(lo[i], 0), hi[i] + 1):
            if i == 0 and j == 0:
                continue
            if j > lo[i]:
                eo, ee = H[i][j - 1] - go - ge, E[i][j - 1] - ge
                E[i][j], eext[i, j] = (ee, True) if ee >= eo else (eo, False)
            if i > 0 and lo[i - 1] <= j <= hi[i - 1]:
                fo, fe = H[i - 1][j] - go - ge, F[i - 1][j] - ge
                F[i][j], fext[i, j] = (fe, True) if fe >= fo else (fo, False)
            d = NEG
            if i > 0 and j > 0 and lo[i - 1] <= j - 1 <= hi[i - 1] and H[i - 1][j - 1] > NEG // 2:
                d = H[i - 1][j - 1] + (0 if q[i - 1] > 3 or t[j - 1] > 3 else (sm if q[i - 1] == t[j - 1] else sx))
            H[i][j], src[i, j] = max((d, 0), (E[i][j], 1), (F[i][j], 2), key=lambda v: (v[0], -v[1]))
    i, j, state, sites, ts, tv, events, bases = n, m, 0, 0, 0, 0, 0, 0
    while i > 0 or j > 0:
        if state == 0:
            state = src[i, j]
            if state == 0:
                if q[i - 1] < 4 and t[j - 1] < 4:
                    sites += 1
                    ts += q[i - 1] != t[j - 1] and (q[i - 1] ^ t[j - 1]) == 2
                    tv += q[i - 1] != t[j - 1] and (q[i - 1] ^ t[j - 1]) != 2
                i, j = i - 1, j - 1
            else:
                events += 1
        elif state == 1:
            state = 1 if eext[i, j] else 0
            j -= 1
            bases += 1
        else:
            state = 2 if fext[i, j] else 0
            i -= 1
            bases += 1
    return sites, ts, tv, events, bases


def _random_pairs(rng, count):
    for it in range(count):
        q = rng.integers(0, 5 if it % 3 == 0 else 4, int(rng.integers(1, 40))).astype(np.uint8)
        t = q.copy()
        hit = rng.random(t.size) < rng.random() * 0.3
        t[hit] = rng.integers(0, 4, hit.sum())
        p, L = int(rng.integers(0, t.size + 1)), int(rng.integers(1, 15))
        t = np.concatenate([t[:p], rng.integers(0, 4, L), t[p:]]) if it % 2 else np.concatenate([t[:p], t[p + L:]])
        t = t.astype(np.uint8)
        if t.size == 0:
            continue
        lo, hi = np.zeros(q.size + 1, np.int64), np.full(q.size + 1, t.size, np.int64)
        if it % 4 == 1:
            r = align.edlib.align(align.ALPH[q].tobytes(), align.ALPH[t].tobytes(), mode="NW", task="path")
            m = align._CIGAR.findall(r["cigar"].encode())
            lo, hi = align._band(q.size, t.size, np.array([int(x) for x, _ in m], np.int64),
                                 np.array([align._OPCODE[o[0]] for _, o in m], np.int64), int(rng.integers(1, 6)))
        yield it, q, t, lo, hi


def test_gotoh_band_matches_a_plain_reference():
    """The banded affine kernel (sites, transitions, transversions of its traceback) equals a plain full-matrix
    implementation with the same tie rules, on short random pairs with N bases, substitutions, indels and bands."""
    for it, q, t, lo, hi in _random_pairs(np.random.default_rng(15), 150):
        assert tuple(align._gotoh_band(q, t, lo, hi, *cst.SCORES)[:3]) == _gotoh_reference(q, t, lo, hi, *cst.SCORES)[:3], it


def test_gotoh_band_gap_counts_match_a_plain_reference():
    """The kernel's gap events (runs of gap columns in one sequence, leading and trailing gaps included) and gap bases
    equal the plain implementation's traceback."""
    for it, q, t, lo, hi in _random_pairs(np.random.default_rng(16), 200):
        got = align._gotoh_band(q, t, lo, hi, *cst.SCORES)
        assert tuple(got[:3]) + tuple(got[4:]) == _gotoh_reference(q, t, lo, hi, *cst.SCORES), it


def _ns(seqs, orient, l5=None, l3=None):
    """The fields of a library the pair counts read: codes, offs, orient, l5, l3 (LTR record offsets, -1 = none)."""
    chunks = [codes(x) for x in seqs]
    n = len(seqs)
    return SimpleNamespace(codes=np.concatenate(chunks), offs=np.concatenate(([0], np.cumsum([c.size for c in chunks]))),
                           orient=np.array(orient, np.int8),
                           l5=np.full((n, 2), -1, np.int64) if l5 is None else np.array(l5, np.int64),
                           l3=np.full((n, 2), -1, np.int64) if l3 is None else np.array(l3, np.int64))


def _rc(x):
    return x.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def test_pair_counts_unknown_strand_flips_orientation():
    rng = np.random.default_rng(2)
    ltr5, inner, ltr3 = ("".join(rng.choice(list("ACGT"), n)) for n in (120, 600, 120))
    a = ltr5 + inner + ltr3
    lib = _ns([a, _rc(a)], [0, 0], l5=[[0, 120], [0, 120]], l3=[[720, 840], [720, 840]])
    assert [f.tolist() for f in align.pairs(lib, [[0, 1]], mode="ltr", threads=1)[:3]] == [[240], [0], [0]]


def test_whole_counts_orientation():
    """Whole-element counts: a reverse-complemented record with strand '-' aligns as is (records are in element
    orientation); with unknown strand the better orientation is chosen."""
    rng = np.random.default_rng(5)
    q = "".join(rng.choice(list("ACGT"), 800))
    t = sg.mutate(rng, q, 8, 4)
    lib = _ns([q, t, _rc(t)], [1, -1, 0])
    c = align.pairs(lib, [[0, 1], [0, 2]], mode="whole", threads=1)
    assert c.sites.tolist() == [800, 800] and c.ts.tolist() == [8, 8] and c.tv.tolist() == [4, 4]


def test_orient():
    """+1 keeps y; -1 when the other strand (by default y's reverse complement) is strictly closer by edit distance."""
    rng = np.random.default_rng(9)
    x = "".join(rng.choice(list("ACGT"), 200))
    y = sg.mutate(rng, x, 4, 2)
    assert align.orient(codes(x), codes(y)) == 1 and align.orient(codes(x), codes(_rc(y))) == -1
    assert align.orient(codes(x), codes(y), codes(_rc(y))) == 1
    assert align.orient(codes("ACGT"), codes("ACGT")) == 1                   # a palindrome: a tie keeps y


def _frozen_pairs(z, prefix):
    """The (x, y) code arrays of frozen sequence pairs (letters, record offsets)."""
    out = []
    for side in ("x", "y"):
        seq, offs = z[f"{prefix}_{side}_seq"], z[f"{prefix}_{side}_offs"]
        out.append([encode(seq[offs[i]:offs[i + 1]].tobytes()) for i in range(offs.size - 1)])
    return list(zip(*out))


def test_orient_ltrmix_rule():
    """LTRmix's rule as orient's keywords: edlib infix (HW) distances of x's first k bases (k = min(len x, len y,
    head)) against the first k + slack of y and of its reverse complement; a tie keeps y. LTRmix's own choices on
    60 frozen pairs (reverse complements, indels, unrelated pairs, sequences longer than the head)."""
    with np.load(GOLDEN / "ltrmix_route.npz") as z:
        pairs, flips = _frozen_pairs(z, "or"), z["or_flip"].tolist()
    got = [align.orient(x, y, mode="HW", head=cst.ORIENT_HEAD, slack=cst.ORIENT_SLACK) == -1 for x, y in pairs]
    assert got == flips and 0 < sum(flips) < len(flips) and max(x.size for x, _ in pairs) > cst.ORIENT_HEAD


def test_unit_cost_path_and_position_map():
    """edlib's unit-cost global path as runs (lengths; 0 aligned column, 1 query base against a gap, 2 target base
    against a gap) and the target position of each query base (-1 at a gap): LTRmix's CIGAR map on 80 frozen pairs; an
    empty sequence aligns as one gap run."""
    path = align.unit_cost_path(codes("AACCGGTTA"), codes("AACCTGGTA"))   # a T inserted in t, a T deleted from it
    assert path.dtype == np.int64 and path.tolist() == [[4, 1, 3, 1, 1], [0, 2, 0, 1, 0]]
    assert align.position_map(path).tolist() == [0, 1, 2, 3, 5, 6, 7, -1, 8]
    with np.load(GOLDEN / "ltrmix_route.npz") as z:
        pairs, want = _frozen_pairs(z, "cg"), z["cg_map"]
    got = np.concatenate([align.position_map(align.unit_cost_path(x, y)) for x, y in pairs])
    assert got.dtype == np.int64 and np.array_equal(got, want)
    for q, t in (("", "ACGT"), ("ACG", ""), ("", "")):
        p = align.unit_cost_path(codes(q), codes(t))
        assert p.shape[0] == 2 and int(p[0].sum()) == len(q) + len(t)
        assert align.position_map(p).tolist() == [-1] * len(q)


def _mixed_lib(rng):
    """40 copies of one element with mutated LTRs and internal regions: every fifth stored reverse-complemented with
    its strand unknown, every seventh (from the fourth) without LTR bounds."""
    seqs, l5, l3, orient = [], [], [], []
    base5, inner = "".join(rng.choice(list("ACGT"), 150)), "".join(rng.choice(list("ACGT"), 400))
    for k in range(40):
        rec = sg.mutate(rng, base5, 3, 1) + sg.mutate(rng, inner, 6, 3) + sg.mutate(rng, base5, 2, 2)
        flip = k % 5 == 4                                  # stored reverse-complemented, strand unknown
        seqs.append(_rc(rec) if flip else rec)
        orient.append(0 if flip else 1)
        has = k % 7 != 3                                   # some copies without LTR bounds
        l5.append([0, 150] if has else [-1, -1])
        l3.append([len(rec) - 150, len(rec)] if has else [-1, -1])
    return _ns(seqs, orient, l5, l3)


@pytest.mark.parametrize("threads", [1, 3])
def test_pairs_modes(threads):
    """mode ltr: 5'-5' plus 3'-3' LTR counts, whole records when a copy has no LTR bounds; whole: whole records;
    self: each copy's own 5' against its 3' LTR. Unknown strands take the closer orientation. Rows in input order,
    whatever the threads and chunks."""
    rng = np.random.default_rng(4)
    lib = _mixed_lib(rng)
    P = rng.integers(0, 40, (300, 2))
    for mode in ("ltr", "whole"):
        c = align.pairs(lib, P, mode=mode, threads=threads)
        assert len(c) == 5 and all(f.shape == (300,) and f.dtype == np.int64 for f in c)
        for n_, (a, b) in enumerate(P[:25]):
            ra, rb = (lib.codes[lib.offs[i]:lib.offs[i + 1]] for i in (a, b))
            if mode == "ltr" and lib.l5[a, 0] >= 0 and lib.l5[b, 0] >= 0:
                x5, x3 = ra[lib.l5[a, 0]:lib.l5[a, 1]], ra[lib.l3[a, 0]:lib.l3[a, 1]]
                y5, y3 = rb[lib.l5[b, 0]:lib.l5[b, 1]], rb[lib.l3[b, 0]:lib.l3[b, 1]]
                if (lib.orient[a] == 0 or lib.orient[b] == 0) and align.orient(x5, y5, align._rc(y3)) < 0:
                    y5, y3 = align._rc(y3), align._rc(y5)
                want = tuple(u + v for u, v in zip(align.global_counts(x5, y5), align.global_counts(x3, y3)))
            else:
                if (lib.orient[a] == 0 or lib.orient[b] == 0) and align.orient(ra, rb) < 0:
                    rb = align._rc(rb)
                want = align.global_counts(ra, rb)
            assert tuple(int(f[n_]) for f in c) == want, (mode, n_)
    own = np.flatnonzero(lib.l5[:, 0] >= 0)
    s = align.pairs(lib, own, mode="self", threads=threads)
    for n_, i in enumerate(own[:10]):
        r = lib.codes[lib.offs[i]:lib.offs[i + 1]]
        assert tuple(int(f[n_]) for f in s) == align.global_counts(r[lib.l5[i, 0]:lib.l5[i, 1]], r[lib.l3[i, 0]:lib.l3[i, 1]])
    empty = align.pairs(lib, np.zeros((0, 2), np.int64), mode="whole", threads=threads)
    assert all(f.shape == (0,) for f in empty)
    with pytest.raises(ValueError, match="mode"):
        align.pairs(lib, P, mode="other", threads=threads)


def test_pairs_align_every_chunk_in_a_worker(monkeypatch):
    """Several chunks and threads: this process aligns no item (the kernels compile here on a short slice of the
    library, before the fork), every chunk runs in a worker, and the counts equal one thread's."""
    rng = np.random.default_rng(5)
    seqs = ["".join(rng.choice(list("ACGT"), 300)) for _ in range(30)]
    lib = _ns(seqs, [0 if k % 3 == 0 else 1 for k in range(30)])
    P = rng.integers(0, 30, (3 * cst.ALIGN_CHUNK["whole"], 2))
    want = align.pairs(lib, P, mode="whole", threads=1)
    here, counts = [], align._counts
    parent = os.getpid()

    def spy(lib_, x, mode):
        if os.getpid() == parent:
            here.append(x)
        return counts(lib_, x, mode)
    monkeypatch.setattr(align, "_counts", spy)
    got = align.pairs(lib, P, mode="whole", threads=3)
    assert not here
    assert all(np.array_equal(g, w) for g, w in zip(got, want))


PROBE = """import sys
import numpy as np
from types import SimpleNamespace
from numba.core.dispatcher import Dispatcher
from ltrtrace.seq import align
z = np.load(sys.argv[1])
lib = SimpleNamespace(**{k: z[k] for k in z.files})
def signatures():
    return {f"{m.__name__}.{n}": len(f.signatures) for m in list(sys.modules.values())
            if getattr(m, "__name__", "").startswith("ltrtrace") for n, f in vars(m).items()
            if isinstance(f, Dispatcher)}
align._compile(lib)
print(signatures())
for mode in ("ltr", "whole"):
    align.pairs(lib, lib.P, mode=mode, threads=1)
align.pairs(lib, np.flatnonzero(lib.l5[:, 0] >= 0), mode="self", threads=1)
print(signatures())
"""


def test_compile_gives_every_signature_the_workers_use(tmp_path):
    """numba compiles one specialization per argument layout (a record slice, a reverse-complement view), and
    signatures are per process: in a fresh process, after the compile step alone, the counts of all three modes add
    no signature to any kernel, so no forked worker compiles."""
    import ast
    import subprocess
    import sys
    rng = np.random.default_rng(4)
    lib = _mixed_lib(rng)
    np.savez(tmp_path / "lib.npz", codes=lib.codes, offs=lib.offs, orient=lib.orient, l5=lib.l5, l3=lib.l3,
             P=rng.integers(0, 40, (300, 2)))
    (tmp_path / "probe.py").write_text(PROBE)
    r = subprocess.run([sys.executable, str(tmp_path / "probe.py"), str(tmp_path / "lib.npz")], check=False,
                       capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stderr[-2000:]
    compiled, used = (ast.literal_eval(line) for line in r.stdout.splitlines()[-2:])
    assert sum(compiled.values()) > 0 and used == compiled
