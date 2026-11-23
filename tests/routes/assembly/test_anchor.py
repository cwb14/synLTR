"""Anchoring (routes.assembly.anchor): the flank-mm2 route (flank pieces, the minimap2 index, piece mapping), the
PAF route and `anchor_pair` with its cache."""
from __future__ import annotations

import contextlib
import gzip
import random
import shutil
from itertools import pairwise
from typing import ClassVar

import builders
import numpy as np
import pytest
from builders import make_element, mutate, pair_context, rand_seq, read_calls, write_fasta
from markers import needs_lastz

import synltr.routes.assembly.anchor as anchor_mod
from synltr.base.checkpoint import code_files
from synltr.base.constants import FLANK_PIECE_MIN, SPLIT_INDEL
from synltr.io.elements import Element
from synltr.io.samples import Sample
from synltr.io.seqio import Fasta, n_runs
from synltr.routes.assembly.anchor import FlankMm2Route, PafRoute, PairContext, Piece, anchor_pair, flank_pieces
from synltr.routes.assembly.blocks import BLOCK_DTYPE
from synltr.routes.assembly.pair import genotype_pair


def el(start, end, level=0):
    """Element at 1-based inclusive start..end on chr1."""
    return Element(genome="g", eid=f"g|chr1:{start}-{end}", name="e", seqid="chr1", start=start, end=end,
                   strand="+", ltr_left=(start, start + 9), ltr_right=(end - 9, end), tsd=None, k=None,
                   family=None, superfamily=None, clade=None, confidence=None, level=level)


@pytest.fixture
def fa_n(tmp_path):
    """30 kb chr1 of random bases with an N run at 9000-9500 and one at 13000-13100."""
    s = list(builders.rand_seq(30_000, random.Random(3)))
    s[9000:9500] = "N" * 500
    s[13000:13100] = "N" * 100
    with Fasta(builders.write_fasta(tmp_path / "g.fa", {"chr1": "".join(s)})) as fa:
        yield fa


# --- flank pieces -----------------------------------------------------------------------------------------


def test_pieces_skip_elements_and_ns(fa_n):
    els = [el(10_001, 12_000), el(10_501, 11_000, level=1), el(12_501, 12_800)]  # child; neighbor
    pcs = flank_pieces(els, fa_n, cap=5000, target=100_000)  # target above the supply: every free piece
    assert pcs and pcs == sorted(pcs, key=lambda p: (p.start, p.end))
    spans = [(e.start - 1, e.end) for e in els]
    ns = n_runs(fa_n.fetch("chr1", 0, 30_000))
    for p in pcs:
        assert p.end - p.start >= FLANK_PIECE_MIN and p.seqid == "chr1"
        assert all(p.end <= s or p.start >= t for s, t in spans + ns), p
    assert len(set(pcs)) == len(pcs)
    assert all(p.start >= 10_000 - 5000 and p.end <= 12_800 + 5000 for p in pcs)  # within cap of an element


def test_pieces_nearest_first_and_capped(fa_n):
    pcs = flank_pieces([el(20_001, 21_000)], fa_n, cap=8000, target=1500)
    left = [p for p in pcs if p.end <= 20_000]
    right = [p for p in pcs if p.start >= 21_000]
    assert (left[-1].end, right[0].start) == (20_000, 21_000)  # nearest piece abuts the element
    for side in (left, right):
        assert 1500 <= sum(p.end - p.start for p in side)
    assert sum(p.end - p.start for p in left[1:]) < 1500  # no piece taken after the target was reached
    capped = flank_pieces([el(20_001, 21_000)], fa_n, cap=1000, target=5000)
    assert min(p.start for p in capped) >= 19_000 and max(p.end for p in capped) <= 22_000
    assert isinstance(capped[0], Piece)


def _naive_pieces(elements, fasta, cap, target):
    """A plain quadratic implementation of flank_pieces (reference)."""
    lengths = fasta.lengths()
    out = set()
    for e in elements:
        if e.level != 0:
            continue
        s0, t0, n = e.start - 1, e.end, lengths[e.seqid]
        for lo, hi, left in ((max(0, s0 - cap), s0, True), (t0, min(n, t0 + cap), False)):
            if hi <= lo:
                continue
            cuts = [(max(s, lo) - lo, min(t, hi) - lo) for s, t in ((x.start - 1, x.end) for x in elements
                                                                    if x.seqid == e.seqid) if s < hi and t > lo]
            cuts += n_runs(fasta.fetch(e.seqid, lo, hi))
            free, pos = [], 0
            for s, t in sorted(cuts):
                if s - pos >= FLANK_PIECE_MIN:
                    free.append((pos, s))
                pos = max(pos, t)
            if hi - lo - pos >= FLANK_PIECE_MIN:
                free.append((pos, hi - lo))
            total = 0
            for s, t in (reversed(free) if left else free):
                if total >= target:
                    break
                out.add(Piece(e.seqid, lo + s, lo + t))
                total += t - s
    return sorted(out, key=lambda p: (p.start, p.end))


def _random_elements(rng, n, length, nested=0.2):
    els = []
    for _ in range(n):
        s = rng.randrange(1, length - 3000)
        e = s + rng.randrange(50, 2500)
        els.append(el(s, e, level=0))
    for _ in range(int(n * nested)):  # nested children inside a random parent
        p = rng.choice(els)
        if p.end - p.start > 40:
            s = rng.randrange(p.start, p.end - 20)
            els.append(el(s, rng.randrange(s + 10, p.end + 1), level=1))
    return els


def test_pieces_match_reference_on_random_case(fa_n):
    rng = random.Random(11)
    els = _random_elements(rng, 60, 30_000, nested=0.5)
    for cap, target in ((3000, 1000), (10_000, 2000)):
        assert flank_pieces(els, fa_n, cap, target) == _naive_pieces(els, fa_n, cap, target)


def test_pieces_many_elements(tmp_path):
    """20,000 elements: the bisect span lookup keeps flank_pieces near linear (a few seconds; a scan of every
    span per element side would take about 100 s). The genome is drawn with random.choices, so the fixture does
    not dominate the test's time."""
    import time
    rng = random.Random(2)
    n_bp = 40_000_000
    with open(tmp_path / "big.fa", "w") as fh:
        fh.write(">chr1\n")
        fh.writelines("".join(rng.choices("ACGT", k=1_000_000)) + "\n" for i in range(0, n_bp, 1_000_000))
    els = [el(s, s + 1500) for s in range(2000, n_bp - 5000, 2000)]  # 20,000 elements, 500-bp gaps
    with Fasta(tmp_path / "big.fa") as fa:
        t0 = time.perf_counter()
        pcs = flank_pieces(els, fa)
        assert time.perf_counter() - t0 < 30
    assert pcs and len({(p.start, p.end) for p in pcs}) == len(pcs)


def test_pieces_at_sequence_ends(fa_n):
    """An element at either end of its sequence has a flank on its inner side only, clipped to the sequence; a
    flank of N runs and free stretches under FLANK_PIECE_MIN yields no piece."""
    first, last = el(1, 1000), el(29_001, 30_000)
    pcs = flank_pieces([first, last], fa_n, cap=2000, target=2000)
    assert [(p.start, p.end) for p in pcs] == [(1000, 3000), (27_000, 29_000)]
    # left: the N run at 9000-9500 fills the 200-bp cap; right: 100 free bp, under FLANK_PIECE_MIN, then N
    assert flank_pieces([el(9_501, 12_900)], fa_n, cap=200, target=200) == []


# --- the flank-mm2 route ----------------------------------------------------------------------------------


@pytest.fixture(scope="module")
def tpair():
    return builders.build_pair(random.Random(7))


def _states(calls, pair):
    names = {e.eid: e.name for e in pair.elements}
    return {names[k]: calls[k]["state"] for k in pair.expect}, {names[k]: v for k, v in pair.expect.items()}


@needs_lastz
def test_flank_mm2_resolves_toy(tpair, tmp_path):
    with pair_context(tmp_path, tpair.a, tpair.b, tpair.elements) as ctx:
        calls = read_calls(genotype_pair(ctx, FlankMm2Route()))
    got, want = _states(calls, tpair)
    for name, st in want.items():
        if name.split(".")[0] in ("P", "E", "D_left", "D_both"):
            assert got[name] == st, (name, calls)


def test_flank_mm2_cache_differs_by_cap(tpair, tmp_path):
    with pair_context(tmp_path, tpair.a, tpair.b, tpair.elements) as ctx:
        a, b = FlankMm2Route(cap=2000), FlankMm2Route(cap=5000)
        assert a.name == b.name == "flank-mm2" and a.reads_elements
        assert a.params(ctx) != b.params(ctx)
        anchor_pair(ctx, a)
        anchor_pair(ctx, b)


def test_map_pieces_bounds_in_flight(tpair, tmp_path, monkeypatch):
    from concurrent.futures import ThreadPoolExecutor

    from synltr.base.constants import PIECES_QUEUED
    peak = []

    class Spy(ThreadPoolExecutor):
        live: ClassVar[list] = []

        def submit(self, fn, *args):
            fut = super().submit(fn, *args)
            self.live.append(fut)
            peak.append(sum(not f.done() for f in self.live))
            return fut

    monkeypatch.setattr(anchor_mod, "ThreadPoolExecutor", Spy)
    with pair_context(tmp_path, tpair.a, tpair.b, tpair.elements, threads=2) as ctx:
        pcs = flank_pieces(ctx.a_elements, ctx.a_fa)
        raw = anchor_mod.map_pieces(ctx, pcs)
        assert len(pcs) > 4 * 2 * PIECES_QUEUED and len(raw)
    assert max(peak) <= PIECES_QUEUED * 2


# =100 X2 =50 I5 =300 | I4800 | =1000 | D90 | =500 from A[10,6767) to B[1000,3042): the 5-bp insertion stays
# inside a piece, the long indels split
_CIGAR = [[100, 7], [2, 8], [50, 7], [5, 1], [300, 7], [4800, 1], [1000, 7], [90, 2], [500, 7]]
_PIECES = {  # (a0, a1, b0, b1, matches, alen) in A order; on '-' the CIGAR runs along B with A descending
    1: [(10, 467, 1000, 1452, 450, 457), (5267, 6267, 1452, 2452, 1000, 1000),
        (6267, 6767, 2542, 3042, 500, 500)],
    -1: [(10, 510, 2542, 3042, 500, 500), (510, 1510, 1452, 2452, 1000, 1000),
         (6310, 6767, 1000, 1452, 450, 457)],
}


@pytest.mark.parametrize("strand", [1, -1])
def test_pieces_split_at_long_indels(strand):
    assert 5 < SPLIT_INDEL <= 90
    assert anchor_mod._pieces(_CIGAR, 10, 1000, 3042, strand) == _PIECES[strand]
    assert anchor_mod._pieces([[700, 7], [3, 8]], 0, 0, 703, strand) == [(0, 703, 0, 703, 700, 703)]
    with pytest.raises(ValueError, match="CIGAR"):
        anchor_mod._pieces([[700, 0]], 0, 0, 700, strand)  # M: matches cannot be counted


# --- B's minimap2 index -----------------------------------------------------------------------------------


@pytest.fixture(scope="module")
def toy(tmp_path_factory):
    """300 kb of shared backbone (chr1 200 kb, chr2 100 kb). A carries 10 A-only elements (6 + 4, 4.8 kb
    each); B holds the backbones mutated 1%, in the other order, chr2's copy reverse-complemented.

    Returns the A and B FASTA paths and the shared (backbone) A intervals per A seqid.
    """
    from synltr.io.seqio import revcomp

    rng = random.Random(15)
    a_recs, b_seqs, shared = {}, [], {}
    for name, size, n_el in (("chr1", 200_000, 6), ("chr2", 100_000, 4)):
        bb = rand_seq(size, rng)
        cuts = sorted(rng.sample(range(5_000, size - 5_000, 1000), n_el))
        parts, ivs, prev, pos = [], [], 0, 0
        for c in cuts:
            elem, _, _ = make_element(rng, 400, 4000)
            parts += [bb[prev:c], elem]
            ivs.append((pos, pos + c - prev))
            pos += c - prev + len(elem)
            prev = c
        parts.append(bb[prev:])
        ivs.append((pos, pos + size - prev))
        a_recs[name] = "".join(parts)
        shared[name] = ivs
        b_seqs.append(mutate(bb, 0.01, rng))
    d = tmp_path_factory.mktemp("toy")
    a = write_fasta(d / "A.fa", a_recs)
    b = write_fasta(d / "B.fa", {"ctgX": revcomp(b_seqs[1]), "ctgY": b_seqs[0]})
    return a, b, shared


@pytest.fixture
def open_ctx(tmp_path):
    """make(a_label, a_fasta, b_label, b_fasta, outdir=None, elements=()) -> PairContext (handles closed at
    teardown; outdir defaults to tmp_path / "out")."""
    with contextlib.ExitStack() as stack:
        def make(a, a_path, b, b_path, outdir=None, elements=()):
            return PairContext(
                a=Sample(a, a_path, None, None, None), b=Sample(b, b_path, None, None, None),
                a_fa=stack.enter_context(Fasta(a_path)), b_fa=stack.enter_context(Fasta(b_path)),
                a_elements=list(elements), b_elements=None, outdir=outdir or tmp_path / "out", threads=2)
        yield make


def _toy_elements(shared):
    """The toy's A-only elements: the A spans between its shared intervals, 1-based closed."""
    out = []
    for seqid, ivs in shared.items():
        for (_, e0), (s1, _) in pairwise(ivs):
            out.append(Element("A", f"A|{seqid}:{e0 + 1}-{s1}", "e", seqid, e0 + 1, s1, "+", (e0 + 1, e0 + 400),
                               (s1 - 399, s1), None, None, None, None, None, None))
    return out


def test_index_built_once_with_or_without_minimap2(toy, open_ctx, tmp_path, monkeypatch):
    a, b, shared = toy
    route = FlankMm2Route()
    elements = _toy_elements(shared)
    ctx = open_ctx("A", a, "B", b, outdir=tmp_path / "one", elements=elements)
    first = route.raw_blocks(ctx)
    assert len(first) and (first["chain"] == -1).all()
    mmi = ctx.workdir / "anchors" / "B.asm20.mmi"
    built = mmi.stat().st_mtime_ns
    assert route.raw_blocks(ctx).tobytes() == first.tobytes()
    assert mmi.stat().st_mtime_ns == built  # reused, not rebuilt
    monkeypatch.setattr(shutil, "which", lambda *args, **kw: None)  # no minimap2 binary: mappy indexes B
    ctx2 = open_ctx("A", a, "B", b, outdir=tmp_path / "two", elements=elements)
    assert route.raw_blocks(ctx2).tobytes() == first.tobytes()
    assert (ctx2.workdir / "anchors" / "B.asm20.mmi").read_bytes() == mmi.read_bytes()


def test_index_key_holds_only_the_index_code():
    """B's index is keyed on the code that builds it: an edit to the rest of the route (chaining, classification)
    never rebuilds an index, while an edit to the anchoring or the pair does reach their own stages."""
    index = code_files("routes.assembly.index")
    assert "routes/assembly/index.py" in index
    assert not [f for f in index if f.startswith(("routes/", "seq/", "io/", "calls/")) and
                f != "routes/assembly/index.py"]
    assert {"routes/assembly/chain.py", "routes/assembly/index.py"} <= set(code_files("routes.assembly.anchor"))
    assert {"routes/assembly/classify.py", "seq/junction.py"} <= set(code_files("routes.assembly.pair"))


# --- PAF route --------------------------------------------------------------------------------------------

_PAF = [  # A (a1) is the query
    "a1\t20000\t1000\t9000\t+\tb1\t30000\t2000\t10000\t7900\t8000\t60",
    "a1\t20000\t10000\t19000\t-\tb2\t15000\t1000\t10000\t8900\t9000\t60",
]


def _swapped(line):
    """The same PAF record with query and target exchanged."""
    f = line.split("\t")
    return "\t".join(f[5:9] + [f[4]] + f[0:4] + f[9:])


def test_paf_route(tmp_path, open_ctx):
    rng = random.Random(3)
    a = write_fasta(tmp_path / "A.fa", {"a1": rand_seq(20_000, rng)})
    b = write_fasta(tmp_path / "B.fa", {"b1": rand_seq(30_000, rng), "b2": rand_seq(15_000, rng)})
    fwd, rev = tmp_path / "fwd", tmp_path / "rev"
    fwd.mkdir()
    rev.mkdir()
    (fwd / "A__B.paf").write_text("\n".join(_PAF) + "\n")
    with gzip.open(rev / "B__A.paf.gz", "wt") as fh:
        fh.write("\n".join(map(_swapped, _PAF)) + "\n")

    ctx = open_ctx("A", a, "B", b)
    got = PafRoute(fwd).raw_blocks(ctx)
    cols = ["a_chrom", "a_start", "a_end", "b_chrom", "b_start", "b_end", "strand", "matches", "alen"]
    assert got[cols].tolist() == [(0, 1000, 9000, 0, 2000, 10000, 1, 7900, 8000),
                                  (0, 10000, 19000, 1, 1000, 10000, -1, 8900, 9000)]
    assert PafRoute(rev).raw_blocks(ctx).tobytes() == got.tobytes()  # B__A.paf.gz read with swap

    ba = open_ctx("B", b, "A", a)  # B as A: A__B.paf read with swap
    bl, an, bn, d = anchor_pair(ba, PafRoute(fwd))
    assert (an, bn) == (["b1", "b2"], ["a1"])
    assert bl[cols[:7]].tolist() == [(0, 2000, 10000, 0, 1000, 9000, 1),
                                     (1, 1000, 10000, 0, 10000, 19000, -1)]
    assert d == pytest.approx(1 - (7900 + 8900) / (8000 + 9000))
    with open(fwd / "A__B.paf", "a") as fh:  # the PAF changes: the netted cache is rebuilt
        fh.write("a1\t20000\t19000\t20000\t+\tb1\t30000\t20000\t21000\t1000\t1000\t60\n")
    assert len(anchor_pair(ba, PafRoute(fwd))[0]) == 3

    with pytest.raises(FileNotFoundError, match=r"A__B\.paf\[\.gz\].*B__A\.paf\[\.gz\]"):
        PafRoute(tmp_path / "none").raw_blocks(ctx)


# --- anchor_pair ------------------------------------------------------------------------------------------


class _FixedRoute:
    """Returns one fixed block (or none); counts its calls; reads elements."""

    name = "fixed"
    reads_elements = True

    def __init__(self, empty=False):
        self.calls, self.version, self.empty = 0, 1, empty

    def params(self, ctx):
        return {"version": self.version}

    def raw_blocks(self, ctx):
        self.calls += 1
        out = np.zeros(0 if self.empty else 1, dtype=BLOCK_DTYPE)
        if not self.empty:
            out[0] = (0, 0, 5000, 0, 0, 5000, 1, 4950, 5000, -1, 0.0, np.nan)
        return out


def test_anchor_pair_cache_key(tmp_path, open_ctx):
    rng = random.Random(4)
    a = write_fasta(tmp_path / "A.fa", {"a1": rand_seq(6000, rng)})
    b = write_fasta(tmp_path / "B.fa", {"b1": rand_seq(6000, rng)})
    b2 = write_fasta(tmp_path / "B2.fa", {"b1": rand_seq(6000, rng)})
    route = _FixedRoute()
    ctx = open_ctx("A", a, "B", b)
    bl, an, bn, d = anchor_pair(ctx, route)
    assert route.calls == 1 and len(bl) == 1 and d == pytest.approx(0.01)
    ctx.threads = 8  # not part of the key
    again = anchor_pair(ctx, route)
    assert route.calls == 1 and again[0].tobytes() == bl.tobytes() and again[1:] == (an, bn, d)
    ctx.a_elements = [Element("A", "A|a1:101-600", "e", "a1", 101, 600, "+", (101, 150), (551, 600),
                              None, None, None, None, None, None)]
    anchor_pair(ctx, route)  # elements changed, and this route reads them
    assert route.calls == 2
    route.version = 2  # a route parameter changed
    anchor_pair(ctx, route)
    assert route.calls == 3
    anchor_pair(open_ctx("A", a, "B", b2), route)  # another B FASTA under the same label
    assert route.calls == 4

    bl, an, bn, d = anchor_pair(open_ctx("A", a, "B", b, outdir=tmp_path / "empty"), _FixedRoute(empty=True))
    assert len(bl) == 0 and (an, bn) == (["a1"], ["b1"]) and np.isnan(d)


def test_anchors_from_the_flat_net_are_recomputed(tmp_path, open_ctx, monkeypatch):
    # a checkpoint whose key's chain parameters lack "net" (a flat-net checkpoint) is never reused
    rng = random.Random(4)
    a = write_fasta(tmp_path / "A.fa", {"a1": rand_seq(6000, rng)})
    b = write_fasta(tmp_path / "B.fa", {"b1": rand_seq(6000, rng)})
    route, ctx = _FixedRoute(), open_ctx("A", a, "B", b)
    assert anchor_mod._CHAIN_PARAMS["net"] == "syntenic"
    with monkeypatch.context() as mp:
        mp.setattr(anchor_mod, "_CHAIN_PARAMS", {k: v for k, v in anchor_mod._CHAIN_PARAMS.items() if k != "net"})
        anchor_pair(ctx, route)  # the old checkpoint
    assert route.calls == 1
    anchor_pair(ctx, route)
    assert route.calls == 2  # recomputed
    anchor_pair(ctx, route)
    assert route.calls == 2  # the new checkpoint is reused
