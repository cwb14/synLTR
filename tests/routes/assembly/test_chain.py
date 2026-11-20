"""Chaining DP (chain_blocks) and the syntenic 1:1 net (net_blocks) of routes.assembly.chain."""
import logging
import math
import time
from collections import Counter, defaultdict

import numpy as np
import pytest

from ltrtrace.base.constants import ANCHOR_MIN, CHAIN_INDEL_COST, CHAIN_MAX_GAP, CHAIN_OVERLAP, MIN_CHAIN_SCORE
from ltrtrace.routes.assembly import chain
from ltrtrace.routes.assembly.blocks import BLOCK_DTYPE, a_to_b
from ltrtrace.routes.assembly.chain import chain_blocks, net_blocks

K = 1000
LOG = "ltrtrace.routes.assembly.chain"


def _blocks(rows):
    """Blocks from (a_chrom, a_start, a_end, b_chrom, b_start, b_end, strand, matches); alen = A length."""
    out = np.zeros(len(rows), dtype=BLOCK_DTYPE)
    for i, (ac, a0, a1, bc, b0, b1, st, m) in enumerate(rows):
        out[i] = (ac, a0, a1, bc, b0, b1, st, m, a1 - a0, -1, 0.0, np.nan)
    return out


def _spans(out):
    """Sorted (a_start, a_end, b_start, b_end) of every row."""
    return sorted((int(r["a_start"]), int(r["a_end"]), int(r["b_start"]), int(r["b_end"])) for r in out)


def _random_blocks(n, seed, n_chrom=2, frac_ortho=0.7):
    """Seeded synthetic blocks, shuffled: collinear runs A chrom c -> B chrom c (odd c inverted, B descending)
    plus random paralogous hits anywhere. Block lengths 20-5000 bp, B length within 5% of A."""
    rng = np.random.default_rng(seed)
    a_len = rng.integers(20, 5000, n)  # ~1.6% shorter than ANCHOR_MIN
    b_len = np.maximum(1, a_len + rng.integers(-(a_len // 20), a_len // 20 + 1))
    a_chrom, b_chrom = np.zeros(n, np.int32), np.zeros(n, np.int32)
    a_start, b_start, strand = np.zeros(n, np.int64), np.zeros(n, np.int64), np.ones(n, np.int8)
    n_orth = int(n * frac_ortho)
    span = 0
    for c in range(n_chrom):
        idx = np.arange(c, n_orth, n_chrom)
        a_start[idx] = np.cumsum(a_len[idx] + rng.exponential(500, len(idx)).astype(np.int64)) - a_len[idx]
        bpos = np.cumsum(b_len[idx] + rng.exponential(500, len(idx)).astype(np.int64)) - b_len[idx]
        top = int(bpos[-1] + b_len[idx][-1])
        b_start[idx] = top - bpos - b_len[idx] if c % 2 else bpos
        strand[idx] = -1 if c % 2 else 1
        a_chrom[idx] = b_chrom[idx] = c
        span = max(span, int(a_start[idx][-1] + a_len[idx][-1]), top)
    par = np.arange(n_orth, n)
    a_chrom[par], b_chrom[par] = rng.integers(0, n_chrom, len(par)), rng.integers(0, n_chrom, len(par))
    a_start[par], b_start[par] = rng.integers(0, span, len(par)), rng.integers(0, span, len(par))
    strand[par] = rng.choice(np.array([-1, 1], np.int8), len(par))
    ident = np.where(np.arange(n) < n_orth, rng.uniform(0.9, 1.0, n), rng.uniform(0.7, 1.0, n))
    out = np.zeros(n, dtype=BLOCK_DTYPE)
    out["a_chrom"], out["a_start"], out["a_end"] = a_chrom, a_start, a_start + a_len
    out["b_chrom"], out["b_start"], out["b_end"] = b_chrom, b_start, b_start + b_len
    out["strand"], out["matches"], out["alen"] = strand, (a_len * ident).astype(np.int32), np.maximum(a_len, b_len)
    out["chain"], out["margin"] = -1, np.nan
    return out[rng.permutation(n)]


def _genome_blocks(seed, n_chrom=3, n_sites=150):
    """Seeded synthetic blocks of a genome pair, shuffled. Per chromosome pair c (odd c inverted in B): a
    collinear run of 2-20 kb orthologous blocks, each followed by an insertion site (none, A-only, B-only or
    both; 0.5-30 kb). Then blocks at 75-95% identity, either strand: fills in place (A and B inside one site's
    insertions), paralogs (an A insertion paired with a B insertion elsewhere), leaks (B past every run),
    two-block chains with one block in an A insertion near a run's end and one past it (top level), and an
    unplaced A scaffold (A sequence n_chrom) paired with B insertions or with an unplaced B scaffold."""
    rng = np.random.default_rng(seed)
    rows, sites, ends = [], [], {}
    for c in range(n_chrom):
        a = b = 0
        run, gaps = [], []
        for _ in range(n_sites):
            n = int(rng.integers(2000, 20000))
            run.append((a, a + n, b, b + n, int(n * rng.uniform(0.95, 1.0))))
            a, b = a + n, b + n
            kind = int(rng.integers(4))  # no insertion, A-only, B-only, both
            ia = int(rng.integers(500, 30000)) if kind in (1, 3) else int(rng.integers(0, 50))
            ib = int(rng.integers(500, 30000)) if kind in (2, 3) else int(rng.integers(0, 50))
            gaps.append((a, a + ia, b, b + ib))
            a, b = a + ia, b + ib
        flip, ends[c] = c % 2 == 1, (a, b)
        for a0, a1, b0, b1, m in run:
            rows.append((c, a0, a1, c, *((b - b1, b - b0) if flip else (b0, b1)), -1 if flip else 1, m))
        sites += [(c, a0, a1, *((b - b1, b - b0) if flip else (b0, b1))) for a0, a1, b0, b1 in gaps]

    def inside(lo, hi, n):
        s = int(rng.integers(lo, hi - n + 1))
        return s, s + n

    def pick(cands):
        return cands[int(rng.integers(len(cands)))]

    def add(ac, a, bc, b, strand, ident=None):
        m = int((a[1] - a[0]) * (ident or rng.uniform(0.75, 0.95)))
        rows.append((ac, a[0], a[1], bc, b[0], b[1], strand, m))

    def strand():
        return int(rng.choice([-1, 1]))

    big_a = [s for s in sites if s[2] - s[1] >= 1500]
    big_b = [s for s in sites if s[4] - s[3] >= 1500]
    b_end = max(e[1] for e in ends.values())
    for _ in range(40):  # in place
        c, a0, a1, b0, b1 = pick([s for s in big_a if s[4] - s[3] >= 1500])
        n = int(rng.integers(1200, min(a1 - a0, b1 - b0) + 1))
        add(c, inside(a0, a1, n), c, inside(b0, b1, n), strand())
    for _ in range(80):  # paralogs
        sa, sb = pick(big_a), pick(big_b)
        n = int(rng.integers(1200, min(sa[2] - sa[1], sb[4] - sb[3]) + 1))
        add(sa[0], inside(sa[1], sa[2], n), sb[0], inside(sb[3], sb[4], n), strand())
    for _ in range(20):  # leaks
        sa = pick(big_a)
        n = int(rng.integers(1200, sa[2] - sa[1] + 1))
        b0 = b_end + int(rng.integers(0, 10**6))
        add(sa[0], inside(sa[1], sa[2], n), int(rng.integers(n_chrom)), (b0, b0 + n), strand())
    for _ in range(20):  # one block in an A insertion near a run's end, one past it
        c = int(rng.integers(n_chrom))
        sa = pick([s for s in big_a if s[0] == c and s[1] > ends[c][0] - 300_000])
        n = int(rng.integers(1200, sa[2] - sa[1] + 1))
        b0, a2 = b_end + 2 * 10**6 + int(rng.integers(0, 10**7)), ends[c][0] + int(rng.integers(0, 100_000))
        add(c, inside(sa[1], sa[2], n), n_chrom + 1, (b0, b0 + n), 1, 0.97)
        add(c, (a2, a2 + 5000), n_chrom + 1, (b0 + n + 1000, b0 + n + 6000), 1, 0.97)
    for _ in range(20):  # the unplaced A scaffold
        sb = pick(big_b)
        n, a0 = int(rng.integers(1200, sb[4] - sb[3] + 1)), int(rng.integers(0, 10**6))
        if rng.random() < 0.5:
            add(n_chrom, (a0, a0 + n), sb[0], inside(sb[3], sb[4], n), strand())
        else:
            add(n_chrom, (a0, a0 + n), n_chrom + 2, (a0, a0 + n), strand())
    return _blocks(rows)[rng.permutation(len(rows))]


def _competing_blocks(n, seed):
    """Seeded single-block chains that compete but never nest: each block has an A or a B sequence of its own
    (A sequence 0 with B sequence 1 + i, or B sequence 0 with A sequence 1 + i), so nothing chains, and the
    blocks sharing A sequence 0 (or B sequence 0) overlap there at random. Lengths 500-5000 bp."""
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n):
        size, pos = int(rng.integers(500, 5000)), int(rng.integers(0, 200_000))
        own, m, st = int(rng.integers(0, 10**6)), int(size * rng.uniform(0.7, 1.0)), int(rng.choice([-1, 1]))
        if i % 2:
            rows.append((0, pos, pos + size, 1 + i, own, own + size, st, m))
        else:
            rows.append((1 + i, own, own + size, 0, pos, pos + size, st, m))
    return _blocks(rows)


def _reference_net(ch, syntenic=True, places=None):
    """`net_blocks` written plainly: the module's competition credit and trimming, then the synteny rule with
    each parent found by scanning every kept chain (`_brute_place`). syntenic=False: the flat 1:1 net, without
    the test. `places` (a Counter) collects each tested piece's (A, B) places, "chain" where its chain's
    top-level piece supports it."""
    srt = ch[np.lexsort((ch["a_end"], ch["a_start"], ch["chain"], -ch["chain_score"]))]
    claimed_a, claimed_b = defaultdict(chain._Claims), defaultdict(chain._Claims)
    kept, rows, comp = [], [], []
    for c in np.split(srt, np.flatnonzero(np.diff(srt["chain"]) != 0) + 1):
        own_a, own_b, pieces = defaultdict(chain._Claims), defaultdict(chain._Claims), []
        for blk in c:
            ac, a0, a1, bc, b0, b1, strand, matches, alen, cid, score, _ = blk.item()
            hits_a, hits_b = claimed_a[ac].overlapping(a0, a1), claimed_b[bc].overlapping(b0, b1)
            chain._credit(hits_a, a0, a1, score, comp)
            chain._credit(hits_b, b0, b1, score, comp)
            for p0, p1, q0, q1 in chain._trim(blk, hits_a, hits_b, own_a[ac], own_b[bc]):
                frac = (p1 - p0) / (a1 - a0)
                pieces.append((ac, p0, p1, bc, q0, q1, strand, round(matches * frac), round(alen * frac), cid,
                               score, np.nan))
        if syntenic and sum(p[7] for p in pieces) >= MIN_CHAIN_SCORE:
            pl = [(_brute_place(kept, p, (0, 1, 2), (3, 4, 5)), _brute_place(kept, p, (3, 4, 5), (0, 1, 2)))
                  for p in pieces]
            pl = [("chain" if a == "out" and any(x == "top" for x, _ in pl) else a,
                   "chain" if b == "out" and any(y == "top" for _, y in pl) else b) for a, b in pl]
            if places is not None:
                places.update(pl)
            pieces = [p for p, (a, b) in zip(pieces, pl) if "out" not in (a, b)]
        if sum(p[7] for p in pieces) < MIN_CHAIN_SCORE:
            continue
        for p in pieces:
            claimed_a[p[0]].add(p[1], p[2], len(rows))
            claimed_b[p[3]].add(p[4], p[5], len(rows))
            rows.append(p)
            comp.append(0.0)
        kept.append(pieces)
    out = np.array(rows, dtype=BLOCK_DTYPE)
    best, margin = np.array(comp), np.ones(len(rows))
    margin[best > 0] = 1 - best[best > 0] / out["chain_score"][best > 0]
    out["margin"] = margin
    return out[np.lexsort((out["a_start"], out["a_chrom"]))]


def _brute_place(kept, p, own, other):
    """p's place in the frame whose (sequence, start, end) fields are `own`: "top" if no kept chain's extent
    holds p's span, else, in the innermost such chain (the first kept on ties), "syn"/"inv" if p's span in the
    other genome lies inside the gap between the parent's pieces flanking p, else "out"."""
    seq, lo, hi = own
    holders = []
    for k in kept:
        if k[0][seq] == p[seq]:
            e0, e1 = min(x[lo] for x in k), max(x[hi] for x in k)
            if e0 <= p[lo] and p[hi] <= e1:
                holders.append((e1 - e0, len(holders), k))
    if not holders:
        return "top"
    par = min(holders)[2]
    oseq, olo, ohi = other
    if par[0][oseq] != p[oseq]:
        return "out"
    left = max((x for x in par if x[hi] <= p[lo]), key=lambda x: x[hi])
    right = min((x for x in par if x[lo] >= p[hi]), key=lambda x: x[lo])
    g0, g1 = (left[ohi], right[olo]) if par[0][6] > 0 else (right[ohi], left[olo])
    if g0 <= p[olo] and p[ohi] <= g1:
        return "syn" if p[6] == par[0][6] else "inv"
    return "out"


def _check_net(ch, out):
    """Invariants of any net of the chained blocks `ch`: sorted, 1:1 on A and B, pieces >= ANCHOR_MIN, chains
    with >= MIN_CHAIN_SCORE matches, margins in [0, 1], each piece inside a block of its chain, on its map."""
    key = out["a_chrom"].astype(np.int64) * 10**12 + out["a_start"]
    assert (np.diff(key) > 0).all()
    for c in set(out["a_chrom"]):
        o = out[out["a_chrom"] == c]
        assert (o["a_start"][1:] >= o["a_end"][:-1]).all()
    for c in set(out["b_chrom"]):
        o = np.sort(out[out["b_chrom"] == c], order="b_start")
        assert (o["b_start"][1:] >= o["b_end"][:-1]).all()
    assert (out["a_end"] - out["a_start"] >= ANCHOR_MIN).all() and (out["b_end"] > out["b_start"]).all()
    for cid in set(out["chain"]):
        assert out[out["chain"] == cid]["matches"].sum() >= MIN_CHAIN_SCORE
    assert np.isfinite(out["margin"]).all() and ((out["margin"] >= 0) & (out["margin"] <= 1)).all()
    for r in out:  # every piece sits inside a block of its chain and follows that block's linear map
        src = ch[(ch["chain"] == r["chain"]) & (ch["a_start"] <= r["a_start"]) & (ch["a_end"] >= r["a_end"])
                 & (ch["b_start"] <= r["b_start"]) & (ch["b_end"] >= r["b_end"])]
        assert len(src) > 0 and (src["strand"] == r["strand"]).all()
        lo, hi = int(r["a_start"]), int(r["a_end"]) - 1
        ends = [(a_to_b(s, lo), a_to_b(s, hi) + 1) if r["strand"] > 0 else (a_to_b(s, hi), a_to_b(s, lo) + 1)
                for s in src]
        assert any(abs(b0 - r["b_start"]) <= 2 and abs(b1 - r["b_end"]) <= 2 for b0, b1 in ends)


def _trimmed(ch, out):
    """How many kept pieces are not whole input blocks."""
    spans = ("a_start", "a_end", "b_start", "b_end")
    whole = set(zip(*(ch[k].tolist() for k in spans)))
    return sum(k not in whole for k in zip(*(out[k].tolist() for k in spans)))


def _net_info(caplog):
    """The numbers of net_blocks' "kept" DEBUG line: kept, top, syn, inv, chain-supported, dropped, dropped A
    bp."""
    [rec] = [r for r in caplog.records if r.name == LOG and "syntenic net: kept" in r.getMessage()]
    return rec.args[1:]  # after the label prefix


# --- chaining ---------------------------------------------------------------------------------------------


def test_chain_bridges_insertion():
    # A carries a 10-kb insertion: A[0,5000)->B[0,5000) and A[15000,20000)->B[5000,10000) stay one chain
    bl = _blocks([(0, 0, 5000, 0, 0, 5000, 1, 4900), (0, 15000, 20000, 0, 5000, 10000, 1, 4900)])
    out = chain_blocks(bl)
    assert len(out) == 2 and len(set(out["chain"])) == 1
    assert out["chain_score"][0] == pytest.approx(4900 + 4900 - 0.05 * 10000)


@pytest.mark.parametrize("second, score", [
    ((1100, 2100, 1100, 2100), 2000 - math.log2(101)),  # a 100-bp gap on both genomes: log2(1 + 100)
    ((1100, 2100, 1300, 2300), 2000 - 0.05 * 200 - math.log2(101)),  # 200 bp more on B: CHAIN_INDEL_COST per bp
    ((960, 1960, 960, 1960), 2000 - 40),  # 40 bp of overlap on both: not counted twice, no gap term
])
def test_chain_dp_terms(second, score):
    """The DP's gap, indel and overlap terms, pinned on two blocks of 1000 matches."""
    a0, a1, b0, b1 = second
    out = chain_blocks(_blocks([(0, 0, 1000, 0, 0, 1000, 1, 1000), (0, a0, a1, 0, b0, b1, 1, 1000)]))
    assert set(out["chain"]) == {0} and out["chain_score"][0] == pytest.approx(score)
    assert CHAIN_INDEL_COST == 0.05


def test_chain_lookback_bounds_the_predecessors(monkeypatch):
    """Only the CHAIN_LOOKBACK blocks before a block (in A order) are tried as predecessors;
    a block between two collinear ones that cannot link them breaks the chain at a lookback of 1."""
    rows = [(0, 0, 1000, 0, 0, 1000, 1, 1000), (0, 500, 1500, 0, 900_000, 901_000, 1, 1000),
            (0, 2000, 3000, 0, 2000, 3000, 1, 1000)]
    assert len(set(chain_blocks(_blocks(rows))["chain"])) == 2  # the middle block is its own chain
    monkeypatch.setattr(chain, "CHAIN_LOOKBACK", 1)
    assert len(set(chain_blocks(_blocks(rows))["chain"])) == 3


@pytest.mark.parametrize("gap, n_chains", [(2_000_000, 2), (500_000, 1)])
def test_chain_breaks_on_far_gap(gap, n_chains):
    # second block `gap` bp further on both genomes: beyond CHAIN_MAX_GAP it starts a new chain
    bl = _blocks([(0, 0, 5000, 0, 0, 5000, 1, 4900),
                  (0, 5000 + gap, 10000 + gap, 0, 5000 + gap, 10000 + gap, 1, 4900)])
    assert len(set(chain_blocks(bl)["chain"])) == n_chains


def test_minus_strand_chain():
    # three '-' blocks with decreasing B are collinear; the same blocks on '+' are not
    rows = [(0, 0, 2000, 0, 10000, 12000, -1, 1900),
            (0, 2100, 4100, 0, 7900, 9900, -1, 1900),
            (0, 4200, 6200, 0, 5800, 7800, -1, 1900)]
    assert len(set(chain_blocks(_blocks(rows))["chain"])) == 1
    plus = [r[:6] + (1,) + r[7:] for r in rows]
    assert len(set(chain_blocks(_blocks(plus))["chain"])) == 3


@pytest.mark.parametrize("rows, spans", [
    # full-length LTR-RT in A vs solo LTR in B: each 40-kb flank block runs through a 1.5-kb LTR, so the two
    # flank alignments overlap by the LTR on B (B[40000,41500))
    ([(0, 0, 41_500, 0, 0, 41_500, 1, 40_670), (0, 46_500, 88_000, 0, 40_000, 81_500, 1, 40_670)],
     [(0, 41_500, 0, 41_500), (48_000, 88_000, 41_500, 81_500)]),
    # the same with an upstream block that both flanks could chain to
    ([(0, 0, 30_000, 0, 0, 30_000, 1, 29_400), (0, 30_000, 71_500, 0, 30_000, 71_500, 1, 40_670),
      (0, 76_500, 118_000, 0, 70_000, 111_500, 1, 40_670)],
     [(0, 30_000, 0, 30_000), (30_000, 71_500, 30_000, 71_500), (78_000, 118_000, 71_500, 111_500)]),
])
def test_chain_links_solo_ltr_flanks(rows, spans):
    ch = chain_blocks(_blocks(rows))
    assert len(set(ch["chain"])) == 1
    out = net_blocks(ch)
    assert _spans(out) == spans  # the shared LTR on B is netted once, to the first flank
    assert (out["margin"] == 1.0).all()


@pytest.mark.parametrize("paralog", [
    (0, 30, 10_020, 0, 10_000, 19_990, 1, 9000),  # tandem copy on B: A span nearly identical, B adjacent
    (0, 10_000, 19_990, 0, 30, 10_020, 1, 9000),  # tandem copy on A: B span nearly identical, A adjacent
])
def test_chain_does_not_link_tandem_paralog(paralog):
    ch = chain_blocks(_blocks([(0, 0, 10_000, 0, 0, 10_000, 1, 10_000), paralog]))
    assert len(set(ch["chain"])) == 2
    out = net_blocks(ch)
    assert _spans(out) == [(0, 10_000, 0, 10_000)]
    assert out["margin"][0] == pytest.approx(1 - 9000 / 10_000, abs=1e-6)


@pytest.mark.parametrize("n_chunks", [2, 10])
def test_chain_links_same_diagonal_overlapping_pieces(n_chunks):
    # one alignment reported in 100-kb chunks that overlap by 10 kb on both genomes (chunked mapping)
    end = 100_000 * n_chunks
    spans = [(s, min(s + 110_000, end)) for s in range(0, end, 100_000)]
    rows = [(0, s, e, 0, s, e, 1, int((e - s) * 0.98)) for s, e in spans]
    ch = chain_blocks(_blocks(rows))
    assert len(set(ch["chain"])) == 1
    out = net_blocks(ch)
    assert (out["margin"] == 1.0).all()
    for lo, hi in (("a_start", "a_end"), ("b_start", "b_end")):  # the union, covered exactly once
        assert out[lo][0] == 0 and out[hi][-1] == end and (out[lo][1:] == out[hi][:-1]).all()


def test_chain_keeps_weak_chains_and_orders_ids():
    # ids by decreasing score; a tie goes to the lower input index of the chain's first block in A order
    bl = _blocks([
        (0, 1000, 2000, 0, 1000, 2000, 1, 1000),   # input 0: second block of a two-block chain (score 2000)
        (0, 50_000, 52_000, 1, 0, 2000, 1, 2000),  # input 1: single block, score 2000
        (0, 0, 1000, 0, 0, 1000, 1, 1000),         # input 2: first block of the two-block chain
        (0, 60_000, 65_000, 2, 0, 5000, 1, 5000),  # input 3: score 5000
        (0, 70_000, 70_500, 3, 0, 500, 1, 400),    # input 4: below MIN_CHAIN_SCORE, still chained
    ])
    out = chain_blocks(bl)
    got = [(int(r["chain"]), int(r["a_start"]), float(r["chain_score"])) for r in out]
    assert got == [(0, 60_000, 5000.0), (1, 50_000, 2000.0), (2, 0, 2000.0), (2, 1000, 2000.0),
                   (3, 70_000, 400.0)]
    assert np.isnan(out["margin"]).all()


def test_chain_tail_off_used_block_scores_its_own_blocks():
    # the second block links to the first across a 550-kb B shift (f = 50000 + 3000 - 27500); extracted after
    # the first is used, it scores its own block (3000), not the marginal f[end] - f[stop] = -24500
    bl = _blocks([(0, 0, 50_000, 0, 0, 50_000, 1, 50_000), (0, 50_000, 53_000, 0, 600_000, 603_000, 1, 3000)])
    ch = chain_blocks(bl)
    assert list(ch["chain"]) == [0, 1]
    assert list(ch["chain_score"]) == [50_000, 3000]
    out = net_blocks(ch)
    assert len(out) == 2 and (out["margin"] == 1.0).all()


def test_net_tandem_paralog_branch_lowers_ortholog_margin():
    # A[12000,18000) also aligns to a tandem copy ~190 kb downstream on B (two collinear blocks P, P2). P's DP
    # predecessor is the orthologous block X1 (f[P] = 10000 + 2900 - 0.05 * 188000 - log2(2001)), so the branch is
    # extracted off a used block: marginal score f[P2] - f[X1] < 0, own score 2900 + 2900 = 5800. Scored by its
    # own blocks it challenges X2, whose A it overlaps: margin 1 - 5800/20000.
    bl = _blocks([(0, 0, 10000, 0, 0, 10000, 1, 10000), (0, 10000, 20000, 0, 10000, 20000, 1, 10000),
                  (0, 12000, 15000, 0, 200_000, 203_000, 1, 2900),
                  (0, 15000, 18000, 0, 203_000, 206_000, 1, 2900)])
    ch = chain_blocks(bl)
    assert list(ch["chain"]) == [0, 0, 1, 1]
    assert ch["chain_score"][2] == pytest.approx(5800)
    out = net_blocks(ch)
    assert _spans(out) == [(0, 10000, 0, 10000), (10000, 20000, 10000, 20000)]
    assert out["margin"][0] == 1.0
    assert out["margin"][1] == pytest.approx(1 - 5800 / 20000, abs=1e-6)


def test_chain_collapses_exact_duplicates(caplog):
    # one alignment reported three times: kept once, the copy with the most matches (the first on ties)
    blk = (0, 0, 5000, 0, 0, 5000, 1, 4900)
    bl = _blocks([blk, blk[:7] + (4800,), blk])
    bl["alen"] = [5001, 5002, 5003]  # tells the copies apart
    with caplog.at_level(logging.DEBUG, logger=LOG):
        ch = chain_blocks(bl)
    assert len(ch) == 1 and ch["matches"][0] == 4900 and ch["alen"][0] == 5001
    assert "2 exact duplicate" in caplog.text
    assert not [r for r in caplog.records if r.levelno >= logging.WARNING]
    out = net_blocks(ch)
    assert _spans(out) == [(0, 5000, 0, 5000)] and out["margin"][0] == 1.0
    # inside a chain, a collapsed duplicate does not compete with its own copy
    rows = [(0, 0, 50_000, 0, 0, 50_000, 1, 49_000), (0, 52_000, 57_000, 0, 52_000, 57_000, 1, 4900)]
    out = net_blocks(chain_blocks(_blocks(rows + rows[1:])))
    assert len(out) == 2 and len(set(out["chain"])) == 1 and (out["margin"] == 1.0).all()


@pytest.mark.parametrize("strand", [1, -1])
@pytest.mark.parametrize("shift", [0, CHAIN_OVERLAP])
def test_chain_collapses_contained_same_diagonal_copy(strand, shift, caplog):
    # chunked mapping reports A[96000,106000) whole and, from the next chunk (starting at 100000), its
    # left-truncated copy A[100000,106000) on the same diagonal (within CHAIN_OVERLAP); they cannot chain
    # (equal A ends), so the copy is dropped instead of competing with its own container
    x = (0, 96_000, 106_000, 0, 10_000, 20_000, strand, 9800)
    if strand > 0:  # diagonal b_start - a_start
        y = (0, 100_000, 106_000, 0, 14_000 + shift, 20_000, strand, 5880)
    else:  # diagonal b_end + a_start
        y = (0, 100_000, 106_000, 0, 10_000, 16_000 - shift, strand, 5880)
    with caplog.at_level(logging.DEBUG, logger=LOG):
        ch = chain_blocks(_blocks([y, x]))
    assert len(ch) == 1 and ch["a_start"][0] == 96_000
    assert "1 block(s) contained" in caplog.text
    out = net_blocks(ch)
    assert _spans(out) == [(96_000, 106_000, 10_000, 20_000)]
    assert len(set(out["chain"])) == 1 and out["margin"][0] == 1.0


@pytest.mark.parametrize("strand", [1, -1])
@pytest.mark.parametrize("shift", [CHAIN_OVERLAP + 1, 4000])
def test_chain_keeps_contained_copy_on_other_diagonal(strand, shift):
    # a tandem copy inside the orthologous block on A and on B, but off its diagonal: real competition
    x = (0, 0, 10_000, 0, 0, 10_000, strand, 10_000)
    if strand > 0:
        y = (0, 1000, 4000, 0, 1000 + shift, 4000 + shift, strand, 2900)
    else:  # the same-diagonal image of A[1000,4000) would be B[6000,9000)
        y = (0, 1000, 4000, 0, 6000 - shift, 9000 - shift, strand, 2900)
    ch = chain_blocks(_blocks([x, y]))
    assert len(ch) == 2 and len(set(ch["chain"])) == 2
    out = net_blocks(ch)
    assert _spans(out) == [(0, 10_000, 0, 10_000)]
    assert out["margin"][0] == pytest.approx(1 - 2900 / 10_000, abs=1e-6)


def test_chain_skips_empty_blocks(caplog):
    bl = _blocks([(0, 0, 5000, 0, 0, 5000, 1, 5000), (0, 6000, 6000, 0, 6000, 6100, 1, 0)])
    with caplog.at_level(logging.WARNING, logger=LOG):
        out = chain_blocks(bl)
    assert len(out) == 1 and out["a_end"][0] == 5000
    assert "1 block" in caplog.text


def test_empty_input():
    empty = np.zeros(0, dtype=BLOCK_DTYPE)
    for out in (chain_blocks(empty), net_blocks(empty)):
        assert out.dtype == BLOCK_DTYPE and len(out) == 0


# --- net --------------------------------------------------------------------------------------------------


def test_net_requires_chained_blocks():
    with pytest.raises(ValueError, match="chain_blocks"):
        net_blocks(_blocks([(0, 0, 5000, 0, 0, 5000, 1, 5000)]))


def test_net_drops_paralog():
    # chain 1: A[0,20000)->B chr0, score 20000; block A[5000,6000)->B chr1, matches 900
    bl = _blocks([(0, 0, 20000, 0, 0, 20000, 1, 20000), (0, 5000, 6000, 1, 300_000, 301_000, 1, 900)])
    out = net_blocks(chain_blocks(bl))
    assert set(out["b_chrom"]) == {0}
    m = out[(out["a_start"] <= 5000) & (out["a_end"] >= 6000)]["margin"][0]
    assert abs(m - (1 - 900 / 20000)) < 1e-3


def test_net_trims_b_overlap():
    # two A regions (distinct A chromosomes, so unchainable) claim B[5000,10000): the lower-scoring one
    # keeps only A[55000,60000) -> B[10000,15000)
    bl = _blocks([(0, 0, 10000, 0, 0, 10000, 1, 10000), (1, 50_000, 60_000, 0, 5000, 15000, 1, 9000)])
    out = net_blocks(chain_blocks(bl))
    assert _spans(out) == [(0, 10000, 0, 10000), (55_000, 60_000, 10000, 15000)]
    win, lose = out[out["a_start"] == 0][0], out[out["a_start"] == 55_000][0]
    assert win["margin"] == pytest.approx(1 - 9000 / 10000, abs=1e-6)
    assert lose["margin"] == 1.0 and lose["matches"] == 4500 and lose["alen"] == 5000


def test_net_minus_strand_chains_compete_for_b():
    # two '-' chains claim B[5000,10000); the loser keeps the A end whose B is free
    bl = _blocks([(0, 0, 10000, 0, 0, 10000, -1, 10000), (0, 100_000, 110_000, 0, 5000, 15000, -1, 8000)])
    out = net_blocks(chain_blocks(bl))
    assert _spans(out) == [(0, 10000, 0, 10000), (100_000, 105_000, 10000, 15000)]
    assert (out["strand"] == -1).all()
    assert out[out["a_start"] == 0]["margin"][0] == pytest.approx(1 - 8000 / 10000, abs=1e-6)
    lose = out[out["a_start"] == 100_000][0]
    assert lose["margin"] == 1.0 and lose["matches"] == 4000
    assert a_to_b(lose, 100_000) == 14_999 and a_to_b(lose, 104_999) == 10_000


def test_net_minus_strand_trim_on_a():
    # a '+' chain claims A[0,10000); a '-' block on A[5000,15000) keeps A[10000,15000), whose B image is B[0,5000)
    bl = _blocks([(0, 0, 10000, 1, 0, 10000, 1, 10000), (0, 5000, 15000, 0, 0, 10000, -1, 9000)])
    out = net_blocks(chain_blocks(bl))
    assert _spans(out) == [(0, 10000, 0, 10000), (10000, 15000, 0, 5000)]
    assert out[out["a_start"] == 0]["margin"][0] == pytest.approx(1 - 9000 / 10000, abs=1e-6)


@pytest.mark.parametrize("side", ["A", "B"])
@pytest.mark.parametrize("ov", [50, 150])
def test_net_edge_overlap_within_slop_is_not_competition(side, ov):
    # two abutting orthologous chains overlap by `ov` bp on one genome; only > CHAIN_OVERLAP counts as competition
    first = (0, 0, 10000, 0, 0, 10000, 1, 10000)
    if side == "A":
        second = (0, 10000 - ov, 20000, 1, 0, 10000 + ov, 1, 9000)
    else:
        second = (1, 0, 10000 + ov, 0, 10000 - ov, 20000, 1, 9000)
    out = net_blocks(chain_blocks(_blocks([first, second])))
    assert len(out) == 2 and len(set(out["chain"])) == 2
    expect = 1.0 if ov <= CHAIN_OVERLAP else 1 - 9000 / 10000
    assert out[out["matches"] == 10000]["margin"][0] == pytest.approx(expect, abs=1e-6)
    assert out[out["matches"] != 10000]["margin"][0] == 1.0


def test_net_trims_own_chain_overlap():
    # consecutive blocks of one chain overlap by 80 bp on A and B: trimmed to stay 1:1, no self-competition
    bl = _blocks([(0, 0, 5000, 0, 0, 5000, 1, 5000), (0, 4920, 10000, 0, 4920, 10000, 1, 5080)])
    ch = chain_blocks(bl)
    assert len(set(ch["chain"])) == 1
    out = net_blocks(ch)
    assert _spans(out) == [(0, 5000, 0, 5000), (5000, 10000, 5000, 10000)]
    assert list(out["margin"]) == [1.0, 1.0] and list(out["matches"]) == [5000, 5000]


def test_net_drops_blocks_below_anchor_min():
    # the 80-bp second block of a chain yields no piece, whether unchallenged or with 1 bp of its B claimed
    # by a better chain
    chain = [(0, 0, 5000, 0, 0, 5000, 1, 5000), (0, 5100, 5180, 0, 5100, 5180, 1, 80)]
    out = net_blocks(chain_blocks(_blocks(chain)))
    assert _spans(out) == [(0, 5000, 0, 5000)] and out["margin"][0] == 1.0
    out = net_blocks(chain_blocks(_blocks(chain + [(1, 0, 20_000, 0, 5179, 25_179, 1, 20_000)])))
    assert _spans(out) == [(0, 5000, 0, 5000), (0, 20_000, 5179, 25_179)]


def test_net_drops_chain_below_min_score():
    weak = chain_blocks(_blocks([(0, 0, 900, 0, 0, 900, 1, 800)]))
    assert len(weak) == 1 and len(net_blocks(weak)) == 0
    # a 1500-match chain trimmed to A[10000,10500) keeps 500 matches: dropped, yet it challenged the winner
    bl = _blocks([(0, 0, 10000, 0, 0, 10000, 1, 10000), (0, 9000, 10500, 1, 0, 1500, 1, 1500)])
    out = net_blocks(chain_blocks(bl))
    assert _spans(out) == [(0, 10000, 0, 10000)]
    assert out["margin"][0] == pytest.approx(1 - 1500 / 10000, abs=1e-6)


@pytest.mark.parametrize("seed, frac_ortho", [(1, 0.5), (2, 0.3)])
def test_net_is_one_to_one_on_random_blocks(seed, frac_ortho):
    bl = _random_blocks(3000, seed=seed, frac_ortho=frac_ortho)
    ch = chain_blocks(bl)
    assert len(ch) == len(bl)
    by_id = np.zeros(ch["chain"].max() + 1)
    by_id[ch["chain"]] = ch["chain_score"]
    assert set(ch["chain"]) == set(range(len(by_id))) and (np.diff(by_id) <= 0).all()

    out = net_blocks(ch)
    # the random paralogs land anywhere: nearly all are fills out of place in both genomes and are dropped,
    # so few chains stay besides the collinear runs. As rivals they still lower the runs' margins, and the few
    # kept are trimmed; test_net_without_drops_equals_the_flat_net has competition between many kept chains.
    assert len(set(out["chain"])) >= 3 and _trimmed(ch, out) >= 3 and (out["margin"] < 1).any()
    _check_net(ch, out)


def test_net_matches_the_rule_written_plainly():
    # the indexed parent lookups agree with a scan of every kept chain, on fixtures covering every place
    places = Counter()
    fixtures = [_genome_blocks(seed, n_chrom=3 + seed % 3) for seed in range(4)]
    fixtures += [_random_blocks(3000, seed=seed, n_chrom=2 + seed, frac_ortho=0.3 + 0.1 * seed)
                 for seed in range(3)]
    for bl in fixtures:
        ch = chain_blocks(bl)
        out = net_blocks(ch)
        assert out.tobytes() == _reference_net(ch, places=places).tobytes()
        _check_net(ch, out)
    kept = [k for k in places if "out" not in k]
    assert {"top", "syn", "inv", "chain"} <= {p for k in kept for p in k}  # every kind of support
    assert sum(n for k, n in places.items() if "out" in k) > 100  # and drops


def test_net_without_drops_equals_the_flat_net(caplog):
    # chains that compete (credit, trimming) but never nest: nothing is out of place, so the syntenic net is
    # the flat 1:1 net, margins included
    ch = chain_blocks(_competing_blocks(400, seed=5))
    assert len(set(ch["chain"])) == 400  # nothing chains
    with caplog.at_level(logging.DEBUG, logger=LOG):
        out = net_blocks(ch)
    assert _net_info(caplog)[5:] == (0, 0)  # no piece dropped
    assert len(set(out["chain"])) > 10 and _trimmed(ch, out) > 20  # it exercises competition and trimming
    assert (out["margin"] < 1).sum() > 20 and (out["margin"] == 1).any()
    assert out.tobytes() == _reference_net(ch, syntenic=False).tobytes()
    _check_net(ch, out)


# Parent chain P on A0 -> B0: L1 A[0,50k) -> B[0,50k); M A[60k,110k) -> B[50k,100k) after an A-only insertion
# A[50k,60k); R A[110k,160k) -> B[110k,160k) after a B-only insertion B[100k,110k). Its A gap A[50k,60k)
# faces an empty B gap at B 50k; its B gap B[100k,110k) faces an empty A gap at A 110k.
_P = [(0, 0, 50 * K, 0, 0, 50 * K, 1, 50 * K), (0, 60 * K, 110 * K, 0, 50 * K, 100 * K, 1, 50 * K),
      (0, 110 * K, 160 * K, 0, 110 * K, 160 * K, 1, 50 * K)]
_P_SPANS = [(0, 50 * K, 0, 50 * K), (60 * K, 110 * K, 50 * K, 100 * K), (110 * K, 160 * K, 110 * K, 160 * K)]


def _net(rows, caplog):
    """net_blocks(chain_blocks(rows)) and its INFO numbers."""
    caplog.clear()
    with caplog.at_level(logging.DEBUG, logger=LOG):
        out = net_blocks(chain_blocks(_blocks(rows)))
    return out, _net_info(caplog)


@pytest.mark.parametrize("strand", [1, -1])
def test_net_drops_paralog_fill_out_of_place_in_both_genomes(strand, caplog):
    # an A-only TE copy inside P's A gap paired with a B-only copy inside P's B gap (8 kb at 85%): a chain of
    # its own, in P's gaps in both genomes and out of place in both
    out, info = _net(_P + [(0, 51 * K, 59 * K, 0, 101 * K, 109 * K, strand, 6800)], caplog)
    assert _spans(out) == _P_SPANS and (out["margin"] == 1.0).all()  # it overlapped nothing
    assert info == (3, 3, 0, 0, 0, 1, 8 * K)


@pytest.mark.parametrize("fill, kind", [
    # same strand, A and B inside the parent's gaps, off its diagonal so that the chain does not take it
    ((0, 52 * K, 55 * K, 0, 93 * K, 96 * K, 1, 3000), "syn"),
    ((0, 60 * K, 65 * K, 0, 60 * K, 65 * K, -1, 5000), "inv"),  # a small inversion
])
def test_net_keeps_fill_in_its_parents_gap(fill, kind, caplog):
    # parent A[0,50k) -> B[0,50k), A[100k,150k) -> B[100k,150k): one chain with a 50-kb gap on both genomes
    parent = [(0, 0, 50 * K, 0, 0, 50 * K, 1, 50 * K), (0, 100 * K, 150 * K, 0, 100 * K, 150 * K, 1, 50 * K)]
    assert len(set(chain_blocks(_blocks(parent + [fill]))["chain"])) == 2
    out, info = _net(parent + [fill], caplog)
    assert len(out) == 3 and _spans(out)[1] == fill[1:3] + fill[4:6]
    assert info == (3, 2, int(kind == "syn"), int(kind == "inv"), 0, 0, 0)


@pytest.mark.parametrize("strand, kind", [(-1, "syn"), (1, "inv")])
def test_net_keeps_fill_in_a_minus_strand_parents_gap(strand, kind, caplog):
    """A minus-strand parent, L = A[0,50k) -> B[100k,150k) and R = A[100k,150k) -> B[0,50k);
    the fill A[52k,55k) -> B[53k,56k) lies in its matching B gap: syn on the parent's strand, inv on the other."""
    parent = [(0, 0, 50 * K, 0, 100 * K, 150 * K, -1, 50 * K), (0, 100 * K, 150 * K, 0, 0, 50 * K, -1, 50 * K)]
    fill = (0, 52 * K, 55 * K, 0, 53 * K, 56 * K, strand, 3000)
    assert len(set(chain_blocks(_blocks(parent + [fill]))["chain"])) == 2
    out, info = _net(parent + [fill], caplog)
    assert len(out) == 3 and (52 * K, 55 * K, 53 * K, 56 * K) in _spans(out)
    assert info == (3, 2, int(kind == "syn"), int(kind == "inv"), 0, 0, 0)


def test_net_keeps_chain_top_level_in_both_genomes(caplog):
    # past P's A extent and past its B extent: a top-level chain (or a top-level paralog; the rule cannot
    # tell them apart)
    out, info = _net(_P + [(0, 200 * K, 210 * K, 0, 300 * K, 310 * K, -1, 9000)], caplog)
    assert _spans(out) == _P_SPANS + [(200 * K, 210 * K, 300 * K, 310 * K)]
    assert info == (4, 4, 0, 0, 0, 0, 0)


def test_net_keeps_crept_in_inversion_by_chain_support(caplog):
    # A 2-Mb inversion (A[100k,2.1M) -> B = 2.2M - A, 100-kb blocks at 90%) between forward chains. The better
    # right chain (A[2.1M,4.1M) -> B[2.1M,4.1M)) starts with a 5-kb forward piece A[1.9M,1.905M) ->
    # B[1.95M,1.955M) inside the inversion, so its A extent covers the inversion's A[1.905M,2.1M) and its B
    # extent the B image of A[100k,245k). Those four inverted pieces are out of place in one genome each, but
    # their chain is top level in both: kept. (Top level in B alone is no escape: see the leak test.)
    rows = [(0, 0, 100 * K, 0, 0, 100 * K, 1, 100 * K),
            (0, 1900 * K, 1905 * K, 0, 1950 * K, 1955 * K, 1, 5 * K)]
    for k in range(0, 2000 * K, 100 * K):
        rows.append((0, 100 * K + k, 200 * K + k, 0, 2000 * K - k, 2100 * K - k, -1, 90 * K))
        rows.append((0, 2100 * K + k, 2200 * K + k, 0, 2100 * K + k, 2200 * K + k, 1, 100 * K))
    ch = chain_blocks(_blocks(rows))
    assert [len(ch[ch["chain"] == c]) for c in range(3)] == [21, 20, 1]  # right chain, inversion, left chain
    out, info = _net(rows, caplog)
    inv = out[out["strand"] < 0]
    assert len(inv) == 21 and int((inv["a_end"] - inv["a_start"]).sum()) == 1990 * K  # 5 kb lost on A, 5 on B
    crept = {(1905 * K, 2000 * K, 200 * K, 295 * K), (2000 * K, 2100 * K, 100 * K, 200 * K),
             (100 * K, 200 * K, 2000 * K, 2100 * K), (200 * K, 245 * K, 1955 * K, 2000 * K)}
    assert crept <= set(_spans(inv))
    assert info[4:] == (4, 0, 0)  # those four are chain-supported; nothing is dropped


@pytest.mark.parametrize("b_seq, b0", [
    (0, 2000 * K),  # B past P's B extent: top level in B only (a leak)
    (1, 0),  # an unplaced B scaffold: the named cost for an ortholog there whose A side is in a chain's gap
])
def test_net_drops_leak_into_uncovered_b(b_seq, b0, caplog):
    # A inside P's A gap, B in B sequence that no chain covers: top level in B, but its chain has no top-level
    # piece in A, so nothing places it: dropped
    out, info = _net(_P + [(0, 51 * K, 59 * K, b_seq, b0, b0 + 8 * K, 1, 7000)], caplog)
    assert _spans(out) == _P_SPANS
    assert info == (3, 3, 0, 0, 0, 1, 8 * K)


def test_net_drops_chain_top_level_in_a_only(caplog):
    # the mirror: past P's A extent (top level in A) but inside P's B gap B[100k,110k), whose A gap is empty:
    # out of place in B, and its chain has no top-level piece in B: dropped (a top-level paralog)
    out, info = _net(_P + [(0, 200 * K, 208 * K, 0, 101 * K, 109 * K, 1, 7000)], caplog)
    assert _spans(out) == _P_SPANS
    assert info == (3, 3, 0, 0, 0, 1, 8 * K)


def test_net_drops_short_transposition(caplog):
    # the rule's named cost: A[50k,60k) moved 50 kb to B[100k,110k) (identical; 10 kb < CHAIN_MAX_GAP) sits
    # in P's gaps in both genomes like a paralog: dropped. An element in it is windowed by P's flanks, missing
    # in B.
    assert 10 * K < CHAIN_MAX_GAP
    out, info = _net(_P + [(0, 50 * K, 60 * K, 0, 100 * K, 110 * K, 1, 10 * K)], caplog)
    assert _spans(out) == _P_SPANS
    assert info[5:] == (1, 10 * K)


def test_net_dropped_piece_still_lowers_margin(caplog):
    # A paralog A[48k,58k) -> B[101k,111k) overlaps P's L1 by 2 kb on A and R by 1 kb on B (> CHAIN_OVERLAP),
    # and its free part A[50k,57k) -> B[103k,110k) is out of place in both genomes: dropped. It is still a
    # rival alignment, so L1's and R's margins fall to 1 - 9000 / 149000; M, unchallenged, keeps 1.
    rows = _P + [(0, 48 * K, 58 * K, 0, 101 * K, 111 * K, 1, 9000)]
    ch = chain_blocks(_blocks(rows))
    assert {int(c): float(s) for c, s in zip(ch["chain"], ch["chain_score"])} == {0: 149_000, 1: 9000}
    out, info = _net(rows, caplog)
    assert _spans(out) == _P_SPANS and info[5:] == (1, 7 * K)
    assert list(out["margin"]) == pytest.approx([1 - 9000 / 149_000, 1.0, 1 - 9000 / 149_000])


def test_net_rejects_mixed_chain():
    bl = _blocks([(0, 0, 5000, 0, 0, 5000, 1, 5000), (0, 6000, 11000, 1, 6000, 11000, 1, 5000)])
    bl["chain"], bl["chain_score"] = 0, 10_000.0  # one chain over two B sequences: not chain_blocks output
    with pytest.raises(ValueError, match="one A sequence, one B sequence and one strand"):
        net_blocks(bl)


def test_net_logs_one_debug_line(caplog):
    with caplog.at_level(logging.DEBUG, logger=LOG):
        net_blocks(chain_blocks(_blocks(_P + [(0, 51 * K, 59 * K, 0, 101 * K, 109 * K, 1, 6800)])))
    info = [r.getMessage() for r in caplog.records if r.getMessage().startswith("syntenic net: kept")]
    assert info == [("syntenic net: kept 3 pieces (3 top in both genomes, 0 syn, 0 inv, 0 chain-supported); "
                     "dropped 1 out of place, 8000 bp of A in total")]
    assert "out/out 1" in caplog.text and "1 chains lost pieces" in caplog.text
    caplog.clear()
    with caplog.at_level(logging.INFO, logger=LOG):  # not at INFO: the default log is minimal
        net_blocks(chain_blocks(_blocks(_P)))
    assert not caplog.records
    with caplog.at_level(logging.DEBUG, logger=LOG):  # a route nets each pair under its label
        net_blocks(chain_blocks(_blocks(_P)), label="A vs B (flank-mm2)")
    assert any(r.getMessage().startswith("A vs B (flank-mm2): syntenic net: kept 3 pieces") for r in caplog.records)


def test_net_50k_blocks_in_seconds():
    # about genome scale per sequence pair: 60k blocks over 5000 sequence pairs, ~8k chains kept
    chain_blocks(_random_blocks(100, seed=0))  # JIT warm-up (compile or load the cache)
    ch = chain_blocks(_random_blocks(60_000, seed=7, n_chrom=5000))
    t0 = time.perf_counter()
    out = net_blocks(ch)
    elapsed = time.perf_counter() - t0
    assert len(set(out["chain"])) > 5000
    assert elapsed < 5, f"net of 60k blocks took {elapsed:.1f} s"


@pytest.mark.slow
def test_chain_and_net_1e6_blocks_under_60s():
    chain_blocks(_random_blocks(100, seed=0))  # JIT warm-up (compile or load the cache)
    bl = _random_blocks(1_000_000, seed=7, n_chrom=10)
    t0 = time.perf_counter()
    out = net_blocks(chain_blocks(bl))
    elapsed = time.perf_counter() - t0
    assert len(out) > 0
    assert elapsed < 60, f"chain + net on 1e6 blocks took {elapsed:.1f} s"
