"""Inverted blocks (analyses.geneflow.inversions): the run-merge rule on one sequence group, blocks of a
simulated run with an inversion, the inside mask and inversions.tsv."""
from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest
from analyses.geneflow.catalog_toy import loc, simulate_run, write_run

from ltrtrace.analyses.geneflow.inversions import (
    INV_COLUMNS,
    InvBlock,
    _find,
    inside,
    inversion_rows,
    inverted_blocks,
    syntenic_pairs,
    write_inversions,
)


def reversed_blocks(y):
    """_find's blocks as (row indexes, n_discordant)."""
    return [(rows, nd) for rows, nd, _, _ in _find(y)]


X = 1001 + np.arange(1500) * 998_000 // 1500 + 250        # 1,500 evenly spaced row midpoints on 1 Mb

def _inv(y, a, b):
    y = y.copy(); m = (X >= a) & (X <= b); y[m] = a + b - X[m]; return y

def make_case(case: str) -> np.ndarray:
    y0 = _inv(X.astype(float), 300_001, 700_000)
    mid = np.flatnonzero((X >= 300_001) & (X <= 700_000))
    if case == "clean inversion":
        return y0
    if case.startswith("one locus moved"):
        y = y0.copy(); y[mid[len(mid) // 2]] += 50_000 if "+" in case else -50_000; return y
    if case == "two adjacent loci moved":
        y = y0.copy(); y[mid[300]] += 50_000; y[mid[301]] += 60_000; return y
    if case == "two separate inversions on one seqid":
        return _inv(_inv(X.astype(float), 100_001, 200_000), 600_001, 800_000)
    if case == "30-locus nested re-inversion":
        y = y0.copy(); s = (X >= 480_001) & (X <= 500_000); y[s] = 500_001 + (X[s] - 480_001); return y
    if case == "collinear":
        return X.astype(float)
    if case == "wholly reversed seqid":
        return (2_000_000 - X).astype(float)
    y = X.astype(float); y[700:704] = y[700:704][::-1]; return y          # "4 reversed loci"

@pytest.mark.parametrize("case, expected", [         # expected = n_discordant of each block, in x order
    ("clean inversion", [0]), ("one locus moved +50 kb", [1]), ("one locus moved -50 kb", [1]),
    ("two adjacent loci moved", [2]), ("two separate inversions on one seqid", [0, 0]),
    ("30-locus nested re-inversion", [29]), ("collinear", []), ("wholly reversed seqid", []), ("4 reversed loci", [])])
def test_reversed_blocks_cases(case, expected):
    blocks = reversed_blocks(make_case(case))
    assert [d for _, d in blocks] == expected
    if case == "clean inversion":
        assert len(blocks[0][0]) == 601
    if case == "two separate inversions on one seqid":
        assert [len(b) for b, _ in blocks] == [150, 301]

TOY_LENGTH = 150  # the toy's element length: a locus starting at or before a breakpoint ends up to this far past it


def test_block_found(tmp_path, toy_fastas, seed=2):  # seeds 1-7 all pass (1.5 s each)
    cat = simulate_run(tmp_path, seed=seed, inversion=("S", "chr2", 300_001, 700_000), fasta_dir=toy_fastas).load()
    blocks = inverted_blocks(cat)
    pair = lambda a, b: [blk for blk in blocks if {blk.x, blk.y} == {a, b}]
    for other in ("R", "D", "O"):                    # every genome against S carries the inversion, S~O included
        (blk,) = pair("S", other)
        assert blk.n_strand == blk.n_loci and blk.n_discordant == 0
        for name in (other, "S"):                    # both axes of the pair
            on_x = blk.x == name
            seq = blk.seq_x if on_x else blk.seq_y
            lo, hi = (blk.start_x, blk.end_x) if on_x else (blk.start_y, blk.end_y)
            left, right = (blk.bp_left_x, blk.bp_right_x) if on_x else (blk.bp_left_y, blk.bp_right_y)
            assert cat.seqids[name][seq] == "chr2" and 300_001 <= lo and hi <= 700_000 + TOY_LENGTH
            # each true breakpoint lies inside its interval, and the interval is no wider than the gap
            # between the flanking loci resolved in both genomes of the pair
            g, go = cat.g(name), cat.g(other if name == "S" else "S")
            ok = (cat.cls(name) >= 0) & (cat.cls(other if name == "S" else "S") >= 0) & (cat.seq[:, g] == seq) \
                & (cat.seq[:, go] >= 0)
            order = np.argsort(cat.start[ok, g])
            st, en = cat.start[ok, g][order], cat.end[ok, g][order]
            for (a, b), bp in ((left, 300_001), (right, 700_000)):
                assert (a if bp == 300_001 else a - TOY_LENGTH) <= bp <= b
                before = en[st < bp][-1]
                after = st[st > bp][0] if bp == 300_001 else st[st >= bp][0]
                assert b - a <= after - before
    assert not pair("R", "D") and not pair("R", "O") and not pair("D", "O")


def test_displaced_locus_keeps_y_extent(tmp_path, toy_fastas):
    cat = simulate_run(tmp_path, seed=2, inversion=("S", "chr2", 300_001, 700_000), fasta_dir=toy_fastas).load()
    (clean,) = [b for b in inverted_blocks(cat) if (b.x, b.y) == ("S", "D")]
    s, d = cat.g("S"), cat.g("D")
    chr2 = cat.seqids["S"].index("chr2")
    both = np.flatnonzero((cat.cls("S") >= 0) & (cat.cls("D") >= 0) & (cat.seq[:, s] == chr2) & (cat.seq[:, d] >= 0)
                          & (cat.start[:, s] > 450_000) & (cat.start[:, s] < 550_000))
    i = both[len(both) // 2]
    cat.start[i, d] += 5_000_000   # one locus transposed 5 Mb away on D
    cat.end[i, d] += 5_000_000
    (blk,) = [b for b in inverted_blocks(cat) if (b.x, b.y) == ("S", "D")]
    assert blk.n_discordant == 1 and (blk.start_y, blk.end_y) == (clean.start_y, clean.end_y)
    assert blk.end_y < 700_000 + TOY_LENGTH
    got = inside([blk], "D", {"S"}, np.array([chr2, chr2]), np.array([500_000, 5_000_500]))
    assert got.tolist() == [True, False]


def test_dominant_flip_finding(tmp_path, toy_fastas):
    cat = simulate_run(tmp_path, seed=2, inversion=("S", "chr2", 300_001, 700_000), fasta_dir=toy_fastas).load()
    assert not any("dominant y direction" in t for _, t in inversion_rows(inverted_blocks(cat), cat)[1].items())
    big = simulate_run(tmp_path / "big", seed=2, inversion=("S", "chr2", 100_001, 900_000), fasta_dir=toy_fastas).load()
    blocks = inverted_blocks(big)
    assert any(b.dominant < 0 for b in blocks)   # an inversion over 80% of the group: the remainder is reported
    assert any("dominant y direction" in t for _, t in inversion_rows(blocks, big)[1].items())


def test_inside(tmp_path, toy_fastas):
    cat = simulate_run(tmp_path, seed=2, inversion=("S", "chr2", 300_001, 700_000), fasta_dir=toy_fastas).load()
    blocks = inverted_blocks(cat)
    chr2 = cat.seqids["R"].index("chr2")
    got = inside(blocks, "R", {"S"}, np.array([chr2, chr2, chr2, 0]), np.array([500_000, 200_000, 800_000, 500_000]))
    assert got.tolist() == [True, False, False, False]
    assert inside(blocks, "R", {"D"}, np.array([chr2]), np.array([500_000])).tolist() == [False]

def test_write_inversions(tmp_path):
    blk = InvBlock("R", "S", 1, 1, 300_001, 699_850, 300_150, 700_000, (299_400, 300_001), (699_850, 700_400),
                   None, (700_000, 700_600), 601, 601, 0)
    cat = SimpleNamespace(seqids={"R": ("chr1", "chr2"), "S": ("chr1", "chr2")})
    write_inversions([blk], cat, tmp_path / "inv.tsv")
    head, row = (tmp_path / "inv.tsv").read_text().splitlines()
    assert head.split("\t") == list(INV_COLUMNS)
    assert row.split("\t") == ["R", "chr2", "300001", "699850", "299400-300001", "699850-700400",
                               "S", "chr2", "300150", "700000", ".", "700000-700600", "601", "601", "0"]


def test_inversion_rows(tmp_path, toy_fastas):
    cat = simulate_run(tmp_path, seed=2, inversion=("S", "chr2", 300_001, 700_000), fasta_dir=toy_fastas).load()
    rows, f = inversion_rows(inverted_blocks(cat), cat)
    got = {(r.scope, r.statistic): r.value for r in rows}
    assert got[("R~S", "n_blocks")] == 1 and got[("R~D", "n_blocks")] == 0
    assert got[("R~S", "n_blocks_1Mb")] == 0 and got[("R~S", "bp_inverted_x")] > 300_000
    assert any("INV_MIN_LOCI" in t for _, t in f.items())


# ---------------------------------------------------------------- syntenic partners only

PAIR = [("R", "tR", None, True), ("S", "tS", None, False)]     # R discovers; S genotyped from R's calls
SLOT = 4_000                                                    # collinear loci every 4 kb, 2 kb apart from the rest


def _backbone(seqids=("chr1", "chr2", "chr3"), n=200):
    """Collinear loci, R and S at the same place: 10,001 + 4 kb x j on each sequence."""
    return [loc(10_001 + SLOT * j, "R:P* S:P", seqid=s) for s in seqids for j in range(n)]


def _cat(tmp_path, toy_fastas, loci, **kw):
    return write_run(tmp_path, sorted(loci, key=lambda x: (x.seqid, x.pos)), genomes=PAIR,
                     seqids=kw.pop("seqids", ("chr1", "chr2", "chr3")), fasta_dir=toy_fastas, **kw).load()


def _blocks(cat):
    return [(cat.seqids["R"][b.seq_x], cat.seqids["S"][b.seq_y], b.start_x, b.end_x, b.n_loci, b.n_discordant)
            for b in inverted_blocks(cat)]


def test_translocated_handful_is_no_block(tmp_path, toy_fastas):
    """A translocation artifact: 11 loci of R's chr1, scattered over 760 kb, sit together on S's chr3, 6 in
    order and 5 reversed. Grouped by sequence pair alone they make a 5-locus block; chr1 and chr3 are not
    syntenic partners (each holds 200 loci with its own namesake), so they are not searched. The real
    inversion on chr2, a homologous pair, is the one block."""
    ys = [12_001 + SLOT * k for k in (50, 51, 52, 53, 54, 55, 65, 64, 63, 62, 61)]     # S chr3, in R's x order
    assert [len(r) for r, _ in reversed_blocks(np.array(ys, dtype=float))] == [5]   # the trap is in the toy
    handful = [loc(12_001 + 76_000 * k, "R:P* S:P", seqid="chr1", at={"S": ("chr3", y)}) for k, y in enumerate(ys)]
    cat = _cat(tmp_path, toy_fastas, _backbone() + handful, inversion=("S", "chr2", 300_001, 700_000))
    assert syntenic_pairs(cat, "R", "S") == {(0, 0): 200, (1, 1): 200, (2, 2): 200}
    assert _blocks(cat) == [("chr2", "chr2", 302_001, 698_150, 100, 0)]   # 302,001-698,001 (+150 bp element)
    rows, f = inversion_rows(inverted_blocks(cat), cat)
    assert next(r for r in rows if r.statistic == "n_blocks").n == 600            # the searched loci
    line = next(t for _, t in f.items() if t.startswith("R~S: "))
    assert "11 of the 611 loci resolved in both lie on 1 sequence pair(s) that are not syntenic partners" in line


def test_scattered_misplacements_in_a_homologous_pair(tmp_path, toy_fastas):
    """Six loci of R's chr1, spread along it, sit together and reversed in 30 kb of S's chr1, among chr1's
    collinear loci: each is one discordant row between collinear neighbors, so no block forms (reciprocal
    dominance leaves no sparse artifact here)."""
    moved = [loc(12_001 + 120_000 * k, "R:P* S:P", seqid="chr1", at={"S": ("chr1", 924_001 - SLOT * k)})
             for k in range(6)]
    cat = _cat(tmp_path, toy_fastas, _backbone() + moved)
    assert syntenic_pairs(cat, "R", "S") == {(0, 0): 206, (1, 1): 200, (2, 2): 200}
    assert _blocks(cat) == []


def _reciprocal():
    """R chr1 = S chr1 slots 0-69 + S chr2 slots 40-99; R chr2 = S chr2 slots 0-39 + S chr1 slots 70-99 (100
    slots per S sequence, 6 kb apart). Inverted stretches: R chr1 rows 20-39 (in its S chr1 piece) and 90-109
    (in its S chr2 piece); R chr2 rows 50-59 (in its S chr1 piece, the smaller share of both sequences)."""
    pieces = {"chr1": [("chr1", j) for j in range(70)] + [("chr2", j) for j in range(40, 100)],
              "chr2": [("chr2", j) for j in range(40)] + [("chr1", j) for j in range(70, 100)]}
    for r, (a, b) in (("chr1", (20, 39)), ("chr1", (90, 109)), ("chr2", (50, 59))):
        pieces[r][a:b + 1] = pieces[r][a:b + 1][::-1]
    return [loc(10_001 + 6_000 * k, "R:P* S:P", seqid=r, at={"S": (s, 10_001 + 6_000 * j)})
            for r, slots in pieces.items() for k, (s, j) in enumerate(slots)]


def test_reciprocal_translocation_keeps_both_partners(tmp_path, toy_fastas):
    """In a reciprocal translocation R chr1 holds the larger share of both S chr1 (70 of 100 loci) and S chr2
    (60 of 100): both are its syntenic partners, and an inversion inside each piece is found. R chr2's share of
    S chr1 (30) is the smaller of both sequences' shares: that pair is not searched, so its inversion is missed
    (the rule's stated cost)."""
    cat = _cat(tmp_path, toy_fastas, _reciprocal(), seqids=("chr1", "chr2"))
    assert syntenic_pairs(cat, "R", "S") == {(0, 0): 70, (0, 1): 60, (1, 1): 40}
    assert _blocks(cat) == [("chr1", "chr1", 10_001 + 6_000 * 20, 10_150 + 6_000 * 39, 20, 0),
                            ("chr1", "chr2", 10_001 + 6_000 * 90, 10_150 + 6_000 * 109, 20, 0)]


def _scaffold_toy():
    """S also has a 100-kb scaffold scf1 whose 10 loci lie in R spread along chr1 (300-660 kb), among chr1's
    collinear loci, 5 in order and 5 reversed: scf1's dominant partner is R chr1, so the pair is a syntenic one,
    and the 5 reversed loci make a run spanning 160 kb of R chr1 among 40 collinear loci."""
    ys = [2_001 + SLOT * k for k in (1, 2, 3, 4, 5, 15, 14, 13, 12, 11)]          # scf1, in R's x order
    moved = [loc(300_001 + 40_000 * k, "R:P* S:P", seqid="chr1", at={"S": ("scf1", y)}) for k, y in enumerate(ys)]
    return ys, moved


def test_scaffold_spread_along_a_chromosome_is_no_block(tmp_path, toy_fastas):
    """The artifact reciprocal dominance leaves: a scaffold whose few loci lie spread along a chromosome
    keeps its pair, and its reversed run would be a 160-kb block on 5 loci; the 40 collinear chr1 loci inside
    that stretch of R map to S chr1, so the block is not contiguous on R and is dropped (INV_FOREIGN_MAX)."""
    ys, moved = _scaffold_toy()
    assert [len(r) for r, _ in reversed_blocks(np.array(ys, dtype=float))] == [5]   # the trap is in the toy
    cat = _cat(tmp_path, toy_fastas, _backbone() + moved, seqids=("chr1", "chr2", "chr3", "scf1"),
               seq_lens={"scf1": 100_000})
    assert syntenic_pairs(cat, "R", "S") == {(0, 0): 200, (0, 3): 10, (1, 1): 200, (2, 2): 200}
    assert _blocks(cat) == []


def test_alt_haplotig_does_not_break_a_block(tmp_path, toy_fastas):
    """S's assembly holds R chr1's 300-700 kb region twice: on chr1, inverted, and on the alternate
    haplotig alt1, collinear (of the region's loci, one in three placed on chr1, two on alt1). alt1 is a
    syntenic partner of R chr1 (all its loci are R chr1's) and chr1~chr1 is mutually dominant, so the 64 alt1
    loci inside the inversion's extent do not count against its 33 loci: it is kept. The scaffold toy above, a
    pair kept from one side only, is dropped."""
    loci = []
    for j in range(200):
        pos = 10_001 + SLOT * j
        on_alt = 300_001 <= pos <= 700_000 and j % 3 != 0
        loci.append(loc(pos, "R:P* S:P", seqid="chr1", at={"S": ("alt1", pos - 290_000)} if on_alt else {}))
    loci += [loc(10_001 + SLOT * j, "R:P* S:P", seqid="chr2") for j in range(200)]
    cat = _cat(tmp_path, toy_fastas, loci, seqids=("chr1", "chr2", "alt1"), seq_lens={"alt1": 500_000},
               inversion=("S", "chr1", 300_001, 700_000))
    assert syntenic_pairs(cat, "R", "S") == {(0, 0): 133, (0, 2): 67, (1, 1): 200}
    assert _blocks(cat) == [("chr1", "chr1", 310_001, 694_150, 33, 0)]
