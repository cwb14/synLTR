"""Alignment blocks (routes.assembly.blocks): PAF records to blocks, split at long indels; the linear A-to-B map."""
import gzip
import logging

import numpy as np
import pytest

from ltrtrace.routes.assembly import blocks as blk
from ltrtrace.routes.assembly.blocks import BLOCK_DTYPE, a_to_b, load_blocks, read_paf, save_blocks

LOG = "ltrtrace.routes.assembly.blocks"
LINE = "a1\t1000\t100\t200\t+\tb1\t5000\t1000\t1100\t95\t100\t60\n"


def test_read_paf_and_map(tmp_path):
    p = tmp_path / "x.paf"
    p.write_text(LINE)
    bl = read_paf(p, {"a1": 0}, {"b1": 0})
    assert bl.dtype == BLOCK_DTYPE
    assert a_to_b(bl[0], 150) == 1050
    assert bl[0]["chain"] == -1 and bl[0]["chain_score"] == 0.0 and np.isnan(bl[0]["margin"])
    assert (bl[0]["matches"], bl[0]["alen"], bl[0]["strand"]) == (95, 100, 1)


def test_minus_strand_mapping(tmp_path):
    p = tmp_path / "x.paf"
    p.write_text(LINE.replace("\t+\t", "\t-\t"))
    bl = read_paf(p, {"a1": 0}, {"b1": 0})
    assert bl[0]["strand"] == -1
    assert a_to_b(bl[0], 100) == 1099
    assert a_to_b(bl[0], 199) == 1000
    assert a_to_b(bl[0], 5000) == 1000  # clamped


def test_swap_and_unknown_names(tmp_path, caplog):
    p = tmp_path / "x.paf.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("# comment\nshort\tline\n")
        fh.write(LINE)
        fh.write(LINE.replace("a1", "zz"))
        fh.write(LINE.replace("b1", "yy"))
    with caplog.at_level(logging.WARNING, logger=LOG):
        bl = read_paf(p, {"b1": 0}, {"a1": 0}, swap=True)
    assert len(bl) == 1
    assert (bl[0]["a_start"], bl[0]["a_end"]) == (1000, 1100)
    assert (bl[0]["b_start"], bl[0]["b_end"]) == (100, 200)
    warns = [r.getMessage() for r in caplog.records if r.name == LOG]
    assert warns[0] == f"{p}:2: PAF record skipped, malformed: 2 fields, PAF has 12"  # "short\tline"
    assert "skipped 2 PAF records with seqids not in the index" in warns[-1] and len(warns) == 3
    assert bl[0]["matches"] == 95 and bl[0]["strand"] == 1


def test_malformed_paf_records_are_skipped(tmp_path, caplog):
    """A PAF record with a field that is not an integer or an empty span is
    skipped with a warning naming the file and line, not an abort; a_to_b of an empty A span does not divide
    by zero."""
    p = tmp_path / "m.paf"
    p.write_text(LINE.replace("\t100\t200\t", "\t1OO\t200\t") + LINE + LINE.replace("\t1100\t", "\t1000\t"))
    with caplog.at_level(logging.WARNING, logger=LOG):
        bl = read_paf(p, {"a1": 0}, {"b1": 0})
    assert len(bl) == 1 and (bl[0]["a_start"], bl[0]["b_end"]) == (100, 1100)
    warns = [r.getMessage() for r in caplog.records]
    assert warns[0].startswith(f"{p}:1: PAF record skipped, malformed: invalid literal")
    assert warns[1] == f"{p}:3: PAF record skipped, malformed: empty or negative span: query 100-200, target 1000-1000"
    assert warns[2] == f"skipped 2 malformed PAF record(s) ({p})"
    empty = bl.copy()
    empty["a_end"] = empty["a_start"]
    assert a_to_b(empty[0], 100) == 1000


def test_save_load_roundtrip(tmp_path):
    p = tmp_path / "x.paf"
    p.write_text(LINE)
    bl = read_paf(p, {"a1": 0}, {"b1": 0})
    out = tmp_path / "b.npz"
    save_blocks(out, bl, ["a1"], ["b1"])
    bl2, an, bn = load_blocks(out)
    assert bl2.dtype == BLOCK_DTYPE and an == ["a1"] and bn == ["b1"]
    assert bl2.tobytes() == bl.tobytes()


def test_multiple_batches(tmp_path, monkeypatch):
    monkeypatch.setattr(blk, "PAF_BATCH_LINES", 3)
    p = tmp_path / "x.paf"
    p.write_text("".join(
        f"a1\t1000\t{i}\t{i + 10}\t+\tb1\t5000\t{i * 2}\t{i * 2 + 10}\t9\t10\t60\n" for i in range(8)))
    bl = read_paf(p, {"a1": 0}, {"b1": 0})
    assert list(bl["a_start"]) == list(range(8))
    assert list(bl["b_start"]) == [i * 2 for i in range(8)]


# a1[1000,10301) on b1[2000,11496), =/X CIGAR: short indels (10I, 5D) and mismatches inside, a 500-bp deletion
# and a 300-bp insertion between the sub-blocks
_CG = "1000=1X2000=10I1000=5D990=500D2000=300I1500=2X498="
_CG_M = "3001M10I1000M5D990M500D2000M300I2000M"  # the same alignment with M ops
_REC = "a1\t20000\t1000\t10301\t{}\tb1\t20000\t2000\t11496\t8988\t9806\t60\ttp:A:P"
_COLS = ["a_chrom", "a_start", "a_end", "b_chrom", "b_start", "b_end", "strand", "matches", "alen"]


def _swapped(bl):
    """The blocks with the A and B fields exchanged (strand and counts unchanged)."""
    out = bl.copy()
    for s in ("chrom", "start", "end"):
        out[f"a_{s}"], out[f"b_{s}"] = bl[f"b_{s}"], bl[f"a_{s}"]
    return out


@pytest.mark.parametrize("strand, cg, want", [
    ("+", _CG, [(0, 1000, 6001, 0, 2000, 6996, 1, 4990, 5006), (0, 6001, 8001, 0, 7496, 9496, 1, 2000, 2000),
                (0, 8301, 10301, 0, 9496, 11496, 1, 1998, 2000)]),
    # '-': the CIGAR runs along B while A descends from its end
    ("-", _CG, [(0, 1000, 3000, 0, 9496, 11496, -1, 1998, 2000),
                (0, 3300, 5300, 0, 7496, 9496, -1, 2000, 2000),
                (0, 5300, 10301, 0, 2000, 6996, -1, 4990, 5006)]),
    # M ops: the record's matches prorated by sub-block alen
    ("+", _CG_M, [(0, 1000, 6001, 0, 2000, 6996, 1, 4996, 5006),
                  (0, 6001, 8001, 0, 7496, 9496, 1, 1996, 2000),
                  (0, 8301, 10301, 0, 9496, 11496, 1, 1996, 2000)]),
])
def test_read_paf_splits_at_long_indels(tmp_path, strand, cg, want):
    p = tmp_path / "x.paf"
    p.write_text(_REC.format(strand) + f"\tcg:Z:{cg}\n")
    bl = read_paf(p, {"a1": 0}, {"b1": 0})
    assert bl[_COLS].tolist() == want
    assert (bl["chain"] == -1).all() and np.isnan(bl["margin"]).all()
    assert read_paf(p, {"b1": 0}, {"a1": 0}, swap=True).tobytes() == _swapped(bl).tobytes()
    whole = [(0, 1000, 10301, 0, 2000, 11496, 1 if strand == "+" else -1, 8988, 9806)]
    assert read_paf(p, {"a1": 0}, {"b1": 0}, split_indel=501)[_COLS].tolist() == whole  # no indel that long
    p.write_text(_REC.format(strand) + "\n")  # no CIGAR: the record stays one block
    assert read_paf(p, {"a1": 0}, {"b1": 0})[_COLS].tolist() == whole


def test_read_paf_short_indels_and_bad_cigar(tmp_path, caplog):
    p = tmp_path / "x.paf"
    p.write_text("a1\t20000\t0\t3010\t+\tb1\t20000\t0\t3079\t3000\t3089\t60\tcg:Z:1000=10I1000=79D1000=\n"
                 # the CIGAR spans 1700 bp of B, the record 1400: kept whole, with a warning
                 "a1\t20000\t5000\t5800\t+\tb1\t20000\t5000\t6400\t800\t1400\t60\tcg:Z:400=900D400=\n")
    with caplog.at_level(logging.WARNING, logger=LOG):
        bl = read_paf(p, {"a1": 0}, {"b1": 0})
    assert bl[_COLS].tolist() == [(0, 0, 3010, 0, 0, 3079, 1, 3000, 3089),
                                  (0, 5000, 5800, 0, 5000, 6400, 1, 800, 1400)]
    assert "1 PAF record(s)" in caplog.text and "CIGAR" in caplog.text
