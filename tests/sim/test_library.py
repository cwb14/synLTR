"""The library writer on a toy genome with one plus- and one minus-strand element, and the screen reading what it
wrote."""
import logging
from pathlib import Path

import numpy as np
import pytest

from synltr.sim.pop.library import write_library

SEQ = "TTTTGATCCAGGTACCATGCAAGTCCGATTACAGGCATTTCCGGAAACTGTTGACCGTAA"  # chrT, 60 bp
GFF3 = """##gff-version 3
chrT\tLTRpav\tLTR_retrotransposon\t5\t24\t.\t+\t.\tID=e1;Name=e1;tsd=ACGTA;K2P_d=0.0123;superfamily=Gypsy
chrT\tLTRpav\tlong_terminal_repeat\t5\t9\t.\t+\t.\tID=e1.lLTR;Parent=e1
chrT\tLTRpav\tlong_terminal_repeat\t20\t24\t.\t+\t.\tID=e1.rLTR;Parent=e1
chrT\tLTRpav\tLTR_retrotransposon\t31\t54\t.\t-\t.\tID=e2;Name=e2
chrT\tLTRpav\tlong_terminal_repeat\t31\t35\t.\t-\t.\tID=e2.lLTR;Parent=e2
chrT\tLTRpav\tlong_terminal_repeat\t48\t54\t.\t-\t.\tID=e2.rLTR;Parent=e2
"""
PLUS = "GATCCAGGTACCATGCAAGT"        # chrT 5-24 as is
MINUS = "TCAACAGTTTCCGGAAATGCCTGT"   # chrT 31-54 (ACAGGCATTTCCGGAAACTGTTGA) reverse-complemented
NAMES = ("chrT:5-24#LTR/Gypsy", "chrT:31-54#LTR/unknown")


def _toy(tmp_path: Path, gff3: str = GFF3) -> tuple[Path, Path, Path]:
    gff, fa = tmp_path / "toy.gff3", tmp_path / "toy.fa"
    gff.write_text(gff3)
    fa.write_text(f">chrT\n{SEQ[:30]}\n{SEQ[30:]}\n")   # wrapped, so the .fai offsets are exercised
    return gff, fa, tmp_path / "lib" / "T"


def test_library(tmp_path, caplog):
    gff, fa, prefix = _toy(tmp_path)
    with caplog.at_level(logging.INFO, logger="synltr"):
        assert write_library(gff, fa, prefix) == 2
    assert "2 elements (1 minus-strand)" in caplog.text
    assert Path(f"{prefix}.fa").read_text() == f">{NAMES[0]}\n{PLUS}\n>{NAMES[1]}\n{MINUS}\n"
    assert Path(f"{prefix}.tsv").read_text().splitlines() == [
        "seq_id\tk2p\ttsd\tstrand\tltr5_len\tltr3_len",
        f"{NAMES[0]}\t0.0123\tACGTA\t+\t5\t5",
        f"{NAMES[1]}\tnan\t.\t-\t5\t7",       # ltr5 = genome-left LTR (31-35), ltr3 = genome-right (48-54)
    ]
    assert not list(prefix.parent.glob(".*.tmp"))


@pytest.mark.parametrize("bad, message", [
    (GFF3.replace("chrT\tLTRpav\tLTR_retrotransposon\t31", "chrX\tLTRpav\tLTR_retrotransposon\t31")
     .replace("chrT\tLTRpav\tlong_terminal_repeat\t31", "chrX\tLTRpav\tlong_terminal_repeat\t31")
     .replace("chrT\tLTRpav\tlong_terminal_repeat\t48", "chrX\tLTRpav\tlong_terminal_repeat\t48"),
     "seqid chrX is not in"),
    (GFF3.replace("\t31\t54\t", "\t31\t64\t").replace("\t48\t54\t", "\t58\t64\t"), "ends past its sequence"),
])
def test_fails_fast(tmp_path, bad, message):
    gff, fa, prefix = _toy(tmp_path, bad)
    with pytest.raises(ValueError, match=message):
        write_library(gff, fa, prefix)
    assert not Path(f"{prefix}.fa").exists() and not Path(f"{prefix}.tsv").exists()
    assert not list(prefix.parent.glob(".*.tmp"))


def test_screen_reads_library(tmp_path):
    """The screen reads the written library: names, coordinates, strand, K2P, canonical TSD, the LTR lengths
    (genome-left first) and the class in each name."""
    from synltr.io.library import read_library
    from synltr.routes.library import Records

    gff, fa, prefix = _toy(tmp_path)
    write_library(gff, fa, prefix)
    lib = read_library(f"{prefix}.fa,{prefix}.tsv", "T")
    rec = Records.of("T", lib)
    assert rec.names.tolist() == list(NAMES)
    assert (rec.chrom.tolist(), rec.start.tolist(), rec.end.tolist()) == (["chrT", "chrT"], [5, 31], [24, 54])
    assert rec.strand.tolist() == ["+", "-"]
    assert rec.k2p[0] == 0.0123 and np.isnan(rec.k2p[1])
    assert rec.tsd.tolist() == ["ACGTA", ""]                # canonical: min(t, reverse complement of t)
    assert rec.ltr_len.tolist() == [[5, 5], [5, 7]]
    assert list(lib.cls) == ["LTR/Gypsy", "LTR/unknown"]
    assert [rec.seq(i).size for i in range(2)] == [len(PLUS), len(MINUS)]
