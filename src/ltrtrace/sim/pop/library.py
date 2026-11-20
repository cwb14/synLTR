"""Write a simulated genome's LTR-RT GFF3 (as `ltrtrace sim pop` writes it) as a generic LTR-RT library, so the
library route can be scored on the same truth:
  PREFIX.fa   each element's LTR-to-LTR span, header "seqid:start-end#LTR/<superfamily or unknown>" (1-based
              closed), minus-strand elements reverse-complemented (as LTRquest writes them), one line each
  PREFIX.tsv  seq_id k2p tsd strand ltr5_len ltr3_len: k2p from K2P_d (nan if absent, %.6g), tsd '.' unless
              ACGT only, strand '.' unless + or -, ltr5_len the genome-left LTR, ltr3_len the genome-right LTR
The GFF3 is read by the package's own reader, so the library and the assembly route see the same elements. Each
span is fetched through the FASTA's .fai index; the FASTA is never read whole. Fails fast, writing nothing, on an
element off its sequence, a seqid missing from the FASTA, or a duplicate name.
"""
from __future__ import annotations

import logging
import os
import re
from pathlib import Path

import pysam

from ltrtrace.io.annotation.gff3 import load_gff3

log = logging.getLogger(__name__)

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
_ACGT = re.compile(r"[ACGTacgt]+")


def write_library(gff3: Path, fasta: Path, prefix: Path) -> int:
    """Write PREFIX.fa and PREFIX.tsv for the elements of `gff3` in `fasta`; the number of elements written.
    ValueError on an unreadable FASTA, an element off its sequence, or a duplicate element."""
    elements, skips = load_gff3(gff3, "lib")
    if skips:
        log.debug("%s: skipped %s", gff3, dict(skips))
    try:
        fa = pysam.FastaFile(str(fasta))
    except (OSError, ValueError) as e:
        raise ValueError(f"cannot open {fasta}: {e}") from None
    lengths = dict(zip(fa.references, fa.lengths, strict=True))
    prefix.parent.mkdir(parents=True, exist_ok=True)
    fa_tmp = prefix.with_name(f".{prefix.name}.fa.{os.getpid()}.tmp")
    tb_tmp = prefix.with_name(f".{prefix.name}.tsv.{os.getpid()}.tmp")
    seen, n_minus = set(), 0
    try:
        with fa, open(fa_tmp, "w") as out_fa, open(tb_tmp, "w") as out_tb:
            out_tb.write("seq_id\tk2p\ttsd\tstrand\tltr5_len\tltr3_len\n")
            for e in elements:
                if e.seqid not in lengths:
                    raise ValueError(f"{gff3}: seqid {e.seqid} is not in {fasta}")
                if e.end > lengths[e.seqid]:
                    raise ValueError(f"{gff3}: element {e.seqid}:{e.start}-{e.end} ends past its sequence "
                                     f"({lengths[e.seqid]} bp)")
                name = f"{e.seqid}:{e.start}-{e.end}#LTR/{e.superfamily or 'unknown'}"
                if name in seen:
                    raise ValueError(f"{gff3}: duplicate element {name}")
                seen.add(name)
                s = fa.fetch(e.seqid, e.start - 1, e.end)
                if e.strand == "-":
                    s = s.translate(_COMP)[::-1]
                    n_minus += 1
                l5 = e.ltr_left[1] - e.ltr_left[0] + 1
                l3 = e.ltr_right[1] - e.ltr_right[0] + 1
                k2p = f"{e.k:.6g}" if e.k is not None else "nan"
                tsd = e.tsd if e.tsd and _ACGT.fullmatch(e.tsd) else "."
                strand = e.strand if e.strand in ("+", "-") else "."
                out_fa.write(f">{name}\n{s}\n")
                out_tb.write(f"{name}\t{k2p}\t{tsd}\t{strand}\t{l5}\t{l3}\n")
        os.replace(fa_tmp, f"{prefix}.fa")
        os.replace(tb_tmp, f"{prefix}.tsv")
    finally:
        fa_tmp.unlink(missing_ok=True)
        tb_tmp.unlink(missing_ok=True)
    log.info("%s: %d elements (%d minus-strand)", prefix, len(elements), n_minus)
    return len(elements)
