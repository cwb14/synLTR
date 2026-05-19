#!/usr/bin/env python3
"""Convert minimap2 PAF (run with -c) to the mmseqs format-mode-4 m8 that
filter_solo_ltr.py consumes, so the filter is aligner-agnostic.

Emitted columns (tab, with a header line, identical to mmseqs --format-mode 4
--format-output query,target,qstart,qend,qlen,tstart,tend,tlen,alnlen,pident,
mismatch,gapopen,evalue,bits,qcov):

  query target qstart qend qlen tstart tend tlen alnlen pident mismatch
  gapopen evalue bits qcov

Conventions made identical to the mmseqs path so the strand-aware filter needs
no change:
  * coords 1-based inclusive (PAF is 0-based half-open)
  * MINUS-strand hit -> query coords reversed (qstart > qend), target coords
    ascending  (exactly the mmseqs --search-type 3 convention)
  * bits  := minimap2 DP alignment score (AS:i: tag) so merge_loci's
    "best hit = max bits" picks the strongest alignment
  * evalue:= minimap2 gap-compressed divergence (de:f: / dv:f:) if present
    else 0.0 (informational only; the filter does not gate on it)
  * pident:= 100 * nmatch / alnblocklen   (PAF col10 / col11)

Reads stdin or a file; writes stdout or -o FILE.
"""
import argparse
import sys

COLS = ("query\ttarget\tqstart\tqend\tqlen\ttstart\ttend\ttlen\talnlen\t"
        "pident\tmismatch\tgapopen\tevalue\tbits\tqcov")


def tag(fields, key, cast, default):
    pref = key + ":"
    for f in fields:
        if f.startswith(pref):
            try:
                return cast(f.split(":", 2)[2])
            except (ValueError, IndexError):
                return default
    return default


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("paf", nargs="?", default="-")
    ap.add_argument("-o", "--out", default="-")
    a = ap.parse_args()
    fin = sys.stdin if a.paf == "-" else open(a.paf)
    fout = sys.stdout if a.out == "-" else open(a.out, "w")
    fout.write(COLS + "\n")
    n = 0
    for line in fin:
        if not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 12:
            continue
        qname, qlen, qs, qe, strand, tname, tlen, ts, te, nmatch, ablen, _mapq \
            = f[:12]
        qlen = int(qlen); qs = int(qs); qe = int(qe)
        tlen = int(tlen); ts = int(ts); te = int(te)
        nmatch = int(nmatch); ablen = int(ablen) or (qe - qs) or 1
        tags = f[12:]
        AS = tag(tags, "AS", int, nmatch)            # DP score -> "bits"
        de = tag(tags, "de", float, tag(tags, "dv", float, 0.0))
        pident = round(100.0 * nmatch / ablen, 3)
        qcov = round((qe - qs) / qlen, 4) if qlen else 0.0
        if strand == "+":
            qS, qE = qs + 1, qe                      # 1-based inclusive
        else:                                        # minus: reverse query
            qS, qE = qe, qs + 1                      # (qstart > qend)
        tS, tE = ts + 1, te                          # target always ascending
        fout.write(
            f"{qname}\t{tname}\t{qS}\t{qE}\t{qlen}\t{tS}\t{tE}\t{tlen}\t"
            f"{ablen}\t{pident}\t0\t0\t{de:.3g}\t{AS}\t{qcov}\n")
        n += 1
    if fin is not sys.stdin:
        fin.close()
    if fout is not sys.stdout:
        fout.close()
    sys.stderr.write(f"[paf2m8] wrote {n} alignments\n")


if __name__ == "__main__":
    main()
