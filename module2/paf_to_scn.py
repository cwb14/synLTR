#!/usr/bin/env python3
import sys

def parse_paf_line(line):
    # PAF required fields (first 12):
    # 0:qname 1:qlen 2:qstart 3:qend 4:strand 5:tname 6:tlen 7:tstart 8:tend 9:nmatch 10:alen 11:mapq ...
    f = line.rstrip("\n").split("\t")
    if len(f) < 12:
        return None
    try:
        qname = f[0]
        qlen  = int(f[1])
        qstart = int(f[2])
        qend   = int(f[3])
        strand = f[4]
        tname = f[5]
        tlen  = int(f[6])
        tstart = int(f[7])
        tend   = int(f[8])
    except ValueError:
        return None
    return {
        "qname": qname, "qlen": qlen, "qstart": qstart, "qend": qend,
        "strand": strand, "tname": tname, "tlen": tlen, "tstart": tstart, "tend": tend
    }

def scn_from_pair(chrom, a_start, a_end, b_start, b_end):
    # Left = smaller start; Right = the other
    if (a_start, a_end) <= (b_start, b_end):
        Ls, Le = a_start, a_end
        Rs, Re = b_start, b_end
    else:
        Ls, Le = b_start, b_end
        Rs, Re = a_start, a_end

    span = Re - Ls
    llen = Le - Ls
    rlen = Re - Rs
    # SCN columns (per your example):
    # [0]span_start [1]span_end [2]span_len [3]left_start [4]left_end [5]left_len
    # [6]right_start [7]right_end [8]right_len [9]NA [10]NA [11]chrom
    return [str(Ls), str(Re), str(span),
            str(Ls), str(Le), str(llen),
            str(Rs), str(Re), str(rlen),
            "NA", "NA", chrom]

def main():
    import argparse
    p = argparse.ArgumentParser(description="Convert PAF LTR-RT candidate alignments to SCN format, deduplicating reciprocal pairs.")
    p.add_argument("paf", nargs="?", help="PAF file (default: stdin)")
    args = p.parse_args()

    fh = sys.stdin if not args.paf or args.paf == "-" else open(args.paf, "r")

    seen = set()
    out = sys.stdout

    for line in fh:
        if not line.strip() or line.startswith("#"):
            continue
        rec = parse_paf_line(line)
        if rec is None:
            continue

        # Only consider intra-chromosomal pairs (your example uses same chr)
        if rec["qname"] != rec["tname"]:
            continue

        chrom = rec["qname"]
        a = (rec["qstart"], rec["qend"])
        b = (rec["tstart"], rec["tend"])

        # Canonicalize an unordered pair key to drop reciprocal duplicates
        key = (chrom, tuple(sorted([a, b])))
        if key in seen:
            continue
        seen.add(key)

        row = scn_from_pair(chrom, a[0], a[1], b[0], b[1])
        out.write(" ".join(row) + "\n")

    if fh is not sys.stdin:
        fh.close()

if __name__ == "__main__":
    main()
