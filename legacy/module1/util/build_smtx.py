#!/usr/bin/env python3
"""
Build a concatenated supermatrix (FASTA) from many gene alignments and
write IQ-TREE partition files (per-gene and per-codon).

Assumptions:
- Each input FASTA has headers like '>Bdact002821' where the taxon is the
  (non-digit) prefix before the first digit.
- Sequences are back-translated codon alignments (length should be multiple of 3).
- Files share the same set of taxa (if some are missing in a file, they'll be padded with gaps).
"""

import argparse, re, sys
from pathlib import Path

def parse_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("-i", "--input_glob", default="*.codon.trim.fna",
                   help="Glob pattern for input FASTA files")
    p.add_argument("-o", "--out_prefix", default="supermatrix",
                   help="Output prefix")
    p.add_argument("--taxon_regex", default=r"^>([A-Za-z_][A-Za-z0-9_]*?)\d",
                   help="Regex to capture taxon name from FASTA header; group(1) is the taxon")
    p.add_argument("--strict", action="store_true",
                   help="Abort if any gene length is not a multiple of 3")
    return p.parse_args()

def read_fasta(path):
    """Return list of (header, seq) tuples; header includes initial '>'."""
    records = []
    h, parts = None, []
    with open(path, "r") as fh:
        for line in fh:
            line = line.rstrip()
            if not line: 
                continue
            if line.startswith(">"):
                if h is not None:
                    records.append((h, "".join(parts)))
                h, parts = line, []
            else:
                parts.append(line.replace(" ", ""))
        if h is not None:
            records.append((h, "".join(parts)))
    return records

def main():
    args = parse_args()
    files = sorted(Path(".").glob(args.input_glob),
                   key=lambda p: [int(x) if x.isdigit() else x.lower()
                                  for x in re.findall(r"\d+|[A-Za-z]+", p.name)])
    if not files:
        sys.exit(f"No files matched: {args.input_glob}")

    taxon_pat = re.compile(args.taxon_regex)

    # First pass: collect taxa set and per-file (gene) sequences
    genes = []  # list of dicts: {"name": <basename>, "len": L, "seqs": {taxon: seq}}
    all_taxa = set()

    for f in files:
        recs = read_fasta(f)
        if not recs:
            sys.exit(f"Empty FASTA: {f}")
        seqs = {}
        for h, s in recs:
            m = taxon_pat.match(h)
            if not m:
                sys.exit(f"Could not parse taxon from header '{h}' in {f}\n"
                         f"Adjust --taxon_regex if needed.")
            taxon = m.group(1)
            seqs.setdefault(taxon, s)
            all_taxa.add(taxon)
        # sanity: all sequences same length
        lengths = {len(s) for s in seqs.values()}
        if len(lengths) != 1:
            sys.exit(f"Non-uniform sequence lengths in {f}: {sorted(lengths)}")
        L = lengths.pop()
        if L % 3 != 0:
            msg = f"WARNING: {f.name} length {L} not divisible by 3."
            if args.strict:
                sys.exit(msg)
            else:
                print(msg, file=sys.stderr)
        genes.append({"name": f.stem, "len": L, "seqs": seqs})

    # Stable taxon order (alphabetical)
    taxa = sorted(all_taxa)

    # Concatenate sequences per taxon, padding gaps where missing
    concat = {t: [] for t in taxa}
    bounds = []  # (gene_name, start, end) 1-based inclusive
    pos = 1
    for g in genes:
        L = g["len"]
        for t in taxa:
            seq = g["seqs"].get(t, "-" * L)
            if len(seq) != L:
                sys.exit(f"Internal error: {g['name']} taxon {t} length {len(seq)} != {L}")
            concat[t].append(seq)
        bounds.append((g["name"], pos, pos + L - 1))
        pos += L

    # Write supermatrix FASTA
    out_fa = Path(f"{args.out_prefix}.fna")
    with open(out_fa, "w") as out:
        for t in taxa:
            s = "".join(concat[t])
            out.write(f">{t}\n")
            # wrap to 80 cols
            for i in range(0, len(s), 80):
                out.write(s[i:i+80] + "\n")

    # Write per-gene partition file (RAxML/IQ-TREE format)
    part_gene = Path(f"{args.out_prefix}.partitions.genes.txt")
    with open(part_gene, "w") as out:
        for name, start, end in bounds:
            out.write(f"DNA, {name} = {start}-{end}\n")

    # Write codon-position partitions across each gene
    part_codon = Path(f"{args.out_prefix}.partitions.codon.txt")
    with open(part_codon, "w") as out:
        for name, start, end in bounds:
            # IQ-TREE/RAxML step pattern: "start-end\3" means every 3rd site
            out.write(f"DNA, {name}_pos1 = {start}-{end}\\3\n")
            out.write(f"DNA, {name}_pos2 = {start+1}-{end}\\3\n")
            out.write(f"DNA, {name}_pos3 = {start+2}-{end}\\3\n")

    # Also drop a simple report
    report = Path(f"{args.out_prefix}.summary.tsv")
    with open(report, "w") as out:
        out.write("gene\tlength\tstart\tend\n")
        for name, start, end in bounds:
            out.write(f"{name}\t{end-start+1}\t{start}\t{end}\n")

    total_len = sum(g["len"] for g in genes)
    print(f"Wrote supermatrix: {out_fa} ({len(taxa)} taxa × {total_len} sites)")
    print(f"Per-gene partitions: {part_gene}")
    print(f"Codon-position partitions: {part_codon}")
    print(f"Gene summary: {report}")
    print("Tip: If any taxa were missing for some genes, those regions were filled with gaps ('-').")

if __name__ == "__main__":
    main()
