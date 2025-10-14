#!/usr/bin/env python3
import sys, argparse, gzip, re

def open_maybe_gzip(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)

def load_fasta_sequences(fasta_path):
    """
    Load a (possibly gzipped) FASTA into a dict: {header_without_gt: uppercased_sequence}.
    Concatenates wrapped lines.
    """
    seqs = {}
    name = None
    buf = []
    with open_maybe_gzip(fasta_path, "rt") as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf).upper()
                name = line[1:].strip().split()[0]  # take first token as the ID
                buf = []
            else:
                buf.append(line.strip())
        if name is not None:
            seqs[name] = "".join(buf).upper()
    return seqs

def parse_locus(locus):
    """
    Parse 'Chr1:269749-278390' -> ('Chr1', 269749, 278390)
    Coordinates are expected 1-based inclusive.
    """
    m = re.match(r'^([^:\s]+):(\d+)-(\d+)$', locus)
    if not m:
        raise ValueError(f"Cannot parse locus '{locus}' (expected Chr:START-END).")
    chrom, start, end = m.group(1), int(m.group(2)), int(m.group(3))
    if end < start:
        raise ValueError(f"End < start in locus '{locus}'.")
    return chrom, start, end

def extract(seqs, chrom, start_1based, end_1based):
    """
    Return subsequence using 1-based inclusive coordinates.
    """
    if chrom not in seqs:
        raise KeyError(f"Chromosome '{chrom}' not found in FASTA.")
    s = seqs[chrom]
    # convert to 0-based slice indices
    start0 = start_1based - 1
    end_excl = end_1based  # because Python slices are end-exclusive on 0-based
    if start0 < 0 or end_excl > len(s):
        raise IndexError(
            f"Requested range {chrom}:{start_1based}-{end_1based} is out of bounds "
            f"(sequence length {len(s)})."
        )
    return s[start0:end_excl]

def main():
    ap = argparse.ArgumentParser(
        description="Extract sequences for loci in column 1 of a TSV and print FASTA to stdout."
    )
    ap.add_argument("fasta", help="Reference genome FASTA (optionally .gz)")
    ap.add_argument("tsv", help="Input TSV with loci like 'Chr1:269749-278390' in column 1 (optionally .gz)")
    ap.add_argument("--col", type=int, default=1,
                    help="1-based column index with loci (default: 1)")
    args = ap.parse_args()

    try:
        seqs = load_fasta_sequences(args.fasta)
    except Exception as e:
        print(f"[error] Failed to read FASTA: {e}", file=sys.stderr)
        sys.exit(1)

    col_idx0 = args.col - 1

    n_ok = n_bad = 0
    with open_maybe_gzip(args.tsv, "rt") as fh:
        for ln, line in enumerate(fh, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()  # splits on any whitespace (tabs or spaces)
            if col_idx0 >= len(fields):
                print(f"[warn] Line {ln}: not enough columns, skipping.", file=sys.stderr)
                n_bad += 1
                continue
            locus = fields[col_idx0]
            try:
                chrom, start, end = parse_locus(locus)
                seq = extract(seqs, chrom, start, end)
            except Exception as e:
                print(f"[warn] Line {ln} ({locus}): {e}", file=sys.stderr)
                n_bad += 1
                continue
            # Write FASTA to stdout
            sys.stdout.write(f">{chrom}:{start}-{end}\n")
            # wrap at 60 columns for readability
            for i in range(0, len(seq), 60):
                sys.stdout.write(seq[i:i+60] + "\n")
            n_ok += 1

    print(f"[info] Extracted {n_ok} sequences; {n_bad} skipped.", file=sys.stderr)

if __name__ == "__main__":
    main()
