#!/usr/bin/env python3
import sys
from collections import defaultdict

def parse_args():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} REF_ACCESSION INPUT_FILE", file=sys.stderr)
        sys.exit(1)
    return sys.argv[1], sys.argv[2]

def chr_sort_key(chr_name):
    """
    Given something like "Zmays_chr10" or "Salte_chrX", extract the part after the last '_'
    (e.g. "chr10" or "chrX"), strip the "chr" prefix, and try to convert to int.
    If it fails, fall back to sorting lexicographically after numeric ones.
    """
    try:
        suffix = chr_name.rsplit('_', 1)[1]  # e.g. "chr10" or "chrX"
        if suffix.lower().startswith("chr"):
            num_part = suffix[3:]
        else:
            num_part = suffix
        num = int(num_part)
        return (0, num)   # numeric chromosomes come first, sorted by number
    except Exception:
        return (1, chr_name)

def main():
    ref_accession, infile = parse_args()

    # For each (query_acc, query_chr), store its best match as (percent_overlap, reference_chr, strand)
    best_match = defaultdict(dict)

    # Track, for each query accession, which reference chromosomes ever appeared
    ref_chrs_seen = defaultdict(set)

    with open(infile) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split()
            # Expecting 10 columns:
            # 0:seq1 1:len1 2:sum1 3:frac1 4:seq2 5:len2 6:sum2 7:frac2 8:frac_plus 9:frac_minus
            if len(parts) < 10:
                continue

            seq1 = parts[0]
            seq2 = parts[4]

            # Robust numeric parsing
            try:
                frac1 = float(parts[3])
                frac2 = float(parts[7])
                frac_plus = float(parts[8])
                frac_minus = float(parts[9])
            except ValueError:
                continue

            # Choose strand by which orientation has more synteny (ties -> '+')
            strand = '+' if frac_plus >= frac_minus else '-'

            acc1 = seq1.rsplit('_', 1)[0]
            acc2 = seq2.rsplit('_', 1)[0]

            # Determine which side is reference and which is query
            if acc1 == ref_accession and acc2 != ref_accession:
                reference_chr = seq1
                query_acc = acc2
                query_chr = seq2
                pct_for_query = frac2  # fraction of seq2 covered by seq1
            elif acc2 == ref_accession and acc1 != ref_accession:
                reference_chr = seq2
                query_acc = acc1
                query_chr = seq1
                pct_for_query = frac1  # fraction of seq1 covered by seq2
            else:
                # Either both are the ref accession or neither; skip
                continue

            # Remember that this reference_chr was seen for this query accession
            ref_chrs_seen[query_acc].add(reference_chr)

            # Update best match for (query_acc, query_chr) if this pct is higher
            prev = best_match[query_acc].get(query_chr)
            if (prev is None) or (pct_for_query > prev[0]):
                best_match[query_acc][query_chr] = (pct_for_query, reference_chr, strand)

    # Now produce output blocks, one per query accession (in sorted order of query_acc)
    for i, query_acc in enumerate(sorted(ref_chrs_seen.keys())):
        # Group by reference_chr: each ref_chr maps to a list of (query_chr, strand, pct)
        grouped_by_ref = defaultdict(list)
        for query_chr, (pct, ref_chr, strand) in best_match[query_acc].items():
            grouped_by_ref[ref_chr].append((query_chr, strand, pct))

        # Ensure every reference chromosome that ever appeared is present
        for ref_chr in ref_chrs_seen[query_acc]:
            grouped_by_ref.setdefault(ref_chr, [])

        # Sort the reference chromosomes (e.g. Zmays_chr1, Zmays_chr2, …)
        sorted_ref_chrs = sorted(ref_chrs_seen[query_acc], key=chr_sort_key)

        # Print each line: if there are no hits, print "<ref_chr>\tNA"
        # Otherwise: "<ref_chr>\t<q_chr1>\t<strand1>\t<q_chr2>\t<strand2> …" (best-first by pct)
        for ref_chr in sorted_ref_chrs:
            hits = grouped_by_ref[ref_chr]
            if not hits:
                print(f"{ref_chr}\tNA")
            else:
                hits.sort(key=lambda x: x[2], reverse=True)  # by pct desc
                out = [ref_chr]
                for (q_chr, strand, _pct) in hits:
                    out.append(q_chr)
                    out.append(strand)
                print("\t".join(out))

        # Blank line between query‐accession blocks (but not after the last one)
        if i < len(ref_chrs_seen) - 1:
            print()

if __name__ == "__main__":
    main()
