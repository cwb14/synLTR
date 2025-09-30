#!/usr/bin/env python3
"""
Script to polish TSV alignments by removing total overlaps and trimming partial overlaps.
Usage:
    python polish_alignments.py input.tsv output.tsv
"""
import argparse
from collections import defaultdict

def parse_alignment(aln_str):
    # Parse strings like 'Accession_chrNUM:start..end'
    acc_chr, coords = aln_str.split(':')
    start_str, end_str = coords.split('..')
    return acc_chr, int(start_str), int(end_str)

def format_alignment(acc_chr, start, end):
    return f"{acc_chr}:{start}..{end}"

def main():
    parser = argparse.ArgumentParser(description="Polish TSV alignments")
    parser.add_argument('input', help='Input TSV file')
    parser.add_argument('output', help='Output TSV file')
    args = parser.parse_args()

    # Read alignments
    alignments = []
    with open(args.input) as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            ref_str, qry_str, strand = line.split('\t')
            ref_acc, ref_start, ref_end = parse_alignment(ref_str)
            qry_acc, qry_start, qry_end = parse_alignment(qry_str)
            alignments.append({
                'line': line,
                'ref_str': ref_str, 'qry_str': qry_str,
                'strand': strand,
                'ref_acc': ref_acc, 'ref_start': ref_start, 'ref_end': ref_end,
                'qry_acc': qry_acc, 'qry_start': qry_start, 'qry_end': qry_end
            })

    # Group by reference, query, and strand
    groups = defaultdict(list)
    for aln in alignments:
        key = (aln['ref_acc'], aln['qry_acc'], aln['strand'])
        groups[key].append(aln)

    polished = []
    # Process each group separately
    for key, group in groups.items():
        g = list(group)
        removed_indices = set()

        # 1) TOTAL OVERLAP: remove contained alignments
        for i in range(len(g)):
            for j in range(len(g)):
                if i == j or j in removed_indices:
                    continue
                ai = g[i]
                aj = g[j]
                # Both share exact key, so just check containment
                if (ai['ref_start'] <= aj['ref_start'] <= aj['ref_end'] <= ai['ref_end'] and
                    ai['qry_start'] <= aj['qry_start'] <= aj['qry_end'] <= ai['qry_end']):
                    print(f"TOTAL OVERLAP: {ai['line']} totally encompasses {aj['line']}. {aj['line']} removed.")
                    removed_indices.add(j)

        # Build reduced group after removals
        reduced = [g[k] for k in range(len(g)) if k not in removed_indices]

        # 2) PARTIAL OVERLAP: trim the later-start alignment
        for i in range(len(reduced)):
            for j in range(len(reduced)):
                if i == j:
                    continue
                ai = reduced[i]
                aj = reduced[j]
                # Identify partial overlap for ref and query
                if (ai['ref_start'] < aj['ref_start'] < ai['ref_end'] < aj['ref_end'] and
                    ai['qry_start'] < aj['qry_start'] < ai['qry_end'] < aj['qry_end']):
                    # ai is earlier alignment, aj is later
                    new_ref_start = ai['ref_end']
                    new_qry_start = ai['qry_end']
                    old_line = aj['line']
                    # Create trimmed strings
                    new_ref_str = format_alignment(aj['ref_acc'], new_ref_start, aj['ref_end'])
                    new_qry_str = format_alignment(aj['qry_acc'], new_qry_start, aj['qry_end'])
                    new_line = f"{new_ref_str}\t{new_qry_str}\t{aj['strand']}"
                    print(f"PARTIAL OVERLAP: {ai['line']} partially overlaps with {old_line}. {old_line} trimmed to {new_line}")
                    # Update aj in reduced
                    aj['ref_str']    = new_ref_str
                    aj['qry_str']    = new_qry_str
                    aj['ref_start']  = new_ref_start
                    aj['qry_start']  = new_qry_start
                    aj['line']       = new_line

        # Collect polished lines for this group
        for aln in reduced:
            polished.append(aln['line'])

    # Write the polished alignments
    with open(args.output, 'w') as outf:
        for ln in polished:
            outf.write(ln + '\n')

if __name__ == '__main__':
    main()
