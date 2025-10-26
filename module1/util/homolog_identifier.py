#!/usr/bin/env python3
import argparse
import os
import re
import subprocess
from collections import defaultdict

# =============================================================================
# Developer notes (strand fraction toggle)
# -----------------------------------------------------------------------------
# By default, we compute plus/minus strand fractions on the LEFT (Seq1) side.
# To switch back to the RIGHT (Seq2) side (original v1.3 behavior), just change:
#
#     MEASURE_STRAND_ON = 'left'   ->   MEASURE_STRAND_ON = 'right'
#
# The denominator used for plus/minus is the UNIQUE bp covered by (+) and (-)
# on the chosen side (i.e., union_plus + union_minus). If you prefer to use the
# total unique coverage on that side (including segments not represented in +/-
# lists), you can replace the denominator with `first_cov` or `second_cov`
# depending on which side you're measuring (see the comment where computed).
# =============================================================================
MEASURE_STRAND_ON = 'left'  # options: 'left' or 'right'

# --------- Helpers ---------

INTERVAL_RE = re.compile(r'^([^:]+):(\d+)(?:\.\.|\-)(\d+)$')

def parse_interval(s):
    """
    Parse 'Name:start..end' or 'Name:start-end' into (name, start, end)
    Returns start, end as integers with the convention that length = end - start (half-open).
    """
    m = INTERVAL_RE.match(s.strip())
    if not m:
        raise ValueError(f"Bad interval format: {s!r}")
    name, a, b = m.group(1), int(m.group(2)), int(m.group(3))
    # Normalize to half-open [min, max) regardless of input ordering
    start, end = (a, b) if a <= b else (b, a)
    return name, start, end

def accession_of(seqname):
    """
    Extract accessionID as the substring before the first underscore.
    Example: 'Bdact_chr1' -> 'Bdact'
    """
    return seqname.split('_', 1)[0]

def ensure_fai(genome_fa):
    """
    Ensure a .fai index exists for genome_fa. If missing, run samtools faidx.
    """
    fai = genome_fa + ".fai"
    if not os.path.exists(fai):
        # run samtools faidx genome_fa
        try:
            subprocess.run(["samtools", "faidx", genome_fa], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except FileNotFoundError:
            raise RuntimeError("samtools not found in PATH. Please install samtools or add it to PATH.")
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"samtools faidx failed for {genome_fa}: {e.stderr.decode(errors='ignore')}")
    return fai

def load_lengths_for_accession(genome_dir, accession):
    """
    Load {seq_name: length} from accession.fa.fai, generating it if needed.
    Caches results across calls via function attribute.
    """
    if not hasattr(load_lengths_for_accession, "_cache"):
        load_lengths_for_accession._cache = {}
    cache = load_lengths_for_accession._cache

    if accession in cache:
        return cache[accession]

    fa_path = os.path.join(genome_dir, f"{accession}.fa")
    if not os.path.exists(fa_path):
        raise FileNotFoundError(f"Genome FASTA not found: {fa_path}")

    fai_path = ensure_fai(fa_path)
    lens = {}
    with open(fai_path, 'r') as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            seqname = parts[0]
            length = int(parts[1])
            lens[seqname] = length

    cache[accession] = lens
    return lens

def union_length(intervals):
    """
    Given a list of (start, end) half-open intervals (start < end),
    return the total length of their union.
    """
    if not intervals:
        return 0
    # sort by start, then end
    intervals = sorted(intervals)
    merged = []
    cur_s, cur_e = intervals[0]
    for s, e in intervals[1:]:
        if s <= cur_e:  # overlap or touch
            if e > cur_e:
                cur_e = e
        else:
            merged.append((cur_s, cur_e))
            cur_s, cur_e = s, e
    merged.append((cur_s, cur_e))
    return sum(e - s for s, e in merged)

# --------- Main processing ---------

def compute_coverage(input_path, genome_dir):
    """
    Reads the synteny file and computes coverage per pair of sequences.
    Returns a list of output rows as tuples of strings/numbers.
    """
    # For each pair key = (seq1, seq2), collect intervals on each side
    per_pair_first = defaultdict(list)    # (seq1, seq2) -> list[(s1, e1)]
    per_pair_second = defaultdict(list)   # (seq1, seq2) -> list[(s2, e2)]

    # Strand-specific intervals on EACH side (so we can choose where to measure)
    per_pair_first_plus = defaultdict(list)    # '+' intervals on seq1
    per_pair_first_minus = defaultdict(list)   # '-' intervals on seq1
    per_pair_second_plus = defaultdict(list)   # '+' intervals on seq2
    per_pair_second_minus = defaultdict(list)  # '-' intervals on seq2

    # We'll also collect all unique sequence names encountered to fetch lengths
    seqnames_needed = set()

    with open(input_path, 'r') as f:
        for ln, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = re.split(r'\t+', line)
            if len(parts) < 3:
                raise ValueError(f"Line {ln}: expected 3 tab-separated fields, got: {line}")
            col1, col2, strand = parts[0], parts[1], parts[2].strip()
            if strand not in ('+', '-'):
                raise ValueError(f"Line {ln}: strand must be '+' or '-', got: {strand!r}")

            name1, s1, e1 = parse_interval(col1)
            name2, s2, e2 = parse_interval(col2)

            key = (name1, name2)
            per_pair_first[key].append((s1, e1))
            per_pair_second[key].append((s2, e2))

            # Record strand-specific intervals on BOTH sides so we can choose later
            if strand == '+':
                per_pair_first_plus[key].append((s1, e1))
                per_pair_second_plus[key].append((s2, e2))
            else:
                per_pair_first_minus[key].append((s1, e1))
                per_pair_second_minus[key].append((s2, e2))

            seqnames_needed.add(name1)
            seqnames_needed.add(name2)

    # Resolve lengths for all needed sequence names via their accessions
    seq_len = {}
    # Group needed sequence names by accession
    by_acc = defaultdict(list)
    for seq in seqnames_needed:
        by_acc[accession_of(seq)].append(seq)

    # Load lens from each accession and populate seq_len
    for acc, seqs in by_acc.items():
        lens = load_lengths_for_accession(genome_dir, acc)
        for seq in seqs:
            if seq not in lens:
                raise KeyError(f"Sequence {seq!r} not found in {acc}.fa.fai")
            seq_len[seq] = lens[seq]

    # Build results
    rows = []
    for (seq1, seq2) in sorted(per_pair_first.keys()):
        first_intervals = per_pair_first[(seq1, seq2)]
        second_intervals = per_pair_second[(seq1, seq2)]

        first_cov = union_length(first_intervals)
        second_cov = union_length(second_intervals)

        # Percentages relative to full sequence lengths
        len1 = seq_len[seq1]
        len2 = seq_len[seq2]
        pct1 = (first_cov / len1) if len1 > 0 else 0.0
        pct2 = (second_cov / len2) if len2 > 0 else 0.0

        # ---------------------------------------------------------------------
        # Strand breakdown (toggle via MEASURE_STRAND_ON above)
        # We compute unique bp for '+' and '-' on the chosen side and divide by
        # (plus_cov + minus_cov) on that same side.
        # If you prefer to use total unique coverage on that side as denom,
        # replace 'denom = plus_cov + minus_cov' with:
        #   denom = first_cov  (when MEASURE_STRAND_ON == 'left')
        #   denom = second_cov (when MEASURE_STRAND_ON == 'right')
        # ---------------------------------------------------------------------
        if MEASURE_STRAND_ON == 'left':
            plus_cov = union_length(per_pair_first_plus[(seq1, seq2)])
            minus_cov = union_length(per_pair_first_minus[(seq1, seq2)])
            denom = plus_cov + minus_cov  # or: denom = first_cov
        elif MEASURE_STRAND_ON == 'right':
            plus_cov = union_length(per_pair_second_plus[(seq1, seq2)])
            minus_cov = union_length(per_pair_second_minus[(seq1, seq2)])
            denom = plus_cov + minus_cov  # or: denom = second_cov
        else:
            raise ValueError(f"Invalid MEASURE_STRAND_ON: {MEASURE_STRAND_ON!r}; use 'left' or 'right'.")

        pct_plus = (plus_cov / denom) if denom > 0 else 0.0
        pct_minus = (minus_cov / denom) if denom > 0 else 0.0

        rows.append((
            seq1, len1, first_cov, pct1,
            seq2, len2, second_cov, pct2,
            pct_plus, pct_minus
        ))
    return rows

def main():
    ap = argparse.ArgumentParser(description="""
Compute per-chromosome-pair syntenic coverage from a block file.

Input format (tab-separated, 3 columns per line):
  col1: 'Seq1:start..end' (Seq1 is always considered plus strand)
  col2: 'Seq2:start..end'
  col3: strand of col2 ('+' or '-')

Assumptions:
  * Coordinates are treated as half-open intervals: length = end - start.
  * Unique positions only: overlaps within a pair are deduplicated via interval unions.
  * Strand percentages (cols 9 & 10) are based on unique coverage on the
    side chosen by MEASURE_STRAND_ON (developer toggle at top).

Genome discovery:
  * For sequence 'Bdact_chr1', accessionID = 'Bdact' (substring before first '_').
  * Script looks for '[accessionID].fa' in --genome-dir (default: current dir).
  * Ensures .fai (runs 'samtools faidx' if missing) and reads lengths from it.
""")
    ap.add_argument('-i', '--input', required=True, help='Syntenic block TSV file')
    ap.add_argument('-g', '--genome-dir', default='.', help='Directory with [accessionID].fa files (default: current directory)')
    args = ap.parse_args()

    rows = compute_coverage(args.input, args.genome_dir)
    # Output header (optional; comment it if you want headerless output)
    # print("\t".join(["seq1","len1","sum1","pct1","seq2","len2","sum2","pct2","pct_plus","pct_minus"]))
    for r in rows:
        # Format: 10 columns, tab-separated
        # Columns:
        # 1 seq1, 2 len1, 3 sum1, 4 pct1, 5 seq2, 6 len2, 7 sum2, 8 pct2, 9 pct_plus, 10 pct_minus
        # Use full precision for integers; 10 decimal places for fractions (as in your example)
        out = [
            str(r[0]),
            str(r[1]),
            str(r[2]),
            f"{r[3]:.10f}",
            str(r[4]),
            str(r[5]),
            str(r[6]),
            f"{r[7]:.10f}",
            f"{r[8]:.10f}",
            f"{r[9]:.10f}",
        ]
        print("\t".join(out))

if __name__ == "__main__":
    main()
