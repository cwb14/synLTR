#!/usr/bin/env python3
import sys

if len(sys.argv) < 3:
    sys.stderr.write(f"Usage: {sys.argv[0]} file1 file2 > output.tsv\n")
    sys.exit(1)

file1, file2 = sys.argv[1], sys.argv[2]

# Build key -> (tsd, motif) from file2 where key is "chrom:start-end"
lookup = {}
with open(file2, 'r', encoding='utf-8') as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split()  # any whitespace (tabs or spaces)
        if len(parts) < 7:
            continue  # skip malformed lines
        chrom, start, end = parts[0], parts[1], parts[4]
        tsd, motif = parts[5], parts[6]
        key = f"{chrom}:{start}-{end}"
        lookup[key] = (tsd, motif)

# Stream file1, append TSD & motif when key matches column 1
with open(file1, 'r', encoding='utf-8') as f:
    for raw in f:
        line = raw.rstrip('\n')
        stripped = line.strip()
        if not stripped:
            print(line)
            continue
        key = stripped.split()[0]  # first field of file1 (e.g., Chr2:5627107-5629580)
        if key in lookup:
            tsd, motif = lookup[key]
            print(f"{line}\t{tsd}\t{motif}")
        else:
            print(line)
