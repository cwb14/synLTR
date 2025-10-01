#!/usr/bin/env python3
import sys
import argparse

def process_file(infile: str, threshold: int, stitch_gaps: bool):
    """
    Process the anchors coordinate file to merge lines based on overlapping
    or touching (butt heads) sequence ranges, per (pair1_id, pair2_id, strand) bin.

    A pair is considered "pass" if both lengths >= threshold; otherwise "fail".
    Only lines where at least one side is "fail" are eligible to merge.
    After merging, the merged record's status becomes "pass".

    If --stitch-gaps is set, after merging we insert synthetic lines to fill
    gaps between consecutive records within each bin (same pair1_id, pair2_id, strand)
    whenever both sequences have a positive gap. Touching intervals (no gap) are not stitched.
    """
    # ---------- parsing ----------
    bins = {}

    def parse_range(token: str):
        # e.g., "9311v2_chr1:18291..25404" -> ("9311v2_chr1", 18291, 25404)
        id_part, coords = token.split(':', 1)
        start_str, end_str = coords.split('..', 1)
        return id_part, int(start_str), int(end_str)

    with open(infile, 'r') as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                continue  # Skip malformed lines quietly

            pair1_id, pair1_start, pair1_end = parse_range(parts[0])
            pair2_id, pair2_start, pair2_end = parse_range(parts[1])
            strand = parts[2].strip()

            # Calculate lengths (half-open style; matches original behavior)
            len1 = pair1_end - pair1_start
            len2 = pair2_end - pair2_start

            threshold_status = "pass" if (len1 >= threshold and len2 >= threshold) else "fail"
            bin_key = f"{pair1_id}_{pair2_id}_{'plus' if strand == '+' else 'minus'}"

            new_line = [
                pair1_id, pair1_start, pair1_end,   # 0..2
                pair2_id, pair2_start, pair2_end,   # 3..5
                strand,                             # 6
                len1, len2,                         # 7..8
                bin_key,                            # 9
                threshold_status                    # 10
            ]
            bins.setdefault(bin_key, []).append(new_line)

    # ---------- merging ----------
    def ranges_touch_or_overlap(a_start, a_end, b_start, b_end):
        # Allow touching: end == start
        return (a_start <= b_end) and (a_end >= b_start)

    final_bins = {}

    for bin_key, lines in bins.items():
        # Sort to make merging deterministic and efficient
        lines.sort(key=lambda x: (x[1], x[4]))  # by pair1_start then pair2_start
        merged_list = []

        while lines:
            current_line = lines.pop(0)
            merged = False

            for i, line in enumerate(lines):
                cond1 = ranges_touch_or_overlap(current_line[1], current_line[2], line[1], line[2])
                cond2 = ranges_touch_or_overlap(current_line[4], current_line[5], line[4], line[5])

                # Merge if both sequences overlap or touch AND at least one record is "fail"
                if cond1 and cond2 and (current_line[10] == "fail" or line[10] == "fail"):
                    merged_line = [
                        current_line[0],
                        min(current_line[1], line[1]),
                        max(current_line[2], line[2]),
                        current_line[3],
                        min(current_line[4], line[4]),
                        max(current_line[5], line[5]),
                        current_line[6],
                        # these len fields aren’t used for output, but keep them coherent
                        max(current_line[2], line[2]) - min(current_line[1], line[1]),
                        max(current_line[5], line[5]) - min(current_line[4], line[4]),
                        current_line[9],
                        "pass"  # merged becomes pass
                    ]
                    lines[i] = merged_line
                    merged = True
                    break

            if not merged:
                merged_list.append(current_line)

        # Keep per-bin results for optional stitching
        final_bins[bin_key] = merged_list

    # ---------- optional stitching ----------
    if stitch_gaps:
        for bin_key, lines in final_bins.items():
            if not lines:
                continue
            # Sort again to ensure chronological order for stitching
            lines.sort(key=lambda x: (x[1], x[4]))

            stitched = []
            for prev, nxt in zip(lines, lines[1:]):
                stitched.append(prev)
                # Compute gaps on both sequences (positive gap only)
                gap1_start, gap1_end = prev[2], nxt[1]
                gap2_start, gap2_end = prev[5], nxt[4]

                if (gap1_end > gap1_start) and (gap2_end > gap2_start):
                    # Insert a synthetic line covering the gap on both sequences
                    stitched.append([
                        prev[0], gap1_start, gap1_end,
                        prev[3], gap2_start, gap2_end,
                        prev[6],
                        gap1_end - gap1_start,
                        gap2_end - gap2_start,
                        bin_key,
                        "stitched"  # internal marker; output format ignores this
                    ])
            # Append the last original record
            stitched.append(lines[-1])
            final_bins[bin_key] = stitched

    # ---------- output ----------
    # Flatten bins in insertion order; within bin keep sorted order for determinism
    for bin_key in final_bins:
        for line in sorted(final_bins[bin_key], key=lambda x: (x[1], x[4])):
            print(f"{line[0]}:{line[1]}..{line[2]}\t{line[3]}:{line[4]}..{line[5]}\t{line[6]}")

def main():
    parser = argparse.ArgumentParser(
        description="Merge anchor coordinate lines by overlapping/touching ranges within (pair1_id, pair2_id, strand) bins, with optional gap stitching."
    )
    parser.add_argument(
        "infile",
        help="Path to the input anchors coordinate file (e.g., all.recip.anchors.coords)"
    )
    parser.add_argument(
        "-t", "--threshold",
        type=int,
        default=1_000_000,
        help="Length threshold for each side to be considered 'pass' (default: 1000000)"
    )
    parser.add_argument(
        "--stitch-gaps",
        action="store_true",
        help="After merging, insert synthetic lines to fill positive gaps between consecutive records within each (pair1, pair2, strand) bin."
    )
    args = parser.parse_args()

    if args.threshold < 0:
        print("Threshold must be non-negative.", file=sys.stderr)
        sys.exit(2)

    process_file(args.infile, args.threshold, args.stitch_gaps)

if __name__ == "__main__":
    main()

