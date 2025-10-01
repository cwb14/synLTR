#!/usr/bin/env python3
import sys
import argparse

def process_file(infile: str, threshold: int, stitch_gaps: bool):
    """
    Process the anchors coordinate file to merge lines based on overlapping
    or touching (butt heads) sequence ranges, per (pair1_id, pair2_id, strand) bin.

    A pair is considered "pass" if both lengths >= threshold; otherwise "fail".
    Only lines where at least one side is "fail" are eligible to merge.
    After merging, the merged record's status becomes "pass" (as in original script).

    If stitch_gaps is True, post-process each bin to insert synthetic lines that
    fill gaps between consecutive anchors where BOTH pair1 and pair2 have a positive gap.
    """
    bins = {}
    merged_lines = []

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
                # Skip malformed lines quietly
                continue

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

    # Helper: overlap or touch? (butt heads allowed)
    def ranges_touch_or_overlap(a_start, a_end, b_start, b_end):
        # Touching means end == start across intervals
        return (a_start <= b_end) and (a_end >= b_start)

    # Merge within each bin
    for bin_key, lines in bins.items():
        # Sort to make merging deterministic
        lines.sort(key=lambda x: (x[1], x[4]))  # sort by pair1_start then pair2_start
        while lines:
            current_line = lines.pop(0)
            merged = False
            for i, line in enumerate(lines):
                cond1 = ranges_touch_or_overlap(current_line[1], current_line[2], line[1], line[2])
                cond2 = ranges_touch_or_overlap(current_line[4], current_line[5], line[4], line[5])

                # Eligible if intervals overlap OR touch on BOTH sequences, and
                # at least one of the two records is "fail"
                if cond1 and cond2 and (current_line[10] == "fail" or line[10] == "fail"):
                    merged_line = [
                        current_line[0],
                        min(current_line[1], line[1]),
                        max(current_line[2], line[2]),
                        current_line[3],
                        min(current_line[4], line[4]),
                        max(current_line[5], line[5]),
                        current_line[6],
                        min(current_line[7], line[7]),
                        max(current_line[8], line[8]),
                        current_line[9],
                        "pass"  # keep behavior: merged becomes "pass"
                    ]
                    lines[i] = merged_line
                    merged = True
                    break

            if not merged:
                merged_lines.append(current_line)

    # Optional: stitch gaps inside each bin
    if stitch_gaps:
        # regroup merged_lines by bin for ordered stitching
        groups = {}
        for rec in merged_lines:
            groups.setdefault(rec[9], []).append(rec)

        stitched_output = []
        # Sort bins to keep output stable (optional)
        for bin_key in sorted(groups.keys()):
            grp = groups[bin_key]
            # Sort records in ascending coordinate order (works for + and - examples)
            grp.sort(key=lambda x: (x[1], x[4]))
            if not grp:
                continue

            # Walk consecutive pairs; insert a gap record when BOTH sides have a positive gap
            stitched_grp = []
            for i in range(len(grp) - 1):
                cur = grp[i]
                nxt = grp[i + 1]
                stitched_grp.append(cur)

                gap1_start, gap1_end = cur[2], nxt[1]
                gap2_start, gap2_end = cur[5], nxt[4]
                g1 = gap1_end - gap1_start
                g2 = gap2_end - gap2_start

                # Insert only if a real gap on BOTH sequences (>0). Butt heads (==0) are ignored.
                if g1 > 0 and g2 > 0:
                    gap_rec = [
                        cur[0], gap1_start, gap1_end,   # pair1 gap
                        cur[3], gap2_start, gap2_end,   # pair2 gap
                        cur[6],                         # strand
                        g1, g2,                         # lengths for completeness
                        cur[9],                         # same bin_key
                        "gap"                           # status label (not used in output)
                    ]
                    stitched_grp.append(gap_rec)

            # add final record
            stitched_grp.append(grp[-1])
            stitched_output.extend(stitched_grp)

        merged_lines = stitched_output

    # Output in original 3-column format
    for line in merged_lines:
        print(f"{line[0]}:{line[1]}..{line[2]}\t{line[3]}:{line[4]}..{line[5]}\t{line[6]}")

def main():
    parser = argparse.ArgumentParser(
        description="Merge anchor coordinate lines by overlapping or touching ranges within (pair1_id, pair2_id, strand) bins. Optionally, stitch gaps between consecutive anchors."
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
        help="After consolidation, insert synthetic lines that fill positive gaps between consecutive anchors within each (pair1_id, pair2_id, strand) bin."
    )
    args = parser.parse_args()

    if args.threshold < 0:
        print("Threshold must be non-negative.", file=sys.stderr)
        sys.exit(2)

    process_file(args.infile, args.threshold, args.stitch_gaps)

if __name__ == "__main__":
    main()
