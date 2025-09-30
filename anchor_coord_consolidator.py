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

    If stitch_gaps is True, perform a post-pass that inserts synthetic lines
    that fill strict gaps between adjacent merged blocks within each bin:
      - '+' strand: fill when both sequences have a forward gap.
      - '-' strand: fill when seq1 has a forward gap and seq2 has a reverse-direction gap.
    """
    bins = {}
    merged_lines_global = []  # fallback if we don't stitch

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
    # Touching means end == start across intervals; we use <= / >= to allow that.
    def ranges_touch_or_overlap(a_start, a_end, b_start, b_end):
        return (a_start <= b_end) and (a_end >= b_start)

    # Merge within each bin
    merged_by_bin = {}
    for bin_key, lines in bins.items():
        # Sort to make merging deterministic and slightly more efficient
        lines.sort(key=lambda x: (x[1], x[4]))  # sort by pair1_start then pair2_start
        merged_lines = []
        work = lines[:]  # copy

        while work:
            current_line = work.pop(0)
            merged = False

            for i, line in enumerate(work):
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
                    work[i] = merged_line
                    merged = True
                    break

            if not merged:
                merged_lines.append(current_line)

        merged_by_bin[bin_key] = merged_lines
        merged_lines_global.extend(merged_lines)

    # Optional: stitch gaps inside each bin
    if stitch_gaps:
        stitched_output = []
        for bin_key, lines in merged_by_bin.items():
            if not lines:
                continue
            # Sort final merged blocks for deterministic adjacency
            # (By seq1 start; for '-' this still makes seq1 increasing.)
            lines_sorted = sorted(lines, key=lambda x: x[1])

            strand = lines_sorted[0][6]
            is_plus = (strand == '+')

            i = 0
            while i < len(lines_sorted) - 1:
                cur = lines_sorted[i]
                nxt = lines_sorted[i + 1]
                stitched_output.append(cur)

                # seq1 gap (forward direction)
                gap1 = (cur[2] < nxt[1])  # cur.end1 < next.start1 (strict gap)

                if is_plus:
                    # seq2 gap must also be forward
                    gap2 = (cur[5] < nxt[4])  # cur.end2 < next.start2
                    if gap1 and gap2:
                        # Insert a synthetic line spanning the gaps
                        stitched_output.append([
                            cur[0], cur[2], nxt[1],   # seq1: [end1 .. next.start1]
                            cur[3], cur[5], nxt[4],   # seq2: [end2 .. next.start2]
                            strand,
                            nxt[1] - cur[2],          # len1
                            nxt[4] - cur[5],          # len2
                            bin_key,
                            "pass"                    # status doesn't affect output format
                        ])
                else:
                    # '-' strand: seq2 runs in reverse across blocks;
                    # adjacency is next.end2 .. current.start2 (increasing coords)
                    gap2 = (nxt[5] < cur[4])  # next.end2 < cur.start2
                    if gap1 and gap2:
                        stitched_output.append([
                            cur[0], cur[2], nxt[1],   # seq1 forward gap
                            cur[3], nxt[5], cur[4],   # seq2 reverse-direction gap as increasing interval
                            strand,
                            nxt[1] - cur[2],
                            cur[4] - nxt[5],
                            bin_key,
                            "pass"
                        ])
                i += 1

            # push the last block
            stitched_output.append(lines_sorted[-1])

        # Emit stitched result
        for line in stitched_output:
            print(f"{line[0]}:{line[1]}..{line[2]}\t{line[3]}:{line[4]}..{line[5]}\t{line[6]}")
    else:
        # Emit un-stitched merged result
        for line in merged_lines_global:
            print(f"{line[0]}:{line[1]}..{line[2]}\t{line[3]}:{line[4]}..{line[5]}\t{line[6]}")

def main():
    parser = argparse.ArgumentParser(
        description="Merge anchor coordinate lines by overlapping/touching ranges within (pair1_id, pair2_id, strand) bins; optionally stitch inter-block gaps."
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
        help="After merging, insert synthetic lines to fill strict gaps between adjacent blocks within each bin."
    )
    args = parser.parse_args()

    if args.threshold < 0:
        print("Threshold must be non-negative.", file=sys.stderr)
        sys.exit(2)

    process_file(args.infile, args.threshold, args.stitch_gaps)

if __name__ == "__main__":
    main()
