#!/usr/bin/env python3
import sys
import argparse

def load_fai_lengths(path):
    """
    Read a .fai (samtools faidx) and return {seq_name: length}.
    """
    lengths = {}
    with open(path, 'r') as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            # .fai: name, length, offset, line_bases, line_width
            if len(parts) >= 2:
                lengths[parts[0]] = int(parts[1])
    return lengths

def process_file(infile: str, threshold: int, stitch_gaps: bool,
                 stitch_terminals: bool, fai1_path: str, fai2_path: str, coord_base: int):
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

    If stitch_terminals is True, also stitch leading/trailing gaps to chromosome bounds
    using lengths from .fai files. Terminal stitching is applied **only to '+' bins**
    to keep coordinates parallel to 0/1-based starts as per your example.
    """
    if stitch_terminals:
        if not (fai1_path and fai2_path):
            print("Error: --stitch-terminals requires --fai1 and --fai2.", file=sys.stderr)
            sys.exit(2)
        len1 = load_fai_lengths(fai1_path)
        len2 = load_fai_lengths(fai2_path)
    else:
        len1 = {}
        len2 = {}

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
                continue

            p1_id, p1_s, p1_e = parse_range(parts[0])
            p2_id, p2_s, p2_e = parse_range(parts[1])
            strand = parts[2].strip()

            # Length as half-open: end - start
            l1 = p1_e - p1_s
            l2 = p2_e - p2_s
            status = "pass" if (l1 >= threshold and l2 >= threshold) else "fail"

            key = f"{p1_id}_{p2_id}_{'plus' if strand == '+' else 'minus'}"
            bins.setdefault(key, []).append([
                p1_id, p1_s, p1_e,
                p2_id, p2_s, p2_e,
                strand,
                l1, l2,
                key,
                status
            ])

    def ranges_touch_or_overlap(a_start, a_end, b_start, b_end):
        return (a_start <= b_end) and (a_end >= b_start)

    merged_by_bin = {}
    all_merged = []

    # First pass: merge by overlap/touch with "fail in either" rule
    for key, lines in bins.items():
        lines.sort(key=lambda x: (x[1], x[4]))
        out = []
        work = lines[:]
        while work:
            cur = work.pop(0)
            merged = False
            for i, nxt in enumerate(work):
                cond1 = ranges_touch_or_overlap(cur[1], cur[2], nxt[1], nxt[2])
                cond2 = ranges_touch_or_overlap(cur[4], cur[5], nxt[4], nxt[5])
                if cond1 and cond2 and (cur[10] == "fail" or nxt[10] == "fail"):
                    m = [
                        cur[0], min(cur[1], nxt[1]), max(cur[2], nxt[2]),
                        cur[3], min(cur[4], nxt[4]), max(cur[5], nxt[5]),
                        cur[6],
                        min(cur[7], nxt[7]), max(cur[8], nxt[8]),
                        cur[9],
                        "pass"
                    ]
                    work[i] = m
                    merged = True
                    break
            if not merged:
                out.append(cur)
        merged_by_bin[key] = out
        all_merged.extend(out)

    # Optional: stitch internal gaps
    def stitch_internal(bin_key, blocks):
        if not stitch_gaps or not blocks:
            return blocks[:]
        strand = blocks[0][6]
        is_plus = (strand == '+')
        blocks_sorted = sorted(blocks, key=lambda x: x[1])

        stitched = []
        for i in range(len(blocks_sorted) - 1):
            cur = blocks_sorted[i]
            nxt = blocks_sorted[i + 1]
            stitched.append(cur)
            gap1 = (cur[2] < nxt[1])  # seq1 forward gap
            if is_plus:
                gap2 = (cur[5] < nxt[4])  # seq2 forward gap
                if gap1 and gap2:
                    stitched.append([
                        cur[0], cur[2], nxt[1],
                        cur[3], cur[5], nxt[4],
                        strand,
                        nxt[1] - cur[2],
                        nxt[4] - cur[5],
                        bin_key,
                        "pass"
                    ])
            else:
                gap2 = (nxt[5] < cur[4])  # seq2 reverse-direction gap
                if gap1 and gap2:
                    stitched.append([
                        cur[0], cur[2], nxt[1],
                        cur[3], nxt[5], cur[4],
                        strand,
                        nxt[1] - cur[2],
                        cur[4] - nxt[5],
                        bin_key,
                        "pass"
                    ])
        stitched.append(blocks_sorted[-1])
        return stitched

    # Optional: stitch terminal gaps to chromosome bounds (PLUS only)
    def stitch_terminals_plus(bin_key, blocks):
        if not stitch_terminals or not blocks:
            return blocks
        p1 = blocks[0][0]
        p2 = blocks[0][3]

        # Require lengths for both chroms
        if p1 not in len1 or p2 not in len2:
            return blocks

        base = coord_base  # 0 or 1
        end1 = len1[p1]
        end2 = len2[p2]

        # Sort by seq1 start (increasing)
        b = sorted(blocks, key=lambda x: x[1])
        out = []

        # Leading gap: [base .. first_start] on both chromosomes
        first = b[0]
        if first[1] > base and first[4] > base:
            out.append([
                p1, base, first[1],
                p2, base, first[4],
                '+',
                first[1] - base,
                first[4] - base,
                bin_key,
                "pass"
            ])

        out.extend(b)

        # Trailing gap: [last_end .. chrom_end] on both chromosomes
        last = b[-1]
        if last[2] < end1 and last[5] < end2:
            out.append([
                p1, last[2], end1,
                p2, last[5], end2,
                '+',
                end1 - last[2],
                end2 - last[5],
                bin_key,
                "pass"
            ])
        return out

    final = []
    for key, blocks in merged_by_bin.items():
        # stitch internal gaps first
        blocks2 = stitch_internal(key, blocks)

        # try terminal stitch (PLUS bins only)
        if blocks2 and blocks2[0][6] == '+':
            blocks2 = stitch_terminals_plus(key, blocks2)

        final.extend(sorted(blocks2, key=lambda x: (x[0], x[3], x[6], x[1], x[4])))

    # Emit
    for line in final:
        print(f"{line[0]}:{line[1]}..{line[2]}\t{line[3]}:{line[4]}..{line[5]}\t{line[6]}")

def main():
    ap = argparse.ArgumentParser(
        description="Merge anchors by overlap/touch within (pair1_id, pair2_id, strand) bins; "
                    "optionally stitch internal gaps and chromosome terminal gaps."
    )
    ap.add_argument("infile", help="Path to the input anchors coordinate file (e.g., all.recip.anchors.coords)")
    ap.add_argument("-t", "--threshold", type=int, default=1_000_000,
                    help="Length threshold for each side to be considered 'pass' (default: 1000000)")
    ap.add_argument("--stitch-gaps", action="store_true",
                    help="Insert synthetic lines to fill strict internal gaps between adjacent blocks within each bin.")
    ap.add_argument("--stitch-terminals", action="store_true",
                    help="Also stitch leading/trailing gaps to chromosome bounds (PLUS strand bins only). Requires --fai1 and --fai2.")
    ap.add_argument("--fai1", help="FAI index for genome 1 (for chromosome lengths). Required if --stitch-terminals.")
    ap.add_argument("--fai2", help="FAI index for genome 2 (for chromosome lengths). Required if --stitch-terminals.")
    ap.add_argument("--coord-base", type=int, choices=[0, 1], default=0,
                    help="Coordinate base for starts (0 or 1). Default: 0")

    args = ap.parse_args()
    if args.threshold < 0:
        print("Threshold must be non-negative.", file=sys.stderr)
        sys.exit(2)
    if args.stitch_terminals and (not args.fai1 or not args.fai2):
        print("Error: --stitch-terminals requires --fai1 and --fai2.", file=sys.stderr)
        sys.exit(2)

    process_file(args.infile, args.threshold, args.stitch_gaps,
                 args.stitch_terminals, args.fai1, args.fai2, args.coord_base)

if __name__ == "__main__":
    main()
