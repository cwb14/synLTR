#!/usr/bin/env python3
"""
hardmask_genome.py

Hard-mask a genome FASTA based on feature coordinates in a TSV or a features FASTA.

Behavior:
- Coordinates come from either:
  * --tsv: TSV col1 like contig:start-end#anything
  * --features-fasta: FASTA headers like contig:start-end#anything (sequences ignored for coords)
- The feature interval (start-end) is always hard-masked (default: N; configurable via --feature-character)
- Everything farther than --distance bp away from ANY feature is also masked
  (default: same as feature mask; configurable separately via --far-character)
- Only the flanking bases within --distance bp on either side of each feature remain unmasked,
  BUT flanks are truncated early if they would run into another feature interval.

Coordinates:
- Interprets start-end as 1-based, inclusive (common in genomics).
  Example: chr1:100-200 with distance=5 keeps bases 95-99 and 201-205 unmasked
  (clipped to contig bounds). The feature itself (100-200) is always masked.

Masking characters:
- --feature-character masks the feature interval itself
- --far-character masks all other non-kept bases (farther than distance from any feature)

NEW behavior:
- If using --features-fasta, the feature sequence itself is scanned.
  Any NON-ATCG character (case-insensitive) in the feature sequence is assumed to have been
  previously masked and will be preserved EXACTLY at the corresponding genome position
  within the feature interval in the output (e.g., 'X' stays 'X', 'N' stays 'N').
  ATCG positions inside the feature interval are still converted to --feature-character.

Critical rule (fixed):
- Masking ALWAYS takes precedence over keeping flanks.
  Flanks never unmask any base that lies in ANY feature interval.
  Flanks extend up to --distance, but stop early at neighboring feature boundaries.
"""

import argparse
import re
import sys
from typing import Dict, List, Tuple, Iterator, Optional

COORD_RE = re.compile(r"^([^:]+):(\d+)-(\d+)(?:#.*)?$")


def fasta_iter(fp) -> Iterator[Tuple[str, str]]:
    """Yield (header, sequence) from a FASTA file handle."""
    header = None
    seq_chunks: List[str] = []
    for line in fp:
        line = line.rstrip("\n")
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(seq_chunks)
            header = line[1:].split()[0]
            seq_chunks = []
        else:
            seq_chunks.append(line.strip())
    if header is not None:
        yield header, "".join(seq_chunks)


def parse_tsv_features(tsv_path: str) -> Dict[str, List[Tuple[int, int]]]:
    """
    Return dict: contig -> list of feature intervals as 0-based, half-open [start0, end0)
    Parsed from TSV col1: contig:start-end#...
    """
    feats: Dict[str, List[Tuple[int, int]]] = {}
    with open(tsv_path, "r", encoding="utf-8") as f:
        for _ln, line in enumerate(f, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if not parts:
                continue

            col1 = parts[0]
            coord_part = col1.split("#", 1)[0]

            m = COORD_RE.match(coord_part)
            if not m:
                continue

            contig, s_str, e_str = m.group(1), m.group(2), m.group(3)
            s = int(s_str)
            e = int(e_str)
            if e < s:
                s, e = e, s

            # Convert 1-based inclusive [s,e] to 0-based half-open [s-1, e)
            start0 = s - 1
            end0 = e
            feats.setdefault(contig, []).append((start0, end0))

    return feats


def parse_features_fasta_with_overrides(features_fa: str) -> Tuple[
    Dict[str, List[Tuple[int, int]]],
    Dict[str, Dict[int, str]]
]:
    """
    Return:
      feats: contig -> list of feature intervals as 0-based, half-open [start0, end0)
      overrides: contig -> { absolute_0based_pos : char_to_preserve }

    Parsed from FASTA headers like: contig:start-end#...
    The feature sequence is scanned; any NON-ATCG character is preserved at that position
    in the final output (only within the feature interval).
    """
    feats: Dict[str, List[Tuple[int, int]]] = {}
    overrides: Dict[str, Dict[int, str]] = {}

    with open(features_fa, "r", encoding="utf-8") as f:
        for header, seq in fasta_iter(f):
            m = COORD_RE.match(header)
            if not m:
                continue

            contig, s_str, e_str = m.group(1), m.group(2), m.group(3)
            s = int(s_str)
            e = int(e_str)
            if e < s:
                s, e = e, s

            start0 = s - 1
            end0 = e  # half-open
            feats.setdefault(contig, []).append((start0, end0))

            expected_len = end0 - start0  # == e-s+1
            seq = seq.strip()
            if not seq:
                continue

            if len(seq) != expected_len:
                print(
                    f"[warn] Feature seq length != interval length for {header}: "
                    f"seq={len(seq)} expected={expected_len}. Mapping min(len, expected).",
                    file=sys.stderr,
                )

            mlen = min(len(seq), expected_len)
            od = overrides.setdefault(contig, {})

            # Preserve any NON-ATCG (case-insensitive) character as an override
            # at the corresponding genome position within the feature interval.
            for i in range(mlen):
                ch = seq[i]
                up = ch.upper()
                if up not in ("A", "T", "C", "G"):
                    pos0 = start0 + i
                    od[pos0] = ch  # keep exact char (case preserved)

    return feats, overrides


def merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Merge overlapping/adjacent half-open intervals."""
    if not intervals:
        return []
    intervals.sort()
    merged = [intervals[0]]
    for s, e in intervals[1:]:
        ps, pe = merged[-1]
        if s <= pe:  # overlap/adjacent
            merged[-1] = (ps, max(pe, e))
        else:
            merged.append((s, e))
    return merged


def wrap_fasta(seq: str, width: int) -> str:
    if width <= 0:
        return seq + "\n"
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width)) + "\n"


def clip_intervals(intervals: List[Tuple[int, int]], n: int) -> List[Tuple[int, int]]:
    """Clip half-open intervals to [0, n) and drop empties."""
    out: List[Tuple[int, int]] = []
    for s, e in intervals:
        s = max(0, min(n, s))
        e = max(0, min(n, e))
        if e > s:
            out.append((s, e))
    return out


def build_keep_flanks_truncated(
    feats_merged: List[Tuple[int, int]],
    n: int,
    dist: int
) -> List[Tuple[int, int]]:
    """
    Build keep-intervals = left/right flanks around each feature, truncated so they
    never enter any feature interval.

    Given merged, sorted feature intervals feats_merged:
      left flank of feature i: [max(prev_end, s - dist), s)
      right flank of feature i: [e, min(next_start, e + dist))

    This ensures:
    - feature masking always wins
    - flanks stop early at neighboring features
    """
    if dist <= 0 or not feats_merged:
        return []

    keep: List[Tuple[int, int]] = []
    k = len(feats_merged)

    for idx, (s, e) in enumerate(feats_merged):
        prev_end = feats_merged[idx - 1][1] if idx > 0 else 0
        next_start = feats_merged[idx + 1][0] if idx + 1 < k else n

        # left flank, truncated by previous feature end
        ls = max(0, s - dist, prev_end)
        le = s
        if le > ls:
            keep.append((ls, le))

        # right flank, truncated by next feature start
        rs = e
        re = min(n, e + dist, next_start)
        if re > rs:
            keep.append((rs, re))

    return merge_intervals(keep)


def mask_contig(
    seq: str,
    features: List[Tuple[int, int]],
    feature_char: str,
    far_char: str,
    dist: int,
    override_map: Optional[Dict[int, str]] = None,
) -> str:
    """
    Mask with two characters:
      - feature interval itself: feature_char (or override char if specified for that position)
      - keep flanks within dist of features: keep original bases, but only outside ALL features
      - everything else (farther than dist from any feature): far_char

    Priority is ALWAYS: feature > keep > far
    """
    n = len(seq)
    if n == 0:
        return seq

    # Clip + merge feature intervals
    feats = merge_intervals(clip_intervals(features, n))

    # If no features: everything is "far"
    if not feats:
        return far_char * n

    # Keep intervals are truncated flanks that never enter any feature interval
    keep = build_keep_flanks_truncated(feats, n, dist)

    out_chunks: List[str] = []

    i = 0
    kp = 0  # keep pointer
    fp = 0  # feature pointer

    # Prepare override scanning
    if override_map:
        ov_items = [(p, c) for p, c in override_map.items() if 0 <= p < n]
        ov_items.sort(key=lambda x: x[0])
        ov_pos = [p for p, _c in ov_items]
        ov_chr = [c for _p, c in ov_items]
    else:
        ov_pos = []
        ov_chr = []
    op = 0  # override pointer

    def in_interval(intervals: List[Tuple[int, int]], p: int, pos: int) -> bool:
        return p < len(intervals) and intervals[p][0] <= pos < intervals[p][1]

    while i < n:
        # Advance pointers if current position is beyond interval end
        while kp < len(keep) and i >= keep[kp][1]:
            kp += 1
        while fp < len(feats) and i >= feats[fp][1]:
            fp += 1
        while op < len(ov_pos) and ov_pos[op] < i:
            op += 1

        in_feat = in_interval(feats, fp, i)
        in_keep = in_interval(keep, kp, i)

        # Determine the next boundary where state could change
        next_b = n
        if fp < len(feats):
            fs, fe = feats[fp]
            if i < fs:
                next_b = min(next_b, fs)
            else:
                next_b = min(next_b, fe)
        if kp < len(keep):
            ks, ke = keep[kp]
            if i < ks:
                next_b = min(next_b, ks)
            else:
                next_b = min(next_b, ke)

        # Emit chunk with priority: FEATURE > KEEP > FAR
        if in_feat:
            # Feature interval: mostly feature_char, but preserve prior non-ATCG mask chars if provided
            if not ov_pos:
                out_chunks.append(feature_char * (next_b - i))
            else:
                start = i
                while op < len(ov_pos) and ov_pos[op] < next_b:
                    p = ov_pos[op]
                    if p > start:
                        out_chunks.append(feature_char * (p - start))
                    out_chunks.append(ov_chr[op])
                    start = p + 1
                    op += 1
                if start < next_b:
                    out_chunks.append(feature_char * (next_b - start))

        elif in_keep:
            out_chunks.append(seq[i:next_b])

        else:
            out_chunks.append(far_char * (next_b - i))

        i = next_b

    return "".join(out_chunks)


def main():
    ap = argparse.ArgumentParser(
        description="Hardmask genome based on TSV/FASTA coordinates; keep only flanks within distance."
    )
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--tsv", help="TSV file with col1 like contig:start-end#...")
    src.add_argument(
        "--features-fasta",
        help="Multi-FASTA where each header is contig:start-end#... (sequences used to preserve prior non-ATCG masks)",
    )
    ap.add_argument("--genome", required=True, help="Genome FASTA")
    ap.add_argument(
        "--feature-character",
        default="N",
        help="Mask character for feature intervals themselves (default: N)",
    )
    ap.add_argument(
        "--far-character",
        default=None,
        help="Mask character for bases farther than distance from any feature "
             "(default: same as --feature-character)",
    )
    ap.add_argument("--distance", type=int, default=0, help="Flank distance to keep unmasked (bp)")
    ap.add_argument("--wrap", type=int, default=60, help="FASTA line wrap width (default: 60; 0 = no wrap)")
    args = ap.parse_args()

    if len(args.feature_character) != 1:
        ap.error("--feature-character must be a single character (e.g., N)")
    if args.far_character is None:
        args.far_character = args.feature_character
    if len(args.far_character) != 1:
        ap.error("--far-character must be a single character (e.g., X)")
    if args.distance < 0:
        ap.error("--distance must be >= 0")

    overrides_by_contig: Dict[str, Dict[int, str]] = {}

    if args.tsv:
        feats_by_contig = parse_tsv_features(args.tsv)
    else:
        feats_by_contig, overrides_by_contig = parse_features_fasta_with_overrides(args.features_fasta)

    with open(args.genome, "r", encoding="utf-8") as gf:
        for header, seq in fasta_iter(gf):
            feats = feats_by_contig.get(header, [])
            ov = overrides_by_contig.get(header, None)
            masked = mask_contig(
                seq=seq,
                features=feats,
                feature_char=args.feature_character,
                far_char=args.far_character,
                dist=args.distance,
                override_map=ov,
            )
            sys.stdout.write(f">{header}\n")
            sys.stdout.write(wrap_fasta(masked, args.wrap))


if __name__ == "__main__":
    main()
