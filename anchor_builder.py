#!/usr/bin/env python3
import sys
import argparse
from typing import Iterable, List, Tuple

TRI = '###'

def parse_pairs(lines: Iterable[str]) -> Tuple[List[str], List[str]]:
    """
    Returns (cleaned_lines, transposed_lines) with internal '###' collapsed.
    cleaned_lines:  'A\tB'
    transposed_lines: 'B\tA'
    """
    cleaned, transposed = [], []
    last_was_tri_clean = False
    last_was_tri_trans = False

    for raw in lines:
        s = raw.strip()
        if not s:
            continue

        if s == TRI:
            if cleaned and not last_was_tri_clean:
                cleaned.append(TRI)
            if transposed and not last_was_tri_trans:
                transposed.append(TRI)
            last_was_tri_clean = True if cleaned else False
            last_was_tri_trans = True if transposed else False
            continue

        parts = s.split()
        if len(parts) < 2:
            print(f"Warning: Line '{s}' does not have at least two parts. Skipping.", file=sys.stderr)
            # After a non-### line, reset tri flags
            last_was_tri_clean = False
            last_was_tri_trans = False
            continue

        a, b = parts[0], parts[1]
        cleaned.append(f"{a}\t{b}")
        transposed.append(f"{b}\t{a}")
        last_was_tri_clean = False
        last_was_tri_trans = False

    # Trim trailing '###' inside each list
    if cleaned and cleaned[-1] == TRI:
        cleaned.pop()
    if transposed and transposed[-1] == TRI:
        transposed.pop()

    return cleaned, transposed

def merge_sections(all_output: List[str], section: List[str]) -> None:
    """
    Append `section` to `all_output` while avoiding consecutive '###'.
    """
    if not section:
        return
    if all_output and all_output[-1] == TRI and section[0] == TRI:
        all_output.extend(section[1:])
    else:
        all_output.extend(section)

def main():
    ap = argparse.ArgumentParser(
        description="Build anchor lists with optional transposed pairs and clean '###' boundaries."
    )
    ap.add_argument("files", nargs="+", help="Input text files")
    ap.add_argument("--transpose", action="store_true",
                    help="Append a '###' divider and the transposed pairs (B\\tA) after each file's cleaned section")
    args = ap.parse_args()

    all_output: List[str] = [TRI]  # start with boundary as required

    for path in args.files:
        try:
            with open(path, "r", encoding="utf-8") as f:
                cleaned, transposed = parse_pairs(f)
        except FileNotFoundError:
            print(f"Error: File not found: {path}", file=sys.stderr)
            continue
        except OSError as e:
            print(f"Error: Cannot open {path}: {e}", file=sys.stderr)
            continue

        # Build this file's section (ensure it is wrapped with TRI boundaries)
        section: List[str] = []
        section.append(TRI)

        # cleaned block
        if cleaned:
            section.extend(cleaned)
        # Always end the cleaned block with TRI to separate from next block/file
        section.append(TRI)

        # optional transposed block
        if args.transpose and transposed:
            # Ensure single TRI between cleaned and transposed
            if section and section[-1] != TRI:
                section.append(TRI)
            # Append transposed items and then close with TRI
            section.extend(transposed)
            section.append(TRI)

        # merge into global output without double TRI
        merge_sections(all_output, section)

    # Ensure final ends with exactly one TRI
    if not all_output or all_output[0] != TRI:
        all_output.insert(0, TRI)
    if all_output[-1] != TRI:
        all_output.append(TRI)

    for line in all_output:
        print(line)

if __name__ == "__main__":
    main()
