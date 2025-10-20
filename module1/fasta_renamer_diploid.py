#!/usr/bin/env python3
import re
import os
import subprocess
from argparse import ArgumentParser
from itertools import combinations
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import defaultdict, Counter


# -----------------------------
# Helpers for chr parsing/naming
# -----------------------------

CHR_PREFIX_RE = re.compile(
    r'(Chr|chr|Chro|chro|Chrom|chrom|Chromosome|chromosome|CHROMOSOME|CHR|CHRO|CHROM)'
    r'[ _-]*:? ?'                       # optional separators/colon
    r'(\d+|[XYZW])'                     # number or X/Y/Z/W
    r'([A-Za-z]+)?'                     # optional existing letter(s), e.g. A, B, AA
)

def index_to_letters(n: int) -> str:
    """
    Convert 1 -> 'A', 2 -> 'B', ... 26 -> 'Z', 27 -> 'AA', etc. (Excel-style)
    Assumes n >= 1.
    """
    letters = []
    while n > 0:
        n, r = divmod(n - 1, 26)
        letters.append(chr(65 + r))
    return ''.join(reversed(letters))


def parse_chr_from_header(header: str):
    """
    Parse 'chr-like' label from a FASTA header line (including '>').
    Returns (chrom_number_str, letter_or_None) or (None, None) if no chr match.
    - Removes leading zeros from numeric chromosomes.
    - Normalizes letter suffix to uppercase (if present).
    """
    m = CHR_PREFIX_RE.search(header)
    if not m:
        return None, None
    chrom = m.group(2)
    if chrom.isdigit():
        chrom = str(int(chrom))  # normalize leading zeros
    letter = m.group(3).upper() if m.group(3) else None
    return chrom, letter


def make_chr_header(prefix: str, chrom: str, letter: str | None) -> str:
    """
    Construct the new header (with '>') using prefix, chrom number, and optional letter.
    """
    if letter:
        return f'>{prefix}_chr{chrom}{letter}'
    return f'>{prefix}_chr{chrom}'


# -----------------------------
# Core processing
# -----------------------------

def process_file(file_path, pass_files, out_dir, out_suffix):
    file_prefix = os.path.splitext(os.path.basename(file_path))[0]
    if not os.path.isfile(file_path):
        print(f"[ERROR] Genome file not found: {file_path}")
        return

    with open(file_path, 'r') as fasta_file:
        lines = fasta_file.readlines()

    # Parse FASTA
    sequences = {}
    current_header = None
    for line in lines:
        if line.startswith('>'):
            current_header = line.strip()
            sequences[current_header] = []
        elif current_header:
            sequences[current_header].append(line.strip())

    sequence_lengths = {header: sum(len(seq) for seq in seqs) for header, seqs in sequences.items()}

    # Partition into chr-like and non-chr (sca) sequences using the new parser
    chr_sequences = {}
    sca_sequences = {}
    for header, seqs in sequences.items():
        chrom, letter = parse_chr_from_header(header)
        if chrom is not None:
            chr_sequences[header] = (seqs, chrom, letter)
        else:
            sca_sequences[header] = seqs

    # Sort sca sequences by size (desc) for stable fallback numbering
    sorted_sca_sequences = sorted(sca_sequences.items(), key=lambda item: -sequence_lengths[item[0]])

    new_lines = []
    used_headers = set()
    header_mapping = []

    # -----------------------------
    # Handle chr-groups (by chrom number)
    #   - Preserve existing letters where present.
    #   - If duplicates exist (group size > 1 OR collisions), assign letters
    #     starting from 'A' to unlettered entries and to any colliding extras.
    #   - If a chrom number is unique and has no letter, keep it unlettered.
    # -----------------------------
    # Group entries: chrom_number -> list of items
    # item: (old_header, seqs, chrom, letter)
    groups = defaultdict(list)
    for old_header, (seqs, chrom, letter) in chr_sequences.items():
        groups[chrom].append((old_header, seqs, chrom, letter))

    for chrom, items in groups.items():
        # Track letters present; handle collisions where the same letter appears multiple times
        # Preserve the first occurrence of an existing letter; extras will be reassigned.
        letter_counts = Counter([it[3] for it in items if it[3] is not None])
        kept_letters = set()
        collisions = []   # items that need new letter due to collision
        unlettered = []   # items with no letter originally

        # Assign buckets while preserving order of appearance
        prepared = []  # (old_header, seqs, chrom, letter_or_None or 'COLLISION')
        seen_for_letter = Counter()
        for old_header, seqs, chromN, letter in items:
            if letter is None:
                unlettered.append((old_header, seqs, chromN))
                prepared.append((old_header, seqs, chromN, None))
            else:
                seen_for_letter[letter] += 1
                if seen_for_letter[letter] == 1:
                    kept_letters.add(letter)
                    prepared.append((old_header, seqs, chromN, letter))
                else:
                    # collision on same letter; will need reassignment
                    collisions.append((old_header, seqs, chromN))
                    prepared.append((old_header, seqs, chromN, 'COLLISION'))

        # Determine if we need to assign letters in this group
        need_letters = (
            len(items) > 1  # multiple entries for this chrom number
            or any(p[3] == 'COLLISION' for p in prepared)
        )

        # Build an iterator of next available letters starting at 'A', skipping kept ones
        def next_letter_generator():
            idx = 1  # A
            while True:
                cand = index_to_letters(idx)
                if cand not in kept_letters:
                    yield cand
                idx += 1

        letter_iter = next_letter_generator()

        # Assign letters to unlettered (only if needed) and to collisions
        assigned_letters = {}
        if need_letters:
            # First handle collisions so their letters change away from duplicates
            for old_header, seqs, chromN in collisions:
                ltr = next(letter_iter)
                kept_letters.add(ltr)
                assigned_letters[old_header] = ltr
            # Then handle truly unlettered ones
            for old_header, seqs, chromN in unlettered:
                ltr = next(letter_iter)
                kept_letters.add(ltr)
                assigned_letters[old_header] = ltr

        # Emit sequences in original group order with their final headers
        for old_header, seqs, chromN, letter in prepared:
            if letter == 'COLLISION':
                new_letter = assigned_letters[old_header]
            elif letter is None:
                # If not needed, leave unlettered (unique number case).
                # If needed, retrieve assigned letter.
                if need_letters:
                    new_letter = assigned_letters[old_header]
                else:
                    new_letter = None
            else:
                new_letter = letter  # preserved

            new_header_with_gt = make_chr_header(file_prefix, chromN, new_letter)
            # Safety: ensure global uniqueness (extremely unlikely after our logic, but guard anyway)
            base = new_header_with_gt
            if base in used_headers:
                # rare: global collision across chroms; bump letters
                bump_iter = next_letter_generator()
                while base in used_headers:
                    bump = next(bump_iter)
                    base = make_chr_header(file_prefix, chromN, bump)
                new_header_with_gt = base

            new_lines.append(f'{new_header_with_gt}\n')
            used_headers.add(new_header_with_gt)
            old_header_no_gt = old_header[1:]
            new_header_no_gt = new_header_with_gt[1:]
            header_mapping.append(f'{old_header_no_gt}\t{new_header_no_gt}\n')
            new_lines.append(''.join(seqs) + '\n')

    # -----------------------------
    # Handle non-chr (sca) sequences
    # -----------------------------
    fallback_counter = 1
    for header, seqs in sorted_sca_sequences:
        # ensure unique sca header
        while f'>{file_prefix}_sca{fallback_counter}' in used_headers:
            fallback_counter += 1
        new_header_full = f'>{file_prefix}_sca{fallback_counter}'
        # Historical constraint: exclude overly long 'sca' headers; chr headers are exempt
        if len(new_header_full) > 14:
            print(f"Excluding sequence with header {new_header_full} due to length > 13 characters.")
            break
        old_header_no_gt = header.strip()[1:]
        new_lines.append(f'{new_header_full}\n')
        used_headers.add(new_header_full)
        fallback_counter += 1
        header_mapping.append(f'{old_header_no_gt}\t{new_header_full[1:]}\n')
        new_lines.append(''.join(seqs) + '\n')

    # -----------------------------
    # Pipe through bioawk for clean FASTA formatting
    # -----------------------------
    bioawk_command = ["bioawk", "-c", "fastx", '{print ">"$name; print $seq}']
    try:
        process = subprocess.Popen(
            bioawk_command,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=False
        )
    except FileNotFoundError:
        raise RuntimeError("bioawk not found on PATH. Please install bioawk and ensure it is accessible.")

    fasta_content = ''.join(new_lines).encode('utf-8')
    bioawk_output, bioawk_err = process.communicate(fasta_content)
    if process.returncode != 0:
        err_text = bioawk_err.decode('utf-8', errors='replace') if bioawk_err else ''
        raise RuntimeError(f"bioawk failed for {file_path} with exit code {process.returncode}.\n{err_text}")

    # -----------------------------
    # Write outputs
    # -----------------------------
    suffix = "" if out_suffix == "disable" else out_suffix
    fasta_filename = f"{file_prefix}{suffix}.fa"
    output_file_path = os.path.join(out_dir if out_dir else ".", fasta_filename)

    with open(output_file_path, 'wb') as output_file:
        output_file.write(bioawk_output)
    print(f"Processed {file_path}, output written to {output_file_path}")

    # Header mapping (always *_chrIDs.txt)
    mapping_filename = f"{file_prefix}_chrIDs.txt"
    mapping_file_path = os.path.join(out_dir if out_dir else ".", mapping_filename)
    with open(mapping_file_path, 'w') as mapping_file:
        mapping_file.writelines(header_mapping)
    print(f"Header mapping for {file_path} written to {mapping_file_path}")

    # Pass file update if present
    pass_file = f'{file_prefix}.pass.list'
    if pass_file in pass_files:
        process_pass_file(pass_file, mapping_file_path, out_dir, out_suffix)


def process_pass_file(pass_file, mapping_file_path, out_dir, out_suffix):
    if not os.path.isfile(pass_file):
        print(f"[WARN] Pass file not found: {pass_file}")
        return
    if not os.path.isfile(mapping_file_path):
        print(f"[WARN] Mapping file not found for pass processing: {mapping_file_path}")
        return

    mappings = {}
    with open(mapping_file_path, 'r') as mapping_file:
        for line in mapping_file:
            parts = line.strip().rsplit('\t', 1)
            if len(parts) != 2:
                continue
            old = parts[0]
            new = parts[1]
            mappings[old] = new

    with open(pass_file, 'r') as pf:
        pass_lines = pf.readlines()

    new_pass_lines = []
    for line in pass_lines:
        columns = line.strip().split('\t')
        if not columns:
            continue
        old_id_part = columns[0].split(':')[0]
        if old_id_part in mappings:
            columns[0] = columns[0].replace(old_id_part, mappings[old_id_part], 1)
        new_pass_lines.append('\t'.join(columns) + '\n')

    suffix = "" if out_suffix == "disable" else out_suffix
    prefix = os.path.splitext(os.path.splitext(os.path.basename(pass_file))[0])[0]
    output_pass_filename = f"{prefix}{suffix}.pass.list"
    output_pass_file = os.path.join(out_dir if out_dir else ".", output_pass_filename)

    with open(output_pass_file, 'w') as opf:
        opf.writelines(new_pass_lines)
    print(f"Processed {pass_file}, output written to {output_pass_file}")


def generate_jcvi_list(genome_prefixes, out_dir):
    path = os.path.join(out_dir if out_dir else ".", 'jcvi_list.txt')
    with open(path, 'w') as jcvi_file:
        for pair in combinations(genome_prefixes, 2):
            jcvi_file.write('\t'.join(pair) + '\n')
    print(f"Generated {path} with all pairwise relationships.")


def write_prefix_list(filename, prefixes, out_dir):
    path = os.path.join(out_dir if out_dir else ".", filename)
    with open(path, 'w') as file:
        for prefix in prefixes:
            file.write(prefix + '\n')
    print(f"Wrote {path}")


def main(genomes, pass_files, processes, out_dir, out_suffix):
    # Safety: disallow disabling suffix unless writing to a separate output directory
    if out_suffix == "disable" and not out_dir:
        raise SystemExit(
            "[FATAL] -out_suffix disable is only allowed when writing to a separate directory via -out_dir. "
            "Otherwise output could overwrite input."
        )

    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    genome_prefixes = [os.path.splitext(os.path.basename(genome))[0] for genome in genomes]
    pass_prefixes = [os.path.splitext(os.path.splitext(os.path.basename(pass_file))[0])[0] for pass_file in pass_files]

    errors = []
    if processes < 1:
        processes = 1

    with ProcessPoolExecutor(max_workers=processes) as executor:
        future_map = {
            executor.submit(process_file, genome_file, pass_files, out_dir, out_suffix): genome_file
            for genome_file in genomes
        }
        for fut in as_completed(future_map):
            genome_file = future_map[fut]
            try:
                fut.result()
            except Exception as e:
                errors.append((genome_file, str(e)))
                print(f"[ERROR] Failed processing {genome_file}: {e}")

    generate_jcvi_list(genome_prefixes, out_dir)
    write_prefix_list('genome_list.txt', genome_prefixes, out_dir)
    write_prefix_list('pass_list.txt', pass_prefixes, out_dir)

    if errors:
        print("\nThe following genomes failed to process:")
        for gf, msg in errors:
            print(f" - {gf}: {msg}")


if __name__ == "__main__":
    parser = ArgumentParser(
        description=(
            "Rename FASTA headers and update pass files based on chromosome numbers, "
            "with optional parallel processing and flexible output control."
        )
    )
    parser.add_argument("-genomes", nargs='+', help="FASTA files to process.", required=True)
    parser.add_argument("-pass_files", nargs='*', help="Pass files to update based on header mappings.", default=[])
    parser.add_argument("-processes", type=int, default=1, help="Max number of genomes to process concurrently (default: 1).")
    parser.add_argument(
        "-out_dir",
        type=str,
        default=None,
        help="Directory to write all output files. If it doesn't exist, it will be created."
    )
    parser.add_argument(
        "-out_suffix",
        type=str,
        default="_mod",
        help=(
            "Suffix to append to output FASTA and pass.list filenames (default: _mod). "
            "Use '_' to append an underscore, etc. Use 'disable' to omit the suffix entirely "
            "(ONLY allowed when using -out_dir to avoid overwriting inputs)."
        )
    )

    args = parser.parse_args()
    main(args.genomes, args.pass_files, args.processes, args.out_dir, args.out_suffix)
