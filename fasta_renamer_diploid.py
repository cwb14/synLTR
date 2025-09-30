#!/usr/bin/env python3
import re
import os
import subprocess
from argparse import ArgumentParser
from itertools import combinations
from concurrent.futures import ProcessPoolExecutor, as_completed


def rename_header(header, file_prefix):
    # chrom labels: chr/Chrom/Chromosome/CHR/etc with optional separators and optional colon+space
    chrom_regex = re.compile(r'(Chr|chr|Chro|chro|Chrom|chrom|Chromosome|chromosome|CHROMOSOME|CHR|CHRO|CHROM)[ _-]*:? ?(\d+|[XYZW])')
    match = chrom_regex.search(header)
    if match:
        # Convert chromosome number to integer to remove leading zeros
        chrom_number = match.group(2)
        if chrom_number.isdigit():
            chrom_number = str(int(chrom_number))  # Remove leading zeros
        return f'>{file_prefix}_chr{chrom_number}', f'{file_prefix}_chr{chrom_number}'
    return None, None


def process_file(file_path, pass_files, out_dir, out_suffix):
    file_prefix = os.path.splitext(os.path.basename(file_path))[0]
    if not os.path.isfile(file_path):
        print(f"[ERROR] Genome file not found: {file_path}")
        return

    with open(file_path, 'r') as fasta_file:
        lines = fasta_file.readlines()

    sequences = {}
    current_header = None
    for line in lines:
        if line.startswith('>'):
            current_header = line.strip()
            sequences[current_header] = []
        elif current_header:
            sequences[current_header].append(line.strip())

    sequence_lengths = {header: sum(len(seq) for seq in seqs) for header, seqs in sequences.items()}

    # Separate sequences into those with chr designations and those without
    chr_sequences = {header: seqs for header, seqs in sequences.items() if rename_header(header, file_prefix)[0]}
    sca_sequences = {header: seqs for header, seqs in sequences.items() if not rename_header(header, file_prefix)[0]}

    # Sort sca sequences by size
    sorted_sca_sequences = sorted(sca_sequences.items(), key=lambda item: -sequence_lengths[item[0]])

    new_lines = []
    used_headers = set()
    fallback_counter = 1
    header_mapping = []

    # Process chr sequences without sorting
    for header, seqs in chr_sequences.items():
        new_header, new_header_full = rename_header(header, file_prefix)
        old_header = header.strip()[1:]
        if new_header and new_header not in used_headers:
            used_headers.add(new_header)
            new_lines.append(f'{new_header}\n')
            header_mapping.append(f'{old_header}\t{new_header.strip(">")}\n')
        new_lines.append(''.join(seqs) + '\n')

    # Process sorted sca sequences with length check
    for header, seqs in sorted_sca_sequences:
        while f'>{file_prefix}_sca{fallback_counter}' in used_headers:
            fallback_counter += 1
        new_header_full = f'>{file_prefix}_sca{fallback_counter}'
        if len(new_header_full) > 14:
            print(f"Excluding sequence with header {new_header_full} due to length > 13 characters.")
            break
        old_header = header.strip()[1:]
        new_lines.append(f'{new_header_full}\n')
        used_headers.add(new_header_full)
        fallback_counter += 1
        header_mapping.append(f'{old_header}\t{new_header_full.strip(">")}\n')
        new_lines.append(''.join(seqs) + '\n')

    # Feed the new_lines output directly to bioawk via stdin
    bioawk_command = ["bioawk", "-c", "fastx", '{print ">"$name; print $seq}']
    try:
        process = subprocess.Popen(bioawk_command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=False)
    except FileNotFoundError:
        raise RuntimeError("bioawk not found on PATH. Please install bioawk and ensure it is accessible.")

    fasta_content = ''.join(new_lines).encode('utf-8')
    bioawk_output, bioawk_err = process.communicate(fasta_content)
    if process.returncode != 0:
        err_text = bioawk_err.decode('utf-8', errors='replace') if bioawk_err else ''
        raise RuntimeError(f"bioawk failed for {file_path} with exit code {process.returncode}.\n{err_text}")

    # Prepare output paths
    # - fasta and pass list honor out_suffix (or disabled when allowed)
    # - mapping file keeps *_chrIDs.txt (no suffix), but goes to out_dir if provided
    suffix = "" if out_suffix == "disable" else out_suffix
    fasta_filename = f"{file_prefix}{suffix}.fa"
    output_file_path = os.path.join(out_dir if out_dir else ".", fasta_filename)

    # Write the bioawk processed output directly to the final file
    with open(output_file_path, 'wb') as output_file:
        output_file.write(bioawk_output)
    print(f"Processed {file_path}, output written to {output_file_path}")

    # Writing the header mapping file (always named *_chrIDs.txt)
    mapping_filename = f"{file_prefix}_chrIDs.txt"
    mapping_file_path = os.path.join(out_dir if out_dir else ".", mapping_filename)
    with open(mapping_file_path, 'w') as mapping_file:
        mapping_file.writelines(header_mapping)
    print(f"Header mapping for {file_path} written to {mapping_file_path}")

    # Check if a pass file exists for this file prefix and process it if so
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

    # Prepare output path for pass list (respects out_suffix setting)
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

    # Create out_dir if provided and doesn't exist
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    genome_prefixes = [os.path.splitext(os.path.basename(genome))[0] for genome in genomes]
    pass_prefixes = [os.path.splitext(os.path.splitext(os.path.basename(pass_file))[0])[0] for pass_file in pass_files]

    # Parallel processing of genome files with up to `processes` concurrent workers
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

    # Generate lists after all processing completes (these live in out_dir if provided)
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

    # New options
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
