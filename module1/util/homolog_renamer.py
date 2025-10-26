#!/usr/bin/env python3
"""
Rename + restrand query chromosomes based on a 'best homolog' TSV.

Adds --allow-sca to optionally treat "_sca" as a chromosome tag (like "_chr").
"""

import os
import sys
import argparse
import shutil
from collections import defaultdict

def parse_args():
    p = argparse.ArgumentParser(
        description="Rename + restrand query chromosomes based on a 'best homolog' TSV."
    )
    p.add_argument("tsv", help="Path to the TSV file (e.g. best_homolog2.2.tsv).")
    p.add_argument("genomes_dir", help="Directory containing genome FASTA files named [accession].fa")
    p.add_argument(
        "--allow-sca",
        action="store_true",
        help="Also treat IDs containing '_sca' as valid chromosome identifiers (e.g. Etef_sca8).",
    )
    return p.parse_args()

def reverse_complement(seq: str) -> str:
    comp = {
        "A": "T", "T": "A", "C": "G", "G": "C",
        "a": "t", "t": "a", "c": "g", "g": "c",
        "N": "N", "n": "n"
    }
    return "".join(comp.get(b, b) for b in reversed(seq))

def wrap_fasta(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))

def split_acc_chr(full_id: str, allowed_tags):
    """
    Try to split an identifier like 'Zmays_chr1' or 'Etef_sca8' into (acc, chr).
    Returns (acc, chr) if any allowed tag is present, else None.
    """
    for tag in allowed_tags:
        if tag in full_id:
            acc, chr_part = full_id.rsplit(tag, 1)
            return acc, chr_part
    return None

def read_best_homolog_tsv(path, allowed_tags):
    """
    Parse the TSV. Build mapping records and sets of accessions.
    """
    mapping_records = []
    ref_accessions = set()
    query_accessions = set()

    with open(path, "r") as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line:
                continue
            tokens = line.split("\t")
            if len(tokens) < 2:
                continue

            ref_full = tokens[0].strip()   # e.g. "Zmays_chr1" or "Etef_sca8"
            ref_parts = split_acc_chr(ref_full, allowed_tags)
            if not ref_parts:
                sys.stderr.write(
                    f"WARNING: expected one of {allowed_tags} in reference '{ref_full}', skipping line.\n"
                )
                continue
            ref_acc, ref_chr = ref_parts
            ref_accessions.add(ref_acc)

            if tokens[1].upper() == "NA":
                continue

            remainder = tokens[1:]
            if len(remainder) % 2 != 0:
                sys.stderr.write(
                    "ERROR: line does not have an even number of query/strand fields:\n"
                    f"    {line}\n"
                )
                sys.exit(1)

            for i in range(0, len(remainder), 2):
                query_full = remainder[i].strip()   # e.g. "Bdhap_chr3" or "Etef_sca8"
                strand = remainder[i+1].strip()

                q_parts = split_acc_chr(query_full, allowed_tags)
                if not q_parts:
                    sys.stderr.write(
                        f"WARNING: expected one of {allowed_tags} in query '{query_full}', skipping.\n"
                    )
                    continue
                query_acc, query_chr = q_parts
                query_accessions.add(query_acc)

                rec = {
                    "ref_acc": ref_acc,
                    "ref_chr": ref_chr,
                    "query_acc": query_acc,
                    "query_chr": query_chr,
                    "orig_qid": query_full,
                    "strand": strand,
                }
                mapping_records.append(rec)

    return mapping_records, ref_accessions, query_accessions

def assign_suffixes(mapping_records):
    group_totals = defaultdict(int)
    for rec in mapping_records:
        key = (rec["query_acc"], rec["ref_acc"], rec["ref_chr"])
        group_totals[key] += 1

    seen_counts = defaultdict(int)
    mapping_info = defaultdict(dict)

    for rec in mapping_records:
        qacc = rec["query_acc"]
        racc = rec["ref_acc"]
        rchr = rec["ref_chr"]
        orig_qid = rec["orig_qid"]
        strand = rec["strand"]

        key = (qacc, racc, rchr)
        total = group_totals[key]

        if total == 1:
            suffix = ""
        else:
            seen_counts[key] += 1
            idx = seen_counts[key]
            if idx > 26:
                sys.stderr.write(f"ERROR: More than 26 mappings for {key}; cannot assign unique letter.\n")
                sys.exit(1)
            suffix = chr(ord("A") + idx - 1)

        new_qid = f"{qacc}_chr{rchr}{suffix}"
        mapping_info[qacc][orig_qid] = {"new_qid": new_qid, "strand": strand}

    return mapping_info

def copy_reference_fastas(ref_accessions, genomes_dir):
    for ref_acc in sorted(ref_accessions):
        src = os.path.join(genomes_dir, f"{ref_acc}.fa")
        dest = f"{ref_acc}_mod.fa"
        if not os.path.exists(src):
            sys.stderr.write(f"WARNING: Reference FASTA '{src}' not found. Skipping copy.\n")
            continue
        shutil.copyfile(src, dest)

def process_query_fasta(query_acc, mapping_for_acc, genomes_dir):
    src = os.path.join(genomes_dir, f"{query_acc}.fa")
    if not os.path.exists(src):
        sys.stderr.write(f"WARNING: Query FASTA '{src}' not found. Skipping {query_acc}.\n")
        return

    out_path = f"{query_acc}_mod.fa"
    out_fh = open(out_path, "w")

    with open(src, "r") as fh:
        curr_header = None
        curr_seq_parts = []
        for raw_line in fh:
            line = raw_line.rstrip("\n")
            if line.startswith(">"):
                if curr_header is not None:
                    header_id = curr_header.split()[0]
                    if header_id in mapping_for_acc:
                        info = mapping_for_acc[header_id]
                        seq = "".join(curr_seq_parts)
                        if info["strand"] == "-":
                            seq = reverse_complement(seq)
                        out_fh.write(f">{info['new_qid']}\n{wrap_fasta(seq)}\n")
                curr_header = line[1:]
                curr_seq_parts = []
            else:
                curr_seq_parts.append(line)

        if curr_header is not None:
            header_id = curr_header.split()[0]
            if header_id in mapping_for_acc:
                info = mapping_for_acc[header_id]
                seq = "".join(curr_seq_parts)
                if info["strand"] == "-":
                    seq = reverse_complement(seq)
                out_fh.write(f">{info['new_qid']}\n{wrap_fasta(seq)}\n")

    out_fh.close()

def main():
    args = parse_args()
    tsv = args.tsv
    genomes_dir = args.genomes_dir.rstrip("/")

    if not os.path.exists(tsv):
        sys.stderr.write(f"ERROR: TSV file '{tsv}' does not exist.\n")
        sys.exit(1)
    if not os.path.isdir(genomes_dir):
        sys.stderr.write(f"ERROR: genomes_dir '{genomes_dir}' is not a directory.\n")
        sys.exit(1)

    allowed_tags = ["_chr"] + (["_sca"] if args.allow_sca else [])

    mapping_records, ref_accessions, query_accessions = read_best_homolog_tsv(tsv, allowed_tags)
    mapping_info = assign_suffixes(mapping_records)
    copy_reference_fastas(ref_accessions, genomes_dir)

    for qacc in sorted(query_accessions):
        if qacc not in mapping_info:
            continue
        process_query_fasta(qacc, mapping_info[qacc], genomes_dir)

    print("Done. Created *_mod.fa files for references (copied) and queries (renamed/RC’d).")

if __name__ == "__main__":
    main()
