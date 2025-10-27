#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import sys
import re
import random
from collections import defaultdict

def parse_args():
    ap = argparse.ArgumentParser(
        description="Build orthogroup FASTAs (pep & cds) from an orthogroups file."
    )
    ap.add_argument("--orthogroups", required=True,
                    help="TSV/whitespace-delimited file; each line is an orthogroup.")
    ap.add_argument("--peps", nargs="+", required=True,
                    help="Protein FASTA file(s). Shell globs are fine (expanded by the shell).")
    ap.add_argument("--cds", nargs="+", required=True,
                    help="CDS FASTA file(s). Shell globs are fine (expanded by the shell).")
    ap.add_argument("--seed", type=int, default=None,
                    help="Seed for random selection of duplicate species within an orthogroup.")
    ap.add_argument("--outdir", default="orthogroups",
                    help="Base output directory (default: ./orthogroups).")
    ap.add_argument("--width", type=int, default=60,
                    help="FASTA line width (default: 60).")
    return ap.parse_args()

def parse_fasta_many(paths):
    seqs = {}
    seen = set()
    for p in paths:
        if not os.path.exists(p):
            print(f"[WARN] FASTA not found: {p}", file=sys.stderr)
            continue
        with open(p, "r") as fh:
            header = None
            chunks = []
            for line in fh:
                if not line:
                    continue
                if line.startswith(">"):
                    if header is not None:
                        sid = header.split()[0]
                        if sid in seqs:
                            if sid not in seen:
                                print(f"[WARN] duplicate sequence ID '{sid}' encountered; last one wins.",
                                      file=sys.stderr)
                                seen.add(sid)
                        seqs[sid] = "".join(chunks).replace(" ", "").replace("\t", "")
                    header = line[1:].strip()
                    chunks = []
                else:
                    chunks.append(line.strip())
            if header is not None:
                sid = header.split()[0]
                if sid in seqs:
                    if sid not in seen:
                        print(f"[WARN] duplicate sequence ID '{sid}' encountered; last one wins.",
                              file=sys.stderr)
                        seen.add(sid)
                seqs[sid] = "".join(chunks).replace(" ", "").replace("\t", "")
    return seqs

_species_re = re.compile(r"^([A-Za-z]+)")

def species_from_id(gene_id):
    """
    Species prefix = leading letters before the first digit.
    E.g., 'Bdact000001_2' -> 'Bdact'
    """
    m = _species_re.match(gene_id)
    return m.group(1) if m else None

def wrap(seq, width):
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))

def ensure_dirs(base):
    pep_dir = os.path.join(base, "pep")
    cds_dir = os.path.join(base, "cds")
    os.makedirs(pep_dir, exist_ok=True)
    os.makedirs(cds_dir, exist_ok=True)
    return pep_dir, cds_dir

def load_orthogroups(path):
    if not os.path.exists(path):
        sys.exit(f"[ERROR] Orthogroups file not found: {path}")
    groups = []
    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # split on tabs or any whitespace
            genes = re.split(r"\s+", line)
            groups.append(genes)
    return groups

def pick_representatives(genes, rng):
    """
    For a list of gene IDs in one orthogroup, return dict species->chosen_id.
    If multiple genes share the same species prefix, choose one at random (using rng).
    """
    by_species = defaultdict(list)
    for g in genes:
        sp = species_from_id(g)
        if not sp:
            print(f"[WARN] could not parse species prefix from '{g}'. Skipping this entry.", file=sys.stderr)
            continue
        by_species[sp].append(g)

    chosen = {}
    for sp, ids in by_species.items():
        if len(ids) == 1:
            chosen[sp] = ids[0]
        else:
            chosen[sp] = rng.choice(ids)
    return chosen

def write_fasta(path, entries, seqs, width):
    """
    entries: list of (gene_id, label_for_header)
    seqs: dict id->sequence
    """
    wrote_any = False
    with open(path, "w") as out:
        for gid, label in entries:
            if gid not in seqs:
                print(f"[WARN] missing sequence for '{gid}' in {os.path.basename(path)}; skipping.", file=sys.stderr)
                continue
            out.write(f">{label}\n")
            out.write(wrap(seqs[gid], width) + "\n")
            wrote_any = True
    return wrote_any

def main():
    args = parse_args()

    rng = random.Random(args.seed)
    pep_dir, cds_dir = ensure_dirs(args.outdir)

    print("[INFO] Reading protein FASTAs...", file=sys.stderr)
    pep_seqs = parse_fasta_many(args.peps)
    print(f"[INFO] Loaded {len(pep_seqs)} protein sequences.", file=sys.stderr)

    print("[INFO] Reading CDS FASTAs...", file=sys.stderr)
    cds_seqs = parse_fasta_many(args.cds)
    print(f"[INFO] Loaded {len(cds_seqs)} CDS sequences.", file=sys.stderr)

    groups = load_orthogroups(args.orthogroups)
    print(f"[INFO] Loaded {len(groups)} orthogroups.", file=sys.stderr)

    n_written_pep = n_written_cds = 0

    for idx, genes in enumerate(groups, start=1):
        og_id = f"OG{idx:05d}"
        # choose one gene per species (consistent for pep & cds)
        chosen = pick_representatives(genes, rng)

        # Order headers deterministically by species name to keep files tidy
        ordered_species = sorted(chosen.keys())

        pep_path = os.path.join(pep_dir, f"{og_id}.fa")
        cds_path = os.path.join(cds_dir, f"{og_id}.fa")

        # We keep the gene IDs as the FASTA headers, as in your example
        pep_entries = [(chosen[sp], chosen[sp]) for sp in ordered_species]
        cds_entries = [(chosen[sp], chosen[sp]) for sp in ordered_species]

        wrote_pep = write_fasta(pep_path, pep_entries, pep_seqs, args.width)
        wrote_cds = write_fasta(cds_path, cds_entries, cds_seqs, args.width)

        if wrote_pep:
            n_written_pep += 1
        else:
            # If nothing was written, remove the empty file to avoid clutter
            try:
                os.remove(pep_path)
            except OSError:
                pass

        if wrote_cds:
            n_written_cds += 1
        else:
            try:
                os.remove(cds_path)
            except OSError:
                pass

    print(f"[INFO] Done. Wrote {n_written_pep} protein orthogroup FASTAs to {pep_dir}/", file=sys.stderr)
    print(f"[INFO] Done. Wrote {n_written_cds} CDS orthogroup FASTAs to {cds_dir}/", file=sys.stderr)

if __name__ == "__main__":
    main()
