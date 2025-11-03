#!/usr/bin/env python3
"""
synLTR multi-genome pipeline wrapper

Steps automated:
(1) Clean genomes with fasta_renamer_diploid.py
(2) cd into DIR_NAME produced in (1)
(3) Liftover proteins for each genome with liftover.py
(4) Generate pairwise anchors with jcvi_diploid.py
(5) For each anchors file, run:
    (a)  anchor_builder.py
    (b)  gene_coords_extractor_all4.py -mcscan
    (c)  anchor_coord_subtracter.py
    (d)  anchor_coord_subtracter.py (second pass)
    (e)  anchor_coord_consolidator.py -t MIN_BLOCK_SIZE --stitch-gaps
(6) Cleanup intermediate per-pair files:
       *.anchors.clean
       *.anchors.coords
       *.anchors.coords.polished
       *.anchors.coords.polished2

Features:
- Resume by detecting existing outputs (uses sentinel .done files for major steps and existence checks for consolidated outputs).
- --overwrite to force re-running steps.
- Parallel per-pair anchors processing, capped by --threads.
- Verbose mode to stream underlying tool output; quiet mode prints high-level progress.

Author: Chris Benson & ChatGPT
"""

import argparse
import itertools
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

# --------------------------- Utilities ---------------------------

def echo(msg, verbose=False):
    print(msg, flush=True)

def run_cmd(cmd, cwd=None, verbose=False):
    """
    Run a command list `cmd`. In quiet mode, suppress child output but still raise on failure.
    """
    if verbose:
        proc = subprocess.run(cmd, cwd=cwd, check=True)
    else:
        with open(os.devnull, "wb") as devnull:
            proc = subprocess.run(cmd, cwd=cwd, stdout=devnull, stderr=devnull, check=True)
    return proc.returncode

def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def stem_without_ext(p: Path):
    # Strip only the last suffix (.fa, .fasta, .fa.gz is not expected here)
    return p.name.rsplit(".", 1)[0]

def sentinel(path: Path) -> Path:
    return path.with_suffix(path.suffix + ".done") if path.suffix else path.with_name(path.name + ".done")

# --------------------------- Pipeline Steps ---------------------------

def step1_clean_genomes(genome_paths, script_dir: Path, out_dir: Path, threads: int, overwrite: bool, verbose: bool):
    """
    Run fasta_renamer_diploid.py to create cleaned copies in out_dir.
    Creates sentinel: out_dir/.clean.done
    """
    done = out_dir / ".clean.done"
    if done.exists() and not overwrite:
        echo(f"[1/6] Clean genomes: already done -> {done}", verbose)
        return

    ensure_dir(out_dir)
    fas = [str(p) for p in genome_paths]
    cmd = [
        sys.executable, str(script_dir / "fasta_renamer_diploid.py"),
        "-genomes", *fas,
        "-processes", str(threads),
        "-out_dir", str(out_dir),
        "-out_suffix", "disable",
    ]
    echo(f"[1/6] Cleaning genomes -> {out_dir}", verbose)
    run_cmd(cmd, cwd=script_dir, verbose=verbose)
    done.touch()

def step3_liftover_per_genome(clean_genome_names, script_dir: Path, work_dir: Path,
                              ref_protein_fa: Path, miniprot_outn: int, threads: int,
                              overwrite: bool, verbose: bool):
    """
    Run liftover.py once for ALL genomes in one command, as per your example.
    Creates sentinel: work_dir/.liftover.done
    """
    done = work_dir / ".liftover.done"
    if done.exists() and not overwrite:
        echo("[3/6] Liftover: already done", verbose)
        return

    cmd = [
        sys.executable, str(script_dir / "liftover.py"),
        "--genome", *clean_genome_names,
        "--reference", str(ref_protein_fa),
        "--outn", str(miniprot_outn),
        "--outs", "0.95",
        "--outc", "0.9",
        "--threads", str(threads),
        "--cdhit",
        "--outputs", "pep", "bed",
#        "--TEsorter",
    ]
    echo("[3/6] Running liftover.py", verbose)
    run_cmd(cmd, cwd=work_dir, verbose=verbose)
    done.touch()

def compute_jcvi_cpus_and_inputs(threads: int, n_inputs: int):
    """
    Implements:
      if THREADS < NUMBER_OF_GENOME_INPUTS:
          THREADS_DIVIDED_BY_NUMBER_OF_GENOME_INPUTS=1
          NUMBER_OF_GENOME_INPUTS=THREADS
      else:
          THREADS_DIVIDED_BY_NUMBER_OF_GENOME_INPUTS = THREADS // NUMBER_OF_GENOME_INPUTS (min 1)
    """
    if threads < n_inputs:
        return threads, 1  # p = threads, cpus=1
    else:
        cpus = max(threads // n_inputs, 1)
        return n_inputs, cpus

def step4_jcvi(work_dir: Path, script_dir: Path, genome_stems, threads: int, cscore: float,
               overwrite: bool, verbose: bool):
    """
    Run jcvi_diploid.py to produce pairwise anchors.
    Creates sentinel: work_dir/.jcvi.done
    """
    done = work_dir / ".jcvi.done"
    if done.exists() and not overwrite:
        echo("[4/6] jcvi_diploid: already done", verbose)
        return

    n_inputs_req = len(genome_stems)
    p, cpus = compute_jcvi_cpus_and_inputs(threads, n_inputs_req)

    cmd = [
        sys.executable, str(script_dir / "jcvi_diploid.py"),
        "-p", str(p),
        "--cpus", str(cpus),
        "--keep", "pdf", "anchors",
        "--prot",
        "--cscore", str(cscore),
    ]
    echo(f"[4/6] Running jcvi_diploid.py with -p {p} --cpus {cpus} --cscore {cscore}", verbose)
    run_cmd(cmd, cwd=work_dir, verbose=verbose)

    # Basic existence check: at least one anchors file should exist
    pairs = list(itertools.combinations(genome_stems, 2))
    expected_any = work_dir / f"{pairs[0][0]}.{pairs[0][1]}.anchors"
    if not expected_any.exists():
        raise RuntimeError("jcvi_diploid.py finished but no anchors files were found.")

    # Collect JCVI dotplot PDFs into a dedicated folder
    dotplot_dir = work_dir / "dotplots"
    ensure_dir(dotplot_dir)
    for pdf in work_dir.glob("*.pdf"):
        dest = dotplot_dir / pdf.name
        try:
            if dest.exists():
                # Overwrite to reflect latest run (consistent with --overwrite semantics elsewhere)
                dest.unlink()
            shutil.move(str(pdf), str(dest))
            if verbose:
                echo(f"    moved {pdf.name} -> dotplots/", verbose)
        except Exception as e:
            echo(f"    warning: could not move {pdf.name} to dotplots/: {e}", verbose=True)
    done.touch()

def anchors_pairs(genome_stems):
    return list(itertools.combinations(genome_stems, 2))

def consolidated_output_for_pair(work_dir: Path, g1: str, g2: str) -> Path:
    return work_dir / f"{g1}.{g2}.anchors.coords.polished2.consolidated"

def step5_process_pair(work_dir: Path, script_dir: Path, g1: str, g2: str, min_block_size: int,
                       overwrite: bool, verbose: bool):
    """
    Process a single pair through (a)–(e).
    """
    base = f"{g1}.{g2}.anchors"
    anchors = work_dir / f"{base}"
    if not anchors.exists():
        raise FileNotFoundError(f"Missing anchors file: {anchors}")

    out_e = consolidated_output_for_pair(work_dir, g1, g2)
    if out_e.exists() and not overwrite:
        if verbose:
            echo(f"[5/6] Pair {g1}.{g2}: consolidated exists, skipping", verbose)
        return out_e

    # (a) anchor_builder
    a_out = work_dir / f"{base}.clean"
    cmd_a = [sys.executable, str(script_dir / "anchor_builder.py"), str(anchors)]
    echo(f"[5/6] (a) anchor_builder -> {a_out.name}", verbose)
    with open(a_out, "w") as fout:
        if verbose:
            subprocess.run(cmd_a, cwd=work_dir, check=True, stdout=fout)
        else:
            with open(os.devnull, "wb") as devnull:
                subprocess.run(cmd_a, cwd=work_dir, check=True, stdout=fout, stderr=devnull)

    # (b) gene_coords_extractor_all4 -mcscan
    b_out = work_dir / f"{base}.coords"
    cmd_b = [sys.executable, str(script_dir / "gene_coords_extractor_all4.py"), "-mcscan", str(a_out)]
    echo(f"[5/6] (b) gene_coords_extractor_all4 -> {b_out.name}", verbose)
    with open(b_out, "w") as fout:
        if verbose:
            subprocess.run(cmd_b, cwd=work_dir, check=True, stdout=fout)
        else:
            with open(os.devnull, "wb") as devnull:
                subprocess.run(cmd_b, cwd=work_dir, check=True, stdout=fout, stderr=devnull)

    # (c) anchor_coord_subtracter
    c_out = work_dir / f"{base}.coords.polished"
    cmd_c = [sys.executable, str(script_dir / "anchor_coord_subtracter.py"), str(b_out), str(c_out)]
    echo(f"[5/6] (c) anchor_coord_subtracter -> {c_out.name}", verbose)
    run_cmd(cmd_c, cwd=work_dir, verbose=verbose)

    # (d) anchor_coord_subtracter (second pass)
    d_out = work_dir / f"{base}.coords.polished2"
    cmd_d = [sys.executable, str(script_dir / "anchor_coord_subtracter.py"), str(c_out), str(d_out)]
    echo(f"[5/6] (d) anchor_coord_subtracter (2nd) -> {d_out.name}", verbose)
    run_cmd(cmd_d, cwd=work_dir, verbose=verbose)

    # (e) anchor_coord_consolidator -t MIN_BLOCK_SIZE --stitch-gaps
    cmd_e = [
        sys.executable, str(script_dir / "anchor_coord_consolidator.py"),
        "-t", str(min_block_size),
        str(d_out),
        "--stitch-gaps", # Seems to work now, but if funky results, check here.
    ]
    echo(f"[5/6] (e) anchor_coord_consolidator -> {out_e.name}", verbose)
    with open(out_e, "w") as fout:
        if verbose:
            subprocess.run(cmd_e, cwd=work_dir, check=True, stdout=fout)
        else:
            with open(os.devnull, "wb") as devnull:
                subprocess.run(cmd_e, cwd=work_dir, check=True, stdout=fout, stderr=devnull)

    return out_e

def step5_parallel(work_dir: Path, script_dir: Path, genome_stems, min_block_size: int,
                   threads: int, overwrite: bool, verbose: bool):
    """
    Process all pairs with up to `threads` concurrent workers.
    """
    pairs = anchors_pairs(genome_stems)
    if not pairs:
        raise RuntimeError("Need at least two genomes to build anchors.")

    max_workers = max(1, min(threads, len(pairs)))
    echo(f"[5/6] Processing {len(pairs)} anchor pairs with up to {max_workers} workers", verbose)

    results = {}
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        fut2pair = {
            ex.submit(step5_process_pair, work_dir, script_dir, g1, g2, min_block_size, overwrite, verbose): (g1, g2)
            for (g1, g2) in pairs
        }
        for fut in as_completed(fut2pair):
            g1, g2 = fut2pair[fut]
            try:
                out = fut.result()
                results[(g1, g2)] = out
                if verbose:
                    echo(f"    ✓ Finished {g1}.{g2}", verbose)
            except Exception as e:
                echo(f"    ✗ Pair {g1}.{g2} failed: {e}", verbose=True)
                raise
    return results

def step6_cleanup(work_dir: Path, genome_stems, verbose: bool):
    """
    Remove intermediate files per pair:
      *.anchors.clean
      *.anchors.coords
      *.anchors.coords.polished
      *.anchors.coords.polished2
    """
    echo("[6/6] Cleanup intermediates", verbose)
    patterns = [
        ".anchors.clean",
        ".anchors.coords",
        ".anchors.coords.polished",
        ".anchors.coords.polished2",
    ]
    for g1, g2 in anchors_pairs(genome_stems):
        base = work_dir / f"{g1}.{g2}"
        for suff in patterns:
            p = Path(str(base) + suff)
            if p.exists():
                try:
                    p.unlink()
                    if verbose:
                        echo(f"    removed {p.name}", verbose)
                except Exception as e:
                    echo(f"    warning: could not remove {p.name}: {e}", verbose=True)

# --------------------------- Main ---------------------------

def parse_args():
    ap = argparse.ArgumentParser(description="Wrapper for synLTR multi-genome pipeline.")
    ap.add_argument("--overwrite", action="store_true",
                    help="Overwrite existing outputs; otherwise resume/skip finished steps.")
    ap.add_argument("--script_dir", required=True,
                    help="Directory containing synLTR scripts (fasta_renamer_diploid.py, liftover.py, etc.).")
    ap.add_argument("--verbose", action="store_true", help="Stream tool output.")
    ap.add_argument("--genomes", nargs="+", required=True,
                    help="One or more genomic FASTA files.")
    ap.add_argument("--threads", type=int, default=4, help="Total threads to use (also caps per-pair concurrency).")
    ap.add_argument("--dir_name", required=True, help="Output directory name produced by renamer; subsequent steps run inside it.")
    ap.add_argument("--protein_fa", required=True, help="Trusted proteome FASTA to use for liftover.")
    ap.add_argument("--miniprot_outn", type=int, default=1,
                    help="miniprot -n value (polyploids may need >1). Default: 1")
    ap.add_argument("--c-score", type=float, default=0.99,
                    help="cscore for jcvi_diploid; 0.99 ~ RBH-like strictness. Loosen for polyploids.")
    ap.add_argument("--min-block-size", type=int, default=15000,
                    help="Minimum block size for consolidation; smaller blocks may be stitched.")
    return ap.parse_args()

def main():
    args = parse_args()

    # Resolve paths
    script_dir = Path(args.script_dir).resolve()
    out_dir = Path(args.dir_name).resolve()
    ref_protein_fa = Path(args.protein_fa).resolve()
    threads = max(1, int(args.threads))

    # Validate inputs
    genome_paths = [Path(g).resolve() for g in args.genomes]
    if len(genome_paths) < 2:
        raise SystemExit("Error: provide at least two genomes via --genomes")

    for p in genome_paths:
        if not p.exists():
            raise SystemExit(f"Missing genome file: {p}")

    if not (script_dir / "fasta_renamer_diploid.py").exists():
        echo("Warning: fasta_renamer_diploid.py not found in script_dir; ensure path is correct.", verbose=True)
    if not ref_protein_fa.exists():
        raise SystemExit(f"Missing --protein_fa file: {ref_protein_fa}")

    # (1) Clean genomes
    step1_clean_genomes(genome_paths, script_dir, out_dir, threads, args.overwrite, args.verbose)

    # (2) Work directory is out_dir created above
    work_dir = out_dir
    if not work_dir.exists():
        raise SystemExit(f"Expected output directory not found: {work_dir}")

    # Cleaned genome filenames are copied into work_dir with out_suffix disabled (same base names).
    # We'll refer to them by filename (not path) for downstream steps executed in work_dir.
    cleaned_names = [p.name for p in genome_paths]
    genome_stems = [stem_without_ext(Path(n)) for n in cleaned_names]

    # (3) Liftover
    step3_liftover_per_genome(cleaned_names, script_dir, work_dir, ref_protein_fa,
                              args.miniprot_outn, threads, args.overwrite, args.verbose)

    # (4) jcvi anchors
    step4_jcvi(work_dir, script_dir, genome_stems, threads, args.c_score,
               args.overwrite, args.verbose)

    # (5) Per-pair processing in parallel (up to `threads` concurrent pairs)
    step5_parallel(work_dir, script_dir, genome_stems, args.min_block_size,
                   threads, args.overwrite, args.verbose)

    # (6) Cleanup intermediates
    step6_cleanup(work_dir, genome_stems, args.verbose)

    echo("Pipeline complete.", verbose=True)

if __name__ == "__main__":
    main()

