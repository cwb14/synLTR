#!/usr/bin/env python3
"""solo_ltr.py - end-to-end solo-LTR detection pipeline.

Wraps the three stages into one command:

  1. mask  : hardmask full-length LTR-RT loci (from internal-FASTA headers)
  2. blast : makeblastdb + blastn (consensus AND internal vs masked genome)
  3. filter: TSD-gated, internal-proximity-filtered solo-LTR calls

Inputs
------
  --consensus    LTRs.alns.consensus.fa   (IUPAC LTR consensus, 1 per LTR-RT)
  --internal     LTRs.alns.internal.fa    (headers = full-length LTR-RT coords)
  --genome       reference genome FASTA

Output
------
  <outdir>/solo_ltr.bed                   predicted solo-LTR loci

Stages are cached: re-running with the same inputs skips the mask/db/blast
work and only re-does the (fast) filter.  Use --threads for blastn.

Defaults are the grid-search F1 optimum on the PrinTE benchmark
(blastn word_size 11, pident>=85, qcov>=95, TSD 5bp slop 2, internal flank 25).
"""
import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

import solo_ltr_core as core

HERE = Path(__file__).resolve().parent
SEARCH = HERE / "solo_ltr_search.py"


def log(msg, t0=None):
    if t0 is not None:
        sys.stderr.write(f"[solo_ltr] {msg} ({time.time()-t0:.1f}s)\n")
    else:
        sys.stderr.write(f"[solo_ltr] {msg}\n")
    sys.stderr.flush()


def run(cmd):
    sys.stderr.write("$ " + " ".join(str(x) for x in cmd) + "\n")
    sys.stderr.flush()
    rc = subprocess.call([str(x) for x in cmd])
    if rc != 0:
        raise SystemExit(f"command failed (rc={rc}): {cmd}")


def blast_out(out_dir, query, task, ws, dust, evalue):
    """Recreate solo_ltr_search.blast_cache_path for a known config."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("ss", SEARCH)
    ss = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ss)
    return ss.blast_cache_path(out_dir, query, str(Path(out_dir).parent / "genome.masked.fasta"),
                               task, ws, dust, evalue, "")


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--consensus", required=True)
    p.add_argument("--internal", required=True)
    p.add_argument("--genome", required=True)
    p.add_argument("--outdir", required=True)

    p.add_argument("--task", default="blastn")
    p.add_argument("--word-size", type=int, default=11)
    p.add_argument("--dust", default="no")
    p.add_argument("--evalue", default="10")
    p.add_argument("--threads", type=int, default=8)

    # filter params (defaults = benchmark F1 optimum: F1=0.581, P=0.95, R=0.42)
    # High-precision alternative: --min-qcov-hsp 95  (F1=0.580, P=0.97, 2 FP)
    p.add_argument("--min-pident", type=float, default=85.0)
    p.add_argument("--min-qcov-hsp", type=float, default=90.0)
    p.add_argument("--min-length", type=int, default=100)
    p.add_argument("--max-evalue", type=float, default=1e-10)
    p.add_argument("--no-tsd", action="store_true")
    p.add_argument("--tsd-k", type=int, default=5)
    p.add_argument("--tsd-slop", type=int, default=1)
    p.add_argument("--int-min-pident", type=float, default=80.0)
    p.add_argument("--int-min-qcov-hsp", type=float, default=0.0)
    p.add_argument("--int-min-length", type=int, default=100)
    p.add_argument("--int-max-evalue", type=float, default=1e-10)
    p.add_argument("--flank", type=int, default=25)
    p.add_argument("--no-internal-filter", action="store_true")
    p.add_argument("--merge-slop", type=int, default=50)

    p.add_argument("--truth", default=None, help="PrinTE BED -> report P/R/F1")
    p.add_argument("--force", action="store_true", help="Redo cached mask/blast")
    p.add_argument("--no-nested", action="store_true",
                   help="Disable intra-LTR-RT (nested) solo-LTR recovery "
                        "(consensus vs internal FASTA). Default: enabled.")
    p.add_argument("--nested-min-pident", type=float, default=95.0,
                   help="pident floor for nested (intra-LTR-RT) calls "
                        "(default 95; the main path uses --min-pident)")
    p.add_argument("--no-nested-guard", action="store_true",
                   help="Disable the nested-intact guard (which drops nested calls "
                        "abutting a different internal region's edge)")
    p.add_argument("--nested-guard-flank", type=int, default=25,
                   help="Abutment distance (bp) for the nested-intact guard")
    args = p.parse_args()

    t0 = time.time()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    masked = outdir / "genome.masked.fasta"
    full_bed = outdir / "full_length_ltr.bed"
    blast_dir = outdir / "blast"

    # 1. mask
    log("stage 1/3: hardmask full-length LTR-RTs")
    mask_cmd = [sys.executable, str(SEARCH), "mask",
                "--internal-fa", args.internal, "--genome", args.genome,
                "--bed-out", str(full_bed), "--fasta-out", str(masked)]
    if args.force:
        mask_cmd.append("--force")
    run(mask_cmd)

    # 2. db + blast
    log("stage 2/3: makeblastdb + blastn")
    db_cmd = [sys.executable, str(SEARCH), "db", "--fasta", str(masked)]
    if args.force:
        db_cmd.append("--force")
    run(db_cmd)

    def do_blast(query, db):
        cmd = [sys.executable, str(SEARCH), "blast",
               "--query", query, "--db", str(db), "--out-dir", str(blast_dir),
               "--task", args.task, "--word-size", str(args.word_size),
               "--dust", args.dust, "--evalue", str(args.evalue),
               "--threads", str(args.threads)]
        if args.force:
            cmd.append("--force")
        out = subprocess.check_output([str(x) for x in cmd]).decode().strip().splitlines()
        return out[-1]  # last line = path

    cons_tsv = do_blast(args.consensus, masked)
    int_tsv = do_blast(args.internal, masked)

    # 3. filter
    log("stage 3/3: TSD + internal filter")
    genome = core.load_genome(args.genome)
    cons = core.load_blast(cons_tsv)
    intn = core.load_blast(int_tsv)
    cands = core.make_candidates(
        cons, genome,
        min_pident=args.min_pident, min_qcov=args.min_qcov_hsp,
        min_length=args.min_length, max_evalue=args.max_evalue,
        use_tsd=not args.no_tsd, tsd_k=args.tsd_k, tsd_slop=args.tsd_slop,
        merge_slop=args.merge_slop)
    if not args.no_internal_filter:
        cands = core.filter_internal(
            cands, intn,
            min_pident=args.int_min_pident, min_qcov=args.int_min_qcov_hsp,
            min_length=args.int_min_length, max_evalue=args.int_max_evalue,
            flank=args.flank)

    # 3b. nested (intra-LTR-RT) solo recovery: consensus vs internal FASTA
    if not args.no_nested:
        log("stage 3b: nested (intra-LTR-RT) solo recovery")
        internal_db = blast_dir / "internal_db"
        db_cmd = [sys.executable, str(SEARCH), "db",
                  "--fasta", args.internal, "--out", str(internal_db),
                  "--no-parse-seqids"]
        if args.force:
            db_cmd.append("--force")
        run(db_cmd)

        nested_tsv = do_blast(args.consensus, internal_db)

        # restrict orientation detection to subjects that actually got a hit
        hit_ids = set()
        with open(nested_tsv) as fh:
            for line in fh:
                if not line or line[0] == "#":
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) > core.COL["sseqid"]:
                    hit_ids.add(f[core.COL["sseqid"]])

        orient_map = core.build_orient_map(args.internal, genome, only=hit_ids)
        nested_hsps = core.load_internal_blast(nested_tsv, orient_map)
        nested = core.make_candidates(
            nested_hsps, genome,
            min_pident=args.nested_min_pident, min_qcov=args.min_qcov_hsp,
            min_length=args.min_length, max_evalue=args.max_evalue,
            use_tsd=not args.no_tsd, tsd_k=args.tsd_k, tsd_slop=args.tsd_slop,
            merge_slop=args.merge_slop)
        if not args.no_nested_guard:
            intervals = core.load_internal_intervals(args.internal)
            nested = core.drop_ltr_of_nested_intact(
                nested, intervals, args.nested_guard_flank)
        nb = core.write_bed(nested, str(outdir / "nested_solo.bed"))
        log(f"nested recovery: {nb} calls -> {outdir / 'nested_solo.bed'}")
        cands = core.union_candidates(cands, nested, merge_slop=args.merge_slop)

    out_bed = outdir / "solo_ltr.bed"
    n = core.write_bed(cands, str(out_bed))
    log(f"wrote {n} solo-LTR calls -> {out_bed}", t0)

    if args.truth:
        solo_idx = core.load_truth_solo(args.truth)
        m = core.score(cands, solo_idx)
        log(f"vs truth: detected {m['detected']}/{m['n_solo']}  "
            f"TP={m['tp']} FP={m['fp']} FN={m['fn']}  "
            f"P={m['precision']:.4f} R={m['recall']:.4f} F1={m['f1']:.4f}")


if __name__ == "__main__":
    main()
