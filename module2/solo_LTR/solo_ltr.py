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
  <outdir>/solo_ltr.tsv      predicted solo-LTR loci (rich TSV; cols 1-3 = BED)
  <outdir>/nested_solo.tsv   intra-LTR-RT (nested) solo-LTR subset

Intermediates (masked genome, BLAST db, blast/, full_length_ltr.bed) are
deleted after a successful run; pass --keep-intermediates to retain them.

Stages are cached: re-running with the same inputs skips the mask/db/blast
work and only re-does the (fast) filter.  Use --threads for blastn.

Blast output is bounded at the source: each blastn runs with the path's own
downstream cutoffs as -evalue/-perc_identity/-qcov_hsp_perc (a touch looser, so
the Python filter stays authoritative), so the on-disk TSV holds ~survivors,
not every permissive hit -- a genome-scale permissive blastn can otherwise
reach terabytes.  --loose-blast restores one permissive blast for filter
sweeps.

Defaults are the grid-search F1 optimum on the PrinTE benchmark
(blastn word_size 11, pident>=85, qcov>=95, TSD 5bp slop 2, internal flank 25).

# A thought:
  I can Blast the flanking sequence of the candidate solo-LTRs against the flanking of all other candidates, similarly to the internal blast.
  The flanking sequence of true solo-LTRs should be unique.
  Imagine a genome full of SINEs (formatted '[SINE_tRNA][SINE_tail]').
  Adjacent SINEs will sometimes be falsely labeled as LTR-RT (eg, '[SINE_tRNA][SINE_tail][intergenic][SINE_tRNA][SINE_tail]' instead of '[LTR][internal][LTR]').
  LTRharvest may have erroneously caught a fragment ([SINE_tRNA]) as LTR or it may have caught the whole thing ([SINE_tRNA][SINE_tail]).
  If the candidate solo-LTR is the fragment ('[SINE_tRNA]') then the flanks would contain [SINE_tail]. So, if I blast those flanks agains each other and see matches, its likely [SINE_tail], indicating a false positive. 
  If the whole thing ([SINE_tRNA][SINE_tail]) is labeled as solo-LTR, then this approach wont work and it may be harder to find.
  Im thinking of the dog genome where many false positive LTR-RTs are fragments of SINEs.
  Probably easier to filter these out at the LTR-RT annotation step than the solo-step.
"""
import argparse
import glob
import os
import shutil
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


def cleanup_intermediates(masked, full_bed, blast_dir):
    """Delete the bulky intermediates a successful run no longer needs.

    Removes the masked genome + its BLAST db sidecars (glob masked*),
    full_length_ltr.bed, and the blast/ directory. The original --genome is
    never touched. Cleanup failures warn but do not fail the run.
    """
    for path in glob.glob(str(masked) + "*") + [str(full_bed)]:
        try:
            p = Path(path)
            if p.is_file():
                p.unlink()
        except OSError as ex:
            log(f"cleanup: could not remove {path}: {ex}")
    try:
        if blast_dir.exists():
            shutil.rmtree(blast_dir)
    except OSError as ex:
        log(f"cleanup: could not remove {blast_dir}: {ex}")


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
    p.add_argument("--evalue", default="10",
                   help="blastn -evalue; used ONLY with --loose-blast. Default "
                        "mode derives the per-path blast e-value from "
                        "--max-evalue / --int-max-evalue so the TSV stays small.")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--chunk-threads", type=int, default=4,
                   help="Threads per blastn process; the worker count is "
                        "(threads // chunk-threads) blastn processes run "
                        "concurrently against the shared DB, scaling past "
                        "blastn's internal threading ceiling on many-core hosts. "
                        "Set >= --threads to run a single blast.")
    p.add_argument("--chunk-oversub", type=int, default=16,
                   help="Split the query into (threads // chunk-threads) * "
                        "chunk-oversub length-balanced chunks fed through a pool "
                        "of (threads // chunk-threads) blastn workers. >1 lets a "
                        "finished worker pull queued work instead of idling on a "
                        "straggler, so the step ends with one small chunk rather "
                        "than the single slowest of one-chunk-per-worker. "
                        "1 = static one chunk per worker. Default 16.")

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
    p.add_argument("--loose-blast", action="store_true",
                   help="Disable the blast-level prefilter and run one permissive "
                        "blastn (-evalue from --evalue, no pident/qcov cap). For "
                        "parameter sweeps that vary the Python cutoffs on a single "
                        "cached blast; NOT for genome-scale runs -- the permissive "
                        "TSV can reach terabytes on disk.")
    p.add_argument("--keep-intermediates", action="store_true",
                   help="Keep the masked genome, BLAST db, blast/ dir, and "
                        "full_length_ltr.bed (default: delete them after a "
                        "successful run)")
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

    def do_blast(query, db, *, min_pident, min_qcov, max_evalue):
        # Push the path's downstream cutoffs into blastn so the on-disk TSV is
        # bounded to ~survivors (disk analogue of the streaming RAM prefilter).
        # --loose-blast restores one permissive blast for filter sweeps.
        if args.loose_blast:
            evalue, perc_id, qcov_hsp = args.evalue, None, None
        else:
            ev, perc_id, qcov_hsp = core.blast_prefilter(
                min_pident, min_qcov, max_evalue)
            evalue = f"{ev:.3e}"
        cmd = [sys.executable, str(SEARCH), "blast",
               "--query", query, "--db", str(db), "--out-dir", str(blast_dir),
               "--task", args.task, "--word-size", str(args.word_size),
               "--dust", args.dust, "--evalue", str(evalue),
               "--threads", str(args.threads),
               "--chunk-threads", str(args.chunk_threads),
               "--chunk-oversub", str(args.chunk_oversub)]
        if perc_id is not None:
            cmd += ["--perc-identity", str(perc_id)]
        if qcov_hsp is not None and qcov_hsp > 0:
            cmd += ["--qcov-hsp-perc", str(qcov_hsp)]
        if args.force:
            cmd.append("--force")
        out = subprocess.check_output([str(x) for x in cmd]).decode().strip().splitlines()
        return out[-1]  # last line = path

    cons_tsv = do_blast(args.consensus, masked,
                        min_pident=args.min_pident, min_qcov=args.min_qcov_hsp,
                        max_evalue=args.max_evalue)
    int_tsv = do_blast(args.internal, masked,
                       min_pident=args.int_min_pident,
                       min_qcov=args.int_min_qcov_hsp,
                       max_evalue=args.int_max_evalue)

    # 3. filter
    log("stage 3/3: TSD + internal filter")
    genome = core.load_genome(args.genome)
    # Consensus hits become the solo-LTR calls, so they're held as HSPs -- but
    # streamed with the path's downstream thresholds at load time (bounds memory
    # by survivors, not file size; lossless -- make_candidates applies the same
    # cutoffs). The internal BLAST, used only as a coverage mask, is handled
    # separately below without ever materializing its HSPs.
    cons = core.load_blast(
        cons_tsv, min_pident=args.min_pident, min_qcov=args.min_qcov_hsp,
        min_length=args.min_length, max_evalue=args.max_evalue)
    cands = core.make_candidates(
        cons, genome,
        min_pident=args.min_pident, min_qcov=args.min_qcov_hsp,
        min_length=args.min_length, max_evalue=args.max_evalue,
        use_tsd=not args.no_tsd, tsd_k=args.tsd_k, tsd_slop=args.tsd_slop,
        merge_slop=args.merge_slop)
    if not args.no_internal_filter:
        # Internal-region mask as a packed-bit genomic coverage track: bounds RAM
        # to ~genome/8 bytes regardless of the (multi-TB on wheat) internal TSV
        # size, vs the old load_blast() which held every HSP's btop in RAM and
        # OOM-killed on maize. Lossless: drop_near_internal_cov reproduces the old
        # prefilter_internal()+drop_near_internal() drops exactly.
        genome_lengths = {c: len(s) for c, s in genome.items()}
        int_cov = core.load_internal_coverage(
            int_tsv, genome_lengths,
            min_pident=args.int_min_pident, min_qcov=args.int_min_qcov_hsp,
            min_length=args.int_min_length, max_evalue=args.int_max_evalue)
        cands = core.drop_near_internal_cov(cands, int_cov, args.flank)

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

        nested_tsv = do_blast(args.consensus, internal_db,
                              min_pident=args.nested_min_pident,
                              min_qcov=args.min_qcov_hsp,
                              max_evalue=args.max_evalue)

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
        nested_hsps = core.load_internal_blast(
            nested_tsv, orient_map,
            min_pident=args.nested_min_pident, min_qcov=args.min_qcov_hsp,
            min_length=args.min_length, max_evalue=args.max_evalue)
        nested = core.make_candidates(
            nested_hsps, genome,
            min_pident=args.nested_min_pident, min_qcov=args.min_qcov_hsp,
            min_length=args.min_length, max_evalue=args.max_evalue,
            use_tsd=not args.no_tsd, tsd_k=args.tsd_k, tsd_slop=args.tsd_slop,
            merge_slop=args.merge_slop, source="nested")
        if not args.no_nested_guard:
            intervals = core.load_internal_intervals(args.internal)
            nested = core.drop_ltr_of_nested_intact(
                nested, intervals, args.nested_guard_flank)
        nb = core.write_tsv(nested, str(outdir / "nested_solo.tsv"),
                            name_prefix="nsolo")
        log(f"nested recovery: {nb} calls -> {outdir / 'nested_solo.tsv'}")
        cands = core.union_candidates(cands, nested, merge_slop=args.merge_slop)

    out_tsv = outdir / "solo_ltr.tsv"
    n = core.write_tsv(cands, str(out_tsv), name_prefix="solo")
    log(f"wrote {n} solo-LTR calls -> {out_tsv}", t0)

    if args.truth:
        solo_idx = core.load_truth_solo(args.truth)
        m = core.score(cands, solo_idx)
        log(f"vs truth: detected {m['detected']}/{m['n_solo']}  "
            f"TP={m['tp']} FP={m['fp']} FN={m['fn']}  "
            f"P={m['precision']:.4f} R={m['recall']:.4f} F1={m['f1']:.4f}")

    if not args.keep_intermediates:
        log("cleaning intermediates (use --keep-intermediates to retain)")
        cleanup_intermediates(masked, full_bed, blast_dir)


if __name__ == "__main__":
    main()
