#!/usr/bin/env python3
"""solo_ltr_search.py - cached BLAST searches for the solo-LTR pipeline.

Three operations, each cached on input + parameter hash so the slow steps
are not repeated across a parameter sweep:

  1. mask :  parse internal-FASTA headers (chr:start-end#...) for the
             full-length LTR-RT coordinates, write a BED, and hardmask
             the genome with bedtools.

  2. db   :  makeblastdb on the masked genome.

  3. blast:  blastn -query <consensus|internal> -db <masked db>
             with a defined task/word_size/dust/evalue. We dump tabular
             results with all the fields the downstream filter needs;
             per-HSP percent-identity and qcov filters are applied
             POST-HOC in solo_ltr_filter.py so one blast run feeds many
             filter combinations.

Heng Li-style: streaming where it matters, no shell-quoting tricks, and
nothing depends on absolute paths.
"""
import argparse
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

# Custom outfmt-6 columns: enough for full filtering + divergence downstream.
BLAST_COLS = (
    "qseqid sseqid pident length qlen slen "
    "qstart qend sstart send evalue bitscore qcovs qcovhsp mismatch gapopen "
    "btop"
).split()

# Standard BED is 0-based half-open. The chr:start-end in our FASTA
# headers comes from BED-derived extractions; we keep that convention.
HEADER_RE = re.compile(r"^([^:\s]+):(\d+)-(\d+)")


def log(msg, t0=None):
    """Minimal, timestamped stderr log line."""
    if t0 is not None:
        sys.stderr.write(f"[solo_ltr_search] {msg} ({time.time()-t0:.1f}s)\n")
    else:
        sys.stderr.write(f"[solo_ltr_search] {msg}\n")
    sys.stderr.flush()


def run(cmd, env=None):
    """Run a subprocess; raise on non-zero. stderr passes through."""
    sys.stderr.write("$ " + " ".join(str(x) for x in cmd) + "\n")
    sys.stderr.flush()
    rc = subprocess.call([str(x) for x in cmd], env=env)
    if rc != 0:
        raise SystemExit(f"command failed (rc={rc}): {cmd}")


# --------------------------------------------------------------------------- #
# Step 1 -- hardmask full-length LTR-RTs out of the genome
# --------------------------------------------------------------------------- #
def headers_to_bed(internal_fa, bed_out):
    """Stream the internal-FASTA, emit a BED of full-length LTR-RT loci.

    The internal-FASTA header is the FULL-LENGTH LTR-RT coordinate (5' LTR
    through 3' LTR), not the internal region. Hardmasking those intervals
    removes the loci that already have a high-confidence intact assignment,
    leaving everywhere else for the consensus search.
    """
    n = 0
    out = open(bed_out, "w")
    opener = gzip.open if str(internal_fa).endswith(".gz") else open
    with opener(internal_fa, "rt") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            m = HEADER_RE.match(line[1:])
            if not m:
                continue
            chrom, s, e = m.group(1), int(m.group(2)), int(m.group(3))
            if e < s:
                s, e = e, s
            out.write(f"{chrom}\t{s}\t{e}\n")
            n += 1
    out.close()
    return n


def cmd_mask(args):
    """Run step 1: write hardmask BED, then bedtools maskfasta -> N's."""
    t0 = time.time()
    log(f"parsing headers from {args.internal_fa}")
    n = headers_to_bed(args.internal_fa, args.bed_out)
    log(f"wrote {n} regions to {args.bed_out}", t0)

    if args.fasta_out is None:
        return

    if Path(args.fasta_out).exists() and not args.force:
        log(f"masked fasta exists: {args.fasta_out} (use --force to redo)")
        return

    t1 = time.time()
    log(f"hardmasking {args.genome} -> {args.fasta_out}")
    # -mc N keeps lower-case off; bedtools default is N for hard mask.
    run(["bedtools", "maskfasta",
         "-fi", args.genome, "-bed", args.bed_out,
         "-fo", args.fasta_out, "-mc", "N"])
    log("hardmask done", t1)


# --------------------------------------------------------------------------- #
# Step 2 -- makeblastdb on the masked genome
# --------------------------------------------------------------------------- #
def cmd_db(args):
    t0 = time.time()
    out_prefix = args.out if getattr(args, "out", None) else args.fasta
    expect = Path(str(out_prefix) + ".nhr")
    if expect.exists() and not args.force:
        log(f"blast db exists for {out_prefix} (use --force to redo)")
        return
    cmd = ["makeblastdb", "-in", args.fasta, "-dbtype", "nucl"]
    if not getattr(args, "no_parse_seqids", False):
        cmd.append("-parse_seqids")
    if getattr(args, "out", None):
        cmd += ["-out", str(args.out)]
    log(f"makeblastdb on {args.fasta} -> {out_prefix}")
    run(cmd)
    log("makeblastdb done", t0)


# --------------------------------------------------------------------------- #
# Step 3 -- blastn with permissive cutoffs (filter post-hoc)
# --------------------------------------------------------------------------- #
def blast_cache_path(out_dir, query_fa, db_fa, task, word_size, dust, evalue, extra):
    """Stable filename for a (query, db, blast-params) tuple."""
    key = json.dumps({
        "query": str(Path(query_fa).resolve()),
        "db": str(Path(db_fa).resolve()),
        "task": task, "word_size": word_size,
        "dust": dust, "evalue": evalue, "extra": extra,
    }, sort_keys=True)
    h = hashlib.sha1(key.encode()).hexdigest()[:12]
    qstem = Path(query_fa).stem
    return Path(out_dir) / f"{qstem}.{task}.ws{word_size}.dust{dust}.e{evalue}.{h}.tsv"


def cmd_blast(args):
    t0 = time.time()
    out_path = blast_cache_path(args.out_dir, args.query, args.db,
                                args.task, args.word_size, args.dust,
                                args.evalue, args.extra)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if out_path.exists() and not args.force:
        log(f"blast result cached: {out_path}")
        print(out_path)
        return

    outfmt = "6 " + " ".join(BLAST_COLS)
    # Write to a .tmp then atomic-rename: a partial file (interrupted run,
    # or a reader checking mid-run) never looks like a complete cache entry.
    tmp_path = out_path.with_suffix(out_path.suffix + ".tmp")
    cmd = ["blastn",
           "-query", args.query,
           "-db", args.db,
           "-task", args.task,
           "-word_size", str(args.word_size),
           "-dust", args.dust,
           "-evalue", str(args.evalue),
           "-outfmt", outfmt,
           "-num_threads", str(args.threads),
           "-mt_mode", str(args.mt_mode),
           "-max_target_seqs", str(args.max_target_seqs),
           "-soft_masking", "false",
           "-out", str(tmp_path)]
    # task-specific extras only added when meaningful
    if args.extra:
        for tok in args.extra.split():
            cmd.append(tok)

    log(f"blastn -> {out_path}")
    run(cmd)
    os.replace(tmp_path, out_path)
    log(f"blast done; {out_path}", t0)
    print(out_path)


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def main():
    p = argparse.ArgumentParser(
        description="Cached BLAST steps for the solo-LTR detection pipeline.")
    sub = p.add_subparsers(dest="cmd", required=True)

    pm = sub.add_parser("mask", help="hardmask full-length LTR-RT loci")
    pm.add_argument("--internal-fa", required=True,
                    help="LTRs.alns.internal.fa (headers carry full-length coords)")
    pm.add_argument("--genome", required=True,
                    help="Reference genome FASTA")
    pm.add_argument("--bed-out", required=True,
                    help="BED of full-length LTR-RT regions (output)")
    pm.add_argument("--fasta-out", default=None,
                    help="Hardmasked genome FASTA (output, optional)")
    pm.add_argument("--force", action="store_true")
    pm.set_defaults(func=cmd_mask)

    pd = sub.add_parser("db", help="makeblastdb on a FASTA")
    pd.add_argument("--fasta", required=True)
    pd.add_argument("--out", default=None,
                    help="DB prefix (default: alongside --fasta)")
    pd.add_argument("--no-parse-seqids", action="store_true",
                    help="Build DB without -parse_seqids (needed when sequence "
                         "IDs exceed BLAST's 50-char local-id limit)")
    pd.add_argument("--force", action="store_true")
    pd.set_defaults(func=cmd_db)

    pb = sub.add_parser("blast", help="run blastn with permissive defaults")
    pb.add_argument("--query", required=True)
    pb.add_argument("--db", required=True,
                    help="BLAST db prefix (= the FASTA path used in makeblastdb)")
    pb.add_argument("--out-dir", required=True,
                    help="Directory for cached blast tabular outputs")
    pb.add_argument("--task", default="blastn",
                    choices=("blastn", "dc-megablast", "megablast", "blastn-short"))
    pb.add_argument("--word-size", type=int, default=11)
    pb.add_argument("--dust", default="no",
                    help="DUST filter: 'no' or 'yes' or 'level window linker'")
    pb.add_argument("--evalue", default="10",
                    help="Permissive default; filter tighter post-hoc")
    pb.add_argument("--extra", default="",
                    help="Extra blastn flags as one string (e.g. for dc-megablast templates)")
    pb.add_argument("--threads", type=int, default=8)
    pb.add_argument("--mt-mode", type=int, default=1,
                    help="1 = split by queries (best for many small queries vs big db)")
    pb.add_argument("--max-target-seqs", type=int, default=5000)
    pb.add_argument("--force", action="store_true")
    pb.set_defaults(func=cmd_blast)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
