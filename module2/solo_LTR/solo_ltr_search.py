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
             results with all the fields the downstream filter needs.
             pident/qcov/evalue cutoffs are pushed into blastn (kept looser
             than the downstream filter, so still lossless) to bound the
             on-disk TSV; the query is split into (threads // chunk-threads)
             length-balanced chunks run concurrently against the shared DB,
             scaling past blastn's internal threading ceiling on many-core
             hosts. Per-HSP filters are also re-applied POST-HOC in
             solo_ltr_core.py so one cached blast feeds many filter runs.

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
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
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
def blast_cache_path(out_dir, query_fa, db_fa, task, word_size, dust, evalue,
                     extra, perc_identity=None, qcov_hsp_perc=None):
    """Stable filename for a (query, db, blast-params) tuple."""
    key = json.dumps({
        "query": str(Path(query_fa).resolve()),
        "db": str(Path(db_fa).resolve()),
        "task": task, "word_size": word_size,
        "dust": dust, "evalue": evalue, "extra": extra,
        "perc_identity": perc_identity, "qcov_hsp_perc": qcov_hsp_perc,
    }, sort_keys=True)
    h = hashlib.sha1(key.encode()).hexdigest()[:12]
    qstem = Path(query_fa).stem
    return Path(out_dir) / f"{qstem}.{task}.ws{word_size}.dust{dust}.e{evalue}.{h}.tsv"


def build_blast_cmd(args, query, out_path, threads):
    """Assemble one blastn command. Shared by the single and chunked paths so
    the prefilter/format flags never drift between them.
    """
    outfmt = "6 " + " ".join(BLAST_COLS)
    cmd = ["blastn",
           "-query", str(query),
           "-db", args.db,
           "-task", args.task,
           "-word_size", str(args.word_size),
           "-dust", args.dust,
           "-evalue", str(args.evalue),
           "-outfmt", outfmt,
           "-num_threads", str(threads),
           "-mt_mode", str(args.mt_mode),
           "-max_target_seqs", str(args.max_target_seqs),
           "-soft_masking", "false",
           "-out", str(out_path)]
    # Blast-level prefilter: bound the on-disk TSV to ~survivors (disk analogue
    # of the downstream streaming filter). Omitted flag == no filter at this
    # threshold; kept looser than the downstream cutoff so it stays lossless.
    if args.perc_identity is not None:
        cmd += ["-perc_identity", str(args.perc_identity)]
    if args.qcov_hsp_perc is not None and args.qcov_hsp_perc > 0:
        cmd += ["-qcov_hsp_perc", str(args.qcov_hsp_perc)]
    # task-specific extras only added when meaningful
    if args.extra:
        cmd += args.extra.split()
    return cmd


def _iter_fasta(path):
    """Stream a FASTA, yielding (header_line, sequence) one record at a time.

    header_line keeps the leading '>' and everything after it (so the query's
    qseqid / '#'-classification survive verbatim); sequence is the joined,
    case-preserved residues. Holds a single record at a time -- never the whole
    library. Supports .gz.
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    header = None
    chunks = []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line.rstrip("\n")
                chunks = []
            else:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def chunk_query_by_length(query_fa, n_chunks, out_dir):
    """Bin-pack query records into <= n_chunks FASTA files balanced by total bp.

    Online greedy: each record is appended to the currently-lightest bin in a
    single streaming pass (one record held at a time -- no whole-library load),
    which keeps the parallel blastn processes finishing at about the same time
    instead of one straggler holding up the rest. Empty bins (when records <
    n_chunks) are dropped. Returns the list of non-empty chunk paths.
    """
    paths = [os.path.join(out_dir, f"chunk_{i}.fa") for i in range(n_chunks)]
    handles = [open(p, "w") for p in paths]
    loads = [0] * n_chunks
    try:
        for header, seq in _iter_fasta(query_fa):
            i = loads.index(min(loads))
            handles[i].write(f"{header}\n{seq}\n")
            loads[i] += len(seq)
    finally:
        for h in handles:
            h.close()
    return [p for p in paths if os.path.getsize(p) > 0]


def _run_blast_chunks(args, chunks, scratch, tmp_path, n_workers):
    """Run blastn over the query chunks through a dynamic pool, concat -> tmp_path.

    A bounded pool of n_workers concurrent blastn processes (chunk-threads each)
    pulls from a shared queue of chunks: the moment a chunk finishes, its worker
    grabs the next pending chunk, so every worker stays busy until the queue
    drains -- no core idles waiting on a straggler. With more chunks than workers
    (see --chunk-oversub) the unavoidable tail shrinks to a single small chunk
    instead of the slowest of n_workers fat ones. blastn runtime tracks hit count,
    not bp, so it is this dynamic pull -- not the bp-balanced split -- that
    actually balances load; chunks are merely submitted largest-first so the big
    ones start early (longest-processing-time scheduling) and the last thing
    running is small.

    Results are identical to a single blast (same query, DB, params -- only HSP
    line order differs, which the order-independent downstream filter ignores).
    """
    # Largest-first submission (LPT): start the slow chunks before the fast ones.
    chunks = sorted(chunks, key=os.path.getsize, reverse=True)
    log(f"blastn pool: {len(chunks)} chunks across {n_workers} workers "
        f"({args.chunk_threads} thread(s) each) -> {tmp_path.with_suffix('')}")

    def run_one(idx, chunk):
        cout = Path(scratch) / f"part_{idx}.tsv"
        cerr_path = Path(scratch) / f"part_{idx}.err"
        cmd = build_blast_cmd(args, chunk, cout, args.chunk_threads)
        with open(cerr_path, "w") as cerr:
            rc = subprocess.call([str(x) for x in cmd], stderr=cerr)
        if rc != 0:
            raise SystemExit(
                f"blast chunk failed (rc={rc}): {cmd}\n"
                f"{cerr_path.read_text()[-2000:]}")
        return cout

    parts = [None] * len(chunks)
    ex = ThreadPoolExecutor(max_workers=n_workers)
    try:
        futs = {ex.submit(run_one, i, c): i for i, c in enumerate(chunks)}
        for fut in as_completed(futs):
            parts[futs[fut]] = fut.result()  # re-raises a failed chunk's error
    except BaseException:
        ex.shutdown(wait=False, cancel_futures=True)  # drop the un-started chunks
        raise
    ex.shutdown(wait=True)

    with open(tmp_path, "wb") as out:
        for cout in parts:
            if cout is not None and os.path.exists(cout):
                with open(cout, "rb") as src:
                    shutil.copyfileobj(src, out)


def cmd_blast(args):
    t0 = time.time()
    out_path = blast_cache_path(args.out_dir, args.query, args.db,
                                args.task, args.word_size, args.dust,
                                args.evalue, args.extra,
                                args.perc_identity, args.qcov_hsp_perc)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if out_path.exists() and not args.force:
        log(f"blast result cached: {out_path}")
        print(out_path)
        return

    # Write to a .tmp then atomic-rename: a partial file (interrupted run,
    # or a reader checking mid-run) never looks like a complete cache entry.
    tmp_path = out_path.with_suffix(out_path.suffix + ".tmp")
    # One worker per (chunk-threads) slice of the thread budget; each worker is a
    # blastn process. Over-decompose the query into chunk-oversub x more chunks
    # than workers so a finished worker always has queued work to pull -- the
    # straggler tail collapses to one small chunk (see _run_blast_chunks).
    n_workers = max(1, args.threads // max(1, args.chunk_threads))

    if n_workers <= 1:
        log(f"blastn -> {out_path}")
        run(build_blast_cmd(args, args.query, tmp_path, args.threads))
        os.replace(tmp_path, out_path)
    else:
        n_chunks = n_workers * max(1, args.chunk_oversub)
        scratch = tempfile.mkdtemp(prefix="blastchunk_", dir=str(out_path.parent))
        try:
            chunks = chunk_query_by_length(args.query, n_chunks, scratch)
            if len(chunks) <= 1:  # query too small to split usefully
                log(f"blastn (query < 2 chunks, not chunked) -> {out_path}")
                run(build_blast_cmd(args, args.query, tmp_path, args.threads))
            else:
                _run_blast_chunks(args, chunks, scratch, tmp_path, n_workers)
            os.replace(tmp_path, out_path)
        finally:
            shutil.rmtree(scratch, ignore_errors=True)

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
    pb.add_argument("--perc-identity", type=float, default=None,
                    help="blastn -perc_identity prefilter (HSP %%-identity floor). "
                         "Set strictly below the downstream pident cutoff so the "
                         "Python filter stays authoritative. Default: no filter.")
    pb.add_argument("--qcov-hsp-perc", type=float, default=None,
                    help="blastn -qcov_hsp_perc prefilter (per-HSP query-coverage "
                         "floor). 0 or None disables it. Default: no filter.")
    pb.add_argument("--extra", default="",
                    help="Extra blastn flags as one string (e.g. for dc-megablast templates)")
    pb.add_argument("--threads", type=int, default=8)
    pb.add_argument("--chunk-threads", type=int, default=4,
                    help="Threads per blastn process when query-chunking. Sets "
                         "the worker count to (threads // chunk-threads): that "
                         "many blastn processes run concurrently against the "
                         "shared DB, scaling past blastn's internal threading "
                         "ceiling on many-core hosts. Set >= --threads (or a "
                         "1-record query) to run a single blast.")
    pb.add_argument("--chunk-oversub", type=int, default=16,
                    help="Over-decomposition factor: split the query into "
                         "(threads // chunk-threads) * chunk-oversub length-"
                         "balanced chunks fed through a pool of (threads // "
                         "chunk-threads) concurrent blastn workers. >1 lets a "
                         "finished worker pull queued work instead of idling on "
                         "a straggler, shrinking the tail to one small chunk. "
                         "1 = one chunk per worker (static; the step ends with "
                         "the single slowest chunk). Default 16.")
    pb.add_argument("--mt-mode", type=int, default=1,
                    help="1 = split by queries (best for many small queries vs big db)")
    pb.add_argument("--max-target-seqs", type=int, default=5000)
    pb.add_argument("--force", action="store_true")
    pb.set_defaults(func=cmd_blast)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
