#!/usr/bin/env python3
import argparse
import os
import sys
import subprocess
import shutil
import tempfile
from collections import defaultdict, namedtuple
from concurrent.futures import ThreadPoolExecutor, as_completed

# Run this after module2.py to remove LTR-RT candidates with too much tandem sequence. 
# python tandem_cleaner.py -fa [module2_in]_Kmer2LTR_TSD_class_purge.fa -d [module2_in]_Kmer2LTR_TSD_class_purge -t 100 -lp 20 -wp 30 -o [module2_in]_Kmer2LTR_TSD_class_purge.tandem

TRF_REPO = "https://github.com/lh3/TRF-mod"
TRF_DIR = "TRF-mod"
TRF_BIN = os.path.join(TRF_DIR, "trf-mod")

SeqResult = namedtuple(
    "SeqResult",
    [
        "seq_id",
        "length",
        "total_tandem",
        "whole_percent",
        "ltr_total",
        "ltr_tandem",
        "ltr_percent",
        "keep",
        "reason",
    ],
)


def run_cmd(cmd, cwd=None):
    """Run a shell command and raise on error."""
    try:
        subprocess.run(cmd, cwd=cwd, check=True)
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"[ERROR] Command failed: {' '.join(cmd)}\n")
        sys.stderr.write(f"        {e}\n")
        sys.exit(1)


def ensure_trf_mod():
    """Clone and compile TRF-mod if not already available."""
    if os.path.isfile(TRF_BIN) and os.access(TRF_BIN, os.X_OK):
        return

    if not os.path.isdir(TRF_DIR):
        sys.stderr.write(f"[INFO] Cloning TRF-mod from {TRF_REPO}\n")
        run_cmd(["git", "clone", TRF_REPO])

    sys.stderr.write("[INFO] Compiling TRF-mod\n")
    run_cmd(["make", "-f", "compile.mak"], cwd=TRF_DIR)

    if not (os.path.isfile(TRF_BIN) and os.access(TRF_BIN, os.X_OK)):
        sys.stderr.write("[ERROR] TRF-mod binary not found after compilation.\n")
        sys.exit(1)


def parse_fasta(fasta_path):
    """Parse FASTA into ordered list and dict of id->sequence."""
    seqs = {}
    order = []
    cur_id = None
    cur_seq_chunks = []

    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if cur_id is not None:
                    seqs[cur_id] = "".join(cur_seq_chunks)
                    order.append(cur_id)
                cur_id = line[1:].split()[0]  # ID up to first whitespace
                cur_seq_chunks = []
            else:
                cur_seq_chunks.append(line)
        if cur_id is not None:
            seqs[cur_id] = "".join(cur_seq_chunks)
            order.append(cur_id)

    return order, seqs


def write_fasta(path, ids, seqs):
    with open(path, "w") as out:
        for sid in ids:
            seq = seqs[sid]
            out.write(f">{sid}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + "\n")


def run_trf_on_fasta(in_fa, out_bed):
    """Run TRF-mod on a single FASTA file."""
    sys.stderr.write(f"[INFO] Running TRF-mod on {in_fa}\n")
    with open(out_bed, "w") as out_fh:
        cmd = [
            TRF_BIN,
            in_fa,
            "-a2",
            "-b7",
            "-g7",
            "-A80",
            "-G10",
            "-s300", # EDTA uses 1000.
            "-p2000",
        ]
        subprocess.run(cmd, stdout=out_fh, check=True)


def split_fasta(order, seqs, threads, tmpdir, base_prefix):
    """Split FASTA into N chunks by sequence count."""
    n = max(1, min(threads, len(order)))
    chunks = [[] for _ in range(n)]
    for i, sid in enumerate(order):
        chunks[i % n].append(sid)

    chunk_files = []
    for i, chunk_ids in enumerate(chunks):
        if not chunk_ids:
            continue
        fname = os.path.join(tmpdir, f"{base_prefix}.chunk{i+1}.fa")
        write_fasta(fname, chunk_ids, seqs)
        chunk_files.append(fname)
    return chunk_files


def run_trf_parallel(in_fa, out_bed, threads, order, seqs):
    """
    Run TRF-mod either single-threaded or by splitting the FASTA into chunks
    and processing them in parallel, then concatenating the BED outputs.
    """
    if threads <= 1:
        run_trf_on_fasta(in_fa, out_bed)
        return

    sys.stderr.write(
        f"[INFO] Running TRF-mod in parallel with {threads} threads (by splitting FASTA)\n"
    )

    tmpdir = tempfile.mkdtemp(prefix="trfmod_tmp_")
    try:
        chunk_files = split_fasta(order, seqs, threads, tmpdir, "candidate_LTR")
        if not chunk_files:
            # no sequences
            open(out_bed, "w").close()
            return

        bed_files = [cf + ".bed" for cf in chunk_files]

        def worker(in_path, out_path):
            run_trf_on_fasta(in_path, out_path)
            return out_path

        futures = []
        with ThreadPoolExecutor(max_workers=len(chunk_files)) as ex:
            for cf, bf in zip(chunk_files, bed_files):
                futures.append(ex.submit(worker, cf, bf))
            for f in as_completed(futures):
                _ = f.result()

        # Concatenate
        with open(out_bed, "w") as out:
            for bf in bed_files:
                if os.path.exists(bf):
                    with open(bf) as in_fh:
                        shutil.copyfileobj(in_fh, out)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def parse_trf_bed(bed_path):
    """
    Parse TRF-mod BED-like output.
    We only use columns 1-3: seq_id, start, end (0-based, end-exclusive).
    Return dict: seq_id -> list[(start, end)].
    """
    repeats = defaultdict(list)
    if not os.path.exists(bed_path):
        return repeats

    with open(bed_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            sid = parts[0]
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue
            if end > start:
                repeats[sid].append((start, end))

    return repeats


def merge_intervals(intervals):
    """Merge potentially overlapping intervals and return merged list."""
    if not intervals:
        return []
    intervals = sorted(intervals, key=lambda x: x[0])
    merged = [list(intervals[0])]
    for s, e in intervals[1:]:
        last_s, last_e = merged[-1]
        if s <= last_e:
            merged[-1][1] = max(last_e, e)
        else:
            merged.append([s, e])
    return [(s, e) for s, e in merged]


def interval_overlap_len(a_start, a_end, b_start, b_end):
    """Length of overlap between [a_start, a_end) and [b_start, b_end)."""
    s = max(a_start, b_start)
    e = min(a_end, b_end)
    return max(0, e - s)


def parse_domain_file(domain_path):
    """
    Parse domain file; return:
      domain_info: dict[id] -> ltr_len (int, from column 2)
      lines: list of raw lines
    """
    domain_info = {}
    lines = []
    with open(domain_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            lines.append(line)
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            sid = parts[0]
            try:
                ltr_len = int(parts[1])
            except ValueError:
                continue
            domain_info[sid] = ltr_len
    return domain_info, lines


def filter_sequences(
    order,
    seqs,
    repeats,
    domain_info,
    lp,
    wp,
):
    """
    Compute tandem repeat coverage for each sequence and decide keep/exclude.

    lp: max allowed percent of LTR bases that can be tandem (None to ignore)
    wp: max allowed percent of whole sequence that can be tandem (None to ignore)
    """
    results = {}
    for sid in order:
        seq = seqs[sid]
        slen = len(seq)
        intervals = merge_intervals(repeats.get(sid, []))

        total_tandem = sum(e - s for s, e in intervals)
        whole_percent = (100.0 * total_tandem / slen) if slen > 0 else 0.0

        ltr_total = 0
        ltr_tandem = 0
        ltr_percent = None

        # LTR-based calculation only if lp is set AND domain info exists for this ID
        if lp is not None and sid in domain_info:
            ltr_len = domain_info[sid]
            # Define 5' and 3' LTR intervals (clipped to [0, slen))
            ltr5_start, ltr5_end = 0, min(ltr_len, slen)
            ltr3_start, ltr3_end = max(slen - ltr_len, 0), slen

            ltr_intervals = [(ltr5_start, ltr5_end), (ltr3_start, ltr3_end)]
            # Merge LTR intervals in case they overlap for very short sequences
            ltr_intervals = merge_intervals(ltr_intervals)

            ltr_total = sum(e - s for s, e in ltr_intervals)
            if ltr_total > 0:
                for ls, le in ltr_intervals:
                    for rs, re in intervals:
                        ltr_tandem += interval_overlap_len(ls, le, rs, re)
                ltr_percent = 100.0 * ltr_tandem / ltr_total
            else:
                ltr_percent = 0.0

        keep = True
        reasons = []

        if wp is not None and whole_percent > wp:
            keep = False
            reasons.append(f"wp>{wp}")

        if lp is not None:
            if sid not in domain_info:
                # Can't evaluate LTR-based filter without domain info.
                # We keep the sequence but record the situation.
                reasons.append("lp_not_evaluated_no_domain")
            else:
                if ltr_percent is not None and ltr_percent > lp:
                    keep = False
                    reasons.append(f"lp>{lp}")

        reason_str = ",".join(reasons) if reasons else "ok"
        results[sid] = SeqResult(
            seq_id=sid,
            length=slen,
            total_tandem=total_tandem,
            whole_percent=whole_percent,
            ltr_total=ltr_total,
            ltr_tandem=ltr_tandem,
            ltr_percent=ltr_percent if ltr_percent is not None else -1,
            keep=keep,
            reason=reason_str,
        )

    return results


def write_filtered_outputs(
    out_prefix,
    order,
    seqs,
    results,
    domain_info,
    domain_lines,
):
    # 1) Filtered FASTA
    fa_out = out_prefix + ".fa"
    keep_ids = [sid for sid in order if results[sid].keep]
    write_fasta(fa_out, keep_ids, seqs)
    sys.stderr.write(f"[INFO] Wrote filtered FASTA: {fa_out}\n")

    # 2) Filtered domain TSV (if domain info provided)
    if domain_lines is not None:
        tsv_out = out_prefix + ".tsv"
        with open(tsv_out, "w") as out:
            for line in domain_lines:
                if line.startswith("#") or not line.strip():
                    out.write(line + "\n")
                    continue
                parts = line.split("\t")
                sid = parts[0]
                # If we have a decision for this ID, apply it; otherwise keep as-is
                if sid in results and not results[sid].keep:
                    continue
                out.write(line + "\n")
        sys.stderr.write(f"[INFO] Wrote filtered domain file: {tsv_out}\n")

    # 3) Report file
    rep_out = out_prefix + ".report"
    with open(rep_out, "w") as out:
        out.write(
            "seq_id\tlength\ttotal_tandem\twhole_percent\t"
            "ltr_total\tltr_tandem\tltr_percent\tkeep\treason\n"
        )
        for sid in order:
            r = results[sid]
            ltr_percent_str = (
                f"{r.ltr_percent:.3f}" if r.ltr_percent >= 0 else "NA"
            )
            out.write(
                f"{r.seq_id}\t{r.length}\t{r.total_tandem}\t{r.whole_percent:.3f}\t"
                f"{r.ltr_total}\t{r.ltr_tandem}\t{ltr_percent_str}\t"
                f"{'KEEP' if r.keep else 'REMOVE'}\t{r.reason}\n"
            )
    sys.stderr.write(f"[INFO] Wrote report: {rep_out}\n")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Run TRF-mod on candidate LTR-RTs and filter sequences based on "
            "tandem repeat content in LTRs and/or whole sequence."
        )
    )
    parser.add_argument(
        "-fa",
        "--fasta",
        required=True,
        help="Multi-seq FASTA of candidate LTR-RTs",
    )
    parser.add_argument(
        "-d",
        "--domain",
        help="Domain file (TSV). Required to use -lp.",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        help="Number of threads for TRF-mod (by splitting FASTA). Default: 1",
    )
    parser.add_argument(
        "-lp",
        "--ltr_percent",
        type=float,
        default=None,
        help=(
            "Max allowed percent of LTR bases that can be tandem. "
            "Requires domain file. Example: 30 means up to 30%%."
        ),
    )
    parser.add_argument(
        "-wp",
        "--whole_percent",
        type=float,
        default=None,
        help=(
            "Max allowed percent of whole sequence bases that can be tandem. "
            "Example: 40 means up to 40%%."
        ),
    )
    parser.add_argument(
        "-o",
        "--out_prefix",
        required=True,
        help="Output prefix for .fa, .tsv, and .report",
    )

    args = parser.parse_args()

    # Sanity checks
    if args.ltr_percent is not None and not args.domain:
        sys.stderr.write(
            "[ERROR] -lp/--ltr_percent was specified but no domain file (-d) was provided.\n"
        )
        sys.exit(1)

    ensure_trf_mod()

    fasta_path = args.fasta
    out_prefix = args.out_prefix
    threads = max(1, args.threads)
    lp = args.ltr_percent
    wp = args.whole_percent

    sys.stderr.write(f"[INFO] Reading FASTA: {fasta_path}\n")
    order, seqs = parse_fasta(fasta_path)
    if not order:
        sys.stderr.write("[ERROR] No sequences found in FASTA.\n")
        sys.exit(1)

    # Run TRF-mod
    trf_bed = fasta_path + ".trf.bed"
    run_trf_parallel(fasta_path, trf_bed, threads, order, seqs)

    # Parse TRF output
    sys.stderr.write(f"[INFO] Parsing TRF output: {trf_bed}\n")
    repeats = parse_trf_bed(trf_bed)

    # Parse domain file if provided
    domain_info = {}
    domain_lines = None
    if args.domain:
        sys.stderr.write(f"[INFO] Reading domain file: {args.domain}\n")
        domain_info, domain_lines = parse_domain_file(args.domain)

    # Compute coverage and filter
    sys.stderr.write("[INFO] Computing tandem repeat coverage & applying filters\n")
    results = filter_sequences(order, seqs, repeats, domain_info, lp, wp)

    # Write outputs
    write_filtered_outputs(out_prefix, order, seqs, results, domain_info, domain_lines)


if __name__ == "__main__":
    main()
