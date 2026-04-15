#!/usr/bin/env python3
"""
(1) Masks genes and non-LTR TEs.
(2) Runs ltrharvest and ltr_finder and merges results.
(3) Identify those with TSD, those with LTR-RT protein, and those with homology to those with either TSD or protein. 
(4) Runs kmer2ltr on the short-list of candidates. 
(5) Purges duplicate LTR-RTs based on kmer2ltr LTR divergence.

# I could consider adding DeepTE or CREATE as an alternative to TEsorter although I need to check their speed and efficiency feasiblility. 
# This script differs from 'ltrharvest3.py' in that it performs TSD searching before TEsorter, then feeds those TSD-containing LTR-RTs into TEsorter for detection durring 2-pass. thus improving annotation of unknown TEs. 

# The bash below runs this script:
bash synLTR/module2/ltrharvest_wrapper.sh --genome burnin.fasta --proteins ../PrinTE/data/TAIR10.pep.fa.gz --threads 20 --out_prefix burnin_ltr
# It runs in 5:17.57 minutes with 20 threads.
python ../PrinTE/util/bedtools.py -pass_scn <( awk 'BEGIN{OFS="\t"} { sub(/#.*/, "", $1); print }' burnin_ltr_r1_kmer2ltr_dedup) -bed <(grep LTR  burnin.bed | grep -v _FRAG) -r 0.9
   Overlapping entries: 4328 (4328 unique) = 4328/4468 = 0.96866 = true positives.
   Entries unique to SCN/PASS file: 170 = 170/4468 = 0.038 = false positive.
   Entries unique to BED file: 140 = 140/4468 = 0.03133 = false negative.
grep LTR  burnin.bed | grep -v _FRAG | wc -l
   4468
   F1 = 0.966

# Compare to LTR_retriever gold-standard.
perl EDTA/EDTA_raw.pl --genome burnin.fasta --type ltr --threads 20
# It runs in 20:28.15 minutes with 20 threads.
python ../PrinTE/util/bedtools.py -pass_scn burnin.fasta.mod.EDTA.raw/LTR/burnin.fasta.mod.pass.list -bed <(grep LTR  burnin.bed | grep -v _FRAG) -r 0.9
   Overlapping entries: 2085 (2085 unique) = 2085/4468 = 0.46665 = true positives.
   Entries unique to SCN/PASS file: 0 = 0/4468 = 0.0 = false positive.
   Entries unique to BED file: 2383 = 2383/4468 = 0.533348 = false negative.
   F1 = 0.64
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Set, Tuple, Iterable, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed


# -----------------------------
# utilities
# -----------------------------

def run(cmd: List[str], cwd: Optional[str] = None, check: bool = False, capture: bool = True, text: bool = True, verbose: bool = False):
    if verbose:
        label = f"  $ {' '.join(cmd)}"
        if cwd:
            label += f"\n    (cwd: {cwd})"
        print(label, file=sys.stderr)
    if capture:
        r = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=text)
    else:
        r = subprocess.run(cmd, cwd=cwd)
    if check and r.returncode != 0:
        raise RuntimeError(f"Command failed ({r.returncode}): {' '.join(cmd)}\nSTDERR:\n{(r.stderr or '').strip()}")
    return r

def which_or_die(prog: str):
    p = shutil.which(prog)
    if not p:
        raise RuntimeError(f"Required tool not found on PATH: {prog}")
    return p

def mkdirp(p: Path):
    p.mkdir(parents=True, exist_ok=True)
    return p

def _clean_path(path, verbose: bool = False):
    """Remove a file or directory quietly (used for incremental --clean purges)."""
    p = Path(path)
    if p.is_dir():
        shutil.rmtree(p, ignore_errors=True)
        if verbose:
            print(f"[CLEAN] removed dir:  {p}", file=sys.stderr)
    elif p.exists():
        try:
            p.unlink()
        except Exception:
            pass
        if verbose:
            print(f"[CLEAN] removed file: {p}", file=sys.stderr)

def _format_hms(seconds: float) -> str:
    if not seconds or seconds == float("inf"):
        return "unknown"
    seconds = int(seconds)
    h = seconds // 3600
    m = (seconds % 3600) // 60
    s = seconds % 60
    if h > 0:
        return f"{h:d}:{m:02d}:{s:02d}"
    else:
        return f"{m:d}:{s:02d}"

# -----------------------------
# FASTA helpers
# -----------------------------

def iter_fasta(path: str) -> Iterable[Tuple[str, str]]:
    name = None
    seq_chunks: List[str] = []
    with open(path, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq_chunks)
                name = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if name is not None:
            yield name, "".join(seq_chunks)

def write_fasta(records: Iterable[Tuple[str, str]], out_path: str, wrap: int = 60):
    with open(out_path, "w") as out:
        for name, seq in records:
            out.write(f">{name}\n")
            for i in range(0, len(seq), wrap):
                out.write(seq[i:i+wrap] + "\n")

def load_fasta_as_dict(path: str) -> Dict[str, str]:
    """Loads entire FASTA into memory. For very large genomes, consider using faidx-based fetch instead."""
    d: Dict[str, str] = {}
    for name, seq in iter_fasta(path):
        d[name] = seq
    return d

def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(comp)[::-1]

def merge_pass2_fastas(user_fa: Optional[str], tsd_fa: Optional[str], out_fa: str) -> str:
    """
    Merge two pass2-classified FASTAs (user provided and TSD-derived),
    deduplicating by header name token (up to whitespace, without '>').
    Returns out_fa path. If both inputs are missing/empty, returns "".
    """
    inputs = []
    if user_fa:
        p = Path(user_fa)
        if p.exists() and p.stat().st_size > 0:
            inputs.append(str(p))
    if tsd_fa:
        p = Path(tsd_fa)
        if p.exists() and p.stat().st_size > 0:
            inputs.append(str(p))

    if not inputs:
        return ""

    seen = set()
    n_written = 0
    with open(out_fa, "w") as out:
        for fa in inputs:
            for name, seq in iter_fasta(fa):
                if name in seen:
                    continue
                seen.add(name)
                out.write(f">{name}\n")
                for i in range(0, len(seq), 60):
                    out.write(seq[i:i+60] + "\n")
                n_written += 1

    if n_written == 0:
        Path(out_fa).touch()
        return ""

    return out_fa


# -----------------------------
# Step 2: miniprot gene masking
# -----------------------------

@dataclass
class Interval:
    start: int  # 0-based inclusive
    end: int    # 0-based exclusive

def parse_feature_intervals_from_gff_text(gff_text: str, feature_types: set) -> Dict[str, List[Interval]]:
    """
    Extract intervals for any requested feature types (e.g., {"mRNA"} or {"CDS"}).
    Coordinates in GFF are 1-based inclusive; we convert to 0-based half-open.
    """
    intervals: Dict[str, List[Interval]] = {}
    for line in gff_text.splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 9:
            continue
        seqid, ftype = parts[0], parts[2]
        if ftype not in feature_types:
            continue
        try:
            s = int(parts[3]) - 1
            e = int(parts[4])
        except ValueError:
            continue
        if e <= s:
            continue
        intervals.setdefault(seqid, []).append(Interval(s, e))
    return intervals

def merge_intervals(intervals: List[Interval]) -> List[Interval]:
    if not intervals:
        return []
    intervals = sorted(intervals, key=lambda x: (x.start, x.end))
    merged = [Interval(intervals[0].start, intervals[0].end)]
    for iv in intervals[1:]:
        cur = merged[-1]
        if iv.start <= cur.end:
            cur.end = max(cur.end, iv.end)
        else:
            merged.append(Interval(iv.start, iv.end))
    return merged

def merge_intervals_per_chrom(d: Dict[str, List[Interval]]) -> Dict[str, List[Interval]]:
    return {chrom: merge_intervals(ivs) for chrom, ivs in d.items()}

def hardmask_fasta_by_intervals(in_fasta: str, intervals: Dict[str, List[Interval]], out_fasta: str, out_bed: str):
    with open(out_bed, "w") as bed:
        for chrom, ivs in intervals.items():
            for iv in ivs:
                bed.write(f"{chrom}\t{iv.start}\t{iv.end}\n")

    out_records = []
    for name, seq in iter_fasta(in_fasta):
        if name not in intervals:
            out_records.append((name, seq))
            continue
        arr = list(seq)
        for iv in intervals[name]:
            s = max(0, iv.start)
            e = min(len(arr), iv.end)
            for i in range(s, e):
                arr[i] = "N"
        out_records.append((name, "".join(arr)))
    write_fasta(out_records, out_fasta)

def run_miniprot_gene_mask(in_fasta: str, protein_faa: str, out_prefix: str,
                           threads: int, miniprot_path: str,
                           outn: int = 1, outs: float = 0.8, outc: float = 0.5,
                           gene_mask: str = "mRNA", verbose: bool = False):

    gff_path = f"{out_prefix}.genic.gff"
    bed_path = f"{out_prefix}.genic.mask.bed"
    masked_fasta = f"{out_prefix}.genic_masked.fa"

    cmd = [
        miniprot_path,
        "--gff-only",
        "-t", str(threads),
        in_fasta,
        protein_faa,
        "-P", out_prefix,
        "--outn", str(outn),
        "--outs", str(outs),
        "--outc", str(outc),
    ]
    if verbose:
        print(f"  $ {' '.join(cmd)}", file=sys.stderr)
    r = run(cmd)
    if r.returncode != 0:
        raise RuntimeError(f"miniprot failed:\n{(r.stderr or '').strip()}")

    gff_text = r.stdout or ""
    Path(gff_path).write_text(gff_text)

    if gene_mask == "mRNA":
        feature_types = {"mRNA"}
    elif gene_mask == "CDS":
        feature_types = {"CDS"}
    else:
        feature_types = {"mRNA", "CDS"}

    intervals = merge_intervals_per_chrom(parse_feature_intervals_from_gff_text(gff_text, feature_types))

    if not intervals:
        shutil.copyfile(in_fasta, masked_fasta)
        Path(bed_path).write_text("")
        return gff_path, bed_path, masked_fasta

    hardmask_fasta_by_intervals(in_fasta, intervals, masked_fasta, bed_path)
    return gff_path, bed_path, masked_fasta


# -----------------------------
# Step 3: non-LTR masking via minimap2 PAF
# -----------------------------

def write_nonltr_library(te_fa: str, out_fa: str):
    kept = []
    for name, seq in iter_fasta(te_fa):
        # (6) case-insensitive match for "LTR"
        if re.search(r"LTR", name, flags=re.IGNORECASE):
            continue
        kept.append((name, seq))
    write_fasta(kept, out_fa)

def paf_to_bed_filtered(paf_path: str, out_bed: str, seq_ident: float, aln_len: int, qcov: float):
    hits = []
    with open(paf_path, "r") as f:
        for line in f:
            if not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 12:
                continue
            try:
                qname = p[0]
                qlen = float(p[1])
                qstart = float(p[2])
                qend = float(p[3])
                tname = p[5]
                tstart = int(p[7])
                tend = int(p[8])
                matches = float(p[9])
                alen_v = float(p[10])
            except ValueError:
                continue
            if qlen <= 0 or alen_v <= 0:
                continue
            ident = matches / alen_v
            qc = (qend - qstart) / qlen
            if ident >= seq_ident and alen_v >= aln_len and qc >= qcov:
                if tend < tstart:
                    tstart, tend = tend, tstart
                hits.append((tname, tstart, tend, qname, ident))

    hits.sort(key=lambda x: (x[0], x[1], x[2]))
    with open(out_bed, "w") as out:
        for tname, s, e, qname, ident in hits:
            out.write(f"{tname}\t{s}\t{e}\t{qname}\t{ident:.6f}\n")


# -----------------------------
# Step 4: chunk genome
# -----------------------------

@dataclass
class ChunkInfo:
    chunk_id: str
    chrom: str
    start0: int
    end0: int
    fasta_path: str

def make_chunks(in_fasta: str, out_dir: str, size: int, overlap: int) -> List[ChunkInfo]:
    outp = mkdirp(Path(out_dir))
    chunks: List[ChunkInfo] = []
    for chrom, seq in iter_fasta(in_fasta):
        L = len(seq)
        if L == 0:
            continue
        step = max(1, size - overlap)
        idx = 0
        for start0 in range(0, L, step):
            end0 = min(L, start0 + size)
            idx += 1
            header = f"{chrom}__chunk{idx}__{start0+1}-{end0}"
            chunk_id = f"{chrom}.chunk{idx}.{start0+1}-{end0}"
            chunk_fa = outp / f"{chunk_id}.fa"
            write_fasta([(header, seq[start0:end0])], str(chunk_fa))
            chunks.append(ChunkInfo(chunk_id=chunk_id, chrom=chrom, start0=start0, end0=end0, fasta_path=str(chunk_fa)))
            if end0 >= L:
                break
    return chunks

def parse_chunk_header(seqid: str) -> Tuple[str, int, int]:
    m = re.match(r"^(.+?)__chunk\d+__(\d+)-(\d+)$", seqid)
    if not m:
        raise ValueError(f"Unrecognized chunk header: {seqid}")
    chrom = m.group(1)
    start1 = int(m.group(2))
    end1 = int(m.group(3))
    return chrom, start1 - 1, end1


# -----------------------------
# Step 5: run ltrharvest per chunk (SCN from stdout, GFF3 to file)
# -----------------------------

def build_suffixerator_index(gt_path: str, chunk_fa: str, indexname: str, verbose: bool = False):
    cmd = [
        gt_path, "suffixerator",
        "-db", chunk_fa,
        "-indexname", indexname,
        "-tis", "-suf", "-lcp", "-des", "-ssp", "-sds",
        "-dna",
    ]
    r = run(cmd, verbose=verbose)
    if r.returncode != 0:
        raise RuntimeError(f"gt suffixerator failed for {chunk_fa}:\n{(r.stderr or '').strip()}")

def run_ltrharvest_scn_and_gff3(
    gt_path: str,
    indexname: str,
    scn_path: str,
    gff3_path: str,
    ltrharvest_args: List[str],
    timeout_s: int = 0,
    verbose: bool = False,
):
    """
    - SCN/tabular output is stdout (with -tabout yes)
    - GFF3 output goes to gff3_path (with -gff3 <file>)
    - If timeout occurs, salvage partial stdout to scn_path and keep any partial gff3 file.
    """
    cmd = [gt_path, "ltrharvest", "-index", indexname] + ltrharvest_args + ["-tabout", "yes", "-gff3", gff3_path]

    if verbose:
        print(f"  $ {' '.join(cmd)}", file=sys.stderr)

    try:
        r = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=(timeout_s if timeout_s and timeout_s > 0 else None),
        )
        if r.returncode != 0:
            raise RuntimeError(f"gt ltrharvest failed for index {indexname}:\n{(r.stderr or '').strip()}")
        Path(scn_path).write_text(r.stdout or "")

    except subprocess.TimeoutExpired as te:
        partial = te.stdout or ""
        Path(scn_path).write_text(partial)
        gp = Path(gff3_path)
        if not gp.exists():
            gp.touch()
        print(f"[WARN] ltrharvest timeout for index {indexname} after {timeout_s}s; salvaged partial SCN/GFF3.", file=sys.stderr)

def process_one_chunk(
    gt_path: str,
    ltrfinder_path: str,
    chunk: ChunkInfo,
    work_dir: str,
    ltrharvest_args: List[str],
    ltrfinder_args: List[str],
    ltr_tools: str,
    ltr_timeout_s: int = 0,
    trfmod_path: Optional[str] = None,
    trf_args: Optional[List[str]] = None,
    trf_timeout_s: int = 0,
    trf_min_copy: float = 0.0,
    verbose: bool = False,
) -> Tuple[ChunkInfo, Optional[str], Optional[str], Optional[str]]:
    """
    Returns:
      (chunk,
       ltrharvest_scn_path or None,
       ltrharvest_gff3_path or None,
       ltrfinder_scn_path (converted) or None)
    """
    wdir = mkdirp(Path(work_dir) / chunk.chunk_id)
    local_fa = str(wdir / Path(chunk.fasta_path).name)
    shutil.copyfile(chunk.fasta_path, local_fa)

    # ---- TRF-mod low-complexity masking (optional) ----
    if trfmod_path:
        trf_bed = str(wdir / f"{chunk.chunk_id}.trf.bed")
        trf_masked_fa = str(wdir / f"{chunk.chunk_id}.trf_masked.fa")

        run_trfmod(
            trfmod_path=trfmod_path,
            chunk_fa=local_fa,
            out_bed=trf_bed,
            trf_args=(trf_args or ["-p", "7", "-s", "35"]),
            timeout_s=trf_timeout_s,
            min_copy=trf_min_copy,
            verbose=verbose,
        )

        hardmask_fasta_by_bed0(local_fa, trf_bed, trf_masked_fa)
        local_fa = trf_masked_fa

    harvest_scn_path = None
    harvest_gff3_path = None
    finder_scn_path = None

    if ltr_tools in ("both", "ltrharvest"):
        indexname = str(wdir / "idx")
        build_suffixerator_index(gt_path, local_fa, indexname, verbose=verbose)

        harvest_scn_path = str(wdir / f"{chunk.chunk_id}.ltrharvest.scn")
        harvest_gff3_path = str(wdir / f"{chunk.chunk_id}.ltrharvest.gff3")
        run_ltrharvest_scn_and_gff3(
            gt_path, indexname, harvest_scn_path, harvest_gff3_path, ltrharvest_args,
            timeout_s=ltr_timeout_s, verbose=verbose,
        )

    if ltr_tools in ("both", "ltr_finder"):
        raw_path = str(wdir / f"{chunk.chunk_id}.ltr_finder.raw.scn")
        conv_path = str(wdir / f"{chunk.chunk_id}.ltr_finder.ltrharvest_like.scn")

        run_ltrfinder_raw(ltrfinder_path, local_fa, raw_path, ltrfinder_args, timeout_s=ltr_timeout_s, verbose=verbose)
        raw_text = Path(raw_path).read_text() if Path(raw_path).exists() else ""
        ltrfinder_w2_to_ltrharvest_scn(raw_text, conv_path)

        finder_scn_path = conv_path

    return chunk, harvest_scn_path, harvest_gff3_path, finder_scn_path

# -----------------------------
# Step 5b: ltr_finder per chunk + convert to ltrharvest-style SCN
# -----------------------------

def run_ltrfinder_raw(ltrfinder_path: str, chunk_fa: str, raw_out_path: str, ltrfinder_args: List[str], timeout_s: int = 0, verbose: bool = False):
    """
    Runs ltr_finder on chunk_fa. Captures stdout to raw_out_path.
    If timeout occurs, salvages partial stdout to raw_out_path.
    NOTE: -w 2 is required (user stated); enforced elsewhere.
    """
    cmd = [ltrfinder_path] + ltrfinder_args + [chunk_fa]

    if verbose:
        print(f"  $ {' '.join(cmd)}", file=sys.stderr)

    try:
        r = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=(timeout_s if timeout_s and timeout_s > 0 else None),
        )
        if r.returncode != 0:
            raise RuntimeError(f"ltr_finder failed:\n{(r.stderr or '').strip()}")
        Path(raw_out_path).write_text(r.stdout or "")

    except subprocess.TimeoutExpired as te:
        partial = te.stdout or ""
        Path(raw_out_path).write_text(partial)
        print(f"[WARN] ltr_finder timeout after {timeout_s}s on {chunk_fa}; salvaged partial stdout.", file=sys.stderr)

def ltrfinder_w2_to_ltrharvest_scn(raw_text: str, out_scn_path: str):
    """
    Convert LTR_FINDER -w 2 output to LTRharvest tabout-like SCN.
    """
    seq_id = -1
    n_written = 0

    with open(out_scn_path, "w") as out:
        out.write("## LTR_FINDER\n")
        out.write("## Converted to LTRharvest-like SCN/tabout\n")
        out.write("## s(ret) e(ret) l(ret) s(lLTR) e(lLTR) l(lLTR) s(rLTR) e(rLTR) l(rLTR) sim(LTRs) seq-nr chr\n")

        for line in raw_text.splitlines():
            if not line:
                continue
            if line.startswith(">"):
                seq_id += 1
                continue
            if not line.startswith("["):
                continue

            line2 = re.sub(r"^\[\s+", "[", line.strip())
            toks = line2.split()
            if len(toks) <= 15:
                continue

            chr_ = toks[1]
            loc = toks[2]
            lens = toks[3]
            ltr_len = toks[4]
            direction = toks[12]
            similarity = toks[15]

            mloc = re.match(r"^(\d+)\-(\d+)$", loc)
            mlens = re.match(r"^(\d+),(\d+)$", lens)
            if not mloc or not mlens:
                continue

            sret = int(mloc.group(1))
            eret = int(mloc.group(2))
            lltr_len = int(mlens.group(1))
            rltr_len = int(mlens.group(2))

            eltr = sret + lltr_len - 1
            srtr = eret - rltr_len + 1

            if sret <= 0 or eret <= 0 or eltr < sret or srtr > eret:
                continue

            sim_val = None
            try:
                if similarity.endswith("%"):
                    sim_val = float(similarity.rstrip("%"))
                else:
                    sim_val = float(similarity)
                    if 0.0 <= sim_val <= 1.0:
                        sim_val *= 100.0
            except Exception:
                continue

            try:
                lret = int(float(ltr_len))
            except Exception:
                lret = eret - sret + 1

            if seq_id < 0:
                seq_id = 0

            out.write(
                f"{sret} {eret} {lret} "
                f"{sret} {eltr} {lltr_len} "
                f"{srtr} {eret} {rltr_len} "
                f"{sim_val:.2f} {seq_id} {chr_}\n"
            )
            n_written += 1

    if n_written == 0:
        Path(out_scn_path).touch()

# -----------------------------
# Step 6: stitch outputs (lift coords and normalize seqid)
# -----------------------------

def stitch_gff3(chunk_triplets: List[Tuple[ChunkInfo, str, str]], stitched_out: str):
    with open(stitched_out, "w") as out:
        out.write("##gff-version 3\n")
        for chunk, _scn, gff_path in chunk_triplets:
            gp = Path(gff_path)
            if not gp.exists() or gp.stat().st_size == 0:
                continue
            with gp.open("r") as f:
                for line in f:
                    if not line.strip():
                        continue
                    if line.startswith("##gff-version"):
                        continue
                    if line.startswith("##sequence-region"):
                        parts = line.rstrip("\n").split()
                        if len(parts) >= 4:
                            seqid = parts[1]
                            try:
                                chrom, start0, end1 = parse_chunk_header(seqid)
                                out.write(f"##sequence-region\t{chrom}\t{start0+1}\t{end1}\n")
                            except Exception:
                                continue
                        continue
                    if line.startswith("#"):
                        continue

                    cols = line.rstrip("\n").split("\t")
                    if len(cols) < 9:
                        continue
                    seqid = cols[0]
                    try:
                        start = int(cols[3])
                        end = int(cols[4])
                    except ValueError:
                        continue

                    try:
                        chrom, chunk_start0, _chunk_end0 = parse_chunk_header(seqid)
                    except Exception:
                        chrom = chunk.chrom
                        chunk_start0 = chunk.start0

                    cols[0] = chrom
                    cols[3] = str(start + chunk_start0)
                    cols[4] = str(end + chunk_start0)
                    out.write("\t".join(cols) + "\n")

def _is_valid_require_run_chars(require_run_chars: Optional[List[str]]) -> List[str]:
    req: List[str] = []
    if not require_run_chars:
        return req
    for c in require_run_chars:
        c = (c or "").strip().upper()
        if not c:
            continue
        if len(c) != 1:
            raise ValueError(f"--require-run-chars expects single characters, got: {c!r}")
        req.append(c)
    return req


def has_nested_run_signature(
    seq: str,
    chars: List[str],
    base_min: int = 800,
    flank_min: int = 80,
) -> bool:
    """
    True iff seq contains an adjacent nested pattern:
      chars[-1]{>=flank} ... chars[1]{>=flank} chars[0]{>=base} chars[1]{>=flank} ... chars[-1]{>=flank}
    """
    if not chars:
        return True

    if base_min <= 0:
        base_min = 1
    if flank_min <= 0:
        flank_min = 1

    s = seq.upper()
    c0 = chars[0]

    n = len(s)
    i = 0

    while i < n:
        if s[i] != c0:
            i += 1
            continue

        j = i
        while j < n and s[j] == c0:
            j += 1

        run_len = j - i
        if run_len >= base_min:
            left_idx = i
            right_idx = j

            ok = True
            for level in range(1, len(chars)):
                c = chars[level]

                l_end = left_idx
                l_start = l_end
                while l_start > 0 and s[l_start - 1] == c:
                    l_start -= 1
                left_run = l_end - l_start
                if left_run < flank_min:
                    ok = False
                    break

                r_start = right_idx
                r_end = r_start
                while r_end < n and s[r_end] == c:
                    r_end += 1
                right_run = r_end - r_start
                if right_run < flank_min:
                    ok = False
                    break

                left_idx = l_start
                right_idx = r_end

            if ok:
                return True

        i = j

    return False


def has_exclude_run_char(seq: str, char: str, min_run: int = 10) -> bool:
    """
    Returns True if seq contains a run of `char` of length >= min_run anywhere.
    Case-insensitive. Used to EXCLUDE candidates containing such a run.
    """
    char = char.upper()
    s = seq.upper()
    count = 0
    for c in s:
        if c == char:
            count += 1
            if count >= min_run:
                return True
        else:
            count = 0
    return False


def stitch_scn(chunk_triplets: List[Tuple[ChunkInfo, str, str]], stitched_out: str):
    """
    SCN/tabular line format (ltrharvest tabout):
      s(ret) e(ret) l(ret) s(lLTR) e(lLTR) l(lLTR) s(rLTR) e(rLTR) l(rLTR) sim(LTRs) seq-nr chrom
    """
    n_written = 0
    with open(stitched_out, "w") as out:
        for chunk, scn_path, _gff in chunk_triplets:
            sp = Path(scn_path)
            if not sp.exists() or sp.stat().st_size == 0:
                continue
            with sp.open("r") as f:
                for line in f:
                    if not line.strip():
                        continue
                    if line.startswith("#"):
                        continue

                    parts = re.split(r"\s+", line.strip())
                    if len(parts) < 11:
                        continue
                    if not (parts[0].isdigit() and parts[1].isdigit()):
                        continue

                    off = chunk.start0
                    for idx in (0, 1, 3, 4, 6, 7):
                        if idx < len(parts) and re.fullmatch(r"\d+", parts[idx]):
                            parts[idx] = str(int(parts[idx]) + off)
                        else:
                            parts = None
                            break
                    if parts is None:
                        continue

                    out.write("  ".join(parts) + f" {chunk.chrom}\n")
                    n_written += 1

    if n_written == 0:
        Path(stitched_out).touch()

def stitch_scn_from_pairs(pairs: List[Tuple[ChunkInfo, str]], stitched_out: str):
    triplets = [(c, scn, "") for (c, scn) in pairs]
    stitch_scn(triplets, stitched_out)

def merge_stitched_scns(stitched_scns: List[str], merged_out: str):
    """
    Merge multiple stitched SCN files into one, removing duplicates.
    """
    seen = set()
    n_written = 0

    with open(merged_out, "w") as out:
        for scn in stitched_scns:
            p = Path(scn)
            if not p.exists() or p.stat().st_size == 0:
                continue
            with p.open("r") as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    parts = re.split(r"\s+", line)
                    if len(parts) < 12:
                        continue
                    chrom = parts[-1]
                    key = (
                        chrom,
                        parts[0], parts[1],
                        parts[3], parts[4],
                        parts[6], parts[7],
                    )
                    if key in seen:
                        continue
                    seen.add(key)
                    out.write("  ".join(parts) + "\n")
                    n_written += 1

    if n_written == 0:
        Path(merged_out).touch()
        
def filter_scn_by_lengths(
    in_scn: str,
    out_scn: str,
    min_ltr_len: int = 0,
    min_ret_len: int = 0,
    max_ret_len: int = 0,
    min_int_len: int = 0,
    max_int_len: int = 0,
) -> Dict[str, int]:
    """
    Filters an LTRharvest tabout-like SCN file by size criteria.
    """
    counts = {
        "kept": 0,
        "total": 0,
        "fail_min_ltr": 0,
        "fail_ret_range": 0,
        "fail_int_range": 0,
        "malformed": 0,
    }

    inp = Path(in_scn)
    if not inp.exists() or inp.stat().st_size == 0:
        Path(out_scn).touch()
        return counts

    with open(in_scn, "r") as fin, open(out_scn, "w") as out:
        for line in fin:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = re.split(r"\s+", line)
            if len(parts) < 12:
                counts["malformed"] += 1
                continue

            counts["total"] += 1

            try:
                lret = int(parts[2])
                lltr = int(parts[5])
                rltr = int(parts[8])
            except ValueError:
                counts["malformed"] += 1
                continue

            if min_ltr_len and (lltr < min_ltr_len or rltr < min_ltr_len):
                counts["fail_min_ltr"] += 1
                continue

            if min_ret_len and lret < min_ret_len:
                counts["fail_ret_range"] += 1
                continue
            if max_ret_len and lret > max_ret_len:
                counts["fail_ret_range"] += 1
                continue

            internal = lret - lltr - rltr
            if min_int_len and internal < min_int_len:
                counts["fail_int_range"] += 1
                continue
            if max_int_len and internal > max_int_len:
                counts["fail_int_range"] += 1
                continue

            out.write("  ".join(parts) + "\n")
            counts["kept"] += 1

    if counts["kept"] == 0:
        Path(out_scn).touch()

    return counts

# -----------------------------
# Step 7: build LTR FASTA from stitched SCN
# -----------------------------

def scn_to_internal_fasta(
    stitched_scn: str,
    genome_fa: str,
    out_fa: str,
    require_run_chars: Optional[List[str]] = None,
    exclude_run_char: Optional[str] = None,
    base_min: int = 800,
    flank_min: int = 80,
):
    """
    Uses e(lLTR) and s(rLTR) (1-based inclusive coords from ltrharvest tabout),
    plus chrom (last column), to extract INTERNAL sequence only (LTRs excluded).
    """
    genome = load_fasta_as_dict(genome_fa)

    req = []
    if require_run_chars:
        for c in require_run_chars:
            c = (c or "").strip().upper()
            if c:
                req.append(c)

    n_written = 0
    with open(out_fa, "w") as out, open(stitched_scn, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = re.split(r"\s+", line)
            if len(parts) < 12:
                continue

            sret, eret = parts[0], parts[1]
            el = parts[4]
            sr = parts[6]
            chrom = parts[-1]

            if not (sret.isdigit() and eret.isdigit() and el.isdigit() and sr.isdigit()):
                continue

            sret1 = int(sret)
            eret1 = int(eret)
            el1   = int(el)
            sr1   = int(sr)

            if eret1 < sret1:
                sret1, eret1 = eret1, sret1

            internal_start1 = el1 + 1
            internal_end1   = sr1 - 1
            if internal_end1 < internal_start1:
                continue

            seq = genome.get(chrom)
            if seq is None:
                continue

            fl_start0 = sret1 - 1
            fl_end0   = eret1
            if fl_start0 < 0 or fl_end0 > len(seq):
                continue

            req = _is_valid_require_run_chars(require_run_chars)

            if req:
                full_len_seq = seq[fl_start0:fl_end0].upper()
                if not has_nested_run_signature(full_len_seq, req, base_min=base_min, flank_min=flank_min):
                    continue

            # Exclude run char filter on FULL LENGTH
            if exclude_run_char:
                full_len_seq = seq[fl_start0:fl_end0].upper()
                if has_exclude_run_char(full_len_seq, exclude_run_char):
                    continue

            start0 = internal_start1 - 1
            end0   = internal_end1
            if start0 < 0 or end0 > len(seq):
                continue

            header = f"{chrom}:{sret1}-{eret1}"
            out.write(f">{header}\n")
            frag = seq[start0:end0]
            for i in range(0, len(frag), 60):
                out.write(frag[i:i+60] + "\n")
            n_written += 1

    if n_written == 0:
        Path(out_fa).touch()

def hardmask_fasta_by_bed0(in_fasta: str, bed_path: str, out_fasta: str) -> int:
    """
    Hardmask regions in a FASTA using a BED file with:
      chrom  start0  end0  ...
    Assumes coordinates are 0-based, end-exclusive (BED convention).
    Returns number of intervals applied.
    """
    intervals: Dict[str, List[Interval]] = {}
    n = 0

    bp = Path(bed_path)
    if not bp.exists() or bp.stat().st_size == 0:
        shutil.copyfile(in_fasta, out_fasta)
        return 0

    with bp.open("r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 3:
                continue
            chrom = cols[0]
            try:
                s = int(cols[1])
                e = int(cols[2])
            except ValueError:
                continue
            if e <= s:
                continue
            intervals.setdefault(chrom, []).append(Interval(s, e))
            n += 1

    intervals = merge_intervals_per_chrom(intervals)
    if not intervals:
        shutil.copyfile(in_fasta, out_fasta)
        return 0

    hardmask_fasta_by_intervals(in_fasta, intervals, out_fasta, out_bed=str(bp))
    return n

def run_trfmod(
    trfmod_path: str,
    chunk_fa: str,
    out_bed: str,
    trf_args: List[str],
    timeout_s: int = 0,
    min_copy: float = 0.0,
    verbose: bool = False,
):
    """
    Runs trf-mod on chunk_fa, filters by copy number, and writes BED0 to out_bed.
    """
    cmd = [trfmod_path] + trf_args + [chunk_fa]

    if verbose:
        print(f"  $ {' '.join(cmd)}", file=sys.stderr)

    def trf_text_to_bed0(text: str) -> str:
        out_lines = []
        for line in text.splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 5:
                continue

            chrom = cols[0]
            try:
                start1 = int(cols[1])
                end1 = int(cols[2])
                copy_num = float(cols[4])
            except ValueError:
                continue

            if min_copy and copy_num < min_copy:
                continue

            if end1 < start1:
                start1, end1 = end1, start1

            start0 = start1 - 1
            end0 = end1

            if end0 <= start0:
                continue

            out_lines.append(f"{chrom}\t{start0}\t{end0}\tcopy={copy_num:.2f}")
        return "\n".join(out_lines) + ("\n" if out_lines else "")

    try:
        r = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=(timeout_s if timeout_s and timeout_s > 0 else None),
        )
        if r.returncode != 0:
            raise RuntimeError(f"trf-mod failed:\n{(r.stderr or '').strip()}")
        Path(out_bed).write_text(trf_text_to_bed0(r.stdout or ""))

    except subprocess.TimeoutExpired as te:
        partial = te.stdout or ""
        Path(out_bed).write_text(trf_text_to_bed0(partial))
        print(f"[WARN] trf-mod timeout after {timeout_s}s on {chunk_fa}; salvaged + filtered BED.", file=sys.stderr)

def scn_to_intact_fasta(
    stitched_scn: str,
    genome_fa: str,
    out_fa: str,
    require_run_chars: Optional[List[str]] = None,
    exclude_run_char: Optional[str] = None,
    base_min: int = 800,
    flank_min: int = 80,
):
    """
    Extract FULL intact LTR-RT sequence using s(ret) and e(ret) (1-based inclusive),
    plus chrom (last column).
    """
    genome = load_fasta_as_dict(genome_fa)

    req = []
    if require_run_chars:
        for c in require_run_chars:
            c = (c or "").strip().upper()
            if c:
                req.append(c)

    n_written = 0
    with open(out_fa, "w") as out, open(stitched_scn, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = re.split(r"\s+", line)
            if len(parts) < 12:
                continue

            sret, eret = parts[0], parts[1]
            chrom = parts[-1]

            if not (sret.isdigit() and eret.isdigit()):
                continue

            sret1 = int(sret)
            eret1 = int(eret)
            if eret1 < sret1:
                sret1, eret1 = eret1, sret1

            seq = genome.get(chrom)
            if seq is None:
                continue

            start0 = sret1 - 1
            end0 = eret1
            if start0 < 0 or end0 > len(seq):
                continue

            frag = seq[start0:end0]

            req = _is_valid_require_run_chars(require_run_chars)

            if req:
                full_len_seq = frag.upper()
                if not has_nested_run_signature(full_len_seq, req, base_min=base_min, flank_min=flank_min):
                    continue

            # Exclude run char filter on FULL LENGTH
            if exclude_run_char:
                if has_exclude_run_char(frag.upper(), exclude_run_char):
                    continue

            header = f"{chrom}:{sret1}-{eret1}"
            out.write(f">{header}\n")
            for i in range(0, len(frag), 60):
                out.write(frag[i:i+60] + "\n")
            n_written += 1

    if n_written == 0:
        Path(out_fa).touch()


def load_tesorter_te_to_annotation(cls_tsv_path: str) -> Dict[str, str]:
    """
    Build mapping:
      TE (e.g. CP002684.1:7065623-7075763) -> "Order/Superfamily/Clade"
    """
    te2ann: Dict[str, str] = {}

    p = Path(cls_tsv_path)
    if not p.exists() or p.stat().st_size == 0:
        return te2ann

    with p.open("r") as f:
        header_seen = False
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue

            if not header_seen:
                header_seen = True
                if line.startswith("TE\t") or line.lower().startswith("te\t"):
                    continue

            cols = line.split("\t")
            if len(cols) < 4:
                continue

            te, order, superfam, clade = cols[0], cols[1], cols[2], cols[3]
            if not te:
                continue

            order = order or "unknown"
            superfam = superfam or "unknown"
            clade = clade or "unknown"

            te2ann[te] = f"{order}/{superfam}/{clade}"

    return te2ann


def load_tesorter_gff3_domains(gff3_path: str) -> Dict[str, list]:
    """Parse a TEsorter .dom.gff3 file into per-TE domain lists.

    Returns:
        Dict mapping 'chrom:start-end' -> list of dicts, each with keys:
        gene, clade, start (int, genomic), end (int, genomic), evalue, coverage.
    """
    result: Dict[str, list] = {}
    p = Path(gff3_path)
    if not p.exists() or p.stat().st_size == 0:
        return result

    with p.open("r") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 9:
                continue

            gff_start = int(cols[3])
            gff_end = int(cols[4])
            attrs_str = cols[8]

            # Parse attributes: ID=chr1:123-456|...; gene=RT; clade=...; ...
            attrs: Dict[str, str] = {}
            for field in attrs_str.split(";"):
                field = field.strip()
                if "=" in field:
                    k, v = field.split("=", 1)
                    attrs[k] = v

            # Extract TE name from ID attribute (before the '|' delimiter)
            id_val = attrs.get("ID", "")
            te_name = id_val.split("|", 1)[0] if "|" in id_val else id_val

            if not te_name:
                continue

            domain_entry = {
                "gene": attrs.get("gene", ""),
                "clade": attrs.get("clade", ""),
                "start": gff_start,
                "end": gff_end,
                "evalue": attrs.get("evalue", ""),
                "coverage": attrs.get("coverage", ""),
            }
            result.setdefault(te_name, []).append(domain_entry)

    return result


# -----------------------------
# Step 8: build kmer2ltr.domain from stitched SCN
# -----------------------------

def scn_to_kmer2ltr_domain(stitched_scn: str, out_domain: str, tesorter_cls_tsv: Optional[str] = None):
    """
    Writes:
      {chrom}:{s(ret)}-{e(ret)}#{Order/Superfamily/Clade} <TAB> max(l(lLTR), l(rLTR))
    """
    te2ann: Dict[str, str] = {}
    if tesorter_cls_tsv:
        te2ann = load_tesorter_te_to_annotation(tesorter_cls_tsv)

    n_written = 0
    with open(stitched_scn, "r") as fin, open(out_domain, "w") as out:
        for line in fin:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            c = line.split()
            if len(c) < 12:
                continue

            sret, eret = c[0], c[1]
            chrom = c[-1]

            if not (sret.isdigit() and eret.isdigit()):
                continue
            try:
                ll = int(c[5])
                rl = int(c[8])
            except ValueError:
                continue

            te_key = f"{chrom}:{sret}-{eret}"
            ann = te2ann.get(te_key, "LTR/unknown/unknown")

            out.write(f"{te_key}#{ann}\t{max(ll, rl)}\n")
            n_written += 1

    if n_written == 0:
        Path(out_domain).touch()

# -----------------------------
# Tool setup (./tools/)
# -----------------------------

def tool_usable_kmer2ltr(k2l_dir: Path) -> bool:
    script = (k2l_dir / "Kmer2LTR.py").resolve()
    if not script.exists():
        return False
    try:
        r = subprocess.run(
            ["python3", str(script), "-h"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        txt = (r.stdout or "") + "\n" + (r.stderr or "")
        return ("--max-win-overdisp" in txt) or ("Kmer2LTR" in txt) or (r.returncode == 0)
    except Exception:
        return False

def ensure_kmer2ltr(tools_dir: Path) -> str:
    tools_dir = mkdirp(tools_dir)
    k2l_dir = tools_dir / "Kmer2LTR"

    if not tool_usable_kmer2ltr(k2l_dir):
        if not k2l_dir.exists():
            run(["git", "clone", "https://github.com/cwb14/Kmer2LTR.git", str(k2l_dir)], check=True)
        if not tool_usable_kmer2ltr(k2l_dir):
            raise RuntimeError(
                f"Kmer2LTR appears unusable in: {k2l_dir}\n"
                f"Try: python3 {k2l_dir}/Kmer2LTR.py -h"
            )

    return str((k2l_dir / "Kmer2LTR.py").resolve())


def tool_usable_generic(bin_path: Path, args: List[str]) -> bool:
    if not bin_path.exists() or not os.access(str(bin_path), os.X_OK):
        return False
    try:
        r = run([str(bin_path)] + args, capture=True)
        return r.returncode == 0
    except Exception:
        return False

def tool_usable_miniprot(bin_path: Path) -> bool:
    if not bin_path.exists() or not os.access(str(bin_path), os.X_OK):
        return False
    try:
        r = run([str(bin_path), "-h"], capture=True)
        txt = (r.stdout or "") + "\n" + (r.stderr or "")
        return ("Usage: miniprot" in txt) or (r.returncode == 0)
    except Exception:
        return False

def tool_usable_minimap2(bin_path: Path) -> bool:
    return tool_usable_generic(bin_path, ["--help"])

def ensure_tools(tools_dir: Path) -> Tuple[str, str]:
    tools_dir = mkdirp(tools_dir)
    mm2_dir = tools_dir / "minimap2"
    mp_dir  = tools_dir / "miniprot"

    mm2_bin = mm2_dir / "minimap2"
    mp_bin  = mp_dir  / "miniprot"

    # ---- minimap2 ----
    if not tool_usable_minimap2(mm2_bin):
        if not mm2_dir.exists():
            run(["git", "clone", "https://github.com/lh3/minimap2", str(mm2_dir)], check=True)

        built = False
        for make_cmd in (["make"], ["make", "CC=gcc"]):
            r = run(make_cmd, cwd=str(mm2_dir))
            if r.returncode == 0 and tool_usable_minimap2(mm2_bin):
                built = True
                break

        if not built:
            sys_mm2 = shutil.which("minimap2")
            if sys_mm2 and tool_usable_minimap2(Path(sys_mm2)):
                mm2_bin = Path(sys_mm2)
            else:
                raise RuntimeError(
                    f"minimap2 build failed and no system minimap2 found on PATH.\n"
                    f"Try: conda install -c bioconda minimap2"
                )

    # ---- miniprot ----
    if not tool_usable_miniprot(mp_bin):
        if not mp_dir.exists():
            run(["git", "clone", "https://github.com/lh3/miniprot", str(mp_dir)], check=True)

        built = False
        for make_cmd in (["make"], ["make", "CC=gcc"]):
            r = run(make_cmd, cwd=str(mp_dir))
            if r.returncode == 0 and tool_usable_miniprot(mp_bin):
                built = True
                break

        if not built:
            sys_mp = shutil.which("miniprot")
            if sys_mp and tool_usable_miniprot(Path(sys_mp)):
                mp_bin = Path(sys_mp)
            else:
                raise RuntimeError(
                    f"miniprot build failed and no system miniprot found on PATH.\n"
                    f"Try: conda install -c bioconda miniprot"
                )

    return str(mm2_bin), str(mp_bin)
    
# -----------------------------
# Step 9: TEsorter
# -----------------------------

def tool_usable_tesorter(te_dir: Path) -> bool:
    if not te_dir.exists():
        return False
    try:
        env = dict(os.environ)
        env["PYTHONPATH"] = str(te_dir)
        r = subprocess.run(
            ["python3", "-m", "TEsorter", "-h"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env,
        )
        txt = (r.stdout or "") + "\n" + (r.stderr or "")
        return ("usage: TEsorter" in txt) or (r.returncode == 0)
    except Exception:
        return False

def ensure_tesorter(tools_dir: Path) -> str:
    tools_dir = mkdirp(tools_dir)
    te_dir = tools_dir / "TEsorter"

    if not tool_usable_tesorter(te_dir):
        if not te_dir.exists():
            run([
                "git", "clone",
                "--branch", "my-new-idea2",
                "--single-branch",
                "https://github.com/cwb14/TEsorter.git",
                str(te_dir)
            ], check=True)

        if not tool_usable_tesorter(te_dir):
            raise RuntimeError(
                f"TEsorter appears unusable in: {te_dir}\n"
                f"Try: PYTHONPATH={te_dir} python3 -m TEsorter -h"
            )

    return str(te_dir.resolve())

def tool_usable_trfmod(bin_path: Path) -> bool:
    if not bin_path.exists() or not os.access(str(bin_path), os.X_OK):
        return False
    try:
        r = subprocess.run(
            [str(bin_path), "-h"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        txt = (r.stdout or "") + "\n" + (r.stderr or "")
        return ("trf" in txt.lower()) or (r.returncode == 0)
    except Exception:
        return False


def ensure_trfmod(tools_dir: Path) -> str:
    tools_dir = mkdirp(tools_dir)
    trf_dir = tools_dir / "TRF-mod"
    trf_bin = trf_dir / "trf-mod"

    if not tool_usable_trfmod(trf_bin):
        if not trf_dir.exists():
            run(["git", "clone", "https://github.com/lh3/TRF-mod", str(trf_dir)], check=True)
        run(["make", "-f", "compile.mak"], cwd=str(trf_dir), check=True)

        if not tool_usable_trfmod(trf_bin):
            raise RuntimeError(f"TRF-mod build failed or binary unusable at: {trf_bin}")

    return str(trf_bin)

def run_tesorter(stitched_fa: str, tesorter_py_path: str, outdir: Path,
                 db: str, cov: int, evalue: str, rule: str, threads: int,
                 pass2_classified_fasta: Optional[str] = None,
                 no_cleanup: bool = False,
                 verbose: bool = False) -> Tuple[str, str, str]:
    """
    Runs TEsorter in outdir so outputs land in {prefix}.work/.
    """
    mkdirp(outdir)

    env = dict(os.environ)

    te_abs = str(Path(tesorter_py_path).resolve())
    env_py = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = te_abs if not env_py else (te_abs + os.pathsep + env_py)

    stitched_fa_abs = str(Path(stitched_fa).resolve())

    cmd = [
        "python3", "-m", "TEsorter",
        stitched_fa_abs,
        "-db", db,
        "-p", str(threads),
        "-cov", str(cov),
        "-eval", str(evalue),
        "-rule", str(rule),
    ]

    if no_cleanup:
        cmd.append("--no-cleanup")

    if pass2_classified_fasta:
        p2 = Path(pass2_classified_fasta).resolve()
        if not p2.exists() or p2.stat().st_size == 0:
            raise RuntimeError(f"--pass2-classified-fasta provided but missing/empty: {p2}")
        cmd += ["--pass2-classified-fasta", str(p2)]

    if verbose:
        label = f"  $ {' '.join(cmd)}\n    (cwd: {outdir})"
        print(label, file=sys.stderr)

    r = subprocess.run(cmd, cwd=str(outdir), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env)
    if r.returncode != 0:
        raise RuntimeError(f"TEsorter failed:\n{(r.stderr or '').strip()}")

    base = Path(stitched_fa_abs).name
    cls_lib = str(outdir / f"{base}.{db}.cls.lib")
    cls_pep = str(outdir / f"{base}.{db}.cls.pep")
    cls_tsv = str(outdir / f"{base}.{db}.cls.tsv")
    return cls_lib, cls_pep, cls_tsv
    
def run_kmer2ltr(kmer2ltr_py: str, in_fa: str, out_prefix: str, outdir: Path,
                threads: int, max_win_overdisp: float, min_retained_fraction: float,
                domain_file: Optional[str] = None,
                wfa_align: bool = False,
                verbose: bool = False) -> str:
    """
    Runs Kmer2LTR in outdir.
    Returns the path to the main TSV output (the file named exactly out_prefix).
    """
    mkdirp(outdir)
    which_or_die("mafft")
    which_or_die("trimal")

    kmer2ltr_py_abs = str(Path(kmer2ltr_py).resolve())

    cmd = [
        "python3", kmer2ltr_py_abs,
        "-i", str(Path(in_fa).resolve()),
        "--max-win-overdisp", str(max_win_overdisp),
        "--min-retained-fraction", str(min_retained_fraction),
        "-p", str(threads),
        "--purge-subdirs",
        "-o", out_prefix,
    ]

    if domain_file:
        df = Path(domain_file).resolve()
        if not df.exists() or df.stat().st_size == 0:
            raise RuntimeError(f"Kmer2LTR domain file requested but missing/empty: {df}")
        cmd += ["-D", str(df)]

    if wfa_align:
        cmd += ["--wfa-align"]

    if verbose:
        label = f"  $ {' '.join(cmd)}\n    (cwd: {outdir})"
        print(label, file=sys.stderr)

    # Let stderr pass through so Kmer2LTR's progress counter is visible
    r = subprocess.run(cmd, cwd=str(outdir), stdout=subprocess.PIPE, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"Kmer2LTR failed (exit {r.returncode})")

    main_out = outdir / out_prefix
    if not main_out.exists():
        raise RuntimeError(f"Kmer2LTR finished but did not create expected output: {main_out}")

    return str(main_out)
    
def _parse_interval_from_kmer2ltr_col1(col1: str) -> Optional[Tuple[str, int, int]]:
    left = col1.split("#", 1)[0]
    if ":" not in left:
        return None
    chrom, rest = left.split(":", 1)
    if "-" not in rest:
        return None
    start_s, end_s = rest.split("-", 1)
    try:
        start = int(start_s)
        end = int(end_s)
    except ValueError:
        return None
    if end < start:
        start, end = end, start
    return chrom, start, end


def _tsd_names_from_fasta(fa_path: str) -> Set[str]:
    """Return set of 'chrom:start-end' keys parsed from TSD-positive FASTA headers.

    Expected header format: >chrom:start-end#classification ...
    """
    names: Set[str] = set()
    p = Path(fa_path)
    if not p.exists() or p.stat().st_size == 0:
        return names
    with open(fa_path) as f:
        for line in f:
            if line.startswith(">"):
                key = line[1:].strip().split("#")[0].split()[0]
                names.add(key)
    return names


def pre_purge_tsd_dominated(
    stitched_scn: str,
    tsd_names: Set[str],
    threshold: float,
    ltr_bounds: Optional[Dict[str, Tuple[int, int, int, int]]] = None,
) -> Set[str]:
    """
    Identify non-TSD SCN candidates whose genomic interval overlaps a TSD+
    candidate at >= *threshold* (fraction of the shorter interval).

    These candidates would be eliminated by Layer-1 of dedup anyway, so
    removing them before TEsorter / Kmer2LTR saves runtime without
    changing the final library.

    Candidates that are fully contained within a TSD+ candidate AND have
    distinct (non-overlapping) LTRs are protected as putative nested TEs.

    Returns the set of 'chrom:start-end' keys to EXCLUDE.
    """
    by_chrom_tsd: Dict[str, List[Tuple[int, int]]] = {}
    by_chrom_nontsd: Dict[str, List[Tuple[int, int, str]]] = {}
    seen: Set[str] = set()

    with open(stitched_scn) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = re.split(r"\s+", line)
            if len(parts) < 12:
                continue
            sret_s, eret_s = parts[0], parts[1]
            chrom = parts[-1]
            if not (sret_s.isdigit() and eret_s.isdigit()):
                continue
            s, e = int(sret_s), int(eret_s)
            if e < s:
                s, e = e, s
            key = f"{chrom}:{s}-{e}"
            if key in seen:
                continue
            seen.add(key)

            if key in tsd_names:
                by_chrom_tsd.setdefault(chrom, []).append((s, e))
            else:
                by_chrom_nontsd.setdefault(chrom, []).append((s, e, key))

    exclude: Set[str] = set()
    for chrom, nontsd_intervals in by_chrom_nontsd.items():
        tsd_intervals = by_chrom_tsd.get(chrom)
        if not tsd_intervals:
            continue
        for s, e, key in nontsd_intervals:
            # Check ALL overlapping TSD+ intervals before deciding.
            # If the non-TSD is contained in any TSD+ interval, keep it
            # (putative nested TE).  If it only partially overlaps, exclude.
            is_contained_in_any = False
            should_exclude = False
            for ts, te in tsd_intervals:
                if _overlap_fraction_of_shorter((s, e), (ts, te)) >= threshold:
                    tsd_key = f"{chrom}:{ts}-{te}"
                    if (_is_contained((s, e), (ts, te)) == "a_in_b"
                            and not _ltrs_shared(key, tsd_key, ltr_bounds)):
                        is_contained_in_any = True
                    else:
                        should_exclude = True
            if is_contained_in_any:
                print(f"[pre_purge] keeping contained candidate "
                      f"{key} (putative nested TE, distinct LTRs)")
            elif should_exclude:
                exclude.add(key)

    return exclude


def _overlap_fraction_of_shorter(a: Tuple[int, int], b: Tuple[int, int]) -> float:
    a0, a1 = a
    b0, b1 = b
    inter = min(a1, b1) - max(a0, b0) + 1
    if inter <= 0:
        return 0.0
    lena = a1 - a0 + 1
    lenb = b1 - b0 + 1
    return inter / min(lena, lenb)


def _is_contained(a: Tuple[int, int], b: Tuple[int, int]) -> str:
    """Check if one interval is strictly contained within the other.

    Returns:
        "a_in_b"  if a is strictly inside b (a shorter than b),
        "b_in_a"  if b is strictly inside a (b shorter than a),
        ""        if neither (partial overlap, identical, or disjoint).
    """
    a0, a1 = a
    b0, b1 = b
    len_a = a1 - a0 + 1
    len_b = b1 - b0 + 1
    if len_a == len_b:
        return ""
    if a0 >= b0 and a1 <= b1:
        return "a_in_b"
    if b0 >= a0 and b1 <= a1:
        return "b_in_a"
    return ""


def _ltrs_shared(
    key_a: str,
    key_b: str,
    ltr_bounds: Optional[Dict[str, Tuple[int, int, int, int]]],
) -> bool:
    """Return True if two elements share (overlap) any LTR.

    An intact nested insertion has its own distinct LTR pair.  If the inner
    and outer share a left or right LTR, it is a truncation / extension
    variant, not true nesting.
    """
    if ltr_bounds is None:
        return False
    ab = ltr_bounds.get(key_a)
    bb = ltr_bounds.get(key_b)
    if ab is None or bb is None:
        return False
    a_lL_s, a_lL_e, a_rL_s, a_rL_e = ab
    b_lL_s, b_lL_e, b_rL_s, b_rL_e = bb
    # Check all four pairwise LTR overlaps (any overlap = shared)
    def _overlaps(s1, e1, s2, e2):
        return s1 <= e2 and s2 <= e1
    if _overlaps(a_lL_s, a_lL_e, b_lL_s, b_lL_e):
        return True
    if _overlaps(a_lL_s, a_lL_e, b_rL_s, b_rL_e):
        return True
    if _overlaps(a_rL_s, a_rL_e, b_lL_s, b_lL_e):
        return True
    if _overlaps(a_rL_s, a_rL_e, b_rL_s, b_rL_e):
        return True
    return False


def load_scn_ltr_boundaries(
    scn_path: str,
) -> Dict[str, Tuple[int, int, int, int]]:
    """Parse a merged SCN file and return per-element LTR boundaries.

    Returns dict: "chrom:s(ret)-e(ret)" -> (s_lLTR, e_lLTR, s_rLTR, e_rLTR).
    """
    result: Dict[str, Tuple[int, int, int, int]] = {}
    with open(scn_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = re.split(r"\s+", line)
            if len(parts) < 12:
                continue
            sret_s, eret_s = parts[0], parts[1]
            chrom = parts[-1]
            if not (sret_s.isdigit() and eret_s.isdigit()):
                continue
            try:
                s_lLTR = int(parts[3])
                e_lLTR = int(parts[4])
                s_rLTR = int(parts[6])
                e_rLTR = int(parts[7])
            except (ValueError, IndexError):
                continue
            s, e = int(sret_s), int(eret_s)
            if e < s:
                s, e = e, s
            key = f"{chrom}:{s}-{e}"
            result[key] = (s_lLTR, e_lLTR, s_rLTR, e_rLTR)
    return result


def dedup_kmer2ltr_tsv(kmer2ltr_tsv: str, out_tsv: str, threshold: float,
                       tsd_names: Optional[Set[str]] = None,
                       gff3_domains: Optional[Dict[str, list]] = None,
                       ltr_bounds: Optional[Dict[str, Tuple[int, int, int, int]]] = None,
                       ) -> None:
    """
    Dedup a Kmer2LTR output TSV (main output file), keeping the lowest p-distance (col7).
    Tie-breaker: if p-distance ties, keep the record with the largest aln_len (col3).

    For "extension" pairs — records sharing a start or end coordinate — one is a
    truncation of the other.  The extension yields slightly higher p-distance but is
    the more complete element.  Resolution order:
      1. If tsd_names is provided and exactly one record has a TSD, keep the TSD record.
      2. Otherwise fall back to lowest p-distance, but print a warning for pairs where
         no TSD information is available so we can study the scale of the problem.
    """
    in_path = Path(kmer2ltr_tsv)
    if not in_path.exists() or in_path.stat().st_size == 0:
        Path(out_tsv).touch()
        return

    tmp_sorted = in_path.parent / (in_path.name + ".sorted.tmp")

    r = run(["sort", "-t:", "-k1,1V", "-k2,2n", str(in_path)], check=True, capture=True)
    tmp_sorted.write_text(r.stdout or "")

    cluster: List[Dict[str, object]] = []
    cluster_chrom: Optional[str] = None

    def flush_cluster(out_handle):
        if not cluster:
            return
        if len(cluster) == 1:
            out_handle.write(cluster[0]["line"])
            return

        # --- Layer 1: TSD tiebreaker (all overlap pairs, not just extension pairs) ---
        # Skip TSD domination for containment pairs — both elements are
        # deferred to Layer 1.5 (nesting diagnostics) instead.
        dominated: Set[int] = set()
        for ii in range(len(cluster)):
            for jj in range(ii + 1, len(cluster)):
                ri, rj = cluster[ii], cluster[jj]
                si, ei = int(ri['s']), int(ri['e'])
                sj, ej = int(rj['s']), int(rj['e'])

                key_i = f"{ri['chrom']}:{si}-{ei}"
                key_j = f"{rj['chrom']}:{sj}-{ej}"

                # Skip TSD domination for true nesting (contained + distinct LTRs)
                if (_is_contained((si, ei), (sj, ej))
                        and not _ltrs_shared(key_i, key_j, ltr_bounds)):
                    continue
                has_i = tsd_names is not None and key_i in tsd_names
                has_j = tsd_names is not None and key_j in tsd_names
                if has_i and not has_j:
                    dominated.add(jj)
                elif has_j and not has_i:
                    dominated.add(ii)
                # both have TSD or both lack TSD → no preference, fall through

        survivors = [cluster[i] for i in range(len(cluster)) if i not in dominated]
        if len(survivors) == 1:
            out_handle.write(survivors[0]["line"])
            return

        # --- Layer 1.5: Containment / nesting detection + diagnostics ---
        # When one candidate fully contains another, keep BOTH and print
        # diagnostic metadata.  Resolution deferred to a future version.
        def _matching_bases(rec):
            return float(rec["aln"]) * (1.0 - float(rec["p"]))

        def _format_domains(te_key):
            """Format GFF3 domain entries for a TE as a concise string."""
            if gff3_domains is None:
                return "no_gff3"
            entries = gff3_domains.get(te_key, [])
            if not entries:
                return "none"
            parts = []
            for d in entries:
                parts.append(f"{d['gene']}|{d['clade']}@{d['start']}-{d['end']}")
            return " ".join(parts)

        def _outer_domains_inside_inner(outer_key, inner_s, inner_e):
            """Count how many of the outer's GFF3 domains fall within inner bounds."""
            if gff3_domains is None:
                return "no_gff3", 0, 0
            entries = gff3_domains.get(outer_key, [])
            if not entries:
                return "no_domains", 0, 0
            inside = []
            total = len(entries)
            for d in entries:
                if d["start"] >= inner_s and d["end"] <= inner_e:
                    inside.append(d["gene"])
            return " ".join(inside) if inside else "none", len(inside), total

        nest_involved: Set[int] = set()
        n_surv = len(survivors)
        if n_surv >= 2:
            for ii in range(n_surv):
                for jj in range(ii + 1, n_surv):
                    ri, rj = survivors[ii], survivors[jj]
                    si, ei = int(ri['s']), int(ri['e'])
                    sj, ej = int(rj['s']), int(rj['e'])
                    containment = _is_contained((si, ei), (sj, ej))
                    if not containment:
                        continue

                    if containment == "a_in_b":
                        inner_rec, outer_rec = ri, rj
                        inner_s, inner_e = si, ei
                        outer_s, outer_e = sj, ej
                    else:
                        inner_rec, outer_rec = rj, ri
                        inner_s, inner_e = sj, ej
                        outer_s, outer_e = si, ei

                    inner_key = f"{inner_rec['chrom']}:{inner_s}-{inner_e}"
                    outer_key = f"{outer_rec['chrom']}:{outer_s}-{outer_e}"

                    # True nesting requires distinct LTR pairs
                    if _ltrs_shared(inner_key, outer_key, ltr_bounds):
                        continue

                    nest_involved.add(ii)
                    nest_involved.add(jj)

                    inner_tsd = "yes" if (tsd_names and inner_key in tsd_names) else "no"
                    outer_tsd = "yes" if (tsd_names and outer_key in tsd_names) else "no"
                    inner_mb = _matching_bases(inner_rec)
                    outer_mb = _matching_bases(outer_rec)

                    inside_genes, n_inside, n_total = _outer_domains_inside_inner(
                        outer_key, inner_s, inner_e)

                    print(f"[dedup] NESTED_PAIR:")
                    print(f"  outer: {outer_key} | aln={outer_rec['aln']} "
                          f"p={float(outer_rec['p']):.6f} mbases={outer_mb:.1f} | "
                          f"TSD={outer_tsd}")
                    print(f"    domains: {_format_domains(outer_key)}")
                    print(f"  inner: {inner_key} | aln={inner_rec['aln']} "
                          f"p={float(inner_rec['p']):.6f} mbases={inner_mb:.1f} | "
                          f"TSD={inner_tsd}")
                    print(f"    domains: {_format_domains(inner_key)}")
                    print(f"  outer_domains_inside_inner: "
                          f"{n_inside}/{n_total} ({inside_genes})")
                    print(f"  BOTH_RETAINED")

        if nest_involved:
            # Write all elements that participated in containment pairs.
            for i in sorted(nest_involved):
                out_handle.write(survivors[i]["line"])
            # Also write any non-involved survivors through Layer 2 below.
            remaining = [survivors[i] for i in range(n_surv)
                         if i not in nest_involved]
            if not remaining:
                return
            survivors = remaining
            if len(survivors) == 1:
                out_handle.write(survivors[0]["line"])
                return

        # --- Layer 2: Matching-bases score with ratio guard ---
        RATIO_CAP = 2.5

        best_score = max(survivors, key=lambda x: _matching_bases(x))
        best_pdist = min(survivors, key=lambda x: (float(x["p"]), -int(x["aln"])))

        if best_score is best_pdist:
            best = best_score
        else:
            p_hi  = max(float(best_score["p"]), 1e-10)
            p_lo  = max(float(best_pdist["p"]), 1e-10)
            aln_hi = max(int(best_score["aln"]), 1)
            aln_lo = max(int(best_pdist["aln"]), 1)

            if p_hi / p_lo > RATIO_CAP * (aln_hi / aln_lo):
                best = best_pdist
            else:
                best = best_score

        out_handle.write(best["line"])

    with open(tmp_sorted, "r") as fin, open(out_tsv, "w") as out:
        for raw in fin:
            if not raw.strip():
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 7:
                out.write(raw)
                continue

            parsed = _parse_interval_from_kmer2ltr_col1(parts[0])
            if parsed is None:
                out.write(raw)
                continue
            chrom, s, e = parsed

            try:
                p = float(parts[6])
            except ValueError:
                out.write(raw)
                continue

            try:
                aln = int(float(parts[2]))
            except ValueError:
                aln = 0

            rec = {"chrom": chrom, "s": s, "e": e, "p": p, "aln": aln, "line": raw}

            if not cluster:
                cluster = [rec]
                cluster_chrom = chrom
                continue

            if chrom != cluster_chrom:
                flush_cluster(out)
                cluster = [rec]
                cluster_chrom = chrom
                continue

            belongs = False
            for rrec in cluster:
                frac = _overlap_fraction_of_shorter((s, e), (int(rrec["s"]), int(rrec["e"])))  # type: ignore
                if frac >= threshold:
                    belongs = True
                    break

            if belongs:
                cluster.append(rec)
            else:
                flush_cluster(out)
                cluster = [rec]
                cluster_chrom = chrom

        flush_cluster(out)

    try:
        tmp_sorted.unlink()
    except Exception:
        pass


def subset_fasta_by_name_set(in_fa: str, out_fa: str, keep_names: set) -> None:
    """
    Keep FASTA records whose header token (up to whitespace, without '>') is in keep_names.
    """
    n_written = 0
    with open(out_fa, "w") as out:
        for name, seq in iter_fasta(in_fa):
            if name in keep_names:
                out.write(f">{name}\n")
                for i in range(0, len(seq), 60):
                    out.write(seq[i:i+60] + "\n")
                n_written += 1
    if n_written == 0:
        Path(out_fa).touch()


def names_from_kmer2ltr_dedup(dedup_tsv: str) -> set:
    keep = set()
    with open(dedup_tsv, "r") as f:
        for line in f:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if cols:
                keep.add(cols[0])
    return keep


def build_tesorter_full_length_ltr_fasta_from_cls_tsv(cls_tsv_path: str, genome_fa: str, out_fa: str):
    """
    Build full-length LTR-RT FASTA using genome + TEsorter cls.tsv.
    """
    genome = load_fasta_as_dict(genome_fa)

    def parse_te_interval(te: str) -> Optional[Tuple[str, int, int]]:
        if ":" not in te:
            return None
        chrom, rng = te.rsplit(":", 1)
        m = re.match(r"^(\d+)-(\d+)$", rng)
        if not m:
            return None
        s1 = int(m.group(1))
        e1 = int(m.group(2))
        if e1 < s1:
            s1, e1 = e1, s1
        return chrom, s1, e1

    n_written = 0
    with open(cls_tsv_path, "r") as fin, open(out_fa, "w") as out:
        header_seen = False
        for line in fin:
            line = line.rstrip("\n")
            if not line:
                continue

            if not header_seen:
                header_seen = True
                if line.startswith("TE\t") or line.lower().startswith("te\t"):
                    continue

            cols = line.split("\t")
            if len(cols) < 7:
                continue

            te, order, superfam, clade, _complete, strand, _domains = cols[:7]
            if order != "LTR":
                continue

            parsed = parse_te_interval(te)
            if parsed is None:
                continue
            chrom, s1, e1 = parsed

            seq = genome.get(chrom)
            if seq is None:
                continue

            start0 = s1 - 1
            end0 = e1
            if start0 < 0 or end0 > len(seq):
                continue

            frag = seq[start0:end0]
            if strand == "-":
                frag = revcomp(frag)

            hdr = f"{te}#LTR/{superfam}/{clade}"
            out.write(f">{hdr}\n")
            for i in range(0, len(frag), 60):
                out.write(frag[i:i+60] + "\n")
            n_written += 1

    if n_written == 0:
        Path(out_fa).touch()


def _has_exact_tsd(left: str, right: str, min_len: int = 5) -> bool:
    """
    Returns True if left and right share ANY exact substring of length >= min_len,
    BUT rejects low-complexity kmers (must contain >=2 distinct bases).
    """
    left = left.upper()
    right = right.upper()

    if len(left) < min_len or len(right) < min_len:
        return False

    maxk = min(len(left), len(right))

    for k in range(maxk, min_len - 1, -1):
        for i in range(0, len(left) - k + 1):
            mer = left[i:i+k]

            if len(set(mer)) < 2:
                continue

            if mer in right:
                return True

    return False

def wfa_guided_tsd_names(
    kmer2ltr_tsv: str,
    genome_fa: str,
    existing_tsd_names: Set[str],
    min_len: int = 5,
) -> Set[str]:
    """
    For kmer2ltr candidates WITHOUT an existing TSD, use WFA left/right trim
    columns (cols 13-14, 1-indexed) to adjust element boundaries inward and
    retry TSD search at the shifted positions.

    The WFA -k trim tells us how many alignment columns at each end of the
    LTR-pair alignment are unreliable (no k consecutive matches).  This
    corresponds to boundary over-extension:
      - left_trim  → 5' LTR outer boundary shifted left into flanking DNA
      - right_trim → 3' LTR outer boundary shifted right into flanking DNA

    Returns set of newly recovered 'chrom:start-end' keys.
    """
    tsv_path = Path(kmer2ltr_tsv)
    if not tsv_path.exists() or tsv_path.stat().st_size == 0:
        return set()

    genome = load_fasta_as_dict(genome_fa)
    recovered: Set[str] = set()

    with open(kmer2ltr_tsv) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 14:
                continue

            parsed = _parse_interval_from_kmer2ltr_col1(parts[0])
            if parsed is None:
                continue
            chrom, s, e = parsed
            key = f"{chrom}:{s}-{e}"

            if key in existing_tsd_names:
                continue

            try:
                left_trim = int(parts[12])
                right_trim = int(parts[13])
            except (ValueError, IndexError):
                continue

            if left_trim == 0 and right_trim == 0:
                continue

            seq = genome.get(chrom)
            if seq is None:
                continue

            adj_s = s + left_trim   # shift 5' outer boundary inward (1-based)
            adj_e = e - right_trim  # shift 3' outer boundary inward (1-based)

            if adj_s >= adj_e or adj_s < 1 or adj_e > len(seq):
                continue

            # Extract flanks at adjusted boundaries (same window as tsd_positive_full_length_from_scn)
            fl0 = adj_s - 1  # 0-based left boundary
            fr0 = adj_e      # 0-based one-past-end

            left_start0  = max(0, fl0 - 6)
            left_end0    = min(len(seq), fl0 + 1)
            right_start0 = max(0, fr0 - 1)
            right_end0   = min(len(seq), (fr0 - 1) + 7)

            left_flank  = seq[left_start0:left_end0]
            right_flank = seq[right_start0:right_end0]

            if _has_exact_tsd(left_flank, right_flank, min_len=min_len):
                recovered.add(key)

    return recovered


def rescue_nonautonomous_by_tsd_from_scn(
    stitched_scn: str,
    genome_fa: str,
    tesorter_retained_te_keys: set,
    out_fa: str,
    min_len: int = 5,
    require_run_chars: Optional[List[str]] = None,
    exclude_run_char: Optional[str] = None,
    base_min: int = 800,
    flank_min: int = 80,
) -> int:
    """
    For each SCN candidate NOT present in tesorter_retained_te_keys, scan for TSD.
    Writes full-length rescued sequences.
    Returns number rescued.
    """
    genome = load_fasta_as_dict(genome_fa)

    rescued = 0
    seen = set()

    with open(out_fa, "w") as out, open(stitched_scn, "r") as fin:
        for line in fin:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = re.split(r"\s+", line)
            if len(parts) < 12:
                continue

            sret_s, eret_s = parts[0], parts[1]
            chrom = parts[-1]
            if not (sret_s.isdigit() and eret_s.isdigit()):
                continue

            sret1 = int(sret_s)
            eret1 = int(eret_s)
            if eret1 < sret1:
                sret1, eret1 = eret1, sret1

            te_key = f"{chrom}:{sret1}-{eret1}"
            if te_key in tesorter_retained_te_keys:
                continue
            if te_key in seen:
                continue
            seen.add(te_key)

            seq = genome.get(chrom)
            if seq is None:
                continue

            fl0 = sret1 - 1
            fr0 = eret1
            if fl0 < 0 or fr0 > len(seq):
                continue

            left_start0 = max(0, fl0 - 6)
            left_end0   = min(len(seq), fl0 + 1)
            right_start0 = max(0, fr0 - 1)
            right_end0   = min(len(seq), (fr0 - 1) + 7)

            left = seq[left_start0:left_end0]
            right = seq[right_start0:right_end0]

            if not _has_exact_tsd(left, right, min_len=min_len):
                continue

            frag = seq[fl0:fr0]

            req = _is_valid_require_run_chars(require_run_chars)
            if req:
                if not has_nested_run_signature(frag.upper(), req, base_min=base_min, flank_min=flank_min):
                    continue

            if exclude_run_char:
                if has_exclude_run_char(frag.upper(), exclude_run_char):
                    continue

            hdr = f"{te_key}#LTR/unknown/unknown"
            out.write(f">{hdr}\n")
            for i in range(0, len(frag), 60):
                out.write(frag[i:i+60] + "\n")
            rescued += 1

    if rescued == 0:
        Path(out_fa).touch()
    return rescued

def tsd_positive_full_length_from_scn(
    stitched_scn: str,
    genome_fa: str,
    out_fa: str,
    min_len: int = 5,
    require_run_chars: Optional[List[str]] = None,
    exclude_run_char: Optional[str] = None,
    base_min: int = 800,
    flank_min: int = 80,
) -> int:
    """
    Scan ALL SCN candidates for a TSD near intact boundaries and write those candidates
    as full-length FASTA suitable for TEsorter --pass2-classified-fasta.
    Returns number written.
    """
    genome = load_fasta_as_dict(genome_fa)

    written = 0
    seen = set()

    req = _is_valid_require_run_chars(require_run_chars)

    with open(out_fa, "w") as out, open(stitched_scn, "r") as fin:
        for line in fin:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = re.split(r"\s+", line)
            if len(parts) < 12:
                continue

            sret_s, eret_s = parts[0], parts[1]
            chrom = parts[-1]
            if not (sret_s.isdigit() and eret_s.isdigit()):
                continue

            sret1 = int(sret_s)
            eret1 = int(eret_s)
            if eret1 < sret1:
                sret1, eret1 = eret1, sret1

            te_key = f"{chrom}:{sret1}-{eret1}"
            if te_key in seen:
                continue
            seen.add(te_key)

            seq = genome.get(chrom)
            if seq is None:
                continue

            fl0 = sret1 - 1
            fr0 = eret1
            if fl0 < 0 or fr0 > len(seq):
                continue

            left_start0  = max(0, fl0 - 6)
            left_end0    = min(len(seq), fl0 + 1)
            right_start0 = max(0, fr0 - 1)
            right_end0   = min(len(seq), (fr0 - 1) + 7)

            left = seq[left_start0:left_end0]
            right = seq[right_start0:right_end0]

            if not _has_exact_tsd(left, right, min_len=min_len):
                continue

            frag = seq[fl0:fr0]

            if req:
                if not has_nested_run_signature(frag.upper(), req, base_min=base_min, flank_min=flank_min):
                    continue

            if exclude_run_char:
                if has_exclude_run_char(frag.upper(), exclude_run_char):
                    continue

            hdr = f"{te_key}#LTR/unknown/unknown"
            out.write(f">{hdr}\n")
            for i in range(0, len(frag), 60):
                out.write(frag[i:i+60] + "\n")
            written += 1

    if written == 0:
        Path(out_fa).touch()
    return written


# -----------------------------
# main
# -----------------------------

def main():
    ap = argparse.ArgumentParser(
        description="miniprot genic mask + minimap2 non-LTR mask + chunk + ltrharvest parallel + stitch (GFF3+SCN) + stitched FASTA",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    ap.add_argument("--genome", required=True, help="Input genome FASTA")
    ap.add_argument("--proteins", default=None, help="Protein FASTA for miniprot (optional; if omitted, skip gene masking)")
    ap.add_argument("--te-library", default=None, help="TE library FASTA (optional; if omitted, skip non-LTR TE masking)")
    ap.add_argument("--out-prefix", required=True, help="Output prefix")
    ap.add_argument("--threads", type=int, default=8, help="Threads for miniprot/minimap2 and chunk parallelism")
    ap.add_argument("--workdir", default=None, help="Working directory (default: {out_prefix}.work)")
    ap.add_argument("--gt", default="gt", help="Path to GenomeTools 'gt' executable")

    ap.add_argument("--clean", action="store_true", help="Remove {out_prefix}.work and ./tools after success")

    ap.add_argument("--verbose", action="store_true",
                    help="Print subcommands and additional progress details to stderr")

    # miniprot tuning
    ap.add_argument("--outn", type=int, default=1000, help="miniprot --outn")
    ap.add_argument("--outs", type=float, default=0.99, help="miniprot --outs")
    ap.add_argument("--outc", type=float, default=0.1, help="miniprot --outc")
    ap.add_argument("--gene-mask", default="CDS", choices=["mRNA", "CDS", "both"], help="Which miniprot GFF features to mask in the genome")

    # minimap2 defaults
    ap.add_argument("--mm2-p", type=float, default=0.0, help="minimap2 -p")
    ap.add_argument("--mm2-N", type=int, default=5000, help="minimap2 -N")
    ap.add_argument("--mm2-k", type=int, default=13, help="minimap2 -k")
    ap.add_argument("--mm2-w", type=int, default=5, help="minimap2 -w")

    # PAF filter defaults
    ap.add_argument("--seq-ident", type=float, default=0.60, help="Min identity (matches/aln_len) for masking hits")
    ap.add_argument("--aln-len", type=int, default=10, help="Min alignment length for masking hits")
    ap.add_argument("--qcov", type=float, default=0.01, help="Min query coverage for masking hits")
    
    # chunking
    ap.add_argument("--size", type=int, default=5_000_000, help="Chunk size (bp)")
    ap.add_argument("--overlap", type=int, default=30_000, help="Chunk overlap (bp)")

    # SCN size-based filtering
    ap.add_argument("--scn-min-ltr-len", type=int, default=0,
                    help="Min LTR length filter (exclude if l(lLTR) or l(rLTR) < this; 0 disables)")
    ap.add_argument("--scn-min-ret-len", type=int, default=0,
                    help="Min LTR-RT length filter on l(ret); 0 disables")
    ap.add_argument("--scn-max-ret-len", type=int, default=0,
                    help="Max LTR-RT length filter on l(ret); 0 disables")
    ap.add_argument("--scn-min-int-len", type=int, default=0,
                    help="Min internal length filter on (l(ret)-l(lLTR)-l(rLTR)); 0 disables")
    ap.add_argument("--scn-max-int-len", type=int, default=0,
                    help="Max internal length filter on (l(ret)-l(lLTR)-l(rLTR)); 0 disables")

    # Kmer2LTR + dedup controls
    ap.add_argument("--kmer2ltr-max-win-overdisp", type=float, default=1000.0,
                    help="Kmer2LTR --max-win-overdisp")
    ap.add_argument("--kmer2ltr-min-retained-fraction", type=float, default=0.01,
                    help="Kmer2LTR --min-retained-fraction")
    ap.add_argument("--dedup-threshold", type=float, default=0.80,
                    help="Overlap threshold (fraction of shorter interval) for dedup (default: 0.80)")
    ap.add_argument("--wfa-align", action="store_true", dest="wfa_align",
                    help="Pass --wfa-align to Kmer2LTR: use WFA instead of mafft for pairwise LTR alignment (~30-50x faster)")
    
    ap.add_argument(
        "--kmer2ltr-domains",
        dest="kmer2ltr_domains",
        action="store_true",
        default=True,
        help="Pass {out_prefix}.kmer2ltr.domain to Kmer2LTR via -D"
    )
    ap.add_argument(
        "--no-kmer2ltr-domains",
        dest="kmer2ltr_domains",
        action="store_false",
        help="Do not pass a domain file to Kmer2LTR"
    )

    # TEsorter
    ap.add_argument(
        "--tesorter",
        dest="use_tesorter",
        action="store_true",
        default=True,
        help="Use TEsorter filtering + cls annotation (default: enabled)."
    )
    ap.add_argument(
        "--no-tesorter",
        dest="use_tesorter",
        action="store_false",
        help="Skip TEsorter entirely; build intact FASTA from SCN and feed directly to Kmer2LTR."
    )

    ap.add_argument(
        "--tesorter-use-ret",
        action="store_true",
        help="For the FASTA fed into TEsorter, extract full-length LTR-RT using s(ret)/e(ret) instead of internal (default: internal)."
    )

    ap.add_argument("--tesorter-db", default="rexdb-plant",
                    help="TEsorter HMM database (-db)")
    ap.add_argument("--tesorter-cov", type=int, default=20,
                    help="TEsorter min coverage (-cov)")
    ap.add_argument("--tesorter-eval", default="1e-2",
                    help="TEsorter max E-value (-eval)")
    ap.add_argument("--tesorter-rule", default="80-80-80",
                    help="TEsorter pass2 rule (-rule)")
    ap.add_argument(
        "--tesorter-no-cleanup",
        action="store_true",
        help="Pass --no-cleanup to TEsorter (keep intermediate files)."
    )

    ap.add_argument(
        "--pass2-classified-fasta",
        default=None,
        help="Optional FASTA of previously-classified elements to augment TEsorter pass-2 database. "
             "Headers must be like: >id#Order/Superfamily/Clade"
    )

    ap.add_argument(
        "--trf",
        dest="use_trf",
        action="store_true",
        default=True,
        help="Run TRF-mod on each chunk and hardmask low-complexity sequence before LTR calling (default: enabled)."
    )
    ap.add_argument(
        "--no-trf",
        dest="use_trf",
        action="store_false",
        help="Disable TRF-mod chunk masking."
    )
    ap.add_argument(
        "--trf-args",
        default="-p 7 -s 35",
        help="Quoted args appended to trf-mod (default: '-p 7 -s 35')."
    )
    ap.add_argument(
        "--trf-timeout",
        type=int,
        default=0,
        help="Max seconds allowed per chunk for trf-mod (0 = no timeout)."
    )
    ap.add_argument(
        "--trf-min-copy",
        type=float,
        default=0.0,
        help="Only retain TRF-mod hits with copyNum >= this value (default: 0 = keep all)."
    )

    ap.add_argument(
        "--tsd-pass2",
        action="store_true",
        help="Scan all SCN candidates for TSDs BEFORE TEsorter and feed TSD+ full-length sequences into TEsorter via --pass2-classified-fasta (default: off)."
    )
    ap.add_argument(
        "--tsd-min-len",
        type=int,
        default=5,
        help="Minimum exact shared k-mer length to call a TSD (default: 5)."
    )

    ap.add_argument("--ltr-tools", default="both", choices=["both", "ltrharvest", "ltr_finder"],
                    help="Which LTR caller(s) to run per chunk")

    ap.add_argument("--ltrfinder", default="ltr_finder",
                    help="Path to ltr_finder v1.07 executable")

    ap.add_argument(
        "--require-run-chars",
        default=None,
        help="Comma-separated chars that must each occur as a run in the FULL LTR-RT "
             "(e.g. 'Y' or 'Y,R'). Optional."
    )
    ap.add_argument("--nested-base-min", type=int, default=800,
                    help="Min length for the first char in --require-run-chars (default: 800)")
    ap.add_argument("--nested-flank-min", type=int, default=80,
                    help="Min length for each flanking nesting char (default: 80)")

    ap.add_argument(
        "--exclude-run-char",
        default=None,
        help="Single character: exclude candidate LTR-RTs whose full-length sequence contains "
             "a run of >=10 of this character anywhere (e.g. '--exclude-run-char V'). Optional."
    )

    ap.add_argument(
        "--ltr-timeout",
        type=int,
        default=0,
        help="Max seconds allowed per chunk per LTR tool call (0 = no timeout). Timed-out chunks salvage partial outputs."
    )

    ap.add_argument(
        "--ltrfinder-args",
        default="-w 2 -C -D 15000 -d 100 -L 7000 -l 100 -p 20 -M 0.00 -S 0.0", 
        help="Quoted string of args appended to ltr_finder (NOTE: -w 2 is required)"
    )

    ap.add_argument(
        "--ltrharvest-args",
        default="-mindistltr 100 -minlenltr 100 -maxlenltr 7000 -mintsd 0 -maxtsd 0 -similar 70 -vic 60 -seed 15 -seqids yes -xdrop 10",
        help="Quoted string of args appended to 'gt ltrharvest -index <idx>'"
    )

    args = ap.parse_args()
    out_prefix = args.out_prefix
    verbose = args.verbose

    # Validate --exclude-run-char
    exclude_run_char = None
    if args.exclude_run_char:
        c = args.exclude_run_char.strip().upper()
        if len(c) != 1:
            ap.error("--exclude-run-char expects a single character")
        exclude_run_char = c

    workdir = Path(args.workdir or f"{out_prefix}.work")
    mkdirp(workdir)

    tools_dir = Path("./tools")
    minimap2_path, miniprot_path = ensure_tools(tools_dir)
    tesorter_py_path = None
    if args.use_tesorter:
        tesorter_py_path = ensure_tesorter(tools_dir)

    kmer2ltr_py = ensure_kmer2ltr(tools_dir)
    trfmod_path = None
 
    trf_args = None
    if args.use_trf:
        trfmod_path = ensure_trfmod(tools_dir)
        trf_args = args.trf_args.strip().split()

    if args.proteins or args.te_library:
        which_or_die("bedtools")

    which_or_die(args.gt)
    if args.ltr_tools in ("both", "ltr_finder"):
        which_or_die(args.ltrfinder)

    lf_args = args.ltrfinder_args.strip().split()
    if args.ltr_tools in ("both", "ltr_finder"):
        if "-w" not in lf_args:
            raise RuntimeError("ltr_finder requires '-w 2' output format; please include '-w 2' in --ltrfinder-args")
        for i, tok in enumerate(lf_args):
            if tok == "-w":
                if i + 1 >= len(lf_args) or lf_args[i + 1] != "2":
                    raise RuntimeError("ltr_finder output format must be '-w 2' for SCN conversion compatibility")

    if verbose:
        print("\n[Config]", file=sys.stderr)
        print(f"  genome:        {args.genome}", file=sys.stderr)
        print(f"  out_prefix:    {out_prefix}", file=sys.stderr)
        print(f"  workdir:       {workdir}", file=sys.stderr)
        print(f"  threads:       {args.threads}", file=sys.stderr)
        print(f"  ltr_tools:     {args.ltr_tools}", file=sys.stderr)
        print(f"  tesorter:      {args.use_tesorter}", file=sys.stderr)
        print(f"  trf masking:   {args.use_trf}", file=sys.stderr)
        if exclude_run_char:
            print(f"  exclude_run:   {exclude_run_char} (>=10bp run)", file=sys.stderr)
        print("", file=sys.stderr)

    # Timing accumulators for core modules
    _t_miniprot = 0.0
    _t_ltr_annotation = 0.0
    _t_tesorter = 0.0
    _t_kmer2ltr = 0.0
    _t_total_start = time.monotonic()

    # Step 2
    _t0 = time.monotonic()
    if args.proteins:
        print("[Step2] miniprot gene masking...")
        gff_path, genic_bed, genic_masked_fa = run_miniprot_gene_mask(
            in_fasta=args.genome,
            protein_faa=args.proteins,
            out_prefix=str(workdir / out_prefix),
            threads=args.threads,
            miniprot_path=miniprot_path,
            outn=args.outn,
            outs=args.outs,
            outc=args.outc,
            gene_mask=args.gene_mask,
            verbose=verbose,
        )
    else:
        print("[Step2] skipping gene masking (no --proteins provided).")
        gff_path = str(workdir / f"{out_prefix}.genic.gff")
        genic_bed = str(workdir / f"{out_prefix}.genic.mask.bed")
        genic_masked_fa = str(workdir / f"{out_prefix}.genic_masked.fa")
        Path(gff_path).write_text("")
        Path(genic_bed).write_text("")
        shutil.copyfile(args.genome, genic_masked_fa)


    # Step 3
    if args.te_library:
        print("[Step3] minimap2 non-LTR masking...")
        nonltr_fa = str(workdir / f"{out_prefix}.nonLTR_TE_library.fa")
        write_nonltr_library(args.te_library, nonltr_fa)

        paf_path = str(workdir / f"{out_prefix}.nonLTR_vs_genome.paf")
        cmd_mm2 = [
            minimap2_path,
            "-c", "--cs=short",
            "-t", str(args.threads),
            "--secondary=yes",
            "-p", str(args.mm2_p),
            "-N", str(args.mm2_N),
            f"-k{args.mm2_k}",
            f"-w{args.mm2_w}",
            genic_masked_fa,
            nonltr_fa,
        ]
        if verbose:
            print(f"  $ {' '.join(cmd_mm2)}", file=sys.stderr)
        r = run(cmd_mm2)
        if r.returncode != 0:
            raise RuntimeError(f"minimap2 failed:\n{(r.stderr or '').strip()}")
        Path(paf_path).write_text(r.stdout or "")

        repeats_bed = str(workdir / f"{out_prefix}.repeats.bed")
        paf_to_bed_filtered(paf_path, repeats_bed, args.seq_ident, args.aln_len, args.qcov)

        merged_bed = str(workdir / f"{out_prefix}.repeats.merged.bed")
        r = run(["bedtools", "merge", "-i", repeats_bed], check=True)
        Path(merged_bed).write_text(r.stdout or "")

        gene_te_masked_fa = str(workdir / f"{out_prefix}.gene_TEmasked.fa")
        run(["bedtools", "maskfasta", "-fi", genic_masked_fa, "-bed", merged_bed, "-fo", gene_te_masked_fa],
        check=True, capture=True)
    else:
        print("[Step3] skipping non-LTR TE masking (no --te-library provided).")
        gene_te_masked_fa = str(workdir / f"{out_prefix}.gene_TEmasked.fa")
        shutil.copyfile(genic_masked_fa, gene_te_masked_fa)

    _t_miniprot = time.monotonic() - _t0

    # Always purge genic_masked.fa; purge other masking intermediates with --clean
    _clean_path(genic_masked_fa, verbose)
    if args.clean:
        for _p in [gff_path, genic_bed]:
            _clean_path(_p, verbose)
        if args.te_library:
            for _p in [nonltr_fa, paf_path, repeats_bed, merged_bed]:
                _clean_path(_p, verbose)

    # Step 4
    _t0 = time.monotonic()
    print("[Step4] chunking masked genome...")
    chunks_dir = str(workdir / "chunks")
    chunks = make_chunks(gene_te_masked_fa, chunks_dir, args.size, args.overlap)
    if not chunks:
        raise RuntimeError("No chunks produced (empty genome?)")
    if verbose:
        print(f"  {len(chunks)} chunks (size={args.size:,}, overlap={args.overlap:,})", file=sys.stderr)

    # Masked genome is now chunked; no longer needed
    _clean_path(gene_te_masked_fa, verbose)

    # Step 5
    print(f"[Step5] LTR calling on {len(chunks)} chunks with {args.threads} workers...")
    ltr_args = args.ltrharvest_args.strip().split()
    lf_args = args.ltrfinder_args.strip().split()

    ltr_work = str(workdir / "ltr_runs")
    mkdirp(Path(ltr_work))

    harvest_pairs: List[Tuple[ChunkInfo, str]] = []
    harvest_gffs: List[Tuple[ChunkInfo, str, str]] = []
    finder_pairs: List[Tuple[ChunkInfo, str]] = []

    failures = 0
    total = len(chunks)
    completed = 0

    start = time.monotonic()
    last_update = 0.0
    update_interval = 1.0  # seconds

    print("", file=sys.stderr)

    with ThreadPoolExecutor(max_workers=max(1, args.threads)) as ex:
        futs = {
            ex.submit(
                process_one_chunk,
                args.gt, args.ltrfinder, c,
                ltr_work, ltr_args, lf_args, args.ltr_tools,
                args.ltr_timeout,
                trfmod_path,
                trf_args,
                args.trf_timeout,
                args.trf_min_copy,
                verbose,
            ): c
            for c in chunks
        }

        for fut in as_completed(futs):
            c = futs[fut]
            completed += 1

            try:
                chunk, h_scn, h_gff, f_scn = fut.result()

                if h_scn:
                    harvest_pairs.append((chunk, h_scn))
                if h_scn and h_gff:
                    harvest_gffs.append((chunk, h_scn, h_gff))
                if f_scn:
                    finder_pairs.append((chunk, f_scn))

            except Exception as e:
                failures += 1
                print(f"\n[WARN] chunk failed: {c.chunk_id} :: {e}", file=sys.stderr)

            # ---- progress display ----
            now = time.monotonic()
            if (now - last_update) >= update_interval or completed == total:
                pct = (completed / total) * 100 if total else 100.0
                msg = f"\r[Step5] LTR calling: {completed}/{total} ({pct:.2f}%)"
                print(msg, end="", file=sys.stderr, flush=True)
                last_update = now

    print("", file=sys.stderr)

    if failures:
        print(f"[WARN] {failures} chunks failed; stitching continues using successes only.")

    # Chunk FASTAs are no longer needed once LTR calling is done
    _clean_path(chunks_dir, verbose)

    # Step 6 stitched outputs (inside workdir)
    stitched_gff3 = str(workdir / f"{out_prefix}.ltrharvest.stitched.gff3")
    stitched_harvest_scn = str(workdir / f"{out_prefix}.ltrharvest.stitched.scn")
    stitched_finder_scn  = str(workdir / f"{out_prefix}.ltrfinder.stitched.scn")

    # Merged SCN stored inside workdir (edit 6)
    merged_scn = str(workdir / f"{out_prefix}.ltrtools.stitched.scn")

    if args.ltr_tools in ("both", "ltrharvest"):
        print(f"[Step6] stitching GFF3 (ltrharvest) -> {stitched_gff3}")
        stitch_gff3(harvest_gffs, stitched_gff3)

        print(f"[Step6] stitching SCN  (ltrharvest) -> {stitched_harvest_scn}")
        stitch_scn_from_pairs(harvest_pairs, stitched_harvest_scn)
    else:
        Path(stitched_gff3).write_text("##gff-version 3\n")
        Path(stitched_harvest_scn).touch()

    if args.ltr_tools in ("both", "ltr_finder"):
        print(f"[Step6] stitching SCN  (ltr_finder) -> {stitched_finder_scn}")
        stitch_scn_from_pairs(finder_pairs, stitched_finder_scn)
    else:
        Path(stitched_finder_scn).touch()

    to_merge = []
    if args.ltr_tools in ("both", "ltrharvest"):
        to_merge.append(stitched_harvest_scn)
    if args.ltr_tools in ("both", "ltr_finder"):
        to_merge.append(stitched_finder_scn)

    print(f"[Step6] merging stitched SCN(s) -> {merged_scn}")
    merge_stitched_scns(to_merge, merged_scn)

    # Per-chunk run dirs are no longer needed after stitching
    _clean_path(ltr_work, verbose)
    if args.clean:
        for _p in [stitched_harvest_scn, stitched_finder_scn, stitched_gff3]:
            _clean_path(_p, verbose)

    # Optional: size-based filtering of merged SCN
    if any([
        args.scn_min_ltr_len,
        args.scn_min_ret_len,
        args.scn_max_ret_len,
        args.scn_min_int_len,
        args.scn_max_int_len,
    ]):
        filtered_scn = str(workdir / f"{out_prefix}.ltrtools.stitched.filtered.scn")
        print(f"[Step6b] filtering SCN by lengths -> {filtered_scn}")

        stats = filter_scn_by_lengths(
            in_scn=merged_scn,
            out_scn=filtered_scn,
            min_ltr_len=args.scn_min_ltr_len,
            min_ret_len=args.scn_min_ret_len,
            max_ret_len=args.scn_max_ret_len,
            min_int_len=args.scn_min_int_len,
            max_int_len=args.scn_max_int_len,
        )

        print(
            f"[Step6b] SCN filter: kept {stats['kept']}/{stats['total']} | "
            f"minLTR_fail={stats['fail_min_ltr']} retRange_fail={stats['fail_ret_range']} "
            f"intRange_fail={stats['fail_int_range']} malformed={stats['malformed']}"
        )

        merged_scn = filtered_scn
    else:
        print("[Step6b] skipping SCN length filtering (no --scn-* filters set).")

    _t_ltr_annotation = time.monotonic() - _t0

    internals_fa = str(workdir / f"{out_prefix}.ltrtools.internals.fa")
    intact_for_tesorter_fa = str(workdir / f"{out_prefix}.ltrtools.intact_for_tesorter.fa")
    intact_fa = f"{out_prefix}.ltrtools.intact.fa"

    req_chars = None
    _t0 = time.monotonic()
    if args.require_run_chars:
        req_chars = [x.strip() for x in args.require_run_chars.split(",") if x.strip()]

    tsd_pass2_fa: Optional[str] = None  # populated below if --tsd-pass2 is active

    if args.use_tesorter:
        if args.tesorter_use_ret:
            tesorter_in_fa = intact_for_tesorter_fa
            print(f"[Step7] building FULL-LENGTH FASTA (sret/eret) for TEsorter -> {tesorter_in_fa}")
            scn_to_intact_fasta(
                merged_scn,
                args.genome,
                tesorter_in_fa,
                require_run_chars=req_chars,
                exclude_run_char=exclude_run_char,
                base_min=args.nested_base_min,
                flank_min=args.nested_flank_min,
            )
        else:
            tesorter_in_fa = internals_fa
            print(f"[Step7] building INTERNALS FASTA (elLTR/srLTR) for TEsorter -> {tesorter_in_fa}")
            scn_to_internal_fasta(
                merged_scn,
                args.genome,
                tesorter_in_fa,
                require_run_chars=req_chars,
                exclude_run_char=exclude_run_char,
                base_min=args.nested_base_min,
                flank_min=args.nested_flank_min,
            )
    else:
        tesorter_in_fa = None
        print(f"[Step7] building INTACT FASTA from SCN -> {intact_fa}")
        scn_to_intact_fasta(
            merged_scn,
            args.genome,
            intact_fa,
            require_run_chars=req_chars,
            exclude_run_char=exclude_run_char,
        )

    # Load LTR boundaries from SCN for nesting validation
    ltr_bounds = load_scn_ltr_boundaries(merged_scn)

    # Step 9: TEsorter
    cls_tsv_path = None
    tesorter_lib_fa = None

    if args.use_tesorter:
        tesorter_outdir = workdir

        pass2_for_tesorter = args.pass2_classified_fasta

        if args.use_tesorter and args.tsd_pass2:
            tsd_pass2_fa = str(workdir / f"{out_prefix}.tsd_pass2.full_length.fa")
            print(f"[Step8.9] scanning SCN candidates for TSDs -> {tsd_pass2_fa}")

            n_tsd = tsd_positive_full_length_from_scn(
                stitched_scn=merged_scn,
                genome_fa=args.genome,
                out_fa=tsd_pass2_fa,
                min_len=args.tsd_min_len,
                require_run_chars=req_chars,
                exclude_run_char=exclude_run_char,
                base_min=args.nested_base_min,
                flank_min=args.nested_flank_min,
            )
            print(f"[Step8.9] TSD pass2: found {n_tsd} TSD+ candidates")

            merged_pass2 = str(workdir / f"{out_prefix}.pass2_classified.merged.fa")
            merged_path = merge_pass2_fastas(args.pass2_classified_fasta, tsd_pass2_fa, merged_pass2)

            pass2_for_tesorter = merged_path if merged_path else None

            # --- Step 8.9b: pre-purge TSD-dominated candidates ---
            # Non-TSD candidates that overlap a TSD+ candidate at >= dedup
            # threshold would be eliminated by Layer-1 of dedup anyway.
            # Removing them now shrinks the TEsorter + Kmer2LTR input.
            # Contained candidates with distinct LTRs are protected
            # (putative nested TEs).
            tsd_names_early = _tsd_names_from_fasta(tsd_pass2_fa)
            if tsd_names_early:
                purge_set = pre_purge_tsd_dominated(
                    merged_scn, tsd_names_early,
                    threshold=args.dedup_threshold,
                    ltr_bounds=ltr_bounds,
                )
                if purge_set:
                    prepurge_fa = str(workdir / f"{out_prefix}.ltrtools.prepurge.fa")
                    n_before = 0
                    n_after = 0
                    with open(prepurge_fa, "w") as out_fh:
                        for name, seq in iter_fasta(tesorter_in_fa):
                            n_before += 1
                            if name in purge_set:
                                continue
                            n_after += 1
                            out_fh.write(f">{name}\n")
                            for i in range(0, len(seq), 60):
                                out_fh.write(seq[i:i+60] + "\n")
                    print(f"[Step8.9b] pre-purge: {n_before} -> {n_after} candidates "
                          f"({n_before - n_after} TSD-dominated removed)")
                    tesorter_in_fa = prepurge_fa

        print(f"[Step9] running TEsorter on stitched FASTA (outputs -> {tesorter_outdir})")

        cls_lib_path, cls_pep_path, cls_tsv_path = run_tesorter(
            stitched_fa=tesorter_in_fa,
            tesorter_py_path=tesorter_py_path,
            outdir=tesorter_outdir,
            db=args.tesorter_db,
            cov=args.tesorter_cov,
            evalue=args.tesorter_eval,
            rule=args.tesorter_rule,
            threads=args.threads,
            pass2_classified_fasta=pass2_for_tesorter,
            no_cleanup=args.tesorter_no_cleanup,
            verbose=verbose,
        )

        tesorter_lib_fa = str(workdir / f"{out_prefix}.ltrharvest.full_length.fa.{args.tesorter_db}.cls.lib.fa")
        print(f"[Step9] building full-length LTR FASTA from cls.tsv -> {tesorter_lib_fa}")
        build_tesorter_full_length_ltr_fasta_from_cls_tsv(cls_tsv_path, args.genome, tesorter_lib_fa)

    _t_tesorter = time.monotonic() - _t0

    # Step 8: kmer2ltr.domain from stitched SCN
    kmer2ltr_domain = str(workdir / f"{out_prefix}.kmer2ltr.domain")
    print(f"[Step8] building kmer2ltr domain -> {kmer2ltr_domain}")
    scn_to_kmer2ltr_domain(merged_scn, kmer2ltr_domain, tesorter_cls_tsv=cls_tsv_path)

    # Step 9b: Kmer2LTR
    _t0 = time.monotonic()
    print("[Step9b] running Kmer2LTR (dedup driver)...")
    k2l_prefix = f"{out_prefix}_kmer2ltr"

    domain_arg = kmer2ltr_domain if args.kmer2ltr_domains else None

    if args.use_tesorter:
        k2l_in_fa = tesorter_lib_fa
    else:
        k2l_in_fa = intact_fa

    k2l_main = run_kmer2ltr(
        kmer2ltr_py=kmer2ltr_py,
        in_fa=k2l_in_fa,
        out_prefix=k2l_prefix,
        outdir=workdir,
        threads=args.threads,
        max_win_overdisp=args.kmer2ltr_max_win_overdisp,
        min_retained_fraction=args.kmer2ltr_min_retained_fraction,
        domain_file=domain_arg,
        wfa_align=args.wfa_align,
        verbose=verbose,
    )

    k2l_dedup_out = f"{out_prefix}_ltr.tsv"
    print(f"[Step9b] deduping Kmer2LTR output -> {k2l_dedup_out}")
    tsd_name_set = _tsd_names_from_fasta(tsd_pass2_fa) if tsd_pass2_fa else set()

    # WFA-guided TSD recovery: use alignment trim info to find TSDs at adjusted boundaries
    wfa_recovered = wfa_guided_tsd_names(k2l_main, args.genome, tsd_name_set,
                                         min_len=args.tsd_min_len)
    if wfa_recovered:
        print(f"[Step9b] WFA-guided TSD recovery: {len(wfa_recovered)} additional TSDs")
    combined_tsd = tsd_name_set | wfa_recovered if (tsd_name_set or wfa_recovered) else None

    # Load GFF3 domain positions for nesting diagnostics
    gff3_dom_data = None
    if cls_tsv_path:
        gff3_path = cls_tsv_path.replace(".cls.tsv", ".dom.gff3")
        if Path(gff3_path).exists():
            gff3_dom_data = load_tesorter_gff3_domains(gff3_path)
            print(f"[Step9b] loaded GFF3 domain data: "
                  f"{len(gff3_dom_data)} TEs with domain annotations")

    dedup_kmer2ltr_tsv(k2l_main, k2l_dedup_out, threshold=args.dedup_threshold,
                       tsd_names=combined_tsd, gff3_domains=gff3_dom_data,
                       ltr_bounds=ltr_bounds)

    _t_kmer2ltr = time.monotonic() - _t0

    keep_names = names_from_kmer2ltr_dedup(k2l_dedup_out)

    if args.use_tesorter:
        tesorter_lib_fa_dedup = f"{out_prefix}_ltr.fa"
        print(f"[Step9b] subsetting full-length FASTA by dedup list -> {tesorter_lib_fa_dedup}")
        subset_fasta_by_name_set(k2l_in_fa, tesorter_lib_fa_dedup, keep_names)
    else:
        intact_dedup_fa = f"{out_prefix}.ltrtools.intact.dedup.fa"
        print(f"[Step9b] subsetting intact FASTA by dedup list -> {intact_dedup_fa}")
        subset_fasta_by_name_set(intact_fa, intact_dedup_fa, keep_names)

    _t_total = time.monotonic() - _t_total_start
    print("\nDone.")
    print("Runtimes:")
    print(f"  Miniprot/masking:   {_format_hms(_t_miniprot)}")
    print(f"  LTR annotation:    {_format_hms(_t_ltr_annotation)}")
    print(f"  TEsorter:          {_format_hms(_t_tesorter)}")
    print(f"  Kmer2LTR:          {_format_hms(_t_kmer2ltr)}")
    print(f"  Total:             {_format_hms(_t_total)}")
    print("")
    print("Primary outputs:")
    print(f"  LTR TSV:   {out_prefix}_ltr.tsv")

    if args.use_tesorter:
        print(f"  LTR FASTA: {out_prefix}_ltr.fa")
        print("")
        print("Workdir outputs:")
        print(f"  Full-length (pre-dedup) FASTA: {tesorter_lib_fa}")
    else:
        print(f"  Intact dedup FASTA:          {out_prefix}.ltrtools.intact.dedup.fa")
        print(f"  Intact (pre-dedup) FASTA:    {out_prefix}.ltrtools.intact.fa")

    print(f"  SCN (merged):                {workdir}/{out_prefix}.ltrtools.stitched.scn")
    print(f"  Kmer2LTR outputs prefix:     {workdir}/{out_prefix}_kmer2ltr*")

    if args.clean:
        print(f"[CLEAN] removing workdir: {workdir}")
        shutil.rmtree(workdir, ignore_errors=True)
        print(f"[CLEAN] removing tools dir: {tools_dir}")
        shutil.rmtree(tools_dir, ignore_errors=True)
    else:
        # Compress FASTA files in workdir in background
        fa_files = list(workdir.glob("*.fa"))
        if fa_files:
            for fa in fa_files:
                subprocess.Popen(["gzip", str(fa)],
                                 stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            print(f"[COMPRESS] gzipping {len(fa_files)} FASTA file(s) in {workdir}/ (background)")

if __name__ == "__main__":
    main()
