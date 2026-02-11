#!/usr/bin/env python3
"""
(1) Masks genes and non-LTR TEs
(2) Runs ltrharvest and ltr_finder and merges results.
(3) Feeds LTR-RT internal sequence to TEsorter to filter out false positives (plus some true positives).
(4) Runs kmer2ltr on remaining pool.
(5) Purges duplicate LTR-RTs based on kmer2ltr LTR divergence.
NOTE: "--tesorter-tree" currently broken. To fix, we'd simply re-run TEsorter on the remaining deduped full length LTR-RTs.

# I could consider adding DeepTE or CREATE as an alternative to TEsorter although I need to check their speed and efficiency feasiblility. 

# Adds TRF-mod for hardmasking low-complexity sequence. 
# LTR-RTs and their flanking sequence contain low-complexity sequence, thus masking will impact the result. 
# Be conservative and only hardmask those with a high copy number. These are the ones that mostly bog down ltrharvest and can bottleneck the run. 
# I recommend (--trf-min-copy 15) but could play with 5, 10, 20, 25, 30:
python ltrharvest3.py --genome Athal.fa --proteins Athal.pep --threads 200 --out-prefix r1_ltr3 --scn-min-ltr-len 100 --scn-min-ret-len 800 --scn-max-ret-len 15000 --scn-min-int-len 500 --scn-max-int-len 14000 --tsd-rescue --trf --trf-args '-p 100 -s 30' --trf-min-copy 15

Benchmarking suggests including the gene protein file is a good idea while non-LTR TEs is optional: 
| Approach              | TP       | FP     | FN      | Precision | Recall    | F1        |
| --------------------- | -------- | ------ | ------- | --------- | --------- | --------- |
| (1) NO prot file      | 3079     | 145    | 692     | 0.955     | 0.816     | 0.880     |
| **(2) NO TE file**    | **3268** | 101    | **503** | 0.970     | **0.867** | **0.915** |
| (3) NO prot / NO TE   | 3179     | 210    | 592     | 0.938     | 0.843     | 0.888     |
| (4) prot + TE file    | 3160     | **49** | 611     | **0.985** | 0.838     | 0.905     |
| (5) EDTA pass + trunc | 2813     | 124    | 958     | 0.958     | 0.746     | 0.839     |
**** This is a friendlier simulation. Durring burn-in, non-LTR-RTs are all fragmented and LTR-RTs are all intact. Durring simulation, only LTR-RTs are mobile ****


| Approach                                | TP       | FP       | FN     | Precision | Recall    | F1        |
| --------------------------------------- | -------- | -------- | ------ | --------- | --------- | --------- |
| (6) EDTA pass + trunc                   | 1255     | **73**   | 806    | **0.945** | 0.609     | 0.741     |
| (7) NO TEsorter                         | **1968** | 1527     | **93** | 0.563     | **0.955** | 0.709     |
| (8) NO TEsorter + SCN filters           | 1928     | 565      | 133    | 0.773     | 0.935     | **0.846** |
| (9) TEsorter                            | 1689     | 405      | 372    | 0.807     | 0.820     | 0.813     |
| **(10) TEsorter + SCN filters**         | 1684     | 285      | 377    | 0.855     | 0.817     | 0.836     |
**** SCN filters are: --scn-min-ltr-len 100 --scn-min-ret-len 800 --scn-max-ret-len 15000 --scn-min-int-len 500 --scn-max-int-len 12000 ****
**** This is a more complicated simulation.  Durring burn-in, all TEs (LTR and non-LTR) are both fragmented and intact. Durring simulation, all TEs (LTR and non-LTR) are mobile ****

# Suggested run:
# Combine #2 and #10. 
# With real data, I suspect TEsorter is required since ltrharvest and ltrfinder parameters are selected to optimize specificity.
# My PriNTE simulations do not test the impact of low-complexity repeats. 
python ltrharvest.py --genome Athal.fa --proteins TAIR10.pep.fa.gz --threads 20 --out-prefix Athal_ltr --tsd-rescue --scn-min-ltr-len 100 --scn-min-ret-len 800 --scn-max-ret-len 15000 --scn-min-int-len 500 --scn-max-int-len 12000

 
# Goldstandard
perl EDTA/EDTA_raw.pl --genome Athal_chr1.fa --type ltr --threads 10
# Finished in 23min.
# Identifes 39 LTR-RTs.
    DB="repbase_RM.nucl"
    Q="Athal.fa.mod.EDTA.raw/LTR/Athal.fa.mod.pass.list.fa"
    PFX="Athal_vs_repbase"
    blastn -task blastn -query "$Q" -db "$DB" -evalue 1e-10 -max_target_seqs 50 -max_hsps 1 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore qcovhsp' -num_threads 16 > "${PFX}.raw.tsv"
    LC_ALL=C sort -k1,1 -k14,14gr -k13,13g "${PFX}.raw.tsv" | awk 'BEGIN{FS=OFS="\t"} !seen[$1]++' > "${PFX}.best.tsv"
    awk 'BEGIN{FS=OFS="\t"}{qcov=100.0*$4/$9; scov=100.0*$4/$12; print $0, sprintf("%.2f",qcov), sprintf("%.2f",scov)}' "${PFX}.best.tsv" > "${PFX}.best.with_cov.tsv"
    { echo -e "qseqid\tsseqid\tpident\taln_len\tmismatch\tgapopen\tqstart\tqend\tqlen\tsstart\tsend\tslen\tevalue\tbitscore\tqcovhsp\tqcov_pct\tscov_pct"; cat "${PFX}.best.with_cov.tsv"; } > "${PFX}.best.with_cov.header.tsv"
    awk 'BEGIN{FS=OFS="\t"} NR==1 || $2 !~ /LTR/ {print}' "${PFX}.best.with_cov.header.tsv" > "${PFX}.no_LTR_hits.tsv"
    awk 'BEGIN{FS="\t"} NR>1 && $2 !~ /LTR/ {n++} END{print n}' "${PFX}.best.with_cov.header.tsv"
# 37 of the 39 have an LTR-RT best-hit in the repbase lib. 2 are LINE.

# This pipeline.
python ltrharvest.py --genome Athal_chr1.fa --proteins ../PrinTE/data/TAIR10.pep.fa.gz --threads 10 --out-prefix ltrharvest --scn-min-ltr-len 100 --scn-min-ret-len 800 --scn-max-ret-len 15000 --scn-min-int-len 500 --scn-max-int-len 12000 --size 500000
# Finished in 30min.
# Identifes 170 LTR-RTs.
# 146 of the 170 have an LTR-RT best-hit in the repbase lib.
# Most of the non LTR-RT alignments are identified using TEsorter 2-pass (LTR/Gypsy/unknown) and very poorly aligned (1-6% query cov). 
# Previous benchmarking suggest TEsorter 2-pass 80-80-80 was useful for detecting real LTR-RTs and lowering 80-80-80 generated artifacts

### DEVELOPING ###
# Nest inserion discovery.
# Round1: Non-nest.
python ltrharvest.py --genome Osati.fa --proteins Osati.pep --threads 100 --out-prefix Osati_r1 --tsd-rescue --scn-min-ltr-len 100 --scn-min-ret-len 800 --scn-max-ret-len 15000 --scn-min-int-len 500 --scn-max-int-len 12000 --size 500000
# Discovers 4919 LTR-RTs.
# Mask LTR-RTs with N. Only flanking sequence is unmasked. 
python mask_ltr.py --features-fasta Osati_r1.ltrharvest.full_length.dedup.fa.rexdb-plant.cls.lib.fa --genome Osati.fa --feature-character N --far-character V --distance 15000 > Osati_r1.fa


# Round2: 1-level nesting.
# Require LTR-RTs to contain runs of N (ie a nested LTR-RT)
# Expand exceptable max LTR-RT lengths since theyre now nested. 
python ltrharvest.py --require-run-chars N --genome Osati_r1.fa --proteins Osati.pep --threads 100 --out-prefix Osati_r2 --tsd-rescue --scn-min-ltr-len 100 --scn-min-ret-len 1000 --scn-max-ret-len 30000 --scn-min-int-len 500 --scn-max-int-len 28000 --ltrharvest-args '-mindistltr 100 -minlenltr 100 -maxlenltr 7000 -mintsd 0 -maxtsd 0 -similar 70 -vic 60 -seed 15 -seqids yes -xdrop 10 -maxdistltr 30000' --ltrfinder-args '-w 2 -C -D 30000 -d 100 -L 7000 -l 100 -p 20 -M 0.00 -S 0.0'
# Discovers 441 1-level nested LTR-RTs. 
# 1-level LTR-RTs masked with N; 2-level LTR-RTs masked with R. Flanking sequence unmasked.
python mask_ltr.py --features-fasta Osati_r2.ltrharvest.full_length.dedup.fa.rexdb-plant.cls.lib.fa --genome Osati.fa --feature-character R --far-character V --distance 15000 > Osati_r2.fa

# Expand exceptable max LTR-RT lengths since theyre now double 2-level nested. 
# Require LTR-RTs to contain runs of N and R (ie 2-level nesting)
python ltrharvest.py --require-run-chars N,R --genome Osati_r2.fa --proteins Osati.pep --threads 100 --out-prefix Osati_r3 --tsd-rescue --scn-min-ltr-len 100 --scn-min-ret-len 1000 --scn-max-ret-len 45000 --scn-min-int-len 500 --scn-max-int-len 43000 --ltrharvest-args '-mindistltr 100 -minlenltr 100 -maxlenltr 7000 -mintsd 0 -maxtsd 0 -similar 70 -vic 60 -seed 15 -seqids yes -xdrop 10 -maxdistltr 45000' --ltrfinder-args '-w 2 -C -D 45000 -d 100 -L 7000 -l 100 -p 20 -M 0.00 -S 0.0'
# Discovers 75  2-level nested LTR-RTs. 
# 1-level LTR-RTs masked with N; 2-level LTR-RTs masked with R. 3-level LTR-RTs masked with D. Flanking sequence unmasked.
python mask_ltr.py --features-fasta Osati_r2.ltrharvest.full_length.dedup.fa.rexdb-plant.cls.lib.fa --genome Osati.fa --feature-character D --far-character V --distance 15000 > Osati_r3.fa

python ltrharvest.py --require-run-chars N,R,D --genome Osati_r3.fa --proteins Osati.pep --threads 100 --out-prefix Osati_r3 --tsd-rescue --scn-min-ltr-len 100 --scn-min-ret-len 1000 --scn-max-ret-len 55000 --scn-min-int-len 500 --scn-max-int-len 54000 --ltrharvest-args '-mindistltr 100 -minlenltr 100 -maxlenltr 7000 -mintsd 0 -maxtsd 0 -similar 70 -vic 60 -seed 15 -seqids yes -xdrop 10 -maxdistltr 55000' --ltrfinder-args '-w 2 -C -D 55000 -d 100 -L 7000 -l 100 -p 20 -M 0.00 -S 0.0'

Needs benchmarked. It may be too strict to expect TEsoerter proteins on higher nesting orders, but this is an outline. 
 Non-nest K2P    :	0.021463
 1-level nest K2P:	0.031816
 2-level nest K2P:	0.060179
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
from typing import Dict, List, Tuple, Iterable, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed


# -----------------------------
# utilities
# -----------------------------

def run(cmd: List[str], cwd: Optional[str] = None, check: bool = False, capture: bool = True, text: bool = True):
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
                           gene_mask: str = "mRNA"):

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

def build_suffixerator_index(gt_path: str, chunk_fa: str, indexname: str):
    cmd = [
        gt_path, "suffixerator",
        "-db", chunk_fa,
        "-indexname", indexname,
        "-tis", "-suf", "-lcp", "-des", "-ssp", "-sds",
        "-dna",
    ]
    r = run(cmd)
    if r.returncode != 0:
        raise RuntimeError(f"gt suffixerator failed for {chunk_fa}:\n{(r.stderr or '').strip()}")

def run_ltrharvest_scn_and_gff3(
    gt_path: str,
    indexname: str,
    scn_path: str,
    gff3_path: str,
    ltrharvest_args: List[str],
    timeout_s: int = 0
):
    """
    - SCN/tabular output is stdout (with -tabout yes)
    - GFF3 output goes to gff3_path (with -gff3 <file>)
    - If timeout occurs, salvage partial stdout to scn_path and keep any partial gff3 file.
    """
    cmd = [gt_path, "ltrharvest", "-index", indexname] + ltrharvest_args + ["-tabout", "yes", "-gff3", gff3_path]

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
        # Salvage partial stdout (SCN/tabout-like)
        partial = te.stdout or ""
        Path(scn_path).write_text(partial)

        # Keep partial gff3 file if it exists; if not, touch it so stitch step doesn't error
        gp = Path(gff3_path)
        if not gp.exists():
            gp.touch()

        # Do NOT raise — allow pipeline to continue and stitch/merge what we got
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
        )


        # Hardmask using TRF BED (BED coords are 0-based, end-exclusive)
        hardmask_fasta_by_bed0(local_fa, trf_bed, trf_masked_fa)

        # Replace local_fa used downstream
        local_fa = trf_masked_fa


    harvest_scn_path = None
    harvest_gff3_path = None
    finder_scn_path = None

    if ltr_tools in ("both", "ltrharvest"):
        indexname = str(wdir / "idx")
        build_suffixerator_index(gt_path, local_fa, indexname)

        harvest_scn_path = str(wdir / f"{chunk.chunk_id}.ltrharvest.scn")
        harvest_gff3_path = str(wdir / f"{chunk.chunk_id}.ltrharvest.gff3")
        run_ltrharvest_scn_and_gff3(
            gt_path, indexname, harvest_scn_path, harvest_gff3_path, ltrharvest_args,
            timeout_s=ltr_timeout_s
        )


    if ltr_tools in ("both", "ltr_finder"):
        raw_path = str(wdir / f"{chunk.chunk_id}.ltr_finder.raw.scn")
        conv_path = str(wdir / f"{chunk.chunk_id}.ltr_finder.ltrharvest_like.scn")

        run_ltrfinder_raw(ltrfinder_path, local_fa, raw_path, ltrfinder_args, timeout_s=ltr_timeout_s)
        raw_text = Path(raw_path).read_text() if Path(raw_path).exists() else ""
        ltrfinder_w2_to_ltrharvest_scn(raw_text, conv_path)

        finder_scn_path = conv_path

    return chunk, harvest_scn_path, harvest_gff3_path, finder_scn_path

# -----------------------------
# Step 5b: ltr_finder per chunk + convert to ltrharvest-style SCN
# -----------------------------

def run_ltrfinder_raw(ltrfinder_path: str, chunk_fa: str, raw_out_path: str, ltrfinder_args: List[str], timeout_s: int = 0):
    """
    Runs ltr_finder on chunk_fa. Captures stdout to raw_out_path.
    If timeout occurs, salvages partial stdout to raw_out_path.
    NOTE: -w 2 is required (user stated); enforced elsewhere.
    """
    cmd = [ltrfinder_path] + ltrfinder_args + [chunk_fa]

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

    Output columns (space-separated):
      s(ret) e(ret) l(ret) s(lLTR) e(lLTR) l(lLTR) s(rLTR) e(rLTR) l(rLTR) sim(LTRs) seq-nr chr

    We mimic the provided Perl logic, but in Python.
    """
    seq_id = -1
    n_written = 0

    with open(out_scn_path, "w") as out:
        # Header lines start with ## so your stitch_scn() will skip them (it skips '#')
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

            # Perl: s/\[\s+/\[/g;
            line2 = re.sub(r"^\[\s+", "[", line.strip())

            toks = line2.split()
            # Need enough tokens to safely index [1,2,3,4,12,15]
            if len(toks) <= 15:
                continue

            chr_ = toks[1]
            loc = toks[2]
            lens = toks[3]      # "lLTRlen,rLTRlen"
            ltr_len = toks[4]   # full-length length
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

            # Compute LTR coordinates
            eltr = sret + lltr_len - 1
            srtr = eret - rltr_len + 1

            # Basic sanity
            if sret <= 0 or eret <= 0 or eltr < sret or srtr > eret:
                continue

            # ltr_finder similarity can be e.g. "0.85" or "85.0%" depending on build/output.
            # ltrharvest "sim(LTRs)" is typically percent-like. We'll normalize:
            sim_val = None
            try:
                if similarity.endswith("%"):
                    sim_val = float(similarity.rstrip("%"))
                else:
                    sim_val = float(similarity)
                    # If it's 0..1, scale to percent
                    if 0.0 <= sim_val <= 1.0:
                        sim_val *= 100.0
            except Exception:
                continue

            # l(ret) should match provided LTR_len token; fallback to eret-sret+1 if needed
            try:
                lret = int(float(ltr_len))
            except Exception:
                lret = eret - sret + 1

            if seq_id < 0:
                # If input had no '>' lines (unlikely), set to 0
                seq_id = 0

            # Direction is not represented in SCN columns; kept only for your debugging if needed.
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

    Example:
      [N] => N{>=base}
      [N,R] => R{>=flank} N{>=base} R{>=flank}
      [N,R,D] => D{>=flank} R{>=flank} N{>=base} R{>=flank} D{>=flank}

    Adjacency is enforced: the flank runs must directly touch the inner run.
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

    # scan for runs of the base char (chars[0])
    while i < n:
        if s[i] != c0:
            i += 1
            continue

        j = i
        while j < n and s[j] == c0:
            j += 1

        run_len = j - i
        if run_len >= base_min:
            # We have a candidate base run [i, j)
            left_idx = i       # boundary just before base run
            right_idx = j      # boundary just after base run

            ok = True
            # require symmetric adjacent flanks for each nesting level
            for level in range(1, len(chars)):
                c = chars[level]

                # left flank must end exactly at left_idx-1
                l_end = left_idx
                l_start = l_end
                while l_start > 0 and s[l_start - 1] == c:
                    l_start -= 1
                left_run = l_end - l_start
                if left_run < flank_min:
                    ok = False
                    break

                # right flank must start exactly at right_idx
                r_start = right_idx
                r_end = r_start
                while r_end < n and s[r_end] == c:
                    r_end += 1
                right_run = r_end - r_start
                if right_run < flank_min:
                    ok = False
                    break

                # move outward for next level
                left_idx = l_start
                right_idx = r_end

            if ok:
                return True

        # continue scanning after this base run
        i = j

    return False


def stitch_scn(chunk_triplets: List[Tuple[ChunkInfo, str, str]], stitched_out: str):
    """
    SCN/tabular line format (ltrharvest tabout):
      s(ret) e(ret) l(ret) s(lLTR) e(lLTR) l(lLTR) s(rLTR) e(rLTR) l(rLTR) sim(LTRs) seq-nr chrom

    We lift 1-based coordinates by adding chunk.start0 to:
      cols 1,2,4,5,7,8   (0,1,3,4,6,7 in 0-based list)
    Then append chromosome as last column.
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

                    off = chunk.start0  # 0-based offset; ltrharvest coords are 1-based inclusive
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
    triplets = [(c, scn, "") for (c, scn) in pairs]  # dummy gff slot
    stitch_scn(triplets, stitched_out)

def merge_stitched_scns(stitched_scns: List[str], merged_out: str):
    """
    Merge multiple stitched SCN files into one, removing duplicates.

    Dedup key:
      (chrom, s(ret), e(ret), s(lLTR), e(lLTR), s(rLTR), e(rLTR))
    Chrom is last column in your stitched format.
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
                    # indices per your SCN convention
                    key = (
                        chrom,
                        parts[0], parts[1],  # sret, eret
                        parts[3], parts[4],  # s(lLTR), e(lLTR)
                        parts[6], parts[7],  # s(rLTR), e(rLTR)
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

    Expected columns (0-based indices):
      0 s(ret)
      1 e(ret)
      2 l(ret)
      3 s(lLTR)
      4 e(lLTR)
      5 l(lLTR)
      6 s(rLTR)
      7 e(rLTR)
      8 l(rLTR)
      9 sim(LTRs)
      10 seq-nr
      11 chrom   (or more; your stitched adds chrom at end)

    Returns counts for reporting.
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

            # Length fields
            try:
                lret = int(parts[2])
                lltr = int(parts[5])
                rltr = int(parts[8])
            except ValueError:
                counts["malformed"] += 1
                continue

            # (1) Min LTR length: exclude if either LTR is too short
            if min_ltr_len and (lltr < min_ltr_len or rltr < min_ltr_len):
                counts["fail_min_ltr"] += 1
                continue

            # (2) Min/Max LTR-RT (lret)
            if min_ret_len and lret < min_ret_len:
                counts["fail_ret_range"] += 1
                continue
            if max_ret_len and lret > max_ret_len:
                counts["fail_ret_range"] += 1
                continue

            # (3) Min/Max internal length
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
    base_min: int = 800,
    flank_min: int = 80,
):
    """
    Uses e(lLTR) and s(rLTR) (1-based inclusive coords from ltrharvest tabout),
    plus chrom (last column), to extract INTERNAL sequence only (LTRs excluded).

    Optional filter:
      If require_run_chars is provided, the FULL LTR-RT (sret..eret) must contain
      a run of `run_len` of EACH requested character somewhere in the full-length
      sequence. Only then do we extract/write the internal.

    Header stays full-length LTR-RT coords:
      >chr1:s(ret)-e(ret)
    """
    genome = load_fasta_as_dict(genome_fa)

    # Normalize required chars: uppercase, strip whitespace, drop empties
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

            # internal is between LTRs: (e(lLTR)+1) .. (s(rLTR)-1)
            internal_start1 = el1 + 1
            internal_end1   = sr1 - 1
            if internal_end1 < internal_start1:
                continue

            seq = genome.get(chrom)
            if seq is None:
                continue

            # bounds check for full-length first
            fl_start0 = sret1 - 1
            fl_end0   = eret1
            if fl_start0 < 0 or fl_end0 > len(seq):
                continue

            req = _is_valid_require_run_chars(require_run_chars)

            # ---- nested signature filter on FULL LENGTH ----
            if req:
                full_len_seq = seq[fl_start0:fl_end0].upper()
                if not has_nested_run_signature(full_len_seq, req, base_min=base_min, flank_min=flank_min):
                    continue

            # Now extract INTERNAL
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

    # reuse your hardmasker (expects Interval start/end already 0-based half-open)
    hardmask_fasta_by_intervals(in_fasta, intervals, out_fasta, out_bed=str(bp))  # keeps bed as-is
    return n

def run_trfmod(
    trfmod_path: str,
    chunk_fa: str,
    out_bed: str,
    trf_args: List[str],
    timeout_s: int = 0,
    min_copy: float = 0.0,
):
    """
    Runs trf-mod on chunk_fa, filters by copy number, and writes BED0 to out_bed.

    TRF-mod stdout format (as you showed):
      ctg start end period copyNum fracMatch fracGap score entropy pattern

    Assumption (matches classic TRF):
      start/end are 1-based inclusive -> convert to BED0 (start0 = start-1, end0 = end)

    If timeout occurs, salvages partial stdout and still applies filtering.
    """
    cmd = [trfmod_path] + trf_args + [chunk_fa]

    def trf_text_to_bed0(text: str) -> str:
        out_lines = []
        for line in text.splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 5:
                # not a data row
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

            # normalize coordinates
            if end1 < start1:
                start1, end1 = end1, start1

            # TRF is (typically) 1-based inclusive; BED is 0-based end-exclusive
            start0 = start1 - 1
            end0 = end1

            if end0 <= start0:
                continue

            # Keep copyNum as an extra column (harmless; mask reader uses first 3 cols)
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
    base_min: int = 800,
    flank_min: int = 80,
):
    """
    Extract FULL intact LTR-RT sequence using s(ret) and e(ret) (1-based inclusive),
    plus chrom (last column).

    Optional filter:
      If require_run_chars is provided, the FULL LTR-RT (sret..eret) must contain
      a run of `run_len` of EACH requested character somewhere in the full-length
      sequence. Only then do we write the intact.

    Header:
      >chr1:s(ret)-e(ret)
    """
    genome = load_fasta_as_dict(genome_fa)

    # Normalize required chars: uppercase, strip whitespace, drop empties
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
            end0 = eret1  # inclusive -> exclusive
            if start0 < 0 or end0 > len(seq):
                continue

            frag = seq[start0:end0]

            req = _is_valid_require_run_chars(require_run_chars)

            # ---- nested signature filter on FULL LENGTH ----
            if req:
                full_len_seq = frag.upper()
                if not has_nested_run_signature(full_len_seq, req, base_min=base_min, flank_min=flank_min):
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
      TE (e.g. CP002684.1:7065623-7075763) -> "Order/Superfamily/Clade" (e.g. "LTR/Gypsy/Retand")

    cls.tsv columns (header):
      TE  Order  Superfamily  Clade  Complete  Strand  Domains
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
                    continue  # skip header

            cols = line.split("\t")
            if len(cols) < 4:
                continue

            te, order, superfam, clade = cols[0], cols[1], cols[2], cols[3]
            if not te:
                continue

            # If any fields are missing, keep placeholders to maintain consistency
            order = order or "unknown"
            superfam = superfam or "unknown"
            clade = clade or "unknown"

            te2ann[te] = f"{order}/{superfam}/{clade}"

    return te2ann

# -----------------------------
# Step 8: build kmer2ltr.domain from stitched SCN
# -----------------------------

def scn_to_kmer2ltr_domain(stitched_scn: str, out_domain: str, tesorter_cls_tsv: Optional[str] = None):
    """
    Writes:
      {chrom}:{s(ret)}-{e(ret)}#{Order/Superfamily/Clade} <TAB> max(l(lLTR), l(rLTR))

    If tesorter_cls_tsv is None or TE not found in mapping, we still write a suffix
    for consistency:
      #{LTR/unknown/unknown}
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
    """
    Checks whether Kmer2LTR is runnable via:
      python3 tools/Kmer2LTR/Kmer2LTR.py -h
    We treat help text presence as usable even if returncode != 0.
    """
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
    """
    Ensures ./tools/Kmer2LTR exists (git clone) and is usable.
    Returns the Kmer2LTR script path: ./tools/Kmer2LTR/Kmer2LTR.py
    """
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
    """
    miniprot uses -h, not --help. Also, some versions may exit nonzero
    while still printing usage text. We treat that as usable.
    """
    if not bin_path.exists() or not os.access(str(bin_path), os.X_OK):
        return False
    try:
        r = run([str(bin_path), "-h"], capture=True)
        txt = (r.stdout or "") + "\n" + (r.stderr or "")
        return ("Usage: miniprot" in txt) or (r.returncode == 0)
    except Exception:
        return False

def tool_usable_minimap2(bin_path: Path) -> bool:
    # minimap2 supports --help and typically returns 0
    return tool_usable_generic(bin_path, ["--help"])

def ensure_tools(tools_dir: Path) -> Tuple[str, str]:
    tools_dir = mkdirp(tools_dir)
    mm2_dir = tools_dir / "minimap2"
    mp_dir = tools_dir / "miniprot"

    mm2_bin = mm2_dir / "minimap2"
    mp_bin = mp_dir / "miniprot"

    need_mm2 = not tool_usable_minimap2(mm2_bin)
    need_mp  = not tool_usable_miniprot(mp_bin)

    if need_mm2:
        if not mm2_dir.exists():
            run(["git", "clone", "https://github.com/lh3/minimap2", str(mm2_dir)], check=True)
        run(["make"], cwd=str(mm2_dir), check=True)
        if not tool_usable_minimap2(mm2_bin):
            raise RuntimeError(f"minimap2 build failed or binary unusable at: {mm2_bin}")

    if need_mp:
        if not mp_dir.exists():
            run(["git", "clone", "https://github.com/lh3/miniprot", str(mp_dir)], check=True)
        run(["make"], cwd=str(mp_dir), check=True)
        if not tool_usable_miniprot(mp_bin):
            raise RuntimeError(f"miniprot build failed or binary unusable at: {mp_bin}")

    return str(mm2_bin), str(mp_bin)

# -----------------------------
# Step 9: TEsorter (+ optional tree)
# -----------------------------

def tool_usable_tesorter(te_dir: Path) -> bool:
    """
    Checks whether TEsorter is runnable via:
      PYTHONPATH=./tools/TEsorter python3 -m TEsorter -h
    We treat "usage: TEsorter" in stdout/stderr as usable even if returncode != 0.
    """
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
    """
    Ensures ./tools/TEsorter exists (git clone) and is usable.
    Returns the TEsorter module dir to be used as PYTHONPATH.
    """
    tools_dir = mkdirp(tools_dir)
    te_dir = tools_dir / "TEsorter"

    if not tool_usable_tesorter(te_dir):
        if not te_dir.exists():
            run(["git", "clone", "https://github.com/cwb14/TEsorter.git", str(te_dir)], check=True)
        # Re-check after clone
        if not tool_usable_tesorter(te_dir):
            raise RuntimeError(
                f"TEsorter appears unusable in: {te_dir}\n"
                f"Try: PYTHONPATH={te_dir} python3 -m TEsorter -h"
            )

    return str(te_dir)

def tool_usable_trfmod(bin_path: Path) -> bool:
    """
    trf-mod prints usage with -h in some builds; others just error nonzero.
    We'll accept either a 0 exit, or any stderr/stdout containing 'trf' / 'TRF'.
    """
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
    """
    Ensures ./tools/TRF-mod exists and is compiled.
    Build command:
      make -f compile.mak
    Returns path to ./tools/TRF-mod/trf-mod
    """
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
                 db: str, cov: int, evalue: str, rule: str, threads: int) -> Tuple[str, str, str]:
    """
    Runs TEsorter in outdir so outputs land in {prefix}.work/,
    but passes stitched_fa as an absolute path so it can be found.
    """
    mkdirp(outdir)

    env = dict(os.environ)
    env["PYTHONPATH"] = tesorter_py_path

    stitched_fa_abs = str(Path(stitched_fa).resolve())  # <<< FIX

    cmd = [
        "python3", "-m", "TEsorter",
        stitched_fa_abs,                 # <<< FIX
        "-db", db,
        "-p", str(threads),
        "-cov", str(cov),
        "-eval", str(evalue),
#        "--no-cleanup",
        "-rule", str(rule),
    ]

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
                domain_file: Optional[str] = None) -> str:
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

    r = subprocess.run(cmd, cwd=str(outdir), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"Kmer2LTR failed:\n{(r.stderr or '').strip()}")

    main_out = outdir / out_prefix
    if not main_out.exists():
        raise RuntimeError(f"Kmer2LTR finished but did not create expected output: {main_out}")

    return str(main_out)
    
def _parse_interval_from_kmer2ltr_col1(col1: str) -> Optional[Tuple[str, int, int]]:
    # col1: Athal_chr1:4406006-4411120#LTR/Copia/Bianca
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


def _overlap_fraction_of_shorter(a: Tuple[int, int], b: Tuple[int, int]) -> float:
    a0, a1 = a
    b0, b1 = b
    inter = min(a1, b1) - max(a0, b0) + 1
    if inter <= 0:
        return 0.0
    lena = a1 - a0 + 1
    lenb = b1 - b0 + 1
    return inter / min(lena, lenb)

def dedup_kmer2ltr_tsv(kmer2ltr_tsv: str, out_tsv: str, threshold: float) -> None:
    """
    Dedup a Kmer2LTR output TSV (main output file), keeping the lowest p-distance (col7).
    Tie-breaker: if p-distance ties, keep the record with the largest aln_len (col3).
    Assumes:
      - p-distance is column 7 (1-based) => parts[6]
      - aln_len is column 3 (1-based)    => parts[2]
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
        # choose lowest p; if tie, choose highest aln (column3)
        best = min(cluster, key=lambda x: (float(x["p"]), -int(x["aln"])))  # type: ignore
        out_handle.write(best["line"])  # type: ignore

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
                p = float(parts[6])  # col7 (0-based 6)
            except ValueError:
                out.write(raw)
                continue

            # Tie-breaker column: aln_len (col3, 0-based 2)
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
    """
    Kmer2LTR first column is the LTR-RT name; your tesorter full-length FASTA uses
    >{TE}#LTR/{Superfamily}/{Clade} so it should match exactly.
    """
    keep = set()
    with open(dedup_tsv, "r") as f:
        for line in f:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if cols:
                keep.add(cols[0])
    return keep


def build_ltr_rt_tree_from_tesorter(stitched_fa: str, tesorter_py_path: str, outdir: Path,
                                   db: str, iqtree3_path: str, ltr_tree_r: str):
    """
    Optional, slow:
      PYTHONPATH=./tools/TEsorter python3 -m TEsorter.modules.concatenate_domains <cls.pep> GAG PROT RH RT INT > <full.aln>
      iqtree3 -s <full.aln> -bb 1000 -nt AUTO
      Rscript <LTR_tree.R> <treefile> <cls.tsv> <out.pdf>

    IMPORTANT:
    - We run with cwd=outdir, so pass filenames (basenames) NOT outdir-prefixed paths.
    """
    env = dict(os.environ)
    env["PYTHONPATH"] = tesorter_py_path

    stitched_fa_abs = str(Path(stitched_fa).resolve())
    base = Path(stitched_fa_abs).name  # e.g. ltrharvest7.ltrharvest.stitched.fa

    # These files are created in outdir by TEsorter, so refer to them by basename while cwd=outdir
    cls_pep_name = f"{base}.{db}.cls.pep"
    cls_tsv_name = f"{base}.{db}.cls.tsv"

    # Keep your original naming convention for the alignment outputs (also in outdir)
    aln_name = f"{base}.rexdb.cls.pep.full.aln"
    aln_path = outdir / aln_name

    # concatenate_domains writes to stdout -> capture and write ourselves
    cmd_cat = ["python3", "-m", "TEsorter.modules.concatenate_domains",
               cls_pep_name, "GAG", "PROT", "RH", "RT", "INT"]

    p = subprocess.run(cmd_cat, cwd=str(outdir), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env)
    if p.returncode != 0:
        raise RuntimeError(f"TEsorter concatenate_domains failed:\n{(p.stderr or '').strip()}")

    aln_path.write_text(p.stdout or "")

    # iqtree3 outputs multiple files next to aln_name (in outdir)
    run([iqtree3_path, "-s", aln_name, "-bb", "1000", "-nt", "AUTO"], cwd=str(outdir), check=True, capture=True)

    treefile_name = aln_name + ".treefile"
    out_pdf = str(outdir / f"{Path(stitched_fa_abs).stem}.TEsorter_tree.pdf")

    run(["Rscript", ltr_tree_r, treefile_name, cls_tsv_name, out_pdf], cwd=str(outdir), check=True, capture=True)
    
def build_tesorter_full_length_ltr_fasta_from_cls_tsv(cls_tsv_path: str, genome_fa: str, out_fa: str):
    """
    Build full-length LTR-RT FASTA using genome + TEsorter cls.tsv.

    cls.tsv columns (with header):
      TE  Order  Superfamily  Clade  Complete  Strand  Domains

    We keep only Order == 'LTR' (column2).
    Header format:
      >{TE}#LTR/{Superfamily}/{Clade}

    If Strand == '-', output reverse-complement of the extracted genomic sequence.
    TE format expected:
      chrom:start-end
    """
    genome = load_fasta_as_dict(genome_fa)

    def parse_te_interval(te: str) -> Optional[Tuple[str, int, int]]:
        # robust split on last ":" in case chrom names contain ":" (rare but safer)
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

            # Handle header row (starts with 'TE\tOrder...' in your example)
            if not header_seen:
                header_seen = True
                if line.startswith("TE\t") or line.lower().startswith("te\t"):
                    continue  # skip header
                # if no header, fall through and parse as data

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

    Windows are tiny (7bp by design), so brute force is fine.
    """
    left = left.upper()
    right = right.upper()

    if len(left) < min_len or len(right) < min_len:
        return False

    maxk = min(len(left), len(right))

    # check longer kmers first so we short-circuit early
    for k in range(maxk, min_len - 1, -1):
        for i in range(0, len(left) - k + 1):
            mer = left[i:i+k]

            # reject low complexity: must contain >=2 distinct bases
            if len(set(mer)) < 2:
                continue

            if mer in right:
                return True

    return False

def rescue_nonautonomous_by_tsd_from_scn(
    stitched_scn: str,
    genome_fa: str,
    tesorter_retained_te_keys: set,
    out_fa: str,
    min_len: int = 5,
    require_run_chars: Optional[List[str]] = None,
    base_min: int = 800,
    flank_min: int = 80,
) -> int:
    """
    For each SCN candidate NOT present in tesorter_retained_te_keys, scan for TSD.

    Window definition (matches your example intent):
      left  = genome[(sret-6) .. sret]  (6 upstream + first base)  => 7bp
      right = genome[eret .. (eret+6)]  (last base + 6 downstream) => 7bp
    We call TSD if left/right share an exact k-mer length >= min_len.

    Writes full-length rescued sequences with header:
      >chrom:sret-eret#LTR/unknown/unknown

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

            # full-length bounds (1-based inclusive -> 0-based half-open)
            fl0 = sret1 - 1
            fr0 = eret1
            if fl0 < 0 or fr0 > len(seq):
                continue

            # Build 7bp windows (clamped)
            left_start0 = max(0, fl0 - 6)
            left_end0   = min(len(seq), fl0 + 1)   # include first base
            right_start0 = max(0, fr0 - 1)         # last base
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
                    
            hdr = f"{te_key}#LTR/unknown/unknown"
            out.write(f">{hdr}\n")
            for i in range(0, len(frag), 60):
                out.write(frag[i:i+60] + "\n")
            rescued += 1

    if rescued == 0:
        Path(out_fa).touch()
    return rescued


# -----------------------------
# main
# -----------------------------

def main():
    ap = argparse.ArgumentParser(
        description="miniprot genic mask + minimap2 non-LTR mask + chunk + ltrharvest parallel + stitch (GFF3+SCN) + stitched FASTA",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,  # (3) show defaults in -h
    )

    ap.add_argument("--genome", required=True, help="Input genome FASTA")
    ap.add_argument("--proteins", default=None, help="Protein FASTA for miniprot (optional; if omitted, skip gene masking)")
    ap.add_argument("--te-library", default=None, help="TE library FASTA (optional; if omitted, skip non-LTR TE masking)")
    ap.add_argument("--out-prefix", required=True, help="Output prefix")
    ap.add_argument("--threads", type=int, default=8, help="Threads for miniprot/minimap2 and chunk parallelism")
    ap.add_argument("--workdir", default=None, help="Working directory (default: {out_prefix}.work)")
    ap.add_argument("--gt", default="gt", help="Path to GenomeTools 'gt' executable")

    # (1) cleaning flag: keep only stitched outputs (+ stitched fasta) and ./tools; remove {prefix}.work
    ap.add_argument("--clean", action="store_true", help="Remove {out_prefix}.work after success (keeps ./tools and stitched outputs)")

    # miniprot tuning
    ap.add_argument("--outn", type=int, default=1000, help="miniprot --outn")
    ap.add_argument("--outs", type=float, default=0.99, help="miniprot --outs")
    ap.add_argument("--outc", type=float, default=0.1, help="miniprot --outc")
    ap.add_argument("--gene-mask", default="CDS", choices=["mRNA", "CDS", "both"], help="Which miniprot GFF features to mask in the genome") # CDS is importablt for simulated data where TEs may insert into genes. Probably less importnat in real data.

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

    # SCN size-based filtering (optional; applied to merged_scn before FASTA/TEsorter/Kmer2LTR)
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
    # default ON
    ap.add_argument(
        "--kmer2ltr-domains",
        dest="kmer2ltr_domains",
        action="store_true",
        default=True,
        help="Pass {out_prefix}.kmer2ltr.domain to Kmer2LTR via -D"
    )
    # allow turning it OFF explicitly
    ap.add_argument(
        "--no-kmer2ltr-domains",
        dest="kmer2ltr_domains",
        action="store_false",
        help="Do not pass a domain file to Kmer2LTR"
    )


    # TEsorter optional toggle (default ON to preserve current behavior)
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

    ap.add_argument("--tesorter-db", default="rexdb-plant",
                    help="TEsorter HMM database (-db)")
    ap.add_argument("--tesorter-cov", type=int, default=20,
                    help="TEsorter min coverage (-cov)")
    ap.add_argument("--tesorter-eval", default="1e-2",
                    help="TEsorter max E-value (-eval)")
    ap.add_argument("--tesorter-rule", default="80-80-80", # Probably dont use 70-30-80. It gives lots. mostly junk.
                    help="TEsorter pass2 rule (-rule)")

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
        "--tsd-rescue",
        action="store_true",
        help="Rescue non-autonomous LTR-RTs that fail TEsorter by scanning for a TSD near intact boundaries (default: off)."
    )
    ap.add_argument(
        "--tsd-min-len",
        type=int,
        default=5,
        help="Minimum exact shared k-mer length to call a TSD (default: 5)."
    )

    # Optional phylogeny (slow!)
    ap.add_argument("--tesorter-tree", action="store_true",
                    help="Build LTR-RT phylogeny from TEsorter domains (VERY slow; adds iqtree3 + Rscript runtime)")
    ap.add_argument("--iqtree3", default="iqtree3", help="Path to iqtree3 (only used if --tesorter-tree)")
    ap.add_argument("--ltr-tree-r", default="./tools/TEsorter/scripts/LTR_tree.R",
                    help="Path to TEsorter LTR_tree.R (only used if --tesorter-tree)")

    ap.add_argument("--ltr-tools", default="both", choices=["both", "ltrharvest", "ltr_finder"],
                    help="Which LTR caller(s) to run per chunk")

    ap.add_argument("--ltrfinder", default="ltr_finder",
                    help="Path to ltr_finder v1.07 executable")

    ap.add_argument(
        "--require-run-chars",
        default=None,
        help="Comma-separated chars that must each occur as a run of 100 in the FULL LTR-RT before internal is kept (e.g. 'Y' or 'Y,R'). Optional."
    )
    ap.add_argument("--nested-base-min", type=int, default=800,
                    help="Min length for the first char in --require-run-chars (default: 800)")
    ap.add_argument("--nested-flank-min", type=int, default=80,
                    help="Min length for each flanking nesting char (default: 80)")

    ap.add_argument(
        "--ltr-timeout",
        type=int,
        default=0,
        help="Max seconds allowed per chunk per LTR tool call (0 = no timeout). Timed-out chunks salvage partial outputs."
    )

    ap.add_argument(
        "--ltrfinder-args",
#        default="-w 2 -C -D 15000 -d 100 -L 7000 -l 100 -p 15 -M 0.75 -S 5.0",
        default="-w 2 -C -D 15000 -d 100 -L 7000 -l 100 -p 20 -M 0.00 -S 0.0", 
        help="Quoted string of args appended to ltr_finder (NOTE: -w 2 is required)"
    )

    # ltrharvest args string
    ap.add_argument(
        "--ltrharvest-args",
#        default="-mindistltr 100 -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -similar 70 -vic 30 -seed 15 -seqids yes -xdrop 10",
        default="-mindistltr 100 -minlenltr 100 -maxlenltr 7000 -mintsd 0 -maxtsd 0 -similar 70 -vic 60 -seed 15 -seqids yes -xdrop 10",
        help="Quoted string of args appended to 'gt ltrharvest -index <idx>'"
    )

    args = ap.parse_args()
    out_prefix = args.out_prefix

    # Workdir is {prefix}.work
    workdir = Path(args.workdir or f"{out_prefix}.work")
    mkdirp(workdir)

    # (4) tools in ./tools (current dir), not inside workdir
    tools_dir = Path("./tools")
    minimap2_path, miniprot_path = ensure_tools(tools_dir)
    tesorter_py_path = ensure_tesorter(tools_dir)
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

    # enforce -w 2 requirement
    lf_args = args.ltrfinder_args.strip().split()
    if args.ltr_tools in ("both", "ltr_finder"):
        if "-w" not in lf_args:
            raise RuntimeError("ltr_finder requires '-w 2' output format; please include '-w 2' in --ltrfinder-args")
        # If user provided '-w' but not '2', catch that too
        for i, tok in enumerate(lf_args):
            if tok == "-w":
                if i + 1 >= len(lf_args) or lf_args[i + 1] != "2":
                    raise RuntimeError("ltr_finder output format must be '-w 2' for SCN conversion compatibility")


    # Step 2    
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

    # Step 4
    print("[Step4] chunking masked genome...")
    chunks_dir = str(workdir / "chunks")
    chunks = make_chunks(gene_te_masked_fa, chunks_dir, args.size, args.overlap)
    if not chunks:
        raise RuntimeError("No chunks produced (empty genome?)")

    # Step 5
    print(f"[Step5] LTR calling on {len(chunks)} chunks with {args.threads} workers...")
    ltr_args = args.ltrharvest_args.strip().split()
    lf_args = args.ltrfinder_args.strip().split()

    ltr_work = str(workdir / "ltr_runs")
    mkdirp(Path(ltr_work))

    harvest_pairs: List[Tuple[ChunkInfo, str]] = []
    harvest_gffs: List[Tuple[ChunkInfo, str, str]] = []  # (chunk, scn, gff) for stitch_gff3
    finder_pairs: List[Tuple[ChunkInfo, str]] = []

    failures = 0
    total = len(chunks)
    completed = 0

    start = time.monotonic()
    last_update = 0.0
    update_interval = 1.0  # seconds

    print("", file=sys.stderr)  # ensure stderr stream exists

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
                elapsed = now - start
                pct = (completed / total) * 100 if total else 100.0
                rate = (completed / elapsed) if elapsed > 0 else 0.0
                remaining = ((total - completed) / rate) if rate > 0 else float("inf")
                eta = _format_hms(remaining)

                msg = (
                    f"\r[Step5] LTR calling: {completed}/{total} "
                    f"({pct:.2f}%) | ETA {eta}"
                )
                print(msg, end="", file=sys.stderr, flush=True)
                last_update = now

    print("", file=sys.stderr)  # newline after progress bar

    if failures:
        print(f"[WARN] {failures} chunks failed; stitching continues using successes only.")

    # Step 6 stitched outputs (NOW INSIDE workdir)
    stitched_gff3 = str(workdir / f"{out_prefix}.ltrharvest.stitched.gff3")
    stitched_harvest_scn = str(workdir / f"{out_prefix}.ltrharvest.stitched.scn")
    stitched_finder_scn  = str(workdir / f"{out_prefix}.ltrfinder.stitched.scn")

    # Keep merged SCN as your "primary" output (unchanged location unless you also want it in workdir)
    merged_scn = f"{out_prefix}.ltrtools.stitched.scn"

    if args.ltr_tools in ("both", "ltrharvest"):
        print(f"[Step6] stitching GFF3 (ltrharvest) -> {stitched_gff3}")
        stitch_gff3(harvest_gffs, stitched_gff3)

        print(f"[Step6] stitching SCN  (ltrharvest) -> {stitched_harvest_scn}")
        stitch_scn_from_pairs(harvest_pairs, stitched_harvest_scn)
    else:
        # If not running ltrharvest, make sure expected outputs exist (optional)
        Path(stitched_gff3).write_text("##gff-version 3\n")
        Path(stitched_harvest_scn).touch()

    if args.ltr_tools in ("both", "ltr_finder"):
        print(f"[Step6] stitching SCN  (ltr_finder) -> {stitched_finder_scn}")
        stitch_scn_from_pairs(finder_pairs, stitched_finder_scn)
    else:
        Path(stitched_finder_scn).touch()

    # Merge whichever stitched SCNs exist into primary merged_scn
    to_merge = []
    if args.ltr_tools in ("both", "ltrharvest"):
        to_merge.append(stitched_harvest_scn)
    if args.ltr_tools in ("both", "ltr_finder"):
        to_merge.append(stitched_finder_scn)

    print(f"[Step6] merging stitched SCN(s) -> {merged_scn}")
    merge_stitched_scns(to_merge, merged_scn)

    # Optional: size-based filtering of merged SCN to reduce downstream runtime
    if any([
        args.scn_min_ltr_len,
        args.scn_min_ret_len,
        args.scn_max_ret_len,
        args.scn_min_int_len,
        args.scn_max_int_len,
    ]):
        filtered_scn = f"{out_prefix}.ltrtools.stitched.filtered.scn"
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

        merged_scn = filtered_scn  # downstream steps use the filtered SCN
    else:
        print("[Step6b] skipping SCN length filtering (no --scn-* filters set).")

    # (2) build LTR FASTA from stitched SCN using s(ret)/e(ret)
    internals_fa = str(workdir / f"{out_prefix}.ltrtools.internals.fa")
    intact_fa = f"{out_prefix}.ltrtools.intact.fa"  # requested path (in ./)

    if args.use_tesorter:
        print(f"[Step7] building INTERNALS FASTA from SCN -> {internals_fa}")
        req_chars = None
        if args.require_run_chars:
            req_chars = [x.strip() for x in args.require_run_chars.split(",") if x.strip()]

        scn_to_internal_fasta(merged_scn, args.genome, internals_fa,
                             require_run_chars=req_chars,
                             base_min=args.nested_base_min,
                             flank_min=args.nested_flank_min)

    else:
        print(f"[Step7] building INTACT FASTA from SCN -> {intact_fa}")

        req_chars = None
        if args.require_run_chars:
            req_chars = [x.strip() for x in args.require_run_chars.split(",") if x.strip()]

        scn_to_intact_fasta(
            merged_scn,
            args.genome,
            intact_fa,
            require_run_chars=req_chars,
        )

    # Step 9: TEsorter on stitched FASTA (required)
    cls_tsv_path = None
    tesorter_lib_fa = None

    if args.use_tesorter:
        tesorter_outdir = workdir  # store ALL TEsorter outputs in {prefix}.work/
        print(f"[Step9] running TEsorter on stitched FASTA (outputs -> {tesorter_outdir})")

        cls_lib_path, cls_pep_path, cls_tsv_path = run_tesorter(
            stitched_fa=internals_fa,
            tesorter_py_path=tesorter_py_path,
            outdir=tesorter_outdir,
            db=args.tesorter_db,
            cov=args.tesorter_cov,
            evalue=args.tesorter_eval,
            rule=args.tesorter_rule,
            threads=args.threads,
        )

        # Build FULL-LENGTH LTR-RT FASTA from genome + cls.tsv (Order == LTR)
        tesorter_lib_fa = str(workdir / f"{out_prefix}.ltrharvest.full_length.fa.{args.tesorter_db}.cls.lib.fa")
        print(f"[Step9] building full-length LTR FASTA from cls.tsv -> {tesorter_lib_fa}")
        build_tesorter_full_length_ltr_fasta_from_cls_tsv(cls_tsv_path, args.genome, tesorter_lib_fa)

        # Optional: rescue non-autonomous LTR-RTs by TSD scan (only for those NOT retained by TEsorter)
        rescued_fa = None
        k2l_merged_fa = tesorter_lib_fa

        if args.tsd_rescue:
            # Build retained set from cls.tsv (Order == LTR)
            retained = set()
            with open(cls_tsv_path, "r") as fin:
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
                    if len(cols) < 2:
                        continue
                    te, order = cols[0], cols[1]
                    if order == "LTR" and te:
                        retained.add(te)

            rescued_fa = str(workdir / f"{out_prefix}.tsd_rescued.full_length.fa")
            print(f"[Step9] TSD rescue enabled: scanning non-TEsorter LTR-RTs -> {rescued_fa}")     
            n_rescued = rescue_nonautonomous_by_tsd_from_scn(
                stitched_scn=merged_scn,
                genome_fa=args.genome,
                tesorter_retained_te_keys=retained,
                out_fa=rescued_fa,
                min_len=args.tsd_min_len,
                require_run_chars=req_chars,
                base_min=args.nested_base_min,
                flank_min=args.nested_flank_min,
            )

            print(f"[Step9] TSD rescue: rescued {n_rescued} candidates")

            # Merge TEsorter full-length + rescued into a single Kmer2LTR input
            k2l_merged_fa = str(workdir / f"{out_prefix}.full_length.protein_plus_tsd.fa")
            print(f"[Step9] merging protein-backed + TSD-rescued full-length FASTAs -> {k2l_merged_fa}")

            seen_hdrs = set()
            with open(k2l_merged_fa, "w") as out:
                # Write TEsorter fasta first (headers already include #LTR/SF/Clade)
                for name, seq in iter_fasta(tesorter_lib_fa):
                    if name in seen_hdrs:
                        continue
                    seen_hdrs.add(name)
                    out.write(f">{name}\n")
                    for i in range(0, len(seq), 60):
                        out.write(seq[i:i+60] + "\n")

                # Then rescued (headers are #LTR/unknown/unknown, matching domain file format)
                if rescued_fa and Path(rescued_fa).exists() and Path(rescued_fa).stat().st_size > 0:
                    for name, seq in iter_fasta(rescued_fa):
                        if name in seen_hdrs:
                            continue
                        seen_hdrs.add(name)
                        out.write(f">{name}\n")
                        for i in range(0, len(seq), 60):
                            out.write(seq[i:i+60] + "\n")

            # IMPORTANT: downstream Kmer2LTR should use the merged file
            tesorter_lib_fa = k2l_merged_fa
    else:
        print("[Step9] skipping TEsorter (--no-tesorter).")

    # Step 8: kmer2ltr.domain from stitched SCN  (MOVE THIS UP before Kmer2LTR)
    kmer2ltr_domain = str(workdir / f"{out_prefix}.kmer2ltr.domain")
    print(f"[Step8] building kmer2ltr domain -> {kmer2ltr_domain}")
    scn_to_kmer2ltr_domain(merged_scn, kmer2ltr_domain, tesorter_cls_tsv=cls_tsv_path)


    # Step 9b: Kmer2LTR on full-length FASTA (for dedup)
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
    )

    # Primary output #2 (in ./): dedupbed full-length FASTA
    k2l_dedup_out = f"{out_prefix}_kmer2ltr_dedup"
    print(f"[Step9b] deduping Kmer2LTR output -> {k2l_dedup_out}")
    dedup_kmer2ltr_tsv(k2l_main, k2l_dedup_out, threshold=args.dedup_threshold)

    keep_names = names_from_kmer2ltr_dedup(k2l_dedup_out)

    if args.use_tesorter:
        tesorter_lib_fa_dedup = f"{out_prefix}.ltrharvest.full_length.dedup.fa.{args.tesorter_db}.cls.lib.fa"
        print(f"[Step9b] subsetting full-length FASTA by dedup list -> {tesorter_lib_fa_dedup}")
        subset_fasta_by_name_set(k2l_in_fa, tesorter_lib_fa_dedup, keep_names)

    else:
        intact_dedup_fa = f"{out_prefix}.ltrtools.intact.dedup.fa"
        print(f"[Step9b] subsetting intact FASTA by dedup list -> {intact_dedup_fa}")
        subset_fasta_by_name_set(intact_fa, intact_dedup_fa, keep_names)


    # Optional: tree building (VERY slow)
    if args.tesorter_tree:
        print("[Step10] building LTR-RT phylogeny from TEsorter (VERY slow)...")
        which_or_die("Rscript")
        which_or_die(args.iqtree3)
        if not Path(args.ltr_tree_r).exists():
            raise RuntimeError(f"LTR_tree.R not found at: {args.ltr_tree_r}")

        build_ltr_rt_tree_from_tesorter(
            stitched_fa=internals_fa,
            tesorter_py_path=tesorter_py_path,
            outdir=tesorter_outdir,
            db=args.tesorter_db,
            iqtree3_path=args.iqtree3,
            ltr_tree_r=args.ltr_tree_r,
        )


    print("\nDone.")
    print("Primary outputs:")
    print(f"  Kmer2LTR dedup TSV:          {out_prefix}_kmer2ltr_dedup")

    if args.use_tesorter:
        print(f"  TEsorter full-length dedup FASTA: {out_prefix}.ltrharvest.full_length.dedup.fa.{args.tesorter_db}.cls.lib.fa")
        print("")
        print("Workdir outputs:")
        print(f"  Full-length (pre-dedup) FASTA: {tesorter_lib_fa}")
    else:
        print(f"  Intact dedup FASTA:          {out_prefix}.ltrtools.intact.dedup.fa")
        print(f"  Intact (pre-dedup) FASTA:    {out_prefix}.ltrtools.intact.fa")

    print(f"  Kmer2LTR outputs prefix:       {workdir}/{out_prefix}_kmer2ltr*")

    if args.tesorter_tree:
        # NOTE: your current code also references stitched_fa here; fix that too (below)
        print(f"  TEsorter tree PDF (in workdir): {Path(internals_fa).stem}.TEsorter_tree.pdf")

    # (1) cleaning: remove {prefix}.work entirely
    if args.clean:
        print(f"[CLEAN] removing workdir: {workdir}")
        shutil.rmtree(workdir, ignore_errors=True)

if __name__ == "__main__":
    main()
