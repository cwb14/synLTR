#!/usr/bin/env python3
"""
LTRharvest pipeline (SCN fixed + minor edits)
LTRharvest parameters are aggressive. Masking genes and non-LTR-RT TEs helps, but still, many false positives with real data. TEsorter is used to clean them up.  

Key fix:
- DO NOT use `-out` for SCN. `-out` produces FASTA sequences.
- For SCN/tabular output (LTR_retriever-style), capture stdout with `-tabout yes`.
- GFF3 is written to a file via `-gff3 <filename>`.

Outputs (outside of {prefix}.work):
- {prefix}.ltrharvest.stitched.gff3  (lifted to genome coords; seqid=chrom)
- {prefix}.ltrharvest.stitched.scn   (lifted to genome coords; appended chrom as last field)
- {prefix}.ltrharvest.stitched.fa    (FASTA extracted from stitched SCN using s(ret)/e(ret))

Tools:
- Uses ./tools/minimap2/minimap2 and ./tools/miniprot/miniprot by default.
- If missing/unusable, auto-downloads and compiles into ./tools/.

Default parameters are benchmarked with PrinTE:
- bash ./PrinTE/PrinTE.sh --cds_percent 17 -itp 5 -sz 168.5Mb -mb ../Athal.ancestral.LTR.bins.freq2 -TsTv 2.0 --ex_LTR -ftp 60 -t 100 -m 2e-9 --burnin_only --TE_lib ../maize_rice_arab_curated_TE.lib.fa -cn 5

- python ltrharvest.py --genome burnin.fasta --proteins PrinTE/data/TAIR10.pep.fa.gz --te-library lib_clean.fa --out-prefix myrun --threads 100

- python ./PrinTE/util/bedtools.py -pass_scn myrun.ltrharvest.stitched.scn -bed <(cat burnin.bed | grep -v FRAG | grep -v gene | grep LTR) -r 0.90 -print unique-bed
  Overlapping entries: 1800 (1793 unique)
  Entries unique to SCN/PASS file: 8
  Entries unique to BED file: 6



"""

import argparse
import os
import re
import shutil
import subprocess
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


# -----------------------------
# Step 2: miniprot gene masking
# -----------------------------

@dataclass
class Interval:
    start: int  # 0-based inclusive
    end: int    # 0-based exclusive

def parse_mrna_intervals_from_gff_text(gff_text: str) -> Dict[str, List[Interval]]:
    intervals: Dict[str, List[Interval]] = {}
    for line in gff_text.splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 9:
            continue
        seqid, ftype = parts[0], parts[2]
        if ftype != "mRNA":
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
                           outn: int = 1, outs: float = 0.8, outc: float = 0.5):
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

    intervals = merge_intervals_per_chrom(parse_mrna_intervals_from_gff_text(gff_text))
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

def run_ltrharvest_scn_and_gff3(gt_path: str, indexname: str, scn_path: str, gff3_path: str, ltrharvest_args: List[str]):
    """
    - SCN/tabular output is stdout (with -tabout yes)
    - GFF3 output goes to gff3_path (with -gff3 <file>)
    """
    cmd = [gt_path, "ltrharvest", "-index", indexname] + ltrharvest_args + ["-tabout", "yes", "-gff3", gff3_path]
    r = run(cmd, capture=True)
    if r.returncode != 0:
        raise RuntimeError(f"gt ltrharvest failed for index {indexname}:\n{(r.stderr or '').strip()}")

    Path(scn_path).write_text(r.stdout or "")

def process_one_chunk(gt_path: str, chunk: ChunkInfo, work_dir: str, ltrharvest_args: List[str]) -> Tuple[ChunkInfo, str, str]:
    wdir = mkdirp(Path(work_dir) / chunk.chunk_id)
    local_fa = str(wdir / Path(chunk.fasta_path).name)
    shutil.copyfile(chunk.fasta_path, local_fa)

    indexname = str(wdir / "idx")
    build_suffixerator_index(gt_path, local_fa, indexname)

    scn_path = str(wdir / f"{chunk.chunk_id}.ltrharvest.scn")
    gff3_path = str(wdir / f"{chunk.chunk_id}.ltrharvest.gff3")

    run_ltrharvest_scn_and_gff3(gt_path, indexname, scn_path, gff3_path, ltrharvest_args)
    return chunk, scn_path, gff3_path


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

# -----------------------------
# Step 7: build LTR FASTA from stitched SCN
# -----------------------------

def scn_to_ltr_fasta(stitched_scn: str, genome_fa: str, out_fa: str):
    """
    Uses s(ret) and e(ret) (1-based inclusive) plus chrom (last column) to extract sequences from genome_fa.

    Output record format:
      >chr1:10006055-10007712
      ACTG...
    """
    genome = load_fasta_as_dict(genome_fa)

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
            s1 = int(sret)
            e1 = int(eret)
            if e1 < s1:
                s1, e1 = e1, s1

            seq = genome.get(chrom)
            if seq is None:
                continue

            start0 = s1 - 1
            end0 = e1  # inclusive -> exclusive
            if start0 < 0 or end0 > len(seq):
                continue

            header = f"{chrom}:{s1}-{e1}"
            out.write(f">{header}\n")
            frag = seq[start0:end0]
            for i in range(0, len(frag), 60):
                out.write(frag[i:i+60] + "\n")
            n_written += 1

    if n_written == 0:
        Path(out_fa).touch()


# -----------------------------
# Step 8: build kmer2ltr.domain from stitched SCN
# -----------------------------

def scn_to_kmer2ltr_domain(stitched_scn: str, out_domain: str):
    """
    Writes: {chrom}:{s(ret)}-{e(ret)} <TAB> max(l(lLTR), l(rLTR))

    Assumes stitched SCN has chromosome appended as last field, like:
      c[0]=sret, c[1]=eret, c[5]=l(lLTR), c[8]=l(rLTR), c[11]=chrom
    """
    n_written = 0
    with open(stitched_scn, "r") as fin, open(out_domain, "w") as out:
        for line in fin:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            c = line.split()
            # Need at least 12 columns to match your one-liner indexing
            if len(c) < 12:
                continue

            sret, eret = c[0], c[1]
            chrom = c[11]

            if not (sret.isdigit() and eret.isdigit()):
                continue
            try:
                ll = int(c[5])
                rl = int(c[8])
            except ValueError:
                continue

            out.write(f"{chrom}:{sret}-{eret}\t{max(ll, rl)}\n")
            n_written += 1

    if n_written == 0:
        Path(out_domain).touch()

# -----------------------------
# Tool setup (./tools/)
# -----------------------------

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

def run_tesorter(stitched_fa: str, tesorter_py_path: str, outdir: Path,
                 db: str, cov: int, evalue: str, rule: str, threads: int) -> Tuple[str, str]:
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
        "-rule", str(rule),
    ]

    r = subprocess.run(cmd, cwd=str(outdir), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env)
    if r.returncode != 0:
        raise RuntimeError(f"TEsorter failed:\n{(r.stderr or '').strip()}")

    base = Path(stitched_fa_abs).name  # <<< FIX: derive base from abs path
    cls_lib = str(outdir / f"{base}.{db}.cls.lib")
    cls_pep = str(outdir / f"{base}.{db}.cls.pep")
    return cls_lib, cls_pep

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
    
def subset_tesorter_cls_lib_to_fasta(cls_lib_path: str, out_fa: str):
    """
    Keep only sequences with 'LTR' in header line, and trim header after the first space.
    Example:
      >Athal_chr5:...#LTR/Copia/Ale Athal_chr5:...
    becomes:
      >Athal_chr5:...#LTR/Copia/Ale
    """
    n_written = 0
    with open(cls_lib_path, "r") as fin, open(out_fa, "w") as out:
        write_it = False
        for line in fin:
            if line.startswith(">"):
                hdr = line[1:].rstrip("\n")
                if "LTR" in hdr:
                    trimmed = hdr.split(" ", 1)[0]
                    out.write(">" + trimmed + "\n")
                    write_it = True
                    n_written += 1
                else:
                    write_it = False
            else:
                if write_it:
                    out.write(line)
    if n_written == 0:
        Path(out_fa).touch()



# -----------------------------
# main
# -----------------------------

def main():
    ap = argparse.ArgumentParser(
        description="miniprot genic mask + minimap2 non-LTR mask + chunk + ltrharvest parallel + stitch (GFF3+SCN) + stitched FASTA",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,  # (3) show defaults in -h
    )

    ap.add_argument("--genome", required=True, help="Input genome FASTA")
    ap.add_argument("--proteins", required=True, help="Protein FASTA for miniprot")
    ap.add_argument("--te-library", required=True, help="TE library FASTA (records containing 'LTR' in header are excluded for step3)")
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

    # TEsorter (required)
    ap.add_argument("--tesorter-db", default="rexdb-plant",
                    help="TEsorter HMM database (-db)")
    ap.add_argument("--tesorter-cov", type=int, default=10,
                    help="TEsorter min coverage (-cov)")
    ap.add_argument("--tesorter-eval", default="1e-2",
                    help="TEsorter max E-value (-eval)")
    ap.add_argument("--tesorter-rule", default="70-30-80",
                    help="TEsorter pass2 rule (-rule)")

    # Optional phylogeny (slow!)
    ap.add_argument("--tesorter-tree", action="store_true",
                    help="Build LTR-RT phylogeny from TEsorter domains (VERY slow; adds iqtree3 + Rscript runtime)")
    ap.add_argument("--iqtree3", default="iqtree3", help="Path to iqtree3 (only used if --tesorter-tree)")
    ap.add_argument("--ltr-tree-r", default="./tools/TEsorter/scripts/LTR_tree.R",
                    help="Path to TEsorter LTR_tree.R (only used if --tesorter-tree)")

    # ltrharvest args string
    ap.add_argument(
        "--ltrharvest-args",
        default="-mindistltr 100 -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -similar 70 -vic 30 -seed 15 -seqids yes -xdrop 10",
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

    which_or_die("bedtools")
    which_or_die(args.gt)

    # Step 2
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
    )

    # Step 3
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

    # (5) store prefix.gene_TEmasked.fa inside prefix.work
    gene_te_masked_fa = str(workdir / f"{out_prefix}.gene_TEmasked.fa")
    run(["bedtools", "maskfasta", "-fi", genic_masked_fa, "-bed", merged_bed, "-fo", gene_te_masked_fa], check=True, capture=True)

    # Step 4
    print("[Step4] chunking masked genome...")
    chunks_dir = str(workdir / "chunks")
    chunks = make_chunks(gene_te_masked_fa, chunks_dir, args.size, args.overlap)
    if not chunks:
        raise RuntimeError("No chunks produced (empty genome?)")

    # Step 5
    print(f"[Step5] ltrharvest on {len(chunks)} chunks with {args.threads} workers...")
    ltr_args = args.ltrharvest_args.strip().split()
    ltr_work = str(workdir / "ltrharvest_runs")
    mkdirp(Path(ltr_work))

    chunk_triplets: List[Tuple[ChunkInfo, str, str]] = []
    failures = 0
    with ThreadPoolExecutor(max_workers=max(1, args.threads)) as ex:
        futs = {ex.submit(process_one_chunk, args.gt, c, ltr_work, ltr_args): c for c in chunks}
        for fut in as_completed(futs):
            c = futs[fut]
            try:
                chunk, scn_path, gff3_path = fut.result()
                chunk_triplets.append((chunk, scn_path, gff3_path))
            except Exception as e:
                failures += 1
                print(f"[WARN] chunk failed: {c.chunk_id} :: {e}")

    if failures:
        print(f"[WARN] {failures} chunks failed; stitching continues using successes only.")

    # Step 6 stitched outputs (outside workdir)
    stitched_gff3 = f"{out_prefix}.ltrharvest.stitched.gff3"
    stitched_scn  = f"{out_prefix}.ltrharvest.stitched.scn"

    print(f"[Step6] stitching GFF3 -> {stitched_gff3}")
    stitch_gff3(chunk_triplets, stitched_gff3)

    print(f"[Step6] stitching SCN  -> {stitched_scn}")
    stitch_scn(chunk_triplets, stitched_scn)

    # (2) build LTR FASTA from stitched SCN using s(ret)/e(ret)
    stitched_fa = f"{out_prefix}.ltrharvest.stitched.fa"
    print(f"[Step7] building FASTA from SCN -> {stitched_fa}")
    scn_to_ltr_fasta(stitched_scn, args.genome, stitched_fa)

    # Step 9: TEsorter on stitched FASTA (required)
    tesorter_outdir = workdir  # store ALL TEsorter outputs in {prefix}.work/
    print(f"[Step9] running TEsorter on stitched FASTA (outputs -> {tesorter_outdir})")
    cls_lib_path, cls_pep_path = run_tesorter(
        stitched_fa=stitched_fa,
        tesorter_py_path=tesorter_py_path,
        outdir=tesorter_outdir,
        db=args.tesorter_db,
        cov=args.tesorter_cov,
        evalue=args.tesorter_eval,
        rule=args.tesorter_rule,
        threads=args.threads,
    )

    # Subset cls.lib into a simple FASTA in ./ (outside workdir)
    tesorter_lib_fa = f"{out_prefix}.ltrharvest.stitched.fa.{args.tesorter_db}.cls.lib.fa"
    print(f"[Step9] subsetting cls.lib -> {tesorter_lib_fa}")
    subset_tesorter_cls_lib_to_fasta(cls_lib_path, tesorter_lib_fa)

    # Optional: tree building (VERY slow)
    if args.tesorter_tree:
        print("[Step10] building LTR-RT phylogeny from TEsorter (VERY slow)...")
        which_or_die("Rscript")
        which_or_die(args.iqtree3)
        if not Path(args.ltr_tree_r).exists():
            raise RuntimeError(f"LTR_tree.R not found at: {args.ltr_tree_r}")

        build_ltr_rt_tree_from_tesorter(
            stitched_fa=stitched_fa,
            tesorter_py_path=tesorter_py_path,
            outdir=tesorter_outdir,
            db=args.tesorter_db,
            iqtree3_path=args.iqtree3,
            ltr_tree_r=args.ltr_tree_r,
        )


    # Step 8: kmer2ltr.domain from stitched SCN
    kmer2ltr_domain = f"{out_prefix}.kmer2ltr.domain"
    print(f"[Step8] building kmer2ltr domain -> {kmer2ltr_domain}")
    scn_to_kmer2ltr_domain(stitched_scn, kmer2ltr_domain)

    print("\nDone.")
    print("Key outputs:")
    print(f"  Stitched GFF3:  {stitched_gff3}")
    print(f"  Stitched SCN:   {stitched_scn}")
    print(f"  Stitched FASTA: {stitched_fa}")
    print(f"  kmer2ltr domain: {kmer2ltr_domain}")
    print(f"  Tools dir:      ./tools/")

    print(f"  TEsorter subset lib FASTA: {tesorter_lib_fa}")
    if args.tesorter_tree:
        print(f"  TEsorter tree PDF (in workdir): {Path(stitched_fa).stem}.TEsorter_tree.pdf")

    # (1) cleaning: remove {prefix}.work entirely
    if args.clean:
        print(f"[CLEAN] removing workdir: {workdir}")
        shutil.rmtree(workdir, ignore_errors=True)


if __name__ == "__main__":
    main()
