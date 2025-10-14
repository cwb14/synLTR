#!/usr/bin/env python3
"""
segment_selfaln.py  (enhanced + TRF-mod + sdust)

Split a genome FASTA into sliding windows and run minimap2 self-alignments
for each segment. Output:
  - <out_prefix>.paf          (raw, '+' strand only, genomic coords)
  - <out_prefix>.filtered.paf (deduped + length filters + optional TSD filter, with 6 extra columns)
  - <out_prefix>.filtered.fa  (multi-seq FASTA of LTR-RTs)
  - <out_prefix>.genome_noLTRs.fa (genome FASTA with those spans excised)

Masking per segment (in this order, merged):
  - longdust (mandatory)
  - TRF-mod (default ON; disable with --disable-trf)
  - sdust   (default ON; disable with --disable-sdust; threshold with --sd-t)

Defaults:
  window=30000; step=5000; threads=3; parallel=1
  minimap2 flags: -X -c --cs=short -k13 -w5 -m25 -r2k
  Filtering defaults:
    --min-internal 100  --max-internal 25000
    --min-ltr 100       --max-ltr 7000
    --tsd-min-len 5     --tsd-max-len 15
    --tsd-max-subs 1    --tsd-search-space 0
    --max-internal-lc 1.0   (skip LC filter unless set < 1.0)
  sdust default threshold: --sd-t 16

Notes:
  - TRF-mod is invoked with: -b5 -g5 -s30 -p200
  - Use --disable-TSD-search to skip TSD detection and filtering (TSD-related args ignored).
  - Use --max-internal-lc to cap the fraction of **uppercase** 'N' (low complexity) in the internal region.
    Lowercase 'n' (protein masking) is ignored in this calculation.
  - Use --keep-temp to retain the <out_prefix>_out working directory.
  - SCN generation and all related parameters were removed in the previous version.
"""
import argparse
import os
import sys
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List, Tuple, Dict, Iterable
import re

# --------------------------- FASTA helpers ---------------------------

def fasta_iter(fa_path):
    """Yield (name, sequence_string) for each record in a FASTA."""
    name = None
    seq_chunks = []
    with open(fa_path, 'r') as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(seq_chunks)
                name = line[1:].strip().split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if name is not None:
            yield name, ''.join(seq_chunks)

def read_fasta_to_dict(fa_path: str) -> Dict[str, str]:
    """Load whole FASTA into memory (contig -> sequence)."""
    return {name: seq for name, seq in fasta_iter(fa_path)}

def write_fasta_record(out_handle, name: str, seq: str, width: int = 60):
    out_handle.write(f">{name}\n")
    for i in range(0, len(seq), width):
        out_handle.write(seq[i:i+width] + "\n")

def segments_for_length(seq_len, window, step):
    """Yield (start, end) 1-based inclusive for windows fully within the sequence."""
    if seq_len < window:
        return
    last_start = seq_len - window + 1
    s = 1
    while s <= last_start:
        e = s + window - 1
        yield s, e
        s += step

# --------------------------- command runner ---------------------------

def run_cmd(cmd, cwd=None, env=None):
    p = subprocess.run(cmd, cwd=cwd, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"Command failed: {' '.join(cmd)}\nSTDOUT:\n{p.stdout}\nSTDERR:\n{p.stderr}")
    return p.stdout

# --------------------------- minimap2 presence / install ---------------------------

def ensure_minimap2_in_cwd() -> tuple[str, bool, Optional[Path]]:
    """
    Ensure there's an executable ./minimap2.
    If not present, clone+build into a hidden temp dir, copy binary to ./minimap2.
    Returns (path_to_mm2, installed_by_us, build_dir).
    """
    here = Path.cwd()
    mm2_path = here / "minimap2"
    if mm2_path.is_file() and os.access(mm2_path, os.X_OK):
        return str(mm2_path), False, None

    build_dir = here / ".mm2_build"
    repo_dir = build_dir / "minimap2"
    build_dir.mkdir(exist_ok=True)
    if not repo_dir.exists():
        run_cmd(["git", "clone", "https://github.com/lh3/minimap2", str(repo_dir)])
    else:
        try:
            run_cmd(["git", "pull"], cwd=str(repo_dir))
        except Exception:
            pass
    run_cmd(["make"], cwd=str(repo_dir))
    built_bin = repo_dir / "minimap2"
    if not (built_bin.is_file() and os.access(built_bin, os.X_OK)):
        raise RuntimeError("Failed to build minimap2 (binary not found).")
    shutil.copy2(built_bin, mm2_path)
    mm2_path.chmod(mm2_path.stat().st_mode | 0o111)
    return str(mm2_path), True, build_dir

# --------------------------- miniprot presence / install ---------------------------

def ensure_miniprot_in_cwd() -> tuple[str, bool, Optional[Path]]:
    here = Path.cwd()
    candidates = [here / "miniprot", here / "miniprot" / "miniprot"]
    for c in candidates:
        if c.is_file() and os.access(c, os.X_OK):
            return str(c), False, None

    build_dir = here / ".mip_build"
    repo_dir = build_dir / "miniprot"
    build_dir.mkdir(exist_ok=True)
    if not repo_dir.exists():
        run_cmd(["git", "clone", "https://github.com/lh3/miniprot", str(repo_dir)])
    else:
        try:
            run_cmd(["git", "pull"], cwd=str(repo_dir))
        except Exception:
            pass
    run_cmd(["make"], cwd=str(repo_dir))
    built_bin = repo_dir / "miniprot"
    if not (built_bin.is_file() and os.access(built_bin, os.X_OK)):
        raise RuntimeError("Failed to build miniprot (binary not found).")
    dest_bin = here / "miniprot"
    shutil.copy2(built_bin, dest_bin)
    dest_bin.chmod(dest_bin.stat().st_mode | 0o111)
    return str(dest_bin), True, build_dir

# --------------------------- longdust presence / install ---------------------------

def ensure_longdust_in_cwd() -> tuple[str, bool, Optional[Path]]:
    here = Path.cwd()
    candidates = [here / "longdust", here / "longdust" / "longdust"]
    for c in candidates:
        if c.is_file() and os.access(c, os.X_OK):
            return str(c), False, None

    build_dir = here / ".ld_build"
    repo_dir = build_dir / "longdust"
    build_dir.mkdir(exist_ok=True)
    if not repo_dir.exists():
        run_cmd(["git", "clone", "https://github.com/lh3/longdust.git", str(repo_dir)])
    else:
        try:
            run_cmd(["git", "pull"], cwd=str(repo_dir))
        except Exception:
            pass
    run_cmd(["make"], cwd=str(repo_dir))
    built_bin = repo_dir / "longdust"
    if not (built_bin.is_file() and os.access(built_bin, os.X_OK)):
        raise RuntimeError("Failed to build longdust (binary not found).")
    dest_bin = here / "longdust"
    shutil.copy2(built_bin, dest_bin)
    dest_bin.chmod(dest_bin.stat().st_mode | 0o111)
    return str(dest_bin), True, build_dir

# --------------------------- TRF-mod presence / install ---------------------------

def ensure_trf_mod_in_cwd() -> tuple[Optional[str], bool, Optional[Path]]:
    """
    Try to ensure ./trf-mod exists. Returns (path or None if disabled/unavailable, installed_by_us, build_dir).
    """
    here = Path.cwd()
    trf_path = here / "trf-mod"
    if trf_path.is_file() and os.access(trf_path, os.X_OK):
        return str(trf_path), False, None

    build_dir = here / ".trf_build"
    repo_dir = build_dir / "TRF-mod"
    try:
        build_dir.mkdir(exist_ok=True)
        if not repo_dir.exists():
            run_cmd(["git", "clone", "https://github.com/lh3/TRF-mod", str(repo_dir)])
        else:
            try:
                run_cmd(["git", "pull"], cwd=str(repo_dir))
            except Exception:
                pass
        # TRF-mod uses a special makefile
        run_cmd(["make", "-f", "compile.mak"], cwd=str(repo_dir))
        # The executable name in repo root is expected to be 'trf-mod'
        built_bin = repo_dir / "trf-mod"
        if not (built_bin.is_file() and os.access(built_bin, os.X_OK)):
            raise RuntimeError("TRF-mod built but binary not found.")
        shutil.copy2(built_bin, trf_path)
        trf_path.chmod(trf_path.stat().st_mode | 0o111)
        return str(trf_path), True, build_dir
    except Exception as e:
        print(f"WARNING: TRF-mod unavailable ({e}). Will continue without it.", file=sys.stderr)
        # Cleanup partial build dir to avoid clutter
        try: shutil.rmtree(build_dir, ignore_errors=True)
        except Exception: pass
        return None, False, None

# --------------------------- barrnap presence / install ---------------------------

def ensure_barrnap_in_cwd() -> tuple[Optional[str], bool, Optional[Path]]:
    """
    Ensure a clone of tseemann/barrnap exists and return path to its bin/barrnap.
    We keep it in a hidden build dir and call the script there (it's a perl wrapper).
    Returns (path_or_None, installed_by_us, build_dir)
    """
    here = Path.cwd()
    build_dir = here / ".bnp_build"
    repo_dir = build_dir / "barrnap"
    exe = repo_dir / "bin" / "barrnap"
    try:
        build_dir.mkdir(exist_ok=True)
        if not repo_dir.exists():
            run_cmd(["git", "clone", "https://github.com/tseemann/barrnap.git", str(repo_dir)])
        else:
            try:
                run_cmd(["git", "pull"], cwd=str(repo_dir))
            except Exception:
                pass
        if not exe.exists():
            print("[WARN] Barrnap executable not found after clone; expected bin/barrnap", file=sys.stderr)
            return None, False, None
        # ensure executable bit (perl script, but let's be nice)
        try:
            os.chmod(exe, (exe.stat().st_mode | 0o111))
        except Exception:
            pass
        return str(exe), True, build_dir
    except Exception as e:
        print(f"WARNING: Barrnap unavailable ({e}). Will continue without it.", file=sys.stderr)
        try: shutil.rmtree(build_dir, ignore_errors=True)
        except Exception: pass
        return None, False, None

def parse_barrnap_gff_to_intervals_by_contig(gff_text: str) -> Dict[str, List[Tuple[int,int]]]:
    """
    Parse Barrnap GFF3 into 1-based inclusive intervals per contig.
    """
    out: Dict[str, List[Tuple[int,int]]] = {}
    for line in gff_text.splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 5:
            continue
        chrom = parts[0]
        try:
            s1 = int(parts[3])  # GFF is 1-based inclusive
            e1 = int(parts[4])
        except Exception:
            continue
        if e1 >= s1:
            out.setdefault(chrom, []).append((s1, e1))
    # normalize/merge per contig
    for c in list(out.keys()):
        out[c] = merge_intervals(out[c])
    return out

def run_barrnap_and_mask(fasta: str, out_prefix: str, threads: int,
                         kingdom: str, reject: float, barrnap_path: str) -> Tuple[str, str, str]:
    """
    Run Barrnap over the (whole) genome FASTA and hard-mask rDNA intervals.
    Returns (gff_path, bed_path, masked_fasta)
    """
    gff_path = f"{out_prefix}.barrnap.gff3"
    bed_path = f"{out_prefix}.barrnap.mask.bed"
    masked_fasta = f"{out_prefix}.barrnap.hardmasked.fa"

    # barrnap prints GFF to stdout
    cmd = [
        barrnap_path,
        "--kingdom", str(kingdom),
        "--reject", str(reject),
        "--threads", str(threads),
        fasta
    ]
    gff_text = run_cmd(cmd)  # raises on nonzero
    with open(gff_path, "w") as gfh:
        gfh.write(gff_text)

    ivals_by_ctg = parse_barrnap_gff_to_intervals_by_contig(gff_text)

    # also emit a 0-based BED for record-keeping
    with open(bed_path, "w") as bo:
        for ctg, lst in ivals_by_ctg.items():
            for s1, e1 in lst:
                bo.write(f"{ctg}\t{s1-1}\t{e1}\n")

    # hard-mask original FASTA (uppercase 'N')
    hardmask_fasta_by_intervals(fasta, ivals_by_ctg, masked_fasta, bed_path)
    return gff_path, bed_path, masked_fasta


# --------------------------- sdust presence / install ---------------------------

def ensure_sdust_in_cwd() -> tuple[Optional[str], bool, Optional[Path]]:
    here = Path.cwd()
    # Accept either a top-level file "./sdust" or "./sdust/sdust" inside a directory
    candidates = [here / "sdust", here / "sdust" / "sdust"]
    for c in candidates:
        if c.is_file() and os.access(c, os.X_OK):
            return str(c), False, None

    build_dir = here / ".sd_build"
    repo_dir = build_dir / "sdust"
    try:
        build_dir.mkdir(exist_ok=True)
        if not repo_dir.exists():
            run_cmd(["git", "clone", "https://github.com/lh3/sdust.git", str(repo_dir)])
        else:
            try:
                run_cmd(["git", "pull"], cwd=str(repo_dir))
            except Exception:
                pass
        run_cmd(["make"], cwd=str(repo_dir))
        built_bin = repo_dir / "sdust"
        if not (built_bin.is_file() and os.access(built_bin, os.X_OK)):
            raise RuntimeError("sdust built but binary not found.")

        # Decide install destination:
        # - if "./sdust" is a directory, install as "./sdust/sdust"
        # - otherwise install as "./sdust" (file)
        dest_dir = here / "sdust"
        if dest_dir.exists() and dest_dir.is_dir():
            dest_file = dest_dir / "sdust"
            dest_dir.mkdir(exist_ok=True)
        else:
            dest_file = here / "sdust"  # top-level file

        shutil.copy2(built_bin, dest_file)
        dest_file.chmod(dest_file.stat().st_mode | 0o111)
        return str(dest_file), True, build_dir
    except Exception as e:
        print(f"WARNING: sdust unavailable ({e}). Will continue without it.", file=sys.stderr)
        try: shutil.rmtree(build_dir, ignore_errors=True)
        except Exception: pass
        return None, False, None



# --------------------------- miniprot masking helpers ---------------------------

def parse_mrna_intervals_from_gff(gff_text: str) -> Dict[str, List[Tuple[int,int]]]:
    intervals: Dict[str, List[Tuple[int,int]]] = {}
    for line in gff_text.splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 5:
            continue
        typ = parts[2]
        if typ != "mRNA":
            continue
        chrom = parts[0]
        try:
            start = int(parts[3])
            end = int(parts[4])
        except ValueError:
            continue
        if start > end:
            start, end = end, start
        intervals.setdefault(chrom, []).append((start, end))
    return intervals

def hardmask_fasta_by_intervals(in_fasta: str, intervals: Dict[str, List[Tuple[int,int]]],
                                out_fasta: str, bed_out: str):
    with open(out_fasta, "w") as fo, open(bed_out, "w") as bo:
        for contig, seq in fasta_iter(in_fasta):
            if contig not in intervals:
                write_fasta_record(fo, contig, seq)
                continue
            seq_list = list(seq)
            L = len(seq_list)
            for s, e in intervals[contig]:
                s0 = max(1, min(s, L))
                e0 = max(1, min(e, L))
                if s0 > e0:
                    s0, e0 = e0, s0
                for i in range(s0 - 1, e0):
                    seq_list[i] = 'N'
                bo.write(f"{contig}\t{s0-1}\t{e0}\n")
            write_fasta_record(fo, contig, ''.join(seq_list))

def run_miniprot_and_mask(fasta: str, protein: str, out_prefix: str,
                          threads: int, outn: int, outs: float, outc: float,
                          miniprot_path: str) -> Tuple[str, str, str]:
    gff_path = f"{out_prefix}.gff"
    bed_path = f"{out_prefix}.mask.bed"
    masked_fasta = f"{out_prefix}.hardmasked.fa"

    cmd = [
        miniprot_path,
        "--gff-only",
        "-t", str(threads),
        fasta,
        protein,
        "-P", out_prefix,
        "--outn", str(outn),
        "--outs", str(outs),
        "--outc", str(outc),
    ]
    gff_text = run_cmd(cmd)
    with open(gff_path, "w") as gfh:
        gfh.write(gff_text)

    intervals = parse_mrna_intervals_from_gff(gff_text)
    if not intervals:
        shutil.copyfile(fasta, masked_fasta)
        with open(bed_path, "w"):
            pass
        return gff_path, bed_path, masked_fasta

    hardmask_fasta_by_intervals(fasta, intervals, masked_fasta, bed_path)
    return gff_path, bed_path, masked_fasta

# --------------------------- BED helpers (0-based -> 1-based inclusive) ---------------------------

def parse_bed_zero_based_to_inclusive1(bed_text: str) -> Dict[str, List[Tuple[int,int]]]:
    """
    Parse generic BED-like (chrom, start0, end0, ...) into 1-based inclusive intervals.
    Ignores malformed lines.
    """
    out: Dict[str, List[Tuple[int, int]]] = {}
    for line in bed_text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.strip().split()
        if len(parts) < 3:
            continue
        name, s, e = parts[0], parts[1], parts[2]
        try:
            s0 = int(s)
            e0 = int(e)
        except ValueError:
            continue
        if e0 <= s0:
            continue
        out.setdefault(name, []).append((s0 + 1, e0))  # convert to 1-based inclusive
    return out

def merge_intervals(intervals: List[Tuple[int,int]]) -> List[Tuple[int,int]]:
    """
    Merge 1-based inclusive intervals; treat touching intervals (end >= next.start - 1) as mergeable.
    """
    if not intervals:
        return []
    intervals = sorted(intervals, key=lambda x: (x[0], x[1]))
    merged = [intervals[0]]
    for s, e in intervals[1:]:
        ms, me = merged[-1]
        if s <= me + 1:
            merged[-1] = (ms, max(me, e))
        else:
            merged.append((s, e))
    return merged

# --------------------------- per-segment masking tool runners ---------------------------

def run_longdust_on_segment(seg_fasta: str, longdust_path: str,
                            k: int, w: int, t: float, e: int,
                            forward_only: bool, approx: bool,
                            bed_dir: Path) -> Dict[str, List[Tuple[int,int]]]:
    bed_dir.mkdir(parents=True, exist_ok=True)
    base = Path(seg_fasta).stem
    out_bed = bed_dir / f"{base}.ld.bed"

    cmd = [longdust_path]
    if k is not None: cmd += ["-k", str(k)]
    if w is not None: cmd += ["-w", str(w)]
    if t is not None: cmd += ["-t", str(t)]
    if e is not None: cmd += ["-e", str(e)]
    if forward_only:  cmd += ["-f"]
    if approx:        cmd += ["-a"]
    cmd += [seg_fasta]

    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"longdust failed for {os.path.basename(seg_fasta)}\nSTDERR:\n{p.stderr.strip()}")

    with open(out_bed, "w") as bo:
        bo.write(p.stdout)

    return parse_bed_zero_based_to_inclusive1(p.stdout)

def run_trfmod_on_segment(seg_fasta: str, trf_path: Optional[str], bed_dir: Path) -> Dict[str, List[Tuple[int,int]]]:
    """
    Run TRF-mod twice if available:
      - microsatellites:          -b5 -g5 -s30  -p200
      - large tandem repeats:     -a2 -b7 -g7 -A80 -G10 -s1000 -p2000
    Both write BED-like to stdout; we parse and merge to 1-based inclusive intervals.
    """
    if not trf_path:
        return {}

    bed_dir.mkdir(parents=True, exist_ok=True)
    base = Path(seg_fasta).stem

    # --- pass 1: microsats
    p1 = subprocess.run([trf_path, "-b4", "-g5", "-s25", "-p500", seg_fasta], # More aggressive.
#   p1 = subprocess.run([trf_path, "-b5", "-g5", "-s30", "-p200", seg_fasta], # Less aggressive.
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p1.returncode != 0:
        print(f"WARNING: TRF-mod pass1 failed on {os.path.basename(seg_fasta)}; skipping TRF mask.\n{p1.stderr.strip()}",
              file=sys.stderr)
        return {}

    out_bed1 = bed_dir / f"{base}.trf.micro.bed"
    with open(out_bed1, "w") as bo:
        bo.write(p1.stdout)
    iv1 = parse_bed_zero_based_to_inclusive1(p1.stdout)

    # --- pass 2: large tandem repeats
    p2 = subprocess.run([trf_path, "-a2", "-b7", "-g7", "-A80", "-G10", "-s1000", "-p2000", seg_fasta],
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p2.returncode != 0:
        print(f"WARNING: TRF-mod pass2 failed on {os.path.basename(seg_fasta)}; continuing with pass1 only.\n{p2.stderr.strip()}",
              file=sys.stderr)
        iv2 = {}
    else:
        out_bed2 = bed_dir / f"{base}.trf.large.bed"
        with open(out_bed2, "w") as bo:
            bo.write(p2.stdout)
        iv2 = parse_bed_zero_based_to_inclusive1(p2.stdout)

    # merge per-seqname lists
    merged: Dict[str, List[Tuple[int,int]]] = {}
    for src in (iv1, iv2):
        for k, lst in src.items():
            merged.setdefault(k, []).extend(lst)
    for k in list(merged.keys()):
        merged[k] = merge_intervals(merged[k])

    return merged

def run_sdust_on_segment(seg_fasta: str, sdust_path: Optional[str], bed_dir: Path, sd_t: int) -> Dict[str, List[Tuple[int,int]]]:
    if not sdust_path:
        return {}
    bed_dir.mkdir(parents=True, exist_ok=True)
    base = Path(seg_fasta).stem
    out_bed = bed_dir / f"{base}.sd.bed"

    try:
        # Added threshold parameter: -t<sd_t>
        p = subprocess.run([sdust_path, f"-t{sd_t}", seg_fasta], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except OSError as e:
        # e.g., when sdust_path points to a directory or isn't executable
        print(f"WARNING: sdust exec failed on {os.path.basename(seg_fasta)} ({e}); skipping sdust mask.",
              file=sys.stderr)
        return {}

    if p.returncode != 0:
        print(f"WARNING: sdust failed on {os.path.basename(seg_fasta)}; skipping sdust mask.\n{p.stderr.strip()}",
              file=sys.stderr)
        return {}

    with open(out_bed, "w") as bo:
        bo.write(p.stdout)
    return parse_bed_zero_based_to_inclusive1(p.stdout)
    
def hardmask_segment_with_intervals(seg_fasta: str,
                                    intervals_by_seqname: Dict[str, List[Tuple[int,int]]],
                                    masked_fa: str, bed_path: str):
    """
    seqname for a segment FASTA is the header before first space (we write like 'contig:start-end').
    We'll map intervals using that seqname.
    """
    # reuse genome-wide function; it keys by contig name
    hardmask_fasta_by_intervals(seg_fasta, intervals_by_seqname, masked_fa, bed_path)

# --------------------------- Work unit ---------------------------

@dataclass
class SegmentJob:
    contig: str
    contig_len: int
    start: int   # 1-based inclusive
    end: int     # 1-based inclusive
    fasta_path: str

def write_segment_fasta(contig, seq, start, end, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    fname = f"{contig}_{start}_{end}.fa"
    fpath = os.path.join(out_dir, fname)
    header = f">{contig}:{start}-{end}\n"
    seg_seq = seq[start-1:end]
    with open(fpath, "w") as out:
        out.write(header)
        out.write(seg_seq + "\n")
    return fpath

def adjust_paf_coords(paf_text, contig, contig_len, offset, keep_strand="+"):
    if not paf_text:
        return ""
    out_lines = []
    for line in paf_text.splitlines():
        if not line.strip():
            continue
        cols = line.split('\t')
        if len(cols) < 12:
            continue
        strand = cols[4]
        if keep_strand in {"+", "-"} and strand != keep_strand:
            continue
        try:
            _qlen = int(cols[1]); qs = int(cols[2]); qe = int(cols[3])
            _tlen = int(cols[6]); ts = int(cols[7]); te = int(cols[8])
        except ValueError:
            continue
        qs_adj = qs + offset
        qe_adj = qe + offset
        ts_adj = ts + offset
        te_adj = te + offset
        cols[0] = contig
        cols[1] = str(contig_len)
        cols[2] = str(qs_adj)
        cols[3] = str(qe_adj)
        cols[5] = contig
        cols[6] = str(contig_len)
        cols[7] = str(ts_adj)
        cols[8] = str(te_adj)
        out_lines.append('\t'.join(cols))
    return '\n'.join(out_lines) + ('\n' if out_lines else '')

def run_minimap2_segment(job: SegmentJob, minimap2, threads, k, w, m, r):
    cmd = [
        minimap2,
        "-X",
        "-c",
        "--cs=short",
        f"-k{k}",
        f"-w{w}",
        f"-m{m}",
        f"-r{r}",
        job.fasta_path,
        job.fasta_path,
        "-t", str(threads),
    ]
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"minimap2 failed for {os.path.basename(job.fasta_path)}\nSTDERR:\n{p.stderr.strip()}")
    offset = job.start - 1
    adjusted = adjust_paf_coords(p.stdout, job.contig, job.contig_len, offset, keep_strand="+")
    return adjusted

def process_segment(job: SegmentJob,
                    minimap2: str, threads: int, k: int, w: int, m: int, r: str,
                    longdust_path: str,
                    ld_k: int, ld_w: int, ld_t: float, ld_e: int,
                    ld_forward: bool, ld_approx: bool,
                    use_trf: bool, trf_path: Optional[str],
                    use_sdust: bool, sdust_path: Optional[str], sd_t: int,
                    seg_mask_dir: Path, mask_bed_dir: Path) -> str:
    """
    For each segment: run longdust (mandatory), plus TRF-mod/sdust if enabled & available.
    Merge their BEDs, hardmask the segment once, then run minimap2.
    """
    base = Path(job.fasta_path).stem
    seg_mask_dir.mkdir(parents=True, exist_ok=True)
    mask_bed_dir.mkdir(parents=True, exist_ok=True)

    # Collect intervals per tool (keyed by the segment sequence name, e.g. "ctg:start-end")
    intervals_all: Dict[str, List[Tuple[int,int]]] = {}

    # longdust
    ld_iv = run_longdust_on_segment(job.fasta_path, longdust_path, ld_k, ld_w, ld_t, ld_e,
                                    ld_forward, ld_approx, bed_dir=mask_bed_dir)
    for kctg, lst in ld_iv.items():
        intervals_all.setdefault(kctg, []).extend(lst)

    # TRF-mod (optional)
    if use_trf and trf_path:
        trf_iv = run_trfmod_on_segment(job.fasta_path, trf_path, bed_dir=mask_bed_dir)
        for kctg, lst in trf_iv.items():
            intervals_all.setdefault(kctg, []).extend(lst)

    # sdust (optional)
    if use_sdust and sdust_path:
        sd_iv = run_sdust_on_segment(job.fasta_path, sdust_path, bed_dir=mask_bed_dir, sd_t=sd_t)
        for kctg, lst in sd_iv.items():
            intervals_all.setdefault(kctg, []).extend(lst)

    # Merge intervals and mask once
    merged_by_seqname: Dict[str, List[Tuple[int,int]]] = {
        k: merge_intervals(v) for k, v in intervals_all.items()
    }

    masked_path = seg_mask_dir / f"{base}.fa"
    bed_out = mask_bed_dir / f"{base}.merged.bed"
    # Write merged BED (0-based half-open) for reference while masking
    with open(bed_out, "w") as bo:
        for seqname, lst in merged_by_seqname.items():
            for s, e in lst:
                bo.write(f"{seqname}\t{s-1}\t{e}\n")

    hardmask_segment_with_intervals(job.fasta_path, merged_by_seqname, str(masked_path), str(bed_out))

    # run minimap2 on masked segment
    job_for_mm2 = SegmentJob(job.contig, job.contig_len, job.start, job.end, str(masked_path))
    return run_minimap2_segment(job_for_mm2, minimap2, threads, k, w, m, r)

# --------------------------- Filtering helpers ---------------------------

def parse_paf_core(line: str):
    """
    Parse essential fields from a PAF line.
    Returns tuple (contig, qs, qe, ts, te) for '+'-strand self alignments.
    Assumes columns >= 9 exist.
    """
    parts = line.rstrip("\n").split("\t")
    if len(parts) < 9:
        return None
    contig_q = parts[0]
    contig_t = parts[5]
    if contig_q != contig_t:
        return None
    try:
        qs = int(parts[2]); qe = int(parts[3])
        ts = int(parts[7]); te = int(parts[8])
    except ValueError:
        return None

    # Canonicalize so "left" is earlier on the genome
    left_s, left_e, right_s, right_e = (qs, qe, ts, te)
    if qs > ts:
        left_s, left_e, right_s, right_e = (ts, te, qs, qe)

    return contig_q, left_s, left_e, right_s, right_e

def ltr_lengths(left_s: int, left_e: int, right_s: int, right_e: int) -> Tuple[int,int,int,int]:
    left_len  = max(0, left_e - left_s)
    right_len = max(0, right_e - right_s)
    internal_len = max(0, right_s - left_e)
    ltr_rt_len = max(0, right_e - left_s)
    return ltr_rt_len, left_len, internal_len, right_len

def hamming(a: str, b: str) -> int:
    return sum(1 for x, y in zip(a, b) if x != y)
    
# ### NEW: masked segment lookup for low-complexity 'N' counting
def _load_masked_segments(seg_mask_dir: str) -> Dict[str, List[Tuple[int,int,str]]]:
    """
    Load masked segment FASTAs from seg_mask_dir into:
      { contig: [(seg_start, seg_end, seq_str), ...] }   (sorted by seg_start)
    Filenames are expected to be '<contig>_<start>_<end>.fa' with headers 'contig:start-end'.
    """
    out: Dict[str, List[Tuple[int,int,str]]] = {}
    p = Path(seg_mask_dir)
    if not p.exists():
        return out
    rx = re.compile(r"^(?P<c>.+)_(?P<s>\d+)_(?P<e>\d+)\.fa$")
    for f in p.glob("*.fa"):
        m = rx.match(f.name)
        if not m:
            continue
        contig = m.group("c")
        s = int(m.group("s")); e = int(m.group("e"))
        # read single-record FASTA quickly
        name = None
        seq_chunks = []
        with open(f, "r") as fh:
            for line in fh:
                if line.startswith(">"):
                    continue
                seq_chunks.append(line.strip())
        seq = "".join(seq_chunks)
        out.setdefault(contig, []).append((s, e, seq))
    for c in out:
        out[c].sort(key=lambda t: t[0])
    return out

def _find_covering_segment(segments: List[Tuple[int,int,str]], start: int, end: int) -> Optional[Tuple[int,int,str]]:
    """
    Return a (s,e,seq) segment such that s <= start and e >= end; else None.
    segments must be sorted by s.
    """
    # small linear scan is fine; segments are overlapping windows
    for s, e, seq in segments:
        if s <= start and e >= end:
            return (s, e, seq)
    return None

def find_tsd(seq: str, left_start: int, right_end: int,
             min_len: int, max_len: int, max_subs: int, search_space: int) -> Optional[Tuple[str, str]]:
    """
    Try to find a TSD flanking the LTR-RT:
      - Left TSD is immediately upstream of left_start.
      - Right TSD is immediately downstream of right_end.
      - Allow +/- search_space bp jitter around both boundaries.
      - Accept the LONGEST L in [min_len, max_len] with <= max_subs mismatches.
    Coordinates are 1-based inclusive; 'seq' is Python 0-based.
    Returns (left_tsd, right_tsd) or None if not found.
    """
    n = len(seq)
    best = None  # (L, subs, |dl|+|dr|, left_seq, right_seq)

    for L in range(max_len, min_len - 1, -1):
        for dl in range(-search_space, search_space + 1):
            ls = left_start - L + dl
            le = left_start - 1 + dl
            if ls < 1 or le > n:
                continue
            left = seq[ls - 1:le]
            for dr in range(-search_space, search_space + 1):
                rs = right_end + 1 + dr
                re = right_end + L + dr
                if rs < 1 or re > n:
                    continue
                right = seq[rs - 1:re]
                subs = hamming(left, right)
                if subs <= max_subs:
                    score = (L, -subs, -(abs(dl)+abs(dr)))
                    if (best is None) or (score > (best[0], -best[1], -best[2])):
                        best = (L, subs, abs(dl)+abs(dr), left, right)
        if best is not None:
            break

    if best is None:
        return None
    return best[3], best[4]

def filter_paf(
    paf_path: str,
    genome_fa: str,
    out_path: str,
    min_internal: int,
    max_internal: int,
    min_ltr: int,
    max_ltr: int,
    tsd_min_len: int,
    tsd_max_len: int,
    tsd_max_subs: int,
    tsd_search_space: int,
    disable_tsd_search: bool = False,
    *,
    # ### NEW: low-complexity filter
    seg_mask_dir: Optional[str] = None,
    max_internal_lc: float = 1.0,
) -> Tuple[int,int,int]:
    """
    Read raw PAF (genomic coords), remove duplicates, apply length (& optional TSD) filters,
    and write filtered lines with six extra columns:
      [total_LTR_RT_len, left_LTR_len, internal_len, right_LTR_len, left_TSD, right_TSD]
    If disable_tsd_search is True, TSD discovery is skipped (and not used for filtering);
    TSD columns are written as '.'.
    Returns (raw_lines, dedup_kept, final_kept).
    """

    # Load original (unmasked) genome to get flanks for TSD
    seqs = read_fasta_to_dict(genome_fa)
    # Preload masked segments only if LC filter is active (< 1.0) and dir provided
    segs_by_ctg: Dict[str, List[Tuple[int,int,str]]] = {}
    if seg_mask_dir and max_internal_lc < 1.0 - 1e-12:
        segs_by_ctg = _load_masked_segments(seg_mask_dir)

    seen_keys = set()    # canonical tuples to drop duplicates
    raw_total = 0
    dedup_kept = 0
    final_kept = 0

    with open(paf_path, "r") as fin, open(out_path, "w") as fout:
        for line in fin:
            if not line.strip():
                continue
            raw_total += 1
            core = parse_paf_core(line)
            if core is None:
                continue
            contig, ls, le, rs, re = core

            # ignore pathological ordering (overlap/inversion)
            if not (ls < le and rs < re and ls < rs):
                continue

            key = (contig, ls, le, rs, re)
            if key in seen_keys:
                continue
            seen_keys.add(key)
            dedup_kept += 1

            lr_len, left_len, internal_len, right_len = ltr_lengths(ls, le, rs, re)

            # length filters
            if not (min_ltr <= left_len  <= max_ltr):  continue
            if not (min_ltr <= right_len <= max_ltr):  continue
            if not (min_internal <= internal_len <= max_internal):  continue

            # ### NEW: internal low-complexity ('N') fraction gate (count uppercase N only)
            if segs_by_ctg and internal_len > 0:
                ctg_segs = segs_by_ctg.get(contig)
                if ctg_segs:
                    # internal region inclusive coordinates
                    is1 = le + 1
                    ie1 = rs - 1
                    seg = _find_covering_segment(ctg_segs, start=is1, end=ie1)
                    if seg is not None:
                        s0, e0, seg_seq = seg
                        a = max(is1, s0) - s0  # 0-based offset into seg_seq
                        b = min(ie1, e0) - s0  # inclusive
                        if b >= a:
                            internal_masked = seg_seq[a:b+1]
                            # count only uppercase 'N' (TRF-mod/sdust/longdust)
                            n_upper = sum(1 for ch in internal_masked if ch == 'N')
                            frac_lc = n_upper / internal_len
                            if frac_lc > max_internal_lc:
                                continue

            left_tsd = "."
            right_tsd = "."
            if not disable_tsd_search:
                seq = seqs.get(contig)
                if seq is None:
                    continue
                tsd = find_tsd(
                    seq, left_start=ls, right_end=re,
                    min_len=tsd_min_len, max_len=tsd_max_len,
                    max_subs=tsd_max_subs, search_space=tsd_search_space
                )
                if tsd is None:
                    continue
                left_tsd, right_tsd = tsd

            # Append six new fields to the original line (maintain original columns first)
            line = line.rstrip("\n")
            extra = f"\t{lr_len}\t{left_len}\t{internal_len}\t{right_len}\t{left_tsd}\t{right_tsd}\n"
            fout.write(line + extra)
            final_kept += 1

    return raw_total, dedup_kept, final_kept

# ---- LTR-RT FASTA + genome-without-LTRs ----

def emit_ltr_fasta_and_collect_intervals(
    filtered_paf_path: str,
    genome_fa: str,
    out_ltr_fa: str
) -> Dict[str, List[Tuple[int,int]]]:
    """
    From filtered PAF, write <out_prefix>.filtered.fa and return intervals per contig to excise.
    Each interval is (left_start, right_end), 1-based inclusive.
    """
    seqs = read_fasta_to_dict(genome_fa)
    intervals: Dict[str, List[Tuple[int,int]]] = {}
    count = 0
    with open(filtered_paf_path, "r") as fin, open(out_ltr_fa, "w") as fof:
        for line in fin:
            if not line.strip():
                continue
            core = parse_paf_core(line)
            if core is None:
                continue
            contig, ls, le, rs, re = core
            if not (ls < le and rs < re and ls < rs):
                continue
            # seq spans [ls, re]
            seq = seqs.get(contig)
            if not seq:
                continue
            subseq = seq[ls-1:re]  # Python slice, inclusive coordinates
            header = f"{contig}_{ls}_{le}_{rs}_{re}"
            write_fasta_record(fof, header, subseq)
            intervals.setdefault(contig, []).append((ls, re))
            count += 1
    print(f"Wrote LTR-RT FASTA with {count} sequences: {out_ltr_fa}", file=sys.stderr)
    return intervals

def write_genome_without_intervals(
    genome_fa: str,
    intervals_by_contig: Dict[str, List[Tuple[int,int]]],
    out_path: str
):
    """
    Excise (remove) provided 1-based inclusive intervals from each contig and write new FASTA.
    Overlapping/touching intervals are merged to avoid double-cutting.
    Contigs that become empty are skipped.
    """
    with open(out_path, "w") as fo:
        for contig, seq in fasta_iter(genome_fa):
            L = len(seq)
            ivals = intervals_by_contig.get(contig, [])
            if not ivals:
                write_fasta_record(fo, contig, seq)
                continue
            merged = merge_intervals([(max(1,s), min(L,e)) for s,e in ivals if 1 <= s <= e])
            if not merged:
                write_fasta_record(fo, contig, seq)
                continue
            # Build complement pieces
            pieces = []
            prev_end = 0
            for s,e in merged:
                # keep (prev_end+1 .. s-1)
                if s-1 > prev_end:
                    pieces.append(seq[prev_end:s-1])
                prev_end = e
            # tail after last cut
            if prev_end < L:
                pieces.append(seq[prev_end:L])
            kept = ''.join(pieces)
            if kept:
                write_fasta_record(fo, contig, kept)
            else:
                # Skip empty contig after excision
                print(f"Note: contig {contig} became empty after excision, skipping.", file=sys.stderr)
    print(f"Wrote genome with LTR-RTs excised: {out_path}", file=sys.stderr)

# --------------------------- Progress logger ---------------------------

def human_eta(seconds):
    if seconds is None or seconds < 0 or not (seconds < 9e18):
        return "unknown"
    m, s = divmod(int(seconds + 0.5), 60)
    h, m = divmod(m, 60)
    if h:
        return f"{h}h {m}m {s}s"
    if m:
        return f"{m}m {s}s"
    return f"{s}s"

def progress_loop(futures, total, start_time, update_interval=0.25):
    done = 0
    next_update = 0
    futures = set(futures)
    while futures:
        done_now = []
        for f in list(futures):
            if f.done():
                done_now.append(f)
                futures.remove(f)
        for f in done_now:
            done += 1
            yield f, done, total
        now = time.time()
        if now >= next_update or not futures:
            elapsed = now - start_time
            pct = (done / total * 100.0) if total else 100.0
            eta = (elapsed / done * (total - done)) if done else None
            eta_text = human_eta(eta)
            msg = (f"\rProcessing {done}/{total}. "
                   f"{pct:.2f}% complete. Estimated time remaining {eta_text}")
            print(msg, end='', file=sys.stderr, flush=True)
            next_update = now + update_interval
        if futures:
            time.sleep(0.05)

# --------------------------- Main pipeline ---------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Sliding-window self-alignment with minimap2; outputs genomic-coordinate PAF. Optional miniprot-guided hard-masking (whole genome) and mandatory per-segment longdust masking, plus optional TRF-mod/sdust. Also writes a filtered PAF (dedup + length + optional TSD) and FASTAs of LTR-RTs and genome-without-LTRs."
    )
    ap.add_argument("fasta", help="Genome FASTA file")
    ap.add_argument("--window", type=int, default=30000, help="Segment window size (default: 30000)")
    ap.add_argument("--step", type=int, default=5000, help="Step/slide size (default: 5000)")
    ap.add_argument("-t", "--threads", type=int, default=3, help="Threads per minimap2 process (default: 3)")
    ap.add_argument("-p", "--parallel", type=int, default=1, help="Number of segments to run simultaneously (default: 1)")

    # minimap2 tuning
    ap.add_argument("-k", "--kmer", type=int, default=13, help="minimap2 k-mer size (default: 13)")
    ap.add_argument("-w", "--minimizer-window", type=int, default=5, help="minimap2 minimizer window size (default: 5)")
    ap.add_argument("-m", "--min-occ-floor", type=int, default=25, help="minimap2 -m (default: 25)")
    ap.add_argument("-r", "--bandwidth", default="2k", help="minimap2 -r bandwidth, e.g. 2k, 5k (default: 2k)")

    # output prefix
    ap.add_argument("--out-prefix", default="aln",
                    help="Prefix for outputs: <prefix>.paf, <prefix>.filtered.paf, <prefix>.filtered.fa, <prefix>.genome_noLTRs.fa; temp dir <prefix>_out (default: aln)")
    ap.add_argument("--keep-temp", action="store_true",
                    help="Keep the temporary working directory <prefix>_out (default: delete)")

    # optional miniprot masking
    ap.add_argument("--protein", help="Protein FASTA; if provided, run miniprot and hard-mask mRNA intervals before segmentation")
    ap.add_argument("--mp-outn", type=int, default=1, help="miniprot --outn (default: 1)")
    ap.add_argument("--mp-outs", type=float, default=0.99, help="miniprot --outs (default: 0.99)")
    ap.add_argument("--mp-outc", type=float, default=0.9, help="miniprot --outc (default: 0.9)")

    # longdust tuning (mandatory masking)
    ap.add_argument("--ld-k", type=int, default=7, help="longdust -k (default: 7)")
    ap.add_argument("--ld-w", type=int, default=5000, help="longdust -w (default: 5000)")
    ap.add_argument("--ld-t", type=float, default=0.6, help="longdust -t (default: 0.6)")
    ap.add_argument("--ld-e", type=int, default=50, help="longdust -e (default: 50)")
    ap.add_argument("--ld-forward", action="store_true", help="longdust -f forward strand only")
    ap.add_argument("--ld-approx", action="store_true", help="longdust -a approximate O(Lw) algorithm")

    # ---------- NEW: TRF-mod / sdust toggles ----------
    ap.add_argument("--disable-trf", action="store_true", help="Disable TRF-mod masking (default: enabled)")
    ap.add_argument("--disable-sdust", action="store_true", help="Disable sdust masking (default: enabled)")
    ap.add_argument("--sd-t", type=int, default=16, help="sdust threshold for -t (default: 16)")

    # optional rDNA masking via Barrnap
    ap.add_argument("--enable-barrnap", action="store_true",
                    help="If set, run Barrnap on the genome and hard-mask rDNA intervals before segmentation (default: off)")
    ap.add_argument("--barrnap-kingdom", default="euk",
                    help="Barrnap --kingdom (bac, arc, mito, euk). Default: euk")
    ap.add_argument("--barrnap-reject", type=float, default=0.1,
                    help="Barrnap --reject (reject hits below this score). Default: 0.1")

    # ---------- Filtering & TSD options ----------
    ap.add_argument("--min-internal", type=int, default=100, help="Minimum internal size (default: 100)")
    ap.add_argument("--max-internal", type=int, default=25000, help="Maximum internal size (default: 25000)")
    ap.add_argument("--min-ltr", type=int, default=100, help="Minimum LTR size (default: 100)")
    ap.add_argument("--max-ltr", type=int, default=7000, help="Maximum LTR size (default: 7000)")
    ap.add_argument("--tsd-min-len", type=int, default=5, help="Minimum TSD length to accept (default: 5)")
    ap.add_argument("--tsd-max-len", type=int, default=15, help="Maximum TSD length to search (default: 15)")
    ap.add_argument("--tsd-max-subs", type=int, default=1, help="Maximum substitutions allowed between TSDs (default: 1)")
    ap.add_argument("--tsd-search-space", type=int, default=0, help="Search +/- bp around left start and right end (default: 0)")
    ap.add_argument("--disable-TSD-search", action="store_true",
                    help="Skip TSD search and filtering; TSD-related args above are ignored when set")
    # ---------- NEW: low-complexity fraction in internal ----------
    ap.add_argument("--max-internal-lc", type=float, default=1.0,
                    help="Max fraction of uppercase 'N' (low-complexity) allowed in internal region; lowercase 'n' from proteins is ignored (default: 1.0 disables this filter)")
    args = ap.parse_args()

    # Ensure minimap2 present (or build locally)
    mm2_path, mm2_installed_by_us, mm2_build_dir = ensure_minimap2_in_cwd()

    # Optional: run miniprot & hard-mask (whole genome, prior to segmentation)
    miniprot_installed_by_us = False
    miniprot_build_dir = None
    fasta_for_run = args.fasta
    if args.protein:
        try:
            mp_path, miniprot_installed_by_us, miniprot_build_dir = ensure_miniprot_in_cwd()
            mp_threads = args.parallel * args.threads
            print("Running miniprot for masking...", file=sys.stderr)
            gff_path, bed_path, masked_fa = run_miniprot_and_mask(
                args.fasta, args.protein, args.out_prefix,
                mp_threads, args.mp_outn, args.mp_outs, args.mp_outc,
                mp_path
            )
            print(f"Miniprot done. GFF: {gff_path}  Mask BED: {bed_path}  Masked FASTA: {masked_fa}", file=sys.stderr)
            fasta_for_run = masked_fa
        except Exception as e:
            print(f"WARNING: miniprot masking skipped due to error:\n{e}", file=sys.stderr)

    # Optional: run barrnap & hard-mask (whole genome, prior to segmentation)
    barrnap_installed_by_us = False
    barrnap_build_dir = None
    if args.enable_barrnap:
        try:
            bn_path, barrnap_installed_by_us, barrnap_build_dir = ensure_barrnap_in_cwd()
            if bn_path:
                bn_threads = args.parallel * args.threads
                print("Running Barrnap for rDNA masking...", file=sys.stderr)
                gff_path_b, bed_path_b, masked_fa_b = run_barrnap_and_mask(
                    fasta_for_run,  # run on whatever genome we're currently using (possibly miniprot-masked)
                    args.out_prefix,
                    bn_threads,
                    args.barrnap_kingdom,
                    args.barrnap_reject,
                    bn_path
                )
                print(f"Barrnap done. GFF: {gff_path_b}  Mask BED: {bed_path_b}  Masked FASTA: {masked_fa_b}", file=sys.stderr)
                fasta_for_run = masked_fa_b
            else:
                print("WARNING: --enable-barrnap set but Barrnap was not available; continuing without rDNA masking.", file=sys.stderr)
        except Exception as e:
            print(f"WARNING: Barrnap masking skipped due to error:\n{e}", file=sys.stderr)

    # Ensure longdust is available (mandatory)
    try:
        longdust_path, longdust_installed_by_us, longdust_build_dir = ensure_longdust_in_cwd()
    except Exception as e:
        raise SystemExit(f"ERROR: longdust is required but unavailable:\n{e}")

    # Ensure TRF-mod and sdust (best-effort; skip if disabled or build fails)
    trf_path = None
    sdust_path = None
    trf_installed_by_us = False
    sdust_installed_by_us = False
    trf_build_dir = None
    sdust_build_dir = None

    use_trf = not args.disable_trf
    use_sdust = not args.disable_sdust

    if use_trf:
        trf_path, trf_installed_by_us, trf_build_dir = ensure_trf_mod_in_cwd()
        if not trf_path:
            use_trf = False  # gracefully disable if unavailable
    if use_sdust:
        sdust_path, sdust_installed_by_us, sdust_build_dir = ensure_sdust_in_cwd()
        if not sdust_path:
            use_sdust = False

    outdir = Path(f"{args.out_prefix}_out").absolute()
    seg_dir = outdir / "segments"
    seg_mask_dir = outdir / "segments_masked"
    mask_bed_dir = outdir / "masks"        # holds .ld.bed / .trf.bed / .sd.bed and merged
    outdir.mkdir(parents=True, exist_ok=True)

    # Final artifact (raw PAF)
    out_paf_path = Path.cwd() / f"{args.out_prefix}.paf"
    with open(out_paf_path, "w"):
        pass

    # Build segment FASTAs and job list
    jobs = []
    print("Indexing FASTA and creating segments...", file=sys.stderr)
    for contig, seq in fasta_iter(fasta_for_run):
        L = len(seq)
        for s, e in segments_for_length(L, args.window, args.step):
            seg_fa = write_segment_fasta(contig, seq, s, e, str(seg_dir))
            jobs.append(SegmentJob(contig=contig, contig_len=L, start=s, end=e, fasta_path=seg_fa))
    total = len(jobs)
    if total == 0:
        print("No segments generated (check window/step vs sequence lengths).", file=sys.stderr)
        # Cleanup
        if mm2_installed_by_us:
            try: os.remove(out_paf_path.parent / "minimap2")
            except Exception: pass
            if mm2_build_dir: shutil.rmtree(mm2_build_dir, ignore_errors=True)
        if miniprot_installed_by_us:
            try: os.remove(Path.cwd() / "miniprot")
            except Exception: pass
            if miniprot_build_dir: shutil.rmtree(miniprot_build_dir, ignore_errors=True)
        if longdust_installed_by_us:
            try: os.remove(Path.cwd() / "longdust")
            except Exception: pass
            if longdust_build_dir: shutil.rmtree(longdust_build_dir, ignore_errors=True)
        if trf_installed_by_us:
            try: os.remove(Path.cwd() / "trf-mod")
            except Exception: pass
            if trf_build_dir: shutil.rmtree(trf_build_dir, ignore_errors=True)
        if sdust_installed_by_us:
            try: os.remove(Path.cwd() / "sdust")
            except Exception: pass
            if sdust_build_dir: shutil.rmtree(sdust_build_dir, ignore_errors=True)
        shutil.rmtree(outdir, ignore_errors=True)
        sys.exit(0)

    print(
        f"Running {total} self-alignments with p={args.parallel}, threads={args.threads} each "
        f"(total nominal threads = {args.parallel * args.threads}); "
        f"minimap2: -X -c --cs=short -k{args.kmer} -w{args.minimizer_window} -m{args.min_occ_floor} -r{args.bandwidth}; "
        f"longdust: -k{args.ld_k} -w{args.ld_w} -t{args.ld_t} -e{args.ld_e}"
        f"{' -f' if args.ld_forward else ''}{' -a' if args.ld_approx else ''}; "
        f"TRF-mod: {'on' if use_trf else 'off'}; sdust: {'on' if use_sdust else 'off'} (t={args.sd_t})",
        file=sys.stderr
    )

    start_time = time.time()
    errors = []

    # Launch parallel alignment runs
    with ThreadPoolExecutor(max_workers=args.parallel) as ex, open(out_paf_path, "a") as paf_out:
        futures = {
            ex.submit(
                process_segment,
                job, mm2_path, args.threads,
                args.kmer, args.minimizer_window, args.min_occ_floor, args.bandwidth,
                longdust_path,
                args.ld_k, args.ld_w, args.ld_t, args.ld_e,
                args.ld_forward, args.ld_approx,
                use_trf, trf_path,
                use_sdust, sdust_path, args.sd_t,
                seg_mask_dir, mask_bed_dir
            )
            for job in jobs
        }

        for f, done, total2 in progress_loop(futures, total, start_time, update_interval=0.25):
            try:
                adj_text = f.result()
                if adj_text:
                    paf_out.write(adj_text)
                    paf_out.flush()
            except Exception as e:
                errors.append(str(e))

    print(file=sys.stderr)
    elapsed = time.time() - start_time
    if errors:
        print("Some jobs failed:", file=sys.stderr)
        for msg in errors[:10]:
            print(f"  - {msg.splitlines()[0]}", file=sys.stderr)
        if len(errors) > 10:
            print(f"  ... and {len(errors)-10} more", file=sys.stderr)

    print(f"Done. Wrote PAF: {out_paf_path}  (elapsed {human_eta(elapsed)})", file=sys.stderr)

    # --------------------------- Filtering pass ---------------------------
    filtered_path = Path.cwd() / f"{args.out_prefix}.filtered.paf"
    print("Filtering PAF: dedup + length " + ("(no TSD)" if args.disable_TSD_search else "+ TSD") + " ...", file=sys.stderr)
    
    print("Filtering PAF: dedup + length "
          + ("(no TSD)" if args.disable_TSD_search else "+ TSD")
          + ("" if args.max_internal_lc >= 1.0 - 1e-12 else f" + internal LC≤{args.max_internal_lc}")
          + " ...", file=sys.stderr)
    raw_total, dedup_kept, final_kept = filter_paf(
        paf_path=str(out_paf_path),
        genome_fa=args.fasta,  # use ORIGINAL genome for TSD extraction
        out_path=str(filtered_path),
        min_internal=args.min_internal,
        max_internal=args.max_internal,
        min_ltr=args.min_ltr,
        max_ltr=args.max_ltr,
        tsd_min_len=args.tsd_min_len,
        tsd_max_len=args.tsd_max_len,
        tsd_max_subs=args.tsd_max_subs,
        tsd_search_space=args.tsd_search_space,
        disable_tsd_search=args.disable_TSD_search,
        seg_mask_dir=str(seg_mask_dir),             # ### NEW: where the masked segments live
        max_internal_lc=float(args.max_internal_lc) # ### NEW
    )
    print(f"Filtering stats: raw={raw_total}, after_dedup={dedup_kept}, kept_after_filters={final_kept}", file=sys.stderr)
    print(f"Wrote filtered PAF: {filtered_path}", file=sys.stderr)

    # ---- LTR-RT FASTA + genome-without-LTRs ----
    ltr_fa_path = Path.cwd() / f"{args.out_prefix}.filtered.fa"
    intervals_by_contig = emit_ltr_fasta_and_collect_intervals(
        filtered_paf_path=str(filtered_path),
        genome_fa=args.fasta,
        out_ltr_fa=str(ltr_fa_path)
    )

    genome_no_ltrs_path = Path.cwd() / f"{args.out_prefix}.genome_noLTRs.fa"
    merged_by_contig = {c: merge_intervals(v) for c, v in intervals_by_contig.items()}
    write_genome_without_intervals(
        genome_fa=args.fasta,
        intervals_by_contig=merged_by_contig,
        out_path=str(genome_no_ltrs_path)
    )

    # Cleanup temp dirs
    if not args.keep_temp:
        try:
            shutil.rmtree(outdir, ignore_errors=True)
        except Exception:
            pass
    else:
        print(f"Keeping temp directory: {outdir}", file=sys.stderr)

    # If we installed minimap2, remove it and its build dir
    if mm2_installed_by_us:
        try: os.remove(Path.cwd() / "minimap2")
        except Exception: pass
        if mm2_build_dir: shutil.rmtree(mm2_build_dir, ignore_errors=True)

    # If we installed miniprot, remove it and its build dir
    if miniprot_installed_by_us:
        try: os.remove(Path.cwd() / "miniprot")
        except Exception: pass
        if miniprot_build_dir: shutil.rmtree(miniprot_build_dir, ignore_errors=True)

    # If we installed longdust, remove it and its build dir
    if longdust_installed_by_us:
        try: os.remove(Path.cwd() / "longdust")
        except Exception: pass
        if longdust_build_dir: shutil.rmtree(longdust_build_dir, ignore_errors=True)

    # If we installed TRF-mod, remove it and its build dir
    if trf_installed_by_us:
        try: os.remove(Path.cwd() / "trf-mod")
        except Exception: pass
        if trf_build_dir: shutil.rmtree(trf_build_dir, ignore_errors=True)

    # If we installed sdust, remove it and its build dir
    if sdust_installed_by_us:
        try: os.remove(Path.cwd() / "sdust")
        except Exception: pass
        if sdust_build_dir: shutil.rmtree(sdust_build_dir, ignore_errors=True)

    # If we installed Barrnap, remove its build dir
    if barrnap_installed_by_us and barrnap_build_dir:
        try: shutil.rmtree(barrnap_build_dir, ignore_errors=True)
        except Exception: pass

if __name__ == "__main__":
    main()

