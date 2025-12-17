#!/usr/bin/env python3
"""
LTR-RT end-to-end pipeline runner

Implements the user-specified steps:

(1a) REQUIRED: LTR_HARVEST_parallel (LTRharvest + LTRfinder)
(1b) OPTIONAL: minimap2_LTR discovery + merge
(2)  WFA filtering
(3)  Kmer2LTR filtering + divergence
(4)  Compile final candidates + TEsorter classification
(5)  Purge overlaps

Arguments:
  (1)  --genome / -g : genome fasta (required)
  (2)  --protein / -p : protein fasta (optional)
  (3)  --threads / -t : threads (default 10)
  (4)  --script-dir / -s : directory containing helper scripts:
       - minimap2_LTR.py
       - paf_to_scn.py
       - wfa_LTR (executable)
       - add_tsd_to_div.py
       - kmer2ltr_fa.py
       - LTR_overlap_purger.py
  (5)  --keep-temp : keep intermediate files (default: delete junk listed)
  (6)  --on-exist : how to handle existing outputs: overwrite|skip (default overwrite)
  (7)  --run-optional : run the optional minimap2 discovery path and merge with LTRharvest

The script will `git clone` these tools if missing:
  - https://github.com/cwb14/LTR_HARVEST_parallel.git
  - https://github.com/cwb14/Kmer2LTR.git
  - https://github.com/cwb14/TEsorter.git

External tools you should have on PATH (the pipeline will check and warn):
  - bedtools, awk, sort, uniq, grep
  - barrnap, sdust, trf, longdust (used inside LTR_HARVEST_parallel if enabled)
  - minimap2 (if you run the optional path)
  - Python 3
  
"""

import argparse
import os
import sys
import shutil
import subprocess
from pathlib import Path
from textwrap import dedent

# --------------------------- helpers ---------------------------------

def run_cmd(cmd, shell=False, check=True, env=None):
    print(f"[cmd] {cmd if isinstance(cmd, str) else ' '.join(cmd)}", flush=True)
    result = subprocess.run(cmd, shell=shell, check=check, env=env)
    return result.returncode


def ensure_repo(url: str, dest: Path):
    if dest.exists():
        print(f"[ok] Repo exists: {dest}")
        return
    print(f"[get] Cloning {url} -> {dest}")
    run_cmd(["git", "clone", url, str(dest)])


def need(prog: str, hard=False):
    p = shutil.which(prog)
    if p:
        print(f"[ok] Found dependency: {prog} -> {p}")
        return True
    msg = f"[warn] Missing dependency on PATH: {prog}"
    if hard:
        print(msg + " (required)"); sys.exit(1)
    else:
        print(msg + " (may be required by some steps)")
        return False


def exists_and_handle(path: Path, on_exist: str):
    if path.exists() and on_exist == "skip":
        print(f"[skip] Exists (on-exist=skip): {path}")
        return True
    # overwrite: nothing to do; downstream tools will overwrite
    return False


# --------------------------- main pipeline ---------------------------

def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("-g", "--genome", required=True, help="Genome fasta")
    ap.add_argument("-p", "--protein", default=None, help="Protein fasta (optional)")
    ap.add_argument("-t", "--threads", type=int, default=10, help="Threads")
    ap.add_argument("-s", "--script-dir", required=True,
                    help="Directory containing helper scripts (minimap2_LTR.py, paf_to_scn.py, wfa_LTR, add_tsd_to_div.py, kmer2ltr_fa.py, LTR_overlap_purger.py)")
    ap.add_argument("--keep-temp", action="store_true", help="Keep intermediate junk files")
    ap.add_argument("--on-exist", choices=["overwrite", "skip"], default="overwrite",
                    help="How to handle existing outputs")
    ap.add_argument("--run-optional", action="store_true",
                    help="Run the OPTIONAL minimap2 discovery and merge with LTRharvest")
    args = ap.parse_args()

    GENOME = Path(args.genome).resolve()
    if not GENOME.exists():
        print(f"[err] Genome not found: {GENOME}"); sys.exit(1)
    PROTEIN = Path(args.protein).resolve() if args.protein else None
    if PROTEIN and not PROTEIN.exists():
        print(f"[err] Protein not found: {PROTEIN}"); sys.exit(1)
    THREADS = str(args.threads)
    SCRIPT_DIR = Path(args.script_dir).resolve()
    if not SCRIPT_DIR.exists():
        print(f"[err] SCRIPT_DIR not found: {SCRIPT_DIR}"); sys.exit(1)

    # Repos to fetch
    repo_ltrharv = Path("LTR_HARVEST_parallel")
    repo_kmer2ltr = Path("Kmer2LTR")
    repo_tesorter = Path("TEsorter")

    ensure_repo("https://github.com/cwb14/LTR_HARVEST_parallel.git", repo_ltrharv)
    ensure_repo("https://github.com/cwb14/Kmer2LTR.git", repo_kmer2ltr)
    ensure_repo("https://github.com/cwb14/TEsorter.git", repo_tesorter)

    # PATH checks (best-effort)
    need("git", hard=True)
    need("bedtools")
    need("awk")
    need("sort")
    need("uniq")
    need("grep")
    if args.run_optional:
        need("minimap2", hard=False)
    # The following are used within LTR_HARVEST_parallel when flags are given
    # We'll warn but not hard fail here.
    for prog in ["barrnap", "sdust", "trf", "longdust"]:
        need(prog, hard=False)

    # Filenames derived from genome basename
    gbase = GENOME.name
    merged_scn = Path(f"{gbase}_merged.scn")
    merged_scn_wfa = Path(f"{gbase}_merged.scn.wfa")
    merged_pass = Path(f"{gbase}_merged.scn.wfa.pass")
    k2ltr_prefix = Path(f"{gbase}_Kmer2LTR")
    k2ltr_tsd = Path(f"{gbase}_Kmer2LTR_TSD")
    k2ltr_tsd_fa = Path(f"{gbase}_Kmer2LTR_TSD.fa")
    cls_prefix = Path(f"{gbase}_Kmer2LTR_TSD.fa.rexdb-plant")

    # ---------------- (1a) REQUIRED: LTRharvest/finder via LTR_HARVEST_parallel ----------------
    print("\n=== (1a) REQUIRED: LTR_HARVEST_parallel (LTRharvest + LTRfinder) ===")

    ltrharvest_py = repo_ltrharv / "LTR_HARVEST_parallel3.py"
    if not ltrharvest_py.exists():
        print(f"[err] Missing script: {ltrharvest_py}"); sys.exit(1)

    ltrharvest_args = (
        '-similar 70 -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 8 -seqids yes -overlaps all -seed 10'
#       '-similar 70 -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 8 -seqids yes -overlaps best -seed 20 -xdrop 5'
    )
    ltrfinder_args = (
        '-w 2 -D 15000 -d 500 -L 7000 -l 100 -S 4 -p 10 -g 150'
#       '-w 2 -D 15000 -d 500 -L 7000 -l 100 -S 6 -p 15 -g 50'
    )
    sweep_file = repo_ltrharv / "sweep.tsv"

    base_cmd = [
        sys.executable, str(ltrharvest_py),
        "-seq", str(GENOME),
        "--overlap", "30000",
        "-t", THREADS,
#       "--barrnap"      
#       "--barrnap", "--sdust", "--trf", "--longdust" # Based on simulation, using all "--barrnap", "--sdust", "--trf", "--longdust", "--protein" leads to loss. I've not tested trying each independelntly...only all or nothing. 
    ]
    # conditionally gene-mask/protein args
    if PROTEIN:
        base_cmd += ["--gene-mask", "--mp-outn", "1000", "--mp-outs", "0.99", "--mp-outc", "0.1", "--protein", str(PROTEIN)]
    base_cmd += [
        "--ltrharvest-args", f"{ltrharvest_args}",
        "--run-ltrfinder", # Tab it out for just harvest.
        "--ltrfinder-args", f"{ltrfinder_args}", # Tab it out for just harvest.
        "--sweep-file", str(sweep_file) # Tab it out to skip sweep.
    ]

    # Skip if main outputs exist and on-exist=skip
    # The LTR_HARVEST_parallel creates {GENOME}.harvest.combine.gff3 and .scn, etc.
    harvest_scn = Path(f"{gbase}.harvest.combine.scn")
    if exists_and_handle(harvest_scn, args.on_exist):
        print("[note] Using existing LTRharvest/finder results.")
    else:
        run_cmd(base_cmd)

    # Clean up junk from (1a)
    if not args.keep_temp:
        junk_1a = [
            f"{gbase}.genic.mask.bed",
            f"{gbase}.genic_masked.fa",
            f"{gbase}.genic.gff",
            f"{gbase}.harvest.combine.gff3"
        ]
        for j in junk_1a:
            p = Path(j)
            if p.exists():
                try:
                    p.unlink()
                    print(f"[rm] {p}")
                except Exception as e:
                    print(f"[warn] Could not remove {p}: {e}")

    # ---------------- (1b) OPTIONAL: Minimap2 discovery + merge ----------------
    if args.run_optional:
        print("\n=== (1b) OPTIONAL: minimap2_LTR discovery + merge ===")
        minimap2_LTR_py = SCRIPT_DIR / "minimap2_LTR.py"
        if not minimap2_LTR_py.exists():
            print(f"[err] Missing script: {minimap2_LTR_py}")
            sys.exit(1)

        mm_prefix = f"{gbase}_mm"
        mm_cmd = [
            sys.executable, str(minimap2_LTR_py), str(GENOME),
            "-t", "1", "-p", THREADS,
#            "--mp-outn", "1000", "--mp-outs", "0.99", "--mp-outc", "0.1",
            "--mp-outn", "1000", "--mp-outs", "0.98", "--mp-outc", "0.1",
            "--out-prefix", mm_prefix,
            "--min-internal", "100", "--max-internal", "30000",
            "--min-ltr", "100", "--max-ltr", "15000",
            "--disable-TSD-search",
            "--sd-t", "12",
            "--enable-barrnap"
        ]
        if PROTEIN:
            mm_cmd += ["--protein", str(PROTEIN)]

        # Skip minimap2_LTR if the merged SCN already exists (since it implies completion)
        if merged_scn.exists() and args.on_exist == "skip":
            print(f"[skip] minimap2_LTR (merged SCN already present: {merged_scn})")
        else:
            mm_filtered_paf = Path(f"{mm_prefix}.filtered.paf")
            if exists_and_handle(mm_filtered_paf, args.on_exist):
                print("[skip] minimap2_LTR (existing filtered PAF).")
            else:
                run_cmd(mm_cmd)

        # Remove junk for the minimap block
        if not args.keep_temp:
            junk_1b = [
                f"{mm_prefix}.gff",
                f"{mm_prefix}.mask.bed",
                f"{mm_prefix}.hardmasked.fa",
                f"{mm_prefix}.barrnap.gff3",
                f"{mm_prefix}.barrnap.mask.bed",
                f"{mm_prefix}.barrnap.hardmasked.fa",
                f"{mm_prefix}.paf",
                f"{mm_prefix}.filtered.fa",
                f"{mm_prefix}.genome_noLTRs.fa",
            ]
            for j in junk_1b:
                p = Path(j)
                if p.exists():
                    try:
                        p.unlink()
                        print(f"[rm] {p}")
                    except Exception as e:
                        print(f"[warn] Could not remove {p}: {e}")

        # Convert PAF -> SCN
        paf_to_scn = SCRIPT_DIR / "paf_to_scn.py"
        if not paf_to_scn.exists():
            print(f"[err] Missing script: {paf_to_scn}"); sys.exit(1)
        mm_scn = Path(f"{mm_prefix}.filtered.scn")

        if not exists_and_handle(mm_scn, args.on_exist):
            # Convert PAF -> SCN (single clean redirect, no stdout clutter)
            run_cmd(f"{sys.executable} {paf_to_scn} {mm_prefix}.filtered.paf > {mm_prefix}.filtered.scn", shell=True)

        # Remove PAF
        p = Path(f"{mm_prefix}.filtered.paf")
        if p.exists() and not args.keep_temp:
            try:
                p.unlink()
                print(f"[rm] {p}")
            except Exception as e:
                print(f"[warn] Could not remove {p}: {e}")

        # Merge harvest SCN + minimap SCN
        if not exists_and_handle(merged_scn, args.on_exist):
            run_cmd(f"cat {gbase}.harvest.combine.scn {mm_prefix}.filtered.scn > {gbase}_merged.scn", shell=True)
    else:
        # Without optional path, use harvest SCN directly
        print("\n=== (1b) OPTIONAL path skipped; using harvest SCN directly ===")
        if not harvest_scn.exists():
            print(f"[err] Expected SCN not found: {harvest_scn}"); sys.exit(1)
        if not exists_and_handle(merged_scn, args.on_exist):
            shutil.copyfile(harvest_scn, merged_scn)
            print(f"[cp] {harvest_scn} -> {merged_scn}")

    # ---------------- (2) WFA filtering ----------------
    print("\n=== (2) WFA filtering ===")
    wfa_exec = SCRIPT_DIR / "wfa_LTR"
    if not wfa_exec.exists():
        print(f"[err] Missing executable: {wfa_exec}"); sys.exit(1)
    if not os.access(wfa_exec, os.X_OK):
        print(f"[err] Not executable: {wfa_exec}"); sys.exit(1)

    # Run WFA on merged SCN
    if not exists_and_handle(merged_scn_wfa, args.on_exist):
        cmd = dedent(f"""
        {wfa_exec} <(grep -v '^#' {merged_scn} | sort | uniq) {GENOME} --win-pairs 34 -x 5 -O 37 -E 7 -o 41 -e 3 -f 45 --vic-out 10 --vic-in 5 --match-thresh 18 --allowed-subs 3 --allowed-gaps 0 -A -w 200 --thresh-high 0.70 --thresh-low 0.5 > {merged_scn_wfa} 2> /dev/null

        # Based on grid searching for best parameters, I think:
        # 'match_thresh 10' to 'match_thresh 18'. Pick 14 or 18. 
        # "allowed_gaps 0" to "allowed_gaps 1". Not sure our simulation is great for gaps, so I should maybe consider raising to 1. 
        # "allowed_subs 1" to "allowed_subs 5".
        # "win_pairs 30" to "win_pairs 34".
        # "thresh_high 0.7" to "thresh_high 0.75".
        # "thresh_low 0.5" to "thresh_low 0.6"
        # "f 40" to "f 45".
        # "x 3" to "x 7".
        # "O 20" to "O 50" & "o 20" to "o 50".
        # "E 3" to "E 9" & "e 1" to "e 5".

        """).strip()
        # Need bash for process substitution
        run_cmd(["bash", "-lc", cmd], check=True)

    # Process alignments -> PASS + basic length filters, columns 3-9
    if not exists_and_handle(merged_pass, args.on_exist):
        cmd = dedent(f"""
        cat {merged_scn_wfa} | grep PASS | tr -s ' ' '\t' | awk -F'\t' '($1!="NA" && $2!="NA" && $3!="NA" && $4!="NA" && $5!="NA" && ($4-$3)>=100 && ($5-$4)>=100)' | sort -u | cut -f 3-9 > {merged_pass}
        """).strip()
        run_cmd(["bash", "-lc", cmd], check=True)

    # ---------------- (3) Kmer2LTR filtering + divergence ----------------
    print("\n=== (3) Kmer2LTR filtering + divergence ===")
    # Make .len and fasta for Kmer2LTR
    pass_len = Path(f"{gbase}_merged.scn.wfa.pass.len")
    pass_bed = Path(f"{gbase}_merged.scn.wfa.pass.bed")
    pass_fa = Path(f"{gbase}_merged.scn.wfa.pass.fa")

    if not exists_and_handle(pass_len, args.on_exist):
        cmd = f"awk '{{d1=$3-$2; d2=$5-$4; printf \"%s:%d-%d\\t%d\\n\", $1, $2, $5, (d1>d2?d1:d2)}}' {merged_pass} > {pass_len}"
        run_cmd(["bash", "-lc", cmd], check=True)

    if not exists_and_handle(pass_bed, args.on_exist):
        cmd = f"awk 'BEGIN{{OFS=\"\\t\"}} {{print $1, $2, $5}}' {merged_pass} > {pass_bed}"
        run_cmd(["bash", "-lc", cmd], check=True)

    if not exists_and_handle(pass_fa, args.on_exist):
        run_cmd(["bedtools", "getfasta", "-fi", str(GENOME), "-bed", str(pass_bed), "-fo", str(pass_fa)], check=True)

    # Remove .bed (junk)
    if pass_bed.exists() and not args.keep_temp:
        try:
            pass_bed.unlink(); print(f"[rm] {pass_bed}")
        except Exception as e:
            print(f"[warn] Could not remove {pass_bed}: {e}")

    # Run Kmer2LTR
    if not exists_and_handle(k2ltr_prefix.with_suffix(".summary"), args.on_exist):
        run_cmd([
            sys.executable, "./Kmer2LTR/Kmer2LTR.py",
#            "-k", # Keep temps
            "-D", str(pass_len),
            "-p", THREADS,
            "--max-win-overdisp", "6",
            "-i", str(pass_fa),
            "-o", str(k2ltr_prefix)
        ], check=True)

    # ---------------- (4) Compile final candidates + TEsorter classification ----------------
    print("\n=== (4) Compile final candidates + classification ===")
    add_tsd_py = SCRIPT_DIR / "add_tsd_to_div.py"
    k2ltr_fa_py = SCRIPT_DIR / "kmer2ltr_fa.py"
    if not add_tsd_py.exists() or not k2ltr_fa_py.exists():
        print(f"[err] Missing one or more scripts in {SCRIPT_DIR}: add_tsd_to_div.py, kmer2ltr_fa.py")
        sys.exit(1)

    # add_tsd_to_div -> *_Kmer2LTR_TSD   (single redirected call; avoid printing lists)
    if not exists_and_handle(k2ltr_tsd, args.on_exist):
        run_cmd(f"{sys.executable} {add_tsd_py} {k2ltr_prefix} {merged_pass} > {k2ltr_tsd}", shell=True)

    # kmer2ltr_fa -> *_Kmer2LTR_TSD.fa
    if not exists_and_handle(k2ltr_tsd_fa, args.on_exist):
        run_cmd(f"{sys.executable} {k2ltr_fa_py} {GENOME} {k2ltr_tsd} > {k2ltr_tsd_fa}", shell=True)

    # Remove Kmer2LTR working files (junk)
    if not args.keep_temp:
        for j in [k2ltr_prefix.with_suffix(".summary"), k2ltr_prefix, Path("kmer2ltr_density.pdf"), Path(f"{gbase}_Kmer2LTR.log")]:
            p = Path(str(j))
            if p.exists():
                try:
                    if p.is_dir():
                        shutil.rmtree(p)
                    else:
                        p.unlink()
                    print(f"[rm] {p}")
                except Exception as e:
                    print(f"[warn] Could not remove {p}: {e}")

    # TEsorter classification
    print("[info] Running TEsorter (rexdb-plant)...")
    env = os.environ.copy()
    env["PYTHONPATH"] = str(repo_tesorter.resolve())
    run_cmd([
        sys.executable, "-m", "TEsorter", str(k2ltr_tsd_fa),
        "-db", "rexdb-plant", "-p", THREADS, "-cov", "10", "-eval", "1e-2", "-rule", "70-30-80"
    ], check=True, env=env)

    # Join classifications to *_Kmer2LTR_TSD -> *_Kmer2LTR_TSD_class
    cls_tsv = Path(f"{k2ltr_tsd_fa}.rexdb-plant.cls.tsv")
    out_class = Path(f"{gbase}_Kmer2LTR_TSD_class")
    awk_join = dedent(f"""
    awk -v OFS="\t" '
    FNR==NR {{
      k = $1
      $1 = ""
      sub(/^[ \t]+/, "", $0)
      map[k] = $0
      cnt = NF - 1
      if (cnt > maxc) maxc = cnt
      next
    }}
    {{
      k = $1
      if (k in map) {{
        n = split(map[k], tmp, /[ \t]+/)
        printf "%s", $0
        for (i=1; i<=n; i++) printf "\t%s", tmp[i]
        for (i=n+1; i<=maxc; i++) printf "\tNA"
        print ""
      }} else {{
        printf "%s", $0
        for (i=1; i<=maxc; i++) printf "\tNA"
        print ""
      }}
    }}
    ' {cls_tsv} {k2ltr_tsd} > {out_class}
    """).strip()
    run_cmd(["bash", "-lc", awk_join], check=True)

    # Remove TEsorter junk (and additional outputs) *after* join
    if not args.keep_temp:
        ts_junk = [
            f"{k2ltr_tsd_fa}.rexdb-plant.domtbl",
            f"{k2ltr_tsd_fa}.rexdb-plant.dom.gff3",
            f"{k2ltr_tsd_fa}.rexdb-plant.dom.faa",
            f"{k2ltr_tsd_fa}.rexdb-plant.dom.tsv",
            "pass1_classified.mmi",
            f"{k2ltr_tsd_fa}.rexdb-plant.cls.lib",
            f"{k2ltr_tsd_fa}.rexdb-plant.cls.lib.rexdb-plant.cls.pep",
            f"{k2ltr_tsd_fa}.rexdb-plant.cls.pep",
            f"{k2ltr_tsd_fa}.rexdb-plant.cls.tsv",
            str(k2ltr_tsd_fa),          # fasta used for classification
            str(k2ltr_tsd),             # table prior to classification join
        ]
        for j in ts_junk:
            p = Path(j)
            if p.exists():
                try:
                    p.unlink(); print(f"[rm] {p}")
                except Exception as e:
                    print(f"[warn] Could not remove {p}: {e}")

    # ---------------- (5) Purge overlaps ----------------
    print("\n=== (5) Purge overlaps ===")
    purger_py = SCRIPT_DIR / "LTR_overlap_purger.py"
    if not purger_py.exists():
        print(f"[err] Missing script: {purger_py}"); sys.exit(1)
    out_purge = Path(f"{gbase}_Kmer2LTR_TSD_class_purge")
    if not exists_and_handle(out_purge, args.on_exist):
        # Single redirected call; avoid printing lists to terminal
        run_cmd(f"{sys.executable} {purger_py} {out_class} --ov-longest > {out_purge}", shell=True)

    # Final FASTA of *purged* calls
    out_purge_fa = Path(f"{out_purge}.fa")
    if not exists_and_handle(out_purge_fa, args.on_exist):
        run_cmd(f"{sys.executable} {k2ltr_fa_py} {GENOME} {out_purge} > {out_purge_fa}", shell=True)

    # Extra cleanup of big intermediates created along the way
    if not args.keep_temp:
        extra_junk = [
            Path(f"{gbase}_merged.scn.wfa.pass.fa"),
            Path(f"{gbase}_merged.scn.wfa.pass.len"),
            Path(f"{gbase}.fai"),
            Path("temp"),  # created by some tools (e.g., Kmer2LTR)
        ]
        for j in extra_junk:
            if j.exists():
                try:
                    if j.is_dir():
                        shutil.rmtree(j)
                    else:
                        j.unlink()
                    print(f"[rm] {j}")
                except Exception as e:
                    print(f"[warn] Could not remove {j}: {e}")

    print("\n[done] Pipeline completed successfully.")
    print(f"[out] Final merged candidates (pre-class): {k2ltr_tsd}")
    print(f"[out] Final classified:                 {out_class}")
    print(f"[out] Final classified (purged):        {out_purge}")
    print(f"[out] Purged FASTA:                     {out_purge_fa}")


if __name__ == "__main__":
    main()
