#!/usr/bin/env bash
set -euo pipefail

# consensus_library.sh
# Cluster intact TEs with mmseqs easy-cluster, align each cluster with MAFFT,
# optionally trim the alignment with trimAl, then generate a consensus per cluster
# either via:
#   (A) hmmbuild + hmmemit (default), or
#   (B) EMBOSS consambig

# Create input using ltrharvest.py.
# Merge all LTR-RTs. 
cat ref*.ltrharvest.full_length.dedup.fa.rexdb-plant.cls.lib.fa | sed '/^>/! { s/[^ATCGatcg]//g; /^$/d }' > ref_all.fa
# Clean headers for this consensus tool. 
sed '/^>/ s/#.*//' ref_all.fa | sed 's/:/_/g' > LTR_intact_TEs.fa
# I merged 'all_intact_TEs.fa' with '$.fa.mod.EDTA.intact.fa' from EDTA to create 'all_intact_TEs.fa'.

# I run it:
bash ./consensus_library.sh -i all_intact_TEs.fa -o TE_consensus_out -p all_intact_TEs -t 20 --min-seq-id 0.80 -c 0.80 --cov-mode 0 --cluster-mode 0 --mafft-mode auto --cons-tool consambig
# '--cons-tool consambig' converts to IUPAC DNA Ambiguity Codes. The HMM approach just picks the best base. 

usage() {
  cat <<'EOF'
Usage:
  consensus_library.sh -i intact.fa -o outdir -p PREFIX -t THREADS [options]

Required:
  -i  Input FASTA of intact LTR-RTs
  -o  Output directory
  -p  Prefix for outputs
Optional:
  -t  Threads (default: all available cores)

MMseqs clustering knobs (optional):
  --min-seq-id FLOAT   (default: 0.80)
  -c FLOAT             coverage fraction (default: 0.80)
  --cov-mode INT       (default: 0)
  -e DOUBLE            E-value (default: 1e-3)
  --cluster-mode INT   (default: 0)

MAFFT knobs (optional):
  --mafft-mode STR     auto|linsi|ginsi|einsi  (default: auto)
  --mafft-extra "..."  extra args passed to mafft (default: "")

Consensus options:
  --cons-tool STR      hmmemit|consambig (default: hmmemit)
  --fancy              (hmmemit only) use hmmemit -C instead of -c
  --minl FLOAT         (hmmemit -C only) hmmemit --minl
  --minu FLOAT         (hmmemit -C only) hmmemit --minu

Trim options:
  --trim               enable trimAl after MAFFT
  --trimal-args "..."  args passed to trimal (default: "-automated1")

Example:
  consensus_ltrrt_library.sh -i intact_LTRRT.fa -o TE_consensus -p LTRRT -t 32 \
    --min-seq-id 0.80 -c 0.80 --cov-mode 0 -e 1e-5 --cluster-mode 0 \
    --mafft-mode auto \
    --trim --trimal-args "-automated1" \
    --cons-tool consambig
EOF
}

# -----------------------------
# Defaults
# -----------------------------
threads="$(nproc || echo 1)"

mm_min_seq_id="0.80"
mm_cov="0.80"
mm_cov_mode="0"
mm_evalue="1e-3"
mm_cluster_mode="0"

mafft_mode="auto"
mafft_extra=""

cons_tool="hmmemit"   # hmmemit|consambig

use_fancy="0"
fancy_minl=""
fancy_minu=""

do_trim="0"
trimal_args="-automated1"

in_fa=""
outdir=""
prefix=""

# -----------------------------
# Arg parsing
# -----------------------------
if [[ $# -eq 0 ]]; then usage; exit 1; fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i) in_fa="$2"; shift 2;;
    -o) outdir="$2"; shift 2;;
    -p) prefix="$2"; shift 2;;
    -t) threads="$2"; shift 2;;

    --min-seq-id) mm_min_seq_id="$2"; shift 2;;
    -c) mm_cov="$2"; shift 2;;
    --cov-mode) mm_cov_mode="$2"; shift 2;;
    -e) mm_evalue="$2"; shift 2;;
    --cluster-mode) mm_cluster_mode="$2"; shift 2;;

    --mafft-mode) mafft_mode="$2"; shift 2;;
    --mafft-extra) mafft_extra="$2"; shift 2;;

    --cons-tool) cons_tool="$2"; shift 2;;

    --fancy) use_fancy="1"; shift 1;;
    --minl) fancy_minl="$2"; shift 2;;
    --minu) fancy_minu="$2"; shift 2;;

    --trim) do_trim="1"; shift 1;;
    --trimal-args) trimal_args="$2"; shift 2;;

    -h|--help) usage; exit 0;;
    *) echo "ERROR: Unknown option: $1" >&2; usage; exit 1;;
  esac
done

# -----------------------------
# Validate inputs + deps
# -----------------------------
[[ -n "$in_fa" && -f "$in_fa" ]] || { echo "ERROR: input FASTA missing: -i" >&2; exit 1; }
[[ -n "$outdir" ]] || { echo "ERROR: output dir missing: -o" >&2; exit 1; }
[[ -n "$prefix" ]] || { echo "ERROR: prefix missing: -p" >&2; exit 1; }

command -v mmseqs >/dev/null 2>&1 || { echo "ERROR: mmseqs not found in PATH" >&2; exit 1; }
command -v mafft  >/dev/null 2>&1 || { echo "ERROR: mafft not found in PATH" >&2; exit 1; }

# Consensus tool deps (conditional)
case "$cons_tool" in
  hmmemit)
    command -v hmmbuild >/dev/null 2>&1 || { echo "ERROR: hmmbuild not found (HMMER)" >&2; exit 1; }
    command -v hmmemit  >/dev/null 2>&1 || { echo "ERROR: hmmemit not found (HMMER)" >&2; exit 1; }
    ;;
  consambig)
    command -v consambig >/dev/null 2>&1 || { echo "ERROR: consambig not found (EMBOSS)" >&2; exit 1; }
    ;;
  *)
    echo "ERROR: --cons-tool must be hmmemit or consambig" >&2
    exit 1
    ;;
esac

# Trim deps (conditional)
if [[ "$do_trim" == "1" ]]; then
  command -v trimal >/dev/null 2>&1 || { echo "ERROR: trimal not found but --trim was requested" >&2; exit 1; }
fi

have_seqkit="0"
if command -v seqkit >/dev/null 2>&1; then
  have_seqkit="1"
else
  echo "WARNING: seqkit not found. Falling back to a slower awk FASTA extractor." >&2
fi

mkdir -p "$outdir"
tmpdir="$outdir/tmp"
mkdir -p "$tmpdir"
workdir="$outdir/work"
mkdir -p "$workdir"
clustdir="$outdir/clusters"
mkdir -p "$clustdir"

cluster_prefix="$outdir/$prefix"
cluster_tsv="${cluster_prefix}_cluster.tsv"
consensus_fa="$outdir/${prefix}.consensus.fa"
log="$outdir/${prefix}.log"
: > "$log"

echo "[INFO] Input:        $in_fa" | tee -a "$log"
echo "[INFO] Outdir:       $outdir" | tee -a "$log"
echo "[INFO] Prefix:       $prefix" | tee -a "$log"
echo "[INFO] Threads:      $threads" | tee -a "$log"
echo "[INFO] MMseqs:       --min-seq-id $mm_min_seq_id  -c $mm_cov  --cov-mode $mm_cov_mode  -e $mm_evalue  --cluster-mode $mm_cluster_mode" | tee -a "$log"
echo "[INFO] MAFFT mode:   $mafft_mode" | tee -a "$log"
echo "[INFO] MAFFT extra:  $mafft_extra" | tee -a "$log"
echo "[INFO] TrimAl:       $([[ "$do_trim" == "1" ]] && echo "ON ($trimal_args)" || echo "OFF")" | tee -a "$log"
echo "[INFO] Consensus:    $cons_tool" | tee -a "$log"
if [[ "$cons_tool" == "hmmemit" ]]; then
  if [[ "$use_fancy" == "1" ]]; then
    echo "[INFO] hmmemit:      -C (fancy) ${fancy_minl:+--minl $fancy_minl} ${fancy_minu:+--minu $fancy_minu}" | tee -a "$log"
  else
    echo "[INFO] hmmemit:      -c (majority-rule consensus)" | tee -a "$log"
  fi
fi

# -----------------------------
# Helper: FASTA extract by ID list
# -----------------------------
extract_fasta_by_ids() {
  local fasta="$1"
  local idlist="$2"
  local outfa="$3"

  if [[ "$have_seqkit" == "1" ]]; then
    seqkit grep -n -f "$idlist" "$fasta" > "$outfa"
  else
    awk -v IDS="$idlist" '
      BEGIN{
        while((getline<IDS)>0){gsub(/\r/,""); if($0!="") id[$0]=1}
      }
      /^>/{
        hdr=$0
        sub(/^>/,"",hdr)
        split(hdr,a,/[\t ]+/)
        keep = (a[1] in id)
      }
      { if(keep) print }
    ' "$fasta" > "$outfa"
  fi
}

# -----------------------------
# Helper: choose MAFFT command
# -----------------------------
run_mafft() {
  local infa="$1"
  local outaln="$2"

  case "$mafft_mode" in
    auto)
      mafft --auto --reorder --thread -1 $mafft_extra "$infa" > "$outaln"
      ;;
    linsi)
      mafft --maxiterate 1000 --localpair --reorder --thread -1 $mafft_extra "$infa" > "$outaln"
      ;;
    ginsi)
      mafft --maxiterate 1000 --globalpair --reorder --thread -1 $mafft_extra "$infa" > "$outaln"
      ;;
    einsi)
      mafft --maxiterate 1000 --genafpair --reorder --thread -1 $mafft_extra "$infa" > "$outaln"
      ;;
    *)
      echo "ERROR: unknown --mafft-mode: $mafft_mode (use auto|linsi|ginsi|einsi)" >&2
      exit 1
      ;;
  esac
}

# -----------------------------
# Helper: optional trimAl
# -----------------------------
maybe_trimal() {
  local inaln="$1"
  local outaln="$2"

  if [[ "$do_trim" == "1" ]]; then
    # trimAl usage: trimal -in <aln> -out <trimmed> [options]
    trimal -in "$inaln" -out "$outaln" $trimal_args
  else
    cp -f "$inaln" "$outaln"
  fi
}

# -----------------------------
# 1) Cluster with mmseqs
# -----------------------------
echo "[INFO] Running mmseqs easy-cluster..." | tee -a "$log"
mmseqs easy-cluster \
  "$in_fa" \
  "$cluster_prefix" \
  "$tmpdir" \
  --threads "$threads" \
  --min-seq-id "$mm_min_seq_id" \
  -c "$mm_cov" \
  --cov-mode "$mm_cov_mode" \
  -e "$mm_evalue" \
  --cluster-mode "$mm_cluster_mode" \
  >> "$log" 2>&1

[[ -s "$cluster_tsv" ]] || { echo "ERROR: expected cluster TSV not found or empty: $cluster_tsv" >&2; exit 1; }
echo "[INFO] mmseqs clusters: $cluster_tsv" | tee -a "$log"

# -----------------------------
# 2) Build per-cluster consensus
# -----------------------------
echo "[INFO] Preparing cluster member lists..." | tee -a "$log"
members_dir="$workdir/members"
mkdir -p "$members_dir"

awk -F'\t' '{
  rep=$1; mem=$2;
  gsub(/\r/,"",rep); gsub(/\r/,"",mem);
  print mem >> (dir "/" rep ".ids")
}' dir="$members_dir" "$cluster_tsv"

for f in "$members_dir"/*.ids; do
  rep="$(basename "$f" .ids)"
  if ! grep -Fxq "$rep" "$f"; then
    echo "$rep" >> "$f"
  fi
done

process_one_cluster() {
  local rep="$1"
  local ids="$members_dir/$rep.ids"
  local cdir="$clustdir/$rep"
  mkdir -p "$cdir"

  local n
  n="$(wc -l < "$ids" | tr -d ' ')"

  local clufa="$cdir/cluster.fa"
  extract_fasta_by_ids "$in_fa" "$ids" "$clufa"

  if [[ ! -s "$clufa" ]]; then
    echo "[WARN] $rep: extracted FASTA empty; skipping" >&2
    return 0
  fi

  # Singleton => consensus is the only member
  if [[ "$n" -le 1 ]]; then
    local only_id
    only_id="$(head -n 1 "$ids")"
    local out="$cdir/consensus.fa"
    extract_fasta_by_ids "$in_fa" <(printf "%s\n" "$only_id") "$out"
    awk -v REP="$rep" -v N="$n" '
      /^>/{print ">" REP "|n=" N; next}
      {print}
    ' "$out" > "$cdir/consensus.renamed.fa"
    return 0
  fi

  # Align
  local aln="$cdir/cluster.aln.fa"
  run_mafft "$clufa" "$aln"

  # Optional trim
  local aln_used="$cdir/cluster.aln.used.fa"
  maybe_trimal "$aln" "$aln_used"

  # Generate consensus
  local cons_raw="$cdir/consensus.raw.fa"

  if [[ "$cons_tool" == "consambig" ]]; then
    # consambig writes to -outseq; -name sets sequence name
    consambig -sequence "$aln_used" -outseq "$cons_raw" -name "${rep}_consensus" >/dev/null 2>&1
  else
    # HMMER route
    local hmm="$cdir/cluster.hmm"
    hmmbuild --dna "$hmm" "$aln_used" >/dev/null 2>&1

    if [[ "$use_fancy" == "1" ]]; then
      if [[ -n "$fancy_minl" || -n "$fancy_minu" ]]; then
        hmmemit -C ${fancy_minl:+--minl "$fancy_minl"} ${fancy_minu:+--minu "$fancy_minu"} "$hmm" > "$cons_raw"
      else
        hmmemit -C "$hmm" > "$cons_raw"
      fi
    else
      hmmemit -c "$hmm" > "$cons_raw"
    fi
  fi

  # Normalize header
  awk -v REP="$rep" -v N="$n" '
    /^>/{print ">" REP "|n=" N; next}
    {print}
  ' "$cons_raw" > "$cdir/consensus.renamed.fa"
}

export -f extract_fasta_by_ids run_mafft maybe_trimal process_one_cluster
export in_fa members_dir clustdir have_seqkit mafft_mode mafft_extra
export cons_tool use_fancy fancy_minl fancy_minu do_trim trimal_args

echo "[INFO] Building per-cluster consensuses..." | tee -a "$log"

# Prevent MAFFT oversubscription when parallelizing clusters
export OMP_NUM_THREADS=1

ls "$members_dir"/*.ids | sed 's#.*/##; s/\.ids$//' \
  | xargs -I{} -P "$threads" bash -c 'process_one_cluster "$@"' _ {}

# -----------------------------
# 3) Concatenate final library
# -----------------------------
echo "[INFO] Writing final consensus library: $consensus_fa" | tee -a "$log"
: > "$consensus_fa"
find "$clustdir" -name 'consensus.renamed.fa' -type f -print0 \
  | sort -z \
  | xargs -0 cat >> "$consensus_fa"

n_cons="$(grep -c '^>' "$consensus_fa" || true)"
echo "[INFO] Done. Consensus sequences: $n_cons" | tee -a "$log"
echo "[INFO] Output: $consensus_fa" | tee -a "$log"
