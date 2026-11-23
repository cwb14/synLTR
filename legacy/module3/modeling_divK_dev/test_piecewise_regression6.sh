#!/usr/bin/env bash
set -euo pipefail

# =========================================
# Batched LTR-RT processing pipeline (no external skipping)
#
# This version **removes all skip/resume logic** from the wrapper.
# Kmer2LTR2.py now implements its own robust skipping/resume behavior,
# so we invoke every step unconditionally from this script.
#
# Still includes:
#   - Up to N pairs in parallel (Bash 5+ required for wait -n)
#   - Failure collection per background job
#
# Updated notes:
#   - Removed old Step 2 (header shortening) — extract_intact_LTR.py now emits short headers.
#   - Step "LTR lengths" uses the .key.tsv from extract_intact_LTR.py.
#   - Pass the resulting .key.tsv.len to Kmer2LTR2 via -D.
# =========================================

# Paths to required scripts (adjust if yours live elsewhere)
EXTRACT_LTR=./extract_intact_LTR.py
KMER2LTR=Kmer2LTR/Kmer2LTR2.py
DIVK_PL=./get_divK_Kmer2LTR.pl
SCATTER_R=./scatterplot_Kmer2LTR2.R

# Concurrency
# Set to how many BEDs you want processed simultaneously.
MAX_JOBS=1

# Quick sanity checks
need() {
  command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found in PATH"; exit 1; }
}
need python
need awk
need grep
need perl
need Rscript

for tool in "$EXTRACT_LTR" "$KMER2LTR" "$DIVK_PL" "$SCATTER_R"; do
  [[ -f "$tool" ]] || { echo "ERROR: Required script not found: $tool"; exit 1; }
done

# Where to record any failures from background jobs
FAIL_FILE="$(mktemp)"
cleanup() { rm -f "$FAIL_FILE"; }
trap cleanup EXIT

shopt -s nullglob

process_pair() {
  local bed="$1"

  local prefix="${bed%.bed}"           # e.g., gen200000_final
  local fasta="${prefix}.fasta"        # e.g., gen200000_final.fasta
  if [[ ! -f "$fasta" ]]; then
    echo "WARN: FASTA missing for $bed; skipping."
    return 0
  fi

  # Derive the base ID without the "_final" suffix for the mut-file
  local id_no_final="${prefix%_final}" # e.g., gen200000
  local mut_file="${id_no_final}_mut.txt"

  # Key paths
  local out_ltr_fa="${prefix}_LTR.fa"               # produced by extract_intact_LTR.py
  local key_tsv="${out_ltr_fa}.key.tsv"             # produced by extract_intact_LTR.py
  local key_len="${out_ltr_fa}.key.tsv.len"         # our derived LTR length file
  local k2l_base="${prefix}_LTR_Kmer2LTR"
  local divk_file="${k2l_base}.divK"
  local pdf_out="${k2l_base}.pdf"

  echo "=== Processing: $prefix ==="

  # 1) Extract LTR-RTs (emits ${out_ltr_fa} and ${key_tsv})
  echo "[1/6] Extracting LTR-RTs -> $out_ltr_fa (+ key TSV)"
  python "$EXTRACT_LTR" --bed "$bed" --genome "$fasta" --out_fasta "$out_ltr_fa"

  # 2) Get LTR length file from the key TSV
  #    key TSV format: <short_id>\t<header/extras...> (must contain LTRlen:NNN in second field)
  echo "[2/6] Generating LTR length file from key TSV -> $key_len"
  if [[ ! -f "$key_tsv" ]]; then
    echo "ERROR: Expected key TSV not found: $key_tsv"
    return 1
  fi
  awk -F'\t' '{if (match($2,/LTRlen:([0-9]+)/,m)) print $1"\t"m[1]}' "$key_tsv" > "$key_len"

  # 3) Calculate LTR-RT age with Kmer2LTR using .key.tsv.len
  #    Kmer2LTR2.py handles its own reuse/skip logic internally.
  echo "[3/6] Running Kmer2LTR2 -> $k2l_base"
  python "$KMER2LTR" \
    --max-win-overdisp 6 \
    --min-retained-fraction 0.6 \
    -i "$out_ltr_fa" \
    -D "$key_len" \
    -o "$k2l_base" \
    -p 200 \
    --reuse-existing \
    -t "/tmp/${prefix}_LTR_Kmer2LTR_temp"

  # 4) Calculate divergence (divK) — overwrite if present
  echo "[4/6] Computing divergence -> $divk_file"
  perl "$DIVK_PL" -b 5 -m 0 -w 0.001 "$k2l_base" > "$divk_file"

  # 5) Truth value for scatterplot (2 * last field from 'Accumulated Mutions / site:' line)
  local truth
  if [[ -f "$mut_file" ]]; then
    truth=$(awk '/Accumulated Mutions \/ site:/ {printf "%.5f\n", $NF*2}' "$mut_file")
    if [[ -z "${truth:-}" ]]; then
      echo "WARN: Pattern not found in $mut_file; defaulting truth to 0.01000"
      truth="0.01000"
    fi
  else
    echo "WARN: $mut_file not found; defaulting truth to 0.01000"
    truth="0.01000"
  fi
  echo "[5/6] Truth value (--truth): $truth"

  # 6) Build the plot PDF
  echo "[6/6] Building scatterplot -> $pdf_out"
  Rscript "$SCATTER_R" "$k2l_base" "$divk_file" no_cancel --truth="$truth" --out="$pdf_out"

  echo "=== Done: $prefix ==="
  echo
}

# Loop over all BED files that match the pattern, up to MAX_JOBS in parallel
for bed in gen*_final.bed; do
  # Launch a subshell for isolation and to trap errors per-job
  (
    set -euo pipefail
    # On any error in this subshell, record the prefix to FAIL_FILE
    # (Pre-resolve the prefix at trap-definition time so it's reliable)
    trap '
      pfx="'"${bed%.bed}"'"
      pfx="${pfx:-unknown}"
      echo "${pfx}" >> "'"$FAIL_FILE"'"
    ' ERR

    process_pair "$bed"
  ) &

  # Throttle to MAX_JOBS concurrent background jobs
  while (( $(jobs -rp | wc -l) >= MAX_JOBS )); do
    # wait for any job to finish (Bash 5+)
    wait -n || true
  done
done

# Wait for all remaining jobs
wait || true

# Report any failures
if [[ -s "$FAIL_FILE" ]]; then
  echo
  echo "Some jobs failed:"
  sort -u "$FAIL_FILE" | sed 's/^/  - /'
  exit 1
fi

echo "All done."
