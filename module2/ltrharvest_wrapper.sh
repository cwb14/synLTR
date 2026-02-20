#!/usr/bin/env bash
set -euo pipefail

# nest_ltr_detector.sh
# Iterative LTR-RT discovery with increasing nesting levels via masking + pass2 library seeding.
#
# Usage:
#   bash nest_ltr_detector.sh --genome genome.fa [--proteins proteins.fa] [--terminate_count 100] \
#        [--script_path ./synLTR/module2/] [--threads 20] [--out_prefix ltrs]
#
# Notes:
# - Round 1 runs on the original genome (no require-run-chars).
# - Each subsequent round masks the original genome with an additional IUPAC ambiguity letter
#   (feature-character), then requires those run-chars in ltrharvest4.py.
# - Terminate after a round if the most recent library has < terminate_count headers.

usage() {
  cat <<'EOF'
bash nest_ltr_detector.sh --genome <genome.fa|fasta> [--proteins <proteins.fa|fasta>]
                         [--terminate_count 100]
                         [--script_path ./synLTR/module2/]
                         [--threads 20]
                         [--out_prefix ltrs]

Only --genome is required.
EOF
}

# -------------------------
# Defaults
# -------------------------
GENOME=""
PROTEINS=""
TERMINATE_COUNT=100
THREADS=20
OUT_PREFIX="ltrs"
SCRIPT_PATH=""   # if empty, infer from this script location

# -------------------------
# Parse args
# -------------------------
if [[ $# -eq 0 ]]; then usage; exit 1; fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    --genome)           GENOME="${2:-}"; shift 2;;
    --proteins)         PROTEINS="${2:-}"; shift 2;;
    --terminate_count)  TERMINATE_COUNT="${2:-}"; shift 2;;
    --threads)          THREADS="${2:-}"; shift 2;;
    --out_prefix)       OUT_PREFIX="${2:-}"; shift 2;;
    --script_path)      SCRIPT_PATH="${2:-}"; shift 2;;
    -h|--help)          usage; exit 0;;
    *) echo "ERROR: Unknown argument: $1" >&2; usage; exit 1;;
  esac
done

# -------------------------
# Validate / infer paths
# -------------------------
if [[ -z "$GENOME" ]]; then
  echo "ERROR: --genome is required." >&2
  usage
  exit 1
fi
if [[ ! -f "$GENOME" ]]; then
  echo "ERROR: Genome file not found: $GENOME" >&2
  exit 1
fi

if [[ -z "$SCRIPT_PATH" ]]; then
  # Default: directory where this wrapper lives
  SCRIPT_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
else
  # Normalize
  SCRIPT_PATH="$(cd "$SCRIPT_PATH" && pwd)"
fi

LTRHARVEST4="${SCRIPT_PATH}/ltrharvest4.py"
MASK_LTR="${SCRIPT_PATH}/mask_ltr.py"

if [[ ! -f "$LTRHARVEST4" ]]; then
  echo "ERROR: Missing dependency: $LTRHARVEST4" >&2
  exit 1
fi
if [[ ! -f "$MASK_LTR" ]]; then
  echo "ERROR: Missing dependency: $MASK_LTR" >&2
  exit 1
fi

if [[ -n "$PROTEINS" && ! -f "$PROTEINS" ]]; then
  echo "ERROR: Proteins file not found: $PROTEINS" >&2
  exit 1
fi

# -------------------------
# IUPAC letters for nesting
# Exclude V (reserved for far-character). N is included and used in round 2 masking.
# Order defines which new letter is introduced each successive round (after round 1).
# Up to 10 rounds total => 9 additional letters beyond round 1.
# -------------------------
IUPAC_LETTERS=(N R Y S W K M B D H)  # 10 letters, no V

MAX_ROUNDS=10

# -------------------------
# Param schedule
# -------------------------
SCN_MIN_LTR_LEN=100
SCN_MIN_RET_LEN=1000
SCN_MIN_INT_LEN=500

# Round 1 baseline (from your example)
BASE_SCN_MAX_RET_LEN=15000
BASE_SCN_MAX_INT_LEN=14000
BASE_MAXDISTLTR=15000
BASE_LTRFINDER_D=15000

# Increment each round AFTER round 1
INC_RET_AND_DIST=15000
INC_INT=14000

# Fixed args (from your example)
LTRHARVEST_ARGS_COMMON="-mindistltr 100 -minlenltr 100 -maxlenltr 7000 -mintsd 0 -maxtsd 0 -similar 70 -vic 60 -seed 15 -seqids yes -xdrop 10"
LTRFINDER_ARGS_COMMON="-w 2 -C -d 100 -L 7000 -l 100 -p 20 -M 0.00 -S 0.0"
SIZE=500000
OVERLAP=25000
TESORTER_RULE="75-80-80"
TESORTER_COV=20
TESORTER_EVAL="1e-2"
NESTED_FLANK_MIN=10
NESTED_BASE_MIN=800

# -------------------------
# Helpers
# -------------------------
count_fasta_headers() {
  local f="$1"
  if [[ ! -f "$f" ]]; then
    echo 0
    return
  fi
  # Count lines that start with '>'
  grep -c '^>' "$f" || true
}

clean_fasta_for_pass2() {
  local in_fa="$1"
  local out_fa="$2"
  # Rewrap to 1-line sequences, then strip non-ATCG from sequence lines
  cat "$in_fa" \
    | awk '/^>/ {printf("\n%s\n",$0); next;} {printf("%s",$0);} END {printf("\n");}' \
    | sed '/^>/! s/[^ATCGatcg]//g' \
    > "$out_fa"
}

run_ltrharvest4() {
  local genome_fa="$1"
  local out_prefix="$2"
  local scn_max_ret="$3"
  local scn_max_int="$4"
  local maxdistltr="$5"
  local ltrfinder_D="$6"
  local require_run_chars="${7:-}"     # optional
  local pass2_fasta="${8:-}"           # optional

  local cmd=(python "$LTRHARVEST4"
    --genome "$genome_fa"
    --threads "$THREADS"
    --out-prefix "$out_prefix"
    --scn-min-ltr-len "$SCN_MIN_LTR_LEN"
    --scn-min-ret-len "$SCN_MIN_RET_LEN"
    --scn-max-ret-len "$scn_max_ret"
    --scn-min-int-len "$SCN_MIN_INT_LEN"
    --scn-max-int-len "$scn_max_int"
    --ltrharvest-args "${LTRHARVEST_ARGS_COMMON} -maxdistltr ${maxdistltr}"
    --ltrfinder-args "${LTRFINDER_ARGS_COMMON} -D ${ltrfinder_D}"
    --no-trf
    --size "$SIZE"
    --overlap "$OVERLAP"
    --tesorter-rule "$TESORTER_RULE"
    --tsd-pass2
    --tesorter-cov "$TESORTER_COV"
    --tesorter-eval "$TESORTER_EVAL"
    --nested-flank-min "$NESTED_FLANK_MIN"
    --nested-base-min "$NESTED_BASE_MIN"
  )

  if [[ -n "$PROTEINS" ]]; then
    cmd+=(--proteins "$PROTEINS")
  fi
  if [[ -n "$pass2_fasta" ]]; then
    cmd+=(--pass2-classified-fasta "$pass2_fasta")
  fi
  if [[ -n "$require_run_chars" ]]; then
    cmd+=(--require-run-chars "$require_run_chars")
  fi

  echo ">>> Running: ${cmd[*]}" >&2
  "${cmd[@]}"
}

# Find most recent library for a given out_prefix pattern (e.g., ltrs, ltrs_r2, ...)
# We know the exact expected name, so just construct it.
library_path_for_prefix() {
  local op="$1"
  echo "${op}.ltrharvest.full_length.dedup.fa.rexdb-plant.cls.lib.fa"
}

# -------------------------
# Main loop
# -------------------------
TEMP_LIB="temp_ltr_2pass_lib.fa"

echo "=== nest_ltr_detector.sh ===" >&2
echo "Genome:          $GENOME" >&2
echo "Proteins:        ${PROTEINS:-<none>}" >&2
echo "Threads:         $THREADS" >&2
echo "Out prefix:      $OUT_PREFIX" >&2
echo "Script path:     $SCRIPT_PATH" >&2
echo "Terminate count: $TERMINATE_COUNT" >&2
echo "Max rounds:      $MAX_ROUNDS" >&2
echo >&2

# Track all library fasta files produced so far for concatenation into TEMP_LIB
declare -a LIBS_SO_FAR=()

# Round 1 (no masking, no pass2 lib, no require-run-chars)
round=1
scn_max_ret=$BASE_SCN_MAX_RET_LEN
scn_max_int=$BASE_SCN_MAX_INT_LEN
maxdistltr=$BASE_MAXDISTLTR
ltrfinder_D=$BASE_LTRFINDER_D

OUTP1="$OUT_PREFIX"
run_ltrharvest4 "$GENOME" "$OUTP1" "$scn_max_ret" "$scn_max_int" "$maxdistltr" "$ltrfinder_D"

LIB1="$(library_path_for_prefix "$OUTP1")"
n1="$(count_fasta_headers "$LIB1")"
echo ">>> Round 1 detected LTR-RTs: $n1 (in $LIB1)" >&2
LIBS_SO_FAR+=("$LIB1")

if (( n1 < TERMINATE_COUNT )); then
  echo ">>> Terminating after round 1 (count $n1 < $TERMINATE_COUNT)." >&2
  rm -f "$TEMP_LIB" || true
  exit 0
fi

# Subsequent rounds
# Round k (k>=2): mask original genome with feature-char = IUPAC_LETTERS[k-2]
# require-run-chars = comma-separated of all used letters so far
# pass2 lib = concat of all prior libraries (cleaned)
for (( round=2; round<=MAX_ROUNDS; round++ )); do
  letter="${IUPAC_LETTERS[$((round-2))]}"

  # Update scheduled params for this round
  scn_max_ret=$(( BASE_SCN_MAX_RET_LEN + (round-1)*INC_RET_AND_DIST ))
  scn_max_int=$(( BASE_SCN_MAX_INT_LEN + (round-1)*INC_INT ))
  maxdistltr=$(( BASE_MAXDISTLTR   + (round-1)*INC_RET_AND_DIST ))
  ltrfinder_D=$(( BASE_LTRFINDER_D + (round-1)*INC_RET_AND_DIST ))

  # Require-run-chars: all letters introduced up to this round (N for round2, N,R for round3, etc.)
  # As comma-separated list
  require_chars="$(IFS=,; echo "${IUPAC_LETTERS[*]:0:$((round-1))}")"

  # Mask genome (always mask the ORIGINAL genome, per your examples)
  masked_genome="${GENOME%.*}_r$((round-1)).fa"
  echo ">>> Masking for round $round using feature-character '${letter}' -> ${masked_genome}" >&2
  python "$MASK_LTR" \
    --features-fasta "${LIBS_SO_FAR[-1]}" \
    --genome "$GENOME" \
    --feature-character "$letter" \
    --far-character V \
    --distance 15000 \
    > "$masked_genome"

  # Build TEMP_LIB from all previous libraries
  echo ">>> Building pass2 library (${#LIBS_SO_FAR[@]} libs) -> $TEMP_LIB" >&2
  cat "${LIBS_SO_FAR[@]}" > "${TEMP_LIB}.raw.fa"
  clean_fasta_for_pass2 "${TEMP_LIB}.raw.fa" "$TEMP_LIB"
  rm -f "${TEMP_LIB}.raw.fa" || true

  # Run ltrharvest4 for this round
  OUTPk="${OUT_PREFIX}_r$round"
  run_ltrharvest4 "$masked_genome" "$OUTPk" "$scn_max_ret" "$scn_max_int" "$maxdistltr" "$ltrfinder_D" "$require_chars" "$TEMP_LIB"

  LIBk="$(library_path_for_prefix "$OUTPk")"
  nk="$(count_fasta_headers "$LIBk")"
  echo ">>> Round $round detected LTR-RTs: $nk (in $LIBk)" >&2
  LIBS_SO_FAR+=("$LIBk")

  if (( nk < TERMINATE_COUNT )); then
    echo ">>> Terminating after round $round (count $nk < $TERMINATE_COUNT)." >&2
    break
  fi
done

# Cleanup
rm -f "$TEMP_LIB" || true
echo "=== Done. ===" >&2
