#!/usr/bin/env bash
# nest_ltr_detector.sh
# Iterative nested LTR-RT detection wrapper for synLTR/module2 pipeline.
#
# Usage:
#   bash nest_ltr_detector.sh --genome genome.fa [--proteins prot.fa] [--terminate_count 100]
#       [--script_path ./synLTR/module2/] [--threads 20] [--out_prefix ltrs]
#
# Notes:
# - Runs Round 1 on the original genome, then masks the ORIGINAL genome each round to build genome_r{N}.fa for next round.
# - Stops after finishing a ltrharvest4.py round if detected LTR-RTs < terminate_count (default 100), or when max rounds reached.
# - Uses IUPAC ambiguity codes for masking letters; excludes 'V' because far-character uses V.
# - At each next round:
#     scn-max-ret-len, -maxdistltr, -D  += 15000
#     scn-max-int-len                  += 14000
# - Builds temp_ltr_2pass_lib.fa by concatenating *all prior* libraries (r1..rK) and cleaning to A/T/C/G only.
# - Deletes temp_ltr_2pass_lib.fa at end.

set -euo pipefail

# ----------------------------
# Defaults
# ----------------------------
GENOME=""
PROTEINS=""
TERMINATE_COUNT=100
SCRIPT_PATH=""
THREADS=20
OUT_PREFIX="ltrs"
RUN_TRF=false
# ----------------------------
# Helpers
# ----------------------------
die() { echo "ERROR: $*" >&2; exit 1; }

usage() {
  cat >&2 <<'EOF'
Usage:
  bash nest_ltr_detector.sh --genome genome.fa [--proteins prot.fa] [--terminate_count 100]
      [--script_path ./synLTR/module2/] [--threads 20] [--out_prefix ltrs]

Required:
  --genome              Genome FASTA (.fa/.fasta)

Optional:
  --proteins            Protein FASTA for ltrharvest4.py
  --terminate_count     Stop if detected LTR-RTs in latest library < this count (default 100)
  --script_path         Path containing ltrharvest4.py and mask_ltr.py (default: same dir as this script)
  --threads             Threads for ltrharvest4.py (default 20)
  --out_prefix          Output prefix (default ltrs)
  --run-trf            Enable TRF instead of default --no-trf
EOF
}

abspath_dir() {
  # portable-ish dirname->abs
  local d
  d="$(cd "$(dirname "$1")" && pwd)"
  echo "$d"
}

count_fasta_headers() {
  local fa="$1"
  if [[ ! -s "$fa" ]]; then
    echo 0
    return
  fi
  # count '>' at line start
  grep -c '^>' "$fa" || true
}

clean_and_concat_libs() {
  # args: output_file, lib1 lib2 ...
  local out="$1"; shift
  local tmpcat
  tmpcat="$(mktemp)"
  cat "$@" > "$tmpcat"
  # single-line sequences and strip non-ATCG from sequence lines
  awk '/^>/ {printf("\n%s\n",$0);next;} {printf("%s",$0);} END {printf("\n");}' "$tmpcat" \
    | sed '/^>/! s/[^ATCGatcg]//g' > "$out"
  rm -f "$tmpcat"
}

# ----------------------------
# Parse args
# ----------------------------
[[ $# -eq 0 ]] && { usage; exit 1; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --genome) GENOME="${2:-}"; shift 2;;
    --proteins) PROTEINS="${2:-}"; shift 2;;
    --terminate_count) TERMINATE_COUNT="${2:-}"; shift 2;;
    --script_path) SCRIPT_PATH="${2:-}"; shift 2;;
    --threads) THREADS="${2:-}"; shift 2;;
    --out_prefix) OUT_PREFIX="${2:-}"; shift 2;;
    --run-trf) RUN_TRF=true; shift;;
    -h|--help) usage; exit 0;;
    *) die "Unknown argument: $1 (use --help)";;
  esac
done

[[ -z "$GENOME" ]] && die "--genome is required"
[[ -f "$GENOME" ]] || die "Genome not found: $GENOME"
if [[ -n "$PROTEINS" ]]; then
  [[ -f "$PROTEINS" ]] || die "Proteins not found: $PROTEINS"
fi

# Resolve SCRIPT_PATH default: same directory as this wrapper
if [[ -z "$SCRIPT_PATH" ]]; then
  SCRIPT_PATH="$(abspath_dir "${BASH_SOURCE[0]}")"
else
  # allow trailing slash
  SCRIPT_PATH="${SCRIPT_PATH%/}"
fi

LTRHARVEST="${SCRIPT_PATH}/ltrharvest4.py"
MASKLTR="${SCRIPT_PATH}/mask_ltr.py"
[[ -f "$LTRHARVEST" ]] || die "Missing: $LTRHARVEST"
[[ -f "$MASKLTR" ]] || die "Missing: $MASKLTR"

# ----------------------------
# Config: IUPAC codes to use for successive rounds (exclude V)
# We start masking with N for round1->round2, then add R, D, Y, S, W, K, M, B, H ...
# N already included first; require-run-chars accumulates.
# ----------------------------
IUPAC_SEQ=(N R D Y S W K M B H)  # 10 codes => up to 10 rounds (r1..r10); nesting levels = rounds+1
MAX_ROUNDS="${#IUPAC_SEQ[@]}"     # 10

# ----------------------------
# Static args (per your pipeline)
# ----------------------------
SCN_MIN_LTR_LEN=100
SCN_MIN_RET_LEN=1000
SCN_MIN_INT_LEN=500
LTR_MAXLENLTR=7000
LTR_MINDISTLTR=100
LTR_MINLENLTR=100
LTR_VIC=60
LTR_SIMILAR=70
LTR_SEED=15
LTR_XDROP=10
LTR_MINTSD=0
LTR_MAXTSD=0

LTRF_W=2
LTRF_C="-C"
LTRF_p=20
LTRF_M="0.00"
LTRF_S="0.0"
LTRF_d=100
LTRF_l=100
LTRF_L=7000

# Chunking / TE sorter / nesting params
# TRF behavior (default = no TRF)
if [[ "$RUN_TRF" == true ]]; then
  TRF_OPTS=(--trf --trf-args "-a 5 -b 30 -g 30 -G 1 -s 150 -p 10" --trf-min-copy 40)
else
  TRF_OPTS=(--no-trf)
fi
SIZE=500000
TESORTER_RULE="75-80-80"
TSD_PASS2="--tsd-pass2"
TESORTER_COV=20
TESORTER_EVAL="1e-2"
NESTED_FLANK_MIN=10
NESTED_BASE_MIN=800
FAR_CHARACTER="V"
MASK_DISTANCE=15000

# ----------------------------
# Derived per-round increments
# Round 1 base values (from your commands)
# ----------------------------
base_scn_max_ret=15000
base_scn_max_int=14000
base_maxdistltr=15000
base_LTRF_D=15000

# Overlap schedule (from your examples; after that, use a simple increasing rule)
# r1: 25000, r2: 40000, r3: 55000, r4: 70000  => +15000 each round starting at r1
overlap_for_round() {
  local r="$1"
  echo $((25000 + (r-1)*15000))
}

# Genome masked file naming: genome_r1.fa used for round2 input, genome_r2.fa for round3, etc.
# We keep the ORIGINAL genome path for masking each time.
orig_genome="$GENOME"

temp_lib="temp_ltr_2pass_lib.fa"

# ----------------------------
# Main loop
# ----------------------------
libs=()   # collected library paths for concatenation
prev_masked_genome=""  # genome used as input to current round

for (( round=1; round<=MAX_ROUNDS; round++ )); do
  round_tag="r${round}"
  out_prefix_round="${OUT_PREFIX}_${round_tag}"

  # Genome input: round1 uses original genome; others use genome_r{round-1}.fa
  if [[ "$round" -eq 1 ]]; then
    genome_in="$orig_genome"
  else
    genome_in="${orig_genome}_r$((round-1)).fa"
    [[ -f "$genome_in" ]] || die "Expected masked genome not found for round ${round}: $genome_in"
  fi

  # Per-round increasing parameters
  scn_max_ret=$((base_scn_max_ret + (round-1)*15000))
  maxdistltr=$((base_maxdistltr + (round-1)*15000))
  ltrf_D=$((base_LTRF_D + (round-1)*15000))
  # round1 scn-max-int-len is 14000; round2 is 28000; etc (+14000 each round)
  scn_max_int=$((base_scn_max_int * round))

  overlap="$(overlap_for_round "$round")"

  # Build ltrharvest args strings
  ltrharvest_args="-mindistltr ${LTR_MINDISTLTR} -minlenltr ${LTR_MINLENLTR} -maxlenltr ${LTR_MAXLENLTR} -mintsd ${LTR_MINTSD} -maxtsd ${LTR_MAXTSD} -similar ${LTR_SIMILAR} -vic ${LTR_VIC} -seed ${LTR_SEED} -seqids yes -xdrop ${LTR_XDROP} -maxdistltr ${maxdistltr}"
  ltrfinder_args="-w ${LTRF_W} ${LTRF_C} -D ${ltrf_D} -d ${LTRF_d} -L ${LTRF_L} -l ${LTRF_l} -p ${LTRF_p} -M ${LTRF_M} -S ${LTRF_S}"

  # pass2 options for nested rounds
  pass2_opts=()
  if [[ "$round" -ge 2 ]]; then
    # require-run-chars: accumulated IUPAC letters used so far (up to round-1 masks, but we include current set too)
    # In your examples:
    #   round2 require N
    #   round3 require N,R
    #   round4 require N,R,D
    # So for round >=2, require chars = IUPAC_SEQ[0..round-2]
    req=""
    for ((i=0; i<=round-2; i++)); do
      req+="${IUPAC_SEQ[$i]},"
    done
    req="${req%,}"  # trim trailing comma

    pass2_opts+=( --pass2-classified-fasta "$temp_lib" --require-run-chars "$req" )
  fi

  # proteins option
  protein_opts=()
  if [[ -n "$PROTEINS" ]]; then
    protein_opts+=( --proteins "$PROTEINS" )
  fi

  echo "============================================================"
  echo "Round ${round} (${round_tag})"
  echo "  genome_in:          $genome_in"
  echo "  out_prefix:         $out_prefix_round"
  echo "  scn-max-ret-len:    $scn_max_ret"
  echo "  scn-max-int-len:    $scn_max_int"
  echo "  -maxdistltr:        $maxdistltr"
  echo "  ltrfinder -D:       $ltrf_D"
  echo "  overlap:            $overlap"
  if [[ "$round" -ge 2 ]]; then
    echo "  pass2 lib:          $temp_lib"
    echo "  require-run-chars:  ${pass2_opts[*]##*--require-run-chars }"
  fi

  # Run ltrharvest4.py (terminate after this step if needed)
  set -x
  python "$LTRHARVEST" \
    --genome "$genome_in" \
    "${protein_opts[@]}" \
    --threads "$THREADS" \
    --out-prefix "$out_prefix_round" \
    --scn-min-ltr-len "$SCN_MIN_LTR_LEN" \
    --scn-min-ret-len "$SCN_MIN_RET_LEN" \
    --scn-max-ret-len "$scn_max_ret" \
    --scn-min-int-len "$SCN_MIN_INT_LEN" \
    --scn-max-int-len "$scn_max_int" \
    --ltrharvest-args "$ltrharvest_args" \
    --ltrfinder-args "$ltrfinder_args" \
    "${TRF_OPTS[@]}" \
    --size "$SIZE" \
    --overlap "$overlap" \
    --tesorter-rule "$TESORTER_RULE" \
    $TSD_PASS2 \
    --tesorter-cov "$TESORTER_COV" \
    --tesorter-eval "$TESORTER_EVAL" \
    --nested-flank-min "$NESTED_FLANK_MIN" \
    --nested-base-min "$NESTED_BASE_MIN" \
    "${pass2_opts[@]}"
  set +x

  # Latest library produced this round
  lib="${out_prefix_round}.ltrharvest.full_length.dedup.fa.rexdb-plant.cls.lib.fa"
  if [[ ! -s "$lib" ]]; then
    die "Expected library not found or empty: $lib"
  fi

  n_hits="$(count_fasta_headers "$lib")"
  echo "Round ${round}: detected ${n_hits} LTR-RTs in ${lib}"

  # Decide whether to terminate (terminate AFTER ltrharvest4.py)
  if [[ "$n_hits" -lt "$TERMINATE_COUNT" ]]; then
    echo "Terminate: ${n_hits} < ${TERMINATE_COUNT}. No further rounds will be run."
    break
  fi

  # If next round would exceed max rounds, stop (can't add more IUPAC letters)
  if [[ "$round" -eq "$MAX_ROUNDS" ]]; then
    echo "Reached MAX_ROUNDS=${MAX_ROUNDS}. Stopping."
    break
  fi

  # Prepare for next round:
  # 1) Mask ORIGINAL genome using this round's features-fasta (the library) with next IUPAC letter as feature-character
  next_feature_char="${IUPAC_SEQ[$((round-1))]}"   # r1 masks with N, r2 masks with R, r3 masks with D, ...
  masked_out="${orig_genome}_r${round}.fa"

  echo "Masking original genome for next round: feature-character=${next_feature_char} -> ${masked_out}"

  # If we have previous libraries (r1..r{round-1}), feed them to --extra-features-fasta
  extra_mask_opts=()
  extra_features_tmp=""

  if (( ${#libs[@]} > 0 )); then
#    extra_features_tmp="$(mktemp --tmpdir "${OUT_PREFIX}.extra_features.XXXXXX.fa")"
    extra_features_tmp="$(mktemp "./${OUT_PREFIX}.extra_features.XXXXXX.fa")"
    cat "${libs[@]}" > "$extra_features_tmp"
    extra_mask_opts+=( --extra-features-fasta "$extra_features_tmp" )
  fi

  echo "Masking original genome for next round: feature-character=${next_feature_char} -> ${masked_out}"
  if (( ${#libs[@]} > 0 )); then
    echo "  extra-features-fasta: r1..r$((round-1)) (${#libs[@]} libs) => $extra_features_tmp"
  else
    echo "  extra-features-fasta: (none; round 1 special-case)"
  fi

  set -x
  python "$MASKLTR" \
    --features-fasta "$lib" \
    --genome "$orig_genome" \
    --feature-character "$next_feature_char" \
    --far-character "$FAR_CHARACTER" \
    --distance "$MASK_DISTANCE" \
    "${extra_mask_opts[@]}" > "$masked_out"
  set +x


  rm -f "$extra_features_tmp" 2>/dev/null || true
  
  [[ -s "$masked_out" ]] || die "Masked genome output is empty: $masked_out"

  # 2) Append library to list and rebuild temp_ltr_2pass_lib.fa for next round (concat ALL libs so far)
  libs+=( "$lib" )
  echo "Rebuilding ${temp_lib} from ${#libs[@]} libraries..."
  clean_and_concat_libs "$temp_lib" "${libs[@]}"
done

# Cleanup
rm -f "$temp_lib" 2>/dev/null || true
echo "Done."

