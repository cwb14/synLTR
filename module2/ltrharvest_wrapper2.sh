#!/usr/bin/env bash
# nest_ltr_detector.sh
# Iterative nested LTR-RT detection wrapper for synLTR/module2 pipeline.
#
# Usage:
#   bash nest_ltr_detector.sh --genome genome.fa [--proteins prot.fa] [--terminate_count 100]
#       [--max-rounds N] [--script_path ./synLTR/module2/] [--threads 20] [--out_prefix ltrs]
#       [--wfa-align] [--ltrharvest5-args "KEY=VALUE ..."] [--ltrharvest5-args-from-round N "KEY=VALUE ..."]
#
# Notes:
# - Runs Round 1 on the original genome, then masks the ORIGINAL genome each round to build genome_r{N}.fa for next round.
# - Stops when: detected LTR-RTs < terminate_count, max-rounds reached, or ltrharvest5.py exits early
#   (e.g. no LTR-RT candidates found) -- early exits are handled gracefully.
# - After the loop, reconcile_nests.py pools all rounds, detects cross-round containment,
#   and writes depth-bucketed {OUT_PREFIX}_depth{N}_ltr.{tsv,fa} files. The raw per-round
#   {OUT_PREFIX}_r{N}_ltr.* files are preserved. depth{N} = elements with N layers of
#   LTR-RT nested inside (depth0 = no inner, i.e. "unnested"; depth1 = single-nested; etc.).
# - Uses IUPAC ambiguity codes for masking letters; excludes 'V' because far-character uses V.
# - At each next round:
#     scn-max-ret-len, -maxdistltr, -D  += 15000
#     scn-max-int-len                  += 14000
# - Builds temp_ltr_2pass_lib.fa by concatenating *all prior* libraries (r1..rK) and cleaning to A/T/C/G only.
# - Deletes temp_ltr_2pass_lib.fa at end.
# - Rounds 2+ automatically add --exclude-run-char V.
#
# Extra ltrharvest5.py args:
#   --ltrharvest5-args "KEY=VALUE [KEY2=VALUE2 ...]"
#       Applied to ALL rounds. Boolean flags: KEY=true (e.g. clean=true).
#       Can be specified multiple times.
#
#   --ltrharvest5-args-from-round N "KEY=VALUE [...]"
#       Applied starting from round N onward. For the same KEY, later from-round
#       values override earlier ones (last-wins per key, per round).
#       Can be specified multiple times.
#
#   Examples:
#     --ltrharvest5-args "clean=true verbose=true"
#     --ltrharvest5-args-from-round 2 "clean=true"
#     --ltrharvest5-args-from-round 1 "tesorter-rule=70-70-80"
#     --ltrharvest5-args-from-round 3 "tesorter-rule=70-40-80"

set -euo pipefail

# Print the exact command for reproducibility
echo "Command: $0 $*"
echo ""

WRAPPER_START=$SECONDS

# ----------------------------
# Defaults
# ----------------------------
GENOME=""
PROTEINS=""
TERMINATE_COUNT=100
MAX_ROUNDS_OVERRIDE=""   # empty = use IUPAC_SEQ length (10)
SCRIPT_PATH=""
THREADS=20
OUT_PREFIX="ltrs"
RUN_TRF=false
WFA_ALIGN=false

# Storage for extra ltrharvest5.py arg directives
# Each entry is tab-separated: "FROMROUND\tKEY\tVALUE\tIS_BOOL"
declare -a EXTRA_ARG_DIRECTIVES=()

# ----------------------------
# Helpers
# ----------------------------
die() { echo "ERROR: $*" >&2; exit 1; }

usage() {
  cat >&2 <<'USAGE_EOF'
Usage:
  bash nest_ltr_detector.sh --genome genome.fa [--proteins prot.fa] [--terminate_count 100]
      [--max-rounds N] [--script_path ./synLTR/module2/] [--threads 20] [--out_prefix ltrs]
      [--wfa-align] [--ltrharvest5-args "KEY=VALUE ..."] [--ltrharvest5-args-from-round N "KEY=VALUE ..."]

Required:
  --genome              Genome FASTA (.fa/.fasta)

Optional:
  --proteins            Protein FASTA for ltrharvest5.py
  --terminate_count     Stop if detected LTR-RTs in latest library < this count (default 100)
  --max-rounds          Maximum number of rounds to run (default: up to 10, limited by IUPAC codes).
                        Use 1 for non-nested only, 2 for single-level nesting, etc.
  --script_path         Path containing ltrharvest5.py and mask_ltr.py (default: same dir as this script)
  --threads             Threads for ltrharvest5.py (default 20)
  --out_prefix          Output prefix (default ltrs)
  --run-trf             Enable TRF instead of default --no-trf
  --wfa-align           Use WFA instead of mafft for Kmer2LTR pairwise alignment (~30-50x faster)

Extra ltrharvest5.py options:
  --ltrharvest5-args "KEY=VALUE [KEY2=VALUE2 ...]"
      Pass additional options to ltrharvest5.py for ALL rounds.
      KEY is the option name without leading '--'. Use '=' to separate key and value.
      Boolean flags (no value): KEY=true  (e.g. clean=true, verbose=true)
      Can be specified multiple times.

  --ltrharvest5-args-from-round N "KEY=VALUE [...]"
      Pass additional options to ltrharvest5.py starting from round N onward.
      For the same KEY, the directive with the highest applicable from-round wins.
      Can be specified multiple times.

  Examples:
      # Always clean workdirs and be verbose
      --ltrharvest5-args "clean=true verbose=true"

      # Clean only from round 2 onward
      --ltrharvest5-args-from-round 2 "clean=true"

      # Use a relaxed tesorter-rule for rounds 1-2, stricter from round 3
      --ltrharvest5-args-from-round 1 "tesorter-rule=70-70-80"
      --ltrharvest5-args-from-round 3 "tesorter-rule=70-40-80"

Notes:
  - Rounds 2+ automatically receive --exclude-run-char V (no need to specify manually).
  - If ltrharvest5.py exits with an error (e.g. no LTR-RT candidates / Kmer2LTR failure),
    the run stops gracefully rather than crashing the whole pipeline.
USAGE_EOF
}

abspath_dir() {
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
  grep -c '^>' "$fa" || true
}

clean_and_concat_libs() {
  local out="$1"; shift
  local tmpcat
  tmpcat="$(mktemp)"
  cat "$@" > "$tmpcat"
  awk '/^>/ {printf("\n%s\n",$0);next;} {printf("%s",$0);} END {printf("\n");}' "$tmpcat" \
    | sed '/^>/! s/[^ATCGatcg]//g' > "$out"
  rm -f "$tmpcat"
}

# Parse a space-separated "KEY=VALUE ..." string and append directives for from_round.
parse_kv_string() {
  local from_round="$1"
  local kv_str="$2"
  local pair key val is_bool

  for pair in $kv_str; do
    if [[ "$pair" != *"="* ]]; then
      die "Bad key-value entry '${pair}': must be KEY=VALUE or KEY=true (use --help)"
    fi
    key="${pair%%=*}"
    val="${pair#*=}"
    is_bool=0
    if [[ "$val" == "true" ]]; then
      is_bool=1
      val=""
    fi
    # Tab-separated: FROMROUND TAB KEY TAB VALUE TAB IS_BOOL
    EXTRA_ARG_DIRECTIVES+=( "$(printf '%s\t%s\t%s\t%s' "$from_round" "$key" "$val" "$is_bool")" )
  done
}

# Build the extra ltrharvest5.py args for a given round.
# Reads EXTRA_ARG_DIRECTIVES; writes result into the array named by $2 (nameref).
# For the same KEY, the directive with the highest from_round that is <= round wins.
build_extra_args_for_round() {
  local round="$1"
  local -n _out_arr="$2"
  _out_arr=()

  # For each key, track best (highest applicable from_round) val and is_bool
  declare -A best_from=()
  declare -A best_val=()
  declare -A best_bool=()
  declare -a key_order=()

  local from_round key val is_bool
  for directive in "${EXTRA_ARG_DIRECTIVES[@]}"; do
    from_round="$(printf '%s' "$directive" | cut -f1)"
    key="$(printf '%s' "$directive" | cut -f2)"
    val="$(printf '%s' "$directive" | cut -f3)"
    is_bool="$(printf '%s' "$directive" | cut -f4)"

    if [[ "$round" -ge "$from_round" ]]; then
      if [[ -z "${best_from[$key]+_}" ]]; then
        key_order+=("$key")
        best_from["$key"]="$from_round"
        best_val["$key"]="$val"
        best_bool["$key"]="$is_bool"
      elif [[ "$from_round" -ge "${best_from[$key]}" ]]; then
        best_from["$key"]="$from_round"
        best_val["$key"]="$val"
        best_bool["$key"]="$is_bool"
      fi
    fi
  done

  for key in "${key_order[@]}"; do
    if [[ "${best_bool[$key]}" == "1" ]]; then
      _out_arr+=( "--${key}" )
    else
      _out_arr+=( "--${key}" "${best_val[$key]}" )
    fi
  done
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
    --max-rounds) MAX_ROUNDS_OVERRIDE="${2:-}"; shift 2;;
    --script_path) SCRIPT_PATH="${2:-}"; shift 2;;
    --threads) THREADS="${2:-}"; shift 2;;
    --out_prefix) OUT_PREFIX="${2:-}"; shift 2;;
    --run-trf) RUN_TRF=true; shift;;
    --wfa-align) WFA_ALIGN=true; shift;;
    --ltrharvest5-args)
      parse_kv_string 1 "${2:-}"
      shift 2;;
    --ltrharvest5-args-from-round)
      _from="${2:-}"
      _kv="${3:-}"
      [[ "$_from" =~ ^[1-9][0-9]*$ ]] || die "--ltrharvest5-args-from-round requires a positive integer round number as the next argument"
      parse_kv_string "$_from" "$_kv"
      shift 3;;
    -h|--help) usage; exit 0;;
    *) die "Unknown argument: $1 (use --help)";;
  esac
done

[[ -z "$GENOME" ]] && die "--genome is required"
[[ -f "$GENOME" ]] || die "Genome not found: $GENOME"
if [[ -n "$PROTEINS" ]]; then
  [[ -f "$PROTEINS" ]] || die "Proteins not found: $PROTEINS"
fi

if [[ -z "$SCRIPT_PATH" ]]; then
  SCRIPT_PATH="$(abspath_dir "${BASH_SOURCE[0]}")"
else
  SCRIPT_PATH="${SCRIPT_PATH%/}"
fi

LTRHARVEST="${SCRIPT_PATH}/ltrharvest5.py"
MASKLTR="${SCRIPT_PATH}/mask_ltr.py"
RECONCILER="${SCRIPT_PATH}/reconcile_nests.py"
[[ -f "$LTRHARVEST" ]] || die "Missing: $LTRHARVEST"
[[ -f "$MASKLTR" ]] || die "Missing: $MASKLTR"
[[ -f "$RECONCILER" ]] || die "Missing: $RECONCILER"

# ----------------------------
# Config: IUPAC codes to use for successive rounds (exclude V)
# ----------------------------
IUPAC_SEQ=(N R D Y S W K M B H)   # 10 codes => up to 10 rounds
IUPAC_MAX="${#IUPAC_SEQ[@]}"       # 10

if [[ -n "$MAX_ROUNDS_OVERRIDE" ]]; then
  [[ "$MAX_ROUNDS_OVERRIDE" =~ ^[1-9][0-9]*$ ]] || die "--max-rounds must be a positive integer"
  if [[ "$MAX_ROUNDS_OVERRIDE" -gt "$IUPAC_MAX" ]]; then
    echo "WARNING: --max-rounds ${MAX_ROUNDS_OVERRIDE} exceeds IUPAC code limit (${IUPAC_MAX}); capping at ${IUPAC_MAX}." >&2
    MAX_ROUNDS="$IUPAC_MAX"
  else
    MAX_ROUNDS="$MAX_ROUNDS_OVERRIDE"
  fi
else
  MAX_ROUNDS="$IUPAC_MAX"
fi

echo "Max rounds set to: ${MAX_ROUNDS}"

# ----------------------------
# Static args
# ----------------------------
SCN_MIN_LTR_LEN=10
SCN_MIN_RET_LEN=80
SCN_MIN_INT_LEN=0
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

if [[ "$RUN_TRF" == true ]]; then
  TRF_OPTS=(--trf --trf-args "-a 5 -b 30 -g 30 -G 1 -s 150 -p 10" --trf-min-copy 40)
else
  TRF_OPTS=(--no-trf)
fi

WFA_OPTS=()
if [[ "$WFA_ALIGN" == true ]]; then
  WFA_OPTS=(--wfa-align)
fi

SIZE=500000
TESORTER_RULE="70-70-80"
TSD_PASS2="--tsd-pass2"
TESORTER_COV=20
TESORTER_EVAL="1e-2"
NESTED_FLANK_MIN=10
NESTED_BASE_MIN=800
FAR_CHARACTER="V"
MASK_DISTANCE=15000

base_scn_max_ret=150000
base_scn_max_int=140000
base_maxdistltr=15000
base_LTRF_D=15000

overlap_for_round() {
  local r="$1"
  echo $((25000 + (r-1)*15000))
}

orig_genome="$GENOME"
temp_lib="temp_ltr_2pass_lib.fa"

# ----------------------------
# Main loop
# ----------------------------
libs=()
# Tracks every round that produced a non-empty lib (used for the reconciler at the end)
completed_round_prefixes=()

for (( round=1; round<=MAX_ROUNDS; round++ )); do
  round_tag="r${round}"
  out_prefix_round="${OUT_PREFIX}_${round_tag}"

  if [[ "$round" -eq 1 ]]; then
    genome_in="$orig_genome"
  else
    genome_in="${orig_genome}_r$((round-1)).fa"
    [[ -f "$genome_in" ]] || die "Expected masked genome not found for round ${round}: $genome_in"
  fi

  scn_max_ret=$((base_scn_max_ret + (round-1)*15000))
  maxdistltr=$((base_maxdistltr + (round-1)*15000))
  ltrf_D=$((base_LTRF_D + (round-1)*15000))
  scn_max_int=$((base_scn_max_int * round))
  overlap="$(overlap_for_round "$round")"

  ltrharvest_args="-mindistltr ${LTR_MINDISTLTR} -minlenltr ${LTR_MINLENLTR} -maxlenltr ${LTR_MAXLENLTR} -mintsd ${LTR_MINTSD} -maxtsd ${LTR_MAXTSD} -similar ${LTR_SIMILAR} -vic ${LTR_VIC} -seed ${LTR_SEED} -seqids yes -xdrop ${LTR_XDROP} -maxdistltr ${maxdistltr}"
  ltrfinder_args="-w ${LTRF_W} ${LTRF_C} -D ${ltrf_D} -d ${LTRF_d} -L ${LTRF_L} -l ${LTRF_l} -p ${LTRF_p} -M ${LTRF_M} -S ${LTRF_S}"

  pass2_opts=()
  if [[ "$round" -ge 2 ]]; then
    req=""
    for ((i=0; i<=round-2; i++)); do
      req+="${IUPAC_SEQ[$i]},"
    done
    req="${req%,}"
    pass2_opts+=(
      --pass2-classified-fasta "$temp_lib"
      --require-run-chars "$req"
      --exclude-run-char "$FAR_CHARACTER"   # Edit (4): always exclude V from round 2+
    )
  fi

  protein_opts=()
  if [[ -n "$PROTEINS" ]]; then
    protein_opts+=( --proteins "$PROTEINS" )
  fi

  # Edit (3): Build per-round extra args
  declare -a extra_round_args=()
  build_extra_args_for_round "$round" extra_round_args

  echo "============================================================"
  echo "Round ${round} / ${MAX_ROUNDS} (${round_tag})"
  echo "  genome_in:          $genome_in"
  echo "  out_prefix:         $out_prefix_round"
  echo "  scn-max-ret-len:    $scn_max_ret"
  echo "  scn-max-int-len:    $scn_max_int"
  echo "  -maxdistltr:        $maxdistltr"
  echo "  ltrfinder -D:       $ltrf_D"
  echo "  overlap:            $overlap"
  if [[ "$round" -ge 2 ]]; then
    echo "  pass2 lib:          $temp_lib"
    echo "  require-run-chars:  $req"
    echo "  exclude-run-char:   $FAR_CHARACTER"
  fi
  if [[ "${#extra_round_args[@]}" -gt 0 ]]; then
    echo "  extra ltrharvest5:  ${extra_round_args[*]}"
  fi

  # Edit (2): Run ltrharvest5.py, catching non-zero exit gracefully
  ltrharvest_exit=0
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
    --tesorter-use-ret \
    --tesorter-rule "$TESORTER_RULE" \
    $TSD_PASS2 \
    --tesorter-cov "$TESORTER_COV" \
    --tesorter-eval "$TESORTER_EVAL" \
    --nested-flank-min "$NESTED_FLANK_MIN" \
    --nested-base-min "$NESTED_BASE_MIN" \
    "${pass2_opts[@]}" \
    "${WFA_OPTS[@]}" \
    "${extra_round_args[@]}" \
    || ltrharvest_exit=$?
  set +x

  # Edit (2): Graceful handling of early ltrharvest5.py failure
  if [[ "$ltrharvest_exit" -ne 0 ]]; then
    echo "" >&2
    echo "============================================================" >&2
    echo "WARNING: ltrharvest5.py exited with code ${ltrharvest_exit} on round ${round}." >&2
    echo "This typically means no LTR-RT candidates were found (e.g. Kmer2LTR" >&2
    echo "reported no usable data). Stopping gracefully after ${round} round(s)." >&2
    echo "============================================================" >&2
    break
  fi

  lib="${out_prefix_round}_ltr.fa"
  if [[ ! -s "$lib" ]]; then
    echo "" >&2
    echo "============================================================" >&2
    echo "WARNING: Expected library not found or empty: $lib" >&2
    echo "Stopping gracefully after round ${round}." >&2
    echo "============================================================" >&2
    break
  fi

  n_hits="$(count_fasta_headers "$lib")"
  echo "Round ${round}: detected ${n_hits} LTR-RTs in ${lib}"

  # Record this successful round for the post-loop reconciler
  completed_round_prefixes+=( "$out_prefix_round" )

  if [[ "$n_hits" -lt "$TERMINATE_COUNT" ]]; then
    echo "Terminate: ${n_hits} < ${TERMINATE_COUNT}. No further rounds will be run."
    break
  fi

  if [[ "$round" -eq "$MAX_ROUNDS" ]]; then
    echo "Reached MAX_ROUNDS=${MAX_ROUNDS}. Stopping."
    break
  fi

  # Prepare for next round: mask original genome
  next_feature_char="${IUPAC_SEQ[$((round-1))]}"
  masked_out="${orig_genome}_r${round}.fa"

  extra_mask_opts=()
  extra_features_tmp=""

  if (( ${#libs[@]} > 0 )); then
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

  libs+=( "$lib" )
  echo "Rebuilding ${temp_lib} from ${#libs[@]} libraries..."
  clean_and_concat_libs "$temp_lib" "${libs[@]}"
done

# ----------------------------
# Reconcile nested-status across rounds into depth-bucketed outputs.
# Each element's inward-chain depth (how many LTR-RT layers are nested inside
# it) determines which {OUT_PREFIX}_depth{N}_ltr.{tsv,fa} file it lands in.
# The raw per-round {OUT_PREFIX}_r{N}_ltr.* files are left untouched.
# ----------------------------
if (( ${#completed_round_prefixes[@]} > 0 )); then
  recon_tsv_args=()
  recon_fa_args=()
  recon_scn_args=()
  for prefix in "${completed_round_prefixes[@]}"; do
    recon_tsv_args+=( "${prefix}_ltr.tsv" )
    recon_fa_args+=(  "${prefix}_ltr.fa"  )
    recon_scn_args+=( "${prefix}.work/${prefix}.ltrtools.stitched.scn" )
  done

  echo ""
  echo "============================================================"
  echo "Reconciling ${#completed_round_prefixes[@]} round(s) into depth-bucketed outputs..."
  set -x
  python "$RECONCILER" \
    --out-prefix "$OUT_PREFIX" \
    --tsv "${recon_tsv_args[@]}" \
    --fa  "${recon_fa_args[@]}" \
    --scn "${recon_scn_args[@]}"
  set +x
fi

# Cleanup
rm -f "$temp_lib" 2>/dev/null || true
rm -rf ./tools/ 2>/dev/null || true

WRAPPER_ELAPSED=$((SECONDS - WRAPPER_START))
WRAPPER_H=$((WRAPPER_ELAPSED / 3600))
WRAPPER_M=$(( (WRAPPER_ELAPSED % 3600) / 60 ))
WRAPPER_S=$((WRAPPER_ELAPSED % 60))
if [[ "$WRAPPER_H" -gt 0 ]]; then
  WRAPPER_TIME="${WRAPPER_H}:$(printf '%02d' $WRAPPER_M):$(printf '%02d' $WRAPPER_S)"
else
  WRAPPER_TIME="${WRAPPER_M}:$(printf '%02d' $WRAPPER_S)"
fi
echo ""
echo "Wrapper total runtime: ${WRAPPER_TIME}"
echo "Done."
