#!/usr/bin/env bash
set -euo pipefail

# Bench wrapper for PrinTE get_divK_Kmer2LTR_benchmarking.pl
# Behavior:
# - Computes in-file UC/SC counts:
#     UC = lines not containing 'shared'
#     SC = lines containing 'shared'
# - max = round_to_step( 10 * min(UC,SC) )
# - low = round_to_step( max / 20 )
# - Randomly sample up to MAX_RUNS (uc,sc) pairs from the implied grid
#   WITHOUT enumerating the full grid
# - Run Perl with BOTH -uc and -sc
# - Parallelized with a simple job slot limiter
#
# Tunables via env:
#   JOBS=N       limit concurrency (default = #CPU cores)
#   MAX_RUNS=N   max random combos per file (default 10000)
#   STEP=N       step size for uc/sc ticks and rounding (default 10)

PERL_SCRIPT="./get_divK_Kmer2LTR_benchmarking.pl"
B_ARG=1
M_ARG=0
W_ARG=0.001

JOBS="${JOBS:-}"                 # export JOBS=N to limit concurrency; default = #CPU cores
MAX_RUNS="${MAX_RUNS:-10000}"    # max random combos PER FILE
STEP="${STEP:-10}"               # step size for grid ticks and rounding

shopt -s nullglob

# Determine default concurrency
if [[ -z "${JOBS}" ]]; then
  if command -v nproc >/dev/null 2>&1; then
    JOBS="$(nproc)"
  else
    JOBS="$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)"
  fi
fi

# Limit concurrent background jobs to $JOBS
wait_for_slot() {
  # shellcheck disable=SC2009
  while (( $(jobs -rp | wc -l) >= JOBS )); do
    # wait for any background job to finish; ignore exit status (we handle per-job)
    wait -n || true
  done
}

# Round to nearest multiple of STEP (standard half-up)
# Usage: round_to_step <integer>
round_to_step() {
  local n="$1"
  local half=$(( STEP / 2 ))
  echo $(( ((n + half) / STEP) * STEP ))
}

process_one() {
  local f="$1" uc="$2" sc="$3"

  local base="${f%_final_LTR_Kmer2LTR}"
  local mut="${base}_mut.txt"
  if [[ ! -f "$mut" ]]; then
    echo "WARNING: Missing mut file for '$f' (expected '$mut'); skipping." >&2
    return 0
  fi

  # Parse the 'Accumulated Mutions / site' value (times 2), tight float formatting
  local actual_value
  if ! actual_value="$(
    awk -F: '
      $1 ~ /^Accumulated Mutions \/ site/ {
        gsub(/[ \t]/,"",$2); val=$2 + 0; printf("%.10f", val*2); found=1
      }
      END { if (!found) exit 1 }
    ' "$mut" 2>/dev/null
  )"; then
    echo "ERROR: Could not find \"Accumulated Mutions / site\" in $mut" >&2
    return 0
  fi

  # Run the Perl script with BOTH -uc and -sc
  local out
  if ! out="$("$PERL_SCRIPT" -b "$B_ARG" -m "$M_ARG" -w "$W_ARG" -uc "$uc" -sc "$sc" "$f")"; then
    echo "ERROR: Perl failed for $f uc=$uc sc=$sc" >&2
    return 0
  fi

  # Extract rescaled S/U ratio from SUMMARY line
  local ratio
  if ! ratio="$(
    awk '
      /^SUMMARY/ && match($0, /rescaled_S\/U=([0-9.]+)/, m) { print m[1]; found=1 }
      END { if (!found) exit 1 }
    ' <<< "$out" 2>/dev/null
  )"; then
    echo "ERROR: Could not parse rescaled_S/U for $f uc=$uc sc=$sc" >&2
    return 0
  fi

  # Grab the first column from the last numeric row as the measured value
  local measured_value
  if ! measured_value="$(
    awk '
      /^[[:space:]]*[0-9]/ { last=$1 }
      END { if (last != "") print last; else exit 1 }
    ' <<< "$out" 2>/dev/null
  )"; then
    echo "ERROR: Could not parse last-row first column for $f uc=$uc sc=$sc" >&2
    return 0
  fi

  printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$measured_value" "$actual_value" "$ratio" "$f" "$uc" "$sc"
}

export -f process_one
export PERL_SCRIPT B_ARG M_ARG W_ARG

# Helpful on Ctrl-C: wait for children so we don't leave orphans
trap 'trap - INT; kill 0 2>/dev/null || true' INT

echo -e "measured_value\tactual_value\tratio\tfile\tuc\tsc"

# Expect files like gen*_final_LTR_Kmer2LTR with a paired *_mut.txt file
for f in gen*_final_LTR_Kmer2LTR; do
  base="${f%_final_LTR_Kmer2LTR}"
  mut="${base}_mut.txt"
  if [[ ! -f "$mut" ]]; then
    echo "WARNING: Missing mut file for '$f' (expected '$mut'); skipping." >&2
    continue
  fi

  # In-file counts
  # UC = lines NOT containing 'shared'; SC = lines containing 'shared'
  UC_COUNT=$(grep -v 'shared' "$f" | wc -l | awk '{print $1}')
  SC_COUNT=$(grep    'shared' "$f" | wc -l | awk '{print $1}')

  # Determine max (10x min of UC/SC), rounded to nearest STEP
  MIN_COUNT=$(( UC_COUNT < SC_COUNT ? UC_COUNT : SC_COUNT ))
  TEN_X=$(( MIN_COUNT * 10 ))
  MAX_RAW=$(round_to_step "$TEN_X")
  if (( MAX_RAW < STEP )); then
    MAX_RAW="$STEP"
  fi

  # Lower bound = max/20, rounded to nearest STEP; guard to >= STEP
  LOW_DIV=$(( MAX_RAW / 20 ))
  LOW_RAW=$(round_to_step "$LOW_DIV")
  if (( LOW_RAW < STEP )); then
    LOW_RAW="$STEP"
  fi

  # Compute ticks-per-axis and total combos for guardrails / logging
  TICKS=$(( (MAX_RAW - LOW_RAW) / STEP + 1 ))
  if (( TICKS <= 0 )); then
    echo "WARNING: Bad tick calculation for '$f' (low=$LOW_RAW, max=$MAX_RAW, step=$STEP); skipping." >&2
    continue
  fi
  TOTAL_COMBOS=$(( TICKS * TICKS ))

  # Build a random sample of unique (uc, sc) pairs without enumerating the full grid
  # Strategy:
  #  - If we need a small fraction of the grid: random unique pairs (O(target))
  #  - If we need a large fraction: systematic stride to cover evenly, then top up randomly
  mapfile -t COMBOS < <(
    awk -v low="$LOW_RAW" -v maxv="$MAX_RAW" -v step="$STEP" -v want="$MAX_RUNS" '
      function emit(u_idx, s_idx) {
        printf("%d\t%d\n", low + u_idx*step, low + s_idx*step)
      }
      BEGIN {
        srand()

        nx = int((maxv - low)/step) + 1
        if (nx < 1) exit 0

        total  = nx * nx
        target = (want < total ? want : total)

        # Heuristic: “large fraction” if > nx * 10 samples
        # This keeps runtime O(target) and avoids many RNG collisions.
        if (target > nx * 10) {
          su = int(nx / sqrt(target/nx)) + 1; if (su < 1) su = 1
          ss = su
          cnt = 0
          for (iu = 0; iu < nx && cnt < target; iu += su) {
            for (is = 0; is < nx && cnt < target; is += ss) {
              emit(iu, is); cnt++
            }
          }
          # Top up randomly if rounding made us short
          while (cnt < target) {
            iu = int(rand()*nx); is = int(rand()*nx)
            key = iu ":" is
            if (!(key in seen)) { seen[key]=1; emit(iu,is); cnt++ }
          }
        } else {
          # Random unique pairs
          cnt = 0
          while (cnt < target) {
            iu = int(rand()*nx); is = int(rand()*nx)
            key = iu ":" is
            if (!(key in seen)) { seen[key]=1; emit(iu,is); cnt++ }
          }
        }
      }'
  )

  if (( ${#COMBOS[@]} == 0 )); then
    echo "WARNING: No parameter combinations produced for '$f' (low=$LOW_RAW, max=$MAX_RAW, step=$STEP); skipping." >&2
    continue
  fi

  # Optional: warn if the implied grid is massive (just informational)
  if (( TOTAL_COMBOS > 50000000 )); then
    echo "NOTE: huge implied grid for '$f' (~${TOTAL_COMBOS} combos across ${TICKS}×${TICKS} ticks); sampling ${#COMBOS[@]} only." >&2
  fi

  # Launch sampled runs
  for line in "${COMBOS[@]}"; do
    uc="${line%%$'\t'*}"
    sc="${line#*$'\t'}"
    wait_for_slot
    process_one "$f" "$uc" "$sc" &
  done
done

# Wait for all background jobs
wait
