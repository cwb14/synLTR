#!/usr/bin/env bash
# nest_ltr_detector.sh
# Iterative nested LTR-RT detection wrapper for synLTR/module2 pipeline.
#
# Usage:
#   bash nest_ltr_detector.sh --genome genome.fa [genome2.fa ...] [--proteins prot.fa] [--terminate_count 100]
#       [--max-rounds N] [--script_path ./synLTR/module2/] [--threads 20] [--out_prefix PREFIX]
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

# FALSE-POSITIVE (FP) FAMILY CORRECTION (post-detection):
#   High-abundance non-LTR repeats (SINEs, DNA/low-complexity) can seed spurious
#   LTR-RT calls. After reconciliation, this wrapper:
#     (1) concatenates the depth-bucketed FASTAs and strips non-ACGT characters
#         (the reconciler's IUPAC-masked nested inners), leaving only the
#         outermost putative LTR-RT of each element -> {OUT_PREFIX}_all_ltr.fa
#     (2) clusters those with Kmer2LTR (consensus + internal clustering, id 0.75)
#     (3) runs flag_fp_families.py to flag FP families and purge their rows from
#         each depth TSV -> {OUT_PREFIX}_depth{N}_clean_ltr.tsv, then writes the
#         matching FP-purged FASTA -> {OUT_PREFIX}_depth{N}_clean_ltr.fa
#   If the FP fraction exceeds --fp-mask-threshold, flag_fp_families.py hard-masks
#   those repeats in the genome ({OUT_PREFIX}_FP_masked.fa) and the ENTIRE pipeline
#   is automatically re-run on the masked genome. Each attempt runs fully isolated
#   in its own staging directory, so first-run files can never be mistaken for the
#   final outputs; the hidden --dev-keep-fp-rounds flag retains every attempt for
#   debugging. See the "FP-correction orchestrator" section below.

# ANNOTATION ADD-ONS (post FP-correction):
#   (a) ltr_annotate.py inserts `strand` and `family` columns into every
#       {OUT_PREFIX}_depth{N}_ltr.tsv and _depth{N}_clean_ltr.tsv, between `tsd`
#       and `domains`. Strand is TEsorter2's call, falling back to protein-domain
#       order and then to the minimap2 pass-2 homology target; elements with no
#       usable evidence keep '.'. Family is the Kmer2LTR/mmseqs consensus
#       cluster, labelled {OUT_PREFIX}_fam00001 onward.
#   (b) ltr_tsv_to_gff3.py pools those tables into
#       {OUT_PREFIX}_all_depth_LTR_cleaned.gff3, and - when round 1 produced a
#       miniprot genic GFF - {OUT_PREFIX}_all_depth_protein_LTR_cleaned.gff3,
#       which additionally carries every miniprot protein alignment. Each
#       element carries its family's clade makeup as family_clades=.

# MULTIPLE GENOMES (shared family nomenclature):
#   --genome accepts several FASTAs. Closely-related species then share one
#   family vocabulary, which is what makes between-species family comparisons
#   meaningful. The pipeline splits into three phases:
#     (1) per genome, sequentially: detection rounds + reconcile_nests.py,
#         producing {P_g}_depth{N}_ltr.{tsv,fa} for each genome. This is the
#         detect-only worker, re-exec'd once per genome.
#     (2) pooled, once: every genome's depth FASTAs are concatenated into
#         {RUN}_all_ltr.fa, clustered by a single Kmer2LTR pass, and
#         FP-corrected by flag_fp_families.py -- so family membership and
#         false-positive calls are computed over the union, not per species.
#     (3) split back out, per genome: ltr_annotate.py and ltr_tsv_to_gff3.py
#         run once per genome against that one pooled cluster table, so
#         {RUN}_fam00001 denotes the same family in every genome's outputs.
#   Naming: each genome keeps its own <basename>_LTRs prefix; --out_prefix
#   names the shared pool (default 'merged'). With a single genome the two are
#   the same thing, so nothing about single-genome runs changes.
#   Sequence IDs must be unique across genomes -- pooled clustering keys
#   elements on 'chrom:start-end', so a shared seqid would cross-assign
#   families and cross-purge real elements. Checked before any work starts.

set -euo pipefail

# Print the exact command for reproducibility
echo "Command: $0 $*"
echo ""

WRAPPER_START=$SECONDS

# Capture the original argv up-front (before parsing consumes it) so the
# FP-correction orchestrator can re-invoke this script as a worker on a masked
# genome without having to reconstruct every flag.
ORIG_ARGS=( "$@" )

# ----------------------------
# Defaults
# ----------------------------
declare -a GENOMES=()
GENOME=""          # GENOMES[0]; the worker only ever sees one genome
N_GENOMES=0
RUN_PREFIX=""      # names the shared family namespace + merged artifacts
PROTEINS=""
TERMINATE_COUNT=100
MAX_ROUNDS_OVERRIDE=""   # empty = use IUPAC_SEQ length (10)
SCRIPT_PATH=""
THREADS=20
OUT_PREFIX=""
RUN_TRF=true
WFA_ALIGN=false
KEEP_WEAK_HMM_PASS2=false
PASS2_ALIGNER="minimap2"
RUN_PLOTS=true

# FP-family correction (post-detection Kmer2LTR stage)
FP_MASK_THRESHOLD="0.10"   # user-tunable: FP-element fraction above which the
                           # genome is hard-masked and the pipeline auto-re-runs
DEV_KEEP_FP_ROUNDS=false   # hidden dev flag: keep every FP-round staging dir
MAX_FP_ROUNDS=10           # hidden safety cap on automatic FP-masking re-runs
                           # (overridable via --dev-max-fp-rounds)

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
  bash nest_ltr_detector.sh --genome genome.fa [genome2.fa ...] [--proteins prot.fa]
      [--terminate_count 100]
      [--max-rounds N] [--script_path ./synLTR/module2/] [--threads 20] [--out_prefix PREFIX]
      [--wfa-align] [--fp-mask-threshold 0.10]
      [--ltrharvest5-args "KEY=VALUE ..."] [--ltrharvest5-args-from-round N "KEY=VALUE ..."]

Required:
  --genome              Genome FASTA (.fa/.fasta, plain or .gz). One or more,
                        space-separated. With 2+ genomes every genome is
                        detected separately, then all detected elements are
                        pooled for a single clustering pass, so the `family`
                        column means the same thing in every genome's output.
                        Sequence IDs must be unique across genomes (checked
                        up front, before any work starts).

Optional:
  --proteins            Protein FASTA for ltrharvest5.py. One file, shared by
                        every genome.
  --terminate_count     Stop if detected LTR-RTs in latest library < this count (default 100)
  --max-rounds          Maximum number of rounds to run (default: up to 10, limited by IUPAC codes).
                        Use 1 for non-nested only, 2 for single-level nesting, etc.
  --script_path         Path containing ltrharvest5.py and mask_ltr.py (default: same dir as this script)
  --threads             Threads for ltrharvest5.py (default 20)
  --out_prefix          With ONE genome: output prefix (default
                        <genome_prefix>_LTRs), which also names the families.
                        With 2+ genomes: names the shared family namespace and
                        the merged cluster tables (default 'merged'); each
                        genome always keeps its own <basename>_LTRs prefix.
  --run-trf             Run TRF tandem-repeat masking (default: on)
  --no-trf              Disable TRF tandem-repeat masking
  --wfa-align           Use WFA instead of mafft for Kmer2LTR pairwise alignment (~30-50x faster)
  --keep-weak-hmm-pass2-matches
                        Include in cls.lib.fa the pass-2 matches whose target is a weak-HMM
                        candidate (had HMM hits but no clade -> LTR/unknown/unknown via augment).
                        By default these matches are dropped. Default: off.
  --pass2-aligner       TEBinSorter pass-2 aligner: minimap2 (default) or blast.
  --fp-mask-threshold   False-positive element fraction above which the genome is
                        hard-masked and the whole pipeline is automatically re-run
                        on the masked genome (default 0.10; higher tolerates more
                        false positives before masking/re-running).
  --no-plots                Skip the post-completion plotting stage
                            (ltrharvest_plots.sh: structure PDFs, summary PDF,
                            TEGV HTML into <out_prefix>_plots/). Default: run it.

Post-detection FP-family correction runs automatically: it clusters the detected
LTR-RTs with Kmer2LTR, purges false-positive families into
<out_prefix>_depth{N}_clean_ltr.{tsv,fa}, and (when FPs are pervasive) masks the
genome and re-runs. See the script header for details.

Annotation add-ons then run automatically: every depth table gains `strand` and
`family` columns, and the tables are pooled into
<out_prefix>_all_depth_LTR_cleaned.gff3 (plus
<out_prefix>_all_depth_protein_LTR_cleaned.gff3 when --proteins was given, which
also carries the miniprot alignments).

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

# Genome prefix = basename minus a trailing .gz and one FASTA extension.
# e.g. Athal_tair10_chr2.fa.gz -> Athal_tair10_chr2
genome_prefix_of() {
  local b
  b="$(basename "$1")"
  b="${b%.gz}"
  b="${b%.*}"
  [[ -n "$b" ]] || die "Could not derive a genome prefix from: $1"
  printf '%s' "$b"
}

# Two genomes whose basenames collapse to the same prefix would overwrite each
# other's outputs in the shared staging directory.
check_unique_basenames() {
  local -A seen=()
  local i p
  for i in "${!OUT_PREFIXES[@]}"; do
    p="${OUT_PREFIXES[$i]}"
    if [[ -n "${seen[$p]+_}" ]]; then
      die "Two genomes derive the same output prefix '${p}':
    ${GENOMES[${seen[$p]}]}
    ${GENOMES[$i]}
  Rename one input file so the basenames differ."
    fi
    seen["$p"]=$i
  done
}

# A seqid shared by two genomes is a correctness hazard, not a cosmetic one:
# flag_fp_families.py purges by bare 'chrom:start-end', so one species' false
# positives would delete the other species' real elements, and ltr_annotate.py's
# by_coord fallback would silently drop the ambiguous coordinates.
check_unique_seqids() {
  local pairs dupes id i
  pairs="$(mktemp "./${RUN_PREFIX}.seqids.XXXXXX")"
  echo "Checking seqid uniqueness across ${N_GENOMES} genomes..."
  for i in "${!GENOMES[@]}"; do
    gzip -cdf "${GENOMES[$i]}" \
      | awk -v tag="${OUT_PREFIXES[$i]}" \
            '/^>/ { split(substr($0,2), a, /[ \t]/); print a[1] "\t" tag }'
  done | LC_ALL=C sort -u > "$pairs"

  dupes="$(cut -f1 "$pairs" | LC_ALL=C uniq -d)"
  if [[ -n "$dupes" ]]; then
    echo "" >&2
    echo "ERROR: these sequence IDs appear in more than one genome:" >&2
    while IFS= read -r id; do
      [[ -z "$id" ]] && continue
      printf '  %-30s %s\n' "$id" \
        "$(awk -F'\t' -v k="$id" '$1==k {printf "%s ", $2}' "$pairs")" >&2
    done <<< "$(printf '%s\n' "$dupes" | head -20)"
    (( $(printf '%s\n' "$dupes" | wc -l) > 20 )) && echo "  ... (truncated)" >&2
    cat >&2 <<'FIXEOF'

Pooled clustering keys elements on 'chrom:start-end', so shared IDs would
cross-assign families AND cross-purge real elements between genomes.

Rename one genome's sequences before re-running, e.g.:
  seqkit replace -p '^' -r 'Aly_' in.fa > out.fa
  awk '/^>/{sub(/^>/,">Aly_")}1' in.fa > out.fa
FIXEOF
    rm -f "$pairs"
    exit 1
  fi
  rm -f "$pairs"
  echo "  seqids unique across all genomes."
}

# Resolve a file or directory to an absolute path (no realpath dependency).
abspath() {
  local p="$1"
  if [[ -d "$p" ]]; then
    (cd "$p" && pwd)
  else
    echo "$(cd "$(dirname "$p")" && pwd)/$(basename "$p")"
  fi
}

# ----------------------------
# FP-correction orchestrator helpers
# ----------------------------
# Ensure the Kmer2LTR clone used by the FP stage is available. ltrharvest5.py
# normally clones it into TOOLS_DIR during the rounds; clone defensively if the
# rounds produced nothing or the layout ever changes.
ensure_kmer2ltr_dir() {
  local k2l="${TOOLS_DIR}/Kmer2LTR"
  if [[ ! -f "${k2l}/Kmer2LTR.py" || ! -f "${k2l}/flag_fp_families.py" ]]; then
    echo "Kmer2LTR not present in ${k2l}; cloning..." >&2
    mkdir -p "$TOOLS_DIR"
    git clone --depth 1 https://github.com/cwb14/Kmer2LTR.git "$k2l" >&2 \
      || die "Failed to clone Kmer2LTR into ${k2l}"
  fi
  [[ -f "${k2l}/Kmer2LTR.py" && -f "${k2l}/flag_fp_families.py" ]] \
    || die "Kmer2LTR clone incomplete in ${k2l}"
}

# Move (default) or copy (--dev-keep-fp-rounds) the final attempt's outputs into
# the user's working directory, dropping the isolation scaffolding and the
# FP-masked trigger file (never a final output). A same-named output from a prior
# run in the destination is removed first so directory moves overwrite cleanly.
promote_fp_outputs() {
  local src="$1" dst="$2" e name p
  for p in "${OUT_PREFIXES[@]}"; do
    rm -f "${src}/${p}.input_genome.fa"* 2>/dev/null || true
    rm -f "${src}/${p}_FP_masked.fa" 2>/dev/null || true
  done

  shopt -s nullglob
  local entries=( "${src}"/* )
  shopt -u nullglob
  for e in "${entries[@]}"; do
    name="$(basename "$e")"
    rm -rf "${dst:?}/${name}"
    if [[ "$DEV_KEEP_FP_ROUNDS" == true ]]; then
      cp -a "$e" "${dst}/"
    else
      mv "$e" "${dst}/"
    fi
  done
}

# ----------------------------
# Merged stage: pooling, clustering, FP correction, annotation, GFF3
# ----------------------------
# Collect a genome's depth tables, excluding the FP-cleaned variants that the
# merged stage itself produces ({p}_depth*_ltr.tsv would otherwise also match
# {p}_depth0_clean_ltr.tsv). Today that only stays safe because every attempt
# starts in a fresh directory; carry-forward breaks that assumption.
depth_tables_for() {
  local prefix="$1" kind="$2"
  local -n _dt_out="$3"
  local f
  _dt_out=()
  shopt -s nullglob
  for f in "${prefix}"_depth*_ltr."${kind}"; do
    [[ "$f" == *_clean_ltr."${kind}" ]] && continue
    _dt_out+=( "$f" )
  done
  shopt -u nullglob
}

# The merged stage needs Kmer2LTR. Reuse a per-genome clone if the detection
# workers left a complete one behind, else clone into a run-level dir.
resolve_merged_tools_dir() {
  local p
  for p in "${OUT_PREFIXES[@]}"; do
    if [[ -f "./${p}_tools/Kmer2LTR/Kmer2LTR.py" \
       && -f "./${p}_tools/Kmer2LTR/flag_fp_families.py" ]]; then
      TOOLS_DIR="./${p}_tools"
      return
    fi
  done
  TOOLS_DIR="./${RUN_PREFIX}_tools"
}

# One pooled flag_fp_families.py call: --domains-tsv takes nargs="+", so a
# single call writes every genome's _clean_ltr.tsv. Stage C (masking) is per
# genome, but its gate -- the FP fraction -- is computed from the cluster tables
# alone (flag_fp_families.py: frac = fp_fraction(len(fp_members), len(member2rep)))
# and is therefore identical for every genome. So parse the fraction from this
# first call and only make the remaining per-genome masking calls when it
# actually exceeded the threshold; below the threshold flag_fp_families.py
# returns before the expensive mmseqs step anyway.
run_fp_stage() {
  local cons="$1" int="$2" cons_fa="$3" flag_fp="$4"; shift 4
  local -a dtsvs=( "$@" )
  local first="${OUT_PREFIXES[0]}"
  local log="${first}_fpcheck.log"
  local i p frac

  # stderr is teed to $log via an fd swap, NOT process substitution: bash does
  # not wait for >(...) to finish, so the log would still be empty when the
  # fraction is parsed out of it a moment later. A pipeline is waited on.
  set -x
  { python "$flag_fp" \
      --consensus-cluster "$cons" \
      --internal-cluster "$int" \
      --ltr-fasta "$cons_fa" \
      --domains-tsv "${dtsvs[@]}" \
      -o "${RUN_PREFIX}_fpcheck" \
      --genome "${first}.input_genome.fa" \
      --masked-out "${first}_FP_masked.fa" \
      --threads "$THREADS" \
      --fp-mask-threshold "$FP_MASK_THRESHOLD" \
      2>&1 1>&3 | tee "$log" >&2; } 3>&1
  set +x

  frac="$(sed -n 's/^\[INFO\] FP fraction: [0-9]*\/[0-9]* = \([0-9.]*\).*/\1/p' "$log" | tail -1)"
  if [[ -z "$frac" ]] || ! awk -v f="$frac" -v t="$FP_MASK_THRESHOLD" \
        'BEGIN{exit !(f>t)}'; then
    echo "FP fraction within --fp-mask-threshold (${FP_MASK_THRESHOLD}); no masking."
    return 0
  fi

  echo "FP fraction ${frac} exceeded ${FP_MASK_THRESHOLD}; masking each genome."
  for (( i=1; i<N_GENOMES; i++ )); do
    p="${OUT_PREFIXES[$i]}"
    set -x
    { python "$flag_fp" \
        --consensus-cluster "$cons" \
        --internal-cluster "$int" \
        --ltr-fasta "$cons_fa" \
        -o "${RUN_PREFIX}_fpcheck" \
        --no-plot \
        --genome "${p}.input_genome.fa" \
        --masked-out "${p}_FP_masked.fa" \
        --threads "$THREADS" \
        --fp-mask-threshold "$FP_MASK_THRESHOLD" \
        2>&1 1>&3 | tee "${p}_fpcheck.log" >&2; } 3>&1
    set +x
  done
}

# Keep only records whose full header equals a surviving element id (cleaned TSV
# column 1; the reconciler writes depth-FASTA headers as exactly that id).
write_clean_fastas() {
  local dtsv clean_tsv depth_fa clean_fa
  for dtsv in "$@"; do
    clean_tsv="${dtsv%_ltr.tsv}_clean_ltr.tsv"
    depth_fa="${dtsv%_ltr.tsv}_ltr.fa"
    clean_fa="${dtsv%_ltr.tsv}_clean_ltr.fa"
    if [[ ! -s "$clean_tsv" || ! -s "$depth_fa" ]]; then
      echo "WARNING: skipping cleaned FASTA for ${dtsv}: missing ${clean_tsv} or ${depth_fa}" >&2
      continue
    fi
    awk -F'\t' 'NR==FNR { if ($0 !~ /^#/ && NF >= 2) keep[$1] = 1; next }
                /^>/ { p = (substr($0, 2) in keep) } p' \
      "$clean_tsv" "$depth_fa" > "$clean_fa"
    echo "Cleaned FASTA: ${clean_fa} ($(count_fasta_headers "$clean_fa") records)"
  done
}

# Every genome reads the one pooled cluster table and labels families from
# RUN_PREFIX, so the family column is directly comparable across species.
# --prefix still drives everything per-genome: TEsorter tables, PAFs,
# {prefix}_LTRRT_00001 ids, miniprot GFF discovery, output GFF3 names.
run_annotation_stage() {
  local p
  local -a fam_opts=() annot_tsvs=()
  # Resolve the pooled cluster table from the filesystem rather than a variable:
  # the orchestrator calls this in a different subshell from run_merged_stage,
  # so anything set there is already gone. Missing table -> no override, and
  # both scripts fall back to their own per-prefix glob.
  shopt -s nullglob
  local -a merged_cons=( "${RUN_PREFIX}"_all_ltr.consensus_id*_cluster.tsv )
  shopt -u nullglob
  if (( ${#merged_cons[@]} == 1 )); then
    fam_opts=( --consensus-cluster "${merged_cons[0]}"
               --family-prefix "$RUN_PREFIX" )
  elif (( ${#merged_cons[@]} > 1 )); then
    die "Expected exactly 1 pooled consensus cluster TSV, found ${#merged_cons[@]}"
  fi

  for p in "${OUT_PREFIXES[@]}"; do
    depth_tables_for "$p" tsv annot_tsvs
    if (( ${#annot_tsvs[@]} == 0 )); then
      echo "WARNING: no ${p}_depth<N>_ltr.tsv present; skipping the" >&2
      echo "strand/family annotation and GFF3 stages for ${p}." >&2
      continue
    fi
    echo ""
    echo "============================================================"
    echo "Annotating ${#annot_tsvs[@]} depth table(s) for ${p} with strand + family..."
    set -x
    python "$ANNOTATOR" --prefix "$p" --indir . "${fam_opts[@]}"
    set +x

    echo ""
    echo "============================================================"
    echo "Writing LTR-RT GFF3 for ${p}..."
    set -x
    python "$GFF3_WRITER" --prefix "$p" --indir . \
      --genome "${p}.input_genome.fa" "${fam_opts[@]}"
    set +x
  done
}

# Pooled post-detection stage: concatenate every genome's depth FASTAs, cluster
# them once, and FP-correct. Runs with cwd = the attempt's staging dir.
# Annotation is NOT done here -- an attempt is only known to be final after the
# FP logs it writes have been read, so the orchestrator calls
# run_annotation_stage separately once it knows.
run_merged_stage() {
  local all_ltr_fa="${RUN_PREFIX}_all_ltr.fa"
  local p cons_fa
  local -a depth_fas=() depth_tsvs=() dfas=() dtsvs=()

  for p in "${OUT_PREFIXES[@]}"; do
    depth_tables_for "$p" fa  dfas
    depth_tables_for "$p" tsv dtsvs
    depth_fas+=( "${dfas[@]}" )
    depth_tsvs+=( "${dtsvs[@]}" )
  done

  if (( ${#depth_fas[@]} == 0 )); then
    echo "WARNING: no depth-bucketed FASTA present; skipping the merged stage." >&2
    return 0
  fi

  echo ""
  echo "============================================================"
  echo "Merged stage: ${#depth_fas[@]} depth FASTA(s) from ${N_GENOMES} genome(s)"

  # (1) Pool, stripping non-ACGT -- the reconciler's IUPAC-masked nested inners
  #     -- so only the outermost putative LTR-RT of each element remains.
  cat "${depth_fas[@]}" \
    | awk '/^>/{printf "\n%s\n",$0;next}{gsub(/[^ACGTacgt]/,"");printf "%s",$0}END{print ""}' \
    > "$all_ltr_fa"

  if [[ ! -s "$all_ltr_fa" ]]; then
    # Degenerate (near-impossible) case: nothing to cluster. The primary depth
    # detection already succeeded, so warn and skip the optional FP add-on rather
    # than aborting and discarding good results.
    echo "WARNING: ${all_ltr_fa} is empty after stripping non-ACGT; skipping FP-family correction." >&2
    return 0
  fi

  resolve_merged_tools_dir
  ensure_kmer2ltr_dir
  local KMER2LTR_PY="${TOOLS_DIR}/Kmer2LTR/Kmer2LTR.py"
  local FLAG_FP_PY="${TOOLS_DIR}/Kmer2LTR/flag_fp_families.py"

  # (2) One clustering pass over the pooled elements -> the shared family basis.
  set -x
  python "$KMER2LTR_PY" \
    -i "$all_ltr_fa" \
    -o "${RUN_PREFIX}_all_ltr" \
    --ltr-cluster --internal-cluster \
    -p "$THREADS" \
    --min-seq-id 0.75
  set +x

  # Resolve the produced cluster tables by glob (robust to id-tag formatting);
  # a single --min-seq-id yields exactly one of each.
  shopt -s nullglob
  local -a cons_cluster=( "${RUN_PREFIX}"_all_ltr.consensus_id*_cluster.tsv )
  local -a int_cluster=(  "${RUN_PREFIX}"_all_ltr.internal_id*_cluster.tsv )
  shopt -u nullglob
  cons_fa="${RUN_PREFIX}_all_ltr.consensus.fa"
  (( ${#cons_cluster[@]} == 1 )) || die "Expected exactly 1 consensus cluster TSV, found ${#cons_cluster[@]}"
  (( ${#int_cluster[@]}  == 1 )) || die "Expected exactly 1 internal cluster TSV, found ${#int_cluster[@]}"
  [[ -s "$cons_fa" ]] || die "Missing consensus FASTA: $cons_fa"

  # (3) Flag FP families, purge their rows from every depth TSV, and -- only if
  #     the FP fraction exceeds --fp-mask-threshold -- write each genome's
  #     masked FASTA, which the orchestrator detects to trigger a re-run.
  run_fp_stage "${cons_cluster[0]}" "${int_cluster[0]}" "$cons_fa" \
               "$FLAG_FP_PY" "${depth_tsvs[@]}"

  # (4) Emit an FP-cleaned FASTA next to each cleaned TSV.
  write_clean_fastas "${depth_tsvs[@]}"
}

# True when this genome's mask actually changed bases.
# flag_fp_families.py writes --masked-out unconditionally once Stage C starts,
# even when it masks 0 bp. Gating on file existence alone therefore re-runs
# forever on an unchanged genome, burning every FP round before warning "not
# converged". Gate on the bp count it reports instead. An unparseable count
# means re-run, which is the safe direction.
genome_needs_rerun() {
  local prefix="$1" log="$2" bp
  [[ -s "./${prefix}_FP_masked.fa" ]] || return 1
  bp="$(sed -n 's/^\[INFO\] masked \([0-9,]*\) bp.*/\1/p' "$log" 2>/dev/null \
        | tail -1 | tr -d ',')"
  if [[ -z "$bp" ]]; then
    echo "WARNING: could not parse masked bp for ${prefix}; assuming it changed." >&2
    return 0
  fi
  (( bp > 0 ))
}

# Reuse a genome's detection outputs when its input did not change. .work dirs
# are hardlinked (bulk data, never rewritten downstream -- only read by
# reconcile_nests.py and the miniprot GFF lookup). Everything else is a real
# copy, because ltr_annotate.py rewrites the depth TSVs IN PLACE and a hardlink
# there would corrupt the previous attempt's directory. The input symlink is
# re-created per attempt, and _clean_/_FP_masked outputs belong to the attempt
# that produced them, so all three are skipped.
carry_forward_genome() {
  local prev="$1" adir="$2" prefix="$3" e name
  shopt -s nullglob
  local entries=( "${prev}/${prefix}"* )
  shopt -u nullglob
  for e in "${entries[@]}"; do
    name="$(basename "$e")"
    # Skip anything the next attempt regenerates or owns itself. The last three
    # patterns only overlap when --out_prefix makes RUN_PREFIX equal to a
    # genome's own prefix, but a stale cluster table or FP log leaking into the
    # next attempt would be a very quiet way to get wrong families.
    case "$name" in
      "${prefix}.input_genome.fa"*|"${prefix}_FP_masked.fa") continue;;
      *_clean_ltr.*|*_fpcheck*|*_all_ltr.*) continue;;
    esac
    if [[ -d "$e" && ( "$name" == *.work || "$name" == *_tools ) ]]; then
      cp -al "$e" "${adir}/"
    else
      cp -a "$e" "${adir}/"
    fi
  done
  echo "  reusing detection outputs for ${prefix} (input unchanged)"
}

# Top-level driver. Runs one or more fully-isolated attempts, each in its own
# staging directory. An attempt detects every genome (sequentially, skipping any
# whose input is unchanged), then runs the pooled merged stage over all of them.
# Genomes that flag_fp_families.py actually masked are re-detected on their
# masked FASTA in the next attempt; when none were masked -- or attempts run out
# -- the attempt is final, gets annotated, and is promoted to the user's
# working directory.
run_fp_orchestrator() {
  local final_dir staging abs_self abs_proteins abs_script
  final_dir="$PWD"
  abs_self="$(abspath "${BASH_SOURCE[0]}")"
  abs_script="$(abspath "$SCRIPT_PATH")"
  [[ -n "$PROTEINS" ]] && abs_proteins="$(abspath "$PROTEINS")"

  staging="${final_dir}/${RUN_PREFIX}_FPstaging"
  rm -rf "$staging"                 # defensive: clear a stale/partial staging dir
  mkdir -p "$staging"

  echo ""
  echo "############################################################"
  echo "FP-correction orchestrator: up to ${MAX_FP_ROUNDS} attempt(s)"
  echo "  genomes:     ${N_GENOMES}"
  echo "  staging dir: ${staging}"
  echo "############################################################"

  local attempt=1 final_adir="" prev_adir="" adir rc i p is_final any_rerun
  local -a cur_genome=() reuse=() abs_genomes=()
  for i in "${!GENOMES[@]}"; do
    abs_genomes+=( "$(abspath "${GENOMES[$i]}")" )
    cur_genome+=( "${abs_genomes[$i]}" )
    reuse+=( false )
  done

  while (( attempt <= MAX_FP_ROUNDS )); do
    adir="${staging}/round${attempt}"
    mkdir -p "$adir"

    echo ""
    echo "======================= FP round ${attempt} ======================="
    echo "  work dir:     ${adir}"

    for i in "${!GENOMES[@]}"; do
      p="${OUT_PREFIXES[$i]}"
      ln -sf "${cur_genome[$i]}" "${adir}/${p}.input_genome.fa"

      if [[ "${reuse[$i]}" == true && -n "$prev_adir" ]]; then
        carry_forward_genome "$prev_adir" "$adir" "$p"
        continue
      fi

      echo "  detecting ${p}: ${cur_genome[$i]}"

      # Replay the original argv, then append overrides so the parser's last-wins
      # semantics pin genome/proteins/out_prefix/script_path to isolated, absolute
      # values. The env var flips the re-exec'd script into worker mode.
      local worker_args=( "${ORIG_ARGS[@]}"
          --genome "${p}.input_genome.fa"
          --out_prefix "$p"
          --script_path "$abs_script" )
      [[ -n "${abs_proteins:-}" ]] && worker_args+=( --proteins "$abs_proteins" )

      rc=0
      ( cd "$adir" && _LTRHARVEST_WRAPPER2_WORKER=1 bash "$abs_self" "${worker_args[@]}" ) || rc=$?
      if (( rc != 0 )); then
        die "Detection failed for genome '${GENOMES[$i]}' (prefix ${p}, FP round ${attempt}, exit ${rc}).
  A genome that yields no LTR-RTs aborts the whole run rather than producing a
  'cross-species' family set built from fewer species than were requested.
  Staging left at: ${adir}"
      fi
    done

    # Annotation is deferred: finality cannot be known until the FP logs exist,
    # and annotating an attempt whose outputs get discarded is pure waste.
    ( cd "$adir" && run_merged_stage )

    any_rerun=false
    for i in "${!GENOMES[@]}"; do
      p="${OUT_PREFIXES[$i]}"
      if ( cd "$adir" && genome_needs_rerun "$p" "${p}_fpcheck.log" ); then
        reuse[$i]=false; any_rerun=true
        echo "FP round ${attempt}: ${p} was masked -> re-detecting next attempt."
      else
        reuse[$i]=true
      fi
    done

    is_final=false
    if [[ "$any_rerun" == false ]]; then
      echo "FP round ${attempt}: no genome was masked -> final pass."
      is_final=true
    elif (( attempt >= MAX_FP_ROUNDS )); then
      echo "" >&2
      echo "WARNING: FP masking still triggered after ${MAX_FP_ROUNDS} attempt(s); not" >&2
      echo "converged. Promoting the last attempt's cleaned outputs anyway." >&2
      is_final=true
    fi

    if [[ "$is_final" == true ]]; then
      ( cd "$adir" && run_annotation_stage )
      final_adir="$adir"
      break
    fi

    for i in "${!GENOMES[@]}"; do
      p="${OUT_PREFIXES[$i]}"
      if [[ "${reuse[$i]}" == false ]]; then
        cur_genome[$i]="$(abspath "${adir}/${p}_FP_masked.fa")"
      fi
    done
    prev_adir="$adir"
    attempt=$((attempt + 1))
  done

  [[ -n "$final_adir" ]] || die "FP orchestrator produced no final attempt (internal error)"

  for p in "${OUT_PREFIXES[@]}"; do rm -rf "${final_adir}/${p}_tools"; done
  rm -rf "${final_adir}/${RUN_PREFIX}_tools"

  echo ""
  echo "Promoting final outputs: ${final_adir} -> ${final_dir}"
  promote_fp_outputs "$final_adir" "$final_dir"

  if [[ "$DEV_KEEP_FP_ROUNDS" == true ]]; then
    echo "[dev] --dev-keep-fp-rounds set: retaining staging dir ${staging}"
  else
    rm -rf "$staging"
  fi

  local elapsed=$((SECONDS - WRAPPER_START))
  printf 'FP-correction complete (%d attempt(s), total %d:%02d).\n' \
    "$attempt" $((elapsed / 60)) $((elapsed % 60))

  if [[ "$RUN_PLOTS" == true ]]; then
    echo ""
    echo "Post-completion plotting (disable with --no-plots)..."
    for i in "${!GENOMES[@]}"; do
      p="${OUT_PREFIXES[$i]}"
      if bash "${SCRIPT_PATH}/ltrharvest_plots.sh" \
           --prefix "$p" --genome "${abs_genomes[$i]}" \
           --indir "$final_dir" --script_path "$SCRIPT_PATH"; then
        echo "Plots written to: ${final_dir}/${p}_plots"
      else
        echo "[WARN] plotting stage reported failures for ${p}; annotation outputs are unaffected." >&2
      fi
    done
  fi

  echo ""
  echo "Done."
}

# ----------------------------
# Parse args
# ----------------------------
[[ $# -eq 0 ]] && { usage; exit 1; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --genome)
      # Variadic: consume every following token that is not a flag. RESETS the
      # list rather than appending, because run_fp_orchestrator replays the
      # original argv and then appends an overriding --genome; last must win.
      shift
      GENOMES=()
      while [[ $# -gt 0 && "$1" != --* ]]; do
        GENOMES+=( "$1" ); shift
      done
      (( ${#GENOMES[@]} > 0 )) || die "--genome requires at least one FASTA"
      ;;
    --proteins) PROTEINS="${2:-}"; shift 2;;
    --terminate_count) TERMINATE_COUNT="${2:-}"; shift 2;;
    --max-rounds) MAX_ROUNDS_OVERRIDE="${2:-}"; shift 2;;
    --script_path) SCRIPT_PATH="${2:-}"; shift 2;;
    --threads) THREADS="${2:-}"; shift 2;;
    --out_prefix) OUT_PREFIX="${2:-}"; shift 2;;
    --run-trf) RUN_TRF=true; shift;;
    --no-trf) RUN_TRF=false; shift;;
    --wfa-align) WFA_ALIGN=true; shift;;
    --keep-weak-hmm-pass2-matches) KEEP_WEAK_HMM_PASS2=true; shift;;
    --pass2-aligner) PASS2_ALIGNER="${2:-}"; shift 2;;
    --fp-mask-threshold) FP_MASK_THRESHOLD="${2:-}"; shift 2;;
    --no-plots) RUN_PLOTS=false; shift;;
    --dev-keep-fp-rounds) DEV_KEEP_FP_ROUNDS=true; shift;;      # hidden dev flag
    --dev-max-fp-rounds) MAX_FP_ROUNDS="${2:-}"; shift 2;;       # hidden dev flag
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

(( ${#GENOMES[@]} > 0 )) || die "--genome is required"
for _g in "${GENOMES[@]}"; do
  [[ -f "$_g" ]] || die "Genome not found: $_g"
done
N_GENOMES=${#GENOMES[@]}
GENOME="${GENOMES[0]}"
if [[ -n "$PROTEINS" ]]; then
  [[ -f "$PROTEINS" ]] || die "Proteins not found: $PROTEINS"
fi

# Validate FP-correction knobs.
[[ "$FP_MASK_THRESHOLD" =~ ^[0-9]*\.?[0-9]+$ ]] \
  || die "--fp-mask-threshold must be a number in [0,1] (got '$FP_MASK_THRESHOLD')"
awk -v x="$FP_MASK_THRESHOLD" 'BEGIN{exit !(x>=0 && x<=1)}' \
  || die "--fp-mask-threshold must be in [0,1] (got '$FP_MASK_THRESHOLD')"
[[ "$MAX_FP_ROUNDS" =~ ^[1-9][0-9]*$ ]] \
  || die "--dev-max-fp-rounds must be a positive integer (got '$MAX_FP_ROUNDS')"

# ----------------------------
# Derived names
# ----------------------------
declare -a GENOME_PREFIXES=() OUT_PREFIXES=()
for _g in "${GENOMES[@]}"; do
  GENOME_PREFIXES+=( "$(genome_prefix_of "$_g")" )
done
genome_prefix="${GENOME_PREFIXES[0]}"

if (( N_GENOMES == 1 )); then
  # Single genome: --out_prefix names this genome's outputs AND the family
  # namespace, exactly as before. tools_base keeps the historic rule so a run
  # without --out_prefix still gets ./<genome_prefix>_tools.
  if [[ -n "$OUT_PREFIX" ]]; then
    tools_base="$OUT_PREFIX"
  else
    tools_base="$genome_prefix"
    OUT_PREFIX="${genome_prefix}_LTRs"
  fi
  OUT_PREFIXES=( "$OUT_PREFIX" )
  RUN_PREFIX="$OUT_PREFIX"
  TOOLS_DIR="./${tools_base}_tools"
else
  # Multi genome: every genome keeps its own auto-derived prefix; --out_prefix
  # names only the shared pool (families, merged cluster tables, staging).
  RUN_PREFIX="${OUT_PREFIX:-merged}"
  for _p in "${GENOME_PREFIXES[@]}"; do
    OUT_PREFIXES+=( "${_p}_LTRs" )
  done
  OUT_PREFIX="${OUT_PREFIXES[0]}"
  TOOLS_DIR="./${RUN_PREFIX}_tools"
fi

if (( N_GENOMES > 1 )); then
  echo "Genomes:         ${N_GENOMES}"
  for _i in "${!GENOMES[@]}"; do
    echo "  [$((_i+1))] ${GENOMES[$_i]}  ->  ${OUT_PREFIXES[$_i]}"
  done
  echo "Run prefix:      ${RUN_PREFIX}"
else
  echo "Output prefix:   ${OUT_PREFIX}"
fi
echo "Tools directory: ${TOOLS_DIR}"

# Preflight guards. Multi-genome only, so the single-genome path is untouched
# and pays no extra genome read. The worker always receives exactly one genome,
# so these only ever fire in the orchestrator.
if (( N_GENOMES > 1 )); then
  check_unique_basenames
  check_unique_seqids
fi

if [[ -z "$SCRIPT_PATH" ]]; then
  SCRIPT_PATH="$(abspath_dir "${BASH_SOURCE[0]}")"
else
  SCRIPT_PATH="${SCRIPT_PATH%/}"
fi

LTRHARVEST="${SCRIPT_PATH}/ltrharvest5.py"
MASKLTR="${SCRIPT_PATH}/mask_ltr.py"
RECONCILER="${SCRIPT_PATH}/reconcile_nests.py"
ANNOTATOR="${SCRIPT_PATH}/ltr_annotate.py"
GFF3_WRITER="${SCRIPT_PATH}/ltr_tsv_to_gff3.py"
[[ -f "$LTRHARVEST" ]] || die "Missing: $LTRHARVEST"
[[ -f "$MASKLTR" ]] || die "Missing: $MASKLTR"
[[ -f "$RECONCILER" ]] || die "Missing: $RECONCILER"
[[ -f "$ANNOTATOR" ]] || die "Missing: $ANNOTATOR"
[[ -f "$GFF3_WRITER" ]] || die "Missing: $GFF3_WRITER"

# ----------------------------
# Orchestrator vs. worker
# ----------------------------
# The top-level invocation orchestrates the (possibly iterated) FP-masking
# re-runs, each isolated in its own staging dir. Every attempt re-execs this
# script with _LTRHARVEST_WRAPPER2_WORKER=1 set, which skips orchestration and
# runs exactly one detection pipeline in the current (isolated) directory.
if [[ -z "${_LTRHARVEST_WRAPPER2_WORKER:-}" ]]; then
  run_fp_orchestrator
  exit $?
fi

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

WEAK_HMM_OPTS=()
if [[ "$KEEP_WEAK_HMM_PASS2" == true ]]; then
  WEAK_HMM_OPTS=(--keep-weak-hmm-pass2-matches)
fi
PASS2_ALIGNER_OPTS=(--pass2-aligner "$PASS2_ALIGNER")

SIZE=500000
TESORTER_RULE="70-75-80"
TSD_PASS2="--tsd-pass2"
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
# Namespaced by prefix: several genomes share one staging directory, and a
# fixed name would let genome 2 clobber genome 1's pass-2 library.
temp_lib="temp_ltr_2pass_lib_${OUT_PREFIX}.fa"

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

  # Same-round nested inners get rewritten with this round's own IUPAC char so
  # the per-round _ltr.fa stays consistent with cross-round nesting (where the
  # inner is already IUPAC-masked via the masked input genome). Without this,
  # same-round outers feed chimeric references (inner seq embedded) into the
  # next round's temp_lib -> TEsorter pass-2.
  this_round_char="${IUPAC_SEQ[$((round-1))]}"

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
    --tools-dir "$TOOLS_DIR" \
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
    --nested-flank-min "$NESTED_FLANK_MIN" \
    --nested-base-min "$NESTED_BASE_MIN" \
    --same-round-inner-char "$this_round_char" \
    "${pass2_opts[@]}" \
    "${WFA_OPTS[@]}" \
    "${WEAK_HMM_OPTS[@]}" \
    "${PASS2_ALIGNER_OPTS[@]}" \
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

# A run that completed no rounds produced no library at all. Without this the
# reconciler below is skipped and the wrapper still exits 0, reporting success
# for a run whose very first round died (e.g. a missing helper tool).
if (( ${#completed_round_prefixes[@]} == 0 )); then
  die "No detection round completed; no LTR library was produced. See the round errors above."
fi

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

# Cleanup. TOOLS_DIR is deliberately left in place: the merged stage runs in the
# orchestrator and needs the Kmer2LTR clone that lives inside it. The
# orchestrator removes every tools dir once the merged stage is done.
rm -f "$temp_lib" 2>/dev/null || true

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
echo "Detection complete for ${OUT_PREFIX} (${WRAPPER_TIME})."
