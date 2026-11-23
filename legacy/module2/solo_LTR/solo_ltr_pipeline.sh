#!/usr/bin/env bash
# =============================================================================
# solo_ltr_pipeline.sh
#
# Identify candidate solo-LTRs by searching per-element LTR consensus sequences
# against a genome in which all intact (full-length) LTR-RTs are hardmasked.
#
#   (1) Hardmask intact LTR-RT element spans in the genome  -> masked genome
#   (2) mmseqs easy-search LTR consensus vs masked genome    -> raw hits (.m8)
#   (3) Filter hits for full-length LTR (solo-LTR rule)      -> solo-LTR loci
#
# Rationale
#   Every intact element (whole element span: both LTRs + internal region) is
#   replaced by N. Any surviving genomic match to an LTR consensus is therefore
#   a NON-intact copy. A match that spans the consensus end-to-end (+/- FLANK bp)
#   means the entire LTR is present  ==>  a bona fide solo-LTR (as opposed to a
#   degraded fragment, which only partially aligns).
#
# Design notes / interrogated mmseqs defaults (mmseqs v17):
#   --strand 2     default is 1 (FORWARD ONLY) -> would miss ~half of solo-LTRs
#   --mask 0       default is 1 (tantan low-complexity masking of target) ->
#                  masks the very repeats we are searching for; disable it
#   --search-type 3   nucleotide vs nucleotide
#   -s 7.5         max sensitivity preset (diverged repeats; small query set)
#   --cov-mode 2 -c 0 no coverage prefilter; the precise +/-FLANK end rule is
#                  applied authoritatively in filter_solo_ltr.py
#
# Coordinate convention (verified against *.fa.info): the intact-LTR FASTA
# headers ">chrom:START-END#LTR/Class/Fam" are 1-based, inclusive
#   (END - START + 1 == element length).  BED = START-1 .. END.
#
# Streaming / memory: bedtools maskfasta streams the genome (no whole-genome
# load); mmseqs and bedtools merge/intersect stream; the Python filter only
# reads the (small) hit table, never the genome.
#
# Resumable: each step is skipped if its output already exists, unless --force.
#
# Usage:
#   bash solo_ltr_pipeline.sh [options]
#   bash solo_ltr_pipeline.sh -h
# =============================================================================
set -Eeuo pipefail

# ----------------------------- defaults --------------------------------------
GENOME="gen5400000_final.fasta"
INTACT_FA="gen5400000_final.fasta_r1_ltr.fa"
CONSENSUS_FA="Kmer2LTR_run.consensus.fa"
OUTDIR="solo_ltr_out"
PREFIX="solo_ltr"

ALIGNER="mmseqs"         # mmseqs | minimap2  (Step 2 / 2b search engine)
MM2_BIN=""               # minimap2 binary (default: <script_dir>/minimap2/
                         #   minimap2, the IUPAC-aware fork; auto-cloned and
                         #   compiled on first use if absent)
MM2_PRESET="map-ont"     # minimap2 -x preset (grid-benchmarked: asm5/10/20
                         #   clip diverged solo-LTRs -> low recall; map-ont
                         #   tolerates the divergence. -N/-p don't matter.)
MM2_N=50                 # minimap2 -N  max secondary alignments/query (repeats
                         #   are high-copy: keep many, unlike the default 5)
MM2_P=0.1                # minimap2 -p  min secondary/primary score ratio (low:
                         #   keep diverged repeat copies)

FLANK=10                 # +/- bp tolerance at each consensus end (solo-LTR rule)
SENS=7.5                 # mmseqs -s sensitivity
EVALUE="1e-3"            # mmseqs -e   (raise to 1e-2/1e-1 for short/diverged LTRs)
MAXSEQS=4000             # mmseqs --max-seqs (avoid capping high-copy families)
MIN_PIDENT=0             # optional post-filter on % identity (0 = off; faithful
                         #   to the user's pure end-position definition)
MERGE_DIST=0             # bedtools merge -d for collapsing loci across consensi
LIB_PURGE=1              # 1 = also search intact LTR-RT library and purge
                         #   fragmented elements (LTR + residual internal seq);
                         #   0 = skip (end-rule only)
IMARGIN=50               # a spanning library hit must extend >= this many bp
                         #   beyond a candidate (=internal seq) to purge it
LCOVER=0.9               # library hit must cover >= this fraction of candidate
LIB_MAXSEQS=10000        # --max-seqs for the library search (avoid missing the
                         #   spanning hit -> would under-purge)
TSD_LEN=5                # target-site-duplication length to require flanking a
                         #   solo-LTR (0 = disable). True solo-LTRs retain the
                         #   insertion TSD; library fragments do not.
TSD_WIN=1                # search +/- this many bp around each LTR boundary
                         #   (W=1 benchmarked precision-optimal vs PrinTE truth;
                         #    W=2 trades ~15pt precision for ~4pt recall)
TSD_MM=0                 # max mismatch between 5'/3' TSD copies (TSDs exact)
THREADS="$(nproc 2>/dev/null || echo 8)"
VERBOSE=0
FORCE=0
KEEP_TMP=0

usage() {
  cat <<EOF
solo_ltr_pipeline.sh - candidate solo-LTR discovery via hardmask + mmseqs

  -g FILE   genome FASTA                      [${GENOME}]
  -i FILE   intact LTR-RT FASTA (coords in header) [${INTACT_FA}]
  -q FILE   LTR consensus FASTA (query)       [${CONSENSUS_FA}]
  -o DIR    output directory                  [${OUTDIR}]
  -p STR    output file prefix                [${PREFIX}]
  -f INT    flank tolerance at each consensus end, bp [${FLANK}]
  -a STR    aligner: mmseqs | minimap2        [${ALIGNER}]
            (minimap2 = IUPAC-aware fork, auto-cloned+built on first use)
  -x STR    minimap2 -x preset                [${MM2_PRESET}]
  -n INT    minimap2 -N max secondary/query   [${MM2_N}]
  -P FLOAT  minimap2 -p min sec/pri ratio     [${MM2_P}]
  -s FLOAT  mmseqs sensitivity (-s)           [${SENS}]
  -e STR    mmseqs E-value (-e)               [${EVALUE}]
  -m INT    mmseqs --max-seqs                 [${MAXSEQS}]
  -d INT    bedtools merge distance for loci  [${MERGE_DIST}]
  -I FLOAT  optional min % identity post-filter (0=off) [${MIN_PIDENT}]
  -L        skip the intact-library fragment purge (end-rule only)
  -M INT    internal margin: lib hit must extend >= this bp past a
            candidate (into internal seq) to purge it [${IMARGIN}]
  -C FLOAT  library hit must cover >= this fraction of candidate [${LCOVER}]
  -T INT    TSD length to require flanking a solo-LTR (0=off) [${TSD_LEN}]
  -W INT    TSD search window +/- bp around each boundary [${TSD_WIN}]
  -X INT    max mismatch between 5'/3' TSD copies [${TSD_MM}]
  -t INT    threads                           [${THREADS}]
  -v        verbose (per-step progress + sanity checks)
  -F        force: recompute even if outputs exist
  -K        keep mmseqs tmp directory
  -h        show this help and exit

Outputs (in DIR):
  <prefix>.intact.bed            merged intact-element intervals (masked)
  <prefix>.masked_genome.fasta   genome with intact LTR-RTs -> N
  <prefix>.raw.m8                LTR-consensus mmseqs hits (BLAST-tab+header)
  <prefix>.intact_lib.m8         intact LTR-RT library mmseqs hits (for purge)
  <prefix>.filtered.tsv          end-rule hits (+ solo_status, + tsd)
  <prefix>.loci.bed              FINAL solo-LTR loci (purged + TSD-required)
  <prefix>.loci.tsv              ALL candidate loci (solo_status, tsd_supported)
  <prefix>.purged_fragments.tsv  loci removed as fragmented LTR-RTs (+evidence)
  <prefix>.summary.txt           counts + per-family breakdown
EOF
}

while getopts ":g:i:q:o:p:f:a:x:n:P:s:e:m:d:I:M:C:T:W:X:Lt:vFKh" opt; do
  case "$opt" in
    g) GENOME="$OPTARG" ;;
    i) INTACT_FA="$OPTARG" ;;
    q) CONSENSUS_FA="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    p) PREFIX="$OPTARG" ;;
    f) FLANK="$OPTARG" ;;
    a) ALIGNER="$OPTARG" ;;
    x) MM2_PRESET="$OPTARG" ;;
    n) MM2_N="$OPTARG" ;;
    P) MM2_P="$OPTARG" ;;
    s) SENS="$OPTARG" ;;
    e) EVALUE="$OPTARG" ;;
    m) MAXSEQS="$OPTARG" ;;
    d) MERGE_DIST="$OPTARG" ;;
    I) MIN_PIDENT="$OPTARG" ;;
    M) IMARGIN="$OPTARG" ;;
    C) LCOVER="$OPTARG" ;;
    T) TSD_LEN="$OPTARG" ;;
    W) TSD_WIN="$OPTARG" ;;
    X) TSD_MM="$OPTARG" ;;
    L) LIB_PURGE=0 ;;
    t) THREADS="$OPTARG" ;;
    v) VERBOSE=1 ;;
    F) FORCE=1 ;;
    K) KEEP_TMP=1 ;;
    h) usage; exit 0 ;;
    \?) echo "ERROR: unknown option -$OPTARG" >&2; usage >&2; exit 2 ;;
    :)  echo "ERROR: option -$OPTARG requires an argument" >&2; exit 2 ;;
  esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
log()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2; }
vlog() { [ "$VERBOSE" -eq 1 ] && log "$*" || true; }
die()  { log "ERROR: $*"; exit 1; }

# --------------------------- minimap2 bootstrap ------------------------------
# The minimap2 aligner is the IUPAC-aware fork (cwb14/minimap2), not stock
# minimap2 (stock is IUPAC-blind in seeding -> drops IUPAC-dense LTR copies).
# It is not a system package, so the pipeline vendors+builds it itself rather
# than expecting a prebuilt binary. Idempotent/resumable like every other
# step: reuse the binary if present, resume an interrupted build if the source
# tree is already there, otherwise shallow-clone then `make`.
ensure_minimap2() {  # bin_path
  local bin="$1"
  local dir; dir="$(dirname "$bin")"
  local blog="${dir%/}.build.log"     # sibling of the source tree

  [ -x "$bin" ] && return 0           # already built -> reuse

  log "minimap2 (IUPAC-aware fork) not found at ${bin}; bootstrapping"

  # build prerequisites: fail fast, clear message (the script can vendor the
  # source but cannot install a toolchain -- that is the user's job).
  for t in git make; do
    command -v "$t" >/dev/null 2>&1 || die \
      "cannot build minimap2: '$t' not on PATH (need git + make + a C compiler).
       Install it, or place a prebuilt minimap2 at ${bin}."
  done
  command -v gcc >/dev/null 2>&1 || command -v cc >/dev/null 2>&1 || die \
    "cannot build minimap2: no C compiler (gcc/cc) on PATH.
     Install one, or place a prebuilt minimap2 at ${bin}."

  if [ ! -e "${dir}/Makefile" ]; then
    if [ -e "$dir" ]; then
      die "${dir} exists but has no Makefile (incomplete clone?).
           Remove it and re-run, or place a prebuilt minimap2 at ${bin}."
    fi
    log "cloning IUPAC-aware minimap2 fork (cwb14/minimap2) -> ${dir}"
    if ! git clone --depth 1 https://github.com/cwb14/minimap2.git "$dir" \
         > "$blog" 2>&1; then
      tail -n 20 "$blog" >&2 || true
      die "git clone of the minimap2 fork failed (full log: ${blog})"
    fi
  else
    log "minimap2 source already present in ${dir}; resuming build"
  fi

  log "compiling minimap2 (make) in ${dir}"
  if ! ( cd "$dir" && make ) > "$blog" 2>&1; then
    tail -n 20 "$blog" >&2 || true
    die "minimap2 'make' failed in ${dir} (full log: ${blog});
         build it manually or place a prebuilt minimap2 at ${bin}"
  fi
  [ "$VERBOSE" -eq 1 ] && { vlog "--- minimap2 build log ---"; cat "$blog" >&2; } || true

  [ -x "$bin" ] || die \
    "minimap2 build finished but ${bin} is missing/not executable (log: ${blog})"
  log "minimap2 bootstrap OK -> ${bin}"
}

# --------------------------- fail-fast checks --------------------------------
for tool in mmseqs bedtools samtools awk python3; do
  command -v "$tool" >/dev/null 2>&1 || die "required tool not found: $tool"
done
for f in "$GENOME" "$INTACT_FA" "$CONSENSUS_FA"; do
  [ -s "$f" ] || die "missing or empty input: $f"
done
[ -s "${SCRIPT_DIR}/filter_solo_ltr.py" ] || \
  die "filter_solo_ltr.py not found next to this script (${SCRIPT_DIR})"

case "$ALIGNER" in
  mmseqs) ;;
  minimap2)
    [ -n "$MM2_BIN" ] || MM2_BIN="${SCRIPT_DIR}/minimap2/minimap2"
    ensure_minimap2 "$MM2_BIN"
    [ -s "${SCRIPT_DIR}/paf2m8.py" ] || \
      die "paf2m8.py not found next to this script (${SCRIPT_DIR})"
    log "aligner=minimap2 ($("$MM2_BIN" --version 2>/dev/null)) preset=${MM2_PRESET} -N ${MM2_N} -p ${MM2_P}"
    ;;
  *) die "unknown --aligner '$ALIGNER' (use: mmseqs | minimap2)" ;;
esac

# Run the chosen aligner: QUERY vs TARGET -> mmseqs format-mode-4-style m8.
# Keeps Step 3 (filter) aligner-agnostic. $4 is the tmp tag/dir.
run_search() {  # query target out_m8 tmptag
  local q="$1" t="$2" out="$3" tag="$4"
  if [ "$ALIGNER" = "mmseqs" ]; then
    local tmp="${OUTDIR}/${PREFIX}.${tag}_tmp"
    [ "$FORCE" -eq 1 ] && rm -rf "$tmp"; mkdir -p "$tmp"
    # --dbtype 2 CRITICAL: IUPAC-rich consensus else auto-typed as protein
    # -> 0 hits, silently. --strand 2 / --mask 0 also non-default on purpose.
    local ms="$MAXSEQS"; [ "$tag" = "lib" ] && ms="$LIB_MAXSEQS"
    local extra=(); [ "$tag" = "lib" ] && extra=(--max-seq-len 100000)
    mmseqs easy-search "$q" "$t" "$out" "$tmp" \
      --dbtype 2 --search-type 3 --strand 2 --mask 0 \
      -s "$SENS" -e "$EVALUE" --min-seq-id 0.0 --cov-mode 2 -c 0.0 \
      --max-seqs "$ms" "${extra[@]}" --format-mode 4 \
      --format-output "query,target,qstart,qend,qlen,tstart,tend,tlen,alnlen,pident,mismatch,gapopen,evalue,bits,qcov" \
      --threads "$THREADS" -v "$([ "$VERBOSE" -eq 1 ] && echo 3 || echo 2)"
    [ "$KEEP_TMP" -eq 1 ] || rm -rf "$tmp"
  else
    # minimap2 (IUPAC-aware fork): -c for NM/AS tags; -N/-p tuned for
    # high-copy repeats; PAF -> m8 so the filter is unchanged. target first.
    local paf="${out%.m8}.paf"
    "$MM2_BIN" -c -x "$MM2_PRESET" -N "$MM2_N" -p "$MM2_P" \
      -t "$THREADS" "$t" "$q" \
      > "$paf" 2> >([ "$VERBOSE" -eq 1 ] && cat >&2 || cat >/dev/null)
    python3 "${SCRIPT_DIR}/paf2m8.py" "$paf" -o "$out"
    [ "$KEEP_TMP" -eq 1 ] || rm -f "$paf"
  fi
}

# Genome .fai (needed for fail-fast bounds check + bedtools); build if absent.
if [ ! -s "${GENOME}.fai" ]; then
  vlog "indexing genome (samtools faidx ${GENOME})"
  samtools faidx "$GENOME"
fi

mkdir -p "$OUTDIR"
INTACT_BED="${OUTDIR}/${PREFIX}.intact.bed"
MASKED="${OUTDIR}/${PREFIX}.masked_genome.fasta"
RAW_M8="${OUTDIR}/${PREFIX}.raw.m8"
LIB_M8="${OUTDIR}/${PREFIX}.intact_lib.m8"
TMPDIR_MM="${OUTDIR}/${PREFIX}.mmseqs_tmp"
TMPDIR_LIB="${OUTDIR}/${PREFIX}.mmseqs_lib_tmp"

log "solo-LTR pipeline start | genome=${GENOME} consensus=${CONSENSUS_FA} flank=${FLANK} lib_purge=${LIB_PURGE}"

# =============================================================================
# STEP 1 - Build intact-element BED and hardmask the genome
# =============================================================================
if [ "$FORCE" -eq 1 ] || [ ! -s "$INTACT_BED" ]; then
  vlog "STEP 1a: parsing intact LTR-RT headers -> BED (1-based incl -> 0-based half-open)"
  # Parse ">chrom:START-END#..." robustly. Fail-fast on coordinate-convention
  # or wrong-genome errors (out-of-range / unknown chrom): these are
  # foundational and must not be silently skipped.
  awk -v FAI="${GENOME}.fai" '
    BEGIN{
      while ((getline line < FAI) > 0){ split(line,a,"\t"); clen[a[1]]=a[2] }
    }
    /^>/{
      h=substr($0,2)
      hash=index(h,"#"); if(hash>0) locus=substr(h,1,hash-1); else locus=h
      ci=match(locus,/:[0-9]+-[0-9]+$/)
      if(ci==0){ print "ERROR: cannot parse coords from header: " h > "/dev/stderr"; exit 3 }
      chrom=substr(locus,1,ci-1)
      coords=substr(locus,ci+1)
      sp=index(coords,"-")
      s=substr(coords,1,sp-1)+0
      e=substr(coords,sp+1)+0
      if(!(chrom in clen)){ print "ERROR: header chrom not in genome .fai: " chrom > "/dev/stderr"; exit 3 }
      if(s<1 || e<s){ print "ERROR: bad interval (expect 1-based incl, s>=1, e>=s): " h > "/dev/stderr"; exit 3 }
      if(e>clen[chrom]){ print "ERROR: end " e " > " chrom " length " clen[chrom] " (wrong coord convention/genome?): " h > "/dev/stderr"; exit 3 }
      printf "%s\t%d\t%d\n", chrom, s-1, e        # BED: 0-based half-open
    }
  ' "$INTACT_FA" | sort -k1,1 -k2,2n > "${INTACT_BED}.unsorted.tmp"

  N_RAW=$(wc -l < "${INTACT_BED}.unsorted.tmp")
  [ "$N_RAW" -gt 0 ] || die "no intervals parsed from ${INTACT_FA}"

  # Merge overlapping / nested intervals -> clean union to mask.
  bedtools merge -i "${INTACT_BED}.unsorted.tmp" > "$INTACT_BED"
  rm -f "${INTACT_BED}.unsorted.tmp"
  N_MERGED=$(wc -l < "$INTACT_BED")
  MASK_BP=$(awk '{s+=$3-$2} END{print s+0}' "$INTACT_BED")
  log "STEP 1: ${N_RAW} intact intervals -> ${N_MERGED} merged; ${MASK_BP} bp to hardmask"
else
  log "STEP 1: reuse existing ${INTACT_BED} (use -F to force)"
fi

if [ "$FORCE" -eq 1 ] || [ ! -s "$MASKED" ]; then
  vlog "STEP 1b: bedtools maskfasta -mc N (streaming hardmask)"
  bedtools maskfasta -fi "$GENOME" -bed "$INTACT_BED" -fo "$MASKED" -mc N
  samtools faidx "$MASKED"
  if [ "$VERBOSE" -eq 1 ]; then
    EXP=$(awk '{s+=$3-$2} END{print s+0}' "$INTACT_BED")
    GOT=$(grep -v '^>' "$MASKED" | tr -cd 'N' | wc -c)
    vlog "SANITY: N's expected (>=, pre-existing N possible) ~${EXP}; observed in masked genome=${GOT}"
  fi
  log "STEP 1: hardmasked genome -> ${MASKED}"
else
  log "STEP 1: reuse existing ${MASKED} (use -F to force)"
fi

# =============================================================================
# STEP 2 - search LTR consensus (query) vs masked genome (target)
#   Aligner-agnostic (mmseqs or IUPAC-aware minimap2); output is always an
#   mmseqs-format-4-style m8 so STEP 3 is unchanged. The consensus carries
#   IUPAC codes: mmseqs needs --dbtype 2 (else auto-typed protein -> 0 hits);
#   the minimap2 fork is IUPAC-aware in seeding+alignment (recovers
#   IUPAC-dense copies the IUPAC-blind seeder drops).
# =============================================================================
if [ "$FORCE" -eq 1 ] || [ ! -s "$RAW_M8" ]; then
  log "STEP 2: ${ALIGNER} search: LTR consensus vs masked genome"
  run_search "$CONSENSUS_FA" "$MASKED" "$RAW_M8" "mmseqs"
  log "STEP 2: $(($(wc -l < "$RAW_M8") - 1)) raw hits -> ${RAW_M8}"
else
  log "STEP 2: reuse existing ${RAW_M8} (use -F to force)"
fi

# =============================================================================
# STEP 2b - search FULL intact LTR-RT library vs masked genome (for the purge)
#   A true solo-LTR carries ONLY the LTR; a truncated element leaves LTR +
#   residual internal sequence, which a single library hit spans AND that
#   reaches into the element's INTERNAL region (.info-based test in STEP 3).
# =============================================================================
if [ "$LIB_PURGE" -eq 1 ]; then
  if [ "$FORCE" -eq 1 ] || [ ! -s "$LIB_M8" ]; then
    log "STEP 2b: ${ALIGNER} search: intact LTR-RT library vs masked genome"
    run_search "$INTACT_FA" "$MASKED" "$LIB_M8" "lib"
    log "STEP 2b: $(($(wc -l < "$LIB_M8") - 1)) library hits -> ${LIB_M8}"
  else
    log "STEP 2b: reuse existing ${LIB_M8} (use -F to force)"
  fi
else
  log "STEP 2b: skipped (-L); fragment purge disabled"
fi

# =============================================================================
# STEP 3 - Filter for solo-LTRs (full-length LTR end rule) + merge + purge
# =============================================================================
log "STEP 3: filtering for solo-LTRs (qstart<=1+${FLANK} AND qend>=qlen-${FLANK})"
LIB_ARGS=()
if [ "$LIB_PURGE" -eq 1 ]; then
  LIB_ARGS=(--intact-lib-m8 "$LIB_M8" \
    --internal-margin "$IMARGIN" --lib-cover-frac "$LCOVER")
  # *.fa.info gives per-element LTR/internal boundaries -> the rigorous
  # internal-overlap purge (distinguishes 'more LTR' from real internal seq).
  if [ -s "${INTACT_FA}.info" ]; then
    LIB_ARGS+=(--intact-lib-info "${INTACT_FA}.info")
  else
    log "WARN: ${INTACT_FA}.info not found -> purge falls back to the coarser"
    log "WARN: genomic-extent proxy (may false-purge true solo-LTRs)."
  fi
fi
TSD_ARGS=(--genome "$GENOME" --tsd-len "$TSD_LEN" \
  --tsd-window "$TSD_WIN" --tsd-mismatch "$TSD_MM")
python3 "${SCRIPT_DIR}/filter_solo_ltr.py" \
  --m8 "$RAW_M8" \
  --intact-bed "$INTACT_BED" \
  --out-prefix "${OUTDIR}/${PREFIX}" \
  --flank "$FLANK" \
  --min-pident "$MIN_PIDENT" \
  --merge-dist "$MERGE_DIST" \
  ${LIB_ARGS[@]+"${LIB_ARGS[@]}"} \
  "${TSD_ARGS[@]}" \
  $([ "$VERBOSE" -eq 1 ] && echo --verbose)

log "DONE. Key outputs:"
log "  filtered hits     : ${OUTDIR}/${PREFIX}.filtered.tsv"
log "  FINAL solo-LTRs   : ${OUTDIR}/${PREFIX}.loci.bed"
[ "$LIB_PURGE" -eq 1 ] && \
log "  purged fragments  : ${OUTDIR}/${PREFIX}.purged_fragments.tsv"
log "  summary           : ${OUTDIR}/${PREFIX}.summary.txt"
