#!/usr/bin/env bash
set -euo pipefail

# Defaults (can be overridden by env or CLI)
PEP_DIR=${PEP_DIR:-orthogroups/pep}
CDS_DIR=${CDS_DIR:-orthogroups/cds}
ALIGN_DIR=${ALIGN_DIR:-align_prot}
CODON_DIR=${CODON_DIR:-codon_align}
THREADS=${THREADS:-4}
PROCESSES=${PROCESSES:-2}
TRIMAL_TRIM_MODE=${TRIMAL_TRIM_MODE:-gappyout}
FORCE=0

usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Options:
  --pep-dir DIR            Protein FASTA dir (default: $PEP_DIR)
  --cds-dir DIR            CDS FASTA dir (default: $CDS_DIR)
  --align-dir DIR          Protein alignment dir (default: $ALIGN_DIR)
  --codon-dir DIR          Codon alignment dir (default: $CODON_DIR)
  --threads N              Threads per MAFFT job (default: $THREADS)
  --processes N            Parallel orthogroups to run (default: $PROCESSES)
  --trim-mode MODE         trimal mode, e.g. gappyout, strict (default: $TRIMAL_TRIM_MODE)
  --force                  Overwrite existing outputs
  -h, --help               Show this help

Environment variables with same names also work (take lower priority than CLI).
EOF
}

# Parse CLI
while [[ $# -gt 0 ]]; do
  case "$1" in
    --pep-dir)        PEP_DIR="$2"; shift 2;;
    --cds-dir)        CDS_DIR="$2"; shift 2;;
    --align-dir)      ALIGN_DIR="$2"; shift 2;;
    --codon-dir)      CODON_DIR="$2"; shift 2;;
    --threads)        THREADS="$2"; shift 2;;
    --processes)      PROCESSES="$2"; shift 2;;
    --trim-mode)      TRIMAL_TRIM_MODE="$2"; shift 2;;
    --force)          FORCE=1; shift;;
    -h|--help)        usage; exit 0;;
    *) echo "Unknown option: $1" >&2; usage; exit 1;;
  esac
done

# Check deps
need() { command -v "$1" >/dev/null 2>&1 || { echo "Error: '$1' not found in PATH" >&2; exit 1; }; }
need mafft
need trimal
need xargs
need bash

# Create outputs
mkdir -p "$ALIGN_DIR" "$CODON_DIR"

# Worker function for a single OG basename (e.g., OG10504.fa)
process_og() {
  local pep_file="$1"
  local pep_base
  pep_base="$(basename "$pep_file")"
  local og="${pep_base%.fa}"
  local cds_file="$CDS_DIR/$og.fa"
  local aln_file="$ALIGN_DIR/$og.aln.fa"
  local codon_file="$CODON_DIR/$og.codon.trim.fna"

  if [[ ! -s "$cds_file" ]]; then
    echo "[skip] $og: missing CDS file: $cds_file" >&2
    return 0
  fi

  if [[ $FORCE -eq 0 && -s "$aln_file" && -s "$codon_file" ]]; then
    echo "[done] $og: outputs exist, skipping (use --force to rerun)"
    return 0
  fi

  echo "[run]  $og: MAFFT -> TRIMAL (threads/job=$THREADS ; trim=$TRIMAL_TRIM_MODE)"

  # MAFFT on protein after stripping trailing "*" from sequences
  mafft --localpair --maxiterate 1000 --thread "$THREADS" \
    <(sed 's/\*$//' "$pep_file") \
    > "$aln_file"

  # Back-translate and trim with trimal
  trimal \
    -in "$aln_file" \
    -backtrans "$cds_file" \
    -ignorestopcodon \
    -"$TRIMAL_TRIM_MODE" \
    -out "$codon_file"
}

export -f process_og
export PEP_DIR CDS_DIR ALIGN_DIR CODON_DIR THREADS TRIMAL_TRIM_MODE FORCE

# Collect protein fasta files
shopt -s nullglob
mapfile -d '' PEP_FILES < <(find "$PEP_DIR" -type f -name '*.fa' -print0 | sort -z)
shopt -u nullglob

if [[ ${#PEP_FILES[@]} -eq 0 ]]; then
  echo "No protein FASTA files found in $PEP_DIR (*.fa)" >&2
  exit 1
fi

# Run in parallel using xargs -P
printf '%s\0' "${PEP_FILES[@]}" \
  | xargs -0 -n1 -P "$PROCESSES" -I{} bash -c 'process_og "$@"' _ {}

echo "All tasks finished."
