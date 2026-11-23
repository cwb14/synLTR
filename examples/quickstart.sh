#!/usr/bin/env bash
# synLTR quick start: toy data with known truth, built by the test suite's builders (tests/builders.py), and every
# command run on it.
#
# Run with the synltr environment active and this checkout installed (README, Install):
#   bash examples/quickstart.sh [WORKDIR]     (default ./synltr_quickstart; WORKDIR must not exist yet)
# Needs lastz, minimap2, samtools and msprime, all in environment.yml. Nothing is downloaded. Uses 4 threads. The
# toys are smoke tests: far too small to measure anything.
set -euo pipefail

REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
W=${1:-synltr_quickstart}
T=4
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1   # numpy's math libraries stay single-threaded

for exe in python synltr lastz minimap2 samtools; do
  command -v "$exe" > /dev/null || { echo "error: $exe not on PATH; activate the environment (mamba activate synltr_pkg) and install this checkout (pip install --no-deps -e .)" >&2; exit 1; }
done
python -c "import msprime" 2> /dev/null || { echo "error: msprime not installed (synltr sim pop needs it); it is in environment.yml" >&2; exit 1; }
if [ -e "$W" ]; then
  echo "error: $W exists; give a new WORKDIR" >&2
  exit 1
fi

mkdir "$W"
cd "$W"
SECONDS=0
step() { printf '\n== [%4ds] %s\n' "$SECONDS" "$*"; }
build() { PYTHONPATH="$REPO/tests" python -c "from pathlib import Path; import builders; $1"; }

# ---------------------------------------------------------------- 1. assemblies: genotype, geneflow, plot, sim score
# A toy backbone (3 x 300 kb with 40 planted LTR-RTs, annotated in GFF3) and four genomes evolved on it by
# `synltr sim pop`: A1, B1 and C1 (annotated, 10% of elements missing from each annotation) and the outgroup O1
# (not annotated). A pulse moves 20% of A's lineages from C, 2,000 generations ago.
step "1. assemblies/: four genomes simulated with known truth (sim pop)"
mkdir backbone
build 'builders.write_backbone(Path("backbone/toy.fa"), Path("backbone/toy.gff3"))'
synltr sim pop --fasta backbone/toy.fa --annotation backbone/toy.gff3 \
  --region s1:1-300000 --region s2:1-300000 --region s3:1-300000 \
  --per-taxon 1 --f 0.2 --ta 2000 --seed 2 -t "$T" -o assemblies
cd assemblies
cat samples.tsv
step "   genotype: each annotated element called in every other genome -> out/"
synltr genotype samples.tsv -o out -t "$T"
step "   geneflow: the genotypes reused, then the analyses of recipient A1, sister B1, donor C1 -> out/"
synltr geneflow samples.tsv -o out -t "$T" --trio A1,B1,C1
step "   plot: the overview figure redrawn"
synltr plot out
step "   sim score: C1's elements called in A1, against the simulated truth -> score/"
synltr sim score --calls out/pairs/C1__A1.calls.tsv.gz --elements out/work/elements/C1.tsv \
  --truth truth/C1__A1.tsv -o score/C1__A1
cd ..

# ---------------------------------------------------------------- 2. libraries: geneflow --screen
# LTR-RT libraries only (FASTA records plus a table of K2P, TSD, strand and LTR lengths) of a trio R, S, D and an
# outgroup O. R carries D's haplotype in three tracts; D carries one inverted segment.
step "2. libraries/: a trio and an outgroup as LTR-RT libraries (geneflow --screen)"
build 'builders.trio_outgroup(Path("libraries"))'
cd libraries
synltr geneflow --screen samples.tsv -o screen -t "$T"
cd ..

# ---------------------------------------------------------------- 3. reads: genotype, source
# A diploid whose first haplotype is the reference: 15x long reads per haplotype, aligned with minimap2.
step "3. reads/: a diploid's long reads on its reference (genotype, source)"
build 'rs = builders.read_set(Path("reads").resolve()); builders.reads_sheet(Path("reads/samples.tsv"), rs.ref, rs.gff, {"sample1": rs.bam})'
cd reads
synltr genotype samples.tsv -o out -t "$T"
synltr source out -t "$T"
cd ..

# ---------------------------------------------------------------- 4. a polyploid: subgenomes
# An allotetraploid's LTR-RT library: progenitors split 2 Myr ago, merged 100 kyr ago, two chromosomes each.
step "4. polyploid/: an allotetraploid's LTR-RT library (sim nucleus, subgenomes)"
mkdir polyploid
build 'import json; Path("polyploid/scenario.json").write_text(json.dumps(builders.SUBGENOMES_TOY))'
cd polyploid
synltr sim nucleus scenario.json -o library --truthdir truth
printf 'label\tlibrary\tannotation\tfasta\ntoy\tlibrary/ltr.fa\tlibrary/annotation.gff3\tlibrary/genome.fai\n' > samples.tsv
synltr subgenomes samples.tsv -o out -t "$T"
cd ..

step "done in ${SECONDS} s; outputs in $PWD"
show() {  # show FILE COLUMN,COLUMN,...: those columns of a TSV with a header line (lines starting # skipped)
  awk -F'\t' -v OFS='\t' -v want="$2" '/^#/ {next} !n {n = split(want, w, ","); for (i = 1; i <= NF; i++) c[$i] = i}
    {s = $c[w[1]]; for (j = 2; j <= n; j++) s = s OFS $c[w[j]]; print s}' "$1"
}
echo "assemblies (a toy too small to show gene flow):"
grep -m1 'D_pav(A1,B1;C1) = ' assemblies/out/findings.txt | cut -c1-100 || true
echo "libraries:"
grep 'Verdict, ' libraries/screen/findings.txt || echo "- no verdict: see libraries/screen/findings.txt"
echo "reads:"
show reads/out/reads.genotypes.tsv element_id,sample,genotype,gq,pav
echo "polyploid:"
show polyploid/out/chromosomes.tsv sequence,subgenome,partner
show polyploid/out/dating.tsv subgenome_a,subgenome_b,event,estimator,reported,years | awk -F'\t' 'NR == 1 || $5 == "yes"'
