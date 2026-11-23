# synLTR

synLTR reads the history of genomes from their LTR retrotransposon (LTR-RT) insertions. An insertion is a derived
state (its absence is ancestral), its target-site duplication (TSD) marks one insertion event, and its two LTRs were
identical when it inserted, so their divergence dates it.

The scripts this repository held before version 0.2.0 (gene-synteny anchors, LTR-RT discovery with nesting, solo
LTRs) are kept, not maintained, in [legacy/](legacy/README.md).

| Command | What it answers | Input |
|---|---|---|
| `synltr genotype` | Which annotated insertions does each genome carry? Every annotated LTR-RT called present or absent in every other assembly, each call checked in the assembly itself (absent only where the empty site with one TSD copy is seen) | 2+ assemblies, 1+ annotated |
| `synltr geneflow` | Did genes flow between genomes, and where? `genotype`, then relationships, D statistics, painting of donor-like stretches, inversions, how insertions were lost, and a figure | 3+ assemblies |
| `synltr geneflow --screen` | The same question from LTR-RT libraries alone: a fast screen | 3+ libraries |
| `synltr genotype` (reads) | What does each sample carry at each reference insertion: the full element, a solo LTR, the empty site, or no sequence, on each chromosome copy | a reference, its annotation, long-read BAMs |
| `synltr source` | Which copy did each new insertion come from? | a finished read genotyping run |
| `synltr subgenomes` | Which chromosomes form each subgenome of a polyploid, which are homoeologs, when the progenitors split and merged, where they exchanged segments | one genome's LTR-RT library and annotation |
| `synltr plot` | Redraws the figure of a finished `geneflow` or `subgenomes` run | a run folder |
| `synltr sim pop`, `sim nucleus`, `sim score` | Genomes with known truth (populations with gene flow; polyploid origins), and calls scored against the truth | a scenario |

How well it works, in brief (measured with each command's earlier standalone version;
[details](#how-well-it-works)):
- `genotype` on assemblies: 99.7% of high-tier P and E calls right in simulations.
- `geneflow`: gene flow found in all three real cases; largest |Z| 3.14 in 60 simulations without it.
- `geneflow --screen`: maize-*mexicana* gene flow found from three LTR-RT libraries alone, Z 4.4-8.4.
- `genotype` on reads: 93-94% (HiFi) and 91-92% (ONT) of genotypes right at 20x.
- `subgenomes`: 90 of 122 natural allotetraploids phased correctly.

## Install

From the repository root, on Linux or macOS, not Windows; tested on Linux.

```bash
mamba env create -f environment.yml    # env "synltr_pkg": Python 3.12, LASTZ, minimap2, samtools, msprime
mamba activate synltr_pkg
pip install --no-deps -e .             # this checkout
synltr --version
```

The conda environment is the tested route. With pip alone (Python 3.12 or later), `pip install -e ".[sim]"`
installs the Python dependencies (pip builds mappy, minimap2's Python binding, from source: it needs a C compiler and
zlib headers); then put LASTZ (assemblies), minimap2 (`source`) and samtools (the quick start) on PATH, e.g. from
bioconda. Tests: `pytest -m "not slow"`; `pytest` adds the slow end-to-end runs (about 50 min); a test skips when a
program it needs (LASTZ, minimap2, samtools, bcftools, poppler) or msprime is missing.

## Quick start

```bash
bash examples/quickstart.sh            # about 2 minutes at 4 threads; writes ./synltr_quickstart
```

It builds four toys with the test suite's builders (`tests/builders.py`) and runs every command on them: four genomes
from `sim pop` (`assemblies/`: A1, B1, C1 annotated, outgroup O1 not, 20% of A's lineages from C), libraries of a trio
R, S, D and an outgroup O (`libraries/`: R carries D's haplotype in three tracts), a diploid's long reads (`reads/`)
and an allotetraploid's library from `sim nucleus` (`polyploid/`). They only check that everything runs.

## Run

These are the quick start's lines, each run in its toy's folder:

```bash
synltr genotype samples.tsv -o out -t 4                    # assemblies: presence/absence calls
synltr geneflow samples.tsv -o out -t 4 --trio A1,B1,C1    # the same calls (reused), then gene flow
synltr plot out                                            # redraw the figure
synltr geneflow --screen samples.tsv -o screen -t 4        # LTR-RT libraries only
synltr genotype samples.tsv -o out -t 4                    # long reads on one reference
synltr source out -t 4                                     # source copies of the new insertions in out/
synltr subgenomes samples.tsv -o out -t 4                  # one polyploid genome's LTR-RT library
```

- `-t` is the total CPU budget across all processes (default 4); `-o` the output folder (default: the current folder;
  `source` writes into its RUN folder, `plot` has no `-o`, the `sim` commands require it); `-v` logs progress and
  sanity checks. Every option: `synltr <command> --help`.
- `geneflow --trio R,S,D` names the recipient, its sister and the donor (default: every trio of ingroup genomes, roles
  from the data); `--identity` adds the LTR divergence of shared insertions (slow). `plot --trio R,S,D` redraws a
  `geneflow` run's figure around another trio the run tested (default: the strongest; on a screen, the first).
  `genotype --pair A,B` runs one ordered pair, to spread a run over jobs. A rerun reuses finished work, keyed on its
  inputs, options and code.
- Exit status: 0 done; 1 an error (one `error:` line, the traceback with `-v`); 2 a usage error (a bad option or
  input file); 3 genotyping finished but the analyses failed; 130 interrupted (Ctrl-C). Warnings start with
  `warning:`.

## Samples sheet

Tab-separated, a header line, one row per genome. `.` or an empty field means none; relative paths resolve against
the sheet's folder; `#` starts a comment line. Labels use `A-Z a-z 0-9 _ . -`, without `__`, without `_` at either
end, and not `all`.

| Input | Columns | Example row |
|---|---|---|
| assemblies | `label fasta annotation`, optional `group role` (role `outgroup`) | `B73  B73.fa.gz  B73.gff3  maize  .` |
| LTR-RT libraries | `label library` (an LTR_peel prefix, a FASTA, or `FASTA,TABLE` with `seq_id` and any of `k2p tsd strand ltr5_len ltr3_len`; plain or gzip; record names `chrom:start-end[#class]`, 1-based closed, and bedtools' 0-based `name::chrom:start-end` converted), optional `role`; a sheet without genome FASTAs runs the screen | `B73  run/B73v5_LTRs  .` |
| long reads | `label fasta annotation role reads`: one row of role `reference` with its FASTA and annotation, then one row per sample with its BAM | `ref  ref.fa  ltr.gff3  reference  .` |
| a polyploid | one row: `label library annotation`, optional `fasta`: the LTR-RT FASTA (several comma-separated), LTR_peel's GFF3 or a K2P table, the genome FASTA with its `.fai` beside it or the `.fai` | `quinoa  LTRs.fa.gz  LTRs.gff3  genome.fa.fai` |

- Genome FASTA: plain or bgzip, not plain gzip; its `.fai` is built beside it, or in `work/` if its folder is read-only.
- Annotation: GFF3 ([LTR_peel](https://github.com/cwb14/LTR_peel), EDTA, LTR_retriever, LTRharvest, DANTE_LTR), BED
  (LTRs found by aligning each element's ends) or an LTR_peel prefix (`<prefix>_depth*_clean_ltr.tsv`). GFF3 and
  LTR_peel tables are 1-based closed, BED 0-based half-open. A genome without annotation is genotyped only; the first
  row sets the VCF's coordinates.
- Reads: coordinate-sorted, indexed BAMs that keep the supplementary alignments (`map-ont` for ONT reads), e.g.
  `minimap2 -ax map-hifi ref.fa reads.fq | samtools sort -o s1.bam - && samtools index s1.bam`. An indexed CRAM
  works too, written against the reference row's FASTA, which decodes it.

## How well it works

Each command began as a standalone program, and the benchmarks were measured with those programs. synLTR reproduces
their outputs on the runs recorded from them (maize chromosome 10, maize LTR-RT libraries, *Cardamine amara* HiFi
reads, five polyploids and a diploid control, simulations), apart from the changes listed in
[legacy/changes_from_original_tools.md](legacy/changes_from_original_tools.md), so the accuracy numbers should carry
over; time and memory are the earlier programs' unless marked as synLTR's. More numbers:
[docs/benchmarks.md](docs/benchmarks.md). The record behind each benchmark:
[legacy/README.md](legacy/README.md#benchmark-records).

**LTR-RT libraries only (`geneflow --screen`).** On maize B73, *mexicana* TIL18 and *parviglumis* TIL01 the four
statistics gave Z = 4.4-8.4 (an outgroup-based SNP test on the same assemblies: Z = 9.2), and chromosomes 5, 8 and 9
ranked lowest in *mexicana*-like signal by every recipient-side statistic, the three Yang et al. 2023 report as
depleted of *mexicana* ancestry (*Science*, https://doi.org/10.1126/science.adg8940). In sweet orange the haplotype
with 12.6% mandarin ancestry gave Z = 4.0-5.3, and every mandarin tract of at least 1.7 Mb in the phased genome of Wu
et al. 2023 (*Horticulture Research*, https://doi.org/10.1093/hr/uhac247) was painted.

Uneven annotation can fake the screen's statistics. False signals (|Z| > 3) in simulated genomes without gene flow
(simulated genome sets of 2026-11-14 and 2026-11-17; records outside this repository):

| Annotation | D_ins | D_kmer | D3 | D_self | D_pav (`geneflow`) |
|---|---|---|---|---|---|
| even (30) | 0 | 0 | 0 | 0 | |
| even, bottleneck (30) | 0 | 0 | 0 | 1 (\|Z\| 3.5) | |
| one sister 10 points less complete (10) | 9 | 5 | 0 | 0 | 0 |
| one sister 30 points less complete (10) | 9 | 9 | 4 | 10 | 0 |

- D_ins and D_kmer never give verdicts: no check separates the 10-point gap from even annotation reliably. They are
  reported as values with these rates.
- D3 and D_self held the null up to the 10-point gap, not beyond, so they give verdicts only if, for every pair X, Y,
  the smaller of r and 1/r is at least 0.889 (0.8/0.9); r: the share of X's old insertions found in Y over Y's in X.
- These rates are printed beside each statistic in `stats.tsv` and `findings.txt`. Held statistics stay in
  `stats.tsv`, flagged `held: annotation-sensitive` (D_ins, D_kmer) or `held: uneven annotation` (D3 and D_self when
  the check fails). Confirm any screen result with `geneflow` on the genome assemblies, whose D_pav gave no false
  signal at either gap. The screen's overview draws the painting panels only; the others say they need assemblies.

**Long reads (`genotype`, `source`).** Truth: haplotype-resolved assemblies of two outcrossing plants; ranges span the
two reference haplotypes.

| | PacBio HiFi (*Cardamine amara*) | ONT (*Medicago sativa* subsp. *falcata*) |
|---|---|---|
| genotypes correct, 10x / 20x / 40x | 89-91% / 93-94% / 94-96% | 88-89% / 91-92% / 93-94% |
| same at 20x: Sniffles2 / kanpig | 54-55% / 46-47% | 50-51% / 34-35% |
| new insertions (F/E) at 20x: recall / precision | 0.96-0.97 / 0.97-0.98 | 0.94-0.95 / 0.93-0.94 |

**Subgenomes (`subgenomes`).** 90 of 122 natural allotetraploids phased correctly; false structure in
7 of 67 control libraries (diploids and others); the right number of subgenomes in 128 of 129 simulations and in 35 of
40 polyploids with progenitor truth (33 perfect); split ages against gene Ks in 78 genomes, Spearman 0.599; mergers of
correctly phased genomes 32 dated, 30 bounded, 26 not datable (both measured one version
earlier).

## Outputs

Coordinates: 1-based closed in the tables, GFF3 and VCF, 0-based half-open in BED and bedGraph files, unless
[docs/outputs.md](docs/outputs.md) says otherwise. A finished `genotype`, `geneflow` or `subgenomes` run writes
`manifest.json` last, so it marks a finished run: versions, command, parameters, and the checksums of the input files
its stages keyed on (on assemblies each genome's FASTA and annotation, and each `--paf` file; for the screen the
library files; on read sets the reference FASTA; none for `subgenomes`, whose stages key on the library as read).
`work/` holds the cache a rerun resumes from, and a `.done.json` beside a file marks it finished, so a rerun skips it.

- `genotype`: `loci.tsv` (each insertion's states; 1-based closed), `pav.vcf.gz` (POS: the base before the insertion).
- `geneflow` (and `--screen`): `findings.txt`, `stats.tsv`, `tracts.<R>.bed` (0-based half-open), `overview.pdf`.
- Reads: `reads.genotypes.tsv`, `reads.sources.tsv` (1-based closed), `reads.vcf.gz` (POS: the base before the element).
- `subgenomes`: `chromosomes.tsv` (subgenomes and homoeologs), `dating.tsv` (split and merger ages), `report.pdf`.
- Before sharing outputs: `manifest.json`, `overview.pdf` and `report.*` record the host name and absolute paths, and
  `findings.txt`, `summary.json`, `reads.vcf.gz` and `sim pop`'s `params.json` record absolute paths.

Every file of every command, with what it holds and its coordinates: [docs/outputs.md](docs/outputs.md).

## Reading the results

- Assembly states: P intact copy, S solo LTR, T truncated (present); E empty site with one TSD copy (absent: the
  insertion never happened in that line); D deletion, O other sequence, G assembly gap, U unresolved (missing data,
  never absence). Tiers: high, medium, low; `annot` where the genome's own annotation holds the insertion.
- D_pav = (n_RD - n_SD) / (n_RD + n_SD), on the donor D's annotated insertions that are resolved in the recipient R
  and its sister S and empty in every outgroup: n_RD are present in R and empty in S, n_SD present in S and empty in
  R. Lineage sorting alone gives 0 on average, gene flow from D into R makes it positive. Z is D_pav over its
  block-jackknife SE; |Z| >= 3 counts, but a Z below 3 does not show the absence of gene flow. Treat every painted
  tract as a lead.
- The screen's statistics, each positive when D shares more with R than with S:
  - D_kmer: single-copy k-mers of D shared with R only, against those shared with S only.
  - D_ins: insertions R shares with D and S lacks, against those S shares with D and R lacks.
  - D3: sequence differences of ortholog triples, positive when R's copy is closer to D's than S's copy is.
  - D_self: ABBA against BABA sites in D's LTRs, each site's ancestral base read from the element's other LTR.

## Development

[AGENTS.md](AGENTS.md) is the guide for developers, human or AI: the layout, the checks to run, the rules the tests
enforce, and how to reuse an approach from `legacy/`.

## Citation and license

There is no paper yet. Please cite the software: Chris Benson, synLTR, version 0.2.0 ([CITATION.cff](CITATION.cff)).
GPL-3.0-or-later ([LICENSE](LICENSE)).
