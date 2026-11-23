# Benchmarks

The README's "How well it works" gives each benchmark's headline numbers; this page gives the rest. They
were measured with the standalone program each command began as, whose outputs synLTR reproduces on the runs
recorded from it (the README lists them), so the accuracy numbers should carry over. Time and memory are the earlier
programs' unless marked as synLTR's. The record behind each benchmark:
[legacy/README.md](../legacy/README.md#benchmark-records).

## Assemblies (`genotype`, `geneflow`)

Measured on 144 simulated genome sets and on real genomes:

- Gene flow found in all three real cases: *mexicana* into maize B73 (D_pav +0.290, Z 7.39), *Zea nicaraguensis*
  into *Z. mays* ssp. *huehuetenangensis* (Z 6.36), mandarin ancestry in sweet orange (Z 16.85). With other
  anchoring routes, flank lengths and jackknife blocks of 1-20 Mb, the maize Z stayed 6.0-8.4.
- In 60 simulated sets without gene flow the largest |Z| was 3.14 (none of the 30 with the default demography above
  3). With one sister annotated 10 or 30 points less completely (20 sets), D_pav gave no false signal.
- Power is limited: Z > 3 in 10 of 36 simulated gene-flow runs on a 50-Mb backbone. A Z below 3 does not show the
  absence of gene flow, and verdicts near 3 are fragile: a positive and a negative maize control gave Z 3.16 and
  3.20. `stats.tsv` repeats the headline D_pav at 1, 2, 5, 10 and 20-Mb blocks, and `findings.txt` warns when the
  verdict changes with the block size.
- Genotypes: 99.7% of high-tier P and E calls right in simulations. Real genomes are harder: 46% of B73's elements
  were resolved in TIL18, against 76% in simulations; elements on unplaced scaffolds end mostly unresolved (70-83%
  U in maize).
- Painting is coarse: it covered 4 of sweet orange's 5 known tracts (1.7-23.8 Mb) but cannot see maize's tracts of
  about 10 kb.
- Resources: four maize genomes (2.2-3.3 Gb) took 2 h 11 min at 32 threads on a busy shared server (an upper bound),
  19.6 GiB at peak across all processes; minimap2 indexes take about 2.4 GB per Gb of assembly (`work/anchors/*.mmi`,
  deletable after the run). Every annotated genome is genotyped against every other, so time grows with the square
  of the genome count (tested up to 7 genomes).

## LTR-RT libraries only (`geneflow --screen`)

- At HPC1, whose highland allele
  entered maize from *mexicana* (Barnes et al. 2022, *PNAS*, https://doi.org/10.1073/pnas.2100036119), the highland
  landrace Palomero Toluqueño was painted *mexicana*-like (with *mexicana* TIL25 as donor) and B73 was not. Inv4m was
  found within 0.03 and 0.3 Mb of its breakpoints. Across 7 trios of 6 genomes LTR-RT D tracked gene-space SNP D with
  an outgroup (k-mers r = 0.99, insertions r = 1.00).
- The every-pair evenness check, on simulated genome sets (2026-11-14 and 2026-11-17; records outside this
  repository; nothing was fitted to them): it passed all 118 evenly annotated runs (lowest value 0.947) and held back
  all 13 runs at the 30-point gap and 10 of 13 at the 10-point gap.
- Speed and memory, synLTR's screen against its earlier standalone version side by side on one server (2026-11-19;
  records outside this repository): the maize trio at 32 threads took 278-281 s and 12.43-12.48 GB at peak against
  503-513 s and 53.08 GB; two simulated trios at 8 threads took 96-117 s and 3.79-3.87 GB against 115-136 s and
  5.16-5.35 GB.

## Long reads (`genotype`, `source`)

The rest of the README's table. Truth: haplotype-resolved assemblies of two outcrossing plants; ranges span the two
reference haplotypes.

| | PacBio HiFi (*Cardamine amara*) | ONT (*Medicago sativa* subsp. *falcata*) |
|---|---|---|
| elements on one haplotype (F/N) at 20x: recall | 0.89-0.92 | 0.87-0.89 |
| calls with GQ >= 20 at 20x: correct (share of sites) | 97-98% (65-70%) | 95% (62-64%) |
| presence/absence in another individual, 10-20x: calls correct | 99.3-99.9% | 97.8-99.5% |
| run time, memory (8 threads) | 1.5 min, 0.14 GB (5,839 elements, 60x) | 16.5 min, 0.33 GB (24,203 elements, 40x) |

- Read GQ as a ranking: calls with GQ >= 20 are 95-98% correct, not 99%.
- Where two haplotypes are identical (inbred lines, identical-by-descent stretches), one copy or two rests on read
  depth and is uncertain below about 30x; highly divergent haplotypes may not align to the reference at all; only
  insertions present in the reference can be genotyped.
- `source` reports the nearest fixed copy, a proxy for the parent (the true parent may be unsampled or lost).

## Subgenomes (`subgenomes`)

Limits:

- When the merger is older than LTR-RT turnover no subgenome-specific lineages survive and the result is "no
  structure" (22 of the 122 benchmark genomes, all woody bamboos among them).
- Large unplaced scaffolds rich in LTR-RTs count as core sequences and can decide the structure as chromosomes do.
- Mergers older than about a third of the split are often undatable.
- Years assume 1.3e-8 substitutions per site per year unless `--mu` is given.
