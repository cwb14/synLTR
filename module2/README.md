A pipeline that finds full-length **LTR retrotransposons (LTR-RTs)**, including
elements nested inside other elements.

---

## 1. Install

```bash
mamba env create -f module2/ltrharvest_env.yml
mamba activate ltrharvest4
```

## 2. Quick start (test data)

Verify your installation using the test data.  It should finish in a few minutes.

```bash
mkdir athal_test && cd athal_test

bash ../module2/ltrharvest_wrapper2.sh \
  --genome   ../test/Athal_tair10_chr2.fa.gz \
  --proteins ../test/Athal.pep.gz \
  --threads  20 --max-rounds 1
```

## 3. Usage

Most runs only need three flags:

| Flag | Meaning |
|---|---|
| `--genome` | Genome FASTA (`.fa`/`.fasta`, plain or `.gz`). **Required.** Accepts several, space-separated — see [§5](#5-multiple-genomes-shared-family-names). |
| `--proteins` | **Can come from any related species** — it does not have to match the genome. Recommended. |
| `--threads` | CPU threads (default 20). |

Occasionally useful:

| Flag | Meaning |
|---|---|
| `--max-rounds 1` | Run a single detection round. Use this if you don't care about nested elements. |
| `--out_prefix` | Output prefix (default: `<genome>_LTRs`). With several genomes it names the shared family namespace instead — see [§5](#5-multiple-genomes-shared-family-names). |
| `--terminate_count` | Stop iterating when a round finds fewer than this many elements (default 100). |
| `--run-sdust` | Drop candidates made mostly of low-complexity sequence, early. Off by default. |

## 4. Primary outputs: `depth<N>_clean_ltr.{tsv,fa}`

Elements are bucketed by nesting depth. `N` =
how many layers of LTR-RT are nested *inside* the element:

```
depth0:  [===============]                     no LTR-RT inside ("un-nested")
depth1:  [=====[depth0]=====]                  one LTR-RT nested inside
depth2:  [===[depth1 [depth0] ]===]            two layers nested inside
...
```

For the test run you'd get, e.g.:

```
Athal_tair10_chr2_LTRs_depth0_clean_ltr.tsv   ← table of un-nested elements
Athal_tair10_chr2_LTRs_depth0_clean_ltr.fa    ← their sequences
Athal_tair10_chr2_LTRs_depth1_clean_ltr.tsv   ← single-nested elements (ie, one LTR-RT nested inside)
Athal_tair10_chr2_LTRs_depth1_clean_ltr.fa
...
```

### 4.1 `depth<N>_clean_ltr.tsv` format

Tab-separated, one row per element, 21 named columns
(`#`-prefixed header line):

| # | Column | Meaning |
|---|---|---|
| 1 | `name` | `chrom:start-end#Class/Superfamily/Clade` — 1-based inclusive genomic coordinates + TEsorter2 classification, e.g. `Chr2:102001-110500#LTR/Gypsy/Reina` |
| 2 | `LTR_len` | Length (bp) of the element's LTR |
| 3 | `aln_len` | Aligned length of the 5′-LTR vs 3′-LTR alignment |
| 4 | `subs` | Substitutions between the two LTRs |
| 5 | `ti` | Transitions |
| 6 | `tv` | Transversions |
| 7 | `raw_d` | Raw (p) distance between the two LTRs |
| 8 | `raw_T` | Age (years) from `raw_d` |
| 9 | `JC69_d` | Jukes–Cantor corrected distance |
| 10 | `JC69_T` | Age (years) from JC69 |
| 11 | `K2P_d` | Kimura 2-parameter distance — **the standard divergence estimate** |
| 12 | `K2P_T` | **Insertion age (years)** = `K2P_d / (2μ)`, default μ = 3×10⁻⁸ subs/site/year |
| 13 | `left_trim` | bp trimmed off the 5′ end by the WFA boundary detector (boundary over-extension) |
| 14 | `right_trim` | bp trimmed off the 3′ end |
| 15 | `ltr5_end` | Last bp of the 5′ LTR (1-based, relative to the element) |
| 16 | `ltr3_start` | First bp of the 3′ LTR (1-based, relative to the element) |
| 17 | `tsd` | Target-site duplication motif, or `.` if none found |
| 18 | `strand` | `+`, `-`, or `.` if undetermined — see 4.1.1 |
| 19 | `family` | `<prefix>_fam00001` … — see section 5 |
| 20 | `domains` | Protein domains with genomic coords: `DOMAIN\|Clade@start-end;...` e.g. `RT\|Bianca@910695-911483;INT\|Bianca@912348-912953`, or `.` |
| 21 | `nest_status` | Nesting relations: `nest-outer:chrom:s-e` (that element is inside me) and/or `nest-inner:chrom:s-e` (I am inside that element), `;`-joined; `.` if un-nested |

### 4.2 `depth<N>_clean_ltr.fa` format — nested regions are masked

Each record is the **full-length LTR-RT**, with every
nested inner element **hard-masked** by an IUPAC letter chosen by the *inner
element's own depth*:

| Inner element's depth | Mask letter |
|---|---|
| 0 | `N` |
| 1 | `R` |
| 2 | `D` |
| 3+ | `Y`, `S`, `W`, `K`, `M`, `B`, `H` |

So a `depth1` record looks like `ACGT...NNNN...ACGT` (its depth0 insert masked
with `N`), and a `depth2` record looks like `ACGT...RRR..NNN..RRR...ACGT`
(its depth1 child masked `R`, whose own depth0 child is masked `N`). 

**To remove the masked inner sequence** and keep only the outer element's sequence:

```bash
awk '/^>/{if(s)print s; s=""; print; next}{gsub(/[^ACGTacgt]/,""); s=s $0}END{if(s)print s}' \
  Athal_tair10_chr2_LTRs_depth1_clean_ltr.fa > depth1_removed_nested.fa
```

`depth0` records contain no masked nest-ins, so this is only needed for depth ≥ 1.

## 5. Multiple genomes: shared family names

Pass several genomes to give closely-related species **one family vocabulary**,
which is what makes between-species family comparisons meaningful:

```bash
bash ../module2/ltrharvest_wrapper2.sh \
  --genome   ../test/Athal_tair10_chr2.fa.gz Alyrata_chr2.fa \
  --proteins ../test/Athal.pep.gz \
  --threads  40
```

Each genome is detected on its own, then **all** detected elements are pooled
for a single clustering pass, so `family` means the same thing everywhere:

```
Athal_tair10_chr2_LTRs_depth0_clean_ltr.tsv     family = merged_fam00001
Alyrata_chr2_LTRs_depth0_clean_ltr.tsv          family = merged_fam00001
Athal_tair10_chr2_LTRs_all_depth_LTR_cleaned.gff3
Alyrata_chr2_LTRs_all_depth_LTR_cleaned.gff3
merged_all_ltr.consensus_id0.75_cluster.tsv     ← the one shared table
```

False-positive families are also called over the pooled set rather than one
species at a time, so a repeat that looks convincing in isolation but wrong
across species gets caught.

Naming: every genome keeps its own `<basename>_LTRs` prefix. `--out_prefix`
names only the shared pool — `--out_prefix Arabidopsis` gives
`Arabidopsis_fam00001` and `Arabidopsis_all_ltr.*`, default `merged`. With a
single genome the two are the same thing, so nothing changes.

One `--proteins` file serves every genome.

> **Sequence IDs must be unique across genomes.** Two files both calling a
> chromosome `Chr2` are rejected before any work starts: pooled clustering keys
> elements on `chrom:start-end`, so a shared ID would cross-assign families
> *and* cross-purge real elements between species. Rename first:
>
> ```bash
> awk '/^>/{sub(/^>/,">Aly_")}1' Alyrata.fa > Alyrata.renamed.fa
> ```

## 6. GFF3 annotation

Two files, pooled across all depths and built from the FP-purged
`_clean_ltr.tsv` set (falling back to the raw set, with a warning, if the FP
stage never ran):

| Output | Contents |
|---|---|
| `<prefix>_all_depth_LTR_cleaned.gff3` | The LTR-RTs |
| `<prefix>_all_depth_protein_LTR_cleaned.gff3` | The same, plus every miniprot protein alignment. Only written when the run had `--proteins` |

Each element is one block:

```
chr  synLTR  LTR_retrotransposon   13031  17307  .  -  .  ID=…_LTRRT_00001;Name=chr:13031-17307;
                                                            classification=LTR/Gypsy/Tekay;superfamily=Gypsy;clade=Tekay;
                                                            family=…_fam00001;family_size=5;depth=0;
                                                            K2P_d=0.026305;K2P_T=438418;strand_source=tesorter
chr  synLTR  long_terminal_repeat  13031  14048  .  -  .  ID=…_LTRRT_00001.lLTR;Parent=…_LTRRT_00001
chr  synLTR  long_terminal_repeat  16290  17307  .  -  .  ID=…_LTRRT_00001.rLTR;Parent=…_LTRRT_00001
chr  synLTR  protein_match         14482  14712  .  -  .  ID=…_LTRRT_00001.CHD.1;Parent=…;Name=CHD;clade=Tekay
###
```

Notes:

- Nested elements are **flat top-level features**, not children of their host:
  GFF3 `Parent` means part-of, and a nested LTR-RT is not a part of the element
  it landed in. Use the `depth` and `nest_status` attributes instead.
- Blocks are coordinate-sorted, and lines within a block are too, but
  **overlapping blocks are not interleaved** — so a nested element's block
  follows its host's in full.

## 7. Plots (`<prefix>_plots/`)

| Output | What it shows |
|---|---|
| `struct/<clade>_average.pdf` | Average element structure per classification clade with bootstrap 95% CI whiskers on feature positions |
| `struct/<clade>_individual.pdf` | Every element of the clade drawn individually, all domains shown |
| `struct/all_elements.pdf` | All elements on one page |
| `<prefix>_summary.pdf` | Multi-page summary |
| `<prefix>_TEGV.html` | Self-contained interactive genome browser (open in web browser) |

## 8. Benchmarks

On a simulated genome (PrinTE) with 4,468 true intact LTR-RTs, 20 threads,
scored at ≥90% reciprocal overlap:

| Tool | Runtime | TP | FP | FN | F1 |
|---|---|---|---|---|---|
| **this pipeline** | **5:18** | 4,328 | 170 | 140 | **0.966** |
| EDTA raw LTR module (LTR_retriever) | 20:28 | 2,085 | 0 | 2,383 | 0.64 |
