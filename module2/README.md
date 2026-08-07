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
| `--genome` | Genome FASTA (`.fa`/`.fasta`, plain or `.gz`). **Required.** |
| `--proteins` | **Can come from any related species** — it does not have to match the genome. Recommended. |
| `--threads` | CPU threads (default 20). |

Occasionally useful:

| Flag | Meaning |
|---|---|
| `--max-rounds 1` | Run a single detection round. Use this if you don't care about nested elements. |
| `--out_prefix` | Output prefix (default: `<genome>_LTRs`). |
| `--terminate_count` | Stop iterating when a round finds fewer than this many elements (default 100). |


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

Tab-separated, one row per element, 19 named columns
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
| 18 | `domains` | Protein domains with genomic coords: `DOMAIN\|Clade@start-end;...` e.g. `RT\|Bianca@910695-911483;INT\|Bianca@912348-912953`, or `.` |
| 19 | `nest_status` | Nesting relations: `nest-outer:chrom:s-e` (that element is inside me) and/or `nest-inner:chrom:s-e` (I am inside that element), `;`-joined; `.` if un-nested |


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

## 5. Families: `*_all_ltr.consensus_id0.75_cluster.tsv`

```
family_representative                 member
Chr2:102001-110500#LTR/Gypsy/Reina    Chr2:102001-110500#LTR/Gypsy/Reina
Chr2:102001-110500#LTR/Gypsy/Reina    Chr2:355010-361200#LTR/Gypsy/Reina
Chr2:900444-905810#LTR/Copia/Ale      Chr2:900444-905810#LTR/Copia/Ale
```

A **family = all rows sharing the same representative**.

## 6. Plots (`<prefix>_plots/`)

| Output | What it shows |
|---|---|
| `struct/<clade>_average.pdf` | Average element structure per classification clade with bootstrap 95% CI whiskers on feature positions |
| `struct/<clade>_individual.pdf` | Every element of the clade drawn individually, all domains shown |
| `struct/all_elements.pdf` | All elements on one page |
| `<prefix>_summary.pdf` | Multi-page summary |
| `<prefix>_TEGV.html` | Self-contained interactive genome browser (open in web browser) |

## 7. Benchmarks

On a simulated genome (PrinTE) with 4,468 true intact LTR-RTs, 20 threads,
scored at ≥90% reciprocal overlap:

| Tool | Runtime | TP | FP | FN | F1 |
|---|---|---|---|---|---|
| **this pipeline** | **5:18** | 4,328 | 170 | 140 | **0.966** |
| EDTA raw LTR module (LTR_retriever) | 20:28 | 2,085 | 0 | 2,383 | 0.64 |
