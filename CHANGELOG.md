# Changelog

## 0.2.0 (2026-11-23)

LTRtrace becomes synLTR: the package moved into the synLTR repository (github.com/cwb14/synLTR) and took its name.
No method or number changed.

- The command `ltrtrace` is now `synltr`, the Python package `ltrtrace` is `synltr`, and the conda environment in
  environment.yml is named `synltr_pkg` (not `synltr`, which differs only by case from the env `synLTR` of the
  earlier pipeline). Options, the samples sheet and every output file and column are unchanged.
- Outputs that name the program now name `synltr` and its version: the GFF3 source column (`synltr-0.2.0`), the
  VCF's `##source` line, the first line of `findings.txt` and `manifest.json`'s `tool`.
- A rerun in a folder written by LTRtrace 0.1.0 recomputes every stage (stage keys hold the version and the code).
- The environment variable `LTRTRACE_LONG_TESTS` (one long test) is `SYNLTR_LONG_TESTS`.
- The quick start writes `./synltr_quickstart` (was `./ltrtrace_quickstart`).
- The hand-built synLTR pipeline that this repository held before (gene-synteny anchors, LTR-RT discovery with
  nesting, solo LTRs, the numbered `ltrharvest_phasing` scripts, and the tarballs of LTRpav, LTRmix and ltrkin) is
  kept under `legacy/`; [legacy/README.md](legacy/README.md) indexes it. The package does not use it.

## 0.1.0 (2026-11-19)

Released as LTRtrace 0.1.0, with the command `ltrtrace`; this entry gives commands under their current name,
`synltr`. First version: four tools merged into one package with one command, `synltr`. On the runs recorded from each tool,
every number it reported is reproduced, apart from the changes below.

| Source | Version | Now |
|---|---|---|
| ltrharvest_phasing10 (`ltrharvest_phasing10.py` and its simulator `simallo.py`, in `legacy/module2/`) | 10.1.0 | `synltr subgenomes`; `synltr sim nucleus` |
| LTRmix (unpublished; `legacy/module2/LTRmix.tar.gz`, sha256 `1dc3439985a2e4c55a5ac72bccb71421a75792bcf6dd3cfd29e7937fcd46a181`) | | `synltr geneflow --screen` |
| LTRpav (`legacy/module2/LTRpav.tar.gz`) | 0.1.0 | `synltr genotype` and `synltr geneflow` on assemblies; `synltr sim pop`, `sim score` |
| ltrkin | 0.2.1 | `synltr genotype` on read sets; `synltr source` |

**Every command**
- One samples sheet (tab-separated, one row per genome) replaces each tool's input options; README.md gives its
  columns for each kind of input.
- `-t` is the total CPU budget across all processes (default 4); numba's thread pool is sized from it. `-o` names the
  output folder (default: the current folder; `source` writes into its RUN folder, `plot` has no `-o`, the `sim`
  commands require it), `-v` logs progress, `--version` prints the version.
- Reruns in the same folder reuse finished work, keyed on the inputs, options and code (LTRmix had no checkpoints).
- A finished `genotype`, `geneflow` or `subgenomes` run writes `manifest.json`: versions, command, parameters, and the
  checksums of the input files its stages keyed on (on assemblies each genome's FASTA and annotation, and each `--paf`
  file; for the screen the library files; on read sets the reference FASTA; none for `subgenomes`, whose stages key on
  the library as read).
- Outputs that name the program that wrote them name `synltr` and its version: the GFF3 source column
  (`synltr-<version>`; was `LTRpav`, `v9` or `simallo`), the VCF's `##source` line (was `LTRpav` or `ltrkin 0.2.1`) and
  the first line of `findings.txt` (was `ltrpav 0.1.0` in `geneflow`'s, `ltrharvest_phasing9 v9.1` in
  `subgenomes`'). Column and field names are unchanged (`##LTRpav_POS`, `v8_exchanges`).
- A figure that fails to draw never fails a run; `synltr plot RUN` redraws it.
- Compressed annotations, library tables and PAFs are told by their first bytes, not their names: a gzip file
  without `.gz` reads as gzip and a plain one named `.gz` as text (these readers went by the name).
- Exit status: 0 done, 1 error, 2 usage error, 3 genotyping finished but the analyses failed, 130 interrupted (Ctrl-C:
  one line, `error: interrupted`). Warnings start with `warning:`.

**ltrharvest_phasing10 -> `synltr subgenomes SHEET`**
- `--ltr_fasta`, `--annotation` and `--fai` become the sheet's `library` (one or more FASTAs, comma-separated),
  `annotation` (one file) and `fasta` (the genome FASTA with its `.fai` beside it, or the `.fai`) columns.
  `--annotation` and `--fai` took several files each; the sheet takes one per genome, so concatenate them first.
  `--outdir` is `-o`; `-n/--n_subgenomes` is `-n`; `--assign`, `--mu`, `--seed`, `-t` and `-v` are unchanged, but a
  negative `--seed` is refused at once.
- The log is the package's: no `[v10 ...s]` prefix, warnings after `warning:`. A fault in the inputs (the library, its
  annotation, the `.fai`, `--assign`, `-n` above what the library can form) is a usage error: one `error:` line naming
  the file or option, exit 2 (was `ERROR:` and exit 1). `--assign` and `-n` are checked against the library before
  the k-mer stage (they were checked after it).
- Removed: `--label` (the sheet's label names the genome), the hidden `--dev_*` options, and the settings that were
  never the default (other split rules, `SPLIT_LTR_UNITS`, `MERGER_PAINT`, `MERGER_CAP="split"`, the off branches of
  `TOP_AGREE` and `AGE_HOMOEOLOGS`).
- New files: `report.png`, `report_numbers.tsv` (every number in the figure), `report_data.json` (what
  `synltr plot` redraws from).
- Fixed: with an LTRquest table as the annotation (not GFF3), minus-strand records whose two LTRs differ in length
  had their LTRs misplaced, because offsets counted along the genome were read as offsets along the record. The
  benchmark used GFF3 and is unaffected.

**LTRmix -> `synltr geneflow --screen SHEET`**
- `--genome LABEL=LIB` becomes one sheet row each (`label library`); `--outdir` is `-o`; `-t` now defaults to 4
  (was 8). `--no-figure` is removed. `ltrmix_plot.py OUTDIR` is `synltr plot RUN`; its `--mu` (an age axis in kyr
  on the insertion-age panel) is removed.
- Verdicts: D_ins and D_kmer never give one; they are reported as values, flagged `held: annotation-sensitive`. D3
  and D_self give verdicts only when the new every-pair evenness check passes (for each pair of genomes, the smaller of
  its ratio of old insertions found and that ratio's inverse must be at least 0.889), else they are flagged
  `held: uneven annotation`. Each statistic's false-signal rate in simulations is printed beside it; the completeness
  ratio is a value only.
- Outputs in the formats `geneflow` writes for assemblies, every number kept: `orthologs.<A>.<B>.tsv.gz` ->
  `pairs/<A>__<B>.calls.tsv.gz` and `pairs/<A>__<B>.div.tsv.gz`; `elements.<R>.tsv` -> `posterior.<R>.tsv`;
  `ancestry.<R>.bed` joins `tracts.<R>.bed` (column `alphabet`: ancestry or recent); new `decoy_posterior.<S>.tsv`;
  `stats.tsv` in the long format under layer `library`. The overview draws the painting panels; the others say they
  need assemblies.
- The k-mer step works in passes, each holding at most 8 GB of k-mer records; the screen keeps every library and its
  pair counters in memory, so its peak grows with the number of genomes. Side by side on one server, the maize trio
  at 32 threads peaked at 12.4 GB (LTRmix: 53 GB) and took 55-67% of LTRmix's time (four runs).
- Fixed: bedtools-style record names (`name::chrom:start-end`, 0-based start) were read as 1-based, so every position
  was 1 bp off.
- Fixed: records without coordinates were written to BED files with start -1; they never reach a BED file now.
- Fixed: LTRmix's README and `findings.txt` said D3 and D_self are robust to annotation dropout; they held the null
  only up to a 10-point completeness gap, and `findings.txt` now says so.

**LTRpav -> `synltr genotype` and `synltr geneflow`**
- `ltrpav run --samples S --outdir O` is `synltr geneflow S -o O`; `synltr genotype S -o O` stops after the
  genotyping. `ltrpav-plot [--trio R,S,D]` is `synltr plot RUN [--trio R,S,D]` (the lead trio, one the run tested);
  `ltrpav-sim pop` and `score` are `synltr sim pop` and `sim score` (`--outdir` and `--out` also as `-o`).
  `--threads` now defaults to 4 (was 8 in `ltrpav run`, 1 in `ltrpav-sim pop`).
- Removed: the single-stage commands `ltrpav-anchor`, `ltrpav-genotype`, `ltrpav-merge`, `ltrpav-stats` and
  `ltrpav-paint`. `synltr genotype --pair A,B` runs one pair job (as `ltrpav-genotype --pair A,B` did), and reruns
  skip finished stages, so `geneflow` with other options reruns only the analyses.
- Removed `--anchor`: one anchoring route remains (minimap2 on the element flanks), and `--paf DIR` (was
  `--anchor paf:<dir>`) takes your own whole-genome alignments. The `mm2`, `fastga`, `flank-blast` and `hybrid` routes
  are gone; on simulated pairs the flank route was as accurate as whole-genome minimap2 at about half its cost
  (LTRpav's docs/METHODS.md).
- Removed `--vcf-ref`: the VCF takes the first sheet row's coordinates (the old default).
- Removed `--block-size`: `stats.tsv` repeats the headline D_pav at 1, 2, 5, 10 and 20-Mb jackknife blocks (test
  `block_scan`), and `findings.txt` warns when the verdict changes with the block size.
- Removed the direction test: it read its prediction in 2 of 8 positive controls and called no direction in 118
  seven-genome runs.
- A bad sheet, annotation or `--pair` label is a usage error (exit 2; LTRpav exited 1), and so is a `sim` fault found
  before simulating or scoring: an unknown `--dropout-genome` label, an `-o` folder that is not empty, a `--f`, `--ta`,
  `--pulse` or `--bottleneck` value the demography refuses, `sim score` files that do not belong together.
- Fixed inside the run, outputs unchanged: an EDTA `ltr_identity`, converted to a JC69 distance, was held as a K2P
  value without a label. It is now labeled JC69 in the run cache (`work/elements/<label>.tsv.done.json`); `loci.tsv`'s
  `<label>_k` columns hold the same values as before, without the label.

**ltrkin -> `synltr genotype` on read sets, and `synltr source RUN`**
- `ltrkin genotype -r REF -a ANNOTATION -o PREFIX BAM...` becomes a sheet: one row of role `reference` (its `fasta`
  and `annotation`), then one row per sample (`reads`: its BAM). Samples are named by their sheet labels. Outputs in
  the `-o` folder: `reads.genotypes.tsv`, `reads.vcf.gz`, new `reads.unplaced.tsv`. `--ploidy`, `--min-mapq`, `-t`
  and `-v` are unchanged.
- `ltrkin source -r -a -g -o` is `synltr source RUN`: it finds the reference and annotation in the run's
  `manifest.json` and writes `reads.sources.tsv`, `reads.copies.tsv` and `reads.locality.tsv` there (or to `-o`).
  `-g` took several genotype tables; `synltr source` reads one run, so genotype all read sets in one sheet.
- Fixed: an individual with one solo LTR and one empty site was reported `absent` although one chromosome carries the
  insertion; it is `solo` now.
- Fixed in the VCF: END disagreed with SVLEN by the TSD length; an element starting at a sequence's first base got
  that base as POS (it now goes to `reads.unplaced.tsv`); contigs and records follow the reference's `.fai` order,
  not their names sorted as text (chr10 before chr2).
- Fixed in `source`: TSDs were compared as raw strings, so an inverted segmental duplicate escaped the duplicate
  check; its alignment cache ignored the code and parameters and could reuse stale alignments.
- Fixed: results came back in completion order, so which copy won an exact tie could change between runs.
