# legacy/

The hand-built synLTR pipeline ("Discovery and characterization of syntenic LTR-RTs for phylogenetic analyses"), kept intact and runnable in place as a source of approaches for future work. It is not the user-facing tool (that is the `synltr` package in `src/synltr`) and it is not covered by the package's tests. All paths below are relative to `legacy/`.

Every file is as it was at the repository root before version 0.2.0, moved here without a change, so `git log --follow --first-parent legacy/<file>` and `git blame` show its history (the move is part of a merge commit, which plain `--follow` does not look into). Usage lines inside the scripts still say `synLTR/module2/...`: read them as `legacy/module2/...`. The wrappers find their helper scripts from their own location or from `--script_dir`, so they run from here. This file is the only one in `legacy/` to keep current: when an approach is ported into the package, set its status to `ported`.

## Map

| Piece | Question it answers |
|---|---|
| `module1.py`, `module1/` | Where are the syntenic blocks between each pair of genomes? Proteins are lifted onto every genome, JCVI finds gene anchors, anchors become cleaned block coordinates. |
| `module1/util/` | Which chromosomes are homologous across genomes (rename and reorient them), and what is the species tree from syntenic orthogroups? |
| `module2.py` + helpers in `module2/` | Where are the full-length LTR-RTs in one genome? LTRharvest + LTR_FINDER candidates, boundary check by WFA alignment, Kmer2LTR divergence, TEsorter classification, overlap purge. |
| `module2/ltrharvest_wrapper2.sh` + `ltrharvest5.py` and helpers | Where are the full-length LTR-RTs, including elements nested inside other elements, with nesting depth, family, strand, GFF3 and plots? Documented in `module2/README.md`. |
| `module2/ltrharvest_phasing*.py`, `module2/simallo.py` | Which chromosomes form each subgenome of a polyploid, from its LTR-RT library alone; a simulator with known truth. |
| `module2/solo_LTR/` | Where are the solo LTRs in a genome? Two search routes (mmseqs/minimap2 and BLAST) and a scorer against simulated truth. |
| `module3/modeling_divK_dev/` | Can a model correct divK, the optimum K that a piecewise regression reads from Kmer2LTR divergences of shared and unique insertions? Development scripts that generate training data on simulated genomes. |
| `module2/ltrkin.tar.gz`, `LTRmix.tar.gz`, `LTRpav.tar.gz` | Full source trees of three of the four tools merged into the package, with their design notes and benchmark records. |
| `test/` | *Arabidopsis thaliana* chromosome 2 and proteome, the test input of `module2/README.md`. |

## How the pieces run

**`module1.py`** (order from the code; sentinel files `.clean.done`, `.liftover.done`, `.jcvi.done` make it resumable):
1. `module1/fasta_renamer_diploid.py` copies each genome into `--dir_name` with headers renamed `<prefix>_chr<N>` / `<prefix>_sca<N>`.
2. `module1/liftover.py --cdhit --TEsorter --outputs pep bed` (in `--dir_name`): cd-hit and TE-protein removal on `--protein_fa`, then miniprot on every genome.
3. `module1/jcvi_diploid.py --prot --cscore 0.99 --keep pdf anchors`: `jcvi.compara.catalog ortholog` for every genome pair; dot plots moved to `dotplots/`.
4. Per pair, in parallel: `anchor_builder.py` -> `gene_coords_extractor_all4.py -mcscan` -> `anchor_coord_subtracter.py` (twice) -> `anchor_coord_consolidator.py -t 15000 --stitch-gaps`. Final file: `<A>.<B>.anchors.coords.polished2.consolidated`.
5. Intermediates (`.anchors.clean`, `.coords`, `.polished`, `.polished2`) are deleted.

Example (`module1/memo`): `python legacy/module1.py --script_dir legacy/module1 --genomes 9311v2.fa DHX2.fa Nipponb.fa --threads 100 --protein_fa Sbico.pep --dir_name module1`

Not called by `module1.py`: `cds_walker.py` (called by `liftover.py` only when the `inframe` output is requested), `riparian.py`, and `module1/util/`.

**`module1/util/`** (run by hand, order from `module1/util/memo`): `module1.py` first pass -> `homolog_identifier.py` per consolidated file -> `best_homolog.py <REF>` -> concatenate -> `homolog_renamer.py --allow-sca` (writes `*_mod.fa`) -> rerun `module1.py` on the renamed genomes -> `jcvi_orthogroups.py` -> `jcvi_OGseqs.py` -> `backtranslate.sh` -> `build_smtx.py` -> `iqtree3`.

**`module2.py`** (order from the code; clones `cwb14/LTR_HARVEST_parallel`, `cwb14/Kmer2LTR`, `cwb14/TEsorter` into the working directory):
1. (1a) `LTR_HARVEST_parallel/LTR_HARVEST_parallel3.py` (LTRharvest + LTR_FINDER; `--gene-mask --protein` added when `-p` is given) -> `<genome>.harvest.combine.scn`.
2. (1b, only with `--run-optional`) `module2/minimap2_LTR.py` -> `module2/paf_to_scn.py` -> concatenated with the harvest SCN -> `<genome>_merged.scn`.
3. (2) `module2/wfa_LTR` on the merged SCN; rows with `PASS` and LTR/internal lengths >= 100 bp kept -> `.scn.wfa.pass`.
4. (3) `bedtools getfasta`, then `Kmer2LTR/Kmer2LTR.py -D <len> --max-win-overdisp 6`.
5. (4) `module2/add_tsd_to_div.py` -> `module2/kmer2ltr_fa.py` -> `TEsorter -db rexdb-plant -cov 10 -eval 1e-2 -rule 70-30-80` -> awk join -> `<genome>_Kmer2LTR_TSD_class`.
6. (5) `module2/LTR_overlap_purger.py --ov-longest` -> `<genome>_Kmer2LTR_TSD_class_purge` and `.fa`.

No example command in the tree. From its argparse: `python legacy/module2.py -g genome.fa [-p proteins.fa] -t 10 -s legacy/module2 [--run-optional]`. `tandem_cleaner.py` is run by hand afterwards (example in its header).

**`module2/ltrharvest_wrapper2.sh`** (order from the code):
1. Per genome, detection rounds (up to 10): `ltrharvest5.py` on the round's genome -> stop if fewer than `--terminate_count` elements -> `mask_ltr.py` masks the original genome for the next round (found elements become one IUPAC letter per round: N, R, D, Y, S, W, K, M, B, H; everything farther than 15 kb from an element becomes V). Search distances grow each round.
2. `reconcile_nests.py` pools the rounds into `<prefix>_depth<N>_ltr.{tsv,fa}` (N = layers of LTR-RT nested inside).
3. Pooled over all genomes: depth FASTAs concatenated and stripped of non-ACGT -> `Kmer2LTR.py --ltr-cluster --internal-cluster --min-seq-id 0.75` -> `flag_fp_families.py` (both from the Kmer2LTR clone) purges false-positive families -> `<prefix>_depth<N>_clean_ltr.{tsv,fa}`.
4. If the false-positive fraction exceeds `--fp-mask-threshold` (default 0.10) the genome is hard-masked and steps 1-3 rerun in a fresh staging directory.
5. Per genome: `ltr_annotate.py` (adds `strand`, `family` columns) -> `ltr_tsv_to_gff3.py`.
6. Unless `--no-plots`: `ltrharvest_plots.sh` -> `ltrharvest_plot_struct.py`, `ltrharvest_plot.py`, `TEGV.py`.

Example (`module2/README.md`, run from the repository root):
```bash
mkdir athal_test && cd athal_test
bash ../legacy/module2/ltrharvest_wrapper2.sh --genome ../legacy/test/Athal_tair10_chr2.fa.gz \
  --proteins ../legacy/test/Athal.pep.gz --threads 20 --max-rounds 1
```

**`module2/ltrharvest_wrapper.sh`** (older wrapper): the same round loop with `ltrharvest4.py` -> `mask_ltr.py`, and nothing after the loop (no reconciliation, false-positive stage, annotation or plots). Example (header of `ltrharvest4.py`): `bash legacy/module2/ltrharvest_wrapper.sh --genome burnin.fasta --proteins ../PrinTE/data/TAIR10.pep.fa.gz --threads 20 --out_prefix burnin_ltr`

**`module2/solo_LTR/`**, two independent routes:
- `solo_ltr_pipeline.sh`: intact-element headers -> BED -> `bedtools maskfasta` -> `mmseqs easy-search` (or minimap2 + `paf2m8.py`) of LTR consensus vs masked genome, and of the intact library vs masked genome -> `filter_solo_ltr.py`. Example (`solo_LTR/memo`): `bash legacy/module2/solo_LTR/solo_ltr_pipeline.sh -g gen5400000_final.fasta -i gen5400000_final.fasta_r1_ltr.fa -q Kmer2LTR_run.consensus.fa -o solo_ltr_test_minimap -p solo_ltr_test_minimap2 -a minimap22 -f 10 -M 50 -C 0.9 -T 5 -W 1 -X 0 -t 132 -v` (`-a minimap22` is as written in the memo; the script's values are `mmseqs` and `minimap2`)
- `solo_ltr.py`: runs `solo_ltr_search.py` (mask, `makeblastdb`, `blastn`) and filters with `solo_ltr_core.py`. Example (`solo_LTR/memo`): `python legacy/module2/solo_LTR/solo_ltr.py --consensus gen5400000_Kmer2LTR.consensus.fa --internal gen5400000_Kmer2LTR.internal.fa --genome ../gen5400000_final.fasta --outdir gen5400000_work --threads 250 --truth ../gen5400000_final.bed`
- `benchmark.py` is standalone; neither route calls it.

**`module3/modeling_divK_dev/`**: `get_training_data_parallel5.sh` runs `get_divK_Kmer2LTR_benchmarking.pl -uc N -sc N` on every `gen*_final_LTR_Kmer2LTR` file over a random sample of (unique, shared) counts. `test_piecewise_regression6.sh` runs `extract_intact_LTR.py` -> `Kmer2LTR/Kmer2LTR2.py` -> `get_divK_Kmer2LTR.pl` -> `scatterplot_Kmer2LTR2.R` per simulated genome; none of those four scripts is in this tree. The PrinTE simulation commands are in `module3/modeling_divK_dev/memo`.

## Catalog

`current` = latest of its line and called by a wrapper. Only rows marked `ported` have a counterpart in the package. "(clones)" = the script runs `git clone` at run time.

### Top level

| File | What it does | External tools | Status |
|---|---|---|---|
| `module1.py` | Genome FASTAs + proteome -> runs the module1 steps above -> consolidated syntenic block files per genome pair | the module1 scripts | current |
| `module2.py` | Genome (+ proteome) -> runs the module2 steps above -> classified, purged LTR-RT table and FASTA | git, bedtools, awk; clones LTR_HARVEST_parallel, Kmer2LTR, TEsorter | current |

### module1/

| File | What it does | External tools | Status |
|---|---|---|---|
| `fasta_renamer_diploid.py` | Genome FASTAs -> parse chromosome labels from headers -> copies with headers `<prefix>_chr<N>[letter]` / `<prefix>_sca<N>` and a header mapping | bioawk | current |
| `liftover.py` | Genomes + reference proteome -> optional cd-hit, optional TE-protein removal (TEsorter, then blastp against the TEsorter hits), miniprot -> per genome `.gff`, `.pep`, `.bed` (BED6, 0-based half-open), `.mRNA`, `.cds`, `.cds.inframe` | miniprot (clones, builds), samtools, bedtools, cd-hit, TEsorter, seqkit, makeblastdb, blastp | current |
| `cds_walker.py` | miniprot CDS + protein -> removes 1-2 bp (cautiously 3-4 bp) frameshift insertions -> CDS that translates to the protein | Biopython | current (called by `liftover.py`) |
| `jcvi_diploid.py` | `*.bed` + `*.pep` (or CDS) in the working directory -> `python -m jcvi.compara.catalog ortholog` for every pair, with retries -> `<A>.<B>.anchors`, dot plot PDFs | JCVI, DIAMOND (`--prot`) or LAST/BLAST | current |
| `anchor_builder.py` | JCVI anchors file -> keeps the two gene columns, collapses repeated `###` block separators -> cleaned anchors on stdout | | current |
| `gene_coords_extractor_all4.py` | Cleaned anchors + the two `<species>.bed` -> for each two consecutive anchors of a block, the span between them in both genomes and their relative orientation -> `seq:start..end<TAB>seq:start..end<TAB>strand` (start/end copied from the BED) | | current |
| `anchor_coord_subtracter.py` | Block TSV -> within each (ref, query, strand) group removes contained blocks and trims partial overlaps -> block TSV | | current |
| `anchor_coord_consolidator.py` | Block TSV -> merges overlapping or touching blocks when one side is shorter than `-t`; `--stitch-gaps` adds synthetic blocks across gaps (length = end - start) -> block TSV | | current |
| `riparian.py` | Block coords + `.fai` files -> stacked ribbon (riparian) plot, ribbons colored by the top genome's chromosome -> PDF, PNG, HTML | numpy, matplotlib | not in the package |
| `memo` | Example `module1.py` command and five to-do notes | | notes |

### module1/util/

| File | What it does | External tools | Status |
|---|---|---|---|
| `homolog_identifier.py` | Block TSV + `[accession].fa` -> union of block intervals per chromosome pair -> 10-column table of syntenic coverage and strand fractions (half-open lengths) | samtools | not in the package |
| `best_homolog.py` | Coverage table + reference accession -> best reference chromosome and strand per query chromosome | | not in the package |
| `homolog_renamer.py` | Best-homolog table + genome directory -> renames and reverse-complements query chromosomes to match the reference -> `*_mod.fa` | | not in the package |
| `jcvi_orthogroups.py` | `*.anchors` in the working directory -> fully reciprocal orthogroups across all species; optional filters by copy count per species (`-f`) and by syntenic block (`--beds`, `--blocks`) -> orthogroup TSV | | not in the package |
| `jcvi_OGseqs.py` | Orthogroups + pep and CDS FASTAs -> one pep and one CDS FASTA per orthogroup (random pick among duplicates, `--seed`) | | not in the package |
| `backtranslate.sh` | Orthogroup pep/CDS FASTAs -> MAFFT L-INS-i protein alignment -> trimAl `-backtrans` codon alignment | mafft, trimal | not in the package |
| `build_smtx.py` | Codon alignments -> concatenated supermatrix FASTA + IQ-TREE partition files (per gene, per codon) | | not in the package |
| `memo` | Commands for chromosome renaming/reorientation and for the orthogroup species tree | (subphaser, iqtree3) | notes |

### module2/

| File | What it does | External tools | Status |
|---|---|---|---|
| `README.md` | User guide for `ltrharvest_wrapper2.sh`: install, test run, the 21-column `depth<N>_clean_ltr.tsv`, mask letters, multi-genome runs, GFF3, plots, one benchmark table | | notes |
| `todo.md` | One line: merge LTRkin, LTRmix, ltrharvest_phasing10.py and LTRpav into one repo | | notes |
| `ltrharvest_env.yml` | Conda env `ltrharvest4` for the wrapper2 pipeline | | current |
| `ltrharvest_wrapper2.sh` | Genome(s) (+ proteome) -> steps above -> depth tables/FASTAs, GFF3, plots | python, git; clones Kmer2LTR | current |
| `ltrharvest_wrapper.sh` | Genome (+ proteome) -> round loop of `ltrharvest4.py` + `mask_ltr.py` -> per-round `<prefix>_r<N>_ltr.fa` | python | superseded by `ltrharvest_wrapper2.sh` |
| `ltrharvest5.py` | Genome (+ proteome) -> mask genes, LTRharvest + LTR_FINDER on chunks, TSD/protein/homology evidence, Kmer2LTR, duplicate purge -> `<prefix>_ltr.{tsv,fa}` (names `chrom:start-end#class`, 1-based closed) | gt, ltr_finder, bedtools, mafft, trimal; clones Kmer2LTR, minimap2, miniprot, TEsorter2, TRF-mod, sdust | current |
| `ltrharvest4.py` | As above, earlier version (see lineage); called by `ltrharvest_wrapper.sh` | gt, ltr_finder, bedtools, mafft, trimal; clones Kmer2LTR, minimap2, miniprot, TEsorter, TRF-mod | superseded by `ltrharvest5.py` |
| `ltrharvest3.py` | Earlier version; see lineage | as v4, plus Rscript for `--tesorter-tree` | superseded by `ltrharvest5.py` |
| `ltrharvest2.py` | Earlier version; see lineage | as v3 without the TRF-mod clone | superseded by `ltrharvest5.py` |
| `ltrharvest.py` | First version; see lineage | as v2 | superseded by `ltrharvest5.py` |
| `mask_ltr.py` | Genome + element FASTA/TSV (`contig:start-end`, 1-based closed) -> elements masked with one letter, sequence farther than `--distance` masked with another, flanks kept; non-ACGT already in an element is preserved -> masked FASTA on stdout | | current |
| `reconcile_nests.py` | Per-round TSV/FASTA/SCN -> containment tree across rounds -> `<prefix>_depth<N>_ltr.{tsv,fa}`; imports helpers from `ltrharvest5.py` | | current |
| `ltr_annotate.py` | Depth tables -> strand (TEsorter2 call, else domain order, else pass-2 homology target) and family (Kmer2LTR/mmseqs consensus cluster, numbered by size) -> tables rewritten in place | | current |
| `ltr_tsv_to_gff3.py` | Annotated depth tables (+ miniprot GFF) -> `<prefix>_all_depth_LTR_cleaned.gff3` and `..._protein_LTR_cleaned.gff3` (1-based closed); imports from `ltr_annotate.py` | | current |
| `ltrharvest_plots.sh` | Wrapper outputs + genome -> runs the three plotters below into `<prefix>_plots/` | python | current |
| `ltrharvest_plot_struct.py` | Depth TSV (`domains` column) -> per-clade average and individual structure plots, optional per-chromosome HTML browser | matplotlib; samtools with `--genome` | current |
| `ltrharvest_plot.py` | `.fai` + element tables -> multi-page PDF: K2P density, chromosome distribution, element / LTR / internal size distributions, optional family pages from an mmseqs cluster TSV | matplotlib | current |
| `TEGV.py` | Gene GFF3 + per-level TE FASTAs, divergence and domain tables (+ optional BAMs) -> self-contained interactive HTML browser of nested TEs and genes | pysam with `--bam` | current |
| `ltrharvest_plot_struct.R` | LTR-RT TSV + SCN + TEsorter domain GFF -> structure plots (the `.py` version needs neither SCN nor GFF) | R optparse | not in the package |
| `ltr_cluster_phylo.py` | Kmer2LTR consensus LTRs + mmseqs cluster TSV + depth tables -> per-cluster alignment, trimming, ML tree with element-structure glyphs at the tips; imports from `ltrharvest_plot_struct.py` | mafft, trimal, iqtree | not in the package |
| `consensus_library.sh` | Intact TE FASTA -> `mmseqs easy-cluster` -> MAFFT per cluster (optional trimAl) -> consensus by `hmmbuild`+`hmmemit` or EMBOSS `consambig` (IUPAC) | mmseqs, mafft, trimal, HMMER, EMBOSS, seqkit (optional) | not in the package |
| `minimap2_LTR.py` | Genome -> sliding 30-kb windows, low-complexity masking, minimap2 self-alignment of each window -> `.paf`, `.filtered.paf` (LTR pair candidates), `.filtered.fa`, genome without LTR-RTs | clones minimap2, miniprot, longdust, TRF-mod, barrnap, sdust | current |
| `paf_to_scn.py` | Filtered PAF -> same-sequence pairs, reciprocal duplicates dropped -> SCN rows (PAF start/end values copied unchanged) | | current |
| `wfa_LTR` | Compiled binary of `wfa_ltr.cpp` (x86-64 Linux ELF, dynamically linked) | | current |
| `wfa_ltr.cpp` | SCN (1-based closed) + genome -> WFA alignment of the two LTR ends of each candidate, boundary trimming, PASS/fail and TSD columns | WFA2-lib (not in the tree) | current |
| `add_tsd_to_div.py` | Kmer2LTR table + WFA pass table -> appends TSD and motif by `chrom:start-end` key | | current |
| `kmer2ltr_fa.py` | Genome + table with `chrom:start-end` in column 1 (read as 1-based closed) -> FASTA on stdout | | current |
| `LTR_overlap_purger.py` | Classified table -> resolves shared-boundary and overlapping candidates (keep longest, lowest divergence, drop all or retain all) -> purged table | | current |
| `tandem_cleaner.py` | `module2.py` FASTA + table -> TRF-mod on each element -> drops candidates above a tandem-repeat percentage in the LTRs (`-lp`) or whole element (`-wp`) | clones TRF-mod | not in the package |
| `ltrharvest_phasing10.py` | LTR-RT FASTA + annotation (GFF3 or K2P table) + optional `.fai` -> subgenomes, homoeolog pairs, split and merger ages, exchanges -> tables, `findings.txt`, `report.pdf` (v10.1.0) | numpy, scipy, numba, edlib, matplotlib | ported: `synltr subgenomes` |
| `ltrharvest_phasing9.py` | v9.1: v8 plus homoeolog pairing and what the pairs show; see lineage | numpy, scipy, numba, edlib, matplotlib | superseded by `ltrharvest_phasing10.py` |
| `ltrharvest_phasing8.py` | v8: v7 with new merger dating; see lineage | numpy, scipy, numba, edlib, matplotlib | superseded by `ltrharvest_phasing10.py` |
| `ltrharvest_phasing7.py` | v7: number of subgenomes searched by insertion-age class; see lineage | numpy, scipy, numba, edlib, matplotlib | superseded by `ltrharvest_phasing10.py` |
| `ltrharvest_phasing6.py` | v6: LTR-RT FASTA + annotation + optional `.fai` -> phasing, structure, dating, ancestry -> tables and `report.md` | numpy, scipy, numba, edlib, matplotlib | superseded by `ltrharvest_phasing10.py` |
| `ltrharvest_phasing5.py` | v5 (5.0): LTR-RT FASTA (+ K2P table, + genome) -> phasing, per-copy ancestry, exchange segments, ages, family tables | numpy, scipy, matplotlib | superseded by `ltrharvest_phasing10.py` |
| `ltrharvest_phasing4.py` | v4 (4.2): first config-free version; LTR-RT FASTA (+ K2P table, + genome) -> subgenome assignment, per-copy lineage calls, merger date, exchange candidates | numpy, scipy, matplotlib | superseded by `ltrharvest_phasing10.py` |
| `ltrharvest_phasing3.py` | LTR library + homoeolog config -> jellyfish k-mer counts per chromosome, differential k-mers, KMeans, bootstrap -> subgenome-specific k-mers, heatmap, PCA | jellyfish, pandas, scipy, scikit-learn | superseded by `ltrharvest_phasing10.py` |
| `ltrharvest_phasing2.py` | LTR library + homoeolog config -> SubPhaser-like phasing, K2P pre/post-polyploid densities, 1-Mb window exchange scan (Fisher, BH), chromosome painting | jellyfish, mmseqs, scipy, scikit-learn | superseded by `ltrharvest_phasing10.py` |
| `ltrharvest_phasing.py` | LTR library + homoeolog config -> SubPhaser-style phasing, K2P annotation, window enrichment scan, chromosome plot | jellyfish, scipy, scikit-learn | superseded by `ltrharvest_phasing10.py` |
| `simallo.py` | Scenario JSON -> forward simulation of LTR-RT families through lineage splits, allopolyploid mergers and homoeologous exchanges -> LTR-RT FASTA, GFF3, K2P table, `.fai` and truth tables | numpy | ported: `synltr sim nucleus` |
| `ltrkin.tar.gz` | ltrkin 0.2.1: source (`src/ltrkin`, 9 files), 7 test files, `benchmark/` (40 scripts and workflow files), `docs/` | | ported: `synltr genotype` (reads), `synltr source` |
| `LTRmix.tar.gz` | LTRmix: `ltrmix.py`, `ltrmix_plot.py`, 3 test files, `benchmark/` (358 files, 29 MB: reports, scripts, result tables, figures), `docs/` | | ported: `synltr geneflow --screen` |
| `LTRpav.tar.gz` | LTRpav 0.1.0: source (`src/ltrpav`, 56 files), 42 test files, `docs/`, `examples/quickstart.sh`, CI workflow, development records | | ported: `synltr genotype`, `synltr geneflow`, `synltr sim pop`, `synltr sim score` |

Notable documents inside the tarballs:
- ltrkin: `docs/DESIGN.md` (spec with decisions and reasons), `docs/lit_review.md` (long-read TE genotyping literature and software), `docs/benchmark_v0.2.md` (haplotype-resolved assemblies as truth), `docs/benchmark.md` and `docs/benchmark_genotype_sim.md` (v0.1, simulated reads vs Sniffles2), `docs/outputs.md` (every output column), `benchmark/README.md` (how to rerun; Snakefile in `benchmark/assembly_truth/`).
- LTRmix: `DESIGN.md` (what LTR-RTs give that SNPs do not; test case maize/*mexicana*/*parviglumis*), `PLAN.md` (method plan), `PROMPT.md` (the original brief), `benchmark/REPORT.md` (full benchmark), `benchmark/dead_ends.md` (rejected approaches and corrected mistakes), `benchmark/memo_ltrmix.md` (run record), `benchmark/figures/SPEC.md` (figure spec), `docs/idea_synteny_input.md` (unbuilt idea: optional synteny input for distant donors). `benchmark/README.md` states the scripts cannot be rerun (inputs not kept, absolute paths).
- LTRpav: `docs/METHODS.md` (genotyping and analyses in full), `docs/plans/` (three implementation plans: genotyping engine, analyses, benchmark).

**Version lineage, `ltrharvest*.py`** (from each file's header):
- `ltrharvest.py`: mask genes and non-LTR TEs; LTRharvest + LTR_FINDER merged; TEsorter on the internal sequence as filter; Kmer2LTR; purge duplicates by LTR divergence. Nesting by manual rounds of `mask_ltr.py` + `--require-run-chars`. Notes `--tesorter-tree` as broken.
- `ltrharvest2.py`: no version note in file. From a diff of the headers: example commands gain `--tsd-rescue`, a to-do on masking low-complexity sequence with TRF-mod, and a third nesting round.
- `ltrharvest3.py`: adds TRF-mod hard-masking of high-copy tandem repeats (`--trf`, `--trf-min-copy`); switches to a custom TEsorter whose second pass uses mmseqs2 instead of blastn and accepts custom libraries.
- `ltrharvest4.py`: TSD search runs before TEsorter and TSD-bearing candidates feed TEsorter's second pass; candidates kept if they have a TSD, an LTR-RT protein, or homology to one that does. Run through `ltrharvest_wrapper.sh`.
- `ltrharvest5.py`: resolves several layers of nesting within one round; header lists two open issues (nested elements classified as one unit by TEsorter; reconciling rounds with single-round nesting). From the code, not a header note: clones TEsorter2 and sdust.

**Version lineage, wrappers** (from the header comments): `ltrharvest_wrapper2.sh` adds to `ltrharvest_wrapper.sh` several genomes with one shared family namespace, `ltrharvest5.py` in place of `ltrharvest4.py`, `reconcile_nests.py` depth buckets, `--run-sdust`, false-positive family correction with automatic rerun, the annotation stage and the plotting stage.

**Version lineage, `ltrharvest_phasing*.py`** (from each file's docstring and closing notes):
- v1 (`ltrharvest_phasing.py`), v2: SubPhaser-style, need a homoeolog config; header says the two "use slightly different approaches and give slightly different results" and both are kept. v2 adds K2P pre/post-polyploid densities and the window exchange scan with chromosome painting.
- v3: "a simpler version", described as possibly more robust on harder data (strawberry); offers a paired t-test, Kruskal-Wallis or Wilcoxon test, and its header says the earlier two scripts would need updating for the latter two.
- v4 (4.2): config-free. Rare k-mers shared by few copies link chromosomes; degree-corrected null; number of subgenomes inferred by split-half replication; per-copy lineage calls; merger dating from K2P; HMM exchange track; `--genome` window painting.
- v5 (5.0): no change note in file. Docstring adds nearest-relative votes to the links, a reported structure tree, per-copy ancestry from two calibrated scores, family and cross-subgenome transposition tables.
- v6: no change note in file. Docstring: exact hashing of canonical 21-mers carried by 2-250 copies; input is the LTR-RT FASTA plus an annotation (GFF3 or K2P table) and optional `.fai`, no genome sequence.
- v7: number of subgenomes searched on all copies and in each LTR-divergence tercile; first file with a "NOTES FOR FUTURE WORK" block at the end.
- v8: merger dating by a genealogy clock on sibling pairs; split-dating order changed; one-page `report.pdf`. Its header names `ltrharvest_phasing8_mod.py`, which is not in the tree.
- v9 / v9.1: homoeolog pairing from inherited insertions; exchanges between homoeologs, intact LTR-RT half-life per subgenome, translocations, reverse-collinear homoeologs; `findings.txt`.
- v10 / v10.1: core-first phasing (unplaced scaffolds placed afterwards), age rule applied to the top split, inherited-insertion split dating, merger-dating cap on one clock; LTRs recovered for records with nested insertions cut out; homoeolog check on age-class sub-splits. Closing notes (line 6299 on) hold benchmark results, dead ends and leads for all versions.

### module2/solo_LTR/

| File | What it does | External tools | Status |
|---|---|---|---|
| `solo_ltr_pipeline.sh` | Genome + intact LTR-RT FASTA (headers 1-based closed, converted to BED) + LTR consensus -> hard-mask intact elements, search consensus and intact library against the masked genome, filter -> solo-LTR loci BED/TSV | mmseqs, bedtools, samtools, awk, python3; minimap2 fork `cwb14/minimap2` (clones, builds) with `-a minimap2` | current |
| `filter_solo_ltr.py` | Search hits (mmseqs format-mode 4) -> keep hits spanning the consensus end to end within `FLANK` bp, purge fragments that extend into internal sequence, require a TSD -> `.loci.bed` and reports | pandas, pyfaidx | current |
| `paf2m8.py` | minimap2 PAF (0-based half-open) -> the mmseqs table above (1-based closed, minus-strand query coordinates reversed) | | current |
| `solo_ltr.py` | LTR consensus + internal FASTA + genome -> mask, BLAST, TSD-gated filter -> `solo_ltr.tsv`, `nested_solo.tsv` (solo LTRs inside other LTR-RTs); `--truth` reports precision/recall/F1 | via `solo_ltr_search.py` | current |
| `solo_ltr_search.py` | Cached mask / `makeblastdb` / chunked `blastn` steps for `solo_ltr.py`; writes header `start`, `end` unchanged into the mask BED | bedtools, makeblastdb, blastn | current |
| `solo_ltr_core.py` | Shared filters for `solo_ltr.py`: HSP cutoffs, TSD at HSP boundaries, locus merging, internal-proximity filter, divergence (0-based half-open throughout) | numpy | current |
| `benchmark.py` | Predicted loci + PrinTE truth BED -> truth classes (clean solo, intact, fragment, disrupted solo) -> precision/recall/F1 by reciprocal overlap | | not in the package |
| `memo` | Commands, parameter comparison table for both routes on one simulated genome, limits noted by the author | (Kmer2LTR, PrinTE) | notes |

### module3/modeling_divK_dev/

| File | What it does | External tools | Status |
|---|---|---|---|
| `get_divK_Kmer2LTR_benchmarking.pl` | Kmer2LTR table (column 1 ID, `shared_` prefix = full site; column 11 divergence) -> piecewise regression scan for the optimum K with bootstraps; `-sc`/`-uc` resample shared and unique counts | Perl `Distributions.pm` (not in the tree) | not in the package |
| `get_training_data_parallel5.sh` | `gen*_final_LTR_Kmer2LTR` files -> random (unique, shared) count pairs -> runs the Perl script in parallel (output layout: not determined) | perl | not in the package |
| `test_piecewise_regression6.sh` | Simulated BED/FASTA pairs -> extract intact LTR-RTs, Kmer2LTR, divK, scatter plot against truth | python, perl, Rscript; four scripts not in the tree | not in the package |
| `memo` | Aim (train a model to correct divK), PrinTE simulation commands, the `$maxdist = 0.1` saturation caveat | (PrinTE) | notes |

### test/

| File | What it does | External tools | Status |
|---|---|---|---|
| `Athal_tair10_chr2.fa.gz` | *A. thaliana* TAIR10 chromosome 2 | | data |
| `Athal.pep.gz` | *A. thaliana* proteome | | data |

## Approaches that exist only here

- **Gene-synteny anchoring with JCVI.** `module1.py`, `module1/liftover.py`, `module1/jcvi_diploid.py`: one trusted proteome is lifted onto every genome with miniprot so all genomes share gene IDs, then JCVI finds pairwise syntenic anchors.
- **Anchors to clean block coordinates.** `module1/gene_coords_extractor_all4.py`, `anchor_coord_subtracter.py`, `anchor_coord_consolidator.py`: consecutive anchors of a block become coordinate intervals with orientation; contained intervals are removed, overlaps trimmed, short blocks merged and gaps stitched.
- **TE-protein removal from a proteome.** `module1/liftover.py` (`run_tesorter_two_pass`): TEsorter flags TE proteins, then blastp against those hits removes a second set.
- **In-frame CDS repair.** `module1/cds_walker.py`: walks the miniprot CDS against its protein and removes frameshift insertions.
- **Synteny-guided chromosome naming and orientation.** `module1/util/homolog_identifier.py`, `best_homolog.py`, `homolog_renamer.py`: syntenic coverage picks each chromosome's best reference homolog and strand; chromosomes are renamed and flipped.
- **Species tree from syntenic orthogroups.** `module1/util/jcvi_orthogroups.py`, `jcvi_OGseqs.py`, `backtranslate.sh`, `build_smtx.py`: fully reciprocal orthogroups from pairwise anchors, filtered by ploidy-aware copy counts and block membership, to a codon supermatrix.
- **Riparian synteny plot.** `module1/riparian.py`.
- **De novo LTR-RT discovery with nesting depth.** `module2/ltrharvest_wrapper2.sh`, `ltrharvest5.py`, `mask_ltr.py`, `reconcile_nests.py`: each round masks found elements with a round-specific IUPAC letter so the next round finds elements that contain them; containment across rounds gives each element a depth. The package reads this output (its README lists `<prefix>_depth*_clean_ltr.tsv` as an annotation input) but does not produce it.
- **False-positive family correction with rerun.** `module2/ltrharvest_wrapper2.sh` (`run_merged_stage`, `run_fp_stage`, `run_fp_orchestrator`): families flagged by Kmer2LTR's `flag_fp_families.py` are purged; above a threshold the repeats are masked and detection reruns.
- **Shared family names across genomes.** `module2/ltrharvest_wrapper2.sh`, `ltr_annotate.py`: elements of all genomes are clustered once so a family label means the same in every genome.
- **Strand and GFF3 for LTR-RTs.** `module2/ltr_annotate.py` (strand from TEsorter2's call, else protein-domain order, else the pass-2 homology target, else unknown; orienting by family members was tried and removed, see its header), `ltr_tsv_to_gff3.py`.
- **WFA check of LTR boundaries.** `module2/wfa_ltr.cpp`, called from `module2.py`: aligns the ends of each candidate's two LTRs, trims over-extended boundaries and reports TSDs. `module2.py` holds the parameter ranges from a grid search (lines 307-317).
- **LTR-RT candidates from minimap2 self-alignment.** `module2/minimap2_LTR.py`: windows of the genome are aligned to themselves and same-strand repeat pairs at LTR-RT spacing become candidates.
- **Overlap and tandem filters.** `module2/LTR_overlap_purger.py`, `module2/tandem_cleaner.py`.
- **Solo-LTR search.** `module2/solo_LTR/solo_ltr_pipeline.sh` (mask intact elements, keep consensus hits that span the whole LTR, purge fragments, require a TSD) and `module2/solo_LTR/solo_ltr.py` (BLAST route; also solo LTRs nested in other elements). `solo_ltr.py`'s header lists untested ideas for both. The package genotypes the solo-LTR state only at annotated full-length elements.
- **divK modeling.** `module3/modeling_divK_dev/`: piecewise regression for the optimum K and a training-data generator on PrinTE simulations; the memo records why the first training set saturates.
- **Consensus library.** `module2/consensus_library.sh`: cluster, align, consensus (HMM or IUPAC).
- **Per-cluster LTR phylogenies with structure glyphs.** `module2/ltr_cluster_phylo.py`.
- **LTR-RT plots and browser.** `module2/ltrharvest_plot_struct.py`, `ltrharvest_plot.py`, `TEGV.py`.
- **Config-based subgenome phasing.** `module2/ltrharvest_phasing.py`, `ltrharvest_phasing2.py`, `ltrharvest_phasing3.py`: differential k-mers between homoeologous chromosome sets given in a config, counted with jellyfish, clustered with KMeans.
- **Genome painting from the assembly.** `module2/ltrharvest_phasing4.py`, `ltrharvest_phasing5.py` (`--genome`): 50-kb windows painted with subgenome-confined k-mers, including unplaced scaffolds. Later versions take no genome sequence; the v10 closing notes record a dropped genome/RepeatMasker mode.
- **Developer options of v10.** `module2/ltrharvest_phasing10.py`: hidden `--dev_*` options, `--label`, and non-default settings that the package CHANGELOG lists as removed.
- **In the tarballs, removed from the package (package CHANGELOG):** LTRpav anchoring routes `mm2`, `fastga`, `flank-blast`, `hybrid` (`src/ltrpav/anchor/`), its direction test (`src/ltrpav/analysis/direction.py`) and single-stage commands; LTRmix's `--mu` age axis in `ltrmix_plot.py`.

## Running legacy code

**Environments.** `module2/ltrharvest_env.yml` (env `ltrharvest4`) is the only env file; it covers `ltrharvest_wrapper2.sh`. `module2/solo_LTR/memo` uses an env named `synLTR` and `module1/riparian.py` an env named `synteny`; neither has a file here. Tools used by scripts and not in the env file: JCVI, DIAMOND, bioawk, cd-hit, seqkit, samtools, the `TEsorter` command (module1); samtools, pandas (solo_LTR); iqtree3, subphaser (`module1/util/memo`); pandas, scipy, scikit-learn, numba, edlib, jellyfish (phasing); iqtree, EMBOSS (`ltr_cluster_phylo.py`, `consensus_library.sh`); pysam (`TEGV.py --bam`); Perl, R (module3, `ltrharvest_plot_struct.R`).

**Commands** (thread counts are the memos'; set your own):
```bash
# module1 (module1/memo)
python legacy/module1.py --script_dir legacy/module1 --genomes 9311v2.fa DHX2.fa Nipponb.fa --threads 100 --protein_fa Sbico.pep --dir_name module1

# module1/util, inside the module1 output directory (module1/util/memo; paths adjusted)
for f in *.anchors.coords.polished2.consolidated; do python ../legacy/module1/util/homolog_identifier.py --input "$f" > "${f}.pairs"; done
python ../legacy/module1/util/homolog_renamer.py --allow-sca all.anchors.coords.polished2.consolidated.pairs.best ./
python ../legacy/module1/util/jcvi_orthogroups.py -f ploidy.tsv --beds ./*bed --blocks "*anchors.coords.polished2.consolidated" > orthogroups.tsv

# module2 discovery with nesting (module2/README.md)
mamba env create -f legacy/module2/ltrharvest_env.yml && mamba activate ltrharvest4
mkdir athal_test && cd athal_test
bash ../legacy/module2/ltrharvest_wrapper2.sh --genome ../legacy/test/Athal_tair10_chr2.fa.gz --proteins ../legacy/test/Athal.pep.gz --threads 20 --max-rounds 1

# subgenome phasing (usage line of ltrharvest_phasing10.py)
python legacy/module2/ltrharvest_phasing10.py --ltr_fasta LTR.fa --annotation ANNOT.gff3 [--fai GENOME.fai] --outdir OUT [-t 8]

# solo LTRs (module2/solo_LTR/memo)
python legacy/module2/solo_LTR/solo_ltr.py --consensus LTRs.alns.consensus.fa --internal LTRs.alns.internal.fa --genome gen5400000_final.fasta --outdir work --threads 48 --truth gen5400000_final.bed --force
```

**Hard-coded paths and machine-specific assumptions.**
- Absolute paths: `module1/riparian.py:44` (`/home/chris/bin/mambaforge/envs/synteny/bin/python`, in an error message; line 45 points to an `environment.yml` beside the script that is not in the tree); `module3/modeling_divK_dev/get_divK_Kmer2LTR_benchmarking.pl:1` (a `cat /home/chris/...` line above the shebang); `module3/modeling_divK_dev/memo:3`; `module2/solo_LTR/memo:52`; provenance-only paths in the closing notes of `module2/ltrharvest_phasing7.py:4332`, `8.py:5011`, `9.py:5868`, `10.py:6301-6302`.
- `git clone` (and `make`) at run time, into the working directory, so network access and a compiler are needed: `module2.py:119-121`; `module1/liftover.py:25-26`; `module2/ltrharvest5.py:1891, 1934, 1956, 2077, 2111, 2142` (and the same in v1-v4); `module2/ltrharvest_wrapper2.sh:429`; `module2/minimap2_LTR.py:113, 140, 168, 199, 236, 331`; `module2/tandem_cleaner.py:51`; `module2/solo_LTR/solo_ltr_pipeline.sh:207`.
- Relative paths that fix the working directory: `module2.py:358` (`./Kmer2LTR/Kmer2LTR.py`); `module3/modeling_divK_dev/get_training_data_parallel5.sh:21`; `module3/modeling_divK_dev/test_piecewise_regression6.sh:22-25`.
- `python` taken from PATH rather than the running interpreter: `module1/jcvi_diploid.py:189`; `module2/ltrharvest_wrapper2.sh:514, 604, 611, 664, 1220, 1309, 1354`; `module2/ltrharvest_wrapper.sh:405, 488`; `module2/ltrharvest_plots.sh:144, 160, 176`. `module2.py:321, 328, 339, 343, 436` run commands through `bash -lc` (a login shell).
- Thread defaults of all cores (`nproc`): `module2/solo_LTR/solo_ltr_pipeline.sh:85`; `module2/consensus_library.sh:67`; `module3/modeling_divK_dev/get_training_data_parallel5.sh:34-35`.
- Dataset-specific defaults: `module2/solo_LTR/solo_ltr_pipeline.sh:45-47` (`gen5400000_final.fasta` and its companions).
- `module2/wfa_LTR` is a prebuilt x86-64 Linux binary. The compile line in `module2/wfa_ltr.cpp:16` needs WFA2-lib and writes `module2/wfa_ltr`, while `module2.py:296` looks for `wfa_LTR`.
- Un-commented example lines that execute or break the script as written: `module2/consensus_library.sh:13, 15, 19` (line 19 calls the script itself, before argument parsing); `module2/ltrharvest_plot_struct.R:4` (a shell command inside the R script); `module3/modeling_divK_dev/get_divK_Kmer2LTR_benchmarking.pl:1`.
- Coordinate hand-offs to check before reuse: `module2/paf_to_scn.py` copies PAF values (0-based half-open) into SCN, which `module2/wfa_ltr.cpp` reads as 1-based closed; `module2.py:342` passes columns 2 and 5 of the `.pass` table unchanged to `bedtools getfasta` (BED) while `module2/kmer2ltr_fa.py` reads the resulting `chrom:start-end` keys as 1-based closed (convention of the `.pass` table: not determined); `module2/solo_LTR/solo_ltr_search.py:48-49, 95` treats `chrom:start-end` headers as 0-based half-open while `module2/solo_LTR/solo_ltr_pipeline.sh:301-322` treats the same header style as 1-based closed.
- `module1/util/memo` passes `*.cds.inframe` files to `jcvi_OGseqs.py`, but `module1.py:112` asks `liftover.py` for `--outputs pep bed` only.
