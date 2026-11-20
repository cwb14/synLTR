# Outputs

Every file each command writes, with its coordinate convention. Coordinates are 1-based closed in the tables, GFF3
and VCF, 0-based half-open in BED and bedGraph files, unless a table below says otherwise. A finished `genotype`,
`geneflow` or `subgenomes` run writes `manifest.json` last, so it marks a finished run: versions, command,
parameters, and the checksums of the input files its stages keyed on (on assemblies each genome's FASTA and
annotation, and each `--paf` file; for the screen the library files; on read sets the reference FASTA; none for
`subgenomes`, whose stages key on the library as read). `work/` holds the cache a rerun resumes from, and a
`.done.json` beside a file marks it finished, so a rerun skips it. GFF3 seqids and attribute values are escaped as
GFF3 requires (a PanSN name `HG002#1#chr1` is written `HG002%231%23chr1`).

**`genotype` on assemblies**

| File | What it holds | Coordinates |
|---|---|---|
| `loci.tsv` | one row per insertion (locus): family, length, TSD, discovering genomes; per genome its state, tier, position, strand, K2P, flags | 1-based closed |
| `pav_matrix.tsv` | loci x genomes, states only | |
| `qc.tsv` | per pair and genome: state fractions, resolved fraction, agreement with the annotation, skipped records | |
| `pav.vcf.gz` (+ `.tbi`) | one site per locus on the first row's genome: GT 1 carries the insertion (P/S/T), 0 the empty site (E), `.` missing; per sample ST (state), TR (tier) | VCF; POS is the base before the insertion |
| `pav.unplaced.tsv` | loci the VCF cannot place, with the reason | |
| `<label>.pav.gff3` | each genome's P, S, T, E and D sites in its own coordinates: state, tier, family, and every genome's state at the locus (`states=`) | 1-based closed |
| `<label>.pav.bed` | the same sites: BED6 plus state, tier and family | 0-based half-open |
| `pairs/<A>__<B>.calls.tsv.gz` | each element of A, its call in B and the evidence | 1-based closed |

**`geneflow`**: `genotype`'s files, plus

| File | What it holds | Coordinates |
|---|---|---|
| `stats.tsv` | every statistic with SE, Z and the loci it counts (`layer test scope stratum statistic value se z n loci note`) | |
| `findings.txt` | the results in plain words, with warnings | |
| `posterior.<R>.tsv`, `decoy_posterior.<S>.tsv` | per insertion of the recipient R (and of its sister S, painted as a decoy): class, posterior of donor-like ancestry, its tract | 1-based closed (the tract string too) |
| `tracts.<R>.bed`, `decoy_tracts.<S>.bed` | donor-like tracts, with q against the decoy | BED 0-based half-open (`start0`) |
| `inversions.tsv` | inverted blocks between genome pairs | 1-based closed |
| `mechanisms.tsv` | what each pair's absences are (empty sites, deletions, replacements) by family and age | |
| `identity.tsv`, `pairs/<A>__<B>.div.tsv.gz` | with `--identity`: the LTR divergence of shared insertions, per locus and per element pair | |
| `overview.png`, `overview.pdf`, `overview_numbers.tsv` | the figure (the PDF adds short and long legends) and every number in it | |

**`geneflow --screen`**

| File | What it holds | Coordinates |
|---|---|---|
| `pairs/<A>__<B>.calls.tsv.gz` | each record of A and its call in B: P (1:1 ortholog), A (inferred absence, never an observed empty site) or U (no collinear block). For P, `b_seqid`, `b_start`, `b_end` are the ortholog (`.` when its name holds no coordinates); for A, the window of B searched for it (between two anchors' midpoints, not a site); `b_strand` the block's sign and `wid` its row among the pair's rows of `blocks.tsv` (0-based); `b_ins` the room between the anchors, not inserted bases; `tsd_a`, `tsd_b` in canonical orientation (the lesser of a TSD and its reverse complement); `tier` and the other evidence columns `.` | 1-based closed |
| `pairs/<A>__<B>.div.tsv.gz` | sites compared and differing for each ortholog pair | |
| `blocks.tsv` | collinear blocks between genome pairs | 1-based (anchor midpoints) |
| `posterior.<R>.tsv`, `decoy_posterior.<S>.tsv` | per record: class and two posteriors (donor-like ancestry, recent shared ancestry) | 1-based closed |
| `tracts.<R>.bed`, `decoy_tracts.<S>.bed` | the tracts of both paintings (column `alphabet`) | BED 0-based half-open |
| `inversions.tsv` | inverted blocks between genome pairs | 1-based closed |
| `stats.tsv`, `findings.txt`, `overview.png`, `overview.pdf`, `overview_numbers.tsv` | as for `geneflow` (`stats.tsv` rows under layer `library`) | |

**`genotype` on read sets, and `source`**

| File | What it holds | Coordinates |
|---|---|---|
| `reads.genotypes.tsv` | one row per element and sample: genotype of F (full), S (solo LTR), E (empty site), N (no sequence) alleles, GQ, `pav` (present, absent, solo, uncertain), reads per allele, flags | 1-based closed |
| `reads.vcf.gz` (+ `.tbi`) | the same calls: REF the element, ALT `<DEL:ME:LTR>` (empty site) and `<DEL:ME:SOLO>` (solo LTR), GT `.` a null haplotype | VCF; POS is the base before the element |
| `reads.unplaced.tsv` | elements the VCF cannot place (one at a sequence's first base), with the reason | |
| `reads.<sample>.evidence.json.gz` | each sample's read evidence, reused by a rerun | |
| `reads.sources.tsv` (`source`) | each new copy with its likely source copy: divergence, ties, age, shared mutations, two-parent test, flags | 1-based closed |
| `reads.copies.tsv` (`source`) | every element: new, fixed or other, LTR divergence, nearest relatives | 1-based closed |
| `reads.locality.tsv` (`source`) | per clade and family: do new copies sit on their source's chromosome, against a null? | |

**`subgenomes`**

| File | What it holds | Coordinates |
|---|---|---|
| `chromosomes.tsv` | each sequence's subgenome, support and homoeolog partner (`core`: decided the structure) | |
| `structure.tsv` | each split tested, and the rule that accepted or rejected it | |
| `dating.tsv` | split and merger ages of each subgenome pair by every estimator (K2P, and years at `--mu`) | |
| `homoeologs.tsv`, `homoeolog_sets.tsv`, `pairing_map.tsv` | homoeolog calls with their evidence, connected sets, and the pairing tests | |
| `homoeolog_anchors.tsv` | every old insertion two sequences share: both copies, whether on the pair's collinear chain | 1-based closed |
| `homoeolog_exchanges.tsv`, `homoeolog_rearrangements.tsv` | near-identical stretches on homoeologs, translocated blocks | 1-based closed |
| `homoeolog_halflife.tsv` | half-lives of intact LTR-RTs, per subgenome pair and direction | |
| `segments.tsv`, `exchanges.gff3`, `exchanges.bed` | painted segments and the called exchanges | 1-based closed; BED 0-based half-open |
| `ancestry_<SG>.bedGraph`, `resolution.bedGraph` | posterior ancestry of each subgenome along each sequence; the smallest detectable exchange | 0-based half-open |
| `elements.tsv` | every LTR-RT: subgenome, clade, family, LTR divergence, age, subgenome posteriors | 1-based closed |
| `clades.tsv`, `insertion_history.tsv`, `nesting.tsv`, `links_oe.tsv` | clades per subgenome, insertions through time, nested pairs, shared-lineage excess between sequences | |
| `findings.txt`, `summary.json` | the results in plain words; every estimate, machine-readable | |
| `report.pdf`, `report.png`, `report_data.json`, `report_numbers.tsv` | the figure, its data and every number in it | |

**`plot`** redraws `overview.*` (`geneflow`; `--trio` names its lead trio) or `report.*` (`subgenomes`) in the run
folder.

**`sim`**

| Command | Files | Coordinates |
|---|---|---|
| `sim pop` | `<label>.fa`, `<label>.gff3` (annotated taxa), `samples.tsv`, `params.json`; `truth/<X>__<Y>.tsv`, `truth/insertions.tsv`, `truth/tracts.<label>.bed`, `truth/fractions.tsv` | GFF3 and `<X>__<Y>.tsv` 1-based closed; `insertions.tsv` `pos` 0-based; BED 0-based half-open |
| `sim nucleus` | `ltr.fa`, `k2p.tsv`, `annotation.gff3`, `genome.fai`; in `--truthdir`: `copies.tsv`, `exchanges.tsv`, `chromosomes.tsv`, `orthologs.tsv`, `paralogs.tsv`, `history.json` | record names and GFF3 1-based closed; `k2p.tsv` LTR offsets 1-based along the copy's genome span from its start, the genome-left LTR first on either strand (as LTRquest writes them) |
| `sim score` | `PREFIX.score.tsv` (precision and recall per state and tier), `PREFIX.confusion.tsv` | |
