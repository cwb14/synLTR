"""Fixed constants, grouped by the module that uses them, each with the reason it is fixed."""

from synltr import __version__

# cli
DEFAULT_THREADS = 4  # total CPU budget when -t is not given; small enough for a shared machine
DEFAULT_PLOIDY = 2  # genotype --ploidy when not given (read sets): chromosome copies per sample, ltrkin's default
DEFAULT_MIN_MAPQ = 20  # genotype --min-mapq when not given: a read placed with MAPQ < 20 (1% error) anchors nothing

# runtime
PR_SET_PDEATHSIG = 1  # <linux/prctl.h>: request a signal when the parent process dies
PARENT_POLL_SECONDS = 1.0  # heuristic: where prctl is missing, how often a worker checks that its parent lives
PAR_CHUNKS = 64  # chunks of a numba parallel accumulation: fixed, so its float sums do not depend on -t
PAR_CHUNK_BYTES = 1_000_000_000  # the chunks' C x C float64 matrices together stay under 1 GB

# checkpoint
CHECKSUM_CACHE = "checksums.json"  # name of the checksum cache file in a work directory
CHECKSUM_BLOCK = 8 << 20  # bytes hashed per read: streams multi-Gb files in constant memory
# heuristic: seconds a tool's version probe may take before the version is recorded as unknown
TOOL_VERSION_TIMEOUT = 10
# executable -> the arguments that print its version (first output line kept)
TOOL_VERSION_ARGS: dict[str, list[str]] = {"lastz": ["--version"], "minimap2": ["--version"], "bcftools": ["--version"]}

# io.seqio
FASTA_BLOCK = 1 << 24  # bytes read at a time when streaming a FASTA: few Python steps per record, little memory

# io.annotation and routes.assembly.blocks
RECORD_WARNINGS = 5  # malformed records warned about one by one per file: a multi-Gb GFF3 or PAF must not flood the log

# io.annotation
TSD_ORDER = (5, 4, 6)  # TSD lengths tried in turn: 5 bp is the modal LTR-RT TSD, then 4 and 6
# heuristic (ltrkin's): an annotated TSD outside 3-8 bp is not a TSD (LTR-RTs make 4-6-bp ones, TSD_DEFAULT)
TSD_MIN_LEN, TSD_MAX_LEN = 3, 8
BED_LTR_EDGE = 20  # heuristic (LTRpav's): a BED element's LTR self-alignment may start or end this far from its ends
BED_LTR_MAX = 6000  # longest LTR searched for in a BED element (LTRs rarely exceed ~5 kb)

# io.library (phasing's rule) and analyses.subgenomes.inputs
# heuristic (phasing's): share of FASTA copies that must be annotated (fail fast otherwise)
MIN_MATCH = 0.95
# io.library and seq.divergence: record-end LTRs of a clean record (nests cut out) are used when their K2P matches the
# annotated LTR K2P as well as 95% of the same genome's copies with genome-mapped LTRs do (z = |difference| / SE; this
# value is the floor: sampling error alone). Per genome because aligners differ: 98-99% of mapped copies pass 1.96 in
# most genomes, 78% in a woody bamboo
CLEAN_Z = 1.96
CLEAN_CAL = 2000  # heuristic (phasing's): mapped copies (random, seed 0) that calibrate that tolerance
# heuristic (phasing's): fewer mapped copies leave the tolerance at its floor (CLEAN_Z)
CLEAN_CAL_MIN = 100
CLEAN_PASS = 0.95  # the tolerance is the mapped copies' 95th z percentile: record-end LTRs pass as 95% of those do
# heuristic: a copy's Poisson K2P s.e. takes K2P >= 1e-6, so identical LTRs get a small s.e., never 0 (which would
# make every comparison with them exact)
K2P_SE_FLOOR = 1e-6
# io.library (the screen's rule, LTRmix's reading): records whose names hold no coordinates are kept, unplaced; a
# library with more than half of them so is not one with LTRquest-style names (fail fast)
UNPLACED_MAX = 0.5

# seq.align (phasing's aligner)
ALIGN_BAND = 32  # affine realignment within +-32 columns of the edlib path (indel runs in LTR-RT pairs are short)
# match, mismatch, gap open, gap extend: BLASTN log-odds scheme (~90% identity, indels ~1/10 of substitutions);
# unit-cost edit distance swaps mismatches for indel pairs (K2P x0.86 at 13.5%)
SCORES = (2, -3, 5, 2)
ALIGN_CHUNK = {"ltr": 500, "whole": 200, "self": 200}  # items per pool job of pair counts: seconds of work each
ALIGN_JOBS = 16  # most pool jobs per thread: several per worker even out uneven alignment lengths
ALIGN_WARM = 50  # library bases aligned before the pool forks, so numba compiles there (any 2 or more bases do)

# seq.kmers (k-mer lineages of a library, record sketches)
K = 21  # k-mer length: unique in plant genomes yet short enough to span one SNP per clade variant
CAP = 250  # largest clade (carrier copies) used for links: the cap the benchmark chose in 42 of 49 runs
MEM_GB = 8.0  # heuristic (phasing's): k-mer extraction holds at most this much (hash, copy) data per pass
KMER_BYTES = 32  # memory a pass plans per k-mer: its (hash, copy) record and bucket-sorted copy (24 bytes) + headroom
KMER_BUDGET = 1.5e8  # k-mers kept: exact up to 150 Mb of LTR-RTs, else 1/ceil(size/150 Mb) (bounds members ~40 M)
# k-mer runs longer than this (occurrences) are skipped: they cannot have <= CAP carriers unless repeated tens of times
# inside each copy (tandem arrays), which carry no lineage signal
RUN_LIM = 12500
CLADE_KEYS_SEED = 20261012  # seeds the per-copy keys that sign carrier sets: hash keys, fixed so clades reproduce
SKETCH_K = 15  # homolog search k-mer: 15-mers stay shared between copies ~15% apart (old splits); 21-mers do not
SKETCH_S = 128  # heuristic (phasing's): bottom-128 MinHash sketch per copy
# k-mers shared across libraries (LTRmix's orthology candidates) count only with <= 4 copies in each library: orthologs
# in young families share low-copy k-mers; it bounds the candidate pairs per k-mer to 16
CMAX_PAIR = 4
PAIR_CHUNKS = 16  # most parallel chunks of the joint scan's pair counters (copies x libraries each): bounds their memory
# the joint scan's int32 copies of the single-copy tallies, one per chunk of a pass, stay under 1 GB in all (past
# it, groups of trios add into the tallies themselves)
TALLY_CHUNK_BYTES = 1_000_000_000

# seq.chain (collinear chains of anchors in rank space; the phasing tool's numbers)
CHAIN_W = 25  # collinear anchors at most 25 copies apart on both sequences (rank space; tolerates local loss)
CHAIN_SD_FLOOR = 0.5  # a chain's z: the permutation null's spread floored at half an anchor (a discrete null)
CHAIN_ALPHA = 0.05  # a chain is significant at the conventional 5%, Bonferroni over the sequence pairs tested
# the chain test's permutations: B = ceil(2 / alpha), so its smallest p, 1 / (B + 1), is about half the level
CHAIN_PERMS_PER_ALPHA = 2.0
NULL_PERMS = 50  # heuristic: partner-position permutations per sequence pair for the null chain length

# seq.relatives (closest relatives of library copies)
BUCKET_CAP = 2000  # sketch hashes carried by more copies are family-wide: they do not single out relatives
REL_BLOCK_NNZ = 20_000_000  # products per block of the relatives self-join: bounds a worker's memory
REL_TOP = 10  # heuristic (phasing's): relatives per focal copy on other sequences
REL_MIN_SHARED = 2  # a relative shares >= 2 sketch hashes or clades with the copy: one shared can be chance

# io.annotation.bed, calls.loci, routes.assembly.classify and routes.assembly.lastz
COV_PRESENT = 0.8  # Wicker et al. 2007 80-80-80 rule: >=80% of a part aligned = present
MIN_SEG = 80  # Wicker's 80-bp minimum: less aligned sequence counts as none

# calls.loci, routes.assembly.classify and routes.reads
TSD_DEFAULT = 5  # a TSD length when none is annotated: most LTR-RTs make 4-6-bp TSDs, 5 bp the common case
# heuristic (LTRpav's): an empty site's gap may lie within about 2 TSD lengths of the element ends (classify rule 4)
SLOP_TSD_MULT = 2

# calls.genotype
RESOLVED_WARN = 0.01  # heuristic (LTRpav's): < 1% of a pair's calls resolved (P/S/T/E), no orthologous alignment; warn

# calls.records
GZ_LEVEL = 6  # zlib's default level: near level 9's size in a fraction of its time

# calls.manifest
MANIFEST = "manifest.json"  # a finished run's record in its output folder: the run's command, inputs and tools

# calls.vcf, calls.tracks, io.annotation.gff3, analyses.geneflow.run, analyses.subgenomes and sim.nucleus: how an output
# names the program that made it (the names of columns and fields, such as ##LTRpav_POS, stay as users' scripts read them)
PROGRAM = f"synltr {__version__}"  # VCF ##source and findings headers: the package and the version that wrote them
GFF3_SOURCE = f"synltr-{__version__}"  # GFF3 source column: the same without a space, which awk-style splitting breaks

# seq.junction, routes.assembly.classify and routes.assembly.lastz
REFINE_PAD = 150  # junction refinement half-width; covers local-aligner end error
SCORE = (1, 4, 6, 2)  # match, mismatch, gap open, gap extend: minimap2 asm20 short-gap scoring

# routes.assembly.classify
SOLO_LEN = (0.8, 1.2)  # a solo-LTR is one LTR's worth of aligned sequence, same 80% rule
IDENT_MIN = 0.8  # flank identity below the 80% rule: site placement not trusted
MARGIN_HIGH = 0.5  # heuristic (LTRpav's): orthologous chain outscores any competitor 2:1
IDENT_FLANK = 1000  # local A-B identity is read over the element +-1 kb: its own neighborhood

# routes.assembly.chain, anchor and window
# FastGA's default minimum alignment length (LTRpav's rule): shorter alignments are not trusted as anchors; also
# chain_blocks' minimum end advance
ANCHOR_MIN = 100

# routes.assembly.chain
CHAIN_MAX_GAP = 1_000_000  # heuristic (LTRpav's): no chaining across >1 Mb without alignment
# heuristic (LTRpav's): cost per bp of diagonal shift; a 30-kb insertion is bridged after > 1.5 kb of chain
CHAIN_INDEL_COST = 0.05
CHAIN_OVERLAP = 100  # an overlap this small between two chains' blocks is edge slop, not competition
CHAIN_LOOKBACK = 200  # heuristic (LTRpav's): predecessors examined per block in the chaining DP
MIN_CHAIN_SCORE = 1000  # heuristic (LTRpav's): a chain must carry >=1 kb of aligned sequence
CLAIM_BUCKET_BITS = 16  # claimed intervals are filed in 64-kb buckets: insert cost stays flat on long chromosomes

# routes.assembly.blocks and anchor
# Wicker's 80-bp minimum; split at longer indels, blocks stay near-collinear: A->B errs < 80 bp per indel
SPLIT_INDEL = 80
PAF_BATCH_LINES = 1_000_000  # PAF lines parsed per structured-array batch: bounds the parser's memory

# routes.assembly.anchor and index
ANCHOR_PRESET = "asm20"  # minimap2's preset for assemblies up to ~5% divergent: covers intraspecific pairs
ANCHOR_BEST_N = 5  # secondaries kept per flank-piece mapping (the rivals margins read): minimap2's generic -N
MM_F_EQX = 0x4000000  # minimap2's MM_F_EQX: CIGAR with =/X, so matches are counted per piece

# routes.assembly.anchor
PIECES_QUEUED = 2  # flank pieces in flight per thread: keeps every mapping thread busy while bounding memory
FLANK_TARGET = 2000  # unique flank bases collected per side (shorter flanks are often repeats)
FLANK_CAP = 2_000  # farthest a flank may reach: longer caps (5-20 kb) cost more, no gain (sims)
FLANK_PIECE_MIN = 200  # shorter flank pieces seed poorly and multi-map

# routes.assembly.window
WINDOW_PAD = 500  # anchored sequence kept at each window end so the local aligner can seed
WINDOW_CAP = 500_000  # runtime guard on window span (A and B)
NCAND = 3  # nearest anchor candidates tried per side: rescues a window whose nearest anchor is out of order (a
#            breakpoint) without searching far

# routes.assembly.lastz
LASTZ_SEED = "12of19"  # LASTZ's default spaced seed, one transition allowed: sensitive at 0-3% divergence

# routes.assembly.pair
CHUNK_WINDOWS = 2000  # most windows per genotyping chunk: bounds its memory and the work a crash loses
CHUNK_MIN = 8  # fewest windows per chunk: fewer add file and marker overhead
CHUNK_TARGET = 256  # chunks per pair: ~4 per worker at 64 threads; not --threads: resumes reuse chunks
U_WARN = 0.99  # heuristic (LTRpav's): > 99% U calls, the pair shares next to no anchored sequence; warn

# routes.library (the screen's evidence: orthology, clock filter, chromosome correspondence and alignment jobs)
# heuristic (LTRmix's): collinear chaining tolerates up to 20 anchors of local disorder, counted in rank space so the
# rule does not depend on anchor density: a local rearrangement or misplaced anchors that small stay in one block
CHAIN_GAP = 20
# heuristic (LTRmix's): floor on a collinear block's anchors; the minimum per chromosome pair is data-derived, one
# more than the longest chain among the same anchors with shuffled partner ranks (NULL_SHUFFLES shuffles)
CHAIN_MIN = 5
NULL_SHUFFLES = 20  # heuristic: shuffles giving that null chain length (its longest of 20: a ~95% bar)
SCREEN_SEED = 0  # seeds each pair's shuffles (one generator per genome pair, used by its chromosome pairs in turn)
# heuristic (LTRmix's): a chromosome pair needs more anchors than independence gives at Poisson p <= 1e-6 / pairs
# tested
CORR_ALPHA = 1e-6
MIN_ALN_FRAC = 0.5  # records whose lengths differ more than twofold are not aligned (orientation or record error)
# a pair's divergence passes the clock when it exceeds the copies' annotated LTR K2P by at most 1.96 sampling SE (or
# one substitution); a record's own two LTRs (self-outgroup) when their K2P is within 1.96 SE of the annotated one
SELF_Z = 1.96
# heuristic (LTRmix's): the clock test compares LTR with LTR when >= 100 LTR columns align, else whole records
CLOCK_LTR_SITES = 100
CLOCK_LTR_FLOOR = 100  # heuristic (LTRmix's): LTR length floor (bp) in the sampling variance of the annotated K2P
# strand of a record pair with an unknown strand: infix edit distances of the first record's first 3 kb against the
# second's first 3 kb + 200 bp and its reverse complement's (enough to tell strands; cheap on 10-kb records)
ORIENT_HEAD = 3000
ORIENT_SLACK = 200
PAIR_CHUNK = 200  # record pairs, or own-LTR checks, per pool job: seconds of work each
TRIPLE_CHUNK = 100  # heuristic (LTRmix's): ortholog triples per pool job (three alignments each)
POOL_JOBS = 32  # most pool jobs per thread: several per worker even out uneven record lengths

# routes.reads.reads (per-read evidence at one element)
ANCHOR = 100  # bp of aligned sequence a read needs beyond a boundary: past alignment end effects
TAU_REL = 0.03  # allele-length tolerance, fraction of the expected allele span: ONT net indel error ~1%, x3
TAU_ABS = 50  # minimum allele-length tolerance (bp): TSD length (4-6 bp) plus breakpoint misalignment
# heuristic (ltrkin's): window of internal region next to an LTR in which a full-junction read needs ANCHOR aligned bp
JUNC_WIN = 300
SLOP = 60  # bp: an alignment end "is at" a boundary within this distance (TSD + aligner shifts)
CLIP_MIN = 500  # soft clip (bp) that marks one piece of a split read: long enough to place uniquely
CLIP_PROBE = 300  # clip bases aligned to the candidate continuation: unique enough to place (ltrkin's design)
SW_MIN_ALN = 150  # heuristic (ltrkin's): aligned columns needed for a clip match
SW_MIN_ID = 0.90  # identity for a clip match: far above chance (~0.5 for unrelated DNA), below read error floors
# clip match scoring (parasail Smith-Waterman): match, mismatch, gap open, gap extend, BLASTN's scheme as in SCORES
SW_SCORE = (2, -3, 5, 2)

# routes.reads.run (the elements the route reads, the run's files) and calls.genotype
READS_TSD_ORDER = (6, 5, 4)  # ltrkin's TSD inference: the longest flank repeat of 6, then 5, then 4 bp
READS_STEM = "reads"  # a reads run's file stem (ltrkin's -o prefix): reads.genotypes.tsv, reads.vcf.gz, reads.<S>.*

# routes.reads.genotype (one read set: evidence over every element, boundary error, flags)
# heuristic (ltrkin's): bp around the expected empty length in which reads inform the annotation's boundary error
RES_WIN = 300
# heuristic (ltrkin's): elements without a TSD, each with >= 2 such reads, needed to estimate that error
BOUND_MIN_SITES = 20
BOUND_K = 3.0  # tolerance = BOUND_K x median |per-site offset|: 95th/50th percentile ratio of |error| for normal or
#               triangular (two uniform ends) errors
FLANK_POINTS = (6000, 3000)  # heuristic (ltrkin's): bp outside the element where local (flank) depth is measured
EVIDENCE_CHUNK = 200  # elements per evidence job: hundreds of jobs on a plant genome keep the workers even
PROGRESS_CHUNKS = 20  # evidence jobs between two -v progress lines (and after the last): ~4,000 elements a line
# heuristic (ltrkin's): informative depth below 1/4 of the expected flags low_depth (a lost or unreadable site),
# above 2.5x high_depth (a collapsed duplication)
DEPTH_LOW = 0.25
DEPTH_HIGH = 2.5
FLAG_FRAC = 0.2  # heuristic (ltrkin's): other or unexplained reads at >= 1/5 of a site's reads (and FLAG_MIN) flag it
FLAG_MIN = 2  # and at least this many such reads: one read alone is no pattern
LOWQ_FRAC = 0.3  # heuristic (ltrkin's): >= 30% of junction reads below --min-mapq flags a site low_mapq
PH_MAX_DIV_NOQUAL = 0.01  # reads without base qualities: phasing only if their median divergence is <= 1% (HiFi)
QUAL_FRAC_MIN = 0.9  # heuristic: a read set has base qualities when >= 90% of its profiled reads carry them
FLANK_LOW_COPIES = 0.5  # flank_depth_low: flank depth nearer one haplotype less than the full ploidy (midpoint)

# routes.reads.model (genotype model)
# heuristic (ltrkin's): read misclassification rate for pass 1; pass 2 uses the per-sample data estimate
EPS_INIT = 0.02
EPS_FLOOR = 0.001  # lower bound of the estimate (keeps likelihoods finite with few reads)
EPS_CAP = 0.1  # heuristic (ltrkin's): upper bound of the estimate
EPS_MIN_READS = 10  # heuristic (ltrkin's): informative reads a homozygous site needs to enter the error estimate
CONF_GQ = 20  # a call this confident (phred 20 = 1% error) enters the error estimate and the depth calibration
# heuristic (ltrkin's): sites with an empty/solo allele needed to calibrate depth at element sites (else genome-wide)
CAL_MIN_SITES = 30
EB_MIN_PER_CLASS = 5  # called sites per genotype class needed to learn the genotype prior (Cochran's rule for
#                      multinomial counts: 50 sites for a diploid); fewer -> flat prior
GQ_CAP = 99  # phred cap of genotype quality: the usual VCF cap (GATK caps GQ at 99)
PAV_P = 0.99  # allele-presence posterior for a presence/absence call (1% error, as GQ 20)
MU_MIN = 0.5  # expected stray informative reads at a site absent from every haplotype (allows mismaps)
# heuristic (ltrkin's): random positions sampled per BAM, for the depth mixture of 1..ploidy copies and the
# read-span distribution
N_DEPTH_POINTS = 3000
# heuristic (ltrkin's): reads whose base qualities and divergence decide whether read phasing is safe
PROFILE_READS = 2000
SEED = 7  # random seed for sampling (reproducible runs)
DEPTH_MIN_POINTS = 20  # heuristic: covered random positions needed to model depth (else the read set is refused)
DEPTH_MIN_READS = 50  # heuristic: and reads, for the read-span distribution
DEPTH_COLLAPSE = 3  # heuristic: depth above 3 x the median is a collapsed repeat, not copies of the genome
DEPTH_GRID = (0.4, 2.5, 46)  # heuristic: haplotype depth searched, 0.4-2.5 x median / ploidy in 46 log steps
NB_SIZE_GRID = (2.0, 1000.0, 16)  # heuristic: negative-binomial size (overdispersion) searched, 2-1000 in 16 log steps
DEPTH_FINE = (-0.045, 0.045, 19)  # the refinement: about +-1 depth grid step (log units) around the best, 19 points
NB_SIZE_FINE = (-0.4, 0.4, 9)  # and about +-1 size grid step, 9 points
MIXTURE_ITERS = 30  # heuristic: EM rounds of the copy-number weights at each grid point
NB_SIZE_MIN, NB_SIZE_MAX = 5.0, 1000.0  # heuristic: the size kept between very overdispersed (ONT) and Poisson-like
EB_ITERS = 20  # heuristic: EM rounds of the empirical-Bayes genotype prior
CAL_MIN_EMPTY = 2  # a calibration site has >= 2 empty or solo reads: one read alone is no pattern (FLAG_MIN)

# routes.reads.phase (distinct haplotypes per allele) and routes.reads.model (their weight)
PH_FLANK = 3000  # bp of unique flank on each side compared between reads: ~30 het sites at 1% heterozygosity
# heuristic (ltrkin's): floor of the minor-allele fraction of a heterozygous site (default max(0.15, 0.5/ploidy))
PH_MIN_AF = 0.15
PH_MAX_SITES = 60  # sites nearest the element used per site (bounds work; more add little)
PH_MIN_Q = 20  # base quality for a phasing allele (error <= 1%): on simulated ONT from identical haplotypes,
#               spurious second haplotypes fell from 13-44% (no filter) to <= 0.4%; real HiFi unaffected
# disagreeing sites that prove two read groups distinct: requiring 3 cost 5 points of detection on real HiFi
# (ltrkin's benchmark)
PH_MIN_DIS = 2
PH_MAX_READS = 80  # reads compared per element (read-name-hash subset): 40 per haplotype of a diploid is plenty;
#                   bounds time and the reads x reads matrices in collapsed repeats
PH_ERR = 0.05  # spurious extra haplotype group: 4.9% of confirmed one-haplotype sites (C. amara HiFi 20x; paralog
#               reads from sequence absent from the reference haplotype)
PH_MIN_READS = 4  # heuristic: fewer reads at an element: no haplotypes counted
PH_MIN_ALLELE_READS = 2  # an allele, or a haplotype group, needs >= 2 reads: one read can be a sequencing error
PH_MIN_LINK_READS = 4  # two sites are tested for linkage on >= 4 reads covering both (one per gamete)
PH_GAMETE_FRAC = 0.1  # heuristic: four-gamete test, linked when the rarest gamete holds <= 10% of those reads
PH_SAME_SITES = 2  # heuristic (ltrkin's): two reads are one haplotype when they share >= 2 sites ...
PH_SAME_DIS = 0.1  # heuristic: ... and disagree at <= 10% of them (sequencing error)
PH_DISTINCT_DIS = 0.7  # heuristic: two read groups are distinct haplotypes when they disagree at >= 70% of shared sites

# analyses.source (ltrkin's validated thresholds; the search is two-stage: a minimap2 screen of every candidate, then an
# exact mappy alignment of the best few)
SOURCE_DIR = "source"  # folder of the search's cached alignments, under the run's work folder
GQ_STATUS = 20  # genotype quality for a status call (phred 20 = 1% error)
MIN_INTERNAL = 1000  # bp of internal region for a candidate: resolves 0.1% divergence (1 difference)
MIN_ALN = 1000  # heuristic (ltrkin's): aligned bp for a relative
MIN_COV = 0.8  # heuristic (ltrkin's): aligned fraction of the shorter internal region for a relative
MIN_COV_LONG = 0.5  # aligned fraction of the longer one: a short internal region must not match inside a long one
TIE_DIFFS = 0.5  # heuristic (ltrkin's): candidates within half a difference of the best are tied
TOP_K = 20  # candidates per new copy re-aligned exactly after the approximate screen (completeness benchmarked)
N_NEAREST = 10  # heuristic (ltrkin's): nearest other relatives used as the local consensus for shared-mutation counts
CONS_MIN = 3  # heuristic (ltrkin's): relatives needed at a site to call the local consensus base
CONS_FRAC = 0.5  # and the consensus base is carried by >= half of the relatives there (a majority)
NEAREST_SPARE = 5  # heuristic: relatives read beyond N_NEAREST, spares for those the checks drop
FLANK = 500  # heuristic (ltrkin's): bp of flank compared between two copies (segmental duplication check)
# heuristic (ltrkin's): flanks this identical (FLANK bp) mark a segmental duplication, not a new insertion: far above
# the identity of unrelated flanks (~0.5 aligned)
FLANK_DUP_ID = 0.9
OLDER_DIFFS = 2  # heuristic (ltrkin's): a source is clearly older with >= 2 more LTR differences (one: an error)
MOSAIC_ALPHA = 0.01  # heuristic (ltrkin's): the two-parent scan's family-wise level, 1% (Bonferroni over pairs)
MOSAIC_COV = 0.9  # heuristic (ltrkin's): a mosaic candidate covers >= 90% of the new copy's internal region
MOSAIC_GAIN = 3  # heuristic (ltrkin's): a mosaic explains >= 3 more differences than the best single candidate
# heuristic (ltrkin's): cumulative source-divergence limits reported in locality
DIV_LIMITS = (0.001, 0.002, 0.005, 0.01, 1.0)
PERM = 5000  # heuristic (ltrkin's): null replicates of the locality tests
LOCALITY_SEED = 11  # random stream of the locality tests
LOCALITY_MIN_NEW = 5  # heuristic (ltrkin's): new copies a family needs for its own locality row
EXACT_CHUNK = 8  # screened new copies per worker task: amortizes the hand-over of results
EXACT_CACHE = 256  # aligner indexes a worker keeps (memory, not results)

# stats
BLOCK_BP = 5_000_000  # jackknife block: >> LD in outcrossing grasses, hundreds of blocks per genome
BLOCK_FRAC = 100  # or 1/100 of the smallest genome, so a small genome still gives ~100 blocks
BLOCK_MIN = 10_000  # heuristic (LTRpav's): smallest block size the rule gives
JK_MIN_BLOCKS = 3  # heuristic (LTRpav's): fewest distinct blocks a delete-one jackknife takes a variance from
Z95 = 1.96  # the normal quantile of a two-sided 95% interval (the convention), as when a model error widens one
RESOLVE_WARN = 0.1  # heuristic: warn when two genomes' resolved fractions differ by > 0.1

# analyses.geneflow (nothing here is tuned to a dataset)
AGE_EDGES = (0.0, 0.01)  # heuristic (LTRpav's): age bins by LTR-LTR K2P, identical, below 1%, older (few bins)
Z_SIG = 3.0  # |Z| >= 3 counts as significant: the conventional bar for D statistics
P_WARN = 0.05  # conventional level of the sister-pair and lineage-sorting symmetry tests
# the analyses import these where they use them; geneflow imports them first, so a broken environment fails before
# any genotyping
ANALYSIS_IMPORTS = ("scipy.optimize", "scipy.special", "scipy.stats")
# the headline D re-jackknifed with blocks from 1 Mb to 20 Mb (stats.tsv test block_scan): the 5-Mb default and up
# to 5x either side, so an se that still grows with the block length (linkage longer than the blocks) shows
BLOCK_SCAN_BP = (1_000_000, 2_000_000, 5_000_000, 10_000_000, 20_000_000)
# heuristic: at a scan length whose informative loci lie in fewer blocks the jackknife se is itself unstable (it moves
# only through the blocks holding informative loci, as stats.report.why counts testability); never warned on
BLOCK_SCAN_MIN_BLOCKS = 30
PI_MAX = 0.5  # donor-like fraction <= 1/2: a recipient more donor-like than sister-like is the donor's sister
PI_MIN = 1e-4  # heuristic (LTRmix's): painting pi's lower clip, the fit's floor where state 1 holds no data
BF_TRACT = 3.0  # tract loci need posterior odds >= 3x prior odds (Kass & Raftery 1995, "positive")
Q_TRACT = 0.05  # a tract is reported donor-like at target-decoy q <= 0.05: the conventional 5% FDR
HMM_PI0, HMM_L0 = 0.1, 1e5  # heuristic (LTRmix's): EM's starting donor-like fraction and switch scale (bp)
# heuristic (LTRmix's): EM's initial state 1, donor-like classes (RD) enriched 5x, the others (RS) depleted 5x
HMM_INIT = (5.0, 0.2)
HMM_L_GRID = (1e4, 5e7, 41)  # heuristic (LTRmix's): mean tract lengths tried by EM, 10 kb-50 Mb in 41 log steps
HMM_ITERS = 100  # heuristic (LTRmix's): EM rounds at most
HMM_TOL = 1e-7  # heuristic (LTRmix's): EM stops when the log-likelihood moves by < 1e-7 of itself
# the library route's painting: an ortholog is near below the 5th percentile of the sister-donor ortholog divergences,
# shared ancestry younger than 95% of sister-donor pairs (lineage sorting makes it equally for both sisters)
NEAR_QUANTILE = 0.05
CALLABLE_FLANK = 20  # heuristic: RS loci on each side of the shortest-callable upper bound
CALLABLE_MAX = 200  # heuristic (LTRpav's): shortest-callable search limit (loci); a longer minimum is not callable
INV_MIN_LOCI = 5  # heuristic: shortest reversed run reported as an inverted block; named in findings.txt
# contiguity: a block's extent on each genome holds fewer loci placed on non-partner sequences than block loci, the
# merge rule's own tolerance; cost: a real inversion whose region's second copy lies on a non-partner sequence (or
# beside a pair kept from one side only) is dropped
INV_FOREIGN_MAX = 1.0
# heuristic: on the library route a chromosome pair needs >= 50 anchors in its chain blocks before its main
# orientation (the sign holding more anchors) is trusted to call the other sign's blocks inverted; named in findings.txt
CHAIN_INV_MIN = 50
# heuristic (LTRpav's): the mixture's post-split gene-flow component has mean divergence <= 0.1%, below the
# splits the isolation component measures (its upper tail holds old paralogs)
MG_MAX = 1e-3
IDENT_MIN_SITES = 200  # LTR sites an ortholog pair needs to enter the split fit (one substitution is 0.5% of 200)
MIX_MIN_PAIRS = 1.0  # mixture: a component with weight x pairs < 1 pair sits at its weight's bound, 0: parameters NA
# heuristic (LTRpav's): divergences where the split's fit check compares observed and fitted
CDF_POINTS = (0.0, 0.001, 0.002, 0.003, 0.005)
IDENTITY_CHUNKS_PER_WORKER = 8  # identity alignment chunks per worker: a slow chunk delays the end by ~1/8 of a share
Q_MIN = 1e-280  # gammaincc nears its underflow (~1e-308) below this: the exact tail sum takes over
TAIL_EPS = 40.0  # tail sum: stop once the terms decrease and the last is < e^-40 of the sum (double precision)
SPLIT_FIT_START = (1e-3, 1e-2)  # heuristic: the isolation fit's start (s, m), a 0.1% shift and a 1% mean
# heuristic (LTRmix's): its Nelder-Mead tolerances (log s, log m) and iterations at most
SPLIT_FIT_TOL, SPLIT_FIT_MAXITER = 1e-6, 4000
# heuristic (LTRmix's): the mixture fit's 18 Nelder-Mead starts around the isolation fit, gene-flow weights 0.5-5%,
# isolation shifts 1-6x the isolation fit's, tail means 5-15%; isolation and tail start at weights 0.95 and 0.03
# beside them
MIX_START_FLOW = (0.005, 0.02, 0.05)
MIX_START_SHIFT = (1, 3, 6)
MIX_START_TAIL = (0.05, 0.15)
MIX_START_WEIGHTS = (0.95, 0.03)
MIX_FIT_TOL, MIX_FIT_MAXITER = 1e-7, 20000  # heuristic (LTRmix's): its tolerances and iteration cap (six parameters)
# heuristic (LTRmix's): simulations of the observed sites behind the split fit's predicted fractions (CDF_POINTS)
CDF_SIMS = 200
# the HMM's emission frequencies take half a count per class (Jeffreys' prior): an unseen class is rare, not
# impossible
EMIT_PSEUDO = 0.5
TRACT_POST = 0.5  # a tract locus: the donor-like state is the more probable one (and BF_TRACT's odds)

# analyses.geneflow.library_stats (the library screen's statistics; the rules they were validated with)
# heuristic (LTRmix's): D_ins age strata by the donor copy's LTR-LTR K2P, identical copies, then bins up to 10%
# (unknown K2P sorts last)
SCREEN_AGE_EDGES = (0.0, 1e-9, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1)
SCREEN_AGE_LABELS = ("0", "(0,0.002]", "(0.002,0.005]", "(0.005,0.01]", "(0.01,0.02]", "(0.02,0.05]", "(0.05,0.1]",
                     ">0.1")
# heuristic (LTRmix's): an age stratum of D_ins is reported with >= 10 informative insertions (|D| has se <= 0.32)
SCREEN_AGE_MIN = 10
SCREEN_CHROM_MIN = 20  # heuristic (LTRmix's): a donor chromosome's D_ins with >= 20 informative insertions
# heuristic (LTRmix's): a chromosome's k-mer share and D_kmer with >= 1000 informative k-mers (linked k-mers count
# far less than their number)
SCREEN_KMER_CHROM_MIN = 1000
ROOM_FRAC = 0.5  # an absence is certain where the syntenic interval is under half the element: too short to hold it
# old insertions: LTR-LTR K2P above the genome's 75th percentile, inserted long before any plausible split and so
# inherited by every genome of the test
OLD_QUANTILE = 0.75
# heuristic (LTRmix's): the completeness ratio needs > 100 of the donor's old insertions resolved in both sisters
COMPLETENESS_MIN = 100
COMPLETENESS_Z = 2.0  # the completeness ratio differs from 1 beyond 2 se (~95%, two-sided)
# heuristic: warn when < 80% of the donor's elements lie in synteny with both sisters, or the two sisters' resolved
# fractions differ by > 0.1 (insertion and ortholog tests then cover only the jointly resolved part)
RESOLVED_BOTH_WARN = 0.8
RESOLVED_GAP_WARN = 0.1
# heuristic (LTRmix's): ortholog divergence bins of the identity spectrum, identical, then doubling from 0.025% to
# 6.4%, then the rest
IDENTITY_BINS = (0.0, 1e-9, 2.5e-4, 5e-4, 1e-3, 2e-3, 4e-3, 8e-3, 1.6e-2, 3.2e-2, 6.4e-2, 1.0)
# findings compare target and mirror orthologs below 0.2% divergence (recent shared ancestry): the identity
# spectrum's bins up to 2e-3 (IDENTITY_BINS), one substitution in 500 LTR sites
IDENTITY_RECENT = 0.002
# the evenness check passes when every pair's ratio, or its inverse, is >= 0.8/0.9: the largest annotation gap (10
# points below 90% complete) at which D3 and D_self held the null in simulation; larger gaps are held back, not judged
EVEN_CUTOFF = 0.889
# heuristic: the evenness ratio's jackknife blocks (its se is printed, never used by the rule)
EVEN_BLOCK_BP = 1_000_000
# false signals (|Z| > 3) per statistic in simulated nulls, (signals, nulls) by annotation; D_pav (geneflow without
# --screen) beside them
FALSE_SIGNALS = {
    "D_kmer": {"even annotation": (0, 30), "even annotation, bottleneck": (0, 30),
               "one sister 10 points less complete": (5, 10), "one sister 30 points less complete": (9, 10)},
    "D_ins": {"even annotation": (0, 30), "even annotation, bottleneck": (0, 30),
              "one sister 10 points less complete": (9, 10), "one sister 30 points less complete": (9, 10)},
    "D3_triples": {"even annotation": (0, 30), "even annotation, bottleneck": (0, 30),
                   "one sister 10 points less complete": (0, 10), "one sister 30 points less complete": (4, 10)},
    "D_self_outgroup": {"even annotation": (0, 30), "even annotation, bottleneck": (1, 30),
                        "one sister 10 points less complete": (0, 10), "one sister 30 points less complete": (10, 10)},
    "D_pav": {"one sister 10 points less complete": (0, 10), "one sister 30 points less complete": (0, 10)},
}
# the evenness check on simulated runs (its lowest pair; nothing fitted to them): (the runs, how many, the range of
# their check values, how many it held back)
EVEN_MEASURED = (
    ("even annotation: 30 nulls and 58 gene-flow runs", 88, (0.948, 0.997), 0),
    ("even annotation: 30 bottleneck nulls", 30, (0.947, 0.994), 0),
    ("one sister 10 points less complete", 13, (0.859, 0.915), 10),
    ("one sister 30 points less complete", 13, (0.636, 0.697), 13),
)

# plot.style, plot.geneflow and plot.subgenomes, and the analyses whose figures they draw (analyses.geneflow.run,
# analyses.subgenomes.run: a run removes the earlier figure when it starts, without importing plot)
FIGURE_SUFFIXES = (".png", ".pdf", "_numbers.tsv")  # a figure's files (plot.style.save): page, page and legends, numbers
OVERVIEW_STEM = "overview"  # the gene-flow overview's files: overview.png, overview.pdf, overview_numbers.tsv
REPORT_STEM = "report"  # the subgenomes page's files: report.png, report.pdf, report_numbers.tsv

# plot.style
FIG_WIDTH_MM = 180.0  # the widest a journal page takes (single column 89 mm; two columns 180 mm)
FIG_DPI = 300  # raster resolution journals ask for at minimum
FONT_SIZE, FONT_SIZE_SMALL, FONT_SIZE_LETTER = 6.5, 6.0, 9.0  # legible at 180 mm; panel letters larger and bold
LETTER_DX_MM, LETTER_DY_MM = 6.0, 0.5  # panel letter sits this far left of and above its axes' top-left corner
KEY_LINE_MM = 3.0  # height of one line of a key row
KEY_GAP_MM = 3.2  # space between key entries
LEGEND_FONT_PT, LEGEND_LEAD_PT, LEGEND_PARA_GAP_PT = 10.0, 12.6, 6.0  # legend page: 10 pt text, 1.26 line spacing
LEGEND_MARGIN_MM = 12.0  # legend page margin
NUMBERS_COLUMNS = ("key", "value", "source")  # the numbers file's header: every plotted value and its source table
GRAY = "#8F8F8F"  # scaffold and context gray of every figure (key outlines, hatching)
KEY_GLYPH_MM = 2.6  # width of a key glyph
KEY_ROLE_GAP_MM = 2.6  # space after a role label in a key row
KEY_TEXT_GAP_MM = 0.9  # space between a key glyph and its text
KEY_ROLE_PAD_MM = 0.3  # padding added to a bold role label's width
KEY_TICK_X_MM, KEY_TICK_HALF_MM, KEY_TICK_LW = 1.2, 0.95, 0.8  # a tick glyph: x offset, half-height (mm), line width
KEY_BAND_HALF_MM, KEY_BAND_H_MM = 1.05, 2.1  # a band glyph: half-height and height (mm)
KEY_HATCH_HALF_MM, KEY_HATCH_H_MM = 0.65, 1.3  # a hatch glyph: half-height and height (mm)
KEY_CMAP_MM, KEY_CMAP_PAD_MM, KEY_CMAP_STEPS = 10.6, 1.6, 30  # a color-scale glyph: width, inset of its bar, steps
KEY_TRACT_LW, KEY_STRONG_LW = 0.9, 1.6  # line widths of the tract glyphs: tract, and tract at the strict q level
DARK = "#303030"  # dark scaffold and tract lines
KEY_HATCH = "//////"  # hatch pattern of the "no data" glyph, dense enough to read at 1.3 mm
LEGEND_LINE_FILL = 0.985  # share of the text width a legend line may use: leaves slack for font-metric rounding
LEGEND_PAD_PT = 4.0  # extra page height at the foot of a legend page
# matplotlib rc of the gene-flow overviews (assemblies and the screen; plot.style adds the fonts): thin lines, small
# ticks, 6-6.5 pt text, text kept as text in PDF (fonttype 42). The subgenomes page keeps its own (FIG_RC): the two
# differ in 19 of their 25 entries (text 6.5 vs 7 pt, lines 0.5 vs 0.6 pt, ticks, spines), so one dict would move
# every pixel of one page or the other
RC_BASE = {
    "font.family": "sans-serif", "font.size": FONT_SIZE, "pdf.fonttype": 42, "ps.fonttype": 42,
    "svg.fonttype": "none", "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2.0, "ytick.major.size": 2.0, "xtick.labelsize": FONT_SIZE_SMALL,
    "ytick.labelsize": FONT_SIZE_SMALL, "xtick.major.pad": 1.5, "ytick.major.pad": 1.5, "hatch.linewidth": 0.35,
    "lines.linewidth": 0.6, "axes.labelsize": FONT_SIZE, "axes.labelpad": 2.0, "legend.frameon": False,
    "mathtext.fontset": "custom"}

# plot.geneflow (display limits and the example's checks; none is a data threshold)
MAX_ROWS = 24  # sequences drawn in panel A: more would be unreadable on one page
MIN_FRAC = 0.02  # sequences shorter than 2% of the longest are not drawn: they would be a speck
BINS = 250  # display bins along the longest drawn sequence (color, hatching): about the page's width in pixels
FLANK_CHECK = 20  # heuristic: bp of each flank compared between R and S's empty site when choosing the example
FLANK_IDENT = 0.9  # heuristic: identity those flanks must reach for the example
LTR_RATIO = 0.9  # heuristic: both LTRs whole when the shorter is at least this share of the longer
FAMILY_MIN = 10  # heuristic: a family with fewer elements is judged by its superfamily's length range
FLANK_SHOW, LTR_SHOW = 10, 10  # bp of flank and of each LTR end drawn in panel C: all that fits a base-pair alignment
RIGHT_MARGIN_MM = 4.0  # right edge of the data areas of panels A, B and D, from the page's edge: one grid
H_PITCH = 2.6  # mm per trio row in panel H: 6-pt labels (about 2.2 mm tall) never touch
H_ROWS = 10  # rows per column of panel H before another column is added (width allowing)
H_FOREST_MIN = 16.0  # mm: the narrowest D axis panel H draws beside its labels; fewer columns when labels are long
ACC, BLUE = "#D55E00", "#0072B2"  # vermillion for donor-sharing (RD), blue for sister-sharing (RS): Okabe-Ito pair
LIGHT, BAND = "#DCDCDC", "#ECECEC"  # light grays: bars and inverted-block bands
LTRC, INNER = "#8C3B00", "#F0B48A"  # dark and light of the element glyph's LTR and internal region
SCREEN_ALPHABET = "recent"  # the painting a library screen's overview draws: near-identity, its headline tracts

# sim (the simulator's defaults and limits; `synltr sim pop --help` prints them without loading the simulator).
# SIM_DEMOG and SIM_RATES are LTRmix's simulation of maize (A), parviglumis (B) and mexicana (C), a generation a
# year, where an entry gives no other source
SIM_DEMOG = {
    "ne_a": 30_000,  # diploid effective size of taxon A
    "ne_b": 100_000,  # of taxon B
    "ne_c": 70_000,  # of taxon C
    "ne_bc": 150_000,  # of the ancestor of B and C
    "ne_o": 150_000,  # heuristic (LTRpav's): of the outgroup O, as of the ancestor of B and C
    "ne_root": 150_000,  # heuristic (LTRpav's): of the ancestor of BC and O, the same
    "t_ab": 9_000,  # generations ago A split from B
    "t_bc": 60_000,  # B split from C
    "t_o": 900_000,  # O split from BC: measured on real genomes, from their outgroup conflict
    "f": 0.0,  # share of the recipient's lineages the pulse moves: no gene flow unless asked
    "ta": 5_000.0,  # pulse time, generations ago
    "pulse": "C>A",  # donor>recipient of the pulse
    "per_taxon": 2,  # heuristic (LTRpav's): genomes sampled per taxon
}
SIM_RATES = {
    "mu": 3e-8,  # substitutions per bp per generation
    "r": 1e-8,  # crossovers per bp per generation
    "r_low": 1e-9,  # crossovers per bp in the low-recombination interval
    "low": (0.3, 0.7),  # that interval, as fractions of each region
    "lam": 1.2e-10,  # insertions per bp per generation
    "delta": 1 / 900_000,  # losses per element per generation
    "indel": 3e-9,  # heuristic (LTRpav's): small indels per bp, mu / 10
    "copy_noise": 0.005,  # substitutions per base of a new copy
    "loss_split": (0.123, 0.164, 0.511, 0.202),  # shares of solo, del_flank, del_partial, replace: measured on real genomes
    "replace_mean": 13_538,  # mean replacement length in bp: measured on real genomes
}
SIM_PULSES = {"C>A": ("C", "A"), "A>C": ("A", "C"), "C>A1": ("C", "A1p"), "C1>A1": ("C1p", "A1p"),
              "A1>C1": ("A1p", "C1p")}  # pulse code -> (donor, recipient) population, forward in time
SIM_MAX_DONORS = 2000  # heuristic (LTRpav's): donor templates sampled from the annotation
SIM_MAX_CHROMS = 99  # seqids c01..c99: the event ids use two digits, and the random streams keep apart below 100
SIM_MAX_SEED = 2**32 - 1  # a seed is one 32-bit word of the random streams' entropy
SIM_FAI_CACHE = "fai"  # the output folder's own cache of the FASTA index, for a read-only FASTA folder
# heuristic (LTRpav's): chance an alive element is missing from its genome's annotation (--dropout's default)
SIM_DROPOUT = 0.1
SIM_JITTER = 0.05  # heuristic (LTRpav's): chance an annotated element's ends move (--jitter's default) ...
SIM_JITTER_BP = 3  # ... by up to 3 bp each: an annotator's boundary error
SIM_GAP_PER_MB = 2.5  # heuristic (LTRpav's): assembly gaps per Mb of a simulated genome
SIM_GAP_LEN = 100  # assembly gaps: 100-bp N runs (the length of most assembly gap placeholders)
SIM_GAP_NEAR = 20  # heuristic (LTRpav's): a near gap starts or ends within 20 bp of an element end or an empty site ...
SIM_GAP_NEAR_DIV = 3  # ... and one gap in 3 is near (they test the calls next to gaps)
# heuristic (LTRpav's): flank deletions take the element plus 50-2000 bp on one or both sides
SIM_FLANK_DEL = (50, 2000)
SIM_PARTIAL = (0.2, 0.6)  # heuristic (LTRpav's): partial deletions take 20-60% of the internal region ...
SIM_PARTIAL_MIN = 1000  # ... of elements with >= 1 kb of it
SIM_MAX_DEPTH = 2  # heuristic (LTRpav's): nested insertions in simulated insertions, at most two deep
SIM_INDEL_P = 1 / 3  # heuristic (LTRpav's): point-mutation indel lengths, geometric with mean 3 bp

# sim nucleus (forward simulation of LTR-RT evolution through splits and polyploid mergers; the scenario file's keys
# override these, and `synltr sim nucleus --help` prints them without loading the simulator): the defaults of
# the phasing tool's benchmark simulator, heuristic where an entry gives no other reason
SIM_NUCLEUS = {
    "mu": 1.3e-8,  # substitutions per site per year: a plant nuclear rate
    "kappa": 2.0,  # transition/transversion rate ratio of the K2P model: plant LTR-RTs
    "root_age": 1.0e7,  # years: start of the root lineage; copies older than about 5 half-lives are rare
    "tree": {"split": 2.0e6, "children": ["A", "B"]},  # the root lineage splits into leaves (subgenome progenitors)
    "mergers": [{"time": 1.0e5, "lineages": ["A", "B"]}],  # years ago, leaves whose nuclei merge
    "chromosomes": {"n": 8, "length": 4.0e7},  # chromosomes per leaf and their length in bp
    "families": {"n": 400, "masters": [1, 3], "branch_rate": 3e-7, "ltr_len_median": 600, "ltr_len_sigma": 0.8,
                 "int_len_median": 4500, "int_len_sigma": 0.4, "size_sigma": 2.0, "gc": 0.40},  # per-family draws
    "activity": {"rate": 0.004, "lineage_sigma": 1.0, "p_off": 0.2, "p_burst": 0.3, "burst_width": 1.0e5,
                 "burst_gain": 10.0},  # insertions per year in all families, and how it varies by lineage
    "target_weight": {},  # leaf -> relative insertion target weight after a merger (default 1)
    "phi": 1.0,  # cross-subgenome insertion preference after a merger (1 = by target weight only)
    "p_copy_source": 0.3,  # share of insertions sourced from a recent copy rather than a master
    "active_age": 1.0e5,  # years: copies younger than this can be sources
    "half_life": 1.5e6,  # years: loss of intact copies (deletion, solo LTR)
    "birth_divergence": 0.0,  # mean substitutions between the two LTRs at insertion
    "gene_conversion": 0.0,  # LTR-LTR conversion events per copy per year
    "conversion_tract": 200,  # bp per conversion event
    "family_rate_sigma": 0.0,  # log-sd of family substitution-rate multipliers (mean 1)
    "exchanges": {"nonreciprocal": 4, "reciprocal": 2, "size_min": 2.0e5, "size_max": 5.0e6},  # homoeologous exchanges
    "exchange_list": [],  # explicit exchanges (override the random ones)
    "singleton_share": 0.14,  # copies written as one-copy families: LTRquest's median over 6 libraries (range 0.07-0.34)
}
SIM_NUCLEUS_CLADES = {"Gypsy": ["Tekay", "CRM", "Athila", "Retand", "Reina", "Ogre", "TatII"],
                      "Copia": ["Ale", "SIRE", "Ivana", "Tork", "Angela", "Bianca", "TAR", "Ikeros"]}  # plant clades
SIM_NUCLEUS_LTR_LEN = (100, 4000)  # heuristic (that simulator's): clip of a family's LTR length, bp
SIM_NUCLEUS_INT_LEN = (1000, 12000)  # heuristic (that simulator's): clip of a family's internal-region length, bp
# heuristic (that simulator's): superfamily names and their draw probabilities
SIM_NUCLEUS_SUPERFAMILIES = (("Gypsy", "Copia", "unknown"), (0.5, 0.4, 0.1))
SIM_NUCLEUS_CACHE = 4096  # heuristic (that simulator's): sequence states kept in memory (least recently used out first)
SIM_NUCLEUS_K2P_RATE = 3e-8  # substitutions per site per year the k2p_time column assumes (LTRquest's default)
SIM_NUCLEUS_SEQID = "chr{:03d}"  # sequence names, numbered in a random order: a name never tells its subgenome

# analyses.subgenomes (phasing a polyploid's LTR-RT library into subgenomes; the phasing tool's numbers)
REPORT_DATA = "report_data.json"  # a subgenomes run's figure data, read by plot.subgenomes and plot.run
# shared by its steps
# 'significantly' above or below another value: by 2 standard errors (~95%), the phasing tool's margin wherever it
# compares an estimate with one (duplicated insertions, sibling pairs, a replication near its bar)
SE_MULT = 2
CI95 = (0.025, 0.975)  # a bootstrap's 95% interval: its 2.5th and 97.5th percentiles (the convention)
LR95 = 1.92  # a likelihood-ratio 95% interval: within half the chi-square(1) 95% quantile (3.84 / 2) of the best
# a zero count taken as half a count (the continuity correction) where a log or a proportion's bounds need one
HALF_COUNT = 0.5
SG_LTR_MIN_SITES = 200  # an LTR-LTR comparison counts with >= 200 aligned sites: one substitution is then 0.5%
PRIOR_POOL_MAX = 20000  # heuristic: copies drawn at most for a g-model age prior (run time)
# inputs
MIN_ELEMENTS = 10  # heuristic: sequences with fewer LTR-RTs are not phased (they cannot be assessed)
# an annotation whose K2P values mostly (> half) reach 0.5 gives identities, not distances: the LTR-LTR K2P of
# intact copies stays far below 0.5 (saturation is near 0.75)
K2P_IDENTITY_LIKE = 0.5
# links
NMAX_KNN = 50  # relatives for nearest-relative votes come from clades of <= 50 copies (benchmark)
K_NEIGH = 3  # each copy votes for the 3 other sequences holding its closest relatives (benchmark)
# heuristic: closest relatives kept per copy and other sequence, spares for the votes a half-library masks
REL_KEEP_MAX = 4
REL_KEEP_CELLS = 5e7  # fewer kept when copies x sequences x kept would pass 5e7 entries (memory)
# heuristic: fixed-point rounds of the degree-corrected null at most (a warning when not converged)
EXPECTED_ITERS = 200
# partition and the number of subgenomes
EXACT_MAX_SETS = 26  # exhaustive two-group search up to 2^25 configurations (Gray code: ~1 s)
MAX_N = 12  # largest number of groups searched (palette and partition search limits)
MIN_GROUP = 2  # a group holds >= 2 sequences: a lineage shared within a group needs two of them
CORE_MIN = 4  # heuristic: with fewer core sequences (a fragmented assembly) every phased sequence is core
LOCAL_ROUNDS = 200  # heuristic: local-search rounds at most from one start (it stops once no move gains)
KMEANS_ITERS = 50  # heuristic: k-means rounds of the spectral start of a partition into more than two groups
PARTITION_RESTARTS = 30  # heuristic: random starts besides the spectral one, beyond the exact two-group search
SPLIT_RESTARTS = 10  # heuristic: random starts of a group's two-way split beyond the exact search
# the structure (N, partition) is decided by 'core' sequences carrying >= 0.1 x the copy-N50 LTR-RTs (within an order
# of magnitude of a typical sequence); smaller ones (unplaced scaffolds) are placed afterwards and never create or veto
# a group
CORE_FRAC = 0.1
MIN_REP = 0.8  # a split is kept when half-libraries reproduce it (adjusted Rand; prediction-strength convention)
MIN_REL = 0.5  # heuristic: a sub-split must be >= this fraction as strong as its parent split on the same copies
# equal-count LTR-K2P terciles searched for splits: the fewest that can show a signal confined to a window of insertion
# ages (absent at one end of the age range)
AGE_STRATA = ("young", "middle", "old")
# a split carried by one LTR-age class (three extra tests beside the pooled copies) must meet the strict end of the
# prediction-strength convention (0.8-0.9) to keep the family-wise false positive rate near the pooled test's
MIN_REP_AGE = 0.9
HALF_PAIRS = 50  # heuristic: complementary half-library pairs (100 halves) for replication, support, partners
# heuristic: half-library pairs drawn at most when replicability straddles MIN_REP (precision, not leniency)
MAX_PAIRS = 200
SUPPORT_OK = 0.9  # support at which a sequence counts as assigned (benchmark: >= 97% correct)
# pairing map: paired sequences share insertions along most of both (homologues); local blocks (segmental
# duplications, exchanges) span a small part of them
PAIR_SPAN = 0.5
# homoeologs
MIN_ANCHORS = 3  # shared insertions that make a chromosome a candidate homoeolog (a chain needs >= 3 points)
CHAIN_MIN_CANDIDATES = 4  # heuristic: chromosome pairs sharing fewer anchors are not tested for a collinear chain
HOMOEOLOG_ALPHA = 0.05  # a homoeolog call's chain test at the conventional 5%, Bonferroni over its candidates
PAIR_SHARED = 0.5  # pairing map: most (> half) of a pair's anchors diverge less than their copies' own LTRs
# what the homoeolog pairs say (each bar fixed before it was scored on simulations and real polyploids)
HB_ALPHA = 0.05  # family-wise level of every test (Bonferroni within its family)
HB_MIN_RUN = 2  # an exchanged stretch needs >= 2 consecutive collinear pairs (one pair can be one converted copy)
HL_FACTOR = 1.5  # copies older than 1.5 x the split are pre-split (retention levels off there in simulations)
HL_MIN_R = 0.06  # retention below 6% leaves the half-life to noise (simulations)
HL_PLATEAU = 0.35  # retention at 2 x the split within 35% of retention at 1.5 x: pre-split copies behave alike
HL_BOOT = 200  # heuristic: bootstrap replicates (chromosomes; inherited pairs behind the split divergence)
SPLIT_MIN = 20  # heuristic: inherited pairs needed for the split divergence the half-life is measured against
HL_BOOT_MIN = HL_BOOT // 2  # heuristic: a half-life interval needs more than half its replicates finite
HL_CHECK_FACTOR = 2  # the plateau check also reads retention at 2 x the split (HL_PLATEAU)
HOMOG_MIN = 5  # heuristic: anchors needed for the empirical null of homogenized pairs (a median, an upper half)
HOMOG_ITERS = 30  # heuristic: rounds at most of that null, refitted on the anchors it leaves unflagged
HOMOG_SCALE_FLOOR = 1e-3  # heuristic: the null's spread (log divergence) never below 1e-3: tied anchors flag nothing
# dating
# heuristic: focal copies per stratum (random; youngest 20% by LTR K2P), ~60k copy pairs for the dating
FOCAL = 3000
YOUNG_MAX = 12000  # heuristic: young focal copies in all for the merger clock; more where the library has them
# heuristic: orthologous insertions needed for the ortholog split (fewer: the profile and clock estimates)
N_ORTH = 20
BOOT = 100  # heuristic: lineage-cluster bootstrap replicates for every dating interval
CLUSTER_D = 0.2  # bootstrap units: copies joined by pairs closer than this (components never mixed families in sims)
# log-scale sd of model error added to each bootstrap interval: the smallest value giving >= 90% coverage on simulation
# replicates 3-5 (12 scenarios), checked on replicates 6-8; inherited: the smallest giving the intervals' nominal 95%
MODEL_ERR = {"merger": 0.24, "orthologs": 0.015, "confinement": 0.10, "profile": 0.08, "inherited": 0.05}
MU_DEFAULT = 1.3e-8  # substitutions/site/year when --mu is not given (rice LTR rate, Ma & Bennetzen 2004; flagged)
# heuristic: sibling pairs teach the genealogy when their divergence means >= 5 substitutions per LTR
GEN_MIN_SUBS = 5.0
# heuristic: with fewer teaching pairs the genealogy prior stays flat (r uniform: one source lineage)
GEN_MIN_PAIRS = 100
BF_DECISIVE = 100.0  # a call needs a 100-fold likelihood ratio, each lineage counted once ('decisive')
BOOT_MIN = BOOT // 2  # heuristic: a bootstrap interval needs at least half its replicates to give an estimate
DATE_MIN_COPIES = 20  # heuristic: phased copies with K2P needed to sample copy pairs for dating
DATE_MIN_PAIRS = 200  # heuristic: copy pairs needed to date at all, and for one pair of subgenomes
YOUNG_QUANTILE = 0.2  # heuristic: the young stratum, the youngest 20% of phased copies by LTR K2P (FOCAL)
PAIR_MIN_SITES = 500  # a copy pair counts with >= 500 aligned sites of the whole records: whole elements (0.2%)
SCALE_NEAR_FACTOR = 2  # heuristic: the LTR / whole-record ratio near an estimate comes from pairs within 2-fold of it
SCALE_NEAR_MIN = 30  # heuristic: with fewer such pairs (or with LTR sites), the genome-wide ratio instead
SCALE_NEAR_MAX = 2000  # heuristic: pairs aligned at most for that ratio (run time)
RATE_RATIO_MAX = 3000  # heuristic: sibling pairs aligned at most for the genome-wide ratio (run time)
# heuristic: the cross-subgenome profile's bins, log10 pair divergence -4 to 0, 10 a decade
PROF_LOG10_EDGES = (-4, 0, 41)
PROF_MIN_PAIRS = 20  # heuristic: a profile bin's cross fraction is read with >= 20 pairs (binomial se <= 0.112)
TROUGH_POOL = 3  # a trough is sought on 3-bin pooled fractions (a bin and its two neighbors): never one noisy bin
YOUNG_CUT_TOL = 0.02  # heuristic: the young cut is the first bin within 2 points of the trough's cross fraction
PLATEAU_ALPHA = 0.05  # the plateau ends at a bin above it at the conventional 5% (one-sided, Bonferroni over bins)
RISE_MIN = 0.05  # heuristic: the profile split needs the monotone cross fraction to climb >= 5 points
CONF_MIN_COPIES = 200  # heuristic: copies with ancestry evidence and K2P needed for the confinement split
CONF_BINS = 50  # heuristic: equal-count LTR K2P bins of the confinement step fit
# heuristic: its step lies below the 99.5th K2P percentile, the oldest 0.5% (saturated) left out
CONF_TOP_QUANTILE = 0.995
CONF_GRID_POINTS = 400  # heuristic: step positions tried, 400 evenly up to that percentile, with the bin edges
CONF_SE_FLOOR = 1e-5  # heuristic: K2P s.e. floor in the step fit (identical LTRs have s.e. 0: a hard edge)
SIB_MIN = 50  # heuristic: sibling pairs, and copies for their age prior, the merger clock needs
CLOCK_MIN_AGES = 3  # heuristic: the merger clock needs >= 3 candidate ages below its cap
EXCESS_FOLD = 2  # datable: a 2-fold excess of cross-subgenome pairs after the merger (phi >= 2 eps), v7's rule
JC_P_CAP = 0.7  # whole-record p-distances capped at 0.7 so their JC69 distance stays finite (diverges at 0.75)
# heuristic: the clock's LTR divergence grid, 0, then 239 log steps from 2e-6 to 0.5
MERGER_AGE_GRID = (2e-6, 0.5, 239)
# heuristic: the age prior (Efron's g-model) with 8 cubic-spline degrees of freedom and ridge 1
G_MODEL_DF, G_MODEL_RIDGE = 8, 1.0
G_MODEL_OFFSET = 1e-4  # heuristic: its splines run on log(K2P + 1e-4): identical LTRs form one end
FIT_ITERS = 100  # Newton steps at most in the clock's two-rate fit (a concave likelihood: a few suffice)
FIT_GTOL = 1e-8  # heuristic: it stops when the gradient falls below 1e-8 x (1 + |log-likelihood|)
# heuristic: the younger copy's age as a share r of its pair's divergence, in 40 steps
GENEALOGY_R_GRID = (0.025, 1.0, 40)
GENEALOGY_G_GRID = (1e-6, 0.3, 79)  # heuristic: the pair divergence grid (LTR units), 0, then 79 log steps
# heuristic: the older copy's divergence grid of a pair's likelihood, 120 log steps
GENEALOGY_D_GRID = (1e-5, 0.6, 120)
GENEALOGY_ITERS = 200  # heuristic: EM rounds that learn the sibling genealogy f(r)
GENEALOGY_CHUNK = 4000  # pairs per block of the genealogy clock's arrays: bounds memory
# derived biology
# heuristic: copies of a subgenome (or of one LTR-RT clade there) for an insertion history (g-model, 8 df)
MIN_HIST = 200
# heuristic: insertion_history.tsv's LTR K2P bins, 0, then 40 log bins from 1e-4 to 0.5
HIST_K2P_GRID = (1e-4, 0.5, 41)
# painting
BF_SUPPORT = 3.0  # a copy supports an origin when its own evidence favors it >= 3-fold ('substantial')
MIN_LINEAGES = 2  # an exchange needs >= 2 independent lineages (one lineage = one piece of evidence)
START_OTHER = 0.01  # heuristic: prior probability that a sequence starts in another subgenome's ancestry
AGE_BINS = 50  # heuristic: equal-count LTR K2P bins for the age weights of painting evidence (clock bins)
PAINT_RES_WIN = 100_000  # resolution track window (bp): a fraction of the shortest LTR-RT exchanges reported (>= 200 kb)
BREADTH_MIN_SEQS = 2  # heuristic: a clade's breadth scores when one group has >= 2 other sequences carrying it
BREADTH_FOLD = 3.0  # heuristic: ... and the groups' breadths differ >= 3-fold
PAINT_FIT_MAX = 200000  # heuristic: copies drawn at most to calibrate the ancestry experts (run time)
AGE_WEIGHT_MIN = 100  # heuristic: copies with ancestry and K2P needed to weigh evidence by age (else all weigh 1)
AGE_CONFINED_MIN = 0.05  # heuristic: age weights apply when some age class is >= 5 points more confined
# heuristic: Baum-Welch's starting switch probability per copy step; the value without phased copies
SWITCH_START = 1e-3
SWITCH_ITERS = 30  # heuristic: Baum-Welch rounds at most
SWITCH_RANGE = (1e-7, 0.2)  # heuristic: the switch probability kept within these (at most one switch per 5 copies)
SWITCH_TOL = 0.01  # heuristic: it stops when the switch probability moves by <= 1% of itself
# report.pdf (display only; none is a data threshold)
FIG_W = 7.09  # figure width (inches): 180 mm, a two-column journal page
FIG_H = 9.45  # figure height (inches): 240 mm (FIG_W x FIG_H fits an A4 page)
FIG_LET = 0.30  # panel letters sit this far (inches) left of each panel's upper-left corner
FIG_MAX_SEQ = 60  # figures draw the 60 longest phased sequences (more cannot be labeled; tables hold all)
FIG_HALVES_MAX = 300  # panel C: half-library link matrices kept up to this many phased sequences
SG_COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9",  # subgenome hues: Okabe-Ito
             "#332288", "#117733", "#882255", "#999933", "#AA4499", "#44AA99"]  # (no yellow), then Tol muted
FIG_GRAY, FIG_LIGHT, FIG_DARK = "#9a9a9a", "#d9d9d9", "#4d4d4d"  # context: copies without the story's property
FIG_EXCH = "#AA3377"  # near-identical inherited pairs (exchanged after the split); outside the subgenome hues
FIG_ACCENT = "#000000"  # the story's entity in the dating figure (pairs across subgenomes, inherited insertions)
FIG_ALPHA = (0.22, 0.45, 0.72, 1.0)  # painting opacity by the posterior of the painted subgenome: <.7, <.9, <.99, more
FIG_SPLIT, FIG_MERGE, FIG_INK = "#e9e9e9", "#f7edd2", "#222222"  # split band, merger band, text and marks
FIG_EST = {"inherited": "inherited", "consensus": "consensus", "orthologs": "orthologs", "confinement": "confinement",
           "profile": "profile", "ltr_clock": "LTR clock", "genealogy_clock": "genealogy clock"}  # dating.tsv's names
FIG_NUM = dict(enumerate(("no", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten",
                          "eleven", "twelve")))  # counts written as words in panel headers
FIG_XL, FIG_XR = 0.22, 3.12  # panel letters of the left and right columns (inches from the left edge)
# matplotlib rc of report.pdf, the subgenomes page (its sans family added from plot.style.setup_fonts): 7-pt text,
# 6-6.5 pt labels and ticks, thin axes and ticks, no top or right spine, text kept as text in PDF (fonttype 42),
# rasters inside the PDF at FIG_DPI; the whole page is drawn and saved under it (the overviews under RC_BASE, above)
FIG_RC = {
    "font.family": "sans-serif", "font.size": 7, "axes.labelsize": 6.5, "xtick.labelsize": 6, "ytick.labelsize": 6,
    "legend.fontsize": 6.5, "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5, "pdf.fonttype": 42, "ps.fonttype": 42, "axes.spines.top": False,
    "axes.spines.right": False, "mathtext.default": "regular", "hatch.linewidth": 0.4, "savefig.dpi": FIG_DPI}
