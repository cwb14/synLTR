#!/usr/bin/env python3
"""ltrharvest_phasing9.py (v9) - subgenomes, their homoeologous chromosomes, their ages and their exchanges in a
polyploid genome from its LTR-retrotransposon library alone (no genome sequence, no progenitors, no homoeolog map).
v9 = v8 (v8 with a one-page report.pdf) + homoeolog pairing (9.0) and what the pairs say (9.1); every v8
table is unchanged.

Inputs: LTR-RT FASTA (headers chrom:start-end[#Class/Superfamily/Clade]; also chrom:start..end, name::chrom:start-end
and strand suffixes), an LTR annotation with each copy's 5'-3' LTR divergence (LTRquest GFF3, or a K2P table), and
optionally the genome .fai (sequence lengths).

Phasing (v5/v6): every canonical 21-mer of the library is hashed exactly; k-mers carried by 2-250 copies on >= 2
sequences define lineages (distinct carrier sets = clades). Sequences are linked by clades (1/(m-1) per pair, each clade
once, same-sequence copies never link) and by nearest-relative votes, compared with a degree-corrected null and split
by maximum within-group excess (exact for two groups).

Number of subgenomes (new in v7): progenitor lineages were confined to one subgenome only between the progenitors'
split and their merger, so a subgenome split lives in a window of insertion ages: copies older than the split link
homoeologues (inherited insertions), copies younger than the merger move freely. Pooled over all ages, one window can
mask another (a recent split hidden by older homoeologue links). Splits are therefore searched on all copies and in
each LTR-divergence tercile (young, middle, old). A split is kept by the v5/v6 rule: disjoint half-libraries reproduce
it (adjusted Rand >= 0.8) and, below the top, it is >= 0.5x as strong as the strongest split of its parent group -
now measured on the same copies, pooled or within one tercile. A split kept only within a tercile must be age-bounded
(no reproducible split of the group in some tercile): structure present at every age is a property of the chromosomes
themselves (e.g. lineages shared by a few chromosomes of one genome), not of a progenitor's isolation.
Support = agreement over the half-libraries of the split that placed each sequence.

Merger dating (new in v8). A sibling pair of copies sits in two subgenomes only if its younger copy was inserted after
the merger (v6/v7's clock). v7 timed the younger copy by its own two LTRs (0.5-1 kb: one substitution ~ 20-80 kyr), so
mergers younger than that were never dated. v8 adds the genealogy clock: a sibling pair's own divergence over the
whole records (5-10 kb) bounds its younger copy's age ten times more finely, and how far below that divergence the
younger copy was inserted is learned from the genome's older sibling pairs, whose LTRs are informative (the sibling
genealogy, f(r), r = younger age / pair divergence); the pair divergence is put on the LTR scale by the measured
LTR / whole-record rate ratio. Both clocks use v7's sibling pairs plus more young focal copies (up to 12,000 in all;
LTR-RT-rich libraries contribute more) with v7's bootstrap units. A merger is datable only when the step in
cross-subgenome pairs has a decisive likelihood ratio (>= 100, the tool's convention) besides v7's two-fold excess.
Reported: the genealogy clock's interval, else the LTR clock's, else the tighter of their bounds; an interval whose lower
limit lies below one substitution per LTR is reported as a bound (younger than its upper limit).
Split dating (v8): inherited insertions (orthologs) are reported when their interval overlaps the lineage-floor
profile's (two largely independent estimators agree); otherwise the profile, then confinement, then orthologs
(v7 reported confinement first).
Ancestry painting and exchanges: as v6 (painting now runs inside the dating step, on the same ancestry: same result).
Homoeologs (new in v9). Most copies were inserted after the progenitors split: they are specific to one subgenome and
spread over all its chromosomes, so they mark subgenomes, not homoeologues. Homoeologous chromosomes are told apart by
the insertions their common ancestor already carried: such an old copy and its ortholog are each other's closest
relative across the two subgenomes and sit in the same order along both chromosomes. Anchors = reciprocal closest
relatives (shared rare k-mer lineages) between copies of the old LTR-divergence tercile in different subgenomes. A
chromosome's homoeolog in another subgenome is the chromosome whose anchors with it form the most collinear chain,
called when that chromosome also shares the most anchors relative to expectation and the chain is significant
(permutations, or the exact longest-monotone-run distribution for <= 25 anchors; Bonferroni over the candidates).
Inherited insertions are deleted over time, so old splits leave too few: those chromosomes are 'unresolved', never
forced into pairs; many-to-one homoeology (fusions, fissions) is reported as found.
What the pairs say (new in v9.1; each kept only after passing a benchmark bar fixed in advance): (1) exchanges -
stretches where called homoeologs carry near-identical copies of inherited insertions (one homoeolog's sequence
replaced the other's after the split), with positions on both homoeologs and their age; (2) the half-life of intact
LTR-RTs in each subgenome (share of one subgenome's pre-split copies whose ortholog in the other is still intact, read
against the split divergence of the inherited pairs) and whether the two lineages differ; (3) blocks of a chromosome
collinear with a chromosome other than its homoeolog (translocations) and sequences with several homoeologs; (4)
homoeologs collinear in reverse (assembly orientation). Stated in plain words in findings.txt. Tested and dropped
(module notes, V9.1): the split age from inherited pairs, inversions, and typing v8's exchange calls.
Output: tables (TSV/BED/GFF3/bedGraph/JSON), findings.txt and one PDF (report.pdf, one page).

usage: ltrharvest_phasing9.py --ltr_fasta LTR.fa --annotation ANNOT.gff3 [--fai GENOME.fai] --outdir OUT [-t 8]

Developer notes (environment, benchmark results vs v5/v6/SubPhaser, per-genome notes, dead ends, gotchas, tests to
rebuild, leads): the NOTES FOR FUTURE WORK block at the end of this file.
"""
from __future__ import annotations

import os
import sys


def _threads_from_argv(argv, default=4):
    """-t/--threads read before numba is imported: numba sizes its thread pool once, at import. Accepts the forms
    argparse accepts here (-t 8, -t8, -t=8, --threads 8, --threads=8; abbreviations are switched off in parse_args)."""
    for i, x in enumerate(argv):
        if x in ("-t", "--threads") and i + 1 < len(argv) and argv[i + 1].isdigit():
            return int(argv[i + 1])
        for pre in ("--threads=", "-t=", "-t"):
            if x.startswith(pre) and x[len(pre):].isdigit():
                return int(x[len(pre):])
    return default


_THREADS0 = max(1, _threads_from_argv(sys.argv[1:]))
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):   # BLAS single-threaded whatever the
    os.environ[_v] = "1"                                                       # shell exports: -t is the true total
os.environ["NUMBA_NUM_THREADS"] = str(min(_THREADS0, os.cpu_count() or 1))     # numba's pool = -t (v6 used every core)

import argparse
import gzip
import hashlib
import json
import math
import re
import time
from collections import defaultdict
from pathlib import Path
from statistics import NormalDist

import numpy as np

VERSION = "9.1.0"
K = 21               # k-mer length: unique in plant genomes yet short enough to span one SNP per clade variant (v4/v5)
CAP = 250            # largest clade (carrier copies) used for links (v5 chose this cap in 42 of 49 runs)
NMAX_KNN = 50        # relatives for nearest-relative votes come from clades of <= 50 copies (v5 benchmark)
K_NEIGH = 3          # each copy votes for the 3 other sequences holding its closest relatives (v5 benchmark)
MIN_ELEMENTS = 10    # sequences with fewer LTR-RTs are not phased (cannot be assessed; v4/v5)
MIN_REP = 0.8        # a split is kept when half-libraries reproduce it (adjusted Rand; prediction-strength convention)
MIN_REL = 0.5        # a sub-split must be >= this fraction as strong as its parent split on the same copies (v4/v5)
AGE_STRATA = ("young", "middle", "old")   # equal-count LTR-K2P terciles searched for splits: the fewest that can show
#                      a signal confined to a window of insertion ages (absent at one end of the age range)
HALF_PAIRS = 50      # complementary half-library pairs (100 halves): auto-N agreement, support, partners
SUPPORT_OK = 0.9     # support at which a sequence counts as assigned (v4 benchmark: >= 97% correct)
MEM_GB = 8.0         # k-mer extraction holds at most this much (hash, copy) data per pass
KMER_BUDGET = 1.5e8  # k-mers kept: exact up to 150 Mb of LTR-RTs, else 1/ceil(size/150 Mb) (bounds members ~40 M)
MAX_N = 12           # largest number of groups searched (palette and partition search limits)
CHAIN_W = 25         # collinear anchors at most 25 copies apart on both sequences (rank space; tolerates local loss)
NULL_PERMS = 50      # partner-position permutations per sequence pair for the null chain length
MAX_PAIRS = 200      # half-library pairs drawn at most when replicability straddles MIN_REP (precision, not leniency)
PAIR_SPAN = 0.5      # paired sequences share insertions along most of both (homologues); local blocks (segmental
#                      duplications, exchanges) span a small part of them
BAND = 32            # affine realignment within +-32 columns of the edlib path (indel runs in LTR-RT pairs are short)
SCORES = (2, -3, 5, 2)  # match, mismatch, gap open, gap extend: BLASTN log-odds scheme (~90% identity, indels ~1/10 of
#                      substitutions); unit-cost edit distance swaps mismatches for indel pairs (K2P x0.86 at 13.5%)
SKETCH_K = 15        # homolog search k-mer: 15-mers stay shared between copies ~15% apart (old splits); 21-mers do not
SKETCH_S = 128       # bottom-128 MinHash sketch per copy
BUCKET_CAP = 2000    # sketch hashes carried by more copies are family-wide: they do not single out relatives
FOCAL = 3000         # focal copies per stratum (random; youngest 20% by LTR K2P): ~60k copy pairs for the dating
YOUNG_MAX = 12000    # v8: young focal copies in all for the merger clock (v7: FOCAL); more where the library has them
MERGER_PAINT = False  # v8: copies inside called exchange segments count as their origin subgenome in the merger clock
GEN_MIN_SUBS = 5.0   # v8: sibling pairs teach the genealogy when their divergence means >= 5 substitutions per LTR
GEN_MIN_PAIRS = 100  # v8: fewer teaching pairs -> the genealogy prior stays flat (r uniform: one source lineage)
SPLIT_RULE = "agree"  # v8: orthologs when their interval overlaps the profile's (two independent estimators agree),
#                      else profile, confinement, orthologs (SPLIT_ORDER['agree']); v7 reported confinement first
SPLIT_LTR_UNITS = False     # v8: orthologs/profile (pair divergences of whole records) put on the LTR scale of the
#                             confinement estimator and the merger clocks, x the LTR/whole K2P ratio of the pairs
REL_TOP = 10         # relatives per focal copy on other sequences
N_ORTH = 20          # orthologous insertions needed for the ortholog split (fewer: the profile and clock estimates)
BOOT = 100           # lineage-cluster bootstrap replicates for every dating interval
CLUSTER_D = 0.2      # bootstrap units: copies joined by pairs closer than this (components never mixed families in sims)
MODEL_ERR = {"merger": 0.24, "orthologs": 0.015, "confinement": 0.10, "profile": 0.08}
#                      log-scale sd of model error added to each bootstrap interval: the smallest value giving >= 90%
#                      coverage on simulation replicates 3-5 (12 scenarios), checked on replicates 6-8 (Plan 4)
MU_DEFAULT = 1.3e-8  # substitutions/site/year when --mu is not given (rice LTR rate, Ma & Bennetzen 2004; flagged)
AGE_GRID = np.concatenate(([0.0], np.geomspace(2e-6, 0.5, 239)))   # LTR divergence grid of the merger clock
MIN_HIST = 200       # copies of a subgenome (or of one LTR-RT clade there) for an insertion history (g-model, 8 df)
HIST_EDGES = np.concatenate(([0.0], np.geomspace(1e-4, 0.5, 41)))   # LTR K2P bins of insertion_history.tsv: 40
#                      log-spaced over 1e-4..0.5 (~11 per decade) plus one for identical LTRs
BF_DECISIVE = 100.0  # an exchange call needs a 100-fold likelihood ratio, each lineage counted once ('decisive', v5)
BF_SUPPORT = 3.0     # a copy supports an origin when its own evidence favours it >= 3-fold ('substantial', v5)
MIN_LINEAGES = 2     # an exchange needs >= 2 independent lineages (one lineage = one piece of evidence, v5)
START_OTHER = 0.01   # prior probability that a sequence starts in another subgenome's ancestry (v5)
AGE_BINS = 50        # equal-count LTR K2P bins for the age weights of painting evidence (v5 clock bins)
RES_WIN = 100_000    # resolution track window (bp): a fraction of the shortest LTR-RT exchanges reported (>= 200 kb)
FIG_W = 7.09         # figure width (inches): 180 mm, a two-column journal page
FIG_LET = 0.30       # panel letters sit this far (inches) left of each panel's upper-left corner
FIG_MAX_SEQ = 60     # figures draw the 60 longest phased sequences (more cannot be labelled; tables hold all)
FIG_HALVES_MAX = 300  # report.pdf panel C: half-library link matrices kept up to this many phased sequences
SG_COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9",   # subgenome hues: Okabe-Ito
             "#332288", "#117733", "#882255", "#999933", "#AA4499", "#44AA99"]  # (no yellow), then Tol muted
FIG_GRAY, FIG_LIGHT, FIG_DARK = "#9a9a9a", "#d9d9d9", "#4d4d4d"   # context: copies without the story's property
FIG_EXCH = "#AA3377"   # v9.1: near-identical inherited pairs (exchanged after the split); outside the subgenome hues
FIG_ACCENT = "#000000"  # the story's entity in the dating figure (pairs across subgenomes, inherited insertions)
T0 = time.time()
VERBOSE = False


def log(msg):
    print(f"[v9 {time.time() - T0:7.1f}s] {msg}", file=sys.stderr, flush=True)


def vlog(msg):
    if VERBOSE:
        log(msg)


def die(msg):
    log("ERROR: " + msg)
    sys.exit(1)


def natural_key(s):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", str(s))]


def parse_args(argv=None):
    ap = argparse.ArgumentParser(allow_abbrev=False,
                                 description="Subgenomes, their ages and their exchanges from an LTR-RT library alone. "
                                 "Writes tables and one PDF (report.pdf).")
    ap.add_argument("--ltr_fasta", nargs="+", required=True, help="LTR-RT FASTA(.gz) file(s)")
    ap.add_argument("--annotation", nargs="+", required=True,
                    help="LTR annotation with each copy's LTR divergence: LTRquest GFF3 or K2P table (TSV)")
    ap.add_argument("--fai", nargs="+", default=None, help="genome .fai (sequence lengths; optional)")
    ap.add_argument("--outdir", required=True, help="output directory")
    ap.add_argument("-n", "--n_subgenomes", type=int, default=None, help="number of subgenomes (default: inferred)")
    ap.add_argument("--assign", default=None, help="fixed partition instead of phasing: 'sequence label' per line (whitespace-separated; "
                    "labels become SG1, SG2... in file order, mapping in summary.json)")
    ap.add_argument("--mu", type=float, default=None, help="substitutions/site/year for ages in years (default 1.3e-8, flagged)")
    ap.add_argument("--label", default=None, help="species name for the PDF (default: FASTA name)")
    ap.add_argument("-t", "--threads", type=int, default=4, help="total threads")
    ap.add_argument("--seed", type=int, default=1, help="random seed")
    ap.add_argument("-v", "--verbose", action="store_true", help="per-step progress and sanity checks")
    ap.add_argument("--dev_scale", type=int, default=0, help=argparse.SUPPRESS)      # k-mer sampling 1/scale (0 = budget)
    ap.add_argument("--dev_skip_dating", action="store_true", help=argparse.SUPPRESS)  # tests of phasing only
    return ap.parse_args(argv)


def setup(a):
    global VERBOSE
    VERBOSE = a.verbose
    for p in list(a.ltr_fasta) + list(a.annotation) + list(a.fai or []):
        if not Path(p).exists():
            die(f"input not found: {p}")
    if a.mu is not None and a.mu <= 0:
        die("--mu must be > 0 (substitutions per site per year)")
    if a.threads < 1:
        die("--threads must be >= 1")
    if a.n_subgenomes is not None and not 2 <= a.n_subgenomes <= MAX_N:
        die(f"-n must be between 2 and {MAX_N} (omit it to infer the structure)")
    if a.assign:
        if not Path(a.assign).exists():
            die(f"input not found: {a.assign}")
        if a.n_subgenomes:
            die("--assign fixes the partition: -n does not apply")
    import numba
    if a.threads > numba.config.NUMBA_NUM_THREADS:
        log(f"WARNING: -t {a.threads} exceeds numba's pool of {numba.config.NUMBA_NUM_THREADS} threads (sized when v9 "
            f"was imported); numba uses {numba.config.NUMBA_NUM_THREADS}, process pools use {a.threads}")
    numba.set_num_threads(min(a.threads, numba.config.NUMBA_NUM_THREADS))
    Path(a.outdir).mkdir(parents=True, exist_ok=True)
    (Path(a.outdir) / ".v9cache").mkdir(exist_ok=True)


SOURCE_SHA = hashlib.sha1(Path(__file__).read_bytes()).hexdigest()[:12]   # code changes invalidate caches


def cache_path(a, stage, files, params):
    """Cache file keyed by version and source, stage, input files (path, size, mtime) and parameters."""
    h = hashlib.sha1(f"{VERSION}|{SOURCE_SHA}|{stage}|{params}".encode())
    for f in files:
        st = os.stat(f)
        h.update(f"{os.path.abspath(f)}|{st.st_size}|{int(st.st_mtime)}".encode())
    return Path(a.outdir) / ".v9cache" / f"{stage}_{h.hexdigest()[:16]}.npz"


# ------------------------------------------------------------------------------------------------ inputs
LUT = np.full(256, 4, np.uint8)
for _i, _b in enumerate(b"ACGT"):
    LUT[_b] = _i
    LUT[_b + 32] = _i
MIN_MATCH = 0.95     # share of FASTA copies that must be annotated (fail fast otherwise)


def parse_header(tok):
    """'chrom:start-end[#class]' (LTRquest/EDTA), 'chrom:start..end', 'name::chrom:start-end' (bedtools -name), with an
    optional '(+)'/'(-)' suffix -> (chrom, start, end, class) or None."""
    core, _, cls = tok.partition("#")
    if core.endswith(("(+)", "(-)")):
        core = core[:-3]
    if "::" in core:
        core = core.rsplit("::", 1)[1]
    chrom, sep, rng = core.rpartition(":")
    if not sep or not chrom:
        return None
    a, _, b = rng.partition("..") if ".." in rng else rng.partition("-")
    if not (a.isdigit() and b.isdigit()):
        return None
    return chrom, int(a), int(b), cls


def _open(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def read_fasta(path):
    """Yield (first header token, sequence bytes)."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rb") as f:
        name, buf = None, []
        for line in f:
            if line[:1] == b">":
                if name is not None:
                    yield name, b"".join(buf)
                tok = line[1:].split()
                name, buf = (tok[0].decode() if tok else ""), []
            else:
                buf.append(line.strip())
        if name is not None:
            yield name, b"".join(buf)


def _num(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else float("nan")
    except (TypeError, ValueError):
        return float("nan")


def _tsd(x):
    """Target site duplication as annotated (upper case); '' when absent ('.', 'NA', empty)."""
    x = (x or "").strip().upper()
    return x if x and x not in (".", "NA") and set(x) <= set("ACGTN") else ""


def _read_gff(path, out):
    elems, ltrs = {}, {}
    with _open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            at = dict(x.split("=", 1) for x in c[8].split(";") if "=" in x)
            if c[2] in ("LTR_retrotransposon", "repeat_region") and ("K2P_d" in at or "k2p" in at or "K2P" in at):
                key = (c[0], int(c[3]), int(c[4]))
                elems[at.get("ID", f"{key}")] = key
                k2p = _num(at.get("K2P_d", at.get("k2p", at.get("K2P"))))
                n_sites = _num(at.get("n_sites", at.get("aln_len")))
                out[key] = dict(k2p=k2p, k2p_se=_num(at.get("k2p_se")), n_sites=n_sites, n_ts=_num(at.get("ti")),
                                n_tv=_num(at.get("tv")), ltr5_len=_num(at.get("ltr5_len")), ltr3_len=_num(at.get("ltr3_len")),
                                tsd=_tsd(at.get("tsd")), family=at.get("family", ""), clade=at.get("clade", ""),
                                superfamily=at.get("superfamily", ""), depth=_num(at.get("depth")),
                                strand=c[6] if c[6] in "+-" else ".")
            elif c[2] == "long_terminal_repeat" and "Parent" in at:
                ltrs.setdefault(at["Parent"], []).append((int(c[3]), int(c[4])))
    for pid, iv in ltrs.items():
        key = elems.get(pid)
        if key in out and len(iv) == 2:
            out[key]["ltr_g"] = sorted(iv)                 # genome intervals (1-based, inclusive)


def _read_table(path, out):
    with _open(path) as f:
        hdr = None
        for line in f:
            if not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if line.startswith("#") or hdr is None:
                if hdr is None:
                    hdr = [h.lstrip("#").strip().lower() for h in p]
                    if "k2p" not in hdr:
                        die(f"{path}: no 'k2p' column in the header (columns found by name)")
                continue
            r = dict(zip(hdr, p))
            x = parse_header(r.get("seq_id", p[0]))
            if x is None:
                continue
            key = x[:3]
            n_sites = _num(r.get("n_sites", r.get("ltr5_len")))
            out[key] = dict(k2p=_num(r["k2p"]), k2p_se=_num(r.get("k2p_se")), n_sites=n_sites, n_ts=_num(r.get("n_ts")),
                            n_tv=_num(r.get("n_tv")), ltr5_len=_num(r.get("ltr5_len")), ltr3_len=_num(r.get("ltr3_len")),
                            tsd=_tsd(r.get("tsd")), family=r.get("family", ""), clade="", superfamily="",
                            depth=float("nan"), strand=r.get("strand", ".") if r.get("strand", ".") in ("+", "-") else ".")
            if all(r.get(c, "").isdigit() for c in ("ltr5_start", "ltr5_end", "ltr3_start", "ltr3_end")):
                out[key]["ltr_r"] = ((int(r["ltr5_start"]) - 1, int(r["ltr5_end"])),
                                     (int(r["ltr3_start"]) - 1, int(r["ltr3_end"])))   # record offsets, half-open


def read_annotation(paths):
    out = {}
    for p in paths:
        with _open(p) as f:
            first = next((l for l in f if l.strip()), "")
        if first.startswith("##gff") or (len(first.split("\t")) == 9 and not first.lower().lstrip("#").startswith("seq_id")):
            _read_gff(p, out)
        else:
            _read_table(p, out)
    return out


def read_fai(paths):
    out = {}
    for p in paths:
        with open(p) as f:
            for line in f:
                c = line.split("\t")
                if len(c) >= 2:
                    out[c[0]] = int(c[1])
    return out


class Library:
    """LTR-RT copies: sequence, coordinates, 2-bit codes (concatenated, offs = copy boundaries) and annotation."""


def record_offsets(ivs, start, end, strand):
    """Genome intervals (1-based, inclusive) of the two LTRs of the element at start..end -> half-open offsets in its
    FASTA record, 5' LTR first. Minus-strand records are reverse complements (LTRquest writes element orientation);
    unknown-strand records are in genome orientation."""
    if strand == "-":
        out = sorted((end - b, end - a + 1) for a, b in ivs)
    else:
        out = sorted((a - start, b - start + 1) for a, b in ivs)
    return out[0], out[1]


def load_library(a):
    ids, chroms, starts, ends, cls, chunks = [], [], [], [], [], []
    bad = 0
    for p in a.ltr_fasta:
        for name, s in read_fasta(p):
            x = parse_header(name)
            if x is None:
                bad += 1
                continue
            ids.append(f"{x[0]}:{x[1]}-{x[2]}")
            chroms.append(x[0])
            starts.append(x[1])
            ends.append(x[2])
            cls.append(x[3])
            chunks.append(LUT[np.frombuffer(s, np.uint8)])
    if bad:
        log(f"WARNING: skipped {bad} FASTA records whose header has no chrom:start-end")
    if not ids:
        die("no usable LTR-RT records in " + ", ".join(map(str, a.ltr_fasta)))
    first = {}
    for i, x in enumerate(ids):
        first.setdefault(x, i)
    n_dup = len(ids) - len(first)
    if n_dup > (1 - MIN_MATCH) * len(ids):
        die(f"{n_dup} of {len(ids)} FASTA records repeat the coordinates of an earlier record")
    if n_dup:
        log(f"WARNING: {n_dup} FASTA records repeat the coordinates of an earlier record; the first is kept")
        keep = sorted(first.values())
        ids, chroms, starts, ends, cls, chunks = ([v[i] for i in keep] for v in (ids, chroms, starts, ends, cls, chunks))
    ann = read_annotation(a.annotation)
    keys = list(zip(chroms, starts, ends))
    hit = [k in ann for k in keys]
    if sum(hit) < MIN_MATCH * len(keys):
        miss = [ids[i] for i, h in enumerate(hit) if not h][:5]
        die(f"only {sum(hit)} of {len(keys)} FASTA copies are annotated (need {MIN_MATCH:.0%}); unmatched e.g. "
            + ", ".join(miss) + " (annotation and FASTA must share chrom:start-end coordinates)")
    lib = Library()
    lib.n, lib.ids, lib.cls = len(ids), ids, cls
    lib.start, lib.end = np.array(starts, np.int64), np.array(ends, np.int64)
    lens = np.array([c.size for c in chunks], np.int64)
    lib.offs = np.concatenate(([0], np.cumsum(lens)))
    lib.codes = np.concatenate(chunks) if chunks else np.empty(0, np.uint8)
    get = lambda f: np.array([ann[k][f] if h else float("nan") for k, h in zip(keys, hit)], float)
    for f in ("k2p", "k2p_se", "n_sites", "n_ts", "n_tv", "ltr5_len", "ltr3_len", "depth"):
        setattr(lib, f, get(f))
    for f in ("family", "clade", "superfamily", "strand"):
        setattr(lib, f, [ann[k][f] if h else "" for k, h in zip(keys, hit)])
    lib.tsd = np.array([ann[k].get("tsd", "") if h else "" for k, h in zip(keys, hit)])
    lib.strand_code = np.array([{"+": 1, "-": -1}.get(x, 0) for x in lib.strand], np.int8)
    lib.l5, lib.l3 = np.full((lib.n, 2), -1, np.int64), np.full((lib.n, 2), -1, np.int64)   # -1: no LTR bounds
    n_bad = 0
    for i, (k, h) in enumerate(zip(keys, hit)):
        x = ann[k] if h else {}
        if "ltr_r" in x:
            o5, o3 = x["ltr_r"]
        elif "ltr_g" in x:
            o5, o3 = record_offsets(x["ltr_g"], k[1], k[2], x["strand"])
        else:
            continue
        L = int(lens[i])
        if 0 <= o5[0] < o5[1] <= o3[0] < o3[1] <= L:
            lib.l5[i], lib.l3[i] = o5, o3
        else:
            n_bad += 1
    if n_bad:
        log(f"WARNING: {n_bad} copies have LTR bounds outside their record; whole elements used for them")
    ok = lib.k2p[np.isfinite(lib.k2p)]
    if ok.size and (np.mean(ok >= 0.5) > 0.5 or np.any(ok < 0)):
        die("the divergence column looks like identity or is negative (K2P distances are expected, mostly < 0.3)")
    if np.isfinite(lib.k2p).sum() < MIN_MATCH * lib.n:
        die(f"K2P missing for {lib.n - np.isfinite(lib.k2p).sum()} of {lib.n} copies")
    lib.fai_given = bool(a.fai)
    if a.fai:
        fai = read_fai(a.fai)
        absent = sorted(set(chroms) - set(fai), key=natural_key)
        if absent:
            die(f"{len(absent)} sequence(s) carrying LTR-RTs are absent from the .fai, e.g. " + ", ".join(absent[:5]))
        lib.seqnames = sorted(fai, key=natural_key)
        lib.seq_len = np.array([fai[s] for s in lib.seqnames], np.int64)
    else:
        lib.seqnames = sorted(set(chroms), key=natural_key)
        last = {}
        for c, e in zip(chroms, ends):
            last[c] = max(last.get(c, 0), e)
        lib.seq_len = np.array([last[s] for s in lib.seqnames], np.int64)
        log("no --fai: sequence lengths = end of their last LTR-RT; LTR-free sequences are not listed")
    six = {s: i for i, s in enumerate(lib.seqnames)}
    lib.seq = np.array([six[c] for c in chroms], np.int32)
    log(f"{lib.n:,} LTR-RTs on {len(set(chroms))} sequences; {int(np.isfinite(lib.k2p).sum()):,} with K2P")
    return lib


# ------------------------------------------------------------------------------------------------ k-mer lineages
import numba  # noqa: E402
from numba import njit, prange  # noqa: E402

numba.set_num_threads(max(1, min(_THREADS0, numba.config.NUMBA_NUM_THREADS)))   # until setup() applies -t

RUN_LIM = 12500      # k-mer runs longer than this (occurrences) are skipped: they cannot have <= CAP carriers
#                      unless repeated tens of times inside each copy (tandem arrays), which carry no lineage signal


@njit(cache=True, inline="always")
def _mix64(x):
    x = x ^ (x >> np.uint64(33))
    x = x * np.uint64(0xFF51AFD7ED558CCD)
    x = x ^ (x >> np.uint64(33))
    x = x * np.uint64(0xC4CEB9FE1A85EC53)
    x = x ^ (x >> np.uint64(33))
    return x


@njit(parallel=True, cache=True)
def _kmer_pass(codes, offs, k, blo, bhi, cut, fill, cnt_or_off, H, E):
    """Canonical k-mer hashes of every copy whose top 16 bits lie in [blo, bhi) and that are <= cut (FracMinHash;
    cut = 2^64-1 keeps all). fill=False: count per copy into cnt_or_off; fill=True: write (hash, copy) from
    cnt_or_off[copy] on."""
    n = offs.size - 1
    mask = (np.uint64(1) << np.uint64(2 * k)) - np.uint64(1)
    sh = np.uint64(2 * (k - 1))
    for i in prange(n):
        f = np.uint64(0)
        r = np.uint64(0)
        run = 0
        o = cnt_or_off[i] if fill else 0
        cnt = 0
        for p in range(offs[i], offs[i + 1]):
            c = codes[p]
            if c > 3:
                run = 0
                f = np.uint64(0)
                r = np.uint64(0)
                continue
            cu = np.uint64(c)
            f = ((f << np.uint64(2)) | cu) & mask
            r = (r >> np.uint64(2)) | ((np.uint64(3) - cu) << sh)
            run += 1
            if run >= k:
                h = _mix64(f if f < r else r)
                b = np.int64(h >> np.uint64(48))
                if h <= cut and b >= blo and b < bhi:
                    if fill:
                        H[o] = h
                        E[o] = i
                        o += 1
                    else:
                        cnt += 1
        if not fill:
            cnt_or_off[i] = cnt


@njit(cache=True)
def _bucket_sort(H, E, blo, nb):
    """(hash, copy) records moved into top-16-bit bucket order (counting sort), and the bucket bounds: each bucket
    is then a contiguous slice that sorts in cache."""
    bounds = np.zeros(nb + 1, np.int64)
    for j in range(H.size):
        bounds[np.int64(H[j] >> np.uint64(48)) - blo + 1] += 1
    for b in range(nb):
        bounds[b + 1] += bounds[b]
    pos = bounds[:-1].copy()
    Hs = np.empty_like(H)
    Es = np.empty_like(E)
    for j in range(H.size):
        b = np.int64(H[j] >> np.uint64(48)) - blo
        Hs[pos[b]] = H[j]
        Es[pos[b]] = E[j]
        pos[b] += 1
    return Hs, Es, bounds


@njit(cache=True)
def _distinct(vals, m, out):
    """Sorted distinct values of vals[:m] written to out; returns their number. In place, no allocation: runs are
    short, so insertion sort (array sort above 32 values)."""
    for i in range(m):
        out[i] = vals[i]
    if m > 32:
        out[:m].sort()
    else:
        for i in range(1, m):
            x = out[i]
            j = i - 1
            while j >= 0 and out[j] > x:
                out[j + 1] = out[j]
                j -= 1
            out[j + 1] = x
    nd = 0
    for t in range(m):
        if t == 0 or out[t] != out[t - 1]:
            out[nd] = out[t]
            nd += 1
    return nd


@njit(parallel=True, cache=True)
def _scan_buckets(H, E, bounds, seqof, rnd, cap, write, nq, nm, qoff, moff, sig, size, nseq, mem, nchunk):
    """Per bucket: sort by hash; runs of one k-mer carried by 2..cap distinct copies on >= 2 sequences qualify.
    write=False: count qualifying runs (nq) and members (nm) per bucket; write=True: write signature (xor of
    per-copy random keys ^ mix(size)), size, n sequences and members at the bucket offsets. Buckets are processed
    in nchunk parallel chunks, each with its own scratch buffers (allocated once per chunk, not per bucket)."""
    nb = bounds.size - 1
    for ch in prange(nchunk):
        cps = np.empty(RUN_LIM, np.int64)
        buf = np.empty(RUN_LIM, np.int64)
        sq = np.empty(RUN_LIM, np.int64)
        sbuf = np.empty(RUN_LIM, np.int64)
        for b in range(nb * ch // nchunk, nb * (ch + 1) // nchunk):
            s, e = bounds[b], bounds[b + 1]
            q = 0
            mm = 0
            if e - s >= 2:
                seg = np.argsort(H[s:e]) + s          # bucket is contiguous (bucket order)
                qo = qoff[b] if write else 0
                mo = moff[b] if write else 0
                L = e - s
                j = 0
                while j < L:
                    t = j + 1
                    hj = H[seg[j]]
                    while t < L and H[seg[t]] == hj:
                        t += 1
                    m = t - j
                    if 2 <= m <= RUN_LIM:
                        for u in range(m):
                            cps[u] = E[seg[j + u]]
                        nd = _distinct(cps, m, buf)
                        if 2 <= nd <= cap:
                            for u in range(nd):
                                sq[u] = seqof[buf[u]]
                            ns = _distinct(sq, nd, sbuf)
                            if ns >= 2:
                                if write:
                                    sg = _mix64(np.uint64(nd))
                                    for u in range(nd):
                                        sg = sg ^ rnd[buf[u]]
                                    sig[qo] = sg
                                    size[qo] = nd
                                    nseq[qo] = ns
                                    for u in range(nd):
                                        mem[mo + u] = buf[u]
                                    qo += 1
                                    mo += nd
                                else:
                                    q += 1
                                    mm += nd
                    j = t
            if not write:
                nq[b] = q
                nm[b] = mm


class Clades:
    """Distinct carrier sets of shared k-mers (2..cap copies on >= 2 sequences). CSR: members of clade i are
    members[indptr[i]:indptr[i+1]]; mult = k-mers with exactly that carrier set; nseq = sequences spanned."""

    def __init__(self, members, size, mult, nseq, n_occ=0, n_multi=0):
        self.members = members.astype(np.int32)
        self.size = size.astype(np.int32)
        self.mult = mult.astype(np.int64)
        self.nseq = nseq.astype(np.int32)
        self.indptr = np.concatenate(([0], np.cumsum(self.size, dtype=np.int64)))
        self.n = int(self.size.size)
        self.n_occ, self.n_multi = int(n_occ), int(n_multi)

    def sid(self):
        return np.repeat(np.arange(self.n), self.size)

    def restrict(self, copy_mask, seq=None):
        """Members inside copy_mask; clades keep >= 2 members (on >= 2 sequences when seq is given). Clades whose
        kept carrier sets coincide merge (mult summed): a lineage reaching removed copies is not counted twice."""
        keep_m = copy_mask[self.members]
        sid = self.sid()[keep_m]
        mem = self.members[keep_m]
        sz = np.bincount(sid, minlength=self.n)
        ok = sz >= 2
        if seq is not None and mem.size:
            u = np.unique(sid.astype(np.int64) * (int(seq.max()) + 1) + seq[mem])
            ns = np.bincount(u // (int(seq.max()) + 1), minlength=self.n)
            ok &= ns >= 2
        else:
            ns = self.nseq.copy()
        km = ok[sid]
        mem, size, mult, nseq = mem[km], sz[ok], self.mult[ok], (ns if seq is not None else self.nseq)[ok]
        if size.size and not copy_mask.all():                    # sets differing only by removed copies coincide now
            sig, size, nseq, mult, mem = _dedupe(_set_signatures(mem, size), size, nseq, mult, mem)
        return Clades(mem, size, mult, nseq, self.n_occ, self.n_multi)


def build_clades(lib, k=K, cap=CAP, scale=1, mem_gb=MEM_GB):
    n = lib.n
    rng = np.random.default_rng(20261012)
    rnd = (rng.integers(1, 2 ** 62, n, dtype=np.int64).astype(np.uint64) << np.uint64(1)) | np.uint64(1)
    cut = np.uint64(2 ** 64 - 1) if scale <= 1 else np.uint64((1 << 64) // scale - 1)
    eH, eE = np.empty(0, np.uint64), np.empty(0, np.int32)
    counts = np.zeros(n, np.int64)
    _kmer_pass(lib.codes, lib.offs, k, 0, 65536, cut, False, counts, eH, eE)
    total = int(counts.sum())
    per_pass = max(1, int(mem_gb * 1e9 / 32))
    n_pass = min(65536, max(1, math.ceil(total / per_pass)))
    edges = np.linspace(0, 65536, n_pass + 1).astype(np.int64)
    seqof = lib.seq.astype(np.int64)
    P_sig, P_size, P_nseq, P_mult, P_mem = [], [], [], [], []
    n_multi = 0
    for p in range(n_pass):
        blo, bhi = int(edges[p]), int(edges[p + 1])
        if bhi <= blo:
            continue
        cnt = np.zeros(n, np.int64)
        _kmer_pass(lib.codes, lib.offs, k, blo, bhi, cut, False, cnt, eH, eE)
        tot = int(cnt.sum())
        if tot == 0:
            continue
        off = np.zeros(n, np.int64)
        off[1:] = np.cumsum(cnt)[:-1]
        H, E = np.empty(tot, np.uint64), np.empty(tot, np.int32)
        _kmer_pass(lib.codes, lib.offs, k, blo, bhi, cut, True, off, H, E)
        H, E, bounds = _bucket_sort(H, E, blo, bhi - blo)
        nb = bhi - blo
        nq, nm = np.zeros(nb, np.int64), np.zeros(nb, np.int64)
        z64, zu, z32 = np.zeros(1, np.int64), np.zeros(1, np.uint64), np.zeros(1, np.int32)
        nchunk = max(1, min(nb, 8 * numba.get_num_threads()))
        _scan_buckets(H, E, bounds, seqof, rnd, cap, False, nq, nm, z64, z64, zu, z32, z32, z64, nchunk)
        qoff = np.concatenate(([0], np.cumsum(nq)[:-1]))
        moff = np.concatenate(([0], np.cumsum(nm)[:-1]))
        Q, M = int(nq.sum()), int(nm.sum())
        n_multi += Q
        if Q == 0:
            continue
        sig, size, nseq, mem = np.empty(Q, np.uint64), np.empty(Q, np.int32), np.empty(Q, np.int32), np.empty(M, np.int64)
        _scan_buckets(H, E, bounds, seqof, rnd, cap, True, nq, nm, qoff, moff, sig, size, nseq, mem, nchunk)
        del H, E
        sig, size, nseq, mult, mem = _dedupe(sig, size, nseq, np.ones(Q, np.int64), mem)
        P_sig.append(sig)
        P_size.append(size)
        P_nseq.append(nseq)
        P_mult.append(mult)
        P_mem.append(mem)
    if not P_sig:
        return Clades(np.empty(0, np.int32), np.empty(0, np.int32), np.empty(0, np.int64), np.empty(0, np.int32), total, 0)
    sig, size, nseq, mult, members = _dedupe(*[np.concatenate(x) for x in (P_sig, P_size, P_nseq, P_mult, P_mem)])
    return Clades(members, size, mult, nseq, total, n_multi)


def _mix64_np(x):
    x = np.asarray(x, np.uint64)
    x = x ^ (x >> np.uint64(33))
    x = x * np.uint64(0xFF51AFD7ED558CCD)
    x = x ^ (x >> np.uint64(33))
    x = x * np.uint64(0xC4CEB9FE1A85EC53)
    return x ^ (x >> np.uint64(33))


def _set_signatures(mem, size):
    """Order-free signature of each carrier set (xor of per-copy keys ^ mix(size)), clades concatenated in mem."""
    keys = _mix64_np(mem.astype(np.uint64) * np.uint64(0x9E3779B97F4A7C15) + np.uint64(1))
    starts = np.concatenate(([0], np.cumsum(size)[:-1])).astype(np.int64)
    return np.bitwise_xor.reduceat(keys, starts) ^ _mix64_np(size)


def _dedupe(sig, size, nseq, mult, mem):
    """Merge entries with the same carrier-set signature (mult summed). Entries keep the order of their first
    occurrence; mem (members concatenated in entry order) keeps the members of each first occurrence."""
    u, first, inv = np.unique(sig, return_index=True, return_inverse=True)
    mult_u = np.bincount(inv, weights=mult, minlength=u.size).astype(np.int64)
    o = np.argsort(first)
    keep = np.zeros(sig.size, bool)
    keep[first] = True
    return u[o], size[first[o]], nseq[first[o]], mult_u[o], mem[np.repeat(keep, size)]


def auto_scale(n_bp):
    """FracMinHash scale keeping at most KMER_BUDGET k-mers (1 = exact)."""
    return max(1, math.ceil(n_bp / KMER_BUDGET))


def load_or_build_clades(a, lib):
    scale = a.dev_scale or auto_scale(int(lib.codes.size))
    f = cache_path(a, "clades", list(a.ltr_fasta), f"k={K}|cap={CAP}|scale={scale}")
    cl = None
    if f.exists():
        try:
            z = np.load(f)
            cl = Clades(z["members"], z["size"], z["mult"], z["nseq"], int(z["n_occ"]), int(z["n_multi"]))
            vlog(f"clades from cache {f.name}")
        except (OSError, ValueError, KeyError, EOFError) as e:
            log(f"WARNING: unreadable cache {f.name} ({type(e).__name__}); rebuilding")
    if cl is None:
        cl = build_clades(lib, scale=scale)
        tmp = f.with_name(f.stem + ".tmp.npz")
        np.savez(tmp, members=cl.members, size=cl.size, mult=cl.mult, nseq=cl.nseq, n_occ=cl.n_occ, n_multi=cl.n_multi)
        os.replace(tmp, f)                                        # never a half-written cache
    cl.scale = scale
    log(f"{cl.n_occ:,} k-mers" + (f" (1/{scale} sampled: library > {KMER_BUDGET / 1e6:.0f} Mb)" if scale > 1 else
                                   " (exact)") + f"; {cl.n:,} clades of 2-{CAP} copies on >= 2 sequences")
    return cl


# ------------------------------------------------------------------------------------------------ links
from scipy import sparse  # noqa: E402


@njit(parallel=True, cache=True)
def _link_kernel(members, indptr, code, C, nchunk):
    """Per chunk of clades (summed by the caller): each clade adds 1/(m-1) to every pair of the m distinct
    sequences spanned by its members with code >= 0 (code = sequence index of copies in the mask, else -1)."""
    n = indptr.size - 1
    out = np.zeros((nchunk, C, C))
    for ch in prange(nchunk):
        seen = np.full(C, -1, np.int64)
        buf = np.empty(C, np.int64)
        for k in range(n * ch // nchunk, n * (ch + 1) // nchunk):
            m = 0
            for p in range(indptr[k], indptr[k + 1]):
                sq = code[members[p]]
                if sq >= 0 and seen[sq] != k:
                    seen[sq] = k
                    buf[m] = sq
                    m += 1
            if m >= 2:
                w = 1.0 / (m - 1)
                for i in range(m):
                    for j in range(i + 1, m):
                        out[ch, buf[i], buf[j]] += w
                        out[ch, buf[j], buf[i]] += w
    return out


def n_chunks(C, n_items):
    """Chunks of a parallel accumulation into C x C matrices: 64 (fixed, so results do not depend on -t), fewer when
    64 copies of the matrix would exceed 1 GB, never more than the items."""
    return int(max(1, min(64, n_items, 1e9 // max(8 * C * C, 1))))


def link_matrix(cl, ek, C, mask):
    """C x C link matrix from the clades' members in mask: each clade adds 1/(m-1) to every pair of the m sequences
    it spans (a clade counts once however many copies or k-mers it has); same-sequence copies never link; copies
    with ek < 0 are ignored."""
    if cl.n == 0:
        return np.zeros((C, C))
    code = np.where(np.asarray(mask, bool), np.asarray(ek, np.int64), -1)
    nchunk = n_chunks(C, cl.n)
    return _link_kernel(cl.members, cl.indptr, code, C, nchunk).sum(0)


def affinity(cl, ek, C):
    return link_matrix(cl, ek, C, np.ones(len(ek), bool))


_EXPECTED_WARNED = False


def expected(W, iters=200):
    """Degree-corrected null for a zero-diagonal matrix: E_ij = b_i b_j with sum_{j!=i} E_ij = row sum. Warns
    (once per run) when the fixed-point iteration has not converged."""
    global _EXPECTED_WARNED
    d = W.sum(1)
    b = d / math.sqrt(max(d.sum(), 1e-300))
    converged = False
    for _ in range(iters):
        nb = np.sqrt(b * d / np.maximum(b.sum() - b, 1e-300))
        if np.allclose(nb, b, rtol=1e-12, atol=0):
            b, converged = nb, True
            break
        b = nb
    E = np.outer(b, b)
    np.fill_diagonal(E, 0)
    if not converged and not _EXPECTED_WARNED:
        _EXPECTED_WARNED = True
        err = float(np.max(np.abs(E.sum(1) - d)) / max(float(d.max()), 1e-300))
        log(f"WARNING: the degree-corrected null did not converge in {iters} iterations (largest row-sum error "
            f"{err:.2g} of the largest row); the last iterate is used")
    return E


def log_oe(W):
    return np.log((W + 1) / (expected(W) + 1))


def _relatives_block(e01):
    """Closest relatives of copies e0..e1-1 (one block of the sparse self-join); run in pool workers."""
    e0, e1 = e01
    AT, DA, ek, C, T = (_G[x] for x in ("AT", "DA", "ek", "C", "T"))
    rel_j = np.full((e1 - e0, C, T), -1, np.int32)
    rel_r = np.zeros((e1 - e0, C, T), np.float32)
    R = (AT[e0:e1] @ DA).tocoo()
    rows, cols, vals = R.row, R.col, R.data
    ch = ek[cols]
    k = (ch != ek[e0 + rows]) & (ch >= 0)
    rows, cols, ch, vals = rows[k], cols[k], ch[k], vals[k]
    if rows.size:
        key = rows.astype(np.int64) * C + ch
        o = np.lexsort((-vals, key))
        key, cols, vals = key[o], cols[o], vals[o]
        st = np.searchsorted(key, key, side="left")
        rank = np.arange(key.size) - st
        m = rank < T
        r_, c_ = key[m] // C, key[m] % C
        rel_j[r_, c_, rank[m]] = cols[m]
        rel_r[r_, c_, rank[m]] = vals[m]
    return e0, rel_j, rel_r


def top_relatives(cl, ek, C, n_el, T, block_nnz=2e7, threads=1):
    """For every copy and every other sequence: its T closest relatives there (copy, relatedness), relatedness =
    shared clades weighted by specificity (sum of mult/(size-1)); -1 = none. Blocks of copies run in a fork pool
    (scipy only, no numba parallel code inside, so the pool is fork-safe)."""
    A = sparse.csr_matrix((np.ones(cl.members.size, np.float32), cl.members, cl.indptr), shape=(cl.n, n_el))
    d = (cl.mult / np.maximum(cl.size - 1, 1)).astype(np.float32)
    DA = sparse.csr_matrix((np.repeat(d, cl.size), cl.members, cl.indptr), shape=(cl.n, n_el))
    cost = np.bincount(cl.members, weights=np.repeat(cl.size, cl.size).astype(float), minlength=n_el)
    cum = np.concatenate(([0], np.cumsum(cost)))
    blocks, e0 = [], 0
    while e0 < n_el:
        e1 = int(np.searchsorted(cum, cum[e0] + block_nnz, side="right"))
        e1 = min(max(e1, e0 + 1), n_el)
        blocks.append((e0, e1))
        e0 = e1
    _G.update(AT=A.T.tocsr(), DA=DA, ek=ek, C=C, T=T)
    rel_j = np.full((n_el, C, T), -1, np.int32)
    rel_r = np.zeros((n_el, C, T), np.float32)
    for e0, bj, br in pool_map(_relatives_block, blocks, threads):
        rel_j[e0:e0 + len(bj)] = bj
        rel_r[e0:e0 + len(br)] = br
    for x in ("AT", "DA", "ek"):
        _G.pop(x, None)
    return rel_j, rel_r


@njit(parallel=True, cache=True)
def _knn_kernel(rel_j, rel_r, ek, mask, C, k, nchunk):
    """Per chunk of copies: each copy in mask (phased) adds 1 to (its sequence, each of the k sequences holding its
    closest relatives within mask); ties go to the lower sequence index."""
    n, _, T = rel_j.shape
    out = np.zeros((nchunk, C, C))
    for ch in prange(nchunk):
        best = np.zeros(C)
        for e in range(n * ch // nchunk, n * (ch + 1) // nchunk):
            if not mask[e] or ek[e] < 0:
                continue
            for c in range(C):
                b = 0.0
                for t in range(T):
                    f = rel_j[e, c, t]
                    if f >= 0 and mask[f] and rel_r[e, c, t] > b:
                        b = rel_r[e, c, t]
                best[c] = b
            for _ in range(k):
                m, arg = 0.0, -1
                for c in range(C):
                    if best[c] > m:
                        m, arg = best[c], c
                if arg < 0:
                    break
                out[ch, ek[e], arg] += 1.0
                best[arg] = 0.0
    return out


def knn_votes(rel_j, rel_r, ek, C, mask):
    """C x C symmetric nearest-relative votes: every copy in mask adds 1 to its sequence x each of the K_NEIGH
    other sequences holding its closest relatives within mask (each copy votes once)."""
    k = min(K_NEIGH, C - 1)
    nchunk = n_chunks(C, rel_j.shape[0])
    W = _knn_kernel(rel_j, rel_r, np.asarray(ek, np.int64), np.asarray(mask, np.bool_), C, k, nchunk).sum(0)
    W = W + W.T
    np.fill_diagonal(W, 0)
    return W


def combine(L, Kv):
    if L.sum() > 0 and Kv.sum() > 0:
        return L + Kv * (L.sum() / Kv.sum())
    return L if L.sum() > 0 else Kv


class LinkModel:
    """Clades restricted to phased copies, and each copy's closest relatives (clades <= NMAX_KNN)."""

    def __init__(self, cl, ek, C, threads=1):
        self.ek, self.C = ek, C
        self.n = ek.size
        phased = ek >= 0
        self.cl = cl.restrict(phased, ek.clip(min=0))
        small = Clades(*_subset(self.cl, self.cl.size <= NMAX_KNN))
        T = int(max(1, min(4, 5e7 // max(self.n * C, 1))))
        self.rel_j, self.rel_r = top_relatives(small, ek, C, self.n, T, threads=threads)

    def links(self, mask):
        return link_matrix(self.cl, self.ek, self.C, mask)

    def matrix(self, mask):
        return combine(self.links(mask), knn_votes(self.rel_j, self.rel_r, self.ek, self.C, mask))


def _subset(cl, keep):
    km = np.repeat(keep, cl.size)
    return cl.members[km], cl.size[keep], cl.mult[keep], cl.nseq[keep], cl.n_occ, cl.n_multi


def pool_map(fn, jobs, threads):
    """Ordered map with a fork process pool (serial if threads <= 1 or fork unavailable). fn must not call numba
    parallel code (the OpenMP threading layer is not fork-safe)."""
    jobs = list(jobs)
    if threads <= 1 or len(jobs) <= 1:
        return [fn(j) for j in jobs]
    import multiprocessing as mp
    try:
        ctx = mp.get_context("fork")
    except ValueError:
        return [fn(j) for j in jobs]
    with ctx.Pool(min(threads, len(jobs))) as pool:
        return pool.map(fn, jobs, chunksize=1)


# ------------------------------------------------------------------------------------------------ partition
from scipy.optimize import linear_sum_assignment  # noqa: E402

EXACT_MAX_SETS = 26  # exhaustive N=2 search up to 2^25 configurations (Gray code: ~1 s)


@njit(cache=True)
def gray_best(M, allsingle):
    """Maximise s^T M s over s in {-1, +1}^V with s_0 = +1, visiting configurations in Gray-code order (one sign
    flips per step; the score and the local fields update in O(V)). With all sets singletons, all-plus (one group)
    is excluded."""
    V = M.shape[0]
    s = np.ones(V)
    h = np.zeros(V)
    q = 0.0
    for i in range(V):
        for w in range(V):
            q += M[i, w]
            if w != i:
                h[i] += M[i, w]
    best = -np.inf
    bs = s.copy()
    if not allsingle:
        best = q
    nneg = 0
    total = 1 << (V - 1)
    for g in range(1, total):
        j = 1
        x = g
        while (x & 1) == 0:
            x >>= 1
            j += 1
        sj = s[j]
        q += -4.0 * sj * h[j]
        s[j] = -sj
        nneg += 1 if sj > 0 else -1
        for w in range(V):
            if w != j:
                h[w] -= 2.0 * M[w, j] * sj
        if (not allsingle or nneg > 0) and q > best:
            best = q
            bs[:] = s
    return bs


def _exhaustive2(B, sets):
    C, V = B.shape[0], len(sets)
    U = np.zeros((C, V))
    for v, st in enumerate(sets):
        U[st[0], v] = 1.0
        if len(st) == 2:
            U[st[1], v] = -1.0
    s = gray_best(U.T @ B @ U, all(len(st) == 1 for st in sets))
    return ((U @ s) < 0).astype(int)


def make_sets(names, config_sets, N):
    """Cannot-link sets (lists of indices); every sequence in exactly one set."""
    ix = {c: i for i, c in enumerate(names)}
    sets, seen = [], set()
    for st in config_sets or []:
        mem = [ix[c] for c in st if c in ix and ix[c] not in seen]
        if len(mem) > N:
            die(f"config set {' '.join(st)} has more sequences than subgenomes (N={N})")
        if mem:
            sets.append(mem)
            seen.update(mem)
    sets += [[i] for i in range(len(names)) if i not in seen]
    return sets


def q_score(B, z):
    return 0.5 * float(B[z[:, None] == z[None, :]].sum())


def _local_search(B, z, sets, N, rng, max_rounds=200):
    z = z.copy()
    C = B.shape[0]
    for _ in range(max_rounds):
        changed = False
        for si in rng.permutation(len(sets)):
            mem = sets[si]
            oh = np.zeros((C, N))
            oh[np.arange(C), z] = 1
            oh[mem] = 0
            A = B[mem] @ oh
            if len(mem) == 1:
                new = [int(np.argmax(A[0]))]
            else:
                r, c = linear_sum_assignment(-A)
                new = list(c[np.argsort(r)])
            for c_, g in zip(mem, new):
                if z[c_] != g:
                    z[c_] = g
                    changed = True
        if not changed:
            break
    return z


def _spectral_init(B, N, rng):
    _, V = np.linalg.eigh(B)
    X = V[:, -max(N - 1, 1):]
    if N == 2:
        return (X[:, -1] < 0).astype(int)
    cen = X[rng.choice(len(X), N, replace=False)]
    for _ in range(50):
        z = np.argmin(((X[:, None, :] - cen[None]) ** 2).sum(-1), 1)
        cen = np.array([X[z == g].mean(0) if (z == g).any() else X[rng.integers(len(X))] for g in range(N)])
    return z


def best_partition(B, N, sets, rng, restarts=30, init=None):
    """Maximum within-group excess (degree-corrected modularity) partition into N groups (exact for N=2)."""
    C = B.shape[0]
    if N == 2 and all(len(s) <= 2 for s in sets) and len(sets) <= EXACT_MAX_SETS:
        return _exhaustive2(B, sets)
    inits = [_spectral_init(B, N, rng)]
    if init is not None:
        inits.insert(0, np.asarray(init))
    inits += [rng.integers(0, N, C) for _ in range(restarts)]
    best, bz = -np.inf, None
    for z0 in inits:
        z = _local_search(B, z0, sets, N, rng)
        if len(set(z.tolist())) < N:
            continue
        q = q_score(B, z)
        if q > best + 1e-12:
            best, bz = q, z
    return bz if bz is not None else _local_search(B, inits[0], sets, N, rng)


def canon_labels(z, mass):
    """Relabel groups so that group 0 holds the most LTR-RTs (ties: lowest index)."""
    groups = sorted(set(z.tolist()), key=lambda g: (-mass[z == g].sum(), np.flatnonzero(z == g)[0]))
    m = {g: i for i, g in enumerate(groups)}
    return np.array([m[g] for g in z])


def align(ref, z, N):
    cont = np.zeros((N, N))
    np.add.at(cont, (z, ref), 1)
    r, c = linear_sum_assignment(-cont)
    m = dict(zip(r, c))
    return np.array([m[g] for g in z])


def ari(a, b):
    """Adjusted Rand index (1 = identical, ~0 = chance)."""
    from scipy.special import comb
    a, b = np.asarray(a), np.asarray(b)
    _, ai = np.unique(a, return_inverse=True)
    _, bi = np.unique(b, return_inverse=True)
    ct = np.zeros((ai.max() + 1, bi.max() + 1))
    np.add.at(ct, (ai, bi), 1)
    s = comb(ct, 2).sum()
    sa, sb = comb(ct.sum(1), 2).sum(), comb(ct.sum(0), 2).sum()
    exp = sa * sb / comb(a.size, 2)
    mx = (sa + sb) / 2
    return 1.0 if mx == exp else float((s - exp) / (mx - exp))


def bipartition(W, idx, rng, sets=None):
    """Best 2-way split of sequences idx from links among them only (null recomputed in the group); sets: cannot-
    link pairs given as local indices into idx (each pair split across the two parts)."""
    Wg = np.asarray(W[np.ix_(idx, idx)], np.float64)
    if Wg.sum() <= 0 or len(idx) < 2:
        return np.zeros(len(idx), int)
    return best_partition(Wg - expected(Wg), 2, local_sets(sets, len(idx)), rng, restarts=10)


def local_sets(pairs, n):
    """Cannot-link pairs (local indices) completed with singletons into sets covering 0..n-1."""
    pairs = [list(p) for p in (pairs or [])]
    seen = {x for p in pairs for x in p}
    return pairs + [[i] for i in range(n) if i not in seen]


def split_contrast(W, idx, z):
    """Mean log(O/E) within the two parts minus between them (null recomputed inside the group)."""
    L = log_oe(W[np.ix_(idx, idx)])
    same = z[:, None] == z[None, :]
    off = ~np.eye(len(idx), dtype=bool)
    if not (same & off).any() or not (~same).any():
        return 0.0
    return float(L[same & off].mean() - L[~same].mean())


def rep_stats(pairs, g, rng, sets=None):
    """Replicability of the split of group g: mean adjusted Rand between the best splits found independently in the
    two halves of each pair, its standard error, and the halves' splits (for per-sequence support)."""
    hz = [(bipartition(Wa, g, rng, sets), bipartition(Wb, g, rng, sets)) for Wa, Wb in pairs]
    v = np.array([ari(x, y) for x, y in hz])
    return float(v.mean()), (float(v.std(ddof=1) / math.sqrt(v.size)) if v.size > 1 else 0.0), hz


def side_support(z, hz):
    """Per member of a split group: fraction of half-library splits (each side-aligned to z) that put it on its side."""
    hits = []
    for pair in hz:
        for h in pair:
            h = np.asarray(h)
            hits.append((h if (h == z).mean() >= 0.5 else 1 - h) == z)
    return np.mean(hits, 0) if hits else np.full(len(z), np.nan)


# ------------------------------------------------------------------------------------------------ number of subgenomes
class Sources:
    """Link matrices of all phased copies ('all') and of each LTR-K2P tercile (AGE_STRATA), with complementary
    half-library pairs drawn from the same random halves for every source; more(s, n) adds n pairs to source s."""

    def __init__(self, lm, k2p, rng):
        self.lm, self.rng = lm, rng
        ph = lm.ek >= 0
        ok = ph & np.isfinite(k2p)
        self.masks = {"all": ph}
        self.edges = None
        if ok.sum() >= len(AGE_STRATA):
            q = np.quantile(k2p[ok], np.linspace(0, 1, len(AGE_STRATA) + 1)[1:-1])
            self.edges = [float(x) for x in q]
            lo = -np.inf
            for name, hi in zip(AGE_STRATA, list(q) + [np.inf]):
                self.masks[name] = ok & (k2p > lo) & (k2p <= hi)
                lo = hi
        self.names = list(self.masks)
        self.W = {s: lm.matrix(m) for s, m in self.masks.items()}
        self.pairs = {s: [] for s in self.names}
        for seed in rng.integers(0, 2 ** 62, HALF_PAIRS):
            u = np.random.default_rng(int(seed)).random(lm.n) < 0.5
            for s, m in self.masks.items():
                self.pairs[s].append(self._pair(u, m))

    def _pair(self, u, m):
        return self.lm.matrix(u & m).astype(np.float32), self.lm.matrix(~u & m).astype(np.float32)

    def more(self, s, n):
        new = [self._pair(np.random.default_rng(int(sd)).random(self.lm.n) < 0.5, self.masks[s])
               for sd in self.rng.integers(0, 2 ** 62, n)]
        self.pairs[s].extend(new)
        return new


def split_candidate(src, s, g, rng, sets=None, min_size=2):
    """Best bipartition of group g on source s's links, its strength, and its replicability in s's half-library pairs
    (more pairs drawn while the mean straddles MIN_REP within 2 s.e., up to MAX_PAIRS); None without links or when a
    side would hold fewer than min_size sequences."""
    W = src.W[s]
    if W[np.ix_(g, g)].sum() <= 0:
        return None
    z = bipartition(W, g, rng, sets)
    if min(np.bincount(z, minlength=2)) < min_size:
        return None
    rep, se, hz = rep_stats(src.pairs[s], g, rng, sets)
    while len(src.pairs[s]) < MAX_PAIRS and abs(rep - MIN_REP) < 2 * se:
        src.more(s, HALF_PAIRS)
        rep, se, hz = rep_stats(src.pairs[s], g, rng, sets)
    return dict(src=s, z=z, con=split_contrast(W, g, z), rep=rep, se=se, n_pairs=len(src.pairs[s]),
                borderline=abs(rep - MIN_REP) < 2 * se, hz=hz)


def parent_info(src, g, z, cands):
    """What children are measured against, on each source: the strongest split of the parent group there (the
    accepted split or that source's own best split, whichever is stronger)."""
    con = {}
    for s in src.names:
        W = src.W[s]
        own = split_contrast(W, g, z) if W[np.ix_(g, g)].sum() > 0 else 0.0
        con[s] = max(own, cands[s]["con"]) if s in cands else own
    return dict(con=con)


def accept_split(cands, parent):
    """The v7 rule. A split passes on some copies when half-libraries of those copies reproduce it (adjusted Rand >=
    MIN_REP), it has positive strength and, below the top, it is >= MIN_REL as strong as the strongest split of its
    parent group on the same copies. Pooled copies decide first (the v5/v6 rule). Otherwise one LTR-K2P tercile can
    carry it, provided the split is age-bounded: some tercile holds no reproducible split of the group (confinement to
    a progenitor has a window; structure present at every age is the chromosomes' own). The top split is decided on
    pooled copies only (no parent to measure strength against). -> (candidate, 'pooled' | 'age') or (None, reason)."""
    def strong(c):
        if c["rep"] < MIN_REP or c["con"] <= 0:
            return False
        return parent is None or (parent["con"][c["src"]] > 0 and c["con"] >= MIN_REL * parent["con"][c["src"]])
    p = cands.get("all")
    if p is not None and strong(p):
        return p, "pooled"
    if parent is None:
        return None, "not reproduced" if p is None or p["rep"] < MIN_REP else "no excess"
    strata = [s for s in AGE_STRATA if s in parent["con"]]
    seen = [s for s in strata if cands.get(s) is not None and cands[s]["rep"] >= MIN_REP]
    ok = [cands[s] for s in seen if strong(cands[s])]
    if ok and len(seen) < len(strata):
        return max(ok, key=lambda c: (round(c["rep"], 2), c["con"] / parent["con"][c["src"]])), "age"
    if ok:
        return None, "present at every age"
    if (p is not None and p["rep"] >= MIN_REP) or seen:
        return None, "weaker than half its parent"
    return None, "not reproduced"


def find_subgenomes(src, rng, n_force=None, min_size=2, max_n=MAX_N, sets_of=None):
    """Divisive search on all copies and on each LTR-K2P tercile. Inferred N: every group is tested once, splits are
    kept by accept_split, N=1 is a possible answer. Forced N (n_force): the group whose best candidate ranks first
    (accepted, then replicability, then strength relative to its parent) is split until there are n_force groups.
    sets_of(g) -> cannot-link pairs (local indices) for the bipartitions of g. Returns labels, every tested split,
    and per-sequence support (the minimum, over the splits that placed it, of its side-support)."""
    C = src.W["all"].shape[0]
    support = np.full(C, np.nan)
    tests = []

    def evaluate(g, parent, path):
        sets = sets_of(g) if sets_of else None
        cands = {}
        if len(g) >= 2 * min_size:
            for s in src.names:
                c = split_candidate(src, s, g, rng, sets, min_size)
                if c is not None:
                    cands[s] = c
        ch, how = accept_split(cands, parent) if cands else (None, "too small" if len(g) < 2 * min_size else "no split")
        t = dict(members=g, path=path, cands=cands, parent=parent, chosen=ch, how=how, accepted=ch is not None,
                 forced=False)
        tests.append(t)
        return t

    def rank(t):
        c = t["chosen"] or max(t["cands"].values(), key=lambda c: (round(c["rep"], 2), c["con"]))
        par = t["parent"]["con"][c["src"]] if t["parent"] else 1.0
        return (t["accepted"], round(c["rep"], 2), c["con"] / max(par, 1e-12)), c

    leaves = [evaluate(np.arange(C), None, "root")]
    final = []
    while leaves:
        if n_force is None:
            t = leaves.pop(0)
            if not t["accepted"] or len(final) + len(leaves) + 2 > max_n:
                final.append(t["members"])
                continue
            c = t["chosen"]
        else:
            open_ = [x for x in leaves if x["cands"]]
            if len(final) + len(leaves) >= n_force or not open_:
                final += [x["members"] for x in leaves]
                break
            t = max(open_, key=lambda x: rank(x)[0])
            leaves = [x for x in leaves if x is not t]              # identity: tests hold arrays
            if not t["accepted"]:
                t["chosen"], t["forced"] = rank(t)[1], True
            c = t["chosen"]
        g, z = t["members"], c["z"]
        sup = side_support(z, c["hz"])
        support[g] = np.fmin(support[g], sup)                      # NaN (not yet placed) gives way
        pi = parent_info(src, g, z, t["cands"])
        leaves += [evaluate(g[z == 0], pi, t["path"] + ".a"), evaluate(g[z == 1], pi, t["path"] + ".b")]
    lab = np.zeros(C, int)
    for k, g in enumerate(sorted(final, key=lambda x: x.min())):
        lab[g] = k
    return lab, tests, support


def support_from_halves(halves, z, N, sets, rng):
    """Fraction of half-library partitions (solved from scratch with the full-data estimator) agreeing with z."""
    hits = []
    for W in halves:
        if W.sum() <= 0:
            continue
        zr = best_partition(W - expected(W), N, sets, rng)
        hits.append(align(z, zr, N) == z)
    return np.mean(hits, 0) if hits else np.full(len(z), np.nan)


def partner_pairs(W, halves, z):
    """Closest partner of each sequence within its own group (most shared lineages, observed/expected), kept when
    mutual; support = fraction of halves where it still is; excess = the pair's O/E over the sequence's median
    within its group (near-identical homologues or progenitors show large excess)."""
    C = W.shape[0]
    same = z[:, None] == z[None, :]
    np.fill_diagonal(same, False)

    def best_of(M):
        E = expected(M)
        oe = np.where(same & (E > 0), M / np.where(E > 0, E, 1.0), np.nan)
        ok = np.isfinite(oe).any(1)
        b = np.where(ok, np.nanargmax(np.where(np.isfinite(oe), oe, -np.inf), 1), -1)
        return b, oe
    b, oe = best_of(W)
    partner = np.array([j if j >= 0 and b[j] == i else -1 for i, j in enumerate(b)])
    support, excess = np.full(C, np.nan), np.full(C, np.nan)
    hb = [best_of(h)[0] for h in halves]
    for i, j in enumerate(partner):
        if j < 0:
            continue
        support[i] = float(np.mean([bb[i] == j and bb[j] == i for bb in hb])) if hb else np.nan
        rest = oe[i][np.isfinite(oe[i])]
        excess[i] = float(oe[i, j] / np.median(rest)) if rest.size > 1 and np.median(rest) > 0 else np.nan
    return partner, support, excess


def chrom_tests(W, z, N):
    """Per sequence: mean O/E to its own group and to the closest other group; one-sided Welch test with the other
    sequences as replicates."""
    from scipy.stats import ttest_ind
    C = len(z)
    E = expected(W)
    L = np.log((W + 1) / (E + 1))
    OE = W / np.where(E > 0, E, np.nan)
    own_oe, oth_oe, pval = (np.full(C, np.nan) for _ in range(3))
    for c in range(C):
        others = np.arange(C) != c
        x = L[c, others & (z == z[c])]
        best, bg = -np.inf, None
        for g in range(N):
            if g != z[c] and (z == g).any() and L[c, z == g].mean() > best:
                best, bg = L[c, z == g].mean(), g
        if bg is None:
            continue
        y = L[c, others & (z == bg)]
        own_oe[c] = np.nanmean(OE[c, others & (z == z[c])]) if x.size else np.nan
        oth_oe[c] = np.nanmean(OE[c, others & (z == bg)])
        if x.size >= 2 and y.size >= 2:
            t, p = ttest_ind(x, y, equal_var=False)
            pval[c] = p / 2 if t > 0 else 1 - p / 2
    return own_oe, oth_oe, pval


# ------------------------------------------------------------------------------------------------ copy pairs
import edlib  # noqa: E402

ALPH = np.frombuffer(b"ACGTN", np.uint8)
RCODE = np.array([3, 2, 1, 0, 4], np.uint8)
_CIGAR = re.compile(rb"(\d+)([=XID])")
_OPCODE = {ord("="): 0, ord("X"): 0, ord("I"): 1, ord("D"): 2}
_G = {}             # state shared with fork-pool workers (set before forking)


@njit(cache=True)
def _band(n, m, lens, ops, B):
    """Per query row i (0..n): target window [lo, hi] spanning every diagonal (j - i) from the lowest to the highest
    among 0, m - n and those the edlib path visits, +- B. A long indel's single-gap path runs between diagonals 0 and
    m - n, which edlib's unit-cost path (many short gaps along a drifting diagonal) can leave by far more than B."""
    i = 0
    j = 0
    dmin = min(0, m - n)
    dmax = max(0, m - n)
    for k in range(lens.size):
        for _ in range(lens[k]):
            if ops[k] == 0:
                i += 1
                j += 1
            elif ops[k] == 1:
                i += 1
            else:
                j += 1
            dmin = min(dmin, j - i)
            dmax = max(dmax, j - i)
    lo = np.empty(n + 1, np.int64)
    hi = np.empty(n + 1, np.int64)
    for r in range(n + 1):
        lo[r] = max(0, r + dmin - B)
        hi[r] = min(m, r + dmax + B)
    lo[0] = 0
    hi[n] = m
    return lo, hi


@njit(cache=True)
def _gotoh_band(q, t, lo, hi, sm, sx, go, ge):
    """Global affine-gap alignment restricted to the band (lo, hi nondecreasing) -> (sites with both bases ACGT,
    transitions, transversions, touched: the best path runs along an inner band edge, so a wider band may score
    better). Ties prefer the diagonal (substitution), then a gap in the query, then in the target. Each row runs in
    three vectorizable passes: gap-in-target and diagonal scores from the previous row, the gap-in-query chain as a
    running maximum (E[j] = max over k < j of Hd[k] - go - ge (j - k), Hd = best without E), then the row's scores.
    int32 scores (no overflow below ~10^8 bp); two rolling rows; only the traceback (1 byte per band cell) is kept."""
    n = q.size
    m = t.size
    NEG = np.int32(-(1 << 29))
    goe = np.int32(go + ge)
    ge = np.int32(ge)
    prof = np.zeros((5, m + 1), np.int32)             # prof[x, j]: score of query base x against t[j - 1]
    for x in range(4):
        for j in range(1, m + 1):
            y = t[j - 1]
            prof[x, j] = 0 if y > 3 else (sm if x == y else sx)
    base = np.zeros(n + 2, np.int64)
    W = 0
    for i in range(n + 1):
        base[i + 1] = base[i] + hi[i] - lo[i] + 1
        W = max(W, hi[i] - lo[i] + 1)
    P = np.empty(base[n + 1], np.uint8)              # bits 0-1: H source (0 diag, 1 E, 2 F); bit 2: E extends; 3: F extends
    Hp = np.full(W + 2, NEG, np.int32)               # rows indexed 1 + (j - lo[i]); padded with NEG on both sides
    Fp = np.full(W + 2, NEG, np.int32)
    Hc = np.full(W + 2, NEG, np.int32)
    Fc = np.full(W + 2, NEG, np.int32)
    D = np.empty(W, np.int32)
    Hd = np.empty(W, np.int32)
    E = np.empty(W, np.int32)
    w0 = hi[0] - lo[0] + 1
    for c in range(w0):                               # row 0: leading gap in the query
        j = lo[0] + c
        if j == 0:
            Hc[1 + c] = 0
            P[base[0] + c] = 0
        else:
            Hc[1 + c] = -(go + ge * j)
            P[base[0] + c] = 1 | (4 if j > 1 else 0)
    wold_c = w0                                       # widths last written into the two row buffers
    wold_p = 0
    for i in range(1, n + 1):
        Hp, Hc = Hc, Hp
        Fp, Fc = Fc, Fp
        wold_c, wold_p = wold_p, wold_c
        l0 = lo[i]
        w = hi[i] - l0 + 1
        sh = l0 - lo[i - 1]
        pr = prof[q[i - 1], l0:l0 + w]
        HpA = Hp[1 + sh:1 + sh + w]                   # previous row at column l0 + c
        HpB = Hp[sh:sh + w]                           # previous row at column l0 + c - 1
        FpA = Fp[1 + sh:1 + sh + w]
        FcV = Fc[1:1 + w]
        HcV = Hc[1:1 + w]
        PV = P[base[i]:base[i] + w]
        Dv = D[:w]
        Hdv = Hd[:w]
        Ev = E[:w]
        for c in range(w):
            fo = HpA[c] - goe
            fe = FpA[c] - ge
            f = max(fe, fo)
            FcV[c] = f
            d = HpB[c] + pr[c]
            Dv[c] = d
            Hdv[c] = max(d, f)
            PV[c] = 8 if fe >= fo else 0
        if l0 == 0:                                   # column 0: no diagonal
            Dv[0] = NEG
            Hdv[0] = FcV[0]
        run = NEG
        Ev[0] = NEG
        for c in range(1, w):
            run = max(run - ge, Hdv[c - 1] - goe)
            Ev[c] = run
        for c in range(w):
            e = Ev[c]
            d = Dv[c]
            f = FcV[c]
            ef = max(e, f)
            HcV[c] = max(d, ef)
            PV[c] = PV[c] | ((d < ef) * (2 - (e >= f)))
        for c in range(1, w):
            PV[c] = PV[c] | (4 * (Ev[c - 1] - ge >= Hdv[c - 1] - goe))
        for c in range(w + 1, wold_c + 1):            # clear what an older, wider row left beyond this one
            Hc[c] = NEG
            Fc[c] = NEG
        wold_c = w
    i = n
    j = m
    state = 0
    sites = 0
    ts = 0
    tv = 0
    touched = False
    while i > 0 or j > 0:
        if (j == lo[i] and lo[i] > 0) or (j == hi[i] and hi[i] < m):
            touched = True
        pb = P[base[i] + j - lo[i]]
        if state == 0:
            arg = pb & 3
            if arg == 0:
                x = q[i - 1]
                y = t[j - 1]
                if x < 4 and y < 4:
                    sites += 1
                    if x != y:
                        if (x ^ y) == 2:
                            ts += 1
                        else:
                            tv += 1
                i -= 1
                j -= 1
            else:
                state = arg
        elif state == 1:
            j -= 1
            if (pb & 4) == 0:
                state = 0
        else:
            i -= 1
            if (pb & 8) == 0:
                state = 0
    return sites, ts, tv, touched


def _rc(codes):
    return RCODE[codes][::-1]


def align_counts(q, t):
    """Global alignment of two code arrays -> (aligned sites, transitions, transversions): the edlib path fixes a
    band (widened while the best path rides its edge), an affine-gap realignment inside it counts the substitutions
    (unit costs undercount them)."""
    if q.size == 0 or t.size == 0:
        return 0, 0, 0
    r = edlib.align(ALPH[q].tobytes(), ALPH[t].tobytes(), mode="NW", task="path")
    m = _CIGAR.findall(r["cigar"].encode())
    lens = np.array([int(x) for x, _ in m], np.int64)
    ops = np.array([_OPCODE[o[0]] for _, o in m], np.int64)
    B = BAND
    while True:                                 # widen 4-fold while the best path rides the band edge
        lo, hi = _band(q.size, t.size, lens, ops, B)
        sites, ts, tv, touched = _gotoh_band(q, t, lo, hi, *SCORES)
        if not touched or (lo[-1] == 0 and hi[0] == t.size):
            return sites, ts, tv
        B *= 4


def _record(lib, i):
    return lib.codes[lib.offs[i]:lib.offs[i + 1]]


def _ltrs(lib, a, b):
    """LTR code arrays (a5, a3, b5, b3) of copies a and b in one orientation; None when either copy has no LTR bounds.
    Stranded records share element orientation; when a strand is unknown, b is reverse-complemented when that brings
    its LTRs closer to a's (smaller edit distance)."""
    if lib.l5[a, 0] < 0 or lib.l5[b, 0] < 0:
        return None
    ra, rb = _record(lib, a), _record(lib, b)
    a5, a3 = ra[lib.l5[a, 0]:lib.l5[a, 1]], ra[lib.l3[a, 0]:lib.l3[a, 1]]
    b5, b3 = rb[lib.l5[b, 0]:lib.l5[b, 1]], rb[lib.l3[b, 0]:lib.l3[b, 1]]
    if lib.strand_code[a] == 0 or lib.strand_code[b] == 0:
        if edlib.align(ALPH[a5].tobytes(), ALPH[_rc(b3)].tobytes(), mode="NW")["editDistance"] < \
                edlib.align(ALPH[a5].tobytes(), ALPH[b5].tobytes(), mode="NW")["editDistance"]:
            b5, b3 = _rc(b3), _rc(b5)
    return a5, a3, b5, b3


def pair_counts(lib, a, b):
    """Between-copy counts of copies a and b: 5'-5' plus 3'-3' LTR alignments (whole records when either copy has no
    LTR bounds). Stranded records share element orientation; when a strand is unknown, the orientation with the
    smaller edit distance is used."""
    lt = _ltrs(lib, a, b)
    if lt is None:
        ra, rb = _record(lib, a), _record(lib, b)
        if lib.strand_code[a] == 0 or lib.strand_code[b] == 0:
            if edlib.align(ALPH[ra].tobytes(), ALPH[_rc(rb)].tobytes(), mode="NW")["editDistance"] < \
                    edlib.align(ALPH[ra].tobytes(), ALPH[rb].tobytes(), mode="NW")["editDistance"]:
                rb = _rc(rb)
        return align_counts(ra, rb)
    a5, a3, b5, b3 = lt
    s5, s3 = align_counts(a5, b5), align_counts(a3, b3)
    return s5[0] + s3[0], s5[1] + s3[1], s5[2] + s3[2]


def _pair_job(chunk):
    lib = _G["lib"]
    out = np.zeros((len(chunk), 3), np.int64)
    for n, (a, b) in enumerate(chunk):
        out[n] = pair_counts(lib, int(a), int(b))
    return out


def align_pairs(lib, pairs, threads):
    """(sites, transitions, transversions) for each copy pair (rows of pairs). A fork pool over chunks: edlib and
    the alignment walker are single-threaded numba-free/non-parallel code, safe after fork."""
    pairs = np.asarray(pairs, np.int64).reshape(-1, 2)
    if len(pairs) == 0:
        return np.zeros((0, 3), np.int64)
    _G["lib"] = lib
    chunks = np.array_split(pairs, max(1, min(len(pairs) // 500, 16 * threads)))
    return np.concatenate(pool_map(_pair_job, chunks, threads))


def whole_counts(lib, a, b):
    """Whole-record counts of copies a and b. Stranded records share element orientation; when a strand is unknown,
    the orientation with the smaller edit distance is used."""
    ra, rb = _record(lib, a), _record(lib, b)
    if lib.strand_code[a] == 0 or lib.strand_code[b] == 0:
        if edlib.align(ALPH[ra].tobytes(), ALPH[_rc(rb)].tobytes(), mode="NW")["editDistance"] < \
                edlib.align(ALPH[ra].tobytes(), ALPH[rb].tobytes(), mode="NW")["editDistance"]:
            rb = _rc(rb)
    return align_counts(ra, rb)


def _whole_job(chunk):
    lib = _G["lib"]
    out = np.zeros((len(chunk), 3), np.int64)
    for n, (a, b) in enumerate(chunk):
        out[n] = whole_counts(lib, int(a), int(b))
    return out


def align_whole(lib, pairs, threads):
    """(sites, transitions, transversions) of whole records for each copy pair; fork pool (edlib and serial njit
    only: no numba parallel code in workers)."""
    pairs = np.asarray(pairs, np.int64).reshape(-1, 2)
    if len(pairs) == 0:
        return np.zeros((0, 3), np.int64)
    _G["lib"] = lib
    align_counts(lib.codes[:min(lib.codes.size, 50)], lib.codes[:min(lib.codes.size, 50)])     # compile in the parent
    chunks = np.array_split(pairs, max(1, min(len(pairs) // 200, 16 * threads)))
    return np.concatenate(pool_map(_whole_job, chunks, threads))


@njit(cache=True)
def _mix_k(x):
    x ^= x >> np.uint64(33)
    x *= np.uint64(0xff51afd7ed558ccd)
    x ^= x >> np.uint64(33)
    x *= np.uint64(0xc4ceb9fe1a85ec53)
    x ^= x >> np.uint64(33)
    return x


@njit(parallel=True, cache=True)
def sketches(codes, offs, k, s):
    """Bottom-s MinHash of each record's canonical k-mers (rows padded with the maximum uint64)."""
    n = offs.size - 1
    pad = np.uint64(0xFFFFFFFFFFFFFFFF)
    out = np.full((n, s), pad, np.uint64)
    mask = np.uint64((1 << (2 * k)) - 1)
    for i in prange(n):
        a0 = offs[i]
        a1 = offs[i + 1]
        m = a1 - a0 - k + 1
        if m <= 0:
            continue
        h = np.empty(m, np.uint64)
        nh = 0
        fw = np.uint64(0)
        rv = np.uint64(0)
        good = 0
        for p in range(a0, a1):
            c = codes[p]
            if c > 3:
                good = 0
                fw = np.uint64(0)
                rv = np.uint64(0)
                continue
            fw = ((fw << np.uint64(2)) | np.uint64(c)) & mask
            rv = (rv >> np.uint64(2)) | (np.uint64(3 - c) << np.uint64(2 * (k - 1)))
            good += 1
            if good >= k:
                h[nh] = _mix_k(min(fw, rv))
                nh += 1
        if nh == 0:
            continue
        h = np.sort(h[:nh])
        j = 0
        last = pad
        for q in range(nh):
            if h[q] != last:
                out[i, j] = h[q]
                last = h[q]
                j += 1
                if j == s:
                    break
    return out


@njit(cache=True)
def _sketch_index(S, use):
    """(hash, copy) of the sketches of copies in use, sorted by hash."""
    n, s = S.shape
    pad = np.uint64(0xFFFFFFFFFFFFFFFF)
    tot = 0
    for i in range(n):
        if use[i]:
            for j in range(s):
                if S[i, j] != pad:
                    tot += 1
    H = np.empty(tot, np.uint64)
    C = np.empty(tot, np.int64)
    k = 0
    for i in range(n):
        if use[i]:
            for j in range(s):
                if S[i, j] != pad:
                    H[k] = S[i, j]
                    C[k] = i
                    k += 1
    o = np.argsort(H)
    return H[o], C[o]


@njit(parallel=True, cache=True)
def _top_by_index(S, seq, focal, Hs, Cs, m, cap, nchunk):
    n, s = S.shape
    nf = focal.size
    J = np.full((nf, m), -1, np.int64)
    V = np.zeros((nf, m), np.int64)
    pad = np.uint64(0xFFFFFFFFFFFFFFFF)
    for ch in prange(nchunk):
        cnt = np.zeros(n, np.int32)
        touched = np.empty(n, np.int64)
        for f in range(nf * ch // nchunk, nf * (ch + 1) // nchunk):
            i = focal[f]
            nt = 0
            for j in range(s):
                h = S[i, j]
                if h == pad:
                    break
                lo = np.searchsorted(Hs, h)
                hi = np.searchsorted(Hs, h, side="right")
                if hi - lo > cap:
                    continue
                for q in range(lo, hi):
                    c = Cs[q]
                    if seq[c] == seq[i]:
                        continue
                    if cnt[c] == 0:
                        touched[nt] = c
                        nt += 1
                    cnt[c] += 1
            for t in range(nt):
                c = touched[t]
                v = cnt[c]
                cnt[c] = 0
                if v < 2 or v < V[f, m - 1] or (v == V[f, m - 1] and J[f, m - 1] >= 0 and c > J[f, m - 1]):
                    continue
                q = m - 1
                while q > 0 and (V[f, q - 1] < v or (V[f, q - 1] == v and J[f, q - 1] > c)):
                    V[f, q] = V[f, q - 1]
                    J[f, q] = J[f, q - 1]
                    q -= 1
                V[f, q] = v
                J[f, q] = c
    return J, V


def top_relatives_sketch(S, seq, focal, m, use=None):
    """For each focal copy, the m copies on other sequences (seq >= 0, in use) sharing the most sketch hashes (>= 2;
    ties to the lower copy index, so results do not depend on threads); -1 = none."""
    use = (seq >= 0) if use is None else use
    Hs, Cs = _sketch_index(S, use)
    return _top_by_index(S, seq, np.asarray(focal, np.int64), Hs, Cs, m, BUCKET_CAP, n_chunks(1, len(focal)))


def k2p_counts(sites, ts, tv):
    """K2P distance and its standard error from aligned sites, transitions and transversions (NaN if saturated)."""
    sites = np.asarray(sites, float)
    n = np.maximum(sites, 1.0)
    P, Q = np.asarray(ts) / n, np.asarray(tv) / n
    w1, w2 = 1 - 2 * P - Q, 1 - 2 * Q
    with np.errstate(all="ignore"):
        d = -0.5 * np.log(w1) - 0.25 * np.log(w2)
        c1, c2 = 1 / w1, 0.5 * (1 / w1 + 1 / w2)
        var = (c1 * c1 * P + c2 * c2 * Q - (c1 * P + c2 * Q) ** 2) / n
    bad = (w1 <= 0) | (w2 <= 0) | (sites == 0)
    return np.where(bad, np.nan, d), np.where(bad, np.nan, np.sqrt(np.maximum(var, 0.0)))


# ------------------------------------------------------------------------------------------------ chromosome pairing

def mutual_best_pairs(lm):
    """Copy pairs that are each other's closest relative on the other's sequence (nearest-relative table), a < b."""
    best = lm.rel_j[:, :, 0]
    e = np.flatnonzero(lm.ek >= 0)
    E = np.repeat(e, lm.C)
    F = best[e].ravel()
    keep = F > E
    E, F = E[keep], F[keep]
    mutual = best[F, lm.ek[E]] == E
    return np.stack([E[mutual], F[mutual]], 1)


@njit(cache=True)
def _chain_dp(x, y, W, sgn):
    """Longest-chain DP over points sorted by x: y increasing (sgn 1) or decreasing (sgn -1), steps <= W in both."""
    n = x.size
    L = np.ones(n, np.int64)
    P = np.full(n, -1, np.int64)
    for i in range(n):
        j = i - 1
        while j >= 0 and x[i] - x[j] <= W:
            dy = (y[i] - y[j]) * sgn
            if x[i] > x[j] and 0 < dy <= W and L[j] + 1 > L[i]:
                L[i] = L[j] + 1
                P[i] = j
            j -= 1
    return L, P


@njit(cache=True)
def _chain_len(x, y, W):
    m = 0
    for sgn in (1, -1):
        L, _ = _chain_dp(x, y, W, sgn)
        if L.size and L.max() > m:
            m = L.max()
    return m


def longest_chain(x, y, W=CHAIN_W):
    """Indices (into x, y; x sorted ascending) of the longest collinear chain, in either orientation."""
    best = None
    for sgn in (1, -1):
        L, P = _chain_dp(np.asarray(x, np.int64), np.asarray(y, np.int64), W, sgn)
        if L.size and (best is None or L.max() > best[0].max()):
            best = (L, P)
    if best is None:
        return np.empty(0, np.int64)
    L, P = best
    i, out = int(np.argmax(L)), []
    while i >= 0:
        out.append(i)
        i = int(P[i])
    return np.array(out[::-1], np.int64)


def chain_significant(n, length, z, n_tests):
    """Chain of `length` among n candidates significant at a Bonferroni level of 0.05 over n_tests sequence pairs.
    When n <= CHAIN_W the window does not constrain chains (any monotone subsequence counts) and the permutation
    null is too discrete for a normal approximation: the exact bound P(chain >= l) <= 2 C(n, l) / l! (expected
    monotone l-subsequences, two orientations) decides. Otherwise z against the permutations decides."""
    from scipy.stats import norm
    alpha = 0.05 / max(n_tests, 1)
    if n <= CHAIN_W:
        lg = math.lgamma
        return math.log(2) + lg(n + 1) - lg(length + 1) - lg(n - length + 1) - lg(length + 1) <= math.log(alpha)
    return z >= float(norm.isf(alpha))


def pairing_map(lib, lm, rng, threads):
    """Sequence pairs joined by collinear chains of mutual-best relatives: insertions shared at matching positions.
    Chain length vs NULL_PERMS permutations of partner positions -> z; significance by chain_significant. Anchors
    of significant chains are aligned; a pair is 'paired' (near-identical sequences: homologues of a young
    polyploid, or haplotypes) when (1) the chain spans >= PAIR_SPAN of both sequences' copies (a local block is a
    segmental duplication or an exchange), (2) most anchors are shared insertions (between-copy K2P below both
    copies' own LTR K2P), and (3) the anchors' median K2P is below the median within-copy K2P of the two sequences'
    copies (they diverged after most of their insertions)."""
    C = lm.C
    pairs = mutual_best_pairs(lm)
    if len(pairs) == 0 or C < 2:
        return []
    rank = np.zeros(lib.n, np.int64)
    for s in range(C):
        idx = np.flatnonzero(lm.ek == s)
        rank[idx[np.argsort(lib.start[idx], kind="stable")]] = np.arange(idx.size)
    sa, sb = lm.ek[pairs[:, 0]], lm.ek[pairs[:, 1]]
    A = np.where(sa <= sb, pairs[:, 0], pairs[:, 1])
    Bp = np.where(sa <= sb, pairs[:, 1], pairs[:, 0])
    key = np.minimum(sa, sb) * C + np.maximum(sa, sb)
    o = np.lexsort((rank[A], key))
    key, A, Bp = key[o], A[o], Bp[o]
    starts = np.r_[0, np.flatnonzero(np.diff(key)) + 1]
    ends = np.r_[starts[1:], key.size]
    n_tests = max(C * (C - 1) // 2, 1)
    n_cp = np.bincount(lm.ek[lm.ek >= 0], minlength=C)
    rows, anchors = [], []
    for s0, e0 in zip(starts, ends):
        if e0 - s0 < 4:
            continue
        x, y = rank[A[s0:e0]], rank[Bp[s0:e0]]
        ch = longest_chain(x, y)
        null = np.array([_chain_len(x, rng.permutation(y), CHAIN_W) for _ in range(NULL_PERMS)], float)
        z = (ch.size - null.mean()) / max(null.std(), 0.5)
        i, j = int(key[s0] // C), int(key[s0] % C)
        xa, yb = x[ch], y[ch]
        rows.append(dict(i=i, j=j, n_candidates=int(e0 - s0), chain=int(ch.size), null_mean=float(null.mean()),
                         z=float(z), significant=chain_significant(int(e0 - s0), int(ch.size), z, n_tests),
                         span_a=float((xa.max() - xa.min() + 1) / n_cp[i]), span_b=float((yb.max() - yb.min() + 1) / n_cp[j]),
                         shared=float("nan"), d_median=float("nan"), k_median=float("nan"), paired=False))
        anchors.append(np.stack([A[s0:e0][ch], Bp[s0:e0][ch]], 1))
    sig = [k for k, r in enumerate(rows) if r["significant"]]
    if sig:
        cnt = align_pairs(lib, np.concatenate([anchors[k] for k in sig]), threads)
        d, _ = k2p_counts(cnt[:, 0], cnt[:, 1], cnt[:, 2])
        at = 0
        for k in sig:
            r, m = rows[k], len(anchors[k])
            dk = d[at:at + m]
            at += m
            within = lib.k2p[(lm.ek == r["i"]) | (lm.ek == r["j"])]
            if np.isfinite(dk).any() and np.isfinite(within).any():
                own = np.minimum(lib.k2p[anchors[k][:, 0]], lib.k2p[anchors[k][:, 1]])
                with np.errstate(invalid="ignore"):
                    r["shared"] = float(np.mean(dk < own))                     # NaN compares False
                r["d_median"], r["k_median"] = float(np.nanmedian(dk)), float(np.nanmedian(within))
                r["paired"] = bool(min(r["span_a"], r["span_b"]) >= PAIR_SPAN and r["shared"] > 0.5
                                   and r["d_median"] < r["k_median"])
    return rows


def paired_sets(rows, C):
    """Cannot-link sets from the pairing map: connected groups of 'paired' sequences of exactly two (larger
    groups are reported only: a bipartition cannot separate three copies of a chromosome)."""
    parent = list(range(C))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for r in rows:
        if r["paired"]:
            parent[find(r["i"])] = find(r["j"])
    comp = defaultdict(list)
    for x in range(C):
        comp[find(x)].append(x)
    return [sorted(v) for v in comp.values() if len(v) == 2], [sorted(v) for v in comp.values() if len(v) > 2]


def mask_pairs(W, pairs):
    """Link matrix without the links between paired sequences: shared insertions join homologues whatever their
    subgenome, so those links carry no subgenome signal."""
    if not pairs:
        return W
    W = W.copy()
    for x, y in pairs:
        W[x, y] = W[y, x] = 0.0
    return W


def pair_splits(src, lab, cl_pairs, rng, support):
    """A final group holding >= 2 paired sequence pairs cannot be one subgenome (a subgenome carries one copy of each
    chromosome): split it with every such pair separated, on pooled links without the paired links, and keep the split
    when half-libraries reproduce it (adjusted Rand >= MIN_REP, more pairs when borderline; v6). Returns the new labels
    and the tests; the support of split sequences is updated in place."""
    lab, tests, k_next = lab.copy(), [], int(lab.max()) + 1
    for k in range(int(lab.max()) + 1):
        g = np.flatnonzero(lab == k)
        inside = [(x, y) for x, y in cl_pairs if lab[x] == k and lab[y] == k]
        if len(inside) < 2:
            continue
        pos = {x: i for i, x in enumerate(g)}
        sets = [[pos[x], pos[y]] for x, y in inside]
        Wm = mask_pairs(src.W["all"], inside)
        hm = [(mask_pairs(a, inside), mask_pairs(b, inside)) for a, b in src.pairs["all"]]
        z = bipartition(Wm, g, rng, sets)
        rep, se, hz = rep_stats(hm, g, rng, sets)
        while len(hm) < MAX_PAIRS and abs(rep - MIN_REP) < 2 * se:
            hm += [(mask_pairs(a, inside), mask_pairs(b, inside)) for a, b in src.more("all", HALF_PAIRS)]
            rep, se, hz = rep_stats(hm, g, rng, sets)
        ok = rep >= MIN_REP
        c = dict(src="all", z=z, con=split_contrast(Wm, g, z), rep=rep, se=se, n_pairs=len(hm),
                 borderline=abs(rep - MIN_REP) < 2 * se, hz=hz)
        tests.append(dict(members=g, path=f"paired{k}", cands={"all": c}, parent=None, chosen=c if ok else None,
                          how="pairs" if ok else "pairs: not reproduced", accepted=ok, forced=False, inside=inside))
        if ok:
            lab[g[z == 1]] = k_next
            k_next += 1
            support[g] = np.fmin(support[g], side_support(z, hz))
    return lab, tests


# ------------------------------------------------------------------------------------------------ homoeologs (v9)
# Homoeologous chromosomes share the insertions their common ancestor carried before the progenitors split: such an
# old copy and its ortholog are each other's closest relative across the two subgenomes, at collinear positions. Most
# copies are younger (subgenome-specific, spread over every chromosome of their subgenome) and say nothing about which
# chromosome pairs with which; the inherited ones are deleted over time (TE turnover), so old splits leave too few.
from functools import lru_cache  # noqa: E402

MIN_ANCHORS = 3      # v9: shared insertions that make a chromosome a candidate homoeolog (a chain needs >= 3 points)


def _old_best_block(e01):
    """Closest relative in every group of copies e0..e1-1 (one block of the old-copy self-join); run in pool workers."""
    e0, e1 = e01
    AT, DA, grp, N = (_G[x] for x in ("h_AT", "h_DA", "h_grp", "h_N"))
    R = (AT[e0:e1] @ DA).tocoo()
    r, c, v = R.row.astype(np.int64), R.col.astype(np.int64), R.data
    g = grp[c]
    k = (c != e0 + r) & (g >= 0)
    r, c, v, g = r[k], c[k], v[k], g[k]
    best = np.full((e1 - e0, N), -1, np.int64)
    if r.size:
        key = r * N + g
        o = np.lexsort((c, -v, key))                           # most related first; ties -> the lower copy index
        key, c = key[o], c[o]
        first = np.r_[True, key[1:] != key[:-1]]
        best.ravel()[key[first]] = c[first]
    return e0, best


def old_relatives(cl, old, grp, N, threads=1, block_nnz=2e7):
    """For every copy in old: its closest relative among the old copies of each group (-1: no shared lineage).
    Relatedness = lineages of <= NMAX_KNN carriers in the whole library shared by the two copies, each weighted by its
    specificity mult/(size-1) (top_relatives' measure), counted over old copies only. Blocks of copies run in the fork
    pool (scipy only, numba-free)."""
    n = old.size
    best = np.full((n, N), -1, np.int64)
    keep = np.repeat(cl.size <= NMAX_KNN, cl.size) & old[cl.members]
    sid = np.repeat(np.arange(cl.n), cl.size)[keep]
    mem = cl.members[keep].astype(np.int64)
    ok = np.bincount(sid, minlength=cl.n)[sid] >= 2                 # lineages carried by >= 2 old copies
    sid, mem = sid[ok], mem[ok]
    if not mem.size:
        return best
    u, row = np.unique(sid, return_inverse=True)
    w = cl.mult[u] / np.maximum(cl.size[u] - 1, 1)
    A = sparse.csr_matrix((np.ones(mem.size), (row, mem)), shape=(u.size, n))
    DA = sparse.csr_matrix((w[row], (row, mem)), shape=(u.size, n))
    cost = np.bincount(mem, weights=np.bincount(row)[row].astype(float), minlength=n)
    cum = np.concatenate(([0], np.cumsum(cost)))
    blocks, e0 = [], 0
    while e0 < n:
        e1 = min(max(int(np.searchsorted(cum, cum[e0] + block_nnz, side="right")), e0 + 1), n)
        blocks.append((e0, e1))
        e0 = e1
    _G.update(h_AT=A.T.tocsr(), h_DA=DA, h_grp=grp, h_N=N)
    for e0, b in pool_map(_old_best_block, blocks, threads):
        best[e0:e0 + len(b)] = b
    for x in ("h_AT", "h_DA", "h_grp", "h_N"):
        _G.pop(x, None)
    return best


def rbh_anchors(best, grp):
    """Reciprocal closest relatives across groups -> pairs (a, b), grp[a] < grp[b]: b is a's closest relative in b's
    group and a is b's in a's group (candidate orthologous insertions)."""
    out = [np.zeros((0, 2), np.int64)]
    for g in range(best.shape[1]):
        for h in range(g + 1, best.shape[1]):
            a = np.flatnonzero((grp == g) & (best[:, h] >= 0))
            b = best[a, h]
            keep = best[b, g] == a
            out.append(np.stack([a[keep], b[keep]], 1))
    return np.concatenate(out)


def _partitions(n, largest=None):
    largest = n if largest is None else largest
    if n == 0:
        yield ()
        return
    for k in range(min(n, largest), 0, -1):
        for rest in _partitions(n - k, k):
            yield (k,) + rest


@lru_cache(maxsize=None)
def lis_tail(n):
    """P(longest increasing subsequence >= l) of a uniformly random permutation of n, for l = 0..n+1, exact: the
    permutations with LIS = l number the sum, over partitions lambda of n with lambda_1 = l, of (f^lambda)^2 (RSK),
    f^lambda by the hook-length formula."""
    cnt = [0] * (n + 2)
    for lam in _partitions(n):
        conj = [sum(1 for r in lam if r > c) for c in range(lam[0])] if lam else []
        hooks = 1
        for i, r in enumerate(lam):
            for j in range(r):
                hooks *= r - j + conj[j] - i - 1
        cnt[lam[0] if lam else 0] += (math.factorial(n) // hooks) ** 2
    tail, acc = [0.0] * (n + 2), 0
    for l in range(n, -1, -1):
        acc += cnt[l]
        tail[l] = acc / math.factorial(n)
    return tuple(tail)


def chain_test(x, y, length, alpha, rng):
    """p of a collinear chain at least `length` long among n anchors (x = 0..n-1 along one chromosome, y their ranks
    along the other). n <= CHAIN_W: the window does not bind, so the chain is the longest monotone subsequence, whose
    tail is exact for one orientation and doubled for the two (a union bound). Larger n: B = ceil(2 / alpha)
    permutations of y, p = (longer-or-equal chains + 1) / (permutations + 1), stopped once p > alpha is certain (v6's
    normal approximation of a 50-permutation z ran 3-7x too liberal on simulated non-homoeologs with 26-60 anchors)."""
    n = x.size
    if n <= CHAIN_W:
        return min(1.0, 2.0 * lis_tail(n)[min(length, n + 1)])
    B, k = int(math.ceil(2.0 / alpha)), 0
    for b in range(1, B + 1):
        k += _chain_len(x, rng.permutation(y), CHAIN_W) >= length
        if (k + 1) / (B + 1) > alpha:
            return (k + 1) / (b + 1)
    return (k + 1) / (B + 1)


def run_homoeologs(a, lib, ph):
    """Homoeologous chromosomes from shared old insertions (module notes). Per chromosome i and other subgenome h: the
    candidate whose anchors with i form the most collinear chain (z against NULL_PERMS permutations) is called i's
    homoeolog when it also holds the most anchors relative to expectation (degree-corrected within the subgenome pair)
    and its chain passes chain_test at Bonferroni 0.05 over the chromosomes of h; else 'unresolved'; 'none' without a
    candidate (>= MIN_ANCHORS anchors). None for one group."""
    N, C = ph["N"], ph["C"]
    edges = ph.get("strata_edges")
    if N < 2 or not edges:
        return None
    z, ek = np.asarray(ph["z"]), np.asarray(ph["ek"])
    old = (ek >= 0) & np.isfinite(lib.k2p) & (lib.k2p > edges[1])
    if VERBOSE and ph.get("strata_n", {}).get("old") not in (None, int(old.sum())):
        log(f"WARNING: old copies {int(old.sum())} != the phasing step's old tercile {ph['strata_n']['old']}")
    grp = np.where(old, z[np.maximum(ek, 0)], -1)
    P = rbh_anchors(old_relatives(ph["clades_all"], old, grp, N, a.threads), grp)
    si, sj = ek[P[:, 0]], ek[P[:, 1]]
    M = np.zeros((C, C))
    np.add.at(M, (si, sj), 1)
    M = M + M.T
    E = np.zeros((C, C))
    for g in range(N):
        for h in range(N):
            ia, ib = np.flatnonzero(z == g), np.flatnonzero(z == h)
            if g != h and ia.size and ib.size:
                B = M[np.ix_(ia, ib)]
                E[np.ix_(ia, ib)] = np.outer(B.sum(1), B.sum(0)) / max(B.sum(), 1e-300)
    S = np.log((M + 1) / (E + 1))
    rng = np.random.default_rng([a.seed, 9])                   # own stream: no v8 draw changes
    key = np.minimum(si, sj) * C + np.maximum(si, sj)
    o = np.argsort(key, kind="stable")
    st = np.r_[0, np.flatnonzero(np.diff(key[o])) + 1] if o.size else np.zeros(0, np.int64)
    en = np.r_[st[1:], o.size] if o.size else np.zeros(0, np.int64)
    pair, in_chain = {}, np.zeros(len(P), bool)
    for s0, e0 in zip(st, en):
        idx = o[s0:e0]
        if idx.size < MIN_ANCHORS:
            continue
        i, j = int(key[idx[0]] // C), int(key[idx[0]] % C)
        ci = np.where(si[idx] == i, P[idx, 0], P[idx, 1])                 # anchor copies on i and on j
        cj = np.where(si[idx] == i, P[idx, 1], P[idx, 0])
        oi = np.argsort(lib.start[ci], kind="stable")
        idx, ci, cj = idx[oi], ci[oi], cj[oi]
        x = np.arange(idx.size, dtype=np.int64)
        y = np.argsort(np.argsort(lib.start[cj], kind="stable"), kind="stable").astype(np.int64)
        ch = longest_chain(x, y)
        null = np.array([_chain_len(x, rng.permutation(y), CHAIN_W) for _ in range(NULL_PERMS)], float)
        zc = (ch.size - null.mean()) / max(null.std(), 0.5)
        in_chain[idx[ch]] = True
        span = lambda c_, s_: float((lib.end[c_].max() - lib.start[c_].min() + 1) / max(lib.seq_len[s_], 1))
        rec = dict(n=int(idx.size), chain=int(ch.size), z=float(zc), x=x, y=y,
                   orientation="+" if ch.size < 2 or y[ch[-1]] > y[ch[0]] else "-")
        si_i, si_j = lib.seq[ci[0]], lib.seq[cj[0]]
        pair[(i, j)] = dict(rec, span=span(ci[ch], si_i), partner_span=span(cj[ch], si_j))
        pair[(j, i)] = dict(rec, span=span(cj[ch], si_j), partner_span=span(ci[ch], si_i))
    calls, tested = {}, {}
    for i in range(C):
        for h in range(N):
            if h == z[i]:
                continue
            js = np.flatnonzero(z == h)
            zs = np.array([pair[(i, j)]["z"] if (i, j) in pair else -np.inf for j in js])
            if not js.size or not np.isfinite(zs).any():
                calls[(i, h)] = dict(status="none", partner=-1, candidate=-1, p=float("nan"))
                continue
            jz, jo = int(js[np.argmax(zs)]), int(js[np.argmax(S[i, js])])
            alpha = 0.05 / js.size
            r = pair[(i, jz)]
            if (min(i, jz), max(i, jz), alpha) not in tested:
                tested[(min(i, jz), max(i, jz), alpha)] = chain_test(r["x"], r["y"], r["chain"], alpha, rng)
            p = tested[(min(i, jz), max(i, jz), alpha)]
            ok = jz == jo and p <= alpha
            calls[(i, h)] = dict(status="homoeolog" if ok else "unresolved", partner=jz if ok else -1, candidate=jz,
                                 p=float(p))
    for (i, h), c in calls.items():
        j = c["partner"]
        c["reciprocal"] = (None if j < 0 else
                           calls.get((j, int(z[i])), {}).get("partner") == i)
    sets = homoeolog_sets(calls, z, N, C)
    called = {(i, c["partner"]) for (i, h), c in calls.items() if c["partner"] >= 0}
    pair_called = np.array([(int(a_), int(b_)) in called or (int(b_), int(a_)) in called for a_, b_ in zip(si, sj)],
                           bool)
    out = dict(P=P, si=si, sj=sj, M=M, E=E, S=S, pair=pair, calls=calls, sets=sets, in_chain=in_chain,
               pair_called=pair_called, old_edge=float(edges[1]), n_old=int(old.sum()))
    for g in range(N):
        for h in range(g + 1, N):
            ng, nh = int((z == g).sum()), int((z == h).sum())
            cg = sum(calls[(i, h)]["status"] == "homoeolog" for i in np.flatnonzero(z == g))
            chh = sum(calls[(j, g)]["status"] == "homoeolog" for j in np.flatnonzero(z == h))
            rc = sum(bool(calls[(i, h)]["reciprocal"]) for i in np.flatnonzero(z == g))
            log(f"homoeologs {group_label(g)}-{group_label(h)}: {cg} of {ng} {group_label(g)} and {chh} of {nh} "
                f"{group_label(h)} sequences paired ({rc} reciprocal pairs)")
    log(f"homoeologs: {len(P):,} shared old insertions (reciprocal closest relatives among {int(old.sum()):,} copies "
        f"with LTR K2P > {edges[1]:.4g}); {len(sets)} homoeolog sets, {sum(s_['complete'] for s_ in sets)} complete")
    return out


def homoeolog_sets(calls, z, N, C):
    """Connected components of the homoeolog calls (either direction), ordered by their first member: members, and
    complete = one member per subgenome with every cross pair called."""
    parent = list(range(C))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    edges = {(min(i, c["partner"]), max(i, c["partner"])) for (i, h), c in calls.items() if c["partner"] >= 0}
    for i, j in edges:
        parent[find(i)] = find(j)
    comp = defaultdict(list)
    for i, j in edges:
        comp[find(i)].append(i)
        comp[find(i)].append(j)
    out = []
    for members in comp.values():
        m = sorted(set(members))
        complete = (sorted(int(z[x]) for x in m) == list(range(N)) and
                    all((min(x, y), max(x, y)) in edges for x in m for y in m if x < y))
        out.append(dict(members=m, complete=bool(complete)))
    return sorted(out, key=lambda s_: (min(int(z[x]) for x in s_["members"]), s_["members"]))


# ------------------------------------------------------------------------------------------------ homoeolog biology
# v9.1: what the homoeolog pairs say about the genome. Each output passed a bar fixed before it was scored
# (LTR_phasing/dev_v91: 147 simulations with exact truth, 36 real polyploids with gene-synteny and exchange truth).
HB_ALPHA = 0.05      # v9.1: family-wise level of every new test (Bonferroni within its family)
HB_MIN_RUN = 2       # v9.1: an exchanged stretch needs >= 2 consecutive collinear pairs (one pair can be a single
                     # converted or misassembled copy)
HL_FACTOR = 1.5      # v9.1: copies older than 1.5 x the split are pre-split (retention levels off there in simulations)
HL_MIN_R = 0.06      # v9.1: retention below 6% leaves the half-life to noise (simulations)
HL_PLATEAU = 0.35    # v9.1: retention at 2 x the split within 35% of retention at 1.5 x: pre-split copies behave alike
HL_BOOT = 200        # v9.1: bootstrap replicates (chromosomes; inherited pairs behind the split divergence)
SPLIT_MIN = 20       # v9.1: inherited pairs needed for the split divergence the half-life is measured against


def flag_homogenized(k2p, sites, limit=0.0, alpha=HB_ALPHA, it=30):
    """Anchors too similar to be inherited from the progenitors' split -> (flags, threshold K2P). Empirical null on
    log divergence (LTR rates vary, so a binomial is too narrow): centre = median, scale = spread of the upper half
    (homogenized pairs cannot reach it); flagged below centre - z(alpha / n) x scale, or at/below the merger's upper
    limit (K2P; 0 = unknown); iterated on the unflagged. NaN divergences are never flagged."""
    k = np.nan_to_num(np.asarray(k2p, float), nan=-1.0)
    x = np.log(np.maximum(np.where(k >= 0, k, np.nan), 0.5 / np.maximum(np.asarray(sites, float), 1.0)))
    ok = np.isfinite(x)
    n = int(ok.sum())
    if n < 5:
        return np.zeros(x.size, bool), float("nan")
    zq = NormalDist().inv_cdf(1 - alpha / n)
    keep = ok & (k > limit)
    if keep.sum() < 5:
        keep = ok.copy()
    thr = float("nan")
    for _ in range(it):
        centre = float(np.median(x[keep]))
        up = x[keep][x[keep] >= centre]
        thr = centre - zq * max(float(np.median(up) - centre) / 0.6745, 1e-3)
        new = ok & (x > thr) & (k > limit)
        if new.sum() < 5 or np.array_equal(new, keep):
            break
        keep = new
    return ok & ~((x > thr) & (k > limit)), float(math.exp(thr))


def runs_of_flags(flags, min_run=HB_MIN_RUN):
    """Runs of >= min_run consecutive flags in chain order -> [(first, last, p)], p = chance of a run at least this
    long anywhere in the chain given the number of flags (Poisson approximation for runs)."""
    f = np.asarray(flags, bool)
    n = f.size
    q = f.sum() / n if n else 0.0
    out, i = [], 0
    while i < n:
        if not f[i]:
            i += 1
            continue
        j = i
        while j + 1 < n and f[j + 1]:
            j += 1
        L = j - i + 1
        if L >= min_run:
            lam = (n - L + 1) * (1 - q) * q ** L + q ** L
            out.append((i, j, float(1 - math.exp(-lam))))
        i = j + 1
    return out


def _merger_k2p(dt, g, h):
    """Upper limit of the merger of groups g and h (K2P; 0 when v8 could not date or bound it)."""
    for r in (dt or {}).get("pairs", []):
        if set(int(x) for x in r["groups"]) == {g, h}:
            mg = r.get("merger")
            return float(mg["hi"]) if mg and mg.get("datable") and mg.get("hi") is not None else 0.0
    return 0.0


def run_homoeolog_biology(a, lib, ph, dt, pa):
    """v9.1 (module notes V9.1): divergence of every shared old insertion; near-identical stretches on called
    homoeologs (exchanges after the split); half-life of intact LTR-RTs in each direction and the two lineages
    compared; blocks collinear with a non-homoeolog (translocations) and sequences paired with several homoeologs;
    homoeologs collinear in reverse. None without homoeolog calls. Reads ph, dt and pa; changes nothing."""
    ho = ph.get("homoeologs")
    if not ho or not len(ho["P"]) or not any(c["partner"] >= 0 for c in ho["calls"].values()):
        return None
    N, C, names = ph["N"], ph["C"], ph["names"]
    z, ek = np.asarray(ph["z"]), np.asarray(ph["ek"])
    P, si, sj = ho["P"], np.asarray(ho["si"]), np.asarray(ho["sj"])
    rng = np.random.default_rng([a.seed, 91])                  # own stream: no earlier draw changes
    pair_counts(lib, int(P[0, 0]), int(P[0, 1]))                # alignment kernels compiled before the pool forks
    cnt = align_pairs(lib, P, a.threads)
    k2p, _ = k2p_counts(cnt[:, 0], cnt[:, 1], cnt[:, 2])
    k2p = np.where(cnt[:, 0] >= 200, k2p, np.nan)
    hot = ho["in_chain"] & ho["pair_called"]
    flags = np.zeros(len(P), bool)
    mu = dt["mu"] if dt else (a.mu or MU_DEFAULT)
    key = np.minimum(si, sj) * C + np.maximum(si, sj)
    gp = np.minimum(z[si], z[sj]) * N + np.maximum(z[si], z[sj])
    pairs = {}
    for g in range(N):
        for h in range(g + 1, N):
            idx = np.flatnonzero(hot & (gp == g * N + h))
            lim = _merger_k2p(dt, g, h)
            f, thr = flag_homogenized(k2p[idx], cnt[idx, 0], lim)
            flags[idx] = f
            use = idx[np.isfinite(k2p[idx]) & ~f]
            pairs[(g, h)] = dict(limit=lim, threshold=thr, use=use, n_flagged=int(f.sum()),
                                 d_s=float(np.median(k2p[use])) if use.size >= SPLIT_MIN else float("nan"))
    called = sorted({(min(i, c["partner"]), max(i, c["partner"])) for (i, h), c in ho["calls"].items()
                     if c["partner"] >= 0})
    v8 = [x for x in (pa or {}).get("segs", []) if x["called"]]
    segs = []
    for i, j in called:
        q = np.flatnonzero(hot & (key == i * C + j))
        ci, cj = np.where(si[q] == i, P[q, 0], P[q, 1]), np.where(si[q] == i, P[q, 1], P[q, 0])
        o = np.argsort(lib.start[ci], kind="stable")
        fin = np.isfinite(k2p[q[o]])
        q, ci, cj = q[o][fin], ci[o][fin], cj[o][fin]
        for s0, s1, p in runs_of_flags(flags[q]):
            if p > HB_ALPHA / len(called):
                continue
            r = slice(s0, s1 + 1)
            seg = dict(i=i, j=j, start_i=int(lib.start[ci[r]].min()), end_i=int(lib.end[ci[r]].max()),
                       start_j=int(lib.start[cj[r]].min()), end_j=int(lib.end[cj[r]].max()), n=s1 - s0 + 1,
                       n_chain=int(q.size), k2p=float(np.median(k2p[q[r]])), p=p)
            lim = pairs[(min(z[i], z[j]), max(z[i], z[j]))]["limit"]
            seg["timing"] = ("after_merger" if seg["k2p"] <= lim else "before_merger") if lim > 0 else "undetermined"
            seg["v8"] = [f"exchange_{k + 1}" for k, x in enumerate(v8)
                         if (x["seq"] == names[i] and x["start"] <= seg["end_i"] and seg["start_i"] <= x["end"]) or
                         (x["seq"] == names[j] and x["start"] <= seg["end_j"] and seg["start_j"] <= x["end"])]
            segs.append(seg)
    # half-life per direction g -> h (h's loss rate, given g's copy survived)
    partner = {(i, h): c["partner"] for (i, h), c in ho["calls"].items() if c["partner"] >= 0}
    sidx = {s_: k_ for k_, s_ in enumerate(lib.seqnames)}
    Lph = np.array([lib.seq_len[sidx[n_]] for n_ in names], float)
    ncop = np.bincount(ek[ek >= 0], minlength=C)
    density = {g: float(ncop[z == g].sum() / max(Lph[z == g].sum() / 1e6, 1e-9)) for g in range(N)}
    fin_k = np.isfinite(lib.k2p)
    by_seq = defaultdict(list)
    for c_ in np.flatnonzero((ek >= 0) & fin_k):
        by_seq[int(ek[c_])].append(int(c_))
    hl, edge = [], float(ho["old_edge"])
    for (g0, h0), pr in pairs.items():
        ds = pr["d_s"]
        if not np.isfinite(ds):
            continue
        thr = max(pr["limit"], pr["threshold"] if np.isfinite(pr["threshold"]) else 0.0)
        inh = np.flatnonzero((gp == g0 * N + h0) & np.isfinite(k2p) & (k2p > thr))
        links = defaultdict(lambda: defaultdict(int))
        for q in inh:
            links[int(P[q, 0])][int(sj[q])] += 1
            links[int(P[q, 1])][int(si[q])] += 1
        chain_k = k2p[pr["use"]]
        for g, h in ((g0, h0), (h0, g0)):
            chroms = [i for i in range(C) if z[i] == g and (i, h) in partner]
            nB = int((z == h).sum())
            if not chroms or nB < 2:
                continue
            per = {}
            for i in chroms:
                j = partner[(i, h)]
                cs = np.array(by_seq[i], np.int64)
                ks = lib.k2p[cs]
                pp = np.array([links[int(c_)].get(j, 0) for c_ in cs], float)
                qq = np.array([sum(n_ for t_, n_ in links[int(c_)].items() if z[t_] == h and t_ != j) for c_ in cs],
                              float)
                o = np.argsort(-ks)
                per[i] = (-ks[o], np.r_[0, np.cumsum(pp[o])], np.r_[0, np.cumsum(qq[o])])

            def ret(t, chs):
                obs = paired = other = 0.0
                for i in chs:
                    nk, cp, cq = per[i]
                    m = int(np.searchsorted(nk, -t, side="left"))
                    obs, paired, other = obs + m, paired + cp[m], other + cq[m]
                return (paired - other / max(nB - 1, 1)) / max(obs, 1.0), obs
            r15, n15 = ret(max(edge, HL_FACTOR * ds), chroms)
            r2, _ = ret(max(edge, 2 * ds), chroms)
            est = math.log(2) * ds / -math.log(r15) if 0 < r15 < 1 else float("nan")
            boot = np.full(HL_BOOT, np.nan)
            for b in range(HL_BOOT):
                chs = list(rng.choice(chroms, len(chroms), replace=True))
                dsb = float(np.median(rng.choice(chain_k, chain_k.size, replace=True)))
                rb = ret(max(edge, HL_FACTOR * dsb), chs)[0]
                boot[b] = math.log(2) * dsb / -math.log(rb) if 0 < rb < 1 else np.nan
            good = np.isfinite(boot).sum() > HL_BOOT // 2
            lo, up = np.nanquantile(boot, [0.025, 0.975]) if good else (np.nan, np.nan)
            qc = dict(plateau=r15 > 0 and abs(r2 - r15) <= HL_PLATEAU * r15, signal=bool(good), range=r15 >= HL_MIN_R)
            hl.append(dict(g=g, h=h, d_s=ds, r15=float(r15), r2=float(r2), n_old=int(n15), est=est, lo=float(lo),
                           hi=float(up), ok=all(qc.values()), qc=qc, boot=boot))
    asym = []
    for e1 in hl:
        e2 = next((e for e in hl if e["g"] == e1["h"] and e["h"] == e1["g"]), None)
        if e2 is None or e1["g"] > e1["h"] or not (e1["ok"] and e2["ok"]):
            continue
        m = np.isfinite(e1["boot"]) & np.isfinite(e2["boot"]) & (e1["boot"] > 0) & (e2["boot"] > 0)
        lr = np.log(e1["boot"][m] / e2["boot"][m])
        lo, up = np.quantile(lr, [0.025, 0.975]) if lr.size else (np.nan, np.nan)
        asym.append(dict(a=e1["h"], b=e2["h"], ratio=e1["est"] / e2["est"], lo=float(np.exp(lo)), hi=float(np.exp(up)),
                         differ=bool(lr.size) and not (lo <= 0 <= up)))
    # translocated blocks: a collinear chain between i and a chromosome other than its homoeolog
    pset = defaultdict(set)
    for (i, h), j in partner.items():
        pset[i].add(j)
        pset[j].add(i)
    called_set = set(called)
    o = np.argsort(key, kind="stable")
    st = np.r_[0, np.flatnonzero(np.diff(key[o])) + 1] if o.size else np.zeros(0, np.int64)
    tests = []
    for s0, e0 in zip(st, np.r_[st[1:], o.size]):
        q = o[s0:e0]
        i, j = int(key[q[0]] // C), int(key[q[0]] % C)
        if q.size < 3 or (i, j) in called_set or not (pset[i] or pset[j]):
            continue
        ci, cj = np.where(si[q] == i, P[q, 0], P[q, 1]), np.where(si[q] == i, P[q, 1], P[q, 0])
        oi = np.argsort(lib.start[ci], kind="stable")
        ci, cj = ci[oi], cj[oi]
        x = np.arange(ci.size, dtype=np.int64)
        y = np.argsort(np.argsort(lib.start[cj], kind="stable"), kind="stable").astype(np.int64)
        tests.append((i, j, ci, cj, x, y, longest_chain(x, y)))
    tl, alpha = [], HB_ALPHA / max(len(tests), 1)
    for i, j, ci, cj, x, y, ch in tests:
        if ch.size < 3:
            continue
        p = chain_test(x, y, int(ch.size), alpha, rng)
        if p <= alpha:
            tl.append(dict(i=i, j=j, start_i=int(lib.start[ci[ch]].min()), end_i=int(lib.end[ci[ch]].max()),
                           start_j=int(lib.start[cj[ch]].min()), end_j=int(lib.end[cj[ch]].max()), n=int(ch.size),
                           n_tested=int(x.size), p=float(p)))
    many = []
    for s_, ps in sorted(pset.items()):
        for h in sorted({int(z[p_]) for p_ in ps}):
            mem = sorted((p_ for p_ in ps if z[p_] == h), key=lambda p_: natural_key(names[p_]))
            if len(mem) >= 2:
                many.append(dict(i=s_, h=h, partners=mem))
    rev = [(i, j) for i, j in called if ho["pair"].get((i, j), {}).get("orientation") == "-"]
    log(f"homoeolog biology: {len(segs)} near-identical stretches, {sum(e['ok'] for e in hl)} of {len(hl)} half-lives "
        f"estimable, {len(tl)} translocated blocks, {len(many)} sequences with several homoeologs, {len(rev)} "
        f"homoeolog pairs collinear in reverse")
    return dict(k2p=k2p, sites=cnt[:, 0], flags=flags, pairs=pairs, segments=segs, halflife=hl, asym=asym,
                density=density, transloc=tl, many=many, reversed=rev, mu=mu, n_called=len(called))


# ------------------------------------------------------------------------------------------------ phasing pipeline
def verdict_text(N, inside, big, tests, source):
    """One-line verdict. inside: paired sequence pairs left in one final group; big: groups of >= 3 near-identical
    sequences (paired map components)."""
    if N >= 2:
        v = "subgenomes found" + (" (borderline)" if any(t["accepted"] and t["chosen"]["borderline"] for t in tests)
                                  else "")
        nf = sum(t["forced"] for t in tests)
        if nf:
            v += f"; {nf} of the requested splits fail the split rule (forced by -n)"
        if inside >= 2:
            v += f"; {inside} pairs of near-identical sequences share a group: more subgenomes than resolved"
        if big:
            v += f"; {len(big)} groups of 3 or more near-identical sequences"
        return v
    if inside >= 2 or big:
        return (f"cannot tell: near-identical copies of {inside + len(big)} chromosomes (homologues of a young "
                f"polyploid, or both haplotypes of a diploid assembly); LTR-RT lineages do not separate them")
    root = tests[0]["cands"].get("all") if tests else None
    if source == "given" or (root is not None and not root["borderline"]):
        return "no structure detected"
    return "cannot tell"


def group_label(g):
    return f"SG{g + 1}"


def read_assign(path, seqnames):
    """--assign file: sequence<whitespace>label per line ('#' comments; further columns ignored) -> {sequence: label}
    in file order. Dies on conflicting labels for one sequence, on a file naming no sequence of the library, and on
    more than MAX_N groups; one-token lines and sequences the library does not hold are reported and skipped; a
    single group is warned about (nothing to date or paint)."""
    out, bad, conflict = {}, 0, []
    with open(path) as f:
        for line in f:
            x = line.split("#", 1)[0].split()
            if len(x) == 1:
                bad += 1
            elif len(x) >= 2:
                if x[0] in out and out[x[0]] != x[1]:
                    conflict.append(f"{x[0]} ({out[x[0]]}, {x[1]})")
                out.setdefault(x[0], x[1])
    if conflict:
        die(f"--assign {path}: {len(conflict)} sequence(s) given two labels: " + ", ".join(conflict[:5]))
    if not out:
        die(f"--assign {path}: no 'sequence label' lines")
    if bad:
        log(f"WARNING: --assign: {bad} line(s) with a single field skipped")
    known = set(seqnames)
    miss = [s for s in out if s not in known]
    if len(miss) == len(out):
        die(f"--assign {path}: none of its sequences is in the library (file: {', '.join(miss[:3])}; library: "
            f"{', '.join(list(seqnames)[:3])})")
    if miss:
        log(f"WARNING: --assign: {len(miss)} sequence(s) not in the library (skipped): " + ", ".join(miss[:10]))
    out = {s: g for s, g in out.items() if s in known}
    n = len(set(out.values()))
    if n > MAX_N:
        die(f"--assign {path}: {n} groups; at most {MAX_N} (is the second column a subgenome label?)")
    if n == 1:
        log("WARNING: --assign: one group only; nothing to date or paint")
    return out


def realized_groups(z, N):
    """Labels renumbered 0..k-1 over the non-empty groups, and k (a forced N can leave groups empty)."""
    u, z2 = np.unique(z, return_inverse=True)
    return z2, int(u.size)


def _shown_candidate(t):
    """The candidate a tested split is reported by: the chosen one, else the pooled one, else the most replicable."""
    if t["chosen"] is not None:
        return t["chosen"]
    if "all" in t["cands"]:
        return t["cands"]["all"]
    return max(t["cands"].values(), key=lambda c: (round(c["rep"], 2), c["con"])) if t["cands"] else None


def tree_row(t, names, sources):
    """structure.tsv / summary.json record of one tested split: the reported candidate and every source's numbers."""
    c = _shown_candidate(t)
    g = t["members"]
    row = dict(level=t["path"], n=int(len(g)), copies=c["src"] if c else None, rule=t["how"],
               accepted=bool(t["accepted"]), forced=bool(t["forced"]),
               part_a=[names[i] for i in g[c["z"] == 0]] if c else [names[i] for i in g],
               part_b=[names[i] for i in g[c["z"] == 1]] if c else [],
               replicability=c["rep"] if c else None, replicability_se=c["se"] if c else None,
               n_pairs=c["n_pairs"] if c else 0, borderline=bool(c["borderline"]) if c else False,
               strength=c["con"] if c else None,
               parent_strength=t["parent"]["con"][c["src"]] if c and t["parent"] else None)
    for s_ in sources:
        x = t["cands"].get(s_)
        row[f"replicability_{s_}"] = x["rep"] if x else None
        row[f"strength_{s_}"] = x["con"] if x else None
        row[f"agrees_{s_}"] = (bool(ari(x["z"], c["z"]) >= MIN_REP) if x and c else None)
        row[f"parent_strength_{s_}"] = t["parent"]["con"][s_] if t["parent"] else None
    return row


def split_log(r):
    """One log line per tested split: replicability/strength on each source (parent strength in brackets)."""
    src = [k[len("replicability_"):] for k in r if k.startswith("replicability_") and k != "replicability_se"]
    parts = []
    for s_ in src:
        if r[f"replicability_{s_}"] is None:
            continue
        par = f" ({r[f'parent_strength_{s_}']:.2f})" if r[f"parent_strength_{s_}"] is not None else ""
        parts.append(f"{s_} {r[f'replicability_{s_}']:.2f}/{r[f'strength_{s_}']:.2f}{par}")
    verdict = (("forced by -n" if r["forced"] else "accepted" + (" (age-bounded, " + r["copies"] + " copies)"
                                                               if r["rule"] == "age" else ""))
               if r["accepted"] or r["forced"] else "rejected: " + r["rule"])
    return (f"split {len(r['part_a'])}|{len(r['part_b'])} [{r['level']}] replicability/strength: " + "; ".join(parts)
            + f" -> {verdict}")


def run_phasing(a, lib, cl):
    rng = np.random.default_rng(a.seed)
    n_el_all = np.bincount(lib.seq, minlength=len(lib.seqnames))
    user = read_assign(a.assign, lib.seqnames) if a.assign else None      # fixed partition: no structure search
    phased = [i for i in range(len(lib.seqnames)) if n_el_all[i] >= MIN_ELEMENTS and
              (user is None or lib.seqnames[i] in user)]
    if user and len(phased) < len(user):
        log(f"WARNING: --assign: {len(user) - len(phased)} sequence(s) carry < {MIN_ELEMENTS} LTR-RTs (unphased)")
    N = a.n_subgenomes
    if len(phased) < max(N or 2, 2):
        die(f"only {len(phased)} sequences carry >= {MIN_ELEMENTS} LTR-RTs; cannot form {max(N or 2, 2)} groups")
    names = [lib.seqnames[i] for i in phased]
    N_requested = N
    C = len(names)
    pos = {s: k for k, s in enumerate(phased)}
    ek = np.array([pos.get(int(s), -1) for s in lib.seq], np.int64)
    log(f"phasing {C} sequences with >= {MIN_ELEMENTS} LTR-RTs" +
        (f" ({len(lib.seqnames) - C} others listed as unphased)" if len(lib.seqnames) > C else ""))
    lm = LinkModel(cl, ek, C, a.threads)
    if lm.cl.n == 0:
        die("no shared k-mers link two sequences; is this an LTR-RT library?")
    informative = np.zeros(lib.n, bool)
    informative[lm.cl.members] = True
    n_inf = np.bincount(ek[(ek >= 0) & informative], minlength=C)
    pairing = pairing_map(lib, lm, np.random.default_rng([a.seed, 3]), a.threads)
    cl_pairs, big = paired_sets(pairing, C)
    log(f"pairing map: {sum(r['significant'] for r in pairing)} sequence pairs with collinear shared insertions, "
        f"{sum(r['paired'] for r in pairing)} paired (near-identical sequences)")
    src = Sources(lm, lib.k2p, rng)
    if src.edges:
        log(f"link sources: all copies, and LTR-K2P terciles young <= {src.edges[0]:.4g} < middle <= "
            f"{src.edges[1]:.4g} < old ({', '.join(f'{s_} {int(m.sum()):,}' for s_, m in src.masks.items())} copies)")
    W = W_raw = src.W["all"]
    n_el = np.bincount(ek[ek >= 0], minlength=C)
    tree, tests, source, use_pairs = [], [], "given" if N else "inferred", []
    support = np.full(C, np.nan)
    if user:
        here = {user[s] for s in names}
        order = [g for g in dict.fromkeys(user.values()) if g in here]    # labels in the file's order
        z, N, source = np.array([order.index(user[s]) for s in names]), len(order), "user"
        log(f"partition given by the user: {N} group(s) ({', '.join(order)} -> "
            f"{', '.join(group_label(g) for g in range(N))})")
    else:
        sets_of = None
        if N and cl_pairs:                                        # forced N: paired sequences kept apart (v6)
            use_pairs = list(cl_pairs)

            def sets_of(g):
                p_ = {x: i for i, x in enumerate(g)}
                return [[p_[x], p_[y]] for x, y in use_pairs if x in p_ and y in p_]
        lab, tests, support = find_subgenomes(src, rng, n_force=N, sets_of=sets_of)
        if N is None and cl_pairs:
            lab, ptests = pair_splits(src, lab, cl_pairs, rng, support)
            tests += ptests
            use_pairs = [p for t in ptests if t["accepted"] for p in t["inside"]]
        if use_pairs:
            log(f"paired sequences: {len(use_pairs)} pairs kept in different groups")
        for t in tests:
            if t["how"] == "too small":                            # groups of < 4 sequences are never split
                continue
            tree.append(tree_row(t, names, src.names))
            log(split_log(tree[-1]))
        z = canon_labels(lab, n_el)
        if N and int(z.max()) + 1 < N:
            log(f"WARNING: only {int(z.max()) + 1} of the {N} requested groups could be formed")
        N = int(z.max()) + 1
        log(f"{'inferred' if source == 'inferred' else 'requested'} number of subgenomes: {N}")
    if use_pairs:                                                  # paired links carry no subgenome signal
        W = mask_pairs(W_raw, use_pairs)
    halves = [mask_pairs(np.asarray(w, np.float64), use_pairs) for p in src.pairs["all"] for w in p]
    one = N == 1
    inside = int(sum(z[x] == z[y] for x, y in cl_pairs))            # paired sequences left in one final group
    verdict = "partition given by the user" if user else verdict_text(N, inside, big, tests, source)
    log(f"verdict: {verdict}")
    user_labels = {group_label(g): order[g] for g in range(len(order))} if user else None
    common = dict(names=names, C=C, ek=ek, tree=tree, W=W, W_raw=W_raw, W_src=dict(src.W), n_el=n_el, n_inf=n_inf,
                  halves=[np.asarray(h, np.float32) for h in halves] if C <= FIG_HALVES_MAX else None,   # report.pdf C
                  source=source, user_labels=user_labels, clades=lm.cl, clades_all=cl, N_requested=N_requested or N,
                  n_clades=lm.cl.n, pairing=pairing, cl_pairs=use_pairs, paired_groups=big, verdict=verdict,
                  n_pairs=len(src.pairs["all"]), strata_edges=src.edges,
                  strata_n={s_: int(m.sum()) for s_, m in src.masks.items()})
    if one:
        partner, psup, pex = partner_pairs(W, halves, z)
        return dict(common, z=z, N=1, support=np.full(C, np.nan), partners=(partner, psup, pex),
                    own_oe=np.full(C, np.nan), oth_oe=np.full(C, np.nan), pval=np.full(C, np.nan))
    if user:
        support = support_from_halves(halves, z, N, make_sets(names, None, N), rng)
    own_oe, oth_oe, pval = chrom_tests(W, z, N)
    partner, psup, pex = partner_pairs(W, halves, z)
    for g in range(N):
        log(f"{group_label(g)}: {int((z == g).sum())} sequences, {int(n_el[z == g].sum()):,} LTR-RTs")
    return dict(common, z=z, N=N, support=support, partners=(partner, psup, pex), own_oe=own_oe, oth_oe=oth_oe,
                pval=pval)


# ------------------------------------------------------------------------------------------------ dating: copy pairs
def k2p_se_copy(lib):
    """Standard error of each copy's LTR K2P (annotation value when given, else Poisson on its aligned LTR sites)."""
    n = np.where(np.isfinite(lib.n_sites) & (lib.n_sites > 0), lib.n_sites,
                 np.nan_to_num(lib.ltr5_len, nan=0.0))
    se = np.sqrt(np.maximum(lib.k2p, 1e-6) / np.maximum(n, 1.0))
    return np.where(np.isfinite(lib.k2p_se) & (lib.k2p_se > 0), lib.k2p_se, se)


def sample_pairs(lib, ek, rng, threads):
    """Copy pairs for dating: FOCAL random phased copies plus FOCAL among the youngest 20% (LTR K2P), each with its
    REL_TOP closest relatives on other phased sequences (shared 15-mer sketch hashes), aligned as whole records.
    Returns dict(P=pairs (focal first), d=pair K2P, se, stratum (0 random, 1 young))."""
    phased = (ek >= 0) & np.isfinite(lib.k2p)
    cand = np.flatnonzero(phased)
    if cand.size < 20:
        return None
    S = sketches(lib.codes, lib.offs, SKETCH_K, SKETCH_S)
    young = cand[lib.k2p[cand] <= np.quantile(lib.k2p[cand], 0.2)]
    fa = rng.choice(cand, min(FOCAL, cand.size), replace=False)
    fy = rng.choice(young, min(FOCAL, young.size), replace=False)
    focal = np.concatenate([fa, fy])
    stratum = np.concatenate([np.zeros(fa.size, np.int8), np.ones(fy.size, np.int8)])
    J, _ = top_relatives_sketch(S, np.where(phased, ek, -1).astype(np.int64), focal, REL_TOP, use=phased)
    r, q = np.nonzero(J >= 0)
    P = np.stack([focal[r], J[r, q]], 1)
    key = np.minimum(P[:, 0], P[:, 1]) * lib.n + np.maximum(P[:, 0], P[:, 1])
    _, first = np.unique(key, return_index=True)
    first.sort()
    P, st = P[first], stratum[r[first]]
    cnt = align_whole(lib, P, threads)
    d, se = k2p_counts(cnt[:, 0], cnt[:, 1], cnt[:, 2])
    ok = np.isfinite(d) & (cnt[:, 0] >= 500)                # >= 500 aligned sites (whole elements)
    return dict(P=P[ok], d=d[ok], se=se[ok], stratum=st[ok], cnt=cnt[ok], sketch=S, young=young, focal_young=fy)


def sample_young_more(lib, ek, pt, rng, threads):
    """v8: more young focal copies for the merger clock. The youngest 20% of phased copies beyond v7's FOCAL young ones,
    up to YOUNG_MAX in all (drawn on their own stream: v7's draws and table are untouched), each with its REL_TOP
    closest relatives on other sequences, aligned like v7's pairs; pairs already in v7's table are not repeated.
    LTR-RT-rich libraries thus contribute more sibling pairs (narrower intervals, younger mergers resolved); the
    budget bounds run time. Returns v7's table extended (stratum 2 = extra young focal copies)."""
    phased = (ek >= 0) & np.isfinite(lib.k2p)
    rest = np.setdiff1d(pt["young"], pt["focal_young"])
    n_extra = int(min(max(YOUNG_MAX - pt["focal_young"].size, 0), rest.size))
    out = dict(P=pt["P"], d=pt["d"], se=pt["se"], stratum=pt["stratum"], cnt=pt["cnt"])
    if n_extra == 0:
        return out
    fy2 = np.sort(rng.choice(rest, n_extra, replace=False))
    J, _ = top_relatives_sketch(pt["sketch"], np.where(phased, ek, -1).astype(np.int64), fy2, REL_TOP, use=phased)
    r, q = np.nonzero(J >= 0)
    P = np.stack([fy2[r], J[r, q]], 1)
    key = np.minimum(P[:, 0], P[:, 1]) * lib.n + np.maximum(P[:, 0], P[:, 1])
    k7 = np.minimum(pt["P"][:, 0], pt["P"][:, 1]) * lib.n + np.maximum(pt["P"][:, 0], pt["P"][:, 1])
    _, first = np.unique(key, return_index=True)
    first = np.sort(first[~np.isin(key[first], k7)])
    P = P[first]
    if not len(P):
        return out
    cnt = align_whole(lib, P, threads)
    d, se = k2p_counts(cnt[:, 0], cnt[:, 1], cnt[:, 2])
    ok = np.isfinite(d) & (cnt[:, 0] >= 500)
    return dict(P=np.concatenate([pt["P"], P[ok]]), d=np.concatenate([pt["d"], d[ok]]),
                se=np.concatenate([pt["se"], se[ok]]), stratum=np.concatenate([pt["stratum"], np.full(ok.sum(), 2, np.int8)]),
                cnt=np.concatenate([pt["cnt"], cnt[ok]]))


def ltr_scale_near(lib, P, d, cnt, est, rng, threads, n_max=2000):
    """LTR / whole-record K2P of the pairs whose whole-record divergence is within a factor of 2 of est (pooled
    counts of their 5'-5' + 3'-3' LTR alignments over those of their whole records); None with < 30 such pairs."""
    sel = np.flatnonzero((d > est / 2) & (d < est * 2))
    if sel.size < 30:
        return None
    sel = rng.choice(sel, min(sel.size, n_max), replace=False)
    lc = align_pairs(lib, P[sel], threads)
    ok = lc[:, 0] >= 200
    if ok.sum() < 30:
        return None
    dl = k2p_counts(lc[ok, 0].sum(), lc[ok, 1].sum(), lc[ok, 2].sum())[0]
    dw = k2p_counts(cnt[sel][ok, 0].sum(), cnt[sel][ok, 1].sum(), cnt[sel][ok, 2].sum())[0]
    return float(dl / dw) if np.isfinite(dl) and np.isfinite(dw) and dw > 0 else None


def ltr_rate_ratio(lib, P, cnt, rng, threads, n_max=3000):
    """LTR rate / whole-record rate of sibling copy pairs: pooled substitutions per site of their 5'-5' + 3'-3' LTR
    alignments over those of their whole-record alignments (up to n_max pairs, both with >= 200 LTR sites and some
    substitutions); 1 when undeterminable. Puts pair divergences on the LTR-age scale of the clock."""
    if not len(P):
        return 1.0
    i = rng.choice(len(P), min(len(P), n_max), replace=False)
    lc = align_pairs(lib, P[i], threads)
    ok = lc[:, 0] >= 200
    pl = (lc[ok, 1] + lc[ok, 2]).sum() / max(lc[ok, 0].sum(), 1)
    pw = (cnt[i][ok, 1] + cnt[i][ok, 2]).sum() / max(cnt[i][ok, 0].sum(), 1)
    return float(pl / pw) if pl > 0 and pw > 0 else 1.0


def exchange_locations(lib, loc, segs):
    """Copy locations for the merger clock: copies inside a called exchange segment take the segment's origin (they
    descend from that subgenome; their relatives there are not moves between subgenomes after the merger)."""
    out = loc.copy()
    si = {s_: i for i, s_ in enumerate(lib.seqnames)}
    for x in segs:
        if x["called"]:
            m = (lib.seq == si[x["seq"]]) & (lib.start >= x["start"]) & (lib.end <= x["end"]) & (loc == x["own"])
            out[m] = x["origin"]
    return out


def units_v7(pt, n):
    """Bootstrap units from v7's pair table: components of copies joined by pairs closer than CLUSTER_D (pure families
    in simulations), and which copies appear in it. Denser sampling must not chain families into one unit, so pairs
    added by v8 take the unit of their focal copy if it is in v7's table, else of their partner."""
    comp = lineage_clusters(pt["P"], pt["d"], n, CLUSTER_D)
    in7 = np.zeros(n, bool)
    in7[pt["P"].ravel()] = True
    return comp, in7


PROF_EDGES = np.linspace(-4, 0, 41)      # log10 pair-divergence bins of the cross-subgenome profile


def cross_profile(d, cross):
    """Pairs and cross-subgenome pairs per log10-divergence bin (pairs at d = 0 go to the first bin)."""
    j = np.clip(np.searchsorted(PROF_EDGES, np.log10(np.maximum(d, 1e-12))) - 1, 0, PROF_EDGES.size - 2)
    return np.bincount(j, minlength=PROF_EDGES.size - 1), np.bincount(j, weights=cross, minlength=PROF_EDGES.size - 1)


def _profile_trough(d, cross):
    """Bins (n, cross fraction f, bins with >= 20 pairs) and the trough: the bin of lowest 3-bin pooled cross
    fraction below the median divergence (the confined plateau between the merger and the split)."""
    n, k = cross_profile(d, cross)
    f = np.where(n >= 20, k / np.maximum(n, 1), np.nan)
    ok = np.flatnonzero(np.isfinite(f))
    if not ok.size:
        return n, k, f, ok, None
    w = np.ones(3)
    fs = np.convolve(k, w, "same") / np.maximum(np.convolve(n, w, "same"), 1)
    low = ok[PROF_EDGES[ok] < np.log10(max(np.median(d), 1e-12))]
    cand = low if low.size else ok
    return n, k, f, ok, int(cand[np.argmin(fs[cand])])


def young_cut(d, cross):
    """Upper edge of the young (post-merger) regime: the first bin from below whose cross fraction is within 0.02 of
    the trough."""
    n, k, f, ok, b0 = _profile_trough(d, cross)
    if b0 is None:
        return 0.0
    b = next(b for b in ok if f[b] <= f[b0] + 0.02)
    return 10 ** PROF_EDGES[b]


def plateau_end(d, cross):
    """End of the confined plateau: scanning up from the trough, bins join the plateau while their cross fraction is
    not significantly above the plateau's pooled fraction (one-sided binomial, Bonferroni over the bins scanned).
    Below it no pair needs to predate the split."""
    from scipy.stats import norm
    n, k, f, ok, b0 = _profile_trough(d, cross)
    if b0 is None:
        return float(np.median(d)) if d.size else 0.0
    up = ok[ok > b0]
    zc = float(norm.isf(0.05 / max(up.size, 1)))
    kk, nn, last = k[b0], n[b0], b0
    for b in up:
        p0 = min(max(kk / nn, 0.5 / nn), 1 - 0.5 / nn)
        if (f[b] - p0) / np.sqrt(p0 * (1 - p0) / n[b]) >= zc:
            break
        kk, nn, last = kk + k[b], nn + n[b], b
    return 10 ** PROF_EDGES[last + 1]


def _iso_increasing(y, w):
    """Pool-adjacent-violators (non-decreasing fit)."""
    blocks = []
    for v, ww in zip(map(float, y), map(float, w)):
        blocks.append([v, ww, 1])
        while len(blocks) > 1 and blocks[-2][0] > blocks[-1][0]:
            a, b = blocks.pop(), blocks.pop()
            W = a[1] + b[1]
            blocks.append([(a[0] * a[1] + b[0] * b[1]) / max(W, 1e-12), W, a[2] + b[2]])
    return np.concatenate([[v] * c for v, _, c in blocks]) if blocks else np.array([])


def half_rise(d, cross):
    """Lineage-floor split: the divergence where the (monotone) cross fraction above the young cut climbs halfway
    from the confined plateau to the ancestral level; NaN when it does not climb by >= 0.05."""
    n, k, f, ok, b0 = _profile_trough(d, cross)
    if b0 is None:
        return float("nan")
    above = ok[ok >= b0]
    iso = _iso_increasing(f[above], n[above])
    if iso.size < 2 or iso[-1] - iso[0] < 0.05:
        return float("nan")
    b = above[np.argmax(iso >= iso[0] + 0.5 * (iso[-1] - iso[0]))]
    return float(10 ** (0.5 * (PROF_EDGES[b] + PROF_EDGES[b + 1])))


def lineage_clusters(P, d, n, max_d):
    """Bootstrap units: connected components of copies joined by pairs closer than max_d (pure families in sims)."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    m = d < max_d
    A = coo_matrix((np.ones(m.sum()), (P[m, 0], P[m, 1])), shape=(n, n))
    return connected_components(A, directed=False)[1]


def widen(est, lo, hi, kind):
    """Bootstrap interval (est, lo, hi) widened by MODEL_ERR[kind] in quadrature on the log scale."""
    me = 1.96 * MODEL_ERR[kind]
    zl = np.log(est / lo) if lo > 0 else np.inf
    zu = np.log(hi / est) if est > 0 else 0.0
    return (float(est * np.exp(-np.hypot(zl, me))) if np.isfinite(zl) else 0.0), float(est * np.exp(np.hypot(zu, me)))


def _cluster_boot(clusters, rng, n_boot):
    u, inv = np.unique(clusters, return_inverse=True)
    idx = [np.flatnonzero(inv == k) for k in range(u.size)]
    for _ in range(n_boot):
        yield np.concatenate([idx[k] for k in rng.integers(0, u.size, u.size)])


def duplicated_insertions(lib, P, d, se):
    """Pairs whose divergence is significantly below both copies' own LTR divergence: one insertion duplicated
    later (inherited by two homoeologues: orthologues; or copied by an exchange or segmental duplication)."""
    ks = k2p_se_copy(lib)
    a, b = P[:, 0], P[:, 1]
    with np.errstate(invalid="ignore"):
        return (d + 2 * se < lib.k2p[a] - 2 * ks[a]) & (d + 2 * se < lib.k2p[b] - 2 * ks[b])


def split_orthologs(d, orth, clusters, rng):
    """Orthologous-insertion split: median divergence of orthologous pairs, lineage-cluster bootstrap interval."""
    v = d[orth]
    if v.size < N_ORTH:
        return None
    bs = [np.median(v[sel]) for sel in _cluster_boot(clusters[orth], rng, BOOT)]
    est = float(np.median(v))
    lo, hi = widen(est, *np.quantile(bs, [0.025, 0.975]), "orthologs")
    return dict(estimate=est, lo=lo, hi=hi, n=int(v.size))


def split_profile(d, cross, clusters, rng):
    """Lineage-floor split (profile half-rise), lineage-cluster bootstrap interval."""
    est = half_rise(d, cross)
    if not np.isfinite(est):
        return None
    bs = [half_rise(d[sel], cross[sel]) for sel in _cluster_boot(clusters, rng, BOOT)]
    bs = np.array([x for x in bs if np.isfinite(x)])
    if bs.size < BOOT // 2:
        return None
    lo, hi = widen(est, *np.quantile(bs, [0.025, 0.975]), "profile")
    return dict(estimate=est, lo=lo, hi=hi, n=int(d.size))


# ------------------------------------------------------------------------------------------------ dating: v5 split
@njit(cache=True)
def _ancestry_scores(members, indptr, mult, seq, zseq, N, T, Cn, nseq_tot):
    """S_ab, n_ab (abundance), S_br, n_br (breadth) per copy. seq: sequence of each copy (-1 unphased: scored on
    every phased carrier, nothing left out); zseq: group of each sequence; T: phased copies per group; Cn: phased
    sequences per group."""
    n = seq.size
    Sa = np.zeros((n, N))
    na = np.zeros(n)
    Sb = np.zeros((n, N))
    nb = np.zeros(n)
    cnt_seq = np.zeros(nseq_tot, np.int64)
    M = np.zeros(N)
    Mc = np.zeros(N)
    m = np.zeros(N)
    sc = np.zeros(N)
    for c in range(indptr.size - 1):
        lo = indptr[c]
        hi = indptr[c + 1]
        M[:] = 0.0
        Mc[:] = 0.0
        for p in range(lo, hi):
            s = seq[members[p]]
            if s >= 0:
                if cnt_seq[s] == 0:
                    Mc[zseq[s]] += 1.0
                cnt_seq[s] += 1
                M[zseq[s]] += 1.0
        w = mult[c]
        for p in range(lo, hi):
            i = members[p]
            s = seq[i]
            g = zseq[s] if s >= 0 else -1                  # unphased copy: scored on all phased carriers
            # abundance: other members, own sequence left out
            tot = 0.0
            for h in range(N):
                m[h] = M[h]
            if g >= 0:
                m[g] -= cnt_seq[s]
            for h in range(N):
                tot += m[h]
            if tot > 0:
                mean = 0.0
                for h in range(N):
                    sc[h] = math.log((m[h] + 0.5) / T[h])
                    mean += sc[h]
                mean /= N
                for h in range(N):
                    Sa[i, h] += w * (sc[h] - mean)
                na[i] += w
            # breadth: other sequences per group carrying the clade, own sequence left out
            mx = 0.0
            for h in range(N):
                m[h] = Mc[h]
            if g >= 0:
                m[g] -= 1.0
            for h in range(N):
                if m[h] > mx:
                    mx = m[h]
            if mx >= 2:
                smax = -1e300
                smin = 1e300
                for h in range(N):
                    avail = Cn[h] - (1.0 if h == g else 0.0)
                    sc[h] = math.log((m[h] + 0.5) / (avail + 1.0))
                    if sc[h] > smax:
                        smax = sc[h]
                    if sc[h] < smin:
                        smin = sc[h]
                if smax - smin >= math.log(3.0):
                    mean = 0.0
                    for h in range(N):
                        mean += sc[h]
                    mean /= N
                    for h in range(N):
                        Sb[i, h] += w * (sc[h] - mean)
                    nb[i] += w
        for p in range(lo, hi):
            s = seq[members[p]]
            if s >= 0:
                cnt_seq[s] = 0
    return Sa, na, Sb, nb


def _calibrate(S, n, y):
    """v5: P(o | copy) = softmax(lam * S / n^beta + b) fitted to the copies' location labels (foreign copies are rare
    label noise that only makes the fit conservative)."""
    from scipy.optimize import minimize
    from scipy.special import logsumexp, expit
    N = S.shape[1]
    ln = np.log(np.maximum(n, 1.0))
    scale = float(np.median(np.abs(S).sum(1))) + 1e-12

    def unpack(par):
        return math.exp(min(par[0], 50.0)) / scale, float(expit(par[1])), np.concatenate(([0.0], par[2:]))

    def nll(par):
        lam, beta, b = unpack(par)
        Z = lam * S * np.exp(-beta * ln)[:, None] + b[None, :]
        return -(Z[np.arange(len(y)), y] - logsumexp(Z, 1)).sum()
    best = min((minimize(nll, np.concatenate(([0.0, b0], np.zeros(N - 1))), method="L-BFGS-B")
                for b0 in (-2.0, 0.0, 2.0)), key=lambda r: r.fun)
    return unpack(best.x)


def _apply_cal(S, n, cal):
    from scipy.special import logsumexp
    lam, beta, b = cal
    Z = lam * S * np.exp(-beta * np.log(np.maximum(n, 1.0)))[:, None] + b[None, :]
    return Z - logsumexp(Z, 1, keepdims=True)


def copy_ancestry(cl, ek, z, N, rng, fit_max=200000):
    """v5 per-copy lineage ancestry on v6 clades: two experts (lineage abundance and breadth on the other sequences,
    own sequence left out), each calibrated to the location labels, averaged, recalibrated. -> log P (n x N), number
    of informative k-mers per copy."""
    own = np.where(ek >= 0, z[np.maximum(ek, 0)], -1)
    T = np.bincount(own[own >= 0], minlength=N).astype(float)
    Cn = np.bincount(z, minlength=N).astype(float)
    Sa, na, Sb, nb = _ancestry_scores(cl.members.astype(np.int64), cl.indptr.astype(np.int64),
                                      cl.mult.astype(np.float64), np.asarray(ek, np.int64), np.asarray(z, np.int64),
                                      N, T, Cn, int(z.size))
    parts = []
    for S, nn in ((Sa, na), (Sb, nb)):
        fit = np.flatnonzero((own >= 0) & (nn > 0))
        if fit.size > fit_max:
            fit = rng.choice(fit, fit_max, replace=False)
        lp = _apply_cal(S, nn, _calibrate(S[fit], nn[fit], own[fit]))
        lp[nn == 0] = -math.log(N)
        parts.append(lp)
    Z = np.mean(parts, 0)
    Z -= Z.mean(1, keepdims=True)
    ninf = np.maximum(na, nb)
    fit = np.flatnonzero((own >= 0) & (ninf > 0))
    if fit.size > fit_max:
        fit = rng.choice(fit, fit_max, replace=False)
    return _apply_cal(Z, np.ones(ek.size), _calibrate(Z[fit], np.ones(fit.size), own[fit])), ninf


def _clock_bins(a, se, y, m, edges):
    """Per K2P bin: median K2P, weighted mean share y, total weight, median s.e. (v5)."""
    bi = np.clip(np.searchsorted(edges, a, side="right") - 1, 0, edges.size - 2)
    nb = edges.size - 1
    W = np.bincount(bi, weights=m, minlength=nb)
    Y = np.bincount(bi, weights=m * y, minlength=nb) / np.maximum(W, 1e-300)
    A = np.array([np.median(a[bi == k]) if (bi == k).any() else np.nan for k in range(nb)])
    S = np.array([np.median(se[bi == k]) if (bi == k).any() else np.nan for k in range(nb)])
    ok = W > 0
    return A[ok], Y[ok], W[ok], S[ok]


def _one_step(A, Y, W, S, grid, young_high):
    """Weighted least-squares step between two levels blurred by K2P s.e. (v5)."""
    from scipy.stats import norm
    w = norm.cdf((grid[:, None] - A[None, :]) / np.maximum(S, 1e-5)[None, :])
    X1, X2 = w, 1 - w
    a11, a12, a22 = (W * X1 * X1).sum(1), (W * X1 * X2).sum(1), (W * X2 * X2).sum(1)
    b1, b2 = (W * X1 * Y).sum(1), (W * X2 * Y).sum(1)
    det = a11 * a22 - a12 ** 2
    det = np.where(np.abs(det) < 1e-12, np.nan, det)
    ly = np.clip((b1 * a22 - b2 * a12) / det, 0, 1)
    lo = np.clip((b2 * a11 - b1 * a12) / det, 0, 1)
    sse = (W[None, :] * (Y[None, :] - X1 * ly[:, None] - X2 * lo[:, None]) ** 2).sum(1)
    sse[~np.isfinite(sse) | ((ly < lo) if young_high else (ly > lo))] = np.inf
    if not np.isfinite(sse).any():
        return None
    j = int(np.argmin(sse))
    return float(grid[j]), float(ly[j]), float(lo[j])


def _divergence_step(A, Y, W, S, grid):
    k = np.convolve(W * Y, np.ones(3), "same") / np.maximum(np.convolve(W, np.ones(3), "same"), 1e-300)
    t = int(np.argmin(k))
    go = grid[grid >= A[t]]
    return _one_step(A[t:], Y[t:], W[t:], S[t:], go, False) if A.size - t >= 2 and go.size else None


def split_confinement(k2p, se, share, m, rng, n_bins=50):
    """v5 tentative split: beyond the confined trough, the share of each copy's lineage evidence pointing to the other
    subgenome rises again (lineages present in both progenitors). Step fit on LTR K2P bins; bootstrap over copies.
    Reported only when the step is supported (rise > 0 in >= 97.5% of resamples)."""
    ok = np.isfinite(k2p) & np.isfinite(se) & (m > 0)
    if ok.sum() < 200:
        return None
    a, s, y, w = k2p[ok], np.maximum(se[ok], 1e-5), share[ok], m[ok]
    hi = float(np.quantile(a, 0.995))
    q = np.unique(np.quantile(a[a > 0], np.linspace(0, 1, n_bins))) if (a > 0).any() else np.array([])
    edges = np.unique(np.concatenate([[0.0, 1e-12], q[(q > 1e-12) & (q < hi)], [hi, np.inf]]))
    grid = np.unique(np.concatenate([[0.0], np.linspace(0, hi, 400), edges[1:-1]]))
    fit = lambda idx: _divergence_step(*_clock_bins(a[idx], s[idx], y[idx], w[idx], edges), grid)
    dv = fit(np.arange(a.size))
    if dv is None:
        return None
    R = [fit(rng.integers(0, a.size, a.size)) for _ in range(BOOT)]
    t = np.array([r[0] for r in R if r is not None])
    rise = np.array([r[2] - r[1] for r in R if r is not None])
    if t.size < BOOT // 2 or np.quantile(rise, 0.025) <= 1e-9:              # rise > 0 beyond rounding
        return None
    lo, hi_ = widen(dv[0], *np.quantile(t, [0.025, 0.975]), "confinement")
    return dict(estimate=dv[0], lo=lo, hi=hi_, n=int(ok.sum()))


# ------------------------------------------------------------------------------------------------ dating: merger clock
def ltr_counts(lib):
    """Substitutions and aligned sites between each copy's two LTRs (annotation counts; else from K2P and length)."""
    n = np.where(np.isfinite(lib.n_sites) & (lib.n_sites > 0), lib.n_sites, np.nan_to_num(lib.ltr5_len, nan=0.0))
    k = lib.n_ts + lib.n_tv
    with np.errstate(invalid="ignore"):
        p = 0.75 * (1 - np.exp(-4 * lib.k2p / 3))                # NaN K2P stays NaN: no age, not age 0
    k = np.where(np.isfinite(k), k, np.round(p * n))
    return k, n


def _age_lik(k, n):
    """Binomial likelihood of k substitutions in n LTR sites at each grid divergence (rows: copies)."""
    from scipy.special import gammaln
    p = np.clip(0.75 * (1 - np.exp(-4 * AGE_GRID / 3)), 1e-12, 1 - 1e-12)
    k, n = np.asarray(k, float)[:, None], np.asarray(n, float)[:, None]
    ll = gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1) + k * np.log(p)[None] + (n - k) * np.log1p(-p)[None]
    return np.exp(ll - ll.max(1, keepdims=True))


def age_prior(k, n, df=8, c=1.0):
    """Smooth distribution of LTR divergence (insertion age) over AGE_GRID, by Efron g-modelling: g proportional to
    exp(Q beta), Q cubic B-splines on log(age + 1e-4), ridge penalty c |beta|^2, from the copies' counts."""
    from scipy.interpolate import BSpline
    from scipy.optimize import minimize
    L = _age_lik(k, n)
    x = np.log(AGE_GRID + 1e-4)
    t = np.linspace(x.min(), x.max(), df - 2)
    Q = BSpline.design_matrix(np.clip(x, t[0], t[-1]), np.concatenate(([t[0]] * 3, t, [t[-1]] * 3)), 3).toarray()
    Q -= Q.mean(0)

    def f(b):
        e = Q @ b
        g = np.exp(e - e.max())
        g /= g.sum()
        m = L @ g
        w = g * (L.T @ (1.0 / m))
        return -(np.log(m).sum() - c * b @ b), -(Q.T @ w - L.shape[0] * (Q.T @ g) - 2 * c * b)
    r = minimize(f, np.zeros(Q.shape[1]), jac=True, method="L-BFGS-B")
    e = Q @ r.x
    g = np.exp(e - e.max())
    return g / g.sum()


@njit(cache=True)
def _two_rates_ll(B, y, phi, eps):
    """Log-likelihood, gradient and Hessian of P(cross) = phi B + eps (1 - B) (concave in phi, eps)."""
    f = 0.0
    g1 = 0.0
    g2 = 0.0
    h11 = 0.0
    h12 = 0.0
    h22 = 0.0
    for i in range(B.size):
        b = B[i]
        p = min(max(phi * b + eps * (1 - b), 1e-12), 1 - 1e-12)
        f += y[i] * np.log(p) + (1 - y[i]) * np.log(1 - p)
        r = y[i] / p - (1 - y[i]) / (1 - p)
        c = y[i] / (p * p) + (1 - y[i]) / ((1 - p) * (1 - p))
        g1 += r * b
        g2 += r * (1 - b)
        h11 -= c * b * b
        h12 -= c * b * (1 - b)
        h22 -= c * (1 - b) * (1 - b)
    return f, g1, g2, h11, h12, h22


@njit(cache=True)
def _phi_eps_fit(B, y):
    """Maximum likelihood (phi, eps) in [1e-9, 1 - 1e-9]^2 with phi >= eps: projected Newton with an active set (a
    parameter at a bound whose gradient points outward stays there) and backtracking (every step raises the
    likelihood). If the unconstrained maximum has phi < eps, the constrained one has phi = eps = mean(y).
    Returns (log-likelihood, phi, eps, converged)."""
    lo = 1e-9
    hi = 1 - 1e-9
    m = 0.0
    for i in range(y.size):
        m += y[i]
    m = min(max(m / max(y.size, 1), lo), hi)
    phi = min(max(2 * m, 0.05), 0.95)
    eps = min(max(m / 2, lo), hi)
    f, g1, g2, h11, h12, h22 = _two_rates_ll(B, y, phi, eps)
    converged = False
    for _ in range(100):
        free1 = not ((phi <= lo and g1 < 0) or (phi >= hi and g1 > 0))
        free2 = not ((eps <= lo and g2 < 0) or (eps >= hi and g2 > 0))
        scale = 1.0 + abs(f)
        if (not free1 or abs(g1) < 1e-8 * scale) and (not free2 or abs(g2) < 1e-8 * scale):
            converged = True
            break
        d1 = 0.0
        d2 = 0.0
        det = h11 * h22 - h12 * h12
        if free1 and free2 and det > 1e-12 * (abs(h11 * h22) + 1e-300):
            d1 = -(h22 * g1 - h12 * g2) / det
            d2 = -(h11 * g2 - h12 * g1) / det
        elif free1 and h11 < 0:
            d1 = -g1 / h11
        elif free2 and h22 < 0:
            d2 = -g2 / h22
        else:                                                    # flat direction: gradient step
            d1 = g1 if free1 else 0.0
            d2 = g2 if free2 else 0.0
        step = 1.0
        moved = False
        while step > 1e-12:
            n1 = min(max(phi + step * d1, lo), hi)
            n2 = min(max(eps + step * d2, lo), hi)
            fn, a1, a2, b11, b12, b22 = _two_rates_ll(B, y, n1, n2)
            if fn >= f - 1e-12 * scale:
                moved = abs(n1 - phi) + abs(n2 - eps) > 1e-15
                phi, eps, f, g1, g2, h11, h12, h22 = n1, n2, fn, a1, a2, b11, b12, b22
                break
            step *= 0.5
        if not moved:
            converged = True
            break
    if phi < eps:                                                # one-sided: no young enrichment
        phi = eps = m
        f = _two_rates_ll(B, y, phi, eps)[0]
    return f, phi, eps, converged


@njit(parallel=True, cache=True)
def _clock_profile(B, y, cols):
    """Log-likelihood, phi, eps, converged at each candidate merger AGE_GRID[cols]; B[i, j] = P(pair i's younger copy
    was inserted less than AGE_GRID[j] ago)."""
    out = np.empty((cols.size, 4))
    for c in prange(cols.size):
        ll, phi, eps, ok = _phi_eps_fit(np.ascontiguousarray(B[:, cols[c]]), y)
        out[c, 0] = ll
        out[c, 1] = phi
        out[c, 2] = eps
        out[c, 3] = 1.0 if ok else 0.0
    return out


def _landing_bound(B, y, cols, phi, clusters):
    """Upper bound on the merger (LTR divergence units) when the clock sees no excess of moves, IF new copies land in
    either subgenome in proportion to its size: phi_i = the other subgenome's genome share seen from pair i's older
    copy. A merger at T predicts P(split) = phi_i P(younger < T) + eps (1 - P(younger < T)), eps = the observed split
    share; per lineage-cluster bootstrap draw (own generator: the dating's draws stay untouched) the oldest T within
    1.92 log-likelihood units of the best (T = 0 included); bound = 97.5% quantile widened by the merger model error
    (merger_clock floors it at the clock's resolution); None when even the oldest candidate survives in > 2.5% of
    draws. Fails when transposition stays within each subgenome after the merger (stated wherever it is reported)."""
    eps = min(max(float(y.mean()), 1e-6), 1 - 1e-6)
    Bm = B[:, cols]
    P = np.clip(phi[:, None] * Bm + eps * (1 - Bm), 1e-12, 1 - 1e-12)
    LL = y[:, None] * np.log(P) + (1 - y[:, None]) * np.log1p(-P)
    l0 = y * math.log(eps) + (1 - y) * math.log1p(-eps)
    rng = np.random.default_rng([int(y.size), int(y.sum())])
    W = np.array([np.bincount(sel, minlength=y.size) for sel in _cluster_boot(clusters, rng, BOOT)], float)
    llT, ll0 = W @ LL, W @ l0
    ok = llT >= (np.maximum(llT.max(1), ll0) - 1.92)[:, None]
    last = np.array([cols[np.flatnonzero(o)[-1]] if o.any() else cols[0] for o in ok])
    if (last == cols[-1]).mean() > 0.025:
        return None
    return float(np.quantile(AGE_GRID[last], 0.975) * math.exp(1.96 * MODEL_ERR["merger"]))


def age_posteriors(k, n, g):
    """Each copy's insertion-age posterior over AGE_GRID (rows) from its LTR counts and the g-model prior g."""
    post = _age_lik(k, n) * g[None]
    return post / np.maximum(post.sum(1, keepdims=True), 1e-300)


def _div_lik(k, n, scale):
    """Binomial likelihood of k substitutions in n sites of a region evolving at 1/scale of the LTR rate, at each LTR
    divergence of AGE_GRID (rows max-normalised; n = 0 rows are flat)."""
    from scipy.special import gammaln
    p = np.clip(0.75 * (1 - np.exp(-4 * AGE_GRID / (3 * scale))), 1e-12, 1 - 1e-12)
    k, n = np.asarray(k, float)[:, None], np.asarray(n, float)[:, None]
    ll = gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1) + k * np.log(p)[None] + (n - k) * np.log1p(-p)[None]
    return np.exp(ll - ll.max(1, keepdims=True))


def younger_cdf(k_u, n_u, k_v, n_v, g, kD=None, nD=None, scale=1.0):
    """B[i, j] = P(the younger copy of pair i was inserted less than AGE_GRID[j] ago).
    LTR ages only (v7): 1 - S_u S_v. With the pair's own substitutions (kD of nD sites, v8): the younger copy cannot
    be older than the pair's divergence (both copies' changes since the younger insertion are in it), so the joint
    age posterior is weighted by max over delta >= min(a_u, a_v) of Binomial(kD | nD, delta / scale), scale = LTR
    rate / rate of the aligned region (measured; LTRs diverge faster than whole elements): no assumption about how
    long before the younger copy the two lineages separated (pairs without counts, nD = 0, are unconstrained)."""
    Fu, Fv = age_posteriors(k_u, n_u, g), age_posteriors(k_v, n_v, g)
    Su = np.cumsum(Fu[:, ::-1], 1)[:, ::-1]                      # P(age >= AGE_GRID[j])
    Sv = np.cumsum(Fv[:, ::-1], 1)[:, ::-1]
    if kD is None:
        return np.ascontiguousarray(1.0 - Su * Sv)
    h = (Fu * Sv + Fv * Su - Fu * Fv) * np.maximum.accumulate(_div_lik(kD, nD, scale)[:, ::-1], 1)[:, ::-1]
    h /= np.maximum(h.sum(1, keepdims=True), 1e-300)             # P(min age = AGE_GRID[j] | LTRs, pair)
    return np.ascontiguousarray(np.clip(np.cumsum(h, 1) - h, 0.0, 1.0))


R_GRID = np.linspace(0.025, 1.0, 40)      # v8: the younger copy's age as a fraction of its pair's divergence
G_GRID = np.concatenate(([0.0], np.geomspace(1e-6, 0.3, 79)))   # v8: pair divergence grid (LTR units)


def _binom_ll(k, n, p):
    from scipy.special import gammaln
    return gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1) + k * np.log(p) + (n - k) * np.log1p(-p)


def _p_of(a):
    """Expected substitution share at divergence a (Jukes-Cantor), clipped away from 0 and 1."""
    return np.clip(0.75 * (1 - np.exp(-4 * np.asarray(a, float) / 3)), 1e-12, 1 - 1e-12)


def sibling_genealogy(k_u, n_u, k_v, n_v, delta, iters=200):
    """f(r) over R_GRID: how far back, as a fraction r of the pair's divergence delta (LTR units), the younger copy of a
    sibling pair was inserted. Learned by EM from sibling pairs whose LTRs carry several substitutions; each pair's
    likelihood is L(r) = 1/2 sum over which copy is the younger of Binom(k_young | r delta) x max over d >= delta of
    Binom(k_older | d). r = 1: copied from its partner; r uniform: both drawn from one source lineage at random times.
    The genome's own genealogy, not an assumed one, then maps a young pair's divergence to its younger copy's age."""
    dgrid = np.geomspace(1e-5, 0.6, 120)
    lik = lambda k, n, a: np.exp(_binom_ll(k[:, None], n[:, None], _p_of(a)))
    Ly_u, Ly_v = lik(k_u, n_u, R_GRID[None] * delta[:, None]), lik(k_v, n_v, R_GRID[None] * delta[:, None])
    above = dgrid[None, :] >= delta[:, None]
    Lo_u = np.where(above, lik(k_u, n_u, dgrid[None]), 0.0).max(1)
    Lo_v = np.where(above, lik(k_v, n_v, dgrid[None]), 0.0).max(1)
    L = 0.5 * (Ly_u * Lo_v[:, None] + Ly_v * Lo_u[:, None]) + 1e-300
    f = np.full(R_GRID.size, 1.0 / R_GRID.size)
    for _ in range(iters):
        post = L * f[None]
        post /= post.sum(1, keepdims=True)
        f = post.mean(0)
    return f


def younger_cdf_genealogy(k_u, n_u, k_v, n_v, kD, nD, scale, f, chunk=4000):
    """B[i, j] = P(the younger copy of pair i was inserted less than AGE_GRID[j] ago): the pair's divergence delta (LTR
    units; its kD of nD whole-record sites at delta / scale, uniform prior over G_GRID), the younger copy's age r delta
    with r from the learned genealogy f updated by both copies' LTR counts (either may be the younger)."""
    n = k_u.size
    B = np.zeros((n, AGE_GRID.size))
    w0 = np.diff(np.concatenate(([0.0], G_GRID)))
    A = R_GRID[None, :] * G_GRID[:, None]                         # younger ages (G, R)
    col = np.searchsorted(AGE_GRID, A, side="right").ravel()       # first grid age above each younger age
    for s in range(0, n, chunk):
        e = min(n, s + chunk)
        m = e - s
        Pd = np.exp(_binom_ll(kD[s:e, None], nD[s:e, None], _p_of(G_GRID[None] / scale))) * w0[None]
        Pd /= np.maximum(Pd.sum(1, keepdims=True), 1e-300)
        lk = lambda k, nn: np.exp(_binom_ll(k[:, None, None], nn[:, None, None], _p_of(A[None])))
        Ou = np.exp(_binom_ll(k_u[s:e, None], n_u[s:e, None], _p_of(G_GRID[None])))[:, :, None]
        Ov = np.exp(_binom_ll(k_v[s:e, None], n_v[s:e, None], _p_of(G_GRID[None])))[:, :, None]
        Pr = f[None, None, :] * 0.5 * (lk(k_u[s:e], n_u[s:e]) * Ov + lk(k_v[s:e], n_v[s:e]) * Ou)
        Pr /= np.maximum(Pr.sum(2, keepdims=True), 1e-300)          # P(r | delta, LTRs)
        Wt = (Pd[:, :, None] * Pr).reshape(m, -1)
        H = np.zeros((m, AGE_GRID.size + 1))
        np.add.at(H, (np.repeat(np.arange(m), col.size), np.tile(col, m)), Wt.ravel())
        B[s:e] = np.cumsum(H[:, :-1], 1)
    return np.ascontiguousarray(np.clip(B, 0.0, 1.0))


def merger_clock(B, cross, clusters, rng, t_max, t_res, phi_land=None):
    """Pairwise sibling clock: a sibling pair sits in two subgenomes only when its younger copy was inserted after the
    merger (eps: exchange or misassignment noise, age-independent). T (LTR divergence units) by profile likelihood
    over AGE_GRID below t_max, given B = P(younger copy inserted < T) per pair (younger_cdf); lineage-cluster bootstrap
    widened by MODEL_ERR["merger"]. Not datable unless cross pairs are >= 2-fold enriched among post-merger pairs
    (phi >= 2 eps) in >= 97.5% of resamples; an upper bound only ("younger than") when the interval reaches the
    clock's resolution t_res (one substitution per LTR for LTR ages, per aligned pair with the pair divergence)."""
    y = cross.astype(np.float64)
    cols = np.flatnonzero((AGE_GRID > 0) & (AGE_GRID < t_max)).astype(np.int64)
    if cols.size < 3 or y.size < 50:
        return None
    prof = _clock_profile(B, y, cols)
    i = int(np.argmax(prof[:, 0]))
    T, phi, eps = AGE_GRID[cols[i]], prof[i, 1], prof[i, 2]
    Tb, gap, bad = [], [], int((prof[:, 3] == 0).sum())
    for sel in _cluster_boot(clusters, rng, BOOT):
        pb = _clock_profile(np.ascontiguousarray(B[sel]), y[sel], cols)
        bad += int((pb[:, 3] == 0).sum())
        j = int(np.argmax(pb[:, 0]))
        Tb.append(AGE_GRID[cols[j]])
        gap.append(pb[j, 1] - 2 * pb[j, 2])
    lo, hi = widen(T, *np.quantile(Tb, [0.025, 0.975]), "merger")
    bound = bool(np.quantile(Tb, 0.025) <= t_res)              # resolution judged on the bootstrap alone
    m_ = min(max(float(y.mean()), 1e-12), 1 - 1e-12)
    log_lr = float(prof[i, 0] - (y * math.log(m_) + (1 - y) * math.log1p(-m_)).sum())   # step vs no step
    datable = bool(phi >= 2 * eps and np.quantile(gap, 0.025) > 0 and log_lr >= math.log(BF_DECISIVE))
    if bad:
        log(f"WARNING: the merger clock's rate fit did not converge in {bad} of {cols.size * (BOOT + 1)} fits")
    lb = _landing_bound(B, y, cols, np.asarray(phi_land, float), clusters) if (
        phi_land is not None and not datable and float(np.mean(phi_land)) >= 2 * float(y.mean())) else None
    if lb is not None:                                     # (random landing must predict >= 2x the background,
        lb = max(lb, t_res)                                #  the clock's rule; resolution floor: one substitution
        lb = lb if lb < t_max else None                    #  per LTR; a bound at or above the split says nothing)
    return dict(landing_bound=lb, estimate=float(T), lo=0.0 if bound else float(lo), hi=float(hi), phi=float(phi),
                eps=float(eps), datable=datable, bound=bound, n=int(y.size), n_cross=int(y.sum()), t_res=float(t_res),
                log_lr=log_lr)


# ------------------------------------------------------------------------------------------------ dating: pipeline
SPLIT_ORDER = {"orthologs": ("orthologs", "confinement", "profile"), "confinement": ("confinement", "orthologs", "profile"),
               "profile": ("profile", "confinement", "orthologs"), "agree": ("profile", "confinement", "orthologs")}


def date_pair(lib, pt, loc, x, y, ancestry, rng, split_rule, sizes=None, ptm=None, units=None, rng_m=None, threads=1,
              locm=None):
    """Split and merger of subgenomes x and y (group indices) from the copy-pair table pt and copy locations loc.
    K2P units (pair divergence = 2 mu T; LTR divergence of a copy = 2 mu age). v8: the merger uses the extended table
    ptm (v7's pairs + extra young focal copies; see merger_pair); the split estimators use v7's table pt."""
    P, d, se = pt["P"], pt["d"], pt["se"]
    la, lb = loc[P[:, 0]], loc[P[:, 1]]
    m0 = np.isin(la, [x, y]) & np.isin(lb, [x, y])
    out = dict(groups=(x, y), split={}, merger=None, n_pairs=int(m0.sum()))
    if m0.sum() < 200:
        return out
    Pm, dm, sm = P[m0], d[m0], se[m0]
    cross = (la[m0] != lb[m0]).astype(np.float64)
    clusters = lineage_clusters(Pm, dm, lib.n, CLUSTER_D)[Pm[:, 0]]
    dup = duplicated_insertions(lib, Pm, dm, sm)
    orth = (cross > 0) & dup & (dm > young_cut(dm, cross))
    out["split"]["orthologs"] = split_orthologs(dm, orth, clusters, rng)
    out["split"]["profile"] = split_profile(dm, cross, clusters, rng)
    lP, ninf = ancestry
    inside = np.isin(loc, [x, y]) & (ninf > 0)
    px, py = np.exp(lP[:, x]), np.exp(lP[:, y])
    share = np.where(loc == x, py, px) / np.maximum(px + py, 1e-300)
    out["split"]["confinement"] = split_confinement(lib.k2p, k2p_se_copy(lib), share, inside.astype(float), rng)
    cap = {k: v["estimate"] for k, v in out["split"].items() if v}       # pair-divergence units (merger clock cap)
    if SPLIT_LTR_UNITS:                          # orthologs/profile: whole-record pair K2P -> LTR K2P (confinement's)
        cr = cross > 0
        for k in ("orthologs", "profile"):
            v = out["split"][k]
            sc = ltr_scale_near(lib, Pm[cr], dm[cr], pt["cnt"][m0][cr], v["estimate"], rng_m or rng, threads) if v \
                else None
            if v and sc:
                out["split"][k] = dict(v, estimate=v["estimate"] * sc, lo=v["lo"] * sc, hi=v["hi"] * sc, ltr_scale=sc)
    rep = next((k for k in SPLIT_ORDER[split_rule] if out["split"][k]), None)
    o_, p_ = out["split"]["orthologs"], out["split"]["profile"]
    if split_rule == "agree" and o_ and p_ and o_["lo"] <= p_["hi"] and p_["lo"] <= o_["hi"]:
        rep = "orthologs"                         # inherited insertions, confirmed by the lineage-floor profile
    out["split_reported"] = rep
    ests = [v for v in out["split"].values() if v and np.isfinite(v["lo"])]
    out["split_agreement"] = ("none" if not ests else "single" if len(ests) == 1 else
                              "agree" if max(v["lo"] for v in ests) <= min(v["hi"] for v in ests) else "disagree")
    cap7 = next((cap[k] for k in SPLIT_ORDER["confinement"] if k in cap), None)     # v7's clock cap (as evaluated)
    out.update(merger_pair(lib, ptm or pt, loc if locm is None else locm, x, y, rng if rng_m is None else rng_m,
                           cap7, sizes, units, threads))
    return out


def merger_pair(lib, pt, loc, x, y, rng, split, sizes, units, threads=1):
    """Merger of subgenomes x and y from sibling copy pairs (MRCA = the older copy's insertion; no duplicated
    insertions) younger than the confined plateau's end and the reported split. Two clocks on the same
    pairs and bootstrap units: the LTR clock, the younger copy's age from its LTR counts (v7); the genealogy clock, the
    pair's divergence over the whole records (one substitution per aligned pair instead of per LTR) read through the
    genome's sibling genealogy f(r) learned from pairs with informative LTRs (sibling_genealogy). Reported:
    merger_report."""
    P, d, se = pt["P"], pt["d"], pt["se"]
    la, lb = loc[P[:, 0]], loc[P[:, 1]]
    m0 = np.isin(la, [x, y]) & np.isin(lb, [x, y])
    out = dict(merger=None, merger_ltr=None, merger_genealogy=None, n_sibling_pairs=0)
    if m0.sum() < 200:
        return out
    Pm, dm, sm, cm = P[m0], d[m0], se[m0], pt["cnt"][m0]
    cross = (la[m0] != lb[m0]).astype(np.float64)
    if units is not None:
        comp, in7 = units
        clusters = np.where(in7[Pm[:, 0]], comp[Pm[:, 0]], comp[Pm[:, 1]])
    else:
        clusters = lineage_clusters(Pm, dm, lib.n, CLUSTER_D)[Pm[:, 0]]
    dup = duplicated_insertions(lib, Pm, dm, sm)
    ks = k2p_se_copy(lib)
    ka, kb = lib.k2p[Pm[:, 0]], lib.k2p[Pm[:, 1]]
    old = np.where(ka >= kb, Pm[:, 0], Pm[:, 1])
    with np.errstate(invalid="ignore"):
        sib = dm <= np.fmax(ka, kb) + 2 * np.sqrt(ks[old] ** 2 + sm ** 2)
    t_max = plateau_end(dm, cross)
    if split:
        t_max = min(t_max, split)
    sel = ~dup & sib & (dm < t_max)
    kc, nc = ltr_counts(lib)
    good = sel & np.isfinite(kc[Pm[:, 0]]) & np.isfinite(kc[Pm[:, 1]]) & (nc[Pm[:, 0]] > 0) & (nc[Pm[:, 1]] > 0)
    out["n_sibling_pairs"] = int(good.sum())
    pool = np.flatnonzero(np.isin(loc, [x, y]) & np.isfinite(kc) & (nc > 0))
    if good.sum() < 50 or pool.size < 50:
        return out
    pool = rng.choice(pool, min(pool.size, 20000), replace=False)
    g = age_prior(kc[pool], nc[pool])
    u, v = Pm[good, 0], Pm[good, 1]
    phi_land = None
    if sizes is not None:                              # random landing: the other subgenome's genome share
        og = loc[old[good]]
        phi_land = np.where(og == x, sizes[y], sizes[x]) / (sizes[x] + sizes[y])
    y_ = cross[good] > 0
    mg = merger_clock(younger_cdf(kc[u], nc[u], kc[v], nc[v], g), y_, clusters[good], rng, t_max,
                      1.0 / np.median(np.concatenate([nc[u], nc[v]])), phi_land=phi_land)
    kD, nD = (cm[good, 1] + cm[good, 2]).astype(float), cm[good, 0].astype(float)
    scale = ltr_rate_ratio(lib, Pm[good], cm[good], rng, threads)
    p_ = np.clip(kD / np.maximum(nD, 1.0), 0, 0.7)
    dlt = -0.75 * np.log1p(-4 * p_ / 3) * scale                  # pair divergence, LTR units
    teach = (dlt >= GEN_MIN_SUBS / np.median(np.minimum(nc[u], nc[v]))) & (dlt < t_max)
    f = (sibling_genealogy(kc[u][teach], nc[u][teach], kc[v][teach], nc[v][teach], dlt[teach])
         if teach.sum() >= GEN_MIN_PAIRS else np.full(R_GRID.size, 1.0 / R_GRID.size))
    mp = merger_clock(younger_cdf_genealogy(kc[u], nc[u], kc[v], nc[v], kD, nD, scale, f), y_, clusters[good], rng,
                      t_max, scale / np.median(nD))
    out.update(ltr_rate_ratio=scale, genealogy=f, n_teach=int(teach.sum()))
    out["merger_ltr"], out["merger_genealogy"] = mg, mp
    out["merger"] = merger_report(mg, mp, 1.0 / float(np.median(nc[np.isfinite(nc) & (nc > 0)])))
    return out


def merger_report(mg, mp, t_ltr):
    """The reported merger of a subgenome pair, by the order fixed on simulations (dev_v8): the genealogy clock's
    two-sided interval; else the LTR clock's; else the tighter of their bounds (the genealogy clock's upper limit is
    usually far tighter: one substitution per aligned pair instead of per LTR); else 'not datable' (the LTR clock's
    record, with its supplementary random-landing bound). A two-sided interval whose lower limit lies below t_ltr (one
    substitution per LTR) is reported as a bound [0, upper limit]: below the LTRs' own resolution a lower limit rests on
    the learned genealogy alone, and simulations put those too old (TE-poor genomes, hexaploids, lineages shared before
    the merger). Each dict gains 'source' (ltr_clock / genealogy_clock)."""
    if mp and mp["datable"] and not mp["bound"]:
        rep = dict(mp, source="genealogy_clock", landing_bound=None)
    elif mg and mg["datable"] and not mg["bound"]:
        rep = dict(mg, source="ltr_clock")
    else:
        bounds = [dict(m, source=k) for k, m in (("genealogy_clock", mp), ("ltr_clock", mg)) if m and m["datable"]]
        if not bounds:
            return dict(mg, source="ltr_clock") if mg else None
        rep = min(bounds, key=lambda m: m["hi"])
        if rep["source"] == "genealogy_clock":
            rep["landing_bound"] = None
    if not rep["bound"] and rep["lo"] < t_ltr:
        rep = dict(rep, bound=True, lo=0.0)
    return rep


def split_text(r, mu):
    """Log text of a subgenome pair's split (years)."""
    sp = r["split"].get(r.get("split_reported")) if r.get("split_reported") else None
    if not sp:
        return "split not datable"
    yr = lambda k: k / (2 * mu) / 1e6
    note = {"single": "one estimator", "agree": "estimators agree", "disagree": "estimators disagree"}
    return (f"split {yr(sp['estimate']):.2f} Myr ({yr(sp['lo']):.2f}-{yr(sp['hi']):.2f}; {r['split_reported']}, "
            f"{note.get(r.get('split_agreement'), '')})")


def merger_text(mg, n_sib, mu):
    """Log text of the merger clock: a date, a bound on the onset of cross-subgenome transposition, or why not."""
    kyr = lambda k: f"{k / (2 * mu) / 1e3:,.1f} kyr"
    if mg is None:
        return f"merger not datable: too few sibling copy pairs ({n_sib})"
    if not mg["datable"]:
        return ("merger not datable: no detectable excess of cross-subgenome sibling pairs (no copies moved between "
                "these subgenomes since the merger, or too few to tell)"
                + (f"; if new copies land in either subgenome at random, the merger is younger than "
                   f"{kyr(mg['landing_bound'])} (wrong where new copies stay in their own subgenome)"
                   if mg.get("landing_bound") else ""))
    src = {"ltr_clock": "LTR clock", "genealogy_clock": "genealogy clock"}.get(mg.get("source"), "clock")
    if mg["bound"]:
        return (f"cross-subgenome transposition began < {kyr(mg['hi'])} ago ({src}; a bound on its onset, not on the "
                f"merger)")
    return (f"merger {kyr(mg['estimate'])} ({kyr(mg['lo'])}-{kyr(mg['hi'])}; {src}, onset of cross-subgenome "
            f"transposition)")


def run_dating(a, lib, ph):
    """Dating of every pair of subgenomes (N >= 2): copy pairs, three split estimators, the merger clock."""
    if ph["N"] < 2:
        return None
    rng = np.random.default_rng([a.seed, 7])
    ek = ph["ek"]
    si = {s_: i for i, s_ in enumerate(lib.seqnames)}
    sizes = np.bincount(ph["z"], weights=[float(lib.seq_len[si[s_]]) for s_ in ph["names"]], minlength=ph["N"])
    loc = np.where(ek >= 0, ph["z"][np.maximum(ek, 0)], -1)
    pt = sample_pairs(lib, ek, rng, a.threads)
    if pt is None or len(pt["P"]) < 200:
        log("dating skipped: too few copy pairs")
        return None
    vlog(f"dating: {len(pt['P']):,} copy pairs aligned")
    ancestry = copy_ancestry(ph["clades"], ek, ph["z"], ph["N"], rng)
    rng_m = np.random.default_rng([a.seed, 80])               # v8's merger draws: v7's stream stays as it was
    ptm = sample_young_more(lib, ek, pt, rng_m, a.threads)
    pt.pop("sketch", None)
    vlog(f"dating: {len(ptm['P']) - len(pt['P']):,} more copy pairs for the merger clock (young focal copies)")
    units = units_v7(pt, lib.n)
    pa = run_painting(a, lib, ph, dict(ancestry=ancestry))     # v8: painted first (same ancestry: same result as v7's)
    locm = exchange_locations(lib, loc, pa["segs"]) if (MERGER_PAINT and pa) else loc
    res = []
    for x in range(ph["N"]):
        for y in range(x + 1, ph["N"]):
            res.append(date_pair(lib, pt, loc, x, y, ancestry, rng, SPLIT_RULE, sizes=sizes, ptm=ptm, units=units,
                                 rng_m=rng_m, threads=a.threads, locm=locm))
    mu = a.mu or MU_DEFAULT
    for r in res:
        x, y = (group_label(g) for g in r["groups"])
        log(f"{x}-{y}: {split_text(r, mu)}; {merger_text(r['merger'], r.get('n_sibling_pairs', 0), mu)} "
            f"[mu={mu:g}{'' if a.mu else ', default'}]")
    return dict(mu=mu, mu_given=bool(a.mu), pairs=res, n_copy_pairs=int(len(pt["P"])), ancestry=ancestry, pt=pt,
                n_copy_pairs_merger=int(len(ptm["P"])), painting=pa)


def dating_summary(dt):
    """Compact dating record for summary.json (years; K2P in dating.tsv)."""
    yr = lambda k: None if k is None else k / (2 * dt["mu"])
    pairs = []
    for r in dt["pairs"]:
        sp = r["split"].get(r.get("split_reported")) if r.get("split_reported") else None
        mg = r["merger"]
        st = "none" if mg is None else ("not_datable" if not mg["datable"] else ("bound" if mg["bound"] else "dated"))
        pairs.append(dict(subgenomes=[group_label(g) for g in r["groups"]],
                          split=dict(estimator=r.get("split_reported"), years=yr(sp["estimate"]),
                                     low=yr(sp["lo"]), high=yr(sp["hi"]), agreement=r.get("split_agreement"))
                          if sp else None,
                          merger=dict(status=st, years=yr(mg["estimate"]) if st == "dated" else None,
                                      low=yr(mg["lo"]) if st == "dated" else (0.0 if st == "bound" else None),
                                      high=yr(mg["hi"]) if st in ("dated", "bound") else None,
                                      clock=mg.get("source") if mg else None,
                                      if_random_landing_below_years=yr(mg["landing_bound"]) if mg and
                                      mg.get("landing_bound") else None),
                          ltr_rate_ratio=r.get("ltr_rate_ratio"),
                          genealogy_mean_r=float(np.dot(R_GRID, r["genealogy"])) if r.get("genealogy") is not None
                          else None, n_sibling_pairs=r.get("n_sibling_pairs")))
    return dict(mu=dt["mu"], mu_given=dt["mu_given"], n_copy_pairs=dt["n_copy_pairs"],
                n_copy_pairs_merger=dt.get("n_copy_pairs_merger"), pairs=pairs)


def write_dating(a, dt):
    od = Path(a.outdir)
    with open(od / "dating.tsv", "w") as f:
        f.write(f"# mu={dt['mu']:g} substitutions/site/year{'' if dt['mu_given'] else ' (default)'}; K2P = 2 mu years\n")
        f.write("subgenome_a\tsubgenome_b\tevent\testimator\treported\tstatus\tk2p\tk2p_low\tk2p_high\tyears\t"
                "years_low\tyears_high\tn\n")
        for r in dt["pairs"]:
            x, y = (group_label(g) for g in r["groups"])
            for k in ("orthologs", "confinement", "profile"):
                v = r["split"].get(k)
                st = "dated" if v else "none"
                row = [x, y, "split", k, "yes" if r.get("split_reported") == k else "no", st]
                vals = [v["estimate"], v["lo"], v["hi"]] if v else [None] * 3
                f.write("\t".join(row + [_f(z_) for z_ in vals] + [_f(z_ / (2 * dt["mu"]) if z_ is not None else None)
                                                                     for z_ in vals] + [str(v["n"]) if v else "0"]) + "\n")
            mg = r["merger"]
            src = mg.get("source") if mg else None
            for k, m_ in (("ltr_clock", r.get("merger_ltr")), ("genealogy_clock", r.get("merger_genealogy"))):
                if k == src:
                    m_ = mg                      # the reported row carries what is reported (a demoted interval: bound)
                st = "none" if m_ is None else ("not_datable" if not m_["datable"] else
                                                ("bound" if m_["bound"] else "dated"))
                vals = ([None if m_["bound"] else m_["estimate"], m_["lo"], m_["hi"]] if m_ and m_["datable"]
                        else [None] * 3)                   # a bound has no point estimate
                f.write("\t".join([x, y, "merger", k, "yes" if src == k else "no", st] + [_f(z_) for z_ in vals] +
                                   [_f(z_ / (2 * dt["mu"]) if z_ is not None else None) for z_ in vals] +
                                   [str(m_["n"]) if m_ else "0"]) + "\n")
            mg = r.get("merger_ltr")
            lb = mg.get("landing_bound") if mg else None             # supplementary: assumes random landing
            vals = [None, 0.0, lb] if lb else [None] * 3
            f.write("\t".join([x, y, "merger_if_random", "landing", "no", "bound" if lb else "none"]
                               + [_f(z_) for z_ in vals] + [_f(z_ / (2 * dt["mu"]) if z_ is not None else None)
                                                            for z_ in vals] + [str(mg["n"]) if mg else "0"]) + "\n")


# ------------------------------------------------------------------------------------------------ derived biology
def clade_names(lib):
    """Annotated LTR-RT clade of each copy (Tekay, CRM, Athila, ...): the GFF clade, else the third field of the FASTA
    header class (LTR/Gypsy/Tekay); 'unknown' when neither names one."""
    out = []
    for c, k in zip(lib.clade, getattr(lib, "cls", [""] * len(lib.clade))):
        if not c:
            f = k.split("/") if k else []
            c = f[2] if len(f) > 2 and f[2] else ""
        out.append(c if c else "unknown")
    return np.array(out)


def merger_cut(r):
    """LTR K2P separating copies inserted before and after the merger of a subgenome pair: the merger estimate; None
    when the merger is only bounded (a bound limits the onset of cross-subgenome transposition, not the merger) or not
    datable."""
    mg = r["merger"]
    if mg is None or not mg["datable"] or mg["bound"]:
        return None
    return float(mg["estimate"])


def after_merger(lib, loc, dt, rng):
    """Per subgenome with a dated merger (the earliest one when several pairs include it): each copy's probability of
    having been inserted after it, from its LTR counts and the subgenome's g-model age distribution (a hard K2P cut
    would count pre-merger copies whose short LTRs carry no substitution) -> {subgenome: array over copies, NaN
    elsewhere}."""
    cut = {}
    for r in (dt["pairs"] if dt else []):
        T = merger_cut(r)
        if T is not None:
            for g in r["groups"]:
                cut[g] = max(cut.get(g, 0.0), T)
    kc, nc = ltr_counts(lib)
    out = {}
    for g, T in sorted(cut.items()):
        idx = np.flatnonzero((loc == g) & np.isfinite(kc) & (nc > 0))
        if idx.size < MIN_HIST:
            continue
        pool = rng.choice(idx, min(idx.size, 20000), replace=False)
        post = _age_lik(kc[idx], nc[idx]) * age_prior(kc[pool], nc[pool])[None]
        p = np.full(lib.n, np.nan)
        p[idx] = post[:, AGE_GRID < T].sum(1) / np.maximum(post.sum(1), 1e-300)
        out[g] = p
    return out


def nesting_check(lib):
    """Elements inside another element (host = the smallest element strictly containing it on the same sequence): the
    inner one inserted later, so its LTR K2P should not exceed its host's. Violation: inner - host > 2 se (each from
    its own LTR sites; pairs without LTR site counts or lengths are not checked) -> list of dicts (inner, host, K2Ps,
    se, violation, inside a host LTR)."""
    se = k2p_se_copy(lib)
    has_se = (np.isfinite(lib.k2p_se) & (lib.k2p_se > 0)) | (np.nan_to_num(lib.n_sites) > 0) | \
        (np.nan_to_num(lib.ltr5_len) > 0)                    # K2P alone carries no standard error
    out = []
    order = np.lexsort((-lib.end, lib.start, lib.seq))
    stack = []
    for i in order:
        while stack and (lib.seq[stack[-1]] != lib.seq[i] or lib.end[stack[-1]] < lib.end[i]):
            stack.pop()
        if stack and lib.start[stack[-1]] <= lib.start[i] and (lib.start[stack[-1]], lib.end[stack[-1]]) != (lib.start[i], lib.end[i]):
            h = stack[-1]
            if np.isfinite(lib.k2p[i]) and np.isfinite(lib.k2p[h]) and has_se[i] and has_se[h]:
                s = float(np.hypot(se[i], se[h]))
                l5 = lib.ltr5_len[h] if np.isfinite(lib.ltr5_len[h]) else 0
                l3 = lib.ltr3_len[h] if np.isfinite(lib.ltr3_len[h]) else 0
                out.append(dict(inner=int(i), host=int(h), k_inner=float(lib.k2p[i]), k_host=float(lib.k2p[h]), se=s,
                                violation=bool(lib.k2p[i] - lib.k2p[h] > 2 * s),
                                in_host_ltr=bool(lib.start[i] < lib.start[h] + l5 or lib.end[i] > lib.end[h] - l3)))
        stack.append(i)
    return out


def hist_mass(g):
    """Mass of the age distribution g (point masses on AGE_GRID) in each HIST_EDGES bin: each positive grid point's
    mass spread over its log-age cell, read at the bin edges, so a bin's mass follows its width rather than how many
    grid points fall in it; the age-0 point goes to the first bin, ages above the last edge to the last."""
    a = AGE_GRID[1:]
    r = np.sqrt(a[1] / a[0])
    x = np.log(np.concatenate(([a[0] / r], a * r)))
    C = np.concatenate(([g[0]], g[0] + np.cumsum(g[1:])))
    Ce = np.interp(np.log(HIST_EDGES[1:]), x, C)
    Ce[-1] = C[-1]
    return np.diff(np.concatenate(([0.0], Ce)))


def insertion_history(lib, loc, N, rng):
    """Copies per LTR K2P bin (HIST_EDGES) per subgenome with >= MIN_HIST copies, all clades and each clade with
    >= MIN_HIST copies there (not one holding the whole subgenome): observed, and deconvolved (g-model age distribution
    of age_prior scaled to the copies; removes the binomial noise of short LTRs, not survival or detection biases; NaN
    without LTR substitution counts or lengths)."""
    kc, nc = ltr_counts(lib)
    fin = np.isfinite(lib.k2p)
    cnt = fin & np.isfinite(kc) & (nc > 0)
    clades = clade_names(lib)
    rows = []
    for g in range(N):
        mg = fin & (loc == g)
        cand = ["all"] + sorted(c for c in set(clades[mg]) if MIN_HIST <= (mg & (clades == c)).sum() < mg.sum())
        for c in cand:
            m = mg if c == "all" else mg & (clades == c)
            if m.sum() < MIN_HIST:
                continue
            idx = np.flatnonzero(m & cnt)
            dens = np.full(HIST_EDGES.size - 1, np.nan)
            if idx.size >= MIN_HIST:
                pool = rng.choice(idx, min(idx.size, 20000), replace=False)
                dens = hist_mass(age_prior(kc[pool], nc[pool])) * idx.size
            obs = np.bincount(np.clip(np.searchsorted(HIST_EDGES, lib.k2p[m], "right") - 1, 0, HIST_EDGES.size - 2),
                              minlength=HIST_EDGES.size - 1)
            for b in range(HIST_EDGES.size - 1):
                rows.append(dict(subgenome=g, clade=c, lo=HIST_EDGES[b], hi=HIST_EDGES[b + 1], observed=int(obs[b]),
                                 deconvolved=float(dens[b])))
    return rows


def clade_table(lib, loc, N, after):
    """Per clade and subgenome: copies, kb, LTR K2P quartiles, and the expected number of copies inserted after the
    subgenome's dated merger (sum of after_merger's probabilities; None without a dated merger)."""
    clades = clade_names(lib)
    L = (lib.end - lib.start + 1).astype(float)
    rows = []
    for g in range(N):
        for c in sorted(set(clades[loc == g])):
            m = (loc == g) & (clades == c)
            k = lib.k2p[m][np.isfinite(lib.k2p[m])]
            q = np.quantile(k, [0.25, 0.5, 0.75]) if k.size else [np.nan] * 3
            rows.append(dict(clade=c, subgenome=g, n=int(m.sum()), kb=L[m].sum() / 1e3, q25=q[0], q50=q[1], q75=q[2],
                             n_after=float(np.nansum(after[g][m])) if g in after else None))
    return rows


def run_biology(a, lib, ph, dt):
    """Derived biology: clade table, insertion history and nesting check (any N)."""
    rng = np.random.default_rng([a.seed, 11])
    N = max(ph["N"], 1)
    ek = ph["ek"]
    loc = np.where(ek >= 0, ph["z"][np.maximum(ek, 0)], -1) if ph["N"] > 1 else np.where(ek >= 0, 0, -1)
    kc, nc = ltr_counts(lib)
    if not (np.isfinite(kc) & (nc > 0)).any():
        log("WARNING: no LTR substitution counts or LTR lengths: no deconvolved insertion history, no nesting check, "
            "no copies-after-merger counts")
    nest = nesting_check(lib)
    nv = sum(r["violation"] for r in nest)
    if nest:
        log(f"nesting check: inner element older than its host in {nv} of {len(nest)} nested pairs "
            f"({100 * nv / len(nest):.1f}%; > 2 se)")
    history = insertion_history(lib, loc, N, rng)                  # drawn first (results as before elements.tsv)
    after = after_merger(lib, loc, dt, rng)
    return dict(nesting=nest, history=history, clades=clade_table(lib, loc, N, after), after=after,
                mu=dt["mu"] if dt else (a.mu or MU_DEFAULT), N=N)


def biology_summary(bio):
    """Compact biology record for summary.json: nesting check."""
    nv = sum(r["violation"] for r in bio["nesting"])
    return dict(nesting=dict(n_nested=len(bio["nesting"]), n_inner_older=nv,
                             rate=nv / len(bio["nesting"]) if bio["nesting"] else None))


def write_biology(a, lib, bio):
    od = Path(a.outdir)
    yr = lambda k: _f(k / (2 * bio["mu"]) if k is not None and np.isfinite(k) else None)
    with open(od / "clades.tsv", "w") as f:
        f.write("clade\tsubgenome\tn_copies\tkb\tk2p_q25\tk2p_median\tk2p_q75\tyears_median\tn_after_merger\n")
        for r in bio["clades"]:
            f.write(f"{r['clade']}\t{group_label(r['subgenome'])}\t{r['n']}\t{r['kb']:.1f}\t{_f(r['q25'])}\t"
                    f"{_f(r['q50'])}\t{_f(r['q75'])}\t{yr(r['q50'])}\t{_f(r['n_after'], '{:.1f}')}\n")
    with open(od / "insertion_history.tsv", "w") as f:
        f.write(f"# mu={bio['mu']:g}; observed = copies per LTR K2P bin; deconvolved = smooth g-model fit (LTR "
                f"counting noise removed; survival and detection biases remain)\n")
        f.write("subgenome\tclade\tk2p_low\tk2p_high\tyears_low\tyears_high\tobserved\tdeconvolved\n")
        for r in bio["history"]:
            f.write(f"{group_label(r['subgenome'])}\t{r['clade']}\t{r['lo']:.4g}\t{r['hi']:.4g}\t{yr(r['lo'])}\t"
                    f"{yr(r['hi'])}\t{r['observed']}\t{_f(r['deconvolved'], '{:.2f}')}\n")
    with open(od / "nesting.tsv", "w") as f:
        f.write("inner\thost\tinner_k2p\thost_k2p\tdifference\tse\tinner_older\tin_host_ltr\n")
        for r in bio["nesting"]:
            f.write(f"{lib.ids[r['inner']]}\t{lib.ids[r['host']]}\t{r['k_inner']:.5f}\t{r['k_host']:.5f}\t"
                    f"{r['k_inner'] - r['k_host']:+.5f}\t{r['se']:.5f}\t{'yes' if r['violation'] else 'no'}\t"
                    f"{'yes' if r['in_host_ltr'] else 'no'}\n")


# ------------------------------------------------------------------------------------------------ painting
def age_weights(lib, own, lP, ninf, N):
    """Weight in [0, 1] of each copy as evidence about its neighborhood's origin (v5): 1 at LTR ages whose copies are
    most confined to their own subgenome, 0 at ages whose copies carry other-subgenome ancestry as often as
    uninformative copies would (1 - prior of their own subgenome: post-merger insertions, lineages older than the
    split); copies without K2P get 1; all 1 when no age class is confined by >= 0.05 or all ages are equal."""
    use = (own >= 0) & (ninf > 0) & np.isfinite(lib.k2p)
    w = np.ones(lib.n)
    if use.sum() < 100:
        return w
    prior = np.bincount(own[own >= 0], minlength=N) / max((own >= 0).sum(), 1)
    share = 1 - np.exp(lP[np.arange(lib.n), np.maximum(own, 0)])
    s_null = float(np.mean(1 - prior[own[use]]))
    k = lib.k2p[use]
    edges = np.unique(np.quantile(k, np.linspace(0, 1, AGE_BINS + 1)))
    if edges.size < 2:                                          # every copy of one age: no age classes to compare
        return w
    bi = np.clip(np.searchsorted(edges, k, "right") - 1, 0, edges.size - 2)
    cnt = np.bincount(bi, minlength=edges.size - 1)
    ok = cnt > 0
    A = np.array([np.median(k[bi == j]) for j in np.flatnonzero(ok)])
    Y = np.bincount(bi, weights=share[use], minlength=edges.size - 1)[ok] / cnt[ok]
    if s_null - Y.min() < 0.05:
        return w
    w = np.clip((s_null - np.interp(np.nan_to_num(lib.k2p, nan=-1.0), A, Y)) / (s_null - Y.min()), 0.0, 1.0)
    w[~np.isfinite(lib.k2p)] = 1.0
    return w
@njit(cache=True)
def _lse(a):
    m = a.max()
    return m + np.log(np.exp(a - m).sum())


@njit(cache=True)
def hmm_posterior(E, l0, sw):
    """Forward-backward over copies (log emissions E, n x N; log initial l0), switch probability sw per step split
    evenly over the other states -> (posteriors n x N, log-likelihood, expected number of switches)."""
    n, N = E.shape
    lsw = np.log(sw / max(N - 1, 1))
    lst = np.log1p(-sw)
    f = np.empty((n, N))
    b = np.zeros((n, N))
    f[0] = l0 + E[0]
    tmp = np.empty(N)
    for t in range(1, n):
        for j in range(N):
            for i in range(N):
                tmp[i] = f[t - 1, i] + (lst if i == j else lsw)
            f[t, j] = _lse(tmp) + E[t, j]
    for t in range(n - 2, -1, -1):
        for i in range(N):
            for j in range(N):
                tmp[j] = (lst if i == j else lsw) + E[t + 1, j] + b[t + 1, j]
            b[t, i] = _lse(tmp)
    ll = _lse(f[n - 1])
    post = np.exp(f + b - ll)
    esw = 0.0
    for t in range(n - 1):
        for i in range(N):
            for j in range(N):
                if i != j:
                    esw += np.exp(f[t, i] + lsw + E[t + 1, j] + b[t + 1, j] - ll)
    return post, ll, esw


@njit(cache=True)
def hmm_viterbi(E, l0, sw):
    """Most likely state path (same model as hmm_posterior)."""
    n, N = E.shape
    lsw = np.log(sw / max(N - 1, 1))
    lst = np.log1p(-sw)
    v = l0 + E[0]
    bp = np.zeros((n, N), np.int64)
    nv = np.empty(N)
    for t in range(1, n):
        for j in range(N):
            best = -1e300
            arg = 0
            for i in range(N):
                s = v[i] + (lst if i == j else lsw)
                if s > best:
                    best = s
                    arg = i
            nv[j] = best + E[t, j]
            bp[t, j] = arg
        v[:] = nv
    path = np.empty(n, np.int64)
    path[n - 1] = np.argmax(v)
    for t in range(n - 1, 0, -1):
        path[t - 1] = bp[t, path[t]]
    return path


def switch_rate(seqs, N):
    """One genome-wide switch probability per copy step by Baum-Welch (emissions fixed; v5)."""
    sw = 1e-3
    for _ in range(30):
        es = tot = 0.0
        for E, l0 in seqs:
            es += hmm_posterior(E, l0, sw)[2]
            tot += E.shape[0] - 1
        new = float(np.clip(es / max(tot, 1.0), 1e-7, 0.2))
        done = abs(new - sw) <= 0.01 * sw
        sw = new
        if done:
            break
    return sw


class Lineages:
    """Copies sharing a k-mer clade belong to one lineage, one piece of evidence (v5): through small clades
    (<= NMAX_KNN carriers) anywhere, and through any clade whose carriers lie mostly (>= half) inside the segment being
    judged (a lineage amplified in one place has no small clade: its k-mers are private or family-wide)."""

    def __init__(self, cl, n):
        cid = np.repeat(np.arange(cl.n), cl.size)
        mem = cl.members.astype(np.int64)
        o = np.argsort(mem, kind="stable")
        self.mem, self.cid, self.size = mem[o], cid[o], cl.size
        self.ptr = np.searchsorted(self.mem, np.arange(n + 1))

    def labels(self, members):
        """Lineage label of each copy in members: connected components through qualifying shared clades."""
        members = np.asarray(members)
        parent = np.arange(members.size)

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        hits = {}
        for p, c in enumerate(members):
            for k in self.cid[self.ptr[c]:self.ptr[c + 1]]:
                hits.setdefault(int(k), []).append(p)
        for k, ps in hits.items():
            if len(ps) < 2 or (self.size[k] > NMAX_KNN and 2 * len(ps) < self.size[k]):
                continue
            r0 = find(ps[0])
            for p in ps[1:]:
                r = find(p)
                if r != r0:
                    parent[r] = r0
        return np.array([find(i) for i in range(members.size)], np.int64)


def resolution_track(mid, n_needed, L):
    """Smallest exchange detectable around each RES_WIN window of a sequence: the shortest span of n_needed consecutive
    informative copies (midpoints, sorted here) that covers the window's center -> (starts, ends, bp); no rows for
    windows centred outside the copies, none at all when the sequence holds fewer copies than needed."""
    mid = np.sort(np.asarray(mid, float))
    starts = np.arange(0, L, RES_WIN, dtype=np.int64)
    ends = np.minimum(starts + RES_WIN, L)
    if mid.size < n_needed:
        return starts[:0], ends[:0], np.zeros(0)
    span = mid[n_needed - 1:] - mid[:mid.size - n_needed + 1]            # copies j .. j + n_needed - 1
    x = (starts + ends) / 2
    j0 = np.maximum(np.searchsorted(mid, x, "left") - n_needed + 1, 0)   # last copy at or after the center
    j1 = np.minimum(np.searchsorted(mid, x, "right") - 1, span.size - 1)  # first copy at or before it
    keep = (x >= mid[0]) & (x <= mid[-1]) & (j1 >= j0)
    out = np.array([span[a:b + 1].min() for a, b in zip(j0[keep], j1[keep])])
    return starts[keep], ends[keep], out


def posterior_tiles(mid, post, L, R):
    """Ancestry posterior tiles of one sequence for the bedGraphs: each informative copy (midpoint order) owns the
    stretch to the midpoints between it and its neighbours, but at most R/2 on either side; longer stretches without
    copies stay unknown (no tile) -> (starts, ends, posteriors)."""
    o = np.argsort(mid, kind="stable")
    m, p = np.asarray(mid, float)[o], post[o]
    half = R / 2
    lo = np.empty(m.size)
    hi = np.empty(m.size)
    if m.size:
        between = (m[1:] + m[:-1]) / 2
        lo[0], hi[-1] = max(0.0, m[0] - half), min(float(L), m[-1] + half)
        lo[1:] = np.maximum(between, m[1:] - half)
        hi[:-1] = np.minimum(between, m[:-1] + half)
    lo, hi = np.floor(lo).astype(np.int64), np.floor(hi).astype(np.int64)
    keep = hi > lo
    return lo[keep], hi[keep], p[keep]


def run_painting(a, lib, ph, dt):
    """LTR-only ancestry along each phased sequence (v5): HMM over informative copies, emissions = calibrated per-copy
    ancestry (log P - log prior) x age weight, one genome-wide switch rate (Baum-Welch); Viterbi segments of another
    subgenome are called exchanges when their likelihood ratio, each lineage counted once, is >= BF_DECISIVE and >=
    MIN_LINEAGES lineages favour the origin >= BF_SUPPORT-fold. Resolution: smallest detectable exchange per window
    (the span of the informative copies whose typical evidence pays the two switches and the call threshold)."""
    if ph["N"] < 2:
        return None
    N, z, ek = ph["N"], ph["z"], ph["ek"]
    ancestry = dt.get("ancestry") if dt else None
    if ancestry is None:
        ancestry = copy_ancestry(ph["clades"], ek, z, N, np.random.default_rng([a.seed, 7]))
    lP, ninf = ancestry
    own = np.where(ek >= 0, z[np.maximum(ek, 0)], -1)
    prior = np.bincount(own[own >= 0], minlength=N) / max((own >= 0).sum(), 1)
    w = age_weights(lib, own, lP, ninf, N)
    inf = (ninf > 0) & (ek >= 0)
    si = {s: i for i, s in enumerate(lib.seqnames)}
    seqs, meta = [], []
    for i, s in enumerate(ph["names"]):
        idx = np.flatnonzero(inf & (ek == i))
        if idx.size < 2:
            continue
        idx = idx[np.argsort(lib.start[idx], kind="stable")]
        E = np.ascontiguousarray((lP[idx] - np.log(prior)[None]) * w[idx][:, None])
        l0 = np.log(np.where(np.arange(N) == z[i], 1 - START_OTHER, START_OTHER / (N - 1)))
        seqs.append((E, l0))
        meta.append((i, s, idx))
    sw = switch_rate(seqs, N) if seqs else 1e-3
    lin = Lineages(ph["clades"], lib.n)
    home = [E[np.arange(E.shape[0]), int(z[i])] - np.max(np.delete(E, int(z[i]), 1), 1)
            for (E, _), (i, _, _) in zip(seqs, meta)]                  # evidence of copies for their own subgenome
    pos = np.concatenate([h[h > 0] for h in home]) if home else np.zeros(0)
    # a segment enters the best path only when its evidence pays both switches (in and out) and reaches the call
    # threshold: copies needed = that log-evidence / the median evidence of a copy for its own subgenome
    need = max(math.log(BF_DECISIVE), 2 * (math.log1p(-sw) - math.log(sw / (N - 1))))
    n_needed = int(max(MIN_LINEAGES, math.ceil(need / np.median(pos)))) if pos.size else MIN_LINEAGES
    segs, post_tracks, res = [], [], []
    for (E, l0), (i, s, idx) in zip(seqs, meta):
        g0 = int(z[i])
        post, _, _ = hmm_posterior(E, l0, sw)
        path = hmm_viterbi(E, l0, sw)
        mid = (lib.start[idx] + lib.end[idx]) // 2
        L = int(lib.seq_len[si[s]])
        post_tracks.append((s, mid, post, L))
        res.append((s,) + resolution_track(mid.astype(float), n_needed, L))
        brk = np.flatnonzero(np.diff(path)) + 1
        for s0, e0 in zip(np.concatenate(([0], brk)), np.concatenate((brk, [idx.size]))):
            g = int(path[s0])
            if g == g0:
                continue
            ii = idx[s0:e0]
            d = E[s0:e0, g] - E[s0:e0, g0]
            lab = lin.labels(ii)
            _, inv = np.unique(lab, return_inverse=True)
            llr = float((np.bincount(inv, weights=d) / np.bincount(inv)).sum())
            sup = d >= math.log(BF_SUPPORT)
            nlin = int(np.unique(lab[sup]).size)
            segs.append(dict(seq=s, start=int(lib.start[ii].min()), end=int(lib.end[ii].max()), own=g0, origin=g,
                             left_flank=int(lib.end[idx[s0 - 1]]) if s0 > 0 else 1,
                             right_flank=int(lib.start[idx[e0]]) if e0 < idx.size else L,
                             n_ltr=int(ii.size), n_support=int(sup.sum()), n_against=int((d <= -math.log(BF_SUPPORT)).sum()),
                             n_lineages=nlin, llr=llr, post=float(post[s0:e0, g].mean()),
                             called=bool(llr >= math.log(BF_DECISIVE) and nlin >= MIN_LINEAGES)))
    nc = sum(x["called"] for x in segs)
    unph = unphased_ancestry(a, lib, ph, prior, N)
    rr = np.concatenate([r[3] for r in res]) if res else np.zeros(0)
    log(f"painting: switch {sw:.2e} per copy; {nc} exchange segments called (likelihood ratio >= {BF_DECISIVE:g} with "
        f"each lineage counted once, >= {MIN_LINEAGES} lineages); smallest detectable exchange: median "
        f"{np.median(rr) / 1e3:,.0f} kb ({n_needed} informative copies)" if rr.size else f"painting: {nc} exchange segments called")
    return dict(segs=segs, switch=sw, n_needed=n_needed, post=post_tracks, res=res, N=N, ancestry=ancestry,
                unphased=unph)


def unphased_ancestry(a, lib, ph, prior, N):
    """Ancestry of each sequence too short or too poor in LTR-RTs to phase, from its informative copies (v5): the
    copies' evidence (log posterior over the prior) summed over lineages, each lineage counted once; called when the
    best subgenome beats the next >= BF_DECISIVE-fold with >= MIN_LINEAGES lineages favoring it >= BF_SUPPORT-fold
    ('decisive'), or >= BF_SUPPORT-fold ('substantial'); otherwise none. Phasing's clades hold phased copies only:
    the unphased copies' ancestry and lineages come from all clades (phased results are untouched)."""
    phased = set(ph["names"])
    ek = np.asarray(ph["ek"])
    if not (ek < 0).any() or "clades_all" not in ph:
        return []
    lP, ninf = copy_ancestry(ph["clades_all"], ek, np.asarray(ph["z"]), N, np.random.default_rng([a.seed, 13]))
    lin = Lineages(ph["clades_all"], lib.n)
    E = lP - np.log(prior)[None]
    o = np.argsort(lib.seq, kind="stable")
    cuts = np.flatnonzero(np.diff(lib.seq[o])) + 1
    out = []
    for grp in np.split(o, cuts):
        if not grp.size or lib.seqnames[lib.seq[grp[0]]] in phased:
            continue
        idx = grp[ninf[grp] > 0]
        row = dict(seq=lib.seqnames[lib.seq[grp[0]]], n_ltr=int(grp.size), n_inf=int(idx.size), n_lin=0, g=-1, fav=-1,
                   lo=0.0, conf="none")
        if idx.size:
            lab = lin.labels(idx)
            _, inv = np.unique(lab, return_inverse=True)
            cnt = np.bincount(inv)
            tot = np.array([(np.bincount(inv, weights=E[idx, g]) / cnt).sum() for g in range(N)])
            best = np.argsort(-tot)
            g, lo = int(best[0]), float(tot[best[0]] - tot[best[1]])
            d = E[idx, g] - np.max(np.delete(E[idx], g, 1), 1)
            nlin = int(np.unique(lab[d >= math.log(BF_SUPPORT)]).size)
            conf = ("decisive" if lo >= math.log(BF_DECISIVE) and nlin >= MIN_LINEAGES else
                    "substantial" if lo >= math.log(BF_SUPPORT) else "none")
            row.update(n_lin=nlin, g=g if conf != "none" else -1, fav=g, lo=lo, conf=conf)
        out.append(row)
    return out


def painting_summary(pa):
    """Compact painting record for summary.json: switch rate, calls, exchanged bp per origin -> recipient."""
    bp = {}
    for x in pa["segs"]:
        if x["called"]:
            k = f"{group_label(x['origin'])}>{group_label(x['own'])}"
            bp[k] = bp.get(k, 0) + x["end"] - x["start"] + 1
    rr = np.concatenate([r[3] for r in pa["res"]]) if pa["res"] else np.zeros(0)
    return dict(switch=pa["switch"], n_segments=len(pa["segs"]), n_called=sum(x["called"] for x in pa["segs"]),
                n_unphased=len(pa.get("unphased", [])),
                n_unphased_called=sum(u["conf"] != "none" for u in pa.get("unphased", [])),
                exchanged_bp=bp, copies_needed=pa["n_needed"],
                resolution_median_bp=float(np.median(rr)) if rr.size else None)


def write_unphased(od, rows):
    """unphased_sequences.tsv: ancestry (called subgenome or NA), favored (the subgenome the odds favor, NA without
    informative copies), log10 odds of the favored over the next subgenome (NA without informative copies)."""
    with open(Path(od) / "unphased_sequences.tsv", "w") as f:
        f.write("sequence\tn_ltr\tn_informative\tn_lineages\tancestry\tfavored\tlog10_odds_vs_next\tconfidence\n")
        for u in rows:
            ev = u["n_inf"] > 0
            f.write(f"{u['seq']}\t{u['n_ltr']}\t{u['n_inf']}\t{u['n_lin']}\t"
                    f"{group_label(u['g']) if u['g'] >= 0 else 'NA'}\t{group_label(u['fav']) if ev else 'NA'}\t"
                    f"{u['lo'] / math.log(10):.2f}\t{u['conf']}\n" if ev else
                    f"{u['seq']}\t{u['n_ltr']}\t0\t0\tNA\tNA\tNA\t{u['conf']}\n")


def write_painting(a, pa):
    od = Path(a.outdir)
    if pa.get("unphased"):
        write_unphased(od, pa["unphased"])
    with open(od / "segments.tsv", "w") as f:
        f.write("sequence\tstart\tend\tlength\tleft_flank\tright_flank\tsubgenome\torigin\tn_ltr\tn_support\t"
                "n_against\tn_lineages\tlog10_lr_lineages\tmean_posterior\tcalled\n")
        for x in pa["segs"]:
            f.write(f"{x['seq']}\t{x['start']}\t{x['end']}\t{x['end'] - x['start'] + 1}\t{x['left_flank']}\t"
                    f"{x['right_flank']}\t{group_label(x['own'])}\t"
                    f"{group_label(x['origin'])}\t{x['n_ltr']}\t{x['n_support']}\t{x['n_against']}\t{x['n_lineages']}\t"
                    f"{x['llr'] / math.log(10):.2f}\t{x['post']:.3f}\t{'yes' if x['called'] else 'no'}\n")
    called = [x for x in pa["segs"] if x["called"]]
    with open(od / "exchanges.bed", "w") as f:
        for k, x in enumerate(called):
            f.write(f"{x['seq']}\t{x['start'] - 1}\t{x['end']}\texchange_{k + 1}:{group_label(x['origin'])}>"
                    f"{group_label(x['own'])}\t{min(1000, int(100 * x['llr'] / math.log(10)))}\t.\n")
    with open(od / "exchanges.gff3", "w") as f:
        f.write("##gff-version 3\n")
        for k, x in enumerate(called):
            f.write(f"{x['seq']}\tv9\tregion\t{x['start']}\t{x['end']}\t{x['llr'] / math.log(10):.2f}\t.\t.\t"
                    f"ID=exchange_{k + 1};origin={group_label(x['origin'])};recipient={group_label(x['own'])};"
                    f"lineages={x['n_lineages']};copies={x['n_ltr']}\n")
    rr = np.concatenate([r[3] for r in pa["res"]]) if pa["res"] else np.zeros(0)
    R = float(np.median(rr)) if rr.size else float(RES_WIN)             # longer copy-free stretches stay unknown
    for g in range(pa["N"]):
        with open(od / f"ancestry_{group_label(g)}.bedGraph", "w") as f:
            f.write(f"track type=bedGraph name=\"P({group_label(g)})\" description=\"posterior ancestry {group_label(g)}"
                    f" per informative LTR-RT (tiles between copy midpoints; copy-free stretches longer than the median"
                    f" resolution, {R / 1e3:,.0f} kb, unknown)\"\n")
            for s, mid, post, L in pa["post"]:
                for p_, q_, v_ in zip(*posterior_tiles(mid, post, L, R)):
                    f.write(f"{s}\t{p_}\t{q_}\t{v_[g]:.4f}\n")
    with open(od / "resolution.bedGraph", "w") as f:
        f.write(f"track type=bedGraph name=\"resolution\" description=\"smallest detectable exchange (bp): span of "
                f"{pa['n_needed']} informative LTR-RTs\"\n")
        for s, st, en, v in pa["res"]:
            for p, q, r in zip(st, en, v):
                f.write(f"{s}\t{p}\t{q}\t{int(r)}\n")


def _favors(lib, ph, anc):
    """Each copy's evidence as the painting reads it (log posterior over the subgenome's share of copies): informative
    copies; those whose lineages favor another subgenome >= BF_SUPPORT-fold over their sequence's own (foreign), and
    their own >= BF_SUPPORT-fold over every other (home); the most favored other subgenome."""
    lP, ninf = anc
    N, ek, z = ph["N"], np.asarray(ph["ek"]), np.asarray(ph["z"])
    own = np.where(ek >= 0, z[np.maximum(ek, 0)], 0)
    prior = np.bincount(own[ek >= 0], minlength=N) / max((ek >= 0).sum(), 1)
    E = lP - np.log(np.maximum(prior, 1e-300))[None]
    r = np.arange(lib.n)
    Eo = E[r, own]
    Ex = E.copy()
    Ex[r, own] = -np.inf
    gbest, d = Ex.argmax(1), Ex.max(1) - Eo
    inf = (ninf > 0) & (ek >= 0)
    return inf, inf & (d >= math.log(BF_SUPPORT)), inf & (-d >= math.log(BF_SUPPORT)), gbest


# ------------------------------------------------------------------------------------------------ report.pdf
# One page, no legend pages: every panel has a one-line header and its own key, so the page reads alone. A: the
# painted karyotype across the top. Left column: B lineage sharing between sequences, C each sequence's placement in
# independent half-libraries, D the two tests behind the number of subgenomes. Right column, on one time axis: E
# transposon births per subgenome, F per main clade, G sharing in young, middle-aged and old copies, H closest
# relatives across subgenomes with every split and merger estimate. Nothing is re-estimated for the figure.
FIG_H = 9.45                  # figure height (inches): 240 mm (FIG_W x FIG_H fits an A4 page)
FIG_ALPHA = (0.22, 0.45, 0.72, 1.0)   # painting opacity by the posterior of the painted subgenome: <.7, <.9, <.99, more
FIG_SPLIT, FIG_MERGE, FIG_INK = "#e9e9e9", "#f7edd2", "#222222"   # split band, merger band, text and marks
FIG_EST = {"orthologs": "orthologs", "confinement": "confinement", "profile": "profile", "ltr_clock": "LTR clock",
           "genealogy_clock": "genealogy clock"}           # dating.tsv's estimators
FIG_NUM = dict(enumerate(("no", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten",
                          "eleven", "twelve")))
FIG_XL, FIG_XR = 0.22, 3.12   # panel letters of the left and right columns (inches from the left edge)


def _fig_style():
    """matplotlib (imported here, when report.pdf is drawn) with Arial if installed, else Liberation Sans from
    the system font directory, else the default sans; 7-pt text; fonttype 42 keeps text as text in PDFs."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager as fm
    names = {f.name for f in fm.fontManager.ttflist}
    family = "Arial" if "Arial" in names else None
    if family is None:
        for d in ("/usr/share/fonts/truetype/liberation2", "/usr/share/fonts/truetype/liberation",
                  "/usr/share/fonts/liberation-sans", "/usr/share/fonts/liberation", "/Library/Fonts"):
            for st in ("Regular", "Bold", "Italic", "BoldItalic"):
                p = os.path.join(d, f"LiberationSans-{st}.ttf")
                if os.path.exists(p):
                    fm.fontManager.addfont(p)
        if "Liberation Sans" in {f.name for f in fm.fontManager.ttflist}:
            family = "Liberation Sans"
    family = family or "DejaVu Sans"
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": [family], "font.size": 7,
                         "axes.labelsize": 6.5, "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6.5,
                         "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                         "xtick.major.size": 2.5, "ytick.major.size": 2.5, "pdf.fonttype": 42, "ps.fonttype": 42,
                         "axes.spines.top": False, "axes.spines.right": False, "mathtext.default": "regular",
                         "hatch.linewidth": 0.4, "savefig.dpi": 300})
    return plt, family


def _ax(fig, x0, y0, w, h):
    """Axes at x0, y0 (inches from the lower left), w x h inches: every panel sits on a grid measured in inches."""
    W, H = fig.get_size_inches()
    return fig.add_axes([x0 / W, y0 / H, w / W, h / H])


def _head(fig, ax, letter, text, xl, dy=0.06):
    """Bold panel letter at xl inches from the left edge and a one-line header after it, dy inches above the axes."""
    W, H = fig.get_size_inches()
    y = ax.get_position().y1 + dy / H
    fig.text(xl / W, y, letter, fontsize=9, fontweight="bold", va="baseline", ha="left")
    fig.text((xl + 0.17) / W, y, text, fontsize=7, va="baseline", ha="left", color=FIG_INK)


def _text_w(fig, s, size, **kw):
    """Width of a text in inches."""
    t = fig.text(0, 0, s, fontsize=size, **kw)
    w = t.get_window_extent(fig.canvas.get_renderer()).width / fig.dpi
    t.remove()
    return w


def _key_row(fig, x, y, items, size=6.2, gap=0.22):
    """On-figure key from x, y (inches, row centre): items (kind, color, label) with kind box, hatch, tri, dens, dot,
    open, num (a homoeolog set number) or frame (an outlined cell); returns the x where it ends."""
    from matplotlib.lines import Line2D
    from matplotlib.patches import Polygon, Rectangle
    W, H = fig.get_size_inches()
    for kind, col, lab in items:
        g = 0.13
        if kind in ("box", "hatch"):
            fig.add_artist(Rectangle((x / W, (y - 0.045) / H), g / W, 0.09 / H, transform=fig.transFigure,
                                     facecolor="white" if kind == "hatch" else col,
                                     edgecolor="#c8c8c8" if kind == "hatch" else "none",
                                     hatch="//////" if kind == "hatch" else None, lw=0))
        elif kind == "tri":
            fig.add_artist(Polygon([((x + 0.01) / W, y / H), ((x + g) / W, (y + 0.05) / H), ((x + g) / W, (y - 0.05) / H)],
                                   closed=True, transform=fig.transFigure, facecolor=col, edgecolor="black", lw=0.3))
        elif kind == "dens":
            for dx_, w_ in ((0.0, 0.05), (0.06, 0.03)):
                fig.add_artist(Rectangle(((x + dx_) / W, (y - 0.05) / H), w_ / W, 0.1 / H, transform=fig.transFigure,
                                         facecolor="#8c8c8c", lw=0))
        elif kind in ("dot", "open"):
            fig.add_artist(Line2D([(x + g / 2) / W], [y / H], marker="o", ms=3.4, lw=0,
                                  mfc=col if kind == "dot" else "white", mec=col, mew=0.7))
        elif kind == "num":
            fig.text((x + g / 2) / W, y / H, "1", ha="center", va="center", fontsize=size)
        elif kind == "frame":
            fig.add_artist(Rectangle(((x + 0.015) / W, (y - 0.045) / H), (g - 0.03) / W, 0.09 / H,
                                     transform=fig.transFigure, fill=False, edgecolor=col, lw=0.6))
        fig.text((x + g + 0.04) / W, y / H, lab, va="center", ha="left", fontsize=size)
        x += g + 0.04 + _text_w(fig, lab, size) + gap
    return x


def _yr_label(v):
    return f"{v / 1e6:g} Myr" if v >= 1e6 else f"{v / 1e3:g} kyr" if v >= 1e3 else f"{v:g} yr"


def _rate(mu):
    """Plain-language substitution rate: 7e-9 -> '7 substitutions per site per billion years'."""
    return f"{mu * 1e9:.3g} substitutions per site per billion years"


def _short_name(s, keep=30):
    """Sequence names longer than keep characters shortened in the middle for figure labels (tables keep them)."""
    return s if len(s) <= keep else s[:keep // 2 - 2] + "…" + s[-(keep // 2 - 1):]


def _ma(t):
    """Readable age from years."""
    return f"{t / 1e6:.2f} Myr" if t >= 1e6 else f"{t / 1e3:,.0f} kyr" if t >= 1e3 else f"{t:,.0f} years"


def _span(lo, est, hi):
    """'4.32 Myr (2.38–18.5)' when the limits share the estimate's unit, else each limit with its own unit."""
    unit = lambda t: "Myr" if t >= 1e6 else "kyr" if t >= 1e3 else "yr"
    num = lambda t: f"{t / 1e6:.3g}" if t >= 1e6 else f"{t / 1e3:.3g}" if t >= 1e3 else f"{t:.3g}"
    if unit(lo) == unit(est) == unit(hi):
        return f"{num(est)} {unit(est)} ({num(lo)}–{num(hi)})"
    return f"{num(est)} {unit(est)} ({num(lo)} {unit(lo)} – {num(hi)} {unit(hi)})"


def _homoeolog_numbers(ph):
    """v9: number of each sequence's homoeolog set (1..K, sets ordered by their first member in their lowest
    subgenome); empty without homoeolog calls."""
    ho = ph.get("homoeologs")
    if not ho or not ho["sets"]:
        return {}
    z = np.asarray(ph["z"])
    rank = {i: r for r, i in enumerate(sorted(range(len(ph["names"])), key=lambda i: natural_key(ph["names"][i])))}

    def first(s_):
        g0 = min(int(z[x]) for x in s_["members"])
        return g0, min(rank[x] for x in s_["members"] if z[x] == g0)
    return {x: k + 1 for k, s_ in enumerate(sorted(ho["sets"], key=first)) for x in s_["members"]}


def _shown(lib, ph):
    """Phased sequences the figures draw, grouped by subgenome: all of them, or the FIG_MAX_SEQ longest. Within a
    subgenome, members of homoeolog sets come first in set order (homoeologs share a position in every group), then
    the rest in natural order -> (indices into ph['names'], lengths in bp of every phased sequence)."""
    si = {s: i for i, s in enumerate(lib.seqnames)}
    L = np.array([lib.seq_len[si[s]] for s in ph["names"]], float)
    keep = np.argsort(-L, kind="stable")[:FIG_MAX_SEQ]
    z = np.asarray(ph["z"])
    num = _homoeolog_numbers(ph)
    return sorted(keep.tolist(), key=lambda i: (z[i], num.get(i, 1 << 30), natural_key(ph["names"][i]))), L


def _oe(W):
    """Observed over expected lineage sharing (degree-corrected null), NaN on the diagonal and without expectation."""
    W = np.asarray(W, float)
    E = expected(W)
    with np.errstate(divide="ignore", invalid="ignore"):
        OE = np.where(E > 0, W / E, np.nan)
    np.fill_diagonal(OE, np.nan)
    return OE


def _split_name(t, zof, gsize):
    """A tested split in final-group terms: 'SG1 | SG2' when each side holds whole groups, else 'SG1: 7 | 3'."""
    sides = []
    for part in (t["part_a"], t["part_b"]):
        gs = sorted({zof[n] for n in part if n in zof})
        whole = bool(gs) and all(sum(zof.get(n) == g for n in part) == gsize[g] for g in gs)
        sides.append(("+".join(group_label(g) for g in gs), whole))
    if all(w for _, w in sides) and t["part_b"]:
        return f"{sides[0][0]} | {sides[1][0]}"
    gs = sorted({zof[n] for n in t["part_a"] + t["part_b"] if n in zof})
    return f"{'+'.join(group_label(g) for g in gs)}: {len(t['part_a'])} | {len(t['part_b'])}"


def _split_words(t, zof, gsize):
    """'SG1 vs SG2' for a split between whole groups, 'within SG1: 5 vs 2' for a split inside one."""
    nm = _split_name(t, zof, gsize)
    if ":" in nm:
        g, rest = nm.split(":", 1)
        return f"within {g.strip()}: {rest.strip().replace(' | ', ' vs ')}"
    return " vs ".join(sorted((x.strip() for x in nm.split("|")), key=natural_key))


def _wilson(k, n, zq=1.96):
    """95% Wilson interval of a binomial share."""
    n = np.maximum(np.asarray(n, float), 1.0)
    p = np.asarray(k, float) / n
    den = 1 + zq ** 2 / n
    c = (p + zq ** 2 / (2 * n)) / den
    h = zq * np.sqrt(p * (1 - p) / n + zq ** 2 / (4 * n ** 2)) / den
    return c - h, c + h


def _repel(ys, gap):
    """Label positions spread so neighbours are >= gap apart (order kept), staying near their anchors."""
    ys = np.asarray(ys, float)
    o = np.argsort(ys, kind="stable")
    y = ys[o].copy()
    for _ in range(200):
        moved = False
        for k in range(1, y.size):
            if y[k] - y[k - 1] < gap:
                push = (gap - (y[k] - y[k - 1])) / 2
                y[k - 1] -= push
                y[k] += push
                moved = True
        if not moved:
            break
    out = np.empty_like(y)
    out[o] = y
    return out


def _fig_data(a, lib, ph, dt, bio, pa):
    """Everything the figure draws, from the run's objects: sequences shown, painting runs, lineage sharing,
    closest-relative profiles, insertion histories and the range of the time axis (years)."""
    N, names, z, ek = ph["N"], list(ph["names"]), np.asarray(ph["z"]), np.asarray(ph["ek"])
    shown, L = _shown(lib, ph)
    mu = dt["mu"] if dt else (a.mu or MU_DEFAULT)
    P = dict(N=N, names=names, z=z, shown=shown, L=L, mu=mu, mu_given=bool(dt["mu_given"] if dt else a.mu), lib=lib,
             dotplot=N > 1 and bool(ph.get("homoeologs")),             # v9: panel B shows the homoeolog evidence
             ph=ph, dt=dt, bio=bio or {}, pa=pa, loc=np.where(ek >= 0, z[np.maximum(ek, 0)], -1),
             si={s: i for i, s in enumerate(lib.seqnames)}, label=a.label or Path(a.ltr_fasta[0]).name.split(".")[0],
             unit="chromosomes" if np.median(L) >= 5e6 else "sequences")
    # painting: posterior tiles merged into runs of one painted subgenome and confidence level; copy-free stretches
    # longer than the median resolution, or than ~2 pixel rows at this scale, stay unknown (bedGraph rule)
    tiles, R = {}, float(RES_WIN)
    if pa and pa.get("post"):
        rr = np.concatenate([r[3] for r in pa["res"]]) if pa.get("res") else np.zeros(0)
        R = float(np.median(rr)) if rr.size else float(RES_WIN)
        R_draw = max(R, float(L[shown].max()) / 400)
        for s, mid, post, Ls in pa["post"]:
            lo, hi, p = posterior_tiles(mid, post, Ls, R_draw)
            if not lo.size:
                continue
            c, lev = p.argmax(1), np.digitize(p.max(1), [0.7, 0.9, 0.99])
            key = c * 4 + lev
            brk = np.flatnonzero((np.diff(key) != 0) | (lo[1:] != hi[:-1])) + 1
            st, en = np.r_[0, brk], np.r_[brk, lo.size]
            tiles[s] = (lo[st], hi[en - 1], c[st], lev[st])
    P.update(tiles=tiles, R=R, calls=[c for c in (pa["segs"] if pa else []) if c["called"]])
    W_src = dict(ph.get("W_src") or {"all": ph["W_raw"]})
    P["oe"] = {s: _oe(W_)[np.ix_(shown, shown)] for s, W_ in W_src.items()}
    prof = []                                          # share of closest-relative pairs in two subgenomes, by divergence
    if dt and dt.get("pt") is not None:
        Pp, d = dt["pt"]["P"], dt["pt"]["d"]
        la, lb = P["loc"][Pp[:, 0]], P["loc"][Pp[:, 1]]
        cen = 10 ** ((PROF_EDGES[:-1] + PROF_EDGES[1:]) / 2)
        for r in dt["pairs"]:
            x, y = r["groups"]
            m0 = np.isin(la, [x, y]) & np.isin(lb, [x, y])
            if not m0.any():
                continue
            n, k = cross_profile(d[m0], (la[m0] != lb[m0]).astype(float))
            ok = n >= 20
            if ok.sum() >= 2:
                lo, hi = _wilson(k[ok], n[ok])
                prof.append(dict(pair=(x, y), t=cen[ok] / (2 * mu), n=n[ok], f=k[ok] / n[ok], lo=lo, hi=hi))
    P["prof"] = prof
    rows = defaultdict(list)                           # insertion history per (subgenome, clade)
    for r in P["bio"].get("history", []):
        rows[(int(r["subgenome"]), str(r["clade"]))].append((r["lo"], r["hi"], r["observed"], r["deconvolved"]))
    P["hist"] = {}
    for key, rr_ in rows.items():
        lo, hi, ob, de = (np.array(v, float) for v in zip(*rr_))
        t = np.where(lo > 0, np.sqrt(np.maximum(lo, 1e-300) * hi), hi / 2) / (2 * mu)
        P["hist"][key] = dict(t=t, y=de if np.isfinite(de).all() else ob, n=float(ob.sum()))
    old = [p_["t"].max() for p_ in prof]
    for r in (dt["pairs"] if dt else []):
        old += [v["hi"] / (2 * mu) for v in r["split"].values() if v and np.isfinite(v.get("hi", np.nan))]
    fin = lib.k2p[np.isfinite(lib.k2p)]
    old.append(float(np.quantile(fin, 0.995)) / (2 * mu) if fin.size else 1e6)
    t0 = 0.4e-4 / (2 * mu)                             # below the identical-LTR bin of the insertion histories
    P["trange"] = (t0, max(1.25 * max(old), 100 * t0))
    return P


# ------------------------------------------------------------------ A: painted karyotype
def _panel_karyotype(fig, ax, P):
    """Every shown sequence a vertical bar to scale, grouped by subgenome, painted by the subgenome its LTR-RT
    lineages point to (opacity: posterior), hatched where copies are too sparse to tell; LTR-RT density in gray to its
    left; called exchanges outlined, with a triangle in the donor's color (the longest >= 1 Mb labelled)."""
    from matplotlib.collections import PatchCollection
    from matplotlib.colors import to_rgba
    from matplotlib.patches import Polygon, Rectangle
    lib, shown, L, z, names = P["lib"], P["shown"], P["L"], P["z"], P["names"]
    xpos, x, prev = {}, 0.0, None
    for i in shown:
        if prev is not None and z[i] != z[prev]:
            x += 0.9
        xpos[i], x, prev = x, x + 1.0, i
    n = len(shown)
    bw = 0.36 if n <= 30 else 0.5
    Lmax = L[shown].max() / 1e6
    num = _homoeolog_numbers(P["ph"])
    ax.set_xlim(-0.75, x - 0.25)
    ax.set_ylim(Lmax * (1.085 if num else 1.03), -Lmax * 0.005)
    ax.spines["bottom"].set_visible(False)
    ax.set_xticks([])
    ax.set_ylabel("position (Mb)")
    pitch = ax.get_position().width * fig.get_size_inches()[0] / (x + 0.5)     # inches per column
    painted = bool(P["tiles"])
    rects, cols, dens = [], [], []
    step = Lmax / 120
    for i in shown:
        s, x0, Lm = names[i], xpos[i], L[i] / 1e6
        if painted:
            ax.add_patch(Rectangle((x0 - bw / 2, 0), bw, Lm, facecolor="white", edgecolor="#c8c8c8", hatch="//////",
                                   lw=0, zorder=1))
        else:
            ax.add_patch(Rectangle((x0 - bw / 2, 0), bw, Lm, facecolor=FIG_LIGHT, edgecolor="none", lw=0, zorder=1))
        if s in P["tiles"]:
            for a_, b_, c_, l_ in zip(*P["tiles"][s]):
                rects.append(Rectangle((x0 - bw / 2, a_ / 1e6), bw, (b_ - a_) / 1e6))
                cols.append(to_rgba(SG_COLORS[c_], FIG_ALPHA[l_]))
        ax.add_patch(Rectangle((x0 - bw / 2, 0), bw, Lm, fill=False, edgecolor=FIG_DARK, lw=0.45, zorder=4))
        if pitch >= 0.2:
            m = lib.seq == P["si"][s]
            e = np.arange(0, Lm + step, step)
            h, _ = np.histogram((lib.start[m] + lib.end[m]) / 2e6, e)
            dens.append((x0, e, h / step))
    if rects:
        ax.add_collection(PatchCollection(rects, facecolors=cols, edgecolors="none", zorder=2,
                                          rasterized=len(rects) > 20000))
    if dens:
        dmax = max(h.max() for _, _, h in dens) or 1.0
        for x0, e, h in dens:
            ax.fill_betweenx(np.repeat(e, 2)[1:-1], x0 - bw / 2 - 0.04 - np.repeat(h, 2) / dmax * 0.3,
                             x0 - bw / 2 - 0.04, color="#8c8c8c", lw=0, zorder=1)
    rot, fs, _ = _name_layout(fig, P)
    for i in shown:
        ax.text(xpos[i], -Lmax * 0.012, _short_name(names[i]), ha="center", va="bottom", rotation=rot, fontsize=fs,
                clip_on=False)
        if i in num:                                                   # v9: homoeolog set number
            ax.text(xpos[i], Lmax * 1.055, str(num[i]), ha="center", va="center", fontsize=fs - 0.4, color=FIG_INK)
    calls = [c for c in P["calls"] if names.index(c["seq"]) in xpos]
    big = sorted((c for c in calls if c["end"] - c["start"] >= 1e6), key=lambda c: c["start"] - c["end"])[:8]
    for c in calls:
        x0, ym = xpos[names.index(c["seq"])], (c["start"] + c["end"]) / 2e6
        h = max((c["end"] - c["start"]) / 1e6, Lmax * 0.006)
        ax.add_patch(Rectangle((x0 - bw / 2, ym - h / 2), bw, h, fill=False, edgecolor="black", lw=0.5, zorder=5))
        ax.add_patch(Polygon([[x0 + bw / 2 + 0.03, ym], [x0 + bw / 2 + 0.2, ym - Lmax * 0.009],
                              [x0 + bw / 2 + 0.2, ym + Lmax * 0.009]], closed=True, facecolor=SG_COLORS[c["origin"]],
                             edgecolor="black", lw=0.3, zorder=5))
        if any(c is b for b in big) and pitch >= 0.3:
            ax.text(x0 + bw / 2 + 0.24, ym, f"{(c['end'] - c['start']) / 1e6:.1f} Mb", fontsize=5.5, va="center",
                    ha="left", zorder=6)
    return xpos


def _name_layout(fig, P):
    """Sequence names over the karyotype: horizontal when few and short, else vertical; 6 pt (5 pt beyond 30
    sequences) -> (rotation, font size, height in inches)."""
    nm = [_short_name(P["names"][i]) for i in P["shown"]]
    fs = 6 if len(nm) <= 30 else 5
    rot = 0 if len(nm) <= 16 and max(map(len, nm)) <= 8 else 90
    return rot, fs, (0.1 if rot == 0 else max(_text_w(fig, x, fs) for x in nm) + 0.02)


def _group_heads(fig, ax, P, xpos, y_in):
    """Each subgenome's columns under a colored rule and 'SG1 · 7 chromosomes · 1,119 Mb' (the name alone where the
    text would not fit)."""
    from matplotlib.lines import Line2D
    W, H = fig.get_size_inches()
    z, L, shown = P["z"], P["L"], P["shown"]
    bb = ax.get_position()
    xl0, xl1 = ax.get_xlim()
    xf = lambda xd: bb.x0 + (xd - xl0) / (xl1 - xl0) * bb.width
    for g in range(P["N"]):
        xs = [xpos[i] for i in shown if z[i] == g]
        if not xs:
            continue
        col = SG_COLORS[g] if P["N"] > 1 else FIG_DARK
        fig.add_artist(Line2D([xf(min(xs) - 0.4), xf(max(xs) + 0.4)], [y_in / H] * 2, color=col, lw=2,
                              solid_capstyle="butt"))
        n_g = int((z == g).sum())
        lab = (f"{group_label(g)} · {n_g} {P['unit']} · {L[z == g].sum() / 1e6:,.0f} Mb" if P["N"] > 1 else
               f"one group · {n_g} {P['unit']} · {L.sum() / 1e6:,.0f} Mb")
        room = (xf(max(xs) + 0.5) - xf(min(xs) - 0.5)) * W
        if _text_w(fig, lab, 7, fontweight="bold") > room:
            lab = group_label(g) if P["N"] > 1 else "one group"
        fig.text(xf((min(xs) + max(xs)) / 2), (y_in + 0.05) / H, lab, ha="center", va="bottom", color=col, fontsize=7,
                 fontweight="bold")


# ------------------------------------------------------------------ B: lineage sharing
def _panel_sharing(fig, ax, P, cmap, norm):
    """Lineage sharing between every pair of shown sequences (observed / expected, all copies), grouped by subgenome
    (strips: subgenome; white lines: group boundaries); names at the right when they fit."""
    from matplotlib.patches import Rectangle
    z, shown = P["z"], P["shown"]
    n = len(shown)
    ax.imshow(np.clip(P["oe"]["all"], norm.vmin, norm.vmax), cmap=cmap, norm=norm, interpolation="nearest")
    s = max(0.6, n / 22)
    for q, i in enumerate(shown):
        col = SG_COLORS[z[i]] if P["N"] > 1 else FIG_GRAY
        ax.add_patch(Rectangle((-0.5 - 1.25 * s, q - 0.5), s, 1, color=col, clip_on=False, lw=0))
        ax.add_patch(Rectangle((q - 0.5, -0.5 - 1.25 * s), 1, s, color=col, clip_on=False, lw=0))
    for b in np.flatnonzero(np.diff(z[shown])) + 0.5:
        ax.axhline(b, color="white", lw=1.0)
        ax.axvline(b, color="white", lw=1.0)
    ax.set_xlim(-0.5 - 1.35 * s, n - 0.5)
    ax.set_ylim(n - 0.5, -0.5 - 1.35 * s)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_xticks([])
    if n <= 30:
        ax.yaxis.tick_right()
        ax.set_yticks(range(n))
        ax.set_yticklabels([_short_name(P["names"][i], 16) for i in shown], fontsize=5.2 if n <= 20 else 4.4)
        ax.tick_params(axis="y", length=0, pad=1.5)
    else:
        ax.set_yticks([])


def _panel_homoeologs(fig, ax, P):
    """v9 panel B: every shared old insertion (reciprocal closest relatives in two subgenomes) as a dot at its
    positions on the two chromosomes, chromosomes to scale in figure order; homoeologs show as collinear diagonals.
    Two subgenomes: SG1 across, SG2 down; more: every shown chromosome on both axes (same-subgenome blocks shaded,
    empty by construction). Dots on the collinear chain of a called pair in black, others gray; called pairs
    framed (dashed: called one way only)."""
    from matplotlib.patches import Rectangle
    shown, z, L, lib, ho = P["shown"], P["z"], P["L"], P["lib"], P["ph"]["homoeologs"]
    two = P["N"] == 2
    xs = [i for i in shown if z[i] == 0] if two else list(shown)
    ys = [i for i in shown if z[i] == 1] if two else list(shown)

    def offsets(seq):
        o, acc = {}, 0.0
        for i in seq:
            o[i], acc = acc, acc + L[i]
        return o, acc
    (ox, tx), (oy, ty) = offsets(xs), offsets(ys)
    ax.set_xlim(0, tx)
    ax.set_ylim(ty, 0)
    for i in xs:
        ax.add_patch(Rectangle((ox[i], -0.034 * ty), L[i], 0.025 * ty, color=SG_COLORS[z[i]], clip_on=False, lw=0))
        ax.axvline(ox[i], color="#dcdcdc", lw=0.3, zorder=1)
    for i in ys:
        ax.add_patch(Rectangle((-0.034 * tx, oy[i]), 0.025 * tx, L[i], color=SG_COLORS[z[i]], clip_on=False, lw=0))
        ax.axhline(oy[i], color="#dcdcdc", lw=0.3, zorder=1)
    if not two:
        for g in range(P["N"]):
            ig = [i for i in xs if z[i] == g]
            if ig:
                a0, a1 = ox[ig[0]], ox[ig[-1]] + L[ig[-1]]
                ax.add_patch(Rectangle((a0, a0), a1 - a0, a1 - a0, color="#f2f2f2", lw=0, zorder=0))
                if a0 > 0:
                    ax.axvline(a0, color="#8c8c8c", lw=0.6, zorder=1)
                    ax.axhline(a0, color="#8c8c8c", lw=0.6, zorder=1)
    ek, Pp, hb = np.asarray(P["ph"]["ek"]), ho["P"], ho.get("bio")
    mid = (lib.start + lib.end) / 2
    hot = ho["in_chain"] & ho["pair_called"]
    for a_, b_ in ((Pp[:, 0], Pp[:, 1]),) + (((Pp[:, 1], Pp[:, 0]),) if not two else ()):
        X = np.array([ox.get(int(c), np.nan) for c in ek[a_]]) + mid[a_]
        Y = np.array([oy.get(int(c), np.nan) for c in ek[b_]]) + mid[b_]
        ok = np.isfinite(X) & np.isfinite(Y)
        ex = hb["flags"] if hb else np.zeros(len(Pp), bool)
        for m_, col, sz, zo in ((ok & ~hot, "#b9b9b9", 0.5, 2), (ok & hot & ~ex, FIG_ACCENT, 0.8, 3),
                                (ok & hot & ex, FIG_EXCH, 2.2, 5)):
            ax.scatter(X[m_], Y[m_], s=sz, c=col, lw=0, zorder=zo, rasterized=True)
    for (i, h), c in ho["calls"].items():
        j = c["partner"]
        for a_, b_ in ((i, j), (j, i)):
            if j >= 0 and a_ in ox and b_ in oy:
                ax.add_patch(Rectangle((ox[a_], oy[b_]), L[a_], L[b_], fill=False, edgecolor=FIG_ACCENT, lw=0.6,
                                       ls="-" if c["reciprocal"] else (0, (2, 1.5)), zorder=4))
    for t in (hb["transloc"] if hb else []):                       # v9.1: block collinear with a non-homoeolog
        for a_, b_, xa, xb, ya, yb in ((t["i"], t["j"], t["start_i"], t["end_i"], t["start_j"], t["end_j"]),
                                       (t["j"], t["i"], t["start_j"], t["end_j"], t["start_i"], t["end_i"])):
            if a_ in ox and b_ in oy:
                pad = 0.006 * max(tx, ty)
                ax.add_patch(Rectangle((ox[a_] + xa - pad, oy[b_] + ya - pad), xb - xa + 2 * pad, yb - ya + 2 * pad,
                                       fill=False, edgecolor=FIG_EXCH, lw=0.6, ls=(0, (1, 1)), zorder=6))
    for sp in ax.spines.values():
        sp.set_visible(False)
    fs = 5.2 if max(len(xs), len(ys)) <= 12 else 4.4
    ax.set_xticks([])
    ax.set_yticks([])
    if len(ys) <= 30:
        ax.yaxis.tick_right()
        ax.set_yticks([oy[i] + L[i] / 2 for i in ys])
        ax.set_yticklabels([_short_name(P["names"][i], 12) for i in ys], fontsize=fs)
        ax.tick_params(axis="y", length=0, pad=1.5)
    if two and len(xs) <= 30:
        ax.set_xticks([ox[i] + L[i] / 2 for i in xs])
        ax.set_xticklabels([_short_name(P["names"][i], 12) for i in xs], fontsize=fs, rotation=90)
        ax.tick_params(axis="x", length=0, pad=1.5)


def _oe_bar(fig, x0, y0, w, norm, cmap, label):
    from matplotlib.cm import ScalarMappable
    from matplotlib.ticker import NullFormatter, NullLocator
    cax = _ax(fig, x0, y0, w, 0.07)
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal",
                      ticks=[0.25, 0.5, 1, 2, 4])
    cb.ax.set_xticklabels(["¼×", "½×", "1×", "2×", "4×"])
    cb.ax.xaxis.set_minor_locator(NullLocator())
    cb.ax.xaxis.set_minor_formatter(NullFormatter())
    cb.ax.tick_params(length=2, pad=1, labelsize=6)
    cb.outline.set_linewidth(0.4)
    if label:
        fig.text((x0 + w / 2) / fig.get_size_inches()[0], (y0 + 0.1) / fig.get_size_inches()[1], label, ha="center",
                 va="bottom", fontsize=5.6, color="#444444")


# ------------------------------------------------------------------ C: placement in half-libraries
def _pca(P):
    """Phased sequences placed by their lineage-sharing profiles: PCA of the rows of the log observed/expected matrix
    the partition maximises (masked links; the missing self-entry takes its column's mean). Each kept half-library's
    profiles are projected onto the same axes -> (positions [C, 2], half-library positions [h, C, 2] or None,
    variance share of the two axes)."""
    ph = P["ph"]

    def M_of(Wm):
        Lm = log_oe(np.asarray(Wm, float))
        np.fill_diagonal(Lm, np.nan)
        return np.where(np.isnan(Lm), np.nanmean(Lm, 0)[None], Lm)
    M = M_of(ph["W"])
    mu_ = M.mean(0)
    _, S, Vt = np.linalg.svd(M - mu_, full_matrices=False)
    k = min(2, Vt.shape[0])
    V = np.zeros((M.shape[1], 2))
    V[:, :k] = Vt[:k].T
    X = (M - mu_) @ V
    z = np.asarray(ph["z"])
    for j in range(2):                                 # the first subgenome on the left; axis 2 toward its largest load
        if (X[z == 0, j].mean() > 0) if j == 0 else (X[np.argmax(np.abs(X[:, j])), j] < 0):
            V[:, j] *= -1
            X[:, j] *= -1
    reps = np.array([(M_of(h) - mu_) @ V for h in ph["halves"]]) if ph.get("halves") else None
    share = np.zeros(2)
    share[:k] = S[:k] ** 2 / max(float((S ** 2).sum()), 1e-300)
    return X, reps, share


def _panel_pca(fig, ax, P):
    """Each phased sequence on the two leading axes of its lineage-sharing profile (dot) and in each half-library
    (cloud), colored by subgenome; names when they fit; the share of sequences kept in their group by >= 90% of the
    half-libraries."""
    z, names, ph = P["z"], P["names"], P["ph"]
    X, reps, share = _pca(P)
    col = lambda i: SG_COLORS[z[i]] if P["N"] > 1 else FIG_DARK
    if reps is not None:
        for i in range(len(names)):
            ax.scatter(reps[:, i, 0], reps[:, i, 1], s=1.5, color=col(i), alpha=0.25, lw=0, rasterized=True, zorder=2)
    for i in range(len(names)):
        ax.scatter([X[i, 0]], [X[i, 1]], s=14, color=col(i), edgecolors="white", linewidths=0.5, zorder=4)
    ax.axhline(0, color="#e2e2e2", lw=0.5, zorder=0)
    ax.axvline(0, color="#e2e2e2", lw=0.5, zorder=0)
    ax.set_xlabel(f"lineage-sharing axis 1 ({100 * share[0]:.0f}% of variance)", fontsize=6.2)
    ax.set_ylabel(f"axis 2 ({100 * share[1]:.0f}%)", fontsize=6.2)
    ax.tick_params(labelsize=5.6)
    ax_ = np.r_[X[:, 0], reps[:, :, 0].ravel()] if reps is not None else X[:, 0]
    ay_ = np.r_[X[:, 1], reps[:, :, 1].ravel()] if reps is not None else X[:, 1]
    xr, yr = (np.ptp(ax_) or 1.0), (np.ptp(ay_) or 1.0)
    ax.set_xlim(ax_.min() - 0.32 * xr, ax_.max() + 0.32 * xr)
    ax.set_ylim(ay_.min() - 0.24 * yr, ay_.max() + 0.25 * yr)
    if P["N"] > 1:
        gs = [g for g in range(P["N"]) if (z == g).any()]
        top = [((reps[:, z == g, 1].max() if reps is not None else X[z == g, 1].max()) + 0.06 * yr) for g in gs]
        per_in = np.diff(ax.get_xlim())[0] / (ax.get_position().width * fig.get_size_inches()[0])
        xs = (_repel([X[z == g, 0].mean() for g in gs], 0.3 * per_in) if P["N"] > 2 else
              [X[z == g, 0].mean() for g in gs])
        if P["N"] == 2:
            top = [ay_.max() + 0.1 * yr] * len(gs)
        for g, x_, y_ in zip(gs, xs, top):
            ax.text(x_, y_, group_label(g), color=SG_COLORS[g], fontsize=6.5, fontweight="bold", ha="center",
                    va="bottom")
    if len(names) <= 24:
        for g in range(max(P["N"], 1)):
            idx = np.flatnonzero(z == g)
            if not idx.size:
                continue
            side = 1 if X[idx, 0].mean() >= np.median(X[:, 0]) else -1
            ys = _repel(X[idx, 1], 0.085 * yr)
            x0 = (X[idx, 0].max() + 0.09 * xr) if side > 0 else (X[idx, 0].min() - 0.09 * xr)
            for i, yy in zip(idx, ys):
                ax.plot([X[i, 0], x0], [X[i, 1], yy], color="#c4c4c4", lw=0.3, zorder=1)
                ax.text(x0 + side * 0.01 * xr, yy, _short_name(names[i], 14), fontsize=5.0, va="center",
                        ha="left" if side > 0 else "right", color="#555555")
    sup = np.asarray(ph.get("support", []), float)
    ok = np.isfinite(sup)
    n_h = len(reps) if reps is not None else 0
    note = (f"clouds: {n_h} half-libraries" if n_h else f"half-libraries not kept (> {FIG_HALVES_MAX} sequences)")
    if ok.any() and P["N"] > 1:
        note += f"\n{int((sup[ok] >= SUPPORT_OK).sum())} of {int(ok.sum())} {P['unit']} kept their group in ≥ " \
                f"{100 * SUPPORT_OK:.0f}%"
    ax.text(0.5, 0.03, note, transform=ax.transAxes, fontsize=5.3, ha="center", va="bottom", color="#444444")


# ------------------------------------------------------------------ D: the two tests behind N
def _panel_tests(fig, axr, axs, P, max_rows=6):
    """The two tests every candidate split must pass, one row per split the search tried (kept first): reproduced by
    disjoint half-libraries (adjusted Rand >= MIN_REP, left) and at least MIN_REL as strong as the split above it
    (right), on the copies the decision used (all copies, or the age class that kept it)."""
    ph, z = P["ph"], P["z"]
    T = [t for t in ph.get("tree", []) if t.get("replicability_all") is not None or t["accepted"] or t["forced"]]
    T = sorted(T, key=lambda t: (not (t["accepted"] or t["forced"]), len(t["level"]), t["level"]))[:max_rows]
    for ax_ in (axr, axs):
        ax_.set_yticks([])
        ax_.spines["left"].set_visible(False)
    if not T:
        for ax_ in (axr, axs):
            ax_.axis("off")
        msg = ("groups given by the user (--assign): no split tested" if ph.get("source") == "user" else
               "no split could be tested")
        axr.text(0.0, 0.5, msg, transform=axr.transAxes, fontsize=6.3, ha="left", va="center", color="#666666")
        return
    zof = {nm: int(z[i]) for i, nm in enumerate(P["names"])}
    gsize = np.bincount(z, minlength=P["N"])
    reason = {"not reproduced": "not reproduced", "weaker than half its parent": "too weak",
              "present at every age": "same at every age", "no excess": "no excess",
              "pairs: not reproduced": "not reproduced"}
    for k, t in enumerate(T):
        kept = bool(t["accepted"] or t["forced"])
        c = FIG_INK if kept else "#9a9a9a"
        s = t["copies"] if (kept and t.get("copies")) else "all"
        r = t.get(f"replicability_{s}")
        st, pst = t.get(f"strength_{s}"), t.get(f"parent_strength_{s}")
        if r is not None:
            v = max(r, -0.05)
            axr.plot([v], [k], "o", ms=4 if kept else 3.2, color=c, clip_on=False, zorder=3)
            axr.text(v - 0.06 if v > 0.45 else v + 0.07, k, f"{r:.2f}", fontsize=5.2, va="center",
                     ha="right" if v > 0.45 else "left", color=c)
        if t["rule"] in ("pairs", "pairs: not reproduced"):
            axs.text(0.6, k, "homologue pairs\nkept apart", fontsize=5.0, va="center", ha="center", color=c)
        elif st is not None and pst:
            rel = st / pst if pst > 0 else np.nan
            if np.isfinite(rel):
                v = min(max(rel, 0.0), 1.2)
                axs.plot([v], [k], "o", ms=4 if kept else 3.2, color=c, clip_on=False, zorder=3)
                axs.text(v + 0.07 if v < 0.8 else v - 0.07, k, f"{100 * rel:.0f}%", fontsize=5.2, va="center",
                         ha="left" if v < 0.8 else "right", color=c)
        else:
            axs.text(0.6, k, "first split", fontsize=5.4, va="center", ha="center", color=c)
        lab = (f"all {P['unit']}: {len(t['part_a'])} vs {len(t['part_b'])}" if P["N"] == 1 else
               _split_words(t, zof, gsize)) + ("" if s == "all" else f" ({s} copies)")
        axr.text(-0.1, k, lab, transform=axr.get_yaxis_transform(), fontsize=5.8, va="center", ha="right", color=c,
                 fontweight="bold" if kept else "normal")
        word = ("forced (-n)" if t["forced"] else "kept") if kept else reason.get(t["rule"], "rejected")
        axs.text(1.04, k, word, transform=axs.get_yaxis_transform(), fontsize=5.8, va="center", ha="left", color=c,
                 fontweight="bold" if kept else "normal")
    for ax_, thr, lim, ttl in ((axr, MIN_REP, (-0.05, 1.05), "reproduced by\nindependent halves"),
                               (axs, MIN_REL, (0.0, 1.2), "as strong as the\nsplit above")):
        ax_.axvspan(thr, lim[1], color="#ececec", lw=0, zorder=0)
        ax_.axvline(thr, color="#888888", lw=0.5, ls=(0, (2, 2)), zorder=1)
        ax_.set_xlim(*lim)
        ax_.set_ylim(len(T) - 0.5, -0.6)
        ax_.tick_params(axis="x", labelsize=5.5, length=2)
        ax_.set_title(ttl, fontsize=5.8, pad=2, color="#444444")
    axr.set_xticks([0, MIN_REP, 1])
    axr.set_xticklabels(["0", f"{MIN_REP:g}", "1"])
    axs.set_xticks([0, MIN_REL, 1])
    axs.set_xticklabels(["0", "½" if MIN_REL == 0.5 else f"{MIN_REL:g}", "1×"])


def _tests_header(P):
    ph, N = P["ph"], P["N"]
    if ph.get("source") == "user":
        return f"{FIG_NUM.get(N, N).capitalize()} groups given by the user (--assign)"
    if any(t["forced"] for t in ph.get("tree", [])):
        return f"{FIG_NUM.get(N, N).capitalize()} groups set with -n: the splits that fail a test were forced"
    if N == 1:
        return "Why one group: no split passes both tests"
    kept = sum(bool(t["accepted"]) for t in ph.get("tree", []))
    if N == 2 and kept == 1:
        return "Why two subgenomes: only the first split passes both tests"
    return f"Why {FIG_NUM.get(N, N)} subgenomes: the kept splits pass both tests"


# ------------------------------------------------------------------ E-H: one time axis
def _time_axis(ax, P, labels=True, xlabel=True):
    """Years ago on a log axis, oldest left, present right; the rate behind it in the axis label."""
    from matplotlib.ticker import FixedFormatter, FixedLocator
    t0, t1 = P["trange"]
    ax.set_xscale("log")
    ax.set_xlim(t1, t0)
    tk = [10.0 ** e for e in range(0, 12) if t0 <= 10.0 ** e <= t1]
    lab = [_yr_label(v) for v in tk]
    if len(tk) > 5:
        lab = [x if (len(tk) - 1 - j) % 2 == 0 else "" for j, x in enumerate(lab)]
    ax.xaxis.set_major_locator(FixedLocator(tk))
    ax.xaxis.set_minor_locator(FixedLocator([m * v for v in [10.0 ** e for e in range(0, 12)] for m in range(2, 10)
                                             if t0 <= m * v <= t1]))
    ax.xaxis.set_major_formatter(FixedFormatter(lab if labels else [""] * len(tk)))
    ax.tick_params(axis="x", which="minor", length=1.5, width=0.4)
    if xlabel:
        ax.set_xlabel(f"years ago ({_rate(P['mu'])}" + ("" if P["mu_given"] else "; default rate") + ")", fontsize=6.2)


def _first_pair_events(P):
    """Reported split and merger of the first subgenome pair, in years: (split dict, merger dict) or Nones."""
    dt = P["dt"]
    if not dt or not dt["pairs"] or P["N"] != 2:
        return None, None
    r = dt["pairs"][0]
    yr = lambda k: k / (2 * P["mu"])
    sp = r["split"].get(r.get("split_reported")) if r.get("split_reported") else None
    mg = r["merger"] if r["merger"] and r["merger"]["datable"] else None
    s = dict(lo=yr(sp["lo"]), est=yr(sp["estimate"]), hi=yr(sp["hi"])) if sp else None
    m = dict(lo=P["trange"][0] if mg["bound"] else yr(mg["lo"]), hi=yr(mg["hi"]), bound=mg["bound"]) if mg else None
    return s, m


def _bands(ax, P):
    """The reported split and merger intervals (two subgenomes) as background bands of a time axis."""
    s, m = _first_pair_events(P)
    if s:
        ax.axvspan(s["lo"], s["hi"], color=FIG_SPLIT, lw=0, zorder=0)
    if m:
        ax.axvspan(m["lo"], m["hi"], color=FIG_MERGE, lw=0, zorder=0)


def _regime_strip(fig, ax, P, h=0.27, pad=0.02):
    """A strip above a time axis naming its periods between the bands (two subgenomes): one ancestral pool, the split,
    the progenitors apart, the merger, one nucleus; words shortened or dropped where a period is too narrow."""
    s, m = _first_pair_events(P)
    if not s and not m:
        return
    W, H = fig.get_size_inches()
    bb = ax.get_position()
    rax = fig.add_axes([bb.x0, bb.y1 + pad / H, bb.width, h / H])
    rax.set_xscale("log")
    rax.set_xlim(ax.get_xlim())
    rax.set_ylim(0, 1)
    rax.axis("off")
    t0, t1 = P["trange"]
    segs = []
    if s:
        segs += [(s["hi"], t1, "white", ["before the split:\none ancestral pool", "before\nsplit"]),
                 (s["lo"], s["hi"], FIG_SPLIT, ["split"]),
                 (m["hi"] if m else t0, s["lo"], "white", ["progenitors apart:\neach amplifies its own lineages",
                                                         "progenitors\napart"])]
    if m:
        segs.append((m["lo"], m["hi"], FIG_MERGE, ["merger\n(bound)", "merger"] if m["bound"] else ["merger"]))
        if not m["bound"]:
            segs.append((t0, m["lo"], "white", ["one nucleus", "after"]))
    span = math.log(t1 / t0)
    for a_, b_, col, words in segs:
        a_, b_ = max(a_, t0), min(b_, t1)
        if b_ <= a_:
            continue
        rax.axvspan(a_, b_, ymin=0, ymax=0.72, color=col, lw=0)
        room = math.log(b_ / a_) / span * bb.width * W
        for txt in words:
            if max(_text_w(fig, x_, 5.8) for x_ in txt.split("\n")) <= room - 0.04:
                rax.text(math.sqrt(a_ * b_), 0.36, txt, ha="center", va="center", fontsize=5.8, color="#333333")
                break
    rax.plot([t0, t1], [0, 0], color="#bbbbbb", lw=0.4)


def _panel_births(ax, axm, P):
    """Transposon births per subgenome: deconvolved copies per log-age bin (one ridge each, on one scale); the margin
    gives each subgenome's copies."""
    gs = [g for g in range(max(P["N"], 1)) if (g, "all") in P["hist"]]
    axm.axis("off")
    if not gs:
        ax.text(0.5, 0.5, "too few copies with LTR ages", transform=ax.transAxes, ha="center", va="center",
                fontsize=6.3, color="#666666")
        return
    vmax = max(P["hist"][(g, "all")]["y"].max() for g in gs) or 1.0
    off = 0.62
    for k, g in enumerate(gs):
        hh = P["hist"][(g, "all")]
        base = (len(gs) - 1 - k) * off
        col = SG_COLORS[g] if P["N"] > 1 else FIG_DARK
        ax.fill_between(hh["t"], base, base + hh["y"] / vmax, color=col, alpha=0.35, lw=0, zorder=2)
        ax.plot(hh["t"], base + hh["y"] / vmax, color=col, lw=0.9, zorder=3)
        tot, big = hh["n"], max(P["hist"][(g_, "all")]["n"] for g_ in gs)
        axm.text(0.0, base + 0.3, group_label(g) if P["N"] > 1 else "all", fontsize=6, fontweight="bold",
                 va="center", ha="left", color=col)
        axm.barh(base + 0.3, 0.85 * tot / big, left=0.62, height=0.22, color=col, lw=0)
        axm.text(0.62 + 0.85 * tot / big + 0.05, base + 0.3, f"{tot:,.0f}", fontsize=5.2, va="center", ha="left",
                 color=col)
    top = off * (len(gs) - 1) + 1.08
    for a_ in (ax, axm):
        a_.set_ylim(0, top)
    axm.set_xlim(0, 2.45)
    axm.text(0.62, top, "copies", fontsize=5.6, va="top", ha="left", color="#555555")
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)


def _panel_clades(ax, axm, P, k=6):
    """One row per main annotated LTR-RT clade (most copies; 'unknown' left out): its copies' births in each subgenome
    (deconvolved, each row scaled to its peak); the margin gives the clade's copies per subgenome."""
    tot = defaultdict(float)
    for r in P["bio"].get("clades", []):
        if str(r["clade"]).lower() not in ("unknown", "na", ""):
            tot[str(r["clade"])] += r["n"]
    gs = range(max(P["N"], 1))
    clades = [c for c, _ in sorted(tot.items(), key=lambda kv: -kv[1]) if any((g, c) in P["hist"] for g in gs)][:k]
    axm.axis("off")
    if not clades:
        ax.text(0.5, 0.5, f"no annotated clade with ≥ {MIN_HIST} copies in a subgenome", transform=ax.transAxes,
                ha="center", va="center", fontsize=6.3, color="#666666")
        ax.set_yticks([])
        return
    cnt = {(int(r["subgenome"]), str(r["clade"])): r["n"] for r in P["bio"].get("clades", [])}
    nmax = max([cnt.get((g, c), 0) for c in clades for g in gs] + [1])
    hb = 0.7 / max(P["N"], 1)
    for j, c in enumerate(clades):
        base = (len(clades) - 1 - j) * 1.0
        vmax = max((P["hist"][(g, c)]["y"].max() for g in gs if (g, c) in P["hist"]), default=1.0) or 1.0
        for g in gs:
            col = SG_COLORS[g] if P["N"] > 1 else FIG_DARK
            if (g, c) in P["hist"]:
                hh = P["hist"][(g, c)]
                ax.fill_between(hh["t"], base, base + hh["y"] / vmax * 0.86, color=col, alpha=0.2, lw=0)
                ax.plot(hh["t"], base + hh["y"] / vmax * 0.86, color=col, lw=0.8)
            n_ = cnt.get((g, c), 0)
            if P["N"] <= 2:
                yy = base + 0.08 + hb * (max(P["N"], 1) - 1 - g) + hb / 2
                axm.barh(yy, n_ / nmax, height=hb * 0.82, color=col, lw=0)
                axm.text(n_ / nmax + 0.04, yy, f"{n_:,}", fontsize=5.1, va="center", ha="left", color=col)
        if P["N"] > 2:                                 # one stacked bar (subgenome shares) and the clade's total
            ns = [cnt.get((g, c), 0) for g in gs]
            tot_c = sum(ns)
            left = 0.0
            for g, n_ in zip(gs, ns):
                axm.barh(base + 0.4, n_ / tot_c, left=left, height=0.38, color=SG_COLORS[g], lw=0)
                left += n_ / tot_c
            axm.text(1.04, base + 0.4, f"{tot_c:,}", fontsize=5.1, va="center", ha="left", color=FIG_DARK)
        ax.text(P["trange"][1] / 1.06, base + 0.12, c, fontsize=6, ha="left", va="bottom", fontstyle="italic")
    for a_ in (ax, axm):
        a_.set_ylim(-0.03, len(clades) + 0.02)
    axm.set_xlim(0, 1.9)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)


def _triangle(ax, OE, xs, norm, cmap):
    """Lineage sharing of every pair of columns as a rhombus below both (a stained-glass triangle)."""
    from matplotlib.collections import PolyCollection
    polys, vals = [], []
    for i in range(len(xs)):
        for j in range(i + 1, len(xs)):
            X, Y = (xs[i] + xs[j]) / 2, -(xs[j] - xs[i]) / 2
            polys.append([(X - 0.5, Y), (X, Y + 0.5), (X + 0.5, Y), (X, Y - 0.5)])
            vals.append(OE[i, j] if np.isfinite(OE[i, j]) else OE[j, i])
    ax.add_collection(PolyCollection(polys, array=np.ma.masked_invalid(np.clip(np.array(vals, float), norm.vmin,
                                                                               norm.vmax)),
                                     cmap=cmap, norm=norm, edgecolors="face", linewidths=0.15))


def _panel_age_classes(fig, ax, P, cmap, norm):
    """B's lineage sharing within young, middle-aged and old copies (LTR-K2P terciles), each triangle centred over
    its copies' age range (bracket beneath; leader lines where triangles had to move apart)."""
    from matplotlib.lines import Line2D
    edges = P["ph"].get("strata_edges") or []
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    if len(edges) != 2 or not all(s in P["oe"] for s in AGE_STRATA):
        ax.text(0.5, 0.5, "no copy-age classes", transform=ax.transAxes, ha="center", va="center", fontsize=6.3,
                color="#666666")
        return
    W, H = fig.get_size_inches()
    yr = lambda k: k / (2 * P["mu"])
    t0, t1 = P["trange"]
    rng = [("young", t0, yr(edges[0])), ("middle", yr(edges[0]), yr(edges[1])), ("old", yr(edges[1]), t1)]
    bb = ax.get_position()
    hin, win = bb.height * H, bb.width * W
    span = math.log(t1 / t0)
    n = len(P["shown"])
    th = hin * 0.6
    tw = min(2 * th * 1.15, win / 3.6)
    want = [(bb.x0 + math.log(t1 / math.sqrt(max(a_, t0) * min(b_, t1))) / span * bb.width) * W for _, a_, b_ in rng]
    lo_, hi_ = bb.x0 * W + tw / 2, (bb.x0 + bb.width) * W - tw / 2
    xs_ = np.clip(_repel(np.clip(want, lo_, hi_), tw + 0.08), lo_, hi_)
    tr = ax.get_xaxis_transform()
    z = P["z"]
    for (s_, a_, b_), x_in, xw in zip(rng, xs_, want):
        sub = fig.add_axes([(x_in - tw / 2) / W, bb.y0 + (hin * 0.38) / H, tw / W, th / H])
        _triangle(sub, P["oe"][s_], list(range(n)), norm, cmap)
        sub.set_xlim(-0.6, n - 0.4)
        sub.set_ylim(-(n - 1) / 2 - 0.55, 0.9)
        sub.axis("off")
        for g in range(P["N"]):
            gi = [k for k, i in enumerate(P["shown"]) if z[i] == g]
            if gi and P["N"] > 1:
                sub.plot([min(gi) - 0.45, max(gi) + 0.45], [0.55, 0.55], color=SG_COLORS[g], lw=1.4,
                         solid_capstyle="butt")
        if abs(x_in - xw) > 0.05:
            fig.add_artist(Line2D([x_in / W, xw / W], [(bb.y0 * H + hin * 0.38) / H, (bb.y0 * H + hin * 0.33) / H],
                                  color="#999999", lw=0.5))
        a2, b2 = max(a_, t0) * 1.05, min(b_, t1) / 1.05
        ax.plot([a2, b2], [0.3, 0.3], color="#777777", lw=0.6, transform=tr, solid_capstyle="butt")
        for e_ in (a2, b2):
            ax.plot([e_, e_], [0.25, 0.35], color="#777777", lw=0.6, transform=tr)
        fig.text(x_in / W, bb.y0 + hin * 0.21 / H, {"young": "young copies", "middle": "middle-aged",
                                                     "old": "old copies"}[s_],
                 ha="center", va="top", fontsize=5.6, color="#444444")


def _panel_profile(ax, P):
    """Share of closest-relative copy pairs that sit in two different subgenomes, by their divergence time: points
    (area: pairs) with 95% intervals, one series per subgenome pair (line in one subgenome's color, dots in the
    other's when there are several pairs)."""
    for pr in P["prof"]:
        x, y = pr["pair"]
        lc, fc = (FIG_INK, FIG_INK) if P["N"] == 2 else (SG_COLORS[x], SG_COLORS[y])
        ax.vlines(pr["t"], 100 * pr["lo"], 100 * pr["hi"], color=lc, lw=0.5, zorder=3)
        ax.plot(pr["t"], 100 * pr["f"], color=lc, lw=0.7, zorder=3)
        ax.scatter(pr["t"], 100 * pr["f"], s=2 + 10 * np.sqrt(pr["n"] / pr["n"].max()), color=fc, zorder=4,
                   edgecolors="white", linewidths=0.3)
    if P["N"] > 2 and P["prof"]:
        ys = _repel([100 * pr["f"][0] for pr in P["prof"]], 0.09 * max(100 * pr["hi"].max() for pr in P["prof"]))
        for pr, yy in zip(P["prof"], ys):
            x, y = pr["pair"]
            ax.text(pr["t"][0] / 1.25, yy, f"{group_label(x)}–{group_label(y)}", fontsize=5.2, va="center",
                    ha="left", color=SG_COLORS[x])
    ax.set_ylabel("closest relative in\nanother subgenome (%)", fontsize=6)
    top = max((100 * pr["hi"].max() for pr in P["prof"]), default=50.0)
    ax.set_ylim(0, max(5.0, 1.12 * top))


def _estimate_rows(P):
    """Rows of the estimate panel: (kind, estimator, estimate dict in K2P units or None, reported, pair prefix).
    Two subgenomes: every estimator; more: each pair's reported split and merger."""
    rows, pairs = [], list(P["dt"]["pairs"])
    if len(pairs) > 6:                                 # the six pairs with the oldest reported splits (all: dating.tsv)
        old = lambda r: (r["split"].get(r.get("split_reported")) or {}).get("estimate", -1.0) if r.get("split_reported") else -1.0
        pairs = sorted(pairs, key=old, reverse=True)[:6]
    for r in pairs:
        pre = f"{group_label(r['groups'][0])}–{group_label(r['groups'][1])} " if P["N"] > 2 else ""
        mg = r["merger"]
        if P["N"] == 2:
            for k in ("orthologs", "confinement", "profile"):
                rows.append(("split", k, r["split"].get(k), r.get("split_reported") == k, pre))
            for k, m_ in (("ltr_clock", r.get("merger_ltr")), ("genealogy_clock", r.get("merger_genealogy"))):
                rep = bool(mg) and mg.get("source") == k
                rows.append(("merger", k, mg if rep else m_, rep, pre))
        else:
            k = r.get("split_reported")
            rows.append(("split", k, r["split"].get(k) if k else None, True, pre))
            rows.append(("merger", mg.get("source") if mg else None, mg, True, pre))
    return rows


def _panel_estimates(fig, ax, axm, P, rows, fs=5.4):
    """Every row on the time axis (an interval with its estimate; a bound runs to the present): its value beside the
    interval on the free side, the estimator (dating.tsv's name) in the margin; bold: reported, gray: not."""
    t0, t1 = P["trange"]
    yr = lambda k: k / (2 * P["mu"])
    r_ = fig.canvas.get_renderer()
    bb = ax.get_window_extent(r_)
    axm.axis("off")
    for y, (kind, key, e, rep, pre) in enumerate(rows):
        c, wt = (FIG_INK, "bold") if rep else ("#8f8f8f", "normal")
        name = f"{pre}{kind}" + (f" · {FIG_EST[key]}" if key and not pre else "")
        axm.text(0.0, y, name, fontsize=5.2, va="center", ha="left", color=c, fontweight=wt)
        if not e or (kind == "merger" and not e["datable"]):
            ax.text(0.5, y, "not datable", transform=ax.get_yaxis_transform(), fontsize=fs, va="center",
                    ha="center", color="#a8a8a8", fontstyle="italic")
            continue
        bound = kind == "merger" and e["bound"]
        lo, hi = (t0 if bound else yr(e["lo"])), yr(e["hi"])
        ax.plot([lo, hi], [y, y], color=c, lw=1.5 if rep else 0.7, solid_capstyle="butt", zorder=3)
        if bound:
            ax.plot([hi], [y], marker="|", ms=5.5, mew=1.1 if rep else 0.7, color=c, zorder=4)
            txt = f"{kind} < {_ma(hi)}"
        else:
            ax.plot([yr(e["estimate"])], [y], "o", ms=3.8 if rep else 2.4, color=c, zorder=4,
                    mfc=c if kind == "split" else "white", mec=c, mew=0.9)
            txt = f"{kind} {_span(lo, yr(e['estimate']), hi)}"
        for side in (("right", "left") if kind == "split" else ("left", "right")):
            t = ax.text(lo / 1.12 if side == "right" else hi * 1.12, y, txt, fontsize=fs, va="center",
                        ha="left" if side == "right" else "right", color=c, fontweight=wt)
            ex = t.get_window_extent(r_)
            if ex.x0 >= bb.x0 - 1 and ex.x1 <= bb.x1 + 1:
                break
            t.remove()
        else:
            ax.text(0.99, y, txt, transform=ax.get_yaxis_transform(), fontsize=fs, va="center", ha="right", color=c,
                    fontweight=wt)
    for a_ in (ax, axm):
        a_.set_ylim(len(rows) - 0.4, -0.6)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)


def _panel_history(fig, P, x0, wt, wm, y_top, cmap, norm):
    """E-H in one column whose time axes all span [x0, x0 + wt] (inches); the margin of width wm at their right
    carries copy counts (E, F) and the estimators' names (H)."""
    dt = P["dt"]
    rows = _estimate_rows(P) if dt and dt["pairs"] else []
    hH = 0.8 if (dt and P["prof"]) else 0.5
    hF = 1.2 if len(rows) <= 6 else max(0.8, 1.2 - 0.09 * (len(rows) - 6))
    hR = 0.11 * len(rows) + 0.1
    gap, strip = 0.33, 0.27 if P["N"] == 2 and dt else 0.0
    y = y_top - strip - 0.04 - 0.62
    axE, axEm = _ax(fig, x0, y, wt, 0.62), _ax(fig, x0 + wt + 0.08, y, wm - 0.08, 0.62)
    _time_axis(axE, P, labels=False, xlabel=False)
    _bands(axE, P)
    _panel_births(axE, axEm, P)
    _regime_strip(fig, axE, P, h=strip or 0.27)
    _head(fig, axE, "E", "When each subgenome's LTR retrotransposons were born" if P["N"] > 1 else
          "When the LTR retrotransposons were born", FIG_XR, dy=strip + 0.07)
    y -= gap + hF
    axF, axFm = _ax(fig, x0, y, wt, hF), _ax(fig, x0 + wt + 0.08, y, wm - 0.08, hF)
    _time_axis(axF, P, labels=False, xlabel=False)
    _bands(axF, P)
    _panel_clades(axF, axFm, P)
    _head(fig, axF, "F", "Births of each main LTR-RT clade" + (", per subgenome" if P["N"] > 1 else ""), FIG_XR)
    y -= gap + 0.6
    axG = _ax(fig, x0, y, wt, 0.6)
    _time_axis(axG, P, labels=False, xlabel=False)
    _bands(axG, P)
    _panel_age_classes(fig, axG, P, cmap, norm)
    if P["dotplot"]:                                                    # v9: B shows homoeologs; the key sits here
        _head(fig, axG, "G", "Lineages shared by each pair of " + P["unit"] + " in young, middle-aged and old copies",
              FIG_XR)
        _oe_bar(fig, x0 + wt + 0.1, y + 0.12, wm - 0.14, norm, cmap, "shared vs random")
    else:
        _head(fig, axG, "G", "B's lineage sharing within young, middle-aged and old copies", FIG_XR)
    y -= gap + hH
    axH = _ax(fig, x0, y, wt, hH)
    _time_axis(axH, P, labels=not rows, xlabel=not rows)
    _bands(axH, P)
    if dt and P["prof"]:
        _panel_profile(axH, P)
        _head(fig, axH, "H", "Closest relatives in another subgenome over time, and every date estimate", FIG_XR)
    else:
        axH.set_yticks([])
        axH.spines["left"].set_visible(False)
        msg = ("one group: no split or merger to date" if P["N"] == 1 else
               "dating skipped" if not dt else "too few copy pairs to date")
        axH.text(0.5, 0.5, msg, transform=axH.transAxes, ha="center", va="center", fontsize=6.3, color="#666666")
        _head(fig, axH, "H", "Closest relatives in another subgenome over time", FIG_XR)
    if rows:
        y -= 0.05 + hR
        axR, axRm = _ax(fig, x0, y, wt, hR), _ax(fig, x0 + wt + 0.08, y, wm - 0.08, hR)
        _time_axis(axR, P)
        _bands(axR, P)
        _panel_estimates(fig, axR, axRm, P, rows)


# ------------------------------------------------------------------ the page
def fig_report(plt, a, lib, ph, dt, bio, pa):
    """report.pdf's one page (see the section comment); returns the figure."""
    from matplotlib.colors import LogNorm
    import matplotlib as mpl
    P = _fig_data(a, lib, ph, dt, bio, pa)
    W, H = FIG_W, FIG_H
    fig = plt.figure(figsize=(W, H))
    cmap, norm = mpl.colormaps["magma"], LogNorm(0.25, 4)
    N = P["N"]
    lab_n = ("groups given by the user" if ph.get("source") == "user" else "subgenomes") if N > 1 else ""
    x = 0.15
    t = fig.text(x / W, (H - 0.2) / H, P["label"], fontsize=9, fontstyle="italic", fontweight="bold", va="baseline")
    x += t.get_window_extent(fig.canvas.get_renderer()).width / fig.dpi + 0.12
    unph = len(lib.seqnames) - len(P["names"])
    verdict = str(ph.get("verdict", "one group"))
    short, _, detail = verdict.partition(":") if N == 1 else (verdict, "", "")
    head = (f"·  {FIG_NUM.get(N, N)} {lab_n}" if N > 1 else f"·  one group ({short.strip()})") + \
        f"  ·  {len(P['names'])} {P['unit']} phased" + (f" ({unph} too short or too poor in LTR-RTs)" if unph else "") \
        + f"  ·  {lib.n:,} LTR retrotransposons"
    fig.text(x / W, (H - 0.2) / H, head, fontsize=7.5, va="baseline", color=FIG_INK)
    if detail.strip():
        fig.text(0.15 / W, (H - 0.33) / H, detail.strip()[:1].upper() + detail.strip()[1:], fontsize=5.8,
                 va="baseline", color="#666666")
    # A: karyotype across the top
    XL, XR = 0.55, W - 0.15
    axA = _ax(fig, XL, 6.3, XR - XL, 8.48 - _name_layout(fig, P)[2] - 6.3)     # names fit under the group heads
    xpos = _panel_karyotype(fig, axA, P)
    _group_heads(fig, axA, P, xpos, 8.56)
    _head(fig, axA, "A", ("Chromosomes painted by the group their LTR retrotransposon lineages point to" if
                          ph.get("source") == "user" else
                          "Chromosomes painted by the progenitor their LTR retrotransposon lineages come from")
          if P["tiles"] else "Chromosomes and their LTR retrotransposons (nothing to paint)", FIG_XL,
          dy=H - 0.47 - axA.get_position().y1 * H)
    items = ([("box", SG_COLORS[g], f"{group_label(g)} lineages") for g in range(N)] if N <= 4 else
             [("box", SG_COLORS[0], "subgenome lineages (group colors)")]) if P["tiles"] else []
    if P["tiles"]:
        items += [("hatch", None, "too few copies to tell")]
    if P["calls"]:
        items += [("tri", SG_COLORS[0], "called exchange (donor color)")]
    if P["dotplot"] and _homoeolog_numbers(ph):
        items += [("num", None, "homoeolog set")]
    items += [("dens", None, "LTR-RT density")]
    room = W - 0.1 - (XL - 0.13)                                       # v9: one more item; tighten when it overflows
    need = lambda gap: sum(0.17 + _text_w(fig, lab, 6.2) + gap for _, _, lab in items)
    _key_row(fig, XL - 0.13, H - 0.61, items, gap=0.22 if need(0.22) <= room else max(0.06, 0.22 - (need(0.22) - room)
                                                                                         / max(len(items) - 1, 1)))
    # left column: B, C, D
    axB = _ax(fig, 0.72, 3.88, 1.62, 1.62)
    if P["dotplot"]:                                                    # v9: homoeolog evidence
        _panel_homoeologs(fig, axB, P)
        _head(fig, axB, "B", f"Old insertions shared by {P['unit']} of two subgenomes", FIG_XL, dy=0.17)
        hb = ph["homoeologs"].get("bio")
        kb = [("dot", FIG_ACCENT, "collinear, homoeologs"), ("dot", "#b9b9b9", "other"), ("frame", FIG_ACCENT, "pair called")]
        if hb and hb["segments"]:
            kb.append(("dot", FIG_EXCH, "near-identical (exchanged)"))
        if hb and hb["transloc"]:
            kb.append(("frame", FIG_EXCH, "translocated"))
        _key_row(fig, 0.62, 3.53, kb[:3], size=5.6, gap=0.1)
        if len(kb) > 3:
            _key_row(fig, 0.62, 3.40, kb[3:], size=5.6, gap=0.1)
    else:
        _panel_sharing(fig, axB, P, cmap, norm)
        _head(fig, axB, "B", f"Transposon lineages shared by each pair of {P['unit']}", FIG_XL, dy=0.17)
        _oe_bar(fig, 0.75, 3.55, 1.5, norm, cmap, "lineages shared vs random")
    axC = _ax(fig, 0.85, 1.78, 1.85, 1.3)
    _panel_pca(fig, axC, P)
    _head(fig, axC, "C", f"Placement of each {P['unit'][:-1]} in independent halves of the library", FIG_XL,
          dy=0.12)
    axDr, axDs = _ax(fig, 1.45, 0.42, 0.72, 0.58), _ax(fig, 2.3, 0.42, 0.55, 0.58)
    _panel_tests(fig, axDr, axDs, P)
    _head(fig, axDr, "D", _tests_header(P), FIG_XL, dy=0.36)
    # right column: E-H on one time axis
    _panel_history(fig, P, 3.75, 2.25, 0.85, 5.8, cmap, norm)
    import socket
    fig.text((W - 0.1) / W, 0.06 / H, f"{Path(__file__).name} {VERSION} · {socket.gethostname()}:"
             f"{Path(a.outdir).resolve() if a.outdir else ''} · {time.strftime('%Y-%m-%d')}", fontsize=4.5,
             ha="right", va="bottom", color="#999999")
    return fig


def write_pdf(a, lib, ph, dt, bio, pa):
    """report.pdf: the one-page figure (v8_mod; v8 added two legend pages), the one PDF of a run."""
    plt, family = _fig_style()
    fig = fig_report(plt, a, lib, ph, dt, bio, pa)
    out = Path(a.outdir) / "report.pdf"
    fig.savefig(out, metadata={"Title": f"{a.label or Path(a.ltr_fasta[0]).name}: subgenomes from LTR retrotransposons",
                               "Creator": f"{Path(__file__).name} {VERSION}"})
    plt.close(fig)
    return out


def write_elements(a, lib, ph, dt, bio, pa):
    """elements.tsv: every LTR-RT with its sequence's subgenome (NA when unphased), clade, family, LTR divergence and
    age, the calibrated ancestry posterior per subgenome for informative copies, a call by the painting's evidence rule
    (foreign: its lineages favor another subgenome >= BF_SUPPORT-fold over its sequence's own; own: the reverse; weak:
    neither; NA for one subgenome or uninformative copies) with the favored other subgenome, and the probability that
    it was inserted after its subgenome's dated merger (NA without a dated merger)."""
    N, ek, z = ph["N"], np.asarray(ph["ek"]), np.asarray(ph["z"])
    sg = np.where(ek >= 0, z[np.maximum(ek, 0)], -1)
    mu = dt["mu"] if dt else (a.mu or MU_DEFAULT)
    clades = clade_names(lib)
    after = bio.get("after", {})
    if N > 1:
        anc = pa["ancestry"] if pa else (dt["ancestry"] if dt else None)
        if anc is None:
            anc = copy_ancestry(ph["clades"], ek, z, N, np.random.default_rng([a.seed, 7]))
        inf, foreign, home, gbest = _favors(lib, ph, anc)
        P = np.exp(anc[0])
    lines = ["element\tsequence\tstart\tend\tstrand\tclade\tfamily\tk2p\tage_years\tsubgenome\tinformative\t"
             + "\t".join(f"P_{group_label(g)}" for g in range(max(N, 1))) + "\tcall\tfavored\tp_after_merger"]
    for i in range(lib.n):
        k = lib.k2p[i]
        if N > 1 and sg[i] >= 0 and inf[i]:
            probs = "\t".join(f"{v:.4f}" for v in P[i])
            call = "foreign" if foreign[i] else "own" if home[i] else "weak"
            fav = group_label(int(gbest[i])) if foreign[i] else "NA"
            info = "yes"
        else:
            probs, call, fav = "\t".join("NA" for _ in range(max(N, 1))), "NA", "NA"
            info = "no" if N > 1 and sg[i] >= 0 else "NA"
        pa_ = after[sg[i]][i] if sg[i] in after else float("nan")
        lines.append(f"{lib.ids[i]}\t{lib.seqnames[lib.seq[i]]}\t{lib.start[i]}\t{lib.end[i]}\t{lib.strand[i] or '.'}\t"
                     f"{clades[i]}\t{lib.family[i] or 'NA'}\t{_f(k, '{:.5f}')}\t"
                     f"{_f(k / (2 * mu) if np.isfinite(k) else None, '{:.0f}')}\t"
                     f"{group_label(int(sg[i])) if sg[i] >= 0 else 'NA'}\t{info}\t{probs}\t{call}\t{fav}\t"
                     f"{_f(pa_, '{:.3f}')}")
    (Path(a.outdir) / "elements.tsv").write_text("\n".join(lines) + "\n")


def _f(x, fmt="{:.4g}"):
    return "NA" if x is None or (isinstance(x, float) and not math.isfinite(x)) else fmt.format(x)


def _jsonable(x):
    if isinstance(x, dict):
        return {k: _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (float, np.floating)):
        return float(x) if math.isfinite(x) else None
    return x


def write_phasing(a, lib, cl, ph, timings, dt=None, bio=None, pa=None):
    od = Path(a.outdir)
    n_el_all = np.bincount(lib.seq, minlength=len(lib.seqnames))
    idx = {s: i for i, s in enumerate(ph["names"])}
    partner, psup, pex = ph["partners"]
    with open(od / "chromosomes.tsv", "w") as f:
        f.write("sequence\tlength\tn_ltr\tn_informative\tsubgenome\tstatus\tsupport\toe_own\toe_other\tpvalue\t"
                "partner\tpartner_support\tpartner_excess\n")
        for s, L, n in zip(lib.seqnames, lib.seq_len, n_el_all):
            i = idx.get(s)
            if i is None:
                f.write(f"{s}\t{L}\t{n}\tNA\tNA\tunphased" + "\tNA" * 7 + "\n")
                continue
            if ph["source"] == "user":
                st = "user"
            elif ph["N"] == 1:
                st = "no_structure" if ph["verdict"] == "no structure detected" else "cannot_tell"
            elif ph["n_inf"][i] < MIN_ELEMENTS:
                st = "insufficient"
            else:
                st = "assigned" if ph["support"][i] >= SUPPORT_OK else "ambiguous"
            pj = partner[i]
            f.write(f"{s}\t{L}\t{n}\t{ph['n_inf'][i]}\t{group_label(ph['z'][i]) if ph['N'] > 1 else 'SG1'}\t{st}\t"
                    f"{_f(ph['support'][i], '{:.3f}')}\t{_f(ph['own_oe'][i])}\t{_f(ph['oth_oe'][i])}\t"
                    f"{_f(ph['pval'][i], '{:.3g}')}\t{ph['names'][pj] if pj >= 0 else 'NA'}\t{_f(psup[i], '{:.2f}')}\t"
                    f"{_f(pex[i], '{:.2f}')}\n")
    with open(od / "structure.tsv", "w") as f:
        src = [k[len("replicability_"):] for k in (ph["tree"][0] if ph["tree"] else {})
               if k.startswith("replicability_") and k != "replicability_se"]
        f.write("level\tn\tn_a\tn_b\tpart_a\tpart_b\tcopies\treplicability\treplicability_se\tn_pairs\tborderline\t"
                "strength\tparent_strength\t" + "\t".join(f"replicability_{x}\tstrength_{x}\tagrees_{x}\t"
                                                            f"parent_strength_{x}" for x in src)
                + "\trule\taccepted\tforced\n")
        yn = lambda v: "NA" if v is None else ("yes" if v else "no")
        for t in ph["tree"]:
            f.write(f"{t['level']}\t{t['n']}\t{len(t['part_a'])}\t{len(t['part_b'])}\t{','.join(t['part_a'])}\t"
                    f"{','.join(t['part_b'])}\t{t['copies'] or 'NA'}\t{_f(t['replicability'], '{:.3f}')}\t"
                    f"{_f(t['replicability_se'], '{:.3f}')}\t{t['n_pairs']}\t{yn(t['borderline'])}\t"
                    f"{_f(t['strength'])}\t{_f(t['parent_strength'])}\t"
                    + "\t".join(f"{_f(t[f'replicability_{x}'], '{:.3f}')}\t{_f(t[f'strength_{x}'])}\t"
                                 f"{yn(t[f'agrees_{x}'])}\t{_f(t[f'parent_strength_{x}'])}"
                                 for x in src)
                    + f"\t{t['rule']}\t{yn(t['accepted'])}\t{yn(t['forced'])}\n")
    with open(od / "pairing_map.tsv", "w") as f:
        f.write("sequence_a\tsequence_b\tn_candidates\tchain\tnull_mean\tz\tsignificant\tspan_a\tspan_b\t"
                "shared_share\tanchor_k2p_median\tcopy_k2p_median\tpaired\n")
        for r in sorted(ph["pairing"], key=lambda r: -r["z"]):
            f.write(f"{ph['names'][r['i']]}\t{ph['names'][r['j']]}\t{r['n_candidates']}\t{r['chain']}\t"
                    f"{r['null_mean']:.2f}\t{r['z']:.2f}\t{'yes' if r['significant'] else 'no'}\t"
                    f"{r['span_a']:.3f}\t{r['span_b']:.3f}\t{_f(r['shared'], '{:.2f}')}\t"
                    f"{_f(r['d_median'], '{:.5f}')}\t{_f(r['k_median'], '{:.5f}')}\t{'yes' if r['paired'] else 'no'}\n")
    Wr = ph["W_raw"]
    E = expected(Wr)
    with open(od / "links_oe.tsv", "w") as f:                    # observed / expected links (v5), NA on the diagonal
        f.write("sequence\t" + "\t".join(ph["names"]) + "\n")
        for i, s in enumerate(ph["names"]):
            f.write(s + "\t" + "\t".join("NA" if i == j or E[i, j] <= 0 else f"{Wr[i, j] / E[i, j]:.4g}"
                                         for j in range(len(ph["names"]))) + "\n")
    sup = ph["support"][np.isfinite(ph["support"])]
    summ = dict(version=VERSION, ltr_fasta=[str(p) for p in a.ltr_fasta], annotation=[str(p) for p in a.annotation],
                k=K, scale=cl.scale, n_ltr=lib.n, n_sequences=len(lib.seqnames), n_phased=ph["C"],
                n_kmers=cl.n_occ, n_multi_kmers=cl.n_multi, n_clades=cl.n, n_clades_phased=ph["n_clades"],
                n_subgenomes=ph["N"], n_subgenomes_requested=ph["N_requested"], n_subgenomes_source=ph["source"],
                user_labels=ph.get("user_labels"), verdict=ph["verdict"],
                min_support=float(sup.min()) if sup.size else None, n_half_pairs=ph["n_pairs"],
                subgenome_sizes={group_label(g): int((ph["z"] == g).sum()) for g in range(ph["N"])},
                pairing=dict(n_tested=len(ph["pairing"]), n_significant=sum(r["significant"] for r in ph["pairing"]),
                             n_paired=sum(r["paired"] for r in ph["pairing"]),
                             cannot_link_used=[[ph["names"][x], ph["names"][y]] for x, y in ph["cl_pairs"]],
                             paired_groups=[[ph["names"][x] for x in g] for g in ph["paired_groups"]]),
                structure_tree=ph["tree"], strata_k2p_edges=ph["strata_edges"], strata_copies=ph["strata_n"],
                homoeologs=homoeolog_summary(ph), timings=timings,
                dating=dating_summary(dt) if dt else None, biology=biology_summary(bio) if bio else None,
                painting=painting_summary(pa) if pa else None)
    (od / "summary.json").write_text(json.dumps(_jsonable(summ), indent=1) + "\n")


def homoeolog_summary(ph):
    """summary.json 'homoeologs': evidence used and, per subgenome pair, sequences paired in each direction."""
    ho = ph.get("homoeologs")
    if not ho:
        return None
    z, N, names = np.asarray(ph["z"]), ph["N"], ph["names"]
    calls, pairs = ho["calls"], []
    for g in range(N):
        for h in range(g + 1, N):
            ig, ih = np.flatnonzero(z == g), np.flatnonzero(z == h)
            stat = lambda rows, k: {s_: sum(calls[(i, k)]["status"] == s_ for i in rows)
                                    for s_ in ("homoeolog", "unresolved", "none")}
            pairs.append(dict(subgenomes=[group_label(g), group_label(h)], sequences=[int(ig.size), int(ih.size)],
                              status=[stat(ig, h), stat(ih, g)],
                              reciprocal=int(sum(bool(calls[(i, h)]["reciprocal"]) for i in ig))))
    hb, bio = ho.get("bio"), None
    if hb:
        yr = lambda k: k / (2 * hb["mu"]) if k is not None and np.isfinite(k) else None
        nf = lambda x: float(x) if x is not None and np.isfinite(x) else None
        bio = dict(exchanges=[dict(a=names[s_["i"]], b=names[s_["j"]], start_a=s_["start_i"], end_a=s_["end_i"],
                                   start_b=s_["start_j"], end_b=s_["end_j"], pairs=s_["n"], years=yr(s_["k2p"]),
                                   timing=s_["timing"], p=s_["p"]) for s_ in hb["segments"]],
                   half_life=[dict(lost_in=group_label(e["h"]), copies_from=group_label(e["g"]), k2p=nf(e["est"]),
                                   low=nf(e["lo"]), high=nf(e["hi"]), years=yr(e["est"]), qc_pass=e["ok"])
                              for e in hb["halflife"]],
                   lineages=[dict(subgenomes=[group_label(s_["a"]), group_label(s_["b"])], ratio=nf(s_["ratio"]),
                                  low=nf(s_["lo"]), high=nf(s_["hi"]), differ=s_["differ"]) for s_ in hb["asym"]],
                   translocated_blocks=len(hb["transloc"]),
                   several_homoeologs=[dict(sequence=names[m["i"]], partners=[names[p_] for p_ in m["partners"]])
                                       for m in hb["many"]],
                   collinear_in_reverse=[[names[i], names[j]] for i, j in hb["reversed"]])
    return dict(old_k2p_edge=ho["old_edge"], old_copies=ho["n_old"], anchors=int(len(ho["P"])), pairs=pairs,
                sets=[{group_label(g): [names[x] for x in s_["members"] if z[x] == g] for g in range(N)}
                      | dict(complete=s_["complete"]) for s_ in ho["sets"]], biology=bio)


def write_homoeologs(a, lib, ph):
    """homoeologs.tsv (each phased sequence x other subgenome: status, the called homoeolog, the best candidate and its
    evidence: anchors, observed/expected, longest collinear chain with its z against NULL_PERMS permutations and its
    p (chain_test: Monte Carlo, resolution alpha/2, or the exact bound for <= CHAIN_W anchors), orientation, spans), homoeolog_sets.tsv (connected sets of calls, one column per subgenome) and homoeolog_anchors.tsv (every
    shared old insertion: both copies, whether on its pair's longest collinear chain, whether the pair is called).
    Nothing is written for one group."""
    ho = ph.get("homoeologs")
    if not ho:
        return
    od, z, N, names, ek = Path(a.outdir), np.asarray(ph["z"]), ph["N"], ph["names"], np.asarray(ph["ek"])
    order = sorted(range(len(names)), key=lambda i: natural_key(names[i]))
    with open(od / "homoeologs.tsv", "w") as f:
        f.write("sequence\tsubgenome\tother_subgenome\tstatus\thomoeolog\treciprocal\tbest_candidate\tanchors\t"
                "anchors_oe\tchain\tchain_z\tchain_p\torientation\tspan\tcandidate_span\n")
        for i in order:
            for h in range(N):
                if h == z[i]:
                    continue
                c = ho["calls"][(i, h)]
                j = c["candidate"]
                r = ho["pair"].get((i, j)) if j >= 0 else None
                yn = "NA" if c["reciprocal"] is None else ("yes" if c["reciprocal"] else "no")
                ev = ("\t".join(["NA"] * 9) if r is None else
                      f"{names[j]}\t{r['n']}\t{ho['M'][i, j] / ho['E'][i, j]:.3g}\t{r['chain']}\t{r['z']:.3g}\t"
                      f"{c['p']:.3g}\t{r['orientation']}\t{r['span']:.3f}\t{r['partner_span']:.3f}")
                f.write(f"{names[i]}\t{group_label(int(z[i]))}\t{group_label(h)}\t{c['status']}\t"
                        f"{names[c['partner']] if c['partner'] >= 0 else 'NA'}\t{yn}\t{ev}\n")
    with open(od / "homoeolog_sets.tsv", "w") as f:
        f.write("set\t" + "\t".join(group_label(g) for g in range(N)) + "\tcomplete\n")
        for k, s_ in enumerate(ho["sets"]):
            cols = [",".join(sorted((names[x] for x in s_["members"] if z[x] == g), key=natural_key)) or "-"
                    for g in range(N)]
            f.write(f"H{k + 1}\t" + "\t".join(cols) + f"\t{'yes' if s_['complete'] else 'no'}\n")
    hb = ho.get("bio")
    with open(od / "homoeolog_anchors.tsv", "w") as f:
        f.write("sequence_a\tstart_a\tend_a\tsubgenome_a\tsequence_b\tstart_b\tend_b\tsubgenome_b\tk2p_a\tk2p_b\t"
                "in_chain\tpair_called" + ("\taligned_sites\tdivergence\tnear_identical" if hb else "") + "\n")
        P = ho["P"]
        for q in np.lexsort((lib.start[P[:, 0]], ek[P[:, 0]])):
            x, y = P[q]
            f.write(f"{names[ek[x]]}\t{lib.start[x]}\t{lib.end[x]}\t{group_label(int(z[ek[x]]))}\t{names[ek[y]]}\t"
                    f"{lib.start[y]}\t{lib.end[y]}\t{group_label(int(z[ek[y]]))}\t{_f(lib.k2p[x], '{:.5f}')}\t"
                    f"{_f(lib.k2p[y], '{:.5f}')}\t{'yes' if ho['in_chain'][q] else 'no'}\t"
                    f"{'yes' if ho['pair_called'][q] else 'no'}" +
                    (f"\t{hb['sites'][q]}\t{_f(hb['k2p'][q], '{:.5f}')}\t{'yes' if hb['flags'][q] else 'no'}"
                     if hb else "") + "\n")
    if hb:
        write_homoeolog_biology(a, ph, hb)


def write_homoeolog_biology(a, ph, hb):
    """v9.1 tables: homoeolog_exchanges.tsv (near-identical stretches on called homoeologs), homoeolog_halflife.tsv
    (half-life of intact LTR-RTs per direction), homoeolog_rearrangements.tsv (translocated blocks; sequences with
    several homoeologs; homoeologs collinear in reverse) and findings.txt (plain-language results)."""
    od, z, names = Path(a.outdir), np.asarray(ph["z"]), ph["names"]
    yr = lambda k: k / (2 * hb["mu"]) if k is not None and np.isfinite(k) else None
    lab = lambda i: f"{names[i]}\t{group_label(int(z[i]))}"
    with open(od / "homoeolog_exchanges.tsv", "w") as f:
        f.write("stretch\tsequence_a\tsubgenome_a\tstart_a\tend_a\tsequence_b\tsubgenome_b\tstart_b\tend_b\t"
                "near_identical_pairs\tchain_pairs\tdivergence\tyears\ttiming\tp\tv8_exchanges\n")
        for k, s_ in enumerate(hb["segments"]):
            f.write(f"X{k + 1}\t{lab(s_['i'])}\t{s_['start_i']}\t{s_['end_i']}\t{lab(s_['j'])}\t{s_['start_j']}\t"
                    f"{s_['end_j']}\t{s_['n']}\t{s_['n_chain']}\t{s_['k2p']:.5f}\t{_f(yr(s_['k2p']))}\t{s_['timing']}\t"
                    f"{s_['p']:.3g}\t{','.join(s_['v8']) or '-'}\n")
    with open(od / "homoeolog_halflife.tsv", "w") as f:
        f.write("lost_in\tcopies_from\tsplit_divergence\tinherited_pairs\tretention_1.5x\tretention_2x\told_copies\t"
                "half_life_k2p\tlow\thigh\thalf_life_years\tyears_low\tyears_high\tcopies_per_mb_lost_in\tqc\n")
        for e in hb["halflife"]:
            bad = ",".join(k for k, v in e["qc"].items() if not v)
            n_use = hb["pairs"][(min(e["g"], e["h"]), max(e["g"], e["h"]))]["use"].size
            f.write(f"{group_label(e['h'])}\t{group_label(e['g'])}\t{e['d_s']:.5f}\t{n_use}\t{e['r15']:.4f}\t"
                    f"{e['r2']:.4f}\t{e['n_old']}\t{_f(e['est'])}\t{_f(e['lo'])}\t{_f(e['hi'])}\t{_f(yr(e['est']))}\t"
                    f"{_f(yr(e['lo']))}\t{_f(yr(e['hi']))}\t{hb['density'][e['h']]:.2f}\t{'pass' if e['ok'] else bad}\n")
    with open(od / "homoeolog_rearrangements.tsv", "w") as f:
        f.write("type\tsequence\tsubgenome\tstart\tend\tpartner\tpartner_subgenome\tpartner_start\tpartner_end\t"
                "collinear_pairs\tpairs_tested\tp\n")
        for t in hb["transloc"]:
            f.write(f"translocated_block\t{lab(t['i'])}\t{t['start_i']}\t{t['end_i']}\t{lab(t['j'])}\t{t['start_j']}\t"
                    f"{t['end_j']}\t{t['n']}\t{t['n_tested']}\t{t['p']:.3g}\n")
        for m in hb["many"]:
            f.write(f"several_homoeologs\t{lab(m['i'])}\tNA\tNA\t{','.join(names[p_] for p_ in m['partners'])}\t"
                    f"{group_label(m['h'])}\tNA\tNA\tNA\tNA\tNA\n")
        for i, j in hb["reversed"]:
            f.write(f"collinear_in_reverse\t{lab(i)}\tNA\tNA\t{lab(j)}\tNA\tNA\tNA\tNA\tNA\n")
    lines = homoeolog_findings(ph, hb)
    with open(od / "findings.txt", "w") as f:
        f.write("\n".join(lines) + "\n")
    for l in lines:
        log(l)


def homoeolog_findings(ph, hb):
    """Plain-language statements of the v9.1 results, each with its numbers and how far the benchmark trusts it."""
    z, names, N = np.asarray(ph["z"]), ph["names"], ph["N"]
    ho = ph["homoeologs"]
    yrs = lambda k: _ma(k / (2 * hb["mu"]))
    nm = lambda i: f"{names[i]} ({group_label(int(z[i]))})"
    out = ["What the homoeologous chromosome pairs say about this genome (ltrharvest_phasing9 v9.1). Each statement "
           "passed a benchmark bar set in advance (see the module notes, V9.1).", ""]
    for g in range(N):
        for h in range(g + 1, N):
            ig, ih = np.flatnonzero(z == g), np.flatnonzero(z == h)
            cg = sum(ho["calls"][(i, h)]["status"] == "homoeolog" for i in ig)
            ch = sum(ho["calls"][(j, g)]["status"] == "homoeolog" for j in ih)
            rest = cg < ig.size or ch < ih.size
            out.append(f"Homoeologs {group_label(g)}-{group_label(h)}: {cg} of {ig.size} {group_label(g)} and {ch} of "
                       f"{ih.size} {group_label(h)} sequences paired" +
                       ("; the others keep too few inherited insertions in collinear order (transposon turnover "
                        "erases them as subgenomes age)." if rest else "."))
    if hb["reversed"]:
        out.append(f"Orientation: {len(hb['reversed'])} homoeolog pair(s) run in opposite directions in this assembly "
                   f"({', '.join(f'{names[i]}-{names[j]}' for i, j in hb['reversed'][:8])}"
                   f"{', ...' if len(hb['reversed']) > 8 else ''}); one member of each is assembled reverse-complemented "
                   f"relative to the other. Benchmark: 252 of 256 called pairs agree with gene order.")
    else:
        out.append("Orientation: every called homoeolog pair runs in the same direction in this assembly.")
    sg_ = hb["segments"]
    if sg_:
        tim = defaultdict(int)
        for s_ in sg_:
            tim[s_["timing"]] += 1
        desc = "; ".join(f"{nm(s_['i'])} {s_['start_i'] / 1e6:.1f}-{s_['end_i'] / 1e6:.1f} Mb = {nm(s_['j'])} "
                         f"{s_['start_j'] / 1e6:.1f}-{s_['end_j'] / 1e6:.1f} Mb ({s_['n']} pairs, "
                         f"{'identical copies: very recent' if s_['k2p'] <= 0 else '~' + yrs(s_['k2p']) + ' ago'})"
                         for s_ in sg_[:6])
        out.append(f"Exchanges: {len(sg_)} stretch(es) where homoeologs carry near-identical copies of inherited "
                   f"insertions, so one homoeolog's sequence replaced the other's after the progenitors split "
                   f"(homoeologous exchange, or introgression before the merger): {desc}"
                   f"{'; ...' if len(sg_) > 6 else ''}. Timing vs the merger: "
                   + ", ".join(f"{v} {k.replace('_', ' ')}" for k, v in sorted(tim.items())) +
                   ". Benchmark: every call checked was real (432 of 432 in simulations, 5 of 5 in real genomes); "
                   "about 9 in 10 simulated replacements >= 1 Mb are found, fewer in real genomes.")
    else:
        out.append("Exchanges: no stretch of near-identical inherited insertions on the called homoeologs (replacing "
                   "exchanges >= 1 Mb are found 9 times in 10 in simulations; smaller or older ones can be missed).")
    for e in hb["halflife"]:
        if e["ok"]:
            out.append(f"Half-life of intact LTR retrotransposons in {group_label(e['h'])}: {e['est']:.4f} LTR "
                       f"divergence [{e['lo']:.4f}-{e['hi']:.4f}] = {yrs(e['est'])} [{yrs(e['lo'])}-{yrs(e['hi'])}] "
                       f"at mu {hb['mu']:.3g} (share of {group_label(e['g'])}'s pre-split copies whose "
                       f"{group_label(e['h'])} ortholog is still intact: {100 * e['r15']:.0f}%); intact copies per Mb "
                       f"in {group_label(e['h'])}: {hb['density'][e['h']]:.1f}.")
        else:
            out.append(f"Half-life in {group_label(e['h'])} (from {group_label(e['g'])}): not estimable "
                       f"({', '.join(k for k, v in e['qc'].items() if not v)} check failed).")
    for s_ in hb["asym"]:
        if s_["differ"]:
            out.append(f"Lineages differ: intact copies persist {s_['ratio']:.2f}x as long in {group_label(s_['a'])} "
                       f"as in {group_label(s_['b'])} [{s_['lo']:.2f}-{s_['hi']:.2f}]. Read this beside each "
                       f"subgenome's intact copies per Mb: a less complete assembly or annotation of one subgenome "
                       f"also shortens its half-life.")
        else:
            out.append(f"Lineages alike: half-lives in {group_label(s_['a'])} and {group_label(s_['b'])} do not differ "
                       f"detectably (ratio {s_['ratio']:.2f} [{s_['lo']:.2f}-{s_['hi']:.2f}]).")
    if hb["halflife"]:
        out.append("Half-life notes: copies the assembly or annotation missed count as lost, so these are lower "
                   "bounds on persistence; in simulations 94% of estimates that pass the checks fall within 25% of the "
                   "truth (splits of 1-5 Myr; very young splits read ~25% short).")
    if hb["transloc"]:
        out.append(f"Translocations: {len(hb['transloc'])} block(s) of a chromosome collinear with a chromosome other "
                   f"than its homoeolog: " + "; ".join(
                       f"{nm(t['i'])} {t['start_i'] / 1e6:.1f}-{t['end_i'] / 1e6:.1f} Mb with {nm(t['j'])} "
                       f"{t['start_j'] / 1e6:.1f}-{t['end_j'] / 1e6:.1f} Mb ({t['n']} pairs)"
                       for t in hb["transloc"][:6]) + ("; ..." if len(hb["transloc"]) > 6 else "") +
                   ". Benchmark: 16 of 19 such blocks in real genomes are supported by gene synteny.")
    if hb["many"]:
        out.append("Sequences with several homoeologs (fusion, fission or translocation): " + "; ".join(
            f"{nm(m['i'])} pairs with {' + '.join(names[p_] for p_ in m['partners'])} ({group_label(m['h'])})"
            for m in hb["many"][:8]) + ("; ..." if len(hb["many"]) > 8 else "") +
            ". Benchmark: 17 of 17 such cases in real genomes are supported by gene synteny.")
    return out


def main(argv=None):
    a = parse_args(argv)
    setup(a)
    t = {}
    t0 = time.time()
    lib = load_library(a)
    t["inputs_s"] = round(time.time() - t0, 1)
    t0 = time.time()
    cl = load_or_build_clades(a, lib)
    t["kmers_s"] = round(time.time() - t0, 1)
    t0 = time.time()
    ph = run_phasing(a, lib, cl)
    t["phasing_s"] = round(time.time() - t0, 1)
    t0 = time.time()
    try:
        ph["homoeologs"] = run_homoeologs(a, lib, ph)                 # v9: nothing else reads or changes it
    except Exception as e:                                               # noqa: BLE001 (never costs the v8 tables)
        ph["homoeologs"] = None
        log(f"WARNING: homoeolog pairing failed ({type(e).__name__}: {e}); every other output is complete")
    t["homoeologs_s"] = round(time.time() - t0, 1)
    t0 = time.time()
    dt = None if a.dev_skip_dating else run_dating(a, lib, ph)
    t["dating_s"] = round(time.time() - t0, 1)
    t0 = time.time()
    bio = run_biology(a, lib, ph, dt)
    t["biology_s"] = round(time.time() - t0, 1)
    t0 = time.time()
    pa = dt["painting"] if dt and dt.get("painting") else run_painting(a, lib, ph, dt)
    t["painting_s"] = round(time.time() - t0, 1)
    t0 = time.time()
    if ph.get("homoeologs"):
        try:
            ph["homoeologs"]["bio"] = run_homoeolog_biology(a, lib, ph, dt, pa)     # v9.1: reads, changes nothing
        except Exception as e:                                           # noqa: BLE001 (never costs the other tables)
            ph["homoeologs"]["bio"] = None
            log(f"WARNING: homoeolog biology failed ({type(e).__name__}: {e}); every other output is complete")
    t["homoeolog_biology_s"] = round(time.time() - t0, 1)
    write_phasing(a, lib, cl, ph, t, dt, bio, pa)
    if dt:
        write_dating(a, dt)
    write_biology(a, lib, bio)
    if pa:
        write_painting(a, pa)
    write_elements(a, lib, ph, dt, bio, pa)
    write_homoeologs(a, lib, ph)
    t0 = time.time()
    try:
        write_pdf(a, lib, ph, dt, bio, pa)
        log(f"report.pdf in {time.time() - t0:.1f} s")
    except Exception as e:                                                   # noqa: BLE001 (a figure never costs a run)
        log(f"WARNING: report.pdf failed ({type(e).__name__}: {e}); the tables are complete")
    log(f"done: {ph['N']} subgenome(s); outputs in {a.outdir} (report.pdf)")


if __name__ == "__main__":
    main()


# ======================================================================================================================
# NOTES FOR FUTURE WORK (v7 notes 2026-10-16; v8 and v8_mod 2026-10-19; v9 2026-10-20; v9.1 2026-10-21). This
# file is the only artifact kept and these notes are self-contained: assume none of the files, runs or data named below
# are reachable. Paths are provenance only (server asc-oulab, under /data2/chris/polyploid_genomes/).
# ======================================================================================================================
#
# V9.1 (2026-10-21): WHAT THE PAIRS SAY. run_homoeolog_biology (after dating and painting; reads, changes nothing;
# rng [seed, 91]); every v8 table and every v9.0 column unchanged (homoeolog_anchors.tsv gains 3 columns at the end).
# - Why: the user asked whether pairing lets the tool tell a one-genome user more biology, and to ship only what
#   survives the benchmark. Five candidates, bars fixed before scoring (harness LTR_phasing/dev_v91, memo.md;
#   design subphaser_LTR_FINAL2/design_v91.md): sims = 129 benchmark + 18 new (simallo + per-lineage half-life:
#   0.75 / 3 Myr, 1 vs 2 Myr); real = 36 polyploids, JCVI gene synteny (coordinates, strand, Ks), hecall PASS
#   exchanges (+ partner intervals), he_parent, 31 literature events; dev half first, held-out half once.
# - Shared: LTR-LTR divergence of every anchor (align_pairs, >= 200 sites). Near-identical anchors: empirical null on
#   log divergence of the chain anchors of called pairs per subgenome pair (median; scale = upper-half spread; flag
#   below median - z(0.05/n) x scale) or <= v8's merger upper limit. (A binomial test is far too narrow on real LTRs.)
# - Shipped: exchanged stretches = runs of >= 2 consecutive near-identical collinear pairs on a called pair (run p,
#   Poisson approximation, Bonferroni over called pairs): sims 432/432 real (chance 24%), recall >= 1 Mb 89%, < 1 Mb
#   58%, pre-merger introgression (N11) 5/5; real 5/5 (Arachis hypogaea 4, held-out 1; chance 7-16%) - real calls are
#   few because real LTR rate variation widens the null. Half-life per direction g -> h = h's loss of intact copies
#   (share of g's copies older than 1.5 x split whose inherited anchor lies on the homoeolog, minus chance; split =
#   median of the unflagged chain anchors, >= 20): 94% of 312 QC-passing sim estimates within 25% (truths 1, 1.5, 2,
#   3 Myr; 0.75 Myr withheld by QC; 0.5-Myr splits read -28%); Musa Silk vs Plantain agree; Arachis hypogaea vs
#   monticola disagree (fragmented monticola). Lineages differ (bootstrap 95% CI of the log ratio excludes 0): 2.0% of
#   148 equal-rate sim pairs flagged (shipped code 4/148); 6/6 asymmetric (2x) recovered (0.48-0.54 vs 0.50) and
#   flagged; real 2-3 of 14
#   flagged (read beside intact copies per Mb). Translocated blocks (collinear chain with a non-homoeolog, chain_test,
#   Bonferroni over tested pairs): gene-supported 8/10 dev, 7/7 held-out (the shipped code: 9/12, 7/7 - borderline
#   permutation p's fall either side of the level); 2 of 147 rearrangement-free sims had one.
#   Several homoeologs: 17/17 gene-supported. Collinear in reverse: 252/256 called pairs agree with gene order.
# - Dropped (dead_ends.md): split age from inherited pairs (sims 1.0% error but CI coverage 47%; dev Spearman vs Ks
#   0.05 vs v8 0.53; held-out MAD 0.058 vs v8 0.189 = a lead needing new held-out genomes); typing v8's exchange calls
#   (sims 97%, real 8/10 dev, 0/4 held-out); inversions (precision 53% / 54%).
# - Outputs: homoeolog_exchanges.tsv, homoeolog_halflife.tsv, homoeolog_rearrangements.tsv, findings.txt (also logged),
#   summary.json homoeologs.biology; report.pdf B: near-identical pairs in FIG_EXCH, translocated blocks framed dotted.
# - Final benchmark (2026-10-21, 188 runs = 170 + 18 half-life sims, 0 failures; 10:58-14:28, budget 96 then 20
#   threads on a shared server): every v8 table and v9.0 column identical in all 188 runs (the 18 new sims differ
#   only in the input paths recorded in summary.json); the shipped outputs reproduce the scrutiny: stretches sims
#   432/432, real 5/5, recall >= 1 Mb 263/293, < 1 Mb 158/273; half-life 294/312 within 25%; translocated blocks
#   real 16/19; several homoeologs 17/17. The benchmarked copy (LTR_phasing/bench_v91/pinned) differs from this file
#   only in the translocation count quoted in findings.txt (15 of 17 -> 16 of 19) and these notes.
# - Tests: tests/test_v9_homoeologs.py (13): + planted near-identical values and runs; a planted 1-Mb replacement found
#   inside its interval and nowhere else; no stretch without exchanges; version 9.1.0.
#
# V9 (2026-10-20): HOMOEOLOG PAIRING. ltrharvest_phasing9.py = v8_mod + run_homoeologs; every v8 table is unchanged
# (byte-identical apart from summary.json 'version'/'homoeologs' and the exchanges.gff3 source column 'v9').
# - Why: v5-v8_mod phased chromosomes into subgenomes but never said which chromosome pairs with which (Poa annua:
#   7 + 7 phased, no 1A-1B pairs). v5/v6 took homoeolog sets as user input; v8's partner = within-subgenome relative.
# - Biology (tested on the dev half, 22 polyploids, MCscan gene-anchor truth): most copies postdate the progenitors'
#   split, are subgenome-specific and spread over every chromosome of their subgenome, so the overall LTR-RT footprint
#   makes a subgenome's chromosomes alike, not a chromosome like its homoeolog. The homoeolog-specific part: insertions
#   inherited from the progenitors' ancestor, at collinear positions on both homoeologs (decisive), and the landscape
#   (TE density / clade mix / age along the chromosome; above chance, never decisive: Oryza minuta 13-15/24 top-1,
#   chance ~2; Aegilops 4-5/14). Inherited insertions are deleted over time, so the signal fades with split age:
#   complete for Aegilops, Arachis, Coffea (gene Ks 0.02-0.07); none for Oryza BC, Brassica B/C, Nicotiana,
#   Brachypodium, A. suecica (Ks 0.08-0.19); simulations (half-life 1.5 Myr): 2-Myr splits 100%, 5 Myr ~85%,
#   8-12 Myr 10-40% recall, no false call.
# - Rule (run_homoeologs): anchors = reciprocal closest relatives between copies of the old LTR-K2P tercile in
#   different subgenomes (old_relatives: shared lineages <= NMAX_KNN carriers in the whole library, weight
#   mult/(size-1), ties -> lower copy index); per cross-subgenome chromosome pair, anchors observed/expected
#   (degree-corrected in the subgenome-pair block) and the longest collinear chain in anchor-rank space (CHAIN_W; z vs
#   NULL_PERMS permutations); per chromosome x other subgenome, candidate = largest z, called 'homoeolog' when it also
#   has the largest O/E and chain_test p <= 0.05 / |candidates| (n <= 25 anchors: 2 x exact P(LIS >= l) by RSK /
#   hook lengths; more: B = ceil(2/alpha) permutations, early stop), else 'unresolved'; 'none' without >= MIN_ANCHORS
#   anchors to any chromosome. Reciprocal flag; sets = components of calls. Nothing forced 1:1 (Brassica, Brachypodium,
#   A. suecica, Camelina are not 1:1).
# - Development (dev half, v9 with dating skipped): 234/234 calls right, recall 234/767 = 0.305 (0.31 incl. old-split
#   genomes that cannot be paired); simulations 2323/2323, recall 0.962. Cost: ~2-25 s (Avena, 463k copies: 22 s).
# - Final benchmark (2026-10-20, end-to-end with dating, 170 runs, 0 failures, 49 min at 96 threads; truth: MCscan
#   gene anchors, a partner holding >= 25% of the chromosome's largest anchor count; rows = chromosome x other
#   subgenome; 36 of 41 polyploids have >= 2 subgenomes): dev 234/234 right, recall 0.305; held-out half (scored once)
#   267/268, recall 0.561 (Solanum 72/72, Musa 22/22 x2, Saccharum 18/18, Prunus 16/16, Glycine 33/40, Sphagnum
#   35/38, Populus 26/38, Poa annua 6/14: 3A-3B, 6A-6B, 7A-7B paired, 1/2/4/5 unresolved); simulations 2323/2323,
#   recall 0.962. The one wrong call: Camelina s011 (chr04SG2) -> s003 (chr04SG1, truth flag 'check'), 0 gene anchors
#   but 9 collinear inherited insertions over 8-23 Mb (possible fusion/truth ambiguity). Recall per subgenome pair
#   follows inherited insertions on homoeologs beyond chance per 1,000 old copies (Spearman 0.84, 51 pairs) more than
#   gene Ks (-0.65); Ks < 0.05: 0.66 paired, 0.05-0.08: 0.47, >= 0.08: 0.05. Every v8 table identical in all 170 runs
#   (bench_v8's dating.tsv predates v8's writer fix on 4 runs: the documented demoted-bound rows). Homoeolog step
#   0.1-12 s per genome (~1% of a run).
# - Rejected (dead ends): landscape profiles; all-copy lineage sharing + one-to-one matching (the v8_mod prototype:
#   Poa 4/7, Oryza minuta 4/12); small-lineage link counts (hub chromosomes: diffuse pre-split paralogs and family
#   co-enrichment, not collinear); clade bootstrap support (279/324 at 0.95) and block bootstrap (168/173, recall 0.22);
#   RBH among all copies (301/308: young C-D movers in Avena, 5 false calls); genome-wide Bonferroni (210/210, recall
#   0.27); v6's normal approximation of the 50-permutation chain z (3-7x too liberal for 26-60 anchors on simulated
#   non-homoeologs); an 'also_collinear' column (flagged a non-homoeolog in 1 of 8 rows of a clean simulation).
# - Outputs: homoeologs.tsv, homoeolog_sets.tsv, homoeolog_anchors.tsv, summary.json 'homoeologs'. report.pdf: B is a
#   dot plot of the anchors (homoeologs = collinear diagonals; called pairs framed), chromosomes ordered so homoeologs
#   share a position in A, B and G, set numbers under A's bars; G now carries the lineage-sharing key. With one group
#   (or a failed homoeolog step) the page is v8_mod's.
# - Tests: subphaser_LTR_FINAL2/tests/test_v9_homoeologs.py (AST guard vs v8_mod; LIS tail vs brute force; chain test
#   level on random anchors; RBH vs brute force with ties, blocks, workers; planted pairs and a shuffled decoy; no
#   shared insertions; e2e: simulated allotetraploid all pairs right + v8 tables identical to v8_mod, 15-Myr split
#   no wrong call, one group, failure isolation; figure with/without homoeologs).
#
# V8_MOD (2026-10-19): REPORT ONLY. ltrharvest_phasing8_mod.py = v8 with a one-page report.pdf; every table is v8's.
# - Why: v8's report (figure + short and long legend pages) showed method internals (copies as ticks, a split-test
#   glyph grid, bare date bars) more than the biology; the user asked for one self-explanatory page, phasing first.
# - Page (180 x 240 mm, one-line header per panel, keys on the page, provenance footer): A every phased sequence to
#   scale, painted by the posterior ancestry tiles of the bedGraphs (opacity by posterior; hatched = copy-free stretch
#   longer than the median resolution or ~2 pixel rows), LTR-RT density at its left, called exchanges outlined with a
#   donor-colored triangle; B lineage sharing (observed/expected, all copies; links_oe's matrix); C PCA of each
#   sequence's log O/E sharing profile (the masked matrix the partition maximises; missing self-entry = column mean)
#   with every half-library projected on the same axes (clouds = placement stability; support >= SUPPORT_OK counted);
#   D the two tests of every tested split (adjusted Rand >= MIN_REP; strength >= MIN_REL x the split above), on the
#   copies the decision used; E-H share one log time axis (K2P / 2 mu, oldest left) and right margin: E deconvolved
#   insertion history per subgenome (+ copies), F the six largest annotated clades (+ copies), G B's sharing within the
#   LTR-K2P terciles placed at their ages, H the closest-relative profile (cross_profile of v7's pair table, Wilson
#   95%) and every split/merger estimate of dating.tsv (bold = reported; > 2 subgenomes: each pair's reported ones,
#   six pairs at most). With two subgenomes, the reported split and merger intervals are bands across E-H and a strip
#   names the periods. N = 1, --assign, -n, no dating and no painting each degrade to messages, never to a crash.
# - The one non-figure change: run_phasing keeps the 100 half-library link matrices (float32; None above
#   FIG_HALVES_MAX = 300 phased sequences) as ph['halves'] for panel C. Nothing else differs: tests/
#   test_v8_mod_report.py compares every top-level definition with v8 (AST) and every output file of v8 and _mod runs
#   on simulated genomes byte for byte (summary.json without timings); the full Poa annua library gave identical
#   tables too. Drawn and inspected on Poa annua (N=2), Solanum demissum (3), Fragaria x ananassa (4) and Actinidia
#   deliciosa (1): fig_v8/final (development harness: subphaser_LTR_FINAL2/fig_v8, memo.md there).
# - Homoeologs (v9 now infers them; see V9): no version (v5-v8) infers homoeolog pairs; v5/v6 took them only as user
#   input (--config). A one-to-one matching of chromosomes across subgenomes on log O/E lineage sharing, with
#   half-library support, recovered the
#   numbered homoeologs in Solanum 36/36, Mimulus 14/14, Coffea 11/11, Musa 11/11, Arachis 10/10, Prunus 8/8, but
#   Fragaria 33/42, Glycine 15/20, Poa 4/7, Oryza minuta 4/12 (support >= 0.7 kept Poa's 2 right pairs). A lead, not
#   implemented (the user kept this edit to the PDF).
# - Gotcha: importing this script under another module name (importlib, e.g. 'p8') makes numba cache its kernels under
#   that name, and every later run then dies with "No module named 'p8'": import it under its own name, or delete the
#   *.nbi/*.nbc files of its __pycache__.
#
# V8 (2026-10-19): DATING ONLY (everything else is v7; non-dating tables byte-identical to v7 on the benchmark)
# - Why: v7 dated 9/41 mergers and none younger than ~20 kyr (each copy timed by its own LTR pair: one substitution
#   per 0.5-1 kb ~ 20-80 kyr); its reported split (confinement) was the least accurate of its three estimators.
# - Genealogy clock (sibling_genealogy, younger_cdf_genealogy): a sibling pair's divergence over the whole records
#   bounds time ~10x more finely, but dates the pair's common source; the younger copy was inserted a fraction r of
#   that ago. f(r) is learned by EM from the genome's older sibling pairs (>= 5 expected substitutions per LTR) and read
#   into every young pair; pair divergence is scaled to LTR units by the measured LTR/whole K2P ratio (ltr_rate_ratio;
#   real genomes median 1.49, IQR 1.25-1.82, 29/31 subgenome pairs > 1; simulations 1.00). Learned genealogies differ
#   between genomes (mean r 0.12 Brachypodium - 0.57 Fragaria; Fragaria alone has a large r ~ 1 copy-of-copy share).
# - More young focal copies (YOUNG_MAX 12,000) on their own stream; bootstrap units from v7's pair table (denser
#   sampling chained families: Arachis largest unit 8% -> 57-65% of focal copies).
# - Datable also needs a decisive likelihood ratio (>= BF_DECISIVE) for the step: without it the clock dated the
#   55-yr synthetic A. magna x longiglumis A-D pair from a few cross pairs with zero background (log LR 0.0).
# - Reported merger (merger_report): genealogy clock two-sided, else LTR clock two-sided, else the tighter bound; a
#   lower limit below one substitution per LTR becomes a bound (simulations put such lower limits too old in TE-poor
#   genomes, hexaploids and after introgression). MODEL_ERR merger 0.24 re-derived for the genealogy clock (r1-r4
#   coverage 0.91; held-apart r5-r8 0.98).
# - Reported split: orthologs when their interval overlaps the profile's, else profile, confinement, orthologs
#   (simulations: utility 286 vs 252 for confinement first; development accessions 19 vs 12 against hecall
#   homoeolog divergence and 19 vs 7 against gene Ks; orthologs are near exact in simulations but contaminated on
#   some real genomes, e.g. A. suecica 0.019 vs gene Ks 0.136).
# - Final benchmark (end-to-end runs; 151 simulated subgenome pairs with exact truth; 41 real polyploids, development
#   and held-out halves; utility +2 narrow / +1 wide or bound / 0 none / -3 wrong):
#   mergers, simulations: v7 utility 147 (dated 59, bound 51, none 36, wrong 5); v8 167 (48, 87, 11, 5). Mergers
#   <= 20 kyr (62): v7 33 bounds, 27 none, 2 wrong; v8 57 bounds (median upper limit 31 kyr, 23 below 20 kyr; 1-kyr
#   mergers < 3-6 kyr, 50-yr synthetic lines < 3-14 kyr), 3 none, 2 wrong. Cost: some correct 50-200 kyr dates became
#   correct bounds (> 20 kyr utility 120 -> 116). Real mergers: held-out half v7 3 bounds / 6 none, v8 8 / 1, no wrong;
#   development v7 5 dated / 6 bound / 11 none, v8 5 / 9 / 7 / 1 wrong (Coffea: < 18 kyr vs 350-610 kyr in the truth
#   table; 10-20 kyr in an earlier estimate). Splits, simulations: wrong 10 -> 3, width 0.17 -> 0.03 log10. Real
#   splits: development half better (utility vs d_hom 12 -> 24, vs gene Ks 7 -> 24), held-out half worse (40 -> 31,
#   30 -> 21; profile ~2x old in Camelina, Poa, Saccharum, Musa); all 49 pairs: wrong 9 -> 7 (d_hom), 12 -> 9 (Ks).
#   SPLIT_RULE = "confinement" restores v7's reported split. Runtime median = v7, max 20 vs 12 min; memory = v7.
# - dating.tsv: the reported merger row carries the reported merger (a demoted lower limit appears as a bound); the
#   benchmarked build wrote the clock's raw interval there (summary.json was right), fixed and rechecked on 11 runs.
# - Dead ends (do not retry): painting-relabelled clock (25 wrong of 85 pairs: exchange-induced cross pairs carry
#   genuine post-merger timing); cross-share profiles along pair divergence with an assumed genealogy (hyperbolic
#   runs young, step old; both dated synthetic Raphanobrassica); pair divergence as an upper bound with the library
#   age prior below it (~2x old); exchange-duplicate lower bound on the merger (synthetic lines carry same-insertion
#   pairs at 3e-3 - 3e-2); targeted ortholog search (picks young duplicates when split-era insertions are gone);
#   interpolated half-rise (keeps a ~15% young bias, narrower intervals: 8/25 missed); split in LTR units.
# - Caveat for every merger: the clocks date the onset of cross-subgenome transposition. Synthetic lines (Saccharum,
#   Raphanobrassica, C. auratus x Cyprinus) carry young sibling pairs shared between subgenomes that no transposition
#   in < 120 years can make (gene flow or TE transfer before the polyploid formed, or misassigned contigs); v8 reports
#   them as bounds containing zero, never as dates, but a natural polyploid's date can be such an onset.
# - Tests: subphaser_LTR_FINAL2/tests/test_v8_dating.py (new functions; V8_SLOW=1 adds a simulated 10-kyr merger end
#   to end) and test_v8_ported.py / test_k_rule_v8.py (v7's suite ported; 11 dating tests updated to v8 behavior);
#   159 tests pass (~21 min with V8_SLOW=1).
#
# ENVIRONMENT
# - python 3.11, numpy >= 2, scipy >= 1.11, numba >= 0.60, matplotlib >= 3.8, python-edlib >= 1.3; pytest for tests;
#   pdfinfo/pdftotext (poppler) only to check report.pdf. Developed and benchmarked with numpy 2.4.6, scipy 1.17.1,
#   numba 0.66.0, matplotlib 3.11.1.
# - Lineage of the code: v4 -> v5 (fast, pooled dating; optional --genome painting, dropped) -> v6 (exact Gray-code
#   2-way search, pairing map, copy-pair dating with three split estimators and a sibling merger clock, calibrated
#   painting) -> v7 = v6's pipeline + the age-resolved number-of-subgenomes rule + one PDF + numba threads = -t.
#
# INPUT FORMATS THE BENCHMARK USED (LTRquest output)
# - FASTA headers chrom:start-end#Class/Superfamily/Clade (e.g. s001:21691-27931#LTR/Copia/Bianca); clade names from
#   DANTE/REXdb (Tekay, CRM, Athila, Retand, Ogre, Ale, Angela, Bianca, Ivana, SIRE, Tork, TAR, ...; "unknown" when no
#   domains; fish/fungal libraries are mostly unknown/Gypsy-unknown/chromo-outgroup).
# - K2P table (TSV, header): seq_id seq_len status ltr5_start ltr5_end ltr3_start ltr3_end ltr5_len ltr3_len flank5_len
#   flank3_len aln_len n_sites n_ts n_tv n_gapcols identity p_dist k2p k2p_se bitscore ... strand family domains
#   nest_status. GFF3 alternative: LTR_retrotransposon/repeat_region features with K2P_d (or k2p), n_sites, k2p_se, ti,
#   tv, ltr5_len, ltr3_len, clade, superfamily and long_terminal_repeat children (Parent=). Both are read by
#   read_annotation(); >= 95% of FASTA copies must match by chrom:start-end.
#
# BENCHMARK DATA AND TRUTH (bench_indep, rebuilt from primary sources; inputs in LTR_phasing/bench_v6/data/<acc>/)
# - Libraries: LTRquest depth0-2 clean LTR-RT libraries of each polyploid alone; sequence names blinded to s001...,
#   family labels removed (they had been pooled with progenitors). 41 polyploid accessions from 38 complexes, 6 extra
#   polyploids, 67 controls.
# - Phasing truth: each chromosome's subgenome from progenitor alignment (cleaned chromosome names, 2026-09-30), with
#   confidence tiers; Rosa x hybrida haplotypes are lineage mosaics (excluded from totals). The token key was rebuilt
#   from sequence (length + md5 of the first and last 20 kb).
# - Exchange truth: hecall v0.4.10 PASS calls (a progenitor-alignment caller; reproduces 30/31 published events) plus 30
#   published events. Trap: Glycine dolichocarpa chr8 carries a rearrangement shared with the diploid parent; any call
#   there is false (v5/v6/v7 make 1). hecall calls its own truth uninterpretable for Rosa x hybrida, Magallana, GR-A3n,
#   A. agadiriana and A. monticola (excluded). Magallana's name labels are incoherent (its "A" haplotype is mostly
#   gigas origin): against hecall's majority origin, k = 1 is right.
# - Date truth: literature (4 searches): deepest published split; youngest published merger (formation of the sequenced
#   polyploid); only events between the accession's own subgenomes (e.g. the 18th-century F. chiloensis x virginiana
#   cross is not one); Digitaria mergers excluded (dated with SubPhaser itself, circular); Prunus keeps the run rate.
#   All tools date as K2P/(2 mu): years were rescaled from mu_run to the literature mu (mu_lit).
#
# PER-ACCESSION RECORD (k true; LTR-RT copies; calls as k_hat/chromosome accuracy for v5, v6, v7 and SubPhaser, which
# was GIVEN k and the homoeolog sets; literature split and merger in Mya; run and literature mu, substitutions/site/yr).
# Extras were scored for v7 only (no v5/v6/SubPhaser or literature rows).
#   accession (k true)                         n_LTR   v5      v6      v7      SubPh.  split  merger  mu_run/mu_lit
#   Actinidia_deliciosa                      2   28141 1/0.50  1/0.50  1/0.50  -         7.1     2.1  9e-10/3.39e-09
#   Aegilops_geniculata                      2  352039 2/1.00  2/1.00  2/1.00  2/1.00   4.21    1.28  1.3e-08/1.3e-08
#   Aegilops_ventricosa                      2  335376 2/1.00  2/1.00  2/1.00  2/1.00    4.6     0.7  1.3e-08/1.3e-08
#   Arabidopsis_kamchatica                   2   11526 2/1.00  2/1.00  2/1.00  2/1.00  0.337     0.1  7e-09/7e-09
#   Arabidopsis_suecica                      2    2903 2/1.00  2/1.00  2/1.00  2/1.00   5.97   0.016  7e-09/7e-09
#   Arabidopsis thaliana x arenosa Allo738   2    3329 2/1.00  2/1.00  2/1.00  2/1.00   5.97 2.3e-05  7e-09/7e-09
#   Arachis_hypogaea                         2   94063 2/1.00  2/1.00  2/1.00  2/1.00   2.16  0.0094  1.64e-08/1.6e-08
#   Arachis_monticola                        2   66749 2/1.00  2/1.00  2/1.00  2/1.00   2.16  0.0117  1.64e-08/1.6e-08
#   Avena_agadiriana                         2  276896 1/0.50  1/0.50  1/0.50  2/0.93   3.11    1.53  1.3e-08/1.3e-08
#   Avena_magna_x_longiglumis                3  462947 2/0.67  2/0.67  3/1.00  3/1.00   10.9 5.5e-05  1.3e-08/1.3e-08
#   extra: Avena maroccana 4x                2  314628 -       -       2/1.00  -           -       -  1.3e-08/-
#   Brachypodium_hybridum                    2    7398 2/1.00  2/1.00  2/1.00  2/1.00   11.5    0.14  1.3e-08/1.3e-08
#   Brassica_carinata                        2   33958 2/1.00  2/1.00  2/1.00  2/1.00    8.5  0.0205  1.5e-08/1.5e-08
#   Brassica_juncea                          2   23203 2/1.00  2/1.00  2/1.00  2/1.00    8.5   0.011  1.5e-08/1.5e-08
#   Brassica_napus                           2   23303 2/1.00  2/1.00  2/1.00  2/1.00    4.6  0.0075  1.5e-08/1.5e-08
#   Camelina_microcarpa_6x_Type1             3   13935 2/0.65  2/0.65  3/1.00  3/1.00    5.2   0.065  8.22e-09/8.22e-09
#   extra: Camelina microcarpa 4x            2    7474 -       -       2/1.00  -           -       -  8.22e-09/-
#   Carassius_GR-A3n_synthetic               3   41828 2/0.82  2/0.82  2/0.82  3/0.62   3.75   3e-06  3.51e-09/3.51e-09
#   Carassius_auratus_x_Cyprinus_carpio      2   21753 2/1.00  2/1.00  2/1.00  2/1.00     11 3.4e-05  3.51e-09/3.51e-09
#   Coffea_arabica                           2   31810 2/1.00  2/1.00  2/1.00  2/1.00    4.9   0.462  1.3e-08/1.3e-08
#   Cryptococcus neoformans x deneoformans   2     137 1/0.12  1/0.12  1/0.12  2/0.77   24.5       1  2e-09/2e-09
#   Digitaria_sanguinalis__4x                2   23249 2/1.00  2/1.00  2/1.00  2/1.00    6.4       -  1.3e-08/1.3e-08
#   Digitaria_sanguinalis__6x                3   31670 3/1.00  3/1.00  3/1.00  3/1.00    6.9       -  1.3e-08/1.3e-08
#   extra: Digitaria milanjiana 4x           2   22664 -       -       2/1.00  -           -       -  1.3e-08/-
#   Fragaria_x_ananassa                      4   19492 4/1.00  4/1.00  4/1.00  4/1.00   6.37     0.8  7e-09/7e-09
#   extra: Fragaria chiloensis 8x            4   23504 -       -       4/1.00  -           -       -  7e-09/-
#   extra: Fragaria virginiana 8x            4   18040 -       -       4/1.00  -           -       -  7e-09/-
#   Glycine_dolichocarpa                     2   48341 2/1.00  2/1.00  2/1.00  2/1.00   4.71    0.38  6.1e-09/1.3e-08
#   Lemna_x_japonica__7182_3x                2   25150 2/1.00  2/1.00  2/1.00  2/1.00    7.5       -  1.1e-08/1.3e-08
#   Lemna_x_japonica__8627_3x                2   25159 2/1.00  2/1.00  2/1.00  2/1.00    7.5       -  1.1e-08/1.3e-08
#   Magallana_gigas_x_Magallana_angulata     2   12758 1/0.50  1/0.50  1/0.50  2/0.55   4.82   8e-06  1.39e-09/2e-09
#   Mimulus_sookensis                        2   13556 2/0.96  2/0.96  2/0.96  2/1.00    0.2   0.071  1.5e-08/1.5e-08
#   Musa_x_paradisiaca__Plantain_AAB_3x      2   30081 2/1.00  2/1.00  2/1.00  2/1.00    4.4 0.00401  4.5e-09/1.3e-08
#   Musa_x_paradisiaca__Silk_AAB_3x          2   31034 2/1.00  2/1.00  2/1.00  2/1.00    4.4   0.007  4.5e-09/1.3e-08
#   Nicotiana_tabacum                        2  131443 2/1.00  2/1.00  2/1.00  2/1.00     10     0.6  1.3e-08/1.3e-08
#   Oryza_malampuzhaensis                    2   19761 2/1.00  2/1.00  2/1.00  2/1.00    6.1    1.73  1.3e-08/1.3e-08
#   Oryza_minuta                             2   22160 2/1.00  2/1.00  2/1.00  2/1.00    6.1    0.71  1.3e-08/1.3e-08
#   Oryza_punctata_4x                        2   18448 2/1.00  2/1.00  2/1.00  2/1.00    6.1   0.296  1.3e-08/1.3e-08
#   Poa_annua                                2   62948 2/1.00  2/1.00  2/1.00  2/1.00   5.64    0.05  1.3e-08/1.3e-08
#   Populus_tomentosa_clone741               2   11318 2/0.84  2/0.84  2/0.84  2/0.47    8.4 5.1e-05  2.5e-09/2.5e-09
#   extra: Populus tomentosa 2x              2   11071 -       -       2/1.00  -           -       -  2.5e-09/-
#   Prunus_cerasus                           2   14547 2/1.00  2/1.00  2/1.00  2/1.00   8.74    1.61  7.77e-09/7.77e-09
#   Rosa_x_hybrida                           4   59434 1/0.25  1/0.25  1/0.25  4/0.39      -       -  4.02e-09/7e-09
#   Saccharum_hybrid                         2   47649 2/1.00  2/1.00  2/1.00  2/1.00   1.56 0.000116  1.3e-08/1.3e-08
#   Solanum_demissum                         3   53955 2/0.61  2/0.67  3/1.00  3/0.97    7.1       -  1.3e-08/1.3e-08
#   Sphagnum_troendelagicum                  2    5005 2/1.00  2/1.00  2/1.00  2/1.00     16    0.04  1.6e-09/9.4e-09
#   xRaphanobrassica_RRCC                    2   26773 2/1.00  2/1.00  2/1.00  2/1.00   8.59 4.3e-05  1.5e-08/1.5e-08
# - Controls (67; v7 k = 1 unless noted): diploid/haploid (55): Actinidia chinensis, A. deliciosa 2x, Aegilops comosa,
#   Ae. tauschii, Ae. umbellulata, Ae. uniaristata, Arabidopsis arenosa, A. halleri, A. lyrata, A. thaliana, Arachis
#   duranensis, A. ipaensis, Avena canariensis, A. hirtula, A. longiglumis, Brachypodium distachyon, B. stacei, Brassica
#   nigra, B. oleracea, B. rapa, Camelina hispida, C. neglecta, Coffea canephora, C. eugenioides, Cryptococcus
#   deneoformans 1x, C. neoformans 1x, Digitaria radicosa, Erythranthe guttata, E. nasuta, Glycine syndetika,
#   G. tomentella, Lemna minor, L. turionifera, Magallana angulata, M. gigas (k=2), Musa acuminata banksii,
#   M. a. malaccensis, M. balbisiana, Nicotiana sylvestris, N. tomentosiformis, Oryza eichingeri (k=2), O. officinalis,
#   O. punctata 2x, Poa infirma (k=2), P. supina, Populus alba, Prunus avium, Raphanus sativus (k=3), Rosa chinensis
#   'Old Blush', R. lucieae, R. multiflora, R. odorata (k=2), Solanum verrucosum, Sphagnum balticum 1x, S. tenellum 1x.
#   Polyploid controls without truth (12): Carassius auratus 4x, C. auratus red var. 4x, C. cuvieri 4x, C. gibelio 6x,
#   Cyprinus carpio 4x, Prunus fruticosa 4x, Rosa gallica 4x (two assemblies: k=2 and k=3), Saccharum officinarum 8x,
#   S. spontaneum (autopolyploid), Solanum acaule 4x (two).
#
# BENCHMARK RESULTS (v7: 114 runs, 8 threads each, 120-thread budget, 23 min wall, 0 failures; v5/v6/SubPhaser from
# the independent benchmark, which ran them on the same inputs)
#                                          v5            v6            v7            SubPhaser (given k, homoeologs)
#   right k (40 non-mosaic polyploids)     32            32            35            given
#   perfectly phased                       30            30            33            33/39
#   mean chromosome accuracy               0.904         0.906         0.931         0.957
#   diploid controls with false structure  5/55          5/55          5/55          no null run (needs >= 2 groups)
#   exchange bp precision / recall         0.531/0.793   0.544/0.791   0.533/0.815   0.35/0.81 (calls 2.4x the truth Mb)
#   exchanges found 1-10 Mb, 0.1-1 Mb      37/90, 58/217 37/90, 61/217 45/90, 64/217 59/90, 70/217
#   published exchanges found (>=50%)      20/28         20/28         22/28         24/28
#   split dated (within 2x of literature)  29 (10)       36 (11)       36 (12)       36 (24; its percentile rule, 3x
#                                                                                    the divergence, see dating)
#   median wall / CPU (polyploids)         1.6 min/0.06h 4.0 min/0.42h 5.0 min/0.52h 5.9 min core; 63-300 GB RSS
# - Scoring rules: chromosome accuracy after the best one-to-one matching of predicted to true labels (Hungarian,
#   length-weighted; unassigned = wrong); exchanges on 10-kb bins (bin-centre rule) against hecall PASS, donor-aware,
#   events "found" when a donor-matched call overlaps; published events found when >= 50% covered; controls false when
#   k_hat > 1; split = the deepest reported pair split, merger = the youngest dated merger, compared on log10(tool/lit)
#   after rescaling to mu_lit; literature bounds scored as consistent or not.
# - v7 changed exactly 3 accessions vs v6, each now perfect: Avena_magna_x_longiglumis (k 2 -> 3; A|D replicated in
#   young copies at 0.76x the parent group's strongest young split, absent in old; pooled only 0.20x its parent),
#   Camelina_microcarpa_6x_Type1 (2 -> 3; SG2|SG3 seen in middle-aged and old copies, not young), Solanum_demissum
#   (2 -> 3; SG1|SG2 only in young copies; pooled, old copies link homoeologues and the best split of SG1+SG2 runs
#   orthogonal to the subgenomes, 6 SG1 + 6 SG2 per side). Painting accuracy there 0.62/0.68/0.63 -> 0.97/0.99/0.95;
#   the published Solanum chr2/chr3 swaps are now found (97%, 92% covered; v5/v6 0%); A. magna x longiglumis exchange
#   precision stays very low (0.04; v6 0.005); Camelina gains 3 false calls (3.1 Mb). Elsewhere v7's tables are
#   byte-identical to v6's (checked on 8 accessions; only chromosomes.tsv partner_support differs: extra half pairs
#   come from the main rng, v6 used a separate one). Extras 6/6 right. One polyploid control without truth (Rosa
#   gallica 4x, GCA_055540455, a 7-sequence haplotype-mosaic assembly) gains a third group from a young-only 3|2 split:
#   the one uncheckable call the added sensitivity makes.
# - Still wrong, and why no k rule fixes them: Actinidia deliciosa (best replicability 0.46 at any age; A/D ~97%
#   identical; AAAADD), Avena agadiriana (pooled replicability 0.76, "cannot tell"; middle-aged copies 0.98 but the top
#   split is never rescued by one age class; SubPhaser gets 0.93 from genome k-mers), Carassius GR-A3n (G and Au, both
#   Carassius, never split reproducibly in any age class, not even with their 15 near-identical pairs kept apart, 0.59
#   agreement; verdict flags "more subgenomes than resolved"; v6 forced to k=3 reached 0.93), Cryptococcus AD (137
#   LTR-RTs, 3 phasable sequences; v5 -n 2 crashed with IndexError there), Magallana (truth incoherent, above).
#   Given the true k, v5/v6 reach only 0.50-0.57 on Actinidia and A. agadiriana: signal, not the rule.
# - Diploid false structure (same five as v6): Oryza eichingeri (6|6, strength 0.58, carried by 5 independent TE clades,
#   strongest in young copies: probably a hybrid/haplotype-mosaic assembly, i.e. real structure), Poa infirma (5|2,
#   strength 0.04 but replicated over 39,521 copies), Raphanus sativus (k=3), Rosa odorata (R. chinensis x gigantea
#   hybrid origin; young-biased 4|3), Magallana gigas (borderline replicability 0.80-0.85). v5 had Rosa chinensis
#   instead of M. gigas; v5's k on near-null genomes depended on FASTA record order (3/67 controls changed on re-sort).
# - Dating: v5/v6 splits are 2.1x young vs the literature, but their divergence equals the homoeolog SNV divergence
#   measured by hecall (ratio 1.04-1.11) and that proxy is itself 2.3x young with the same rates: the offset is the
#   divergence-to-years conversion, not the estimators. Spearman vs literature v5 0.65, v6 0.60. SubPhaser's rule (mean
#   97.5th percentile of subgenome-specific LTR ages) overshoots divergence 3x and lands near published years for the
#   wrong reason. Large misses: A. suecica/Allo738 (0.7-0.8 Myr vs 6 Mya), Prunus 8x young, Mimulus and A. kamchatica
#   3-5x old (the two shallowest published splits).
# - Mergers: no tool dates mergers younger than ~20 kyr (15 accessions; 8 made in the last ~120 years: Allo738,
#   A. magna x longiglumis, RRCC, C. auratus x Cyprinus, GR-A3n, Magallana, Populus 741, Saccharum). LTR clocks cannot:
#   identical LTRs (K2P = 0) reach back 18-550 kyr (half still identical at ln2/(2 mu L)); in those 8 lines,
#   identical-LTR copies at sites empty in the progenitor reference (5,730; 44-97%) are mostly parent-vs-reference
#   polymorphisms, not post-merger insertions (needs the actual parents); v6's "foreign" copies are not enriched there.
#   v6/v7 report "not datable" or bounds. Older mergers: O. minuta 0.72 Myr vs 0.71 Mya (good); Fragaria 3.6 Myr vs
#   0.8 Mya (matches the older 3.0 Mya CD step instead). Scoring quirk: the benchmark compares the youngest DATED merger
#   with the literature's youngest event; v7 dates A. magna's own C-D merger (448 kyr, 273-725) while the ~55-yr
#   synthetic cross is reported, correctly, as not datable (< 23 kyr under random landing).
# - Runtime: dating (copy-pair sketches + alignments + bootstraps) dominates (B. napus 160 of 190 s). v7's age-class
#   link matrices add ~15-60 s (4 sources x 50 half pairs). Largest library (A. magna x longiglumis, 463k LTR-RTs) 12
#   min on a loaded server. v7 used 3.8 cores per 8-thread job (median; max 7.1); v6 burned 4-25 cores at -t 8 (its
#   numba pool was sized to all 256 cores). SubPhaser needs the genome plus a homoeolog config (one column per
#   subgenome); its k-mer step peaks at 200-440 GB on 2.5-11 Gb genomes; its LTR dating costs ~5 thread-min per 20 Mb
#   (ltrharvest); defaults fail on Actinidia (no k-mer with fold > 2) and Cryptococcus (fallback -q 50 recovers only
#   Cryptococcus).
#
# DEAD ENDS (number-of-subgenomes work). Method: an offline harness held, per accession, the link matrices for all
# copies and for each LTR-K2P tercile with 30 copy-half pairs, so a rule could be scored on all 114 accessions in
# ~1 min; "k right /40" below is from that harness (v6's rule scored 33 there, 32 in the full benchmark).
# - Drop the 0.5x-parent rule: 29/40. Splits five correct genomes: Brassica A05/A06/A09 in B. napus AND B. juncea
#   (replicability 1.00; an inherited B. rapa compartment: these three chromosomes share lineages exclusively, at every
#   age, across many TE clades, pericentromeric, median K2P 0.036), Prunus cerasus 5|3, Poa annua 5|2, Arachis
#   hypogaea 5|5. No strength threshold separates true from false sub-splits (true 0.20 and 0.41 of the parent;
#   false 0.19, 0.33, 0.47).
# - Homoeologs from collinear shared insertions (v6 pairing map) to confirm sub-splits: diploids show significant chains
#   too (Rosa odorata 7, Raphanus 9), Camelina only 4. v6 had found homoeolog maps from chains recall 6/21 (B. napus),
#   8/42 (Fragaria); mutual best lineage sharing 9/343.
# - Homoeologs from duplicated insertions (pair K2P below both copies' LTR K2P): MORE frequent within subgenomes than
#   across (B. napus 0.0056 vs 0.0010; Camelina 0.0031 vs 0.0001). Segmental duplicates or LTR-age noise.
# - TE-clade header labels as independent replicates: the Brassica compartment is carried by many clades; fish and
#   fungal libraries are mostly "unknown". Not used.
# - Age classes as extra candidates compared with the parent on another source: fixes Solanum, breaks B. napus.
# - Rescuing the top split from one age class: no gain; A. agadiriana -> k=3.
# - Requiring the parent split to be the best split in the rescuing age class: blocks Avena (its young copies' best top
#   split is A|CD, not the pooled C|AD). Replaced by "strongest split of the parent group on the same copies".
# - Pooled re-partition of the final labels (v6's best_partition with N groups): undoes Solanum (1.00 -> 0.67).
# - From v6's design: LTR-only cross-subgenome flux failed truth; duplicated insertions add nothing to exchange calls;
#   distance-dependent switching no gain; homoeolog maps/ribbons dropped (low recall).
#
# GOTCHAS
# - numba sizes its pool at import: v7 reads -t from sys.argv before importing it (_threads_from_argv). In-process use
#   (tests, imports) keeps the pool of import time; setup() caps numba at it with a warning. argparse abbreviations are
#   off so both parsers agree.
# - numba parallel code is not fork-safe: fork pools (pool_map) must run numba-free functions; compile kernels in the
#   parent before forking workers (else BadZipFile from importlib.metadata in children).
# - list.remove() on dicts holding numpy arrays raises (ambiguous truth value): remove by identity (fixed in -n).
# - `pkill -f <pattern>` kills the calling shell when the pattern is in its own command line: use pgrep + kill PID.
# - Peak RSS summed over a process tree counts the forked workers' shared pages repeatedly; per-process max RSS
#   (/usr/bin/time -v) is the honest number (v6 = v7: 7.3 GB on a 2x Aegilops library).
# - The clade cache (.v7cache in the output dir) is keyed by this file's sha1: any edit rebuilds clades (7-60 s).
# - Pin (copy) the script for benchmark jobs while development continues; launch long jobs detached (setsid nohup) and
#   cap them (1 h per run); run largest libraries first.
#
# TESTS (the test files are not kept; rebuild from this)
# - v6's 132 tests were ported (133 cases) and 13 v7 cases added; all 146 passed on this file (11.7 min). End-to-end
#   tests use a forward simulator of LTR-RT evolution ("simallo"; rewrite it if needed): lineages form a tree with split
#   times, leaves are subgenomes; mergers join leaves at a time. Each TE family starts with 1-3 master elements in the
#   root lineage that mutate and branch; at a split every master and copy is inherited by both daughters at the same
#   coordinates (homoeologues collinear), then evolves independently. Insertions are Poisson per family and lineage
#   (lineage-specific activity, optional bursts); a new copy's two LTRs equal its source's LTR; source = a master of the
#   inserting lineage, or (p 0.3) a copy younger than 100 kyr anywhere in the same nucleus; after a merger a copy lands
#   in subgenome j with probability proportional to target weight x (1 if j is the source's subgenome else phi = 1).
#   Copies are lost with half-life 1.5 Myr; substitutions are Poisson per object (K2P, kappa 2, mu 1.3e-8); optional
#   homoeologous exchanges after the youngest merger. Output: LTRquest-style FASTA + GFF3 + K2P table + .fai and exact
#   truth (copy origins, chromosome subgenomes, homoeologue index). Defaults: 8 chromosomes x 40 Mb, 400 families
#   (LTR length median 600 bp, internal 4.5 kb), activity rate 0.004, 14% of copies written as singletons.
# - Key v7 regression (scenario): root_age 8e6; 5 chromosomes x 6 Mb per subgenome; 80 families; activity rate 0.001,
#   no bursts; no exchanges; tree {split 4.0e6: [{split 1.0e6: [A, B]}, C]}; one merger of A, B, C at 2.0e4 years;
#   seed 1; v7 -t 4 with dating skipped. Expected: 3 subgenomes, all 15 sequences right; the A|B split (root.a) kept
#   on young copies (rule "age") although pooled it is < 0.5x its parent (v6's rule stops at 2). Also true with A/B
#   splits at 0.8 and 1.5 Myr and C at 5 Myr. With A/B at 0.3 Myr the young window is too thin and the near-identical-
#   pair split (pair_splits) finds A|B instead.
# - accept_split unit cases: pooled pass -> "pooled"; weak pooled + young pass + old absent -> "age"; reproduced at
#   every age but < 0.5x -> "present at every age"; top split never rescued; parent strength <= 0 on a source blocks it.
#   find_subgenomes on a planted 3-block matrix ((B,C),A): inferred and -n 3 give the blocks; -n 5 forces 2 flagged
#   splits. Other ported cases: clades vs brute force, Gray-code search vs brute force, link/kNN kernels vs references,
#   HMM vs brute force, header/annotation parsing, outputs and verdicts, report.pdf has 3 pages and its text.
#
# LEADS (not findings)
# - Brassica A05/A06/A09 lineage compartment shared by two independent allopolyploids (B. napus, B. juncea): B. rapa
#   centromeres carry CentBr1 on 8 and CentBr2 on 2 chromosomes (Lim et al. 2005) - which ones was not checked.
# - Old-copy link matrices light up homoeologous chromosome pairs (Solanum report panel E: diagonal lines of bright
#   cells in the between-subgenome blocks): an LTR-only homoeolog map may be possible for close subgenomes.
# - Several "diploid false positives" look like real hybrid-haplotype structure (young-biased, many TE clades).
# - Provenance (not needed to use these notes): independent benchmark LTR_phasing/bench_indep (report.md, results/,
#   truth/), v7 benchmark LTR_phasing/bench_v7 (report.md, memo.md, score/score_v7.py), development harness
#   LTR_phasing/dev_v7, v6 subphaser_LTR/v6 (design_v6.md, VALIDATION.md, tests/simallo.py).
