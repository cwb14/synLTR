#!/usr/bin/env python3
"""v6.py - subgenome phasing, structure, dating and ancestry of polyploid genomes from their LTR-retrotransposon
library alone (no genome sequence, no progenitors).

Inputs: LTR-RT FASTA (headers chrom:start-end[#Class/Superfamily/Clade]; also chrom:start..end, name::chrom:start-end
and strand suffixes), an LTR annotation with each copy's 5'-3' LTR divergence (LTRquest GFF3, or a K2P table), and
optionally the genome .fai (sequence lengths).

Phasing: every canonical 21-mer of the library is hashed exactly; k-mers carried by 2-250 copies on >= 2 sequences
define lineages (distinct carrier sets = clades). Sequences are linked by clades (1/(m-1) per pair, each clade once,
same-sequence copies never link) and by nearest-relative votes, compared with a degree-corrected null and split by
maximum within-group excess (exact for two groups). Subgenomes are found by divisive splitting, a split being kept
when disjoint half-libraries reproduce it; support = agreement over 100 half-libraries.

usage: v6.py --ltr_fasta LTR.fa --annotation ANNOT.gff3 [--fai GENOME.fai] --outdir OUT [-t 8]
"""
from __future__ import annotations

import os

_THREADS0 = int(os.environ.get("NUMBA_NUM_THREADS") or os.environ.get("OMP_NUM_THREADS") or 1)   # default share (no -t)
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):   # BLAS single-threaded whatever the
    os.environ[_v] = "1"                                                       # shell exports: -t is the true total
os.environ["NUMBA_NUM_THREADS"] = str(os.cpu_count() or 1)                    # numba may use every core; -t picks the share

import argparse
import gzip
import hashlib
import json
import math
import re
import socket
import sys
import time
from pathlib import Path

import numpy as np

VERSION = "6.0.0-dev"
K = 21               # k-mer length: unique in plant genomes yet short enough to span one SNP per clade variant (v4/v5)
CAP = 250            # largest clade (carrier copies) used for links (v5 chose this cap in 42 of 49 runs)
NMAX_KNN = 50        # relatives for nearest-relative votes come from clades of <= 50 copies (v5 benchmark)
K_NEIGH = 3          # each copy votes for the 3 other sequences holding its closest relatives (v5 benchmark)
MIN_ELEMENTS = 10    # sequences with fewer LTR-RTs are not phased (cannot be assessed; v4/v5)
MIN_REP = 0.8        # a split is kept when half-libraries reproduce it (adjusted Rand; prediction-strength convention)
MIN_REL = 0.5        # rule R0 (v5): a sub-split must be >= this fraction as strong as its parent (heuristic, from v4)
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
FIG_LET = 0.42       # panel letters sit this far (inches) left of each panel's upper-left corner
FIG_MAX_SEQ = 60     # figures draw the 60 longest phased sequences (more cannot be labelled; tables hold all)
FIG_EXAMPLES = 300   # candidates aligned to pick the typical worked example (a random 300 bound the time, ~1 s)
SG_COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9",   # subgenome hues: Okabe-Ito
             "#332288", "#117733", "#882255", "#999933", "#AA4499", "#44AA99"]  # (no yellow), then Tol muted
FIG_GRAY, FIG_LIGHT, FIG_DARK = "#9a9a9a", "#d9d9d9", "#4d4d4d"   # context: copies without the story's property
FIG_ACCENT = "#000000"  # the story's entity in the dating figure (pairs across subgenomes, inherited insertions)
T0 = time.time()
VERBOSE = False


def log(msg):
    print(f"[v6 {time.time() - T0:7.1f}s] {msg}", file=sys.stderr, flush=True)


def vlog(msg):
    if VERBOSE:
        log(msg)


def die(msg):
    log("ERROR: " + msg)
    sys.exit(1)


def natural_key(s):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", str(s))]


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description="Subgenome phasing, structure, dating and ancestry from an LTR-RT library.")
    ap.add_argument("--ltr_fasta", nargs="+", required=True, help="LTR-RT FASTA(.gz) file(s)")
    ap.add_argument("--annotation", nargs="+", required=True,
                    help="LTR annotation with each copy's LTR divergence: LTRquest GFF3 or K2P table (TSV)")
    ap.add_argument("--fai", nargs="+", default=None, help="genome .fai (sequence lengths; optional)")
    ap.add_argument("--outdir", required=True, help="output directory")
    ap.add_argument("-n", "--n_subgenomes", type=int, default=None, help="number of subgenomes (default: inferred)")
    ap.add_argument("--config", default=None, help="optional hint: sets of sequences known to be in different subgenomes")
    ap.add_argument("--assign", default=None, help="fixed partition instead of phasing: 'sequence label' per line (whitespace-separated; "
                    "labels become SG1, SG2... in file order, mapping in summary.json)")
    ap.add_argument("--features", default=None, help="BED/GFF3 features to give the ancestry painted beneath them")
    ap.add_argument("--mu", type=float, default=None, help="substitutions/site/year for ages in years (default 1.3e-8, flagged)")
    ap.add_argument("--label", default=None, help="species name for figure legends and report.md (default: FASTA name)")
    ap.add_argument("-t", "--threads", type=int, default=4, help="total threads")
    ap.add_argument("--seed", type=int, default=1, help="random seed")
    ap.add_argument("--no_plots", action="store_true", help="skip figures")
    ap.add_argument("-v", "--verbose", action="store_true", help="per-step progress and sanity checks")
    ap.add_argument("--dev_scale", type=int, default=0, help=argparse.SUPPRESS)      # k-mer sampling 1/scale (0 = budget)
    ap.add_argument("--dev_rule", default="R0", help=argparse.SUPPRESS)              # split acceptance rule
    ap.add_argument("--dev_split_rule", default="confinement", help=argparse.SUPPRESS)  # reported split estimator order
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
        if a.n_subgenomes or a.config:
            die("--assign fixes the partition: -n and --config do not apply")
    if a.config:
        if not Path(a.config).exists():
            die(f"input not found: {a.config}")
        cfg = read_config(a.config)
        if not a.n_subgenomes and max(len(x) for x in cfg) < 2:
            die("config sets hold one sequence each: they define no cannot-link constraint (N would be 1)")
    import numba
    if a.threads > numba.config.NUMBA_NUM_THREADS:
        log(f"WARNING: -t {a.threads} exceeds the {numba.config.NUMBA_NUM_THREADS} threads available; using those")
    numba.set_num_threads(min(a.threads, numba.config.NUMBA_NUM_THREADS))
    Path(a.outdir).mkdir(parents=True, exist_ok=True)
    (Path(a.outdir) / ".v6cache").mkdir(exist_ok=True)


SOURCE_SHA = hashlib.sha1(Path(__file__).read_bytes()).hexdigest()[:12]   # code changes invalidate caches


def cache_path(a, stage, files, params):
    """Cache file keyed by v6 version and source, stage, input files (path, size, mtime) and parameters."""
    h = hashlib.sha1(f"{VERSION}|{SOURCE_SHA}|{stage}|{params}".encode())
    for f in files:
        st = os.stat(f)
        h.update(f"{os.path.abspath(f)}|{st.st_size}|{int(st.st_mtime)}".encode())
    return Path(a.outdir) / ".v6cache" / f"{stage}_{h.hexdigest()[:16]}.npz"


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
                                family=at.get("family", ""), clade=at.get("clade", ""),
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
                            family=r.get("family", ""), clade="", superfamily="", depth=float("nan"),
                            strand=r.get("strand", ".") if r.get("strand", ".") in ("+", "-") else ".")
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


def half_pairs(lm, pairs, rng, threads=None):
    """`pairs` complementary random halves of the library -> list of (W_a, W_b) link matrices. Computed in this
    process with the parallel link kernel: numba's OpenMP layer is not fork-safe, so no process pool here."""
    ph = lm.ek >= 0
    out = []
    for seed in rng.integers(0, 2 ** 62, pairs):
        u = np.random.default_rng(int(seed)).random(lm.n) < 0.5
        out.append((lm.matrix(u & ph), lm.matrix(~u & ph)))
    return out


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
    Wg = W[np.ix_(idx, idx)]
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
    two halves of each pair, and its standard error."""
    v = np.array([ari(bipartition(Wa, g, rng, sets), bipartition(Wb, g, rng, sets)) for Wa, Wb in pairs])
    return float(v.mean()), (float(v.std(ddof=1) / math.sqrt(v.size)) if v.size > 1 else 0.0)


def find_subgenomes(W, pairs, rng, rule="R0", min_size=2, max_n=MAX_N, more=None):
    """Divisive search for subgenome structure: a group is split when the best splits found independently in
    disjoint halves of the library agree (mean adjusted Rand >= MIN_REP) and the rule accepts it (R0: below the
    top, >= MIN_REL as strong as the parent split). When the mean is within 2 standard errors of MIN_REP, more
    half-library pairs are drawn (more(n) -> n pairs; up to MAX_PAIRS) and the test is flagged borderline if it still
    straddles. N=1 is a possible answer. Returns labels and every tested split."""
    C = W.shape[0]
    groups = [(np.arange(C), None, "root")]
    final, tests = [], []
    while groups:
        g, parent, path = groups.pop(0)
        if len(g) < 2 * min_size or len(final) + len(groups) + 1 >= max_n:
            final.append(g)
            continue
        z = bipartition(W, g, rng)
        if min(np.bincount(z, minlength=2)) < min_size:
            tests.append(dict(members=g, split=z, rep=0.0, se=0.0, n_pairs=0, borderline=False, contrast=0.0,
                              parent_contrast=parent, accepted=False, path=path, rule=rule))
            final.append(g)
            continue
        con = split_contrast(W, g, z)
        rep, se = rep_stats(pairs, g, rng)
        while more is not None and len(pairs) < MAX_PAIRS and abs(rep - MIN_REP) < 2 * se:
            pairs.extend(more(HALF_PAIRS))
            rep, se = rep_stats(pairs, g, rng)
        if rule == "R0":
            ok = rep >= MIN_REP and con > 0 and (parent is None or con >= MIN_REL * parent)
        else:
            die(f"unknown split rule {rule}")
        tests.append(dict(members=g, split=z, rep=rep, se=se, n_pairs=len(pairs), borderline=abs(rep - MIN_REP) < 2 * se,
                          contrast=con, parent_contrast=parent, accepted=ok, path=path, rule=rule))
        if ok:
            groups += [(g[z == 0], con, path + ".a"), (g[z == 1], con, path + ".b")]
        else:
            final.append(g)
    lab = np.zeros(C, int)
    for k, g in enumerate(sorted(final, key=lambda x: x.min())):
        lab[g] = k
    return lab, tests


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
from collections import defaultdict  # noqa: E402

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


def pair_splits(W, lab, cl_pairs, pairs, rng, more=None):
    """A final group holding >= 2 paired sequence pairs cannot be one subgenome (a subgenome carries one copy of
    each chromosome): split it with every such pair separated, on links without the paired links, and keep the
    split when half-libraries reproduce it (adjusted Rand >= MIN_REP, more pairs when borderline). Returns the new
    labels and the tests."""
    lab, tests, k_next = lab.copy(), [], int(lab.max()) + 1
    for k in range(int(lab.max()) + 1):
        g = np.flatnonzero(lab == k)
        inside = [(x, y) for x, y in cl_pairs if lab[x] == k and lab[y] == k]
        if len(inside) < 2:
            continue
        pos = {x: i for i, x in enumerate(g)}
        sets = [[pos[x], pos[y]] for x, y in inside]
        Wm = mask_pairs(W, inside)
        hm = [(mask_pairs(a, inside), mask_pairs(b, inside)) for a, b in pairs]
        z = bipartition(Wm, g, rng, sets)
        rep, se = rep_stats(hm, g, rng, sets)
        while more is not None and len(hm) < MAX_PAIRS and abs(rep - MIN_REP) < 2 * se:
            new = more(HALF_PAIRS)
            pairs.extend(new)
            hm += [(mask_pairs(a, inside), mask_pairs(b, inside)) for a, b in new]
            rep, se = rep_stats(hm, g, rng, sets)
        ok = rep >= MIN_REP
        tests.append(dict(members=g, split=z, rep=rep, se=se, n_pairs=len(hm), borderline=abs(rep - MIN_REP) < 2 * se,
                          contrast=split_contrast(Wm, g, z), parent_contrast=None, accepted=ok, path=f"paired{k}",
                          rule="pairs", inside=inside))
        if ok:
            lab[g[z == 1]] = k_next
            k_next += 1
    return lab, tests


# ------------------------------------------------------------------------------------------------ phasing pipeline
def verdict_text(N, inside, big, tests, source):
    """One-line verdict. inside: paired sequence pairs left in one final group; big: groups of >= 3 near-identical
    sequences (paired map components)."""
    if N >= 2:
        v = "subgenomes found" + (" (borderline)" if any(t["accepted"] and t["borderline"] for t in tests) else "")
        if inside >= 2:
            v += f"; {inside} pairs of near-identical sequences share a group: more subgenomes than resolved"
        if big:
            v += f"; {len(big)} groups of 3 or more near-identical sequences"
        return v
    if inside >= 2 or big:
        return (f"cannot tell: near-identical copies of {inside + len(big)} chromosomes (homologues of a young "
                f"polyploid, or both haplotypes of a diploid assembly); LTR-RT lineages do not separate them")
    if source == "given" or (tests and tests[0]["n_pairs"] and not tests[0]["borderline"]):
        return "no structure detected"
    return "cannot tell"


def group_label(g):
    return f"SG{g + 1}"


def read_config(path):
    sets = []
    with open(path) as f:
        for line in f:
            line = line.split("#", 1)[0].split()
            if line:
                sets.append(line)
    if not sets:
        die(f"no sequence sets in config {path}")
    return sets


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


def run_phasing(a, lib, cl):
    rng = np.random.default_rng(a.seed)
    n_el_all = np.bincount(lib.seq, minlength=len(lib.seqnames))
    user = read_assign(a.assign, lib.seqnames) if a.assign else None      # fixed partition: no structure search
    phased = [i for i in range(len(lib.seqnames)) if n_el_all[i] >= MIN_ELEMENTS and
              (user is None or lib.seqnames[i] in user)]
    if user and len(phased) < len(user):
        log(f"WARNING: --assign: {len(user) - len(phased)} sequence(s) carry < {MIN_ELEMENTS} LTR-RTs (unphased)")
    cfg = read_config(a.config) if a.config else None
    N = a.n_subgenomes or (max(len(s) for s in cfg) if cfg else None)
    if cfg and not a.n_subgenomes and N < 2:
        die("config sets hold one sequence each: they define no cannot-link constraint (N would be 1)")
    if len(phased) < max(N or 2, 2):
        die(f"only {len(phased)} sequences carry >= {MIN_ELEMENTS} LTR-RTs; cannot form {max(N or 2, 2)} groups")
    names = [lib.seqnames[i] for i in phased]
    if cfg:
        miss = sorted({x for st in cfg for x in st} - set(names), key=natural_key)
        if miss:
            log(f"WARNING: {len(miss)} config sequence(s) absent or with < {MIN_ELEMENTS} LTR-RTs (ignored): "
                + ", ".join(miss[:10]) + (" ..." if len(miss) > 10 else ""))
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
    W = W_raw = lm.matrix(ek >= 0)
    pairs = half_pairs(lm, HALF_PAIRS, rng, a.threads)
    rng_more = np.random.default_rng([a.seed, 5])
    more = lambda n: half_pairs(lm, n, rng_more, a.threads)
    use_pairs = []
    n_el = np.bincount(ek[ek >= 0], minlength=C)
    tree, init, source, tests = [], None, "given" if N else "inferred", []
    if user:
        here = {user[s] for s in names}
        order = [g for g in dict.fromkeys(user.values()) if g in here]    # labels in the file's order
        init, N, source = np.array([order.index(user[s]) for s in names]), len(order), "user"
        log(f"partition given by the user: {N} group(s) ({', '.join(order)} -> "
            f"{', '.join(group_label(g) for g in range(N))})")
    elif N is None:
        lab, tests = find_subgenomes(W, pairs, rng, rule=a.dev_rule, more=more)
        if not cfg and cl_pairs:
            lab, ptests = pair_splits(W, lab, cl_pairs, pairs, rng, more)
            tests += ptests
            use_pairs = [p for t in ptests if t["accepted"] for p in t["inside"]]
            if use_pairs:
                cfg = [[names[x], names[y]] for x, y in use_pairs]
                W = mask_pairs(W, use_pairs)
                pairs[:] = [(mask_pairs(x, use_pairs), mask_pairs(y, use_pairs)) for x, y in pairs]
                log(f"paired sequences: {len(use_pairs)} pairs kept in different groups")
        for t in tests:
            tree.append(dict(level=t["path"], part_a=[names[i] for i in t["members"][t["split"] == 0]],
                             part_b=[names[i] for i in t["members"][t["split"] == 1]], replicability=t["rep"],
                             replicability_se=t["se"], n_pairs=t["n_pairs"], borderline=bool(t["borderline"]),
                             strength=t["contrast"], parent_strength=t["parent_contrast"], rule=t["rule"],
                             accepted=bool(t["accepted"])))
            log(f"split {len(tree[-1]['part_a'])}|{len(tree[-1]['part_b'])}: replicability {t['rep']:.2f} "
                f"(se {t['se']:.3f}, {t['n_pairs']} pairs{', borderline' if t['borderline'] else ''}), strength "
                f"{t['contrast']:.2f} -> {'accepted' if t['accepted'] else 'rejected'}")
        N, init = int(lab.max()) + 1, lab
        log(f"inferred number of subgenomes: {N}")
    elif not cfg and cl_pairs:                                       # forced N: paired sequences kept apart too
        use_pairs = list(cl_pairs)
        cfg = [[names[x], names[y]] for x, y in use_pairs]
        W = mask_pairs(W, use_pairs)
        pairs[:] = [(mask_pairs(x, use_pairs), mask_pairs(y, use_pairs)) for x, y in pairs]
        log(f"paired sequences: {len(use_pairs)} pairs kept in different groups")
    halves = [w for p in pairs for w in p]
    one = N == 1
    if one:
        z = np.zeros(C, int)
    elif user:
        sets = make_sets(names, cfg, N)
        z = init
    else:
        sets = make_sets(names, cfg, N)
        z = canon_labels(best_partition(W - expected(W), N, sets, rng, init=init), n_el)
        z, realized = realized_groups(z, N)
        if realized < N:
            log(f"WARNING: only {realized} of the {N} requested groups are non-empty; {realized} reported")
            N = realized
    inside = int(sum(z[x] == z[y] for x, y in cl_pairs))            # paired sequences left in one final group
    verdict = "partition given by the user" if user else verdict_text(N, inside, big, tests, source)
    log(f"verdict: {verdict}")
    # W (paired links masked) drives the partition, support and per-sequence tests; links_oe.tsv reports W_raw
    user_labels = {group_label(g): order[g] for g in range(len(order))} if user else None
    common = dict(names=names, C=C, ek=ek, tree=tree, W=W, W_raw=W_raw, n_el=n_el, n_inf=n_inf, source=source,
                  user_labels=user_labels,
                  clades=lm.cl, clades_all=cl,
                  N_requested=N_requested or N,
                  n_clades=lm.cl.n, pairing=pairing, cl_pairs=use_pairs, paired_groups=big, verdict=verdict,
                  n_pairs=len(pairs))
    if one:
        partner, psup, pex = partner_pairs(W, halves, z)
        return dict(common, z=z, N=1, support=np.full(C, np.nan), partners=(partner, psup, pex),
                    own_oe=np.full(C, np.nan), oth_oe=np.full(C, np.nan), pval=np.full(C, np.nan))
    support = support_from_halves(halves, z, N, sets, rng)
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
    return dict(P=P[ok], d=d[ok], se=se[ok], stratum=st[ok])


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
def _clock_profile(Su, Sv, y, cols):
    """Log-likelihood, phi, eps, converged at each candidate merger AGE_GRID[cols]: P(younger < T) = 1 - Su Sv."""
    out = np.empty((cols.size, 4))
    for c in prange(cols.size):
        j = cols[c]
        B = 1.0 - Su[:, j] * Sv[:, j]
        ll, phi, eps, ok = _phi_eps_fit(B, y)
        out[c, 0] = ll
        out[c, 1] = phi
        out[c, 2] = eps
        out[c, 3] = 1.0 if ok else 0.0
    return out


def _landing_bound(Su, Sv, y, cols, phi, clusters):
    """Upper bound on the merger (LTR divergence units) when the clock sees no excess of moves, IF new copies land in
    either subgenome in proportion to its size: phi_i = the other subgenome's genome share seen from pair i's older
    copy. A merger at T predicts P(split) = phi_i P(younger < T) + eps (1 - P(younger < T)), eps = the observed split
    share; per lineage-cluster bootstrap draw (own generator: the dating's draws stay untouched) the oldest T within
    1.92 log-likelihood units of the best (T = 0 included); bound = 97.5% quantile widened by the merger model error
    (merger_clock floors it at one substitution per LTR); None when even the oldest candidate survives in > 2.5% of
    draws. Fails when transposition stays within each subgenome after the merger (stated wherever it is reported)."""
    eps = min(max(float(y.mean()), 1e-6), 1 - 1e-6)
    Bm = 1.0 - Su[:, cols] * Sv[:, cols]
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


def merger_clock(k_u, n_u, k_v, n_v, cross, clusters, g, rng, t_max, phi_land=None):
    """Pairwise sibling clock: a sibling pair sits in two subgenomes only when its younger copy was inserted after the
    merger (eps: exchange or misassignment noise, age-independent). T (LTR divergence units) by profile likelihood
    over AGE_GRID below t_max; lineage-cluster bootstrap widened by MODEL_ERR["merger"]. Not datable unless cross pairs are >= 2-fold
    enriched among post-merger pairs (phi >= 2 eps) in >= 97.5% of resamples; an upper bound only ("younger than")
    when the interval reaches the divergence of one substitution per LTR (1 / median LTR sites)."""
    def surv(k, n):
        post = _age_lik(k, n) * g[None]
        post /= np.maximum(post.sum(1, keepdims=True), 1e-300)
        return np.ascontiguousarray(np.cumsum(post[:, ::-1], 1)[:, ::-1])          # P(age >= AGE_GRID[j])
    Su, Sv = surv(k_u, n_u), surv(k_v, n_v)
    y = cross.astype(np.float64)
    cols = np.flatnonzero((AGE_GRID > 0) & (AGE_GRID < t_max)).astype(np.int64)
    if cols.size < 3 or y.size < 50:
        return None
    prof = _clock_profile(Su, Sv, y, cols)
    i = int(np.argmax(prof[:, 0]))
    T, phi, eps = AGE_GRID[cols[i]], prof[i, 1], prof[i, 2]
    Tb, gap, bad = [], [], int((prof[:, 3] == 0).sum())
    for sel in _cluster_boot(clusters, rng, BOOT):
        pb = _clock_profile(Su[sel], Sv[sel], y[sel], cols)
        bad += int((pb[:, 3] == 0).sum())
        j = int(np.argmax(pb[:, 0]))
        Tb.append(AGE_GRID[cols[j]])
        gap.append(pb[j, 1] - 2 * pb[j, 2])
    lo, hi = widen(T, *np.quantile(Tb, [0.025, 0.975]), "merger")
    t_res = 1.0 / np.median(np.concatenate([n_u, n_v]))
    bound = bool(np.quantile(Tb, 0.025) <= t_res)              # resolution judged on the bootstrap alone
    datable = bool(phi >= 2 * eps and np.quantile(gap, 0.025) > 0)
    if bad:
        log(f"WARNING: the merger clock's rate fit did not converge in {bad} of {cols.size * (BOOT + 1)} fits")
    lb = _landing_bound(Su, Sv, y, cols, np.asarray(phi_land, float), clusters) if (
        phi_land is not None and not datable and float(np.mean(phi_land)) >= 2 * float(y.mean())) else None
    if lb is not None:                                     # (random landing must predict >= 2x the background,
        lb = max(lb, t_res)                                #  the clock's rule; resolution floor: one substitution
        lb = lb if lb < t_max else None                    #  per LTR; a bound at or above the split says nothing)
    return dict(landing_bound=lb, estimate=float(T), lo=0.0 if bound else float(lo), hi=float(hi), phi=float(phi),
                eps=float(eps),
                datable=datable, bound=bound, n=int(y.size), n_cross=int(y.sum()))


# ------------------------------------------------------------------------------------------------ dating: pipeline
SPLIT_ORDER = {"orthologs": ("orthologs", "confinement", "profile"), "confinement": ("confinement", "orthologs", "profile")}


def date_pair(lib, pt, loc, x, y, ancestry, rng, split_rule, sizes=None):
    """Split and merger of subgenomes x and y (group indices) from the copy-pair table pt and copy locations loc.
    K2P units (pair divergence = 2 mu T; LTR divergence of a copy = 2 mu age)."""
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
    rep = next((k for k in SPLIT_ORDER[split_rule] if out["split"][k]), None)
    out["split_reported"] = rep
    ests = [v for v in out["split"].values() if v and np.isfinite(v["lo"])]
    out["split_agreement"] = ("none" if not ests else "single" if len(ests) == 1 else
                              "agree" if max(v["lo"] for v in ests) <= min(v["hi"] for v in ests) else "disagree")
    # merger: sibling pairs (MRCA = older copy's insertion) younger than the confined plateau's end and the split
    ks = k2p_se_copy(lib)
    ka, kb = lib.k2p[Pm[:, 0]], lib.k2p[Pm[:, 1]]
    old = np.where(ka >= kb, Pm[:, 0], Pm[:, 1])
    with np.errstate(invalid="ignore"):
        sib = dm <= np.fmax(ka, kb) + 2 * np.sqrt(ks[old] ** 2 + sm ** 2)
    t_max = plateau_end(dm, cross)
    if rep:
        t_max = min(t_max, out["split"][rep]["estimate"])
    sel = ~dup & sib & (dm < t_max)
    kc, nc = ltr_counts(lib)
    good = sel & np.isfinite(kc[Pm[:, 0]]) & np.isfinite(kc[Pm[:, 1]]) & (nc[Pm[:, 0]] > 0) & (nc[Pm[:, 1]] > 0)
    out["n_sibling_pairs"] = int(good.sum())
    out["ev"] = dict(rows=np.flatnonzero(m0), orth=orth, clock=good)        # the evidence drawn by fig_dating
    pool = np.flatnonzero(np.isin(loc, [x, y]) & np.isfinite(kc) & (nc > 0))
    if good.sum() >= 50 and pool.size >= 50:
        pool = rng.choice(pool, min(pool.size, 20000), replace=False)
        g = age_prior(kc[pool], nc[pool])
        u, v = Pm[good, 0], Pm[good, 1]
        phi_land = None
        if sizes is not None:                              # random landing: the other subgenome's genome share
            og = loc[old[good]]
            phi_land = np.where(og == x, sizes[y], sizes[x]) / (sizes[x] + sizes[y])
        out["merger"] = merger_clock(kc[u], nc[u], kc[v], nc[v], cross[good] > 0, clusters[good], g, rng, t_max,
                                     phi_land=phi_land)
    return out


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
    if mg["bound"]:
        return f"cross-subgenome transposition began < {kyr(mg['hi'])} ago (a bound on its onset, not on the merger)"
    return (f"merger {kyr(mg['estimate'])} ({kyr(mg['lo'])}-{kyr(mg['hi'])}; onset of cross-subgenome "
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
    res = []
    for x in range(ph["N"]):
        for y in range(x + 1, ph["N"]):
            res.append(date_pair(lib, pt, loc, x, y, ancestry, rng, a.dev_split_rule, sizes=sizes))
    mu = a.mu or MU_DEFAULT
    for r in res:
        x, y = (group_label(g) for g in r["groups"])
        log(f"{x}-{y}: {split_text(r, mu)}; {merger_text(r['merger'], r.get('n_sibling_pairs', 0), mu)} "
            f"[mu={mu:g}{'' if a.mu else ', default'}]")
    return dict(mu=mu, mu_given=bool(a.mu), pairs=res, n_copy_pairs=int(len(pt["P"])), ancestry=ancestry, pt=pt)


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
                                      if_random_landing_below_years=yr(mg["landing_bound"]) if mg and
                                      mg.get("landing_bound") else None)))
    return dict(mu=dt["mu"], mu_given=dt["mu_given"], n_copy_pairs=dt["n_copy_pairs"], pairs=pairs)


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
            st = "none" if mg is None else ("not_datable" if not mg["datable"] else ("bound" if mg["bound"] else "dated"))
            vals = ([None if mg["bound"] else mg["estimate"], mg["lo"], mg["hi"]] if mg and mg["datable"]
                    else [None] * 3)                       # a bound has no point estimate
            f.write("\t".join([x, y, "merger", "clock", "yes", st] + [_f(z_) for z_ in vals] +
                               [_f(z_ / (2 * dt["mu"]) if z_ is not None else None) for z_ in vals] +
                               [str(mg["n"]) if mg else "0"]) + "\n")
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
            f.write(f"{x['seq']}\tv6\tregion\t{x['start']}\t{x['end']}\t{x['llr'] / math.log(10):.2f}\t.\t.\t"
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


# ------------------------------------------------------------------------------------------------ figures and report
def _fig_style():
    """matplotlib (imported here: --no_plots runs never load it) with Arial if installed, else Liberation Sans from
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
                         "axes.labelsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
                         "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                         "xtick.major.size": 2.5, "ytick.major.size": 2.5, "pdf.fonttype": 42, "ps.fonttype": 42,
                         "axes.spines.top": False, "axes.spines.right": False, "mathtext.default": "regular",
                         "hatch.linewidth": 0.4, "savefig.dpi": 300})
    return plt, family


def _ax(fig, x0, y0, w, h):
    """Axes at x0, y0 (inches from the lower left), w x h inches: every panel sits on a grid measured in inches."""
    W, H = fig.get_size_inches()
    return fig.add_axes([x0 / W, y0 / H, w / W, h / H])


def _letter(fig, ax, s):
    """Bold panel letter at one fixed offset from the axes' upper-left corner (FIG_LET inches left, 0.05 above)."""
    W, H = fig.get_size_inches()
    bb = ax.get_position()
    fig.text(bb.x0 - FIG_LET / W, bb.y1 + 0.05 / H, s, fontsize=9, fontweight="bold", va="bottom", ha="left")


def _key(fig, x0, y0, w, items, draw=True):
    """Color/symbol key on the figure from x0, y0 (inches, row centre), w wide, wrapping to further rows 0.16 in
    lower: items (kind, color, label), kind one of tick, box, hatch, outline, dot, diamond -> rows used (draw=False:
    measure only)."""
    from matplotlib.lines import Line2D
    from matplotlib.patches import Rectangle
    W, H = fig.get_size_inches()
    r = fig.canvas.get_renderer()
    x, y, g, rows = x0, y0, 0.14, 1
    for kind, col, lab in items:
        t = fig.text(0, 0, lab, fontsize=7)
        tw = t.get_window_extent(r).width / fig.dpi
        t.remove()
        if x > x0 and x + g + 0.05 + tw > x0 + w:
            x, y, rows = x0, y - 0.16, rows + 1
        if not draw:
            x += g + 0.05 + tw + 0.28
            continue
        if kind == "tick":
            fig.add_artist(Line2D([x / W, (x + g) / W], [y / H] * 2, color=col, lw=1.2))
        elif kind in ("dot", "diamond"):
            fig.add_artist(Line2D([(x + g / 2) / W], [y / H], marker="o" if kind == "dot" else "D", ms=3.2, color=col,
                                  lw=0))
        else:
            fig.add_artist(Rectangle((x / W, (y - 0.045) / H), g / W, 0.09 / H, transform=fig.transFigure,
                                     facecolor="white" if kind == "hatch" else "none" if kind == "outline" else col,
                                     edgecolor=col if kind != "box" else "none",
                                     hatch="//////" if kind == "hatch" else None,
                                     lw=0.8 if kind == "outline" else 0))
        fig.text((x + g + 0.05) / W, y / H, lab, va="center", ha="left", fontsize=7)
        x += g + 0.05 + tw + 0.28
    return rows


_MARK = re.compile(r"^([(\[\"']*)(\*\*|\*)?(.*?)(\*\*|\*)?([,.;:!?)\]\"']*)$")


def _tokens(markup):
    """'**bold**' and '*italic*' markup -> [(word, weight, style)]; markers may span several words and sit inside
    punctuation ('(*A. suecica*),')."""
    out, bold, ital = [], False, False
    for w in markup.split():
        pre, op, core, cl, post = _MARK.match(w).groups()
        if op:
            bold, ital = (not bold, ital) if op == "**" else (bold, not ital)
        word = pre + core + post
        if word:
            out.append((word, "bold" if bold else "normal", "italic" if ital else "normal"))
        if cl:
            bold, ital = (not bold, ital) if cl == "**" else (bold, not ital)
    return out


def _legend_page(plt, pdf, paragraphs, width_in, family, size=10):
    """One legend page as wide as the figure, height fitted to the text: ragged-right 10-pt paragraphs with *italic*
    and **bold** markup, laid out word by word (text stays text)."""
    lh, margin = size * 1.35 / 72, 0.35
    probe = plt.figure(figsize=(width_in, 2), dpi=720)                 # fine pixels: word widths without rounding
    r = probe.canvas.get_renderer()
    cache = {}

    def wid(word, wt, st):
        if (word, wt, st) not in cache:
            t = probe.text(0, 0, word, fontsize=size, fontweight=wt, fontstyle=st, family=family)
            cache[(word, wt, st)] = t.get_window_extent(r).width / probe.dpi
            t.remove()
        return cache[(word, wt, st)]
    space = wid("a a", "normal", "normal") - 2 * wid("a", "normal", "normal")
    avail = width_in - 2 * margin

    def pieces(tok):                                            # a word wider than the page: break it, keep it all
        word, wt, st = tok
        if wid(word, wt, st) <= avail:
            return [tok]
        out, cur = [], ""
        for chunk in re.split(r"(?<=[/_\-.])", word):
            if cur and wid(cur + chunk, wt, st) > avail:
                out.append(cur)
                cur = ""
            while wid(chunk, wt, st) > avail:
                k = len(chunk)
                while k > 1 and wid(chunk[:k], wt, st) > avail:
                    k -= 1
                out.append(chunk[:k])
                chunk = chunk[k:]
            cur += chunk
        return [(x_, wt, st) for x_ in out + ([cur] if cur else [])]
    laid = []
    for p in paragraphs:
        line, x = [], 0.0
        for tok in (pc for t_ in _tokens(p) for pc in pieces(t_)):
            w = wid(*tok)
            if line and x + space + w > avail:
                laid.append(line)
                line, x = [], 0.0
            line.append((tok, x + space if line else 0.0))
            x = line[-1][1] + w
        laid += [line, None]
    plt.close(probe)
    h = 2 * margin + lh * len(laid)
    fig = plt.figure(figsize=(width_in, h))
    y = h - margin
    for line in laid:
        for (word, wt, st), x in (line or []):
            fig.text((margin + x) / width_in, y / h, word, fontsize=size, fontweight=wt, fontstyle=st, family=family,
                     va="top", ha="left")
        y -= lh
    pdf.savefig(fig)
    plt.close(fig)


def _save_figure(plt, fig, stem, short, long_, family):
    """<stem>.png (figure only, 300 dpi) and <stem>.pdf: figure, short legend, long legend (written from the run's
    numbers; the last line names the run)."""
    from matplotlib.backends.backend_pdf import PdfPages
    fig.savefig(f"{stem}.png", dpi=300)
    with PdfPages(f"{stem}.pdf") as pdf:
        pdf.savefig(fig)
        w = fig.get_size_inches()[0]
        _legend_page(plt, pdf, short, w, family)
        _legend_page(plt, pdf, long_ + [f"Analysis: {socket.gethostname()}:{Path(stem).parent.resolve()} · Script: "
                                        f"{Path(__file__).resolve()} · {time.strftime('%Y-%m-%d')}"], w, family)
    plt.close(fig)


def _yr_label(v):
    return f"{v / 1e6:g} Myr" if v >= 1e6 else f"{v / 1e3:g} kyr" if v >= 1e3 else f"{v:g} yr"


def _log_axis(ax, axis, unit):
    """Log axis labelled in plain units at each power of ten (and x3 when fewer than three powers show): 'years'
    (1 kyr, 10 kyr, 1 Myr) or plain numbers ('Mb': 0.1, 1, 10)."""
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter, NullLocator
    lo, hi = ax.get_xlim() if axis == "x" else ax.get_ylim()
    lo, hi = min(lo, hi), max(lo, hi)
    if lo <= 0:
        return
    e0, e1 = int(np.floor(np.log10(lo))), int(np.ceil(np.log10(hi)))
    for series in ((1,), (1, 3), (1, 2, 5)):                          # the coarsest series giving >= 3 ticks
        ticks = [m * 10.0 ** e for e in range(e0 - 1, e1 + 1) for m in series if lo <= m * 10.0 ** e <= hi]
        if len(ticks) >= 3:
            break
    a = ax.xaxis if axis == "x" else ax.yaxis
    a.set_major_locator(FixedLocator(ticks))
    a.set_major_formatter(FuncFormatter(lambda v, _: _yr_label(v) if unit == "years" else f"{v:,g}"))
    a.set_minor_locator(NullLocator())
    a.set_minor_formatter(NullFormatter())


def _rate(mu):
    """Plain-language substitution rate: 7e-9 -> '7 substitutions per site per billion years'."""
    return f"{mu * 1e9:.3g} substitutions per site per billion years"


def _rows(st, en, gap):
    """IGV 'expanded' rows: each element goes to the first row whose last element ended more than gap before it."""
    ends, row = [], np.zeros(st.size, int)
    for i in np.argsort(st, kind="stable"):
        for r, e in enumerate(ends):
            if e + gap < st[i]:
                ends[r], row[i] = en[i], r
                break
        else:
            ends.append(en[i])
            row[i] = len(ends) - 1
    return row


def _n(k, word):
    """A count with its noun, plural unless k == 1: '1 inherited insertion', '2 inherited insertions'."""
    return f"{k:,} {word}" + ("" if k == 1 else "s")


def _short_name(s, keep=30):
    """Sequence names longer than keep characters shortened in the middle for figure labels (tables keep them)."""
    return s if len(s) <= keep else s[:keep // 2 - 2] + "…" + s[-(keep // 2 - 1):]


_ESTIMATORS = (("orthologs", "inherited", "inherited insertions (their median)"),
               ("confinement", "confinement", "confinement (the age at which lineages stop being confined to one "
                                              "subgenome)"),
               ("profile", "profile", "profile (where pairs across the subgenomes become common)"))


def _estimates_text(r):
    """Long-legend sentence naming the split estimates drawn for a subgenome pair (those with an interval)."""
    have = [(k, lab, d) for k, lab, d in _ESTIMATORS if r["split"].get(k) and np.isfinite(r["split"][k]["lo"])]
    if not have:
        return "No split estimate."
    ds = [d for _, _, d in have]
    joined = ds[0] if len(ds) == 1 else ", ".join(ds[:-1]) + " and " + ds[-1]
    rep = dict((k, lab) for k, lab, _ in have).get(r.get("split_reported"))
    return (f"Above, {['one split estimate', 'two split estimates', 'three split estimates'][len(have) - 1]} with "
            f"95% intervals: {joined}; filled: the one reported" + (f" ({rep}); estimators "
                                                                     f"{r.get('split_agreement')}." if rep else "."))


def _step_text(mg):
    """The merger clock's fitted step in words: phi applies to sibling pairs whose younger copy is younger than the
    merger, eps to older pairs."""
    return (f"the fitted step: {100 * mg['phi']:.0f}% across after the merger, {100 * mg['eps']:.0f}% before it "
            f"(misplaced copies, exchanges)")


def _ma(t):
    """Readable age from years."""
    return f"{t / 1e6:.2f} Myr" if t >= 1e6 else f"{t / 1e3:,.0f} kyr" if t >= 1e3 else f"{t:,.0f} years"


def _shown(lib, ph):
    """Phased sequences the figures draw, grouped by subgenome in natural order: all of them, or the FIG_MAX_SEQ
    longest -> (indices into ph['names'], lengths in bp of every phased sequence)."""
    si = {s: i for i, s in enumerate(lib.seqnames)}
    L = np.array([lib.seq_len[si[s]] for s in ph["names"]], float)
    keep = np.argsort(-L, kind="stable")[:FIG_MAX_SEQ]
    z = np.asarray(ph["z"])
    return sorted(keep.tolist(), key=lambda i: (z[i], natural_key(ph["names"][i]))), L


def _bin_shares(st, en, val, L, nb):
    """Length-weighted mean of tile values (tiles st..en, sorted, disjoint; val n x N) in nb equal bins of 0..L; bins
    partly or wholly outside the tiles keep the uncovered part at 0 -> nb x N."""
    w = L / nb
    edges = np.arange(nb + 1) * w
    xs = np.column_stack([st, en]).ravel().astype(float)
    out = np.zeros((nb, val.shape[1]))
    for g in range(val.shape[1]):
        Fe = np.cumsum((en - st) * val[:, g])
        F = np.column_stack([Fe - (en - st) * val[:, g], Fe]).ravel()
        out[:, g] = np.diff(np.interp(edges, xs, F, left=0.0, right=F[-1] if F.size else 0.0)) / w
    return np.clip(out, 0, 1)


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


def _structure_head(name, N, within, between, source):
    """Headline of fig_structure: inferred subgenomes, a number set by the user (-n) or a partition given by the user
    (--assign; never called subgenomes, and said to be unsupported when sharing is not higher within groups)."""
    rate = (f"chromosomes share LTR retrotransposon lineages at {within:.1f}× the expected rate within a group and "
            f"{between:.1f}× between groups.")
    if source == "user":
        return (f"{name} was split into {N} groups given by the user: {rate}"
                + ("" if within > between else " The transposons do not support this partition: groups do not "
                   "share more lineages within than between them."))
    if source == "given":
        return f"{name} was split into {N} groups (their number set by the user): {rate}"
    return (f"{name} carries {N} subgenomes, each with LTR retrotransposon lineages of its own: chromosomes share "
            f"lineages at {within:.1f}× the expected rate within a subgenome and {between:.1f}× between subgenomes. "
            f"Each subgenome kept the transposons its progenitor amplified before the progenitors met.")


def fig_structure(plt, a, lib, ph, dt, pa, anc, name):
    """Figure 1: A every phased sequence to scale by subgenome, every informative LTR-RT a tick (left: its lineages
    point to the sequence's own subgenome; right: to another, in that subgenome's color; opacity = confidence above
    the uninformative 1/N), called exchanges as boxes on the bars; B lineage sharing between sequences, observed over
    expected; C the dated subgenome tree (without dates: the replicability of every tested split); D per sequence, the
    share of informative copies most likely from another subgenome."""
    from matplotlib.colors import LogNorm, to_rgba
    from matplotlib.patches import Rectangle
    from matplotlib.ticker import NullFormatter, NullLocator
    N, names, z, ek = ph["N"], list(ph["names"]), np.asarray(ph["z"]), np.asarray(ph["ek"])
    shown, L = _shown(lib, ph)
    W, H = FIG_W, 8.6
    XL, XR, CW = 0.62, FIG_W - 0.12, 2.8
    fig = plt.figure(figsize=(W, H))
    inf = ek >= 0
    foreign, home, gbest = np.zeros(lib.n, bool), np.zeros(lib.n, bool), np.zeros(lib.n, int)
    if N > 1:
        inf, foreign, home, gbest = _favors(lib, ph, anc)
    else:
        home = inf.copy()
    axA = _ax(fig, XL, 3.3, XR - XL, 4.0)
    xpos, x, prev = {}, 0.0, None
    for i in shown:
        if prev is not None and z[i] != z[prev]:
            x += 1.0
        xpos[i], x, prev = x, x + 1.0, i
    Lmax = L[shown].max() / 1e6
    bw = 0.2
    axA.set_xlim(-0.7, x - 0.3)
    axA.set_ylim(Lmax * 1.02, 0)
    axA.set_ylabel("position (Mb)")
    axA.spines["bottom"].set_visible(False)
    axA.xaxis.tick_top()
    axA.set_xticks([xpos[i] for i in shown])
    axA.set_xticklabels([_short_name(names[i]) for i in shown], rotation=90, fontsize=6 if len(shown) <= 30 else 5)
    axA.tick_params(axis="x", length=0, pad=2)
    r = fig.canvas.get_renderer()
    names_h = max(t.get_window_extent(r).height for t in axA.get_xticklabels()) / fig.dpi
    top = H - 0.42 - 0.24 - names_h                                           # key strip, subgenome row, names
    axA.set_position([XL / W, 3.3 / H, (XR - XL) / W, (top - 3.3) / H])
    px_mb = (top - 3.3) * 300 / (Lmax * 1.02)                                # PNG pixel rows per Mb
    bw = 0.2
    calls = [c for c in (pa["segs"] if pa else []) if c["called"]]
    for i in shown:
        x0 = xpos[i]
        axA.add_patch(Rectangle((x0 - bw / 2, 0), bw, L[i] / 1e6, facecolor="white", edgecolor=FIG_DARK, lw=0.5,
                                zorder=3))
        m = inf & (ek == i)
        y = (lib.start[m] + lib.end[m]) / 2e6
        f, hm, gb = foreign[m], home[m], gbest[m]
        op = float(np.clip(L[i] / 1e6 * px_mb / max(m.sum(), 1), 0.03, 1.0))   # more copies than pixel rows: ink
        if (~f).any():                                                           # follows the share, not the overlap
            axA.hlines(y[~f], x0 - bw / 2 - 0.26, x0 - bw / 2 - 0.03, lw=0.35, zorder=2, rasterized=True,
                       colors=[to_rgba(FIG_GRAY, op * (0.75 if h_ else 0.3)) for h_ in hm[~f]])
        for g in range(N):
            mk = f & (gb == g)
            if mk.any():
                axA.hlines(y[mk], x0 + bw / 2 + 0.03, x0 + bw / 2 + 0.26, color=to_rgba(SG_COLORS[g], op), lw=0.35,
                           zorder=2, rasterized=True)
        for c in calls:
            if c["seq"] == names[i]:
                axA.add_patch(Rectangle((x0 - bw / 2, c["start"] / 1e6), bw, (c["end"] - c["start"]) / 1e6,
                                        facecolor=SG_COLORS[c["origin"]], edgecolor="none", zorder=4))
    xl0, xl1 = axA.get_xlim()
    xf = lambda xd: (XL + (xd - xl0) / (xl1 - xl0) * (XR - XL)) / W
    yb = (top + names_h + 0.08) / H
    for g in range(N):
        xs = [xpos[i] for i in shown if z[i] == g]
        if xs:
            col = SG_COLORS[g] if N > 1 else FIG_DARK
            fig.add_artist(plt.Line2D([xf(min(xs) - 0.35), xf(max(xs) + 0.35)], [yb, yb], color=col, lw=2,
                                      solid_capstyle="butt"))
            lab = f"{group_label(g)}: {int((z == g).sum())} sequences, {L[z == g].sum() / 1e6:,.0f} Mb"
            t = fig.text(xf((min(xs) + max(xs)) / 2), yb + 0.04 / H, lab, ha="center", va="bottom", color=col,
                         fontsize=7, fontweight="bold")
            if t.get_window_extent(r).width / fig.dpi > (xf(max(xs) + 0.5) - xf(min(xs) - 0.5)) * W:
                t.set_text(group_label(g))                             # narrow group: the name only (tables: rest)
    _letter(fig, axA, "A")
    if N > 1:
        _key(fig, XL, H - 0.2, XR - XL, [("tick", FIG_GRAY, "LTR-RT (left of its chromosome)"),
                                         ("tick", SG_COLORS[1], f"lineages favor another subgenome {BF_SUPPORT:g}-fold "
                                                                "(right, in its color)"),
                                         ("box", SG_COLORS[1], "called exchange (origin's color)")])
    else:
        _key(fig, XL, H - 0.2, XR - XL, [("tick", FIG_GRAY, "LTR-RT")])
    # B lineage sharing, observed over expected (the null of links_oe.tsv), in A's order
    side = 2.35
    axB = _ax(fig, XL, 0.55, side, side)
    Wr = np.asarray(ph["W_raw"], float)
    E = expected(Wr)
    with np.errstate(divide="ignore", invalid="ignore"):
        OE = np.where(E > 0, Wr / E, np.nan)[np.ix_(shown, shown)]
    np.fill_diagonal(OE, np.nan)
    im = axB.imshow(np.clip(OE, 0.25, 4), cmap="cividis", norm=LogNorm(0.25, 4), interpolation="nearest")
    n = len(shown)
    strip = max(0.6, n / 40)
    for k, i in enumerate(shown):
        col = SG_COLORS[z[i]] if N > 1 else FIG_GRAY
        axB.add_patch(Rectangle((-0.5 - 1.3 * strip, k - 0.5), strip, 1.0, color=col, clip_on=False, lw=0))
        axB.add_patch(Rectangle((k - 0.5, -0.5 - 1.3 * strip), 1.0, strip, color=col, clip_on=False, lw=0))
    axB.set_xlim(-0.5 - 1.4 * strip, n - 0.5)
    axB.set_ylim(n - 0.5, -0.5 - 1.4 * strip)
    axB.set_xticks([])
    axB.set_yticks([])
    for sp in axB.spines.values():
        sp.set_visible(False)
    cb = fig.colorbar(im, cax=_ax(fig, XL + side + 0.08, 0.95, 0.07, 1.5), ticks=[0.25, 0.5, 1, 2, 4])
    cb.ax.set_yticklabels(["¼×", "½×", "1×", "2×", "4×"])
    cb.ax.yaxis.set_minor_locator(NullLocator())
    cb.ax.yaxis.set_minor_formatter(NullFormatter())
    cb.outline.set_linewidth(0.4)
    cb.set_label("lineages shared,\nobserved / expected", fontsize=6.5)
    _letter(fig, axB, "B")
    # C dated tree, or the tested splits when nothing is dated; D shares from another subgenome
    xr0 = XR - CW
    axC = _ax(fig, xr0, 2.05, CW, 0.85) if N > 1 else _ax(fig, xr0, 0.55, CW, side)
    tree_txt = _draw_tree(axC, ph, dt, N) if N > 1 and dt else None
    if tree_txt is None:
        axC.cla()
        tree_txt = _draw_splits(axC, ph)
    _letter(fig, axC, "C")
    med_sh = {}
    if N > 1:
        axD = _ax(fig, xr0, 0.55, CW, 1.0)
        top_sh = 0.0
        for g in range(N):
            idx = [i for i in shown if z[i] == g]
            nf = np.array([(foreign & (ek == i)).sum() for i in idx], float)
            ni = np.array([(inf & (ek == i)).sum() for i in idx], float)
            with np.errstate(invalid="ignore", divide="ignore"):
                sh = np.where(ni >= 10, 100 * nf / ni, np.nan)
            if np.isfinite(sh).any():
                med_sh[g] = float(np.nanmedian(sh))
                top_sh = max(top_sh, float(np.nanmax(sh)))
            axD.scatter(g + np.linspace(-0.25, 0.25, len(idx)) if len(idx) > 1 else [g] * len(idx), sh, s=7,
                        color=SG_COLORS[g], lw=0, zorder=3)
        axD.set_ylim(0, max(10.0, top_sh * 1.15))
        axD.set_xlim(-0.6, N - 0.4)
        axD.set_xticks(range(N))
        axD.set_xticklabels([group_label(g) for g in range(N)], rotation=0 if N <= 6 else 90)
        axD.set_ylabel("foreign copies (%)")
        _letter(fig, axD, "D")
    # legend text from the run's numbers
    zz = z[shown]
    pw = [OE[p, q] for p in range(n) for q in range(n) if p != q and zz[p] == zz[q]]
    pb = [OE[p, q] for p in range(n) for q in range(n) if zz[p] != zz[q]]
    within = float(np.nanmedian(pw)) if np.isfinite(pw).any() else float("nan")
    between = float(np.nanmedian(pb)) if pb and np.isfinite(pb).any() else float("nan")
    n_inf, n_for = int((inf & (ek >= 0)).sum()) if N > 1 else 0, int(foreign.sum())
    more = (f" (the {n} longest of {len(names)} phased sequences; all are in chromosomes.tsv)" if n < len(names)
            else f" ({n} phased sequences)")
    unph = len(lib.seqnames) - len(names)
    given = ph.get("source") == "given"
    if N > 1:
        head = _structure_head(name, N, within, between, ph.get("source"))
        short = [head,
                 "**A** Every chromosome with its transposons placed by the subgenome their lineages come from; copies "
                 "of another subgenome's lineages, and called exchanges, mark where the subgenomes mixed.",
                 ("**B** Chromosomes of one subgenome share transposon lineages with each other, not with the rest."
                  if within > between else "**B** Lineage sharing between chromosomes, relative to chance."),
                 f"**C** {tree_txt[0]}",
                 "**D** On each chromosome, the share of transposons whose lineages favor another subgenome: low "
                 "where a chromosome kept its progenitor's transposons, high where it holds exchanges or is poorly "
                 "phased."]
    else:
        verdict = ph.get("verdict", "")
        head = (f"{name} shows no subgenome structure: no grouping of its chromosomes shares LTR retrotransposon "
                f"lineages reproducibly more within groups than between them. One progenitor genome, progenitors "
                f"too alike, or too few transposons to tell." if verdict.startswith("no structure") else
                f"{name}, one group of chromosomes: {verdict}.")
        short = [head,
                 "**A** Every chromosome with its transposons.",
                 "**B** Lineage sharing between every pair of chromosomes, relative to chance.",
                 f"**C** {tree_txt[0]}"]
    long_ = [(f"Terms. **LTR-RT**: an LTR retrotransposon copy. **Lineage**: copies sharing a 21-bp sequence found in "
              f"2-{CAP} copies on more than one sequence (a recent common source). **Informative copy**: one whose "
              f"lineages tell the subgenomes apart ({n_inf:,} of {lib.n:,} copies). Subgenomes are numbered SG1-SG{N} "
              f"by size.") if N > 1 else
             (f"Terms. **LTR-RT**: an LTR retrotransposon copy ({lib.n:,} here). **Lineage**: copies sharing a 21-bp "
              f"sequence found in 2-{CAP} copies on more than one sequence (a recent common source)."),
             f"**A** Sequences to scale{more}, grouped by subgenome"
             + (f"; {unph} sequences too short or too poor in LTR-RTs to phase are not drawn. " if unph else ". ")
             + (f"Each informative LTR-RT is a tick at its position: right of its sequence, in that subgenome's "
                f"color, when its lineages favor another subgenome at least {BF_SUPPORT:g}-fold over the sequence's own "
                f"({n_for:,} copies); left in gray otherwise (dark: favoring the sequence's own subgenome "
                f"{BF_SUPPORT:g}-fold; light: weaker evidence either way). Where a sequence holds more copies than the "
                f"image has pixel rows, every tick is drawn fainter in proportion, so color density shows the share "
                f"of copies, not their overlap. Boxes on the bars: the {len(calls)} called "
                f"exchanges (fig_exchanges), colored by origin." if N > 1 else "Each LTR-RT is a tick at its position."),
             f"**B** Lineage sharing between every pair of sequences, observed over the number expected from their "
             f"copy numbers (1× = random). Median within a subgenome {within:.2f}×"
             + (f", between subgenomes {between:.2f}×" if N > 1 else "") + ". Strips: each sequence's subgenome, in "
             f"A's order.",
             f"**C** {tree_txt[1]}"]
    if N > 1:
        long_.append(f"**D** Per sequence (dots), the share of its informative copies that are foreign: their "
                     f"lineages favor another subgenome at least {BF_SUPPORT:g}-fold over the sequence's own ("
                     + "; ".join(f"{group_label(g)} median {v:.1f}%" for g, v in med_sh.items())
                     + "); sequences with fewer than 10 informative copies are not drawn. Exchanges, insertions after "
                     "the merger, misplaced copies and, with similar subgenomes, chance all add to it.")
    if any(_short_name(names[i]) != names[i] for i in shown):
        long_.append("Sequence names longer than 30 characters are shortened in the middle (…); tables keep them.")
    long_.append(f"Phasing verdict: {ph['verdict']}" + ("; the number of subgenomes was set by the user." if given
                                                        else "; the partition was given by the user (--assign)."
                                                        if ph.get("source") == "user" else "."))
    return fig, short, long_


def _draw_tree(ax, ph, dt, N):
    """Dated subgenome tree on a time axis (million years ago, present at the right): a node for each accepted split
    of the phasing whose sides separate subgenomes (two subgenomes: their one split), at the median of the reported
    pairwise split estimates across it, its bar spanning the medians of their 95% bounds; subgenomes no accepted split
    separates (N set by the user, groups merged after the search) join one unresolved node. Dated mergers as dashed
    connectors between tips. -> (short, long) legend sentences, or None when nothing is dated."""
    from matplotlib.patches import Rectangle
    mu = dt["mu"]
    yr = lambda k: k / (2 * mu) / 1e6
    names, z = list(ph["names"]), np.asarray(ph["z"])
    zof = {nm: int(z[i]) for i, nm in enumerate(names)}
    splits = [(frozenset(zof[n] for n in t["part_a"] if n in zof), frozenset(zof[n] for n in t["part_b"] if n in zof))
              for t in ph.get("tree", []) if t.get("accepted")]
    rep = {}
    for r in dt["pairs"]:
        sp = r["split"].get(r.get("split_reported")) if r.get("split_reported") else None
        if sp:
            rep[frozenset(r["groups"])] = sp

    def node(groups):
        if len(groups) == 1:
            return next(iter(groups))
        sub = [(p, q) for p, q in splits if p | q == groups and p and q and not p & q]   # sides must separate
        if sub or len(groups) == 2:
            p, q = sub[0] if sub else (frozenset(sorted(groups)[:1]), frozenset(sorted(groups)[1:]))
            sides, resolved = [p, q], True
        else:
            sides, resolved = [frozenset([g]) for g in sorted(groups)], False
        est = [rep[frozenset((u, v))] for i_, a_ in enumerate(sides) for b_ in sides[i_ + 1:] for u in a_ for v in b_
               if frozenset((u, v)) in rep]
        f = lambda key: float(np.median([yr(e[key]) for e in est])) if est else float("nan")
        return dict(kids=[node(x) for x in sides], sides=sides, resolved=resolved, age=f("estimate"), lo=f("lo"),
                    hi=f("hi"))
    tree, tips = node(frozenset(range(N))), []

    def collect(t):
        if isinstance(t, int):
            tips.append(t)
        else:
            for k in t["kids"]:
                collect(k)
    collect(tree)
    ys = {g: k for k, g in enumerate(tips)}
    ages, txt, unres = [], [], []
    lab = lambda side: "+".join(group_label(g) for g in sorted(side))

    def draw(t):
        if isinstance(t, int):
            return ys[t], 0.0
        pts = [draw(k) for k in t["kids"]]
        yv = [p_[0] for p_ in pts]
        age, lo, hi = t["age"], t["lo"], t["hi"]
        if not np.isfinite(age):
            return float(np.mean(yv)), max(p_[1] for p_ in pts)
        ages.append(hi)
        ax.add_patch(Rectangle((lo, np.mean(yv) - 0.09), hi - lo, 0.18, facecolor=FIG_LIGHT, edgecolor="none",
                               zorder=1))
        for (yy, xx), sub in zip(pts, t["kids"]):
            ax.plot([xx, age], [yy, yy], color=SG_COLORS[sub] if isinstance(sub, int) else FIG_DARK,
                    lw=2.0 if isinstance(sub, int) else 1.0, solid_capstyle="butt", zorder=2)
        ax.plot([age, age], [min(yv), max(yv)], color=FIG_DARK, lw=1.0, zorder=2,
                ls="solid" if t["resolved"] else (0, (2, 1.5)))
        what = (" from ".join(lab(x) for x in t["sides"]) if t["resolved"] else
                ", ".join(lab(x) for x in t["sides"][:-1]) + " and " + lab(t["sides"][-1]) + " (unresolved)")
        txt.append(f"{what} {_ma(age * 1e6)} ({_ma(lo * 1e6)}-{_ma(hi * 1e6)})")
        if not t["resolved"]:
            unres.append(what)
        return float(np.mean(yv)), age
    draw(tree)
    if not ages:
        return None
    nm, mtxt = 0, []
    for r in dt["pairs"]:
        mg = r["merger"]
        if mg and mg["datable"] and not mg["bound"]:
            u, v = r["groups"]
            t = yr(mg["estimate"])
            for g in (u, v):
                ax.plot([yr(mg["lo"]), yr(mg["hi"])], [ys[g], ys[g]], color=FIG_DARK, lw=5, alpha=0.5,
                        solid_capstyle="butt", zorder=1)
            ax.plot([t, t], [ys[u], ys[v]], color=FIG_ACCENT, lw=0.8, ls=(0, (3, 1.5)), zorder=3)
            ax.text(t, max(ys[u], ys[v]) + 0.2, "merger", ha="center", va="bottom", fontsize=6)
            nm += 1
            mtxt.append(f"{group_label(u)}-{group_label(v)} {_ma(t * 1e6)} ({_ma(yr(mg['lo']) * 1e6)}-"
                        f"{_ma(yr(mg['hi']) * 1e6)})")
    xmax = max(ages) * 1.08
    for g in range(N):
        ax.text(-0.02 * xmax, ys[g], group_label(g), color=SG_COLORS[g], ha="left", va="center", fontweight="bold")
    ax.set_xlim(xmax, -0.16 * xmax)
    ax.spines["bottom"].set_bounds(0, xmax)
    ax.set_xticks([v for v in ax.get_xticks() if 0 <= v <= xmax])
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_ylim(-0.5, N - 0.3)
    ax.set_xlabel(f"million years ago ({_rate(mu)})", fontsize=6.5)
    s = ("The subgenomes' progenitors split: " + "; ".join(txt) + " ago"
         + (f", and they merged {'; '.join(mtxt)} ago." if nm else
            "; when they merged cannot be dated from the transposons."))
    lg = (s + " Nodes: " + ("the split of the two subgenomes" if N == 2 else "accepted splits of the phasing")
          + ", at the median of the reported split estimates between their subgenomes (fig_dating); light bars: the "
          + ("95% interval. " if N == 2 else "medians of their 95% interval bounds. ")
          + (f"Unresolved (dashed node): {'; '.join(unres)}: no accepted split separates these subgenomes (their "
             f"number set by the user, or groups merged after the search); the node sits at the median of all their "
             f"pairwise splits. " if unres else "")
          + ("Dashed connector and dark bands: dated mergers with 95% intervals. " if nm else "")
          + f"Ages assume {_rate(mu)}" + ("" if dt["mu_given"] else " (the default rate, not given for this species)")
          + "; another rate rescales every age.")
    return s, lg


def _draw_splits(ax, ph):
    """Every tested split of the structure search: replicability across disjoint half-libraries (points; filled =
    accepted) against the acceptance threshold -> (short, long) legend sentences."""
    T = ph.get("tree", [])
    if not T:
        ax.axis("off")
        msg = "No split could be tested: too few phased sequences or lineages."
        return msg, msg
    for k, t in enumerate(T):
        acc = t.get("accepted")
        ax.plot([k], [t["replicability"]], "o", ms=4, mfc=FIG_DARK if acc else "white", mec=FIG_DARK, mew=0.8)
    ax.axhline(MIN_REP, color=FIG_GRAY, lw=0.6, ls=(0, (3, 2)))
    ax.text(max(len(T), 1) - 0.45, MIN_REP + 0.03, "threshold", ha="right", va="bottom", fontsize=6, color=FIG_GRAY)
    ax.set_xticks(range(len(T)))
    ax.set_xticklabels([f"{len(t['part_a'])} | {len(t['part_b'])}" for t in T])
    ax.set_xlim(-0.6, max(len(T), 1) - 0.4)
    ax.set_ylim(min(-0.1, min([t["replicability"] for t in T], default=0) - 0.05), 1.05)
    ax.set_xlabel("tested split (sequences on each side)")
    ax.set_ylabel("replicability between\ndisjoint library halves")
    acc = sum(bool(t.get("accepted")) for t in T)
    return ((f"Of {len(T)} tested splits, {acc} {'is' if acc == 1 else 'are'} reproduced by independent halves of "
             f"the transposon library; no dates are available." if acc else
             f"None of the {len(T)} tested splits is reproduced by independent halves of the transposon library."),
            f"Each tested split of the sequences (labels: sequences on each side), with its replicability: agreement "
            f"(adjusted Rand index) between the splits found in disjoint halves of the lineages; filled: accepted "
            f"(>= {MIN_REP:g}, dashed). {acc} of {len(T)} accepted.")


def fig_dating(plt, a, lib, ph, dt, bio, r, name):
    """Figure 2 for one subgenome pair, drawing the inference: A the two clocks as a cartoon; B every sibling copy
    pair of the merger clock at its younger copy's age, within or across subgenomes, the share across per age bin and
    the fitted step; C every cross-subgenome copy pair by divergence, inherited insertions as ticks, the three split
    estimates; D each subgenome's insertion history on one time axis; E one typical inherited insertion:
    differences along its LTRs within each copy and between the two copies."""
    from matplotlib.patches import Rectangle
    from scipy.stats import rankdata
    mu = dt["mu"]
    x, y = r["groups"]
    lx, ly = group_label(x), group_label(y)
    ev = r["ev"]
    P, d = dt["pt"]["P"][ev["rows"]], dt["pt"]["d"][ev["rows"]]
    ek, z = np.asarray(ph["ek"]), np.asarray(ph["z"])
    loc = np.where(ek >= 0, z[np.maximum(ek, 0)], -1)
    cross = loc[P[:, 0]] != loc[P[:, 1]]
    orth, clk = ev["orth"], ev["clock"]
    yrs = lambda k: np.asarray(k, float) / (2 * mu)
    W, H = FIG_W, 8.2
    XL, XR, CW = 0.62, FIG_W - 0.12, 2.8
    xr0 = XR - CW
    fig = plt.figure(figsize=(W, H))
    # A cartoon (flat, schematic)
    axA = _ax(fig, XL, 5.75, CW, 2.1)
    axA.axis("off")
    axA.set_xlim(0, 10)
    axA.set_ylim(0, 10)
    cx, cy = SG_COLORS[x], SG_COLORS[y]
    axA.plot([0.3, 3.2], [5, 5], color=FIG_GRAY, lw=2.5, solid_capstyle="butt")
    axA.plot([3.2, 3.2], [2.2, 7.8], color=FIG_GRAY, lw=2.5, solid_capstyle="butt")
    axA.plot([3.2, 9.0], [7.8, 7.8], color=cy, lw=2.5, solid_capstyle="butt")
    axA.plot([3.2, 9.0], [2.2, 2.2], color=cx, lw=2.5, solid_capstyle="butt")
    axA.plot([6.4, 6.4], [2.2, 7.8], color=FIG_DARK, lw=0.9, ls=(0, (3, 2)))
    axA.text(3.2, 8.4, "split", ha="center", va="bottom", fontsize=6.5)
    axA.text(6.4, 8.4, "merger", ha="center", va="bottom", fontsize=6.5)
    axA.plot([1.6], [5.35], marker="v", color=FIG_ACCENT, ms=5, lw=0)
    for yy in (7.8, 2.2):
        axA.plot([4.6], [yy + 0.35], marker="v", color=FIG_ACCENT, ms=4, lw=0)
    axA.text(1.6, 5.9, "insertion", ha="center", va="bottom", fontsize=6)
    axA.text(4.6, 5.0, "inherited\nby both", ha="center", va="center", fontsize=6)
    axA.plot([7.4], [2.55], marker="v", color=FIG_DARK, ms=5, lw=0)
    axA.annotate("", xy=(8.2, 7.4), xytext=(7.5, 2.9), arrowprops=dict(arrowstyle="-|>", lw=0.7, color=FIG_DARK))
    axA.text(8.45, 4.9, "new copy\nlands\nacross", ha="left", va="center", fontsize=6)
    axA.text(9.15, 7.8, ly, color=cy, fontweight="bold", va="center", fontsize=7)
    axA.text(9.15, 2.2, lx, color=cx, fontweight="bold", va="center", fontsize=7)
    axA.annotate("", xy=(9.0, 0.6), xytext=(0.3, 0.6), arrowprops=dict(arrowstyle="-|>", lw=0.6, color=FIG_GRAY))
    axA.text(4.65, 0.0, "time", ha="center", va="bottom", fontsize=6, color=FIG_GRAY)
    _letter(fig, axA, "A")
    # B merger clock: every sibling pair it used
    axB = _ax(fig, xr0, 5.75, CW, 2.1)
    young = yrs(np.fmin(lib.k2p[P[clk, 0]], lib.k2p[P[clk, 1]]))
    c = cross[clk]
    nsite = lib.n_sites[np.isfinite(lib.n_sites) & (lib.n_sites > 0)]
    nsite = nsite if nsite.size else lib.ltr5_len[np.isfinite(lib.ltr5_len) & (lib.ltr5_len > 0)]
    floor = float(yrs(0.5 / (np.median(nsite) if nsite.size else 500.0)))      # half a substitution per LTR
    pos = young > floor
    x0 = floor * 0.5                                                             # column of identical LTRs (age 0)
    ap = np.where(pos, young, x0)
    jit = np.random.default_rng(1).uniform(-0.06, 0.06, young.size)
    axB.scatter(ap[~c], jit[~c], s=1.5, color=FIG_GRAY, alpha=0.35, lw=0, rasterized=True)
    axB.scatter(ap[c], 1 + jit[c], s=1.5, color=FIG_ACCENT, alpha=0.6, lw=0, rasterized=True)
    hi_age = max(float(young.max()) if young.size else floor * 10, floor * 10)
    edges = np.geomspace(floor, hi_age * 1.01, 13)
    bi = np.clip(np.searchsorted(edges, young) - 1, 0, edges.size - 2)
    for k in range(edges.size - 1):
        mk = (bi == k) & pos
        if mk.sum() >= 10:
            axB.plot(np.sqrt(edges[k] * edges[k + 1]), c[mk].mean(), "o", mfc="white", mec="black", mew=0.8, ms=3.5,
                     zorder=5)
    if (~pos).sum() >= 10:
        axB.plot(x0, c[~pos].mean(), "o", mfc="white", mec="black", mew=0.8, ms=3.5, zorder=5)
    axB.axvline(floor * 0.75, color=FIG_LIGHT, lw=0.6)
    mg = r["merger"]
    dated = bool(mg and mg["datable"] and not mg["bound"])
    if dated:
        T = float(yrs(mg["estimate"]))
        axB.axvspan(float(yrs(mg["lo"])), float(yrs(mg["hi"])), color=FIG_LIGHT, zorder=0, lw=0)
        axB.plot([floor, T, T, edges[-1]], [mg["phi"], mg["phi"], mg["eps"], mg["eps"]], color="black", lw=0.8)
        axB.text(T, 1.16, "merger", ha="center", va="bottom", fontsize=6)
    elif mg and mg["datable"]:
        axB.axvline(float(yrs(mg["hi"])), color="black", lw=0.8, ls=(0, (3, 2)))
        axB.text(float(yrs(mg["hi"])), 1.16, "onset bound", ha="center", va="bottom", fontsize=6)
    axB.set_xscale("log")
    axB.set_xlim(edges[-1] * 1.6, floor * 0.3)
    _log_axis(axB, "x", "years")
    tk = [v for v in axB.get_xticks() if floor <= v <= edges[-1] * 1.6]
    axB.set_xticks([x0] + tk, ["0"] + [_yr_label(v) for v in tk])
    axB.text(edges[-1] * 1.5, 1.07, "across", ha="left", va="bottom", fontsize=6, color=FIG_DARK)
    axB.text(edges[-1] * 1.5, -0.07, "within", ha="left", va="top", fontsize=6, color=FIG_DARK)
    axB.set_ylim(-0.17, 1.17)
    axB.set_yticks([0, 0.5, 1])
    axB.set_yticklabels(["0", "50", "100"])
    axB.spines["left"].set_bounds(0, 1)
    axB.set_ylabel("sibling pairs split\nbetween subgenomes (%)")
    axB.set_xlabel(f"age of the younger copy (years ago)\nat {_rate(mu)}")
    _letter(fig, axB, "B")
    # C split: every cross-subgenome pair by divergence; estimates in a band above, in the same axes
    axC = _ax(fig, XL, 2.85, CW, 2.1)
    tc = yrs(d[cross])
    tpos = tc[tc > 0]
    n_zero = int((tc <= 0).sum())
    sp = r["split"].get(r.get("split_reported")) if r.get("split_reported") else None
    lab_est = {"orthologs": "inherited", "confinement": "confinement", "profile": "profile"}
    if tpos.size:
        lo_b = max(float(np.quantile(tpos, 0.002)), 1e3)
        bins = np.geomspace(lo_b, max(tpos.max(), lo_b * 10) * 1.01, 45)
        cnt, _ = np.histogram(tpos, bins)
        ymax = max(cnt.max(), 1)
        axC.stairs(cnt, bins, fill=True, color=FIG_LIGHT, lw=0)
        to = yrs(d[cross & orth])
        axC.vlines(to[to > 0], ymax * 1.05, ymax * 1.17, color=FIG_ACCENT, lw=0.4, alpha=0.8)
        axC.text(max(to.max() if to.size else bins[-1], bins[0]) * 1.2, ymax * 1.11,
                 f"inherited (n = {int((cross & orth).sum())})", ha="right", va="center", fontsize=6)
        for j, k in enumerate(("orthologs", "confinement", "profile")):
            yy = ymax * (1.31 + 0.13 * j)
            v = r["split"].get(k)
            if v and np.isfinite(v["lo"]):
                l_, h_ = float(yrs(v["lo"])), float(yrs(v["hi"]))
                axC.plot([l_, h_], [yy, yy], color="black", lw=0.9)
                axC.plot([float(yrs(v["estimate"]))], [yy], "o", ms=3.5, lw=0, mec="black", mew=0.7,
                         color=FIG_ACCENT if k == r.get("split_reported") else "white")
                axC.text(h_ * 1.15, yy, lab_est[k], va="center", ha="right", fontsize=6)
        axC.set_xscale("log")
        axC.set_xlim(bins[-1] * 6, bins[0] * 0.8)
        _log_axis(axC, "x", "years")
        axC.set_ylim(0, ymax * 1.62)
        axC.spines["left"].set_bounds(0, ymax)
        axC.set_yticks([t for t in axC.get_yticks() if 0 <= t <= ymax])
    axC.set_xlabel(f"common ancestor of the two copies (years ago)\nat {_rate(mu)}")
    axC.set_ylabel("copy pairs across\nthe subgenomes")
    _letter(fig, axC, "C")
    # D insertion histories per subgenome
    axD = _ax(fig, xr0, 2.85, CW, 2.1)
    peaks, used_obs, curves = {}, False, {}
    for g in (x, y):
        R_ = [h for h in bio["history"] if h["subgenome"] == g and h["clade"] == "all" and h["lo"] > 0]
        if not R_:
            continue
        lo, hi = np.array([h["lo"] for h in R_]), np.array([h["hi"] for h in R_])
        dec = np.array([h["deconvolved"] for h in R_], float)
        if not np.isfinite(dec).all():
            dec, used_obs = np.array([h["observed"] for h in R_], float), True
        yv = dec / np.log10(hi / lo)
        mid = yrs(np.sqrt(lo * hi))
        axD.plot(mid, yv, color=SG_COLORS[g], lw=1.2)
        peaks[g] = float(mid[int(np.nanargmax(yv))])
        curves[g] = (mid, yv)
    for g, (mid, yv) in curves.items():                       # direct label at the peak: above it, or below it
        j = int(np.nanargmax(yv))                               # when the other curve passes above there
        other = [c_[1] for h_, c_ in curves.items() if h_ != g]
        below = bool(other and other[0].size == yv.size and other[0][j] > yv[j])
        axD.annotate(group_label(g), (mid[j], yv[j]), xytext=(0, -4 if below else 3), textcoords="offset points",
                     color=SG_COLORS[g], ha="center", va="top" if below else "bottom", fontsize=6.5, fontweight="bold")
    yt = axD.get_ylim()[1]
    for est, lab in ((sp, "split"), (mg if dated else None, "merger")):
        if est:
            axD.axvspan(float(yrs(est["lo"])), float(yrs(est["hi"])), color=FIG_LIGHT, zorder=0, lw=0)
            axD.text(float(yrs(est["estimate"])), yt * 1.04, lab, ha="center", va="bottom", fontsize=6)
    axD.set_ylim(0, yt * 1.14)
    axD.spines["left"].set_bounds(0, yt)
    axD.set_yticks([t for t in axD.get_yticks() if 0 <= t <= yt])
    axD.set_xscale("log")
    xl = axD.get_xlim()
    axD.set_xlim(max(xl), min(xl))
    _log_axis(axD, "x", "years")
    axD.set_xlabel(f"insertion (years ago)\nat {_rate(mu)}")
    axD.set_ylabel("insertions per\ntenfold age range")
    _letter(fig, axD, "D")
    # E worked example: the most typical inherited insertion in what E shows (difference rate between the copies'
    # 5' LTRs and LTR length nearest the medians; at most FIG_EXAMPLES candidates, drawn at random)
    axE = _ax(fig, XL, 0.45, XR - XL, 1.6)
    inh = np.flatnonzero(cross & orth)
    lts = {}
    for k in inh:
        lt_ = _ltrs(lib, int(P[k, 0]), int(P[k, 1]))
        if lt_ is not None:
            lts[int(k)] = lt_
    n_lt = len(lts)
    cand = sorted(lts)
    if len(cand) > FIG_EXAMPLES:
        cand = sorted(np.random.default_rng(1).choice(cand, FIG_EXAMPLES, replace=False).tolist())
    ex_long = "No inherited insertion with both LTRs annotated in both copies: no worked example."
    if cand:
        rate = [_diff_positions(lts[k][0], lts[k][2]).size / max(lts[k][0].size, 1) for k in cand]
        rd = rankdata(rate) / len(cand)
        rl = rankdata([lts[k][0].size for k in cand]) / len(cand)
        k = cand[int(np.argmin(np.abs(rd - 0.5) + np.abs(rl - 0.5)))]
        ca, cb = int(P[k, 0]), int(P[k, 1])
        a5, a3, b5, b3 = lts[k]
        za, zb = int(loc[ca]), int(loc[cb])
        rows = [(f"copy 1 ({group_label(za)}): 5' vs 3' LTR", a5, a3, SG_COLORS[za]),
                (f"copy 2 ({group_label(zb)}): 5' vs 3' LTR", b5, b3, SG_COLORS[zb]),
                ("copy 1 vs copy 2: 5' LTRs", a5, b5, FIG_ACCENT)]
        Lm = max(rw[1].size for rw in rows)
        nd = []
        for j, (nm, q, tt, col) in enumerate(rows):
            pz = _diff_positions(q, tt)
            nd.append(pz.size)
            yy = 2 - j
            axE.add_patch(Rectangle((0, yy - 0.2), q.size, 0.4, facecolor="#f2f2f2", edgecolor=FIG_GRAY, lw=0.5))
            axE.vlines(pz, yy - 0.2, yy + 0.2, color=col, lw=0.8)
            axE.text(-Lm * 0.015, yy, nm, ha="right", va="center", fontsize=6.5)
            axE.text(q.size + Lm * 0.015, yy, f"{pz.size} differences", ha="left", va="center", fontsize=6.5)
        where = [f"{lib.seqnames[lib.seq[c_]]}:{lib.start[c_]:,}-{lib.end[c_]:,}" for c_ in (ca, cb)]
        axE.set_xlim(-Lm * 0.36, Lm * 1.14)
        axE.set_ylim(-0.55, 2.55)
        axE.spines["bottom"].set_bounds(0, Lm)
        axE.set_xticks([v for v in axE.get_xticks() if 0 <= v <= Lm])
        ex_long = (f"Copy 1 ({group_label(za)}, {where[0]}) and copy 2 ({group_label(zb)}, {where[1]}): their 5' LTRs "
                   f"differ at {nd[2]} sites, fewer than each copy's own two LTRs ({nd[0]} and {nd[1]}). A copy's two "
                   f"LTRs start identical at its insertion; two copies of one insertion start identical at the split. "
                   f"Fewer differences between the copies than within them place this insertion before the split, "
                   f"inherited by both subgenomes. "
                   + ("The only inherited insertion with both LTRs annotated. " if n_lt == 1 else
                      "One of the two inherited insertions with both LTRs annotated. " if n_lt == 2 else
                      f"Chosen as typical: the difference rate between its copies' 5' LTRs and its LTR length lie "
                      f"nearest the medians of the {n_lt:,} inherited insertions with both LTRs annotated"
                      + (f" ({len(cand)} drawn at random)" if len(cand) < n_lt else "") + ". ")
                   + f"Ticks: differences on a "
                   f"unit-cost alignment (display); the dating counts substitutions on an affine realignment of whole "
                   f"elements.")
    else:
        axE.axis("off")
    axE.set_yticks([])
    axE.spines["left"].set_visible(False)
    axE.set_xlabel("position along the LTR (bp)")
    _letter(fig, axE, "E")
    # legends
    ma = lambda k: _ma(float(yrs(k)))
    if mg is None:
        mtxt = f"the merger is not datable: too few sibling copy pairs ({r.get('n_sibling_pairs', 0)})"
    elif not mg["datable"]:
        mtxt = ("the merger is not datable: young sibling copies are not split between the subgenomes more often than "
                "old ones (no copies moved after the merger, or too few to tell)"
                + (f"; if new copies land in either subgenome at random, the merger is younger than "
                   f"{ma(mg['landing_bound'])} (an assumption that fails where transposition stays within each "
                   f"subgenome)" if mg.get("landing_bound") else ""))
    elif mg["bound"]:
        mtxt = f"copies began moving between the subgenomes less than {ma(mg['hi'])} ago (a bound, not a date)"
    else:
        mtxt = f"copies began moving between the subgenomes {ma(mg['estimate'])} ago ({ma(mg['lo'])}-{ma(mg['hi'])})"
    stxt = (f"the progenitors of {lx} and {ly} split {ma(sp['estimate'])} ago ({ma(sp['lo'])}-{ma(sp['hi'])})" if sp
            else f"the split of {lx} and {ly} cannot be dated")
    m_short = (f"they merged {ma(mg['estimate'])} ago ({ma(mg['lo'])}-{ma(mg['hi'])})" if dated else
               f"copies began moving between them less than {ma(mg['hi'])} ago" if mg and mg["datable"] else
               "the merger cannot be dated (too few young copies moving between them)"
               + (f"; if new copies land in either subgenome at random, it is younger than {ma(mg['landing_bound'])}"
                  if mg and mg.get("landing_bound") else ""))
    short = [f"In {name}, {stxt}; {m_short}. Insertions inherited by both subgenomes date the split; new copies "
             f"landing in the other subgenome date the merger.",
             "**A** The two clocks: an insertion made before the split sits in both subgenomes; a copy made after the "
             "merger can land in either.",
             ("**B** Young sibling copies sit in different subgenomes far more often than old ones; the age where that "
              "stops is the merger." if dated else
              "**B** Sibling copies of every age, within or across the subgenomes: young ones show no clear excess of "
              "moves, so no merger date."),
             "**C** Copy pairs across the two subgenomes by age; insertions inherited by both (ticks) mark the "
             "split.",
             "**D** When each subgenome's transposons were active.",
             "**E** One inherited insertion: its two copies differ less from each other than each copy's own LTRs do."]
    n_sib, n_cr = int(clk.sum()), int(c.sum())
    long_ = [f"Terms. Time runs from left (past) to right (present) in every panel. **Age** of a copy: the "
             f"divergence of its two LTRs (corrected for multiple hits) over twice the "
             f"rate, {_rate(mu)}" + ("" if dt["mu_given"] else " (default, not given for this species)") + ". "
             f"**Sibling pair**: two related copies no more diverged than the older copy's own LTRs (the older copy's "
             f"insertion is their ancestor). **Inherited insertion**: a pair across the subgenomes whose divergence is "
             f"significantly below both copies' own LTR divergence (one insertion, older than the split).",
             "**A** Cartoon, not data.",
             f"**B** {n_sib:,} sibling pairs within or across {lx} and {ly}, each a dot at its younger copy's age "
             f"(identical LTRs: 0), top row if the two copies sit in different subgenomes ({n_cr:,}), bottom row if "
             f"not; open circles: percentage across per age bin (bins with >= 10 pairs). "
             + (f"Line: {_step_text(mg)}; band: 95% interval. " if dated else "")
             + f"Result: {mtxt}.",
             f"**C** {int(cross.sum()):,} sampled copy pairs across {lx} and {ly} by divergence (in years: divergence "
             f"over twice the rate; {n_zero} identical pairs not drawn); ticks: the "
             f"{_n(int((cross & orth).sum()), 'inherited insertion')}. {_estimates_text(r)}",
             f"**D** Insertions per age in each subgenome, all LTR-RT clades together; "
             + ("observed LTR ages (no substitution counts to deconvolve)" if used_obs else
                "LTR ages deconvolved for the counting noise of short LTRs")
             + " (peaks: " + "; ".join(f"{group_label(g)} {_ma(v)}" for g, v in peaks.items()) + "). Old copies decay "
             "and are lost: older activity is under-represented. Bands: split and merger intervals.",
             f"**E** {ex_long}"]
    return fig, short, long_


def fig_exchanges(plt, a, lib, ph, pa, anc, name):
    """Figure 3: A every phased sequence painted by ancestry (posterior shares stacked per bin; copy-free stretches
    hatched as unknown), called exchanges boxed; B beside A on the same rows: the smallest exchange the transposons
    could reveal on each sequence and every called exchange's length; C zoom on the typical called exchange (median
    evidence) in one browser panel: posterior trace, the call with its breakpoint intervals, every LTR-RT in stacked
    rows colored by its most likely subgenome, the smallest detectable exchange along the window."""
    from matplotlib.colors import to_rgba
    from matplotlib.patches import ConnectionPatch, Rectangle
    N, names, z, ek = ph["N"], list(ph["names"]), np.asarray(ph["z"]), np.asarray(ph["ek"])
    shown, Larr = _shown(lib, ph)
    Ls = {names[i]: int(Larr[i]) for i in range(len(names))}
    tracks = {t[0]: t for t in pa["post"]}
    res = {r_[0]: r_ for r_ in pa["res"]}
    rr = np.concatenate([r_[3] for r_ in pa["res"]]) if pa["res"] else np.zeros(0)
    R = float(np.median(rr)) if rr.size else float(RES_WIN)
    calls = [c for c in pa["segs"] if c["called"]]
    nrow = len(shown)
    rowh = float(np.clip(0.12 * nrow, 1.8, 4.6))
    XR = FIG_W - 0.12
    fs = 6 if nrow <= 30 else 5
    probe = plt.figure(figsize=(FIG_W, 2))
    rp = probe.canvas.get_renderer()
    wname = max(probe.text(0, 0, _short_name(names[i]), fontsize=fs).get_window_extent(rp).width
                for i in shown) / probe.dpi
    plt.close(probe)
    XL = max(1.05, wname + 0.45)                              # names, subgenome bars and their labels fit
    if calls:
        yA = 3.6
    else:
        yA, rowh = 0.25, max(rowh, 2.75)
    items = [("box", SG_COLORS[g], f"{group_label(g)} ancestry") for g in range(N)] + \
        [("hatch", FIG_GRAY, "unknown (no informative LTR-RTs)")] + \
        ([("outline", "black", "called exchange")] if calls else []) + \
        [("dot", FIG_GRAY, "smallest detectable exchange (B)")] + \
        ([("diamond", FIG_DARK, "called exchange length (B; origin's color)")] if calls else [])
    probe = plt.figure(figsize=(FIG_W, 2))
    krows = _key(probe, XL, 1.0, XR - XL, items, draw=False)
    plt.close(probe)
    H = yA + rowh + 0.62 + 0.16 * krows
    fig = plt.figure(figsize=(FIG_W, H))
    WA = XR - XL - 0.6 - 1.37                                     # A and C share the left column; B, D 1.37 in
    axA = _ax(fig, XL, yA, WA, rowh)
    h = 0.72
    maxL = max(Ls[names[i]] for i in shown) / 1e6
    rowy = {names[i]: nrow - 1 - k for k, i in enumerate(shown)}
    for k, i in enumerate(shown):
        s, y0 = names[i], nrow - 1 - k
        axA.add_patch(Rectangle((0, y0 - h / 2), Ls[s] / 1e6, h, facecolor="white", edgecolor=FIG_GRAY,
                                hatch="//////", lw=0, zorder=0))
        t = tracks.get(s)
        if t is not None:
            st, en, val = posterior_tiles(t[1], t[2], Ls[s], R)
            nb = max(1, int(1500 * Ls[s] / (maxL * 1e6)))
            share = _bin_shares(st, en, val, Ls[s], nb)
            xs = np.arange(nb + 1) * (Ls[s] / nb) / 1e6
            base = np.zeros(nb)
            for g in [int(z[i])] + [g for g in range(N) if g != z[i]]:
                top = base + share[:, g]
                axA.fill_between(xs, y0 - h / 2 + np.append(base, base[-1]) * h,
                                 y0 - h / 2 + np.append(top, top[-1]) * h, step="post", color=SG_COLORS[g], lw=0, zorder=1)
                base = top
        axA.add_patch(Rectangle((0, y0 - h / 2), Ls[s] / 1e6, h, facecolor="none", edgecolor=FIG_DARK, lw=0.4,
                                zorder=3))
        for c in calls:
            if c["seq"] == s:
                axA.add_patch(Rectangle((c["start"] / 1e6, y0 - h / 2 - 0.12), (c["end"] - c["start"]) / 1e6, h + 0.24,
                                        facecolor="none", edgecolor="black", lw=0.8, zorder=4))
    axA.set_yticks(range(nrow))
    axA.set_yticklabels([_short_name(names[i]) for i in shown][::-1], fontsize=fs)
    axA.tick_params(axis="y", length=0, pad=2)
    axA.spines["left"].set_visible(False)
    axA.spines["bottom"].set_visible(False)
    axA.spines["top"].set_visible(True)
    axA.xaxis.tick_top()
    axA.xaxis.set_label_position("top")
    axA.set_xlim(0, maxL * 1.01)
    axA.set_ylim(-0.6, nrow - 0.4)
    axA.set_xlabel("position (Mb)")
    r = fig.canvas.get_renderer()
    wname = max(t.get_window_extent(r).width for t in axA.get_yticklabels()) / fig.dpi
    for g in range(N):
        rows_ = [nrow - 1 - k for k, i in enumerate(shown) if z[i] == g]
        if rows_:
            yy = lambda v: (yA + (v + 0.6) / nrow * rowh) / H
            xb = (XL - wname - 0.1) / FIG_W
            fig.add_artist(plt.Line2D([xb, xb], [yy(min(rows_) - 0.4), yy(max(rows_) + 0.4)], color=SG_COLORS[g],
                                      lw=2.2, solid_capstyle="butt"))
            fig.text(xb - 0.05 / FIG_W, yy((min(rows_) + max(rows_)) / 2), group_label(g), color=SG_COLORS[g],
                     rotation=90, ha="right", va="center", fontweight="bold")
    _letter(fig, axA, "A")
    _key(fig, XL, H - 0.14, XR - XL, items)
    # B smallest detectable exchange and called exchanges, rows shared with A
    xb0 = XL + WA + 0.6
    axB = _ax(fig, xb0, yA, XR - xb0, rowh)
    for k, i in enumerate(shown):
        s, y0 = names[i], nrow - 1 - k
        if s in res and res[s][3].size:
            axB.plot([np.median(res[s][3]) / 1e6], [y0], "o", ms=3, color=FIG_GRAY, lw=0, zorder=2)
        for c in calls:
            if c["seq"] == s:
                axB.plot([(c["end"] - c["start"] + 1) / 1e6], [y0], "D", ms=3.5, color=SG_COLORS[c["origin"]], lw=0,
                         zorder=3)
    axB.set_xscale("log")
    allx = [np.median(res[names[i]][3]) / 1e6 for i in shown if names[i] in res and res[names[i]][3].size] + \
        [(c["end"] - c["start"] + 1) / 1e6 for c in calls]
    if allx:
        axB.set_xlim(min(allx) / 2, max(allx) * 2)
    _log_axis(axB, "x", "Mb")
    axB.set_ylim(-0.6, nrow - 0.4)
    axB.set_yticks([])
    axB.spines["left"].set_visible(False)
    axB.spines["bottom"].set_visible(False)
    axB.spines["top"].set_visible(True)
    axB.xaxis.tick_top()
    axB.xaxis.set_label_position("top")
    for k in range(nrow):
        axB.axhline(k, color="#efefef", lw=0.5, zorder=0)
    axB.set_xlabel("exchange length (Mb)")
    _letter(fig, axB, "B")
    ex_txt = ""
    if calls:
        pool = [c_ for c_ in calls if c_["seq"] in rowy] or calls        # the showcase must sit on a drawn sequence
        lrs = np.array([c_["llr"] for c_ in pool])
        c = pool[int(np.argsort(lrs, kind="stable")[len(lrs) // 2])]
        span, s = c["end"] - c["start"] + 1, c["seq"]
        w0, w1 = max(0, c["left_flank"] - span), min(Ls[s], c["right_flank"] + span)
        i_s = names.index(s)
        m = np.flatnonzero((ek == i_s) & (lib.end >= w0) & (lib.start <= w1))
        lP, ninf = anc
        Pm = np.exp(lP[m])
        bst = Pm.argmax(1)
        inf_m = ninf[m] > 0
        al = np.clip((Pm.max(1) - 1 / N) / (1 - 1 / N), 0.15, 1)
        axC = _ax(fig, XL, 0.45, WA, 2.6)
        st, en = lib.start[m], lib.end[m]
        row = _rows(st, en, (w1 - w0) * 0.002)
        nr = int(row.max()) + 1 if row.size else 1
        yR0, yE0 = 0.0, 2.2                     # tracks bottom-up: resolution, LTR-RTs, call, posterior
        yC = yE0 + nr * 0.8 + 0.6
        yP0 = yC + 1.0
        for k in range(m.size):
            col = to_rgba(SG_COLORS[bst[k]], al[k]) if inf_m[k] else to_rgba(FIG_GRAY, 0.5)
            axC.add_patch(Rectangle((st[k] / 1e6, yE0 + (nr - 1 - row[k]) * 0.8), (en[k] - st[k] + 1) / 1e6, 0.55,
                                    facecolor=col, edgecolor="none"))
        axC.add_patch(Rectangle((c["start"] / 1e6, yC - 0.25), span / 1e6, 0.5, facecolor=SG_COLORS[c["origin"]],
                                edgecolor="none"))
        for p_, q_ in ((c["left_flank"], c["start"]), (c["end"], c["right_flank"])):
            axC.add_patch(Rectangle((p_ / 1e6, yC - 0.25), (q_ - p_) / 1e6, 0.5, facecolor="none", hatch="//////",
                                    edgecolor=SG_COLORS[c["origin"]], lw=0))
        t = tracks[s]
        mm = (t[1] >= w0) & (t[1] <= w1)
        o = np.argsort(t[1][mm])
        for g in range(N):
            axC.plot(t[1][mm][o] / 1e6, yP0 + 1.4 * t[2][mm][o, g], color=SG_COLORS[g], lw=0.9)
        rmax = float(span)
        if s in res:
            rs, re_, rb = res[s][1], res[s][2], res[s][3]
            mk = (re_ >= w0) & (rs <= w1)
            if mk.any():
                rmax = max(rmax, float(rb[mk].max()))
                xs_ = np.column_stack([rs[mk], re_[mk]]).ravel() / 1e6
                axC.plot(xs_, np.repeat(rb[mk], 2) / rmax * 1.4 + yR0, color=FIG_DARK, lw=0.8)
        axC.plot([w0 / 1e6, w1 / 1e6], [yR0 + span / rmax * 1.4] * 2, color=SG_COLORS[c["origin"]], lw=0.6,
                 ls=(0, (3, 2)))
        for yy, lab in ((yP0 + 0.7, "posterior\n(0-1)"), (yC, "call"), (yE0 + nr * 0.4, "LTR-RTs"),
                        (yR0 + 0.7, f"smallest\ndetectable\n(0-{rmax / 1e6:.2g} Mb)")):
            axC.text(-0.015, yy, lab, transform=axC.get_yaxis_transform(), ha="right", va="center", fontsize=6.5)
        axC.set_xlim(w0 / 1e6, w1 / 1e6)
        axC.set_ylim(-0.2, yP0 + 1.55)
        axC.set_yticks([])
        axC.spines["left"].set_visible(False)
        axC.set_xlabel(f"{s} position (Mb)")
        _letter(fig, axC, "C")
        if s in rowy:
            axA.add_patch(Rectangle((w0 / 1e6, rowy[s] - 0.5), (w1 - w0) / 1e6, 1.0, facecolor="none",
                                    edgecolor="black", lw=0.6, ls=(0, (2, 1)), zorder=5))
        axA.plot([w0 / 1e6, w1 / 1e6], [-0.6, -0.6], color="black", lw=2, solid_capstyle="butt", clip_on=False,
                 zorder=6)                                     # the window on A's position axis, joined to C
        for xa in (w0 / 1e6, w1 / 1e6):
            fig.add_artist(ConnectionPatch((xa, -0.6), (xa, yP0 + 1.55), "data", "data", axesA=axA, axesB=axC,
                                           color=FIG_GRAY, lw=0.5, ls=(0, (2, 1))))
        axD = _ax(fig, xb0, 0.45, XR - xb0, 2.6)
        segs = pa["segs"]
        Lseg = np.array([sg_["end"] - sg_["start"] + 1 for sg_ in segs]) / 1e6
        lr = np.array([sg_["llr"] for sg_ in segs]) / math.log(10)
        cd = np.array([sg_["called"] for sg_ in segs])
        axD.scatter(Lseg[~cd], lr[~cd], s=9, facecolors="none", edgecolors=FIG_GRAY, lw=0.6, zorder=2)
        for g in range(N):
            mk = cd & np.array([sg_["origin"] == g for sg_ in segs])
            axD.scatter(Lseg[mk], lr[mk], s=12, marker="D", color=SG_COLORS[g], lw=0, zorder=3)
        jc = segs.index(c)
        axD.scatter([Lseg[jc]], [lr[jc]], s=60, facecolors="none", edgecolors="black", lw=0.8, zorder=4)
        axD.axhline(math.log10(BF_DECISIVE), color=FIG_DARK, lw=0.6, ls=(0, (3, 2)), zorder=1)
        axD.set_xscale("log")
        axD.set_xlim(Lseg.min() / 1.6, Lseg.max() * 1.6)
        _log_axis(axD, "x", "Mb")
        from matplotlib.ticker import FuncFormatter, MaxNLocator
        axD.set_ylim(min(-0.5, lr.min() - 0.5), lr.max() * 1.08 + 0.5)
        axD.yaxis.set_major_locator(MaxNLocator(5, integer=True))
        axD.yaxis.set_major_formatter(FuncFormatter(
            lambda v, _: {0: "1", 1: "10", 2: "100"}.get(int(v), f"$10^{{{int(v)}}}$")))
        axD.set_xlabel("segment length (Mb)")
        axD.set_ylabel("likelihood ratio")
        _letter(fig, axD, "D")
        lr_txt = (f"{math.exp(c['llr']):,.0f}" if c["llr"] < 7 * math.log(10) else
                  f"10^{c['llr'] / math.log(10):.0f}")                   # plain number, or a power of ten
        ex_txt = (f"The typical call (median evidence of the {len(pool)} calls"
                  + (" on drawn sequences" if len(pool) < len(calls) else "")
                  + f"): {s}:{c['start']:,}-{c['end']:,} ({span / 1e6:.2f} Mb), painted {group_label(c['origin'])} on a {group_label(c['own'])} sequence; {c['n_ltr']} "
                  f"informative copies inside, {c['n_lineages']} independent lineages favoring the origin, likelihood "
                  f"ratio {lr_txt} with each lineage counted once. Window: the call and its "
                  f"length again on either side ({m.size} LTR-RTs, {int(inf_m.sum())} informative; uninformative ones "
                  f"gray; opacity = confidence). Posterior: probability of each ancestry along the sequence. Hatched "
                  f"ends of the call: the breakpoints lie between the outermost copies of the origin's ancestry and the "
                  f"nearest informative copies of the sequence's own. Bottom: the smallest exchange detectable along "
                  f"the window; dashed: this call's length.")
    bp = {}
    for c_ in calls:
        k_ = f"{group_label(c_['origin'])} into {group_label(c_['own'])}"
        bp[k_] = bp.get(k_, 0) + c_["end"] - c_["start"] + 1
    above = sum(1 for c_ in calls if c_["seq"] in res and res[c_["seq"]][3].size and
                c_["end"] - c_["start"] + 1 >= np.median(res[c_["seq"]][3]))
    more = (f"The {nrow} longest of {len(names)} phased sequences" if nrow < len(names) else
            f"All {nrow} phased sequences")
    short = [(f"{name} carries {len(calls)} stretches where a chromosome of one subgenome holds another subgenome's "
              f"transposon lineages (" + "; ".join(f"{k_}: {v / 1e6:.1f} Mb" for k_, v in bp.items()) + "): "
              f"homoeologous exchanges or translocations between the subgenomes.") if calls else
             (f"{name} shows no stretch where a chromosome of one subgenome holds another subgenome's transposon "
              f"lineages strongly enough to call; " + (f"exchanges shorter than about {R / 1e6:.1f} Mb would go "
                                                        f"unseen." if rr.size else
                                                        "the smallest detectable exchange could not be measured "
                                                        "(too few informative copies).")),
             "**A** Each chromosome painted by the subgenome its transposon lineages come from.",
             "**B** The smallest exchange the transposons could reveal on each chromosome, and the exchanges found."]
    if calls:
        short.append("**C** One typical exchange up close: the transposons themselves switch ancestry at its borders.")
        short.append("**D** Every candidate exchange by length and strength of evidence: those called against those "
                     "rejected.")
    long_ = [f"Terms. **Ancestry posterior**: the probability of each subgenome's ancestry along a sequence, from a "
             f"hidden Markov model over its informative LTR-RTs (copies whose lineages tell the subgenomes apart; "
             f"switch probability {pa['switch']:.1e} per copy, fitted to the genome). **Called exchange**: a stretch "
             f"painted with another subgenome whose likelihood ratio, each lineage counted once, is >= {BF_DECISIVE:g}, "
             f"with >= {MIN_LINEAGES} lineages favoring the origin >= {BF_SUPPORT:g}-fold. **Smallest detectable "
             f"exchange**: the span of the {pa['n_needed']} informative copies whose typical evidence pays for two "
             f"ancestry switches and the call threshold.",
             f"**A** {more}, to scale; colors stack the posterior shares in bins of 1/1,500 of the longest sequence; "
             f"hatched: stretches without informative copies longer than "
             + (f"the median detectable exchange ({R / 1e6:.2f} Mb)" if rr.size else
                f"{R / 1e6:.2f} Mb (the smallest detectable exchange could not be measured)")
             + f", of unknown ancestry. Black boxes: the {len(calls)} called exchanges"
             + ("; dashed box: the window of C, marked under A and joined to C." if calls else "."),
             f"**B** Per sequence (rows as in A): gray dot, the median smallest detectable exchange; diamonds, called "
             f"exchanges by length, colored by origin ({above} of {len(calls)} at least as long as that median). "
             f"Shorter exchanges, and any in stretches without informative copies, cannot be seen from transposons "
             f"alone; the true number of exchanges is higher."]
    if calls:
        n_rej = len(pa["segs"]) - len(calls)
        long_.append(f"**C** {ex_txt}")
        long_.append(f"**D** Every segment of another subgenome's ancestry on the best path of the model "
                     f"({len(pa['segs'])}): called (diamonds, origin colors; {len(calls)}) or rejected (gray circles; "
                     f"{n_rej}: likelihood ratio below {BF_DECISIVE:g} or fewer than {MIN_LINEAGES} lineages favoring the "
                     f"origin), by length and likelihood ratio with each lineage counted once; dashed: {BF_DECISIVE:g}; "
                     f"ringed: the call in C. Its evidence is the median of the calls: it was chosen as typical, not as "
                     f"the clearest.")
    return fig, short, long_


def fig_quality(plt, a, lib, ph, bio, mu, name):
    """Figure 4, the evidence behind the calls: A every phased sequence's informative LTR-RTs against its support
    (agreement of half-libraries), assigned filled, ambiguous open; B every copy's LTR length against its age, with
    the age of one substitution per LTR (below it, ages are coarse); C every nested pair: the host's age against the
    inner element's (an insertion cannot predate its host)."""
    N, z = ph["N"], np.asarray(ph["z"])
    ek = np.asarray(ph["ek"])
    loc = np.where(ek >= 0, z[np.maximum(ek, 0)], -1) if N > 1 else np.where(ek >= 0, 0, -1)
    fig = plt.figure(figsize=(FIG_W, 3.7))
    XL, gap = 0.62, 0.66
    cw = (FIG_W - 0.12 - XL - 2 * gap) / 3
    axA, axB, axC = (_ax(fig, XL + k * (cw + gap), 0.55, cw, 2.75) for k in range(3))
    yr = lambda k: np.asarray(k, float) / (2 * mu)
    # A support against informative copies
    ninf, sup = np.asarray(ph["n_inf"], float), np.asarray(ph["support"], float)
    n_ok = n_amb = n_ins = 0
    if N > 1 and np.isfinite(sup).any():
        for g in range(N):
            m = (z == g) & np.isfinite(sup)
            ok = m & (sup >= SUPPORT_OK) & (ninf >= MIN_ELEMENTS)
            n_ok += int(ok.sum())
            n_ins += int((m & (ninf < MIN_ELEMENTS)).sum())
            n_amb += int((m & ~ok & (ninf >= MIN_ELEMENTS)).sum())
            axA.scatter(np.maximum(ninf[ok], 1), sup[ok], s=12, color=SG_COLORS[g], lw=0, zorder=3)
            axA.scatter(np.maximum(ninf[m & ~ok], 1), sup[m & ~ok], s=12, facecolors="none", edgecolors=SG_COLORS[g],
                        lw=0.7, zorder=3)
        axA.axhline(SUPPORT_OK, color=FIG_GRAY, lw=0.6, ls=(0, (3, 2)), zorder=1)
        axA.set_xscale("log")
        axA.set_xlim(max(1.0, np.nanmin(ninf) / 1.6), np.nanmax(ninf) * 1.6)
        _log_axis(axA, "x", "n")
        axA.set_ylim(min(0.0, float(np.nanmin(sup))) - 0.05, 1.05)
        axA.text(axA.get_xlim()[1], SUPPORT_OK - 0.02, "assigned ", ha="right", va="top", fontsize=6, color=FIG_GRAY)
        axA.set_ylabel("support (agreement of\nindependent library halves)")
        a_head = f"The calls for {name} rest on {np.nanmedian(ninf):,.0f} informative transposons per chromosome"
        a_short = ("**A** How much transposon evidence each chromosome's subgenome call rests on, and how "
                   "reproducible it is.")
        a_txt = (f"{n_ok} of {len(z)} phased sequences are assigned (filled; support >= {SUPPORT_OK:g}), {n_amb} "
                 f"ambiguous (open) and {n_ins} carry fewer than {MIN_ELEMENTS} informative copies; median "
                 f"{np.nanmedian(ninf):,.0f} informative copies per sequence.")
    else:
        si = {s: i for i, s in enumerate(lib.seqnames)}
        L = np.array([lib.seq_len[si[s]] for s in ph["names"]], float) / 1e6
        cnt = np.array([(ek == i).sum() for i in range(len(ph["names"]))], float)
        axA.scatter(L, cnt, s=12, color=FIG_DARK, lw=0)
        axA.set_xlabel("sequence length (Mb)")
        axA.set_ylabel("LTR-RTs per sequence")
        a_txt = (f"One subgenome: no support to show; {len(ph['names'])} phased sequences with a median of "
                 f"{np.median(cnt):,.0f} LTR-RTs each.")
        a_head = (f"The {len(ph['names'])} phased chromosomes of {name} carry a median of {np.median(cnt):,.0f} "
                  f"LTR retrotransposons each")
        a_short = "**A** How many transposons each chromosome carries."
    if N > 1 and np.isfinite(sup).any():
        axA.set_xlabel("informative LTR-RTs per sequence")
    _letter(fig, axA, "A")
    # B LTR length against age, with the one-substitution floor
    Lt = np.where(np.isfinite(lib.ltr5_len) & (lib.ltr5_len > 0), lib.ltr5_len, np.nan)
    ok = np.isfinite(Lt) & np.isfinite(lib.k2p) & (loc >= 0)
    b_txt = "No LTR lengths in the annotation: not drawn."
    if ok.sum() >= 10:
        age, Lo = yr(lib.k2p[ok]), Lt[ok]
        floor0 = float(yr(1.0 / np.median(Lo)))
        zero = age <= 0
        ya = np.where(zero, floor0 / 4, age)
        al = float(np.clip(3000 / ok.sum(), 0.03, 0.6))
        col = np.array([SG_COLORS[g] if N > 1 else FIG_DARK for g in loc[ok]])
        axB.scatter(Lo, ya, s=1.5, c=col, alpha=al, lw=0, rasterized=True)
        xl = np.geomspace(Lo.min(), Lo.max(), 50)
        axB.plot(xl, yr(1.0 / xl), color="black", lw=0.8, ls=(0, (3, 2)))
        axB.set_xscale("log")
        axB.set_yscale("log")
        axB.set_xlim(Lo.min() / 1.3, Lo.max() * 1.3)
        axB.set_ylim(floor0 / 8, max(float(ya.max()), floor0 * 10) * 2)
        _log_axis(axB, "x", "n")
        _log_axis(axB, "y", "years")
        tk = [v for v in axB.get_yticks() if floor0 / 2 < v]
        axB.set_yticks([floor0 / 4] + tk, ["0"] + [_yr_label(v) for v in tk])
        xt = float(np.sqrt(xl[0] * xl[-1]))
        (p0, q0), (p1, q1) = axB.transData.transform([(xt, float(yr(1.0 / xt))), (xt * 2, float(yr(0.5 / xt)))])
        axB.text(xt, float(yr(1.0 / xt)) / 1.25, "one substitution", ha="left", va="top", fontsize=6,
                 rotation=math.degrees(math.atan2(q1 - q0, p1 - p0)), rotation_mode="anchor")   # along the line
        axB.set_xlabel("LTR length (bp)")
        axB.set_ylabel(f"age (years ago)\nat {_rate(mu)}", fontsize=6.5)
        nz = int(zero.sum())
        b_txt = (f"{int(ok.sum()):,} copies with an annotated LTR length (median {np.median(Lo):,.0f} bp): one "
                 f"substitution in a typical LTR is {_ma(floor0)} (dashed: one substitution in an LTR of each length); "
                 f"{nz:,} copies ({100 * nz / ok.sum():.0f}%) have identical LTRs (age 0, bottom row). Below the line, "
                 f"ages rest on zero or one substitution and are coarse; the dating models this counting noise.")
    else:
        axB.axis("off")
    _letter(fig, axB, "B")
    # C nested pairs: host age against inner age
    nest = bio["nesting"]
    c_txt = "No nested elements with LTR standard errors: not drawn."
    if nest:
        hx = yr([r_["k_host"] for r_ in nest]) / 1e6
        iy = yr([r_["k_inner"] for r_ in nest]) / 1e6
        v = np.array([r_["violation"] for r_ in nest])
        alc = float(np.clip(1500 / len(nest), 0.08, 0.8))                      # dense clouds: ink follows counts
        axC.scatter(hx[~v], iy[~v], s=6, color=FIG_GRAY, lw=0, alpha=alc, zorder=2, rasterized=True)
        axC.scatter(hx[v], iy[v], s=6, color=FIG_ACCENT, lw=0, alpha=alc, zorder=3, rasterized=True)
        top = max(float(hx.max()), float(iy.max()), 1e-3) * 1.05
        axC.plot([0, top], [0, top], color=FIG_DARK, lw=0.6, ls=(0, (3, 2)), zorder=1)
        axC.set_xlim(0, top)
        axC.set_ylim(0, top)
        axC.text(0.95 * top, 0.95 * top, "same age", ha="right", va="bottom", fontsize=6, color=FIG_DARK,
                 rotation=45, rotation_mode="anchor", transform_rotates_text=True)
        axC.set_xlabel("age of the host element (Myr ago)")
        axC.set_ylabel("age of the element inside it\n(Myr ago)")
        n_in = sum(r_["in_host_ltr"] for r_ in nest)
        c_txt = (f"{len(nest):,} elements nested in another (host = the smallest element containing it). An insertion "
                 f"cannot predate its host, so points should fall below the dashed line; {int(v.sum())} "
                 f"({100 * v.mean():.1f}%; black) look older than their host by more than two standard errors: an "
                 f"upper estimate of how often single-copy ages mislead (LTR gene conversion homogenizes host LTRs; "
                 f"chimeric annotations). {n_in:,} sit inside a host LTR.")
    else:
        axC.axis("off")
    _letter(fig, axC, "C")
    rate_v = 100 * np.mean([r_["violation"] for r_ in nest]) if nest else float("nan")
    short = [f"{a_head}; ages from LTR divergence resolve to about one substitution per LTR"
             + (f", and {rate_v:.1f}% of nested elements look older than their host, the error rate of a single copy's "
                f"age." if nest else "."),
             a_short,
             "**B** Short LTRs carry few substitutions, so young ages are coarse; the dating models this counting "
             "noise.",
             "**C** Elements inside another must be younger than their host; those that are not show how often a "
             "single copy's age misleads."]
    long_ = [f"Terms. **Informative copy**: an LTR-RT whose lineages tell the subgenomes apart. **Support**: the share "
             f"of independent half-libraries of the lineages that put the sequence in the same subgenome. **Age**: LTR "
             f"divergence (corrected for multiple hits) over twice the rate, {_rate(mu)}.",
             f"**A** {a_txt}", f"**B** {b_txt}", f"**C** {c_txt}"]
    return fig, short, long_


def _diff_positions(q, t):
    """Positions in q of the substitutions on the edlib path between two code arrays (display only)."""
    r = edlib.align(ALPH[q].tobytes(), ALPH[t].tobytes(), mode="NW", task="path")
    pos, i = [], 0
    for n, op in re.findall(r"(\d+)([=XID])", r["cigar"]):
        n = int(n)
        if op == "X":
            pos += range(i, i + n)
        if op in "=XI":
            i += n
    return np.array(pos, int)


def write_figures(a, lib, ph, dt, bio, pa):
    """fig_structure, fig_dating_<SGx>-<SGy> (each subgenome pair with a dated split or merger), fig_exchanges and
    fig_quality, each <stem>.pdf (figure, short legend, long legend) and <stem>.png. A figure that fails is skipped with
    a warning: the tables are already written."""
    t0 = time.time()
    plt, family = _fig_style()
    name = f"*{a.label}*" if a.label else Path(a.ltr_fasta[0]).name.split(".")[0]
    od = Path(a.outdir)
    anc = pa["ancestry"] if pa else (dt["ancestry"] if dt else None)
    if anc is None and ph["N"] > 1:
        anc = copy_ancestry(ph["clades"], ph["ek"], ph["z"], ph["N"], np.random.default_rng([a.seed, 7]))
    jobs = [("fig_structure", lambda: fig_structure(plt, a, lib, ph, dt, pa, anc, name))]
    for r in (dt["pairs"] if dt else []):
        if "ev" in r and (r.get("split_reported") or (r["merger"] and r["merger"]["datable"])):
            jobs.append((f"fig_dating_{group_label(r['groups'][0])}-{group_label(r['groups'][1])}",
                         (lambda r_: lambda: fig_dating(plt, a, lib, ph, dt, bio, r_, name))(r)))
    if pa:
        jobs.append(("fig_exchanges", lambda: fig_exchanges(plt, a, lib, ph, pa, anc, name)))
    jobs.append(("fig_quality", lambda: fig_quality(plt, a, lib, ph, bio, dt["mu"] if dt else (a.mu or MU_DEFAULT),
                                                    name)))
    done = []
    for stem, fn in jobs:
        try:
            fig, short, long_ = fn()
            _save_figure(plt, fig, str(od / stem), short, long_, family)
            done.append(stem)
        except Exception as e:                                              # noqa: BLE001 (a figure never costs a run)
            log(f"WARNING: figure {stem} failed ({type(e).__name__}: {e}); skipped")
            plt.close("all")
    if done:
        log(f"figures: {', '.join(done)} (.png; .pdf with legend pages) in {time.time() - t0:.1f} s")
    return done


def write_report(a, lib, ph, dt, bio, pa, figs):
    """report.md: the run in plain language from its numbers (structure, dates, exchanges, age quality, files)."""
    N = ph["N"]
    name = a.label or Path(a.ltr_fasta[0]).name.split(".")[0]
    z, ek = np.asarray(ph["z"]), np.asarray(ph["ek"])
    si = {s: i for i, s in enumerate(lib.seqnames)}
    L = np.array([lib.seq_len[si[s]] for s in ph["names"]], float)
    mu = dt["mu"] if dt else (a.mu or MU_DEFAULT)
    out = [f"# v6 report: {name}", "",
           f"{lib.n:,} LTR retrotransposons on {len(lib.seqnames)} sequences; {len(ph['names'])} sequences phased into "
           f"{N} subgenome{'s' if N != 1 else ''}. Verdict: {ph['verdict']}.", "", "## Subgenomes", "",
           "| subgenome | sequences | Mb | LTR-RTs |", "|---|---:|---:|---:|"]
    for g in range(N):
        idx = np.flatnonzero(z == g)
        out.append(f"| {group_label(g)} | {idx.size} | {L[idx].sum() / 1e6:,.1f} | {int(np.isin(ek, idx).sum()):,} |")
    out += ["", "Per-sequence assignments, support and each sequence's closest partner in its own subgenome "
            "(near-identical homologues share far more lineages than the rest): chromosomes.tsv."]
    if dt:
        out += ["", "## Dates", "", f"Ages assume {_rate(mu)}"
                + ("." if dt["mu_given"] else " (the default; give --mu for your species: every age scales with it)."),
                "", "| subgenomes | progenitors split | merger |", "|---|---|---|"]
        for r in dt["pairs"]:
            sp = r["split"].get(r.get("split_reported")) if r.get("split_reported") else None
            stx = (f"{_ma(sp['estimate'] / 2 / mu)} ({_ma(sp['lo'] / 2 / mu)}-{_ma(sp['hi'] / 2 / mu)}; "
                   f"{r['split_reported']}, estimators {r.get('split_agreement')})" if sp else "not datable")
            out.append(f"| {group_label(r['groups'][0])}-{group_label(r['groups'][1])} | {stx} | "
                       f"{merger_text(r['merger'], r.get('n_sibling_pairs', 0), mu)} |")
        out += ["", "A merger date is the onset of transposition between the subgenomes; the clock's bound (\"began < "
                "X ago\") limits that onset, not the merger. The supplementary bound (\"if new copies land in either "
                "subgenome at random\") limits the merger itself, but only under that assumption. All estimates with "
                "their intervals: dating.tsv."]
    elif N > 1:
        out += ["", "## Dates", "", "Not dated (too few copy pairs)."]
    if pa:
        calls = [c for c in pa["segs"] if c["called"]]
        rr = np.concatenate([r_[3] for r_ in pa["res"]]) if pa["res"] else np.zeros(0)
        bp = {}
        for c in calls:
            k = f"{group_label(c['origin'])} into {group_label(c['own'])}"
            bp[k] = bp.get(k, 0) + c["end"] - c["start"] + 1
        out += ["", "## Exchanges", "", f"{len(calls)} exchanges called"
                + (" (" + "; ".join(f"{k}: {v / 1e6:.2f} Mb" for k, v in bp.items()) + ")" if bp else "")
                + ": segments.tsv, exchanges.bed, exchanges.gff3."
                + (f" The smallest exchange the transposons could reveal is {np.median(rr) / 1e6:.2f} Mb in the median "
                   f"(resolution.bedGraph); shorter ones, and any in stretches without informative LTR-RTs, go unseen, "
                   f"so the true number is higher. Boundaries are fuzzy: each call carries the interval its breakpoints "
                   f"lie in (left_flank, right_flank)." if rr.size else "")]
    nest = bio["nesting"]
    out += ["", "## Age quality", ""]
    if nest:
        nv = sum(r_["violation"] for r_ in nest)
        out.append(f"Of {len(nest):,} elements nested inside another, {nv} ({100 * nv / len(nest):.1f}%) look older "
                   f"than their host by more than two standard errors, which an insertion cannot be: an upper estimate "
                   f"of how often single-copy ages mislead here (LTR gene conversion, chimeric annotations). nesting.tsv.")
    else:
        out.append("No nested elements with LTR substitution counts: single-copy ages could not be checked.")
    out += ["", "## Files", ""]
    out += [f"- {f}.pdf (figure, short legend, long legend), {f}.png" for f in figs]
    out += ["- chromosomes.tsv, structure.tsv, pairing_map.tsv, links_oe.tsv: phasing"]
    if dt:
        out.append("- dating.tsv: split and merger estimates with intervals")
    out.append("- elements.tsv: every LTR-RT with its subgenome, ancestry posteriors, call and age")
    out.append("- clades.tsv, insertion_history.tsv, nesting.tsv: LTR-RT clades, insertion histories, age check")
    if pa:
        out.append("- segments.tsv, exchanges.bed, exchanges.gff3, ancestry_<SG>.bedGraph, resolution.bedGraph: "
                   "painting and exchanges" + ("; unphased_sequences.tsv: ancestry of sequences not phased"
                                               if pa.get("unphased") else ""))
    if a.features:
        out.append("- features.tsv: ancestry beneath the --features intervals")
    out.append("- summary.json: the run in machine-readable form")
    Path(a.outdir, "report.md").write_text("\n".join(out) + "\n")


def read_features(path):
    """--features BED (0-based, half-open) or GFF3 (1-based) -> [(sequence, start, end, name)], 1-based inclusive."""
    out, bad = [], 0
    with _open(path) as f:
        for ln in f:
            ln = ln.decode() if isinstance(ln, bytes) else ln
            if not ln.strip() or ln.startswith(("#", "track", "browser")):
                continue
            x = ln.rstrip("\n").split("\t")
            if len(x) >= 9 and x[3].isdigit() and x[4].isdigit() and x[6] in ("+", "-", ".", "?"):
                at = dict(kv.split("=", 1) for kv in x[8].split(";") if "=" in kv)
                out.append((x[0], int(x[3]), int(x[4]), at.get("ID", at.get("Name", f"{x[0]}:{x[3]}-{x[4]}"))))
            elif len(x) >= 3 and x[1].isdigit() and x[2].isdigit():
                out.append((x[0], int(x[1]) + 1, int(x[2]),
                            x[3] if len(x) > 3 and x[3] else f"{x[0]}:{int(x[1]) + 1}-{x[2]}"))
            else:
                bad += 1
    if not out:
        die(f"--features {path}: no BED or GFF3 intervals read")
    if bad:
        log(f"WARNING: --features: {bad} lines are neither BED nor GFF3; skipped")
    return out


def feature_ancestry(a, lib, ph, pa):
    """features.tsv: for each --features interval, the length-weighted ancestry posterior of the painting tiles it
    overlaps ('unknown' in copy-free stretches and on unphased sequences), its sequence's subgenome, and the called
    exchange it overlaps, if any."""
    feats = read_features(a.features)
    N = max(ph["N"], 1)
    sg = {s: int(g) for s, g in zip(ph["names"], ph["z"])}
    rr = np.concatenate([r[3] for r in pa["res"]]) if pa and pa["res"] else np.zeros(0)
    R = float(np.median(rr)) if rr.size else float(RES_WIN)
    tiles = {s: posterior_tiles(mid, post, L, R) for s, mid, post, L in (pa["post"] if pa else [])}
    calls = [c for c in (pa["segs"] if pa else []) if c["called"]]
    n_known = 0
    with open(Path(a.outdir) / "features.tsv", "w") as f:
        f.write("feature\tsequence\tstart\tend\tsubgenome\t" + "\t".join(f"P_{group_label(g)}" for g in range(N))
                + "\tancestry\tin_exchange\n")
        for s, st, en, nm in feats:
            v = np.full(N, np.nan)
            if s in tiles:
                ts, te, tv = tiles[s]
                ov = np.clip(np.minimum(te, en) - np.maximum(ts, st - 1), 0, None)
                if ov.sum() > 0:
                    v = (tv * ov[:, None]).sum(0) / ov.sum()
            if np.isfinite(v).all():
                anc = group_label(int(np.argmax(v)))
            elif ph["N"] == 1 and s in sg:
                anc = "SG1"                                                  # one subgenome: nothing to paint
            else:
                anc = "unknown"
            n_known += anc != "unknown"
            ex = next((f"{group_label(c['origin'])}>{group_label(c['own'])}" for c in calls
                       if c["seq"] == s and c["start"] <= en and c["end"] >= st), "no")
            f.write(f"{nm}\t{s}\t{st}\t{en}\t{group_label(sg[s]) if s in sg else 'unphased'}\t"
                    + "\t".join(_f(x_, "{:.3f}") for x_ in v) + f"\t{anc}\t{ex}\n")
    log(f"features: {len(feats):,} intervals, {n_known:,} with a known ancestry (features.tsv)")


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
        f.write("level\tn_a\tn_b\tpart_a\tpart_b\treplicability\treplicability_se\tn_pairs\tborderline\tstrength\t"
                "parent_strength\trule\taccepted\n")
        for t in ph["tree"]:
            f.write(f"{t['level']}\t{len(t['part_a'])}\t{len(t['part_b'])}\t{','.join(t['part_a'])}\t"
                    f"{','.join(t['part_b'])}\t{t['replicability']:.3f}\t{t['replicability_se']:.3f}\t{t['n_pairs']}\t"
                    f"{'yes' if t['borderline'] else 'no'}\t{t['strength']:.4f}\t{_f(t['parent_strength'])}\t"
                    f"{t['rule']}\t{'yes' if t['accepted'] else 'no'}\n")
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
                structure_tree=ph["tree"], rule=a.dev_rule, timings=timings,
                dating=dating_summary(dt) if dt else None, biology=biology_summary(bio) if bio else None,
                painting=painting_summary(pa) if pa else None)
    (od / "summary.json").write_text(json.dumps(_jsonable(summ), indent=1) + "\n")


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
    dt = None if a.dev_skip_dating else run_dating(a, lib, ph)
    t["dating_s"] = round(time.time() - t0, 1)
    t0 = time.time()
    bio = run_biology(a, lib, ph, dt)
    t["biology_s"] = round(time.time() - t0, 1)
    t0 = time.time()
    pa = run_painting(a, lib, ph, dt)
    t["painting_s"] = round(time.time() - t0, 1)
    write_phasing(a, lib, cl, ph, t, dt, bio, pa)
    if dt:
        write_dating(a, dt)
    write_biology(a, lib, bio)
    if pa:
        write_painting(a, pa)
    write_elements(a, lib, ph, dt, bio, pa)
    if a.features:
        feature_ancestry(a, lib, ph, pa)
    figs = [] if a.no_plots else write_figures(a, lib, ph, dt, bio, pa)
    write_report(a, lib, ph, dt, bio, pa, figs)
    log(f"done: {ph['N']} subgenome(s); outputs in {a.outdir} (report.md)")


if __name__ == "__main__":
    main()
