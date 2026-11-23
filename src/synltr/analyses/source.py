"""synltr source: the likely source (master) copy of each new insertion of a read genotyping run, with the evidence
to judge it.

Input: a finished reads run folder (reads.genotypes.tsv and manifest.json, which names the reference FASTA and its
annotation; read with calls.manifest.read_manifest). Output, in the same folder (whose manifest.json, the read
genotyping run's, always stays as it is) or in -o (which gets its own manifest.json; an earlier one there is removed
once the checks have passed): reads.sources.tsv (one row per new copy and tied source), reads.copies.tsv (every
candidate copy with its status, age and nearest relatives) and reads.locality.tsv (same-chromosome tests per clade
and family). Without minimap2 on PATH the run stops first (an error, exit 1).

Status, from the calls with GQ >= GQ_STATUS: new = an empty-site allele in some sample; fixed = every called sample
F/F (full on every haplotype); other = anything else; a call flagged high_depth or unexplained_clips blocks new and
fixed. A copy whose internal region, with the elements nested inside it cut out, is at least MIN_INTERNAL bp is a
candidate. The nearest fixed copies of a new copy (internal-region divergence) are its candidate sources, and ties
share weight. The search runs per family (copies without a family form one group) in two stages: minimap2 asm5
without base alignment maps every new copy against every target of its group (-N = the group's size, so no target is
crowded out), then mappy aligns the TOP_K targets with the lowest approximate divergence exactly. A copy nested in
or overlapping another is never its source. The segdup flag marks a copy sharing its source's TSD site (the TSD
compared by `canon_tsd`, so either strand) or a flank identity of at least FLANK_DUP_ID.

Coordinates: elements are the reads route's Targets (0-based half-open); the tables print 1-based closed starts.

The two search stages are cached in <outdir>/work/source/ under a key of the candidates (id, status, family,
sequence), MIN_ALN, MIN_COV, MIN_COV_LONG, TOP_K, the minimap2 and mappy versions and the code fingerprint of this
module's import closure, so a change in any of them recomputes the alignments and nothing else does.

Table cells: a bool column prints True/False, an integer column no decimal point, a column holding a NaN prints
floats with the NaN empty.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import logging
import math
import shutil
import subprocess
from collections import Counter
from pathlib import Path

import mappy
import numpy as np
import parasail

from synltr.base.checkpoint import (
    atomic_path,
    atomic_write,
    clear_done,
    code_fingerprint,
    is_done,
    mark_done,
    marker,
    stage_key,
    tool_version,
)
from synltr.base.constants import (
    CONS_FRAC,
    CONS_MIN,
    DIV_LIMITS,
    EXACT_CACHE,
    EXACT_CHUNK,
    FLANK,
    FLANK_DUP_ID,
    GQ_STATUS,
    LOCALITY_MIN_NEW,
    LOCALITY_SEED,
    MANIFEST,
    MIN_ALN,
    MIN_COV,
    MIN_COV_LONG,
    MIN_INTERNAL,
    MOSAIC_ALPHA,
    MOSAIC_COV,
    MOSAIC_GAIN,
    N_NEAREST,
    NEAREST_SPARE,
    OLDER_DIFFS,
    PERM,
    READS_STEM,
    SOURCE_DIR,
    SW_SCORE,
    TIE_DIFFS,
    TOP_K,
)
from synltr.base.log import SynltrError, UsageError
from synltr.base.runtime import pool_map
from synltr.calls.manifest import Finished, manifest_samples, read_manifest, remove_manifest
from synltr.io.annotation.tsd import canon_tsd
from synltr.io.seqio import Fasta, open_fasta
from synltr.routes.reads.reads import Target, sw_identity
from synltr.routes.reads.run import load_targets
from synltr.seq.divergence import k2p

log = logging.getLogger(__name__)

SRC_COLS = ["new_id", "new_chrom", "new_start", "new_end", "family", "clade", "source_id", "source_chrom",
            "source_start", "source_end", "divergence", "aln_len", "est_diffs", "n_tied", "weight", "same_chrom",
            "distance", "nearest_any_id", "nearest_any_status", "nearest_any_div", "new_ltr_diffs",
            "source_ltr_diffs", "source_older", "same_tsd", "flank_identity", "shared_derived", "source_private",
            "new_private", "consensus_sites", "mosaic_left", "mosaic_right", "mosaic_breakpoint", "mosaic_support",
            "mosaic_conflicts", "mosaic_p_adj", "flags"]
COPY_COLS = ["element_id", "chrom", "start", "end", "family", "clade", "status", "n_called", "internal_len",
             "ltr_diffs", "ltr_k2p", "nearest_any_id", "nearest_any_div", "nearest_fixed_id", "nearest_fixed_div"]
LOC_COLS = ["max_divergence", "level", "group", "n_new", "same", "expected", "obs_exp", "p", "perm_expected",
            "perm_obs_exp", "perm_p", "median_distance_mb"]

_SC = parasail.matrix_create("ACGT", SW_SCORE[0], SW_SCORE[1])
_ACGTN_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")  # rc's table: a letter other than ACGTN stays as it is
TS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
NAN = float("nan")
REMEDY = "run synltr genotype on read sets there first"  # what a folder without a finished read genotyping run needs


def rc(s: str) -> str:
    return s.translate(_ACGTN_COMPLEMENT)[::-1]


def same_site(a: str | None, b: str | None) -> bool:
    """True if two TSDs are the same site: both are TSDs and equal in canonical orientation (`canon_tsd`)."""
    ca = canon_tsd(a)
    return bool(ca) and ca == canon_tsd(b)


# ---- the run, its calls and its elements --------------------------------------------------------------------------
def reference_of(run: Path) -> tuple[Path, str]:
    """The reference FASTA and annotation of a reads run folder, from its manifest.json (calls.manifest.read_manifest:
    the row of role reference)."""
    run = Path(run)
    doc = read_manifest(run, REMEDY, absent=f"{run} is not a finished read genotyping run: no {MANIFEST} (synltr "
                                             "genotype on read sets writes it)")
    try:
        (ref,) = [s for s in doc["samples"] if s.get("role") == "reference"]
        return Path(ref["fasta"]), ref["annotation"]
    except (ValueError, KeyError, TypeError, AttributeError):
        raise UsageError(f"{run / MANIFEST}: no sample of role reference with a fasta and an annotation") from None


def read_status(path: Path):
    """Element status from the genotype table at `path` (reads.genotypes.tsv; the chrom, start and end columns
    are optional). Returns (status by element ID, called samples by ID, (chrom, start, end) by ID or None); a call
    counts when its GQ is at least GQ_STATUS."""
    need = {"element_id", "sample", "genotype", "gq", "flags"}
    with open(path, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        if not need <= set(rd.fieldnames or ()):
            raise UsageError(f"{path}: not a genotype table, missing column(s) {sorted(need - set(rd.fieldnames or ()))}")
        has_coords = {"chrom", "start", "end"} <= set(rd.fieldnames)
        agg: dict[str, list] = {}
        coords: dict[str, tuple] = {}
        seen = []
        for r in rd:
            eid = r["element_id"]
            if eid not in coords:
                seen.append(eid)
                coords[eid] = (r["chrom"], int(r["start"]), int(r["end"])) if has_coords else None
            if float(r["gq"]) < GQ_STATUS:
                continue
            gt = r["genotype"]
            a = agg.setdefault(eid, [0, False, True, False])
            a[0] += 1
            a[1] |= "E" in gt
            a[2] &= all(x == "F" for x in gt.split("/"))
            a[3] |= "high_depth" in r["flags"] or "unexplained_clips" in r["flags"]
    status = {}
    for eid in seen:
        _, any_e, all_f, bad = agg.get(eid, (0, False, False, False))
        status[eid] = "new" if any_e and not bad else "fixed" if all_f and not bad and eid in agg else "other"
    ncalled = {k: v[0] for k, v in agg.items()}
    return status, ncalled, (coords if has_coords else None)


def own_internal(fa: Fasta, el: Target) -> str:
    """Internal region of an element with the elements nested inside it excised (its own sequence)."""
    s, e = el.internal
    keep, cur = [], s
    for a, b in sorted(el.inner):
        if b <= s or a >= e:
            continue
        if a > cur:
            keep.append((cur, a))
        cur = max(cur, b)
    if cur < e:
        keep.append((cur, e))
    return "".join(fa.fetch(el.chrom, a, b) for a, b in keep).upper()


def nested_pair(a: Target, b: Target) -> bool:
    """True if one element contains or overlaps the other (a host cannot be its own insertion's source)."""
    return a.chrom == b.chrom and a.start < b.end and b.start < a.end


def ltr_diffs(fa: Fasta, el: Target) -> tuple:
    """Differences between the two LTRs of one element by global alignment:
    (substitutions, transitions, indel events, aligned columns without gaps)."""
    a = fa.fetch(el.chrom, *el.lltr).upper()
    b = fa.fetch(el.chrom, *el.rltr).upper()
    if not a or not b:
        return (NAN, NAN, NAN, 0)
    tb = parasail.nw_trace_scan_16(a, b, SW_SCORE[2], SW_SCORE[3], _SC).traceback
    sub = ts = indel = n = 0
    prev_gap = False
    for x, y in zip(tb.query, tb.ref, strict=False):
        if x == "-" or y == "-":
            if not prev_gap:
                indel += 1
            prev_gap = True
            continue
        prev_gap = False
        n += 1
        if x != y:
            sub += 1
            ts += (x, y) in TS
    return (sub, ts, indel, n)


def ltr_k2p(sub, ts, n) -> float:
    """The K2P distance between an element's LTRs (NaN if undefined), in the product form of K2P (never -0.0)."""
    if not n or math.isnan(sub):
        return NAN
    return float(k2p([n], [ts], [sub - ts], form="product")[0][0])


# ---- alignments --------------------------------------------------------------------------------------------------
def _cs_ops(cs):
    """Parse a minimap2 short cs string into (op, value) pairs."""
    i, n = 0, len(cs)
    while i < n:
        op = cs[i]
        j = i + 1
        if op == ":":
            while j < n and cs[j].isdigit():
                j += 1
            yield ":", int(cs[i + 1:j])
        elif op == "*":
            yield "*", (cs[i + 1], cs[i + 2])
            j = i + 3
        else:  # + or -
            while j < n and cs[j].isalpha():
                j += 1
            yield op, cs[i + 1:j]
        i = j


def bases_on_target(hit, qseq):
    """Map target positions -> query base ('-' deletion) for one mappy hit (forward-strand coordinates)."""
    out = {}
    t = hit.r_st
    for op, val in _cs_ops(hit.cs):
        if op == ":":
            t += val
        elif op == "*":
            out[t] = val[1].upper()
            t += 1
        elif op == "-":
            for k in range(len(val)):
                out[t + k] = "-"
            t += len(val)
    return out   # positions absent from 'out' but inside [r_st, r_en) carry the target base


def shared_mutations(m_seq, n_seq, others):
    """Shared derived substitutions between source M and new copy N, polarized by the consensus of other relatives.
    Returns (shared, source_private, new_private, informative_sites)."""
    al = mappy.Aligner(seq=m_seq, preset="asm5", best_n=1)

    def one(q):
        h = next(iter(al.map(q, cs=True)), None)
        if h is None:
            return None
        return h, bases_on_target(h, q)
    n = one(n_seq)
    if n is None:
        return (NAN,) * 4
    hn, bn = n
    obs = [o for o in (one(s) for s in others) if o is not None]
    if len(obs) < CONS_MIN:
        return (NAN,) * 4
    shared = sprivate = nprivate = informative = 0
    var_sites = set(bn)
    for h, b in obs:
        var_sites |= set(b)
    for t in sorted(var_sites):
        mb = m_seq[t].upper()
        calls = [b.get(t, mb) for h, b in obs if h.r_st <= t < h.r_en]
        if len(calls) < CONS_MIN:
            continue
        cons, k = Counter(calls).most_common(1)[0]
        if k < max(CONS_MIN, CONS_FRAC * len(calls)) or cons == "-":
            continue
        if not (hn.r_st <= t < hn.r_en):
            continue
        nb = bn.get(t, mb)
        if nb == "-" or mb == "-":
            continue
        informative += 1
        if mb != cons:                         # derived in the source lineage
            if nb == mb:
                shared += 1
            elif nb == cons:
                sprivate += 1
        elif nb != mb:
            nprivate += 1
    return shared, sprivate, nprivate, informative


def scan_p(labels):
    """Exact one-breakpoint scan test (MaxChi-like). labels: +1 where the new copy follows c1, -1 where it follows c2,
    ordered along the sequence. Under no recombination the labels are randomly ordered; the statistic is the maximum
    h of the running sum, with P(max >= h) = C(n, n1 - h) / C(n, n1) (reflection principle). Returns (p, h, k)."""
    n, n1 = len(labels), labels.count(1)
    w = h = k = 0
    for i, x in enumerate(labels, 1):
        w += x
        if w > h:
            h, k = w, i
    if n1 - h < 0:
        return 0.0, h, k
    return math.comb(n, n1 - h) / math.comb(n, n1), h, k


def mosaic_test(n_seq, cands):
    """Two-parent (template-switch) test. cands: dict id -> internal sequence of candidate parents (nearest relatives
    of any status and the tied sources). Each candidate is aligned to the new copy N (positions where it differs).
    For every ordered pair (c1, c2), the diagnostic sites (exactly one of c1, c2 differs from N) are labeled by the
    parent N follows there; scan_p tests for one breakpoint with c1 left and c2 right. Bonferroni over pairs. Called
    recombinant if adjusted P <= MOSAIC_ALPHA and the mosaic explains >= MOSAIC_GAIN more differences than the best
    single candidate. Returns (c1, c2, breakpoint, support_left, support_right, conflicts, p_adj) or None."""
    al = mappy.Aligner(seq=n_seq, preset="asm5", best_n=1)
    L = len(n_seq)
    diff = {}
    for cid, q in cands.items():
        h = next(iter(al.map(q, cs=True)), None)
        if h is None or (h.r_en - h.r_st) < MOSAIC_COV * L:
            continue
        diff[cid] = (h.r_st, h.r_en, set(bases_on_target(h, q)))
    if len(diff) < 2:
        return None
    ids = list(diff)
    m_tests = len(ids) * (len(ids) - 1)
    best = None
    for c1 in ids:
        for c2 in ids:
            if c2 == c1:
                continue
            s1, e1, d1 = diff[c1]
            s2, e2, d2 = diff[c2]
            lo, hi = max(s1, s2), min(e1, e2)
            sites = sorted([(x, 1) for x in d2 - d1 if lo <= x < hi] + [(x, -1) for x in d1 - d2 if lo <= x < hi])
            if len(sites) < 2:
                continue
            lab = [v for _, v in sites]
            p, h, k = scan_p(lab)
            if k == 0 or k == len(lab):
                continue
            b = (sites[k - 1][0] + sites[k][0]) // 2 + 1
            single = min(len([x for x in d if lo <= x < hi]) for (s_, e_, d) in diff.values() if s_ <= lo and e_ >= hi)
            mos = sum(1 for x in d1 if lo <= x < b) + sum(1 for x in d2 if b <= x < hi)
            if single - mos < MOSAIC_GAIN:
                continue
            sl, cl = lab[:k].count(1), lab[:k].count(-1)
            sr, cr = lab[k:].count(-1), lab[k:].count(1)
            padj = min(1.0, p * m_tests)
            if best is None or padj < best[6]:
                best = (c1, c2, b, sl, sr, cl + cr, padj)
    return best if best is not None and best[6] <= MOSAIC_ALPHA else None


_X: dict = {}


def exact_div(tseq, qseq, al):
    """Exact divergence of qseq against tseq: (substitutions + indel events) / aligned columns, with the coverage
    rules of a relative; None if they fail."""
    h = next(iter(al.map(qseq, cs=True)), None)
    if h is None:
        return None
    span = min(h.q_en - h.q_st, h.r_en - h.r_st)
    if h.blen < MIN_ALN or span / min(len(qseq), len(tseq)) < MIN_COV or span / max(len(qseq), len(tseq)) < MIN_COV_LONG:
        return None
    d = sum(1 for op, _ in _cs_ops(h.cs) if op in "*+-")
    return d / h.blen, h.blen


def _exact_chunk(chunk):
    """The exact hits of a list of (query, targets) jobs, in job order: (query, target, divergence, aligned bp)."""
    seqs, cache = _X["seqs"], _X["al"]
    out = []
    for q, ts in chunk:
        for t in ts:
            if t == q:
                continue
            if t not in cache:
                if len(cache) > EXACT_CACHE:          # bounded per-worker cache of aligner indexes
                    cache.clear()
                cache[t] = mappy.Aligner(seq=seqs[t], preset="asm5", best_n=1)
            r = exact_div(seqs[t], seqs[q], cache[t])
            if r is not None:
                out.append((q, t, round(r[0], 6), r[1]))
    return out


def screen_paf(path):
    """Approximate hits (no base alignment): (query, target, approximate divergence dv)."""
    with open(path) as fh:
        for line in fh:
            x = line.rstrip("\n").split("\t")
            if x[0] == x[5]:
                continue
            dv = next((float(f[5:]) for f in x[12:] if f.startswith("dv:f:")), 1.0)
            yield x[0], x[5], dv


def screen(work: Path, fams: dict, seqs: dict, status: dict, target: str, threads: int) -> dict:
    """Stage 1 for the target set "fixed" (the fixed copies of each family) or "all" (every candidate of it): minimap2
    asm5 without base alignment of each family's new copies against its targets, -N = the number of targets. Returns
    {query: [(approximate divergence, target)]}."""
    pairs: dict[str, list] = {}
    for fam, lst in sorted(fams.items()):
        new = [e for e in lst if status[e.id] == "new"]
        tg = [e for e in lst if status[e.id] == "fixed"] if target == "fixed" else lst
        if not new or not tg or (target == "all" and len(tg) < 2):
            continue
        qf, tf, pf = work / "q.fa", work / "t.fa", work / "part.paf"
        for path, ss in ((qf, new), (tf, tg)):
            with open(path, "w") as fh:
                fh.writelines(f">{e.id}\n{seqs[e.id]}\n" for e in ss)
        with open(pf, "w") as fo, open(f"{pf}.log", "w") as fe:
            r = subprocess.run(["minimap2", "-x", "asm5", "-p", "0", "-t", str(threads), "-N", str(len(tg)),
                                str(tf), str(qf)], stdout=fo, stderr=fe, stdin=subprocess.DEVNULL, check=False)
        if r.returncode:
            raise SynltrError(f"minimap2 failed (exit {r.returncode}); see {pf}.log")
        for q, t, dv in screen_paf(pf):
            pairs.setdefault(q, []).append((dv, t))
    return pairs


def search(work: Path, fams: dict, seqs: dict, status: dict, target: str, threads: int, out: Path) -> int:
    """Both stages for one target set, the exact hits written to `out` (q, t, de, blen). Returns the number of new
    copies searched."""
    pairs = screen(work, fams, seqs, status, target, threads)
    jobs = [(q, [t for _, t in sorted(v)[:TOP_K]]) for q, v in pairs.items()]
    chunks = [jobs[i:i + EXACT_CHUNK] for i in range(0, len(jobs), EXACT_CHUNK)]
    _X["seqs"], _X["al"] = seqs, {}
    try:
        done = pool_map(_exact_chunk, chunks, threads, warm=False)
    finally:
        _X.clear()
    with atomic_path(out) as tmp, open(tmp, "w") as fo:
        fo.write("q\tt\tde\tblen\n")
        for res in done:
            fo.writelines(f"{q}\t{t}\t{de}\t{blen}\n" for q, t, de, blen in res)
    log.debug("search new vs %s: %d new copies, %d exact alignments", target, len(jobs),
              sum(len(t) for _, t in jobs))
    return len(jobs)


def read_exact(path: Path) -> list[tuple]:
    """The exact hits of `path` as (query, target, divergence, aligned bp), the lowest divergence of each query first
    (ties in file order), one row per query and target."""
    with open(path) as fh:
        rows = [(r["q"], r["t"], float(r["de"]), int(r["blen"])) for r in csv.DictReader(fh, delimiter="\t")]
    rows.sort(key=lambda r: (r[0], r[2]))
    seen, out = set(), []
    for r in rows:
        if (r[0], r[1]) not in seen:
            seen.add((r[0], r[1]))
            out.append(r)
    return out


# ---- the run -----------------------------------------------------------------------------------------------------
def main(args: argparse.Namespace) -> int | Finished:
    """`synltr source RUN`: args.run, args.outdir (None: RUN) and args.threads. Into RUN, returns the exit code:
    RUN's manifest.json, the read genotyping run's, stays as it is. Into another folder, returns the record of the
    source run there (calls.manifest.Finished: the samples of RUN's manifest), which the command line writes as that
    folder's manifest.json."""
    genotype_dir = Path(args.run)
    outdir = genotype_dir if args.outdir is None else Path(args.outdir)
    code = run(genotype_dir, outdir, args.threads)
    if outdir.resolve() == genotype_dir.resolve():
        return code
    samples = manifest_samples(read_manifest(genotype_dir, REMEDY))
    return Finished(outdir, tuple(samples))


def run(genotype_dir: Path, outdir: Path, threads: int) -> int:
    """Likely source copies of the new insertions of the reads run in `genotype_dir`, written to `outdir`, with
    `threads` in all (minimap2 threads, then pool workers: never both at once). Returns 0. A missing minimap2 is an
    SynltrError (exit 1), a folder that holds no finished read genotyping run a UsageError; once these checks pass,
    an `outdir` other than `genotype_dir` loses its earlier manifest.json (the run's own record always stays)."""
    genotype_dir, outdir = Path(genotype_dir), Path(outdir)
    if shutil.which("minimap2") is None:
        raise SynltrError("minimap2 not found on PATH (needed by synltr source); install it with the conda "
                            "environment in environment.yml (bioconda's minimap2)")
    reference, annotation = reference_of(genotype_dir)
    table = genotype_dir / f"{READS_STEM}.genotypes.tsv"
    if not table.is_file():
        raise UsageError(f"{genotype_dir} is not a finished read genotyping run: no {table.name}")
    if outdir.resolve() != genotype_dir.resolve():
        remove_manifest(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    work = outdir / "work" / SOURCE_DIR
    work.mkdir(parents=True, exist_ok=True)
    with open_fasta(reference, "reference", outdir / "work" / "fai") as fa:
        targets, _ = load_targets(annotation, fa, "reference")
        els = {e.id: e for e in targets}
        status, ncalled, coords = read_status(table)
        shared = [k for k in status if k in els]
        if not shared:
            raise UsageError(f"no element ID of {table.name} is in {annotation}: give the annotation used for "
                             "genotyping")
        moved = [] if coords is None else [k for k in shared if (els[k].chrom, els[k].start + 1, els[k].end)
                                           != coords[k]]
        if moved:
            raise UsageError(f"{len(moved)} element(s) have other coordinates in {annotation} than in {table.name}, "
                             f"e.g. {moved[0]}: give the annotation used for genotyping")
        status = {k: status.get(k, "other") for k in els}
        seqs = {e.id: own_internal(fa, e) for e in els.values()}
        cand = [e for e in els.values() if len(seqs[e.id]) >= MIN_INTERNAL]
        log.info("source: %d elements; status %s; candidates (internal >= %d bp) %d", len(els),
                 dict(Counter(status.values())), MIN_INTERNAL, len(cand))
        groups = {"new": [e for e in cand if status[e.id] == "new"],
                  "fixed": [e for e in cand if status[e.id] == "fixed"]}
        if not groups["new"] or not groups["fixed"]:
            log.warning("no source to assign: %d new and %d fixed copies with an internal region >= %d bp (status "
                        "from calls with GQ >= %d); tables left empty", len(groups["new"]), len(groups["fixed"]),
                        MIN_INTERNAL, GQ_STATUS)
            for suf, cols in (("sources", SRC_COLS), ("copies", COPY_COLS), ("locality", LOC_COLS)):
                write_table(outdir / f"{READS_STEM}.{suf}.tsv", cols, [])
            return 0
        fams: dict[str, list] = {}
        for e in cand:
            fams.setdefault(e.family if e.family != "." else "_unlabeled", []).append(e)
        digest = hashlib.sha256(";".join(f"{e.id}:{status[e.id]}:{e.family}:{hashlib.sha1(seqs[e.id].encode()).hexdigest()}"
                                         for e in sorted(cand, key=lambda x: x.id)).encode()).hexdigest()
        params = {"search": "two-stage", "MIN_ALN": MIN_ALN, "MIN_COV": MIN_COV, "MIN_COV_LONG": MIN_COV_LONG,
                  "TOP_K": TOP_K, "minimap2": tool_version("minimap2"), "mappy": mappy.__version__}
        for tgt in ("fixed", "all"):
            out = work / f"new_vs_{tgt}.tsv"
            key = stage_key(f"source.new_vs_{tgt}", {"candidates": digest}, params, code_fingerprint("analyses.source"))
            if is_done(marker(out), key):
                log.debug("search new vs %s: reused", tgt)
                continue
            clear_done(marker(out))
            search(work, fams, seqs, status, tgt, threads, out)
            mark_done(marker(out), key, [out])
        fx = [r for r in read_exact(work / "new_vs_fixed.tsv") if not nested_pair(els[r[0]], els[r[1]])]
        al = [r for r in read_exact(work / "new_vs_all.tsv") if not nested_pair(els[r[0]], els[r[1]])]
        src, cp = tables(fa, els, cand, groups, status, ncalled, seqs, fx, al)
        write_table(outdir / f"{READS_STEM}.sources.tsv", SRC_COLS, src)
        write_table(outdir / f"{READS_STEM}.copies.tsv", COPY_COLS, cp)
        loc = locality(src, fa, cand)
        write_table(outdir / f"{READS_STEM}.locality.tsv", LOC_COLS, loc)
    log.info("source: new copies with a fixed relative: %d; pairs %d", len({r["new_id"] for r in src}), len(src))
    for r in loc:
        log.debug("locality: %s", " ".join(f"{k}={v}" for k, v in r.items()))
    log.info("source: done -> %s.{sources,copies,locality}.tsv in %s", READS_STEM, outdir)
    return 0


def tables(fa, els, cand, groups, status, ncalled, seqs, fx, al):
    """The sources and copies tables as lists of rows, from the exact hits of the new copies against the fixed copies
    (`fx`) and against every candidate (`al`)."""
    al_status = [(q, t, de, blen, status[t]) for q, t, de, blen in al]
    ages = {}
    for eid in {e.id for e in groups["new"]} | {r[1] for r in fx}:
        sub, ts, _, n = ltr_diffs(fa, els[eid])
        ages[eid] = (sub, ts, ltr_k2p(sub, ts, n))
    near_any: dict = {}
    al_q: dict = {}
    for r in al_status:
        near_any.setdefault(r[0], r)
        al_q.setdefault(r[0], []).append(r)
    fx_q: dict = {}
    for r in fx:
        fx_q.setdefault(r[0], []).append(r)
    rows = []
    for e in groups["new"]:
        h = fx_q.get(e.id)
        if not h:
            continue
        best = min(r[2] for r in h)
        tied = [r for r in h if (r[2] - best) * r[3] <= TIE_DIFFS]
        na = near_any.get(e.id)
        oth = al_q.get(e.id, [])[:N_NEAREST + NEAREST_SPARE]
        mcand = {t: seqs[t] for t in [r[1] for r in oth][:N_NEAREST] + [r[1] for r in tied] if t != e.id}
        mos = mosaic_test(seqs[e.id], mcand)
        for _, t, de, blen in tied:
            m = els[t]
            others = [seqs[x[1]] for x in oth if x[1] not in (e.id, m.id)][:N_NEAREST]
            sh = shared_mutations(seqs[m.id], seqs[e.id], others)
            fid = flank_identity(fa, e, m)
            same_tsd = same_site(e.tsd, m.tsd)
            closer = na is not None and na[1] != m.id and (de - na[2]) * blen >= 1 and na[4] != "fixed"
            flags = []
            if same_tsd or fid >= FLANK_DUP_ID:
                flags.append("segdup")
            if closer:
                flags.append("closer_nonfixed_relative")
            if len(tied) > 1:
                flags.append("tied")
            if mos is not None:
                flags.append("recombinant")
            older = ages[m.id][0] - ages[e.id][0] >= OLDER_DIFFS if not math.isnan(ages[m.id][0]) else False
            rows.append({
                "new_id": e.id, "new_chrom": e.chrom, "new_start": e.start + 1, "new_end": e.end,
                "family": e.family, "clade": e.clade,
                "source_id": m.id, "source_chrom": m.chrom, "source_start": m.start + 1, "source_end": m.end,
                "divergence": round(de, 5), "aln_len": blen, "est_diffs": round(de * blen, 1), "n_tied": len(tied),
                "weight": round(1 / len(tied), 4), "same_chrom": e.chrom == m.chrom,
                "distance": abs((e.start + e.end) / 2 - (m.start + m.end) / 2) if e.chrom == m.chrom else NAN,
                "nearest_any_id": na[1] if na is not None else ".",
                "nearest_any_status": na[4] if na is not None else ".",
                "nearest_any_div": float(np.round(na[2], 5)) if na is not None else NAN,  # numpy rounding
                "new_ltr_diffs": ages[e.id][0], "source_ltr_diffs": ages[m.id][0], "source_older": older,
                "same_tsd": same_tsd, "flank_identity": round(fid, 3), "shared_derived": sh[0],
                "source_private": sh[1],
                "new_private": sh[2], "consensus_sites": sh[3],
                "mosaic_left": mos[0] if mos else ".", "mosaic_right": mos[1] if mos else ".",
                "mosaic_breakpoint": mos[2] if mos else NAN, "mosaic_support": f"{mos[3]},{mos[4]}" if mos else ".",
                "mosaic_conflicts": mos[5] if mos else NAN, "mosaic_p_adj": mos[6] if mos else NAN,
                "flags": ",".join(flags) or "PASS"})
    cp = []
    for e in cand:
        na = near_any.get(e.id)
        h = fx_q.get(e.id)
        cp.append({
            "element_id": e.id, "chrom": e.chrom, "start": e.start + 1, "end": e.end, "family": e.family,
            "clade": e.clade, "status": status[e.id], "n_called": int(ncalled.get(e.id, 0)),
            "internal_len": e.internal[1] - e.internal[0],
            "ltr_diffs": ages.get(e.id, (NAN,))[0], "ltr_k2p": ages.get(e.id, (NAN,) * 3)[2],
            "nearest_any_id": na[1] if na is not None else ".", "nearest_any_div": na[2] if na is not None else NAN,
            "nearest_fixed_id": h[0][1] if h else ".", "nearest_fixed_div": h[0][2] if h else NAN})
    return rows, cp


def flank_identity(fa: Fasta, a: Target, b: Target) -> float:
    """Best SW identity between the FLANK-bp flanks of two copies (both orientations)."""
    def fl(e):
        L = fa.length(e.chrom)
        t = len(e.tsd) if e.tsd else 0
        return (fa.fetch(e.chrom, max(0, e.start - t - FLANK), max(0, e.start - t)).upper(),
                fa.fetch(e.chrom, min(L, e.end + t), min(L, e.end + t + FLANK)).upper())
    a1, b1 = fl(a)
    a2, b2 = fl(b)
    best = 0.0
    for x, y in ((a1, a2), (b1, b2), (a1, rc(b2)), (b1, rc(a2))):
        best = max(best, sw_identity(x, y))
    return best


# ---- locality ----------------------------------------------------------------------------------------------------
def locality_test(d, share, rng, n_perm=PERM):
    """Same-chromosome test for one set of (new copy, source) pairs (rows of the sources table) with tie weights.
    bp null: each new copy shares its source's chromosome with probability = its own chromosome's bp share
             (Poisson-binomial; assumes independent new copies).
    permutation null: the new copies' chromosome labels are shuffled among new copies, sources and their offspring
             clusters kept; controls for shared sources and for chromosomes rich in both sources and new copies."""
    ids = sorted({r["new_id"] for r in d})
    pos = {i: k for k, i in enumerate(ids)}
    same = np.zeros(len(ids))
    chrom = [None] * len(ids)
    for r in d:
        k = pos[r["new_id"]]
        same[k] += r["same_chrom"] * r["weight"]
        if chrom[k] is None:
            chrom[k] = r["new_chrom"]
    obs = float(same.sum())
    nul = np.array([share.get(c, 0.0) for c in chrom])
    p_bp = float(np.mean((rng.random((n_perm, len(ids))) < nul).sum(1) >= obs - 1e-9))
    lab = np.array(chrom)
    idx = np.array([pos[r["new_id"]] for r in d])
    sc = np.array([r["source_chrom"] for r in d])
    w = np.array([r["weight"] for r in d])
    null = np.empty(n_perm)
    for k in range(n_perm):
        perm = lab[rng.permutation(len(lab))]
        null[k] = float(np.sum(w * (perm[idx] == sc)))
    p_perm = float(np.mean(null >= obs - 1e-9))
    return {"n_new": len(ids), "same": round(obs, 2), "expected": round(nul.sum(), 2),
            "obs_exp": round(obs / nul.sum(), 2) if nul.sum() else NAN, "p": p_bp,
            "perm_expected": round(float(null.mean()), 2),
            "perm_obs_exp": round(obs / null.mean(), 2) if null.mean() else NAN, "perm_p": p_perm}


def locality(src, fa, cand, seed=LOCALITY_SEED):
    """Same-chromosome locality per clade (and per family with >= LOCALITY_MIN_NEW new copies) for cumulative
    source-divergence limits (a distant nearest fixed copy is not a plausible source). Two nulls: see locality_test.
    Rows are sorted by level, group and limit."""
    if not src:
        return []
    chroms = sorted({e.chrom for e in cand})
    lens = np.array([fa.length(c) for c in chroms], dtype=float)
    share = dict(zip(chroms, lens / lens.sum(), strict=True))
    ok = [dict(r) for r in src if "segdup" not in r["flags"]]
    size = Counter(r["new_id"] for r in ok)
    for r in ok:
        r["weight"] = 1.0 / size[r["new_id"]]
    rng = np.random.default_rng(seed)
    out = []
    for maxd in DIV_LIMITS:
        sel = [r for r in ok if r["divergence"] <= maxd]
        for level in ("clade", "family"):
            for g in sorted({r[level] for r in sel}):
                d = [r for r in sel if r[level] == g]
                res = locality_test(d, share, rng)
                dd = [r["distance"] for r in d if r["same_chrom"]]
                out.append(dict(max_divergence=maxd, level=level, group=g, **res,
                                median_distance_mb=round(np.median(dd) / 1e6, 2) if dd else NAN))
    out = [r for r in out if r["level"] == "clade" or r["n_new"] >= LOCALITY_MIN_NEW]
    out.sort(key=lambda r: (r["level"], r["group"], r["max_divergence"]))
    return out


# ---- tables ------------------------------------------------------------------------------------------------------
def _cells(values: list) -> list[str]:
    """One column's cells: bool as True/False; ints without a decimal point; a column
    with any float (a NaN included) as floats with the NaN empty; text as is."""
    if values and all(isinstance(v, (bool, np.bool_)) for v in values):
        return ["True" if v else "False" for v in values]
    num = (int, float, np.integer, np.floating)
    if values and all(isinstance(v, num) and not isinstance(v, (bool, np.bool_)) for v in values):
        if any(isinstance(v, (float, np.floating)) for v in values):
            return ["" if math.isnan(v) else repr(float(v)) for v in values]
        return [str(int(v)) for v in values]
    return ["" if isinstance(v, float) and math.isnan(v) else str(v) for v in values]


def write_table(path: Path, columns: list[str], rows: list[dict]) -> None:
    """Write `rows` (dicts) as a tab-separated table with `columns`, a header line only when there are none."""
    cols = [_cells([r[c] for r in rows]) for c in columns]
    lines = ["\t".join(columns)] + ["\t".join(c[i] for c in cols) for i in range(len(rows))]
    atomic_write(path, "\n".join(lines) + "\n")
