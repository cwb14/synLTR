"""One read set: evidence at every element, the read model, the calls and its rows of ltrkin's genotypes.tsv.

`genotype_sample` runs a read set: its evidence (`scan`), the annotation boundary error learned from it, the model's
three passes (calls at a fixed error rate, then the sample's own error rate and depth at element sites, then its
learned genotype prior), and its genotypes.tsv lines (one per element on the BAM's sequences, in the order of the
route's elements). It returns them with the calls the reference VCF needs; `write_genotypes` writes the run's table,
every read set's lines in turn, as ltrkin wrote a run of several BAMs.

The evidence is a checkpoint: <outdir>/reads.<S>.evidence.json.gz (the depths at the sampled points, the read spans
and each element's counts; S the read set's label) with its done marker, keyed by the BAM (base.checkpoint.bam_key:
size, modification time and header), the reference FASTA's checksum, --min-mapq, --ploidy, every element (ID,
sequence, extent, LTRs, TSD) and the code of this module's import closure. Workers gather it EVIDENCE_CHUNK elements
at a time (base.runtime.pool_map, each worker with its own BAM and FASTA handles) and it comes back in element order
whatever the number of workers, so the sums over it (the error and depth estimates) and their exact ties do not
depend on which worker finishes first.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import logging
import math
import os
import time
from dataclasses import dataclass
from pathlib import Path
from statistics import median

import pysam

from synltr.base.checkpoint import (
    atomic_path,
    bam_key,
    clear_done,
    code_fingerprint,
    file_checksum,
    is_done,
    mark_done,
    marker,
    stage_key,
)
from synltr.base.constants import (
    BOUND_K,
    BOUND_MIN_SITES,
    CHECKSUM_CACHE,
    DEPTH_HIGH,
    DEPTH_LOW,
    EPS_INIT,
    EVIDENCE_CHUNK,
    FLAG_FRAC,
    FLAG_MIN,
    FLANK_LOW_COPIES,
    FLANK_POINTS,
    LOWQ_FRAC,
    PH_MAX_DIV_NOQUAL,
    PROGRESS_CHUNKS,
    QUAL_FRAC_MIN,
    READS_STEM,
    RES_WIN,
)
from synltr.base.log import UsageError
from synltr.base.runtime import pool_map
from synltr.io.seqio import Fasta
from synltr.routes.reads.model import (
    CLASSES,
    ReadModel,
    allele_posteriors,
    call,
    depth_and_spans,
    estimate_eps,
    genotype_prior,
    genotypes,
    pav_call,
    read_profile,
    sample_points,
    site_calibration,
)
from synltr.routes.reads.phase import haplotypes
from synltr.routes.reads.reads import Site, Target, evidence

log = logging.getLogger(__name__)

COUNT_KEYS = ["F", "F_span", "F_other", "F_left", "F_right", "S", "S_split", "E", "E_split", "P", "O", "U", "lowq"]
PHASE_KEYS = ["hap_F", "hap_S", "hap_E", "phase_sites"]
VCF_IDX = {"F": "0", "E": "1", "S": "2", "N": ".", ".": "."}  # allele -> VCF allele index
TSV_COLUMNS = (["element_id", "chrom", "start", "end", "strand", "family", "clade", "length", "ltr5_len", "ltr3_len",
                "tsd", "n_inner", "n_outer", "sample"] + [f"n_{k}" for k in COUNT_KEYS] + PHASE_KEYS
               + ["other_delta_median", "r_F", "r_S", "r_E", "exp_reads", "flank_depth", "flank_ratio", "genotype",
                  "gq", "vcf_gt", "p_F", "p_E", "p_S", "pav", "flags"])


@dataclass(frozen=True, slots=True)
class Call:
    """One read set's call at one element, as the reference VCF holds it."""

    gt: str  # VCF genotype: 0 full, 1 empty site, 2 solo LTR, . null or no call (e.g. 0/1)
    gq: int
    ad: tuple[int, int, int]  # informative reads: full, empty, solo
    flags: tuple[str, ...]  # per-sample flags; none = PASS


@dataclass(frozen=True)
class SampleCalls:
    """A read set's calls: element index (into the route's elements) -> Call, for the elements on its BAM."""

    label: str
    ploidy: int
    calls: dict[int, Call]


# ---------------------------------------------------------------- evidence (forked workers)

_job: dict = {}  # the read set being scanned (BAM, reference, MAPQ, phasing ploidy): set before the workers fork
_handles: dict = {}  # pid -> (BAM, FASTA) that process opened: pysam handles are never shared across a fork


def _open() -> tuple:
    """This process's BAM and FASTA handles for the read set in `_job`, opened on first use."""
    h = _handles.get(os.getpid())
    if h is None:
        bam = pysam.AlignmentFile(str(_job["bam"]), reference_filename=str(_job["reference"]))
        h = _handles[os.getpid()] = (bam, Fasta(_job["reference"], cache_dir=_job["fai_cache"]))
    return h


def _close() -> None:
    """Close the handles this process opened."""
    h = _handles.pop(os.getpid(), None)
    if h is not None:
        h[0].close()
        h[1].close()


def _point_depth(bam, chrom: str, p: int, q: int) -> int:
    return sum(1 for x in bam.fetch(chrom, p, p + 1)
               if not (x.is_unmapped or x.is_secondary or x.is_supplementary or x.is_duplicate or x.is_qcfail)
               and x.mapping_quality >= q)


def _scan_chunk(chunk: list[tuple[int, Target]]) -> list[tuple[int, dict]]:
    """Evidence, read phasing and flank depth of each (index, element) of `chunk`."""
    bam, fa = _open()
    q = _job["min_mapq"]
    out = []
    for i, el in chunk:
        rd = []
        c = evidence(bam, fa, el, el.tsd_len, q, reads_out=rd)
        h, nsites = haplotypes(rd, fa, el.chrom, el.start, el.end, _job["ploidy"])
        c.update({f"hap_{a}": h[a] for a in CLASSES}, phase_sites=nsites)
        L = fa.length(el.chrom)
        pts = [el.start - d for d in FLANK_POINTS] + [el.end + d for d in FLANK_POINTS]
        fd = [_point_depth(bam, el.chrom, p, q) for p in pts if 0 <= p < L]
        c["flank_depth"] = median(fd) if fd else -1
        out.append((i, c))
    return out


def _elements_digest(targets: list[Target]) -> str:
    """sha256 of every element's ID, sequence, extent, LTRs and TSD, in order: the evidence is indexed by element."""
    h = hashlib.sha256()
    for e in targets:
        h.update(f"{e.id}\t{e.chrom}\t{e.start}\t{e.end}\t{e.lltr}\t{e.rltr}\t{e.tsd}\n".encode())
    return h.hexdigest()


def evidence_path(outdir: Path, label: str) -> Path:
    """The evidence checkpoint of read set `label`: <outdir>/reads.<label>.evidence.json.gz."""
    return Path(outdir) / f"{READS_STEM}.{label}.evidence.json.gz"


def scan(label: str, bam_path: Path, reference: Path, targets: list[Target], threads: int, *, ploidy: int,
         min_mapq: int, outdir: Path) -> tuple[dict[int, dict], list, list]:
    """The evidence of read set `label` at every element on its BAM's sequences: (counts by element index, depths at
    the sampled points, read spans), reused from its checkpoint when that matches, else gathered with `threads`
    workers and checkpointed."""
    outdir = Path(outdir)
    cache = evidence_path(outdir, label)
    done = marker(cache)
    inputs = {"bam": bam_key(bam_path), "reference": file_checksum(reference, outdir / "work" / CHECKSUM_CACHE)}
    params = {"min_mapq": min_mapq, "ploidy": ploidy, "elements": _elements_digest(targets)}
    key = stage_key("reads_evidence", inputs, params, code_fingerprint("routes.reads.genotype"))
    if is_done(done, key):
        try:
            with gzip.open(cache, "rt") as fh:
                d = json.load(fh)
            ev = {int(k): v for k, v in d["evidence"].items()}
            log.info("%s: evidence done, skip (%s)", label, cache)
            return ev, d["depths"], d["spans"]
        except (OSError, EOFError, ValueError, KeyError) as exc:
            log.warning("%s: evidence checkpoint %s is unreadable (%s); rescanning", label, cache, exc)
    clear_done(done)
    ev, depths, spans = _collect(label, bam_path, reference, targets, threads, ploidy=ploidy, min_mapq=min_mapq,
                                 outdir=outdir)
    data = json.dumps({"depths": depths, "spans": spans, "evidence": {str(i): c for i, c in ev.items()}})
    # through a temporary file: a killed run leaves no truncated checkpoint
    with atomic_path(cache) as tmp, open(tmp, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0,
                                                                          filename="") as gz:
        gz.write(data.encode())
    mark_done(done, key, [cache])
    return ev, depths, spans


def _collect(label: str, bam_path: Path, reference: Path, targets: list[Target], threads: int, *, ploidy: int,
             min_mapq: int, outdir: Path) -> tuple[dict[int, dict], list, list]:
    """Sample the read set's depth and read spans, decide whether its reads can be phased, then gather the evidence
    of every element on its sequences (in element order)."""
    with pysam.AlignmentFile(str(bam_path), reference_filename=str(reference)) as bam:
        held = set(bam.references)
        use = [i for i, el in enumerate(targets) if el.chrom in held]
        # sorted: the sampled points do not depend on the order of the BAM's header
        chroms = {c: bam.get_reference_length(c) for c in sorted({targets[i].chrom for i in use})}
        pts = sample_points(chroms)
        depths, spans = depth_and_spans(bam, pts, min_mapq)
        try:
            ReadModel(depths, spans, ploidy)  # fail before the evidence scan, not after it
        except ValueError as err:
            raise UsageError(f"{label} ({bam_path}): {err}") from None
        qfrac, div = read_profile(bam, pts, min_mapq)
    phase_ploidy = ploidy
    if qfrac < QUAL_FRAC_MIN and not (div <= PH_MAX_DIV_NOQUAL):
        phase_ploidy = 1  # haplotypes() is off for ploidy 1
        log.warning("%s: reads lack base qualities and diverge %.3f from the reference: read phasing off "
                    "(sequencing errors would look like haplotype differences)", label, div)
    log.debug("%s: reads with base qualities %.2f, median divergence %.4f", label, qfrac, div)
    work = [(i, targets[i]) for i in use]
    chunks = [work[j:j + EVIDENCE_CHUNK] for j in range(0, len(work), EVIDENCE_CHUNK)]
    _job.update(bam=bam_path, reference=reference, min_mapq=min_mapq, ploidy=phase_ploidy,
                fai_cache=Path(outdir) / "work" / "fai")
    t0 = time.perf_counter()
    finished = [0, 0]  # chunks and elements finished so far

    def progress(k: int) -> None:
        finished[0] += 1
        finished[1] += len(chunks[k])
        if finished[0] % PROGRESS_CHUNKS == 0 or finished[0] == len(chunks):
            log.debug("%s: %d/%d elements (%.0f s)", label, finished[1], len(work), time.perf_counter() - t0)

    try:
        results = pool_map(_scan_chunk, chunks, threads, warm=False, progress=progress)
    finally:
        _close()
        _job.clear()
    return {i: c for res in results for i, c in res}, depths, spans


# ---------------------------------------------------------------- the sample's model and calls

def boundary_tolerance(ev: dict[int, dict], sites: dict[int, Site]) -> tuple[float, int]:
    """Annotation boundary error, learned from the sample: at elements whose ends no TSD confirms, the median offset
    of empty-site and other-length spanning reads from the expected empty length (within RES_WIN bp; >= 2 reads
    per element). Returns (extra tolerance = BOUND_K x median |offset| over elements, elements used); 0 with fewer
    than BOUND_MIN_SITES elements (too little evidence: ends taken as exact)."""
    meds = []
    for i, c in ev.items():
        st = sites[i]
        if st.tsd_found:
            continue
        exp = st.delta["E"]
        r = [d - exp for d in c.get("other_delta", [])] + [d - exp for lab, d in c.get("es_delta", []) if lab == "E"]
        r = [x for x in r if abs(x) <= RES_WIN]
        if len(r) >= 2:
            meds.append(abs(median(r)))
    if len(meds) < BOUND_MIN_SITES:
        return 0.0, len(meds)
    return BOUND_K * float(median(meds)), len(meds)


def relabel_other(ev: dict[int, dict], sites: dict[int, Site], extra: float) -> int:
    """Other-length spanning reads at elements without a TSD that fit the empty or solo length within the widened
    tolerance are relabeled (counts updated in place). Returns the number of reads relabeled."""
    n = 0
    if extra <= 0:
        return 0
    for i, c in ev.items():
        st = sites[i]
        if st.tsd_found or not c.get("other_delta"):
            continue
        keep = []
        for d in c["other_delta"]:
            lab = st.classify_delta(d, extra)
            if lab in ("E", "S"):
                c[lab] += 1
                c["O"] -= 1
                n += 1
            else:
                keep.append(d)
        c["other_delta"] = keep
    return n


def flags_for(c: dict, gt: str, mu: float, flank_ratio: float, ploidy: int) -> list[str]:
    f = []
    T = sum(c[k] for k in CLASSES)
    if mu > 0 and T < DEPTH_LOW * mu:
        f.append("low_depth")
    if mu > 0 and T > DEPTH_HIGH * mu:
        f.append("high_depth")
    # other allele lengths are only visible in reads that span the site: compare with those
    spanning = c["F_span"] + c["S"] + c["E"] + c["O"]
    if c["O"] >= FLAG_MIN and c["O"] >= FLAG_FRAC * spanning:
        f.append("other_allele")
    if c.get("F_other", 0) >= FLAG_MIN and c.get("F_other", 0) >= FLAG_FRAC * spanning:
        f.append("full_other_length")  # a full allele whose internal region differs in length (indel, nesting)
    if c["U"] >= FLAG_MIN and c["U"] >= FLAG_FRAC * (T + c["U"]):
        f.append("unexplained_clips")
    if 0 <= flank_ratio < (ploidy - FLANK_LOW_COPIES) / ploidy:
        f.append("flank_depth_low")
    lq = c.get("lowq", 0)
    if lq >= FLAG_MIN and lq >= LOWQ_FRAC * (lq + c.get("highq", 0)):
        f.append("low_mapq")  # many reads here map ambiguously: depth-based dosage is unreliable
    if T == 0:
        f.append("no_reads")
    return f


def genotype_sample(label: str, bam_path: Path, reference: Path, targets: list[Target], threads: int, *,
                    ploidy: int, min_mapq: int, outdir: Path) -> tuple[SampleCalls, list[str]]:
    """Genotype read set `label` (its BAM aligned to `reference`) at every element of `targets`: its calls and its
    genotypes.tsv lines."""
    outdir = Path(outdir)
    ev, depths, spans = scan(label, bam_path, reference, targets, threads, ploidy=ploidy, min_mapq=min_mapq,
                             outdir=outdir)
    model = ReadModel(depths, spans, ploidy)
    sites = {i: Site(targets[i], targets[i].tsd_len) for i in ev}
    tol, nb = boundary_tolerance(ev, sites)
    moved = relabel_other(ev, sites, tol)
    log.debug("%s: boundary error of TSD-less element ends: tolerance +%.0f bp (%d elements); %d other-length reads "
              "relabeled empty/solo", label, tol, nb, moved)
    rates = {i: model.rates(sites[i]) for i in ev}
    log.debug("%s: genome-wide depth: per haplotype %.1f, copy weights %s, NB size %.0f", label, model.d_hap,
              ",".join(f"{w:.2f}" for w in model.copy_weights), model.k)
    calls = {i: call(ev[i], rates[i], model, EPS_INIT) for i in ev}
    eps = estimate_eps((ev[i], calls[i][0], calls[i][1]) for i in ev) or EPS_INIT
    cal = site_calibration(((ev[i], rates[i], calls[i][0], calls[i][1]) for i in ev), model)
    if cal is not None:
        model.d_hap, model.k = cal[0], cal[1]
        model.depth = model.d_hap * model.ploidy
        log.debug("%s: depth at element sites from %d sites with an empty/solo allele: per haplotype %.1f, NB size "
                  "%.0f", label, cal[2], cal[0], cal[1])
    else:
        log.debug("%s: too few sites with an empty/solo allele to calibrate depth at element sites (genome-wide "
                  "depth kept)", label)
    calls = {i: call(ev[i], rates[i], model, eps) for i in ev}
    gts = genotypes(model.ploidy)
    prior = genotype_prior([calls[i][2] for i in ev if not calls[i][0].startswith(".")], gts)
    if prior is not None:
        calls = {i: call(ev[i], rates[i], model, eps, prior) for i in ev}
    log.debug("%s: depth %.1f (per haplotype %.1f), NB size %.0f, mean read span %.0f bp, error %.4f", label,
              model.depth, model.d_hap, model.k, model.mean_span, eps)
    if prior is None:
        log.debug("%s: too few called sites to learn a genotype prior (flat prior)", label)
    else:
        log.debug("%s: genotype prior (EM over sites): %s", label, ", ".join(
            f"{g} {math.exp(p):.3f}" for g, p in sorted(prior.items(), key=lambda kv: -kv[1]) if math.exp(p) >= 0.01))
    out, lines = {}, []
    for i in sorted(ev):
        row, out[i] = _row(targets[i], label, ev[i], rates[i], calls[i], model)
        lines.append("\t".join(map(str, row)) + "\n")
    counts: dict[str, int] = {}
    for i in ev:
        counts[calls[i][0]] = counts.get(calls[i][0], 0) + 1
    log.info("%s: %d elements: %s", label, len(ev), ", ".join(f"{g} {n}" for g, n in
                                                               sorted(counts.items(), key=lambda x: -x[1])))
    return SampleCalls(label, ploidy, out), lines


def write_genotypes(path: Path, parts: list[list[str]]) -> None:
    """ltrkin's genotypes.tsv at `path`: the header, then each read set's lines (`parts`), in order."""
    with atomic_path(Path(path)) as tmp, open(tmp, "w") as fh:
        fh.write("\t".join(TSV_COLUMNS) + "\n")
        for lines in parts:
            fh.writelines(lines)


def _row(el: Target, label: str, c: dict, r: dict, called: tuple, model: ReadModel) -> tuple[list, Call]:
    """The genotypes.tsv row of element `el` in read set `label`, and its Call."""
    gt, gq, lls, mu = called
    fr = c["flank_depth"] / model.depth if c["flank_depth"] >= 0 else -1
    od = c["other_delta"]
    flags = flags_for(c, gt, mu, fr, model.ploidy)
    vcf_gt = "/".join(VCF_IDX[a] for a in gt.split("/"))
    row = [el.id, el.chrom, el.start + 1, el.end, el.strand, el.family, el.clade, el.length,
           el.lltr[1] - el.lltr[0], el.rltr[1] - el.rltr[0], el.tsd or ".", el.n_inner, el.n_outer, label]
    row += [c.get(k, 0) for k in COUNT_KEYS] + [c.get(k, 0) for k in PHASE_KEYS]
    pa = allele_posteriors(lls)
    row += [int(median(od)) if od else ".", f"{r['F']:.3f}", f"{r['S']:.3f}", f"{r['E']:.3f}",
            f"{mu:.1f}", c["flank_depth"], f"{fr:.2f}", gt, gq, vcf_gt,
            f"{pa['F']:.4f}", f"{pa['E']:.4f}", f"{pa['S']:.4f}", pav_call(pa), ",".join(flags) or "PASS"]
    return row, Call(vcf_gt, gq, (c["F"], c["E"], c["S"]), tuple(flags))
