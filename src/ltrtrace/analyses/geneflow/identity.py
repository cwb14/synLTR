"""Optional sequence-identity layer: divergence between orthologous LTR copies of shared
insertions, the rate-free split estimate and the shared insertions younger than the split. It lies outside the
PAV-only core: it runs only with --identity, its rows are layer "identity", its findings have their own
section, and its per-locus counts go to identity.tsv.

Alignment. The pair files keep no LTR site counts, so copies are realigned with the route's own window aligner,
`routes.assembly.lastz.align_window` (LASTZ, asm20 scoring): the same aligner, one more pass. `ltr_counts` compares
the bases of each gap-free aligned segment over its overlap with x's LTR intervals, so sites and substitutions are
exact counts; a column with N on either side is skipped, and so is one with any other ambiguity code (only A, C, G
and T are compared). `self_counts` aligns a copy's 5' LTR to its 3' LTR the same way.

Loci read, per trio (R, S; D), each state from D's or R's own calls:
  pair loci   D's members where D's own call into Y (R or S) is P: only P copies are compared, intact copies
              holding both LTRs
  counted     the D_pav counted rows, `relations.discover(cat, R, S, D, dist, "E")` (all tiers, every outgroup E)
  RD class    R's RD painting class (`paint.classes`: R's members where R's own calls give D present, S E)

identity.tsv (IDENTITY_COLUMNS), per locus in catalog order:
  self  one per annotated copy of R, S or D at those loci whose two LTRs align over >= MIN_SEG sites: genome
        the copy's, partner ".", div = diffs / sites, u95 = `seq.divergence.garwood_u95(diffs, sites)`; young is
        1 or 0 on D's copy at the counted R+D and S+D loci once s exists, "." elsewhere
  pair  one per pair locus and partner Y with two P copies: genome D, partner Y; div "." below IDENT_MIN_SITES
        LTR sites (such pairs never enter a fit); u95 and young "."
Loci with neither get no row. Self rows are what painting dates tracts with (`paint`'s self_counts: from R's
copy, else from D's copy of the same insertion, which exists only where D is a member).

Pair alignment. x = D's element +- WINDOW_PAD and y = Y's site from D's own call (pseq, pstart, pend) +-
WINDOW_PAD, both fetched on their genome's + strand: the shape of the engine's own windows. D's flanks anchor
LASTZ's chain on the site's diagonal, so an LTR cannot slide onto paralogous sequence inside the partner's
site, as it can with a bare element. Substitutions are still counted in D's LTR columns only: ltr_left and
ltr_right in x coordinates, shifted by the left pad actually fetched (clipped at a sequence start). y is
reverse-complemented where D's own call reverses the window: its b_strand (`Catalog.pbstrand`), "+" or "-" on
every sited call and known whatever D's own strand, "." included; the orientation is never read from the
alignment. A pair whose call has no usable b_strand is skipped and counted, with a WARNING.

Rate-free split (`nll`, `fit`, `fit_mix`, `pred_cdf`). S-D ortholog pairs
coalesced before the split: k ~ Poisson(n (s + X)), X ~ Exp(mean m), with s = 2 mu T_split and X the ancestral
coalescence, so log P(k) = s/m - log(m n) - (k+1) log(1 + 1/(n m)) + log Q(k+1, n (1 + 1/(n m)) s), Q the
regularized upper incomplete gamma; Nelder-Mead on (log s, log m), on the S-D pairs with >= IDENT_MIN_SITES
LTR sites. Q is never clamped: a clamp (at 1e-300, say) would make the
likelihood grow without bound as m -> 0 (s/m grows while log Q stays at the clamp), so a fit to data without
overdispersion runs to m ~ 1e-310 and a likelihood of -inf; where Q < Q_MIN the exact finite sum Q(k+1, x) =
e^-x sum_{j<=k} x^j / j! is used instead, in log space with s/m cancelled, and the likelihood tends to the
Poisson one. Where Q >= Q_MIN, log Q is taken from gammaincc directly. The se of s and m: a delete-one
block jackknife over D's blocks (`stats.jackknife`, blocks of `size` bp on D's element), each refit starting
from the full fit; se is nan, the estimate kept, when the pairs lie in fewer than JK_MIN_BLOCKS blocks or the
pseudovalues do not vary (trio.py's testability rule). No z: s and m are positive by construction.

Guard (`fit_mix`, 18 Nelder-Mead starts): isolation + post-split gene flow (unshifted, mean <= MG_MAX) + an
unshifted upper tail (old paralogs, mis-assigned pairs). Post-split S-D gene flow pulls the isolation s down,
so the mixture's s_mix stands beside s, with its weights, means and AIC (2 nll + 12 against 2 nll + 4). s_mix
is the shift of the isolation component alone: it corrects for gene flow only through the gene-flow weight,
and the tail is an upper tail only when its mean exceeds the isolation component's (s_mix + m_mix); the
findings line states both, since a better AIC alone does not say which component earned it. A component whose
weight sits at its lower bound (weight x pairs < MIX_MIN_PAIRS: the weights are a softmax of free logits, so 0
is approached, never reached) rests on no data: its parameters are written NA. Fit check (CDF_POINTS: 0,
0.001, 0.002, 0.003, 0.005): the observed fraction of S-D pairs with diffs / sites <= x against `pred_cdf`
under the isolation fit.

Younger than the split, on the D_pav counted rows: an R+D (R present, S E) or S+D locus whose D copy has a
self row is young when its u95 < s; n_RD_young, n_SD_young and D_young = (n_RD_young - n_SD_young) /
(n_RD_young + n_SD_young), jackknifed on D's blocks with s held at its estimate (trio.py's testability rule).
sites_needed = u95(0, 1) / s: the identical LTR sites a single insertion needs to qualify. Without gene flow,
an insertion younger than the split cannot be shared by D and one sister; LTR-LTR gene conversion makes some
old insertions look young, in both classes alike. Below the split: the fractions of R-D and S-D pairs (>=
IDENT_MIN_SITES sites) with diffs / sites < s.

Pair-divergence layer (`write_pair_divergences`): the pair rows, every trio's, as the catalog's sparse layer of
LTR sites and differences between orthologous copies, pairs/<D>__<Y>.div.tsv.gz per donor D and partner Y
(calls.records.write_divergences; DIV_COLUMNS: D's element, Y's member element at the locus or none, the sites
and differences, d = differences / sites, none without a site; meta engine "assembly", d_kind "p", the LASTZ
p-distance over aligned LTR sites). A pair aligned in several trios is written once. Each run replaces the
assembly layer of the run before; another route's divergence files stay.

Workers: base.runtime.pool_map over chunks of loci (IDENTITY_CHUNKS_PER_WORKER per worker), every chunk in a
forked worker, all at once (warm=False: no chunk waits for this process), each worker opening its own Fasta (cache
work/fai: handles are never shared across processes) and killing its running lastz child if it is terminated
(lastz.kill_lastz); results come back in locus order. With one thread, all in this process.

Checkpoint: the alignment of each trio, about 20 min on maize, is the one checkpointed analysis (base.checkpoint):
<outdir>/work/geneflow/identity.<R>__<S>__<D>.json.gz and its done marker, keyed by the alignment jobs (every
copy's and window's coordinates), the checksums of the FASTAs they read, LASTZ's version and the code
fingerprint of this module's import closure. A rerun with the same key reuses the counts; everything else reruns.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import logging
import math
import os
import time
from collections import Counter
from pathlib import Path

import numpy as np

from ltrtrace.analyses.geneflow.catalog import STATE_CODE, AnalysisError, Catalog, patterns, site_blocks
from ltrtrace.analyses.geneflow.paint import RD, classes
from ltrtrace.analyses.geneflow.relations import discover
from ltrtrace.analyses.geneflow.route import ROUTE_PROPS
from ltrtrace.analyses.geneflow.trio import NOT_TESTABLE, Trio
from ltrtrace.base.checkpoint import (
    atomic_path,
    clear_done,
    code_fingerprint,
    done_meta,
    file_checksum,
    is_done,
    mark_done,
    marker,
    stage_key,
)
from ltrtrace.base.constants import (
    CDF_POINTS,
    CDF_SIMS,
    CHECKSUM_CACHE,
    IDENT_MIN_SITES,
    IDENTITY_CHUNKS_PER_WORKER,
    MG_MAX,
    MIN_SEG,
    MIX_FIT_MAXITER,
    MIX_FIT_TOL,
    MIX_MIN_PAIRS,
    MIX_START_FLOW,
    MIX_START_SHIFT,
    MIX_START_TAIL,
    MIX_START_WEIGHTS,
    Q_MIN,
    SPLIT_FIT_MAXITER,
    SPLIT_FIT_START,
    SPLIT_FIT_TOL,
    TAIL_EPS,
    WINDOW_PAD,
)
from ltrtrace.base.runtime import PoolError, pool_map
from ltrtrace.calls.records import tsv_field, write_divergences
from ltrtrace.io.elements import Element, read_elements
from ltrtrace.io.samples import PAIR_SEP
from ltrtrace.io.seqio import Fasta, revcomp
from ltrtrace.routes.assembly.lastz import AlnMap, align_window, kill_lastz, lastz_version, sweep_lastz_tmp
from ltrtrace.seq.divergence import garwood_u95 as u95
from ltrtrace.seq.intervals import union
from ltrtrace.stats import report
from ltrtrace.stats.report import Findings, StatRow
from ltrtrace.stats.resample import d_ratio, jackknife

log = logging.getLogger(__name__)

IDENTITY_COLUMNS = ("kind", "lid", "trio", "genome", "partner", "sites", "diffs", "div", "u95", "young")
SECTION = "Sequence-identity layer (outside the PAV-only core)"  # its findings section (stats.report.SECTIONS)
LAYER = "identity"
ALIGNER = "lastz"  # the window aligner of the assembly route, named in the findings
CHECKPOINTS = Path("work") / "geneflow"  # under the run directory: the alignment checkpoints
DIV_META = {"engine": "assembly", "d_kind": "p"}  # the pair-divergence layer's meta: LASTZ p-distance
_ACGT = np.zeros(256, dtype=bool)
_ACGT[np.frombuffer(b"ACGT", dtype=np.uint8)] = True  # the bases compared; N and other codes are skipped


# ---------------------------------------------------------------- LTR counts


def ltr_counts(x: str, y: str, ltrs: list[tuple[int, int]]) -> tuple[int, int]:
    """(sites, substitutions) in x's LTR columns: x aligned to y (already on x's strand) with the route's
    aligner, align_window(x, y); each gap-free segment compared base by base over its overlap with the intervals
    `ltrs` ([a, b), x coordinates; overlapping intervals count once). A column with any base other than A, C, G, T
    (N included) on either side is skipped."""
    return _count(align_window(x, y), x, y, ltrs)


def self_counts(ltr5: str, ltr3: str) -> tuple[int, int]:
    """(sites, substitutions) between a copy's two LTRs: the 5' LTR aligned to the 3' LTR as `ltr_counts`
    does, every column of the 5' LTR an LTR column."""
    return ltr_counts(ltr5, ltr3, [(0, len(ltr5))])


def _count(aln: AlnMap, x: str, y: str, ltrs: list[tuple[int, int]]) -> tuple[int, int]:
    """`ltr_counts` on an alignment of x to y."""
    xb = np.frombuffer(x.upper().encode("ascii"), dtype=np.uint8)
    yb = np.frombuffer(y.upper().encode("ascii"), dtype=np.uint8)
    sites = diffs = 0
    for a0, a1 in union((max(0, int(a)), min(len(xb), int(b))) for a, b in ltrs):  # clipped to x
        lo = np.maximum(aln.seg_a, a0)
        hi = np.minimum(aln.seg_a + aln.seg_len, a1)
        for j in np.flatnonzero(hi > lo).tolist():
            p, q = int(lo[j]), int(hi[j])
            b = int(aln.seg_b[j]) + p - int(aln.seg_a[j])
            xs, ys = xb[p:q], yb[b:b + q - p]
            ok = _ACGT[xs] & _ACGT[ys]
            sites += int(np.count_nonzero(ok))
            diffs += int(np.count_nonzero(ok & (xs != ys)))
    return sites, diffs


# ---------------------------------------------------------------- the rate-free split


def split_nll(params: np.ndarray, k: np.ndarray, n: np.ndarray) -> float:
    """Negative log-likelihood of the isolation model, params = (log s, log m), k substitutions
    over n LTR sites per pair; exact where Q underflows (module docstring), +inf where it is not finite."""
    with np.errstate(all="ignore"):
        s, m = np.exp(np.asarray(params, dtype=np.float64))
        v = -float(_lp(k, n, s, m).sum())
    return v if math.isfinite(v) else math.inf


def fit_split(k: np.ndarray, n: np.ndarray, start: tuple[float, float] = SPLIT_FIT_START) -> tuple[float, float]:
    """(s, m), the isolation model's maximum-likelihood fit (Nelder-Mead on the logs from
    `start`, xatol = fatol = 1e-6, at most 4000 iterations)."""
    from scipy.optimize import minimize
    k, n = np.asarray(k, dtype=np.float64), np.asarray(n, dtype=np.float64)
    with np.errstate(all="ignore"):  # a +inf vertex makes Nelder-Mead's own convergence check warn
        r = minimize(split_nll, np.log(start), args=(k, n), method="Nelder-Mead",
                     options={"xatol": SPLIT_FIT_TOL, "fatol": SPLIT_FIT_TOL, "maxiter": SPLIT_FIT_MAXITER})
    s, m = np.exp(r.x)
    return float(s), float(m)


def fit_mix(k: np.ndarray, n: np.ndarray, s0: float, m0: float
            ) -> tuple[np.ndarray, float, float, float, float, float]:
    """The mixture fit: isolation (shifted exponential) + post-split gene flow (unshifted, mean <= MG_MAX) +
    an unshifted upper tail (old paralogs, mis-assigned pairs), Nelder-Mead from its 18 starts around the
    isolation fit (s0, m0). Returns (weights isolation/gene flow/tail, gene-flow mean, s, m, tail mean, nll)."""
    from scipy.optimize import minimize
    k, n = np.asarray(k, dtype=np.float64), np.asarray(n, dtype=np.float64)
    best = None
    w_iso, w_tail = MIX_START_WEIGHTS
    for fg in MIX_START_FLOW:
        for sm in MIX_START_SHIFT:
            for mh0 in MIX_START_TAIL:
                x0 = [np.log(fg / w_iso), np.log(w_tail / w_iso), 0.0, np.log(max(s0, 1e-5) * sm), np.log(m0),
                      np.log(mh0)]
                with np.errstate(all="ignore"):  # as in fit_split
                    r = minimize(_mix_nll, x0, args=(k, n), method="Nelder-Mead",
                                 options={"xatol": MIX_FIT_TOL, "fatol": MIX_FIT_TOL, "maxiter": MIX_FIT_MAXITER})
                if best is None or r.fun < best.fun:
                    best = r
    w, mg, s, m, mh = _mix_par(best.x)
    return w, float(mg), float(s), float(m), float(mh), float(best.fun)


def pred_cdf(x: float, n: np.ndarray, s: float, m: float, nsim: int = CDF_SIMS, seed: int = 0) -> float:
    """The fraction of pairs with divergence k / n <= x predicted by the isolation fit (s,
    m), simulated `nsim` times at the observed n (generator seeded with `seed`)."""
    rng = np.random.default_rng(seed)
    nn = np.tile(np.asarray(n, dtype=np.float64), nsim)
    kk = rng.poisson(nn * (s + rng.exponential(m, len(nn))))
    return float(np.mean(kk / nn <= x))


def _lp(k: np.ndarray, n: np.ndarray, s: float, m: float) -> np.ndarray:
    """Per-pair log-likelihood (s = 0 is the unshifted exponential). Where Q(k+1, x) < Q_MIN,
    log Q = -x + log sum_{j<=k} x^j / j! with x = n s + s/m, so s/m cancels exactly."""
    from scipy.special import gammaincc
    c = 1.0 + 1.0 / (n * m)
    if s <= 0:
        return -np.log(m * n) - (k + 1) * np.log(c)
    x = c * n * s
    q = gammaincc(k + 1, x)
    tiny = q < Q_MIN
    lp = s / m - np.log(m * n) - (k + 1) * np.log(c) + np.log(np.where(tiny, 1.0, q))
    if tiny.any():
        kt, xt = k[tiny], x[tiny]
        lp[tiny] = -np.log(m * n[tiny]) - (kt + 1) * np.log(c[tiny]) - n[tiny] * s + _log_poisson_sum(kt, xt)
    return lp


def _log_poisson_sum(k: np.ndarray, x: np.ndarray) -> np.ndarray:
    """log sum_{j=0..k} x^j / j! (= x + log Q(k + 1, x)) for integer k >= 0 and x > 0, in log space from the
    j = k term down, so no term overflows."""
    from scipy.special import gammaln
    head = k * np.log(x) - gammaln(k + 1)
    lt = np.zeros_like(x)   # log(term j = k - i - 1 / term j = k)
    acc = np.zeros_like(x)  # log(sum of the terms so far / term j = k)
    for i in range(int(k.max()) if k.size else 0):
        with np.errstate(divide="ignore"):
            lt = lt + np.log(np.maximum(k - i, 0.0) / x)
        acc = np.logaddexp(acc, lt)
        if np.all((lt < acc - TAIL_EPS) & (k - i <= x)):
            break
    return head + acc


def _mix_par(par: np.ndarray) -> tuple[np.ndarray, float, float, float, float]:
    """The mixture's parameters: (weights isolation/gene flow/tail, gene-flow mean, s, m, tail mean)."""
    with np.errstate(over="ignore"):
        w = np.exp(np.r_[0.0, par[0], par[1]])
        w /= w.sum()
        mg = MG_MAX / (1 + np.exp(-par[2]))
        s, m, mh = np.exp(par[3:])
    return w, mg, s, m, mh


def _mix_nll(par: np.ndarray, k: np.ndarray, n: np.ndarray) -> float:
    """Negative log-likelihood of the mixture; +inf where it is not finite."""
    with np.errstate(all="ignore"):
        w, mg, s, m, mh = _mix_par(par)
        lp = np.stack([np.log(w[0]) + _lp(k, n, s, m), np.log(w[1]) + _lp(k, n, 0.0, mg),
                       np.log(w[2]) + _lp(k, n, 0.0, mh)])
        v = -float(np.logaddexp.reduce(lp, axis=0).sum())
    return v if math.isfinite(v) else math.inf


# ---------------------------------------------------------------- workers

_FASTA_SPECS: dict[str, tuple[Path, Path]] = {}  # genome -> (FASTA, index cache), set before the alignments
_FASTAS: dict = {"pid": None, "open": {}}         # the FASTAs this process opened (pysam handles are per process)


def _fetch(label: str, seqid: str, a: int, b: int) -> str:
    """[a, b) of a genome's sequence (0-based, clipped), from this process's own Fasta: a forked worker drops the
    handles it inherited, never reads through them, and opens its own on first use."""
    if _FASTAS["pid"] != os.getpid():
        _FASTAS.update(pid=os.getpid(), open={})
    fa = _FASTAS["open"].get(label)
    if fa is None:
        path, cache = _FASTA_SPECS[label]
        fa = _FASTAS["open"][label] = Fasta(path, cache_dir=cache)
    return fa.fetch(seqid, a, b)


def _close_fastas() -> None:
    if _FASTAS["pid"] == os.getpid():
        for fa in _FASTAS["open"].values():
            fa.close()
    _FASTAS.update(pid=None, open={})


def _locus(job: tuple) -> tuple[tuple, tuple]:
    """One locus: (self counts per copy, (sites, diffs, strand) per pair; strand 0 = skipped, the call has no
    usable b_strand)."""
    selfs, pair = job
    out_self = tuple(self_counts(_fetch(g, sq, a0, a1), _fetch(g, sq, b0, b1))
                     for g, sq, a0, a1, b0, b1 in selfs)
    out_pair = []
    if pair is not None:
        d, dsq, x0, x1, ltrs, ys = pair
        x = _fetch(d, dsq, x0, x1)
        for y_label, ysq, y0, y1, strand in ys:
            if strand == 0:
                out_pair.append((0, 0, 0))
                continue
            y = _fetch(y_label, ysq, y0, y1)
            yy = y if strand > 0 else revcomp(y)
            out_pair.append((*_count(align_window(x, yy), x, yy, ltrs), strand))
    return out_self, tuple(out_pair)


def _chunk(jobs: list[tuple]) -> list[tuple[tuple, tuple]]:
    return [_locus(j) for j in jobs]


def _align(cat: Catalog, jobs: list[tuple], threads: int) -> list[tuple[tuple, tuple]]:
    """Every locus job, in order: base.runtime.pool_map over chunks of them (module docstring)."""
    _close_fastas()
    _FASTA_SPECS.clear()
    _FASTA_SPECS.update({s.label: (Path(s.fasta), cat.outdir / "work" / "fai") for s in cat.samples})
    n_chunks = min(len(jobs), max(threads, 1) * IDENTITY_CHUNKS_PER_WORKER)
    size = -(-len(jobs) // n_chunks)
    chunks = [jobs[a:a + size] for a in range(0, len(jobs), size)]
    sweep_lastz_tmp()
    try:
        return [r for part in pool_map(_chunk, chunks, threads, on_exit=kill_lastz, warm=False) for r in part]
    finally:
        _close_fastas()
        sweep_lastz_tmp()  # the folders of workers that died


def _aligned(cat: Catalog, trio: Trio, jobs: list[tuple], threads: int) -> list[tuple[tuple, tuple]]:
    """Every locus job's counts, in order: from the trio's alignment checkpoint when its key holds (module
    docstring), else aligned (`_align`) and checkpointed. A failure, a missing FASTA or a dead worker included, is
    one RuntimeError naming the identity layer."""
    if not jobs:
        return []
    out = cat.outdir / CHECKPOINTS / f"identity.{PAIR_SEP.join((trio.R, trio.S, trio.D))}.json.gz"
    try:
        sums = cat.outdir / "work" / CHECKSUM_CACHE
        used = sorted({j[0] for js, _ in jobs for j in js} | {p[0] for _, p in jobs if p is not None}
                      | {y[0] for _, p in jobs if p is not None for y in p[5]})
        fastas = {s.label: s.fasta for s in cat.samples}
        inputs = {"jobs": hashlib.sha256(json.dumps(jobs, separators=(",", ":")).encode()).hexdigest(),
                  **{f"fasta:{g}": file_checksum(fastas[g], sums) for g in used}}
        key = stage_key("identity_alignment", inputs, {"trio": trio.scope, "aligner": lastz_version()},
                        code_fingerprint("analyses.geneflow.identity"))
        if is_done(marker(out), key):
            log.info("identity %s: alignment done, skip (%s)", trio.scope, out)
            with gzip.open(out, "rt", encoding="utf-8") as fh:
                return [(tuple(map(tuple, a)), tuple(map(tuple, b))) for a, b in json.load(fh)]
        clear_done(marker(out))
        n_pairs = sum(len(p[5]) for _, p in jobs if p is not None)
        log.info("identity %s: aligning %d ortholog pair(s) and %d copies' LTRs at %d loci (%d thread(s))",
                 trio.scope, n_pairs, sum(len(js) for js, _ in jobs), len(jobs), threads)
        t0 = time.monotonic()
        results = _align(cat, jobs, threads)
        dt = time.monotonic() - t0
        log.info("identity %s: aligned in %.1f s (%.1f pairs/s, %.1f loci/s)", trio.scope, dt,
                 n_pairs / dt if dt else 0, len(jobs) / dt if dt else 0)
        out.parent.mkdir(parents=True, exist_ok=True)
        with atomic_path(out) as tmp, gzip.open(tmp, "wt", encoding="utf-8") as fh:
            json.dump(results, fh, separators=(",", ":"))
        mark_done(marker(out), key, [out], {"trio": trio.scope})
        return results
    except PoolError as exc:
        raise RuntimeError(f"identity layer: {exc}") from exc
    except Exception as exc:  # a worker's own error, or a FASTA gone: one line naming the layer
        raise RuntimeError(f"identity layer: an alignment failed: {type(exc).__name__}: {exc}") from exc


# ---------------------------------------------------------------- the layer


def identity(cat: Catalog, trio: Trio, size: int, threads: int, dist: np.ndarray
             ) -> tuple[list[StatRow], Findings, list[tuple]]:
    """The sequence-identity layer of one trio (module docstring): its stats rows (layer identity), findings
    and identity.tsv rows (IDENTITY_COLUMNS). Skipped with one line when D has no annotation (no LTR bounds)."""
    f = Findings()
    R, S, D, scope = trio.R, trio.S, trio.D, trio.scope
    if not cat.annotated[D]:
        f.add(SECTION, f"Sequence-identity layer of {scope} skipped: {D} has no annotation, so its copies have "
                       "no LTR bounds to count substitutions in")
        return [], f, []
    gd = cat.g(D)
    c = discover(cat, R, S, D, dist, "E")  # D is annotated: it discovers, never None
    partners = [y for y in cat.labels if y in (R, S)]  # sheet order
    pmask = {y: cat.member[:, gd] & (cat.proj[D][:, cat.g(y)] == STATE_CODE["P"]) for y in partners}
    rd_class = np.zeros(len(cat.lid), dtype=bool)
    if cat.annotated[R]:
        # the layer realigns the assembly route's calls, so R's RD class is that route's painting alphabet's
        idx, cls = classes(cat, R, S, D, alphabet=ROUTE_PROPS["assembly"].alphabet[0])
        rd_class[idx[cls == RD]] = True
    reads = c.rows | rd_class
    for m in pmask.values():
        reads = reads | m
    loci = np.flatnonzero(reads)
    copies = [g for g in cat.labels if g in (R, S, D) and cat.annotated[g]]  # sheet order
    elements = {g: _elements(cat, g, loci) for g in copies}
    jobs = [_job(cat, int(i), copies, elements, pmask, D) for i in loci]
    results = _aligned(cat, trio, jobs, threads)
    tab = _Tables(loci, jobs, results, D)
    rows: list[StatRow] = []
    f.add(SECTION, _intro(trio, cat, loci, pmask, c.rows, rd_class, tab))
    if tab.n_skipped:
        f.warn(SECTION, f"{tab.n_skipped} pair(s) skipped: no usable b_strand in {D}'s own P call, so the "
                        "orientation of the partner's site is unknown (a sited call always has one: the pair "
                        "files are not the engine's as written)")
    sd = tab.fit_pairs(S)
    log.info("identity %s: fitting the divergence of %d %s-%s ortholog pair(s) with >= %d LTR sites, then "
             "%s's", scope, sd[0].size, S, D, IDENT_MIN_SITES, R)
    young = np.full(len(cat.lid), -1, dtype=np.int8)  # D's self row at a counted R+D or S+D locus: 1/0
    if not sd[0].size:
        f.add(SECTION, f"Split {S};{D} skipped: no {S}-{D} pair with >= {IDENT_MIN_SITES} LTR sites; so are the "
                       "younger-than-split counts and the fractions below the split")
    else:
        split_rows, s = _split(cat, trio, c, size, sd, f)
        rows += split_rows
        young_rows, young = _young(cat, trio, c, size, s, tab, f)
        rows += young_rows
        rows += _below(trio, s, tab.fit_pairs(R), sd, f)
    return rows, f, tab.out(cat, scope, young)


def write_identity(rows: list[tuple], path: Path) -> None:
    """identity.tsv: a header row of IDENTITY_COLUMNS, then one line per row (`tsv_field` values: None and
    nan as ".")."""
    with open(path, "w", newline="") as fh:
        fh.write("\t".join(IDENTITY_COLUMNS) + "\n")
        for r in rows:
            if len(r) != len(IDENTITY_COLUMNS):
                raise ValueError(f"identity row {r!r}: {len(r)} fields, IDENTITY_COLUMNS has "
                                 f"{len(IDENTITY_COLUMNS)}")
            fh.write("\t".join(tsv_field(v) for v in r) + "\n")


def assembly_divergences(outdir: Path) -> list[Path]:
    """The assembly route's pair-divergence files in <outdir>/pairs (by their done marker's meta), in name order."""
    return [p for p in sorted((Path(outdir) / "pairs").glob("*.div.tsv.gz"))
            if (done_meta(marker(p)) or {}).get("engine") == DIV_META["engine"]]


def write_pair_divergences(cat: Catalog, rows: list[tuple]) -> list[Path]:
    """The pair-divergence layer (module docstring) of the identity layer's rows (IDENTITY_COLUMNS, every trio's):
    pairs/<D>__<Y>.div.tsv.gz per donor D and partner Y, in sheet order of the pair, each row once in the order
    of `rows`; returns the paths. The assembly layer of the run before goes first; another route's divergence
    files stay."""
    pairs_dir = cat.outdir / "pairs"
    for path in assembly_divergences(cat.outdir):
        path.unlink()
        clear_done(marker(path))
    index = {lid: i for i, lid in enumerate(cat.lid.tolist())}
    layer: dict[tuple[str, str], dict[str, dict]] = {}
    for kind, lid, _, d, y, sites, diffs, *_ in rows:
        if kind == "pair":
            i = index[lid]
            layer.setdefault((d, y), {}).setdefault(lid, {
                "eid": cat.eid[i, cat.g(d)], "b_eid": cat.eid[i, cat.g(y)], "n_sites": sites, "n_diff": diffs,
                "d": diffs / sites if sites else None})
    out = []
    for (d, y), found in sorted(layer.items(), key=lambda kv: (cat.g(kv[0][0]), cat.g(kv[0][1]))):
        out.append(pairs_dir / f"{d}{PAIR_SEP}{y}.div.tsv.gz")
        write_divergences(out[-1], found.values(), DIV_META)
    return out


def _elements(cat: Catalog, g: str, loci: np.ndarray) -> dict[str, Element]:
    """Genome g's member elements at `loci` (eid -> Element), from its element table."""
    want = {e for e in cat.eid[loci, cat.g(g)].tolist() if e is not None}
    path = cat.outdir / "work" / "elements" / f"{g}.tsv"
    try:
        found = {e.eid: e for e in read_elements(path) if e.eid in want}
    except (OSError, ValueError, TypeError) as exc:
        raise AnalysisError(f"{path}: cannot read the element table ({exc})") from exc
    lost = want - found.keys()
    if lost:
        raise AnalysisError(f"{path}: {len(lost)} member element(s) of {g} in loci.tsv are not in it (e.g. "
                            f"{min(lost)}): the run's files do not match each other (rerun ltrtrace genotype)")
    return found


def _job(cat: Catalog, i: int, copies: list[str], elements: dict[str, dict[str, Element]],
         pmask: dict[str, np.ndarray], D: str) -> tuple:
    """Locus i's work, in plain values for the workers: each annotated copy's two LTRs (0-based [a, b) on its
    genome's + strand); D's element +- WINDOW_PAD with its LTRs in those x coordinates (shifted by the left pad
    actually fetched: it is clipped at a sequence start); and each partner's site +- WINDOW_PAD with the
    orientation of D's own call into it, b_strand (+1, -1; 0 = none, skipped)."""
    selfs = []
    for g in copies:
        e = cat.eid[i, cat.g(g)]
        if e is not None:
            el = elements[g][e]
            selfs.append((g, el.seqid, el.ltr_left[0] - 1, el.ltr_left[1], el.ltr_right[0] - 1, el.ltr_right[1]))
    pair = None
    ys = []
    gd = cat.g(D)
    for y, m in pmask.items():
        if not m[i]:
            continue
        gy = cat.g(y)
        ys.append((y, cat.seqids[y][int(cat.pseq[D][i, gy])], int(cat.pstart[D][i, gy]) - 1 - WINDOW_PAD,
                   int(cat.pend[D][i, gy]) + WINDOW_PAD, int(cat.pbstrand[D][i, gy])))
    if ys:
        el = elements[D][cat.eid[i, gd]]
        x0 = max(0, el.start - 1 - WINDOW_PAD)
        ltrs = tuple((a - 1 - x0, b - x0) for a, b in (el.ltr_left, el.ltr_right))
        pair = (D, el.seqid, x0, el.end + WINDOW_PAD, ltrs, tuple(ys))
    return tuple(selfs), pair


class _Tables:
    """The aligned counts of one trio's loci: self counts per (locus, genome), pairs per partner."""

    def __init__(self, loci: np.ndarray, jobs: list[tuple], results: list[tuple[tuple, tuple]],
                 D: str):
        self.D = D
        self.loci = loci
        self.self_: dict[tuple[int, str], tuple[int, int]] = {}  # (locus, genome) -> (sites, diffs), >= MIN_SEG
        self.short: Counter = Counter()  # genome -> copies whose LTRs align over < MIN_SEG sites
        self.pairs: dict[str, list[tuple[int, int, int]]] = {}   # partner -> (locus, sites, diffs), in order
        self.n_skipped = 0  # pairs whose call has no usable b_strand
        self.per_locus: list[tuple[list[str], list[str]]] = []  # each locus's self genomes and pair partners
        for i, (selfs, pair), (rs, rp) in zip(loci.tolist(), jobs, results, strict=True):
            gs, ps = [], []
            for (g, *_), (sites, diffs) in zip(selfs, rs, strict=True):
                if sites >= MIN_SEG:
                    self.self_[(i, g)] = (sites, diffs)
                    gs.append(g)
                else:
                    self.short[g] += 1
            if pair is not None:
                for (y, *_), (sites, diffs, o) in zip(pair[5], rp, strict=True):
                    if o == 0:
                        self.n_skipped += 1
                        continue
                    self.pairs.setdefault(y, []).append((i, sites, diffs))
                    ps.append(y)
            self.per_locus.append((gs, ps))
        self._pair_at = {(i, y): (sites, diffs) for y, v in self.pairs.items() for i, sites, diffs in v}

    def fit_pairs(self, y: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """(locus, k, n) of the pairs with partner y and >= IDENT_MIN_SITES LTR sites, in catalog order."""
        v = [t for t in self.pairs.get(y, []) if t[1] >= IDENT_MIN_SITES]
        idx = np.array([t[0] for t in v], dtype=np.int64)
        return idx, np.array([t[2] for t in v], dtype=np.float64), np.array([t[1] for t in v], dtype=np.float64)

    def self_u95(self, g: str, n_loci: int) -> tuple[np.ndarray, np.ndarray]:
        """(n,) mask of the loci where g's copy has a self row, and (n,) u95 there (nan elsewhere)."""
        has = np.zeros(n_loci, dtype=bool)
        bound = np.full(n_loci, np.nan)
        keys = [(i, c) for (i, gg), c in self.self_.items() if gg == g]
        if keys:
            idx = np.array([i for i, _ in keys], dtype=np.int64)
            has[idx] = True
            bound[idx] = u95(np.array([c[1] for _, c in keys]), np.array([c[0] for _, c in keys]))
        return has, bound

    def out(self, cat: Catalog, scope: str, young: np.ndarray) -> list[tuple]:
        """identity.tsv rows (IDENTITY_COLUMNS) in catalog order: each locus's self rows, then its pair rows."""
        rows: list[tuple] = []
        selfs = [(i, g) for i, (gs, _) in zip(self.loci.tolist(), self.per_locus) for g in gs]
        bounds = u95(np.array([self.self_[k][1] for k in selfs]), np.array([self.self_[k][0] for k in selfs])
                     ) if selfs else np.zeros(0)
        b = dict(zip(selfs, bounds.tolist()))
        for i, (gs, ps) in zip(self.loci.tolist(), self.per_locus):
            lid = str(cat.lid[i])
            for g in gs:
                sites, diffs = self.self_[(i, g)]
                yg = int(young[i]) if g == self.D and young[i] >= 0 else None
                rows.append(("self", lid, scope, g, ".", sites, diffs, diffs / sites, b[(i, g)], yg))
            for y in ps:
                sites, diffs = self._pair_at[(i, y)]
                div = diffs / sites if sites >= IDENT_MIN_SITES else None
                rows.append(("pair", lid, scope, self.D, y, sites, diffs, div, None, None))
        return rows


def _intro(trio: Trio, cat: Catalog, loci: np.ndarray, pmask: dict[str, np.ndarray], counted: np.ndarray,
           rd_class: np.ndarray, tab: _Tables) -> str:
    """The section's first line: what was aligned, at which loci, with what result."""
    R, S, D = trio.R, trio.S, trio.D
    n_pl = int(np.count_nonzero(np.logical_or.reduce(list(pmask.values()))))
    pairs = "; ".join(
        f"{y}-{D} {len(tab.pairs.get(y, []))} aligned, {sum(t[1] >= IDENT_MIN_SITES for t in tab.pairs.get(y, []))} "
        f"with >= {IDENT_MIN_SITES} LTR sites" for y in pmask)
    selfs = Counter(g for _, g in tab.self_)
    copies = ", ".join(f"{g} {selfs[g]}" + (f" (and {tab.short[g]} whose LTRs align over < {MIN_SEG} sites)"
                                            if tab.short[g] else "") for g in cat.labels if g in selfs or tab.short[g])
    return (f"Sequence-identity layer of {trio.scope}: {D}'s copies (+- {WINDOW_PAD} bp of flank, the engine's "
            f"window padding) realigned with the engine's aligner (align_window, {ALIGNER}) against their "
            f"sites in {R} and {S} (+- {WINDOW_PAD} bp) from {D}'s own P calls, oriented by those calls' b_strand, "
            f"and each annotated copy's 5' LTR against its 3' LTR; substitutions counted exactly in the LTR columns "
            f"(N and other ambiguity codes skipped), at {len(loci)} loci ({n_pl} with a P call of {D} into {R} "
            f"or {S}, {int(np.count_nonzero(counted))} D_pav counted, {int(np.count_nonzero(rd_class))} in "
            f"{R}'s RD painting class). Pairs: {pairs or 'none'}; self rows (LTR-LTR, >= {MIN_SEG} sites): "
            f"{copies or 'none'}")


def _split(cat: Catalog, trio: Trio, c, size: int, sd: tuple[np.ndarray, np.ndarray, np.ndarray],
           f: Findings) -> tuple[list[StatRow], float]:
    """Rows identity_split (scope S;D): the isolation fit with its jackknife se, the mixture guard, the fit
    check and n_pairs; and s."""
    S, D = trio.S, trio.D
    idx, k, n = sd
    scope, loci = f"{S};{D}", f"disc={D};P={S};sites>={IDENT_MIN_SITES}"
    s, m = fit_split(k, n)
    rows_mask = np.zeros(len(cat.lid), dtype=bool)
    rows_mask[idx] = True
    blocks = site_blocks(cat, c.gstar, D, rows_mask, size)  # catalog order, as idx
    cache: dict[int, tuple[float, float]] = {}

    def refit(keep: np.ndarray) -> tuple[float, float]:
        gone = blocks[~keep]
        if not gone.size:
            return s, m
        b = int(gone[0])
        if b not in cache:
            cache[b] = fit_split(k[keep], n[keep], start=(s, m))
        return cache[b]
    why = report.why({f"{S}-{D} pair": len(n)}, blocks)
    jk_s, why_s = report.guard(jackknife(lambda keep: refit(keep)[0], blocks), why)
    jk_m, why_m = report.guard(jackknife(lambda keep: refit(keep)[1], blocks), why)
    w, mg, s_mix, m_mix, mh, nll_mix = fit_mix(k, n, s, m)
    par, empty = _mix_params(w, mg, s_mix, m_mix, mh, len(n))
    aic_mix, aic_iso = 2 * nll_mix + 12, 2 * split_nll(np.log([s, m]), k, n) + 4
    jk_note = f"{jk_s.nblocks} jackknife blocks of {size:,} bp on {D}'s sites, each refit from the full fit"
    rows = [StatRow(LAYER, "identity_split", scope, "all", "s", s, jk_s.se, n=len(n), loci=loci,
                    note=f"isolation-model split, substitutions per LTR site; se: "
                         f"{jk_note}" + (f"; {NOT_TESTABLE}: {why_s}" if why_s else "")),
            StatRow(LAYER, "identity_split", scope, "all", "m", m, jk_m.se, n=len(n), loci=loci,
                    note=f"mean ancestral excess of S-D divergence over s; se: {jk_note}"
                         + (f"; {NOT_TESTABLE}: {why_m}" if why_m else ""))]
    na = f"; NA: the component is empty (weight x {len(n)} pairs < {MIX_MIN_PAIRS:g})"
    mix = (("s_mix", par["s_mix"], "split of the mixture guard (isolation + post-split gene flow + upper tail)"
            + (na if "isolation" in empty else "")),
           ("m_mix", par["m_mix"], "the mixture's isolation-component mean excess" + (na if "isolation" in empty
                                                                                     else "")),
           ("w_iso", float(w[0]), "the mixture's isolation weight"),
           ("w_flow", float(w[1]), f"the mixture's post-split gene-flow weight (component mean <= {MG_MAX})"),
           ("w_tail", float(w[2]), "the mixture's upper-tail weight (old paralogs, mis-assigned pairs)"),
           ("m_flow", par["m_flow"], "the gene-flow component's mean divergence" + (na if "gene flow" in empty
                                                                                   else "")),
           ("m_tail", par["m_tail"], "the upper-tail component's mean divergence" + (na if "tail" in empty
                                                                                    else "")),
           ("aic_mix", aic_mix, "2 nll + 12 (6 parameters)"),
           ("aic_iso", aic_iso, "the isolation fit's 2 nll + 4 (2 parameters)"))
    rows += [StatRow(LAYER, "identity_split", scope, "all", st, v, n=len(n), loci=loci, note=note)
             for st, v, note in mix]
    obs = {x: float(np.mean(k / n <= x)) for x in CDF_POINTS}
    fit = {x: pred_cdf(x, n, s, m) for x in CDF_POINTS}
    for x in CDF_POINTS:
        rows += [StatRow(LAYER, "identity_split", scope, "all", f"cdf_obs_{x:g}", obs[x], n=len(n), loci=loci,
                         note=f"observed fraction of pairs with diffs/sites <= {x:g}"),
                 StatRow(LAYER, "identity_split", scope, "all", f"cdf_fit_{x:g}", fit[x], n=len(n), loci=loci,
                         note=f"the isolation fit's predicted fraction <= {x:g} (200 simulations at the observed "
                              "sites, seed 0)")]
    rows.append(StatRow(LAYER, "identity_split", scope, "all", "n_pairs", len(n), n=len(n), loci=loci,
                        note=f"{S}-{D} ortholog pairs (P copies in {D}'s own calls) with >= {IDENT_MIN_SITES} "
                             "aligned LTR sites"))
    se = (f"({NOT_TESTABLE}: {why_s})" if why_s else f"± {report.fmt(jk_s.se, '.3g')} (jackknife over "
          f"{jk_s.nblocks} blocks of {size:,} bp on {D}'s sites, each refit from the full fit)")
    f.add(SECTION, f"Split {S};{D}, rate-free (isolation model: {S}-{D} LTR divergence = s + X, X ~ "
                   f"Exp(mean m), in substitutions per LTR site) on {len(n)} {S}-{D} pairs with >= "
                   f"{IDENT_MIN_SITES} LTR sites: s = {report.fmt(s, '.4g')} {se}; m = {report.fmt(m, '.4g')}"
                   + (f" ± {report.fmt(jk_m.se, '.3g')}" if not why_m else ""))
    better = "fits better than isolation alone" if aic_mix < aic_iso else "does not fit better than isolation"
    empty_txt = (f"; empty, its parameters NA (weight x {len(n)} pairs < {MIX_MIN_PAIRS:g}, its weight's lower "
                 f"bound): {', '.join(empty)}" if empty else "")
    tail_txt = ("" if {"isolation", "tail"} & set(empty) else
                f", and the tail is an upper tail only when its mean ({mh:.4g}) exceeds the isolation component's "
                f"(s_mix + m_mix = {s_mix + m_mix:.4g})")
    f.add(SECTION, f"Guard against post-split gene flow (18 starts): isolation + gene flow (mean "
                   f"<= {MG_MAX:g}) + upper tail, weights {w[0]:.3f} / {w[1]:.4f} / {w[2]:.3f}, gene-flow mean "
                   f"{_na(par['m_flow'])}, s_mix = {_na(par['s_mix'])}, m_mix = {_na(par['m_mix'])}, tail mean "
                   f"{_na(par['m_tail'])}{empty_txt}; AIC {aic_mix:.1f} against {aic_iso:.1f}: the mixture "
                   f"{better}. s_mix is the shift of the isolation component alone (the gene-flow and tail "
                   "components are unshifted): it corrects for post-split gene flow only through the gene-flow "
                   f"weight ({w[1]:.4f}){tail_txt}. Post-split {S}-{D} gene flow pulls the isolation s down")
    f.add(SECTION, f"Fit check of the isolation fit (fraction of {S}-{D} pairs with divergence <= x, observed "
                   "against fitted): " + "; ".join(f"{x:g}: {obs[x]:.4f} / {fit[x]:.4f}" for x in CDF_POINTS))
    return rows, s


def _mix_params(w: np.ndarray, mg: float, s: float, m: float, mh: float, n_pairs: int
                ) -> tuple[dict[str, float], list[str]]:
    """The mixture's component parameters {s_mix, m_mix, m_flow, m_tail}, nan for a component at its weight's
    lower bound, and the names of those components (isolation, gene flow, tail). The weights are a softmax of
    free logits, so Nelder-Mead drives an unused component's weight toward 0 without reaching it: a weight
    worth fewer than MIX_MIN_PAIRS of the n pairs (weight x n_pairs < MIX_MIN_PAIRS) is taken as that bound,
    and the component's parameters then rest on no data (the fit can return a tail mean of 6e284 there)."""
    empty = [float(w[j]) * n_pairs < MIX_MIN_PAIRS for j in range(3)]
    par = {"s_mix": math.nan if empty[0] else float(s), "m_mix": math.nan if empty[0] else float(m),
           "m_flow": math.nan if empty[1] else float(mg), "m_tail": math.nan if empty[2] else float(mh)}
    return par, [name for name, e in zip(("isolation", "gene flow", "tail"), empty, strict=True) if e]


def _na(v: float) -> str:
    """A fitted parameter for findings.txt: NA where it is nan (an empty component)."""
    return "NA" if math.isnan(v) else f"{v:.4g}"


def _young(cat: Catalog, trio: Trio, c, size: int, s: float, tab: _Tables, f: Findings
           ) -> tuple[list[StatRow], np.ndarray]:
    """Rows identity_young (trio scope) on the D_pav counted rows, and young per locus: 1 or 0 where D's copy
    at a counted R+D or S+D locus has a self row, -1 elsewhere."""
    R, S, D, scope = trio.R, trio.S, trio.D, trio.scope
    a, b = patterns(cat, c.gstar, R, S)
    has, bound = tab.self_u95(D, len(cat.lid))
    rd, sd = c.rows & a & has, c.rows & b & has
    is_young = np.zeros(len(cat.lid), dtype=bool)
    is_young[has] = bound[has] < s
    young = np.full(len(cat.lid), -1, dtype=np.int8)
    young[rd | sd] = is_young[rd | sd]
    tested = rd | sd
    ay, by = (rd & is_young)[tested], (sd & is_young)[tested]
    n_ry, n_sy = int(np.count_nonzero(ay)), int(np.count_nonzero(by))
    loci = f"{c.loci};self={D}"
    sites_needed = float(u95(np.array([0]), np.array([1.0]))[0]) / s
    rows = [StatRow(LAYER, "identity_young", scope, "all", "n_RD_young", n_ry, n=int(np.count_nonzero(rd)),
                    loci=loci, note=f"{R} present, {S} E, {D}'s copy younger than the split (u95 < s; n: those "
                                    f"with a self row of {D})"),
            StatRow(LAYER, "identity_young", scope, "all", "n_SD_young", n_sy, n=int(np.count_nonzero(sd)),
                    loci=loci, note=f"{S} present, {R} E, {D}'s copy younger than the split (u95 < s; n: those "
                                    f"with a self row of {D})")]
    if tested.any():
        blocks = site_blocks(cat, c.gstar, D, tested, size)
        why = report.why({f"young {R}+{D}": n_ry, f"young {S}+{D}": n_sy}, blocks[ay | by])
        jk, why = report.guard(jackknife(d_ratio(ay, by), blocks), why)
        rows.append(StatRow(LAYER, "identity_young", scope, "all", "D_young", jk.est, jk.se, jk.z, n=n_ry + n_sy,
                            loci=loci, note=f"(n_RD_young - n_SD_young) / (n_RD_young + n_SD_young); "
                                            f"{jk.nblocks} jackknife blocks of {size:,} bp on {D}'s sites, s held "
                                            "at its estimate" + (f"; {NOT_TESTABLE}: {why}" if why else "")))
        ztext = f"{NOT_TESTABLE}: {why}" if why else f"± {report.fmt(jk.se, '.3f')}, Z {report.fmt(jk.z, '.2f')}"
        dline = f"D_young = {report.fmt(jk.est, '+.3f')} ({ztext})"
    else:
        dline = f"D_young not computed: no counted {R}+{D} or {S}+{D} locus whose {D} copy has a self row"
    rows.append(StatRow(LAYER, "identity_young", scope, "all", "sites_needed", sites_needed, loci="none",
                        note="identical LTR sites a single insertion needs for u95 < s: chi2.ppf(0.95, 2) / 2 / s"))
    f.add(SECTION, f"Younger than the split (u95 of {D}'s own LTR-LTR divergence < s = {report.fmt(s, '.4g')}; D_pav "
                   f"counted loci whose {D} copy has a self row): {R}+{D} {n_ry} of {int(np.count_nonzero(rd))}, "
                   f"{S}+{D} {n_sy} of {int(np.count_nonzero(sd))}; {dline}. An insertion needs >= "
                   f"{sites_needed:,.0f} identical LTR sites to qualify. Without gene flow, an insertion younger "
                   f"than the split cannot be shared by {D} and one sister; LTR-LTR gene conversion makes some old "
                   "insertions look young, in both classes alike")
    return rows, young


def _below(trio: Trio, s: float, rdp: tuple[np.ndarray, np.ndarray, np.ndarray],
           sdp: tuple[np.ndarray, np.ndarray, np.ndarray], f: Findings) -> list[StatRow]:
    """Rows identity_below_split (trio scope): the fractions of R-D and S-D pairs with diffs / sites < s."""
    R, S, D, scope = trio.R, trio.S, trio.D, trio.scope
    rows, parts = [], []
    for stat, y, (_, k, n) in (("frac_RD_below", R, rdp), ("frac_SD_below", S, sdp)):
        frac = float(np.mean(k / n < s)) if n.size else math.nan
        rows.append(StatRow(LAYER, "identity_below_split", scope, "all", stat, frac, n=int(n.size),
                            loci=f"disc={D};P={y};sites>={IDENT_MIN_SITES}",
                            note=f"the fraction of {y}-{D} pairs with diffs/sites < s"))
        parts.append(f"{y}-{D} {report.fmt(frac, '.4f')} of {n.size}")
    f.add(SECTION, f"Ortholog pairs below the split (diffs/sites < s, >= {IDENT_MIN_SITES} LTR sites): "
                   + ", ".join(parts))
    return rows
