"""Anchoring: raw blocks of genome A vs genome B from a route, chained and netted to syntenic 1:1 blocks.

A route maps one ordered genome pair to raw, unchained blocks (BLOCK_DTYPE, 0-based half-open). `anchor_pair`
chains and nets them (`chain_blocks`, `net_blocks`) into the orthologous block set that windows are built
from, cached per ordered pair and route under <workdir>/anchors. Two routes:
  flank-mm2  (FlankMm2Route, the default) the unique sequence beside each annotated element of A, cut into
             pieces (`flank_pieces`) and mapped onto B's minimap2 asm20 index (index) in-process with mappy
  paf        (PafRoute, `--paf DIR`) the pair's own whole-genome alignments, <dir>/<A>__<B>.paf[.gz]

`flank_pieces` cuts, for each top-level element of A and each side, the free sequence (no annotated element,
no N run) within FLANK_CAP bp, nearest first until FLANK_TARGET bp are collected. minimap2 aligns through
indels up to its 10-kb chaining gap (an element insertion included), so each hit is cut at its indels of
>= SPLIT_INDEL bp (`blocks.split_cigar`): a TE-sized indel lies between blocks, where chaining bridges it, never
inside one, where the linear A-to-B map and the divergence would absorb it.
"""
from __future__ import annotations

import bisect
import logging
import math
import threading
import time
import zipfile
from collections import deque
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from itertools import accumulate
from pathlib import Path
from typing import Protocol

import mappy
import numpy as np

from synltr.base.checkpoint import (
    atomic_path,
    clear_done,
    code_fingerprint,
    file_checksum,
    is_done,
    mark_done,
    marker,
    stage_key,
)
from synltr.base.constants import (
    ANCHOR_BEST_N,
    ANCHOR_MIN,
    ANCHOR_PRESET,
    CHAIN_INDEL_COST,
    CHAIN_LOOKBACK,
    CHAIN_MAX_GAP,
    CHAIN_OVERLAP,
    CHECKSUM_CACHE,
    FLANK_CAP,
    FLANK_PIECE_MIN,
    FLANK_TARGET,
    MIN_CHAIN_SCORE,
    PIECES_QUEUED,
    SPLIT_INDEL,
)
from synltr.io.elements import Element, elements_digest
from synltr.io.samples import Sample
from synltr.io.seqio import Fasta, n_runs
from synltr.routes.assembly.blocks import (
    BLOCK_DTYPE,
    blocks_from_rows,
    load_blocks,
    read_paf,
    save_blocks,
    split_cigar,
)
from synltr.routes.assembly.chain import chain_blocks, net_blocks
from synltr.routes.assembly.index import index_versions, load_index

log = logging.getLogger(__name__)

# the netted blocks depend on these besides the raw blocks; "net" names the net's rule, so anchors netted
# under another rule are recomputed, never reused
_CHAIN_PARAMS = {"anchor_min": ANCHOR_MIN, "max_gap": CHAIN_MAX_GAP, "indel_cost": CHAIN_INDEL_COST,
                 "overlap": CHAIN_OVERLAP, "lookback": CHAIN_LOOKBACK, "min_score": MIN_CHAIN_SCORE,
                 "net": "syntenic"}


@dataclass
class PairContext:
    """One ordered genome pair: A (whose elements are genotyped) against B."""

    a: Sample
    b: Sample
    a_fa: Fasta
    b_fa: Fasta
    a_elements: list[Element]
    b_elements: list[Element] | None
    outdir: Path
    threads: int

    @property
    def workdir(self) -> Path:
        """Intermediate files and caches: <outdir>/work."""
        return Path(self.outdir) / "work"


class Route(Protocol):
    """An anchoring route: raw (unchained) alignment blocks of A against B. Every member is required (anchor_pair
    reads them directly).

    name: short id, used in cache file names.
    reads_elements: True if the blocks depend on ctx.a_elements / ctx.b_elements (then part of the cache key).
    params(ctx): JSON values the blocks depend on besides the two FASTAs (settings, tool versions, checksums of
        external inputs); part of the cache key.
    raw_blocks(ctx): BLOCK_DTYPE array; A indices follow ctx.a_fa's fai order, B indices ctx.b_fa's.
    """

    name: str
    reads_elements: bool

    def params(self, ctx: PairContext) -> dict: ...

    def raw_blocks(self, ctx: PairContext) -> np.ndarray: ...


def anchor_pair(ctx: PairContext, route: Route) -> tuple[np.ndarray, list[str], list[str], float]:
    """Netted syntenic 1:1 blocks of A vs B by `route`, the A and B seqids (fai order) and the divergence d.

    d = 1 - sum(matches) / sum(alen) over the netted blocks (nan if there are none). The result is cached at
    <workdir>/anchors/<A>__<B>.<route>.npz, keyed by both FASTA checksums, the elements (routes that read
    them), the route's params, the chaining constants with the net's rule (_CHAIN_PARAMS) and the code of this
    module's import closure.
    """
    out = anchors_path(ctx, route)
    anchors = out.parent
    done = marker(out)
    sums = ctx.workdir / CHECKSUM_CACHE
    inputs = {"A": file_checksum(ctx.a_fa.path, sums), "B": file_checksum(ctx.b_fa.path, sums)}
    if route.reads_elements:
        inputs["A_elements"] = elements_digest(ctx.a_elements)
        inputs["B_elements"] = elements_digest(ctx.b_elements)
    key = stage_key("anchor", inputs, {"route": route.name, "route_params": route.params(ctx),
                                       "chain": _CHAIN_PARAMS}, code_fingerprint("routes.assembly.anchor"))
    pair = f"{ctx.a.label} vs {ctx.b.label} ({route.name})"
    if is_done(done, key):
        cached = _load(out, pair)
        if cached is not None:
            blocks, a_names, b_names = cached
            d = divergence(blocks)
            log.info("%s: anchors done, skip (%d netted blocks, d = %.4f)", pair, len(blocks), d)
            return blocks, a_names, b_names, d

    clear_done(done)
    anchors.mkdir(parents=True, exist_ok=True)
    log.debug("%s: anchoring", pair)
    raw = route.raw_blocks(ctx)
    blocks = net_blocks(chain_blocks(raw), label=pair)
    a_names, b_names = list(ctx.a_fa.lengths()), list(ctx.b_fa.lengths())
    _save(out, blocks, a_names, b_names)
    mark_done(done, key, [out])
    d = divergence(blocks)
    anchored = int((blocks["a_end"] - blocks["a_start"]).sum(dtype=np.int64))
    log.debug("%s: %d raw blocks -> %d netted, %.1f Mb of A anchored, d = %.4f",
              pair, len(raw), len(blocks), anchored / 1e6, d)
    return blocks, a_names, b_names, d


def anchors_path(ctx: PairContext, route: Route) -> Path:
    """Where `anchor_pair` caches the netted blocks: <workdir>/anchors/<A>__<B>.<route>.npz."""
    return ctx.workdir / "anchors" / f"{ctx.a.label}__{ctx.b.label}.{route.name}.npz"


def divergence(blocks: np.ndarray) -> float:
    """1 - sum(matches) / sum(alen); nan if the blocks hold no aligned column."""
    alen = int(blocks["alen"].sum(dtype=np.int64))
    return 1 - int(blocks["matches"].sum(dtype=np.int64)) / alen if alen > 0 else math.nan


def _load(path: Path, what: str) -> tuple[np.ndarray, list[str], list[str]] | None:
    """load_blocks(path), or None with a warning when the cached file cannot be read (a damaged file of the
    recorded size passes is_done): the caller recomputes it."""
    try:
        return load_blocks(path)
    except (OSError, ValueError, KeyError, EOFError, zipfile.BadZipFile) as exc:
        log.warning("%s: cached blocks %s are unreadable (%s: %s); recomputing", what, path, type(exc).__name__, exc)
        return None


def _save(path: Path, blocks: np.ndarray, a_names: list[str], b_names: list[str]) -> None:
    """save_blocks to a temp file, then rename: a reader never sees a partial file."""
    with atomic_path(path) as tmp:
        save_blocks(tmp, blocks, a_names, b_names)


# ---------------------------------------------------------------- the flank-mm2 route


@dataclass(frozen=True)
class Piece:
    """A flank piece of A: 0-based half-open."""

    seqid: str
    start: int
    end: int


def flank_pieces(elements: list[Element], fasta: Fasta, cap: int = FLANK_CAP,
                 target: int = FLANK_TARGET) -> list[Piece]:
    """Flank pieces of the top-level `elements` of A, deduplicated, in fai order then start.

    Per element and side the flank region is [start - cap, start) or [end, end + cap) clipped to the sequence.
    Every annotated element (all levels) and every N run is cut out; the free stretches of >= FLANK_PIECE_MIN
    bp are taken nearest-first, whole, until they total >= `target` bp (the last may overshoot).
    """
    lengths = fasta.lengths()
    by_seq: dict[str, list[tuple[int, int]]] = {}
    for e in elements:
        by_seq.setdefault(e.seqid, []).append((e.start - 1, e.end))
    # per seqid: spans sorted by start, their starts, and the running max of ends (non-decreasing, so nested
    # spans cannot hide an overlap): a region [lo, hi) can only meet spans[j:k]
    spans = {}
    for sid, lst in by_seq.items():
        lst.sort()
        spans[sid] = (lst, [s for s, _ in lst], list(accumulate((t for _, t in lst), max)))
    out: set[Piece] = set()
    for e in elements:
        if e.level != 0:
            continue
        s0, t0, n = e.start - 1, e.end, lengths[e.seqid]
        for lo, hi, left in ((max(0, s0 - cap), s0, True), (t0, min(n, t0 + cap), False)):
            if hi <= lo:
                continue
            lst, starts, max_end = spans[e.seqid]
            j, k = bisect.bisect_right(max_end, lo), bisect.bisect_left(starts, hi)
            cuts = [(max(s, lo) - lo, min(t, hi) - lo) for s, t in lst[j:k] if t > lo]
            cuts += n_runs(fasta.fetch(e.seqid, lo, hi))
            free, pos = [], 0
            for s, t in sorted(cuts):
                if s - pos >= FLANK_PIECE_MIN:
                    free.append((pos, s))
                pos = max(pos, t)
            if hi - lo - pos >= FLANK_PIECE_MIN:
                free.append((pos, hi - lo))
            total = 0
            for s, t in (reversed(free) if left else free):
                if total >= target:
                    break
                out.add(Piece(e.seqid, lo + s, lo + t))
                total += t - s
    order = {name: i for i, name in enumerate(lengths)}
    return sorted(out, key=lambda p: (order[p.seqid], p.start, p.end))


class FlankMm2Route:
    """Flank pieces of A mapped onto B's minimap2 index (built once per B FASTA)."""

    name = "flank-mm2"
    reads_elements = True

    def __init__(self, cap: int = FLANK_CAP):
        self.cap = cap

    def params(self, ctx: PairContext) -> dict:
        return {"cap": self.cap, "target": FLANK_TARGET, "piece_min": FLANK_PIECE_MIN, "preset": ANCHOR_PRESET,
                "best_n": ANCHOR_BEST_N, "split_indel": SPLIT_INDEL, **index_versions()}

    def raw_blocks(self, ctx: PairContext) -> np.ndarray:
        return map_pieces(ctx, flank_pieces(ctx.a_elements, ctx.a_fa, self.cap, FLANK_TARGET))


def map_pieces(ctx: PairContext, pieces: list[Piece]) -> np.ndarray:
    """Raw blocks of `pieces` (A) mapped with mappy onto B's asm20 index.

    Pieces are fetched in this thread (pysam handles are not thread-safe) and mapped by a thread pool of
    ctx.threads, one mappy ThreadBuffer per thread, at most PIECES_QUEUED pieces in flight per thread; results
    come back in piece order, so the output is deterministic.
    """
    t0 = time.perf_counter()
    if not pieces:
        return np.zeros(0, dtype=BLOCK_DTYPE)
    aligner = load_index(ctx.b.label, ctx.b_fa, ctx.workdir, ctx.threads)
    a_index = {n: i for i, n in enumerate(ctx.a_fa.lengths())}
    b_index = {n: i for i, n in enumerate(ctx.b_fa.lengths())}
    local = threading.local()

    def work(job: tuple[Piece, str]) -> list[tuple]:
        piece, seq = job
        buf = getattr(local, "buf", None)
        if buf is None:
            buf = local.buf = mappy.ThreadBuffer()
        return map_rows(aligner, buf, seq, a_index[piece.seqid], piece.start, b_index)

    threads = max(1, ctx.threads)
    parts: list[list[tuple]] = []
    pending: deque = deque()
    with ThreadPoolExecutor(threads) as pool:
        for p in pieces:  # fetched just before submission; at most PIECES_QUEUED * threads pieces in flight
            pending.append(pool.submit(work, (p, ctx.a_fa.fetch(p.seqid, p.start, p.end))))
            if len(pending) >= PIECES_QUEUED * threads:  # collected in piece order: deterministic
                parts.append(pending.popleft().result())
        parts.extend(f.result() for f in pending)
    rows = [r for part in parts for r in part]
    raw = blocks_from_rows(rows) if rows else np.zeros(0, dtype=BLOCK_DTYPE)
    log.debug("%s vs %s (flank-mm2): %d pieces, %d raw blocks, %.1f s", ctx.a.label, ctx.b.label, len(pieces),
              len(raw), time.perf_counter() - t0)
    return raw


def map_rows(aligner: mappy.Aligner, buf: mappy.ThreadBuffer, seq: str, ai: int, w0: int,
             b_index: dict[str, int]) -> list[tuple]:
    """Block rows (A seqid index, A start, A end, B seqid index, B start, B end, strand, matches, alen) of the
    hits of `seq` (A[w0, w0 + len(seq)) of seqid index `ai`), each cut at its long indels (`_pieces`)."""
    rows = []
    for h in aligner.map(seq, buf=buf):
        bi = b_index[h.ctg]
        for q0, q1, r0, r1, matches, alen in _pieces(h.cigar, h.q_st, h.r_st, h.r_en, h.strand):
            rows.append((ai, w0 + q0, w0 + q1, bi, r0, r1, h.strand, matches, alen))
    return rows


def _pieces(cigar: list, q_st: int, r_st: int, r_en: int, strand: int) -> list[tuple[int, ...]]:
    """One mappy hit (its CIGAR as mappy's [length, op] pairs) as (q0, q1, r0, r1, matches, alen) pieces in query
    order: `split_cigar` at SPLIT_INDEL. The aligner reports =/X (EQX), so matches are exact '=' counts; an M op
    raises ValueError."""
    ops = np.array([op for _, op in cigar], dtype=np.int64)
    lens = np.array([n for n, _ in cigar], dtype=np.int64)
    pieces, has_m = split_cigar(ops, lens, q_st, r_st, r_en, strand, SPLIT_INDEL)
    if has_m:
        raise ValueError("CIGAR op M: mappy must report =/X ops (EQX)")
    return pieces


# ---------------------------------------------------------------- the PAF route


class PafRoute:
    """Raw blocks read from an external whole-genome PAF per genome pair: <dir>/<A>__<B>.paf[.gz] (A = query, PAF
    column 1), else <dir>/<B>__<A>.paf[.gz] (read swapped)."""

    name = "paf"
    reads_elements = False

    def __init__(self, dirpath: Path):
        self.dirpath = Path(dirpath)

    def find(self, ctx: PairContext) -> tuple[Path, bool]:
        """The pair's PAF and whether its query is B (read swapped); .paf before .paf.gz, A__B before B__A."""
        a, b = ctx.a.label, ctx.b.label
        for query, target, swap in ((a, b, False), (b, a, True)):
            for ext in (".paf", ".paf.gz"):
                path = self.dirpath / f"{query}__{target}{ext}"
                if path.is_file():
                    return path, swap
        raise FileNotFoundError(f"no PAF for {a} vs {b}: neither {self.dirpath / f'{a}__{b}.paf[.gz]'} "
                                f"nor {self.dirpath / f'{b}__{a}.paf[.gz]'} exists")

    def params(self, ctx: PairContext) -> dict:
        path, _ = self.find(ctx)
        return {"paf": path.name, "sha256": file_checksum(path, ctx.workdir / CHECKSUM_CACHE),
                "split_indel": SPLIT_INDEL}

    def raw_blocks(self, ctx: PairContext) -> np.ndarray:
        path, swap = self.find(ctx)
        log.info("%s vs %s (paf): reading %s%s", ctx.a.label, ctx.b.label, path, " as B vs A" if swap else "")
        a_index = {name: i for i, name in enumerate(ctx.a_fa.lengths())}
        b_index = {name: i for i, name in enumerate(ctx.b_fa.lengths())}
        return read_paf(path, a_index, b_index, swap=swap, split_indel=SPLIT_INDEL)
