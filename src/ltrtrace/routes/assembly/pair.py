"""Pair genotyping: the call file of one ordered genome pair (A, whose elements are called; B, the target).

`genotype_pair` runs anchor_pair -> make_windows -> per window: the A and B sequences (B reverse-complemented on a
minus-strand window), align_window (LASTZ, asm20 scoring) and classify_window -> one row per element of A. The
windows, in A fai order, are cut into chunks (`_chunk_size`: about CHUNK_TARGET per pair, whatever the thread
count), genotyped through base.runtime.pool_map, every chunk in a forked worker process, largest first, all at once
(each worker opens both FASTAs itself: pysam handles must not be shared across processes), or all in this process
with one thread. A chunk that fails, or a worker that dies (killed, out of memory, crashed), stops the pair with an
error naming it: pending chunks are canceled, running ones finish first, finished ones stay done. A worker exits when
the parent dies, killing its running lastz child and removing its temp folder (lastz.kill_lastz); the folders of
workers killed outright are swept at the next pair.

Outputs: gzip TSV with a header line, 1-based closed coordinates and reproducible bytes (calls.records):
  <outdir>/pairs/<A>__<B>.calls.tsv.gz       CALL_COLUMNS, every element of A (all levels) exactly once,
                                              sorted by (A fai order, start, -end, level). An element in no
                                              window (its parent is not annotated) is U no_anchor, wid ".".
  <workdir>/windows/<A>__<B>.tsv.gz          WINDOW_COLUMNS, one row per window (wid order).
  <workdir>/calls/<A>__<B>/chunk_<n>.tsv.gz  the rows of chunk n (same columns and order as the pair file).
The pair file and each chunk have a done marker keyed by the anchors, both FASTAs, the elements and the lastz
version (a chunk also by the windows and its wids), and by the code fingerprint of this module's import closure
(every module a pair runs): a pair finished by other genotyping code or another aligner is recomputed. A rerun
skips a finished pair; an interrupted one recomputes only its unfinished chunks, then merges the chunks' lines
again, as they are. The pair's marker also records the engine ("assembly"), the digest of A's elements and the
code fingerprint.
"""
from __future__ import annotations

import heapq
import logging
import math
import os
import time
from collections import Counter
from collections.abc import Callable, Iterator
from dataclasses import dataclass
from pathlib import Path

from ltrtrace.base.checkpoint import clear_done, code_fingerprint, file_checksum, is_done, marker, stage_key
from ltrtrace.base.constants import CHECKSUM_CACHE, CHUNK_MIN, CHUNK_TARGET, CHUNK_WINDOWS, U_WARN
from ltrtrace.base.runtime import PoolError, pool_map
from ltrtrace.calls.records import CALL_COLUMNS, call_line, call_lines, gz_text, write_call_lines
from ltrtrace.calls.states import PAIR_STATES
from ltrtrace.io.elements import Element, elements_digest
from ltrtrace.io.seqio import Fasta, revcomp
from ltrtrace.routes.assembly.anchor import PairContext, Route, anchor_pair, anchors_path
from ltrtrace.routes.assembly.classify import call_row, classify_window
from ltrtrace.routes.assembly.lastz import align_window, kill_lastz, lastz_version, pid_alive, sweep_lastz_tmp
from ltrtrace.routes.assembly.window import Window, make_windows

log = logging.getLogger(__name__)

WINDOW_COLUMNS = ("wid", "a_seqid", "a_start", "a_end", "b_seqid", "b_start", "b_end", "strand", "status",
                  "left_anchor", "right_anchor", "margin", "eids")
ENGINE = "assembly"  # the meta of every call file this route writes (calls.records)


def genotype_pair(ctx: PairContext, route: Route) -> Path:
    """Call every element of A in B; returns the call file <outdir>/pairs/<A>__<B>.calls.tsv.gz.

    `route` anchors the pair (anchor_pair); LASTZ aligns each window (align_window).
    """
    t0 = time.perf_counter()
    pair, label = f"{ctx.a.label}__{ctx.b.label}", f"{ctx.a.label} vs {ctx.b.label}"
    out = Path(ctx.outdir) / "pairs" / f"{pair}.calls.tsv.gz"
    order = {s: i for i, s in enumerate(ctx.a_fa.lengths())}  # A fai order
    elements, children = _index(ctx.a_elements, order)  # fails fast, before the anchoring
    blocks, a_names, b_names, _ = anchor_pair(ctx, route)
    sums = ctx.workdir / CHECKSUM_CACHE
    # the anchors alone miss what they do not hold: element fields (a route may not read elements at all)
    # and sequence outside the blocks
    a_digest = elements_digest(ctx.a_elements)
    inputs = {"anchors": file_checksum(anchors_path(ctx, route), sums), "A": file_checksum(ctx.a_fa.path, sums),
              "B": file_checksum(ctx.b_fa.path, sums), "elements": a_digest}
    params = {"pair": pair, "aligner": lastz_version()}
    code = code_fingerprint("routes.assembly.pair")  # every module the pair runs (its import closure)
    key = stage_key("genotype_pair", inputs, {**params, "route": route.name}, code)
    if is_done(marker(out), key):
        log.info("%s: calls done, skip (%s)", label, out)
        return out
    clear_done(marker(out))

    windows = make_windows(ctx.a_elements, blocks, a_names, b_names, ctx.b_fa.lengths(), ctx.a_fa.lengths())
    wpath = ctx.workdir / "windows" / f"{pair}.tsv.gz"
    for path in (out, wpath):
        _sweep_tmp(path)
    with gz_text(wpath) as fh:
        fh.write(call_line(WINDOW_COLUMNS))
        fh.writelines(map(_window_line, windows))
    chunk_inputs = {**inputs, "windows": file_checksum(wpath, sums)}
    chunks, jobs, size = [], [], _chunk_size(len(windows))
    for n, i in enumerate(range(0, len(windows), size)):
        part = windows[i:i + size]
        path = ctx.workdir / "calls" / pair / f"chunk_{n}.tsv.gz"
        ckey = stage_key("genotype_chunk", chunk_inputs,
                         {**params, "chunk": n, "wids": [part[0].wid, part[-1].wid]}, code)
        chunks.append(path)
        if not is_done(marker(path), ckey):
            clear_done(marker(path))
            jobs.append(_Job(n, part, path, ckey))
    cdir = ctx.workdir / "calls" / pair
    for extra in cdir.glob("chunk_*.tsv.gz*") if cdir.is_dir() else ():  # chunks of an earlier, longer run
        n = extra.name.split(".", 1)[0].removeprefix("chunk_")
        if n.isdigit() and int(n) >= len(chunks):
            extra.unlink(missing_ok=True)
    if jobs:
        for tmp in cdir.glob(".chunk_*.tmp") if cdir.is_dir() else ():  # left by a worker killed mid-write
            tmp.unlink(missing_ok=True)
        # the most work first (LPT): a large chunk submitted last becomes the pair's straggler
        jobs.sort(key=lambda j: -sum(w.a_end - w.a_start + w.b_end - w.b_start for w in j.windows
                                     if w.status == "ok"))
        _run(ctx, jobs, elements, children, label)

    alone = sorted(map(_no_window_line, _in_no_window(windows, elements, children)), key=_sort_key(order))
    meta = {"engine": ENGINE, "key": key, "elements": a_digest, "code": code}
    states = _merge(chunks, alone, out, elements, order, label, meta)

    total = sum(states.values())
    summary = ", ".join(f"{states[s]} {s}" for s in PAIR_STATES if states[s]) or "none"
    log.info("%s: %d windows (%d anchored), %d chunk(s) (%d computed), %d calls (%s) -> %s, %.1f s", label,
             len(windows), sum(w.status == "ok" for w in windows), len(chunks), len(jobs), total, summary, out,
             time.perf_counter() - t0)
    if total and states["U"] > U_WARN * total:
        log.warning("%s: %d of %d calls are U: no orthologous alignment between %s and %s?", label,
                    states["U"], total, ctx.a.label, ctx.b.label)
    return out


def _chunk_size(n_windows: int) -> int:
    """Windows per chunk: max(CHUNK_MIN, min(CHUNK_WINDOWS, ceil(n_windows / CHUNK_TARGET))).

    About CHUNK_TARGET chunks per pair balance the workers; the size never depends on the thread count, so a
    resumed run finds its chunks (keyed by index and wids) at any --threads. The calls do not depend on it.
    """
    return max(CHUNK_MIN, min(CHUNK_WINDOWS, -(-n_windows // CHUNK_TARGET)))


# ---------------------------------------------------------------- workers


@dataclass
class _Job:
    """Chunk n: its windows, its output and its stage key."""

    n: int
    windows: list[Window]
    path: Path
    key: str

    def __repr__(self) -> str:
        return f"chunk {self.n}"


class _Genotyper:
    """Genotypes chunks in one process: both FASTAs open, the elements and children map."""

    def __init__(self, a_fa: Fasta, b_fa: Fasta, elements: dict[str, Element], children: dict[str, list[str]]):
        self.a_fa, self.b_fa = a_fa, b_fa
        self.elements, self.children = elements, children
        self.sort_key = _sort_key({s: i for i, s in enumerate(a_fa.lengths())})

    def rows(self, win: Window) -> list[dict]:
        """The calls of the window's elements as CALL_COLUMNS rows; sequences only for an ok window."""
        a_seq = b_seq = ""
        aln = None
        if win.status == "ok":
            a_seq = self.a_fa.fetch(win.a_seqid, win.a_start, win.a_end)
            b_seq = self.b_fa.fetch(win.b_seqid, win.b_start, win.b_end)
            if win.strand == -1:
                b_seq = revcomp(b_seq)  # into A's orientation
            aln = align_window(a_seq, b_seq)
        calls = classify_window(win, aln, a_seq, b_seq, self.elements, self.children)
        return [call_row(c, self.elements[c.eid], win) for c in calls]

    def chunk(self, job: _Job) -> tuple[int, int]:
        """Write the job's rows, sorted, to its chunk file and mark it done; returns (n, rows)."""
        lines = sorted((call_line(r[c] for c in CALL_COLUMNS) for w in job.windows for r in self.rows(w)),
                       key=self.sort_key)
        write_call_lines(job.path, lines, {"engine": ENGINE, "key": job.key})
        return job.n, len(lines)


def _run(ctx: PairContext, jobs: list[_Job], elements: dict[str, Element], children: dict[str, list[str]],
         label: str) -> None:
    """Genotype the jobs' chunks (largest first, as given) through pool_map with ctx.threads.

    Every chunk goes to the workers at once (warm=False: a worker loads the chunks' numba kernels from numba's cache,
    or compiles them, itself, so no chunk waits for this process); a worker opens both FASTAs itself at its first
    chunk. With one thread or one chunk, they run in this process with its open FASTAs. A chunk is marked done only
    by the process that finished it. A failing chunk raises a RuntimeError naming the pair and the chunk; a worker
    that dies (killed by a signal or the OOM killer, crashed, or never started) raises PoolError naming the pair and
    the chunk(s) it held.
    """
    threads = max(1, ctx.threads)
    if n := sweep_lastz_tmp():
        log.info("%s: removed %d LASTZ temp folder(s) left by killed workers", label, n)
    here, parent = _Genotyper(ctx.a_fa, ctx.b_fa, elements, children), os.getpid()
    opened: dict[int, _Genotyper] = {}  # a worker's genotyper, by pid (each worker fills its own copy)
    files = (ctx.a_fa.path, ctx.a_fa.cache_dir, ctx.b_fa.path, ctx.b_fa.cache_dir)

    def chunk(job: _Job) -> tuple[int, int]:
        try:
            pid = os.getpid()
            if pid == parent:
                g = here
            elif (g := opened.get(pid)) is None:
                a, a_cache, b, b_cache = files
                g = opened[pid] = _Genotyper(Fasta(a, a_cache), Fasta(b, b_cache), elements, children)
            return g.chunk(job)
        except Exception as exc:
            raise RuntimeError(f"{label}: chunk {job.n} failed: {type(exc).__name__}: {exc}") from exc

    try:
        for n, rows in pool_map(chunk, jobs, threads, on_exit=kill_lastz, warm=False):
            log.debug("chunk %d: %d calls", n, rows)
    except PoolError as exc:
        raise PoolError(f"{label}: {exc}") from exc
    finally:
        sweep_lastz_tmp()  # the folders of workers that died in this pool


# ---------------------------------------------------------------- elements, rows and files


def _index(elements: list[Element], order: dict[str, int]) -> tuple[dict[str, Element], dict[str, list[str]]]:
    """eid -> element, and eid -> its children's eids in A order (`order`: seqid -> fai index; start, -end).

    Fails on a duplicate eid or a seqid absent from A: each element gets exactly one row.
    """
    by: dict[str, Element] = {}
    for e in elements:
        if e.seqid not in order:
            raise ValueError(f"element {e.eid}: seqid {e.seqid!r} is not in the A FASTA")
        if e.eid in by:
            raise ValueError(f"duplicate element id {e.eid}")
        by[e.eid] = e
    children: dict[str, list[str]] = {}
    for e in sorted(elements, key=lambda e: (order[e.seqid], e.start, -e.end)):
        if e.parent is not None and e.parent in by:
            children.setdefault(e.parent, []).append(e.eid)
    return by, children


def _in_no_window(windows: list[Window], elements: dict[str, Element],
                  children: dict[str, list[str]]) -> list[Element]:
    """The elements no window reaches: neither top-level nor below an annotated ancestor."""
    seen: set[str] = set()
    todo = [eid for w in windows for eid in w.eids]
    while todo:
        eid = todo.pop()
        if eid not in seen:
            seen.add(eid)
            todo.extend(children.get(eid, ()))
    return [e for eid, e in elements.items() if eid not in seen]


def _no_window_line(e: Element) -> str:
    """U no_anchor for an element in no window: the call of an unresolved window holding only it, wid "."."""
    win = Window(-1, e.seqid, e.start - 1, e.end, None, 0, 0, 0, [e.eid], "no_anchor", 0, 0, math.nan)
    [call] = classify_window(win, None, "", "", {e.eid: e}, {})
    row = call_row(call, e, win)
    row["wid"] = None
    return call_line(row[c] for c in CALL_COLUMNS)


def _window_line(w: Window) -> str:
    """The window as a WINDOW_COLUMNS line: 1-based closed (the anchors: the last anchored A base on the left,
    the first on the right); B fields and anchors "." without a B locus."""
    b = w.b_seqid is not None
    return call_line((w.wid, w.a_seqid, w.a_start + 1, w.a_end, w.b_seqid, w.b_start + 1 if b else None,
                  w.b_end if b else None, {1: "+", -1: "-"}.get(w.strand), w.status,
                  w.left_anchor + 1 if b else None, w.right_anchor + 1 if b else None, w.margin, ",".join(w.eids)))


_SEQID, _START, _END, _LEVEL, _STATE_COL = map(CALL_COLUMNS.index, ("a_seqid", "a_start", "a_end", "level",
                                                                   "state"))
_SPLIT = max(_SEQID, _START, _END, _LEVEL, _STATE_COL) + 1  # a calls line is split only this far


def _sort_key(order: dict[str, int]) -> Callable[[str], tuple[int, int, int, int]]:
    """Key of a calls line: (A fai order, a_start, -a_end, level)."""
    def key(line: str) -> tuple[int, int, int, int]:
        f = line.split("\t", _SPLIT)
        return order[f[_SEQID]], int(f[_START]), -int(f[_END]), int(f[_LEVEL])
    return key


def _merge(chunks: list[Path], alone: list[str], out: Path, elements: dict[str, Element],
           order: dict[str, int], label: str, meta: dict) -> Counter:
    """Merge the sorted chunk files' lines and the no-window lines, as they are, into the call file `out`
    (write_call_lines, with `meta`); returns the count of each state.

    Every element must appear exactly once, or nothing is written and the chunks' markers are cleared (a
    rerun recomputes them).
    """
    states: Counter = Counter()

    def lines() -> Iterator[str]:
        seen: Counter = Counter()
        for line in heapq.merge(*map(call_lines, chunks), alone, key=_sort_key(order)):
            f = line.split("\t", _SPLIT)
            seen[f[0]] += 1
            states[f[_STATE_COL]] += 1
            yield line
        bad = {"element(s) missing from the calls": [x for x in elements if x not in seen],
               "element(s) called more than once": [x for x, n in seen.items() if n > 1],
               "call(s) of no annotated element": [x for x in seen if x not in elements]}
        problems = [f"{len(v)} {what}: {', '.join(v[:5])}{' ...' if len(v) > 5 else ''}"
                    for what, v in bad.items() if v]
        if problems:
            for path in chunks:
                clear_done(marker(path))
            raise RuntimeError(f"{label}: {'; '.join(problems)}; chunks cleared, they are recomputed on rerun")

    write_call_lines(out, lines(), meta)
    return states


def _sweep_tmp(path: Path) -> None:
    """Remove the temp files of `path` (checkpoint.atomic_path's .<name>.<pid>.<thread>.tmp) whose writer no
    longer runs: left by a run killed mid-write."""
    for tmp in path.parent.glob(f".{path.name}.*.tmp") if path.parent.is_dir() else ():
        pid = tmp.name[len(path.name) + 2:-len(".tmp")].split(".")[0]
        if pid.isdigit() and not pid_alive(int(pid)):
            tmp.unlink(missing_ok=True)
