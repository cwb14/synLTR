"""The one process pool (fork, ordered results, workers that die with the parent, a named dead worker, workers that
stop quietly on Ctrl-C) and the sizing of numba's parallel work."""

import contextlib
import ctypes
import multiprocessing
import os
import signal
import sys
import threading
import time
from collections.abc import Callable, Sequence
from concurrent.futures import FIRST_COMPLETED, FIRST_EXCEPTION, ProcessPoolExecutor, wait
from concurrent.futures.process import BrokenProcessPool

from ltrtrace.base import constants
from ltrtrace.base.log import LtrtraceError

shared: dict = {}  # state that forked workers read; fill it before calling pool_map

_job: dict = {}  # fn, items and the per-item worker slots, inherited by the forked workers


class PoolError(LtrtraceError):
    """A worker process died (killed by a signal, out of memory, crashed) or could not start."""


def set_numba_threads(n: int) -> None:
    """Use `n` numba threads, capped at the number the process started with (NUMBA_NUM_THREADS)."""
    import numba

    numba.set_num_threads(max(1, min(n, numba.config.NUMBA_NUM_THREADS)))


def n_chunks(C: int, n_items: int) -> int:
    """Chunks of a numba parallel loop that accumulates into one C x C matrix per chunk: PAR_CHUNKS (fixed, so results do
    not depend on -t), fewer when that many float64 matrices would exceed PAR_CHUNK_BYTES, never more than the items,
    at least one."""
    return int(max(1, min(constants.PAR_CHUNKS, n_items, constants.PAR_CHUNK_BYTES // max(8 * C * C, 1))))


def pool_map[T, R](fn: Callable[[T], R], items: Sequence[T], threads: int, *,
                   on_exit: Callable[[], None] | None = None, warm: bool = True,
                   progress: Callable[[int], None] | None = None) -> list[R]:
    """Return [fn(x) for x in items], using up to `threads` forked workers; results are in input order.

    With `warm` (the default) the first item runs in this process before any fork, so `fn`'s numba kernels compile (or
    load from numba's cache) once, here, and the workers inherit them. With warm=False every item goes to the workers
    at once (none runs here): for an `fn` that runs no numba kernel, or whose kernels each worker then compiles or
    loads itself once this process has run a numba kernel, so that numba's one-time setup, which reads the installed
    packages' metadata through file handles a fork shares, has run here and never runs in several workers at once.
    With one thread or one item everything runs here.

    An exception raised by `fn` propagates unchanged, as soon as it happens: the items no worker has started are
    canceled, the running ones finish first. A worker that dies raises PoolError naming the item(s) it was running.
    `fn` and the items are inherited through the fork, not pickled; only the results travel back, so they must be
    picklable. `on_exit` runs in a worker that is terminated (SIGTERM: the pool's own terminate, or the parent's death
    on Linux) or that finds its parent gone, just before it exits: it ends what the worker started, such as a child
    process. `progress(i)` runs in this process once per item as it finishes, with its index (items that finish
    together in index order): a caller logs progress with it.

    Ctrl-C (SIGINT, which a terminal sends to every process of the run) stops the map: each worker unwinds the item it
    holds, so the item's own cleanup runs (a temporary file removed), then runs `on_exit` and ends without printing
    anything or taking another item; a SIGINT sent to this process alone is passed on to the workers. KeyboardInterrupt
    then propagates here once the workers are gone. A process started with SIGINT ignored (a script's background job,
    `ltrtrace ... &`) keeps ignoring it, workers included, and the map runs to its end.
    """
    if not items:
        return []
    results: list = [None] * len(items)
    if threads <= 1 or len(items) == 1:
        for i, item in enumerate(items):
            results[i] = fn(item)
            if progress is not None:
                progress(i)
        return results
    first = 1 if warm else 0
    if warm:
        results[0] = fn(items[0])
        if progress is not None:
            progress(0)
    running = multiprocessing.get_context("fork").RawArray("l", len(items))  # pid of the worker per item
    _job.update(fn=fn, items=items, running=running, on_exit=on_exit)
    pool = ProcessPoolExecutor(min(threads, len(items) - first), mp_context=multiprocessing.get_context("fork"),
                               initializer=_init_worker,
                               initargs=(os.getpid(), on_exit, signal.getsignal(signal.SIGINT) is signal.SIG_IGN))
    try:
        futures = {pool.submit(_run_item, i): i for i in range(first, len(items))}
        pending = set(futures)
        while pending:  # the first failure stops the map; the finally block cancels what has not started
            done, pending = wait(pending, return_when=FIRST_EXCEPTION if progress is None else FIRST_COMPLETED)
            for fut in sorted(done, key=lambda f: (f.exception() is None, futures[f])):  # a failure first
                i = futures[fut]
                try:
                    results[i] = fut.result()
                except BrokenProcessPool as exc:
                    lost = [j for j in futures.values() if running[j] > 0 and not _done(futures, j)]
                    names = ", ".join(f"{j} ({_short(items[j])})" for j in lost) or f"{i} ({_short(items[i])})"
                    raise PoolError(f"a worker process died or could not start (killed by a signal or out of "
                                    f"memory, or crashed); item {names}") from exc
                if progress is not None:
                    progress(i)
    except KeyboardInterrupt:
        for worker in multiprocessing.active_children():  # a SIGINT to this process alone reaches them too
            with contextlib.suppress(ProcessLookupError):
                os.kill(worker.pid, signal.SIGINT)
        raise
    finally:
        pool.shutdown(wait=True, cancel_futures=True)
        _job.clear()
    return results


def _done(futures: dict, i: int) -> bool:
    """True if item `i`'s future finished without an error."""
    fut = next(f for f, j in futures.items() if j == i)
    return fut.done() and not fut.cancelled() and fut.exception() is None


def _short(item) -> str:
    """The item's repr, cut to 80 characters."""
    text = repr(item)
    return text if len(text) <= 80 else text[:77] + "..."


def _run_item(i: int):
    """Run item `i` in a worker, recording which worker has it. An item interrupted (SystemExit, raised by
    _interrupted) has unwound when it reaches here: the worker then ends, after `on_exit`, without taking another."""
    _job["running"][i] = os.getpid()
    try:
        return _job["fn"](_job["items"][i])
    except SystemExit as exc:
        _exit(exc.code if isinstance(exc.code, int) else 1, _job.get("on_exit"))


def _interrupted(signum: int, frame) -> None:
    """A worker's SIGINT (Ctrl-C): SystemExit, which unwinds the item the worker holds (_run_item then ends the
    worker) or ends an idle worker without a word; a second SIGINT is ignored, so it cannot cut that cleanup short."""
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    raise SystemExit(128 + signum)


def _init_worker(parent: int, on_exit: Callable[[], None] | None = None, sigint_ignored: bool = False) -> None:
    """Make this worker exit when the parent dies; SIGTERM from the pool's shutdown also ends it. `on_exit` runs
    first either way (pool_map). SIGINT stops it quietly (_interrupted), unless the parent ignores SIGINT
    (`sigint_ignored`): then the worker keeps the ignore it inherited."""
    if not sigint_ignored:
        signal.signal(signal.SIGINT, _interrupted)
    if on_exit is not None:
        signal.signal(signal.SIGTERM, lambda signum, frame: _exit(128 + signum, on_exit))
    if not _pdeathsig():
        threading.Thread(target=_watch_parent, args=(parent, on_exit), name="watch-parent", daemon=True).start()
    if os.getppid() != parent:  # the parent died before this ran
        _exit(1, on_exit)


def _exit(status: int, on_exit: Callable[[], None] | None) -> None:
    """Run `on_exit` (if any), then end this process at once with `status`."""
    try:
        if on_exit is not None:
            on_exit()
    finally:
        os._exit(status)


def _pdeathsig() -> bool:
    """Ask the kernel for SIGTERM when this process's parent dies (Linux prctl); False where it cannot."""
    if not sys.platform.startswith("linux"):
        return False
    try:
        prctl = ctypes.CDLL(None, use_errno=True).prctl
    except (OSError, AttributeError):
        return False
    zero = ctypes.c_ulong(0)
    return prctl(ctypes.c_int(constants.PR_SET_PDEATHSIG), ctypes.c_ulong(signal.SIGTERM),
                 zero, zero, zero) == 0


def _watch_parent(parent: int, on_exit: Callable[[], None] | None = None) -> None:
    """Exit this process (after `on_exit`) once its parent is no longer `parent`."""
    while os.getppid() == parent:
        time.sleep(constants.PARENT_POLL_SECONDS)
    _exit(1, on_exit)
