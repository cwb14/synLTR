import contextlib
import os
import signal
import subprocess
import sys
import time

import pytest

from ltrtrace.base.log import LtrtraceError
from ltrtrace.base.runtime import PoolError, pool_map, shared


def square(x):
    return x * x


def test_pool_map_order_and_parent_first():
    items = list(range(12))
    assert pool_map(square, items, 4) == [square(x) for x in items]
    pids = pool_map(lambda _x: os.getpid(), items, 4)
    assert pids[0] == os.getpid()
    assert set(pids[1:]) != {os.getpid()}


def test_pool_map_cold_runs_nothing_in_the_parent():
    items = list(range(8))
    assert pool_map(square, items, 3, warm=False) == [square(x) for x in items]
    pids = pool_map(lambda _x: os.getpid(), items, 3, warm=False)
    assert os.getpid() not in pids
    assert pool_map(lambda _x: os.getpid(), [1, 2, 3], 1, warm=False) == [os.getpid()] * 3  # one thread: here


def test_pool_map_cold_exception_and_dead_worker():
    def fn(x):
        if x == 0:
            raise ValueError("boom 0")
        return x

    with pytest.raises(ValueError, match="boom 0"):
        pool_map(fn, [0, 1, 2], 3, warm=False)

    def die(x):
        if x == "die":
            os._exit(9)
        time.sleep(0.2)
        return x

    with pytest.raises(PoolError, match="die"):
        pool_map(die, ["a", "die", "c"], 2, warm=False)


def _slow_first(x):
    if x == 0:
        time.sleep(1.0)
    return x * x


def test_pool_map_progress_reports_each_item_as_it_finishes():
    """progress(i) runs in this process once per item as it finishes: in completion order with workers (the slow
    first item last), in item order when everything runs here, the warm item first; results stay in item order."""
    items, want = list(range(6)), [x * x for x in range(6)]
    seen = []
    assert pool_map(_slow_first, items, 3, warm=False, progress=seen.append) == want
    assert sorted(seen) == items and seen[-1] == 0
    seen.clear()
    assert pool_map(square, items, 1, progress=seen.append) == want and seen == items
    seen.clear()
    assert pool_map(square, items, 3, progress=seen.append) == want and seen[0] == 0 and sorted(seen) == items


def test_pool_map_serial_and_empty():
    assert pool_map(square, [], 4) == []
    assert pool_map(lambda _x: os.getpid(), [1, 2, 3], 1) == [os.getpid()] * 3


def test_pool_map_reads_shared_state():
    shared["offset"] = 100
    try:
        assert pool_map(lambda x: x + shared["offset"], [1, 2, 3, 4], 3) == [101, 102, 103, 104]
    finally:
        shared.clear()


def test_pool_map_exception_propagates():
    def fn(x):
        if x == 3:
            raise ValueError("boom 3")
        return x

    with pytest.raises(ValueError, match="boom 3"):
        pool_map(fn, [0, 1, 2, 3, 4], 3)


def test_pool_map_dead_worker():
    def fn(x):
        if x == "die":
            os._exit(9)
        time.sleep(0.2)
        return x

    start = time.monotonic()
    with pytest.raises(PoolError, match="die") as info:
        pool_map(fn, ["a", "b", "die", "c"], 2)
    assert time.monotonic() - start < 30
    assert isinstance(info.value, LtrtraceError)


_PARENT = r"""
import os, sys, time
from pathlib import Path
from ltrtrace.base import runtime

root = Path(sys.argv[1])
if sys.argv[2] == "poll":
    runtime._pdeathsig = lambda: False  # as without prctl (macOS): the workers poll their parent


def on_exit():
    (root / f"cleaned.{os.getpid()}").touch()


def work(x):
    if os.getpid() != parent:
        (root / f"busy.{os.getpid()}").touch()
        time.sleep(60)
    return x


parent = os.getpid()
runtime.pool_map(work, [0, 1, 2], 2, on_exit=on_exit)
"""


@pytest.mark.skipif(not sys.platform.startswith("linux"), reason="lists processes through /proc")
@pytest.mark.parametrize("guard", ["prctl", "poll"])
def test_pool_map_on_exit_runs_when_the_parent_dies(tmp_path, guard):
    """A worker whose parent is killed outright (SIGKILL, the OOM killer) runs `on_exit` before it exits: on Linux
    when prctl's SIGTERM arrives, elsewhere ("poll") when it sees its parent gone."""
    script = tmp_path / "parent.py"
    script.write_text(_PARENT)
    proc = subprocess.Popen([sys.executable, str(script), str(tmp_path), guard])
    try:
        t0 = time.monotonic()
        while len(list(tmp_path.glob("busy.*"))) < 2:
            assert proc.poll() is None and time.monotonic() - t0 < 30
            time.sleep(0.05)
        busy = sorted(int(p.name.split(".")[1]) for p in tmp_path.glob("busy.*"))
        proc.kill()
        proc.wait()
        t0 = time.monotonic()
        while len(list(tmp_path.glob("cleaned.*"))) < 2 and time.monotonic() - t0 < 10:
            time.sleep(0.05)
        assert sorted(int(p.name.split(".")[1]) for p in tmp_path.glob("cleaned.*")) == busy
    finally:
        if proc.poll() is None:
            proc.kill()
            proc.wait()
        for p in tmp_path.glob("busy.*"):
            with contextlib.suppress(ProcessLookupError):
                os.kill(int(p.name.split(".")[1]), signal.SIGKILL)


def test_pool_map_on_exit_runs_when_a_worker_is_terminated(tmp_path):
    """SIGTERM to a worker (the pool's own terminate) runs `on_exit`, then the worker exits with 128 + SIGTERM."""
    def fn(x):
        if x:
            os.kill(os.getpid(), signal.SIGTERM)
        return x

    with pytest.raises(PoolError):
        pool_map(fn, [0, 1], 2, on_exit=lambda: (tmp_path / f"cleaned.{os.getpid()}").touch())
    assert len(list(tmp_path.glob("cleaned.*"))) == 1


def test_pool_map_raises_at_the_first_failure(tmp_path):
    """A failing item stops the map at once, without waiting for slower items submitted before it: the items
    not yet handed to the workers are canceled (the executor queues up to workers + 1 ahead, so a few still
    run)."""
    def fn(x):
        if x == 1:
            time.sleep(2)
        elif x == 2:
            raise ValueError("boom 2")
        elif x > 2:
            (tmp_path / f"ran.{x}").touch()
            time.sleep(0.3)
        return x

    with pytest.raises(ValueError, match="boom 2"):
        pool_map(fn, list(range(20)), 2)
    assert len(list(tmp_path.glob("ran.*"))) <= 5  # waiting in submission order, all 17 later items would run


# a parent that maps slow items over two workers; argv: root folder, "group" (Ctrl-C reaches every process of the
# group) or "parent" (a SIGINT to this process alone)
_CTRL_C = r"""
import os, sys, time
from pathlib import Path
from ltrtrace.base import runtime

root = Path(sys.argv[1])


def work(x):
    try:
        (root / f"busy.{os.getpid()}.{x}").touch()
        time.sleep(60)
    finally:
        (root / f"cleaned.{os.getpid()}.{x}").touch()  # the item's own cleanup, as an atomic write's
    return x


try:
    runtime.pool_map(work, list(range(6)), 2, warm=False, on_exit=lambda: (root / f"on_exit.{os.getpid()}").touch())
except KeyboardInterrupt:
    print("interrupted", file=sys.stderr)
    sys.exit(130)
"""


@pytest.mark.skipif(not sys.platform.startswith("linux"), reason="lists processes through /proc")
@pytest.mark.parametrize("who", ["group", "parent"])
def test_pool_map_ctrl_c_ends_the_workers_quietly(tmp_path, who):
    """Ctrl-C while pool_map works (SIGINT to the whole process group, as a terminal sends it; or to this process
    alone, which passes it on to its workers): each busy worker unwinds its item, so the item's cleanup runs, then
    runs `on_exit` and ends without printing anything, and takes no other item; pool_map raises KeyboardInterrupt
    here at once."""
    (tmp_path / "parent.py").write_text(_CTRL_C)
    proc = subprocess.Popen([sys.executable, str(tmp_path / "parent.py"), str(tmp_path), who],
                            stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True, start_new_session=True)
    try:
        t0 = time.monotonic()
        while len(list(tmp_path.glob("busy.*"))) < 2:  # both workers inside an item
            assert proc.poll() is None and time.monotonic() - t0 < 30, proc.communicate()[1]
            time.sleep(0.05)
        busy = sorted(p.name.removeprefix("busy.") for p in tmp_path.glob("busy.*"))
        t0 = time.monotonic()
        if who == "group":
            os.killpg(proc.pid, signal.SIGINT)
        else:
            os.kill(proc.pid, signal.SIGINT)
        err = proc.communicate(timeout=30)[1]
        assert time.monotonic() - t0 < 10, "the workers finished their items instead of stopping"
    finally:
        if proc.poll() is None:
            os.killpg(proc.pid, signal.SIGKILL)
            proc.wait()
    assert proc.returncode == 130 and err == "interrupted\n", err
    assert sorted(p.name.removeprefix("cleaned.") for p in tmp_path.glob("cleaned.*")) == busy
    assert sorted(p.name.removeprefix("busy.") for p in tmp_path.glob("busy.*")) == busy  # no other item started
    assert len(list(tmp_path.glob("on_exit.*"))) == 2


# a parent that ignores SIGINT, as a script's background job (`ltrtrace ... &`) starts with it ignored, maps short
# items over two workers and prints the results; argv: root folder
_IGNORED = r"""
import os, signal, sys, time
from pathlib import Path
from ltrtrace.base import runtime

signal.signal(signal.SIGINT, signal.SIG_IGN)
root = Path(sys.argv[1])


def work(x):
    (root / f"busy.{os.getpid()}.{x}").touch()
    time.sleep(0.5)
    (root / f"done.{x}").touch()
    return x * x


print(runtime.pool_map(work, list(range(6)), 2, warm=False))
"""


@pytest.mark.parametrize("who", ["group", "parent"])
def test_pool_map_keeps_an_inherited_sigint_ignore(tmp_path, who):
    """A run started with SIGINT ignored (a script's background job, `ltrtrace ... &`) ignores Ctrl-C in its workers
    too: SIGINT to the whole process group, or to this process alone, while the workers hold items changes nothing,
    and the map finishes with every item's output and exit 0."""
    (tmp_path / "parent.py").write_text(_IGNORED)
    proc = subprocess.Popen([sys.executable, str(tmp_path / "parent.py"), str(tmp_path)], stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, text=True, start_new_session=True)
    try:
        t0 = time.monotonic()
        while len(list(tmp_path.glob("busy.*"))) < 2:  # both workers inside an item
            assert proc.poll() is None and time.monotonic() - t0 < 30, proc.communicate()[1]
            time.sleep(0.05)
        if who == "group":
            os.killpg(proc.pid, signal.SIGINT)
        else:
            os.kill(proc.pid, signal.SIGINT)
        out, err = proc.communicate(timeout=30)
    finally:
        if proc.poll() is None:
            os.killpg(proc.pid, signal.SIGKILL)
            proc.wait()
    assert proc.returncode == 0 and err == "", err
    assert out == f"{[x * x for x in range(6)]}\n"
    assert sorted(p.name for p in tmp_path.glob("done.*")) == [f"done.{x}" for x in range(6)]
