"""Pair genotyping (routes.assembly.pair): the toy pair (builders.build_pair) end to end through anchoring,
windows, LASTZ window alignment and classification in worker processes, with chunk checkpoints and resume."""
from __future__ import annotations

import contextlib
import logging
import math
import os
import random
import signal
import subprocess
import sys
import time
from collections import Counter
from dataclasses import replace
from pathlib import Path

import builders as toy
import mappy
import pytest
from builders import pair_context, read_calls, read_rows
from markers import needs_lastz

import synltr.routes.assembly.pair as genotyping
from synltr.base import runtime
from synltr.base.runtime import PoolError
from synltr.io.samples import Sample
from synltr.io.seqio import Fasta, revcomp
from synltr.routes.assembly.anchor import FlankMm2Route, PairContext
from synltr.routes.assembly.blocks import load_blocks, save_blocks
from synltr.routes.assembly.pair import genotype_pair

pytestmark = needs_lastz

SEED = 7
SEEDS = (7, 0, 4, 6, 9, 11, 14)  # test_classify's seeds: chance microhomology at the junctions
NAMES = [n for n in toy.CASES if n != "U_flank"]  # build_pair's default cases, in its order
OLD_NS = 10**18  # an mtime (2001-09-09) that no file written by a test run has


def run_pair(d, a_seq, b_seq, elements, threads=2):
    with pair_context(d, a_seq, b_seq, elements, threads) as ctx:
        return genotype_pair(ctx, FlankMm2Route())


def marker(path):
    return path.with_name(path.name + ".done.json")


@contextlib.contextmanager
def deadline(seconds):
    """Fail instead of hanging: TimeoutError if the block still runs after `seconds`."""
    def expire(signum, frame):
        raise TimeoutError(f"still running after {seconds} s")

    old = signal.signal(signal.SIGALRM, expire)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)


def flags(row):
    return set() if row["flags"] == "." else set(row["flags"].split(","))


def assert_states(calls, pair, eids=None):
    names = {e.eid: e.name for e in pair.elements}
    for eid in pair.expect if eids is None else eids:
        assert calls[eid]["state"] == pair.expect[eid], (names[eid], calls[eid])


@pytest.fixture(scope="module")
def pair():
    return toy.build_pair(random.Random(SEED))


@pytest.fixture(scope="module")
def base(pair, tmp_path_factory):
    """The toy pair genotyped with 2 worker processes: its pair file."""
    return run_pair(tmp_path_factory.mktemp("base"), pair.a, pair.b, pair.elements)


@pytest.mark.parametrize("seed", SEEDS)
def test_pair_end_to_end(seed, pair, base, tmp_path):
    if seed != SEED:
        pair = toy.build_pair(random.Random(seed))
        base = run_pair(tmp_path, pair.a, pair.b, pair.elements)
    calls = read_calls(base)
    assert set(calls) == set(pair.expect)
    assert_states(calls, pair)  # exact, every case
    for eid, want in pair.flags.items():
        assert want <= flags(calls[eid]), calls[eid]
    for eid, unwanted in pair.noflags.items():
        assert not unwanted & flags(calls[eid]), calls[eid]
    for eid, tier in pair.tiers.items():
        assert calls[eid]["tier"] == tier, calls[eid]
    for eid, reason in pair.reasons.items():
        assert calls[eid]["reason"] == reason, calls[eid]


def test_pair_file_has_each_element_once_and_threads_do_not_matter(pair, base, tmp_path):
    rows = read_rows(base)
    assert Counter(r["eid"] for r in rows) == Counter(e.eid for e in pair.elements)  # each exactly once
    keys = [(int(r["a_start"]), -int(r["a_end"]), int(r["level"])) for r in rows]
    assert keys == sorted(keys)  # A order (one seqid)
    out = base.parent.parent
    work = out / "work"
    assert base == out / "pairs" / "A__B.calls.tsv.gz" and marker(base).is_file()
    windows = read_rows(work / "windows" / "A__B.tsv.gz", genotyping.WINDOW_COLUMNS)
    assert [int(w["wid"]) for w in windows] == list(range(len(windows)))
    n = math.ceil(len(windows) / genotyping._chunk_size(len(windows)))  # the toy's 20 windows: 3 chunks of <= 8
    chunks = [work / "calls" / "A__B" / f"chunk_{i}.tsv.gz" for i in range(n)]
    assert sorted(chunks[0].parent.iterdir()) == sorted(chunks + [marker(p) for p in chunks])  # no temp left
    assert sorted(p.name for p in base.parent.iterdir()) == [base.name, marker(base).name]
    top = [e for e in pair.elements if e.level == 0]
    assert sorted(x for w in windows for x in w["eids"].split(",")) == sorted(e.eid for e in top)
    e = next(e for e in top if e.name == "E.e")
    w = next(w for w in windows if w["eids"] == e.eid)
    assert (w["status"], w["strand"], w["a_seqid"]) == ("ok", "+", "chr1")
    assert int(w["a_start"]) < e.start and e.end < int(w["a_end"])  # 1-based inclusive, around the element
    for p in (base, *chunks, work / "windows" / "A__B.tsv.gz"):
        head = p.read_bytes()[:8]
        assert head[3] == 0 and head[4:8] == bytes(4)  # gzip without a file name and with mtime 0

    one = run_pair(tmp_path, pair.a, pair.b, pair.elements, threads=1)  # in-process, no worker pool
    assert one.read_bytes() == base.read_bytes()


@pytest.mark.parametrize("n, size", [(0, 8), (1, 8), (100, 8), (5_600, 22), (60_000, 235), (10**6, 2000)])
def test_chunk_size(n, size):
    # max(CHUNK_MIN, min(CHUNK_WINDOWS, ceil(n / CHUNK_TARGET))): about 256 chunks of 8 to 2000 windows
    assert genotyping._chunk_size(n) == size


def test_calls_do_not_depend_on_threads_or_chunks(pair, base, tmp_path, monkeypatch):
    # 1 and 3 threads (3 chunks of the toy's 20 windows), then chunks of 1 window and of all 20: the same calls
    one = run_pair(tmp_path / "t1", pair.a, pair.b, pair.elements, threads=1)
    three = run_pair(tmp_path / "t3", pair.a, pair.b, pair.elements, threads=3)
    assert one.read_bytes() == three.read_bytes() == base.read_bytes()
    for size in (1, 2000):
        monkeypatch.setattr(genotyping, "CHUNK_MIN", size)
        calls = run_pair(tmp_path / f"c{size}", pair.a, pair.b, pair.elements, threads=3)
        assert calls.read_bytes() == base.read_bytes()


def test_resume_at_another_thread_count_reuses_chunks(pair, base, tmp_path):
    # the chunk size does not depend on --threads: a run resumed with 1 thread finds a 3-thread run's chunks
    route = FlankMm2Route()
    with pair_context(tmp_path, pair.a, pair.b, pair.elements, threads=3) as ctx:
        path = genotype_pair(ctx, route)
        chunks = sorted((ctx.workdir / "calls" / "A__B").glob("chunk_*.tsv.gz"))
        assert len(chunks) == 3
        for p in chunks:
            os.utime(p, ns=(OLD_NS, OLD_NS))
        for p in (path, marker(path)):  # interrupted before the merge
            p.unlink()
        ctx.threads = 1
        assert genotype_pair(ctx, route).read_bytes() == base.read_bytes()
        assert all(p.stat().st_mtime_ns == OLD_NS for p in chunks)  # reused, not recomputed


def test_pair_lowercase_target(pair, base, tmp_path):
    # a soft-masked (lowercase) B: the same calls, byte for byte
    assert run_pair(tmp_path, pair.a, pair.b.lower(), pair.elements).read_bytes() == base.read_bytes()


def test_pair_minus_strand_target(pair, base, tmp_path):
    # B holds the toy on its minus strand: each B window is reverse-complemented to A's orientation before
    # it is aligned, so the states stay; the B sites are those of the plus-strand run, flipped
    b = revcomp(pair.b)
    calls, fwd = read_calls(run_pair(tmp_path, pair.a, b, pair.elements)), read_calls(base)
    assert_states(calls, pair)
    sites = 0
    for eid, r in calls.items():
        if r["b_seqid"] == ".":
            continue
        sites += 1
        f = fwd[eid]
        assert r["b_strand"] == "-" and "inv" in flags(r) and f["b_strand"] == "+", (r, f)
        assert (int(r["b_start"]), int(r["b_end"])) == (len(b) + 1 - int(f["b_end"]),
                                                       len(b) + 1 - int(f["b_start"])), (r, f)
    assert sites >= 15


def _with_decoy(pair):
    """B with a copy of the E case's A (L T e T R: its flanks with the element inside) pasted into the middle
    of the spacer nearest to 100 kb past the E site: a paralog of the site, out of collinear order.

    Returns (B, decoy start, the E site's B interval)."""
    rng = random.Random(SEED)
    cases = [toy.case(n, rng) for n in NAMES]  # build_pair draws the cases first, in order
    starts, ob = {}, toy.SPACER
    for n, c in zip(NAMES, cases):
        starts[n] = ob
        ob += len(c.b) + toy.SPACER
    ec = cases[NAMES.index("E")]
    site = (starts["E"], starts["E"] + len(ec.b))
    assert pair.b[site[0]:site[1]] == ec.b
    mids = [starts[n] + len(c.b) + toy.SPACER // 2 for n, c in zip(NAMES, cases)]
    at = min(mids, key=lambda m: abs(m - site[0] - 100_000))
    return pair.b[:at] + ec.a + pair.b[at:], at, site, ec.a


def test_pair_repetitive_flank_decoy(pair, tmp_path):
    b, at, site, region = _with_decoy(pair)
    assert abs(at - site[0] - 100_000) < 20_000
    # the trap is real: the E case's A region aligns best, in full, to the decoy
    best = next(h for h in mappy.Aligner(seq=b, preset="asm20").map(region) if h.is_primary)
    assert best.r_st == at and best.q_en - best.q_st == len(region) and best.mlen == len(region)
    calls = read_calls(run_pair(tmp_path, pair.a, b, pair.elements))
    assert_states(calls, pair)  # E.e stays E (the collinear ortholog), not P (the paralog)
    e = next(e for e in pair.elements if e.name == "E.e")
    assert site[0] < int(calls[e.eid]["b_start"]) <= int(calls[e.eid]["b_end"]) <= site[1]


def test_pair_unrelated_genomes(pair, tmp_path, caplog):
    b = toy.rand_seq(300_000, random.Random(99))
    with caplog.at_level(logging.INFO, logger="synltr"):
        calls = read_calls(run_pair(tmp_path, pair.a, b, pair.elements))
    assert set(calls) == set(pair.expect)
    assert {(r["state"], r["reason"]) for r in calls.values()} == {("U", "no_anchor")}
    warn = [r for r in caplog.records if r.levelno == logging.WARNING and r.name == "synltr.routes.assembly.pair"]
    assert len(warn) == 1 and "no orthologous alignment between A and B?" in warn[0].getMessage()


def test_resume(pair, base, tmp_path, monkeypatch, caplog):
    monkeypatch.setattr(genotyping, "CHUNK_MIN", 6)  # the toy's 20 windows in 4 chunks
    route = FlankMm2Route()
    with pair_context(tmp_path, pair.a, pair.b, pair.elements) as ctx:
        path = genotype_pair(ctx, route)
        assert path.read_bytes() == base.read_bytes()  # chunking changes nothing
        cdir = ctx.workdir / "calls" / "A__B"
        chunks = [cdir / f"chunk_{n}.tsv.gz" for n in range(4)]
        assert sorted(cdir.iterdir()) == sorted(chunks + [marker(p) for p in chunks])
        before = {p: (p.read_bytes(), marker(p).read_bytes()) for p in chunks}
        for p in chunks:
            os.utime(p, ns=(OLD_NS, OLD_NS))
        # a run interrupted before chunk_0 and the pair file were finished
        for p in (chunks[0], marker(chunks[0]), path, marker(path)):
            p.unlink()
        assert genotype_pair(ctx, route) == path
        assert path.read_bytes() == base.read_bytes()
        assert chunks[0].stat().st_mtime_ns != OLD_NS and chunks[0].read_bytes() == before[chunks[0]][0]
        for p in chunks[1:]:  # only chunk_0 was recomputed
            assert p.stat().st_mtime_ns == OLD_NS and (p.read_bytes(), marker(p).read_bytes()) == before[p]

        # interrupted after the last chunk, before the pair file: only the merge runs again
        for p in chunks:
            os.utime(p, ns=(OLD_NS, OLD_NS))
        for p in (path, marker(path)):
            p.unlink()
        assert genotype_pair(ctx, route).read_bytes() == base.read_bytes()
        assert all(p.stat().st_mtime_ns == OLD_NS for p in chunks)

        # a finished pair is skipped entirely
        windows = ctx.workdir / "windows" / "A__B.tsv.gz"
        kept = [path, marker(path), windows, *chunks]
        for p in kept:
            os.utime(p, ns=(OLD_NS, OLD_NS))
        with caplog.at_level(logging.INFO, logger="synltr.routes.assembly.pair"):
            assert genotype_pair(ctx, route) == path
        assert "skip" in caplog.text
        assert all(p.stat().st_mtime_ns == OLD_NS for p in kept)


def test_calls_key_follows_the_aligner_not_the_anchor_file_time(pair, base, tmp_path, monkeypatch, caplog):
    """Anchors saved again (as after deleting work/) are the same bytes: an .npz's zip members carry no write
    time, so the calls' key, which holds the anchors' checksum, keeps the calls. Another aligner version
    recomputes them."""
    route = FlankMm2Route()
    with pair_context(tmp_path, pair.a, pair.b, pair.elements) as ctx:
        path = genotype_pair(ctx, route)
        os.utime(path, ns=(OLD_NS, OLD_NS))
        npz = ctx.workdir / "anchors" / "A__B.flank-mm2.npz"
        before = npz.read_bytes()
        monkeypatch.setattr(time, "time", lambda: 4e9)  # saved again, in 2096
        save_blocks(npz, *load_blocks(npz))
        monkeypatch.undo()
        assert npz.read_bytes() == before
        with caplog.at_level(logging.INFO, logger="synltr"):
            assert genotype_pair(ctx, route) == path
        assert "calls done, skip" in caplog.text and path.stat().st_mtime_ns == OLD_NS
        monkeypatch.setattr(genotyping, "lastz_version", lambda: "lastz 9.99")
        assert genotype_pair(ctx, route).read_bytes() == base.read_bytes()
        assert path.stat().st_mtime_ns != OLD_NS  # recomputed


def test_rerun_after_annotation_change(pair, tmp_path):
    route = FlankMm2Route()
    by = {e.name: e for e in pair.elements}
    with pair_context(tmp_path, pair.a, pair.b, pair.elements) as ctx:
        genotype_pair(ctx, route)
        # a field the windows do not hold (the anchors do not read elements): recomputed, not reused
        e = by["E.e"]
        ctx.a_elements = [replace(x, tsd="GATTA") if x is e else x for x in pair.elements]
        calls = read_calls(genotype_pair(ctx, route))
        assert calls[e.eid]["tsd_a"] == "GATTA" and calls[e.eid]["state"] == "E"
        # the host removed: its child, in no window, is U no_anchor
        h, c = by["host_E.h"], by["host_E.c"]
        ctx.a_elements = [x for x in pair.elements if x is not h]
        calls = read_calls(genotype_pair(ctx, route))
        assert set(calls) == {x.eid for x in ctx.a_elements}
        assert [calls[c.eid][k] for k in ("state", "tier", "reason", "wid", "b_seqid")] == ["U", ".", "no_anchor",
                                                                                       ".", "."]
        assert_states(calls, pair, [x.eid for x in ctx.a_elements if x is not c])


FAULTS = {  # fault -> the error genotype_pair must raise, and its message
    "killed": (PoolError, r"A vs B: a worker process died or could not start.*chunk 1"),  # SIGKILL (OOM killer)
    "init_raises": (PoolError, r"A vs B: a worker process died or could not start"),  # the pool initializer
    "fasta_fails": (RuntimeError, r"A vs B: chunk \d+ failed: OSError: cannot open the FASTA"),  # a worker's open
    "raises": (RuntimeError, r"A vs B: chunk 1 failed: ValueError: boom"),  # an error inside a chunk
}


@pytest.mark.parametrize("fault", list(FAULTS))
def test_worker_failure_fails_fast_and_resumes(pair, base, tmp_path, monkeypatch, fault):
    """A failing or dying worker makes genotype_pair raise promptly (never hang), naming the pair; the failed
    chunk is never marked done; a rerun reuses the finished chunks and writes the same pair file."""
    monkeypatch.setattr(genotyping, "CHUNK_MIN", 6)  # the toy's 20 windows in 4 chunks; wid 8 in chunk 1
    parent, real_rows = os.getpid(), genotyping._Genotyper.rows

    def rows(self, win):
        if win.wid == 8:
            assert os.getpid() != parent  # in a worker, never in the test process
            if fault == "killed":
                os.kill(os.getpid(), signal.SIGKILL)
            raise ValueError("boom")
        return real_rows(self, win)

    def fail(*args, **kwargs):
        raise OSError("cannot open the FASTA")

    route = FlankMm2Route()
    with pair_context(tmp_path, pair.a, pair.b, pair.elements) as ctx:
        with pytest.MonkeyPatch.context() as mp:
            if fault in ("killed", "raises"):
                mp.setattr(genotyping._Genotyper, "rows", rows)
            else:
                mp.setattr(*((runtime, "_init_worker") if fault == "init_raises" else (genotyping, "Fasta")), fail)
            error, message = FAULTS[fault]
            with deadline(30), pytest.raises(error, match=message) as caught:
                genotype_pair(ctx, route)
            assert type(caught.value) is error
        out = ctx.outdir / "pairs" / "A__B.calls.tsv.gz"
        assert not out.exists() and not marker(out).exists()
        cdir = ctx.workdir / "calls" / "A__B"
        chunks = [cdir / f"chunk_{n}.tsv.gz" for n in range(4)]
        done = [p for p in chunks if marker(p).exists()]
        if fault in ("killed", "raises"):
            assert chunks[1] not in done
        for p in done:
            os.utime(p, ns=(OLD_NS, OLD_NS))
        with deadline(30):
            assert genotype_pair(ctx, route).read_bytes() == base.read_bytes()
        assert all(p.stat().st_mtime_ns == OLD_NS for p in done)  # finished chunks reused
        assert sorted(cdir.iterdir()) == sorted(chunks + [marker(p) for p in chunks])  # no temp files left


# the genotyping parent for test_workers_exit_when_the_parent_dies: argv = tests dir, guard, run dir
_PARENT = r"""
import os, random, sys, time
from pathlib import Path

tests, guard, root = sys.argv[1], sys.argv[2], Path(sys.argv[3])
sys.path.insert(0, tests)
import builders
import synltr.routes.assembly.pair as genotyping
from synltr.base import runtime
from synltr.io.samples import Sample
from synltr.io.seqio import Fasta
from synltr.routes.assembly.anchor import FlankMm2Route, PairContext

if guard == "poll":
    runtime._pdeathsig = lambda: False  # as without prctl (macOS): the workers poll their parent
genotyping.CHUNK_MIN = 1  # one window per chunk
pair = builders.build_pair(random.Random(7), ["E", "P", "T"])
parent, real_rows = os.getpid(), genotyping._Genotyper.rows


def rows(self, win):
    if os.getpid() != parent:  # the first chunk a worker takes runs slowly; the other worker goes idle
        try:
            os.close(os.open(root / "slow", os.O_CREAT | os.O_EXCL))
        except FileExistsError:
            return real_rows(self, win)
        time.sleep(120)
    return real_rows(self, win)


genotyping._Genotyper.rows = rows
root.mkdir(parents=True, exist_ok=True)
fastas = []
for name, seq in (("A", pair.a), ("B", pair.b)):
    fastas.append(root / f"{name}.fa")
    builders.write_fasta(fastas[-1], {"chr1": seq})
with Fasta(fastas[0]) as fa, Fasta(fastas[1]) as fb:
    ctx = PairContext(Sample("A", fastas[0], None, None, None), Sample("B", fastas[1], None, None, None), fa, fb,
                      pair.elements, None, root / "out", 2)
    genotyping.genotype_pair(ctx, FlankMm2Route())
"""

def _stat(pid):
    """The fields of /proc/<pid>/stat after the command name (state first), or None without that process."""
    try:
        text = Path(f"/proc/{pid}/stat").read_text()
    except OSError:
        return None
    return text[text.rindex(")") + 2:].split()


def _children(pid):
    """(pid, start time) of the live (non-zombie) children of process `pid`."""
    out = []
    for d in Path("/proc").iterdir():
        st = _stat(d.name) if d.name.isdigit() else None
        if st is not None and st[0] != "Z" and int(st[1]) == pid:
            out.append((int(d.name), st[19]))
    return out


def _alive(pid, start):
    st = _stat(pid)
    return st is not None and st[0] != "Z" and st[19] == start  # the same process, not a reused pid


@pytest.mark.skipif(not sys.platform.startswith("linux"), reason="lists processes through /proc")
@pytest.mark.parametrize("guard", ["prctl", "poll"])
def test_workers_exit_when_the_parent_dies(tmp_path, guard):
    """A genotyping parent killed mid-run (SIGKILL, as by the OOM killer) takes its workers with it, the busy
    one and the idle one: on Linux by prctl(PR_SET_PDEATHSIG), elsewhere ("poll": prctl off) by polling."""
    (tmp_path / "parent.py").write_text(_PARENT)
    root, err = tmp_path / "run", tmp_path / "parent.err"
    with open(err, "w") as fh:
        proc = subprocess.Popen([sys.executable, str(tmp_path / "parent.py"), str(Path(__file__).parents[2]), guard,
                                 str(root)], stdout=subprocess.DEVNULL, stderr=fh)
    workers = []
    try:
        with deadline(30):
            while not (root / "slow").exists():  # a worker is inside the slow chunk
                assert proc.poll() is None, err.read_text()
                time.sleep(0.05)
            workers = _children(proc.pid)
            assert len(workers) == 2, workers
            proc.kill()
            proc.wait()
            t0 = time.monotonic()
            while any(_alive(*w) for w in workers) and time.monotonic() - t0 < 5:
                time.sleep(0.05)
        assert not [w for w in workers if _alive(*w)], f"workers alive 5 s after the parent died: {workers}"
    finally:
        if proc.poll() is None:
            proc.kill()
            proc.wait()
        for pid, start in workers:
            if _alive(pid, start):
                os.kill(pid, signal.SIGKILL)


def test_bad_elements_fail_before_anchoring(tmp_path):
    small = toy.build_pair(random.Random(3), ["E"])
    [e] = small.elements
    for bad, msg in (([e, e], "duplicate element id"), ([replace(e, seqid="chr9")], "not in the A FASTA")):
        with pair_context(tmp_path, small.a, small.b, bad) as ctx:
            with pytest.raises(ValueError, match=msg):
                genotype_pair(ctx, FlankMm2Route())
            assert not ctx.workdir.exists()  # nothing anchored


@pytest.mark.parametrize("fault", ["missing", "duplicate"])
def test_pair_file_needs_each_element_once(tmp_path, monkeypatch, fault):
    small = toy.build_pair(random.Random(3), ["E", "host_E"])
    real = genotyping.classify_window

    def broken(*args):
        calls = real(*args)
        return calls[:-1] if fault == "missing" else calls + calls[:1]

    monkeypatch.setattr(genotyping, "classify_window", broken)
    route = FlankMm2Route()
    with pair_context(tmp_path, small.a, small.b, small.elements, threads=1) as ctx:
        with pytest.raises(RuntimeError, match="missing" if fault == "missing" else "more than once"):
            genotype_pair(ctx, route)
        out = ctx.outdir / "pairs" / "A__B.calls.tsv.gz"
        assert not out.exists() and not marker(out).exists()
        monkeypatch.undo()  # fixed: the chunks were cleared, so a rerun recomputes them
        calls = read_calls(genotype_pair(ctx, route))
    assert {eid: r["state"] for eid, r in calls.items()} == small.expect


def test_stale_chunks_and_temp_files_are_swept(pair, base, tmp_path, monkeypatch):
    """A rerun removes the chunk files an earlier run with more chunks left, and the temp
    files of the pair file and the windows file whose writer is dead (a live writer's stay)."""
    route = FlankMm2Route()
    with pair_context(tmp_path, pair.a, pair.b, pair.elements) as ctx:
        monkeypatch.setattr(genotyping, "CHUNK_MIN", 1)  # one window per chunk
        out = genotype_pair(ctx, route)
        cdir = ctx.workdir / "calls" / "A__B"
        assert len(list(cdir.glob("chunk_*.tsv.gz"))) > 3
        dead = subprocess.Popen(["true"])
        dead.wait()
        stale = [out.with_name(f".{out.name}.{dead.pid}.tmp"),
                 ctx.workdir / "windows" / f".A__B.tsv.gz.{dead.pid}.tmp"]
        live = out.with_name(f".{out.name}.{os.getppid()}.tmp")  # a running process's
        for p in [*stale, live]:
            p.write_text("partial")
        for p in (out, marker(out)):  # rerun with chunks of 3 windows: only the merge's inputs change
            p.unlink()
        monkeypatch.setattr(genotyping, "CHUNK_MIN", 3)
        assert genotype_pair(ctx, route).read_bytes() == base.read_bytes()
        n = math.ceil(len(read_rows(ctx.workdir / "windows" / "A__B.tsv.gz", genotyping.WINDOW_COLUMNS)) / 3)
        assert sorted(p.name for p in cdir.iterdir()) == sorted(
            x for i in range(n) for x in (f"chunk_{i}.tsv.gz", f"chunk_{i}.tsv.gz.done.json"))
        assert not any(p.exists() for p in stale) and live.exists()


def test_every_chunk_runs_in_a_worker(pair, base, tmp_path, monkeypatch):
    """Several chunks and threads: this process genotypes no chunk, every chunk runs in a forked worker (as soon as
    one is free, none waiting for this process), and the pair file equals the one-thread run's."""
    monkeypatch.setattr(genotyping, "CHUNK_MIN", 1)  # one window per chunk: the toy's 20 windows in 20 chunks
    ran, real, parent = tmp_path / "ran", genotyping._Genotyper.chunk, os.getpid()

    def chunk(self, job):
        with open(ran, "a") as fh:  # the process that genotypes each chunk
            fh.write(f"{os.getpid()}\n")
        return real(self, job)

    monkeypatch.setattr(genotyping._Genotyper, "chunk", chunk)
    calls = run_pair(tmp_path / "run", pair.a, pair.b, pair.elements, threads=3)
    pids = [int(x) for x in ran.read_text().split()]
    assert len(pids) == 20 and parent not in pids
    assert calls.read_bytes() == base.read_bytes()


def test_large_chunks_are_submitted_first(pair, tmp_path, monkeypatch):
    """Chunks go to the workers largest first (summed window spans), so a large chunk is not
    the pair's straggler; chunk numbers, files and calls are unchanged."""
    seen = []
    real = genotyping._run

    def spy(ctx, jobs, *args):
        seen.extend(sum(w.a_end - w.a_start + w.b_end - w.b_start for w in j.windows if w.status == "ok")
                    for j in jobs)
        return real(ctx, jobs, *args)

    monkeypatch.setattr(genotyping, "_run", spy)
    monkeypatch.setattr(genotyping, "CHUNK_MIN", 1)
    run_pair(tmp_path, pair.a, pair.b, pair.elements)
    assert len(seen) > 3 and seen == sorted(seen, reverse=True) and seen[0] > seen[-1]


def test_unreadable_cached_anchors_are_recomputed(pair, base, tmp_path, caplog):
    """A cached anchors file damaged in place (same size: is_done cannot tell) is
    recomputed with a warning, not a traceback."""
    route = FlankMm2Route()
    with pair_context(tmp_path, pair.a, pair.b, pair.elements) as ctx:
        out = genotype_pair(ctx, route)
        npz = ctx.workdir / "anchors" / "A__B.flank-mm2.npz"
        data = bytearray(npz.read_bytes())
        data[len(data) // 2:len(data) // 2 + 64] = bytes(64)
        npz.write_bytes(bytes(data))
        for p in (out, marker(out)):
            p.unlink()
        with caplog.at_level(logging.WARNING, logger="synltr.routes.assembly.anchor"):
            assert genotype_pair(ctx, route).read_bytes() == base.read_bytes()
        assert "unreadable" in caplog.text and "recomputing" in caplog.text


def test_calls_follow_a_fai_order(tmp_path):
    """Calls are sorted by A's FASTA index order, here chr2 before chr10 (lexicographic
    order would put chr10 first). The two cases have their own seeds: flank pieces alone (the flank-mm2 route)
    cannot tell apart two sequences built from the same seed."""
    parts = {}
    for seqid, case, seed in (("chr2", "E", 3), ("chr10", "P", 4)):
        p = toy.build_pair(random.Random(seed), [case])
        parts[seqid] = (p.a, p.b, [replace(e, seqid=seqid, eid=e.eid.replace("chr1", seqid)) for e in p.elements])
    d = tmp_path
    for name, k in (("A", 0), ("B", 1)):
        toy.write_fasta(d / f"{name}.fa", {seqid: part[k] for seqid, part in parts.items()})
    els = [e for part in parts.values() for e in part[2]]
    with Fasta(d / "A.fa") as fa, Fasta(d / "B.fa") as fb:
        ctx = PairContext(Sample("A", d / "A.fa", None, None, None), Sample("B", d / "B.fa", None, None, None), fa,
                          fb, els, None, d / "out", 1)
        rows = read_rows(genotype_pair(ctx, FlankMm2Route()))
    assert [r["a_seqid"] for r in rows] == ["chr2", "chr10"] and [r["state"] for r in rows] == ["E", "P"]
