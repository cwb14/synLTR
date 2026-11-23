"""Inputs of the subgenomes analysis: the LTR-RT library (FASTA records named seqid:start-end) joined to its annotation
(LTRquest GFF3 or a K2P table) by io.library's phasing rule, the analysis's input checks and the genome's .fai; the
run's log helpers. Coordinates follow io.library: copies 1-based closed, LTR offsets 0-based half-open within the
record, 5' LTR first.

The run logs through the package's logger (base.log): `log` a milestone (INFO, shown by default), `warn` a warning,
`vlog` per-step progress and `verbose()` whether the sanity checks run (DEBUG: -v). An input that fails a check is a
UsageError naming the file or option (exit 2)."""
from __future__ import annotations

import contextlib
import logging
import math
import re

import numpy as np

from synltr.base.constants import K2P_IDENTITY_LIKE, MIN_MATCH
from synltr.base.log import UsageError
from synltr.io.annotation import load_annotations
from synltr.io.library import load_library as read_library
from synltr.io.seqio import read_fai

_log = logging.getLogger(__package__)


def log(msg):
    _log.info(msg)


def warn(msg):
    _log.warning(msg)


def vlog(msg):
    _log.debug(msg)


def verbose() -> bool:
    """Whether -v's per-step progress and sanity checks are on."""
    return _log.isEnabledFor(logging.DEBUG)


@contextlib.contextmanager
def quiet():
    """Silence the run's log inside the block."""
    was, _log.disabled = _log.disabled, True
    try:
        yield
    finally:
        _log.disabled = was


def natural_key(s):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", str(s))]


def group_label(g):
    return f"SG{g + 1}"


def _f(x, fmt="{:.4g}"):
    return "NA" if x is None or (isinstance(x, float) and not math.isfinite(x)) else fmt.format(x)


class Library:
    """LTR-RT copies: the shared store (`store`, io.library.Library) under the names the analysis reads: 2-bit codes
    (concatenated, offs = copy boundaries), LTR record offsets l5 and l3 (0-based half-open, -1 = none), strand_code
    (+1, -1, 0), coordinates (start, end: 1-based closed; ids seqid:start-end; cls: the record name's class), the
    annotation's columns, and the genome's sequences (seqnames, seq_len; seq: each copy's sequence index)."""


# phasing's annotation readers: LTRquest GFF3 (LTR_retrotransposon and repeat_region hosts with a K2P key, kept without
# an LTR pair) or a K2P table
READER = dict(ltr_rule="optional", host_rule="k2p", table_rule="phasing", tsd_rule="acgtn",
              k_keys=("K2P_d", "k2p", "K2P"))


def read_annotation(paths, genome=""):
    """The annotation files (each LTRquest GFF3 or a K2P table, read under READER's rules) as one ElementTable, nested
    as a whole (io.annotation.load_annotations)."""
    return load_annotations(paths, genome, **READER)


def load_library(a):
    """The run's LTR-RT copies: the records of a.ltr_fasta joined to a.annotation (io.library's phasing rule), checked:
    K2P values that look like distances (not identities, none negative), K2P for at least MIN_MATCH of the copies,
    and every sequence carrying a copy in the .fai (a.fai) when given; without one, a sequence's length is the end
    of its last copy. An input that fails a check is a UsageError naming the file (io.library's and io.annotation's
    own errors pass unchanged)."""
    genome = getattr(a, "label", None) or ""
    store = read_library(a.ltr_fasta, read_annotation(a.annotation, genome), rule="phasing",
                         threads=getattr(a, "threads", 1))
    e = store.elements
    chroms, starts, ends = list(e.column("seqid")), e.column("start"), e.column("end")
    lib = Library()
    lib.store, lib.n = store, store.n
    lib.ids = [f"{c}:{s}-{t}" for c, s, t in zip(chroms, starts.tolist(), ends.tolist())]
    lib.cls = [c or "" for c in store.cls]
    lib.start, lib.end = starts, ends
    lib.offs, lib.codes = store.offs, store.codes
    lib.k2p, lib.k2p_se = e.column("k"), e.column("k_se")
    lib.n_sites, lib.n_ts, lib.n_tv = e.column("n_sites"), e.column("n_ts"), e.column("n_tv")
    lib.ltr5_len, lib.ltr3_len = e.column("ltr_left_len"), e.column("ltr_right_len")
    for f in ("family", "clade", "superfamily"):
        setattr(lib, f, [v or "" for v in e.column(f)])
    lib.strand = list(e.column("strand"))
    lib.tsd = np.array([t or "" for t in e.column("tsd")])
    lib.strand_code, lib.l5, lib.l3 = store.orient, store.l5, store.l3
    ok = lib.k2p[np.isfinite(lib.k2p)]
    annotation = ", ".join(map(str, a.annotation))
    if ok.size and (np.mean(ok >= K2P_IDENTITY_LIKE) > 0.5 or np.any(ok < 0)):     # most look like identities
        raise UsageError(f"{annotation}: the divergence column looks like identity or is negative (K2P distances are "
                         "expected, mostly < 0.3)")
    if np.isfinite(lib.k2p).sum() < MIN_MATCH * lib.n:
        raise UsageError(f"{annotation}: K2P missing for {lib.n - np.isfinite(lib.k2p).sum()} of {lib.n} copies")
    lib.fai_given = bool(a.fai)
    if a.fai:
        fai = {}
        for p in a.fai:
            try:
                fai.update(read_fai(p))
            except ValueError as exc:  # names the file and the line
                raise UsageError(str(exc)) from None
        absent = sorted(set(chroms) - set(fai), key=natural_key)
        if absent:
            raise UsageError(f"{', '.join(map(str, a.fai))}: {len(absent)} sequence(s) carrying LTR-RTs are absent from "
                             "the .fai, e.g. " + ", ".join(absent[:5]))
        lib.seqnames = sorted(fai, key=natural_key)
        lib.seq_len = np.array([fai[s] for s in lib.seqnames], np.int64)
    else:
        lib.seqnames = sorted(set(chroms), key=natural_key)
        last = {}
        for c, t in zip(chroms, ends.tolist()):
            last[c] = max(last.get(c, 0), t)
        lib.seq_len = np.array([last[s] for s in lib.seqnames], np.int64)
        log("no .fai: sequence lengths = end of their last LTR-RT; LTR-free sequences are not listed")
    six = {s: i for i, s in enumerate(lib.seqnames)}
    lib.seq = np.array([six[c] for c in chroms], np.int32)
    log(f"{lib.n:,} LTR-RTs on {len(set(chroms))} sequences; {int(np.isfinite(lib.k2p).sum()):,} with K2P")
    return lib
