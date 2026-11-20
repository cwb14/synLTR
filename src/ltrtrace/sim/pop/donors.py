"""Inputs the simulators read from a real annotation: donor templates for new insertions (pick_donors), element
records moved to a backbone segment's coordinates (shift_element), and the backbone elements of one segment of
`ltrtrace sim pop` (backbone_elements); the donor record and the constants of assembly gaps and deletions that
pop's events share; and how an annotation reports an element's ends (jitter_bounds, read_tsd)."""
from __future__ import annotations

import logging
from collections import Counter
from dataclasses import dataclass, replace

import numpy as np

from ltrtrace.base.constants import SIM_MAX_DONORS, TSD_DEFAULT
from ltrtrace.io.elements import Element
from ltrtrace.io.seqio import Fasta, revcomp
from ltrtrace.sim.pop.genome import SimElement

log = logging.getLogger(__name__)

MAX_DONORS = SIM_MAX_DONORS


@dataclass
class Donor:
    """A template element in its own orientation (5' LTR first) and its LTR lengths."""

    seq: str
    ltr_left: int
    ltr_right: int
    family: str | None = None
    superfamily: str | None = None


def shift_element(e: Element, off: int) -> Element:
    """e with every coordinate moved `off` bases down (a segment starting at base off + 1 begins at base 1)."""
    return replace(e, start=e.start - off, end=e.end - off,
                   ltr_left=(e.ltr_left[0] - off, e.ltr_left[1] - off),
                   ltr_right=(e.ltr_right[0] - off, e.ltr_right[1] - off))


def pick_donors(elements: list[Element], fa: Fasta, seed: int) -> list[Donor]:
    """Up to MAX_DONORS top-level, childless, high-confidence elements (any seqid), sampled with the seed, in
    their own orientation; those with an N or without an internal region are skipped. An annotation without
    confidence values (EDTA, BED) gives every top-level, childless element."""
    parents = {e.parent for e in elements if e.parent is not None}
    rated = any(e.confidence is not None for e in elements)
    pool = [e for e in elements
            if e.level == 0 and e.eid not in parents and (e.confidence == "high" or not rated)]
    rng = np.random.default_rng([seed, 0])
    pick = sorted(rng.choice(len(pool), size=min(MAX_DONORS, len(pool)), replace=False).tolist()) \
        if pool else []
    out, skipped = [], Counter()
    for i in pick:
        e = pool[i]
        s = fa.fetch(e.seqid, e.start - 1, e.end)
        left, right = e.ltr_left[1] - e.ltr_left[0] + 1, e.ltr_right[1] - e.ltr_right[0] + 1
        if "N" in s:
            skipped["N"] += 1
        elif not (left > 0 and right > 0 and left + right < len(s)):
            skipped["no_internal"] += 1
        else:
            out.append(Donor(revcomp(s), right, left, e.family, e.superfamily) if e.strand == "-"
                       else Donor(s, left, right, e.family, e.superfamily))
    if skipped:
        log.warning("skipped %d donor candidates (%s)", sum(skipped.values()),
                    ", ".join(f"{k}={v}" for k, v in sorted(skipped.items())))
    return out


def backbone_elements(elements: list[Element], seqid: str, start: int, end: int) -> list[SimElement]:
    """The elements lying wholly inside seqid:start-end (1-based inclusive) as origin-backbone SimElements on the
    segment (0-based half-open, segment coordinates), in (start, -end) order; an element crossing either end is
    left out. Each is named by its source position, "<seqid>:<start>-<end>" (1-based inclusive): unique across
    disjoint segments of one annotation, whose loader drops duplicate spans. Parents follow Element.parent when the
    parent is a backbone element too, else None; records whose LTRs do not lie in order inside the element are
    skipped."""
    length = end - start + 1
    names: dict[str, str] = {}
    out = []
    for e in sorted((e for e in elements if e.seqid == seqid and e.start >= start and e.end <= end),
                    key=lambda x: (x.start, -x.end)):
        name = f"{seqid}:{e.start}-{e.end}"
        e = shift_element(e, start - 1)
        s, t = e.start - 1, e.end
        (a, b), (c, d) = (e.ltr_left[0] - 1, e.ltr_left[1]), (e.ltr_right[0] - 1, e.ltr_right[1])
        if not 0 <= s <= a < b <= c < d <= t <= length:
            continue
        out.append(SimElement(name, s, t, (a, b), (c, d), e.strand, e.tsd or "", e.k, e.family, e.superfamily,
                              "backbone", names.get(e.parent) if e.parent else None))
        names[e.eid] = name
    return out


def jitter_bounds(start: int, end: int, ltr_left: tuple[int, int], ltr_right: tuple[int, int], length: int,
                  jit: bool, ds: int, dt: int) -> tuple[int, int, int, int]:
    """(start, end, left LTR start, right LTR end), 0-based half-open, of an element reported with jitter: when
    `jit`, start moves by ds and end by dt (the caller draws them), each LTR's outer bound with them, clamped so the
    element stays on the sequence of `length` and both LTRs stay non-empty (their inner ends stay); else unmoved."""
    (a, b), (c, d) = ltr_left, ltr_right
    ds = min(max(ds, -start), b - 1 - a) if jit else 0
    dt = max(min(dt, length - end), c + 1 - d) if jit else 0
    return start + ds, end + dt, a + ds, d + dt


def read_tsd(seq, start: int, end: int) -> str | None:
    """The TSD an annotator reads beside reported ends [start, end) (0-based half-open) of `seq` (bytes-like): the
    TSD_DEFAULT bases before start, when they equal the TSD_DEFAULT after end and hold no N; else None."""
    left, right = bytes(seq[start - TSD_DEFAULT:start]), bytes(seq[end:end + TSD_DEFAULT])
    return left.decode() if start >= TSD_DEFAULT and left == right and b"N" not in left else None
