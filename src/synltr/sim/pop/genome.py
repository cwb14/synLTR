"""Building blocks of the simulators: SimElement (an LTR-RT of a simulated genome, 0-based half-open), reverse
complement of bytes, and base substitution."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from synltr.io.seqio import COMP_FROM, COMP_TO

_COMP = bytes.maketrans(COMP_FROM.encode(), COMP_TO.encode())
_BASES = np.frombuffer(b"ACGT", dtype=np.uint8)
_CODE = np.full(256, 255, dtype=np.uint8)  # A C G T -> 0 1 2 3; anything else (N) is never substituted
_CODE[_BASES] = np.arange(4, dtype=np.uint8)


@dataclass
class SimElement:
    """An LTR-RT of a simulated genome. Coordinates 0-based half-open, LTR spans in genome orientation.

    origin: "backbone" (every element of the simulated backbone). parent: eid of the smallest element holding it
    (the annotation's nesting), None at the top level.
    """

    eid: str
    start: int
    end: int
    ltr_left: tuple[int, int]
    ltr_right: tuple[int, int]
    strand: str
    tsd: str
    k: float | None
    family: str | None
    superfamily: str | None
    origin: str
    parent: str | None = None


def revcomp_bytes(seq: bytes | bytearray) -> bytes:
    return bytes(seq).translate(_COMP)[::-1]


def substitute(seq: bytes, rate: float, rng: np.random.Generator) -> bytes:
    """Each A/C/G/T replaced with probability `rate` by one of the other three."""
    a = np.frombuffer(seq, dtype=np.uint8).copy()
    hit = np.flatnonzero(rng.random(len(a)) < rate)
    code = _CODE[a[hit]]
    ok = code < 4
    a[hit[ok]] = _BASES[(code[ok] + rng.integers(1, 4, int(ok.sum()))) % 4]
    return a.tobytes()
