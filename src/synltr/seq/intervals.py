"""Half-open intervals [lo, hi) on one sequence: their union."""

from __future__ import annotations

from collections.abc import Iterable


def union(spans: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    """The union of half-open intervals as sorted, disjoint (lo, hi) pairs: overlapping or touching intervals merge,
    empty ones (hi <= lo) drop out."""
    out: list[list[int]] = []
    for lo, hi in sorted(spans):
        if hi <= lo:
            continue
        if out and lo <= out[-1][1]:
            out[-1][1] = max(out[-1][1], hi)
        else:
            out.append([lo, hi])
    return [(lo, hi) for lo, hi in out]
