"""Nesting of elements: which element contains which (genome coordinates, 1-based closed): LTRpav's parent and level
(`assign_nesting`) and ltrkin's containment counts (`containment`)."""
from __future__ import annotations

import numpy as np

from synltr.io.elements import ElementTable


def assign_nesting(table: ElementTable) -> ElementTable:
    """Set `parent` (eid of the innermost containing element) and `level` on each element; the rows come back
    sorted by (seqid, start, -end).

    Sweeps per seqid in that order with a stack. Partially overlapping elements stay on the stack but are never
    parents.
    """
    seqid, start, end, eid = (table.column(c) for c in ("seqid", "start", "end", "eid"))
    order = sorted(range(len(table)), key=lambda i: (seqid[i], start[i], -end[i]))
    parent = np.empty(len(order), dtype=object)
    level = np.zeros(len(order), dtype=np.int64)
    stack: list[tuple[int, str, int]] = []  # (end, eid, level) of the elements swept so far
    current = None
    for pos, i in enumerate(order):
        if seqid[i] != current:
            stack, current = [], seqid[i]
        while stack and stack[-1][0] < start[i]:
            stack.pop()
        up = next((p for p in reversed(stack) if p[0] >= end[i]), None)
        parent[pos], level[pos] = (up[1], up[2] + 1) if up else (None, 0)
        stack.append((int(end[i]), eid[i], int(level[pos])))
    return table.select(order).replace(parent=parent, level=level)


def containment(table: ElementTable) -> tuple[list[int], list[int], list[list[tuple[int, int]]]]:
    """ltrkin's nesting counts of each element, in table order: n_inner (the elements it holds, at any depth),
    n_outer (the elements that hold it) and inner (the spans (start, end), 1-based closed, of the elements it holds,
    in sweep order).

    Sweeps each seqid in (start, -end) order, keeping the elements that still reach the current start: an earlier
    element holds a later one that ends at or before its end. Partial overlaps stay in the sweep and count for
    neither side, so an element inside two partly overlapping ones has n_outer 2, where `assign_nesting` gives it
    one parent and level 1.
    """
    seqid, start, end = table.column("seqid"), table.column("start").tolist(), table.column("end").tolist()
    order = sorted(range(len(table)), key=lambda i: (seqid[i], start[i], -end[i]))
    n_inner, n_outer = [0] * len(table), [0] * len(table)
    inner: list[list[tuple[int, int]]] = [[] for _ in range(len(table))]
    reach: list[int] = []  # the elements of the current seqid that still reach the current start
    current = None
    for i in order:
        if seqid[i] != current:
            reach, current = [], seqid[i]
        reach = [o for o in reach if end[o] >= start[i]]
        for o in reach:
            if end[i] <= end[o]:
                n_inner[o] += 1
                inner[o].append((start[i], end[i]))
                n_outer[i] += 1
        reach.append(i)
    return n_inner, n_outer, inner
