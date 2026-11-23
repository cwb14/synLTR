"""Union-find forests: the root of a node, the path to it halved on the way."""

from __future__ import annotations


def find(parent, x):
    """The root of node `x` in the forest `parent` (a list or array of parent indexes, a root its own parent); each node
    passed on the way is pointed to its grandparent (path halving), so later finds are shorter."""
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x
