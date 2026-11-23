"""Chromosome correspondence and collinear blocks of one genome pair, learned from its anchors (LTRmix's).

A chromosome pair (x of a, y of b) is kept when it holds at least CHAIN_MIN anchors and more than independence
predicts: Poisson upper tail P(>= n | n_x n_y / N) <= CORR_ALPHA / (chromosome pairs with an anchor), N anchors in
all. Its anchors, placed at their records' midpoints, then give collinear blocks (seq.chain.chain_blocks, steps <=
CHAIN_GAP ranks, both orientations) of at least max(CHAIN_MIN, 1 + the longest chain among the same anchors with their
partners' ranks shuffled, seq.chain.shuffled_chain, NULL_SHUFFLES shuffles) anchors. Chromosome pairs are taken in the
order their first anchor appears, and the shuffles of every kept pair come in turn from the one generator passed in.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from synltr.base.constants import CHAIN_GAP, CHAIN_MIN, CORR_ALPHA, NULL_SHUFFLES
from synltr.seq.chain import chain_blocks, shuffled_chain


@dataclass(eq=False)
class Block:
    """A collinear block: chromosomes, orientation (+1, -1), the extents of its anchors' midpoints on each genome
    (1-based) and its anchors' midpoints in chain order (ax on a, bx on b)."""

    chrom_a: str
    chrom_b: str
    sign: int
    a_start: int
    a_end: int
    b_start: int
    b_end: int
    n_anchor: int
    ax: np.ndarray
    bx: np.ndarray


def collinear_blocks(la, lb, a, b, anc, rng) -> tuple[list[Block], np.ndarray, int]:
    """The blocks of the pair (la, lb: Records) from its anchors `anc` (indices into the candidate pairs a, b), with
    `rng` for the shuffles: (blocks, per candidate pair whether it is an anchor inside a block, chromosome pairs
    kept)."""
    from scipy.stats import poisson  # a command's import, not the command line's

    pairs: dict[tuple, list[int]] = {}
    for i, (x, y) in enumerate(zip(la.chrom[a[anc]], lb.chrom[b[anc]])):
        pairs.setdefault((x, y), []).append(int(anc[i]))
    na: dict = {}
    nb: dict = {}
    for (x, y), ids in pairs.items():
        na[x] = na.get(x, 0) + len(ids)
        nb[y] = nb.get(y, 0) + len(ids)
    total = max(len(anc), 1)
    alpha = CORR_ALPHA / max(len(pairs), 1)
    blocks, in_block, kept = [], np.zeros(len(a), bool), 0
    for (x, y), ids in pairs.items():
        ids = np.array(ids)
        if len(ids) < CHAIN_MIN or poisson.sf(len(ids) - 1, na[x] * nb[y] / total) > alpha:
            continue
        kept += 1
        px, py = la.mid[a[ids]], lb.mid[b[ids]]
        minlen = max(CHAIN_MIN, shuffled_chain(px, rng, n=NULL_SHUFFLES, W=CHAIN_GAP) + 1)
        for chain, sign in chain_blocks(px, py, minlen, W=CHAIN_GAP):
            sel = ids[chain]
            in_block[sel] = True
            ax, bx = la.mid[a[sel]], lb.mid[b[sel]]
            blocks.append(Block(x, y, sign, int(ax.min()), int(ax.max()), int(bx.min()), int(bx.max()), len(sel),
                                ax, bx))
    return blocks, in_block, kept
