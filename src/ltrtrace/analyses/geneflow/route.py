"""Route properties: the internal parameters through which the one gene-flow engine reproduces each evidence route's
validated numbers. They are fixed per route and documented, never user options. The engine honors each one, and
refuses a value it does not implement (`require`) rather than compute another statistic. Where each is honored:
  roles           run.analyses (trio.trios) and run.screen_analyses (library_stats.roles)
  alphabet        run's paintings, one per alphabet in this order; paint.classes and paint.paint (paint.ALPHABETS)
  decoy_scale, donor_state, tract_position, tract_end   paint.paint (and paint.call_tracts for the tract end)
  posterior_mean  paint.posterior_mean, read by the painting rows of both routes and by the screen's overview
  block           route.block_bp and library_stats.block_bp
  inversions      run.analyses (inversions.inverted_blocks) and run.screen_analyses (inversions.chain_inversions)

The assembly route (ROUTE_PROPS["assembly"]):
  roles           "shared_derived": the sister pair shares the most insertions the third genome lacks (merged
                  states, no outgroup present), with a two-sided binomial check of the top count against the
                  second; the recipient is the sister with the larger donor-shared count, by the sign of D_pav
                  (trio.trios)
  alphabet        ("insertion",): one painting, its classes RD (donor present, sister E) and RS (the mirror), loci
                  an outgroup carries set aside (paint.classes)
  decoy_scale     "informative": the target-decoy scale is the recipient's informative (RD, RS) loci over the
                  decoy's
  donor_state     "rd_up": the fitted state 1 is donor-like when its RD emission exceeds state 0's
  tract_position  "site_start": a locus sits at its site's start on the painted genome (the insertion point)
  tract_end       "largest_end": a tract ends at the largest end of its loci
  posterior_mean  "unweighted": posterior means count each painted locus once
  block           "fai": jackknife blocks of BLOCK_BP, or 1/BLOCK_FRAC of the smallest genome by its FASTA index
                  (at least BLOCK_MIN; stats.resample.block_size)
  inversions      "reversed_runs": inverted blocks are reversed runs of resolved loci between syntenic partner
                  sequences (inversions.inverted_blocks)

The library route (ROUTE_PROPS["library"]; the screen's evidence: routes.library):
  roles           "single_copy_kmers": the sister pair shares the most single-copy k-mers that the third genome
                  lacks; the recipient is the sister by the sign of D_kmer on the donor's records
  alphabet        ("ancestry", "recent"): two paintings, each its own HMM: "ancestry" with the classes near_donor,
                  donor_only, near_sister, sister_only ("near": divergence below the 5th percentile of the
                  sister-donor orthologs' divergence) and "recent" with near_donor and near_sister
  decoy_scale     "all": the target-decoy scale is the recipient's elements over the decoy's
  donor_state     "near_donor_up_near_sister_down": the fitted state 1 is donor-like when its near_donor emissions
                  exceed state 0's and its near_sister emissions fall below them
  tract_position  "midpoint": an element sits at its midpoint on the painted genome
  tract_end       "last_end": a tract ends at its last element's end
  posterior_mean  "length_weighted": posterior means weight each element by its length
  block           "placed": jackknife blocks of BLOCK_BP, or 1/BLOCK_FRAC of the smallest genome's placed span (the
                  sum over its sequences of their last element's end; a record whose name holds no coordinates is
                  its own sequence, ending at its length), at least BLOCK_MIN, every genome of the sheet counted
                  (library_stats.block_bp)
  inversions      "chain_blocks": inverted blocks are the collinear chain blocks whose sign is opposite to their
                  chromosome pair's main orientation (blocks.tsv)
The library route's other properties belong to its evidence and are documented where they live: the record sequence
and the element-vs-element divergence (routes.library), inferred absences (A, never an observed empty site; the
statistics that accept them: library_stats.ABSENCES) and the loci of the insertion D (each sister's own library:
library_stats, D_ins).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from ltrtrace.stats.resample import block_size

if TYPE_CHECKING:
    from ltrtrace.analyses.geneflow.catalog import Catalog


@dataclass(frozen=True)
class RouteProps:
    """One evidence route's properties (module docstring)."""

    roles: str
    alphabet: tuple[str, ...]  # the painting alphabets, one painting each, in this order
    decoy_scale: str
    donor_state: str
    tract_position: str
    tract_end: str
    posterior_mean: str
    block: str
    inversions: str


ROUTE_PROPS: dict[str, RouteProps] = {
    "assembly": RouteProps(roles="shared_derived", alphabet=("insertion",), decoy_scale="informative",
                           donor_state="rd_up", tract_position="site_start", tract_end="largest_end",
                           posterior_mean="unweighted", block="fai", inversions="reversed_runs"),
    "library": RouteProps(roles="single_copy_kmers", alphabet=("ancestry", "recent"), decoy_scale="all",
                          donor_state="near_donor_up_near_sister_down", tract_position="midpoint",
                          tract_end="last_end", posterior_mean="length_weighted", block="placed",
                          inversions="chain_blocks"),
}


def require(name: str, value: str, implemented: tuple[str, ...]) -> None:
    """A ValueError when a route property's `value` is not one the engine implements."""
    if value not in implemented:
        raise ValueError(f"{name} {value!r} is not implemented by the gene-flow engine ({', '.join(implemented)})")


def block_bp(table: Catalog, props: RouteProps) -> int:
    """The jackknife block length of a run under its route's block rule (module docstring)."""
    require("jackknife block rule", props.block, ("fai",))
    return block_size([int(np.sum(v, dtype=np.int64)) for v in table.lengths.values()])
