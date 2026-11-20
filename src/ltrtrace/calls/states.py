"""Call states and the class of a genotype.

A call holds one allele letter (ALLELES). The class of a set of alleles is `call_class`: present when any is P, S
or T; absent when none is present and at least one is an absence the caller accepts; missing otherwise. With every
absence accepted, observed or inferred, the class is `call_class(absent=OBSERVED_ABSENT | INFERRED_ABSENT)`. The
default, E alone, serves the statistics that were validated on E (an empty site seen with its TSD); a statistic accepts
A only by passing the larger set. Every class decision on these letters goes through this module's sets (PRESENT,
OBSERVED_ABSENT, INFERRED_ABSENT) or `call_class`: the pair calls of both engines, the loci and the catalog, the
catalog VCF and tracks, the statistics, the figures and sim's scores. The read genotyper decides otherwise, on
alleles of its own (F, S, E, N): its reads VCF codes them itself, and its presence/absence column
(routes.reads.model.pav_call: present, solo, absent or uncertain) reads the alleles' posteriors.

FLAGS is the one registry of call flags: those a window call carries and those the catalog adds to a locus.
"""
from __future__ import annotations

from collections.abc import Iterable

ALLELES = "PSTEADOGNU"  # P intact copy, S solo LTR, T truncated, E empty site observed, A absent inferred,
# D deletion, O other sequence, G assembly gap, N null haplotype (reads), U unresolved
PAIR_STATES = tuple(a for a in ALLELES if a != "N")  # what a pair call holds: N is a reads-route state
PRESENT = frozenset("PST")
OBSERVED_ABSENT = frozenset("E")
INFERRED_ABSENT = frozenset("A")
CLASSES = ("present", "absent", "missing")

FLAGS = {
    "adjacent": "an E or O gap reaches past the element end into neighbouring deleted sequence",
    "small_ins": "an O call whose B insertion is only slightly larger than the TSD tolerance",
    "filler": "a D call whose junction holds more than the tolerated B sequence",
    "tsd_match": "an E call whose B site matches an A TSD copy within the allowance",
    "nestB": "a P call with a B insertion inside the element",
    "solo_tsd": "an S call with B's bases on both sides matching A's TSD",
    "inv": "the B window lies on A's reverse strand",
    "host_empty": "an E inherited from the host element's empty site",
    "inherit": "a D, O, G or U inherited from the host element",
    "conflict": "top-tier calls of one locus disagree on the genotype class",
    "no_call": "no member element of the locus was called in this genome",
    "annot_conflict": "another genome's member element calls this member genome empty, deleted or other",
    "multi_member": "a genome holds more than one member element of the locus",
}


def call_class(alleles: Iterable[str], *, absent: frozenset[str] | set[str] = OBSERVED_ABSENT) -> str:
    """"present", "absent" or "missing" for a set of alleles.

    present: any P, S or T. absent: no present allele and at least one allele in `absent` (default E alone;
    OBSERVED_ABSENT | INFERRED_ABSENT adds A). missing: everything else. A letter outside ALLELES is a ValueError.
    """
    seen = set(alleles)
    unknown = seen - set(ALLELES)
    if unknown:
        raise ValueError(f"unknown allele(s) {''.join(sorted(unknown))!r}; alleles are {ALLELES}")
    if seen & PRESENT:
        return "present"
    return "absent" if seen & absent else "missing"
