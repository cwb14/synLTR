"""Route properties (analyses.geneflow.route): the assembly and library routes' values, fixed and never options,
and the engine refusing a value it does not implement."""
from __future__ import annotations

import dataclasses

import pytest
from analyses.geneflow.catalog_toy import loc, write_run

from ltrtrace.analyses.geneflow.paint import classes, paint
from ltrtrace.analyses.geneflow.route import ROUTE_PROPS, RouteProps, block_bp
from ltrtrace.analyses.geneflow.trio import Trio


def test_assembly_route_props():
    props = ROUTE_PROPS["assembly"]
    assert dataclasses.asdict(props) == {
        "roles": "shared_derived", "alphabet": ("insertion",), "decoy_scale": "informative", "donor_state": "rd_up",
        "tract_position": "site_start", "tract_end": "largest_end", "posterior_mean": "unweighted",
        "block": "fai", "inversions": "reversed_runs"}
    with pytest.raises(dataclasses.FrozenInstanceError):
        props.alphabet = "ancestry"


def test_library_route_props():
    props = ROUTE_PROPS["library"]
    assert dataclasses.asdict(props) == {
        "roles": "single_copy_kmers", "alphabet": ("ancestry", "recent"), "decoy_scale": "all",
        "donor_state": "near_donor_up_near_sister_down", "tract_position": "midpoint", "tract_end": "last_end",
        "posterior_mean": "length_weighted", "block": "placed", "inversions": "chain_blocks"}
    with pytest.raises(dataclasses.FrozenInstanceError):
        props.block = "fai"


def test_block_bp_from_the_fasta_indexes(tmp_path, toy_fastas):
    """The assembly route's jackknife block: BLOCK_BP, or 1/BLOCK_FRAC of the smallest genome by its FASTA index
    (at least BLOCK_MIN): 40 kb for the toy's four 1-Mb sequences, the 10-kb floor for one 500-kb sequence."""
    cat = write_run(tmp_path / "a", [loc(10_001, "R:P* S:E D:P* O:E")], fasta_dir=toy_fastas).load()
    assert block_bp(cat, ROUTE_PROPS["assembly"]) == 40_000
    small = write_run(tmp_path / "b", [loc(10_001, "R:P* S:E D:P* O:E")], seqids=("chr1",), seq_len=500_000,
                      fasta_dir=toy_fastas).load()
    assert block_bp(small, ROUTE_PROPS["assembly"]) == 10_000
    with pytest.raises(ValueError, match="jackknife block rule 'placed'"):
        block_bp(cat, dataclasses.replace(ROUTE_PROPS["assembly"], block="placed"))


def test_engine_refuses_properties_it_does_not_implement(tmp_path, toy_fastas):
    """A value the engine does not implement, an alphabet of the other input (the library route's records, not a
    catalog's loci), and a donor-state check reading classes the alphabet lacks are ValueErrors before anything is
    painted."""
    cat = write_run(tmp_path, [loc(1000 * i, "R:P* S:E D:P O:E") for i in range(1, 4)]
                    + [loc(10_000 + 1000 * i, "R:P* S:P D:E O:E") for i in range(3)], fasta_dir=toy_fastas).load()
    with pytest.raises(ValueError, match="painting alphabet 'ancestry' does not paint a run's catalog"):
        classes(cat, "R", "S", "D", alphabet="ancestry")
    for field, value in (("decoy_scale", "decoy"), ("donor_state", "near_donor_up"),
                         ("tract_position", "site_end"), ("tract_end", "first_end"),
                         ("posterior_mean", "median")):
        props = dataclasses.replace(ROUTE_PROPS["assembly"], **{field: value})
        with pytest.raises(ValueError, match=f"{field.replace('_', ' ')} '{value}'"):
            paint(cat, Trio("R", "S", "D", True), props)
    with pytest.raises(ValueError, match="painting alphabet 'words'"):
        paint(cat, Trio("R", "S", "D", True), dataclasses.replace(ROUTE_PROPS["assembly"], alphabet=("words",)))
    with pytest.raises(ValueError, match="does not paint a run's catalog"):
        paint(cat, Trio("R", "S", "D", True), dataclasses.replace(ROUTE_PROPS["assembly"], alphabet=("recent",)))
    with pytest.raises(ValueError, match="reads the class\\(es\\) near_donor, near_sister"):
        paint(cat, Trio("R", "S", "D", True),
              dataclasses.replace(ROUTE_PROPS["assembly"], donor_state="near_donor_up_near_sister_down"))
    assert isinstance(ROUTE_PROPS["assembly"], RouteProps) and set(ROUTE_PROPS) == {"assembly", "library"}
