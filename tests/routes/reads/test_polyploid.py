"""Polyploid genotyping (--ploidy 4) on a synthetic autotetraploid with planted allele dosages."""
import csv

import builders
import pytest

from synltr import cli


@pytest.fixture(scope="module")
def tetra(tmp_path_factory):
    d = tmp_path_factory.mktemp("tetra")
    rs = builders.read_set(d, 4)
    sheet = builders.reads_sheet(d / "samples.tsv", rs.ref, rs.gff, {"p4": rs.bam})
    assert cli.main(["genotype", str(sheet), "-o", str(d / "out"), "-t", "2", "--ploidy", "4"]) == 0
    with open(d / "out" / "reads.genotypes.tsv") as fh:
        return {r["element_id"]: r for r in csv.DictReader(fh, delimiter="\t")}


def test_tetraploid_dosage(tetra):
    for eid, st in builders.TETRA_DOSES.items():
        r = tetra[eid]
        want = "/".join(sorted(st, key="FSEN".index))
        if "N" in st:  # a missing haplotype rests on depth: only its presence of E/S matters
            assert "E" not in r["genotype"] and "S" not in r["genotype"], (eid, r["genotype"])
            continue
        assert r["genotype"] == want, (eid, r["genotype"], r["gq"], r["n_F"], r["n_E"], r["n_S"], r["hap_F"])


def test_tetraploid_pav(tetra):
    assert tetra["D1F3E"]["pav"] == "present" and float(tetra["D1F3E"]["p_E"]) > 0.99
    assert tetra["D2F2S"]["pav"] == "present" and float(tetra["D2F2S"]["p_S"]) > 0.99
