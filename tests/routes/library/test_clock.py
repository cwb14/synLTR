"""The library route's alignment jobs (routes.library.clock): whole-record and LTR counts of record pairs (the
clock filter's), ortholog triples projected onto the donor's record, and the self-LTR check; LTRmix's outputs on
three small LTR libraries exactly (stored outputs), at one and three threads."""
import math
from pathlib import Path

import numpy as np
import pytest

from synltr.io.library import read_library
from synltr.routes.library import Records, clock

GOLDEN = Path(__file__).resolve().parents[2] / "golden"


@pytest.fixture(scope="module")
def jobs(tmp_path_factory):
    """The frozen libraries R, S and D (names chr1:start-end, LTRmix's columns) loaded under the screen rule, and the
    frozen inputs and outputs of the jobs."""
    d = tmp_path_factory.mktemp("jobs")
    z = dict(np.load(GOLDEN / "ltrmix_route.npz"))
    recs = {}
    for x in "RSD":
        seq, offs = z[f"jb_{x}_seq"], z[f"jb_{x}_offs"]
        fa, tb = d / f"{x}.fa", d / f"{x}.tsv"
        names = []
        with open(fa, "w") as f, open(tb, "w") as t:
            t.write("seq_id\tk2p\ttsd\tstrand\tltr5_len\tltr3_len\n")
            for i in range(offs.size - 1):
                s = i * 10_000 + 1
                name = f"chr1:{s}-{s + offs[i + 1] - offs[i] - 1}"
                names.append(name)
                f.write(f">{name}\n{seq[offs[i]:offs[i + 1]].tobytes().decode()}\n")
                k = float(z[f"jb_{x}_k2p"][i])
                t.write(f"{name}\t{'NA' if math.isnan(k) else repr(k)}\t.\t{z[f'jb_{x}_strand'][i]}\t"
                        f"{z[f'jb_{x}_l5'][i]}\t{z[f'jb_{x}_l3'][i]}\n")
        recs[x] = Records.of(x, read_library(f"{fa},{tb}", x))
    return recs, z


@pytest.mark.parametrize("threads", [1, 3])
def test_pair_counts_match_ltrmix(jobs, threads):
    """Per pair: aligned columns and mismatches of the two N-stripped records (the second in the first's orientation:
    as given when both strands are known, else the closer by LTRmix's rule; none below half the longer length), then
    the same over the first record's LTRs (none without its bounds)."""
    recs, z = jobs
    pairs = z["jb_pairs"]
    got = clock.pair_counts(recs["R"], recs["D"], pairs[:, 0], pairs[:, 1], threads)
    want = z["jb_pair_out"][:, [0, 1, 3, 4]]
    assert got.dtype == np.int64 and np.array_equal(got, want)
    assert (want[:, 0] == 0).sum() > 10 and (want[:, 2] > 0).sum() > 10


@pytest.mark.parametrize("threads", [1, 3])
def test_triple_counts_match_ltrmix(jobs, threads):
    """Per (donor, recipient, sister) triple: sites aligned in all three on the donor's record, the three pairs'
    mismatches, then ABBA, BABA and BBAA at LTR sites polarized by the donor's other LTR (only where the donor's own
    LTRs pass the self check)."""
    recs, z = jobs
    got = clock.triple_counts(recs["D"], recs["R"], recs["S"], z["jb_triples"], z["jb_self_ok"], threads)
    assert got.dtype == np.int64 and np.array_equal(got, z["jb_triple_out"])
    assert z["jb_triple_out"][:, 4:].sum() > 10


@pytest.mark.parametrize("threads", [1, 3])
def test_self_check_matches_ltrmix(jobs, threads):
    """1 where a record's own two LTRs (at its record ends) give a K2P within SELF_Z standard errors (or one
    substitution) of its annotated K2P: for every record, and for a subset (0 elsewhere)."""
    recs, z = jobs
    n = recs["D"].n
    assert np.array_equal(clock.self_ltr_ok(recs["D"], np.arange(n), threads), z["jb_self_out"])
    got = clock.self_ltr_ok(recs["D"], np.arange(0, n, 2), threads)
    assert got.dtype == np.int8 and np.array_equal(got, z["jb_self_ltr_ok"])
    assert 0 < z["jb_self_out"].sum() < n
