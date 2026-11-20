"""The library route's orthology (routes.library.orthology): LTRmix's ortholog tests on its synthetic trio, and
LTRmix's own arrays for every pair of that trio, unplaced records included (stored outputs)."""
import re
from pathlib import Path

import numpy as np
import pytest

from ltrtrace.calls.records import read_calls
from ltrtrace.io.elements import read_elements
from ltrtrace.io.library import read_library
from ltrtrace.routes.library import Records
from ltrtrace.routes.library.orthology import orthology

GOLDEN = Path(__file__).resolve().parents[2] / "golden"


def _states(run, a, b):
    """{element name: state} of a's elements called in b, with each element's start and chromosome number."""
    names = {e.eid: e for e in read_elements(run / "work" / "elements" / f"{a}.tsv")}
    out = {}
    for c in read_calls(run / "pairs" / f"{a}__{b}.calls.tsv.gz"):
        e = names[c["eid"]]
        out[e.name] = (int(re.search(r"chr(\d+)", e.seqid)[1]), e.start, c["state"])
    return out


def test_orthology_precision(admixed):
    """R's elements called present in D (state P) against the truth, decoy copies (no orthologs by construction)
    left out: precision > 0.95 and recall > 0.8."""
    sc, run = admixed
    calls = {(c, s): st for name, (c, s, st) in _states(run, "R", "D").items() if "fam9" not in name}
    truth = [t for t in sc.truth if "R" in t["carriers"]]
    hits = [(calls[t["chrom"], t["pos"]] == "P", "D" in t["carriers"]) for t in truth
            if (t["chrom"], t["pos"]) in calls]
    called, true = sum(c for c, _ in hits), sum(t for _, t in hits)
    both = sum(c and t for c, t in hits)
    prec, rec = both / max(called, 1), both / max(true, 1)
    assert len(hits) > 600 and prec > 0.95 and rec > 0.8, (prec, rec)


def test_decoy_not_orthologous(admixed):
    """The decoy family's copies (near-identical, at random positions) are rarely called present in D."""
    _, run = admixed
    decoys = [st for name, (_, _, st) in _states(run, "R", "D").items() if "fam9" in name]
    assert len(decoys) == 40 and sum(st == "P" for st in decoys) / len(decoys) < 0.2


FIELDS = ("orth_b", "orth_a", "status_a", "status_b", "room_a", "room_b")


def _same(po, z, key):
    for f in FIELDS:
        assert np.array_equal(getattr(po, f), z[f"{key}_{f}"]), f
    assert np.array_equal(po.div_a, z[f"{key}_div_a"], equal_nan=True)
    blocks = [(b.sign, b.a_start, b.a_end, b.b_start, b.b_end, b.n_anchor) for b in po.blocks]
    assert blocks == [tuple(r) for r in z[f"{key}_blocks"].tolist()]
    assert [f"{b.chrom_a}\t{b.chrom_b}" for b in po.blocks] == z[f"{key}_block_chroms"].tolist()


@pytest.fixture(scope="module")
def trio_records(admixed):
    sc, _ = admixed
    return {g: Records.of(g, read_library(sc.libraries[g], g)) for g in "RSD"}


@pytest.mark.parametrize("pair", ["RS", "RD", "SD"])
def test_orthology_matches_ltrmix(trio_records, pair):
    """Each pair of LTRmix's trio: ortholog of every record both ways, status (2 present, 1 absent inside a
    collinear window, 0 unresolved), whole-record divergence, the room between flanking anchors and the collinear
    blocks, exactly LTRmix's (frozen), at 2 threads."""
    a, b = pair
    po = orthology(trio_records[a], trio_records[b], threads=2)
    with np.load(GOLDEN / "ltrmix_route.npz") as z:
        _same(po, z, f"tr_{pair}")
    assert (po.status_a == 2).sum() > 300 and len(po.blocks) >= 3


def test_unplaced_records_as_ltrmix(admixed, trio_records, tmp_path):
    """Records whose names hold no coordinates (copies of R records, so they share k-mers) take part in the k-mers
    and anchors, on a chromosome named after themselves: every array of the R-D pair is LTRmix's, theirs unresolved."""
    sc, _ = admixed
    fa, tb = sc.libraries["R"].split(",")
    with np.load(GOLDEN / "ltrmix_route.npz") as z:
        seq, offs = z["tr_unplaced_seq"], z["tr_unplaced_offs"]
        extra = "".join(f">contig_{k}\n{seq[offs[k]:offs[k + 1]].tobytes().decode()}\n" for k in range(offs.size - 1))
        (tmp_path / "Ru.fa").write_text(Path(fa).read_text() + extra)
        lib = read_library(f"{tmp_path / 'Ru.fa'},{tb}", "R")
        assert (~lib.placed).sum() == 4
        po = orthology(Records.of("R", lib), trio_records["D"], threads=2)
        _same(po, z, "tr_uD")
    assert po.status_a[-4:].tolist() == [0, 0, 0, 0]


def test_one_generator_per_genome_pair(trio_records, monkeypatch):
    """The shuffles of all a pair's chromosome pairs come in turn from one generator, default_rng(SCREEN_SEED):
    LTRmix's blocks (frozen) on anchors where the shuffled-null minimum decides which blocks are kept; another seed, or
    a fresh generator per chromosome pair, keeps other blocks. orthology hands each genome pair one new generator."""
    from types import SimpleNamespace

    from ltrtrace.base.constants import SCREEN_SEED
    from ltrtrace.routes.library import correspondence
    from ltrtrace.routes.library import orthology as module

    with np.load(GOLDEN / "ltrmix_route.npz") as z:
        la = SimpleNamespace(chrom=z["nb_chrom_a"].astype(object), mid=z["nb_mid_a"])
        lb = SimpleNamespace(chrom=z["nb_chrom_b"].astype(object), mid=z["nb_mid_b"])
        anc, want = z["nb_anc"], ([tuple(b) for b in z["nb_blocks"].tolist()], z["nb_block_chroms"].tolist())
    ab = np.arange(la.mid.size, dtype=np.int64)

    def blocks(rng):
        found, _, _ = correspondence.collinear_blocks(la, lb, ab, ab, anc, rng)
        return ([(b.sign, b.a_start, b.a_end, b.b_start, b.b_end, b.n_anchor) for b in found],
                [f"{b.chrom_a}\t{b.chrom_b}" for b in found])

    assert blocks(np.random.default_rng(SCREEN_SEED)) == want
    assert blocks(np.random.default_rng(SCREEN_SEED + 1)) != want
    real = correspondence.shuffled_chain
    with monkeypatch.context() as m:
        m.setattr(correspondence, "shuffled_chain", lambda pos, rng, **kw: real(pos, np.random.default_rng(0), **kw))
        assert blocks(np.random.default_rng(SCREEN_SEED)) != want
    states, run = [], module.collinear_blocks
    monkeypatch.setattr(module, "collinear_blocks",
                        lambda *a: states.append(a[-1].bit_generator.state) or run(*a))
    module.orthology(trio_records["R"], trio_records["S"], threads=2)
    module.orthology(trio_records["R"], trio_records["D"], threads=2)
    fresh = np.random.default_rng(SCREEN_SEED).bit_generator.state
    assert states == [fresh, fresh]
