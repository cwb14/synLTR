"""Collinear chains of shared insertions: the pairing map of near-identical sequences and the homoeologs."""
from types import SimpleNamespace

import numpy as np
import pytest
from analyses.subgenomes import sgtools as sg

from synltr.analyses.subgenomes import homoeologs, inputs, lineages, links, run
from synltr.base import constants as cst
from synltr.seq import kmers, relatives


def _pairing_rows(tmp_path, split):
    """Pairing map of a simulated allotetraploid (direct calls, no pipeline): paired sequence pairs and the true
    homologue pairs (same chromosome number in the two progenitors)."""
    d, truth = sg.simulate(tmp_path, tree={"split": split, "children": ["A", "B"]},
                         mergers=[{"time": min(1e4, split / 5), "lineages": ["A", "B"]}])
    hom = {l.split("\t")[0]: l.split("\t")[2] for l in (tmp_path / "truth" / "chromosomes.tsv").read_text().splitlines()[1:]}
    a = sg.parse_args(["--ltr_fasta", str(d / "ltr.fa"), "--annotation", str(d / "annotation.gff3"), "--fai",
                       str(d / "genome.fai"), "--outdir", str(tmp_path / "o"), "-t", "2"])
    run.setup(a)
    lib = inputs.load_library(a)
    cl = lineages.load_or_build_clades(a, lib)
    lm = links.LinkModel(cl, lib.seq.astype(np.int64), len(lib.seqnames))       # every simulated sequence is phased
    rows = homoeologs.pairing_map(lib, lm, np.random.default_rng(1), 2)
    n = lib.seqnames
    paired = {frozenset((n[r["i"]], n[r["j"]])) for r in rows if r["paired"]}
    true_pairs = {frozenset((x, y)) for x in truth for y in truth if x < y and hom[x] == hom[y] and truth[x] != truth[y]}
    return paired, true_pairs


@pytest.mark.slow
@pytest.mark.parametrize("split,expect_paired", [(5e4, True), (2e6, False)])
def test_pairing_map_recent_vs_old_split(tmp_path, split, expect_paired):
    """Homologous chromosomes of subgenomes split 50 kyr ago share most insertions at matching positions: every
    homologue pair is paired and nothing else; after a 2-Myr split no pair diverged after most insertions."""
    paired, true_pairs = _pairing_rows(tmp_path, split)
    assert paired == (true_pairs if expect_paired else set())


def _plant_block(d, src, dst, n_copies):
    """Copy n_copies consecutive LTR-RTs of sequence src onto sequence dst at the same coordinates (a segmental
    duplication between two non-homologous chromosomes): identical records and annotation."""
    fa = (d / "ltr.fa").read_text().split(">")[1:]
    recs = [(r.split("\n", 1)[0], r.split("\n", 1)[1]) for r in fa]
    src_recs = sorted((r for r in recs if r[0].split(":")[0] == src), key=lambda r: int(r[0].split(":")[1].split("-")[0]))
    mid = len(src_recs) // 2 - n_copies // 2
    block = src_recs[mid:mid + n_copies]
    names = {h.split("#")[0] for h, _ in block}
    gff = (d / "annotation.gff3").read_text().splitlines()
    ids, new_gff = {}, []
    for l in gff:
        f = l.split("\t")
        if len(f) == 9 and f[2] == "LTR_retrotransposon" and f"Name={src}:{f[3]}-{f[4]};" in f[8] + ";" \
                and f"{src}:{f[3]}-{f[4]}" in names:
            ids[f[8].split("ID=")[1].split(";")[0]] = True
    for l in gff:
        f = l.split("\t")
        if len(f) == 9 and f[0] == src:
            fid = f[8].split("ID=")[1].split(";")[0] if "ID=" in f[8] else ""
            par = f[8].split("Parent=")[1].split(";")[0] if "Parent=" in f[8] else ""
            if fid in ids or par in ids:
                attr = f[8].replace(f"ID={fid}", f"ID={fid}_dup") if fid else f[8]
                attr = attr.replace(f"Parent={par}", f"Parent={par}_dup") if par else attr
                attr = attr.replace(f"Name={src}:", f"Name={dst}:")
                new_gff.append("\t".join([dst] + f[1:8] + [attr]))
    (d / "annotation.gff3").write_text("\n".join(gff + new_gff) + "\n")
    with open(d / "ltr.fa", "a") as fh:
        for h, s in block:
            fh.write(">" + dst + h[len(src):] + "\n" + s)
    return len(block)


@pytest.mark.slow
def test_pairing_map_rejects_local_duplicated_block(tmp_path):
    """A diploid whose chromosome 2 carries a copied 40-element block of chromosome 1: the block is a significant
    collinear chain of shared insertions, but it spans a small part of the chromosomes, so they are not paired."""
    d, truth = sg.simulate(tmp_path, tree="A", mergers=[])
    s1, s2 = sorted(truth)[:2]
    assert _plant_block(d, s1, s2, 40) == 40
    a = sg.parse_args(["--ltr_fasta", str(d / "ltr.fa"), "--annotation", str(d / "annotation.gff3"), "--fai",
                       str(d / "genome.fai"), "--outdir", str(tmp_path / "o"), "-t", "2"])
    run.setup(a)
    lib = inputs.load_library(a)
    cl = lineages.load_or_build_clades(a, lib)
    lm = links.LinkModel(cl, lib.seq.astype(np.int64), len(lib.seqnames))
    rows = homoeologs.pairing_map(lib, lm, np.random.default_rng(1), 2)
    n = lib.seqnames
    blk = [r for r in rows if {n[r["i"]], n[r["j"]]} == {s1, s2}]
    assert blk and blk[0]["significant"] and blk[0]["chain"] >= 30
    assert not any(r["paired"] for r in rows)


def _clades(members_list, mult=None):
    size = np.array([len(m) for m in members_list])
    mem = np.concatenate([np.array(m) for m in members_list])
    mult = np.ones(len(members_list), np.int64) if mult is None else np.array(mult)
    return kmers.Clades(mem, size, mult, np.full(len(members_list), 2))


def _brute(cl, old, grp, N):
    n = old.size
    R = np.zeros((n, n))
    for k in range(cl.n):
        if cl.size[k] > cst.NMAX_KNN:
            continue
        m = [x for x in cl.members[cl.indptr[k]:cl.indptr[k + 1]] if old[x]]
        if len(m) < 2:
            continue
        w = cl.mult[k] / max(cl.size[k] - 1, 1)
        for a_ in m:
            for b_ in m:
                if a_ != b_:
                    R[a_, b_] += w
    best = np.full((n, N), -1)
    for e in range(n):
        for h in range(N):
            cand = [f for f in range(n) if grp[f] == h and R[e, f] > 0]
            if cand:
                best[e, h] = max(cand, key=lambda f: (R[e, f], -f))  # most related, ties -> lower index
    return best


def test_old_relatives_and_reciprocal_best_hits_match_brute_force():
    rng = np.random.default_rng(1)
    n, N = 60, 3
    grp_all = rng.integers(0, N, n)
    old = rng.random(n) < 0.7
    grp = np.where(old, grp_all, -1)
    cl = _clades([sorted(rng.choice(n, rng.integers(2, 6), replace=False).tolist()) for _ in range(120)]
                 + [[0, 1], [0, 2]],                                    # a tie for copy 0
                 mult=list(rng.integers(1, 4, 120)) + [1, 1])
    want = _brute(cl, old, grp, N)
    for threads, block in ((1, 2e7), (3, 7)):                        # blocks and workers change nothing
        got = homoeologs.old_relatives(cl, old, grp, N, threads=threads, block_nnz=block)
        assert np.array_equal(got, want)
    P = relatives.rbh_anchors(want, grp)
    for a_, b_ in P:
        assert grp[a_] < grp[b_] and want[a_, grp[b_]] == b_ and want[b_, grp[a_]] == a_
    exp = {(a_, b_) for a_ in range(n) for b_ in range(n) if grp[a_] >= 0 and grp[b_] > grp[a_]
           and want[a_, grp[b_]] == b_ and want[b_, grp[a_]] == a_}
    assert {tuple(map(int, p)) for p in P} == exp


def _homoeologs_on(lib, ph, seed=1):
    return homoeologs.run_homoeologs(SimpleNamespace(seed=seed, threads=2), lib, ph)


def _planted(n_pairs=3, per=60, shuffle=None, seed=2):
    """Two subgenomes of n_pairs chromosomes; homoeologs (i, i + n_pairs) share `per` old insertions at collinear
    positions (one pair-private lineage each), plus random young copies; optionally one chromosome's copies shuffled."""
    rng = np.random.default_rng(seed)
    C = 2 * n_pairs
    seq, start, k2p, mem = [], [], [], []
    for p in range(n_pairs):
        for q in range(per):
            pos = q * 50_000 + 1000
            for s in (p, p + n_pairs):
                seq.append(s)
                start.append(pos + int(rng.integers(0, 2000)))
                k2p.append(0.1 + rng.random() * 0.05)
            mem.append([len(seq) - 2, len(seq) - 1])
    for s in range(C):                                                  # young copies, no shared lineages
        for q in range(per):
            seq.append(s)
            start.append(int(rng.integers(0, per * 50_000)))
            k2p.append(rng.random() * 0.01)
    seq, start, k2p = np.array(seq), np.array(start), np.array(k2p)
    if shuffle is not None:
        m = np.flatnonzero(seq == shuffle)
        start[m] = rng.permutation(start[m])
    lib = SimpleNamespace(n=seq.size, seq=seq.astype(np.int32), start=start.astype(np.int64),
                          end=(start + 6000).astype(np.int64), k2p=k2p, seq_len=np.full(C, per * 50_000 + 10_000))
    z = np.array([0] * n_pairs + [1] * n_pairs)
    ph = dict(N=2, C=C, z=z, ek=seq.astype(np.int64), names=[f"c{i}" for i in range(C)],
              strata_edges=[0.005, 0.05], clades_all=_clades(mem), strata_n={})
    return lib, ph


def test_planted_homoeologs_and_shuffled_decoy():
    lib, ph = _planted()
    ho = _homoeologs_on(lib, ph)
    for i in range(6):
        c = ho["calls"][(i, 1 - ph["z"][i])]
        assert c["status"] == "homoeolog" and c["partner"] == (i + 3) % 6 and c["reciprocal"]
    assert len(ho["sets"]) == 3 and all(s["complete"] for s in ho["sets"])
    lib, ph = _planted(shuffle=1)                                       # same anchors, positions scrambled on c1
    ho = _homoeologs_on(lib, ph)
    assert ho["calls"][(1, 1)]["status"] == "unresolved" and ho["calls"][(4, 0)]["status"] == "unresolved"
    assert ho["calls"][(0, 1)]["partner"] == 3 and ho["calls"][(2, 1)]["partner"] == 5


def test_no_shared_insertions_no_calls():
    lib, ph = _planted()
    ph["clades_all"] = _clades([[0, 2]])                               # one lineage within a subgenome
    ho = _homoeologs_on(lib, ph)
    assert {c["status"] for c in ho["calls"].values()} == {"none"} and not ho["sets"]
    assert _homoeologs_on(lib, dict(ph, N=1)) is None
