"""k-mer lineages (seq.kmers): carrier sets (clades) of shared canonical k-mers, exact against brute force and against
the arrays phasing pinned on the builders' polyploid scenario; the one 64-bit mixer; bottom-s sketches; the clade
checkpoint; k-mers shared across libraries (candidate record pairs, single-copy tallies), exact against brute force
and against LTRmix's own outputs."""
import json
from collections import defaultdict
from itertools import combinations
from pathlib import Path
from types import SimpleNamespace

import numba
import numpy as np
import pytest
from analyses.subgenomes import sgtools as sg

from ltrtrace.io.seqio import encode
from ltrtrace.seq import kmers

GOLDEN = Path(__file__).resolve().parent.parent / "golden"
TINY = GOLDEN / "nucleus_tiny" / "lib"
PINNED = GOLDEN / "clades_nucleus_tiny.npz"
SETTINGS = {"default": {"k": 21, "cap": 250, "scale": 1, "mem_gb": 8.0},     # phasing's run settings
            "passes": {"k": 15, "cap": 40, "scale": 1, "mem_gb": 1e-4},      # many bounded passes
            "sampled": {"k": 21, "cap": 250, "scale": 3, "mem_gb": 8.0}}     # FracMinHash, 1/3 of the k-mers
FIELDS = ("members", "size", "mult", "nseq")
M64 = (1 << 64) - 1


def _brute(seqs, seqof, k, cap, min_seq=2):
    rc = str.maketrans("ACGT", "TGCA")
    d = defaultdict(set)
    for i, s in enumerate(seqs):
        for p in range(len(s) - k + 1):
            w = s[p:p + k]
            if "N" in w:
                continue
            d[min(w, w.translate(rc)[::-1])].add(i)
    cl = defaultdict(int)
    for cs in d.values():
        if 2 <= len(cs) <= cap and len({seqof[c] for c in cs}) >= min_seq:
            cl[frozenset(cs)] += 1
    return dict(cl)


def _as_dict(cl):
    return {frozenset(cl.members[cl.indptr[i]:cl.indptr[i + 1]].tolist()): int(cl.mult[i]) for i in range(cl.n)}


def _fmix64(x):
    """MurmurHash3's 64-bit finalizer on Python integers."""
    x ^= x >> 33
    x = (x * 0xFF51AFD7ED558CCD) & M64
    x ^= x >> 33
    x = (x * 0xC4CEB9FE1A85EC53) & M64
    return x ^ (x >> 33)


@pytest.fixture(scope="module")
def tiny(tmp_path_factory):
    """The builders' polyploid scenario pinned in tests/golden/nucleus_tiny (136 copies on six sequences), loaded as
    the subgenomes analysis loads it: the shared store and each copy's sequence index."""
    from ltrtrace.analyses.subgenomes import inputs, run

    a = run.options([TINY / "ltr.fa"], [TINY / "annotation.gff3"], tmp_path_factory.mktemp("tiny") / "o",
                    fai=[TINY / "genome.fai"], threads=1)
    return inputs.load_library(a)


@pytest.mark.parametrize("name", list(SETTINGS))
def test_clades_match_phasing(tiny, name):
    """The phasing analysis's clades of the tiny polyploid, pinned in tests/golden/clades_nucleus_tiny.npz, for its run
    settings, many bounded passes and FracMinHash sampling: same values, same dtypes, same order."""
    cl = kmers.clades(tiny.store, tiny.seq, threads=2, **SETTINGS[name])
    with np.load(PINNED) as z:
        for f in FIELDS:
            want = z[f"{name}_{f}"]
            assert getattr(cl, f).dtype == want.dtype and np.array_equal(getattr(cl, f), want), f
        assert (cl.n_occ, cl.n_multi) == (int(z[f"{name}_n_occ"]), int(z[f"{name}_n_multi"]))
    assert cl.scale == SETTINGS[name]["scale"]


def test_one_mixer_is_fmix64_on_scalars_and_arrays():
    """One mixer for numba kernels (uint64 scalars) and numpy code (uint64 arrays): MurmurHash3's fmix64."""
    rng = np.random.default_rng(5)
    xs = [0, 1, 2, 3, 0x9E3779B97F4A7C15, M64, 1 << 63] + [int(v) for v in rng.integers(0, 2 ** 63, 20, np.int64)]
    want = [_fmix64(x) for x in xs]
    assert [int(kmers.mix64(np.uint64(x))) for x in xs] == want
    got = kmers.mix64(np.array(xs, np.uint64))
    assert got.dtype == np.uint64 and got.tolist() == want


def test_clades_match_brute_force():
    rng = np.random.default_rng(1)
    seqs = sg.random_family(rng) + sg.random_family(rng, 30)
    seqs[3] = seqs[3][:100] + "N" + seqs[3][101:]
    seqof = rng.integers(0, 4, len(seqs)).tolist()
    cl = kmers.clades(sg.lib_from_seqs(seqs, seqof), seqof, k=15, cap=25, threads=1)
    assert _as_dict(cl) == _brute(seqs, seqof, 15, 25) and cl.n > 10


def test_clades_min_seq():
    """min_seq sets how many sequences a carrier set must span (phasing: 2)."""
    rng = np.random.default_rng(6)
    seqs = sg.random_family(rng, 24, 200)
    seqof = rng.integers(0, 4, len(seqs)).tolist()
    lib = sg.lib_from_seqs(seqs, seqof)
    for m in (1, 2, 3):
        assert _as_dict(kmers.clades(lib, seqof, k=13, cap=30, min_seq=m, threads=1)) == _brute(seqs, seqof, 13, 30, m)


def test_clades_pass_invariance():
    rng = np.random.default_rng(2)
    seqs = sg.random_family(rng, 60)
    seqof = rng.integers(0, 5, len(seqs)).tolist()
    lib = sg.lib_from_seqs(seqs, seqof)
    one = _as_dict(kmers.clades(lib, seqof, k=15, cap=40, mem_gb=10.0, threads=1))
    many = _as_dict(kmers.clades(lib, seqof, k=15, cap=40, mem_gb=1e-6, threads=1))      # forces many passes
    assert one == many


def test_clades_skip_huge_runs():
    rep = "ACGTTGCAAGT" * 6000          # each 11-mer occurs ~18,000 times over three copies (> RUN_LIM): skipped
    uniq = "GATTACAGGCATTAGCCTAGGATCCAGT"
    seqs = [rep, rep, rep, uniq, uniq]
    cl = kmers.clades(sg.lib_from_seqs(seqs, [0, 1, 2, 0, 1]), [0, 1, 2, 0, 1], k=11, cap=250, threads=1)
    assert set(_as_dict(cl)) == {frozenset({3, 4})}
    assert frozenset({0, 1, 2}) in _brute(seqs, [0, 1, 2, 0, 1], 11, 250)     # brute force keeps the tandem set


def test_clades_restrict():
    rng = np.random.default_rng(3)
    seqs = sg.random_family(rng, 20)
    seqof = list(range(4)) * 5
    cl = kmers.clades(sg.lib_from_seqs(seqs, seqof), seqof, k=15, cap=40, threads=1)
    keep = np.zeros(len(seqs), bool)
    keep[:10] = True
    r = cl.restrict(keep)
    assert r.n > 0 and set(r.members.tolist()) <= set(range(10))
    assert all(r.size[i] >= 2 for i in range(r.n))


def test_auto_scale():
    assert kmers.auto_scale(2_000_000) == 1 and kmers.auto_scale(150_000_000) == 1
    assert kmers.auto_scale(153_000_000) == 2 and kmers.auto_scale(4_640_000_000) == 31


def test_restrict_merges_sets_equal_after_masking():
    """Carrier sets that differ only by excluded copies become one clade (mult summed), so a lineage reaching
    unphased scaffolds is not counted twice."""
    cl = kmers.Clades(np.array([0, 1, 5, 0, 1, 6, 2, 3]), np.array([3, 3, 2]), np.array([4, 7, 1]), np.array([2, 2, 2]))
    keep = np.ones(7, bool)
    keep[[5, 6]] = False
    r = cl.restrict(keep, np.array([0, 1, 0, 1, 0, 2, 2]))
    got = {frozenset(r.members[r.indptr[i]:r.indptr[i + 1]].tolist()): int(r.mult[i]) for i in range(r.n)}
    assert got == {frozenset({0, 1}): 11, frozenset({2, 3}): 1}


def test_clades_auto_scale_and_threads(tiny):
    """scale "auto" samples only above KMER_BUDGET bases (exact here); the result does not depend on threads."""
    a = kmers.clades(tiny.store, tiny.seq, threads=1)
    b = kmers.clades(tiny.store, tiny.seq, threads=3)
    assert a.scale == 1 and all(np.array_equal(getattr(a, f), getattr(b, f)) for f in FIELDS)


def test_clades_checkpoint(tiny, tmp_path, monkeypatch, caplog):
    """With a cache folder the clades are a checkpoint stage: written once with its done marker, loaded while the
    key holds (library, parameters, code), rebuilt for other parameters or another library, and rebuilt with a
    warning when the file is damaged."""
    first = kmers.clades(tiny.store, tiny.seq, threads=1, cache=tmp_path)
    out = tmp_path / "clades.npz"
    assert out.is_file() and (tmp_path / "clades.npz.done.json").is_file()
    assert json.loads((tmp_path / "clades.npz.done.json").read_text())["outputs"][0]["path"] == "clades.npz"
    build = kmers._build
    monkeypatch.setattr(kmers, "_build", lambda *a, **kw: pytest.fail("rebuilt"))
    again = kmers.clades(tiny.store, tiny.seq, threads=1, cache=tmp_path)
    assert all(np.array_equal(getattr(first, f), getattr(again, f)) for f in FIELDS)
    assert (again.n_occ, again.n_multi, again.scale) == (first.n_occ, first.n_multi, first.scale)
    calls = []
    monkeypatch.setattr(kmers, "_build", lambda *a, **kw: calls.append(kw) or build(*a, **kw))
    kmers.clades(tiny.store, tiny.seq, k=15, threads=1, cache=tmp_path)               # other parameters
    other = np.where(tiny.seq == 0, 1, tiny.seq)                                       # another library partition
    kmers.clades(tiny.store, other, k=15, threads=1, cache=tmp_path)
    assert len(calls) == 2
    kmers.clades(tiny.store, other, k=15, threads=1, cache=tmp_path)
    assert len(calls) == 2
    out.write_bytes(b"\0" * out.stat().st_size)                                          # damaged, same size
    with caplog.at_level("WARNING", logger="ltrtrace"):
        kmers.clades(tiny.store, other, k=15, threads=1, cache=tmp_path)
    assert len(calls) == 3 and "unreadable" in caplog.text


def test_sketch_is_bottom_s_of_mixed_canonical_kmers():
    """Each row: the s smallest distinct fmix64 values of the record's canonical k-mers (2-bit codes, k-mers with
    another base skipped), padded with the largest uint64."""
    rng = np.random.default_rng(8)
    seqs = ["".join(rng.choice(list("ACGT"), n)) for n in (60, 9, 40, 300)]
    seqs[2] = seqs[2][:15] + "N" + seqs[2][16:]
    lib = sg.lib_from_seqs(seqs, [0, 1, 2, 3])
    k, s = 7, 16
    S = kmers.sketch(lib, k=k, s=s)
    code = {"A": 0, "C": 1, "G": 2, "T": 3}
    for i, x in enumerate(seqs):
        vals = set()
        for p in range(len(x) - k + 1):
            w = x[p:p + k]
            if "N" in w:
                continue
            f = sum(code[b] << 2 * (k - 1 - j) for j, b in enumerate(w))
            r = sum((3 - code[b]) << 2 * j for j, b in enumerate(w))
            vals.add(_fmix64(min(f, r)))
        want = sorted(vals)[:s]
        assert S[i].tolist() == want + [M64] * (s - len(want)), i


# ---- k-mers shared across libraries: LTRmix's candidate record pairs and single-copy tallies ----
LTRMIX = GOLDEN / "ltrmix_route.npz"   # LTRmix's own outputs on seeded inputs, stored


def _store(z, prefix):
    """A library store (codes, offs) of records frozen as their letters."""
    return SimpleNamespace(codes=encode(z[f"{prefix}_seq"].tobytes()), offs=z[f"{prefix}_offs"])


def _occurrences(lib, k):
    """canonical k-mer -> the record of each occurrence (brute force)."""
    rc = str.maketrans("ACGT", "TGCA")
    occ = defaultdict(list)
    for r in range(len(lib.offs) - 1):
        s = "".join("ACGTN"[c] for c in lib.codes[lib.offs[r]:lib.offs[r + 1]])
        for p in range(len(s) - k + 1):
            w = s[p:p + k]
            if "N" not in w:
                occ[min(w, w.translate(rc)[::-1])].append(r)
    return occ


def test_cross_pairs_match_ltrmix():
    """Record pairs of two libraries sharing k-mers with at most cmax occurrences in each, weighted by the sum of
    1/(ca cb) per shared occurrence pair (float32 terms, float64 sum), sorted by (a, b): LTRmix's candidate pairs
    exactly (soft-masked letters, N runs, an IUPAC letter, a k-mer twice in one record, an empty record, one shorter
    than k and a family above cmax), in one pass or many bounded ones, at any thread count."""
    with np.load(LTRMIX) as z:
        libs = {x: _store(z, f"km_{x}") for x in "ABC"}
        for x, y in ("AB", "AC", "BC"):
            want = [z[f"km_pairs_{x}{y}_{c}"] for c in "abw"]
            for mem_gb, threads in ((8.0, 1), (1e-5, 3)):
                got = kmers.cross_pairs(libs[x], libs[y], mem_gb=mem_gb, threads=threads)
                assert isinstance(got, kmers.PairWeights)
                assert [v.dtype for v in got] == [np.dtype(np.int64)] * 2 + [np.dtype(np.float64)]
                assert all(np.array_equal(g, w) for g, w in zip(got, want)), (x, y, mem_gb)


def test_cross_pairs_brute_force():
    rng = np.random.default_rng(4)
    fam = "".join(rng.choice(list("ACGT"), 200))
    seqs = [[fam[:120], fam[60:], "ACGT" * 20], [fam[:100] + "N" + fam[101:], fam[150:] * 2, fam[30:90]]]
    a, b = (sg.lib_from_seqs(s, [0] * len(s)) for s in seqs)
    k, cmax = 15, 2
    oa, ob = _occurrences(a, k), _occurrences(b, k)
    want = defaultdict(float)
    for km, xs in oa.items():
        ys = ob.get(km, [])
        if ys and len(xs) <= cmax and len(ys) <= cmax:
            for x in xs:
                for y in ys:
                    want[x, y] += float(np.float32(1.0 / (len(xs) * len(ys))))
    got = kmers.cross_pairs(a, b, k=k, cmax=cmax)
    assert dict(zip(zip(got.a.tolist(), got.b.tolist()), got.w.tolist())) == dict(want) and len(want) > 3
    assert len(set(want.values())) > 2


def test_single_copy_tallies_match_ltrmix():
    """Per record of the focal library: its k-mers that occur once in it and at most once in each other library,
    counted by the other libraries carrying them (keys: tuples of their labels in the libraries' order): LTRmix's
    tallies exactly, in one pass or many."""
    with np.load(LTRMIX) as z:
        libs = {x: _store(z, f"km_{x}") for x in "ABC"}
        for x in "ABC":
            o1, o2 = z[f"km_tally_{x}_others"].tolist()
            want = z[f"km_tally_{x}"]
            for mem_gb, threads in ((8.0, 1), (1e-5, 2)):
                got = kmers.single_copy_tallies(libs, focal=x, mem_gb=mem_gb, threads=threads)
                assert list(got) == [(), (o1,), (o2,), (o1, o2)]
                assert np.array_equal(np.stack(list(got.values()), 1), want), (x, mem_gb)


def test_cross_library_counts_are_exact_only():
    lib = sg.lib_from_seqs(["ACGT" * 10], [0])
    with pytest.raises(ValueError, match="exact"):
        kmers.cross_pairs(lib, lib, scale=2)
    with pytest.raises(ValueError, match="exact"):
        kmers.single_copy_tallies({"a": lib, "b": lib}, focal="a", scale=2)


def _random_libraries(rng, n_libs=4):
    """Small random libraries (four unless asked) sharing two families (mutated copies, N runs, repeats, a high-copy
    family)."""
    fams = [rng.integers(0, 4, 400) for _ in range(3)]
    libs = {}
    for g in range(n_libs):
        seqs = []
        for r in range(int(rng.integers(12, 20))):
            s = fams[r % 3].copy()
            hit = rng.random(s.size) < [0.002, 0.01, 0.04][r % 3]
            s[hit] = (s[hit] + rng.integers(1, 4, int(hit.sum()))) % 4
            t = "".join("ACGT"[c] for c in s[int(rng.integers(0, 50)):])
            if r % 4 == 1:
                p = int(rng.integers(0, len(t) - 20))
                t = t[:p] + "N" * 10 + t[p + 10:]
            if r % 5 == 2:
                t = t + t[:60]
            seqs.append(t)
        seqs += [seqs[0][:200]] * (5 if g % 2 else 1)            # above cmax in some libraries only
        libs["L" + "ABCDEFGH"[g]] = sg.lib_from_seqs(seqs, [0] * len(seqs))
    return libs


def test_shared_kmers_tallies_split_over_every_thread(monkeypatch):
    """The trios' single-copy tallies are counted by work items, groups of trios (each counting only its own rows) by
    chunks of the buckets, about as many as threads: on 6 libraries (20 trios), in many passes, with the per-chunk
    tally copies fitting TALLY_CHUNK_BYTES (one group by as many chunks as threads, int32 copies per chunk), fitting two
    chunks only (groups by two chunks: 3 uneven groups at 6 threads, 4 groups at 8) and forced past it (one chunk per
    group, each adding into the int64 tallies), every trio's tallies equal single_copy_tallies on its three
    libraries."""
    libs = _random_libraries(np.random.default_rng(5), n_libs=6)
    want = {trio: {f: kmers.single_copy_tallies({t: libs[t] for t in trio}, focal=f, k=15) for f in trio}
            for trio in combinations(libs, 3)}
    width = sum(len(libs[t].offs) - 1 for trio in want for t in trio)     # tally rows: every trio's members' records
    real, splits = kmers._joint_tallies, []

    def spy(H, E, bounds, lo, trios, base, gstart, out, nchunk):
        splits.append((gstart.size - 1, nchunk, out.dtype))
        return real(H, E, bounds, lo, trios, base, gstart, out, nchunk)

    monkeypatch.setattr(kmers, "_joint_tallies", spy)
    cases = ([(kmers.TALLY_CHUNK_BYTES, t, (1, t)) for t in (1, 2, 7)]     # (copies' byte limit, threads, split)
             + [(32 * width, 6, (3, 2)), (32 * width, 8, (4, 2))]
             + [(1, t, (t, 1)) for t in (1, 2, 7)])
    for limit, threads, split in cases:
        monkeypatch.setattr(kmers, "TALLY_CHUNK_BYTES", limit)
        splits.clear()
        got = kmers.shared_kmers(libs, k=15, mem_gb=1e-5, threads=threads)
        assert got.passes > 5 and len(got.tallies) == 20 == len(want)
        for trio, by_focal in got.tallies.items():
            for focal, tally in by_focal.items():
                ref = want[trio][focal]
                assert list(tally) == list(ref) and all(np.array_equal(tally[p], ref[p]) for p in ref)
        if threads <= numba.config.NUMBA_NUM_THREADS:                       # fewer threads give a smaller split
            assert set(splits) == {(*split, np.dtype(np.int64 if split[1] == 1 else np.int32))}


def test_shared_kmers_equal_the_pair_and_trio_references(caplog):
    """One bounded scan over every library gives each pair's candidates (cross_pairs on that pair) and each trio's
    single-copy tallies for each member (single_copy_tallies on that trio's libraries alone), at any pass count and
    thread count; LTRmix's own outputs on its three frozen libraries; the passes are logged."""
    libs = _random_libraries(np.random.default_rng(12))
    labels = list(libs)
    for mem_gb, threads in ((8.0, 1), (1e-5, 3)):
        with caplog.at_level("DEBUG", logger="ltrtrace"):
            got = kmers.shared_kmers(libs, k=15, mem_gb=mem_gb, threads=threads)
        assert got.passes == (1 if mem_gb > 1 else got.passes) and (mem_gb > 1 or got.passes > 5)
        assert f"in {got.passes} pass" in caplog.text
        assert list(got.pairs) == [(x, y) for i, x in enumerate(labels) for y in labels[i + 1:]]
        for (x, y), pw in got.pairs.items():
            want = kmers.cross_pairs(libs[x], libs[y], k=15)
            assert all(np.array_equal(g, w) for g, w in zip(pw, want)) and pw.a.size > 0
        assert len(got.tallies) == 4
        for trio, by_focal in got.tallies.items():
            for focal, tally in by_focal.items():
                want = kmers.single_copy_tallies({t: libs[t] for t in trio}, focal=focal, k=15)
                assert list(tally) == list(want) and all(np.array_equal(tally[p], want[p]) for p in want)
    with np.load(LTRMIX) as z:
        frozen = {x: _store(z, f"km_{x}") for x in "ABC"}
        got = kmers.shared_kmers(frozen, threads=2)
        for x, y in ("AB", "AC", "BC"):
            assert all(np.array_equal(g, z[f"km_pairs_{x}{y}_{c}"]) for g, c in zip(got.pairs[x, y], "abw"))
        for x in "ABC":
            tally = got.tallies["A", "B", "C"][x]
            assert np.array_equal(np.stack(list(tally.values()), 1), z[f"km_tally_{x}"])
    two = kmers.shared_kmers({x: libs[x] for x in labels[:2]})
    assert list(two.pairs) == [tuple(labels[:2])] and two.tallies == {}
    with pytest.raises(ValueError, match="exact"):
        kmers.shared_kmers(libs, scale=2)
