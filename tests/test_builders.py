import json
import random
from pathlib import Path

import builders
import pytest
from builders import (
    LTRQUEST_COLUMNS,
    ElementSpec,
    _parse,
    _read_fasta,
    build_scene,
    make_ltr_element,
    mutate,
    pav_toy,
    plant,
    plant_scaffolds,
    rand_seq,
    simulate,
    truth_subgenomes,
    write_fasta,
    write_gff3,
    write_ltrquest,
)

from synltr.io.annotation import load_annotation
from synltr.io.annotation.gff3 import load_gff3
from synltr.io.annotation.ltrquest import load_ltrquest
from synltr.io.elements import make_eid
from synltr.io.seqio import Fasta, n_runs, revcomp


def scene(seed=3):
    """A plus element, a minus host with unequal LTRs (250, 320) and a nested child, one more plus element."""
    child = ElementSpec(tsd="GAGAG", left=200, right=200, internal=500)
    host = ElementSpec(strand="-", tsd="CTCTC", left=250, right=320, internal=1800, nests=[(700, child)])
    return build_scene(random.Random(seed), [ElementSpec(), host, ElementSpec(strand="-", tsd="TTAGA")])


def key(e):
    return (e.seqid, e.start, e.end, e.strand, e.ltr_left, e.ltr_right, e.tsd, e.level)


def test_primitives_are_deterministic():
    assert rand_seq(50, random.Random(1)) == rand_seq(50, random.Random(1))
    s = rand_seq(500, random.Random(1))
    assert mutate(s, 0.0, random.Random(2)) == s and mutate(s, 1.0, random.Random(2)) != s
    a, b = make_ltr_element(random.Random(1), 100, 100, 200), make_ltr_element(random.Random(1), 100, 100, 200)
    assert a == b


def test_plant_places_the_element_after_the_tsd():
    g = "A" * 10 + "C" * 10
    out = plant(g, "GGGG", "TT", 10)
    assert out == "A" * 10 + "TT" + "GGGG" + "TT" + "C" * 10
    assert out[10 + 2:10 + 2 + 4] == "GGGG"  # 0-based pos + len(tsd), so 1-based pos + len(tsd) + 1


def test_unequal_ltrs_are_copies_of_each_other():
    seq, nl, nr = make_ltr_element(random.Random(1), 100, 140, 300)
    assert (nl, nr, len(seq)) == (100, 140, 540) and seq[:100] == seq[-140:-40]
    seq, nl, nr = make_ltr_element(random.Random(1), 140, 100, 300)
    assert seq[:100] == seq[-100:]


def test_scene_coordinates_match_the_sequence():
    records, els = scene()
    chrom = records["chr1"]
    assert len(els) == 4 and [e.level for e in els] == [0, 0, 1, 0]
    for e in els:
        s, t = e.start, e.end
        assert chrom[s - 1 - len(e.tsd):s - 1] == e.tsd == chrom[t:t + len(e.tsd)]  # TSDs flank the 1-based span
        left, right = chrom[e.ltr_left[0] - 1:e.ltr_left[1]], chrom[e.ltr_right[0] - 1:e.ltr_right[1]]
        n = min(len(left), len(right))
        assert left[:n] == right[:n] and e.ltr_left[0] == s and e.ltr_right[1] == t
    host = next(e for e in els if e.strand == "-" and e.level == 0 and e.ltr_left_len == 250)
    assert (host.ltr_left_len, host.ltr_right_len) == (250, 320)
    child = next(e for e in els if e.level == 1)
    assert child.parent == host.eid and host.start < child.start and child.end < host.end


def test_fasta_reads_back_through_fasta(tmp_path):
    records, _ = scene()
    path = write_fasta(tmp_path / "sub" / "g.fa", records)
    with Fasta(path) as fa:
        assert fa.lengths() == {"chr1": len(records["chr1"])}
        assert fa.fetch("chr1", 100, 400) == records["chr1"][100:400]


def test_gff3_reads_back_with_the_same_coordinates(tmp_path):
    _, els = scene()
    path = write_gff3(tmp_path / "g.gff3", els)
    got, skipped = load_gff3(path, "G")
    assert not skipped and sorted(key(e)[:6] for e in got) == sorted(key(e)[:6] for e in els)
    assert sorted(e.tsd for e in got) == sorted(e.tsd for e in els)
    assert {e.eid for e in got} == {e.eid for e in els} and all(e.k == 0.01 for e in got)
    table = load_annotation(str(path), "G")
    assert [(e.start, e.end, e.parent, e.level) for e in table] == [(e.start, e.end, e.parent, e.level) for e in els]
    assert all(e.ltr_left_len is None for e in got)


def test_gff3_lengths_are_the_ltrs_own(tmp_path):
    """With `lengths` the GFF3 carries ltr5_len and ltr3_len, read back as the LTRs' own lengths."""
    _, els = scene()
    got, _ = load_gff3(write_gff3(tmp_path / "g.gff3", els, lengths=True), "G")
    own = sorted((e.start, e.ltr_left_len, e.ltr_right_len) for e in els)
    assert sorted((e.start, e.ltr_left_len, e.ltr_right_len) for e in got) == own and own[0][1] > 0


@pytest.mark.parametrize("clean", ["none", "cut", "masked"])
def test_ltrquest_headers_and_record_sizes(tmp_path, clean):
    records, els = scene()
    paths = write_ltrquest(str(tmp_path / "X"), els, records, clean=clean)
    assert [p.name for p in paths] == ["X_depth0_clean_ltr.tsv", "X_depth1_clean_ltr.tsv"]
    fa = {}
    for p in paths:
        with Fasta(p.with_suffix(".fa")) as f:
            fa.update({n: f.fetch(n, 0, L) for n, L in f.lengths().items()})
    host = next(e for e in els if e.ltr_left_len == 250)
    child = next(e for e in els if e.level == 1)
    span, nested = host.end - host.start + 1, child.end - child.start + 1
    rec = fa[f"{host.seqid}:{host.start}-{host.end}#LTR/Gypsy/Ale"]
    assert host.strand == "-"  # the record is the reverse complement of the cleaned genome span
    whole = records["chr1"][host.start - 1:host.end]
    a = child.start - host.start
    if clean == "none":
        assert rec == revcomp(whole)
    elif clean == "masked":
        assert rec == revcomp(whole[:a] + "N" * nested + whole[a + nested:])
        assert len(rec) == span and rec.count("N") == nested and n_runs(rec) == [(span - a - nested, span - a)]
    else:
        assert len(rec) == span - nested and "N" not in rec
        assert rec == revcomp(whole[:a] + whole[a + nested:])
    # a plus-strand record is the genome span itself, whole
    assert child.strand == "+"
    assert fa[f"{child.seqid}:{child.start}-{child.end}#LTR/Gypsy/Ale"] == records["chr1"][child.start - 1:child.end]


@pytest.mark.parametrize("clean", ["none", "cut", "masked"])
def test_ltrquest_reads_back_with_the_same_coordinates(tmp_path, clean):
    records, els = scene()
    write_ltrquest(str(tmp_path / "X"), els, records, clean=clean)
    got, skipped = load_ltrquest(str(tmp_path / "X"), "G")
    assert not skipped and len(got) == len(els)
    by_eid = {e.eid: e for e in els}
    for e in got:
        w = by_eid[e.eid]
        assert (e.seqid, e.start, e.end, e.strand) == (w.seqid, w.start, w.end, w.strand)
        assert (e.ltr_left, e.ltr_right, e.tsd, e.k) == (w.ltr_left, w.ltr_right, w.tsd, w.k)
        assert (e.ltr_left_len, e.ltr_right_len) == (w.ltr_left_len, w.ltr_right_len)
        assert (e.superfamily, e.clade, e.family) == ("Gypsy", "Ale", "fam1")
    table = load_annotation(str(tmp_path / "X"), "G")
    assert [(e.start, e.end, e.level) for e in table] == [(e.start, e.end, e.level) for e in els]


def test_ltrquest_cut_table_keeps_span_offsets(tmp_path):
    records, els = scene()
    write_ltrquest(str(tmp_path / "X"), els, records, clean="cut")
    got, skipped = load_ltrquest(str(tmp_path / "X"), "G")
    assert not skipped and len(got) == len(els)
    by_eid = {e.eid: e for e in els}
    for e in got:  # the planted genome coordinates, with no shift for the shortened host record
        w = by_eid[e.eid]
        assert (e.start, e.end, e.ltr_left, e.ltr_right) == (w.start, w.end, w.ltr_left, w.ltr_right)
    host = next(e for e in els if e.ltr_left_len == 250)
    child = next(e for e in els if e.level == 1)
    row = next(r for r in (tmp_path / "X_depth0_clean_ltr.tsv").read_text().splitlines() if f":{host.start}-" in r)
    col = dict(zip(LTRQUEST_COLUMNS, row.split("\t"), strict=True))
    span = host.end - host.start + 1
    assert int(col["seq_len"]) == span - (child.end - child.start + 1)  # seq_len is the shortened record
    assert int(col["ltr3_end"]) == span  # span offsets, not record offsets
    assert (int(col["ltr5_len"]), int(col["ltr3_len"])) == (250, 320)  # unequal LTRs, host on the minus strand
    assert col["strand"] == "-"


def test_ltrquest_rejects_an_unknown_clean_kind(tmp_path):
    records, els = scene()
    with pytest.raises(ValueError, match="clean"):
        write_ltrquest(str(tmp_path / "X"), els, records, clean="trimmed")


def test_pav_toy(tmp_path):
    toy = pav_toy(tmp_path)
    assert len(toy.elements) == len(toy.expect) and toy.a != toy.b
    with Fasta(toy.fasta_a) as fa, Fasta(toy.fasta_b) as fb:
        assert fa.fetch("chr1", 0, len(toy.a)) == toy.a and fb.fetch("chr1", 0, len(toy.b)) == toy.b
        for e in toy.elements:  # the 1-based closed element is flanked by its TSD in A
            assert fa.fetch("chr1", e.start - 1 - len(e.tsd), e.start - 1) == e.tsd
            assert fa.fetch("chr1", e.end, e.end + len(e.tsd)) == e.tsd
    got = load_annotation(str(toy.gff3), "toy")
    assert [(e.start, e.end, e.parent) for e in got] == [(e.start, e.end, e.parent) for e in toy.elements]


# ---------------------------------------------------------------- the genotype toy cases (case, build_pair)

ELEM_LEN = {"e": 4800, "e1": 4800, "e2": 4800, "h": 7810, "c": 3000}
LTR_LEN = {"e": 400, "e1": 400, "e2": 400, "h": 400, "c": 300}
HOST = ("nested_E", "host_E", "host_S", "host_G")
STATES = {
    "P": {"e": "P"}, "P_nestB": {"e": "P"}, "S": {"e": "S"}, "S_notsd": {"e": "S"}, "T": {"e": "T"},
    "E": {"e": "E"}, "E_mm1": {"e": "E"}, "E_mm2": {"e": "E"}, "D_left": {"e": "D"}, "D_both": {"e": "D"},
    "O": {"e": "O"}, "G": {"e": "G"}, "G_near": {"e": "G"}, "G_far": {"e": "E"}, "U_flank": {"e": "U"},
    "inv_site": {"e": "U"}, "nested_E": {"h": "P", "c": "E"}, "host_E": {"h": "E", "c": "E"},
    "host_S": {"h": "S", "c": "D"}, "host_G": {"h": "G", "c": "G"}, "adjacent_E": {"e1": "E", "e2": "E"},
}
FLAGS = {("P_nestB", "e"): {"nestB"}, ("S", "e"): {"solo_tsd"}, ("E", "e"): {"tsd_match"},
         ("E_mm1", "e"): {"tsd_match"}, ("host_E", "c"): {"host_empty"}, ("host_G", "c"): {"inherit"},
         ("adjacent_E", "e1"): {"adjacent"}, ("adjacent_E", "e2"): {"adjacent"}}
NOFLAGS = {("S_notsd", "e"): {"solo_tsd"}, ("E_mm2", "e"): {"tsd_match"}}
TIERS = {("E", "e"): "high", ("E_mm2", "e"): "medium"}
REASONS = {("U_flank", "e"): "flank_unaligned", ("inv_site", "e"): "inv_site"}


def _hamming(x, y):
    return sum(p != q for p, q in zip(x, y))


def _expected_b(name, a):
    """B built from A's pieces, for the cases without random filler or mutation."""
    L, T, R = a[:3000], a[3000:3005], a[-3000:]
    e = a[3005:-3005]  # the element (host for host cases)
    return {
        "S": L + T + e[:400] + T + R,
        "S_notsd": L + T + e[:400] + R,
        "T": L + T + e[:400] + e[400:2000] + e[-400:] + T + R,
        "E": L + T + R,
        "D_left": L[:-500] + R,
        "D_both": L[:-300] + R[400:],
        "G": L + T + "N" * 100 + R,
        "G_near": L[:-70] + "N" * 50 + L[-20:] + T + R,
        "G_far": L[:-550] + "N" * 50 + L[-500:] + T + R,
        "inv_site": L + T + revcomp(e) + T + R,
        "nested_E": L + T + e[:2405] + e[5410:] + T + R,  # h without c and one GAGAG
        "host_E": L + T + R,
        "host_S": L + T + e[:400] + T + R,
        "host_G": L + T + "N" * 100 + R,
        "adjacent_E": L + T + "GAGAG" + R,
    }.get(name)


def test_make_element():
    seq, ll, rl = builders.make_element(random.Random(1), 400, 4000)
    assert (len(seq), ll, rl) == (4800, 400, 400) and seq[:400] == seq[-400:]


def test_cases_list():
    assert builders.CASES == list(STATES) and len(builders.CASES) == 21
    with pytest.raises(ValueError):
        builders.case("nope", random.Random(1))


@pytest.mark.parametrize("name", builders.CASES)
def test_toy_case_elements(name):
    c = builders.case(name, random.Random(7))
    by_name = {e.name: e for e in c.elements}
    assert set(by_name) == set(STATES[name]) and len(by_name) == len(c.elements)
    for e in c.elements:
        assert (e.genome, e.seqid, e.strand, e.k) == ("toy", "a", "+", 0.0)
        assert e.eid == make_eid("toy", "a", e.start, e.end)
        assert e.end - e.start + 1 == ELEM_LEN[e.name]
        n = LTR_LEN[e.name]
        assert e.ltr_left == (e.start, e.start + n - 1) and e.ltr_right == (e.end - n + 1, e.end)
        assert c.a[e.start - 1:e.start - 1 + n] == c.a[e.end - n:e.end]  # identical LTRs at the bounds
        assert e.tsd == ("GAGAG" if e.name in ("c", "e2") else "CTCTC")
        k = len(e.tsd)
        assert c.a[e.start - 1 - k:e.start - 1] == e.tsd == c.a[e.end:e.end + k]
    eid = {n: e.eid for n, e in by_name.items()}
    assert c.expect == {eid[n]: s for n, s in STATES[name].items()}
    assert c.flags == {eid[n]: f for (cn, n), f in FLAGS.items() if cn == name}
    assert c.noflags == {eid[n]: f for (cn, n), f in NOFLAGS.items() if cn == name}
    assert c.tiers == {eid[n]: t for (cn, n), t in TIERS.items() if cn == name}
    assert c.reasons == {eid[n]: r for (cn, n), r in REASONS.items() if cn == name}


@pytest.mark.parametrize("name", builders.CASES)
def test_toy_case_b(name):
    c = builders.case(name, random.Random(7))
    a, b = c.a, c.b
    want = _expected_b(name, a)
    if want is not None:
        assert b == want
    elif name == "P":
        assert len(b) == len(a) and 0.005 <= _hamming(a, b) / len(a) <= 0.015
    elif name in ("E_mm1", "E_mm2"):
        assert b[:3000] == a[:3000] and b[3005:] == a[-3000:] and len(b) == 6005
        assert _hamming(b[3000:3005], "CTCTC") == (1 if name == "E_mm1" else 2)
    elif name == "P_nestB":
        assert len(b) == len(a) + 3000 and b[:5405] == a[:5405] and b[8405:] == a[5405:]
    elif name == "O":
        assert len(b) == 8010 and b[:3005] == a[:3005] and b[-3005:] == a[-3005:]
    else:
        assert name == "U_flank" and len(b) == 6005 and b[:3005] == a[:3005] and b[3005:] != a[-3000:]


def test_toy_host_nesting():
    for name in HOST:
        c = builders.case(name, random.Random(8))
        h, ch = sorted(c.elements, key=lambda e: e.level)
        assert (h.name, h.level, h.parent) == ("h", 0, None)
        assert (ch.name, ch.level, ch.parent) == ("c", 1, h.eid)
        assert h.start + 2405 == ch.start and ch.end + 2405 == h.end
    for name in set(builders.CASES) - set(HOST):
        assert all(e.level == 0 and e.parent is None for e in builders.case(name, random.Random(8)).elements)


def _name_of(c, eid):
    return next(e.name for e in c.elements if e.eid == eid)


def test_build_pair_offsets():
    names = ["E", "host_E", "adjacent_E", "P_nestB"]
    pair = builders.build_pair(random.Random(5), names)
    rng = random.Random(5)  # build_pair draws the cases first, in order, then the spacers
    cases = [builders.case(n, rng) for n in names]
    sp = builders.SPACER
    assert len(pair.a) == sum(len(c.a) for c in cases) + (len(names) + 1) * sp
    assert len(pair.b) == sum(len(c.b) for c in cases) + (len(names) + 1) * sp
    oa = ob = sp
    for n, c in zip(names, cases):
        assert pair.a[oa - sp:oa] == pair.b[ob - sp:ob]  # same spacer before the case in A and B
        assert pair.a[oa:oa + len(c.a)] == c.a and pair.b[ob:ob + len(c.b)] == c.b
        new = {e.name: e for e in pair.elements if e.name.startswith(n + ".")}
        assert len(new) == len(c.elements)
        for e in c.elements:
            p = new[f"{n}.{e.name}"]
            assert (p.genome, p.seqid, p.tsd, p.k, p.level) == ("toy", "chr1", e.tsd, 0.0, e.level)
            assert (p.start, p.end) == (e.start + oa, e.end + oa)
            assert p.ltr_left == (e.ltr_left[0] + oa, e.ltr_left[1] + oa)
            assert p.ltr_right == (e.ltr_right[0] + oa, e.ltr_right[1] + oa)
            assert p.eid == make_eid("toy", "chr1", p.start, p.end)
            assert pair.a[p.start - 1:p.end] == c.a[e.start - 1:e.end]
            assert pair.expect[p.eid] == c.expect[e.eid]
            for d_pair, d_case in ((pair.flags, c.flags), (pair.noflags, c.noflags), (pair.tiers, c.tiers),
                                   (pair.reasons, c.reasons)):
                assert d_pair.get(p.eid) == d_case.get(e.eid)
            if e.parent is not None:
                assert p.parent == new[f"{n}.{_name_of(c, e.parent)}"].eid
        oa += len(c.a) + sp
        ob += len(c.b) + sp
    assert pair.a[-sp:] == pair.b[-sp:]
    assert len(pair.expect) == len(pair.elements) == sum(len(c.elements) for c in cases)


def test_build_pair_default_cases():
    pair = builders.build_pair(random.Random(6))
    assert {e.name.split(".")[0] for e in pair.elements} == set(builders.CASES) - {"U_flank"}
    assert all(e.seqid == "chr1" for e in pair.elements)
    assert sorted(e.start for e in pair.elements) == [e.start for e in pair.elements]


# ---- polyploid scenario ----
TINY = json.loads((Path(__file__).parent / "golden" / "nucleus_tiny" / "scenario.json").read_text())


@pytest.fixture(scope="module")
def poly(tmp_path_factory):
    d = tmp_path_factory.mktemp("poly")
    out, truth = simulate(TINY, d / "lib", d / "truth", seed=7)
    return d, out, truth


def test_simulate_is_sim_nucleus(poly):
    d, out, _ = poly
    gold = Path(__file__).parent / "golden" / "nucleus_tiny" / "lib"
    assert (out / "ltr.fa").read_bytes() == (gold / "ltr.fa").read_bytes()
    assert json.loads((d / "scenario_7.json").read_text()) == TINY


def test_simulate_failure_raises(tmp_path):
    """A scenario `sim nucleus` refuses (a merger naming no leaf of the tree: a usage error, exit 2) raises."""
    with pytest.raises(RuntimeError, match="exited 2"):
        simulate({**TINY, "mergers": [{"time": 1e5, "lineages": ["A", "Z"]}]}, tmp_path / "l", tmp_path / "t", 1)


def test_truth_subgenomes(poly):
    tsg = truth_subgenomes(poly[2])
    assert sorted(tsg) == [f"chr{i:03d}" for i in range(1, 7)]
    assert sorted(set(tsg.values())) == ["A", "C", "D"]


def test_plant_scaffolds_and_decoys(poly, tmp_path):
    out = poly[1]
    new, parent, decoys = plant_scaffolds(out, tmp_path / "planted", n_decoys=3, per_decoy=5, seed=2)
    assert decoys == ["decoy1", "decoy2", "decoy3"]
    assert sorted(parent.values()) == sorted(c for c in truth_subgenomes(poly[2]) for _ in range(3))
    fai = {ln.split("\t")[0]: int(ln.split("\t")[1]) for ln in (new / "genome.fai").read_text().splitlines()}
    assert set(fai) == set(truth_subgenomes(poly[2])) | set(parent) | set(decoys)
    recs = dict(_read_fasta(new / "ltr.fa"))
    for h in recs:
        chrom, a, b, _ = _parse(h)
        assert 1 <= a < b <= fai[chrom]
    assert sum(_parse(h)[0] in decoys for h in recs) == 15
    table = (new / "k2p.tsv").read_text().splitlines()
    assert [ln.split("\t")[0] for ln in table[1:]] == list(recs)
    placed = {h: q for h, q in _read_fasta(out / "ltr.fa")}
    moved = [h for h in recs if _parse(h)[0] in parent]
    assert moved and all(q in placed.values() for q in (recs[h] for h in moved))


def test_plant_scaffolds_is_deterministic(poly, tmp_path):
    a = plant_scaffolds(poly[1], tmp_path / "a", seed=5)
    b = plant_scaffolds(poly[1], tmp_path / "b", seed=5)
    assert (a[0] / "ltr.fa").read_bytes() == (b[0] / "ltr.fa").read_bytes() and a[1:] == b[1:]


def test_trio_is_ltrmix_synth(tmp_path):
    """builders.trio writes LTRmix's synthetic trio (its tests/synth.py make_trio, seed 3) byte for byte (LTRmix's
    sha256s, stored); the outgroup variant leaves R, S and D as they are and adds O, role outgroup, carrying old
    insertions only and its own young ones."""
    import hashlib

    import numpy as np

    from synltr.io.samples import read_samples

    with np.load(Path(__file__).resolve().parent / "golden" / "ltrmix_route.npz") as z:
        want = dict(zip(z["tr_files"].tolist(), z["tr_sha256"].tolist()))
    for d, make in ((tmp_path / "trio", builders.trio), (tmp_path / "og", builders.trio_outgroup)):
        sc = make(d, seed=3)
        assert {f: hashlib.sha256((d / f).read_bytes()).hexdigest() for f in want} == want
    sheet = read_samples(sc.sheet)
    assert [(s.label, s.role) for s in sheet] == [("R", None), ("S", None), ("D", None), ("O", "outgroup")]
    o = [t for t in sc.truth if "O" in t["carriers"]]
    assert {t["cls"] for t in o} == {"old", "pO"} and sum(t["cls"] == "old" for t in o) > 50
    assert (tmp_path / "og" / "O.fa").read_text().count(">") == len(o)
