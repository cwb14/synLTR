"""Locus catalog (ltrpav.loci): synthetic pair files in the CALL_COLUMNS format (gzip TSV, "." for none) of
genomes A and B (annotated) and C (unannotated), all on seqid chr1 unless stated; build_loci's links, states,
attributes, lids and QC rows."""
from __future__ import annotations

import gzip
import logging
import math
import random
from collections import Counter

import pytest

from synltr.calls.loci import GenomeState, build_loci
from synltr.calls.records import CALL_COLUMNS, tsv_field, write_calls
from synltr.calls.states import INFERRED_ABSENT, OBSERVED_ABSENT, PAIR_STATES
from synltr.io.elements import Element, make_eid
from synltr.io.samples import Sample

STATES = tuple(s for s in PAIR_STATES if s != "A")  # the states LTRpav's QC rows list

ORDER = {"chr1": 0, "chr2": 1}  # every genome's fai order unless a test says otherwise
NO_SITE = (None, None, None, None)


def el(genome, start, end, seqid="chr1", strand="+", ltr=100, **kw):
    """An element of `genome` (1-based inclusive) with an `ltr`-bp LTR at each end; `kw` replaces fields."""
    fields = {"genome": genome, "eid": make_eid(genome, seqid, start, end), "name": f"{genome}{start}",
              "seqid": seqid, "start": start, "end": end, "strand": strand,
              "ltr_left": (start, start + ltr - 1), "ltr_right": (end - ltr + 1, end), "tsd": "ACGTA",
              "k": 0.01, "family": None, "superfamily": None, "clade": None, "confidence": None}
    return Element(**{**fields, **kw})


def call(e, b, state, site=None, tier=None, strand="+", flags=None, reason=None):
    """((e's genome, b), the CALL_COLUMNS row of e's call in b): `tier` defaults to high ("." for G and U);
    `site` = (seqid, start, end) in b or None."""
    if tier is None:
        tier = "." if state in ("G", "U") else "high"
    seqid, start, end = site if site else (None, None, None)
    row = {"eid": e.eid, "a_seqid": e.seqid, "a_start": e.start, "a_end": e.end, "parent": e.parent,
           "level": e.level, "state": state, "tier": tier, "flags": flags, "reason": reason, "b_seqid": seqid,
           "b_start": start, "b_end": end, "b_strand": strand if site else None, "wid": 0}
    return (e.genome, b), row


def build(d, elements, calls=(), labels="ABC", pairs=None, order=None, rng=None, skipped=None, absent=OBSERVED_ABSENT):
    """build_loci over genomes `labels` (sheet order; annotated: the keys of `elements`) with pair files in d.

    calls: call() results; an element with no call in a written pair gets U no_alignment there. pairs: the
    ordered pairs written (default: every annotated A with every other genome). order: seq_order (default
    ORDER for every genome). rng: shuffles the element lists, the rows of each pair file and call_paths.
    skipped: build_loci's skip counts. absent: build_loci's absent alleles.
    """
    d.mkdir(parents=True, exist_ok=True)
    elements = {g: list(v) for g, v in elements.items()}
    given: dict[tuple[str, str], dict[str, dict]] = {}
    for pair, row in calls:
        given.setdefault(pair, {})[row["eid"]] = row
    if pairs is None:
        pairs = [(a, b) for a in elements for b in labels if b != a]
    assert set(given) <= set(pairs), "a test call in a pair that is not written"
    paths = {}
    for a, b in pairs:
        rows = [given.get((a, b), {}).get(e.eid) or call(e, b, "U", reason="no_alignment")[1]
                for e in elements[a]]
        if rng is not None:
            rng.shuffle(rows)
        paths[(a, b)] = d / f"{a}__{b}.calls.tsv.gz"
        write_calls(paths[(a, b)], rows, {"engine": "assembly"})
    if rng is not None:
        for v in elements.values():
            rng.shuffle(v)
        items = list(paths.items())
        rng.shuffle(items)
        paths = dict(items)
    samples = [Sample(x, d / f"{x}.fa", f"{x}.gff3" if x in elements else None, None, None) for x in labels]
    seq_order = order if order is not None else {x: dict(ORDER) for x in labels}
    return build_loci(samples, elements, paths, seq_order, skipped=skipped, absent=absent)


def by_member(loci):
    """{first member eid (sheet order): locus}."""
    return {next(iter(x.members.values())): x for x in loci}


# ---------------------------------------------------------------- links, states and conflicts


def test_link_and_states(tmp_path):
    a1, b1 = el("A", 1001, 2000, k=0.02), el("B", 5001, 6000, strand="-", k=0.03)
    loci, _ = build(tmp_path, {"A": [a1], "B": [b1]}, [
        call(a1, "B", "P", ("chr1", 5011, 6005), strand="-"),  # reciprocal overlap 0.99: links a1 and b1
        call(a1, "C", "E", ("chr1", 3001, 3005), flags="tsd_match"),
        call(b1, "C", "E", ("chr1", 3002, 3006), strand="-", flags="tsd_match"),
    ])
    [x] = loci
    assert (x.lid, x.members, x.extra, x.flags) == ("L0000001", {"A": a1.eid, "B": b1.eid}, {}, [])
    assert x.discovery == ["A", "B"]
    assert list(x.states) == ["A", "B", "C"]
    assert x.states["A"] == GenomeState("P", "annot", "chr1", 1001, 2000, "+", 0.02, [])
    assert x.states["B"] == GenomeState("P", "annot", "chr1", 5001, 6000, "-", 0.03, [])
    # both C calls E/high: the first candidate in sheet order (a1's) gives the site and flags
    assert x.states["C"] == GenomeState("E", "high", "chr1", 3001, 3005, "+", None, ["tsd_match"])


def _site(state, start):
    return None if state in ("G", "U") else ("chr1", start, start + 4)


@pytest.mark.parametrize("a, b, want", [
    (("E", "high"), ("D", "high"), ("U", "high", ["conflict"], None)),  # absent vs missing
    (("P", "high"), ("E", "high"), ("U", "high", ["conflict"], None)),  # present vs absent
    (("T", "medium"), ("D", "medium"), ("U", "medium", ["conflict"], None)),  # present vs missing
    (("E", "high"), ("D", "medium"), ("E", "high", [], 3001)),  # the highest tier decides alone
    (("D", "low"), ("E", "medium"), ("E", "medium", [], 3101)),
    (("T", "low"), ("U", "."), ("T", "low", [], 3001)),
    (("P", "high"), ("T", "high"), ("P", "high", [], 3001)),  # one class, present: a tie goes to P
    (("T", "high"), ("P", "high"), ("P", "high", [], 3101)),  # ... at the first P candidate (b1's call)
    (("T", "low"), ("S", "low"), ("S", "low", [], 3101)),  # S > T
    (("O", "low"), ("D", "low"), ("D", "low", [], 3101)),  # one class, missing: D > O
    (("G", "."), ("U", "."), ("G", ".", [], None)),  # G > U (no site)
    (("U", "."), ("G", "."), ("G", ".", [], None)),
    (("E", "high"), ("E", "high"), ("E", "high", [], 3001)),
    (("E", "low"), ("E", "high"), ("E", "high", [], 3101)),
])
def test_conflict_rules(tmp_path, a, b, want):
    """Top-tier calls of different genotype classes (present P/S/T, absent E, missing D/O/G/U) conflict;
    within one class the most frequent state wins (a tie: P > S > T, D > O > G > U)."""
    a1, b1 = el("A", 1001, 2000), el("B", 5001, 6000)
    loci, _ = build(tmp_path, {"A": [a1], "B": [b1]}, [
        call(a1, "B", "P", ("chr1", 5001, 6000)),
        call(a1, "C", a[0], _site(a[0], 3001), tier=a[1]),
        call(b1, "C", b[0], _site(b[0], 3101), tier=b[1]),
    ])
    [x] = loci
    c = x.states["C"]
    assert (c.state, c.tier, c.flags, c.start) == want
    if c.start is None:
        assert (c.seqid, c.start, c.end, c.strand) == NO_SITE


def test_conflict_rules_first_candidate_in_sheet_order(tmp_path):
    a1, b1 = el("A", 1001, 2000), el("B", 5001, 6000)
    loci, _ = build(tmp_path, {"A": [a1], "B": [b1]}, [
        call(a1, "B", "P", ("chr1", 5001, 6000)),
        call(a1, "C", "E", _site("E", 3001)),
        call(b1, "C", "E", _site("E", 3101)),
    ], labels="BAC")
    [x] = loci
    assert (x.states["C"].start, x.discovery, list(x.members)) == (3101, ["B", "A"], ["B", "A"])


@pytest.mark.parametrize("state, conflict", [("E", True), ("D", True), ("O", True), ("U", False),
                                             ("G", False), ("T", False)])
def test_annotation_conflict(tmp_path, state, conflict):
    a2, b2 = el("A", 7001, 8000), el("B", 9001, 10000)
    site = None if state in ("G", "U") else ("chr1", 20001, 20100)  # T's image elsewhere links nothing
    loci, _ = build(tmp_path, {"A": [a2], "B": [b2]}, [
        call(b2, "A", "P", ("chr1", 7001, 8000)),  # links b2 and a2
        call(a2, "B", state, site),
    ])
    [x] = loci
    assert x.members == {"A": a2.eid, "B": b2.eid}
    assert x.states["B"] == GenomeState("P", "annot", "chr1", 9001, 10000, "+", 0.01,
                                        ["annot_conflict"] if conflict else [])
    assert x.states["A"].flags == []


def test_unannotated_genome(tmp_path):
    a1, b1 = el("A", 1001, 2000), el("B", 1001, 2000)
    loci, _ = build(tmp_path, {"A": [a1], "B": [b1]}, [
        call(a1, "C", "P", ("chr1", 4001, 5000)),
        call(b1, "C", "P", ("chr1", 4001, 5000)),  # the same C copy: C holds no element, so nothing links
    ])
    assert len(loci) == 2
    for x in loci:
        assert "C" not in x.members and "C" not in x.extra and "C" not in x.discovery
        assert x.states["C"] == GenomeState("P", "high", "chr1", 4001, 5000, "+", None, [])


def test_shared_seqid_names(tmp_path):
    a, b = el("A", 1000, 2000), el("B", 1000, 2000)  # one seqid name and span in two genomes
    loci, _ = build(tmp_path, {"A": [a], "B": [b]}, [call(a, "B", "P", ("chr1", 50001, 51000))])
    assert [x.members for x in loci] == [{"A": a.eid}, {"B": b.eid}]
    assert loci[0].states["B"] == GenomeState("P", "high", "chr1", 50001, 51000, "+", None, [])
    assert loci[1].states["A"] == GenomeState("U", ".", None, None, None, None, None, [], "no_alignment")


def _world():
    """Five loci on chr1 and chr2: {a1, a1x, b1} multi-member with a conflict in C, {a2, b2} with an
    annotation conflict in B, {a3} (A only), {b3} (B only, no A site), {b4} (B only, an A empty site)."""
    a1 = el("A", 1001, 2000)
    a1x = el("A", 1011, 2000, level=1, parent=a1.eid)  # a second annotation of a1's insertion
    a2, a3 = el("A", 7001, 8000), el("A", 501, 1500, seqid="chr2")
    b1, b2 = el("B", 5001, 6000), el("B", 9001, 10000)
    b3, b4 = el("B", 3001, 4000, seqid="chr2"), el("B", 100, 1000)
    elements = {"A": [a1, a1x, a2, a3], "B": [b1, b2, b3, b4]}
    calls = [
        call(a1, "B", "P", ("chr1", 5001, 6000)),
        call(a1x, "B", "P", ("chr1", 5011, 6000)),
        call(b1, "A", "P", ("chr1", 1001, 2000)),
        call(a1, "C", "E", ("chr1", 3001, 3005)),
        call(a1x, "C", "D", ("chr1", 3001, 3010)),  # an extra member's call is a candidate too
        call(b2, "A", "P", ("chr1", 7001, 8000)),
        call(a2, "B", "E", ("chr1", 9001, 9005)),
        call(a3, "C", "P", ("chr2", 101, 1100)),
        call(b3, "C", "T", ("chr2", 4001, 4500), tier="low"),
        call(b4, "A", "E", ("chr1", 30001, 30005)),
    ]
    return elements, calls


def test_lid_deterministic(tmp_path):
    elements, calls = _world()
    want = build(tmp_path / "ref", elements, calls)
    assert [x.lid for x in want[0]] == [f"L{i:07d}" for i in range(1, 6)]
    for seed in (1, 2, 3):
        assert build(tmp_path / str(seed), elements, calls, rng=random.Random(seed)) == want


def test_qc_counters(tmp_path):
    elements, calls = _world()
    _, qc = build(tmp_path, elements, calls)
    assert qc[:4] == [("all", "n_loci", 5), ("all", "n_conflict", 1), ("all", "n_multi_member", 1),
                      ("all", "n_annot_conflict", 1)]
    assert all(scope != "all" for scope, _, _ in qc[4:])


# ---------------------------------------------------------------- members, states and attributes


def test_multi_member(tmp_path):
    elements, calls = _world()
    a1, a1x, b1 = elements["A"][0], elements["A"][1], elements["B"][0]
    x = by_member(build(tmp_path, elements, calls)[0])[a1.eid]
    assert (x.members, x.extra, x.flags) == ({"A": a1.eid, "B": b1.eid}, {"A": [a1x.eid]}, ["multi_member"])
    assert x.discovery == ["A", "B"]
    assert x.states["A"] == GenomeState("P", "annot", "chr1", 1001, 2000, "+", 0.01, [])  # the first member
    assert x.states["C"] == GenomeState("U", "high", None, None, None, None, None, ["conflict"])


def test_no_call(tmp_path, caplog):
    caplog.set_level(logging.WARNING, logger="ltrpav.loci")
    a1, a2, b2 = el("A", 1001, 2000), el("A", 7001, 8000), el("B", 9001, 10000)
    loci, _ = build(tmp_path, {"A": [a1, a2], "B": [b2]}, [
        call(b2, "A", "P", ("chr1", 7001, 8000)),
        call(b2, "C", "E", ("chr1", 501, 505)),
    ], pairs=[("A", "B"), ("B", "A"), ("B", "C")])  # A -> C was not run
    x = by_member(loci)
    assert x[a1.eid].states["C"] == GenomeState("U", ".", None, None, None, None, None, ["no_call"])
    assert x[a2.eid].states["C"] == GenomeState("E", "high", "chr1", 501, 505, "+", None, [])
    assert "no calls for 1 pair(s) (A__C)" in caplog.text


def test_u_reason(tmp_path):
    a1 = el("A", 1001, 2000)
    loci, _ = build(tmp_path, {"A": [a1]}, [call(a1, "C", "U", reason="flank_unaligned")], labels="AC")
    assert loci[0].states["C"] == GenomeState("U", ".", None, None, None, None, None, [], "flank_unaligned")


@pytest.mark.parametrize("strand, b_strand, want", [("+", "+", "+"), ("+", "-", "-"), ("-", "+", "-"),
                                                    ("-", "-", "+"), (".", "-", ".")])
def test_projected_strand(tmp_path, strand, b_strand, want):
    """A projected copy's strand is the element's own, turned when the B site is reversed (b_strand -)."""
    a1 = el("A", 1001, 2000, strand=strand)
    loci, _ = build(tmp_path, {"A": [a1]}, [call(a1, "C", "P", ("chr1", 4001, 5000), strand=b_strand)],
                    labels="AC")
    assert loci[0].states["C"].strand == want


def _three(tmp_path, a1, b1, d1, labels="ABCD", calls=()):
    """a1, b1, d1 (genomes A, B, D: all annotated) linked by a1's P calls onto b1 and d1; `calls` added."""
    calls = [call(a1, "B", "P", ("chr1", b1.start, b1.end)), call(a1, "D", "P", ("chr1", d1.start, d1.end)),
             *calls]
    [x] = build(tmp_path, {"A": [a1], "B": [b1], "D": [d1]}, calls, labels=labels)[0]
    assert x.members == {"A": a1.eid, "B": b1.eid, "D": d1.eid}
    return x


@pytest.mark.parametrize("states, want", [
    ("TTP", ("T", 3001)),  # the most frequent state beats the tie order
    ("PTT", ("T", 3101)),  # ... at its first candidate
    ("DOO", ("O", 3101)),
])
def test_most_frequent_state(tmp_path, states, want):
    a1, b1, d1 = el("A", 1001, 2000), el("B", 5001, 6000), el("D", 2001, 3000)
    calls = [call(e, "C", st, _site(st, start))
             for e, st, start in zip((a1, b1, d1), states, (3001, 3101, 3201), strict=True)]
    x = _three(tmp_path, a1, b1, d1, calls=calls)
    assert (x.states["C"].state, x.states["C"].tier, x.states["C"].flags, x.states["C"].start) == (
        want[0], "high", [], want[1])


def test_majority_family(tmp_path):
    x = _three(tmp_path, el("A", 1001, 2000, family="Ji", superfamily="Copia"),
               el("B", 5001, 6000, family="Opie", superfamily=None, clade="Retrofit"),
               el("D", 2001, 3000, family="Opie", superfamily="Copia"))
    assert (x.family, x.superfamily, x.clade) == ("Opie", "Copia", "Retrofit")  # None never votes


@pytest.mark.parametrize("labels, family", [("ABC", "Ji"), ("BAC", "Opie")])
def test_majority_family_tie(tmp_path, labels, family):
    a1, b1 = el("A", 1001, 2000, family="Ji"), el("B", 5001, 6000, family="Opie")
    calls = [call(a1, "B", "P", ("chr1", 5001, 6000))]
    [x] = build(tmp_path, {"A": [a1], "B": [b1]}, calls, labels=labels)[0]
    assert x.family == family  # a tie goes to the first member in sheet order


def test_median_length(tmp_path):
    x = _three(tmp_path, el("A", 1001, 2000, ltr=100), el("B", 5001, 6200, ltr=150),
               el("D", 1, 5000, ltr=400))
    assert (x.length, x.ltr_len) == (1200, 150)


def test_median_length_even(tmp_path):
    a1 = el("A", 1001, 2000)  # 1000 bp, LTRs 100 and 100
    b1 = el("B", 5001, 6001, ltr_left=(5001, 5101), ltr_right=(5900, 6001))  # 1001 bp, LTRs 101 and 102
    [x] = build(tmp_path, {"A": [a1], "B": [b1]}, [call(a1, "B", "P", ("chr1", 5001, 6001))])[0]
    assert (x.length, x.ltr_len) == (1001, 101)  # medians 1000.5 and 100.75, rounded half up


def test_tsd_first_member_with_one(tmp_path):
    a1, b1 = el("A", 1001, 2000, tsd=None), el("B", 5001, 6000, tsd="TTAGC")
    [x] = build(tmp_path, {"A": [a1], "B": [b1]}, [call(a1, "B", "P", ("chr1", 5001, 6000))])[0]
    assert x.tsd == "TTAGC"


# ---------------------------------------------------------------- lids


@pytest.mark.parametrize("a_order, want", [
    ({"chr1": 0, "chr2": 1}, ["a1", "a2", "b4", "a3", "b3"]),
    ({"chr2": 0, "chr1": 1}, ["a3", "a1", "a2", "b4", "b3"]),  # fai order, not seqid names, sorts
])
def test_lid_order(tmp_path, a_order, want):
    """Loci sort by the first sheet genome with coordinates (b4 by its empty site in A, b3 by B), its fai
    seqid index, then start."""
    elements, calls = _world()
    a1, _, a2, a3 = elements["A"]
    _, _, b3, b4 = elements["B"]
    names = {a1.eid: "a1", a2.eid: "a2", a3.eid: "a3", b3.eid: "b3", b4.eid: "b4"}
    order = {"A": a_order, "B": dict(ORDER), "C": dict(ORDER)}
    loci, _ = build(tmp_path, elements, calls, order=order)
    assert [names[next(iter(x.members.values()))] for x in loci] == want
    assert [x.lid for x in loci] == [f"L{i:07d}" for i in range(1, 6)]


# ---------------------------------------------------------------- linking


@pytest.mark.parametrize("span, image, linked", [
    ((5001, 6000), (5001, 5800), True),  # 800 bp: 0.80 of the element, all of the call
    ((5001, 6000), (5001, 5790), False),  # 0.79 of the element
    ((5001, 5800), (5001, 6000), True),  # 0.80 of the call
    ((5001, 5790), (5001, 6000), False),  # 0.79 of the call
    ((5001, 6000), (5201, 6200), True),  # shifted: 800 bp, 0.80 of each
    ((5001, 6000), (5211, 6210), False),  # 790 bp, 0.79 of each
    ((5201, 6200), (5001, 6000), True),  # the element starts after the call: 0.80 of each
    ((5211, 6210), (5001, 6000), False),  # 0.79 of each
])
def test_reciprocal_overlap_threshold(tmp_path, span, image, linked):
    a1, b1 = el("A", 1001, 2000), el("B", *span)
    loci, _ = build(tmp_path, {"A": [a1], "B": [b1]}, [call(a1, "B", "P", ("chr1", *image))])
    assert len(loci) == (1 if linked else 2)


@pytest.mark.parametrize("state, linked", [("P", True), ("S", True), ("T", True), ("E", False), ("D", False),
                                           ("O", False)])
def test_link_states(tmp_path, state, linked):
    a1, b1 = el("A", 1001, 2000), el("B", 5001, 6000)
    loci, _ = build(tmp_path, {"A": [a1], "B": [b1]}, [call(a1, "B", state, ("chr1", 5001, 6000))])
    assert len(loci) == (1 if linked else 2)


def test_link_needs_same_seqid(tmp_path):
    a1, b1 = el("A", 1001, 2000), el("B", 5001, 6000)
    loci, _ = build(tmp_path, {"A": [a1], "B": [b1]}, [call(a1, "B", "P", ("chr2", 5001, 6000))])
    assert len(loci) == 2


def test_link_nested_b_element(tmp_path):
    """B elements of every level are candidates; reciprocal overlap keeps the host out."""
    a1 = el("A", 1001, 2000)
    host = el("B", 5001, 20000)
    inner = el("B", 8001, 9000, level=1, parent=host.eid)
    loci, _ = build(tmp_path, {"A": [a1], "B": [host, inner]}, [call(a1, "B", "P", ("chr1", 8001, 9000))])
    x = by_member(loci)[a1.eid]
    assert x.members == {"A": a1.eid, "B": inner.eid}
    assert len(loci) == 2


# ---------------------------------------------------------------- input checks


def _inputs(tmp_path, rows=None, header=CALL_COLUMNS):
    """(samples, elements, call_paths, seq_order, a1): A (a1 only) and C, one pair file A__C of `rows`
    (default: a1 U) under `header`."""
    a1 = el("A", 1001, 2000)
    path = tmp_path / "A__C.calls.tsv.gz"
    with gzip.open(path, "wt", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(header) + "\n")
        for r in [call(a1, "C", "U")[1]] if rows is None else rows:
            fh.write("\t".join(tsv_field(r.get(c)) for c in CALL_COLUMNS) + "\n")
    samples = [Sample("A", tmp_path / "A.fa", "A.gff3", None, None),
               Sample("C", tmp_path / "C.fa", None, None, None)]
    return samples, {"A": [a1]}, {("A", "C"): path}, {"A": dict(ORDER), "C": dict(ORDER)}, a1


def test_inputs_ok(tmp_path):
    samples, elements, paths, order, a1 = _inputs(tmp_path)
    loci, _ = build_loci(samples, elements, paths, order)
    assert [x.members for x in loci] == [{"A": a1.eid}]


def test_bad_header(tmp_path):
    with pytest.raises(ValueError, match="header"):
        build_loci(*_inputs(tmp_path, header=CALL_COLUMNS[:-1])[:4])


def test_call_of_unknown_element(tmp_path):
    rows = [call(el("A", 1001, 2000), "C", "U")[1], call(el("A", 5001, 6000), "C", "U")[1]]
    with pytest.raises(ValueError, match="not an element of A"):
        build_loci(*_inputs(tmp_path, rows=rows)[:4])


def test_element_without_call(tmp_path):
    with pytest.raises(ValueError, match="no call"):
        build_loci(*_inputs(tmp_path, rows=[])[:4])


def test_element_called_twice(tmp_path):
    row = call(el("A", 1001, 2000), "C", "U")[1]
    with pytest.raises(ValueError, match="more than once"):
        build_loci(*_inputs(tmp_path, rows=[row, row])[:4])


def test_b_seqid_not_in_fai(tmp_path):
    row = call(el("A", 1001, 2000), "C", "E", ("chrX", 1, 5))[1]
    with pytest.raises(ValueError, match="B seqid 'chrX' is not in C's seq_order"):
        build_loci(*_inputs(tmp_path, rows=[row])[:4])


@pytest.mark.parametrize("state, site, message", [
    ("P", None, "a P call without a B site"),
    ("E", ("chr1", 9, 5), r"B site chr1:9-5 is not two integers with 1 <= start <= end"),
    ("D", ("chr1", 0, 1), r"B site chr1:0-1 is not"),
    ("S", ("chr1", "x", 5), r"B site chr1:x-5 is not"),
])
def test_malformed_b_site(tmp_path, state, site, message):
    """A B site that is not 1 <= start <= end, or a present call without one, fails."""
    row = call(el("A", 1001, 2000), "C", state, site)[1]
    with pytest.raises(ValueError, match=message):
        build_loci(*_inputs(tmp_path, rows=[row])[:4])


@pytest.mark.parametrize("field, value", [("state", "Q"), ("tier", "top"), ("b_strand", "?")])
def test_unknown_state_tier_or_strand(tmp_path, field, value):
    row = {**call(el("A", 1001, 2000), "C", "E", ("chr1", 1, 5))[1], field: value}
    with pytest.raises(ValueError, match="unknown state, tier or b_strand"):
        build_loci(*_inputs(tmp_path, rows=[row])[:4])


def test_short_line(tmp_path):
    samples, elements, paths, order, _ = _inputs(tmp_path)
    path = paths[("A", "C")]
    with gzip.open(path, "wt") as fh:
        fh.write("\t".join(CALL_COLUMNS) + "\n" + "x\ty\n")
    with pytest.raises(ValueError, match=r"line 2 has 2 column\(s\); a pair file has 29"):
        build_loci(samples, elements, paths, order)


def test_bad_pairs_and_genomes(tmp_path):
    samples, elements, paths, order, a1 = _inputs(tmp_path)
    with pytest.raises(ValueError, match="genome not in the samples"):
        build_loci(samples, elements, {("A", "Z"): paths[("A", "C")]}, order)
    with pytest.raises(ValueError, match="a genome against itself"):
        build_loci(samples, elements, {("A", "A"): paths[("A", "C")]}, order)
    with pytest.raises(ValueError, match="duplicate sample label"):
        build_loci(samples + samples[:1], elements, paths, order)
    with pytest.raises(ValueError, match=r"elements of genome\(s\) not in the samples: Z"):
        build_loci(samples, {**elements, "Z": [a1]}, paths, order)
    with pytest.raises(ValueError, match="seq_order lacks genome"):
        build_loci(samples, elements, paths, {"A": dict(ORDER)})


def test_corrupt_pair_file(tmp_path):
    samples, elements, paths, order, _ = _inputs(tmp_path)
    path = paths[("A", "C")]
    path.write_bytes(path.read_bytes()[:-10])  # truncated: no end-of-stream marker
    with pytest.raises(ValueError, match=r"A__C\.calls\.tsv\.gz: corrupt pair file"):
        build_loci(samples, elements, paths, order)


def test_pair_from_unannotated_genome(tmp_path):
    samples, elements, paths, order, _ = _inputs(tmp_path)
    with pytest.raises(ValueError, match="no annotation"):
        build_loci(samples, elements, {("C", "A"): paths[("A", "C")]}, order)


def test_duplicate_element(tmp_path):
    samples, _, paths, order, a1 = _inputs(tmp_path)
    with pytest.raises(ValueError, match="duplicate element id"):
        build_loci(samples, {"A": [a1, a1]}, paths, order)


def test_element_seqid_not_in_fai(tmp_path):
    samples, elements, paths, _, _ = _inputs(tmp_path)
    with pytest.raises(ValueError, match="seqid 'chr1' is not in A's seq_order"):
        build_loci(samples, elements, paths, {"A": {"chr2": 0}, "C": dict(ORDER)})


# ---------------------------------------------------------------- QC metrics


def qc_map(qc):
    """{(scope, metric): value} of QC rows (each key once)."""
    out = {(s, m): v for s, m, v in qc}
    assert len(out) == len(qc)
    return out


def test_qc_state_fractions(tmp_path):
    a = [el("A", 1001 + 10000 * i, 2000 + 10000 * i) for i in range(4)]
    _, qc = build(tmp_path, {"A": a}, [
        call(a[0], "C", "P", ("chr1", 5001, 6000)),  # tier high (the default)
        call(a[1], "C", "E", ("chr1", 7001, 7005), tier="medium"),
        call(a[2], "C", "G"),
    ], labels="AC")  # a[3]: U
    q = qc_map(qc)
    want = {"P": 0.25, "S": 0, "T": 0, "E": 0.25, "D": 0, "O": 0, "G": 0.25, "U": 0.25}
    assert {s: q[("A>C", f"frac_{s}")] for s in STATES} == want
    assert (q[("A>C", "resolved_frac")], q[("A>C", "high_frac")]) == (0.5, 0.25)
    assert ("A>C", "annot_agree") not in q and ("A>C", "e_in_annot") not in q  # C has no annotation
    assert math.isnan(q[("A", "gap_rate")]) and q[("C", "gap_rate")] == 0.25  # no calls into A


def test_qc_annotation_agreement(tmp_path):
    """annot_agree over the calls with a B element at the B site, e_in_annot over the E calls whose site meets a
    B element; an element in the locus of an ancestor of the called element does not count (a nested element's
    site lies inside its host's ortholog), and one spanning the site is not at it (annot_agree)."""
    h = el("A", 1001, 9000)
    c = el("A", 3001, 4000, level=1, parent=h.eid)
    c2 = el("A", 3201, 3500, level=2, parent=c.eid)
    a2, a3, a4, a5, a6 = (el("A", s, s + 999) for s in (20001, 30001, 40001, 50001, 60001))
    hb, b2, b5 = el("B", 101001, 108000), el("B", 120001, 121000), el("B", 140001, 141000)
    _, qc = build(tmp_path, {"A": [h, c, c2, a2, a3, a4, a5, a6], "B": [hb, b2, b5]}, [
        call(h, "B", "P", ("chr1", 101001, 108000)),  # links h and hb; agrees
        call(c, "B", "E", ("chr1", 103001, 103005)),  # inside hb, its host's ortholog: not counted
        call(c2, "B", "E", ("chr1", 103201, 103205)),  # inside hb, its grandparent's ortholog: not counted
        call(a2, "B", "E", ("chr1", 120501, 120505)),  # deep inside b2 (it spans the site): an E in annotation
        call(a3, "B", "D", ("chr1", 140000, 140001)),  # the junction meets b5 by 1 bp: b5 at the site, disagrees
        call(a4, "B", "O", ("chr1", 120001, 121000)),  # O has no site
        call(a5, "B", "E", ("chr1", 150001, 150005)),  # meets no B element
        call(a6, "B", "E", ("chr1", 140998, 141002)),  # overlaps b5's end: b5 at the site, disagrees
    ], labels="AB")
    q = qc_map(qc)
    assert q[("A>B", "annot_agree")] == 1 / 3  # h of h, a3, a6 (b2 spans a2's site)
    assert q[("A>B", "e_in_annot")] == 0.4  # a2 and a6 of the 5 E calls
    # B's calls in A are all U (no site): nothing to agree with, no E call
    assert math.isnan(q[("B>A", "annot_agree")]) and q[("B>A", "e_in_annot")] == 0


def test_qc_e_in_annot_exempts_the_hosts_present_image(tmp_path):
    """e_in_annot: a nested E inside a B element hb that its host's P call places (the host's B image overlaps
    hb), though hb is not in the host's locus (reciprocal overlap under COV_PRESENT: B's annotation gives the
    host other boundaries), is not counted: the host's present call is direct evidence that hb holds the host's
    ortholog. An E inside an unrelated B element is still counted, a nested one included. annot_agree exempts
    only the locus of an ancestor."""
    h = el("A", 1001, 9000)
    c = el("A", 3001, 4000, level=1, parent=h.eid)
    c2 = el("A", 5001, 6000, level=1, parent=h.eid)
    a2 = el("A", 20001, 21000)
    hb = el("B", 101001, 104000)  # 3 kb of the host's 8-kb image: overlap 0.375
    b2, b3 = el("B", 120001, 121000), el("B", 130001, 131000)
    loci, qc = build(tmp_path, {"A": [h, c, c2, a2], "B": [hb, b2, b3]}, [
        call(h, "B", "P", ("chr1", 101001, 108000)),
        call(c, "B", "E", ("chr1", 103001, 103005)),  # inside hb, which the host's image overlaps: exempt
        call(c2, "B", "E", ("chr1", 130501, 130505)),  # inside b3, away from the host's image: counted
        call(a2, "B", "E", ("chr1", 120501, 120505)),  # inside b2, no ancestor: counted
    ], labels="AB")
    host = next(x for x in loci if x.members.get("A") == h.eid)
    assert "B" not in host.members  # hb is not in the host's locus
    q = qc_map(qc)
    assert q[("A>B", "e_in_annot")] == 2 / 3  # c2 and a2 of the 3 E calls
    assert q[("A>B", "annot_agree")] == 1  # h; hb, b3 and b2 span the E sites of c, c2 and a2


def test_qc_e_in_annot_needs_a_present_host(tmp_path):
    """Only a present host (P, S or T) exempts through its B image: a host called E inside hb (its empty site
    in an element of B's annotation) leaves its own E and its nested element's inherited E counted."""
    h = el("A", 1001, 9000)
    c = el("A", 3001, 4000, level=1, parent=h.eid)
    hb = el("B", 101001, 104000)
    _, qc = build(tmp_path, {"A": [h, c], "B": [hb]}, [
        call(h, "B", "E", ("chr1", 103001, 103005)),
        call(c, "B", "E", ("chr1", 103001, 103005), flags="host_empty"),
    ], labels="AB")
    assert qc_map(qc)[("A>B", "e_in_annot")] == 1.0


def test_qc_site_inside_long_element(tmp_path):
    """A site deep inside a long B element, past a nested element that starts after it: it meets the long one
    (e_in_annot), which spans it (no B element at the site: annot_agree has nothing to count)."""
    a1 = el("A", 1001, 2000)
    host = el("B", 100001, 200000)
    inner = el("B", 150001, 151000, level=1, parent=host.eid)
    _, qc = build(tmp_path, {"A": [a1], "B": [host, inner]}, [call(a1, "B", "E", ("chr1", 180001, 180005))],
                  labels="AB")
    q = qc_map(qc)
    assert math.isnan(q[("A>B", "annot_agree")]) and q[("A>B", "e_in_annot")] == 1


def test_qc_annot_agree_needs_a_b_element_at_the_site(tmp_path):
    """annot_agree counts a call only where a B element is at its site. An
    element reaching more than w = 2 x 5 bp (two TSD lengths) past both site ends spans it (a host annotated
    only in B) and is set aside; one ending within w of the site, or inside it, is at it. e_in_annot counts
    the E calls overlapping any B element (spanning ones included), partly overlapping ones too."""
    a1, a2, a3, a6, a7, a8 = (el("A", s, s + 999) for s in (1001, 11001, 21001, 31001, 41001, 51001))
    a4, a5 = el("A", 61001, 62000), el("A", 71001, 72000)
    bh = el("B", 1001, 9000)  # spans a1's site by kb on both sides
    b2 = el("B", 19995, 21000)  # starts 6 bp before a2's site: at the site
    b3 = el("B", 29991, 31000)  # starts w = 10 bp before a3's site: at the site
    b6 = el("B", 32990, 34000)  # starts w + 1 bp before a6's site and ends far after it: spans it
    b8 = el("B", 90003, 91000)  # starts inside a8's site: at the site, a partial overlap
    b4, bh4 = el("B", 40001, 41000), el("B", 35001, 50000)  # a4's ortholog and a host spanning it
    bh5 = el("B", 55001, 70000)  # spans a5's present site, which has no B element at it
    _, qc = build(tmp_path, {"A": [a1, a2, a3, a6, a7, a8, a4, a5], "B": [bh, b2, b3, b6, b8, b4, bh4, bh5]}, [
        call(a1, "B", "E", ("chr1", 5001, 5005)),
        call(a2, "B", "E", ("chr1", 20001, 20005)),
        call(a3, "B", "E", ("chr1", 30001, 30005)),
        call(a6, "B", "E", ("chr1", 33001, 33005)),
        call(a7, "B", "E", ("chr1", 80001, 80005)),  # meets no B element
        call(a8, "B", "E", ("chr1", 90001, 90005)),
        call(a4, "B", "P", ("chr1", 40001, 41000)),
        call(a5, "B", "T", ("chr1", 60001, 61000)),
    ], labels="AB")
    q = qc_map(qc)
    assert q[("A>B", "annot_agree")] == 1 / 4  # a4 of a2, a3, a8, a4
    assert q[("A>B", "e_in_annot")] == 5 / 6  # a1, a2, a3, a6, a8 of the 6 E calls (not a7)


def test_qc_reciprocity(tmp_path):
    """Reciprocity is the share of links resolved both ways (both calls P, S, T or
    E) whose two calls are both present; link_symmetry, the share of links seen both ways, also drops a link
    whose other call is unresolved or misses the reciprocal overlap (a boundary difference)."""
    a1, a2, a3, a4, a5 = (el("A", s, s + 999) for s in (1001, 5001, 9001, 13001, 17001))
    b1, b2, b3, b4, b5 = (el("B", s, s + 999) for s in (21001, 25001, 29001, 33001, 37001))
    _, qc = build(tmp_path, {"A": [a1, a2, a3, a4, a5], "B": [b1, b2, b3, b4, b5]}, [
        call(a1, "B", "P", ("chr1", 21001, 22000)),  # a1-b1: seen both ways, both present: agree
        call(b1, "A", "S", ("chr1", 1001, 2000)),
        call(a2, "B", "T", ("chr1", 25001, 26000)),  # a2-b2: seen from A only, b2 calls A E: disagree
        call(b2, "A", "E", ("chr1", 5001, 5005)),
        call(b3, "A", "P", ("chr1", 9001, 10000)),  # a3-b3: seen from B only, a3 calls B U: not resolved
        call(a4, "B", "P", ("chr1", 33001, 34000)),  # a4-b4: seen from A only, b4 present in A at an image
        call(b4, "A", "P", ("chr1", 13501, 14400)),  # overlapping a4 by 500 of 900 bp: agree
        call(a5, "B", "P", ("chr1", 37001, 38000)),  # a5-b5: b5 calls A D: not resolved
        call(b5, "A", "D", ("chr1", 17000, 17001)),
    ])
    q = qc_map(qc)
    assert q[("A~B", "reciprocity")] == pytest.approx(2 / 3)  # a1-b1 and a4-b4 of a1-b1, a2-b2, a4-b4
    assert q[("A~B", "link_symmetry")] == pytest.approx(1 / 5)  # a1-b1 of the 5 links
    assert [k for k in q if k[1] in ("reciprocity", "link_symmetry")] == [("A~B", "reciprocity"),
                                                                         ("A~B", "link_symmetry")]  # C: none


def test_qc_reciprocity_needs_both_directions(tmp_path):
    elements = {"A": [el("A", 1001, 2000)], "B": [el("B", 1001, 2000)]}
    q = qc_map(build(tmp_path / "both", elements, labels="AB")[1])
    assert math.isnan(q[("A~B", "reciprocity")]) and math.isnan(q[("A~B", "link_symmetry")])  # no link
    q = qc_map(build(tmp_path / "one", elements, labels="AB", pairs=[("A", "B")])[1])
    assert ("A~B", "reciprocity") not in q and ("A~B", "link_symmetry") not in q
    q = qc_map(build(tmp_path / "unresolved", elements, [call(elements["A"][0], "B", "P", ("chr1", 1001, 2000))],
                     labels="AB")[1])
    assert math.isnan(q[("A~B", "reciprocity")]) and q[("A~B", "link_symmetry")] == 0  # B's call is U


def test_qc_annot_sensitivity_and_gap_rate(tmp_path):
    a1, a2, a3, a4 = (el("A", s, s + 999) for s in (1001, 5001, 9001, 13001))
    b1 = el("B", 21001, 22000)
    _, qc = build(tmp_path, {"A": [a1, a2, a3, a4], "B": [b1]}, [
        call(a1, "B", "P", ("chr1", 21001, 22000)),  # high; B's annotation has it (b1)
        call(b1, "A", "P", ("chr1", 1001, 2000)),
        call(a2, "B", "P", ("chr1", 31001, 32000)),  # high; B's annotation missed it
        call(a3, "B", "P", ("chr1", 41001, 42000), tier="medium"),  # not high: not counted
        call(a4, "B", "E", ("chr1", 51001, 51005)),  # not P
        call(a1, "C", "G"),
    ])
    q = qc_map(qc)
    assert q[("B", "annot_sensitivity")] == 0.5  # a1's locus has b1, a2's has no B member
    assert q[("A", "annot_sensitivity")] == 1.0  # b1's P in A lands on a1
    assert ("C", "annot_sensitivity") not in q  # unannotated
    assert (q[("A", "gap_rate")], q[("B", "gap_rate")], q[("C", "gap_rate")]) == (0, 0, 0.2)  # C: 1 G of 5


def test_qc_annot_sensitivity_without_projection(tmp_path):
    a1, b1 = el("A", 1001, 2000), el("B", 1001, 2000)
    q = qc_map(build(tmp_path, {"A": [a1], "B": [b1]}, [call(a1, "B", "E", ("chr1", 5001, 5005))],
                     labels="AB")[1])
    assert math.isnan(q[("A", "annot_sensitivity")]) and math.isnan(q[("B", "annot_sensitivity")])


def test_qc_skipped_counts(tmp_path):
    _, qc = build(tmp_path, {"A": [el("A", 1001, 2000)]}, labels="AC",
                  skipped={"A": Counter(no_ltr_pair=3, duplicate=1), "C": Counter()})
    assert [r for r in qc if r[1].startswith("skipped_")] == [("A", "skipped_duplicate", 1),
                                                             ("A", "skipped_no_ltr_pair", 3)]


def test_qc_row_order(tmp_path):
    elements, calls = _world()
    _, qc = build(tmp_path, elements, calls, skipped={"B": Counter(no_ltrs=2)})
    assert list(dict.fromkeys(s for s, _, _ in qc)) == ["all", "A>B", "A>C", "B>A", "B>C", "A~B", "A", "B", "C"]
    pair = [m for s, m, _ in qc if s == "A>B"]
    assert pair == [*(f"frac_{x}" for x in STATES), "resolved_frac", "high_frac", "annot_agree", "e_in_annot"]
    assert [m for s, m, _ in qc if s == "A>C"] == pair[:-2]  # C unannotated
    assert [m for s, m, _ in qc if s == "B"] == ["annot_sensitivity", "gap_rate", "skipped_no_ltrs"]
    assert [m for s, m, _ in qc if s == "C"] == ["gap_rate"]


# ---------------------------------------------------------------- inferred absence (A)


def test_inferred_absence_conflicts_with_observed_absence_unless_accepted(tmp_path):
    """A call A is kept as state A; against E in the top tier it conflicts, as present against absent does,
    unless the caller accepts A as an absence."""
    a1, b1 = el("A", 1001, 2000), el("B", 5001, 6000)
    calls = [call(a1, "B", "P", ("chr1", 5001, 6000)), call(a1, "C", "A"), call(b1, "C", "E", ("chr1", 3001, 3005))]
    loci, qc = build(tmp_path / "x", {"A": [a1], "B": [b1]}, calls)
    assert loci[0].states["C"].state == "U" and loci[0].states["C"].flags == ["conflict"]
    q = {(s, m): v for s, m, v in qc}
    assert q[("A>C", "frac_A")] == 1 and ("A>B", "frac_A") not in q  # a column only where an A occurs
    loci, _ = build(tmp_path / "y", {"A": [a1], "B": [b1]}, calls, absent=OBSERVED_ABSENT | INFERRED_ABSENT)
    st = loci[0].states["C"]
    assert (st.state, st.flags) == ("E", [])  # one class; E precedes A in a tie of the absent class


def test_inferred_absence_alone_is_called_a(tmp_path):
    a1 = el("A", 1001, 2000)
    loci, _ = build(tmp_path, {"A": [a1]}, [call(a1, "C", "A")], labels="AC")
    assert loci[0].states["C"].state == "A"
