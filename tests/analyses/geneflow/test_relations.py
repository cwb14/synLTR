"""Relationships (analyses.geneflow.relations): neighbor joining and rooting, PAV distances, shared derived
insertions, the sharing counts that pick the sister pair, the KKSC symmetry test read from the third genome's
own calls, the one discovery function and its four-point fallback, outgroup conflicts and the run summary,
on toy run directories (catalog_toy.py)."""
from __future__ import annotations

import math
from types import SimpleNamespace

import numpy as np
import pytest
from analyses.geneflow.catalog_toy import GENOMES, dmat, loc, write_run
from scipy.stats import binomtest

from ltrtrace.analyses.geneflow.relations import (
    Node,
    discover,
    eligible_fallback,
    four_point_ok,
    neighbor_joining,
    pav_distances,
    relations,
    root_on,
    run_summary,
    shared_derived,
    sister_pair,
)

D5 = np.array([[0, 5, 9, 9, 8], [5, 0, 10, 10, 9], [9, 10, 0, 8, 7], [9, 10, 8, 0, 3], [8, 9, 7, 3, 0]], float)


def leaf_distances(t: Node, labels: list[str]) -> np.ndarray:
    """Path lengths between leaves on the unrooted tree."""
    edges, stack = {}, [t]
    while stack:
        node = stack.pop()
        for c in node.children:
            edges.setdefault(id(node), []).append((id(c), c.length)); edges.setdefault(id(c), []).append((id(node), c.length))
            stack.append(c)
    leaf = {n.name: id(n) for n in _walk(t) if not n.children}
    def dist_from(a):
        seen, todo = {a: 0.0}, [a]
        while todo:
            u = todo.pop()
            for v, w in edges.get(u, ()):
                if v not in seen:
                    seen[v] = seen[u] + w; todo.append(v)
        return seen
    return np.array([[dist_from(leaf[a])[leaf[b]] for b in labels] for a in labels])


def splits(t: Node) -> set[frozenset[str]]:
    """Leaf sets below every inner node."""
    return {n.leaves() for n in _walk(t) if n.children}


def _walk(t):
    yield t
    for c in t.children:
        yield from _walk(c)


def test_nj_recovers_additive_tree():
    t = neighbor_joining(D5, list("abcde"))
    assert leaf_distances(t, list("abcde")) == pytest.approx(D5)
    assert {frozenset("ab"), frozenset("de")} <= splits(t) | {frozenset("abcde") - s for s in splits(t)}


def test_root_on_outgroup():
    t = root_on(neighbor_joining(D5, list("abcde")), "e")
    assert t.leaves() == frozenset("abcde")
    assert {c.leaves() for c in t.children} == {frozenset("e"), frozenset("abcd")} and t.newick().endswith(";")


R_S_SISTER = ([loc(1000 * i, "R:P* S:P D:E O:E") for i in range(1, 13)]          # R+S x12
              + [loc(20_000 + 1000 * i, "R:P S:E D:P* O:E") for i in range(2)]      # R+D x2, D members
              + [loc(30_000, "R:E S:P D:P* O:E")]                                   # S+D x1, D member
              + [loc(40_000 + 1000 * i, "R:P* S:P D:E O:P") for i in range(3)]      # outgroup present x3
              + [loc(50_000, "R:P* S:U D:E O:E")])                                  # S unresolved


def test_sharing_and_sister(tmp_path, toy_fastas):
    cat = write_run(tmp_path, R_S_SISTER, fasta_dir=toy_fastas).load()
    sc = sister_pair(cat, ("R", "S", "D"), pav_distances(cat)[0])
    assert sc.pair == ("R", "S") and sc.donor == "D" and sc.counts == {"R+S": 12, "R+D": 2, "S+D": 1}
    assert sc.p_top == pytest.approx(binomtest(12, 14).pvalue) and not sc.ambiguous          # 0.0129
    assert sc.p_kksc == pytest.approx(1.0) and sc.kksc_loci == "disc=D;present=D;res=R,S;outgroup=not_present"


def test_kksc_reads_the_third_genomes_calls(tmp_path, toy_fastas):
    # sister pair (R, S), donor D: D's own calls are read. One more R+D locus whose merged S is E (R, a member, calls it E)
    # while D's own call into S is U: the merged sharing counts it, p_kksc does not
    loci = R_S_SISTER + [loc(60_000, "R:P* S:E D:P* O:E", calls={("D", "S"): "U"})]
    cat = write_run(tmp_path / "a", loci, fasta_dir=toy_fastas).load()
    sc = sister_pair(cat, ("R", "S", "D"), pav_distances(cat)[0])
    assert sc.pair == ("R", "S") and sc.counts["R+D"] == 3 and not sc.ambiguous             # binomtest(12, 15): 0.035
    assert sc.p_kksc == pytest.approx(binomtest(2, 3).pvalue)                               # n_R+D 2, n_S+D 1
    unannotated_d = [GENOMES[0], GENOMES[1], ("D", "taxD", None, False), GENOMES[3]]
    loci2 = [loc(1000 * i, "R:P* S:P D:E O:E") for i in range(1, 13)] + [loc(30_000, "R:P* S:E D:P O:E")]
    cat2 = write_run(tmp_path / "b", loci2, genomes=unannotated_d, fasta_dir=toy_fastas).load()
    sc2 = sister_pair(cat2, ("R", "S", "D"), pav_distances(cat2)[0])
    assert math.isnan(sc2.p_kksc) and sc2.kksc_loci == "none"


def test_discover_rules(tmp_path, toy_fastas):
    loci = [loc(10_000, "R:P S:E D:P* O:E"), loc(20_000, "R:P S:E D:P* O:U"), loc(30_000, "R:P S:E D:P* O:P"),
            loc(40_000, "R:P* S:E D:P O:E"), loc(50_000, "R:P S:U D:P* O:E")]
    cat = write_run(tmp_path, loci, fasta_dir=toy_fastas).load()
    dist = pav_distances(cat)[0]
    e, np_ = discover(cat, "R", "S", "D", dist, "E"), discover(cat, "R", "S", "D", dist, "not_present")
    assert e.cand.tolist() == [True, False, False, False, True] and e.rows.tolist() == [True, False, False, False, False]
    assert np_.cand.tolist() == [True, True, False, False, True]
    assert e.gstar.tolist() == [2, 2, 2, -1, 2] and e.disc == ("D",)
    assert e.loci == "disc=D;present=D;res=R,S;outgroup=E" and np_.loci.endswith("outgroup=not_present")


def test_sister_tie_warns(tmp_path, toy_fastas):
    loci = ([loc(1000 * i, "R:P* S:P D:E O:E") for i in range(1, 4)]
            + [loc(10_000 + 1000 * i, "R:P* S:E D:P O:E") for i in range(3)])
    cat = write_run(tmp_path, loci, fasta_dir=toy_fastas).load()
    _, f, _, dist = relations(cat)
    sc = sister_pair(cat, ("R", "S", "D"), dist)
    assert sc.pair == ("R", "S") and sc.ambiguous
    assert any(t.startswith("WARNING:") and "R+S" in t for _, t in f.items())


def test_outgroup_conflict(tmp_path, toy_fastas):
    # counted from the discoverers' own calls, not from the merged states (which give 1/8 of 8): R's seven
    # member loci, P by construction, do not enter; S discovers the E locus and the two U loci, where S's own call
    # into R is U
    loci = ([loc(1000 * i, "R:P* S:P D:P O:P") for i in range(1, 8)] + [loc(9000, "R:E S:P* D:P O:P")]
            + [loc(10_000 + 1000 * i, "R:U S:P* D:P O:P") for i in range(2)])
    cat = write_run(tmp_path, loci, fasta_dir=toy_fastas).load()
    rows = relations(cat)[0]
    r = next(x for x in rows if x.test == "outgroup_conflict" and x.scope == "R")
    assert r.value == pytest.approx(1.0) and r.n == 1 and "upper bound" in r.note
    assert r.loci == "disc=S,D;res=R;outgroup=present"
    s = next(x for x in rows if x.test == "outgroup_conflict" and x.scope == "S")
    assert (s.value, s.n) == (0.0, 7)                    # R discovers S's seven loci: S P in R's own calls


def test_distances_pairwise_complete(tmp_path, toy_fastas):
    loci = [loc(1000, "R:P* S:P D:E O:E"), loc(2000, "R:P* S:E D:E O:E"), loc(3000, "R:E S:P* D:E O:E"),
            loc(4000, "R:P* S:U D:E O:E")]
    cat = write_run(tmp_path, loci, fasta_dir=toy_fastas).load()
    d, n = pav_distances(cat)
    assert d[cat.g("R"), cat.g("S")] == pytest.approx(2 / 3) and n[cat.g("R"), cat.g("S")] == 3


def test_shared_derived(tmp_path, toy_fastas):
    loci = ([loc(1000 * i, "R:P* S:P D:E O:E") for i in range(1, 4)]
            + [loc(5000 + 1000 * i, "R:P* S:P D:P O:E") for i in range(2)] + [loc(9000, "R:P* S:P D:E O:P")])
    cat = write_run(tmp_path, loci, fasta_dir=toy_fastas).load()
    assert shared_derived(cat)[cat.g("R"), cat.g("S")] == 3


def test_four_point():
    cat = SimpleNamespace(g=["R", "S", "D", "X"].index)
    near_d = np.array([[0, .2, .5, .5], [.2, 0, .5, .5], [.5, .5, 0, .1], [.5, .5, .1, 0]])
    near_r = np.array([[0, .2, .5, .1], [.2, 0, .5, .25], [.5, .5, 0, .5], [.1, .25, .5, 0]])
    assert four_point_ok(near_d, cat, "R", "S", "X", "D")          # 0.3 < 1.0
    assert not four_point_ok(near_r, cat, "R", "S", "X", "D")      # 0.7 > 0.6


def test_run_summary_quotes_qc(tmp_path, toy_fastas):
    cat = write_run(tmp_path, [loc(1000, "R:P* S:E D:P* O:E")], fasta_dir=toy_fastas).load()
    rows, f = run_summary(cat)
    lines = [t for s, t in f.items() if s == "Run"]
    assert any("annot_sensitivity" in t and "0.9" in t for t in lines)
    assert any("skipped_malformed_coordinates 2" in t for t in lines)
    assert next(r for r in rows if r.statistic == "resolved_frac" and r.scope == "S").value == 1.0


# ---------------------------------------------------------------- edge cases


def test_newick_and_nj_details():
    t = Node(None, (Node("a", length=0.1234567), Node(None, (Node("b", length=2.0), Node("c", length=1e-7)),
                                                         length=0.5)))
    assert t.newick() == "(a:0.123457,(b:2,c:1e-07):0.5);" and t.leaves() == frozenset("abc")
    two = neighbor_joining(np.array([[0, 0.4], [0.4, 0]]), ["a", "b"])
    assert [(c.name, c.length) for c in two.children] == [("a", pytest.approx(0.2)),
                                                          ("b", pytest.approx(0.2))]
    star = neighbor_joining(np.ones((4, 4)) - np.eye(4), list("abcd"))     # every Q ties: a, b join first
    assert frozenset("ab") in splits(star) and len(star.children) == 3
    neg = neighbor_joining(np.array([[0, 1, 1], [1, 0, 3], [1, 3, 0]], float), list("abc"))
    assert [c.length for c in neg.children] == [0.0, 1.5, 1.5]             # a's -0.5 set to 0
    with pytest.raises(ValueError, match="finite"):
        neighbor_joining(np.array([[0, np.nan], [np.nan, 0]]), ["a", "b"])


def test_root_on_keeps_paths_and_merges_old_root():
    rooted = Node(None, (Node(None, (Node("a", length=1), Node("b", length=1)), length=1),
                         Node(None, (Node("c", length=1), Node("d", length=1)), length=1)))
    t = root_on(rooted, "c")
    assert t.newick() == "(c:0.5,((a:1,b:1):2,d:1):0.5);"
    assert leaf_distances(t, list("abcd")) == pytest.approx(leaf_distances(rooted, list("abcd")))
    two = neighbor_joining(np.array([[0, 0.4], [0.4, 0]]), ["a", "b"])
    assert root_on(two, "b").newick() == "(b:0.2,a:0.2);"                # the old root merged away
    with pytest.raises(ValueError, match="x"):
        root_on(rooted, "x")


def test_fallback_discovery_reasons(tmp_path, toy_fastas):
    genomes = [("R", "taxR", None, True), ("S", "taxS", None, True), ("D", "taxD", None, False),
               ("X", "taxD", None, True), ("Y", "taxR", None, True), ("O", "taxO", "outgroup", False)]
    cat = write_run(tmp_path, [loc(1000, "R:P S:E D:P X:P* Y:E O:E"), loc(2000, "R:P* S:E D:P X:P Y:E O:E")],
                    genomes=genomes, fasta_dir=toy_fastas).load()
    dist = dmat({("R", "S"): 0.2, ("D", "X"): 0.1, ("R", "Y"): 0.1, ("S", "Y"): 0.25}, 0.5, cat.labels)
    assert eligible_fallback(cat, "R", "S", "D", dist) == (("X",), {"Y": "groups with R"})
    c = discover(cat, "R", "S", "D", dist, "E")
    assert c.disc == ("X",) and c.excluded == {"Y": "groups with R"} and c.gstar.tolist() == [cat.g("X"), -1]
    assert c.loci == "disc=X;present=D;res=R,S;outgroup=E" and c.rows.tolist() == [True, False]
    dist[cat.g("D"), cat.g("X")] = dist[cat.g("X"), cat.g("D")] = np.nan
    excluded = eligible_fallback(cat, "R", "S", "D", dist)[1]
    assert excluded == {"X": "distance unresolved", "Y": "groups with R"}
    assert discover(cat, "R", "S", "D", dist, "E") is None
    with pytest.raises(ValueError, match="outgroup"):
        discover(cat, "R", "S", "D", dist, "absent")


def test_relations_skips_and_rooting(tmp_path, toy_fastas):
    # two ingroup genomes and an outgroup: a rooted three-leaf tree, no trio
    two = write_run(tmp_path / "a", [loc(1000 * i, "R:P* S:E O:E") for i in range(1, 4)]
                    + [loc(5000, "R:P* S:P O:P")], genomes=[GENOMES[0], GENOMES[1], GENOMES[3]],
                    fasta_dir=toy_fastas).load()
    rows, f, tree, _ = relations(two)
    assert tree.children[0].name == "O" and next(r for r in rows if r.test == "tree").value.startswith("(O:")
    assert not [r for r in rows if r.test == "sharing"]
    assert any("fewer than three ingroup genomes" in t for _, t in f.items())
    # no outgroup: an unrooted tree; D unannotated with no fallback: the KKSC test is skipped with a WARNING
    no_out = [GENOMES[0], GENOMES[1], ("D", "taxD", None, False)]
    cat = write_run(tmp_path / "b", [loc(1000 * i, "R:P* S:P D:E") for i in range(1, 4)]
                    + [loc(9000, "R:P* S:E D:P")], genomes=no_out, fasta_dir=toy_fastas).load()
    rows, f, tree, _ = relations(cat)
    assert len(tree.children) == 3 and not [r for r in rows if r.test == "outgroup_conflict"]
    kksc = next(r for r in rows if r.statistic == "p_kksc")
    assert math.isnan(kksc.value) and kksc.loci == "none"
    assert any(t.startswith("WARNING:") and "skipped" in t and "KKSC" in t for _, t in f.items())
    # a pair never resolved together: no distance, no tree
    gap = write_run(tmp_path / "c", [loc(1000, "R:P* S:U D:P* O:E"), loc(2000, "R:U S:P* D:P* O:E")],
                    fasta_dir=toy_fastas).load()
    rows, f, tree, dist = relations(gap)
    assert tree is None and math.isnan(dist[gap.g("R"), gap.g("S")])
    assert any(t.startswith("Tree skipped") and "R" in t and "S" in t for _, t in f.items())


def test_sister_without_informative_loci_is_ambiguous(tmp_path, toy_fastas):
    cat = write_run(tmp_path, [loc(1000, "R:P* S:P D:P* O:E")], fasta_dir=toy_fastas).load()
    sc = sister_pair(cat, ("D", "R", "S"), pav_distances(cat)[0])           # any order: sheet order inside
    assert sc.pair == ("R", "S") and sc.donor == "D" and math.isnan(sc.p_top) and sc.ambiguous
    with pytest.raises(ValueError, match="O"):
        sister_pair(cat, ("R", "S", "O"), pav_distances(cat)[0])


# ---------------------------------------------------------------- discoverers' own calls, fallbacks, two outgroups


def test_outgroup_conflict_ignores_annotation_depth(tmp_path, toy_fastas):
    # the same eight S-discovered loci with O present (S's own call puts R E at one), annotated in R shallowly
    # or deeply: deep also makes R a member where it is present and adds eight loci only R discovers
    shallow = [loc(1000 * i, "R:P S:P* D:P O:P") for i in range(1, 8)] + [loc(9000, "R:E S:P* D:P O:P")]
    deep = ([loc(1000 * i, "R:P* S:P* D:P O:P") for i in range(1, 8)] + [loc(9000, "R:E S:P* D:P O:P")]
            + [loc(20_000 + 1000 * i, "R:P* S:P D:P O:P") for i in range(8)])
    got, merged = [], []
    for name, loci in (("shallow", shallow), ("deep", deep)):
        cat = write_run(tmp_path / name, loci, fasta_dir=toy_fastas).load()
        r = next(x for x in relations(cat)[0] if x.test == "outgroup_conflict" and x.scope == "R")
        got.append((r.value, r.n))
        m = (cat.cls("O") == 1) & (cat.cls("R") >= 0)                 # the merged-state rule, for contrast
        merged.append(float((m & (cat.cls("R") == 0)).sum() / m.sum()))
    assert got == [(pytest.approx(1 / 8), 8), (pytest.approx(1 / 8), 8)]
    assert merged == [pytest.approx(1 / 8), pytest.approx(1 / 16)]         # diluted by R's own member loci


def test_rows_declare_their_loci(tmp_path, toy_fastas):
    cat = write_run(tmp_path, R_S_SISTER, fasta_dir=toy_fastas).load()
    rows = run_summary(cat)[0] + relations(cat)[0]
    assert all(r.loci != "." for r in rows)
    loci = {(r.scope, r.statistic): r.loci for r in rows if r.test == "catalog"}
    assert loci[("R", "annotated")] == loci[("R", "skipped_malformed_coordinates")] == "none"
    assert loci[("R", "annot_sensitivity")] == "disc=any_but_R;present=R"


def test_fallback_first_member_and_exclusions(tmp_path, toy_fastas):
    genomes = [("R", "taxR", None, True), ("S", "taxS", None, True), ("D", "taxD", None, False),
               ("X1", "taxD", None, True), ("X2", "taxD", None, True), ("Y", "taxS", None, True),
               ("O", "taxO", "outgroup", False)]
    loci = [loc(1000, "R:P S:E D:P X1:P* X2:P* Y:E O:E"),                                 # both: gstar X1
            loc(2000, "R:P S:E D:P X1:E X2:P* Y:E O:E"),                                  # X2 only: gstar X2
            loc(3000, "R:P S:E D:P X1:P* X2:P* Y:E O:E", calls={("X1", "D"): "U"})]       # X1's U drops it
    cat = write_run(tmp_path, loci, genomes=genomes, fasta_dir=toy_fastas).load()
    dist = dmat({("R", "S"): 0.2, ("D", "X1"): 0.1, ("D", "X2"): 0.1, ("X1", "X2"): 0.1, ("S", "Y"): 0.1,
                 ("R", "Y"): 0.25}, 0.5, cat.labels)
    c = discover(cat, "R", "S", "D", dist, "E")
    assert c.disc == ("X1", "X2") and c.excluded == {"Y": "groups with S"}
    assert c.gstar.tolist() == [cat.g("X1"), cat.g("X2"), cat.g("X1")]
    assert c.cand.tolist() == [True, True, False]


def test_kksc_line_names_the_fallback_exclusions(tmp_path, toy_fastas):
    # PAV distances from the loci: X splits with the unannotated D from R and S, Y groups with S
    genomes = [("R", "taxR", None, True), ("S", "taxS", None, True), ("D", "taxD", None, False),
               ("X", "taxD", None, True), ("Y", "taxS", None, True), ("O", "taxO", "outgroup", False)]
    loci = ([loc(1000 * i, "R:P* S:P D:E X:E Y:P O:E") for i in range(1, 4)]
            + [loc(5000 + 1000 * i, "R:E S:P* D:E X:E Y:P O:E") for i in range(3)]
            + [loc(9000 + 1000 * i, "R:E S:E D:P X:P* Y:E O:E") for i in range(3)]
            + [loc(13_000, "R:P S:E D:P X:P* Y:E O:E")]
            + [loc(14_000 + 1000 * i, "R:P* S:P D:P X:P Y:P O:E") for i in range(2)])
    cat = write_run(tmp_path, loci, genomes=genomes, fasta_dir=toy_fastas).load()
    _, f, _, dist = relations(cat)
    assert eligible_fallback(cat, "R", "S", "D", dist) == (("X",), {"Y": "groups with S"})
    line = next(t for _, t in f.items() if t.startswith("KKSC symmetry of R,S;D"))
    assert "on loci X discovers" in line and "not eligible: Y: groups with S" in line


def test_kksc_reads_a_member_as_the_discoverer_calls_it(tmp_path, toy_fastas):
    # R is a member (merged P/annot), D's own call into R is E (annot_conflict): the merged sharing counts an
    # R+D locus; the KKSC test reads D's E, so neither of its patterns
    loci = R_S_SISTER + [loc(60_000, "R:P* S:E D:P* O:E", calls={("D", "R"): "E"})]
    cat = write_run(tmp_path, loci, fasta_dir=toy_fastas).load()
    sc = sister_pair(cat, ("R", "S", "D"), pav_distances(cat)[0])
    assert cat.cls("R")[-1] == 1 and cat.pcls("D", "R")[-1] == 0
    assert sc.counts["R+D"] == 3 and sc.p_kksc == pytest.approx(binomtest(2, 3).pvalue)


def test_two_outgroups(tmp_path, toy_fastas):
    genomes = [*GENOMES[:3], ("O1", "taxO", "outgroup", False), ("O2", "taxO", "outgroup", False)]
    loci = ([loc(1000 * i, "R:P* S:P D:E O1:E O2:E") for i in range(1, 5)]
            + [loc(6000 + 1000 * i, "R:P S:E D:P* O1:E O2:E") for i in range(2)]
            + [loc(9000 + 1000 * i, "R:P* S:P D:P O1:E O2:E") for i in range(6)]
            + [loc(16_000 + 1000 * i, "R:P* S:P D:P O1:P O2:P") for i in range(3)]
            + [loc(20_000, "R:E S:P* D:P O1:P O2:E"),                 # only O1 present: counts, R E
               loc(21_000, "R:P S:P* D:P O1:U O2:P"),                 # only O2 present: counts, R P
               loc(22_000, "R:E S:P* D:E O1:E O2:E")])                # no outgroup present: left out
    run = write_run(tmp_path / "clade", loci, genomes=genomes, fasta_dir=toy_fastas)
    rows, _, tree, _ = relations(run.load())
    assert [c.leaves() for c in tree.children] == [frozenset({"O1", "O2"}), frozenset("RSD")]
    r = next(x for x in rows if x.test == "outgroup_conflict" and x.scope == "R")
    assert (r.value, r.n) == (pytest.approx(0.5), 2)
    # O2 carries what R carries: the outgroups form no clade; the tree is rooted on the first, with a WARNING
    loci = ([loc(1000 * i, "R:P* S:P D:E O1:E O2:P") for i in range(1, 5)]
            + [loc(6000 + 1000 * i, "R:P S:E D:P* O1:E O2:P") for i in range(2)]
            + [loc(9000 + 1000 * i, "R:P* S:P D:P O1:P O2:P") for i in range(3)]
            + [loc(13_000 + 1000 * i, "R:E S:P* D:P O1:E O2:E") for i in range(2)])
    run = write_run(tmp_path / "split", loci, genomes=genomes, fasta_dir=toy_fastas)
    _, f, tree, _ = relations(run.load())
    assert tree.children[0].name == "O1"
    assert any(t.startswith("WARNING:") and "do not form a clade" in t for _, t in f.items())


def test_run_summary_names_a_foreign_qc_tsv(tmp_path, toy_fastas):
    """A qc.tsv with another header quotes nothing and says so in one Run line, naming the command that writes a
    qc.tsv (ltrtrace genotype)."""
    run = write_run(tmp_path, [loc(1000, "R:P* S:E D:P* O:E")], fasta_dir=toy_fastas)
    (run.outdir / "qc.tsv").write_text("a\tb\tc\n")
    rows, f = run_summary(run.load())
    lines = [t for s, t in f.items() if s == "Run"]
    assert lines[-1] == (f"{run.outdir / 'qc.tsv'}: not a qc.tsv of ltrtrace genotype (header not scope, metric, value); "
                         "nothing quoted")
    assert not [r for r in rows if r.statistic == "annot_sensitivity"]
