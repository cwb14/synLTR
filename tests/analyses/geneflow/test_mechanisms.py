"""Mechanisms and missingness (analyses.geneflow.mechanisms): states by family and age with O and D split by
their junctions, absence fractions and events, missingness and Dollo placement, on toy run directories
(catalog_toy.py)."""
from __future__ import annotations

import gzip

import pytest
from analyses.geneflow.catalog_toy import GENOMES, loc, write_run

from synltr.analyses.geneflow.catalog import AnalysisError
from synltr.analyses.geneflow.mechanisms import (
    MECH_COLUMNS,
    STATE_METRICS,
    _absence,
    _dollo,
    d_kind,
    mechanisms,
    o_kind,
    pair_calls,
    write_mechanisms,
)
from synltr.analyses.geneflow.relations import Node
from synltr.calls.records import CALL_COLUMNS

TREE = Node(None, (Node("O"), Node(None, (Node(None, (Node("R"), Node("S"))), Node("D")))))  # ((R,S),D),O


def _table(rows, table, scope):
    """{(superfamily, family, age_bin, metric): value} of one table and scope."""
    return {(sf, fam, age, m): v for t, sc, sf, fam, age, m, v in rows if (t, sc) == (table, scope)}


def _lines(f):
    return [t for _, t in f.items()]


def _rewrite(path, edit):
    """Rewrite pair file `path` through edit(rows) -> rows, each a dict CALL_COLUMNS -> field."""
    with gzip.open(path, "rt") as fh:
        lines = fh.read().splitlines()
    rows = edit([dict(zip(CALL_COLUMNS, x.split("\t"))) for x in lines[1:]])
    with gzip.open(path, "wt") as fh:
        fh.write("\n".join(["\t".join(CALL_COLUMNS), *("\t".join(r[c] for c in CALL_COLUMNS) for r in rows)]) + "\n")


def test_o_kind_and_d_kind():
    assert o_kind(3, -2, "ACGTA") == "O_site" and o_kind(400, 0, "ACGTA") == "O_complex"
    assert o_kind(10, 10, None) == "O_site" and o_kind(11, 0, None) == "O_complex"    # TSD_DEFAULT 5 -> w 10
    assert d_kind(-12, 0, "ACGTA") == "D_inside" and d_kind(300, 0, "ACGTA") == "D_beyond"


def test_state_counts_by_family_age(tmp_path, toy_fastas):
    loci = [loc(1000, "R:P* S:E D:E O:E", k=0.0), loc(2000, "R:P* S:E D:E O:E", k=0.02),
            loc(3000, "R:P* S:D D:E O:E", k=0.02, family="famB", superfamily="Copia", evidence={"S": (300, 0, 0)}),
            loc(4000, "R:P* S:O D:E O:E", k=0.004, family="famB", superfamily="Copia", evidence={"S": (2, -1, 500)}),
            loc(5000, "R:P* S:P D:E O:E", k=0.0)]
    cat = write_run(tmp_path, loci, fasta_dir=toy_fastas).load()
    st = _table(mechanisms(cat, None)[0], "states", "R>S")
    assert st[("Gypsy", "famA", "K=0", "n_E")] == 1 and st[("Gypsy", "famA", "K>=0.01", "n_E")] == 1
    assert st[("Copia", "famB", "K>=0.01", "n_D_beyond")] == 1 and st[("Copia", "famB", "0<K<0.01", "n_O_site")] == 1
    assert st[("all", "all", "all", "n_E")] == 2 and st[("all", "all", "all", "n_P")] == 1
    assert st[("Gypsy", "all", "K=0", "n_P")] == 1


def test_absence_fractions(tmp_path, toy_fastas):
    specs = ([("E", (0, 0, 0))] * 4 + [("D", (-12, 0, 0))] + [("D", (300, 0, 0))] * 2
             + [("O", (2, -1, 500))] + [("O", (400, 0, 900))] * 2)
    loci = [loc(1000 * (i + 1), f"R:P* S:{s} D:E O:E", evidence={"S": ev}) for i, (s, ev) in enumerate(specs)]
    cat = write_run(tmp_path, loci, fasta_dir=toy_fastas).load()
    ab = _table(mechanisms(cat, None)[0], "absence", "R>S")
    got = [ab[("all", "all", "all", m)] for m in ("n_absent", "frac_E", "frac_D_inside", "frac_D_beyond",
                                                  "frac_O_site", "frac_O_complex")]
    assert got == [10, 0.4, 0.1, 0.2, 0.1, 0.2]


def test_absence_events(tmp_path, toy_fastas):
    loci = [loc(10_000, "R:P* S:E D:E O:E", length=6000),
            loc(11_000, "R:P* S:E D:E O:E", parent=0, flags={"S": ["host_empty"]}),
            loc(12_000, "R:P* S:E D:E O:E", parent=0, flags={"S": ["host_empty"]}),
            loc(30_000, "R:P* S:D D:E O:E", length=6000),
            loc(31_000, "R:P* S:D D:E O:E", parent=3, flags={"S": ["inherit"]}),
            loc(50_000, "R:P* S:E D:E O:E"),
            loc(60_000, "R:P* S:E D:E O:E", flags={"S": ["adjacent"]}, wid=7),
            loc(60_400, "R:P* S:E D:E O:E", flags={"S": ["adjacent"]}, wid=7)]
    cat = write_run(tmp_path, loci, fasta_dir=toy_fastas).load()
    ev = {m: v for (sf, fam, age, m), v in _table(mechanisms(cat, None)[0], "events", "R>S").items()}
    assert (ev["n_E_single"], ev["n_E_cluster"], ev["n_E_cluster_elements"]) == (1, 2, 5)
    assert (ev["n_D_single"], ev["n_D_multi"], ev["n_D_multi_elements"]) == (0, 1, 2)


def test_dollo_branches(tmp_path, toy_fastas):
    loci = [loc(1000, "R:P* S:P D:E O:E"),      # insertion on R+S
            loc(2000, "R:S S:P* D:E O:E"),      # insertion on R+S, solo on R
            loc(3000, "R:D S:P* D:P O:E"),      # insertion on R+S+D, deletion on R
            loc(4000, "R:D S:P* D:E O:E"),      # insertion on S; R outside the clade: missing, no deletion
            loc(5000, "R:P* S:E D:P O:E")]      # incongruent: S E inside the R+S+D clade
    cat = write_run(tmp_path, loci, fasta_dir=toy_fastas).load()
    tree = Node(None, (Node("O"), Node(None, (Node(None, (Node("R"), Node("S"))), Node("D")))))
    rows, *_ = _dollo(cat, tree)
    br = lambda scope, m: _table(rows, "branches", scope).get(("all", "all", "all", m), 0)
    assert (br("R+S", "n_insertion"), br("R+S+D", "n_insertion"), br("S", "n_insertion")) == (2, 1, 1)
    assert (br("R", "n_solo"), br("R", "n_deletion"), br("incongruent", "n_loci")) == (1, 1, 1)


def test_missingness(tmp_path, toy_fastas):
    loci = ([loc(1000 * (i + 1), "R:P* S:U D:E O:E") for i in range(3)]
            + [loc(5000, "R:P* S:G D:E O:E"), loc(6000, "R:P* S:E D:E O:E")])
    cat = write_run(tmp_path, loci, fasta_dir=toy_fastas).load()
    mi = _table(mechanisms(cat, None)[0], "missingness", "R>S")
    assert mi[("all", "all", "all", "U_no_anchor")] == pytest.approx(3 / 5) and mi[("all", "all", "all", "G")] == pytest.approx(1 / 5)


def test_write_mechanisms(tmp_path):
    rows = [("states", "R>S", "all", "all", "all", "n_E", 2), ("absence", "R>S", "all", "all", "K=0", "frac_E", 0.5)]
    write_mechanisms(rows, tmp_path / "m.tsv")
    lines = (tmp_path / "m.tsv").read_text().splitlines()
    assert lines[0].split("\t") == list(MECH_COLUMNS)
    assert lines[1:] == ["absence\tR>S\tall\tall\tK=0\tfrac_E\t0.5", "states\tR>S\tall\tall\tall\tn_E\t2"]


# ---------------------------------------------------------------- geometry, events, damaged inputs, tables, Dollo losses


def test_inherited_calls_are_undetermined(tmp_path, toy_fastas):
    """A D or O passed from its host (rule 0, flag inherit) has no junction of its own (dl, dr "."): never
    shown derived, never O at the site."""
    loci = [loc(10_000, "R:P* S:D D:E O:E", length=6000, evidence={"S": (-12, 0, 0)}),          # D_inside
            loc(11_000, "R:P* S:D D:E O:E", parent=0, flags={"S": ["inherit"]},
                evidence={"S": (None, None, None)}),                                             # inherited D
            loc(20_000, "R:P* S:O D:E O:E", length=6000, evidence={"S": (2, -1, 500)}),         # O_site
            loc(21_000, "R:P* S:O D:E O:E", parent=2, flags={"S": ["inherit"]},
                evidence={"S": (None, None, 500)})]                                              # inherited O
    cat = write_run(tmp_path, loci, fasta_dir=toy_fastas).load()
    st = _table(mechanisms(cat, None)[0], "states", "R>S")
    assert [st[("all", "all", "all", m)] for m in ("n_D_inside", "n_D_beyond", "n_O_site", "n_O_inside",
                                                   "n_O_complex")] == [1, 1, 1, 0, 1]


def test_o_inside_is_derived(tmp_path, toy_fastas):
    """An O call whose gap end lies more than w inside the element: B kept part of it, the geometry of
    D_inside, so it is derived and shown, and the lower bound adds it."""
    assert o_kind(-12, 0, "ACGTA") == "O_inside" and o_kind(400, -11, None) == "O_inside"   # inside wins
    assert o_kind(-10, 10, None) == "O_site" and o_kind(-10, 11, None) == "O_complex"       # -w is at the site
    specs = [("E", (0, 0, 0))] * 2 + [("D", (-12, 0, 0)), ("O", (-12, 300, 500)), ("O", (400, 0, 900))]
    loci = [loc(1000 * (i + 1), f"R:P* S:{s} D:E O:E", evidence={"S": ev}) for i, (s, ev) in enumerate(specs)]
    cat = write_run(tmp_path, loci, fasta_dir=toy_fastas).load()
    rows, stat, f = mechanisms(cat, None)
    st = _table(rows, "states", "R>S")
    assert [st[("all", "all", "all", m)] for m in ("n_O_site", "n_O_inside", "n_O_complex")] == [0, 1, 1]
    ab = _table(rows, "absence", "R>S")
    assert [ab[("all", "all", "all", m)] for m in ("frac_E", "frac_D_inside", "frac_O_inside", "frac_D_beyond",
                                                   "frac_O_site", "frac_O_complex")] == [0.4, 0.2, 0.2, 0, 0, 0.2]
    note = next(r.note for r in stat if (r.test, r.scope, r.statistic) == ("mechanisms", "R>S", "frac_O_inside"))
    assert note.startswith("derived, shown")
    line = next(t for t in _lines(f) if t.startswith("R>S: "))
    assert "frac_D_inside + frac_O_inside = 0.2 + 0.2 = 0.4 (a lower bound" in line
    assert "1 - frac_E = 0.6" in line and "the bounds bracket the derived share; they do not measure it" in line


def test_inside_stretch_with_a_child_is_unshown(tmp_path, toy_fastas):
    """A gap end inside the element whose cut-off stretch holds an annotated child: the stretch is
    not the element's own sequence, so the call is inside but unshown, never shown derived."""
    inherit = {"flags": {"S": ["inherit"]}, "evidence": {"S": (None, None, None)}}
    loci = [loc(10_000, "R:P* S:D D:E O:E", length=6000, evidence={"S": (-40, 0, 0)}),     # child 20 bp in
            loc(10_020, "R:P* S:D D:E O:E", parent=0, **inherit),
            loc(20_000, "R:P* S:D D:E O:E", length=6000, evidence={"S": (-40, 0, 0)}),     # child 2 kb in: shown
            loc(22_000, "R:P* S:D D:E O:E", parent=2, **inherit),
            loc(30_000, "R:P* S:O D:E O:E", length=6000, evidence={"S": (300, -40, 500)}), # child 10 bp from the end
            loc(35_900, "R:P* S:O D:E O:E", length=90, parent=4, flags={"S": ["inherit"]},
                evidence={"S": (None, None, 500)})]
    cat = write_run(tmp_path, loci, fasta_dir=toy_fastas).load()
    rows, _, f = mechanisms(cat, None)
    st = _table(rows, "states", "R>S")
    assert [st[("all", "all", "all", m)] for m in ("n_D_inside", "n_D_inside_unshown", "n_D_beyond", "n_O_inside",
                                                   "n_O_inside_unshown", "n_O_complex")] == [1, 1, 2, 0, 1, 1]
    ab = _table(rows, "absence", "R>S")
    assert (ab[("all", "all", "all", "frac_D_inside_unshown")], ab[("all", "all", "all", "frac_O_inside_unshown")]) \
        == (1 / 6, 1 / 6)
    assert any(t.startswith("R>S: ") and "= 0.167 + 0 = 0.167 (a lower bound" in t for t in _lines(f))



def test_bounds_line_keeps_small_fractions():
    """A nonzero fraction in the bounds line never prints as 0 (3 fixed decimals would write
    24 of 76,418 absences, 0.0003, as 0.000): 3 significant digits."""
    counts = dict.fromkeys(STATE_METRICS, 0) | {"n_E": 9000, "n_D_inside": 990, "n_O_inside": 3, "n_D_beyond": 7}
    line = _absence({("all", "all", "all"): [counts[m] for m in STATE_METRICS]}, "R", "S", "R>S")[2]
    assert "frac_D_inside + frac_O_inside = 0.099 + 0.0003 = 0.0993 (a lower bound" in line
    assert "1 - frac_E = 0.1 (an upper bound" in line and "O at the site (0)" in line

def test_long_inside_stretch_is_unshown(tmp_path, toy_fastas):
    """A gap end MIN_SEG or more inside the element: rule 4 is reached with fewer than MIN_SEG own
    bases aligned, so so long a stretch cannot be aligned own sequence: unshown. One shown end is enough."""
    assert d_kind(-79, 0, "ACGTA") == "D_inside" and d_kind(-80, 0, "ACGTA") == "D_inside_unshown"
    assert o_kind(300, -80, None) == "O_inside_unshown" and o_kind(-79, -500, None) == "O_inside"
    assert d_kind(-40, 0, None, own=(39, 0)) == "D_inside_unshown" and d_kind(-40, 0, None, own=(40, 0)) == "D_inside"
    specs = [("D", (-500, 0, 0)), ("D", (-79, 0, 0)), ("O", (-2, -3000, 900))]
    loci = [loc(10_000 * (i + 1), f"R:P* S:{s} D:E O:E", length=6000, evidence={"S": ev})
            for i, (s, ev) in enumerate(specs)]
    cat = write_run(tmp_path, loci, fasta_dir=toy_fastas).load()
    st = _table(mechanisms(cat, None)[0], "states", "R>S")
    assert [st[("all", "all", "all", m)] for m in ("n_D_inside", "n_D_inside_unshown", "n_O_inside",
                                                   "n_O_inside_unshown")] == [1, 1, 0, 1]


def test_w_is_the_elements_own(tmp_path, toy_fastas):
    """w = SLOP_TSD_MULT x the element's TSD length, TSD_DEFAULT without one, as the classifier's."""
    assert d_kind(-9, 0, "ACGT") == "D_inside" and d_kind(-11, 0, "ACGTAC") == "D_beyond"        # w 8, w 12
    assert d_kind(-11, 0, None) == "D_inside" and d_kind(-11, 0, "") == "D_inside"              # w 10
    assert o_kind(12, -12, "ACGTAC") == "O_site" and o_kind(9, 0, "ACGT") == "O_complex"
    loci = [loc(1000 * (i + 1), "R:P* S:D D:E O:E", evidence={"S": (-11, 0, 0)}) for i in range(3)]
    run = write_run(tmp_path, loci, fasta_dir=toy_fastas)
    cat = run.load()
    tsd = {"1000": "ACGTAC", "2000": "."}  # a 6-bp TSD (w 12), none (w 10); the third keeps the toy's ACGTA (w 10)
    _rewrite(run.outdir / "pairs" / "R__S.calls.tsv.gz",
             lambda rows: [{**r, "tsd_a": tsd.get(r["a_start"], r["tsd_a"])} for r in rows])
    st = _table(mechanisms(cat, None)[0], "states", "R>S")
    assert (st[("all", "all", "all", "n_D_inside")], st[("all", "all", "all", "n_D_beyond")]) == (2, 1)


def test_events_group_adjacent_calls_per_window_and_state(tmp_path, toy_fastas):
    """adjacent calls join per (state, wid) only, and host links are transitive."""
    adj = {"S": ["adjacent"]}
    loci = [loc(10_000, "R:P* S:E D:E O:E", flags=adj, wid=100), loc(10_400, "R:P* S:E D:E O:E", flags=adj, wid=100),
            loc(20_000, "R:P* S:E D:E O:E", flags=adj, wid=101), loc(20_400, "R:P* S:E D:E O:E", flags=adj, wid=101),
            loc(30_000, "R:P* S:D D:E O:E", flags=adj, wid=100),                                # a D in wid 100
            loc(40_000, "R:P* S:E D:E O:E", length=9000),                                       # host
            loc(41_000, "R:P* S:E D:E O:E", length=3000, parent=5, flags={"S": ["host_empty"]}),  # its child
            loc(42_000, "R:P* S:E D:E O:E", parent=6, flags={"S": ["host_empty"]})]              # a grandchild
    cat = write_run(tmp_path, loci, fasta_dir=toy_fastas).load()
    ev = {m: v for (_, _, _, m), v in _table(mechanisms(cat, None)[0], "events", "R>S").items()}
    assert (ev["n_E_single"], ev["n_E_cluster"], ev["n_E_cluster_elements"]) == (0, 3, 7)
    assert (ev["n_D_single"], ev["n_D_multi"], ev["n_D_multi_elements"]) == (0, 1, 1)


@pytest.mark.parametrize("damage, match", [
    ("unknown_state", "unknown state 'X'"),
    ("text_dl", "are not integers"),
    ("duplicate", "called more than once"),
    ("missing", "have no call"),
    ("truncated", "cannot read the pair file"),
])
def test_damaged_pair_file_fails_fast(tmp_path, toy_fastas, damage, match):
    """Each fail-fast path of the pair-file read names the file."""
    run = write_run(tmp_path, [loc(1000, "R:P* S:D D:E O:E"), loc(2000, "R:P* S:E D:E O:E")], fasta_dir=toy_fastas)
    cat = run.load()
    path = run.outdir / "pairs" / "R__S.calls.tsv.gz"
    if damage == "truncated":
        data = path.read_bytes()
        path.write_bytes(data[:len(data) // 2])
    else:
        _rewrite(path, {"unknown_state": lambda rows: [{**rows[0], "state": "X"}, *rows[1:]],
                        "text_dl": lambda rows: [{**rows[0], "dl": "x"}, *rows[1:]],      # rows[0]: the D call
                        "duplicate": lambda rows: [*rows, rows[-1]],
                        "missing": lambda rows: rows[:-1]}[damage])
    with pytest.raises(AnalysisError, match=match) as err:
        mechanisms(cat, None)
    assert "R__S.calls.tsv.gz" in str(err.value)


def test_pair_tables_stats_and_findings(tmp_path, toy_fastas):
    loci = [loc(1000, "R:P* S:E D:P* O:E", k=0.0), loc(2000, "R:P* S:D D:U O:E", evidence={"S": (-12, 0, 0)}),
            loc(3000, "R:P* S:U D:E O:E", k=None)]
    cat = write_run(tmp_path, loci, fasta_dir=toy_fastas).load()
    rows, stat, f = mechanisms(cat, None)
    assert {r[1] for r in rows} == {f"{a}>{b}" for a in "RSD" for b in "RSDO" if b != a}   # every A annotated
    assert {r.test for r in stat} == {"mechanisms", "absence_events", "missingness"}
    assert _table(rows, "states", "S>R")[("all", "all", "all", "n_E")] == 0              # S has no element here
    assert [t for t in _lines(f) if t.startswith("S>R")] == ["S>R: no call (the element table of S is empty)"]
    st = _table(rows, "states", "R>S")   # a family: nonzero counts only; a rollup: every metric; no K: all only
    assert ("Gypsy", "famA", "K=0", "n_U") not in st and st[("Gypsy", "famA", "all", "n_U")] == 1
    assert st[("Gypsy", "all", "K=0", "n_U")] == 0 and st[("all", "all", "K=0", "n_G")] == 0
    got = {(r.stratum, r.statistic): r for r in stat if (r.test, r.scope) == ("mechanisms", "R>S")}
    assert got[("all", "frac_D_inside")].value == 0.5 and got[("all", "frac_D_inside")].n == 2
    assert got[("all", "n_absent")].n == 3 and ("K=0", "frac_E") in got and ("0<K<0.01", "frac_E") not in got
    miss = {r.statistic: r.value for r in stat if (r.test, r.scope) == ("missingness", "R>S")}
    assert miss["resolved_S"] == pytest.approx(1 / 3) and miss["U_no_anchor"] == pytest.approx(1 / 3)
    line = next(t for t in _lines(f) if t.startswith("R>S: "))
    assert "1 - frac_E = 0.5 (an upper bound on derived absences)" in line
    assert "frac_D_inside + frac_O_inside = 0.5 + 0 = 0.5 (a lower bound" in line
    assert any(t.startswith("R>S absence events") for t in _lines(f))
    assert any(t.startswith("R>S missing: U 0.333") and "no_anchor" in t for t in _lines(f))


def test_strata_named_all_count_once(tmp_path, toy_fastas):
    loci = [loc(1000, "R:P* S:E D:E O:E", family="all", superfamily="all"), loc(2000, "R:P* S:P D:E O:E")]
    cat = write_run(tmp_path, loci, fasta_dir=toy_fastas).load()
    st = _table(mechanisms(cat, None)[0], "states", "R>S")
    assert (st[("all", "all", "all", "n_E")], st[("all", "all", "all", "n_P")]) == (1, 1)
    rows, stat, _ = _dollo(cat, TREE)
    assert _table(rows, "branches", "R")[("all", "all", "all", "n_insertion")] == 1
    assert {r.n for r in stat if r.scope != "incongruent" and r.stratum == "all"} == {2}


def test_age_strata_without_k(tmp_path, toy_fastas):
    loci = [loc(1000 * (i + 1), "R:P* S:E D:E O:E", k=None) for i in range(3)]
    cat = write_run(tmp_path, loci, fasta_dir=toy_fastas).load()
    rows, _, f = mechanisms(cat, None)
    ab = _table(rows, "absence", "R>S")
    assert {age for _, _, age, _ in ab} == {"all"} and ab[("all", "all", "all", "frac_E")] == 1.0
    assert any(t.startswith("R>S: ") and "no K2P" in t for t in _lines(f))


def test_dollo_losses(tmp_path, toy_fastas):
    loci = [loc(1000, "R:S S:S D:P* O:E"),                        # insertion R+S+D; one solo on R+S
            loc(2000, "R:T S:U D:P* O:E"),                        # insertion R+S+D; truncation on R+S (S missing)
            loc(3000, "R:D S:G D:P* O:P", superfamily="Copia"),   # insertion before_root; deletion on R+S
            loc(4000, "R:D S:S D:P* O:E")]                        # insertion R+S+D; solo on R+S (R's D is
    #                                                               missing for S events), deletion on R
    cat = write_run(tmp_path, loci, fasta_dir=toy_fastas).load()
    rows, stat, _ = _dollo(cat, TREE)
    assert _table(rows, "branches", "R+S+D")[("all", "all", "all", "n_insertion")] == 3
    assert _table(rows, "branches", "R+S+D")[("Gypsy", "all", "all", "n_insertion")] == 3
    assert _table(rows, "branches", "before_root")[("Copia", "all", "all", "n_insertion")] == 1
    rs = _table(rows, "branches", "R+S")
    assert [rs[("all", "all", "all", m)] for m in ("n_insertion", "n_solo", "n_truncation", "n_deletion")] \
        == [0, 2, 1, 1]
    assert _table(rows, "branches", "R")[("all", "all", "all", "n_deletion")] == 1
    assert _table(rows, "branches", "incongruent") == {("all", "all", "all", "n_loci"): 0}
    assert sorted((r.scope, r.stratum, r.statistic, r.value) for r in stat if r.test == "branches") \
        == sorted((sc, sf, m, v) for t, sc, sf, _, _, m, v in rows)
    _, stat, f = mechanisms(cat, TREE)
    assert any(t.startswith("Dollo placement") and "R+S 0/2/1/1" in t for t in _lines(f))
    assert {r.statistic for r in stat if r.test == "branches_resolution"} == {f"resolved_{g}" for g in "RSDO"}


def test_dollo_needs_outgroup_and_rooted_tree(tmp_path, toy_fastas):
    cat = write_run(tmp_path / "a", [loc(1000, "R:P* S:P D:E O:E")], fasta_dir=toy_fastas).load()
    rows, _, f = mechanisms(cat, None)
    assert not _table(rows, "branches", "R+S") and any("Dollo placement skipped: no tree" in t for t in _lines(f))
    with pytest.raises(ValueError, match="rooted"):
        _dollo(cat, Node(None, (Node("O"), Node("D"), Node(None, (Node("R"), Node("S"))))))
    with pytest.raises(ValueError, match="lacks genome"):
        _dollo(cat, Node(None, (Node("O"), Node(None, (Node("R"), Node("S"))))))
    genomes = [g for g in GENOMES if g[0] != "O"]
    cat = write_run(tmp_path / "b", [loc(1000, "R:P* S:P D:E")], genomes=genomes, fasta_dir=toy_fastas).load()
    tree = Node(None, (Node("D"), Node(None, (Node("R"), Node("S")))))
    rows, _, f = mechanisms(cat, tree)
    assert not _table(rows, "branches", "R+S") and any("Dollo placement skipped: no outgroup" in t
                                                       for t in _lines(f))
    with pytest.raises(ValueError, match="outgroup"):
        _dollo(cat, tree)


def test_pair_calls(tmp_path):
    good, bad, short = (tmp_path / f"{x}.calls.tsv.gz" for x in ("good", "bad", "short"))
    row = [str(i) for i in range(len(CALL_COLUMNS))]
    for path, text in ((good, ["\t".join(CALL_COLUMNS), "\t".join(row)]), (bad, ["\t".join(CALL_COLUMNS[:-1])]),
                       (short, ["\t".join(CALL_COLUMNS), "\t".join(row[:-1])])):
        with gzip.open(path, "wt") as fh:
            fh.write("\n".join(text) + "\n")
    assert list(pair_calls(good)) == [dict(zip(CALL_COLUMNS, row))]
    with pytest.raises(AnalysisError, match="CALL_COLUMNS"):
        list(pair_calls(bad))
    with pytest.raises(AnalysisError, match="line 2"):
        list(pair_calls(short))
    with pytest.raises(AnalysisError, match="cannot read"):
        list(pair_calls(tmp_path / "missing.calls.tsv.gz"))


def test_element_table_must_match_the_calls(tmp_path, toy_fastas):
    run = write_run(tmp_path, [loc(1000, "R:P* S:E D:E O:E"), loc(2000, "R:P* S:E D:E O:E")],
                    fasta_dir=toy_fastas)
    cat = run.load()
    path = run.outdir / "work" / "elements" / "R.tsv"
    path.write_text("".join(path.read_text().splitlines(keepends=True)[:-1]))   # one element fewer
    with pytest.raises(AnalysisError, match="not an element of R"):
        mechanisms(cat, None)
