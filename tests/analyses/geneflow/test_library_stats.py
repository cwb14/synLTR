"""The library screen's statistics (analyses.geneflow.library_stats): LTRmix's trio tests on its synthetic trios (the
signal of gene flow when admixed, none without), the f-ratio with a second donor (D's sister E, never an outgroup),
every statistic's numbers against LTRmix's own (stored outputs, on seeded inputs and on whole trios, records without
coordinates included), the every-pair evenness check, the gate, and inferred absences never counted as empty sites."""
import math
import re
import shutil
from dataclasses import replace
from pathlib import Path

import builders
import numpy as np
import pytest

from ltrtrace.analyses.geneflow import library_stats as ls
from ltrtrace.analyses.geneflow.library_stats import LibraryRun, compute, evenness, gate
from ltrtrace.analyses.geneflow.route import ROUTE_PROPS
from ltrtrace.base.constants import EVEN_CUTOFF, FALSE_SIGNALS
from ltrtrace.base.log import AnalysisError
from ltrtrace.calls.records import read_divergences, stream_calls, write_calls, write_divergences
from ltrtrace.calls.states import OBSERVED_ABSENT
from ltrtrace.io.elements import NO_LTR, Element, ElementTable
from ltrtrace.io.samples import Sample, read_samples
from ltrtrace.routes.library.run import screen_pairs
from ltrtrace.stats.report import StatRow
from ltrtrace.stats.resample import block_size

GOLDEN = Path(__file__).resolve().parents[2] / "golden" / "ltrmix_stats.npz"
PROPS = ROUTE_PROPS["library"]
TRIO = ("R", "S", "D")
NEW_TESTS = {"evenness"}  # rows LTRmix never wrote
EXTRA = {"triples", "BBAA", "median_d", "n_target_recent", "n_mirror_recent"}  # statistics LTRmix printed, not tabled
ZERO_DEN = ("its denominator is 0: among the donor's records resolved in the recipient, the sister and the second donor, "
            "those present in the second donor and not the sister are as many as those present in the sister and not "
            "the second donor")
REPLICATE = "deleting one jackknife block leaves the denominator 0"


@pytest.fixture(scope="module")
def stored():
    with np.load(GOLDEN) as z:
        return {k: z[k] for k in z.files}


def _run(sc, d):
    screen_pairs(read_samples(sc.sheet), d / "run", 2)
    return LibraryRun(d / "run", tuple(read_samples(sc.sheet)), 2, {})


def _lib(screened) -> LibraryRun:
    sc, run = screened
    return LibraryRun(run, tuple(read_samples(sc.sheet)), 2, {})


@pytest.fixture(scope="module")
def admixed_run(admixed):
    return _lib(admixed)


@pytest.fixture(scope="module")
def null_run(null_screen):
    return _lib(null_screen)


@pytest.fixture(scope="module")
def four_run(four_screen):
    return _lib(four_screen)


@pytest.fixture(scope="module")
def outgroup_run(tmp_path_factory):
    d = tmp_path_factory.mktemp("outgroup")
    return _run(builders.trio_outgroup(d / "in", seed=3), d)


@pytest.fixture(scope="module")
def unplaced_run(unplaced_screen):
    """The admixed trio with 4 records whose names hold no coordinates (mutated copies of placed records, no table
    row) appended to R's library and 4 to D's."""
    return _lib(unplaced_screen)


@pytest.fixture(scope="module")
def unplaced_k2p_run(tmp_path_factory):
    """The admixed trio (seed 3) with 40 records whose names hold no coordinates appended to S's library and 40 to
    D's (random sequences), each with a table row carrying a K2P above the genomes' usual range (0.2-0.3)."""
    d = tmp_path_factory.mktemp("unplaced_k2p")
    sc = builders.trio(d / "in", seed=3, admix=True)
    rng = np.random.default_rng(29)
    lines = ["label\tlibrary"]
    for g in TRIO:
        fa, tb = sc.libraries[g].split(",")
        if g != "R":
            names = [f"{g}_ctg{k}" for k in range(40)]
            seqs = ["".join("ACGT"[c] for c in rng.integers(0, 4, int(rng.integers(1000, 3000)))) for _ in names]
            fa2, tb2 = d / "in" / f"{g}k.fa", d / "in" / f"{g}k.tsv"
            fa2.write_text(Path(fa).read_text() + "".join(f">{n}\n{s}\n" for n, s in zip(names, seqs)))
            tb2.write_text(Path(tb).read_text() + "".join(f"{n}\t{rng.uniform(0.2, 0.3):.5f}\t.\t+\t250\t250\n"
                                                          for n in names))
            fa, tb = fa2, tb2
        lines.append(f"{g}\t{fa},{tb}")
    sheet = d / "in" / "unplaced_k2p.tsv"
    sheet.write_text("\n".join(lines) + "\n")
    screen_pairs(read_samples(sheet), d / "run", 2)
    return LibraryRun(d / "run", tuple(read_samples(sheet)), 2, {})


def _row(rows, test, stratum="all", statistic="D"):
    found = [r for r in rows if (r.test, r.stratum, r.statistic) == (test, stratum, statistic)]
    assert len(found) == 1, (test, stratum, statistic)
    return found[0]


def _equal(x, y) -> bool:
    return (math.isnan(x) and math.isnan(y)) if (math.isnan(x) or math.isnan(y)) else x == y


def _ltrmix_values(rows) -> dict:
    """{(test, stratum, key): value} of the rows in LTRmix's stats.tsv form: each statistic under its own key, the
    se of an estimate (D, ratio, the insertion f-ratio's alpha) as SE and the z of a D as Z."""
    out = {}
    for r in rows:
        if r.test in NEW_TESTS or r.statistic in EXTRA:
            continue
        out[r.test, r.stratum, r.statistic] = r.value
        if r.statistic in ("D", "ratio") or (r.test, r.statistic) == ("f_ratio_ins", "alpha"):
            out[r.test, r.stratum, "SE"] = r.se
        if r.statistic == "D":
            out[r.test, r.stratum, "Z"] = r.z
    return out


def _stored_rows(stored, case) -> dict:
    keys = zip(*(stored[f"e2e_{case}_{x}"].tolist() for x in ("test", "stratum", "key", "value")))
    return {(t, s, k): v for t, s, k, v in keys}


# ---------------------------------------------------------------- LTRmix's trio tests


def test_D_positive_when_admixed(admixed_run):
    """R carries D's haplotype in three tracts: D_kmer and D_ins are positive and significant."""
    rows = compute(admixed_run, TRIO, PROPS)
    for test in ("D_kmer", "D_ins"):
        r = _row(rows, test)
        assert r.value > 0 and r.z > 3, (test, r.value, r.z)


def test_D_null_when_not_admixed(null_run):
    """Without gene flow neither D_kmer nor D_ins departs from 0."""
    rows = compute(null_run, TRIO, PROPS)
    for test in ("D_kmer", "D_ins"):
        r = _row(rows, test)
        assert abs(r.z) < 3, (test, r.value, r.z)


def test_four_genomes_fratio(four_run, outgroup_run):
    """A fourth genome on D's side of the tree (E, D's sister) is the second donor: the f-ratio rows are reported for
    it, alpha finite. An outgroup (closer to none of the trio) is never a second donor: no f-ratio row."""
    rows = compute(four_run, TRIO, PROPS)
    f = [r for r in rows if r.test.startswith("f_ratio")]
    assert {(r.test, r.scope, r.stratum) for r in f} == {("f_ratio_ins", "R,S;D;E", "second donor E"),
                                                         ("f_ratio_identity", "R,S;D;E", "second donor E")}
    assert np.isfinite(_row(f, "f_ratio_ins", "second donor E", "alpha").value)
    assert not [r for r in compute(outgroup_run, TRIO, PROPS) if r.test.startswith("f_ratio")]


def test_outgroup_is_never_a_second_donor(four_screen, four_run):
    """A genome whose sheet role is outgroup is no second donor even where its 1:1 orthologs meet the rule (the trios
    leave outgroups out alike): with E, D's sister, as the outgroup the trio gets no f-ratio row for E, while E in the
    ingroup, the same genome in the same position, gets them."""
    sc, run_dir = four_screen
    sheet = tuple(replace(s, role="outgroup") if s.label == "E" else s for s in read_samples(sc.sheet))
    out = LibraryRun(run_dir, sheet, 2, {})
    assert ls.second_donors(four_run, *TRIO) == ["E"] and ls.second_donors(out, *TRIO) == []
    assert not [r for r in compute(out, TRIO, PROPS) if r.test.startswith("f_ratio")]
    assert {r.scope for r in compute(four_run, TRIO, PROPS) if r.test.startswith("f_ratio")} == {"R,S;D;E"}


# ---------------------------------------------------------------- LTRmix's numbers


@pytest.mark.parametrize("case", ["admixed", "null", "four", "unplaced"])
def test_rows_match_ltrmix(case, request, stored):
    """Every row LTRmix wrote for the trio (its Trio.patterns, divergences and fratio through its own orthology),
    value for value and in its order, from the screen's files and the trio LTRmix inferred; the jackknife block
    from every genome on the sheet. `unplaced` holds records without coordinates in R and D, which count in the
    arrays' denominators, in their own jackknife blocks and in the placed span (LTRmix's block size moves)."""
    run = request.getfixturevalue(f"{case}_run")
    rows = compute(run, tuple(stored[f"e2e_{case}_roles"].tolist()), PROPS)
    mine, want = _ltrmix_values(rows), _stored_rows(stored, case)
    assert set(mine) == set(want)
    assert [k for k in want if not _equal(mine[k], want[k])] == []
    groups = list(dict.fromkeys((r.test, r.stratum) for r in rows if r.test not in NEW_TESTS))
    assert groups == list(dict.fromkeys(k[:2] for k in want))
    assert ls.block_bp(run, PROPS) == int(stored[f"e2e_{case}_block_bp"])


def test_block_bp_counts_records_without_coordinates(stored):
    """The placed span sums each sequence's last element end; a record whose name holds no coordinates is its own
    sequence, ending at its length: A's 600-bp contig makes B the smallest genome (LTRmix's set_block_size)."""
    spans = [ls.placed_span(stored[f"bs_{g}_chrom"].astype(object), stored[f"bs_{g}_end"]) for g in "AB"]
    assert spans == [1_000_600, 1_000_500]
    assert block_size(spans) == int(stored["bs_out"]) == 10_005
    placed = ls.placed_span(stored["bs_A_chrom"][:2].astype(object), stored["bs_A_end"][:2])
    assert block_size([placed, spans[1]]) == int(stored["bs_out_placed_only"]) == 10_000


def test_d_stat_matches_ltrmix(stored):
    """The D-form delete-one block jackknife, (sum a - sum b) / (sum a + sum b), on LTRmix's text block labels: its
    estimate, se, z and block count, on 17 seeded cases (fewer than 3 blocks, nothing informative, one-sided counts,
    every informative row in one block among them)."""
    for k in range(int(stored["jk_n"])):
        a, b, blocks = stored[f"jk{k}_a"], stored[f"jk{k}_b"], stored[f"jk{k}_blocks"]
        jk = ls.d_stat(a, b, blocks)
        est, se, z, n_a, n_b, nb = stored[f"jk{k}_out"].tolist()
        assert (_equal(jk.est, est), _equal(jk.se, se), _equal(jk.z, z), jk.nblocks) == (True, True, True, nb), k
        assert (a.sum(), b.sum()) == (n_a, n_b)


def test_completeness_matches_ltrmix(stored):
    """The completeness ratio (present in R over present in S, at least 1 below) and its jackknife se."""
    for k in range(int(stored["cp_n"])):
        jk = ls.completeness(stored[f"cp{k}_pr"], stored[f"cp{k}_ps"], stored[f"cp{k}_blocks"])
        assert (jk.est, jk.se) == tuple(stored[f"cp{k}_out"].tolist()), k


def test_f_ratio_matches_ltrmix(stored):
    """The f-ratio from insertion patterns (alpha, se, num, den, blocks) and from identity-spectrum excesses. A zero
    denominator gives alpha NaN with its reason and se NaN (LTRmix's alpha was NaN too; its se came from the
    pseudovalues, NaN in this case as a block's own denominator is 0, finite only when the denominator cancels across
    blocks none of whose own denominators is 0: test_f_ratio_cancelling_denominator); a jackknife replicate whose
    denominator is 0 gives se NaN with its reason (LTRmix's se NaN: its pseudovalues hold that replicate's inf)."""
    for k in range(int(stored["fr_n"])):
        st = [stored[f"fr{k}_{x}"] for x in ("st_dr", "st_ds", "st_d2")]
        d = [stored[f"fr{k}_{x}"] for x in ("dR", "dS", "d2")]
        f = ls.f_ratio(*st, stored[f"fr{k}_blocks"], *d)
        alpha, se, num, den, nb, alpha_id, e_r, e_2 = stored[f"fr{k}_out"].tolist()
        assert (f.num, f.den, f.ins.nblocks) == (num, den, nb), k
        assert (_equal(f.identity, alpha_id), f.excess_r, f.excess_2) == (True, e_r, e_2), k
        if den:
            assert f.ins.est == alpha and _equal(f.ins.se, se), k
            assert f.why == (None if math.isfinite(se) else REPLICATE), k
        else:
            assert math.isnan(alpha) and math.isnan(f.ins.est) and math.isnan(f.ins.se)
            assert f.why == ZERO_DEN


def test_f_ratio_zero_denominator():
    """D2 present exactly where S is: the denominator is 0, so alpha is NaN with the reason (never inf or 0); a
    denominator held by one jackknife block alone leaves alpha but no se, with the reason."""
    st_r, st_s = np.array([2, 2, 1, 2, 1, 2]), np.array([1, 2, 2, 1, 2, 2])
    blocks, d = np.array([0, 0, 1, 1, 2, 2]), np.full(6, np.nan)
    f = ls.f_ratio(st_r, st_s, st_s, blocks, d, d, d)
    assert (f.num, f.den) == (0, 0) and math.isnan(f.ins.est) and math.isnan(f.ins.se) and math.isnan(f.identity)
    assert f.why == ZERO_DEN
    st_2 = st_s.copy()
    st_2[0] = 2
    f = ls.f_ratio(st_r, st_s, st_2, blocks, d, d, d)
    assert (f.num, f.den, f.ins.est) == (0, 1, 0.0) and math.isnan(f.ins.se) and f.why == REPLICATE


def test_f_ratio_cancelling_denominator():
    """The denominator cancels across blocks (+1, +1, -2: the second donor and the sister differ in both directions,
    equally often): alpha NaN and its se NaN, a standard error of an undefined estimate (LTRmix printed one from its
    pseudovalues, 1.333 here, as no block's own denominator is 0; with one that is, its se was NaN too), with the
    reason true for this case too."""
    st_r, st_s, st_2 = np.array([2, 2, 2, 2]), np.array([1, 1, 2, 2]), np.array([2, 2, 1, 1])
    f = ls.f_ratio(st_r, st_s, st_2, np.array([0, 1, 2, 2]), *[np.full(4, np.nan)] * 3)
    assert (f.num, f.den, f.ins.nblocks) == (2, 0, 3) and math.isnan(f.ins.est) and math.isnan(f.ins.se)
    assert f.why == ZERO_DEN


def test_f_ratio_note_says_what_is_nan():
    """The f-ratio row's note: alpha NaN gives ": NaN, <reason>"; alpha finite with its se NaN gives "; se NaN,
    <reason>"; otherwise the definition alone."""
    jk = ls.Jk(math.nan, math.nan, math.nan, 3)
    assert ls._ratio_note("alpha", ls.FRatio(jk, 2, 0, ZERO_DEN, math.nan, 0.0, 0.0)) == f"alpha: NaN, {ZERO_DEN}"
    jk = ls.Jk(0.25, math.nan, math.nan, 3)
    assert ls._ratio_note("alpha", ls.FRatio(jk, 1, 4, REPLICATE, 1.0, 0.1, 0.1)) == f"alpha; se NaN, {REPLICATE}"
    jk = ls.Jk(0.25, 0.1, 2.5, 3)
    assert ls._ratio_note("alpha", ls.FRatio(jk, 1, 4, None, 1.0, 0.1, 0.1)) == "alpha"


def test_findings_give_the_f_ratio_reason(admixed_run):
    """The f-ratio's findings line gives the reason in full: after alpha NaN, and after a finite alpha whose se alone
    is NaN (se NaN, then the reason)."""
    rows = compute(admixed_run, TRIO, PROPS)
    for value, se, note, said in ((math.nan, math.nan, f"definition: NaN, {ZERO_DEN}", f"; {ZERO_DEN})"),
                                  (0.25, math.nan, f"definition; se NaN, {REPLICATE}", f"; se NaN, {REPLICATE})"),
                                  (0.25, 0.1, "definition", "(insertion patterns), ")):
        extra = [StatRow("library", "f_ratio_ins", "R,S;D;E", "second donor E", "alpha", value, se, note=note),
                 StatRow("library", "f_ratio_identity", "R,S;D;E", "second donor E", "alpha", 0.5)]
        lines = [t for _, t in ls.findings(rows + extra).items() if t.startswith("Admixture fraction")]
        assert len(lines) == 1 and said in lines[0], (note, lines)


def test_age_bins_on_the_edges():
    """D_ins's age strata: a K2P on an edge falls in the bin below it (right-closed bins: 0, (0, 0.002], ...), above
    0.1 or unknown in the last."""
    k = np.array([0.0, 1e-9, 2e-9, 0.002, 0.0020001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.1000001, math.nan])
    assert ls.age_bins(k).tolist() == [0, 0, 1, 1, 2, 2, 3, 4, 5, 6, 7, 7]


def test_completeness_needs_more_than_100_old_insertions():
    """The completeness ratio's old insertions: K2P above the 0.75 quantile of every record's finite K2P, resolved in
    both sisters; reported only with more than COMPLETENESS_MIN (100): exactly 100 give none, 101 give all."""
    for n, want in ((400, 0), (404, 101)):
        k = np.linspace(0.001, 0.2, n)
        k_old, old = ls.completeness_set(k, np.full(n, 2), np.full(n, 1))
        assert old.size == want and (k[old] > k_old).all(), n
    k_old, old = ls.completeness_set(np.full(10, math.nan), np.full(10, 2), np.full(10, 2))
    assert k_old == math.inf and old.size == 0


def test_absence_certain_needs_room_under_half_the_element():
    """An absence is certain where the syntenic interval is shorter than half the element: room equal to half is not;
    no room (-1) is not."""
    room = np.array([-1, 0, 499, 500, 501])
    assert ls.absence_certain(room, np.full(5, 1000)).tolist() == [False, True, True, False, False]


def _counts_run(counts: dict[str, int], outgroups: str = "") -> LibraryRun:
    """A run of genomes R, S, D, E whose 1:1 ortholog counts are given ("XY" -> count), those in `outgroups` with
    role outgroup, no files."""
    run = LibraryRun(Path("unused"), tuple(Sample(g, None, None, None, "outgroup" if g in outgroups else None,
                                                  library="x") for g in "RSDE"), 1, {})
    run._orthologs.update({frozenset(k): v for k, v in counts.items()})
    return run


def test_second_donors_rule():
    """A second donor has strictly more 1:1 orthologs with D than D has with R or S, and than it has itself with R or
    S: an exact tie fails, and so does a candidate that fails only on its own side. An outgroup never is one."""
    trio = {"RS": 500, "RD": 400, "SD": 380}
    assert ls.second_donors(_counts_run({**trio, "DE": 600, "RE": 390, "SE": 370}), "R", "S", "D") == ["E"]
    assert ls.second_donors(_counts_run({**trio, "DE": 400, "RE": 390, "SE": 370}), "R", "S", "D") == []
    assert ls.second_donors(_counts_run({**trio, "DE": 600, "RE": 600, "SE": 370}), "R", "S", "D") == []
    assert ls.second_donors(_counts_run({**trio, "DE": 600, "RE": 390, "SE": 650}), "R", "S", "D") == []
    assert ls.second_donors(_counts_run({**trio, "DE": 600, "RE": 390, "SE": 370}, "E"), "R", "S", "D") == []


def test_identity_spectrum_matches_ltrmix(stored):
    """Per-bin counts and calibrated ratios, the admixed-fraction lower bound, the calibration, the median and the
    recent counts of the identity spectrum."""
    for k in range(int(stored["is_n"])):
        sp = ls.spectrum(stored[f"is{k}_t"], stored[f"is{k}_m"])
        assert np.array_equal(sp.n_target, stored[f"is{k}_ht"]) and np.array_equal(sp.n_mirror, stored[f"is{k}_hm"])
        assert np.array_equal(sp.ratio, stored[f"is{k}_ratio"], equal_nan=True)
        f, nt, nm, cal, med, lo_t, lo_m = stored[f"is{k}_out"].tolist()
        assert (sp.f, len(stored[f"is{k}_t"]), len(stored[f"is{k}_m"])) == (f, nt, nm)
        assert (sp.calibration, sp.median, sp.recent_target, sp.recent_mirror) == (cal, med, lo_t, lo_m)


def test_trio_second_donor_is_ltrmix_E(tmp_path, stored):
    """builders.trio_second_donor writes LTRmix's fourth genome E (D's library lightly mutated, renamed) byte for
    byte, LTRmix's sha256s stored."""
    import hashlib

    sc = builders.trio_second_donor(tmp_path, seed=5)
    got = [hashlib.sha256((tmp_path / f).read_bytes()).hexdigest() for f in ("E.fa", "E.tsv")]
    assert got == stored["e2e_four_E_sha256"].tolist()
    assert [s.label for s in read_samples(sc.sheet)] == ["R", "S", "D", "E"]


# ---------------------------------------------------------------- the pair files


def test_partner_without_coordinates_stops_the_run(admixed_run, tmp_path):
    """A P call whose ortholog has no coordinates (b_eid ".") or names no placed record of the partner stops the
    statistics with one error naming the pair and the call."""
    run = LibraryRun(tmp_path / "run", admixed_run.samples, 2, admixed_run.records)
    shutil.copytree(admixed_run.outdir, run.outdir)
    path = run.outdir / "pairs" / "R__D.div.tsv.gz"
    rows = list(read_divergences(path))
    eid = rows[3]["eid"]
    for bad, why in ((None, "no coordinates"), ("D|chrX:1-2", "no placed record")):
        rows[3]["b_eid"] = bad
        write_divergences(path, rows, {"engine": "library", "d_kind": "p_whole_unit"})
        with pytest.raises(AnalysisError, match=f"R__D.*{re.escape(eid)}.*{why}"):
            ls.pair_calls(run, "R", "D")


def test_pair_arrays_cover_every_record(unplaced_run):
    """The arrays run over every record in library order: a record whose name holds no coordinates has no call and
    stays status 0, ortholog -1, room -1, divergence NaN."""
    c = ls.pair_calls(unplaced_run, "R", "D")
    r = unplaced_run.records["R"]
    u = ~r.placed
    assert u.sum() == 4 and c.orth.size == r.n
    assert (c.status()[u] == 0).all() and (c.orth[u] == -1).all() and (c.room[u] == -1).all()
    assert np.isnan(c.div[u]).all() and (c.status()[r.placed] > 0).mean() > 0.9


# ---------------------------------------------------------------- absences


def test_inferred_absence_never_counts_as_e(admixed_run, monkeypatch):
    """The screen's absences are inferred (A). A statistic validated on observed empty sites (E) finds no absence it
    accepts among them: its value is NaN and its note says why. The statistics validated on A are unchanged."""
    c = ls.pair_calls(admixed_run, "R", "S")
    assert (c.status(OBSERVED_ABSENT) == 1).sum() == 0 and (c.status() == 1).sum() > 100
    base = compute(admixed_run, TRIO, PROPS)
    monkeypatch.setitem(ls.ABSENCES, "D_ins", OBSERVED_ABSENT)
    rows = compute(admixed_run, TRIO, PROPS)
    d_ins = [r for r in rows if r.test == "D_ins"]
    assert d_ins and all(math.isnan(r.value) and math.isnan(r.se) for r in d_ins)
    assert all("inferred absences (A)" in r.note and "never count as observed empty sites (E)" in r.note
               for r in d_ins)
    assert _row(rows, "D_kmer").value == _row(base, "D_kmer").value
    assert _row(rows, "D3_triples").value == _row(base, "D3_triples").value
    monkeypatch.setitem(ls.ABSENCES, "evenness", OBSERVED_ABSENT)  # the check itself: NaN, so D3 and D_self are held
    rows = compute(admixed_run, TRIO, PROPS)
    ratios = [r for r in rows if r.test == "evenness" and r.statistic == "r"]
    assert len(ratios) == 3 and all(math.isnan(r.value) and "never count as observed empty sites" in r.note
                                    for r in ratios)
    assert gate(rows)["D3_triples"] == "held: uneven annotation"


# ---------------------------------------------------------------- the evenness check and the gate


def _elements(label, ks, mids):
    return ElementTable.from_rows(
        Element(genome=label, eid=f"{label}|c1:{m - 50}-{m + 50}", name=f"e{i}", seqid="c1", start=m - 50,
                end=m + 50, strand="+", ltr_left=NO_LTR, ltr_right=NO_LTR, tsd=None, k=k, family=None,
                superfamily=None, clade=None, confidence=None)
        for i, (k, m) in enumerate(zip(ks, mids)))


def test_evenness_by_hand(tmp_path):
    """Two genomes, 16 elements each with K2P (one unknown) and 4 old ones (K2P above the 0.75 quantile), one per
    1-Mb block. X's old ones in Y: P P A U (2 of 3 resolved present); Y's in X: P A A P (2 of 4). r = (2/3)/(2/4);
    se of log r from each genome's own delete-one block jackknife, the variances added."""
    young = [0.001 * i for i in range(1, 12)] + [math.nan]
    ks = young + [0.05, 0.06, 0.07, 0.08]
    mids = [100 + 1000 * i for i in range(12)] + [500_000, 1_500_000, 2_500_000, 3_500_000]
    els = {"X": _elements("X", ks, mids), "Y": _elements("Y", ks, mids)}
    states = {"X": ["P"] * 12 + ["P", "P", "A", "U"], "Y": ["A"] * 12 + ["P", "A", "A", "P"]}
    pairs = tmp_path / "pairs"
    for x, y in (("X", "Y"), ("Y", "X")):
        write_calls(pairs / f"{x}__{y}.calls.tsv.gz",
                    [{"eid": e, "state": s} for e, s in zip(els[x].column("eid"), states[x])], {"engine": "library"})
    rows = evenness(pairs, els)

    def jk_var(present):  # delete-one jackknife variance of log(share), one resolved element per block
        p = np.array([math.log(np.delete(present, i).mean()) for i in range(len(present))])
        return (len(p) - 1) / len(p) * ((p - p.mean()) ** 2).sum()

    r = _row(rows, "evenness", "old insertions", "r")
    assert r.scope == "X,Y" and r.value == pytest.approx((2 / 3) / (2 / 4), rel=1e-12) and r.n == 7
    log_r = _row(rows, "evenness", "old insertions", "log_r")
    assert log_r.value == pytest.approx(math.log(4 / 3), rel=1e-12)
    assert log_r.se == pytest.approx(math.sqrt(jk_var(np.array([1, 1, 0.])) + jk_var(np.array([1, 0, 0, 1.]))),
                                     rel=1e-12)
    assert (_row(rows, "evenness", "old insertions", "share_X_in_Y").value,
            _row(rows, "evenness", "old insertions", "share_Y_in_X").value) == pytest.approx((2 / 3, 0.5))
    check = _row(rows, "evenness", "all", "check")
    assert check.value == pytest.approx(0.75) and "X,Y" in check.note


def _write_states(pairs, els, states):
    for x, y in (("X", "Y"), ("Y", "X")):
        write_calls(pairs / f"{x}__{y}.calls.tsv.gz",
                    [{"eid": e, "state": s} for e, s in zip(els[x].column("eid"), states[x])], {"engine": "library"})


def test_evenness_quantile_tie(tmp_path):
    """An old insertion's K2P must be above the 0.75 quantile: X's five K2P put the quantile exactly on 0.04, so only
    0.05 is old (present in Y): share 1 of 1, not 1 of 2 with the absent 0.04 one. Y's old one is absent in X."""
    ks = [0.01, 0.02, 0.03, 0.04, 0.05]
    mids = [500_000 + 1_000_000 * i for i in range(5)]
    els = {g: _elements(g, ks, mids) for g in "XY"}
    _write_states(tmp_path, els, {"X": ["P", "P", "P", "A", "P"], "Y": ["P", "P", "P", "P", "A"]})
    rows = evenness(tmp_path, els)
    share = _row(rows, "evenness", "old insertions", "share_X_in_Y")
    assert (share.value, share.n) == (1.0, 1) and _row(rows, "evenness", "old insertions", "share_Y_in_X").n == 1


def test_evenness_note_without_k2p(tmp_path):
    """A genome without any known K2P has no old insertion: r NaN, and the note says no K2P is known."""
    mids = [500_000 + 1_000_000 * i for i in range(5)]
    els = {"X": _elements("X", [math.nan] * 5, mids), "Y": _elements("Y", [0.01, 0.02, 0.03, 0.04, 0.05], mids)}
    _write_states(tmp_path, els, {"X": ["P"] * 5, "Y": ["P"] * 5})
    r = _row(evenness(tmp_path, els), "evenness", "old insertions", "r")
    assert math.isnan(r.value) and r.note.endswith("NaN: no K2P known for X")


def _pair_evenness_share(run, x, y, every):
    """pair_evenness.py's share of x's old insertions present in y and its count: old = finite K2P above the 0.75
    quantile of x's finite K2P over every record (`every`) or over its placed records; resolved = P or A."""
    rec = run.records[x]
    k = rec.k2p if every else rec.k2p[rec.placed]
    k_old = np.quantile(k[np.isfinite(k)], 0.75)
    old = {e for e, kk in zip(rec.lib.elements.column("eid").tolist(), rec.k2p) if e is not None and kk > k_old}
    path = run.outdir / "pairs" / f"{x}__{y}.calls.tsv.gz"
    states = [s for _, (e, s) in stream_calls(path, ("eid", "state")) if e in old and s in ("P", "A")]
    return sum(s == "P" for s in states) / len(states), len(states)


def test_evenness_quantile_over_every_record(unplaced_k2p_run):
    """The old-insertion quantile is taken over every record's finite K2P, records whose names hold no coordinates
    included (pair_evenness.py's rule on LTRmix's orthologs tables): r and n of each pair as that rule gives them on
    the same calls; over the placed records alone the counts would differ."""
    run = unplaced_k2p_run
    rows = compute(run, TRIO, PROPS)
    for x, y in (("R", "S"), ("R", "D"), ("S", "D")):
        (s1, n1), (s2, n2) = (_pair_evenness_share(run, a, b, True) for a, b in ((x, y), (y, x)))
        r = next(q for q in rows if (q.test, q.scope, q.statistic) == ("evenness", f"{x},{y}", "r"))
        assert r.value == pytest.approx(s1 / s2, rel=1e-12) and r.n == n1 + n2, (x, y)
    placed = [_pair_evenness_share(run, g, h, False)[1] for g, h in (("S", "D"), ("D", "S"))]
    assert placed != [_pair_evenness_share(run, g, h, True)[1] for g, h in (("S", "D"), ("D", "S"))]


def _even(r):
    return [StatRow("library", "evenness", p, "old insertions", "r", v)
            for p, v in (("A,B", r), ("A,C", 1.01), ("B,C", 0.99))]


def test_gate_holds():
    """D_ins and D_kmer are always held (annotation-sensitive); D3 and D_self are held when any pair's ratio or its
    inverse is below EVEN_CUTOFF (0.889), or when the check has no value."""
    sensitive = {t: "held: annotation-sensitive" for t in
                 ("D_kmer", "D_kmer_by_chrom", "D_ins", "D_ins_by_age", "D_ins_by_chrom")}
    uneven = {"D3_triples": "held: uneven annotation", "D_self_outgroup": "held: uneven annotation"}
    assert gate(_even(0.66)) == sensitive | uneven
    assert gate(_even(1 / 0.66)) == sensitive | uneven
    assert gate(_even(0.95)) == sensitive
    assert gate(_even(EVEN_CUTOFF)) == sensitive
    assert gate(_even(math.nan)) == sensitive | uneven
    assert gate([]) == sensitive | uneven


def test_hold_replaces_an_earlier_flag():
    """hold() on rows already carrying a flag: the gate's flag replaces it, and a statistic the gate no longer holds
    loses it."""
    rows = _even(0.95) + [StatRow("library", "D3_triples", "R,S;D", "all", "D", 0.1, note="held: uneven annotation; d"),
                          StatRow("library", "D_kmer", "R,S;D", "all", "D", 0.1, note="held: uneven annotation; k"),
                          StatRow("library", "D_ins", "R,S;D", "all", "D", 0.1, note="held: annotation-sensitive")]
    notes = [r.note for r in ls.hold(rows)[3:]]
    assert notes == ["d", "held: annotation-sensitive; k", "held: annotation-sensitive"]


def test_held_statistics_flagged_with_their_rates(admixed_run):
    """compute's rows: the evenness check of the trio's three pairs; every row of a held statistic carries its gate
    note; each headline statistic carries its simulated false-signal rates; the completeness ratio is a value only."""
    rows = compute(admixed_run, TRIO, PROPS)
    assert {r.scope for r in rows if r.test == "evenness" and r.statistic == "r"} == {"R,S", "R,D", "S,D"}
    notes = gate(rows)
    for r in rows:
        if r.test in notes:
            assert r.note.startswith(notes[r.test]), r
        elif r.test in ("D3_triples", "D_self_outgroup"):
            assert not r.note.startswith("held"), r
    for test, stratum in (("D_kmer", "all"), ("D_ins", "all"), ("D3_triples", "all"),
                          ("D_self_outgroup", "LTR sites")):
        note = _row(rows, test, stratum).note
        for condition, (n, of) in FALSE_SIGNALS[test].items():
            assert f"{condition} {n}/{of}" in note, (test, condition)
    ratio = _row(rows, "completeness_ratio", next(r.stratum for r in rows if r.test == "completeness_ratio"),
                 "ratio")
    assert "value only" in ratio.note and "10-point" in ratio.note


def test_findings_give_verdicts_only_for_unheld_statistics(admixed_run):
    """Findings: each statistic with its value and its false-signal rates; D_ins and D_kmer never give a verdict; D3
    and D_self give one when the evenness check passes, else say why not and point to geneflow without --screen."""
    rows = compute(admixed_run, TRIO, PROPS)
    even = gate(rows).get("D3_triples") is None
    text = "\n".join(t for _, t in ls.findings(rows).items())
    assert "D_ins and D_kmer give no verdict" in text and "9/10" in text
    uneven = ls.hold([replace(r, value=0.66) if (r.test, r.statistic) == ("evenness", "r") else r for r in rows])
    held = "\n".join(t for _, t in ls.findings(uneven).items())
    assert "D3 and D_self give no verdict" in held and "geneflow` without --screen" in held
    assert ("D3 and D_self give no verdict" in text) is not even
    assert not [r for r in uneven if r.test in ("D3_triples", "D_self_outgroup") and not r.note.startswith("held")]


def test_duplicate_divergence_row_named(admixed_run, tmp_path):
    """A divergence row repeated for one P call is reported as a duplicate."""
    run = LibraryRun(tmp_path / "run", admixed_run.samples, 2, admixed_run.records)
    shutil.copytree(admixed_run.outdir, run.outdir)
    path = run.outdir / "pairs" / "R__D.div.tsv.gz"
    rows = list(read_divergences(path))
    write_divergences(path, rows[:4] + rows[3:], {"engine": "library", "d_kind": "p_whole_unit"})
    with pytest.raises(AnalysisError, match=f"duplicate divergence row of {re.escape(rows[3]['eid'])}"):
        ls.pair_calls(run, "R", "D")


def test_findings_no_verdict_without_z(admixed_run):
    """A D3 or D_self without a z (not computed) gets no verdict, not "no signal"."""
    rows = [replace(r, z=math.nan) if (r.test, r.statistic) == ("D3_triples", "D") else r
            for r in compute(admixed_run, TRIO, PROPS)]
    assert gate(rows).get("D3_triples") is None
    text = "\n".join(t for _, t in ls.findings(rows).items())
    assert "Verdict, D3: no verdict (z not computed)." in text


def test_calls_read_only_from_finished_library_files(admixed_run, tmp_path):
    """A pair whose call file is missing, or is not a finished library screen's file, stops the statistics with an
    error naming the file."""
    run = LibraryRun(tmp_path / "run", admixed_run.samples, 2, admixed_run.records)
    shutil.copytree(admixed_run.outdir, run.outdir)
    path = run.outdir / "pairs" / "S__R.calls.tsv.gz"
    path.with_name(path.name + ".done.json").unlink()
    with pytest.raises(AnalysisError, match="S__R.calls"):
        ls.pair_calls(run, "S", "R")
    path.unlink()
    with pytest.raises(AnalysisError, match="S__R.calls"):
        ls.pair_calls(run, "S", "R")
