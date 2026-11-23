"""Derived biology: nesting check, insertion histories, clade table, copies after the merger."""
import numpy as np

from synltr.analyses.subgenomes import biology, inputs, merger


def _nest_lib(spans, k2p, seq=None):
    lib = inputs.Library()
    lib.n = len(spans)
    lib.start = np.array([s for s, _ in spans], np.int64)
    lib.end = np.array([e for _, e in spans], np.int64)
    lib.seq = np.zeros(lib.n, np.int64) if seq is None else np.array(seq, np.int64)
    lib.k2p = np.array(k2p, float)
    lib.k2p_se = np.full(lib.n, 0.002)
    lib.n_sites = np.full(lib.n, 500.0)
    lib.ltr5_len = np.full(lib.n, 500.0)
    lib.ltr3_len = np.full(lib.n, 500.0)
    return lib


def test_nesting_check_planted():
    """Host 1000-12000 holds a younger element (no violation) and, inside that, an older one (violation, its host is
    the smallest container); an overlapping-but-not-contained element, a copy at identical coordinates and an element
    on another sequence are not nested; an inner element within a host LTR is flagged."""
    spans = [(1000, 12000), (3000, 9000), (4000, 6000), (11000, 15000), (1000, 12000), (3500, 5000), (1100, 1300)]
    lib = _nest_lib(spans, [0.05, 0.01, 0.03, 0.02, 0.05, 0.001, 0.01], seq=[0, 0, 0, 0, 0, 1, 0])
    res = {(r["inner"], r["host"]): r for r in biology.nesting_check(lib)}
    assert set(res) == {(1, 4), (2, 1), (6, 4)} or set(res) == {(1, 0), (2, 1), (6, 0)}
    assert not res[next(k for k in res if k[0] == 1)]["violation"] and res[(2, 1)]["violation"]
    assert res[next(k for k in res if k[0] == 6)]["in_host_ltr"] and not res[(2, 1)]["in_host_ltr"]


def test_nesting_check_ignores_overlaps_and_missing_ages():
    """Abutting or overlapping elements are not nested; an inner element without K2P is not checked."""
    lib = _nest_lib([(100, 900), (900, 2000), (1500, 3000), (200, 600)], [0.01, 0.02, 0.03, np.nan])
    assert biology.nesting_check(lib) == []


def test_insertion_history_counts():
    """Observed and deconvolved copies per bin each sum to the copies of that subgenome (and clade)."""
    rng = np.random.default_rng(22)
    n = 900
    lib = inputs.Library()
    lib.n = n
    a = 10 ** rng.uniform(-4, -1, n)
    k = rng.binomial(500, 0.75 * (1 - np.exp(-4 * a / 3)))
    lib.n_ts, lib.n_tv, lib.n_sites = k.astype(float), np.zeros(n), np.full(n, 500.0)
    lib.ltr5_len = np.full(n, 500.0)
    lib.k2p = -0.5 * np.log(1 - 2 * k / 500)
    lib.clade = np.array(["Tekay"] * 600 + ["CRM"] * 300)
    loc = np.array([0, 1, 1] * 300)
    rows = biology.insertion_history(lib, loc, 2, rng)
    for g in (0, 1):
        for c in ("all", "Tekay"):
            R = [r for r in rows if r["subgenome"] == g and r["clade"] == c]
            m = (loc == g) & ((lib.clade == c) if c != "all" else True)
            assert sum(r["observed"] for r in R) == m.sum() and abs(sum(r["deconvolved"] for r in R) - m.sum()) < 1e-6
    assert not [r for r in rows if r["clade"] == "CRM" and r["subgenome"] == 0]          # 100 CRM copies: < MIN_HIST


def test_clade_table_counts():
    """Copies, kb, LTR K2P quartiles per clade and subgenome; expected copies after the subgenome's merger (sum of
    probabilities); a missing clade is 'unknown'."""
    lib = _nest_lib([(1, 1000), (1, 2000), (1, 3000), (1, 4000)], [0.01, 0.02, 0.03, 0.04])
    lib.clade = ["Tekay", "Tekay", "", "Tekay"]
    after = {0: np.array([0.9, 0.2, 0.5, np.nan])}            # P(after the merger) of subgenome 0's copies
    rows = {(r["clade"], r["subgenome"]): r for r in biology.clade_table(lib, np.array([0, 0, 0, 1]), 2, after)}
    t0 = rows[("Tekay", 0)]
    assert t0["n"] == 2 and t0["kb"] == 3.0 and t0["q50"] == 0.015 and abs(t0["n_after"] - 1.1) < 1e-12
    assert rows[("unknown", 0)]["n"] == 1 and rows[("Tekay", 1)]["n_after"] is None


def test_merger_cut_needs_a_dated_merger():
    """Copies after the merger need a dated merger: a bound limits the onset of cross-subgenome transposition, not the
    merger, and a not-datable pair has none."""
    assert biology.merger_cut(dict(merger=None)) is None
    assert biology.merger_cut(dict(merger=dict(datable=False, bound=False, estimate=0.01, lo=0.0, hi=0.02))) is None
    assert biology.merger_cut(dict(merger=dict(datable=True, bound=True, estimate=0.001, lo=0.0, hi=0.003))) is None
    assert biology.merger_cut(dict(merger=dict(datable=True, bound=False, estimate=0.002, lo=0.001, hi=0.004))) == 0.002


def test_hist_mass_has_no_grid_aliasing():
    """An age distribution flat in log age gives equal mass in every log-spaced history bin (no bin gains mass from
    holding one more AGE_GRID point than its neighbors); masses sum to 1."""
    g = np.where(merger.AGE_GRID > 0, 1.0, 0.0)
    g /= g.sum()
    m = biology.hist_mass(g)
    assert m.size == biology.HIST_EDGES.size - 1 and abs(m.sum() - 1) < 1e-9
    inner = m[1:-1]
    assert inner.max() / inner.min() < 1.01


def _after_lib(rng, n=3000, T=0.004, n_sites=600):
    p = lambda a: 0.75 * (1 - np.exp(-4 * a / 3))
    age = 10 ** rng.uniform(-4.5, np.log10(0.05), n)
    lib = inputs.Library()
    lib.n = n
    k = rng.binomial(n_sites, p(age))
    lib.n_ts, lib.n_tv, lib.n_sites = k.astype(float), np.zeros(n), np.full(n, float(n_sites))
    lib.ltr5_len, lib.ltr3_len, lib.k2p_se = np.full(n, float(n_sites)), np.full(n, float(n_sites)), np.full(n, np.nan)
    lib.k2p = -0.5 * np.log(1 - 2 * k / n_sites)
    lib.start, lib.end = np.arange(n) * 10000 + 1, np.arange(n) * 10000 + 5000
    lib.clade = ["Tekay"] * n
    lib.cls = [""] * n
    return lib, age


def test_after_merger_is_the_expected_count():
    """Copies after a dated merger: the sum of each copy's probability of being younger than it (its LTR counts, the
    subgenome's age distribution) matches the planted count; a hard K2P cut would count zero-substitution pre-merger
    copies too. A bound gives no count."""
    rng = np.random.default_rng(31)
    lib, age = _after_lib(rng)
    loc = np.array([0, 1] * (lib.n // 2))
    T = 0.004
    mg = dict(datable=True, bound=False, estimate=T, lo=T / 2, hi=2 * T)
    after = biology.after_merger(lib, loc, dict(pairs=[dict(groups=(0, 1), merger=mg)]), rng)
    rows = {r["subgenome"]: r for r in biology.clade_table(lib, loc, 2, after)}
    for g in (0, 1):
        true = int(((age < T) & (loc == g)).sum())
        assert abs(rows[g]["n_after"] - true) < 0.1 * true, (rows[g]["n_after"], true)
    bound = dict(datable=True, bound=True, estimate=T, lo=0.0, hi=T)
    assert biology.after_merger(lib, loc, dict(pairs=[dict(groups=(0, 1), merger=bound)]), rng) == {}


def test_clade_names_from_fasta_class():
    """Without a GFF clade (K2P-table input), the clade is the third field of the FASTA header class
    (LTR/Gypsy/Tekay); 'unknown' when neither names one."""
    lib = inputs.Library()
    lib.clade = ["", "", "Tekay", ""]
    lib.cls = ["LTR/Gypsy/CRM", "LTR/Copia/unknown", "LTR/Gypsy/Athila", ""]
    assert biology.clade_names(lib).tolist() == ["CRM", "unknown", "Tekay", "unknown"]


def test_insertion_history_skips_a_clade_that_is_the_whole_subgenome():
    """A clade holding every copy of a subgenome (e.g. all 'unknown') would repeat the 'all' rows: not written."""
    rng = np.random.default_rng(32)
    lib, _ = _after_lib(rng, n=600)
    rows = biology.insertion_history(lib, np.zeros(lib.n, int), 1, rng)
    assert {r["clade"] for r in rows} == {"all"}


def test_k2p_only_input_history_and_nesting():
    """K2P without LTR substitution counts or lengths: the observed history is written (deconvolved NA) and nested
    pairs are not checked (no standard error), instead of reporting 0% violations."""
    rng = np.random.default_rng(33)
    lib, _ = _after_lib(rng, n=600)
    lib.n_ts = lib.n_tv = lib.n_sites = lib.ltr5_len = lib.ltr3_len = np.full(lib.n, np.nan)
    rows = biology.insertion_history(lib, np.zeros(lib.n, int), 1, rng)
    assert sum(r["observed"] for r in rows if r["clade"] == "all") == lib.n
    assert all(not np.isfinite(r["deconvolved"]) for r in rows)
    nest = _nest_lib([(1000, 12000), (3000, 9000)], [0.01, 0.05])
    nest.k2p_se = nest.n_sites = nest.ltr5_len = nest.ltr3_len = np.full(2, np.nan)
    assert biology.nesting_check(nest) == []
