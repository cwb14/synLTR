"""report.pdf of a subgenomes run: drawn from report_data.json, which keeps every value the page reads; report.png and
report_numbers.tsv beside it."""
import csv
import datetime
import json
import re
import shutil
import subprocess
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from analyses.subgenomes import sgtools as sg

from ltrtrace import cli
from ltrtrace.analyses.subgenomes import report
from ltrtrace.base.log import UsageError
from ltrtrace.calls.manifest import write_manifest
from ltrtrace.plot import run as plot_run
from ltrtrace.plot import subgenomes as plotsg
from ltrtrace.plot.style import Numbers


def test_report_with_and_without_homoeologs(tmp_path):
    """The one-page report draws from stub run objects: the plain page when there are no homoeologs, the B dot plot
    and set numbers when there are."""
    import matplotlib
    matplotlib.use("Agg")
    for with_ho in (False, True):
        a, lib, ph, dt, bio, pa = sg.figure_stub(2, 6)
        if with_ho:
            z, ek = np.asarray(ph["z"]), np.asarray(ph["ek"])
            P_ = np.array([[e, f] for e in range(lib.n) for f in range(e + 1, lib.n)
                           if ek[f] == ek[e] + 1 and z[ek[e]] == 0 and e % 40 == f % 40][:200])
            calls = {(i, 1 - int(z[i])): dict(status="homoeolog", partner=i + 1 if z[i] == 0 else i - 1,
                                              candidate=i + 1 if z[i] == 0 else i - 1, reciprocal=True, p=0.001)
                     for i in range(6)}
            ph["homoeologs"] = dict(P=P_, calls=calls, in_chain=np.ones(len(P_), bool),
                                    pair_called=np.ones(len(P_), bool),
                                    sets=[dict(members=[0, 1], complete=True), dict(members=[2, 3], complete=True),
                                          dict(members=[4, 5], complete=True)])
        a.outdir = str(tmp_path / ("ho" if with_ho else "plain"))
        Path(a.outdir).mkdir()
        report.write_report_data(a, lib, ph, dt, bio, pa)
        plotsg.report(Path(a.outdir))
        out = Path(a.outdir) / "report.pdf"
        assert out.exists() and out.stat().st_size > 10_000
        if with_ho:
            assert report._homoeolog_numbers(ph) == {0: 1, 1: 1, 2: 2, 3: 2, 4: 3, 5: 3}
            shown, _ = report._shown(lib, ph)
            assert [shown.index(x) for x in (0, 2, 4)] == [0, 1, 2]       # homoeologs share a position per group


def test_report_data_keeps_types_and_values(tmp_path):
    """Arrays keep dtype, shape and every value (NaN and infinities too, float32 exactly); tuples, dicts keyed by
    tuples or numbers, namespaces, numpy scalars and None come back as they went in."""
    x = dict(a=np.array([[1.5, np.nan], [np.inf, -0.0]]), b=np.arange(4, dtype=np.int32), c=np.zeros(0, bool),
             d=np.array([0.1, 1 / 3, 7e-9], np.float32), t=(1, "x", None), k={(0, "Tekay"): [1, 2], 3: "three"},
             ns=SimpleNamespace(n=2, seqnames=["s1", "s2"]), s=np.int64(7), f=np.float64(0.25), g=np.bool_(True),
             e=[], m={"nested": {(1, 2): np.array([1, 2], np.int64)}})
    path = tmp_path / "d.json"
    report.dump_data(path, x)
    json.loads(path.read_text(), parse_constant=lambda c: c)               # one JSON document
    y = report.load_data(path)
    assert set(y) == set(x)
    for k in ("a", "b", "c", "d"):
        assert y[k].dtype == x[k].dtype and y[k].shape == x[k].shape
        assert np.array_equal(y[k], x[k], equal_nan=True)
    assert np.signbit(y["a"][1, 1]) and y["t"] == (1, "x", None) and y["k"] == {(0, "Tekay"): [1, 2], 3: "three"}
    assert y["ns"] == x["ns"] and y["s"] == 7 and type(y["s"]) is int and y["f"] == 0.25 and y["g"] is True
    assert y["e"] == [] and np.array_equal(y["m"]["nested"][(1, 2)], [1, 2])


def test_figure_data_carries_what_the_page_reads(tmp_path):
    """The figure data holds the run objects cut to the fields the page reads, and no library sequence."""
    a, lib, ph, dt, bio, pa = sg.figure_stub(3, 12, dating=True)
    a.outdir = str(tmp_path)
    P = report.figure_data(a, lib, ph, dt, bio, pa)
    assert P["pa"] is None and "loc" not in P
    assert set(vars(P["lib"])) == {"n", "seqnames", "seq", "start", "end"}
    assert set(P["dt"]) == {"pairs"} and set(P["bio"]) == {"clades"}
    assert {"names", "z", "ek", "W", "halves", "support", "tree", "source", "verdict", "strata_edges",
            "homoeologs"} == set(P["ph"])


def test_plot_command_redraws_a_subgenomes_run(tmp_path):
    """`ltrtrace plot RUN` on a subgenomes run folder (its manifest.json says so) draws report.pdf again from
    report_data.json."""
    a, lib, ph, dt, bio, pa = sg.figure_stub(2, 8, dating=True)
    a.outdir = str(tmp_path)
    report.write_report_data(a, lib, ph, dt, bio, pa)
    now = datetime.datetime.now(datetime.UTC)
    write_manifest(tmp_path / "manifest.json", ["ltrtrace", "subgenomes", "sheet.tsv"], {"command": "subgenomes"},
                   [], {}, now, now)
    assert plot_run.plot(tmp_path) == 0 and (tmp_path / "report.pdf").stat().st_size > 10_000
    (tmp_path / "report.pdf").unlink()
    assert cli.main(["plot", str(tmp_path)]) == 0 and (tmp_path / "report.pdf").exists()


def test_plot_command_names_a_folder_without_a_run(tmp_path):
    """A folder without manifest.json holds no finished run, figure data or not: a usage error."""
    assert cli.main(["plot", str(tmp_path)]) == 2
    with pytest.raises(UsageError, match="no finished"):
        plot_run.plot(tmp_path)
    a, lib, ph, dt, bio, pa = sg.figure_stub(2, 8, dating=True)
    a.outdir = str(tmp_path)
    report.write_report_data(a, lib, ph, dt, bio, pa)
    assert cli.main(["plot", str(tmp_path)]) == 2 and not (tmp_path / "report.pdf").exists()


def _pdf_pages(path):
    if shutil.which("pdfinfo"):
        r = subprocess.run(["pdfinfo", str(path)], capture_output=True, text=True, check=True)
        return int(re.search(r"^Pages:\s+(\d+)", r.stdout, re.MULTILINE).group(1))
    return len(re.findall(rb"/Type /Page\b(?!s)", path.read_bytes()))


def _numbers(path):
    """report_numbers.tsv -> (header, {key: (value, source)})."""
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    return rows[0], {k: (v, src) for k, v, src in rows[1:]}


@pytest.mark.parametrize("N,S,dating,source", [(2, 8, True, "inferred"), (3, 12, True, "inferred"),
                                               (12, 24, False, "inferred"), (1, 10, False, "inferred"),
                                               (2, 8, False, "user")])
def test_report_renders_one_page(tmp_path, N, S, dating, source):
    """Two to twelve groups, one group, dated or not, a user partition: one page, no exception."""
    a, lib, ph, dt, bio, pa = sg.figure_stub(N, S, dating=dating, source=source)
    a.outdir = str(tmp_path)
    report.write_report_data(a, lib, ph, dt, bio, pa)
    plotsg.report(tmp_path)
    assert _pdf_pages(tmp_path / "report.pdf") == 1


def test_report_keeps_every_label_on_the_page(tmp_path):
    """Scaffold names of 35 characters, three dated subgenomes: every text of the figure stays on the page."""
    a, lib, ph, dt, bio, pa = sg.figure_stub(3, 12, prefix="Bnapus_ZS11_chrA01_RagTag_polished_", dating=True)
    a.outdir = str(tmp_path)
    P = report.load_data(report.write_report_data(a, lib, ph, dt, bio, pa))
    with plotsg._rc():
        fig = plotsg.fig_report(P, Numbers())
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        texts = list(fig.texts) + [t for ax in fig.axes for t in ax.texts + ax.get_yticklabels() + ax.get_xticklabels()
                                   if t.get_visible() and t.get_text()]
        w, h = fig.canvas.get_width_height()
        box = [t.get_window_extent(r) for t in texts]
    assert min(b.x0 for b in box) >= -1 and min(b.y0 for b in box) >= -1
    assert max(b.x1 for b in box) <= w + 1 and max(b.y1 for b in box) <= h + 1


def test_pca_places_groups_apart_and_projects_halves():
    """Panel C's PCA: two planted groups separate on axis 1 (first group on the left); each half-library projects
    one point per sequence."""
    ph = sg.figure_stub(2, 10)[2]
    X, reps, share = plotsg._pca(dict(ph=ph))
    z = ph["z"]
    assert X[z == 0, 0].max() < X[z == 1, 0].min() and share[0] > 0.5
    assert reps.shape == (len(ph["halves"]), 10, 2)


def test_report_png_and_numbers(tmp_path):
    """report.png (the page at 300 dpi) and report_numbers.tsv: every result number printed on the page, as printed,
    with the report_data.json key it comes from (or the constant behind a printed threshold)."""
    import matplotlib.image as mpimg

    a, lib, ph, dt, bio, pa = sg.figure_stub(2, 8, dating=True)
    a.outdir = str(tmp_path)
    report.write_report_data(a, lib, ph, dt, bio, pa)
    plotsg.report(tmp_path)
    assert mpimg.imread(tmp_path / "report.png").shape[:2] == (round(9.45 * 300), round(7.09 * 300))
    pdf = (tmp_path / "report.pdf").read_bytes()      # the page's fonts embedded as TrueType, never as Type 3 glyphs
    assert b"/FontFile2" in pdf and b"/Subtype /Type3" not in pdf
    head, nums = _numbers(tmp_path / "report_numbers.tsv")
    assert head == ["key", "value", "source"]
    assert all(src.startswith(("report_data.json:", "base.constants:")) for _, src in nums.values())
    assert nums["groups"] == ("two", "report_data.json:N") and nums["phased"][0] == "8"
    assert nums["ltr_rts"][0] == "320" and "unphased" not in nums
    for g in (0, 1):                                  # copies with an LTR age per subgenome (panel E's margin)
        n = sum(r["observed"] for r in bio["history"] if r["subgenome"] == g and r["clade"] == "all")
        assert nums[f"E_SG{g + 1}_copies"][0] == f"{n:,.0f}"
        assert nums[f"F_Tekay_SG{g + 1}_copies"][0] == f"{300 + 10 * g:,}"
    assert nums["D_split1"][0] == "SG1 vs SG2" and nums["D_split1_replicability"][0] == "0.95"
    assert nums["D_split1_result"][0] == "kept"
    assert nums["H_split_profile"][0] == "split 1.92 Myr (1.15\u20133.85)"           # K2P / (2 mu), mu 1.3e-8
    assert nums["H_split_inherited"][0] == "not datable"
    assert nums["H_merger_genealogy_clock"][0] == "merger < 77 kyr"                  # a bound: hi only
    assert nums["time_axis_rate"][0] == "13"


_AGE = re.compile(r"([\d.,]+)(?: (Myr|kyr|yr))?")


def _years(text):
    """'split 1.99 Myr (1.8\u20132.2)' -> [estimate, low, high] in years (a limit without a unit takes the
    estimate's)."""
    vals = _AGE.findall(text.split(" ", 1)[1])
    unit = vals[0][1]
    scale = {"Myr": 1e6, "kyr": 1e3, "yr": 1.0}
    return [float(v.replace(",", "")) * scale[u or unit] for v, u in vals]


@pytest.mark.slow
def test_report_numbers_on_the_builders_polyploid(tmp_path):
    """The builders' dated allotetraploid: report_numbers.tsv holds the page's counts and dates as the run's tables give
    them (phased sequences and LTR-RTs from summary.json, copies per subgenome from insertion_history.tsv, the
    reported split from dating.tsv)."""
    d, _ = sg.simulate(tmp_path, tree={"split": 2e6, "children": ["A", "B"]},
                       mergers=[{"time": 1e5, "lineages": ["A", "B"]}])
    sg.run_tool(d, tmp_path / "out", "--mu", "1.3e-8", "-t", "8")
    out = tmp_path / "out"
    _, nums = _numbers(out / "report_numbers.tsv")
    s = json.loads((out / "summary.json").read_text())
    assert nums["groups"][0] == "two" and nums["phased"][0] == str(s["n_phased"])
    assert nums["ltr_rts"][0] == f"{s['n_ltr']:,}"
    H = [ln.split("\t") for ln in (out / "insertion_history.tsv").read_text().splitlines() if not ln.startswith("#")]
    for g in ("SG1", "SG2"):
        n = sum(int(r[H[0].index("observed")]) for r in H[1:] if r[0] == g and r[1] == "all")
        assert nums[f"E_{g}_copies"][0] == f"{n:,}"
    sp = next(r for r in sg.dating_rows(out) if r["event"] == "split" and r["reported"] == "yes")
    shown = _years(nums[f"H_split_{sp['estimator']}"][0])
    want = [float(sp[k]) for k in ("years", "years_low", "years_high")]
    assert shown == pytest.approx(want, rel=0.01)                      # 3 significant digits on the page
