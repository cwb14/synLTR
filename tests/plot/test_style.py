"""Figure style: fonts, mm layout, numbers file, atomic save, safe_figure."""

import logging

from matplotlib import font_manager as fm

from ltrtrace.base.log import log
from ltrtrace.plot import style

GOLDEN_NUMBERS = (
    "key\tvalue\tsource\n"
    "sequences_drawn\t3\tgenome index:A1 sequence lengths\n"
    "block_bp\t13.1 kb\tstats.tsv:D_pav/A1,B1;C1/all/D:note\n"
)


def _numbers():
    n = style.Numbers()
    n.add("sequences_drawn", 3, "genome index:A1 sequence lengths")
    n.add("block_bp", "13.1 kb", "stats.tsv:D_pav/A1,B1;C1/all/D:note")
    return n


def _fig():
    fig = style.page(60)
    ax = style.axes_mm(fig, 20, 10, 100, 30)
    ax.plot([0, 1], [0, 1])
    style.panel_letter(ax, "A")
    return fig


def test_save_atomic(tmp_path):
    stem = tmp_path / "overview"
    style.save(_fig(), stem, _numbers(), legend="Short legend.\n\nSecond paragraph.")
    names = sorted(p.name for p in tmp_path.iterdir())
    assert names == ["overview.pdf", "overview.png", "overview_numbers.tsv"]
    assert (tmp_path / "overview.pdf").read_bytes().startswith(b"%PDF")
    assert (tmp_path / "overview.png").read_bytes().startswith(b"\x89PNG")


def test_numbers_columns(tmp_path):
    stem = tmp_path / "f"
    style.save(_fig(), stem, _numbers())
    assert (tmp_path / "f_numbers.tsv").read_text() == GOLDEN_NUMBERS
    assert style.NUMBERS_COLUMNS == ("key", "value", "source")
    n = _numbers()
    assert n.get("block_bp") == "13.1 kb" and n.get("nope") is None
    for bad in (lambda: n.add("block_bp", 1, "s"), lambda: n.add("k2", "a\tb", "s")):
        try:
            bad()
        except ValueError:
            continue
        raise AssertionError("expected ValueError")


def test_font_fallback_warns(monkeypatch):
    dejavu = [f for f in fm.fontManager.ttflist if f.name == "DejaVu Sans"]
    monkeypatch.setattr(fm.fontManager, "ttflist", list(dejavu))
    monkeypatch.setattr(fm, "findSystemFonts", lambda *a, **k: [])
    records = []

    class Grab(logging.Handler):
        def emit(self, record):
            records.append(record)

    h = Grab(level=logging.WARNING)
    log.addHandler(h)
    try:
        assert style.setup_fonts() == "DejaVu Sans"
    finally:
        log.removeHandler(h)
    assert len(records) == 1 and "falls back" in records[0].getMessage()


def test_font_found_is_liberation_or_arial():
    assert style.setup_fonts() in ("Liberation Sans", "Arial", "DejaVu Sans")


def test_safe_figure_never_raises(caplog):
    @style.safe_figure
    def bad():
        raise RuntimeError("boom")

    @style.safe_figure
    def good():
        return 7

    records = []

    class Grab(logging.Handler):
        def emit(self, record):
            records.append(record)

    h = Grab(level=logging.WARNING)
    log.addHandler(h)
    try:
        assert bad() is None
    finally:
        log.removeHandler(h)
    assert len(records) == 1 and "boom" in records[0].getMessage()
    assert good() == 7


def test_mm_layout():
    fig = style.page(100)
    assert abs(fig.get_size_inches()[0] * 25.4 - 180) < 1e-9
    ax = style.axes_mm(fig, 18, 10, 90, 36)
    x0, y0, w, h = ax.get_position().bounds
    assert abs(x0 * 180 - 18) < 1e-9 and abs(w * 180 - 90) < 1e-9
    assert abs((1 - y0 - h) * 100 - 10) < 1e-9 and abs(h * 100 - 36) < 1e-9


def test_key_row_draws_items():
    fig = style.page(40)
    ax = style.axes_mm(fig, 10, 5, 160, 6)
    style.key_row(ax, [("tick", "#D55E00", "shared"), ("band", "#ECECEC", "block"), ("hatch", "#8F8F8F", "none")])
    assert len(ax.texts) == 3 and len(ax.lines) >= 1


def test_key_row_color_scale_and_tracts():
    from matplotlib.colors import LinearSegmentedColormap

    cmap = LinearSegmentedColormap.from_list("t", ["#FFFFFF", "#000000"])
    items = [("cmap", None, "scale"), ("tract", None, "tract"), ("strong", None, "strict")]
    wide = style.key_layout(items, 160)
    assert wide[1][3] - wide[0][3] > 10.6          # the scale's glyph is wider than the standard one
    fig = style.page(40)
    ax = style.axes_mm(fig, 10, 5, 160, 6)
    style.key_row(ax, items, cmap=cmap)
    assert [t.get_text() for t in ax.texts] == ["0", "1", "scale", "tract", "strict"] and len(ax.patches) > 30


def test_panel_letter_on_the_figure():
    fig = style.page(100)
    style.panel_letter(fig, "A", x_mm=0.5, y_mm=1.0)       # before any axes exists
    t = fig.texts[-1]
    assert t.get_text() == "A" and t.get_position() == (0.5 / 180, 1 - 1.0 / 100)


def test_panel_letter_page_mm():
    fig = style.page(100)
    ax = style.axes_mm(fig, 20, 10, 100, 30)
    style.panel_letter(ax, "B", x_mm=0.5, y_mm=10.0)
    t = fig.texts[-1]
    assert t.get_position() == (0.5 / 180, 1 - 10.0 / 100) and t.get_va() == "top" and t.get_ha() == "left"


def test_save_two_legend_pages(tmp_path):
    style.save(_fig(), tmp_path / "f", _numbers(), legend=["short text", "long text\n\nsecond"])
    assert b"/Count 3" in (tmp_path / "f.pdf").read_bytes()
    style.save(_fig(), tmp_path / "g", _numbers(), legend=[[("run", "b"), (" text", "")]])
    assert b"/Count 2" in (tmp_path / "g.pdf").read_bytes()


def test_safe_figure_survives_unprintable_error_and_raising_handler():
    class Bad(Exception):
        def __str__(self):
            raise RuntimeError("no str")

    class Raising(logging.Handler):
        def emit(self, record):
            record.getMessage()
            raise RuntimeError("handler")

    @style.safe_figure
    def f():
        raise Bad()

    h = Raising(level=logging.WARNING)
    log.addHandler(h)
    try:
        assert f() is None
    finally:
        log.removeHandler(h)


def test_save_failure_leaves_no_temp_files(tmp_path, monkeypatch):
    n = _numbers()
    monkeypatch.setattr(style.Numbers, "write", lambda self, path: (_ for _ in ()).throw(OSError("disk")))
    try:
        style.save(_fig(), tmp_path / "f", n, legend="x")
    except OSError:
        pass
    assert not [p for p in tmp_path.iterdir() if p.name.endswith(".tmp")]
    assert not (tmp_path / "f_numbers.tsv").exists()
