"""Figure scaffolding shared by every figure: fonts, a page laid out in millimetres, panel letters, key rows,
the numbers file, legend pages and an atomic save. Coordinates on a page are mm from its top-left corner."""

from __future__ import annotations

import contextlib
import functools
import os
import re
from collections.abc import Callable, Sequence
from pathlib import Path

from synltr.base import constants as C
from synltr.base.checkpoint import atomic_path
from synltr.base.log import log

NUMBERS_COLUMNS = C.NUMBERS_COLUMNS
FS, FS_SMALL, FS_LETTER = C.FONT_SIZE, C.FONT_SIZE_SMALL, C.FONT_SIZE_LETTER
WIDTH_MM = C.FIG_WIDTH_MM
DPI = C.FIG_DPI

_FONTS = {"sans": "DejaVu Sans", "mono": "DejaVu Sans Mono"}


# ---------------------------------------------------------------- fonts


def setup_fonts() -> str:
    """The sans family of the figures: Liberation Sans or Arial, else DejaVu Sans (bundled with matplotlib);
    Liberation Mono, else DejaVu Sans Mono for the mono family. One WARNING names any fallback (conda-forge has
    no Liberation font package: the host provides it). The font files come from matplotlib's font list and the
    system's font folders, and the list is rebuilt from them in one fixed order, the newest version of each font
    first. Matplotlib takes the first of equal matches, and its cached font list and the system's folders come
    in no fixed order, so two installed versions would otherwise make the figure depend on the matplotlib
    version, its cache and the run."""
    from matplotlib import font_manager as fm

    want = ("Liberation Sans", "Arial", "Liberation Mono")
    files = {f.fname for f in fm.fontManager.ttflist if f.name in want}
    files.update(p for p in fm.findSystemFonts()
                 if os.path.basename(p).startswith(("LiberationSans-", "LiberationMono-", "Arial.", "arial.", "Arial-",
                                                    "Arial ")))
    if files:
        fm.fontManager.ttflist[:] = [f for f in fm.fontManager.ttflist if f.name not in want]
    for p in sorted(files, key=lambda p: (tuple(-v for v in _font_version(p)), p)):
        try:
            fm.fontManager.addfont(p)
        except (OSError, RuntimeError, ValueError) as exc:  # an unreadable font file: skip it
            log.debug("figure: font %s skipped (%s)", p, exc)
    have = {f.name for f in fm.fontManager.ttflist}
    sans = next((f for f in want[:2] if f in have), None)
    mono = "Liberation Mono" if "Liberation Mono" in have else None
    fallback = ([] if sans else ["DejaVu Sans for Liberation Sans or Arial"]) + (
        [] if mono else ["DejaVu Sans Mono for Liberation Mono"])
    if fallback:
        log.warning("figure: falls back to %s, which this host lacks (install the Liberation fonts, "
                    "e.g. the fonts-liberation package; conda-forge has none)", " and ".join(fallback))
    _FONTS["sans"], _FONTS["mono"] = sans or "DejaVu Sans", mono or "DejaVu Sans Mono"
    return _FONTS["sans"]


def _font_version(path: str) -> tuple[int, ...]:
    """A font file's version (its name table's version string, e.g. "Version 2.1.5" -> (2, 1, 5)); () when it
    cannot be read."""
    from matplotlib.ft2font import FT2Font

    try:
        names = FT2Font(path).get_sfnt()
    except (OSError, RuntimeError, ValueError):
        return ()
    for (platform, _, _, name_id), raw in sorted(names.items()):
        if name_id == 5:
            text = raw.decode("utf-16-be" if platform in (0, 3) else "latin-1", errors="ignore")
            mm = re.search(r"\d+(?:\.\d+)*", text)
            if mm:
                return tuple(int(x) for x in mm.group(0).split("."))
    return ()


def rc_params() -> dict:
    """matplotlib rc settings of every figure, for the fonts chosen by `setup_fonts` (call it first); use as
    `with matplotlib.rc_context(rc_params()):`."""
    sans, mono = _FONTS["sans"], _FONTS["mono"]
    return {**C.RC_BASE, "font.sans-serif": [sans, "DejaVu Sans"], "mathtext.rm": sans,
            "mathtext.it": f"{sans}:italic", "mathtext.bf": f"{sans}:bold", "mathtext.sf": sans,
            "mathtext.cal": sans, "mathtext.tt": mono}


# ---------------------------------------------------------------- the page


def page(height_mm: float, width_mm: float = WIDTH_MM):
    """An empty white figure `width_mm` wide and `height_mm` tall."""
    from matplotlib.figure import Figure

    fig = Figure(figsize=(width_mm / 25.4, height_mm / 25.4))
    fig.patch.set_facecolor("white")
    return fig


def _size_mm(fig) -> tuple[float, float]:
    w, h = fig.get_size_inches()
    return w * 25.4, h * 25.4


def axes_mm(fig, x: float, y: float, w: float, h: float, **kw):
    """Axes with their top-left corner at (x, y) mm from the figure's top-left corner, `w` x `h` mm."""
    W, H = _size_mm(fig)
    ax = fig.add_axes((x / W, 1 - (y + h) / H, w / W, h / H), **kw)
    ax.patch.set_alpha(0)
    return ax


def panel_letter(ax, letter: str, *, x_mm: float | None = None, y_mm: float | None = None) -> None:
    """A bold panel letter, top-left anchored. By default above-left of `ax`'s top-left corner at the same offset
    for every panel; with `x_mm` and `y_mm` (mm from the page's top-left corner) it sits exactly there, and the
    first argument may be the figure itself."""
    fig = ax if x_mm is not None and y_mm is not None and not hasattr(ax, "get_position") else ax.figure
    W, H = _size_mm(fig)
    if x_mm is None or y_mm is None:
        box = ax.get_position()
    x = box.x0 - C.LETTER_DX_MM / W if x_mm is None else x_mm / W
    y = box.y1 + C.LETTER_DY_MM / H if y_mm is None else 1 - y_mm / H
    fig.text(x, y, letter, fontsize=FS_LETTER, fontweight="bold", ha="left", va="top")


def bare(ax) -> None:
    """Hide every spine, tick and tick label of `ax`."""
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])


_TTP: list = []


def text_mm(s: str, size: float = FS) -> float:
    """Width of a string in the figure's sans font, in mm."""
    from matplotlib.font_manager import FontProperties
    from matplotlib.textpath import TextToPath

    if not _TTP:
        _TTP.append(TextToPath())
    w = _TTP[0].get_text_width_height_descent(s, FontProperties(family=_FONTS["sans"], size=size),
                                              ismath=False)[0]
    return w * 25.4 / 72


# ---------------------------------------------------------------- key rows



def _glyph_mm(kind: str, color: str | None) -> float:
    """Width of a key glyph: a role's bold label, a color scale, else the standard glyph."""
    if kind == "role":
        return text_mm(color or "", FS) + C.KEY_ROLE_PAD_MM
    return C.KEY_CMAP_MM if kind == "cmap" else C.KEY_GLYPH_MM


def key_layout(items: list[tuple[str, str | None, str]], width: float) -> list[tuple[str, str | None, str, float, int]]:
    """(kind, color, text, x, line) of each key item, wrapped to `width` mm. A "role" item's color field holds its
    bold label and it starts on its own line(s); roles come first."""
    out = []
    x, line, prev = 0.0, 0, "role"
    for kind, color, text in items:
        glyph = _glyph_mm(kind, color)
        need = glyph + C.KEY_TEXT_GAP_MM + text_mm(text)
        if (x > 0 and x + need > width) or (prev == "role" and kind != "role"):
            x, line = 0.0, line + 1
        prev = kind
        out.append((kind, color, text, x, line))
        x += need + (C.KEY_ROLE_GAP_MM if kind == "role" else C.KEY_GAP_MM)
    return out


def key_row(ax, items: list[tuple[str, str | None, str]], cmap=None) -> None:
    """Draw a key in `ax` (axes of an `axes_mm` page; it is turned into a mm frame, y downward): items are
    (kind, color, text) with kind "tick" (a vertical tick), "band" (a filled bar), "hatch" (a hatched bar) or
    "role" (bold label in the color field), wrapped to the axes' width line by line. With a colormap `cmap`, the
    kinds "cmap" (its scale from 0 to 1), "tract" (a bar at its 0.35 with a dark line on top) and "strong" (the
    same with a heavier black line and a mark above) are drawn too."""
    from matplotlib.patches import Rectangle

    W, H = _size_mm(ax.figure)
    box = ax.get_position()
    width, height = box.width * W, box.height * H
    ax.set_xlim(0, width)
    ax.set_ylim(height, 0)
    bare(ax)
    for kind, color, text, x, line in key_layout(items, width):
        yr = line * C.KEY_LINE_MM + C.KEY_LINE_MM / 2
        glyph = _glyph_mm(kind, color)
        if kind == "role":
            ax.text(x, yr, color, va="center", ha="left", fontsize=FS, fontweight="bold")
        elif kind == "tick":
            ax.plot([x + C.KEY_TICK_X_MM] * 2, [yr - C.KEY_TICK_HALF_MM, yr + C.KEY_TICK_HALF_MM], color=color,
                    lw=C.KEY_TICK_LW)
        elif kind == "cmap" and cmap is not None:
            g0, gw = x + C.KEY_CMAP_PAD_MM, glyph - 2 * C.KEY_CMAP_PAD_MM
            for k in range(C.KEY_CMAP_STEPS):
                ax.add_patch(Rectangle((g0 + k * gw / C.KEY_CMAP_STEPS, yr - C.KEY_HATCH_HALF_MM),
                                       gw / C.KEY_CMAP_STEPS + 0.01, C.KEY_HATCH_H_MM,
                                       fc=cmap(k / (C.KEY_CMAP_STEPS - 1)), ec="none"))
            ax.add_patch(Rectangle((g0, yr - C.KEY_HATCH_HALF_MM), gw, C.KEY_HATCH_H_MM, fill=False, ec=C.GRAY,
                                   lw=0.3))
            ax.text(g0 - 0.4, yr, "0", ha="right", va="center", fontsize=FS_SMALL)
            ax.text(g0 + gw + 0.4, yr, "1", ha="left", va="center", fontsize=FS_SMALL)
        elif kind in ("tract", "strong") and cmap is not None:
            ax.add_patch(Rectangle((x, yr - C.KEY_HATCH_HALF_MM), glyph, C.KEY_HATCH_H_MM, fc=cmap(0.35),
                                   ec="none"))
            ax.plot([x, x + glyph], [yr - C.KEY_HATCH_HALF_MM] * 2, color=C.DARK if kind == "tract" else "black",
                    lw=C.KEY_TRACT_LW if kind == "tract" else C.KEY_STRONG_LW, solid_capstyle="butt")
            if kind == "strong":
                ax.scatter([x + glyph / 2], [yr - 1.35], marker="v", s=6, color="black", lw=0)
        elif kind == "band":
            ax.add_patch(Rectangle((x, yr - C.KEY_BAND_HALF_MM), glyph, C.KEY_BAND_H_MM, fc=color, ec="none"))
        elif kind == "hatch":
            ax.add_patch(Rectangle((x, yr - C.KEY_HATCH_HALF_MM), glyph, C.KEY_HATCH_H_MM, fc="white", ec=C.GRAY,
                                   lw=0.0, hatch=C.KEY_HATCH))
        else:
            raise ValueError(f"unknown key item kind {kind!r}")
        ax.text(x + glyph + C.KEY_TEXT_GAP_MM, yr, text, va="center", ha="left", fontsize=FS)


# ---------------------------------------------------------------- the numbers file


class Numbers:
    """Every number written on a figure or in its legend: (key, the text as written, its source)."""

    def __init__(self) -> None:
        self.rows: list[tuple[str, str, str]] = []
        self._keys: set[str] = set()

    def add(self, key: str, value, source: str) -> str:
        """Record `value` (as the text written) under `key` with its source ("<file>:<column>" or
        "<file>:<test>/<scope>/<stratum>/<statistic>"); returns the text."""
        text = str(value)
        if key in self._keys:
            raise ValueError(f"numbers key {key!r} added twice")
        if "\t" in text or "\n" in text or "\t" in source:
            raise ValueError(f"numbers row {key!r}: a tab or newline in {text!r} or {source!r}")
        self._keys.add(key)
        self.rows.append((key, text, source))
        return text

    def get(self, key: str) -> str | None:
        """The text registered under key, or None."""
        return next((text for k, text, _ in self.rows if k == key), None)

    def write(self, path: Path) -> None:
        with open(path, "w", encoding="utf-8", newline="") as fh:
            fh.write("\t".join(NUMBERS_COLUMNS) + "\n")
            fh.writelines("\t".join(r) + "\n" for r in self.rows)


# ---------------------------------------------------------------- legend pages

Paragraphs = list[list[tuple[str, str]]]


def _paragraphs(text: str | Paragraphs) -> Paragraphs:
    if isinstance(text, str):
        return [[(p.strip(), "")] for p in re.split(r"\n\s*\n", text) if p.strip()]
    return text


def text_page(paras: Paragraphs, family: str | None = None):
    """A legend page as wide as the figure, its height fitted to the text: 10 pt, ragged right, each run's style
    ("" normal, "b" bold, "i" italic) kept, every line a few text objects (searchable)."""
    from matplotlib.figure import Figure
    from matplotlib.font_manager import FontProperties
    from matplotlib.textpath import TextToPath

    family = family or _FONTS["sans"]
    size, lead, para_gap = C.LEGEND_FONT_PT, C.LEGEND_LEAD_PT, C.LEGEND_PARA_GAP_PT
    width_pt = WIDTH_MM / 25.4 * 72
    margin = C.LEGEND_MARGIN_MM / 25.4 * 72
    usable = width_pt - 2 * margin
    props = {"": FontProperties(family=family, size=size),
             "b": FontProperties(family=family, size=size, weight="bold"),
             "i": FontProperties(family=family, size=size, style="italic")}
    ttp = TextToPath()
    cache: dict[tuple[str, str], float] = {}

    def width(s: str, st: str) -> float:
        k = (s, st)
        if k not in cache:
            cache[k] = ttp.get_text_width_height_descent(s, props[st], ismath=False)[0]
        return cache[k]

    room = usable * C.LEGEND_LINE_FILL

    def pieces(w: str, st: str) -> list[tuple[str, bool]]:
        """A word as (piece, glued to the one before): whole when it fits a line, else broken after "/" or "_"
        (a piece still too wide is cut by characters)."""
        if width(w, st) <= room:
            return [(w, False)]
        out: list[str] = []
        for part in re.findall(r"[^/_]*[/_]|[^/_]+$", w):
            while width(part, st) > room and len(part) > 1:
                k = len(part) - 1
                while k > 1 and width(part[:k], st) > room:
                    k -= 1
                out.append(part[:k])
                part = part[k:]
            out.append(part)
        return [(x, k > 0) for k, x in enumerate(out) if x]

    lines: list[list[tuple[float, str, str]] | None] = []
    for para in paras:
        words = [(piece, st, glue) for text, st in para for w in text.split() for piece, glue in pieces(w, st)]
        line: list[tuple[float, str, str]] = []
        x = 0.0
        for w, st, glue in words:
            ww = width(w, st)
            sp = 0.0 if glue or not line else width(" ", st)
            if line and x + sp + ww > room:
                lines.append(line)
                line, x, sp = [], 0.0, 0.0
            if line and line[-1][2] == st:
                x0, txt, _ = line[-1]
                line[-1] = (x0, f"{txt}{'' if glue else ' '}{w}", st)
            else:
                line.append((x + sp, w, st))
            x += sp + ww
        if line:
            lines.append(line)
        lines.append(None)
    if lines and lines[-1] is None:
        lines.pop()
    height = 2 * margin + sum(lead if ln is not None else para_gap for ln in lines) + C.LEGEND_PAD_PT
    fig = Figure(figsize=(width_pt / 72, height / 72))
    fig.patch.set_facecolor("white")
    y = margin + size
    for ln in lines:
        if ln is None:
            y += para_gap
            continue
        for x, txt, st in ln:
            fig.text((margin + x) / width_pt, 1 - y / height, txt, fontproperties=props[st], va="baseline",
                     ha="left")
        y += lead
    return fig


def legend_page(pdf, text: str | Paragraphs) -> None:
    """Append a legend page to an open `PdfPages`. `text` is paragraphs separated by blank lines, or a list of
    paragraphs, each a list of (text, style) runs with style "" / "b" / "i"."""
    pdf.savefig(text_page(_paragraphs(text)))


# ---------------------------------------------------------------- saving


def save(fig, stem: Path, numbers: Numbers, legend: str | Paragraphs | Sequence[str | Paragraphs] | None = None) -> None:
    """Write `<stem>.png` (figure only), `<stem>.pdf` (the figure, then one text page per entry of `legend`) and
    `<stem>_numbers.tsv`, all three to temporary files first, renamed into place together at the end: a drawing or a
    write that fails leaves the figure's files as they were. `legend` is one page's text or a list of pages (a page:
    paragraphs separated by blank lines, or a list of paragraphs of (text, style) runs)."""
    from matplotlib.backends.backend_pdf import PdfPages

    if legend is None or isinstance(legend, str) or (legend and isinstance(legend[0], list)
                                                     and legend[0] and isinstance(legend[0][0], tuple)):
        pages = [legend] if legend else []
    else:
        pages = list(legend)
    with contextlib.ExitStack() as stack:  # each temporary file is renamed once the block ends without an error
        png, pdf, tsv = (stack.enter_context(atomic_path(path)) for path in figure_files(stem))
        fig.savefig(png, dpi=DPI, format="png")
        with PdfPages(pdf, metadata={"Title": Path(stem).name, "Creator": "synltr"}) as pp:
            pp.savefig(fig)
            for text in pages:
                legend_page(pp, text)
        numbers.write(tsv)


def figure_files(stem: Path) -> list[Path]:
    """The files `save` writes for `stem`: <stem>.png, <stem>.pdf and <stem>_numbers.tsv (C.FIGURE_SUFFIXES)."""
    stem = Path(stem)
    return [stem.with_name(f"{stem.name}{end}") for end in C.FIGURE_SUFFIXES]


def safe_figure(fn: Callable) -> Callable:
    """Decorator: a figure never fails a run. Any exception is logged as one warning and the call returns None."""

    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        try:
            return fn(*args, **kwargs)
        except Exception as exc:  # noqa: BLE001 - a figure failure must never fail the run
            try:
                log.warning("figure %s failed: %s: %s", getattr(fn, "__name__", "?"), type(exc).__name__, exc)
            except Exception:  # noqa: BLE001, S110 - an unprintable error or a failing handler stays silent
                pass
            return None

    return wrapper
