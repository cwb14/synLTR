"""report.pdf of a subgenomes run, drawn from its report_data.json (analyses.subgenomes.report.figure_data). One page,
no legend pages: every panel has a one-line header and its own key, so the page reads alone. A: the painted karyotype
across the top. Left column: B lineage sharing between sequences (or the homoeolog dot plot), C each sequence's
placement in independent half-libraries, D the two tests behind the number of subgenomes. Right column, on one time
axis: E transposon births per subgenome, F per main clade, G sharing in young, middle-aged and old copies, H closest
relatives across subgenomes with every split and merger estimate. Nothing is re-estimated for the figure.

The page is drawn in plot.style's fonts (setup_fonts) under FIG_RC and saved by plot.style.save as report.pdf,
report.png (the page at FIG_DPI) and report_numbers.tsv: every result number printed on the page (values, counts, ages
and their intervals, test results, and the thresholds and rate printed with them; no tick label), as printed, with the
report_data.json key it comes from (or the constant behind a printed threshold)."""
from __future__ import annotations

import math
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

from ltrtrace.analyses.subgenomes.inputs import group_label, natural_key
from ltrtrace.analyses.subgenomes.links import log_oe
from ltrtrace.analyses.subgenomes.report import REPORT_DATA, _homoeolog_numbers, _ma, load_data
from ltrtrace.base.constants import (
    AGE_STRATA,
    FIG_ACCENT,
    FIG_ALPHA,
    FIG_DARK,
    FIG_EST,
    FIG_EXCH,
    FIG_GRAY,
    FIG_H,
    FIG_HALVES_MAX,
    FIG_INK,
    FIG_LIGHT,
    FIG_MERGE,
    FIG_NUM,
    FIG_RC,
    FIG_SPLIT,
    FIG_W,
    FIG_XL,
    FIG_XR,
    MIN_HIST,
    MIN_REL,
    MIN_REP,
    PROGRAM,
    REPORT_STEM,
    SG_COLORS,
    SUPPORT_OK,
)
from ltrtrace.plot import style

STEM = REPORT_STEM  # the page's files: report.png, report.pdf and report_numbers.tsv (plot.style.save)


def report(outdir: Path) -> None:
    """report.pdf (one page), report.png and report_numbers.tsv in the subgenomes run folder outdir, drawn from its
    report_data.json and saved by plot.style.save."""
    outdir = Path(outdir)
    P = load_data(outdir / REPORT_DATA)
    nums = style.Numbers()
    with _rc():
        style.save(fig_report(P, nums, outdir), outdir / STEM, nums, legend=None)


def _rc():
    """The rc context report.pdf is drawn and saved in: FIG_RC with the sans family of plot.style.setup_fonts
    (Liberation Sans or Arial, else DejaVu Sans)."""
    import matplotlib as mpl
    return mpl.rc_context({**FIG_RC, "font.sans-serif": [style.setup_fonts()]})


def _ax(fig, x0, y0, w, h):
    """Axes at x0, y0 (inches from the lower left), w x h inches: every panel sits on a grid measured in inches."""
    W, H = fig.get_size_inches()
    return fig.add_axes([x0 / W, y0 / H, w / W, h / H])


def _head(fig, ax, letter, text, xl, dy=0.06):
    """Bold panel letter at xl inches from the left edge and a one-line header after it, dy inches above the axes."""
    W, H = fig.get_size_inches()
    y = ax.get_position().y1 + dy / H
    fig.text(xl / W, y, letter, fontsize=9, fontweight="bold", va="baseline", ha="left")
    fig.text((xl + 0.17) / W, y, text, fontsize=7, va="baseline", ha="left", color=FIG_INK)


def _text_w(fig, s, size, **kw):
    """Width of a text in inches."""
    t = fig.text(0, 0, s, fontsize=size, **kw)
    w = t.get_window_extent(fig.canvas.get_renderer()).width / fig.dpi
    t.remove()
    return w


def _key_row(fig, x, y, items, size=6.2, gap=0.22):
    """On-figure key from x, y (inches, row center): items (kind, color, label) with kind box, hatch, tri, dens, dot,
    open, num (a homoeolog set number) or frame (an outlined cell); returns the x where it ends."""
    from matplotlib.lines import Line2D
    from matplotlib.patches import Polygon, Rectangle
    W, H = fig.get_size_inches()
    for kind, col, lab in items:
        g = 0.13
        if kind in ("box", "hatch"):
            fig.add_artist(Rectangle((x / W, (y - 0.045) / H), g / W, 0.09 / H, transform=fig.transFigure,
                                     facecolor="white" if kind == "hatch" else col,
                                     edgecolor="#c8c8c8" if kind == "hatch" else "none",
                                     hatch="//////" if kind == "hatch" else None, lw=0))
        elif kind == "tri":
            fig.add_artist(Polygon([((x + 0.01) / W, y / H), ((x + g) / W, (y + 0.05) / H), ((x + g) / W, (y - 0.05) / H)],
                                   closed=True, transform=fig.transFigure, facecolor=col, edgecolor="black", lw=0.3))
        elif kind == "dens":
            for dx_, w_ in ((0.0, 0.05), (0.06, 0.03)):
                fig.add_artist(Rectangle(((x + dx_) / W, (y - 0.05) / H), w_ / W, 0.1 / H, transform=fig.transFigure,
                                         facecolor="#8c8c8c", lw=0))
        elif kind in ("dot", "open"):
            fig.add_artist(Line2D([(x + g / 2) / W], [y / H], marker="o", ms=3.4, lw=0,
                                  mfc=col if kind == "dot" else "white", mec=col, mew=0.7))
        elif kind == "num":
            fig.text((x + g / 2) / W, y / H, "1", ha="center", va="center", fontsize=size)
        elif kind == "frame":
            fig.add_artist(Rectangle(((x + 0.015) / W, (y - 0.045) / H), (g - 0.03) / W, 0.09 / H,
                                     transform=fig.transFigure, fill=False, edgecolor=col, lw=0.6))
        fig.text((x + g + 0.04) / W, y / H, lab, va="center", ha="left", fontsize=size)
        x += g + 0.04 + _text_w(fig, lab, size) + gap
    return x


def _yr_label(v):
    return f"{v / 1e6:g} Myr" if v >= 1e6 else f"{v / 1e3:g} kyr" if v >= 1e3 else f"{v:g} yr"


def _rate(mu, nums=None):
    """Plain-language substitution rate: 7e-9 -> '7 substitutions per site per billion years'; the number as printed
    is recorded in nums."""
    v = f"{mu * 1e9:.3g}"
    if nums is not None:
        v = nums.add("time_axis_rate", v, "report_data.json:mu x 1e9 (substitutions per site per billion years)")
    return f"{v} substitutions per site per billion years"


def _short_name(s, keep=30):
    """Sequence names longer than keep characters shortened in the middle for figure labels (tables keep them)."""
    return s if len(s) <= keep else s[:keep // 2 - 2] + "…" + s[-(keep // 2 - 1):]


def _span(lo, est, hi):
    """'4.32 Myr (2.38–18.5)' when the limits share the estimate's unit, else each limit with its own unit."""
    unit = lambda t: "Myr" if t >= 1e6 else "kyr" if t >= 1e3 else "yr"
    num = lambda t: f"{t / 1e6:.3g}" if t >= 1e6 else f"{t / 1e3:.3g}" if t >= 1e3 else f"{t:.3g}"
    if unit(lo) == unit(est) == unit(hi):
        return f"{num(est)} {unit(est)} ({num(lo)}–{num(hi)})"
    return f"{num(est)} {unit(est)} ({num(lo)} {unit(lo)} – {num(hi)} {unit(hi)})"


def _split_name(t, zof, gsize):
    """A tested split in final-group terms: 'SG1 | SG2' when each side holds whole groups, else 'SG1: 7 | 3'."""
    sides = []
    for part in (t["part_a"], t["part_b"]):
        gs = sorted({zof[n] for n in part if n in zof})
        whole = bool(gs) and all(sum(zof.get(n) == g for n in part) == gsize[g] for g in gs)
        sides.append(("+".join(group_label(g) for g in gs), whole))
    if all(w for _, w in sides) and t["part_b"]:
        return f"{sides[0][0]} | {sides[1][0]}"
    gs = sorted({zof[n] for n in t["part_a"] + t["part_b"] if n in zof})
    return f"{'+'.join(group_label(g) for g in gs)}: {len(t['part_a'])} | {len(t['part_b'])}"


def _split_words(t, zof, gsize):
    """'SG1 vs SG2' for a split between whole groups, 'within SG1: 5 vs 2' for a split inside one."""
    nm = _split_name(t, zof, gsize)
    if ":" in nm:
        g, rest = nm.split(":", 1)
        return f"within {g.strip()}: {rest.strip().replace(' | ', ' vs ')}"
    return " vs ".join(sorted((x.strip() for x in nm.split("|")), key=natural_key))


def _repel(ys, gap):
    """Label positions spread so neighbors are >= gap apart (order kept), staying near their anchors."""
    ys = np.asarray(ys, float)
    o = np.argsort(ys, kind="stable")
    y = ys[o].copy()
    for _ in range(200):
        moved = False
        for k in range(1, y.size):
            if y[k] - y[k - 1] < gap:
                push = (gap - (y[k] - y[k - 1])) / 2
                y[k - 1] -= push
                y[k] += push
                moved = True
        if not moved:
            break
    out = np.empty_like(y)
    out[o] = y
    return out


def _panel_karyotype(fig, ax, P, nums):
    """Every shown sequence a vertical bar to scale, grouped by subgenome, painted by the subgenome its LTR-RT
    lineages point to (opacity: posterior), hatched where copies are too sparse to tell; LTR-RT density in gray to its
    left; called exchanges outlined, with a triangle in the donor's color (the longest >= 1 Mb labeled)."""
    from matplotlib.collections import PatchCollection
    from matplotlib.colors import to_rgba
    from matplotlib.patches import Polygon, Rectangle
    lib, shown, L, z, names = P["lib"], P["shown"], P["L"], P["z"], P["names"]
    xpos, x, prev = {}, 0.0, None
    for i in shown:
        if prev is not None and z[i] != z[prev]:
            x += 0.9
        xpos[i], x, prev = x, x + 1.0, i
    n = len(shown)
    bw = 0.36 if n <= 30 else 0.5
    Lmax = L[shown].max() / 1e6
    num = _homoeolog_numbers(P["ph"])
    ax.set_xlim(-0.75, x - 0.25)
    ax.set_ylim(Lmax * (1.085 if num else 1.03), -Lmax * 0.005)
    ax.spines["bottom"].set_visible(False)
    ax.set_xticks([])
    ax.set_ylabel("position (Mb)")
    pitch = ax.get_position().width * fig.get_size_inches()[0] / (x + 0.5)     # inches per column
    painted = bool(P["tiles"])
    rects, cols, dens = [], [], []
    step = Lmax / 120
    for i in shown:
        s, x0, Lm = names[i], xpos[i], L[i] / 1e6
        if painted:
            ax.add_patch(Rectangle((x0 - bw / 2, 0), bw, Lm, facecolor="white", edgecolor="#c8c8c8", hatch="//////",
                                   lw=0, zorder=1))
        else:
            ax.add_patch(Rectangle((x0 - bw / 2, 0), bw, Lm, facecolor=FIG_LIGHT, edgecolor="none", lw=0, zorder=1))
        if s in P["tiles"]:
            for a_, b_, c_, l_ in zip(*P["tiles"][s]):
                rects.append(Rectangle((x0 - bw / 2, a_ / 1e6), bw, (b_ - a_) / 1e6))
                cols.append(to_rgba(SG_COLORS[c_], FIG_ALPHA[l_]))
        ax.add_patch(Rectangle((x0 - bw / 2, 0), bw, Lm, fill=False, edgecolor=FIG_DARK, lw=0.45, zorder=4))
        if pitch >= 0.2:
            m = lib.seq == P["si"][s]
            e = np.arange(0, Lm + step, step)
            h, _ = np.histogram((lib.start[m] + lib.end[m]) / 2e6, e)
            dens.append((x0, e, h / step))
    if rects:
        ax.add_collection(PatchCollection(rects, facecolors=cols, edgecolors="none", zorder=2,
                                          rasterized=len(rects) > 20000))
    if dens:
        dmax = max(h.max() for _, _, h in dens) or 1.0
        for x0, e, h in dens:
            ax.fill_betweenx(np.repeat(e, 2)[1:-1], x0 - bw / 2 - 0.04 - np.repeat(h, 2) / dmax * 0.3,
                             x0 - bw / 2 - 0.04, color="#8c8c8c", lw=0, zorder=1)
    rot, fs, _ = _name_layout(fig, P)
    for i in shown:
        ax.text(xpos[i], -Lmax * 0.012, _short_name(names[i]), ha="center", va="bottom", rotation=rot, fontsize=fs,
                clip_on=False)
        if i in num:                                                   # homoeolog set number
            ax.text(xpos[i], Lmax * 1.055, nums.add(f"A_set_{names[i]}", num[i],
                                                    "report_data.json:ph.homoeologs.sets (numbered from 1, ordered by "
                                                    "their first member in their lowest subgenome)"),
                    ha="center", va="center", fontsize=fs - 0.4, color=FIG_INK)
    calls = [c for c in P["calls"] if names.index(c["seq"]) in xpos]
    big = sorted((c for c in calls if c["end"] - c["start"] >= 1e6), key=lambda c: c["start"] - c["end"])[:8]
    row = {id(c): k for k, c in enumerate(P["calls"])}
    n_lab = 0
    for c in calls:
        x0, ym = xpos[names.index(c["seq"])], (c["start"] + c["end"]) / 2e6
        h = max((c["end"] - c["start"]) / 1e6, Lmax * 0.006)
        ax.add_patch(Rectangle((x0 - bw / 2, ym - h / 2), bw, h, fill=False, edgecolor="black", lw=0.5, zorder=5))
        ax.add_patch(Polygon([[x0 + bw / 2 + 0.03, ym], [x0 + bw / 2 + 0.2, ym - Lmax * 0.009],
                              [x0 + bw / 2 + 0.2, ym + Lmax * 0.009]], closed=True, facecolor=SG_COLORS[c["origin"]],
                             edgecolor="black", lw=0.3, zorder=5))
        if any(c is b for b in big) and pitch >= 0.3:
            n_lab += 1
            mb = nums.add(f"A_exchange{n_lab}", f"{(c['end'] - c['start']) / 1e6:.1f} Mb",
                          f"report_data.json:calls[{row[id(c)]}] ({c['seq']}), end - start")
            ax.text(x0 + bw / 2 + 0.24, ym, mb, fontsize=5.5, va="center", ha="left", zorder=6)
    return xpos


def _name_layout(fig, P):
    """Sequence names over the karyotype: horizontal when few and short, else vertical; 6 pt (5 pt beyond 30
    sequences) -> (rotation, font size, height in inches)."""
    nm = [_short_name(P["names"][i]) for i in P["shown"]]
    fs = 6 if len(nm) <= 30 else 5
    rot = 0 if len(nm) <= 16 and max(map(len, nm)) <= 8 else 90
    return rot, fs, (0.1 if rot == 0 else max(_text_w(fig, x, fs) for x in nm) + 0.02)


def _group_heads(fig, ax, P, xpos, y_in, nums):
    """Each subgenome's columns under a colored rule and 'SG1 · 7 chromosomes · 1,119 Mb' (the name alone where the
    text would not fit)."""
    from matplotlib.lines import Line2D
    W, H = fig.get_size_inches()
    z, L, shown = P["z"], P["L"], P["shown"]
    bb = ax.get_position()
    xl0, xl1 = ax.get_xlim()
    xf = lambda xd: bb.x0 + (xd - xl0) / (xl1 - xl0) * bb.width
    for g in range(P["N"]):
        xs = [xpos[i] for i in shown if z[i] == g]
        if not xs:
            continue
        col = SG_COLORS[g] if P["N"] > 1 else FIG_DARK
        fig.add_artist(Line2D([xf(min(xs) - 0.4), xf(max(xs) + 0.4)], [y_in / H] * 2, color=col, lw=2,
                              solid_capstyle="butt"))
        n_g = int((z == g).sum())
        mb = f"{L[z == g].sum() / 1e6:,.0f}" if P["N"] > 1 else f"{L.sum() / 1e6:,.0f}"
        lab = (f"{group_label(g)} · {n_g} {P['unit']} · {mb} Mb" if P["N"] > 1 else
               f"one group · {n_g} {P['unit']} · {mb} Mb")
        room = (xf(max(xs) + 0.5) - xf(min(xs) - 0.5)) * W
        if _text_w(fig, lab, 7, fontweight="bold") > room:
            lab = group_label(g) if P["N"] > 1 else "one group"
        else:
            gl = group_label(g) if P["N"] > 1 else "all"
            nums.add(f"A_{gl}_sequences", n_g, f"report_data.json:z (sequences of {gl})")
            nums.add(f"A_{gl}_Mb", mb, f"report_data.json:L (sum over {gl}, Mb)")
        fig.text(xf((min(xs) + max(xs)) / 2), (y_in + 0.05) / H, lab, ha="center", va="bottom", color=col, fontsize=7,
                 fontweight="bold")


def _panel_sharing(fig, ax, P, cmap, norm):
    """Lineage sharing between every pair of shown sequences (observed / expected, all copies), grouped by subgenome
    (strips: subgenome; white lines: group boundaries); names at the right when they fit."""
    from matplotlib.patches import Rectangle
    z, shown = P["z"], P["shown"]
    n = len(shown)
    ax.imshow(np.clip(P["oe"]["all"], norm.vmin, norm.vmax), cmap=cmap, norm=norm, interpolation="nearest")
    s = max(0.6, n / 22)
    for q, i in enumerate(shown):
        col = SG_COLORS[z[i]] if P["N"] > 1 else FIG_GRAY
        ax.add_patch(Rectangle((-0.5 - 1.25 * s, q - 0.5), s, 1, color=col, clip_on=False, lw=0))
        ax.add_patch(Rectangle((q - 0.5, -0.5 - 1.25 * s), 1, s, color=col, clip_on=False, lw=0))
    for b in np.flatnonzero(np.diff(z[shown])) + 0.5:
        ax.axhline(b, color="white", lw=1.0)
        ax.axvline(b, color="white", lw=1.0)
    ax.set_xlim(-0.5 - 1.35 * s, n - 0.5)
    ax.set_ylim(n - 0.5, -0.5 - 1.35 * s)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_xticks([])
    if n <= 30:
        ax.yaxis.tick_right()
        ax.set_yticks(range(n))
        ax.set_yticklabels([_short_name(P["names"][i], 16) for i in shown], fontsize=5.2 if n <= 20 else 4.4)
        ax.tick_params(axis="y", length=0, pad=1.5)
    else:
        ax.set_yticks([])


def _panel_homoeologs(fig, ax, P):
    """Panel B with homoeologs: every shared old insertion (reciprocal closest relatives in two subgenomes) as a dot at its
    positions on the two chromosomes, chromosomes to scale in figure order; homoeologs show as collinear diagonals.
    Two subgenomes: SG1 across, SG2 down; more: every shown chromosome on both axes (same-subgenome blocks shaded,
    empty by construction). Dots on the collinear chain of a called pair in black, others gray; called pairs
    framed (dashed: called one way only)."""
    from matplotlib.patches import Rectangle
    shown, z, L, lib, ho = P["shown"], P["z"], P["L"], P["lib"], P["ph"]["homoeologs"]
    two = P["N"] == 2
    xs = [i for i in shown if z[i] == 0] if two else list(shown)
    ys = [i for i in shown if z[i] == 1] if two else list(shown)

    def offsets(seq):
        o, acc = {}, 0.0
        for i in seq:
            o[i], acc = acc, acc + L[i]
        return o, acc
    (ox, tx), (oy, ty) = offsets(xs), offsets(ys)
    ax.set_xlim(0, tx)
    ax.set_ylim(ty, 0)
    for i in xs:
        ax.add_patch(Rectangle((ox[i], -0.034 * ty), L[i], 0.025 * ty, color=SG_COLORS[z[i]], clip_on=False, lw=0))
        ax.axvline(ox[i], color="#dcdcdc", lw=0.3, zorder=1)
    for i in ys:
        ax.add_patch(Rectangle((-0.034 * tx, oy[i]), 0.025 * tx, L[i], color=SG_COLORS[z[i]], clip_on=False, lw=0))
        ax.axhline(oy[i], color="#dcdcdc", lw=0.3, zorder=1)
    if not two:
        for g in range(P["N"]):
            ig = [i for i in xs if z[i] == g]
            if ig:
                a0, a1 = ox[ig[0]], ox[ig[-1]] + L[ig[-1]]
                ax.add_patch(Rectangle((a0, a0), a1 - a0, a1 - a0, color="#f2f2f2", lw=0, zorder=0))
                if a0 > 0:
                    ax.axvline(a0, color="#8c8c8c", lw=0.6, zorder=1)
                    ax.axhline(a0, color="#8c8c8c", lw=0.6, zorder=1)
    ek, Pp, hb = np.asarray(P["ph"]["ek"]), ho["P"], ho.get("bio")
    mid = (lib.start + lib.end) / 2
    hot = ho["in_chain"] & ho["pair_called"]
    for a_, b_ in ((Pp[:, 0], Pp[:, 1]),) + (((Pp[:, 1], Pp[:, 0]),) if not two else ()):
        X = np.array([ox.get(int(c), np.nan) for c in ek[a_]]) + mid[a_]
        Y = np.array([oy.get(int(c), np.nan) for c in ek[b_]]) + mid[b_]
        ok = np.isfinite(X) & np.isfinite(Y)
        ex = hb["flags"] if hb else np.zeros(len(Pp), bool)
        for m_, col, sz, zo in ((ok & ~hot, "#b9b9b9", 0.5, 2), (ok & hot & ~ex, FIG_ACCENT, 0.8, 3),
                                (ok & hot & ex, FIG_EXCH, 2.2, 5)):
            ax.scatter(X[m_], Y[m_], s=sz, c=col, lw=0, zorder=zo, rasterized=True)
    for (i, h), c in ho["calls"].items():
        j = c["partner"]
        for a_, b_ in ((i, j), (j, i)):
            if j >= 0 and a_ in ox and b_ in oy:
                ax.add_patch(Rectangle((ox[a_], oy[b_]), L[a_], L[b_], fill=False, edgecolor=FIG_ACCENT, lw=0.6,
                                       ls="-" if c["reciprocal"] else (0, (2, 1.5)), zorder=4))
    for t in (hb["transloc"] if hb else []):                       # block collinear with a non-homoeolog
        for a_, b_, xa, xb, ya, yb in ((t["i"], t["j"], t["start_i"], t["end_i"], t["start_j"], t["end_j"]),
                                       (t["j"], t["i"], t["start_j"], t["end_j"], t["start_i"], t["end_i"])):
            if a_ in ox and b_ in oy:
                pad = 0.006 * max(tx, ty)
                ax.add_patch(Rectangle((ox[a_] + xa - pad, oy[b_] + ya - pad), xb - xa + 2 * pad, yb - ya + 2 * pad,
                                       fill=False, edgecolor=FIG_EXCH, lw=0.6, ls=(0, (1, 1)), zorder=6))
    for sp in ax.spines.values():
        sp.set_visible(False)
    fs = 5.2 if max(len(xs), len(ys)) <= 12 else 4.4
    ax.set_xticks([])
    ax.set_yticks([])
    if len(ys) <= 30:
        ax.yaxis.tick_right()
        ax.set_yticks([oy[i] + L[i] / 2 for i in ys])
        ax.set_yticklabels([_short_name(P["names"][i], 12) for i in ys], fontsize=fs)
        ax.tick_params(axis="y", length=0, pad=1.5)
    if two and len(xs) <= 30:
        ax.set_xticks([ox[i] + L[i] / 2 for i in xs])
        ax.set_xticklabels([_short_name(P["names"][i], 12) for i in xs], fontsize=fs, rotation=90)
        ax.tick_params(axis="x", length=0, pad=1.5)


def _oe_bar(fig, x0, y0, w, norm, cmap, label):
    from matplotlib.cm import ScalarMappable
    from matplotlib.ticker import NullFormatter, NullLocator
    cax = _ax(fig, x0, y0, w, 0.07)
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal",
                      ticks=[0.25, 0.5, 1, 2, 4])
    cb.ax.set_xticklabels(["¼×", "½×", "1×", "2×", "4×"])
    cb.ax.xaxis.set_minor_locator(NullLocator())
    cb.ax.xaxis.set_minor_formatter(NullFormatter())
    cb.ax.tick_params(length=2, pad=1, labelsize=6)
    cb.outline.set_linewidth(0.4)
    if label:
        fig.text((x0 + w / 2) / fig.get_size_inches()[0], (y0 + 0.1) / fig.get_size_inches()[1], label, ha="center",
                 va="bottom", fontsize=5.6, color="#444444")


def _pca(P):
    """Phased sequences placed by their lineage-sharing profiles: PCA of the rows of the log observed/expected matrix
    the partition maximizes (masked links; the missing self-entry takes its column's mean). Each kept half-library's
    profiles are projected onto the same axes -> (positions [C, 2], half-library positions [h, C, 2] or None,
    variance share of the two axes)."""
    ph = P["ph"]

    def M_of(Wm):
        Lm = log_oe(np.asarray(Wm, float))
        np.fill_diagonal(Lm, np.nan)
        return np.where(np.isnan(Lm), np.nanmean(Lm, 0)[None], Lm)
    M = M_of(ph["W"])
    mu_ = M.mean(0)
    _, S, Vt = np.linalg.svd(M - mu_, full_matrices=False)
    k = min(2, Vt.shape[0])
    V = np.zeros((M.shape[1], 2))
    V[:, :k] = Vt[:k].T
    X = (M - mu_) @ V
    z = np.asarray(ph["z"])
    for j in range(2):                                 # the first subgenome on the left; axis 2 toward its largest load
        if (X[z == 0, j].mean() > 0) if j == 0 else (X[np.argmax(np.abs(X[:, j])), j] < 0):
            V[:, j] *= -1
            X[:, j] *= -1
    reps = np.array([(M_of(h) - mu_) @ V for h in ph["halves"]]) if ph.get("halves") else None
    share = np.zeros(2)
    share[:k] = S[:k] ** 2 / max(float((S ** 2).sum()), 1e-300)
    return X, reps, share


def _panel_pca(fig, ax, P, nums):
    """Each phased sequence on the two leading axes of its lineage-sharing profile (dot) and in each half-library
    (cloud), colored by subgenome; names when they fit; the share of sequences kept in their group by >= 90% of the
    half-libraries."""
    z, names, ph = P["z"], P["names"], P["ph"]
    X, reps, share = _pca(P)
    col = lambda i: SG_COLORS[z[i]] if P["N"] > 1 else FIG_DARK
    if reps is not None:
        for i in range(len(names)):
            ax.scatter(reps[:, i, 0], reps[:, i, 1], s=1.5, color=col(i), alpha=0.25, lw=0, rasterized=True, zorder=2)
    for i in range(len(names)):
        ax.scatter([X[i, 0]], [X[i, 1]], s=14, color=col(i), edgecolors="white", linewidths=0.5, zorder=4)
    ax.axhline(0, color="#e2e2e2", lw=0.5, zorder=0)
    ax.axvline(0, color="#e2e2e2", lw=0.5, zorder=0)
    pca = "report_data.json:ph.W (PCA of the log observed/expected rows), variance share of axis"
    v1 = nums.add("C_axis1_variance", f"{100 * share[0]:.0f}%", f"{pca} 1")
    v2 = nums.add("C_axis2_variance", f"{100 * share[1]:.0f}%", f"{pca} 2")
    ax.set_xlabel(f"lineage-sharing axis 1 ({v1} of variance)", fontsize=6.2)
    ax.set_ylabel(f"axis 2 ({v2})", fontsize=6.2)
    ax.tick_params(labelsize=5.6)
    ax_ = np.r_[X[:, 0], reps[:, :, 0].ravel()] if reps is not None else X[:, 0]
    ay_ = np.r_[X[:, 1], reps[:, :, 1].ravel()] if reps is not None else X[:, 1]
    xr, yr = (np.ptp(ax_) or 1.0), (np.ptp(ay_) or 1.0)
    ax.set_xlim(ax_.min() - 0.32 * xr, ax_.max() + 0.32 * xr)
    ax.set_ylim(ay_.min() - 0.24 * yr, ay_.max() + 0.25 * yr)
    if P["N"] > 1:
        gs = [g for g in range(P["N"]) if (z == g).any()]
        top = [((reps[:, z == g, 1].max() if reps is not None else X[z == g, 1].max()) + 0.06 * yr) for g in gs]
        per_in = np.diff(ax.get_xlim())[0] / (ax.get_position().width * fig.get_size_inches()[0])
        xs = (_repel([X[z == g, 0].mean() for g in gs], 0.3 * per_in) if P["N"] > 2 else
              [X[z == g, 0].mean() for g in gs])
        if P["N"] == 2:
            top = [ay_.max() + 0.1 * yr] * len(gs)
        for g, x_, y_ in zip(gs, xs, top):
            ax.text(x_, y_, group_label(g), color=SG_COLORS[g], fontsize=6.5, fontweight="bold", ha="center",
                    va="bottom")
    if len(names) <= 24:
        for g in range(max(P["N"], 1)):
            idx = np.flatnonzero(z == g)
            if not idx.size:
                continue
            side = 1 if X[idx, 0].mean() >= np.median(X[:, 0]) else -1
            ys = _repel(X[idx, 1], 0.085 * yr)
            x0 = (X[idx, 0].max() + 0.09 * xr) if side > 0 else (X[idx, 0].min() - 0.09 * xr)
            for i, yy in zip(idx, ys):
                ax.plot([X[i, 0], x0], [X[i, 1], yy], color="#c4c4c4", lw=0.3, zorder=1)
                ax.text(x0 + side * 0.01 * xr, yy, _short_name(names[i], 14), fontsize=5.0, va="center",
                        ha="left" if side > 0 else "right", color="#555555")
    sup = np.asarray(ph.get("support", []), float)
    ok = np.isfinite(sup)
    n_h = len(reps) if reps is not None else 0
    if n_h:
        note = f"clouds: {nums.add('C_half_libraries', n_h, 'report_data.json:ph.halves (count)')} half-libraries"
    else:
        most = nums.add("C_half_libraries_max", FIG_HALVES_MAX, "base.constants:FIG_HALVES_MAX")
        note = f"half-libraries not kept (> {most} sequences)"
    if ok.any() and P["N"] > 1:
        kept = nums.add("C_kept_group", int((sup[ok] >= SUPPORT_OK).sum()),
                        "report_data.json:ph.support (>= SUPPORT_OK)")
        of = nums.add("C_kept_group_of", int(ok.sum()), "report_data.json:ph.support (finite)")
        level = nums.add("C_kept_group_level", f"{100 * SUPPORT_OK:.0f}%", "base.constants:SUPPORT_OK")
        note += f"\n{kept} of {of} {P['unit']} kept their group in ≥ {level}"
    ax.text(0.5, 0.03, note, transform=ax.transAxes, fontsize=5.3, ha="center", va="bottom", color="#444444")


def _panel_tests(fig, axr, axs, P, nums, max_rows=6):
    """The two tests every candidate split must pass, one row per split the search tried (kept first): reproduced by
    disjoint half-libraries (adjusted Rand >= MIN_REP, left) and at least MIN_REL as strong as the split above it
    (right), on the copies the decision used (all copies, or the age class that kept it)."""
    ph, z = P["ph"], P["z"]
    at = {id(t): i for i, t in enumerate(ph.get("tree", []))}               # each split's place in ph.tree
    T = [t for t in ph.get("tree", []) if t.get("replicability_all") is not None or t["accepted"] or t["forced"]]
    T = sorted(T, key=lambda t: (not (t["accepted"] or t["forced"]), len(t["level"]), t["level"]))[:max_rows]
    for ax_ in (axr, axs):
        ax_.set_yticks([])
        ax_.spines["left"].set_visible(False)
    if not T:
        for ax_ in (axr, axs):
            ax_.axis("off")
        msg = ("groups given by the user (--assign): no split tested" if ph.get("source") == "user" else
               "no split could be tested")
        axr.text(0.0, 0.5, msg, transform=axr.transAxes, fontsize=6.3, ha="left", va="center", color="#666666")
        return
    zof = {nm: int(z[i]) for i, nm in enumerate(P["names"])}
    gsize = np.bincount(z, minlength=P["N"])
    reason = {"not reproduced": "not reproduced", "weaker than half its parent": "too weak",
              "present at every age": "same at every age", "no excess": "no excess",
              "pairs: not reproduced": "not reproduced"}
    for k, t in enumerate(T):
        kept = bool(t["accepted"] or t["forced"])
        c = FIG_INK if kept else "#9a9a9a"
        s = t["copies"] if (kept and t.get("copies")) else "all"
        src, row = f"report_data.json:ph.tree[{at[id(t)]}]", f"D_split{k + 1}"
        lab = (f"all {P['unit']}: {len(t['part_a'])} vs {len(t['part_b'])}" if P["N"] == 1 else
               _split_words(t, zof, gsize)) + ("" if s == "all" else f" ({s} copies)")
        nums.add(row, lab, f"{src}.part_a, part_b (sizes and groups), copies")
        r = t.get(f"replicability_{s}")
        st, pst = t.get(f"strength_{s}"), t.get(f"parent_strength_{s}")
        if r is not None:
            v = max(r, -0.05)
            axr.plot([v], [k], "o", ms=4 if kept else 3.2, color=c, clip_on=False, zorder=3)
            axr.text(v - 0.06 if v > 0.45 else v + 0.07, k,
                     nums.add(f"{row}_replicability", f"{r:.2f}", f"{src}.replicability_{s}"), fontsize=5.2,
                     va="center", ha="right" if v > 0.45 else "left", color=c)
        if t["rule"] in ("pairs", "pairs: not reproduced"):
            axs.text(0.6, k, "homologue pairs\nkept apart", fontsize=5.0, va="center", ha="center", color=c)
        elif st is not None and pst:
            rel = st / pst if pst > 0 else np.nan
            if np.isfinite(rel):
                v = min(max(rel, 0.0), 1.2)
                axs.plot([v], [k], "o", ms=4 if kept else 3.2, color=c, clip_on=False, zorder=3)
                axs.text(v + 0.07 if v < 0.8 else v - 0.07, k,
                         nums.add(f"{row}_strength", f"{100 * rel:.0f}%", f"{src}.strength_{s} / parent_strength_{s}"),
                         fontsize=5.2, va="center", ha="left" if v < 0.8 else "right", color=c)
        else:
            axs.text(0.6, k, "first split", fontsize=5.4, va="center", ha="center", color=c)
        axr.text(-0.1, k, lab, transform=axr.get_yaxis_transform(), fontsize=5.8, va="center", ha="right", color=c,
                 fontweight="bold" if kept else "normal")
        word = ("forced (-n)" if t["forced"] else "kept") if kept else reason.get(
            t["rule"], "no homoeologs\nacross it"
            if str(t["rule"]).endswith("without homoeologs across it") else "rejected")
        nums.add(f"{row}_result", " ".join(word.split("\n")), f"{src}.accepted, forced, rule")
        axs.text(1.04, k, word, transform=axs.get_yaxis_transform(), fontsize=5.8, va="center", ha="left", color=c,
                 fontweight="bold" if kept else "normal")
    for ax_, thr, lim, ttl in ((axr, MIN_REP, (-0.05, 1.05), "reproduced by\nindependent halves"),
                               (axs, MIN_REL, (0.0, 1.2), "as strong as the\nsplit above")):
        ax_.axvspan(thr, lim[1], color="#ececec", lw=0, zorder=0)
        ax_.axvline(thr, color="#888888", lw=0.5, ls=(0, (2, 2)), zorder=1)
        ax_.set_xlim(*lim)
        ax_.set_ylim(len(T) - 0.5, -0.6)
        ax_.tick_params(axis="x", labelsize=5.5, length=2)
        ax_.set_title(ttl, fontsize=5.8, pad=2, color="#444444")
    axr.set_xticks([0, MIN_REP, 1])
    axr.set_xticklabels(["0", f"{MIN_REP:g}", "1"])
    axs.set_xticks([0, MIN_REL, 1])
    axs.set_xticklabels(["0", "½" if MIN_REL == 0.5 else f"{MIN_REL:g}", "1×"])


def _tests_header(P):
    ph, N = P["ph"], P["N"]
    if ph.get("source") == "user":
        return f"{FIG_NUM.get(N, N).capitalize()} groups given by the user (--assign)"
    if any(t["forced"] for t in ph.get("tree", [])):
        return f"{FIG_NUM.get(N, N).capitalize()} groups set with -n: the splits that fail a test were forced"
    if N == 1:
        return "Why one group: no split passes both tests"
    kept = sum(bool(t["accepted"]) for t in ph.get("tree", []))
    if N == 2 and kept == 1:
        return "Why two subgenomes: only the first split passes both tests"
    return f"Why {FIG_NUM.get(N, N)} subgenomes: the kept splits pass both tests"


def _time_axis(ax, P, labels=True, xlabel=True, nums=None):
    """Years ago on a log axis, oldest left, present right; the rate behind it in the axis label (recorded in nums)."""
    from matplotlib.ticker import FixedFormatter, FixedLocator
    t0, t1 = P["trange"]
    ax.set_xscale("log")
    ax.set_xlim(t1, t0)
    tk = [10.0 ** e for e in range(0, 12) if t0 <= 10.0 ** e <= t1]
    lab = [_yr_label(v) for v in tk]
    if len(tk) > 5:
        lab = [x if (len(tk) - 1 - j) % 2 == 0 else "" for j, x in enumerate(lab)]
    ax.xaxis.set_major_locator(FixedLocator(tk))
    ax.xaxis.set_minor_locator(FixedLocator([m * v for v in [10.0 ** e for e in range(0, 12)] for m in range(2, 10)
                                             if t0 <= m * v <= t1]))
    ax.xaxis.set_major_formatter(FixedFormatter(lab if labels else [""] * len(tk)))
    ax.tick_params(axis="x", which="minor", length=1.5, width=0.4)
    if xlabel:
        ax.set_xlabel(f"years ago ({_rate(P['mu'], nums)}" + ("" if P["mu_given"] else "; default rate") + ")",
                      fontsize=6.2)


def _first_pair_events(P):
    """Reported split and merger of the first subgenome pair, in years: (split dict, merger dict) or Nones."""
    dt = P["dt"]
    if not dt or not dt["pairs"] or P["N"] != 2:
        return None, None
    r = dt["pairs"][0]
    yr = lambda k: k / (2 * P["mu"])
    sp = r["split"].get(r.get("split_reported")) if r.get("split_reported") else None
    mg = r["merger"] if r["merger"] and r["merger"]["datable"] else None
    s = dict(lo=yr(sp["lo"]), est=yr(sp["estimate"]), hi=yr(sp["hi"])) if sp else None
    m = dict(lo=P["trange"][0] if mg["bound"] else yr(mg["lo"]), hi=yr(mg["hi"]), bound=mg["bound"]) if mg else None
    return s, m


def _bands(ax, P):
    """The reported split and merger intervals (two subgenomes) as background bands of a time axis."""
    s, m = _first_pair_events(P)
    if s:
        ax.axvspan(s["lo"], s["hi"], color=FIG_SPLIT, lw=0, zorder=0)
    if m:
        ax.axvspan(m["lo"], m["hi"], color=FIG_MERGE, lw=0, zorder=0)


def _regime_strip(fig, ax, P, h=0.27, pad=0.02):
    """A strip above a time axis naming its periods between the bands (two subgenomes): one ancestral pool, the split,
    the progenitors apart, the merger, one nucleus; words shortened or dropped where a period is too narrow."""
    s, m = _first_pair_events(P)
    if not s and not m:
        return
    W, H = fig.get_size_inches()
    bb = ax.get_position()
    rax = fig.add_axes([bb.x0, bb.y1 + pad / H, bb.width, h / H])
    rax.set_xscale("log")
    rax.set_xlim(ax.get_xlim())
    rax.set_ylim(0, 1)
    rax.axis("off")
    t0, t1 = P["trange"]
    segs = []
    if s:
        segs += [(s["hi"], t1, "white", ["before the split:\none ancestral pool", "before\nsplit"]),
                 (s["lo"], s["hi"], FIG_SPLIT, ["split"]),
                 (m["hi"] if m else t0, s["lo"], "white", ["progenitors apart:\neach amplifies its own lineages",
                                                         "progenitors\napart"])]
    if m:
        segs.append((m["lo"], m["hi"], FIG_MERGE, ["merger\n(bound)", "merger"] if m["bound"] else ["merger"]))
        if not m["bound"]:
            segs.append((t0, m["lo"], "white", ["one nucleus", "after"]))
    span = math.log(t1 / t0)
    for a_, b_, col, words in segs:
        a_, b_ = max(a_, t0), min(b_, t1)
        if b_ <= a_:
            continue
        rax.axvspan(a_, b_, ymin=0, ymax=0.72, color=col, lw=0)
        room = math.log(b_ / a_) / span * bb.width * W
        for txt in words:
            if max(_text_w(fig, x_, 5.8) for x_ in txt.split("\n")) <= room - 0.04:
                rax.text(math.sqrt(a_ * b_), 0.36, txt, ha="center", va="center", fontsize=5.8, color="#333333")
                break
    rax.plot([t0, t1], [0, 0], color="#bbbbbb", lw=0.4)


def _panel_births(ax, axm, P, nums):
    """Transposon births per subgenome: deconvolved copies per log-age bin (one ridge each, on one scale); the margin
    gives each subgenome's copies."""
    gs = [g for g in range(max(P["N"], 1)) if (g, "all") in P["hist"]]
    axm.axis("off")
    if not gs:
        ax.text(0.5, 0.5, "too few copies with LTR ages", transform=ax.transAxes, ha="center", va="center",
                fontsize=6.3, color="#666666")
        return
    vmax = max(P["hist"][(g, "all")]["y"].max() for g in gs) or 1.0
    off = 0.62
    for k, g in enumerate(gs):
        hh = P["hist"][(g, "all")]
        base = (len(gs) - 1 - k) * off
        col = SG_COLORS[g] if P["N"] > 1 else FIG_DARK
        ax.fill_between(hh["t"], base, base + hh["y"] / vmax, color=col, alpha=0.35, lw=0, zorder=2)
        ax.plot(hh["t"], base + hh["y"] / vmax, color=col, lw=0.9, zorder=3)
        tot, big = hh["n"], max(P["hist"][(g_, "all")]["n"] for g_ in gs)
        axm.text(0.0, base + 0.3, group_label(g) if P["N"] > 1 else "all", fontsize=6, fontweight="bold",
                 va="center", ha="left", color=col)
        axm.barh(base + 0.3, 0.85 * tot / big, left=0.62, height=0.22, color=col, lw=0)
        gl = group_label(g) if P["N"] > 1 else "all"
        n_ = nums.add(f"E_{gl}_copies", f"{tot:,.0f}", f"report_data.json:hist[{g}, all].n (copies with an LTR age)")
        axm.text(0.62 + 0.85 * tot / big + 0.05, base + 0.3, n_, fontsize=5.2, va="center", ha="left", color=col)
    top = off * (len(gs) - 1) + 1.08
    for a_ in (ax, axm):
        a_.set_ylim(0, top)
    axm.set_xlim(0, 2.45)
    axm.text(0.62, top, "copies", fontsize=5.6, va="top", ha="left", color="#555555")
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)


def _panel_clades(ax, axm, P, nums, k=6):
    """One row per main annotated LTR-RT clade (most copies; 'unknown' left out): its copies' births in each subgenome
    (deconvolved, each row scaled to its peak); the margin gives the clade's copies per subgenome."""
    tot = defaultdict(float)
    for r in P["bio"].get("clades", []):
        if str(r["clade"]).lower() not in ("unknown", "na", ""):
            tot[str(r["clade"])] += r["n"]
    gs = range(max(P["N"], 1))
    clades = [c for c, _ in sorted(tot.items(), key=lambda kv: -kv[1]) if any((g, c) in P["hist"] for g in gs)][:k]
    axm.axis("off")
    if not clades:
        ax.text(0.5, 0.5, f"no annotated clade with ≥ {nums.add('F_min_copies', MIN_HIST, 'base.constants:MIN_HIST')} "
                "copies in a subgenome", transform=ax.transAxes, ha="center", va="center", fontsize=6.3,
                color="#666666")
        ax.set_yticks([])
        return
    cnt = {(int(r["subgenome"]), str(r["clade"])): r["n"] for r in P["bio"].get("clades", [])}
    nmax = max([cnt.get((g, c), 0) for c in clades for g in gs] + [1])
    hb = 0.7 / max(P["N"], 1)
    for j, c in enumerate(clades):
        base = (len(clades) - 1 - j) * 1.0
        vmax = max((P["hist"][(g, c)]["y"].max() for g in gs if (g, c) in P["hist"]), default=1.0) or 1.0
        for g in gs:
            col = SG_COLORS[g] if P["N"] > 1 else FIG_DARK
            if (g, c) in P["hist"]:
                hh = P["hist"][(g, c)]
                ax.fill_between(hh["t"], base, base + hh["y"] / vmax * 0.86, color=col, alpha=0.2, lw=0)
                ax.plot(hh["t"], base + hh["y"] / vmax * 0.86, color=col, lw=0.8)
            n_ = cnt.get((g, c), 0)
            if P["N"] <= 2:
                yy = base + 0.08 + hb * (max(P["N"], 1) - 1 - g) + hb / 2
                axm.barh(yy, n_ / nmax, height=hb * 0.82, color=col, lw=0)
                gl = group_label(g) if P["N"] > 1 else "all"
                axm.text(n_ / nmax + 0.04, yy, nums.add(f"F_{c}_{gl}_copies", f"{n_:,}",
                                                         f"report_data.json:bio.clades (subgenome {g}, clade {c}).n"),
                         fontsize=5.1, va="center", ha="left", color=col)
        if P["N"] > 2:                                 # one stacked bar (subgenome shares) and the clade's total
            ns = [cnt.get((g, c), 0) for g in gs]
            tot_c = sum(ns)
            left = 0.0
            for g, n_ in zip(gs, ns):
                axm.barh(base + 0.4, n_ / tot_c, left=left, height=0.38, color=SG_COLORS[g], lw=0)
                left += n_ / tot_c
            axm.text(1.04, base + 0.4, nums.add(f"F_{c}_copies", f"{tot_c:,}",
                                                f"report_data.json:bio.clades (clade {c}).n, summed over subgenomes"),
                     fontsize=5.1, va="center", ha="left", color=FIG_DARK)
        ax.text(P["trange"][1] / 1.06, base + 0.12, c, fontsize=6, ha="left", va="bottom", fontstyle="italic")
    for a_ in (ax, axm):
        a_.set_ylim(-0.03, len(clades) + 0.02)
    axm.set_xlim(0, 1.9)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)


def _triangle(ax, OE, xs, norm, cmap):
    """Lineage sharing of every pair of columns as a rhombus below both (a stained-glass triangle)."""
    from matplotlib.collections import PolyCollection
    polys, vals = [], []
    for i in range(len(xs)):
        for j in range(i + 1, len(xs)):
            X, Y = (xs[i] + xs[j]) / 2, -(xs[j] - xs[i]) / 2
            polys.append([(X - 0.5, Y), (X, Y + 0.5), (X + 0.5, Y), (X, Y - 0.5)])
            vals.append(OE[i, j] if np.isfinite(OE[i, j]) else OE[j, i])
    ax.add_collection(PolyCollection(polys, array=np.ma.masked_invalid(np.clip(np.array(vals, float), norm.vmin,
                                                                               norm.vmax)),
                                     cmap=cmap, norm=norm, edgecolors="face", linewidths=0.15))


def _panel_age_classes(fig, ax, P, cmap, norm):
    """B's lineage sharing within young, middle-aged and old copies (LTR-K2P terciles), each triangle centered over
    its copies' age range (bracket beneath; leader lines where triangles had to move apart)."""
    from matplotlib.lines import Line2D
    edges = P["ph"].get("strata_edges") or []
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    if len(edges) != 2 or not all(s in P["oe"] for s in AGE_STRATA):
        ax.text(0.5, 0.5, "no copy-age classes", transform=ax.transAxes, ha="center", va="center", fontsize=6.3,
                color="#666666")
        return
    W, H = fig.get_size_inches()
    yr = lambda k: k / (2 * P["mu"])
    t0, t1 = P["trange"]
    rng = [("young", t0, yr(edges[0])), ("middle", yr(edges[0]), yr(edges[1])), ("old", yr(edges[1]), t1)]
    bb = ax.get_position()
    hin, win = bb.height * H, bb.width * W
    span = math.log(t1 / t0)
    n = len(P["shown"])
    th = hin * 0.6
    tw = min(2 * th * 1.15, win / 3.6)
    want = [(bb.x0 + math.log(t1 / math.sqrt(max(a_, t0) * min(b_, t1))) / span * bb.width) * W for _, a_, b_ in rng]
    lo_, hi_ = bb.x0 * W + tw / 2, (bb.x0 + bb.width) * W - tw / 2
    xs_ = np.clip(_repel(np.clip(want, lo_, hi_), tw + 0.08), lo_, hi_)
    tr = ax.get_xaxis_transform()
    z = P["z"]
    for (s_, a_, b_), x_in, xw in zip(rng, xs_, want):
        sub = fig.add_axes([(x_in - tw / 2) / W, bb.y0 + (hin * 0.38) / H, tw / W, th / H])
        _triangle(sub, P["oe"][s_], list(range(n)), norm, cmap)
        sub.set_xlim(-0.6, n - 0.4)
        sub.set_ylim(-(n - 1) / 2 - 0.55, 0.9)
        sub.axis("off")
        for g in range(P["N"]):
            gi = [k for k, i in enumerate(P["shown"]) if z[i] == g]
            if gi and P["N"] > 1:
                sub.plot([min(gi) - 0.45, max(gi) + 0.45], [0.55, 0.55], color=SG_COLORS[g], lw=1.4,
                         solid_capstyle="butt")
        if abs(x_in - xw) > 0.05:
            fig.add_artist(Line2D([x_in / W, xw / W], [(bb.y0 * H + hin * 0.38) / H, (bb.y0 * H + hin * 0.33) / H],
                                  color="#999999", lw=0.5))
        a2, b2 = max(a_, t0) * 1.05, min(b_, t1) / 1.05
        ax.plot([a2, b2], [0.3, 0.3], color="#777777", lw=0.6, transform=tr, solid_capstyle="butt")
        for e_ in (a2, b2):
            ax.plot([e_, e_], [0.25, 0.35], color="#777777", lw=0.6, transform=tr)
        fig.text(x_in / W, bb.y0 + hin * 0.21 / H, {"young": "young copies", "middle": "middle-aged",
                                                     "old": "old copies"}[s_],
                 ha="center", va="top", fontsize=5.6, color="#444444")


def _panel_profile(ax, P):
    """Share of closest-relative copy pairs that sit in two different subgenomes, by their divergence time: points
    (area: pairs) with 95% intervals, one series per subgenome pair (line in one subgenome's color, dots in the
    other's when there are several pairs)."""
    for pr in P["prof"]:
        x, y = pr["pair"]
        lc, fc = (FIG_INK, FIG_INK) if P["N"] == 2 else (SG_COLORS[x], SG_COLORS[y])
        ax.vlines(pr["t"], 100 * pr["lo"], 100 * pr["hi"], color=lc, lw=0.5, zorder=3)
        ax.plot(pr["t"], 100 * pr["f"], color=lc, lw=0.7, zorder=3)
        ax.scatter(pr["t"], 100 * pr["f"], s=2 + 10 * np.sqrt(pr["n"] / pr["n"].max()), color=fc, zorder=4,
                   edgecolors="white", linewidths=0.3)
    if P["N"] > 2 and P["prof"]:
        ys = _repel([100 * pr["f"][0] for pr in P["prof"]], 0.09 * max(100 * pr["hi"].max() for pr in P["prof"]))
        for pr, yy in zip(P["prof"], ys):
            x, y = pr["pair"]
            ax.text(pr["t"][0] / 1.25, yy, f"{group_label(x)}–{group_label(y)}", fontsize=5.2, va="center",
                    ha="left", color=SG_COLORS[x])
    ax.set_ylabel("closest relative in\nanother subgenome (%)", fontsize=6)
    top = max((100 * pr["hi"].max() for pr in P["prof"]), default=50.0)
    ax.set_ylim(0, max(5.0, 1.12 * top))


def _estimate_rows(P):
    """Rows of the estimate panel: (kind, estimator, estimate dict in K2P units or None, reported, pair prefix, the
    estimate's report_data.json key). Two subgenomes: every estimator; more: each pair's reported split and merger."""
    rows, pairs = [], list(enumerate(P["dt"]["pairs"]))
    if len(pairs) > 6:                                 # the six pairs with the oldest reported splits (all: dating.tsv)
        old = lambda r: (r["split"].get(r.get("split_reported")) or {}).get("estimate", -1.0) if r.get("split_reported") else -1.0
        pairs = sorted(pairs, key=lambda pr: old(pr[1]), reverse=True)[:6]
    for p, r in pairs:
        pre = f"{group_label(r['groups'][0])}–{group_label(r['groups'][1])} " if P["N"] > 2 else ""
        mg = r["merger"]
        if P["N"] == 2:
            for k in ("inherited", "consensus", "orthologs", "confinement", "profile"):
                rows.append(("split", k, r["split"].get(k), r.get("split_reported") == k, pre,
                             f"dt.pairs[{p}].split.{k}"))
            for k, m_, own in (("ltr_clock", r.get("merger_ltr"), "merger_ltr"),
                               ("genealogy_clock", r.get("merger_genealogy"), "merger_genealogy")):
                rep = bool(mg) and mg.get("source") == k
                rows.append(("merger", k, mg if rep else m_, rep, pre, f"dt.pairs[{p}].{'merger' if rep else own}"))
        else:
            k = r.get("split_reported")
            rows.append(("split", k, r["split"].get(k) if k else None, True, pre,
                         f"dt.pairs[{p}].split.{k}" if k else f"dt.pairs[{p}].split_reported"))
            rows.append(("merger", mg.get("source") if mg else None, mg, True, pre, f"dt.pairs[{p}].merger"))
    return rows


def _panel_estimates(fig, ax, axm, P, rows, nums, fs=5.4):
    """Every row on the time axis (an interval with its estimate; a bound runs to the present): its value beside the
    interval on the free side, the estimator (dating.tsv's name) in the margin; bold: reported, gray: not."""
    t0, t1 = P["trange"]
    yr = lambda k: k / (2 * P["mu"])
    r_ = fig.canvas.get_renderer()
    bb = ax.get_window_extent(r_)
    axm.axis("off")
    for y, (kind, key, e, rep, pre, src) in enumerate(rows):
        c, wt = (FIG_INK, "bold") if rep else ("#8f8f8f", "normal")
        name = f"{pre}{kind}" + (f" · {FIG_EST[key]}" if key and not pre else "")
        axm.text(0.0, y, name, fontsize=5.2, va="center", ha="left", color=c, fontweight=wt)
        row = f"H_{pre.strip() + '_' if pre else ''}{kind}_{key or 'none'}"
        if not e or (kind == "merger" and not e["datable"]):
            nums.add(row, "not datable", f"report_data.json:{src}" + (".datable" if e else " (none)"))
            ax.text(0.5, y, "not datable", transform=ax.get_yaxis_transform(), fontsize=fs, va="center",
                    ha="center", color="#a8a8a8", fontstyle="italic")
            continue
        bound = kind == "merger" and e["bound"]
        lo, hi = (t0 if bound else yr(e["lo"])), yr(e["hi"])
        ax.plot([lo, hi], [y, y], color=c, lw=1.5 if rep else 0.7, solid_capstyle="butt", zorder=3)
        if bound:
            ax.plot([hi], [y], marker="|", ms=5.5, mew=1.1 if rep else 0.7, color=c, zorder=4)
            txt = nums.add(row, f"{kind} < {_ma(hi)}", f"report_data.json:{src}.hi / (2 mu), years")
        else:
            ax.plot([yr(e["estimate"])], [y], "o", ms=3.8 if rep else 2.4, color=c, zorder=4,
                    mfc=c if kind == "split" else "white", mec=c, mew=0.9)
            txt = nums.add(row, f"{kind} {_span(lo, yr(e['estimate']), hi)}",
                           f"report_data.json:{src}.estimate, lo, hi / (2 mu), years")
        for side in (("right", "left") if kind == "split" else ("left", "right")):
            t = ax.text(lo / 1.12 if side == "right" else hi * 1.12, y, txt, fontsize=fs, va="center",
                        ha="left" if side == "right" else "right", color=c, fontweight=wt)
            ex = t.get_window_extent(r_)
            if ex.x0 >= bb.x0 - 1 and ex.x1 <= bb.x1 + 1:
                break
            t.remove()
        else:
            ax.text(0.99, y, txt, transform=ax.get_yaxis_transform(), fontsize=fs, va="center", ha="right", color=c,
                    fontweight=wt)
    for a_ in (ax, axm):
        a_.set_ylim(len(rows) - 0.4, -0.6)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)


def _panel_history(fig, P, x0, wt, wm, y_top, cmap, norm, nums):
    """E-H in one column whose time axes all span [x0, x0 + wt] (inches); the margin of width wm at their right
    carries copy counts (E, F) and the estimators' names (H)."""
    dt = P["dt"]
    rows = _estimate_rows(P) if dt and dt["pairs"] else []
    hH = 0.8 if (dt and P["prof"]) else 0.5
    hF = 1.2 if len(rows) <= 6 else max(0.8, 1.2 - 0.09 * (len(rows) - 6))
    hR = 0.11 * len(rows) + 0.1
    gap, strip = 0.33, 0.27 if P["N"] == 2 and dt else 0.0
    y = y_top - strip - 0.04 - 0.62
    axE, axEm = _ax(fig, x0, y, wt, 0.62), _ax(fig, x0 + wt + 0.08, y, wm - 0.08, 0.62)
    _time_axis(axE, P, labels=False, xlabel=False)
    _bands(axE, P)
    _panel_births(axE, axEm, P, nums)
    _regime_strip(fig, axE, P, h=strip or 0.27)
    _head(fig, axE, "E", "When each subgenome's LTR retrotransposons were born" if P["N"] > 1 else
          "When the LTR retrotransposons were born", FIG_XR, dy=strip + 0.07)
    y -= gap + hF
    axF, axFm = _ax(fig, x0, y, wt, hF), _ax(fig, x0 + wt + 0.08, y, wm - 0.08, hF)
    _time_axis(axF, P, labels=False, xlabel=False)
    _bands(axF, P)
    _panel_clades(axF, axFm, P, nums)
    _head(fig, axF, "F", "Births of each main LTR-RT clade" + (", per subgenome" if P["N"] > 1 else ""), FIG_XR)
    y -= gap + 0.6
    axG = _ax(fig, x0, y, wt, 0.6)
    _time_axis(axG, P, labels=False, xlabel=False)
    _bands(axG, P)
    _panel_age_classes(fig, axG, P, cmap, norm)
    if P["dotplot"]:                                                    # B shows homoeologs; the key sits here
        _head(fig, axG, "G", "Lineages shared by each pair of " + P["unit"] + " in young, middle-aged and old copies",
              FIG_XR)
        _oe_bar(fig, x0 + wt + 0.1, y + 0.12, wm - 0.14, norm, cmap, "shared vs random")
    else:
        _head(fig, axG, "G", "B's lineage sharing within young, middle-aged and old copies", FIG_XR)
    y -= gap + hH
    axH = _ax(fig, x0, y, wt, hH)
    _time_axis(axH, P, labels=not rows, xlabel=not rows, nums=nums)
    _bands(axH, P)
    if dt and P["prof"]:
        _panel_profile(axH, P)
        _head(fig, axH, "H", "Closest relatives in another subgenome over time, and every date estimate", FIG_XR)
    else:
        axH.set_yticks([])
        axH.spines["left"].set_visible(False)
        msg = ("one group: no split or merger to date" if P["N"] == 1 else
               "dating skipped" if not dt else "too few copy pairs to date")
        axH.text(0.5, 0.5, msg, transform=axH.transAxes, ha="center", va="center", fontsize=6.3, color="#666666")
        _head(fig, axH, "H", "Closest relatives in another subgenome over time", FIG_XR)
    if rows:
        y -= 0.05 + hR
        axR, axRm = _ax(fig, x0, y, wt, hR), _ax(fig, x0 + wt + 0.08, y, wm - 0.08, hR)
        _time_axis(axR, P, nums=nums)
        _bands(axR, P)
        _panel_estimates(fig, axR, axRm, P, rows, nums)


def fig_report(P, nums, outdir=None):
    """report.pdf's one page (see the module docstring) from the figure data P, its printed numbers recorded in nums
    (plot.style.Numbers); returns the figure (Agg canvas). Draw it inside _rc()."""
    import matplotlib as mpl
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.colors import LogNorm
    from matplotlib.figure import Figure
    lib, ph = P["lib"], P["ph"]
    W, H = FIG_W, FIG_H
    fig = Figure(figsize=(W, H))
    FigureCanvasAgg(fig)
    cmap, norm = mpl.colormaps["magma"], LogNorm(0.25, 4)
    N = P["N"]
    lab_n = ("groups given by the user" if ph.get("source") == "user" else "subgenomes") if N > 1 else ""
    x = 0.15
    t = fig.text(x / W, (H - 0.2) / H, P["label"], fontsize=9, fontstyle="italic", fontweight="bold", va="baseline")
    x += t.get_window_extent(fig.canvas.get_renderer()).width / fig.dpi + 0.12
    unph = len(lib.seqnames) - len(P["names"])
    verdict = str(ph.get("verdict", "one group"))
    short, _, detail = verdict.partition(":") if N == 1 else (verdict, "", "")
    groups = nums.add("groups", FIG_NUM.get(N, N) if N > 1 else "one", "report_data.json:N")
    if N == 1:
        nums.add("verdict", short.strip(), "report_data.json:ph.verdict")
    head = (f"·  {groups} {lab_n}" if N > 1 else f"·  {groups} group ({short.strip()})") + \
        f"  ·  {nums.add('phased', len(P['names']), 'report_data.json:names (count)')} {P['unit']} phased" + \
        (f" ({nums.add('unphased', unph, 'report_data.json:lib.seqnames minus names (count)')} too short or too "
         f"poor in LTR-RTs)" if unph else "") + \
        f"  ·  {nums.add('ltr_rts', f'{lib.n:,}', 'report_data.json:lib.n')} LTR retrotransposons"
    fig.text(x / W, (H - 0.2) / H, head, fontsize=7.5, va="baseline", color=FIG_INK)
    if detail.strip():
        fig.text(0.15 / W, (H - 0.33) / H, nums.add("verdict_detail", detail.strip()[:1].upper() + detail.strip()[1:],
                                                    "report_data.json:ph.verdict"), fontsize=5.8,
                 va="baseline", color="#666666")
    # A: karyotype across the top
    XL, XR = 0.55, W - 0.15
    axA = _ax(fig, XL, 6.3, XR - XL, 8.48 - _name_layout(fig, P)[2] - 6.3)     # names fit under the group heads
    xpos = _panel_karyotype(fig, axA, P, nums)
    _group_heads(fig, axA, P, xpos, 8.56, nums)
    _head(fig, axA, "A", ("Chromosomes painted by the group their LTR retrotransposon lineages point to" if
                          ph.get("source") == "user" else
                          "Chromosomes painted by the progenitor their LTR retrotransposon lineages come from")
          if P["tiles"] else "Chromosomes and their LTR retrotransposons (nothing to paint)", FIG_XL,
          dy=H - 0.47 - axA.get_position().y1 * H)
    items = ([("box", SG_COLORS[g], f"{group_label(g)} lineages") for g in range(N)] if N <= 4 else
             [("box", SG_COLORS[0], "subgenome lineages (group colors)")]) if P["tiles"] else []
    if P["tiles"]:
        items += [("hatch", None, "too few copies to tell")]
    if P["calls"]:
        items += [("tri", SG_COLORS[0], "called exchange (donor color)")]
    if P["dotplot"] and _homoeolog_numbers(ph):
        items += [("num", None, "homoeolog set")]
    items += [("dens", None, "LTR-RT density")]
    room = W - 0.1 - (XL - 0.13)                                       # tighten the key when it overflows
    need = lambda gap: sum(0.17 + _text_w(fig, lab, 6.2) + gap for _, _, lab in items)
    _key_row(fig, XL - 0.13, H - 0.61, items, gap=0.22 if need(0.22) <= room else max(0.06, 0.22 - (need(0.22) - room)
                                                                                         / max(len(items) - 1, 1)))
    # left column: B, C, D
    axB = _ax(fig, 0.72, 3.88, 1.62, 1.62)
    if P["dotplot"]:                                                    # homoeolog evidence
        _panel_homoeologs(fig, axB, P)
        _head(fig, axB, "B", f"Old insertions shared by {P['unit']} of two subgenomes", FIG_XL, dy=0.17)
        hb = ph["homoeologs"].get("bio")
        kb = [("dot", FIG_ACCENT, "collinear, homoeologs"), ("dot", "#b9b9b9", "other"), ("frame", FIG_ACCENT, "pair called")]
        if hb and hb["segments"]:
            kb.append(("dot", FIG_EXCH, "near-identical (exchanged)"))
        if hb and hb["transloc"]:
            kb.append(("frame", FIG_EXCH, "translocated"))
        _key_row(fig, 0.62, 3.53, kb[:3], size=5.6, gap=0.1)
        if len(kb) > 3:
            _key_row(fig, 0.62, 3.40, kb[3:], size=5.6, gap=0.1)
    else:
        _panel_sharing(fig, axB, P, cmap, norm)
        _head(fig, axB, "B", f"Transposon lineages shared by each pair of {P['unit']}", FIG_XL, dy=0.17)
        _oe_bar(fig, 0.75, 3.55, 1.5, norm, cmap, "lineages shared vs random")
    axC = _ax(fig, 0.85, 1.78, 1.85, 1.3)
    _panel_pca(fig, axC, P, nums)
    _head(fig, axC, "C", f"Placement of each {P['unit'][:-1]} in independent halves of the library", FIG_XL,
          dy=0.12)
    axDr, axDs = _ax(fig, 1.45, 0.42, 0.72, 0.58), _ax(fig, 2.3, 0.42, 0.55, 0.58)
    _panel_tests(fig, axDr, axDs, P, nums)
    _head(fig, axDr, "D", _tests_header(P), FIG_XL, dy=0.36)
    # right column: E-H on one time axis
    _panel_history(fig, P, 3.75, 2.25, 0.85, 5.8, cmap, norm, nums)
    import socket
    fig.text((W - 0.1) / W, 0.06 / H, f"{PROGRAM} · {socket.gethostname()}:"
             f"{Path(outdir).resolve() if outdir else ''} · {time.strftime('%Y-%m-%d')}", fontsize=4.5,
             ha="right", va="bottom", color="#999999")
    return fig
