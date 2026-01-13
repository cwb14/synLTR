#!/usr/bin/env python3
"""
Takes genome FAI and the dedup results of the FAI to plot LTR-RT annotation features.
(1) age density.
(2) LTR-RT size dist.
(3) LTR size dist.
(4) internal size dist.
(5) Chrom locations.

python3 ltrharvest_plot.py --k2p-xmax 0.2 --bin 200 --legend-page --stacked
python3 ltrharvest_plot.py --k2p-xmax 0.2 --bin 200 --legend-page
# Large genomes need --chr-rasterize.
python ltrharvest_plot.py --species B73 --aln-suffix _kmer2ltr_dedup --outpdf B73_ltr3.pdf --bin 100 --k2p-xmax 0.1 --legend-page --timing --chr-rasterize --chr-merge-gap 50

LTR-RT multi-page PDF plots across species.

Per species, the script can produce 5 plot types:
  1) K2P density (ONE overall density; area partitioned by TE type; uses column11)
  2) Chromosome distribution (karyotype-style; TE intervals colored by type; uses FAI + column1)
     Optional: overlay dots for any user-specified TE type(s)
  3) Full-length LTR-RT size distribution (stacked bars; 50 bp bins; from col1 start/end)
  4) LTR size distribution (stacked bars; 50 bp bins; col3)
  5) Internal size distribution (stacked bars; 50 bp bins; full - 2*ltr)

Legend behavior:
- If --legend-page is set:
    * A single legend-only page is added (colors are global & consistent)
    * NO legends are embedded in any plots.
- If --legend-page is NOT set:
    * Each page includes an embedded legend (fit into the figure using a side panel).
    * This applies to both --stacked and non-stacked modes.

Stacked panels mode (--stacked):
- All K2P density panels on a single page (shared x/y axes)
- All chromosome distribution panels on a single page (shared x axis)
- All size panels on a single page for each size type (shared x/y axes)

Density correctness:
- Per species, KDE per type is computed using SAME bandwidth as global KDE, then weighted by n_type/n_total.
  Stacked fill sums to the global density estimate, avoiding disproportionate bumps from tiny-count types.

Assumptions:
- Alignment file is whitespace-delimited with >= 11 columns
- Column 1: "Chr:start-end#TYPE"
- Column 3: LTR length (bp)
- Column 11: K2P divergence
- Full-length LTR-RT length is (end - start) (NOT +1), per your example (100-200 => 100)
- Ignore scaffolds: any chromosome name containing "sca" (case-insensitive) excluded in FAI and TE intervals
"""

import re
import math
import argparse
import sys
import time
from datetime import datetime
from collections import defaultdict, OrderedDict

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Patch


COL1_RE = re.compile(r"^(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)#(?P<typ>.+)$")


# -----------------------------
# Timing / logging
# -----------------------------
def _now_str():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def tlog(msg, enabled=True):
    """Timestamped log to stderr."""
    if enabled:
        print(f"[{_now_str()}] {msg}", file=sys.stderr, flush=True)


class Timer:
    """Context manager for timing blocks."""
    def __init__(self, label, enabled=True):
        self.label = label
        self.enabled = enabled
        self.t0 = None

    def __enter__(self):
        if self.enabled:
            tlog(f"START: {self.label}", enabled=True)
            self.t0 = time.perf_counter()
        return self

    def __exit__(self, exc_type, exc, tb):
        if self.enabled and self.t0 is not None:
            dt = time.perf_counter() - self.t0
            tlog(f"END:   {self.label}  (elapsed {dt:.3f}s)", enabled=True)


# -----------------------------
# IO / parsing
# -----------------------------
def read_fai_lengths(fai_path):
    """Return OrderedDict {chrom: length}, excluding seqnames with 'sca' (case-insensitive)."""
    chrom2len = {}
    with open(fai_path, "r") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            chrom = parts[0]
            if "sca" in chrom.lower():
                continue
            try:
                length = int(parts[1])
            except ValueError:
                continue
            chrom2len[chrom] = length

    def sort_key(name):
        m = re.search(r"(\d+)$", name)
        return (name[:m.start()] if m else name, int(m.group(1)) if m else 10**18, name)

    return OrderedDict((k, chrom2len[k]) for k in sorted(chrom2len, key=sort_key))


def parse_alignment_file(path):
    """
    Parse alignment file:
      col1 => chrom, start, end, type
      col3 => LTR length
      col11 => K2P
    """
    records = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 11:
                continue

            m = COL1_RE.match(parts[0])
            if not m:
                continue

            chrom = m.group("chrom")
            start = int(m.group("start"))
            end = int(m.group("end"))
            te_type = m.group("typ")

            # per your example: 100-200 => 100
            full_len = max(0, end - start)

            try:
                ltr_len = int(float(parts[2]))  # column3
            except ValueError:
                ltr_len = None

            internal_len = None
            if ltr_len is not None:
                internal_len = full_len - 2 * ltr_len

            try:
                k2p = float(parts[10])  # column11
            except ValueError:
                k2p = None

            records.append(
                {
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "type": te_type,
                    "k2p": k2p,
                    "full_len": full_len,
                    "ltr_len": ltr_len,
                    "internal_len": internal_len,
                }
            )
    return records


def compute_chrom_density(records, chrom2len, window_bp=1_000_000, step_bp=500_000):
    """
    Sliding-window density per chromosome per TE type.
    Returns:
      dens[chrom][type] = np.array(density_fraction_per_window)
      centers_bp[chrom] = np.array(window_center_positions_bp)
    Density is fraction of bases in window covered by that TE type (0..1).
    """
    dens = {}
    centers = {}

    if window_bp <= 0 or step_bp <= 0:
        return dens, centers

    # index intervals by chrom
    by_chrom = defaultdict(list)
    for r in records:
        chrom = r["chrom"]
        if "sca" in chrom.lower():
            continue
        if chrom not in chrom2len:
            continue
        start = max(0, int(r["start"]))
        end = min(int(r["end"]), int(chrom2len[chrom]))
        if end <= start:
            continue
        by_chrom[chrom].append((start, end, r["type"]))

    for chrom, L in chrom2len.items():
        if chrom not in by_chrom:
            continue

        # windows: [i*step, i*step+window)
        # include last window whose start < L
        nwin = int(math.floor((max(L - 1, 0)) / step_bp)) + 1
        w_starts = np.arange(0, nwin * step_bp, step_bp, dtype=np.int64)
        w_ends = w_starts + window_bp
        centers_bp = w_starts + window_bp / 2.0

        # accumulate covered bp per window per type
        acc = defaultdict(lambda: np.zeros(nwin, dtype=np.float64))

        for start, end, typ in by_chrom[chrom]:
            # window indices that could overlap this interval
            i0 = max(0, int((start - window_bp) // step_bp))
            i1 = min(nwin - 1, int(end // step_bp))
            for i in range(i0, i1 + 1):
                ws = int(w_starts[i])
                we = int(w_ends[i])
                ov = max(0, min(end, we) - max(start, ws))
                if ov:
                    acc[typ][i] += ov

        # convert bp-covered to density fraction in the window (cap at 1)
        dens[chrom] = {}
        for typ, covered in acc.items():
            d = covered / float(window_bp)
            d[d < 0] = 0
            d[d > 1] = 1
            dens[chrom][typ] = d

        centers[chrom] = centers_bp

    return dens, centers

def compute_chrom_k2p_window_stats(records, chrom2len, window_bp=1_000_000, step_bp=500_000):
    """
    Sliding-window mean K2P per chromosome (all types pooled).

    Each TE contributes its K2P to every window whose [start, end) contains the TE midpoint.
    (This matches the idea of sliding windows with overlap without double-counting by bp.)

    Returns:
      mean_k2p[chrom] = np.array(mean per window; NaN where n==0)
      err_k2p[chrom]  = np.array(95% CI half-width; NaN where n<2)
      n_k2p[chrom]    = np.array(count per window)
      centers_bp[chrom] = np.array(window centers in bp)
    """
    mean_k2p = {}
    err_k2p = {}
    n_k2p = {}
    centers = {}

    if window_bp <= 0 or step_bp <= 0:
        return mean_k2p, err_k2p, n_k2p, centers

    # index per chrom
    by_chrom = defaultdict(list)
    for r in records:
        chrom = r["chrom"]
        if "sca" in chrom.lower():
            continue
        if chrom not in chrom2len:
            continue
        k = r.get("k2p", None)
        if k is None or (not np.isfinite(k)) or k < 0:
            continue
        start = max(0, int(r["start"]))
        end = min(int(r["end"]), int(chrom2len[chrom]))
        if end <= start:
            continue
        mid = 0.5 * (start + end)
        by_chrom[chrom].append((mid, float(k)))

    for chrom, L in chrom2len.items():
        if chrom not in by_chrom:
            continue

        nwin = int(math.floor((max(L - 1, 0)) / step_bp)) + 1
        w_starts = np.arange(0, nwin * step_bp, step_bp, dtype=np.int64)
        w_ends = w_starts + window_bp
        centers_bp = w_starts + window_bp / 2.0

        s = np.zeros(nwin, dtype=np.float64)
        ss = np.zeros(nwin, dtype=np.float64)
        c = np.zeros(nwin, dtype=np.int64)

        # For each midpoint, add to all windows where ws <= mid < we.
        # Window i has ws = i*step, we = ws+window.
        # Solve for i:
        #   i*step <= mid < i*step + window
        #   (mid - window) < i*step <= mid
        # So:
        #   i0 = floor((mid - window)/step) + 1  (smallest integer i satisfying i*step > mid-window)
        #   i1 = floor(mid/step)
        for mid, k in by_chrom[chrom]:
            i1 = int(math.floor(mid / step_bp))
            i0 = int(math.floor((mid - window_bp) / step_bp) + 1)
            if i1 < 0 or i0 > (nwin - 1):
                continue
            i0 = max(i0, 0)
            i1 = min(i1, nwin - 1)
            if i1 < i0:
                continue
            # vectorized update on slice
            s[i0:i1+1] += k
            ss[i0:i1+1] += k * k
            c[i0:i1+1] += 1

        mean = np.full(nwin, np.nan, dtype=np.float64)
        err = np.full(nwin, np.nan, dtype=np.float64)

        ok = c > 0
        mean[ok] = s[ok] / c[ok]

        # sample SD only defined for n>=2; CI uses SE = sd/sqrt(n)
        ok2 = c >= 2
        # var = (sumsq - sum^2/n) / (n-1)
        var = np.zeros(nwin, dtype=np.float64)
        var[ok2] = (ss[ok2] - (s[ok2] ** 2) / c[ok2]) / (c[ok2] - 1)
        var[var < 0] = 0.0
        sd = np.sqrt(var)
        se = np.zeros(nwin, dtype=np.float64)
        se[ok2] = sd[ok2] / np.sqrt(c[ok2])
        err[ok2] = 1.96 * se[ok2]  # 95% CI half-width

        mean_k2p[chrom] = mean
        err_k2p[chrom] = err
        n_k2p[chrom] = c
        centers[chrom] = centers_bp

    return mean_k2p, err_k2p, n_k2p, centers


def filter_density_types(dens_by_chrom, threshold_percent=None, top_n=14):
    """
    dens_by_chrom: dens[chrom][type] = np.array(0..1)
    Returns list of kept types sorted by max density desc.
    """
    max_by_type = defaultdict(float)
    for chrom, m in dens_by_chrom.items():
        for t, arr in m.items():
            if arr.size:
                mx = float(np.nanmax(arr))
                if np.isfinite(mx):
                    max_by_type[t] = max(max_by_type[t], mx)

    types = list(max_by_type.keys())
    # threshold
    if threshold_percent is not None:
        thr = float(threshold_percent) / 100.0
        types = [t for t in types if max_by_type[t] > thr]

    # sort by max density desc
    types.sort(key=lambda t: max_by_type[t], reverse=True)

    # cap
    if top_n is not None and int(top_n) > 0 and len(types) > int(top_n):
        types = types[: int(top_n)]

    return types


def style_minimal_axes(ax):
    # simple "theme_minimal-ish" without global matplotlib styles
    ax.grid(True, alpha=0.2)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_alpha(0.25)
    ax.tick_params(labelsize=8)


def plot_chrom_density_facets(
    pdf, sp, chrom2len, records, type_colors,
    window_bp=1_000_000, step_bp=500_000,
    threshold_percent=None, top_n=14,
    legend_page=False, embedded_legend_cols=2,
    timing=False,
):
    """
    One PDF page for a species: facet grid of chromosomes; line per type.
    X in Mb, Y in percent.
    """
    with Timer(f"[{sp}] compute_chrom_density(window={window_bp}, step={step_bp})", enabled=timing):
        dens, centers = compute_chrom_density(records, chrom2len, window_bp=window_bp, step_bp=step_bp)

    if not dens:
        fig, ax = plt.subplots(figsize=(11, 6))
        ax.text(0.5, 0.5, f"{sp}: No density data", ha="center", va="center")
        ax.set_axis_off()
        pdf.savefig(fig)
        plt.close(fig)
        return

    with Timer(f"[{sp}] filter_density_types(thr={threshold_percent}, top_n={top_n})", enabled=timing):
        kept_types = filter_density_types(dens, threshold_percent=threshold_percent, top_n=top_n)

    if not kept_types:
        fig, ax = plt.subplots(figsize=(11, 6))
        ax.text(0.5, 0.5, f"{sp}: No types pass filters", ha="center", va="center")
        ax.set_axis_off()
        pdf.savefig(fig)
        plt.close(fig)
        return

    chroms = [c for c in chrom2len.keys() if c in dens]  # preserve sorted chrom order
    r, c = panel_grid(len(chroms))

    with Timer(f"[{sp}] plot_chrom_density_facets(render + save)", enabled=timing):
        fig, axes = plt.subplots(r, c, figsize=(11, 8.5), squeeze=False)
        axes_flat = axes.ravel()

        # shared y max across facets for readability
        y_max = 0.0
        for chrom in chroms:
            for t in kept_types:
                arr = dens[chrom].get(t)
                if arr is not None and arr.size:
                    y_max = max(y_max, float(np.nanmax(arr) * 100.0))
        y_max = y_max * 1.05 if y_max > 0 else None

        for i, chrom in enumerate(chroms):
            ax = axes_flat[i]
            x_mb = centers[chrom] / 1e6

            for t in kept_types:
                y = dens[chrom].get(t)
                if y is None:
                    continue
                ax.plot(x_mb, y * 100.0, linewidth=1.2, alpha=0.9, color=type_colors.get(t, (0, 0, 0, 1)))

            ax.set_title(chrom, fontsize=9)
            ax.set_xlim(0, float(chrom2len[chrom]) / 1e6)
            if y_max is not None:
                ax.set_ylim(0, y_max)
            style_minimal_axes(ax)

            # label only left/bottom panels
            if i % c == 0:
                ax.set_ylabel("% in window", fontsize=9)
            else:
                ax.set_ylabel("")
            if i >= (r - 1) * c:
                ax.set_xlabel("Position (Mb)", fontsize=9)
            else:
                ax.set_xlabel("")

        # turn off unused
        for j in range(len(chroms), len(axes_flat)):
            axes_flat[j].set_axis_off()

        step_mb = step_bp / 1e6
        win_mb = window_bp / 1e6
        title = f"{sp}: Chromosome distribution (density-style)  |  window={win_mb:g}Mb step={step_mb:g}Mb"
        if threshold_percent is not None:
            title += f"  |  thr>{float(threshold_percent):g}%"
        if top_n and int(top_n) > 0:
            title += f"  |  top={int(top_n)}"
        fig.suptitle(title, fontsize=12)

        # legend handling matches your existing philosophy
        if not legend_page:
            add_embedded_legend(fig, {t: type_colors[t] for t in kept_types if t in type_colors}, ncol=embedded_legend_cols, fontsize=7)
            fig.tight_layout(rect=[0, 0, 0.72, 0.96])
        else:
            fig.tight_layout(rect=[0, 0, 1, 0.96])

        pdf.savefig(fig)
        plt.close(fig)

def plot_chrom_k2p_mean_facets(
    pdf, sp, chrom2len, records,
    window_bp=1_000_000, step_bp=500_000,
    legend_page=False,  # kept for API symmetry; no legend needed here
    timing=False,
):
    """
    One PDF page per species: facet grid of chromosomes;
    plot windowed mean K2P with 95% CI error bars.
    """
    with Timer(f"[{sp}] compute_chrom_k2p_window_stats(window={window_bp}, step={step_bp})", enabled=timing):
        mean_k2p, err_k2p, n_k2p, centers = compute_chrom_k2p_window_stats(
            records, chrom2len, window_bp=window_bp, step_bp=step_bp
        )

    if not mean_k2p:
        fig, ax = plt.subplots(figsize=(11, 6))
        ax.text(0.5, 0.5, f"{sp}: No K2P window stats", ha="center", va="center")
        ax.set_axis_off()
        pdf.savefig(fig)
        plt.close(fig)
        return

    chroms = [c for c in chrom2len.keys() if c in mean_k2p]
    r, c = panel_grid(len(chroms))

    with Timer(f"[{sp}] plot_chrom_k2p_mean_facets(render + save)", enabled=timing):
        fig, axes = plt.subplots(r, c, figsize=(11, 8.5), squeeze=False)
        axes_flat = axes.ravel()

        # shared y-lims
        y_min = np.inf
        y_max = -np.inf
        for chrom in chroms:
            m = mean_k2p[chrom]
            e = err_k2p[chrom]

            ok = np.isfinite(m)
            if np.any(ok):
                y_min = min(y_min, float(np.nanmin(m[ok])))
                y_max = max(y_max, float(np.nanmax(m[ok])))

            ok2 = ok & np.isfinite(e)
            if np.any(ok2):
                y_min = min(y_min, float(np.nanmin(m[ok2] - e[ok2])))
                y_max = max(y_max, float(np.nanmax(m[ok2] + e[ok2])))


        if not np.isfinite(y_min) or not np.isfinite(y_max):
            y_min, y_max = 0.0, 1.0
        pad = 0.05 * (y_max - y_min) if (y_max > y_min) else 0.01
        y_min = max(0.0, y_min - pad)
        y_max = y_max + pad

        for i, chrom in enumerate(chroms):
            ax = axes_flat[i]
            x_mb = centers[chrom] / 1e6
            m = mean_k2p[chrom]
            e = err_k2p[chrom]
            n = n_k2p[chrom]

            ok = np.isfinite(m)
            if np.any(ok):
                # only draw error bars where defined (n>=2)
                ok_err = ok & np.isfinite(e)
                if np.any(ok_err):
                    ax.errorbar(
                        x_mb[ok_err],
                        m[ok_err],
                        yerr=e[ok_err],
                        fmt="o",
                        markersize=2.5,
                        linewidth=0.8,
                        capsize=2,
                    )
                # also show mean-only points where n==1 (no CI)
                ok_noerr = ok & (~np.isfinite(e))
                if np.any(ok_noerr):
                    ax.plot(x_mb[ok_noerr], m[ok_noerr], "o", markersize=2.0, alpha=0.7)

            ax.set_title(chrom, fontsize=9)
            ax.set_xlim(0, float(chrom2len[chrom]) / 1e6)
            ax.set_ylim(y_min, y_max)
            style_minimal_axes(ax)

            # label only left/bottom panels
            if i % c == 0:
                ax.set_ylabel("Mean K2P", fontsize=9)
            else:
                ax.set_ylabel("")
            if i >= (r - 1) * c:
                ax.set_xlabel("Position (Mb)", fontsize=9)
            else:
                ax.set_xlabel("")

            # tiny annotation: n summary (optional but helpful)
            if n is not None and n.size:
                nn = int(np.nanmax(n)) if np.any(n > 0) else 0
                ax.text(0.98, 0.02, f"max n={nn}", transform=ax.transAxes,
                        ha="right", va="bottom", fontsize=7, alpha=0.7)

        for j in range(len(chroms), len(axes_flat)):
            axes_flat[j].set_axis_off()

        step_mb = step_bp / 1e6
        win_mb = window_bp / 1e6
        fig.suptitle(
            f"{sp}: Windowed mean K2P (95% CI)  |  window={win_mb:g}Mb step={step_mb:g}Mb",
            fontsize=12
        )
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        pdf.savefig(fig)
        plt.close(fig)


# -----------------------------
# Color mapping + legend utils
# -----------------------------
def make_type_colors(all_types):
    """Fixed, consistent mapping from TE type -> RGBA (stable across all pages)."""
    types_sorted = sorted(all_types)
    cmaps = [plt.get_cmap("tab20"), plt.get_cmap("tab20b"), plt.get_cmap("tab20c")]
    palette = []
    for cm in cmaps:
        palette.extend([cm(i) for i in range(cm.N)])

    colors = {}
    for i, t in enumerate(types_sorted):
        colors[t] = palette[i % len(palette)]
    return colors


def build_legend_handles(type_colors):
    types_sorted = sorted(type_colors.keys())
    return [Patch(facecolor=type_colors[t], edgecolor="none", label=t, alpha=0.85) for t in types_sorted]


def make_legend_page(pdf, type_colors, ncol=2, title="TE type color legend"):
    """Legend-only page."""
    handles = build_legend_handles(type_colors)
    n = len(handles)
    rows = max(1, math.ceil(n / ncol))
    fig_h = min(20, max(6, rows * 0.3))
    fig, ax = plt.subplots(figsize=(11, fig_h))
    ax.set_axis_off()
    ax.set_title(title, fontsize=14, pad=10)
    ax.legend(handles=handles, loc="upper left", frameon=False, ncol=ncol, fontsize=8)
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def add_embedded_legend(fig, type_colors, ncol=2, fontsize=7):
    """
    Add an embedded legend into the figure without overlapping axes by reserving a right-side margin.
    This works for both single-plot pages and multi-panel pages.
    """
    handles = build_legend_handles(type_colors)
    # Reserve space on right for legend
    fig.subplots_adjust(right=0.72)
    fig.legend(
        handles=handles,
        loc="center left",
        bbox_to_anchor=(0.74, 0.5),
        frameon=False,
        ncol=ncol,
        fontsize=fontsize,
    )


# -----------------------------
# KDE utilities
# -----------------------------
def silverman_bandwidth(x):
    """Silverman's rule-of-thumb bandwidth."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    n = x.size
    if n <= 1:
        return 0.01
    std = np.std(x, ddof=1)
    if not np.isfinite(std) or std <= 0:
        std = 1.0
    h = 1.06 * std * (n ** (-1 / 5))
    if not np.isfinite(h) or h <= 0:
        h = max(1e-6, std * 0.1)
    return h


def gaussian_kde_numpy(x, grid, bandwidth):
    """Gaussian KDE in numpy, using a provided bandwidth."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    n = x.size
    if n == 0:
        return np.zeros_like(grid, dtype=float)

    h = float(bandwidth)
    if not np.isfinite(h) or h <= 0:
        h = 1e-6

    diffs = (grid[:, None] - x[None, :]) / h
    dens = np.exp(-0.5 * diffs**2).sum(axis=1) / (n * h * np.sqrt(2 * np.pi))
    return dens


# -----------------------------
# Computation for plots
# -----------------------------
def compute_species_k2p_stack(records, k2p_xmax=None, grid_n=700):
    """
    Compute (xgrid, dens_global, dens_weighted_by_type) for one species.
    dens_weighted_by_type[t] sums to dens_global over types (numerically).
    """
    per_type = defaultdict(list)
    all_vals = []
    for r in records:
        k = r["k2p"]
        if k is None or not np.isfinite(k) or k < 0:
            continue
        per_type[r["type"]].append(k)
        all_vals.append(k)

    if len(all_vals) == 0:
        return None

    all_vals = np.asarray(all_vals, dtype=float)
    data_max = float(np.nanmax(all_vals))
    if not np.isfinite(data_max):
        data_max = 0.1

    xmin = 0.0
    if k2p_xmax is not None:
        xmax = float(k2p_xmax)
    else:
        pad = 0.05 * data_max if data_max > 0 else 0.01
        xmax = max(data_max + pad, 0.05)

    xgrid = np.linspace(xmin, xmax, grid_n)
    bw = silverman_bandwidth(all_vals)
    dens_global = gaussian_kde_numpy(all_vals, xgrid, bw)

    n_total = float(len(all_vals))
    dens_weighted = {}
    for t, vals in per_type.items():
        vals = np.asarray(vals, dtype=float)
        n_t = float(vals.size)
        if n_t == 0:
            continue
        dens_t = gaussian_kde_numpy(vals, xgrid, bw)
        dens_weighted[t] = dens_t * (n_t / n_total)

    return xgrid, dens_global, dens_weighted


def compute_size_by_type(records, key):
    out = defaultdict(list)
    for r in records:
        v = r.get(key, None)
        if v is None or not np.isfinite(v) or v < 0:
            continue
        out[r["type"]].append(float(v))
    return out


def compute_hist_counts(values_by_type, bins):
    counts_by_type = {}
    for t, vals in values_by_type.items():
        vv = np.asarray([v for v in vals if np.isfinite(v) and v >= 0], dtype=float)
        if vv.size == 0:
            continue
        counts, _ = np.histogram(vv, bins=bins)
        counts_by_type[t] = counts
    return counts_by_type


# -----------------------------
# Plotting (no legends here; legend handled by caller)
# -----------------------------
def plot_k2p_panel(ax, sp_label, k2p_data, type_colors, y_max=None, show_global_line=False):
    ax.set_title(sp_label, fontsize=10)

    if k2p_data is None:
        ax.text(0.5, 0.5, "No K2P", ha="center", va="center", fontsize=9)
        ax.set_axis_off()
        return

    xgrid, dens_global, dens_weighted = k2p_data
    types_order = sorted(
        dens_weighted.keys(),
        key=lambda t: float(np.trapz(dens_weighted[t], xgrid)),
        reverse=True,
    )

    base = np.zeros_like(xgrid)
    for t in types_order:
        d = dens_weighted[t]
        ax.fill_between(xgrid, base, base + d, color=type_colors[t], alpha=0.65, linewidth=0.0)
        base = base + d

    if show_global_line:
        ax.plot(xgrid, dens_global, linewidth=1.0)

    ax.set_xlim(float(xgrid[0]), float(xgrid[-1]))
    if y_max is not None and np.isfinite(y_max) and y_max > 0:
        ax.set_ylim(0, y_max)

    ax.grid(True, alpha=0.2)
    ax.tick_params(labelsize=8)


from collections import defaultdict

def merge_ranges(xranges, gap=0):
    """Merge sorted (start, width) into non-overlapping (start, width)."""
    if not xranges:
        return []
    xranges = sorted(xranges, key=lambda x: x[0])
    merged = []
    s0, w0 = xranges[0]
    e0 = s0 + w0
    for s, w in xranges[1:]:
        e = s + w
        if s <= e0 + gap:   # overlap / near-touch
            e0 = max(e0, e)
        else:
            merged.append((s0, e0 - s0))
            s0, e0 = s, e
    merged.append((s0, e0 - s0))
    return merged

def plot_chrom_panel(ax, sp_label, chrom2len, records, type_colors, dot_types, x_max=None,
                    merge_gap_bp=0, rasterize=True):
    ax.set_title(sp_label, fontsize=10)

    if not chrom2len:
        ax.text(0.5, 0.5, "No chromosomes", ha="center", va="center", fontsize=9)
        ax.set_axis_off()
        return

    # intervals_by_chrom_type[chrom][type] = [(start, width), ...]
    intervals_by_chrom_type = defaultdict(lambda: defaultdict(list))
    dots = defaultdict(list)
    dot_set = set(dot_types) if dot_types else set()

    for r in records:
        chrom = r["chrom"]
        if "sca" in chrom.lower():
            continue
        if chrom not in chrom2len:
            continue

        start = max(0, int(r["start"]))
        end = min(int(r["end"]), chrom2len[chrom])
        if end <= start:
            continue

        t = r["type"]
        intervals_by_chrom_type[chrom][t].append((start, end - start))

        if dot_set and t in dot_set:
            dots[chrom].append((start + end) / 2.0)

    chrom_names = list(chrom2len.keys())
    n = len(chrom_names)
    ystep = 10
    height = 6
    y_positions = {c: (n - 1 - i) * ystep for i, c in enumerate(chrom_names)}

    # chromosome backbones (cheap)
    for chrom in chrom_names:
        y = y_positions[chrom]
        L = chrom2len[chrom]
        ax.broken_barh([(0, L)], (y, height),
                       facecolors=(0, 0, 0, 0.08),
                       edgecolors=(0, 0, 0, 0.12),
                       linewidth=0.5)

    # TE intervals (batched + optionally merged)
    for chrom, type_map in intervals_by_chrom_type.items():
        y = y_positions[chrom]
        for t, xranges in type_map.items():
            if merge_gap_bp is not None and merge_gap_bp >= 0:
                xranges = merge_ranges(xranges, gap=int(merge_gap_bp))
            coll = ax.broken_barh(
                xranges,
                (y, height),
                facecolors=type_colors.get(t, (0, 0, 0, 1)),
                edgecolors="none",
                alpha=0.85,
            )
            if rasterize:
                coll.set_rasterized(True)  # huge PDF speedup

    # dots (usually not the bottleneck)
    if dot_set:
        for chrom, mids in dots.items():
            if not mids:
                continue
            y = y_positions[chrom] + height + 0.8
            ax.scatter(mids, [y] * len(mids), s=6, rasterized=rasterize)

    max_len = max(chrom2len.values())
    xmax = x_max if x_max is not None else max_len
    ax.set_xlim(0, xmax * 1.02)
    ax.set_ylim(-ystep, n * ystep + 6)
    ax.set_yticks([y_positions[c] + height / 2 for c in chrom_names])
    ax.set_yticklabels(chrom_names, fontsize=7)
    ax.grid(True, axis="x", alpha=0.2)
    ax.tick_params(labelsize=8)

def plot_size_panel(ax, sp_label, counts_by_type, bins, type_colors, y_max=None):
    ax.set_title(sp_label, fontsize=10)

    if not counts_by_type:
        ax.text(0.5, 0.5, "No values", ha="center", va="center", fontsize=9)
        ax.set_axis_off()
        return

    types_order = sorted(counts_by_type.keys(), key=lambda t: int(np.sum(counts_by_type[t])), reverse=True)

    bottom = np.zeros(len(bins) - 1, dtype=float)
    for t in types_order:
        c = counts_by_type[t]
        ax.bar(
            bins[:-1],
            c,
            width=(bins[1] - bins[0]),
            bottom=bottom,
            align="edge",
            alpha=0.8,
            edgecolor="none",
            color=type_colors[t],
        )
        bottom += c

    ax.set_xlim(float(bins[0]), float(bins[-1]))
    if y_max is not None and np.isfinite(y_max) and y_max > 0:
        ax.set_ylim(0, y_max)
    ax.grid(True, axis="y", alpha=0.2)
    ax.tick_params(labelsize=8)


# -----------------------------
# Panel layout helpers
# -----------------------------
def panel_grid(n):
    if n <= 1:
        return 1, 1
    if n == 2:
        return 2, 1
    if n <= 4:
        return 2, 2
    if n <= 6:
        return 3, 2
    if n <= 9:
        return 3, 3
    return math.ceil(n / 3), 3


def main():
    ap = argparse.ArgumentParser(description="Generate multi-page LTR-RT plots PDF across species.")
    ap.add_argument(
        "--species",
        nargs="+",
        default=["Aaren", "Ahall", "Alyra", "Athal", "Chisp"],
        help="Species prefixes (default: Aaren Ahall Alyra Athal Chisp)",
    )
    ap.add_argument("--fai-suffix", default=".fa.fai", help="Suffix for FAI files (default: .fa.fai)")
    ap.add_argument(
        "--aln-suffix",
        default="_ltrharvest6_kmer2ltr_dedup",
        help="Suffix for alignment files (default: _ltrharvest6_kmer2ltr_dedup)",
    )
    ap.add_argument("--outpdf", default="ltrrt_summary_plots.pdf", help="Output PDF filename")
    ap.add_argument("--bin", type=int, default=50, help="Bin size for size hists (default: 50 bp)")
    ap.add_argument(
        "--k2p-xmax",
        type=float,
        default=None,
        help="If set, fix K2P density x-axis max to this value (x-axis will be [0, xmax])",
    )
    ap.add_argument("--chr-merge-gap", type=int, default=0,
                help="Merge TE intervals if separated by <= this many bp (default: 0).")
    ap.add_argument("--chr-rasterize", action="store_true",
                help="Rasterize TE rectangles in chromosome plot for speed (recommended).")

    ap.add_argument(
        "--dot-type",
        action="append",
        default=[],
        help='Overlay dots for this TE type on chromosome plots (repeatable). Example: --dot-type "LTR/Gypsy/CRM"',
    )
    ap.add_argument(
        "--stacked",
        action="store_true",
        help="Stack like plots onto single pages for easy cross-species comparison (shared axes).",
    )
    ap.add_argument(
        "--legend-page",
        action="store_true",
        help="Add a legend-only page and omit legends from all plots.",
    )
    ap.add_argument("--legend-cols", type=int, default=2, help="Columns for legend (default: 2)")
    ap.add_argument(
        "--embedded-legend-cols",
        type=int,
        default=2,
        help="Columns for embedded legend when --legend-page is NOT used (default: 2)",
    )
    ap.add_argument(
        "--density-global-line",
        action="store_true",
        help="Draw a thin outline of the global density on top of the stacked fill (verification).",
    )
    ap.add_argument("--dens-window", type=int, default=1_000_000, help="Window size for chrom density (bp). Default: 1,000,000")
    ap.add_argument("--dens-step", type=int, default=500_000, help="Step size for chrom density (bp). Default: 500,000")
    ap.add_argument("--dens-threshold", type=float, default=None, help="Exclude types whose max density percent never exceeds this value (0-100).")
    ap.add_argument("--dens-top-types", type=int, default=14, help="Keep only top N types by max density (default: 14). Use 0 for no cap.")
    ap.add_argument(
        "--timing",
        action="store_true",
        help="Print timestamped timing logs to stderr to identify bottlenecks.",
    )

    args = ap.parse_args()
    timing = bool(args.timing)

    with Timer("Load species: read FAI + parse alignment + collect global TE types", enabled=timing):
        # Load species and build global TE type set for consistent colors
        species_data = {}
        all_types = set()
        for sp in args.species:
            fai = f"{sp}{args.fai_suffix}"
            aln = f"{sp}{args.aln_suffix}"

            with Timer(f"[{sp}] read_fai_lengths({fai})", enabled=timing):
                chrom2len = read_fai_lengths(fai)

            with Timer(f"[{sp}] parse_alignment_file({aln})", enabled=timing):
                records = parse_alignment_file(aln)

            with Timer(f"[{sp}] collect types + store species_data", enabled=timing):
                for r in records:
                    all_types.add(r["type"])
                species_data[sp] = {"chrom2len": chrom2len, "records": records}

            tlog(f"[{sp}] parsed records: {len(species_data[sp]['records'])}", enabled=timing)
            tlog(f"[{sp}] chromosomes kept (non-sca): {len(species_data[sp]['chrom2len'])}", enabled=timing)

    with Timer("Build global type_colors", enabled=timing):
        type_colors = make_type_colors(all_types)
        tlog(f"Global TE types: {len(type_colors)}", enabled=timing)

    with Timer(f"Open PDF for writing: {args.outpdf}", enabled=timing):
        with PdfPages(args.outpdf) as pdf:
            # Legend-only page
            if args.legend_page:
                with Timer("Write legend-only page", enabled=timing):
                    make_legend_page(pdf, type_colors, ncol=args.legend_cols)

            # Helper: embed legend unless legend-page mode
            def maybe_embed_legend(fig):
                if not args.legend_page:
                    add_embedded_legend(fig, type_colors, ncol=args.embedded_legend_cols, fontsize=7)

            if args.stacked:
                # ---------- Precompute shared axes limits ----------
                with Timer("STACKED precompute: K2P cache + global x/y limits", enabled=timing):
                    k2p_cache = {}
                    global_k2p_xmax = 0.0
                    global_k2p_ymax = 0.0
                    for sp in args.species:
                        with Timer(f"[{sp}] compute_species_k2p_stack(initial)", enabled=timing):
                            kdat = compute_species_k2p_stack(species_data[sp]["records"], k2p_xmax=args.k2p_xmax)
                        k2p_cache[sp] = kdat
                        if kdat is None:
                            continue
                        xgrid, dens_global, _ = kdat
                        global_k2p_xmax = max(global_k2p_xmax, float(xgrid[-1]))
                        global_k2p_ymax = max(global_k2p_ymax, float(np.nanmax(dens_global)))

                    # If not fixed, recompute each species to use global xmax so x-axes match
                    if args.k2p_xmax is None and global_k2p_xmax > 0:
                        for sp in args.species:
                            with Timer(f"[{sp}] compute_species_k2p_stack(recompute global xmax={global_k2p_xmax:.4g})", enabled=timing):
                                k2p_cache[sp] = compute_species_k2p_stack(species_data[sp]["records"], k2p_xmax=global_k2p_xmax)

                with Timer("STACKED precompute: global chromosome x max", enabled=timing):
                    global_chr_xmax = 0
                    for sp in args.species:
                        c2l = species_data[sp]["chrom2len"]
                        if c2l:
                            global_chr_xmax = max(global_chr_xmax, max(c2l.values()))
                    tlog(f"Global chromosome max length (bp): {global_chr_xmax}", enabled=timing)

                # Size shared bins and y-max
                bin_size = int(args.bin)

                def global_bins_for_key(key):
                    vmax = 0.0
                    for sp in args.species:
                        vals = []
                        for r in species_data[sp]["records"]:
                            v = r.get(key, None)
                            if v is None or not np.isfinite(v) or v < 0:
                                continue
                            vals.append(float(v))
                        if vals:
                            vmax = max(vmax, float(np.nanmax(vals)))
                    vmax = max(vmax, float(bin_size))
                    return np.arange(0, vmax + bin_size, bin_size)

                with Timer("STACKED precompute: global bins for size plots", enabled=timing):
                    bins_full = global_bins_for_key("full_len")
                    bins_ltr = global_bins_for_key("ltr_len")
                    bins_int = global_bins_for_key("internal_len")

                size_cache = {"full_len": {}, "ltr_len": {}, "internal_len": {}}
                y_full = y_ltr = y_int = 0.0

                def stacked_ymax(counts_by_type):
                    if not counts_by_type:
                        return 0.0
                    mat = np.vstack(list(counts_by_type.values()))
                    return float(np.sum(mat, axis=0).max()) if mat.size else 0.0

                with Timer("STACKED precompute: hist counts per species + y-max", enabled=timing):
                    for sp in args.species:
                        recs = species_data[sp]["records"]

                        with Timer(f"[{sp}] hist full_len", enabled=timing):
                            cb = compute_hist_counts(compute_size_by_type(recs, "full_len"), bins_full)
                        size_cache["full_len"][sp] = cb
                        y_full = max(y_full, stacked_ymax(cb))

                        with Timer(f"[{sp}] hist ltr_len", enabled=timing):
                            cb = compute_hist_counts(compute_size_by_type(recs, "ltr_len"), bins_ltr)
                        size_cache["ltr_len"][sp] = cb
                        y_ltr = max(y_ltr, stacked_ymax(cb))

                        with Timer(f"[{sp}] hist internal_len", enabled=timing):
                            cb = compute_hist_counts(compute_size_by_type(recs, "internal_len"), bins_int)
                        size_cache["internal_len"][sp] = cb
                        y_int = max(y_int, stacked_ymax(cb))

                # ---------- Page builder: multi-panel ----------
                def save_multipanel_page(title, plotter, xlabel=None, ylabel=None, embed_legend=True):
                    with Timer(f"Render+save multipanel page: {title}", enabled=timing):
                        r, c = panel_grid(len(args.species))
                        fig, axes = plt.subplots(r, c, figsize=(11, 8.5), squeeze=False)
                        axes_flat = axes.ravel()

                        for i, sp in enumerate(args.species):
                            ax = axes_flat[i]
                            plotter(ax, sp)

                            if ylabel is not None:
                                if i % c == 0:
                                    ax.set_ylabel(ylabel, fontsize=9)
                                else:
                                    ax.set_ylabel("")
                            if xlabel is not None:
                                if i >= (r - 1) * c:
                                    ax.set_xlabel(xlabel, fontsize=9)
                                else:
                                    ax.set_xlabel("")

                        for j in range(len(args.species), len(axes_flat)):
                            axes_flat[j].set_axis_off()

                        fig.suptitle(title, fontsize=12)

                        if embed_legend and (not args.legend_page):
                            add_embedded_legend(fig, type_colors, ncol=args.embedded_legend_cols, fontsize=7)
                            fig.tight_layout(rect=[0, 0, 0.72, 0.96])
                        else:
                            fig.tight_layout(rect=[0, 0, 1, 0.96])

                        pdf.savefig(fig)
                        plt.close(fig)

                # K2P multipanel
                def k2p_plotter(ax, sp):
                    plot_k2p_panel(
                        ax=ax,
                        sp_label=sp,
                        k2p_data=k2p_cache[sp],
                        type_colors=type_colors,
                        y_max=global_k2p_ymax * 1.02 if global_k2p_ymax > 0 else None,
                        show_global_line=args.density_global_line,
                    )

                save_multipanel_page(
                    title="LTR-RT K2P density (shared axes; area partitioned by type)",
                    plotter=k2p_plotter,
                    xlabel="K2P",
                    ylabel="Density",
                    embed_legend=True,
                )

                # Chromosome multipanel
                dot_note = ("  |  dots: " + ", ".join(args.dot_type)) if args.dot_type else ""
                def chr_plotter(ax, sp):
                    plot_chrom_panel(
                        ax=ax,
                        sp_label=sp,
                        chrom2len=species_data[sp]["chrom2len"],
                        records=species_data[sp]["records"],
                        type_colors=type_colors,
                        dot_types=args.dot_type,
                        x_max=global_chr_xmax if global_chr_xmax > 0 else None,
                    )

                save_multipanel_page(
                    title=f"Chromosome distribution of LTR-RTs (shared x-axis){dot_note}",
                    plotter=chr_plotter,
                    xlabel="bp",
                    ylabel="Chromosomes",
                    embed_legend=True,
                )

                # Density-style chromosome distribution (one page per species; faceted by chromosome)
                for sp in args.species:
                    plot_chrom_density_facets(
                        pdf=pdf,
                        sp=sp,
                        chrom2len=species_data[sp]["chrom2len"],
                        records=species_data[sp]["records"],
                        type_colors=type_colors,
                        window_bp=args.dens_window,
                        step_bp=args.dens_step,
                        threshold_percent=args.dens_threshold,
                        top_n=args.dens_top_types,
                        legend_page=args.legend_page,
                        embedded_legend_cols=args.embedded_legend_cols,
                        timing=timing,
                    )

                    plot_chrom_k2p_mean_facets(
                            pdf=pdf,
                            sp=sp,
                            chrom2len=species_data[sp]["chrom2len"],
                            records=species_data[sp]["records"],
                            window_bp=args.dens_window,
                            step_bp=args.dens_step,
                            legend_page=args.legend_page,
                            timing=timing,
                        )

                # Size pages
                def make_size_page(key, bins, y_max, title, xlabel):
                    def size_plotter(ax, sp):
                        plot_size_panel(
                            ax=ax,
                            sp_label=sp,
                            counts_by_type=size_cache[key][sp],
                            bins=bins,
                            type_colors=type_colors,
                            y_max=y_max * 1.05 if y_max > 0 else None,
                        )
                    save_multipanel_page(
                        title=title,
                        plotter=size_plotter,
                        xlabel=xlabel,
                        ylabel="Count",
                        embed_legend=True,
                    )

                make_size_page("full_len", bins_full, y_full, f"Full-length LTR-RT size distribution (bin={bin_size} bp; shared axes)", "Full length (bp)")
                make_size_page("ltr_len", bins_ltr, y_ltr, f"LTR size distribution (bin={bin_size} bp; shared axes)", "LTR length (bp)")
                make_size_page("internal_len", bins_int, y_int, f"Internal size distribution (bin={bin_size} bp; shared axes)", "Internal length (bp)")

            else:
                # Non-stacked: 5 pages/species
                for sp in args.species:
                    chrom2len = species_data[sp]["chrom2len"]
                    records = species_data[sp]["records"]

                    # 1) K2P density
                    with Timer(f"[{sp}] page: K2P density compute+render+save", enabled=timing):
                        with Timer(f"[{sp}] compute_species_k2p_stack", enabled=timing):
                            kdat = compute_species_k2p_stack(records, k2p_xmax=args.k2p_xmax)

                        fig, ax = plt.subplots(figsize=(11, 6))
                        plot_k2p_panel(
                            ax=ax,
                            sp_label=f"{sp}: LTR-RT K2P density",
                            k2p_data=kdat,
                            type_colors=type_colors,
                            y_max=None,
                            show_global_line=args.density_global_line,
                        )
                        ax.set_xlabel("K2P divergence")
                        ax.set_ylabel("Density")
                        maybe_embed_legend(fig)
                        if not args.legend_page:
                            fig.tight_layout(rect=[0, 0, 0.72, 1])
                        else:
                            fig.tight_layout()
                        pdf.savefig(fig)
                        plt.close(fig)

                    # 2) Chromosome distribution
                    with Timer(f"[{sp}] page: Chromosome distribution render+save", enabled=timing):
                        fig_h = max(6, 0.35 * max(1, len(chrom2len)))
                        fig, ax = plt.subplots(figsize=(11, fig_h))
                        dot_note = (" | dots: " + ", ".join(args.dot_type)) if args.dot_type else ""
                        plot_chrom_panel(
                            ax=ax,
                            sp_label=f"{sp}: Chromosome distribution{dot_note}",
                            chrom2len=chrom2len,
                            records=records,
                            type_colors=type_colors,
                            dot_types=args.dot_type,
                            x_max=None,
                            merge_gap_bp=args.chr_merge_gap,
                            rasterize=args.chr_rasterize,
                        )

                        ax.set_xlabel("Position (bp)")
                        maybe_embed_legend(fig)
                        if not args.legend_page:
                            fig.tight_layout(rect=[0, 0, 0.72, 1])
                        else:
                            fig.tight_layout()
                        pdf.savefig(fig)
                        plt.close(fig)

                    # 2b) Chromosome distribution (density-style; faceted)
                    plot_chrom_density_facets(
                        pdf=pdf,
                        sp=sp,
                        chrom2len=chrom2len,
                        records=records,
                        type_colors=type_colors,
                        window_bp=args.dens_window,
                        step_bp=args.dens_step,
                        threshold_percent=args.dens_threshold,
                        top_n=args.dens_top_types,
                        legend_page=args.legend_page,
                        embedded_legend_cols=args.embedded_legend_cols,
                        timing=timing,
                    )

                    plot_chrom_k2p_mean_facets(
                            pdf=pdf,
                            sp=sp,
                            chrom2len=species_data[sp]["chrom2len"],
                            records=species_data[sp]["records"],
                            window_bp=args.dens_window,
                            step_bp=args.dens_step,
                            legend_page=args.legend_page,
                            timing=timing,
                        )


                    # 3/4/5) Size distributions
                    for key, title, xlabel in [
                        ("full_len", "Full-length LTR-RT size distribution", "Full length (bp)"),
                        ("ltr_len", "LTR size distribution", "LTR length (bp)"),
                        ("internal_len", "Internal size distribution", "Internal length (bp)"),
                    ]:
                        with Timer(f"[{sp}] page: {key} size dist compute+render+save", enabled=timing):
                            with Timer(f"[{sp}] compute_size_by_type({key})", enabled=timing):
                                vb = compute_size_by_type(records, key)

                            # per-species bins in non-stacked mode
                            with Timer(f"[{sp}] bins + hist({key})", enabled=timing):
                                all_vals = []
                                for vals in vb.values():
                                    all_vals.extend([v for v in vals if np.isfinite(v) and v >= 0])
                                vmax = max(all_vals) if all_vals else float(args.bin)
                                bins = np.arange(0, max(vmax, float(args.bin)) + args.bin, args.bin)
                                counts_by_type = compute_hist_counts(vb, bins)

                            fig, ax = plt.subplots(figsize=(11, 6))
                            plot_size_panel(
                                ax=ax,
                                sp_label=f"{sp}: {title} (bin={args.bin} bp)",
                                counts_by_type=counts_by_type,
                                bins=bins,
                                type_colors=type_colors,
                                y_max=None,
                            )
                            ax.set_xlabel(xlabel)
                            ax.set_ylabel("Count")
                            maybe_embed_legend(fig)
                            if not args.legend_page:
                                fig.tight_layout(rect=[0, 0, 0.72, 1])
                            else:
                                fig.tight_layout()
                            pdf.savefig(fig)
                            plt.close(fig)

    print(f"Wrote: {args.outpdf}")


if __name__ == "__main__":
    main()
