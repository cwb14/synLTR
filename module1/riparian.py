#!/usr/bin/env python3
"""Syntenic ribbon (riparian) plot from pairwise anchor coords + FAI indexes.

Genomes are stacked as horizontal tracks in a user-given order. Ribbons connect
ADJACENT tracks only (the riparian convention), and every ribbon is coloured by
the chromosome of the FIRST genome in the stack -- the colour is carried down the
stack layer by layer, so a colour tracks one ancestral chromosome all the way to
the bottom even where it fissions or fuses.

Inputs
  coords : 3-column TSV, one syntenic block per line
             seq:start..end <TAB> seq:start..end <TAB> strand
  fai    : samtools faidx index; column 1 = seq name, column 2 = length

Outputs: <prefix>.pdf, <prefix>.png, <prefix>.html (interactive zoom/pan).

python riparian.py --coords Ahall.Aare.anchors.coords Athal.Aare.anchors.coords Athal.Ahall.anchors.coords --fai Athal.fa.fai Ahall.fa.fai Aare.fa.fai --order Athal,Ahall,Aare -o riparian --scale bp
"""

import argparse
import html
import os
import re
import sys
from collections import OrderedDict

try:
    import numpy as np
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patheffects as pe
    from matplotlib.path import Path
    from matplotlib.patches import PathPatch, Rectangle
except ImportError as exc:
    sys.exit(
        f"ERROR: missing dependency ({exc.name}).\n"
        f"       riparian.py needs numpy + matplotlib, and is running under\n"
        f"       {sys.executable}\n"
        f"       which does not have them. Either activate an env that does:\n"
        f"           mamba activate synteny\n"
        f"       or call that interpreter directly:\n"
        f"           /home/chris/bin/mambaforge/envs/synteny/bin/python riparian.py ...\n"
        f"       (a fresh env can be built from the environment.yml beside this script)"
    )


# -----------------------------
# Logging
# -----------------------------

def tlog(msg, enabled=True):
    if enabled:
        sys.stderr.write(msg + "\n")


def die(msg):
    sys.stderr.write("ERROR: " + msg + "\n")
    sys.exit(1)


# -----------------------------
# Parsing
# -----------------------------

def read_fai(path):
    """name -> length, in file order."""
    d = OrderedDict()
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 2:
                continue
            d[f[0]] = int(f[1])
    return d


_NAT_RE = re.compile(r"(\d+)")


def sort_chrom_names(names):
    """Natural sort: chr2 before chr10."""
    def key(n):
        return [int(p) if p.isdigit() else p.lower() for p in _NAT_RE.split(n)]
    return sorted(names, key=key)


_LOCUS_RE = re.compile(r"^(.+):(\d+)\.\.(\d+)$")


def _parse_locus(tok):
    m = _LOCUS_RE.match(tok)
    if not m:
        return None
    name, s, e = m.group(1), int(m.group(2)), int(m.group(3))
    if s > e:
        s, e = e, s
    return name, s, e


def parse_coords(path, seq2genome):
    """Read a coords file into blocks tagged with their genome of origin.

    Genome identity comes from the sequence names (looked up in the FAI union),
    not from the filename or column order, so a mislabelled file cannot silently
    transpose the plot.
    """
    blocks = []
    bad = 0
    with open(path) as fh:
        for ln, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 2:
                bad += 1
                continue
            a, b = _parse_locus(f[0]), _parse_locus(f[1])
            if a is None or b is None:
                bad += 1
                continue
            strand = f[2].strip() if len(f) > 2 and f[2].strip() in ("+", "-") else "+"
            ga, gb = seq2genome.get(a[0]), seq2genome.get(b[0])
            if ga is None or gb is None:
                bad += 1
                continue
            blocks.append({"ga": ga, "ca": a[0], "sa": a[1], "ea": a[2],
                           "gb": gb, "cb": b[0], "sb": b[1], "eb": b[2],
                           "strand": strand})
    if bad:
        tlog(f"  [{os.path.basename(path)}] skipped {bad} unparseable/unknown-seq line(s)")
    if not blocks:
        die(f"no usable blocks in {path} (are the FAI files for the same assemblies?)")
    genomes = {blk["ga"] for blk in blocks} | {blk["gb"] for blk in blocks}
    if len(genomes) != 2:
        die(f"{path}: expected exactly 2 genomes, found {sorted(genomes)}")
    return blocks


# -----------------------------
# Chromosome detection
# -----------------------------
#
# A mixed assembly holds two populations: chromosomes (few, long, usually sharing
# a header convention) and everything else (many, short, a different convention).
# Length and header turn out to be the same cliff-detection problem at two
# granularities, so a single routine covers both:
#
#   header granularity  group sequences by header signature (digit runs collapsed
#                       to '#'), then cut between groups. CM023363.1 and
#                       CM023401.1 share 'CM#.#'; JAAQRD020000065.1 does not.
#   length granularity  when headers carry no class signal they collapse to one
#                       group; every sequence then becomes its own group and the
#                       identical code path degenerates to a pure length cut.
#
# Headers are tried first because they are the only signal that survives size
# classes overlapping -- an unplaced scaffold can be longer than a microchromosome,
# so the chromosome set is not always a length-rank prefix. Length is the fallback
# for assemblies whose headers say nothing (e.g. chromosomes still named
# 'scaffold_N').
#
# Cuts are RANKED (not argmax'd) by Otsu between-class variance on log10(length),
# and the best cut that survives the guards wins. Ranking is what makes this
# robust: a spurious cliff gets skipped rather than obeyed.

_DIGIT_RUN_RE = re.compile(r"\d+")

# Guard defaults, tuned against a sweep of synthetic assembly shapes.
CHROM_MAX_COUNT = 200     # chromosomes are FEW
CHROM_MIN_FRAC = 0.30     # chromosomes ARE the genome (rejects fragmented assemblies)
CHROM_MIN_FOLD = 2.0      # the size drop must be a cliff, not a gradient
CHROM_PURE_RATIO = 0.10   # no junk class at all: everything is chromosome-scale


def header_signature(name):
    """Collapse digit runs, keeping the characters that mark the sequence class."""
    return _DIGIT_RUN_RE.sub("#", name)


def _rank_cuts(medians, counts):
    """Rank cut points by Otsu between-class variance on log10(median length).

    Cut j assigns groups[:j+1] to the chromosome class. Scoring on a log scale
    makes this scale-free; weighting by group size keeps a 500-member scaffold
    class from being outvoted by a 1-member chromosome class. Returns
    [(score, j), ...] best first.
    """
    x = np.log10(np.maximum(np.asarray(medians, dtype=float), 1.0))
    w = np.asarray(counts, dtype=float)
    cw = np.cumsum(w)
    cx = np.cumsum(w * x)
    total_w, total_x = cw[-1], cx[-1]

    scored = []
    for j in range(len(medians) - 1):
        w1, w2 = cw[j], total_w - cw[j]
        if w1 <= 0 or w2 <= 0:
            continue
        m1, m2 = cx[j] / w1, (total_x - cx[j]) / w2
        scored.append(((w1 * w2 / (total_w * total_w)) * (m1 - m2) ** 2, j))

    scored.sort(reverse=True)
    return scored


def _best_partition(name2len, by_header, max_chroms, min_frac, min_fold):
    """Best guard-passing chromosome set at one granularity, else None."""
    buckets = OrderedDict()
    for name in name2len:
        key = header_signature(name) if by_header else name
        buckets.setdefault(key, []).append(name)

    groups = sorted(
        ((names, float(np.median([name2len[n] for n in names]))) for names in buckets.values()),
        key=lambda g: -g[1],
    )
    if len(groups) < 2:
        return None

    medians = [m for _, m in groups]
    sizes = [len(names) for names, _ in groups]
    cum_count = np.cumsum(sizes)
    cum_mass = np.cumsum([sum(name2len[n] for n in names) for names, _ in groups])
    total = float(sum(name2len.values()))

    for _, j in _rank_cuts(medians, sizes):
        fold = medians[j] / medians[j + 1] if medians[j + 1] > 0 else float("inf")
        if fold < min_fold:
            continue
        if not (2 <= cum_count[j] <= max_chroms):
            continue
        mass = cum_mass[j] / total
        if mass < min_frac:
            continue
        return {
            "chroms": [n for names, _ in groups[: j + 1] for n in names],
            "fold": float(fold),
            "mass": float(mass),
            "n_groups": len(groups),
        }
    return None


def detect_chromosomes(
    name2len,
    max_chroms=CHROM_MAX_COUNT,
    min_frac=CHROM_MIN_FRAC,
    min_fold=CHROM_MIN_FOLD,
    pure_ratio=CHROM_PURE_RATIO,
):
    """Split an assembly into chromosomes and everything else.

    Returns (chrom_names, info). An empty chrom_names means the assembly is not
    chromosome-level, and coordinate-based plots should be skipped.
    """
    info = {"rule": "not chromosome-level", "n_seqs": len(name2len),
            "fold": None, "mass": None, "agree": None, "n_groups": None}
    if len(name2len) < 2:
        return [], info

    lengths = sorted(name2len.values(), reverse=True)

    # No junk class at all: every sequence is within pure_ratio of the largest.
    if lengths[-1] >= pure_ratio * lengths[0] and len(lengths) <= max_chroms:
        info.update(rule="pure-scale", mass=1.0)
        return list(name2len), info

    by_header = _best_partition(name2len, True, max_chroms, min_frac, min_fold)
    by_length = _best_partition(name2len, False, max_chroms, min_frac, min_fold)

    chosen = by_header or by_length
    if chosen is None:
        return [], info

    if by_header and by_length:
        # Independent corroboration when both granularities land on the same set.
        info["agree"] = set(by_header["chroms"]) == set(by_length["chroms"])

    info.update(
        rule="header-class" if by_header else "length-cliff",
        fold=chosen["fold"],
        mass=chosen["mass"],
        n_groups=chosen["n_groups"],
    )
    return chosen["chroms"], info


def select_chromosomes(fai_path, chrom_regex=None, min_chrom_len=None):
    """Read a FAI and decide which of its sequences are chromosomes.

    Returns (chrom2len, all2len, info). chrom2len is empty when the assembly is
    not chromosome-level. Explicit overrides bypass detection entirely.
    """
    all2len = read_fai(fai_path)
    if not all2len:
        return OrderedDict(), all2len, {"rule": "empty FAI", "n_seqs": 0,
                                        "fold": None, "mass": None, "agree": None}

    if chrom_regex:
        rx = re.compile(chrom_regex)
        chroms = [n for n in all2len if rx.search(n)]
        info = {"rule": f"--chrom-regex {chrom_regex!r}", "n_seqs": len(all2len),
                "fold": None, "mass": None, "agree": None}
    elif min_chrom_len:
        chroms = [n for n, L in all2len.items() if L >= int(min_chrom_len)]
        info = {"rule": f"--min-chrom-len {int(min_chrom_len)}", "n_seqs": len(all2len),
                "fold": None, "mass": None, "agree": None}
    else:
        chroms, info = detect_chromosomes(all2len)

    if chroms:
        total = float(sum(all2len.values()))
        info["mass"] = sum(all2len[n] for n in chroms) / total if total else 0.0

    chrom2len = OrderedDict((n, all2len[n]) for n in sort_chrom_names(chroms))
    return chrom2len, all2len, info


def report_chromosome_call(sp, chrom2len, all2len, info):
    """Always say out loud what was decided -- a misfire must be visible, not silent."""
    n_all = len(all2len)
    if not chrom2len:
        tlog(f"[{sp}] NOT chromosome-level ({n_all} sequences, rule: {info['rule']}). "
             f"Skipping chromosome-based plots.")
        return

    bits = [f"{len(chrom2len)}/{n_all} sequences", f"rule: {info['rule']}"]
    if info.get("mass") is not None:
        bits.append(f"{info['mass']:.1%} of assembly bp")
    if info.get("fold"):
        bits.append(f"{info['fold']:,.1f}x size cliff")
    if info.get("agree") is True:
        bits.append("header + length agree")
    elif info.get("agree") is False:
        bits.append("WARNING: header and length disagree -- check --chrom-regex")

    tlog(f"[{sp}] chromosome-level: " + "  |  ".join(bits))
    kept = list(chrom2len)
    preview = (", ".join(kept) if len(kept) <= 5
               else ", ".join(kept[:4]) + f", ... , {kept[-1]}")
    tlog(f"[{sp}]   chromosomes: {preview}")


# -----------------------------
# Colours (colourblind-safe)
# -----------------------------

# Okabe-Ito, minus black (unusable as a translucent ribbon fill).
_OKABE_ITO = ["#0072B2", "#E69F00", "#009E73", "#D55E00",
              "#CC79A7", "#56B4E9", "#F0E442", "#999999"]
UNPAINTED = "#D9D9D9"   # ribbon whose ancestry could not be traced


def build_palette(names):
    if len(names) <= len(_OKABE_ITO):
        cols = _OKABE_ITO[:len(names)]
    else:
        cmap = plt.get_cmap("tab20")
        cols = [matplotlib.colors.to_hex(cmap(i % 20)) for i in range(len(names))]
    return OrderedDict(zip(names, cols))


# -----------------------------
# Layout
# -----------------------------

class Track:
    """One genome laid out along x, with per-chromosome offsets and flips.

    Span is independent of chromosome order and of flips, so `denom` and `shift`
    can be fixed once, up front, and stay valid while the optimiser reshuffles.
    """

    def __init__(self, genome, chrom2len):
        self.genome = genome
        self.lens = chrom2len
        self.order = list(chrom2len)          # left-to-right chromosome order
        self.flip = {c: False for c in self.order}
        self.offset = {}                      # bp offset of each chromosome's left edge
        self.shift = 0.0                      # centring shift, bp
        self.span = 0.0
        self.denom = 1.0                      # bp that map onto one unit of x

    def relayout(self, gap_bp):
        pos = 0.0
        for c in self.order:
            self.offset[c] = pos
            pos += self.lens[c] + gap_bp
        self.span = max(pos - gap_bp, 0.0)

    def x(self, chrom, coord):
        """bp coordinate on a chromosome -> bp position along this track."""
        L = self.lens[chrom]
        c = min(max(coord, 0), L)
        local = (L - c) if self.flip[chrom] else c
        return self.offset[chrom] + local + self.shift

    def xn(self, chrom, coord):
        """...and the same position normalised into the shared [0, 1] plot space."""
        return self.x(chrom, coord) / self.denom


def _weighted_cov(xs, ys, ws):
    xs, ys, ws = np.asarray(xs, float), np.asarray(ys, float), np.asarray(ws, float)
    W = ws.sum()
    if W <= 0:
        return 0.0
    mx, my = (ws * xs).sum() / W, (ws * ys).sum() / W
    return float((ws * (xs - mx) * (ys - my)).sum() / W)


def optimise_track(track, prev_track, blocks, do_order=True, do_flip=True):
    """Order and orient one track against the track already placed above it.

    Order: each chromosome goes to the length-weighted mean x of its partners on
    the track above, so ribbons run as vertically as possible.
    Flip: sign of the length-weighted covariance between the partner's *laid-out*
    x (which already absorbs that partner's own flip) and this chromosome's bp
    coordinate. Negative covariance means the chromosome reads backwards.
    """
    by_chrom = {}
    for b in blocks:
        by_chrom.setdefault(b["cb"], []).append(b)

    if do_flip:
        for c, bl in by_chrom.items():
            # Judge orientation against the single partner carrying the most bp:
            # pooling partners would mix unrelated x-offsets and blur the sign.
            mass = {}
            for b in bl:
                mass[b["ca"]] = mass.get(b["ca"], 0) + (b["eb"] - b["sb"])
            dom = max(mass, key=mass.get)
            xs, ys, ws = [], [], []
            for b in bl:
                if b["ca"] != dom:
                    continue
                w = b["eb"] - b["sb"]
                xs.append(prev_track.xn(b["ca"], (b["sa"] + b["ea"]) / 2.0))
                ys.append((b["sb"] + b["eb"]) / 2.0)
                ws.append(w)
            if len(xs) >= 2 and _weighted_cov(xs, ys, ws) < 0:
                track.flip[c] = True

    if do_order:
        anchor = {}
        for c, bl in by_chrom.items():
            xs = [prev_track.xn(b["ca"], (b["sa"] + b["ea"]) / 2.0) for b in bl]
            ws = [float(b["eb"] - b["sb"]) for b in bl]
            W = sum(ws)
            anchor[c] = sum(x * w for x, w in zip(xs, ws)) / W if W > 0 else None

        placed = [c for c in track.order if anchor.get(c) is not None]
        orphan = [c for c in track.order if anchor.get(c) is None]
        placed.sort(key=lambda c: anchor[c])
        track.order = placed + orphan   # unanchored chromosomes park on the right


def uncovered_frac(cov, chrom, s, e):
    """Fraction of [s, e] on `chrom` not already covered by intervals in `cov`."""
    span = e - s
    if span <= 0:
        return 0.0
    clipped = sorted((max(s, a), min(e, b))
                     for a, b in cov.get(chrom, ()) if min(e, b) > max(s, a))
    if not clipped:
        return 1.0

    covered, cur_s, cur_e = 0, clipped[0][0], clipped[0][1]
    for a, b in clipped[1:]:
        if a > cur_e:
            covered += cur_e - cur_s
            cur_s, cur_e = a, b
        else:
            cur_e = max(cur_e, b)
    covered += cur_e - cur_s
    return max(0.0, 1.0 - covered / span)


def triage_blocks(blocks, max_len_ratio, rescue_uncovered, rule="both"):
    """Split blocks into clean, rescued and dropped.

    A length-skewed block is not automatically wrong -- it is only *redundant* if
    clean collinear blocks already explain the same territory. So skew alone does
    not condemn a block; skew plus redundancy does. Anything skewed that lands
    where nothing clean reaches is the only evidence for that locus, and dropping
    it would punch a hole in the figure.

    Novelty must hold on BOTH sides by default. Scoring on either side alone does
    not discriminate: a degenerate block is skewed precisely because one side is
    tiny, and a tiny window almost always falls in some gap between clean blocks,
    so its "novelty" is near 1 no matter how redundant the long side is. On this
    data `either` rescued 148/161 degenerate blocks, the 206x worst offender
    included; `both` rescues 21 and drops the rest. A rescued block should fill a
    genuine two-sided hole, not duplicate territory one track already explains.

    Degenerate blocks are considered least-skewed first, and each rescued block
    joins the coverage, so a single hole attracts one rescuer instead of a pile.
    """
    score = min if rule == "both" else max
    clean, degen = [], []
    for b in blocks:
        la, lb = b["ea"] - b["sa"], b["eb"] - b["sb"]
        ratio = max(la, lb) / min(la, lb) if min(la, lb) > 0 else float("inf")
        (degen if max_len_ratio > 0 and ratio > max_len_ratio else clean).append((ratio, b))

    cov_up, cov_dn = {}, {}
    for _, b in clean:
        cov_up.setdefault(b["ca"], []).append((b["sa"], b["ea"]))
        cov_dn.setdefault(b["cb"], []).append((b["sb"], b["eb"]))

    rescued, dropped = [], []
    for ratio, b in sorted(degen, key=lambda rb: rb[0]):
        u_up = uncovered_frac(cov_up, b["ca"], b["sa"], b["ea"])
        u_dn = uncovered_frac(cov_dn, b["cb"], b["sb"], b["eb"])
        novelty = score(u_up, u_dn)
        if rescue_uncovered <= novelty:
            b = dict(b, rescued=True, ratio=ratio, novelty=novelty)
            rescued.append(b)
            cov_up.setdefault(b["ca"], []).append((b["sa"], b["ea"]))
            cov_dn.setdefault(b["cb"], []).append((b["sb"], b["eb"]))
        else:
            dropped.append((ratio, novelty, b))

    return [b for _, b in clean], rescued, dropped


def overlap_paint(paint_chrom, s, e):
    """Trace [s, e] back to a reference chromosome.

    Prefer the greatest bp overlap. Failing that, fall back to the nearest painted
    segment on the same chromosome: filtering upstream leaves gaps in the paint,
    and a block sitting in one of those gaps still descends from whatever flanks
    it. Nearest (not chromosome-majority) so that a fission chromosome keeps the
    right ancestry on each of its arms. Returns None only for a chromosome with
    no paint at all.
    """
    best, best_ov = None, 0
    for (ps, pe, ref) in paint_chrom:
        ov = min(e, pe) - max(s, ps)
        if ov > best_ov:
            best, best_ov = ref, ov
    if best is not None:
        return best

    best_d = None
    for (ps, pe, ref) in paint_chrom:
        d = ps - e if ps > e else (s - pe if s > pe else 0)
        if best_d is None or d < best_d:
            best_d, best = d, ref
    return best


# -----------------------------
# Geometry -> renderer-neutral primitives
# -----------------------------

def bezier_ribbon(a1, a2, b1, b2, y_top, y_bot, slack=0.5):
    """Corner points + control points for one ribbon.

    (a1,a2) are the block's endpoints on the upper track, (b1,b2) the matching
    endpoints below -- already strand-swapped by the caller, so a '-' block
    yields a twisted (bow-tie) ribbon, which is exactly how an inversion should
    read.
    """
    dy = (y_bot - y_top) * slack
    return {
        "a1": a1, "a2": a2, "b1": b1, "b2": b2,
        "y_top": y_top, "y_bot": y_bot,
        "cy1": y_top + dy, "cy2": y_bot - dy,
    }


def paint_segments(intervals, length, nbins=2000):
    """Collapse overlapping painted intervals into contiguous segments over [0, L].

    The paint arrives as one interval per syntenic block, so it overlaps itself,
    leaves gaps, and disagrees with itself where two reference chromosomes both
    hit the same place. Binning resolves all three at once: each bin takes the
    reference with the most bp in it, unpainted bins inherit the nearest painted
    one (matching overlap_paint's fallback), and runs of equal bins merge.

    Returns [(start_bp, end_bp, ref_chrom), ...] tiling the whole chromosome.
    """
    if not intervals or length <= 0:
        return []

    refs = sorted({r for _, _, r in intervals})
    edges = np.linspace(0.0, float(length), nbins + 1)
    cov = np.zeros((len(refs), nbins))

    for k, ref in enumerate(refs):
        for s, e, r in intervals:
            if r != ref or e <= s:
                continue
            b0 = max(0, min(nbins - 1, int(s / length * nbins)))
            b1 = max(b0 + 1, min(nbins, int(np.ceil(e / length * nbins))))
            cov[k, b0:b1] += (np.minimum(e, edges[b0 + 1:b1 + 1])
                              - np.maximum(s, edges[b0:b1]))

    win = np.argmax(cov, axis=0)
    painted = cov.sum(axis=0) > 0
    if not painted.any():
        return []

    # Unpainted bins inherit whichever painted bin is closer.
    idx = np.flatnonzero(painted)
    bins = np.arange(nbins)
    pos = np.searchsorted(idx, bins)
    left = idx[np.clip(pos - 1, 0, len(idx) - 1)]
    right = idx[np.clip(pos, 0, len(idx) - 1)]
    nearest = np.where(np.abs(bins - left) <= np.abs(right - bins), left, right)
    win = np.where(painted, win, win[nearest])

    segs, start = [], 0
    for b in range(1, nbins + 1):
        if b == nbins or win[b] != win[start]:
            segs.append((edges[start], edges[b], refs[win[start]]))
            start = b
    return segs


def build_scene(order, tracks, pair_blocks, palette, ref_genome, scale, bp_per_x,
                bar_h=0.10, track_gap=1.0):
    """Everything both renderers need: bars, ribbons, labels -- in unit x, track y."""
    X = lambda t, chrom, coord: t.xn(chrom, coord)

    y_of = {g: i * track_gap for i, g in enumerate(order)}

    ribbons = []
    # paint[genome][chrom] = [(start, end, ref_chrom), ...] -- ancestry carried down
    paint = {order[0]: {c: [(0, tracks[order[0]].lens[c], c)]
                        for c in tracks[order[0]].order}}

    for i in range(len(order) - 1):
        up, dn = order[i], order[i + 1]
        tu, td = tracks[up], tracks[dn]
        y_top = y_of[up] + bar_h
        y_bot = y_of[dn]
        paint.setdefault(dn, {})

        for b in pair_blocks[(up, dn)]:
            ref = overlap_paint(paint[up].get(b["ca"], []), b["sa"], b["ea"])
            colour = palette.get(ref, UNPAINTED) if ref else UNPAINTED

            a1, a2 = X(tu, b["ca"], b["sa"]), X(tu, b["ca"], b["ea"])
            if b["strand"] == "+":
                b1, b2 = X(td, b["cb"], b["sb"]), X(td, b["cb"], b["eb"])
            else:
                b1, b2 = X(td, b["cb"], b["eb"]), X(td, b["cb"], b["sb"])

            r = bezier_ribbon(a1, a2, b1, b2, y_top, y_bot)
            r["color"] = colour
            r["strand"] = b["strand"]
            r["rescued"] = bool(b.get("rescued"))
            r["ratio"] = b.get("ratio")
            r["ref"] = ref or "unassigned"
            r["up"] = f'{b["ca"]}:{b["sa"]:,}-{b["ea"]:,}'
            r["dn"] = f'{b["cb"]}:{b["sb"]:,}-{b["eb"]:,}'
            r["len_up"] = b["ea"] - b["sa"]
            r["len_dn"] = b["eb"] - b["sb"]
            ribbons.append(r)

            if ref:
                paint[dn].setdefault(b["cb"], []).append((b["sb"], b["eb"], ref))

    # Draw widest ribbons first so slivers land on top and stay visible.
    ribbons.sort(key=lambda r: -(r["len_up"] + r["len_dn"]))

    # Bars come last: a non-reference chromosome is filled with the mosaic of the
    # reference chromosomes it descends from, so a fission reads straight off the
    # bar (Ahall_chr3 lands half green, half yellow) without tracing a ribbon.
    bars = []
    for g in order:
        t = tracks[g]
        for c in t.order:
            L = t.lens[c]
            x0, x1 = X(t, c, 0), X(t, c, L)

            if g == ref_genome:
                intervals = [(0, L, c)]
            else:
                intervals = paint.get(g, {}).get(c, [])
            segs = []
            for s, e, ref in paint_segments(intervals, L):
                sx0, sx1 = X(t, c, s), X(t, c, e)
                segs.append({"x0": min(sx0, sx1), "x1": max(sx0, sx1),
                             "color": palette.get(ref, UNPAINTED), "ref": ref})

            bars.append({
                "genome": g, "chrom": c,
                "x0": min(x0, x1), "x1": max(x0, x1),
                "y": y_of[g], "h": bar_h,
                "len": L, "flip": t.flip[c], "segs": segs,
            })

    return {
        "bars": bars, "ribbons": ribbons, "y_of": y_of,
        "scale": scale, "bp_per_x": bp_per_x, "bar_h": bar_h,
        "height": (len(order) - 1) * track_gap + bar_h,
        "order": order, "palette": palette, "ref_genome": ref_genome,
    }


def nice_scale_bar(bp_per_x):
    """A round Mb value spanning roughly a tenth of the plot width.

    Only meaningful under --scale bp; under --scale fit each genome has its own
    bp-per-pixel, so there is no single ruler to draw.
    """
    target = bp_per_x / 10.0
    for mb in (1, 2, 5, 10, 20, 25, 50, 100, 200, 500):
        if mb * 1e6 >= target:
            return mb * 1e6 / bp_per_x, f"{mb} Mb"
    return 1e9 / bp_per_x, "1 Gb"


# -----------------------------
# Renderer: matplotlib (PDF / PNG)
# -----------------------------

def render_mpl(scene, out_paths, width, height, dpi, alpha, linewidth):
    fig, ax = plt.subplots(figsize=(width, height))

    for r in scene["ribbons"]:
        verts = [
            (r["a1"], r["y_top"]),
            (r["a2"], r["y_top"]),
            (r["a2"], r["cy1"]), (r["b2"], r["cy2"]), (r["b2"], r["y_bot"]),
            (r["b1"], r["y_bot"]),
            (r["b1"], r["cy2"]), (r["a1"], r["cy1"]), (r["a1"], r["y_top"]),
            (r["a1"], r["y_top"]),
        ]
        codes = [Path.MOVETO, Path.LINETO,
                 Path.CURVE4, Path.CURVE4, Path.CURVE4,
                 Path.LINETO,
                 Path.CURVE4, Path.CURVE4, Path.CURVE4,
                 Path.CLOSEPOLY]
        # Rescued (length-skewed) blocks get a dashed outline: they are the only
        # evidence at their locus, but they are still degenerate, and the static
        # figure must not present them with the same confidence as a clean block.
        if r["rescued"]:
            ax.add_patch(PathPatch(Path(verts, codes), facecolor=r["color"],
                                   edgecolor=r["color"], alpha=alpha * 0.8,
                                   linewidth=0.5, linestyle=(0, (2, 1.5))))
        else:
            ax.add_patch(PathPatch(Path(verts, codes), facecolor=r["color"],
                                   edgecolor="none", alpha=alpha, linewidth=0))

    # Plain rectangles: x is in unit-genome space and y in track units, so any
    # corner rounding expressed in data coordinates comes out as an ellipse.
    for bar in scene["bars"]:
        w = max(bar["x1"] - bar["x0"], 1e-6)
        for seg in bar["segs"]:
            ax.add_patch(Rectangle(
                (seg["x0"], bar["y"]), max(seg["x1"] - seg["x0"], 1e-9), bar["h"],
                facecolor=seg["color"], edgecolor="none", zorder=3))
        # Outline drawn over the mosaic so the seams never break the chromosome border.
        ax.add_patch(Rectangle(
            (bar["x0"], bar["y"]), w, bar["h"],
            facecolor="none" if bar["segs"] else "#FFFFFF",
            edgecolor="#333333", linewidth=linewidth, zorder=3.5))

        label = re.sub(r"^.*?_", "", bar["chrom"])
        if bar["flip"]:
            label += "′"      # prime mark = reverse-complemented for display
        ax.text((bar["x0"] + bar["x1"]) / 2, bar["y"] + bar["h"] / 2, label,
                ha="center", va="center", fontsize=6.5, zorder=4, color="#000000",
                path_effects=[pe.withStroke(linewidth=1.8, foreground="white")])

    for g, y in scene["y_of"].items():
        ax.text(-0.012, y + scene["bar_h"] / 2, g, ha="right", va="center",
                fontsize=10, fontstyle="italic", fontweight="bold")

    # The y-axis is inverted, so va="top" is what makes the caption hang BELOW the
    # rule; va="bottom" would ride up and strike through it.
    y_sb = scene["height"] + 0.12
    if scene["scale"] == "bp":
        bw, bar_lab = nice_scale_bar(scene["bp_per_x"])
        ax.plot([0, bw], [y_sb, y_sb], color="#333333", lw=1.2, solid_capstyle="butt")
        ax.plot([0, 0], [y_sb - 0.015, y_sb + 0.015], color="#333333", lw=1.2)
        ax.plot([bw, bw], [y_sb - 0.015, y_sb + 0.015], color="#333333", lw=1.2)
        ax.text(bw / 2, y_sb + 0.045, bar_lab, ha="center", va="top", fontsize=7.5)
    else:
        ax.text(0, y_sb + 0.03, "genomes scaled to equal width — "
                                "lengths are not comparable between tracks",
                ha="left", va="top", fontsize=7, color="#777777")

    ax.set_xlim(-0.09, 1.02)
    ax.set_ylim(scene["height"] + 0.32, -0.10)   # inverted: first genome on top
    ax.axis("off")
    fig.tight_layout(pad=0.4)

    for p in out_paths:
        fig.savefig(p, dpi=dpi, bbox_inches="tight",
                    transparent=False, facecolor="white")
        tlog(f"  wrote {p}")
    plt.close(fig)


# -----------------------------
# Renderer: self-contained interactive SVG/HTML
# -----------------------------

_HTML_TMPL = """<!doctype html>
<meta charset="utf-8">
<title>{title}</title>
<style>
  html,body{{margin:0;height:100%;font:13px/1.45 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif;background:#fff;color:#222}}
  #wrap{{position:relative;width:100%;height:100%;overflow:hidden}}
  svg{{width:100%;height:100%;display:block;cursor:grab;touch-action:none}}
  svg.grabbing{{cursor:grabbing}}
  #hud{{position:absolute;top:10px;left:12px;background:rgba(255,255,255,.94);border:1px solid #ddd;border-radius:6px;padding:8px 10px;max-width:270px}}
  #hud b{{font-size:12px}}
  .key{{display:flex;align-items:center;gap:6px;margin:3px 0;cursor:pointer;user-select:none;font-size:11px;padding:1px 3px;border-radius:3px}}
  .key:hover{{background:#f0f0f0}}
  .key.off{{opacity:.32}}
  .sw{{width:13px;height:13px;border-radius:3px;border:1px solid #999;flex:none}}
  #tip{{position:absolute;pointer-events:none;background:rgba(20,20,20,.92);color:#fff;border-radius:5px;padding:6px 9px;font-size:11.5px;line-height:1.5;white-space:nowrap;opacity:0;transition:opacity .08s;z-index:9}}
  #tip .s{{color:#9fe6b0}}
  #tip .w{{color:#ffcf7a}}
  #hint{{position:absolute;bottom:10px;left:12px;font-size:11px;color:#888}}
  button{{font:11px inherit;margin-top:6px;padding:3px 9px;border:1px solid #ccc;background:#fafafa;border-radius:4px;cursor:pointer}}
  button:hover{{background:#f0f0f0}}
  .rib{{transition:opacity .1s}}
  .dim{{opacity:.04 !important}}
  .hot{{opacity:.95 !important;stroke:#111;stroke-width:.4}}
</style>
<div id="wrap">
<svg id="svg" viewBox="{vb}" preserveAspectRatio="xMidYMid meet">
  <g id="scene">
{body}
  </g>
</svg>
<div id="hud">
  <b>{ref} chromosomes</b>
  <div id="legend">
{legend}
  </div>
  <button id="reset">reset view</button>
</div>
<div id="tip"></div>
<div id="hint">scroll = zoom &nbsp;&middot;&nbsp; drag = pan &nbsp;&middot;&nbsp; click a colour to isolate</div>
</div>
<script>
(function(){{
  var svg=document.getElementById('svg'), tip=document.getElementById('tip'),
      wrap=document.getElementById('wrap'), HOME=[{vbx},{vby},{vbw},{vbh}], vb=HOME.slice();

  function apply(){{ svg.setAttribute('viewBox', vb.join(' ')); }}

  // Cursor-anchored zoom: the bp under the pointer must not move.
  svg.addEventListener('wheel', function(e){{
    e.preventDefault();
    var r=svg.getBoundingClientRect(),
        fx=(e.clientX-r.left)/r.width, fy=(e.clientY-r.top)/r.height,
        k=Math.exp(e.deltaY*0.0016),
        nw=Math.min(Math.max(vb[2]*k, HOME[2]/600), HOME[2]*1.6),
        s=nw/vb[2];
    vb[0]+=(vb[2]-nw)*fx; vb[1]+=(vb[3]-vb[3]*s)*fy; vb[2]=nw; vb[3]*=s;
    apply();
  }}, {{passive:false}});

  var drag=null;
  svg.addEventListener('pointerdown', function(e){{
    drag={{x:e.clientX, y:e.clientY, vx:vb[0], vy:vb[1]}};
    svg.setPointerCapture(e.pointerId); svg.classList.add('grabbing');
  }});
  svg.addEventListener('pointermove', function(e){{
    if(!drag) return;
    var r=svg.getBoundingClientRect();
    vb[0]=drag.vx-(e.clientX-drag.x)*vb[2]/r.width;
    vb[1]=drag.vy-(e.clientY-drag.y)*vb[3]/r.height;
    apply();
  }});
  function endDrag(){{ drag=null; svg.classList.remove('grabbing'); }}
  svg.addEventListener('pointerup', endDrag);
  svg.addEventListener('pointercancel', endDrag);

  document.getElementById('reset').onclick=function(){{ vb=HOME.slice(); apply(); }};

  var ribs=[].slice.call(document.querySelectorAll('.rib'));
  ribs.forEach(function(p){{
    p.addEventListener('mouseenter', function(){{
      p.classList.add('hot');
      tip.innerHTML = p.dataset.up + '<br>' + p.dataset.dn +
                      '<br><span class="s">strand ' + p.dataset.strand +
                      ' &middot; ' + p.dataset.size + ' &middot; ' + p.dataset.ref + '</span>' +
                      (p.dataset.note ? '<br><span class="w">' + p.dataset.note + '</span>' : '');
      tip.style.opacity=1;
    }});
    p.addEventListener('mousemove', function(e){{
      var r=wrap.getBoundingClientRect(), x=e.clientX-r.left+14, y=e.clientY-r.top+14;
      if(x+tip.offsetWidth > r.width) x -= tip.offsetWidth+26;
      if(y+tip.offsetHeight > r.height) y -= tip.offsetHeight+26;
      tip.style.left=x+'px'; tip.style.top=y+'px';
    }});
    p.addEventListener('mouseleave', function(){{
      p.classList.remove('hot'); tip.style.opacity=0;
    }});
  }});

  var solo=null;
  [].slice.call(document.querySelectorAll('.key')).forEach(function(k){{
    k.onclick=function(){{
      var c=k.dataset.chrom;
      solo = (solo===c) ? null : c;
      ribs.forEach(function(p){{
        p.classList.toggle('dim', solo!==null && p.dataset.ref!==solo);
      }});
      [].slice.call(document.querySelectorAll('.key')).forEach(function(o){{
        o.classList.toggle('off', solo!==null && o.dataset.chrom!==solo);
      }});
    }};
  }});
  apply();
}})();
</script>
"""


def _fmt_bp(n):
    if n >= 1e6:
        return f"{n / 1e6:.2f} Mb"
    if n >= 1e3:
        return f"{n / 1e3:.1f} kb"
    return f"{n} bp"


def render_html(scene, out_path, alpha, linewidth, px_w=1500, px_track=250):
    """Map unit-space geometry onto a pixel viewBox, then emit SVG + vanilla JS."""
    M_L, M_R, M_T, M_B = 105, 40, 55, 45
    inner_w = px_w - M_L - M_R
    px_h = M_T + M_B + int(scene["height"] * px_track)

    def PX(x):
        return M_L + x * inner_w

    def PY(y):
        return M_T + y * px_track

    o = []
    o.append('    <rect x="0" y="0" width="%d" height="%d" fill="#fff"/>' % (px_w, px_h))

    o.append('    <g id="ribbons">')
    for r in scene["ribbons"]:
        a1, a2 = PX(r["a1"]), PX(r["a2"])
        b1, b2 = PX(r["b1"]), PX(r["b2"])
        yt, yb = PY(r["y_top"]), PY(r["y_bot"])
        c1, c2 = PY(r["cy1"]), PY(r["cy2"])
        d = (f'M{a1:.2f},{yt:.2f} L{a2:.2f},{yt:.2f} '
             f'C{a2:.2f},{c1:.2f} {b2:.2f},{c2:.2f} {b2:.2f},{yb:.2f} '
             f'L{b1:.2f},{yb:.2f} '
             f'C{b1:.2f},{c2:.2f} {a1:.2f},{c1:.2f} {a1:.2f},{yt:.2f} Z')
        size = f'{_fmt_bp(r["len_up"])} / {_fmt_bp(r["len_dn"])}'
        if r["rescued"]:
            style = (f'fill-opacity="{alpha * 0.8}" stroke="{r["color"]}" '
                     f'stroke-width="0.8" stroke-dasharray="3 2"')
            note = f'rescued &#183; {r["ratio"]:.0f}x length-skew, only evidence here'
        else:
            style = f'fill-opacity="{alpha}" stroke="none"'
            note = ""
        o.append(
            f'      <path class="rib" d="{d}" fill="{r["color"]}" {style} '
            f'data-up="{html.escape(r["up"])}" data-dn="{html.escape(r["dn"])}" '
            f'data-strand="{r["strand"]}" data-size="{size}" '
            f'data-ref="{html.escape(r["ref"])}" data-note="{note}"/>'
        )
    o.append("    </g>")

    o.append('    <g id="bars">')
    for bar in scene["bars"]:
        x0, x1 = PX(bar["x0"]), PX(bar["x1"])
        y, h = PY(bar["y"]), bar["h"] * px_track
        w = max(x1 - x0, 0.6)

        # Square corners (not rx-rounded) so the mosaic segments meet the border
        # exactly, and so HTML matches the PDF/PNG pixel for pixel.
        for seg in bar["segs"]:
            sx0, sx1 = PX(seg["x0"]), PX(seg["x1"])
            o.append(
                f'      <rect x="{sx0:.2f}" y="{y:.2f}" width="{max(sx1 - sx0, 0.3):.2f}" '
                f'height="{h:.2f}" fill="{seg["color"]}" stroke="none"/>'
            )
        fill = "none" if bar["segs"] else "#FFFFFF"
        o.append(
            f'      <rect x="{x0:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" '
            f'fill="{fill}" stroke="#333" stroke-width="{linewidth}">'
            f'<title>{html.escape(bar["chrom"])} &#183; {bar["len"]:,} bp'
            f'{" &#183; displayed reverse-complemented" if bar["flip"] else ""}</title></rect>'
        )
        label = html.escape(re.sub(r"^.*?_", "", bar["chrom"]) + ("′" if bar["flip"] else ""))
        o.append(
            f'      <text x="{(x0 + x1) / 2:.2f}" y="{y + h / 2:.2f}" text-anchor="middle" '
            f'dominant-baseline="central" font-size="10" pointer-events="none" '
            f'fill="#000" stroke="#fff" stroke-width="2.6" paint-order="stroke">{label}</text>'
        )
    o.append("    </g>")

    o.append('    <g id="labels">')
    for g, y in scene["y_of"].items():
        o.append(
            f'      <text x="{M_L - 14}" y="{PY(y) + bar_px(scene, px_track) / 2:.2f}" '
            f'text-anchor="end" dominant-baseline="central" font-size="15" '
            f'font-style="italic" font-weight="600">{html.escape(g)}</text>'
        )
    sy = PY(scene["height"]) + 26
    if scene["scale"] == "bp":
        bw, bar_lab = nice_scale_bar(scene["bp_per_x"])
        sx0, sx1 = PX(0), PX(bw)
        o.append(f'      <line x1="{sx0:.1f}" y1="{sy}" x2="{sx1:.1f}" y2="{sy}" '
                 f'stroke="#333" stroke-width="1.4"/>')
        o.append(f'      <text x="{(sx0 + sx1) / 2:.1f}" y="{sy + 15}" text-anchor="middle" '
                 f'font-size="11" fill="#333">{bar_lab}</text>')
    else:
        o.append(f'      <text x="{PX(0):.1f}" y="{sy + 6}" font-size="11" fill="#888">'
                 f'genomes scaled to equal width &#8212; lengths are not comparable '
                 f'between tracks</text>')
    o.append("    </g>")

    legend = []
    for c, col in scene["palette"].items():
        legend.append(
            f'    <div class="key" data-chrom="{html.escape(c)}">'
            f'<span class="sw" style="background:{col}"></span>{html.escape(c)}</div>'
        )

    doc = _HTML_TMPL.format(
        title=f'Riparian: {" / ".join(scene["order"])}',
        vb=f"0 0 {px_w} {px_h}",
        vbx=0, vby=0, vbw=px_w, vbh=px_h,
        body="\n".join(o),
        legend="\n".join(legend),
        ref=html.escape(scene["ref_genome"]),
    )
    with open(out_path, "w") as fh:
        fh.write(doc)
    tlog(f"  wrote {out_path}")


def bar_px(scene, px_track):
    return scene["bar_h"] * px_track


# -----------------------------
# Main
# -----------------------------

def main(argv=None):
    p = argparse.ArgumentParser(
        description="Syntenic ribbon (riparian) plot from pairwise anchor coords + FAI.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--coords", nargs="+", required=True,
                   help="pairwise coords files: 'seq:start..end<TAB>seq:start..end<TAB>strand'")
    p.add_argument("--fai", nargs="+", required=True,
                   help="samtools faidx indexes, one per genome (e.g. Athal.fa.fai)")
    p.add_argument("--order", default=None,
                   help="comma-separated genome stack order, top to bottom. "
                        "The FIRST genome supplies the ribbon colours. "
                        "Default: order the --fai files were given.")
    p.add_argument("-o", "--out-prefix", default="riparian", help="output path prefix")
    p.add_argument("--formats", default="pdf,png,html", help="comma-separated: pdf,png,html")
    p.add_argument("--min-block-len", type=int, default=10000,
                   help="drop syntenic blocks shorter than this on either side (bp; 0 = keep all)")
    p.add_argument("--max-len-ratio", type=float, default=5.0,
                   help="a block whose two sides differ in length by more than this factor "
                        "is called degenerate: an anchor chain dragged out by one distant "
                        "outlier, spanning megabases on one genome and kilobases on the "
                        "other. Degenerate blocks are not dropped outright -- see "
                        "--rescue-uncovered. (0 = treat none as degenerate)")
    p.add_argument("--rescue-uncovered", type=float, default=0.5,
                   help="keep a degenerate block if at least this fraction of its span is "
                        "territory no clean block already covers. Skew alone does not condemn "
                        "a block -- skew plus redundancy does; a degenerate block over an "
                        "otherwise-empty locus is the only evidence there is. "
                        "(0 = keep every degenerate block; >1 = drop them all)")
    p.add_argument("--rescue-rule", choices=("both", "either"), default="both",
                   help="'both': a degenerate block must open new territory on BOTH sides to "
                        "be rescued. 'either': one side suffices -- but a degenerate block's "
                        "short side is nearly always novel, so this rescues almost everything "
                        "and defeats the filter. Use 'both' unless you know why not")
    p.add_argument("--chrom-regex", default=None,
                   help="override chromosome detection: keep sequences matching this regex")
    p.add_argument("--min-chrom-len", type=int, default=None,
                   help="override chromosome detection: keep sequences at least this long (bp)")
    p.add_argument("--scale", choices=("fit", "bp"), default="fit",
                   help="'fit': every genome spans the full width (the GENESPACE look; "
                        "ribbons run near-vertical, but lengths are NOT comparable between "
                        "tracks). 'bp': one shared bp scale, shorter genomes centred "
                        "(lengths comparable; a scale bar is drawn)")
    p.add_argument("--no-optimise", action="store_true",
                   help="keep chromosomes in FAI order instead of minimising ribbon crossings")
    p.add_argument("--no-flip", action="store_true",
                   help="never reverse-complement a chromosome for display")
    p.add_argument("--width", type=float, default=11.0, help="figure width (inches)")
    p.add_argument("--track-height", type=float, default=1.9,
                   help="height per genome track (inches)")
    p.add_argument("--dpi", type=int, default=300, help="raster DPI for PNG")
    p.add_argument("--alpha", type=float, default=0.62, help="ribbon opacity")
    p.add_argument("--linewidth", type=float, default=0.6, help="chromosome outline width")
    p.add_argument("--gap-frac", type=float, default=0.011,
                   help="gap between chromosomes, as a fraction of the widest genome")
    p.add_argument("-v", "--verbose", action="store_true", help="per-step progress")
    a = p.parse_args(argv)

    fmts = {f.strip().lower() for f in a.formats.split(",") if f.strip()}
    bad = fmts - {"pdf", "png", "html"}
    if bad:
        die(f"unknown format(s): {sorted(bad)}")

    # --- genomes from FAI ---
    tlog("Reading FAI indexes")
    genome_chroms, seq2genome, dupes = OrderedDict(), {}, []
    for path in a.fai:
        if not os.path.exists(path):
            die(f"missing FAI: {path}")
        g = os.path.basename(path).split(".")[0]
        chrom2len, all2len, info = select_chromosomes(path, a.chrom_regex, a.min_chrom_len)
        report_chromosome_call(g, chrom2len, all2len, info)
        if not chrom2len:
            die(f"{g} is not chromosome-level; cannot draw a riparian plot. "
                f"Force with --chrom-regex or --min-chrom-len.")
        genome_chroms[g] = chrom2len
        for name in all2len:
            if name in seq2genome and seq2genome[name] != g:
                dupes.append(name)
            seq2genome[name] = g
    if dupes:
        die(f"sequence names are not unique across genomes (e.g. {dupes[:3]}); "
            f"prefix them per genome so blocks can be attributed unambiguously.")

    # --- stack order ---
    order = ([s.strip() for s in a.order.split(",") if s.strip()]
             if a.order else list(genome_chroms))
    unknown = [g for g in order if g not in genome_chroms]
    if unknown:
        die(f"--order names unknown genome(s) {unknown}; known: {list(genome_chroms)}")
    if len(order) < 2:
        die("need at least 2 genomes to draw ribbons")
    ref_genome = order[0]
    tlog(f"Stack: {' -> '.join(order)}   (colours from {ref_genome})")

    # --- coords ---
    tlog("Reading coords")
    by_pair, used = {}, {}
    for path in a.coords:
        if not os.path.exists(path):
            die(f"missing coords: {path}")
        blocks = parse_coords(path, seq2genome)
        pair = frozenset((blocks[0]["ga"], blocks[0]["gb"]))
        by_pair[pair] = blocks
        used[pair] = path

    pair_blocks = {}
    for i in range(len(order) - 1):
        up, dn = order[i], order[i + 1]
        key = frozenset((up, dn))
        if key not in by_pair:
            die(f"no coords file for adjacent pair {up} / {dn}")
        survivors, n_raw, n_off, n_small = [], 0, 0, 0
        for b in by_pair[key]:
            n_raw += 1
            if b["ga"] == dn:      # transpose so 'a' is always the upper genome
                b = {"ga": b["gb"], "ca": b["cb"], "sa": b["sb"], "ea": b["eb"],
                     "gb": b["ga"], "cb": b["ca"], "sb": b["sa"], "eb": b["ea"],
                     "strand": b["strand"]}
            if b["ca"] not in genome_chroms[up] or b["cb"] not in genome_chroms[dn]:
                n_off += 1
                continue
            if min(b["ea"] - b["sa"], b["eb"] - b["sb"]) < a.min_block_len:
                n_small += 1
                continue
            survivors.append(b)

        clean, rescued, dropped = triage_blocks(survivors, a.max_len_ratio,
                                                a.rescue_uncovered, a.rescue_rule)
        kept = clean + rescued
        if not kept:
            die(f"no blocks survive filtering for {up} / {dn}")
        pair_blocks[(up, dn)] = kept

        tlog(f"  {up} <-> {dn}: {len(kept):,}/{n_raw:,} blocks kept "
             f"[{os.path.basename(used[key])}]")
        tlog(f"      {len(clean):,} clean, {len(rescued):,} rescued "
             f"(length-skewed > {a.max_len_ratio:g}x but filling territory no clean "
             f"block reaches), {len(dropped):,} dropped as redundant", a.verbose)
        tlog(f"      also dropped: {n_off:,} off-chromosome, {n_small:,} < "
             f"{a.min_block_len:,} bp", a.verbose)
        for b in sorted(rescued, key=lambda x: -x["ratio"])[:8]:
            tlog(f"        rescued {b['ratio']:6.1f}x  {b['novelty']:.0%} new  "
                 f"{b['ca']}:{b['sa']:,}-{b['ea']:,} <-> "
                 f"{b['cb']}:{b['sb']:,}-{b['eb']:,}", a.verbose)

    for pair, path in used.items():
        if not any(frozenset(k) == pair for k in pair_blocks):
            tlog(f"  note: {os.path.basename(path)} is unused "
                 f"({' / '.join(sorted(pair))} are not adjacent in --order)")

    # --- layout ---
    tlog("Laying out tracks", a.verbose)
    tracks = {g: Track(g, genome_chroms[g]) for g in order}
    totals = {g: sum(genome_chroms[g].values()) for g in order}
    widest = max(totals.values())

    # Under 'fit' each genome gets gaps proportional to its own length, so the
    # gaps stay visually equal once the genomes are stretched to a common width.
    for g in order:
        gap = a.gap_frac * (totals[g] if a.scale == "fit" else widest)
        tracks[g].relayout(gap)

    # Span is order- and flip-invariant, so normalisation can be pinned now and
    # stays correct through the optimiser below.
    max_span = max(t.span for t in tracks.values())
    for g in order:
        t = tracks[g]
        if a.scale == "fit":
            t.denom, t.shift = t.span, 0.0
        else:
            t.denom, t.shift = max_span, (max_span - t.span) / 2.0
    bp_per_x = max_span   # bp spanned by one unit of x (meaningful only for 'bp')

    for i in range(len(order) - 1):
        up, dn = order[i], order[i + 1]
        optimise_track(tracks[dn], tracks[up], pair_blocks[(up, dn)],
                       do_order=not a.no_optimise, do_flip=not a.no_flip)
        gap = a.gap_frac * (totals[dn] if a.scale == "fit" else widest)
        tracks[dn].relayout(gap)        # re-lay: order and flips have changed
        if a.verbose:
            flipped = [c for c in tracks[dn].order if tracks[dn].flip[c]]
            tlog(f"  {dn}: {' '.join(tracks[dn].order)}"
                 + (f"   (flipped: {', '.join(flipped)})" if flipped else ""))

    palette = build_palette(list(genome_chroms[ref_genome]))
    scene = build_scene(order, tracks, pair_blocks, palette, ref_genome,
                        a.scale, bp_per_x)

    n_grey = sum(1 for r in scene["ribbons"] if r["ref"] == "unassigned")
    tlog(f"Ribbons: {len(scene['ribbons']):,}"
         + (f"  ({n_grey:,} could not be traced to a {ref_genome} chromosome)" if n_grey else ""))

    # --- render ---
    tlog("Rendering")
    height = a.track_height * len(order)
    raster = [f"{a.out_prefix}.{f}" for f in ("pdf", "png") if f in fmts]
    if raster:
        render_mpl(scene, raster, a.width, height, a.dpi, a.alpha, a.linewidth)
    if "html" in fmts:
        render_html(scene, f"{a.out_prefix}.html", a.alpha, a.linewidth)

    tlog("Done")


if __name__ == "__main__":
    main()
