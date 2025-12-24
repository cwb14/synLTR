#!/usr/bin/env python3
"""
plot.py

python ltrharvest_plot_struct.py --scn <(awk '{if (NF>12) {for(i=1;i<=NF;i++) if(i!=NF-1) printf "%s%s",$i,(i==NF?"":" "); print ""} else print}' Bdact_ltr.ltrtools.stitched.scn) --gff  Bdact_ltr.work/Bdact_ltr.ltrtools.internals.fa.rexdb-plant.dom.gff3 --tsv Bdact_ltr_kmer2ltr_dedup --out_dir Bdact_struct --min_presence 0.5 --max_individual 100000000000

Usage:
  python plot.py --tsv ltrrt.tsv --scn file.scn --gff hmm.gff --out_dir dir_name

Outputs (PDF):
  <out_dir>/<safe_family>_average.pdf
  <out_dir>/<safe_family>_individual.pdf
  <out_dir>/all_elements.pdf

Key features:
- Strand-normalizes all elements for plotting/averaging (flips '-' into '+').
- Average plots EXCLUDE rare proteins unless present in >= --min_presence fraction of elements.
- Average plots include whiskers showing 95% CI for start/end of each feature (bootstrap CI).
- Prints per-family summary stats + 95% CI to terminal.
"""

import argparse
import os
import re
import random
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional

import matplotlib.pyplot as plt
from matplotlib.patches import Patch


# -----------------------------
# Feature colors (edit freely)
# -----------------------------
FEATURE_COLORS = {
    "LTR":  "#333333",
    "GAG":  "#fc850e",
    "PROT": "#f80bfb",
    "RT":   "#0808f7",
    "RH":   "#fc1413",
    "INT":  "#05fc09",
    "CH":   "#20b5f4",
    "CHD":  "#20b5f4",
    "CHDCR":"#7a001a",
    "ARH":  "#70403c",
    "ENDO": "#a67c52",   # add ENDO explicitly so it has a stable color if present
}
DEFAULT_COLOR = "#AAAAAA"

FEATURE_ORDER = ["LTR", "GAG", "PROT", "RT", "RH", "INT", "CH", "CHD", "CHDCR", "ARH", "ENDO"]


# -----------------------------
# Data structures
# -----------------------------
@dataclass
class Feature:
    name: str
    start: int
    end: int

    def clamp(self, lo: int, hi: int) -> "Feature":
        s = max(lo, min(self.start, hi))
        e = max(lo, min(self.end, hi))
        if e < s:
            s, e = e, s
        return Feature(self.name, s, e)


@dataclass
class Element:
    element_id: str  # chrom:start-end
    chrom: str
    start: int
    end: int
    length: int      # end - start (matches your example)
    ltr_len: int
    cls: str
    superfamily: str
    family: str
    shift: int = 0

    strand: str = "?"
    strand_counts: Dict[str, int] = field(default_factory=lambda: {"+": 0, "-": 0})

    proteins: List[Feature] = field(default_factory=list)

    k2p: Optional[float] = None   # <-- NEW (TSV column 11, 0-based index 10)


# -----------------------------
# Utilities
# -----------------------------
def safe_name(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", s).strip("_") or "NA"


def median_int(vals: List[int]) -> int:
    if not vals:
        return 0
    vals2 = sorted(vals)
    return int(vals2[len(vals2)//2])


def percentile(vals: List[float], p: float) -> float:
    """p in [0,100]. Simple linear interpolation percentile."""
    if not vals:
        return float("nan")
    xs = sorted(vals)
    if len(xs) == 1:
        return float(xs[0])
    k = (len(xs) - 1) * (p / 100.0)
    f = int(k)
    c = min(f + 1, len(xs) - 1)
    if c == f:
        return float(xs[f])
    d = k - f
    return float(xs[f] * (1 - d) + xs[c] * d)


def bootstrap_ci(vals: List[int], boot: int, alpha: float = 0.05, stat="median") -> Tuple[int, int, int]:
    """
    Returns (center, lo, hi) as ints.
    center uses median (default) or mean.
    CI is bootstrap percentile CI of the chosen stat.
    """
    if not vals:
        return (0, 0, 0)

    if len(vals) == 1:
        v = int(vals[0])
        return (v, v, v)

    if stat == "mean":
        center = sum(vals) / len(vals)
        def f(sample): return sum(sample) / len(sample)
    else:
        center = median_int(vals)
        def f(sample): return median_int(sample)

    stats = []
    n = len(vals)
    for _ in range(boot):
        samp = [vals[random.randrange(n)] for _ in range(n)]
        stats.append(f(samp))

    lo = percentile(stats, 100 * (alpha / 2))
    hi = percentile(stats, 100 * (1 - alpha / 2))

    return (int(round(center)), int(round(lo)), int(round(hi)))


# -----------------------------
# Parsing
# -----------------------------
def parse_tsv(tsv_path: str) -> Dict[str, Element]:
    elements: Dict[str, Element] = {}
    print(f"[INFO] Reading TSV: {tsv_path}")

    with open(tsv_path, "r", encoding="utf-8") as f:
        for ln, line in enumerate(f, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                print(f"[WARN] TSV line {ln}: expected >=3 columns, got {len(parts)}. Skipping.")
                continue

            raw_id = parts[0]
            try:
                ltr1 = int(float(parts[1]))
                ltr2 = int(float(parts[2]))
            except ValueError:
                print(f"[WARN] TSV line {ln}: cannot parse LTR lengths from '{parts[1]}' '{parts[2]}'. Skipping.")
                continue

            # Column 11 (1-based) = index 10 (0-based): K2P
            k2p = None
            if len(parts) >= 11:
                try:
                    k2p = float(parts[10])
                except ValueError:
                    k2p = None


            if "#" in raw_id:
                elem_part, class_part = raw_id.split("#", 1)
            else:
                elem_part, class_part = raw_id, "NA/NA/NA"

            m = re.match(r"^(.+):(\d+)-(\d+)$", elem_part)
            if not m:
                print(f"[WARN] TSV line {ln}: cannot parse element coords from '{elem_part}'. Skipping.")
                continue
            chrom, s_str, e_str = m.group(1), m.group(2), m.group(3)
            s0, e0 = int(s_str), int(e_str)
            length = e0 - s0
            if length <= 0:
                print(f"[WARN] TSV line {ln}: non-positive length for '{elem_part}'. Skipping.")
                continue

            class_bits = class_part.split("/")
            cls = class_bits[0] if len(class_bits) > 0 else "NA"
            sup = class_bits[1] if len(class_bits) > 1 else "NA"
            fam = class_bits[2] if len(class_bits) > 2 else "NA"

            element_id = elem_part
            if element_id in elements:
                print(f"[WARN] TSV line {ln}: duplicate element_id '{element_id}'. Keeping first, ignoring this line.")
                continue

            elements[element_id] = Element(
                element_id=element_id,
                chrom=chrom,
                start=s0,
                end=e0,
                length=length,
                ltr_len=max(ltr1, ltr2),
                cls=cls,
                superfamily=sup,
                family=fam,
                k2p=k2p,   # <-- NEW
            )

    print(f"[INFO] TSV loaded: {len(elements)} elements")
    for el in list(elements.values())[:3]:
        print(f"[DEBUG] TSV sample: {el.element_id} len={el.length} ltr={el.ltr_len} fam={el.superfamily}/{el.family}")
    return elements


def parse_scn(scn_path: str) -> Dict[str, int]:
    shifts: Dict[str, int] = {}
    print(f"[INFO] Reading SCN: {scn_path}")

    with open(scn_path, "r", encoding="utf-8") as f:
        for ln, line in enumerate(f, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = re.split(r"\s+", line)
            if len(parts) < 12:
                print(f"[WARN] SCN line {ln}: expected >=12 columns, got {len(parts)}. Skipping.")
                continue
            try:
                s0 = int(parts[0])
                e0 = int(parts[1])
                shift = int(float(parts[5]))
                chrom = parts[11]
            except ValueError:
                print(f"[WARN] SCN line {ln}: cannot parse numeric fields. Skipping.")
                continue

            element_id = f"{chrom}:{s0}-{e0}"
            if element_id in shifts:
                if shifts[element_id] != shift:
                    print(f"[WARN] SCN line {ln}: duplicate '{element_id}' with different shift "
                          f"{shifts[element_id]} vs {shift}. Keeping first.")
                continue
            shifts[element_id] = shift

    print(f"[INFO] SCN loaded: {len(shifts)} element shifts")
    for k, v in list(shifts.items())[:5]:
        print(f"[DEBUG] SCN sample: {k} shift={v}")
    return shifts


def extract_attr(col9: str, key: str) -> Optional[str]:
    m = re.search(rf"(?:^|;){re.escape(key)}=([^;]+)", col9)
    return m.group(1) if m else None


def normalize_protein_name(raw: str) -> str:
    p = raw.strip().upper()
    p = re.sub(r"[^A-Z0-9]+", "", p)

    if p == "INTEGRASE":
        return "INT"
    if p == "PROTEASE":
        return "PROT"
    if p == "RNASEH":
        return "RH"
    if p in ("ARH", "ARNASEH"):
        return "ARH"
    if p.startswith("CH"):
        return p
    return p


def parse_gff_domains(gff_path: str) -> List[Tuple[str, str, int, int, str]]:
    out: List[Tuple[str, str, int, int, str]] = []
    print(f"[INFO] Reading GFF: {gff_path}")

    with open(gff_path, "r", encoding="utf-8") as f:
        for ln, line in enumerate(f, start=1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                print(f"[WARN] GFF line {ln}: expected 9 columns, got {len(parts)}. Skipping.")
                continue

            strand = parts[6].strip()
            if strand not in ("+", "-"):
                strand = "?"

            try:
                gff_start = int(parts[3])
                gff_end = int(parts[4])
            except ValueError:
                print(f"[WARN] GFF line {ln}: cannot parse start/end '{parts[3]}' '{parts[4]}'. Skipping.")
                continue

            attr = parts[8]
            id_val = extract_attr(attr, "ID")
            if not id_val:
                print(f"[WARN] GFF line {ln}: no ID= in col9. Skipping.")
                continue

            element_part = id_val.split("|", 1)[0].split(";", 1)[0]
            if not re.match(r"^.+:\d+-\d+$", element_part):
                print(f"[WARN] GFF line {ln}: cannot parse element_part '{element_part}'. Skipping.")
                continue

            gene = extract_attr(attr, "gene")
            if gene:
                prot = normalize_protein_name(gene)
            else:
                tail = id_val.split(":")[-1]
                tail = re.sub(r"^TY\d+-", "", tail, flags=re.IGNORECASE)
                tail = re.sub(r"^TY3-", "", tail, flags=re.IGNORECASE)
                tail = re.sub(r"^TY1-", "", tail, flags=re.IGNORECASE)
                tail = re.sub(r"^.*-", "", tail)
                prot = normalize_protein_name(tail)

            out.append((element_part, prot, gff_start, gff_end, strand))

    print(f"[INFO] GFF loaded: {len(out)} protein features")
    for tup in out[:6]:
        print(f"[DEBUG] GFF sample: element={tup[0]} prot={tup[1]} {tup[2]}-{tup[3]} strand={tup[4]}")
    return out


# -----------------------------
# Join TSV + SCN + GFF
# -----------------------------
def attach_shifts(elements: Dict[str, Element], shifts: Dict[str, int]) -> None:
    hit = 0
    for eid, el in elements.items():
        if eid in shifts:
            el.shift = shifts[eid]
            hit += 1
    print(f"[INFO] Shift attachment: matched {hit}/{len(elements)} TSV elements to SCN shifts")
    missing = [k for k in shifts.keys() if k not in elements]
    if missing:
        print(f"[WARN] SCN has {len(missing)} elements not present in TSV. First 5:")
        for k in missing[:5]:
            print(f"       - {k}")


def attach_proteins_and_strand(elements: Dict[str, Element], gff_feats: List[Tuple[str, str, int, int, str]]) -> None:
    matched = 0
    unmatched = 0

    for (eid, prot, gs, ge, strand) in gff_feats:
        el = elements.get(eid)
        if el is None:
            unmatched += 1
            continue

        if strand in ("+", "-"):
            el.strand_counts[strand] = el.strand_counts.get(strand, 0) + 1

        adj_start = gs + el.shift
        adj_end = ge + el.shift

        rel_s = adj_start - el.start + 1
        rel_e = adj_end - el.start + 1

        feat = Feature(prot, int(rel_s), int(rel_e)).clamp(1, el.length)
        el.proteins.append(feat)
        matched += 1

    print(f"[INFO] Protein attachment: matched {matched} GFF features to TSV elements; unmatched={unmatched}")

    # decide strand
    plus_major = minus_major = unknown = 0
    for el in elements.values():
        p = el.strand_counts.get("+", 0)
        m = el.strand_counts.get("-", 0)
        if p == 0 and m == 0:
            el.strand = "?"
            unknown += 1
        elif p >= m:
            el.strand = "+"
            plus_major += 1
        else:
            el.strand = "-"
            minus_major += 1

    print(f"[INFO] Strand inference (majority across domains): +={plus_major}, -={minus_major}, ?={unknown}")

    shown = 0
    for eid, el in elements.items():
        if el.proteins:
            prots = ", ".join([f"{p.name}:{p.start}-{p.end}" for p in sorted(el.proteins, key=lambda x: x.start)[:10]])
            print(f"[DEBUG] Proteins for {eid} (strand={el.strand}, votes={el.strand_counts}, shift={el.shift}, len={el.length}, ltr={el.ltr_len}): {prots}")
            shown += 1
        if shown >= 3:
            break


def flip_features_for_minus_strand(el: Element) -> None:
    if el.strand != "-":
        return
    L = el.length
    flipped = []
    for f in el.proteins:
        new_s = L - f.end + 1
        new_e = L - f.start + 1
        flipped.append(Feature(f.name, new_s, new_e).clamp(1, L))
    el.proteins = flipped
    el.strand = "+"
    print(f"[DEBUG] Flipped element to '+' orientation for plotting: {el.element_id}")


def normalize_all_elements_strand(elements: Dict[str, Element]) -> None:
    n_flip = n_plus = n_unknown = 0
    for el in elements.values():
        if el.strand == "-":
            flip_features_for_minus_strand(el)
            n_flip += 1
        elif el.strand == "+":
            n_plus += 1
        else:
            n_unknown += 1
    print(f"[INFO] Strand normalization: flipped={n_flip}, already_plus={n_plus}, unknown_left_as_is={n_unknown}")


# -----------------------------
# Plotting helpers
# -----------------------------
def build_legend_handles(features_present: List[str]) -> List[Patch]:
    ordered = [f for f in FEATURE_ORDER if f in features_present]
    extras = sorted([f for f in features_present if f not in ordered])
    final = ordered + extras

    handles = []
    for f in final:
        col = FEATURE_COLORS.get(f, DEFAULT_COLOR)
        handles.append(Patch(facecolor=col, edgecolor="black", label=f))
    return handles


def draw_element(ax, y: float, el_len: int, ltr_len: int, proteins: List[Feature],
                 label: Optional[str] = None, height: float = 0.6):
    ax.add_patch(plt.Rectangle((0, y - height/2), el_len, height, fill=False, linewidth=1.0))

    # LTRs
    ltr_col = FEATURE_COLORS["LTR"]
    ax.add_patch(plt.Rectangle((0, y - height/2), ltr_len, height,
                               facecolor=ltr_col, edgecolor="black", linewidth=0.6))
    ax.add_patch(plt.Rectangle((max(0, el_len - ltr_len), y - height/2), ltr_len, height,
                               facecolor=ltr_col, edgecolor="black", linewidth=0.6))

    for p in proteins:
        w = max(1, p.end - p.start + 1)
        x = max(0, p.start - 1)
        col = FEATURE_COLORS.get(p.name, DEFAULT_COLOR)
        ax.add_patch(plt.Rectangle((x, y - height/2), w, height,
                                   facecolor=col, edgecolor="black", linewidth=0.6))

    if label:
        ax.text(-0.01 * el_len, y, label, ha="right", va="center", fontsize=8)


def draw_whisker(ax, y: float, lo: int, hi: int, color: str, lw: float = 1.2):
    """
    A simple CI whisker: horizontal line from lo..hi with small end caps.
    """
    ax.plot([lo, hi], [y, y], linewidth=lw, color=color)
    cap = 0.06
    ax.plot([lo, lo], [y - cap, y + cap], linewidth=lw, color=color)
    ax.plot([hi, hi], [y - cap, y + cap], linewidth=lw, color=color)


# -----------------------------
# Average computation w/ rarity filter + CI
# -----------------------------
@dataclass
class AvgFeatureCI:
    name: str
    start_center: int
    start_lo: int
    start_hi: int
    end_center: int
    end_lo: int
    end_hi: int


@dataclass
class FamilyAverages:
    n_elements: int
    len_center: int
    len_lo: int
    len_hi: int
    ltr_center: int
    ltr_lo: int
    ltr_hi: int
    features: List[AvgFeatureCI]  # proteins only (LTR handled separately)
    features_present: List[str]   # for legend


def compute_family_averages_with_ci(
    fam_elements: List[Element],
    min_presence: float,
    boot: int,
    alpha: float = 0.05,
) -> FamilyAverages:
    """
    Strand-normalization must already be done before calling this.
    Uses bootstrap percentile CI on median by default.
    """

    n = len(fam_elements)
    lengths = [e.length for e in fam_elements]
    ltrs = [e.ltr_len for e in fam_elements]

    len_center, len_lo, len_hi = bootstrap_ci(lengths, boot=boot, alpha=alpha, stat="median")
    ltr_center, ltr_lo, ltr_hi = bootstrap_ci(ltrs, boot=boot, alpha=alpha, stat="median")

    # Count presence (protein appears in an element if at least one hit of that name exists)
    prot_presence: Dict[str, int] = {}
    for e in fam_elements:
        present = set(p.name for p in e.proteins)
        for name in present:
            prot_presence[name] = prot_presence.get(name, 0) + 1

    # Filter by min_presence
    keep = []
    for name, cnt in prot_presence.items():
        frac = cnt / n if n else 0.0
        if frac >= min_presence:
            keep.append(name)
    keep = sorted(keep, key=lambda x: (FEATURE_ORDER.index(x) if x in FEATURE_ORDER else 999, x))

    # For each kept protein, gather per-element “representative” start/end
    # If multiple hits exist, use earliest start and latest end (covers fragmented hits).
    features_ci: List[AvgFeatureCI] = []
    for prot in keep:
        starts = []
        ends = []
        for e in fam_elements:
            hits = [p for p in e.proteins if p.name == prot]
            if not hits:
                continue
            s = min(h.start for h in hits)
            en = max(h.end for h in hits)
            starts.append(int(s))
            ends.append(int(en))

        # CI on median start/end
        s_center, s_lo, s_hi = bootstrap_ci(starts, boot=boot, alpha=alpha, stat="median")
        e_center, e_lo, e_hi = bootstrap_ci(ends, boot=boot, alpha=alpha, stat="median")

        features_ci.append(AvgFeatureCI(
            name=prot,
            start_center=s_center, start_lo=s_lo, start_hi=s_hi,
            end_center=e_center, end_lo=e_lo, end_hi=e_hi
        ))

    # Legend includes LTR + kept proteins
    features_present = ["LTR"] + [f.name for f in features_ci]
    return FamilyAverages(
        n_elements=n,
        len_center=len_center, len_lo=len_lo, len_hi=len_hi,
        ltr_center=ltr_center, ltr_lo=ltr_lo, ltr_hi=ltr_hi,
        features=features_ci,
        features_present=features_present
    )


def print_family_summary(family_key: str, avg: FamilyAverages):
    print(f"\n[SUMMARY] {family_key}  (N={avg.n_elements})")
    print(f"  Average LTR-RT len: {avg.len_center} (95% CI: {avg.len_lo}-{avg.len_hi})")
    print(f"  Average LTR len:    {avg.ltr_center} (95% CI: {avg.ltr_lo}-{avg.ltr_hi})")
    for feat in avg.features:
        print(f"  Average {feat.name} start: {feat.start_center} (95% CI: {feat.start_lo}-{feat.start_hi})")
        print(f"  Average {feat.name} end:   {feat.end_center} (95% CI: {feat.end_lo}-{feat.end_hi})")


# -----------------------------
# Plotting: average with whiskers
# -----------------------------
def plot_family_average_with_whiskers(family_key: str, avg: FamilyAverages, out_path: str):
    print(f"[INFO] Plot average (CI whiskers): {family_key} -> {out_path}")

    legend_handles = build_legend_handles(avg.features_present)

    fig, ax = plt.subplots(figsize=(12.5, 3.2))

    # Backbone outline for "average length" (center)
    L = max(1, avg.len_center)
    ltr = max(1, min(avg.ltr_center, L))

    # Draw the average element cartoon (rectangles at centers)
    # Build features at center positions
    center_feats = []
    for f in avg.features:
        s = max(1, min(f.start_center, L))
        e = max(1, min(f.end_center, L))
        if e < s:
            s, e = e, s
        center_feats.append(Feature(f.name, s, e))

    draw_element(ax, y=1.0, el_len=L, ltr_len=ltr, proteins=sorted(center_feats, key=lambda x: x.start),
                 label=None, height=0.55)

    # Draw whiskers for LTR-RT length CI (as a line under the backbone)
    # We show the CI span on x-axis at a lower y
    draw_whisker(ax, y=0.55, lo=max(1, avg.len_lo), hi=max(1, avg.len_hi), color="black", lw=1.4)
    ax.text(L * 0.01, 0.40, "Element length 95% CI", fontsize=8, ha="left", va="center")

    # Draw whiskers for LTR length CI (left and right)
    # Left LTR CI shows length span from 1..ltr_len
    # We'll render whisker for "LTR length" as span on x-axis from lo..hi near left,
    # and also mirrored near right.
    ltr_lo = max(1, min(avg.ltr_lo, L))
    ltr_hi = max(1, min(avg.ltr_hi, L))
    # Left
    draw_whisker(ax, y=1.33, lo=0, hi=ltr_hi, color=FEATURE_COLORS["LTR"], lw=1.2)
    draw_whisker(ax, y=1.33, lo=0, hi=ltr_lo, color=FEATURE_COLORS["LTR"], lw=1.2)
    # Right (display as span [L-ltr_hi .. L] and [L-ltr_lo .. L])
    draw_whisker(ax, y=1.33, lo=max(0, L - ltr_hi), hi=L, color=FEATURE_COLORS["LTR"], lw=1.2)
    draw_whisker(ax, y=1.33, lo=max(0, L - ltr_lo), hi=L, color=FEATURE_COLORS["LTR"], lw=1.2)
    ax.text(L * 0.5, 1.43, "LTR length 95% CI (both ends)", fontsize=8, ha="center", va="center")

    # For each protein: whiskers for start and end CIs
    # Put start whisker slightly above, end whisker slightly below the feature bar
    for f in avg.features:
        col = FEATURE_COLORS.get(f.name, DEFAULT_COLOR)
        # start CI whisker
        draw_whisker(ax, y=1.18, lo=max(1, min(f.start_lo, L)), hi=max(1, min(f.start_hi, L)), color=col, lw=1.2)
        # end CI whisker
        draw_whisker(ax, y=0.82, lo=max(1, min(f.end_lo, L)), hi=max(1, min(f.end_hi, L)), color=col, lw=1.2)

    ax.set_title(f"{family_key} : Structure of average element (rare proteins filtered; whiskers=95% CI)")
    ax.set_xlim(0, max(avg.len_hi, L) * 1.02)
    ax.set_ylim(0.25, 1.65)
    ax.set_yticks([])
    ax.set_xlabel("Position [bp] (relative)")

    ax.legend(handles=legend_handles, loc="upper right", frameon=False, ncol=min(6, len(legend_handles)))
    plt.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


# -----------------------------
# Individual plots + all-elements
# -----------------------------
def plot_family_individual(family_key: str, fam_elements: List[Element], out_path: str, max_elements: int = 200):
    print(f"[INFO] Plot individuals: {family_key} -> {out_path}")
    fam_elements2 = fam_elements
    if len(fam_elements2) > max_elements:
        print(f"[WARN] Too many elements in {family_key} ({len(fam_elements2)}). Plotting first {max_elements}.")
        fam_elements2 = fam_elements2[:max_elements]

    max_len = max(e.length for e in fam_elements2)

    prots_seen = sorted(set(p.name for e in fam_elements2 for p in e.proteins))
    features_present = ["LTR"] + prots_seen
    legend_handles = build_legend_handles(features_present)

    fig_h = max(3.0, 0.18 * len(fam_elements2) + 1.8)
    fig, ax = plt.subplots(figsize=(14, fig_h))

    for i, e in enumerate(fam_elements2):
        y = len(fam_elements2) - i
        prots = sorted(e.proteins, key=lambda x: x.start)
        k2p_txt = "NA" if e.k2p is None else f"{e.k2p * 100:.1f}%"
        draw_element(ax, y=y, el_len=e.length, ltr_len=e.ltr_len, proteins=prots,
                     label=f"{e.element_id}  {k2p_txt}", height=0.65)

    ax.set_title(f"{family_key} : Structure of individual elements (n={len(fam_elements2)} shown)")
    ax.set_xlim(0, max_len * 1.02)
    ax.set_ylim(0, len(fam_elements2) + 1)
    ax.set_yticks([])
    ax.set_xlabel("Position [bp] (relative)")

    # bottom-right as requested
    ax.legend(handles=legend_handles, loc="lower right", frameon=False, ncol=min(6, len(legend_handles)))
    plt.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def plot_all_elements_average(family_avg_map: Dict[str, FamilyAverages], out_path: str):
    fam_keys = sorted(family_avg_map.keys())
    print(f"[INFO] Plot all_elements averages -> {out_path}")
    print(f"[DEBUG] Families included: {len(fam_keys)}")

    # union of proteins included in averages (after rarity filter)
    all_prots_seen = set()
    max_len_hi = 1
    for fk in fam_keys:
        avg = family_avg_map[fk]
        max_len_hi = max(max_len_hi, avg.len_hi)
        for feat in avg.features:
            all_prots_seen.add(feat.name)

    features_present = ["LTR"] + sorted(all_prots_seen)
    legend_handles = build_legend_handles(features_present)

    fig_h = max(3.0, 0.35 * len(fam_keys) + 2.0)
    fig, ax = plt.subplots(figsize=(14, fig_h))

    for i, fk in enumerate(fam_keys):
        avg = family_avg_map[fk]
        y = len(fam_keys) - i

        L = max(1, avg.len_center)
        ltr = max(1, min(avg.ltr_center, L))
        feats = []
        for f in avg.features:
            s = max(1, min(f.start_center, L))
            e = max(1, min(f.end_center, L))
            if e < s:
                s, e = e, s
            feats.append(Feature(f.name, s, e))


        draw_element(ax, y=y, el_len=L, ltr_len=ltr, proteins=sorted(feats, key=lambda x: x.start),
                     label=fk, height=0.65)

    ax.set_title("all_elements : Average structure per family (rare proteins filtered)")
    ax.set_xlim(0, max_len_hi * 1.02)
    ax.set_ylim(0, len(fam_keys) + 1)
    ax.set_yticks([])
    ax.set_xlabel("Position [bp] (relative)")

    ax.legend(handles=legend_handles, loc="upper right", frameon=False, ncol=min(6, len(legend_handles)))
    plt.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


# -----------------------------
# Main
# -----------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", required=True, help="LTR-RT TSV file")
    ap.add_argument("--scn", required=True, help="SCN file")
    ap.add_argument("--gff", required=True, help="Domain HMM gff file")
    ap.add_argument("--out_dir", required=True, help="Output directory")

    ap.add_argument("--max_individual", type=int, default=200, help="Max elements per-family in individual plot")
    ap.add_argument("--min_presence", type=float, default=0.10,
                    help="Min fraction of elements a protein must appear in to be included in the AVERAGE plot (default 0.10)")
    ap.add_argument("--boot", type=int, default=2000,
                    help="Bootstrap iterations for 95%% CI (default 2000). Increase for smoother CI.")
    ap.add_argument("--seed", type=int, default=1, help="Random seed for bootstrapping (default 1)")

    args = ap.parse_args()
    random.seed(args.seed)

    os.makedirs(args.out_dir, exist_ok=True)

    elements = parse_tsv(args.tsv)
    shifts = parse_scn(args.scn)
    gff_feats = parse_gff_domains(args.gff)

    attach_shifts(elements, shifts)
    attach_proteins_and_strand(elements, gff_feats)

    # Strand normalize (critical for averages + ordering)
    normalize_all_elements_strand(elements)

    # Group by superfamily/family
    family_to_elements: Dict[str, List[Element]] = {}
    for el in elements.values():
        fam_key = f"{el.superfamily}/{el.family}"
        family_to_elements.setdefault(fam_key, []).append(el)

    for fk in family_to_elements:
        family_to_elements[fk].sort(key=lambda e: e.length, reverse=True)

    print(f"[INFO] Grouped into {len(family_to_elements)} families")

    # Compute averages w/ rarity filter + CI, print summaries, and plot
    family_avg_map: Dict[str, FamilyAverages] = {}

    for fk, fam_elements in sorted(family_to_elements.items(), key=lambda x: x[0]):
        # Averages
        avg = compute_family_averages_with_ci(
            fam_elements=fam_elements,
            min_presence=args.min_presence,
            boot=args.boot
        )
        family_avg_map[fk] = avg

        # Print summary to terminal
        # Also print which proteins got included/excluded
        present_counts = {}
        for e in fam_elements:
            for p in set(pp.name for pp in e.proteins):
                present_counts[p] = present_counts.get(p, 0) + 1
        excluded = []
        for p, c in sorted(present_counts.items(), key=lambda x: (-x[1], x[0])):
            if p not in [ff.name for ff in avg.features]:
                frac = c / len(fam_elements)
                excluded.append(f"{p} ({c}/{len(fam_elements)}={frac:.3f})")
        print_family_summary(fk, avg)
        if excluded:
            print(f"  Excluded rare proteins (<{args.min_presence:.2f} presence): " + ", ".join(excluded[:12]) +
                  (" ..." if len(excluded) > 12 else ""))

        # Plot average (with whiskers)
        safe_fk = safe_name(fk)
        out_avg = os.path.join(args.out_dir, f"{safe_fk}_average.pdf")
        plot_family_average_with_whiskers(fk, avg, out_avg)

        # Plot individuals
        out_ind = os.path.join(args.out_dir, f"{safe_fk}_individual.pdf")
        plot_family_individual(fk, fam_elements, out_ind, max_elements=args.max_individual)

    # Global plot of averages
    out_all = os.path.join(args.out_dir, "all_elements.pdf")
    plot_all_elements_average(family_avg_map, out_all)

    print("\n[INFO] Done.")
    print(f"[INFO] Output directory: {args.out_dir}")
    print(f"[INFO] Average plots exclude proteins with presence < {args.min_presence:.2f} and show 95% bootstrap CI whiskers.")


if __name__ == "__main__":
    main()
