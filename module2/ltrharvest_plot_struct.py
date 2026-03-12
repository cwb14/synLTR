#!/usr/bin/env python3
"""
plot.py

Usage:
  python plot.py --tsv ltrrt.tsv --scn file.scn --gff hmm.gff --out_dir dir_name [--fai genome.fai] [--html]
  python plot.py --scn <(awk '{if (NF>12) {for(i=1;i<=NF;i++) if(i!=NF-1) printf "%s%s",$i,(i==NF?"":" "); print ""} else print}' Col0_ltr_r1.work/Col0_ltr_r1.ltrfinder.stitched.scn) --gff   Col0_ltr_r1.work/Col0_ltr_r1.ltrtools.intact_for_tesorter.fa.rexdb-plant.dom.gff3 --tsv Col0_ltr_r1_kmer2ltr_dedup --out_dir plot_struct --fai ../Col0.fa.fai --html
Outputs:
  <out_dir>/<safe_family>_average.pdf   (or embedded in HTML if --html)
  <out_dir>/<safe_family>_individual.pdf (or embedded in HTML if --html)
  <out_dir>/all_elements.pdf            (or embedded in HTML if --html)
  <out_dir>/chrom_plots/<chrom>.html    (always HTML, if --fai provided)
  <out_dir>/index.html                  (if --html, master page with everything)

Key features:
- Strand-normalizes all elements for plotting/averaging (flips '-' into '+').
- Auto-detects whether GFF coordinates are genomic or internal-to-element (shift vs no-shift).
- Average plots EXCLUDE rare proteins unless present in >= --min_presence fraction of elements.
- Individual plots always show ALL proteins for each element.
- Average plots include whiskers showing 95% CI for start/end of each feature (bootstrap CI).
- Per-chromosome positional plots are interactive IGV-like HTML viewers.
- Prints per-family summary stats + 95% CI to terminal.
- Flags potential false positives: elements with outlier-long lengths that lack expected proteins.
"""

import argparse
import base64
import io
import json
import math
import os
import re
import random
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional, Set

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.backends.backend_pdf import PdfPages


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
    "ENDO": "#a67c52",
}
DEFAULT_COLOR = "#AAAAAA"
FLAGGED_COLOR = "#ff4444"

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
    element_id: str
    chrom: str
    start: int
    end: int
    length: int
    ltr_len: int
    cls: str
    superfamily: str
    family: str
    shift: int = 0
    ltr_len_raw: int = 0   # column 2 of TSV (LTR_LEN)
    aln_len_raw: int = 0   # column 3 of TSV (ALN_LEN)

    strand: str = "?"
    strand_counts: Dict[str, int] = field(default_factory=lambda: {"+": 0, "-": 0})

    proteins: List[Feature] = field(default_factory=list)

    k2p: Optional[float] = None


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


def fig_to_base64_png(fig) -> str:
    """Convert a matplotlib figure to a base64-encoded PNG string."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode("ascii")
    buf.close()
    return b64


# -----------------------------
# Outlier detection & false-positive flagging
# -----------------------------
def _distributional_fence(lengths: List[float], n: int, iqr_mul: float = 1.5) -> float:
    """
    Compute a distributional upper fence for outlier detection.

    For N >= 15: IQR method (Q3 + iqr_mul*IQR).  Default iqr_mul=1.5 (Tukey).
        Lower values (e.g. 1.2–1.4) tighten the fence; higher values loosen it.
    For 5 <= N < 15: MAD-based method, which is more robust for small samples.
        MAD = median(|xi - median(x)|)
        Scaled MAD (consistent estimator of sigma) = MAD / 0.6745
        Upper fence = median + (iqr_mul * 2) * scaled_MAD
        (The MAD multiplier is scaled proportionally to iqr_mul so tuning
        has a consistent effect across both methods.)

    The MAD approach is preferred for small samples because IQR estimates
    become unreliable when quartiles are computed from few data points.
    The factor 0.6745 makes MAD a consistent estimator of the standard
    deviation for normally distributed data.
    """
    if n >= 15:
        q1 = percentile(lengths, 25.0)
        q3 = percentile(lengths, 75.0)
        iqr = q3 - q1
        return q3 + iqr_mul * iqr
    else:
        med = percentile(lengths, 50.0)
        abs_devs = sorted(abs(x - med) for x in lengths)
        mad = percentile(abs_devs, 50.0)
        if mad < 1e-9:
            return med * 1.15
        scaled_mad = mad / 0.6745
        # Scale MAD multiplier proportionally: default iqr_mul=1.5 → 3.0σ
        mad_mul = iqr_mul * 2.0
        return med + mad_mul * scaled_mad


def _gap_fence(lengths: List[float], n: int) -> float:
    """
    Detect a natural breakpoint in the upper portion of the length distribution.

    When outliers constitute a large fraction of a small sample (e.g. 2 of 9),
    they inflate IQR/MAD and hide themselves.  A biologist looking at the plot
    spots these as a visible *gap* between the main cluster and the long tail.

    Algorithm:
      1. Sort lengths, compute ALL consecutive gaps to establish baseline
         spacing (median of all gaps).
      2. Only consider gaps whose *lower* value is at or above the median
         length as candidate breakpoints (we don't want to split the main
         body of the distribution).
      3. A gap qualifies as a breakpoint if it is:
           a) > 2× the median gap (across ALL gaps), AND
           b) > 1000 bp absolute (avoids triggering on noise).
      4. Among qualifying gaps, pick the largest one.
      5. The gap fence is placed at the value just below that gap, plus a
         small tolerance (half the median gap) so elements right at the
         boundary aren't caught by rounding.

    Returns float('inf') if no qualifying gap is found.
    """
    if n < 5:
        return float("inf")

    sorted_L = sorted(lengths)
    med = percentile(sorted_L, 50.0)

    # Compute ALL consecutive gaps for baseline spacing
    all_gaps = []
    for i in range(len(sorted_L) - 1):
        all_gaps.append(sorted_L[i + 1] - sorted_L[i])

    if not all_gaps:
        return float("inf")

    median_gap = percentile(sorted(all_gaps), 50.0)

    # Candidate gaps: only in the upper half (lower value >= median length)
    upper_gaps = []  # (gap_size, index_of_lower_value)
    for i in range(len(sorted_L) - 1):
        if sorted_L[i] >= med:
            gap = sorted_L[i + 1] - sorted_L[i]
            upper_gaps.append((gap, i))

    if not upper_gaps:
        return float("inf")

    # Qualifying gaps: large relative to typical spacing AND absolute minimum
    min_abs_gap = 1000.0
    min_relative_factor = 2.0
    threshold = max(min_abs_gap, median_gap * min_relative_factor)

    qualifying = [(g, idx) for g, idx in upper_gaps if g > threshold]
    if not qualifying:
        return float("inf")

    # Pick the largest qualifying gap
    qualifying.sort(key=lambda x: x[0], reverse=True)
    _best_gap, best_idx = qualifying[0]

    # Fence just above the value below the gap (small tolerance)
    tolerance = median_gap * 0.5
    return sorted_L[best_idx] + tolerance


def _compute_upper_fence(lengths: List[float], n: int, iqr_mul: float = 1.5) -> float:
    """
    Final upper fence = min(distributional fence, gap fence).

    The distributional fence (IQR or MAD) works well when outliers are a
    small fraction of the data.  The gap fence catches cases where outliers
    inflate the spread estimate and hide themselves — it looks for a natural
    discontinuity in the sorted lengths.  Taking the minimum of the two
    ensures both mechanisms can contribute.
    """
    d_fence = _distributional_fence(lengths, n, iqr_mul=iqr_mul)
    g_fence = _gap_fence(lengths, n)
    return min(d_fence, g_fence)


def _is_dubious_clade(family_key: str) -> bool:
    """Check if a family key indicates a dubious clade (mixture or unknown)."""
    parts = family_key.lower().replace("/", " ").split()
    return any(p in ("mixture", "unknown") for p in parts)


def _is_unknown_clade(family_key: str) -> bool:
    """Check if a family key indicates an 'unknown' clade specifically."""
    parts = family_key.lower().replace("/", " ").split()
    return "unknown" in parts


def detect_false_positives(
    fam_elements: List[Element],
    expected_proteins: List[str],
    family_key: str,
    min_family_size: int = 5,
    protein_recovery_frac: float = 0.50,
    iqr_mul: float = 0.5,
    min_ltr_aln: int = 120,
    dubious_k2p: float = 0.15,
    recovery_k2p: float = 0.10,
    skip_unknown_iqr: bool = False,
) -> Set[str]:
    """
    Identify potential false positives within a family.

    Two independent mechanisms:

    1. LENGTH OUTLIER DETECTION (all families with N >= min_family_size,
       unless skip_unknown_iqr is True for unknown clades):
       Uses IQR for N >= 15 or MAD for 5 <= N < 15, supplemented by gap
       analysis.  Among outliers, elements with sufficient expected proteins
       are "recovered" (considered real).  When recovery requires only 1
       protein, an additional K2P < recovery_k2p constraint is applied.
       Protein-rich elements (more distinct proteins than the minimum
       threshold) are recovered unconditionally.

    2. DUBIOUS-CLADE FLAGGING (mixture/unknown families only):
       Any element with K2P > dubious_k2p, or LTR_LEN < min_ltr_aln,
       or ALN_LEN < min_ltr_aln is flagged regardless of element length.

    Returns a set of element_ids flagged as potential false positives.
    """
    flagged: Set[str] = set()
    n = len(fam_elements)
    dubious = _is_dubious_clade(family_key)
    unknown = _is_unknown_clade(family_key)

    # --- Mechanism 2: dubious-clade quality flags ---
    if dubious:
        for el in fam_elements:
            reasons = []
            if el.k2p is not None and el.k2p > dubious_k2p:
                reasons.append(f"k2p>{dubious_k2p*100:.0f}%")
            if el.ltr_len_raw < min_ltr_aln:
                reasons.append(f"LTR_LEN<{min_ltr_aln}")
            if el.aln_len_raw < min_ltr_aln:
                reasons.append(f"ALN_LEN<{min_ltr_aln}")
            if reasons:
                flagged.add(el.element_id)

    # --- Mechanism 1: length outlier detection ---
    # Skip IQR for unknown clades if requested
    if unknown and skip_unknown_iqr:
        return flagged

    if n < min_family_size:
        return flagged  # too few for outlier analysis, return dubious flags only

    lengths = sorted(float(e.length) for e in fam_elements)
    upper_fence = _compute_upper_fence(lengths, n, iqr_mul=iqr_mul)

    # Identify outlier-long elements
    outlier_long = [e for e in fam_elements if e.length > upper_fence]
    if not outlier_long:
        return flagged

    n_expected = len(expected_proteins)

    # Compute the minimum number of proteins required for recovery
    # This is the smallest integer h such that h / n_expected > protein_recovery_frac
    if n_expected > 0:
        min_proteins_for_recovery = math.floor(n_expected * protein_recovery_frac) + 1
        # Edge case: if frac * n_expected is an exact integer, we still need +1
        # because the comparison is strictly greater than
        if n_expected * protein_recovery_frac == math.floor(n_expected * protein_recovery_frac):
            min_proteins_for_recovery = int(n_expected * protein_recovery_frac) + 1
    else:
        min_proteins_for_recovery = 0

    for el in outlier_long:
        if n_expected == 0:
            # Non-autonomous: no proteins to check, flag purely on length
            flagged.add(el.element_id)
        else:
            # Check how many expected proteins this element has
            el_prots = set(p.name for p in el.proteins)
            hits = sum(1 for ep in expected_proteins if ep in el_prots)
            total_distinct = len(el_prots)

            if hits >= min_proteins_for_recovery:
                # Has enough expected proteins.
                # Additional recovery layer: if the element has more total
                # distinct proteins than the minimum required, it's well-
                # supported and recovered unconditionally.  K2P stringency
                # only applies when the element barely meets the threshold.
                if total_distinct > min_proteins_for_recovery:
                    pass  # recovered — protein-rich element
                elif min_proteins_for_recovery == 1:
                    # Barely meets threshold with 1 protein — apply K2P gate
                    if el.k2p is None or el.k2p >= recovery_k2p:
                        flagged.add(el.element_id)
                    # else: recovered (has protein AND low K2P)
                # else: recovered (multiple expected proteins required and met)
            else:
                # Does NOT have enough expected proteins.
                # But check if it's protein-rich overall — an element with
                # many proteins is likely real even if they don't match the
                # family's typical set (e.g. a clade where avg only has 1
                # expected protein but the element has 5).
                if total_distinct > min_proteins_for_recovery:
                    pass  # recovered — protein-rich despite missing expected set
                else:
                    flagged.add(el.element_id)

    return flagged


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
                continue

            raw_id = parts[0]
            try:
                ltr1 = int(float(parts[1]))
                ltr2 = int(float(parts[2]))
            except ValueError:
                continue

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
                continue
            chrom, s_str, e_str = m.group(1), m.group(2), m.group(3)
            s0, e0 = int(s_str), int(e_str)
            length = e0 - s0
            if length <= 0:
                continue

            class_bits = class_part.split("/")
            cls = class_bits[0] if len(class_bits) > 0 else "NA"
            sup = class_bits[1] if len(class_bits) > 1 else "NA"
            fam = class_bits[2] if len(class_bits) > 2 else "NA"

            element_id = elem_part
            if element_id in elements:
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
                k2p=k2p,
                ltr_len_raw=ltr1,
                aln_len_raw=ltr2,
            )

    print(f"[INFO] TSV loaded: {len(elements)} elements")
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
                continue
            try:
                s0 = int(parts[0])
                e0 = int(parts[1])
                shift = int(float(parts[5]))
                chrom = parts[11]
            except ValueError:
                continue

            element_id = f"{chrom}:{s0}-{e0}"
            if element_id in shifts:
                continue
            shifts[element_id] = shift

    print(f"[INFO] SCN loaded: {len(shifts)} element shifts")
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
                continue

            strand = parts[6].strip()
            if strand not in ("+", "-"):
                strand = "?"

            try:
                gff_start = int(parts[3])
                gff_end = int(parts[4])
            except ValueError:
                continue

            attr = parts[8]
            id_val = extract_attr(attr, "ID")
            if not id_val:
                continue

            element_part = id_val.split("|", 1)[0].split(";", 1)[0]
            if not re.match(r"^.+:\d+-\d+$", element_part):
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
    return out


# -----------------------------
# Auto-detect GFF coordinate mode
# -----------------------------
def detect_gff_coordinate_mode(
    elements: Dict[str, Element],
    gff_feats: List[Tuple[str, str, int, int, str]],
) -> str:
    """
    Determine if GFF coordinates are:
      'genomic'  — true genomic coordinates (fall within element start..end)
      'internal' — internal-to-element (need shift to map into element coords)

    Strategy: for each GFF feature that maps to a known element, check if the
    GFF start..end falls within element.start..element.end.  If most do,
    coordinates are genomic.  Otherwise, they are internal.
    """
    genomic_hits = 0
    internal_hits = 0
    tested = 0

    for (eid, _prot, gs, ge, _strand) in gff_feats:
        el = elements.get(eid)
        if el is None:
            continue
        tested += 1
        if tested > 5000:
            break  # sample enough

        # Check if GFF coords fall within [el.start, el.end]
        if el.start <= gs <= el.end and el.start <= ge <= el.end:
            genomic_hits += 1
        else:
            internal_hits += 1

    if tested == 0:
        print("[WARN] No GFF features matched TSV elements for coordinate-mode detection; defaulting to 'internal'.")
        return "internal"

    frac_genomic = genomic_hits / tested
    if frac_genomic >= 0.8:
        mode = "genomic"
    else:
        mode = "internal"

    print(f"[INFO] GFF coordinate mode auto-detected: '{mode}' "
          f"(tested={tested}, genomic_hits={genomic_hits} ({frac_genomic:.1%}), internal_hits={internal_hits})")
    return mode


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


def attach_proteins_and_strand(
    elements: Dict[str, Element],
    gff_feats: List[Tuple[str, str, int, int, str]],
    coord_mode: str,
) -> None:
    """
    coord_mode: 'genomic' or 'internal'

    If 'genomic', GFF coords are true genomic — convert to element-relative by
    subtracting element.start (no shift needed).

    If 'internal', GFF coords are internal and need el.shift applied, then
    converted to element-relative.
    """
    matched = 0
    unmatched = 0

    for (eid, prot, gs, ge, strand) in gff_feats:
        el = elements.get(eid)
        if el is None:
            unmatched += 1
            continue

        if strand in ("+", "-"):
            el.strand_counts[strand] = el.strand_counts.get(strand, 0) + 1

        if coord_mode == "genomic":
            # GFF coords are genomic; convert to element-relative [1..length]
            rel_s = gs - el.start + 1
            rel_e = ge - el.start + 1
        else:
            # Internal coords: apply shift then convert
            adj_start = gs + el.shift
            adj_end = ge + el.shift
            rel_s = adj_start - el.start + 1
            rel_e = adj_end - el.start + 1

        feat = Feature(prot, int(rel_s), int(rel_e)).clamp(1, el.length)
        el.proteins.append(feat)
        matched += 1

    print(f"[INFO] Protein attachment ({coord_mode} mode): matched {matched} GFF features; unmatched={unmatched}")

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

    print(f"[INFO] Strand inference: +={plus_major}, -={minus_major}, ?={unknown}")


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
    print(f"[INFO] Strand normalization: flipped={n_flip}, already_plus={n_plus}, unknown={n_unknown}")


# -----------------------------
# Plotting helpers
# -----------------------------
def build_legend_handles(features_present: List[str], include_flagged: bool = False) -> List[Patch]:
    ordered = [f for f in FEATURE_ORDER if f in features_present]
    extras = sorted([f for f in features_present if f not in ordered])
    final = ordered + extras
    handles = []
    for f in final:
        col = FEATURE_COLORS.get(f, DEFAULT_COLOR)
        handles.append(Patch(facecolor=col, edgecolor="black", label=f))
    if include_flagged:
        handles.append(Patch(facecolor="none", edgecolor=FLAGGED_COLOR, linewidth=2.0,
                             label="Potential FP"))
    return handles


def draw_element(ax, y: float, el_len: int, ltr_len: int, proteins: List[Feature],
                 label: Optional[str] = None, height: float = 0.6,
                 flagged: bool = False):
    outline_color = FLAGGED_COLOR if flagged else "black"
    outline_lw = 2.0 if flagged else 1.0
    ax.add_patch(plt.Rectangle((0, y - height/2), el_len, height, fill=False,
                               edgecolor=outline_color, linewidth=outline_lw))
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
        label_color = FLAGGED_COLOR if flagged else "black"
        ax.text(-0.01 * el_len, y, label, ha="right", va="center", fontsize=8,
                color=label_color, fontweight="bold" if flagged else "normal")


def draw_whisker(ax, y: float, lo: int, hi: int, color: str, lw: float = 1.2):
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
    features: List[AvgFeatureCI]
    features_present: List[str]


def compute_family_averages_with_ci(
    fam_elements: List[Element],
    min_presence: float,
    boot: int,
    alpha: float = 0.05,
) -> FamilyAverages:
    n = len(fam_elements)
    lengths = [e.length for e in fam_elements]
    ltrs = [e.ltr_len for e in fam_elements]

    len_center, len_lo, len_hi = bootstrap_ci(lengths, boot=boot, alpha=alpha, stat="median")
    ltr_center, ltr_lo, ltr_hi = bootstrap_ci(ltrs, boot=boot, alpha=alpha, stat="median")

    prot_presence: Dict[str, int] = {}
    for e in fam_elements:
        present = set(p.name for p in e.proteins)
        for name in present:
            prot_presence[name] = prot_presence.get(name, 0) + 1

    keep = []
    for name, cnt in prot_presence.items():
        frac = cnt / n if n else 0.0
        if frac >= min_presence:
            keep.append(name)
    keep = sorted(keep, key=lambda x: (FEATURE_ORDER.index(x) if x in FEATURE_ORDER else 999, x))

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

        s_center, s_lo, s_hi = bootstrap_ci(starts, boot=boot, alpha=alpha, stat="median")
        e_center, e_lo, e_hi = bootstrap_ci(ends, boot=boot, alpha=alpha, stat="median")

        features_ci.append(AvgFeatureCI(
            name=prot,
            start_center=s_center, start_lo=s_lo, start_hi=s_hi,
            end_center=e_center, end_lo=e_lo, end_hi=e_hi
        ))

    features_present = ["LTR"] + [f.name for f in features_ci]
    return FamilyAverages(
        n_elements=n,
        len_center=len_center, len_lo=len_lo, len_hi=len_hi,
        ltr_center=ltr_center, ltr_lo=ltr_lo, ltr_hi=ltr_hi,
        features=features_ci,
        features_present=features_present
    )


def print_family_summary(family_key: str, avg: FamilyAverages, flagged_ids: Set[str] = None):
    print(f"\n[SUMMARY] {family_key}  (N={avg.n_elements})")
    print(f"  Average LTR-RT len: {avg.len_center} (95% CI: {avg.len_lo}-{avg.len_hi})")
    print(f"  Average LTR len:    {avg.ltr_center} (95% CI: {avg.ltr_lo}-{avg.ltr_hi})")
    for feat in avg.features:
        print(f"  Average {feat.name} start: {feat.start_center} (95% CI: {feat.start_lo}-{feat.start_hi})")
        print(f"  Average {feat.name} end:   {feat.end_center} (95% CI: {feat.end_lo}-{feat.end_hi})")
    if flagged_ids is not None:
        n_flagged = len(flagged_ids)
        if n_flagged > 0:
            print(f"  ** Potential false positives (outlier-long, lacking proteins): {n_flagged}")
            for eid in sorted(flagged_ids):
                print(f"     - {eid}")
        else:
            print(f"  No potential false positives detected.")


# -----------------------------
# Plotting: average with whiskers
# -----------------------------
def plot_family_average_with_whiskers(family_key: str, avg: FamilyAverages, out_path: str):
    legend_handles = build_legend_handles(avg.features_present)

    fig, ax = plt.subplots(figsize=(12.5, 3.2))
    L = max(1, avg.len_center)
    ltr = max(1, min(avg.ltr_center, L))

    center_feats = []
    for f in avg.features:
        s = max(1, min(f.start_center, L))
        e = max(1, min(f.end_center, L))
        if e < s:
            s, e = e, s
        center_feats.append(Feature(f.name, s, e))

    draw_element(ax, y=1.0, el_len=L, ltr_len=ltr,
                 proteins=sorted(center_feats, key=lambda x: x.start), label=None, height=0.55)

    draw_whisker(ax, y=0.55, lo=max(1, avg.len_lo), hi=max(1, avg.len_hi), color="black", lw=1.4)
    ax.text(L * 0.01, 0.40, "Element length 95% CI", fontsize=8, ha="left", va="center")

    ltr_lo = max(1, min(avg.ltr_lo, L))
    ltr_hi = max(1, min(avg.ltr_hi, L))
    draw_whisker(ax, y=1.33, lo=0, hi=ltr_hi, color=FEATURE_COLORS["LTR"], lw=1.2)
    draw_whisker(ax, y=1.33, lo=0, hi=ltr_lo, color=FEATURE_COLORS["LTR"], lw=1.2)
    draw_whisker(ax, y=1.33, lo=max(0, L - ltr_hi), hi=L, color=FEATURE_COLORS["LTR"], lw=1.2)
    draw_whisker(ax, y=1.33, lo=max(0, L - ltr_lo), hi=L, color=FEATURE_COLORS["LTR"], lw=1.2)
    ax.text(L * 0.5, 1.43, "LTR length 95% CI (both ends)", fontsize=8, ha="center", va="center")

    for f in avg.features:
        col = FEATURE_COLORS.get(f.name, DEFAULT_COLOR)
        draw_whisker(ax, y=1.18, lo=max(1, min(f.start_lo, L)), hi=max(1, min(f.start_hi, L)), color=col, lw=1.2)
        draw_whisker(ax, y=0.82, lo=max(1, min(f.end_lo, L)), hi=max(1, min(f.end_hi, L)), color=col, lw=1.2)

    ax.set_title(f"{family_key} : Average element structure (rare proteins filtered; whiskers=95% CI)")
    ax.set_xlim(0, max(avg.len_hi, L) * 1.02)
    ax.set_ylim(0.25, 1.65)
    ax.set_yticks([])
    ax.set_xlabel("Position [bp] (relative)")
    ax.legend(handles=legend_handles, loc="upper right", frameon=False, ncol=min(6, len(legend_handles)))
    plt.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


def plot_family_individual(family_key: str, fam_elements: List[Element], out_path: str,
                           flagged_ids: Set[str] = None):
    """Plot ALL individual elements (no max cap). Always shows all proteins.
    Elements in flagged_ids are highlighted with red outlines."""
    if flagged_ids is None:
        flagged_ids = set()

    n = len(fam_elements)
    max_len = max(e.length for e in fam_elements) if fam_elements else 1

    prots_seen = sorted(set(p.name for e in fam_elements for p in e.proteins))
    features_present = ["LTR"] + prots_seen
    has_flagged = any(e.element_id in flagged_ids for e in fam_elements)
    legend_handles = build_legend_handles(features_present, include_flagged=has_flagged)

    fig_h = max(3.0, 0.18 * n + 1.8)
    fig, ax = plt.subplots(figsize=(14, fig_h))

    for i, e in enumerate(fam_elements):
        y = n - i
        prots = sorted(e.proteins, key=lambda x: x.start)
        k2p_txt = "NA" if e.k2p is None else f"{e.k2p * 100:.1f}%"
        is_flagged = e.element_id in flagged_ids
        label = f"{e.element_id}  {k2p_txt}"
        if is_flagged:
            label += "  [FP?]"
        draw_element(ax, y=y, el_len=e.length, ltr_len=e.ltr_len, proteins=prots,
                     label=label, height=0.65, flagged=is_flagged)

    n_flagged = sum(1 for e in fam_elements if e.element_id in flagged_ids)
    title = f"{family_key} : Individual elements (n={n})"
    if n_flagged > 0:
        title += f"  [{n_flagged} potential FP highlighted]"
    ax.set_title(title)
    ax.set_xlim(0, max_len * 1.02)
    ax.set_ylim(0, n + 1)
    ax.set_yticks([])
    ax.set_xlabel("Position [bp] (relative)")
    ax.legend(handles=legend_handles, loc="lower right", frameon=False, ncol=min(6, len(legend_handles)))
    plt.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


def plot_all_elements_average(family_avg_map: Dict[str, FamilyAverages], out_path: str):
    fam_keys = sorted(family_avg_map.keys())

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

    ax.set_title("All elements : Average structure per family (rare proteins filtered)")
    ax.set_xlim(0, max_len_hi * 1.02)
    ax.set_ylim(0, len(fam_keys) + 1)
    ax.set_yticks([])
    ax.set_xlabel("Position [bp] (relative)")
    ax.legend(handles=legend_handles, loc="upper right", frameon=False, ncol=min(6, len(legend_handles)))
    plt.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


# -----------------------------
# FAI parsing
# -----------------------------
def parse_fai(fai_path: str) -> Dict[str, int]:
    chrom_lens: Dict[str, int] = {}
    print(f"[INFO] Reading FAI: {fai_path}")
    with open(fai_path, "r", encoding="utf-8") as f:
        for ln, line in enumerate(f, start=1):
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            chrom = parts[0]
            try:
                L = int(parts[1])
            except ValueError:
                continue
            chrom_lens[chrom] = L
    print(f"[INFO] FAI loaded: {len(chrom_lens)} contigs")
    return chrom_lens


# -----------------------------
# IGV-like HTML chromosome viewer
# -----------------------------
def generate_chrom_html(
    chrom: str,
    chrom_len: int,
    elements_on_chrom: List[Element],
    feature_colors: Dict[str, str],
    default_color: str,
) -> str:
    """Generate a self-contained IGV-like HTML viewer for one chromosome."""

    # Serialize element data to JSON for the JS viewer
    els_json = []
    for el in sorted(elements_on_chrom, key=lambda e: (e.start, e.end)):
        prots = []
        for p in el.proteins:
            prots.append({
                "name": p.name,
                "start": p.start,
                "end": p.end,
            })
        els_json.append({
            "id": el.element_id,
            "chrom": el.chrom,
            "start": el.start,
            "end": el.end,
            "length": el.length,
            "ltr_len": el.ltr_len,
            "superfamily": el.superfamily,
            "family": el.family,
            "strand": el.strand,
            "k2p": el.k2p,
            "proteins": prots,
        })

    data_js = json.dumps({
        "chrom": chrom,
        "chrom_len": chrom_len,
        "elements": els_json,
        "feature_colors": feature_colors,
        "default_color": default_color,
    })

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{chrom} — LTR-RT Genome Browser</title>
<style>
* {{ margin: 0; padding: 0; box-sizing: border-box; }}
body {{
    font-family: 'Segoe UI', -apple-system, BlinkMacSystemFont, sans-serif;
    background: #0d1117; color: #c9d1d9;
}}
#header {{
    background: #161b22; border-bottom: 1px solid #30363d;
    padding: 10px 18px; display: flex; align-items: center; gap: 16px; flex-wrap: wrap;
}}
#header h1 {{ font-size: 16px; color: #58a6ff; font-weight: 600; white-space: nowrap; }}
#header .info {{ font-size: 12px; color: #8b949e; }}
#nav {{
    background: #161b22; border-bottom: 1px solid #30363d;
    padding: 8px 18px; display: flex; align-items: center; gap: 10px; flex-wrap: wrap;
}}
#nav label {{ font-size: 12px; color: #8b949e; }}
#nav input[type=text] {{
    background: #0d1117; border: 1px solid #30363d; color: #c9d1d9;
    padding: 4px 8px; border-radius: 4px; font-size: 13px; width: 220px;
}}
#nav button {{
    background: #21262d; border: 1px solid #30363d; color: #c9d1d9;
    padding: 4px 10px; border-radius: 4px; cursor: pointer; font-size: 12px;
}}
#nav button:hover {{ background: #30363d; }}
#nav .sep {{ width: 1px; height: 20px; background: #30363d; }}
#zoom-info {{ font-size: 11px; color: #8b949e; min-width: 160px; }}
#overview-wrap {{
    background: #0d1117; border-bottom: 1px solid #30363d; padding: 6px 18px;
}}
#overview {{ width: 100%; height: 30px; cursor: pointer; display: block; }}
#main-canvas-wrap {{
    padding: 0 18px; overflow: hidden; position: relative;
}}
#main-canvas {{
    display: block; width: 100%; cursor: grab;
}}
#main-canvas.grabbing {{ cursor: grabbing; }}
#tooltip {{
    display: none; position: fixed; background: #1c2128; border: 1px solid #30363d;
    border-radius: 6px; padding: 10px 14px; font-size: 12px; color: #c9d1d9;
    pointer-events: none; z-index: 100; max-width: 380px; line-height: 1.5;
    box-shadow: 0 4px 12px rgba(0,0,0,0.4);
}}
#tooltip .tt-label {{ color: #8b949e; }}
#tooltip .tt-val {{ color: #f0f6fc; font-weight: 500; }}
#legend {{
    background: #161b22; border-top: 1px solid #30363d;
    padding: 8px 18px; display: flex; gap: 14px; flex-wrap: wrap; align-items: center;
}}
#legend .leg-item {{
    display: flex; align-items: center; gap: 5px; font-size: 11px; color: #c9d1d9;
}}
#legend .leg-swatch {{
    width: 14px; height: 14px; border: 1px solid #484f58; border-radius: 2px;
}}
#search-results {{
    font-size: 11px; color: #58a6ff;
}}
</style>
</head>
<body>

<div id="header">
    <h1>LTR-RT Browser</h1>
    <span class="info" id="chrom-info"></span>
</div>

<div id="nav">
    <label>Region:</label>
    <input type="text" id="region-input" placeholder="e.g. 1000000-2000000" />
    <button id="go-btn">Go</button>
    <div class="sep"></div>
    <button id="zoom-in">Zoom In</button>
    <button id="zoom-out">Zoom Out</button>
    <button id="zoom-all">Full View</button>
    <div class="sep"></div>
    <label>Search:</label>
    <input type="text" id="search-input" placeholder="family or element ID" />
    <button id="search-btn">Find</button>
    <span id="search-results"></span>
    <div style="flex:1"></div>
    <span id="zoom-info"></span>
</div>

<div id="overview-wrap">
    <canvas id="overview"></canvas>
</div>

<div id="main-canvas-wrap">
    <canvas id="main-canvas"></canvas>
</div>

<div id="legend"></div>

<div id="tooltip"></div>

<script>
const DATA = {data_js};

const chromLen = DATA.chrom_len;
const elements = DATA.elements;
const featureColors = DATA.feature_colors;
const defaultColor = DATA.default_color;

// State
let viewStart = 1;
let viewEnd = chromLen;
let canvasW = 0;
let canvasH = 0;
let trackH = 22;
let trackGap = 3;
let headerH = 26;
let rows = []; // assigned rows for non-overlapping layout
let hoveredEl = null;
let isDragging = false;
let dragStartX = 0;
let dragViewStart = 0;
let searchHighlight = null;

const overviewCanvas = document.getElementById('overview');
const overviewCtx = overviewCanvas.getContext('2d');
const mainCanvas = document.getElementById('main-canvas');
const mainCtx = mainCanvas.getContext('2d');
const tooltip = document.getElementById('tooltip');

function init() {{
    document.getElementById('chrom-info').textContent =
        DATA.chrom + '  |  ' + chromLen.toLocaleString() + ' bp  |  ' +
        elements.length + ' elements';

    buildLegend();
    assignRows();
    resize();

    window.addEventListener('resize', resize);
    mainCanvas.addEventListener('wheel', onWheel, {{ passive: false }});
    mainCanvas.addEventListener('mousedown', onMouseDown);
    mainCanvas.addEventListener('mousemove', onMouseMove);
    mainCanvas.addEventListener('mouseleave', () => {{ hoveredEl = null; tooltip.style.display = 'none'; }});
    window.addEventListener('mouseup', onMouseUp);
    overviewCanvas.addEventListener('click', onOverviewClick);

    document.getElementById('zoom-in').addEventListener('click', () => zoom(0.5));
    document.getElementById('zoom-out').addEventListener('click', () => zoom(2.0));
    document.getElementById('zoom-all').addEventListener('click', () => {{ viewStart = 1; viewEnd = chromLen; render(); }});
    document.getElementById('go-btn').addEventListener('click', goToRegion);
    document.getElementById('region-input').addEventListener('keydown', (e) => {{ if (e.key === 'Enter') goToRegion(); }});
    document.getElementById('search-btn').addEventListener('click', doSearch);
    document.getElementById('search-input').addEventListener('keydown', (e) => {{ if (e.key === 'Enter') doSearch(); }});

    render();
}}

function buildLegend() {{
    const leg = document.getElementById('legend');
    // Collect all protein names present
    const names = new Set(['LTR']);
    elements.forEach(el => el.proteins.forEach(p => names.add(p.name)));
    const ordered = ["LTR","GAG","PROT","RT","RH","INT","CH","CHD","CHDCR","ARH","ENDO"];
    const sorted_names = ordered.filter(n => names.has(n)).concat(
        [...names].filter(n => !ordered.includes(n)).sort()
    );
    sorted_names.forEach(n => {{
        const col = featureColors[n] || defaultColor;
        const item = document.createElement('div');
        item.className = 'leg-item';
        item.innerHTML = '<div class="leg-swatch" style="background:' + col + '"></div>' + n;
        leg.appendChild(item);
    }});
}}

function assignRows() {{
    // Interval scheduling: pack elements into non-overlapping rows
    elements.sort((a, b) => a.start - b.start || a.end - b.end);
    const rowEnds = []; // end coord of last element in each row
    rows = [];
    elements.forEach(el => {{
        let placed = false;
        for (let r = 0; r < rowEnds.length; r++) {{
            if (el.start > rowEnds[r] + Math.max(100, (el.end - el.start) * 0.02)) {{
                rowEnds[r] = el.end;
                rows.push(r);
                placed = true;
                break;
            }}
        }}
        if (!placed) {{
            rows.push(rowEnds.length);
            rowEnds.push(el.end);
        }}
    }});
}}

function resize() {{
    const dpr = window.devicePixelRatio || 1;
    const wrap = document.getElementById('main-canvas-wrap');
    canvasW = wrap.clientWidth;
    const numRows = rows.length > 0 ? Math.max(...rows) + 1 : 1;
    canvasH = headerH + numRows * (trackH + trackGap) + 20;
    canvasH = Math.max(canvasH, 120);
    // Cap visible height, allow scroll
    const maxVisH = Math.min(canvasH, window.innerHeight - 180);
    wrap.style.height = maxVisH + 'px';
    wrap.style.overflowY = canvasH > maxVisH ? 'auto' : 'hidden';

    mainCanvas.width = canvasW * dpr;
    mainCanvas.height = canvasH * dpr;
    mainCanvas.style.width = canvasW + 'px';
    mainCanvas.style.height = canvasH + 'px';
    mainCtx.setTransform(dpr, 0, 0, dpr, 0, 0);

    // Overview
    const ow = overviewCanvas.parentElement.clientWidth;
    overviewCanvas.width = ow * dpr;
    overviewCanvas.height = 30 * dpr;
    overviewCanvas.style.width = ow + 'px';
    overviewCanvas.style.height = '30px';
    overviewCtx.setTransform(dpr, 0, 0, dpr, 0, 0);

    render();
}}

function bpToX(bp) {{
    return (bp - viewStart) / (viewEnd - viewStart) * canvasW;
}}

function xToBp(x) {{
    return viewStart + x / canvasW * (viewEnd - viewStart);
}}

function render() {{
    renderOverview();
    renderMain();
    updateZoomInfo();
}}

function renderOverview() {{
    const w = overviewCanvas.width / (window.devicePixelRatio || 1);
    const h = 30;
    overviewCtx.clearRect(0, 0, w, h);

    // Background
    overviewCtx.fillStyle = '#161b22';
    overviewCtx.fillRect(0, 0, w, h);

    // Element density as ticks
    overviewCtx.fillStyle = '#58a6ff';
    elements.forEach(el => {{
        const x1 = (el.start / chromLen) * w;
        const x2 = Math.max(x1 + 1, (el.end / chromLen) * w);
        overviewCtx.fillRect(x1, 4, x2 - x1, h - 8);
    }});

    // Viewport indicator
    const vx1 = (viewStart / chromLen) * w;
    const vx2 = (viewEnd / chromLen) * w;
    overviewCtx.strokeStyle = '#f0883e';
    overviewCtx.lineWidth = 2;
    overviewCtx.strokeRect(vx1, 1, Math.max(2, vx2 - vx1), h - 2);
}}

function renderMain() {{
    mainCtx.clearRect(0, 0, canvasW, canvasH);

    // Background
    mainCtx.fillStyle = '#0d1117';
    mainCtx.fillRect(0, 0, canvasW, canvasH);

    // Ruler
    drawRuler();

    // Elements
    elements.forEach((el, idx) => {{
        const row = rows[idx];
        const x1 = bpToX(el.start);
        const x2 = bpToX(el.end);
        if (x2 < -50 || x1 > canvasW + 50) return; // off screen

        const y = headerH + row * (trackH + trackGap) + trackGap;
        const w = Math.max(2, x2 - x1);
        const th = trackH;

        const isHighlighted = (searchHighlight && (
            el.id === searchHighlight ||
            el.family.toLowerCase().includes(searchHighlight) ||
            el.superfamily.toLowerCase().includes(searchHighlight)
        ));
        const isHovered = (hoveredEl === idx);

        // Element outline
        mainCtx.strokeStyle = isHighlighted ? '#f0883e' : isHovered ? '#58a6ff' : '#30363d';
        mainCtx.lineWidth = isHighlighted ? 2 : isHovered ? 1.5 : 0.8;
        mainCtx.strokeRect(x1, y, w, th);

        // LTRs
        const ltrPx = Math.max(1, (el.ltr_len / (el.end - el.start)) * w);
        const ltrCol = featureColors['LTR'] || '#333';
        mainCtx.fillStyle = ltrCol;
        mainCtx.fillRect(x1, y, Math.min(ltrPx, w), th);
        mainCtx.fillRect(Math.max(x1, x2 - ltrPx), y, Math.min(ltrPx, w), th);

        // Proteins (relative coords mapped to genomic)
        el.proteins.forEach(p => {{
            const pStart = el.start + (p.start - 1);
            const pEnd = el.start + (p.end - 1);
            const px1 = bpToX(pStart);
            const px2 = bpToX(pEnd);
            const pw = Math.max(1, px2 - px1);
            mainCtx.fillStyle = featureColors[p.name] || defaultColor;
            mainCtx.fillRect(px1, y, pw, th);
        }});

        // Label if zoomed in enough
        if (w > 60) {{
            const label = el.superfamily + '/' + el.family;
            mainCtx.fillStyle = '#c9d1d9';
            mainCtx.font = '10px sans-serif';
            mainCtx.fillText(label, x1 + 3, y + th - 4, w - 6);
        }}
    }});
}}

function drawRuler() {{
    const span = viewEnd - viewStart;
    // Choose nice tick interval
    const raw = span / 8;
    const mag = Math.pow(10, Math.floor(Math.log10(raw)));
    let step = mag;
    if (raw / mag >= 5) step = mag * 5;
    else if (raw / mag >= 2) step = mag * 2;

    mainCtx.fillStyle = '#8b949e';
    mainCtx.strokeStyle = '#21262d';
    mainCtx.lineWidth = 1;
    mainCtx.font = '10px sans-serif';

    const first = Math.ceil(viewStart / step) * step;
    for (let bp = first; bp <= viewEnd; bp += step) {{
        const x = bpToX(bp);
        mainCtx.beginPath();
        mainCtx.moveTo(x, 0);
        mainCtx.lineTo(x, headerH - 4);
        mainCtx.stroke();

        let label;
        if (bp >= 1e6) label = (bp / 1e6).toFixed(bp % 1e6 === 0 ? 0 : 2) + ' Mb';
        else if (bp >= 1e3) label = (bp / 1e3).toFixed(bp % 1e3 === 0 ? 0 : 1) + ' kb';
        else label = bp.toString();
        mainCtx.fillText(label, x + 3, headerH - 8);
    }}

    // Baseline
    mainCtx.strokeStyle = '#30363d';
    mainCtx.beginPath();
    mainCtx.moveTo(0, headerH);
    mainCtx.lineTo(canvasW, headerH);
    mainCtx.stroke();
}}

function updateZoomInfo() {{
    const span = viewEnd - viewStart;
    let txt;
    if (span >= 1e6) txt = (span / 1e6).toFixed(2) + ' Mb';
    else if (span >= 1e3) txt = (span / 1e3).toFixed(1) + ' kb';
    else txt = span + ' bp';
    document.getElementById('zoom-info').textContent =
        viewStart.toLocaleString() + ' – ' + Math.round(viewEnd).toLocaleString() + '  (' + txt + ')';
    document.getElementById('region-input').placeholder =
        Math.round(viewStart) + '-' + Math.round(viewEnd);
}}

// --- Interaction ---
function zoom(factor, centerBp) {{
    const span = viewEnd - viewStart;
    if (!centerBp) centerBp = (viewStart + viewEnd) / 2;
    const frac = (centerBp - viewStart) / span;
    const newSpan = Math.max(500, Math.min(chromLen, span * factor));
    viewStart = Math.max(1, centerBp - frac * newSpan);
    viewEnd = Math.min(chromLen, viewStart + newSpan);
    viewStart = Math.max(1, viewEnd - newSpan);
    render();
}}

function onWheel(e) {{
    e.preventDefault();
    const rect = mainCanvas.getBoundingClientRect();
    const x = e.clientX - rect.left;
    const bp = xToBp(x);
    const factor = e.deltaY > 0 ? 1.3 : 1 / 1.3;
    zoom(factor, bp);
}}

function onMouseDown(e) {{
    isDragging = true;
    dragStartX = e.clientX;
    dragViewStart = viewStart;
    mainCanvas.classList.add('grabbing');
}}

function onMouseUp() {{
    isDragging = false;
    mainCanvas.classList.remove('grabbing');
}}

function onMouseMove(e) {{
    if (isDragging) {{
        const dx = e.clientX - dragStartX;
        const bpPerPx = (viewEnd - viewStart) / canvasW;
        let newStart = dragViewStart - dx * bpPerPx;
        const span = viewEnd - viewStart;
        newStart = Math.max(1, Math.min(chromLen - span, newStart));
        viewStart = newStart;
        viewEnd = newStart + span;
        render();
        return;
    }}

    // Hover detection
    const rect = mainCanvas.getBoundingClientRect();
    const mx = e.clientX - rect.left;
    const my = e.clientY - rect.top;
    let found = null;

    for (let idx = 0; idx < elements.length; idx++) {{
        const el = elements[idx];
        const row = rows[idx];
        const x1 = bpToX(el.start);
        const x2 = bpToX(el.end);
        const y = headerH + row * (trackH + trackGap) + trackGap;
        if (mx >= x1 && mx <= x2 && my >= y && my <= y + trackH) {{
            found = idx;
            break;
        }}
    }}

    if (found !== hoveredEl) {{
        hoveredEl = found;
        renderMain(); // re-render for hover highlight
    }}

    if (found !== null) {{
        const el = elements[found];
        const k2pTxt = el.k2p !== null ? (el.k2p * 100).toFixed(2) + '%' : 'NA';
        const protsHtml = el.proteins.map(p =>
            '<span style="color:' + (featureColors[p.name] || defaultColor) + ';font-weight:600">' +
            p.name + '</span> ' + p.start + '-' + p.end
        ).join(', ') || '<em>none</em>';

        tooltip.innerHTML =
            '<div><span class="tt-label">ID:</span> <span class="tt-val">' + el.id + '</span></div>' +
            '<div><span class="tt-label">Position:</span> <span class="tt-val">' + el.start.toLocaleString() + ' – ' + el.end.toLocaleString() + '</span></div>' +
            '<div><span class="tt-label">Length:</span> <span class="tt-val">' + el.length.toLocaleString() + ' bp</span></div>' +
            '<div><span class="tt-label">LTR:</span> <span class="tt-val">' + el.ltr_len + ' bp</span></div>' +
            '<div><span class="tt-label">Family:</span> <span class="tt-val">' + el.superfamily + '/' + el.family + '</span></div>' +
            '<div><span class="tt-label">K2P:</span> <span class="tt-val">' + k2pTxt + '</span></div>' +
            '<div><span class="tt-label">Domains:</span> ' + protsHtml + '</div>';
        tooltip.style.display = 'block';
        tooltip.style.left = Math.min(e.clientX + 14, window.innerWidth - 400) + 'px';
        tooltip.style.top = (e.clientY + 14) + 'px';
    }} else {{
        tooltip.style.display = 'none';
    }}
}}

function onOverviewClick(e) {{
    const rect = overviewCanvas.getBoundingClientRect();
    const x = e.clientX - rect.left;
    const w = rect.width;
    const bp = (x / w) * chromLen;
    const span = viewEnd - viewStart;
    viewStart = Math.max(1, Math.min(chromLen - span, bp - span / 2));
    viewEnd = viewStart + span;
    render();
}}

function goToRegion() {{
    const val = document.getElementById('region-input').value.trim();
    const m = val.match(/(\\d[\\d,]*)\\s*[-:]\\s*(\\d[\\d,]*)/);
    if (m) {{
        const s = parseInt(m[1].replace(/,/g, ''));
        const e = parseInt(m[2].replace(/,/g, ''));
        if (s > 0 && e > s) {{
            viewStart = Math.max(1, s);
            viewEnd = Math.min(chromLen, e);
            render();
        }}
    }}
}}

function doSearch() {{
    const q = document.getElementById('search-input').value.trim().toLowerCase();
    const res = document.getElementById('search-results');
    if (!q) {{ searchHighlight = null; res.textContent = ''; render(); return; }}

    searchHighlight = q;
    // Find first match and navigate
    const idx = elements.findIndex(el =>
        el.id.toLowerCase().includes(q) ||
        el.family.toLowerCase().includes(q) ||
        el.superfamily.toLowerCase().includes(q)
    );
    const count = elements.filter(el =>
        el.id.toLowerCase().includes(q) ||
        el.family.toLowerCase().includes(q) ||
        el.superfamily.toLowerCase().includes(q)
    ).length;

    if (idx >= 0) {{
        const el = elements[idx];
        const pad = Math.max(5000, (el.end - el.start) * 3);
        viewStart = Math.max(1, el.start - pad);
        viewEnd = Math.min(chromLen, el.end + pad);
        res.textContent = count + ' match' + (count > 1 ? 'es' : '') + ' found';
    }} else {{
        res.textContent = 'No matches';
    }}
    render();
}}

init();
</script>
</body>
</html>"""
    return html


# -----------------------------
# Master HTML with embedded plots
# -----------------------------
def generate_master_html(
    family_plots: Dict[str, Dict[str, str]],  # fk -> {"average": base64, "individual": base64}
    all_elements_b64: str,
    chrom_html_files: List[str],
    out_path: str,
):
    """Generate a single index.html with all plots embedded as base64 PNGs and links to chrom HTMLs."""

    sections = []

    # All elements average
    sections.append(f"""
    <div class="section">
        <h2>All Elements — Average Structure per Family</h2>
        <img src="data:image/png;base64,{all_elements_b64}" class="plot-img" />
    </div>""")

    # Per-family
    for fk in sorted(family_plots.keys()):
        plots = family_plots[fk]
        sections.append(f"""
    <div class="section">
        <h2>{fk}</h2>
        <h3>Average</h3>
        <img src="data:image/png;base64,{plots['average']}" class="plot-img" />
        <h3>Individual Elements</h3>
        <img src="data:image/png;base64,{plots['individual']}" class="plot-img" />
    </div>""")

    # Chrom links
    chrom_section = ""
    if chrom_html_files:
        links = "\n".join(
            f'        <li><a href="{os.path.relpath(cf, os.path.dirname(out_path))}">{os.path.basename(cf)}</a></li>'
            for cf in sorted(chrom_html_files)
        )
        chrom_section = f"""
    <div class="section">
        <h2>Per-Chromosome Genome Browsers</h2>
        <ul>
{links}
        </ul>
    </div>"""

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>LTR-RT Structure Plots</title>
<style>
body {{
    font-family: 'Segoe UI', -apple-system, sans-serif;
    background: #f6f8fa; color: #1f2328; max-width: 1200px; margin: 0 auto; padding: 24px;
}}
h1 {{ font-size: 22px; border-bottom: 2px solid #d0d7de; padding-bottom: 8px; margin-bottom: 20px; }}
h2 {{ font-size: 17px; color: #0969da; margin-top: 32px; margin-bottom: 8px; }}
h3 {{ font-size: 14px; color: #656d76; margin: 12px 0 4px; }}
.section {{ margin-bottom: 36px; }}
.plot-img {{ max-width: 100%; border: 1px solid #d0d7de; border-radius: 6px; margin: 4px 0 16px; }}
ul {{ list-style: none; padding: 0; }}
li {{ margin: 4px 0; }}
a {{ color: #0969da; text-decoration: none; }}
a:hover {{ text-decoration: underline; }}
</style>
</head>
<body>
<h1>LTR-RT Structure Report</h1>
{"".join(sections)}
{chrom_section}
</body>
</html>"""

    with open(out_path, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"[INFO] Master HTML written: {out_path}")


# -----------------------------
# Main
# -----------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", required=True, help="LTR-RT TSV file")
    ap.add_argument("--scn", required=True, help="SCN file")
    ap.add_argument("--gff", required=True, help="Domain HMM gff file")
    ap.add_argument("--out_dir", required=True, help="Output directory")

    ap.add_argument("--min_presence", type=float, default=0.50,
                    help="Min fraction of elements a protein must appear in to be included "
                         "in the AVERAGE plot only (default 0.50)")
    ap.add_argument("--boot", type=int, default=2000,
                    help="Bootstrap iterations for 95%% CI (default 2000)")
    ap.add_argument("--fai", default=None,
                    help="Optional FASTA .fai index to enable per-chromosome interactive HTML plots")
    ap.add_argument("--html", action="store_true",
                    help="Generate a master index.html with all outputs embedded; "
                         "otherwise only chromosome plots are HTML and family plots are PDF")
    ap.add_argument("--seed", type=int, default=1, help="Random seed for bootstrapping (default 1)")

    # --- False-positive detection parameters ---
    ap.add_argument("--no_fp", action="store_true",
                    help="Skip false-positive detection entirely; just plot LTR structure")
    ap.add_argument("--min_family_size", type=int, default=5,
                    help="Minimum number of elements in a family to perform outlier detection (default 5)")
    ap.add_argument("--fp_recovery_frac", type=float, default=0.50,
                    help="Fraction of expected proteins an outlier-long element must have "
                         "to be recovered from false-positive flagging (default 0.50)")
    ap.add_argument("--iqr_mul", type=float, default=0.5,
                    help="IQR multiplier for outlier fence: Q3 + iqr_mul*IQR (default 0.5). "
                         "Lower values tighten the fence and flag more elements; "
                         "higher values loosen it. Also scales the MAD multiplier for small families.")
    ap.add_argument("--min_ltr_aln", type=int, default=100,
                    help="Minimum LTR_LEN and ALN_LEN (bp) for dubious-clade quality filter "
                         "(default 120). Elements in mixture/unknown clades with LTR_LEN or "
                         "ALN_LEN below this value are flagged.")
    ap.add_argument("--dubious_k2p", type=float, default=0.15,
                    help="K2P divergence cutoff for dubious clades (mixture/unknown). "
                         "Elements above this threshold are flagged (default 0.15 = 15%%).")
    ap.add_argument("--recovery_k2p", type=float, default=0.10,
                    help="K2P divergence cutoff for single-protein recovery. When an outlier-"
                         "long element has exactly 1 protein and that barely meets the threshold, "
                         "it must have K2P below this value to be recovered (default 0.10 = 10%%).")
    ap.add_argument("--skip_unknown_iqr", action="store_true",
                    help="Skip IQR-based length outlier detection for 'unknown' clades. "
                         "If provided, unknown-clade elements are flagged solely based on "
                         "K2P and LTR_LEN/ALN_LEN quality filters.")

    args = ap.parse_args()
    random.seed(args.seed)

    os.makedirs(args.out_dir, exist_ok=True)

    elements = parse_tsv(args.tsv)
    shifts = parse_scn(args.scn)
    gff_feats = parse_gff_domains(args.gff)

    attach_shifts(elements, shifts)

    # --- Edit (1): auto-detect GFF coordinate mode ---
    coord_mode = detect_gff_coordinate_mode(elements, gff_feats)
    attach_proteins_and_strand(elements, gff_feats, coord_mode)

    # Strand normalize
    normalize_all_elements_strand(elements)

    # Group by superfamily/family
    family_to_elements: Dict[str, List[Element]] = {}
    for el in elements.values():
        fam_key = f"{el.superfamily}/{el.family}"
        family_to_elements.setdefault(fam_key, []).append(el)

    for fk in family_to_elements:
        family_to_elements[fk].sort(key=lambda e: e.length, reverse=True)

    print(f"[INFO] Grouped into {len(family_to_elements)} families")

    # Compute averages, detect false positives, print summaries, plot
    family_avg_map: Dict[str, FamilyAverages] = {}
    family_flagged: Dict[str, Set[str]] = {}
    family_plots_b64: Dict[str, Dict[str, str]] = {}  # for --html mode
    total_flagged = 0

    for fk, fam_elements in sorted(family_to_elements.items()):
        avg = compute_family_averages_with_ci(
            fam_elements=fam_elements,
            min_presence=args.min_presence,
            boot=args.boot,
        )
        family_avg_map[fk] = avg

        # --- False-positive detection ---
        if args.no_fp:
            flagged_ids: Set[str] = set()
        else:
            # Expected proteins are those that passed the min_presence filter in the average
            expected_proteins = [f.name for f in avg.features]
            flagged_ids = detect_false_positives(
                fam_elements=fam_elements,
                expected_proteins=expected_proteins,
                family_key=fk,
                min_family_size=args.min_family_size,
                protein_recovery_frac=args.fp_recovery_frac,
                iqr_mul=args.iqr_mul,
                min_ltr_aln=args.min_ltr_aln,
                dubious_k2p=args.dubious_k2p,
                recovery_k2p=args.recovery_k2p,
                skip_unknown_iqr=args.skip_unknown_iqr,
            )
        family_flagged[fk] = flagged_ids
        total_flagged += len(flagged_ids)

        # Excluded proteins info
        present_counts = {}
        for e in fam_elements:
            for p in set(pp.name for pp in e.proteins):
                present_counts[p] = present_counts.get(p, 0) + 1
        excluded = []
        for p, c in sorted(present_counts.items(), key=lambda x: (-x[1], x[0])):
            if p not in [ff.name for ff in avg.features]:
                frac = c / len(fam_elements)
                excluded.append(f"{p} ({c}/{len(fam_elements)}={frac:.3f})")
        print_family_summary(fk, avg, flagged_ids)
        if excluded:
            print(f"  Excluded rare proteins (<{args.min_presence:.2f} presence): " +
                  ", ".join(excluded[:12]) + (" ..." if len(excluded) > 12 else ""))

        safe_fk = safe_name(fk)

        # Average plot
        out_avg = os.path.join(args.out_dir, f"{safe_fk}_average.pdf")
        plot_family_average_with_whiskers(fk, avg, out_avg)

        # Individual plot (shows ALL proteins, no max cap) — with FP highlighting
        out_ind = os.path.join(args.out_dir, f"{safe_fk}_individual.pdf")
        plot_family_individual(fk, fam_elements, out_ind, flagged_ids=flagged_ids)

        # If --html, also generate base64 PNGs
        if args.html:
            # Average as PNG
            fig_avg, ax_avg = plt.subplots(figsize=(12.5, 3.2))
            L = max(1, avg.len_center)
            ltr = max(1, min(avg.ltr_center, L))
            center_feats = []
            for f in avg.features:
                s = max(1, min(f.start_center, L))
                e = max(1, min(f.end_center, L))
                if e < s: s, e = e, s
                center_feats.append(Feature(f.name, s, e))
            draw_element(ax_avg, y=1.0, el_len=L, ltr_len=ltr,
                         proteins=sorted(center_feats, key=lambda x: x.start), height=0.55)
            draw_whisker(ax_avg, y=0.55, lo=max(1, avg.len_lo), hi=max(1, avg.len_hi), color="black", lw=1.4)
            for f in avg.features:
                col = FEATURE_COLORS.get(f.name, DEFAULT_COLOR)
                draw_whisker(ax_avg, y=1.18, lo=max(1, min(f.start_lo, L)), hi=max(1, min(f.start_hi, L)), color=col)
                draw_whisker(ax_avg, y=0.82, lo=max(1, min(f.end_lo, L)), hi=max(1, min(f.end_hi, L)), color=col)
            ax_avg.set_title(f"{fk} : Average (n={avg.n_elements})")
            ax_avg.set_xlim(0, max(avg.len_hi, L)*1.02)
            ax_avg.set_ylim(0.25, 1.65)
            ax_avg.set_yticks([])
            ax_avg.set_xlabel("Position [bp]")
            ax_avg.legend(handles=build_legend_handles(avg.features_present), loc="upper right", frameon=False, ncol=6)
            plt.tight_layout()
            avg_b64 = fig_to_base64_png(fig_avg)
            plt.close(fig_avg)

            # Individual as PNG — with FP highlighting
            n = len(fam_elements)
            max_len = max(e.length for e in fam_elements) if fam_elements else 1
            prots_seen = sorted(set(p.name for e in fam_elements for p in e.proteins))
            fig_h = max(3.0, 0.18 * n + 1.8)
            fig_ind, ax_ind = plt.subplots(figsize=(14, fig_h))
            has_flagged = any(e.element_id in flagged_ids for e in fam_elements)
            for i, e in enumerate(fam_elements):
                y = n - i
                prots = sorted(e.proteins, key=lambda x: x.start)
                k2p_txt = "NA" if e.k2p is None else f"{e.k2p*100:.1f}%"
                is_flagged = e.element_id in flagged_ids
                label = f"{e.element_id}  {k2p_txt}"
                if is_flagged:
                    label += "  [FP?]"
                draw_element(ax_ind, y=y, el_len=e.length, ltr_len=e.ltr_len, proteins=prots,
                             label=label, height=0.65, flagged=is_flagged)
            title = f"{fk} : Individual (n={n})"
            n_flagged_fam = sum(1 for e in fam_elements if e.element_id in flagged_ids)
            if n_flagged_fam > 0:
                title += f"  [{n_flagged_fam} potential FP]"
            ax_ind.set_title(title)
            ax_ind.set_xlim(0, max_len*1.02)
            ax_ind.set_ylim(0, n+1)
            ax_ind.set_yticks([])
            ax_ind.set_xlabel("Position [bp]")
            ax_ind.legend(handles=build_legend_handles(["LTR"]+prots_seen, include_flagged=has_flagged),
                          loc="lower right", frameon=False, ncol=6)
            plt.tight_layout()
            ind_b64 = fig_to_base64_png(fig_ind)
            plt.close(fig_ind)

            family_plots_b64[fk] = {"average": avg_b64, "individual": ind_b64}

    # Global all-elements plot
    out_all = os.path.join(args.out_dir, "all_elements.pdf")
    plot_all_elements_average(family_avg_map, out_all)

    all_b64 = ""
    if args.html:
        # Regenerate as PNG for embedding
        fam_keys = sorted(family_avg_map.keys())
        all_prots = set()
        max_lh = 1
        for fk in fam_keys:
            a = family_avg_map[fk]
            max_lh = max(max_lh, a.len_hi)
            for feat in a.features:
                all_prots.add(feat.name)
        fig_h = max(3.0, 0.35*len(fam_keys)+2.0)
        fig_all, ax_all = plt.subplots(figsize=(14, fig_h))
        for i, fk in enumerate(fam_keys):
            a = family_avg_map[fk]
            y = len(fam_keys)-i
            L = max(1, a.len_center)
            ltr = max(1, min(a.ltr_center, L))
            feats = []
            for f in a.features:
                s = max(1, min(f.start_center, L))
                e = max(1, min(f.end_center, L))
                if e < s: s, e = e, s
                feats.append(Feature(f.name, s, e))
            draw_element(ax_all, y=y, el_len=L, ltr_len=ltr,
                         proteins=sorted(feats, key=lambda x: x.start), label=fk, height=0.65)
        ax_all.set_title("All elements: Average structure per family")
        ax_all.set_xlim(0, max_lh*1.02)
        ax_all.set_ylim(0, len(fam_keys)+1)
        ax_all.set_yticks([])
        ax_all.set_xlabel("Position [bp]")
        ax_all.legend(handles=build_legend_handles(["LTR"]+sorted(all_prots)),
                      loc="upper right", frameon=False, ncol=6)
        plt.tight_layout()
        all_b64 = fig_to_base64_png(fig_all)
        plt.close(fig_all)

    # Per-chromosome interactive HTML plots (always HTML, edit 4)
    chrom_html_files = []
    if args.fai:
        chrom_lens = parse_fai(args.fai)

        chrom_to_elements: Dict[str, List[Element]] = {}
        for el in elements.values():
            chrom_to_elements.setdefault(el.chrom, []).append(el)

        chrom_outdir = os.path.join(args.out_dir, "chrom_plots")
        os.makedirs(chrom_outdir, exist_ok=True)

        for chrom, clen in chrom_lens.items():
            out_html = os.path.join(chrom_outdir, f"{safe_name(chrom)}.html")
            els = chrom_to_elements.get(chrom, [])
            html_content = generate_chrom_html(
                chrom=chrom,
                chrom_len=clen,
                elements_on_chrom=els,
                feature_colors=FEATURE_COLORS,
                default_color=DEFAULT_COLOR,
            )
            with open(out_html, "w", encoding="utf-8") as fh:
                fh.write(html_content)
            chrom_html_files.append(out_html)
            print(f"[INFO] Chromosome HTML: {out_html} ({len(els)} elements)")

        missing_contigs = sorted([c for c in chrom_to_elements if c not in chrom_lens])
        if missing_contigs:
            print(f"[WARN] {len(missing_contigs)} contigs with elements not in FAI. First 10:")
            for c in missing_contigs[:10]:
                print(f"       - {c}")

    # Master HTML (edit 5)
    if args.html:
        master_path = os.path.join(args.out_dir, "index.html")
        generate_master_html(
            family_plots=family_plots_b64,
            all_elements_b64=all_b64,
            chrom_html_files=chrom_html_files,
            out_path=master_path,
        )

    # Write flagged elements to a TSV and filtered input TSV (excluding FPs)
    all_flagged_ids: Set[str] = set()
    if not args.no_fp:
        for fids in family_flagged.values():
            all_flagged_ids.update(fids)

        flagged_tsv_path = os.path.join(args.out_dir, "flagged_false_positives.tsv")
        with open(flagged_tsv_path, "w", encoding="utf-8") as fout:
            fout.write("#element_id\tfamily\tlength\tltr_len_raw\taln_len_raw\tk2p\t"
                       "num_proteins\texpected_proteins\treason\n")
            for fk in sorted(family_flagged.keys()):
                flagged_ids = family_flagged[fk]
                if not flagged_ids:
                    continue
                expected_proteins = [f.name for f in family_avg_map[fk].features]
                dubious = _is_dubious_clade(fk)
                unknown = _is_unknown_clade(fk)

                # Recompute upper fence for reason annotation
                fam_els = family_to_elements[fk]
                n_fam = len(fam_els)
                skip_iqr_here = unknown and args.skip_unknown_iqr
                if n_fam >= args.min_family_size and not skip_iqr_here:
                    lengths = sorted(float(e.length) for e in fam_els)
                    upper_fence = _compute_upper_fence(lengths, n_fam, iqr_mul=args.iqr_mul)
                else:
                    upper_fence = float("inf")

                n_expected = len(expected_proteins)

                for el in fam_els:
                    if el.element_id not in flagged_ids:
                        continue
                    el_prots = sorted(set(p.name for p in el.proteins))
                    k2p_str = f"{el.k2p:.4f}" if el.k2p is not None else "NA"
                    reasons = []

                    # Check dubious-clade reasons
                    if dubious:
                        if el.k2p is not None and el.k2p > args.dubious_k2p:
                            reasons.append(f"dubious_clade;k2p>{args.dubious_k2p*100:.0f}%")
                        if el.ltr_len_raw < args.min_ltr_aln:
                            reasons.append(f"dubious_clade;LTR_LEN<{args.min_ltr_aln}")
                        if el.aln_len_raw < args.min_ltr_aln:
                            reasons.append(f"dubious_clade;ALN_LEN<{args.min_ltr_aln}")

                    # Check length-outlier reason
                    if el.length > upper_fence:
                        if n_expected == 0:
                            reasons.append("outlier_long;non_autonomous")
                        else:
                            n_hits = sum(1 for ep in expected_proteins if ep in el_prots)
                            reasons.append(f"outlier_long;proteins={n_hits}/{n_expected}")

                    reason_str = "|".join(reasons) if reasons else "flagged"
                    fout.write(f"{el.element_id}\t{fk}\t{el.length}\t"
                               f"{el.ltr_len_raw}\t{el.aln_len_raw}\t{k2p_str}\t"
                               f"{len(el_prots)}\t{','.join(expected_proteins) or 'none'}\t"
                               f"{reason_str}\n")
        print(f"[INFO] Flagged false positives written: {flagged_tsv_path} ({total_flagged} elements)")

    # Write filtered TSV (input TSV minus flagged elements)
    if not args.no_fp:
        filtered_tsv_path = os.path.join(args.out_dir, "filtered.tsv")
        n_kept = 0
        n_removed = 0
        with open(args.tsv, "r", encoding="utf-8") as fin, \
             open(filtered_tsv_path, "w", encoding="utf-8") as fout:
            for line in fin:
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    fout.write(line)
                    continue
                parts = stripped.split("\t")
                if len(parts) < 3:
                    fout.write(line)
                    continue
                raw_id = parts[0]
                if "#" in raw_id:
                    elem_part = raw_id.split("#", 1)[0]
                else:
                    elem_part = raw_id
                if elem_part in all_flagged_ids:
                    n_removed += 1
                else:
                    fout.write(line)
                    n_kept += 1
        print(f"[INFO] Filtered TSV written: {filtered_tsv_path} "
              f"(kept={n_kept}, removed={n_removed})")

    print(f"\n[INFO] Done.")
    print(f"[INFO] Output directory: {args.out_dir}")
    print(f"[INFO] Average plots exclude proteins with presence < {args.min_presence:.2f} "
          f"and show 95% bootstrap CI whiskers.")
    if args.no_fp:
        print(f"[INFO] False-positive detection: DISABLED (--no_fp)")
    else:
        print(f"[INFO] False-positive detection: IQR (N>=15) / MAD (N<15) + gap analysis, "
              f"iqr_mul={args.iqr_mul:.2f}, min_family_size={args.min_family_size}, "
              f"protein_recovery_frac={args.fp_recovery_frac:.2f}, "
              f"recovery_k2p={args.recovery_k2p:.0%}, "
              f"dubious_k2p={args.dubious_k2p:.0%}, min_ltr_aln={args.min_ltr_aln}")
        if args.skip_unknown_iqr:
            print(f"[INFO] IQR skipped for 'unknown' clades (--skip_unknown_iqr)")
        print(f"[INFO] Total potential false positives flagged: {total_flagged}")
    if args.html:
        print(f"[INFO] Master HTML: {os.path.join(args.out_dir, 'index.html')}")
    if chrom_html_files:
        print(f"[INFO] Chromosome browsers: {os.path.join(args.out_dir, 'chrom_plots')}/")


if __name__ == "__main__":
    main()
