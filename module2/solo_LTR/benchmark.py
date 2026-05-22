#!/usr/bin/env python3
"""benchmark.py - score predicted solo-LTR loci against a PrinTE ground truth BED.

The ground-truth labels mirror PrinTE / solo_intact_count.py:

  clean_solo    : feature_id has _SOLO, no _FRAG, te_class == LTR,
                  AND is NOT disrupted (no CUT_BY and not split by a
                  nested insertion within +/-100 records).
  intact_ltr    : feature_id has neither _SOLO nor _FRAG, te_class == LTR,
                  not disrupted.
  frag_ltr      : feature_id has _FRAG, te_class == LTR.
  disrupted_solo: feature_id has _SOLO, te_class == LTR, IS disrupted.
  other         : anything else (genes, non-LTR TEs, etc).

A predicted interval is bucketed by the BEST overlap it has with any
ground-truth interval (reciprocal-overlap rule, --min-overlap-frac).
A clean-solo truth is detected if ANY prediction meets the rule.

Reports F1 against clean_solo as the headline metric, plus breakdowns
that explain where the false positives went (intact? frag? novel?).
"""
import argparse
import os
import re
import sys
from collections import defaultdict, deque

NEIGHBORHOOD = 100  # records, mirrors solo_intact_count.py


# --------------------------------------------------------------------------- #
# Ground-truth parsing (byte-identical-in-spirit to solo_intact_count.py)
# --------------------------------------------------------------------------- #
class Rec:
    __slots__ = ("chrom", "start", "end", "name", "tsd", "strand",
                 "feature_id", "supp")

    def __init__(self, chrom, start, end, name, tsd, strand):
        self.chrom, self.start, self.end = chrom, start, end
        self.name, self.tsd, self.strand = name, tsd, strand
        if ";" in name:
            parts = name.split(";")
            self.feature_id, self.supp = parts[0], parts[1:]
        else:
            self.feature_id, self.supp = name, []


def te_class_super(fid):
    """'te_name#TE_class/TE_super[~junk]' -> (class, super)."""
    if "#" not in fid:
        return None, None
    try:
        _, rest = fid.split("#", 1)
        cls, sup = rest.split("/", 1)
        return cls, sup.split("~")[0]
    except ValueError:
        return None, None


def is_disrupted(rec, window, center):
    if not rec.supp:
        return False
    if "CUT_BY" in rec.supp[0]:
        return True
    lo = max(0, center - NEIGHBORHOOD)
    hi = min(len(window), center + NEIGHBORHOOD + 1)
    for j in range(lo, hi):
        if j == center:
            continue
        other = window[j]
        if other.feature_id.startswith("gene"):
            continue
        if other.tsd == rec.tsd and other.strand == rec.strand:
            if rec.name.startswith(other.name) or other.name.startswith(rec.name):
                return True
    return False


def classify(rec, window, center):
    fid = rec.feature_id
    if fid.startswith("gene"):
        return "gene"
    if "_FRAG" in fid:
        cls, _ = te_class_super(fid)
        return "frag_ltr" if cls == "LTR" else "frag_other"
    cls, _ = te_class_super(fid)
    if cls != "LTR":
        return "non_ltr"
    disrupted = is_disrupted(rec, window, center)
    if "_SOLO" in fid:
        return "disrupted_solo" if disrupted else "clean_solo"
    return "fragmented_ltr" if disrupted else "intact_ltr"


def stream_truth(path):
    """Stream truth BED with PrinTE classification. Yields (Rec, category)."""
    buf = deque()
    cur = 0
    src = open(path)
    line_iter = iter(src)
    primed = False
    pending = []  # (Rec, center idx within buf at classify time)
    try:
        while True:
            rec = None
            for raw in line_iter:
                if not raw.strip() or raw.startswith("#"):
                    continue
                parts = raw.rstrip("\n").split("\t")
                if len(parts) < 6:
                    continue
                rec = Rec(parts[0], int(parts[1]), int(parts[2]),
                          parts[3], parts[4], parts[5])
                break
            if rec is not None:
                buf.append(rec)
            if not primed:
                if rec is not None and len(buf) <= NEIGHBORHOOD:
                    continue
                primed = True
            if cur < len(buf):
                cat = classify(buf[cur], buf, cur)
                yield buf[cur], cat
                cur += 1
            while cur > NEIGHBORHOOD and len(buf) > 2 * NEIGHBORHOOD + 1:
                buf.popleft(); cur -= 1
            if rec is None and cur >= len(buf):
                break
    finally:
        src.close()


def load_truth(path):
    """Return dict cat -> list of (chrom, start, end)."""
    by_cat = defaultdict(list)
    for rec, cat in stream_truth(path):
        by_cat[cat].append((rec.chrom, rec.start, rec.end))
    for k in by_cat:
        by_cat[k].sort()
    return by_cat


def load_pred(path):
    out = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            out.append((f[0], int(f[1]), int(f[2])))
    out.sort()
    return out


# --------------------------------------------------------------------------- #
# Interval overlap utilities
# --------------------------------------------------------------------------- #
def index_by_chrom(intervals):
    """chrom -> sorted list of (start, end). Sorting kept stable."""
    by = defaultdict(list)
    for c, s, e in intervals:
        by[c].append((s, e))
    for c in by:
        by[c].sort()
    return by


def best_overlap(pred, truth_idx):
    """Return (best_chrom_truth_idx_pos, ovl_frac_min) for the prediction.

    ovl_frac_min := min( ovl / len_pred , ovl / len_truth ) -- the
    reciprocal-overlap measure used by --min-overlap-frac.
    """
    chrom, ps, pe = pred
    plen = max(1, pe - ps)
    chrom_list = truth_idx.get(chrom)
    if not chrom_list:
        return None, 0.0

    # Linear scan -- list is small per chrom (10**4-ish even in mammalian
    # genomes), and we are not in a tight loop.
    best_idx, best_score = None, 0.0
    for i, (ts, te) in enumerate(chrom_list):
        if te <= ps:
            continue
        if ts >= pe:
            break
        ovl = min(pe, te) - max(ps, ts)
        if ovl <= 0:
            continue
        score = min(ovl / plen, ovl / max(1, te - ts))
        if score > best_score:
            best_score = score
            best_idx = i
    return best_idx, best_score


# --------------------------------------------------------------------------- #
# Reporting
# --------------------------------------------------------------------------- #
PRIORITY = ["clean_solo", "intact_ltr", "fragmented_ltr",
            "frag_ltr", "disrupted_solo", "non_ltr", "frag_other", "gene"]


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--truth", required=True, help="PrinTE BED (gen<N>_final.bed)")
    p.add_argument("--pred", required=True, help="Predicted solo-LTR BED")
    p.add_argument("--min-overlap-frac", type=float, default=0.5,
                   help="Reciprocal overlap fraction required for a match")
    p.add_argument("--tsv", default=None,
                   help="Optional: append a 1-line TSV summary to this file")
    p.add_argument("--tag", default="",
                   help="Free-form tag included in the TSV summary line")
    p.add_argument("--verbose", action="store_true")
    args = p.parse_args()

    sys.stderr.write(f"[benchmark] loading truth from {args.truth}\n")
    truth = load_truth(args.truth)
    truth_idx = {cat: index_by_chrom(v) for cat, v in truth.items()}

    n_clean_solo = len(truth.get("clean_solo", []))
    n_intact = len(truth.get("intact_ltr", []))
    n_disrupted_solo = len(truth.get("disrupted_solo", []))
    n_frag_ltr = len(truth.get("frag_ltr", []))
    sys.stderr.write(
        f"[benchmark] truth counts: clean_solo={n_clean_solo}, "
        f"intact_ltr={n_intact}, frag_ltr={n_frag_ltr}, "
        f"disrupted_solo={n_disrupted_solo}\n")

    preds = load_pred(args.pred)
    n_pred = len(preds)
    sys.stderr.write(f"[benchmark] predictions: {n_pred}\n")

    # 1. Predictions -> bucket by best matching ground-truth category.
    pred_buckets = defaultdict(int)         # cat -> #predictions
    pred_unmatched = 0
    cats_present = [c for c in PRIORITY if truth.get(c)]
    for pred in preds:
        best_cat, best_score = None, 0.0
        for cat in cats_present:
            _, s = best_overlap(pred, truth_idx[cat])
            if s >= args.min_overlap_frac and s > best_score:
                best_cat, best_score = cat, s
        if best_cat is None:
            pred_unmatched += 1
        else:
            pred_buckets[best_cat] += 1

    # 2. Detected clean_solo truths: any prediction with reciprocal overlap.
    detected_clean = 0
    if "clean_solo" in truth_idx:
        pred_idx = index_by_chrom(preds)
        for chrom, s, e in truth.get("clean_solo", []):
            tlen = max(1, e - s)
            chrom_list = pred_idx.get(chrom, [])
            hit = False
            for ps, pe in chrom_list:
                if pe <= s: continue
                if ps >= e: break
                ovl = min(pe, e) - max(ps, s)
                if ovl <= 0: continue
                score = min(ovl / tlen, ovl / max(1, pe - ps))
                if score >= args.min_overlap_frac:
                    hit = True
                    break
            if hit:
                detected_clean += 1

    # 3. Headline metrics: clean_solo vs all (intact/frag/disrupted/other = FP).
    TP = pred_buckets.get("clean_solo", 0)
    FP = n_pred - TP
    FN = n_clean_solo - detected_clean

    precision = TP / n_pred if n_pred else 0.0
    recall = detected_clean / n_clean_solo if n_clean_solo else 0.0
    f1 = (2 * precision * recall / (precision + recall)
          if (precision + recall) > 0 else 0.0)

    # 4. Print headline + breakdown.
    print(f"=== solo-LTR benchmark ===")
    print(f"truth: clean_solo={n_clean_solo}  intact={n_intact}  "
          f"frag_ltr={n_frag_ltr}  disrupted_solo={n_disrupted_solo}")
    print(f"predictions: {n_pred}")
    print(f"detected clean solos: {detected_clean} / {n_clean_solo}")
    print(f"TP={TP}  FP={FP}  FN={FN}")
    print(f"precision={precision:.4f}  recall={recall:.4f}  F1={f1:.4f}")
    print(f"prediction breakdown (best-overlap category):")
    for cat in PRIORITY:
        if cat in pred_buckets:
            print(f"  {cat:<18} {pred_buckets[cat]}")
    if pred_unmatched:
        print(f"  {'(no overlap)':<18} {pred_unmatched}")

    if args.tsv:
        new = not os.path.exists(args.tsv)
        with open(args.tsv, "a") as fh:
            if new:
                fh.write("tag\tn_pred\tdetected_clean\tn_clean_solo\t"
                         "TP\tFP\tFN\tprecision\trecall\tF1\t"
                         "fp_intact\tfp_frag_ltr\tfp_disrupted_solo\t"
                         "fp_no_overlap\tmin_overlap_frac\n")
            fh.write("\t".join(str(x) for x in [
                args.tag, n_pred, detected_clean, n_clean_solo,
                TP, FP, FN,
                f"{precision:.6f}", f"{recall:.6f}", f"{f1:.6f}",
                pred_buckets.get("intact_ltr", 0),
                pred_buckets.get("frag_ltr", 0) + pred_buckets.get("fragmented_ltr", 0),
                pred_buckets.get("disrupted_solo", 0),
                pred_unmatched, args.min_overlap_frac,
            ]) + "\n")


if __name__ == "__main__":
    main()
