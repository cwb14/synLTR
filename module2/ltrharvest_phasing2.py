#!/usr/bin/env python3
"""
LTR-only SubPhaser-like subgenome phasing (with K2P integration + binned density PDF)
+ genomic windows exchange scan + optional chromosome painting with FAI.

python ./test7.2.py --ltr_fasta Oalta.ltr.lib.fa --homoeolog_config sg.config --outdir ltr_phaser_out -k 15 --hash_size 100000000 --dump_min 3 --ratio_f 2 --q_min 200 --threads 8 --boot_reps 1000 --k2p_file  <(cat *_dedup) --k2p_col 11 --k2p_bin 0.0005 --genome_fai Oalta.fa.fai

Adds (existing from your version):
  - --k2p_file: TSV where col1 is LTR id (exact header string) and K2P is column 11 (1-based)
  - Adds "k2p" column to ltr_enrichment_and_age.tsv
  - Outputs binned density distributions (NOT KDE) of K2P for:
      (1) post_polyploid
      (2) pre_polyploid
    using bin size 0.0005 by default.
  - New output: k2p_density_pre_vs_post.pdf

Adds (requested now):
  - Window scan (default 1 Mb, non-overlapping) over each chromosome:
      Aggregates subgenome-specific k-mer hits across all LTRs in each window.
      Tests enrichment for each SG in each window vs the rest of the SAME chromosome (Fisher greater),
      BH-corrects across all windows*SG,
      and flags windows whose best enriched SG != chromosome SG assignment as candidate exchanges/assembly issues.
  - Outputs:
      window_enrichment.tsv
      window_exchange_candidates.tsv (subset where mismatch and significant)
  - Optional chromosome painting if --genome_fai provided:
      chromosome_painting_exchanges.pdf

Original outputs still produced:
  - chromosome_phasing.tsv
  - differential_kmers.txt
  - subgenome_specific_kmers.<SG>.txt
  - ltr_enrichment_and_age.tsv
  - pca_chromosomes.pdf
  - heatmap_differential_kmers.pdf
"""

import argparse
import os
import re
import sys
import shutil
import tempfile
import subprocess
import math
from collections import defaultdict, Counter
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import scipy.sparse as sp
from scipy.stats import ttest_ind, fisher_exact
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
from scipy.optimize import linear_sum_assignment

from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

from statsmodels.stats.multitest import multipletests


# -------------------------
# FASTA parsing
# -------------------------
def fasta_iter(path):
    header = None
    seq_chunks = []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if header is not None:
            yield header, "".join(seq_chunks)


def parse_ltr_header(h):
    """
    Expected:
      CHR:START-END#CLASS/SUPERFAM/FAM
    Return dict with chrom,start,end,klass,superfamily,family.
    """
    chrom_part, rest = (h.split("#", 1) + [""])[:2]
    chrom = chrom_part.split(":", 1)[0]
    start, end = (None, None)
    m = re.search(r":(\d+)-(\d+)", chrom_part)
    if m:
        start, end = int(m.group(1)), int(m.group(2))

    klass, superfam, fam = ("", "", "")
    if rest:
        parts = rest.split("/")
        if len(parts) >= 1:
            klass = parts[0]
        if len(parts) >= 2:
            superfam = parts[1]
        if len(parts) >= 3:
            fam = parts[2]
    return {
        "chrom": chrom,
        "start": start,
        "end": end,
        "klass": klass,
        "superfamily": superfam,
        "family": fam,
    }


# -------------------------
# Config parsing
# -------------------------
def read_homoeolog_config(path):
    """
    Lines: tokens separated by whitespace. '#' comments allowed.
    Returns list of lists, each list is a homoeolog set of chromosome names.
    """
    sets = []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.split("#", 1)[0].strip()
            if not line:
                continue
            toks = line.split()
            if toks:
                sets.append(toks)
    return sets


def infer_n_subgenomes(hsets):
    # Use the modal size among non-singleton sets, else max size, else 1.
    sizes = [len(s) for s in hsets if len(s) >= 2]
    if not sizes:
        return 1
    c = Counter(sizes)
    modal = c.most_common(1)[0][0]
    return modal


# -------------------------
# K2P parsing (optional)
# -------------------------
def read_k2p_map(path, k2p_col_1based=11, sep=None):
    """
    Returns: dict {ltr_id: k2p_float}
    Keeps first occurrence; warns on malformed lines.
    """
    idx = k2p_col_1based - 1
    d = {}
    bad = 0
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for ln, raw in enumerate(f, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split() if sep is None else line.split(sep)
            if len(parts) <= idx:
                bad += 1
                continue
            ltr_id = parts[0]
            try:
                k2p = float(parts[idx])
            except ValueError:
                bad += 1
                continue
            if ltr_id not in d:
                d[ltr_id] = k2p
    if bad:
        print(f"[WARN] K2P: skipped {bad} malformed lines", file=sys.stderr)
    return d


# -------------------------
# FAI parsing (optional)
# -------------------------
def read_fai_lengths(path):
    """
    FASTA index (.fai): col1=seqname col2=length ...
    Returns dict {chrom: length_int}
    """
    d = {}
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            if not raw.strip() or raw.startswith("#"):
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            name = parts[0]
            try:
                length = int(parts[1])
            except ValueError:
                continue
            d[name] = length
    return d


# -------------------------
# Jellyfish helpers
# -------------------------
def run(cmd, check=True):
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if check and p.returncode != 0:
        raise RuntimeError(
            f"Command failed ({p.returncode}): {' '.join(cmd)}\nSTDERR:\n{p.stderr}"
        )
    return p


def jellyfish_count_and_dump(fasta_path, out_prefix, k, hash_size, dump_min, canonical=True):
    """
    Creates:
      out_prefix.jf
      out_prefix.dump.txt  (kmer count per line: "KMER COUNT")
    """
    jf_path = out_prefix + ".jf"
    dump_path = out_prefix + ".dump.txt"

    cmd = ["jellyfish", "count", "-m", str(k), "-s", str(hash_size), "-o", jf_path]
    if canonical:
        cmd.append("--canonical")
    cmd.append(fasta_path)
    run(cmd)

    cmd = ["jellyfish", "dump", "-c", "-L", str(dump_min), jf_path]
    p = run(cmd)
    with open(dump_path, "w", encoding="utf-8") as out:
        out.write(p.stdout)

    return dump_path


def parse_jellyfish_dump(dump_path):
    """
    dump -c format: "kmer count" per line
    Returns dict {kmer: count}, and total_count (sum counts).
    """
    d = {}
    total = 0
    with open(dump_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            kmer, cnt = line.split()
            c = int(cnt)
            d[kmer] = c
            total += c
    return d, total


# -------------------------
# Matrix building & differential kmers
# -------------------------
def build_sparse_count_matrix(kmer_counts_by_chrom, chroms):
    """
    kmer_counts_by_chrom: dict chrom -> dict{kmer:count}
    chroms: ordered list of chrom names
    Returns:
      counts_csr (m x n), kmer_index_to_seq (list of kmers), col_totals (n,), row_totals (m,)
    """
    kmer_to_row = {}
    rows = []
    cols = []
    data = []

    for j, chrom in enumerate(chroms):
        d = kmer_counts_by_chrom.get(chrom, {})
        for kmer, cnt in d.items():
            i = kmer_to_row.get(kmer)
            if i is None:
                i = len(kmer_to_row)
                kmer_to_row[kmer] = i
            rows.append(i)
            cols.append(j)
            data.append(cnt)

    m = len(kmer_to_row)
    n = len(chroms)
    counts = sp.coo_matrix((data, (rows, cols)), shape=(m, n), dtype=np.int32).tocsr()
    kmer_index_to_seq = [None] * m
    for kmer, i in kmer_to_row.items():
        kmer_index_to_seq[i] = kmer

    col_totals = np.asarray(counts.sum(axis=0)).ravel().astype(np.float64)
    row_totals = np.asarray(counts.sum(axis=1)).ravel().astype(np.int64)
    return counts, kmer_index_to_seq, col_totals, row_totals


def compute_ratio_matrix(counts_csr, col_totals):
    """
    M0: ratios per kmer per chromosome (counts / total kmers of chrom).
    Returns float32 CSC for fast column slicing.
    """
    inv = np.zeros_like(col_totals, dtype=np.float64)
    nz = col_totals > 0
    inv[nz] = 1.0 / col_totals[nz]

    M0 = counts_csr.astype(np.float32).tocsc(copy=True)
    for j in range(M0.shape[1]):
        if inv[j] != 0:
            start, end = M0.indptr[j], M0.indptr[j + 1]
            M0.data[start:end] *= inv[j]
    return M0


def differential_kmer_mask(M0_csc, row_totals, chroms, homoeolog_sets, ratio_f=2.0, q_min=200):
    """
    For each homoeolog set, sort ratios high->low for each kmer: r0, r1...
    f = r0/r1 >= ratio_f must be true in ALL non-singleton sets
    AND q (row_totals) >= q_min
    """
    m, n = M0_csc.shape
    ok = (row_totals >= q_min)

    chrom_to_idx = {c: i for i, c in enumerate(chroms)}

    for hset in homoeolog_sets:
        if len(hset) < 2:
            continue
        idxs = [chrom_to_idx[x] for x in hset if x in chrom_to_idx]
        if len(idxs) < 2:
            continue

        sub = M0_csc[:, idxs].toarray().astype(np.float32)  # (m, g)

        if sub.shape[1] == 2:
            r0 = np.maximum(sub[:, 0], sub[:, 1])
            r1 = np.minimum(sub[:, 0], sub[:, 1])
        else:
            part = np.partition(sub, kth=sub.shape[1] - 2, axis=1)
            r0 = part[:, -1]
            r1 = part[:, -2]

        f = np.zeros_like(r0, dtype=np.float32)
        mask_r1 = (r1 > 0)
        f[mask_r1] = r0[mask_r1] / r1[mask_r1]
        f[~mask_r1] = np.where(r0[~mask_r1] > 0, np.inf, 0.0)

        ok &= (f >= ratio_f)

    return ok


def zscale_rows_dense(X):
    means = X.mean(axis=1, keepdims=True)
    stds = X.std(axis=1, keepdims=True)
    stds[stds == 0] = 1.0
    return (X - means) / stds


# -------------------------
# Clustering with bootstrap stability
# -------------------------
def hungarian_relabel(base_labels, new_labels, n_clusters):
    cm = np.zeros((n_clusters, n_clusters), dtype=np.int64)
    for a, b in zip(base_labels, new_labels):
        cm[a, b] += 1
    cost = cm.max() - cm
    r, c = linear_sum_assignment(cost)
    mapping = {c_i: r_i for r_i, c_i in zip(r, c)}
    relabeled = np.array([mapping.get(x, x) for x in new_labels], dtype=np.int32)
    return relabeled


def kmeans_with_bootstrap(M1, n_clusters, seed, boot_frac, boot_reps):
    rng = np.random.default_rng(seed)
    n_kmers, n_chrom = M1.shape
    X = M1.T  # (n_chrom, n_kmers)

    km = KMeans(n_clusters=n_clusters, n_init="auto", random_state=seed)
    base_labels = km.fit_predict(X).astype(np.int32)

    same_counts = np.zeros(n_chrom, dtype=np.int64)
    sample_size = max(1, int(round(boot_frac * n_kmers)))

    for _ in range(boot_reps):
        idx = rng.choice(n_kmers, size=sample_size, replace=False)
        Xb = M1[idx, :].T  # (n_chrom, sample_size)
        km_b = KMeans(
            n_clusters=n_clusters,
            n_init="auto",
            random_state=int(rng.integers(0, 2**31 - 1)),
        )
        lb = km_b.fit_predict(Xb).astype(np.int32)
        lb2 = hungarian_relabel(base_labels, lb, n_clusters)
        same_counts += (lb2 == base_labels)

    stability = (same_counts / float(boot_reps)) * 100.0
    return base_labels, stability


# -------------------------
# Visualization
# -------------------------
def save_pca_pdf(M1, chroms, labels, out_pdf):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    X = M1.T
    pca = PCA(n_components=2, random_state=0)
    coords = pca.fit_transform(X)

    plt.figure()
    for i, c in enumerate(chroms):
        plt.scatter(coords[i, 0], coords[i, 1], s=40)
        plt.text(coords[i, 0], coords[i, 1], f" {c} (SG{labels[i]+1})", fontsize=8, va="center")
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    plt.title("Chromosome PCA on differential k-mers (LTR-only)")
    plt.tight_layout()
    plt.savefig(out_pdf)
    plt.close()


def save_heatmap_pdf(M1, chroms, out_pdf, max_kmers=10000, seed=0):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rng = np.random.default_rng(seed)
    m, n = M1.shape
    take = min(max_kmers, m)
    idx = rng.choice(m, size=take, replace=False) if take < m else np.arange(m)
    X = M1[idx, :]

    dcol = pdist(X.T, metric="euclidean")
    lcol = linkage(dcol, method="average")
    col_order = leaves_list(lcol)

    drow = pdist(X, metric="euclidean")
    lrow = linkage(drow, method="average")
    row_order = leaves_list(lrow)

    Xo = X[row_order, :][:, col_order]
    chroms_o = [chroms[i] for i in col_order]

    plt.figure(figsize=(max(6, n * 0.5), 8))
    plt.imshow(Xo, aspect="auto", interpolation="nearest")
    plt.xticks(np.arange(n), chroms_o, rotation=90, fontsize=8)
    plt.yticks([])
    plt.title(f"Heatmap (Z-scaled) of {take} sampled differential k-mers")
    plt.tight_layout()
    plt.savefig(out_pdf)
    plt.close()


def save_k2p_binned_density_pdf(k2p_pre, k2p_post, out_pdf, bin_size=0.0005, title=None):
    """
    Density bins (histogram density=True) for pre vs post.
    Not KDE, not smoothing.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    k2p_pre = np.asarray(k2p_pre, dtype=float)
    k2p_post = np.asarray(k2p_post, dtype=float)

    # Remove NaNs
    k2p_pre = k2p_pre[np.isfinite(k2p_pre)]
    k2p_post = k2p_post[np.isfinite(k2p_post)]

    if k2p_pre.size == 0 and k2p_post.size == 0:
        raise RuntimeError("No finite K2P values for either pre_polyploid or post_polyploid.")

    xmin = np.inf
    xmax = -np.inf
    for arr in (k2p_pre, k2p_post):
        if arr.size:
            xmin = min(xmin, float(arr.min()))
            xmax = max(xmax, float(arr.max()))
    if not np.isfinite(xmin) or not np.isfinite(xmax) or xmin == xmax:
        xmin = 0.0 if not np.isfinite(xmin) else xmin - bin_size
        xmax = xmin + 10 * bin_size

    # Build bins aligned to bin_size
    left = math.floor(xmin / bin_size) * bin_size
    right = math.ceil(xmax / bin_size) * bin_size
    if right <= left:
        right = left + bin_size

    nbins = int(round((right - left) / bin_size))
    nbins = max(1, nbins)
    bins = left + bin_size * np.arange(nbins + 1, dtype=float)

    plt.figure(figsize=(9, 5))
    if k2p_post.size:
        plt.hist(
            k2p_post,
            bins=bins,
            density=True,
            histtype="step",
            linewidth=2,
            label=f"post_polyploid (n={k2p_post.size})",
        )
    if k2p_pre.size:
        plt.hist(
            k2p_pre,
            bins=bins,
            density=True,
            histtype="step",
            linewidth=2,
            label=f"pre_polyploid (n={k2p_pre.size})",
        )

    plt.xlabel(f"K2P (bin size {bin_size})")
    plt.ylabel("Density")
    plt.title(title or "K2P binned density: pre_polyploid vs post_polyploid")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_pdf)
    plt.close()


def save_chromosome_painting_pdf(chrom_lengths, chrom_order, chrom_sg_map, windows, n_sg, out_pdf,
                                window_size, title="Chromosome painting (SG + exchanges)"):
    """
    chrom_lengths: dict chrom->len
    chrom_order: list chrom
    chrom_sg_map: dict chrom->sg (0-based) or None
    windows: list of dicts with keys:
        chrom, start, end, best_sg (or None), best_q, mismatch_bool, significant_bool
    Paints each chromosome as a horizontal bar partitioned into windows.
      - default color: chromosome SG
      - exchange windows: color = best_sg
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    # Use tab10/20 deterministically
    cmap = plt.get_cmap("tab10" if n_sg <= 10 else "tab20")

    # Index windows by chrom and sort
    by_chrom = defaultdict(list)
    for w in windows:
        by_chrom[w["chrom"]].append(w)
    for c in by_chrom:
        by_chrom[c].sort(key=lambda x: (x["start"], x["end"]))

    # Layout
    fig_h = max(4.0, 0.35 * len(chrom_order) + 1.5)
    plt.figure(figsize=(12, fig_h))
    ax = plt.gca()

    y = 0
    yticks = []
    ylabels = []

    maxlen = max(chrom_lengths.get(c, 0) for c in chrom_order) if chrom_order else 1

    for chrom in chrom_order:
        L = chrom_lengths.get(chrom, 0)
        if L <= 0:
            continue

        # Baseline bar (light gray)
        ax.add_patch(Rectangle((0, y - 0.18), L, 0.36, facecolor="0.9", edgecolor="0.8", linewidth=0.5))

        base_sg = chrom_sg_map.get(chrom, None)

        # Paint windows
        for w in by_chrom.get(chrom, []):
            ws, we = w["start"], w["end"]
            if we <= ws:
                continue
            # choose color
            color_sg = None
            if w["significant"] and w["best_sg"] is not None:
                # if mismatch window: color by best_sg (exchange)
                if w["mismatch"]:
                    color_sg = w["best_sg"]
                else:
                    # match window: paint by base SG if known, else best_sg
                    color_sg = base_sg if base_sg is not None else w["best_sg"]
            else:
                # not significant: paint by base SG if known, else leave gray
                color_sg = base_sg

            if color_sg is None:
                continue

            fc = cmap(int(color_sg) % cmap.N)
            ax.add_patch(Rectangle((ws, y - 0.18), we - ws, 0.36, facecolor=fc, edgecolor=None, linewidth=0))

            # Outline exchange windows to pop visually
            if w["significant"] and w["mismatch"]:
                ax.add_patch(Rectangle((ws, y - 0.18), we - ws, 0.36, fill=False, edgecolor="k", linewidth=0.7))

        yticks.append(y)
        ylabels.append(chrom)
        y += 1

    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=8)
    ax.set_xlim(0, maxlen)
    ax.set_xlabel("Genomic position (bp)")
    ax.set_title(f"{title}\nWindow size: {window_size} bp")
    ax.grid(False)

    # Legend
    handles = []
    labels = []
    for sg in range(n_sg):
        handles.append(Rectangle((0, 0), 1, 1, facecolor=cmap(sg % cmap.N), edgecolor="none"))
        labels.append(f"SG{sg+1}")
    ax.legend(handles, labels, loc="upper right", ncol=min(6, n_sg), frameon=False, fontsize=8)

    plt.tight_layout()
    plt.savefig(out_pdf)
    plt.close()


# -------------------------
# Subgenome-specific kmers & LTR enrichment
# -------------------------
def one_tailed_ttest_greater(a, b):
    t, p2 = ttest_ind(a, b, equal_var=False, nan_policy="omit")
    if np.isnan(t):
        return np.nan
    if t > 0:
        return p2 / 2.0
    else:
        return 1.0 - (p2 / 2.0)


def call_subgenome_specific_kmers(M0_csc, diff_idx, chroms, chrom_labels, p_thresh=0.05, test="ttest"):
    n_sg = int(chrom_labels.max() + 1)
    chrom_by_sg = [[] for _ in range(n_sg)]
    for j, lab in enumerate(chrom_labels):
        chrom_by_sg[int(lab)].append(j)

    sg_specific = {sg: set() for sg in range(n_sg)}
    sub = M0_csc[diff_idx, :].toarray().astype(np.float32)

    all_cols = np.arange(len(chroms), dtype=int)

    for i_local, i_global in enumerate(diff_idx):
        row = sub[i_local, :]
        for sg in range(n_sg):
            A_idx = np.array(chrom_by_sg[sg], dtype=int)
            if A_idx.size == 0:
                continue
            B_mask = np.ones(len(chroms), dtype=bool)
            B_mask[A_idx] = False
            B_idx = all_cols[B_mask]

            A = row[A_idx]
            B = row[B_idx]
            if len(A) < 2 or len(B) < 2:
                continue

            p = one_tailed_ttest_greater(A, B)
            if p is not None and not np.isnan(p) and p < p_thresh:
                sg_specific[sg].add(int(i_global))

    return sg_specific


def ltr_kmer_hits(seq, k, kmer_sets_by_sg):
    seq = seq.upper()
    total = max(0, len(seq) - k + 1)
    hits = [0] * len(kmer_sets_by_sg)
    if total == 0:
        return hits, 0
    for i in range(total):
        kmer = seq[i:i + k]
        for sg, sset in enumerate(kmer_sets_by_sg):
            if kmer in sset:
                hits[sg] += 1
    return hits, total


def fisher_one_tailed_greater(c00, c01, c10, c11):
    table = np.array([[c00, c01], [c10, c11]], dtype=np.int64)
    try:
        _, p = fisher_exact(table, alternative="greater")
    except TypeError:
        odds, p2 = fisher_exact(table)
        if odds > 1:
            p = p2 / 2.0
        else:
            p = 1.0 - (p2 / 2.0)
    return p


# -------------------------
# Window enrichment scan (new)
# -------------------------
def infer_chrom_lengths_from_ltrs(ltrs):
    """
    If no FAI, infer chromosome length as max(end) among LTR headers (or max(start)+1).
    """
    mx = defaultdict(int)
    for _header, _seq, meta in ltrs:
        c = meta["chrom"]
        e = meta["end"]
        s = meta["start"]
        if e is not None:
            mx[c] = max(mx[c], int(e))
        elif s is not None:
            mx[c] = max(mx[c], int(s) + 1)
    # If ends are 1-based inclusive, the "length" here is approximate.
    # We just need a window tiling coordinate system; using max(end) is fine.
    return dict(mx)


def build_window_assignments(ltrs, window_size):
    """
    Returns dict chrom -> list of (window_index, ltr_index)
    LTR is assigned to window based on midpoint coordinate.
    """
    by_chrom = defaultdict(list)
    for i, (_h, _seq, meta) in enumerate(ltrs):
        c = meta["chrom"]
        s = meta["start"]
        e = meta["end"]
        if s is None or e is None:
            continue
        mid = (int(s) + int(e)) // 2
        widx = int(mid // window_size)
        by_chrom[c].append((widx, i))
    return by_chrom


def compute_window_enrichment(per_ltr_counts, ltrs, chrom_to_sg, n_sg,
                              window_size, enrich_q=0.05):
    """
    Enrichment for each window vs rest of chromosome (Fisher greater), BH across all windows*SG.

    per_ltr_counts[i] = (hits_by_sg(list len n_sg), total_k)
    ltrs[i] = (header, seq, meta)
    chrom_to_sg: dict chrom->sg (0-based)
    Returns:
      windows: list of dict with per-window summary
    """
    # Precompute chromosome totals for background (sum of LTRs on that chrom)
    chrom_hits = defaultdict(lambda: np.zeros(n_sg, dtype=np.int64))
    chrom_total = defaultdict(int)
    chrom_ltr_n = defaultdict(int)

    for i, (_h, _seq, meta) in enumerate(ltrs):
        c = meta["chrom"]
        hits, tot = per_ltr_counts[i]
        chrom_total[c] += int(tot)
        chrom_ltr_n[c] += 1
        for sg in range(n_sg):
            chrom_hits[c][sg] += int(hits[sg])

    # Assign LTRs to windows by midpoint
    by_chrom_w = build_window_assignments(ltrs, window_size)

    # Build raw window aggregates
    window_rows = []  # each: dict, with pvals placeholder
    tests = []        # (win_row_index, sg, pval)
    for chrom, pairs in by_chrom_w.items():
        if not pairs:
            continue
        # aggregate by window index
        agg_hits = defaultdict(lambda: np.zeros(n_sg, dtype=np.int64))
        agg_tot = defaultdict(int)
        agg_n = defaultdict(int)

        for widx, ltr_i in pairs:
            hits, tot = per_ltr_counts[ltr_i]
            agg_tot[widx] += int(tot)
            agg_n[widx] += 1
            for sg in range(n_sg):
                agg_hits[widx][sg] += int(hits[sg])

        for widx in sorted(agg_tot.keys()):
            win_tot = agg_tot[widx]
            win_n = agg_n[widx]
            win_hits = agg_hits[widx]
            if win_tot <= 0:
                continue

            start = int(widx * window_size)
            end = int(start + window_size)

            row = {
                "chrom": chrom,
                "start": start,
                "end": end,
                "ltr_count": int(win_n),
                "total_kmers": int(win_tot),
                "hits_by_sg": win_hits.copy(),
                "chrom_sg": chrom_to_sg.get(chrom, None),
                "best_sg": None,
                "best_q": 1.0,
                "enriched_sgs": [],
                "significant": False,
                "mismatch": False,
            }
            row_idx = len(window_rows)
            window_rows.append(row)

            # Enrichment vs rest of same chromosome
            chr_tot = int(chrom_total.get(chrom, 0))
            chr_hits = chrom_hits.get(chrom, np.zeros(n_sg, dtype=np.int64))
            rest_tot = chr_tot - win_tot
            if chr_tot <= 0 or rest_tot <= 0:
                # can't test; keep p=1 for all
                for sg in range(n_sg):
                    tests.append((row_idx, sg, 1.0))
                continue

            for sg in range(n_sg):
                c00 = int(win_hits[sg])
                c01 = int(win_tot - win_hits[sg])
                c10 = int(chr_hits[sg] - win_hits[sg])
                c11 = int((chr_tot - win_tot) - (chr_hits[sg] - win_hits[sg]))
                if (c00 + c01 + c10 + c11) <= 0:
                    p = 1.0
                else:
                    p = fisher_one_tailed_greater(c00, c01, c10, c11)
                tests.append((row_idx, sg, p))

    if not tests:
        return []

    pvals = np.array([t[2] for t in tests], dtype=np.float64)
    _, qvals, _, _ = multipletests(pvals, alpha=enrich_q, method="fdr_bh")

    # Attach qvals
    q_by_win = defaultdict(dict)
    for (win_i, sg, _p), q in zip(tests, qvals):
        q_by_win[win_i][sg] = float(q)

    for wi, row in enumerate(window_rows):
        qs = [q_by_win[wi].get(sg, 1.0) for sg in range(n_sg)]
        enriched = [sg for sg in range(n_sg) if qs[sg] <= enrich_q]
        row["enriched_sgs"] = enriched

        best_sg = None
        best_q = 1.0
        if enriched:
            # tie-break by (q, -hits)
            best_sg = min(enriched, key=lambda sg: (qs[sg], -int(row["hits_by_sg"][sg])))
            best_q = qs[best_sg]

        row["best_sg"] = best_sg
        row["best_q"] = float(best_q)
        row["significant"] = (best_sg is not None and best_q <= enrich_q)

        chrom_sg = row["chrom_sg"]
        row["mismatch"] = bool(row["significant"] and chrom_sg is not None and best_sg != chrom_sg)

    return window_rows


# -------------------------
# Main
# -------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ltr_fasta", required=True, help="Intact LTR-RT FASTA")
    ap.add_argument("--homoeolog_config", required=True, help="Homoeolog chromosome sets file")
    ap.add_argument("--outdir", required=True)

    ap.add_argument("--k2p_file", default="", help="TSV/space file with col1=LTR id and K2P in column 11 (1-based)")
    ap.add_argument("--k2p_col", type=int, default=11, help="1-based K2P column in --k2p_file (default 11)")
    ap.add_argument("--k2p_bin", type=float, default=0.0005, help="Bin size for K2P binned density (default 0.0005)")

    # Window exchange scan + painting
    ap.add_argument("--window_size", type=int, default=1000000, help="Window size in bp (default 1000000)")
    ap.add_argument("--window_step", type=int, default=0,
                    help="Window step in bp (default 0 => equals window_size, i.e. non-overlapping)")
    ap.add_argument("--genome_fai", default="", help="Optional genome .fai for chromosome lengths & painting")
    ap.add_argument("--paint_pdf", default="", help="Optional output PDF name (default chromosome_painting_exchanges.pdf)")
    ap.add_argument("--exchange_q", type=float, default=0.05, help="BH q-value threshold for window exchange calls (default 0.05)")
    ap.add_argument("--min_ltrs_per_window", type=int, default=1, help="Min LTR count per window to report (default 1)")

    ap.add_argument("-k", type=int, default=15, help="k-mer size (default 15)")
    ap.add_argument("--hash_size", default="100000000", help="jellyfish -s (default 100000000)")
    ap.add_argument("--dump_min", type=int, default=3, help="jellyfish dump -L (default 3)")
    ap.add_argument("--ratio_f", type=float, default=2.0, help="Differential k-mer r0/r1 threshold (default 2)")
    ap.add_argument("--q_min", type=int, default=200, help="Differential k-mer total genome count q (default 200)")
    ap.add_argument("--n_subgenomes", type=int, default=0, help="Override N subgenomes; else inferred from config")
    ap.add_argument("--boot_frac", type=float, default=0.5, help="Bootstrap fraction of kmers (default 0.5)")
    ap.add_argument("--boot_reps", type=int, default=1000, help="Bootstrap replicates (default 1000)")
    ap.add_argument("--threads", type=int, default=4, help="Parallel jellyfish jobs (default 4)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--specific_p", type=float, default=0.05, help="p-value for sg-specific kmers (default 0.05)")
    ap.add_argument("--enrich_p", type=float, default=0.05, help="BH q-value for LTR enrichment (default 0.05)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    hsets = read_homoeolog_config(args.homoeolog_config)
    N = args.n_subgenomes if args.n_subgenomes > 0 else infer_n_subgenomes(hsets)

    if args.window_step <= 0:
        args.window_step = args.window_size
    if args.window_step != args.window_size:
        print("[WARN] This implementation reports non-overlapping windows based on window_size; "
              "window_step != window_size is accepted but used only for painting tiling. "
              "If you truly need sliding-window stats, say so and I’ll adapt it.",
              file=sys.stderr)

    # K2P map (optional)
    k2p_map = {}
    if args.k2p_file:
        k2p_map = read_k2p_map(args.k2p_file, k2p_col_1based=args.k2p_col)
        print(f"[INFO] Loaded K2P for {len(k2p_map)} LTRs from {args.k2p_file}", file=sys.stderr)

    tmp = tempfile.mkdtemp(prefix="ltr_subphaser_")
    try:
        chrom_fastas = {}
        chrom_to_handle = {}
        chrom_seen = set()

        def get_handle(chrom):
            if chrom not in chrom_to_handle:
                path = os.path.join(tmp, f"{chrom}.fa")
                chrom_fastas[chrom] = path
                chrom_to_handle[chrom] = open(path, "w", encoding="utf-8")
            return chrom_to_handle[chrom]

        ltrs = []  # (header, seq, meta)
        for header, seq in fasta_iter(args.ltr_fasta):
            meta = parse_ltr_header(header)
            chrom = meta["chrom"]
            h = get_handle(chrom)
            h.write(f">{header}\n")
            for i in range(0, len(seq), 80):
                h.write(seq[i:i + 80] + "\n")
            chrom_seen.add(chrom)
            ltrs.append((header, seq, meta))

        for h in chrom_to_handle.values():
            h.close()

        chroms = sorted(chrom_seen)
        if len(chroms) < 2:
            raise RuntimeError("Need >=2 chromosomes (from headers) to phase.")

        # Jellyfish per chromosome (parallel)
        kmer_counts_by_chrom = {}
        col_totals = {}

        def do_one(chrom):
            fa = chrom_fastas[chrom]
            out_prefix = os.path.join(tmp, f"{chrom}.k{args.k}")
            dump = jellyfish_count_and_dump(
                fasta_path=fa,
                out_prefix=out_prefix,
                k=args.k,
                hash_size=args.hash_size,
                dump_min=args.dump_min,
                canonical=True,
            )
            d, total = parse_jellyfish_dump(dump)
            return chrom, d, total

        with ThreadPoolExecutor(max_workers=max(1, args.threads)) as ex:
            futs = [ex.submit(do_one, c) for c in chroms]
            for fu in as_completed(futs):
                chrom, d, total = fu.result()
                kmer_counts_by_chrom[chrom] = d
                col_totals[chrom] = total

        # Build count matrix, M0 ratios
        counts, kmer_list, col_sums, row_sums = build_sparse_count_matrix(kmer_counts_by_chrom, chroms)
        M0 = compute_ratio_matrix(counts, col_sums)

        # Differential kmers
        mask = differential_kmer_mask(M0, row_sums, chroms, hsets, ratio_f=args.ratio_f, q_min=args.q_min)
        diff_idx = np.where(mask)[0].astype(np.int64)
        if diff_idx.size == 0:
            raise RuntimeError(
                "No differential k-mers found. Consider lowering --q_min or --ratio_f, or --dump_min."
            )

        diff_kmers_path = os.path.join(args.outdir, "differential_kmers.txt")
        with open(diff_kmers_path, "w", encoding="utf-8") as out:
            for i in diff_idx:
                out.write(kmer_list[int(i)] + "\n")

        # M1 (Z-scaled) on differential kmers only
        M0_diff = M0[diff_idx, :].toarray().astype(np.float32)
        M1 = zscale_rows_dense(M0_diff).astype(np.float32)

        # KMeans + bootstrap
        labels, stability = kmeans_with_bootstrap(
            M1, n_clusters=N, seed=args.seed, boot_frac=args.boot_frac, boot_reps=args.boot_reps
        )

        # Outputs: chromosome phasing
        chrom_tsv = os.path.join(args.outdir, "chromosome_phasing.tsv")
        with open(chrom_tsv, "w", encoding="utf-8") as out:
            out.write("chromosome\tsubgenome\tbootstrap_percent\tltr_kmer_total\n")
            for j, c in enumerate(chroms):
                out.write(f"{c}\tSG{int(labels[j])+1}\t{stability[j]:.2f}\t{int(col_sums[j])}\n")

        # PCA + heatmap
        save_pca_pdf(M1, chroms, labels, os.path.join(args.outdir, "pca_chromosomes.pdf"))
        save_heatmap_pdf(M1, chroms, os.path.join(args.outdir, "heatmap_differential_kmers.pdf"), seed=args.seed)

        # Subgenome-specific kmers
        sg_specific_idx = call_subgenome_specific_kmers(
            M0_csc=M0,
            diff_idx=diff_idx,
            chroms=chroms,
            chrom_labels=labels,
            p_thresh=args.specific_p,
            test="ttest",
        )

        sg_kmer_sets = []
        for sg in range(N):
            kmers = set(kmer_list[i] for i in sg_specific_idx.get(sg, set()))
            sg_kmer_sets.append(kmers)
            outp = os.path.join(args.outdir, f"subgenome_specific_kmers.SG{sg+1}.txt")
            with open(outp, "w", encoding="utf-8") as out:
                for kmer in sorted(kmers):
                    out.write(kmer + "\n")

        # Enrichment per LTR (Fisher one-tailed + BH)
        chrom_to_sg = {c: int(labels[i]) for i, c in enumerate(chroms)}

        all_tests = []          # (ltr_i, sg, pval)
        per_ltr_counts = []     # (hits_by_sg, total_k)
        for i, (header, seq, meta) in enumerate(ltrs):
            hits, total_k = ltr_kmer_hits(seq, args.k, sg_kmer_sets)
            per_ltr_counts.append((hits, total_k))

        bg_hits_by_sg = np.zeros(N, dtype=np.int64)
        bg_total = 0
        for hits, total_k in per_ltr_counts:
            bg_total += total_k
            for sg in range(N):
                bg_hits_by_sg[sg] += hits[sg]

        for i, (hits, total_k) in enumerate(per_ltr_counts):
            for sg in range(N):
                c00 = hits[sg]
                c01 = total_k - hits[sg]
                c10 = int(bg_hits_by_sg[sg] - hits[sg])
                c11 = int((bg_total - total_k) - (bg_hits_by_sg[sg] - hits[sg]))
                if total_k == 0 or (c00 + c01 + c10 + c11) <= 0:
                    p = 1.0
                else:
                    p = fisher_one_tailed_greater(c00, c01, c10, c11)
                all_tests.append((i, sg, p))

        pvals = np.array([x[2] for x in all_tests], dtype=np.float64)
        _, qvals, _, _ = multipletests(pvals, alpha=args.enrich_p, method="fdr_bh")

        # Collate per LTR results (+K2P)
        ltr_out = os.path.join(args.outdir, "ltr_enrichment_and_age.tsv")
        with open(ltr_out, "w", encoding="utf-8") as out:
            out.write(
                "ltr_id\tchrom\tstart\tend\tclass\tsuperfamily\tfamily\tchrom_subgenome\t"
                "enriched_subgenomes\tbest_enriched_subgenome\tbest_q\t"
                "hits_by_sg\ttotal_kmers\tcategory\tk2p\n"
            )

            q_by = defaultdict(dict)
            for (i, sg, _p), q in zip(all_tests, qvals):
                q_by[i][sg] = q

            # We'll also collect K2P values for density plot
            k2p_pre = []
            k2p_post = []

            for i, (header, seq, meta) in enumerate(ltrs):
                chrom = meta["chrom"]
                chrom_sg = chrom_to_sg.get(chrom, None)
                chrom_sg_str = f"SG{chrom_sg+1}" if chrom_sg is not None else "NA"

                hits, total_k = per_ltr_counts[i]
                qs = [q_by[i].get(sg, 1.0) for sg in range(N)]
                enriched = [sg for sg in range(N) if qs[sg] <= args.enrich_p]

                best_sg = None
                best_q = 1.0
                if enriched:
                    best_sg = min(enriched, key=lambda sg: (qs[sg], -hits[sg]))
                    best_q = qs[best_sg]

                enriched_str = ",".join([f"SG{sg+1}" for sg in enriched]) if enriched else ""
                best_str = f"SG{best_sg+1}" if best_sg is not None else ""

                if chrom_sg is None:
                    if len(enriched) == 1:
                        category = "ambiguous"
                    else:
                        category = "post_polyploid"
                else:
                    if len(enriched) == 1 and best_sg == chrom_sg:
                        category = "pre_polyploid"
                    else:
                        category = "post_polyploid"

                k2p = k2p_map.get(header, float("nan"))

                # collect for density plot (only pre vs post, ignore ambiguous)
                if np.isfinite(k2p):
                    if category == "pre_polyploid":
                        k2p_pre.append(k2p)
                    elif category == "post_polyploid":
                        k2p_post.append(k2p)

                out.write(
                    f"{header}\t{chrom}\t{meta['start']}\t{meta['end']}\t{meta['klass']}\t{meta['superfamily']}\t{meta['family']}\t"
                    f"{chrom_sg_str}\t{enriched_str}\t{best_str}\t{best_q:.3g}\t"
                    f"{','.join(str(x) for x in hits)}\t{total_k}\t{category}\t"
                    f"{k2p if np.isfinite(k2p) else 'NA'}\n"
                )

        # Density bins PDF (pre vs post) if we have a K2P file
        if args.k2p_file:
            pdf_out = os.path.join(args.outdir, "k2p_density_pre_vs_post.pdf")
            save_k2p_binned_density_pdf(
                k2p_pre=k2p_pre,
                k2p_post=k2p_post,
                out_pdf=pdf_out,
                bin_size=args.k2p_bin,
                title=f"K2P binned density (bin={args.k2p_bin})",
            )

        # -------------------------
        # Window exchange scan (new)
        # -------------------------
        window_rows = compute_window_enrichment(
            per_ltr_counts=per_ltr_counts,
            ltrs=ltrs,
            chrom_to_sg=chrom_to_sg,
            n_sg=N,
            window_size=args.window_size,
            enrich_q=args.exchange_q,
        )

        # Write window_enrichment.tsv
        win_out = os.path.join(args.outdir, "window_enrichment.tsv")
        with open(win_out, "w", encoding="utf-8") as out:
            out.write(
                "chrom\twindow_start\twindow_end\tchrom_subgenome\tltr_count\t"
                "total_kmers\thits_by_sg\tenriched_subgenomes\tbest_enriched_subgenome\tbest_q\t"
                "significant\tmismatch\n"
            )
            for w in window_rows:
                if w["ltr_count"] < args.min_ltrs_per_window:
                    continue
                chrom_sg = w["chrom_sg"]
                chrom_sg_str = f"SG{chrom_sg+1}" if chrom_sg is not None else "NA"
                enriched_str = ",".join([f"SG{sg+1}" for sg in w["enriched_sgs"]]) if w["enriched_sgs"] else ""
                best_str = f"SG{w['best_sg']+1}" if w["best_sg"] is not None else ""
                hits_str = ",".join(str(int(w["hits_by_sg"][sg])) for sg in range(N))
                out.write(
                    f"{w['chrom']}\t{w['start']}\t{w['end']}\t{chrom_sg_str}\t{w['ltr_count']}\t"
                    f"{w['total_kmers']}\t{hits_str}\t{enriched_str}\t{best_str}\t{w['best_q']:.3g}\t"
                    f"{1 if w['significant'] else 0}\t{1 if w['mismatch'] else 0}\n"
                )

        # Write window_exchange_candidates.tsv
        exch_out = os.path.join(args.outdir, "window_exchange_candidates.tsv")
        n_exch = 0
        with open(exch_out, "w", encoding="utf-8") as out:
            out.write(
                "chrom\twindow_start\twindow_end\tchrom_subgenome\tbest_enriched_subgenome\tbest_q\t"
                "ltr_count\ttotal_kmers\thits_by_sg\tenriched_subgenomes\n"
            )
            for w in window_rows:
                if w["ltr_count"] < args.min_ltrs_per_window:
                    continue
                if not (w["significant"] and w["mismatch"]):
                    continue
                chrom_sg = w["chrom_sg"]
                chrom_sg_str = f"SG{chrom_sg+1}" if chrom_sg is not None else "NA"
                best_str = f"SG{w['best_sg']+1}" if w["best_sg"] is not None else ""
                enriched_str = ",".join([f"SG{sg+1}" for sg in w["enriched_sgs"]]) if w["enriched_sgs"] else ""
                hits_str = ",".join(str(int(w["hits_by_sg"][sg])) for sg in range(N))
                out.write(
                    f"{w['chrom']}\t{w['start']}\t{w['end']}\t{chrom_sg_str}\t{best_str}\t{w['best_q']:.3g}\t"
                    f"{w['ltr_count']}\t{w['total_kmers']}\t{hits_str}\t{enriched_str}\n"
                )
                n_exch += 1

        # -------------------------
        # Optional chromosome painting (new)
        # -------------------------
        if args.genome_fai:
            chrom_lengths = read_fai_lengths(args.genome_fai)
            if not chrom_lengths:
                print(f"[WARN] genome_fai provided but no lengths parsed: {args.genome_fai}", file=sys.stderr)
            else:
                # Prefer plotting in the same chrom order as phasing output, but keep only those present in FAI.
                chrom_order = [c for c in chroms if c in chrom_lengths]
                # include other FAI chroms that have any window rows
                w_chroms = sorted({w["chrom"] for w in window_rows if w["chrom"] in chrom_lengths})
                for c in w_chroms:
                    if c not in chrom_order:
                        chrom_order.append(c)

                paint_name = args.paint_pdf.strip() if args.paint_pdf.strip() else "chromosome_painting_exchanges.pdf"
                paint_pdf = os.path.join(args.outdir, paint_name)
                save_chromosome_painting_pdf(
                    chrom_lengths=chrom_lengths,
                    chrom_order=chrom_order,
                    chrom_sg_map=chrom_to_sg,
                    windows=window_rows,
                    n_sg=N,
                    out_pdf=paint_pdf,
                    window_size=args.window_size,
                    title="Chromosome painting: SG assignment + exchange candidate windows",
                )
        else:
            # If no FAI, we still can infer lengths (useful for sanity), but we won't plot.
            _ = infer_chrom_lengths_from_ltrs(ltrs)

        # Final prints
        print("DONE")
        print(f"Outdir: {args.outdir}")
        print(f"Differential kmers: {diff_kmers_path}")
        print(f"Chromosome phasing: {chrom_tsv}")
        print(f"PCA: {os.path.join(args.outdir, 'pca_chromosomes.pdf')}")
        print(f"Heatmap: {os.path.join(args.outdir, 'heatmap_differential_kmers.pdf')}")
        print(f"LTR enrichment & age TSV (+K2P): {ltr_out}")
        if args.k2p_file:
            print(f"K2P binned density PDF: {os.path.join(args.outdir, 'k2p_density_pre_vs_post.pdf')}")
        print(f"Window enrichment TSV: {win_out}")
        print(f"Exchange candidate windows TSV: {exch_out} (n={n_exch})")
        if args.genome_fai:
            paint_name = args.paint_pdf.strip() if args.paint_pdf.strip() else "chromosome_painting_exchanges.pdf"
            print(f"Chromosome painting PDF: {os.path.join(args.outdir, paint_name)}")

    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    main()
