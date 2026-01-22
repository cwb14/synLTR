#!/usr/bin/env python3
"""
This is a simpler version of the ltr_phasing, but may be more robust to more complicated data (strawberry). 
The previous two versions were used with easy Oryza alta. 
This version gets the A and B genomes of strawberry and mixes C and D. Subphaser gets C and D, but im ok missing them since these are the complicated ones. 
Its also 
# kruskal,wilcoxon seems better than ttest.
# Use low '--q_min'.
# --kmer 17 or --kmer 15.
python ltr_phasing3.py --fasta Fanan.ltr.lib.fa --config Fanan.config --outdir temp5 --threads 200 --canonical -L 3 --f_ratio 2 --q_min 40 --test kruskal --kmer 17

ltr_phasing3.py

Identify subgenome-specific LTR-RT k-mers (default k=15) from a polyploid LTR-RT library,
using homoeologous chromosome sets from a config file.

Pipeline (as you described):
1) Split LTR-RT library FASTA into per-chromosome FASTA using headers like:
   >chrom:start-end#class/superfamily/family
2) For each chromosome FASTA: jellyfish count + dump k-mers (canonical optional; dump -L mincount)
3) Build m x n matrix of k-mer counts (m species, n chromosomes)
4) Normalize to ratios per chromosome => M0
5) Differential k-mer filter:
   For each homoeologous set, sort ratios r0>=r1>=...
   If (r0/r1 >= f) holds for ALL sets and total genome count q >= q_min => keep k-mer
6) Subset M0 to differential k-mers; Z-scale across chromosomes => M1
7) KMeans on chromosomes (samples) using M1 (features), N = number of subgenomes (from config columns)
   Bootstrap by resampling k-mers (rows) with fraction p for B replicates; report stability
8) Visualization: hierarchical clustering heatmap + PCA scatter (PDFs)
9) Subgenome-specific k-mers (post-clustering):
   Group ratios by assigned subgenome; also do per-homoeologous-set paired tests:
   - ttest (paired) OR
   - Kruskal-Wallis across subgenomes OR
   - Wilcoxon signed-rank (paired)
   For each k-mer, “top subgenome” = highest mean ratio across subgenomes.

Requirements:
- jellyfish in PATH
- Python: numpy, pandas, scipy, scikit-learn, matplotlib

Example:
  ./subgenome_ltr_kmers.py \
    --fasta LTRlib.fa \
    --config homeologs.tsv \
    --outdir out \
    --threads 16

Author: (you)
"""

import argparse
import collections
import gzip
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from typing import Dict, List, Tuple, Iterable, Optional

import numpy as np
import pandas as pd

from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
from scipy.stats import ttest_rel, kruskal, wilcoxon

from sklearn.cluster import KMeans
from sklearn.decomposition import PCA


HEADER_RE = re.compile(r"^>([^:]+):(\d+)-(\d+)#(.+)$")  # chrom, start, end, class/superfam/fam


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def check_exe(name: str) -> None:
    if shutil.which(name) is None:
        raise SystemExit(f"ERROR: required executable not found in PATH: {name}")


def open_maybe_gz(path: str, mode: str = "rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def fasta_iter(path: str) -> Iterable[Tuple[str, str]]:
    """Yield (header_without_>, sequence) from FASTA (supports .gz)."""
    with open_maybe_gz(path, "rt") as f:
        header = None
        seq_chunks = []
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


def parse_chrom_from_header(header: str) -> str:
    """
    Expected header form:
      chrom:start-end#class/superfamily/family
    Returns chrom.
    """
    m = HEADER_RE.match(">" + header if not header.startswith(">") else header)
    if not m:
        # try forgiving parse: split at ':' first
        if header.startswith(">"):
            header = header[1:]
        chrom = header.split(":", 1)[0]
        return chrom
    return m.group(1)


def write_per_chrom_fastas(fasta_path: str, outdir: str) -> List[str]:
    """
    Create per-chrom FASTA files containing all LTR-RT sequences from that chromosome.
    Returns list of chromosome names encountered (in sorted order).
    """
    os.makedirs(outdir, exist_ok=True)
    handles = {}
    chroms = set()

    try:
        for hdr, seq in fasta_iter(fasta_path):
            chrom = parse_chrom_from_header(hdr)
            chroms.add(chrom)
            if chrom not in handles:
                handles[chrom] = open(os.path.join(outdir, f"{chrom}.fa"), "wt")
            handles[chrom].write(f">{hdr}\n")
            # wrap 60
            for i in range(0, len(seq), 60):
                handles[chrom].write(seq[i:i+60] + "\n")
    finally:
        for h in handles.values():
            h.close()

    return sorted(chroms)


def run(cmd: List[str]) -> None:
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        raise SystemExit(
            "ERROR running command:\n"
            f"  {' '.join(cmd)}\n\nSTDOUT:\n{p.stdout}\n\nSTDERR:\n{p.stderr}\n"
        )


def jellyfish_count_and_dump(
    fasta_path: str,
    jf_out: str,
    dump_out: str,
    k: int,
    threads: int,
    canonical: bool,
    hash_size: str,
    min_dump_count: int,
):
    """
    jellyfish count -m k -s hash_size -t threads [-C] -o jf_out fasta
    jellyfish dump -c -L min_dump_count jf_out > dump_out
    """
    cmd = ["jellyfish", "count", "-m", str(k), "-s", str(hash_size), "-t", str(threads)]
    if canonical:
        cmd.append("-C")
    cmd += ["-o", jf_out, fasta_path]
    run(cmd)

    cmd = ["jellyfish", "dump", "-c", "-L", str(min_dump_count), jf_out]
    with open(dump_out, "wt") as out:
        p = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)
        if p.returncode != 0:
            raise SystemExit(
                "ERROR running command:\n"
                f"  {' '.join(cmd)}\n\nSTDERR:\n{p.stderr}\n"
            )


def parse_jf_dump_counts(dump_path: str) -> Dict[str, int]:
    """
    jellyfish dump -c format: "<kmer> <count>" per line
    """
    d = {}
    with open(dump_path, "rt") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) != 2:
                continue
            kmer, c = parts
            try:
                d[kmer] = int(c)
            except ValueError:
                continue
    return d


def read_homeolog_config(path: str) -> List[List[str]]:
    """
    Each line: homoeologous chromosome set, tab or space separated.
    Returns list of lists.
    """
    sets = []
    with open(path, "rt") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = re.split(r"\s+", line)
            parts = [p for p in parts if p]
            if len(parts) < 2:
                continue
            sets.append(parts)
    if not sets:
        raise SystemExit(f"ERROR: no homoeolog sets found in config: {path}")
    # sanity: all rows same length?
    widths = {len(x) for x in sets}
    if len(widths) != 1:
        eprint("WARNING: config lines have varying numbers of columns; "
               "subgenome count will be set to the maximum width.")
    return sets


def build_count_matrix(chrom_to_counts: Dict[str, Dict[str, int]]) -> pd.DataFrame:
    """
    Returns DataFrame with rows=k-mers, cols=chromosomes, values=counts (int).
    """
    chroms = sorted(chrom_to_counts.keys())
    # union of kmers
    all_kmers = set()
    for d in chrom_to_counts.values():
        all_kmers.update(d.keys())
    all_kmers = sorted(all_kmers)

    mat = np.zeros((len(all_kmers), len(chroms)), dtype=np.int64)
    kmer_index = {k: i for i, k in enumerate(all_kmers)}
    for j, chrom in enumerate(chroms):
        d = chrom_to_counts[chrom]
        for kmer, c in d.items():
            mat[kmer_index[kmer], j] = c

    df = pd.DataFrame(mat, index=all_kmers, columns=chroms)
    return df


def counts_to_ratios(count_df: pd.DataFrame, eps: float = 0.0) -> pd.DataFrame:
    """
    M0: ratios per chromosome = count / total_kmers_in_chromosome
    """
    totals = count_df.sum(axis=0).astype(float)
    if (totals == 0).any():
        zero_cols = totals[totals == 0].index.tolist()
        eprint(f"WARNING: chromosomes with zero counted kmers (after dump filter) -> {zero_cols}")
    ratios = count_df.astype(float).div(totals + eps, axis=1)
    return ratios


def zscale_rows(df: pd.DataFrame, eps: float = 1e-12) -> pd.DataFrame:
    """
    Z-scale per k-mer across chromosomes: (x - mean)/std
    """
    x = df.values.astype(float)
    mu = x.mean(axis=1, keepdims=True)
    sd = x.std(axis=1, ddof=0, keepdims=True)
    z = (x - mu) / (sd + eps)
    return pd.DataFrame(z, index=df.index, columns=df.columns)


def differential_filter(
    M0: pd.DataFrame,
    homeolog_sets: List[List[str]],
    f_ratio: float,
    q_min: int,
    count_df: pd.DataFrame,
    eps: float = 1e-15
) -> List[str]:
    """
    Keep k-mer if:
      - total genome count q >= q_min
      - for every homeolog set: (r0/r1) >= f_ratio, where r0>=r1 are top two ratios in that set
    Missing chromosomes in a set are treated as 0 ratio.
    """
    chroms_all = list(M0.columns)
    # total genome count per kmer
    q = count_df.sum(axis=1).astype(int)

    keep = []
    for kmer in M0.index:
        if q.loc[kmer] < q_min:
            continue
        ok = True
        for s in homeolog_sets:
            vals = []
            for chrom in s:
                if chrom in M0.columns:
                    vals.append(M0.at[kmer, chrom])
                else:
                    vals.append(0.0)
            # need at least 2 values to compare; if set length is 1, skip
            if len(vals) < 2:
                continue
            vals_sorted = sorted(vals, reverse=True)
            r0, r1 = vals_sorted[0], vals_sorted[1]
            ratio = (r0 + eps) / (r1 + eps)
            if ratio < f_ratio:
                ok = False
                break
        if ok:
            keep.append(kmer)
    return keep


def hungarian_match_labels(ref_labels: np.ndarray, new_labels: np.ndarray, n_clusters: int) -> np.ndarray:
    """
    Match new cluster labels to reference labels via maximum overlap (Hungarian algorithm).
    Returns relabeled new_labels in reference label space.
    """
    from scipy.optimize import linear_sum_assignment

    # contingency matrix: rows=ref, cols=new
    C = np.zeros((n_clusters, n_clusters), dtype=int)
    for r, n in zip(ref_labels, new_labels):
        C[r, n] += 1
    # maximize overlap => minimize negative overlap
    row_ind, col_ind = linear_sum_assignment(-C)
    mapping = {c: r for r, c in zip(row_ind, col_ind)}
    relabeled = np.array([mapping[x] for x in new_labels], dtype=int)
    return relabeled


@dataclass
class KMeansBootstrapResult:
    ref_labels: np.ndarray
    stability: pd.Series  # per chromosome
    label_counts: pd.DataFrame  # chrom x cluster (counts)
    assignments: pd.Series  # chrom -> cluster (int)


def kmeans_with_bootstrap(
    M1: pd.DataFrame,
    n_clusters: int,
    bootstrap_fraction: float,
    bootstrap_reps: int,
    random_state: int,
) -> KMeansBootstrapResult:
    """
    Fit KMeans on chromosomes using M1 (k-mers x chromosomes) => samples=chromosomes.
    Bootstrap by resampling k-mers (rows) without replacement.
    """
    X = M1.T.values  # chrom x kmer
    chroms = M1.columns.tolist()

    km = KMeans(n_clusters=n_clusters, n_init="auto", random_state=random_state)
    ref_labels = km.fit_predict(X)

    rng = np.random.default_rng(random_state)
    n_kmers = M1.shape[0]
    bsize = max(2, int(round(n_kmers * bootstrap_fraction)))

    counts = np.zeros((len(chroms), n_clusters), dtype=int)

    for b in range(bootstrap_reps):
        idx = rng.choice(n_kmers, size=bsize, replace=False)
        Xb = M1.iloc[idx, :].T.values

        km_b = KMeans(n_clusters=n_clusters, n_init="auto", random_state=rng.integers(0, 2**31-1))
        lab_b = km_b.fit_predict(Xb)

        # align labels to reference
        lab_b = hungarian_match_labels(ref_labels, lab_b, n_clusters=n_clusters)

        for i, lab in enumerate(lab_b):
            counts[i, lab] += 1

    label_counts = pd.DataFrame(counts, index=chroms, columns=[f"cluster{c}" for c in range(n_clusters)])
    modal = label_counts.values.argmax(axis=1)
    stability = label_counts.max(axis=1) / bootstrap_reps

    assignments = pd.Series(ref_labels, index=chroms, name="cluster")
    return KMeansBootstrapResult(
        ref_labels=ref_labels,
        stability=stability.rename("bootstrap_stability"),
        label_counts=label_counts,
        assignments=assignments
    )


def plot_heatmap(M1: pd.DataFrame, out_pdf: str, assignments: pd.Series):
    import matplotlib.pyplot as plt

    # order chromosomes by hierarchical clustering
    X = M1.T.values
    dist = pdist(X, metric="euclidean")
    Z = linkage(dist, method="average")
    order = leaves_list(Z)
    chrom_order = [M1.columns[i] for i in order]

    data = M1[chrom_order].values  # kmers x ordered chroms

    plt.figure(figsize=(max(8, len(chrom_order)*0.35), 10))
    plt.imshow(data, aspect="auto", interpolation="nearest")
    plt.colorbar(label="Z-scaled k-mer ratio (M1)")
    plt.yticks([])

    xticks = np.arange(len(chrom_order))
    plt.xticks(xticks, chrom_order, rotation=90, fontsize=7)

    # annotate cluster assignment on top as text line
    top = plt.gca()
    for i, chrom in enumerate(chrom_order):
        top.text(i, -0.5, str(int(assignments.loc[chrom])), ha="center", va="bottom", fontsize=8)

    plt.title("Differential LTR-RT k-mers (M1) heatmap\n(top numbers = KMeans cluster IDs)")
    plt.tight_layout()
    plt.savefig(out_pdf)
    plt.close()


def plot_pca(M1: pd.DataFrame, out_pdf: str, assignments: pd.Series, stability: pd.Series):
    import matplotlib.pyplot as plt

    X = M1.T.values  # chrom x features
    pca = PCA(n_components=2, random_state=0)
    coords = pca.fit_transform(X)
    df = pd.DataFrame(coords, index=M1.columns, columns=["PC1", "PC2"])
    df["cluster"] = assignments.astype(int)
    df["stability"] = stability

    plt.figure(figsize=(8, 6))
    # no explicit colors (per instruction style); matplotlib will cycle defaults
    for cl in sorted(df["cluster"].unique()):
        sub = df[df["cluster"] == cl]
        plt.scatter(sub["PC1"], sub["PC2"], label=f"cluster{cl}", s=60, alpha=0.85)
        for chrom, row in sub.iterrows():
            plt.text(row["PC1"], row["PC2"], chrom, fontsize=7, ha="left", va="bottom")

    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    plt.title("PCA of differential LTR-RT k-mers (M1)")
    plt.legend(loc="best", fontsize=8)
    plt.tight_layout()
    plt.savefig(out_pdf)
    plt.close()


def per_kmer_subgenome_tests(
    M0_diff: pd.DataFrame,
    homeolog_sets: List[List[str]],
    assignments: pd.Series,
    method: str,
    eps: float = 1e-15,
) -> pd.DataFrame:
    """
    For each k-mer:
      - compute mean ratio per subgenome (cluster)
      - top_subgenome = argmax(mean)
      - statistical test using values pulled per homeolog set:
            For each set, take ratio in top-subgenome chromosome within that set (if any),
            and ratio of the best non-top chromosome within the set (max of others).
        This produces paired vectors across sets suitable for paired t-test / wilcoxon.
      - For Kruskal: within each set, we collect the ratio for each subgenome if present;
        then pool across sets to form groups for kruskal.

    Returns DataFrame indexed by k-mer with columns:
      top_cluster, mean_cluster0.., pvalue, n_sets_used
    """
    clusters = sorted(assignments.unique().tolist())
    chrom_to_cluster = assignments.to_dict()

    # Precompute: for each homoeolog set, map cluster -> chromosomes in that set
    set_cluster_chroms = []
    for s in homeolog_sets:
        m = collections.defaultdict(list)
        for chrom in s:
            if chrom in chrom_to_cluster:
                m[chrom_to_cluster[chrom]].append(chrom)
        set_cluster_chroms.append(m)

    rows = []
    for kmer, row in M0_diff.iterrows():
        # mean ratios by cluster across all chromosomes assigned to cluster
        means = {}
        for cl in clusters:
            chroms = assignments.index[assignments == cl].tolist()
            if chroms:
                means[cl] = float(row[chroms].mean())
            else:
                means[cl] = 0.0
        top = max(means.items(), key=lambda x: x[1])[0]

        pval = np.nan
        n_used = 0

        if method.lower() in {"ttest", "t", "paired_t"}:
            a = []
            b = []
            for m in set_cluster_chroms:
                top_chroms = m.get(top, [])
                other_chroms = []
                for cl, chroms in m.items():
                    if cl != top:
                        other_chroms.extend(chroms)
                if not top_chroms or not other_chroms:
                    continue
                # if multiple per cluster appear in set, use max as “most abundant”
                top_val = max([row[c] for c in top_chroms] + [0.0])
                other_val = max([row[c] for c in other_chroms] + [0.0])
                a.append(top_val)
                b.append(other_val)
            n_used = len(a)
            if n_used >= 3:
                stat, pval = ttest_rel(a, b, alternative="greater")
        elif method.lower() in {"wilcoxon", "w"}:
            a = []
            b = []
            for m in set_cluster_chroms:
                top_chroms = m.get(top, [])
                other_chroms = []
                for cl, chroms in m.items():
                    if cl != top:
                        other_chroms.extend(chroms)
                if not top_chroms or not other_chroms:
                    continue
                top_val = max([row[c] for c in top_chroms] + [0.0])
                other_val = max([row[c] for c in other_chroms] + [0.0])
                a.append(top_val)
                b.append(other_val)
            n_used = len(a)
            if n_used >= 3:
                # wilcoxon signed-rank on paired diffs; alt = greater
                try:
                    stat, pval = wilcoxon(np.array(a) - np.array(b), alternative="greater", zero_method="wilcox")
                except ValueError:
                    pval = np.nan
        elif method.lower() in {"kruskal", "kw"}:
            groups = {cl: [] for cl in clusters}
            for m in set_cluster_chroms:
                # per set, take max value for each cluster present (or skip if absent)
                for cl in clusters:
                    chroms = m.get(cl, [])
                    if not chroms:
                        continue
                    groups[cl].append(max([row[c] for c in chroms] + [0.0]))
            # require at least 2 groups with >=3 points
            usable = [v for v in groups.values() if len(v) >= 3]
            n_used = min([len(v) for v in groups.values() if len(v) > 0], default=0)
            if len(usable) >= 2:
                stat, pval = kruskal(*usable)
        else:
            raise SystemExit(f"ERROR: unknown test method: {method}")

        out = {
            "top_cluster": int(top),
            "pvalue": pval,
            "n_sets_used": int(n_used),
        }
        for cl in clusters:
            out[f"mean_cluster{cl}"] = means[cl]
        rows.append((kmer, out))

    df = pd.DataFrame.from_dict(dict(rows), orient="index")
    df.index.name = "kmer"
    # also compute a simple effect size: top_mean / second_mean
    mean_cols = [c for c in df.columns if c.startswith("mean_cluster")]
    if mean_cols:
        mean_vals = df[mean_cols].values
        sorted_means = np.sort(mean_vals, axis=1)[:, ::-1]
        df["top_over_second_mean"] = (sorted_means[:, 0] + eps) / (sorted_means[:, 1] + eps)
    return df


def main():
    ap = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Identify subgenome-specific LTR-RT k-mers from polyploid LTR library + homoeolog config."
    )
    ap.add_argument("--fasta", required=True, help="LTR-RT library FASTA (.fa/.fasta, optionally .gz)")
    ap.add_argument("--config", required=True, help="Homoeologous chromosome sets config (tab/space separated)")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("-k", "--kmer", type=int, default=15, help="k-mer length")
    ap.add_argument("--canonical", action="store_true", help="Use canonical k-mers (-C) in jellyfish")
    ap.add_argument("-L", "--min_dump_count", type=int, default=3, help="Minimum k-mer count for jellyfish dump (-L)")
    ap.add_argument("--hash_size", default="200M", help="Jellyfish hash size (-s). Increase if many kmers.")
    ap.add_argument("--threads", type=int, default=4, help="Threads for jellyfish")
    ap.add_argument("--f_ratio", type=float, default=2.0, help="Differential threshold f = r0/r1 per homoeolog set")
    ap.add_argument("--q_min", type=int, default=200, help="Minimum total genome count q for a k-mer")
    ap.add_argument("--bootstrap_fraction", type=float, default=0.5, help="Fraction of k-mers to resample per bootstrap")
    ap.add_argument("--bootstrap_reps", type=int, default=1000, help="Bootstrap replicates")
    ap.add_argument("--seed", type=int, default=1, help="Random seed")
    ap.add_argument("--test", choices=["ttest", "kruskal", "wilcoxon"], default="ttest",
                    help="Statistical test for subgenome-specific k-mers (paired-by-homeolog-set where applicable)")
    ap.add_argument("--min_p", type=float, default=0.01, help="Report subgenome-specific k-mers with p <= min_p")
    ap.add_argument("--keep_tmp", action="store_true", help="Keep temporary working directory (for debugging)")
    args = ap.parse_args()

    check_exe("jellyfish")

    os.makedirs(args.outdir, exist_ok=True)
    homeolog_sets = read_homeolog_config(args.config)
    n_subgenomes = max(len(s) for s in homeolog_sets)

    # temp workspace
    tmpdir = tempfile.mkdtemp(prefix="subgenome_ltr_kmers_")
    try:
        eprint(f"[INFO] tempdir: {tmpdir}")

        per_chr_dir = os.path.join(tmpdir, "per_chrom_fastas")
        eprint("[INFO] splitting LTR library into per-chromosome FASTAs...")
        chroms_seen = write_per_chrom_fastas(args.fasta, per_chr_dir)
        eprint(f"[INFO] chromosomes with LTR-RT entries: {len(chroms_seen)}")

        # run jellyfish per chromosome
        chrom_to_counts: Dict[str, Dict[str, int]] = {}
        jf_dir = os.path.join(tmpdir, "jellyfish")
        os.makedirs(jf_dir, exist_ok=True)

        eprint("[INFO] counting and dumping k-mers with jellyfish...")
        for chrom in chroms_seen:
            fa = os.path.join(per_chr_dir, f"{chrom}.fa")
            jf = os.path.join(jf_dir, f"{chrom}.jf")
            dump = os.path.join(jf_dir, f"{chrom}.dump.txt")

            jellyfish_count_and_dump(
                fasta_path=fa,
                jf_out=jf,
                dump_out=dump,
                k=args.kmer,
                threads=args.threads,
                canonical=args.canonical,
                hash_size=args.hash_size,
                min_dump_count=args.min_dump_count,
            )
            chrom_to_counts[chrom] = parse_jf_dump_counts(dump)

        eprint("[INFO] building k-mer count matrix...")
        count_df = build_count_matrix(chrom_to_counts)
        count_df.to_csv(os.path.join(args.outdir, "kmer_counts.tsv"), sep="\t")

        eprint("[INFO] normalizing to ratios (M0)...")
        M0 = counts_to_ratios(count_df, eps=0.0)
        M0.to_csv(os.path.join(args.outdir, "M0_ratios.tsv"), sep="\t")

        eprint("[INFO] filtering differential k-mers...")
        diff_kmers = differential_filter(
            M0=M0,
            homeolog_sets=homeolog_sets,
            f_ratio=args.f_ratio,
            q_min=args.q_min,
            count_df=count_df,
        )
        eprint(f"[INFO] differential k-mers retained: {len(diff_kmers)}")
        with open(os.path.join(args.outdir, "differential_kmers.txt"), "wt") as out:
            for kmer in diff_kmers:
                out.write(kmer + "\n")

        if len(diff_kmers) < 10:
            raise SystemExit("ERROR: too few differential k-mers after filtering. "
                             "Consider lowering --q_min or --f_ratio, or lowering -L.")

        M0_diff = M0.loc[diff_kmers, :]
        M0_diff.to_csv(os.path.join(args.outdir, "M0_diff.tsv"), sep="\t")

        eprint("[INFO] Z-scaling across chromosomes (M1)...")
        M1 = zscale_rows(M0_diff)
        M1.to_csv(os.path.join(args.outdir, "M1_zscaled.tsv"), sep="\t")

        eprint(f"[INFO] KMeans clustering chromosomes into N={n_subgenomes} subgenomes...")
        boot = kmeans_with_bootstrap(
            M1=M1,
            n_clusters=n_subgenomes,
            bootstrap_fraction=args.bootstrap_fraction,
            bootstrap_reps=args.bootstrap_reps,
            random_state=args.seed,
        )

        # outputs: assignments + stability
        assign_df = pd.DataFrame({
            "chromosome": boot.assignments.index,
            "cluster": boot.assignments.values.astype(int),
            "bootstrap_stability": boot.stability.loc[boot.assignments.index].values
        })
        assign_df.to_csv(os.path.join(args.outdir, "chromosome_subgenome_assignments.tsv"),
                         sep="\t", index=False)

        boot.label_counts.to_csv(os.path.join(args.outdir, "bootstrap_cluster_counts.tsv"), sep="\t")

        eprint("[INFO] plotting heatmap + PCA...")
        plot_heatmap(M1, os.path.join(args.outdir, "M1_heatmap_hclust.pdf"), boot.assignments)
        plot_pca(M1, os.path.join(args.outdir, "M1_PCA.pdf"), boot.assignments, boot.stability)

        eprint(f"[INFO] testing subgenome-specific k-mers using: {args.test}")
        test_df = per_kmer_subgenome_tests(
            M0_diff=M0_diff,
            homeolog_sets=homeolog_sets,
            assignments=boot.assignments,
            method=args.test,
        )
        test_df.to_csv(os.path.join(args.outdir, "kmer_subgenome_tests.tsv"), sep="\t")

        sig = test_df[(~test_df["pvalue"].isna()) & (test_df["pvalue"] <= args.min_p)].copy()
        sig.sort_values(["pvalue", "top_over_second_mean"], ascending=[True, False], inplace=True)
        sig.to_csv(os.path.join(args.outdir, "subgenome_specific_kmers.tsv"), sep="\t")
        eprint(f"[INFO] significant subgenome-specific k-mers (p <= {args.min_p}): {sig.shape[0]}")

        eprint("[DONE] Outputs written to:", args.outdir)

    finally:
        if args.keep_tmp:
            eprint(f"[INFO] keeping tempdir: {tmpdir}")
        else:
            shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == "__main__":
    main()
