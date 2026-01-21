#!/usr/bin/env python3
"""
Note 'ltrharvest_phasing.py' and 'ltrharvest_phasing2.py' use slightly different approaches and give slightly different results. Im not sure which is better, so theyre both included. 

# Get non-nest.
python ltrharvest.py --genome Oalta.fa --proteins Osati.pep --threads 20 --out-prefix Oalta_ltr --scn-min-ltr-len 100 --scn-min-ret-len 800 --scn-max-ret-len 15000 --scn-min-int-len 500 --scn-max-int-len 12000

# Get 1-level nesters. 
python v2/mask_ltr.py --features-fasta Oalta.ltrharvest.full_length.dedup.fa.rexdb-plant.cls.lib.fa --genome Oalta.fa --feature-character N --far-character V --distance 15000 > Oalta_r1.fa
python ./ltrharvest.py --require-run-chars N --genome Oalta_r1.fa --proteins Osati.pep --threads 100 --out-prefix Oalta_r2 --scn-min-ltr-len 100 --scn-min-ret-len 1000 --scn-max-ret-len 30000 --scn-min-int-len 200 --scn-max-int-len 28000 --ltrharvest-args '-mindistltr 100 -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -similar 70 -vic 30 -seed 15 -seqids yes -xdrop 10 -maxdistltr 30000' --ltrfinder-args '-w 2 -C -D 30000 -d 100 -L 7000 -l 100 -p 20 -M 0.00 -S 0.0'

# Merge non-nest and 1-level nesters. 
sed '/^>/! { s/[^ATCGatcg]//g; /^$/d }' Oalta_r2.ltrharvest.full_length.dedup.fa.rexdb-plant.cls.lib.fa >temp.fa
cat temp.fa Oalta.ltrharvest.full_length.dedup.fa.rexdb-plant.cls.lib.fa > Oalta.ltr.lib.fa

cat sg.config 
CC1	DD1
CC2	DD2
CC3	DD3
CC4	DD4
CC5	DD5
CC6	DD6
CC7	DD7
CC8	DD8
CC9	DD9
CC10	DD10
CC11	DD11
CC12	DD12

# Now Phase. 
python ltrharvest_phasing.py --ltr_fasta Oalta.ltr.lib.fa --homoeolog_config sg.config --outdir ltr_phased --n_subgenomes 2 --canonical -k 15 -L 3 --threads 200 --bootstrap 1000 --k2p_tsv <(cat *_dedup) --genome_fai Oalta.fa.fai

LTR-only SubPhaser-style subgenome phasing + LTR subphasing + K2P annotation +
K2P window-density plot + 1 Mb genomic-window enrichment scan +
(OPTIONAL) chromosome plot using genome .fai.

New (chromosome plot):
  - If --genome_fai is provided, emits chromosome_plot.pdf
  - Each chromosome is a horizontal bar colored by chromosome-level subgenome assignment
  - Mismatch windows (potential exchanges/switch errors) are overlaid as segments colored by
    the window's best enriched subgenome, outlined in black.

Window enrichment is still LTR-driven (windows are evaluated using LTRs that fall into them).
"""

from __future__ import annotations

import argparse
import gzip
import os
import re
import shutil
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional

import numpy as np

# Optional heavy deps: sklearn/scipy/matplotlib. Fail with clear message if absent.
try:
    from sklearn.cluster import KMeans
    from sklearn.decomposition import PCA
except ImportError as e:
    raise SystemExit("ERROR: scikit-learn is required (pip install scikit-learn).") from e

try:
    from scipy.stats import ttest_ind, fisher_exact
    from scipy.cluster.hierarchy import linkage, leaves_list
except ImportError as e:
    raise SystemExit("ERROR: scipy is required (pip install scipy).") from e

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
except ImportError as e:
    raise SystemExit("ERROR: matplotlib is required (pip install matplotlib).") from e


HEADER_RE = re.compile(
    r"^>(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)#(?P<class>[^/\s#]+)/(?P<superfamily>[^/\s#]+)/(?P<family>[^\s#]+)\s*$"
)


def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(comp)[::-1]


def canonical_kmer(kmer: str) -> str:
    rc = revcomp(kmer)
    return rc if rc < kmer else kmer


def bh_adjust(pvals: List[float]) -> List[float]:
    """Benjamini-Hochberg FDR adjustment (returns q-values)."""
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = np.empty(n, dtype=float)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        idx = order[i]
        rank = i + 1
        q = min(prev, pvals[idx] * n / rank)
        ranked[idx] = q
        prev = q
    return ranked.tolist()


@dataclass
class LTRRecord:
    chrom: str
    start: int
    end: int
    klass: str
    superfamily: str
    family: str
    seq: str

    def key(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end}#{self.klass}/{self.superfamily}/{self.family}"


def read_ltr_fasta(path: str) -> Tuple[List[LTRRecord], Dict[str, List[int]]]:
    records: List[LTRRecord] = []
    chrom_to_indices: Dict[str, List[int]] = defaultdict(list)

    with open(path, "r") as f:
        header = None
        seq_chunks: List[str] = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    rec = parse_record(header, "".join(seq_chunks))
                    chrom_to_indices[rec.chrom].append(len(records))
                    records.append(rec)
                header = line
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header is not None:
            rec = parse_record(header, "".join(seq_chunks))
            chrom_to_indices[rec.chrom].append(len(records))
            records.append(rec)

    return records, chrom_to_indices


def parse_record(header: str, seq: str) -> LTRRecord:
    m = HEADER_RE.match(header)
    if not m:
        raise ValueError(
            f"Bad header format:\n  {header}\nExpected:\n  >Chr:start-end#class/superfamily/family"
        )
    chrom = m.group("chrom")
    start = int(m.group("start"))
    end = int(m.group("end"))
    klass = m.group("class")
    sf = m.group("superfamily")
    fam = m.group("family")
    seq = seq.upper().replace("U", "T")
    return LTRRecord(chrom, start, end, klass, sf, fam, seq)


def read_homoeolog_config(path: str) -> List[List[str]]:
    sets: List[List[str]] = []
    with open(path, "r") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if "#" in line:
                line = line.split("#", 1)[0].strip()
            parts = re.split(r"\s+", line)
            parts = [p for p in parts if p]
            if parts:
                sets.append(parts)
    return sets


def read_k2p_map(path: str, k2p_col_1based: int = 11) -> Dict[str, float]:
    if k2p_col_1based < 2:
        raise ValueError("k2p_col_1based must be >=2 (since col1 is the key).")
    idx = k2p_col_1based - 1  # 0-based

    m: Dict[str, float] = {}
    with open(path, "r") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = re.split(r"\s+", line)
            if len(parts) <= idx:
                continue
            key = parts[0].lstrip(">")
            try:
                val = float(parts[idx])
            except ValueError:
                continue
            m[key] = val
    return m


def read_fai_lengths(path: str) -> Dict[str, int]:
    """
    .fai format: name length offset linebases linewidth
    """
    lengths: Dict[str, int] = {}
    with open(path, "r") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                parts = line.split()
            if len(parts) < 2:
                continue
            name = parts[0]
            try:
                L = int(parts[1])
            except ValueError:
                continue
            lengths[name] = L
    return lengths


def ensure_dir(p: str) -> None:
    os.makedirs(p, exist_ok=True)


def which(exe: str) -> Optional[str]:
    return shutil.which(exe)


def write_chrom_concat_fastas(records: List[LTRRecord],
                              chrom_to_indices: Dict[str, List[int]],
                              outdir: str,
                              sep_n: int = 50) -> Dict[str, str]:
    chrom_fas: Dict[str, str] = {}
    for chrom, idxs in chrom_to_indices.items():
        fp = os.path.join(outdir, f"{chrom}.ltr_concat.fa")
        with open(fp, "w") as out:
            out.write(f">{chrom}\n")
            pieces = [records[i].seq for i in idxs]
            seq = ("N" * sep_n).join(pieces)
            for j in range(0, len(seq), 80):
                out.write(seq[j:j+80] + "\n")
        chrom_fas[chrom] = fp
    return chrom_fas


def compute_total_kmers_per_chrom(records: List[LTRRecord],
                                  chrom_to_indices: Dict[str, List[int]],
                                  k: int) -> Dict[str, int]:
    totals: Dict[str, int] = {}
    for chrom, idxs in chrom_to_indices.items():
        total = 0
        for i in idxs:
            L = len(records[i].seq)
            if L >= k:
                total += (L - k + 1)
        totals[chrom] = total
    return totals


def jellyfish_count_and_dump(chrom_fa: str, k: int, size: int, canonical: bool, Lmin: int,
                             outprefix: str, threads: int) -> str:
    jf = f"{outprefix}.jf"
    dump = f"{outprefix}.dump"
    cmd_count = ["jellyfish", "count", "-m", str(k), "-s", str(size), "-o", jf, "-t", str(threads)]
    if canonical:
        cmd_count.append("--canonical")
    cmd_count.append(chrom_fa)

    subprocess.run(cmd_count, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    cmd_dump = ["jellyfish", "dump", "-c", "-L", str(Lmin), jf]
    with open(dump, "w") as out:
        subprocess.run(cmd_dump, check=True, stdout=out, stderr=subprocess.DEVNULL)
    return dump


def python_kmer_count_dump(chrom_fa: str, k: int, canonical: bool, Lmin: int, dump_path: str) -> None:
    seq: List[str] = []
    with open(chrom_fa, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq.append(line.upper())
    s = "".join(seq)
    counts = defaultdict(int)
    for i in range(0, len(s) - k + 1):
        kmer = s[i:i+k]
        if "N" in kmer:
            continue
        if canonical:
            kmer = canonical_kmer(kmer)
        counts[kmer] += 1
    with open(dump_path, "w") as out:
        for kmer, c in counts.items():
            if c >= Lmin:
                out.write(f"{kmer} {c}\n")


def build_kmer_count_matrix(chroms: List[str],
                            chrom_fas: Dict[str, str],
                            total_kmers: Dict[str, int],
                            k: int,
                            jf_size: int,
                            canonical: bool,
                            Lmin: int,
                            threads: int,
                            workdir: str) -> Tuple[List[str], np.ndarray]:
    have_jf = which("jellyfish") is not None
    kmer_to_row: Dict[str, int] = {}
    rows: List[np.ndarray] = []
    n = len(chroms)

    def get_row(kmer: str) -> int:
        if kmer in kmer_to_row:
            return kmer_to_row[kmer]
        idx = len(rows)
        kmer_to_row[kmer] = idx
        rows.append(np.zeros(n, dtype=np.int32))
        return idx

    for j, chrom in enumerate(chroms):
        if total_kmers.get(chrom, 0) <= 0:
            continue

        outprefix = os.path.join(workdir, chrom)
        dump_path = f"{outprefix}.dump"

        if have_jf:
            dump_path = jellyfish_count_and_dump(
                chrom_fas[chrom],
                k=k,
                size=jf_size,
                canonical=canonical,
                Lmin=Lmin,
                outprefix=outprefix,
                threads=threads,
            )
        else:
            python_kmer_count_dump(chrom_fas[chrom], k=k, canonical=canonical, Lmin=Lmin, dump_path=dump_path)

        with open(dump_path, "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split()
                if len(parts) != 2:
                    continue
                kmer, c = parts[0], int(parts[1])
                r = get_row(kmer)
                rows[r][j] = c

    kmers = [None] * len(rows)
    for kmer, idx in kmer_to_row.items():
        kmers[idx] = kmer

    counts = np.vstack(rows) if rows else np.zeros((0, n), dtype=np.int32)
    return kmers, counts


def identify_differential_kmers(kmers: List[str],
                               counts: np.ndarray,
                               chroms: List[str],
                               total_kmers: Dict[str, int],
                               homoeolog_sets: List[List[str]],
                               fmin: float,
                               qmin: int) -> np.ndarray:
    chrom_index = {c: i for i, c in enumerate(chroms)}
    denoms = np.array([max(1, total_kmers.get(c, 1)) for c in chroms], dtype=np.float64)
    ratios = counts.astype(np.float64) / denoms[None, :]

    q = counts.sum(axis=1).astype(np.int64)
    keep = (q >= qmin)

    for hs in homoeolog_sets:
        hs = [c for c in hs if c in chrom_index]
        if len(hs) < 2:
            continue
        cols = np.array([chrom_index[c] for c in hs], dtype=np.int32)

        sub = ratios[:, cols]
        r_sorted = np.sort(sub, axis=1)[:, ::-1]
        r0 = r_sorted[:, 0]
        r1 = r_sorted[:, 1]
        fold = np.where((r1 > 0), (r0 / r1), np.where(r0 > 0, np.inf, 0.0))
        keep &= (fold >= fmin)

    var = ratios.var(axis=1)
    keep &= (var > 0)

    return keep


def zscore_rows(mat: np.ndarray) -> np.ndarray:
    mu = mat.mean(axis=1, keepdims=True)
    sd = mat.std(axis=1, keepdims=True)
    sd = np.where(sd == 0, 1.0, sd)
    return (mat - mu) / sd


def bootstrap_kmeans(chrom_features: np.ndarray,
                     n_clusters: int,
                     n_boot: int,
                     resample_frac: float,
                     seed: int) -> Tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)

    km = KMeans(n_clusters=n_clusters, n_init="auto", random_state=seed)
    consensus = km.fit_predict(chrom_features)

    matches = np.zeros((n_boot, chrom_features.shape[0]), dtype=np.int8)

    n_feat = chrom_features.shape[1]
    m = max(2, int(round(resample_frac * n_feat)))

    for b in range(n_boot):
        feat_idx = rng.choice(n_feat, size=m, replace=False)
        Xb = chrom_features[:, feat_idx]
        km_b = KMeans(n_clusters=n_clusters, n_init="auto", random_state=seed + b + 1)
        lab_b = km_b.fit_predict(Xb)

        cont = np.zeros((n_clusters, n_clusters), dtype=np.int32)
        for i in range(len(consensus)):
            cont[lab_b[i], consensus[i]] += 1

        mapping = [-1] * n_clusters
        used = set()
        for src in range(n_clusters):
            best = np.argsort(cont[src])[::-1]
            for tgt in best:
                if tgt not in used:
                    mapping[src] = tgt
                    used.add(tgt)
                    break
            if mapping[src] == -1:
                mapping[src] = int(best[0])

        lab_aligned = np.array([mapping[x] for x in lab_b], dtype=np.int32)
        matches[b, :] = (lab_aligned == consensus).astype(np.int8)

    stability = matches.mean(axis=0)
    return consensus, stability


def call_subgenome_specific_kmers(diff_ratios: np.ndarray,
                                 chrom_labels: np.ndarray,
                                 pval: float) -> Dict[int, np.ndarray]:
    N = int(chrom_labels.max()) + 1
    specific = {g: np.zeros(diff_ratios.shape[0], dtype=bool) for g in range(N)}

    idx_by_g = {g: np.where(chrom_labels == g)[0] for g in range(N)}
    all_idx = np.arange(diff_ratios.shape[1])

    for g in range(N):
        idx_g = idx_by_g[g]
        idx_o = np.array([i for i in all_idx if i not in set(idx_g)], dtype=np.int32)
        if len(idx_g) < 2 or len(idx_o) < 2:
            continue

        Xg = diff_ratios[:, idx_g]
        Xo = diff_ratios[:, idx_o]

        mean_g = Xg.mean(axis=1)
        mean_o = Xo.mean(axis=1)

        try:
            pv = ttest_ind(Xg, Xo, axis=1, equal_var=False, alternative="greater").pvalue
        except TypeError:
            pv2 = ttest_ind(Xg, Xo, axis=1, equal_var=False).pvalue
            pv = np.where(mean_g > mean_o, pv2 / 2.0, 1.0)

        specific[g] = (pv < pval)

    return specific


def plot_pca(chrom_features: np.ndarray, chroms: List[str], labels: np.ndarray, out_pdf: str) -> None:
    pca = PCA(n_components=2, random_state=1)
    coords = pca.fit_transform(chrom_features)
    plt.figure(figsize=(7, 6))
    for g in sorted(set(labels.tolist())):
        idx = np.where(labels == g)[0]
        plt.scatter(coords[idx, 0], coords[idx, 1], label=f"subg{g}", s=45)
        for i in idx:
            plt.text(coords[i, 0], coords[i, 1], chroms[i], fontsize=8)
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(out_pdf)
    plt.close()


def plot_heatmap(M1: np.ndarray, chroms: List[str], out_pdf: str, max_kmers: int, seed: int) -> None:
    rng = np.random.default_rng(seed)
    R = M1
    if R.shape[0] > max_kmers:
        idx = rng.choice(R.shape[0], size=max_kmers, replace=False)
        R = R[idx, :]

    C = np.corrcoef(R.T)
    dist = 1.0 - C
    iu = np.triu_indices(dist.shape[0], k=1)
    condensed = dist[iu]
    Z = linkage(condensed, method="average")
    order = leaves_list(Z)
    R2 = R[:, order]
    chroms2 = [chroms[i] for i in order]

    plt.figure(figsize=(max(6, 0.35 * len(chroms2)), 8))
    plt.imshow(R2, aspect="auto", interpolation="nearest")
    plt.xticks(range(len(chroms2)), chroms2, rotation=90, fontsize=8)
    plt.yticks([])
    plt.colorbar(label="Z-scaled k-mer ratio")
    plt.tight_layout()
    plt.savefig(out_pdf)
    plt.close()


def count_ltr_kmer_hits(seq: str, k: int, canonical: bool, subg_sets: Dict[int, set]) -> Dict[int, int]:
    hits = {g: 0 for g in subg_sets.keys()}
    if len(seq) < k:
        return hits
    for i in range(0, len(seq) - k + 1):
        kmer = seq[i:i+k]
        if "N" in kmer:
            continue
        if canonical:
            kmer = canonical_kmer(kmer)
        for g, s in subg_sets.items():
            if kmer in s:
                hits[g] += 1
    return hits


def plot_k2p_window_density(k2p_by_group: Dict[str, List[float]],
                            out_pdf: str,
                            window: float = 0.0005) -> None:
    groups = [
        "POST_polyploidy_shared_or_mixed",
        "PRE_polyploidy_subgenome_specific",
        "EXCHANGE_OR_ERROR",
    ]

    all_vals: List[float] = []
    for g in groups:
        all_vals.extend([v for v in k2p_by_group.get(g, []) if v is not None and np.isfinite(v)])

    if len(all_vals) == 0:
        plt.figure(figsize=(8, 5))
        plt.text(0.5, 0.5, "No K2P values available for plotting.", ha="center", va="center")
        plt.axis("off")
        plt.tight_layout()
        plt.savefig(out_pdf)
        plt.close()
        return

    vmin = min(all_vals)
    vmax = max(all_vals)

    start = np.floor(vmin / window) * window
    end = np.ceil(vmax / window) * window
    if end <= start:
        end = start + window

    bins = np.arange(start, end + window, window)

    plt.figure(figsize=(9, 5.5))
    plotted_any = False
    for g in groups:
        vals = [v for v in k2p_by_group.get(g, []) if v is not None and np.isfinite(v)]
        if len(vals) == 0:
            continue
        plt.hist(vals, bins=bins, density=True, alpha=0.35, label=f"{g} (n={len(vals)})")
        plotted_any = True

    if not plotted_any:
        plt.text(0.5, 0.5, "No K2P values available for plotting.", ha="center", va="center")
        plt.axis("off")
    else:
        plt.xlabel("K2P divergence")
        plt.ylabel("Density (area=1 per group)")
        plt.legend(frameon=False)
        plt.tight_layout()

    plt.savefig(out_pdf)
    plt.close()


def window_enrichment_scan(records: List[LTRRecord],
                           subg_sets: Dict[int, set],
                           chrom_to_subg: Dict[str, int],
                           k: int,
                           canonical: bool,
                           window_size: int,
                           window_alpha: float,
                           min_ltrs_per_window: int,
                           min_total_hits_per_window: int,
                           k2p_map: Dict[str, float]) -> Tuple[List[dict], List[dict], int]:
    subg_ids = sorted(subg_sets.keys())

    acc = defaultdict(lambda: defaultdict(lambda: {
        "sum_hits": {g: 0 for g in subg_ids},
        "n_ltrs": 0,
        "k2p_vals": [],
        "min_pos": None,
        "max_pos": None,
    }))

    missing_k2p_used = 0

    for rec in records:
        chrom_subg = chrom_to_subg.get(rec.chrom, -1)
        if chrom_subg == -1:
            continue

        win_id = rec.start // window_size
        a = acc[rec.chrom][win_id]

        hits = count_ltr_kmer_hits(rec.seq, k=k, canonical=canonical, subg_sets=subg_sets)
        for g in subg_ids:
            a["sum_hits"][g] += int(hits.get(g, 0))

        a["n_ltrs"] += 1
        a["min_pos"] = rec.start if a["min_pos"] is None else min(a["min_pos"], rec.start)
        a["max_pos"] = rec.end if a["max_pos"] is None else max(a["max_pos"], rec.end)

        kv = k2p_map.get(rec.key(), float("nan"))
        if np.isfinite(kv):
            a["k2p_vals"].append(float(kv))
        else:
            missing_k2p_used += 1

    all_windows: List[dict] = []
    mismatches: List[dict] = []

    bg_sizes = {g: len(subg_sets[g]) for g in subg_ids}

    for chrom, bywin in acc.items():
        chrom_subg = chrom_to_subg.get(chrom, -1)
        if chrom_subg == -1:
            continue

        for win_id, a in sorted(bywin.items(), key=lambda x: x[0]):
            n_ltrs = a["n_ltrs"]
            sum_hits = a["sum_hits"]
            total_hits = sum(sum_hits.values())

            win_start = int(win_id * window_size)
            win_end = int(win_start + window_size)

            if n_ltrs < min_ltrs_per_window or total_hits < min_total_hits_per_window:
                row = {
                    "chrom": chrom,
                    "win_start": win_start,
                    "win_end": win_end,
                    "chrom_subgenome": f"subg{chrom_subg}",
                    "best_enriched_subgenome": "NA",
                    "best_p": np.nan,
                    "best_q": np.nan,
                    "call": "UNINFORMATIVE",
                    "n_ltrs": n_ltrs,
                    "total_hits": total_hits,
                    "mean_k2p": (np.mean(a["k2p_vals"]) if a["k2p_vals"] else np.nan),
                    "min_ltr_start": a["min_pos"] if a["min_pos"] is not None else "NA",
                    "max_ltr_end": a["max_pos"] if a["max_pos"] is not None else "NA",
                    "sum_hits": {g: sum_hits[g] for g in subg_ids},
                }
                all_windows.append(row)
                continue

            pvals: List[float] = []
            for g in subg_ids:
                c00 = sum_hits[g]
                c01 = total_hits - c00
                bg_g = bg_sizes[g]
                bg_o = sum(bg_sizes[x] for x in subg_ids if x != g)
                _, p = fisher_exact([[c00, c01], [bg_g, bg_o]], alternative="greater")
                pvals.append(p)

            qvals = bh_adjust(pvals)
            best_i = int(np.argmin(qvals))
            best_g = subg_ids[best_i]
            best_p = pvals[best_i]
            best_q = qvals[best_i]

            if best_q < window_alpha and best_g != chrom_subg:
                call = "POTENTIAL_EXCHANGE_OR_SWITCH_ERROR"
            elif best_q < window_alpha and best_g == chrom_subg:
                call = "MATCHES_CHROMOSOME_ASSIGNMENT"
            else:
                call = "UNINFORMATIVE"

            row = {
                "chrom": chrom,
                "win_start": win_start,
                "win_end": win_end,
                "chrom_subgenome": f"subg{chrom_subg}",
                "best_enriched_subgenome": f"subg{best_g}",
                "best_p": best_p,
                "best_q": best_q,
                "call": call,
                "n_ltrs": n_ltrs,
                "total_hits": total_hits,
                "mean_k2p": (np.mean(a["k2p_vals"]) if a["k2p_vals"] else np.nan),
                "min_ltr_start": a["min_pos"] if a["min_pos"] is not None else "NA",
                "max_ltr_end": a["max_pos"] if a["max_pos"] is not None else "NA",
                "sum_hits": {g: sum_hits[g] for g in subg_ids},
            }
            all_windows.append(row)
            if call == "POTENTIAL_EXCHANGE_OR_SWITCH_ERROR":
                mismatches.append(row)

    return all_windows, mismatches, missing_k2p_used


def write_window_tables(all_windows: List[dict],
                        mismatches: List[dict],
                        out_all: str,
                        out_mis: str,
                        subg_ids: List[int]) -> None:
    def write_one(rows: List[dict], path: str) -> None:
        with open(path, "w") as out:
            header = [
                "chrom", "win_start", "win_end",
                "chrom_subgenome",
                "best_enriched_subgenome", "best_p", "best_q",
                "call",
                "n_ltrs", "total_hits",
                "mean_k2p",
                "min_ltr_start", "max_ltr_end",
            ] + [f"sum_hits_subg{g}" for g in subg_ids]
            out.write("\t".join(header) + "\n")

            for r in rows:
                best_p = r["best_p"]
                best_q = r["best_q"]
                mp = r["mean_k2p"]
                row = [
                    str(r["chrom"]),
                    str(r["win_start"]),
                    str(r["win_end"]),
                    str(r["chrom_subgenome"]),
                    str(r["best_enriched_subgenome"]),
                    f"{best_p:.3e}" if isinstance(best_p, (float, np.floating)) and np.isfinite(best_p) else "NA",
                    f"{best_q:.3e}" if isinstance(best_q, (float, np.floating)) and np.isfinite(best_q) else "NA",
                    str(r["call"]),
                    str(r["n_ltrs"]),
                    str(r["total_hits"]),
                    f"{mp:.6f}" if isinstance(mp, (float, np.floating)) and np.isfinite(mp) else "NA",
                    str(r["min_ltr_start"]),
                    str(r["max_ltr_end"]),
                ]
                for g in subg_ids:
                    row.append(str(int(r["sum_hits"][g])))
                out.write("\t".join(row) + "\n")

    write_one(all_windows, out_all)
    write_one(mismatches, out_mis)


def get_subg_color_map(n_subgenomes: int):
    """
    Deterministic color assignment using matplotlib's tab20 cycle.
    Returns dict: subg index -> RGBA
    """
    cmap = plt.get_cmap("tab20")
    colors = {}
    for g in range(n_subgenomes):
        colors[g] = cmap(g % 20)
    return colors


def plot_chromosomes_with_exchanges(genome_lengths: Dict[str, int],
                                   chroms: List[str],
                                   chrom_to_subg: Dict[str, int],
                                   mismatches: List[dict],
                                   n_subgenomes: int,
                                   out_pdf: str,
                                   max_chroms_per_page: int = 60) -> None:
    """
    Plots chromosomes as horizontal bars colored by chrom-level subgenome assignment,
    overlays mismatch windows as colored segments using the window's best enriched subgenome.

    - chroms controls the order (we plot only those present in chroms unless lengths contain extras)
    - if there are many chromosomes, auto-scales figure height
    """
    subg_colors = get_subg_color_map(n_subgenomes)

    # group mismatch windows by chrom
    mis_by_chrom = defaultdict(list)
    for r in mismatches:
        mis_by_chrom[r["chrom"]].append(r)

    # decide plotting order: keep input chrom order, but only those with fai lengths
    plot_chroms = [c for c in chroms if c in genome_lengths]
    if not plot_chroms:
        # fall back: intersect of fai and chrom_to_subg
        plot_chroms = [c for c in genome_lengths.keys() if c in chrom_to_subg]

    # sort mismatch windows within chrom
    for c in mis_by_chrom:
        mis_by_chrom[c].sort(key=lambda r: (int(r["win_start"]), int(r["win_end"])))

    # handle very large lists by paging (rare, but nice)
    pages = [plot_chroms[i:i+max_chroms_per_page] for i in range(0, len(plot_chroms), max_chroms_per_page)]

    # Make multi-page PDF if needed
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(out_pdf) as pdf:
        for page_idx, chrom_list in enumerate(pages, start=1):
            fig_h = max(6.0, 0.18 * len(chrom_list) + 2.0)
            plt.figure(figsize=(12, fig_h))

            y_positions = np.arange(len(chrom_list))[::-1]  # top to bottom
            max_len = max(genome_lengths[c] for c in chrom_list)

            # Draw base chromosomes
            for y, chrom in zip(y_positions, chrom_list):
                L = genome_lengths[chrom]
                subg = chrom_to_subg.get(chrom, -1)
                base_color = (0.7, 0.7, 0.7, 1.0) if subg == -1 else subg_colors[subg]

                # base bar
                plt.plot([0, L], [y, y], linewidth=10, solid_capstyle="butt", color=base_color, zorder=1)

                # overlay mismatch windows
                for r in mis_by_chrom.get(chrom, []):
                    ws = int(r["win_start"])
                    we = int(r["win_end"])
                    # clamp to chromosome length
                    ws = max(0, min(ws, L))
                    we = max(0, min(we, L))
                    if we <= ws:
                        continue

                    # best enriched subg index
                    bes = r["best_enriched_subgenome"]
                    if isinstance(bes, str) and bes.startswith("subg"):
                        try:
                            best_g = int(bes.replace("subg", ""))
                        except ValueError:
                            best_g = -1
                    else:
                        best_g = -1

                    seg_color = (0.0, 0.0, 0.0, 0.8) if best_g == -1 else subg_colors.get(best_g, (0.0, 0.0, 0.0, 0.8))

                    # colored segment (exchange)
                    plt.plot([ws, we], [y, y], linewidth=10, solid_capstyle="butt",
                             color=seg_color, zorder=2)
                    # outline to emphasize
                    plt.plot([ws, we], [y, y], linewidth=12, solid_capstyle="butt",
                             color=(0, 0, 0, 0.6), zorder=1.9)

            plt.yticks(y_positions, chrom_list, fontsize=9)
            plt.xticks(fontsize=9)
            plt.xlim(0, max_len * 1.02)
            plt.xlabel("Genomic position (bp)")
            title = "Chromosome subgenome assignments with mismatch (exchange/switch-error) windows"
            if len(pages) > 1:
                title += f" (page {page_idx}/{len(pages)})"
            plt.title(title)

            # legend: subgenomes
            legend_items = []
            for g in range(n_subgenomes):
                legend_items.append(Patch(facecolor=subg_colors[g], edgecolor="none", label=f"chrom=subg{g}"))

            # exchange note
            legend_items.append(Patch(facecolor=(1, 1, 1, 0), edgecolor="black", label="mismatch windows overlaid (colored by best-enriched subg)"))

            plt.legend(handles=legend_items, frameon=False, ncol=4, fontsize=9, loc="upper right")
            plt.tight_layout()
            pdf.savefig()
            plt.close()


def main():
    ap = argparse.ArgumentParser(
        description="LTR-only SubPhaser-style subgenome phasing + LTR subphasing (+K2P) + window scan + optional chrom plot."
    )
    ap.add_argument("--ltr_fasta", required=True, help="Intact LTR-RT fasta.")
    ap.add_argument("--homoeolog_config", required=True, help="TSV/space-separated homoeologous chromosome sets.")
    ap.add_argument("--k2p_tsv", required=True, help="K2P table; col1=LTR header key, col11=K2P.")
    ap.add_argument("--outdir", required=True, help="Output directory.")
    ap.add_argument("--n_subgenomes", type=int, required=True, help="Number of subgenomes (K for k-means).")

    ap.add_argument("-k", type=int, default=15, help="k-mer size (default: 15).")
    ap.add_argument("--canonical", action="store_true", help="Use canonical kmers (default off; turn on to match SubPhaser).")
    ap.add_argument("-L", type=int, default=3, help="Min k-mer count to dump from jellyfish (default: 3).")
    ap.add_argument("--jf_size", type=int, default=100000000, help="Jellyfish -s hash size (default: 100000000).")
    ap.add_argument("--threads", type=int, default=4, help="Threads for jellyfish (default: 4).")

    ap.add_argument("--fmin", type=float, default=2.0, help="Differential k-mer fold threshold r0/r1 (default: 2).")
    ap.add_argument("--qmin", type=int, default=200, help="Min total genome count for k-mer (default: 200).")

    ap.add_argument("--bootstrap", type=int, default=1000, help="Bootstrap replicates (default: 1000).")
    ap.add_argument("--resample_frac", type=float, default=0.5, help="Fraction of kmers to resample per bootstrap (default: 0.5).")
    ap.add_argument("--seed", type=int, default=1, help="Random seed (default: 1).")

    ap.add_argument("--subg_kmer_p", type=float, default=0.05, help="p-value cutoff for subgenome-specific k-mers (default: 0.05).")
    ap.add_argument("--ltr_alpha", type=float, default=0.05, help="q-value cutoff for LTR enrichment calls (default: 0.05).")

    ap.add_argument("--heatmap_kmers", type=int, default=10000, help="Max kmers to plot in heatmap (default: 10000).")

    ap.add_argument("--k2p_col", type=int, default=11, help="1-based column index for K2P in --k2p_tsv (default: 11).")
    ap.add_argument("--k2p_window", type=float, default=0.0005, help="K2P window/bin width for density hist (default: 0.0005).")

    # Window scan parameters
    ap.add_argument("--window_size", type=int, default=1_000_000, help="Genomic window size in bp (default: 1000000).")
    ap.add_argument("--window_alpha", type=float, default=0.05, help="q-value cutoff for window enrichment (default: 0.05).")
    ap.add_argument("--min_ltrs_per_window", type=int, default=3, help="Min LTRs required in a window to test enrichment (default: 3).")
    ap.add_argument("--min_total_hits_per_window", type=int, default=100, help="Min summed k-mer hits in a window to test (default: 100).")

    # NEW: genome fai for chromosome plot
    ap.add_argument("--genome_fai", default=None, help="Genome .fai file (optional). If provided, outputs chromosome_plot.pdf")

    args = ap.parse_args()

    ensure_dir(args.outdir)
    workdir = os.path.join(args.outdir, "work")
    ensure_dir(workdir)

    records, chrom_to_indices = read_ltr_fasta(args.ltr_fasta)
    homoeolog_sets = read_homoeolog_config(args.homoeolog_config)
    k2p_map = read_k2p_map(args.k2p_tsv, k2p_col_1based=args.k2p_col)

    chroms = sorted(chrom_to_indices.keys())
    if len(chroms) < 2:
        raise SystemExit("ERROR: Need LTRs from at least 2 chromosomes to phase.")

    chrom_fas = write_chrom_concat_fastas(records, chrom_to_indices, workdir, sep_n=50)
    total_kmers = compute_total_kmers_per_chrom(records, chrom_to_indices, args.k)

    kmers, counts = build_kmer_count_matrix(
        chroms=chroms,
        chrom_fas=chrom_fas,
        total_kmers=total_kmers,
        k=args.k,
        jf_size=args.jf_size,
        canonical=args.canonical,
        Lmin=args.L,
        threads=args.threads,
        workdir=workdir,
    )

    if counts.shape[0] == 0:
        raise SystemExit("ERROR: No k-mers passed the dump threshold. Try lowering -L or using smaller k.")

    denoms = np.array([max(1, total_kmers.get(c, 1)) for c in chroms], dtype=np.float64)
    ratios = counts.astype(np.float64) / denoms[None, :]

    diff_mask = identify_differential_kmers(
        kmers=kmers,
        counts=counts,
        chroms=chroms,
        total_kmers=total_kmers,
        homoeolog_sets=homoeolog_sets,
        fmin=args.fmin,
        qmin=args.qmin,
    )
    diff_idx = np.where(diff_mask)[0]
    if diff_idx.size < 50:
        raise SystemExit(
            f"ERROR: Too few differential k-mers ({diff_idx.size}). "
            "Try lowering --qmin, lowering --fmin, or lowering -L, or increasing input LTRs."
        )

    diff_kmers = [kmers[i] for i in diff_idx]
    diff_ratios = ratios[diff_idx, :]

    M1 = zscore_rows(diff_ratios)

    chrom_features = M1.T
    labels, stability = bootstrap_kmeans(
        chrom_features=chrom_features,
        n_clusters=args.n_subgenomes,
        n_boot=args.bootstrap,
        resample_frac=args.resample_frac,
        seed=args.seed,
    )

    chrom_out = os.path.join(args.outdir, "chromosome_phasing.tsv")
    with open(chrom_out, "w") as out:
        out.write("chromosome\tsubgenome\tbootstrap_stability\ttotal_LTR_kmers\n")
        for c, lab, st in zip(chroms, labels.tolist(), stability.tolist()):
            out.write(f"{c}\tsubg{lab}\t{st:.3f}\t{total_kmers.get(c,0)}\n")

    pca_pdf = os.path.join(args.outdir, "PCA.pdf")
    heat_pdf = os.path.join(args.outdir, "heatmap.pdf")
    plot_pca(chrom_features, chroms, labels, pca_pdf)
    plot_heatmap(M1, chroms, heat_pdf, args.heatmap_kmers, args.seed)

    spec_masks = call_subgenome_specific_kmers(diff_ratios, labels, pval=args.subg_kmer_p)

    subg_sets: Dict[int, set] = {}
    for g, m in spec_masks.items():
        s = set(np.array(diff_kmers, dtype=object)[m].tolist())
        subg_sets[g] = s
        outk = os.path.join(args.outdir, f"subgenome_specific_kmers.subg{g}.txt.gz")
        with gzip.open(outk, "wt") as gz:
            for kmer in sorted(s):
                gz.write(kmer + "\n")

    chrom_to_subg = {c: int(lab) for c, lab in zip(chroms, labels.tolist())}
    subg_ids = sorted(subg_sets.keys())

    ltr_out = os.path.join(args.outdir, "ltr_subphasing.tsv")
    k2p_by_group: Dict[str, List[float]] = defaultdict(list)

    with open(ltr_out, "w") as out:
        header = [
            "chrom", "start", "end", "class", "superfamily", "family",
            "chrom_subgenome",
            "best_enriched_subgenome", "best_p", "best_q",
            "call",
            "K2P",
        ]
        for g in subg_ids:
            header.append(f"hits_subg{g}")
        out.write("\t".join(header) + "\n")

        missing_k2p = 0

        for rec in records:
            hits = count_ltr_kmer_hits(rec.seq, args.k, args.canonical, subg_sets)

            pvals = []
            for g in subg_ids:
                c00 = hits[g]
                c01 = sum(hits[x] for x in subg_ids if x != g)
                bg_g = len(subg_sets[g])
                bg_o = sum(len(subg_sets[x]) for x in subg_ids if x != g)
                _, p = fisher_exact([[c00, c01], [bg_g, bg_o]], alternative="greater")
                pvals.append(p)

            qvals = bh_adjust(pvals)
            best_i = int(np.argmin(qvals))
            best_g = subg_ids[best_i]
            best_p = pvals[best_i]
            best_q = qvals[best_i]

            chrom_subg = chrom_to_subg.get(rec.chrom, -1)

            if best_q < args.ltr_alpha and best_g == chrom_subg:
                call = "PRE_polyploidy_subgenome_specific"
            elif best_q < args.ltr_alpha and best_g != chrom_subg and chrom_subg != -1:
                call = "EXCHANGE_OR_ERROR"
            else:
                call = "POST_polyploidy_shared_or_mixed"

            key = rec.key()
            k2p_val = k2p_map.get(key, float("nan"))
            if not np.isfinite(k2p_val):
                missing_k2p += 1
            else:
                k2p_by_group[call].append(k2p_val)

            row = [
                rec.chrom, str(rec.start), str(rec.end),
                rec.klass, rec.superfamily, rec.family,
                f"subg{chrom_subg}" if chrom_subg != -1 else "NA",
                f"subg{best_g}",
                f"{best_p:.3e}",
                f"{best_q:.3e}",
                call,
                f"{k2p_val:.6f}" if np.isfinite(k2p_val) else "NA",
            ]
            for g in subg_ids:
                row.append(str(hits[g]))
            out.write("\t".join(row) + "\n")

    k2p_pdf = os.path.join(args.outdir, "K2P_density_windows.pdf")
    plot_k2p_window_density(k2p_by_group, k2p_pdf, window=args.k2p_window)

    all_windows, mismatches, missing_k2p_used = window_enrichment_scan(
        records=records,
        subg_sets=subg_sets,
        chrom_to_subg=chrom_to_subg,
        k=args.k,
        canonical=args.canonical,
        window_size=args.window_size,
        window_alpha=args.window_alpha,
        min_ltrs_per_window=args.min_ltrs_per_window,
        min_total_hits_per_window=args.min_total_hits_per_window,
        k2p_map=k2p_map,
    )

    win_all_tsv = os.path.join(args.outdir, "window_enrichment.tsv")
    win_mis_tsv = os.path.join(args.outdir, "window_enrichment.mismatches.tsv")
    write_window_tables(all_windows, mismatches, win_all_tsv, win_mis_tsv, subg_ids=subg_ids)

    chrom_plot_pdf = None
    if args.genome_fai:
        genome_lengths = read_fai_lengths(args.genome_fai)
        chrom_plot_pdf = os.path.join(args.outdir, "chromosome_plot.pdf")
        plot_chromosomes_with_exchanges(
            genome_lengths=genome_lengths,
            chroms=chroms,
            chrom_to_subg=chrom_to_subg,
            mismatches=mismatches,
            n_subgenomes=args.n_subgenomes,
            out_pdf=chrom_plot_pdf,
        )

    print("DONE")
    print(f"- {chrom_out}")
    print(f"- {ltr_out}")
    print(f"- {pca_pdf}")
    print(f"- {heat_pdf}")
    print(f"- {k2p_pdf}")
    print(f"- {win_all_tsv}")
    print(f"- {win_mis_tsv}   (n={len(mismatches)})")
    if chrom_plot_pdf:
        print(f"- {chrom_plot_pdf}")

    if missing_k2p_used:
        print(f"NOTE: {missing_k2p_used} LTRs were missing K2P during window scan (ignored for mean_k2p).")
    if missing_k2p:
        print(f"NOTE: {missing_k2p} LTRs had no matching K2P entry (K2P=NA). "
              f"Expected key format: chrom:start-end#class/superfamily/family")


if __name__ == "__main__":
    main()
