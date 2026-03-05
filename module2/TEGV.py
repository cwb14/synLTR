#!/usr/bin/env python3
"""
te_viewer.py — Interactive nested TE + gene genome browser
Generates a self-contained HTML file with D3.js visualization.

Now supports --bam and --bam-labels to add per-sample:
  • Coverage barplots over exon/TE positions
  • Collapsed transcript tracks (spliced alignments shown as exon boxes + intron lines)

Usage:
    python TEGV.py \
        --gff3 genes.gff3 \
        --te-fastas level0.fa level1.fa level2.fa \
        --ltr-divergence level0.div level1.div level2.div \
        --domains level0.domains level1.domains level2.domains \
        --bam wt.bam met1.bam ddm1.bam \
        --bam-labels wt met1 ddm1 \
        --output viewer.html

python TEGV.py --gff3 Col0_ltr_r1.work/Col0_ltr_r1.genic.gff  --te-fastas Col0_ltr_r1.ltrharvest.full_length.dedup.fa.rexdb-plant.cls.lib.fa Col0_ltr_r2.ltrharvest.full_length.dedup.fa.rexdb-plant.cls.lib.fa Col0_ltr_r3.ltrharvest.full_length.dedup.fa.rexdb-plant.cls.lib.fa --ltr-divergence Col0_ltr_r1_kmer2ltr_dedup Col0_ltr_r2_kmer2ltr_dedup  Col0_ltr_r3_kmer2ltr_dedup --domains  Col0_ltr_r1.work/Col0_ltr_r1.ltrtools.intact_for_tesorter.fa.rexdb-plant.dom.gff3 Col0_ltr_r2.work/Col0_ltr_r2.ltrtools.intact_for_tesorter.fa.rexdb-plant.dom.gff3 Col0_ltr_r3.work/Col0_ltr_r3.ltrtools.intact_for_tesorter.fa.rexdb-plant.dom.gff3 --bam ltr_rt_validation/bam/wt.bam ltr_rt_validation/bam/met1.bam ltr_rt_validation/bam/ddm1.bam --bam-labels wt met1 ddm1 --output Col0_viewer_new.html --bin-size 200 --threads 200
"""

import argparse
import json
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path


# ─────────────────────────────────────────────
#  Data classes
# ─────────────────────────────────────────────

@dataclass
class ProteinDomain:
    start: int
    end: int
    name: str
    gene: str
    evalue: float
    probability: float
    strand: str


@dataclass
class TEElement:
    chrom: str
    start: int
    end: int
    family: str
    superfamily: str
    te_class: str
    level: int
    ltr_len: int = 0
    k2p_dist: float = 0.0
    k2p_time: float = 0.0
    has_divergence: bool = False
    domains: list = field(default_factory=list)
    children: list = field(default_factory=list)
    visible_frags: list = field(default_factory=list)

    @property
    def key(self):
        return f"{self.chrom}:{self.start}-{self.end}#{self.superfamily}/{self.family}"

    @property
    def length(self):
        return self.end - self.start


@dataclass
class Exon:
    start: int
    end: int
    feature_type: str


@dataclass
class Gene:
    chrom: str
    start: int
    end: int
    name: str
    strand: str
    gene_id: str
    exons: list = field(default_factory=list)
    targets: list = field(default_factory=list)


# ─────────────────────────────────────────────
#  Parsers
# ─────────────────────────────────────────────

def parse_fasta_headers(fasta_path, level):
    elements = []
    with open(fasta_path, buffering=1 << 20) as f:
        for line in f:
            if not line or line[0] != '>':
                continue
            header = line[1:].rstrip().split()[0]
            if '#' not in header:
                print(f"Warning: skipping malformed header: {header}", file=sys.stderr)
                continue
            coord_part, classif = header.split('#', 1)
            m = re.match(r'^(.+):(\d+)-(\d+)$', coord_part)
            if not m:
                print(f"Warning: cannot parse coords in: {header}", file=sys.stderr)
                continue
            chrom = m.group(1)
            start = int(m.group(2))
            end = int(m.group(3))
            parts = classif.split('/')
            te_class = parts[0] if len(parts) > 0 else 'Unknown'
            superfamily_parts = parts[:2] if len(parts) >= 2 else parts
            superfamily = '/'.join(superfamily_parts)
            family = parts[-1] if len(parts) >= 2 else parts[0]
            elements.append(TEElement(
                chrom=chrom, start=start, end=end,
                family=family, superfamily=superfamily,
                te_class=te_class, level=level
            ))
    return elements


def parse_ltr_divergence(div_path, te_map):
    with open(div_path, buffering=1 << 20) as f:
        for line in f:
            if not line or line[0] == '#':
                continue
            cols = line.split()
            if not cols or len(cols) < 12:
                continue
            key = cols[0]
            try:
                ltr_len = int(cols[1])
                k2p_dist = float(cols[10])
                k2p_time = float(cols[11])
            except (ValueError, IndexError):
                continue
            if key in te_map:
                te_map[key].ltr_len = ltr_len
                te_map[key].k2p_dist = k2p_dist
                te_map[key].k2p_time = k2p_time
                te_map[key].has_divergence = True


def parse_domains(domain_path, te_map):
    coord_index = {(te.chrom, te.start, te.end): te for te in te_map.values()}
    _id_re   = re.compile(r'ID=([^;\t]+)')
    _gene_re = re.compile(r'gene=([^;\t]+)')
    _name_re = re.compile(r'Name=([^;\t]+)')
    _eval_re = re.compile(r'evalue=([^;\t]+)')
    _prob_re = re.compile(r'probability=([^;\t]+)')
    _coord_re = re.compile(r':(\d+)-(\d+)[|#]')

    with open(domain_path, buffering=1 << 20) as f:
        for line in f:
            if not line or line[0] == '#':
                continue
            cols = line.split('\t')
            if len(cols) < 9:
                continue
            chrom      = cols[0]
            feat_start = int(cols[3])
            feat_end   = int(cols[4])
            strand     = cols[6]
            attrs      = cols[8]

            id_m = _id_re.search(attrs)
            if not id_m:
                continue
            coord_m = _coord_re.search(id_m.group(1))
            if not coord_m:
                continue
            te_start = int(coord_m.group(1))
            te_end   = int(coord_m.group(2))

            te = coord_index.get((chrom, te_start, te_end))
            if te is None:
                continue

            gene_m = _gene_re.search(attrs)
            name_m = _name_re.search(attrs)
            eval_m = _eval_re.search(attrs)
            prob_m = _prob_re.search(attrs)

            gene_name = (gene_m.group(1) if gene_m else
                         name_m.group(1) if name_m else 'domain')
            try:
                evalue = float(eval_m.group(1)) if eval_m else 0.0
                prob   = float(prob_m.group(1)) if prob_m else 0.0
            except ValueError:
                evalue = prob = 0.0

            te.domains.append(ProteinDomain(
                start=feat_start, end=feat_end,
                name=gene_name, gene=gene_name,
                evalue=evalue, probability=prob, strand=strand
            ))


def parse_gff3(gff3_path):
    genes = {}
    mrna_to_gene = {}

    with open(gff3_path, buffering=1 << 20) as f:
        for line in f:
            if not line or line[0] == '#':
                continue
            cols = line.split('\t')
            if len(cols) < 9:
                continue
            chrom, source, feature, start, end, score, strand, phase, attrs = cols
            start, end = int(start), int(end)
            attr_dict = {}
            for attr in attrs.split(';'):
                eq = attr.find('=')
                if eq > 0:
                    attr_dict[attr[:eq].strip()] = attr[eq+1:].strip()

            gid    = attr_dict.get('ID', '')
            parent = attr_dict.get('Parent', '')
            target_raw  = attr_dict.get('Target', '')
            target_name = target_raw.split()[0] if target_raw else ''
            name = target_name or attr_dict.get('Name', '') or gid

            if feature == 'gene':
                g = Gene(chrom=chrom, start=start, end=end,
                         name=name, strand=strand, gene_id=gid)
                if target_name:
                    g.targets = [target_name]
                genes[gid] = g
            elif feature == 'mRNA':
                if parent and parent in genes:
                    mrna_to_gene[gid] = parent
                    if target_name and target_name not in genes[parent].targets:
                        genes[parent].targets.append(target_name)
                        if not genes[parent].name or genes[parent].name == genes[parent].gene_id:
                            genes[parent].name = target_name
                elif parent:
                    mrna_to_gene[gid] = parent
                    g = Gene(chrom=chrom, start=start, end=end,
                             name=name, strand=strand, gene_id=parent)
                    g.targets = [target_name] if target_name else []
                    genes[parent] = g
                else:
                    mrna_to_gene[gid] = gid
                    g = Gene(chrom=chrom, start=start, end=end,
                             name=name, strand=strand, gene_id=gid)
                    g.targets = [target_name] if target_name else []
                    genes[gid] = g
            elif feature in ('exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR',
                              '5UTR', '3UTR', 'UTR', 'stop_codon'):
                target_gene = mrna_to_gene.get(parent, parent)
                if target_gene in genes and feature != 'stop_codon':
                    genes[target_gene].exons.append(
                        Exon(start=start, end=end, feature_type=feature)
                    )

    raw = list(genes.values())
    return _collapse_overlapping_genes(raw)


def _collapse_overlapping_genes(genes):
    from collections import defaultdict
    buckets = defaultdict(list)
    for g in genes:
        buckets[(g.chrom, g.strand)].append(g)

    merged = []
    for (chrom, strand), group in buckets.items():
        group.sort(key=lambda g: g.start)
        clusters = []
        cur_cluster = []
        cur_end = -1
        for g in group:
            if cur_cluster and g.start > cur_end:
                clusters.append(cur_cluster)
                cur_cluster = []
                cur_end = -1
            cur_cluster.append(g)
            if g.end > cur_end:
                cur_end = g.end
        if cur_cluster:
            clusters.append(cur_cluster)

        for cluster in clusters:
            rep = cluster[0]
            all_targets = []
            seen_targets = set()
            all_exons = []
            seen_exons = set()
            for g in cluster:
                for t in g.targets:
                    if t not in seen_targets:
                        all_targets.append(t)
                        seen_targets.add(t)
                for ex in g.exons:
                    key = (ex.start, ex.end, ex.feature_type)
                    if key not in seen_exons:
                        all_exons.append(ex)
                        seen_exons.add(key)
            merged_gene = Gene(
                chrom=chrom,
                start=min(g.start for g in cluster),
                end=max(g.end   for g in cluster),
                name=all_targets[0] if all_targets else rep.name,
                strand=strand,
                gene_id=rep.gene_id,
                exons=all_exons,
                targets=all_targets,
            )
            merged.append(merged_gene)

    return merged


# ─────────────────────────────────────────────
#  Nesting logic
# ─────────────────────────────────────────────

def build_nesting(all_elements):
    by_level = defaultdict(list)
    for el in all_elements:
        by_level[el.level].append(el)

    max_level = max(by_level.keys()) if by_level else 0

    children_by_chrom = defaultdict(lambda: defaultdict(list))
    for level, elems in by_level.items():
        for el in elems:
            children_by_chrom[level][el.chrom].append(el)

    for parent_level in range(1, max_level + 1):
        child_level = parent_level - 1
        for parent in by_level[parent_level]:
            for child in children_by_chrom[child_level][parent.chrom]:
                if child.start >= parent.start and child.end <= parent.end:
                    parent.children.append(child)

    for el in all_elements:
        compute_visible_frags(el)

    return all_elements


def compute_visible_frags(el):
    if not el.children:
        el.visible_frags = [(el.start, el.end)]
        return
    children_sorted = sorted(el.children, key=lambda c: c.start)
    frags = []
    cursor = el.start
    for child in children_sorted:
        if child.start > cursor:
            frags.append((cursor, child.start))
        cursor = max(cursor, child.end)
    if cursor < el.end:
        frags.append((cursor, el.end))
    el.visible_frags = frags if frags else []


# ─────────────────────────────────────────────
#  BAM processing  (fast, parallel)
# ─────────────────────────────────────────────

def _merge_regions(regions_by_chrom):
    """Merge overlapping / nearby regions per chromosome."""
    merged = {}
    for chrom, regs in regions_by_chrom.items():
        regs.sort()
        m = []
        for s, e in regs:
            if m and s <= m[-1][1] + 1000:
                m[-1] = (m[-1][0], max(m[-1][1], e))
            else:
                m.append((s, e))
        merged[chrom] = m
    return merged


def _process_bam_chrom(bam_path, chrom, regions, bin_size):
    """
    Worker function: process one (BAM, chrom) pair.
    Each worker opens its own file handle (required for multiprocessing).

    Returns (chrom, cov_bins, tx_clusters).

    FAST coverage strategy:
      • One count_coverage() call per merged region (can span 100 kb+).
      • numpy sums four base-count arrays → per-base depth vector.
      • reshape/mean gives binned depth in one vectorised step.

    FAST transcript strategy:
      • read.get_blocks() (C-level) instead of Python CIGAR walk.
      • Tuples instead of dicts for intermediate storage.
    """
    import pysam
    import numpy as np

    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
    except Exception:
        return (chrom, [], [])

    if chrom not in set(bam.references):
        bam.close()
        return (chrom, [], [])

    chrom_bins = []
    all_raw_reads = []   # list of (start, end, strand, blocks_tuple)

    for reg_start, reg_end in regions:
        reg_len = reg_end - reg_start

        # ── Coverage: single call, numpy binning ──
        try:
            counts = bam.count_coverage(
                chrom, reg_start, reg_end, read_callback='all'
            )
            # counts is a tuple of 4 array.array objects (A, C, G, T)
            # array.array typecode varies by platform ('L' = 4 or 8 bytes)
            # so use np.array() which handles any type correctly
            depth = np.zeros(reg_len, dtype=np.float64)
            for arr in counts:
                depth += np.array(arr, dtype=np.float64)

            # Bin: trim to multiple of bin_size, reshape, mean
            n_full_bins = reg_len // bin_size
            if n_full_bins > 0:
                trimmed = depth[:n_full_bins * bin_size]
                binned = trimmed.reshape(n_full_bins, bin_size).mean(axis=1)
                for i in range(n_full_bins):
                    d = float(binned[i])
                    if d > 0.01:
                        chrom_bins.append({
                            "s": reg_start + i * bin_size,
                            "e": reg_start + (i + 1) * bin_size,
                            "d": round(d, 2)
                        })
            # Remainder bin
            remainder = reg_len - n_full_bins * bin_size
            if remainder > 0:
                d = float(depth[n_full_bins * bin_size:].mean())
                if d > 0.01:
                    chrom_bins.append({
                        "s": reg_start + n_full_bins * bin_size,
                        "e": reg_end,
                        "d": round(d, 2)
                    })
        except Exception:
            pass

        # ── Transcripts: extract spliced reads ──
        try:
            for read in bam.fetch(chrom, reg_start, reg_end):
                if read.flag & 0x904:          # unmapped | secondary | supplementary
                    continue
                if read.mapping_quality < 5:
                    continue
                blocks = read.get_blocks()     # C-level, returns list of (start, end)
                if not blocks:
                    continue
                strand = '-' if read.is_reverse else '+'
                all_raw_reads.append((blocks[0][0], blocks[-1][1], strand, blocks))
        except Exception:
            pass

    bam.close()

    # ── Collapse transcripts ──
    tx_clusters = _collapse_transcripts_fast(all_raw_reads, chrom)

    return (chrom, chrom_bins, tx_clusters)


def _collapse_transcripts_fast(raw_reads, chrom):
    """
    Collapse overlapping spliced reads into representative transcript clusters.
    Input: list of (start, end, strand, blocks) tuples.
    Returns: list of dicts ready for JSON serialisation.
    """
    if not raw_reads:
        return []

    # Group by strand
    by_strand = defaultdict(list)
    for start, end, strand, blocks in raw_reads:
        by_strand[strand].append((start, end, blocks))

    results = []
    for strand, reads in by_strand.items():
        reads.sort()  # by (start, end, ...)
        # Sweep-line clustering
        clusters = []
        cur_reads = []
        cur_end = -1
        for start, end, blocks in reads:
            if cur_reads and start > cur_end:
                clusters.append(cur_reads)
                cur_reads = []
                cur_end = -1
            cur_reads.append(blocks)
            if end > cur_end:
                cur_end = end
        if cur_reads:
            clusters.append(cur_reads)

        for cluster in clusters:
            depth = len(cluster)
            # Merge all exon blocks
            all_blocks = []
            for blocks in cluster:
                all_blocks.extend(blocks)
            all_blocks.sort()
            merged = []
            for s, e in all_blocks:
                if merged and s <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], e))
                else:
                    merged.append((s, e))
            results.append({
                "c": chrom,
                "s": merged[0][0],
                "e": merged[-1][1],
                "st": strand,
                "b": [[s, e] for s, e in merged],
                "d": depth
            })

    return results


def compute_coverage_and_transcripts(bam_paths, bam_labels, genes, all_elements,
                                     bin_size=50, threads=1):
    """
    Parallel BAM processing.  Each (BAM, chrom) pair is dispatched as an
    independent job to a ProcessPoolExecutor.

    Returns:
      coverage_data:   { label: { chrom: [ {s, e, d}, ... ] } }
      transcript_data: { label: [ {c, s, e, st, b, d}, ... ] }
    """
    from concurrent.futures import ProcessPoolExecutor, as_completed
    import time

    # Collect regions of interest per chrom
    regions_by_chrom = defaultdict(list)
    for g in genes:
        regions_by_chrom[g.chrom].append((g.start, g.end))
    for te in all_elements:
        regions_by_chrom[te.chrom].append((te.start, te.end))

    merged_regions = _merge_regions(regions_by_chrom)

    total_region_bp = sum(
        sum(e - s for s, e in regs) for regs in merged_regions.values()
    )
    print(f"  Regions of interest: {len(merged_regions)} chroms, "
          f"{sum(len(r) for r in merged_regions.values())} regions, "
          f"{total_region_bp/1e6:.1f} Mb total",
          file=sys.stderr)

    coverage_data = {}
    transcript_data = {}

    # Build job list: (bam_path, label, chrom, regions)
    jobs = []
    for bam_path, label in zip(bam_paths, bam_labels):
        for chrom, regions in merged_regions.items():
            jobs.append((bam_path, label, chrom, regions))

    print(f"  Dispatching {len(jobs)} jobs across {threads} workers",
          file=sys.stderr)

    # Initialise result containers
    for label in bam_labels:
        coverage_data[label] = {}
        transcript_data[label] = []

    t0 = time.time()
    completed = 0

    if threads <= 1:
        # Sequential — avoids pickling overhead for single-thread
        for bam_path, label, chrom, regions in jobs:
            ch, cov_bins, tx_clusters = _process_bam_chrom(
                bam_path, chrom, regions, bin_size
            )
            if cov_bins:
                coverage_data[label][ch] = cov_bins
            if tx_clusters:
                transcript_data[label].extend(tx_clusters)
            completed += 1
            elapsed = time.time() - t0
            print(f"\r  [{completed}/{len(jobs)}] {label} {chrom} "
                  f"({len(cov_bins)} cov bins, {len(tx_clusters)} tx) "
                  f"[{elapsed:.0f}s]",
                  end='', file=sys.stderr)
        print(file=sys.stderr)
    else:
        with ProcessPoolExecutor(max_workers=threads) as pool:
            future_map = {}
            for bam_path, label, chrom, regions in jobs:
                fut = pool.submit(
                    _process_bam_chrom, bam_path, chrom, regions, bin_size
                )
                future_map[fut] = (label, chrom)

            for fut in as_completed(future_map):
                label, chrom = future_map[fut]
                try:
                    ch, cov_bins, tx_clusters = fut.result()
                except Exception as ex:
                    print(f"\n  ERROR {label}/{chrom}: {ex}", file=sys.stderr)
                    continue
                if cov_bins:
                    coverage_data[label][ch] = cov_bins
                if tx_clusters:
                    transcript_data[label].extend(tx_clusters)
                completed += 1
                elapsed = time.time() - t0
                print(f"\r  [{completed}/{len(jobs)}] {label} {chrom} "
                      f"({len(cov_bins)} cov bins, {len(tx_clusters)} tx) "
                      f"[{elapsed:.0f}s]",
                      end='', file=sys.stderr)
            print(file=sys.stderr)

    for label in bam_labels:
        n_cov = sum(len(v) for v in coverage_data[label].values())
        n_tx = len(transcript_data[label])
        print(f"  {label}: {n_cov} coverage bins, {n_tx} transcript clusters",
              file=sys.stderr)

    return coverage_data, transcript_data


# ─────────────────────────────────────────────
#  Serialisation for JS
# ─────────────────────────────────────────────

def te_to_dict(el):
    return {
        "chrom": el.chrom,
        "start": el.start,
        "end": el.end,
        "family": el.family,
        "superfamily": el.superfamily,
        "te_class": el.te_class,
        "level": el.level,
        "ltr_len": el.ltr_len,
        "k2p_dist": el.k2p_dist,
        "k2p_time": el.k2p_time,
        "has_divergence": el.has_divergence,
        "visible_frags": el.visible_frags,
        "domains": [
            {
                "start": d.start, "end": d.end,
                "name": d.name, "gene": d.gene,
                "evalue": d.evalue, "probability": d.probability,
                "strand": d.strand
            } for d in el.domains
        ],
        "key": el.key,
    }


def gene_to_dict(g):
    return {
        "chrom": g.chrom,
        "start": g.start,
        "end": g.end,
        "name": g.name,
        "strand": g.strand,
        "gene_id": g.gene_id,
        "exons": [{"start": e.start, "end": e.end, "type": e.feature_type}
                  for e in g.exons],
        "targets": g.targets,
    }


# ─────────────────────────────────────────────
#  HTML template
# ─────────────────────────────────────────────

HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>TE Genome Viewer</title>
<script src="https://cdnjs.cloudflare.com/ajax/libs/d3/7.8.5/d3.min.js"></script>
<style>
  @import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500;600&family=IBM+Plex+Sans:wght@300;400;500&display=swap');

  :root {
    --bg: #0e1117;
    --surface: #161b27;
    --surface2: #1e2535;
    --border: #2a3348;
    --text: #c8d3e8;
    --text-dim: #5a6a88;
    --text-bright: #e8f0ff;
    --accent: #4a9eff;
    --gene-color: #5b9bd5;
    --gene-exon: #7bb8f0;
    --gene-utr: #3a6a9a;
    --gene-cds: #9acfff;
  }

  * { box-sizing: border-box; margin: 0; padding: 0; }

  body {
    background: var(--bg);
    color: var(--text);
    font-family: 'IBM Plex Sans', sans-serif;
    font-size: 13px;
    height: 100vh;
    display: flex;
    flex-direction: column;
    overflow: hidden;
  }

  header {
    background: var(--surface);
    border-bottom: 1px solid var(--border);
    padding: 10px 20px;
    display: flex;
    align-items: center;
    gap: 20px;
    flex-wrap: wrap;
    flex-shrink: 0;
  }

  .brand {
    font-family: 'IBM Plex Mono', monospace;
    font-weight: 600;
    font-size: 14px;
    color: var(--accent);
    letter-spacing: 0.05em;
    white-space: nowrap;
  }

  .controls {
    display: flex;
    align-items: center;
    gap: 10px;
    flex-wrap: wrap;
  }

  label { color: var(--text-dim); font-size: 11px; text-transform: uppercase; letter-spacing: 0.08em; }

  select, input[type=number], input[type=text] {
    background: var(--surface2);
    border: 1px solid var(--border);
    color: var(--text-bright);
    font-family: 'IBM Plex Mono', monospace;
    font-size: 12px;
    padding: 5px 9px;
    border-radius: 4px;
    outline: none;
    transition: border-color 0.15s;
  }
  select:focus, input:focus { border-color: var(--accent); }
  input[type=number] { width: 100px; }

  button {
    background: var(--surface2);
    border: 1px solid var(--border);
    color: var(--text);
    font-family: 'IBM Plex Sans', sans-serif;
    font-size: 12px;
    padding: 5px 12px;
    border-radius: 4px;
    cursor: pointer;
    transition: background 0.15s, border-color 0.15s;
  }
  button:hover { background: var(--border); border-color: var(--accent); color: var(--text-bright); }
  button.primary { background: var(--accent); border-color: var(--accent); color: #fff; font-weight: 500; }
  button.primary:hover { background: #6ab4ff; }

  .sep { width: 1px; height: 24px; background: var(--border); }

  #minimap-container {
    background: var(--surface);
    border-bottom: 1px solid var(--border);
    padding: 6px 20px;
    flex-shrink: 0;
  }
  #minimap { width: 100%; height: 28px; cursor: crosshair; }

  #viewer-wrap {
    flex: 1;
    overflow: hidden;
    position: relative;
  }

  #viewer-svg { width: 100%; height: 100%; cursor: grab; }
  #viewer-svg:active { cursor: grabbing; }

  #tooltip {
    position: fixed;
    background: var(--surface2);
    border: 1px solid var(--border);
    border-left: 3px solid var(--accent);
    color: var(--text);
    font-family: 'IBM Plex Mono', monospace;
    font-size: 11px;
    padding: 8px 12px;
    border-radius: 4px;
    pointer-events: none;
    opacity: 0;
    transition: opacity 0.1s;
    z-index: 1000;
    max-width: 380px;
    line-height: 1.7;
  }
  #tooltip.visible { opacity: 1; }
  #tooltip .tt-title { font-weight: 600; color: var(--text-bright); font-size: 12px; margin-bottom: 4px; }
  #tooltip .tt-row { display: flex; gap: 12px; }
  #tooltip .tt-label { color: var(--text-dim); min-width: 70px; }

  #legend {
    background: var(--surface);
    border-top: 1px solid var(--border);
    padding: 6px 20px;
    display: flex;
    flex-wrap: wrap;
    gap: 14px;
    align-items: center;
    flex-shrink: 0;
    font-size: 11px;
  }
  .legend-item { display: flex; align-items: center; gap: 5px; cursor: pointer; opacity: 1; transition: opacity 0.2s; }
  .legend-item.hidden { opacity: 0.3; }
  .legend-swatch { width: 14px; height: 10px; border-radius: 2px; flex-shrink: 0; }
  .legend-label { color: var(--text); font-family: 'IBM Plex Mono', monospace; }

  .scale-bar text { fill: var(--text-dim); font-family: 'IBM Plex Mono', monospace; font-size: 10px; }
  .scale-bar line { stroke: var(--text-dim); }
  .track-label { fill: var(--text-dim); font-family: 'IBM Plex Mono', monospace; font-size: 10px; text-transform: uppercase; letter-spacing: 0.06em; }
  .ruler text { fill: var(--text-dim); font-family: 'IBM Plex Mono', monospace; font-size: 10px; }
  .ruler line { stroke: var(--border); }
  .ruler .major { stroke: var(--text-dim); }

  #statusbar {
    position: absolute;
    bottom: 8px;
    right: 12px;
    font-family: 'IBM Plex Mono', monospace;
    font-size: 10px;
    color: var(--text-dim);
    pointer-events: none;
  }

  .no-data-msg {
    fill: var(--text-dim);
    font-family: 'IBM Plex Sans', sans-serif;
    font-size: 14px;
  }
</style>
</head>
<body>

<header>
  <div class="brand">⬡ TE VIEWER</div>
  <div class="controls">
    <label>Chromosome</label>
    <select id="chrom-select"></select>
  </div>
  <div class="sep"></div>
  <div class="controls">
    <label>Start</label>
    <input type="number" id="coord-start" placeholder="0">
    <label>End</label>
    <input type="number" id="coord-end" placeholder="1000000">
    <button class="primary" id="btn-go">Go</button>
  </div>
  <div class="sep"></div>
  <div class="controls">
    <button id="btn-zoom-in">＋ Zoom</button>
    <button id="btn-zoom-out">－ Zoom</button>
    <button id="btn-reset">↺ Reset</button>
  </div>
  <div class="sep"></div>
  <div class="controls">
    <label>Search gene</label>
    <input type="text" id="gene-search" placeholder="gene name…" style="width:130px">
    <button id="btn-gene-go">Find</button>
  </div>
</header>

<div id="minimap-container">
  <svg id="minimap"></svg>
</div>

<div id="viewer-wrap">
  <svg id="viewer-svg"></svg>
  <div id="statusbar"></div>
</div>

<div id="legend"></div>
<div id="tooltip"></div>

<script>
// ═══════════════════════════════════════════════════════
//  DATA (injected by Python)
// ═══════════════════════════════════════════════════════
const RAW_TES  = __TE_DATA__;
const RAW_GENES = __GENE_DATA__;
const RAW_COVERAGE = __COVERAGE_DATA__;
const RAW_TRANSCRIPTS = __TRANSCRIPT_DATA__;
const BAM_LABELS = __BAM_LABELS__;

// ═══════════════════════════════════════════════════════
//  CONSTANTS & CONFIG
// ═══════════════════════════════════════════════════════
const RULER_H       = 28;
const GENE_H        = 18;
const GENE_ROW_H    = 26;
const TE_ROW_H      = 22;
const TE_BAR_H      = 14;
const LTR_SHADE     = 0.6;
const DOMAIN_H      = 8;
const PAD_LEFT      = 70;
const PAD_RIGHT     = 20;
const PAD_TOP       = 14;

// BAM track dimensions
const COV_TRACK_H   = 50;   // height of each coverage barplot
const COV_GAP       = 4;
const TX_ROW_H      = 16;   // height of each transcript row
const TX_BAR_H      = 10;
const TX_TRACK_PAD  = 6;

// Colors for BAM samples
const BAM_COLORS = [
  "#e06c75", "#61afef", "#98c379", "#d19a66", "#c678dd",
  "#56b6c2", "#e5c07b", "#ff6b81"
];

const FAMILY_PALETTE = [
  "#4a9eff","#ff6b6b","#51cf66","#ffd43b","#cc5de8",
  "#ff922b","#20c997","#f06595","#74c0fc","#a9e34b",
  "#e599f7","#ffec99","#63e6be","#ffb3c6","#d0bfff",
  "#a5d8ff","#b2f2bb","#ffcccc","#ffe8cc","#c5f6fa"
];

// ═══════════════════════════════════════════════════════
//  COLOUR SYSTEM
// ═══════════════════════════════════════════════════════
const allFamilies = [...new Set(RAW_TES.map(d => d.family))].sort();
const familyColor = new Map();
allFamilies.forEach((f, i) => {
  familyColor.set(f, FAMILY_PALETTE[i % FAMILY_PALETTE.length]);
});

function ltrColor(family) {
  const c = d3.color(familyColor.get(family) || "#888");
  if (c) { c.opacity = 1; return c.darker(0.9).formatHex(); }
  return "#555";
}
function internalColor(family) { return familyColor.get(family) || "#888"; }

const GENE_BODY_COLOR  = "#3a6a9a";
const GENE_EXON_COLOR  = "#7bb8f0";
const GENE_CDS_COLOR   = "#9acfff";
const GENE_UTR_COLOR   = "#3a6a9a";
const GENE_ARROW_COLOR = "#5b9bd5";

// ═══════════════════════════════════════════════════════
//  STATE
// ═══════════════════════════════════════════════════════
const chroms = [...new Set([
  ...RAW_TES.map(d => d.chrom),
  ...RAW_GENES.map(d => d.chrom)
])].sort();

let currentChrom = chroms[0] || '';
let viewStart = 0;
let viewEnd = 1;
let hiddenFamilies = new Set();
let hiddenSamples = new Set();

// ═══════════════════════════════════════════════════════
//  CHROM DATA
// ═══════════════════════════════════════════════════════
function getChromExtent(chrom) {
  const tes = RAW_TES.filter(d => d.chrom === chrom);
  const genes = RAW_GENES.filter(d => d.chrom === chrom);
  const all = [...tes.map(d => [d.start, d.end]), ...genes.map(d => [d.start, d.end])];
  if (!all.length) return [0, 1000000];
  return [0, Math.max(...all.map(d => d[1])) + 1000];
}

function setChrom(chrom) {
  currentChrom = chrom;
  const ext = getChromExtent(chrom);
  viewStart = ext[0];
  viewEnd = ext[1];
  renderAll();
}

// ═══════════════════════════════════════════════════════
//  LAYOUT CALCULATOR
// ═══════════════════════════════════════════════════════
function computeLayout(chrom, vStart, vEnd) {
  const tes = RAW_TES.filter(d =>
    d.chrom === chrom && d.end > vStart && d.start < vEnd && !hiddenFamilies.has(d.family)
  );
  const genes = RAW_GENES.filter(d =>
    d.chrom === chrom && d.end > vStart && d.start < vEnd
  );

  const maxLevel = tes.length ? Math.max(...tes.map(d => d.level)) : 0;
  const teRowCount = maxLevel + 1;
  const geneRows = packRows(genes);

  // Compute visible BAM labels (not hidden)
  const visibleLabels = BAM_LABELS.filter(l => !hiddenSamples.has(l));

  // Coverage bins per visible label
  const covByLabel = {};
  for (const label of visibleLabels) {
    const chromCov = (RAW_COVERAGE[label] || {})[chrom] || [];
    covByLabel[label] = chromCov.filter(b => b.e > vStart && b.s < vEnd);
  }

  // Transcripts per visible label
  const txByLabel = {};
  for (const label of visibleLabels) {
    const allTx = (RAW_TRANSCRIPTS[label] || []).filter(t =>
      t.c === chrom && t.e > vStart && t.s < vEnd
    );
    txByLabel[label] = packRows(allTx.map(t => ({start: t.s, end: t.e, ...t})));
  }

  return { tes, genes, geneRows, teRowCount, maxLevel, visibleLabels, covByLabel, txByLabel };
}

function packRows(features) {
  const rows = [];
  const sorted = [...features].sort((a, b) => a.start - b.start);
  for (const f of sorted) {
    let placed = false;
    for (const row of rows) {
      if (row[row.length - 1].end <= f.start) {
        row.push(f);
        placed = true;
        break;
      }
    }
    if (!placed) rows.push([f]);
  }
  return rows;
}

// ═══════════════════════════════════════════════════════
//  SVG DIMENSIONS
// ═══════════════════════════════════════════════════════
function svgDims() {
  const wrap = document.getElementById('viewer-wrap');
  return { W: wrap.clientWidth, H: wrap.clientHeight };
}
function trackWidth(W) { return W - PAD_LEFT - PAD_RIGHT; }

// ═══════════════════════════════════════════════════════
//  MAIN RENDER
// ═══════════════════════════════════════════════════════
function renderAll() {
  const { W, H } = svgDims();
  const TW = trackWidth(W);

  const layout = computeLayout(currentChrom, viewStart, viewEnd);
  const { tes, genes, geneRows, teRowCount, visibleLabels, covByLabel, txByLabel } = layout;

  const xScale = d3.scaleLinear()
    .domain([viewStart, viewEnd])
    .range([PAD_LEFT, PAD_LEFT + TW]);

  // compute total height
  const geneBlockH = geneRows.length * GENE_ROW_H;
  const teBlockH   = teRowCount * TE_ROW_H;

  // BAM tracks height
  let bamBlockH = 0;
  for (const label of visibleLabels) {
    bamBlockH += COV_TRACK_H + COV_GAP;  // coverage
    const txRows = txByLabel[label] || [];
    const txH = Math.max(0, txRows.length) * TX_ROW_H + TX_TRACK_PAD;
    bamBlockH += txH + 8;  // transcripts + gap
  }

  const totalContentH = PAD_TOP + RULER_H + geneBlockH + 10 + teBlockH + 20 + bamBlockH + 40;
  const svgH = Math.max(H, totalContentH);

  const svg = d3.select('#viewer-svg');
  svg.selectAll('*').remove();
  svg.attr('viewBox', `0 0 ${W} ${svgH}`)
     .attr('width', W).attr('height', svgH);

  let yOffset = PAD_TOP;

  // ── Ruler ──
  drawRuler(svg, xScale, yOffset, W, TW);
  yOffset += RULER_H;

  // ── Gene track ──
  if (geneRows.length > 0) {
    svg.append('text').attr('class','track-label')
      .attr('x', 4).attr('y', yOffset + 12).text('GENES');
    drawGenes(svg, geneRows, xScale, yOffset);
    yOffset += geneBlockH + 10;
  }

  // ── TE track ──
  svg.append('text').attr('class','track-label')
    .attr('x', 4).attr('y', yOffset + 12).text('TEs');

  for (let lvl = 0; lvl <= layout.maxLevel; lvl++) {
    const rowTEs = tes.filter(d => d.level === lvl);
    const rowY = yOffset + lvl * TE_ROW_H;
    drawTERow(svg, rowTEs, xScale, rowY, lvl);
  }
  drawNestingConnectors(svg, tes, xScale, yOffset);
  yOffset += teBlockH + 20;

  // ── BAM tracks (coverage + transcripts per sample) ──
  for (let si = 0; si < visibleLabels.length; si++) {
    const label = visibleLabels[si];
    const color = BAM_COLORS[si % BAM_COLORS.length];

    // Separator line
    svg.append('line')
      .attr('x1', PAD_LEFT).attr('x2', PAD_LEFT + TW)
      .attr('y1', yOffset - 2).attr('y2', yOffset - 2)
      .attr('stroke', 'var(--border)').attr('stroke-width', 0.5);

    // ── Coverage barplot ──
    svg.append('text').attr('class','track-label')
      .attr('x', 4).attr('y', yOffset + 10)
      .text(label.toUpperCase());
    svg.append('text')
      .attr('x', 4).attr('y', yOffset + 20)
      .attr('fill','var(--text-dim)')
      .attr('font-family','IBM Plex Mono,monospace')
      .attr('font-size','8px')
      .text('coverage');

    drawCoverageTrack(svg, covByLabel[label] || [], xScale, yOffset, TW, color);
    yOffset += COV_TRACK_H + COV_GAP;

    // ── Transcript track ──
    svg.append('text')
      .attr('x', 4).attr('y', yOffset + 9)
      .attr('fill','var(--text-dim)')
      .attr('font-family','IBM Plex Mono,monospace')
      .attr('font-size','8px')
      .text('transcripts');

    const txRows = txByLabel[label] || [];
    drawTranscriptTrack(svg, txRows, xScale, yOffset, color);
    const txH = Math.max(0, txRows.length) * TX_ROW_H + TX_TRACK_PAD;
    yOffset += txH + 8;
  }

  // ── Interaction ──
  setupInteraction(svg, xScale, W, TW, svgH);

  document.getElementById('statusbar').textContent =
    `${currentChrom}  ${fmt(viewStart)}–${fmt(viewEnd)}  (${fmt(viewEnd - viewStart)} bp)`;

  renderMinimap(tes, genes);
}

// ═══════════════════════════════════════════════════════
//  RULER
// ═══════════════════════════════════════════════════════
function drawRuler(svg, xScale, y, W, TW) {
  const g = svg.append('g').attr('class','ruler');
  const span = viewEnd - viewStart;
  const step = niceStep(span, Math.floor(TW / 80));

  g.append('line')
    .attr('x1', PAD_LEFT).attr('x2', PAD_LEFT + TW)
    .attr('y1', y + RULER_H - 1).attr('y2', y + RULER_H - 1)
    .attr('stroke','var(--border)').attr('stroke-width',1);

  const ticks = d3.range(
    Math.ceil(viewStart / step) * step, viewEnd, step
  );
  ticks.forEach(t => {
    const x = xScale(t);
    if (x < PAD_LEFT || x > PAD_LEFT + TW) return;
    g.append('line').attr('class','major')
      .attr('x1',x).attr('x2',x)
      .attr('y1',y+RULER_H-8).attr('y2',y+RULER_H)
      .attr('stroke-width',1);
    g.append('text')
      .attr('x',x).attr('y',y+RULER_H-11)
      .attr('text-anchor','middle')
      .text(fmt(t));
  });
}

function niceStep(span, targetTicks) {
  const raw = span / targetTicks;
  const mag = Math.pow(10, Math.floor(Math.log10(raw)));
  const norm = raw / mag;
  if (norm < 1.5) return mag;
  if (norm < 3.5) return 2 * mag;
  if (norm < 7.5) return 5 * mag;
  return 10 * mag;
}

// ═══════════════════════════════════════════════════════
//  GENE DRAWING
// ═══════════════════════════════════════════════════════
function drawGenes(svg, geneRows, xScale, yOffset) {
  const g = svg.append('g').attr('class','gene-track');

  geneRows.forEach((row, ri) => {
    const ry = yOffset + ri * GENE_ROW_H + (GENE_ROW_H - GENE_H) / 2;
    row.forEach(gene => {
      const gx = xScale(gene.start);
      const gw = Math.max(2, xScale(gene.end) - xScale(gene.start));
      const gg = g.append('g').attr('class','gene-feat');

      gg.append('line')
        .attr('x1', gx).attr('x2', gx + gw)
        .attr('y1', ry + GENE_H / 2).attr('y2', ry + GENE_H / 2)
        .attr('stroke', GENE_BODY_COLOR).attr('stroke-width', 1.5);

      const arrowSpacing = 40;
      const arrowN = Math.floor(gw / arrowSpacing);
      for (let i = 0; i <= arrowN; i++) {
        const ax = gx + (i / Math.max(arrowN, 1)) * gw;
        const dir = gene.strand === '-' ? -1 : 1;
        const aw = 5;
        gg.append('path')
          .attr('d', `M${ax},${ry + GENE_H/2 - 3} L${ax + dir*aw},${ry + GENE_H/2} L${ax},${ry + GENE_H/2 + 3}`)
          .attr('fill','none').attr('stroke', GENE_ARROW_COLOR)
          .attr('stroke-width',1).attr('opacity',0.5);
      }

      gene.exons.forEach(ex => {
        const ex_x = xScale(ex.start);
        const ex_w = Math.max(1, xScale(ex.end) - xScale(ex.start));
        let color = GENE_EXON_COLOR;
        let h = GENE_H - 4;
        let dy = 2;
        if (ex.type === 'CDS') { color = GENE_CDS_COLOR; h = GENE_H; dy = 0; }
        if (ex.type.includes('UTR') || ex.type.includes('utr')) { color = GENE_UTR_COLOR; h = GENE_H - 6; dy = 3; }
        gg.append('rect')
          .attr('x', ex_x).attr('y', ry + dy)
          .attr('width', ex_w).attr('height', h)
          .attr('fill', color).attr('rx', 1);
      });

      if (gene.exons.length === 0) {
        gg.append('rect')
          .attr('x', gx).attr('y', ry)
          .attr('width', gw).attr('height', GENE_H)
          .attr('fill', GENE_EXON_COLOR).attr('rx', 2);
      }

      if (gw > 30) {
        gg.append('text')
          .attr('x', gx + gw / 2).attr('y', ry - 2)
          .attr('text-anchor','middle')
          .attr('fill','var(--text-dim)')
          .attr('font-family','IBM Plex Mono, monospace')
          .attr('font-size','9px')
          .text(gene.name.length > 16 ? gene.name.slice(0,14)+'…' : gene.name);
      }

      gg.selectAll('rect,line').on('mousemove', (event) => {
        showTooltip(event, buildGeneTooltip(gene));
      }).on('mouseleave', hideTooltip);
    });
  });
}

// ═══════════════════════════════════════════════════════
//  TE DRAWING
// ═══════════════════════════════════════════════════════
function drawTERow(svg, rowTEs, xScale, rowY, level) {
  const g = svg.append('g').attr('class', `te-row te-level-${level}`);
  const barY = rowY + (TE_ROW_H - TE_BAR_H) / 2;

  svg.append('text')
    .attr('x', PAD_LEFT - 4)
    .attr('y', rowY + TE_ROW_H / 2 + 4)
    .attr('text-anchor','end')
    .attr('fill','var(--text-dim)')
    .attr('font-family','IBM Plex Mono,monospace')
    .attr('font-size','9px')
    .text(`L${level}`);

  rowTEs.forEach(te => {
    const isLTR = te.te_class === 'LTR';
    const fragG = g.append('g').attr('class','te-element');

    te.visible_frags.forEach(([fs, fe]) => {
      const fx = xScale(fs);
      const fw = Math.max(1, xScale(fe) - xScale(fs));

      if (isLTR && te.ltr_len > 0) {
        const ltrLeft  = { s: te.start, e: te.start + te.ltr_len };
        const ltrRight = { s: te.end - te.ltr_len, e: te.end };
        drawLTRFragment(fragG, fs, fe, te, ltrLeft, ltrRight, xScale, barY);
      } else {
        fragG.append('rect')
          .attr('x', fx).attr('y', barY)
          .attr('width', fw).attr('height', TE_BAR_H)
          .attr('fill', internalColor(te.family))
          .attr('rx', 2).attr('opacity', 0.88);
      }
    });

    if (te.domains && te.domains.length > 0) {
      te.domains.forEach(dom => {
        if (dom.start >= te.start && dom.end <= te.end) {
          const dx = xScale(dom.start);
          const dw = Math.max(2, xScale(dom.end) - xScale(dom.start));
          const domY = barY + (TE_BAR_H - DOMAIN_H) / 2;
          fragG.append('rect')
            .attr('x', dx).attr('y', domY)
            .attr('width', dw).attr('height', DOMAIN_H)
            .attr('fill', 'rgba(255,255,255,0.25)')
            .attr('rx', 1);
          if (dw > 14) {
            fragG.append('text')
              .attr('x', dx + dw/2).attr('y', domY + DOMAIN_H - 1)
              .attr('text-anchor','middle')
              .attr('fill','rgba(255,255,255,0.85)')
              .attr('font-family','IBM Plex Mono,monospace')
              .attr('font-size','7px')
              .text(dom.name);
          }
        }
      });
    }

    const teX = xScale(te.start);
    const teW = xScale(te.end) - xScale(te.start);
    if (te.k2p_dist > 0 && teW > 40 && te.visible_frags.length > 0) {
      const firstFrag = te.visible_frags[0];
      const lx = xScale(firstFrag[0]) + 2;
      fragG.append('text')
        .attr('x', lx).attr('y', barY - 2)
        .attr('fill','var(--text-dim)')
        .attr('font-family','IBM Plex Mono,monospace')
        .attr('font-size','8px')
        .text(`K2P:${te.k2p_dist.toFixed(3)}`);
    }

    if (te.visible_frags.length > 0) {
      const allFragsX1 = xScale(Math.min(...te.visible_frags.map(f=>f[0])));
      const allFragsX2 = xScale(Math.max(...te.visible_frags.map(f=>f[1])));
      fragG.append('rect')
        .attr('x', allFragsX1).attr('y', barY - 3)
        .attr('width', Math.max(1, allFragsX2 - allFragsX1))
        .attr('height', TE_BAR_H + 6)
        .attr('fill','transparent')
        .on('mousemove', (event) => showTooltip(event, buildTETooltip(te)))
        .on('mouseleave', hideTooltip);
    }
  });
}

function drawLTRFragment(g, fs, fe, te, ltrLeft, ltrRight, xScale, barY) {
  const segments = [];
  const ll_s = Math.max(fs, ltrLeft.s);
  const ll_e = Math.min(fe, ltrLeft.e);
  if (ll_e > ll_s) segments.push({ s: ll_s, e: ll_e, type: 'ltr' });
  const int_s = Math.max(fs, ltrLeft.e);
  const int_e = Math.min(fe, ltrRight.s);
  if (int_e > int_s) segments.push({ s: int_s, e: int_e, type: 'internal' });
  const rl_s = Math.max(fs, ltrRight.s);
  const rl_e = Math.min(fe, ltrRight.e);
  if (rl_e > rl_s) segments.push({ s: rl_s, e: rl_e, type: 'ltr' });

  if (segments.length === 0) {
    g.append('rect')
      .attr('x', xScale(fs)).attr('y', barY)
      .attr('width', Math.max(1, xScale(fe) - xScale(fs)))
      .attr('height', TE_BAR_H)
      .attr('fill', internalColor(te.family))
      .attr('rx', 2).attr('opacity', 0.88);
    return;
  }

  segments.forEach((seg, idx) => {
    const sx = xScale(seg.s);
    const sw = Math.max(1, xScale(seg.e) - xScale(seg.s));
    const color = seg.type === 'ltr' ? ltrColor(te.family) : internalColor(te.family);
    const h = seg.type === 'ltr' ? TE_BAR_H : TE_BAR_H - 2;
    const dy = seg.type === 'ltr' ? 0 : 1;
    g.append('rect')
      .attr('x', sx).attr('y', barY + dy)
      .attr('width', sw).attr('height', h)
      .attr('fill', color)
      .attr('rx', (idx === 0 || idx === segments.length - 1) ? 2 : 0)
      .attr('opacity', 0.9);
  });
}

// ═══════════════════════════════════════════════════════
//  NESTING CONNECTORS
// ═══════════════════════════════════════════════════════
function drawNestingConnectors(svg, tes, xScale, teTrackY) {
  const g = svg.append('g').attr('class','connectors').attr('opacity',0.25);
  tes.forEach(parent => {
    if (parent.level === 0) return;
    const childLevel = parent.level - 1;
    const children = tes.filter(c =>
      c.level === childLevel &&
      c.chrom === parent.chrom &&
      c.start >= parent.start && c.end <= parent.end
    );
    children.forEach(child => {
      const cy_child  = teTrackY + childLevel  * TE_ROW_H + (TE_ROW_H - TE_BAR_H) / 2;
      const cy_parent = teTrackY + parent.level * TE_ROW_H + (TE_ROW_H - TE_BAR_H) / 2;
      const cx = xScale((child.start + child.end) / 2);
      g.append('line')
        .attr('x1', cx).attr('y1', cy_child + TE_BAR_H)
        .attr('x2', cx).attr('y2', cy_parent)
        .attr('stroke','var(--text-dim)').attr('stroke-width',0.5)
        .attr('stroke-dasharray','2,2');
    });
  });
}

// ═══════════════════════════════════════════════════════
//  COVERAGE BARPLOT
// ═══════════════════════════════════════════════════════
function drawCoverageTrack(svg, bins, xScale, yOffset, TW, color) {
  const g = svg.append('g').attr('class','cov-track');

  if (!bins.length) {
    g.append('text')
      .attr('x', PAD_LEFT + TW / 2).attr('y', yOffset + COV_TRACK_H / 2)
      .attr('text-anchor','middle')
      .attr('fill','var(--text-dim)')
      .attr('font-size','10px')
      .text('no coverage in view');
    return;
  }

  const maxDepth = Math.max(...bins.map(b => b.d), 1);
  const yScale = d3.scaleLinear()
    .domain([0, maxDepth])
    .range([yOffset + COV_TRACK_H, yOffset + 4]);

  // Draw axis ticks (max value)
  g.append('text')
    .attr('x', PAD_LEFT - 4).attr('y', yOffset + 8)
    .attr('text-anchor','end')
    .attr('fill','var(--text-dim)')
    .attr('font-family','IBM Plex Mono,monospace')
    .attr('font-size','8px')
    .text(maxDepth >= 1000 ? (maxDepth/1000).toFixed(1)+'k' : Math.round(maxDepth));

  g.append('text')
    .attr('x', PAD_LEFT - 4).attr('y', yOffset + COV_TRACK_H)
    .attr('text-anchor','end')
    .attr('fill','var(--text-dim)')
    .attr('font-family','IBM Plex Mono,monospace')
    .attr('font-size','8px')
    .text('0');

  // Baseline
  g.append('line')
    .attr('x1', PAD_LEFT).attr('x2', PAD_LEFT + TW)
    .attr('y1', yOffset + COV_TRACK_H).attr('y2', yOffset + COV_TRACK_H)
    .attr('stroke','var(--border)').attr('stroke-width', 0.5);

  // Draw bars
  bins.forEach(bin => {
    const bx = xScale(bin.s);
    const bw = Math.max(1, xScale(bin.e) - xScale(bin.s));
    const by = yScale(bin.d);
    const bh = yOffset + COV_TRACK_H - by;
    g.append('rect')
      .attr('x', bx).attr('y', by)
      .attr('width', bw).attr('height', Math.max(0.5, bh))
      .attr('fill', color).attr('opacity', 0.7);
  });

  // Tooltip on the coverage region
  g.append('rect')
    .attr('x', PAD_LEFT).attr('y', yOffset)
    .attr('width', TW).attr('height', COV_TRACK_H)
    .attr('fill', 'transparent')
    .on('mousemove', (event) => {
      const [mx] = d3.pointer(event);
      const genomicPos = xScale.invert(mx);
      // Find closest bin
      let closest = null;
      let minDist = Infinity;
      for (const bin of bins) {
        const mid = (bin.s + bin.e) / 2;
        const dist = Math.abs(mid - genomicPos);
        if (dist < minDist) { minDist = dist; closest = bin; }
      }
      if (closest) {
        showTooltip(event,
          `<div class="tt-title">Coverage</div>
           <div class="tt-row"><span class="tt-label">Position</span><span>${fmt(closest.s)}-${fmt(closest.e)}</span></div>
           <div class="tt-row"><span class="tt-label">Depth</span><span>${closest.d.toFixed(1)}×</span></div>`
        );
      }
    })
    .on('mouseleave', hideTooltip);
}

// ═══════════════════════════════════════════════════════
//  TRANSCRIPT TRACK
// ═══════════════════════════════════════════════════════
function drawTranscriptTrack(svg, txRows, xScale, yOffset, color) {
  const g = svg.append('g').attr('class','tx-track');

  txRows.forEach((row, ri) => {
    const ry = yOffset + ri * TX_ROW_H + (TX_ROW_H - TX_BAR_H) / 2;
    row.forEach(tx => {
      const txG = g.append('g');

      // Intron line spanning full extent
      const tx_x1 = xScale(tx.s);
      const tx_x2 = xScale(tx.e);
      txG.append('line')
        .attr('x1', tx_x1).attr('x2', tx_x2)
        .attr('y1', ry + TX_BAR_H / 2).attr('y2', ry + TX_BAR_H / 2)
        .attr('stroke', color).attr('stroke-width', 1).attr('opacity', 0.5);

      // Strand arrows
      const dir = tx.st === '-' ? -1 : 1;
      const gw = tx_x2 - tx_x1;
      const arrowN = Math.max(1, Math.floor(gw / 50));
      for (let i = 0; i <= arrowN; i++) {
        const ax = tx_x1 + (i / Math.max(arrowN, 1)) * gw;
        const aw = 4;
        txG.append('path')
          .attr('d', `M${ax},${ry + TX_BAR_H/2 - 2} L${ax + dir*aw},${ry + TX_BAR_H/2} L${ax},${ry + TX_BAR_H/2 + 2}`)
          .attr('fill','none').attr('stroke', color)
          .attr('stroke-width', 0.8).attr('opacity', 0.4);
      }

      // Exon blocks
      (tx.b || []).forEach(block => {
        const bs = block[0], be = block[1];
        const bx = xScale(bs);
        const bw = Math.max(1, xScale(be) - xScale(bs));
        txG.append('rect')
          .attr('x', bx).attr('y', ry)
          .attr('width', bw).attr('height', TX_BAR_H)
          .attr('fill', color).attr('opacity', 0.75).attr('rx', 1);
      });

      // Depth label if space
      if (gw > 30 && tx.d > 1) {
        txG.append('text')
          .attr('x', tx_x1 + gw / 2).attr('y', ry - 1)
          .attr('text-anchor','middle')
          .attr('fill','var(--text-dim)')
          .attr('font-family','IBM Plex Mono,monospace')
          .attr('font-size','7px')
          .text(`${tx.d}×`);
      }

      // Tooltip hit area
      txG.append('rect')
        .attr('x', tx_x1).attr('y', ry - 2)
        .attr('width', Math.max(1, tx_x2 - tx_x1))
        .attr('height', TX_BAR_H + 4)
        .attr('fill','transparent')
        .on('mousemove', (event) => {
          const nBlocks = (tx.b || []).length;
          const totalExon = (tx.b || []).reduce((a, b) => a + (b[1] - b[0]), 0);
          showTooltip(event,
            `<div class="tt-title">Transcript cluster</div>
             <div class="tt-row"><span class="tt-label">Span</span><span>${fmt(tx.s)}-${fmt(tx.e)}</span></div>
             <div class="tt-row"><span class="tt-label">Strand</span><span>${tx.st}</span></div>
             <div class="tt-row"><span class="tt-label">Exon blocks</span><span>${nBlocks}</span></div>
             <div class="tt-row"><span class="tt-label">Exon bp</span><span>${fmt(totalExon)}</span></div>
             <div class="tt-row"><span class="tt-label">Read depth</span><span>${tx.d} reads</span></div>`
          );
        })
        .on('mouseleave', hideTooltip);
    });
  });
}

// ═══════════════════════════════════════════════════════
//  MINIMAP
// ═══════════════════════════════════════════════════════
function renderMinimap(tes, genes) {
  const container = document.getElementById('minimap-container');
  const W = container.clientWidth - 40;
  const H = 28;
  const ext = getChromExtent(currentChrom);
  const chromLen = ext[1];

  const svg = d3.select('#minimap');
  svg.attr('width', W + 40).attr('height', H);
  svg.selectAll('*').remove();

  const mx = d3.scaleLinear().domain([0, chromLen]).range([PAD_LEFT, PAD_LEFT + W - PAD_LEFT]);

  svg.append('rect').attr('x',PAD_LEFT).attr('y',4).attr('width',W-PAD_LEFT).attr('height',H-8)
     .attr('fill','var(--surface2)').attr('rx',2);

  const allChromTEs = RAW_TES.filter(d => d.chrom === currentChrom && !hiddenFamilies.has(d.family));
  allChromTEs.forEach(te => {
    const tx = mx(te.start);
    const tw = Math.max(1, mx(te.end) - mx(te.start));
    svg.append('rect').attr('x', tx).attr('y', 4 + H*0.25)
       .attr('width', tw).attr('height', H * 0.5)
       .attr('fill', internalColor(te.family)).attr('opacity', 0.4);
  });

  const allChromGenes = RAW_GENES.filter(d => d.chrom === currentChrom);
  allChromGenes.forEach(g => {
    const gx = mx(g.start);
    svg.append('line').attr('x1',gx).attr('x2',gx).attr('y1',4).attr('y2',H-4)
       .attr('stroke', GENE_ARROW_COLOR).attr('stroke-width',1).attr('opacity',0.6);
  });

  const vx = mx(viewStart);
  const vw = Math.max(4, mx(viewEnd) - mx(viewStart));
  const vp = svg.append('rect')
    .attr('x', vx).attr('y', 2)
    .attr('width', vw).attr('height', H - 4)
    .attr('fill','rgba(74,158,255,0.12)')
    .attr('stroke','var(--accent)').attr('stroke-width',1)
    .attr('rx',2).attr('cursor','ew-resize');

  const drag = d3.drag().on('drag', (event) => {
    const clickX = Math.max(PAD_LEFT, Math.min(PAD_LEFT + W - PAD_LEFT, event.x));
    const genomicPos = mx.invert(clickX);
    const span = viewEnd - viewStart;
    viewStart = Math.max(0, genomicPos - span / 2);
    viewEnd = viewStart + span;
    renderAll();
  });
  vp.call(drag);

  svg.on('click', (event) => {
    const [px] = d3.pointer(event);
    const genomicPos = mx.invert(px);
    const span = viewEnd - viewStart;
    viewStart = Math.max(0, genomicPos - span / 2);
    viewEnd = viewStart + span;
    renderAll();
  });
}

// ═══════════════════════════════════════════════════════
//  PAN / ZOOM INTERACTION
// ═══════════════════════════════════════════════════════
function setupInteraction(svg, xScale, W, TW, svgH) {
  let dragStartX   = null;
  let dragStartView = null;
  let isDragging   = false;

  svg.on('mousedown', (event) => {
    if (event.button !== 0) return;
    dragStartX    = event.clientX;
    dragStartView = [viewStart, viewEnd];
    isDragging    = false;
    event.preventDefault();
  });

  window.addEventListener('mousemove', (event) => {
    if (dragStartX === null) return;
    const dx = event.clientX - dragStartX;
    if (!isDragging && Math.abs(dx) > 2) {
      isDragging = true;
      document.getElementById('viewer-svg').style.cursor = 'grabbing';
    }
    if (!isDragging) return;
    const span = dragStartView[1] - dragStartView[0];
    const genomicDx = -dx / TW * span;
    viewStart = dragStartView[0] + genomicDx;
    viewEnd   = dragStartView[1] + genomicDx;
    renderAll();
  });

  window.addEventListener('mouseup', () => {
    dragStartX = null;
    isDragging = false;
    document.getElementById('viewer-svg').style.cursor = 'grab';
  });

  svg.on('wheel', (event) => {
    event.preventDefault();
    const factor = event.deltaY > 0 ? 1.25 : 0.8;
    const span = viewEnd - viewStart;
    const mx = event.offsetX;
    const genomicMx = xScale.invert(mx);
    const ratio = (genomicMx - viewStart) / span;
    const newSpan = Math.max(100, span * factor);
    viewStart = genomicMx - ratio * newSpan;
    viewEnd   = viewStart + newSpan;
    renderAll();
  }, { passive: false });
}

// ═══════════════════════════════════════════════════════
//  TOOLTIP
// ═══════════════════════════════════════════════════════
function buildTETooltip(te) {
  const domains = te.domains.map(d => `${d.name} [${fmt(d.start)}-${fmt(d.end)}] e=${d.evalue.toExponential(1)}`).join('<br>');
  return `<div class="tt-title">${te.family} <span style="opacity:0.5;font-size:10px">${te.superfamily}</span></div>
<div class="tt-row"><span class="tt-label">Coords</span><span>${te.chrom}:${fmt(te.start)}-${fmt(te.end)}</span></div>
<div class="tt-row"><span class="tt-label">Length</span><span>${fmt(te.end - te.start)} bp</span></div>
${te.ltr_len ? `<div class="tt-row"><span class="tt-label">LTR len</span><span>${te.ltr_len} bp</span></div>` : ''}
${te.has_divergence ? `<div class="tt-row"><span class="tt-label">K2P dist</span><span>${te.k2p_dist.toFixed(4)}</span></div>` : ''}
${te.has_divergence ? `<div class="tt-row"><span class="tt-label">K2P time</span><span>${(te.k2p_time/1e6).toFixed(2)} Mya</span></div>` : ''}
<div class="tt-row"><span class="tt-label">Level</span><span>${te.level}</span></div>
${domains ? `<div style="margin-top:4px;color:var(--text-dim);font-size:10px">Domains:<br>${domains}</div>` : ''}`;
}

function buildGeneTooltip(gene) {
  const targets = (gene.targets && gene.targets.length > 1)
    ? `<div style="margin-top:4px;color:var(--text-dim);font-size:10px">All targets (${gene.targets.length}):<br>${gene.targets.join('<br>')}</div>`
    : '';
  return `<div class="tt-title">${gene.name}</div>
<div class="tt-row"><span class="tt-label">Coords</span><span>${gene.chrom}:${fmt(gene.start)}-${fmt(gene.end)}</span></div>
<div class="tt-row"><span class="tt-label">Strand</span><span>${gene.strand}</span></div>
<div class="tt-row"><span class="tt-label">Length</span><span>${fmt(gene.end - gene.start)} bp</span></div>
<div class="tt-row"><span class="tt-label">CDS blocks</span><span>${gene.exons.length}</span></div>
${gene.targets && gene.targets.length > 1 ? `<div class="tt-row"><span class="tt-label">Merged</span><span>${gene.targets.length} overlapping alignments</span></div>` : ''}
${targets}`;
}

function showTooltip(event, html) {
  const tt = document.getElementById('tooltip');
  tt.innerHTML = html;
  tt.classList.add('visible');
  positionTooltip(event);
}

function hideTooltip() {
  document.getElementById('tooltip').classList.remove('visible');
}

function positionTooltip(event) {
  const tt = document.getElementById('tooltip');
  const margin = 12;
  let x = event.clientX + margin;
  let y = event.clientY + margin;
  if (x + 350 > window.innerWidth)  x = event.clientX - 350 - margin;
  if (y + 200 > window.innerHeight) y = event.clientY - 200 - margin;
  tt.style.left = x + 'px';
  tt.style.top  = y + 'px';
}

// ═══════════════════════════════════════════════════════
//  LEGEND
// ═══════════════════════════════════════════════════════
function buildLegend() {
  const container = document.getElementById('legend');
  container.innerHTML = '';

  const geneItem = document.createElement('div');
  geneItem.className = 'legend-item';
  geneItem.innerHTML = `<div class="legend-swatch" style="background:${GENE_EXON_COLOR}"></div><span class="legend-label">Gene</span>`;
  container.appendChild(geneItem);

  allFamilies.forEach(family => {
    const item = document.createElement('div');
    item.className = 'legend-item' + (hiddenFamilies.has(family) ? ' hidden' : '');
    item.innerHTML = `<div class="legend-swatch" style="background:${familyColor.get(family)}"></div><span class="legend-label">${family}</span>`;
    item.addEventListener('click', () => {
      if (hiddenFamilies.has(family)) hiddenFamilies.delete(family);
      else hiddenFamilies.add(family);
      item.classList.toggle('hidden');
      renderAll();
    });
    container.appendChild(item);
  });

  const ltrItem = document.createElement('div');
  ltrItem.className = 'legend-item';
  ltrItem.innerHTML = `<div style="display:flex;gap:1px;align-items:center">
    <div style="width:5px;height:10px;background:#555;border-radius:2px 0 0 2px"></div>
    <div style="width:14px;height:8px;background:#666;margin-top:1px"></div>
    <div style="width:5px;height:10px;background:#555;border-radius:0 2px 2px 0"></div>
  </div><span class="legend-label" style="margin-left:5px">LTR | internal | LTR</span>`;
  container.appendChild(ltrItem);

  // BAM sample legend items
  BAM_LABELS.forEach((label, i) => {
    const color = BAM_COLORS[i % BAM_COLORS.length];
    const item = document.createElement('div');
    item.className = 'legend-item' + (hiddenSamples.has(label) ? ' hidden' : '');
    item.innerHTML = `<div class="legend-swatch" style="background:${color}"></div><span class="legend-label">${label} (RNA)</span>`;
    item.addEventListener('click', () => {
      if (hiddenSamples.has(label)) hiddenSamples.delete(label);
      else hiddenSamples.add(label);
      item.classList.toggle('hidden');
      renderAll();
    });
    container.appendChild(item);
  });
}

// ═══════════════════════════════════════════════════════
//  CONTROLS
// ═══════════════════════════════════════════════════════
function initControls() {
  const chromSel = document.getElementById('chrom-select');
  chroms.forEach(c => {
    const opt = document.createElement('option');
    opt.value = c; opt.text = c;
    chromSel.appendChild(opt);
  });
  chromSel.value = currentChrom;
  chromSel.addEventListener('change', () => setChrom(chromSel.value));

  document.getElementById('btn-go').addEventListener('click', () => {
    const s = parseInt(document.getElementById('coord-start').value);
    const e = parseInt(document.getElementById('coord-end').value);
    if (!isNaN(s) && !isNaN(e) && e > s) {
      viewStart = s; viewEnd = e; renderAll();
    }
  });

  document.getElementById('btn-zoom-in').addEventListener('click', () => {
    const center = (viewStart + viewEnd) / 2;
    const span = (viewEnd - viewStart) * 0.5;
    viewStart = center - span / 2;
    viewEnd = center + span / 2;
    renderAll();
  });

  document.getElementById('btn-zoom-out').addEventListener('click', () => {
    const center = (viewStart + viewEnd) / 2;
    const span = (viewEnd - viewStart) * 2;
    viewStart = center - span / 2;
    viewEnd = center + span / 2;
    renderAll();
  });

  document.getElementById('btn-reset').addEventListener('click', () => {
    setChrom(currentChrom);
  });

  document.getElementById('btn-gene-go').addEventListener('click', () => {
    const query = document.getElementById('gene-search').value.trim().toLowerCase();
    const match = RAW_GENES.find(g =>
      g.chrom === currentChrom &&
      (g.name.toLowerCase().includes(query) || g.gene_id.toLowerCase().includes(query))
    );
    if (match) {
      const pad = Math.max(5000, (match.end - match.start) * 2);
      viewStart = match.start - pad;
      viewEnd = match.end + pad;
      renderAll();
    }
  });

  document.getElementById('gene-search').addEventListener('keydown', e => {
    if (e.key === 'Enter') document.getElementById('btn-gene-go').click();
  });
  document.getElementById('coord-start').addEventListener('keydown', e => {
    if (e.key === 'Enter') document.getElementById('btn-go').click();
  });
  document.getElementById('coord-end').addEventListener('keydown', e => {
    if (e.key === 'Enter') document.getElementById('btn-go').click();
  });

  window.addEventListener('resize', renderAll);
}

// ═══════════════════════════════════════════════════════
//  UTILITIES
// ═══════════════════════════════════════════════════════
function fmt(n) {
  if (n >= 1e6) return (n/1e6).toFixed(2) + ' Mb';
  if (n >= 1e3) return (n/1e3).toFixed(1) + ' kb';
  return Math.round(n).toLocaleString();
}

// ═══════════════════════════════════════════════════════
//  INIT
// ═══════════════════════════════════════════════════════
initControls();
buildLegend();
if (currentChrom) setChrom(currentChrom);
</script>
</body>
</html>
"""


# ─────────────────────────────────────────────
#  Main
# ─────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Generate an interactive nested TE + gene genome browser.')
    parser.add_argument('--gff3', required=False, default=None,
                        help='Gene annotation GFF3 file')
    parser.add_argument('--te-fastas', nargs='+', required=False, default=[],
                        metavar='FASTA',
                        help='TE FASTA files, ordered level0 (innermost) first')
    parser.add_argument('--ltr-divergence', nargs='+', required=False, default=[],
                        metavar='DIV',
                        help='LTR divergence files, same order as --te-fastas')
    parser.add_argument('--domains', nargs='+', required=False, default=[],
                        metavar='DOM',
                        help='Protein domain GFF files, same order as --te-fastas')
    parser.add_argument('--bam', nargs='+', required=False, default=[],
                        metavar='BAM',
                        help='BAM files (indexed) for coverage + transcript tracks')
    parser.add_argument('--bam-labels', nargs='+', required=False, default=[],
                        metavar='LABEL',
                        help='Labels for each BAM file (same order as --bam)')
    parser.add_argument('--bin-size', type=int, default=50,
                        help='Coverage bin size in bp (default: 50)')
    parser.add_argument('--threads', type=int, default=1,
                        help='Number of parallel workers for BAM processing (default: 1)')
    parser.add_argument('--output', default='viewer.html',
                        help='Output HTML file (default: viewer.html)')
    args = parser.parse_args()

    # Validate BAM args
    if args.bam and not args.bam_labels:
        args.bam_labels = [Path(b).stem for b in args.bam]
    if args.bam and len(args.bam_labels) != len(args.bam):
        print("Error: --bam-labels must match number of --bam files", file=sys.stderr)
        sys.exit(1)

    # ── Parse TEs ──
    all_elements = []
    te_map = {}
    for level, fasta_path in enumerate(args.te_fastas):
        print(f"Parsing TE fasta level {level}: {fasta_path}", file=sys.stderr)
        elements = parse_fasta_headers(fasta_path, level)
        all_elements.extend(elements)
        for el in elements:
            te_map[el.key] = el

    # ── Attach divergence data ──
    for level, div_path in enumerate(args.ltr_divergence):
        if level < len(args.te_fastas):
            print(f"Parsing divergence level {level}: {div_path}", file=sys.stderr)
            parse_ltr_divergence(div_path, te_map)

    # ── Attach domain data ──
    for level, dom_path in enumerate(args.domains):
        if level < len(args.te_fastas):
            print(f"Parsing domains level {level}: {dom_path}", file=sys.stderr)
            parse_domains(dom_path, te_map)

    # ── Build nesting ──
    if all_elements:
        print("Building nesting tree...", file=sys.stderr)
        build_nesting(all_elements)

    # ── Parse genes ──
    genes = []
    if args.gff3:
        print(f"Parsing GFF3: {args.gff3}", file=sys.stderr)
        genes = parse_gff3(args.gff3)

    print(f"Loaded {len(all_elements)} TE elements, {len(genes)} genes", file=sys.stderr)

    # ── Process BAMs ──
    coverage_data = {}
    transcript_data = {}
    bam_labels = args.bam_labels

    if args.bam:
        try:
            import pysam
        except ImportError:
            print("Error: pysam is required for BAM processing. Install with: pip install pysam",
                  file=sys.stderr)
            sys.exit(1)
        coverage_data, transcript_data = compute_coverage_and_transcripts(
            args.bam, args.bam_labels, genes, all_elements,
            bin_size=args.bin_size, threads=args.threads
        )

    # ── Serialise ──
    te_data   = json.dumps([te_to_dict(el) for el in all_elements], separators=(',', ':'))
    gene_data = json.dumps([gene_to_dict(g) for g in genes], separators=(',', ':'))
    cov_data  = json.dumps(coverage_data, separators=(',', ':'))
    tx_data   = json.dumps(transcript_data, separators=(',', ':'))
    labels_data = json.dumps(bam_labels, separators=(',', ':'))

    html = HTML_TEMPLATE
    html = html.replace('__TE_DATA__', te_data)
    html = html.replace('__GENE_DATA__', gene_data)
    html = html.replace('__COVERAGE_DATA__', cov_data)
    html = html.replace('__TRANSCRIPT_DATA__', tx_data)
    html = html.replace('__BAM_LABELS__', labels_data)

    out_path = Path(args.output)
    out_path.write_text(html, encoding='utf-8')
    print(f"\n✓ Viewer written to: {out_path.resolve()}", file=sys.stderr)
    print(f"  TEs:          {len(all_elements)}", file=sys.stderr)
    print(f"  Genes:        {len(genes)}", file=sys.stderr)
    print(f"  BAM samples:  {len(bam_labels)}", file=sys.stderr)
    for label in bam_labels:
        n_cov = sum(len(v) for v in coverage_data.get(label, {}).values())
        n_tx  = len(transcript_data.get(label, []))
        print(f"    {label}: {n_cov} coverage bins, {n_tx} transcript clusters", file=sys.stderr)
    kb = out_path.stat().st_size / 1024
    print(f"  Size:  {kb:.1f} KB", file=sys.stderr)


if __name__ == '__main__':
    main()
