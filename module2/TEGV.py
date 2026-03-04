#!/usr/bin/env python3
"""
te_viewer.py — Interactive nested TE + gene genome browser
Generates a self-contained HTML file with D3.js visualization.

Usage:
    python TEGV.py \
        --gff3 genes.gff3 \
        --te-fastas level0.fa level1.fa level2.fa \
        --ltr-divergence level0.div level1.div level2.div \
        --domains level0.domains level1.domains level2.domains \
        --output viewer.html

python TEGV.py --gff3 cotton_ltr_r1.work/cotton_ltr_r1.genic.gff  --te-fastas cotton_ltr_r1.ltrharvest.full_length.dedup.fa.rexdb-plant.cls.lib.fa cotton_ltr_r2.ltrharvest.full_length.dedup.fa.rexdb-plant.cls.lib.fa cotton_ltr_r3.ltrharvest.full_length.dedup.fa.rexdb-plant.cls.lib.fa cotton_ltr_r4.ltrharvest.full_length.dedup.fa.rexdb-plant.cls.lib.fa cotton_ltr_r5.ltrharvest.full_length.dedup.fa.rexdb-plant.cls.lib.fa --ltr-divergence cotton_ltr_r1_kmer2ltr_dedup cotton_ltr_r2_kmer2ltr_dedup  cotton_ltr_r3_kmer2ltr_dedup cotton_ltr_r4_kmer2ltr_dedup cotton_ltr_r5_kmer2ltr_dedup --domains cotton_ltr_r1.work/cotton_ltr_r1.ltrtools.intact_for_tesorter.fa.rexdb-plant.dom.gff3 cotton_ltr_r2.work/cotton_ltr_r2.ltrtools.intact_for_tesorter.fa.rexdb-plant.dom.gff3 cotton_ltr_r3.work/cotton_ltr_r3.ltrtools.intact_for_tesorter.fa.rexdb-plant.dom.gff3 cotton_ltr_r4.work/cotton_ltr_r4.ltrtools.intact_for_tesorter.fa.rexdb-plant.dom.gff3 cotton_ltr_r5.work/cotton_ltr_r5.ltrtools.intact_for_tesorter.fa.rexdb-plant.dom.gff3
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
    family: str        # e.g. "Ivana"
    superfamily: str   # e.g. "LTR/Copia"
    te_class: str      # e.g. "LTR"
    level: int         # 0=innermost
    ltr_len: int = 0
    k2p_dist: float = 0.0
    k2p_time: float = 0.0
    has_divergence: bool = False   # True if a divergence row was matched
    domains: list = field(default_factory=list)
    children: list = field(default_factory=list)
    visible_frags: list = field(default_factory=list)  # [(start,end),...]

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
    feature_type: str  # exon / CDS / UTR


@dataclass
class Gene:
    chrom: str
    start: int
    end: int
    name: str
    strand: str
    gene_id: str
    exons: list = field(default_factory=list)
    targets: list = field(default_factory=list)  # all Target= IDs merged into this gene


# ─────────────────────────────────────────────
#  Parsers
# ─────────────────────────────────────────────

def parse_fasta_headers(fasta_path, level):
    """Parse TE elements from FASTA headers like >Nem_chr1:4000-5000#LTR/Copia/Ivana"""
    elements = []
    with open(fasta_path, buffering=1 << 20) as f:
        for line in f:
            if not line or line[0] != '>':
                continue
            header = line[1:].rstrip().split()[0]
            # split coord part and classification
            if '#' not in header:
                print(f"Warning: skipping malformed header: {header}", file=sys.stderr)
                continue
            coord_part, classif = header.split('#', 1)
            # parse chrom:start-end
            m = re.match(r'^(.+):(\d+)-(\d+)$', coord_part)
            if not m:
                print(f"Warning: cannot parse coords in: {header}", file=sys.stderr)
                continue
            chrom = m.group(1)
            start = int(m.group(2))
            end = int(m.group(3))
            # parse classification e.g. LTR/Copia/Ivana or DNA/TIR/Helitron
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
    """
    Parse LTR divergence file and attach ltr_len + k2p stats to matching elements.
    Key format: Nem_chr1:4000-5000#LTR/Copia/Ivana
    Columns: LTR-RT  LTR_LEN  ALN_LEN  subs  transitions  transversions
             p-dist  p-time  JC69-dist  JC69-time  K2P-dist  K2P-time
    """
    with open(div_path, buffering=1 << 20) as f:
        for line in f:
            if not line or line[0] == '#':
                continue
            cols = line.split()
            if not cols:
                continue
            if len(cols) < 12:
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
    """
    Parse TEsorter-style domain GFF.
    Uses (chrom, start, end) index for O(1) TE lookup; pre-compiled
    regexes; buffered I/O; fast attr splitting. Previously O(N) per line.
    """
    # O(1) lookup index: (chrom, start, end) -> TEElement
    coord_index = {(te.chrom, te.start, te.end): te for te in te_map.values()}

    # Compile patterns once
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
    """Parse GFF3 for genes/mRNAs, exons, CDSs, UTRs.
    Handles both standard GFF3 (with gene features) and miniprot-style
    GFF3 where mRNA is the top-level feature (no parent gene line).
    Gene label priority: Target= > Name= > ID=
    After parsing, overlapping genes on the same chrom+strand are collapsed
    into a single representative gene; all Target IDs are preserved in .targets.
    """
    genes = {}         # gene_id -> Gene
    mrna_to_gene = {}  # mRNA_id -> gene_id

    with open(gff3_path, buffering=1 << 20) as f:
        for line in f:
            if not line or line[0] == '#':
                continue
            cols = line.split('\t')
            if len(cols) < 9:
                continue
            chrom, source, feature, start, end, score, strand, phase, attrs = cols
            start, end = int(start), int(end)
            # Fast attr parse: find '=' index directly, no strip needed for well-formed GFF3
            attr_dict = {}
            for attr in attrs.split(';'):
                eq = attr.find('=')
                if eq > 0:
                    attr_dict[attr[:eq].strip()] = attr[eq+1:].strip()

            gid    = attr_dict.get('ID', '')
            parent = attr_dict.get('Parent', '')

            # Determine display name: prefer Target (first token), then Name, then ID
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
                    # miniprot: mRNA is top-level, no parent gene line
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
    """
    Merge genes that overlap on the same chrom+strand into one representative.
    The merged gene spans the union of all overlapping genes.
    Its exons are the union of all CDS/exon blocks (deduplicated).
    All Target IDs are pooled into .targets; .name = first target ID.
    """
    from collections import defaultdict

    # group by (chrom, strand)
    buckets = defaultdict(list)
    for g in genes:
        buckets[(g.chrom, g.strand)].append(g)

    merged = []
    for (chrom, strand), group in buckets.items():
        # sort by start
        group.sort(key=lambda g: g.start)
        # Sweep line: O(N log N) — track running max-end of current cluster
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
    """
    Assign children and compute visible_frags for each element.
    Children of level-N elements are level-(N-1) elements whose coords
    fall entirely within the parent on the same chromosome.
    """
    by_level = defaultdict(list)
    for el in all_elements:
        by_level[el.level].append(el)

    max_level = max(by_level.keys()) if by_level else 0

    # Index children by (level, chrom) to skip cross-chromosome comparisons
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

    # compute visible fragments for each element
    for el in all_elements:
        compute_visible_frags(el)

    return all_elements


def compute_visible_frags(el):
    """
    Visible fragments = el's span minus its direct children's spans.
    """
    if not el.children:
        el.visible_frags = [(el.start, el.end)]
        return

    # sort children by start
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

  /* ── Header ── */
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

  /* ── Minimap ── */
  #minimap-container {
    background: var(--surface);
    border-bottom: 1px solid var(--border);
    padding: 6px 20px;
    flex-shrink: 0;
  }

  #minimap {
    width: 100%;
    height: 28px;
    cursor: crosshair;
  }

  /* ── Main canvas area ── */
  #viewer-wrap {
    flex: 1;
    overflow: hidden;
    position: relative;
  }

  #viewer-svg {
    width: 100%;
    height: 100%;
    cursor: grab;
  }
  #viewer-svg:active { cursor: grabbing; }

  /* ── Tooltip ── */
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
    max-width: 340px;
    line-height: 1.7;
  }
  #tooltip.visible { opacity: 1; }
  #tooltip .tt-title { font-weight: 600; color: var(--text-bright); font-size: 12px; margin-bottom: 4px; }
  #tooltip .tt-row { display: flex; gap: 12px; }
  #tooltip .tt-label { color: var(--text-dim); min-width: 70px; }

  /* ── Legend ── */
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

  /* ── Scale bar ── */
  .scale-bar text { fill: var(--text-dim); font-family: 'IBM Plex Mono', monospace; font-size: 10px; }
  .scale-bar line { stroke: var(--text-dim); }

  /* ── Track labels ── */
  .track-label { fill: var(--text-dim); font-family: 'IBM Plex Mono', monospace; font-size: 10px; text-transform: uppercase; letter-spacing: 0.06em; }

  /* ── Coord ruler ── */
  .ruler text { fill: var(--text-dim); font-family: 'IBM Plex Mono', monospace; font-size: 10px; }
  .ruler line { stroke: var(--border); }
  .ruler .major { stroke: var(--text-dim); }

  /* ── Status bar ── */
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

// ═══════════════════════════════════════════════════════
//  CONSTANTS & CONFIG
// ═══════════════════════════════════════════════════════
const RULER_H     = 28;
const GENE_H      = 18;
const GENE_ROW_H  = 26;
const TE_ROW_H    = 22;
const TE_BAR_H    = 14;
const LTR_SHADE   = 0.6;   // darken factor for LTR ends
const DOMAIN_H    = 8;
const PAD_LEFT    = 70;
const PAD_RIGHT   = 20;
const PAD_TOP     = 14;

// Family → color palette (qualitative, colorblind-friendly extended)
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

function internalColor(family) {
  return familyColor.get(family) || "#888";
}

function domainColor() { return "#fff"; }

// Gene colors
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
    d.chrom === chrom &&
    d.end > vStart && d.start < vEnd &&
    !hiddenFamilies.has(d.family)
  );
  const genes = RAW_GENES.filter(d =>
    d.chrom === chrom && d.end > vStart && d.start < vEnd
  );

  // max nesting level present
  const maxLevel = tes.length ? Math.max(...tes.map(d => d.level)) : 0;

  // Build rows: we lay out interleaved genes and TE stacks
  // TE rows: one row per nesting level (level 0 at top, maxLevel at bottom)
  // Genes: placed in their own rows using a greedy non-overlapping packer

  const teRowCount = maxLevel + 1;

  // Pack genes into non-overlapping rows
  const geneRows = packRows(genes);

  return { tes, genes, geneRows, teRowCount, maxLevel };
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
  const { tes, genes, geneRows, teRowCount } = layout;

  const xScale = d3.scaleLinear()
    .domain([viewStart, viewEnd])
    .range([PAD_LEFT, PAD_LEFT + TW]);

  // compute total height
  const geneBlockH = geneRows.length * GENE_ROW_H;
  const teBlockH   = teRowCount * TE_ROW_H;
  const totalContentH = PAD_TOP + RULER_H + geneBlockH + 10 + teBlockH + 40;
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

  // ── TE track label ──
  svg.append('text').attr('class','track-label')
    .attr('x', 4).attr('y', yOffset + 12).text('TEs');

  // ── TE rows: innermost (level 0) at top ──
  for (let lvl = 0; lvl <= layout.maxLevel; lvl++) {
    const rowTEs = tes.filter(d => d.level === lvl);
    const rowY = yOffset + lvl * TE_ROW_H;
    drawTERow(svg, rowTEs, xScale, rowY, lvl);
  }

  // ── Nesting connectors ──
  drawNestingConnectors(svg, tes, xScale, yOffset);

  // ── Viewport pan/zoom ──
  setupInteraction(svg, xScale, W, TW, svgH);

  // ── Status bar ──
  document.getElementById('statusbar').textContent =
    `${currentChrom}  ${fmt(viewStart)}–${fmt(viewEnd)}  (${fmt(viewEnd - viewStart)} bp)`;

  // ── Minimap ──
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
    Math.ceil(viewStart / step) * step,
    viewEnd,
    step
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

      // intron line (gene body)
      gg.append('line')
        .attr('x1', gx).attr('x2', gx + gw)
        .attr('y1', ry + GENE_H / 2).attr('y2', ry + GENE_H / 2)
        .attr('stroke', GENE_BODY_COLOR).attr('stroke-width', 1.5);

      // direction arrows along body
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

      // exons / CDS / UTR
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

      // gene body outline (if no exons, draw a box)
      if (gene.exons.length === 0) {
        gg.append('rect')
          .attr('x', gx).attr('y', ry)
          .attr('width', gw).attr('height', GENE_H)
          .attr('fill', GENE_EXON_COLOR).attr('rx', 2);
      }

      // gene label
      if (gw > 30) {
        gg.append('text')
          .attr('x', gx + gw / 2).attr('y', ry - 2)
          .attr('text-anchor','middle')
          .attr('fill','var(--text-dim)')
          .attr('font-family','IBM Plex Mono, monospace')
          .attr('font-size','9px')
          .text(gene.name.length > 16 ? gene.name.slice(0,14)+'…' : gene.name);
      }

      // tooltip interaction
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

  // level indicator
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
        // determine which part of the element this fragment covers
        const ltrLeft  = { s: te.start, e: te.start + te.ltr_len };
        const ltrRight = { s: te.end - te.ltr_len, e: te.end };

        // draw sub-segments of this fragment
        drawLTRFragment(fragG, fs, fe, te, ltrLeft, ltrRight, xScale, barY);
      } else {
        // plain bar
        fragG.append('rect')
          .attr('x', fx).attr('y', barY)
          .attr('width', fw).attr('height', TE_BAR_H)
          .attr('fill', internalColor(te.family))
          .attr('rx', 2).attr('opacity', 0.88);
      }
    });

    // protein domains (drawn on innermost fragments at level 0)
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

    // K2P label if visible
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

    // invisible hit target for tooltip
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
  // Split fragment into up to 3 segments: ltr_left_part | internal | ltr_right_part
  const segments = [];
  // left LTR overlap
  const ll_s = Math.max(fs, ltrLeft.s);
  const ll_e = Math.min(fe, ltrLeft.e);
  if (ll_e > ll_s) segments.push({ s: ll_s, e: ll_e, type: 'ltr' });
  // internal
  const int_s = Math.max(fs, ltrLeft.e);
  const int_e = Math.min(fe, ltrRight.s);
  if (int_e > int_s) segments.push({ s: int_s, e: int_e, type: 'internal' });
  // right LTR overlap
  const rl_s = Math.max(fs, ltrRight.s);
  const rl_e = Math.min(fe, ltrRight.e);
  if (rl_e > rl_s) segments.push({ s: rl_s, e: rl_e, type: 'ltr' });

  if (segments.length === 0) {
    // fallback: plain bar
    g.append('rect')
      .attr('x', xScale(fs)).attr('y', barY)
      .attr('width', Math.max(1, xScale(fe) - xScale(fs)))
      .attr('height', TE_BAR_H)
      .attr('fill', internalColor(te.family))
      .attr('rx', 2).attr('opacity', 0.88);
    return;
  }

  const isFirstSeg = idx => idx === 0;
  const isLastSeg  = (idx, arr) => idx === arr.length - 1;

  segments.forEach((seg, idx) => {
    const sx = xScale(seg.s);
    const sw = Math.max(1, xScale(seg.e) - xScale(seg.s));
    const color = seg.type === 'ltr' ? ltrColor(te.family) : internalColor(te.family);
    const h = seg.type === 'ltr' ? TE_BAR_H : TE_BAR_H - 2;
    const dy = seg.type === 'ltr' ? 0 : 1;
    // rounded corners only on outer ends
    g.append('rect')
      .attr('x', sx).attr('y', barY + dy)
      .attr('width', sw).attr('height', h)
      .attr('fill', color)
      .attr('rx', (isFirstSeg(idx) || isLastSeg(idx, segments)) ? 2 : 0)
      .attr('opacity', 0.9);
  });
}

// ═══════════════════════════════════════════════════════
//  NESTING CONNECTORS
// ═══════════════════════════════════════════════════════
function drawNestingConnectors(svg, tes, xScale, teTrackY) {
  // Draw subtle vertical lines connecting parent fragments to child elements
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
      // connect left edge of child to parent
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

  // background
  svg.append('rect').attr('x',PAD_LEFT).attr('y',4).attr('width',W-PAD_LEFT).attr('height',H-8)
     .attr('fill','var(--surface2)').attr('rx',2);

  // TE density blobs
  const allChromTEs = RAW_TES.filter(d => d.chrom === currentChrom && !hiddenFamilies.has(d.family));
  allChromTEs.forEach(te => {
    const tx = mx(te.start);
    const tw = Math.max(1, mx(te.end) - mx(te.start));
    svg.append('rect').attr('x', tx).attr('y', 4 + H*0.25)
       .attr('width', tw).attr('height', H * 0.5)
       .attr('fill', internalColor(te.family)).attr('opacity', 0.4);
  });

  // gene marks
  const allChromGenes = RAW_GENES.filter(d => d.chrom === currentChrom);
  allChromGenes.forEach(g => {
    const gx = mx(g.start);
    svg.append('line').attr('x1',gx).attr('x2',gx).attr('y1',4).attr('y2',H-4)
       .attr('stroke', GENE_ARROW_COLOR).attr('stroke-width',1).attr('opacity',0.6);
  });

  // viewport highlight
  const vx = mx(viewStart);
  const vw = Math.max(4, mx(viewEnd) - mx(viewStart));
  const vp = svg.append('rect')
    .attr('x', vx).attr('y', 2)
    .attr('width', vw).attr('height', H - 4)
    .attr('fill','rgba(74,158,255,0.12)')
    .attr('stroke','var(--accent)').attr('stroke-width',1)
    .attr('rx',2).attr('cursor','ew-resize');

  // drag viewport on minimap
  const drag = d3.drag().on('drag', (event) => {
    const clickX = Math.max(PAD_LEFT, Math.min(PAD_LEFT + W - PAD_LEFT, event.x));
    const genomicPos = mx.invert(clickX);
    const span = viewEnd - viewStart;
    viewStart = Math.max(0, genomicPos - span / 2);
    viewEnd = viewStart + span;
    renderAll();
  });
  vp.call(drag);

  // click to jump
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
  // Drag state lives outside all event handlers so it survives
  // cursor moving over child elements (TE bars, gene rects, etc.)
  // which previously fired mouseout and killed the drag.
  let dragStartX   = null;
  let dragStartView = null;
  let isDragging   = false;

  // Mousedown on the SVG starts a potential drag
  svg.on('mousedown', (event) => {
    if (event.button !== 0) return;
    dragStartX    = event.clientX;
    dragStartView = [viewStart, viewEnd];
    isDragging    = false;
    event.preventDefault();
  });

  // Mousemove and mouseup are on window so they fire even when the
  // cursor slides over child SVG elements or leaves the SVG entirely.
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
    // zoom toward the mouse position on the genome
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

  // gene item
  const geneItem = document.createElement('div');
  geneItem.className = 'legend-item';
  geneItem.innerHTML = `<div class="legend-swatch" style="background:${GENE_EXON_COLOR}"></div><span class="legend-label">Gene</span>`;
  container.appendChild(geneItem);

  // TE families
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

  // LTR indicator
  const ltrItem = document.createElement('div');
  ltrItem.className = 'legend-item';
  ltrItem.innerHTML = `<div style="display:flex;gap:1px;align-items:center">
    <div style="width:5px;height:10px;background:#555;border-radius:2px 0 0 2px"></div>
    <div style="width:14px;height:8px;background:#666;margin-top:1px"></div>
    <div style="width:5px;height:10px;background:#555;border-radius:0 2px 2px 0"></div>
  </div><span class="legend-label" style="margin-left:5px">LTR | internal | LTR</span>`;
  container.appendChild(ltrItem);
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
    parser.add_argument('--output', default='viewer.html',
                        help='Output HTML file (default: viewer.html)')
    args = parser.parse_args()

    # ── Parse TEs ──
    all_elements = []
    te_map = {}  # key → TEElement
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

    # ── Serialise ──
    te_data   = json.dumps([te_to_dict(el) for el in all_elements], separators=(',', ':'))
    gene_data = json.dumps([gene_to_dict(g) for g in genes], separators=(',', ':'))

    html = HTML_TEMPLATE.replace('__TE_DATA__', te_data).replace('__GENE_DATA__', gene_data)

    out_path = Path(args.output)
    out_path.write_text(html, encoding='utf-8')
    print(f"\n✓ Viewer written to: {out_path.resolve()}", file=sys.stderr)
    print(f"  TEs:   {len(all_elements)}", file=sys.stderr)
    print(f"  Genes: {len(genes)}", file=sys.stderr)
    kb = out_path.stat().st_size / 1024
    print(f"  Size:  {kb:.1f} KB", file=sys.stderr)


if __name__ == '__main__':
    main()
