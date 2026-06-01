#!/usr/bin/env python3
"""Per-cluster LTR-RT phylogenies with structure glyphs at the tips.

Inputs : Kmer2LTR consensus LTR FASTA, an mmseqs cluster TSV, and one or more
         depth*_ltr.tsv annotation files. Output: one tree+glyph figure per
         non-singleton cluster. See docs/superpowers/specs/2026-05-29-*.md.


Pipeline : 
# De novo LTR-RT annotation.
bash synLTR/module2/ltrharvest_wrapper2.sh \
  --genome   synLTR/test/Athal_tair10_chr2.fa.gz \
  --proteins synLTR/test/Athal.pep.gz \
  --threads  250 \
  --terminate_count 10

# Merge nest-level fastas.
cat Athal_tair10_chr2_LTRs_depth*_ltr.fa > Athal_tair10_chr2_LTRs_all.fa

# Classify LTR region and cluster consensus LTRs. 
python Kmer2LTR/Kmer2LTR.py \
  -i Athal_tair10_chr2_LTRs_all.fa \
  --ltr-consensus --ltr-cluster --internal-fasta

# MSA (mafft), trimming (trimal), ML phylogeny (IQTREE).
# Unrooted trees of LTR clusters. Each cluster is its own treefile. 
python synLTR/module2/ltr_cluster_phylo.py \
  --consensus Athal_tair10_chr2_LTRs_all.LTRs.alns.consensus.fa \
  --clusters  Athal_tair10_chr2_LTRs_all.LTRs.alns.consensus_id0.70_cluster.tsv \
  --annot     Athal_tair10_chr2_LTRs_depth*_ltr.tsv \
  --out-dir   cluster_trees_id0.70 \
  --jobs 8 --threads 1 \
  -v
"""
import argparse
import math
import os
import re
import subprocess
import sys
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

# Reuse the existing structure-plot helpers (import-safe; top-level code guarded).
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ltrharvest_plot_struct import (  # noqa: E402
    Feature, FEATURE_COLORS, DEFAULT_COLOR, normalize_protein_name,
    build_legend_handles,
)

# Okabe-Ito colorblind-safe palette overriding the module's domain colors.
OKABE_ITO = {
    "LTR":  "#333333",
    "GAG":  "#E69F00",  # orange
    "PROT": "#CC79A7",  # reddish purple
    "RT":   "#0072B2",  # blue
    "RH":   "#D55E00",  # vermillion
    "INT":  "#009E73",  # bluish green
    "CH":   "#56B4E9",  # sky blue
    "CHD":  "#56B4E9",
    "CHDCR":"#882255",
    "ARH":  "#F0E442",  # yellow
    "ENDO": "#999999",
}
BACKBONE_COLOR = "#DDDDDD"

CANONICAL_ORDER = {  # internal domain order 5'->3' on the + strand
    "Gypsy": ["GAG", "PROT", "RT", "RH", "INT"],
    "Copia": ["GAG", "PROT", "INT", "RT", "RH"],
}


@dataclass
class Cluster:
    index: int
    rep: str
    members: List[str]

    @property
    def size(self) -> int:
        return len(self.members)


@dataclass
class ElementAnnot:
    eid: str
    chrom: str
    start: int          # 1-based genomic, inclusive
    end: int            # 1-based genomic, inclusive
    ltr_len: int
    k2p: Optional[float]
    order: str
    superfamily: str
    clade: str
    domains: List[Tuple[str, int, int]]   # (gene, gStart, gEnd) genomic 1-based
    inners: List[Tuple[int, int]]         # (iStart, iEnd) genomic 1-based
    is_inner: bool


@dataclass
class Glyph:
    el_len: int
    ltr_len: int
    domains: List[Feature]                 # element-relative, 1-based start..end
    gaps: List[Tuple[int, int]]            # element-relative, 0-based [s, e)
    strand: str


# ---------------------------------------------------------------------------
# Labels
# ---------------------------------------------------------------------------
def sanitize_label(s: str) -> str:
    """Make a label safe for Newick (':' is reserved; also strip '#','/')."""
    return re.sub(r"[^A-Za-z0-9_.\-]", "_", s)


class LabelMap:
    """Bidirectional element-id <-> safe-label map with collision handling."""

    def __init__(self, ids: List[str]):
        self._to_safe: Dict[str, str] = {}
        self._to_full: Dict[str, str] = {}
        seen: Dict[str, int] = {}
        for i in ids:
            base = sanitize_label(i)
            if base in seen:
                seen[base] += 1
                safe = f"{base}__{seen[base]}"
            else:
                seen[base] = 0
                safe = base
            self._to_safe[i] = safe
            self._to_full[safe] = i

    def to_safe(self, full: str) -> str:
        return self._to_safe[full]

    def to_full(self, safe: str) -> str:
        return self._to_full[safe]


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------
_ID_RE = re.compile(r"^(?P<chrom>.+):(?P<start>\d+)-(?P<end>\d+)#(?P<cls>.+)$")


def parse_element_id(eid: str) -> Tuple[str, int, int, str, str, str]:
    m = _ID_RE.match(eid)
    if not m:
        raise ValueError(f"unparseable element id: {eid!r}")
    cls_parts = m.group("cls").split("/")
    order = cls_parts[0] if len(cls_parts) > 0 else "unknown"
    sf = cls_parts[1] if len(cls_parts) > 1 else "unknown"
    clade = cls_parts[2] if len(cls_parts) > 2 else "unknown"
    return (m.group("chrom"), int(m.group("start")), int(m.group("end")),
            order, sf, clade)


def parse_clusters(path: str, min_size: int = 2) -> List[Cluster]:
    """Read mmseqs rep<TAB>member TSV. Drop clusters smaller than min_size
    (singletons are size 1). Returned largest-first with a stable 0-based index."""
    members: Dict[str, List[str]] = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            rep, _, mem = line.partition("\t")
            if not mem:
                raise ValueError(f"cluster TSV line is not 2-column: {line!r}")
            members.setdefault(rep, [])
            if mem not in members[rep]:
                members[rep].append(mem)
    kept = [(rep, mem) for rep, mem in members.items() if len(mem) >= min_size]
    kept.sort(key=lambda rm: (-len(rm[1]), rm[0]))
    return [Cluster(index=i, rep=rep, members=mem)
            for i, (rep, mem) in enumerate(kept)]


def parse_domains(field: str) -> List[Tuple[str, int, int]]:
    field = field.strip()
    if field in (".", ""):
        return []
    out = []
    for tok in field.split(";"):
        tok = tok.strip()
        if not tok:
            continue
        gene = tok.split("|", 1)[0]
        coords = tok.split("@", 1)[1]
        gs, ge = coords.split("-")
        out.append((gene, int(gs), int(ge)))
    return out


def parse_nest(field: str) -> Tuple[List[Tuple[int, int]], bool]:
    field = field.strip()
    if field in (".", ""):
        return [], False
    inners: List[Tuple[int, int]] = []
    is_inner = False
    for tok in field.split(";"):
        tok = tok.strip()
        if not tok:
            continue
        role, _, locus = tok.partition(":")
        if role == "nest-outer":
            # locus = "chrom:start-end" (no #class suffix); reuse the id parser.
            _, s, e, *_ = parse_element_id(
                locus if "#" in locus else locus + "#LTR/x/y")
            inners.append((s, e))
        elif role == "nest-inner":
            is_inner = True
    return inners, is_inner


def parse_annot(paths: List[str]) -> Dict[str, ElementAnnot]:
    """Concatenate depth*_ltr.tsv files into {id: ElementAnnot}.

    Header names 17 columns but data rows have 19 (two unnamed end5p/start3p
    before tsd/domains/nest_status), so the tail is read by NEGATIVE index:
    nest_status=cols[-1], domains=cols[-2], tsd=cols[-3]. Named positional:
    name=cols[0], LTR_len=cols[1], K2P_d=cols[10].
    """
    out: Dict[str, ElementAnnot] = {}
    for path in paths:
        with open(path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 17:
                    continue  # malformed row: warn-and-skip
                eid = cols[0]
                chrom, start, end, order, sf, clade = parse_element_id(eid)
                try:
                    k2p = float(cols[10])
                except (ValueError, IndexError):
                    k2p = None
                inners, is_inner = parse_nest(cols[-1])
                out[eid] = ElementAnnot(
                    eid=eid, chrom=chrom, start=start, end=end,
                    ltr_len=int(cols[1]), k2p=k2p, order=order,
                    superfamily=sf, clade=clade,
                    domains=parse_domains(cols[-2]),
                    inners=inners, is_inner=is_inner)
    return out


def infer_strand(a: ElementAnnot) -> str:
    """'+'/'-' from genomic domain order vs the superfamily canonical order.
    Falls back to '+' when undeterminable (logged by caller)."""
    canon = CANONICAL_ORDER.get(a.superfamily)
    if not canon:
        return "+"
    first_start: Dict[str, int] = {}
    for gene, gs, _ge in a.domains:
        name = normalize_protein_name(gene)
        first_start.setdefault(name, gs)
    present = [d for d in canon if d in first_start]
    if len(present) < 2:
        return "+"
    lo, hi = present[0], present[-1]   # lo precedes hi in canonical 5'->3'
    return "+" if first_start[lo] < first_start[hi] else "-"


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------
def element_geometry(a: ElementAnnot, strand: str) -> Glyph:
    el_len = a.end - a.start + 1
    domains: List[Feature] = []
    for gene, gs, ge in a.domains:
        name = normalize_protein_name(gene)
        s_rel = gs - a.start + 1            # 1-based element-relative
        e_rel = ge - a.start + 1
        domains.append(Feature(name, s_rel, e_rel))
    gaps: List[Tuple[int, int]] = []
    for i_start, i_end in a.inners:
        gs = i_start - a.start              # 0-based left edge
        ge = i_end - a.start + 1            # exclusive right edge
        gaps.append((max(0, gs), min(el_len, ge)))
    gaps.sort()
    if strand == "-":
        domains = [Feature(d.name, el_len - d.end + 1, el_len - d.start + 1)
                   for d in domains]
        gaps = sorted((el_len - e, el_len - s) for s, e in gaps)
    return Glyph(el_len=el_len, ltr_len=a.ltr_len, domains=domains,
                 gaps=gaps, strand=strand)


# ---------------------------------------------------------------------------
# Sequence IO
# ---------------------------------------------------------------------------
def fetch_sequences(consensus_fa: str, ids: List[str]) -> Dict[str, str]:
    """Random-access fetch by header via pyfastx (never loads the whole file)."""
    import pyfastx
    fa = pyfastx.Fasta(consensus_fa, build_index=True)
    out: Dict[str, str] = {}
    for i in ids:
        try:
            out[i] = fa[i].seq
        except KeyError:
            sys.stderr.write(f"[warn] id not found in consensus FASTA: {i}\n")
    return out


def write_fasta(seqs: Dict[str, str], label_map: "LabelMap", path: str) -> None:
    with open(path, "w") as fh:
        for full, seq in seqs.items():
            fh.write(f">{label_map.to_safe(full)}\n{seq.upper()}\n")


# ---------------------------------------------------------------------------
# Alignment + trimming
# ---------------------------------------------------------------------------
def run_mafft(in_fa: str, out_fa: str, threads: int = 1) -> None:
    with open(out_fa, "w") as out:
        subprocess.run(
            ["mafft", "--auto", "--preservecase", "--thread", str(threads), in_fa],
            check=True, stdout=out, stderr=subprocess.DEVNULL)
    _uppercase_fasta(out_fa)


def run_trimal(in_fa: str, out_fa: str, mode: str = "gappyout") -> None:
    if mode == "none":
        import shutil as _sh
        _sh.copyfile(in_fa, out_fa)
        return
    flag = {"gappyout": "-gappyout", "automated1": "-automated1"}[mode]
    subprocess.run(["trimal", "-in", in_fa, "-out", out_fa, flag], check=True,
                   stderr=subprocess.DEVNULL)


def _uppercase_fasta(path: str) -> None:
    with open(path) as fh:
        lines = [l if l.startswith(">") else l.upper() for l in fh]
    with open(path, "w") as fh:
        fh.writelines(lines)


# ---------------------------------------------------------------------------
# Trees
# ---------------------------------------------------------------------------
def jc_pdistance(a: str, b: str) -> float:
    """Jukes-Cantor distance over comparable ACGT columns; IUPAC/gaps ignored."""
    acgt = set("ACGT")
    n = diff = 0
    for x, y in zip(a.upper(), b.upper()):
        if x in acgt and y in acgt:
            n += 1
            if x != y:
                diff += 1
    if n == 0:
        return 0.0
    p = diff / n
    if p >= 0.75:
        return 3.0  # saturated; cap
    return -0.75 * math.log(1 - 4.0 / 3.0 * p)


def write_two_tip_tree(safe_a: str, safe_b: str, dist: float, out_nwk: str) -> None:
    half = dist / 2.0
    with open(out_nwk, "w") as fh:
        fh.write(f"({safe_a}:{half:g},{safe_b}:{half:g});\n")


def _read_alignment(path: str) -> List[Tuple[str, str]]:
    recs, name, seq = [], None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(seq)))
                name = line[1:].strip().split()[0]
                seq = []
            else:
                seq.append(line.strip())
    if name is not None:
        recs.append((name, "".join(seq)))
    return recs


def run_iqtree(aln: str, prefix: str, n: int, model: str, ufboot: int,
               threads: int) -> Tuple[str, str]:
    cmd = ["iqtree", "-s", aln, "-m", model, "-nt", str(threads),
           "-keep-ident", "-redo", "-quiet", "-pre", prefix]
    support = "none"
    if n >= 4:
        cmd += ["-bb", str(ufboot)]
        support = "UFBoot"
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL)
    return prefix + ".treefile", support


def run_veryfasttree(aln: str, out_nwk: str, threads: int) -> Tuple[str, str]:
    common = ["-gamma", "-quiet", "-nopr", "-threads", str(threads)]
    for extra in (["-gtr"], []):   # GTR, then JC fallback
        try:
            with open(out_nwk, "w") as out:
                subprocess.run(["VeryFastTree", "-nt"] + extra + common + [aln],
                               check=True, stdout=out, stderr=subprocess.DEVNULL)
            if os.path.getsize(out_nwk) > 0:
                return out_nwk, "SH-like"
        except subprocess.CalledProcessError:
            continue
    raise RuntimeError(f"VeryFastTree failed on {aln}")


def build_tree(aln: str, prefix: str, n: int, method: str, max_ml_tips: int,
               model: str, ufboot: int, threads: int) -> Tuple[str, str]:
    if n == 2:
        (na, sa), (nb, sb) = _read_alignment(aln)
        out = prefix + ".nwk"
        write_two_tip_tree(na, nb, jc_pdistance(sa, sb), out)
        return out, "none"
    use_vft = method == "veryfasttree" or (method == "auto" and n > max_ml_tips)
    if use_vft:
        if method == "auto":
            sys.stderr.write(
                f"[info] cluster n={n} > max-ml-tips={max_ml_tips}: "
                f"using VeryFastTree (IUPAC-blind) instead of IQ-TREE\n")
        return run_veryfasttree(aln, prefix + ".nwk", threads)
    return run_iqtree(aln, prefix, n, model, ufboot, threads)


# ---------------------------------------------------------------------------
# Glyph drawing
# ---------------------------------------------------------------------------
def _solid_segments(el_len: int, gaps: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    segs, cur = [], 0
    for s, e in sorted(gaps):
        if s > cur:
            segs.append((cur, s))
        cur = max(cur, e)
    if cur < el_len:
        segs.append((cur, el_len))
    return segs or [(0, el_len)]


def draw_tip_glyph(ax, y: float, glyph: Glyph, scale: float = 1.0,
                   height: float = 0.62) -> None:
    import matplotlib.pyplot as plt
    bb_h = height * 0.45

    def X(v):
        return v * scale

    # backbone thick bar on solid segments
    for s, e in _solid_segments(glyph.el_len, glyph.gaps):
        ax.add_patch(plt.Rectangle((X(s), y - bb_h / 2), X(e - s), bb_h,
                                   facecolor=BACKBONE_COLOR, edgecolor="none",
                                   zorder=1))
    # nested inners as thin "intron" lines
    for s, e in glyph.gaps:
        ax.plot([X(s), X(e)], [y, y], color="black", linewidth=0.8, zorder=2)
    # LTR boxes (both ends)
    ltr_col = FEATURE_COLORS["LTR"]
    ax.add_patch(plt.Rectangle((X(0), y - height / 2), X(glyph.ltr_len), height,
                               facecolor=ltr_col, edgecolor="black",
                               linewidth=0.4, zorder=3))
    ax.add_patch(plt.Rectangle((X(glyph.el_len - glyph.ltr_len), y - height / 2),
                               X(glyph.ltr_len), height, facecolor=ltr_col,
                               edgecolor="black", linewidth=0.4, zorder=3))
    # domains
    for p in glyph.domains:
        col = FEATURE_COLORS.get(p.name, DEFAULT_COLOR)
        x = X(max(0, p.start - 1))
        w = X(max(1, p.end - p.start + 1))
        ax.add_patch(plt.Rectangle((x, y - height / 2), w, height,
                                   facecolor=col, edgecolor="black",
                                   linewidth=0.4, zorder=4))


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
def layout_phylogram(tree):
    """Return (xs, ys): clade -> root-distance x and leaf-order y."""
    xs = tree.depths()
    if not xs or max(xs.values()) == 0:
        xs = tree.depths(unit_branch_lengths=True)
    ys = {}
    for i, lf in enumerate(tree.get_terminals()):
        ys[lf] = float(i)

    def assign(clade):
        if clade.is_terminal():
            return ys[clade]
        cy = [assign(c) for c in clade.clades]
        ys[clade] = sum(cy) / len(cy)
        return ys[clade]

    assign(tree.root)
    return xs, ys


def _draw_tree(ax, tree):
    xs, ys = layout_phylogram(tree)

    def recur(clade):
        x0 = xs[clade]
        y0 = ys[clade]
        for c in clade.clades:
            x1 = xs[c]
            y1 = ys[c]
            ax.plot([x0, x1], [y1, y1], color="black", linewidth=0.8)  # horizontal
            recur(c)
        if clade.clades:
            child_ys = [ys[c] for c in clade.clades]
            ax.plot([x0, x0], [min(child_ys), max(child_ys)],
                    color="black", linewidth=0.8)                      # vertical
            if clade.confidence is not None:
                ax.text(x0, y0 + 0.12, f"{int(round(float(clade.confidence)))}",
                        fontsize=6, ha="center", va="bottom", color="#555555")

    recur(tree.root)
    return xs, ys


def render_cluster_figure(nwk_path: str, label_map: "LabelMap",
                          glyphs: Dict[str, Glyph],
                          k2p: Dict[str, Optional[float]],
                          caption: str, out_pdf: str, out_svg: str,
                          show_k2p: bool = True, glyph_scale: str = "true") -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    from Bio import Phylo

    # apply colorblind-safe palette to the reused module dict
    FEATURE_COLORS.update(OKABE_ITO)

    tree = Phylo.read(nwk_path, "newick")
    try:
        tree.root_at_midpoint()
    except (UnboundLocalError, ValueError):
        # Bio.Phylo's midpoint rooting fails on a zero-divergence tree (all
        # branch lengths ~0, e.g. a cluster of near-identical recent copies).
        # Leave the tree as the builder rooted it.
        pass
    tree.ladderize()
    for t in tree.get_terminals():
        t.name = label_map.to_full(t.name)

    n = len(tree.get_terminals())
    fig = plt.figure(figsize=(12, max(2.2, 0.5 * n + 1.3)), dpi=300)
    # columns: tree | tip labels | glyphs (with a 5'-3' K2P value column at right)
    gs = GridSpec(1, 3, width_ratios=[1.0, 1.35, 1.9], wspace=0.04, figure=fig)
    ax_tree = fig.add_subplot(gs[0, 0])
    ax_mid = fig.add_subplot(gs[0, 1], sharey=ax_tree)
    ax_gly = fig.add_subplot(gs[0, 2], sharey=ax_tree)

    _xs, ys = _draw_tree(ax_tree, tree)
    yvals = {t.name: ys[t] for t in tree.get_terminals()}

    def _short(full: str) -> str:
        # "chrom:start-end#order/super/clade" -> "chrom:start-end  super/clade"
        try:
            chrom, s, e, _o, sf, cl = parse_element_id(full)
            return f"{chrom}:{s}-{e}  {sf}/{cl}"
        except ValueError:
            return full

    ax_mid.set_xlim(0, 1)
    max_len = max(g.el_len for g in glyphs.values())
    k2p_x = max_len * 1.07          # x of the K2P value column (just right of glyphs)
    for full, y in yvals.items():
        g = glyphs[full]
        scale = 1.0 if glyph_scale == "true" else (max_len / g.el_len)
        draw_tip_glyph(ax_gly, y, g, scale=scale)
        ax_mid.text(0.0, y, _short(full), ha="left", va="center", fontsize=6)
        if show_k2p:
            kv = k2p.get(full)
            ax_gly.text(k2p_x, y, f"{kv:.4f}" if kv is not None else "NA",
                        ha="left", va="center", fontsize=6, family="monospace")

    ax_tree.set_ylim(-0.9, n - 0.1)
    ax_tree.set_xlabel("substitutions / site", fontsize=7)
    ax_tree.tick_params(labelsize=6)
    ax_tree.spines[["top", "right", "left"]].set_visible(False)
    ax_tree.set_yticks([])
    ax_mid.axis("off")
    ax_gly.set_xlim(-0.02 * max_len, max_len * (1.22 if show_k2p else 1.04))
    ax_gly.set_yticks([])
    ax_gly.tick_params(labelsize=6)
    ax_gly.spines[["top", "right", "left"]].set_visible(False)
    ax_gly.set_xlabel("element length (bp)", fontsize=7)
    if show_k2p:
        ax_gly.text(k2p_x, n - 0.45, "5'-3' K2P", ha="left", va="bottom",
                    fontsize=6, fontweight="bold")

    # domain legend (floated just above the glyph panel so it never covers a tip)
    feats = sorted({d.name for g in glyphs.values() for d in g.domains} | {"LTR"})
    ax_gly.legend(handles=build_legend_handles(feats), loc="lower right",
                  bbox_to_anchor=(1.0, 1.01), fontsize=6, frameon=False, ncol=4)

    fig.text(0.99, 0.005, caption, fontsize=6, color="#777777",
             ha="right")  # corner annotation, not a title
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def cluster_basename(c: Cluster) -> str:
    _, s, e, _o, sf, clade = parse_element_id(c.rep)
    tag = sanitize_label(f"{sf}_{clade}_{s}-{e}")
    return f"c{c.index:02d}_{tag}_n{c.size}"


def process_cluster(c: Cluster, consensus_fa: str,
                    annot: Dict[str, ElementAnnot], out_dir: str,
                    opts: argparse.Namespace) -> Optional[dict]:
    base = cluster_basename(c)
    prefix = os.path.join(out_dir, base)
    out_pdf, out_svg = prefix + ".pdf", prefix + ".svg"
    if os.path.exists(out_pdf):                      # resumable
        return {"cluster": base, "status": "skipped-exists", "pdf": out_pdf}
    members = [m for m in c.members if m in annot]
    if len(members) < opts.min_cluster_size:
        sys.stderr.write(
            f"[warn] {base}: <{opts.min_cluster_size} annotated members; skip\n")
        return None
    lm = LabelMap(members)
    fa = prefix + ".members.fa"
    write_fasta(fetch_sequences(consensus_fa, members), lm, fa)
    aln, trim = prefix + ".aln.fa", prefix + ".trim.fa"
    run_mafft(fa, aln, threads=opts.threads)
    run_trimal(aln, trim, mode=opts.trim)
    nwk, support = build_tree(trim, prefix, n=len(members), method=opts.tree,
                              max_ml_tips=opts.max_ml_tips, model=opts.iqtree_model,
                              ufboot=opts.ufboot, threads=opts.threads)
    glyphs, k2p = {}, {}
    for m in members:
        a = annot[m]
        glyphs[m] = element_geometry(a, infer_strand(a))
        k2p[m] = a.k2p
    render_cluster_figure(nwk, lm, glyphs, k2p, f"{base}  support={support}",
                          out_pdf, out_svg, show_k2p=opts.show_k2p,
                          glyph_scale=opts.glyph_scale)
    kvals = [v for v in k2p.values() if v is not None]
    return {"cluster": base, "rep": c.rep, "size": c.size, "method": opts.tree,
            "support": support,
            "mean_k2p": (sum(kvals) / len(kvals) if kvals else ""),
            "pdf": out_pdf}


def _write_manifest(rows: List[dict], path: str) -> None:
    cols = ["cluster", "rep", "size", "method", "support", "mean_k2p", "pdf"]
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


def _write_combined_pdf(rows: List[dict], path: str) -> None:
    # Concatenate the per-cluster PDFs (vector, no re-render) via pypdf if
    # available; otherwise warn and skip.
    try:
        from pypdf import PdfWriter
    except ImportError:
        sys.stderr.write("[warn] pypdf not available; skipping combined PDF\n")
        return
    w = PdfWriter()
    for r in rows:
        p = r.get("pdf")
        if p and os.path.exists(p):
            w.append(p)
    with open(path, "wb") as fh:
        w.write(fh)


def _write_memo(opts: argparse.Namespace, rows: List[dict]) -> None:
    path = os.path.join(opts.out_dir, "memo.md")
    n_made = sum(1 for r in rows if r.get("status") != "skipped-exists")
    with open(path, "w") as fh:
        fh.write("# memo: ltr_cluster_phylo\n\n")
        fh.write("date: 2026-05-29\n")
        fh.write(f"consensus: {opts.consensus}\n")
        fh.write(f"clusters: {opts.clusters}\n")
        fh.write(f"annot: {' '.join(opts.annot)}\n")
        fh.write(f"min_cluster_size: {opts.min_cluster_size}\n")
        fh.write(f"tree: {opts.tree} (max_ml_tips={opts.max_ml_tips}, "
                 f"model={opts.iqtree_model}, ufboot={opts.ufboot})\n")
        fh.write(f"trim: {opts.trim}; glyph_scale: {opts.glyph_scale}; "
                 f"show_k2p: {opts.show_k2p}\n")
        fh.write(f"clusters_processed: {len(rows)} (figures: {n_made})\n\n")
        fh.write("```\n")
        fh.write("python ltr_cluster_phylo.py --consensus <fa> --clusters <tsv> "
                 "--annot <depth0_ltr.tsv> [...] --out-dir <dir>\n")
        fh.write("```\n")


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(
        description="Per-cluster LTR-RT phylogenies with structure glyphs at tips.")
    p.add_argument("--consensus", required=True, help="Kmer2LTR consensus LTR FASTA")
    p.add_argument("--clusters", required=True, help="mmseqs rep<TAB>member TSV")
    p.add_argument("--annot", required=True, nargs="+",
                   help="depth*_ltr.tsv annotation file(s)")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--min-cluster-size", type=int, default=2)
    p.add_argument("--tree", choices=["auto", "iqtree", "veryfasttree"], default="auto")
    p.add_argument("--max-ml-tips", type=int, default=500)
    p.add_argument("--iqtree-model", default="GTR+G")
    p.add_argument("--ufboot", type=int, default=1000)
    p.add_argument("--trim", choices=["gappyout", "automated1", "none"],
                   default="gappyout")
    p.add_argument("--glyph-scale", choices=["true", "normalized"], default="true")
    p.add_argument("--show-k2p", dest="show_k2p", action="store_true", default=True,
                   help="annotate each tip with its 5'-3' LTR K2P value (default)")
    p.add_argument("--no-k2p", dest="show_k2p", action="store_false",
                   help="omit the per-element K2P value")
    p.add_argument("--combined-pdf", dest="combined_pdf", action="store_true",
                   default=True)
    p.add_argument("--no-combined-pdf", dest="combined_pdf", action="store_false")
    p.add_argument("--jobs", type=int, default=1)
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("-v", "--verbose", action="store_true")
    opts = p.parse_args(argv)

    os.makedirs(opts.out_dir, exist_ok=True)
    clusters = parse_clusters(opts.clusters, min_size=opts.min_cluster_size)
    annot = parse_annot(opts.annot)
    if opts.verbose:
        sys.stderr.write(f"[info] {len(clusters)} clusters >= size "
                         f"{opts.min_cluster_size}; {len(annot)} annotated elements\n")

    rows: List[dict] = []
    if opts.jobs > 1:
        from concurrent.futures import ProcessPoolExecutor, as_completed
        with ProcessPoolExecutor(max_workers=opts.jobs) as ex:
            futs = {ex.submit(process_cluster, c, opts.consensus, annot,
                              opts.out_dir, opts): c for c in clusters}
            for fut in as_completed(futs):
                r = fut.result()
                if r:
                    rows.append(r)
    else:
        for c in clusters:
            r = process_cluster(c, opts.consensus, annot, opts.out_dir, opts)
            if r:
                rows.append(r)
            if opts.verbose:
                sys.stderr.write(f"[info] done {cluster_basename(c)}\n")

    rows.sort(key=lambda r: r.get("cluster", ""))
    _write_manifest(rows, os.path.join(opts.out_dir, "manifest.tsv"))
    if opts.combined_pdf:
        _write_combined_pdf(rows, os.path.join(opts.out_dir, "all_clusters_trees.pdf"))
    _write_memo(opts, rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
