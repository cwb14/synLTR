"""`ltrtrace sim pop`: haploid genomes A1..An, B1..Bn, C1..Cn and the outgroup O1 evolved on a real
backbone, with their annotations and the truth of every element in every other genome.

Each region of the real FASTA is one chromosome, simulated alone: seqids c01, c02, ... in the order given, the
backbone its region's sequence and its elements those lying wholly inside it (donors.backbone_elements). Per
chromosome c: popmodel.ancestry (the genealogy), popevents.draw_events (insertions, losses, SNPs, small indels; its
`chrom` numbers the iids, "i<c:02d><n:06d>") and popbuild.build_chrom (every sample's sequence, Copy records and
lift; its `chrom` is the index c - 1 into PopGenome.seqs, which Copy.chrom holds). The parent then joins each
sample's chromosomes into one PopGenome (sample order: popmodel.sample_labels, which is the sheet order) and adds
assembly gaps (popbuild.add_gaps), annotations (poptruth.annotate), truth (poptruth.pair_truth), lifted tracts
(poptruth.tract_bed) and the insertion table (poptruth.insertion_rows).

Random streams (numpy seed sequences, one per purpose):
  donors             [seed, 0]       the donor sample (donors.pick_donors)
  chromosome c       [seed, c]       c = 1, 2, ...; its first draw is msprime's random_seed for ancestry; then
                                     draw_events (which draws its own msprime seeds first) and build_chrom read it
  gaps of genome g   [seed, 100, g]  g its sheet index, 0-based
  annotation of g    [seed, 200, g]
The streams stay apart only while c < 100 (MAX_CHROMS) and the seed is one 32-bit word (MAX_SEED): numpy's
SeedSequence zero-pads short entropy, so [seed, 100] would equal [seed, 100, 0], and a seed of 2^32 or more becomes
two words, so [x + 100 * 2^32, c] would equal [x, 100, c].
Results do not depend on `threads`: a chromosome is one task, in a worker of a fork pool of min(threads,
chromosomes) processes (inline for one), and builds its stream from [seed, c]; no random state crosses a process.

Memory: a worker returns its chromosome's per-sample sequences, copies and lifts, its kept insertions, event
counts and tracts. The parent keeps one copy of each genome's sequence (its chromosomes' bytes, never joined); gaps
are added one genome at a time and a FASTA is written in 1-Mb blocks.

Event counts (params.json "events"; per chromosome, per genome for gaps, and totals):
  requested   what the processes drew: insertions (primary and nested, before the drop rules), losses, SNP sites,
              small indels; gaps round(gap_per_mb x genome length / 1 Mb)
  realized    what the genomes hold: kept insertions (nested among them), kept losses by kind, SNP sites at ACGT
              bases, small indels at least one genome holds (build_chrom's `held`), gap runs written
  dropped     draw_events' counters (ins_n, ins_end, snp_n, ins_after_loss, nested_after_loss, loss_after_loss),
              build_chrom's small indels in no genome (indel_in_no_genome: a deletion an overlapping removal above
              it pre-empts, inserted bases every carrier loses with both junction bases; popbuild's small-indel
              rule) and add_gaps' near runs that fell off a chromosome (gap_off_chromosome); requested = realized +
              dropped
"""
from __future__ import annotations

import json
import logging
import math
import numbers
import os
from collections import Counter
from dataclasses import asdict, dataclass, field
from itertools import pairwise
from pathlib import Path

import numpy as np

from ltrtrace import __version__
from ltrtrace.base.checkpoint import sha256_file
from ltrtrace.base.constants import (
    SIM_DROPOUT,
    SIM_FAI_CACHE,
    SIM_GAP_PER_MB,
    SIM_JITTER,
    SIM_MAX_CHROMS,
    SIM_MAX_SEED,
    WINDOW_PAD,
)
from ltrtrace.base.log import UsageError
from ltrtrace.base.runtime import PoolError, pool_map
from ltrtrace.io.annotation import annotation_files, load_annotation
from ltrtrace.io.annotation.gff3 import write_gff3
from ltrtrace.io.elements import Element
from ltrtrace.io.seqio import Fasta
from ltrtrace.sim.pop.build import PopGenome, add_gaps, build_chrom
from ltrtrace.sim.pop.donors import Donor, backbone_elements, pick_donors
from ltrtrace.sim.pop.events import KINDS, Ins, draw_events
from ltrtrace.sim.pop.genome import SimElement
from ltrtrace.sim.pop.model import PULSES, Demog, Rates, _msprime, ancestry, sample_labels, tracts
from ltrtrace.sim.pop.truth import INSERTION_COLUMNS, alive, annotate, insertion_rows, pair_truth, tract_bed
from ltrtrace.sim.score import TRUTH_COLUMNS

log = logging.getLogger(__name__)

MAX_CHROMS = SIM_MAX_CHROMS
MAX_SEED = SIM_MAX_SEED
_GAPS, _ANNOT = 100, 200  # the gap and annotation streams' second word
assert MAX_CHROMS < _GAPS < _ANNOT, "chromosome streams [seed, c] must stay apart from [seed, 100, g]"
DROPPED = ("ins_n", "ins_end", "snp_n", "ins_after_loss", "nested_after_loss", "loss_after_loss")
INDEL_DROPPED = "indel_in_no_genome"  # build_chrom's small indels no genome holds
CACHE = SIM_FAI_CACHE
SIMPLIFICATIONS = (  # every approximation the popmodel, popevents, popbuild and poptruth docstrings name
    ("each region recombines at r_low between its fractions low[0] and low[1] (by default its middle 40%) and at r "
     "elsewhere: a per-chromosome rule applied to every region"),
    "replacement lengths are exponential with mean replace_mean (an assumed shape; the mean is measured)",
    "replacing sequence and the solo LTR of a simulated insertion receive no insertions",
    "an event is discarded when the lineage of any one of its carriers lacks its sequence (any-carrier thinning)",
    "a backbone element can take later losses until older losses have removed its whole span",
    "replacing sequence and inserted indel bases take no substitutions",
    "a primary insertion's TSD bases at insertion time assume one mutation per SNP site",
    ("a small deletion that overlaps a removal with the same carriers is taken as younger than the removal and is "
     "not made (indel times are not kept; a few bp per event)"),
    ("a backbone element's annotated K is read as its K at the root, so each sample's copy exceeds it by about "
     "2 x mu x T_root (T_root, the root's age, is at least t_o: 0.054 or more at the defaults)"),
    ("truth: a copy Y holds without an image (deleted), or an absent one whose empty-site window is broken, has no "
     "site in Y, so it is never G"),
    "truth: O (another insertion at the target in Y) does not ask whether X carries that insertion too",
    ("truth: a tract end inside a stretch the sample deleted lifts to the next kept base, so a payload anchored "
     "there falls on the wrong side of the end"),
)
_FA_BLOCK = 80 * 12_500  # FASTA written in 1-Mb blocks of 80-column lines


def _fraction(name: str, v, below_one: bool) -> float:
    if isinstance(v, bool) or not isinstance(v, numbers.Real) or not math.isfinite(v) or not 0 <= v <= 1 \
            or (below_one and v == 1):
        raise ValueError(f"{name} must be in [0, {'1)' if below_one else '1]'}: {v!r}")
    return float(v)


@dataclass(frozen=True)
class PopParams:
    """The population, its rates, and how the genomes are annotated and gapped. Checked on construction: a bad
    value, or dropout_genome naming the outgroup or a genome the demography does not sample, raises ValueError."""

    demog: Demog
    rates: Rates
    dropout: float = SIM_DROPOUT  # annotation noise defaults
    jitter: float = SIM_JITTER
    dropout_genome: tuple[tuple[str, float], ...] = ()  # (label, dropout) of genomes annotated at their own rate
    gap_per_mb: float = SIM_GAP_PER_MB

    def __post_init__(self):
        if not isinstance(self.demog, Demog) or not isinstance(self.rates, Rates):
            raise ValueError("PopParams needs a Demog and a Rates")  # noqa: TRY004 (bad input)
        _fraction("dropout", self.dropout, below_one=True)
        _fraction("jitter", self.jitter, below_one=False)
        if isinstance(self.gap_per_mb, bool) or not isinstance(self.gap_per_mb, numbers.Real) \
                or not math.isfinite(self.gap_per_mb) or self.gap_per_mb < 0:
            raise ValueError(f"gap_per_mb must be a finite number >= 0: {self.gap_per_mb!r}")
        groups = dict(sample_labels(self.demog))
        annotated = [lab for lab, grp in groups.items() if grp != "O"]
        seen: set[str] = set()
        pairs = []
        for item in self.dropout_genome:
            label, x = item
            if label not in groups:
                raise ValueError(f"dropout_genome: unknown genome {label!r}; the annotated genomes are "
                                 f"{', '.join(annotated)}")
            if groups[label] == "O":
                raise ValueError(f"dropout_genome: {label} is the outgroup, which has no annotation")
            if label in seen:
                raise ValueError(f"dropout_genome: {label} given twice")
            seen.add(label)
            pairs.append((label, _fraction(f"dropout_genome {label}", x, below_one=True)))
        object.__setattr__(self, "dropout_genome", tuple(pairs))

    def dropout_of(self, label: str) -> float:
        return dict(self.dropout_genome).get(label, self.dropout)


@dataclass
class PopResult:
    params: PopParams
    seed: int
    regions: list[tuple[str, int, int]]  # 1-based inclusive, chromosome c is regions[c - 1]
    genomes: list[PopGenome]  # sheet order
    tracts: dict[str, list[tuple[str, int, int]]]  # label -> lifted tracts (recipient-taxon samples)
    truth: dict[tuple[str, str], list[dict]]  # (X, Y) -> TRUTH_COLUMNS rows
    annotations: dict[str, list[Element]]  # ingroup labels
    insertions: list[dict]
    meta: dict = field(default_factory=dict)  # inputs, versions, trees, events, truth counts, fractions,
    #                                             simplifications: params.json's entries after seed, params, regions


@dataclass
class _ChromOut:
    """What a worker returns for chromosome c."""

    c: int
    samples: list  # per sample (sequence, copies, lift), ts.samples() order
    ins: list[Ins]
    requested: dict[str, int]
    realized: dict[str, int]
    dropped: dict[str, int]
    tracts: dict[int, list[tuple[int, int]]]
    trees: int


class RegionError(UsageError, ValueError):
    """A backbone region (--region) the simulation cannot use: malformed, overlapping another, off its FASTA
    sequence or too short (exit 2)."""


def _region_text(r: tuple[str, int, int]) -> str:
    return f"{r[0]}:{r[1]}-{r[2]}"


def check_regions(regions: list[tuple[str, int, int]], lengths: dict[str, int] | None = None,
                  fasta: Path | None = None) -> None:
    """RegionError, naming --region, unless there are 1 to MAX_CHROMS regions, each (seqid, start, end) with 1 <= start
    <= end (and, given the FASTA's lengths, on a sequence of the FASTA `fasta` and inside it), no two overlapping."""
    if not 1 <= len(regions) <= MAX_CHROMS:
        raise RegionError(f"--region given {len(regions)} times: need 1 to {MAX_CHROMS} regions (one chromosome each)")
    for r in regions:
        seqid, start, end = r
        if isinstance(start, bool) or isinstance(end, bool) or not isinstance(start, numbers.Integral) \
                or not isinstance(end, numbers.Integral) or not 1 <= start <= end:
            raise RegionError(f"--region {r!r}: expected (seqid, start, end), 1 <= start <= end (1-based inclusive)")
        if lengths is not None:
            if seqid not in lengths:
                raise RegionError(f"--region {_region_text(r)}: sequence {seqid!r} is not in "
                                  f"{fasta if fasta is not None else 'the FASTA'}")
            if end > lengths[seqid]:
                raise RegionError(f"--region {_region_text(r)} is outside {seqid} (1-{lengths[seqid]})")
    ordered = sorted(regions, key=lambda r: (r[0], r[1], r[2]))
    for a, b in pairwise(ordered):
        if a[0] == b[0] and b[1] <= a[2]:
            raise RegionError(f"--region {_region_text(a)} and {_region_text(b)} overlap; each region is one "
                              "chromosome and must be disjoint from the others")


def _versions() -> dict[str, str]:
    """ltrtrace's, msprime's and tskit's versions; RuntimeError(MSPRIME_HELP) without msprime."""
    msprime = _msprime()
    import tskit  # msprime depends on it

    return {"ltrtrace": __version__, "msprime": msprime.__version__, "tskit": tskit.__version__}


# ---------------------------------------------------------------- one chromosome (a worker)

def _chrom(c: int, backbone: bytes, elements: list[SimElement], donors: list[Donor], demog: Demog, rates: Rates,
           seed: int) -> _ChromOut:
    """Chromosome c (1-based): its genealogy, events and built samples, from the stream [seed, c]."""
    rng = np.random.default_rng([seed, c])
    ts = ancestry(demog, rates, len(backbone), int(rng.integers(1, 2**32 - 1)))
    ev = draw_events(ts, backbone, elements, donors, rates, rng, chrom=c)
    held: set[int] = set()
    samples = build_chrom(c - 1, backbone, elements, ts, ev, donors, rates, rng, held=held)
    dropped = {k: int(ev.dropped.get(k, 0)) for k in DROPPED}
    dropped[INDEL_DROPPED] = len(ev.indel_pos) - len(held)
    kinds = Counter(x.kind for x in ev.losses)
    n_ins = sum(bool(ev.indel_ins[k]) for k in held)
    realized = {"insertions": len(ev.ins), "nested": sum(x.depth > 0 for x in ev.ins), "losses": len(ev.losses),
                **{k: kinds[k] for k in KINDS}, "snps": len(ev.snp_pos), "indels": len(held),
                "indel_insertions": n_ins, "indel_deletions": len(held) - n_ins}
    requested = {"insertions": len(ev.ins) + dropped["ins_n"] + dropped["ins_end"] + dropped["ins_after_loss"]
                 + dropped["nested_after_loss"], "losses": len(ev.losses) + dropped["loss_after_loss"],
                 "snps": len(ev.snp_pos) + dropped["snp_n"], "indels": len(ev.indel_pos)}
    return _ChromOut(c, samples, ev.ins, requested, realized, dropped, tracts(ts, demog), ts.num_trees)


def _run_chroms(jobs: list[tuple[int, bytes, list[SimElement]]], donors: list[Donor], demog: Demog, rates: Rates,
                seed: int, threads: int, names: dict[int, str]) -> list[_ChromOut]:
    """Every chromosome's _ChromOut in chromosome order; a failure raises RuntimeError naming the chromosome."""
    def one(job: tuple[int, bytes, list[SimElement]]) -> _ChromOut:
        c, backbone, elements = job
        try:
            r = _chrom(c, backbone, elements, donors, demog, rates, seed)
        except Exception as exc:
            raise RuntimeError(f"{names[c]} failed: {type(exc).__name__}: {exc}") from exc
        log.info("%s done: %d trees, %d insertions, %d losses, %d SNPs", names[c], r.trees, r.realized["insertions"],
                 r.realized["losses"], r.realized["snps"])
        return r

    try:
        return pool_map(one, jobs, min(threads, len(jobs)), warm=False)  # no numba here: all chromosomes at once
    except PoolError as exc:
        raise RuntimeError(f"a worker process died or could not start (killed by a signal or out of memory, or "
                           f"crashed): {exc}") from exc


# ---------------------------------------------------------------- the population


def simulate_pop(regions: list[tuple[str, int, int]], fasta: Path, annotation: str, params: PopParams, seed: int,
                 threads: int, *, cache_dir: Path | None = None) -> PopResult:
    """Simulate the population of `params` on the backbone `regions` of `fasta` (1-based inclusive; chromosome c is
    regions[c - 1]), the annotation giving the backbone elements and the donors (module docstring). Fails fast,
    before anything is simulated: without msprime (RuntimeError(MSPRIME_HELP)); on bad or overlapping regions, or a
    region shorter than 2 x WINDOW_PAD plus the longest donor (RegionError, a ValueError); with no donor while lam >
    0, or backbone element names that repeat across chromosomes (ValueError). `cache_dir` holds the FASTA's index if
    its directory is read-only and it has none."""
    versions = _versions()
    if isinstance(seed, bool) or not isinstance(seed, numbers.Integral) or not 0 <= seed <= MAX_SEED:
        raise ValueError(f"seed must be an integer in [0, {MAX_SEED}]: {seed!r}")
    if isinstance(threads, bool) or not isinstance(threads, numbers.Integral) or threads < 1:
        raise ValueError(f"threads must be an integer >= 1: {threads!r}")
    regions = [tuple(r) for r in regions]
    check_regions(regions)
    demog, rates = params.demog, params.rates
    with Fasta(fasta, cache_dir=cache_dir) as fa:
        check_regions(regions, fa.lengths(), fasta)
        log.debug("loading the annotation %s", annotation)
        elements = list(load_annotation(annotation, "ref", fa))  # rows built once: the scans below pass several times
        donors = pick_donors(elements, fa, seed)
        if rates.lam > 0 and not donors:
            raise ValueError("no donor elements to insert (lam > 0): the annotation has no top-level, childless "
                             "element (high-confidence where it rates confidence) without an N")
        longest = max((len(d.seq) for d in donors), default=0)
        need = 2 * WINDOW_PAD + longest
        for r in regions:
            if r[2] - r[1] + 1 < need:
                raise RegionError(f"--region {_region_text(r)} is {r[2] - r[1] + 1} bp, shorter "
                                  f"than 2 x {WINDOW_PAD} bp plus the longest donor ({longest} bp) = "
                                  f"{need} bp")
        jobs, where = [], {}
        for c, (seqid, start, end) in enumerate(regions, 1):
            els = backbone_elements(elements, seqid, start, end)
            for e in els:
                if e.eid in where:
                    raise ValueError(f"backbone element {e.eid} is on chromosomes c{where[e.eid]:02d} and c{c:02d}: "
                                     "backbone names must be unique (they are the truth tables' names)")
                where[e.eid] = c
            jobs.append((c, fa.fetch(seqid, start - 1, end).encode(), els))
    n_annotation = len(elements)
    del elements, where
    names = {c: f"c{c:02d} ({_region_text(r)})" for c, r in enumerate(regions, 1)}
    labels = sample_labels(demog)
    bp = sum(len(b) for _, b, _ in jobs)
    log.info("ltrtrace sim pop: %d genomes (%d per taxon), %d chromosomes (%d bp), %d backbone elements, %d donors, "
             "seed %d, pulse %s f %g at %g generations, %d threads", len(labels), demog.per_taxon, len(jobs), bp,
             sum(len(e) for _, _, e in jobs), len(donors), seed, demog.pulse, demog.f, demog.ta, threads)
    meta = {"fasta": {"path": str(Path(fasta).resolve()), "sha256": sha256_file(Path(fasta))},
            "annotation": {"path": os.path.abspath(annotation),
                           "files": [{"path": str(p.resolve()), "sha256": sha256_file(p)}
                                     for p in annotation_files(annotation)]},
            "versions": versions,
            "backbone": {"annotation_elements": n_annotation, "donors": len(donors), "longest_donor": longest,
                         "chromosomes": {f"c{c:02d}": {"region": _region_text(regions[c - 1]), "length": len(b),
                                                       "elements": len(e)} for c, b, e in jobs}}}
    out = _run_chroms(jobs, donors, demog, rates, seed, threads, names)
    del jobs

    genomes = _join(out, labels)
    requested = {"per_chromosome": {f"c{r.c:02d}": r.requested for r in out}, "per_genome": {}}
    realized = {"per_chromosome": {f"c{r.c:02d}": r.realized for r in out}, "per_genome": {}}
    dropped = {"per_chromosome": {f"c{r.c:02d}": r.dropped for r in out}, "per_genome": {}}
    for i, g in enumerate(genomes):
        off = add_gaps([g], params.gap_per_mb, [np.random.default_rng([seed, _GAPS, i])])[0]
        requested["per_genome"][g.label] = {"gaps": len(g.gaps) + off}
        realized["per_genome"][g.label] = {"gaps": len(g.gaps)}
        dropped["per_genome"][g.label] = {"gap_off_chromosome": off}
    for d in (requested, realized, dropped):
        total: Counter = Counter()
        for part in ("per_chromosome", "per_genome"):
            for v in d[part].values():
                total.update(v)
        d["total"] = dict(sorted(total.items()))
    dropped["total"] = {k: dropped["total"].get(k, 0) for k in (*DROPPED, INDEL_DROPPED, "gap_off_chromosome")}

    annotations = {g.label: annotate(g, params.dropout_of(g.label), params.jitter,
                                     np.random.default_rng([seed, _ANNOT, i]))
                   for i, g in enumerate(genomes) if g.group != "O"}
    truth = {(x.label, y.label): pair_truth(x, y) for x in genomes if x.label in annotations
             for y in genomes if y is not x}
    recipient = PULSES[demog.pulse][1][0]  # the recipient population's taxon: A or C
    lifted = {g.label: tract_bed(g, {r.c - 1: r.tracts.get(s, []) for r in out}) for s, g in enumerate(genomes)}
    fractions = []
    for s, g in enumerate(genomes):
        tract_bp = sum(hi - lo for r in out for lo, hi in r.tracts.get(s, []))
        lifted_bp = sum(b - a for _, a, b in lifted[g.label])
        length = sum(len(x) for x in g.seqs)
        fractions.append({"label": g.label, "group": g.group, "backbone_bp": bp, "tract_bp": tract_bp,
                          "fraction": tract_bp / bp, "genome_bp": length, "lifted_bp": lifted_bp,
                          "lifted_fraction": lifted_bp / length if length else 0.0})
    insertions = insertion_rows(genomes, {r.c - 1: r.ins for r in out})
    meta.update({
        "trees": {f"c{r.c:02d}": r.trees for r in out},
        "events": {"requested": requested, "realized": realized, "dropped": dropped},
        "genomes": {g.label: {"group": g.group, "length": sum(len(x) for x in g.seqs),
                              "alive": sum(alive(c) for c in g.copies.values()),
                              "annotated": len(annotations.get(g.label, ())),
                              "dropout": params.dropout_of(g.label) if g.label in annotations else None,
                              "gaps": len(g.gaps)} for g in genomes},
        "truth": {x: {y: dict(sorted(Counter(row["state"] for row in truth[(x, y)]).items()))
                      for (x2, y) in truth if x2 == x} for x in annotations},
        "fractions": fractions,
        "simplifications": list(SIMPLIFICATIONS),
    })
    log.info("truth: %d ordered pairs; %d insertions; tracts in %s", len(truth), len(insertions),
             ", ".join(f"{lab} {sum(b - a for _, a, b in v)} bp" for lab, v in lifted.items() if v) or "no genome")
    return PopResult(params, int(seed), regions, genomes,
                     {g.label: lifted[g.label] for g in genomes if g.group == recipient}, truth, annotations,
                     insertions, meta)


def _join(out: list[_ChromOut], labels: list[tuple[str, str]]) -> list[PopGenome]:
    """Each sample's chromosomes as one PopGenome, in sample order. The contract checked here: chromosome c's
    copies all have Copy.chrom == c - 1 (its index in seqs and seqids), and no element id is on two chromosomes."""
    seqids = [f"c{r.c:02d}" for r in out]
    if [r.c for r in out] != list(range(1, len(out) + 1)):
        raise RuntimeError(f"chromosomes out of order: {[r.c for r in out]}")
    genomes = []
    for s, (label, group) in enumerate(labels):
        seqs, lifts, copies = [], [], {}
        for k, r in enumerate(out):
            if len(r.samples) != len(labels):
                raise RuntimeError(f"{seqids[k]}: {len(r.samples)} samples built, {len(labels)} expected")
            seq, cps, lift = r.samples[s]
            r.samples[s] = None  # the genome holds it now
            bad = next((e for e, cp in cps.items() if cp.chrom != k), None)
            if bad is not None:
                raise RuntimeError(f"{seqids[k]}: copy {bad} has chrom {cps[bad].chrom}, not its index {k}")
            if not copies.keys().isdisjoint(cps):
                raise RuntimeError(f"{seqids[k]}: element id(s) {sorted(copies.keys() & cps.keys())[:3]} already on "
                                   "another chromosome")
            seqs.append(seq)
            lifts.append(lift)
            copies.update(cps)
        genomes.append(PopGenome(label, group, list(seqids), seqs, copies, [], lifts))
    return genomes


# ---------------------------------------------------------------- outputs


def _write_fasta(path: Path, g: PopGenome) -> None:
    with open(path, "wb") as fh:
        for seqid, s in zip(g.seqids, g.seqs):
            fh.write(f">{seqid}\n".encode())
            for i in range(0, len(s), _FA_BLOCK):
                block = b"\n".join(s[j:j + 80] for j in range(i, min(i + _FA_BLOCK, len(s)), 80))
                fh.write(block + b"\n")


def _write_tsv(path: Path, columns: list[str] | tuple[str, ...], rows: list[dict]) -> None:
    with open(path, "w") as fh:
        fh.write("\t".join(columns) + "\n")
        fh.writelines("\t".join(str(r[c]) for c in columns) + "\n" for r in rows)


def check_outdir(outdir: Path) -> None:
    """ValueError unless `outdir` is missing or a directory holding nothing but its own cache (CACHE, a
    directory): a rerun into an old outdir would leave its stale truth files beside the new ones."""
    out = Path(outdir)
    if not out.exists():
        return
    if not out.is_dir():
        raise ValueError(f"{out} exists and is not a directory")
    other = sorted(p.name for p in out.iterdir() if not (p.name == CACHE and p.is_dir()))
    if other:
        more = f" and {len(other) - 3} more" if len(other) > 3 else ""
        raise ValueError(f"{out} is not empty ({', '.join(other[:3])}{more}): give a new or empty directory, "
                         f"so no stale file of an earlier run stays beside the new outputs")


def write_pop(outdir: Path, result: PopResult) -> None:
    """<label>.fa per genome; <label>.gff3 per annotated genome; samples.tsv (label, fasta, annotation, group,
    role; paths relative to it, the outgroup unannotated with role outgroup); truth/<X>__<Y>.tsv per annotated X and
    other Y (start, end 1-based inclusive in X); truth/insertions.tsv (pos 0-based, as in BED: the copy goes before
    base pos of its chromosome's backbone, 0 at its region's start, or of its host's donor body when host is an iid,
    so pos is a coordinate of no output genome); truth/tracts.<label>.bed per recipient-taxon genome (0-based
    half-open); truth/fractions.tsv (introgressed share per genome: backbone bp in tracts, and the lifted tracts' bp
    in the genome); params.json (seed, params, regions and the result's meta). An outdir holding anything but its own
    cache is refused first (check_outdir: ValueError)."""
    out = Path(outdir)
    check_outdir(out)
    (out / "truth").mkdir(parents=True, exist_ok=True)
    for g in result.genomes:
        _write_fasta(out / f"{g.label}.fa", g)
        if g.label in result.annotations:
            write_gff3(result.annotations[g.label], out / f"{g.label}.gff3")
    lines = ["label\tfasta\tannotation\tgroup\trole\n"]
    for g in result.genomes:
        ann = f"{g.label}.gff3" if g.label in result.annotations else "."
        lines.append(f"{g.label}\t{g.label}.fa\t{ann}\t{g.group}\t{'outgroup' if g.group == 'O' else '.'}\n")
    (out / "samples.tsv").write_text("".join(lines))
    for (x, y), rows in result.truth.items():
        _write_tsv(out / "truth" / f"{x}__{y}.tsv", TRUTH_COLUMNS, rows)
    labels = [g.label for g in result.genomes]
    _write_tsv(out / "truth" / "insertions.tsv", [*INSERTION_COLUMNS, *(f"fate_{lab}" for lab in labels)],
               result.insertions)
    for label, recs in result.tracts.items():
        (out / "truth" / f"tracts.{label}.bed").write_text("".join(f"{s}\t{a}\t{b}\n" for s, a, b in recs))
    fractions = [{**f, "fraction": f"{f['fraction']:.6f}", "lifted_fraction": f"{f['lifted_fraction']:.6f}"}
                 for f in result.meta.get("fractions", ())]
    _write_tsv(out / "truth" / "fractions.tsv", ("label", "group", "backbone_bp", "tract_bp", "fraction",
                                                  "genome_bp", "lifted_bp", "lifted_fraction"), fractions)
    info = {"seed": result.seed, "params": asdict(result.params), "regions": [list(r) for r in result.regions],
            **result.meta}
    (out / "params.json").write_text(json.dumps(info, indent=2, sort_keys=True) + "\n")
