"""Toy run directories for the analysis tests, written through the engine's own writers and merged by its own
merge: work/elements/<label>.tsv (`write_elements`), the pair files pairs/<A>__<B>.calls.tsv.gz
(CALL_COLUMNS, `tsv_field` values), loci.tsv (`build_loci` over those pair files, then `write_loci_tsv`),
qc.tsv (`write_qc`), the samples sheet and N-only FASTAs.

The toy writes only what the engine writes (classify.py, genotype.py, loci.py):
  tiers     a call's tier is placement + its state's evidence: P, S (with solo_tsd) and E (with tsd_match)
            high; T and O medium (they never have the evidence a high tier needs); D medium too (the
            engine writes D high when its junction is refined; medium is engine-possible and the toy
            never refines); G and U ".".
  sites     P and T at the element; S at its solo LTR (pos .. pos + ltr - 1); E at its TSD copy pos - 5 ..
            pos - 1; D at the junction pos - 1 .. pos + b_ins (two flank bases when b_ins is 0); O at the
            sequence B holds in the element's
            place, pos .. pos + b_ins - 1; G and U none.
  strands   b_strand is the window's orientation, "+" or "-" for every call with a site: "-" where exactly
            one genome of the pair is mirrored by the inversion, flipped again where both genomes' given
            strands are known and differ. A call's element strand is the caller's turned by it, so it is
            unknown only where the caller's own strand is: "." is an annotated member's strand (a BED
            annotation), never a non-member's.
  inv       on every call from a reversed window, G included; U is no_anchor with no inv (an
            unresolved window has no strand) but keeps its window id, as the engine gives every top-level
            element a window.
  loci.tsv  loci.py's merge of the pair files: each state, tier, site, strand and flag is the deciding
            call's, a member is P/annot. Lids are in list order (the engine orders them by position); each
            locus's states must be what its calls merge to, and its members must be linked by P/S/T calls
            (the engine splits a locus whose members never find each other).
Not modeled: rule 0 (a host's E, D, O, G or U passed to what is nested in it, with the host's site and tier).
A nested element's calls are written as given, so a test can nest a P inside a D host (test_load_roundtrip
does); inherited calls (flags host_empty, inherit) carry no evidence flag, as in the engine.

Exact runs come from hand-written loci (`loc`, `write_run`); simulated runs (`simulate_run`) draw insertion
patterns, members and missing calls from a seeded generator: the engine's behavior in miniature.
Coordinates are 1-based inclusive. A locus sits at the same place (seqid, pos) in every genome unless `at`
places it elsewhere in one genome (a translocation: the engine anchors a window wherever the flanks lie) or
`inversion` mirrors it in one.
"""
from __future__ import annotations

import gzip
import hashlib
import os
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ltrtrace.analyses.geneflow.catalog import STATE_CODE, Catalog, load_catalog
from ltrtrace.analyses.geneflow.route import ROUTE_PROPS, block_bp
from ltrtrace.calls.loci import Locus, build_loci
from ltrtrace.calls.records import CALL_COLUMNS, tsv_field
from ltrtrace.calls.tables import write_loci_tsv, write_qc
from ltrtrace.io.elements import Element, make_eid, write_elements
from ltrtrace.io.samples import PAIR_SEP, Sample, read_samples

STATES = tuple(STATE_CODE)  # the states the assembly route writes (no A: an inferred absence is the library route's)

# label, group, role, annotated
GENOMES = [("R", "taxR", None, True), ("S", "taxS", None, True), ("D", "taxD", None, True),
           ("O", "taxO", "outgroup", False)]
TRACTS = [("chr1", 100_001, 400_000), ("chr2", 500_001, 800_000), ("chr3", 200_001, 500_000)]  # R intervals

TSD = "ACGTA"  # every element's TSD
E_SITE = 5  # an empty site is the TSD copy pos - 5 .. pos - 1
EVIDENCE = {"E": (0, 0, 0), "D": (300, 0, 0), "O": (2, -1, 500)}  # a call's default (dl, dr, b_ins)
TIER = {"P": "high", "S": "high", "T": "medium", "E": "high", "D": "medium", "O": "medium",  # a call's tier
        "G": ".", "U": "."}
_PROOF = {"S": "solo_tsd", "E": "tsd_match"}  # the evidence that puts these states at tier high
_INHERITED = ("host_empty", "inherit")  # flags of a call passed from a host (rule 0): no evidence of its own
_DERIVED = ("inv", "conflict", "no_call", "annot_conflict")  # flags the engine sets itself, never given
_FLIP = {"+": "-", "-": "+", ".": "."}
_LINE = 80  # FASTA line width


@dataclass
class ToyLocus:
    seqid: str
    pos: int
    states: dict[str, str]                                             # loci.tsv (merged) states
    members: frozenset[str] = frozenset()
    k: float | None = 0.01
    length: int = 150
    ltr: int = 40                                                      # LTRs = the first and last `ltr` bp
    family: str = "famA"
    superfamily: str = "Gypsy"
    parent: int | None = None                                          # index of the host locus in the list
    flags: dict[str, list[str]] = field(default_factory=dict)          # genome -> extra flags
    evidence: dict[str, tuple[int, int, int]] = field(default_factory=dict)   # genome -> (dl, dr, b_ins) in pair files
    calls: dict[tuple[str, str], str] = field(default_factory=dict)    # (A, B) -> state of A's call into B
    strand: dict[str, str] = field(default_factory=dict)               # genome -> strand ("+" by default)
    wid: int | None = None                                             # pair-file wid; None: the locus index
    at: dict[str, tuple[str, int]] = field(default_factory=dict)       # genome -> (seqid, pos) where it holds
                                                                       # the locus (default: seqid, pos)


@dataclass
class ToyRun:
    """A written run: its outdir, the samples sheet (path and parsed rows), the loci it holds (one per
    loci.tsv row, in order), and a simulated run's truth."""

    outdir: Path
    sheet: Path
    samples: list[Sample]
    loci: list[ToyLocus]
    truth: dict = field(default_factory=dict)

    def load(self) -> Catalog:
        return load_catalog(self.outdir, self.samples, self.sheet)


def loc(pos: int, spec: str, seqid: str = "chr1", **kw) -> ToyLocus:
    """A locus at seqid:pos. `spec` is space-separated `G:X`, X a state letter with `*` after it marking a
    member (always P): "R:P* S:E D:P* O:E". Other keywords go to ToyLocus."""
    states, members = {}, set()
    for word in spec.split():
        label, sep, x = word.partition(":")
        member = x.endswith("*")
        x = x[:-1] if member else x
        if not sep or not label or x not in STATES or label in states:
            raise ValueError(f"loc: bad word {word!r} in {spec!r}: G:X, X one of {''.join(STATES)}, "
                             "* for a member")
        if member and x != "P":
            raise ValueError(f"loc: member {label} must be P (a member is P/annot in loci.tsv): {spec!r}")
        states[label] = x
        if member:
            members.add(label)
    return ToyLocus(seqid=seqid, pos=pos, states=states, members=frozenset(members), **kw)


def block_size(cat: Catalog) -> int:
    """The assembly route's jackknife block of a toy run (route.block_bp): 1/100 of its smallest genome, at least
    10 kb (the default toy's four 1-Mb sequences: 40 kb)."""
    return block_bp(cat, ROUTE_PROPS["assembly"])


def dmat(pairs: dict[tuple[str, str], float], default: float, labels) -> np.ndarray:
    """A symmetric (G, G) distance matrix over `labels`: 0 on the diagonal, pairs[(a, b)] where given, else
    `default`."""
    labels = list(labels)
    d = np.full((len(labels), len(labels)), float(default))
    np.fill_diagonal(d, 0.0)
    for (a, b), v in pairs.items():
        i, j = labels.index(a), labels.index(b)
        d[i, j] = d[j, i] = v
    return d


def read_tsv(path: Path) -> list[dict[str, str]]:
    """The rows of a TSV as dicts keyed by its header (a leading # on the header is dropped); blank lines
    skipped."""
    with open(path, encoding="utf-8", newline="") as fh:
        lines = [line for line in fh.read().splitlines() if line]
    if not lines:
        return []
    header = lines[0].removeprefix("#").split("\t")
    return [dict(zip(header, line.split("\t"), strict=True)) for line in lines[1:]]


def write_run(tmp: Path, loci: list[ToyLocus], genomes=GENOMES, seqids=("chr1", "chr2", "chr3", "chr4"),
              seq_len: int = 1_000_000, seq_lens: dict[str, int] | None = None, fasta_dir: Path | None = None,
              fastas: dict[str, Path] | None = None, inversion: tuple[str, str, int, int] | None = None,
              prefix: bool = False) -> ToyRun:
    """Write a finished run of `loci` (lids L0000001... in list order) into <tmp>/run, its sheet
    <tmp>/samples.tsv. A ValueError names any locus the engine could not have written (module docstring).

    genomes: (label, group, role, annotated) in sheet order. seqids, seq_len, seq_lens (seqid -> length):
    every genome's sequences, N only, written once per layout as <fasta_dir>/<layout>/<label>.fa (fasta_dir:
    default <tmp>/fasta); fastas: label -> a FASTA with real sequence instead. prefix: name each genome's
    seqids <label>_<seqid>. A locus's `at` places it at another (seqid, pos) in the genomes it names.
    inversion (label, seqid, a, b): that genome's sites of the loci it holds starting in [a, b] (its own
    placement) are mirrored (start a + b - end, end a + b - start) and its element strands there reversed.

    Each member has an element in work/elements/<label>.tsv (LTRs its first and last `ltr` bp, TSD ACGTA,
    parent the host locus's element there). Each annotated A's pair file into every other B holds one call
    per A element, in the engine's row order: by default P into a member, else B's state in `states`
    (`calls` overrides), with the tier, site, b_strand and flags of the module docstring, the given flags
    of B after the evidence flag; reason no_anchor for U; dl dr b_ins from `evidence`, else EVIDENCE by
    state; tsd_a ACGTA, other evidence ".". loci.tsv is build_loci's merge of these files. qc.tsv:
    annot_sensitivity 0.9 per annotated genome and skipped_malformed_coordinates 2 for the first.
    """
    tmp = Path(tmp)
    labels = [g[0] for g in genomes]
    annotated = [g[0] for g in genomes if g[3]]
    lengths = {s: (seq_lens or {}).get(s, seq_len) for s in seqids}
    _check(loci, labels, set(annotated), lengths, inversion)
    names = {x: [f"{x}_{s}" if prefix else s for s in seqids] for x in labels}
    place = _Places(loci, seqids, names, lengths, inversion)
    elements = _elements(loci, annotated, place)
    calls = {(a, b): _calls(a, b, loci, elements[a], place) for a in annotated for b in labels if b != a}
    layout = hashlib.sha1(repr((tuple(seqids), tuple(lengths.values()), prefix)).encode()).hexdigest()[:12]
    fdir = Path(fasta_dir) if fasta_dir is not None else tmp / "fasta"
    paths = {x: Path(fastas[x]) if fastas and x in fastas
             else _n_fasta(fdir / layout / f"{x}.fa", names[x], list(lengths.values())) for x in labels}
    sheet = _write_sheet(tmp / "samples.tsv", genomes, paths)
    samples = read_samples(sheet)
    outdir = tmp / "run"
    for d in (outdir / "pairs", outdir / "work" / "elements"):
        d.mkdir(parents=True, exist_ok=True)
    for a in annotated:
        write_elements(list(elements[a].values()), outdir / "work" / "elements" / f"{a}.tsv")
    call_paths = {}
    for (a, b), lines in calls.items():
        call_paths[(a, b)] = path = outdir / "pairs" / f"{a}{PAIR_SEP}{b}.calls.tsv.gz"
        with open(path, "wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=1,
                                                    mtime=0) as gz:
            gz.write("".join(lines).encode())
    merged = _merge(loci, labels, elements, call_paths, samples, names)
    write_loci_tsv(merged, samples, outdir / "loci.tsv")
    qc = [(x, "annot_sensitivity", 0.9) for x in annotated]
    if annotated:
        qc.insert(1, (annotated[0], "skipped_malformed_coordinates", 2))
    write_qc(qc, outdir / "qc.tsv")
    return ToyRun(outdir, sheet, samples, list(loci))


# simulate_run: insertion patterns over (R, S, D) with their probabilities, the two K values (each drawn with
# probability 1/2) and whether O carries the insertion; "old" is present in all four genomes
_PATTERNS = (((1, 1, 1), 0.24, (0.02, 0.02), False),
             ((1, 1, 0), 0.15, (0.01, 0.02), False),
             ((1, 0, 1), 0.06, (0.004, 0.02), False),
             ((0, 1, 1), 0.06, (0.004, 0.02), False),
             ((1, 0, 0), 0.11, (0.0, 0.004), False),
             ((0, 1, 0), 0.11, (0.0, 0.004), False),
             ((0, 0, 1), 0.11, (0.0, 0.004), False),
             ((1, 1, 1), 0.16, (0.03, 0.03), True))
_D_PRIVATE = 6  # index of pattern 001 in _PATTERNS
_SIM_SEQIDS, _SIM_SPAN, _SIM_LENGTH = ("chr1", "chr2", "chr3", "chr4"), 998_000, 150
_TRUE = np.array(list("PST"))  # a present genome's true state, drawn with probabilities _TRUE_P
_TRUE_P = (0.9, 0.05, 0.05)


def simulate_run(tmp: Path, *, seed: int, tracts=(), inversion: tuple[str, str, int, int] | None = None,
                 p_missing: float | dict[str, float] = 0.15, p_member: float | dict[str, float] = 0.9,
                 n_per_seq: int = 1500, k_known: bool = True, prefix: bool = False,
                 fasta_dir: Path | None = None) -> ToyRun:
    """A simulated run of genomes GENOMES (R, S, D annotated; O the outgroup), written with write_run.

    n_per_seq loci per seqid on chr1-chr4 (1 Mb each) at 1001 + j * 998_000 // n_per_seq, 150 bp with 40-bp
    LTRs; the spacing must be >= 2 x the element length (n_per_seq <= 3000 keeps 332 bp), else ValueError.
    Each locus draws a pattern from _PATTERNS (its K one of two values, probability 1/2 each). Inside a tract
    (seqid, a, b), an R interval, R's presence := D's presence; a locus with no genome present is dropped. A
    present genome's true state is P, S or T (0.9, 0.05, 0.05), an absent genome's E. Each present annotated
    genome is a member with probability p_member[g], its true state then P; a locus with no member is
    dropped. Each member A's call into each other genome B (members included) is B's true state with
    probability 1 - p_missing[B], else U. Members are linked as the engine links them, by a resolved call of
    either into the other (a P call on its element): members never linked are separate loci, one per linked
    set (sheet order), as loci.py writes them. Each such locus merges like loci.py: a member is P/annot;
    another genome takes its members' resolved call (they agree), else U. k_known=False writes no K. A dict
    for p_missing or p_member sets the genomes it lists; the others keep the defaults (0.15, 0.9). The draws
    do not depend on tracts, inversion or prefix. truth: {"tracts": tracts, "n_rd_tract": the insertions
    written whose D-private presence the tract copied into R}.
    """
    spacing = _SIM_SPAN // n_per_seq
    if spacing < 2 * _SIM_LENGTH:
        raise ValueError(f"simulate_run: locus spacing {spacing} bp ({_SIM_SPAN} // {n_per_seq}) is below "
                         f"2 x the {_SIM_LENGTH}-bp element length; use n_per_seq <= 3000")
    labels = [g[0] for g in GENOMES]
    miss = np.array([_param(p_missing, x, 0.15) for x in labels])
    memb = np.array([_param(p_member, x, 0.9) for x in labels[:3]])
    rng = np.random.default_rng(seed)
    pos = np.tile(1001 + np.arange(n_per_seq) * _SIM_SPAN // n_per_seq, len(_SIM_SEQIDS))
    seqid = np.repeat(np.array(_SIM_SEQIDS), n_per_seq)
    n = len(pos)
    pat = rng.choice(len(_PATTERNS), size=n, p=[p for _, p, _, _ in _PATTERNS])
    k = np.array([kk for _, _, kk, _ in _PATTERNS])[pat, (rng.random(n) < 0.5).astype(int)]
    present = np.zeros((n, 4), dtype=bool)
    present[:, :3] = np.array([pt for pt, _, _, _ in _PATTERNS], dtype=bool)[pat]
    present[:, 3] = np.array([old for *_, old in _PATTERNS])[pat]
    inside = np.zeros(n, dtype=bool)
    for sid, a, b in tracts:
        inside |= (seqid == sid) & (pos >= a) & (pos <= b)
    present[inside, 0] = present[inside, 2]
    copied = inside & (pat == _D_PRIVATE)
    true = _TRUE[rng.choice(3, size=(n, 4), p=_TRUE_P)]
    true[~present] = "E"
    member = (rng.random((n, 3)) < memb) & present[:, :3]
    true[:, :3][member] = "P"
    resolved = rng.random((n, 3, 4)) >= miss  # [locus, A (R, S, D), B]: A's call into B resolves
    loci, n_rd = [], 0
    for i in np.flatnonzero(member.any(axis=1)):
        for linked in _linked([a for a in range(3) if member[i, a]], resolved[i]):
            calls = {(labels[a], labels[b]): str(true[i, b]) if resolved[i, a, b] else "U"
                     for a in linked for b in range(4) if b != a}
            states = {}
            for b, x in enumerate(labels):
                if b in linked:
                    states[x] = "P"
                else:
                    first = next((a for a in linked if resolved[i, a, b]), None)
                    states[x] = "U" if first is None else str(true[i, b])
            loci.append(ToyLocus(seqid=str(seqid[i]), pos=int(pos[i]), states=states,
                                 members=frozenset(labels[a] for a in linked),
                                 k=float(k[i]) if k_known else None, length=_SIM_LENGTH, ltr=40, calls=calls))
        n_rd += bool(copied[i])
    run = write_run(tmp, loci, inversion=inversion, prefix=prefix, fasta_dir=fasta_dir)
    run.truth = {"tracts": list(tracts), "n_rd_tract": n_rd}
    return run


# ---------------------------------------------------------------- helpers


def _param(p: float | dict[str, float], label: str, default: float) -> float:
    return p.get(label, default) if isinstance(p, dict) else p


def _linked(members: list[int], resolved: np.ndarray) -> list[list[int]]:
    """The sets of members (indexes into R, S, D) the engine links: two are linked when either's call into
    the other resolves (resolved[a, b]: a P call on the other's element). Sheet order inside and across."""
    sets: list[list[int]] = []
    for a in members:
        touch = [s for s in sets if any(resolved[a, b] or resolved[b, a] for b in s)]
        sets = [s for s in sets if s not in touch] + [sorted({a, *(b for s in touch for b in s)})]
    return sorted(sets)


def _check(loci: list[ToyLocus], labels: list[str], annotated: set[str], lengths: dict[str, int],
           inversion: tuple[str, str, int, int] | None) -> None:
    """Fail on a toy the engine could not have written: a locus without every genome's state or without a
    member, a member without annotation, a call from a non-member, a position outside its sequence, an
    unknown strand on a non-member, a flag only the engine sets, a bad parent or inversion. Sites are
    checked where they are placed (_Places.site)."""
    if len(set(labels)) != len(labels):
        raise ValueError(f"duplicate genome label in {labels}")
    if inversion is not None and (inversion[0] not in labels or inversion[1] not in lengths
                                  or not 1 <= inversion[2] <= inversion[3]):
        raise ValueError(f"inversion {inversion}: not (label, seqid, a, b) with 1 <= a <= b")
    for i, x in enumerate(loci):
        where = f"locus {i} ({x.seqid}:{x.pos})"
        if sorted(x.states) != sorted(labels):
            raise ValueError(f"{where}: states for {sorted(x.states)}, the genomes are {labels}")
        if not x.members:
            raise ValueError(f"{where}: no member; every locus is some annotated genome's element")
        if not x.members <= annotated:
            raise ValueError(f"{where}: member(s) {sorted(x.members - annotated)} without annotation")
        if x.seqid not in lengths or x.pos <= E_SITE or x.pos + x.length - 1 > lengths[x.seqid]:
            raise ValueError(f"{where}: not inside a sequence of {lengths}, empty site included")
        for g, (s, p) in x.at.items():
            if g not in labels or s not in lengths or p <= E_SITE or p + x.length - 1 > lengths[s]:
                raise ValueError(f"{where}: at[{g!r}] = {s}:{p}: not a genome, or not inside a sequence of "
                                 f"{lengths}, empty site included")
        if not 0 < 2 * x.ltr <= x.length:
            raise ValueError(f"{where}: LTRs of {x.ltr} bp do not fit a {x.length}-bp element")
        if x.parent is not None and not (0 <= x.parent < len(loci) and x.parent != i):
            raise ValueError(f"{where}: parent {x.parent} is not another locus of the list")
        for (a, b), s in x.calls.items():
            if a not in x.members or b not in labels or b == a or s not in STATES:
                raise ValueError(f"{where}: call {a}->{b} {s!r}: A must be a member, B another genome")
        for g in (*x.flags, *x.evidence, *x.strand):
            if g not in labels:
                raise ValueError(f"{where}: {g!r} is not a genome")
        for g, s in x.strand.items():
            if s not in _FLIP or (s == "." and g not in x.members):
                raise ValueError(f"{where}: {g}'s strand {s!r}: +, - or ., and . only for a member (its own "
                                 "annotation's unknown strand; a call's strand is the caller's, turned)")
        for g, words in x.flags.items():
            if any(w in _DERIVED for w in words):
                raise ValueError(f"{where}: {g}'s flags {words}: {', '.join(_DERIVED)} are set by the engine "
                                 "(inv by a reversed window: use inversion or strands)")


class _Places:
    """Where each genome holds each locus: its seqid name, the site of each state (inside its sequence), the
    element's strand and the orientation of each pair's window, with the inversion applied."""

    def __init__(self, loci: list[ToyLocus], seqids, names: dict[str, list[str]], lengths: dict[str, int],
                 inversion: tuple[str, str, int, int] | None):
        self.loci, self.names, self.lengths, self.inv = loci, names, lengths, inversion
        self.index = {s: j for j, s in enumerate(seqids)}

    def where(self, i: int, g: str) -> tuple[str, int]:
        """(seqid, pos) where genome g holds locus i, before any inversion: `at`, else the locus's own."""
        x = self.loci[i]
        return x.at.get(g, (x.seqid, x.pos))

    def mirrored(self, i: int, g: str) -> bool:
        (seqid, pos), inv = self.where(i, g), self.inv
        return inv is not None and g == inv[0] and seqid == inv[1] and inv[2] <= pos <= inv[3]

    def strand(self, i: int, g: str) -> str:
        """Genome g's element strand at locus i (a member's own)."""
        s = self.loci[i].strand.get(g, "+")
        return _FLIP[s] if self.mirrored(i, g) else s

    def reversed(self, i: int, a: str, b: str) -> bool:
        """Whether the window of A's element at locus i in B runs on B's reverse strand: exactly one of them
        mirrored, flipped again where both given strands are known and differ."""
        x = self.loci[i]
        sa, sb = x.strand.get(a, "+"), x.strand.get(b, "+")
        return (self.mirrored(i, a) != self.mirrored(i, b)) != ("." not in (sa, sb) and sa != sb)

    def site(self, i: int, g: str, state: str, b_ins: int | None = None) -> tuple[str | None, int | None,
                                                                                    int | None]:
        """(seqid name, start, end) of genome g's site for `state` at locus i (module docstring); a
        ValueError if it leaves the sequence."""
        x = self.loci[i]
        seqid, pos = self.where(i, g)
        if state in ("P", "T"):
            s, e = pos, pos + x.length - 1
        elif state == "S":
            s, e = pos, pos + x.ltr - 1
        elif state == "E":
            s, e = pos - E_SITE, pos - 1
        elif state == "D":
            s, e = pos - 1, pos + (b_ins or 0)
        elif state == "O":
            if b_ins is None or b_ins < 1:
                raise ValueError(f"locus {i}: an O call needs B's inserted sequence (b_ins >= 1), not "
                                 f"{b_ins}")
            s, e = pos, pos + b_ins - 1
        else:
            return None, None, None
        if self.mirrored(i, g):
            a, b = self.inv[2], self.inv[3]
            s, e = a + b - e, a + b - s
        name = self.names[g][self.index[seqid]]
        if not 1 <= s <= e <= self.lengths[seqid]:
            raise ValueError(f"locus {i} ({seqid}:{pos}): {g}'s {state} site {name}:{s}-{e} leaves the "
                             f"sequence (1-{self.lengths[seqid]})")
        return name, s, e


def _elements(loci: list[ToyLocus], annotated: list[str], place: _Places) -> dict[str, dict[int, Element]]:
    """Each annotated genome's elements, by locus index: one per locus it is a member of."""
    out: dict[str, dict[int, Element]] = {}
    for g in annotated:
        mine = [i for i, x in enumerate(loci) if g in x.members]
        eids = {}
        for i in mine:
            eids[i] = make_eid(g, *place.site(i, g, "P"))
        if len(set(eids.values())) != len(eids):
            raise ValueError(f"{g}: two member loci at the same place")
        out[g] = {}
        for i in mine:
            x = loci[i]
            seqid, s, e = place.site(i, g, "P")
            host = x.parent if x.parent is not None and x.parent in eids else None
            out[g][i] = Element(genome=g, eid=eids[i], name=f"{g}_toy{i + 1}", seqid=seqid, start=s, end=e,
                                strand=place.strand(i, g), ltr_left=(s, s + x.ltr - 1),
                                ltr_right=(e - x.ltr + 1, e), tsd=TSD, k=x.k, family=x.family,
                                superfamily=x.superfamily, clade=None, confidence=None,
                                parent=None if host is None else eids[host], level=_level(loci, i, eids))
    return out


def _level(loci: list[ToyLocus], i: int, mine) -> int:
    """Nesting depth of locus i's element among the genome's member loci `mine` (a container of locus
    indexes; 0 = top level)."""
    level, seen = 0, {i}
    while loci[i].parent is not None and loci[i].parent in mine:
        i = loci[i].parent
        if i in seen:
            raise ValueError(f"locus {i}: its parent links form a cycle")
        seen.add(i)
        level += 1
    return level


def _calls(a: str, b: str, loci: list[ToyLocus], mine: dict[int, Element], place: _Places) -> list[str]:
    """A's pair file into B: the header, then one CALL_COLUMNS line per A element, sorted as the engine
    sorts them (A's fai order, start, -end, level)."""
    rows = []
    for i, el in mine.items():
        x = loci[i]
        state = x.calls.get((a, b)) or ("P" if b in x.members else x.states[b])
        dl, dr, b_ins = x.evidence.get(b) or EVIDENCE.get(state, (None, None, None))
        seqid, s, e = place.site(i, b, state, b_ins)
        rev = state != "U" and place.reversed(i, a, b)  # U is no_anchor: its window is unresolved, no inv
        given = x.flags.get(b, [])
        proof = [_PROOF[state]] if state in _PROOF and not any(f in _INHERITED for f in given) else []
        flags = list(dict.fromkeys([*proof, *given, *(["inv"] if rev else [])]))
        row = {"eid": el.eid, "a_seqid": el.seqid, "a_start": el.start, "a_end": el.end, "parent": el.parent,
               "level": el.level, "state": state, "tier": TIER[state], "flags": ",".join(flags) or None,
               "reason": "no_anchor" if state == "U" else None, "b_seqid": seqid, "b_start": s, "b_end": e,
               "b_strand": None if seqid is None else "-" if rev else "+",
               "wid": i if x.wid is None else x.wid, "dl": dl, "dr": dr,
               "b_ins": b_ins, "tsd_a": TSD}
        key = (place.index[place.where(i, a)[0]], el.start, -el.end, el.level)
        rows.append((key, "\t".join([_field(row.get(c)) for c in CALL_COLUMNS]) + "\n"))
    return ["\t".join(CALL_COLUMNS) + "\n", *(line for _, line in sorted(rows))]


def _field(v) -> str:
    """tsv_field(v), with its two commonest cases (None, a string) first."""
    return "." if v is None else v if type(v) is str else tsv_field(v)


def _merge(loci: list[ToyLocus], labels: list[str], elements: dict[str, dict[int, Element]],
           call_paths: dict[tuple[str, str], Path], samples: list[Sample],
           names: dict[str, list[str]]) -> list[Locus]:
    """build_loci over the written pair files, its loci in list order (lids renumbered); a ValueError where
    the engine merges two toy loci, splits one, or merges a locus's calls to states other than its own."""
    order = {x: {name: j for j, name in enumerate(names[x])} for x in labels}
    merged, _ = build_loci(samples, {a: list(els.values()) for a, els in elements.items()}, call_paths, order)
    owner = {(g, el.eid): i for g, els in elements.items() for i, el in els.items()}
    rows: list[Locus | None] = [None] * len(loci)
    for locus in merged:
        found = {owner[(g, eid)] for g, eid in locus.members.items()}
        found |= {owner[(g, eid)] for g, eids in locus.extra.items() for eid in eids}
        if len(found) > 1:
            raise ValueError(f"loci {sorted(found)}: the engine merges them into one (a P, S or T call of "
                             "one lands on another's element)")
        i = found.pop()
        if rows[i] is not None:
            raise ValueError(f"locus {i}: the engine splits it (no P, S or T call links its members "
                             f"{sorted(loci[i].members)})")
        rows[i] = locus
    for i, (x, locus) in enumerate(zip(loci, rows, strict=True)):
        got = {g: locus.states[g].state for g in labels}
        if got != x.states:
            raise ValueError(f"locus {i}: its calls merge to {got} (loci.py), not its states {x.states}")
        locus.lid = f"L{i + 1:07d}"
    return rows


def _write_sheet(path: Path, genomes, fastas: dict[str, Path]) -> Path:
    """samples.tsv with paths relative to its directory; annotated genomes get ann/<label>.gff3 (never
    opened)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = ["label\tfasta\tannotation\tgroup\trole"]
    for label, group, role, ann in genomes:
        rows.append("\t".join((label, os.path.relpath(fastas[label], path.parent),
                               f"ann/{label}.gff3" if ann else ".", group or ".", role or ".")))
    path.write_text("\n".join(rows) + "\n")
    return path


def _n_fasta(path: Path, names: list[str], lengths: list[int]) -> Path:
    """An N-only FASTA of the given sequences at `path`, written once (atomically) and reused after."""
    if not path.is_file():
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
        with open(tmp, "w") as fh:
            for name, n in zip(names, lengths, strict=True):
                full, rest = divmod(n, _LINE)
                fh.write(f">{name}\n" + ("N" * _LINE + "\n") * full + ("N" * rest + "\n" if rest else ""))
        os.replace(tmp, path)
    return path
