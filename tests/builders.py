"""Small synthetic-sequence builders shared by tests (seeded through the caller's Random)."""
from __future__ import annotations

import contextlib
import csv
import gzip
import json
import math
import random
import re
import shutil
import subprocess
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import NamedTuple

import numpy as np

from synltr.cli import main as _cli_main
from synltr.io.annotation.nesting import assign_nesting
from synltr.io.elements import Element, ElementTable, make_eid
from synltr.io.seqio import revcomp


def write_fasta(path: Path, records: dict[str, str]) -> Path:
    """A FASTA of `records` (name -> sequence) in 80-column lines, with the parent directory made; returns
    the path."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for name, seq in records.items():
            fh.write(f">{name}\n")
            fh.writelines(seq[i:i + 80] + "\n" for i in range(0, len(seq), 80))
    return path


def rand_seq(n: int, rng: random.Random) -> str:
    """Uniform random ACGT of length n."""
    return "".join(rng.choice("ACGT") for _ in range(n))


def mutate(seq: str, rate: float, rng: random.Random) -> str:
    """Substitute each base independently with probability `rate` by one of the other three bases."""
    return "".join(rng.choice([c for c in "ACGT" if c != b]) if rng.random() < rate else b for b in seq)


def make_element(rng: random.Random, ltr_len: int, int_len: int) -> tuple[str, int, int]:
    """LTR + INT + LTR with identical LTRs: (sequence, left LTR length, right LTR length)."""
    ltr = rand_seq(ltr_len, rng)
    return ltr + rand_seq(int_len, rng) + ltr, ltr_len, ltr_len


FLANK = 3000  # length of the flanks L and R
SPACER = 20_000  # build_pair spacer length
T = "CTCTC"  # TSD of the top-level elements
T2 = "GAGAG"  # TSD of the host's child and of adjacent_E's second element
T_MM1 = "CTGTC"  # T with 1 substitution (E_mm1)
T_MM2 = "CACAC"  # T with 2 substitutions (E_mm2)
INT_SPLIT = 2000  # the host's internal region holds its child (or its empty site) after this many bp


def _nested(elements: list[Element]) -> list[Element]:
    """The elements with `parent` and `level` set, sorted by (seqid, start, -end)."""
    return list(assign_nesting(ElementTable.from_rows(elements)))


class _Label(NamedTuple):
    state: str
    flags: tuple[str, ...] = ()  # must be present
    noflags: tuple[str, ...] = ()  # must be absent
    tier: str | None = None
    reason: str | None = None  # U reason


# case -> element name -> expected call
_LABELS: dict[str, dict[str, _Label]] = {
    "P": {"e": _Label("P")},
    "P_nestB": {"e": _Label("P", flags=("nestB",))},
    "S": {"e": _Label("S", flags=("solo_tsd",))},
    "S_notsd": {"e": _Label("S", noflags=("solo_tsd",))},
    "T": {"e": _Label("T")},
    "E": {"e": _Label("E", flags=("tsd_match",), tier="high")},
    "E_mm1": {"e": _Label("E", flags=("tsd_match",))},
    "E_mm2": {"e": _Label("E", noflags=("tsd_match",), tier="medium")},
    "D_left": {"e": _Label("D")},
    "D_both": {"e": _Label("D")},
    "O": {"e": _Label("O")},
    "G": {"e": _Label("G")},
    "G_near": {"e": _Label("G")},
    "G_far": {"e": _Label("E")},
    "U_flank": {"e": _Label("U", reason="flank_unaligned")},
    "inv_site": {"e": _Label("U", reason="inv_site")},
    "nested_E": {"h": _Label("P"), "c": _Label("E")},
    "host_E": {"h": _Label("E"), "c": _Label("E", flags=("host_empty",))},
    "host_S": {"h": _Label("S"), "c": _Label("D")},
    "host_G": {"h": _Label("G"), "c": _Label("G", flags=("inherit",))},
    "adjacent_E": {"e1": _Label("E", flags=("adjacent",)), "e2": _Label("E", flags=("adjacent",))},
}
CASES: list[str] = list(_LABELS)
_HOST = ("nested_E", "host_E", "host_S", "host_G")


@dataclass
class ToyCase:
    """One A/B window pair; the dicts are keyed by element eid."""

    a: str
    b: str
    elements: list[Element]
    expect: dict[str, str]  # expected state
    flags: dict[str, set[str]] = field(default_factory=dict)  # flags that must be present
    noflags: dict[str, set[str]] = field(default_factory=dict)  # flags that must be absent
    tiers: dict[str, str] = field(default_factory=dict)  # expected tier
    reasons: dict[str, str] = field(default_factory=dict)  # expected U reason


@dataclass
class ToyPair(ToyCase):
    """Toy cases laid end to end on seqid chr1 (same fields as ToyCase)."""


def _el(name: str, off: int, length: int, ltr: int, tsd: str) -> Element:
    """Toy element at 0-based offset `off` of A (the record is 1-based inclusive)."""
    s, e = off + 1, off + length
    return Element(genome="toy", eid=make_eid("toy", "a", s, e), name=name, seqid="a", start=s, end=e,
                   strand="+", ltr_left=(s, s + ltr - 1), ltr_right=(e - ltr + 1, e), tsd=tsd, k=0.0,
                   family=None, superfamily=None, clade=None, confidence=None)


def _finish(name: str, a: str, b: str, elements: list[Element]) -> ToyCase:
    out = ToyCase(a, b, _nested(elements), {})
    for e in out.elements:
        lab = _LABELS[name][e.name]
        out.expect[e.eid] = lab.state
        if lab.flags:
            out.flags[e.eid] = set(lab.flags)
        if lab.noflags:
            out.noflags[e.eid] = set(lab.noflags)
        if lab.tier:
            out.tiers[e.eid] = lab.tier
        if lab.reason:
            out.reasons[e.eid] = lab.reason
    return out


def _single(name: str, L: str, R: str, rng: random.Random) -> ToyCase:
    """A = L T e T R with e = LTR(400) + INT(4000) + LTR."""
    e, n, _ = make_element(rng, 400, 4000)
    ltr, inner = e[:n], e[n:-n]
    b = {  # lambdas: only the chosen case draws filler from rng
        "P": lambda: mutate(L + T + e + T + R, 0.01, rng),
        "P_nestB": lambda: L + T + e[:2400] + rand_seq(3000, rng) + e[2400:] + T + R,
        "S": lambda: L + T + ltr + T + R,
        "S_notsd": lambda: L + T + ltr + R,
        "T": lambda: L + T + ltr + inner[:1600] + ltr + T + R,
        "E": lambda: L + T + R,
        "E_mm1": lambda: L + T_MM1 + R,
        "E_mm2": lambda: L + T_MM2 + R,
        "D_left": lambda: L[:-500] + R,
        "D_both": lambda: L[:-300] + R[400:],
        "O": lambda: L + T + rand_seq(2000, rng) + T + R,
        "G": lambda: L + T + "N" * 100 + R,
        "G_near": lambda: L[:-70] + "N" * 50 + L[-20:] + T + R,
        "G_far": lambda: L[:-550] + "N" * 50 + L[-500:] + T + R,
        "U_flank": lambda: L + T + rand_seq(3000, rng),
        "inv_site": lambda: L + T + revcomp(e) + T + R,
    }[name]()
    return _finish(name, L + T + e + T + R, b, [_el("e", len(L) + len(T), len(e), n, T)])


def _host(name: str, L: str, R: str, rng: random.Random) -> ToyCase:
    """A = L T h T R with h = LTRh(400) + INTh[:INT_SPLIT] + T2 + c + T2 + INTh[INT_SPLIT:] + LTRh, INTh 4 kb,
    c = LTRc(300) + INTc(2400) + LTRc."""
    hseq, n, _ = make_element(rng, 400, 4000)
    c, m, _ = make_element(rng, 300, 2400)
    ltr, inner = hseq[:n], hseq[n:-n]
    h = ltr + inner[:INT_SPLIT] + T2 + c + T2 + inner[INT_SPLIT:] + ltr
    b = {
        "nested_E": L + T + ltr + inner[:INT_SPLIT] + T2 + inner[INT_SPLIT:] + ltr + T + R,
        "host_E": L + T + R,
        "host_S": L + T + ltr + T + R,
        "host_G": L + T + "N" * 100 + R,
    }[name]
    off = len(L) + len(T)
    els = [_el("h", off, len(h), n, T), _el("c", off + n + INT_SPLIT + len(T2), len(c), m, T2)]
    return _finish(name, L + T + h + T + R, b, els)


def _adjacent(L: str, R: str, rng: random.Random) -> ToyCase:
    """A = L T e1 T T2 e2 T2 R with two distinct elements; B = L T T2 R."""
    e1, n1, _ = make_element(rng, 400, 4000)
    e2, n2, _ = make_element(rng, 400, 4000)
    o1 = len(L) + len(T)
    o2 = o1 + len(e1) + len(T) + len(T2)
    els = [_el("e1", o1, len(e1), n1, T), _el("e2", o2, len(e2), n2, T2)]
    return _finish("adjacent_E", L + T + e1 + T + T2 + e2 + T2 + R, L + T + T2 + R, els)


def case(name: str, rng: random.Random) -> ToyCase:
    """Toy case `name` (one of CASES); draws L and R, then the element(s), then any filler."""
    if name not in _LABELS:
        raise ValueError(f"unknown toy case {name!r}")
    L, R = rand_seq(FLANK, rng), rand_seq(FLANK, rng)
    if name in _HOST:
        return _host(name, L, R, rng)
    if name == "adjacent_E":
        return _adjacent(L, R, rng)
    return _single(name, L, R, rng)


def build_pair(rng: random.Random, cases: list[str] | None = None) -> ToyPair:
    """Toy cases end to end on seqid chr1: spacer, case, spacer, ..., case, spacer.

    Each spacer is a fresh 20-kb random sequence, the same in A and B. The default cases are CASES
    without U_flank (its U only arises when a window ends inside the case). Draw order: the cases,
    in order, then the spacers.
    """
    names = [n for n in CASES if n != "U_flank"] if cases is None else list(cases)
    built = [case(n, rng) for n in names]
    spacers = [rand_seq(SPACER, rng) for _ in range(len(built) + 1)]
    out = ToyPair("", "", [], {})
    a_parts, b_parts, els = [spacers[0]], [spacers[0]], []
    off = SPACER  # A offset of the current case
    for n, c, sp in zip(names, built, spacers[1:]):
        eid = {}
        for e in c.elements:
            s, t = e.start + off, e.end + off
            eid[e.eid] = make_eid("toy", "chr1", s, t)
            els.append(replace(e, eid=eid[e.eid], name=f"{n}.{e.name}", seqid="chr1", start=s, end=t,
                               ltr_left=(e.ltr_left[0] + off, e.ltr_left[1] + off),
                               ltr_right=(e.ltr_right[0] + off, e.ltr_right[1] + off)))
        for src, dst in ((c.expect, out.expect), (c.flags, out.flags), (c.noflags, out.noflags),
                         (c.tiers, out.tiers), (c.reasons, out.reasons)):
            dst.update((eid[k], v) for k, v in src.items())
        a_parts += [c.a, sp]
        b_parts += [c.b, sp]
        off += len(c.a) + SPACER
    out.a, out.b, out.elements = "".join(a_parts), "".join(b_parts), _nested(els)
    return out


# ---- locus-catalog fixtures: genomes A (reference), B and C on chr1 and chr2

ABC_ORDER = {g: {"chr1": 0, "chr2": 1} for g in "ABC"}  # every genome's fai order
ABC_LENGTHS = {g: {"chr1": 1000, "chr2": 500} for g in "ABC"}


def abc_samples():
    """The sheet of genomes A, B and C (FASTA and annotation paths are placeholders)."""
    from synltr.io.samples import Sample

    return [Sample(g, Path(f"/x/{g}.fa"), f"/x/{g}", None, None) for g in "ABC"]


def genome_state(state, tier="high", seqid=None, start=None, end=None, strand=None, k=None, flags=(), reason=None):
    """A locus's state in one genome (coordinates 1-based closed)."""
    from synltr.calls.loci import GenomeState

    return GenomeState(state, tier, seqid, start, end, strand, k, list(flags), reason)


def make_locus(lid, a, b, c, family="RLG_1", **kw):
    """A locus with member element A:e and the states `a`, `b`, `c` of genomes A, B and C; `kw` replaces fields."""
    from synltr.calls.loci import Locus

    fields = {"lid": lid, "members": {"A": "A:e"}, "discovery": ["A"], "family": family, "superfamily": "Gypsy",
              "clade": "Ale", "length": 5000, "ltr_len": 300, "tsd": "ACGTA", "states": {"A": a, "B": b, "C": c}}
    return Locus(**{**fields, **kw})


# ---- planted LTR-retrotransposons and the formats that carry them

LTRQUEST_COLUMNS = ["seq_id", "seq_len", "status", "ltr5_start", "ltr5_end", "ltr3_start", "ltr3_end", "ltr5_len", "ltr3_len", "flank5_len", "flank3_len", "aln_len", "n_sites", "n_ts", "n_tv", "n_gapcols", "identity", "p_dist", "k2p", "k2p_se", "bitscore", "flank_margin_bits", "cigar", "motif", "k2p_time", "orientation", "tsd", "tsd_offset", "tsd_input", "strand", "family", "confidence", "domains", "nest_status"]
CLEAN_KINDS = ("none", "cut", "masked")


def plant(genome: str, element: str, tsd: str, pos: int) -> str:
    """`genome` with `tsd + element + tsd` inserted before the 0-based offset `pos`. `element` is the sequence as
    it lies in the genome (whatever its strand label), so the element starts at
    0-based `pos + len(tsd)`, which is 1-based `pos + len(tsd) + 1`."""
    return genome[:pos] + tsd + element + tsd + genome[pos:]


def make_ltr_element(rng: random.Random, left_len: int, right_len: int, int_len: int) -> tuple[str, int, int]:
    """LTR + INT + LTR in genome orientation (sequence, left LTR length, right LTR length). The right LTR is the
    left one cut to `right_len`, or extended with random bases when it is longer, so the LTRs can differ in length."""
    left = rand_seq(left_len, rng)
    right = left[:right_len] + rand_seq(max(0, right_len - left_len), rng)
    return left + rand_seq(int_len, rng) + right, left_len, right_len


@dataclass
class ElementSpec:
    """An element to plant: LTR lengths in genome orientation, internal length, strand label, TSD, its class and
    k, and `nests`, a list of (offset into the internal region, ElementSpec) planted inside it with their own TSDs."""

    left: int = 300
    right: int = 300
    internal: int = 1500
    strand: str = "+"
    tsd: str = "ACGTA"
    nests: list[tuple[int, ElementSpec]] = field(default_factory=list)
    family: str | None = "fam1"
    superfamily: str | None = "Gypsy"
    clade: str | None = "Ale"
    k: float | None = 0.01


def _render(spec: ElementSpec, rng: random.Random):
    """(sequence, records): each record is (spec, start, end, left LTR (s, e), right LTR (s, e)), 0-based half-open
    and relative to the sequence's first base, the element itself first."""
    seq, nl, nr = make_ltr_element(rng, spec.left, spec.right, spec.internal)
    body, out, last, records = seq[:nl], seq[nl:-nr], 0, []
    inner = [body]
    shift = nl
    for off, child in sorted(spec.nests, key=lambda n: n[0]):
        cseq, crec = _render(child, rng)
        inner.append(out[last:off])
        shift += off - last + len(child.tsd)
        records += [(s, a + shift, b + shift, (c0 + shift, c1 + shift), (d0 + shift, d1 + shift))
                    for s, a, b, (c0, c1), (d0, d1) in crec]
        inner.append(child.tsd + cseq + child.tsd)
        shift += len(cseq) + len(child.tsd)
        last = off
    inner.append(out[last:])
    full = "".join(inner) + seq[-nr:]
    return full, [(spec, 0, len(full), (0, nl), (len(full) - nr, len(full)))] + records


def build_scene(rng: random.Random, specs: list[ElementSpec], seqid: str = "chr1", genome: str = "G",
                flank: int = 300, gap: int = 300) -> tuple[dict[str, str], list]:
    """One chromosome `seqid`: random flank, then each top-level spec as `tsd + element + tsd` (nests inside it)
    separated by random gaps, then a flank. Returns ({seqid: sequence}, elements), the elements (1-based closed,
    nests included in a host's extent) with `parent` and `level` set, sorted by (start, -end)."""
    chrom, elements = rand_seq(flank, rng), []
    for spec in specs:
        seq, records = _render(spec, rng)
        off = len(chrom) + len(spec.tsd)  # 0-based start of the element in the chromosome
        for s, a, b, (c0, c1), (d0, d1) in records:
            start, end = off + a + 1, off + b
            elements.append(Element(
                genome=genome, eid=make_eid(genome, seqid, start, end), name=f"{seqid}:{start}-{end}", seqid=seqid,
                start=start, end=end, strand=s.strand, ltr_left=(off + c0 + 1, off + c1),
                ltr_right=(off + d0 + 1, off + d1), tsd=s.tsd, k=s.k, family=s.family, superfamily=s.superfamily,
                clade=s.clade, confidence="high", ltr_left_len=c1 - c0, ltr_right_len=d1 - d0, tsd_source="annot"))
        chrom = plant(chrom, seq, spec.tsd, len(chrom)) + rand_seq(gap, rng)
    chrom += rand_seq(flank - gap, rng) if flank > gap else ""
    return {seqid: chrom}, _nested(elements)


def write_gff3(path: Path, elements: list, lengths: bool = False) -> Path:
    """LTRquest-convention GFF3 of `elements` (1-based closed, as the elements): an LTR_retrotransposon line and two
    long_terminal_repeat lines (Parent = its ID) each. With `lengths` the host line also carries ltr5_len and
    ltr3_len, the LTRs' own lengths (their extents' when the element does not give them). Returns the path."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["##gff-version 3\n"]
    for e in elements:
        opt = [("Name", e.name), ("tsd", e.tsd), ("K2P_d", None if e.k is None else repr(e.k)),
               ("family", e.family), ("superfamily", e.superfamily), ("clade", e.clade),
               ("confidence", e.confidence)]
        if lengths:
            opt += [("ltr5_len", str(e.ltr_left_len or e.ltr_left[1] - e.ltr_left[0] + 1)),
                    ("ltr3_len", str(e.ltr_right_len or e.ltr_right[1] - e.ltr_right[0] + 1))]
        attrs = f"ID={e.eid}" + "".join(f";{k}={v}" for k, v in opt if v is not None)
        head, tail = f"{e.seqid}\tbuilder\t", f"\t.\t{e.strand}\t.\t"
        lines.append(f"{head}LTR_retrotransposon\t{e.start}\t{e.end}{tail}{attrs}\n")
        for tag, (s, t) in (("lLTR", e.ltr_left), ("rLTR", e.ltr_right)):
            lines.append(f"{head}long_terminal_repeat\t{s}\t{t}{tail}ID={e.eid}.{tag};Parent={e.eid}\n")
    path.write_text("".join(lines))
    return path


def _clean_record(span: str, e, children: list, clean: str) -> str:
    """The element's record (genome orientation) under the `clean` convention: `cut` removes the direct nests'
    spans, `masked` replaces them with N, `none` keeps them."""
    holes = sorted((c.start - e.start, c.end - e.start + 1) for c in children)
    if clean == "none":
        return span
    if clean == "masked":
        for a, b in holes:
            span = span[:a] + "N" * (b - a) + span[b:]
        return span
    out, last = [], 0
    for a, b in holes:
        out.append(span[last:a])
        last = b
    out.append(span[last:])
    return "".join(out)


def write_ltrquest(prefix: str, elements: list, records: dict[str, str], clean: str = "none") -> list[Path]:
    """LTRquest library files `{prefix}_depth<level>_clean_ltr.tsv` and `.fa`, one pair per nesting depth, from
    `elements` (with `level` and `parent` set) cut out of the genome `records` (seqid -> sequence).

    Coordinates: the record header `seqid:start-end#LTR/<superfamily>/<clade>` is 1-based closed. The record is
    the cleaned span in the element's own orientation, so a minus-strand element's record is the reverse
    complement of the genome span (as LTRquest writes it). The LTR columns are 1-based offsets from the span's
    left end in genome orientation on either strand (ltr5 is the genome-left LTR), in every `clean` mode and
    whatever the record's orientation or length: a cut record is shorter than the span the offsets count along.
    `clean` says what a host's record holds in place of its direct nests: "none" the nests, "masked" N of the same length
    (the record is as long as the span), "cut" nothing (the record is shorter than the span). The clean step
    happens before any reverse complement. `seq_len` is the record's length."""
    if clean not in CLEAN_KINDS:
        raise ValueError(f"clean must be one of {CLEAN_KINDS}, not {clean!r}")
    kids: dict[str, list] = {}
    for e in elements:
        if e.parent is not None:
            kids.setdefault(e.parent, []).append(e)
    paths = []
    for depth in sorted({e.level for e in elements}):
        rows, fa = [], []
        for e in (x for x in elements if x.level == depth):
            span = records[e.seqid][e.start - 1:e.end]
            rec = _clean_record(span, e, kids.get(e.eid, []), clean)
            if e.strand == "-":
                rec = revcomp(rec)
            header = f"{e.seqid}:{e.start}-{e.end}#LTR/{e.superfamily or 'unknown'}/{e.clade or 'unknown'}"
            c = dict.fromkeys(LTRQUEST_COLUMNS, "NA")
            c.update(seq_id=header, seq_len=len(rec), status="pass", ltr5_start=e.ltr_left[0] - e.start + 1,
                     ltr5_end=e.ltr_left[1] - e.start + 1,
                     ltr3_start=e.ltr_right[0] - e.start + 1,
                     ltr3_end=e.ltr_right[1] - e.start + 1,
                     ltr5_len=e.ltr_left[1] - e.ltr_left[0] + 1, ltr3_len=e.ltr_right[1] - e.ltr_right[0] + 1,
                     k2p="NA" if e.k is None else e.k, tsd=e.tsd or "NA", strand=e.strand, family=e.family or "NA",
                     confidence=e.confidence or "NA", nest_status="." if e.parent is None else f"nest-inner:{e.parent}")
            rows.append("\t".join(str(c[k]) for k in LTRQUEST_COLUMNS))
            fa.append(f">{header}\n" + "".join(rec[i:i + 60] + "\n" for i in range(0, len(rec), 60)))
        table = Path(f"{prefix}_depth{depth}_clean_ltr.tsv")
        table.parent.mkdir(parents=True, exist_ok=True)
        table.write_text("#" + "\t".join(LTRQUEST_COLUMNS) + "\n" + "".join(r + "\n" for r in rows))
        table.with_suffix(".fa").write_text("".join(fa))
        paths.append(table)
    return paths


def pav_toy(tmp_path: Path, seed: int = 1):
    """LTRpav's toy pair: `build_pair` over its default cases, written under `tmp_path` as A.fa and B.fa (seqid
    chr1) and A.gff3 (A's elements). Returns the ToyPair with `fasta_a`, `fasta_b` and `gff3` set."""
    pair = build_pair(random.Random(seed))
    pair.fasta_a = write_fasta(tmp_path / "A.fa", {"chr1": pair.a})
    pair.fasta_b = write_fasta(tmp_path / "B.fa", {"chr1": pair.b})
    pair.gff3 = write_gff3(tmp_path / "A.gff3", pair.elements)
    return pair


def write_backbone(fasta: Path, gff3: Path) -> int:
    """Three 300 kb random sequences with 13 planted LTR retrotransposons each (and one nested in the first of
    s1), as a FASTA and a GFF3 (LTR_retrotransposon with two long_terminal_repeat children, 1-based closed).
    Returns the number of elements. `synltr sim pop` evolves the quickstart's four genomes on it."""
    rng = random.Random(7)
    starts = [10_000 + 22_000 * k for k in range(13)]
    n = 0
    with open(fasta, "w") as fa, open(gff3, "w") as gff:
        gff.write("##gff-version 3\n")
        for j, seqid in enumerate(("s1", "s2", "s3")):
            r = random.Random(j)
            seq = bytearray("".join(r.choices("ACGT", k=300_000)).encode())
            spans = [(s, rng.randint(200, 600), rng.randint(2_000, 6_000)) for s in starts]
            if seqid == "s1":
                spans.append((starts[0] + spans[0][1] + 500, 150, 800))
            for i, (s, ltr, internal) in enumerate(sorted(spans)):  # s: 0-based start
                e = s + 2 * ltr + internal  # 0-based, half-open end
                lt = "".join(r.choices("ACGT", k=ltr)).encode()
                seq[s:e] = lt + "".join(r.choices("ACGT", k=internal)).encode() + lt
                seq[e:e + 5] = seq[s - 5:s]  # target-site duplication
                n += 1
                eid, strand = f"{seqid}e{i}", "-" if i % 3 == 2 else "+"
                col = f"{seqid}\ttoy\t{{}}\t{{}}\t{{}}\t.\t{strand}\t.\t{{}}\n"
                gff.write(col.format("LTR_retrotransposon", s + 1, e,
                                     f"ID={eid};tsd={seq[s - 5:s].decode()};K2P_d=0.01;family=fam{i % 4};"
                                     "superfamily=Gypsy;confidence=high"))
                gff.write(col.format("long_terminal_repeat", s + 1, s + ltr, f"Parent={eid}"))
                gff.write(col.format("long_terminal_repeat", e - ltr + 1, e, f"Parent={eid}"))
            fa.write(f">{seqid}\n" + "\n".join(seq[k:k + 80].decode() for k in range(0, len(seq), 80)) + "\n")
    return n


# ---- the assembly route: a genome pair's context and its call files

@contextlib.contextmanager
def pair_context(d: Path, a_seq: str, b_seq: str, elements: list, threads: int = 2):
    """PairContext (routes.assembly.anchor) of A (a_seq) vs B (b_seq), each one sequence chr1 written to d; outdir
    d / "out"."""
    from synltr.io.samples import Sample
    from synltr.io.seqio import Fasta
    from synltr.routes.assembly.anchor import PairContext

    d.mkdir(parents=True, exist_ok=True)
    a, b = write_fasta(d / "A.fa", {"chr1": a_seq}), write_fasta(d / "B.fa", {"chr1": b_seq})
    with Fasta(a) as fa, Fasta(b) as fb:
        yield PairContext(a=Sample("A", a, None, None, None), b=Sample("B", b, None, None, None), a_fa=fa,
                          b_fa=fb, a_elements=list(elements), b_elements=None, outdir=d / "out", threads=threads)


def read_rows(path: Path, columns: tuple[str, ...] | None = None) -> list[dict]:
    """The rows of a gzip TSV with a header line (CALL_COLUMNS unless `columns` is given), values as written."""
    from synltr.calls.records import CALL_COLUMNS

    with gzip.open(path, "rt") as fh:
        header = tuple(fh.readline().rstrip("\n").split("\t"))
        rows = [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh]
    assert header == (CALL_COLUMNS if columns is None else columns)
    return rows


def read_calls(path: Path) -> dict[str, dict]:
    """{eid: row} of a pair call file, values as written ("." for none)."""
    rows = read_rows(path)
    calls = {r["eid"]: r for r in rows}
    assert len(calls) == len(rows)
    return calls


# ---- the reads route: synthetic genomes with planted LTR-RT alleles, simulated long reads and their alignments

def _ltr_family(rng: random.Random, ltr_len: int, int_len: int) -> dict[str, str]:
    """A family's LTR (TG...CA) and internal region."""
    return {"ltr": "TG" + rand_seq(ltr_len - 4, rng) + "CA", "int": rand_seq(int_len, rng)}


def _ltr_copy(rng: random.Random, fam: dict[str, str], div: float) -> dict[str, str]:
    """A copy of a family, each base substituted with probability `div`."""
    return {"ltr": mutate(fam["ltr"], div, rng), "int": mutate(fam["int"], div, rng)}


def _allele_seq(cp: dict[str, str], inner=None, inner_state: str = "F") -> tuple[str, tuple | None]:
    """The full element of copy `cp`; `inner` = (offset in the internal region, copy) is planted inside it with its
    own 5-bp TSD when `inner_state` is F. Returns (sequence, (offset, length, LTR length) of the inner element or
    None)."""
    internal, info = cp["int"], None
    if inner is not None and inner_state == "F":
        q, icp = inner
        iseq = icp["ltr"] + icp["int"] + icp["ltr"]
        internal = internal[:q + 5] + iseq + internal[q:q + 5] + internal[q + 5:]
        info = (len(cp["ltr"]) + q + 5, len(iseq), len(icp["ltr"]))
    return cp["ltr"] + internal + cp["ltr"], info


def _haplotype(anc: str, events: list[dict], hap: str, ann: list | None = None, chrom: str = "chr1") -> str:
    """One haplotype of ancestral sequence `anc`: each event's allele on `hap` at anc[p] with a 5-bp TSD
    (T = anc[p:p+5]): F full (T + element + T), Fx full with an internal deletion, S solo LTR (T + LTR + T), E empty
    (anc unchanged), X anc[x0:x1] deleted (the site and its flanks missing). With `ann`, the elements of the
    haplotype are appended to it as dicts (0-based half-open coordinates on this haplotype)."""
    out, cur, n = [], 0, 0
    for ev in sorted(events, key=lambda e: e["p"]):
        st, p = ev["state"][hap], ev["p"]
        if st == "X":
            x0, x1 = ev["xdel"]
            out.append(anc[cur:x0])
            n += x0 - cur
            cur = x1
            continue
        out.append(anc[cur:p + 5])
        n += p + 5 - cur
        cur = p + 5
        tsd = anc[p:p + 5]
        if st == "F":
            seq, info = _allele_seq(ev["copy"], ev.get("inner"), ev.get("inner_state", {}).get(hap, "F"))
            if ann is not None:
                full, ltr = len(seq), len(ev["copy"]["ltr"])
                ann.append({"id": ev["id"], "chrom": chrom, "s": n, "e": n + full, "l5": (n, n + ltr),
                            "l3": (n + full - ltr, n + full), "fam": ev["fam"]})
                if info:
                    off, il, ll = info
                    ann.append({"id": ev["id"] + "_inner", "chrom": chrom, "s": n + off, "e": n + off + il,
                                "l5": (n + off, n + off + ll), "l3": (n + off + il - ll, n + off + il),
                                "fam": ev["inner_fam"]})
            out.append(seq + tsd)
            n += len(seq) + 5
        elif st == "Fx":
            seq, _ = _allele_seq(ev["copy"])
            ltr = len(ev["copy"]["ltr"])
            d0, d1 = ev["int_del"]
            seq = seq[:ltr + d0] + seq[ltr + d1:]
            out.append(seq + tsd)
            n += len(seq) + 5
        elif st == "S":
            out.append(ev["copy"]["ltr"] + tsd)
            n += len(ev["copy"]["ltr"]) + 5
    out.append(anc[cur:])
    return "".join(out)


def _write_read_gff3(path: Path, ann: list[dict]) -> Path:
    """GFF3 of the planted elements (an LTR_retrotransposon with two long_terminal_repeat children each; the dicts
    are 0-based half-open, the file 1-based closed)."""
    lines = ["##gff-version 3\n"]
    for a in sorted(ann, key=lambda r: (r["chrom"], r["s"])):
        lines.append(f"{a['chrom']}\tsynth\tLTR_retrotransposon\t{a['s'] + 1}\t{a['e']}\t.\t+\t.\t"
                     f"ID={a['id']};family={a['fam']};clade=test\n")
        lines += [f"{a['chrom']}\tsynth\tlong_terminal_repeat\t{x + 1}\t{y}\t.\t+\t.\tID={a['id']}.{tag};"
                  f"Parent={a['id']}\n" for tag, (x, y) in (("lLTR", a["l5"]), ("rLTR", a["l3"]))]
    path.write_text("".join(lines))
    return path


def _simulate_reads(rng: random.Random, seqs: dict[str, str], depth: float, mean_len: int, err: float, fh,
                    tag: str) -> int:
    """FASTQ reads drawn uniformly from each sequence (lognormal lengths, at least 1.5 kb), errors split equally
    into substitutions, insertions and deletions, half of them reverse-complemented. Returns the read count."""
    n_tot = 0
    for name, s in seqs.items():
        for k in range(int(depth * len(s) / mean_len)):
            ln = max(1500, min(len(s), int(rng.lognormvariate(math.log(mean_len) - 0.06, 0.35))))
            st = rng.randrange(0, len(s) - ln + 1)
            r = s[st:st + ln]
            if err > 0:
                o = []
                for b in r:
                    u = rng.random()
                    if u < err / 3:
                        o.append(rng.choice([c for c in "ACGT" if c != b]))
                    elif u < 2 * err / 3:
                        o.append(b + rng.choice("ACGT"))
                    elif u < err:
                        continue
                    else:
                        o.append(b)
                r = "".join(o)
            if rng.random() < 0.5:
                r = revcomp(r)
            fh.write(f"@{tag}_{name}_{k}\n{r}\n+\n{'I' * len(r)}\n")
            n_tot += 1
    return n_tot


def _align_reads(ref: Path, fastq: Path, bam: Path, primary: Path | None = None) -> None:
    """minimap2 (map-hifi) alignments of `fastq` on `ref`, sorted and indexed into `bam`; with `primary`, also the
    primary alignments only (no secondary or supplementary records), indexed. Without minimap2 or samtools the
    calling test skips (pytest is imported only then, so the quick start runs without it)."""
    if shutil.which("minimap2") is None or shutil.which("samtools") is None:
        import pytest

        pytest.skip("minimap2 and samtools (to align the simulated reads) not on PATH")
    cmd = (f"minimap2 -ax map-hifi -t 2 {ref} {fastq} 2>/dev/null | samtools sort -o {bam} - 2>/dev/null && "
           f"samtools index {bam}")
    if primary is not None:
        cmd += f" && samtools view -b -F 0x904 -o {primary} {bam} && samtools index {primary}"
    subprocess.run(cmd, shell=True, check=True)


@dataclass
class ReadSet:
    """A reads-route scenario: reference FASTA, its elements (GFF3), the read alignments (`bam`: every record;
    `primary`: primary records only, ploidy 2), and the true genotype of each element."""

    ref: Path
    gff: Path
    bam: Path
    primary: Path | None
    truth: dict[str, str]
    ploidy: int


TETRA_DOSES = {  # element -> allele on haplotypes h1..h4 (h1 is the reference)
    "D4F": "FFFF", "D3F1E": "FFFE", "D2F2E": "FFEE", "D1F3E": "FEEE", "D2F2S": "FFSS", "D3F1N": "FFFN"}


def read_set(tmp_path: Path, ploidy: int = 2, *, seed: int = 1, het: float = 0.0) -> ReadSet:
    """ltrkin's synthetic read sets, written under `tmp_path` (ref.fa, elements.gff3, reads.fq, full.bam and, for
    ploidy 2, primary.bam). Skips the calling test without minimap2 or samtools.

    Ploidy 2: a diploid whose haplotype h1 is the reference (chr1 160 kb, chr2 60 kb, chr3 40 kb), with the
    elements E1_het_empty F/E, E2_het_solo F/S, E3_hom F/F, E4_hemi F/N (h2 lacks 16 kb around it),
    E5_nested_outer F/F holding E5_nested_outer_inner F/E, E7_hom and E8_decoy_hom F/F, and E9_internal_del F/F (h2
    lacks 1.5 kb of its internal region); 15x HiFi-like reads (0.5% error) per haplotype; `het` substitutes bases of
    h2 at that rate. Ploidy 4: an autotetraploid with the allele dosages of TETRA_DOSES on four haplotypes 0.4%
    apart, 10x per haplotype (`seed` and `het` unused)."""
    tmp_path.mkdir(parents=True, exist_ok=True)
    if ploidy == 4:
        return _tetraploid(tmp_path)
    if ploidy != 2:
        raise ValueError(f"read_set builds ploidy 2 or 4, not {ploidy}")
    rng = random.Random(seed)
    fam_a, fam_b = _ltr_family(rng, 400, 4000), _ltr_family(rng, 600, 6000)
    a1, a2 = rand_seq(160_000, rng), rand_seq(60_000, rng)
    a3 = rand_seq(40_000, random.Random(seed + 7))  # its own stream: the other sequences stay as they were

    def cp(fam):
        return _ltr_copy(rng, fam, 0.002)

    ev1 = [
        {"id": "E1_het_empty", "p": 20000, "copy": cp(fam_a), "fam": "famA", "state": {"h1": "F", "h2": "E"}},
        {"id": "E2_het_solo", "p": 45000, "copy": cp(fam_a), "fam": "famA", "state": {"h1": "F", "h2": "S"}},
        {"id": "E3_hom", "p": 70000, "copy": cp(fam_b), "fam": "famB", "state": {"h1": "F", "h2": "F"}},
        {"id": "E4_hemi", "p": 95000, "copy": cp(fam_a), "fam": "famA", "state": {"h1": "F", "h2": "X"},
         "xdel": (87000, 103000)},
        {"id": "E5_nested_outer", "p": 125000, "copy": cp(fam_b), "fam": "famB", "state": {"h1": "F", "h2": "F"},
         "inner": (3000, cp(fam_a)), "inner_fam": "famA", "inner_state": {"h1": "F", "h2": "E"}},
    ]
    ev2 = [{"id": "E7_hom", "p": 15000, "copy": cp(fam_b), "fam": "famB", "state": {"h1": "F", "h2": "F"}},
           {"id": "E8_decoy_hom", "p": 40000, "copy": cp(fam_b), "fam": "famB", "state": {"h1": "F", "h2": "F"}}]
    ev3 = [{"id": "E9_internal_del", "p": 18000, "copy": _ltr_copy(random.Random(seed + 8), fam_a, 0.002),
            "fam": "famA", "state": {"h1": "F", "h2": "Fx"}, "int_del": (1000, 2500)}]
    ann: list[dict] = []
    h1 = {"chr1": _haplotype(a1, ev1, "h1", ann, "chr1"), "chr2": _haplotype(a2, ev2, "h1", ann, "chr2"),
          "chr3": _haplotype(a3, ev3, "h1", ann, "chr3")}
    h2 = {"chr1": _haplotype(a1, ev1, "h2"), "chr2": _haplotype(a2, ev2, "h2"), "chr3": _haplotype(a3, ev3, "h2")}
    if het > 0:
        h2 = {c: mutate(v, het, random.Random(seed + 100)) for c, v in h2.items()}
    truth = {"E1_het_empty": "F/E", "E2_het_solo": "F/S", "E3_hom": "F/F", "E4_hemi": "F/N",
             "E5_nested_outer": "F/F", "E5_nested_outer_inner": "F/E", "E7_hom": "F/F", "E8_decoy_hom": "F/F",
             "E9_internal_del": "F/F"}
    ref = tmp_path / "ref.fa"
    ref.write_text("".join(f">{c}\n" + "\n".join(s[i:i + 80] for i in range(0, len(s), 80)) + "\n"
                           for c, s in h1.items()))
    gff = _write_read_gff3(tmp_path / "elements.gff3", ann)
    fq = tmp_path / "reads.fq"
    with open(fq, "w") as fh:
        _simulate_reads(rng, h1, 15, 12000, 0.005, fh, "h1")
        _simulate_reads(rng, h2, 15, 12000, 0.005, fh, "h2")
    _align_reads(ref, fq, tmp_path / "full.bam", tmp_path / "primary.bam")
    return ReadSet(ref, gff, tmp_path / "full.bam", tmp_path / "primary.bam", truth, 2)


def _tetraploid(d: Path) -> ReadSet:
    """ltrkin's autotetraploid (read_set with ploidy 4)."""
    rng = random.Random(9)
    fam = _ltr_family(rng, 450, 4500)
    anc = rand_seq(40000 * len(TETRA_DOSES) + 20000, rng)
    events = []
    for k, (eid, st) in enumerate(TETRA_DOSES.items()):
        p = 20000 + 40000 * k
        ev = {"id": eid, "p": p, "copy": _ltr_copy(rng, fam, 0.003), "fam": "f",
              "state": {f"h{i + 1}": ("X" if a == "N" else a) for i, a in enumerate(st)}}
        if "N" in st:
            ev["xdel"] = (p - 8000, p + 14000)
        events.append(ev)
    ann: list[dict] = []
    haps = {}
    for i in range(4):
        h = _haplotype(anc, events, f"h{i + 1}", ann if i == 0 else None, "chr1")
        haps[f"h{i + 1}"] = h if i == 0 else mutate(h, 0.004, random.Random(100 + i))  # 4 distinct haplotypes
    ref = d / "ref.fa"
    ref.write_text(">chr1\n" + haps["h1"] + "\n")
    gff = _write_read_gff3(d / "el.gff3", ann)
    fq = d / "reads.fq"
    with open(fq, "w") as fh:
        for name, h in haps.items():
            _simulate_reads(rng, {"chr1": h}, 10, 12000, 0.005, fh, name)
    _align_reads(ref, fq, d / "reads.bam")
    truth = {eid: "/".join(sorted(st, key="FSEN".index)) for eid, st in TETRA_DOSES.items()}
    return ReadSet(ref, gff, d / "reads.bam", None, truth, 4)


def reads_sheet(path: Path, ref: Path, annotation: Path, reads: dict[str, Path], ref_label: str = "ref") -> Path:
    """A samples sheet for the reads route: a reference row (`ref_label`, FASTA and annotation) and one row per
    read set (label -> BAM). Returns the path."""
    rows = [f"{ref_label}\t{ref}\t{annotation}\treference\t."] + [f"{k}\t.\t.\t.\t{v}" for k, v in reads.items()]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("label\tfasta\tannotation\trole\treads\n" + "".join(r + "\n" for r in rows))
    return path


# ---- polyploid scenario: a `sim nucleus` library with planted unplaced scaffolds and a decoy compartment ----
# `sim nucleus` allotetraploid, two chromosomes per subgenome: split 2 Myr, merger 100 kyr. The insertion rate keeps
# the library small (928 LTR-RTs) while every chromosome finds its homoeolog and the split is dated from
# >= SPLIT_MIN inherited pairs.
SUBGENOMES_TOY = {"root_age": 6e6, "chromosomes": {"n": 2, "length": 4e6}, "families": {"n": 60},
                  "activity": {"rate": 0.000225, "p_burst": 0.0}, "exchanges": {"nonreciprocal": 0, "reciprocal": 0},
                  "tree": {"split": 2e6, "children": ["A", "B"]}, "mergers": [{"time": 1e5, "lineages": ["A", "B"]}]}


def simulate(scenario, outdir, truthdir, seed):
    """Run `synltr sim nucleus` in this process on `scenario` (a dict of the scenario file's keys): the library in
    `outdir` (ltr.fa, k2p.tsv, annotation.gff3, genome.fai), the truth in `truthdir`. The scenario file is written
    next to `outdir`. Returns (outdir, truthdir)."""
    outdir, truthdir = Path(outdir), Path(truthdir)
    sc = outdir.parent / f"scenario_{seed}.json"
    sc.write_text(json.dumps(scenario))
    code = _cli_main(["sim", "nucleus", str(sc), "-o", str(outdir), "--truthdir", str(truthdir), "--seed", str(seed)])
    if code:
        raise RuntimeError(f"sim nucleus exited {code}")
    return outdir, truthdir


def _read_fasta(p):
    name, seq = None, []
    with open(p) as f:
        for line in f:
            if line.startswith(">"):
                if name:
                    yield name, "".join(seq)
                name, seq = line[1:].strip(), []
            else:
                seq.append(line.strip())
    if name:
        yield name, "".join(seq)


def _parse(h):
    loc, _, cls = h.partition("#")
    chrom, _, se = loc.rpartition(":")
    s, e = (int(x) for x in se.split("-"))
    return chrom, s, e, cls


def plant_scaffolds(outdir, newdir, scaffold_frac=0.1, n_scaffolds=3, n_decoys=8, per_decoy=15, decoy_div=0.03,
                    seed=1):
    """Copy a simulation into newdir with (1) the last scaffold_frac of every chromosome cut into n_scaffolds unplaced
    scaffolds (copies renamed and shifted; scaffold -> parent chromosome recorded), (2) n_decoys small sequences each
    holding per_decoy copies of one family's young elements carrying a shared private mutation set (decoy_div) plus
    rare own mutations. Annotation = the K2P table (record-relative LTR coordinates need no change).
    Returns (newdir, {scaffold: parent chromosome}, [decoy names])."""
    rnd = random.Random(seed)
    outdir, newdir = Path(outdir), Path(newdir)
    newdir.mkdir(parents=True, exist_ok=True)
    with open(outdir / "genome.fai") as fh:
        fai = {r[0]: int(r[1]) for r in csv.reader(fh, delimiter="\t")}
    with open(outdir / "k2p.tsv") as fh:
        k2p_rows = list(csv.reader(fh, delimiter="\t"))
    head, body = k2p_rows[0], {r[0]: r for r in k2p_rows[1:]}
    seqs = dict(_read_fasta(outdir / "ltr.fa"))
    parent, newfai, recs, krows = {}, {}, [], []
    for c, L in fai.items():
        cut = int(L * (1 - scaffold_frac))
        newfai[c] = cut
        w = (L - cut) // n_scaffolds
        for k in range(n_scaffolds):
            nm = f"{c}_scf{k + 1}"
            parent[nm] = c
            newfai[nm] = w if k < n_scaffolds - 1 else L - cut - w * (n_scaffolds - 1)
    for h, s in seqs.items():
        c, a, b, cls = _parse(h)
        L = fai[c]
        cut = newfai[c]
        if a <= cut and b <= cut:
            nh = h
        else:
            if a <= cut:                                      # straddles the cut: drop (assemblies break there)
                continue
            w = (L - cut) // n_scaffolds
            k = min((a - cut - 1) // w, n_scaffolds - 1)
            off = cut + k * w
            nm = f"{c}_scf{k + 1}"
            if b - off > newfai[nm]:
                continue
            nh = f"{nm}:{a - off}-{b - off}#{cls}"
        recs.append((nh, s))
        r = list(body[h])
        r[0] = nh
        krows.append(r)
    # decoy compartment: one family's youngest copies, a shared private mutation set, rare own mutations
    fam_col = head.index("family")
    k_col = head.index("k2p")
    by_fam = {}
    for h, r in body.items():
        by_fam.setdefault(r[fam_col], []).append(h)
    fam = max(by_fam, key=lambda f: len(by_fam[f]))
    young = sorted(by_fam[fam], key=lambda h: float(body[h][k_col]))[:max(per_decoy * n_decoys, 20)]
    ref = seqs[young[0]]
    shared = {i: rnd.choice([x for x in "ACGT" if x != ref[i].upper()]) for i in range(len(ref))
              if rnd.random() < decoy_div}
    decoys = []
    for d in range(n_decoys):
        nm = f"decoy{d + 1}"
        decoys.append(nm)
        pos = 1000
        for j in range(per_decoy):
            h0 = young[(d * per_decoy + j) % len(young)]
            s = list(seqs[young[0]])
            for i, x in shared.items():
                s[i] = x
            for i in range(len(s)):
                if rnd.random() < 0.001:
                    s[i] = rnd.choice("ACGT")
            s = "".join(s)
            cls = _parse(h0)[3]
            nh = f"{nm}:{pos}-{pos + len(s) - 1}#{cls}"
            recs.append((nh, s))
            r = list(body[young[0]])
            r[0] = nh
            krows.append(r)
            pos += len(s) + 5000
        newfai[nm] = pos + 1000
    with open(newdir / "ltr.fa", "w") as f:
        for h, s in recs:
            f.write(f">{h}\n{s}\n")
    with open(newdir / "k2p.tsv", "w") as f:
        f.write("\t".join(head) + "\n")
        for r in krows:
            f.write("\t".join(r) + "\n")
    with open(newdir / "genome.fai", "w") as f:
        for c, L in newfai.items():
            f.write(f"{c}\t{L}\t0\t0\t0\n")
    return newdir, parent, decoys


def truth_subgenomes(truthdir):
    """{chromosome token: subgenome} of a `sim nucleus` truth folder."""
    with open(Path(truthdir) / "chromosomes.tsv") as fh:
        return {r["token"]: r["subgenome"] for r in csv.DictReader(fh, delimiter="\t")}


# ---- the library route: three genomes ((R,S),D) as LTR-RT libraries (FASTA records and tables) with known truth ----
TRIO_GENOMES = "RSD"
_TRIO_ACGT = np.array(list(b"ACGT"), np.uint8)
TRIO_OUT_SPLIT = 0.03  # divergence at which the outgroup splits from ((R,S),D): older insertions are its too


def _trio_mut(seq, d, rng):
    """`seq` (codes 0-3) with Binomial(len, min(d, 0.75)) distinct positions each changed to another base."""
    s = seq.copy()
    n = rng.binomial(len(s), min(d, 0.75))
    if n:
        pos = rng.choice(len(s), n, replace=False)
        s[pos] = (s[pos] + rng.integers(1, 4, n).astype(np.uint8)) % 4
    return s


def _trio_element(fam, ltr, rng):
    """A new insertion: the family master with 2% copy-specific noise, its two LTRs identical at insertion."""
    core = _trio_mut(fam, 0.02, rng)
    left = core[:ltr].copy()
    return np.concatenate([left, core[ltr:], left])


@dataclass
class TrioScenario:
    """A library scenario written by `trio`: the samples sheet (columns label, library = "FASTA,TABLE", role), each
    genome's library, the truth per locus (chrom number, position in R and S, class, carriers, in_tract), the planted
    tracts (chrom number, start, end in R) and the inversion (chrom number, start, end in R; D holds it reversed)."""

    sheet: Path
    libraries: dict[str, str]
    truth: list[dict]
    tracts: list[tuple[int, int, int]]
    inversion: tuple[int, int, int] | None


def trio(outdir: Path, *, seed: int = 1, admix: bool = True, outgroup: bool = False, n_chrom: int = 3,
         chrom_len: int = 3_000_000, spacing: int = 8000, inversion: bool = True, decoy: bool = True,
         t_rs: float = 0.0015, t_split: float = 0.006, t_deep: float = 0.02, ltr: int = 250,
         elen: int = 2000) -> TrioScenario:
    """Tiny LTR-RT libraries of three genomes with species tree ((R,S),D) and known truth, in `outdir`:
    <G>.fa (records `<G>_chr<c>:<start>-<end>#LTR/Gypsy/fam<f>`, 1-based closed) and <G>.tsv (seq_id, k2p, tsd,
    strand, ltr5_len, ltr3_len), and samples.tsv. Divergences are substitutions per site.

    Loci share positions across genomes (D's are shifted by 100 kb and one D segment is inverted). Classes: old
    (inserted before every split: R, S, D), rs (R-S branch), rd and sd (incomplete lineage sorting, equally many:
    D = 0 without gene flow), private (one genome). `admix`: inside three tracts R carries a recent copy of D's
    haplotype (D's elements at near-identity, none of its own R-S insertions). `decoy`: a high-copy family of
    near-identical copies at random positions in every genome (no orthology), family 9. `outgroup` adds genome O
    (role outgroup in the sheet), split from ((R,S),D) at divergence TRIO_OUT_SPLIT: it carries the old insertions
    older than that and its own private ones, positions shifted by 200 kb, drawn from a second stream so R, S and D
    stay as without it. Paths in the sheet and in `libraries` are absolute."""
    outdir = Path(outdir).resolve()
    rng = np.random.default_rng(seed)
    outdir.mkdir(parents=True, exist_ok=True)
    fams = [rng.integers(0, 4, elen).astype(np.uint8) for _ in range(6)]
    loci = []   # (chrom, pos, class, family)
    classes = ["old"] * 30 + ["rs"] * 25 + ["rd"] * 6 + ["sd"] * 6 + ["pR"] * 11 + ["pS"] * 11 + ["pD"] * 11
    for c in range(1, n_chrom + 1):
        for p in range(50_000, chrom_len - 50_000, spacing):
            loci.append((c, p + int(rng.integers(-2000, 2000)), classes[rng.integers(len(classes))],
                         int(rng.integers(5))))
    tracts = [(1, 900_000, 1_700_000), (2, 2_200_000, 2_600_000), (3, 400_000, 1_000_000)] if admix else []
    inv = (2, 1_000_000, 1_800_000) if inversion else None
    recs = {g: [] for g in TRIO_GENOMES}
    truth, roots = [], []
    for (c, p, cl, fi) in loci:
        tsd = "".join("ACGT"[x] for x in rng.integers(0, 4, 5))
        in_tract = any(c == t[0] and t[1] <= p <= t[2] for t in tracts)
        if cl == "old":
            age, carriers = t_deep + rng.exponential(0.01), "RSD"
        elif cl == "rs":
            age, carriers = rng.uniform(t_rs, t_split), "RS"
        elif cl in ("rd", "sd"):
            age, carriers = t_split + rng.exponential(0.005), ("RD" if cl == "rd" else "SD")
        else:
            age, carriers = rng.uniform(0, t_rs), cl[1]
        base = _trio_element(fams[fi], ltr, rng)
        # LTR-LTR divergence at the root = age (each LTR mutates age/2 independently)
        root = base.copy()
        root[:ltr] = _trio_mut(base[:ltr], age / 2, rng)
        root[-ltr:] = _trio_mut(base[-ltr:], age / 2, rng)
        root[ltr:-ltr] = _trio_mut(base[ltr:-ltr], age / 2, rng)
        # coalescence of carriers: R-S at t_rs, (R,S)-D deep; each carrier mutates from the shared node
        node_rs = _trio_mut(root, max(age - t_split, 0) / 2, rng) if "D" in carriers else root
        seqs = {}
        if "D" in carriers:
            seqs["D"] = _trio_mut(node_rs, min(age, t_split + 0.004) / 2, rng)
        rs_node = _trio_mut(node_rs, (t_split - t_rs) / 2, rng) if ("R" in carriers or "S" in carriers) else None
        for g in "RS":
            if g in carriers:
                seqs[g] = _trio_mut(rs_node, t_rs / 2, rng)
        if in_tract:
            # R carries D's haplotype here: D's element at near-identity; R's own R-S insertions are absent
            seqs.pop("R", None)
            if "D" in seqs:
                seqs["R"] = _trio_mut(seqs["D"], 0.0004, rng)
        for g, s in seqs.items():
            recs[g].append((c, p, s, tsd, age, fi))
        truth.append({"chrom": c, "pos": p, "cls": cl, "carriers": "".join(sorted(seqs)), "in_tract": in_tract})
        roots.append((root, tsd, age))
    if decoy:
        dec = rng.integers(0, 4, elen).astype(np.uint8)
        for g in TRIO_GENOMES:
            for _ in range(40):
                c = int(rng.integers(1, n_chrom + 1))
                p = int(rng.integers(60_000, chrom_len - 60_000))
                recs[g].append((c, p, _trio_mut(dec, 0.001, rng), "".join("ACGT"[x] for x in rng.integers(0, 4, 5)),
                                0.001, 9))
    genomes, shift = list(TRIO_GENOMES), {"D": 100_000}
    if outgroup:
        genomes.append("O")
        shift["O"] = 200_000
        recs["O"] = _outgroup_records(np.random.default_rng([seed, 1]), loci, roots, truth, fams, ltr, chrom_len)
    paths = {}
    for g in genomes:
        fa, tb = outdir / f"{g}.fa", outdir / f"{g}.tsv"
        with open(fa, "w") as f1, open(tb, "w") as f2:
            f2.write("seq_id\tk2p\ttsd\tstrand\tltr5_len\tltr3_len\n")
            for (c, p, s, tsd, age, fi) in sorted(recs[g], key=lambda r: (r[0], r[1])):
                pos = p + shift.get(g, 0)
                if g == "D" and inv and c == inv[0] and inv[1] <= p <= inv[2]:
                    pos = inv[1] + inv[2] - p + 100_000
                name = f"{g}_chr{c}:{pos}-{pos + len(s) - 1}#LTR/Gypsy/fam{fi}"
                f1.write(f">{name}\n{bytes(_TRIO_ACGT[s]).decode()}\n")
                f2.write(f"{name}\t{age:.5f}\t{tsd}\t+\t{ltr}\t{ltr}\n")
        paths[g] = f"{fa},{tb}"
    sheet = outdir / "samples.tsv"
    sheet.write_text("label\tlibrary\trole\n" + "".join(
        f"{g}\t{paths[g]}\t{'outgroup' if g == 'O' else '.'}\n" for g in genomes))
    return TrioScenario(sheet, paths, truth, tracts, inv)


def _outgroup_records(rng, loci, roots, truth, fams, ltr, chrom_len):
    """The outgroup's records: each old locus older than TRIO_OUT_SPLIT (its copy diverged from the insertion's root
    copy since the split), plus one private young insertion per nine loci, at a random position (its own TSD).
    The truth rows of the old loci it carries gain O."""
    out = []
    for (c, p, cl, fi), (root, tsd, age), row in zip(loci, roots, truth):
        if cl == "old" and age > TRIO_OUT_SPLIT:
            out.append((c, p, _trio_mut(root, min(age, TRIO_OUT_SPLIT + 0.004) / 2, rng), tsd, age, fi))
            row["carriers"] = "".join(sorted(row["carriers"] + "O"))
    for c, p, _, _ in loci[::9]:
        fi = int(rng.integers(5))
        age = float(rng.uniform(0, 0.0015))
        s = _trio_element(fams[fi], ltr, rng)
        s[:ltr] = _trio_mut(s[:ltr], age / 2, rng)
        s[-ltr:] = _trio_mut(s[-ltr:], age / 2, rng)
        q = min(p + 4000, chrom_len - 60_000)
        out.append((c, q, s, "".join("ACGT"[x] for x in rng.integers(0, 4, 5)), age, fi))
        truth.append({"chrom": c, "pos": q, "cls": "pO", "carriers": "O", "in_tract": False})
    return out


def trio_outgroup(outdir: Path, *, seed: int = 1, admix: bool = True) -> TrioScenario:
    """The trio-plus-outgroup library scenario: `trio` with genome O as the outgroup (samples sheet role outgroup)."""
    return trio(outdir, seed=seed, admix=admix, outgroup=True)


def trio_second_donor(outdir: Path, *, seed: int = 5, admix: bool = True) -> TrioScenario:
    """`trio` plus genome E, a sister of the donor D: D's library with light extra mutation (in each record
    max(1, length // 400) distinct positions, each changed to the next base in ACGT order; default_rng(9)), records
    and table rows renamed D_ -> E_ (E.fa, E.tsv). The samples sheet lists R, S, D, E; the truth is `trio`'s."""
    sc = trio(outdir, seed=seed, admix=admix)
    outdir = Path(outdir).resolve()
    rng = np.random.default_rng(9)
    fa, tb = sc.libraries["D"].split(",")
    with open(fa) as fi, open(outdir / "E.fa", "w") as fo:
        for rec in fi.read().split(">")[1:]:
            name, seq = rec.split("\n", 1)
            seq = list(seq.replace("\n", ""))
            for i in rng.choice(len(seq), max(1, len(seq) // 400), replace=False):
                seq[i] = "ACGT"[("ACGT".index(seq[i]) + 1) % 4]
            fo.write(f">{name.replace('D_', 'E_')}\n{''.join(seq)}\n")
    with open(tb) as fi, open(outdir / "E.tsv", "w") as fo:
        fo.writelines(line.replace("D_", "E_") for line in fi)
    libraries = {**sc.libraries, "E": f"{outdir / 'E.fa'},{outdir / 'E.tsv'}"}
    sheet = outdir / "samples.tsv"
    sheet.write_text("label\tlibrary\trole\n" + "".join(f"{g}\t{libraries[g]}\t.\n" for g in "RSDE"))
    return TrioScenario(sheet, libraries, sc.truth, sc.tracts, sc.inversion)


def screen_table_variants(outdir: Path) -> dict[str, str]:
    """Library specs ("FASTA,TABLE") of one tiny trio library (R of `trio`, seed 5, one 400-kb chromosome, no
    decoys) whose table varies within LTRmix's table format ("seq_id plus any of k2p, tsd, strand, ltr5_len,
    ltr3_len"): as is; '#' before the header; gzip; CRLF line ends; no K2P column; K2P only as K2P_d; K2P under the
    header K2P; a header-only table; a repeated seq_id (the last row differs); non-integer LTR lengths; an unplaced
    record (no coordinates in its name) with a row; seq_ids without the record names' #class (same intervals, other
    names); odd values (k2p inf, strand '?', TSD 'acgtn', an LTR length 0); a FASTA alone; k2p "NA" in every other row
    beside a K2P_d column holding a value in every row."""
    outdir = Path(outdir).resolve()
    sc = trio(outdir / "base", seed=5, n_chrom=1, chrom_len=400_000, decoy=False)
    fa0, tb0 = sc.libraries["R"].split(",")
    lines = Path(fa0).read_text().splitlines()
    recs = [(lines[i][1:], lines[i + 1]) for i in range(0, len(lines), 2)]
    head, *rows = (r.split("\t") for r in Path(tb0).read_text().splitlines())

    def write(name, records, header, table_rows, gz=False, eol="\n"):
        d = outdir / name
        d.mkdir(parents=True, exist_ok=True)
        fa, tb = d / ("R.fa.gz" if gz else "R.fa"), d / ("R.tsv.gz" if gz else "R.tsv")
        fa_text = "".join(f">{n}\n{s}\n" for n, s in records)
        tb_text = "\t".join(header) + eol + "".join("\t".join(r) + eol for r in table_rows)
        for path, text in ((fa, fa_text), (tb, tb_text)):
            if gz:
                with gzip.open(path, "wt", newline="") as fh:
                    fh.write(text)
            else:
                with open(path, "w", newline="") as fh:
                    fh.write(text)
        return f"{fa},{tb}"

    ki = head.index("k2p")
    dup = list(rows[0])
    dup[1:4] = ["0.33333", "GGGGG", "-"]
    frac = [list(r) for r in rows]
    frac[0][4:6], frac[1][4:6] = ["250.0", "250.0"], ["250.7", "249.2"]
    odd = [list(r) for r in rows]
    odd[0][1], odd[1][3], odd[2][2], odd[3][4] = "inf", "?", "acgtn", "0"
    kd_head = [*head[:ki + 1], "K2P_d", *head[ki + 1:]]
    kd_rows = [[*r[:ki], "NA" if j % 2 == 0 else r[ki], "0.0222" if j % 2 == 0 else "0.9", *r[ki + 1:]]
               for j, r in enumerate(rows)]
    alone = outdir / "fasta_alone"
    alone.mkdir(parents=True, exist_ok=True)
    shutil.copy(fa0, alone / "R.fa")
    return {
        "as_is": write("as_is", recs, head, rows),
        "hash_header": write("hash_header", recs, ["#" + head[0], *head[1:]], rows),
        "gz": write("gz", recs, head, rows, gz=True),
        "crlf_table": write("crlf_table", recs, head, rows, eol="\r\n"),
        "no_k2p_column": write("no_k2p_column", recs, [c for c in head if c != "k2p"],
                               [[v for j, v in enumerate(r) if j != ki] for r in rows]),
        "only_K2P_d": write("only_K2P_d", recs, ["K2P_d" if c == "k2p" else c for c in head], rows),
        "upper_K2P": write("upper_K2P", recs, ["K2P" if c == "k2p" else c for c in head], rows),
        "header_only": write("header_only", recs, head, []),
        "repeated_seq_id": write("repeated_seq_id", recs, head, [*rows, dup]),
        "non_integer_ltr_len": write("non_integer_ltr_len", recs, head, frac),
        "unplaced_with_row": write("unplaced_with_row", [*recs, ("contig_X", recs[3][1])], head,
                                   [*rows, ["contig_X", "0.01234", "ACGTA", "+", "250", "250"]]),
        "seq_id_without_class": write("seq_id_without_class", recs, head,
                                      [[r[0].split("#")[0], *r[1:]] for r in rows]),
        "odd_values": write("odd_values", recs, head, odd),
        "fasta_alone": str(alone / "R.fa"),
        "k2p_na_beside_K2P_d": write("k2p_na_beside_K2P_d", recs, kd_head, kd_rows),
    }


_GFF3_SOURCE = re.compile(rb"^([^#\t\n][^\t\n]*\tsynltr)-[^\t\n]*\t", re.MULTILINE)


def unversioned_gff3(data: bytes) -> bytes:
    """GFF3 bytes with `synltr` alone in the source column of each feature line (written `synltr-<version>`): how
    a pinned GFF3 file is compared, so that a version bump fails no pin."""
    return _GFF3_SOURCE.sub(rb"\1\t", data)


def subgenomes_text(path: Path) -> bytes:
    """A file of a subgenomes run as tests/golden/subgenomes_toy pins it: findings.txt without its header line (it
    names the version), a GFF3 file with `synltr` alone in its source column; every other byte as written."""
    data = path.read_bytes()
    if path.name == "findings.txt":
        return data.split(b"\n", 1)[1]
    return unversioned_gff3(data) if path.suffix == ".gff3" else data
