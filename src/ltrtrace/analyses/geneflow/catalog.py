"""Run directory -> Catalog: a finished run's outputs as numpy arrays, read as data; age bins;
nested-cluster keys.

Two kinds of state sit side by side, one row per locus (loci.tsv order), one column per genome (sheet order):
  merged  loci.tsv's state of each genome, merged across the discovering genomes by loci.py. A member genome
          is P (tier annot) whatever the calls say, and a member's own calls give the other genomes extra
          chances to be resolved, so merged states let annotation depth into a count: they serve the
          any-discovery descriptions only.
  own     each annotated genome A's own calls, pairs/<A>__<B>.calls.tsv.gz, joined to loci.tsv on A's member
          eid (the locus's first member element of A; A's other member elements are ignored): proj[A],
          ptier[A], pseq[A], pstart[A], pend[A], pbstrand[A]. A's own column holds its member
          element (P, tier high); rows where A has no member hold no call (0). A statistic with a discovery rule
          reads these (relations.discover).

Reading: loci.tsv in two passes (count the rows, then fill preallocated arrays); each pair file streamed
once, its evidence columns never parsed; sequence names and lengths from the FASTA indexes (cache work/fai);
the nesting (`top`) from work/elements/<label>.tsv. A missing run, or files that do not match the samples
sheet or each other, fail fast with an AnalysisError naming the file. The files are data: no genotyping
checkpoint is checked, so pair files without done markers are read as they stand. Nothing is written into the
run directory, which may be read-only.

Codes: state STATE_CODE (0 = no call); tier TIER_CODE; seq the genome's fai index (-1 = no coordinates);
start, end 1-based inclusive (0 = none); strand +1, -1, 0 (unknown or none); flags a FLAG_BIT mask. A call's
window orientation, b_strand, is kept as pbstrand (+1, -1; 0 for a call without a site), known whatever A's
own strand (a BED annotation writes "." for it).

Scale: proj, ptier, pseq, pstart, pend and pbstrand are dense (n, G) arrays per annotated genome A, 23
bytes per locus x genome x A: about 89 MB at 322,538 loci, 4 genomes, 3 annotated, about 2.3 GB at one
million loci and 10 genomes all annotated, the largest run this version targets. A pangenome-scale run needs
each A's calls stored on A's member rows only. One load reads every file (about 11 s at 322,538 loci): load
once per run.

Helpers of the statistics with a discovery rule (`check_trio`, `first_member`, `patterns`, `site_blocks`): the
trio check, each row's discovering genome, the R+D / S+D patterns and the jackknife blocks on the donor's sites.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ltrtrace.base.constants import AGE_EDGES
from ltrtrace.base.log import AnalysisError as AnalysesFailed
from ltrtrace.calls.records import stream_calls
from ltrtrace.calls.states import FLAGS, INFERRED_ABSENT, PAIR_STATES, call_class
from ltrtrace.io.elements import read_elements
from ltrtrace.io.samples import PAIR_SEP, Sample
from ltrtrace.io.seqio import Fasta
from ltrtrace.stats.resample import block_ids

# the states of the assembly route's pair calls (every pair state but the inferred absence), coded 1.. in the
# registry's order: P 1, S 2, T 3, E 4, D 5, O 6, G 7, U 8; 0 = no call
STATE_CODE = {x: c for c, x in enumerate((x for x in PAIR_STATES if x not in INFERRED_ABSENT), 1)}
TIER_CODE = {".": 0, "low": 1, "medium": 2, "high": 3, "annot": 4}
FLAG_BIT = {f: 1 << i for i, f in enumerate(FLAGS)}   # the call flag registry's words; other words are ignored
AGE_LABELS = ("K=0", "0<K<0.01", "K>=0.01")

# state code -> genotype class (calls.states.call_class, E the one absence): 1 present (P/S/T), 0 absent (E),
# -1 missing (D/O/G/U) or no call
_CLASS = {"present": 1, "absent": 0, "missing": -1}
_CLS = np.array([-1, *(_CLASS[call_class(x)] for x in STATE_CODE)], dtype=np.int8)
_STRAND = {"+": 1, "-": -1, ".": 0}
_CALL_TIER = {t: c for t, c in TIER_CODE.items() if t != "annot"}  # "annot" is a member's loci.tsv tier
_ATTRS = ("lid", "family", "superfamily", "length", "tsd")  # the loci.tsv attributes read
_PER_GENOME = ("eid", "state", "tier", "seqid", "start", "end", "strand", "k", "flags")  # <label>_<column>
_CALL = ("eid", "state", "tier", "b_seqid", "b_start", "b_end", "b_strand")  # columns read (evidence: never)
_FIRST = "run ltrtrace genotype first"
_REMERGE = "the run's files do not match each other (rerun ltrtrace genotype)"


class AnalysisError(AnalysesFailed, ValueError):
    """A run that is missing, or whose files do not match the samples sheet or each other; or analyses asked for
    genomes the run does not hold. The analyses failed after genotyping (exit status 3)."""


@dataclass(eq=False, repr=False)
class Catalog:
    """A finished run as arrays (module docstring): n loci in loci.tsv order, G genomes in sheet order."""

    outdir: Path
    samples: list[Sample]
    sheet: Path                         # the samples sheet's path
    labels: tuple[str, ...]             # sheet order
    ingroup: tuple[str, ...]            # sheet order: role other than "outgroup"
    outgroups: tuple[str, ...]          # sheet order: role "outgroup"
    annotated: dict[str, bool]
    groups: dict[str, str | None]
    seqids: dict[str, tuple[str, ...]]  # fai order, from seqio.Fasta (work/fai cache)
    lengths: dict[str, np.ndarray]      # int64, fai order
    # loci.tsv, merged across discoverers: any-discovery descriptions only
    lid: np.ndarray                     # (n,) str, loci.tsv order
    state: np.ndarray                   # (n, G) uint8 STATE_CODE
    tier: np.ndarray                    # (n, G) uint8 TIER_CODE
    seq: np.ndarray                     # (n, G) int32 fai index; -1 = no coordinates
    start: np.ndarray                   # (n, G) int64, 1-based inclusive; 0 = none
    end: np.ndarray                     # (n, G) int64, 1-based inclusive; 0 = none
    strand: np.ndarray                  # (n, G) int8: +1, -1, 0 (unknown or none)
    k: np.ndarray                       # (n, G) float64: a member's K2P, else nan
    flags: np.ndarray                   # (n, G) uint32 FLAG_BIT mask
    member: np.ndarray                  # (n, G) bool: the genome's annotation holds the locus (its discovery)
    eid: np.ndarray                     # (n, G) object: the member eid, else None
    family: np.ndarray                  # (n,) object, None = unknown (loci.tsv columns)
    superfamily: np.ndarray             # (n,) object, None = unknown
    tsd: np.ndarray                     # (n,) object, None = unknown
    length: np.ndarray                  # (n,) int64 (loci.tsv `length`)
    top: dict[str, dict[str, str]]      # genome -> member eid -> eid of its outermost annotated ancestor
                                        # (itself at level 0)
    # each annotated genome A's own calls: pairs/<A>__<B>.calls.tsv.gz joined on eid[:, g(A)] (extra members
    # ignored); A's own column holds its member element (P, tier high); rows where A has no member: no call
    proj: dict[str, np.ndarray]         # A -> (n, G) uint8 STATE_CODE of A's call into each genome; 0 = none
    ptier: dict[str, np.ndarray]        # same shape: TIER_CODE of those calls
    pseq: dict[str, np.ndarray]         # same shape: int32 fai index of the call's site there; -1 = none
    pstart: dict[str, np.ndarray]       # same shape: int64, 1-based inclusive; 0 = none
    pend: dict[str, np.ndarray]         # same shape: int64, 1-based inclusive; 0 = none
    pbstrand: dict[str, np.ndarray]     # same shape: int8, the call's b_strand: +1 the site runs on A's strand,
                                        # -1 reversed, 0 no site; A's own column +1 at its members

    def __post_init__(self) -> None:
        self._gi = {x: g for g, x in enumerate(self.labels)}

    def __repr__(self) -> str:
        return f"Catalog({str(self.outdir)!r}: {len(self.lid)} loci, genomes {', '.join(self.labels)})"

    def g(self, label: str) -> int:
        """The genome's column: its sheet index."""
        g = self._gi.get(label)
        if g is None:
            raise AnalysisError(f"{label!r} is not a genome of the run ({', '.join(self.labels)})")
        return g

    def cls(self, label: str) -> np.ndarray:
        """(n,) int8 from the merged state: 1 = P/S/T, 0 = E, -1 = D/O/G/U."""
        return _CLS[self.state[:, self.g(label)]]

    def mid(self, label: str) -> np.ndarray:
        """(n,) float64 (start + end) / 2 of the merged site; nan without one."""
        g = self.g(label)
        return np.where(self.seq[:, g] >= 0, (self.start[:, g] + self.end[:, g]) / 2, np.nan)

    def pcls(self, disc: str, label: str) -> np.ndarray:
        """(n,) int8 from disc's own calls into label: 1 P/S/T, 0 E, -1 D/O/G/U or no call (disc has no
        member there; everywhere for an unannotated disc)."""
        gl, gd = self.g(label), self.g(disc)
        p = self.proj.get(self.labels[gd])
        return np.full(len(self.lid), -1, dtype=np.int8) if p is None else _CLS[p[:, gl]]

    def pcls_by(self, gstar: np.ndarray, label: str) -> np.ndarray:
        """(n,) int8: row i's pcls(labels[gstar[i]], label); -1 where gstar[i] < 0."""
        gstar, gl = self._gstar(gstar), self.g(label)
        out = np.full(len(self.lid), -1, dtype=np.int8)
        for a, p in self.proj.items():
            rows = gstar == self._gi[a]
            out[rows] = _CLS[p[rows, gl]]
        return out

    def ptier_by(self, gstar: np.ndarray, label: str) -> np.ndarray:
        """(n,) uint8: row i's TIER_CODE of genome gstar[i]'s own call into label; 0 where gstar[i] < 0."""
        gstar, gl = self._gstar(gstar), self.g(label)
        out = np.zeros(len(self.lid), dtype=np.uint8)
        for a, t in self.ptier.items():
            rows = gstar == self._gi[a]
            out[rows] = t[rows, gl]
        return out

    def proj_site(self, gstar: np.ndarray, label: str) -> tuple[np.ndarray, np.ndarray]:
        """(fai index int64, mid float64) of label's site in row i's call from genome gstar[i] (in its own
        column: its member element); (-1, nan) where there is none."""
        gstar, gl = self._gstar(gstar), self.g(label)
        seq = np.full(len(self.lid), -1, dtype=np.int64)
        mid = np.full(len(self.lid), np.nan)
        for a, ps in self.pseq.items():
            rows = np.flatnonzero(gstar == self._gi[a])
            rows = rows[ps[rows, gl] >= 0]
            seq[rows] = ps[rows, gl]
            mid[rows] = (self.pstart[a][rows, gl] + self.pend[a][rows, gl]) / 2
        return seq, mid

    def _gstar(self, gstar: np.ndarray) -> np.ndarray:
        gstar = np.asarray(gstar)
        if gstar.shape != self.lid.shape:
            raise ValueError(f"gstar has shape {gstar.shape}; the catalog has {len(self.lid)} loci")
        return gstar


def age_bin(k: np.ndarray) -> np.ndarray:
    """int8 age bin of each K2P (AGE_LABELS): 0 K == 0; 1 0 < K < AGE_EDGES[1]; 2 K >= AGE_EDGES[1]; -1 nan
    (and a negative K, which no K2P is)."""
    k = np.asarray(k, dtype=np.float64)
    out = np.full(k.shape, -1, dtype=np.int8)
    out[k == AGE_EDGES[0]] = 0
    out[(k > AGE_EDGES[0]) & (k < AGE_EDGES[1])] = 1
    out[k >= AGE_EDGES[1]] = 2
    return out


def cluster_keys(cat: Catalog, gstar: np.ndarray, rows: np.ndarray) -> np.ndarray:
    """One key per selected row (rows: a boolean mask over the n loci; catalog order): (label of gstar[i], the
    top eid of that genome's member element there), so a host and the elements nested in it share a key. An
    index array for `rows`, or a selected row whose gstar genome has no member there, is a caller error."""
    gstar, rows = cat._gstar(gstar), np.asarray(rows)
    if rows.dtype != bool or rows.shape != cat.lid.shape:
        raise ValueError(f"cluster_keys: rows must be a boolean mask over the {len(cat.lid)} loci, not "
                         f"{rows.dtype} of shape {rows.shape} (an index array?)")
    idx = np.flatnonzero(rows)
    out = np.empty(len(idx), dtype=object)
    for j, i in enumerate(idx):
        g = int(gstar[i])
        e = cat.eid[i, g] if g >= 0 else None
        if e is None:
            raise ValueError(f"cluster_keys: locus {cat.lid[i]} has no member element of its discovering "
                             f"genome (gstar {g})")
        label = cat.labels[g]
        out[j] = (label, cat.top[label][e])
    return out


def load_catalog(outdir: Path, samples: list[Sample], sheet: Path) -> Catalog:
    """The finished run in `outdir`, made with the samples sheet at `sheet` (`samples`: its rows), as a
    Catalog (module docstring)."""
    outdir, sheet = Path(outdir), Path(sheet)
    labels = tuple(s.label for s in samples)
    if len(set(labels)) != len(labels):
        raise AnalysisError(f"duplicate genome label in {', '.join(labels)}")
    loci_path = outdir / "loci.tsv"
    if not loci_path.is_file():
        raise AnalysisError(f"no loci.tsv in {outdir}: {_FIRST}")
    width, attrs, gcols = _loci_columns(loci_path, labels, sheet)
    annotated = [s.label for s in samples if s.annotation is not None]
    elements = {a: outdir / "work" / "elements" / f"{a}.tsv" for a in annotated}
    pairs = {(a, b): outdir / "pairs" / f"{a}{PAIR_SEP}{b}.calls.tsv.gz"
             for a in annotated for b in labels if b != a}
    gone = [p.name for p in elements.values() if not p.is_file()]
    if gone:
        raise AnalysisError(f"element table(s) {', '.join(gone)} missing in {outdir / 'work' / 'elements'}: "
                            f"{_FIRST}")
    gone = [f"{a}{PAIR_SEP}{b}" for (a, b), p in pairs.items() if not p.is_file()]
    if gone:
        raise AnalysisError(f"{len(gone)} of {len(pairs)} pair file(s) missing in {outdir / 'pairs'}: "
                            f"{', '.join(gone)}; {_FIRST}")
    seqids, lengths = _fasta_lengths(samples, outdir / "work" / "fai")
    m = _read_loci(loci_path, width, attrs, gcols, labels, seqids)
    top: dict[str, dict[str, str]] = {x: {} for x in labels}
    for a in annotated:
        top[a] = _tops(elements[a], a, m["eid"][:, labels.index(a)])
    own = _read_own_calls(pairs, labels, annotated, m, seqids)
    return Catalog(outdir=outdir, samples=list(samples), sheet=sheet, labels=labels,
                   ingroup=tuple(s.label for s in samples if s.role != "outgroup"),
                   outgroups=tuple(s.label for s in samples if s.role == "outgroup"),
                   annotated={s.label: s.annotation is not None for s in samples},
                   groups={s.label: s.group for s in samples}, seqids=seqids, lengths=lengths, top=top,
                   **m, **own)


# ---------------------------------------------------------------- loci.tsv


def _loci_columns(path: Path, labels: tuple[str, ...], sheet: Path
                  ) -> tuple[int, dict[str, int], list[tuple[int, ...]]]:
    """loci.tsv's width, its _ATTRS columns and each genome's _PER_GENOME columns, after checking that its
    genomes are the sheet's, in sheet order."""
    with open(path, encoding="utf-8", newline="") as fh:
        header = fh.readline().rstrip("\r\n").split("\t")
    found = [c[:-len("_state")] for c in header if c.endswith("_state")]
    if found != list(labels):
        raise AnalysisError(f"{path}: its genomes ({', '.join(found) or 'none'}) are not those of the "
                            f"samples sheet {sheet} ({', '.join(labels)}): use the sheet the run was made "
                            "with")
    index = {c: i for i, c in enumerate(header)}
    lacking = [c for c in (*_ATTRS, *(f"{x}_{c}" for x in labels for c in _PER_GENOME)) if c not in index]
    if lacking:
        raise AnalysisError(f"{path}: column(s) {', '.join(lacking)} missing: not a loci.tsv of ltrtrace genotype")
    gcols = [tuple(index[f"{x}_{c}"] for c in _PER_GENOME) for x in labels]
    return len(header), {c: index[c] for c in _ATTRS}, gcols


def _count_rows(path: Path) -> tuple[int, int]:
    """loci.tsv's data rows and the width of its widest lid (the first pass)."""
    n = width = 0
    with open(path, "rb") as fh:
        fh.readline()
        for line in fh:
            n += 1
            width = max(width, line.find(b"\t"))
    return n, width


def _read_loci(path: Path, width: int, attrs: dict[str, int], gcols: list[tuple[int, ...]],
               labels: tuple[str, ...], seqids: dict[str, tuple[str, ...]]) -> dict[str, np.ndarray]:
    """loci.tsv as the Catalog's merged arrays: the rows counted first, then read into preallocated arrays."""
    n, lid_width = _count_rows(path)
    G = len(labels)
    m = {"lid": np.empty(n, dtype=f"<U{max(lid_width, 1)}"),
         "state": np.zeros((n, G), dtype=np.uint8), "tier": np.zeros((n, G), dtype=np.uint8),
         "seq": np.full((n, G), -1, dtype=np.int32), "start": np.zeros((n, G), dtype=np.int64),
         "end": np.zeros((n, G), dtype=np.int64), "strand": np.zeros((n, G), dtype=np.int8),
         "k": np.full((n, G), np.nan), "flags": np.zeros((n, G), dtype=np.uint32),
         "member": np.zeros((n, G), dtype=bool), "eid": np.full((n, G), None, dtype=object),
         "family": np.full(n, None, dtype=object), "superfamily": np.full(n, None, dtype=object),
         "tsd": np.full(n, None, dtype=object), "length": np.zeros(n, dtype=np.int64)}
    lid, family, superfamily, tsd, length = (
        m[c] for c in ("lid", "family", "superfamily", "tsd", "length"))
    state, tier, seq, start, end, strand, k, flags, member, eid = (
        m[c] for c in ("state", "tier", "seq", "start", "end", "strand", "k", "flags", "member", "eid"))
    c_lid, c_fam, c_sf, c_len, c_tsd = (attrs[c] for c in _ATTRS)
    orders = [{name: j for j, name in enumerate(seqids[x])} for x in labels]
    names: dict[str, str] = {}  # each distinct family, superfamily and TSD string, held once
    masks: dict[str, int] = {}  # flags field -> FLAG_BIT mask
    rows = 0
    with open(path, encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)
        next(reader)  # the header, checked by _loci_columns
        for i, row in enumerate(reader):
            line = reader.line_num
            if i >= n:
                raise AnalysisError(f"{path} changed while it was read")
            if len(row) != width:
                raise AnalysisError(f"{path}: line {line} has {len(row)} column(s); the header has {width}")
            lid[i] = row[c_lid]
            family[i], superfamily[i], tsd[i] = (None if v == "." else names.setdefault(v, v)
                                                 for v in (row[c_fam], row[c_sf], row[c_tsd]))
            try:
                length[i] = int(row[c_len])
            except ValueError:
                raise AnalysisError(f"{path}: line {line}: length is not an integer") from None
            for g, (ce, cs, ct, cq, cb, cn, cd, ck, cf) in enumerate(gcols):
                st, tr, sd = STATE_CODE.get(row[cs]), TIER_CODE.get(row[ct]), _STRAND.get(row[cd])
                if st is None or tr is None or sd is None:
                    raise AnalysisError(f"{path}: line {line}: {labels[g]}: unknown state, tier or strand "
                                        f"({row[cs]!r}, {row[ct]!r}, {row[cd]!r})")
                state[i, g], tier[i, g], strand[i, g] = st, tr, sd
                if row[ce] != ".":
                    eid[i, g], member[i, g] = row[ce], True
                if row[cq] != ".":
                    j = orders[g].get(row[cq])
                    if j is None:
                        raise AnalysisError(f"{path}: line {line}: {labels[g]}'s seqid {row[cq]!r} is not in "
                                            "its FASTA: the samples sheet names another FASTA than the run's")
                    try:
                        s, e = int(row[cb]), int(row[cn])
                    except ValueError:
                        s = e = 0
                    if not 1 <= s <= e:
                        raise AnalysisError(f"{path}: line {line}: {labels[g]}'s site {row[cq]}:{row[cb]}-"
                                            f"{row[cn]} is not two integers with 1 <= start <= end")
                    seq[i, g], start[i, g], end[i, g] = j, s, e
                if row[ck] != ".":
                    try:
                        k[i, g] = float(row[ck])
                    except ValueError:
                        raise AnalysisError(f"{path}: line {line}: {labels[g]}'s k {row[ck]!r} is not a "
                                            "number") from None
                mask = masks.get(row[cf])
                if mask is None:
                    mask = masks[row[cf]] = _flag_mask(row[cf])
                flags[i, g] = mask
            rows += 1
    if rows != n:
        raise AnalysisError(f"{path} changed while it was read")
    return m


def _flag_mask(field: str) -> int:
    mask = 0
    for word in field.split(","):
        mask |= FLAG_BIT.get(word, 0)
    return mask


# ---------------------------------------------------------------- FASTA indexes, elements


def _fasta_lengths(samples: list[Sample], cache: Path
                   ) -> tuple[dict[str, tuple[str, ...]], dict[str, np.ndarray]]:
    """Each genome's sequence names and lengths in fai order, from its index (built in `cache` if its
    directory is read-only); the sequence itself is never read."""
    seqids, lengths = {}, {}
    for s in samples:
        try:
            with Fasta(s.fasta, cache_dir=cache) as fa:
                lens = fa.lengths()
        except Exception as exc:  # pysam raises several kinds on a malformed or unindexable FASTA
            raise AnalysisError(f"{s.label}: cannot read the FASTA index of {s.fasta}: {exc}") from exc
        seqids[s.label] = tuple(lens)
        lengths[s.label] = np.fromiter(lens.values(), dtype=np.int64, count=len(lens))
    return seqids, lengths


def _tops(path: Path, label: str, members: np.ndarray) -> dict[str, str]:
    """Each member eid of genome `label` (members: its eid column, None where it has none) -> the eid of its
    outermost annotated ancestor, by the `parent` links of its element table at `path`."""
    try:
        parents = {e.eid: e.parent for e in read_elements(path)}
    except (ValueError, TypeError) as exc:
        raise AnalysisError(f"{path}: not an element table of ltrtrace ({exc})") from exc
    memo: dict[str, str] = {}
    top = {}
    for e in members:
        if e is None:
            continue
        if e not in parents:
            raise AnalysisError(f"{path}: {label}'s element {e}, a member in loci.tsv, is not in it: "
                                f"{_REMERGE}")
        top[e] = _outermost(e, parents, memo)
    return top


def _outermost(eid: str, parents: dict[str, str | None], memo: dict[str, str]) -> str:
    """The outermost annotated ancestor of `eid` (itself without a parent); a parent that is not annotated,
    or a cycle, ends the walk."""
    chain, x, seen = [], eid, {eid}
    while x not in memo:
        chain.append(x)
        p = parents.get(x)
        if p is None or p not in parents or p in seen:
            found = x
            break
        seen.add(p)
        x = p
    else:
        found = memo[x]
    for y in chain:
        memo[y] = found
    return found


# ---------------------------------------------------------------- each annotated genome's own calls


def _read_own_calls(pairs: dict[tuple[str, str], Path], labels: tuple[str, ...], annotated: list[str],
                    m: dict[str, np.ndarray], seqids: dict[str, tuple[str, ...]]
                    ) -> dict[str, dict[str, np.ndarray]]:
    """proj, ptier, pseq, pstart, pend and pbstrand: each annotated A's calls into every other
    genome, one pair file at a time, and A's own column."""
    n, G = m["state"].shape
    out: dict[str, dict[str, np.ndarray]] = {c: {} for c in ("proj", "ptier", "pseq", "pstart", "pend",
                                                             "pbstrand")}
    orders = {x: {name: j for j, name in enumerate(seqids[x])} for x in labels}
    for a in annotated:
        ga = labels.index(a)
        arr = {"proj": np.zeros((n, G), dtype=np.uint8), "ptier": np.zeros((n, G), dtype=np.uint8),
               "pseq": np.full((n, G), -1, dtype=np.int32), "pstart": np.zeros((n, G), dtype=np.int64),
               "pend": np.zeros((n, G), dtype=np.int64),
               "pbstrand": np.zeros((n, G), dtype=np.int8)}
        rows = {e: i for i, e in enumerate(m["eid"][:, ga].tolist()) if e is not None}
        mem = m["member"][:, ga]
        for gb, b in enumerate(labels):
            if gb == ga:
                continue
            path = pairs[(a, b)]
            _read_pair(path, b, rows, orders[b], *(arr[c][:, gb] for c in arr))
            lost = np.flatnonzero(mem & (arr["proj"][:, gb] == 0))
            if lost.size:
                some = ", ".join(m["eid"][lost[:5], ga])
                raise AnalysisError(f"{path}: {lost.size} member element(s) of {a} in loci.tsv have no call "
                                    f"({some}{' ...' if lost.size > 5 else ''}): {_REMERGE}")
        arr["proj"][mem, ga] = STATE_CODE["P"]
        arr["ptier"][mem, ga] = TIER_CODE["high"]
        for c, src in (("pseq", "seq"), ("pstart", "start"), ("pend", "end")):
            arr[c][mem, ga] = m[src][mem, ga]
        arr["pbstrand"][mem, ga] = 1  # a member's own site runs on its own strand
        for c, v in arr.items():
            out[c][a] = v
    return out


def _read_pair(path: Path, b: str, rows: dict[str, int], order: dict[str, int],
               state: np.ndarray, tier: np.ndarray, seq: np.ndarray, start: np.ndarray, end: np.ndarray,
               bstrand: np.ndarray) -> None:
    """Stream one pair file (A's calls into B; calls.records.stream_calls, each line split only as far as its
    b_strand) into the column views: the call of an element that loci.tsv lists as its locus's member (`rows`:
    eid -> row) fills that row; A's other elements are skipped."""
    try:
        for line, (eid, st_, tr_, sq, sb, se, sd_) in stream_calls(path, _CALL):
            i = rows.get(eid)
            if i is None:
                continue  # not its locus's first member element of A (multi_member): ignored
            st, tr, sd = STATE_CODE.get(st_), _CALL_TIER.get(tr_), _STRAND.get(sd_)
            if st is None or tr is None or sd is None:
                raise AnalysisError(f"{path}: line {line}: unknown state, tier or b_strand ({st_!r}, {tr_!r}, "
                                    f"{sd_!r})")
            if state[i]:
                raise AnalysisError(f"{path}: line {line}: {eid} is called more than once")
            state[i], tier[i] = st, tr
            if sq == ".":
                if st <= STATE_CODE["T"]:
                    raise AnalysisError(f"{path}: line {line}: a {st_} call without a site: not a pair file of "
                                        "ltrtrace genotype")
                continue
            j = order.get(sq)
            if j is None:
                raise AnalysisError(f"{path}: line {line}: {b}'s seqid {sq!r} is not in its FASTA: the samples "
                                    "sheet names another FASTA than the run's")
            try:
                s, e = int(sb), int(se)
            except ValueError:
                s = e = 0
            if not 1 <= s <= e:
                raise AnalysisError(f"{path}: line {line}: site {sq}:{sb}-{se} is not two integers with 1 <= start "
                                    "<= end")
            seq[i], start[i], end[i], bstrand[i] = j, s, e, sd
    except AnalysisError:
        raise
    except ValueError as exc:  # the reader's: a header that is not CALL_COLUMNS, a short line, a corrupt file
        raise AnalysisError(str(exc)) from exc


# ---------------------------------------------------------------- statistics with a discovery rule


def check_trio(cat: Catalog, labels: tuple[str, ...], who: str) -> None:
    """AnalysisError unless `labels` are distinct ingroup genomes of the run."""
    for x in labels:
        cat.g(x)  # an unknown label raises, naming the run's genomes
    out = [x for x in labels if x in cat.outgroups]
    if out or len(set(labels)) != len(labels):
        raise AnalysisError(f"{who}: {', '.join(labels)} must be distinct ingroup genomes"
                            + (f" ({', '.join(out)}: outgroup)" if out else ""))


def first_member(cat: Catalog, disc: tuple[str, ...] | list[str]) -> np.ndarray:
    """gstar, (n,) int16: per row, the first genome of `disc`, in order, that is a member there; else -1."""
    gstar = np.full(len(cat.lid), -1, dtype=np.int16)
    for w in disc:
        g = cat.g(w)
        gstar[(gstar < 0) & cat.member[:, g]] = g
    return gstar


def site_blocks(cat: Catalog, gstar: np.ndarray, D: str, rows: np.ndarray, size: int) -> np.ndarray:
    """Jackknife block ids of the selected rows (catalog order): D's site in each row's discoverer call, in
    blocks of `size` bp (stats.resample.block_ids)."""
    seq, mid = cat.proj_site(gstar, D)
    seq, mid = seq[rows], mid[rows]
    if (seq < 0).any():  # D is present in every counted row's call, and a present call has a site
        raise ValueError(f"{np.count_nonzero(seq < 0)} counted locus/loci have no site of {D} in their "
                         f"discoverer's call: D must be present at every counted locus")
    return block_ids(seq, mid, size)


def patterns(cat: Catalog, gstar: np.ndarray, R: str, S: str) -> tuple[np.ndarray, np.ndarray]:
    """(n,) masks from gstar's own calls: R present and S E (R+D, when D is present), and the mirror."""
    pr, ps = cat.pcls_by(gstar, R), cat.pcls_by(gstar, S)
    return (pr == 1) & (ps == 0), (ps == 1) & (pr == 0)
