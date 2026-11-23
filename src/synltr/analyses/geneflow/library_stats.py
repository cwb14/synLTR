"""The library screen's statistics of a trio, from the library route's pair files; the every-pair evenness check; and
the gate that holds back the statistics uneven annotation can fake.

Input: a finished screen (`LibraryRun`; routes.library.run.screen_pairs): for each ordered genome pair X, Y the call
file pairs/<X>__<Y>.calls.tsv.gz and the divergence file pairs/<X>__<Y>.div.tsv.gz (opened by name from the sheet's
pairs), each genome's records (routes.library.Records) and the k-mer stage's single-copy tallies
(routes.library.run.read_tallies). Each call file becomes arrays over every record of X in library order
(`pair_calls`): the ortholog's record index in Y (a P call's divergence row names it by eid among Y's placed records),
the call's state, its room (b_ins) and the divergence n_diff / n_sites (NaN without an aligned column). A record whose
name holds no coordinates has no call: it lies in no collinear block, so it has no ortholog, no room and no
divergence and is neither present nor absent; it still counts in denominators, in the jackknife blocks (its own
sequence) and in the placed span (`block_bp`).

Absences (`ABSENCES`): each statistic that reads absences counts those it was validated on, the screen's inferred
absences (A; calls.states.INFERRED_ABSENT). A statistic validated on observed empty sites (E) never counts an A: when
the calls hold absences it does not accept and none it does, its rows are NaN and their note says why.

Roles are given: R the recipient, S its sister, D the donor. Rows, layer "library", scope "R,S;D" ("R,S;D;D2" for a
second donor's), one row per statistic, an estimate's se and z on its own row:
  sharing             single_copy_kmers  "X+Y" (labels sorted) per pair of the trio: single-copy k-mers the pair
                                         shares and the third genome lacks
  D_kmer              all                D = (nA - nB) / (nA + nB) with se and z; nA, nB: single-copy k-mers D's
                                         records share with R only, with S only; nblocks (blocks on D's records)
  donor_resolved      "<D> elements inside synteny"  resolved_R, resolved_S, resolved_both: the share of D's records
                                         resolved (present or absent) in R, in S, in both; n, D's records
  D_ins               all                D over insertions: nA R+D copies absent from S, nB S+D copies absent from R,
                                         each counted only when the donor copy is also absent from the other sister
                                         (symmetric resolution); blocks on the donor copies
                      absence_certain    the same, where the absent sister's syntenic interval is shorter than
                                         ROOM_FRAC of the element (too short to hold it)
  D_ins_by_age        "K2P <bin>"        by the donor copy's K2P (SCREEN_AGE_EDGES; an unknown K2P sorts into the last
                                         bin), bins with >= SCREEN_AGE_MIN insertions
  D_ins_by_chrom      <D's sequence>     per donor sequence with >= SCREEN_CHROM_MIN insertions
  completeness_ratio  "<R>/<S> on <D> K2P><k>"  ratio (se): of D's old insertions (K2P above its OLD_QUANTILE
                                         quantile k) resolved in both sisters, those present in R over those present in
                                         S (at least 1 below); n; only with > COMPLETENESS_MIN of them; a value only
  R_share_D_by_chrom  <R's sequence>     of R's single-copy k-mers and insertions shared with exactly one other genome,
                                         the shares with D: kmer_share, kmers_D, kmers_S, ins_share, ins_D, ins_S;
                                         sequences with >= SCREEN_KMER_CHROM_MIN such k-mers
  D_kmer_by_chrom     <D's sequence>     D_kmer per donor sequence with >= SCREEN_KMER_CHROM_MIN informative k-mers
  identity_spectrum   "d in [a,b)"       per IDENTITY_BINS bin: n_target (R-D ortholog divergences whose donor copy is
                                         resolved in S), n_mirror (S-D, resolved in R), ratio_calibrated
  f_lower_bound       identity_spectrum  f: the target's excess over the calibrated mirror per bin, summed, over the
                                         target's orthologs: a lower bound on R's admixed fraction; n_target, n_mirror;
                                         calibration (mirror over target counts above the pooled median divergence);
                                         median_d; n_target_recent, n_mirror_recent (below IDENTITY_RECENT, the mirror
                                         calibrated)
  D3_triples          all                D from ortholog triples (R's, S's and D's copies, both sisters projected onto
                                         the donor's): (nA - nB) / (nA + nB) with nA the S-D and nB the R-D mismatches,
                                         positive when R is closer to D; triples; blocks on R's copies
  D_self_outgroup     LTR sites          D at the donor element's LTR sites polarized by its other LTR (copies whose own
                                         LTRs pass the own-LTR check): nA ABBA (R shares D's derived allele), nB BABA;
                                         BBAA; blocks on R's copies
  f_ratio_ins         "second donor <D2>"  alpha (se) = ([R+D, not S] - [S+D, not R]) / ([D2+D, not S] - [S+D, not
                                         D2]) over D's records resolved in R, S and D2; num, den, nblocks; NaN with the
                                         reason when den is 0, its se NaN when a jackknife replicate's den is 0
  f_ratio_identity    "second donor <D2>"  alpha = excess_R / excess_D2: the identity-spectrum excesses of R's and of
                                         D2's orthologs of D's records over S's
A second donor D2 is each other ingroup genome of the sheet (sheet order; never an outgroup) closer to D, and D to
it, than either is to R or S by 1:1 orthologs (P calls). Then the evenness rows of the trio's pairs (`evenness`).

The jackknife is stats.resample.jackknife (delete-one block): blocks of `block_bp` bp named <sequence>:<index>, taken
in the text order of those names (it fixes the order of the pseudovalues, and so the last bits of each se). The
D-form statistics are stats.resample.d_ratio on summed counts; the completeness ratio's denominator is at least 1.

Annotation and the gate. D_ins and D_kmer read annotation completeness as gene flow: in simulated nulls with one sister
10 points less completely annotated they gave false signals (|Z| > 3) in 9 and 5 of 10, and no check separates that
gap reliably from even annotation; every row of theirs is flagged "held: annotation-sensitive" and gives no verdict.
D3 and D_self held the null up to that 10-point gap but not at 30 points (4 and 10 false signals in 10): they are
flagged "held: uneven annotation" unless the every-pair evenness check passes (`gate`, `hold`). The completeness
ratio compares the sisters only and gates nothing: it is a value. Each headline row (D_kmer, D_ins, D3, D_self)
carries its false-signal rates (base.constants.FALSE_SIGNALS) in its note; `findings` prints them beside each value
and gives verdicts only for the statistics not held.
"""
from __future__ import annotations

import hashlib
import logging
import math
import re
import time
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np

from synltr.analyses.geneflow.route import RouteProps, require
from synltr.base.checkpoint import atomic_path, clear_done, code_fingerprint, is_done, mark_done, marker, stage_key
from synltr.base.constants import (
    COMPLETENESS_MIN,
    COMPLETENESS_Z,
    EVEN_BLOCK_BP,
    EVEN_CUTOFF,
    EVEN_MEASURED,
    FALSE_SIGNALS,
    IDENTITY_BINS,
    IDENTITY_RECENT,
    OLD_QUANTILE,
    RESOLVED_BOTH_WARN,
    RESOLVED_GAP_WARN,
    ROOM_FRAC,
    SCREEN_AGE_EDGES,
    SCREEN_AGE_LABELS,
    SCREEN_AGE_MIN,
    SCREEN_CHROM_MIN,
    SCREEN_KMER_CHROM_MIN,
    Z_SIG,
)
from synltr.base.log import AnalysisError
from synltr.calls.records import read_calls_meta, read_divergences, stream_calls
from synltr.calls.states import INFERRED_ABSENT, OBSERVED_ABSENT, PAIR_STATES, call_class
from synltr.io.elements import ElementTable
from synltr.io.library import read_library
from synltr.io.samples import PAIR_SEP, Sample
from synltr.routes.library import Records
from synltr.routes.library.clock import self_ltr_ok, triple_counts
from synltr.routes.library.run import library_digest, pair_files, read_tallies
from synltr.stats.report import Findings, StatRow, fmt
from synltr.stats.resample import Jk, block_ids, block_size, d_ratio, jackknife

log = logging.getLogger(__name__)

LAYER = "library"
SECTION = "Trio tests"  # the findings section (stats.report.SECTIONS)
TRIPLES = Path("work") / "geneflow" / "triples"  # under the screen's folder: each trio's alignments (`triples`)
HELD_SENSITIVE = "held: annotation-sensitive"
HELD_UNEVEN = "held: uneven annotation"
SENSITIVE = ("D_kmer", "D_kmer_by_chrom", "D_ins", "D_ins_by_age", "D_ins_by_chrom")  # always held
EVEN_GATED = ("D3_triples", "D_self_outgroup")  # held unless the every-pair evenness check passes
# the absences each statistic that reads them accepts: the screen's inferred absences, on which each was validated
ABSENCES: dict[str, frozenset[str]] = dict.fromkeys(
    ("donor_resolved", "D_ins", "completeness_ratio", "R_share_D_by_chrom", "identity_spectrum", "f_ratio_ins",
     "evenness", "painting"), INFERRED_ABSENT)
_STATE = {s: i for i, s in enumerate(PAIR_STATES)}
NO_CALL = len(PAIR_STATES)  # the state code of a record without a call
_ANY_ABSENT = OBSERVED_ABSENT | INFERRED_ABSENT
_HELD = re.compile(f"^(?:{re.escape(HELD_SENSITIVE)}|{re.escape(HELD_UNEVEN)})(?:; |$)")


@dataclass(frozen=True, eq=False)
class LibraryRun:
    """A finished library screen: its output folder (routes.library.run.screen_pairs), the samples sheet's samples
    in sheet order (each with its library), the threads to use, and label -> Records of the libraries already read
    (screen_pairs(..., records=) fills it); a library missing there is read when a statistic needs it, once. Each
    direction's pair calls (`pair_calls`) and each genome's sequence names (`sequences`) are read once per run: build
    it after the screen has finished, as its caches assume the screen's files no longer change."""

    outdir: Path
    samples: Sequence[Sample]
    threads: int
    records: dict[str, Records] = field(default_factory=dict)
    _orthologs: dict[frozenset[str], int] = field(default_factory=dict, init=False, repr=False)
    _calls: dict[tuple[str, str], PairCalls] = field(default_factory=dict, init=False, repr=False)
    _sequences: dict[str, tuple[tuple[str, ...], np.ndarray]] = field(default_factory=dict, init=False, repr=False)

    @property
    def labels(self) -> list[str]:
        return [s.label for s in self.samples]

    def rec(self, label: str) -> Records:
        """Genome `label`'s Records, read from its library the first time."""
        if label not in self.records:
            sample = next((s for s in self.samples if s.label == label), None)
            if sample is None or not sample.library:
                raise AnalysisError(f"library statistics: {label} is not a genome of the sheet with a library")
            self.records[label] = Records.of(label, read_library(sample.library, label, threads=self.threads))
        return self.records[label]

    def sequences(self, label: str) -> tuple[tuple[str, ...], np.ndarray]:
        """Genome `label`'s sequence names in text order (each record's chromosome; a record whose name holds no
        coordinates is its own sequence) and each record's sequence code, the rank of its name there (int64)."""
        if label not in self._sequences:
            names, codes = np.unique(self.rec(label).chrom.astype(str), return_inverse=True)
            self._sequences[label] = (tuple(names.tolist()), codes.astype(np.int64))
        return self._sequences[label]

    @property
    def seqids(self) -> Mapping[str, tuple[str, ...]]:
        """label -> its sequence names in text order (`sequences`), read on first use: the names a sequence code
        stands for in the painting's files and the inverted blocks."""
        return _Names(self)


class _Names(Mapping):
    """LibraryRun.seqids: label -> the genome's sequence names, from LibraryRun.sequences."""

    def __init__(self, run: LibraryRun) -> None:
        self._run = run

    def __getitem__(self, label: str) -> tuple[str, ...]:
        if label not in self._run.labels:
            raise KeyError(label)
        return self._run.sequences(label)[0]

    def __iter__(self):
        return iter(self._run.labels)

    def __len__(self) -> int:
        return len(self._run.labels)


# ---------------------------------------------------------------- the pair files


@dataclass(frozen=True, eq=False)
class PairCalls:
    """X's records called in Y (`pair_calls`), arrays over every record of X in library order."""

    x: str
    y: str
    state: np.ndarray  # uint8: the call's state, its index in calls.states.PAIR_STATES; NO_CALL without a call
    orth: np.ndarray  # int64: Y's record index of a P call's ortholog; -1 otherwise
    room: np.ndarray  # int64: a P or A call's room (b_ins: Y's span between the flanking anchors); -1 otherwise
    div: np.ndarray  # float64: a P call's n_diff / n_sites; NaN otherwise and without an aligned column

    def status(self, absent: frozenset[str] = INFERRED_ABSENT) -> np.ndarray:
        """int8 per record: 2 present, 1 absent by an absence in `absent`, 0 otherwise (`classes`)."""
        lut = np.zeros(NO_CALL + 1, np.int8)
        for s, c in classes(absent).items():
            lut[_STATE[s]] = c
        return lut[self.state]

    def absences(self) -> set[str]:
        """The absence states (observed or inferred) among the calls."""
        held = {PAIR_STATES[i] for i in np.unique(self.state).tolist() if i < NO_CALL}
        return {s for s in held if call_class(s, absent=_ANY_ABSENT) == "absent"}


def classes(absent: frozenset[str]) -> dict[str, int]:
    """Each pair state's class code under the absences `absent` (calls.states.call_class): 2 present, 1 absent, 0
    missing."""
    code = {"present": 2, "absent": 1, "missing": 0}
    return {s: code[call_class(s, absent=absent)] for s in PAIR_STATES}


def _finished(path: Path) -> None:
    """An AnalysisError naming `path` unless it is a finished library screen's file."""
    if not path.is_file():
        raise AnalysisError(f"{path}: missing: run the library screen first")
    try:
        meta = read_calls_meta(path)
    except ValueError:
        raise AnalysisError(f"{path}: unfinished (no done marker): rerun the library screen") from None
    if meta.get("engine") != "library":
        raise AnalysisError(f"{path}: written by the {meta.get('engine')} engine, not by the library screen")


def _pair(run: LibraryRun, x: str, y: str) -> tuple[Path, Path]:
    """(calls, divergences) of x's records in y, by name from the sheet's pair (its first genome first)."""
    labels = run.labels
    for g in (x, y):
        if g not in labels:
            raise AnalysisError(f"library statistics: {g} is not a genome of the sheet ({', '.join(labels)})")
    first = labels.index(x) < labels.index(y)
    files = pair_files(run.outdir, *((x, y) if first else (y, x)))
    return (files[0], files[2]) if first else (files[1], files[3])


def pair_calls(run: LibraryRun, x: str, y: str) -> PairCalls:
    """x's records called in y, over every record of x in library order (module docstring), from the screen's call
    and divergence files, read once per run (cached in `run`; the arrays are read-only). Each P call's ortholog is
    found among y's placed records by its divergence row's b_eid; a
    P call whose b_eid is "." (an ortholog whose name holds no coordinates) or names no placed record of y, a file
    that is missing or unfinished, and files that do not match the libraries are AnalysisErrors naming the pair or
    the file."""
    if (x, y) in run._calls:
        return run._calls[x, y]
    calls, divs = _pair(run, x, y)
    for path in (calls, divs):
        _finished(path)
    rx, ry = run.rec(x), run.rec(y)
    at = {e: i for i, e in enumerate(rx.lib.elements.column("eid").tolist()) if e is not None}
    state = np.full(rx.n, NO_CALL, np.uint8)
    room = np.full(rx.n, -1, np.int64)
    for line, (eid, st, b_ins) in stream_calls(calls, ("eid", "state", "b_ins")):
        i, code = at.get(eid), _STATE.get(st)
        if i is None or code is None:
            raise AnalysisError(f"{calls}: line {line}: {eid} {st} is no call of a record of {x}: the screen's files "
                                "do not match its libraries (rerun the screen)")
        state[i] = code
        if b_ins != ".":
            room[i] = int(b_ins)
    in_y = {e: j for j, e in enumerate(ry.lib.elements.column("eid").tolist()) if e is not None}
    orth = np.full(rx.n, -1, np.int64)
    div = np.full(rx.n, np.nan)
    p = state == _STATE["P"]
    pair = f"{x}{PAIR_SEP}{y}"
    for row in read_divergences(divs):
        i = at.get(row["eid"])
        if i is not None and orth[i] >= 0:
            raise AnalysisError(f"{divs}: a duplicate divergence row of {row['eid']}: the screen's files do not match "
                                "(rerun the screen)")
        if i is None or not p[i]:
            raise AnalysisError(f"{divs}: the divergence row of {row['eid']} matches no P call of {calls.name}: the "
                                "screen's files do not match (rerun the screen)")
        if row["b_eid"] is None:
            raise AnalysisError(f"pair {pair}: the ortholog of P call {row['eid']} has no coordinates (b_eid '.' in "
                                f"{divs.name}): its record in {y} is unknown, so no statistic can use the pair")
        j = in_y.get(row["b_eid"])
        if j is None:
            raise AnalysisError(f"pair {pair}: P call {row['eid']}: its ortholog {row['b_eid']} is no placed record "
                                f"of {y}'s library")
        orth[i] = j
        div[i] = row["n_diff"] / row["n_sites"] if row["n_sites"] else math.nan
    if (p & (orth < 0)).any():
        raise AnalysisError(f"{divs}: {int((p & (orth < 0)).sum())} P call(s) of {calls.name} have no divergence row: "
                            "the screen's files do not match (rerun the screen)")
    for a in (state, orth, room, div):  # shared by every trio of the run: read-only
        a.setflags(write=False)
    run._calls[x, y] = PairCalls(x, y, state, orth, room, div)
    return run._calls[x, y]


def orthologs(run: LibraryRun, x: str, y: str) -> int:
    """The 1:1 orthologs of genomes x and y: the P calls of their pair (either direction holds as many)."""
    key = frozenset((x, y))
    if key not in run._orthologs:
        calls = _pair(run, x, y)[0]
        _finished(calls)
        run._orthologs[key] = sum(st == "P" for _, (st,) in stream_calls(calls, ("state",)))
    return run._orthologs[key]


def second_donors(run: LibraryRun, R: str, S: str, D: str) -> list[str]:
    """The ingroup genomes of the sheet outside the trio (sheet order; a genome with role outgroup never is one, as no
    trio holds one) closer to D, and D to them, than either is to R or S: more 1:1 orthologs with D than D has with R
    or S, and than they have with R or S."""
    out = []
    for g in (s.label for s in run.samples if s.role != "outgroup"):
        if g in (R, S, D):
            continue
        nd = orthologs(run, D, g)
        if all(nd > orthologs(run, x, y) for x in (D, g) for y in (R, S)):
            out.append(g)
    return out


# ---------------------------------------------------------------- jackknife blocks


def placed_span(chrom: np.ndarray, end: np.ndarray) -> int:
    """A genome's placed span: the sum over its sequences of their last element's end (a record whose name holds no
    coordinates is its own sequence, ending at its length)."""
    if not len(chrom):
        return 0
    _, inv = np.unique(np.asarray(chrom, object), return_inverse=True)
    last = np.zeros(int(inv.max()) + 1, np.int64)
    np.maximum.at(last, inv, np.asarray(end, np.int64))
    return int(last.sum())


def block_bp(run: LibraryRun, props: RouteProps) -> int:
    """The jackknife block length of a screen under its route's block rule, "placed": BLOCK_BP, or 1/BLOCK_FRAC of the
    smallest placed span among every genome of the sheet, at least BLOCK_MIN (stats.resample.block_size)."""
    require("jackknife block rule", props.block, ("placed",))
    return block_size([placed_span(run.rec(g).chrom, run.rec(g).end) for g in run.labels])


def _blocks(rec: Records, idx, size: int) -> np.ndarray:
    """The jackknife blocks of records `idx` of `rec`: <sequence>:<midpoint // size>, as int64 ranks of those names
    in text order. The blocks are stats.resample.block_ids' blocks, but not its ids: the jackknife sums its
    pseudovalues in id order, and LTRmix's text order is what keeps the screen's standard errors bit-identical to
    LTRmix's."""
    idx = np.asarray(idx, np.int64)
    if not idx.size:
        return np.zeros(0, np.int64)
    names = np.array([f"{c}:{p // size}" for c, p in zip(rec.chrom[idx].tolist(), rec.mid[idx].tolist())], object)
    return np.unique(names, return_inverse=True)[1].astype(np.int64)


def chrom_key(c) -> tuple:
    """Sequences with a trailing number first, by that number, then the others; ties by name."""
    m = re.search(r"(\d+)$", str(c))
    return (0, int(m.group(1)), str(c)) if m else (1, 0, str(c))


# ---------------------------------------------------------------- the statistics


def age_bins(k) -> np.ndarray:
    """D_ins's age stratum of each K2P: 0 for 0, then the right-closed bins of SCREEN_AGE_EDGES ((0, 0.002], ...), the
    last for K2P above 0.1 or unknown (SCREEN_AGE_LABELS order)."""
    return np.digitize(np.asarray(k, float), SCREEN_AGE_EDGES[1:], right=True)


def absence_certain(room, length) -> np.ndarray:
    """An absence is certain where its syntenic interval (room, -1 without one) is shorter than ROOM_FRAC of the
    element: too short to hold it."""
    room = np.asarray(room)
    return (room >= 0) & (room < ROOM_FRAC * np.asarray(length))


def completeness_set(k2p, st_r, st_s) -> tuple[float, np.ndarray]:
    """The completeness ratio's old insertions of the donor: k_old, the OLD_QUANTILE quantile of every record's finite
    K2P (inf without one), and the records with K2P above it resolved (status > 0) in both sisters; none unless more
    than COMPLETENESS_MIN."""
    k2p = np.asarray(k2p, float)
    k_old = float(np.nanquantile(k2p, OLD_QUANTILE)) if np.isfinite(k2p).any() else math.inf
    old = np.flatnonzero((k2p > k_old) & (np.asarray(st_r) > 0) & (np.asarray(st_s) > 0))
    return k_old, old if old.size > COMPLETENESS_MIN else old[:0]


def d_stat(a, b, blocks) -> Jk:
    """D = (sum a - sum b) / (sum a + sum b), with its delete-one block jackknife (NaN without a count)."""
    return jackknife(d_ratio(np.asarray(a), np.asarray(b)), np.asarray(blocks))


def completeness(present_r, present_s, blocks) -> Jk:
    """The share present in R over the share present in S of the same old insertions: the count present in R over
    the count present in S (at least 1), with its delete-one block jackknife."""
    pr, ps = np.asarray(present_r, float), np.asarray(present_s, float)

    def stat(keep: np.ndarray) -> float:
        return float(pr[keep].sum()) / max(float(ps[keep].sum()), 1)
    return jackknife(stat, np.asarray(blocks))


@dataclass(frozen=True)
class FRatio:
    """An f-ratio with a second donor (`f_ratio`)."""

    ins: Jk  # alpha from insertion patterns, jackknifed
    num: int
    den: int
    why: str | None  # why alpha, or its se, is NaN (a zero denominator); None otherwise
    identity: float  # alpha from identity-spectrum excesses
    excess_r: float
    excess_2: float


def f_ratio(st_r, st_s, st_2, blocks, d_r, d_s, d_2) -> FRatio:
    """The f-ratio over the donor's records from their statuses (2 present, 1 absent, 0 unresolved) in R, S and the
    second donor D2 and their jackknife `blocks`: alpha = sum num / sum den over the records resolved in all three,
    num = [R present, S not] - [S present, R not], den = [D2 present, S not] - [S present, D2 not]. When den sums to
    0 (no record differs, or the two directions are equally many), alpha is NaN and so is its se (a standard error of
    an undefined estimate), with the reason; when deleting one jackknife block leaves den 0, that replicate has no
    value and the se is NaN, with the reason. And from their divergences to R's, S's and D2's orthologs (NaN without
    one): each one's histogram over IDENTITY_BINS per donor record, its excess over S's (clipped at 0) summed,
    alpha = excess_R / excess_D2 (NaN when excess_D2 is 0)."""
    st_r, st_s, st_2 = (np.asarray(s) for s in (st_r, st_s, st_2))
    ok = (st_r > 0) & (st_s > 0) & (st_2 > 0)
    xr, xs, x2 = st_r == 2, st_s == 2, st_2 == 2
    num = ((xr & ~xs).astype(np.int64) - (xs & ~xr))[ok]
    den = ((x2 & ~xs).astype(np.int64) - (xs & ~x2))[ok]

    def stat(keep: np.ndarray) -> float:
        d = den[keep].sum()
        return float(num[keep].sum() / d) if d else math.nan
    blocks = np.asarray(blocks)[ok]
    jk = jackknife(stat, blocks)
    total = int(den.sum())
    why = None if total else ("its denominator is 0: among the donor's records resolved in the recipient, the sister and "
                              "the second donor, those present in the second donor and not the sister are as many as "
                              "those present in the sister and not the second donor")
    if total and jk.nblocks:
        _, inv = np.unique(blocks, return_inverse=True)
        if (np.bincount(inv, den, jk.nblocks) == total).any():
            jk = Jk(jk.est, math.nan, math.nan, jk.nblocks)
            why = "deleting one jackknife block leaves the denominator 0"
    bins = np.asarray(IDENTITY_BINS)
    n = len(d_r)
    h = [np.histogram(d[np.isfinite(d)], bins)[0] / n for d in (np.asarray(d_r), np.asarray(d_s), np.asarray(d_2))]
    e_r, e_2 = np.clip(h[0] - h[1], 0, None).sum(), np.clip(h[2] - h[1], 0, None).sum()
    return FRatio(jk, int(num.sum()), total, why, float(e_r / e_2) if e_2 > 0 else math.nan, float(e_r), float(e_2))


def _ratio_note(base: str, f: FRatio) -> str:
    """The insertion f-ratio row's note: its definition, then what is NaN and why (alpha, or its se alone)."""
    if math.isnan(f.ins.est):
        return f"{base}: NaN, {f.why}" if f.why else f"{base}: NaN"
    return f"{base}; se NaN, {f.why}" if f.why else base


@dataclass(frozen=True)
class Spectrum:
    """The identity spectrum of target and mirror ortholog divergences (`spectrum`)."""

    n_target: np.ndarray  # int64 per IDENTITY_BINS bin
    n_mirror: np.ndarray
    ratio: np.ndarray  # float64 per bin: n_target x calibration / n_mirror (NaN without a mirror count)
    f: float  # the admixed-fraction lower bound
    calibration: float
    median: float
    recent_target: int  # target orthologs below IDENTITY_RECENT
    recent_mirror: float  # mirror orthologs below IDENTITY_RECENT, calibrated


def spectrum(t, m) -> Spectrum:
    """Target divergences `t` against mirror divergences `m`: the calibration (mirror over target counts above the
    pooled median: equal totals in the deep half), per-bin counts and calibrated ratios, and the lower bound on the
    admixed fraction: the target's per-bin excess over the calibrated mirror, clipped at 0, summed, over the target's
    count (at least 1). Empty inputs give NaN."""
    t, m = np.asarray(t, float), np.asarray(m, float)
    bins = np.asarray(IDENTITY_BINS)
    both = np.concatenate([t, m])
    med = np.median(both) if both.size else math.nan
    with np.errstate(divide="ignore", invalid="ignore"):
        cal = (m > med).sum() / max((t > med).sum(), 1)
        ht, _ = np.histogram(t, bins)
        hm, _ = np.histogram(m, bins)
        ratio = np.array([ht[i] * cal / hm[i] if hm[i] else np.nan for i in range(len(bins) - 1)], float)
        f_low = np.clip(ht - hm / cal, 0, None).sum() / max(len(t), 1)
        recent_m = (m < IDENTITY_RECENT).sum() / cal
    return Spectrum(ht, hm, ratio, float(f_low), float(cal), float(med), int((t < IDENTITY_RECENT).sum()),
                    float(recent_m))


# ---------------------------------------------------------------- one trio


def rates(test: str) -> str:
    """The false-signal rates of a headline statistic, for its note and findings."""
    return "false signals (|Z| > 3) in simulated nulls: " + ", ".join(
        f"{condition} {n}/{of}" for condition, (n, of) in FALSE_SIGNALS[test].items())


def _join(*notes: str) -> str:
    return "; ".join(n for n in notes if n and n != ".") or "."


def _jk_rows(test: str, scope: str, stratum: str, est: float, jk: Jk, z: float, n_a: int, n_b: int, notes: tuple,
             where: str, size: int) -> list[StatRow]:
    """A D statistic's rows: D (se, z), nA, nB and nblocks; notes = (D's, nA's, nB's)."""
    return [StatRow(LAYER, test, scope, stratum, "D", float(est), float(jk.se), float(z), note=notes[0]),
            StatRow(LAYER, test, scope, stratum, "nA", int(n_a), note=notes[1]),
            StatRow(LAYER, test, scope, stratum, "nB", int(n_b), note=notes[2]),
            StatRow(LAYER, test, scope, stratum, "nblocks", jk.nblocks,
                    note=f"jackknife blocks of {size:,} bp on {where}")]


def _refused(rows: list[StatRow], why: str | None) -> list[StatRow]:
    """The rows unchanged, or, with `why`, every value NaN and the note `why`."""
    if why is None:
        return rows
    return [replace(r, value=math.nan, se=math.nan, z=math.nan, note=why) for r in rows]


def _unaccepted(name: str, *calls: PairCalls) -> str | None:
    """Why statistic `name` has no absence to count (ABSENCES): the calls hold absences, none of a kind it accepts;
    None otherwise."""
    return _why_unaccepted(name, set().union(*(c.absences() for c in calls)))


def _why_unaccepted(name: str, seen: set[str]) -> str | None:
    """`_unaccepted` from the absence states `seen` among the calls."""
    accepted = ABSENCES[name]
    if not seen or seen & accepted:
        return None
    why = (f"no absence {name} accepts ({', '.join(sorted(accepted))}) among the calls, which hold "
           f"{', '.join(sorted(seen))}")
    if seen <= INFERRED_ABSENT and accepted <= OBSERVED_ABSENT:
        why += ": inferred absences (A) never count as observed empty sites (E)"
    return why + ": not computed"


@dataclass
class _Trio:
    run: LibraryRun
    R: str
    S: str
    D: str
    size: int

    def __post_init__(self) -> None:
        self.L = {g: self.run.rec(g) for g in (self.R, self.S, self.D)}
        self.calls = {(x, y): pair_calls(self.run, x, y) for x in self.L for y in self.L if x != y}
        self.order = [g for g in self.run.labels if g in self.L]
        self.tallies = {x: read_tallies(self.run.outdir, tuple(self.order), x) for x in self.order}
        self.scope = f"{self.R},{self.S};{self.D}"
        self.dblocks = _blocks(self.L[self.D], np.arange(self.L[self.D].n), self.size)

    def c(self, x: str, y: str) -> PairCalls:
        return self.calls[x, y]

    def st(self, name: str, x: str, y: str) -> np.ndarray:
        return self.calls[x, y].status(ABSENCES[name])

    # ---- single-copy k-mers
    def kmers(self) -> list[StatRow]:
        R, S, D, scope = self.R, self.S, self.D, self.scope
        share = sharing(self.order, self.tallies)
        rows = [StatRow(LAYER, "sharing", scope, "single_copy_kmers", "+".join(sorted(k)), v,
                        note=f"single-copy k-mers {' and '.join(sorted(k))} share and "
                             f"{next(g for g in self.order if g not in k)} lacks")
                for k, v in share.items()]
        a, b = self.tallies[D][(R,)], self.tallies[D][(S,)]
        jk = d_stat(a, b, self.dblocks)
        rows += _jk_rows("D_kmer", scope, "all", jk.est, jk, jk.z, a.sum(), b.sum(),
                         (f"D = (nA - nB) / (nA + nB): single-copy k-mers on {D}'s records; {rates('D_kmer')}",
                          f"single-copy k-mers {D} shares with {R} only",
                          f"single-copy k-mers {D} shares with {S} only"),
                         f"{D}'s records", self.size)
        return rows

    # ---- insertion patterns
    def patterns(self) -> list[StatRow]:
        R, S, D, scope, size = self.R, self.S, self.D, self.scope, self.size
        lD, lR, lS = self.L[D], self.L[R], self.L[S]
        rd, sd, rs = self.c(R, D).orth, self.c(S, D).orth, self.c(R, S).orth
        dr_, ds_ = self.st("donor_resolved", D, R), self.st("donor_resolved", D, S)
        cr, cs, cb = float((dr_ > 0).mean()), float((ds_ > 0).mean()), float(((dr_ > 0) & (ds_ > 0)).mean())
        resolved = f"of {D}'s records, the share resolved (present or absent) in"
        rows = _refused([StatRow(LAYER, "donor_resolved", scope, f"{D} elements inside synteny", stat, v, note=note)
                         for stat, v, note in (("resolved_R", cr, f"{resolved} {R}"),
                                               ("resolved_S", cs, f"{resolved} {S}"),
                                               ("resolved_both", cb, f"{resolved} both {R} and {S}"),
                                               ("n", int(lD.n), f"{D}'s records"))],
                        _unaccepted("donor_resolved", self.c(D, R), self.c(D, S)))
        # D_ins: R+D copies absent from S, S+D absent from R; the donor copy absent from the other sister too
        st_rs, st_sr = self.st("D_ins", R, S), self.st("D_ins", S, R)
        st_ds, st_dr = self.st("D_ins", D, S), self.st("D_ins", D, R)
        abba = np.flatnonzero((rd >= 0) & (st_rs == 1))
        baba = np.flatnonzero((sd >= 0) & (st_sr == 1))
        abba = abba[st_ds[rd[abba]] == 1]
        baba = baba[st_dr[sd[baba]] == 1]
        dA, dB = rd[abba], sd[baba]
        ins = []
        what = (f"{R}+{D} insertions absent from {S}", f"{S}+{D} insertions absent from {R}")
        on = f"{D}'s copies"

        def d_ins(stratum: str, ma: np.ndarray, mb: np.ndarray, test: str = "D_ins", note: str = "") -> list[StatRow]:
            x = np.concatenate([np.ones(int(ma.sum())), np.zeros(int(mb.sum()))])
            jk = d_stat(x, 1 - x, _blocks(lD, np.concatenate([dA[ma], dB[mb]]), size))
            return _jk_rows(test, scope, stratum, jk.est, jk, jk.z, ma.sum(), mb.sum(),
                            (note or "D = (nA - nB) / (nA + nB) over insertions", *what), on, size)

        every_a, every_b = np.ones(len(dA), bool), np.ones(len(dB), bool)
        ins += d_ins("all", every_a, every_b, note=f"D = (nA - nB) / (nA + nB) over insertions, each absence "
                                                   f"inferred (A) inside synteny; {rates('D_ins')}")
        room_r, room_s = self.c(R, S).room, self.c(S, R).room
        len_r, len_s = lR.end - lR.start + 1, lS.end - lS.start + 1
        ok_a, ok_b = absence_certain(room_r[abba], len_r[abba]), absence_certain(room_s[baba], len_s[baba])
        ins += d_ins("absence_certain", ok_a, ok_b, note=f"absences whose syntenic interval is shorter than "
                                                          f"{ROOM_FRAC:g} of the element")
        c_a, c_b = age_bins(lD.k2p[dA]), age_bins(lD.k2p[dB])
        for c, label in enumerate(SCREEN_AGE_LABELS):
            m1, m2 = c_a == c, c_b == c
            if m1.sum() + m2.sum() >= SCREEN_AGE_MIN:
                unknown = " or unknown" if c == len(SCREEN_AGE_LABELS) - 1 else ""
                ins += d_ins(f"K2P {label}", m1, m2, "D_ins_by_age",
                             f"insertions whose {D} copy has LTR K2P {label}{unknown}")
        ch_a, ch_b = lD.chrom[dA], lD.chrom[dB]
        for c in sorted(set(ch_a.tolist()) | set(ch_b.tolist()), key=chrom_key):
            m1, m2 = ch_a == c, ch_b == c
            if m1.sum() + m2.sum() >= SCREEN_CHROM_MIN:
                ins += d_ins(c, m1, m2, "D_ins_by_chrom", f"insertions whose {D} copy lies on {c}")
        rows += _refused(ins, _unaccepted("D_ins", self.c(R, S), self.c(S, R), self.c(D, S), self.c(D, R)))
        rows += self._completeness()
        rows += self._shares(rd, rs)
        rows += self._kmers_by_chrom()
        return rows

    def _completeness(self) -> list[StatRow]:
        R, S, D = self.R, self.S, self.D
        lD = self.L[D]
        dr_c, ds_c = self.st("completeness_ratio", D, R), self.st("completeness_ratio", D, S)
        k_old, old = completeness_set(lD.k2p, dr_c, ds_c)
        if not old.size:
            return []
        jk = completeness(dr_c[old] == 2, ds_c[old] == 2, _blocks(lD, old, self.size))
        stratum = f"{R}/{S} on {D} K2P>{k_old:.3g}"
        note = (f"of {D}'s old insertions (LTR K2P above its {OLD_QUANTILE:g} quantile) resolved in both sisters, the "
                f"count present in {R} over that in {S}; a value only, it gates nothing: it reads {R}'s annotation "
                f"completeness relative to {S}'s (or old polymorphisms carried by gene flow), and one-sided dropout "
                f"biases D_kmer and D_ins toward the more complete sister; D3 and D_self held the null in simulation "
                f"up to a 10-point completeness gap, not beyond (the every-pair evenness check gates them)")
        return _refused([StatRow(LAYER, "completeness_ratio", self.scope, stratum, "ratio", float(jk.est),
                                 float(jk.se), note=note),
                         StatRow(LAYER, "completeness_ratio", self.scope, stratum, "n", int(old.size),
                                 note=f"{D}'s old insertions resolved in both sisters")],
                        _unaccepted("completeness_ratio", self.c(D, R), self.c(D, S)))

    def _shares(self, rd: np.ndarray, rs: np.ndarray) -> list[StatRow]:
        """R's single-copy k-mers and insertions shared with exactly one other genome, per R sequence."""
        R, S, D = self.R, self.S, self.D
        st_rs, st_rd = self.st("R_share_D_by_chrom", R, S), self.st("R_share_D_by_chrom", R, D)
        rd_only = ((rd >= 0) & (st_rs == 1)).astype(float)
        rs_only = ((rs >= 0) & (st_rd == 1)).astype(float)
        names, inv = np.unique(self.L[R].chrom, return_inverse=True)
        sums = [np.bincount(inv, w, len(names)) for w in (self.tallies[R][(D,)], self.tallies[R][(S,)], rd_only,
                                                           rs_only)]
        rows = []
        for k in sorted(range(len(names)), key=lambda k: chrom_key(names[k])):
            kd, ks, nd, ns = (int(s[k]) for s in sums)
            if kd + ks < SCREEN_KMER_CHROM_MIN:
                continue
            for stat, v, note in (
                    ("kmer_share", kd / (kd + ks), (f"of {R}'s single-copy k-mers shared with exactly one of {S} and "
                                                    f"{D}, the share with {D}")),
                    ("kmers_D", kd, f"single-copy k-mers {R} shares with {D} only"),
                    ("kmers_S", ks, f"single-copy k-mers {R} shares with {S} only"),
                    ("ins_share", nd / max(nd + ns, 1), (f"of {R}'s insertions shared with exactly one of {S} and "
                                                         f"{D}, the share with {D}")),
                    ("ins_D", nd, f"{R}+{D} insertions absent from {S}"),
                    ("ins_S", ns, f"{R}+{S} insertions absent from {D}")):
                rows.append(StatRow(LAYER, "R_share_D_by_chrom", self.scope, names[k], stat, v, note=note))
        return _refused(rows, _unaccepted("R_share_D_by_chrom", self.c(R, S), self.c(R, D)))

    def _kmers_by_chrom(self) -> list[StatRow]:
        R, S, D = self.R, self.S, self.D
        a, b = self.tallies[D][(R,)], self.tallies[D][(S,)]
        names, inv = np.unique(self.L[D].chrom, return_inverse=True)
        sa, sb = np.bincount(inv, a, len(names)), np.bincount(inv, b, len(names))
        rows = []
        for k in sorted(range(len(names)), key=lambda k: chrom_key(names[k])):
            if sa[k] + sb[k] < SCREEN_KMER_CHROM_MIN:
                continue
            m = inv == k
            jk = d_stat(a[m], b[m], self.dblocks[m])
            rows += _jk_rows("D_kmer_by_chrom", self.scope, names[k], jk.est, jk, jk.z, a[m].sum(), b[m].sum(),
                             (f"D_kmer on {D}'s records on {names[k]}", f"single-copy k-mers {D} shares with {R} only",
                              f"single-copy k-mers {D} shares with {S} only"), f"{D}'s records", self.size)
        return rows

    # ---- ortholog divergences: the identity spectrum, triples
    def divergences(self) -> list[StatRow]:
        R, S, D, scope, size = self.R, self.S, self.D, self.scope, self.size
        rd, sd, rs = self.c(R, D).orth, self.c(S, D).orth, self.c(R, S).orth
        d_rd, d_sd = self.c(R, D).div, self.c(S, D).div
        ds_, dr_ = self.st("identity_spectrum", D, S), self.st("identity_spectrum", D, R)
        okt = np.isfinite(d_rd) & (rd >= 0)
        okt[okt] = ds_[rd[okt]] > 0
        okm = np.isfinite(d_sd) & (sd >= 0)
        okm[okm] = dr_[sd[okm]] > 0
        sp = spectrum(d_rd[okt], d_sd[okm])
        bins = IDENTITY_BINS
        rows = []
        for i in range(len(bins) - 1):
            stratum = f"d in [{bins[i]:.2g},{bins[i + 1]:.2g})"
            rows += [StatRow(LAYER, "identity_spectrum", scope, stratum, "n_target", int(sp.n_target[i]),
                             note=f"{R}-{D} orthologs whose {D} copy is resolved in {S}"),
                     StatRow(LAYER, "identity_spectrum", scope, stratum, "n_mirror", int(sp.n_mirror[i]),
                             note=f"{S}-{D} orthologs whose {D} copy is resolved in {R}"),
                     StatRow(LAYER, "identity_spectrum", scope, stratum, "ratio_calibrated", float(sp.ratio[i]),
                             note="n_target x calibration / n_mirror")]
        rows += [StatRow(LAYER, "f_lower_bound", scope, "identity_spectrum", stat, v, note=note) for stat, v, note in (
            ("f", sp.f, (f"lower bound on {R}'s admixed fraction: the excess of {R}-{D} over calibrated {S}-{D} "
                         f"orthologs per bin, summed, over {R}'s {D}-orthologs")),
            ("n_target", int(okt.sum()), f"{R}-{D} orthologs"), ("n_mirror", int(okm.sum()), f"{S}-{D} orthologs"),
            ("calibration", sp.calibration, "mirror over target orthologs above the pooled median divergence"),
            ("median_d", sp.median, "the pooled median divergence"),
            ("n_target_recent", sp.recent_target, f"{R}-{D} orthologs below {IDENTITY_RECENT:g} divergence"),
            ("n_mirror_recent", sp.recent_mirror, f"{S}-{D} orthologs below {IDENTITY_RECENT:g}, calibrated"))]
        rows = _refused(rows, _unaccepted("identity_spectrum", self.c(D, S), self.c(D, R)))
        # triples: the R-D, R-S and S-D orthologs agree
        tri = np.flatnonzero((rd >= 0) & (rs >= 0))
        tri = tri[sd[rs[tri]] == rd[tri]]
        t0 = time.perf_counter()
        self_ok, res, reused = triples(self.run, R, S, D, np.stack([rd[tri], tri, rs[tri]], 1))
        log.debug("%s: %d ortholog triples %s in %.1f s; own-LTR outgroup usable for %d", scope, len(tri),
                  "reused" if reused else "aligned", time.perf_counter() - t0, int(self_ok[rd[tri]].sum()))
        blk = _blocks(self.L[R], tri, size)
        d3 = d_stat(res[:, 1], res[:, 2], blk)
        rows += _jk_rows("D3_triples", scope, "all", -d3.est, d3, -d3.z if math.isfinite(d3.z) else math.nan,
                         res[:, 2].sum(), res[:, 1].sum(),
                         ((f"D = (nA - nB) / (nA + nB) over the triples' mismatches, positive when {R} is closer "
                           f"to {D}; {rates('D3_triples')}"), f"{S}-{D} mismatches", f"{R}-{D} mismatches"),
                         f"{R}'s copies", size)
        rows.append(StatRow(LAYER, "D3_triples", scope, "all", "triples", len(tri),
                            note=f"ortholog triples ({R}, {S} and {D} copies called orthologous in each pair)"))
        ds = d_stat(res[:, 4], res[:, 5], blk)
        rows += _jk_rows("D_self_outgroup", scope, "LTR sites", ds.est, ds, ds.z, res[:, 4].sum(), res[:, 5].sum(),
                         ((f"D = (ABBA - BABA) / (ABBA + BABA) at {D}'s LTR sites polarized by its other LTR; "
                           f"{rates('D_self_outgroup')}"), f"ABBA: {R} shares {D}'s derived allele",
                          f"BABA: {S} shares {D}'s derived allele"), f"{R}'s copies", size)
        rows.append(StatRow(LAYER, "D_self_outgroup", scope, "LTR sites", "BBAA", int(res[:, 6].sum()),
                            note=f"sites where {R} and {S} share an allele {D} lacks"))
        return rows

    # ---- a second donor
    def fratio(self, d2: str) -> list[StatRow]:
        R, S, D = self.R, self.S, self.D
        c2 = pair_calls(self.run, D, d2)
        name = "f_ratio_ins"
        f = f_ratio(self.st(name, D, R), self.st(name, D, S), c2.status(ABSENCES[name]), self.dblocks,
                    self.c(D, R).div, self.c(D, S).div, c2.div)
        scope, stratum = f"{self.scope};{d2}", f"second donor {d2}"
        alpha = _ratio_note(f"([{R}+{D}, not {S}] - [{S}+{D}, not {R}]) / ([{d2}+{D}, not {S}] - [{S}+{D}, not "
                            f"{d2}]) over {D}'s records resolved in {R}, {S} and {d2}", f)
        rows = _refused([StatRow(LAYER, name, scope, stratum, "alpha", float(f.ins.est), float(f.ins.se), note=alpha),
                         StatRow(LAYER, name, scope, stratum, "num", f.num, note="the numerator, summed"),
                         StatRow(LAYER, name, scope, stratum, "den", f.den, note="the denominator, summed"),
                         StatRow(LAYER, name, scope, stratum, "nblocks", f.ins.nblocks,
                                 note=f"jackknife blocks of {self.size:,} bp on {D}'s records")],
                        _unaccepted(name, self.c(D, R), self.c(D, S), c2))
        return rows + [
            StatRow(LAYER, "f_ratio_identity", scope, stratum, "alpha", f.identity,
                    note="excess_R / excess_D2" + ("" if math.isfinite(f.identity) else ": NaN, excess_D2 is 0")),
            StatRow(LAYER, "f_ratio_identity", scope, stratum, "excess_R", f.excess_r,
                    note=f"the excess of {D}'s records' {R} orthologs over its {S} orthologs per divergence bin"),
            StatRow(LAYER, "f_ratio_identity", scope, stratum, "excess_D2", f.excess_2,
                    note=f"the excess of {D}'s records' {d2} orthologs over its {S} orthologs per divergence bin")]


def triples(run: LibraryRun, R: str, S: str, D: str, trip: np.ndarray) -> tuple[np.ndarray, np.ndarray, bool]:
    """The trio's alignment work: clock.self_ltr_ok over D's records at the triples' donor copies and
    clock.triple_counts of the triples `trip` (rows: D's, R's and S's record index), with whether they were reused.
    Checkpointed under TRIPLES (base.checkpoint), keyed by the three libraries' files, the triples, the roles and the
    code of the library route (routes.library and every module it imports: the record reader, the records' LTR
    bounds and the alignments, the code the screen's pair files are keyed on), so a rerun after a kill redoes no
    finished trio and an edit to that code realigns every trio."""
    trip = np.ascontiguousarray(trip, dtype=np.int64)
    path = Path(run.outdir) / TRIPLES / f"{PAIR_SEP.join((R, S, D))}.npz"
    spec = {s.label: s.library for s in run.samples}
    inputs = {g: library_digest(spec[g], run.outdir) for g in (R, S, D)}
    inputs["triples"] = hashlib.sha256(trip.tobytes()).hexdigest()
    key = stage_key("screen_triples", inputs, {"roles": [R, S, D], "n": len(trip)}, code_fingerprint("routes.library"))
    if is_done(marker(path), key):
        with np.load(path) as z:
            return z["self_ok"], z["counts"], True
    clear_done(marker(path))
    self_ok = self_ltr_ok(run.rec(D), trip[:, 0], run.threads)
    res = triple_counts(run.rec(D), run.rec(R), run.rec(S), trip, self_ok, run.threads)
    path.parent.mkdir(parents=True, exist_ok=True)
    with atomic_path(path) as tmp, open(tmp, "wb") as fh:
        np.savez(fh, self_ok=self_ok, counts=res)
    mark_done(marker(path), key, [path])
    return self_ok, res, False


def sharing(order: Sequence[str], tallies: dict[str, dict[tuple[str, ...], np.ndarray]]) -> dict[frozenset[str], int]:
    """The single-copy k-mers each pair of the trio `order` shares and the third genome lacks, from each member's
    tallies (read_tallies), in the order {a,b}, {a,c}, {b,c}: each member in turn writes its pairs' counts, the last
    one's kept (the counts are symmetric)."""
    share: dict[frozenset[str], int] = {}
    for x in order:
        for o in (y for y in order if y != x):
            share[frozenset((x, o))] = int(tallies[x][(o,)].sum())
    return share


@dataclass(frozen=True)
class Roles:
    """A trio's roles by the single-copy k-mers (`roles`): recipient R, sister S, donor D; the sharing counts
    (`sharing`); D_kmer on D's records, nA k-mers D shares with R only and nB with S only; `tie`: the top sharing count
    was tied, so the first such pair in the order {a,b}, {a,c}, {b,c} was taken."""

    R: str
    S: str
    D: str
    share: dict[frozenset[str], int]
    d: Jk
    n_a: int
    n_b: int
    tie: bool

    @property
    def source(self) -> str:
        """The roles row's source: "data", or "sheet order" when a tie or D_kmer (nan or 0) left the choice to the
        sheet's order."""
        return "sheet order" if self.tie or not (math.isfinite(self.d.est) and self.d.est != 0) else "data"


def pick_roles(order: Sequence[str], tallies: dict[str, dict[tuple[str, ...], np.ndarray]],
               blocks: Callable[[str], np.ndarray]) -> Roles:
    """The roles of the trio `order` (sheet order) from each member's single-copy tallies and the jackknife blocks of
    a genome's records (`blocks(label)`): the sisters are the pair sharing the most single-copy k-mers the third
    lacks (`sharing`; a tie keeps the first pair in that order), in sheet order, and the donor the third; D_kmer =
    (nA - nB) / (nA + nB) on the donor's records, nA its k-mers shared with the first sister only, nB with the
    second; the recipient is the first sister when D_kmer >= 0, else the second (a nan D_kmer gives the second
    sister), with D_kmer taken from the recipient's side (>= 0)."""
    share = sharing(order, tallies)
    top = max(share.values())
    pair = next(k for k, v in share.items() if v == top)
    sis = tuple(sorted(pair, key=list(order).index))
    don = next(x for x in order if x not in sis)
    a, b = tallies[don][(sis[0],)], tallies[don][(sis[1],)]
    blk = blocks(don)
    d = d_stat(a, b, blk)
    R, S = (sis[0], sis[1]) if d.est >= 0 else (sis[1], sis[0])
    if d.est < 0:
        a, b = b, a
        d = d_stat(a, b, blk)
    return Roles(R, S, don, share, d, int(a.sum()), int(b.sum()), sum(v == top for v in share.values()) > 1)


def roles(lib_run: LibraryRun, genomes: Sequence[str], props: RouteProps) -> Roles:
    """The roles of three genomes of the screen by the route's roles rule, "single_copy_kmers" (`pick_roles`), from
    the k-mer stage's tallies (routes.library.run.read_tallies) and blocks of `block_bp` on the donor's records. An
    AnalysisError for genomes that are not three distinct genomes of the sheet."""
    require("roles", props.roles, ("single_copy_kmers",))
    order = [g for g in lib_run.labels if g in genomes]
    if len(order) != 3 or len(set(genomes)) != 3:
        raise AnalysisError(f"roles: the trio must be three distinct genomes of the sheet ({', '.join(lib_run.labels)}),"
                            f" not {', '.join(map(str, genomes)) or 'none'}")
    tallies = {x: read_tallies(lib_run.outdir, tuple(order), x) for x in order}
    size = block_bp(lib_run, props)
    return pick_roles(order, tallies, lambda g: _blocks(lib_run.rec(g), np.arange(lib_run.rec(g).n), size))


def role_rows(r: Roles, order: Sequence[str]) -> tuple[list[StatRow], Findings]:
    """The roles row (layer library, test roles, scope R,S;D, statistic source: "data" or "sheet order") and the
    findings lines of the roles `r` of the trio `order` (sheet order); a WARNING when a tie or D_kmer left the choice
    to the sheet's order."""
    f = Findings()
    trio = ",".join(order)
    sis = [g for g in order if g in (r.R, r.S)]
    pair = f"{sis[0]}+{sis[1]}"
    sign = f"D_kmer({r.R},{r.S};{r.D}) = {fmt(r.d.est, '+.4f')}"
    top = r.share[frozenset(sis)]
    note = (f"sister pair {pair} (the pair sharing the most single-copy k-mers the third genome lacks: {top}); "
            f"recipient {r.R} by the sign of {sign} on {r.D}'s records")
    if r.tie:
        tied = [f"{'+'.join(g for g in order if g in k)}" for k, v in r.share.items() if v == top]
        f.warn(SECTION, f"trio {trio}: the sister pair is a tie ({', '.join(tied)} each share {top} single-copy "
                        f"k-mers the third genome lacks); {pair}, the first in sheet order, is taken")
    if not math.isfinite(r.d.est):
        f.warn(SECTION, f"trio {trio}: the recipient is not identified ({sign}: no single-copy k-mer of {r.D} is "
                        f"shared with one sister only); {r.R}, the second of the sister pair in sheet order, is taken "
                        "as R")
    elif r.d.est == 0:
        f.warn(SECTION, f"trio {trio}: the recipient is not identified ({sign}: a tie); {r.R}, the first of the "
                        "sister pair in sheet order, is taken as R")
    f.add(SECTION, f"Trio {trio}: roles from the single-copy k-mers: sisters {sis[0]} and {sis[1]} (the pair sharing "
                   f"the most single-copy k-mers the third genome lacks), donor {r.D}; recipient {r.R}, the sister "
                   f"sharing more of them with {r.D} ({sign})")
    return [StatRow(LAYER, "roles", f"{r.R},{r.S};{r.D}", "all", "source", r.source, note=note)], f


def compute(lib_run: LibraryRun, trio: tuple[str, str, str], props: RouteProps) -> list[StatRow]:
    """The library statistics of `trio` (R, S, D: three distinct genomes of the sheet) on a finished screen (module
    docstring): the trio's rows, the f-ratio rows of every second donor, the evenness rows of the trio's pairs, with
    the gate's notes (`hold`). Jackknife blocks by the route's block rule (`block_bp`). An AnalysisError for another
    trio, or pair files that do not match the libraries."""
    trio = tuple(trio)
    if len(trio) != 3 or len(set(trio)) != 3 or not set(trio) <= set(lib_run.labels):
        raise AnalysisError(f"library statistics: the trio must be three distinct genomes of the sheet "
                            f"({', '.join(lib_run.labels)}), not {', '.join(map(str, trio)) or 'none'}")
    t0 = time.perf_counter()
    R, S, D = trio
    t = _Trio(lib_run, R, S, D, block_bp(lib_run, props))
    rows = t.kmers() + t.patterns() + t.divergences()
    for d2 in second_donors(lib_run, R, S, D):
        rows += t.fratio(d2)
    rows += evenness(Path(lib_run.outdir) / "pairs", {g: t.L[g].lib.elements for g in t.order})
    log.info("library statistics of %s: %d rows in %.1f s", t.scope, len(rows), time.perf_counter() - t0)
    return hold(rows)


# ---------------------------------------------------------------- the every-pair evenness check


def _old_states(path: Path, el: ElementTable) -> tuple[np.ndarray, np.ndarray, set[str], bool]:
    """For the genome's old insertions (K2P above the OLD_QUANTILE quantile of its table's finite K2P) resolved in
    the other genome (present or absent there, by an absence `evenness` accepts), whether each is present (0/1) and
    its EVEN_BLOCK_BP block on the genome; the absence states (of any kind) among the calls; whether any K2P is
    known."""
    _finished(path)
    k = el.column("k")
    known = bool(np.isfinite(k).any())
    k_old = np.nanquantile(k, OLD_QUANTILE) if known else np.inf
    at = {e: i for i, e in enumerate(el.column("eid").tolist()) if e is not None}
    code = classes(ABSENCES["evenness"])
    absent = classes(_ANY_ABSENT)
    st = np.zeros(len(el), np.int8)
    seen = set()
    for line, (eid, s) in stream_calls(path, ("eid", "state")):
        if eid not in at or s not in code:
            raise AnalysisError(f"{path}: line {line}: {eid} {s} is no call of the genome's elements")
        st[at[eid]] = code[s]
        if absent[s] == 1:
            seen.add(s)
    keep = (k > k_old) & (st > 0)
    seqs = np.unique(el.column("seqid")[keep], return_inverse=True)[1]
    mid = (el.column("start")[keep] + el.column("end")[keep]) // 2
    return (st[keep] == 2).astype(float), block_ids(seqs, mid, EVEN_BLOCK_BP), seen, known


def _log_share(present: np.ndarray, blocks: np.ndarray) -> Jk:
    """log of the share present, with its delete-one block jackknife (NaN without a present one)."""
    def stat(keep: np.ndarray) -> float:
        n, p = int(keep.sum()), float(present[keep].sum())
        return math.log(p / n) if n and p else math.nan
    return jackknife(stat, blocks)


def evenness(pairs_dir: Path, elements: dict[str, ElementTable]) -> list[StatRow]:
    """The every-pair evenness check of the genomes of `elements` (label -> its element table, with each copy's K2P
    and position), from their pair call files in `pairs_dir` (<X>__<Y>.calls.tsv.gz). Per pair X, Y (in `elements`'
    order), scope "X,Y", stratum "old insertions": r = the share of X's old insertions (LTR K2P above the
    OLD_QUANTILE quantile of X's table's finite K2P) present (P) in Y among those resolved in Y (P or A), over the
    same share of Y's old insertions in X: Y's annotation completeness relative to X's (n: the resolved old
    insertions of both); log_r with se = sqrt(var_X + var_Y), each the delete-one jackknife variance of the log share
    over EVEN_BLOCK_BP blocks on its own genome; share_<X>_in_<Y> and share_<Y>_in_<X> (n each). r is NaN when either
    share is 0 or has no resolved old insertion (or no K2P is known), or when the calls hold absences none of a kind
    evenness accepts (ABSENCES). Then the check (scope the labels, stratum "all"): the lowest min(r, 1/r) over the
    pairs; it passes at EVEN_CUTOFF for every pair (`gate`).

    The quantile is taken over every row of a genome's table. Pass each library's whole table (Records.lib.elements,
    as `compute` does): its records whose names hold no coordinates have no eid and no call, and count in the
    quantile only, as in the completeness ratio. A table of placed records alone (elements.tsv) takes the quantile
    over those."""
    labels = list(elements)
    pairs_dir = Path(pairs_dir)
    old = {(x, y): _old_states(pairs_dir / f"{x}{PAIR_SEP}{y}.calls.tsv.gz", elements[x])
           for x in labels for y in labels if x != y}
    rows, checks = [], []
    for i, x in enumerate(labels):
        for y in labels[i + 1:]:
            (p1, b1, seen1, known1), (p2, b2, seen2, known2) = old[x, y], old[y, x]
            refused = _why_unaccepted("evenness", seen1 | seen2)
            s1 = float(p1.mean()) if p1.size and not refused else math.nan
            s2 = float(p2.mean()) if p2.size and not refused else math.nan
            r = s1 / s2 if s1 > 0 and s2 > 0 else math.nan
            j1, j2 = _log_share(p1, b1), _log_share(p2, b2)
            se = math.sqrt(j1.se ** 2 + j2.se ** 2) if not refused else math.nan
            pair = f"{x},{y}"
            what = (f"share of {x}'s old insertions present in {y} among those resolved there, over the same share of "
                    f"{y}'s in {x}: {y}'s annotation completeness relative to {x}'s")
            unknown = " and ".join(g for g, k in ((x, known1), (y, known2)) if not k)
            if math.isnan(r):
                what += (f"; NaN: {refused}" if refused else f"; NaN: no K2P known for {unknown}" if unknown else
                         "; NaN: no old insertion resolved" if not (p1.size and p2.size) else
                         "; NaN: no old insertion present")
            rows += [StatRow(LAYER, "evenness", pair, "old insertions", "r", r, n=int(p1.size + p2.size), note=what),
                     StatRow(LAYER, "evenness", pair, "old insertions", "log_r", math.log(r) if r > 0 else math.nan,
                             se, note=f"se: delete-one jackknife over {EVEN_BLOCK_BP:,}-bp blocks on each genome, the "
                                      "two variances added (printed only)"),
                     StatRow(LAYER, "evenness", pair, "old insertions", f"share_{x}_in_{y}", s1, n=int(p1.size),
                             note=f"{x}'s old insertions (LTR K2P above its {OLD_QUANTILE:g} quantile) present in {y} "
                                  "among those resolved there"),
                     StatRow(LAYER, "evenness", pair, "old insertions", f"share_{y}_in_{x}", s2, n=int(p2.size),
                             note=f"{y}'s old insertions present in {x} among those resolved there")]
            checks.append((min(r, 1 / r) if r > 0 else math.nan, pair))
    if checks:
        value, pair = min(checks, key=lambda c: (not math.isnan(c[0]), c[0]))
        measured = "; ".join(f"{name} {lo:.3f}-{hi:.3f} ({held}/{n} held)" for name, n, (lo, hi), held in EVEN_MEASURED)
        rows.append(StatRow(LAYER, "evenness", ",".join(labels), "all", "check", value,
                            note=f"the lowest min(r, 1/r), pair {pair}; the check passes when every pair reaches "
                                 f"{EVEN_CUTOFF} (simulated runs: {measured})"))
    return rows


# ---------------------------------------------------------------- the gate


def _value(r: StatRow) -> bool:
    """The row's value is a positive finite number."""
    return isinstance(r.value, (int, float)) and math.isfinite(r.value) and r.value > 0


def _even(r: StatRow) -> bool:
    return _value(r) and min(r.value, 1 / r.value) >= EVEN_CUTOFF


def gate(rows: Iterable[StatRow]) -> dict[str, str]:
    """The held statistics' notes by test: D_ins and D_kmer (with their strata) always "held: annotation-sensitive";
    D3 and D_self "held: uneven annotation" unless every evenness pair's r, or its inverse, is at least EVEN_CUTOFF
    (a pair without a value, or no evenness row, fails the check). `rows` are one trio's (`compute` gates each trio
    on its own): every evenness row given counts, so the rows of several trios at once would share one check."""
    ratios = [r for r in rows if r.test == "evenness" and r.statistic == "r"]
    notes = dict.fromkeys(SENSITIVE, HELD_SENSITIVE)
    if not ratios or not all(_even(r) for r in ratios):
        notes.update(dict.fromkeys(EVEN_GATED, HELD_UNEVEN))
    return notes


def hold(rows: list[StatRow]) -> list[StatRow]:
    """The rows with the gate's notes (`gate`): every row of a held statistic's note begins with its flag; an earlier
    flag is replaced, and removed from a statistic no longer held. `rows` are one trio's, as for `gate`."""
    notes = gate(rows)
    out = []
    for r in rows:
        base = _HELD.sub("", r.note) or "."
        note = _join(notes[r.test], base) if r.test in notes else base
        out.append(r if note == r.note else replace(r, note=note))
    return out


# ---------------------------------------------------------------- findings


def findings(rows: list[StatRow]) -> Findings:
    """Statements for findings.txt (section "Trio tests") from one trio's rows (`compute`): each statistic with its
    value, the four with simulated rates beside them; the every-pair evenness check; then the verdicts: D_ins and
    D_kmer never give one; D3 and D_self give one when they are not held, else say why and point to `synltr
    geneflow` without --screen."""
    f = Findings()
    get = {(r.test, r.stratum, r.statistic): r for r in rows}
    head = next((r for r in rows if r.test == "D_kmer" and r.statistic == "D"), None)
    if head is None:
        return f
    scope = head.scope
    (R, S), D = scope.split(";")[0].split(","), scope.split(";")[1]

    def v(test, stratum, stat, spec=".4g"):
        r = get.get((test, stratum, stat))
        return "nan" if r is None else (fmt(float(r.value), spec) if isinstance(r.value, float) else str(r.value))

    f.add(SECTION, f"D(single-copy k-mers; {R},{S};{D}) = {v('D_kmer', 'all', 'D', '+.4f')} +- "
                   f"{fmt(head.se, '.4f')} (Z = {fmt(head.z, '.1f')}, {v('D_kmer', 'all', 'nblocks')} jackknife blocks "
                   f"on {D}); {rates('D_kmer')}.")
    stratum = f"{D} elements inside synteny"
    cr, cs, cb = (float(get[("donor_resolved", stratum, s)].value) for s in ("resolved_R", "resolved_S",
                                                                            "resolved_both"))
    line = (f"{D} elements inside synteny with {R}: {fmt(cr, '.2f')}, with {S}: {fmt(cs, '.2f')}, with both: "
            f"{fmt(cb, '.2f')}.")
    if cb < RESOLVED_BOTH_WARN or abs(cr - cs) > RESOLVED_GAP_WARN:
        f.warn(SECTION, line + " Little or uneven shared synteny (a distant donor or sparse annotation): insertion and "
                               "ortholog tests cover only the jointly resolved part and can be unreliable (heuristic "
                               f"levels {RESOLVED_BOTH_WARN:g} and {RESOLVED_GAP_WARN:g}).")
    else:
        f.add(SECTION, line)
    d_ins = get.get(("D_ins", "all", "D"))
    if d_ins is not None:
        f.add(SECTION, f"D(insertions; {R},{S};{D}) = {v('D_ins', 'all', 'D', '+.4f')} +- {fmt(d_ins.se, '.4f')} (Z = "
                       f"{fmt(d_ins.z, '.1f')}; {v('D_ins', 'all', 'nA')} {R}+{D} vs {v('D_ins', 'all', 'nB')} {S}+{D} "
                       f"insertions absent from the third); {rates('D_ins')}.")
    ratio = next((r for r in rows if r.test == "completeness_ratio" and r.statistic == "ratio"), None)
    if ratio is not None:
        k = ratio.stratum.rsplit("K2P>", 1)[1]
        differs = math.isfinite(ratio.se) and abs(ratio.value - 1) > COMPLETENESS_Z * ratio.se
        f.add(SECTION, f"Annotation completeness: {D}'s old insertions (K2P > {k}) present in {R} vs {S}: ratio "
                       f"{fmt(float(ratio.value), '.3f')} +- {fmt(ratio.se, '.3f')}"
                       + (f", beyond {COMPLETENESS_Z:g} se from 1" if differs else "")
                       + ". A value only: one-sided annotation dropout (or old polymorphisms carried by gene flow) "
                         "biases D_kmer and D_ins toward the more complete sister; D3 and D_self held the null in "
                         "simulation up to a 10-point completeness gap, not beyond, so the every-pair evenness check "
                         "below decides whether they give a verdict.")
    fl = "f_lower_bound", "identity_spectrum"
    if (*fl, "f") in get:
        lo_t, lo_m = get[(*fl, "n_target_recent")].value, float(get[(*fl, "n_mirror_recent")].value)
        f.add(SECTION, f"Identity spectrum: {R}-{D} vs mirror {S}-{D} orthologs (equal totals in the deep half): "
                       f"ratio at d >= {v(*fl, 'median_d', '.3g')} = 1 by calibration; at d < {IDENTITY_RECENT:g}: "
                       f"{lo_t} vs {fmt(lo_m, '.0f')} ({fmt(lo_t / max(lo_m, 1), '.2f')}x). Excess = "
                       f"admixture-fraction lower bound {v(*fl, 'f', '.3f')} (fraction of {R}'s {D}-orthologs).")
    d3, ds = get.get(("D3_triples", "all", "D")), get.get(("D_self_outgroup", "LTR sites", "D"))
    if d3 is not None and ds is not None:
        f.add(SECTION, f"D3 (ortholog triples, divergence; sign as D): {v('D3_triples', 'all', 'D', '+.4f')} (Z = "
                       f"{fmt(d3.z, '.1f')}; {int(get[('D3_triples', 'all', 'triples')].value):,} triples); "
                       f"{rates('D3_triples')}. "
                       f"D(self-outgroup: sister LTR polarizes) = {v('D_self_outgroup', 'LTR sites', 'D', '+.4f')} "
                       f"(Z = {fmt(ds.z, '.1f')}; ABBA {v('D_self_outgroup', 'LTR sites', 'nA')}, BABA "
                       f"{v('D_self_outgroup', 'LTR sites', 'nB')}, BBAA {v('D_self_outgroup', 'LTR sites', 'BBAA')}); "
                       f"{rates('D_self_outgroup')}.")
    for r in rows:
        if r.test == "f_ratio_ins" and r.statistic == "alpha":
            d2 = r.stratum.removeprefix("second donor ")
            ident = get[("f_ratio_identity", r.stratum, "alpha")].value
            why = (r.note.split(": NaN, ", 1)[1] if ": NaN, " in r.note else
                   "se NaN, " + r.note.split("; se NaN, ", 1)[1] if "; se NaN, " in r.note else None)
            f.add(SECTION, f"Admixture fraction of {R} from {D}'s population, with {d2} as a second donor genome: "
                           f"{fmt(float(r.value), '.3f')} +- {fmt(r.se, '.3f')} (insertion patterns"
                           + (f"; {why}" if why else "")
                           + f"), {fmt(float(ident), '.3f')} (identity-spectrum excess ratio).")
    pairs = [r for r in rows if r.test == "evenness" and r.statistic == "r"]
    check = next((r for r in rows if r.test == "evenness" and r.statistic == "check"), None)
    even = not any(gate(rows).get(t) for t in EVEN_GATED)
    if pairs and check is not None:
        se = {r.scope: r.se for r in rows if r.test == "evenness" and r.statistic == "log_r"}
        f.add(SECTION, "Every-pair evenness check (share of each genome's old insertions present in the other, X's "
                       "over Y's): " + ", ".join(f"{r.scope} r = {fmt(float(r.value), '.3f')} (se of log r "
                                                  f"{fmt(se.get(r.scope, math.nan), '.3f')})" for r in pairs)
                       + f"; lowest min(r, 1/r) {fmt(float(check.value), '.3f')}: "
                       + ("passes" if even else "fails") + f" (at least {EVEN_CUTOFF} for every pair).")
    f.add(SECTION, f"D_ins and D_kmer give no verdict: annotation-sensitive (false signals in "
                   f"{'/'.join(map(str, FALSE_SIGNALS['D_ins']['one sister 10 points less complete']))} and "
                   f"{'/'.join(map(str, FALSE_SIGNALS['D_kmer']['one sister 10 points less complete']))} simulated "
                   "nulls with one sister 10 points less completely annotated, a gap no check separates reliably "
                   "from even annotation); for a verdict, run `synltr geneflow` without --screen on the genome "
                   "assemblies.")
    if even:
        for name, r in (("D3", d3), ("D_self", ds)):
            if r is None:
                continue
            if math.isnan(r.z):
                f.add(SECTION, f"Verdict, {name}: no verdict (z not computed).")
                continue
            verdict = (f"a signal of gene flow between {R if r.z > 0 else S} and {D}" if abs(r.z) >= Z_SIG
                       else "no signal of gene flow")
            f.add(SECTION, f"Verdict, {name}: {verdict} (Z = {fmt(r.z, '.1f')}; |Z| >= {Z_SIG:g} counts).")
    else:
        low = min(pairs, key=lambda r: min(r.value, 1 / r.value) if _value(r) else -1.0, default=None)
        why = (f"the every-pair evenness check failed (lowest pair {low.scope}, min(r, 1/r) "
               f"{fmt(float(check.value), '.3f')}, below {EVEN_CUTOFF})" if low is not None and _value(low)
               else "the every-pair evenness check has no value")
        gap = "one sister 30 points less complete"
        f.add(SECTION, f"D3 and D_self give no verdict: {why}, and uneven annotation can fake them (false signals in "
                       f"{'/'.join(map(str, FALSE_SIGNALS['D3_triples'][gap]))} and "
                       f"{'/'.join(map(str, FALSE_SIGNALS['D_self_outgroup'][gap]))} simulated nulls with one sister "
                       "30 points less complete); run `synltr geneflow` without --screen on the genome assemblies, "
                       "whose D_pav gave no false signal at either gap ("
                       + ", ".join(f"{n}/{of}" for n, of in FALSE_SIGNALS["D_pav"].values()) + ").")
    return f
