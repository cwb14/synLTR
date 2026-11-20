"""Alignment blocks: PAF -> compact numpy structured array (0-based half-open), split at long indels."""

from __future__ import annotations

import logging
import re
from pathlib import Path

import numpy as np

from ltrtrace.base.constants import PAF_BATCH_LINES, RECORD_WARNINGS, SPLIT_INDEL
from ltrtrace.io.files import open_input
from ltrtrace.seq.cigar import ALIGNMENT_OPS, OPS, parse, walk

log = logging.getLogger(__name__)

_INDEL = re.compile(r"(\d+)[IDN]")  # a PAF CIGAR's indel lengths, read before the CIGAR is parsed
_INDELS = (1, 2, 3)  # I, D, N

BLOCK_DTYPE = np.dtype([
    ("a_chrom", "i4"), ("a_start", "i8"), ("a_end", "i8"),
    ("b_chrom", "i4"), ("b_start", "i8"), ("b_end", "i8"),
    ("strand", "i1"), ("matches", "i4"), ("alen", "i4"),
    ("chain", "i4"), ("chain_score", "f8"), ("margin", "f4"),
])  # chain=-1, margin=nan until set

_ROW = [("a_chrom", "i4"), ("a_start", "i8"), ("a_end", "i8"), ("b_chrom", "i4"),
        ("b_start", "i8"), ("b_end", "i8"), ("strand", "i1"), ("matches", "i4"), ("alen", "i4")]


def blocks_from_rows(rows: list[tuple]) -> np.ndarray:
    """Blocks from (a_chrom, a_start, a_end, b_chrom, b_start, b_end, strand, matches, alen) tuples."""
    tmp = np.array(rows, dtype=_ROW)
    out = np.zeros(len(tmp), dtype=BLOCK_DTYPE)
    for name, _ in _ROW:
        out[name] = tmp[name]
    out["chain"] = -1
    out["margin"] = np.nan
    return out


def read_paf(path: Path, a_index: dict[str, int], b_index: dict[str, int],
             swap: bool = False, split_indel: int = SPLIT_INDEL) -> np.ndarray:
    """Read a PAF (plain or gzip, told by its first bytes) into blocks. A = query, B = target (swapped if `swap`).

    A record carrying a cg:Z: CIGAR is split at its I/D ops of >= split_indel bp into sub-blocks with exact
    coordinates on both genomes (`split_cigar`, strand-aware), so a long indel lies between blocks, never
    inside one. A sub-block's matches are its '=' columns when the CIGAR uses =/X, else the record's matches
    prorated by sub-block alen (its aligned and shorter-indel columns). A record without a CIGAR, or without
    an indel that long, stays one block with the record's matches and alignment length; so does a record whose
    CIGAR is malformed or disagrees with its coordinates, counted in one warning. Lines starting with '#' or
    blank are skipped; a record with < 12 fields, or with coordinates or counts that are not integers or spans
    that are empty, is malformed: skipped, with a warning naming the file and line for the first
    RECORD_WARNINGS and a count in one warning. Records with a seqid missing from its index are skipped and
    counted in one warning.
    """
    path = Path(path)
    qi, ti = (b_index, a_index) if swap else (a_index, b_index)
    parts: list[np.ndarray] = []
    rows: list[tuple] = []
    skipped = bad_cigar = malformed = 0
    with open_input(path) as fh:
        for n, line in enumerate(fh, 1):
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\r\n").split("\t")
            try:
                if len(f) < 12:
                    raise ValueError(f"{len(f)} fields, PAF has 12")
                qs, qe, ts, te, nm, al = (int(f[i]) for i in (2, 3, 7, 8, 9, 10))
                if not (0 <= qs < qe and 0 <= ts < te):
                    raise ValueError(f"empty or negative span: query {qs}-{qe}, target {ts}-{te}")
            except ValueError as exc:
                malformed += 1
                if malformed <= RECORD_WARNINGS:
                    log.warning("%s:%d: PAF record skipped, malformed: %s", path, n, exc)
                continue
            q, t = qi.get(f[0]), ti.get(f[5])
            if q is None or t is None:
                skipped += 1
                continue
            strand = 1 if f[4] == "+" else -1
            pieces = [(qs, qe, ts, te, nm, al)]
            cg = next((x for x in f[12:] if x.startswith("cg:Z:")), None)
            if cg is not None and any(int(n) >= split_indel for n in _INDEL.findall(cg)):
                try:
                    pieces = _paf_pieces(cg[5:], qs, qe, ts, te, strand, nm, split_indel)
                except ValueError:
                    bad_cigar += 1
            for q0, q1, r0, r1, m, n in pieces:
                if swap:
                    rows.append((t, r0, r1, q, q0, q1, strand, m, n))
                else:
                    rows.append((q, q0, q1, t, r0, r1, strand, m, n))
            if len(rows) >= PAF_BATCH_LINES:
                parts.append(blocks_from_rows(rows))
                rows = []
    if rows:
        parts.append(blocks_from_rows(rows))
    if malformed:
        log.warning("skipped %d malformed PAF record(s) (%s)", malformed, path)
    if skipped:
        log.warning("skipped %d PAF records with seqids not in the index (%s)", skipped, path)
    if bad_cigar:
        log.warning("kept %d PAF record(s) whole: CIGAR malformed or disagreeing with the coordinates (%s)",
                    bad_cigar, path)
    return np.concatenate(parts) if parts else np.zeros(0, dtype=BLOCK_DTYPE)


def _paf_pieces(cg: str, qs: int, qe: int, ts: int, te: int, strand: int, matches: int,
                split_indel: int) -> list[tuple[int, ...]]:
    """Sub-blocks (q0, q1, r0, r1, matches, alen) of one PAF record from its CIGAR string `cg`.

    ValueError if the CIGAR is malformed, disagrees with the record's coordinates or aligns nothing. With M
    ops the record's matches are prorated by sub-block alen (cumulative rounding keeps their total).
    """
    ops, lens = parse(cg, "sam")
    q, r = walk(ops, lens)
    if q[-1] != qe - qs or r[-1] != te - ts:
        raise ValueError("CIGAR disagrees with the record's coordinates")
    pieces, has_m = split_cigar(ops, lens, qs, ts, te, strand, split_indel)
    if not pieces:
        raise ValueError("CIGAR aligns no column")
    if not has_m:
        return pieces
    total = sum(p[5] for p in pieces)
    out, cum, prev = [], 0, 0
    for q0, q1, r0, r1, _, alen in pieces:
        cum += alen
        m = round(matches * cum / total)
        out.append((q0, q1, r0, r1, m - prev, alen))
        prev = m
    return out


def split_cigar(ops, lens, q_st: int, r_st: int, r_en: int, strand: int,
                split_indel: int) -> tuple[list[tuple[int, ...]], bool]:
    """One alignment cut at each indel (I, D or N) of >= split_indel bp: ([(q0, q1, r0, r1, eq, alen)], has_m).

    `ops` and `lens` hold the CIGAR (seq.cigar) in reference order, as in SAM, PAF cg:Z: and mappy: on the
    minus strand the query runs backward, so the walk takes the ops in reverse, from r_en down, and the pieces
    come out in query order. Coordinates are forward-strand, 0-based half-open. A piece spans the ops between two
    long indels (or an alignment end), shorter indels at its edges included, and is kept when it holds query and
    reference bases. eq counts '=' columns, alen the M, =, X and shorter-indel columns; has_m tells whether M ops
    hide the match count. S, H and P ops raise ValueError.
    """
    plus = strand > 0
    ops = np.asarray(ops, dtype=np.int64)
    lens = np.asarray(lens, dtype=np.int64)
    if not plus:
        ops, lens = ops[::-1], lens[::-1]
    bad = np.flatnonzero(~np.isin(ops, ALIGNMENT_OPS))
    if len(bad):
        op = int(ops[bad[0]])
        raise ValueError(f"CIGAR op {OPS[op] if 0 <= op < len(OPS) else op} not allowed")
    q, d = walk(ops, lens)  # d: reference bases walked, along the walk (descending positions on '-')
    eq = np.concatenate(([0], np.cumsum(np.where(ops == 7, lens, 0))))
    cols = np.concatenate(([0], np.cumsum(lens)))
    cuts = np.flatnonzero(np.isin(ops, _INDELS) & (lens >= split_indel)).tolist()
    out = []
    for i0, i1 in zip([0, *(c + 1 for c in cuts)], [*cuts, len(ops)]):
        q0, q1, d0, d1 = q_st + int(q[i0]), q_st + int(q[i1]), int(d[i0]), int(d[i1])
        if q1 > q0 and d1 != d0:
            r0, r1 = (r_st + d0, r_st + d1) if plus else (r_en - d1, r_en - d0)
            out.append((q0, q1, r0, r1, int(eq[i1] - eq[i0]), int(cols[i1] - cols[i0])))
    return out, bool((ops == 0).any())


def a_to_b(block: np.void, apos: int) -> int:
    """Map an A position to B by linear scaling, clamped to [b_start, b_end - 1]; a block with an empty A span
    maps every position to its first B base on its strand."""
    a0, a1 = int(block["a_start"]), int(block["a_end"])
    b0, b1 = int(block["b_start"]), int(block["b_end"])
    off = round((apos - a0) * (b1 - b0) / (a1 - a0)) if a1 > a0 else 0
    pos = b0 + off if block["strand"] > 0 else b1 - 1 - off
    return max(b0, min(b1 - 1, pos))


def save_blocks(path: Path, blocks: np.ndarray, a_names: list[str], b_names: list[str]) -> None:
    """Write blocks and both name lists to one compressed .npz."""
    with open(path, "wb") as fh:
        np.savez_compressed(fh, blocks=blocks, a_names=np.array(a_names, dtype=np.str_),
                            b_names=np.array(b_names, dtype=np.str_))


def load_blocks(path: Path) -> tuple[np.ndarray, list[str], list[str]]:
    with np.load(path, allow_pickle=False) as z:
        return (z["blocks"].astype(BLOCK_DTYPE, copy=False),
                [str(x) for x in z["a_names"]], [str(x) for x in z["b_names"]])
