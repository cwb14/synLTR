#!/usr/bin/env python3
"""solo_ltr_core.py - shared logic for the solo-LTR detection pipeline.

Both the single-run CLI (solo_ltr_filter.py) and the parameter sweep
(grid_search.py) import these functions, so the filtering logic lives in
exactly one place.

Pipeline (after the cached BLAST in solo_ltr_search.py):

  consensus HSPs --[pident/qcov/len/evalue]--> kept HSPs
                 --[TSD at HSP boundaries]----> solo-like HSPs
                 --[merge within slop]---------> candidate loci
                 --[drop if internal HSP within flank]--> solo-LTR calls

The TSD test is the key discriminator: a true solo-LTR retains the target
site duplication (a short direct repeat flanking the element), whereas a
truncated LTR fragment does not.  Because the BLAST alignment boundaries
(sstart/send) sit at or very near the LTR termini, we look for an identical
k-mer just-outside the 5' boundary and just-outside the 3' boundary, allowing
a few bp of slop for boundary wobble.

Coordinates are 0-based half-open (BED) throughout.
"""
import bisect
import collections
import re
import sys

# Tabular BLAST columns produced by solo_ltr_search.py
BLAST_COLS = (
    "qseqid sseqid pident length qlen slen "
    "qstart qend sstart send evalue bitscore qcovs qcovhsp mismatch gapopen"
).split()
COL = {c: i for i, c in enumerate(BLAST_COLS)}


_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq):
    """Reverse complement; N maps to N. Case-preserving."""
    return seq.translate(_COMP)[::-1]


# Header: '>FULLid \t chrom:start-end' where start-end is 1-based inclusive.
_INTERNAL_COORD_RE = re.compile(r"^(.+):(\d+)-(\d+)$")


def build_orient_map(internal_fa, genome, only=None, min_ident=0.6):
    """Map internal-FASTA record id -> (chrom, g0, g1, strand).

    The id is the first whitespace-delimited token of the header (matches the
    blastn `sseqid` after makeblastdb -parse_seqids). The second token is the
    internal-region coordinate chrom:start-end (1-based inclusive); the genomic
    window is g0 = start-1, g1 = end (0-based half-open).

    Orientation is detected by comparing the stored sequence to genome[g0:g1]
    over NON-N positions: '+' if the forward orientation matches at least as
    many non-N positions as the reverse complement, else '-'. A record is
    SKIPPED (warned, omitted from the map) when its header lacks an internal
    coordinate, its chrom is absent, the window length disagrees with the stored
    length, it has no non-N positions, or the best orientation matches
    < min_ident of non-N positions.

    `only`: optional set/collection of ids to build the map for (others skipped
    without genome comparison) -- lets a caller restrict work to ids that
    actually received a BLAST hit.
    """
    only = set(only) if only is not None else None
    out = {}
    cur_id = None
    cur_coord = None
    buf = []

    def _flush():
        if cur_id is None:
            return
        if only is not None and cur_id not in only:
            return
        m = _INTERNAL_COORD_RE.match(cur_coord) if cur_coord else None
        if not m:
            sys.stderr.write(f"[solo_ltr_core] orient: no internal coordinate in header for {cur_id}; skip\n")
            return
        chrom, start, end = m.group(1), int(m.group(2)), int(m.group(3))
        g = genome.get(chrom)
        if g is None:
            sys.stderr.write(f"[solo_ltr_core] orient: chrom {chrom} not in genome; skip {cur_id}\n")
            return
        g0, g1 = start - 1, end
        win = g[g0:g1].upper()
        seq = "".join(buf).upper()
        if len(win) != len(seq):
            sys.stderr.write(f"[solo_ltr_core] orient: length mismatch for {cur_id}; skip\n")
            return
        n_nonN = sum(1 for c in seq if c != "N")
        if n_nonN == 0:
            return
        fwd = sum(1 for sc, gc in zip(seq, win) if sc != "N" and sc == gc)
        rev = sum(1 for sc, gc in zip(reverse_complement(seq), win) if sc != "N" and sc == gc)
        best = max(fwd, rev)
        if best < min_ident * n_nonN:
            sys.stderr.write(f"[solo_ltr_core] orient: unresolved orientation for {cur_id}; skip\n")
            return
        out[cur_id] = (chrom, g0, g1, "+" if fwd >= rev else "-")

    with open(internal_fa) as fh:
        for line in fh:
            if line.startswith(">"):
                _flush()
                parts = line[1:].rstrip("\n").split()
                cur_id = parts[0] if parts else None
                cur_coord = parts[1] if len(parts) > 1 else None
                buf = []
            else:
                buf.append(line.strip())
    _flush()
    return out


# --------------------------------------------------------------------------- #
# Genome (for TSD lookups)
# --------------------------------------------------------------------------- #
def load_genome(path):
    """Load a FASTA into a dict chrom -> uppercase sequence string.

    Whole-genome in RAM is intentional: TSD checks do millions of tiny random
    slices, and per-slice pyfaidx/samtools calls would dominate runtime.  A
    fungal genome (~200 Mb) is a few hundred MB as a Python str -- fine.  For
    very large genomes, swap in pyfaidx and accept the slowdown.
    """
    genome = {}
    name = None
    chunks = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    genome[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        genome[name] = "".join(chunks).upper()
    return genome


# --------------------------------------------------------------------------- #
# BLAST loading
# --------------------------------------------------------------------------- #
def load_blast(path):
    """Parse tabular BLAST -> dict chrom -> list of HSP tuples.

    Each tuple: (s0, e1, pident, qcov_hsp, length, evalue) where (s0, e1) is a
    0-based half-open SUBJECT interval.  Lists are NOT sorted here.
    """
    by = collections.defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < len(BLAST_COLS):
                continue
            chrom = f[COL["sseqid"]]
            s = int(f[COL["sstart"]]); e = int(f[COL["send"]])
            if s > e:
                s, e = e, s
            by[chrom].append((s - 1, e,
                              float(f[COL["pident"]]),
                              float(f[COL["qcovhsp"]]),
                              int(f[COL["length"]]),
                              float(f[COL["evalue"]])))
    return by


def load_internal_blast(path, orient_map):
    """Parse a consensus-vs-internal blastn tsv and map SUBJECT coordinates to
    GENOMIC coordinates, returning the same shape as load_blast():
        dict chrom -> list of (s0, e1, pident, qcovhsp, length, evalue)
    with (s0, e1) a 0-based half-open genomic interval, so the result feeds
    make_candidates() directly.

    orient_map: id -> (chrom, g0, g1, strand) from build_orient_map(). A hit
    whose sseqid is not in orient_map is skipped.

    Subject mapping (validated): with lo=min(sstart,send), hi=max(sstart,send):
        '+' : [ g0 + (lo - 1),  g0 + hi )
        '-' : [ g1 - hi,        g1 - (lo - 1) )
    """
    by = collections.defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < len(BLAST_COLS):
                continue
            rec = orient_map.get(f[COL["sseqid"]])
            if rec is None:
                continue
            chrom, g0, g1, strand = rec
            ss = int(f[COL["sstart"]]); se = int(f[COL["send"]])
            lo, hi = (ss, se) if ss <= se else (se, ss)
            if strand == "+":
                s0, e1 = g0 + (lo - 1), g0 + hi
            else:
                s0, e1 = g1 - hi, g1 - (lo - 1)
            by[chrom].append((s0, e1,
                              float(f[COL["pident"]]),
                              float(f[COL["qcovhsp"]]),
                              int(f[COL["length"]]),
                              float(f[COL["evalue"]])))
    return by


# --------------------------------------------------------------------------- #
# TSD detection
# --------------------------------------------------------------------------- #
def has_tsd(seq, s, e, k, slop):
    """True if an identical k-mer flanks the [s, e) interval within +/-slop.

    Looks for seq[s5-k : s5] == seq[e3 : e3+k] over boundary offsets
    o5, o3 in [-slop, slop].  Ns disqualify a comparison.
    """
    Lc = len(seq)
    for o5 in range(-slop, slop + 1):
        s5 = s + o5
        if s5 - k < 0 or s5 > Lc:
            continue
        a = seq[s5 - k:s5]
        if "N" in a or len(a) < k:
            continue
        for o3 in range(-slop, slop + 1):
            e3 = e + o3
            if e3 < 0 or e3 + k > Lc:
                continue
            b = seq[e3:e3 + k]
            if "N" in b or len(b) < k:
                continue
            if a == b:
                return True
    return False


# --------------------------------------------------------------------------- #
# Candidate construction
# --------------------------------------------------------------------------- #
def make_candidates(cons_by_chrom, genome, *, min_pident, min_qcov,
                    min_length, max_evalue, use_tsd, tsd_k, tsd_slop,
                    merge_slop):
    """Filter consensus HSPs, optionally TSD-gate, then merge into loci.

    Returns dict chrom -> sorted list of (start, end) candidate intervals.
    """
    cands = {}
    for chrom, arr in cons_by_chrom.items():
        seq = genome.get(chrom)
        kept = []
        for s, e, p, q, L, ev in arr:
            if p < min_pident or q < min_qcov or L < min_length or ev > max_evalue:
                continue
            if use_tsd:
                if seq is None or not has_tsd(seq, s, e, tsd_k, tsd_slop):
                    continue
            kept.append((s, e))
        if not kept:
            continue
        kept.sort()
        merged = [list(kept[0])]
        for s, e in kept[1:]:
            if s <= merged[-1][1] + merge_slop:
                if e > merged[-1][1]:
                    merged[-1][1] = e
            else:
                merged.append([s, e])
        cands[chrom] = [(s, e) for s, e in merged]
    return cands


def prefilter_internal(internal_by_chrom, *, min_pident, min_qcov,
                       min_length, max_evalue):
    """Filter + sort internal HSPs once -> dict chrom -> sorted [(s, e), ...].

    Split out from filter_internal so a parameter sweep can cache this (the
    internal cutoffs rarely change) instead of re-scanning millions of HSPs
    for every flank value.
    """
    internal = {}
    for chrom, arr in internal_by_chrom.items():
        v = [(s, e) for s, e, p, q, L, ev in arr
             if p >= min_pident and q >= min_qcov
             and L >= min_length and ev <= max_evalue]
        v.sort()
        internal[chrom] = v
    return internal


def drop_near_internal(cands, internal, flank):
    """Drop candidates with any pre-filtered internal interval within flank bp.

    `internal` is the output of prefilter_internal().
    """
    out = {}
    for chrom, arr in cands.items():
        iarr = internal.get(chrom, [])
        if not iarr:
            out[chrom] = list(arr)
            continue
        istart = [x[0] for x in iarr]
        keep = []
        for s, e in arr:
            j = bisect.bisect_left(istart, s - flank - 50_000)
            near = False
            while j < len(iarr):
                isx, iex = iarr[j]
                if isx > e + flank:
                    break
                if iex <= s:
                    gap = s - iex
                elif isx >= e:
                    gap = isx - e
                else:
                    gap = 0
                if gap <= flank:
                    near = True
                    break
                j += 1
            if not near:
                keep.append((s, e))
        out[chrom] = keep
    return out


def filter_internal(cands, internal_by_chrom, *, min_pident, min_qcov,
                    min_length, max_evalue, flank):
    """Convenience: prefilter_internal + drop_near_internal in one call.

    Used by the single-run CLIs.  A sweep should call the two pieces
    separately and cache prefilter_internal().
    """
    internal = prefilter_internal(
        internal_by_chrom, min_pident=min_pident, min_qcov=min_qcov,
        min_length=min_length, max_evalue=max_evalue)
    return drop_near_internal(cands, internal, flank)


def load_internal_intervals(internal_fa):
    """Parse internal-FASTA headers -> {chrom: sorted [(g0, g1), ...]} of the
    INTERNAL-region genomic intervals (the 2nd, tab-separated header coordinate,
    1-based inclusive -> 0-based half-open). Used by the nested-intact guard.
    """
    by = collections.defaultdict(list)
    with open(internal_fa) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            parts = line[1:].rstrip("\n").split()
            if len(parts) < 2:
                continue
            m = _INTERNAL_COORD_RE.match(parts[1])
            if not m:
                continue
            chrom, start, end = m.group(1), int(m.group(2)), int(m.group(3))
            by[chrom].append((start - 1, end))
    for c in by:
        by[c].sort()
    return by


def drop_ltr_of_nested_intact(cands, internal_intervals, flank):
    """Remove nested candidates that look like one LTR of a nested INTACT element.

    A true nested solo lies strictly inside a single host internal region (both
    flanks are the same element). An LTR of an intact element X nested in host H
    abuts the boundary of X's OWN internal region (a separate record). So drop a
    candidate if some internal interval is within `flank` bp of it WITHOUT cleanly
    containing it (margin >= flank on both sides); keep it if every nearby internal
    interval is a clean container.
    """
    out = {}
    for chrom, arr in cands.items():
        ivs = internal_intervals.get(chrom, [])
        keep = []
        for s, e in arr:
            drop = False
            for is_, ie in ivs:
                if ie <= s - flank:
                    continue
                if is_ >= e + flank:
                    break
                if is_ <= s and e <= ie and (s - is_) >= flank and (ie - e) >= flank:
                    continue  # clean host container
                if max(is_ - e, s - ie, 0) <= flank:
                    drop = True
                    break
            if not drop:
                keep.append((s, e))
        out[chrom] = keep
    return out


def write_bed(cands, path):
    n = 0
    with open(path, "w") as fh:
        for chrom in sorted(cands):
            for s, e in cands[chrom]:
                fh.write(f"{chrom}\t{s}\t{e}\n")
                n += 1
    return n


def union_candidates(*cand_dicts, merge_slop=0):
    """Combine candidate dicts {chrom: [(s, e), ...]} into one, merging
    intervals that overlap or sit within merge_slop bp. Returns a sorted,
    merged dict. Inputs are not mutated.
    """
    chroms = set()
    for d in cand_dicts:
        chroms.update(d)
    out = {}
    for c in sorted(chroms):
        ivs = []
        for d in cand_dicts:
            ivs.extend(d.get(c, []))
        if not ivs:
            continue
        ivs.sort()
        merged = [list(ivs[0])]
        for s, e in ivs[1:]:
            if s <= merged[-1][1] + merge_slop:
                if e > merged[-1][1]:
                    merged[-1][1] = e
            else:
                merged.append([s, e])
        out[c] = [(s, e) for s, e in merged]
    return out


# --------------------------------------------------------------------------- #
# Scoring against truth (mirrors benchmark.py headline metric)
# --------------------------------------------------------------------------- #
def score(cands, solo_idx, min_overlap_frac=0.5):
    """Return dict with TP/FP/FN/precision/recall/F1 vs clean-solo truth.

    solo_idx: dict chrom -> sorted list of (start, end) clean-solo intervals.
    A prediction is a TP if it reciprocally overlaps a solo by >= frac.
    A solo is detected if any prediction reciprocally overlaps it by >= frac.
    """
    n_solo = sum(len(v) for v in solo_idx.values())
    frac = min_overlap_frac

    detected = 0
    for chrom, recs in solo_idx.items():
        cl = sorted(cands.get(chrom, []))
        for ss, se in recs:
            tlen = max(1, se - ss)
            for ps, pe in cl:
                if pe <= ss:
                    continue
                if ps >= se:
                    break
                ovl = min(pe, se) - max(ps, ss)
                if ovl > 0 and min(ovl / tlen, ovl / max(1, pe - ps)) >= frac:
                    detected += 1
                    break

    n_pred = sum(len(v) for v in cands.values())
    tp = 0
    for chrom, arr in cands.items():
        recs = solo_idx.get(chrom, [])
        for s, e in arr:
            plen = max(1, e - s)
            for ts, te in recs:
                if te <= s:
                    continue
                if ts >= e:
                    break
                ovl = min(te, e) - max(ts, s)
                if ovl > 0 and min(ovl / plen, ovl / max(1, te - ts)) >= frac:
                    tp += 1
                    break
    fp = n_pred - tp
    fn = n_solo - detected
    precision = tp / n_pred if n_pred else 0.0
    recall = detected / n_solo if n_solo else 0.0
    f1 = (2 * precision * recall / (precision + recall)
          if (precision + recall) > 0 else 0.0)
    return dict(tp=tp, fp=fp, fn=fn, detected=detected, n_solo=n_solo,
                n_pred=n_pred, precision=precision, recall=recall, f1=f1)


def load_truth_solo(bed_path):
    """clean-solo intervals from a PrinTE BED, via benchmark.stream_truth."""
    from benchmark import stream_truth
    idx = collections.defaultdict(list)
    for rec, cat in stream_truth(bed_path):
        if cat == "clean_solo":
            idx[rec.chrom].append((rec.start, rec.end))
    for c in idx:
        idx[c].sort()
    return idx
