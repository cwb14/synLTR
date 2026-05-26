#!/usr/bin/env python3
"""solo_ltr_core.py - shared logic for the solo-LTR detection pipeline.

The end-to-end CLI (solo_ltr.py) is the sole importer of this module, so the
parsing, filtering, and divergence logic lives in exactly one place.

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
import math
import re
import sys

# Tabular BLAST columns produced by solo_ltr_search.py
BLAST_COLS = (
    "qseqid sseqid pident length qlen slen "
    "qstart qend sstart send evalue bitscore qcovs qcovhsp mismatch gapopen "
    "btop"
).split()
COL = {c: i for i, c in enumerate(BLAST_COLS)}

# Per-HSP record. (s0, e1) is a 0-based half-open SUBJECT/GENOMIC interval.
# strand is the genomic strand ('+'/'-'); host is None for the main path and
# the internal-record sseqid for the nested path.
HSP = collections.namedtuple(
    "HSP",
    "s0 e1 pident qcovhsp length evalue bitscore mismatch gapopen "
    "strand qseqid btop host")

# One merged solo-LTR locus. name is assigned at write time (not stored here).
Candidate = collections.namedtuple(
    "Candidate",
    "chrom start end strand best_consensus classification pident aln_len "
    "qcov mismatch evalue bitscore tsd tsd_len pdist k2p source nested_in")


_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq):
    """Reverse complement; N maps to N. Case-preserving."""
    return seq.translate(_COMP)[::-1]


_PURINES = frozenset("AG")
_PYRIMIDINES = frozenset("CT")
_ACGT = frozenset("ACGT")


def _is_transition(a, b):
    """True if {a,b} is a purine-purine or pyrimidine-pyrimidine substitution."""
    return (a in _PURINES and b in _PURINES) or (a in _PYRIMIDINES and b in _PYRIMIDINES)


def btop_distances(btop):
    """Parse a BLAST BTOP string -> (pdist, k2p) as formatted strings.

    BTOP grammar: a run of digits = that many consecutive matches; a 2-char
    token = one alignment column where char1 = query (consensus) base and
    char2 = subject (genomic) base; '-' marks a gap.

    We count substitutions over A/C/G/T columns only:
      transition   : both purines (A/G) or both pyrimidines (C/T)
      transversion : any other A/C/G/T substitution
    Gap columns and IUPAC-ambiguous query/subject columns are excluded from the
    effective length L = match + ts + tv. With P = ts/L, Q = tv/L:
      pdist = (ts + tv) / L
      K2P   = -0.5*ln(1 - 2P - Q) - 0.25*ln(1 - 2Q)
    Returns ('NA','NA') when L == 0 (or btop empty); K2P alone is 'NA' on
    saturation (a log argument <= 0). Values are rounded to 4 dp.
    """
    if not btop:
        return "NA", "NA"
    match = ts = tv = 0
    i, n = 0, len(btop)
    while i < n:
        c = btop[i]
        if c.isdigit():
            j = i
            while j < n and btop[j].isdigit():
                j += 1
            match += int(btop[i:j])
            i = j
            continue
        if i + 1 >= n:
            break  # malformed trailing single char; ignore (excluded from L)
        a = btop[i].upper()
        b = btop[i + 1].upper()
        i += 2
        if a == "-" or b == "-":
            continue  # gap column: excluded from L
        if a not in _ACGT or b not in _ACGT:
            continue  # IUPAC-ambiguous: cannot classify ts/tv
        if a == b:
            match += 1  # defensive; BTOP should not emit equal pairs (counted as match to preserve L)
        elif _is_transition(a, b):
            ts += 1
        else:
            tv += 1
    L = match + ts + tv
    if L == 0:
        return "NA", "NA"
    pdist = (ts + tv) / L
    if ts == 0 and tv == 0:
        return f"{pdist:.4f}", "0.0000"
    P, Q = ts / L, tv / L
    a_term, b_term = 1 - 2 * P - Q, 1 - 2 * Q
    if a_term <= 0 or b_term <= 0:
        return f"{pdist:.4f}", "NA"
    k2p = -0.5 * math.log(a_term) - 0.25 * math.log(b_term)
    return f"{pdist:.4f}", f"{k2p:.4f}"


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
    """Parse tabular BLAST -> dict chrom -> list of HSP (main path; host=None).

    (s0, e1) is a 0-based half-open SUBJECT interval. strand is '+' when the
    subject alignment runs forward (sstart < send) else '-'. Lists are NOT
    sorted here. Warns if every data line is too narrow (a stale pre-btop
    cache), which would otherwise silently yield zero candidates.
    """
    by = collections.defaultdict(list)
    n_lines = n_narrow = 0
    with open(path) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            n_lines += 1
            f = line.rstrip("\n").split("\t")
            if len(f) < len(BLAST_COLS):
                n_narrow += 1
                continue
            chrom = f[COL["sseqid"]]
            s = int(f[COL["sstart"]]); e = int(f[COL["send"]])
            strand = "+" if s < e else "-"
            if s > e:
                s, e = e, s
            by[chrom].append(HSP(
                s - 1, e,
                float(f[COL["pident"]]), float(f[COL["qcovhsp"]]),
                int(f[COL["length"]]), float(f[COL["evalue"]]),
                float(f[COL["bitscore"]]), int(f[COL["mismatch"]]),
                int(f[COL["gapopen"]]), strand,
                f[COL["qseqid"]], f[COL["btop"]], None))
    if n_lines and n_narrow == n_lines:
        sys.stderr.write(
            f"[solo_ltr_core] WARNING: {path} had {n_lines} data lines but all "
            f"were too narrow ({len(BLAST_COLS)} cols expected) -- stale/"
            f"incompatible BLAST cache? re-run with --force\n")
    return by


def load_internal_blast(path, orient_map):
    """Parse consensus-vs-internal blastn -> dict chrom -> list of HSP with
    SUBJECT coords mapped to GENOMIC coords, host = internal-record sseqid.

    orient_map: id -> (chrom, g0, g1, strand) from build_orient_map(). A hit
    whose sseqid is absent is skipped. Subject mapping (lo=min, hi=max of
    sstart/send): '+' -> [g0+(lo-1), g0+hi); '-' -> [g1-hi, g1-(lo-1)). The
    genomic strand combines the record orientation with the hit orientation.
    """
    by = collections.defaultdict(list)
    n_lines = n_narrow = 0
    with open(path) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            n_lines += 1
            f = line.rstrip("\n").split("\t")
            if len(f) < len(BLAST_COLS):
                n_narrow += 1
                continue
            sseqid = f[COL["sseqid"]]
            rec = orient_map.get(sseqid)
            if rec is None:
                continue
            chrom, g0, g1, orient = rec
            ss = int(f[COL["sstart"]]); se = int(f[COL["send"]])
            lo, hi = (ss, se) if ss <= se else (se, ss)
            if orient == "+":
                s0, e1 = g0 + (lo - 1), g0 + hi
            else:
                s0, e1 = g1 - hi, g1 - (lo - 1)
            genomic_strand = "+" if (orient == "+") == (ss <= se) else "-"
            by[chrom].append(HSP(
                s0, e1,
                float(f[COL["pident"]]), float(f[COL["qcovhsp"]]),
                int(f[COL["length"]]), float(f[COL["evalue"]]),
                float(f[COL["bitscore"]]), int(f[COL["mismatch"]]),
                int(f[COL["gapopen"]]), genomic_strand,
                f[COL["qseqid"]], f[COL["btop"]], sseqid))
    if n_lines and n_narrow == n_lines:
        sys.stderr.write(
            f"[solo_ltr_core] WARNING: {path} had {n_lines} data lines but all "
            f"were too narrow ({len(BLAST_COLS)} cols expected) -- stale/"
            f"incompatible BLAST cache? re-run with --force\n")
    return by


# --------------------------------------------------------------------------- #
# TSD detection
# --------------------------------------------------------------------------- #
def find_tsd(seq, s, e, k, slop):
    """Return the k-mer flanking [s, e) within +/-slop, or None.

    Looks for seq[s5-k:s5] == seq[e3:e3+k] over boundary offsets o5, o3 in
    [-slop, slop]; Ns disqualify a comparison. Returns the matched 5'-flank
    k-mer (the TSD sequence) on the first match found, else None.
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
                return a
    return None


# --------------------------------------------------------------------------- #
# Candidate construction
# --------------------------------------------------------------------------- #
def _build_candidate(chrom, start, end, hsp, tsd, source):
    """Assemble a Candidate from a locus span + its representative HSP."""
    qseqid = hsp.qseqid
    classification = qseqid.split("#", 1)[1] if "#" in qseqid else "."
    pdist, k2p = btop_distances(hsp.btop)
    return Candidate(
        chrom=chrom, start=start, end=end, strand=hsp.strand,
        best_consensus=qseqid, classification=classification,
        pident=hsp.pident, aln_len=hsp.length, qcov=hsp.qcovhsp,
        mismatch=hsp.mismatch, evalue=hsp.evalue, bitscore=hsp.bitscore,
        tsd=(tsd if tsd else "."), tsd_len=(len(tsd) if tsd else 0),
        pdist=pdist, k2p=k2p, source=source,
        nested_in=(hsp.host if hsp.host else "."))


def make_candidates(hsps_by_chrom, genome, *, min_pident, min_qcov,
                    min_length, max_evalue, use_tsd, tsd_k, tsd_slop,
                    merge_slop, source="main"):
    """Filter HSPs, optionally TSD-gate, then merge into loci.

    Returns dict chrom -> sorted list of Candidate. For each merged locus the
    representative is the highest-bitscore kept HSP; its metadata (best
    consensus, strand, pident, divergence, TSD, host) populates the Candidate.
    `source` labels the provenance ('main' or 'nested').
    """
    cands = {}
    for chrom, arr in hsps_by_chrom.items():
        seq = genome.get(chrom)
        kept = []  # list of (hsp, tsd_or_None)
        for h in arr:
            if (h.pident < min_pident or h.qcovhsp < min_qcov
                    or h.length < min_length or h.evalue > max_evalue):
                continue
            tsd = None
            if use_tsd:
                if seq is None:
                    continue
                tsd = find_tsd(seq, h.s0, h.e1, tsd_k, tsd_slop)
                if tsd is None:
                    continue
            kept.append((h, tsd))
        if not kept:
            continue
        kept.sort(key=lambda ht: (ht[0].s0, ht[0].e1))
        groups = [[kept[0]]]
        cur_end = kept[0][0].e1
        for ht in kept[1:]:
            if ht[0].s0 <= cur_end + merge_slop:
                groups[-1].append(ht)
                if ht[0].e1 > cur_end:
                    cur_end = ht[0].e1
            else:
                groups.append([ht])
                cur_end = ht[0].e1
        out = []
        for grp in groups:
            start = min(h.s0 for h, _ in grp)
            end = max(h.e1 for h, _ in grp)
            rep_h, rep_tsd = max(grp, key=lambda ht: ht[0].bitscore)
            out.append(_build_candidate(chrom, start, end, rep_h, rep_tsd, source))
        cands[chrom] = out
    return cands


def prefilter_internal(internal_by_chrom, *, min_pident, min_qcov,
                       min_length, max_evalue):
    """Filter + sort internal HSPs once -> dict chrom -> sorted [(s, e), ...]."""
    internal = {}
    for chrom, arr in internal_by_chrom.items():
        v = [(h.s0, h.e1) for h in arr
             if h.pident >= min_pident and h.qcovhsp >= min_qcov
             and h.length >= min_length and h.evalue <= max_evalue]
        v.sort()
        internal[chrom] = v
    return internal


def drop_near_internal(cands, internal, flank):
    """Drop Candidates with any pre-filtered internal interval within flank bp."""
    out = {}
    for chrom, arr in cands.items():
        iarr = internal.get(chrom, [])
        if not iarr:
            out[chrom] = list(arr)
            continue
        istart = [x[0] for x in iarr]
        keep = []
        for c in arr:
            s, e = c.start, c.end
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
                keep.append(c)
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
    """Remove nested Candidates that look like one LTR of a nested INTACT element.

    Keep a candidate only if every nearby internal interval cleanly contains it
    (margin >= flank both sides); drop it if some internal interval is within
    flank bp without cleanly containing it (it abuts that element's own edge).
    """
    out = {}
    for chrom, arr in cands.items():
        ivs = internal_intervals.get(chrom, [])
        keep = []
        for c in arr:
            s, e = c.start, c.end
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
                keep.append(c)
        out[chrom] = keep
    return out


TSV_HEADER = ("chrom start end name strand source best_consensus classification "
              "pident aln_len qcov mismatch evalue bitscore tsd tsd_len pdist "
              "k2p nested_in").split()


def write_tsv(cands, path, name_prefix="solo"):
    """Write Candidates as a TSV with a single '#'-commented header line.

    Columns 1-3 (chrom/start/end) are a valid BED interval, so
    `grep -v '^#' file.tsv | cut -f1-3` yields a BED and bedtools skips the
    header. Rows are sorted by chrom then start; `name` is assigned in that
    order as <name_prefix>_NNNN. Returns the row count.
    """
    n = 0
    with open(path, "w") as fh:
        fh.write("#" + "\t".join(TSV_HEADER) + "\n")
        for chrom in sorted(cands):
            for c in sorted(cands[chrom], key=lambda c: (c.start, c.end)):
                n += 1
                name = f"{name_prefix}_{n:04d}"
                # explicitly ordered to match TSV_HEADER, NOT Candidate._fields
                # (source/nested_in sit at the tuple end but are emitted early).
                fh.write("\t".join(str(x) for x in (
                    c.chrom, c.start, c.end, name, c.strand, c.source,
                    c.best_consensus, c.classification,
                    f"{c.pident:.3f}", c.aln_len, f"{c.qcov:.1f}", c.mismatch,
                    f"{c.evalue:.2e}", f"{c.bitscore:.1f}",
                    c.tsd, c.tsd_len, c.pdist, c.k2p, c.nested_in)) + "\n")
    return n


def _merge_candidate_group(grp):
    """Collapse overlapping Candidates into one: highest-bitscore representative,
    union span, nested provenance preserved if any member is nested."""
    if len(grp) == 1:
        return grp[0]
    rep = max(grp, key=lambda c: c.bitscore)
    start = min(c.start for c in grp)
    end = max(c.end for c in grp)
    # first non-'.' host wins (grp is start-sorted by union_candidates); in
    # practice one merged locus comes from a single host element.
    nested_in = next((c.nested_in for c in grp if c.nested_in != "."), ".")
    source = "nested" if any(c.source == "nested" for c in grp) else rep.source
    return rep._replace(start=start, end=end, nested_in=nested_in, source=source)


def union_candidates(*cand_dicts, merge_slop=0):
    """Combine Candidate dicts, merging intervals that overlap or sit within
    merge_slop bp. Returns a sorted, merged dict. Inputs are not mutated.
    """
    chroms = set()
    for d in cand_dicts:
        chroms.update(d)
    out = {}
    for chrom in sorted(chroms):
        members = []
        for d in cand_dicts:
            members.extend(d.get(chrom, []))
        if not members:
            continue
        members.sort(key=lambda c: (c.start, c.end))
        groups = [[members[0]]]
        cur_end = members[0].end
        for c in members[1:]:
            if c.start <= cur_end + merge_slop:
                groups[-1].append(c)
                if c.end > cur_end:
                    cur_end = c.end
            else:
                groups.append([c])
                cur_end = c.end
        out[chrom] = [_merge_candidate_group(g) for g in groups]
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
        cl = sorted((c.start, c.end) for c in cands.get(chrom, []))
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
        for c in arr:
            s, e = c.start, c.end
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
