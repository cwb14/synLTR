#!/usr/bin/env python3
"""
filter_solo_ltr.py

Filter mmseqs easy-search hits (LTR consensus -> hardmasked genome) down to
candidate solo-LTRs, then collapse them into unique genomic loci.

Solo-LTR true-positive rule (user definition)
---------------------------------------------
A hit is a candidate solo-LTR iff the alignment spans the LTR *consensus*
(the query) end-to-end within FLANK bp at BOTH ends, i.e. the entire LTR is
present (not a degraded fragment):

    qstart <= 1 + FLANK            (starts near the consensus start)
    qend   >= qlen - FLANK         (ends   near the consensus end)

mmseqs keeps query coordinates forward-oriented even for minus-strand target
hits (those flip tstart/tend instead), so this rule is strand-safe.

Because all intact full-length elements were hardmasked upstream, any surviving
hit is a non-intact copy; the end rule then keeps only full-length LTRs.

Inputs
------
--m8         mmseqs result, --format-mode 4 (header line) with columns:
             query target qstart qend qlen tstart tend tlen alnlen
             pident mismatch gapopen evalue bits qcov
--intact-bed merged intact-element BED (defensive overlap removal)

Two precision filters refine the full-length candidates:

Fragment purge (--intact-lib-m8): a solo-LTR carries ONLY the LTR. A
truncated/partially-excised element leaves an LTR PLUS adjacent internal
sequence. Given an mmseqs search of the FULL intact LTR-RT library vs the same
masked genome, a candidate is a 'fragment' if some library hit SPANS it
(>= --lib-cover-frac) AND that hit's element coordinates reach >=
--internal-margin bp into the element's INTERNAL region (.info).

TSD filter (--tsd-len, requires --genome): a real LTR-RT insertion duplicates
a short target site; recombination to a solo-LTR keeps that flanking direct
repeat, but library fragments do not. A candidate is kept only if an exact
L-bp direct repeat sits immediately 5' and 3' of the LTR (searched +/-
--tsd-window bp, the consensus boundary being approximate). Benchmarked vs
PrinTE truth this is the dominant precision lever (e.g. 0.06 -> ~0.8).

Outputs (<out-prefix>.*)
------------------------
  .filtered.tsv          every hit passing the end-rule (+ solo_status, tsd)
  .loci.tsv              every candidate locus (+ solo_status, lib + tsd ev.)
  .loci.bed              FINAL solo-LTR loci (purged + TSD-required), BED6
  .purged_fragments.tsv  loci removed as fragmented LTR-RTs (with evidence)
  .purged_fragments.bed  same, BED6                  [only with --intact-lib-m8]
  .summary.txt           counts + per-family breakdown

The genome is read only via indexed pyfaidx slices (TSD flanks); never loaded
whole. Only the (small) hit table is held in memory.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

COLS = [
    "query", "target", "qstart", "qend", "qlen", "tstart", "tend", "tlen",
    "alnlen", "pident", "mismatch", "gapopen", "evalue", "bits", "qcov",
]
NUMERIC = ["qstart", "qend", "qlen", "tstart", "tend", "tlen", "alnlen",
           "pident", "mismatch", "gapopen", "evalue", "bits", "qcov"]


def eprint(*a):
    print(*a, file=sys.stderr)


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Filter mmseqs hits for full-length solo-LTRs and merge to loci.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--m8", required=True,
                   help="mmseqs easy-search result (--format-mode 4, header line)")
    p.add_argument("--intact-bed", required=True,
                   help="merged intact-element BED (for defensive overlap removal)")
    p.add_argument("--out-prefix", required=True,
                   help="output path prefix (writes .filtered.tsv/.loci.bed/.loci.tsv/.summary.txt)")
    p.add_argument("--flank", type=int, default=10,
                   help="bp tolerance at EACH consensus end for the solo-LTR rule")
    p.add_argument("--min-pident", type=float, default=0.0,
                   help="optional minimum %% identity (0 = off; keeps the rule "
                        "purely positional as defined)")
    p.add_argument("--merge-dist", type=int, default=0,
                   help="merge passing hits into one locus if within this many bp "
                        "(same chrom & strand)")
    p.add_argument("--intact-lib-m8", default=None,
                   help="optional mmseqs m8 of the FULL intact LTR-RT library vs "
                        "the masked genome. If given, candidate loci that a "
                        "library hit spans AND that carries adjacent internal "
                        "element sequence are purged as truncated/fragmented "
                        "elements (not solo-LTRs).")
    p.add_argument("--intact-lib-info", default=None,
                   help="the *.fa.info file (name, ltr_start, ltr_end, "
                        "total_len; 2 LTR rows/element). Enables the rigorous "
                        "test: the library hit's ELEMENT coordinates must reach "
                        ">= --internal-margin bp into that element's INTERNAL "
                        "region (between its two LTRs) -- distinguishes 'more "
                        "LTR than the consensus' (kept) from 'real internal "
                        "sequence' (purged). Without it, a coarser genomic-"
                        "extent proxy is used (less precise; warns).")
    p.add_argument("--internal-margin", type=int, default=50,
                   help="library hit must carry >= this many bp of the "
                        "element's INTERNAL region (with --intact-lib-info), or "
                        "extend >= this many bp beyond the candidate (proxy "
                        "fallback), to trigger a fragment purge")
    p.add_argument("--lib-cover-frac", type=float, default=0.9,
                   help="a library hit must cover >= this fraction of the "
                        "candidate locus to be considered spanning it")
    p.add_argument("--genome", default=None,
                   help="genome FASTA (indexed) for de novo TSD detection")
    p.add_argument("--tsd-len", type=int, default=0,
                   help="target-site-duplication length to require flanking a "
                        "solo-LTR (0 = disabled; the pipeline sets 5). A true "
                        "solo-LTR retains the insertion TSD; library fragments "
                        "do not. Requires --genome when > 0.")
    p.add_argument("--tsd-window", type=int, default=2,
                   help="search +/- this many bp around each LTR boundary for "
                        "the flanking direct repeat (the consensus alignment "
                        "boundary is approximate)")
    p.add_argument("--tsd-mismatch", type=int, default=0,
                   help="max mismatches between the 5' and 3' TSD copies "
                        "(0 = exact; PrinTE-style TSDs are exact duplications)")
    p.add_argument("--verbose", action="store_true",
                   help="per-step progress and sanity checks")
    return p.parse_args(argv)


def read_m8(path: str) -> pd.DataFrame:
    """Read mmseqs format-mode-4 output (header line) or headerless BLAST-tab."""
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame(columns=COLS)
    with path.open() as fh:
        first = fh.readline().rstrip("\n")
    has_header = first.split("\t")[0].strip().lower() == "query"
    df = pd.read_csv(
        path, sep="\t",
        header=0 if has_header else None,
        names=None if has_header else COLS,
        dtype=str, comment=None,
    )
    # Normalise column names (mmseqs header uses the requested names already).
    df.columns = [str(c).strip() for c in df.columns]
    missing = [c for c in COLS if c not in df.columns]
    if missing:
        raise SystemExit(
            f"ERROR: {path} missing expected columns {missing}. "
            f"Re-run mmseqs with the documented --format-output.")
    for c in NUMERIC:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def family_of(query: str) -> str:
    """'chr2:..-..#LTR/Gypsy/mixture' -> 'LTR/Gypsy/mixture'; else 'NA'."""
    return query.split("#", 1)[1] if "#" in query else "NA"


def load_intact(path: str) -> dict[str, list[tuple[int, int]]]:
    """chrom -> sorted list of (start0, end) intervals (already merged upstream)."""
    by_chrom: dict[str, list[tuple[int, int]]] = {}
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return by_chrom
    with p.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            by_chrom.setdefault(f[0], []).append((int(f[1]), int(f[2])))
    for c in by_chrom:
        by_chrom[c].sort()
    return by_chrom


def overlaps_intact(chrom, s0, e, intact) -> bool:
    """True if [s0, e) overlaps any intact interval on chrom (linear scan;
    interval counts here are small)."""
    for a, b in intact.get(chrom, ()):
        if a < e and s0 < b:
            return True
        if a >= e:
            break
    return False


def annotate_tsd(df: pd.DataFrame, genome_path: str, L: int, W: int,
                 mm: int, verbose: bool) -> pd.DataFrame:
    """Add tsd_ok / tsd_seq per hit by de-novo scanning the genome flanks.

    A real LTR-RT insertion duplicates an L-bp target site; recombination to
    a solo-LTR keeps that flanking direct repeat, while library fragments do
    not. The consensus-alignment boundary is approximate, so search +/- W bp:
    a hit is TSD-supported if some boundary pair has identical (<= mm mismatch)
    L-mers immediately 5' and 3' of the LTR.
    """
    from pyfaidx import Fasta
    g = Fasta(genome_path, as_raw=True, sequence_always_upper=True)
    clen = {c: len(g[c]) for c in g.keys()}

    def scan(chrom, gs, ge):
        cl = clen.get(chrom)
        if cl is None or gs - W - L < 0 or ge + W + L > cl:
            return False, ""
        up = g[chrom][gs - W - L: gs + W]      # 2W+L bp around 5' boundary
        dn = g[chrom][ge - W: ge + W + L]      # 2W+L bp around 3' boundary
        for j in range(0, 2 * W + 1):
            t5 = up[j:j + L]
            if len(t5) < L or "N" in t5:
                continue
            for k in range(0, 2 * W + 1):
                t3 = dn[k:k + L]
                if len(t3) < L or "N" in t3:
                    continue
                if mm == 0:
                    if t5 == t3:
                        return True, t5
                elif sum(a != b for a, b in zip(t5, t3)) <= mm:
                    return True, t5
        return False, ""

    res = [scan(c, int(s), int(e))
           for c, s, e in zip(df["chrom"], df["gstart"], df["gend"])]
    df = df.copy()
    df["tsd_ok"] = [r[0] for r in res]
    df["tsd_seq"] = [r[1] for r in res]
    if verbose:
        n = len(df)
        eprint(f"[filter] TSD (L={L} W={W} mm={mm}): "
               f"{int(df['tsd_ok'].sum())}/{n} hits carry a flanking TSD")
    return df


def merge_loci(df: pd.DataFrame, merge_dist: int) -> pd.DataFrame:
    """Collapse passing hits into unique loci per (chrom, strand).

    Representative per locus = hit with the highest bit score. Also records
    support breadth (#hits, #distinct consensi, distinct families)."""
    if df.empty:
        return pd.DataFrame(columns=[
            "chrom", "start", "end", "strand", "n_hits", "n_consensi",
            "families", "best_query", "best_pident", "best_evalue",
            "best_bits", "tsd_supported", "tsd_seq"])
    has_tsd = "tsd_ok" in df.columns
    rows = []
    df = df.sort_values(["chrom", "strand", "gstart", "gend"])
    for (chrom, strand), g in df.groupby(["chrom", "strand"], sort=False):
        cs = ce = None
        bucket = []

        def flush(bucket, cs, ce):
            best = max(bucket, key=lambda r: r["bits"])
            # a locus is TSD-supported if ANY of its hits has a flanking TSD
            tsd_hits = [r for r in bucket if has_tsd and r.get("tsd_ok")]
            return {
                "chrom": chrom, "start": cs, "end": ce, "strand": strand,
                "n_hits": len(bucket),
                "n_consensi": len({r["query"] for r in bucket}),
                "families": ";".join(sorted({r["family"] for r in bucket})),
                "best_query": best["query"],
                "best_pident": round(best["pident"], 3),
                "best_evalue": best["evalue"],
                "best_bits": best["bits"],
                "tsd_supported": bool(tsd_hits),
                "tsd_seq": (tsd_hits[0]["tsd_seq"] if tsd_hits else ""),
            }

        for r in g.to_dict("records"):
            if cs is None:
                cs, ce, bucket = r["gstart"], r["gend"], [r]
            elif r["gstart"] <= ce + merge_dist:
                ce = max(ce, r["gend"])
                bucket.append(r)
            else:
                rows.append(flush(bucket, cs, ce))
                cs, ce, bucket = r["gstart"], r["gend"], [r]
        if bucket:
            rows.append(flush(bucket, cs, ce))
    return pd.DataFrame(rows).sort_values(["chrom", "start", "end"])


def load_internal(path: str) -> dict:
    """Parse *.fa.info -> {element_name: (internal_start, internal_end)} in
    1-based-inclusive ELEMENT coordinates.

    Each element has two LTR rows (name, ltr_start, ltr_end, total_len). The
    internal region is between them: [min-LTR.end + 1, max-LTR.start - 1].
    Elements without exactly two valid LTR rows are omitted -> their library
    hits can never trigger an internal-based purge (conservative: keep)."""
    rows = {}
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return {}
    with p.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            try:
                a, b = int(f[1]), int(f[2])
            except ValueError:
                continue
            rows.setdefault(f[0], []).append((min(a, b), max(a, b)))
    internal = {}
    for name, ltrs in rows.items():
        if len(ltrs) != 2:
            continue
        ltrs.sort()
        ia, ib = ltrs[0][1] + 1, ltrs[1][0] - 1     # between the two LTRs
        if ia <= ib:
            internal[name] = (ia, ib)
    return internal


def library_intervals(path: str, internal: dict | None):
    """Per-chromosome arrays of intact-library hits for fast overlap queries.

    In this search the library element is the QUERY and the genome chromosome
    the TARGET. Returns {chrom: (gs, ge, internal_bp, qname)} sorted by gs:
      gs, ge       genomic interval, 0-based half-open
      internal_bp  bp of the hit that lies in the element's INTERNAL region
                   (0 if --intact-lib-info absent or element boundary unknown)
      qname        library element (for purge evidence)
    """
    import numpy as np
    df = read_m8(path)
    out = {}
    if df.empty:
        return out
    df = df.dropna(subset=["tstart", "tend", "qstart", "qend"]).copy()
    gs = df[["tstart", "tend"]].min(axis=1).astype(int) - 1   # genomic 0-based
    ge = df[["tstart", "tend"]].max(axis=1).astype(int)       # exclusive
    qlo = df[["qstart", "qend"]].min(axis=1).astype(int)      # element coords
    qhi = df[["qstart", "qend"]].max(axis=1).astype(int)
    qn = df["query"].astype(str)
    if internal:
        # element unknown -> sentinel (1, 0): internal overlap evaluates to 0
        ia = qn.map(lambda x: internal.get(x, (1, 0))[0]).to_numpy(dtype=int)
        ib = qn.map(lambda x: internal.get(x, (1, 0))[1]).to_numpy(dtype=int)
        # bp of [qlo,qhi] (incl.) overlapping the element's internal [ia,ib]
        ibp = np.clip(np.minimum(qhi.to_numpy(), ib)
                      - np.maximum(qlo.to_numpy(), ia) + 1, 0, None)
    else:
        ibp = np.zeros(len(df), dtype=int)
    tmp = pd.DataFrame({"chrom": df["target"].astype(str),
                        "s": gs.values, "e": ge.values,
                        "ibp": np.asarray(ibp, dtype=int),
                        "q": qn.values})
    for chrom, gg in tmp.groupby("chrom", sort=False):
        gg = gg.sort_values("s")
        out[chrom] = (gg["s"].to_numpy(), gg["e"].to_numpy(),
                      gg["ibp"].to_numpy(), gg["q"].to_numpy())
    return out


def classify_loci(loci: pd.DataFrame, lib, margin: int, cover_frac: float,
                   use_info: bool, verbose: bool) -> pd.DataFrame:
    """Tag each candidate locus 'solo' or 'fragment'.

    A locus is a 'fragment' (truncated / partially-excised element, NOT a
    solo-LTR) if some intact-library hit SPANS the candidate (covers
    >= cover_frac of it) AND, on that same hit, carries >= `margin` bp of the
    element's INTERNAL region (with --intact-lib-info; the biologically exact
    test) -- or, as a coarser fallback without .info, extends >= `margin` bp
    beyond the candidate genomically. A bona fide solo-LTR has only the LTR,
    so library hits there carry no internal element sequence.
    """
    import numpy as np
    n = len(loci)
    status = ["solo"] * n
    ev_bp = [0] * n
    ev_q = [""] * n
    for i, (_, L) in enumerate(loci.reset_index(drop=True).iterrows()):
        arrs = lib.get(L["chrom"])
        if arrs is None:
            continue
        s, e, ibp, q = arrs
        ls, le = int(L["start"]), int(L["end"])
        clen = max(le - ls, 1)
        hi = int(np.searchsorted(s, le, side="right"))   # start < locus end
        if hi == 0:
            continue
        ss, ee, ii, qq = s[:hi], e[:hi], ibp[:hi], q[:hi]
        ov = np.minimum(le, ee) - np.maximum(ls, ss)             # overlap bp
        spans = ov >= cover_frac * clen                           # covers cand.
        if not spans.any():
            continue
        if use_info:
            evidence = ii                       # internal bp carried by hit
        else:
            evidence = np.maximum(ls - ss, ee - le)   # genomic overhang proxy
        trig = spans & (evidence >= margin)
        if trig.any():
            j = int(np.argmax(np.where(trig, evidence, -1)))
            status[i] = "fragment"
            ev_bp[i] = int(evidence[j])
            ev_q[i] = str(qq[j])
    loci = loci.copy()
    loci["solo_status"] = status
    loci["lib_internal_bp" if use_info else "lib_extend_bp"] = ev_bp
    loci["lib_evidence_query"] = ev_q
    if verbose:
        nf = status.count("fragment")
        mode = "internal-overlap (.info)" if use_info else "genomic-extent PROXY"
        eprint(f"[filter] library purge [{mode}]: {nf}/{n} loci reclassified "
               f"as fragmented LTR-RT (margin={margin}, cover>={cover_frac})")
    return loci


def main(argv=None):
    args = parse_args(argv)
    flank = args.flank
    df = read_m8(args.m8)
    n_raw = len(df)
    if args.verbose:
        eprint(f"[filter] raw hits: {n_raw}")

    out = Path(args.out_prefix)
    out.parent.mkdir(parents=True, exist_ok=True)

    if n_raw == 0:
        for ext in ("filtered.tsv", "loci.tsv"):
            (out.parent / f"{out.name}.{ext}").write_text("")
        (out.parent / f"{out.name}.loci.bed").write_text("")
        (out.parent / f"{out.name}.summary.txt").write_text(
            "No mmseqs hits to filter (0 raw hits).\n")
        eprint("[filter] no hits; wrote empty outputs")
        return 0

    df = df.dropna(subset=["qstart", "qend", "qlen", "tstart", "tend"]).copy()
    for c in ("qstart", "qend", "qlen", "tstart", "tend"):
        df[c] = df[c].astype(int)

    # mmseqs 'pident' may be a fraction (0-1) or a percentage (0-100) depending
    # on build; normalise to percentage so --min-pident is always in percent.
    pmax = df["pident"].max()
    if pd.notna(pmax) and pmax <= 1.0:
        if args.verbose:
            eprint("[filter] pident looks fractional (max<=1) -> scaling x100")
        df["pident"] = df["pident"] * 100.0

    # mmseqs --search-type 3 encodes a minus-strand hit by REVERSING the QUERY
    # coordinates (qstart > qend) while keeping target coords ascending
    # (verified empirically: 0/444 rows had tstart>tend; 236/444 had
    # qstart>qend; a reverse-complement probe reproduced it). So the consensus
    # span covered by a hit is [min(qstart,qend), max(qstart,qend)] regardless
    # of orientation; applying the end-rule to raw qstart/qend would silently
    # reject EVERY minus-strand solo-LTR.
    qlo = df[["qstart", "qend"]].min(axis=1)
    qhi = df[["qstart", "qend"]].max(axis=1)

    # ---- solo-LTR rule: full-length over the consensus (query), strand-agnostic
    rule = (qlo <= 1 + flank) & (qhi >= df["qlen"] - flank)
    passed = df[rule].copy()
    if args.verbose:
        eprint(f"[filter] pass end-rule (qlo<=1+{flank} & qhi>=qlen-{flank}, "
               f"strand-agnostic): {len(passed)} / {n_raw}")

    if args.min_pident > 0:
        before = len(passed)
        passed = passed[passed["pident"] >= args.min_pident].copy()
        if args.verbose:
            eprint(f"[filter] pass min-pident>={args.min_pident}: "
                   f"{len(passed)} / {before}")

    # Derived genomic coords (BED 0-based half-open) + strand. Strand is encoded
    # by QUERY coordinate order (qstart>qend => minus); target coords are always
    # ascending so the genomic interval is [tstart-1, tend) (min/max kept for
    # robustness against any build that flips them).
    passed["chrom"] = passed["target"].astype(str)
    passed["strand"] = ["-" if qs > qe else "+"
                        for qs, qe in zip(passed["qstart"], passed["qend"])]
    passed["gstart"] = passed[["tstart", "tend"]].min(axis=1) - 1   # 0-based
    passed["gend"] = passed[["tstart", "tend"]].max(axis=1)         # exclusive
    passed["family"] = passed["query"].map(family_of)

    # Defensive: drop anything overlapping an intact (masked) element. Masked
    # bases are N so this should be ~empty; it catches edge bleed at element
    # boundaries and proves the masking did its job.
    intact = load_intact(args.intact_bed)
    if intact:
        keep = [not overlaps_intact(c, s, e, intact)
                for c, s, e in zip(passed["chrom"], passed["gstart"], passed["gend"])]
        n_drop = len(keep) - sum(keep)
        if n_drop and args.verbose:
            eprint(f"[filter] dropped {n_drop} hit(s) overlapping intact mask "
                   f"(expected ~0)")
        passed = passed[keep].copy()

    # ---- de-novo TSD annotation (precision filter) ----
    tsd_enabled = args.tsd_len > 0
    if tsd_enabled and not args.genome:
        raise SystemExit("ERROR: --tsd-len > 0 requires --genome FASTA "
                         "(use --tsd-len 0 to disable the TSD filter)")
    if tsd_enabled and not passed.empty:
        passed = annotate_tsd(passed, args.genome, args.tsd_len,
                              args.tsd_window, args.tsd_mismatch, args.verbose)

    # ---- merge into unique loci ----
    loci = merge_loci(passed, args.merge_dist)

    # ---- optional: purge fragmented/truncated elements using the intact
    #      LTR-RT library search (a true solo-LTR carries ONLY the LTR) ----
    use_lib = bool(args.intact_lib_m8)
    internal = load_internal(args.intact_lib_info) if args.intact_lib_info else {}
    use_info = bool(internal)
    ev_col = "lib_internal_bp" if use_info else "lib_extend_bp"
    if use_lib and not internal and args.intact_lib_m8:
        eprint("[filter] WARNING: no --intact-lib-info (or unparyable) -> "
               "using the coarser genomic-extent proxy; this can FALSE-purge "
               "true solo-LTRs whose genomic LTR is longer than the consensus. "
               "Pass the *.fa.info file for the exact internal-overlap test.")
    if use_lib and not loci.empty:
        lib = library_intervals(args.intact_lib_m8, internal)
        loci = classify_loci(loci, lib, args.internal_margin,
                             args.lib_cover_frac, use_info, args.verbose)
    else:
        loci = loci.copy()
        loci["solo_status"] = "solo" if not loci.empty else []
        loci[ev_col] = 0 if not loci.empty else []
        loci["lib_evidence_query"] = "" if not loci.empty else []

    # Assign each passing hit to its locus to carry solo_status onto
    # filtered.tsv. Loci are non-overlapping unions of these hits within a
    # (chrom,strand), sorted by start -> bisect for an O(log n) exact lookup.
    import bisect as _bisect

    def assign_status(p):
        if p.empty:
            p = p.copy(); p["solo_status"] = []; return p
        idx = {}
        for _, L in loci.iterrows():
            k = (L["chrom"], L["strand"])
            idx.setdefault(k, []).append((L["start"], L["end"],
                                          L["solo_status"]))
        for k in idx:
            idx[k].sort()
        st = []
        for c, s, e, sd in zip(p["chrom"], p["gstart"], p["gend"],
                               p["strand"]):
            lst = idx.get((c, sd))
            if not lst:
                st.append("solo"); continue
            starts = [t[0] for t in lst]
            j = _bisect.bisect_right(starts, s) - 1   # last start <= s
            st.append(lst[j][2] if j >= 0 and lst[j][1] >= e else "solo")
        p = p.copy(); p["solo_status"] = st
        return p
    passed = assign_status(passed)

    # ---- write filtered hits (+ solo_status, + tsd if enabled) ----
    filt_cols = ["query", "family", "chrom", "gstart", "gend", "strand",
                 "qstart", "qend", "qlen", "tstart", "tend", "alnlen",
                 "pident", "mismatch", "gapopen", "evalue", "bits", "qcov",
                 "solo_status"]
    if tsd_enabled:
        filt_cols += ["tsd_ok", "tsd_seq"]
    passed[filt_cols].sort_values(["chrom", "gstart", "gend"]).to_csv(
        f"{args.out_prefix}.filtered.tsv", sep="\t", index=False)

    # ---- loci tables: all (annotated), final solo-only BED, purged BED ----
    loci.to_csv(f"{args.out_prefix}.loci.tsv", sep="\t", index=False)
    if loci.empty:
        solo = frag = no_tsd = loci
    else:
        is_solo = loci["solo_status"] == "solo"
        frag = loci[~is_solo]
        if tsd_enabled:
            solo = loci[is_solo & loci["tsd_supported"]]
            no_tsd = loci[is_solo & ~loci["tsd_supported"]]
        else:
            solo = loci[is_solo]
            no_tsd = loci.iloc[0:0]

    def write_bed(sub, path, extra=None):
        if sub.empty:
            Path(path).write_text(""); return
        b = sub.copy()
        b["name"] = b["best_query"] + "|" + b["families"]
        if extra == "frag":
            b["name"] = (b["name"] + "|" + ev_col + "=" + b[ev_col].astype(str)
                         + "|" + b["lib_evidence_query"])
        b["score"] = b["best_bits"].round().astype(int).clip(0, 1000)
        b[["chrom", "start", "end", "name", "score", "strand"]].to_csv(
            path, sep="\t", index=False, header=False)

    write_bed(solo, f"{args.out_prefix}.loci.bed")          # FINAL solo-LTRs
    if use_lib:
        write_bed(frag, f"{args.out_prefix}.purged_fragments.bed", "frag")
        frag.to_csv(f"{args.out_prefix}.purged_fragments.tsv",
                    sep="\t", index=False)

    # ---- summary ----
    lines = [
        "solo-LTR candidate filtering summary",
        "=" * 38,
        f"m8 input              : {args.m8}",
        f"raw hits              : {n_raw}",
        f"flank (each end, bp)  : {flank}",
        f"min %identity         : {args.min_pident if args.min_pident>0 else 'off'}",
        f"merge distance (bp)   : {args.merge_dist}",
        f"hits passing rule     : {len(passed)}",
        f"distinct consensi     : {passed['query'].nunique() if len(passed) else 0}",
        f"candidate loci        : {len(loci)}",
    ]
    if use_lib:
        lines += [
            f"intact-lib m8         : {args.intact_lib_m8}",
            f"purge mode            : "
            f"{'internal-overlap (.info)' if use_info else 'genomic-extent PROXY'}",
            f"internal margin (bp)  : {args.internal_margin}",
            f"lib cover fraction    : {args.lib_cover_frac}",
            f"purged as fragments   : {len(frag)}",
        ]
    if tsd_enabled:
        lines += [
            f"TSD filter            : L={args.tsd_len} W={args.tsd_window} "
            f"mm={args.tsd_mismatch}",
            f"solo, no TSD (removed): {len(no_tsd)}",
        ]
    lines.append(f"FINAL solo-LTR loci   : {len(solo)}"
                 + ("" if (use_lib or tsd_enabled)
                    else " (no library purge / TSD applied)"))
    lines += [
        "",
        "final solo-LTR loci per family (class/family of best hit):",
    ]
    if not solo.empty:
        fam = (solo["best_query"].map(family_of)
               .value_counts().sort_values(ascending=False))
        for k, v in fam.items():
            lines.append(f"  {k:<32} {v}")
    else:
        lines.append("  (none)")
    lines.append("")
    Path(f"{args.out_prefix}.summary.txt").write_text("\n".join(lines) + "\n")
    eprint("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
