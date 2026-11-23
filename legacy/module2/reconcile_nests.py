#!/usr/bin/env python3
"""reconcile_nests.py

Reconcile nested-LTR-RT calls across multiple rounds of ltrharvest5.py
into depth-bucketed libraries.

Each surviving element is assigned an inward-chain depth:
    chain_inward(x) = 0 if no LTR-RT is strictly inside x, else
                      1 + max(chain_inward(direct children of x)).

Direct children: y is a direct child of x iff x strictly contains y AND
no other z in the pool satisfies x strictly-contains z strictly-contains y.

Outputs {out_prefix}_depth{N}_ltr.tsv and {out_prefix}_depth{N}_ltr.fa
for every observed depth N (shadows the raw per-round files, which are
left untouched).

Usage:
  python reconcile_nests.py \
      --out-prefix mafft_update \
      --tsv mafft_update_r1_ltr.tsv mafft_update_r2_ltr.tsv ... \
      --fa  mafft_update_r1_ltr.fa  mafft_update_r2_ltr.fa  ... \
      --scn mafft_update_r1.work/mafft_update_r1.ltrtools.stitched.scn ...
"""
from __future__ import annotations

import argparse
import os
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Import shared helpers from ltrharvest5
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ltrharvest5 import (  # noqa: E402
    _is_contained,
    _ltrs_shared,
    _sub_dedup_shared_ltr_group,
    load_scn_ltr_boundaries,
    iter_fasta,
)

COORD_RE = re.compile(r"^([^:]+):(\d+)-(\d+)")

# Depth -> IUPAC mask char. Must match ltrharvest_wrapper2.sh's IUPAC_SEQ so a
# cross-round N from a round-1 element (already embedded in the extracted outer
# sequence) lines up with depth0's char here, depth1 -> R, depth2 -> D, etc.
# V is reserved for the wrapper's far-character and is NOT in this list.
IUPAC_DEPTH_SEQ = ("N", "R", "D", "Y", "S", "W", "K", "M", "B", "H")


def parse_tsv(path: str, round_idx: int) -> Tuple[Optional[str], List[dict]]:
    """Return (header_line_or_None, list_of_record_dicts)."""
    header: Optional[str] = None
    recs: List[dict] = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                if header is None:
                    header = line
                continue
            cols = line.split("\t")
            if not cols:
                continue
            m = COORD_RE.match(cols[0])
            if not m:
                continue
            chrom, s, e = m.group(1), int(m.group(2)), int(m.group(3))
            key = f"{chrom}:{s}-{e}"
            try:
                p = float(cols[6])
            except (ValueError, IndexError):
                p = 0.0
            try:
                aln = int(float(cols[2]))
            except (ValueError, IndexError):
                aln = 0
            recs.append({
                "key": key,
                "chrom": chrom,
                "s": s,
                "e": e,
                "round": round_idx,
                "col1": cols[0],
                "line": line,
                "cols": cols,
                "p": p,
                "aln": aln,
            })
    return header, recs


def build_all_in(survivors: List[dict],
                 ltr_bounds: Dict[str, Tuple[int, int, int, int]]
                 ) -> Dict[str, List[str]]:
    """For each survivor, list all survivors strictly contained by it
    (distinct LTRs required). Returns dict key -> list of descendant keys.
    """
    all_in: Dict[str, List[str]] = defaultdict(list)
    by_chr: Dict[str, List[dict]] = defaultdict(list)
    for e in survivors:
        by_chr[e["chrom"]].append(e)

    for chrom, arr in by_chr.items():
        # Sort by start asc, then end desc (outer-first on ties)
        arr.sort(key=lambda x: (x["s"], -x["e"]))
        n = len(arr)
        for i in range(n):
            x = arr[i]
            xe = x["e"]
            for j in range(i + 1, n):
                y = arr[j]
                if y["s"] > xe:
                    break  # no further overlap possible
                if y["e"] > xe:
                    continue  # partial overlap, not containment
                if _is_contained((y["s"], y["e"]), (x["s"], x["e"])) != "a_in_b":
                    continue
                if _ltrs_shared(x["key"], y["key"], ltr_bounds):
                    continue
                all_in[x["key"]].append(y["key"])
    return all_in


def build_direct_children(all_in: Dict[str, List[str]]) -> Dict[str, List[str]]:
    """y is a direct child of x iff y in all_in[x] and no z in all_in[x]
    has y in all_in[z]."""
    children: Dict[str, List[str]] = defaultdict(list)
    for xk, descendants in all_in.items():
        desc_set = set(descendants)
        for y in descendants:
            is_direct = True
            for z in descendants:
                if z == y:
                    continue
                if y in all_in.get(z, ()):
                    is_direct = False
                    break
            if is_direct:
                children[xk].append(y)
    return children


def compute_chain_inward(keys: List[str],
                         children: Dict[str, List[str]]) -> Dict[str, int]:
    """Longest chain from each node going inward (downward)."""
    memo: Dict[str, int] = {}

    def walk(k: str) -> int:
        if k in memo:
            return memo[k]
        kids = children.get(k, [])
        if not kids:
            memo[k] = 0
            return 0
        m = 1 + max(walk(c) for c in kids)
        memo[k] = m
        return m

    for k in keys:
        walk(k)
    return memo


def build_updated_nest_status(key: str,
                              all_in: Dict[str, List[str]]) -> str:
    """Pairwise nest_status string: all strict containments for `key`
    across the entire pool (both inners and outers)."""
    rels: List[Tuple[str, str]] = []
    for y in all_in.get(key, []):
        rels.append(("nest-outer", y))
    for xk, descendants in all_in.items():
        if key in descendants:
            rels.append(("nest-inner", xk))
    if not rels:
        return "."
    seen = set()
    uniq: List[Tuple[str, str]] = []
    for r in rels:
        if r not in seen:
            seen.add(r)
            uniq.append(r)
    return ";".join(f"{role}:{k}" for role, k in uniq)


def apply_depth_masking(outer_seq: str,
                        outer_rec: dict,
                        direct_children: Dict[str, List[str]],
                        rec_by_key: Dict[str, dict],
                        depth_map: Dict[str, int]) -> str:
    """Return outer_seq with every descendant region overwritten by the IUPAC
    char corresponding to that descendant's depth.

    The traversal is parent-first / grandchild-after so deeper descendants
    overwrite the shallower marks placed by their ancestors, producing the
    depth-indexed pattern:
        depth-1 outer  -> ...N... (direct child at depth 0)
        depth-2 outer  -> ...R..N..R... (direct child at depth 1 with its own
                                         depth-0 grandchild inside)
    """
    outer_s = outer_rec["s"]
    outer_chrom = outer_rec["chrom"]
    outer_len = len(outer_seq)
    chars = list(outer_seq)

    def paint(parent_key: str) -> None:
        for child_key in direct_children.get(parent_key, []):
            if child_key == parent_key:
                continue
            child = rec_by_key.get(child_key)
            if child is None or child["chrom"] != outer_chrom:
                continue
            cd = depth_map.get(child_key, 0)
            if 0 <= cd < len(IUPAC_DEPTH_SEQ):
                ch = IUPAC_DEPTH_SEQ[cd]
            else:
                ch = "X"
            # 1-based inclusive coords -> 0-based half-open relative to outer
            rel_s = max(0, child["s"] - outer_s)
            rel_e = min(outer_len, child["e"] - outer_s + 1)
            if rel_e <= rel_s:
                continue
            for p in range(rel_s, rel_e):
                chars[p] = ch
            paint(child_key)

    paint(outer_rec["key"])
    return "".join(chars)


def cross_round_dedup(pool: List[dict],
                      ltr_bounds: Dict[str, Tuple[int, int, int, int]]
                      ) -> Tuple[List[dict], int]:
    """Union-find across the pool: any two records (from different rounds)
    that share LTR boundaries are merged; keep the best per group via
    _sub_dedup_shared_ltr_group. Returns (survivors, n_merged_groups).
    """
    n = len(pool)
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    by_chr: Dict[str, List[int]] = defaultdict(list)
    for i, e in enumerate(pool):
        by_chr[e["chrom"]].append(i)

    for chrom, idxs in by_chr.items():
        m = len(idxs)
        for ii in range(m):
            i = idxs[ii]
            for jj in range(ii + 1, m):
                j = idxs[jj]
                if pool[i]["round"] == pool[j]["round"]:
                    continue
                if _ltrs_shared(pool[i]["key"], pool[j]["key"], ltr_bounds):
                    union(i, j)

    groups: Dict[int, List[int]] = defaultdict(list)
    for i in range(n):
        groups[find(i)].append(i)

    survivors: List[dict] = []
    n_merged = 0
    for members in groups.values():
        if len(members) == 1:
            survivors.append(pool[members[0]])
        else:
            n_merged += 1
            recs = [pool[m] for m in members]
            best = _sub_dedup_shared_ltr_group(recs)
            survivors.append(best)
    return survivors, n_merged


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Reconcile nested LTR-RT calls across rounds into "
                    "depth-bucketed libraries."
    )
    ap.add_argument("--tsv", nargs="+", required=True,
                    help="Per-round _ltr.tsv files (in round order).")
    ap.add_argument("--fa", nargs="+", required=True,
                    help="Per-round _ltr.fa files matching --tsv.")
    ap.add_argument("--scn", nargs="*", default=[],
                    help="Per-round merged SCN files (ltrtools.stitched.scn) "
                         "for LTR-boundary-based shared-LTR detection. "
                         "If omitted, shared-LTR check is skipped (falls back "
                         "to pure coordinate containment).")
    ap.add_argument("--out-prefix", required=True,
                    help="Output prefix for _depth{N}_ltr.{tsv,fa} files.")
    args = ap.parse_args()

    if len(args.tsv) != len(args.fa):
        ap.error("--tsv and --fa must have the same number of entries")
    if args.scn and len(args.scn) != len(args.tsv):
        ap.error("--scn count must match --tsv when provided")

    # 1. Load pool
    pool: List[dict] = []
    header: Optional[str] = None
    for r_idx, tsv in enumerate(args.tsv, start=1):
        h, recs = parse_tsv(tsv, r_idx)
        if header is None and h is not None:
            header = h
        pool.extend(recs)
    print(f"[reconcile] loaded {len(pool)} records from {len(args.tsv)} rounds",
          file=sys.stderr)

    # 2. Load LTR boundaries (union across rounds)
    ltr_bounds: Dict[str, Tuple[int, int, int, int]] = {}
    for scn in args.scn:
        if not Path(scn).exists():
            print(f"[reconcile] WARNING: SCN not found: {scn}", file=sys.stderr)
            continue
        lb = load_scn_ltr_boundaries(scn)
        for k, v in lb.items():
            ltr_bounds.setdefault(k, v)
    print(f"[reconcile] LTR boundaries loaded: {len(ltr_bounds)} keys",
          file=sys.stderr)

    # 3. Cross-round dedup (shared-LTR collapse)
    survivors, n_merged = cross_round_dedup(pool, ltr_bounds)
    if n_merged:
        print(f"[reconcile] cross-round shared-LTR merges: {n_merged} groups",
              file=sys.stderr)
    else:
        print("[reconcile] no cross-round shared-LTR merges", file=sys.stderr)

    # 4. Containment graph
    all_in = build_all_in(survivors, ltr_bounds)
    children = build_direct_children(all_in)

    # 5. Inward-chain depth
    keys = [e["key"] for e in survivors]
    depth = compute_chain_inward(keys, children)

    rec_by_key: Dict[str, dict] = {e["key"]: e for e in survivors}

    # 6. Bucket survivors; rewrite nest_status col with cross-round view
    buckets: Dict[int, List[dict]] = defaultdict(list)
    for e in survivors:
        d = depth.get(e["key"], 0)
        new_cols = list(e["cols"])
        new_cols[-1] = build_updated_nest_status(e["key"], all_in)
        e["line_out"] = "\t".join(new_cols)
        buckets[d].append(e)

    # 7. Load FASTA records keyed by full header (col1 == "chrom:s-e#class")
    fa_seqs: Dict[str, str] = {}
    for fa in args.fa:
        for h, seq in iter_fasta(fa):
            if h not in fa_seqs or len(seq) > len(fa_seqs[h]):
                fa_seqs[h] = seq

    # 8. Write outputs
    out_prefix = args.out_prefix
    any_written = False
    n_repainted = 0
    for d in sorted(buckets):
        recs = buckets[d]
        tsv_out = f"{out_prefix}_depth{d}_ltr.tsv"
        fa_out = f"{out_prefix}_depth{d}_ltr.fa"
        with open(tsv_out, "w") as tout:
            if header:
                tout.write(header + "\n")
            for r in recs:
                tout.write(r["line_out"] + "\n")
        with open(fa_out, "w") as fout:
            for r in recs:
                seq = fa_seqs.get(r["col1"])
                if seq is None:
                    continue
                if d > 0:
                    seq = apply_depth_masking(seq, r, children, rec_by_key, depth)
                    n_repainted += 1
                fout.write(f">{r['col1']}\n")
                for i in range(0, len(seq), 60):
                    fout.write(seq[i:i + 60] + "\n")
        any_written = True
        print(f"[reconcile] depth{d}: {len(recs)} -> {tsv_out}, {fa_out}",
              file=sys.stderr)
    if n_repainted:
        print(f"[reconcile] depth-indexed IUPAC masking applied to {n_repainted} "
              f"outer sequence(s)", file=sys.stderr)
    if not any_written:
        print("[reconcile] no records to write", file=sys.stderr)

    # Summary line for logs
    total = sum(len(v) for v in buckets.values())
    parts = [f"depth{d}={len(buckets[d])}" for d in sorted(buckets)]
    print(f"[reconcile] total={total} ({', '.join(parts)})", file=sys.stderr)


if __name__ == "__main__":
    main()
