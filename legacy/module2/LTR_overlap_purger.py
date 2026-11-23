#!/usr/bin/env python3
import sys
import argparse
import re
from collections import defaultdict, namedtuple, deque

Row = namedtuple("Row", ["idx","raw","chrom","start","end","ltr_len","div","int_start","int_end"])

coord_re = re.compile(r"^(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)$")
mixture_re = re.compile(r"\bmixture\b", re.IGNORECASE)

def parse_row(line, idx):
    # Split on tabs; tolerate extra whitespace
    parts = line.rstrip("\n").split("\t")
    if len(parts) < 2:
        raise ValueError(f"Line {idx+1}: expected at least 2 tab-separated columns")
    coord = parts[0].strip()
    m = coord_re.match(coord)
    if not m:
        raise ValueError(f"Line {idx+1}: bad coordinate '{coord}'")
    chrom = m.group("chrom")
    start = int(m.group("start"))
    end   = int(m.group("end"))
    try:
        ltr_len = int(parts[1].strip())
    except Exception:
        raise ValueError(f"Line {idx+1}: column2 (LTR length) not an integer: '{parts[1]}'")
    # divergence (column 7) may be missing; default high so it loses on --ov-lowest-div
    div = None
    if len(parts) >= 7 and parts[6].strip() != "":
        try:
            div = float(parts[6].strip())
        except Exception:
            div = None
    if div is None:
        div = float("inf")
    # internal interval
    int_start = start + ltr_len
    int_end   = end   - ltr_len
    return Row(idx, line.rstrip("\n"), chrom, start, end, ltr_len, div, int_start, int_end)

def shared_groups(rows, key):
    """Group rows by a key on a per-chrom basis (start or end)."""
    groups = defaultdict(list)
    for r in rows:
        groups[(r.chrom, key(r))].append(r)
    # Keep only groups with size > 1
    return {k:v for k,v in groups.items() if len(v) > 1}

def mark_shared_start_end(rows):
    """Apply rule (1). Return sets: keep_ids, drop_ids and a reason dict."""
    keep, drop, reason = set(), set(), {}
    # Shared start
    for (_, _), group in shared_groups(rows, key=lambda r: r.start).items():
        max_len = max(r.ltr_len for r in group)
        max_group = [r for r in group if r.ltr_len == max_len]
        if len(max_group) == 1:
            winner = max_group[0].idx
            if winner not in drop:
                keep.add(winner)
                reason[winner] = "keep_shared_start_longest"
            for r in group:
                if r.idx != winner and r.idx not in keep:
                    drop.add(r.idx)
                    reason[r.idx] = "drop_shared_start_shorter"
        else:
            # tie for longest → drop all
            for r in group:
                drop.add(r.idx)
                reason[r.idx] = "drop_shared_start_tie_longest"

    # Shared end
    for (_, _), group in shared_groups(rows, key=lambda r: r.end).items():
        max_len = max(r.ltr_len for r in group)
        max_group = [r for r in group if r.ltr_len == max_len]
        if len(max_group) == 1:
            winner = max_group[0].idx
            if winner not in drop:
                keep.add(winner)
                reason.setdefault(winner, "keep_shared_end_longest")
            for r in group:
                if r.idx != winner and r.idx not in keep:
                    drop.add(r.idx)
                    reason.setdefault(r.idx, "drop_shared_end_shorter")
        else:
            for r in group:
                drop.add(r.idx)
                reason.setdefault(r.idx, "drop_shared_end_tie_longest")

    return keep, drop, reason

def interval_overlap(a_start, a_end, b_start, b_end):
    return not (a_end <= b_start or b_end <= a_start)

def is_nested_inner(inner: Row, outer: Row):
    """True if inner.outer interval fully inside outer.internal interval."""
    return (inner.chrom == outer.chrom and
            inner.start >= outer.int_start and
            inner.end   <= outer.int_end)

def find_protected_nested(rows, excluded_ids):
    """Rule (2): identify pairs that are valid nestings; protect all involved."""
    protected = set()
    by_chrom = defaultdict(list)
    for r in rows:
        if r.idx in excluded_ids:
            continue
        by_chrom[r.chrom].append(r)
    for chrom, lst in by_chrom.items():
        lst.sort(key=lambda r: r.start)
        for i in range(len(lst)):
            a = lst[i]
            for j in range(i+1, len(lst)):
                b = lst[j]
                if b.start > a.end:
                    break  # nothing further overlaps a
                # Check nesting both ways
                if is_nested_inner(a, b) or is_nested_inner(b, a):
                    protected.add(a.idx)
                    protected.add(b.idx)
    return protected

def build_type3_components(rows, excluded_ids, protected_ids):
    """Build components of type-3 overlaps (non-nested, non-shared)."""
    candidates = [r for r in rows if r.idx not in excluded_ids]
    by_chrom = defaultdict(list)
    for r in candidates:
        by_chrom[r.chrom].append(r)
    edges = defaultdict(set)
    for chrom, lst in by_chrom.items():
        lst.sort(key=lambda r: r.start)
        active = []
        for r in lst:
            active = [a for a in active if a.end > r.start]
            for a in active:
                if interval_overlap(a.start, a.end, r.start, r.end):
                    if a.start == r.start or a.end == r.end:
                        continue  # shared start/end handled earlier
                    if is_nested_inner(a, r) or is_nested_inner(r, a):
                        continue  # nesting handled separately
                    edges[a.idx].add(r.idx)
                    edges[r.idx].add(a.idx)
            active.append(r)
    seen = set()
    components = []
    for node in edges:
        if node in seen:
            continue
        comp = set()
        dq = deque([node])
        seen.add(node)
        while dq:
            u = dq.popleft()
            comp.add(u)
            for v in edges[u]:
                if v not in seen:
                    seen.add(v)
                    dq.append(v)
        components.append(comp)
    components = [c for c in components if len(c) >= 2]
    filtered = []
    for comp in components:
        c2 = {i for i in comp if i not in protected_ids}
        if len(c2) >= 2:
            filtered.append(c2)
    return filtered

def resolve_component(comp_ids, rows_by_id, policy):
    """Return sets (keep, drop) within the component according to policy."""
    if policy == "retain":
        return set(comp_ids), set()
    if policy == "drop":
        return set(), set(comp_ids)
    comp = [rows_by_id[i] for i in comp_ids]
    if policy == "longest":
        max_len = max(r.ltr_len for r in comp)
        keep = {r.idx for r in comp if r.ltr_len == max_len}
        drop = set(comp_ids) - keep
        return keep, drop
    if policy == "lowest_div":
        min_div = min(r.div for r in comp)
        keep = {r.idx for r in comp if r.div == min_div}
        drop = set(comp_ids) - keep
        return keep, drop
    raise ValueError("unknown policy")

def main():
    ap = argparse.ArgumentParser(
        description="Purge LTR-RT duplicates/false positives by coordinate overlap."
    )
    ap.add_argument("infile", nargs="?", default="-",
                    help="Input TSV (default: stdin)")
    ap.add_argument("--ov-longest", dest="ov", action="store_const", const="longest",
                    help="For non-nested, non-shared overlaps keep longest (column2).")
    ap.add_argument("--ov-lowest-div", dest="ov", action="store_const", const="lowest_div",
                    help="For non-nested, non-shared overlaps keep lowest divergence (column7).")
    ap.add_argument("--ov-drop-all", dest="ov", action="store_const", const="drop",
                    help="For non-nested, non-shared overlaps drop all.")
    ap.add_argument("--ov-retain-all", dest="ov", action="store_const", const="retain",
                    help="For non-nested, non-shared overlaps retain all.")
    ap.add_argument("--log", default=None, metavar="FILE",
                    help="Optional: write a drop/keep log with reasons.")
    ap.set_defaults(ov="longest")
    args = ap.parse_args()

    # Read raw lines
    if args.infile == "-" or args.infile == "/dev/stdin":
        raw_lines = [l for l in sys.stdin if l.strip() != ""]
    else:
        with open(args.infile) as fh:
            raw_lines = [l for l in fh if l.strip() != ""]

    # (0) Pre-filter: drop lines containing the word 'mixture' (case-insensitive)
    # This happens BEFORE all other logic.
    mixture_dropped_idx = set()
    clean_lines = []
    for i, line in enumerate(raw_lines):
        if line.lstrip().startswith("#"):
            clean_lines.append(line)  # keep comments
            continue
        if mixture_re.search(line):
            mixture_dropped_idx.add(i)
        else:
            clean_lines.append(line)

    # Parse remaining (non-mixture) data lines
    rows = []
    idx_map = {}  # map new sequential index to original raw line index (for logging consistency)
    new_idx = 0
    for i, line in enumerate(clean_lines):
        if line.lstrip().startswith("#"):
            continue
        row = parse_row(line, new_idx)
        rows.append(row)
        idx_map[new_idx] = i
        new_idx += 1

    rows_by_id = {r.idx: r for r in rows}

    # (1) Shared start / end
    keep1, drop1, reason = mark_shared_start_end(rows)

    # (2) Valid nesting → keep both
    protected = find_protected_nested(rows, excluded_ids=drop1)

    # (3) Non-nested, non-shared overlaps with policy
    comps = build_type3_components(rows, excluded_ids=drop1, protected_ids=protected)

    keep3 = set()
    drop3 = set()
    for comp in comps:
        k, d = resolve_component(comp, rows_by_id, policy=args.ov)
        keep3 |= k
        drop3 |= d
        for i in d:
            reason.setdefault(i, f"drop_type3_{args.ov}")
        for i in k:
            reason.setdefault(i, f"keep_type3_{args.ov}")

    # Collate final sets
    dropped = set()
    dropped |= drop1
    dropped |= drop3
    kept = {r.idx for r in rows if r.idx not in dropped}
    kept |= protected  # ensure protected kept
    for i in protected:
        reason.setdefault(i, "keep_nested")

    # Output kept rows in the current (parsed) order
    for i in sorted(kept):
        print(rows_by_id[i].raw)

    # Optional log (includes mixture-dropped and comments-awareness)
    if args.log:
        with open(args.log, "w") as lf:
            print("#row_seq\tstatus\treason\tchrom\tstart\tend\tLTRlen\tdiv\tline_example", file=lf)
            # First, report mixture-dropped lines
            for i, line in enumerate(raw_lines):
                if i in mixture_dropped_idx and not line.lstrip().startswith("#"):
                    # Try to parse minimal fields for context
                    chrom = start = end = ltr_len = div = "NA"
                    try:
                        parts = line.rstrip("\n").split("\t")
                        m = coord_re.match(parts[0].strip())
                        if m:
                            chrom = m.group("chrom")
                            start = m.group("start")
                            end   = m.group("end")
                        if len(parts) >= 2:
                            ltr_len = parts[1].strip()
                        if len(parts) >= 7:
                            div = parts[6].strip()
                    except Exception:
                        pass
                    print(f"{i}\tdropped\tdrop_contains_mixture\t{chrom}\t{start}\t{end}\t{ltr_len}\t{div}\t{line.strip()}", file=lf)

            # Then, report all parsed rows
            for r in sorted(rows, key=lambda x: x.idx):
                status = "kept" if r.idx in kept else "dropped"
                rsn = reason.get(r.idx, "kept_default" if status=="kept" else "drop_default")
                print(f"{r.idx}\t{status}\t{rsn}\t{r.chrom}\t{r.start}\t{r.end}\t{r.ltr_len}\t{r.div}\t{rows_by_id[r.idx].raw}", file=lf)

if __name__ == "__main__":
    main()
