#!/usr/bin/env python3
"""ltr_annotate.py

Add `strand` and `family` columns to the depth-bucketed LTR-RT tables written by
reconcile_nests.py (and FP-purged by flag_fp_families.py).

Both {prefix}_depth{N}_ltr.tsv and {prefix}_depth{N}_clean_ltr.tsv are rewritten
in place, gaining two columns immediately after `tsd`:

    ... 16 ltr3_start | 17 tsd | 18 strand | 19 family | 20 domains | 21 nest_status

The columns go *before* `domains`/`nest_status` on purpose. ltrharvest_plot_struct.py
and TEGV.py both read those two fields from the end of the row (parts[-2], parts[-1]);
appending would silently feed them the new columns instead, emptying the domain track
in the structure PDFs and the TEGV viewer. This placement also leaves every positional
index intact for consumers that count from the front (cols[10] == K2P_d, etc).

Strand is resolved by a four-tier cascade, most authoritative first:

  1. tesorter     TEsorter2's own call, column 6 of {prefix}_r*.work/*.cls.tsv
  2. domain_order inferred from the genomic order of protein domains
  3. pass2        inherited from the minimap2 pass-2 homology target
  4. .            unknown

Family labels ({prefix}_fam00001, ...) come from the Kmer2LTR/mmseqs consensus
cluster table, numbered by descending family size.

Usage:
  python ltr_annotate.py --prefix Athal_tair10_chr2_LTRs [--indir .] [-v]

Missing inputs degrade rather than abort: with no .work directories every element
gets strand '.', with no cluster table every element gets family '.'. A non-zero
exit means a real fault. Re-running is safe -- existing strand/family columns are
replaced, not duplicated.
"""
from __future__ import annotations

import argparse
import glob
import os
import re
import stat
import sys
import tempfile
from collections import defaultdict
from typing import Dict, List, NamedTuple, Optional, Sequence, Tuple

# Columns this script manages, and the column they are inserted after.
STRAND_COL = "strand"
FAMILY_COL = "family"
ANCHOR_COL = "tsd"

UNKNOWN = "."

# The wrapper hard-caps detection at 10 rounds, so depth cannot exceed 10.
# Probing numerically (rather than globbing) keeps depth0, depth1, ... in
# numeric order regardless of glob lexicography.
MAX_DEPTH = 20

# Pass-2 strand transfer is iterated so a chain of homology transfers resolves.
# In practice it converges after one hop; the cap is a safety net.
MAX_PASS2_ITERS = 10

FLIP = {"+": "-", "-": "+"}

ELEMENT_RE = re.compile(r"^(.+):(\d+)-(\d+)$")
# 'GENE|clade@gStart-gEnd' -- the token grammar of the TSV `domains` column.
DOMAIN_TOKEN_RE = re.compile(r"^([^|@;]+)\|([^@;]*)@(\d+)-(\d+)$")

# Canonical 5'->3' internal-domain order per superfamily. Mirrors
# ltrharvest_plot_struct.py's _CANON_RANK_* tables; the two must stay in sync so
# that the strand written here and the strand the structure plots infer agree.
_CANON_RANK_COPIA = {"GAG": 0, "PROT": 1, "INT": 2, "RT": 3, "RH": 4}
_CANON_RANK_GYPSY = {"GAG": 0, "PROT": 1, "RT": 2, "RH": 3, "INT": 4}
_CANON_RANK_CORE = {"GAG": 0, "PROT": 1, "RT": 2, "RH": 3}


def warn(msg: str) -> None:
    print(f"[ltr_annotate] WARNING: {msg}", file=sys.stderr)


def log(msg: str) -> None:
    print(f"[ltr_annotate] {msg}")


# -----------------------------
# Depth-table discovery and I/O
# -----------------------------
class DepthTable(NamedTuple):
    path: str
    depth: int
    variant: str  # "clean" or "raw"


def discover_depth_tables(prefix: str, indir: str = ".",
                          variants: Sequence[str] = ("clean", "raw")
                          ) -> List[DepthTable]:
    """Every existing, non-empty depth table, numeric depth order, clean first."""
    found: List[DepthTable] = []
    for variant in variants:
        suffix = "_clean_ltr.tsv" if variant == "clean" else "_ltr.tsv"
        for depth in range(MAX_DEPTH + 1):
            path = os.path.join(indir, f"{prefix}_depth{depth}{suffix}")
            if os.path.isfile(path) and os.path.getsize(path) > 0:
                found.append(DepthTable(path, depth, variant))
    return found


def select_annotation_set(prefix: str, indir: str = "."
                          ) -> Tuple[str, List[DepthTable]]:
    """The set a consumer should read: FP-purged if any exists, else raw.

    Mirrors ltrharvest_plots.sh. A depth whose _clean_ table is absent while
    other _clean_ tables exist was fully purged by flag_fp_families.py and is
    correctly omitted -- never back-filled from the raw table.
    """
    for variant in ("clean", "raw"):
        tables = discover_depth_tables(prefix, indir, variants=(variant,))
        if tables:
            return variant, tables
    return "", []


def read_table(path: str) -> Tuple[Optional[List[str]], List[List[str]]]:
    """Return (header fields or None, data rows). The header keeps its '#'."""
    header: Optional[List[str]] = None
    rows: List[List[str]] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                if header is None:
                    header = line.split("\t")
                continue
            rows.append(line.split("\t"))
    return header, rows


def _target_mode(path: str) -> int:
    """Permissions the rewritten file should end up with.

    mkstemp deliberately creates at 0600 and os.replace carries that mode onto
    the target, so an atomic rewrite would otherwise strip group/world read from
    tables that were world-readable -- and leave them inconsistent with their
    untouched .fa siblings. Preserve what was there; fall back to the umask
    default for a file that does not exist yet.
    """
    try:
        return stat.S_IMODE(os.stat(path).st_mode)
    except OSError:
        umask = os.umask(0)
        os.umask(umask)
        return 0o666 & ~umask


def write_table(path: str, header: Optional[Sequence[str]],
                rows: Sequence[Sequence[str]]) -> None:
    """Rewrite `path` atomically, so an interrupted run cannot truncate it."""
    directory = os.path.dirname(os.path.abspath(path))
    mode = _target_mode(path)
    fd, tmp = tempfile.mkstemp(prefix=os.path.basename(path) + ".", suffix=".tmp",
                               dir=directory)
    try:
        with os.fdopen(fd, "w") as out:
            if header is not None:
                out.write("\t".join(header) + "\n")
            for row in rows:
                out.write("\t".join(row) + "\n")
        os.chmod(tmp, mode)
        os.replace(tmp, path)
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


def header_names(header: Optional[Sequence[str]]) -> List[str]:
    """Header fields with the leading '#' of the first column stripped."""
    if header is None:
        return []
    return [f.lstrip("#") for f in header]


def element_key(name: str) -> Optional[str]:
    """'chrom:s-e#Order/Super/Clade' -> 'chrom:s-e'. None if unparseable."""
    coord = name.split("#", 1)[0]
    return coord if ELEMENT_RE.match(coord) else None


def superfamily_of(name: str) -> str:
    """Middle field of the '#Order/Superfamily/Clade' suffix, or ''."""
    if "#" not in name:
        return ""
    bits = name.split("#", 1)[1].split("/")
    return bits[1] if len(bits) > 1 else ""


def parse_domains_field(field: str) -> List[Tuple[str, int, int]]:
    """Parse the `domains` column into [(gene, genomic_start, genomic_end), ...].

    Returns [] for '.', for empty, and for any field that is not a domains field
    at all -- callers may hand us a mis-indexed column and must not be poisoned
    by it.
    """
    if not field or field == UNKNOWN:
        return []
    out: List[Tuple[str, int, int]] = []
    for token in field.split(";"):
        m = DOMAIN_TOKEN_RE.match(token.strip())
        if not m:
            return []
        out.append((m.group(1), int(m.group(3)), int(m.group(4))))
    return out


# -----------------------------
# Strand: tier 1 (TEsorter2)
# -----------------------------
def load_tesorter_strands(prefix: str, indir: str = ".",
                          verbose: bool = False,
                          warn_missing: bool = True) -> Dict[str, str]:
    """element key -> '+', '-' or '.' from the TEsorter2 cls tables.

    The glob matches both the combined {base}.cls.tsv and the per-database
    {base}.rexdb.cls.tsv; both share the 7-column layout and column 6 is Strand.
    An explicit +/- from any round beats a '?' from another.
    """
    pattern = os.path.join(indir, f"{prefix}_r*.work", "*.cls.tsv")
    paths = sorted(glob.glob(pattern))
    if not paths:
        if warn_missing:
            warn(f"no TEsorter cls tables matched {pattern}; strand will be "
                 f"'{UNKNOWN}' unless recovered by domain order")
        return {}

    strands: Dict[str, str] = {}
    n_conflict = 0
    for path in paths:
        with open(path) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line or line.startswith("#"):
                    continue
                cols = line.split("\t")
                if len(cols) < 6 or not cols[0] or cols[0].lower() == "te":
                    continue
                key, value = cols[0], cols[5]
                if value in FLIP:
                    prev = strands.get(key)
                    if prev in FLIP and prev != value:
                        n_conflict += 1
                    strands[key] = value
                else:
                    strands.setdefault(key, UNKNOWN)
    if n_conflict:
        warn(f"{n_conflict} element(s) had conflicting strands across rounds; "
             f"kept the last one seen")
    if verbose:
        resolved = sum(1 for v in strands.values() if v in FLIP)
        log(f"tesorter: {resolved}/{len(strands)} stranded from "
            f"{len(paths)} cls table(s)")
    return strands


# -----------------------------
# Strand: tier 2 (domain order)
# -----------------------------
def infer_strand_from_domains(domains: Sequence[Tuple[str, int, int]],
                              superfamily: str) -> str:
    """Infer strand from the genomic order of domains against canonical order.

    Pairwise concordance vote over the ranked domain types present: a pair whose
    genomic order matches canonical order votes '+', otherwise '-'. Majority
    wins; a tie, or fewer than two ranked types, is unknown.

    Mirrors ltrharvest_plot_struct.py:infer_strand_from_domains, but returns '.'
    rather than '?' for unknown to match this table's placeholder convention.
    """
    sf = (superfamily or "").lower()
    if sf == "copia":
        rank = _CANON_RANK_COPIA
    elif sf == "gypsy":
        rank = _CANON_RANK_GYPSY
    else:
        rank = _CANON_RANK_CORE

    # 5'-most genomic start per ranked domain type.
    pos: Dict[str, int] = {}
    for gene, start, _end in domains:
        if gene not in rank:
            continue
        if gene not in pos or start < pos[gene]:
            pos[gene] = start
    if len(pos) < 2:
        return UNKNOWN

    genes = list(pos)
    plus = minus = 0
    for i in range(len(genes)):
        for j in range(i + 1, len(genes)):
            gi, gj = genes[i], genes[j]
            if rank[gi] == rank[gj] or pos[gi] == pos[gj]:
                continue
            if (rank[gi] < rank[gj]) == (pos[gi] < pos[gj]):
                plus += 1
            else:
                minus += 1
    if plus > minus:
        return "+"
    if minus > plus:
        return "-"
    return UNKNOWN


# -----------------------------
# Strand: tier 3 (pass-2 homology)
# -----------------------------
def load_pass2_links(prefix: str, indir: str = ".", verbose: bool = False,
                     warn_missing: bool = True
                     ) -> Tuple[Dict[str, str], Dict[Tuple[str, str], str]]:
    """Return (query -> pass-2 target, (query, target) -> alignment orientation).

    Only `pass` rows of pass2.tsv contribute a link. Orientation comes from the
    highest-scoring (largest matching-base count) PAF record for that pair.
    An empty orientation map means pass-2 transfer is unavailable, which is the
    normal state under --pass2-aligner blast.
    """
    target_of: Dict[str, str] = {}
    tsv_paths = sorted(glob.glob(
        os.path.join(indir, f"{prefix}_r*.work", "minimap2_pass2", "pass2.tsv")))
    for path in tsv_paths:
        with open(path) as fh:
            for line in fh:
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 6 or cols[1] != "pass":
                    continue
                target_of.setdefault(cols[0], cols[5])

    orientation: Dict[Tuple[str, str], str] = {}
    best_score: Dict[Tuple[str, str], int] = {}
    paf_paths = sorted(glob.glob(
        os.path.join(indir, f"{prefix}_r*.work", "minimap2_pass2", "pass2.paf")))
    for path in paf_paths:
        with open(path) as fh:
            for line in fh:
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 12 or cols[4] not in FLIP:
                    continue
                try:
                    score = int(cols[9])
                except ValueError:
                    continue
                pair = (cols[0], cols[5])
                if score > best_score.get(pair, -1):
                    best_score[pair] = score
                    orientation[pair] = cols[4]

    if target_of and not orientation and warn_missing:
        warn("pass-2 links found but no usable pass2.paf (expected with "
             "--pass2-aligner blast); skipping pass-2 strand transfer")
    if verbose:
        log(f"pass2: {len(target_of)} homology link(s) from {len(tsv_paths)} "
            f"table(s), {len(orientation)} oriented pair(s)")
    return target_of, orientation


# -----------------------------
# Strand cascade
# -----------------------------
class ElementInfo(NamedTuple):
    name: str          # full 'chrom:s-e#Order/Super/Clade' as written in column 1
    superfamily: str
    domains: List[Tuple[str, int, int]]


def collect_elements(loaded) -> Dict[str, ElementInfo]:
    """Element key -> ElementInfo, from already-stripped (header, rows) pairs.

    `loaded` is an iterable of (header, rows) as returned by load_unannotated.
    Taking the union across every table means the raw and clean variants of a
    depth can never disagree about an element.
    """
    elements: Dict[str, ElementInfo] = {}
    for header, rows in loaded:
        names = header_names(header)
        tail = _tail_width(names) if names else 2
        for row in rows:
            key = element_key(row[0]) if row else None
            if key is None or key in elements:
                continue
            # `domains` is the first of the trailing pair, wherever the row ends.
            domain_field = row[len(row) - tail] if len(row) >= tail else ""
            elements[key] = ElementInfo(row[0], superfamily_of(row[0]),
                                        parse_domains_field(domain_field))
    return elements


def resolve_strands(elements: Dict[str, ElementInfo],
                    tesorter: Dict[str, str],
                    target_of: Dict[str, str],
                    orientation: Dict[Tuple[str, str], str],
                    verbose: bool = False
                    ) -> Tuple[Dict[str, str], Dict[str, str]]:
    """Run the four-tier cascade.

    `elements` maps element key -> ElementInfo for the rows being annotated; it
    drives tier 2, which needs per-element evidence. Tiers 1 and 3 operate on the
    whole pool, so a pass-2 target that is not itself an output element can still
    donate its strand.

    Returns (key -> '+'/'-', key -> tier name).
    """
    strand: Dict[str, str] = {}
    source: Dict[str, str] = {}

    for key, value in tesorter.items():
        if value in FLIP:
            strand[key] = value
            source[key] = "tesorter"

    for key, info in elements.items():
        if key in strand:
            continue
        inferred = infer_strand_from_domains(info.domains, info.superfamily)
        if inferred in FLIP:
            strand[key] = inferred
            source[key] = "domain_order"

    for _ in range(MAX_PASS2_ITERS):
        added = 0
        for query, target in target_of.items():
            if query in strand:
                continue
            donor = strand.get(target)
            orient = orientation.get((query, target))
            if donor in FLIP and orient in FLIP:
                strand[query] = donor if orient == "+" else FLIP[donor]
                source[query] = "pass2"
                added += 1
        if not added:
            break

    if verbose:
        tiers = defaultdict(int)
        for key in elements:
            tiers[source.get(key, "unknown")] += 1
        log("strand tiers: " + ", ".join(
            f"{name}={tiers[name]}"
            for name in ("tesorter", "domain_order", "pass2", "unknown")))
    return strand, source


# -----------------------------
# Family
# -----------------------------
class Family(NamedTuple):
    family_id: str
    representative: str
    size: int


def load_families(prefix: str, indir: str = ".", verbose: bool = False,
                  warn_missing: bool = True
                  ) -> Tuple[Dict[str, Family], Dict[str, Family]]:
    """Return (by full element name, by bare coordinate key).

    Reads the Kmer2LTR/mmseqs consensus cluster table, whose two columns are
    representative and member; a family is every row sharing a representative.
    Families are numbered by descending size, ties broken on representative
    name, so labels are stable across reruns of the same data.
    """
    pattern = os.path.join(indir, f"{prefix}_all_ltr.consensus_id*_cluster.tsv")
    paths = sorted(glob.glob(pattern))
    if not paths:
        if warn_missing:
            warn(f"no consensus cluster table matched {pattern}; family will be "
                 f"'{UNKNOWN}'")
        return {}, {}
    if len(paths) > 1:
        warn(f"{len(paths)} consensus cluster tables matched {pattern}; "
             f"using {os.path.basename(paths[0])}")

    members_by_rep: Dict[str, List[str]] = defaultdict(list)
    seen: Dict[str, str] = {}
    with open(paths[0]) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 2 or not cols[0] or not cols[1]:
                continue
            rep, member = cols[0], cols[1]
            if member in seen:
                continue
            seen[member] = rep
            members_by_rep[rep].append(member)

    order = sorted(members_by_rep, key=lambda r: (-len(members_by_rep[r]), r))
    width = max(5, len(str(len(order))))

    by_name: Dict[str, Family] = {}
    by_coord: Dict[str, Family] = {}
    ambiguous_coords = set()
    for index, rep in enumerate(order, start=1):
        family = Family(f"{prefix}_fam{index:0{width}d}", rep,
                        len(members_by_rep[rep]))
        for member in members_by_rep[rep]:
            by_name[member] = family
            coord = element_key(member)
            if coord is None:
                continue
            if coord in by_coord and by_coord[coord].family_id != family.family_id:
                ambiguous_coords.add(coord)
            else:
                by_coord[coord] = family
    for coord in ambiguous_coords:
        by_coord.pop(coord, None)

    if verbose:
        log(f"family: {len(by_name)} member(s) in {len(order)} families from "
            f"{os.path.basename(paths[0])}")
    return by_name, by_coord


# -----------------------------
# Table rewriting
# -----------------------------
def _tail_width(names: Sequence[str]) -> int:
    """How many columns sit after `tsd` -- i.e. where to insert, counted from
    the end of each row. Falls back to 2 (domains, nest_status) when the anchor
    is missing, which keeps the trailing pair trailing whatever the layout."""
    if ANCHOR_COL in names:
        return len(names) - 1 - names.index(ANCHOR_COL)
    return 2


def strip_annotation_columns(header: Optional[List[str]], rows: List[List[str]]
                             ) -> Tuple[Optional[List[str]], List[List[str]]]:
    """Remove previously-written strand/family columns so a rerun replaces them.

    Indices are converted to offsets from the end of the row before deletion, so
    rows that are shorter than the header still lose the right fields.
    """
    names = header_names(header)
    targets = [names.index(c) for c in (STRAND_COL, FAMILY_COL) if c in names]
    if not header or not targets:
        return header, rows

    from_end = sorted((len(names) - i for i in targets), reverse=True)
    new_header = list(header)
    for offset in from_end:
        del new_header[len(new_header) - offset]

    new_rows = []
    for row in rows:
        trimmed = list(row)
        for offset in from_end:
            index = len(trimmed) - offset
            if 0 <= index < len(trimmed):
                del trimmed[index]
        new_rows.append(trimmed)
    return new_header, new_rows


def lookup_family(name: str, key: Optional[str],
                  family_by_name: Dict[str, Family],
                  family_by_coord: Dict[str, Family]) -> Optional[Family]:
    """Family for an element, matching on the full name first.

    The cluster table's members are full 'chrom:s-e#Order/Super/Clade' strings.
    The bare-coordinate fallback covers the case where a classification string
    drifted between clustering and the final table.
    """
    found = family_by_name.get(name)
    if found is None and key is not None:
        found = family_by_coord.get(key)
    return found


def load_unannotated(path: str) -> Tuple[Optional[List[str]], List[List[str]]]:
    """Read a depth table with any previously-written strand/family removed.

    Every consumer of a table's field offsets must go through this, so that a
    rerun sees the same layout as a first run. Reading the raw file instead
    would put `tsd`'s tail width at 4 and make `domains` resolve to the `strand`
    column.
    """
    return strip_annotation_columns(*read_table(path))


def write_annotated_table(table: DepthTable,
                          header: Optional[List[str]],
                          rows: List[List[str]],
                          strand: Dict[str, str],
                          family_by_name: Dict[str, Family],
                          family_by_coord: Dict[str, Family],
                          verbose: bool = False) -> Tuple[int, int]:
    """Insert strand and family into already-stripped rows and rewrite the file.

    Returns (rows written, rows whose element name was unparseable).
    """
    names = header_names(header)
    tail = _tail_width(names) if names else 2

    if header is not None:
        insert_at = len(header) - tail
        header = header[:insert_at] + [STRAND_COL, FAMILY_COL] + header[insert_at:]

    skipped = 0
    for row in rows:
        name = row[0] if row else ""
        key = element_key(name)
        if key is None:
            skipped += 1
        family = lookup_family(name, key, family_by_name, family_by_coord)
        values = [
            strand.get(key, UNKNOWN) if key else UNKNOWN,
            family.family_id if family else UNKNOWN,
        ]
        insert_at = max(0, len(row) - tail)
        row[insert_at:insert_at] = values

    write_table(table.path, header, rows)
    if verbose:
        log(f"  {os.path.basename(table.path)}: {len(rows)} row(s)"
            + (f", {skipped} unparseable name(s)" if skipped else ""))
    return len(rows), skipped


# -----------------------------
# Driver
# -----------------------------
def annotate(prefix: str, indir: str = ".", verbose: bool = False) -> int:
    tables = discover_depth_tables(prefix, indir)
    if not tables:
        print(f"[ltr_annotate] ERROR: no {prefix}_depth<N>[_clean]_ltr.tsv found "
              f"in {indir}", file=sys.stderr)
        return 1

    # Read every table once, stripped of any prior annotation, and keep the rows
    # for the write pass below.
    loaded = [(table,) + load_unannotated(table.path) for table in tables]
    elements = collect_elements((header, rows) for _t, header, rows in loaded)

    tesorter = load_tesorter_strands(prefix, indir, verbose)
    target_of, orientation = load_pass2_links(prefix, indir, verbose)
    strand, source = resolve_strands(elements, tesorter, target_of, orientation,
                                     verbose)
    family_by_name, family_by_coord = load_families(prefix, indir, verbose)

    total_rows = 0
    total_skipped = 0
    for table, header, rows in loaded:
        written, skipped = write_annotated_table(
            table, header, rows, strand, family_by_name, family_by_coord, verbose)
        total_rows += written
        total_skipped += skipped

    tiers = defaultdict(int)
    for key in elements:
        tiers[source.get(key, "unknown")] += 1
    log(f"strand: {tiers['tesorter']} tesorter, {tiers['domain_order']} "
        f"domain_order, {tiers['pass2']} pass2, {tiers['unknown']} unknown "
        f"({len(elements)} elements)")

    hits = [lookup_family(info.name, key, family_by_name, family_by_coord)
            for key, info in elements.items()]
    assigned = [f for f in hits if f is not None]
    n_families = len({f.family_id for f in assigned})
    log(f"family: {len(assigned)}/{len(elements)} assigned across "
        f"{n_families} families")

    if total_skipped:
        warn(f"{total_skipped} row(s) had an unparseable element name and got "
             f"'{UNKNOWN}' for both columns")
    log(f"annotated {len(tables)} table(s): "
        + ", ".join(os.path.basename(t.path) for t in tables))
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Add strand and family columns to the depth-bucketed "
                    "LTR-RT tables produced by ltrharvest_wrapper2.sh.")
    parser.add_argument("--prefix", required=True,
                        help="Wrapper output prefix, e.g. Athal_tair10_chr2_LTRs")
    parser.add_argument("--indir", default=".",
                        help="Directory holding the wrapper outputs (default: .)")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Per-step progress and per-file counts")
    args = parser.parse_args(argv)

    if not os.path.isdir(args.indir):
        print(f"[ltr_annotate] ERROR: not a directory: {args.indir}",
              file=sys.stderr)
        return 1
    return annotate(args.prefix, args.indir, args.verbose)


if __name__ == "__main__":
    sys.exit(main())
