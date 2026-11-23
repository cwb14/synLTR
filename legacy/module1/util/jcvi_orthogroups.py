#!/usr/bin/env python3
import argparse
import glob
import os
import re
from collections import defaultdict, Counter

gene_re = re.compile(r'^([A-Za-z]+)(\d+)(?:_(\d+))?$')
# Parses "Bdact016340" or "Etef739571_2" -> ("Bdact","016340","2"/None)

chrom_sp_re = re.compile(r'^([A-Za-z]+)_(chr|sca)')  # e.g., Bdact_chr2, Etef_sca12

def parse_gene(tok):
    m = gene_re.match(tok)
    if not m:
        raise ValueError(f"Unrecognized gene token: {tok}")
    sp, base, dup = m.group(1), m.group(2), m.group(3)
    return sp, base, dup  # species, base numeric string, copy suffix or None

def read_anchors():
    """
    Returns:
      species: sorted list of species appearing in filenames
      edges: dict keyed by (A,B) -> dict A_base -> set(B_base)
      copies_seen: copies_seen[sp][base] = set of full ids (with suffixes) observed
      pair_fulls: pair_fulls[(A,B)][(A_base,B_base)] -> set of B FULL ids that pair with ANY copy of A_base
                  (used to output duplicate copy IDs per species)
    """
    edges = dict()
    pair_fulls = dict()
    copies_seen = defaultdict(lambda: defaultdict(set))
    species_set = set()

    for path in glob.glob("*.anchors"):
        base = os.path.basename(path)
        parts = base.split(".")
        if len(parts) < 3 or parts[-1] != "anchors":
            continue
        A, B = parts[0], parts[1]
        species_set.update([A, B])

        keyAB = (A, B)
        keyBA = (B, A)
        if keyAB not in edges:
            edges[keyAB] = defaultdict(set)
        if keyBA not in edges:
            edges[keyBA] = defaultdict(set)
        if keyAB not in pair_fulls:
            pair_fulls[keyAB] = defaultdict(set)
        if keyBA not in pair_fulls:
            pair_fulls[keyBA] = defaultdict(set)

        with open(path) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line or line.startswith("#"):
                    continue
                cols = line.split("\t")
                if len(cols) < 2:
                    continue
                gA_full, gB_full = cols[0], cols[1]

                spA, baseA, dupA = parse_gene(gA_full)
                spB, baseB, dupB = parse_gene(gB_full)

                copies_seen[spA][baseA].add(gA_full)
                copies_seen[spB][baseB].add(gB_full)

                edges[keyAB][baseA].add(baseB)
                edges[keyBA][baseB].add(baseA)

                pair_fulls[keyAB][(baseA, baseB)].add(gB_full)
                pair_fulls[keyBA][(baseB, baseA)].add(gA_full)

    species = sorted(species_set)
    return species, edges, copies_seen, pair_fulls

def choose_reference_species(species, prefer=None):
    if prefer is not None:
        if prefer not in species:
            raise SystemExit(f"Requested reference species '{prefer}' not found among: {', '.join(species)}")
        return prefer
    return species[0]

def find_orthogroups(species, edges, copies_seen, pair_fulls, ref_sp):
    """
    Finds fully reciprocal, across-all-species orthogroups at the *base-gene* level.
    OUTPUT: the reference species outputs FULL IDs (with suffixes) for all reference duplicates
            that pair with at least one chosen base in the other species. These ref copies are first.
    """
    others = [s for s in species if s != ref_sp]

    needed_pairs = set()
    for i in range(len(species)):
        for j in range(i+1, len(species)):
            needed_pairs.add((species[i], species[j]))
            needed_pairs.add((species[j], species[i]))
    for a,b in needed_pairs:
        if (a,b) not in edges:
            raise SystemExit(f"Missing anchors for pair {a}.{b}.anchors (or equivalent); cannot enforce reciprocity.")

    ref_bases = set(copies_seen[ref_sp].keys())
    emitted_keys = set()

    for ref_base in sorted(ref_bases):
        candidates = dict()
        for sp in others:
            neigh = edges[(ref_sp, sp)].get(ref_base, set())
            candidates[sp] = set(neigh) if neigh else set()

        if any(len(candidates[sp]) == 0 for sp in others):
            continue

        order = sorted(others, key=lambda s: len(candidates[s]))
        selection = dict()
        solutions = []

        def ok_pair(si, bi, sj, bj):
            return (bj in edges[(si, sj)].get(bi, set())) and (bi in edges[(sj, si)].get(bj, set()))

        def dfs(idx):
            if idx == len(order):
                solutions.append(selection.copy())
                return
            sp = order[idx]
            for base_choice in sorted(candidates[sp]):
                good = True
                for sp2, base2 in selection.items():
                    if not ok_pair(sp, base_choice, sp2, base2):
                        good = False
                        break
                if not good:
                    continue
                selection[sp] = base_choice
                dfs(idx+1)
                del selection[sp]

        dfs(0)
        if not solutions:
            continue

        for sol in solutions:
            key_tuple = (ref_sp, ref_base) + tuple((sp, sol[sp]) for sp in sorted(sol.keys()))
            if key_tuple in emitted_keys:
                continue
            emitted_keys.add(key_tuple)

            out_tokens = []

            # Collect all reference FULL copies that actually pair with the chosen bases
            ref_fulls_union = set()
            for sp in sorted(others):
                chosen_base = sol[sp]
                # from (sp -> ref_sp): which REF FULL ids pair with chosen_base?
                ref_fulls_union |= pair_fulls[(sp, ref_sp)].get((chosen_base, ref_base), set())

            if not ref_fulls_union:
                ref_fulls_union = set(copies_seen[ref_sp][ref_base])

            for f in sorted(ref_fulls_union):
                out_tokens.append(f)

            # Then, for each other species, append their FULL copies that pair with the ref_base
            for sp in sorted(others):
                chosen_base = sol[sp]
                fulls = pair_fulls[(ref_sp, sp)].get((ref_base, chosen_base), set())
                if not fulls:
                    fulls = copies_seen[sp][chosen_base]
                for f in sorted(fulls):
                    out_tokens.append(f)

            yield out_tokens

def load_filter_counts(path):
    """
    Reads a two-column file specifying per-species copy-count constraints.
    Format per line:
      <Species><whitespace><CountSpec>
    Where CountSpec can be:
      N   (exactly N)
      N-  (≤ N)
      N+  (≥ N)

    Lines starting with '#' or blank lines are ignored.
    Unspecified species are unconstrained.

    Returns:
      dict: species -> (op, value) where op ∈ {'==', '<=', '>='}
    """
    filt = {}
    pat = re.compile(r'^(\d+)([+-])?$')  # captures integer and optional +/- suffix
    with open(path) as fh:
        for line in fh:
            raw = line.strip()
            if not raw or raw.startswith("#"):
                continue
            parts = raw.split()
            if len(parts) < 2:
                raise SystemExit(f"Malformed filter line (need 2 columns): {line.rstrip()}")
            sp = parts[0].strip()
            m = pat.match(parts[1])
            if not m:
                raise SystemExit(f"Invalid count spec for species '{sp}': {parts[1]} (use N, N+, or N-)")
            val = int(m.group(1))
            suf = m.group(2)
            if suf == '+':
                op = '>='
            elif suf == '-':
                op = '<='
            else:
                op = '=='
            filt[sp] = (op, val)
    if not filt:
        raise SystemExit("Filter file parsed but contained no rules.")
    return filt

def passes_filter(tokens, filt_counts):
    """
    tokens: list of full IDs (e.g., 'Bdact001522', 'Cdact003658_1', ...)
    filt_counts: dict species -> (op, value), where op ∈ {'==', '<=', '>='}
    Returns True if each specified species satisfies its constraint.
    """
    sp_counts = Counter()
    for t in tokens:
        sp, _, _ = parse_gene(t)
        sp_counts[sp] += 1

    for sp, (op, val) in filt_counts.items():
        c = sp_counts.get(sp, 0)
        if op == '==' and c != val:
            return False
        if op == '<=' and c > val:
            return False
        if op == '>=' and c < val:
            return False
    return True

# ----------------------
# Syntenic-block filter
# ----------------------

def _extract_species_from_chrom(chrom):
    """
    chrom examples: 'Bdact_chr2', 'Cdact_chr2A', 'Etef_sca12'
    returns 'Bdact', 'Cdact', 'Etef' respectively (based on prefix before _chr/_sca)
    """
    m = chrom_sp_re.match(chrom)
    if not m:
        # fallback: take alpha prefix
        return re.match(r'^([A-Za-z]+)', chrom).group(1)
    return m.group(1)

def load_beds(bed_patterns):
    """
    bed_patterns: list of patterns or file paths; we will glob-expand each
    Returns: dict full_gene_id -> (species, chrom, start, end)
    """
    gene2pos = {}
    files = []
    for p in bed_patterns:
        files.extend(glob.glob(p))
    for fp in files:
        with open(fp) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):  # tolerant
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 4:
                    continue
                chrom, s, e, name = cols[0], cols[1], cols[2], cols[3]
                try:
                    start = int(s)
                    end = int(e)
                except ValueError:
                    continue
                sp = _extract_species_from_chrom(chrom)
                gene2pos[name] = (sp, chrom, start, end)
    return gene2pos

def parse_region(tok):
    """
    tok like 'Bdact_chr1:13911..75338' -> (species, chrom, 13911, 75338)
    """
    chrom, coords = tok.split(":")
    s, e = coords.split("..")
    start, end = int(s), int(e)
    sp = _extract_species_from_chrom(chrom)
    return sp, chrom, start, end

def load_blocks(block_patterns):
    """
    block_patterns: list of patterns or file paths (glob-expanded)
    Returns:
      blocks[(spA, spB)] = list of (A_chrom, A_s, A_e, B_chrom, B_s, B_e)
    We store both directions, so lookups are easy either way.
    """
    blocks = defaultdict(list)
    files = []
    for p in block_patterns:
        files.extend(glob.glob(p))
    for fp in files:
        with open(fp) as fh:
            for line in fh:
                raw = line.strip()
                if not raw or raw.startswith("#"):
                    continue
                parts = raw.split()
                if len(parts) < 2:
                    continue
                # accept tab or space separation; first two fields are the regions
                a_tok = parts[0]
                b_tok = parts[1]
                try:
                    spA, A_chrom, A_s, A_e = parse_region(a_tok)
                    spB, B_chrom, B_s, B_e = parse_region(b_tok)
                except Exception:
                    continue
                # store both directions
                blocks[(spA, spB)].append((A_chrom, A_s, A_e, B_chrom, B_s, B_e))
                blocks[(spB, spA)].append((B_chrom, B_s, B_e, A_chrom, A_s, A_e))
    return blocks

def _gene_within(chrom, s, e, block_chrom, bs, be):
    """Return True if chrom matches and [s,e] is fully inside [bs,be]."""
    if chrom != block_chrom:
        return False
    return s >= bs and e <= be

def _pair_supported_by_block(geneA_pos, geneB_pos, block_list):
    """
    geneA_pos: (chromA, sA, eA); geneB_pos: (chromB, sB, eB)
    block_list: list of (A_chrom, A_s, A_e, B_chrom, B_s, B_e) for a specific (spA, spB)
    Return True if any block contains BOTH geneA and geneB.
    """
    chromA, sA, eA = geneA_pos
    chromB, sB, eB = geneB_pos
    for A_chrom, A_s, A_e, B_chrom, B_s, B_e in block_list:
        if _gene_within(chromA, sA, eA, A_chrom, A_s, A_e) and _gene_within(chromB, sB, eB, B_chrom, B_s, B_e):
            return True
    return False

def passes_synteny(tokens, ref_sp, gene2pos, blocks):
    """
    tokens: list of full IDs for one orthogroup line
    ref_sp: reference species (first in chain)
    gene2pos: dict full_id -> (species, chrom, start, end)
    blocks: blocks[(spA, spB)] -> list of block tuples
    Strategy:
      - Build species -> list of full IDs present in this orthogroup.
      - Form a "chain" of species: [ref, *sorted(others)].
      - For each adjacent pair in the chain (si, sj):
            require that there exists at least one (gene_i, gene_j) pair
            that is covered by a syntenic block in blocks[(si, sj)].
        If no block data exists for (si, sj), treat this pair as "unknown" and pass it.
      - If any required pair fails, the orthogroup fails.
    """
    # group genes by species
    sp2genes = defaultdict(list)
    for g in tokens:
        sp, _, _ = parse_gene(g)
        sp2genes[sp].append(g)

    # ensure we have positions for every gene; if not, fail
    for g in tokens:
        if g not in gene2pos:
            return False

    species_in_group = list(sp2genes.keys())
    others = sorted([s for s in species_in_group if s != ref_sp])
    chain = [ref_sp] + others

    for i in range(len(chain) - 1):
        si, sj = chain[i], chain[i+1]
        block_key = (si, sj)
        if block_key not in blocks or len(blocks[block_key]) == 0:
            # No info to verify this adjacency; allow it.
            continue
        ok = False
        # try all combinations; success if any pair shares a block
        for gi in sp2genes[si]:
            _, chrom_i, s_i, e_i = gene2pos[gi]
            for gj in sp2genes[sj]:
                _, chrom_j, s_j, e_j = gene2pos[gj]
                if _pair_supported_by_block(
                    (chrom_i, s_i, e_i),
                    (chrom_j, s_j, e_j),
                    blocks[block_key]
                ):
                    ok = True
                    break
            if ok:
                break
        if not ok:
            return False
    return True

def main():
    ap = argparse.ArgumentParser(
        description="Build fully reciprocal, across-all-species orthogroups from *.anchors files, "
                    "optionally filtering by (A) exact per-species copy counts OR (B) syntenic block matching."
    )
    ap.add_argument("--ref", help="Reference species to key output lines by (default: lexicographically first species).")
    ap.add_argument("-o", "--out", help="Output TSV file (default: stdout).")

    # Filter A: per-species exact copy counts
    ap.add_argument("-f", "--filter",
                    help="Two-column file: <Species> <CountSpec>. CountSpec: N (exact), N- (≤N), N+ (≥N). "
                         "Only keep orthogroups whose per-species copy counts satisfy the provided constraints. "
                         "Unspecified species are unconstrained.")

    # Filter B: syntenic block matching
    ap.add_argument("--beds", nargs="+",
                    help="BED files (glob ok). Each row: chrom  start  end  geneID  ... "
                         "Chrom must be like '[species]_chr#' or '[species]_sca#'.")
    ap.add_argument("--blocks", nargs="+",
                    help="Syntenic block files (glob ok). Lines like: "
                         "'Bdact_chr1:13911..75338\\tCdact_chr1A:411902..478939\\t+'. "
                         "We only use the first two columns.")

    args = ap.parse_args()

    # Read anchors & build orthogroups
    species, edges, copies_seen, pair_fulls = read_anchors()
    if not species:
        raise SystemExit("No *.anchors files found.")

    ref_sp = choose_reference_species(species, args.ref)
    lines = list(find_orthogroups(species, edges, copies_seen, pair_fulls, ref_sp))

    # Decide which filters (if any) are active
    use_count_filter = bool(args.filter)
    use_block_filter = bool(args.beds and args.blocks)

    # Apply per-species copy-count filter first (if requested)
    if use_count_filter:
        filt_counts = load_filter_counts(args.filter)
        lines = [toks for toks in lines if passes_filter(toks, filt_counts)]

    # Then apply syntenic-block filter (if requested)
    if use_block_filter:
        gene2pos = load_beds(args.beds)
        blocks = load_blocks(args.blocks)
        filtered = []
        for toks in lines:
            if passes_synteny(toks, ref_sp, gene2pos, blocks):
                filtered.append(toks)
        lines = filtered

    # Optional syntenic-block filter
    elif use_block_filter:
        gene2pos = load_beds(args.beds)
        blocks = load_blocks(args.blocks)
        filtered = []
        for toks in lines:
            if passes_synteny(toks, ref_sp, gene2pos, blocks):
                filtered.append(toks)
        lines = filtered

    # else: no filtering

    if args.out:
        with open(args.out, "w") as out:
            for toks in lines:
                out.write("\t".join(toks) + "\n")
    else:
        for toks in lines:
            print("\t".join(toks))

if __name__ == "__main__":
    main()
