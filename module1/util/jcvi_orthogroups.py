#!/usr/bin/env python3
import argparse
import glob
import os
import re
from collections import defaultdict, Counter

gene_re = re.compile(r'^([A-Za-z]+)(\d+)(?:_(\d+))?$')
# Parses "Bdact016340" or "Etef739571_2" -> ("Bdact","016340","2"/None)

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
    OUTPUT CHANGE (already present): the reference species now outputs FULL IDs (with suffixes)
    for all reference duplicates that actually pair with at least one chosen
    base in the other species. These reference copies are placed first on the line.
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
    Reads a two-column file specifying exact copy counts per species.
    Format: <Species><whitespace><Count>
    Lines starting with '#' or blank lines are ignored.
    Unspecified species are unconstrained.
    """
    filt = {}
    with open(path) as fh:
        for line in fh:
            raw = line.strip()
            if not raw or raw.startswith("#"):
                continue
            parts = raw.split()
            if len(parts) < 2:
                raise SystemExit(f"Malformed filter line (need 2 columns): {line.rstrip()}")
            sp = parts[0].strip()
            try:
                cnt = int(parts[1])
            except ValueError:
                raise SystemExit(f"Non-integer count in filter for species '{sp}': {parts[1]}")
            if cnt < 0:
                raise SystemExit(f"Negative count in filter for species '{sp}': {cnt}")
            filt[sp] = cnt
    if not filt:
        raise SystemExit("Filter file parsed but contained no rules.")
    return filt

def passes_filter(tokens, filt_counts):
    """
    tokens: list of full IDs (e.g., 'Bdact001522', 'Cdact003658_1', ...)
    filt_counts: dict species -> exact count required
    Returns True if for each specified species, the number of tokens with that prefix equals the required count.
    """
    # Count species occurrences in this orthogroup line
    sp_counts = Counter()
    for t in tokens:
        sp, _, _ = parse_gene(t)
        sp_counts[sp] += 1

    # Check exact matches for only the specified species
    for sp, required in filt_counts.items():
        if sp_counts.get(sp, 0) != required:
            return False
    return True

def main():
    ap = argparse.ArgumentParser(
        description="Build fully reciprocal, across-all-species orthogroups from *.anchors files, "
                    "optionally filtering by exact per-species copy counts."
    )
    ap.add_argument("--ref", help="Reference species to key output lines by (default: lexicographically first species).")
    ap.add_argument("-o", "--out", help="Output TSV file (default: stdout).")
    ap.add_argument("-f", "--filter", help="Two-column file: <Species> <Count>. Only keep orthogroups whose per-species "
                                           "copy counts EXACTLY match the provided counts. Unspecified species are unconstrained.")
    args = ap.parse_args()

    species, edges, copies_seen, pair_fulls = read_anchors()
    if not species:
        raise SystemExit("No *.anchors files found.")

    ref_sp = choose_reference_species(species, args.ref)
    lines = list(find_orthogroups(species, edges, copies_seen, pair_fulls, ref_sp))

    # Optional filtering step
    if args.filter:
        filt_counts = load_filter_counts(args.filter)
        lines = [toks for toks in lines if passes_filter(toks, filt_counts)]

    if args.out:
        with open(args.out, "w") as out:
            for toks in lines:
                out.write("\t".join(toks) + "\n")
    else:
        for toks in lines:
            print("\t".join(toks))

if __name__ == "__main__":
    main()
