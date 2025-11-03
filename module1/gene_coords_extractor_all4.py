import argparse
import os
import sys

def read_bed_data(bed_files):
    merged_data = []
    for bed_file in bed_files:
        try_paths = [bed_file]
        if not os.path.exists(bed_file):
            raise FileNotFoundError(f"BED not found: {bed_file}")
        with open(bed_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith(('#', 'track', 'browser')):
                    continue
                merged_data.append(line)
    return merged_data

def get_gene_coords(merged_data):
    gene_coords = {}
    for line in merged_data:
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        chrom, start, end, geneID = parts[0], parts[1], parts[2], parts[3]
        try:
            gene_coords[geneID] = (chrom, int(start), int(end))
        except ValueError:
            continue
    return gene_coords

def read_gene_ids_file(gene_ids_filename):
    clusters = []
    with open(gene_ids_filename, 'r') as f:
        cluster = []
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s == "###":
                if cluster:
                    clusters.append(cluster)
                    cluster = []
            else:
                fields = s.split()
                if len(fields) >= 2:
                    cluster.append(fields[:2])
    if cluster:
        clusters.append(cluster)
    return clusters

def process_clusters(clusters, gene_coords):
    output_strings = []
    for cluster in clusters:
        for i in range(len(cluster) - 1):
            g1a, g1b = cluster[i][0], cluster[i][1]
            g2a, g2b = cluster[i + 1][0], cluster[i + 1][1]

            missing = [g for g in (g1a, g1b, g2a, g2b) if g not in gene_coords]
            if missing:
                print(f"WARNING: missing {len(missing)} gene(s) in BED: {', '.join(missing)}", file=sys.stderr)
                continue

            coords_1 = gene_coords[g1a]
            coords_2 = gene_coords[g1b]
            coords_3 = gene_coords[g2a]
            coords_4 = gene_coords[g2b]

            direction_1 = "+" if coords_1[1] < coords_3[1] else "-"
            direction_2 = "+" if coords_2[1] < coords_4[1] else "-"
            directionality = "+" if direction_1 == direction_2 else "-"

            left_span = f"{coords_1[0]}:{min(coords_1[1], coords_3[1])}..{max(coords_1[2], coords_3[2])}"
            right_span = f"{coords_2[0]}:{min(coords_2[1], coords_4[1])}..{max(coords_2[2], coords_4[2])}"

            output_strings.append(f"{left_span}\t{right_span}\t{directionality}")
    return output_strings

def parse_species_from_mcscan(mcscan_filename):
    """
    Expect basename like 'SpeciesA.SpeciesB.anchors' or 'SpeciesA.SpeciesB.anchors.clean'
    Return ('SpeciesA', 'SpeciesB').
    """
    base = os.path.basename(mcscan_filename)
    toks = base.split('.')
    if len(toks) < 3:
        raise ValueError(f"Cannot parse species names from anchors file name: {base}")
    # First two tokens are species
    return toks[0], toks[1]

def resolve_bed_path(species, mcscan_path, bed_dir=None):
    """
    Try a series of plausible locations for <species>.bed.
    """
    tried = []
    candidates = []

    if bed_dir:
        candidates.append(os.path.join(bed_dir, f"{species}.bed"))

    anchors_dir = os.path.dirname(os.path.abspath(mcscan_path))
    if anchors_dir:
        candidates.append(os.path.join(anchors_dir, f"{species}.bed"))

    # common convention in this repo layout
    candidates.append(os.path.join(os.getcwd(), "results", f"{species}.bed"))
    # fallback to CWD
    candidates.append(os.path.join(os.getcwd(), f"{species}.bed"))

    for p in candidates:
        tried.append(p)
        if os.path.exists(p):
            return p
    raise FileNotFoundError(
        f"Could not locate BED for '{species}'. Tried:\n  - " + "\n  - ".join(tried)
    )

def derive_bed_files(mcscan_filename, bed_dir=None, bed1_override=None, bed2_override=None):
    if bed1_override and bed2_override:
        return [bed1_override, bed2_override]

    sp1, sp2 = parse_species_from_mcscan(mcscan_filename)
    bed1 = bed1_override or resolve_bed_path(sp1, mcscan_filename, bed_dir=bed_dir)
    bed2 = bed2_override or resolve_bed_path(sp2, mcscan_filename, bed_dir=bed_dir)
    return [bed1, bed2]

def main():
    parser = argparse.ArgumentParser(description='Extract gene coordinates and determine directionality.')
    parser.add_argument('-mcscan', required=True, help='Input anchors file (e.g., SpeciesA.SpeciesB.anchors[.clean])')
    parser.add_argument('--bed-dir', default=None, help='Directory containing <species>.bed files (optional).')
    parser.add_argument('--bed1', default=None, help='Explicit path to the first BED (overrides inference).')
    parser.add_argument('--bed2', default=None, help='Explicit path to the second BED (overrides inference).')
    args = parser.parse_args()

    try:
        bed_files = derive_bed_files(args.mcscan, bed_dir=args.bed_dir, bed1_override=args.bed1, bed2_override=args.bed2)
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(2)

    print(f"# Using BEDs:\n#   1) {bed_files[0]}\n#   2) {bed_files[1]}", file=sys.stderr)

    merged_data = read_bed_data(bed_files)
    gene_coords = get_gene_coords(merged_data)
    clusters = read_gene_ids_file(args.mcscan)
    output_strings = process_clusters(clusters, gene_coords)

    for output_str in output_strings:
        print(output_str)

if __name__ == "__main__":
    main()
