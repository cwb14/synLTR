import argparse
import os
import sys

def read_bed_data(bed_files):
    merged_data = []
    for bed_file in bed_files:
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
            # skip malformed lines
            continue
        chrom, start, end, geneID = parts[0], parts[1], parts[2], parts[3]
        try:
            gene_coords[geneID] = (chrom, int(start), int(end))
        except ValueError:
            # skip lines with non-integer coords
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
                # anchor lines: geneA  geneB  score
                # we keep the first two tokens (gene IDs); any trailing fields ignored
                fields = s.split()
                if len(fields) >= 2:
                    cluster.append(fields[:2])
    if cluster:  # flush last cluster if file doesn't end with ###
        clusters.append(cluster)
    return clusters

def process_clusters(clusters, gene_coords):
    output_strings = []
    for cluster in clusters:
        # iterate over adjacent anchor pairs within the cluster
        for i in range(len(cluster) - 1):
            g1a, g1b = cluster[i][0], cluster[i][1]
            g2a, g2b = cluster[i + 1][0], cluster[i + 1][1]

            # verify presence in gene_coords
            missing = [g for g in (g1a, g1b, g2a, g2b) if g not in gene_coords]
            if missing:
                print(f"WARNING: missing {len(missing)} gene(s) in BED: {', '.join(missing)}", file=sys.stderr)
                # skip this pair if any gene is missing
                continue

            coords_1 = gene_coords[g1a]
            coords_2 = gene_coords[g1b]
            coords_3 = gene_coords[g2a]
            coords_4 = gene_coords[g2b]

            # direction by start coordinate ordering
            direction_1 = "+" if coords_1[1] < coords_3[1] else "-"
            direction_2 = "+" if coords_2[1] < coords_4[1] else "-"
            directionality = "+" if direction_1 == direction_2 else "-"

            # span each syntenic side from min(start) to max(end)
            left_span = f"{coords_1[0]}:{min(coords_1[1], coords_3[1])}..{max(coords_1[2], coords_3[2])}"
            right_span = f"{coords_2[0]}:{min(coords_2[1], coords_4[1])}..{max(coords_2[2], coords_4[2])}"

            output_strings.append(f"{left_span}\t{right_span}\t{directionality}")
    return output_strings

def derive_bed_files(mcscan_filename):
    name_parts = mcscan_filename.split('.')[0:2]
    return [f"{name_parts[0]}.bed", f"{name_parts[1]}.bed"]

def main():
    parser = argparse.ArgumentParser(description='Extract gene coordinates and determine directionality.')
    parser.add_argument('-mcscan', required=True, help='Input anchors file (e.g., SpeciesA.SpeciesB.anchors)')
    args = parser.parse_args()

    bed_files = derive_bed_files(args.mcscan)
    merged_data = read_bed_data(bed_files)
    gene_coords = get_gene_coords(merged_data)
    clusters = read_gene_ids_file(args.mcscan)
    output_strings = process_clusters(clusters, gene_coords)

    for output_str in output_strings:
        print(output_str)

if __name__ == "__main__":
    main()
