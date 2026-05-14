#!/bin/bash
#SBATCH --job-name="extractinterbaits"
#SBATCH --export=ALL
#SBATCH --mem=200G  
#SBATCH --partition=medium



conda activate seqkit_env



BED="intersection_GERP_HyPhy.bed"
FASTA="/home/tmichel/projects/rbge/tmichel/reference_genomes/Hannah_Begonia_baits_edited.fasta"
OUTPUT="intersected_baits.fasta"

# Extract bait names from BED (first column)
awk '{print $1}' "$BED" | sort -u > bait_names.txt

# Extract matching sequences, ignoring the "Hannah_Begonia_baits-" prefix
awk '
BEGIN {
    while ((getline < "bait_names.txt") > 0) names[$1]=1
}
/^>/ {
    # Remove ">" and prefix "Hannah_Begonia_baits-"
    hdr = substr($0,2)
    sub(/^Hannah_Begonia_baits-/, "", hdr)
    split(hdr, id, " ")   # take first word
    keep = (id[1] in names)
}
keep { print }
' "$FASTA" > "$OUTPUT"

echo "Sequences extracted to $OUTPUT"

