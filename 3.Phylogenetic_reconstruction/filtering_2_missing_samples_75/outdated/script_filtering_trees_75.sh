#!/bin/bash
#SBATCH --job-name=check_gene_trees
#SBATCH --mem=8G
#SBATCH --partition=short

#!/usr/bin/env bash

TREE_DIR="/home/tmichel/scratch/Deleterious_alleles_PNG/CAPTUS/CAPTUS-PNG/05_phylogeny_FastTreeAstral/genetree_from_filtered_MSA"



# Get list of .tree files with full path
TREE_FILES=("$TREE_DIR"/*.tree)

if [ ${#TREE_FILES[@]} -eq 0 ]; then
    echo "No .tree files found in: $TREE_DIR"
    exit 1
fi

# Extract all taxa (specimen names) from all trees
ALL_SPS=$(grep -ho "[^():,;]*" "${TREE_FILES[@]}" \
    | awk 'NF' \
    | sort -u)

SPS_COUNT=$(echo "$ALL_SPS" | wc -l)
GENE_COUNT=${#TREE_FILES[@]}

# 75% threshold
THRESHOLD=$(printf "%.0f" "$(echo "$SPS_COUNT * 0.75" | bc -l)")

# Output files (written to current directory!)
GOOD_LIST="good_genes.txt"
BAD_LIST="bad_genes.txt"
SUMMARY="gene_tree_summary.txt"

> "$GOOD_LIST"
> "$BAD_LIST"
> "$SUMMARY"

# Evaluate each tree
for TREE in "${TREE_FILES[@]}"; do
    TAXA_IN_TREE=$(grep -ho "[^():,;]*" "$TREE" \
        | awk 'NF' \
        | sort -u \
        | wc -l)

    if (( TAXA_IN_TREE >= THRESHOLD )); then
        echo "$TREE" >> "$GOOD_LIST"
    else
        echo "$TREE" >> "$BAD_LIST"
    fi
done

GOOD_COUNT=$(wc -l < "$GOOD_LIST")
BAD_COUNT=$(wc -l < "$BAD_LIST")

# Write summary file
cat <<EOF > "$SUMMARY"
Total unique specimens (Sps): $SPS_COUNT
Total number of genes (trees): $GENE_COUNT
75% of specimens: $THRESHOLD
Genes with >= 75% specimens: $GOOD_COUNT
Genes with < 75% specimens: $BAD_COUNT

Good genes list: $GOOD_LIST
Bad genes list: $BAD_LIST
EOF

echo "Done!"
echo "Summary written to $SUMMARY"
echo "Good genes: $GOOD_COUNT"
echo "Bad genes: $BAD_COUNT"

