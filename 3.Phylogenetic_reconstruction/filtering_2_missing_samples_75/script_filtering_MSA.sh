#!/bin/bash
#SBATCH --job-name=second_filtering
#SBATCH --mem=8G
#SBATCH --partition=short




#!/bin/bash

# Directory containing the MSA files
MSA_DIR="/home/tmichel/scratch/Deleterious_alleles_PNG/CAPTUS/CAPTUS-PNG/05_phylogeny_FastTreeAstral/filtering_1_missing_data_70/filtered_MSA_70"

# Output files in the working directory
SUMMARY_FILE="summary_filtering.txt"
GOOD_GENES="good_genes.csv"
BAD_GENES="bad_genes.csv"

# Step 1: Get total unique specimens across all files
echo "Collecting unique specimen names..."
ALL_SPECIMENS=$(grep -h "^>" "$MSA_DIR"/*.fna | sed 's/^>//' | sort -u)
TOTAL_SPECIMENS=$(echo "$ALL_SPECIMENS" | wc -l)

# Step 2: Get total number of genes (MSA files)
TOTAL_GENES=$(ls "$MSA_DIR"/*.fna | wc -l)

# Step 3: Calculate 75% of specimens (rounded down)
THRESHOLD=$(echo "$TOTAL_SPECIMENS * 0.75" | bc | awk '{print int($1)}')

# Initialize counts and clean any existing output files
GOOD_COUNT=0
BAD_COUNT=0
> "$GOOD_GENES"
> "$BAD_GENES"

echo "Processing each gene file..."
for FILE in "$MSA_DIR"/*.fna; do
    GENE=$(basename "$FILE" .fna)
    NUM_SPECIMENS=$(grep -c "^>" "$FILE")
    if [ "$NUM_SPECIMENS" -ge "$THRESHOLD" ]; then
        echo "$GENE" >> "$GOOD_GENES"
        ((GOOD_COUNT++))
    else
        echo "$GENE" >> "$BAD_GENES"
        ((BAD_COUNT++))
    fi
done

# Step 4: Write summary
{
    echo "Summary of filtering results"
    echo "============================"
    echo "Total unique specimens: $TOTAL_SPECIMENS"
    echo "Total number of genes: $TOTAL_GENES"
    echo "75% of specimens (threshold): $THRESHOLD"
    echo "Genes with >= 75% specimens: $GOOD_COUNT"
    echo "Genes with < 75% specimens: $BAD_COUNT"
    echo
    echo "Lists written to:"
    echo "  $(pwd)/$GOOD_GENES"
    echo "  $(pwd)/$BAD_GENES"
} > "$SUMMARY_FILE"

echo "Done. Summary written to $(pwd)/$SUMMARY_FILE"

