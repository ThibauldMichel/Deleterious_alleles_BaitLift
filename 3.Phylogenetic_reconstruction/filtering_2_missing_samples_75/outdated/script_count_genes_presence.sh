#!/bin/bash
#SBATCH --partition=short


# Directory containing the .fna alignments
DIR=~/scratch/Deleterious_alleles_PNG/CAPTUS/CAPTUS-PNG/05_phylogeny_FastTreeAstral/filtering_gene_trees

# Output files
OUT_MULTIPLE="$DIR/genes_in_multiple_samples.txt"
OUT_SINGLE="$DIR/genes_in_single_sample.txt"

# Empty the output files
> "$OUT_MULTIPLE"
> "$OUT_SINGLE"

# Loop through all .fna files
for f in "$DIR"/*.fna; do
    # Count the number of sample headers in the FASTA file
    n=$(grep -c "^>" "$f")
    
    # Get gene name (filename without path and extension)
    gene=$(basename "$f" .fna)
    
    # Sort into the appropriate list
    if [ "$n" -ge 2 ]; then
        echo "$gene" >> "$OUT_MULTIPLE"
    elif [ "$n" -eq 1 ]; then
        echo "$gene" >> "$OUT_SINGLE"
    fi
done

# Print summary
echo "✅ Genes present in ≥2 samples: $(wc -l < "$OUT_MULTIPLE")"
echo "✅ Genes present in only 1 sample: $(wc -l < "$OUT_SINGLE")"
echo "Results saved to:"
echo "  - $OUT_MULTIPLE"
echo "  - $OUT_SINGLE"

