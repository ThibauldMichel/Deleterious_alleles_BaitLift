#!/bin/bash
#SBATCH --job-name="ASTRAL_species_tree"
#SBATCH --export=ALL
#SBATCH --mem=64G
#SBATCH --partition=long

source /mnt/apps/users/tmichel/conda/etc/profile.d/conda.sh
conda activate FastTree_env

Path="/home/tmichel/scratch/Deleterious_alleles_PNG/CAPTUS/CAPTUS-PNG/05_phylogeny_FastTreeAstral"
JAR="/home/tmichel/ASTRAL/Astral/astral.5.7.8.jar"
GOOD_GENES="$Path/filtering_2_missing_samples_75/good_genes.csv"
CONCAT_TREE="$Path/ASTRAL/gene_trees_good.tre"
SPECIES_TREE="$Path/ASTRAL/species_tree_astral.tre"

# Make a folder for species tree
mkdir -p "$Path/ASTRAL"

# Concatenate only the good gene trees
> "$CONCAT_TREE"  # empty or create the output file

while read -r gene; do
    TREE_FILE="$Path/genetree_from_filtered_MSA/${gene}_FastTree.tree"
    if [[ -f "$TREE_FILE" ]]; then
        cat "$TREE_FILE" >> "$CONCAT_TREE"
    else
        echo "Warning: Tree file not found for $gene" >&2
    fi
done < "$GOOD_GENES"

# Quick sanity check: count trees (one ';' per Newick tree)
echo "Number of gene trees:" $(grep -c ';' "$CONCAT_TREE")

# Run ASTRAL
java -Xmx32G -jar "$JAR" \
    -i "$CONCAT_TREE" \
    -o "$SPECIES_TREE" \
    -t 1 \
    > "$Path/ASTRAL/astral_run.log" 2>&1



# java -Xmx32G  Allocates 32 GB of memory for the Java Virtual Machine running ASTRAL. Adjust this based on your available resources.
# -jar /path/to/astral.5.7.8.jar        Specifies the ASTRAL executable JAR file to run.
# -i ./gene_trees_all.tre       Input file containing all individual gene trees (one per line or concatenated). ASTRAL will use these to infer the species tree.
# -o ./species_tree_astral.tre  Output file where the inferred species tree will be written.
# -t 3  Requests that ASTRAL output additional information — this is the “output option flag”, where:
        # -t 1 → only species tree (default)
        # -t 2 → species tree + local posterior probabilities (LPPs)
        # -t 3 → species tree + branch quartet scores + local posterior probabilities (the most detailed output)
# > astral_run.log 2>&1 Redirects both standard output (stdout) and standard error (stderr) into the file astral_run.log, so the run’s full log is saved.




