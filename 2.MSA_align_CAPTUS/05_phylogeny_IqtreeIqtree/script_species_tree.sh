#!/bin/bash
#SBATCH --job-name=iqtree_species
#SBATCH --export=ALL
#SBATCH --mem=32G
#SBATCH --partition=medium
#SBATCH --output=logs/iqtree_species_%j.out
#SBATCH --error=logs/iqtree_species_%j.err
#SBATCH --cpus-per-task=8

# --- Activate environment ---
source /mnt/apps/users/tmichel/conda/etc/profile.d/conda.sh
conda activate captus_env

# Create output directory if it doesn't exist
mkdir -p ./species_tree

# Define directories
GENE_TREES_DIR="/home/tmichel/scratch/Deleterious_alleles_PNG/hyphy_iqtree_OUTDATED/gene_trees"
OUTPUT_DIR="./species_tree"
mkdir -p $OUTPUT_DIR logs

# Concatenate all gene trees into a single file
cat $GENE_TREES_DIR/*.treefile > $OUTPUT_DIR/all_gene_trees.trees

# Run IQ-TREE to estimate species tree
iqtree3 \
    -t $OUTPUT_DIR/all_gene_trees.trees \
    -con 0.5\
    -m MFP+MERGE \
    -pre $OUTPUT_DIR/species_tree \
    -T AUTO
