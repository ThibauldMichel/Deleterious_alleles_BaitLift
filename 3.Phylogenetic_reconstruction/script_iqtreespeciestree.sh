#!/bin/bash
#SBATCH --job-name=iqtree_speciesGT
#SBATCH --export=ALL
#SBATCH --mem=200G
#SBATCH --partition=long
#SBATCH --cpus-per-task=16


# -------------------------
# 1. Activate conda environment
# -------------------------
source /mnt/apps/users/tmichel/conda/etc/profile.d/conda.sh
conda activate captus_env

# -------------------------
# 2. Set directories
# -------------------------
# Directory containing all good gene MSAs (*.fna)
MSA_DIR=~/scratch/Deleterious_alleles_PNG/CAPTUS/05_phylogeny_FastTreeAstral/good_genes_MSA

# Output directory
OUTDIR=~/scratch/Deleterious_alleles_PNG/CAPTUS/05_phylogeny_FastTreeAstral/iqtree_species_tree

# -------------------------
# 3. Prepare directories and files
# -------------------------


# Create output directory if it doesn't exist
mkdir -p "$OUTDIR"

# Move into output directory
cd "$OUTDIR"

echo "🔹 Concatenating all .fna MSAs into supermatrix.fna ..."

# Concatenate all MSA files (same taxa order assumed)
cat "$MSA_DIR"/*.fna > supermatrix.fna

# -------------------------
# 4. Make the tree
# -------------------------



echo "🔹 Running IQ-TREE on concatenated supermatrix..."

iqtree3 \
    -s "$MSA_DIR" \
    -m MFP \
    -bb 1000 \
    -alrt 1000 \
    -nt AUTO \
    -pre species_tree \
    -T AUTO

echo "✅ Species tree ML inference complete!"
echo "   Results are in: $OUTDIR"

# -s supermatrix.fna — Input alignment
# -m MFP — ModelFinder Plus (automatic model selection)
# -bb 1000 — Ultrafast bootstrap (UFBoot)
# -alrt 1000 — SH-aLRT branch test (Runs 1,000 replicates of the SH-aLRT test)
# -nt AUTO — Automatic detection of CPU threads
# -pre species_tree — Output prefix




