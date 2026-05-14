#!/bin/bash
#SBATCH --job-name=compare_trees
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --partition=short

# --- Activate environment ---
source /mnt/apps/users/tmichel/conda/etc/profile.d/conda.sh
conda activate snakemake_env






# Define paths to your trees
TREE1="/home/tmichel/scratch/Deleterious_alleles_PNG/CAPTUS/CAPTUS-PNG/05_phylogeny_FastTreeAstral/ASTRAL/species_tree_fasttreeastral_clean.tre"
TREE2="/home/tmichel/scratch/Deleterious_alleles_PNG/CAPTUS/CAPTUS-PNG/05_phylogeny_iqTreeASTRAL/ASTRAL/species_tree_iqtreeastral_clean.tre"

# Create and run the Python comparison inline
python << EOF
from ete3 import Tree

t1_path = "$TREE1"
t2_path = "$TREE2"

print("Loading trees:")
print(t1_path)
print(t2_path)

# Load trees
t1 = Tree(t1_path)
t2 = Tree(t2_path)

# Compare topologies
rf, max_rf, common_leaves, *_ = t1.robinson_foulds(t2)

print("\\n=== Topology Comparison ===")
print(f"RF distance: {rf}")
print(f"Normalized RF distance: {rf / max_rf:.4f}")
print(f"Common leaves: {len(common_leaves)}")
print(f"Same topology? {'YES' if rf == 0 else 'NO'}")
EOF
