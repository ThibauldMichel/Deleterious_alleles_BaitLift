#!/bin/bash
#SBATCH --job-name=print_trees
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --partition=short

# --- Activate environment ---
source /mnt/apps/users/tmichel/conda/etc/profile.d/conda.sh
conda activate snakemake_env




# Run the R script inline
Rscript - << 'EOF'
# Load required libraries
if(!requireNamespace("ape", quietly = TRUE)) install.packages("ape", repos="https://cloud.r-project.org")
if(!requireNamespace("phytools", quietly = TRUE)) install.packages("phytools", repos="https://cloud.r-project.org")

library(ape)
library(phytools)

# Paths to your cleaned species trees
tree1_path <- "/home/tmichel/scratch/Deleterious_alleles_PNG/CAPTUS/CAPTUS-PNG/05_phylogeny_FastTreeAstral/ASTRAL/species_tree_fasttreeastral_clean.tre"
tree2_path <- "/home/tmichel/scratch/Deleterious_alleles_PNG/CAPTUS/CAPTUS-PNG/05_phylogeny_iqTreeASTRAL/ASTRAL/species_tree_iqtreeastral_clean.tre"

# Read the trees
tree1 <- read.tree(tree1_path)
tree2 <- read.tree(tree2_path)

# Keep only common taxa
common_tips <- intersect(tree1$tip.label, tree2$tip.label)
tree1 <- drop.tip(tree1, setdiff(tree1$tip.label, common_tips))
tree2 <- drop.tip(tree2, setdiff(tree2$tip.label, common_tips))

# Force small positive branch lengths
tree1$edge.length[tree1$edge.length == 0] <- 1e-8
tree2$edge.length[tree2$edge.length == 0] <- 1e-8

# Root trees at midpoint
tree1 <- midpoint.root(tree1)
tree2 <- midpoint.root(tree2)

# Create a cophylo object for tanglegram
cophy <- cophylo(tree1, tree2, rotate=TRUE)

# Plot the tanglegram to a PNG
png("species_trees_tanglegram.png", width=2400, height=1600, res=200)
plot(cophy, 
     link.type = "curved",
     link.lwd = 1,
     link.col = "grey50",
     fsize = 0.6,
     mar=c(5,5,5,5))
title("Tanglegram: FastTree+ASTRAL vs IQ-TREE+ASTRAL")
dev.off()

cat("Tanglegram saved to species_trees_tanglegram.png\n")
EOF
