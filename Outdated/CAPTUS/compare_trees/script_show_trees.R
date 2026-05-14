library(ape)
library(phytools)

# Paths to your cleaned species trees
tree1_path <- "/home/thibauld/Documents/Bioinformatics/Deleterious_alleles_pipeline/Deleterious_alleles_PNG/CAPTUS/CAPTUS-PNG/05_phylogeny_FastTreeAstral/ASTRAL/species_tree_fasttreeastral_clean.tre"
tree2_path <- "/home/thibauld/Documents/Bioinformatics/Deleterious_alleles_pipeline/Deleterious_alleles_PNG/CAPTUS/CAPTUS-PNG/05_phylogeny_iqTreeASTRAL/ASTRAL/species_tree_iqtreeastral_clean.tre"

# Read the trees
tree1 <- read.tree(tree1_path)
tree2 <- read.tree(tree2_path)

# Replace NaN and zero branch lengths with a tiny positive value
tree1$edge.length[is.na(tree1$edge.length) | tree1$edge.length == 0] <- 1e-8
tree2$edge.length[is.na(tree2$edge.length) | tree2$edge.length == 0] <- 1e-8

# Keep only common taxa
common_tips <- intersect(tree1$tip.label, tree2$tip.label)
tree1 <- drop.tip(tree1, setdiff(tree1$tip.label, common_tips))
tree2 <- drop.tip(tree2, setdiff(tree2$tip.label, common_tips))

# Root trees at the node including Hillebrandia and B.antsiranensis175
root_taxa <- c("Hillebrandia", "B.antsiranensis175")
tree1 <- root(tree1, node = getMRCA(tree1, root_taxa))
tree2 <- root(tree2, node = getMRCA(tree2, root_taxa))

# Create a cophylo object for tanglegram
cophy <- cophylo(tree1, tree2, rotate=TRUE)

# Plot tanglegram
png("/home/thibauld/Documents/Bioinformatics/Deleterious_alleles_pipeline/Deleterious_alleles_PNG/CAPTUS/CAPTUS-PNG/compare_trees/species_trees_tanglegram.png", width=2000, height=8000, res=200)
plot(cophy, 
     link.type = "curved",
     link.lwd = 1,
     link.col = "grey50",
     fsize = 0.5,
     mar=c(5,5,5,5))
title("Tanglegram: FastTree+ASTRAL vs IQ-TREE+ASTRAL")
dev.off()
