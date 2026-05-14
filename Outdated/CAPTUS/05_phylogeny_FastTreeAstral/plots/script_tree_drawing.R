#install.packages("ape")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("treeio")
BiocManager::install("ggtree")














library(ape)
library(treeio)
library(ggtree)




path <- "/home/thibauld/Documents/Bioinformatics/Deleterious_alleles_pipeline/Deleterious_alleles_PNG/CAPTUS/05_phylogeny_FastTreeAstral"

# Loading tree
tree <- read.tree(file.path(path, "ASTRAL/species_tree_before_rescale.treefile"))

# Reboot graphics
dev.off()

# Preview of the tree
plot(tree)

# Look at the labels
tree$tip.label

# Show node labels (support values)
tree$node.label

# Store node labels separately
node_support <- tree$node.label

# Remove node labels temporarily so no NAs interfere
tree$node.label <- NULL

# Replace any NA branch lengths
tree$edge.length[is.na(tree$edge.length)] <- 0.001

# Plot tree
plot(tree, cex = 0.6)

# Add node labels manually
nodelabels(node_support, cex = 0.5)




# Root and ladderize
taxa_to_root <- c("Hillebrandia", "B.antsiranensis175")
rooted_tree <- root(tree, outgroup = taxa_to_root, resolve.root = TRUE)
#rooted_tree <- ladderize(rooted_tree)

# Prune problematic samples
tips_to_remove <- c("B.S.Pet.93", "B.brachybotrys147")
pruned_tree <- drop.tip(rooted_tree, tip = tips_to_remove)

# Scale branch lengths proportionally
#scale_factor <- 1
#scaled_tree <- rooted_tree
#scaled_tree$edge.length <- scaled_tree$edge.length * scale_factor

# Set x-lim to match scaled tree
#xlim_max <- max(node.depth.edgelength(scaled_tree)) * 1.05


# Save to PNG
png(file.path(path, "plots/rooted_tree_scaled_branches.png"),
    width = 2000,
    height = ntips * 30,
    res = 300)

#plot.phylo(pruned_tree)

plot.phylo(pruned_tree,
           type = "phylogram",
           cex = 0.5,
           direction = "rightwards",
           show.tip.label = TRUE,
           no.margin = FALSE,
#           x.lim = c(0, 1) 
)




library(treeio)
library(ggtree)

# Read tree with bootstrap values
tree <- read.tree("your_tree_file.tree")  














# Add bootstrap values
# Assume your tree has 'node.label' with support values
nodelabels(pruned_tree$node.label, 
           adj = c(1, -0.2),  # adjust position: c(horizontal, vertical)
           frame = "none",    # no box around label
           cex = 0.6,         # size of numbers
           col = "blue")      # color of bootstrap numbers
dev.off()
