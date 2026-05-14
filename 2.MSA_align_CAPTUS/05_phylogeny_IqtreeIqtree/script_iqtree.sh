#!/bin/bash

#SBATCH --job-name="iqtree_PNG"
#SBATCH --export=ALL
#SBATCH --mem=500G
#SBATCH --partition=himem


# To make individual genes trees and create them in the gene_trees directory

mkdir gene_trees

for aln in cleaned_msa/*.fna; do
    basename=$(basename "$aln" .fna)
    basename=${basename%_NOSTOP}
    iqtree2 --prefix "gene_trees/${basename}" -s "${aln}" -m MFP -B 1000 --alrt 1000 --bnni --safe --seed 42 --threads AUTO

done

#iqtree2
#Calls the IQ-TREE version 2 executable.

#--prefix "gene_trees/${basename}"
#Sets the prefix for output files. All results will be saved in files starting with gene_trees/${basename}.
#basename is likely a variable holding the name of the alignment file without its extension.

#-s "${aln}"
#Specifies the input alignment file, where ${aln} is a variable pointing to the path of the alignment (in FASTA, PHYLIP, etc.).

#-m MFP
# Tells IQ-TREE to automatically find the best-fit model of sequence evolution using ModelFinder Plus (MFP).

# -B 1000
# Performs 1000 ultrafast bootstrap replicates, a fast method for assessing node support in the tree.

# --alrt 1000
# Performs 1000 SH-aLRT (approximate likelihood-ratio test) replicates for branch support.

# --bnni
# Optimizes the ultrafast bootstrap trees using branch length and nearest-neighbor interchange to reduce overestimation of support.

# --safe
# Enables a safe mode that automatically activates various sanity checks to avoid common mistakes (e.g., wrong file formats, wrong tree rooting, etc.).

# --seed 42
# Sets the random seed for reproducibility. Using the same seed ensures the same random processes (like bootstrap sampling) produce the same results.

# --threads AUTO
# Tells IQ-TREE to automatically detect and use an appropriate number of CPU threads for parallel computation.

################################################################################################################

# Tree Estimation (Maximum Likelihood Tree)

# IQ-TREE starts by estimating a single maximum likelihood (ML) tree from the alignment. This is not a set of 1000 ML trees. Only one best-scoring ML tree is inferred using the optimal model selected by -m MFP.
# Then, it estimates branch support using resampling methods.


# Bootstrap Support: -B 1000

# This tells IQ-TREE to compute 1000 ultrafast bootstrap replicates.
# - Ultrafast bootstrapping (UFBoot) is a fast approximation of traditional bootstrapping.
# - It generates 1000 pseudo-replicate datasets by resampling alignment columns.
# - For each, a tree is built quickly using fast heuristics.
# - The branch support values are derived from how consistently a given clade appears across these 1000 replicate trees.
# - UFBoot is designed to be statistically efficient and computationally fast.



# SH-aLRT Support: --alrt 1000

# This requests 1000 replicates of the SH-like approximate likelihood ratio test for branch support. It is independent of bootstrapping and provides a different statistical measure of support.
# How SH-aLRT Works:
# - For each internal branch in the ML tree, the method compares the likelihood of:
# - The best ML tree topology (with that branch present)
# - Alternative topologies with that branch removed or rearranged
# - It does this 1000 times (by generating test data sets with parametric resampling).
# - The result is a p-value-like support value between 0–100%, showing how likely that branch is real.
