#!/bin/bash
#SBATCH --job-name="re-calibrate gene trees"
#SBATCH --mem=8G
#SBATCH --partition=short
#SBATCH --export=ALL


# Load conda into this shell session
source /mnt/apps/users/tmichel/conda/etc/profile.d/conda.sh

conda activate snakemake_env



############################################################################
# GENES TREES
# We will recalibrate gene trees for HyPhy as it is sensitive to gene-specific substitution scaling. Per-gene ML branch-length optimization.

ASTRAL="/home/tmichel/scratch/Deleterious_alleles_PNG/CAPTUS/CAPTUS-PNG/05_phylogeny_FastTreeAstral/ASTRAL/species_tree_astral.treefile"
FNA="/home/tmichel/scratch/Deleterious_alleles_PNG/CAPTUS/CAPTUS-PNG/04_alignments/03_trimmed/06_informed/01_coding_NUC/02_NT"

# remove node supports in the ASTRAL tree (numbers before colon), IQ-TREE choke on it.
ASTR_BASENAME=$(basename "$ASTRAL" .treefile)   # species_tree_astral
ASTR_DIR=$(dirname "$ASTRAL")                    # /home/.../ASTRAL
ASTR_OUT="${ASTR_DIR}/${ASTR_BASENAME}_no_support.treefile"

nw_condense "$ASTRAL" \
| sed -E 's/\)([0-9.]+):/\):/g' \
| sed 's/:0\.0/:1e-6/g' \
> "$ASTR_OUT"

################################################################################################



mkdir -p pruned_ASTRAL
mkdir -p recalibrated_gene_trees_from_ASTRAL

# full ASTRAL topology (full taxon set)
for aln in "$FNA"/*.fna; do
  base=$(basename "$aln") ;
  # get base name of gene alignment
  taxa=$(grep ">" "$aln" | sed 's/>//' | sed 's/ .*//') ;
  # create list of taxa in the gene alignment
  nw_prune $ASTR_OUT $taxa  > pruned_ASTRAL/${base}_tree.tre
    # run IQ-TREE on this gene, using pruned tree as starting tree and disabling tree search
  # Note: -te or -t used as starting tree; check your iqtree2 version for flags to disable topology search
  iqtree2 \
     -s "$aln" \
     -t pruned_ASTRAL/${base}_tree.tre \
     -m MFP \
     -redo \
     -pre recalibrated_gene_trees_from_ASTRAL/${base}
  # -s input alignment file
  # -t starting tree (here, pruned ASTRAL tree)
  # -m MFP tells IQ-TREE to use ModelFinder Plus, which tests a set of nucleotide substitution models (JC, HKY, GTR, TN, etc.) with and without rate heterogeneity (+G, +I, etc.) and selects the best fit according to AICc or BIC.
  # -redo allows overwriting previous results
  # -pre specifies the prefix for output files
done
