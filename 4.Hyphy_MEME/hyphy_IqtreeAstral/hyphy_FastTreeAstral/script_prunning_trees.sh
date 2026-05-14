#!/bin/bash

#SBATCH --job-name="prune"
#SBATCH --export=ALL
#SBATCH --partition=short

source activate newick_env


# iqtree -t gene_trees/ACmerged_contig_10072.treefile -prune Hannah_Begonia_baits__ref -pre gene_trees/pruned



# nw_prune your_treefile.tree sample_to_remove > pruned_tree.tree

#nw_prune gene_trees/ACmerged_contig_10072.treefile Hannah_Begonia_baits__ref > gene_trees/ACmerged_contig_10072_pruned.treefile

mkdir gene_trees

while read -r locus; do
  pathtree="gene_trees/${locus}_pruned.treefile"
  pathmsa="cleaned_msa/${locus}_NOSTOP_prunned.fna"
  pathout="output_meme/${locus}_MEME.fna"

nw_prune /home/tmichel/scratch/Deleterious_alleles_PNG/CAPTUS/CAPTUS-PNG/05_phylogeny_FastTreeAstral/genes_tree_FastTree/"$locus".tree Hannah_Begonia_baits__ref > "$pathtree"
done < mfa_files.txt
