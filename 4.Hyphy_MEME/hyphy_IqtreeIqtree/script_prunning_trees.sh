#!/bin/bash

#SBATCH --job-name="prune"
#SBATCH --export=ALL
#SBATCH --partition=short

source activate newick_env


# iqtree -t gene_trees/ACmerged_contig_10072.treefile -prune Hannah_Begonia_baits__ref -pre gene_trees/pruned



# nw_prune your_treefile.tree sample_to_remove > pruned_tree.tree

#nw_prune gene_trees/ACmerged_contig_10072.treefile Hannah_Begonia_baits__ref > gene_trees/ACmerged_contig_10072_pruned.treefile



while read -r locus; do
  pathtree="gene_trees/${locus}_pruned.treefile"
  pathmsa="cleaned_msa/${locus}_NOSTOP_prunned.fna"
  pathout="output_meme/${locus}_MEME.fna"

nw_prune gene_trees/"$locus".treefile Hannah_Begonia_baits__ref > "$pathtree"
#  hyphy meme \
#	--alignment "$pathmsa" \
#	--tree "$pathtree" \
#	--branches All \
#	--pvalue 0.1 \
#	--output "$pathout" \
#	--full-model Yes ;
done < mfa_files.txt
