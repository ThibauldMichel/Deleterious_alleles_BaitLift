#!/bin/bash
for method in iqTreeASTRAL FastTreeAstral; do
  echo "=== $method ==="
  Path="/home/tmichel/scratch/Deleterious_alleles_PNG/CAPTUS/CAPTUS-PNG/05_phylogeny_${method}"
  nw_labels "$Path/ASTRAL/species_tree_astral.tre" | sort > "$Path/ASTRAL/tree_taxa.txt"
  grep ">" "$Path/MSA_renamed/concatenated_MSA.fna" | sed 's/>//' | sort > "$Path/ASTRAL/msa_taxa.txt"
  n_tree=$(wc -l < "$Path/ASTRAL/tree_taxa.txt")
  n_msa=$(wc -l < "$Path/ASTRAL/msa_taxa.txt")
  echo "Tree taxa: $n_tree"
  echo "MSA taxa: $n_msa"
  diff_count=$(comm -3 "$Path/ASTRAL/tree_taxa.txt" "$Path/ASTRAL/msa_taxa.txt" | wc -l)
  echo "Taxa mismatch: $diff_count"
  echo
done

