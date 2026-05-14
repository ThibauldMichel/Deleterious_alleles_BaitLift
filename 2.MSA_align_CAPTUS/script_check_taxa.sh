#!/bin/bash




# Step 1 — Extract taxa from each ASTRAL tree

# For IQ-TREE-based ASTRAL
Path_iq="/home/tmichel/scratch/Deleterious_alleles_PNG/CAPTUS/CAPTUS-PNG/05_phylogeny_iqTreeASTRAL"
nw_labels "$Path_iq/ASTRAL/species_tree_astral.tre" | sort > "$Path_iq/ASTRAL/tree_taxa.txt"
echo "IQ-TREE ASTRAL tree taxa count:"
wc -l < "$Path_iq/ASTRAL/tree_taxa.txt"

# For FastTree-based ASTRAL
Path_ft="/home/tmichel/scratch/Deleterious_alleles_PNG/CAPTUS/CAPTUS-PNG/05_phylogeny_FastTreeAstral"
nw_labels "$Path_ft/ASTRAL/species_tree_astral.tre" | sort > "$Path_ft/ASTRAL/tree_taxa.txt"
echo "FastTree ASTRAL tree taxa count:"
wc -l < "$Path_ft/ASTRAL/tree_taxa.txt"


# Step 2 — Extract taxa from concatenated MSAs

# IQ-TREE concatenated alignment
grep ">" "$Path_iq/ASTRAL/../MSA_renamed/concatenated_MSA.fna" | sed 's/>//' | sort > "$Path_iq/ASTRAL/msa_taxa.txt"
echo "IQ-TREE concatenated MSA taxa count:"
wc -l < "$Path_iq/ASTRAL/msa_taxa.txt"

# FastTree concatenated alignment
grep ">" "$Path_ft/ASTRAL/../MSA_renamed/concatenated_MSA.fna" | sed 's/>//' | sort > "$Path_ft/ASTRAL/msa_taxa.txt"
echo "FastTree concatenated MSA taxa count:"
wc -l < "$Path_ft/ASTRAL/msa_taxa.txt"


# Step 3 — Compare taxa sets directly

# Compare for IQ-TREE
echo "Differences between IQ-TREE ASTRAL tree and MSA:"
comm -3 "$Path_iq/ASTRAL/tree_taxa.txt" "$Path_iq/ASTRAL/msa_taxa.txt"

# Compare for FastTree
echo "Differences between FastTree ASTRAL tree and MSA:"
comm -3 "$Path_ft/ASTRAL/tree_taxa.txt" "$Path_ft/ASTRAL/msa_taxa.txt"


echo "IQ-TREE mismatch count:"
comm -3 "$Path_iq/ASTRAL/tree_taxa.txt" "$Path_iq/ASTRAL/msa_taxa.txt" | wc -l

echo "FastTree mismatch count:"
comm -3 "$Path_ft/ASTRAL/tree_taxa.txt" "$Path_ft/ASTRAL/msa_taxa.txt" | wc -l

