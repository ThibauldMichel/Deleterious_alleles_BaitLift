#!/bin/bash
#SBATCH --job-name=iqtree_rescale
#SBATCH --export=ALL
#SBATCH --mem=128G
#SBATCH --partition=long
#SBATCH --output=logs/iqtree_%A_%a.out
#SBATCH --error=logs/iqtree_%A_%a.err
#SBATCH --cpus-per-task=16


# --- Activate environment ---
source /mnt/apps/users/tmichel/conda/etc/profile.d/conda.sh
conda activate captus_env

#----------------------------:q

# 1. Prepare a concatenated alignment (supermatrix) of all genes, or at least all taxa present in the ASTRAL tree.

ASTRAL_TREE="/home/tmichel/scratch/Deleterious_alleles_PNG/CAPTUS/CAPTUS-PNG/05_phylogeny_FastTreeAstral/ASTRAL/species_tree_astral.tre"

ALIGNMENTS_DIR="/home/tmichel/scratch/Deleterious_alleles_PNG/CAPTUS/CAPTUS-PNG/05_phylogeny_FastTreeAstral/good_genes_MSA"

CONCAT_ALIGNEMENT="./concatenated_MSA.fna"

nw_labels $ASTRAL_TREE > taxa_in_tree.txt

AMAS.py concat \
 -i $ALIGNMENTS_DIR/*.fna \
 -f fasta \
 -d dna \
 -p partitions.txt \
 -t concatenated_MSA.fna

# -i → input alignments (*.fna)
# -f fasta → input format
# -d dna → data type (DNA)
# -u fasta → output format (can also use phylip or nexus)
# -t → output concatenated alignment file
#-p → output partition file (useful for IQ-TREE)

# 1.5. Reformat the partition file for iqtree
sed -i 's/^/DNA, /' partitions.txt


# 2. Use IQ-TREE to rescale branch lengths of your species tree while keeping its topology:

iqtree3 -s concatenated_MSA.fna \
        -te $ASTRAL_TREE \
        -p partitions.txt \
        -m MFP+MERGE \
        --tree-fix \
        -nt AUTO \
        -pre species_tree_rescaled

# -te: provide your species tree (topology).
# --tree-fix: do not change topology, only optimize branch lengths.
# Output will be one rescaled tree: species_tree_rescaled.treefile.
# -m MFP+MERGE → IQ-TREE automatically finds the best model for each partition and merges similar ones (avoiding overfitting).
# --tree-fix → keeps your ASTRAL topology fixed while optimizing branch lengths.
# Output (species_tree_rescaled.treefile) will be your rescaled ASTRAL tree.
