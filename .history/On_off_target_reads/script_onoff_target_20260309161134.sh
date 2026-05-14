#!/bin/bash
#SBATCH --job-name="blast"
#SBATCH --export=ALL
#SBATCH --partition=short
#SBATCH --mem=16G

# Load conda into this shell session
source /mnt/apps/users/tmichel/conda/etc/profile.d/conda.sh

conda activate blast_env


# Paths
SUBJ="/home/tmichel/projects/rbge/tmichel/reference_genomes/Bmas.fa"
QUERY="intersected_baits.fasta"
OUTPUT="results_blast.csv"
FILT="results_filtered.csv"

# Step 1:Run blastn

blastn \
  -subject "$SUBJ" \
  -query "$QUERY" \
  -out "$OUTPUT" \
  -outfmt '6 qseqid sseqid pident length sstart send slen evalue bitscore' \
  -evalue 1e-20 \
  -max_target_seqs 50






# Raw results of the blasts:
/home/tmichel/scratch/Deleterious_alleles_PNG_baits/blast_masoniana/results_blast.csv

