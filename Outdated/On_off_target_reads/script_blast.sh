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
QUERY="/home/tmichel/projects/rbge/tmichel/reference_genomes/Hannah_Begonia_baits_edited.fasta"
OUTPUT="results_blast.csv"
FILT="results_filtered.csv"
BAM="/home/tmichel/scratch/ROH-pipeline_PNG/mapped"

# Step 1:Run blastn

blastn \
  -subject "$SUBJ" \
  -query "$QUERY" \
  -out "$OUTPUT" \
  -outfmt '6 qseqid sseqid pident length sstart send slen evalue bitscore' \
  -evalue 1e-20 \
  -max_target_seqs 50 > "$OUTPUT"


# Step 2: Filter results for matches with >= 98% identity

awk '$3 >= 90' "$OUTPUT" > "$FILT"


# We filter at 90% as we want all the baits mapped (1200 with this score) and we wish to see what is "off target". So even if baits are matching different locations, we count these locations as "on target" as any reads NOT coming from these sites is surely off-target.
