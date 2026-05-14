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

awk '$3 >= 98' "$OUTPUT" > "$FILT"

# Step 3: Extract on-target and off-target reads

# Et up environment
source /mnt/apps/users/tmichel/conda/etc/profile.d/conda.sh
conda activate msa_env

SAMTOOLS=$(which samtools)
MAFFT=$(which mafft)
BCFTOOLS=$(which bcftools)

mkdir -p on_off_target_bam

bedtools intersect \
    -abam aligned_reads.bam \
    -b baits_coords.bed > on_target.bam
bedtools intersect \
    -abam aligned_reads.bam \
    -b baits_coords.bed \
    -v > off_target.bam

