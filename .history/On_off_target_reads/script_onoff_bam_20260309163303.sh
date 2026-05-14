#!/bin/bash
#SBATCH --job-name="blast"
#SBATCH --export=ALL
#SBATCH --partition=short
#SBATCH --mem=16G

# Load conda into this shell session
source /mnt/apps/users/tmichel/conda/etc/profile.d/conda.sh
conda activate msa_env

# Set up environment
SAMTOOLS=$(which samtools)
MAFFT=$(which mafft)
BCFTOOLS=$(which bcftools)


# Paths
SUBJ="/home/tmichel/projects/rbge/tmichel/reference_genomes/Bmas.fa"
QUERY="/home/tmichel/projects/rbge/tmichel/reference_genomes/Hannah_Begonia_baits_edited.fasta"
OUTPUT="results_blast.csv"
FILT="results_filtered.csv"
BAM="/home/tmichel/scratch/ROH-pipeline_PNG/mapped"

# Step 3: Extract on-target and off-target reads

mkdir -p on_off_target_bam

bedtools intersect \
    -abam aligned_reads.bam \
    -b baits_coords.bed > on_target.bam
    
bedtools intersect \
    -abam aligned_reads.bam \
    -b baits_coords.bed \
    -v > off_target.bam
