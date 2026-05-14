#!/bin/bash
#SBATCH --job-name="blast"
#SBATCH --export=ALL
#SBATCH --partition=short
#SBATCH --mem=16G

# Load conda into this shell session
source /mnt/apps/users/tmichel/conda/etc/profile.d/conda.sh
conda activate msa_env

# Set up environment
BEDTOOLS=$(which bedtools)
SAMTOOLS=$(which samtools)
MAFFT=$(which mafft)
BCFTOOLS=$(which bcftools)


# Paths
SUBJ="/home/tmichel/projects/rbge/tmichel/reference_genomes/Bmas.fa"
QUERY="/home/tmichel/projects/rbge/tmichel/reference_genomes/Hannah_Begonia_baits_edited.fasta"
OUTPUT="results_blast.csv"
FILT="results_filtered.csv"
# Directory containing BAM files
BAM_DIR="/home/tmichel/scratch/ROH-pipeline_PNG/mapped"
# Bait coordinates BED file
BAITS_BED="/home/tmichel/scratch/ROH-pipeline_PNG/baits_coords.bed"
# Output directories
OUT_DIR="${BAM_DIR}/on_off_target"

# Step 3: Extract on-target and off-target reads

mkdir -p on_target_bam
mkdir -p off_target_bam

$BEDTOOLS intersect \
    -abam aligned_reads.bam \
    -b baits_coords.bed > on_target.bam
    
$BEDTOOLS intersect \
    -abam aligned_reads.bam \
    -b baits_coords.bed \
    -v > off_target.bam
    
# Count the reads    
$SAMTOOLS view -c on_target.bam
$SAMTOOLS view -c off_target.bam


# Loop through all BAM files
for BAM in "$BAM_DIR"/*.sorted.bam; do
    BASENAME=$(basename "$BAM" .sorted.bam)
    echo "Processing $BASENAME"

    # 1️⃣ Extract on-target and off-target reads
    bedtools intersect -abam "$BAM" -b "$BAITS_BED" > "$OUT_DIR/on_target_bam/${BASENAME}_on.bam"
    bedtools intersect -abam "$BAM" -b "$BAITS_BED" -v > "$OUT_DIR/off_target_bam/${BASENAME}_off.bam"

    # 2️⃣ Count reads
    ON_READS=$(samtools view -c "$OUT_DIR/on_target_bam/${BASENAME}_on.bam")
    OFF_READS=$(samtools view -c "$OUT_DIR/off_target_bam/${BASENAME}_off.bam")
    echo "On-target reads: $ON_READS"
    echo "Off-target reads: $OFF_READS"




############################################

# On-target variants
bedtools intersect -a variants.vcf -b baits_coords.bed > on_target_variants.vcf

# Off-target variants
bedtools intersect -a variants.vcf -b baits_coords.bed -v > off_target_variants.vcf