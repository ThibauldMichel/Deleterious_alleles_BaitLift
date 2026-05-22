#!/bin/bash
#SBATCH --job-name="Filtering snpeff"
#SBATCH --mem=40G
#SBATCH --partition=medium
#SBATCH --export=ALL

# Load conda environment
# ------------------------------------------------------------------------------
source /mnt/apps/users/tmichel/conda/etc/profile.d/conda.sh
conda activate bcftools_env

# ------------------------------------------------------------------------------
# Input paths — edit these as needed
# ------------------------------------------------------------------------------
# snpEff installation directory (contains snpEff.jar and snpEff.config)
SNPEFF_DIR="$HOME/software/snpEff"
SNPEFF_JAR="$SNPEFF_DIR/snpEff.jar"
SNPEFF_CFG="$SNPEFF_DIR/snpEff.config"

# SIFT4G_Annotator JAR
SIFT4G_JAR="$HOME/software/SIFT4G_Annotator/SIFT4G_Annotator.jar"

# Genome name — must match the entry in snpEff.config
GENOME_NAME="Bmas"

SNPEFF_DATA_DIR="$SNPEFF_DIR/data"

# ------------------------------------------------------------------------------
# Output directory: use $1 if supplied, otherwise current directory
# ------------------------------------------------------------------------------
OUTPUT_DIR="${1:-.}"
mkdir -p "$OUTPUT_DIR"
 
 #==============================================================================
# Remove WARNING_REF_DOES_NOT_MATCH_GENOME variants
# ==============================================================================
echo ""
echo "script 5: Filtering WARNING_REF_DOES_NOT_MATCH_GENOME variants..."

SNPEFF="${OUTPUT_DIR}/all_annotated_snpeff.vcf"

SNPEFF_VCF_FILTERED="${OUTPUT_DIR}/all_annotated_snpeff.filtered.vcf"

bcftools filter \
    -e 'INFO/ANN ~ "WARNING_REF_DOES_NOT_MATCH_GENOME"' \
    -O v \
    -o "$SNPEFF_VCF_FILTERED" \
    "$SNPEFF_VCF"

if [ $? -ne 0 ]; then
    echo "[!] Filtering step failed." >&2
    exit 1
fi

# Quick sanity check
N_BEFORE=$(bcftools view -H "$SNPEFF_VCF" | wc -l)
N_AFTER=$(bcftools view -H "$SNPEFF_VCF_FILTERED" | wc -l)
echo "    Variants before filtering : $N_BEFORE"
echo "    Variants after filtering  : $N_AFTER"
echo "    Removed                   : $((N_BEFORE - N_AFTER))"
