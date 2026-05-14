#!/bin/bash
#SBATCH --job-name="SIFT_annotation"
#SBATCH --mem=8G
#SBATCH --partition=medium
#SBATCH --export=ALL

# ==============================================================================
# SIFT4G annotation pipeline
#
# Usage (direct execution or sbatch):
#   sbatch run_sift_pipeline.sh [output_dir]
#
# All input paths are hardcoded below. OUTPUT_DIR defaults to the
# current directory if not supplied as $1.
#
# Steps:
#   1. Liftover VCF from bait coordinates to genome coordinates
#   2. Sort, index, and validate the lifted VCF
#   3. Build SIFT4G genomic database from reference FASTA + GFF
#   4. Annotate with SIFT4G
# ==============================================================================

# ------------------------------------------------------------------------------
# Load conda environment
# ------------------------------------------------------------------------------
source /mnt/apps/users/tmichel/conda/etc/profile.d/conda.sh
conda activate bcftools_env

# ------------------------------------------------------------------------------
# Hardcoded input paths — edit these as needed
# ------------------------------------------------------------------------------
VCF_INPUT="/home/tmichel/scratch/Deleterious_alleles_PNG_baits/6.Variant_calling_annotations/all.vcf.gz"
MAPPING_TSV="/home/tmichel/scratch/Deleterious_alleles_PNG_baits/1.Baits_to_genome/bait_locus_annotation.tsv"
BMAS_FASTA="/home/tmichel/projects/rbge/tmichel/reference_genomes/Bmas.fa"
BMAS_GFF="/home/tmichel/projects/rbge/tmichel/reference_genomes/Bmas.gff"
SIFT4G="/home/tmichel/sift4g//bin/sift4g"

# ------------------------------------------------------------------------------
# Output directory: use $1 if supplied, otherwise current directory
# ------------------------------------------------------------------------------
OUTPUT_DIR="${1:-.}"

# ------------------------------------------------------------------------------
# Validate inputs before doing any work
# ------------------------------------------------------------------------------
missing=0
for f in "$VCF_INPUT" "$MAPPING_TSV" "$BMAS_FASTA" "$BMAS_GFF" "$SIFT4G_JAR"; do
    if [ ! -f "$f" ]; then
        echo "[!] Missing required file: $f" >&2
        missing=1
    fi
done
if [ "$missing" -eq 1 ]; then
    echo "[!] One or more required files are missing. Aborting." >&2
    exit 1
fi

# Verify the input VCF is tabix-indexed (required for bcftools operations)
if [ ! -f "${VCF_INPUT}.tbi" ] && [ ! -f "${VCF_INPUT}.csi" ]; then
    echo "[!] Input VCF has no tabix/CSI index. Run: bcftools index $VCF_INPUT" >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# Resolve the directory containing this script so Python script is found
# regardless of where the job is launched from
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Convenience variables for intermediate/final files
LIFTED_VCF="$OUTPUT_DIR/all_renamed.genome_coords.vcf.gz"
SORTED_VCF="$OUTPUT_DIR/all_renamed.genome_coords.sorted.vcf.gz"
FAILED_VCF="$OUTPUT_DIR/failed_liftover.vcf"
LIFTED_CHROMS="$OUTPUT_DIR/lifted_chroms.txt"
SIFT_DB="$OUTPUT_DIR/sift4g_bmas"
SIFT_RESULTS="$OUTPUT_DIR/sift_results"

# ------------------------------------------------------------------------------
# Step 1 — Liftover: bait coordinates → genome coordinates
# ------------------------------------------------------------------------------
echo ""
echo "[1/4] Lifting VCF from bait to genome coordinates..."
python3 "./liftover_bait_to_genome.py" \
    "$VCF_INPUT" \
    "$MAPPING_TSV" \
    "$LIFTED_VCF" \
    "$FAILED_VCF"

if [ $? -ne 0 ]; then
    echo "[!] Liftover failed. Aborting." >&2
    exit 1
fi

# ------------------------------------------------------------------------------
# Step 2 — Sort, index, and validate lifted VCF
# ------------------------------------------------------------------------------
echo ""
echo "[2/4] Sorting and indexing lifted VCF..."
bcftools sort -O z -o "$SORTED_VCF" "$LIFTED_VCF"

if [ $? -ne 0 ]; then
    echo "[!] bcftools sort failed. Aborting." >&2
    exit 1
fi

bcftools index "$SORTED_VCF"

# Validate: check chromosome names in the lifted VCF
echo "[*] Validating chromosome names in lifted VCF..."
bcftools view -H "$SORTED_VCF" | awk '{print $1}' | sort -u > "$LIFTED_CHROMS"
echo "    Chromosomes present in lifted VCF:"
sed 's/^/      /' "$LIFTED_CHROMS"

# Cross-check against reference FASTA sequence names
REF_CHROMS="$OUTPUT_DIR/ref_chroms.txt"
grep '^>' "$BMAS_FASTA" | sed 's/^>//' | awk '{print $1}' | sort -u > "$REF_CHROMS"

UNMATCHED=$(comm -23 "$LIFTED_CHROMS" "$REF_CHROMS" | wc -l)
if [ "$UNMATCHED" -gt 0 ]; then
    echo "[!] WARNING: $UNMATCHED chromosome(s) in lifted VCF not found in reference FASTA:" >&2
    comm -23 "$LIFTED_CHROMS" "$REF_CHROMS" | sed 's/^/      /' >&2
fi

# ------------------------------------------------------------------------------
# Step 3 — Build SIFT4G genomic database
# ------------------------------------------------------------------------------
echo ""
echo "[3/4] Building SIFT4G database..."
"$SIFT4G" SIFT4G_Create_Genomic_DB \
    -genome "$BMAS_FASTA" \
    -dbName "$SIFT_DB" \
    -gff "$BMAS_GFF"

if [ $? -ne 0 ]; then
    echo "[!] SIFT4G database build failed. Aborting." >&2
    exit 1
fi

# ------------------------------------------------------------------------------
# Step 4 — Annotate with SIFT4G
# ------------------------------------------------------------------------------
echo ""
echo "[4/4] Running SIFT4G annotation..."
mkdir -p "$SIFT_RESULTS"

java -jar "$SIFT4G_JAR" SIFT4G_Annotator \
    -c "$SORTED_VCF" \
    -d "$SIFT_DB" \
    -r "$SIFT_RESULTS"

if [ $? -ne 0 ]; then
    echo "[!] SIFT4G annotation failed." >&2
    exit 1
fi

# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------
echo ""
echo "[+] SIFT4G annotation complete!"
echo ""
echo "    Key output files:"
echo "      $SORTED_VCF"
echo "        Lifted + sorted VCF in genome coordinates"
echo "      $FAILED_VCF"
echo "        Variants that could not be lifted (inspect for QC)"
echo "      $SIFT_RESULTS/SIFT_results.txt"
echo "        SIFT4G scores for all annotated variants"
echo "      $SIFT_RESULTS/*_exonicSIFT.vcf"
echo "        Annotated VCF (exonic variants only)"
