#!/bin/bash
#SBATCH --job-name="SIFT_annotation"
#SBATCH --mem=8G
#SBATCH --partition=medium
#SBATCH --export=ALL

# Load conda into this shell session
source /mnt/apps/users/tmichel/conda/etc/profile.d/conda.sh

conda activate bcftools_env




# Liftover VCF from bait coordinates to genome coordinates, then prepare for SIFT4G annotation

VCF_INPUT="${1}"
MAPPING_TSV="${2}"
BMAS_FASTA="${3}"
BMAS_GFF="${4}"
OUTPUT_DIR="${5:-.}"

if [ -z "$VCF_INPUT" ] || [ -z "$MAPPING_TSV" ]; then
    cat <<EOF
Usage: $0 <input.vcf.gz> <mapping.tsv> <Bmas.fa> <Bmas.gff> [output_dir]

Steps:
  1. Lift VCF from bait coordinates to genome coordinates
  2. Validate lifted VCF (sanity checks)
  3. Build SIFT4G database
  4. Annotate with SIFT4G
EOF
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "[1/4] Lifting VCF from bait to genome coordinates..."
python3 "$SCRIPT_DIR/liftover_bait_to_genome.py" \
    "$VCF_INPUT" \
    "$MAPPING_TSV" \
    "$OUTPUT_DIR/all_renamed.genome_coords.vcf.gz"

echo "[2/4] Sorting and validating lifted VCF..."
bcftools sort -O z -o "$OUTPUT_DIR/all_renamed.genome_coords.sorted.vcf.gz" \
    "$OUTPUT_DIR/all_renamed.genome_coords.vcf.gz"
bcftools index "$OUTPUT_DIR/all_renamed.genome_coords.sorted.vcf.gz"

# Quick validation: check that chromosome names match reference
echo "[*] Validating chromosome names..."
bcftools view -H "$OUTPUT_DIR/all_renamed.genome_coords.sorted.vcf.gz" | \
    awk '{print $1}' | sort -u > "$OUTPUT_DIR/lifted_chroms.txt"

echo "Chromosomes in lifted VCF:"
head "$OUTPUT_DIR/lifted_chroms.txt"

echo ""
echo "[3/4] Building SIFT4G database..."
SIFT4G_Create_Genomic_DB \
    -genome "$BMAS_FASTA" \
    -dbName "$OUTPUT_DIR/sift4g_bmas" \
    -gff "$BMAS_GFF"

echo ""
echo "[4/4] Running SIFT4G annotation..."
SIFT4G_Annotator \
    -c "$OUTPUT_DIR/all_renamed.genome_coords.sorted.vcf.gz" \
    -d "$OUTPUT_DIR/sift4g_bmas" \
    -r "$OUTPUT_DIR/sift_results/"

echo ""
echo "[+] SIFT4G annotation complete!"
echo "    Results: $OUTPUT_DIR/sift_results/"
echo ""
echo "    Key files:"
echo "      - $OUTPUT_DIR/all_renamed.genome_coords.sorted.vcf.gz (lifted VCF)"
echo "      - $OUTPUT_DIR/sift_results/SIFT_results.txt"
echo "      - $OUTPUT_DIR/sift_results/*_exonicSIFT.vcf (annotated VCF)"
