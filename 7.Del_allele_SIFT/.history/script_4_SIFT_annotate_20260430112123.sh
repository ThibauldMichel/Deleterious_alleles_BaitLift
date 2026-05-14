#!/bin/bash
#SBATCH --job-name="SIFT_snpEff_annotate"
#SBATCH --mem=40G
#SBATCH --partition=medium
#SBATCH --export=ALL

# ==============================================================================
# Script 3c — SIFT4G + snpEff annotation
#
# Usage:
#   sbatch script_3c_annotate.sh [output_dir]
#   OUTPUT_DIR defaults to current directory if not supplied as $1.
#
# Prerequisites:
#   script_3a_build_snpeff_db.sh must have completed successfully.
#   script_3b_build_regions.sh must have completed successfully.
#   The sorted, lifted VCF from script_1 must exist at:
#     $OUTPUT_DIR/all_renamed.genome_coords.sorted.vcf.gz
#
# What this script does:
#
#   Step 4 — SIFT4G annotation (SIFT4G_Annotator.jar)
#   -------------------------------------------------
#   Annotates each variant in the VCF with a SIFT score and prediction
#   (TOLERATED / DELETERIOUS) by looking up the per-scaffold .regions files
#   built in script_3b.
#
#   SIFT4G_Annotator.jar does not accept bgzipped (.vcf.gz) input, so the
#   sorted VCF is first decompressed to a plain .vcf. The uncompressed copy
#   is removed after annotation to save space.
#
#   Flags used:
#     -c        headless / command-line mode (no GUI)
#     -i <vcf>  input VCF (must be plain, not bgzipped)
#     -d <dir>  snpEff data/Bmas directory containing the .regions files
#     -r <dir>  output directory for annotated VCF
#     -t        multithreaded mode
#
#   Output file: $SIFT_RESULTS/*_SIFTannotations.vcf
#   Added INFO fields: SIFT_SCORE, SIFT_PRED (per transcript)
#
#   Step 5 — snpEff functional annotation (snpEff ann)
#   ---------------------------------------------------
#   Adds ANN fields to each variant describing the functional consequence of
#   the mutation at the transcript level: gene name, transcript ID, effect
#   (missense_variant, synonymous_variant, splice_region_variant, etc.), and
#   predicted impact (HIGH / MODERATE / LOW / MODIFIER).
#
#   Run on the SIFT-annotated VCF from Step 4 so that both SIFT scores and
#   snpEff ANN fields are present in the final output.
#
#   Flags used:
#     -v              verbose (progress and summary statistics)
#     -c <cfg>        path to snpEff.config
#     -stats <html>   write HTML summary report with variant counts by
#                     consequence type
#     -Xmx32g         Java heap size — increase if snpEff runs out of memory
#                     (keep below the SBATCH --mem allocation)
#
#   Output file: $OUTPUT_DIR/all_annotated_snpeff.vcf
#   Added INFO fields: ANN (pipe-delimited per-transcript consequence records)
#
# Output files:
#   $SIFT_RESULTS/*_SIFTannotations.vcf   — VCF with SIFT scores
#   $OUTPUT_DIR/all_annotated_snpeff.vcf  — final VCF with SIFT + ANN fields
#   $OUTPUT_DIR/snpeff_summary.html       — snpEff HTML summary report
#
# Run order:
#   script_3a_build_snpeff_db.sh   ← build snpEff database (run once)
#   script_3b_build_regions.sh     ← build SIFT4G .regions files (run once)
#   script_3c_annotate.sh          ← this script (run per dataset)
# ==============================================================================

# ------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------
# Derived paths
# ------------------------------------------------------------------------------
SIFT_DB="$SNPEFF_DATA_DIR/$GENOME_NAME"

# Input VCF — produced by script_1 (lifted + sorted)
SORTED_VCF="$OUTPUT_DIR/all_renamed.genome_coords.sorted.vcf.gz"

# Temporary uncompressed VCF for SIFT4G_Annotator (removed after Step 4)
SORTED_VCF_PLAIN="$OUTPUT_DIR/all_renamed.genome_coords.sorted.vcf"

SIFT_RESULTS="$OUTPUT_DIR/sift_results"
SNPEFF_VCF="$OUTPUT_DIR/all_annotated_snpeff.vcf"
SNPEFF_STATS="$OUTPUT_DIR/snpeff_summary.html"

# ------------------------------------------------------------------------------
# Validate prerequisites
# ------------------------------------------------------------------------------
missing=0
for f in "$SNPEFF_JAR" "$SNPEFF_CFG" "$SIFT4G_JAR" "$SORTED_VCF"; do
    if [ ! -f "$f" ]; then
        echo "[!] Missing required file: $f" >&2
        missing=1
    fi
done
if [ "$missing" -eq 1 ]; then
    echo "[!] One or more required files are missing. Aborting." >&2
    exit 1
fi

if [ ! -f "$SIFT_DB/snpEffectPredictor.bin" ]; then
    echo "[!] snpEff database not found. Run script_3a_build_snpeff_db.sh first." >&2
    exit 1
fi

REGIONS_COUNT=$(find "$SIFT_DB" -name "*.SIFTprediction" 2>/dev/null | wc -l)
if [ "$REGIONS_COUNT" -gt 0 ]; then
    echo "[3b] SIFT4G .SIFTprediction files already exist ($REGIONS_COUNT found) — skipping scorer."
    exit 0
fi

# ==============================================================================
# Step 4 — SIFT4G annotation
# ==============================================================================
echo ""
echo "[4/5] Running SIFT4G annotation..."
mkdir -p "$SIFT_RESULTS"

# Decompress the bgzipped VCF — SIFT4G_Annotator requires plain .vcf input
echo "    Decompressing VCF for SIFT4G_Annotator..."
bcftools view -O v -o "$SORTED_VCF_PLAIN" "$SORTED_VCF"
if [ $? -ne 0 ]; then
    echo "[!] Failed to decompress VCF. Aborting." >&2
    exit 1
fi

java -jar "$SIFT4G_JAR" \
    -c \
    -i "$SORTED_VCF_PLAIN" \
    -d "$SIFT_DB" \
    -r "$SIFT_RESULTS" \
    -t
SIFT4G_EXIT=$?

# Remove the uncompressed copy regardless of outcome to save space
rm -f "$SORTED_VCF_PLAIN"

if [ $SIFT4G_EXIT -ne 0 ]; then
    echo "[!] SIFT4G annotation failed (exit code $SIFT4G_EXIT)." >&2
    exit 1
fi

# SIFT4G_Annotator writes the output with the suffix _SIFTannotations.vcf
SIFT_ANNOTATED_VCF=$(find "$SIFT_RESULTS" -name "*_SIFTannotations.vcf" | head -1)
if [ -z "$SIFT_ANNOTATED_VCF" ]; then
    echo "[!] Could not find SIFT-annotated VCF in $SIFT_RESULTS" >&2
    echo "    Contents of results directory:" >&2
    ls "$SIFT_RESULTS" >&2
    exit 1
fi

echo "    SIFT-annotated VCF: $SIFT_ANNOTATED_VCF"

# ==============================================================================
# Step 5 — snpEff functional annotation
# ==============================================================================
echo ""
echo "[5/5] Running snpEff functional annotation..."

java -Xmx32g -jar "$SNPEFF_JAR" ann \
    -v \
    -c "$SNPEFF_CFG" \
    -stats "$SNPEFF_STATS" \
    "$GENOME_NAME" \
    "$SIFT_ANNOTATED_VCF" \
    > "$SNPEFF_VCF"

if [ $? -ne 0 ]; then
    echo "[!] snpEff annotation failed." >&2
    exit 1
fi


# Line near the end — the success check
REGIONS_COUNT=$(find "$SIFT_DB" -name "*.SIFTprediction" 2>/dev/null | wc -l)
if [ "$REGIONS_COUNT" -eq 0 ]; then
    echo "[!] sift4g completed but no .SIFTprediction files were written to $SIFT_DB." >&2
    exit 1
fi


# ==============================================================================
# Summary
# ==============================================================================
echo ""
echo "[+] Annotation pipeline complete!"
echo ""
echo "    Key output files:"
echo ""
echo "      $SORTED_VCF"
echo "        Lifted + sorted VCF in genome coordinates (from script_1)"
echo ""
echo "      $SIFT_ANNOTATED_VCF"
echo "        VCF with SIFT scores (SIFT_SCORE, SIFT_PRED fields)"
echo ""
echo "      $SNPEFF_VCF"
echo "        Final VCF with both SIFT scores and snpEff ANN fields"
echo ""
echo "      $SNPEFF_STATS"
echo "        snpEff HTML summary report (variant counts by consequence)"
