#!/bin/bash
#SBATCH --job-name="SIFT_annotation"
#SBATCH --mem=32G
#SBATCH --partition=medium
#SBATCH --export=ALL

# ==============================================================================
# SIFT4G + snpEff annotation pipeline
#
# Usage:
#   sbatch run_sift_pipeline.sh [output_dir]
#   OUTPUT_DIR defaults to current directory if not supplied as $1.
#
# Steps:
#   1. Liftover VCF from bait coordinates to genome coordinates
#   2. Sort, index, and validate the lifted VCF
#   3. Build snpEff database for Bmas (required once; skip if already built)
#   4. Annotate with SIFT4G (uses the snpEff database directory)
#   5. Annotate with snpEff (adds functional consequence labels: ANN field)
#
# ---- Tool clarification -----------------------------------------------------
#
#   snpEff build         — builds the species database from FASTA + GFF.
#                          The resulting data/Bmas/ directory is consumed by
#                          BOTH snpEff annotate AND SIFT4G_Annotator.jar.
#
#   SIFT4G_Annotator.jar — annotates variants with SIFT scores. It reads the
#                          same snpEff data/Bmas/ directory via -d. Flags:
#                            -c        command-line / headless mode (boolean)
#                            -i <vcf>  input VCF
#                            -d <dir>  snpEff database directory (data/Bmas)
#                            -r <dir>  output results directory
#                            -t        use multiple threads
#
#   snpEff ann           — adds ANN fields (gene name, transcript consequence,
#                          impact) to the VCF. Run after SIFT4G so both
#                          annotations are present in the final VCF.
#
# ==============================================================================

# ------------------------------------------------------------------------------
# Load conda environment
# ------------------------------------------------------------------------------
source /mnt/apps/users/tmichel/conda/etc/profile.d/conda.sh
conda activate bcftools_env

# ------------------------------------------------------------------------------
# Input paths — edit these as needed
# ------------------------------------------------------------------------------
VCF_INPUT="/home/tmichel/scratch/Deleterious_alleles_PNG_baits/6.Variant_calling_annotations/all.vcf.gz"
MAPPING_TSV="/home/tmichel/scratch/Deleterious_alleles_PNG_baits/1.Baits_to_genome/bait_locus_annotation.tsv"
BMAS_FASTA="/home/tmichel/projects/rbge/tmichel/reference_genomes/Bmas.fa"
BMAS_GFF="/home/tmichel/projects/rbge/tmichel/reference_genomes/Bmas.gff"

# snpEff installation directory (contains snpEff.jar and snpEff.config)
SNPEFF_DIR="$HOME/software/snpEff"
SNPEFF_JAR="$SNPEFF_DIR/snpEff.jar"
SNPEFF_CFG="$SNPEFF_DIR/snpEff.config"

# SIFT4G_Annotator JAR (MO-SIFT-21Apr14.jar or similar)
SIFT4G_JAR="$HOME/software/SIFT4G_Annotator/SIFT4G_Annotator.jar"

# Genome name used in snpEff — must match the entry added to snpEff.config
GENOME_NAME="Bmas"

# snpEff data directory where the Bmas database will live
# snpEff expects: $SNPEFF_DATA_DIR/$GENOME_NAME/{sequences.fa, genes.gff}
SNPEFF_DATA_DIR="$SNPEFF_DIR/data"

# # ------------------------------------------------------------------------------
# Output directory: use $1 if supplied, otherwise current directory
# ------------------------------------------------------------------------------
OUTPUT_DIR="${1:-.}"

# ---- ADD THIS BLOCK ----------------------------------------------------------
# Intermediate / output file paths  (must match script_1 exactly)
LIFTED_VCF="$OUTPUT_DIR/all_renamed.genome_coords.vcf.gz"
SORTED_VCF="$OUTPUT_DIR/all_renamed.genome_coords.sorted.vcf.gz"
FAILED_VCF="$OUTPUT_DIR/failed_liftover.vcf"
SIFT_RESULTS="$OUTPUT_DIR/sift_results"
SNPEFF_VCF="$OUTPUT_DIR/all_annotated_snpeff.vcf"
SNPEFF_STATS="$OUTPUT_DIR/snpeff_summary.html"
SIFT_DB="$SNPEFF_DATA_DIR/$GENOME_NAME"
# ------------------------------------------------------------------------------


# ==============================================================================
# Step 3 — Build snpEff database for Bmas
#
# This step is only needed once. If the database already exists (i.e.
# $SIFT_DB/snpEffectPredictor.bin is present) this step is skipped.
#
# What happens here:
#   a) The Bmas genome entry is appended to snpEff.config (if not already there)
#   b) sequences.fa and genes.gff are copied into $SNPEFF_DATA_DIR/Bmas/
#   c) snpEff build compiles the binary database
#
# The resulting $SNPEFF_DATA_DIR/Bmas/ directory is used by BOTH
# SIFT4G_Annotator (-d flag) and snpEff ann in Step 5.
# ==============================================================================
SIFT_DB="$SNPEFF_DATA_DIR/$GENOME_NAME"

echo ""
if [ -f "$SIFT_DB/snpEffectPredictor.bin" ]; then
    echo "[3/5] snpEff database already exists — skipping build."
    echo "      $SIFT_DB/snpEffectPredictor.bin"
else
    echo "[3/5] Building snpEff database for $GENOME_NAME..."

    # a) Add genome entry to snpEff.config if not already present
    if ! grep -q "^${GENOME_NAME}.genome" "$SNPEFF_CFG"; then
        echo "" >> "$SNPEFF_CFG"
        echo "# ${GENOME_NAME} — Begonia masoniana" >> "$SNPEFF_CFG"
        echo "${GENOME_NAME}.genome : ${GENOME_NAME}" >> "$SNPEFF_CFG"
        echo "    Added $GENOME_NAME entry to $SNPEFF_CFG"
    else
        echo "    $GENOME_NAME already present in $SNPEFF_CFG"
    fi

    # b) Populate the data directory
    SIFT_DB="$SNPEFF_DATA_DIR/$GENOME_NAME"
    mkdir -p "$SIFT_DB"
    cp "$BMAS_FASTA" "$SIFT_DB/sequences.fa"
    cp "$BMAS_GFF"   "$SIFT_DB/genes.gff"

    # c) Build the binary database
    #
    #    -gff3           annotation format (GLEAN output is GFF3)
    #    -noCheckCds     skip CDS sequence validation — we have no cds.fa
    #    -noCheckProtein skip protein sequence validation — we have no protein.fa
    #
    #    NOTE: snpEff exits with code 1 when cds.fa / protein.fa are absent,
    #    even though the database was written correctly. We therefore ignore
    #    the exit code and instead test whether snpEffectPredictor.bin was
    #    actually produced.
    #
    #    WARNING_GENE_NOT_FOUND is expected for GLEAN annotations, which lack
    #    explicit 'gene' feature lines. snpEff infers gene boundaries from mRNA
    #    features and the database is still valid.
    java -Xmx32g -jar "$SNPEFF_JAR" build \
        -c "$SNPEFF_CFG" \
        -gff3 \
        -noCheckCds \
        -noCheckProtein \
        -v \
        "$GENOME_NAME" \
        > "$OUTPUT_DIR/snpeff_build.stdout" \
        2> "$OUTPUT_DIR/snpeff_build.stderr"

    # Ignore exit code — snpEff returns 1 when optional check files are absent
    # even when the database was written successfully. Test the output directly.
    if [ ! -f "$SIFT_DB/snpEffectPredictor.bin" ]; then
        echo "[!] snpEff database build failed — snpEffectPredictor.bin not found." >&2
        echo "    Check: $OUTPUT_DIR/snpeff_build.stderr" >&2
        exit 1
    fi

    echo "    Database built: $SIFT_DB/snpEffectPredictor.bin"
fi

# ==============================================================================
# Step 4 — SIFT4G annotation
#
# SIFT4G_Annotator.jar reads the snpEff database directory directly via -d.
# Flags:
#   -c        headless / command-line mode (boolean, no argument)
#   -i <vcf>  input VCF
#   -d <dir>  snpEff data/Bmas directory (same one built in Step 3)
#   -r <dir>  results output directory
#   -t        multithreaded mode
# ==============================================================================
echo ""
echo "[4/5] Running SIFT4G annotation..."
mkdir -p "$SIFT_RESULTS"

java -jar "$SIFT4G_JAR" \
    -c \
    -i "$SORTED_VCF" \
    -d "$SIFT_DB" \
    -r "$SIFT_RESULTS" \
    -t

if [ $? -ne 0 ]; then
    echo "[!] SIFT4G annotation failed." >&2
    exit 1
fi

# Locate the SIFT-annotated VCF — SIFT4G_Annotator appends _SIFTannotations.vcf
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
#
# Adds ANN fields to each variant: gene name, transcript ID, consequence
# (missense, synonymous, splice, etc.), and predicted impact (HIGH/MODERATE/
# LOW/MODIFIER). Run on the SIFT-annotated VCF so both annotations are merged
# in the final output.
#
# -v            verbose (prints progress and stats)
# -c            path to snpEff.config
# -stats        write HTML summary report
# -Xmx32g       heap size — increase if snpEff runs out of memory
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

# ==============================================================================
# Summary
# ==============================================================================
echo ""
echo "[+] Annotation pipeline complete!"
echo ""
echo "    Key output files:"
echo ""
echo "      $SORTED_VCF"
echo "        Lifted + sorted VCF in genome coordinates"
echo ""
echo "      $FAILED_VCF"
echo "        Variants that could not be lifted (inspect for QC)"
echo ""
echo "      $SIFT_ANNOTATED_VCF"
echo "        VCF with SIFT scores (SIFT_SCORE, SIFT_PRED fields)"
echo ""
echo "      $SNPEFF_VCF"
echo "        Final VCF with both SIFT scores and snpEff ANN fields"
echo ""
echo "      $SNPEFF_STATS"
echo "        snpEff HTML summary report (variant counts by consequence)"
