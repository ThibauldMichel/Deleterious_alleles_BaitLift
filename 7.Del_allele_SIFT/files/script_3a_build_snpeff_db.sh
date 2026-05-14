#!/bin/bash
#SBATCH --job-name="snpEff_build"
#SBATCH --mem=40G
#SBATCH --partition=medium
#SBATCH --export=ALL

# ==============================================================================
# Script 3a — Build snpEff database for Bmas
#
# Usage:
#   sbatch script_3a_build_snpeff_db.sh [output_dir]
#   OUTPUT_DIR defaults to current directory if not supplied as $1.
#
# What this script does:
#   This step is only needed ONCE. If the binary database already exists
#   ($SIFT_DB/snpEffectPredictor.bin) the script exits immediately.
#
#   a) Appends the Bmas genome entry to snpEff.config (if not already present).
#      snpEff requires a line of the form:
#          Bmas.genome : Bmas
#      in its config file before it will build or query a custom genome.
#
#   b) Copies the reference FASTA and GFF into the snpEff data directory:
#          $SNPEFF_DATA_DIR/Bmas/sequences.fa   ← reference genome
#          $SNPEFF_DATA_DIR/Bmas/genes.gff      ← gene annotation
#      snpEff expects exactly these filenames in a subdirectory named after
#      the genome (GENOME_NAME).
#
#   c) Runs `snpEff build` to compile the binary database
#      (snpEffectPredictor.bin). This binary is consumed by both snpEff ann
#      (Step 5) and SIFT4G_Annotator.jar (Step 4).
#
#      Flags used:
#        -gff3           input annotation is GFF3 format
#        -noCheckCds     skip CDS sequence validation (no cds.fa available)
#        -noCheckProtein skip protein sequence validation (no protein.fa)
#        -v              verbose output for debugging
#
#      NOTE: snpEff exits with code 1 when cds.fa / protein.fa are absent,
#      even when the database is written correctly. The exit code is therefore
#      ignored; success is tested by checking whether snpEffectPredictor.bin
#      was actually produced.
#
#      WARNING_GENE_NOT_FOUND is expected for GLEAN annotations, which omit
#      explicit 'gene' feature lines. snpEff infers gene boundaries from mRNA
#      features and the database is still valid.
#
# Output:
#   $SNPEFF_DATA_DIR/Bmas/snpEffectPredictor.bin  — binary snpEff database
#   $OUTPUT_DIR/snpeff_build.stdout               — snpEff build stdout log
#   $OUTPUT_DIR/snpeff_build.stderr               — snpEff build stderr log
#
# Run order:
#   script_3a_build_snpeff_db.sh   ← this script (run once)
#   script_3b_build_regions.sh     ← build SIFT4G .regions files (run once)
#   script_3c_annotate.sh          ← SIFT4G + snpEff annotation (run per dataset)
# ==============================================================================

# ------------------------------------------------------------------------------
# Load conda environment
# ------------------------------------------------------------------------------
source /mnt/apps/users/tmichel/conda/etc/profile.d/conda.sh
conda activate bcftools_env

# ------------------------------------------------------------------------------
# Input paths — edit these as needed
# ------------------------------------------------------------------------------
BMAS_FASTA="/home/tmichel/projects/rbge/tmichel/reference_genomes/Bmas.fa"
BMAS_GFF="/home/tmichel/projects/rbge/tmichel/reference_genomes/Bmas.gff"

# snpEff installation directory (contains snpEff.jar and snpEff.config)
SNPEFF_DIR="$HOME/software/snpEff"
SNPEFF_JAR="$SNPEFF_DIR/snpEff.jar"
SNPEFF_CFG="$SNPEFF_DIR/snpEff.config"

# Genome name — must match the entry added to snpEff.config
GENOME_NAME="Bmas"

# snpEff data directory; snpEff expects:
#   $SNPEFF_DATA_DIR/$GENOME_NAME/sequences.fa
#   $SNPEFF_DATA_DIR/$GENOME_NAME/genes.gff
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

# ------------------------------------------------------------------------------
# Validate required files and tools
# ------------------------------------------------------------------------------
missing=0
for f in "$BMAS_FASTA" "$BMAS_GFF" "$SNPEFF_JAR" "$SNPEFF_CFG"; do
    if [ ! -f "$f" ]; then
        echo "[!] Missing required file: $f" >&2
        missing=1
    fi
done
if [ "$missing" -eq 1 ]; then
    echo "[!] One or more required files are missing. Aborting." >&2
    exit 1
fi

# ==============================================================================
# Step 3a — Build snpEff database
# ==============================================================================
echo ""
if [ -f "$SIFT_DB/snpEffectPredictor.bin" ]; then
    echo "[3a] snpEff database already exists — skipping build."
    echo "     $SIFT_DB/snpEffectPredictor.bin"
    exit 0
fi

echo "[3a] Building snpEff database for $GENOME_NAME..."

# a) Add genome entry to snpEff.config if not already present
if ! grep -q "^${GENOME_NAME}.genome" "$SNPEFF_CFG"; then
    echo "" >> "$SNPEFF_CFG"
    echo "# ${GENOME_NAME} — Begonia masoniana" >> "$SNPEFF_CFG"
    echo "${GENOME_NAME}.genome : ${GENOME_NAME}" >> "$SNPEFF_CFG"
    echo "    Added $GENOME_NAME entry to $SNPEFF_CFG"
else
    echo "    $GENOME_NAME already present in $SNPEFF_CFG"
fi

# b) Populate the data directory with the required filenames
mkdir -p "$SIFT_DB"
cp "$BMAS_FASTA" "$SIFT_DB/sequences.fa"
cp "$BMAS_GFF"   "$SIFT_DB/genes.gff"
echo "    Copied FASTA and GFF into $SIFT_DB"

# c) Build the binary database
#    Exit code is intentionally ignored — see header note above.
java -Xmx32g -jar "$SNPEFF_JAR" build \
    -c "$SNPEFF_CFG" \
    -gff3 \
    -noCheckCds \
    -noCheckProtein \
    -v \
    "$GENOME_NAME" \
    > "$OUTPUT_DIR/snpeff_build.stdout" \
    2> "$OUTPUT_DIR/snpeff_build.stderr"

# Test for the output binary directly rather than trusting the exit code
if [ ! -f "$SIFT_DB/snpEffectPredictor.bin" ]; then
    echo "[!] snpEff database build failed — snpEffectPredictor.bin not found." >&2
    echo "    Check: $OUTPUT_DIR/snpeff_build.stderr" >&2
    exit 1
fi

echo ""
echo "[+] snpEff database built successfully."
echo "    $SIFT_DB/snpEffectPredictor.bin"
echo "    Build logs: $OUTPUT_DIR/snpeff_build.stdout / .stderr"
