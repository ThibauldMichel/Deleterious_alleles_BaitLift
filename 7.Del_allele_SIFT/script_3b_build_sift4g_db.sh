#!/bin/bash
#SBATCH --job-name="SIFT4G_build_db"
#SBATCH --mem=40G
#SBATCH --partition=medium
#SBATCH --export=ALL

# ==============================================================================
# Script 3b — Build SIFT4G genomic database using SIFT4G_Create_Genomic_DB
#
# Usage:
#   sbatch script_3b_build_sift4g_db.sh [output_dir]
#   OUTPUT_DIR defaults to current directory if not supplied as $1.
#
# Prerequisite:
#   script_2_build_snpeff_db.sh must have completed successfully.
#
# What this script does:
#   This step is only needed ONCE per genome. It replaces the manual sift4g
#   scorer call in the original script_3_build_regions.sh, which produced
#   per-gene .SIFTprediction files but NOT the scaffold-level .regions and
#   .gz files that SIFT4G_Annotator.jar requires.
#
#   This script uses the SIFT4G_Create_Genomic_DB pipeline (pauline-ng,
#   https://github.com/pauline-ng/SIFT4G_Create_Genomic_DB) which handles
#   the full database build: protein extraction, sift4g scoring, and
#   coordinate-indexed .regions + .gz file generation.
#
#   Steps performed:
#   a) Converts the GFF annotation to GTF format using gffread, then
#      compresses both the FASTA and GTF with gzip, as SIFT4G_Create_Genomic_DB
#      expects gzip-compressed inputs.
#
#   b) Builds the required directory structure:
#        <SIFT4G_DB_PARENT>/chr-src/          ← genome FASTA (.fa.gz)
#        <SIFT4G_DB_PARENT>/gene-annotation-src/ ← annotation (.gtf.gz)
#
#   c) Writes a config file for make-SIFT-db-all.pl, setting all required
#      parameters for a non-Ensembl local-file build.
#
#   d) Runs make-SIFT-db-all.pl, which:
#        - Extracts protein sequences per transcript
#        - Runs sift4g against UniRef90 to score each protein
#        - Produces per-scaffold .regions and .gz files in
#          <SIFT4G_DB_PARENT>/<ORG_VERSION>/
#      This is the directory that must be passed to SIFT4G_Annotator.jar
#      via the -d flag in script_4_SIFT_annotate.sh.
#
#   NOTE: This step can take many hours depending on genome size and the
#   number of annotated transcripts. The UniRef90 database must be
#   decompressed (~100 GB) before running — the script checks for this
#   and will decompress it if needed.
#
# Output:
#   <SIFT4G_DB_PARENT>/<ORG_VERSION>/         ← SIFT4G database directory
#     scaffold1.gz, scaffold1.regions, ...    ← per-scaffold score files
#     CHECK_GENES.LOG                         ← summary of SIFT coverage
#
# Run order:
#   script_1_lift_VCF.sh              ← liftover (run per dataset)
#   script_2_build_snpeff_db.sh       ← build snpEff database (run once)
#   script_3b_build_sift4g_db.sh      ← this script (run once)
#   script_3_build_regions.sh         ← now superseded by this script
#   script_4_SIFT_annotate.sh         ← annotation (run per dataset)
#                                        NOTE: update -d path (see below)
#
# IMPORTANT — update script_4_SIFT_annotate.sh:
#   The -d flag passed to SIFT4G_Annotator.jar must point to the output of
#   THIS script, not the snpEff data/Bmas directory. Change SIFT_DB in
#   script_4 to:
#     SIFT_DB="$SIFT4G_DB_PARENT/$ORG_VERSION"
# ==============================================================================

# ------------------------------------------------------------------------------
# Load conda environment
# ------------------------------------------------------------------------------
source /mnt/apps/users/tmichel/conda/etc/profile.d/conda.sh
conda activate perl_env

# ------------------------------------------------------------------------------
# Input paths — edit these as needed
# ------------------------------------------------------------------------------
BMAS_FASTA="/home/tmichel/projects/rbge/tmichel/reference_genomes/Bmas.fa"
BMAS_GFF="/home/tmichel/projects/rbge/tmichel/reference_genomes/Bmas.gff"

# Compressed UniRef90 database
UNIREF_DB="/mnt/shared/datasets/databases/uniprot/uniref/uniref90/uniref90.fasta.gz"

# Destination for decompressed UniRef90 — requires ~100 GB free
UNIREF_UNZIPPED="/home/tmichel/scratch/Deleterious_alleles_PNG_baits/7.Del_allele_SIFT/uniref90.fasta"

# Path to the SIFT4G_Create_Genomic_DB scripts directory
# Clone with: git clone https://github.com/pauline-ng/SIFT4G_Create_Genomic_DB.git
SIFT4G_CREATE_DB_DIR="$HOME/software/SIFT4G_Create_Genomic_DB"

# Path to the sift4g executable
SIFT4G_BIN="$(which sift4g)"

# Organism name and version — used to name the output database directory
ORG="Bmas"
ORG_VERSION="Bmas_v1"

# Parent directory for the SIFT4G database build
# Final database will be at: $SIFT4G_DB_PARENT/$ORG_VERSION/
SIFT4G_DB_PARENT="/home/tmichel/scratch/Deleterious_alleles_PNG_baits/7.Del_allele_SIFT/sift4g_db"

# Genetic code table (1 = Standard; plants typically use 1)
GENETIC_CODE_TABLE=1
GENETIC_CODE_TABLENAME="Standard"

# Mitochondrial genetic code (for sequences named Mt/chrM/Mito)
# Set to 1 to use Standard for all — change if Bmas has annotated mito genes
MITO_GENETIC_CODE_TABLE=1
MITO_GENETIC_CODE_TABLENAME="Standard"

# Number of threads for sift4g
SIFT4G_THREADS=8

# ------------------------------------------------------------------------------
# Output directory: use $1 if supplied, otherwise current directory
# ------------------------------------------------------------------------------
OUTPUT_DIR="${1:-.}"
mkdir -p "$OUTPUT_DIR"

# ------------------------------------------------------------------------------
# Validate prerequisites
# ------------------------------------------------------------------------------
missing=0
for f in "$BMAS_FASTA" "$BMAS_GFF"; do
    if [ ! -f "$f" ]; then
        echo "[!] Missing required file: $f" >&2
        missing=1
    fi
done
if [ ! -d "$SIFT4G_CREATE_DB_DIR" ]; then
    echo "[!] SIFT4G_Create_Genomic_DB not found at $SIFT4G_CREATE_DB_DIR" >&2
    echo "    Clone with: git clone https://github.com/pauline-ng/SIFT4G_Create_Genomic_DB.git $SIFT4G_CREATE_DB_DIR" >&2
    missing=1
fi
if [ ! -f "$SIFT4G_CREATE_DB_DIR/make-SIFT-db-all.pl" ]; then
    echo "[!] make-SIFT-db-all.pl not found in $SIFT4G_CREATE_DB_DIR" >&2
    missing=1
fi
if [ -z "$SIFT4G_BIN" ] || [ ! -x "$SIFT4G_BIN" ]; then
    echo "[!] sift4g executable not found in PATH." >&2
    missing=1
fi
if [ "$missing" -eq 1 ]; then
    echo "[!] One or more prerequisites are missing. Aborting." >&2
    exit 1
fi

# Check if the database already exists
FINAL_DB="$SIFT4G_DB_PARENT/$ORG_VERSION"
if [ -d "$FINAL_DB" ] && [ "$(find "$FINAL_DB" -name "*.regions" 2>/dev/null | wc -l)" -gt 0 ]; then
    echo "[3b] SIFT4G database already exists at $FINAL_DB — skipping build."
    REGIONS_COUNT=$(find "$FINAL_DB" -name "*.regions" | wc -l)
    echo "     Found $REGIONS_COUNT .regions files."
    exit 0
fi

echo ""
echo "[3b] Building SIFT4G genomic database — this will take several hours..."

# ==============================================================================
# a) Convert GFF to GTF and compress inputs
# ==============================================================================
echo ""
echo "    [3b-a] Converting GFF to GTF format..."

WORK_DIR="$OUTPUT_DIR/sift4g_build_work"
mkdir -p "$WORK_DIR"

GTF_FILE="$WORK_DIR/${ORG}.gene.gtf"

gffread "$BMAS_GFF" -T -o "$GTF_FILE"
if [ $? -ne 0 ] || [ ! -f "$GTF_FILE" ]; then
    echo "[!] gffread GFF->GTF conversion failed." >&2
    exit 1
fi
echo "    GTF written to $GTF_FILE"

# Compress inputs — SIFT4G_Create_Genomic_DB expects .fa.gz and .gtf.gz
echo "    Compressing FASTA and GTF..."
FASTA_GZ="$WORK_DIR/${ORG}.fa.gz"
GTF_GZ="$WORK_DIR/${ORG}.gene.gtf.gz"

gzip -c "$BMAS_FASTA" > "$FASTA_GZ"
gzip -c "$GTF_FILE"   > "$GTF_GZ"

if [ ! -f "$FASTA_GZ" ] || [ ! -f "$GTF_GZ" ]; then
    echo "[!] Failed to compress FASTA or GTF." >&2
    exit 1
fi

# ==============================================================================
# b) Build the SIFT4G_Create_Genomic_DB directory structure
# ==============================================================================
echo ""
echo "    [3b-b] Setting up SIFT4G_Create_Genomic_DB directory structure..."

CHR_SRC="$SIFT4G_DB_PARENT/chr-src"
GENE_SRC="$SIFT4G_DB_PARENT/gene-annotation-src"

mkdir -p "$CHR_SRC" "$GENE_SRC" "$SIFT4G_DB_PARENT/dbSNP"

cp "$FASTA_GZ" "$CHR_SRC/"
cp "$GTF_GZ"   "$GENE_SRC/"

echo "    Files placed:"
echo "      $CHR_SRC/$(basename $FASTA_GZ)"
echo "      $GENE_SRC/$(basename $GTF_GZ)"

# ==============================================================================
# c) Write the config file for make-SIFT-db-all.pl
# ==============================================================================
echo ""
echo "    [3b-c] Writing config file..."

CONFIG_FILE="/home/tmichel/scratch/Deleterious_alleles_PNG_baits/7.Del_allele_SIFT/sift4g_${ORG}.config.txt"

cat > "$CONFIG_FILE" << EOF
SIFT4G_PATH=$SIFT4G_BIN
PROTEIN_DB=$UNIREF_UNZIPPED
PARENT_DIR=$SIFT4G_DB_PARENT
ORG=$ORG
ORG_VERSION=$ORG_VERSION
GENETIC_CODE_TABLE=$GENETIC_CODE_TABLE
GENETIC_CODE_TABLENAME=$GENETIC_CODE_TABLENAME
MITO_GENETIC_CODE_TABLE=$MITO_GENETIC_CODE_TABLE
MITO_GENETIC_CODE_TABLENAME=$MITO_GENETIC_CODE_TABLENAME
EOF

echo "    Config written to $CONFIG_FILE"
echo ""
echo "    Config contents:"
sed 's/^/      /' "$CONFIG_FILE"

# ==============================================================================
# d) Decompress UniRef90 if needed
#    make-SIFT-db-all.pl invokes sift4g which requires plain FASTA input.
# ==============================================================================
echo ""
if [ ! -f "$UNIREF_UNZIPPED" ]; then
    if [ ! -f "$UNIREF_DB" ]; then
        echo "[!] UniRef90 database not found at $UNIREF_DB" >&2
        echo "    Either update UNIREF_DB above or download with:" >&2
        echo "    wget https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz" >&2
        exit 1
    fi
    echo "    [3b-d] Decompressing UniRef90 — requires ~100 GB and will take a while..."
    echo "    Checking available space:"
    df -h "$(dirname "$UNIREF_UNZIPPED")"
    gunzip -c "$UNIREF_DB" > "$UNIREF_UNZIPPED"
    if [ $? -ne 0 ]; then
        echo "[!] Failed to decompress UniRef90." >&2
        exit 1
    fi
    echo "    Decompression complete: $UNIREF_UNZIPPED"
else
    echo "    [3b-d] Uncompressed UniRef90 already exists — skipping decompression."
fi

# ==============================================================================
# e) Run make-SIFT-db-all.pl
# ==============================================================================
echo ""
echo "    [3b-e] Running make-SIFT-db-all.pl..."
echo "    This can take 1-24+ hours depending on genome size."
echo "    Monitor progress with:"
echo "      ls -lt $SIFT4G_DB_PARENT/singleRecords/"
echo "      ls $SIFT4G_DB_PARENT/SIFT_predictions/"

cd "$SIFT4G_CREATE_DB_DIR" || exit 1

perl make-SIFT-db-all.pl -config "$CONFIG_FILE"
MAKE_EXIT=$?

# The script prints "All done!" on success; also verify .regions files exist
REGIONS_COUNT=$(find "$FINAL_DB" -name "*.regions" 2>/dev/null | wc -l)

if [ $MAKE_EXIT -ne 0 ] || [ "$REGIONS_COUNT" -eq 0 ]; then
    echo "[!] make-SIFT-db-all.pl failed or produced no .regions files." >&2
    echo "    Exit code: $MAKE_EXIT" >&2
    echo "    .regions files found: $REGIONS_COUNT" >&2
    echo "    Check the contents of $SIFT4G_DB_PARENT for partial output." >&2
    exit 1
fi

# ==============================================================================
# Summary
# ==============================================================================
echo ""
echo "[+] SIFT4G genomic database built successfully."
echo ""
echo "    Database directory (use this as -d in script_4_SIFT_annotate.sh):"
echo "      $FINAL_DB"
echo ""
echo "    $REGIONS_COUNT .regions files produced."
echo ""
echo "    Database coverage summary:"
if [ -f "$FINAL_DB/CHECK_GENES.LOG" ]; then
    tail -5 "$FINAL_DB/CHECK_GENES.LOG"
else
    echo "    (CHECK_GENES.LOG not found — check $SIFT4G_DB_PARENT)"
fi
echo ""
echo "    IMPORTANT: Update SIFT_DB in script_4_SIFT_annotate.sh to:"
echo "      SIFT_DB=\"$FINAL_DB\""
