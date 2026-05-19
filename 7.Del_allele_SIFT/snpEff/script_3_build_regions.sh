#!/bin/bash
#SBATCH --job-name="SIFT4G_regions"
#SBATCH --mem=40G
#SBATCH --partition=medium
#SBATCH --export=ALL

# ==============================================================================
# Script 3b — Build SIFT4G per-scaffold .regions files
#
# Usage:
#   sbatch script_3b_build_regions.sh [output_dir]
#   OUTPUT_DIR defaults to current directory if not supplied as $1.
#
# Prerequisite:
#   script_3a_build_snpeff_db.sh must have completed successfully.
#   The snpEff database ($SIFT_DB/snpEffectPredictor.bin) must exist.
#
# What this script does:
#   This step is only needed ONCE per genome. If any .regions files already
#   exist in $SIFT_DB the script exits immediately.
#
#   a) Extracts protein sequences from the reference genome and GFF using
#      gffread. The output (Bmas_proteins.fa) is used as the query for sift4g.
#      A brief sanity check is printed: sequence count and the first few lines,
#      which should begin with M (methionine) for valid protein sequences.
#      This step is also skipped if Bmas_proteins.fa already exists.
#
#   b) Decompresses the UniRef90 database to a plain FASTA if not already done.
#      sift4g does not support gzip-compressed input and validates the file
#      path before opening it, so named pipes and process substitutions are
#      not viable workarounds. The uncompressed file is large (~100 GB) so
#      ensure sufficient scratch space before running:
#          df -h /home/tmichel/scratch/
#
#   c) Runs `sift4g` to align each protein against UniRef90 and write one
#      .regions file per scaffold into $SIFT_DB. These files encode the
#      per-position substitution tolerance scores that SIFT4G_Annotator.jar
#      looks up in script_3c.
#
#      Flags used:
#        -q <fasta>   query protein sequences (extracted in step a)
#        -d <fasta>   UniRef90 database (uncompressed)
#        --out <dir>  output directory — must match $SIFT_DB used by annotator
#        -t <n>       number of threads
#
# Output:
#   $OUTPUT_DIR/Bmas_proteins.fa          — extracted protein sequences
#   $SIFT_DB/*.regions                    — per-scaffold SIFT4G score files
#                                           (consumed by script_3c)
#
# Run order:
#   script_3a_build_snpeff_db.sh   ← build snpEff database (run once)
#   script_3b_build_regions.sh     ← this script (run once)
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

# Compressed UniRef90 database (as stored on shared storage)
UNIREF_DB="/mnt/shared/datasets/databases/uniprot/uniref/uniref90/uniref90.fasta.gz"

# Destination for the decompressed UniRef90 — needs ~100 GB free in scratch
UNIREF_UNZIPPED="/home/tmichel/scratch/Deleterious_alleles_PNG_baits/7.Del_allele_SIFT/uniref90.fasta"

# snpEff installation directory
SNPEFF_DIR="$HOME/software/snpEff"
GENOME_NAME="Bmas"
SNPEFF_DATA_DIR="$SNPEFF_DIR/data"

# Number of threads for sift4g
SIFT4G_THREADS=8

# ------------------------------------------------------------------------------
# Output directory: use $1 if supplied, otherwise current directory
# ------------------------------------------------------------------------------
OUTPUT_DIR="${1:-.}"
mkdir -p "$OUTPUT_DIR"

# ------------------------------------------------------------------------------
# Derived paths
# ------------------------------------------------------------------------------
SIFT_DB="$SNPEFF_DATA_DIR/$GENOME_NAME"
PROTEINS_FA="$OUTPUT_DIR/Bmas_proteins.fa"

# ------------------------------------------------------------------------------
# Validate prerequisites
# ------------------------------------------------------------------------------
if [ ! -f "$SIFT_DB/snpEffectPredictor.bin" ]; then
    echo "[!] snpEff database not found at $SIFT_DB/snpEffectPredictor.bin" >&2
    echo "    Run script_3a_build_snpeff_db.sh first." >&2
    exit 1
fi

missing=0
for f in "$BMAS_FASTA" "$BMAS_GFF"; do
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
# Step 3b — Build SIFT4G .regions files
# ==============================================================================
echo ""
REGIONS_COUNT=$(find "$SIFT_DB" -name "*.regions" 2>/dev/null | wc -l)
if [ "$REGIONS_COUNT" -gt 0 ]; then
    echo "[3b] SIFT4G .regions files already exist ($REGIONS_COUNT found) — skipping scorer."
    exit 0
fi

echo "[3b] Running SIFT4G scorer — this may take several hours..."

# ------------------------------------------------------------------------------
# a) Extract protein sequences with gffread
# ------------------------------------------------------------------------------
if [ ! -f "$PROTEINS_FA" ]; then
    echo "    Extracting protein sequences with gffread..."
    gffread "$BMAS_GFF" \
        -g "$BMAS_FASTA" \
        -y "$PROTEINS_FA"

    if [ $? -ne 0 ] || [ ! -f "$PROTEINS_FA" ]; then
        echo "[!] gffread failed to extract protein sequences." >&2
        exit 1
    fi
else
    echo "    Protein FASTA already exists — skipping gffread extraction."
fi

# Sanity check: print sequence count and first few lines.
# Valid protein sequences should begin with M (methionine).
echo "    Protein sequence count:"
grep -c '^>' "$PROTEINS_FA"
echo "    First 4 lines of $PROTEINS_FA:"
head -4 "$PROTEINS_FA"

# ------------------------------------------------------------------------------
# b) Decompress UniRef90 if not already done
#    sift4g requires a plain (uncompressed) FASTA — it validates the file path
#    before opening it and rejects pipes, FIFOs, and process substitutions.
#    The uncompressed database is ~100 GB; check available space first:
#        df -h $(dirname "$UNIREF_UNZIPPED")
# ------------------------------------------------------------------------------
if [ ! -f "$UNIREF_UNZIPPED" ]; then
    if [ ! -f "$UNIREF_DB" ]; then
        echo "[!] UniRef90 database not found at $UNIREF_DB" >&2
        echo "    Either update UNIREF_DB above or download with:" >&2
        echo "    wget https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz" >&2
        exit 1
    fi
    echo "    Decompressing UniRef90 — this will take a while and requires ~100 GB..."
    gunzip -c "$UNIREF_DB" > "$UNIREF_UNZIPPED"
    if [ $? -ne 0 ]; then
        echo "[!] Failed to decompress UniRef90 database." >&2
        exit 1
    fi
    echo "    Decompression complete: $UNIREF_UNZIPPED"
else
    echo "    Uncompressed UniRef90 already exists — skipping decompression."
fi

# ------------------------------------------------------------------------------
# c) Run sift4g scorer
#    Writes one .regions file per scaffold into $SIFT_DB.
#    These files are what SIFT4G_Annotator.jar (script_3c) looks up per variant.
# ------------------------------------------------------------------------------
echo "    Running sift4g..."
sift4g \
    -q "$PROTEINS_FA" \
    -d "$UNIREF_UNZIPPED" \
    --out "$SIFT_DB" \
    -t "$SIFT4G_THREADS"
SIFT4G_EXIT=$?

if [ $SIFT4G_EXIT -ne 0 ]; then
    echo "[!] SIFT4G scorer failed (exit code $SIFT4G_EXIT)." >&2
    exit 1
fi

REGIONS_COUNT=$(find "$SIFT_DB" -name "*.regions" 2>/dev/null | wc -l)
if [ "$REGIONS_COUNT" -eq 0 ]; then
    echo "[!] sift4g completed but no .regions files were written to $SIFT_DB." >&2
    exit 1
fi

echo ""
echo "[+] SIFT4G scorer complete."
echo "    $REGIONS_COUNT .regions files written to $SIFT_DB"
