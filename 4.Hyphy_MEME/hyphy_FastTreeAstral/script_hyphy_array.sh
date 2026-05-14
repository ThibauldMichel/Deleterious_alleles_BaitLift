#!/bin/bash
#SBATCH --job-name="Hyphy_PNG"
#SBATCH --export=ALL
#SBATCH --mem=64G
#SBATCH --partition=medium
#SBATCH --array=1-1239%10   # 10 jobs in parallel max
#SBATCH --output=logs/%A_%a.out
#SBATCH --error=logs/%A_%a.err

# Load conda
source /mnt/apps/users/tmichel/conda/etc/profile.d/conda.sh
conda activate hyphy_env

# Path setup
pathtree="/home/tmichel/scratch/Deleterious_alleles_PNG/CAPTUS/CAPTUS-PNG/05_phylogeny_FastTreeAstral/ASTRAL/species_tree_astral.treefile"

# Make sure logs and output directories exist
mkdir -p logs output_meme synced_trees synced_msa

# Generate mfa_files.txt (if not already)
if [ ! -f mfa_files.txt ]; then
    for file in synced_msa/*.fna; do
        basename "$file" _synced.fna
    done > mfa_files.txt
fi

# Get locus name for this array task
locus=$(sed -n "${SLURM_ARRAY_TASK_ID}p" mfa_files.txt)

# Define paths for this locus
pathtree="synced_trees/${locus}_synced.treefile"
pathmsa="synced_msa/${locus}_synced.fna"
pathout="output_meme/${locus}.json"
logfile="logs/${locus}.log"

# Skip if already processed
if [[ -f "$pathout" ]]; then
    echo "Skipping $locus, already processed."
    exit 0
fi

echo "[$(date)] Processing $locus..."

hyphy meme \
    --alignment "$pathmsa" \
    --tree "$pathtree" \
    --branches All \
    --pvalue 0.1 \
    --output "$pathout" \
    --full-model Yes \
    > "$logfile" 2>&1

echo "[$(date)] Finished $locus."




