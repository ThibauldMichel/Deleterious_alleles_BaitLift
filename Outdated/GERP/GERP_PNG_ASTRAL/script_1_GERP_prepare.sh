#!/bin/bash
#SBATCH --job-name="GERP_prepare"
#SBATCH --mem=2G

# Define paths
pathalig="$HOME/scratch/Deleterious_alleles_PNG/CAPTUS/04_alignments/03_trimmed/03_informed_w_refs/01_coding_NUC/02_NT"

# Ensure mfa_dir exists and copy files
rm -rf mfa_dir
mkdir -p mfa_dir

for file in "$pathalig"/*.fna; do
    filename=$(basename "$file")
    rsync -az "$file" "mfa_dir/${filename%.fna}.mfa"
done

# Make list of mfa files
for file in mfa_dir/*.mfa; do basename "$file" .mfa; done > mfa_files.txt

# Submit job array (number of lines = number of files)
nfiles=$(wc -l < mfa_files.txt)
sbatch --array=1-${nfiles}%10 script_2_GERP_array.sh

