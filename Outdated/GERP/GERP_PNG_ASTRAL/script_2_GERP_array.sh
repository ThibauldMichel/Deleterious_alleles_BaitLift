#!/bin/bash
#SBATCH --job-name="GERP++"
#SBATCH --output=GERP_%A_%a.out
#SBATCH --mem=8G
#SBATCH --partition=medium
#SBATCH --export=ALL



# Paths
pathgerp="$HOME/GERPplusplus"
pathtree="/home/tmichel/scratch/Deleterious_alleles_PNG/CAPTUS/05_phylogeny_FastTreeAstral/species_tree_rescaled.treefile"
pathrefe="$HOME/projects/rbge/tmichel/reference_genomes/Hannah_Begonia_baits.fasta"

# Get file name corresponding to this array task
file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" mfa_files.txt)
echo "Processing file: $file"

# Clean header
sed -E 's/^(>[^ ]+).*/\1/' "mfa_dir/${file}.mfa" > "mfa_dir/${file}_edited.mfa"

# Run gerpcol
"$pathgerp"/gerpcol -a -v -e Hannah_Begonia_baits__ref \
    -f "mfa_dir/${file}_edited.mfa" \
    -t "$pathtree" \
    -x ".rates"

# Run gerpelem (only if not done)
infile="mfa_dir/${file}_edited.mfa.rates"
outfile="${infile}.elems"

if [ ! -f "$outfile" ]; then
    echo "Running gerpelem on $infile"
    "$pathgerp"/gerpelem -f "$infile"
else
    echo "Skipping $file, $outfile exists."
fi

