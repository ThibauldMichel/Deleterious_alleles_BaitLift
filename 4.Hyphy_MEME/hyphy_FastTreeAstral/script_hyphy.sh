#!/bin/bash

#SBATCH --job-name="Hyphy_PNG"
#SBATCH --export=ALL
#SBATCH --mem=64G
#SBATCH --partition=long

# Load conda into this shell session
source /mnt/apps/users/tmichel/conda/etc/profile.d/conda.sh

conda activate hyphy_env



# Make sure hyphy does not stop from previous data
#rm -fr raw_msa
#rm -fr cleaned_msa
#rm -fr output_meme

# From https://github.com/veg/hyphy


#hyphy <method_name> --alignment <path_to_alignment_file> 

# Path to the tree file
#pathtree="$HOME/scratch/CAPTUS/CAPTUS-PNG/05_phylogeny/concat.treefile"
pathtree="/home/tmichel/scratch/Deleterious_alleles_PNG/CAPTUS/CAPTUS-PNG/05_phylogeny_FastTreeAstral/ASTRAL/species_tree_astral.treefile"

# Path to the alignment files
pathalig="$HOME/scratch/Deleterious_alleles_PNG/CAPTUS/CAPTUS-PNG/04_alignments/03_trimmed/06_informed/01_coding_NUC/02_NT"



# Create a list of the mfa files without path and extension.
#for file in "$pathalig"/*.fna; do basename  "$file" .fna; done > mfa_files.txt

# Create a directory for msa with clean headers
#mkdir raw_msa

# Remove the extra square brackets from the header of the fna files.
#cat mfa_files.txt | while read file;
#do
#sed -E 's/^(>[^ ]+).*/\1/' "$pathalig/${file}.fna"  > "raw_msa/${file}_edited.fna";
#done

# Create a folder to store files without STOP codons
#mkdir cleaned_msa

# Remove the STOP codons from alignments with an unnecessary complicated command line
#while IFS= read -r file; do
#    hyphy rmv Universal "raw_msa/${file}_edited.fna" "Disallow stops" "./cleaned_msa/${file}_NOSTOP.fna"
#done < mfa_files.txt


# put outputs of MEME analysis
#mkdir output_meme

# Make the analysis of the cleaned alignments
#while read file;
#do
#hyphy meme --alignment ./cleaned_msa/"$file"_NOSTOP.fna --tree "$pathtree" --branches "All" --pvalue 0.1 --output ./output_meme/ --full-model Yes ;
#    done < mfa_files.txt



#hyphy meme --alignment ./cleaned_msa/ACmerged_contig_1011_NOSTOP.fna --tree "$pathtree" --branches "All" --pvalue 0.1 --output ./output_meme/ACmerged_contig_1011_MEME.fna --full-model Yes


# Remake the mfa_files.txt for synced data
for file in synced_msa/*.fna; do basename  "$file" _synced.fna; done > mfa_files.txt


# Loop over loci
while read -r locus; do
  pathtree="synced_trees/${locus}_synced.treefile"
  pathmsa="synced_msa/${locus}_synced.fna"
  pathout="output_meme/${locus}.json"
  logdir="logs"


  # Skip if output file already exists
  if [[ -f "$pathout" ]]; then
    echo "Skipping $locus, already processed."
    continue
  fi

  echo "Processing $locus..."
  hyphy meme \
    --alignment "$pathmsa" \
    --tree "$pathtree" \
    --branches All \
    --pvalue 0.1 \
    --output "$pathout" \
    --full-model Yes \
    > "${logdir}/${locus}.log" 2>&1

done < mfa_files.txt


#hyphy meme --help

#Available analysis command line options
#---------------------------------------
#Use --option VALUE syntax to invoke
#If a [reqired] option is not provided on the command line, the analysis will #prompt for its value
#[conditionally required] options may or not be required based on the values of other options

#code
#	Which genetic code should be used
#	default value: Universal

#alignment [required]
#	An in-frame codon alignment in one of the formats supported by HyPhy

#tree [conditionally required]
#	A phylogenetic tree (optionally annotated with {})
#	applies to: Please select a tree file for the data:

#branches
#	Branches to test
#	default value: All

#pvalue
#	The p-value threshold to use when testing for selection
#	default value: 0.1

#output
#	Write the resulting JSON to this file (default is to save to the same path as the alignment file + 'MEME.json')
#	default value: meme.codon_data_info[terms.json.json] [computed at run time]

#full-model
#	Perform branch length re-optimization under the full codon model
#	default value: Yes




