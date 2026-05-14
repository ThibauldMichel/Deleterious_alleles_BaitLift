#!/bin/bash
#SBATCH --job-name="missing_baits"
#SBATCH --export=ALL
#SBATCH --mem=8G  
#SBATCH --partition=medium



BAITS=/home/tmichel/projects/rbge/tmichel/reference_genomes/Hannah_Begonia_baits_edited.fasta

# Get bait names from FASTA
grep ">" $BAITS | sed 's/>//' | sort > baits_all.txt

# Get bait names from BED
cut -f4 bait_hits.bed | sort -u > baits_mapped.txt

# Find the missing ones
comm -23 baits_all.txt baits_mapped.txt > baits_unmapped.txt
wc -l baits_unmapped.txt
cat baits_unmapped.txt
