#!/bin/bash
#SBATCH --job-name="check"
#SBATCH --export=ALL
#SBATCH --mem=8G  
#SBATCH --partition=medium

# Checking that the baits don't match in parts very distant in the genomes (on different contigs) which cannot be accounted as result of the intron splicing.

# Count how many distinct genomic regions each bait maps to
cut -f4 bait_hits.bed | sort | uniq -c | sort -rn | head -20

# Count baits with more than one genomic locus
cut -f4 bait_hits.bed | sort -u | wc -l
# compare to total number of baits in your FASTA
grep -c ">" $BAITS
