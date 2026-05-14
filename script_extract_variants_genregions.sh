#!/bin/bash
#SBATCH --job-name="compressVCF"
#SBATCH --export=ALL
#SBATCH --mem=200G  
#SBATCH --partition=medium

# Step 1: Extract chromosome names
bcftools view -h calls/all.vcf.gz | \
grep "^##contig" | \
awk -F'[=,>]' '{print $3}' > chromosomes.txt

# Step 2: Create rename mapping file
awk 'BEGIN{OFS="\t"} {print $1, gensub(/Hannah_Begonia_baits-/, "", 1, $1)}' chromosomes.txt > rename_chrs.txt

# Step 3: Apply renaming
bcftools annotate --rename-chrs rename_chrs.txt calls/all.vcf.gz -Oz -o calls/all_renamed.vcf.gz

# Index the new file
tabix calls/all_renamed.vcf.gz

# Extract the variants matching intersection of HyPhy and GERP
bcftools view -R intersection_GERP_HyPhy.bed calls/all_renamed.vcf.gz > calls/bait_variants.vcf
