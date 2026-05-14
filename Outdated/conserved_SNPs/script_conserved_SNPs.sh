#!/bin/bash
#SBATCH --job-name="bcftools"
#SBATCH --export=ALL
#SBATCH --mem=64G  
#SBATCH --partition=short

source /mnt/apps/users/tmichel/conda/etc/profile.d/conda.sh

conda activate bcftools_env


#bgzip /home/tmichel/scratch/Deleterious_alleles_PNG/SIFT/all_annotated.vcf
bgzip -c /home/tmichel/scratch/Deleterious_alleles_PNG/SIFT/all_annotated.vcf > \
/home/tmichel/scratch/Deleterious_alleles_PNG/SIFT/all_annotated.vcf.gz


bcftools index /home/tmichel/scratch/Deleterious_alleles_PNG/SIFT/all_annotated.vcf.gz



bcftools view \
  -R /home/tmichel/scratch/Deleterious_alleles_PNG/Merge_pipeline/Scaffold_conserved_loci.bed \
  -o filtered_in_regions.vcf \
  /home/tmichel/scratch/Deleterious_alleles_PNG/SIFT/all_annotated.vcf.gz

