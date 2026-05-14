#!/bin/bash
#SBATCH --job-name="SIFT_annotation"
#SBATCH --mem=8G
#SBATCH --partition=medium
#SBATCH --export=ALL

VCF=/home/tmichel/scratch/Deleterious_alleles_PNG_baits/6.Variant_calling_annotations/all_renamed.vcf.gz
TSV=/home/tmichel/scratch/Deleterious_alleles_PNG_baits/1.Baits_to_genome/bait_locus_annotation.tsv
FAS=/home/tmichel/projects/rbge/tmichel/reference_genomes/Bmas.fa
GFF=/home/tmichel/projects/rbge/tmichel/reference_genomes/Bmas.gff
OUT=./


bash sift_annotation_workflow.sh $VCF $TSV $FAS $GFF $OUT