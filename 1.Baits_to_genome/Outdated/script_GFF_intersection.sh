#!/bin/bash
#SBATCH --job-name="GFF intersection"
#SBATCH --export=ALL
#SBATCH --mem=8G  
#SBATCH --partition=medium

source /mnt/apps/users/tmichel/conda/etc/profile.d/conda.sh

conda activate bcftools_env

GFF=/home/tmichel/projects/rbge/tmichel/reference_genomes/Bmas.gff

bedtools intersect \
  -a bait_hits.bed \
  -b $GFF \
  -wa -wb \
  > bait_locus_annotation.tsv
  
  
  
  
  # -a — the "query" file. Your bait exon blocks BED. Each line here is tested for overlap against the -b file.
 # -b — the "database" file. Your GFF annotation. Every feature in the GFF (gene, mRNA, CDS, exon, UTR...) is a potential target for overlap.
 # -wa — write the original -a entry in full when an overlap is found. Without this flag, bedtools only reports the overlapping portion (the intersection coordinates). With -wa you keep your original bait coordinates intact in the output.
 # -wb — write the original -b entry in full alongside it. This is what gives you the GFF columns — chromosome, source, feature type, start, end, score, strand, frame, and the attributes field containing gene ID, transcript ID etc.
