#!/bin/bash
#SBATCH --job-name="minimap2"
#SBATCH --export=ALL
#SBATCH --mem=16G  
#SBATCH --partition=medium

source /mnt/apps/users/tmichel/conda/etc/profile.d/conda.sh

conda activate bcftools_env


BAITS=/home/tmichel/projects/rbge/tmichel/reference_genomes/Hannah_Begonia_baits_edited.fasta
REFERENCE=/home/tmichel/projects/rbge/tmichel/reference_genomes/Bmas.fa

minimap2 -ax splice \
         --cs \
         -C5 \
         -u b \
         $REFERENCE \
         $BAITS \
         | samtools sort -o bait_hits.bam

samtools index bait_hits.bam

# -a  Output in SAM format (instead of the default PAF format). Required if you want to pipe into samtools.

# -x splice, Uses a preset for spliced alignment.
# Originally designed for RNA-seq (aligning transcripts to genome with introns). Implies a bundle of parameters optimized for long gaps (introns) and exon-intron structure. In our case (baits → genome), this means it allows large gaps, which might or might not be desirable depending on whether our baits span introns.

# -C5 which reduces the penalty for non-canonical splice sites (important because splice sites may differ slightly between B. luzhaiensis and B. masoniana)

# --cs  adds a cs tag to each alignment. This tag encodes the base-level alignment differences (mismatches, insertions, deletions). It captures exact matches, substitutions, and indels. Useful to indicate fine-scale mismatch info and downstream parsing of sequence variation


# Remove secondary and supplementary alignments
samtools view -F 2308 -b bait_hits.bam > bait_hits_primary.bam
samtools index bait_hits_primary.bam

# -split reports each exon block separately — correct for GFF intersection
bedtools bamtobed -i bait_hits_primary.bam -split > bait_hits.bed

# The flag -F 2308 filters out unmapped (4), secondary (256), and supplementary (2048) alignments. The -split flag in bedtools tells it to report each exon block as a separate BED line but within the same record — which is what you want for downstream bedtools intersect with the GFF, since it will correctly identify which exons the bait covers rather than spanning across the intron.


# Results are given following this format in the BED file:
# chromosome   start   end   name   score   strand
# Most of the Mapping quality score (MAPQ from minimap2) are 60, indicating very high confidence, typically unique mapping
