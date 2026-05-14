#!/bin/bash
#SBATCH --job-name="SIFT"
#SBATCH --export=ALL
#SBATCH --mem=32G



#java -jar \
#	SIFT4G_Annotator.jar \
#	-c \
#	-i all.vcf \
#	-d ~/scratch/Deleterious_alleles_PNG/SIFT/snpEff/data/Bmas \
#	-r results \
#	-t


# -c	Command-line mode (as opposed to GUI mode). Essential for headless environments.
# -i	Input file: path to the input VCF file (in your case, all.vcf).
# -d	Database directory: location of the SIFT4G database directory (here, CADRE.23). This folder should contain the pre-built SIFT4G database for the organism.
# -r	Results directory: directory where the annotation results will be saved (results). SIFT will generate an annotated VCF and possibly other output files here.
# -t	Use multiple threads (multithreading for faster annotation). Optional but recommended for speed.



java -Xmx32g -jar ./snpEff/snpEff.jar -v Bmas all.vcf > all_annotated.vcf

