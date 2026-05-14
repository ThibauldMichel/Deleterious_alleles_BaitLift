

vcf_path="/home/thibauld/Documents/Bioinformatics/Deleterious_alleles_pipeline/Deleterious_alleles_PNG/conserved_SNPs/filtered_in_regions.vcf"



bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/ANN\t[%GT]\n' $vcf_path > variants.tsv




