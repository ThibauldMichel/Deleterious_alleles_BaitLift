
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Rsamtools")
# install.packages('Rsamtools', repos = c('https://bioc.r-universe.dev', 'https://cloud.r-project.org'))

BiocManager::install("VariantAnnotation")

library(Rsamtools)
library(VariantAnnotation)

path_vcf <- "/home/thibauld/Documents/Bioinformatics/Deleterious_alleles_pipeline/Deleterious_alleles_PNG/SIFT/all_annotated_sample.vcf"
vcf <- scanVcf(path_vcf)


