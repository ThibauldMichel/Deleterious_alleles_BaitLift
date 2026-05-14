#!/bin/bash
#SBATCH --job-name="goodgenesPNG"
#SBATCH --export=ALL
#SBATCH --mem=8G  
#SBATCH --partition=short



#conda create -n FastTree_env -c bioconda -c conda-forge fasttree
conda activate FastTree



mkdir -p good_genes_MSA


while read f;
do
	rsync -azvh /home/tmichel/scratch/Deleterious_alleles_PNG/CAPTUS/CAPTUS-PNG/05_phylogeny_FastTreeAstral/filtering_1_missing_data_70/filtered_MSA_70/"$f".fna good_genes_MSA/ ;
done < /home/tmichel/scratch/Deleterious_alleles_PNG/CAPTUS/CAPTUS-PNG/05_phylogeny_FastTreeAstral/filtering_2_missing_samples_75/good_genes.csv




# -nt	Nucleotide sequences	Tells FastTree to treat the alignment as nucleotide data (not amino acids). It uses nucleotide substitution models.
# -gtr	General Time Reversible model	Uses the full GTR model, which estimates six substitution rate parameters (A↔C, A↔G, A↔T, C↔G, C↔T, G↔T) and empirical base frequencies — the most general reversible model commonly used in phylogenetics.
# -gamma	Gamma-distributed rate heterogeneity	Models variation in substitution rates among sites with a discrete gamma distribution. This accounts for some sites evolving faster than others — an important realism that improves branch length estimates.
