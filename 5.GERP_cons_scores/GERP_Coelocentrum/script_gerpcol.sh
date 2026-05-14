#!/bin/bash

#SBATCH --job-name="GERP++"
#SBATCH --export=ALL
#SBATCH --mem=200G
#SBATCH --partition=short



# Path of GERPplusplus
pathgerp="$HOME/scratch/GERPplusplus"

# Path to the tree file
pathtree="$HOME/projects/rbge/tmichel/CAPTUS/CAPTUS-Coelocentrum/05_phylogeny/concat.treefile"

# Path to the alignment files
pathalig="$HOME/projects/rbge/tmichel/CAPTUS/CAPTUS-Coelocentrum/04_alignments/03_trimmed/03_informed_w_refs/01_coding_NUC/02_NT"

# Path to the reference sequence
#pathrefe="$HOME/projects/rbge/tmichel/reference_genomes/Hannah_Begonia_baits.fasta"

# Copy the alignment files to the current directory and rename them as mfa.
mkdir -p mfa_dir  # Ensure the destination directory exists

# Copy all .fna files to mfa_dir and rename them to .mfa
for file in "$pathalig"/*.fna; do
    filename=$(basename "$file") ;
    rsync -azvh "$file" "mfa_dir/${filename%.fna}.mfa"
done 

# Create a list of the mfa files without path and extension.
for file in mfa_dir/*.mfa; do basename "$file" .mfa; done > mfa_files.txt

# Remove the extra squarre brackets from the header of the mfa files.
cat mfa_files.txt | while read file;
do
sed -E 's/^(>[^ ]+).*/\1/' "mfa_dir/${file}.mfa"  > "mfa_dir/${file}_edited.mfa";
done





# Calculate the RS score for each locus and position.
# The debug of gerpcol:
# https://www.biostars.org/p/207518/ 
cat mfa_files.txt | while read file;
do
"$pathgerp"/gerpcol -a -v -e Hannah_Begonia_baits__ref -f "mfa_dir/${file}_edited.mfa" -t "$pathtree" -x ".rates" ;
done

# -v verbose mode
# -a 	 alignment in mfa format [default = false]
# -e <reference seq>
#    	 name of reference sequence
# -f <filename>
#    	 alignment filename
# -t <tree filename>
#    	 evolutionary tree
# -x <suffix>
#    	 suffix for naming output files [default = ".rates"]


# From manual:
#The tree file should contain exactly one line specifying the tree in the standard nested parentheses/Newick format, including branch lengths.


#These output files contain one line for each position of the alignment. Each #line consists of the neutral rate N (from step 2) and RS score S (from step 4), #separated by a tab character, for that alignment position.  Sample output may #look like this:

#1.05	-1.1
#1.05	1.05
#1.23	1.23
#1.23	-1.98
#1.23	-1.23
#1.23	0.288

# Make the list of elements constrained found for each locus.
cat mfa_files.txt | while read file;
do
"$pathgerp"/gerpelem -f "mfa_dir/${file}_edited.mfa.rates"
done

# The output file(s) will contain one line for each constrained element found, # # listed in order of increasing p-value.  We report the following values for each # element, separated by spaces:

# start	end	length	     RS-score	p-value


cat mfa_files.txt | while read file;
do
 "mfa_dir/${file}_edited.mfa.rates.elems"
done


cat mfa_files.txt | while read file;
do
    awk -v name="$file" '{print name, $0}' OFS='\t' "mfa_dir/${file}_edited.mfa.rates.elems" 
done > all.mfa.rates.elems
