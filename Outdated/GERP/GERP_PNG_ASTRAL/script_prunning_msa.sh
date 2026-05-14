#!/bin/bash

#SBATCH --job-name="prunemsa"
#SBATCH --export=ALL
#SBATCH --partition=short

source activate newick_env


#python3 script_prunning_msa.py

sample_to_remove="Hannah_Begonia_baits__ref"

while read -r locus; do
  input_file="cleaned_msa/${locus}_NOSTOP.fna"
  output_file="cleaned_msa/${locus}_NOSTOP_prunned.fna"

  echo "Processing $locus"

  python3 - <<EOF
from Bio import SeqIO

input_file = "${input_file}"
output_file = "${output_file}"
sample_to_remove = "${sample_to_remove}"

with open(input_file) as infile, open(output_file, "w") as outfile:
    for record in SeqIO.parse(infile, "fasta"):
        if record.id != sample_to_remove:
            SeqIO.write(record, outfile, "fasta")
EOF

done < mfa_files.txt
