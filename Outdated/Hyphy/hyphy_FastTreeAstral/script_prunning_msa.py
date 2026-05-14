from Bio import SeqIO

input_file = "cleaned_msa/ACmerged_contig_10072_NOSTOP.fna"
output_file = "cleaned_msa/ACmerged_contig_10072_NOSTOP_prunned.fna"
sample_to_remove = "Hannah_Begonia_baits__ref"

with open(input_file) as infile, open(output_file, "w") as outfile:
    for record in SeqIO.parse(infile, "fasta"):
        if record.id != sample_to_remove:
            SeqIO.write(record, outfile, "fasta")
