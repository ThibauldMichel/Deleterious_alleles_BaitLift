from Bio import SeqIO
import os

# -----------------------------
# PARAMETERS
# -----------------------------
input_dir = "~/scratch/Deleterious_alleles_PNG/CAPTUS/CAPTUS-PNG/05_phylogeny_FastTreeAstral/MSA_renamed"
missing_threshold = 0.7  # 70% missing data

# Expand ~ to home directory
input_dir = os.path.expanduser(input_dir)

# Output directory: current working directory
output_dir = os.path.join(os.getcwd(), "filtered_MSA_70")
os.makedirs(output_dir, exist_ok=True)

# Summary file in current working directory
summary_file = os.path.join(os.getcwd(), "filter_summary.txt")

print("Input directory:", input_dir)
print("Output directory:", output_dir)
print("Missing data threshold:", missing_threshold)

# -----------------------------
# START FILTERING
# -----------------------------
with open(summary_file, "w") as summary:
    summary.write("MSA_file\tkept_sequences\tremoved_sequences\tremoved_sample_names\n")

    try:
        files = os.listdir(input_dir)
    except Exception as e:
        print("ERROR: Cannot list files in input_dir:", e)
        exit(1)

    if not files:
        print("WARNING: No files found in the input directory!")

    for fasta_file in files:
        # Skip non-FASTA files
        if not (fasta_file.endswith(".fa") or fasta_file.endswith(".fasta") or fasta_file.endswith(".fna")):
            print(f"Skipping non-FASTA file: {fasta_file}")
            continue

        in_path = os.path.join(input_dir, fasta_file)
        out_path = os.path.join(output_dir, fasta_file)

        kept_seqs = []
        removed_seqs = []

        print(f"\nProcessing file: {fasta_file}")

        try:
            seq_count = 0
            for record in SeqIO.parse(in_path, "fasta"):
                seq_count += 1
                seq = str(record.seq)
                if len(seq) == 0:
                    print(f"WARNING: Sequence {record.id} is empty, skipping")
                    removed_seqs.append(record.id)
                    continue

                missing_fraction = seq.count('-') / len(seq)

                # Debug info
                if missing_fraction > missing_threshold:
                    print(f"Removing {record.id}: {missing_fraction*100:.1f}% missing")
                    removed_seqs.append(record.id)
                else:
                    kept_seqs.append(record)

            print(f"Total sequences in file: {seq_count}")
            print(f"Kept: {len(kept_seqs)}, Removed: {len(removed_seqs)}")

            # Write filtered sequences
            SeqIO.write(kept_seqs, out_path, "fasta")

            # Update summary
            summary.write(f"{fasta_file}\t{len(kept_seqs)}\t{len(removed_seqs)}\t{','.join(removed_seqs)}\n")

        except Exception as e:
            print(f"ERROR processing file {fasta_file}:", e)

print("\nFiltering complete! Summary written to", summary_file)

