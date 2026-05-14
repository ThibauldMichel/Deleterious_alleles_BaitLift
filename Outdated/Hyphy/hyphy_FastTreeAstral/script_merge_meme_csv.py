import os
import pandas as pd

# Directory containing the csv files
csv_dir = './csv_extracted_MEME'

# Read the list of file names from mfa_files.txt
with open('mfa_files.txt', 'r') as f:
    samples = [line.strip() for line in f]

# Create a list to store all dataframes
all_dfs = []

for sample in samples:
    csv_file = os.path.join(csv_dir, f"{sample}.csv")
    
    if os.path.isfile(csv_file):
        # Read the csv file
        df = pd.read_csv(csv_file)
        
        # Insert new column with sample name at the first position
        df.insert(0, 'sample', sample)
        
        # Append dataframe to the list
        all_dfs.append(df)
    else:
        print(f"Warning: File {csv_file} not found. Skipping.")

# Concatenate all dataframes
if all_dfs:
    merged_df = pd.concat(all_dfs, ignore_index=True)
    
    # Write the merged dataframe to a new CSV file inside ./csv_extracted_MEME
    output_file = os.path.join(csv_dir, 'all_loci.csv')
    merged_df.to_csv(output_file, index=False)
    print(f"Successfully created {output_file}")
else:
    print("No CSV files were processed.")

