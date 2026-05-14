import json
import csv
import sys
import os
import numpy as np

# Benjamini-Hochberg correction for q-values
def benjamini_hochberg(pvalues):
    pvalues = np.array(pvalues)
    n = len(pvalues)
    sorted_indices = np.argsort(pvalues)
    sorted_pvalues = pvalues[sorted_indices]
    qvalues = np.empty(n)
    min_coeff = 1
    for i in reversed(range(n)):
        rank = i + 1
        coeff = n / rank
        qvalue = coeff * sorted_pvalues[i]
        min_coeff = min(min_coeff, qvalue)
        qvalues[sorted_indices[i]] = min_coeff
    qvalues[qvalues > 1] = 1
    return qvalues.tolist()

# Check if the user provided the input argument
if len(sys.argv) != 2:
    print("Usage: python3 script_meme_logs_parser.py PATTERN.json")
    sys.exit(1)

# Get the input filename from command-line argument
input_file = sys.argv[1]

# Extract the pattern name (without .json)
basename = os.path.basename(input_file)
pattern_name = os.path.splitext(basename)[0]

# Create output directory if it doesn't exist
output_dir = "./csv_extracted_MEME"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Full output path
output_file = os.path.join(output_dir, pattern_name + ".csv")

# Load JSON file
with open(input_file, "r") as f:
    data = json.load(f)

# Access content
rows = data["MLE"]["content"]["0"]

# Collect p-values for FDR correction (index 6)
p_values = [row[6] for row in rows]
q_values = benjamini_hochberg(p_values)

# Column names (with codon first, then omega)
columns = [
    "codon", "omega", "alpha", "beta-", "p-", "beta+", "p+", "LRT", "p-value", "q-value", "branches_with_selection", "LogL", "AICc"
]

# Prepare data rows
table = []
for idx, (row, q) in enumerate(zip(rows, q_values)):
    site_number = idx + 1
    alpha = row[0]
    beta_plus = row[3]
    # Compute omega safely
    if alpha != 0:
        omega = beta_plus / alpha
    else:
        omega = "NA"
    output_row = [
        site_number,  # codon
        omega,        # omega
        alpha,        # alpha
        row[1],       # beta-
        row[2],       # p-
        beta_plus,    # beta+
        row[4],       # p+
        row[5],       # LRT
        row[6],       # p-value
        q,            # q-value
        row[7],       # branches_with_selection
        row[9],       # LogL
        row[10]       # AICc
    ]
    table.append(output_row)

# Write CSV output
with open(output_file, "w", newline='') as f:
    writer = csv.writer(f)
    writer.writerow(columns)
    writer.writerows(table)

print(f"✅ Successfully parsed '{input_file}' to '{output_file}'")



