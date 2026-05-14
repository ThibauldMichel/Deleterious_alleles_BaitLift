#!/usr/bin/env bash

# Directory containing .mfa.rates files
DIR="./mfa_dir"

# Output CSV file
OUT="all_chromosomes_rates_nonzero.csv"

# Header
echo "chromosome,neutral_rate,observed_rate,RS" > "$OUT"

# Loop through rate files
for f in "$DIR"/*.mfa.rates; do
    chrom=$(basename "$f" .mfa.rates)

    # Remove trailing "_edited" if present
    chrom=${chrom%_edited}

    while read -r neutral observed; do
        RS=$(awk -v n="$neutral" -v o="$observed" 'BEGIN {print n - o}')

        # Only print lines where RS != 0
        awk -v r="$RS" 'BEGIN {exit !(r != 0)}' && \
        echo "${chrom},${neutral},${observed},${RS}" >> "$OUT"

    done < "$f"
done

