#!/usr/bin/env bash

# Loop through all loci
for msa in ./synced_msa/*_synced.fna; do
  # Extract locus name (assuming filenames match before _NOSTOP_prunned)
  locus=$(basename "$msa" _synced.fna)
  tree="./synced_trees/${locus}_synced.treefile"

  if [[ ! -f "$tree" ]]; then
    echo "$locus    MISSING_TREE"
    continue
  fi

  # Count unique taxa in tree
  n_tree=$(tr '(),;' '\n' < "$tree" \
    | sed 's/:[^[:space:]]\+//g' \
    | sed '/^\s*$/d' \
    | grep -E '[A-Za-z]' \
    | sed 's/^[ \t]*//; s/[ \t]*$//' \
    | sort -u \
    | wc -l)

  # Count unique taxa in MSA
  n_msa=$(grep '^>' "$msa" \
    | sed 's/^>//; s/[[:space:]]*$//' \
    | sort -u \
    | wc -l)

  # Print side by side
  printf "%s\tTree:%s\tMSA:%s\n" "$locus" "$n_tree" "$n_msa"
done
