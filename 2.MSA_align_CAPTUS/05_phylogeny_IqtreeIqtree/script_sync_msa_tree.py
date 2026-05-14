#!/usr/bin/env python3
"""
sync_msa_tree_normalize.py

For each locus listed (one per line) in mfa_files.txt this script:
 - reads cleaned_msa/<locus>_NOSTOP_prunned.fna (FASTA)
 - reads gene_trees/<locus>_pruned.treefile (Newick)
 - normalizes tree tip names (dots -> underscores, strip quotes/spaces)
 - finds the intersection of MSA names and normalized tree names
 - writes synced_msa/<locus>_synced.fna  (only shared sequences)
 - writes synced_trees/<locus>_synced.treefile  (tree pruned to same taxa)

Original files are NOT overwritten.
"""
import os
import sys
from Bio import SeqIO, Phylo

msa_dir = "cleaned_msa"
tree_dir = "gene_trees"
out_msa_dir = "synced_msa"
out_tree_dir = "synced_trees"
os.makedirs(out_msa_dir, exist_ok=True)
os.makedirs(out_tree_dir, exist_ok=True)

def normalize_name(s):
    if s is None:
        return None
    # basic normalization: dots -> underscores, strip quotes and spaces
    norm = s.replace('.', '_').replace(' ', '_').strip().strip('"').strip("'")
    return norm

def process_locus(locus):
    msa_path = os.path.join(msa_dir, f"{locus}_NOSTOP_prunned.fna")
    tree_path = os.path.join(tree_dir, f"{locus}_pruned.treefile")
    if not os.path.exists(msa_path) or not os.path.exists(tree_path):
        print(f"Skipping {locus}: missing MSA or tree", file=sys.stderr)
        return

    # Read MSA headers
    msa_records = list(SeqIO.parse(msa_path, "fasta"))
    msa_ids = [rec.id for rec in msa_records]  # typical FASTA id (before whitespace)
    msa_set = set(msa_ids)

    # Read tree
    tree = Phylo.read(tree_path, "newick")

    # Build mapping original_tree_name -> normalized_name
    orig_to_norm = {}
    for term in tree.get_terminals():
        if term.name is None:
            continue
        orig_to_norm[term.name] = normalize_name(term.name)

    # Build reverse mapping normalized -> [originals]
    norm_to_orig = {}
    for orig, norm in orig_to_norm.items():
        norm_to_orig.setdefault(norm, []).append(orig)

    # Intersection between msa ids and normalized tree names
    common_norm = sorted(set(msa_set).intersection(norm_to_orig.keys()))
    if len(common_norm) == 0:
        print(f"{locus}: NO common taxa after normalization (0).", file=sys.stderr)
        return

    # Filter MSA records (keep only those in common_norm)
    common_set = set(common_norm)
    filtered = [rec for rec in msa_records if rec.id in common_set]
    out_msa = os.path.join(out_msa_dir, f"{locus}_synced.fna")
    SeqIO.write(filtered, out_msa, "fasta")
    print(f"{locus}: wrote {len(filtered)} sequences to {out_msa}")

    # Rename terminals in the tree to their normalized names
    for term in tree.get_terminals():
        if term.name is None:
            continue
        term.name = orig_to_norm.get(term.name, term.name)

    # Prune any terminal whose normalized name is NOT in common_set
    term_names = [t.name for t in tree.get_terminals() if t.name]
    for name in term_names:
        if name not in common_set:
            try:
                tree.prune(target=name)
            except Exception as e:
                # best-effort: print warning but continue
                print(f"Warning pruning {name} in {locus}: {e}", file=sys.stderr)

    out_tree = os.path.join(out_tree_dir, f"{locus}_synced.treefile")
    Phylo.write(tree, out_tree, "newick")
    print(f"{locus}: wrote pruned tree with {len(tree.get_terminals())} terminals to {out_tree}")


if __name__ == "__main__":
    if not os.path.exists("mfa_files.txt"):
        print("mfa_files.txt not found in current dir", file=sys.stderr)
        sys.exit(1)
    with open("mfa_files.txt") as fh:
        loci = [line.strip() for line in fh if line.strip()]
    for locus in loci:
        process_locus(locus)

