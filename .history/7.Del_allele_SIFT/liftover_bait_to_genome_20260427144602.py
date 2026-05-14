#!/usr/bin/env python3
"""
Liftover VCF from bait-relative coordinates to genome coordinates.

Input:
  - VCF file with variants called against bait sequences (from CAPTUS/HybPiper)
  - TSV mapping table (bait_locus_annotation.tsv) produced in Stage 1

Output:
  - VCF file with variants remapped to genome coordinates

Key logic:
  The TSV table has one row per GFF feature (CDS exon) overlapping each bait.
  All rows for the same bait share identical bait-level alignment fields
  (t_start, t_end, q_start, q_end, strand), while gff_start/gff_end differ
  per exon.

  For a variant at VCF position P (1-based) on bait B:
    1. Look up the bait's alignment: q_start, t_start, t_end, strand
    2. offset = (P - 1) - q_start          # convert to 0-based, then subtract bait start
    3. if + strand: genome_pos = t_start + offset + 1   (back to 1-based)
       if - strand: genome_pos = t_end   - offset       (t_end is 0-based exclusive -> already 1-based)
    4. Optionally annotate with whichever GFF exon spans genome_pos
"""

import sys
import gzip
import re
from collections import defaultdict


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def parse_mapping_table(mapping_file):
    """
    Parse bait_locus_annotation.tsv.

    Column order (0-based index):
      0  chr
      1  t_start       (0-based, start of bait alignment on genome)
      2  t_end         (0-based exclusive, end of bait alignment on genome)
      3  bait_id
      4  mapq
      5  strand
      6  q_start       (0-based, start of bait query used in alignment)
      7  q_end         (0-based exclusive)
      8  q_len
      9  gff_seqid
      10 gff_source
      11 gff_type
      12 gff_start     (1-based GFF coordinate, start of CDS exon)
      13 gff_end       (1-based GFF coordinate, end of CDS exon, inclusive)
      14 gff_score
      15 gff_strand
      16 gff_phase
      17 gff_attributes

    Returns:
        dict: {bait_id -> {
                'chr':     str,
                't_start': int,   # 0-based
                't_end':   int,   # 0-based exclusive
                'q_start': int,   # 0-based
                'q_end':   int,   # 0-based exclusive
                'q_len':   int,
                'strand':  str,
                'mapq':    int,
                'exons':   [ {'gff_start': int,   # 1-based inclusive
                               'gff_end':   int,   # 1-based inclusive
                               'gff_type':  str,
                               'gff_strand':str,
                               'gff_phase': int,
                               'parent':    str }, ... ]
               }
              }
    """
    baits = {}

    with open(mapping_file, 'r') as fh:
        fh.readline()  # skip header

        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue

            f = line.split('\t')
            if len(f) < 17:
                # Malformed row — skip with a warning
                print(f"[!] Skipping malformed row (only {len(f)} fields): {line[:80]}",
                      file=sys.stderr)
                continue

            bait_id = f[3]

            # ------------------------------------------------------------------
            # Bait-level alignment fields (identical for every row of same bait)
            # ------------------------------------------------------------------
            chr_id  = f[0]
            t_start = int(f[1])
            t_end   = int(f[2])
            mapq    = int(f[4])
            strand  = f[5]
            q_start = int(f[6])
            q_end   = int(f[7])
            q_len   = int(f[8])

            # ------------------------------------------------------------------
            # GFF exon-level fields (differ per row)
            # ------------------------------------------------------------------
            gff_type   = f[11]
            gff_start  = int(f[12])   # 1-based inclusive
            gff_end    = int(f[13])   # 1-based inclusive
            gff_strand = f[15]
            gff_phase  = int(f[16]) if f[16] not in ('.', '') else 0
            gff_attrs  = f[17] if len(f) > 17 else ''

            # Parse Parent= from attributes
            parent_match = re.search(r'Parent=([^;]+)', gff_attrs)
            parent = parent_match.group(1) if parent_match else '.'

            # Register bait-level info once
            if bait_id not in baits:
                baits[bait_id] = {
                    'chr':     chr_id,
                    't_start': t_start,
                    't_end':   t_end,
                    'q_start': q_start,
                    'q_end':   q_end,
                    'q_len':   q_len,
                    'strand':  strand,
                    'mapq':    mapq,
                    'exons':   [],
                }

            baits[bait_id]['exons'].append({
                'gff_start':  gff_start,
                'gff_end':    gff_end,
                'gff_type':   gff_type,
                'gff_strand': gff_strand,
                'gff_phase':  gff_phase,
                'parent':     parent,
            })

    print(f"[*] Loaded {len(baits)} bait sequences from mapping table", file=sys.stderr)
    return baits


# ---------------------------------------------------------------------------
# Coordinate lifting
# ---------------------------------------------------------------------------

def liftover_position(bait_id, vcf_pos, baits):
    """
    Convert a 1-based VCF position on a bait to a 1-based genome position.

    Args:
        bait_id : CHROM field from the bait-relative VCF
        vcf_pos : POS field (1-based) from the bait-relative VCF
        baits   : dict from parse_mapping_table()

    Returns:
        dict with keys (chr, pos, strand, gff_type, parent, exon_matched)
        or None if the bait is not in the table or position is out of range.
    """
    if bait_id not in baits:
        return None

    b = baits[bait_id]
    q_start = b['q_start']   # 0-based
    q_end   = b['q_end']     # 0-based exclusive
    t_start = b['t_start']   # 0-based
    t_end   = b['t_end']     # 0-based exclusive
    strand  = b['strand']

    # Convert VCF 1-based pos → 0-based
    pos_0 = vcf_pos - 1

    # Check the position falls within the aligned portion of the bait
    if not (q_start <= pos_0 < q_end):
        return None

    offset = pos_0 - q_start  # distance from start of aligned bait region

    if strand == '+':
        # Genome position (0-based) = t_start + offset
        genome_pos_1based = t_start + offset + 1
    else:
        # Reverse strand: the first bait base maps to t_end-1 (0-based)
        # t_end is 0-based exclusive, so last included base = t_end - 1
        genome_pos_1based = t_end - offset  # t_end - offset gives 1-based directly

    # Identify which GFF exon the genome position falls in (for annotation)
    matched_exon = None
    for exon in b['exons']:
        if exon['gff_start'] <= genome_pos_1based <= exon['gff_end']:
            matched_exon = exon
            break

    return {
        'chr':          b['chr'],
        'pos':          genome_pos_1based,
        'strand':       strand,
        'gff_type':     matched_exon['gff_type']   if matched_exon else 'intergenic',
        'parent':       matched_exon['parent']     if matched_exon else '.',
        'exon_matched': matched_exon is not None,
    }


# ---------------------------------------------------------------------------
# VCF processing
# ---------------------------------------------------------------------------

def rewrite_contig_headers(header_lines, baits):
    """
    Replace ##contig lines that reference bait IDs with genome scaffold contigs.
    Deduplicates by scaffold so each appears only once.
    """
    genome_contigs_seen = set()
    new_headers = []

    for line in header_lines:
        if line.startswith('##contig='):
            # Extract ID from ##contig=<ID=xxx,...>
            m = re.search(r'ID=([^,>]+)', line)
            if m:
                bait_id = m.group(1)
                if bait_id in baits:
                    scaffold = baits[bait_id]['chr']
                    if scaffold not in genome_contigs_seen:
                        genome_contigs_seen.add(scaffold)
                        new_headers.append(f'##contig=<ID={scaffold}>\n')
                    # Skip the original bait contig line
                    continue
        new_headers.append(line)

    return new_headers


def liftover_vcf(vcf_file, mapping_file, output_file, failed_file=None):
    """
    Main lifting function.

    Args:
        vcf_file    : input VCF (plain or .gz)
        mapping_file: TSV from Stage 1
        output_file : output VCF (.gz)
        failed_file : optional path to write un-liftable variants for QC
    """
    baits = parse_mapping_table(mapping_file)

    opener = gzip.open if vcf_file.endswith('.gz') else open

    lifted_count  = 0
    failed_count  = 0
    skipped_count = 0   # variants outside aligned bait region

    # Collect header lines first so we can rewrite ##contig entries
    header_lines   = []
    variant_lines  = []

    with opener(vcf_file, 'rt') as vcf_in:
        for line in vcf_in:
            if line.startswith('#'):
                header_lines.append(line)
            else:
                variant_lines.append(line)

    # Rewrite ##contig headers to use genome scaffold IDs
    header_lines = rewrite_contig_headers(header_lines, baits)

    fail_out = None
    if failed_file:
        fail_out = open(failed_file, 'w')

    with gzip.open(output_file, 'wt') as vcf_out:

        # Write (rewritten) headers
        for line in header_lines:
            vcf_out.write(line)

        for line in variant_lines:
            line = line.rstrip('\n')
            if not line:
                continue

            fields = line.split('\t')
            bait_id = fields[0]
            vcf_pos = int(fields[1])

            lifted = liftover_position(bait_id, vcf_pos, baits)

            if lifted is None:
                # Bait not in table OR position outside aligned region
                failed_count += 1
                if fail_out:
                    fail_out.write(line + '\n')
                continue

            # Update CHROM and POS
            fields[0] = lifted['chr']
            fields[1] = str(lifted['pos'])

            # Append provenance to INFO
            info = fields[7]
            provenance = (
                f"LIFTOVER_BAIT={bait_id}"
                f";LIFTOVER_BAIT_POS={vcf_pos}"
                f";LIFTOVER_STRAND={lifted['strand']}"
                f";LIFTOVER_GFFTYPE={lifted['gff_type']}"
                f";LIFTOVER_PARENT={lifted['parent']}"
            )
            fields[7] = '.' if info == '.' else info + ';'
            fields[7] = (provenance if info == '.'
                         else info + ';' + provenance)

            vcf_out.write('\t'.join(fields) + '\n')
            lifted_count += 1

    if fail_out:
        fail_out.close()

    print(f"[+] Lifted    : {lifted_count} variants", file=sys.stderr)
    print(f"[-] Failed    : {failed_count} variants (bait not in table / out of range)",
          file=sys.stderr)
    if failed_file:
        print(f"    Failed variants written to: {failed_file}", file=sys.stderr)
    print(f"[+] Output    : {output_file}", file=sys.stderr)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    if len(sys.argv) not in (4, 5):
        print(
            f"Usage: {sys.argv[0]} <input.vcf[.gz]> <mapping.tsv> "
            f"<output.vcf.gz> [failed_variants.vcf]",
            file=sys.stderr,
        )
        sys.exit(1)

    vcf_in   = sys.argv[1]
    map_tsv  = sys.argv[2]
    vcf_out  = sys.argv[3]
    failed   = sys.argv[4] if len(sys.argv) == 5 else None

    liftover_vcf(vcf_in, map_tsv, vcf_out, failed_file=failed)
