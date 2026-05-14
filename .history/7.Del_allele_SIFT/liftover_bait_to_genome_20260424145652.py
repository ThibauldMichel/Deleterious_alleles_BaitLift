#!/usr/bin/env python3
"""
Liftover VCF from bait-relative coordinates to genome coordinates.

Input:
  - VCF file with variants in bait coordinates (from CAPTUS/HybPiper)
  - TSV mapping table from Stage 1 (bait_id, q_start, q_end, t_start, t_end, strand, gff features)

Output:
  - VCF file with variants in genome coordinates
  - Lifted positions map to exon features in the reference GFF

Key logic:
  For each variant at position P in bait Q:
    1. Find all exon blocks for bait Q in the mapping table
    2. Sort exons by query coordinate (q_start, q_end)
    3. Accumulate query bases until reaching position P
    4. Map to the corresponding target (genome) exon block
    5. Emit genome coordinate (chr, pos)
"""

import sys
import gzip
from collections import defaultdict
from pathlib import Path

def parse_mapping_table(mapping_file):
    """
    Parse the enriched coordinate mapping table.
    
    Returns:
        dict: {bait_id -> [(q_start, q_end, q_len, t_start, t_end, chr, strand, gff_type, gff_start, gff_end), ...]}
    """
    bait_exons = defaultdict(list)
    
    with open(mapping_file, 'r') as f:
        # Skip header
        header = f.readline()
        
        for line in f:
            fields = line.strip().split('\t')
            
            # Parse fields
            chr_id = fields[0]
            t_start = int(fields[1])
            t_end = int(fields[2])
            bait_id = fields[3]
            mapq = int(fields[4])
            strand = fields[5]
            q_start = int(fields[6])
            q_end = int(fields[7])
            q_len = int(fields[8])
            gff_seqid = fields[9]
            gff_type = fields[11]
            gff_start = int(fields[12])
            gff_end = int(fields[13])
            gff_strand = fields[15]
            gff_phase = int(fields[16]) if fields[16] != '.' else 0
            
            # Store exon info indexed by bait_id
            # Key: (q_start, q_end) to sort properly
            bait_exons[bait_id].append({
                'q_start': q_start,
                'q_end': q_end,
                'q_len': q_len,
                't_start': t_start,
                't_end': t_end,
                'chr': chr_id,
                'strand': strand,
                'gff_type': gff_type,
                'gff_start': gff_start,
                'gff_end': gff_end,
                'gff_strand': gff_strand,
                'gff_phase': gff_phase,
                'mapq': mapq,
            })
    
    # Sort exons by query coordinate for each bait
    for bait_id in bait_exons:
        bait_exons[bait_id].sort(key=lambda x: x['q_start'])
    
    return bait_exons


def liftover_position(bait_id, query_pos, bait_exons):
    """
    Lift a position from bait query coordinates to genome coordinates.
    
    Args:
        bait_id: identifier of the bait/contig
        query_pos: 1-based position in the bait sequence
        bait_exons: mapping dict from parse_mapping_table()
    
    Returns:
        (chr, genome_pos, strand, gff_type, exon_info) or None if unmappable
    """
    
    if bait_id not in bait_exons:
        return None
    
    exons = bait_exons[bait_id]
    
    # Find which exon block contains this query position
    # The exons may not be contiguous in query space (they could overlap or have gaps)
    # But typically for a spliced alignment they form a continuous path
    
    # Accumulate query bases as we walk through exons in order
    cumulative_query_base = 0
    
    for exon in exons:
        q_start = exon['q_start']
        q_end = exon['q_end']
        exon_length = q_end - q_start  # 0-based half-open
        
        # Check if query_pos falls within this exon
        # query_pos is 1-based from VCF; convert to 0-based for comparison
        query_pos_0based = query_pos - 1
        
        if q_start <= query_pos_0based < q_end:
            # Found the exon
            offset_in_exon = query_pos_0based - q_start
            
            # Map to target (genome) coordinate
            t_start = exon['t_start']
            t_end = exon['t_end']
            strand = exon['strand']
            
            if strand == '+':
                # Forward strand: offset in query maps directly to offset in target
                genome_pos = t_start + offset_in_exon
            else:
                # Reverse strand: offset in query maps to reverse direction in target
                # For - strand, the target coordinates are still in forward direction
                # but the alignment is reversed
                genome_pos = t_end - offset_in_exon - 1
            
            return {
                'chr': exon['chr'],
                'pos': genome_pos + 1,  # Convert back to 1-based for VCF
                'strand': strand,
                'gff_type': exon['gff_type'],
                'exon': exon,
            }
    
    # Position not found in any exon
    return None


def liftover_vcf(vcf_file, mapping_file, output_file):
    """
    Lift VCF from bait coordinates to genome coordinates.
    
    Args:
        vcf_file: input VCF (gzipped or plain)
        mapping_file: TSV mapping table from Stage 1
        output_file: output VCF (gzipped)
    """
    
    print(f"[*] Parsing mapping table from {mapping_file}", file=sys.stderr)
    bait_exons = parse_mapping_table(mapping_file)
    print(f"[*] Loaded {len(bait_exons)} bait sequences", file=sys.stderr)
    
    # Detect if VCF is gzipped
    if vcf_file.endswith('.gz'):
        vcf_open = gzip.open
    else:
        vcf_open = open
    
    lifted_count = 0
    failed_count = 0
    
    with vcf_open(vcf_file, 'rt') as vcf_in, gzip.open(output_file, 'wt') as vcf_out:
        
        for line in vcf_in:
            # Write header lines unchanged
            if line.startswith('##') or line.startswith('#CHROM'):
                vcf_out.write(line)
                continue
            
            if line.startswith('#'):
                continue
            
            # Parse VCF variant line
            fields = line.strip().split('\t')
            chrom = fields[0]  # This is the bait_id in bait-relative VCF
            pos = int(fields[1])  # 1-based position in bait
            ref = fields[3]
            alt = fields[4]
            
            # Try to liftover
            lifted = liftover_position(chrom, pos, bait_exons)
            
            if lifted is None:
                # Could not lift this variant
                failed_count += 1
                # Optionally write to a failed variants file for debugging
                continue
            
            # Update VCF fields with genome coordinates
            fields[0] = lifted['chr']  # chromosome
            fields[1] = str(lifted['pos'])  # 1-based genome position
            
            # Add INFO fields for provenance
            info = fields[7]
            if info == '.':
                info = ''
            else:
                info += ';'
            
            info += f"LIFTOVER=bait:{chrom}:{pos};GFFTYPE={lifted['gff_type']};STRAND={lifted['strand']}"
            fields[7] = info
            
            # Write lifted variant
            vcf_out.write('\t'.join(fields) + '\n')
            lifted_count += 1
    
    print(f"[+] Lifted {lifted_count} variants", file=sys.stderr)
    print(f"[-] Failed to lift {failed_count} variants", file=sys.stderr)
    print(f"[+] Output: {output_file}", file=sys.stderr)


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <input.vcf> <mapping.tsv> <output.vcf.gz>", file=sys.stderr)
        sys.exit(1)
    
    vcf_file = sys.argv[1]
    mapping_file = sys.argv[2]
    output_file = sys.argv[3]
    
    liftover_vcf(vcf_file, mapping_file, output_file)
