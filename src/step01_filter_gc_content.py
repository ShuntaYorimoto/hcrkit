#!/usr/bin/env python3
"""
Generate probe candidate regions with GC content filtering
"""

import argparse
import os
from Bio import SeqIO

def calculate_gc_content(sequence):
    """Calculate GC content percentage"""
    gc_count = sequence.upper().count('G') + sequence.upper().count('C')
    return (gc_count / len(sequence)) * 100

parser = argparse.ArgumentParser(description='Generate probe candidates with GC filtering')
parser.add_argument('-i', '--input', required=True, help='Input FASTA file')
parser.add_argument('-p', '--prefix', required=True, help='Prefix for output files and directory')
parser.add_argument('--min_gc', type=float, default=45, help='Min GC percentage (default: 45)')
parser.add_argument('--max_gc', type=float, default=55, help='Max GC percentage (default: 55)')
parser.add_argument('-o', '--outdir', required=True, help='Output directory')
args = parser.parse_args()

# Create output directory and filename
outdir = args.outdir
gc_suffix = f"gc{int(args.min_gc)}-{int(args.max_gc)}"
outfile = os.path.join(outdir, f"{args.prefix}_probe_candidates_{gc_suffix}.fasta")

probe_count = 0
total_count = 0

with open(outfile, 'w') as out_f:
    for record in SeqIO.parse(args.input, 'fasta'):
        sequence = str(record.seq).upper()
        seq_len = len(sequence)
        
        # Generate sliding window probes (52bp)
        for i in range(seq_len - 52 + 1):
            probe_seq = sequence[i:i + 52]
            
            # Split into P1 (first 25bp) and P2 (last 25bp, skip 2bp spacer)
            p1_seq = probe_seq[0:25]
            p2_seq = probe_seq[27:52]
            
            p1_gc = calculate_gc_content(p1_seq)
            p2_gc = calculate_gc_content(p2_seq)
            total_count += 1
            
            # Filter: both P1 and P2 must meet GC criteria
            if (args.min_gc <= p1_gc <= args.max_gc and 
                args.min_gc <= p2_gc <= args.max_gc):
                probe_id = f"{record.id}_probe_{i+1}_{i+52}"
                out_f.write(f">{probe_id}\n{probe_seq}\n")
                probe_count += 1

print(f"Generated {probe_count}/{total_count} probe candidates -> {outfile}")