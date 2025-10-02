#!/usr/bin/env python3
"""
Generate probe candidate regions with GC content filtering
"""

import argparse
import os
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Generate probe candidates with GC filtering', add_help=False)
parser.add_argument('-h', '--help', action='help', help='Print help information')
parser.add_argument('-i', '--input', required=True, metavar='FILE', help='Input FASTA file')
parser.add_argument('-p', '--prefix', required=True, metavar='STRING', help='Prefix for output files and directory')
parser.add_argument('--min_gc', type=float, default=45, metavar='FLOAT', help='Min GC percentage (default: 45)')
parser.add_argument('--max_gc', type=float, default=55, metavar='FLOAT', help='Max GC percentage (default: 55)')
parser.add_argument('-o', '--outdir', required=True, metavar='STRING', help='Output directory')
args = parser.parse_args()

def calculate_gc_content(sequence):
    """Calculate GC content percentage"""
    gc_count = sequence.upper().count('G') + sequence.upper().count('C')
    return (gc_count / len(sequence)) * 100

def filter_duplicated_probes(probe_sequences):
    """Retain unique sequences after filtering duplicated sequences"""
    seen_sequences = set()
    unique_probes = {}
    for probe_id, sequence in probe_sequences.items():
        if sequence not in seen_sequences:
            seen_sequences.add(sequence)
            unique_probes[probe_id] = sequence
    return unique_probes

# Create output directory and filename
outdir = args.outdir
temp_dir = os.path.join(outdir, 'temp')
os.makedirs(temp_dir, exist_ok=True)
gc_suffix = f"gc{int(args.min_gc)}-{int(args.max_gc)}"
outfile = os.path.join(temp_dir, f"{args.prefix}_probe_candidates_{gc_suffix}.fasta")

probe_count = 0
total_count = 0
gc_filtered_probes ={}

with open(outfile, 'w') as out_f:
    for record in SeqIO.parse(args.input, 'fasta'):
        sequence = str(record.seq).upper()
        seq_len = len(sequence)
        
        # Generate sliding window probes (52 nt)
        for i in range(seq_len - 51):
            probe_seq = sequence[i:i + 52]
            
            # Split into P1 (first 25 nt) and P2 (last 25 nt, skip 2 nt spacer)
            p1_seq, p2_seq = probe_seq[:25], probe_seq[27:]
            total_count += 1
            
            # Filter by GC percent
            if (args.min_gc <= calculate_gc_content(p1_seq) <= args.max_gc and 
                args.min_gc <= calculate_gc_content(p2_seq) <= args.max_gc):
                probe_id = f"{record.id}_probe_{i+1}_{i+52}"
                gc_filtered_probes[probe_id] = probe_seq
    
    # Filter duplicated probes
    unique_probes = filter_duplicated_probes(gc_filtered_probes)

    # Write unique probes to file
    for probe_id, probe_seq in unique_probes.items():
        out_f.write(f">{probe_id}\n{probe_seq}\n")
        probe_count += 1

print(f"Generated {len(gc_filtered_probes)} probe candidates after GC filtering")
print(f"Kept {probe_count} unique probe candidates -> {outfile}")