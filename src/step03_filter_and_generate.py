#!/usr/bin/env python3
"""
Filter probes by specificity, select non-overlapping probes, and generate probe pairs
"""

import argparse
import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

parser = argparse.ArgumentParser(description='Filter and generate HCR probe pairs', add_help=False)
parser.add_argument('-h', '--help', action='help', help='Print help information')
parser.add_argument('-i','--input', required=True, metavar='FILE', help='Input FASTA file')
parser.add_argument('--blast_results', required=True, metavar='FILE', help='BLAST results TSV file')
parser.add_argument('--probe_fasta', required=True, metavar='FILE', help='Probe candidates FASTA file')
parser.add_argument('-p', '--prefix', required=True, metavar='STRING', help='Prefix for output files')
parser.add_argument('--initiator_id', required=True, metavar='ID', 
                    help='HCR initiator ID. Predefined: S23, S41, S45, S72, S73, A161. Custom IDs available with --initiator_custom')
parser.add_argument('--initiator_custom', metavar='FILE', help='Custom initiators CSV file')
parser.add_argument('--target_ids', metavar='FILE', help='On-target mRNA ID list file')
parser.add_argument('--min_gc', type=float, default=45, metavar='FLOAT', help='Min GC percentage (default: 45)')
parser.add_argument('--max_gc', type=float, default=55, metavar='FLOAT', help='Max GC percentage (default: 55)')
parser.add_argument('-o', '--outdir', required=True, metavar='STRING', help='Output directory')
args = parser.parse_args()

# Predefined HCR initiator sequences
INITIATORS = {
    'S23': 'GGGTGGTCGTCGAAGTCGTAT',
    'S41': 'GCTCGACGTTCCTTTGCAACA', 
    'S45': 'CCTCCACGTTCCATCTAAGCT',
    'S72': 'CGGTGGAGTGGCAAGTAGGAT',
    'S73': 'CGGTCAGGTGGCTAGTATGGA',
    'A161': 'GGTACGCGAAGGTAGGTGTAA'
}

def load_target_ids(target_ids_file, input_fasta_file):
    """Load on-target IDs from file or extract from original input FASTA"""
    if target_ids_file:
        with open(target_ids_file, 'r') as f:
            return set(line.strip() for line in f if line.strip())
    else:
        target_ids = set()
        for record in SeqIO.parse(input_fasta_file, 'fasta'):
            clean_id = record.id.split()[0]
            target_ids.add(clean_id)
        return target_ids

def load_probe_sequences(probe_fasta_file):
    """Load probe sequences from FASTA file"""
    probe_sequences = {}
    for record in SeqIO.parse(probe_fasta_file, 'fasta'):
        probe_sequences[record.id] = str(record.seq)
    return probe_sequences

def filter_by_blast_results(blast_file, target_ids):
    """Filter probes based on BLAST results and calculate off-target coverage"""
    # Initialize all probes as valid
    probe_coverage = {
        probe_id: {'max_offtarget': 0, 'valid': True}
        for probe_id in probe_sequences.keys()
    }
    
    with open(blast_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            query_id, subject_id, length = fields[0], fields[1], int(fields[3])
            coverage = length / 52
            
            # Check if off-target hit
            is_ontarget = any(target_id in subject_id for target_id in target_ids)
            if not is_ontarget:
                probe_coverage[query_id]['max_offtarget'] = max(
                    probe_coverage[query_id]['max_offtarget'], coverage * 100
                )
                if coverage >= 0.5:
                    probe_coverage[query_id]['valid'] = False
    
    return probe_coverage

def filter_overlapping_probes(valid_probes):
    """Select non-overlapping probes using greedy algorithm"""
    probe_list = []
    for probe_id in valid_probes:
        parts = probe_id.split('_')
        start = int(parts[-2])
        end = int(parts[-1])
        probe_list.append((start, end, probe_id))
    
    probe_list.sort()
    selected_probes = []
    last_end = 0

    for start, end, probe_id in probe_list:
        if start > last_end:
            selected_probes.append(probe_id)
            last_end = end
    
    return selected_probes

def create_probe_pairs_and_summary(final_probes, probe_sequences, probe_coverage, 
                                 initiator_seq, prefix, initiator_id):
    """Create probe pairs and summary data"""
    p1_init = initiator_seq[:9] + "aa"
    p2_init = "aa" + initiator_seq[9:]
    probe_pairs = []
    summary_data = []
    
    for probe_id in final_probes:
        probe_seq = probe_sequences[probe_id]
        start_pos = int(probe_id.split('_')[-2])
        
        p1_region = probe_seq[0:25]
        p2_region = probe_seq[27:52]
        p1_seq = p1_init + str(Seq(p1_region).reverse_complement())
        p2_seq = str(Seq(p2_region).reverse_complement()) + p2_init
        
        probe_set_name = f"{prefix}_{initiator_id}_s{start_pos}"
        probe_pairs.extend([
            [f"{probe_set_name}_P1", p1_seq],
            [f"{probe_set_name}_P2", p2_seq]
        ])
        
        max_coverage = probe_coverage.get(probe_id, {}).get('max_offtarget', 0)
        summary_data.append([probe_set_name, p1_seq, p2_seq, f"{max_coverage:.0f}"])
    
    return probe_pairs, summary_data

# Load initiators
if args.initiator_custom:
    initiators_df = pd.read_csv(args.initiator_custom, header=None, names=['ID', 'Sequence'])
    initiator_dict = dict(zip(initiators_df['ID'], initiators_df['Sequence']))
    print(f"Using custom initiators from: {args.initiator_custom}")
else:
    initiator_dict = INITIATORS
    print("Using predefined initiators")

if args.initiator_id not in initiator_dict:
    print(f"Error: Initiator '{args.initiator_id}' not found")
    exit(1)

initiator_seq = initiator_dict[args.initiator_id]
print(f"Selected initiator: {args.initiator_id} ({initiator_seq})")

# Filter valid probes
target_ids = load_target_ids(args.target_ids, args.input)
probe_sequences = load_probe_sequences(args.probe_fasta)
probe_coverage = filter_by_blast_results(args.blast_results, target_ids)
valid_probes = [probe_id for probe_id, data in probe_coverage.items() if data['valid']]
print(f"Valid probes after BLAST filtering: {len(valid_probes)}")

# filter overlapped probes
final_probes = filter_overlapping_probes(valid_probes)
print(f"Final non-overlapping probes: {len(final_probes)}")

probe_pairs, summary_data = create_probe_pairs_and_summary(
    final_probes, probe_sequences, probe_coverage, initiator_seq, args.prefix, args.initiator_id
)

# Output files
outdir = args.outdir
gc_suffix = f"gc{int(args.min_gc)}-{int(args.max_gc)}"
csv_out = os.path.join(outdir, f"{args.prefix}_{args.initiator_id}_probe_pairs_{gc_suffix}.csv")
probe_pairs_df = pd.DataFrame(probe_pairs, columns=['Oligoname', 'Sequence'])
probe_pairs_df.to_csv(csv_out, index=False)

summary_out = os.path.join(outdir, f"{args.prefix}_{args.initiator_id}_probe_summary_{gc_suffix}.txt")
with open(summary_out, 'w') as f:
    f.write("Probe set name\tP1\tP2\tMax off-target coverage (%)\n")
    for row in summary_data:
        f.write('\t'.join(row) + '\n')

print(f"Generated {len(probe_pairs)} split probe sequences -> {csv_out}")
print(f"Summary written to -> {summary_out}")