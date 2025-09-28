#!/usr/bin/env python3
"""
HCR Probe Design Pipeline
"""

import subprocess
import argparse
import os

def run_step(step_name, cmd):
    """Run pipeline step with error handling"""
    print(f"\n--- {step_name} ---")
    try:
        subprocess.run(cmd, check=True)
        print(f"{step_name} completed!")
    except subprocess.CalledProcessError as e:
        print(f"Error in {step_name}: {e}")
        exit(1)

parser = argparse.ArgumentParser(description='HCR Probe Design Pipeline')
parser.add_argument('-i', '--input', required=True, help='Input FASTA file')
parser.add_argument('-d', '--database', required=True, help='BLAST database or FASTA file')
parser.add_argument('-p', '--prefix', required=True, help='Prefix for output files')
parser.add_argument('--initiator_id', required=True, help='HCR initiator ID')
parser.add_argument('--initiator_custom', help='Custom initiators CSV file')
parser.add_argument('--target_ids', help='On-target mRNA ID list file')
parser.add_argument('--make_db', action='store_true', help='Create BLAST database')
parser.add_argument('--min_gc', type=float, default=45, help='Min GC percentage')
parser.add_argument('--max_gc', type=float, default=55, help='Max GC percentage')
parser.add_argument('-t', '--threads', type=int, default=1, help='CPU threads for BLAST')

args = parser.parse_args()

# Set output directory
outdir = f'{args.prefix}_out'
os.makedirs(outdir, exist_ok=True)
gc_suffix = f"gc{int(args.min_gc)}-{int(args.max_gc)}"

# Step 1: GC content filtering
cmd1 = ['python3', 'src/step01_filter_gc_content.py', 
        '-i', args.input,
        '-p', args.prefix,
        '--min_gc', str(args.min_gc), 
        '--max_gc', str(args.max_gc), 
        '-o', outdir]
run_step("Step 1: GC Content Filtering", cmd1)

# Step 2: BLAST search
probe_fasta = os.path.join(outdir, f"{args.prefix}_probe_candidates_{gc_suffix}.fasta")
cmd2 = ['python3', 'src/step02_blast_search.py',
        '--probe_fasta', probe_fasta, 
        '-d', args.database, 
        '-p', args.prefix,
        '-t', str(args.threads), 
        '--min_gc', str(args.min_gc), 
        '--max_gc', str(args.max_gc), 
        '-o', outdir]
if args.make_db: cmd2.append('--make_db')
run_step("Step 2: BLAST Search", cmd2)

# Step 3: Filter and generate probes
blast_results = os.path.join(outdir, f"{args.prefix}_blast_results_{gc_suffix}.tsv")
cmd3 = ['python3', 'src/step03_filter_and_generate.py',
        '-i', args.input,
        '--blast_results', blast_results, 
        '--probe_fasta', probe_fasta, 
        '-p', args.prefix, 
        '--initiator_id', args.initiator_id, 
        '--min_gc', str(args.min_gc), 
        '--max_gc', str(args.max_gc), 
        '-o', outdir]
if args.target_ids: cmd3.extend(['--target_ids', args.target_ids])
if args.initiator_custom: cmd3.extend(['--initiator_custom', args.initiator_custom])
run_step("Step 3: Filter and Generate Probes", cmd3)

print("\nhcrkit completed!")