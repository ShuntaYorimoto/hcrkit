#!/usr/bin/env python3
"""
HCR Probe Design Pipeline
"""

# Version information
__version__ = "1.0.0"

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

parser = argparse.ArgumentParser(description='HCR Probe Design Pipeline', add_help=False)
parser.add_argument('-h', '--help', action='help', help='Print help information')
parser.add_argument('-i', '--input', required=True, metavar='FILE', help='Input FASTA file')
parser.add_argument('-d', '--database', required=True, metavar='FILE', help='BLAST database or FASTA file')
parser.add_argument('-p', '--prefix', required=True, metavar='STRING', help='Prefix for output files')
parser.add_argument('--initiator_id', required=True, metavar='ID', 
                    help='HCR initiator ID. Predefined: S23, S41, S45, S72, S73, A161. Custom IDs available with --initiator_custom')
parser.add_argument('--initiator_custom', metavar='FILE', help='Custom initiators CSV file')
parser.add_argument('--initiator_split', type=int, default=9, 
                    help='Position to split initiator sequence between P1 and P2 (default: 9)')
parser.add_argument('--target_ids', metavar='FILE', help='On-target mRNA ID list file')
parser.add_argument('--make_db', metavar='STRING', help='Create BLAST database with specified name')
parser.add_argument('--min_gc', type=float, default=45, metavar='FLOAT', help='Min GC percentage (default: 45)')
parser.add_argument('--max_gc', type=float, default=55, metavar='FLOAT', help='Max GC percentage (default: 55)')
parser.add_argument('-t', '--threads', type=int, default=1, metavar='INT', help='CPU threads for BLAST')
parser.add_argument('-v', '--version', action='version', version=f'hcrkit {__version__}', 
                    help='Print version information')

args = parser.parse_args()

# Set output directory
outdir = f'{args.prefix}_out'
os.makedirs(outdir, exist_ok=True)
gc_suffix = f"gc{int(args.min_gc)}-{int(args.max_gc)}"

# Step 1: GC content filtering
cmd1 = ['step01_filter_gc_content.py', 
        '-i', args.input,
        '-p', args.prefix,
        '--min_gc', str(args.min_gc), 
        '--max_gc', str(args.max_gc), 
        '-o', outdir]
run_step("Step 1: GC Content Filtering", cmd1)

# Step 2: BLAST search
temp_dir = os.path.join(outdir, 'temp')
probe_fasta = os.path.join(temp_dir, f"{args.prefix}_probe_candidates_{gc_suffix}.fasta")
cmd2 = ['step02_blast_search.py',
        '--probe_fasta', probe_fasta, 
        '-d', args.database, 
        '-p', args.prefix,
        '-t', str(args.threads), 
        '--min_gc', str(args.min_gc), 
        '--max_gc', str(args.max_gc), 
        '-o', outdir]
if args.make_db: cmd2.extend(['--make_db', args.make_db])
run_step("Step 2: BLAST Search", cmd2)

# Step 3: Filter and generate probes
blast_results = os.path.join(temp_dir, f"{args.prefix}_blast_results_{gc_suffix}.tsv")
cmd3 = ['step03_filter_and_generate.py',
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
if args.initiator_split: cmd3.extend(['--initiator_split', str(args.initiator_split)])
run_step("Step 3: Filter and Generate Probes", cmd3)

print("\nhcrkit completed!")