#!/usr/bin/env python3
"""
BLAST search for probe candidate regions
"""

import argparse
import subprocess
import os

parser = argparse.ArgumentParser(description='BLAST probe candidates against database', add_help=False)
parser.add_argument('-h', '--help', action='help', help='Print help information')
parser.add_argument('--probe_fasta', required=True, metavar='FILE', help='Input probe candidates FASTA file')
parser.add_argument('-d', '--database', required=True, metavar='FILE', help='BLAST database or FASTA file')
parser.add_argument('-p', '--prefix', required=True, metavar='STRING', help='Prefix for output files')
parser.add_argument('--make_db', metavar='STRING', help='Create BLAST database with specified name')
parser.add_argument('-t', '--threads', type=int, default=1, metavar='INT', help='Number of CPU threads (default: 1)')
parser.add_argument('--min_gc', type=float, default=45, metavar='FLOAT', help='Min GC percentage (default: 45)')
parser.add_argument('--max_gc', type=float, default=55, metavar='FLOAT', help='Max GC percentage (default: 55)')
parser.add_argument('-o', '--outdir', required=True, metavar='STRING', help='Output directory')
args = parser.parse_args()

# Create output directory and filename
outdir = args.outdir
temp_dir = os.path.join(outdir, 'temp')
os.makedirs(temp_dir, exist_ok=True)
gc_suffix = f"gc{int(args.min_gc)}-{int(args.max_gc)}"
blast_out = os.path.join(temp_dir, f"{args.prefix}_blast_results_{gc_suffix}.tsv")

# Create BLAST database if requested
if args.make_db:
    db_dir = f"{args.make_db}_blastdb"
    os.makedirs(db_dir, exist_ok=True)
    db_path = os.path.join(db_dir, args.make_db)
    print(f"Creating BLAST database: {args.make_db}")
    cmd = [
        'makeblastdb',
        '-in', args.database,
        '-dbtype', 'nucl',
        '-out', db_path
    ]
    subprocess.run(cmd, check=True)
    database = db_path
    print(f"BLAST database created: {db_path}")
else:
    database = args.database

# Run BLASTN
print("Running BLASTN...")
cmd = [
    'blastn',
    '-query', args.probe_fasta,
    '-db', database,
    '-task', 'blastn-short',
    '-outfmt', '6',
    '-evalue', '1.0e-2',
    '-num_threads', str(args.threads),
    '-out', blast_out
]

subprocess.run(cmd, check=True)
print(f"BLAST completed -> {blast_out}")