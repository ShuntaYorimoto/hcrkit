#!/usr/bin/env python3
"""
BLAST search for probe candidate regions
"""

import argparse
import subprocess
import os

parser = argparse.ArgumentParser(description='BLAST probe candidates against database')
parser.add_argument('--probe_fasta', required=True, help='Input probe candidates FASTA file')
parser.add_argument('-d', '--database', required=True, help='BLAST database or FASTA file')
parser.add_argument('-p', '--prefix', required=True, help='Prefix for output files')
parser.add_argument('--make_db', action='store_true', help='Create BLAST database from FASTA file')
parser.add_argument('-t', '--threads', type=int, default=1, help='Number of CPU threads (default: 1)')
parser.add_argument('--min_gc', type=float, default=45, help='Min GC percentage (default: 45)')
parser.add_argument('--max_gc', type=float, default=55, help='Max GC percentage (default: 55)')
parser.add_argument('-o', '--outdir', required=True, help='Output directory')
args = parser.parse_args()

# Create output directory and filename
outdir = args.outdir
gc_suffix = f"gc{int(args.min_gc)}-{int(args.max_gc)}"
blast_out = os.path.join(outdir, f"{args.prefix}_blast_results_{gc_suffix}.tsv")

# Create BLAST database if requested
if args.make_db:
    db_name = os.path.splitext(args.database)[0]
    print("Creating BLAST database...")
    cmd = [
        'makeblastdb',
        '-in', args.database,
        '-dbtype', 'nucl',
        '-out', db_name
    ]
    subprocess.run(cmd, check=True)
    database = db_name
    print(f"BLAST database created: {db_name}")
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