#!/usr/bin/env python3
"""
Extract mRNA IDs from GFF3 file based on exact gene name match
"""

import argparse
import os
import re

parser = argparse.ArgumentParser(description='Extract mRNA IDs from GFF3 based on exact gene match')
parser.add_argument('-i', '--input', required=True, help='Input GFF3 file')
parser.add_argument('-s', '--search', required=True, help='Exact search term for gene name')
parser.add_argument('-p', '--prefix', required=True, help='Prefix for output files')
args = parser.parse_args()

# Create output directory
outdir = f'{args.prefix}_out'
os.makedirs(outdir, exist_ok=True)

# Output filenames
id_file = os.path.join(outdir, f"{args.prefix}_target_ids.txt")
desc_file = os.path.join(outdir, f"{args.prefix}_target_description.tsv")

mrna_ids = []
descriptions = []
partial_matches = set()

with open(args.input, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        
        fields = line.strip().split('\t')
        if len(fields) >= 9 and fields[2] == 'mRNA':
            # Extract gene value from attributes
            if 'gene=' in fields[8]:
                gene_match = re.search(r'gene=([^;]+)', fields[8])
                if gene_match:
                    gene_value = gene_match.group(1)
                    
                    # Exact match check
                    if gene_value.lower() == args.search.lower():
                        # Extract ID and clean rna- prefix
                        id_part = fields[8].split('ID=')[1].split(';')[0]
                        clean_id = id_part.replace('rna-', '')
                        mrna_ids.append(clean_id)
                        descriptions.append(line.strip())
                    
                    # Partial match check (for suggestions)
                    elif args.search.lower() in gene_value.lower():
                        partial_matches.add(gene_value)

# Check results
if mrna_ids:
    # Success: exact matches found
    with open(id_file, 'w') as f:
        f.write('\n'.join(mrna_ids) + '\n')

    with open(desc_file, 'w') as f:
        f.write('#seqid\tsource\ttype\tstart\tend\tscore\tstrand\tphase\tattributes\n')
        f.write('\n'.join(descriptions) + '\n')

    print(f"Found {len(mrna_ids)} mRNA entries -> {id_file}, {desc_file}")

else:
    # No exact matches found
    if partial_matches:
        print(f"Error: No exact match found for '{args.search}'")
        print("Did you mean one of these?")
        for match in sorted(partial_matches):
            print(f"  - {match}")
    else:
        print(f"Error: No mRNA entries found containing '{args.search}'")
    
    exit(1)