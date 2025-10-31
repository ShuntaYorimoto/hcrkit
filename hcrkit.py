#!/usr/bin/env python3
"""
HCR Probe Design Pipeline
"""

__version__ = "2.0.0"

import argparse
import sys
from core import HCRProbeDesigner


def main():
    parser = argparse.ArgumentParser(description='HCR Probe Design Pipeline', add_help=False)
    
    parser.add_argument('-h', '--help', action='help', help='Print help')
    parser.add_argument('-i', '--input', required=True, metavar='PATH', help='Input FASTA file')
    parser.add_argument('-d', '--database', required=True, metavar='PATH', help='BLAST database')
    parser.add_argument('-p', '--prefix', required=True, metavar='STR', help='Output prefix')
    parser.add_argument('--initiator_id', required=True, metavar='ID',
                       help='Initiator ID (S23, S41, S45, S72, S73, A161)')
    
    target_group = parser.add_mutually_exclusive_group()
    target_group.add_argument('--target_ids', metavar='PATH', help='Target ID list file')
    target_group.add_argument('--gff3', metavar='PATH', help='GFF3 file (requires --gene_name)')
    
    parser.add_argument('--gene_name', metavar='STR', help='Gene name for GFF3 extraction')
    parser.add_argument('--initiator_custom', metavar='PATH', help='Custom initiators CSV')
    parser.add_argument('--initiator_split', type=int, metavar='INT', 
                       help='Initiator split position (default: 9)')
    parser.add_argument('--min_gc', type=float, metavar='FLOAT', help='Min GC%% (default: 45.0)')
    parser.add_argument('--max_gc', type=float, metavar='FLOAT', help='Max GC%% (default: 55.0)')
    parser.add_argument('-t', '--threads', type=int, metavar='INT', help='BLAST threads (default: 1)')
    parser.add_argument('-v', '--version', action='version', version=f'hcrkit {__version__}')
    
    args = parser.parse_args()
    
    if args.gff3 and not args.gene_name:
        parser.error("--gff3 requires --gene_name")
    
    if args.min_gc is not None and args.max_gc is not None:
        if args.min_gc >= args.max_gc:
            parser.error("--min_gc must be less than --max_gc")
    
    outdir = f'{args.prefix}_out'
    
    print("=" * 60)
    print("HCR Probe Design Pipeline")
    print("=" * 60)
    
    try:
        designer = HCRProbeDesigner(args.input, args.prefix, outdir)
        
        print("\n" + "=" * 60)
        print("Step 1: Load Sequences")
        print("=" * 60)
        designer.load_sequences()
        
        print("\n" + "=" * 60)
        print("Step 2: Load Target IDs")
        print("=" * 60)
        
        if args.gff3:
            if not designer.extract_target_ids_from_gff3(args.gff3, args.gene_name):
                sys.exit(1)
        elif args.target_ids:
            designer.load_target_ids_from_file(args.target_ids)
        else:
            designer.load_target_ids_from_sequences()
        
        print("\n" + "=" * 60)
        print("Step 3: Generate Probe Region Candidates")
        print("=" * 60)
        
        gc_kwargs = {}
        if args.min_gc is not None:
            gc_kwargs['min_gc'] = args.min_gc
        if args.max_gc is not None:
            gc_kwargs['max_gc'] = args.max_gc
        
        region_fasta = designer.generate_probe_candidates(**gc_kwargs)
        
        print("\n" + "=" * 60)
        print("Step 4: BLAST Search")
        print("=" * 60)
        
        blast_kwargs = {'database': args.database, 'region_fasta': region_fasta}
        if args.threads is not None:
            blast_kwargs['threads'] = args.threads
        
        blast_results = designer.run_blast_search(**blast_kwargs)
        
        print("\n" + "=" * 60)
        print("Step 5: Filter by Specificity")
        print("=" * 60)
        designer.parse_blast_results(blast_results)
        
        print("\n" + "=" * 60)
        print("Step 6: Select Non-overlapping Probe Regions")
        print("=" * 60)
        final_region_ids = designer.select_non_overlapping_probes()
        
        if not final_region_ids:
            raise ValueError("No valid probe regions found")
        
        print("\n" + "=" * 60)
        print("Step 7: Generate Probe Pairs")
        print("=" * 60)
        
        initiator_seq = designer.load_initiators(args.initiator_id, args.initiator_custom)
        
        pairs_kwargs = {
            'final_region_ids': final_region_ids,
            'initiator_id': args.initiator_id,
            'initiator_seq': initiator_seq
        }
        if args.initiator_split is not None:
            pairs_kwargs['split_pos'] = args.initiator_split
        
        probe_pairs_df, summary_data = designer.generate_probe_pairs(**pairs_kwargs)
        
        print("\n" + "=" * 60)
        print("Step 8: Write Outputs")
        print("=" * 60)
        designer.write_outputs(probe_pairs_df, summary_data, args.initiator_id)
        
        print("\n" + "=" * 60)
        print("COMPLETED SUCCESSFULLY")
        print("=" * 60)
    
    except KeyboardInterrupt:
        print("\n\nInterrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\nERROR: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()