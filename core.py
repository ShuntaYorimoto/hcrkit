#!/usr/bin/env python3
"""
HCR Probe Design - Core Functions
"""

import os
import re
import subprocess
from typing import Dict, List, Set, Tuple, Optional
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd

# Predefined HCR initiator sequences
PREDEFINED_INITIATORS = {
    'S23': 'GGGTGGTCGTCGAAGTCGTAT',
    'S41': 'GCTCGACGTTCCTTTGCAACA', 
    'S45': 'CCTCCACGTTCCATCTAAGCT',
    'S72': 'CGGTGGAGTGGCAAGTAGGAT',
    'S73': 'CGGTCAGGTGGCTAGTATGGA',
    'A161': 'GGTACGCGAAGGTAGGTGTAA'
}


class HCRProbeDesigner:
    """HCR probe design pipeline"""
    
    def __init__(self, input_fasta: str, prefix: str, outdir: str):
        """Initialize designer with configuration"""
        self.input_fasta = input_fasta
        self.prefix = prefix
        self.outdir = outdir
        self.temp_dir = os.path.join(outdir, 'temp')
        
        os.makedirs(self.outdir, exist_ok=True)
        os.makedirs(self.temp_dir, exist_ok=True)
        
        # Data storage
        self.sequences: Dict[str, str] = {}
        self.regions: Dict[str, str] = {}
        self.target_ids: Set[str] = set()
        self.region_coverage: Dict[str, Dict] = {}
        self.initiators: Dict[str, str] = {}
        
        # Store GC parameters for file naming
        self.min_gc_used: Optional[float] = None
        self.max_gc_used: Optional[float] = None
    
    def load_sequences(self) -> None:
        """Load single sequence from FASTA file"""
        print(f"\nLoading sequence from {self.input_fasta}")
        
        records = list(SeqIO.parse(self.input_fasta, 'fasta'))
        
        if len(records) == 0:
            raise ValueError("Input FASTA file is empty")
        if len(records) > 1:
            raise ValueError(
                f"Found {len(records)} sequences. "
                "Only single-sequence FASTA files are supported."
            )
        
        record = records[0]
        self.sequences[record.id] = str(record.seq).upper()
        print(f"Loaded 1 sequence: {record.id} ({len(record.seq)} bp)")
    
    def load_target_ids_from_file(self, target_ids_file: str) -> None:
        """Load target IDs from file"""
        with open(target_ids_file, 'r') as f:
            self.target_ids = set(line.strip() for line in f if line.strip())
        print(f"\nLoaded {len(self.target_ids)} target IDs")
    
    def load_target_ids_from_sequences(self) -> None:
        """Use loaded sequence ID as target"""
        self.target_ids = set(self.sequences.keys())
        seq_id = list(self.target_ids)[0]
        print(f"\nUsing sequence ID as on-target: {seq_id}")
    
    def extract_target_ids_from_gff3(self, gff3_file: str, gene_name: str) -> bool:
        """
        Extract target mRNA IDs from GFF3 by gene name
        
        Searches for exact gene name matches (case-insensitive).
        Also reports similar gene names for user awareness.
        """
        mrna_ids = []
        descriptions = []
        exact_matches = set()
        partial_matches = set()
        
        with open(gff3_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 9 or fields[2] != 'mRNA':
                    continue
                
                gene_match = re.search(r'gene=([^;]+)', fields[8])
                if not gene_match:
                    continue
                
                gene_value = gene_match.group(1)
                search_lower = gene_name.lower()
                gene_lower = gene_value.lower()
                
                if gene_lower == search_lower:
                    id_part = fields[8].split('ID=')[1].split(';')[0]
                    clean_id = id_part.replace('rna-', '')
                    mrna_ids.append(clean_id)
                    descriptions.append(line.strip())
                    exact_matches.add(gene_value)
                elif search_lower in gene_lower:
                    partial_matches.add(gene_value)
        
        partial_matches -= exact_matches
        
        if mrna_ids:
            self.target_ids = set(mrna_ids)
            
            # Write to temp directory
            id_file = os.path.join(self.temp_dir, f"{self.prefix}_target_ids.txt")
            desc_file = os.path.join(self.temp_dir, f"{self.prefix}_target_description.tsv")
            
            with open(id_file, 'w') as f:
                f.write('\n'.join(mrna_ids) + '\n')
            
            with open(desc_file, 'w') as f:
                f.write('#seqid\tsource\ttype\tstart\tend\tscore\tstrand\tphase\tattributes\n')
                f.write('\n'.join(descriptions) + '\n')
            
            print(f"\nFound {len(mrna_ids)} mRNA entries: {', '.join(sorted(exact_matches))}")
            
            # Show similar genes to avoid missing targets
            if partial_matches:
                print(f"\nNote: Found {len(partial_matches)} similar gene(s) (not included):")
                for match in sorted(partial_matches):
                    print(f"  - {match}")
            
            return True
        
        # No exact matches found
        if partial_matches:
            print(f"Error: No exact match for '{gene_name}'")
            print(f"Did you mean: {', '.join(sorted(partial_matches))}?")
        else:
            print(f"Error: No genes found matching '{gene_name}'")
        
        return False
    
    def calculate_gc_content(self, sequence: str) -> float:
        """Calculate GC percentage"""
        gc_count = sequence.count('G') + sequence.count('C')
        return (gc_count / len(sequence)) * 100
    
    def _filter_duplicates(self, region_dict: Dict[str, str]) -> Dict[str, str]:
        """Remove duplicate sequences, keeping first occurrence"""
        seen_sequences = set()
        unique_regions = {}
        for region_id, sequence in region_dict.items():
            if sequence not in seen_sequences:
                seen_sequences.add(sequence)
                unique_regions[region_id] = sequence
        return unique_regions
    
    def generate_probe_candidates(self, min_gc: float = 45.0, max_gc: float = 55.0) -> str:
        """
        Generate probe region candidates with GC filtering
        
        Uses sliding window (52 nt) to generate candidates.
        Each region is split into P1 (25 nt) and P2 (25 nt) with 2 nt spacer.
        Both P1 and P2 must meet GC criteria.
        """
        self.min_gc_used = min_gc
        self.max_gc_used = max_gc
        
        print(f"\nGenerating probe region candidates (GC: {min_gc}-{max_gc}%)")
        
        gc_suffix = f"gc{int(min_gc)}-{int(max_gc)}"
        outfile = os.path.join(self.temp_dir, f"{self.prefix}_region_candidates_{gc_suffix}.fasta")
        
        gc_filtered_regions = {}
        total_count = 0
        
        for seq_id, sequence in self.sequences.items():
            for i in range(len(sequence) - 51):
                region_seq = sequence[i:i + 52]
                total_count += 1
                
                # Split: P1 (0-24), spacer (25-26), P2 (27-51)
                p1_seq = region_seq[:25]
                p2_seq = region_seq[27:]
                
                if (min_gc <= self.calculate_gc_content(p1_seq) <= max_gc and 
                    min_gc <= self.calculate_gc_content(p2_seq) <= max_gc):
                    region_id = f"{seq_id}_region_{i+1}_{i+52}"
                    gc_filtered_regions[region_id] = region_seq
        
        self.regions = self._filter_duplicates(gc_filtered_regions)
        
        with open(outfile, 'w') as f:
            for region_id, region_seq in self.regions.items():
                f.write(f">{region_id}\n{region_seq}\n")
        
        print(f"Total windows: {total_count}, After GC filter: {len(gc_filtered_regions)}, Unique: {len(self.regions)}")
        return outfile
    
    def run_blast_search(self, database: str, region_fasta: str, threads: int = 1) -> str:
        """Run BLASTN search against database"""
        print(f"\nRunning BLASTN (threads: {threads})")
        
        gc_suffix = f"gc{int(self.min_gc_used)}-{int(self.max_gc_used)}"
        blast_out = os.path.join(self.temp_dir, f"{self.prefix}_blast_results_{gc_suffix}.tsv")
        
        cmd = [
            'blastn', '-query', region_fasta, '-db', database,
            '-task', 'blastn-short', '-outfmt', '6', '-evalue', '1.0e-2',
            '-num_threads', str(threads), '-out', blast_out
        ]
        
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(f"BLAST completed -> {blast_out}")
        return blast_out
    
    def parse_blast_results(self, blast_file: str) -> None:
        """
        Parse BLAST results and filter by specificity
        
        Marks regions as invalid if they have ≥50% coverage to off-target sequences.
        Tracks maximum off-target coverage for each region.
        """
        print(f"\nParsing BLAST results")
        
        self.region_coverage = {
            region_id: {'max_offtarget': 0.0, 'valid': True}
            for region_id in self.regions.keys()
        }
        
        with open(blast_file, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                query_id = fields[0]
                subject_id = fields[1]
                length = int(fields[3])
                coverage = length / 52
                
                is_ontarget = any(target_id in subject_id for target_id in self.target_ids)
                
                if not is_ontarget:
                    current_max = self.region_coverage[query_id]['max_offtarget']
                    self.region_coverage[query_id]['max_offtarget'] = max(
                        current_max, coverage * 100
                    )
                    if coverage >= 0.5:
                        self.region_coverage[query_id]['valid'] = False
        
        valid_count = sum(1 for data in self.region_coverage.values() if data['valid'])
        print(f"Valid probe regions: {valid_count}")
    
    def select_non_overlapping_probes(self) -> List[str]:
        """
        Select non-overlapping probe regions using greedy algorithm
        
        Sorts regions by position and selects the first available non-overlapping region.
        """
        print(f"\nSelecting non-overlapping probe regions")
        
        valid_regions = [
            region_id for region_id, data in self.region_coverage.items() 
            if data['valid']
        ]
        
        # Extract positions and sort by start position
        region_list = []
        for region_id in valid_regions:
            parts = region_id.split('_')
            start = int(parts[-2])
            end = int(parts[-1])
            region_list.append((start, end, region_id))
        region_list.sort()
        
        # Greedy selection: pick first available non-overlapping region
        selected_regions = []
        last_end = 0
        for start, end, region_id in region_list:
            if start > last_end:
                selected_regions.append(region_id)
                last_end = end
        
        print(f"Selected {len(selected_regions)} non-overlapping probe regions")
        return selected_regions
    
    def load_initiators(self, initiator_id: str, 
                       custom_file: Optional[str] = None) -> str:
        """Load initiator sequences from predefined set or custom file"""
        if custom_file:
            df = pd.read_csv(custom_file, header=None, names=['ID', 'Sequence'])
            self.initiators = dict(zip(df['ID'], df['Sequence']))
            print(f"\nUsing custom initiators from: {custom_file}")
        else:
            self.initiators = PREDEFINED_INITIATORS
            print("\nUsing predefined initiators")
        
        if initiator_id not in self.initiators:
            available = ', '.join(self.initiators.keys())
            raise ValueError(f"Initiator '{initiator_id}' not found. Available: {available}")
        
        initiator_seq = self.initiators[initiator_id]
        print(f"Selected: {initiator_id} ({initiator_seq})")
        return initiator_seq
    
    def generate_probe_pairs(self, final_region_ids: List[str], 
                            initiator_id: str, initiator_seq: str,
                            split_pos: int = 9) -> Tuple[pd.DataFrame, List[List[str]]]:
        """
        Generate probe pairs with initiator sequences
        
        Creates split probes (P1 and P2) by:
        1. Taking reverse complement of P1 and P2 regions
        2. Adding initiator sequences (split at specified position)
        """
        print(f"\nGenerating probe pairs")
        
        if split_pos >= len(initiator_seq):
            raise ValueError(
                f"Split position ({split_pos}) >= initiator length ({len(initiator_seq)})"
            )
        
        # Split initiator and add spacer
        p1_init = initiator_seq[:split_pos] + "aa"
        p2_init = "aa" + initiator_seq[split_pos:]
        
        probe_pairs = []
        summary_data = []
        
        for region_id in final_region_ids:
            region_seq = self.regions[region_id]
            start_pos = int(region_id.split('_')[-2])
            
            # Extract P1 and P2 regions
            p1_region = region_seq[0:25]
            p2_region = region_seq[27:52]
            
            # Create probe sequences: initiator + reverse_complement(region)
            p1_seq = p1_init + str(Seq(p1_region).reverse_complement())
            p2_seq = str(Seq(p2_region).reverse_complement()) + p2_init
            
            probe_set_name = f"{self.prefix}_{initiator_id}_s{start_pos}"
            
            probe_pairs.extend([
                [f"{probe_set_name}_P1", p1_seq],
                [f"{probe_set_name}_P2", p2_seq]
            ])
            
            max_coverage = self.region_coverage.get(region_id, {}).get('max_offtarget', 0.0)
            summary_data.append([
                probe_set_name, p1_seq, p2_seq, f"{max_coverage:.0f}"
            ])
        
        probe_pairs_df = pd.DataFrame(probe_pairs, columns=['Oligoname', 'Sequence'])
        print(f"Generated {len(probe_pairs)} probes ({len(summary_data)} pairs)")
        
        return probe_pairs_df, summary_data
    
    def write_outputs(self, probe_pairs_df: pd.DataFrame, 
                     summary_data: List[List[str]], 
                     initiator_id: str) -> Tuple[str, str]:
        """Write probe pairs and summary files"""
        gc_suffix = f"gc{int(self.min_gc_used)}-{int(self.max_gc_used)}"
        
        csv_out = os.path.join(self.outdir, 
                               f"{self.prefix}_{initiator_id}_probe_pairs_{gc_suffix}.csv")
        probe_pairs_df.to_csv(csv_out, index=False)
        
        summary_out = os.path.join(self.outdir, 
                                   f"{self.prefix}_{initiator_id}_probe_summary_{gc_suffix}.txt")
        with open(summary_out, 'w') as f:
            f.write("Probe set name\tP1\tP2\tMax off-target coverage (%)\n")
            for row in summary_data:
                f.write('\t'.join(row) + '\n')
        
        print(f"\nOutputs:\n  {csv_out}\n  {summary_out}")
        return csv_out, summary_out