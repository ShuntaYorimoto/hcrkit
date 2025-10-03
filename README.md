# hcrkit
Automated pipeline for HCR (Hybridization Chain Reaction) probe design with customizable parameters

## Table of Contents
- [Overview](#overview)
- [Installation](#installation)
  - [Dependencies](#dependencies)
  - [Setup](#setup)
  - [Download transcript and gff files from NCBI (Optional)](#download-transcript-and-gff-files-from-ncbi-optional)
- [Usage](#usage)
  - [Preparation (Optional)](#preparation-optional)
  - [Basic Usage](#basic-usage)
  - [Parameters](#parameters)
  - [Output Files](#output-files)
- [Tips](#tips)
- [Citation](#citation)

## Overview
### What is isHCR (_in situ_ Hybridization Chain Reaction)?
isHCR is a powerful _in situ_ hybridization technique that allows visualization of RNA molecules in tissues with high sensitivity and specificity. Unlike traditional methods, HCR uses a cascade amplification system where probe binding triggers a chain reaction of DNA hybridization, dramatically enhancing signal strength.

### How does hcrkit work?
hcrkit automates the complex process of designing isHCR probe pairs by:

1. **Isoform detection**: Extracting isoform IDs from GFF3 annotations for comprehensive on-target identification
2. **Probe candidate generation**: Scanning your target sequence to identify potential 52-nucleotide probe regions
3. **Quality filtering**: Selecting probes based on GC content (45-55% recommended) and filtering duplicated sequences
4. **Specificity screening**: Using BLAST (-task blastn-short -evalue 1.0e-2) to eliminate probes that might bind to off-target transcripts (>50% coverage)
5. **Optimal selection**: Choosing non-overlapping probes for maximum coverage
6. **Split-probe generation**: Creating reverse-complemented probe pairs (P1/P2) with initiator sequences

## Installation

### Dependencies
- Python 3.7+
- BLAST+
- Python packages: biopython, pandas

**Developed and tested with:**
Python 3.13.3, BLAST+ 2.16.0, BioPython 1.85, pandas 2.3.2

### Setup
1. Clone the repository:
```bash
git clone https://github.com/ShuntaYorimoto/hcrkit.git
cd hcrkit
```
2. Create and activate conda environment:
```bash
conda env create -f environment.yml
conda activate hcrkit
```

3. Verify installation: <br>
Installation suceeded when the help guidance is shown for all commands.
```bash
hcrkit.py --help
extract_target_ids.py --help
blastn -help
```

### Download transcript and gff files from NCBI (Optional)
For comprehensive specificity screening, download the complete transcriptome of your organism:

1. Go to NCBI Genome database and search for your organism
2. Download `GCF_*_genomic.gff.gz` and `GCF_*_rna.fna.gz` files from ftp server
3. Decompress downloaded files

Example:
```bash
# Download using curl
curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/005/508/785/GCF_005508785.2_pea_aphid_22Mar2018_4r6ur_v2/GCF_005508785.2_pea_aphid_22Mar2018_4r6ur_v2_genomic.gff.gz
curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/005/508/785/GCF_005508785.2_pea_aphid_22Mar2018_4r6ur_v2/GCF_005508785.2_pea_aphid_22Mar2018_4r6ur_v2_rna.fna.gz

# Decompress files
gunzip *.gz
```

## Usage
### Preparation (Optional)
#### Prepare on-target isoform IDs list
If your gene has multiple isoforms, create a list of transcript IDs that should be considered as valid targets (not off-targets). The `extract_target_ids.py` script helps extract all isoforms of a specific gene from NCBI GFF3 files.

```bash
# Syntax
extract_target_ids.py -i FILE -s STRING -p STRING

# Example: Extract all isoforms of LOC100216490 (nanos-like protein) in pea aphid
extract_target_ids.py -i GCF_005508785.2_pea_aphid_22Mar2018_4r6ur_v2_genomic.gff -s LOC100216490 -p ApNos1
```

This creates `ApNos1_target_ids.txt` containing:
```
XM_016802089.2
XM_016802088.2
XM_003242475.4
```

#### Automatic on-target detection
When `--target_ids` is not specified, hcrkit automatically uses the input FASTA sequence ID as on-target reference. This is suitable for:
- Single-isoform gene
- Pre-selected isoform sequences with non-redundant databases (e.g., longest isoform representatives)

### Basic Usage
```bash
# Syntax
hcrkit.py -i FILE -d FILE -p STRING --initiator_id ID [OPTIONS]

# First run: Create BLAST database
hcrkit.py -i ApDll.fasta -d GCF_005508785.2_pea_aphid_22Mar2018_4r6ur_v2_rna.fna -p ApDll --initiator_id A161 --make_db Apis -t 4

# Subsequent runs: Use existing database
hcrkit.py -i ApVas1.fasta -d Apis_blastdb/Apis -p ApVas1 --initiator_id S73 -t 4

# With target IDs: Input a list of transcript IDs for multi-isogorm genes
hcrkit.py -i ApNos1.fasta -d Apis_blastdb/Apis -p ApNos1 --initiator_id A161 --target_ids ApNos1_out/ApNos1_target_ids.txt -t 4
```

### Parameters
#### Required Parameters
`-i, --input FILE`: Input transcript FASTA file containing your target sequence<br>

Example:  `ApNos1.fasta`
```fasta
>XM_003242475.4 PREDICTED: Acyrthosiphon pisum nanos-like protein (LOC100216490), transcript variant X3, mRNA
AAAAAAATGTATGTTTTGAATTTGAAGTTTAAATTTAAGTGGGTGTTTTCAAGAGGACCCTTATTGCCAAGATATCACAA
TGAGTGGCCGTAGCAGATTCGCGTTCGGCAGCCAGACCCAGGCGCCCCAGTACCAGAACCGATGCGCTCTGTGCTTGACC
...
```

`-d, --database FILE`: BLAST database or FASTA file
- Use your organism's complete transcriptome
- Can be a pre-built BLAST database or FASTA file

`-p, --prefix STRING`: Prefix for output files and directories
- Organizes outputs in `{prefix}_out/` directory to prevent file clutter

`--initiator_id ID`: HCR initiator ID
- Predefined options: S23, S41, S45, S72, S73, A161 (Nepagene sequences)
- Input different initiator IDs for custom initiators

#### Optional Parameters

`--target_ids FILE`: File listing transcript IDs that should be considered as valid targets
- Prevents isoforms of the same gene from being flagged as off-targets
- Essential for genes with multiple splice variants

`--make_db STRING`: Create BLAST database with specified name
- Database files are saved in `{STRING}_blastdb/` directory for reuse

`--min_gc`, `--max_gc`: GC content range (default: 45-55%)
- The GC content should be with 40–60% (45-55% is recommended)
- Default GC content range (45-55%) is the most strict condition

`--initiator_custom`: Custom initiators CSV file
- Example: `custom_initiators.csv`
```csv
B1,GAGGAGGGCAGCAAACGG
B2,CCTCGTAAATCCTCATCA
B3,GTCCCTGCCTCTATATCT
```

`--initiator_split INT`: Position to split initiator between P1 and P2 (default: 9)
- For custom initiators

`-t, --threads`: Number of CPU threads for BLAST (default: 1)

`-h, --help`: Print help information

`-v, --version`: Print version information

### Output Files
hcrkit creates a structured output directory with the following files:

#### Main Output Files
`{prefix}_out/{prefix}_{initiator}_probe_pairs_gc{min}-{max}.csv`
- Ready-to-order probe sequences in CSV format (eurofins)
- Example: `ApNos1_out/ApNos1_A161_probe_pairs_gc45-55.csv`
```
Oligoname,Sequence
ApNos1_A161_s38_P1,GGTACGCGAaaAAGGGTCCTCTTGAAAACACCCACT
ApNos1_A161_s38_P2,GGCCACTCATTGTGATATCTTGGCAaaAGGTAGGTGTAA
ApNos1_A161_s363_P1,GGTACGCGAaaAACTCGGCCAGTGTGTAGTAGTTGG
ApNos1_A161_s363_P2,TGTTGCATGCGGTTCATGCGCAACAaaAGGTAGGTGTAA
...
```

`{prefix}_out/{prefix}_{initiator}_probe_summary_gc{min}-{max}.txt`
- Probe summary with specificity information
- Example: `ApNos1_out/ApNos1_A161_probe_summary_gc45-55.txt`
```tsv
Probe set name    P1                                    P2                                       Max off-target coverage (%)
ApNos1_A161_s38   GGTACGCGAaaAAGGGTCCTCTTGAAAACACCCACT  GGCCACTCATTGTGATATCTTGGCAaaAGGTAGGTGTAA  0
ApNos1_A161_s363  GGTACGCGAaaAACTCGGCCAGTGTGTAGTAGTTGG  TGTTGCATGCGGTTCATGCGCAACAaaAGGTAGGTGTAA  0
...
```
- Max off-target coverage: Percentage of probe sequence that matches off-target transcripts.
- Probes with >50% coverage are automatically removed
- 0%: Highly specific probes (recommended)

#### Intermediate Files (in temp/ directory)
`{prefix}_out/temp/{prefix}_probe_candidates_gc{min}-{max}.fasta`
- All potential probe sequences before specificity filtering
- Example: `ApNos1_out/temp/ApNos1_probe_candidates_gc45-55.fasta`
```
>XM_003242475.4_probe_38_89
AGTGGGTGTTTTCAAGAGGACCCTTATTGCCAAGATATCACAATGAGTGGCC
>XM_003242475.4_probe_39_90
GTGGGTGTTTTCAAGAGGACCCTTATTGCCAAGATATCACAATGAGTGGCCG
...
```

`{prefix}_out/temp/{prefix}_blast_results_gc{min}-{max}.tsv`
- Raw BLASTN results (outfmt 6)
- Example: `ApNos1_out/temp/ApNos1_blast_results_gc45-55.tsv`
```
XM_003242475.4_probe_38_89	XM_003242475.4	100.000	52	0	0	1	52	38	89	1.68e-22	103
XM_003242475.4_probe_38_89	XM_016802088.2	100.000	21	0	0	32	52	231	251	5.30e-04	42.1
...
```

## Tips
### Appropreate Probe Number
- **Start with 10 probe pairs** as a baseline for most targets
- **If no signal is detected**, increase the number of probe pairs
- **Fewer pairs may be sufficient** for highly expressed genes (signals were detected with only 2 pairs in some cases)

### If Few probes are designed
- Relax GC content requirements:

```bash
# Default (strict: 45-55%)
hcrkit.py -i target.fasta -d database -p output --initiator_id A161

# Relaxed (40-60%)
hcrkit -i target.fasta -d database -p output --initiator_id A161 --min_gc 40 --max_gc 60
```

- Check the on-target list:
  - Verify that --target_ids file contains all relevant transcript IDs
  - Missing isoforms may be incorrectly flagged as off-targets

## Citation
Preparing...
