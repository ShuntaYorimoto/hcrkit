# hcrkit
Automated pipeline for HCR (Hybridization Chain Reaction) probe design with customizable parameters

## Features

- **Automated pipeline** from sequence input to final probe pairs
- **Flexible GC content filtering** with customizable range (default 45-55%)
- **BLAST-based specificity screening** to minimize off-target hybridization
- **Multiple initiator support** (predefined Nepagene sequences + custom CSV)
- **Isoform detection** from GFF3 annotations for comprehensive on-target identification

## Installation

### Requirements
- Python 3.7+
- BLAST+
- Python packages: biopython, pandas

**Developed and tested with:**
- Python 3.13.3, BLAST+ 2.16.0, BioPython 1.85, pandas 2.3.2

### Setup
1. Clone the repository:
```
git clone https://github.com/ShuntaYorimoto/hcrkit.git
cd hcrkit
```
2. Create and activate conda environment:
```
conda create -n hcrkit python=3.13 blast -c conda-forge -c bioconda
conda activate hcrkit
```

3. Install hcrkit and Python dependencies:
```
pip install -e .
```

4. Verify installation:
```
hcrkit.py -h
extract_target_ids.py -h
blastn -help
```

## Usage
### Quick Usage
```
hcrkit.py -i gene.fasta -d transcripts.fasta -p MyGene --initiator_id S73 --make_db
```

### Command Syntax

```
hcrkit.py -i INPUT_FASTA -d DATABASE -p PREFIX --initiator_id INITIATOR_ID [OPTIONS]
```

### Parameters
Use `hcrkit.py -h` or `hcrkit.py --help` for full parameter list.
| Parameter | Description | Default |
|-----------|-------------|---------|
| `-i, --input` | Target gene FASTA file | Required |
| `-d, --database` | Transcriptome FASTA or BLAST database | Required |
| `-p, --prefix` | Output file prefix | Required |
| `--initiator_id` | HCR initiator sequence ID | Required |
| `--min_gc` | Minimum GC percentage | 45 |
| `--max_gc` | Maximum GC percentage | 55 |
| `-t, --threads` |	CPU threads for BLAST |	1 |
| `--make_db` |	Create BLAST database from FASTA | Flag |
| `--target_ids` | On-target sequence ID list	| Optional |
| `--initiator_custom` | Custom initiator CSV file | Optional |

### Advanced Usage
1. Extract Target IDs from GFF3  
For genes with multiple isoforms, first extract target sequence IDs by gene name:

```
# Extract all DDX4 isoforms
extract_target_ids.py -i human_annotation.gff3 -s "DDX4" -p HsDDX4

# Use extracted IDs in probe design
hcrkit.py -i HsDDX4.fasta -d human_transcriptome.fasta -p HsDDX4 \
  --initiator_id A161 --target_ids HsDDX4_out/HsDDX4_target_ids.txt
```

2. Custom GC Content Ranges  
Design probes with different GC stringency:

```
# Relaxed GC range (more probe candidates)
hcrkit.py -i gene.fasta -d db.fasta -p test --initiator_id S23 \
  --min_gc 40 --max_gc 60
```

3. Custom Initiator Sequences  
Create a CSV file with custom initiators (no header) (_Example sequences from Choi et al., 2018_):  
```
B1,GAGGAGGGCAGCAAACGG
B2,CCTCGTAAATCCTCATCA
B3,GTCCCTGCCTCTATATCT
```

```
hcrkit.py -i gene.fasta -d db.fasta -p test \
  --initiator_id B1 --initiator_custom my_initiators.csv
```

## Algorithm and Implementation Details
### Workflow Overview
HCRKit implements a 3-step automated pipeline:

1. GC Content Filtering (`step01_filter_gc_content.py`)
- Generates 52bp probe candidates using sliding window (1bp step)
- Filters candidates based on GC content of both P1 (first 25bp) and P2 (last 25bp) regions
- Default GC range: 45-55% (customizable via --min_gc and --max_gc)

2. BLAST-based Specificity Screening (`step02_blast_search.py`)
- Performs BLASTN search against provided transcriptome database
- BLAST parameters: -task blastn-short -evalue 1.0e-2 -outfmt 6
- Optimized for short sequence alignment with relaxed e-value threshold

3. Probe Selection and Pair Generation (`step03_filter_and_generate.py`)
- Filters probes with ≥50% coverage against off-target sequences
- Selects non-overlapping probes using greedy algorithm
- Removes duplicate sequences from repetitive regions
- Generates final HCR probe pairs with specified initiator sequences

### On-target Sequence Handling
**Automatic on-target detection:** When `--target_ids` is not specified, hcrkit automatically uses the input FASTA sequence ID as on-target references. This is suitable for:
- Single-isoform gene
- Pre-selected isoform gene with non-redundant database (e.g., longest isoform representatives)

**Manual on-target specification:** Use `extract_target_ids.py` to identify all isoforms of a target gene from from GFF3 annotations. This is essential when:
- Working with NCBI annotations where isoforms have unrelated mRNA IDs (e.g., XM_001234567, XM_007654321 for the same gene)
- Gene name-based isoform discovery is needed to ensure all variants are recognized as on-target
- Comprehensive isoform coverage is required for accurate specificity filtering

### BLAST Specificity Criteria
- **Coverage threshold:** Probes showing ≥50% alignment coverage to off-target sequences are excluded (≥26bp out of 52bp total length)
- **On-target hits:** Alignments to sequences specified in --target_ids (or input FASTA ID) are allowed regardless of coverage
- **E-value:** Uses relaxed threshold (1.0e-2) to detect potential cross-hybridization

## Output Files
All outputs are saved in `{prefix}_out/` directory:

| File | Description |
|------|-------------|
| `{prefix}_{initiator}_probe_pairs_gc{min}-{max}.csv` | Ready-to-order sequences |
| `{prefix}_{initiator}_probe_summary_gc{min}-{max}.txt` |	Summary includes max off-target coverage for each probe set |
| `{prefix}_probe_candidates_gc{min}-{max}.fasta` |	Initial probe candidates |
| `{prefix}_blast_results_gc{min}-{max}.tsv` |	BLAST search results |

## Predefined Initiators

| ID | Sequence | Source |
|----|----------|---------|
| S23 | GGGTGGTCGTCGAAGTCGTAT | Nepagene |
| S41 | GCTCGACGTTCCTTTGCAACA | Nepagene |
| S45 | CCTCCACGTTCCATCTAAGCT | Nepagene |
| S72 | CGGTGGAGTGGCAAGTAGGAT | Nepagene |
| S73 | CGGTCAGGTGGCTAGTATGGA | Nepagene |
| A161 | GGTACGCGAAGGTAGGTGTAA | Nepagene |

## Citation
Preparing...
