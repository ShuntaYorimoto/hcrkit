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

### Developed and tested with:
- Python 3.13.3
- BLAST+ 2.16.0
- BioPython 1.85
- pandas 2.3.2

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
hcrkit --help
extract_target_ids --help
blastn -help
```

## Usage
### Basic Usage
```
# Usage
hcrkit -i INPUT_FASTA -d DATABASE -p PREFIX --initiator_id INITIATOR_ID [OPTIONS]

# Example
hcrkit -i ApDll.fasta -d Apisum_transcripts.fasta -p ApDll --initiator_id S73
```

#### Options
Use `hcrkit -h` or `hcrkit --help` for full parameter list.
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
extract_target_ids -i human_annotation.gff3 -s "DDX4" -p HsDDX4

# Use extracted IDs in probe design
hcrkit -i HsDDX4.fasta -d human_transcriptome.fasta -p HsDDX4 \
  --initiator_id A161 --target_ids HsDDX4_out/HsDDX4_target_ids.txt
```

2. Custom GC Content Ranges
Design probes with different GC stringency:

```
# Relaxed GC range (more probe candidates)
hcrkit -i gene.fasta -d db.fasta -p test --initiator_id S23 \
  --min_gc 40 --max_gc 60
```

3. Custom Initiator Sequences
Create a CSV file with custom initiators (no header):
```
B1,GAGGAGGGCAGCAAACGG
B2,CCTCGTAAATCCTCATCA
B3,GTCCCTGCCTCTATATCT
```
_Example sequences from Choi et al., 2018_

Use custom initiators:
```
hcrkit -i gene.fasta -d db.fasta -p test \
  --initiator_id B1 --initiator_custom my_initiators.csv
```

## Output Files
All outputs are saved in {prefix}_out/ directory:

| File | Description |
|------|-------------|
| `{prefix}_{initiator}_probe_pairs_gc{min}-{max}.csv` | Final probe sequences for ordering |
| `{prefix}_{initiator}_probe_summary_gc{min}-{max}.txt` |	Summary with off-target coverage |
| `{prefix}_probe_candidates_gc{min}-{max}.fasta` |	Initial probe candidates |
| `{prefix}_blast_results_gc{min}-{max}.tsv` |	BLAST search results |

## Ciation

