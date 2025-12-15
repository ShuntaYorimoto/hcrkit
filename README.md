# hcrkit<!-- omit in toc -->

hcrkit is an automated probe design tool for EC-isHCR.</br>
**For beginners, please check [Overview](OVERVIEW.md).**

## Table of Contents<!-- omit in toc -->

- [Build hcrkit environment](#build-hcrkit-environment)
  - [Dependencies](#dependencies)
  - [Installation](#installation)
  - [File requirements](#file-requirements)
- [Running hcrkit](#running-hcrkit)
  - [Syntax](#syntax)
  - [Quick Examples](#quick-examples)
  - [Parameters](#parameters)
- [Output Files](#output-files)
- [Tips](#tips)
- [Citation](#citation)

## Build `hcrkit` environment

### Dependencies

- Python 3.7+
- BLAST+
- Python packages: BioPython, pandas

**Developed and tested with**:
- Python 3.13.3
- BLAST+ 2.16.0
- BioPython 1.85, pandas 2.3.2

### Installation

1. Prerequisites

    It is recommended to prepare environment for using `conda`(https://conda-forge.org/download/).</br>
    If you use a Windows machine, Windows Subsystem for Linux and Ubuntu are also required.

2. Clone the repository

    ```bash
    git clone https://github.com/ShuntaYorimoto/hcrkit.git
    cd hcrkit
    ```

3. Create and activate conda environment

    ```bash
    conda env create -f environment.yml
    conda activate hcrkit
    ```

4. Verify installation

    To confirm a successful installation, run the command below and ensure that the help messages are displayed.

    ```bash
    hcrkit.py -h
    blastn -help
    ```

### File requirements
`hcrkit` uses `BLAST` against a reference transcript sequences to remove probes that may hybridize with off-target transcripts, ensuring probe specificity.

#### Required files
1. Target transcript sequence (FASTA format)
2. Reference transcript sequences (FASTA format)
    - From NCBI or your own assembly
3. BLAST database (created from the FASTA file)

#### Optional files (for genes with multiple isoforms)
4. GFF3 annotation (NCBI-format GFF3 only)
    - Required when: Your target gene has multiple isoforms AND you want `hcrkit` to automatically identify all on-target transcript IDs.

#### Download reference and annotation files from NCBI

> [!Note]
> New to NCBI downloads?</br>
> See our [NCBI Download Guide](docs/ncbi_download_guide.md) for step-by-step instructions with screenshots.

1. Go to NCBI Genome database and search for your organism
2. Download `GCF_*_genomic.gff.gz` and `GCF_*_rna.fna.gz` files from FTP server
3. Decompress downloaded files

#### Create BLAST database
Before running `hcrkit`, create a BLAST database from the reference RNA sequences:

Syntax:
```bash
makeblastdb -in <PATH> -dbtype nucl -out <STR>
```

Example:
```bash
makeblastdb -in GCF_005508785.2_pea_aphid_22Mar2018_4r6ur_v2_rna.fna -dbtype nucl -out Apis
```

## Running hcrkit

### Syntax

```bash
hcrkit.py [-h] -i <PATH> -d <PATH> -p <STR> --initiator_id <ID> [--target_ids <PATH> | --gff3 <PATH>] [--gene_name <STR>] [--initiator_custom <PATH>] [--initiator_split <INT>] [--min_gc <FLOAT>] [--max_gc <FLOAT>] [-t <INT>] [-v]
```

### Quick Examples
Basic usage
```bash
hcrkit.py \
  -i ApDll.fasta \
  -d Apis \
  -p ApDll \
  --initiator_id A161
```

With GFF3 for target ID extraction
```bash
hcrkit.py \
  -i ApNos1.fasta \
  -d Apis \
  -p ApNos1 \
  --initiator_id A161 \
  --gff3 GCF_005508785.2_pea_aphid_22Mar2018_4r6ur_v2_genomic.gff \
  --gene_name LOC100216490
```

### Parameters

#### Required Parameters

- `-i, --input PATH`: Path to FASTA file of the target RNA sequence

- `-d, --database PATH`: Path to BLAST database
  - Specify the database name without extension (e.g., `Apis`)
  - Users should create a BLAST database before running `hcrkit`

- `-p, --prefix STR`: String used as a prefix for output files and a directory (e.g. if `-p ApNos1`, then file: `ApNos1_[...].csv`, dir:`ApNos1/`)

- `--initiator_id ID`: The name of fluorescently-labeled hairpin DNA

  - Specify the hairpin DNAs of ISHpalette™ Short hairpin amplifier ([Nepagene](https://nepagene.jp/en/products/fluorescent-tissue-staining-in-situ-hcr/ishpalette-2)): S23, S41, S45, S72, S73, A161
  - Users also specify hairpin DNAs of Molecular Instruments. In this case, `--initiator_custom` parameter should also be used.

#### Optional Parameters

- `--gff3 PATH`: Path to GFF3 file for automatic target ID extraction

  - Requires `--gene_name` parameter

- `--gene_name STR`: Gene name to search in GFF3 file

- `--target_ids PATH`: Path to file of a custom list of on-target transctipts IDs

- `--min_gc FLOAT`: Minimum GC percentage (default: 45.0)
- `--max_gc FLOAT`: Maximum GC percentage (default: 55.0)

  - The GC content can be changed with 40–60%

- `--initiator_custom`: Path to CSV file of custom initiator sequences

  - CSV format: initiator_name,sequence (no header)
  - Example:
    ```csv
    B1,GAGGAGGGCAGCAAACGG
    B2,CCTCGTAAATCCTCATCA
    B3,GTCCCTGCCTCTATATCT
    ```

- `--initiator_split INT`: Split position of initiator sequence (default: 9)

  - The initiator sequence is split, and each split portion is conjugated to the P1 and P2 probes.
  - Specify the split position from 5′ end (default: 9 for ISHpalette™).

- `-t, --threads`: Number of CPU threads for BLAST (default: 1)

- `-h, --help`: Print help information

- `-v, --version`: Print version information

## Output Files

hcrkit creates the outputs as a tree structure:

```
{prefix}_out/
├── {prefix}_{initiator}_probe_pairs_gc{min}-{max}.csv
├── {prefix}_{initiator}_probe_summary_gc{min}-{max}.txt`
└── temp/
    ├── {prefix}_region_candidates_gc{min}-{max}.fasta
    ├── {prefix}_blast_results_gc{min}-{max}.tsv
    ├── {prefix}_target_ids.txt (if using --gff3)
    └── {prefix}_target_description.tsv (if using --gff3)
```

- `{prefix}`: string specified by `-p` parameter
- `{initiator}`: initiator ID
- `{min}`, `{max}`: GC content range

Example: Running with `-p ApNos1 --initiator_id A161` creates:
```
ApNos1_out/
├── ApNos1_A161_probe_pairs_gc45-55.csv
├── ApNos1_A161_probe_summary_gc45-55.txt
└── temp/
...
```

**Details of main output files**</br>
- `{prefix}_{initiator}_probe_pairs_gc{min}-{max}.csv`:

  - Probe sequences formatted for ordering (e.g., Eurofins Genomics)
  - Example:
  
    ```csv
    Oligoname,Sequence
    ApNos1_A161_s38_P1,GGTACGCGAaaAAGGGTCCTCTTGAAAACACCCACT
    ApNos1_A161_s38_P2,GGCCACTCATTGTGATATCTTGGCAaaAGGTAGGTGTAA
    ApNos1_A161_s363_P1,GGTACGCGAaaAACTCGGCCAGTGTGTAGTAGTTGG
    ApNos1_A161_s363_P2,TGTTGCATGCGGTTCATGCGCAACAaaAGGTAGGTGTAA
    ...
    ```

- `{prefix}_{initiator}_probe_summary_gc{min}-{max}.txt`:

  - Summary information for each probe set
  - `Max off-target coverage (%)` indicates the maximum percentage of bases matching off-target transcripts
  - Example:
  
    ```
    Probe set name    P1                                    P2                                       Max off-target coverage (%)
    ApNos1_A161_s38   GGTACGCGAaaAAGGGTCCTCTTGAAAACACCCACT  GGCCACTCATTGTGATATCTTGGCAaaAGGTAGGTGTAA  0
    ApNos1_A161_s363  GGTACGCGAaaAACTCGGCCAGTGTGTAGTAGTTGG  TGTTGCATGCGGTTCATGCGCAACAaaAGGTAGGTGTAA  0
    ...
    ```

**Details of intermediate files (in `temp` directory)**</br>

- `{prefix}_region_candidates_gc{min}-{max}.fasta`:

  - All candidates of probe regions before BLAST evaluation
  - Example:

    ```
    >XM_003242475.4_region_38_89
    AGTGGGTGTTTTCAAGAGGACCCTTATTGCCAAGATATCACAATGAGTGGCC
    >XM_003242475.4_region_39_90
    GTGGGTGTTTTCAAGAGGACCCTTATTGCCAAGATATCACAATGAGTGGCCG
    ...
    ```

- `{prefix}_blast_results_gc{min}-{max}.tsv`:

  - Raw BLAST results in tabular output format 6
  - Example:
 
    ```
    XM_003242475.4_region_38_89	XM_003242475.4	100.000	52	0	0	1	52	38	89	1.68e-22	103
    XM_003242475.4_region_38_89	XM_016802088.2	100.000	21	0	0	32	52	231	251	5.30e-04	42.1
    ...
    ```

- `{prefix}_target_ids.txt` (when using `--gff3` and `--gene_name`):

  - List of on-target transcript IDs
  - Example:
  
    ```
    XM_016802089.2
    XM_016802088.2
    XM_003242475.4
    ```

- `{prefix}_target_description.tsv` (when using `--gff3` and `--gene_name`):

  - List of on-target transcript IDs and the description

## Tips

Although several mRNAs were detected with 2 probe sets, using 10 or more probe sets was standard for most mRNAs.</br>

**If only a few probe sets are obtained...**

- Relax GC content requirements from default (45–55%) to 40–60%:
  
  ```bash
  hcrkit -i target.fasta -d database -p output --initiator_id A161 --min_gc 40 --max_gc 60
  ```

- Since missing on-target transcript IDs cause incorrect filtering, review the on-target list:

    - When using `--gff3` and `--gene_name`: Check the list of on-target IDs (`{prefix}_target_ids.txt`)
    - When using `--target_ids`: Verify that the file contains all on-target transcript IDs

## Citing `hcrkit`

If you use `hcrkit` in your research, please cite the following paper:
> Manuscript in preparation
