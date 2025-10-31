# hcrkit<!-- omit in toc -->

hcrkit is an automated pipeline for design of *in situ* Hybridization Chain Reaction (isHCR) probes. **For beginners, please check [Overview](#overview) section.**

## Table of Contents<!-- omit in toc -->

- [Installation](#installation)
  - [Dependencies](#dependencies)
  - [Setup](#setup)
- [Design isHCR probe](#design-ishcr-probe)
  - [Download reference files from NCBI](#download-reference-files-from-ncbi)
    - [Procedure for download files via NCBI web page](#procedure-for-download-files-via-ncbi-web-page)
    - [Procedure for download files via CLI](#procedure-for-download-files-via-cli)
  - [Prepare BLAST database](#prepare-blast-database)
  - [Running hcrkit](#running-hcrkit)
    - [Syntax](#syntax)
    - [Quick Examples](#quick-examples)
    - [Parameters](#parameters)
- [Tips](#tips)
- [Overview](#overview)
  - [What is isHCR (_in situ_ Hybridization Chain Reaction)?](#what-is-ishcr-in-situ-hybridization-chain-reaction)
  - [Guidline for probe design](#guidline-for-probe-design)
  - [Workflow \& algorithm](#workflow--algorithm)
- [Citation](#citation)

## Installation

### Dependencies

- Python 3.7+
- BLAST+
- Python packages: BioPython, pandas

**Developed and tested with**:
- Python 3.13.3
- BLAST+ 2.16.0
- BioPython 1.85, pandas 2.3.2

### Setup

1. Prepare for conda
   
   hcrkit can be easily set up using conda, so it is recommended to prepare environment for using conda. If you use a Windows machine, Windows Subsystem for Linux and Ubuntu are also required.

2. Clone the repository

   ```bash
   git clone https://github.com/ShuntaYorimoto/hcrkit.git
   cd hcrkit
   ```

4. Create and activate conda environment

   ```bash
   conda env create -f environment.yml
   conda activate hcrkit
   ```

5. Verify installation
   
   To verify successful installation, run the command below and check whether the help guidance is shown.

   ```bash
   hcrkit.py -h
   blastn -help
   ```

## Design isHCR probe

### Download reference files from NCBI

It is necessary to download a reference RNA sequence file and a gff3-formated annotation file. `hcrkit` automatically filters candidate probes that might bind off-target transcripts using BLAST (Details are shown in [Workflow \& algorithm](#workflow--algorithm) section). It is recommended to download files from NCBI to match transcript IDs between a reference RNA sequence file and a gff file.

#### Procedure for download files via NCBI web page

> [!Note]
> New to NCBI downloads? See our [NCBI Download Guide](docs/ncbi_download_guide.md) for step-by-step instructions with screenshots.

1. Go to NCBI Genome database and search for your organism
2. Download `GCF_*_genomic.gff.gz` and `GCF_*_rna.fna.gz` files from FTP server
3. Decompress downloaded files

#### Procedure for download files via CLI

1. Download files with curl to any directory.

   Example:
   
   ```bash
   curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/005/508/785/GCF_005508785.2_pea_aphid_22Mar2018_4r6ur_v2/GCF_005508785.2_pea_aphid_22Mar2018_4r6ur_v2_genomic.gff.gz
   curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/005/508/785/GCF_005508785.2_pea_aphid_22Mar2018_4r6ur_v2/GCF_005508785.2_pea_aphid_22Mar2018_4r6ur_v2_rna.fna.gz
   ```

3. Decompress files.
  
   ```bash
   gunzip *.gz
   ```

### Prepare BLAST database
Before running `hcrkit`, create a BLAST database from the reference RNA sequences:

```bash
makeblastdb -in transcripts.fasta -dbtype nucl -out database_name
```

Example:

```bash
makeblastdb -in GCF_005508785.2_pea_aphid_22Mar2018_4r6ur_v2_rna.fna -dbtype nucl -out Apis
```

### Running hcrkit

#### Syntax

```bash
hcrkit.py -i PATH -d PATH -p STR --initiator_id ID [OPTIONS]
```

#### Quick Examples
Basic usage (automatic on-target detection)
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
  --gene_name LOC10021649
```

With custom GC range and multiple threads
```bash
hcrkit.py \
  -i ApVas1.fasta \
  -d Apis \
  -p ApVas1 \
  --initiator_id S73 \
  --min_gc 40 \
  --max_gc 60 \
  -t 4
```

With pre-extracted target IDs
```bash
hcrkit.py \
  -i ApNos1.fasta \
  -d Apis \
  -p ApNos1 \
  --initiator_id A161 \
  --target_ids ApNos1_target_ids.txt
```

#### Parameters

**Required Parameters**

`-i, --input PATH`: Path to FASTA file of the target RNA sequence

- Example of FASTA file content:
  ```fasta
  >XM_003242475.4 PREDICTED: Acyrthosiphon pisum nanos-like protein (LOC100216490), transcript variant X3, mRNA
  AAAAAAATGTATGTTTTGAATTTGAAGTTTAAATTTAAGTGGGTGTTTTCAAGAGGACCCTTATTGCCAAGATATCACAA
  TGAGTGGCCGTAGCAGATTCGCGTTCGGCAGCCAGACCCAGGCGCCCCAGTACCAGAACCGATGCGCTCTGTGCTTGACC
  ...
  ```

`-d, --database PATH`: Path to BLAST database

- Specify the database name without extension (e.g., `Apis`)
- The database must be created beforehand using `makeblastdb` (see [Prepare BLAST database](#prepare-blast-database))

`-p, --prefix STR`: String used as a prefix for output files and a directory

`--initiator_id ID`: The name of fluorescently-labeled hairpin DNA

- The names which can be specified as options are only products lineup of ISHpalette™ Short hairpin amplifier (production of [Nepagene](https://nepagene.jp/en/products/fluorescent-tissue-staining-in-situ-hcr/ishpalette-2)): S23, S41, S45, S72, S73, A161
- The other names of products can be used, if you specify `--initiator_custom` parameter

**Optional Parameters**

`--target_ids PATH`: Path to file containing list of on-target transctipts IDs

- When omitted, `hcrkit` extracts the transcript ID from the FASTA header of `-i, --input PATH` file, taking the string up to the first whitespace

`--gff3 PATH`: Path to GFF3 file for automatic target ID extraction

- Requires `--gene_name` parameter
- Useful when the target gene has multiple isoforms or when transcript IDs are unknown
- Mutually exclusive with `--target_ids`

`--gene_name STR`: Gene name to search in GFF3 file

- Required when using `--gff3`
- Case-insensitive exact match

`--min_gc FLOAT`: Minimum GC percentage (default: 45.0)
`--max_gc FLOAT`: Maximum GC percentage (default: 55.0)

- The GC content can be changed with 40–60% (see [Guidline for probe design](#guidline-for-probe-design) section)

`--initiator_custom`: Path to CSV file of custom initiator sequences

- CSV format: initiator_name,sequence (no header)
- Example:
  
  ```csv
  B1,GAGGAGGGCAGCAAACGG
  B2,CCTCGTAAATCCTCATCA
  B3,GTCCCTGCCTCTATATCT
  ```

`--initiator_split INT`: Split position of initiator sequence (default: 9)

- The initiator sequence is split at this position for P1 and P2 probes
- Default value (9) is suitable for ISHpalette™ Short hairpin amplifiers

`-t, --threads`: Number of CPU threads for BLAST (default: 1)

`-h, --help`: Print help information

`-v, --version`: Print version information

**Output Files**

hcrkit creates the following output structure:

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

- `{prefix}`: string specified by `-p` parameter (e.g., if `-p ApNos1`, then `ApNos1_out/`)
- `{initiator}` = initiator ID (e.g., `A161`, `S73`)
- `{min}`, `{max}` = GC content range (e.g., `gc45-55`)

Example: Running with `-p ApNos1 --initiator_id A161` creates:
```
ApNos1_out/
├── ApNos1_A161_probe_pairs_gc45-55.csv
├── ApNos1_A161_probe_summary_gc45-55.txt
└── temp/
...
```

Main output files:</br>
`{prefix}_{initiator}_probe_pairs_gc{min}-{max}.csv`:

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

`{prefix}_{initiator}_probe_summary_gc{min}-{max}.txt`:

- Summary information for each probe set
- `Max off-target coverage (%)` indicates the maximum percentage of bases matching off-target transcripts
- Example:
  
  ```
  Probe set name    P1                                    P2                                       Max off-target coverage (%)
  ApNos1_A161_s38   GGTACGCGAaaAAGGGTCCTCTTGAAAACACCCACT  GGCCACTCATTGTGATATCTTGGCAaaAGGTAGGTGTAA  0
  ApNos1_A161_s363  GGTACGCGAaaAACTCGGCCAGTGTGTAGTAGTTGG  TGTTGCATGCGGTTCATGCGCAACAaaAGGTAGGTGTAA  0
  ...
  ```

Intermediate files (in `temp` directory):</br>

`{prefix}_region_candidates_gc{min}-{max}.fasta`:

- All probe region candidates before specificity filtering
- Example:

  ```
  >XM_003242475.4_region_38_89
  AGTGGGTGTTTTCAAGAGGACCCTTATTGCCAAGATATCACAATGAGTGGCC
  >XM_003242475.4_region_39_90
  GTGGGTGTTTTCAAGAGGACCCTTATTGCCAAGATATCACAATGAGTGGCCG
  ...
  ```

`{prefix}_blast_results_gc{min}-{max}.tsv`:

- Raw BLAST results in tabular output format 6
- Example:
  
  ```
  XM_003242475.4_region_38_89	XM_003242475.4	100.000	52	0	0	1	52	38	89	1.68e-22	103
  XM_003242475.4_region_38_89	XM_016802088.2	100.000	21	0	0	32	52	231	251	5.30e-04	42.1
  ...
  ```

`{prefix}_target_ids.txt` (when using `--gff3`):

- List of on-target transcript IDs extracted from GFF3
- Example:
  
  ```
  XM_016802089.2
  XM_016802088.2
  XM_003242475.4
  ```

`{prefix}_target_description.tsv` (when using `--gff3`):

- Detailed GFF3 entries for extracted target transcripts

## Tips

The number of probes can affect RNA detection. Starting with a 10-probe set is recommended. A larger number of probes may enable detection of low-copy RNAs.</br>

**If the number of probes were few...**

- Relax GC content requirements from default (45-55%) to 40-60%:
  
  ```bash
  hcrkit -i target.fasta -d database -p output --initiator_id A161 --min_gc 40 --max_gc 60
  ```

- Review the on-target list:

  - When using `--gff3`: Check the console output for "Found X mRNA entries" and similar genes
  - When using `--target_ids`: Verify that the file contains all relevant transcript IDs
  - Missing on-target transcripts may cause incorrect filtering

## Overview

### What is isHCR (_in situ_ Hybridization Chain Reaction)?

isHCR is a powerful technique for RNA detection *in situ*. Unlike enzyme-based methods, isHCR visualize RNAs though formation of polymors composed of fluorescently-labeled oligonucleotides.</br>

![isHCR](images/isHCR.png)
<p align="center"> Figure 1 </p>

The probe contains an initiator sequence (Figure 1A). The initiator hybridizes with hairpin DNAs, triggering a hybridization chain reaction that produces fluorescently labeled polymers (Figure 1B). To suppress background signals, the probe is split (P1 and P2), and the pair of probes can efficiently trigger amplification.

### Guidline for design of probes

-	GC content: 40–60% (45–55% recomended)
-	For each probe pair, sequence identity with off-target RNAs: ≤ 50% (verify with BLAST)
-	Design region in target RNA: CDS + UTR (CDS recomended)

### Workflow & algorithm

`hcrkit` automates the probe design process through th following workflow:

#### Step 1. Load sequence
`hcrkit` loads the target RNA sequence from the input FASTA file (`-i` parameter).

#### Step 2. Load Target IDs
`hcrkit` determines which transcripts should be considered as on-target for specificity filtering:

- **Automatic mode (default)**: Extracts transcript ID from the FASTA header
  - Example: `>XM_003242475.4 description...` → on-target ID: `XM_003242475.4`

- **GFF3 mode (`--gff3` + `--gene_name`)**: Extracts all transcript IDs for the specified gene
  - Searches the GFF3 file for mRNA entries where `gene=` matches the specified gene name (case-insensitive)
  - Example: `--gene_name LOC10021649` finds all nanos isoforms
  - Useful when:
    - The target gene has multiple isoforms
    - You don't know the exact transcript ID
    - You want to avoid off-target hits to other isoforms of the same gene

- **Manual mode (`--target_ids`)**: Uses IDs from the provided file
  - Useful for your own data not registered in NCBI
  - Gives precise control over on-target definition

#### Step 3. Generate Probe Region Candidates
![load](images/step1.png)

**Generate 52-nt sliding windows**

The probe region for a split probe pair is 52 nt. `hcrkit` divides the target transcript into overlapping 52-nt fragments using a sliding window (moving 1 nt at a time from 5′ to 3′ end).

**GC content filtering**

Each 52-nt probe region is split into:
- P1 binding site: first 25 nt (positions 1-25)
- Spacer: 2 nt (positions 26-27)
- P2 binding site: last 25 nt (positions 28-52)

`hcrkit` calculates the GC content for both P1 and P2 binding sites. Only regions where **both** binding sites meet the GC criteria (default: 45-55%) are retained.</br>

**Remove duplicates**

If multiple probe regions have identical sequences, `hcrkit` keeps only one and removes the duplicates.

#### Step 4. BLAST Search
![step2](images/step2.png)

`hcrkit` performs BLASTN search against the reference transcriptome database to identify potential off-target binding sites for each probe region candidate.

#### Step 5. Filter by Specificity

For each BLAST hit, `hcrkit`:
1. Checks if the hit is to an on-target transcript (based on the IDs determined in step 2)
2. For off-target hits, calculates coverage: `(alignment length) / 52 * 100`
3. Tracks the maximum off-target coverage for each probe region

`hcrkit` removes probe regions where the maximum off-target coverage is ≥50%. This ensures that each probe has sufficient specificity to the target transcript.

#### 6. Select Non-overlapping Probe Regions
![step3](images/step3.png)

After specificity filtering, many valid probe regions may overlap with each other. To maximize probe coverage across the transcript, `hcrkit` selects non-overlapping regions using a greedy algorithm:

1. Sort all valid probe regions by start position (5′ to 3′)
2. Select the first region
3. Skip any regions that overlap with the selected region
4. Select the next non-overlapping region
5. Repeat until all regions are processed

This approach ensures that selected probe regions do not overlap, maximizing the number of probes that can be used simultaneously.

#### 7. Generate Probe Pairs
For each selected probe region, `hcrkit` generates a split probe pair (P1 and P2):

**Reverse complementation**

Since probe regions are extracted from the sense strand, `hcrkit` converts them to antisense sequences (which will bind to the target mRNA) by reverse complementation.

**Add initiator sequences**

The initiator sequence is split at the specified position (default: 9 nt):
- P1 probe: `[initiator (1-9)] + aa + [reverse_complement(P1 region)]`
- P2 probe: `[reverse_complement(P2 region)] + aa + [initiator (10-21)]`

The "aa" spacer sequences are added at the junction between initiator and binding region.

#### 8. Write Outputs

The final probe sequences are formatted as:
- CSV file: for ordering oligonucleotides
- Summary file: with probe set names and off-target coverage information

This workflow is implemented across `hcrkit.py` (main pipeline) and `core.py` (core functions), working together to automate the entire probe design process.

## Citation

Preparing...