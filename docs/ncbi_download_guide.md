# How to Download GFF and Reference RNA Sequences Files 

## I. Note
This document was created on October 3, 2025. <br>
The NCBI (National Center of Biotechnology Information) interface may change over time. 

## II. Scope of this page
- **GFF files** are required to identify transcript IDs of your target genes
- **Reference transcript sequences** are used to create BLAST detabases for specificity screening
- Our programs are designed for NCBI-format flles (may not compatible with species-specific databases like FlyBase)

## III. Download
1. **Access NCBI homepage** at [https://www.ncbi.nlm.nih.gov](https://www.ncbi.nlm.nih.gov) and enter your species name in the search bar. 
![Figure 1](images/fig1.png)
<br>

2. **Click "Assembly/Genomes"**. 
![Figure 2](images/fig2.png)
<br>

3. **Select an appropriate assembly** for your experiments. In most cases, choose the assembly marked as the reference genome (indicated by a green checkmark).
![Figure 3](images/fig3.png)
<br>

4. **Click "FTP"**. 
![Figure 4](images/fig4.png)
<br>   

5. **Download the required files**.
- `GCF_*_genomic.gff.gz`: gene annotations
- `GCF_*_rna.fna.gz`: transcript sequences
![Figure 5](images/fig5.png)
<br>

> [!IMPORTANT]
These files can exceed 10 GB. Ensure sufficient storage capacity before downloading.


