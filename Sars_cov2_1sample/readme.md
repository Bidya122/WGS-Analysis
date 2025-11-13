# WGS_Analysis
This repository contains scripts, workflows, and documentation for performing Whole Genome Sequencing (WGS) analysis . All steps were implemented hands-on, including buildup of end-to-end workflow curation on nextflow, raw data processing, alignment, QC, and downstream analysis, demonstrating practical proficiency in WGS-Seq workflows.


# Introduction  
For this analysis, I chose the SARS-CoV-2 genome. The sequencing data was obtained from a human sample available in the SRA database (NCBI). The entire workflow was implemented on WSL1, covering end-to-end WGS analysis including quality control, alignment, variant calling, and annotation.  
The main goal of this project was to learn and demonstrate a complete WGS pipeline, preparing me to scale this workflow to larger datasets or different organisms in future analyses.  
- **Sample**:  SRX9620660  [Find Dataset](https://www.ncbi.nlm.nih.gov/sra/?term=SRX9620660)  
<img width="723" height="598" alt="image" src="https://github.com/user-attachments/assets/256f0d1b-efcf-4e62-8284-d605f791d9bc" />

```bash
#to download the file directly
wget -r -nH --cut-dirs=5 -A "*.fastq.gz" ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR131/025/SRR13182925/
```
(OR)
```bash
# Install SRA Toolkit 
sudo apt install sra-toolkit

# Download and convert to FASTQ
prefetch SRR13182925
fasterq-dump SRR13182925

```


- **Reference**: NCBI > Assembly database> search for GCF_009858895.2 > click FTP to download the files
<img width="822" height="635" alt="image" src="https://github.com/user-attachments/assets/fb784126-e9dd-4e71-b4aa-46fc5075dfd2" />  

```bash
#to download the file directly
curl –O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna.gz
```

```bash
#to download the file directly
gunzip GCF_009858895.2_ASM985889v3_genomic.fna.gz
```
```bash
#to rename the file
mv GCF_009858895.2_ASM985889v3_genomic.fna NC_045512.2.fasta
```
