# WGS_Analysis
This repository contains scripts, workflows, and documentation for performing Whole Genome Sequencing (WGS) analysis . All steps were implemented hands-on, including buildup of end-to-end workflow curation on **nextflow**, raw data processing, alignment, QC, and downstream analysis, demonstrating practical proficiency in WGS-Seq workflows. The goal of this project was to get fimiliar with working with nextflow and understand variant calling and analysis.


# Introduction  
**Project Overview**  
This project is a Nextflow DSL2 workflow for processing SARS-CoV-2 sequencing data. It performs:  
-Quality trimming with Trimmomatic  
-Quality control with FastQC  
-Alignment to the SARS-CoV-2 reference genome using BWA  
-Sorting and indexing BAMs  
-Merging BAMs, counting total reads, and computing coverage  

**Learning Context**  
During this exercise, I ran Trimmomatic in single-end (SE) mode on a dataset that was originally paired-end.
I made this choice by mistake initially, but I continued the workflow to learn Nextflow, FASTQC, BWA alignment, and BAM handling. The purpose of this project is educational, not to produce publication-grade results.
The workflow is designed as a learning exercise to practice building reproducible bioinformatics pipelines.  
The main goal of this project was to learn and demonstrate a complete WGS pipeline, preparing me to scale this workflow to larger datasets or different organisms in future analyses. This workflow can be easily adapted to true paired-end processing by modifying the Trimmomatic and alignment steps to handle read pairs.   

**Key Notes About the Workflow**

Trimming -   
Each FASTQ file was trimmed independently using SE Trimmomatic.  
For this small viral genome, SE trimming is sufficient to achieve good coverage.  

Alignment -  
Trimmed reads were aligned to the SARS-CoV-2 reference genome (NC_045512.2).  
Although the original data were paired-end, alignment as SE reads is reasonable because the genome is small (~30 kb) and amplicon-based.  

MERGE_BAMS -   
Sorted BAMs for each sample are merged into a single BAM.  
The merged BAM is indexed, total reads counted, and per-base coverage computed.  
Treating paired-end reads as SE does not significantly affect coverage or variant detection for a viral genome of this size.  
For larger genomes or structural variant detection, proper PE alignment would be necessary.  

# Sample
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

# Reference
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

# CONDA_Installation
```bash
#downlaod miniconda installer
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```
```bash
#Run the installer
bash Miniconda3-latest-Linux-x86_64.sh
```
```bash
#list the installation directory
ls ~/miniconda3
```
```bash
#add conda to path
echo 'export PATH=$HOME/miniconda3/bin:$PATH' >> ~/.bashrc
```
```bash
#reload bash configuration
source ~/.bashrc
```
```bash
conda –version
```
# TOOLS INSTALLATION
## FASTQC
