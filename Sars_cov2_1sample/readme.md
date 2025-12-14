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
conda --version
```
# TOOLS INSTALLATION
## FASTQC   
```bash
sudo apt install fastqc -y
```
```bash
fastqc --version
```
## TRIMMOMATIC
(on home directory)
```bash
mkdir -p ~/tools
```
```bash
cd ~/tools
```
```bash
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
```
```bash
unzip Trimmomatic-0.39.zip
```
```bash
ls
```
```bash
cd Trimmomatic-0.39
```
## BWA INDEX
```bash
#Using APT (Ubuntu/WSL)
sudo apt update
sudo apt install bwa
```
```bash
#Using Conda (Bioconda)
conda create -n bioinfo -y
conda activate bioinfo
conda install -c bioconda bwa
```
## BCFTOOLS
```bash
#Creates a new Conda environment named variant_calling and installs bcftools → for manipulating VCF/BCF files, bedtools → for genomic interval operations, snpeff → for variant annotation
conda create -n variant_calling -c bioconda -c conda-forge bcftools bedtools snpeff
```
```bash
#Tabix is required to index compressed VCF files (.vcf.gz). It also installs bgzip, which is needed for compressing VCFs.
sudo apt install tabix
```
```bash
bcftools --version
```
```bash
bedtools --version
```
```bash
tabix --version
```
## SNPEFF
```bash
conda install -c bioconda snpeff
```
```bash
snpEff --version
```
```bash
conda list | grep -i snp 
```
```bash
find ~/miniconda3/ -type f -iname "snpEff*.jar"
```

# NEXTFLOW.CONFIG 
nextflow.config is the core configuration file that tells Nextflow how to run the pipeline. I have uploaded the config file [Find config file](https://github.com/Bidya122/Bulk_RNAseq_Analysis/tree/main/Data)

<img width="720" height="170" alt="image" src="https://github.com/user-attachments/assets/0f8f91bf-6c28-4d09-ba4e-447bf174dbae" />  

Nextflow will run jobs locally (same machine). Each process will use 2 CPUs and 1 GB memory. It will use the Conda environment described in environment.yml. 

<img width="341" height="72" alt="image" src="https://github.com/user-attachments/assets/0926f94a-5058-4153-8662-50b3691602cf" />  

All pipeline output will go into a folder called results. It can also be override while running:
```bash
nextflow run main.nf -profile conda --outdir output_folder
```
<img width="406" height="82" alt="image" src="https://github.com/user-attachments/assets/7de0d25e-fda1-4494-98ef-f47095c919f7" />  

This is optional
enabled = true → Nextflow knows to activate conda.
useMamba = true → Faster and smarter installation solver. 

<img width="455" height="132" alt="image" src="https://github.com/user-attachments/assets/02df77b1-0049-4b10-a6b5-b9ac02924acd" />  

withName: ".*" means apply this to all processes (because .* matches everything).  
Sets:  
1 CPU per process  
1 GB memory 
This overrides earlier settings unless a specific process has its own resources defined. Why have two process blocks?  
The first process block gives global defaults.  
The second one refines or overrides them for all processes matched by the pattern.  

<img width="439" height="219" alt="image" src="https://github.com/user-attachments/assets/a59c3221-425a-44d8-ab0b-b6667e618df6" />  

Runs pipeline locally without using Conda environments.  
Gives 2 CPUs and 1 GB to each process.  
Saves output to results.  

# Main1.nf NEXTFLOW WORKFLOW
In the earlier sections, I have described the configuration file (nextflow.config), including profiles and execution settings. From this point onward, the focus is on the Nextflow workflow file (.nf).Since this was my first time working with Nextflow, I chose to build the pipeline incrementally, starting with one tool and one process at a time rather than designing a full workflow at once. This approach helped me:  
- Understand the structure of a Nextflow pipeline  
- Learn how process, workflow, and channels interact  
- Debug errors more effectively    
- Ensure each step works independently before scaling  
The file main1.nf therefore represents a minimal, functional workflow, intended as a learning and validation step before extending the pipeline further.
Starting with Trimming and QC check of the samples.

<img width="940" height="166" alt="image" src="https://github.com/user-attachments/assets/ff3f10e9-381a-4e08-9e90-3f170f0348d3" />  

## Quality Control 






