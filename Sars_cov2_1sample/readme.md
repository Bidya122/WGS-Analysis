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

I have given the process of the nf file below but you can also find the nf file attached in this folder.

```bash
process TRIMMOMATIC {
publishDir "${params.outdir}/trimmed", mode: 'copy'
cpus = 1
    memory = '2 GB'
    input:
    path sample_file

    output:
    path "trimmed_${sample_file.name}"

    script:
    """
    echo "Running Trimmomatic on $sample_file"
    java -jar /home/bidya122/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 2 \
        $sample_file trimmed_${sample_file.name} SLIDINGWINDOW:4:20 MINLEN:50
    """
}

process FASTQC {
publishDir "${params.outdir}/fastqc", mode: 'copy'
cpus = 1
    memory = '1 GB'

    input:
    path trimmed_file

    output:
    path "*_fastqc.html"
    path "*_fastqc.zip"

    script:
    """
    echo "Running FastQC on ${trimmed_file}"
    fastqc ${trimmed_file} -o .
    """
}
```

## Trimming and Quality Control 
<img width="608" height="332" alt="image" src="https://github.com/user-attachments/assets/b8c6bc83-c0a4-4f0f-87ca-0be941edc179" />  
<img width="940" height="574" alt="image" src="https://github.com/user-attachments/assets/c86282fa-2f37-4fa3-a47a-dbabf463832c" />  
<img width="940" height="708" alt="image" src="https://github.com/user-attachments/assets/e1632d68-9a40-41cb-8653-27ce0298b2cc" />  
<img width="940" height="559" alt="image" src="https://github.com/user-attachments/assets/62f99647-e7b8-4fcf-b8f7-e1a13000d3ee" />  
<img width="940" height="580" alt="image" src="https://github.com/user-attachments/assets/7f1f8435-2eec-4092-8b20-32311f3fb9e2" />  

The sample on its quality check passed Basic Statistics, Per Base Sequence Quality, Per sequence Quality Scores, Per Base N Content and Adapter Content, which means data file looks normal - correct read length, GC %, etc. No sequences were flagged as poor quality, Read lengths fall within the expected range (50–301 bp) after trimming, GC content (39%) is within a biologically reasonable range, Encoding format is consistent with Illumina sequencing (Sanger / Illumina 1.9), A sufficient number of reads (3,059,010) were retained for downstream analysis. The blue line shows the average quality at each base position across all reads in the per base sequence quality graph.   
The Per Sequence Quality Scores shows that No reads have an average quality below 30. In other words: Every single read in your dataset has average Phred ≥ 30. Which means every read is 99.9% accurate or better. N content usually means "an unknown base" or the sequencer puts an N, when machine can’t decide what base is there — for example, the signal was too weak or mixed and the last one shows there are no adapter remnannts. 

<img width="769" height="575" alt="image" src="https://github.com/user-attachments/assets/d3ea9d80-a648-48a1-b964-1a0e1845f9f4" />  
<img width="722" height="366" alt="image" src="https://github.com/user-attachments/assets/941a039f-8538-474a-8ead-16c118fe8d38" />  

Although FastQC flagged the sequence length distribution with a warning, this reflects expected post-trimming variability rather than data quality issues. No corrective action was required, and downstream processing was continued. Overrepresented sequences are sequences (usually short reads or parts of reads) that occur much more frequently than expected in the dataset. FastQC says “No Hit” for the overrepresented sequences, that means:  
•	Those repeated sequences don’t match any known adapter, primer, or contaminant in its built-in database, and
•	They’re just repeated short reads that occur more often than expected.

<img width="761" height="372" alt="image" src="https://github.com/user-attachments/assets/3244748b-87ab-423a-bf40-e7532b712593" />  
<img width="841" height="646" alt="image" src="https://github.com/user-attachments/assets/feda942d-0eda-4b7e-b012-3415969349c3" />    
<img width="940" height="686" alt="image" src="https://github.com/user-attachments/assets/c2427c62-9ccd-4d90-b560-864bbdd07490" />  

In Per base sequence content That pattern = sequence data is AT-rich — meaning the genome or transcript that was sequenced naturally contains more adenines (A) and thymines (T) than guanines (G) and cytosines (C).
As it’s SARS-CoV-2, which is indeed AT-biased (~62% A+T). So that exact separation (A/T lines above, G/C lines below) is expected and biologically correct for this virus.
In per sequnce GC content, for random genomic DNA, the plot should look like:
•	A smooth, single bell-shaped curve (like a hill),
•	Centered around the organism’s expected GC content.  
o	e.g. Human ~40–50%  
o	E. coli ~51%  
o	SARS-CoV-2 ~38%  
So the peak of the curve should roughly match the expected GC % for your species or genome and it does!!  
Although FastQC flagged sequence duplication levels as a failure, this is expected for SARS-CoV-2 sequencing due to the small viral genome and high sequencing depth. The observed duplication reflects biological and experimental enrichment rather than poor data quality. Therefore, this QC result did not prevent downstream analysis.  

Following successful quality control, reads were aligned to the indexed reference genome using BWA. Before sequencing reads can be aligned to a reference genome, the reference must be indexed. Indexing does not cut or modify the genome. Instead, it preprocesses the reference genome into a searchable data structure that allows BWA to locate short DNA sequences quickly and efficiently.   
You can think of indexing like preparing a table of contents for a very large book. Rather than scanning the entire book every time, BWA uses this index to jump directly to the most likely locations where a read may align.    
Sequencing reads are very short compared to the full genome.
Indexing is required because, Without indexing, BWA would need to scan the entire genome for every read, which would be extremely slow. By indexing the reference genome:  
- Read alignment becomes fast
- Memory usage is efficient
- Millions of reads can be aligned in a practical time
I have given the process below but, you can also find it on the main1.nf file that I have provided

```bash
process BWA_INDEX {
    publishDir "${params.outdir}/bwa_index", mode: 'copy'

    input:
    path ref_file

    output:
    path "bwa_index", emit: indexed_ref_dir

    script:
    """
    mkdir bwa_index
    cp ${ref_file} bwa_index/
    cd bwa_index
    bwa index `basename ${ref_file}`
    """
}
```
<img width="940" height="147" alt="image" src="https://github.com/user-attachments/assets/e681398f-76d5-4288-bca0-a4409237b9ed" />    

This step aligns sequencing reads to a pre-indexed reference genome using BWA-MEM and generates sorted and indexed BAM files. The alignment output is converted from SAM to BAM format, sorted by genomic coordinates, and indexed to enable efficient downstream analysis.  
These processed BAM files are essential for tasks such as variant calling, coverage analysis, and visualization. The sorted and indexed BAM files can be directly loaded into IGV (Integrative Genomics Viewer) for manual inspection of read alignments and mapping quality. The resulting alignment files are suitable for versioning and sharing as reproducible outputs of the pipeline.  

```bash
process ALIGN_READS {
    publishDir "${params.outdir}/alignment", mode: 'copy'

    input:
    tuple path(reads), path(indexed_ref_dir)

    output:
    path "*.sorted.bam"

    script:
    """
    echo "Aligning reads: ${reads} using reference inside ${indexed_ref_dir}"

    # Get the FASTA file from inside the indexed reference directory
    ref_fasta=\$(ls ${indexed_ref_dir}/*.fasta)

    echo "Using reference: \${ref_fasta}"

    bwa mem \${ref_fasta} ${reads} > \$(basename ${reads} .fastq.gz).sam
    samtools view -S -b \$(basename ${reads} .fastq.gz).sam > \$(basename ${reads} .fastq.gz).bam
    samtools sort \$(basename ${reads} .fastq.gz).bam -o \$(basename ${reads} .fastq.gz).sorted.bam
    samtools index \$(basename ${reads} .fastq.gz).sorted.bam
    """
}
```
<img width="940" height="260" alt="image" src="https://github.com/user-attachments/assets/1e4d0819-b7b7-4c0b-8ddc-2f285f4da751" />  





  
















