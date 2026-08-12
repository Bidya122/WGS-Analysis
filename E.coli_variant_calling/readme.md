# E. coli Phage Coevolution — WGS Variant Analysis  
## Introduction  
Bacteriophages, or phages, are viruses that infect bacteria. When Escherichia coli populations are exposed to phages, susceptible bacterial cells can be infected and eliminated, creating selective pressure for bacteria that acquire or already possess mutations that provide a survival advantage.  
However, this interaction is not one-sided. As bacteria evolve resistance to phages, phages can also evolve to overcome bacterial resistance. This continuous process of adaptation between bacteria and phages is known as phage–bacteria coevolution and can be viewed as an evolutionary arms race.  
For this project, I am using a publicly available whole-genome sequencing (WGS) dataset of E. coli K-12 obtained from the NCBI Sequence Read Archive (SRA). The dataset comes from a study investigating bacterial–phage interactions under different phage exposure and coevolution conditions.  

## Dataset Details  
Experiment: SRX32505351
Run: SRR37624945
Organism: Escherichia coli K-12
Sample: T2-First Late Day 6 R5.3
Sequencing: Illumina NextSeq 2000
Strategy: Whole-genome sequencing (WGS)
Layout: Single-end

The purpose of using this dataset is not to reproduce the original study, but to use it as a practice dataset for learning and documenting a complete WGS variant-analysis workflow.

The basic biological question underlying the analysis is:  
What genetic variants can be identified in an evolved E. coli population following phage-associated selection, and what might these variants tell us about bacterial adaptation? The analysis will therefore take us from raw sequencing reads → genomic variants → biological interpretation.    

## Dataset Acquisition
The dataset was obtained from the NCBI Sequence Read Archive (SRA), which provides publicly available sequencing data. The dataset was identified by searching the NCBI SRA database using the experiment accession: SRX32505351    
The SRA record provided information about the sequencing experiment, including the organism, sample, sequencing platform, library strategy, and associated sequencing run.    
For this project, the relevant records are:    
Study       → SRP683594  
Sample      → SRS28398370  
Experiment  → SRX32505351  
Run         → SRR37624945  

[Dataset Overview](https://www.ncbi.nlm.nih.gov/sra/?term=SRX32505351) 

<img width="1107" height="602" alt="image" src="https://github.com/user-attachments/assets/a2d4d7d4-c597-4ca9-9bdb-04c9bfb5279b" />    

Then I Clicked on the SRR to download the FASTQ file for my analysis. 

<img width="1462" height="733" alt="image" src="https://github.com/user-attachments/assets/993396d7-ee71-4f1a-b27b-8a768ea72ef9" />    

The raw sequencing data for run SRR37624945 were obtained from the NCBI SRA in FASTQ format. Since the analysis was performed using WSL (Ubuntu), the downloaded FASTQ file was first unzipped and then moved into the Ubuntu home directory so that it could be accessed easily during the downstream analysis. The resulting raw sequencing file was `SRR37624945.fastq`.

But if you directly want to uncompress using the terminal then below is the command
```
gunzip SRR37624945.fastq.gz
```
## Installation of Tools
**SRA Toolkit Installation - necessary to access, download, and convert sequencing data from the NCBI Sequence Read Archive (SRA) into usable formats for analysis.**  
```bash
sudo apt update && sudo apt upgrade -y
```
```bash
sudo apt install sra-toolkit -y
```
```bash
which prefetch
```
```bash
which fastq-dump
```
**Fastq Installation.**  
```bash
sudo apt install fastqc -y
```
```bash
fastqc --version
```
**Trimmomatic Installation**
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
You should see something like this 

<img width="707" height="48" alt="image" src="https://github.com/user-attachments/assets/23194517-1525-43e9-b683-16d0a4251759" />

## Data Visualization and Quality Control

Before performing quality control, I first inspected the raw FASTQ file to understand its basic structure. A FASTQ file stores sequencing reads along with the quality information associated with each base. Each read is represented by four lines:

`@Read_ID`  
`SEQUENCE`  
`+`  
`QUALITY_SCORE`  

For example:
@SRR37624945.1  
ACGTACGTACGT...  
+  
IIIIIIIIIIII...  
What the four lines represent  
- Read identifier (@): Contains information used to identify the sequencing read.  
- Nucleotide sequence: Contains the DNA bases determined by the sequencer (A, T, G, C, and sometimes N for an undetermined base).  
- Separator (+): Marks the beginning of the quality-score line. It may optionally repeat the read identifier.  
- Quality scores: Contains one character for each base in the sequence. These characters encode the Phred quality score, which represents the confidence of the sequencer's base call.  
 
```bash
#To check the raw fastq file
head SRR37624945.fastq
```
```bash
#To check the first twenty lines of the file
head -n 20 SRR37624945.fastq
```
<img width="940" height="295" alt="image" src="https://github.com/user-attachments/assets/80052e8d-b5ca-4b57-9ebb-26af6295f225" />

FastQC was performed on the downloaded FASTQ file to assess the quality of the sequencing reads. It provides information on read quality scores, GC content, sequence duplication, and other quality-related parameters to identify any potential issues before further analysis.

```bash
#To check the quality of the fastq file
fastqc SRR37624945.fastq
```
[Fastqc Output](https://github.com/Bidya122/WGS-Analysis/blob/main/E.coli_variant_calling/SRR37624945_fastqc.html) 

**FastQC Report**

<img width="1256" height="612" alt="image" src="https://github.com/user-attachments/assets/2022f271-de7b-45c2-bc67-63098133f8b8" />

FastQC was performed to assess the quality and characteristics of the sequencing reads before downstream analysis. The dataset contained 110,951 sequences with approximately 231.3 Mbp of total sequence data and a GC content of 50%. No sequences were flagged as poor quality. Some FastQC modules showed warnings or failures, which were further examined to identify potential quality or sequence-composition issues.

<img width="1120" height="783" alt="image" src="https://github.com/user-attachments/assets/6a69fadd-5b56-4f8c-b830-b4071b9d141f" />

The reads showed relatively poor and variable quality at the beginning, as indicated by the lower blue line and wider yellow boxes. The quality gradually improved with increasing read length and became good and relatively consistent later in the reads. The blue line represents the mean (average) quality score, while the red line represents the median quality score within each position. When these values fall into the lower-quality region, it indicates poorer-quality bases. Hence, it needs trimming for the initial few bases. 

<img width="1005" height="753" alt="image" src="https://github.com/user-attachments/assets/66ab2d2d-348d-4689-8f24-99940c7adfc2" />

The Y-axis represents the number of reads, while the X-axis represents the mean quality score of the reads. The majority of reads are concentrated around Q38–Q40, indicating that most reads have a high average sequencing quality.

<img width="1000" height="731" alt="image" src="https://github.com/user-attachments/assets/ba4747c8-7d65-47d9-b3f8-2373d1fca951" />

The sequence content shows substantial variation in nucleotide composition at the beginning of the reads, resulting in a FastQC failure. At the first few positions (1–9 bp), the nucleotide percentages fluctuate considerably. For example, C reaches ~47% at position 2, while G reaches ~46% at position 4. After the initial region, the four nucleotide lines become very close and stable around ~25% each. So, after the initial bases, the A, T, G and C contents become stable and remain close to 25% each, indicating a balanced and consistent nucleotide composition across the majority of the reads.












