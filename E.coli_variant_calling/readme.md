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

**BWA Installation**

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

<img width="1002" height="747" alt="image" src="https://github.com/user-attachments/assets/aa478848-502c-4f93-988f-88c12e14a471" />

The per-sequence GC content shows that most reads have a GC content of approximately 52–53%. Although the observed distribution is broadly close to the theoretical distribution, some deviation is evident, particularly in the shape and width of the distribution. Therefore, FastQC has flagged this parameter as failed. However, the overall GC content is reasonably consistent with that expected for E. coli, and this result alone does not necessarily indicate a major quality problem or contamination.    

<img width="1022" height="727" alt="image" src="https://github.com/user-attachments/assets/8b12dcee-bdb5-4609-bce9-ca522b820a12" />

No unknown bases detected.

<img width="990" height="722" alt="image" src="https://github.com/user-attachments/assets/f1b3a08c-6caa-49b0-9d27-8b2aca1fbc58" />

The sequence length distribution shows that most reads are concentrated within the shorter sequence-length range, with a sharp decrease in read counts as sequence length increases. A small number of longer reads are also present, indicating variability in read lengths. FastQC has therefore assigned a warning for this parameter.   

The remaining FastQC parameters, including sequence duplication levels, overrepresented sequences, and adapter content, showed no significant abnormalities. The absence of excessive sequence duplication suggests good library complexity, while no overrepresented sequences or adapter contamination indicates that the reads are not substantially affected by unwanted sequence enrichment or adapter contamination. Overall, these parameters indicate good sequencing quality.    

**Trimming**

Following the initial quality assessment using FastQC, the sequencing reads were subjected to quality trimming. Trimming was performed to remove low-quality bases and any residual adapter sequences, thereby improving the overall quality of the reads for downstream analysis. The trimmed reads were subsequently used for further analysis. 

```bash
java -jar /home/bidya/tools/Trimmomatic-0.39/trimmomatic-0.39.jar SE -phread33 \
SRR37624945.fastq \
SRR37624945_trimmed.fastq \
LEADING:20 SLIDINGWINDOW:4:20 MINLEN:36
```
Trimmomatic was used to remove low-quality bases and improve the overall quality of the sequencing reads. The LEADING:20 parameter was used to remove low-quality bases from the beginning of reads, while SLIDINGWINDOW:4:20 removed regions where the average quality score within a 4-base window fell below 20. Reads shorter than 36 bp after trimming were discarded using MINLEN:36, as very short reads may not be suitable for reliable downstream analysis. The Phred+33 quality score encoding was specified using the -phred33 option.

<img width="1463" height="163" alt="image" src="https://github.com/user-attachments/assets/c5988817-e01d-4f09-b20b-437e59380237" />

Following trimming, FastQC was performed again on the trimmed FASTQ file to assess the effectiveness of the quality-control step. The post-trimming FastQC results were compared with the initial quality assessment to confirm the removal of low-quality bases and to check for any remaining quality issues before proceeding with downstream analysis.

```bash
fastqc SRR37624945_trimmed.fastq
```
<img width="851" height="528" alt="image" src="https://github.com/user-attachments/assets/3b0d24da-e415-4b57-a952-2b8ab69361c6" /> 

<img width="1320" height="746" alt="image" src="https://github.com/user-attachments/assets/bf2e9cf9-f7ac-4fdd-8d37-40dd33c838bf" />

<img width="1322" height="743" alt="image" src="https://github.com/user-attachments/assets/107afa0b-f7d2-498e-9317-10e6c527d903" />

<img width="1070" height="776" alt="image" src="https://github.com/user-attachments/assets/f569f5c7-a744-4c89-bff2-ed2ef8f3ecdb" />

Additional trimming was not performed because the post-trimming FastQC report showed a substantial improvement in overall sequence quality. The Per Base Sequence Quality module showed consistently high Phred scores, predominantly above Q40, across most of the read length, indicating excellent base-level sequencing quality. However, the Per Base Sequence Content module remained flagged due to variations in the proportions of A, T, G, and C at different positions. Since this dataset represents whole-genome sequencing of E. coli, variations in nucleotide composition may reflect the underlying genomic sequence composition as well as technical factors introduced during library preparation. Further quality trimming was therefore not considered necessary, as it was unlikely to correct the observed sequence-content variation and could result in unnecessary loss of useful genomic sequence. After trimming, the sequence length distribution showed a predominance of shorter reads, with most sequences concentrated below approximately 300 bp. The presence of a smaller number of longer reads indicates variation in read length following trimming. FastQC assigned a warning to this parameter, which is consistent with the variable read lengths produced during quality trimming.    


## Reference Genome Download    
The dataset was identified as Escherichia coli K-12. Based on the strain information, E. coli K-12 MG1655 was selected as the reference strain. The corresponding reference genome assembly, ASM584v2, was obtained from the NCBI genome database in FASTA format. This reference genome was subsequently used for alignment of the trimmed sequencing reads.  

Go to NCBI > Search "E coli K12 NCBI Genome" > Select of the correct strain and substrain and download the reference genome.    
[Click to Download Reference Genome](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000005845.2/)    





























