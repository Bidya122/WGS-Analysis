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

**SAMTOOLS Installation**
```bash
sudo apt install samtools
samtools --version
```

**BCFTOOLS Installation**
```bash
sudo apt install bcftools
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

Go to NCBI > Search "E coli K12 NCBI Genome" > Select of the correct strain and substrain and download the reference genome > Click FTP > Click on "GCF_000005845.2_ASM584v2_genomic.fna.gz" > Decompress it Take this file to the working directory on Ubuntu  

[Click to Download Reference Genome](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000005845.2/)    

## Indexing the Reference Genome
Before performing read alignment, the downloaded E. coli K-12 substr. MG1655 reference genome (ASM584v2) was indexed using BWA. Indexing generates auxiliary files that allow BWA to efficiently search and align sequencing reads against the reference genome. The reference genome was indexed using the following command:

```bash
bwa index  GCF_000005845.2_ASM584v2_genomic.fna
```
<img width="1020" height="232" alt="image" src="https://github.com/user-attachments/assets/8fb3d98c-10d7-436d-8bc9-3400d8235b86" />

<img width="702" height="227" alt="image" src="https://github.com/user-attachments/assets/a72ae888-08ac-4d49-a550-359546397359" />

Indexing the reference genome using BWA generated several auxiliary index files with the extensions .amb, .ann, .bwt, .pac, and .sa. These files contain different indexing information required by BWA to efficiently search the reference genome during read alignment. The original FASTA reference genome was retained, while the generated index files were used by BWA during the subsequent alignment step.

## Alignment of trimmed fastq to the Reference Genome
Following reference genome indexing, the trimmed sequencing reads were aligned to the Escherichia coli K-12 substr. MG1655 reference genome (ASM584v2) using BWA-MEM. The purpose of alignment was to determine the genomic position of each sequencing read relative to the reference genome. This allows the reads to be mapped to their corresponding regions of the E. coli genome and provides the basis for subsequent analyses such as genome coverage and variant identification. Since the dataset consists of single-end reads, the trimmed FASTQ file was aligned against the indexed reference genome using the bwa mem command. The alignment output was initially generated in SAM (Sequence Alignment/Map) format, which contains information about how each read aligns to the reference genome.  

```bash
bwa mem GCF_000005845.2_ASM584v2_genomic.fna SRR37624945_trimmed.fastq  > aligned.sam
```
<img width="1456" height="210" alt="image" src="https://github.com/user-attachments/assets/0ad5d4f8-7511-41aa-9324-5aff2d029d56" />

<img width="1455" height="537" alt="image" src="https://github.com/user-attachments/assets/d209c4fc-afc7-4fd6-a8ba-bbac810126d5" />

**SAM File – Alignment Output**  
After alignment using BWA-MEM, the reads were saved in SAM (Sequence Alignment/Map) format.  

The SAM file contains:  
Header information – describes the reference genome and the software used for alignment.  
@SQ → gives information about the reference sequence, including its name and length.  
@PG → records the program and version used to generate the alignment.  
Alignment records – each subsequent line represents one sequencing read and its alignment to the reference genome.  

Important SAM fields to understand:  

QNAME	Name/ID of the sequencing read  
FLAG	Indicates the alignment properties, such as strand information  
RNAME	Reference sequence/chromosome to which the read mapped  
POS	Starting position of the alignment on the reference  
MAPQ	Mapping quality/confidence  
CIGAR	Describes how the read aligns to the reference  
SEQ	Sequence of the read  
QUAL	Base quality scores  
NM	Number of differences between the read and reference  
MD	Information about mismatches/deletions  
AS	Alignment score  
XS	Secondary alignment score  

For example, in the alignment output:  

NC_000913.3 represents the E. coli reference chromosome, and a MAPQ value of 60 indicates a highly confident alignment. A CIGAR value such as 66M indicates that 66 bases of the read were aligned.  

The SAM file is a human-readable alignment format. It will later be converted to BAM, which is the binary and more compact version used for efficient downstream analysis.  

The SAM alignment file was converted into BAM (Binary Alignment/Map) format using samtools. BAM is the compressed binary version of SAM and is more efficient for storage and downstream analysis. The resulting aligned.bam file contains the same alignment information as the SAM file but in a compressed binary format.    

```bash
#To convert .sam to .bam
samtools view -S -b aligned.sam > aligned.bam
```
After converting the SAM file to BAM format, the BAM file was sorted based on the position of the reads on the reference genome using Samtools. Sorting the BAM file is important for further steps such as indexing, coverage analysis, and variant calling.

```bash
#To sort the .bam file
samtools sort aligned.bam > aligned_sorted.bam
```
After sorting the BAM file, the next step was indexing the .bam file. Indexing creates a .bai index file that acts like a map of the sorted BAM file. A sorted BAM contains reads arranged according to their genomic coordinates, and the index allows tools such as samtools and genome browsers like IGV to efficiently locate and retrieve reads from particular regions. 

```bash
#To sort the bam file
samtools index alignment_sorted.bam
```

## IGV Visualization
After indexing the sorted BAM file, I loaded it into IGV to visually inspect the alignment of sequencing reads against the E. coli reference genome and I have recorded my observations below. 

<img width="1424" height="994" alt="image" src="https://github.com/user-attachments/assets/341bf547-0deb-4bff-a7a0-a108e5b95511" />

<img width="1420" height="989" alt="image" src="https://github.com/user-attachments/assets/068972da-08af-4b56-8907-bcdf4bd8a95b" />

<img width="1912" height="936" alt="image" src="https://github.com/user-attachments/assets/21b8c9aa-958f-4302-bb80-22319e83c145" />

<img width="1899" height="959" alt="image" src="https://github.com/user-attachments/assets/0f93940a-3bcf-4800-b60d-4dc41d1065c3" />



At genomic position NC_000913.3:101,770, the reference sequence contains C, while one aligned read contains G, indicating a C → G mismatch. The selected read has a MAPQ of 60 and a base quality of QV 36, indicating high-confidence mapping and a high-quality base call. The selected read is 265 bp long, has a MAPQ of 60, QV 36, CIGAR 265M, and NM=1. The read is aligned to the reverse strand (FLAG 16). I also saw an Insertion of one of the base C which is again worth investigating because its a variant and IGV shows possible differences from the reference genome. However, the mismatch is observed in only one of the reads covering this position, while the remaining reads agree with the reference. Therefore, this observation is recorded as a read-level mismatch and not considered a confirmed variant. Examining multiple reads at the same genomic position is important before interpreting a mismatch as a potential variant.  

To better understand the different components displayed in IGV, I created the following labelled diagram.

<img width="1536" height="1024" alt="image" src="https://github.com/user-attachments/assets/8e120fdd-377d-4424-a38b-220065b40ecf" />

## Variant Calling
After visualizing the aligned reads in IGV, the next step was to identify potential variants in the E. coli genome. The samtools mpileup command was used to examine the bases present at each genomic position based on the aligned reads in the sorted BAM file. The resulting mpileup file contains position-wise information about the reference base, read depth, bases observed in the reads, and base-quality information. This information provides the evidence required for identifying potential SNPs and indels during the subsequent variant-calling step. So basically the purpose of this step is to summarize the sequencing reads at each genomic position and generate the input evidence for variant calling.

```bash
samtools mplieup -f GCF_000005845.2_ASM584v2_genomic.fna  aligned_sorted.bam > output.mpileup
```
<img width="1778" height="510" alt="image" src="https://github.com/user-attachments/assets/e0df540f-8bf8-4ce7-ba8e-339270e9fbd6" />

The generated mpileup file was inspected using head to verify the output. The file contains position-wise information for the reference sequence, including the genomic position, reference nucleotide, read depth, bases observed in the aligned reads, and base-quality scores. The . and , symbols indicate bases matching the reference on the forward and reverse strands, respectively. This output provides the read-level evidence used for subsequent variant calling.

| Column | Example from output      | What it represents                                            |
| ------ | ------------------------ | ------------------------------------------------------------- |
| 1      | `NC_000913.3`            | Reference sequence/chromosome                                 |
| 2      | `1`                      | Genomic position                                              |
| 3      | `A`                      | Reference nucleotide at that position                         |
| 4      | `2`                      | Read depth/coverage — number of reads covering the position   |
| 5      | `^].^].`                 | Bases observed in the reads, along with alignment information |
| 6      | `EE`                     | Base-quality scores encoded using ASCII characters            |

| Symbol    | Meaning                                                 |
| --------- | ------------------------------------------------------- |
| `.`       | Read base matches the reference on the forward strand   |
| `,`       | Read base matches the reference on the reverse strand   |
| `A/C/G/T` | Base observed in a read that differs from the reference |
| `^`       | Beginning of a read/alignment segment                   |
| `$`       | End of a read                                           |
| `+`       | Insertion relative to the reference                     |
| `-`       | Deletion relative to the reference                      |

So, 
| Position | Reference | Depth | Observed bases | Interpretation                    |
| -------: | --------- | ----: | -------------- | --------------------------------- |
|        1 | A         |     2 | `^].^].`       | Read-start information is present |
|        2 | G         |     2 | `.,`           | Both reads match the reference    |
|        3 | C         |     2 | `.,`           | Both reads match the reference    |
|        4 | T         |     2 | `.,`           | Both reads match the reference    |
|        5 | T         |     2 | `.,`           | Both reads match the reference    |

After alignment and visualization of the reads in IGV, variant calling was performed using BCFtools. The bcftools mpileup command was used to generate genotype likelihoods from the aligned reads, and bcftools call was then used to identify potential variants such as SNPs and indels. Here, the reference genome was provided using -f, while the sorted BAM file containing the aligned reads was used as input. The output from bcftools mpileup was passed directly to bcftools call using the pipe (|).    
The -m option enables the multiallelic variant caller, while -v reports variant sites only. The -Ov option specifies VCF as the output format. The resulting variants were saved in variants.vcf. The generated VCF file contains the candidate variants identified from the sequencing data and will be used for subsequent variant filtering and validation using IGV.    

<img width="1852" height="97" alt="image" src="https://github.com/user-attachments/assets/23440e6b-88da-4024-9707-467861b8ce88" />

After aligning the sequencing reads to the E. coli reference genome and visualizing the alignments in IGV, the next step is to identify genomic positions that differ from the reference sequence. These differences may include single nucleotide variants (SNVs/SNPs) and small insertions or deletions (indels). 

```bash
bcftools mpileup -f GCF_000005845.2_ASM584v2_genomic.fna aligned_sorted.bam | bcftools call -mv -Ov > variants.vcf
cat variants.vcf
```
<img width="1905" height="641" alt="image" src="https://github.com/user-attachments/assets/3ebe050b-a835-40d2-b792-da54c9720886" />

| Column               | Meaning                                               | Example from your file |
| -------------------- | ----------------------------------------------------- | ---------------------- |
| `CHROM`              | Reference chromosome/sequence                         | `NC_000913.3`          |
| `POS`                | Genomic position of the variant                       | `49348`                |
| `ID`                 | Variant identifier, if available                      | `.`                    |
| `REF`                | Base/sequence present in the reference                | `.` / `G` / `C` etc.   |
| `ALT`                | Alternative base/sequence detected                    | `G`, `C`, `T`, etc.    |
| `QUAL`               | Quality/confidence score assigned to the variant call | `110.415`              |
| `FILTER`             | Filtering status                                      | `.`                    |
| `INFO`               | Additional information about the variant              | `DP=4;...`             |
| `FORMAT`             | Describes the sample-specific fields                  | `GT:PL`                |
| `aligned_sorted.bam` | Your sample's genotype information                    | `1/1:36,3,0`           |


The generated VCF file was inspected to understand the called variants. Each variant record contains information about its genomic position, reference allele (REF), alternative allele (ALT), quality score (QUAL), read depth and other supporting information. The results included both single-nucleotide variants (SNPs) and small insertions/deletions (indels). Since the initial variant calls may contain low-confidence calls, further filtering based on parameters such as read depth and variant quality, followed by visual validation in IGV, is required before considering them high-confidence variants.

So, to filter out the data 
```bash
 bcftools filter -i 'DP>=10 && QUAL>=20' filtered_variants.vcf
```
```bash
#to check how many variants are recorded in .vcf file
bcftools view -H filtered_variants.vcf | wc -l
```
<img width="1045" height="101" alt="image" src="https://github.com/user-attachments/assets/65efde4f-8cea-41ec-b7d8-92e1684dab42" />

Initial variant filtering using QUAL ≥20 and DP ≥10 resulted in no retained variants because all candidate calls had INFO-level DP below the predefined threshold. To avoid prematurely excluding potentially informative calls, variants meeting the QUAL criterion were retained as an exploratory candidate set. Twenty-one candidates were subsequently evaluated using read-level evidence in IGV, including coverage, alternate-allele support, strand representation, mapping quality, and local alignment context.

## Inspection of the output data
Because the filtering removed all the potential variants I wanted to check if any of them were actually worth reporting. The initial variant records were manually inspected to examine genomic position, reference and alternate alleles, variant quality, and read depth.

```bash
bcftools view -H variants.vcf | head
```
<img width="1456" height="485" alt="image" src="https://github.com/user-attachments/assets/2dad3d5c-2d0c-456f-9456-511396ec77fe" /> 

Here it showed depth to be the main reason why from the 124 potential variants none passed the filter. Then to get an idea about the whole data. Mean sequencing depth was calculated using samtools depth. The dataset showed an average depth of approximately 3.16×, indicating relatively shallow sequencing coverage. 

```bash
#To check the average sequencing depth. sum of all depths ÷ number of positions
 samtools depth aligned_sorted.bam | awk '{sum+=$3; count++} END {print "Average depth:", sum/count}'
```
<img width="1453" height="75" alt="image" src="https://github.com/user-attachments/assets/f7a09a73-6da1-452f-9159-d3d090962b31" />

Mean sequencing depth was calculated using samtools depth. The dataset showed an average depth of approximately 3.16×, indicating relatively shallow sequencing coverage. As the mean depth was so low I wanted to check the alignment statistics of the dataset and vcf file.

```bash
#Alignment QC
samtools flagstat aligned_sorted.bam
```
<img width="954" height="394" alt="image" src="https://github.com/user-attachments/assets/8965f2d2-b1b9-4937-9eef-4b2edecb93ae" />

Alignment statistics were assessed using samtools flagstat. Of the 74,716 reads, 74,559 (99.79%) mapped to the reference genome, indicating excellent alignment efficiency. No duplicate reads were detected. As this came out fine, I wanted to check the Coverage Statistics

```bash
#Coverage QC
samtools coverage aligned_sorted.bam
```

<img width="1453" height="67" alt="image" src="https://github.com/user-attachments/assets/752b8822-d9b5-479b-8dcb-d3874f0467d4" />

Genome-wide coverage was assessed using samtools coverage. Approximately 93.09% of the reference genome was covered, with a mean depth of 2.94×. The mean base quality (41.7) and mapping quality (58.5) were high, indicating good read and alignment quality despite the relatively low sequencing depth. Then as the coverage and alignment showed no problem I wanted to check the distribution of sequencing depth. 

```bash
 samtools depth aligned_sorted.bam | awk '{print $3}' | sort -n | uniq -c
```
<img width="1315" height="465" alt="image" src="https://github.com/user-attachments/assets/8dde0056-7ff0-4bdd-b918-8867504f6828" />

So, this made the whole picture much clearer.   
773,753 positions have 1 read covering them.   
1,034,427 positions have 2 reads  
955,324 positions have 3 reads  
691,252 positions have 4 reads  
This shows that most of your genome has very shallow coverage. The distribution of per-base sequencing depth was examined to assess coverage uniformity. Most genomic positions were covered at approximately 1–4× depth, confirming that the dataset provides broad but shallow genome coverage.   

Then in the total number of 124 variants I wanted to know how many of them are SNPs and Indels.   
```bash
#to count SNP variants
bcftools view -H -v snps variants.vcf | wc -l

#To count indels
bcftools view -H -v indels variants.vcf | wc -l
```

<img width="1018" height="95" alt="image" src="https://github.com/user-attachments/assets/0ec4777c-8580-4465-b823-fb84b152bbc4" />

The 124 candidate variants consisted of 87 SNPs (70.2%) and 37 indels (29.8%), indicating that SNPs represented the majority of the preliminary variant calls. Then I performed a Variant-level quality assessment.   
After identifying 124 candidate variants, their sequencing depth (DP) and variant quality (QUAL) were assessed to understand the reliability and confidence of the preliminary variant calls. The mean and range of these values were calculated to evaluate the overall support for the detected variants and the variation in quality among individual calls.  

The mean sequencing depth and mean QUAL score were calculated to obtain an overall assessment of the read support and quality of the 124 candidate variants.

```bash
 bcftools query -f '%QUAL\t%INFO/DP\n' variants.vcf |
awk '{sumDP+=$2; sumQ+=$1; n++} END {print "Variants:",n, "Mean DP:",sumDP/n, "Mean QUAL:",sumQ/n}'
```
<img width="1226" height="70" alt="image" src="https://github.com/user-attachments/assets/f9f5c613-b085-4e65-9ee9-afd1d594c6f7" />

Result: Mean DP = 2.83×, Mean QUAL = 18.23

The minimum and maximum DP and QUAL values were then calculated to determine the variation in sequencing support and confidence across individual candidate variants.

```bash
bcftools query -f '%QUAL\t%INFO/DP\n' variants.vcf |
awk 'BEGIN{minDP=999999; maxDP=0; minQ=999999; maxQ=0}
{
if($2<minDP)minDP=$2;
if($2>maxDP)maxDP=$2;
if($1<minQ)minQ=$1;
if($1>maxQ)maxQ=$1
}
END{
print "DP range:",minDP,"-",maxDP;
print "QUAL range:",minQ,"-",maxQ
}'
```

<img width="1069" height="328" alt="image" src="https://github.com/user-attachments/assets/b0ecaa83-78f2-474e-9420-5ea8f1841d1f" />

Result: DP = 1–12×, QUAL = 3.18–153.42

The candidate variants showed relatively low read support, with a mean depth of 2.83× and a range of 1–12×. Variant quality also varied considerably, with QUAL scores ranging from 3.18 to 153.42. This assessment was useful for determining the confidence of the preliminary calls and evaluating the effect of subsequent filtering criteria indicating that many candidate calls were supported by relatively few reads.

Since none of the 124 preliminary variant calls satisfied both the DP ≥10 and QUAL ≥20 criteria, the candidate calls were further examined to determine whether any variants met either criterion individually. This was done to assess whether the absence of filtered variants was due to both thresholds being simultaneously restrictive, particularly given the low overall sequencing depth of the dataset.

```bash
#with depth as filter
bcftools query -f '%QUAL\t%INFO/DP\n' variants.vcf |
awk '$2 >= 10 {n++} END {print "Variants with DP >=10:",n}'

#with qual as filter
bcftools query -f '%QUAL\t%INFO/DP\n' variants.vcf |
awk '$1 >= 20 {n++} END {print "Variants with QUAL >=20:",n}'

#with both as filter
bcftools query -f '%QUAL\t%INFO/DP\n' variants.vcf |
awk '$2 >= 10 && $1 >= 20 {n++} END {print "Variants with DP >=10 AND QUAL >=20:", n+0}'
```
Results: DP ≥10: 1 candidate
QUAL ≥20: 21 candidates
Both: 0 candidates

<img width="1092" height="278" alt="image" src="https://github.com/user-attachments/assets/52f29915-0247-47cb-b632-455ebe2107e8" />

A total of 21 variants were identified during variant calling. Since reporting every detected variant was unnecessary, variants were prioritized based on variant quality (QUAL), read depth (DP), alternate-allele support (DP4), mapping quality (MQ), and strand balance. Variants with higher confidence, adequate read support, high mapping quality, and alternate alleles supported on both sequencing strands were prioritized for downstream IGV validation. Variants located within repetitive homopolymer regions or showing limited/strand-biased read support were considered lower-confidence candidates. Based on these criteria, variants at positions 70289, 705013, 2844010, 4602509, and 4616669 were selected for detailed visualization in IGV because they showed relatively high QUAL scores, high mapping quality (MQ = 60), complete or strong alternate-allele support, and support from both forward and reverse reads. 

```bash
 { echo -e "CHROM\tPOS\tREF\tALT\tQUAL\tDP\tDP4\tMQ"; bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/DP\t%INFO/DP4\t%INFO/MQ\n' variants.vcf | awk '$5 >= 20'; } > qual20_candidates.tsv
```

<img width="861" height="530" alt="image" src="https://github.com/user-attachments/assets/ee9656ff-fddc-4ab1-b855-15a94555c6f7" />


I then checked the variants on IGV for the interpretation.     

**1. Variant One:**

<img width="860" height="54" alt="image" src="https://github.com/user-attachments/assets/489d3e6c-0c54-4ac2-ad14-4390e496e52d" />

<img width="1683" height="883" alt="image" src="https://github.com/user-attachments/assets/b0d449b8-2cff-480a-a1d7-b9a9f69cdd2f" />

The variant at NC_000913.3:70289 (G>T) had a variant quality (QUAL) of 110.415 and a read depth (DP) of 4. IGV inspection showed that all four reads covering the position supported the alternate T allele. The alternate allele was observed on both forward and reverse strands (2 reads on each strand). The mapping quality was 60. The variant therefore showed consistent read-level support in IGV, although the low coverage (DP=4) limits the confidence of the call. So its supported by the available reads, but requires caution due to low sequencing depth.    


**2. Variant two:**

<img width="863" height="50" alt="image" src="https://github.com/user-attachments/assets/7a93b0a0-a9d0-4350-9ebd-07e2d6c57748" />

<img width="1907" height="992" alt="image" src="https://github.com/user-attachments/assets/94111319-f5e6-46ef-bc24-f2009a090919" />

The variant at NC_000913.3:705013 (T>C) had a QUAL score of 115.415, read depth of 4, and mapping quality of 60. IGV inspection showed that all four reads covering the position supported the alternate C allele. The alternate allele was observed on both strands, with 3 forward-strand and 1 reverse-strand reads. No reads supporting the reference T allele were observed. The variant shows consistent read-level support with high mapping quality and representation on both strands. However, the low coverage (DP=4) limits confidence in the call. It was therefore considered a promising candidate requiring further validation. 


**3. Variant three:**

<img width="864" height="53" alt="image" src="https://github.com/user-attachments/assets/1e8b3374-6fb4-492c-b6c2-d6424bc853aa" />

<img width="1917" height="989" alt="image" src="https://github.com/user-attachments/assets/db9347e2-d8bc-4136-9e83-17c8d1568319" />

NC_000913.3:2,844,010 G>T. The variant was supported by all four reads covering the position (DP=4), with three forward-strand and one reverse-strand reads supporting the alternate T allele (DP4=0,0,3,1). The variant showed a high QUAL score of 99.42 and a mapping quality of 60, indicating good alignment quality. The variant was therefore selected for further validation using IGV.

**4. Variant four:**

<img width="859" height="49" alt="image" src="https://github.com/user-attachments/assets/7526e51a-779e-40c6-8f1b-31fe7e633276" />

<img width="1896" height="967" alt="image" src="https://github.com/user-attachments/assets/5c9e106e-8b23-49a7-9eb9-6fd0f2ba2e51" />

NC_000913.3:4,602,509 C>T. The variant was supported by all four reads covering the position (DP=4), with three forward-strand and one reverse-strand reads supporting the alternate T allele (DP4=0,0,3,1). The variant showed a high QUAL score of 113.42 and a mapping quality of 60, indicating strong variant-calling and read-mapping evidence. Therefore, this variant was selected for further validation using IGV.

**5. Variant five:**

<img width="859" height="54" alt="image" src="https://github.com/user-attachments/assets/090fa41c-a24b-421b-b1b6-60370b604a19" />

<img width="1900" height="966" alt="image" src="https://github.com/user-attachments/assets/c9e792ec-dd6f-47ef-9389-ba48d41eae86" />

NC_000913.3:4,616,669 G>T. The variant showed the highest QUAL score among the prioritized candidates (QUAL=153.42) and was supported by all five reads covering the position (DP=5). The alternate T allele was supported by three forward-strand and two reverse-strand reads (DP4=0,0,3,2), indicating balanced strand support. The mapping quality was 60, indicating high-quality read alignment. This variant was therefore prioritized for IGV validation.































