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
⚠️ Note on Read Layout (Learning Context)   
During this exercise, I ran Trimmomatic in single-end (SE) mode on a dataset that was originally paired-end.
I made this choice by mistake initially, but I continued the workflow to learn Nextflow, FASTQC, BWA alignment, and BAM handling. The purpose of this project is educational, not to produce publication-grade results.
The workflow is designed as a learning exercise to practice building reproducible bioinformatics pipelines.  
The main goal of this project was to learn and demonstrate a complete WGS pipeline, preparing me to scale this workflow to larger datasets or different organisms in future analyses. This workflow can be easily adapted to true paired-end processing by modifying the Trimmomatic and alignment steps to handle read pairs.  

🚀 Pipeline Overview

The workflow performs the following steps:  
Quality trimming – Trimmomatic  
Quality control – FastQC  
Reference indexing – BWA  
Read alignment – BWA-MEM  
BAM processing – sorting, indexing, merging  
Coverage & read statistics – samtools  
Variant calling – bcftools  
Variant annotation – SnpEff  
Consensus genome generation – bcftools consensus 
Clade & lineage assignment – Nextclade  

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
## NEXTCLADE
```bash
conda activate base   # or a tools env
```
```bash
conda install -c bioconda -c conda-forge nextclade
```
```bash
nextclade –version
```
```bash
#dataset download
nextclade dataset get --name sars-cov-2 --output-dir nextclade_dataset
```
<img width="940" height="105" alt="image" src="https://github.com/user-attachments/assets/d94310c2-e7c7-41ab-996d-49e4ed0dcc97" />  


# NEXTFLOW.CONFIG 
nextflow.config is the core configuration file that tells Nextflow how to run the pipeline. I have uploaded the config file [Find config file](https://github.com/Bidya122/WGS-Analysis/tree/main/Sars_cov2_1sample/Nextflow)

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
 [Find nf file](https://github.com/Bidya122/WGS-Analysis/tree/main/Sars_cov2_1sample/Nextflow)

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

**Let's see actually how does my file looks like on IGV:**

Open the IGV > On TOP right click on "Genomes" > Load genome from file > Select the .fasta file (reference genome), in this case its NC_045512.2 > You can see that on the bottom of the page > Then again on the top right click on "File" > Load from file > Select the .bam file. And this is what I see: 

<img width="1365" height="728" alt="image" src="https://github.com/user-attachments/assets/2ef64cde-604c-432f-aae0-25f45f533f47" />  

The aligned sequencing reads generated from the WGS Nextflow pipeline were visualized using Integrative Genomics Viewer (IGV) to assess alignment quality, coverage, and sequence variation.
*Reference Genome*  
The alignment was performed against the NC_045512.2 reference genome (SARS-CoV-2). The same FASTA file used during alignment was loaded into IGV to ensure coordinate and chromosome name consistency.  

*Input Files*  
The following files were used for visualization:  
Sorted BAM file (trimmed_SRR13182925_1.sorted.bam) containing aligned reads  
BAM index file (.bai) enabling rapid random access  
Reference genome FASTA file (NC_045512.2.fasta)  

*Coverage Track*  
The coverage track displays read depth across the genome, indicating the number of reads aligned at each genomic position. Continuous and uniform coverage across most regions suggests successful whole-genome sequencing and alignment.  

*Alignment Track*
Individual sequencing reads are shown as horizontal bars aligned to the reference genome. Stacked reads indicate regions of high coverage. Each read represents a trimmed sequencing fragment mapped to the reference.  

*Variant Representation*
Colored bases within reads represent mismatches between the aligned read and the reference sequence, potentially indicating single nucleotide variants (SNVs). Insertion and deletion events are visualized as small symbols within the reads, reflecting gaps or additional bases relative to the reference. If many reads show the same color at the same position, that is a real variant.  

*Genome-wide View*  
At the genome-wide scale (~29 kb), reads appear densely packed. Zooming into smaller genomic regions allows base-level inspection of alignments and clearer identification of sequence variations.  

<img width="1363" height="725" alt="image" src="https://github.com/user-attachments/assets/a264c900-27a3-47b1-8571-6cd096210595" />    

Base-level IGV view showing a single nucleotide variant (SNV). Multiple aligned reads consistently support an A substitution relative to the NC_045512.2 reference, indicating a high-confidence variant.  

<img width="211" height="224" alt="image" src="https://github.com/user-attachments/assets/2bcee49c-e80f-4914-b6fc-e28c51ebbdfe" />    

When the green A nucleotide is clicked in IGV, a base count summary is displayed for the selected genomic position (NC_045512.2:21,846). This summary reports a total read depth of 2252, with 2251 reads supporting an A substitution. Strand-specific counts show balanced support from both forward (1094+) and reverse (1157−) orientations, indicating the variant is not strand biased. Minimal support for other bases and low deletion counts suggest high confidence in this single nucleotide variant.

Although the reads were initially processed in single-end mode due to an earlier oversight, the resulting BAM files were merged to generate a complete sample-level alignment and to understand downstream viral genomics analyses.

```bash
process MERGE_BAMS {
    publishDir "${params.outdir}/merged", mode: 'copy'

    input:
    path bam_files  // This will be a list of BAMs

    output:
    tuple path("merged.bam"), path("merged.bam.bai"), path("total_reads.txt"), path("coverage.txt"), emit: merged_info

    script:
    """
    # Merge BAMs
    samtools merge merged.bam ${bam_files.join(' ')}
    samtools index merged.bam

    # Count total reads
    samtools view -c merged.bam > total_reads.txt

    # Compute coverage
    samtools depth merged.bam > coverage.txt
    """
}
```

<img width="940" height="402" alt="image" src="https://github.com/user-attachments/assets/a41cf695-03de-4a8f-b3ed-ced309ebec1f" />  

<img width="863" height="300" alt="image" src="https://github.com/user-attachments/assets/2d352856-db94-440f-8554-e8966334cd3e" />


<img width="940" height="354" alt="image" src="https://github.com/user-attachments/assets/31fd44a0-fbb7-453b-8be1-3104e212cab8" />  


| Metric                                  | Value   | Interpretation                                                               |
| --------------------------------------- | ------- | ---------------------------------------------------------------------------- |
| Total reads                             | 291,246 | Total number of reads in the merged BAM file                                 |
| Mapped reads                            | 265,884 | Reads successfully aligned to the reference genome                           |
| Mapping rate                            | 91.29%  | High alignment efficiency, suitable for downstream analysis                  |
| Supplementary alignments                | 665     | Likely split or partially aligned reads, common in viral amplicon sequencing |
| Paired reads                            | 0       | Dataset is treated as single-end                                             |
| Properly paired reads                   | N/A     | Not applicable for single-end data                                           |
| Reads with mate on different chromosome | 0       | Expected for a single-chromosome viral genome                                |  


The alignment results indicate high-quality mapping of sequencing reads to the SARS-CoV-2 reference genome (NC_045512.2), with over 91% of reads successfully aligned. This level of alignment accuracy is sufficient for reliable detection of single nucleotide variants (SNVs) and small indels.  
Paired-end–specific metrics are reported as 0 or N/A because the BAM file does not contain read-pair information and is treated as single-end. This does not affect variant calling, as variant detection relies primarily on:  
- Read depth  
- Base quality  
- Mapping quality  
- Consistency of variant support across reads  
Supplementary alignments likely arise from primer boundaries or fragmented reads and are expected in amplicon-based viral sequencing workflows.

<img width="940" height="471" alt="image" src="https://github.com/user-attachments/assets/4861afb8-74e9-4251-a235-04da3e1caa86" />  

```bash
process VARIANT_CALLING {
    tag "$sample_id"

    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
    tuple path(bam), path(bai), path(reference), val(sample_id)

    output:
    path "${sample_id}.variants.vcf.gz"

    script:
    """
    echo "Running variant calling for $sample_id"
    ls -lh  #show what’s inside work dir for debugging
    echo "Reference file is: $reference"
    bcftools mpileup -Ou -f $reference $bam | \
    bcftools call -mv -Oz -o ${sample_id}.variants.vcf.gz

    """
}
```
The resulting VCF files serve as the basis for downstream variant annotation, consensus genome generation, and lineage/clade assignment. Below is a representation of the results.  

<img width="940" height="234" alt="image" src="https://github.com/user-attachments/assets/1eb7d0dc-06c8-4f59-957a-6eef62091012" />

<img width="940" height="176" alt="image" src="https://github.com/user-attachments/assets/8d7fe042-bc9e-4983-9e5c-c390799da905" />


[Find the vcf file here](https://github.com/Bidya122/WGS-Analysis/tree/main/Sars_cov2_1sample/Output)

The generated VCF files capture high-confidence genomic variants supported by read depth and quality metrics, enabling downstream annotation and biological interpretation. The variant calling step produces a VCF (Variant Call Format) file containing genomic positions where the sample sequence differs from the reference genome.

**VCF Columns Explained**  

| Column     | Description                                                               |
| ---------- | ------------------------------------------------------------------------- |
| **CHROM**  | Reference sequence / chromosome name (e.g., *NC_045512.2* for SARS-CoV-2) |
| **POS**    | Genomic position of the variant (1-based coordinate)                      |
| **REF**    | Reference nucleotide at this position                                     |
| **ALT**    | Alternate nucleotide observed in the sample                               |
| **QUAL**   | Phred-scaled quality score indicating confidence in the variant call      |
| **FILTER** | Filter status (`.` indicates the variant passed all filters)              |
| **INFO**   | Additional annotations describing read depth, quality, and allele support |

**Key INFO Field Annotations (Commonly Used)**  

| INFO Tag | Meaning                                                                |
| -------- | ---------------------------------------------------------------------- |
| **DP**   | Total read depth at the variant position                               |
| **DP4**  | Read counts supporting REF and ALT alleles (forward & reverse strands) |
| **MQ**   | Root mean square mapping quality of reads                              |
| **MQSB** | Mapping quality strand bias                                            |
| **VDB**  | Variant distance bias (lower values indicate better support)           |
| **SGB**  | Segregation bias (used to detect alignment artefacts)                  |
| **AC**   | Alternate allele count                                                 |
| **AN**   | Total number of alleles considered                                     |
| **AF**   | Allele frequency of the variant                                        |
| **HOB**  | Homozygosity bias (useful in viral genomics)                           |
| **MQ0F** | Fraction of reads with zero mapping quality                            |

Variant calling identified multiple high-confidence single nucleotide variants across the viral genome, all supported by high read depth (DP > 80) and mapping quality (MQ = 60). Most variants showed allele frequencies close to 1, indicating fixation within the sample population. These variants were subsequently used for consensus genome generation and lineage assignment, suggesting that the sample represents a well-resolved viral genome suitable for downstream epidemiological and functional analysis.  

After obtaining high-confidence variants supported by sufficient read depth and quality, the next logical step is to translate raw nucleotide changes into biological context. Variant annotation enables:  
- Mapping variants to specific genes or ORFs  
- Predicting amino acid changes caused by coding variants  
- Classifying variants based on predicted functional impact  
This transformation is essential for downstream interpretation, such as functional analysis, lineage characterization, and comparative genomics.

```bash
process SNPEFF_ANNOTATION {
    publishDir "${params.outdir}/annotation", mode: 'copy'

    conda 'bioconda::snpeff=5.4.0a'

    input:
    tuple path(vcf), val(sample_id)

    output:
    path "merged.ann.vcf"
    path "*.html"
    path "*.txt"

    script:
    """
    echo "Running SnpEff annotation for $sample_id"

    java -Xmx4g -jar /home/bidya122/miniconda3/share/snpeff-5.4.0a-0/snpEff.jar \
        -stats snpeff_stats.html \
        NC_045512.2 \
        $vcf > merged.ann.vcf
    """
}
```

<img width="940" height="471" alt="image" src="https://github.com/user-attachments/assets/42a286ee-1ecf-42c0-806a-082d085528a8" />  

SnpEff annotation revealed that the majority of detected variants were synonymous or missense mutations within coding regions of the viral genome. Several missense variants were identified in functionally relevant genes, suggesting potential implications for protein structure or viral fitness. All annotated variants were supported by high sequencing depth and quality metrics.  

<img width="2174" height="801" alt="image" src="https://github.com/user-attachments/assets/36296fbc-8841-452c-9178-bd16193d58f9" />  

SnpEff-based functional annotation revealed that the detected variants were distributed across multiple viral genes, with the highest number observed in ORF1ab, consistent with its large genomic size. Most variants were classified as low-impact or modifier mutations, indicating synonymous or non-coding changes. Notably, several moderate-impact (missense) variants were identified in functionally relevant genes, including the Spike (S) protein and ORF3a, suggesting potential effects on protein function or viral phenotype. Overall, the variant profile reflects a well-resolved viral genome suitable for downstream functional and epidemiological analyses.

[Find the snpeff output here](https://github.com/Bidya122/Bulk_RNAseq_Analysis/tree/main/Data)

Variant calling and annotation identify individual nucleotide differences relative to a reference genome, but many downstream analyses require a complete genome sequence that represents the sample itself. To enable such analyses, a consensus genome was generated by integrating high-confidence variants into the reference sequence.   
The consensus genome provides a single, continuous representation of the sample’s genome, reflecting all validated variants detected during variant calling. This step is essential because:  
- Lineage and clade assignment tools require a full genome sequence  
- Comparative and phylogenetic analyses operate on FASTA sequences, not VCFs  
- It allows direct comparison of the sample genome with other reference or public genomes
Low coverage regions were masked, This approach ensures that the consensus genome reflects only well-supported genomic positions, improving reliability for downstream interpretation. Low-coverage regions (<10×) were masked to avoid unreliable base calls, resulting in a high-quality consensus genome sequence.

<img width="940" height="352" alt="image" src="https://github.com/user-attachments/assets/09a485ff-e6d1-4215-b664-a4f31ac07f08" />  

<img width="940" height="194" alt="image" src="https://github.com/user-attachments/assets/c692bd6b-5ff6-4abb-a88a-0302f26bbbc2" />  

The consensus genome length was verified to ensure it does not exceed the reference genome length, confirming correct consensus generation.  
```bash
#Checks that the FASTA file contains a single sequence header, ensuring the genome was not accidentally concatenated.
grep "^>" merged.consensus.fasta
```
```bash
#Counts the total number of nucleotide bases in the consensus genome by excluding the FASTA header and removing line breaks.
grep -v ">" merged.consensus.fasta | tr -d '\n' | wc -c
```
```bash
#Calculates the exact genome length by summing all nucleotide characters while ignoring the FASTA header lines.
awk 'BEGIN{len=0} !/^>/{len+=length($0)} END{print len}' merged.consensus.fasta
```
Has (2675 / 29903) × 100 ≈ 8.94%   
Approximately 8.94% of the consensus genome consists of masked (N) bases, reflecting low-coverage regions that were intentionally excluded to maintain sequence reliability.

<img width="940" height="120" alt="image" src="https://github.com/user-attachments/assets/e553630e-eb86-4804-8a8a-9b9db5511be1" />

Consensus genomes were generated using bcftools consensus with regions having <10× coverage masked to Ns. The final consensus genome was 29,903 bp in length with ~9% ambiguous bases, reflecting uneven sequencing coverage.  

After generating a high-confidence consensus genome, the next step is to contextualize the sample within known viral diversity. While variant calling and annotation describe mutations at the nucleotide and gene level, they do not provide information about lineage assignment, clade classification, or overall genome quality.To address this, the consensus genome was analyzed using Nextclade, a widely used tool for SARS-CoV-2 quality control and clade assignment.  
  
```bash
process NEXTCLADE_QC {

    publishDir "${params.outdir}/nextclade", mode: 'copy'

    input:
    path consensus_fasta
    path nextclade_dataset

    output:
    path "nextclade.tsv"
    path "nextclade.json"
    path "nextclade.aligned.fasta"
    
    script:
    """
    nextclade run \
        --input-dataset ${nextclade_dataset} \
        ${consensus_fasta} \
        --output-tsv nextclade.tsv \
        --output-json nextclade.json \
        --output-fasta nextclade.aligned.fasta
    """
}
```
<img width="940" height="363" alt="image" src="https://github.com/user-attachments/assets/2c99257f-10b9-41c0-8aac-b914ba248afe" />  

| Feature                      | Observation              | Interpretation                                                                              |
| ---------------------------- | ------------------------ | ------------------------------------------------------------------------------------------- |
| Reference genome             | NC_045512.2 (Wuhan-Hu-1) | Standard reference genome for SARS-CoV-2 comparative analysis.                              |
| Clade (Nextstrain / WHO)     | 20A / B.1                | Early pandemic clade corresponding to the B.1 lineage.                                      |
| Overall QC score             | 77.37 (mediocre)         | Genome quality is moderate due to missing regions, but sufficient for clade assignment.     |
| Total substitutions          | 7                        | Limited divergence from the reference genome.                                               |
| Total deletions / insertions | 0                        | No indels detected; genome structure is preserved.                                          |
| Total missing bases          | 2,675 (~8.94%)           | Missing regions reflect low coverage and may impact fine-scale mutation analysis.           |
| Total non-ACGTN bases        | 0                        | No ambiguous base calls in covered regions, indicating reliable sequencing data.            |
| Amino acid substitutions     | 5                        | Includes substitutions in E, ORF1b, ORF3a, and Spike proteins (e.g., S:D614G).              |
| Frame shifts / stop codons   | 0                        | No disruptive mutations detected; coding sequences remain intact.                           |
| Private mutations            | 4 labeled substitutions  | Mutations not shared with the broader clade, possibly reflecting sample-specific variation. |
| Missing regions              | Multiple ORF ranges      | Partial CDS coverage may affect gene-level mutation interpretation.                         |
| Pangolin lineage             | B.1                      | Consistent with Nextstrain clade assignment.                                                |

[Find the Nextclade output here](https://github.com/Bidya122/WGS-Analysis/tree/main/Sars_cov2_1sample/Output)

1.	The sequenced genome belongs to clade 20A / lineage B.1, an early SARS-CoV-2 variant.  
2.	Overall sequence quality is moderate, with ~9% missing data, but no major structural mutations (frameshifts or stop codons) were found.  
3.	Detected nucleotide substitutions are few and consistent with early B.1 lineage markers (C3037T, C14408T, A23403G).  
4.	Observed amino acid changes in spike (S:T95N, S:D614G) and other ORFs may have functional relevance, consistent with lineage characteristics.  
5.	No insertions, deletions, or frameshifts were detected, indicating a largely intact genome.  
✅ Overall, the genome is a moderately complete early B.1 lineage SARS-CoV-2 sequence, with a few private mutations but no major disruptive changes.  
I identified 7 high-confidence SNPs, including lineage-defining mutations such as Spike D614G and ORF1b P314L, which place this genome within clade 20A (B.1). These variants are consistent with globally circulating strains and explain the observed phylogenetic placement.  

Although Nextclade provides placement information, a full phylogenetic tree was not generated in this project because only a single consensus genome was analyzed. Phylogenetic reconstruction is inherently comparative and requires multiple sequences to infer evolutionary relationships reliably. In future work, I plan to process larger cohorts of genomes and integrate phylogenetic tree construction (e.g., using IQ-TREE, UShER, or Nextstrain’s Augur) to uncover evolutionary dynamics among samples.



  
















