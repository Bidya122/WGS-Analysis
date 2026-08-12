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




