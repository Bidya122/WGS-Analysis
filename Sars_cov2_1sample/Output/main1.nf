#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ---------------- PARAMETERS ----------------
params.samples = 'data/*.fastq.gz'
params.ref     = 'data/NC_045512.2.fasta'
params.outdir  = 'results'
params.snpEff_jar = "snpEff" 
params.nextclade_dataset = "nextclade_dataset"

  // 1️⃣ Define your input channels
    samples_ch = Channel.fromPath(params.samples)
    ref_ch     = Channel.fromPath(params.ref)
    nextclade_dataset_ch = Channel.fromPath(params.nextclade_dataset)

// ---------------- PROCESSES ----------------

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
    ls -lh  # ✅ show what’s inside work dir for debugging
    echo "Reference file is: $reference"
    bcftools mpileup -Ou -f $reference $bam | \
    bcftools call -mv -Oz -o ${sample_id}.variants.vcf.gz

    """
}

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

process CONSENSUS_GENOME {

    publishDir "${params.outdir}/consensus", mode: 'copy'

    input:
    tuple path(bam), path(bai), path(vcf), path(reference), val(sample_id)

    output:
    path "${sample_id}.consensus.fasta"

    script:
    """
    echo "Generating consensus genome for ${sample_id}"

    # Ensure VCF is indexed (required by bcftools consensus)
    bcftools index -t ${vcf}

    # Create BED of low coverage regions (<10x)
    samtools depth ${bam} | \
        awk '\$3 < 10 {print \$1"\t"\$2-1"\t"\$2}' > lowcov.bed

    # Generate consensus genome
    bcftools consensus \
        -f ${reference} \
        -m lowcov.bed \
        ${vcf} > ${sample_id}.consensus.fasta
    """
}

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





// ---------------- WORKFLOW ----------------


reference_ch = Channel.value(file("data/NC_045512.2.fasta"))

workflow {
    trimmed_ch = TRIMMOMATIC(samples_ch)
    FASTQC(trimmed_ch)

    indexed_ref_ch = BWA_INDEX(Channel.value(file("data/NC_045512.2.fasta")))

    align_input_ch = trimmed_ch.combine(indexed_ref_ch)
    aligned_ch = ALIGN_READS(align_input_ch)
    aligned_bams_ch = aligned_ch.map { it instanceof Tuple ? it[0] : it }

    merged_bams_ch = MERGE_BAMS(aligned_bams_ch)

    variant_input_ch = merged_bams_ch.map { bam, bai, total, cov ->
        tuple(bam, bai, file("data/NC_045512.2.fasta"), "merged")
    }

    // Call VARIANT_CALLING once
    variant_vcf_ch = VARIANT_CALLING(variant_input_ch)

    // Prepare for SnpEff
    annotated_ch = variant_vcf_ch.map { vcf -> tuple(vcf, "merged") }

    SNPEFF_ANNOTATION(annotated_ch)
   
consensus_input_ch = variant_vcf_ch.map { vcf ->
    tuple(
        file("results/merged/merged.bam"),
        file("results/merged/merged.bam.bai"),
        vcf,
        file(params.ref),  // Pass reference as input
        "merged"
    )
}

// Run consensus genome generation

consensus_ch = CONSENSUS_GENOME(consensus_input_ch)

NEXTCLADE_QC(
    consensus_ch,
    nextclade_dataset_ch
)


}






















