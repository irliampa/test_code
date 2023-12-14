#!/usr/bin/env nextflow


/*
 * Workflow parameters
 */

params.reads = "$projectDir/test_data/data"
params.pedigree = "$projectDir/test_data/data/pedigree.ped"
params.outDir = "$projectDir/test_data/results"
params.indexDir = "$projectDir/test_data/ref"
params.genome = "$projectDir/test_data/ref/hg19_chr8.fasta"

 log.info """\
    VARIANT CALLING - N F   P I P E L I N E
    =======================================
    readsDir     : ${params.reads}
    outDir       : ${params.outDir}
    indexDir     : ${params.indexDir}
    genome       : ${params.genome}
    pedigree     : ${params.pedigree}
    """
    .stripIndent()


/*
 * Trim reads with trimmomatic
 */

process TRIMMING {
    publishDir "${params.outDir}/trimmed", mode:'copy'
	cpus 2

    input:
    tuple val(sampleID), path(reads)

    output:
    tuple val(sampleID), path('*_trim*')

    script:
    fq_1_paired = sampleID + '_1_trim.fastq.gz'
    fq_1_unpaired = sampleID + '_1un_trim.fastq.gz'
    fq_2_paired = sampleID + '_2_trim.fastq.gz'
    fq_2_unpaired = sampleID + '_2un_trim.fastq.gz'
    """
    mkdir trimmed
    trimmomatic PE -threads ${task.cpus} -phred33 \
                ${reads[0]} ${reads[1]} \
                $fq_1_paired $fq_1_unpaired \
                $fq_2_paired $fq_2_unpaired \
                HEADCROP:5 SLIDINGWINDOW:4:30 MINLEN:25
    """

}


/*
 * Check quality of reads files using fastQC 
 */

process FASTQC {
	publishDir "${params.outDir}/qc", mode: 'copy'
	cpus 2

	input:
	tuple val(sampleID), path(reads), path(trimmed)

	output:
	path "fastqc_${sampleID}"

	script:
    """
	mkdir fastqc_${sampleID}
	fastqc -q ${reads} -t ${task.cpus} -o fastqc_${sampleID}
    fastqc -q ${trimmed} -t ${task.cpus} -o fastqc_${sampleID}
	"""
}


/*
 * Collect fastQC results in one report using multiQC
 */

process MULTIQC {
    publishDir "${params.outDir}/qc", mode:'copy'
    
    input:
    path '*'
    
    output:
    path 'multiqc_report.html'
    path 'multiqc_data'
    
    script:
    """
    multiqc . 
    """
} 


/*
 * Align reads to reference FASTA with BWA MEM
 */

process ALIGNMENT {
    publishDir "${params.outDir}/aligned", mode:'copy'
    cpus 4

    input:
    tuple val(sampleID), path(trimmed)
    path(refGenome)
    path(indexDir)

    output:
    tuple val(sampleID), path('*.sam')

    script:
    """
    mkdir aligned
    bwa mem -M -P -t ${task.cpus} ref/${refGenome} ${trimmed[0]} ${trimmed[2]}  > ${sampleID}_aln.sam
    """
}


/*
 * Convert SAM files to sorted BAM files using samtools
 */

process SAM_TO_BAM {
    publishDir "${params.outDir}/aligned", mode:'copy'

    input:
    tuple val(sampleID), path(aligned)

    output:
    tuple val(sampleID), path('*_sorted.bam')

    script:
    """
    samtools view -Sb ${aligned} > ${aligned.baseName}.bam
    samtools sort ${aligned.baseName}.bam -O bam > ${aligned.baseName}_sorted.bam
    """
}


/*
 * Add Read Group information to BAM files
 */

process ADD_READGROUPS {
    publishDir "${params.outDir}/aligned", mode:'copy'

    input:
    tuple val(sampleID), path(sorted)
    path(refGenome)

    output:
    tuple val(sampleID), path('*_rg.bam') 

    script:
    """ 
    /gatk/gatk AddOrReplaceReadGroups \
        --INPUT ${sorted} \
        --OUTPUT ${sampleID}_sorted_rg.bam \
        --RGID ${sampleID} \
        --RGLB ${sampleID} \
        --RGPL ILLUMINA \
        --RGPU unit \
        --RGSM ${sampleID} \
        --SORT_ORDER coordinate \
        -R ${refGenome}
    """
}


process MARK_DUPLICATES {
    publishDir "${params.outDir}/aligned", mode:'copy'

    input:
    tuple val(sampleID), path(sampleFile)

    output:
    tuple val(sampleID), path('*_dupl.bam'), path('*dupl.bai'), path('*_dupl_metrics.txt')

    script:
    """ 
    /gatk/gatk MarkDuplicates \
        --INPUT ${sampleID}_sorted_rg.bam \
        --OUTPUT ${sampleID}_dupl.bam \
        --METRICS_FILE ${sampleID}_dupl_metrics.txt \
        --ASSUME_SORT_ORDER coordinate \
        --CREATE_INDEX true \
        --java-options '-Xmx6g'
    """
}


println "Execution starts!"

workflow {

    reads_ch = Channel
                        .fromFilePairs("${params.reads}/*_{,R}{1,2}{,_001}.{fq,fastq}{,.gz}")
	    
    refGenome = file("${params.genome}")
	indexDir = Channel.value("${params.indexDir}")

    trim_ch = TRIMMING(reads_ch)
    trim_ch.view()
    fastqc_ch = FASTQC(reads_ch.join(trim_ch))
    MULTIQC(fastqc_ch.collect())
    align_ch = ALIGNMENT(trim_ch, refGenome, indexDir)
    convert_ch = SAM_TO_BAM(align_ch)
    rg_ch = ADD_READGROUPS(convert_ch, refGenome)
    dupl_ch = MARK_DUPLICATES(rg_ch)

}
