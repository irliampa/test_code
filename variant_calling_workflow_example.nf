#!/usr/bin/env nextflow


/*
 * Workflow parameters
 */

params.reads = "$projectDir/data"
params.pedigree = "$projectDir/data/pedigree.ped"
params.outDir = "$projectDir/results"
params.indexDir = "$projectDir/ref"
params.genome = "$projectDir/ref/Homo_sapiens_assembly19.fasta"

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
                HEADCROP:15 SLIDINGWINDOW:4:20 MINLEN:25
    """

}


/*
 * Check quality of reads files using fastQC 
 */

process FASTQC {
	publishDir "${params.outDir}/qc", mode: 'copy'
	cpus 2

	input:
	tuple val(sampleID), path(reads)

	output:
	path "fastqc_${sampleID}"

	script:
    """
	mkdir fastqc_${sampleID}
	fastqc -q ${reads} -t ${task.cpus} -o fastqc_${sampleID}
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


println "Execution starts!"

workflow {

    reads_ch = Channel
                        .fromFilePairs("${params.reads}/*_{,R}{1,2}{,_001}.{fq,fastq}{,.gz}")
	    
    refGenome = file("${params.genome}")
	indexDir = Channel.value("${params.indexDir}")

    trim_ch = TRIM(reads_ch)
    fastqc_ch = FASTQC(reads_ch.join(trim_ch))
    align_ch = ALIGN(trim_ch, refGenome, index_ch)
    MULTIQC(fastqc_ch.collect())
    

}
