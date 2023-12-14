#!/usr/bin/env nextflow


/*
 * Workflow parameters
 */

params.reads = "$projectDir/test_data/data"
params.pedigree = "$projectDir/test_data/data/pedigree.ped"
params.outDir = "$projectDir/test_data/results"
params.indexDir = "$projectDir/test_data/ref"
params.genome = "$projectDir/test_data/ref/hg19_chr8.fasta"
params.knownSites = "$projectDir/test_data/ref/1000G_phase1.snps.high_confidence.b37.chr18.vcf.gz"


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


/*
 * Base recalibration using GATK4
*/

process BASERECALIBRATION {
    publishDir "${params.outDir}/GATK_recalibration", mode:'copy'

    input:
    tuple val(sampleID), path(sampleFile), path(sampleIndex), path(sampleMetrics)
    path(indexDir)
    path(knownSites)

    output: 
    tuple val(sampleID), file('*_recal_data.table')

    script:
    """
    /gatk/gatk BaseRecalibrator \
        --input ${sampleID}_dupl.bam \
        --output ${sampleID}_recal_data.table \
        --reference ${indexDir}/*.fasta \
        --known-sites ${indexDir}/${knownSites}
    """
}


process BQSR {
    publishDir "${params.outDir}/GATK_recalibration", mode:'copy'

    input:
    tuple val(sampleID), path(sampleFile), path(sampleIndex), path(sampleMetrics), path(sampleTable)
    path(indexDir)
    
    output:
    tuple val(sampleID), file('*_recal.bam'), file('*_recal.bai')

    script:
    """
    /gatk/gatk ApplyBQSR \
        --bqsr-recal-file ${sampleID}_recal_data.table \
        --input ${sampleID}_dupl.bam \
        --output ${sampleID}_recal.bam \
        --reference ${indexDir}/*.fasta \
        --create-output-bam-index true
    """
}


process INDEX_BAM {
    publishDir "${params.outDir}/GATK_recalibration", mode:'copy'

    input:
    tuple val(sampleID), path(bamFile)

    output:
    tuple val(sampleID), file('*_recal_sorted.bam'), file('*_recal_sorted.bai')

    script:
    """
    samtools sort ${sampleID}_recal.bam -O bam > ${sampleID}_recal_sorted.bam 
    samtools index ${sampleID}_recal_sorted.bam > ${sampleID}_recal_sorted.bai
    """
}


/* 
 * Pedigree-based Variant Calling info using GATK HaplotypeCaller
*/

process FAMILY_CALL {
    publishDir "${params.outDir}/GATK_variant_calling", mode:'copy'

    input:
    tuple val(sampleID), path(sampleFile), path(sampleIndex), path(sampleTable)
    path(pedFile)
    path(indexDir)

    output:
    tuple val(sampleID), file('*_familyvariants.g.vcf'), file('*_familyvariants.g.vcf.idx')

    script:
    """
    /gatk/gatk  HaplotypeCaller  \
        --input ${sampleID}_recal.bam \
        --output ${sampleID}_familyvariants.g.vcf \
        --pedigree ${pedFile} \
        --reference ${indexDir}/*.fasta \
        -ERC GVCF \
        --create-output-bam-index true \
        --create-output-variant-index true \
        --java-options "-Xmx8G"
    """
}


/*
 * Create cohort VCF 
 */ 

process COHORT_COMB {
    publishDir "${params.outDir}/GATK_variant_calling", mode:'copy'

    input:
    tuple val(sampleID), path(gVCFFile), path(gVCFIndex)
    path(indexDir)

    output:
    file('cohort.g.vcf.gz')
    file('cohort.g.vcf.gz.tbi')
    file('cohort.vcf.gz')
    file('cohort.vcf.gz.tbi')

    script:
    """
    /gatk/gatk  CombineGVCFs \
        --variant ${sampleID}_familyvariants.g.vcf \
        --reference ${indexDir}/*.fasta \
        --output cohort.g.vcf.gz
    /gatk/gatk GenotypeGVCFs \
        --reference ${indexDir}/*.fasta \
        --variant cohort.g.vcf.gz \
        --output cohort.vcf.gz
    """
}


/* 
 * Variant Calling using parallel FREEBAYES
 */

process FREEBAYES { 
    publishDir "${params.outDir}/FREEBAYES_variant_calling", mode:'copy'

    input:
    tuple val(sampleID), file(sampleFile), file(sampleIndex)
    path(indexDir)

    output:
    tuple val(sampleID), file('*_variants.vcf')

    script:
    """
    freebayes-parallel <(fasta_generate_regions.py ${indexDir}/*.fai 1000000) 5 \
            -f ${indexDir}/*.fasta ${sampleID}_recal.bam > ${sampleID}_variants.vcf
    """
}

println "Execution starts!"

workflow {

    reads_ch = Channel
                        .fromFilePairs("${params.reads}/*_{,R}{1,2}{,_001}.{fq,fastq}{,.gz}")
	    
    refGenome = file("${params.genome}")
	indexDir = Channel.value("${params.indexDir}")
    knownSites = file("${params.knownSites}")
    pedFile = file("${params.pedigree}")

    trim_ch = TRIMMING(reads_ch)
    fastqc_ch = FASTQC(reads_ch.join(trim_ch))
    MULTIQC(fastqc_ch.collect())
    align_ch = ALIGNMENT(trim_ch, refGenome, indexDir)
    convert_ch = SAM_TO_BAM(align_ch)
    rg_ch = ADD_READGROUPS(convert_ch, refGenome)
    dupl_ch = MARK_DUPLICATES(rg_ch)
    recal_ch = BASERECALIBRATION(dupl_ch, indexDir, knownSites)
    corr_ch = BQSR(dupl_ch.join(recal_ch), indexDir)
    fam_ch = FAMILY_CALL(recal_ch.join(corr_ch), pedFile, indexDir)
    cohort_ch = COHORT_COMB(fam_ch, indexDir)

    vcFB_ch = FREEBAYES(corr_ch, indexDir)
}
