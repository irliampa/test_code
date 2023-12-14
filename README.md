# test_code  
This is a basic implementation of GATK best practices variant calling workflow in nextflow.  

## Benchmarking data  
A small family trio dataset was retrieved from the Galaxy communinty tutorial [**Exome sequencing data analysis for diagnosing a genetic disease**](https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/exome-seq/tutorial.html).  

The dataset contains partial paired-end whole exome sequencing FASTQ files from a family trio (father - mother - proband) and the respective pedigree file and reference files for chromosome 18 (hg18) to run the pipeline.  

## Workflow overview  
![alt text](https://github.com/irliampa/test_code/blob/main/test_data/pipeline_runtime/workflow.png)
This is the DAG of the workflow, automatically created when running the pipeline.  

## Workflow information  
This implementation is based on the basic [Germline short variant discovery (SNPs + Indels) pipeline](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-), but designed for family trio data and without the variants postprocessing steps, as described in the figure from the aforementioned post: 
![alt text](https://github.com/irliampa/test_code/blob/main/assets/01.png)  

In addition it contains variant calling using [freebayes variant caller](https://github.com/freebayes/freebayes), enabling comparison of the retrieved variants from both tools, i.e. GATK HaplotypeCaller and breebayes.  

## Toolset  
The workflow runs community maintained dockerized tools, from public container repositories.  

## How to run  
```
nextflow run -c nextflow.config -profile local variant_calling_workflow_example.nf --reads test_data/data/ --outDir test_data/results/ 
```
The runtime reports were created by adding to the previous command:  
```
-with-report test_data/pipeline_runtime/report.html -with-trace test_data/pipeline_runtime/trace.txt -with-timeline test_data/pipeline_runtime/timeline.html -with-dag test_data/pipeline_runtime/workflow.png
```
