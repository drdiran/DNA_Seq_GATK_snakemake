########################################################
##### README file for DNA-Seq GATK workflow ############
##### Developer: Oyediran Akinrinade, PhD ##############
########################################################

The pipeline has two parts each with one submission script
- 1. GATK_workflow_V1.1.snakefile (submit_batch_GATKworkflow_V1.1_50jobs.sh)
- 2. GATK_GenotypeGVCF_MultiSample_Parallel_V1.1.snakefile (submit_GATK_MultiSample_genotypeGVCF_parallel_V1.1.sh)

## GATK_workflow_V1.1.snakefile / GATK_workflow_BamInput_V1.1.snakefile
This accepts preprocessed fastq files as input. The other variant (GATK_workflow_BamInput_V1.1.snakefile)
accepts bam files as input. The workflow is based on GATK/3.7.0. It aligns paired end reads to build 37
of the human genome (hg19) and follows the standard GATK best practices workflow.
This script generates gvcf file for cohort analysis. Some parameters need to be changed if
generating a multi-sample vcf file is not our goal. The snakefile generates script for each sample
and submit the job to the cluster.
# How to use: Run
bash path_to_script/submit_batch_GATKworkflow_V1.1_50jobs.sh &

to run 50 jobs at a time. The parameter 'j' in the 'submit_batch_GATKworkflow_V1.1_50jobs.sh' script
could be modified to submit more or less jobs at a time.
NB: Check the file 'GATK_workflow_V1.1.log' in the CWD for progress/updates on the jobs.
- Please rename this file if you have to resume the job at any point.

 ## GATK_GenotypeGVCF_MultiSample_Parallel_V1.1.snakefile
 This file takes a gvcf list file as input and runs GATK genotypeGVCFs in parallel per chromosome to generate
 vcf files for each chromosome. The chromosome level vcf files are then merged to generate a genome-wide
 vcf file for variants quality score recalibration (VQSR) and variants failing VQSR
 are subsequently flagged/filtered before annotation with any tool of choice (VEP, snpeff, annovar).
 (see VEP_annotation_snakemake folder for annotation script).
