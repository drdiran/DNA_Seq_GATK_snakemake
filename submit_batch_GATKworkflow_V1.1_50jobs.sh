#!/bin/sh
module load python/3.5.2
module load R/3.2.0

#
work_dir=`pwd`
snakefile="GATK_workflow_V1.1.snakefile"
## Make project result directories
mkdir FASTQ
mkdir BAM
mkdir LOGFILES
mkdir Raw_gVCF
mkdir Recal_Data

chmod -R 770 $work_dir

## Execute SnakeMake , GATK pipeline

snakemake \
 --nolock \
 --jobname 'SM.{rulename}.{jobid}' \
 --directory $work_dir \
 --snakefile $snakefile \
 -k -p -w 10 \
 --rerun-incomplete \
 --stats $work_dir/GATK_workflow_V1.stats \
 -j 50 \
 --timestamp \
 --cluster "qsub -l vmem=30g -l nodes=1:ppn=8 -l walltime=48:00:00 \
 -e LOGFILES -o LOGFILES" >& GATK_workflow_V1.1.log
