#!/bin/sh
module load python/3.5.2
module load R/3.2.0

#
work_dir=`pwd`
snakefile="/hpf/largeprojects/smital/TOOLS/GATK_GenotypeGVCF_MultiSample_Parallel_V1.1.snakefile"
## Make project result directories
#mkdir FASTQ
#mkdir BAM
mkdir LOGFILES
#mkdir Raw_gVCF
#mkdir Recal_Data

chmod -R 770 $work_dir

## Execute SnakeMake , GATK pipeline

snakemake \
 --nolock \
 --jobname 'SM.{rulename}.{jobid}' \
 --directory $work_dir \
 --snakefile $snakefile \
 -k -p -w 10 \
 --rerun-incomplete \
 --stats $work_dir/GATK_genotypeGVCF_MultiSample.stats \
 -j 25 \
 --timestamp \
 --cluster "qsub -l vmem=40g -l nodes=1:ppn=8 -l walltime=48:00:00 \
 -e LOGFILES -o LOGFILES" >& GenotypeGVCF_MultiSample_parallel_V1.1.log

