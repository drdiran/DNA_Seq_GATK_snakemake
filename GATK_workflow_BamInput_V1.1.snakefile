import glob, os


# Declare the resources being used
Human_Ref = "~/references/homo_sapiens/v37_decoy/gatk/human_g1k_v37_decoy_mod.fasta"
Human_Ref_bwa = "~/references/homo_sapiens/v37_decoy/bwa/human_g1k_v37_decoy_mod.fasta"
INDELS = "~/references/homo_sapiens/v37_decoy/gatk/1000G_phase1.indels.b37_mod.vcf"
MILLS = "~/references/homo_sapiens/v37_decoy/gatk/Mills_and_1000G_gold_standard.indels.b37_mod.vcf"
DBSNP = "~/references/homo_sapiens/v37_decoy/gatk/dbsnp_138.b37_mod.vcf"
HAPMAP = "~/references/homo_sapiens/v37_decoy/gatk/hapmap_3.3.b37_mod.vcf"
OMNI = "~/references/homo_sapiens/v37_decoy/gatk/1000G_omni2.5.b37_mod.vcf"
SNPS = "~/references/homo_sapiens/v37_decoy/gatk/1000G_phase1.snps.high_confidence.b37_mod.vcf"
EXCLUDED_CONTIGS = "~/TOOLS/Contigs_Excluded_b37.list"
# Tools directory
picard_dir = "~TOOLS/picard-tools/2.5.0"
gatk_dir = "~TOOLS/gatk/3.7.0"

# generate dataset
# This workflow assumes the bam files are located in the cwd
# if the bam files are in diff directories, one could create symbolic links of
# the bam files in the intended working directory.

SAMP = [os.path.basename(sample) for sample in glob.glob('*.bam')]
SAMPLES = [sample.split('_sorted.bam')[0] for sample in SAMP]


rule all:
     input: expand("BAM_dupMarked/{sample}_dupMarked.bam.bai", sample=SAMPLES),
            expand("Raw_gVCF/{sample}_raw.g.vcf.gz", sample=SAMPLES)

# This rule converts the bam file to fastq files in PE mode
rule samToFastq:
     input:"{sample}_sorted.bam"
     output:R1="FASTQ/{sample}/{sample}_R1.fastq",
            R2="FASTQ/{sample}/{sample}_R2.fastq",
            unpaired="FASTQ/{sample}/{sample}_unpaired.fastq"
     log: "LOGFILES/samToFastq/{sample}.log"
     shell:"""
           module load picard-tools/2.5.0
           java -Xmx15G -jar $picard_dir/picard.jar SamToFastq \
           I={input} \
           FASTQ={output.R1} \
           SECOND_END_FASTQ={output.R2} \
           UNPAIRED_FASTQ={output.unpaired} 2> {log}
           """
# This rule runs bwa alignment on the fastq files produced by rule SamToFastq
rule bwaAlignment:
     input: fwd="FASTQ/{sample}/{sample}_R1.fastq",
            rev="FASTQ/{sample}/{sample}_R2.fastq"
     params: bwa_ref=Human_Ref_bwa,
             tmp="FASTQ/{sample}/{sample}_tmp",
             rg="-R '@RG\tID:bwa\tSM:{sample}\tPU:Illumina\tLB:lib1\tPL:ILLUMINA'"
     output: "BAM/{sample}_sorted.bam"
     threads: 8
     log: "LOGFILES/bwaAlignment/{sample}.log"
     shell:"""
           module load bwa/0.7.15
           module load samtools/1.3.1
           (bwa mem {params.rg} -t {threads} {params.bwa_ref} {input.fwd} {input.rev} | \
           samtools view -S -b - | \
           samtools sort -T {params.tmp} -O bam -o {output} - ) 2> {log}
           rm -rf {input.fwd}
           rm -rf {input.rev}
           """

# Removed as it takes longer time. This has been included in the bwa mem .
#rule addOrReplaceReadGroups:
#     input:"BAM/{sample}_sorted.bam"
#     output:"BAM/{sample}_rg.bam"
#     params:"RGLB=lib1 RGPL=illumina RGPU={sample} RGSM={sample}"
#     log: "LOGFILES/addOrReplaceReadGroups/{sample}.log"
#     shell:"""
#           module load picard-tools/2.5.0
#           java -Xmx10G -jar /hpf/tools/centos6/picard-tools/2.5.0/picard.jar AddOrReplaceReadGroups \
#           {params} I={input} O={output} &> {log}
#           rm -rf {input}
#          """

rule markDuplicates:
     input: "BAM/{sample}_sorted.bam"
     output: bam="BAM_dupMarked/{sample}_dupMarked.bam",
             bai="BAM_dupMarked/{sample}_dupMarked.bam.bai",
             metrics="BAM_dupMarked/METRICS/{sample}_markedDupMetrics.txt"
     params: tmp="BAM/{sample}_tmp"
     log: "LOGFILES/markDuplicates/{sample}.log"
     shell:"""
           module load picard-tools/2.5.0
           module load samtools/1.3.1
           java -Xmx20G -jar $picard_dir/picard.jar MarkDuplicates \
           I={input} \
           O={output.bam} \
           M={output.metrics} TMP_DIR={params.tmp} 2> {log}
           samtools index {output.bam}
           rm -rf {input}
           rm -rf {params.tmp}
           """

####################### NOTE ################################################################
## I didn't include indel realignment step as it is not required with/for haplotype caller ##
#############################################################################################

# Baserecalibration is a 4-step process splitted into step1-step4 here in this workflow
# 1 -  Analyze patterns of covariation in the sequence dataset
rule baseRecalibrator_step1:
     input: "BAM_dupMarked/{sample}_dupMarked.bam"
     output: "Recal_Data/{sample}_recal_data.table"
     params: Ref=Human_Ref,
             dbsnp=DBSNP,
             kg_indels=INDELS,
             mills=MILLS
     log: "LOGFILES/baseRecalibrator_step1/{sample}.log"
     threads: 8
     shell:"""
           module load gatk/3.7.0
           java -Xmx15g -jar $gatk_dir/GenomeAnalysisTK.jar \
           -T BaseRecalibrator -R {params.Ref} -I {input} -knownSites {params.dbsnp} \
           -knownSites {params.kg_indels} -knownSites {params.mills} -o {output} -nct {threads} 2> {log}
           """

# 2 - Second pass to analyze covariation remaining after recalibration
rule baseRecalibrator_step2:
     input: bam="BAM_dupMarked/{sample}_dupMarked.bam",
            table="Recal_Data/{sample}_recal_data.table"
     output: "Recal_Data/{sample}_post_recal_data.table"
     params: Ref=Human_Ref,
             dbsnp=DBSNP,
             kg_indels=INDELS,
             mills=MILLS
     log: "LOGFILES/baseRecalibrator_step2/{sample}.log"
     threads: 8
     shell:"""
           module load gatk/3.7.0
           java -Xmx15g -jar $gatk_dir/GenomeAnalysisTK.jar \
           -T BaseRecalibrator -R {params.Ref} -I {input} -knownSites {params.dbsnp} \
           -knownSites {params.kg_indels} -knownSites {params.mills} -BQSR {input.table} \
           -o {output} -nct {threads} 2> {log}
           """

###################################################################################
# 3 - Generate before/after plots
# I have removed this part as it has some bugs. The step 4 doesn't depend on it

#rule analyzeCovariates_step3:
#     input: before="BAM/{sample}/{sample}_recal_data.table",
#            after="BAM/{sample}/{sample}_post_recal_data.table"
#     output: "BAM/{sample}/{sample}_recalibration_plots.pdf"
#     params: Ref=Human_Ref
#     log: "LOGFILES/AnalyzeCovariates/{sample}.log"
#     shell:"""
#           module load gatk/3.7.0
#           module load R/3.3.2
#           java -Xmx6g -jar /hpf/tools/centos6/gatk/3.7.0/GenomeAnalysisTK.jar \
#           -T AnalyzeCovariates -R {params.Ref} -before {input.before} \
#           -after {input.after} -plots {output} 2> {log}
#           """
###################################################################################
# 4 - Apply the recalibration to the sequence data
rule printReads_step4:
     input: bam="BAM_dupMarked/{sample}_dupMarked.bam",
            table="Recal_Data/{sample}_recal_data.table"
     output: "bwa_BAM/{sample}_recal.bam"
     params: Ref=Human_Ref
     log: "LOGFILES/PrintReads/{sample}.log"
     threads: 8
     shell:"""
           module load gatk/3.7.0
           module load R/3.3.2
           java -Xmx15g -jar $gatk_dir/GenomeAnalysisTK.jar \
           -T PrintReads -R {params.Ref} -I {input.bam} \
           -BQSR {input.table} -o {output} -nct {threads} 2> {log}
           rm -rf {input.bam}
           """

# This generates the gVCF file for cohort analysis workflow
rule HaplotypeCaller:
     input: "bwa_BAM/{sample}_recal.bam"
     output: "Raw_gVCF/{sample}_raw.g.vcf.gz"
     params: Ref=Human_Ref,
             dbsnp=DBSNP,
             call_conf=30,
             other_contigs=EXCLUDED_CONTIGS,
             maxALTalleles=2
     log: "LOGFILES/HaplotypeCaller/{sample}.log"
     threads: 8
     shell:"""
           module load gatk/3.7.0
           java -Xmx18g -jar $gatk_dir/GenomeAnalysisTK.jar \
           -T HaplotypeCaller -R {params.Ref} -I {input} --dbsnp {params.dbsnp} --genotyping_mode DISCOVERY \
           -stand_call_conf {params.call_conf} -XL {params.other_contigs} --emitRefConfidence GVCF -o {output} -nct {threads} --max_alternate_alleles {params.maxALTalleles} 2> {log}
           """
