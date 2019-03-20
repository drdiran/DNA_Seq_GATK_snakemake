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
dbsnp_147 = "~/references/homo_sapiens/hg19/dbSNP_147/All_20160601.vcf.gz"
dbsnp_150 = "~/references/homo_sapiens/hg19/dbSNP_150/GATK/All_20170403.vcf.gz"
# Tools directory
picard_dir = "~TOOLS/picard-tools/2.5.0"
gatk_dir = "~TOOLS/gatk/3.7.0"
snpEff="~/TOOLS/snpEff/4.3"
#gvcfList ="~/PCS2_GATK/cohortGVCF.list"

# generate dataset
# This workflow assumes the raw gVCF files are located in the cwd (e.g Raw_gVCF dir)
# if the gVCF files are in diff directories, one could create symbolic links of
# the gVCFs (without necessarily copying them) in the intended working directory.

SAMP = [os.path.basename(sample) for sample in glob.glob('*_raw.g.vcf.gz')]
SAMPLES = [sample.split('_raw')[0] for sample in SAMP]

CHROMOSOMES = 'chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM'.split()

rule all:
     input: expand("GenotypeGVCF/Combined_multisampleGenotype_{chr}.vcf.gz",chr=CHROMOSOMES),
            "MultiSampleVCF/Combined_multisampleGenotype_recalSNPs_recalINDELs_passed.vcf.gz",
            "MultiSampleVCF/Recalibration/recalibrate_SNP.plots.R",
            "MultiSampleVCF/Recalibration/recalibrate_INDEL.plots.R"


rule genotypeGVCFs:
     input: "cohortGVCF.list"
     output: "GenotypeGVCF/Combined_multisampleGenotype_{chr}.vcf.gz"
     params: Ref=Human_Ref,
             dbsnp=DBSNP,
             chrm="{chr}"
     log: "LOGFILES/GenotypeGVCFs.log"
     shell:"""
           module load gatk/3.7.0
           java -Xmx30g -jar $gatk_dir/GenomeAnalysisTK.jar \
           -T GenotypeGVCFs -R {params.Ref} --dbsnp {params.dbsnp} \
           -V {input} -o {output} --max_alternate_alleles 6 -L {params.chrm} 2> {log}
           """

rule joinVCF:
     input: expand("GenotypeGVCF/Combined_multisampleGenotype_{chr}.vcf.gz",chr=CHROMOSOMES)
     output: "MultiSampleVCF/Combined_multisampleGenotype_snpsift.vcf.gz"
     log: "joinVCF.log"
     shell: """
          module load snpEff/4.3
          module load java/1.8.0_91
          module load htslib/1.4.1
          java -jar -Xmx20g $snpEff/SnpSift.jar split -j {input} | bgzip -c > {output}
          tabix {output}
           """


# VARIANT QUALITY SCORE RECALIBRATOR
# Done separately for SNPs and INDELs

############## VQSR SNP  #################################
rule variantRecalibratorSNP:
     input: "MultiSampleVCF/Combined_multisampleGenotype_snpsift.vcf.gz"
     output: recalFile="MultiSampleVCF/Recalibration/recalibrate_SNP.recal",
             tranchesFile="MultiSampleVCF/Recalibration/recalibrate_SNP.tranches",
             rscriptFile="MultiSampleVCF/Recalibration/recalibrate_SNP.plots.R"
     params: Ref=Human_Ref,
             hapmap=HAPMAP,
             omni=OMNI,
             KG_snps=SNPS,
             dbsnp=DBSNP,
             maxGau=4
     threads: 8
     log: "LOGFILES/variantRecalibratorSNP.log"
     shell:"""
           module load gatk/3.7.0
           module load R/3.3.2
           java -Xmx20g -jar $gatk_dir/GenomeAnalysisTK.jar \
           -T VariantRecalibrator -R {params.Ref} -input {input} \
           -resource:hapmap,known=false,training=true,truth=true,prior=15.0 {params.hapmap} \
           -resource:omni,known=false,training=true,truth=true,prior=12.0 {params.omni} \
           -resource:1000G,known=false,training=true,truth=false,prior=10.0 {params.KG_snps} \
           -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {params.dbsnp} \
           -an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum \
           -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
           -recalFile {output.recalFile} -tranchesFile {output.tranchesFile} \
           -rscriptFile {output.rscriptFile} -nt {threads} --maxGaussians {params.maxGau} 2> {log}
           """

############## APPLY RECAL: SNP  #################################
rule applyRecalibrationSNP:
     input: raw_vcf="MultiSampleVCF/Combined_multisampleGenotype_snpsift.vcf.gz",
            recalFile="MultiSampleVCF/Recalibration/recalibrate_SNP.recal",
            tranchesFile="MultiSampleVCF/Recalibration/recalibrate_SNP.tranches"
     output: "MultiSampleVCF/Combined_multisampleGenotype_recalSNPs_rawINDELs.vcf.gz"
     params: Ref=Human_Ref,
             level=99.5
     threads: 8
     log: "LOGFILES/applyRecalibrationSNP.log"
     shell:"""
           module load gatk/3.7.0
           module load R/3.3.2
           java -Xmx20g -jar $gatk_dir/GenomeAnalysisTK.jar \
           -T ApplyRecalibration -R {params.Ref} -input {input.raw_vcf} \
           -mode SNP -recalFile {input.recalFile} -tranchesFile {input.tranchesFile} \
           --ts_filter_level {params.level} -o {output} -nt {threads} 2> {log}
           rm -rf {input.raw_vcf}
           """

###############################################################################
############## VQSR INDEL  #################################
## Input here is the gVCF file with recalibrated SNPs and raw indels
###############################################################################
rule variantRecalibratorINDEL:
     input: "MultiSampleVCF/Combined_multisampleGenotype_recalSNPs_rawINDELs.vcf.gz"
     output: recalFile="MultiSampleVCF/Recalibration/recalibrate_INDEL.recal",
             tranchesFile="MultiSampleVCF/Recalibration/recalibrate_INDEL.tranches",
             rscriptFile="MultiSampleVCF/Recalibration/recalibrate_INDEL.plots.R"
     params: mill=MILLS,
             dbsnp=DBSNP,
             Ref=Human_Ref,
             maxGau=4
     threads: 8
     log: "LOGFILES/variantRecalibratorINDEL.log"
     shell:"""
           module load gatk/3.7.0
           module load R/3.3.2
           java -Xmx20g -jar $gatk_dir/GenomeAnalysisTK.jar \
           -T VariantRecalibrator -R {params.Ref} -input {input} \
           -resource:mills,known=false,training=true,truth=true,prior=12.0 {params.mill} \
           -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {params.dbsnp} \
           -an DP -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum \
           -mode INDEL --maxGaussians {params.maxGau} -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
           -recalFile {output.recalFile} -tranchesFile {output.tranchesFile} \
           -rscriptFile {output.rscriptFile} -nt {threads} 2> {log}
           """

###############################################################################
############## APPLY RECAL: INDEL  ############################################
###############################################################################
rule applyRecalibrationINDEL:
     input: recalSNP_rawINDEL="MultiSampleVCF/Combined_multisampleGenotype_recalSNPs_rawINDELs.vcf.gz",
            recalFile="MultiSampleVCF/Recalibration/recalibrate_INDEL.recal",
            tranchesFile="MultiSampleVCF/Recalibration/recalibrate_INDEL.tranches"
     output: "MultiSampleVCF/Combined_multisampleGenotype_recalSNPs_recalINDELs.vcf.gz"
     params: Ref=Human_Ref,
             level=99.0
     threads: 8
     log: "LOGFILES/applyRecalibrationINDEL.log"
     shell:"""
           module load gatk/3.7.0
           module load R/3.3.2
           java -Xmx30g -jar $gatk_dir/GenomeAnalysisTK.jar \
           -T ApplyRecalibration -R {params.Ref} -input {input.recalSNP_rawINDEL} \
           -mode INDEL -recalFile {input.recalFile} -tranchesFile {input.tranchesFile} \
           --ts_filter_level {params.level} -o {output} -nt {threads} 2> {log}
           rm -rf {input.recalSNP_rawINDEL}
           """
###############################################################################

# 1. Remove filtered variants using snpEff too
# 2. Annotate with dbsnp rsID (dbsnp150)


rule removeFilteredVariants:
     input: "MultiSampleVCF/Combined_multisampleGenotype_recalSNPs_recalINDELs.vcf.gz"
     output: "MultiSampleVCF/Combined_multisampleGenotype_recalSNPs_recalINDELs_passed.vcf.gz"
     params: dbsnp150=dbsnp_150
     log:"LOGFILES/removeFilteredVariants.log"
     shell:"""
           module load snpEff/4.3
           (zcat {input} | java -Xmx10g -jar $snpEff/SnpSift.jar filter "( na FILTER ) | (FILTER = 'PASS')" \
           | awk -v OFS='\t' '/^chr/{{$3="."}} {{print}}' | java -Xmx10g -jar $snpEff/SnpSift.jar annotate -id -noInfo {params.dbsnp150} - ) > {output} 2> {log}
           rm -rf {input}
           """
