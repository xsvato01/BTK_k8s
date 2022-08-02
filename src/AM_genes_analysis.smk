####################################
## Andrea Mareckova gene pipeline ##
####################################

# running under lymphopanel conda environment

import re
import datetime

## REFERENCES AND BEDS
REF="/mnt/ssd/ssd_3/references/homsap/GRCh37-p13/index/BWA/GRCh37-p13"
# for coverage
GFF= "/mnt/nfs/shared/MedGen/BTK/beds/BTK_PLCG2_corr.bed"
# for vardict
BED= "/mnt/nfs/shared/MedGen/BTK/beds/BTK_PLCG2_corr_toVardict.bed"

# Tools and scripts
SAMTOOLS= "samtools"
BEDTOOLS= "/usr/bin/bedtools"
BWA = "/mnt/ssd/ssd_2/install/dir/anaconda/envs/lymphopanel/bin/bwa"
CUTADAPT = "/mnt/ssd/ssd_2/install/dir/anaconda/envs/lymphopanel/bin/cutadapt"
VARSCAN = "varscan"
VARDICT = "vardict"
BCFTOOLS = "/mnt/ssd/ssd_2/install/dir/anaconda/envs/lymphopanel/bin/bcftools"
GATK = "/mnt/nfs/shared/999993-Bioda/projects/TP53/resources/GATK_3.6/GenomeAnalysisTK.jar"
VEP = "variant_effect_predictor.pl"
BIOPET = "biopet"

VCF_SIMPLIFY= "/mnt/nfs/shared/MedGen/BTK/src/vcf-simplify.py"
CREATE_TABLE= "/mnt/nfs/shared/MedGen/BTK/src/rearange_table_BTK_AMareckova.R"
COV_STAT= "/mnt/nfs/shared/MedGen/BTK/src/coverage_final_table.R"

VEP_DATA = "/mnt/ssd/ssd_3/references/homsap/GRCh37-p13_r87/annot/vep/"

def getFasta(wildcards):
    files = config["Samples"][wildcards.sample]
    return files

def covlist(wildcards):
    covlist = []
    samples = config["Samples"]
    run = config["Run"]
    for sample in samples:
        covlist = covlist + ["../data/" + run + "/coverage/" + sample + ".perexon_stat.txt"]
    return covlist

def variantlist(wildcards):
     variantlist = []
     samples = config["Samples"]
     run = config["Run"]
     for sample in samples:
         variantlist = variantlist + ["../data/" + run + "/annot/" + sample + ".sample.merged.anot.txt"]
     return variantlist

def bamlist(wildcards):
    bamlist = []
    samples = config["Samples"]
    run = config["Run"]
    for sample in samples:
        bamlist = bamlist + ["../data/" + run + "/bam/" + sample + ".sorted.bam"]
    return bamlist


def all_variants(wildcards):
    all_variants = []
    run = config["Run"]
    sample = config["Samples"[0]]
    all_variants = all_variants + ["../data/" + run + "/annot/" + sample + ".sample.merged.anot.txt"]
    return all_variants

rule all:
    output:
        "../run_" + str(datetime.datetime.now()).replace(" ","_")
    input:
        bamlist,
        covlist,
        variantlist
    shell:
        "echo {input}; touch {output}"

rule coverage_stat:
    output:
        expand("../data/{run}/coverage/{sample}.perexon_stat.txt", run = config["Run"], sample = config["Samples"])
    input:
        expand("../data/{run}/coverage/{sample}.PBcoverage.txt", run = config["Run"], sample = config["Samples"])
    run:
        run_ID = config["Run"]
        command = "Rscript --vanilla " + COV_STAT \
            + " ../data/" + run_ID + "/coverage " + run_ID
        shell(command)

rule create_final_table:
    output:
        expand("../data/{run}/annot/{sample}.sample.merged.anot.txt", run = config["Run"], sample = config["Samples"])
    input:
        expand("../data/{run}/annot/{sample}.norm.merged.annot.normVEP.txt", run = config["Run"],sample = config["Samples"])
    run:
        run_ID = config["Run"]
        command = "Rscript --vanilla " + CREATE_TABLE \
            + " ../data/" + run_ID + "/annot " + run_ID
        shell(command)

rule create_txt:
    output:
        "../data/{run}/annot/{sample}.norm.merged.annot.normVEP.txt"
    input:
        "../data/{run}/annot/{sample}.norm.merged.annot.normVEP.vcf"
    run:
        command = "python3.6 " + VCF_SIMPLIFY + " SimplifyVCF -toType table" \
            + " -inVCF " + str(input) \
            + " -out " + str(output)
        shell(command)

rule vep_normalizer:
    output:
        "../data/{run}/annot/{sample}.norm.merged.annot.normVEP.vcf"
    input:
        "../data/{run}/annot/{sample}.norm.merged.annot.vcf"
    run:
        command = BIOPET + " tool VepNormalizer" \
            + " -I " + str(input) \
            + " -O " + str(output) + " -m explode"
        shell(command)

rule annotate:
    output:
        "../data/{run}/annot/{sample}.norm.merged.annot.vcf"
    input:
        "../data/{run}/vcf/merged/{sample}.allcallers.merged.norm.vcf"

    run:
        command = VEP + " -i " + str(input) \
            + " --cache --cache_version 87 --dir_cache " + VEP_DATA \
            + " --fasta " + "/mnt/ssd/ssd_3/references/homsap/GRCh37-p13/seq/GRCh37-p13.fa" \
            + " --merged --offline --vcf --hgvs --sift b --polyphen b" \
            + " -o " + str(output)
        shell(command)

rule normalize_merged_variants:
    output:
        "../data/{run}/vcf/merged/{sample}.allcallers.merged.norm.vcf"
    input:
        "../data/{run}/vcf/merged/{sample}.allcallers.merged.vcf"
    run:
        command = BCFTOOLS + " norm -f " + "/mnt/ssd/ssd_3/references/homsap/GRCh37-p13/seq/GRCh37-p13.fa" + " " + str(input) \
            + " -o " + str(output)
        shell(command)

rule merge:
    output:
        "../data/{run}/vcf/merged/{sample}.allcallers.merged.vcf"
    input:
        varscan_snv = "../data/{run}/vcf/varscan2/{sample}.varscan.snv.norm.vcf",
        varscan_indel = "../data/{run}/vcf/varscan2/{sample}.varscan.indel.norm.vcf",
        vardict = "../data/{run}/vcf/vardict/{sample}.vardict.norm.vcf"
    run:
        command = "java -jar " + GATK + " -T CombineVariants" \
            + " --variant:vardict " + str(input.vardict) \
            + " --variant:varscan_SNV " + str(input.varscan_snv) \
            + " --variant:varscan_INDEL " + str(input.varscan_indel) \
            + " -R " + "/mnt/ssd/ssd_3/references/homsap/GRCh37-p13/seq/GRCh37-p13.fa" + " -genotypeMergeOptions UNSORTED" \
            + " --disable_auto_index_creation_and_locking_when_reading_rods" \
            + " -o " + str(output)
        shell(command)

rule normalize_variants:
    output:
        varscan_snv_norm = "../data/{run}/vcf/varscan2/{sample}.varscan.snv.norm.vcf",
        varscan_indel_norm = "../data/{run}/vcf/varscan2/{sample}.varscan.indel.norm.vcf",
        vardict_norm = "../data/{run}/vcf/vardict/{sample}.vardict.norm.vcf"
    input:
        varscan_snv = "../data/{run}/vcf/varscan2/{sample}.varscan.snv.vcf",
        varscan_indel = "../data/{run}/vcf/varscan2/{sample}.varscan.indel.vcf",
        vardict = "../data/{run}/vcf/vardict/{sample}.vardict.vcf"
    run:
        command1 = BCFTOOLS + " norm -f " + "/mnt/ssd/ssd_3/references/homsap/GRCh37-p13/seq/GRCh37-p13.fa" + " " + str(input.varscan_snv) \
            + " -o " + str(output.varscan_snv_norm)
        command2 = BCFTOOLS + " norm -f " + "/mnt/ssd/ssd_3/references/homsap/GRCh37-p13/seq/GRCh37-p13.fa" + " " + str(input.varscan_indel) \
            + " -o " + str(output.varscan_indel_norm)
        command3 = BCFTOOLS + " norm -f " + "/mnt/ssd/ssd_3/references/homsap/GRCh37-p13/seq/GRCh37-p13.fa" + " " +  str(input.vardict) \
            + " -o " + str(output.vardict_norm) 
#           + " | " + BCFTOOLS + " sort - -o " + str(output.vardict_norm)
        shell(command1)
        shell(command2)
        shell(command3)

rule vardict:
    output:
        "../data/{run}/vcf/vardict/{sample}.vardict.vcf"
    input:
        "../data/{run}/bam/{sample}.sorted.bam"
    run:
        command = VARDICT + " -G " + "/mnt/ssd/ssd_3/references/homsap/GRCh37-p13/seq/GRCh37-p13.fa" + " -f 0.0005 -b" \
            + " " + str(input) \
            + " -c 1 -S 2 -E 3 -r 8 -Q 1 -q 25 -P 2 -m 8" \
            + " " + BED + " | teststrandbias.R | var2vcf_valid.pl -f 0.0005 -d 50 -c 5 -p 2 -q 25 -Q 1 -v 8 -m 8 -N vardict -" \
            + " > " + str(output)
        shell(command)

rule varscan:
    output:
        snv = "../data/{run}/vcf/varscan2/{sample}.varscan.snv.vcf",
        indel = "../data/{run}/vcf/varscan2/{sample}.varscan.indel.vcf"
    input:
        "../data/{run}/vcf/varscan2/{sample}.mpileup"
    run:
        command_snv = VARSCAN + " mpileup2snp " + str(input) \
            + " --strand-filter 0 --p-value 0.95 --min-coverage 50 --min-reads2 8 --min-avg-qual 25 --min-var-freq 0.0005" \
            + " --output-vcf > " + str(output.snv)
        command_indel = VARSCAN + " mpileup2indel " + str(input) \
            + " --strand-filter 0 --p-value 0.95 --min-coverage 50 --min-reads2 8 --min-avg-qual 25 --min-var-freq 0.0005" \
            + " --output-vcf > " + str(output.indel)
        shell(command_snv)
        shell(command_indel)
        shell("sed -i 's/Sample1/varscan_SNV/g' {output.snv}")
        shell("sed -i 's/Sample1/varscan_INDEL/g' {output.indel}")

rule mpileup:
    output:
        "../data/{run}/vcf/varscan2/{sample}.mpileup"
    input:
        "../data/{run}/bam/{sample}.sorted.bam"
    run:
        command = SAMTOOLS + " mpileup -x -B -Q 25 -d 999999 -L 999999 -F 0.0005" \
            + " -f " + "/mnt/ssd/ssd_3/references/homsap/GRCh37-p13/seq/GRCh37-p13.fa" + " " + str(input) \
            + " > " + str(output)
        shell(command)

rule coverage:
    output:
        "../data/{run}/coverage/{sample}.PBcoverage.txt"
    input:
        "../data/{run}/bam/{sample}.sorted.bam"
    run:
        command = BEDTOOLS + " coverage -abam " + GFF \
            + " -b " + str(input) + " -d > " + str(output)
        shell(command)

rule align:
    output:
        "../data/{run}/bam/{sample}.sorted.bam"
    input:
        fwd = "../data/{run}/fastq_trimmed/{sample}.R1.trimmed.fastq.gz",
        rev = "../data/{run}/fastq_trimmed/{sample}.R2.trimmed.fastq.gz"
    log:
        run = "../data/{run}/bam/{sample}.log"
    run:
        id = wildcards.sample
        sm = wildcards.sample
        command = BWA + " mem -t 30 -R '@RG\\tID:" + id + "\\tSM:" + sm + "\\tPL:illumina' -v 1 " \
            + REF + " " + str(input.fwd) + " " + str(input.rev) + " 2>>" \
            + " " + log.run + " | " + SAMTOOLS + " view -bS - | " + SAMTOOLS + " sort -T - -o" \
            + " " + str(output)
        shell(command)
        shell("samtools index {output}")

rule trim2:
    output:
        fwd = "../data/{run}/fastq_trimmed/{sample}.R1.trimmed.fastq.gz",
        rev = "../data/{run}/fastq_trimmed/{sample}.R2.trimmed.fastq.gz"
    input:
        fwd = "../data/{run}/fastq_trimmed/{sample}.R1.tmp.fastq.gz",
        rev = "../data/{run}/fastq_trimmed/{sample}.R2.tmp.fastq.gz"
    run:
        command = CUTADAPT + " -a CAAGGGGGACTGTAGATGGG...TAGGATCTGACTGCGGCTCC" \
            + " -A GGAGCCGCAGTCAGATCCTA...CCCATCTACAGTCCCCCTTG" \
            + " -a ACAACGTTCTGGTAAGGACAX -A TGTCCTTACCAGAACGTTGTX --overlap 4" \
            + " -o " + str(output.fwd) + " -p " + str(output.rev) \
            + " " + str(input.fwd) + " " + str(input.rev)
        shell(command)

rule trim1:
    output:
        fwd = temp("../data/{run}/fastq_trimmed/{sample}.R1.tmp.fastq.gz"),
        rev = temp("../data/{run}/fastq_trimmed/{sample}.R2.tmp.fastq.gz")
    input:
        unpack(getFasta)
    run:
        command = CUTADAPT + " -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT" \
            + " -o " + str(output.fwd) + " -p " + str(output.rev) \
            + " " + str(input.fwd) + " " + str(input.rev)
        shell(command)
