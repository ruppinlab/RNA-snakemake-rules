from os.path import join, exists
import json


# Directories
ENV_DIR = join("..", "envs")
STAR_OUTPUT_DIR = join("output", "star", "{patient}-{sample}-{plate}-{cell}")
STAR_PE_OUTPUT_DIR = join(STAR_OUTPUT_DIR, "_STARpe")
STAR_SE_OUTPUT_DIR = join(STAR_OUTPUT_DIR, "_STARse")
# Files
STAR_ENV_FILE = join(ENV_DIR, "star.yml")

# STAR output directories
STAR_GENOME_INDEX = join("output", "star-index")

SE_FASTQ_FILE = join("FASTQ", "trimmed", "{patient}-{sample}-{cell}_3.fastq.gz")

STAR_PE_BAM_FILE = join(STAR_PE_OUTPUT_DIR, "Aligned.out.bam")
STAR_PE_READCOUNT_FILE = join(STAR_PE_OUTPUT_DIR, "ReadsPerGene.out.tab")
STAR_SE_BAM_FILE = join(STAR_SE_OUTPUT_DIR, "Aligned.out.bam")
STAR_SE_READCOUNT_FILE = join(STAR_SE_OUTPUT_DIR, "ReadsPerGene.out.tab")

def get_fq(wildcards):
    if config["trimming"]["skip"]:
        # no trimming, use raw reads
        dir = join("FASTQ", "raw")
    else:
        dir = join("FASTQ", "trimmed")
    return {
        'fq1': join(dir, "{wildcards.patient}-{wildcards.sample}-{wildcards.cell}_1.fastq.gz".format(wildcards=wildcards)),
        'fq2': join(dir, "{wildcards.patient}-{wildcards.sample}-{wildcards.cell}_2.fastq.gz".format(wildcards=wildcards))
        }


rule create_star_index:
    conda:
        STAR_ENV_FILE
    input:
        genome = config["ref"]["genome"],
        gtf = config["ref"]["annotation"],
    output:
        directory(STAR_GENOME_INDEX)
    threads:
        16
    shell:
        "mkdir '{output}' && STAR "
        "--runThreadN {threads} "
        "--runMode genomeGenerate "
        "--genomeDir '{output}' "
        "--genomeFastaFiles '{input[genome]}' "
        "--sjdbGTFfile '{input[gtf]}'"


rule STAR_PE:
    group:
        "FASTQ"
    conda:
        STAR_ENV_FILE
    input:
        STAR_GENOME_INDEX,
        unpack(get_fq)
    params:
        odir = join(STAR_PE_OUTPUT_DIR, "")
    output:
        STAR_PE_BAM_FILE,
        STAR_PE_READCOUNT_FILE
    threads:
        16
    benchmark:
        "benchmarks/{patient}-{sample}-{plate}-{cell}.STAR_PE.benchmark.txt"
    shell:
        "STAR "
        "--runThreadN {threads} "
        "--readFilesIn '{input.fq1}' '{input.fq2}' "
        "--genomeDir '{input[0]}' "
        "--outSAMunmapped Within "
        "--readFilesCommand zcat "
        "--outSAMtype BAM Unsorted "
        "--outSAMattrRGline ID:{wildcards.cell} "
        "PL:illumina SM:{wildcards.cell} LB:RNA "
        "--outFileNamePrefix '{params.odir}' "
        "--quantMode GeneCounts "

rule STAR_SE:
    group:
        "FASTQ"
    conda:
        STAR_ENV_FILE
    input:
        STAR_GENOME_INDEX,
        SE_FASTQ_FILE
    params:
        odir = join(STAR_SE_OUTPUT_DIR, "")
    output:
        STAR_SE_BAM_FILE,
        STAR_SE_READCOUNT_FILE
    threads:
        16
    benchmark:
        "benchmarks/{patient}-{sample}-{plate}-{cell}.STAR_SE.benchmark.txt"
    shell:
        "STAR "
        "--runThreadN {threads} "
        "--readFilesIn '{input[1]}' "
        "--genomeDir '{input[0]}' "
        "--outSAMunmapped Within "
        "--readFilesCommand zcat "
        "--outSAMtype BAM Unsorted "
        "--outSAMattrRGline ID:{wildcards.cell} "
        "PL:illumina SM:{wildcards.cell} LB:RNA "
        "--outFileNamePrefix '{params.odir}' "
        "--quantMode GeneCounts "
