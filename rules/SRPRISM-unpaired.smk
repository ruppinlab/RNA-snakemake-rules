from os.path import join
import pandas as pd


include: "SRPRISM.smk"

# input files
SRPRISM_UNPAIRED_INPUT_FQ = join("output", "SRPRISM", "{patient}", "{identifier}", "unaligned_3.fq")

# sam files
SRPRISM_UNPAIRED_SAM = join("output", "SRPRISM", "{patient}", "{identifier}", "{genome}-unpaired.sam")
SRPRISM_UNPAIRED_PRIMARY_SAM = join("output", "SRPRISM", "{patient}", "{identifier}", "{genome}-unpaired.primary.sam")
SRPRISM_UNPAIRED_PRIMARY_BAM = join("output", "SRPRISM", "{patient}", "{identifier}", "{genome}-unpaired.primary.bam")
SRPRISM_UNPAIRED_PRIMARY_SORTED_BAM = join("output", "SRPRISM", "{patient}", "{identifier}", "{genome}-unpaired.primary.sorted.bam")
SRPRISM_UNPAIRED_PRIMARY_SORTED_BAI = join("output", "SRPRISM", "{patient}", "{identifier}", "{genome}-unpaired.primary.sorted.bam.bai")


localrules: extract_primary_alignment_from_unpaired_bam, convert_unpaired_sam_to_bam, sort_unpaired_bam, index_unpaired_bam


rule index_unpaired_bam:
    group:
        "SRPRISM"
    input:
        SRPRISM_UNPAIRED_PRIMARY_SORTED_BAM
    output:
        SRPRISM_UNPAIRED_PRIMARY_SORTED_BAI
    shell:
        "module load samtools && "
        "samtools index {input}"

# sort bam by genomic coordinates so we can easily filter by region(s)
rule sort_unpaired_bam:
    group:
        "SRPRISM"
    input:
        SRPRISM_UNPAIRED_PRIMARY_BAM
    output:
        SRPRISM_UNPAIRED_PRIMARY_SORTED_BAM
    shell:
        "module load samtools && "
        "samtools sort -o {output} {input}"

rule convert_unpaired_sam_to_bam:
    group:
        "SRPRISM"
    input:
        SRPRISM_UNPAIRED_PRIMARY_SAM
    output:
        SRPRISM_UNPAIRED_PRIMARY_BAM
    shell:
        "module load samtools && "
        "samtools view -h -b -o {output} {input}"


# -F means exclude and 256 means "not primary alignment"
rule extract_primary_alignment_from_unpaired_bam:
    group:
        "SRPRISM"
    input:
        SRPRISM_UNPAIRED_SAM
    output:
        SRPRISM_UNPAIRED_PRIMARY_SAM
    shell:
        "module load samtools && "
        "samtools view -h -F 256 -o {output} {input}"

rule map_SRPRISM_genome_unpaired:
    group:
        "SRPRISM"
    params:
        GENOME_DB_PREFIX
    input:
        SRPRISM_UNPAIRED_INPUT_FQ,
        GENOME_DB_FILE
    output:
        SRPRISM_UNPAIRED_SAM
    shell:
        "srprism search -I {params} -i {input[0]} -F fastq -p false -o {output} --sam-header true"
