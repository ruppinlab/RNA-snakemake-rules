from os.path import join
import pandas as pd

include: "SRPRISM.smk"

# input files
SRPRISM_INPUT_FQ1 = join("output", "SRPRISM", "{patient}", "{identifier}", "unaligned_1.fq")
SRPRISM_INPUT_FQ2 = join("output", "SRPRISM", "{patient}", "{identifier}", "unaligned_2.fq")


# sam files
SRPRISM_PAIRED_SAM = join("output", "SRPRISM", "{patient}", "{identifier}", "{genome}-paired.sam")
SRPRISM_PAIRED_PRIMARY_SAM = join("output", "SRPRISM", "{patient}", "{identifier}", "{genome}-paired.primary.sam")
SRPRISM_PROPER_PAIRED_PRIMARY_SAM = join("output", "SRPRISM", "{patient}", "{identifier}", "{genome}-proper-paired.primary.sam")
SRPRISM_PROPER_PAIRED_PRIMARY_BAM = join("output", "SRPRISM", "{patient}", "{identifier}", "{genome}-proper-paired.primary.bam")
SRPRISM_PROPER_PAIRED_PRIMARY_SORTED_BAM = join("output", "SRPRISM", "{patient}", "{identifier}", "{genome}-paired.primary.sorted.bam")
SRPRISM_PROPER_PAIRED_PRIMARY_SORTED_BAI = join("output", "SRPRISM", "{patient}", "{identifier}", "{genome}-paired.primary.sorted.bam.bai")


localrules: extract_primary_alignment, convert_to_bam, sort_bam, index_bam, exclude_non_proper_pairs


# index bam by genomic coordinates so we can easily filter by region(s)
rule index_bam:
    group:
        "index_bam"
    input:
        SRPRISM_PROPER_PAIRED_PRIMARY_SORTED_BAM
    output:
        SRPRISM_PROPER_PAIRED_PRIMARY_SORTED_BAI
    shell:
        "module load samtools && "
        "samtools index {input}"

# sort bam by genomic coordinates so we can index
rule sort_bam:
    group:
        "sort_bam"
    input:
        SRPRISM_PROPER_PAIRED_PRIMARY_BAM
    output:
        SRPRISM_PROPER_PAIRED_PRIMARY_SORTED_BAM
    shell:
        "module load samtools && "
        "samtools sort -o {output} {input}"

rule convert_to_bam:
    group:
        "convert_to_bam"
    input:
        SRPRISM_PROPER_PAIRED_PRIMARY_SAM
    output:
        temp(SRPRISM_PROPER_PAIRED_PRIMARY_BAM)
    shell:
        "module load samtools && "
        "samtools view -h -b -o {output} {input}"

# -f means include and 2 means "read mapped in proper pair"
rule exclude_non_proper_pairs:
    group:
        "exclude_non_proper_pairs"
    input:
        SRPRISM_PAIRED_PRIMARY_SAM
    output:
        temp(SRPRISM_PROPER_PAIRED_PRIMARY_SAM)
    shell:
        "module load samtools && "
        "samtools view -h -f 2 -o {output} {input}"

# -F means exclude and 256 means "not primary alignment"
rule extract_primary_alignment:
    group:
        "extract_primary_alignment"
    input:
        SRPRISM_PAIRED_SAM
    output:
        temp(SRPRISM_PAIRED_PRIMARY_SAM)
    shell:
        "module load samtools && "
        "samtools view -h -F 256 -o {output} {input}"

rule map_SRPRISM_genome_paired:
    group:
        "map_SRPRISM_genome_paired"
    params:
        GENOME_DB_PREFIX
    input:
        SRPRISM_INPUT_FQ1,
        SRPRISM_INPUT_FQ2,
        GENOME_DB_FILE
    output:
        temp(SRPRISM_PAIRED_SAM)
    shell:
        "srprism search -I {params} -i {input[0]},{input[1]} -F fastq -p true -o {output} --sam-header true"
