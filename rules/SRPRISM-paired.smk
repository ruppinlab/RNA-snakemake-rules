from os.path import join
import pandas as pd

# files
GENOME_FA = join("raw", "{genome}.fa")
SRPRISM_DB_DIR = join("output", "genomes_db")
GENOME_DB_PREFIX = join(SRPRISM_DB_DIR, "{genome}")
GENOME_DB_FILE = join(SRPRISM_DB_DIR, "{genome}.idx")

# sam files
SRPRISM_PAIRED_SAM = join("output", "SRPRISM", "{patient}", "{sample}-{plate}-{cell}", "{genome}-paired.sam")
SRPRISM_PAIRED_PRIMARY_SAM = join("output", "SRPRISM", "{patient}", "{sample}-{plate}-{cell}", "{genome}-paired.primary.sam")
SRPRISM_PROPER_PAIRED_PRIMARY_SAM = join("output", "SRPRISM", "{patient}", "{sample}-{plate}-{cell}", "{genome}-proper-paired.primary.sam")
SRPRISM_PROPER_PAIRED_PRIMARY_BAM = join("output", "SRPRISM", "{patient}", "{sample}-{plate}-{cell}", "{genome}-proper-paired.primary.bam")
SRPRISM_PROPER_PAIRED_PRIMARY_SORTED_BAM = join("output", "SRPRISM", "{patient}", "{sample}-{plate}-{cell}", "{genome}-paired.primary.sorted.bam")
SRPRISM_PROPER_PAIRED_PRIMARY_SORTED_BAI = join("output", "SRPRISM", "{patient}", "{sample}-{plate}-{cell}", "{genome}-paired.primary.sorted.bam.bai")

# SRPRISM_COUNT = join("output", "{genome}-read_count.tsv")

localrules: extract_primary_alignment, convert_to_bam, sort_bam, index_bam, exclude_non_proper_pairs


# rule all:
#     input:
#         expand(SRPRISM_COUNT, genome="nucleatum")
#
# rule count_nreads:
#     conda:
#         "../envs/pysam-env.yaml"
#     params:
#         cells=df["cell"]
#     input:
#         expand(SRPRISM_PROPER_PAIRED_PRIMARY_SORTED_BAM, genome="{genome}", cell=df["cell"])
#     output:
#         SRPRISM_COUNT
#     script:
#         "src/count_nreads.py"
#
# rule index_bam:
#     input:
#         SRPRISM_PROPER_PAIRED_PRIMARY_SORTED_BAM
#     output:
#         SRPRISM_PROPER_PAIRED_PRIMARY_SORTED_BAI
#     shell:
#         "module load samtools && "
#         "samtools index {input}"
#
# # sort bam by genomic coordinates so we can easily filter by region(s)
# rule sort_bam:
#     input:
#         SRPRISM_PROPER_PAIRED_PRIMARY_BAM
#     output:
#         SRPRISM_PROPER_PAIRED_PRIMARY_SORTED_BAM
#     shell:
#         "module load samtools && "
#         "samtools sort -o {output} {input}"

rule convert_to_bam:
    input:
        SRPRISM_PROPER_PAIRED_PRIMARY_SAM
    output:
        SRPRISM_PROPER_PAIRED_PRIMARY_BAM
    shell:
        "module load samtools && "
        "samtools view -h -b -o {output} {input}"

# -f means include and 2 means "read mapped in proper pair"
rule exclude_non_proper_pairs:
    input:
        SRPRISM_PAIRED_PRIMARY_SAM
    output:
        SRPRISM_PROPER_PAIRED_PRIMARY_SAM
    shell:
        "module load samtools && "
        "samtools view -h -f 2 -o {output} {input}"

# -F means exclude and 256 means "not primary alignment"
rule extract_primary_alignment:
    input:
        SRPRISM_PAIRED_SAM
    output:
        SRPRISM_PAIRED_PRIMARY_SAM
    shell:
        "module load samtools && "
        "samtools view -h -F 256 -o {output} {input}"

rule map_SRPRISM_genome_paired:
    params:
        GENOME_DB_PREFIX
    input:
        SRPRISM_INPUT_FQ1,
        SRPRISM_INPUT_FQ2,
        GENOME_DB_FILE
    output:
        SRPRISM_PAIRED_SAM
    shell:
        "srprism search -I {params} -i {input[0]},{input[1]} -F fastq -p true -o {output} --sam-header true"

rule make_reference_DB:
    params:
        GENOME_DB_PREFIX
    input:
        GENOME_FA
    output:
        GENOME_DB_FILE
    shell:
        "srprism mkindex -i {input} -o {params}"
