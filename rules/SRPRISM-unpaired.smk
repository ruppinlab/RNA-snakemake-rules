from os.path import join
import pandas as pd

wildcard_constraints:
    genome="SL1344"

SL1344_genome_URL = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/210/855/GCF_000210855.2_ASM21085v2/GCF_000210855.2_ASM21085v2_genomic.fna.gz"

# files
GENOME_FA = join("raw", "{genome}.fa")
SRPRISM_DB_DIR = join("output", "genomes_db")
GENOME_DB_PREFIX = join(SRPRISM_DB_DIR, "{genome}")
GENOME_DB_FILE = join(SRPRISM_DB_DIR, "{genome}.idx")

#FASTQ_PREFIX = "/data/Robinson-SB/scRNA-seq-microbe-identification/Chung2017/FASTQ/"

#Input_FQ1 = join(FASTQ_PREFIX, "raw", "{cell}_1.fastq.gz")
#Input_FQ2 = join(FASTQ_PREFIX, "raw", "{cell}_2.fastq.gz")

# sam files
SRPRISM_UNPAIRED_SAM = join("output", "SRPRISM", "{patient}", "{sample}-{cell}", "{genome}-unpaired.sam")
SRPRISM_UNPAIRED_PRIMARY_SAM = join("output", "SRPRISM", "{patient}", "{sample}-{cell}", "{genome}-paired.primary.sam")
SRPRISM_PROPER_PAIRED_PRIMARY_SAM = join("output", "SRPRISM", "{patient}", "{sample}-{cell}", "{genome}-proper-paired.primary.sam")
SRPRISM_PROPER_PAIRED_PRIMARY_BAM = join("output", "SRPRISM", "{patient}", "{sample}-{cell}", "{genome}-proper-paired.primary.bam")
SRPRISM_PROPER_PAIRED_PRIMARY_SORTED_BAM = join("output", "SRPRISM", "{patient}", "{sample}-{cell}", "{genome}-paired.primary.sorted.bam")
SRPRISM_PROPER_PAIRED_PRIMARY_SORTED_BAI = join("output", "SRPRISM", "{patient}", "{sample}-{cell}", "{genome}-paired.primary.sorted.bam.bai")

SRPRISM_COUNT = join("output", "{genome}-read_count.tsv")

localrules: extract_primary_alignment, convert_to_bam, sort_bam, index_bam

#
# rule all:
#     input:
#         expand(SRPRISM_COUNT, genome="nucleatum")

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
        SRPRISM_UNPAIRED_PRIMARY_SAM
    output:
        SRPRISM_UNPAIRED_PRIMARY_BAM
    shell:
        "module load samtools && "
        "samtools view -h -b -o {output} {input}"


# -F means exclude and 256 means "not primary alignment"
rule extract_primary_alignment:
    input:
        SRPRISM_UNPAIRED_SAM
    output:
        SRPRISM_UNPAIRED_PRIMARY_SAM
    shell:
        "module load samtools && "
        "samtools view -h -F 256 -o {output} {input}"

rule map_SRPRISM_genome_unpaired:
    params:
        GENOME_DB_PREFIX
    input:
        SRPRISM_INPUT_FQ1,
        GENOME_DB_FILE
    output:
        SRPRISM_UNPAIRED_SAM
    shell:
        "srprism search -I {params} -i {input[0]} -F fastq -p false -o {output} --sam-header true"

rule make_reference_DB:
    params:
        GENOME_DB_PREFIX
    input:
        GENOME_FA
    output:
        GENOME_DB_FILE
    shell:
        "srprism mkindex -i {input} -o {params}"

rule download_relevant_genomes:
    wildcard_constraints:
        genome="SL1344"
    params:
        SL1344_genome_URL
    output:
        GENOME_FA
    shell:
        "wget -O - {params[0]} | gunzip -c > {output}"
