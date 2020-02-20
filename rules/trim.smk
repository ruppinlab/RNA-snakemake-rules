from os.path import join

# directories
ENV_DIR = join("..", "envs")
FASTQ_DIR = "FASTQ"
DATA_DIR = "data"
TRIMMED_DIR = "trimmed"
# input files
FASTQ1_FILE = join(FASTQ_DIR, "raw", "{patient}-{sample}_1.fastq.gz")
FASTQ2_FILE = join(FASTQ_DIR, "raw", "{patient}-{sample}_2.fastq.gz")

# output files
FASTP_JSON_REPORT = join(FASTQ_DIR, "trimmed", "{patient}-{sample}-report.json")
FASTP_HTML_REPORT = join(FASTQ_DIR, "trimmed", "{patient}-{sample}-report.html")
TRIMMED_FASTQ1_FILE = join(FASTQ_DIR, "trimmed", "{patient}-{sample}_1.fastq.gz")
TRIMMED_FASTQ2_FILE = join(FASTQ_DIR, "trimmed", "{patient}-{sample}_2.fastq.gz")


# fastp
rule run_fastp:
    conda:
        join(ENV_DIR, "fastp.yml")
    input:
        FASTQ1_FILE,
        FASTQ2_FILE
    output:
        temp(TRIMMED_FASTQ1_FILE),
        temp(TRIMMED_FASTQ2_FILE),
        FASTP_JSON_REPORT,
        FASTP_HTML_REPORT
    threads:
        6
    shell:
        # -l 50 requires reads be of length 50 or greater
        # quality filtering is enabled (default); default phred quality >= 15
        # to be qualified; at least 40% of bases must be qualified adapter
        # trimming happens for paired end reads using a pre-read overlap
        # analysis per Jared, we want to cut tail  with window length 4 and
        # phred score >= 25
        "fastp -w {threads} --cut_tail --cut_tail_window_size 4 "
        "--cut_tail_mean_quality 25 -l 50 --in1 '{input[0]}' "
        "--in2 '{input[1]}' --out1 '{output[0]}' --out2 '{output[1]}' "
        "-j '{output[2]}' -h '{output[3]}'"
