from os.path import join

# directories
ENV_DIR = join("..", "envs")
FASTQ_DIR = "FASTQ"
DATA_DIR = "data"
TRIMMED_DIR = "trimmed"
# input files
FASTQ1_FILE = join(FASTQ_DIR, "raw", "{patient}-{sample}-{cell}_1.fastq.gz")
FASTQ2_FILE = join(FASTQ_DIR, "raw", "{patient}-{sample}-{cell}_2.fastq.gz")

# output files
FASTP_JSON_REPORT = join(FASTQ_DIR, "trimmed", "{patient}-{sample}-{cell}-report.json")
FASTP_HTML_REPORT = join(FASTQ_DIR, "trimmed", "{patient}-{sample}-{cell}-report.html")
TRIMMED_FASTQ1_FILE = join(FASTQ_DIR, "trimmed", "{patient}-{sample}-{cell}_1.fastq.gz")
TRIMMED_FASTQ2_FILE = join(FASTQ_DIR, "trimmed", "{patient}-{sample}-{cell}_2.fastq.gz")
TRIMMED_UNPAIRED_FILE = join(FASTQ_DIR, "trimmed", "{patient}-{sample}-{cell}_3.fastq.gz")
FAILED_READS_FILE = join(FASTQ_DIR, "trimmed", "{patient}-{sample}-{cell}_failed.fastq.gz")

# fastp
rule run_fastp:
    group:
        "run_fastp"
    conda:
        join(ENV_DIR, "fastp.yml")
    input:
        FASTQ1_FILE,
        FASTQ2_FILE
    output:
        TRIMMED_FASTQ1_FILE,
        TRIMMED_FASTQ2_FILE,
        TRIMMED_UNPAIRED_FILE,
        FAILED_READS_FILE,
        FASTP_JSON_REPORT,
        FASTP_HTML_REPORT
    threads:
        6
    benchmark:
        "benchmarks/{patient}-{sample}-{cell}.run_fastp.benchmark.txt"
    shell:
        "fastp -w {threads} "
        "--unqualified_percent_limit 40 " # filter reads where 40% of bases have phred quality < 15
        "--cut_tail " # use defaults --cut_window_size 4 --cut_mean_quality 20
        "--low_complexity_filter " # filter reads with less than 30% complexity (30% of the bases are different from the preceeding base)
        "--trim_poly_x " # trim poly X's - useful for trimming polyA tails
        "-i {input[0]} -I {input[1]} -o {output[0]} -O {output[1]} "
        "--unpaired1 {output[2]} --unpaired2 {output[2]} --failed_out {output[3]} "
        "-j {output[4]} -h {output[5]} "
        + config["params"]["fastp"]
