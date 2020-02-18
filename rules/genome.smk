from os.path import join

# URLs for downloading reference sequence and annotation files
GENCODE_HUMAN_BASE_URL = (
    "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32")
GENCODE_TRANSCRIPT_FASTA_URL = join(GENCODE_HUMAN_BASE_URL,
                                    "gencode.v32.transcripts.fa.gz")
GENCODE_GENOME_FASTA_URL = join(GENCODE_HUMAN_BASE_URL,
                                "GRCh38.p13.genome.fa.gz")
GENCODE_GTF_URL = join(GENCODE_HUMAN_BASE_URL,
                       "gencode.v32.annotation.gtf.gz")
GDC_DATA_BASE_URL = "https://api.gdc.cancer.gov/data"
GDC_GENOME_FASTA_URL = join(GDC_DATA_BASE_URL,
                            "254f697d-310d-4d7d-a27b-27fbf767a834")
GDC_GTF_URL = join(GDC_DATA_BASE_URL, "25aa497c-e615-4cb7-8751-71f744f9691f")

# files
GDC_GENOME_FASTA_TARGZ_FILE = join("genome", "raw", "GRCh38.d1.vd1.fa.tar.gz")
GDC_GENOME_FASTA_FILE = join("genome", "raw", "GRCh38.d1.vd1.fa")
GDC_GTF_FILE = join("genome", "raw", "gencode.v22.annotation.gtf")
GENCODE_TRANSCRIPT_FASTA_FILE = join("genome", "raw", "gencode.v32.transcripts.fa.gz")
GENCODE_GENOME_FASTA_FILE = join("genome", "raw", "GRCh38.p13.genome.fa")
GENCODE_GTF_FILE = join("genome", "raw", "gencode.v32.annotation.gtf")


# rules for downloading reference sequence and annoation files
rule uncompress_gdc_genome_fasta:
    input:
        GDC_GENOME_FASTA_TARGZ_FILE
    output:
        GDC_GENOME_FASTA_FILE
    shell:
        "tar -zxvf {input}"

rule download_gdc_genome_fasta:
    params:
        GDC_GENOME_FASTA_URL
    output:
        GDC_GENOME_FASTA_TARGZ_FILE
    shell:
        "wget -O {output} {params}"  # | tar -xvz > {output}

rule download_gdc_gtf:
    params:
        GDC_GTF_URL
    output:
        GDC_GTF_FILE
    shell:
        "wget -O - {params} | gunzip -c > {output}"

rule download_gencode_transcript_fasta:
    params:
        GENCODE_TRANSCRIPT_FASTA_URL
    output:
        GENCODE_TRANSCRIPT_FASTA_FILE
    shell:
        "wget {params} -O {output}"

rule download_gencode_genome_fasta:
    params:
        GENCODE_GENOME_FASTA_URL
    output:
        GENCODE_GENOME_FASTA_FILE
    shell:
        "wget -O - {params} | gunzip -c > {output}"

rule download_gencode_gtf:
    params:
        GENCODE_GTF_URL
    output:
        GENCODE_GTF_FILE
    shell:
        "wget -O - {params} | gunzip -c > {output}"
