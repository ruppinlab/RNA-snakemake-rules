from os.path import join


# URLs
GENCODE_MOUSE_BASE_URL = (
    "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24")

GENCODE_MOUSE_GENOME_FASTA_URL = join(GENCODE_MOUSE_BASE_URL,
                                "GRCm38.p6.genome.fa.gz")
GENCODE_MOUSE_GTF_URL = join(GENCODE_MOUSE_BASE_URL,
                       "gencode.vM24.annotation.gtf.gz")
# output files
GENCODE_MOUSE_GENOME_FASTA_FILE = join("raw", "genome", "mouse", "GRCm38.p6.genome.fa")
GENCODE_MOUSE_GTF_FILE = join("raw", "genome", "mouse", "gencode.vM24.annotation.gtf")


rule download_gencode_genome_fasta:
    params:
        GENCODE_MOUSE_GENOME_FASTA_URL
    output:
        GENCODE_MOUSE_GENOME_FASTA_FILE
    shell:
        "wget -O - {params} | gunzip -c > {output}"

rule download_gencode_gtf:
    params:
        GENCODE_MOUSE_GTF_URL
    output:
        GENCODE_MOUSE_GTF_FILE
    shell:
        "wget -O - {params} | gunzip -c > {output}"
