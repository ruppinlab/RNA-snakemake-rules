from os.path import join

# GDC genomes are used in GDC for all sequencing and array-based analysis (including TCGA)
# They include
## GCA_000001405.15_GRCh38_no_alt_analysis_set
## Sequence Decoys (GenBank Accession GCA_000786075)
## Virus Sequences (https://gdc.cancer.gov/files/public/file/GRCh83.d1.vd1_virus_decoy.txt)

# There are multiple versions of the genome from GENCODE
## Genome Sequence - contains chromosomes, scaffolds, assembly patches and haplotypes (no EBV)
## Genome Sequence, primary assembly - contains chromosomes and scaffolds (no EBV)
## Gene annotation (CHR) - only for genes on the reference chromosomes
## Gene annotation (ALL) - gene annotations on the reference chromosomes, scaffolds, assembly patches and alternate loci

# URLs for downloading reference sequence and annotation files
# Links for downloading GENCODE files
GENCODE_HUMAN_BASE_URL = (
    "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34")
GENCODE_CHR_GTF_URL = join(GENCODE_HUMAN_BASE_URL, "gencode.v34.annotation.gtf.gz")
GENCODE_PRI_GTF_URL = join(GENCODE_HUMAN_BASE_URL, "gencode.v34.primary_assembly.annotation.gtf.gz")
GENCODE_ALL_GTF_URL = join(GENCODE_HUMAN_BASE_URL, "gencode.v34.chr_patch_hapl_scaff.annotation.gtf.gz")
# ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh38.p13.genome.fa.gz
GENCODE_PRIMARY_GENOME_FASTA_URL = join(GENCODE_HUMAN_BASE_URL, "GRCh38.primary_assembly.genome.fa.gz")
GENCODE_GENOME_FASTA_URL = join(GENCODE_HUMAN_BASE_URL, "GRCh38.p13.genome.fa.gz")

# Links for downloading GDC files
GDC_DATA_BASE_URL = "https://api.gdc.cancer.gov/data"
GDC_GENOME_FASTA_URL = join(GDC_DATA_BASE_URL,
                            "254f697d-310d-4d7d-a27b-27fbf767a834")
GDC_GTF_URL = join(GDC_DATA_BASE_URL, "25aa497c-e615-4cb7-8751-71f744f9691f")

# Directories
GENOME_DIR = join("raw", "genome")

# files
# GENCODE files
GENCODE_CHR_GTF = join(GENOME_DIR, "gencode.v34.annotation.gtf")
GENCODE_PRI_GTF = join(GENOME_DIR, "gencode.v34.primary_assembly.annotation.gtf")
GENCODE_ALL_GTF = join(GENOME_DIR, "gencode.v34.chr_patch_hapl_scaff.annotation.gtf")
GENCODE_GENOME_FASTA = join(GENOME_DIR, "GRCh38.p13.genome.fa")
GENCODE_PRIMARY_GENOME_FASTA = join(GENOME_DIR, "GRCh38.primary_assembly.genome.fa")

# GDC files
GDC_GENOME_FASTA_FILE = join(GENOME_DIR, "GRCh38.d1.vd1.fa")
GDC_GTF_FILE = join(GENOME_DIR, "gencode.v22.annotation.gtf")



# rules for downloading GDC reference sequence and annotation files
rule download_gdc_genome_fasta:
    params:
        url=GDC_GENOME_FASTA_URL,
        odir=GENOME_DIR
    output:
        GDC_GENOME_FASTA_FILE
    shell:
        "curl {params.url} | tar xvz -C {params.odir}"

rule download_gdc_gtf:
    params:
        GDC_GTF_URL
    output:
        GDC_GTF_FILE
    shell:
        "wget -O - {params} | gunzip -c > {output}"

# rules for downloading GENCODE files
rule download_gencode_genome_fasta:
    params:
        GENCODE_GENOME_FASTA_URL
    output:
        GENCODE_GENOME_FASTA
    shell:
        "wget -O - {params} | gunzip -c > {output}"

rule download_gencode_primary_genome_fasta:
    params:
        GENCODE_PRIMARY_GENOME_FASTA_URL
    output:
        GENCODE_PRIMARY_GENOME_FASTA
    shell:
        "wget -O - {params} | gunzip -c > {output}"

rule download_gencode_CHR_gtf:
    params:
        GENCODE_CHR_GTF_URL
    output:
        GENCODE_CHR_GTF
    shell:
        "wget -O - {params} | gunzip -c > {output}"

rule download_gencode_ALL_gtf:
    params:
        GENCODE_ALL_GTF_URL
    output:
        GENCODE_ALL_GTF
    shell:
        "wget -O - {params} | gunzip -c > {output}"

rule download_gencode_PRI_gtf:
    params:
        GENCODE_PRI_GTF_URL
    output:
        GENCODE_PRI_GTF
    shell:
        "wget -O - {params} | gunzip -c > {output}"
