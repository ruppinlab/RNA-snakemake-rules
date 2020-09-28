from os.path import join


# Directories
ENV_DIR = join("..", "envs")
STAR_OUTPUT_DIR = join("output", "star", "{patient}-{sample}-{plate}")
STAR_PE_OUTPUT_DIR = join("output", "star", "_STARpe", "{patient}-{sample}-{plate}")
STAR_SE_OUTPUT_DIR = join("output", "star", "_STARse", "{patient}-{sample}-{plate}")
STAR_GENOME_INDEX = join("output", "star-index")

# Files
STAR_ENV_FILE = join(ENV_DIR, "star.yml")
PE_MANIFEST_FILE = join("output", "manifest", "{patient}-{sample}-{plate}-pe-manifest.tsv")
SE_MANIFEST_FILE = join("output", "manifest", "{patient}-{sample}-{plate}-se-manifest.tsv")

# output files
STAR_PE_BAM_FILE = join(STAR_OUTPUT_DIR, "_STARpe", "Aligned.out.bam")
STAR_SE_BAM_FILE = join(STAR_OUTPUT_DIR, "_STARse", "Aligned.out.bam")

localrules: generate_se_manifest_file, generate_pe_manifest_file

rule generate_se_manifest_file:
    input:
        config["units"]
    output:
        SE_MANIFEST_FILE
    script:
        "../src/generate_se_manifest_file.py"

rule generate_pe_manifest_file:
    input:
        config["units"]
    output:
        PE_MANIFEST_FILE
    script:
        "../src/generate_pe_manifest_file.py"


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


def get_pe_fq_files(wildcards):
    c = cells.loc[(cells.patient == wildcards.patient) & (cells["sample"] == wildcards.sample) & (cells.plate == wildcards.plate)]
    return {
        "FQ1": expand(TRIMMED_FASTQ1_FILE, zip, patient=c["patient"], sample=c["sample"], cell=c["cell"]),
        "FQ2": expand(TRIMMED_FASTQ2_FILE, zip, patient=c["patient"], sample=c["sample"], cell=c["cell"]),
    }

rule STAR_manifest_PE:
    group:
        "STAR_2pass"
    conda:
        STAR_ENV_FILE
    input:
        STAR_GENOME_INDEX,
        PE_MANIFEST_FILE,
        unpack(get_pe_fq_files)
    params:
        odir = join(STAR_PE_OUTPUT_DIR, ""),
    output:
        STAR_PE_BAM_FILE
    threads:
        16
    shell:
        "STAR "
        "--runThreadN {threads} "
        "--readFilesManifest {input[1]} "
        "--genomeDir '{input[0]}' "
        "--outSAMunmapped Within "
        "--readFilesCommand zcat "
        "--outSAMtype BAM Unsorted "
        "--outFileNamePrefix '{params.odir}' "


def get_se_fq_files(wildcards):
    c = cells.loc[(cells.patient == wildcards.patient) & (cells["sample"] == wildcards.sample) & (cells.plate == wildcards.plate)]
    return {
        "FQ1": expand(TRIMMED_UNPAIRED_FILE, zip, patient=c["patient"], sample=c["sample"], cell=c["cell"]),
    }

rule STAR_manifest_SE:
    group:
        "STAR_2pass"
    conda:
        STAR_ENV_FILE
    input:
        STAR_GENOME_INDEX,
        SE_MANIFEST_FILE,
        unpack(get_se_fq_files)
    params:
        odir = join(STAR_SE_OUTPUT_DIR, ""),
    output:
        STAR_SE_BAM_FILE
    threads:
        16
    shell:
        "STAR "
        "--runThreadN {threads} "
        "--readFilesManifest {input[1]} "
        "--genomeDir '{input[0]}' "
        "--outSAMunmapped Within "
        "--readFilesCommand zcat "
        "--outSAMtype BAM Unsorted "
        "--outFileNamePrefix '{params.odir}' "
