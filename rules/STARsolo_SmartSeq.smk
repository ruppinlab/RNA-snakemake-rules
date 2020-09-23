from os.path import join
import json


# Directories
ENV_DIR = join("..", "envs")
STAR_OUTPUT_DIR = join("output", "star", "{patient}-{plate}")
STAR_GENOME_INDEX = join("output", "star-index")

# Files
STAR_ENV_FILE = join(ENV_DIR, "star.yml")
STAR_BAM_FILE = join(STAR_OUTPUT_DIR, "Aligned.out.bam")
PE_MANIFEST_FILE = join("output", "{patient}-{plate}")

rule generate_manifest_file:
    input:
        config["units"]
    output:
        PE_MANIFEST_FILE
    script:
        "src/generate_manifest_file.py"

def get_SRA_pe_fq_files(wildcards):
    cells = cells.loc[(cells.patient == wildcards.patient) and (cells.plate == wildcards.plate)]
    return {
        "FQ1": expand(TRIMMED_FASTQ1_FILE, patient=cells["patient"], sample=cells["sample"], cell=cells["cell"]),
        "FQ2": expand(TRIMMED_FASTQ2_FILE, patient=cells["patient"], sample=cells["sample"], cell=cells["cell"]),
    }

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

rule STARsolo_smartseq2_PE:
    group:
        "STAR_2pass"
    conda:
        STAR_ENV_FILE
    input:
        index = STAR_GENOME_INDEX,
        manifest = PE_MANIFEST_FILE,
        unpack(get_input_pe_fastq_files)
    params:
        odir = join(STAR_OUTPUT_DIR, ""),
    output:
        STAR_BAM_FILE
    threads:
        16
    shell:
        "STAR "
        "--runThreadN {threads} "
        "--soloType SmartSeq "
        "--readFilesManifest {input.manifest} "
        "--soloUMIdedup Exact --soloStrand Unstranded "
        "--genomeDir '{input.index}' "
        "--outSAMunmapped Within "
        "--readFilesCommand zcat "
        "--outSAMtype BAM Unsorted "
