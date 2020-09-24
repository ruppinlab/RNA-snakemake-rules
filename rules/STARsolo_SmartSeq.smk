from os.path import join


# Directories
ENV_DIR = join("..", "envs")
STAR_OUTPUT_DIR = join("output", "star", "{patient}-{sample}-{plate}")
IND_STAR_OUTPUT_DIR = join("output", "star", "{patient}", "{sample}-{cell}")
STAR_GENOME_INDEX = join("output", "star-index")

# Files
STAR_ENV_FILE = join(ENV_DIR, "star.yml")
STAR_BAM_FILE = join(STAR_OUTPUT_DIR, "Aligned.out.bam")
IND_STAR_BAM_FILE = join(IND_STAR_OUTPUT_DIR, "Aligned.out.bam")
PE_MANIFEST_FILE = join("output", "{patient}-{sample}-{plate}-manifest.tsv")

rule generate_manifest_file:
    input:
        config["units"]
    output:
        PE_MANIFEST_FILE
    script:
        "../src/generate_manifest_file.py"

def get_pe_fq_files(wildcards):
    c = cells.loc[(cells.patient == wildcards.patient) & (cells.plate == wildcards.plate)]
    return {
        "FQ1": expand(TRIMMED_FASTQ1_FILE, zip, patient=c["patient"], sample=c["sample"], cell=c["cell"]),
        "FQ2": expand(TRIMMED_FASTQ2_FILE, zip, patient=c["patient"], sample=c["sample"], cell=c["cell"]),
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

rule STAR_PE:
    group:
        "STAR_2pass"
    conda:
        STAR_ENV_FILE
    input:
        STAR_GENOME_INDEX,
        TRIMMED_FASTQ1_FILE,
        TRIMMED_FASTQ2_FILE
    params:
        odir = join(IND_STAR_OUTPUT_DIR, ""),
    output:
        IND_STAR_BAM_FILE
    threads:
        16
    shell:
        "STAR "
        "--runThreadN {threads} "
        "--readFilesIn '{input[1]}' '{input[2]}' "
        "--genomeDir '{input[0]}' "
        "--outSAMunmapped Within "
        "--readFilesCommand zcat "
        "--outSAMtype BAM Unsorted "
        "--outFileNamePrefix '{params.odir}' "
        "--outSAMattrRGline ID:{wildcards.cell}"

rule STARsolo_smartseq2_PE:
    group:
        "STAR_2pass"
    conda:
        STAR_ENV_FILE
    input:
        STAR_GENOME_INDEX,
        PE_MANIFEST_FILE,
        unpack(get_pe_fq_files)
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
        "--readFilesManifest {input[1]} "
        "--soloUMIdedup Exact --soloStrand Unstranded "
        "--genomeDir '{input[0]}' "
        "--outSAMunmapped Within "
        "--readFilesCommand zcat "
        "--outSAMtype BAM Unsorted "
        "--outFileNamePrefix '{params.odir}' "
