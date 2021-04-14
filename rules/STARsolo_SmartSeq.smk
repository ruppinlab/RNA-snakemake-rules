from os.path import join


# Directories
ENV_DIR = join("..", "envs")
STAR_OUTPUT_DIR = join("output", "star", "{patient}-{sample}-{plate}")
STAR_PE_OUTPUT_DIR = join(STAR_OUTPUT_DIR, "_STARpe")
STAR_SE_OUTPUT_DIR = join(STAR_OUTPUT_DIR, "_STARse")
STAR_GENOME_INDEX = join("output", "star-index")

# Files
STAR_ENV_FILE = join(ENV_DIR, "star.yml")
PE_MANIFEST_FILE = join("output", "manifest", "{patient}-{sample}-{plate}-pe-manifest.tsv")
SE_MANIFEST_FILE = join("output", "manifest", "{patient}-{sample}-{plate}-se-manifest.tsv")

# output files
STAR_PE_BAM_FILE = join(STAR_PE_OUTPUT_DIR, "Aligned.out.bam")
STAR_SE_BAM_FILE = join(STAR_SE_OUTPUT_DIR, "Aligned.out.bam")
STAR_PE_READCOUNT_FILE = join(STAR_PE_OUTPUT_DIR, "ReadsPerGene.out.tab")

localrules: generate_se_manifest_file, generate_pe_manifest_file

# function
def get_SAMattrRGline(wildcards):
    c = cells.loc[(cells.patient == wildcards.patient) & (cells["sample"] == wildcards.sample) & (cells.plate == wildcards.plate)]
    return " , ".join(["ID:{}".format(x) for x in c["cell"]])

def get_pe_fq_files(wildcards):
    c = cells.loc[(cells.patient == wildcards.patient) & (cells["sample"] == wildcards.sample) & (cells.plate == wildcards.plate)]
    return {
        "FQ1": expand(TRIMMED_FASTQ1_FILE, zip, patient=c["patient"], sample=c["sample"], cell=c["cell"]),
        "FQ2": expand(TRIMMED_FASTQ2_FILE, zip, patient=c["patient"], sample=c["sample"], cell=c["cell"]),
    }

def get_se_fq_files(wildcards):
    c = cells.loc[(cells.patient == wildcards.patient) & (cells["sample"] == wildcards.sample) & (cells.plate == wildcards.plate)]
    return {
        "FQ1": expand(TRIMMED_UNPAIRED_FILE, zip, patient=c["patient"], sample=c["sample"], cell=c["cell"]),
    }

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


rule STAR_manifest_PE:
    conda:
        STAR_ENV_FILE
    input:
        STAR_GENOME_INDEX,
        PE_MANIFEST_FILE,
        unpack(get_pe_fq_files)
    params:
        odir = join(STAR_PE_OUTPUT_DIR, ""),
        SAMattrRGline = get_SAMattrRGline
    output:
        STAR_PE_BAM_FILE
    threads:
        48
    benchmark:
        "benchmarks/{patient}-{sample}-{plate}.STAR_manifest_PE.benchmark.txt"
    shell:
        "STAR "
        "--soloType SmartSeq --soloUMIdedup Exact --soloStrand Unstranded "
        "--runThreadN {threads} "
        "--readFilesManifest {input[1]} "
        "--genomeDir '{input[0]}' "
        "--outSAMunmapped Within "
        "--readFilesCommand zcat "
        "--outSAMtype BAM Unsorted "
        "--outSAMattrRGline {params.SAMattrRGline} "
        "--outFileNamePrefix '{params.odir}' "
        "--quantMode GeneCounts "


rule STAR_manifest_SE:
    conda:
        STAR_ENV_FILE
    input:
        STAR_GENOME_INDEX,
        SE_MANIFEST_FILE,
        unpack(get_se_fq_files)
    params:
        odir = join(STAR_SE_OUTPUT_DIR, ""),
        SAMattrRGline = get_SAMattrRGline
    output:
        STAR_SE_BAM_FILE
    threads:
        16
    benchmark:
        "benchmarks/{patient}-{sample}-{plate}.STAR_manifest_SE.benchmark.txt"
    shell:
        "STAR "
        "--soloType SmartSeq --soloUMIdedup Exact --soloStrand Unstranded "
        "--runThreadN {threads} "
        "--readFilesManifest {input[1]} "
        "--genomeDir '{input[0]}' "
        "--outSAMunmapped Within "
        "--readFilesCommand zcat "
        "--outSAMtype BAM Unsorted "
        "--outSAMattrRGline {params.SAMattrRGline} "
        "--outFileNamePrefix '{params.odir}' "
        "--quantMode GeneCounts "
