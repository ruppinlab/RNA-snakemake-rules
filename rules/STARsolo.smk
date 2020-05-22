from os.path import join

# Directories
ENV_DIR = join("..", "envs")
STAR_GENOME_INDEX = config["STARsolo"]["genome_index"]
STAR_SOLO_OUTPUT = config["STARsolo"]["output_dir"]

# Files
STAR_ENV_FILE = join(ENV_DIR, "star.yml")


rule create_star_index:
    conda:
        STAR_ENV_FILE
    input:
        config["ref"]["genome"]
    output:
        directory(STAR_GENOME_INDEX)
    threads:
        16
    shell:
        "mkdir '{output}' && STAR "
        "--runThreadN {threads} "
        "--runMode genomeGenerate "
        "--genomeDir '{output}' "
        "--genomeFastaFiles '{input}'"

# STARsolo
rule run_STARsolo_10x:
    conda:
        STAR_ENV_FILE
    input:
        barcode_whitelist = config["ref"]["barcode_whitelist"],
        cDNA_fq = cDNA_FASTQ_FILE,
        CellBarcode_fq = CellBarcode_FASTQ_FILE,
        index = STAR_GENOME_INDEX,
        gtf = config["ref"]["annotation"]
    output:
        join(STAR_SOLO_OUTPUT, "Aligned.sortedByCoord.out.bam"),
        join(STAR_SOLO_OUTPUT, "Solo.out", "Gene", "filtered", "matrix.mtx")
    params:
        odir = join(STAR_SOLO_OUTPUT, "")
    threads:
        48
    shell:
        "STAR "
        "--runThreadN {threads} "
        "--genomeDir '{input.index}' "
        "--genomeLoad NoSharedMemory "
        "--outFileNamePrefix '{params.odir}' "
        "--readFilesIn '{input.cDNA_fq}' '{input.CellBarcode_fq}' "
        "--soloType CB_UMI_Simple "
        "--soloCBwhitelist {input.barcode_whitelist} "
        "--soloFeatures Gene "
        "--outSAMattrRGline ID:{wildcards.patient}.{wildcards.sample} "
        "PL:illumina SM:{wildcards.patient}.{wildcards.sample} LB:RNA "
        # CR/UR is raw (uncorrected) CellBarcode/UMI while CB/UB is corrected
        "--outSAMattributes NH HI AS nM NM ch CR UR CB UB "
        "--outSAMtype BAM SortedByCoordinate "  # required for CB/UB
        "--outSAMunmapped Within "  # include the unmapped reads in the output file
        "--readFilesCommand zcat "  # used for reading gz files
        "--sjdbGTFfile '{input.gtf}' "
        + config["params"]["STARsolo"]
