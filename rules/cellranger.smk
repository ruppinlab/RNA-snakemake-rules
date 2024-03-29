from os.path import join

PATIENT_FASTQ_DIR = join("FASTQ", "raw", "{patient}")

# cellranger complains when you pass directory as --id
CR_SAMPLE_ODIR = "{patient}-{sample}"

CB_FASTQ_FILE = join(PATIENT_FASTQ_DIR, "{sample}_{lane}_R1_001.fastq.gz")
cDNA_FASTQ_FILE = join(PATIENT_FASTQ_DIR, "{sample}_{lane}_R2_001.fastq.gz")
CR_BAM_FILE = join(CR_SAMPLE_ODIR, "outs", "possorted_genome_bam.bam")
CR_BARCODES_FILE = join(CR_SAMPLE_ODIR, "outs", "filtered_feature_bc_matrix", "barcodes.tsv.gz") 
CR_RAW_BARCODES_FILE = join(CR_SAMPLE_ODIR, "outs", "raw_feature_bc_matrix", "barcodes.tsv.gz")

def get_cellranger_fq_files(wildcards):
    lanes_of_interest = lanes.loc[(wildcards.patient, wildcards.sample), "lane"]
    return {
        "CB": expand(CB_FASTQ_FILE, patient=wildcards.patient, sample=wildcards.sample, lane=lanes_of_interest),
        "cDNA": expand(cDNA_FASTQ_FILE, patient=wildcards.patient, sample=wildcards.sample, lane=lanes_of_interest),
    }

# expected input format for FASTQ file
rule cellranger_count:
    input:
        unpack(get_cellranger_fq_files)
    params:
        PATIENT_FASTQ_DIR,
        CR_SAMPLE_ODIR,
        config["CellRanger"]["genome_dir"],
        config["CellRanger"]["chemistry"]
    output:
        CR_BAM_FILE,
        CR_BARCODES_FILE,
        CR_RAW_BARCODES_FILE
    shell:
        "module load CellRanger/6.0.1 && "
        # snakemake auto creates directories for output files but cellranger expects existing directories to pipestance directory
        "rm -rf {params[1]} && "
        "cellranger count --id={params[1]} "
        "--fastqs={params[0]} " # this is the path to the directory containing the FASTQ files
        "--sample={wildcards.sample} " # this is the sample to use
        "--transcriptome=refdata-gex-GRCh38-2020-A "
        "--localcores=4 "
        "--chemistry={params[3]} "
        "--localmem=60"
