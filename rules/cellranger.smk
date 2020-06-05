from os.path import join

PATIENT_FASTQ_DIR = join("FASTQ", "raw", "{patient}")
CR_SAMPLE_ODIR = join("output", "CellRanger", "{patient}-{sample}")

CB_FASTQ_FILE = join(PATIENT_FASTQ_DIR, "{sample}_{lane}_R1_001.fastq.gz")
cDNA_FASTQ_FILE = join(PATIENT_FASTQ_DIR, "{sample}_{lane}_R2_001.fastq.gz")
IDX_FASTQ_FILE = join(PATIENT_FASTQ_DIR, "{sample}_{lane}_I1_001.fastq.gz")
CR_BAM_FILE = join(CR_SAMPLE_ODIR, "outs", "possorted_genome_bam.bam")

def get_cellranger_fq_files(wildcards):
    lanes_of_interest = lanes.loc[(wildcards.patient, wildcards.sample), "lane"]
    return {
        "CB": expand(CB_FASTQ_FILE, patient=wildcards.patient, sample=wildcards.sample, lane=lanes_of_interest),
        "cDNA": expand(cDNA_FASTQ_FILE, patient=wildcards.patient, sample=wildcards.sample, lane=lanes_of_interest),
        "IDX": expand(IDX_FASTQ_FILE, patient=wildcards.patient, sample=wildcards.sample, lane=lanes_of_interest),
    }

# expected input format for FASTQ file
rule cellranger_count:
    input:
        unpack(get_cellranger_fq_files)
    params:
        PATIENT_FASTQ_DIR,
        CR_SAMPLE_ODIR
    output:
        CR_BAM_FILE
    shell:
        "module load cellranger && "
        "cellranger count --id={params[1]} "
        "--fastqs={params[0]} " # this is the path to the directory containing the FASTQ files
        "--sample={wildcards.sample} " # this is the sample to use
        "--transcriptome=$CELLRANGER_REF300/GRCh38 "
        "--localcores=$SLURM_CPUS_PER_TASK "
        "--chemistry=SC3Pv2 "
        "--localmem=60"
