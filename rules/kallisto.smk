from os.path import join

KALLISTO_ENV_FILE = join("..", "envs", "kallisto.yaml")
KALLISTO_INDEX_FILE = join("results", "kallisto", "transcripts.idx")

def get_fq(wildcards):
    if config["trimming"]["skip"]:
        # no trimming, use raw reads
        dir = "data"
    else:
        dir = "trimmed"
    return {
        'fq1': join(dir, "{wildcards.patient}.{wildcards.sample}_1.fastq.gz".format(wildcards=wildcards)),
        'fq2': join(dir, "{wildcards.patient}.{wildcards.sample}_2.fastq.gz".format(wildcards=wildcards))
        }

def kallisto_params(wildcards, input):
    extra = config["params"]["kallisto"]
    if "strandedness" in samples.columns:
        strandedness = samples.loc[(wildcards.patient, wildcards.sample),
                                   "strandedness"]
        if strandedness == "reverse":
            extra += " --rf-stranded"
        if strandedness == "yes":
            extra += " --fr-stranded"
    return extra


rule create_kallisto_index:
    conda:
        KALLISTO_ENV_FILE
    log:
        "results/logs/kallisto/index.log"
    input:
        config["ref"]["transcriptome"]
    output:
        KALLISTO_INDEX_FILE
    shell:
        "kallisto index -i {output} {input} 2> {log}"

rule run_kallisto:
    conda:
        KALLISTO_ENV_FILE
    log:
        "results/logs/kallisto/quant/{sample}-{unit}.log"
    input:
        unpack(get_fq),
        idx = KALLISTO_INDEX_FILE
    output:
        directory("results/kallisto/{patient}-{sample}")
    params:
        extra = kallisto_params
    shell:
        "kallisto quant -i '{input.idx}' -o '{output}' {params.extra} "
        "'{input.fq1}' '{input.fq2}' 2> {log}"
