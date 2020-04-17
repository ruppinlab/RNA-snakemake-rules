
ERCC92_URL = "https://tools.thermofisher.com/content/sfs/manuals/ERCC92.zip"
ERCC92_FA = join("data", "ERCC92.fa")
ERCC92_GTF = join("data", "ERCC92.gtf")

localrules: download_ERCC92_spike_ins, combine_human_ERCC92_GTF, combine_human_ERCC92_fasta

rule combine_human_ERCC92_fasta:
   input:
       config["human_ref"]["genome"],
       ERCC92_FA
   output:
       config["ref"]["genome"]
   shell:
       "cat {input[0]} {input[1]} > {output}"

rule combine_human_ERCC92_GTF:
   input:
       config["human_ref"]["annotation"],
       ERCC92_GTF
   output:
       config["ref"]["annotation"]
   shell:
       "cat {input[0]} {input[1]} > {output}"

rule download_ERCC92_spike_ins:
    params:
        url=ERCC92_URL,
        odir="data"
    output:
        zip=temp("ERCC92.zip"),
        gtf=ERCC92_GTF,
        fa=ERCC92_FA
    shell:
        "wget {params.url} && unzip {output.zip} -d {params.odir}"
