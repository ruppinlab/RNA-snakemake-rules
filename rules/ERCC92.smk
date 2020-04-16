
ERCC92_URL = "https://tools.thermofisher.com/content/sfs/manuals/ERCC92.zip"
ERCC92_FA = join("data", "ERCC92.fa")
ERCC92_GTF = join("data", "ERCC92.gtf")
localrules: download_ERCC92_spike_ins

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
