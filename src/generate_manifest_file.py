import pandas as pd

df = pd.read_csv(snakemake.input[0], sep="\t")

df["fq1"] = df.apply(lambda x: "FASTQ/trimmed/{}-{}-{}_1.fastq.gz".format(x["patient"], x["sample"], x["cell"]), axis=1)
df["fq2"] = df.apply(lambda x: "FASTQ/trimmed/{}-{}-{}_1.fastq.gz".format(x["patient"], x["sample"], x["cell"]), axis=1)
df[["fq1", "fq2", "cell"]].to_csv(snakemake.output[0], sep="\t", index=False, header=False)
