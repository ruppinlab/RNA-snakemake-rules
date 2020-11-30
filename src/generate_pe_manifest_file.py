import pandas as pd

df = pd.read_csv(snakemake.input[0], sep="\t", dtype={"patient": "str", "sample": "str", "plate": "str"})
df = df.loc[(df.plate == snakemake.wildcards["plate"]) & (df["sample"] == snakemake.wildcards["sample"]) & (df.patient == snakemake.wildcards["patient"])]
df["fq1"] = df.apply(lambda x: "FASTQ/trimmed/{}-{}-{}_1.fastq.gz".format(x["patient"], x["sample"], x["cell"]), axis=1)
df["fq2"] = df.apply(lambda x: "FASTQ/trimmed/{}-{}-{}_2.fastq.gz".format(x["patient"], x["sample"], x["cell"]), axis=1)
df["RG_ID"] = df["cell"].apply(lambda x: "ID:{}".format(x))
df["RG_PL"] = "PL:illumina"
df["RG_SM"] = df["cell"].apply(lambda x: "SM:{}".format(x))
df["RG_LB"] = "LB:RNA"
df[["fq1", "fq2", "RG_ID", "RG_PL", "RG_SM", "RG_LB"]].to_csv(snakemake.output[0], sep="\t", index=False, header=False)
