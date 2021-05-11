from os.path import join

# files
GENOME_FA = join("raw", "{genome}.fa")
SRPRISM_DB_DIR = join("output", "SRPRISM", "genomes_db")
GENOME_DB_PREFIX = join(SRPRISM_DB_DIR, "{genome}")
GENOME_DB_FILE = join(SRPRISM_DB_DIR, "{genome}.idx")

rule make_reference_DB:
    params:
        GENOME_DB_PREFIX
    input:
        GENOME_FA
    output:
        GENOME_DB_FILE
    shell:
        "srprism mkindex -i {input} -o {params}"
