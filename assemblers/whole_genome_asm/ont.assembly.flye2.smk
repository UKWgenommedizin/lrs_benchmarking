# ************************************************************************************************
# Flye ONT whole-genome assembly
# ************************************************************************************************

#################
# Include shared assembler header
include: "../../header_assembler.smk"

#################
# Flye version
FLYE_VERSION = "2.9.6"
DOCKER_FLYE = "nicolasardila1/lrs-flye2:" + FLYE_VERSION

print("Flye version: " + FLYE_VERSION)

#####################
#Discover datasets and create wildcards
DATASETS_FASTQ, = glob_wildcards(CWD + r"/fastq/{dataset,[A-Za-z0-9._-]+}.fastq.gz")
DATASETS = [
    dataset
    for dataset in DATASETS_FASTQ
    if ".ont." in dataset.lower()]


##############
#Targets

OUTPUT = []

OUTPUT = OUTPUT + expand(CWD + "/assemblies/flye/{dataset}/assembly.fasta",zip,dataset=DATASETS)

rule all:
    input:
        OUTPUT

print ("Discover datasets and create wildcards")
print(OUTPUT)

################
# Flye resource requirements

def get_flye_memory(wildcards):
    dataset = wildcards.dataset.lower()

    if ".1k" in dataset:
        return 8

    if ".chr21." in dataset:
        return 8

    return 450

################
#Rules

rule flye_assemble:
    input:
        fastq = CWD + "/fastq/{dataset}.fastq.gz"

    output:
        fasta = CWD + "/assemblies/flye/{dataset}/assembly.fasta"

    params:
        outdir = CWD + "/assemblies/flye/{dataset}"
    
    log:
        CWD + "/assemblies/flye/{dataset}/flye.log"

    message:
        "executing {rule} with output {output} and input {input}"

    threads: 12

    resources:
        mem_gb = get_flye_memory

    shell:
        """
        mkdir -p "$(dirname "{log}")"

        (
            set -eo pipefail

            echo "[$(date -Is)] START flye_assemble {wildcards.dataset}"
            echo "Dataset: {wildcards.dataset}"
            echo "Read technology: ONT"
            echo "Read mode: --nano-hq"
            echo "Threads: {threads}"
            echo "Memory: {resources.mem_gb} GB"

            mkdir -p "{CWD}/assemblies/flye"

            docker run --rm \
                --workdir /tmp \
                -u $UID:$(id -g) \
                --cpus {threads} \
                -m {resources.mem_gb}g \
                -v {CWD}:{CWD} \
                {DOCKER_FLYE} \
                flye \
                --nano-hq \
                {input.fastq} \
                --out-dir {params.outdir} \
                --threads {threads}

            [[ -s "{output.fasta}" ]] || {{
                echo "ERROR: Flye assembly.fasta is missing or empty"
                exit 101;
            }}

            grep -q '^>' "{output.fasta}" || {{
                echo "ERROR: Flye output does not appear to be a valid FASTA"
                exit 101;
            }}

            echo "[$(date -Is)] END flye_assemble {wildcards.dataset}"

        ) > "{log}" 2>&1
        """