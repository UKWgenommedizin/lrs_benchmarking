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
# Discover datasets and create wildcards

DATASETS_FASTQ, = glob_wildcards(
    CWD + r"/fastq/{dataset,[A-Za-z0-9._-]+}.fastq.gz"
)

# Only ONT whole-genome datasets are included.
# Local 1k and chromosome 21 validation datasets are intentionally excluded. 
DATASETS = [
    dataset
    for dataset in DATASETS_FASTQ
    if ".ont." in dataset.lower()
    and ".1k" not in dataset.lower()
    and ".chr21." not in dataset.lower()
]


##############
# Targets

OUTPUT = []

OUTPUT += expand(
    CWD + "/assemblies/flye/{dataset}/assembly.fasta",
    zip,
    dataset=DATASETS
)

rule all:
    input:
        OUTPUT

print("Discover datasets and create wildcards")
print(OUTPUT)


################
# Prevent local test datasets from being requested explicitly

wildcard_constraints:
    dataset = r"(?=.*\.ont\.)(?!.*\.1k(?:\.|$))(?!.*\.chr21\.)[A-Za-z0-9._-]+"


################
# Flye resource requirements

# The server provides approximately 128 CPU threads and 256 GB RAM total.
# A maximum of 240 GB is used for Flye to leave memory available
# for the operating system and other server processes.
#
# Important:
# Flye documentation reports that human ONT WGS assemblies may require 450GB RAM depending on coverage and read characteristics.
def get_flye_memory(wildcards):
    return 240000


################
# Rules

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

    # Flye recommendations described 30-50 threads.
    # The existing server configuration assigns 32 threads to Flye.
    threads: 32

    resources:
        mem_mb = get_flye_memory

    shell:
        """
        mkdir -p "$(dirname "{log}")"

        (
            set -euo pipefail

            echo "[$(date -Is)] START flye_assemble {wildcards.dataset}"
            echo "Dataset: {wildcards.dataset}"
            echo "Read technology: ONT"
            echo "Read mode: --nano-hq"
            echo "Threads: {threads}"
            echo "Memory: {resources.mem_mb} MB"
            echo "Container hostname: flye-{wildcards.dataset}"

            mkdir -p "{params.outdir}"

            docker run --rm \
                --hostname flye-{wildcards.dataset} \
                --workdir /tmp \
                -u $UID:$(id -g) \
                --cpus {threads} \
                -m {resources.mem_mb}m \
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