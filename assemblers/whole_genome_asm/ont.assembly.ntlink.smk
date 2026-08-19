# ************************************************************************************************
# ntLink ONT whole-genome scaffolding
# ************************************************************************************************

#################
# Include shared assembler header
include: "../../header_assembler.smk"

#################
# ntLink version
NTLINK_VERSION = "1.3.11"
DOCKER_NTLINK = "nicolasardila1/lrs-ntlink:" + NTLINK_VERSION

print("ntLink version: " + NTLINK_VERSION)

#####################
# Discover datasets and create wildcards

DATASETS_FASTQ, = glob_wildcards(CWD + r"/fastq/{dataset,[A-Za-z0-9._-]+}.fastq.gz")

# Only ONT whole-genome datasets are included.
# Local 1k and chromosome 21 validation datasets are intentionally excluded.
DATASETS = [
    dataset
    for dataset in DATASETS_FASTQ
    if ".ont." in dataset.lower()
    and ".1k" not in dataset.lower()
    and ".chr21." not in dataset.lower()
    and ".localtest." not in dataset.lower()]

##############
# Targets

OUTPUT = []

OUTPUT += expand(CWD + "/assemblies/ntlink/{dataset}/assembly.fasta",zip,dataset=DATASETS)

rule all:
    input:
        OUTPUT

print("Discover datasets and create wildcards")
print(OUTPUT)


################
# Prevent local test datasets from being requested explicitly

wildcard_constraints:
    dataset = r"(?=.*\.ont\.)(?!.*\.1k(?:\.|$))(?!.*\.chr21\.)(?!.*\.localtest\.)[A-Za-z0-9._-]+"

################
# ntLink resource requirements

# The server has a maximum of 128 CPU threads and 256 GB RAM total.
# The existing server configuration assigns 16 threads and 32 GB RAM to ntLink
def get_ntlink_memory(wildcards):
    return 32000

################
# Rules

rule ntLink_scaffold:
    input:
        draft = CWD + "/assemblies/goldrush/{dataset}/assembly.fasta",
        fastq = CWD + "/fastq/{dataset}.fastq.gz"

    output:
        fasta = CWD + "/assemblies/ntlink/{dataset}/assembly.fasta"

    params:
        outdir = CWD + "/assemblies/ntlink/{dataset}"

    log:
        CWD + "/assemblies/ntlink/{dataset}/ntlink.log"

    message:
        "executing {rule} with output {output} and input {input}"

    threads: 16

    resources:
        mem_mb = get_ntlink_memory

    shell:
        r"""
        mkdir -p "$(dirname "{log}")"

        (
            set -euo pipefail

            echo "[$(date -Is)] START ntLink_scaffold {wildcards.dataset}"
            echo "Dataset: {wildcards.dataset}"
            echo "Read technology: ONT"
            echo "Draft assembly: {input.draft}"
            echo "Long reads: {input.fastq}"
            echo "Threads: {threads}"
            echo "Memory: {resources.mem_mb} MB"
            echo "Container hostname: ntLink-{wildcards.dataset}"

            mkdir -p "{params.outdir}"

            # Prepare the draft assembly inside the ntLink working directory.
            cp "{input.draft}" "{params.outdir}/draft.fasta"

            # Prepare the long-read FASTQ inside the ntLink working directory
            READS_FASTQ="{params.outdir}/reads.fastq"

            echo "Preparing uncompressed ntLink input..."
            gzip -cd "{input.fastq}" > "$READS_FASTQ.tmp"
            mv "$READS_FASTQ.tmp" "$READS_FASTQ"

            docker run --rm \
            --hostname ntLink-{wildcards.dataset} \
            --workdir {params.outdir} \
            -u $UID:$(id -g) \
            --cpus {threads} \
            -m {resources.mem_mb}m \
            -v {CWD}:{CWD} \
            {DOCKER_NTLINK} \
            ntLink scaffold \
                target=draft.fasta \
                reads='reads.fastq' \
                t={threads} \
                k=32 \
                w=100 \
                z=1000

            NTLINK_ASSEMBLY="{params.outdir}/draft.fasta.k32.w100.z1000.ntLink.scaffolds.fa"

            [[ -s "$NTLINK_ASSEMBLY" ]] || {{
                echo "ERROR: ntLink final assembly is missing or empty"
                echo "Expected: $NTLINK_ASSEMBLY"
                exit 101;
            }}

            cp "$NTLINK_ASSEMBLY" "{output.fasta}"

            [[ -s "{output.fasta}" ]] || {{
                echo "ERROR: ntLink assembly.fasta is missing or empty"
                exit 101;
            }}

            grep -q '^>' "{output.fasta}" || {{
                echo "ERROR: ntLink output does not appear to be a valid FASTA"
                exit 101;
            }}

            echo "[$(date -Is)] END ntLink_scaffold {wildcards.dataset}"

        ) > "{log}" 2>&1
        """