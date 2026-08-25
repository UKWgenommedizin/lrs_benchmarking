# ************************************************************************************************
# GoldRush ONT whole-genome assembly
# ************************************************************************************************

#################
# Include shared assembler header
include: "../../header_assembler.smk"

#################
# GoldRush version
GOLDRUSH_VERSION = "1.2.2-ntlinkfix"
DOCKER_GOLDRUSH = "nicolasardila1/lrs-goldrush:" + GOLDRUSH_VERSION

print("GoldRush version: " + GOLDRUSH_VERSION)

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

OUTPUT += expand(CWD + "/assemblies/goldrush/{dataset}/assembly.fasta",zip,dataset=DATASETS)

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
# GoldRush resource requirements

# The server has a maximum of 128 CPU threads and 256 GB RAM total.
# A maximum of 64 GB is used for GoldRush.
#
# Important:
# GoldRush documentation reports for human WGS assembly with approximately
# 60x coverage requires at least 64 GB RAM and recommends at least 48 threads.
# The existing server configuration assigns 32 threads and 64 GB RAM to GoldRush.
def get_goldrush_memory(wildcards):
    return 64000

################
# Rules

rule goldrush_assemble:
    input:
        fastq = CWD + "/fastq/{dataset}.fastq.gz"

    output:
        fasta = CWD + "/assemblies/goldrush/{dataset}/assembly.fasta"

    params:
        outdir = CWD + "/assemblies/goldrush/{dataset}",
        genome_size = "3e9",
        prefix = "{dataset}_goldrush",
        shm_size = "8g"

    log:
        CWD + "/assemblies/goldrush/{dataset}/goldrush.log"

    message:
        "executing {rule} with output {output} and input {input}"

    # GoldRush recommends at least 48 threads for human WGS.
    # The existing server configuration assigns 32 threads to GoldRush.
    threads: 32

    resources:
        mem_mb = get_goldrush_memory

    shell:
        """
        mkdir -p "$(dirname "{log}")"

        # Preserve the previous GoldRush log before starting a new execution.
        if [[ -s "{log}" ]]; then
        PREVIOUS_LOG="{log}.previous.$(date +%Y%m%dT%H%M%S)"
        cp -- "{log}" "$PREVIOUS_LOG"
        echo "Previous GoldRush log saved to: $PREVIOUS_LOG"
        fi

        (
            set -euo pipefail

            echo "[$(date -Is)] START goldrush_assemble {wildcards.dataset}"
            echo "Dataset: {wildcards.dataset}"
            echo "Read technology: ONT"
            echo "Genome size: {params.genome_size}"
            echo "Threads: {threads}"
            echo "Memory: {resources.mem_mb} MB"
            echo "Shared memory: {params.shm_size}"
            echo "Container hostname: goldrush-{wildcards.dataset}"

            mkdir -p "{params.outdir}"

            # GoldRush requires an uncompressed FASTQ file.
            # The reads argument must be provided without the .fastq suffix.
            READS_FASTQ="{params.outdir}/{wildcards.dataset}.fastq"

            # Prepare the uncompressed FASTQ when it does not exist or when the
            # compressed source FASTQ is newer than the prepared copy.
            if [[ ! -s "$READS_FASTQ" || "{input.fastq}" -nt "$READS_FASTQ" ]]; then
                echo "Preparing uncompressed GoldRush input..."

                rm -f -- "$READS_FASTQ.tmp"

                gzip -cd "{input.fastq}" > "$READS_FASTQ.tmp"

            [[ -s "$READS_FASTQ.tmp" ]] || {{
                echo "ERROR: decompressed GoldRush FASTQ is missing or empty"
                exit 102
            }}

                mv "$READS_FASTQ.tmp" "$READS_FASTQ"
            fi

            # -------------------------------------------------------------------------
            # Reset GoldRush internal workflow state.
            #
            # GoldRush manages an internal Make/ntLink workflow whose intermediate
            # files are not declared as individual Snakemake outputs. Files left from
            # a failed or previous run may therefore be reused by GoldRush/ntLink,
            # including stale or incomplete ntLink checkpoints.
            #
            # The prepared uncompressed FASTQ is stored outside this directory and
            # is intentionally preserved.
            # -------------------------------------------------------------------------
            INTERMEDIATE_DIR="{params.outdir}/goldrush_intermediate_files"

            if [[ -d "$INTERMEDIATE_DIR" ]]; then
                echo "Removing previous GoldRush internal state:"
                echo "$INTERMEDIATE_DIR"
                rm -rf -- "$INTERMEDIATE_DIR"
            fi

            rm -f -- "{output.fasta}"

            docker run --rm \
                --hostname goldrush-{wildcards.dataset} \
                --workdir "{params.outdir}" \
                -u $UID:$(id -g) \
                --cpus {threads} \
                -m {resources.mem_mb}m \
                --shm-size {params.shm_size} \
                -v "{CWD}:{CWD}" \
                {DOCKER_GOLDRUSH} \
                goldrush run \
                reads={wildcards.dataset} \
                G={params.genome_size} \
                t={threads} \
                m=5000 \
                P=0 \
                p={params.prefix} \
                track_time=1

            FINAL_ASSEMBLY=$(find \
                "{params.outdir}/goldrush_intermediate_files" \
                -type f \
                -name "*ntLink-5rounds.polished.fa" \
                -print \
                | sort \
                | tail -n 1)

            [[ -n "$FINAL_ASSEMBLY" && -s "$FINAL_ASSEMBLY" ]] || {{
                echo "ERROR: GoldRush final assembly is missing or empty"
                echo "Expected under: {params.outdir}/goldrush_intermediate_files"
                exit 101;
            }}

            cp "$FINAL_ASSEMBLY" "{output.fasta}"

            [[ -s "{output.fasta}" ]] || {{
                echo "ERROR: GoldRush assembly.fasta is missing or empty"
                exit 101;
            }}

            grep -q '^>' "{output.fasta}" || {{
                echo "ERROR: GoldRush output does not appear to be a valid FASTA"
                exit 101;
            }}

            echo "[$(date -Is)] END goldrush_assemble {wildcards.dataset}"
            echo "Assembly size: $(du -h "{output.fasta}" | cut -f1)"

        ) > "{log}" 2>&1
        """