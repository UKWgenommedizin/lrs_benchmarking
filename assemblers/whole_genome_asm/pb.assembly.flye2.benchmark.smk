# ************************************************************************************************
# Flye PacBio HiFi whole-genome assembly -- BENCHMARK-ONLY variant
#
# Same Docker image, thread count (32) and memory ceiling as
# pb.assembly.flye2.smk, so the measured numbers describe the same
# execution conditions as the assembly already reported in Figure 4 /
# assembly_benchmark_30x.tsv -- this file does not change any parameter
# from the production rule, only what is measured and kept.
#
# See ont.assembly.flye2.benchmark.smk's header comment for why
# /usr/bin/time -v inside the container is used instead of Snakemake's
# own `benchmark:` directive, and why the regenerated assembly.fasta is
# discarded rather than kept.
# ************************************************************************************************

import os

CWD = os.getcwd()
print("Current working directory: " + CWD)

try:
    DATASET_FILTER = config["dataset_filter"]
except (KeyError, NameError):
    DATASET_FILTER = None

#################
# Flye version (identical to pb.assembly.flye2.smk)
FLYE_VERSION = "2.9.6"
DOCKER_FLYE = "nicolasardila1/lrs-flye2:" + FLYE_VERSION

print("Flye version: " + FLYE_VERSION)

#####################
# Discover datasets and create wildcards (identical filter to pb.assembly.flye2.smk)
DATASETS_FASTQ, = glob_wildcards(CWD + r"/fastq/{dataset,[A-Za-z0-9._-]+}.fastq.gz")
DATASETS = [
    dataset
    for dataset in DATASETS_FASTQ
    if ".pb." in dataset.lower()
    and ".1k" not in dataset.lower()
    and ".chr21." not in dataset.lower()
    and "localtest" not in dataset.lower()
    and "smoke" not in dataset.lower()]

##############
# Targets

OUTPUT = []

OUTPUT = OUTPUT + expand(
    CWD + "/assemblers_benchmark/flye/{dataset}.perf.tsv", zip, dataset=DATASETS
)

rule all:
    input:
        OUTPUT

print("Discover datasets and create wildcards")
print(OUTPUT)

################
# Flye resource requirements (identical to pb.assembly.flye2.smk)

def get_flye_memory(wildcards):
    dataset = wildcards.dataset.lower()

    if ".1k" in dataset:
        return 8

    if ".chr21." in dataset:
        return 8

    return 240

################
# Rules

rule flye_benchmark:
    input:
        fastq = CWD + "/fastq/{dataset}.fastq.gz"

    output:
        perf = CWD + "/assemblers_benchmark/flye/{dataset}.perf.tsv"

    params:
        scratch = CWD + "/assemblers_benchmark/flye/{dataset}.scratch",
        time_file = CWD + "/assemblers_benchmark/flye/{dataset}.time"

    log:
        CWD + "/assemblers_benchmark/flye/{dataset}.benchmark.log"

    message:
        "executing {rule} with output {output} and input {input}"

    threads: 32

    resources:
        mem_gb = get_flye_memory

    shell:
        """
        mkdir -p "$(dirname "{log}")"

        (
            set -eo pipefail

            echo "[$(date -Is)] START flye_benchmark {wildcards.dataset}"
            echo "Dataset: {wildcards.dataset}"
            echo "Read technology: PacBio HiFi"
            echo "Read mode: --pacbio-hifi"
            echo "Threads: {threads}"
            echo "Memory: {resources.mem_gb} GB"
            echo "Container hostname: flye-benchmark-{wildcards.dataset}"

            mkdir -p "{params.scratch}"
            mkdir -p "$(dirname "{output.perf}")"

            echo "Checking /usr/bin/time -v is available inside {DOCKER_FLYE} ..."
            docker run --rm --entrypoint sh "{DOCKER_FLYE}" \
                -c 'command -v /usr/bin/time' || {{
                echo "ERROR: /usr/bin/time is not available inside {DOCKER_FLYE}"
                echo "Cannot measure peak RAM for this container without it -- aborting"
                echo "before running the multi-hour assembly, rather than after."
                exit 103;
            }}

            docker run --rm \
                --tmpfs /tmp:size=50g,exec \
                --hostname flye-benchmark-{wildcards.dataset} \
                --workdir /tmp \
                -u $UID:$(id -g) \
                --cpus {threads} \
                -m {resources.mem_gb}g \
                -v {CWD}:{CWD} \
                --entrypoint /usr/bin/time \
                {DOCKER_FLYE} \
                -v -o "{params.time_file}" \
                flye \
                --pacbio-hifi \
                {input.fastq} \
                --out-dir {params.scratch} \
                --threads {threads}

            [[ -s "{params.scratch}/assembly.fasta" ]] || {{
                echo "ERROR: Flye assembly.fasta is missing or empty"
                exit 101;
            }}

            grep -q '^>' "{params.scratch}/assembly.fasta" || {{
                echo "ERROR: Flye output does not appear to be a valid FASTA"
                exit 101;
            }}

            [[ -s "{params.time_file}" ]] || {{
                echo "ERROR: /usr/bin/time output ({params.time_file}) is missing or empty"
                exit 102;
            }}

            python3 "{CWD}/assembly_analysis/scripts/metrics/parse_assembler_time_v.py" \
                --assembler flye \
                --sample "{wildcards.dataset}" \
                --technology pb \
                --threads {threads} \
                --time-file "{params.time_file}" \
                --output "{output.perf}"

            echo "Discarding regenerated assembly (already have QUAST metrics for it):"
            echo "  {params.scratch}"
            rm -rf "{params.scratch}"
            rm -f "{params.time_file}"

            echo "[$(date -Is)] END flye_benchmark {wildcards.dataset}"

        ) > "{log}" 2>&1
        """
