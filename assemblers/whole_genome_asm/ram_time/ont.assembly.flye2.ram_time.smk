# ************************************************************************************************
# Flye ONT whole-genome assembly -- RUNTIME / PEAK-RAM MEASUREMENT ONLY
#
# Re-runs Flye with the same Docker image and parameters as
# ont.assembly.flye2.smk, but at RAM_TIME_THREADS (default 64), only to
# record wall-clock time and peak RAM. Nothing else is kept.
#
# Rules (per dataset):
#   flye_run      -> runs Flye under /usr/bin/time -v inside the container,
#                    writing everything to a temp() scratch directory
#   flye_record   -> parses the time file into {dataset}.ram_time.tsv
#                    (the only result this workflow keeps)
#
# The scratch directory is a temp() output, so Snakemake deletes it automatically
# as soon as the record rule succeeds (and also removes it if the run fails).
# Every rule writes its own log ({run,record}.log) per CONSTITUTION VI.3/VI.4.
#
# /usr/bin/time -v runs *inside* the container because Snakemake's
# `benchmark:` directive samples the host process tree and cannot see the
# container's real memory use.
#
# Configurable (--config key=value), no paths hard-coded in the rules;
# both directories must stay inside the repository root (CONSTITUTION I.1):
#   ram_time_threads      (default 64)
#   ram_time_dir          (default assemblers_ram_time)          -- ram_time TSVs + logs
#   ram_time_scratch_dir  (default assemblers_ram_time/scratch)  -- deleted automatically after each run
# Relative paths are resolved against the repository root (the working directory).
# ************************************************************************************************

import os

CWD = os.getcwd()
print("Current working directory: " + CWD)

try:
    DATASET_FILTER = config["dataset_filter"]
except (KeyError, NameError):
    DATASET_FILTER = None

#################
# RAM / time measurement settings

RAM_TIME_THREADS = int(config.get("ram_time_threads", 64))
RAM_TIME_DIR = os.path.join(CWD, config.get("ram_time_dir", "assemblers_ram_time"), "flye")
SCRATCH_DIR = os.path.join(CWD, config.get("ram_time_scratch_dir", "assemblers_ram_time/scratch"), "flye")

for _path in (RAM_TIME_DIR, SCRATCH_DIR):
    if os.path.commonpath([CWD, os.path.realpath(_path)]) != os.path.realpath(CWD):
        raise ValueError("Constitution I.1: ram_time paths must be inside the repository root: " + _path)

print("RAM/time threads: " + str(RAM_TIME_THREADS))
print("RAM/time results: " + RAM_TIME_DIR)
print("RAM/time scratch: " + SCRATCH_DIR)

#################
# Flye version (identical to ont.assembly.flye2.smk)
FLYE_VERSION = "2.9.6"
DOCKER_FLYE = "nicolasardila1/lrs-flye2:" + FLYE_VERSION

print("Flye version: " + FLYE_VERSION)

#####################
# Discover datasets and create wildcards (identical filter to ont.assembly.flye2.smk)
DATASETS_FASTQ, = glob_wildcards(CWD + r"/fastq/{dataset,[A-Za-z0-9._-]+}.fastq.gz")
DATASETS = [
    dataset
    for dataset in DATASETS_FASTQ
    if ".ont." in dataset.lower()
    and ".1k" not in dataset.lower()
    and ".chr21." not in dataset.lower()
    and "localtest" not in dataset.lower()
    and "smoke" not in dataset.lower()]

##############
# Targets

OUTPUT = []

OUTPUT += expand(RAM_TIME_DIR + "/{dataset}.ram_time.tsv", dataset=DATASETS)

rule all:
    input:
        OUTPUT

print("Discover datasets and create wildcards")
print(OUTPUT)

################
# Flye resource requirements (identical to ont.assembly.flye2.smk)

def get_flye_memory(wildcards):
    dataset = wildcards.dataset.lower()

    if ".1k" in dataset:
        return 8

    if ".chr21." in dataset:
        return 8

    return 240

################
# Rules

rule flye_run:
    input:
        fastq = CWD + "/fastq/{dataset}.fastq.gz"

    output:
        scratch = temp(directory(SCRATCH_DIR + "/{dataset}"))

    log:
        RAM_TIME_DIR + "/{dataset}.run.log"

    message:
        "executing {rule} with output {output} and input {input}"

    threads: RAM_TIME_THREADS

    resources:
        mem_gb = get_flye_memory

    shell:
        """
        mkdir -p "$(dirname "{log}")"

        (
            set -eo pipefail

            echo "[$(date -Is)] START flye_run {wildcards.dataset}"
            echo "Dataset: {wildcards.dataset}"
            echo "Read technology: ONT"
            echo "Read mode: --nano-hq"
            echo "Threads: {threads}"
            echo "Memory: {resources.mem_gb} GB"
            echo "Scratch: {output.scratch}"

            [[ {threads} -eq {RAM_TIME_THREADS} ]] || {{
                echo "ERROR: Snakemake granted {threads} threads, expected {RAM_TIME_THREADS}."
                echo "Run with --cores {RAM_TIME_THREADS} (or more) so the measurement is comparable."
                exit 104;
            }}

            mkdir -p "{output.scratch}"

            echo "Checking /usr/bin/time -v is available inside {DOCKER_FLYE} ..."
            docker run --rm --entrypoint sh "{DOCKER_FLYE}" \
                -c 'command -v /usr/bin/time && command -v python3' || {{
                echo "ERROR: /usr/bin/time or python3 (needed by the record rule) is not available inside {DOCKER_FLYE}"
                exit 103;
            }}

            docker run --rm \
                --tmpfs /tmp:size=50g,exec \
                --hostname flye-ram-time-{wildcards.dataset} \
                --workdir /tmp \
                -u $UID:$(id -g) \
                --cpus {threads} \
                -m {resources.mem_gb}g \
                -v {CWD}:{CWD} \
                --entrypoint /usr/bin/time \
                {DOCKER_FLYE} \
                -v -o "{output.scratch}/time_v.txt" \
                flye \
                --nano-hq \
                {input.fastq} \
                --out-dir {output.scratch} \
                --threads {threads}

            [[ -s "{output.scratch}/assembly.fasta" ]] && grep -q '^>' "{output.scratch}/assembly.fasta" || {{
                echo "ERROR: Flye assembly.fasta is missing, empty or not FASTA -- run did not complete"
                exit 101;
            }}

            [[ -s "{output.scratch}/time_v.txt" ]] || {{
                echo "ERROR: /usr/bin/time output is missing or empty"
                exit 102;
            }}

            echo "[$(date -Is)] END flye_run {wildcards.dataset}"

        ) > "{log}" 2>&1
        """

rule flye_record:
    input:
        scratch = SCRATCH_DIR + "/{dataset}"

    output:
        ram_time = RAM_TIME_DIR + "/{dataset}.ram_time.tsv"

    log:
        RAM_TIME_DIR + "/{dataset}.record.log"

    message:
        "executing {rule} with output {output} and input {input}"

    shell:
        """
        (
            set -eo pipefail

            echo "[$(date -Is)] START flye_record {wildcards.dataset}"
            echo "Container hostname: flye-record-{wildcards.dataset}"

            # Parser is stdlib-only Python, run inside the assembler's own pinned image.
            docker run --rm \
                --tmpfs /tmp:size=50g,exec \
                --hostname flye-record-{wildcards.dataset} \
                -u $UID:$(id -g) \
                -v {CWD}:{CWD} \
                --entrypoint python3 \
                {DOCKER_FLYE} \
                "{CWD}/assembly_analysis/scripts/metrics/parse_assembler_time_v.py" \
                --assembler flye \
                --sample "{wildcards.dataset}" \
                --technology ont \
                --threads {RAM_TIME_THREADS} \
                --time-file "{input.scratch}/time_v.txt" \
                --output "{output.ram_time}"

            [[ $(wc -l < "{output.ram_time}") -eq 2 ]] || {{
                echo "ERROR: {output.ram_time} must contain exactly one header and one data row"
                exit 101;
            }}

            echo "[$(date -Is)] END flye_record {wildcards.dataset}"

        ) > "{log}" 2>&1
        """
