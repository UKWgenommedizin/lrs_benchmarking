# ************************************************************************************************
# GoldRush ONT whole-genome assembly -- RUNTIME / PEAK-RAM MEASUREMENT ONLY
#
# Re-runs GoldRush with the same Docker image, memory ceiling and
# `goldrush run` parameters as ont.assembly.goldrush.smk, but at
# RAM_TIME_THREADS (default 64), only to record wall-clock time and peak
# RAM. Nothing else is kept.
#
# GoldRush self-instruments every internal Makefile stage via
# `command time -v -o <stage>.time` when track_time=1 is passed, so no extra
# wrapper is needed: the per-stage *.time files are summed (wall-clock) and
# max'd (peak RSS) by parse_assembler_time_v.py.
#
# Rules (per dataset):
#   goldrush_run      -> runs GoldRush in a temp() scratch directory
#   goldrush_record   -> parses the stage *.time files into {dataset}.ram_time.tsv
#                        (the only result this workflow keeps)
#
# The scratch directory is a temp() output, so Snakemake deletes it automatically
# (decompressed reads, assembly, intermediates, *.time files) as soon as the
# record rule succeeds, and also removes it if the run fails.
# Every rule writes its own log ({run,record}.log) per CONSTITUTION VI.3/VI.4.
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
RAM_TIME_DIR = os.path.join(CWD, config.get("ram_time_dir", "assemblers_ram_time"), "goldrush")
SCRATCH_DIR = os.path.join(CWD, config.get("ram_time_scratch_dir", "assemblers_ram_time/scratch"), "goldrush")

for _path in (RAM_TIME_DIR, SCRATCH_DIR):
    if os.path.commonpath([CWD, os.path.realpath(_path)]) != os.path.realpath(CWD):
        raise ValueError("Constitution I.1: ram_time paths must be inside the repository root: " + _path)

print("RAM/time threads: " + str(RAM_TIME_THREADS))
print("RAM/time results: " + RAM_TIME_DIR)
print("RAM/time scratch: " + SCRATCH_DIR)

#################
# GoldRush version (identical to ont.assembly.goldrush.smk)
GOLDRUSH_VERSION = "1.2.2-ntlinkfix"
DOCKER_GOLDRUSH = "nicolasardila1/lrs-goldrush:" + GOLDRUSH_VERSION

print("GoldRush version: " + GOLDRUSH_VERSION)

#####################
# Discover datasets and create wildcards (identical filter to ont.assembly.goldrush.smk)

DATASETS_FASTQ, = glob_wildcards(CWD + r"/fastq/{dataset,[A-Za-z0-9._-]+}.fastq.gz")

DATASETS = [
    dataset
    for dataset in DATASETS_FASTQ
    if ".ont." in dataset.lower()
    and ".1k" not in dataset.lower()
    and ".chr21." not in dataset.lower()
    and ".localtest." not in dataset.lower()
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
# Prevent local test datasets from being requested explicitly

wildcard_constraints:
    dataset = r"(?=.*\.ont\.)(?!.*\.1k(?:\.|$))(?!.*\.chr21\.)(?!.*\.localtest\.)[A-Za-z0-9._-]+"


################
# GoldRush resource requirements (identical to ont.assembly.goldrush.smk)

def get_goldrush_memory(wildcards):
    return 64000

################
# Rules

rule goldrush_run:
    input:
        fastq = CWD + "/fastq/{dataset}.fastq.gz"

    output:
        scratch = temp(directory(SCRATCH_DIR + "/{dataset}"))

    params:
        genome_size = "3e9",
        prefix = "{dataset}_goldrush",
        min_length = 5000,
        shm_size = "8g"

    log:
        RAM_TIME_DIR + "/{dataset}.run.log"

    message:
        "executing {rule} with output {output} and input {input}"

    threads: RAM_TIME_THREADS

    resources:
        mem_mb = get_goldrush_memory

    shell:
        """
        mkdir -p "$(dirname "{log}")"

        (
            set -euo pipefail

            echo "[$(date -Is)] START goldrush_run {wildcards.dataset}"
            echo "Dataset: {wildcards.dataset}"
            echo "Read technology: ONT"
            echo "Genome size: {params.genome_size}"
            echo "Threads: {threads}"
            echo "Memory: {resources.mem_mb} MB"
            echo "Shared memory: {params.shm_size}"
            echo "Scratch: {output.scratch}"

            [[ {threads} -eq {RAM_TIME_THREADS} ]] || {{
                echo "ERROR: Snakemake granted {threads} threads, expected {RAM_TIME_THREADS}."
                echo "Run with --cores {RAM_TIME_THREADS} (or more) so the measurement is comparable."
                exit 104;
            }}

            mkdir -p "{output.scratch}"

            echo "Checking gzip and python3 are available inside {DOCKER_GOLDRUSH} ..."
            docker run --rm --entrypoint sh "{DOCKER_GOLDRUSH}" \
                -c 'command -v gzip && command -v python3' || {{
                echo "ERROR: gzip or python3 (needed by the record rule) is not available inside {DOCKER_GOLDRUSH}"
                exit 103;
            }}

            # GoldRush requires an uncompressed FASTQ next to its working directory.
            echo "Preparing uncompressed GoldRush input..."
            docker run --rm \
                --tmpfs /tmp:size=50g,exec \
                --hostname goldrush-gunzip-{wildcards.dataset} \
                -u $UID:$(id -g) \
                -v "{CWD}:{CWD}" \
                --entrypoint sh \
                {DOCKER_GOLDRUSH} \
                -c 'gzip -cd "{input.fastq}" > "{output.scratch}/{wildcards.dataset}.fastq"'

            [[ -s "{output.scratch}/{wildcards.dataset}.fastq" ]] || {{
                echo "ERROR: decompressed GoldRush FASTQ is missing or empty"
                exit 102
            }}

            docker run --rm \
                --tmpfs /tmp:size=50g,exec \
                --hostname goldrush-ram-time-{wildcards.dataset} \
                --workdir "{output.scratch}" \
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
                m={params.min_length} \
                P=0 \
                p={params.prefix} \
                track_time=1

            FINAL_ASSEMBLY=$(find \
                "{output.scratch}/goldrush_intermediate_files" \
                -type f \
                -name "*ntLink-5rounds.polished.fa" \
                -print \
                | sort \
                | tail -n 1)

            [[ -n "$FINAL_ASSEMBLY" && -s "$FINAL_ASSEMBLY" ]] || {{
                echo "ERROR: GoldRush final assembly is missing or empty -- run did not"
                echo "complete, so its partial *.time files cannot be trusted as a total."
                exit 101;
            }}

            echo "[$(date -Is)] END goldrush_run {wildcards.dataset}"

        ) > "{log}" 2>&1
        """

rule goldrush_record:
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

            echo "[$(date -Is)] START goldrush_record {wildcards.dataset}"
            echo "Container hostname: goldrush-record-{wildcards.dataset}"

            mapfile -t TIME_FILES < <(find "{input.scratch}/goldrush_intermediate_files" -maxdepth 1 -type f -name "*.time" | sort)

            [[ ${{#TIME_FILES[@]}} -gt 0 ]] || {{
                echo "ERROR: no *.time files found under {input.scratch}/goldrush_intermediate_files"
                echo "track_time=1 should have produced one per internal stage."
                exit 102;
            }}

            echo "Found ${{#TIME_FILES[@]}} stage timing file(s)"

            # Parser is stdlib-only Python, run inside the assembler's own pinned image.
            docker run --rm \
                --tmpfs /tmp:size=50g,exec \
                --hostname goldrush-record-{wildcards.dataset} \
                -u $UID:$(id -g) \
                -v {CWD}:{CWD} \
                --entrypoint python3 \
                {DOCKER_GOLDRUSH} \
                "{CWD}/assembly_analysis/scripts/metrics/parse_assembler_time_v.py" \
                --assembler goldrush \
                --sample "{wildcards.dataset}" \
                --technology ont \
                --threads {RAM_TIME_THREADS} \
                --time-file "${{TIME_FILES[@]}}" \
                --output "{output.ram_time}"

            [[ $(wc -l < "{output.ram_time}") -eq 2 ]] || {{
                echo "ERROR: {output.ram_time} must contain exactly one header and one data row"
                exit 101;
            }}

            echo "[$(date -Is)] END goldrush_record {wildcards.dataset}"

        ) > "{log}" 2>&1
        """
