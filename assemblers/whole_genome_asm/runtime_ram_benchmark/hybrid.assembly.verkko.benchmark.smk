# ************************************************************************************************
# Verkko hybrid whole-genome assembly -- RUNTIME / PEAK-RAM MEASUREMENT ONLY
#
# Re-runs Verkko with the same Docker image, memory ceiling and parameters
# as hybrid.assembly.verkko.smk, but at BENCHMARK_THREADS (default 64), only
# to record wall-clock time and peak RAM. Nothing else is kept.
#
# Rules (per sample):
#   verkko_run      -> runs Verkko under /usr/bin/time -v inside the container,
#                      writing everything to a temp() scratch directory
#   verkko_record   -> parses the time file into {sample}.perf.tsv
#                      (the only result this workflow keeps)
#
# The scratch directory is a temp() output, so Snakemake deletes it automatically
# as soon as the record rule succeeds (and also removes it if the run fails).
# Every rule writes its own log ({run,record}.log) per CONSTITUTION VI.3/VI.4.
#
# See ont.assembly.flye2.benchmark.smk for why /usr/bin/time -v runs inside
# the container instead of using Snakemake's `benchmark:` directive.
#
# Configurable (--config key=value), no paths hard-coded in the rules;
# both directories must stay inside the repository root (CONSTITUTION I.1):
#   benchmark_threads      (default 64)
#   benchmark_dir          (default assemblers_benchmark)          -- perf TSVs + logs
#   benchmark_scratch_dir  (default assemblers_benchmark/scratch)  -- deleted automatically after each run
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
# Benchmark settings

BENCHMARK_THREADS = int(config.get("benchmark_threads", 64))
BENCHMARK_DIR = os.path.join(CWD, config.get("benchmark_dir", "assemblers_benchmark"), "verkko")
SCRATCH_DIR = os.path.join(CWD, config.get("benchmark_scratch_dir", "assemblers_benchmark/scratch"), "verkko")

for _path in (BENCHMARK_DIR, SCRATCH_DIR):
    if os.path.commonpath([CWD, os.path.realpath(_path)]) != os.path.realpath(CWD):
        raise ValueError("Constitution I.1: benchmark paths must be inside the repository root: " + _path)

print("Benchmark threads: " + str(BENCHMARK_THREADS))
print("Benchmark results: " + BENCHMARK_DIR)
print("Benchmark scratch: " + SCRATCH_DIR)

#################
# Verkko version (identical to hybrid.assembly.verkko.smk)
VERKKO_VERSION = "2.3.2"
DOCKER_VERKKO = "nicolasardila1/lrs-verkko2:" + VERKKO_VERSION

print("Verkko version: " + VERKKO_VERSION)

#####################
# Discover samples and create wildcards (identical to hybrid.assembly.verkko.smk)

ONT_SAMPLES, = glob_wildcards(CWD + r"/fastq/{sample,[A-Za-z0-9_-]+}.ont.30x.fastq.gz")

PB_SAMPLES, = glob_wildcards(CWD + r"/fastq/{sample,[A-Za-z0-9_-]+}.pb.30x.fastq.gz")

TEST_SAMPLE_MARKERS = ("SMOKE", "LOCALTEST")

def is_production_sample(sample):
    upper = sample.upper()
    return not any(marker in upper for marker in TEST_SAMPLE_MARKERS)

ONT_SAMPLE_SET = {sample for sample in ONT_SAMPLES if is_production_sample(sample)}
PB_SAMPLE_SET = {sample for sample in PB_SAMPLES if is_production_sample(sample)}

SAMPLES = sorted(ONT_SAMPLE_SET & PB_SAMPLE_SET)

######################
# Input samples and unpaired control checkpoint
MISSING_PB = sorted(ONT_SAMPLE_SET - PB_SAMPLE_SET)
MISSING_ONT = sorted(PB_SAMPLE_SET - ONT_SAMPLE_SET)

if MISSING_PB:
    raise ValueError("Missing PacBio HiFi input for samples: " + ", ".join(MISSING_PB))

if MISSING_ONT:
    raise ValueError("Missing ONT input for samples: " + ", ".join(MISSING_ONT))

if not SAMPLES:
    raise ValueError("No paired Verkko WGS inputs were found. " "Expected files such as "
    f"{CWD}/fastq/HG002.pb.30x.fastq.gz and " f"{CWD}/fastq/HG002.ont.30x.fastq.gz")

##############
# Targets

OUTPUT = []

OUTPUT += expand(BENCHMARK_DIR + "/{sample}.perf.tsv", sample=SAMPLES)

rule all:
    input:
        OUTPUT

print("Discover samples and create wildcards")
print(OUTPUT)


################
# Prevent local test sample from being requested explicitly

wildcard_constraints:
    sample = r"[A-Za-z0-9_-]+"


################
# Verkko resource requirements (identical to hybrid.assembly.verkko.smk)

def get_verkko_memory(wildcards):
    return 72000

################
# Rules

rule verkko_run:
    input:
        hifi = CWD + "/fastq/{sample}.pb.30x.fastq.gz",
        ont  = CWD + "/fastq/{sample}.ont.30x.fastq.gz"

    output:
        scratch = temp(directory(SCRATCH_DIR + "/{sample}"))

    params:
        local_memory_gb = 64

    log:
        BENCHMARK_DIR + "/{sample}.run.log"

    message:
        "executing {rule} with output {output} and input {input}"

    threads: BENCHMARK_THREADS

    resources:
        mem_mb = get_verkko_memory

    shell:
        """
        mkdir -p "$(dirname "{log}")"

        (
            set -euo pipefail

            echo "[$(date -Is)] START verkko_run {wildcards.sample}"
            echo "Sample: {wildcards.sample}"
            echo "Read technologies: PacBio HiFi + ONT"
            echo "HiFi input: {input.hifi}"
            echo "ONT input: {input.ont}"
            echo "Threads: {threads}"
            echo "Memory: {resources.mem_mb} MB"
            echo "Scratch: {output.scratch}"

            [[ {threads} -eq {BENCHMARK_THREADS} ]] || {{
                echo "ERROR: Snakemake granted {threads} threads, expected {BENCHMARK_THREADS}."
                echo "Run with --cores {BENCHMARK_THREADS} (or more) so the measurement is comparable."
                exit 104;
            }}

            mkdir -p "{output.scratch}"

            echo "Checking /usr/bin/time -v is available inside {DOCKER_VERKKO} ..."
            docker run --rm --entrypoint sh "{DOCKER_VERKKO}" \
                -c 'command -v /usr/bin/time && command -v python3' || {{
                echo "ERROR: /usr/bin/time or python3 (needed by the record rule) is not available inside {DOCKER_VERKKO}"
                exit 103;
            }}

            docker run --rm \
                --hostname verkko-benchmark-{wildcards.sample} \
                --cpus {threads} \
                -m {resources.mem_mb}m \
                --tmpfs /tmp:size=50g,exec \
                -u $UID:$(id -g) \
                -e HOME=/tmp \
                -e TMPDIR=/tmp \
                --workdir {CWD} \
                -v {CWD}:{CWD} \
                --entrypoint /usr/bin/time \
                {DOCKER_VERKKO} \
                -v -o "{output.scratch}/time_v.txt" \
                verkko \
                    -d "{output.scratch}/work" \
                    --hifi "{input.hifi}" \
                    --nano "{input.ont}" \
                    --local \
                    --local-cpus {threads} \
                    --local-memory {params.local_memory_gb}

            [[ -s "{output.scratch}/work/assembly.fasta" ]] && grep -q '^>' "{output.scratch}/work/assembly.fasta" || {{
                echo "ERROR: Verkko assembly.fasta is missing, empty or not FASTA -- run did not complete"
                exit 101;
            }}

            [[ -s "{output.scratch}/time_v.txt" ]] || {{
                echo "ERROR: /usr/bin/time output is missing or empty"
                exit 102;
            }}

            echo "[$(date -Is)] END verkko_run {wildcards.sample}"

        ) > "{log}" 2>&1
        """

rule verkko_record:
    input:
        scratch = SCRATCH_DIR + "/{sample}"

    output:
        perf = BENCHMARK_DIR + "/{sample}.perf.tsv"

    log:
        BENCHMARK_DIR + "/{sample}.record.log"

    message:
        "executing {rule} with output {output} and input {input}"

    shell:
        """
        (
            set -eo pipefail

            echo "[$(date -Is)] START verkko_record {wildcards.sample}"
            echo "Container hostname: verkko-record-{wildcards.sample}"

            # Parser is stdlib-only Python, run inside the assembler's own pinned image.
            docker run --rm \
                --tmpfs /tmp:size=50g,exec \
                --hostname verkko-record-{wildcards.sample} \
                -u $UID:$(id -g) \
                -v {CWD}:{CWD} \
                --entrypoint python3 \
                {DOCKER_VERKKO} \
                "{CWD}/assembly_analysis/scripts/metrics/parse_assembler_time_v.py" \
                --assembler verkko \
                --sample "{wildcards.sample}" \
                --technology hybrid \
                --threads {BENCHMARK_THREADS} \
                --time-file "{input.scratch}/time_v.txt" \
                --output "{output.perf}"

            [[ $(wc -l < "{output.perf}") -eq 2 ]] || {{
                echo "ERROR: {output.perf} must contain exactly one header and one data row"
                exit 101;
            }}

            echo "[$(date -Is)] END verkko_record {wildcards.sample}"

        ) > "{log}" 2>&1
        """
