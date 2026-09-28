# ************************************************************************************************
# Verkko hybrid whole-genome assembly -- BENCHMARK-ONLY variant
#
# Same Docker image, thread count (32) and memory ceiling (72 GB) as
# hybrid.assembly.verkko.smk, so the measured numbers describe the same
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

OUTPUT += expand(CWD + "/assemblers_benchmark/verkko/{sample}.perf.tsv", sample=SAMPLES)

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

rule verkko_benchmark:
    input:
        hifi = CWD + "/fastq/{sample}.pb.30x.fastq.gz",
        ont  = CWD + "/fastq/{sample}.ont.30x.fastq.gz"

    output:
        perf = CWD + "/assemblers_benchmark/verkko/{sample}.perf.tsv"

    params:
        outdir = CWD + "/assemblers_benchmark/verkko/{sample}.scratch/work",
        scratch = CWD + "/assemblers_benchmark/verkko/{sample}.scratch",
        time_file = CWD + "/assemblers_benchmark/verkko/{sample}.time",
        local_memory_gb = 64

    log:
        CWD + "/assemblers_benchmark/verkko/{sample}.benchmark.log"

    message:
        "executing {rule} with output {output} and input {input}"

    threads: 32

    resources:
        mem_mb = get_verkko_memory

    shell:
        """
        mkdir -p "$(dirname "{log}")"

        (
            set -euo pipefail

            echo "[$(date -Is)] START verkko_benchmark {wildcards.sample}"
            echo "Sample: {wildcards.sample}"
            echo "Read technologies: PacBio HiFi + ONT"
            echo "HiFi input: {input.hifi}"
            echo "ONT input: {input.ont}"
            echo "Threads: {threads}"
            echo "Memory: {resources.mem_mb} MB"
            echo "Container hostname: verkko-benchmark-{wildcards.sample}"

            mkdir -p "{params.outdir}"
            mkdir -p "$(dirname "{output.perf}")"

            echo "Checking /usr/bin/time -v is available inside {DOCKER_VERKKO} ..."
            docker run --rm --entrypoint sh "{DOCKER_VERKKO}" \
                -c 'command -v /usr/bin/time' || {{
                echo "ERROR: /usr/bin/time is not available inside {DOCKER_VERKKO}"
                echo "Cannot measure peak RAM for this container without it -- aborting"
                echo "before running the multi-hour assembly, rather than after."
                exit 103;
            }}

            docker run --rm \
                --hostname verkko-benchmark-{wildcards.sample} \
                --cpus {threads} \
                -m {resources.mem_mb}m \
                --tmpfs /tmp:size=50g,exec \
                --user "$(id -u):$(id -g)" \
                -e HOME=/tmp \
                -e TMPDIR=/tmp \
                --workdir {CWD} \
                -v {CWD}:{CWD} \
                --entrypoint /usr/bin/time \
                {DOCKER_VERKKO} \
                -v -o "{params.time_file}" \
                verkko \
                    -d "{params.outdir}" \
                    --hifi "{input.hifi}" \
                    --nano "{input.ont}" \
                    --local \
                    --local-cpus {threads} \
                    --local-memory {params.local_memory_gb}

            [[ -s "{params.outdir}/assembly.fasta" ]] || {{
                echo "ERROR: Verkko assembly.fasta is missing or empty"
                echo "Expected: {params.outdir}/assembly.fasta"
                exit 101;
            }}

            grep -q '^>' "{params.outdir}/assembly.fasta" || {{
                echo "ERROR: Verkko output does not appear to be a valid FASTA"
                exit 101;
            }}

            [[ -s "{params.time_file}" ]] || {{
                echo "ERROR: /usr/bin/time output ({params.time_file}) is missing or empty"
                exit 102;
            }}

            python3 "{CWD}/assembly_analysis/scripts/metrics/parse_assembler_time_v.py" \
                --assembler verkko \
                --sample "{wildcards.sample}" \
                --technology hybrid \
                --threads {threads} \
                --time-file "{params.time_file}" \
                --output "{output.perf}"

            echo "Discarding regenerated assembly (already have QUAST metrics for it):"
            echo "  {params.scratch}"
            rm -rf "{params.scratch}"
            rm -f "{params.time_file}"

            echo "[$(date -Is)] END verkko_benchmark {wildcards.sample}"

        ) > "{log}" 2>&1
        """
