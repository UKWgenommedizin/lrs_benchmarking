# ************************************************************************************************
# GoldRush PacBio HiFi whole-genome assembly -- BENCHMARK-ONLY variant
#
# Same Docker image, thread count (32), memory ceiling (64 GB) and
# goldrush run parameters as pb.assembly.goldrush.smk, so the measured
# numbers describe the same execution conditions as the assembly already
# reported in Figure 4 / assembly_benchmark_30x.tsv -- this file does not
# change any parameter from the production rule, only what is measured
# and kept.
#
# See ont.assembly.goldrush.benchmark.smk's header comment for why
# GoldRush's own internal track_time=1 timing files are collected and
# aggregated here (sum of wall-clock, max of peak RSS across stages)
# instead of adding a new instrumentation mechanism, and why the
# regenerated assembly/intermediates are discarded rather than kept.
# ************************************************************************************************

import os

CWD = os.getcwd()
print("Current working directory: " + CWD)

try:
    DATASET_FILTER = config["dataset_filter"]
except (KeyError, NameError):
    DATASET_FILTER = None

#################
# GoldRush version (identical to pb.assembly.goldrush.smk)

GOLDRUSH_VERSION = "1.2.2-ntlinkfix"
DOCKER_GOLDRUSH = "nicolasardila1/lrs-goldrush:" + GOLDRUSH_VERSION

print("GoldRush version: " + GOLDRUSH_VERSION)

#####################
# Discover datasets and create wildcards (identical filter to pb.assembly.goldrush.smk)

DATASETS_FASTQ, = glob_wildcards(CWD + r"/fastq/{dataset,[A-Za-z0-9._-]+}.fastq.gz")

DATASETS = [
    dataset
    for dataset in DATASETS_FASTQ
    if ".pb." in dataset.lower()
    and ".1k" not in dataset.lower()
    and ".chr21." not in dataset.lower()
    and ".localtest." not in dataset.lower()
    and "smoke" not in dataset.lower()]


##############
# Targets

OUTPUT = []

OUTPUT += expand(
    CWD + "/assemblers_benchmark/goldrush/{dataset}.perf.tsv", zip, dataset=DATASETS
)

rule all:
    input:
        OUTPUT

print("Discover datasets and create wildcards")
print(OUTPUT)


################
# Prevent local test datasets from being requested explicitly

wildcard_constraints:
    dataset = r"(?=.*\.pb\.)(?!.*\.1k(?:\.|$))(?!.*\.chr21\.)(?!.*\.localtest\.)[A-Za-z0-9._-]+"


################
# GoldRush resource requirements (identical to pb.assembly.goldrush.smk)

def get_goldrush_memory(wildcards):
    return 64000

################
# Rules

rule goldrush_benchmark:
    input:
        fastq = CWD + "/fastq/{dataset}.fastq.gz"

    output:
        perf = CWD + "/assemblers_benchmark/goldrush/{dataset}.perf.tsv"

    params:
        outdir = CWD + "/assemblers_benchmark/goldrush/{dataset}.scratch",
        genome_size = "3e9",
        prefix = "{dataset}_goldrush",
        shm_size = "8g"

    log:
        CWD + "/assemblers_benchmark/goldrush/{dataset}.benchmark.log"

    message:
        "executing {rule} with output {output} and input {input}"

    threads: 32

    resources:
        mem_mb = get_goldrush_memory

    shell:
        """
        mkdir -p "$(dirname "{log}")"

        (
            set -euo pipefail

            echo "[$(date -Is)] START goldrush_benchmark {wildcards.dataset}"
            echo "Dataset: {wildcards.dataset}"
            echo "Read technology: PB"
            echo "Genome size: {params.genome_size}"
            echo "Threads: {threads}"
            echo "Memory: {resources.mem_mb} MB"
            echo "Shared memory: {params.shm_size}"
            echo "Container hostname: goldrush-benchmark-{wildcards.dataset}"

            mkdir -p "{params.outdir}"
            mkdir -p "$(dirname "{output.perf}")"

            # GoldRush requires an uncompressed FASTQ file.
            READS_FASTQ="{params.outdir}/{wildcards.dataset}.fastq"

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

            # Reset GoldRush internal workflow state (same rationale as the
            # production rule: stale/incomplete ntLink checkpoints from a
            # previous run must not be silently reused).
            INTERMEDIATE_DIR="{params.outdir}/goldrush_intermediate_files"

            if [[ -d "$INTERMEDIATE_DIR" ]]; then
                echo "Removing previous GoldRush internal state:"
                echo "$INTERMEDIATE_DIR"
                rm -rf -- "$INTERMEDIATE_DIR"
            fi

            docker run --rm \
                --tmpfs /tmp:size=50g,exec \
                --hostname goldrush-benchmark-{wildcards.dataset} \
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
                m=10000 \
                P=0 \
                p={params.prefix} \
                track_time=1

            FINAL_ASSEMBLY=$(find \
                "$INTERMEDIATE_DIR" \
                -type f \
                -name "*ntLink-5rounds.polished.fa" \
                -print \
                | sort \
                | tail -n 1)

            [[ -n "$FINAL_ASSEMBLY" && -s "$FINAL_ASSEMBLY" ]] || {{
                echo "ERROR: GoldRush final assembly is missing or empty -- run did not"
                echo "complete, so its partial *.time files cannot be trusted as a total."
                echo "Expected under: $INTERMEDIATE_DIR"
                exit 101;
            }}

            mapfile -t TIME_FILES < <(find "$INTERMEDIATE_DIR" -maxdepth 1 -type f -name "*.time" | sort)

            [[ ${{#TIME_FILES[@]}} -gt 0 ]] || {{
                echo "ERROR: no *.time files found under $INTERMEDIATE_DIR"
                echo "track_time=1 should have produced one per internal stage."
                exit 102;
            }}

            echo "Found ${{#TIME_FILES[@]}} stage timing file(s):"
            printf '  %s\\n' "${{TIME_FILES[@]}}"

            python3 "{CWD}/assembly_analysis/scripts/metrics/parse_assembler_time_v.py" \
                --assembler goldrush \
                --sample "{wildcards.dataset}" \
                --technology pb \
                --threads {threads} \
                --time-file "${{TIME_FILES[@]}}" \
                --output "{output.perf}"

            echo "Discarding regenerated assembly and intermediates (already have"
            echo "QUAST metrics for this assembly):"
            echo "  {params.outdir}"
            rm -rf "{params.outdir}"

            echo "[$(date -Is)] END goldrush_benchmark {wildcards.dataset}"

        ) > "{log}" 2>&1
        """
