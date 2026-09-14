##
# ont.read_mapping.minimap2.smk
# Read mapping workflow for Oxford Nanopore data using minimap2.
# mapper_tag: mm2-ont
# Constitution: Articles I-VIII
##

include: "header_mapper.smk"

import os

####################
# Containers

MINIMAP2_VERSION = "2.28"
SAMTOOLS_VERSION = "1.24"

DOCKER_MINIMAP2 = "quay.io/biocontainers/minimap2:2.28--h577a1d6_4"
DOCKER_SAMTOOLS = "quay.io/biocontainers/samtools:1.24--h9dcdb79_1"

print("minimap2 version: " + MINIMAP2_VERSION)
print("minimap2 Docker image: " + DOCKER_MINIMAP2)
print("samtools version: " + SAMTOOLS_VERSION)
print("samtools Docker image: " + DOCKER_SAMTOOLS)

####################
# Reference

REFERENCE = "hg38"
MAPPER_TAG = "mm2-ont"

LOCAL_REFERENCE = os.path.join(
    CWD,
    "reference",
    "GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta",
)

RAW_REFERENCE = config.get("reference", LOCAL_REFERENCE)
REF = os.path.expanduser(RAW_REFERENCE)
if not os.path.isabs(REF):
    REF = os.path.join(CWD, REF)
REF = os.path.abspath(REF)
REFERENCE_DIR = os.path.dirname(REF)

print("Reference genome: " + REF)

####################
# Discover FASTQ inputs

FASTQ_DIR = "fastq"
DATASETS, = glob_wildcards(FASTQ_DIR + "/{dataset}.fastq.gz")
DATASETS = [d for d in DATASETS if ".ont." in d.lower()]

if DATASET_FILTER:
    DATASETS = [d for d in DATASETS if DATASET_FILTER in d]

####################
# Targets

rule all:
    input:
        expand(
            "cram/{dataset}.{ref}.{tag}.cram",
            dataset=DATASETS,
            ref=REFERENCE,
            tag=MAPPER_TAG,
        ),
        expand(
            "cram/{dataset}.{ref}.{tag}.cram.crai",
            dataset=DATASETS,
            ref=REFERENCE,
            tag=MAPPER_TAG,
        ),
        expand(
            "cram/{dataset}.{ref}.{tag}.cram.idxstats",
            dataset=DATASETS,
            ref=REFERENCE,
            tag=MAPPER_TAG,
        ),
        expand(
            "cram/{dataset}.{ref}.{tag}.cram.stats",
            dataset=DATASETS,
            ref=REFERENCE,
            tag=MAPPER_TAG,
        ),

####################
# Mapping

rule minimap2_ont_map_sort:
    input:
        fastq=FASTQ_DIR + "/{dataset}.fastq.gz",
        ref=lambda wildcards: REF,
    output:
        cram="cram/{dataset}." + REFERENCE + "." + MAPPER_TAG + ".cram",
        crai="cram/{dataset}." + REFERENCE + "." + MAPPER_TAG + ".cram.crai",
    log:
        "cram/{dataset}." + REFERENCE + "." + MAPPER_TAG + ".map_sort.log",
    threads: 64
    resources:
        mem_mb=131072,
    params:
        map_threads=48,
        sort_threads=16,
        map_mem_gb=64,
        sort_mem_gb=64,
    shell:
        r"""
        (
            set -eo pipefail

            echo "[$(date -Is)] START minimap2_ont_map_sort {wildcards.dataset}"
            echo "Dataset: {wildcards.dataset}"
            echo "Reference: {input.ref}"
            echo "Mapper: minimap2 {MINIMAP2_VERSION}"
            echo "Preset: map-ont"
            echo "Mapper container hostname: mm2-ont-{wildcards.dataset}"
            echo "Sort container hostname: samtools-sort-{wildcards.dataset}"
            echo "Threads: {threads}"
            echo "Memory resource: {resources.mem_mb} MB"

            mkdir -p "{CWD}/cram/tmp"

            docker run --rm \
                --hostname "mm2-ont-{wildcards.dataset}" \
                --tmpfs /tmp:size=50g,exec \
                -u $UID:$(id -g) \
                --cpus {params.map_threads} \
                -m {params.map_mem_gb}g \
                --workdir "{CWD}" \
                -v "{CWD}:{CWD}" \
                -v "{REFERENCE_DIR}:{REFERENCE_DIR}:ro" \
                --entrypoint minimap2 \
                "{DOCKER_MINIMAP2}" \
                -t {params.map_threads} \
                -R '@RG\tID:{wildcards.dataset}\tSM:{wildcards.dataset}' \
                -a \
                -y \
                -x map-ont \
                -L \
                --cs \
                --MD \
                "{input.ref}" \
                "{CWD}/{input.fastq}" \
            | docker run --rm -i \
                --hostname "samtools-sort-{wildcards.dataset}" \
                --tmpfs /tmp:size=50g,exec \
                -u $UID:$(id -g) \
                --cpus {params.sort_threads} \
                -m {params.sort_mem_gb}g \
                --workdir "{CWD}" \
                -v "{CWD}:{CWD}" \
                -v "{REFERENCE_DIR}:{REFERENCE_DIR}:ro" \
                --entrypoint samtools \
                "{DOCKER_SAMTOOLS}" \
                sort \
                -@ {params.sort_threads} \
                -m 3G \
                --reference "{input.ref}" \
                --no-PG \
                -O CRAM \
                -T "{CWD}/cram/tmp/{wildcards.dataset}.{REFERENCE}.{MAPPER_TAG}" \
                -o "{CWD}/{output.cram}" \
                -

            if [[ $(du -b "{output.cram}" | cut -f 1) -le 64 ]]; then
                echo "ERROR: CRAM is missing or truncated"
                exit 101
            fi

            docker run --rm \
                --hostname "samtools-index-{wildcards.dataset}" \
                --tmpfs /tmp:size=50g,exec \
                -u $UID:$(id -g) \
                --cpus 8 \
                -m 16g \
                --workdir "{CWD}" \
                -v "{CWD}:{CWD}" \
                -v "{REFERENCE_DIR}:{REFERENCE_DIR}:ro" \
                --entrypoint samtools \
                "{DOCKER_SAMTOOLS}" \
                index \
                -@ 8 \
                -o "{CWD}/{output.crai}" \
                "{CWD}/{output.cram}"

            if [[ ! -s "{output.crai}" ]]; then
                echo "ERROR: CRAI is missing or empty"
                exit 101
            fi

            echo "[$(date -Is)] END minimap2_ont_map_sort {wildcards.dataset}"
        ) > "{log}" 2>&1
        """

####################
# idxstats

rule minimap2_ont_idxstats:
    input:
        cram="cram/{dataset}." + REFERENCE + "." + MAPPER_TAG + ".cram",
        crai="cram/{dataset}." + REFERENCE + "." + MAPPER_TAG + ".cram.crai",
        ref=lambda wildcards: REF,
    output:
        idxstats="cram/{dataset}." + REFERENCE + "." + MAPPER_TAG + ".cram.idxstats",
    log:
        "cram/{dataset}." + REFERENCE + "." + MAPPER_TAG + ".idxstats.log",
    threads: 2
    resources:
        mem_mb=4096,
    shell:
        r"""
        (
            set -eo pipefail

            echo "[$(date -Is)] START minimap2_ont_idxstats {wildcards.dataset}"
            echo "Container hostname: samtools-idxstats-{wildcards.dataset}"

            docker run --rm \
                --hostname "samtools-idxstats-{wildcards.dataset}" \
                --tmpfs /tmp:size=50g,exec \
                -u $UID:$(id -g) \
                --cpus {threads} \
                -m 4g \
                --workdir "{CWD}" \
                -v "{CWD}:{CWD}" \
                -v "{REFERENCE_DIR}:{REFERENCE_DIR}:ro" \
                --entrypoint samtools \
                "{DOCKER_SAMTOOLS}" \
                idxstats \
                "{CWD}/{input.cram}" \
                > "{output.idxstats}"

            if [[ ! -s "{output.idxstats}" ]]; then
                echo "ERROR: idxstats output is missing or empty"
                exit 101
            fi

            echo "[$(date -Is)] END minimap2_ont_idxstats {wildcards.dataset}"
        ) > "{log}" 2>&1
        """

####################
# samtools stats

rule minimap2_ont_stats:
    input:
        cram="cram/{dataset}." + REFERENCE + "." + MAPPER_TAG + ".cram",
        crai="cram/{dataset}." + REFERENCE + "." + MAPPER_TAG + ".cram.crai",
        ref=lambda wildcards: REF,
    output:
        stats="cram/{dataset}." + REFERENCE + "." + MAPPER_TAG + ".cram.stats",
    log:
        "cram/{dataset}." + REFERENCE + "." + MAPPER_TAG + ".stats.log",
    threads: 16
    resources:
        mem_mb=32768,
    shell:
        r"""
        (
            set -eo pipefail

            echo "[$(date -Is)] START minimap2_ont_stats {wildcards.dataset}"
            echo "Container hostname: samtools-stats-{wildcards.dataset}"

            docker run --rm \
                --hostname "samtools-stats-{wildcards.dataset}" \
                --tmpfs /tmp:size=50g,exec \
                -u $UID:$(id -g) \
                --cpus {threads} \
                -m 32g \
                --workdir "{CWD}" \
                -v "{CWD}:{CWD}" \
                -v "{REFERENCE_DIR}:{REFERENCE_DIR}:ro" \
                --entrypoint samtools \
                "{DOCKER_SAMTOOLS}" \
                stats \
                -@ {threads} \
                --reference "{input.ref}" \
                --remove-overlaps \
                "{CWD}/{input.cram}" \
                > "{output.stats}"

            if [[ $(du -b "{output.stats}" | cut -f 1) -lt 5000 ]]; then
                echo "ERROR: samtools stats output is too small"
                exit 101
            fi

            echo "[$(date -Is)] END minimap2_ont_stats {wildcards.dataset}"
        ) > "{log}" 2>&1
        """
