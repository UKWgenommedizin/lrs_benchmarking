##
# ont.read_mapping.vacmap.smk
# Read mapping workflow for ONT data using VACmap
# mapper_tag: vacmap-ont
# Constitution: Articles I–VIII
#
# VACmap CLI (from github.com/micahvista/VACmap, v1.0.0):
#   vacmap -ref <ref.fasta> -read <reads.fastq.gz> -mode H -t <threads> \
#          --rg-id <ID> --rg-sm <SM> | samtools sort ...
#
#   -mode H  = high error rate long reads (ONT / PacBio CLR)  ← correct for ONT
#   -mode L  = low error rate (HiFi)
#   -mode S  = increased sensitivity for small variants
#
#   Read-group flags (all separate, unlike minimap2's single --rg string):
#     --rg-id <string>   RG ID tag
#     --rg-sm <string>   RG SM (sample) tag
#
#   Output: SAM to stdout by default; pipe to samtools sort → CRAM.
#
# NOTE: VACmap is a Python tool. The Docker image must have Python + VACmap
#       installed (conda env or pip install from source), plus samtools/bcftools.
##

include: "header_mapper.smk"

####################
# Docker image
# Must contain: vacmap (Python), samtools, bcftools (Art. VII.1)

DOCKER_VACMAP = "schimar/lrs-vacmap:v1.2.0"

####################
# Reference

REF = (CWD + "/ref/"
    "GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta")

MAPPER_TAG = "vacmap-ont"
REFERENCE  = "hg38"

####################
# Discover inputs

FASTQ_DIR = "fastq"
DATASETS, = glob_wildcards(FASTQ_DIR + "/{dataset}.fastq.gz")
if DATASET_FILTER:
    DATASETS = [d for d in DATASETS if DATASET_FILTER in d]

####################
# Targets

rule all:
    input:
        expand(
            "cram/{dataset}.{ref}.{tag}.cram",
            dataset=DATASETS, ref=REFERENCE, tag=MAPPER_TAG
        ),
        expand(
            "cram/{dataset}.{ref}.{tag}.cram.crai",
            dataset=DATASETS, ref=REFERENCE, tag=MAPPER_TAG
        ),
        expand(
            "cram/{dataset}.{ref}.{tag}.cram.idxstats",
            dataset=DATASETS, ref=REFERENCE, tag=MAPPER_TAG
        ),
        expand(
            "cram/{dataset}.{ref}.{tag}.cram.stats",
            dataset=DATASETS, ref=REFERENCE, tag=MAPPER_TAG
        ),

####################
# Rules

rule vacmap_map_sort:
    input:
        fastq = FASTQ_DIR + "/{dataset}.fastq.gz",
        ref   = REF,
    output:
        cram  = "cram/{dataset}." + REFERENCE + "." + MAPPER_TAG + ".cram",
        crai  = "cram/{dataset}." + REFERENCE + "." + MAPPER_TAG + ".cram.crai",
    log:
        "cram/{dataset}." + REFERENCE + "." + MAPPER_TAG + ".map_sort.log",
    threads: 16
    shell:
        """
        (
        echo "[$(date -Is)] START vacmap_map_sort {wildcards.dataset}" >&2

        # VACmap → SAM stdout → samtools sort → CRAM
        # -mode H: high error rate ONT reads
        # --rg-id / --rg-sm: per-field RG tags (VACmap has no single --rg string)
        docker run --rm \
            --workdir /tmp \
            -u $UID:$(id -g) \
            --cpus {threads} \
            -m 24g \
            -v {CWD}:{CWD} \
            -v {input.ref}:{input.ref}:ro \
            --entrypoint vacmap \
            {DOCKER_VACMAP} \
            -ref {input.ref} \
            -read {CWD}/{input.fastq} \
            -mode H \
            -t {threads} \
            --rg-id {wildcards.dataset} \
            --rg-sm {wildcards.dataset} \
        | docker run --rm \
            --workdir /tmp \
            -u $UID:$(id -g) \
            --cpus 4 \
            -m 16g \
            -v {CWD}:{CWD} \
            -v {input.ref}:{input.ref}:ro \
            --entrypoint samtools \
            {DOCKER_VACMAP} \
            sort \
            -@ 4 \
            -O CRAM \
            --reference {input.ref} \
            -o {CWD}/{output.cram} \
            -

        docker run --rm \
            --workdir /tmp \
            -u $UID:$(id -g) \
            --cpus 4 \
            -m 8g \
            -v {CWD}:{CWD} \
            -v {input.ref}:{input.ref}:ro \
            --entrypoint samtools \
            {DOCKER_VACMAP} \
            index {CWD}/{output.cram}

        [[ $(du -b {output.cram} | cut -f 1) -le 64 ]] && exit 101

        echo "[$(date -Is)] END vacmap_map_sort {wildcards.dataset}" >&2
        ) > {log} 2>&1
        """

rule vacmap_idxstats:
    input:
        cram = "cram/{dataset}." + REFERENCE + "." + MAPPER_TAG + ".cram",
        crai = "cram/{dataset}." + REFERENCE + "." + MAPPER_TAG + ".cram.crai",
        ref  = REF,
    output:
        idxstats = "cram/{dataset}." + REFERENCE + "." + MAPPER_TAG + ".cram.idxstats",
    log:
        "cram/{dataset}." + REFERENCE + "." + MAPPER_TAG + ".idxstats.log",
    threads: 2
    shell:
        """
        (
        echo "[$(date -Is)] START vacmap_idxstats {wildcards.dataset}" >&2

        docker run --rm \
            --workdir /tmp \
            -u $UID:$(id -g) \
            --cpus {threads} \
            -m 4g \
            -v {CWD}:{CWD} \
            -v {input.ref}:{input.ref}:ro \
            --entrypoint samtools \
            {DOCKER_VACMAP} \
            idxstats {CWD}/{input.cram} \
            > {output.idxstats}

        [[ $(du -b {output.idxstats} | cut -f 1) -lt 5000 ]] && exit 101

        echo "[$(date -Is)] END vacmap_idxstats {wildcards.dataset}" >&2
        ) > {log} 2>&1
        """

rule vacmap_stats:
    input:
        cram = "cram/{dataset}." + REFERENCE + "." + MAPPER_TAG + ".cram",
        crai = "cram/{dataset}." + REFERENCE + "." + MAPPER_TAG + ".cram.crai",
        ref  = REF,
    output:
        stats = "cram/{dataset}." + REFERENCE + "." + MAPPER_TAG + ".cram.stats",
    log:
        "cram/{dataset}." + REFERENCE + "." + MAPPER_TAG + ".stats.log",
    threads: 4
    shell:
        """
        (
        echo "[$(date -Is)] START vacmap_stats {wildcards.dataset}" >&2

        docker run --rm \
            --workdir /tmp \
            -u $UID:$(id -g) \
            --cpus {threads} \
            -m 8g \
            -v {CWD}:{CWD} \
            -v {input.ref}:{input.ref}:ro \
            --entrypoint samtools \
            {DOCKER_VACMAP} \
            stats \
            --reference {input.ref} \
            {CWD}/{input.cram} \
            > {output.stats}

        [[ $(du -b {output.stats} | cut -f 1) -lt 5000 ]] && exit 101

        echo "[$(date -Is)] END vacmap_stats {wildcards.dataset}" >&2
        ) > {log} 2>&1
        """
