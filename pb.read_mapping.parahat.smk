##
# pb.read_mapping.parahat.smk
# Read mapping workflow for PacBio HiFi data using ParaHAT
# mapper_tag: parahat-pb
# Constitution: Articles I–VIII
##

include: "header_mapper.smk"

####################
# Docker image

DOCKER_PARAHAT = "schimar/lrs-parahat:v1.0.0-cuda"

####################
# Reference

REF = (CWD + "/ref/"
    "GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta")

MAPPER_TAG = "parahat-pb"
REFERENCE  = "hg38"
PARAHAT_INDEX_DIR = CWD + "/parahat_index"

####################
# Discover inputs

FASTQ_DIR = "fastq"
DATASETS, = glob_wildcards(FASTQ_DIR + "/{dataset}.fastq.gz")
DATASETS = [d for d in DATASETS if ".pb." in d]
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

rule parahat_index:
    """Build the ParaHAT hash index for the reference (run once)."""
    input:
        ref = REF,
    output:
        sentinel = PARAHAT_INDEX_DIR + "/parahat_index.done",
    log:
        PARAHAT_INDEX_DIR + "/parahat_index.log",
    threads: 1
    shell:
        """
        mkdir -p {PARAHAT_INDEX_DIR}
        (
        echo "[$(date -Is)] START parahat_index" >&2

        docker run --rm \
            --workdir /tmp \
            -u $UID:$(id -g) \
            --cpus {threads} \
            -m 32g \
            -v {CWD}:{CWD} \
            -v {input.ref}:{input.ref}:ro \
            --entrypoint ParaHAT-indexer \
            {DOCKER_PARAHAT} \
            -k 13 \
            {PARAHAT_INDEX_DIR} \
            {input.ref}

        touch {output.sentinel}
        echo "[$(date -Is)] END parahat_index" >&2
        ) > {log} 2>&1
        """

rule parahat_map_sort:
    input:
        fastq   = FASTQ_DIR + "/{dataset}.fastq.gz",
        ref     = REF,
        idx_done = PARAHAT_INDEX_DIR + "/parahat_index.done",
    output:
        cram = "cram/{dataset}." + REFERENCE + "." + MAPPER_TAG + ".cram",
        crai = "cram/{dataset}." + REFERENCE + "." + MAPPER_TAG + ".cram.crai",
    log:
        "cram/{dataset}." + REFERENCE + "." + MAPPER_TAG + ".map_sort.log",
    threads: 16
    shell:
        """
        (
        set -eo pipefail
        echo "[$(date -Is)] START parahat_map_sort {wildcards.dataset}" >&2
        mkdir -p "{CWD}/cram/tmp"

        # ParaHAT-aligner → raw SAM file (no RG tag)
        # Decompress FASTQ first: ParaHAT does not handle .gz
        TMP_FASTQ="{CWD}/cram/tmp/{wildcards.dataset}.fastq"
        TMP_SAM="{CWD}/cram/tmp/{wildcards.dataset}.{REFERENCE}.{MAPPER_TAG}.raw.sam"

        zcat {CWD}/{input.fastq} > "$TMP_FASTQ"

        docker run --rm \
            --tmpfs /tmp:size=50g,exec \
            --shm-size 1g \
            --ulimit stack=-1 \
            -u $UID:$(id -g) \
            --cpus {threads} \
            -m 48g \
            -v {CWD}:{CWD} \
            -v {input.ref}:{input.ref}:ro \
            --entrypoint mpirun \
            {DOCKER_PARAHAT} \
            -n 1 ParaHAT-aligner \
            -t {threads} \
            {PARAHAT_INDEX_DIR} \
            "$TMP_FASTQ" \
            {input.ref} \
            > "$TMP_SAM"

        rm -f "$TMP_FASTQ"

        [[ -s "$TMP_SAM" ]] || {{ echo "ERROR: ParaHAT produced empty/missing $TMP_SAM"; exit 101; }}

        # samtools addreplacerg (inject @RG) → sort → CRAM
        docker run --rm \
            --workdir /tmp \
            -u $UID:$(id -g) \
            --cpus 4 \
            -m 8g \
            -v {CWD}:{CWD} \
            -v {input.ref}:{input.ref}:ro \
            --entrypoint samtools \
            {DOCKER_PARAHAT} \
            addreplacerg \
            -r "@RG\\tID:{wildcards.dataset}\\tSM:{wildcards.dataset}" \
            -o "{CWD}/cram/tmp/{wildcards.dataset}.{REFERENCE}.{MAPPER_TAG}.rg.sam" \
            "$TMP_SAM"

        rm -f "$TMP_SAM"

        docker run --rm \
            --workdir /tmp \
            -u $UID:$(id -g) \
            --cpus 4 \
            -m 16g \
            -v {CWD}:{CWD} \
            -v {input.ref}:{input.ref}:ro \
            --entrypoint samtools \
            {DOCKER_PARAHAT} \
            sort \
            -@ 4 \
            -O CRAM \
            --reference {input.ref} \
            -o {CWD}/{output.cram} \
            "{CWD}/cram/tmp/{wildcards.dataset}.{REFERENCE}.{MAPPER_TAG}.rg.sam"

        rm -f "{CWD}/cram/tmp/{wildcards.dataset}.{REFERENCE}.{MAPPER_TAG}.rg.sam"

        docker run --rm \
            --workdir /tmp \
            -u $UID:$(id -g) \
            --cpus 4 \
            -m 8g \
            -v {CWD}:{CWD} \
            -v {input.ref}:{input.ref}:ro \
            --entrypoint samtools \
            {DOCKER_PARAHAT} \
            index {CWD}/{output.cram}

        [[ $(du -b {output.cram} | cut -f 1) -le 64 ]] && exit 101

        echo "[$(date -Is)] END parahat_map_sort {wildcards.dataset}" >&2
        ) > {log} 2>&1
        """

rule parahat_idxstats:
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
        set -eo pipefail
        echo "[$(date -Is)] START parahat_idxstats {wildcards.dataset}" >&2

        docker run --rm \
            --workdir /tmp \
            -u $UID:$(id -g) \
            --cpus {threads} \
            -m 4g \
            -v {CWD}:{CWD} \
            -v {input.ref}:{input.ref}:ro \
            --entrypoint samtools \
            {DOCKER_PARAHAT} \
            idxstats {CWD}/{input.cram} \
            > {output.idxstats}

        [[ $(du -b {output.idxstats} | cut -f 1) -lt 5000 ]] && exit 101

        echo "[$(date -Is)] END parahat_idxstats {wildcards.dataset}" >&2
        ) > {log} 2>&1
        """

rule parahat_stats:
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
        set -eo pipefail
        echo "[$(date -Is)] START parahat_stats {wildcards.dataset}" >&2

        docker run --rm \
            --workdir /tmp \
            -u $UID:$(id -g) \
            --cpus {threads} \
            -m 16g \
            -v {CWD}:{CWD} \
            -v {input.ref}:{input.ref}:ro \
            --entrypoint samtools \
            {DOCKER_PARAHAT} \
            stats \
            --reference {input.ref} \
            {CWD}/{input.cram} \
            > {output.stats}

        [[ $(du -b {output.stats} | cut -f 1) -lt 5000 ]] && exit 101

        echo "[$(date -Is)] END parahat_stats {wildcards.dataset}" >&2
        ) > {log} 2>&1
        """