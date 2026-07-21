##
# ont.read_mapping.vg.smk
# Read mapping workflow for ONT data using VG Giraffe (long-read mode)
# mapper_tag: vg-ont
# Constitution: Articles I–VIII
#
# VG Giraffe long-read CLI (vgteam/vg v1.63.0+):
#   vg giraffe -Z <gbz> -b r10 -f <fastq> -o SAM -t N -R <rg> -N <sample>
#
# Required indexes (pre-build with: vg autoindex --workflow lr-giraffe ...):
#   vg_index/hg38.giraffe.gbz
#   vg_index/hg38.dist
#   vg_index/hg38.longread.withzip.min
#   vg_index/hg38.longread.zipcodes
##

include: "header_mapper.smk"

####################
# Docker image (must be v1.63.0+ for long-read mode)

DOCKER_VG = "schimar/lrs-vg:v1.73.0"

####################
# Reference & index paths

REF = (CWD + "/ref/"
    "GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta")

VG_INDEX_DIR = CWD + "/vg_index"
VG_GBZ       = VG_INDEX_DIR + "/hg38.gbz"
VG_DIST      = VG_INDEX_DIR + "/hg38.dist"
VG_MIN       = VG_INDEX_DIR + "/hg38.longread.withzip.min"
VG_ZIPCODES  = VG_INDEX_DIR + "/hg38.longread.zipcodes"

MAPPER_TAG = "vg-ont"
REFERENCE  = "hg38"

####################
# Discover inputs

FASTQ_DIR = "fastq"
DATASETS, = glob_wildcards(FASTQ_DIR + "/{dataset}.fastq.gz")
DATASETS = [d for d in DATASETS if ".ont." in d]
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

rule vg_map_sort:
    input:
        fastq    = FASTQ_DIR + "/{dataset}.fastq.gz",
        ref      = REF,
        gbz      = VG_GBZ,
        dist     = VG_DIST,
        min_idx  = VG_MIN,
        zipcodes = VG_ZIPCODES,
    output:
        cram  = "cram/{dataset}." + REFERENCE + "." + MAPPER_TAG + ".cram",
        crai  = "cram/{dataset}." + REFERENCE + "." + MAPPER_TAG + ".cram.crai",
    log:
        "cram/{dataset}." + REFERENCE + "." + MAPPER_TAG + ".map_sort.log",
    threads: 16
    shell:
        """
        (
        set -eo pipefail
        echo "[$(date -Is)] START vg_map_sort {wildcards.dataset}" >&2
        mkdir -p "{CWD}/cram/tmp"

        TMP_SAM="{CWD}/cram/tmp/{wildcards.dataset}.{REFERENCE}.{MAPPER_TAG}.sam"

        # vg giraffe → SAM file
        docker run --rm \
            --tmpfs /tmp:size=50g,exec \
            -u $UID:$(id -g) \
            --cpus {threads} \
            -m 64g \
            -v {CWD}:{CWD} \
            -v {input.ref}:{input.ref}:ro \
            --entrypoint vg \
            {DOCKER_VG} \
            giraffe \
            -t {threads} \
            -Z {CWD}/{input.gbz} \
            -d {CWD}/{input.dist} \
            -m {CWD}/{input.min_idx} \
            -z {CWD}/{input.zipcodes} \
            -b r10 \
            -f {CWD}/{input.fastq} \
            -o SAM \
            -R "ID:{wildcards.dataset}\tSM:{wildcards.dataset}" \
            -N {wildcards.dataset} \
            > "$TMP_SAM"

        [[ -s "$TMP_SAM" ]] || {{ echo "ERROR: vg giraffe produced empty/missing $TMP_SAM"; exit 101; }}

        # samtools sort SAM → CRAM
        docker run --rm \
            --workdir /tmp \
            -u $UID:$(id -g) \
            --cpus 4 \
            -m 16g \
            -v {CWD}:{CWD} \
            -v {input.ref}:{input.ref}:ro \
            --entrypoint samtools \
            {DOCKER_VG} \
            sort \
            -@ 4 \
            -O CRAM \
            --reference {input.ref} \
            -o {CWD}/{output.cram} \
            "$TMP_SAM"

        rm -f "$TMP_SAM"

        docker run --rm \
            --workdir /tmp \
            -u $UID:$(id -g) \
            --cpus 4 \
            -m 8g \
            -v {CWD}:{CWD} \
            -v {input.ref}:{input.ref}:ro \
            --entrypoint samtools \
            {DOCKER_VG} \
            index {CWD}/{output.cram}

        [[ $(du -b {output.cram} | cut -f 1) -le 64 ]] && exit 101

        echo "[$(date -Is)] END vg_map_sort {wildcards.dataset}" >&2
        ) > {log} 2>&1
        """

####################
# Rule: build VG long-read indexes (run once)

rule build_vg_index:
    """Build VG giraffe long-read indexes from the reference FASTA."""
    input:
        ref = REF,
    output:
        gbz      = VG_GBZ,
        dist     = VG_DIST,
        min_idx  = VG_MIN,
        zipcodes = VG_ZIPCODES,
    log:
        VG_INDEX_DIR + "/build_index.log",
    threads: 16
    resources:
        mem_mb = 225280,
    shell:
        """
        mkdir -p {VG_INDEX_DIR}
        (
        set -eo pipefail
        echo "[$(date -Is)] START build_vg_index" >&2

        docker run --rm \
            --tmpfs /tmp:size=50g,exec \
            -u $UID:$(id -g) \
            --cpus {threads} \
            -m 220g \
            -v {CWD}:{CWD} \
            -v {input.ref}:{input.ref}:ro \
            --entrypoint vg \
            {DOCKER_VG} \
            autoindex \
            --workflow lr-giraffe \
            --ref {input.ref} \
            --prefix {VG_INDEX_DIR}/hg38 \
            --target-mem 200G \
            --threads {threads}
        # Large indexes (.gbz ~2.3GB) may not be visible immediately after docker exits
        # due to fsync/buffering on bind mounts. Poll with backoff.
        for attempt in $(seq 1 30); do
            missing=0
            for f in {VG_GBZ} {VG_DIST} {VG_MIN} {VG_ZIPCODES}; do
                [[ -s "$f" ]] || missing=1
            done
            [[ $missing -eq 0 ]] && break
            echo "Waiting for index files to flush (attempt $attempt)..." >&2
            sleep 2
        done
        for f in {VG_GBZ} {VG_DIST} {VG_MIN} {VG_ZIPCODES}; do
            [[ -s "$f" ]] || {{ echo "Missing/empty after 30 retries: $f"; exit 101; }}
        done

        echo "[$(date -Is)] END build_vg_index" >&2
        ) > {log} 2>&1
        """

rule vg_idxstats:
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
        echo "[$(date -Is)] START vg_idxstats {wildcards.dataset}" >&2

        docker run --rm \
            --workdir /tmp \
            -u $UID:$(id -g) \
            --cpus {threads} \
            -m 4g \
            -v {CWD}:{CWD} \
            -v {input.ref}:{input.ref}:ro \
            --entrypoint samtools \
            {DOCKER_VG} \
            idxstats {CWD}/{input.cram} \
            > {output.idxstats}

        [[ $(du -b {output.idxstats} | cut -f 1) -lt 5000 ]] && exit 101

        echo "[$(date -Is)] END vg_idxstats {wildcards.dataset}" >&2
        ) > {log} 2>&1
        """

rule vg_stats:
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
        echo "[$(date -Is)] START vg_stats {wildcards.dataset}" >&2

        docker run --rm \
            --workdir /tmp \
            -u $UID:$(id -g) \
            --cpus {threads} \
            -m 16g \
            -v {CWD}:{CWD} \
            -v {input.ref}:{input.ref}:ro \
            --entrypoint samtools \
            {DOCKER_VG} \
            stats \
            --reference {input.ref} \
            {CWD}/{input.cram} \
            > {output.stats}

        [[ $(du -b {output.stats} | cut -f 1) -lt 5000 ]] && exit 101

        echo "[$(date -Is)] END vg_stats {wildcards.dataset}" >&2
        ) > {log} 2>&1
        """