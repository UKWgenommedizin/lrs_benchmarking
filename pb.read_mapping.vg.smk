##
# pb.read_mapping.vg.smk
# Read mapping workflow for PacBio HiFi data using VG (vg giraffe / vg map)
# mapper_tag: vg-pb
# Constitution: Articles I–VIII
#
# VG requires a graph index (GBZ/XG/GCSA2/dist) in addition to the
# linear reference. This workflow uses `vg giraffe` for speed on
# ONT reads. If you prefer `vg map`, replace the giraffe block with
# the vg map equivalent and provide the appropriate index files.
#
# Required index files (place under vg_index/ in CWD):
#   vg_index/hg38.giraffe.gbz
#   vg_index/hg38.dist
#   vg_index/hg38.min
##

include: "header_mapper.smk"

####################
# Docker image

DOCKER_VG = "schimar/lrs-vg:v1.73.0"

####################
# Reference / index paths

REF = (CWD + "/ref/"
    "GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta")

VG_INDEX_DIR = CWD + "/vg_index"
VG_GBZ       = VG_INDEX_DIR + "/hg38.giraffe.gbz"
VG_DIST      = VG_INDEX_DIR + "/hg38.dist"
VG_MIN       = VG_INDEX_DIR + "/hg38.longread.withzip.min"
VG_ZIPCODES  = VG_INDEX_DIR + "/hg38.longread.zipcodes"

MAPPER_TAG = "vg-pb"
REFERENCE  = "hg38"

####################
# Samtools utility image (bundled in vg image)

DOCKER_SAMTOOLS = DOCKER_VG

####################
# Discover inputs

FASTQ_DIR = "fastq"
DATASETS, = glob_wildcards(FASTQ_DIR + "/{dataset}.fastq.gz")
if DATASET_FILTER:
    DATASETS = [d for d in DATASETS if DATASET_FILTER in d]

####################
# Rule: build VG index if not present

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
    shell:
        """
        (
        echo "[$(date -Is)] START build_vg_index" >&2
        mkdir -p {VG_INDEX_DIR}

        docker run --rm \
            --workdir /tmp \
            -u $UID:$(id -g) \
            --cpus {threads} \
            -m 64g \
            -v {CWD}:{CWD} \
            -v {input.ref}:{input.ref}:ro \
            --entrypoint vg \
            {DOCKER_VG} \
            autoindex \
            --workflow lr-giraffe \
            --ref {input.ref} \
            --prefix {VG_INDEX_DIR}/hg38 \
            --threads {threads}

        for f in {VG_GBZ} {VG_DIST} {VG_MIN} {VG_ZIPCODES}; do
            [[ -s "$f" ]] || {{ echo "Missing/empty: $f"; exit 101; }}
        done

        echo "[$(date -Is)] END build_vg_index" >&2
        ) > {log} 2>&1
        """

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
        fastq = FASTQ_DIR + "/{dataset}.fastq.gz",
        ref   = REF,
        gbz   = VG_GBZ,
        dist  = VG_DIST,
        min_idx = VG_MIN,
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
        echo "[$(date -Is)] START vg_map_sort {wildcards.dataset}" >&2

        # vg giraffe → surject to linear reference → SAM → samtools sort → CRAM
        docker run --rm \
            --workdir /tmp \
            -u $UID:$(id -g) \
            --cpus {threads} \
            -m 64g \
            -v {CWD}:{CWD} \
            -v {input.ref}:{input.ref}:ro \
            --entrypoint vg \
            {DOCKER_VG} \
            giraffe \
            -t {threads} \
            -Z {input.gbz} \
            -d {input.dist} \
            -m {input.min_idx} \
            -z {input.zipcodes} \
            -b hifi \
            -f {CWD}/{input.fastq} \
            --read-group "ID:{wildcards.dataset}\tSM:{wildcards.dataset}" \
            --sample {wildcards.dataset} \
            --output-format SAM \
        | docker run --rm \
            --workdir /tmp \
            -u $UID:$(id -g) \
            --cpus 4 \
            -m 16g \
            -v {CWD}:{CWD} \
            -v {input.ref}:{input.ref}:ro \
            --entrypoint samtools \
            {DOCKER_SAMTOOLS} \
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
            {DOCKER_SAMTOOLS} \
            index {CWD}/{output.cram}

        # Validate CRAM size
        [[ $(du -b {output.cram} | cut -f 1) -le 64 ]] && exit 101

        echo "[$(date -Is)] END vg_map_sort {wildcards.dataset}" >&2
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
        echo "[$(date -Is)] START vg_idxstats {wildcards.dataset}" >&2

        docker run --rm \
            --workdir /tmp \
            -u $UID:$(id -g) \
            --cpus {threads} \
            -m 4g \
            -v {CWD}:{CWD} \
            -v {input.ref}:{input.ref}:ro \
            --entrypoint samtools \
            {DOCKER_SAMTOOLS} \
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
        echo "[$(date -Is)] START vg_stats {wildcards.dataset}" >&2

        docker run --rm \
            --workdir /tmp \
            -u $UID:$(id -g) \
            --cpus {threads} \
            -m 8g \
            -v {CWD}:{CWD} \
            -v {input.ref}:{input.ref}:ro \
            --entrypoint samtools \
            {DOCKER_SAMTOOLS} \
            stats \
            --reference {input.ref} \
            {CWD}/{input.cram} \
            > {output.stats}

        [[ $(du -b {output.stats} | cut -f 1) -lt 5000 ]] && exit 101

        echo "[$(date -Is)] END vg_stats {wildcards.dataset}" >&2
        ) > {log} 2>&1
        """
