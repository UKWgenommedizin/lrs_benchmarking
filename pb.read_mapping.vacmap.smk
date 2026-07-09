##
# pb.read_mapping.vacmap.smk
# Read mapping workflow for PacBio HiFi data using VACmap
# mapper_tag: vacmap-pb
# Constitution: Articles I–VIII
##

include: "header.smk"

####################
# Docker image

DOCKER_VACMAP = "schimar/lrs-vacmap:v1.2.0"

####################
# Reference

REF = os.path.expanduser(
    "~/smb/Analyses/Reference_sequence/hg38_KGGM/"
    "GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta"
)

MAPPER_TAG = "vacmap-pb"
REFERENCE  = "hg38"

####################
# Discover inputs

FASTQ_DIR = "fastq"
DATASETS, = glob_wildcards(FASTQ_DIR + "/{dataset}.fastq.gz")

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

        # Map with VACmap, pipe to samtools sort → CRAM
        docker run --rm \
            --workdir /tmp \
            -u $UID:$(id -g) \
            --cpus {threads} \
            -m 48g \
            -v {CWD}:{CWD} \
            -v {input.ref}:{input.ref}:ro \
            --entrypoint vacmap \
            {DOCKER_VACMAP} \
            -t {threads} \
            -x hifi \
            --rg "@RG\\tID:{wildcards.dataset}\\tSM:{wildcards.dataset}" \
            {input.ref} \
            {CWD}/{input.fastq} \
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

        # Validate CRAM size
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
