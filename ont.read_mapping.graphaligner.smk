##
# ont.read_mapping.graphaligner.smk
# Read mapping workflow for ONT data using GraphAligner
# mapper_tag: ga-ont
# Constitution: Articles I–VIII
#
# GraphAligner CLI (github.com/maickrau/GraphAligner):
#   GraphAligner -g <graph.gfa> -f <reads.fastq.gz> -a <out.gaf> -x vg -t N
#
#   There is NO --read-group / --rg flag in GraphAligner.
#   Output is GAF (Graph Alignment Format), not SAM/BAM/CRAM.
#
# Pipeline (file-based, no docker pipes):
#   Step 1  vg convert  : hg38.fasta → hg38.gfa  (once per ref)
#   Step 2  GraphAligner: reads + hg38.gfa → GAF file
#   Step 3  vg surject  : GAF + hg38.gfa → SAM file
#   Step 4  samtools    : addreplacerg → sort → CRAM → index
#
# Required Docker images:
#   DOCKER_GRAPHALIGNER  must contain: GraphAligner binary
#   DOCKER_VG            must contain: vg binary + samtools + bcftools
##

include: "header_mapper.smk"

####################
# Docker images

DOCKER_GRAPHALIGNER = "schimar/lrs-graphaligner:v1.0.20"
DOCKER_VG           = "schimar/lrs-vg:v1.73.0"

####################
# Reference

REF = (CWD + "/ref/"
    "GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta")

# GFA derived from REF — built once by rule ref_to_gfa below
GRAPH_GFA = CWD + "/graphaligner_graph/hg38.gfa"

MAPPER_TAG = "ga-ont"
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

rule ref_to_gfa:
    """Convert the linear hg38 FASTA to a single-path GFA (run once)."""
    input:
        ref = REF,
    output:
        gfa = GRAPH_GFA,
    log:
        CWD + "/graphaligner_graph/ref_to_gfa.log",
    threads: 4
    shell:
        """
        (
        set -eo pipefail
        echo "[$(date -Is)] START ref_to_gfa" >&2
        mkdir -p {CWD}/graphaligner_graph

        docker run --rm \
            --workdir /tmp \
            -u $UID:$(id -g) \
            --cpus {threads} \
            -m 16g \
            -v {CWD}:{CWD} \
            -v {input.ref}:{input.ref}:ro \
            --entrypoint vg \
            {DOCKER_VG} \
            construct \
            -r {input.ref} \
            --gfa-out \
            {input.ref} \
            > {output.gfa}

        [[ ! -s {output.gfa} ]] && exit 101

        echo "[$(date -Is)] END ref_to_gfa" >&2
        ) > {log} 2>&1
        """

rule graphaligner_map:
    """Align ONT reads to the hg38 GFA → GAF file."""
    input:
        fastq = FASTQ_DIR + "/{dataset}.fastq.gz",
        gfa   = GRAPH_GFA,
    output:
        gaf   = temp("cram/tmp/{dataset}." + REFERENCE + "." + MAPPER_TAG + ".gaf"),
    log:
        "cram/tmp/{dataset}." + REFERENCE + "." + MAPPER_TAG + ".graphaligner.log",
    threads: 16
    shell:
        """
        (
        set -eo pipefail
        echo "[$(date -Is)] START graphaligner_map {wildcards.dataset}" >&2
        mkdir -p cram/tmp

        docker run --rm \
            --workdir /tmp \
            -u $UID:$(id -g) \
            --cpus {threads} \
            -m 64g \
            -v {CWD}:{CWD} \
            --entrypoint GraphAligner \
            {DOCKER_GRAPHALIGNER} \
            -g {CWD}/{input.gfa} \
            -f {CWD}/{input.fastq} \
            -a {CWD}/{output.gaf} \
            -x vg \
            -t {threads} \
            --precise-clipping 0.75

        [[ ! -s {output.gaf} ]] && exit 101

        echo "[$(date -Is)] END graphaligner_map {wildcards.dataset}" >&2
        ) > {log} 2>&1
        """

rule graphaligner_surject_sort:
    """Surject GAF → SAM → addreplacerg → sort → CRAM → index (file-based)."""
    input:
        gaf = "cram/tmp/{dataset}." + REFERENCE + "." + MAPPER_TAG + ".gaf",
        gfa = GRAPH_GFA,
        ref = REF,
    output:
        cram = "cram/{dataset}." + REFERENCE + "." + MAPPER_TAG + ".cram",
        crai = "cram/{dataset}." + REFERENCE + "." + MAPPER_TAG + ".cram.crai",
    log:
        "cram/{dataset}." + REFERENCE + "." + MAPPER_TAG + ".surject_sort.log",
    threads: 8
    shell:
        """
        (
        set -eo pipefail
        echo "[$(date -Is)] START graphaligner_surject_sort {wildcards.dataset}" >&2
        mkdir -p "{CWD}/cram/tmp"

        TMP_SAM="{CWD}/cram/tmp/{wildcards.dataset}.{REFERENCE}.{MAPPER_TAG}.raw.sam"
        TMP_RG="{CWD}/cram/tmp/{wildcards.dataset}.{REFERENCE}.{MAPPER_TAG}.rg.sam"

        # vg surject: GAF → SAM
        docker run --rm \
            --workdir /tmp \
            -u $UID:$(id -g) \
            --cpus {threads} \
            -m 32g \
            -v {CWD}:{CWD} \
            -v {input.ref}:{input.ref}:ro \
            --entrypoint vg \
            {DOCKER_VG} \
            surject \
            -x {CWD}/{input.gfa} \
            -G \
            -s \
            {CWD}/{input.gaf} \
            > "$TMP_SAM"

        [[ -s "$TMP_SAM" ]] || {{ echo "ERROR: vg surject produced empty/missing $TMP_SAM"; exit 101; }}

        # samtools addreplacerg (inject @RG)
        docker run --rm \
            --workdir /tmp \
            -u $UID:$(id -g) \
            --cpus 2 \
            -m 4g \
            -v {CWD}:{CWD} \
            -v {input.ref}:{input.ref}:ro \
            --entrypoint samtools \
            {DOCKER_VG} \
            addreplacerg \
            -r "@RG\\tID:{wildcards.dataset}\\tSM:{wildcards.dataset}" \
            -o "$TMP_RG" \
            "$TMP_SAM"

        rm -f "$TMP_SAM"

        # samtools sort → CRAM
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
            "$TMP_RG"

        rm -f "$TMP_RG"

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

        echo "[$(date -Is)] END graphaligner_surject_sort {wildcards.dataset}" >&2
        ) > {log} 2>&1
        """

rule graphaligner_idxstats:
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
        echo "[$(date -Is)] START graphaligner_idxstats {wildcards.dataset}" >&2

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

        echo "[$(date -Is)] END graphaligner_idxstats {wildcards.dataset}" >&2
        ) > {log} 2>&1
        """

rule graphaligner_stats:
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
        echo "[$(date -Is)] START graphaligner_stats {wildcards.dataset}" >&2

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

        echo "[$(date -Is)] END graphaligner_stats {wildcards.dataset}" >&2
        ) > {log} 2>&1
        """