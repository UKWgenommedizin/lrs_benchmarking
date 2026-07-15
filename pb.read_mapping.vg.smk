##
# pb.read_mapping.vg.smk
# Read mapping workflow for PacBio HiFi data using VG Giraffe
# mapper_tag: vg-pb
# Constitution: Articles I–VIII
#
# VG Giraffe long-read CLI (vgteam/vg v1.63.0+):
#   vg giraffe -Z <gbz> -b hifi -f <fastq> -o SAM -t N -R <rg> -N <sample>
#
# Required indexes (pre-build with: vg autoindex --workflow lr-giraffe ...):
#   vg_index/hg38.giraffe.gbz
#   vg_index/hg38.dist
#   vg_index/hg38.longread.withzip.min
#   vg_index/hg38.longread.zipcodes
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
VG_GBZ       = VG_INDEX_DIR + "/hg38.gbz"
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
    resources:
        mem_mb = 131072,
    shell:
        """
        mkdir -p {VG_INDEX_DIR}
        exec 3>&1 4>&2 1> >(tee -a {log}) 2>&1
        echo "[$(date -Is)] START build_vg_index"

        docker run --rm \
            --workdir /tmp \
            -u $UID:$(id -g) \
            --cpus {threads} \
            -m 120g \
            --memory-swap 120g \
            -v {CWD}:{CWD} \
            -v {input.ref}:{input.ref}:ro \
            --entrypoint vg \
            {DOCKER_VG} \
            autoindex \
            --workflow lr-giraffe \
            --ref {input.ref} \
            --prefix {VG_INDEX_DIR}/hg38 \
            --target-mem 100G \
            --threads {threads}
        DOCKER_EXIT=$?
        echo "[$(date -Is)] docker exit code: $DOCKER_EXIT"

        if [[ $DOCKER_EXIT -ne 0 ]]; then
            echo "ERROR: vg autoindex failed with exit code $DOCKER_EXIT"
            exit $DOCKER_EXIT
        fi

        for f in {VG_GBZ} {VG_DIST} {VG_MIN} {VG_ZIPCODES}; do
            if [[ -s "$f" ]]; then
                echo "OK: $f ($(du -h "$f" | cut -f 1))"
            else
                echo "MISSING/EMPTY: $f"
                exit 101
            fi
        done

        echo "[$(date -Is)] END build_vg_index"
        """