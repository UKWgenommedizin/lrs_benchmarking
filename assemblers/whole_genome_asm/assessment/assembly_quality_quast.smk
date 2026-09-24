#
# assembly_quality_quast.smk
#
# Whole-genome assembly quality assessment with QUAST-LG.
#
# Evaluates existing Flye, GoldRush and Verkko assemblies.

import os


# ************************************************************************************************
# Working repository
# ************************************************************************************************

CWD = os.path.abspath(os.getcwd())

YU_ROOT = "/data/genmedbfx/yu_j/lrs_benchmarking"
STOIBER_ROOT = "/home/stoiber_l/smbshare/lrs_benchmarking"

print("Current working directory: " + CWD)
print("Yu repository: " + YU_ROOT)
print("Stoiber assembly repository: " + STOIBER_ROOT)


# ************************************************************************************************
# QUAST configuration
# ************************************************************************************************

QUAST_VERSION = "5.3.0"

DOCKER_QUAST = (
    "quay.io/biocontainers/"
    "quast:5.3.0--py313pl5321h5ca1c30_2"
)

print("QUAST version: " + QUAST_VERSION)
print("QUAST Docker image: " + DOCKER_QUAST)


# ************************************************************************************************
# Genome reference
# ************************************************************************************************

DEFAULT_REFERENCE = (
    "/data/genmedbfx/schilling_m/repos/lrs_benchmarking/ref/"
    "GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_"
    "MAP2K3_KMT2C_KCNJ18.fasta"
)

RAW_REFERENCE = config.get("reference", DEFAULT_REFERENCE)

REFERENCE = os.path.abspath(
    os.path.expanduser(RAW_REFERENCE)
)

if not os.path.isfile(REFERENCE):
    raise ValueError(
        "Reference genome does not exist: " + REFERENCE
    )

REFERENCE_DIR = os.path.dirname(REFERENCE)

print("Reference genome: " + REFERENCE)


# ************************************************************************************************
# Assembly locations
# ************************************************************************************************

ASSEMBLY_ROOTS = {
    "flye": os.path.join(
        YU_ROOT,
        "assemblies",
        "flye"
    ),
    "goldrush": os.path.join(
        STOIBER_ROOT,
     #   "assemblies",
        "goldrush"
    ),
    "verkko": os.path.join(
        STOIBER_ROOT,
     #   "assemblies",
        "verkko"
    ),
}

for assembler, root in ASSEMBLY_ROOTS.items():
    print(f"{assembler} assembly root: {root}")


# ************************************************************************************************
# Production filtering
# ************************************************************************************************

TEST_MARKERS = (
    ".1k",
    ".chr21",
    ".localtest",
    "smoke",
)

PRODUCTION_SAMPLES = {
    "hg002",
    "hg003",
    "hg004",
}


def is_production_assembly(assembler, dataset):

    assembler = assembler.lower()
    dataset = dataset.lower()

    if any(marker in dataset for marker in TEST_MARKERS):
        return False

    # Verkko is hybrid and should have one assembly per GIAB sample.
    if assembler == "verkko":
        return dataset in PRODUCTION_SAMPLES

    # Flye and GoldRush are technology-specific:
    # HG002.ont.30x, HG002.pb.30x, etc.
    if assembler in {"flye", "goldrush"}:
        return ".30x" in dataset

    return False


# ************************************************************************************************
# Discover completed assembly FASTAs
# ************************************************************************************************

ASSEMBLIES = []

for assembler, root in ASSEMBLY_ROOTS.items():

    if not os.path.isdir(root):
        print(
            "WARNING: assembly directory does not exist: "
            + root
        )
        continue

    datasets, = glob_wildcards(
        os.path.join(
            root,
            "{dataset}",
            "assembly.fasta"
        )
    )

    for dataset in datasets:

        if is_production_assembly(
            assembler,
            dataset
        ):
            ASSEMBLIES.append(
                (assembler, dataset)
            )


ASSEMBLIES = sorted(set(ASSEMBLIES))


if not ASSEMBLIES:
    raise ValueError(
        "No completed production Flye, GoldRush or Verkko "
        "assembly FASTAs were discovered."
    )


print("Discovered production assemblies:")

for assembler, dataset in ASSEMBLIES:
    print(
        "  "
        + assembler
        + "\t"
        + dataset
    )


# ************************************************************************************************
# Resolve assembly FASTA
# ************************************************************************************************

def assembly_path(wildcards):

    root = ASSEMBLY_ROOTS[
        wildcards.assembler
    ]

    path = os.path.join(
        root,
        wildcards.dataset,
        "assembly.fasta"
    )

    if not os.path.isfile(path):
        raise ValueError(
            "Assembly FASTA does not exist: "
            + path
        )

    return path


# ************************************************************************************************
# QUAST targets
# ************************************************************************************************

OUTPUT = [
    f"assembly_quality/quast/{assembler}/{dataset}/report.tsv"
    for assembler, dataset in ASSEMBLIES
]


print("QUAST targets:")

for target in OUTPUT:
    print("  " + target)


rule all:
    input:
        OUTPUT


# ************************************************************************************************
# QUAST-LG
# ************************************************************************************************

rule quast_assembly:
    input:
        assembly=assembly_path,
        reference=REFERENCE

    output:
        quast_tsv="assembly_quality/quast/{assembler}/{dataset}/report.tsv"

    params:
        outdir="assembly_quality/quast/{assembler}/{dataset}"

    log:
        "assembly_quality/quast/{assembler}/{dataset}/quast.log"

    threads: 16

    resources:
        mem_gb=128

    message:
        "Evaluating {wildcards.assembler} {wildcards.dataset} with QUAST-LG"

    shell:
        r"""
        set -eo pipefail

        mkdir -p "{params.outdir}"
        : > "{log}"

        docker run --rm \
            --cpus {threads} \
            -m {resources.mem_gb}g \
            --tmpfs /tmp:size=50g,exec \
            -u $UID:$(id -g) \
            -v "{YU_ROOT}:{YU_ROOT}" \
            -v "{STOIBER_ROOT}:{STOIBER_ROOT}:ro" \
            -v "{REFERENCE_DIR}:{REFERENCE_DIR}:ro" \
            --workdir "{CWD}" \
            {DOCKER_QUAST} \
            /bin/bash -c '
                set -eo pipefail

                mkdir -p "{params.outdir}"

                printf "Container hostname:\t"
                hostname

                printf "Start time:\t"
                date -Is

                echo "Assembler: {wildcards.assembler}"
                echo "Dataset: {wildcards.dataset}"
                echo "Assembly: {input.assembly}"
                echo "Reference: {input.reference}"
                echo "Threads: {threads}"
                echo "Memory: {resources.mem_gb} GB"

                echo
                echo "QUAST version:"
                quast.py --version

                echo
                echo "Running QUAST-LG"

                quast.py \
                    --large \
                    --threads {threads} \
                    --min-contig 500 \
                    --reference "{input.reference}" \
                    --output-dir "{params.outdir}" \
                    "{input.assembly}"

                echo
                printf "End time:\t"
                date -Is

                [[ -s "{output.quast_tsv}" ]] || {{
                    echo "ERROR: QUAST report.tsv is missing or empty"
                    exit 101
                }}

                grep -q "^N50" "{output.quast_tsv}" || {{
                    echo "ERROR: N50 not found in QUAST report"
                    exit 101
                }}
            ' >> "{log}" 2>&1
        """