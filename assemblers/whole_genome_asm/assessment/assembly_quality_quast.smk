#
# assembly_quality_quast.smk
#
# Whole-genome assembly quality assessment with QUAST-LG.
#
# Evaluates existing Flye, GoldRush, Verkko and ntLink outputs.
# ntLink is treated as a scaffolder, not as an assembler.
#
# Constitution: Articles I, II, III, V, VI and VII
#
# ************************************************************************************************

import os


# ************************************************************************************************
# Repository root
# ************************************************************************************************

CWD = os.path.abspath(os.getcwd())

print("Current working directory: " + CWD)


# ************************************************************************************************
# QUAST configuration
# ************************************************************************************************

QUAST_VERSION = "5.3.0"

DOCKER_QUAST = (
    "quay.io/biocontainers/"
    "quast:5.3.0--py313pl5321h5ca1c30_2")

print("QUAST version: " + QUAST_VERSION)
print("QUAST Docker image: " + DOCKER_QUAST)


# ************************************************************************************************
# Reference configuration
# ************************************************************************************************

RAW_REFERENCE = config.get("reference")

if not RAW_REFERENCE:
    raise ValueError(
        "Missing reference genome. Provide it with "
        "--config reference=/absolute/path/to/reference.fasta")


REFERENCE = os.path.abspath(os.path.expanduser(RAW_REFERENCE))

if not os.path.isfile(REFERENCE):
    raise ValueError(
        "Reference genome does not exist: " + REFERENCE)

REFERENCE_DIR = os.path.dirname(REFERENCE)

print("Reference genome: " + REFERENCE)


# ************************************************************************************************
# Assembly-root configuration
# ************************************************************************************************

ASSEMBLIES_ROOT = os.path.join(CWD, "assemblies")

if not os.path.isdir(ASSEMBLIES_ROOT):
    raise ValueError(
        "Assemblies directory does not exist: " + ASSEMBLIES_ROOT)

print("Assemblies root: " + ASSEMBLIES_ROOT)


# ************************************************************************************************
# Discover completed assembly FASTAs
# ************************************************************************************************

ASSEMBLERS_FOUND, DATASETS_FOUND = glob_wildcards(
    ASSEMBLIES_ROOT
    + r"/{assembler}/{dataset}/assembly.fasta")

print("Raw assemblers found:", ASSEMBLERS_FOUND)
print("Raw datasets found:", DATASETS_FOUND)


ALLOWED_OUTPUTS = {
    "flye",
    "goldrush",
    "verkko",
    "ntlink",}

TEST_MARKERS = (
    ".1k",
    ".chr21",
    ".localtest",
    "smoke",)

PRODUCTION_SAMPLES = {
    "hg002",
    "hg003",
    "hg004",}


def is_production_assembly(assembler, dataset):

    assembler = assembler.lower()
    dataset = dataset.lower()

    if assembler not in ALLOWED_OUTPUTS:
        return False

    if any(marker in dataset for marker in TEST_MARKERS):
        return False

    # Verkko is hybrid and uses sample-level names:
    # HG002, HG003, HG004
    if assembler == "verkko":
        return dataset in PRODUCTION_SAMPLES

    # Flye / GoldRush / ntLink use technology-specific 30x datasets
    return ".30x" in dataset


ASSEMBLIES = sorted(
    {
        (assembler, dataset)
        for assembler, dataset in zip(
            ASSEMBLERS_FOUND,
            DATASETS_FOUND
        )
        if is_production_assembly(
            assembler,
            dataset)})


if not ASSEMBLIES:
    raise ValueError(
        "No completed production assembly FASTAs were discovered under "
        + ASSEMBLIES_ROOT
        + "/{assembler}/{dataset}/assembly.fasta")


print("Discovered assemblies:")

for assembler, dataset in ASSEMBLIES:
    print("  " + assembler + "\t" + dataset)


# ************************************************************************************************
# QUAST targets
# ************************************************************************************************

OUTPUT = [
    f"assembly_quality/quast/{assembler}/{dataset}/report.tsv"
    for assembler, dataset in ASSEMBLIES]


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
        assembly = "assemblies/{assembler}/{dataset}/assembly.fasta",
        reference = REFERENCE

    output: quast_tsv = "assembly_quality/quast/{assembler}/{dataset}/report.tsv"

    params: outdir = "assembly_quality/quast/{assembler}/{dataset}"

    log: "assembly_quality/quast/{assembler}/{dataset}/quast.log"

    threads: 16

    resources: mem_gb = 128

    message: "Evaluating {wildcards.assembler} {wildcards.dataset} with QUAST-LG"

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
            -v "{CWD}:{CWD}" \
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
