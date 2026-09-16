# ************************************************************************************************
#
# assembly_quality_busco.smk
#
# Whole-genome assembly completeness assessment with BUSCO.
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

ASSEMBLIES_ROOT = os.path.join(CWD, "assemblies")

if not os.path.isdir(ASSEMBLIES_ROOT):
    raise ValueError(
        "Assemblies directory does not exist: "
        + ASSEMBLIES_ROOT)

print("Current working directory: " + CWD)
print("Assemblies root: " + ASSEMBLIES_ROOT)


# ************************************************************************************************
# BUSCO configuration
# ************************************************************************************************

BUSCO_VERSION = "6.1.0"

DOCKER_BUSCO = (
    "quay.io/biocontainers/"
    "busco:6.1.0--pyhdfd78af_2")

print("BUSCO version: " + BUSCO_VERSION)
print("BUSCO Docker image: " + DOCKER_BUSCO)


# ************************************************************************************************
# BUSCO lineage configuration
# ************************************************************************************************

# ************************************************************************************************
# BUSCO lineage configuration
# ************************************************************************************************

DEFAULT_BUSCO_LINEAGE = (
    "/data/genmedbfx/ref/busco/lineages/"
    "primates_odb12.2")

RAW_BUSCO_LINEAGE = config.get(
    "busco_lineage",
    DEFAULT_BUSCO_LINEAGE)

BUSCO_LINEAGE = os.path.expanduser(RAW_BUSCO_LINEAGE)

if not os.path.isabs(BUSCO_LINEAGE):
    raise ValueError(
        "BUSCO lineage path must be absolute: "
        + BUSCO_LINEAGE)

BUSCO_LINEAGE = os.path.abspath(BUSCO_LINEAGE)

if not os.path.isdir(BUSCO_LINEAGE):
    raise ValueError(
        "BUSCO lineage directory does not exist: "
        + BUSCO_LINEAGE)

BUSCO_DATASET_CONFIG = os.path.join(
    BUSCO_LINEAGE,
    "dataset.cfg")

if not os.path.isfile(BUSCO_DATASET_CONFIG):
    raise ValueError(
        "BUSCO lineage does not contain dataset.cfg: "
        + BUSCO_LINEAGE)

BUSCO_LINEAGE_NAME = os.path.basename(
    BUSCO_LINEAGE.rstrip(os.sep))

print("BUSCO lineage: " + BUSCO_LINEAGE)
print("BUSCO lineage name: " + BUSCO_LINEAGE_NAME)


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

    # Verkko is a hybrid assembler with sample-level production outputs.
    if assembler == "verkko":
        return dataset in PRODUCTION_SAMPLES

    # Flye, GoldRush and ntLink use technology-specific 30x outputs.
    return ".30x" in dataset


ASSEMBLIES = sorted({
    (assembler, dataset)
    for assembler, dataset in zip(
        ASSEMBLERS_FOUND,
        DATASETS_FOUND)
    if is_production_assembly(assembler, dataset)})

if not ASSEMBLIES:
    raise ValueError(
        "No completed production assembly FASTAs were discovered under "
        "assemblies/{assembler}/{dataset}/assembly.fasta")

print("Discovered production assemblies:")

for assembler, dataset in ASSEMBLIES:
    print("  " + assembler + "\t" + dataset)


# ************************************************************************************************
# BUSCO targets
# ************************************************************************************************

OUTPUT = [
    (
        f"assembly_quality/busco/{assembler}/{dataset}/"
        f"busco/run_{BUSCO_LINEAGE_NAME}/short_summary.json"
    )
    for assembler, dataset in ASSEMBLIES]

print("BUSCO targets:")

for target in OUTPUT:
    print("  " + target)


rule all:
    input:
        OUTPUT


# ************************************************************************************************
# BUSCO genome-mode assessment
# ************************************************************************************************

rule busco_assembly:
    input:
        assembly = "assemblies/{assembler}/{dataset}/assembly.fasta"

    output:
        summary_json = (
            "assembly_quality/busco/"
            "{assembler}/{dataset}/"
            "busco/run_" + BUSCO_LINEAGE_NAME + "/short_summary.json")

    params:
        outdir = "assembly_quality/busco/{assembler}/{dataset}",
        run_name = "busco"

    log:
        "assembly_quality/busco/{assembler}/{dataset}/busco.log"

    threads:
        16

    resources:
        mem_gb = 64

    message:
        "Evaluating {wildcards.assembler} {wildcards.dataset} with BUSCO"

    shell:
        r"""
        set -eo pipefail

        mkdir -p "{params.outdir}"
        rm -rf "{params.outdir}/{params.run_name}"
        : > "{log}"

        echo "Assembler: {wildcards.assembler}" >> "{log}"
        echo "Dataset: {wildcards.dataset}" >> "{log}"
        echo "Assembly: {input.assembly}" >> "{log}"
        echo "BUSCO lineage: {BUSCO_LINEAGE}" >> "{log}"
        echo "Threads: {threads}" >> "{log}"
        echo "Memory: {resources.mem_gb} GB" >> "{log}"
        echo "Start time: $(date -Is)" >> "{log}"
        echo >> "{log}"

        docker run --rm \
            --cpus {threads} \
            -m {resources.mem_gb}g \
            --tmpfs /tmp:size=50g,exec \
            -u $UID:$(id -g) \
            -v "{CWD}:{CWD}" \
            -v "{BUSCO_LINEAGE}:{BUSCO_LINEAGE}:ro" \
            --workdir "{CWD}" \
            {DOCKER_BUSCO} \
            /bin/bash -c '
                set -eo pipefail

                printf "Container hostname:\t"
                hostname

                echo
                echo "BUSCO version:"
                busco --version

                echo
                echo "Running BUSCO"

                busco \
                    --in "{input.assembly}" \
                    --mode genome \
                    --lineage_dataset "{BUSCO_LINEAGE}" \
                    --cpu {threads} \
                    --out "{params.run_name}" \
                    --out_path "{params.outdir}" \
                    --offline \
                    --miniprot \
                    --opt-out-run-stats
            ' >> "{log}" 2>&1

        echo >> "{log}"
        echo "End time: $(date -Is)" >> "{log}"

        [[ -s "{output.summary_json}" ]] || {{
            echo "ERROR: BUSCO short_summary.json is missing or empty" \
                >> "{log}"
            exit 101
        }}
        """