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

CWD = os.getcwd()

print("Current working directory: " + CWD)


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

RAW_BUSCO_LINEAGE = config.get("busco_lineage")

if not RAW_BUSCO_LINEAGE:
    raise ValueError(
        "BUSCO lineage directory is required. "
        "Provide it with --config busco_lineage=/path/to/primates_odb12.2")

BUSCO_LINEAGE = os.path.expanduser(RAW_BUSCO_LINEAGE)

if not os.path.isabs(BUSCO_LINEAGE):
    raise ValueError(
        "BUSCO lineage path must be absolute after expansion: "
        + BUSCO_LINEAGE)

BUSCO_LINEAGE = os.path.abspath(BUSCO_LINEAGE)

if not os.path.isdir(BUSCO_LINEAGE):
    raise ValueError(
        "BUSCO lineage directory does not exist: " + BUSCO_LINEAGE)

BUSCO_DATASET_CONFIG = os.path.join(
    BUSCO_LINEAGE,
    "dataset.cfg")

if not os.path.isfile(BUSCO_DATASET_CONFIG):
    raise ValueError(
        "BUSCO lineage does not contain dataset.cfg: " + BUSCO_LINEAGE)

BUSCO_LINEAGE_DIR = os.path.dirname(BUSCO_LINEAGE)
BUSCO_LINEAGE_NAME = os.path.basename(
    BUSCO_LINEAGE.rstrip(os.sep))

print("BUSCO lineage: " + BUSCO_LINEAGE)
print("BUSCO lineage name: " + BUSCO_LINEAGE_NAME)


# ************************************************************************************************
# Discover completed assembly FASTAs
# ************************************************************************************************

ASSEMBLERS_FOUND, DATASETS_FOUND = glob_wildcards(
    CWD + r"/assemblies/{assembler}/{dataset}/assembly.fasta")

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
    "smoke",
)

ASSEMBLIES = sorted(
    {
        (assembler, dataset)
        for assembler, dataset in zip(
            ASSEMBLERS_FOUND,
            DATASETS_FOUND
        )
        if assembler.lower() in ALLOWED_OUTPUTS
        and ".30x" in dataset.lower()
        and not any(
            marker in dataset.lower()
            for marker in TEST_MARKERS)})

if not ASSEMBLIES:
    raise ValueError(
        "No completed production assembly FASTAs were discovered under "
        "assemblies/{assembler}/{dataset}/assembly.fasta")

print("Discovered assemblies:")

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

        echo "Assembler: {wildcards.assembler}" > "{log}"
        echo "Dataset: {wildcards.dataset}" >> "{log}"
        echo "Assembly: {input.assembly}" >> "{log}"
        echo "BUSCO lineage: {BUSCO_LINEAGE}" >> "{log}"
        echo "Threads: {threads}" >> "{log}"
        echo "Memory: {resources.mem_gb} GB" >> "{log}"
        echo "Start time: $(date -Is)" >> "{log}"
        echo >> "{log}"

        docker run --rm --cpus {threads} -m {resources.mem_gb}g --tmpfs /tmp:size=50g,exec -u $UID:$(id -g) -v "{CWD}:{CWD}" -v "{BUSCO_LINEAGE_DIR}:{BUSCO_LINEAGE_DIR}:ro" --workdir "{CWD}" {DOCKER_BUSCO} /bin/bash -c '
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
            echo "ERROR: BUSCO short_summary.json is missing or empty" >> "{log}"
            exit 101
        }}
        """