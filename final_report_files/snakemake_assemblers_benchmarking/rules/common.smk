#####################################################################################
# SHARED ASSEMBLER WORKFLOW CONFIGURATION
#
# Responsibilities:
#   - Read sample metadata from samples.tsv.
#   - Keep machine-specific paths OUT of samples.tsv.
#   - Resolve source FASTQs from input_root automatically.
#   - Support local pre-extracted Chr21 testing and server WGS mode.
#   - Define shared assembler parameters and expected workflow outputs.
#####################################################################################

from pathlib import Path
import csv


# -----------------------------------------------------------------------------
# General workflow configuration
# -----------------------------------------------------------------------------

PROJECT_DIR = Path(workflow.basedir).resolve()

SAMPLE_SHEET = config.get("sample_sheet", "samples.tsv")
ACTIVE_ASSEMBLERS = config.get("active_assemblers", ["flye"])
GENOME_SIZE = config.get("genome_size", 46709983)
DEFAULT_THREADS = int(config.get("threads", 4))

INPUT_MODE = config.get("input_mode")

if INPUT_MODE not in {"preextracted", "whole_genome"}:
    raise ValueError(
        "Missing or invalid input_mode. "
        "Run the workflow with either config/local.yaml "
        "or config/server.yaml. Expected input_mode to be "
        "'preextracted' or 'whole_genome'."
    )


# -----------------------------------------------------------------------------
# Input root
# -----------------------------------------------------------------------------

RAW_INPUT_ROOT = config.get("input_root")

if not RAW_INPUT_ROOT:
    raise ValueError(
        "Missing input_root. "
        "Set input_root in config/local.yaml or config/server.yaml."
    )

INPUT_ROOT = Path(RAW_INPUT_ROOT).expanduser()

if not INPUT_ROOT.is_absolute():
    INPUT_ROOT = PROJECT_DIR / INPUT_ROOT

INPUT_ROOT = INPUT_ROOT.resolve(strict=False)


# -----------------------------------------------------------------------------
# Input/output naming
# -----------------------------------------------------------------------------

WHOLE_GENOME_PATTERN = config.get(
    "whole_genome_pattern",
    "{sample}.{technology}.30x.fastq.gz"
)

PREEXTRACTED_PATTERN = config.get(
    "preextracted_pattern",
    "{sample}.{technology}.fastq.gz"
)

CHR21_OUTPUT_DIR = config.get(
    "chr21_output_dir",
    "results/chr21"
)

# Separate validation markers prevent a successful local validation from being
# accidentally reused during a server run, or vice versa.
VALIDATION_OK = (
    f"results/logs/validate_inputs.{INPUT_MODE}.ok"
)


# -----------------------------------------------------------------------------
# Conda environments
# -----------------------------------------------------------------------------

# Absolute paths prevent included rule files from looking for rules/envs/.
ENV_FLYE = str(PROJECT_DIR / "envs" / "flye.yaml")
ENV_GOLDRUSH = str(PROJECT_DIR / "envs" / "goldrush.yaml")
ENV_NTLINK = str(PROJECT_DIR / "envs" / "ntlink.yaml")
ENV_VERKKO = str(PROJECT_DIR / "envs" / "verkko.yaml")
ENV_ASSESSMENT = str(PROJECT_DIR / "envs" / "assessment.yaml")


# -----------------------------------------------------------------------------
# Supported assemblers
# -----------------------------------------------------------------------------

SUPPORTED_ASSEMBLERS = {
    "flye",
    "goldrush",
    "ntlink",
    "verkko",
}

UNKNOWN_ASSEMBLERS = (
    set(ACTIVE_ASSEMBLERS)
    - SUPPORTED_ASSEMBLERS
)

if UNKNOWN_ASSEMBLERS:
    raise ValueError(
        "Unsupported assembler names in configuration: "
        f"{sorted(UNKNOWN_ASSEMBLERS)}"
    )


# -----------------------------------------------------------------------------
# Read portable sample table
# -----------------------------------------------------------------------------

def read_sample_sheet(path):
    """
    Read sample metadata.

    samples.tsv deliberately contains NO file paths.

    Required columns:
        sample
        technology

    Example:
        HG002   ont
        HG002   pb
    """

    path = Path(path)

    if not path.is_absolute():
        path = PROJECT_DIR / path

    if not path.exists():
        raise FileNotFoundError(
            f"Missing sample sheet: {path}"
        )

    rows = []
    seen_pairs = set()

    with path.open(newline="") as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t"
        )

        required_columns = {
            "sample",
            "technology",
        }

        observed_columns = set(
            reader.fieldnames or []
        )

        missing_columns = (
            required_columns
            - observed_columns
        )

        if missing_columns:
            raise ValueError(
                "Sample sheet is missing columns: "
                f"{sorted(missing_columns)}"
            )

        for row in reader:
            sample = row["sample"].strip()
            technology = row["technology"].strip()

            if not sample or not technology:
                raise ValueError(
                    f"Invalid empty value in sample sheet row: {row}"
                )

            if technology not in {"ont", "pb"}:
                raise ValueError(
                    "Unsupported sequencing technology: "
                    f"{technology}"
                )

            pair = (
                sample,
                technology,
            )

            if pair in seen_pairs:
                raise ValueError(
                    "Duplicate sample/technology combination: "
                    f"{sample} {technology}"
                )

            seen_pairs.add(pair)

            rows.append(
                {
                    "sample": sample,
                    "technology": technology,
                }
            )

    if not rows:
        raise ValueError(
            f"Sample sheet is empty: {path}"
        )

    return rows


SAMPLE_ROWS = read_sample_sheet(
    SAMPLE_SHEET
)

SAMPLES = list(
    dict.fromkeys(
        row["sample"]
        for row in SAMPLE_ROWS
    )
)

TECHNOLOGIES = list(
    dict.fromkeys(
        row["technology"]
        for row in SAMPLE_ROWS
    )
)

SAMPLE_TECH_PAIRS = list(
    dict.fromkeys(
        (
            row["sample"],
            row["technology"],
        )
        for row in SAMPLE_ROWS
    )
)

SAMPLE_TECH_SET = set(
    SAMPLE_TECH_PAIRS
)


# -----------------------------------------------------------------------------
# Source input path construction
# -----------------------------------------------------------------------------

def _check_sample_technology(
    sample,
    technology
):
    key = (
        sample,
        technology,
    )

    if key not in SAMPLE_TECH_SET:
        raise ValueError(
            "No dataset defined for "
            f"sample={sample}, "
            f"technology={technology}"
        )


def source_fastq_for_values(
    sample,
    technology
):
    """
    Return the ORIGINAL source FASTQ.

    Local mode:
        input_root/HG002.ont.fastq.gz

    Server mode:
        input_root/HG002.ont.30x.fastq.gz

    The filename is derived automatically.
    It is never stored in samples.tsv.
    """

    _check_sample_technology(
        sample,
        technology
    )

    if INPUT_MODE == "preextracted":
        pattern = PREEXTRACTED_PATTERN

    elif INPUT_MODE == "whole_genome":
        pattern = WHOLE_GENOME_PATTERN

    else:
        raise ValueError(
            f"Unsupported input_mode: {INPUT_MODE}"
        )

    filename = pattern.format(
        sample=sample,
        technology=technology,
    )

    return str(
        INPUT_ROOT / filename
    )


def source_fastq_for(wildcards):
    """
    Return the source FASTQ for a Snakemake wildcard object.

    This helper will be used by the future Chr21 extraction rule.
    """

    return source_fastq_for_values(
        wildcards.sample,
        wildcards.technology,
    )


# -----------------------------------------------------------------------------
# Canonical chromosome-21 paths
# -----------------------------------------------------------------------------

def generated_chr21_fastq_for_values(
    sample,
    technology
):
    """
    Canonical Chr21 FASTQ generated by server preprocessing.
    """

    _check_sample_technology(
        sample,
        technology
    )

    return (
        f"{CHR21_OUTPUT_DIR}/"
        f"{sample}.{technology}.fastq.gz"
    )


# -----------------------------------------------------------------------------
# Shared assembler input helpers
# -----------------------------------------------------------------------------

def assembler_fastq_for_values(
    sample,
    technology
):
    """
    Return the FASTQ that assemblers should consume.

    LOCAL:
        directly use Nicolas's pre-extracted real Chr21 30x reads.

    SERVER:
        use the canonical Chr21 FASTQ generated automatically from WGS.
    """

    _check_sample_technology(
        sample,
        technology
    )

    if INPUT_MODE == "preextracted":
        return source_fastq_for_values(
            sample,
            technology
        )

    if INPUT_MODE == "whole_genome":
        return generated_chr21_fastq_for_values(
            sample,
            technology
        )

    raise ValueError(
        f"Unsupported input_mode: {INPUT_MODE}"
    )


# Preserve the existing helper API used by Flye, GoldRush and ntLink.
def fastq_for(wildcards):
    return assembler_fastq_for_values(
        wildcards.sample,
        wildcards.technology,
    )


# Preserve the existing helper API used by Verkko.
def ont_fastq_for_sample(wildcards):
    return assembler_fastq_for_values(
        wildcards.sample,
        "ont",
    )


def pb_fastq_for_sample(wildcards):
    return assembler_fastq_for_values(
        wildcards.sample,
        "pb",
    )


# -----------------------------------------------------------------------------
# Technology-specific assembler parameters
# -----------------------------------------------------------------------------

def flye_read_option(wildcards):
    if wildcards.technology == "ont":
        return "--nano-raw"

    if wildcards.technology == "pb":
        return "--pacbio-hifi"

    raise ValueError(
        "Unsupported sequencing technology: "
        f"{wildcards.technology}"
    )


# -----------------------------------------------------------------------------
# Expected workflow outputs
# -----------------------------------------------------------------------------

EXPECTED_FINAL_OUTPUTS = []


# Flye
if "flye" in ACTIVE_ASSEMBLERS:
    EXPECTED_FINAL_OUTPUTS.extend(
        [
            (
                f"results/assemblies/"
                f"{sample}.{technology}.flye/"
                f"assembly.fasta"
            )
            for sample, technology
            in SAMPLE_TECH_PAIRS
        ]
    )


# GoldRush
if "goldrush" in ACTIVE_ASSEMBLERS:
    EXPECTED_FINAL_OUTPUTS.extend(
        [
            (
                f"results/smoke_test/goldrush/"
                f"{sample}.{technology}/assembly.fasta"
            )
            for sample, technology
            in SAMPLE_TECH_PAIRS
        ]
    )


# ntLink
if "ntlink" in ACTIVE_ASSEMBLERS:
    EXPECTED_FINAL_OUTPUTS.extend(
        [
            (
                f"results/smoke_test/ntlink/"
                f"{sample}.{technology}/assembly.fasta"
            )
            for sample, technology
            in SAMPLE_TECH_PAIRS
        ]
    )


# Verkko
if "verkko" in ACTIVE_ASSEMBLERS:
    EXPECTED_FINAL_OUTPUTS.extend(
        [
            (
                f"results/smoke_test/verkko/"
                f"{sample}/assembly.fasta"
            )
            for sample
            in SAMPLES
        ]
    )


# -----------------------------------------------------------------------------
# Wildcard constraints
# -----------------------------------------------------------------------------

wildcard_constraints:
    sample="|".join(SAMPLES),
    technology="|".join(TECHNOLOGIES)
