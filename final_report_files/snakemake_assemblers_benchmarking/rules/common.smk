#####################################################################################
# SHARED ASSEMBLER WORKFLOW CONFIGURATION
# This module reads the sample table, validates configuration values, defines
# helper functions, and generates the expected outputs used by the workflow.
#####################################################################################

# Import packages used to verify paths and read tab-separated sample information
from pathlib import Path
import csv


# -----------------------------------------------------------------------------
# General workflow configuration
# -----------------------------------------------------------------------------

# Define the sample table and general assembler settings
SAMPLE_SHEET = config.get("sample_sheet", "samples.tsv")
ACTIVE_ASSEMBLERS = config.get("active_assemblers", ["flye"])
GENOME_SIZE = config.get("genome_size", "3.1g")
DEFAULT_THREADS = int(config.get("threads", 4))

# Marker file created after all required inputs have been validated
VALIDATION_OK = "results/logs/validate_inputs.ok"

# Conda environments used by the individual assembler and assessment modules
ENV_FLYE = "envs/flye.yaml"
ENV_GOLDRUSH = "envs/goldrush.yaml"
ENV_NTLINK = "envs/ntlink.yaml"
ENV_VERKKO = "envs/verkko.yaml"
ENV_ASSESSMENT = "envs/assessment.yaml"

# Define all assembler and scaffolding names accepted by this workflow
SUPPORTED_ASSEMBLERS = {
    "flye",
    "goldrush",
    "ntlink",
    "verkko",
}

# Stop immediately if config.yaml contains an unknown tool name
UNKNOWN_ASSEMBLERS = set(ACTIVE_ASSEMBLERS) - SUPPORTED_ASSEMBLERS

if UNKNOWN_ASSEMBLERS:
    raise ValueError(
        f"Unsupported assembler names in config.yaml: "
        f"{sorted(UNKNOWN_ASSEMBLERS)}"
    )


# -----------------------------------------------------------------------------
# Read the sample table
# -----------------------------------------------------------------------------

# Read samples.tsv, verify its required columns, reject missing values, and
# return the sample name, sequencing technology, and FASTQ path for each row
def read_sample_sheet(path):
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(
            f"Missing sample sheet: {path}"
        )

    rows = []

    with path.open(newline="") as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t"
        )

        required_columns = {
            "sample",
            "technology",
            "fastq",
        }

        observed_columns = set(
            reader.fieldnames or []
        )

        missing_columns = (
            required_columns - observed_columns
        )

        if missing_columns:
            raise ValueError(
                f"Sample sheet is missing columns: "
                f"{sorted(missing_columns)}"
            )

        for row in reader:
            sample = row["sample"].strip()
            technology = row["technology"].strip()
            fastq = row["fastq"].strip()

            if not sample or not technology or not fastq:
                raise ValueError(
                    f"Invalid empty value in sample sheet row: {row}"
                )

            if technology not in {"ont", "pb"}:
                raise ValueError(
                    f"Unsupported sequencing technology: {technology}"
                )

            rows.append(
                {
                    "sample": sample,
                    "technology": technology,
                    "fastq": fastq,
                }
            )

    if not rows:
        raise ValueError(
            f"Sample sheet is empty: {path}"
        )

    return rows


# Transform samples.tsv into structures used to generate Snakemake jobs
SAMPLE_ROWS = read_sample_sheet(SAMPLE_SHEET)

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

FASTQ_BY_SAMPLE_TECH = {
    (
        row["sample"],
        row["technology"],
    ): row["fastq"]
    for row in SAMPLE_ROWS
}


# -----------------------------------------------------------------------------
# Shared input helper functions
# -----------------------------------------------------------------------------

# Return the FASTQ associated with one sample and sequencing technology
def fastq_for(wildcards):
    key = (
        wildcards.sample,
        wildcards.technology,
    )

    if key not in FASTQ_BY_SAMPLE_TECH:
        raise ValueError(
            f"No FASTQ defined for "
            f"sample={wildcards.sample}, "
            f"technology={wildcards.technology}"
        )

    return FASTQ_BY_SAMPLE_TECH[key]


# Return the ONT FASTQ associated with one sample
def ont_fastq_for_sample(wildcards):
    key = (
        wildcards.sample,
        "ont",
    )

    if key not in FASTQ_BY_SAMPLE_TECH:
        raise ValueError(
            f"No ONT FASTQ defined for "
            f"sample={wildcards.sample}"
        )

    return FASTQ_BY_SAMPLE_TECH[key]


# Return the PacBio HiFi FASTQ associated with one sample
def pb_fastq_for_sample(wildcards):
    key = (
        wildcards.sample,
        "pb",
    )

    if key not in FASTQ_BY_SAMPLE_TECH:
        raise ValueError(
            f"No PacBio HiFi FASTQ defined for "
            f"sample={wildcards.sample}"
        )

    return FASTQ_BY_SAMPLE_TECH[key]


# -----------------------------------------------------------------------------
# Technology-specific parameter functions
# -----------------------------------------------------------------------------

# Select the correct Flye input option for ONT or PacBio HiFi reads
def flye_read_option(wildcards):
    if wildcards.technology == "ont":
        return "--nano-raw"

    if wildcards.technology == "pb":
        return "--pacbio-hifi"

    raise ValueError(
        f"Unsupported sequencing technology: "
        f"{wildcards.technology}"
    )


# -----------------------------------------------------------------------------
# Expected workflow outputs
# -----------------------------------------------------------------------------

# Build the final output collection according to the tools enabled in config.yaml
EXPECTED_FINAL_OUTPUTS = []

# Flye generates one independent assembly for each sample and technology
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

# GoldRush outputs will be added after its command and installation are validated

# ntLink outputs will be added after selecting which draft assemblies it will
# scaffold using the corresponding long-read datasets

# Verkko outputs will be added after validating its hybrid ONT and PacBio inputs


# -----------------------------------------------------------------------------
# Wildcard constraints
# -----------------------------------------------------------------------------

# Restrict wildcards to the samples and technologies listed in samples.tsv
wildcard_constraints:
    sample="|".join(SAMPLES),
    technology="|".join(TECHNOLOGIES)
