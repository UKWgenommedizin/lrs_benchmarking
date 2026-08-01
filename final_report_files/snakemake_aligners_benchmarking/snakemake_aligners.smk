#####################################################################################
# ALIGNMENT PIPELINE: MINIMAP2, PBMM2, VACMAP, AND VG GIRAFFE
# This workflow benchmarks long-read aligners using the same sample information, reference genome, and output structure.
#####################################################################################

# Import packages to verify file existence and read tab-separated sample tables
from pathlib import Path
import csv

configfile: "config.yaml"

# -----------------------------------------------------------------------------
# Configuration of the environment
# -----------------------------------------------------------------------------

# Variables containing the sample table and reference genome
SAMPLE_SHEET = config.get("sample_sheet", "samples.tsv")
REFERENCE = config.get("reference", "data/reference/genome.fasta")

# Define which aligners will be included in the benchmark
ACTIVE_ALIGNERS = config.get(
    "active_aligners",
    ["minimap2", "pbmm2", "vacmap", "vg_giraffe"]
)

# Store the command used to run VACmap
VACMAP_BIN = config.get("vacmap_bin", "conda run -n vacmap_env vacmap")

# Parameters required by vg Giraffe because it uses graph indexes instead of a linear FASTA reference like the other aligners
VG_GBZ = config.get("vg_gbz", "data/graph/graph.gbz")  # Store the path to the GBZ graph file
VG_MIN = config.get("vg_min", "data/graph/graph.min")  # Store the path to the minimizer index
VG_DIST = config.get("vg_dist", "data/graph/graph.dist")  # Store the path to the distance index used by vg Giraffe
VG_ZIPCODES = config.get("vg_zipcodes", "")  # Optional zipcodes file for vg Giraffe

# Marker file created after the environment validation step passes; this supports reproducibility
VALIDATION_OK = "results/logs/validate_environment.ok"

# Conda environment used for alignment, statistics, and plotting rules
ENV_ALIGNMENT = "envs/alignment.yaml"

# Conda environment used specifically for vg Giraffe
ENV_VG = "envs/vg.yaml"


# -----------------------------------------------------------------------------
# Read the sample sheet containing HG002, HG003, and HG004, with ONT and PacBio HiFi technologies
# -----------------------------------------------------------------------------

# Read samples.tsv, verify required columns and missing values, and return sample, technology, and FASTQ paths
def read_sample_sheet(path): 
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(f"Missing sample sheet: {path}")

    rows = []

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

    rows = []

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required_columns = {"sample", "technology", "fastq"}
        observed_columns = set(reader.fieldnames or [])

        missing_columns = required_columns - observed_columns
        if missing_columns:
            raise ValueError(
                f"Sample sheet is missing columns: {sorted(missing_columns)}"
            )

        for row in reader:
            sample = row["sample"].strip()
            technology = row["technology"].strip()
            fastq = row["fastq"].strip()

            if not sample or not technology or not fastq:
                raise ValueError(f"Invalid empty value in sample sheet row: {row}")

            rows.append(
                {
                    "sample": sample,
                    "technology": technology,
                    "fastq": fastq,
                }
            )

    if not rows:
        raise ValueError(f"Sample sheet is empty: {path}")

    return rows


# Transform the sample sheet into lists and dictionaries to generate Snakemake jobs automatically
SAMPLE_ROWS = read_sample_sheet(SAMPLE_SHEET)

# Define sample and technology lists from the sample sheet
SAMPLES = list(dict.fromkeys(row["sample"] for row in SAMPLE_ROWS))
TECHNOLOGIES = list(dict.fromkeys(row["technology"] for row in SAMPLE_ROWS))

SAMPLE_TECH_PAIRS = list(
    dict.fromkeys((row["sample"], row["technology"]) for row in SAMPLE_ROWS)
)

FASTQ_BY_SAMPLE_TECH = {
    (row["sample"], row["technology"]): row["fastq"]
    for row in SAMPLE_ROWS
}

# Function used to define which FASTQ file belongs to each sample and technology combination
def fastq_for(wildcards):
    key = (wildcards.sample, wildcards.technology)

    if key not in FASTQ_BY_SAMPLE_TECH:
        raise ValueError(
            f"No FASTQ defined for sample={wildcards.sample}, "
            f"technology={wildcards.technology}"
        )

    return FASTQ_BY_SAMPLE_TECH[key]


# Define the expected samtools stats files for each sample, technology, and aligner combination
EXPECTED_STATS = [
    f"results/stats/{sample}.{technology}.{aligner}.stats.txt"
    for sample, technology in SAMPLE_TECH_PAIRS
    for aligner in ACTIVE_ALIGNERS
]

# Define the expected verification marker for every requested alignment
EXPECTED_ALIGNMENT_CHECKS = [
    f"results/checkpoints/{sample}.{technology}.{aligner}.verified"
    for sample, technology in SAMPLE_TECH_PAIRS
    for aligner in ACTIVE_ALIGNERS
]

# Define the expected flagstat files for every requested alignment
EXPECTED_FLAGSTATS = [
    f"results/stats/{sample}.{technology}.{aligner}.flagstat.txt"
    for sample, technology in SAMPLE_TECH_PAIRS
    for aligner in ACTIVE_ALIGNERS
]

# Define the expected idxstat files for every requested alignment
EXPECTED_IDXSTATS = [
    f"results/stats/{sample}.{technology}.{aligner}.idxstat.txt"
    for sample, technology in SAMPLE_TECH_PAIRS
    for aligner in ACTIVE_ALIGNERS
]

# Restrict wildcards to valid sample names, technologies, and aligners
wildcard_constraints:
    sample="|".join(SAMPLES),
    technology="|".join(TECHNOLOGIES),
    aligner="|".join(ACTIVE_ALIGNERS)


# -----------------------------------------------------------------------------
# Preset / mode functions
# -----------------------------------------------------------------------------

# Select the correct preset or mode based on the sequencing technology
def minimap2_preset(wildcards):
    if wildcards.technology == "ont":
        return "map-ont"
    if wildcards.technology == "pb":
        return "map-hifi"
    raise ValueError(f"Unknown technology: {wildcards.technology}")


def pbmm2_preset(wildcards):
    if wildcards.technology == "ont":
        return "SUBREAD"
    if wildcards.technology == "pb":
        return "HIFI"
    raise ValueError(f"Unknown technology: {wildcards.technology}")


def vacmap_mode(wildcards):
    if wildcards.technology == "ont":
        return "H"
    if wildcards.technology == "pb":
        return "L"
    raise ValueError(f"Unknown technology: {wildcards.technology}")


def vg_giraffe_preset(wildcards):
    if wildcards.technology == "ont":
        return "r10"
    if wildcards.technology == "pb":
        return "hifi"
    raise ValueError(f"Unknown technology: {wildcards.technology}")


def vg_zipcodes_arg(wildcards):
    if VG_ZIPCODES:
        return f"-z {VG_ZIPCODES}"
    return ""


# -----------------------------------------------------------------------------
# Final target rule
# -----------------------------------------------------------------------------

# Define the expected final files produced by the workflow
rule all:
    input:
        VALIDATION_OK,
        "results/tables/alignment_summary.tsv",
        "results/tables/quality_report.tsv",
        "results/tables/plot_data_clean.tsv",
        "results/tables/plot_data_pivot.tsv",
        "results/tables/plot_arrays.npz",
        "results/tables/paired_alignment_tests.tsv",
        "results/figures/final_alignment_error_rate.png",
        "results/figures/final_alignment_error_rate.pdf",
        "results/figures/mapped_bases.png",
        "results/figures/mapped_bases.pdf",
        "results/figures/mapped_bases_cigar.png",
        "results/figures/mapped_bases_cigar.pdf",
        "results/checkpoints/workflow_complete.ok"


# -----------------------------------------------------------------------------
# Check environment and required files before running the workflow
# -----------------------------------------------------------------------------

# Check the required input files, Conda environments, reference genome, FASTQ files, graph indexes, and Python scripts before running the benchmark
# Prevent the workflow from failing in later steps due to missing files or incomplete setup
rule validate_environment:
    output:
        VALIDATION_OK
    run:
        Path(output[0]).parent.mkdir(parents=True, exist_ok=True)

        errors = []

        for required_file in [SAMPLE_SHEET, ENV_ALIGNMENT, ENV_VG, REFERENCE]:
            if not Path(required_file).exists():
                errors.append(f"Missing required file: {required_file}")

        for row in SAMPLE_ROWS:
            fastq = Path(row["fastq"])
            if not fastq.exists():
                errors.append(f"Missing FASTQ file: {fastq}")

        if "vg_giraffe" in ACTIVE_ALIGNERS:
            for graph_file in [VG_GBZ, VG_MIN, VG_DIST]:
                if not Path(graph_file).exists():
                    errors.append(f"Missing vg Giraffe graph/index file: {graph_file}")

            if VG_ZIPCODES and not Path(VG_ZIPCODES).exists():
                errors.append(f"Missing vg Giraffe zipcodes file: {VG_ZIPCODES}")

        required_scripts = [
            "scripts/alignment_summary.py",
            "scripts/quality_check.py",
            "scripts/prepare_plot_data.py",
            "scripts/paired_t_tests.py",
            "scripts/alignment_metrics.py",
        ]

        for script in required_scripts:
            if not Path(script).exists():
                errors.append(f"Missing Python script: {script}")

        if errors:
            message = "\n".join(errors)
            raise RuntimeError(
                "\nENVIRONMENT VALIDATION FAILED\n"
                "Fix these problems before running the workflow:\n\n"
                f"{message}\n"
            )

        Path(output[0]).write_text("Environment validation passed.\n")


# -----------------------------------------------------------------------------
# Rule 1: minimap2 alignment
# -----------------------------------------------------------------------------

# Each aligner uses a FASTQ file and the required reference input, creates a sorted BAM file, and indexes it with samtools
rule align_minimap2:
    input:
        validation=VALIDATION_OK,
        fastq=fastq_for,
        ref=REFERENCE
    output:
        bam="results/bam/{sample}.{technology}.minimap2.sorted.bam",
        bai="results/bam/{sample}.{technology}.minimap2.sorted.bam.bai"
    conda:
        ENV_ALIGNMENT
    params:
        preset=minimap2_preset
    threads: 4
    shell:
        r"""
        mkdir -p results/bam

        minimap2 -ax {params.preset} -t {threads} {input.ref} {input.fastq} \
        | samtools sort -@ {threads} -o {output.bam}

        samtools index {output.bam}
        """


# -----------------------------------------------------------------------------
# Rule 2: pbmm2 alignment
# -----------------------------------------------------------------------------


rule align_pbmm2:
    input:
        validation=VALIDATION_OK,
        fastq=fastq_for,
        ref=REFERENCE
    output:
        bam="results/bam/{sample}.{technology}.pbmm2.sorted.bam",
        bai="results/bam/{sample}.{technology}.pbmm2.sorted.bam.bai"
    conda:
        ENV_ALIGNMENT
    params:
        preset=pbmm2_preset
    threads: 4
    shell:
        r"""
        mkdir -p results/bam

        pbmm2 align \
            --preset {params.preset} \
            --sort \
            -j {threads} \
            {input.ref} \
            {input.fastq} \
            {output.bam}

        samtools index {output.bam}
        """


# -----------------------------------------------------------------------------
# Rule 3: VACmap alignment
# -----------------------------------------------------------------------------

rule align_vacmap:
    input:
        validation=VALIDATION_OK,
        fastq=fastq_for,
        ref=REFERENCE
    output:
        bam="results/bam/{sample}.{technology}.vacmap.sorted.bam",
        bai="results/bam/{sample}.{technology}.vacmap.sorted.bam.bai"
    conda:
        ENV_ALIGNMENT
    params:
        mode=vacmap_mode,
        vacmap_bin=VACMAP_BIN
    threads: 4
    shell:
        r"""
        mkdir -p results/bam

        {params.vacmap_bin} \
            -ref {input.ref} \
            -read {input.fastq} \
            -mode {params.mode} \
            -t {threads} \
            --force \
            -o {output.bam}

        samtools index {output.bam}
        """


# -----------------------------------------------------------------------------
# Rule 4: vg Giraffe alignment
# -----------------------------------------------------------------------------

# Unlike the other aligners, vg Giraffe uses graph indexes instead of a linear reference genome: GBZ graph, minimizer index, and distance index
rule align_vg_giraffe:
    input:
        validation=VALIDATION_OK,
        fastq=fastq_for,
        gbz=VG_GBZ,
        minimizer=VG_MIN,
        dist=VG_DIST
    output:
        bam="results/bam/{sample}.{technology}.vg_giraffe.sorted.bam",
        bai="results/bam/{sample}.{technology}.vg_giraffe.sorted.bam.bai"
    conda:
        ENV_VG
    params:
        preset=vg_giraffe_preset,
        zipcodes=vg_zipcodes_arg
    threads: 4
    shell:
        r"""
        mkdir -p results/bam

        vg giraffe \
            -Z {input.gbz} \
            -m {input.minimizer} \
            -d {input.dist} \
            {params.zipcodes} \
            -f {input.fastq} \
            -t {threads} \
            -b {params.preset} \
            -o BAM \
        | samtools sort -@ {threads} -o {output.bam}

        samtools index {output.bam}
        """


# -----------------------------------------------------------------------------
# Checkpoint 1: Verify each completed alignment
# -----------------------------------------------------------------------------

# Verify that each BAM and BAI file exists, is non-empty, has a readable
# header, passes samtools quickcheck, and can be queried through its index.
rule verify_alignment:
    input:
        bam="results/bam/{sample}.{technology}.{aligner}.sorted.bam",
        bai="results/bam/{sample}.{technology}.{aligner}.sorted.bam.bai"
    output:
        verified="results/checkpoints/{sample}.{technology}.{aligner}.verified"
    log:
        "results/logs/{sample}.{technology}.{aligner}.verification.log"
    conda:
        ENV_ALIGNMENT
    shell:
        r"""
        set -euo pipefail

        mkdir -p results/checkpoints
        mkdir -p results/logs

        {{
            echo "ALIGNMENT VERIFICATION"
            echo "Sample: {wildcards.sample}"
            echo "Technology: {wildcards.technology}"
            echo "Aligner: {wildcards.aligner}"
            echo "BAM: {input.bam}"
            echo "BAI: {input.bai}"
            echo "Started: $(date --iso-8601=seconds)"
        }} > {log}

        test -s {input.bam}
        test -s {input.bai}

        samtools quickcheck -v {input.bam} >> {log} 2>&1
        samtools view -H {input.bam} > /dev/null
        samtools idxstats {input.bam} > /dev/null

        {{
            echo "Finished: $(date --iso-8601=seconds)"
            echo "STATUS: PASSED"
        }} >> {log}

        {{
            echo "ALIGNMENT VERIFIED"
            echo "Sample: {wildcards.sample}"
            echo "Technology: {wildcards.technology}"
            echo "Aligner: {wildcards.aligner}"
            echo "BAM: {input.bam}"
            echo "Verified: $(date --iso-8601=seconds)"
        }} > {output.verified}
        """


# -----------------------------------------------------------------------------
# Checkpoint 2: Confirm that every requested alignment passed verification
# -----------------------------------------------------------------------------

# Create one stage marker only after all requested BAM files have passed the
# individual verification rule.
rule checkpoint_alignments_complete:
    input:
        verified=EXPECTED_ALIGNMENT_CHECKS
    output:
        checkpoint="results/checkpoints/alignments_complete.ok"
    run:
        checkpoint_path = Path(output.checkpoint)
        checkpoint_path.parent.mkdir(parents=True, exist_ok=True)

        invalid_markers = [
            str(marker)
            for marker in input.verified
            if not Path(marker).exists() or Path(marker).stat().st_size == 0
        ]

        if invalid_markers:
            raise RuntimeError(
                "Alignment checkpoint failed. Missing or empty verification "
                "markers:\n" + "\n".join(invalid_markers)
            )

        checkpoint_path.write_text(
            "ALIGNMENT CHECKPOINT PASSED\n"
            f"Verified alignments: {len(input.verified)}\n"
            "All requested BAM and BAI files passed samtools verification.\n"
        )


# -----------------------------------------------------------------------------
# Rule 5: Samtools quality statistics
# -----------------------------------------------------------------------------

# Calculate alignment statistics for each BAM file using samtools flagstat, idxstat, and stats
rule samtools_quality:
    input:
        validation=VALIDATION_OK,
        alignment_checkpoint="results/checkpoints/alignments_complete.ok",
        verified="results/checkpoints/{sample}.{technology}.{aligner}.verified",
        bam="results/bam/{sample}.{technology}.{aligner}.sorted.bam",
        bai="results/bam/{sample}.{technology}.{aligner}.sorted.bam.bai"
    output:
        flagstat="results/stats/{sample}.{technology}.{aligner}.flagstat.txt",
        idxstat="results/stats/{sample}.{technology}.{aligner}.idxstat.txt",
        stats="results/stats/{sample}.{technology}.{aligner}.stats.txt"
    conda:
        ENV_ALIGNMENT
    shell:
        r"""
        mkdir -p results/stats

        samtools flagstat {input.bam} > {output.flagstat}
        samtools idxstats {input.bam} > {output.idxstat}
        samtools stats {input.bam} > {output.stats}
        """


# -----------------------------------------------------------------------------
# Checkpoint 3: Verify all samtools statistics files
# -----------------------------------------------------------------------------

# Verify that flagstat, idxstat, and samtools stats outputs exist and contain
# data before the workflow creates the combined alignment summary.
rule checkpoint_statistics_complete:
    input:
        alignment_checkpoint="results/checkpoints/alignments_complete.ok",
        flagstats=EXPECTED_FLAGSTATS,
        idxstats=EXPECTED_IDXSTATS,
        stats=EXPECTED_STATS
    output:
        checkpoint="results/checkpoints/statistics_complete.ok"
    run:
        checkpoint_path = Path(output.checkpoint)
        checkpoint_path.parent.mkdir(parents=True, exist_ok=True)

        files_to_check = (
            list(input.flagstats)
            + list(input.idxstats)
            + list(input.stats)
        )

        invalid_files = [
            str(file_path)
            for file_path in files_to_check
            if not Path(file_path).exists()
            or Path(file_path).stat().st_size == 0
        ]

        if invalid_files:
            raise RuntimeError(
                "Statistics checkpoint failed. Missing or empty files:\n"
                + "\n".join(invalid_files)
            )

        checkpoint_path.write_text(
            "STATISTICS CHECKPOINT PASSED\n"
            f"Flagstat files: {len(input.flagstats)}\n"
            f"Idxstat files: {len(input.idxstats)}\n"
            f"Samtools stats files: {len(input.stats)}\n"
            f"Total statistics files: {len(files_to_check)}\n"
        )


# -----------------------------------------------------------------------------
# Rule 6: Build alignment summary table
# -----------------------------------------------------------------------------

# Combine all samtools stats files into one alignment summary table
rule alignment_summary:
    input:
        statistics_checkpoint="results/checkpoints/statistics_complete.ok",
        stats=EXPECTED_STATS
    output:
        summary="results/tables/alignment_summary.tsv"
    conda:
        ENV_ALIGNMENT
    shell:
        r"""
        mkdir -p results/tables

        python scripts/alignment_summary.py \
            --stats {input.stats} \
            --out {output.summary}
        """


# -----------------------------------------------------------------------------
# Rule 7: Quality check of the alignment summary
# -----------------------------------------------------------------------------

# Check the alignment summary for duplicated conditions, missing values, inconsistent labels, and other data-quality problems
rule quality_check:
    input:
        summary="results/tables/alignment_summary.tsv"
    output:
        quality_report="results/tables/quality_report.tsv"
    conda:
        ENV_ALIGNMENT
    shell:
        r"""
        python scripts/quality_check.py \
            --input {input.summary} \
            --out {output.quality_report}
        """


# -----------------------------------------------------------------------------
# Rule 8: Prepare data for plotting
# -----------------------------------------------------------------------------

# Clean and reshape the summary table into plot-ready data structures
rule prepare_plot_data:
    input:
        summary="results/tables/alignment_summary.tsv",
        quality_report="results/tables/quality_report.tsv"
    output:
        clean="results/tables/plot_data_clean.tsv",
        pivot="results/tables/plot_data_pivot.tsv",
        arrays="results/tables/plot_arrays.npz"
    conda:
        ENV_ALIGNMENT
    shell:
        r"""
        python scripts/prepare_plot_data.py \
            --input {input.summary} \
            --quality-report {input.quality_report} \
            --cleaned-data {output.clean} \
            --pivoted-data {output.pivot} \
            --arrays {output.arrays}
        """


# -----------------------------------------------------------------------------
# Checkpoint 4: Verify summary and plotting data files
# -----------------------------------------------------------------------------

# Verify that the summary, quality report, cleaned data, pivot data, and NumPy
# arrays have all been created before the statistical tests are performed.
rule checkpoint_tables_complete:
    input:
        statistics_checkpoint="results/checkpoints/statistics_complete.ok",
        summary="results/tables/alignment_summary.tsv",
        quality="results/tables/quality_report.tsv",
        clean="results/tables/plot_data_clean.tsv",
        pivot="results/tables/plot_data_pivot.tsv",
        arrays="results/tables/plot_arrays.npz"
    output:
        checkpoint="results/checkpoints/tables_complete.ok"
    run:
        checkpoint_path = Path(output.checkpoint)
        checkpoint_path.parent.mkdir(parents=True, exist_ok=True)

        files_to_check = [
            input.summary,
            input.quality,
            input.clean,
            input.pivot,
            input.arrays,
        ]

        invalid_files = [
            str(file_path)
            for file_path in files_to_check
            if not Path(file_path).exists()
            or Path(file_path).stat().st_size == 0
        ]

        if invalid_files:
            raise RuntimeError(
                "Table checkpoint failed. Missing or empty files:\n"
                + "\n".join(invalid_files)
            )

        checkpoint_path.write_text(
            "TABLE CHECKPOINT PASSED\n"
            "Alignment summary, quality report, cleaned data, pivot data, "
            "and NumPy arrays are present and non-empty.\n"
        )


# -----------------------------------------------------------------------------
# Rule 9: Paired t-tests
# -----------------------------------------------------------------------------

# Perform paired statistical tests using the prepared pivot table
rule paired_t_tests:
    input:
        tables_checkpoint="results/checkpoints/tables_complete.ok",
        pivot="results/tables/plot_data_pivot.tsv"
    output:
        tests="results/tables/paired_alignment_tests.tsv"
    conda:
        ENV_ALIGNMENT
    shell:
        r"""
        python scripts/paired_t_tests.py \
            --input {input.pivot} \
            --out {output.tests}
        """


# -----------------------------------------------------------------------------
# Rule 10: Final alignment figures
# -----------------------------------------------------------------------------

# Generate the final benchmark figures
rule alignment_figures:
    input:
        clean="results/tables/plot_data_clean.tsv",
        pivot="results/tables/plot_data_pivot.tsv",
        arrays="results/tables/plot_arrays.npz",
        tests="results/tables/paired_alignment_tests.tsv"
    output:
        error_png="results/figures/final_alignment_error_rate.png",
        error_pdf="results/figures/final_alignment_error_rate.pdf",
        mapped_png="results/figures/mapped_bases.png",
        mapped_pdf="results/figures/mapped_bases.pdf",
        cigar_png="results/figures/mapped_bases_cigar.png",
        cigar_pdf="results/figures/mapped_bases_cigar.pdf"
    conda:
        ENV_ALIGNMENT
    shell:
        r"""
        mkdir -p results/figures

        python scripts/alignment_metrics.py \
            --cleaned-data {input.clean} \
            --pivoted-data {input.pivot} \
            --arrays {input.arrays} \
            --tests {input.tests} \
            --error-png {output.error_png} \
            --error-pdf {output.error_pdf} \
            --mapped-png {output.mapped_png} \
            --mapped-pdf {output.mapped_pdf} \
            --cigar-png {output.cigar_png} \
            --cigar-pdf {output.cigar_pdf}
        """

# -----------------------------------------------------------------------------
# Checkpoint 5: Confirm that the complete workflow finished successfully
# -----------------------------------------------------------------------------

# Verify that all stage checkpoints, statistical tests, and final figures are
# present before marking the complete workflow as successful.
rule checkpoint_workflow_complete:
    input:
        alignment_checkpoint="results/checkpoints/alignments_complete.ok",
        statistics_checkpoint="results/checkpoints/statistics_complete.ok",
        tables_checkpoint="results/checkpoints/tables_complete.ok",
        tests="results/tables/paired_alignment_tests.tsv",
        error_png="results/figures/final_alignment_error_rate.png",
        error_pdf="results/figures/final_alignment_error_rate.pdf",
        mapped_png="results/figures/mapped_bases.png",
        mapped_pdf="results/figures/mapped_bases.pdf",
        cigar_png="results/figures/mapped_bases_cigar.png",
        cigar_pdf="results/figures/mapped_bases_cigar.pdf"
    output:
        checkpoint="results/checkpoints/workflow_complete.ok"
    run:
        checkpoint_path = Path(output.checkpoint)
        checkpoint_path.parent.mkdir(parents=True, exist_ok=True)

        files_to_check = [
            input.tests,
            input.error_png,
            input.error_pdf,
            input.mapped_png,
            input.mapped_pdf,
            input.cigar_png,
            input.cigar_pdf,
        ]

        invalid_files = [
            str(file_path)
            for file_path in files_to_check
            if not Path(file_path).exists()
            or Path(file_path).stat().st_size == 0
        ]

        if invalid_files:
            raise RuntimeError(
                "Final workflow checkpoint failed. Missing or empty files:\n"
                + "\n".join(invalid_files)
            )

        checkpoint_path.write_text(
            "WORKFLOW CHECKPOINT PASSED\n"
            "All alignments, statistics, tables, statistical tests, and "
            "figures completed successfully.\n"
        )

