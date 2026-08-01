#!/usr/bin/env bash
set -euo pipefail

echo "============================================================"
echo " UPDATE WORKFLOW TO USE samples.tsv"
echo " DRY-RUN ONLY: HG002 + HG003 + HG004"
echo " Aligners: minimap2 + pbmm2 + VACmap + vg Giraffe"
echo "============================================================"

WORKDIR="$HOME/lrs_benchmarking/final_report_files/snakemake_aligners_benchmarking"
PROJECT_ROOT="$HOME/lrs_benchmarking"
FASTQ_SOURCE="$PROJECT_ROOT/samples_try"
REFERENCE_SOURCE="$PROJECT_ROOT/reference/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta"

cd "$WORKDIR"

echo "[1/8] Creating folders..."
mkdir -p data/samples data/reference data/graph .dryrun_inputs/graph envs results/logs results/bam results/stats results/tables results/figures

echo "[2/8] Backing up current files..."
STAMP="$(date +%Y%m%d_%H%M%S)"
for f in snakemake_aligners.smk config.yaml config.dryrun.yaml samples.tsv; do
    if [ -f "$f" ]; then
        cp "$f" "${f}.backup_${STAMP}"
        echo "Backup created: ${f}.backup_${STAMP}"
    fi
done

echo "[3/8] Linking real FASTQ files into data/samples..."
MISSING=0

for sample in HG002 HG003 HG004; do
    for tech in ont pb; do
        SRC="$FASTQ_SOURCE/${sample}.${tech}.1k.fastq.gz"
        DEST="data/samples/${sample}.${tech}.1k.fastq.gz"

        if [ -f "$SRC" ]; then
            ln -sf "$SRC" "$DEST"
            echo "Linked: $DEST -> $SRC"
        else
            echo "MISSING FASTQ: $SRC"
            MISSING=1
        fi
    done
done

echo "[4/8] Linking reference genome..."
if [ -f "$REFERENCE_SOURCE" ]; then
    ln -sf "$REFERENCE_SOURCE" data/reference/genome.fasta
    echo "Linked: data/reference/genome.fasta -> $REFERENCE_SOURCE"
else
    echo "MISSING REFERENCE: $REFERENCE_SOURCE"
    MISSING=1
fi

if [ "$MISSING" -eq 1 ]; then
    echo
    echo "STOP: Some FASTQ/reference files are missing."
    echo "The workflow files were not updated because paths must be fixed first."
    exit 1
fi

echo "[5/8] Writing samples.tsv..."
cat > samples.tsv <<'TSV'
sample	technology	fastq
HG002	ont	data/samples/HG002.ont.1k.fastq.gz
HG002	pb	data/samples/HG002.pb.1k.fastq.gz
HG003	ont	data/samples/HG003.ont.1k.fastq.gz
HG003	pb	data/samples/HG003.pb.1k.fastq.gz
HG004	ont	data/samples/HG004.ont.1k.fastq.gz
HG004	pb	data/samples/HG004.pb.1k.fastq.gz
TSV

echo
echo "samples.tsv:"
column -t -s $'\t' samples.tsv || cat samples.tsv
echo

echo "[6/8] Writing config files..."

cat > config.yaml <<'YAML'
sample_sheet: samples.tsv
reference: data/reference/genome.fasta

active_aligners:
  - minimap2
  - pbmm2
  - vacmap
  - vg_giraffe

vacmap_bin: "conda run -n vacmap_env vacmap"

vg_gbz: data/graph/graph.gbz
vg_min: data/graph/graph.min
vg_dist: data/graph/graph.dist
vg_zipcodes: ""
YAML

# Dry-run config uses dummy vg graph placeholders so the DAG can be tested
# without real vg graph indexes.
touch .dryrun_inputs/graph/graph.gbz
touch .dryrun_inputs/graph/graph.min
touch .dryrun_inputs/graph/graph.dist

cat > config.dryrun.yaml <<'YAML'
sample_sheet: samples.tsv
reference: data/reference/genome.fasta

active_aligners:
  - minimap2
  - pbmm2
  - vacmap
  - vg_giraffe

vacmap_bin: "conda run -n vacmap_env vacmap"

vg_gbz: .dryrun_inputs/graph/graph.gbz
vg_min: .dryrun_inputs/graph/graph.min
vg_dist: .dryrun_inputs/graph/graph.dist
vg_zipcodes: ""
YAML

echo "[7/8] Writing environment YAML files..."

cat > envs/alignment.yaml <<'YAML'
channels:
  - bioconda
  - conda-forge
  - nodefaults

dependencies:
  - python=3.11
  - minimap2
  - pbmm2
  - samtools
  - pandas
  - numpy
  - scipy
  - matplotlib
  - pip
YAML

cat > envs/vg.yaml <<'YAML'
channels:
  - bioconda
  - conda-forge
  - nodefaults

dependencies:
  - vg
  - samtools
YAML

echo "[8/8] Writing updated snakemake_aligners.smk..."

cat > snakemake_aligners.smk <<'SNAKE'
#####################################################################################
# ALIGNMENT PIPELINE: MINIMAP2, PBMM2, VACMAP, AND VG GIRAFFE
# Sample information is read from samples.tsv through config.yaml.
#####################################################################################

from pathlib import Path
import csv


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

SAMPLE_SHEET = config.get("sample_sheet", "samples.tsv")
REFERENCE = config.get("reference", "data/reference/genome.fasta")

ACTIVE_ALIGNERS = config.get(
    "active_aligners",
    ["minimap2", "pbmm2", "vacmap", "vg_giraffe"]
)

VACMAP_BIN = config.get("vacmap_bin", "conda run -n vacmap_env vacmap")

VG_GBZ = config.get("vg_gbz", "data/graph/graph.gbz")
VG_MIN = config.get("vg_min", "data/graph/graph.min")
VG_DIST = config.get("vg_dist", "data/graph/graph.dist")
VG_ZIPCODES = config.get("vg_zipcodes", "")

PREFLIGHT_OK = "results/logs/preflight.ok"

ENV_ALIGNMENT = "envs/alignment.yaml"
ENV_VG = "envs/vg.yaml"


# -----------------------------------------------------------------------------
# Read sample sheet
# -----------------------------------------------------------------------------

def read_sample_sheet(path):
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(f"Missing sample sheet: {path}")

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


SAMPLE_ROWS = read_sample_sheet(SAMPLE_SHEET)

SAMPLES = list(dict.fromkeys(row["sample"] for row in SAMPLE_ROWS))
TECHNOLOGIES = list(dict.fromkeys(row["technology"] for row in SAMPLE_ROWS))

SAMPLE_TECH_PAIRS = list(
    dict.fromkeys((row["sample"], row["technology"]) for row in SAMPLE_ROWS)
)

FASTQ_BY_SAMPLE_TECH = {
    (row["sample"], row["technology"]): row["fastq"]
    for row in SAMPLE_ROWS
}


def fastq_for(wildcards):
    key = (wildcards.sample, wildcards.technology)

    if key not in FASTQ_BY_SAMPLE_TECH:
        raise ValueError(
            f"No FASTQ defined for sample={wildcards.sample}, "
            f"technology={wildcards.technology}"
        )

    return FASTQ_BY_SAMPLE_TECH[key]


EXPECTED_STATS = [
    f"results/stats/{sample}.{technology}.{aligner}.stats.txt"
    for sample, technology in SAMPLE_TECH_PAIRS
    for aligner in ACTIVE_ALIGNERS
]


wildcard_constraints:
    sample="|".join(SAMPLES),
    technology="|".join(TECHNOLOGIES),
    aligner="|".join(ACTIVE_ALIGNERS)


# -----------------------------------------------------------------------------
# Preset / mode functions
# -----------------------------------------------------------------------------

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
# Final target
# -----------------------------------------------------------------------------

rule all:
    input:
        PREFLIGHT_OK,
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
        "results/figures/mapped_bases_cigar.pdf"


# -----------------------------------------------------------------------------
# Preflight check
# -----------------------------------------------------------------------------

rule preflight:
    output:
        PREFLIGHT_OK
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
                "\nPREFLIGHT CHECK FAILED\n"
                "Fix these problems before running the workflow:\n\n"
                f"{message}\n"
            )

        Path(output[0]).write_text("Preflight check passed.\n")


# -----------------------------------------------------------------------------
# Rule 1: minimap2 alignment
# -----------------------------------------------------------------------------

rule align_minimap2:
    input:
        preflight=PREFLIGHT_OK,
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
        preflight=PREFLIGHT_OK,
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
        preflight=PREFLIGHT_OK,
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

rule align_vg_giraffe:
    input:
        preflight=PREFLIGHT_OK,
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
# Rule 5: Samtools quality statistics
# -----------------------------------------------------------------------------

rule samtools_quality:
    input:
        preflight=PREFLIGHT_OK,
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
# Rule 6: Build alignment summary table
# -----------------------------------------------------------------------------

rule alignment_summary:
    input:
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

rule quality_check:
    input:
        summary="results/tables/alignment_summary.tsv"
    output:
        report="results/tables/quality_report.tsv"
    conda:
        ENV_ALIGNMENT
    shell:
        r"""
        python scripts/quality_check.py \
            --input {input.summary} \
            --out {output.report}
        """


# -----------------------------------------------------------------------------
# Rule 8: Prepare data for plotting
# -----------------------------------------------------------------------------

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
# Rule 9: Paired t-tests
# -----------------------------------------------------------------------------

rule paired_t_tests:
    input:
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
SNAKE

echo
echo "============================================================"
echo " FILES UPDATED"
echo "============================================================"
echo "Updated:"
echo "  samples.tsv"
echo "  config.yaml"
echo "  config.dryrun.yaml"
echo "  envs/alignment.yaml"
echo "  envs/vg.yaml"
echo "  snakemake_aligners.smk"
echo
echo "Now running dry-run only..."
echo

snakemake -s snakemake_aligners.smk --configfile config.dryrun.yaml --sdm conda --cores 4 -np --forceall | tee dryrun_sample_sheet.log

echo
echo "============================================================"
echo " DRY-RUN COMPLETE"
echo "============================================================"
echo "Expected alignment jobs:"
echo "  3 samples × 2 technologies × 4 aligners = 24 planned alignment jobs"
echo
echo "Dry-run log saved at:"
echo "  dryrun_sample_sheet.log"
echo
echo "No real alignment was executed."
echo "No BAM files were created."
