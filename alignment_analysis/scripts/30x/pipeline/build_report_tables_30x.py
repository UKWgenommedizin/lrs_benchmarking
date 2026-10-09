#!/usr/bin/env python3
"""Build every data table of the F2 thesis report from the 30x study tables.

Inputs (both 24 rows = 3 samples x 2 technologies x 4 aligners):
  tables/30x/source/alignment_summary_30x_Samtools_Christian.tsv
      samtools stats counts (reads, bases, CIGAR bases, error rate), written by
      build_alignment_summary_30x_table.py
  tables/30x/final/alignment_benchmark_30x.tsv
      canonical benchmark; supplies reads_mq0, threads, wall-clock runtime and peak RAM

Outputs, written next to the Christian table in tables/30x/source/:
  report_per_run_values_30x.tsv                            unrounded per-run values behind every table
  report_input_sequencing_data_30x.tsv                     tab:input_data
  report_thread_hours_30x.tsv                              tab:thread_hours
  report_alignment_summary_30x.tsv                         tab:alignment_summary
  report_cigar_yield_error_correlation_diagnostics_30x.tsv tab:supp_cigar_correlation_diagnostics

Definitions:
- Verified raw input = raw_total_sequences / total_length shared by the aligners that keep
  unmapped reads (reads_unmapped > 0). VACmap drops unmapped reads, so its own totals are
  smaller; every aligner is divided by the verified input instead.
- Coverage = input bases / 3.1 Gb; mean read length = input bases / input reads.
- Mapped reads % = reads_mapped / input reads; mapped bases % = bases_mapped / input bases;
  CIGAR yield % = bases_mapped_cigar / input bases; MQ0 % = reads_mq0 / reads_mapped.
- Thread-hours = allocated threads x wall-clock hours.
- Summary tables give mean +/- sample SD (ddof=1) over HG002, HG003 and HG004.
- Correlation diagnostics (per technology, n = 12): Pearson and Spearman between CIGAR yield %
  (x) and error_percent (y), linear regression y ~ x, Shapiro-Wilk on the regression residuals.
  error_percent is the Christian table's value rounded to 4 decimals, as the report used.

Usage: python3 alignment_analysis/scripts/30x/pipeline/build_report_tables_30x.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from scipy import stats

PROJECT = next(
    parent for parent in Path(__file__).resolve().parents
    if (parent / "alignment_analysis").is_dir()
)
TABLE_DIR = PROJECT / "alignment_analysis" / "tables" / "30x"
CHRISTIAN_TSV = TABLE_DIR / "source" / "alignment_summary_30x_Samtools_Christian.tsv"
BENCHMARK_TSV = TABLE_DIR / "final" / "alignment_benchmark_30x.tsv"
OUTPUT_DIR = CHRISTIAN_TSV.parent

GENOME_SIZE = 3.1e9
SAMPLE_ORDER = ["HG002", "HG003", "HG004"]
TECHNOLOGY_ORDER = ["ONT", "PacBio"]
TECHNOLOGY_TITLES = {"ONT": "ONT", "PacBio": "PacBio HiFi"}
ALIGNER_ORDER = ["minimap2", "pbmm2", "VACmap", "VG Giraffe"]
KEYS = ["sample", "read_technology", "aligner", "configuration"]

# alignment summary column -> decimals shown in the report
SUMMARY_METRICS = {
    "error_percent": 2,
    "mapped_reads_percent": 2,
    "mapped_bases_percent": 2,
    "cigar_yield_percent": 2,
    "mq0_reads_percent": 2,
    "runtime_hours": 2,
    "peak_ram_gb": 1,
}

#Read the parameters of Christian fiinal tsv files
def load_runs() -> pd.DataFrame:
    christian = pd.read_csv(CHRISTIAN_TSV, sep="\t")
    benchmark = pd.read_csv(BENCHMARK_TSV, sep="\t")
    benchmark["aligner"] = benchmark["aligner"].replace({"VACMap": "VACmap"})
    runs = christian[KEYS + [
        "raw_total_sequences", "reads_mapped", "reads_unmapped", "total_length",
        "bases_mapped", "bases_mapped_cigar", "error_percent",
    ]].merge(
        benchmark[KEYS + ["reads_mq0", "threads", "wallclock_runtime_hours", "peak_ram_gb",
                          "runtime_configuration"]],
        on=KEYS, how="inner", validate="one_to_one",
    )
    if len(runs) != 24:
        raise ValueError(f"Expected 24 runs after merging the Christian and benchmark tables, found {len(runs)}.")

    # verified raw input per sample x technology
    retaining = runs.loc[runs["reads_unmapped"] > 0]
    agreed = retaining.groupby(["sample", "read_technology"])[["raw_total_sequences", "total_length"]].nunique()
    if (agreed != 1).any().any():
        raise ValueError(f"Aligners retaining unmapped reads disagree on the input:\n{agreed}")
    verified = (retaining.groupby(["sample", "read_technology"])[["raw_total_sequences", "total_length"]].first()
                .rename(columns={"raw_total_sequences": "input_reads", "total_length": "input_bases"}))
    runs = runs.merge(verified, left_on=["sample", "read_technology"], right_index=True, how="left")

    runs["mapped_reads_percent"] = 100 * runs["reads_mapped"] / runs["input_reads"]
    runs["mapped_bases_percent"] = 100 * runs["bases_mapped"] / runs["input_bases"]
    runs["cigar_yield_percent"] = 100 * runs["bases_mapped_cigar"] / runs["input_bases"]
    runs["mq0_reads_percent"] = 100 * runs["reads_mq0"] / runs["reads_mapped"]
    runs["runtime_hours"] = runs["wallclock_runtime_hours"]
    runs["thread_hours"] = runs["threads"] * runs["runtime_hours"]

    order = {name: i for i, name in enumerate(TECHNOLOGY_ORDER + ALIGNER_ORDER + SAMPLE_ORDER)}
    return runs.sort_values(["read_technology", "aligner", "sample"], key=lambda s: s.map(order)).reset_index(drop=True)

#Check the values but in case are not present define n.m.
def mean_sd(values: pd.Series, decimals: int) -> str:
    values = values.dropna()
    if values.empty:
        return "n.m."
    return f"{values.mean():.{decimals}f} $\\pm$ {values.std(ddof=1):.{decimals}f}"


def groups(runs: pd.DataFrame):
    """Yield (technology title shown once per block, aligner, 3-sample group)."""
    for technology in TECHNOLOGY_ORDER:
        for index, aligner in enumerate(ALIGNER_ORDER):
            group = runs.loc[(runs["read_technology"] == technology) & (runs["aligner"] == aligner)]
            if len(group) != 3:
                raise ValueError(f"Expected 3 samples for {technology}/{aligner}, found {len(group)}.")
            if group["threads"].nunique() != 1:
                raise ValueError(f"Mixed thread counts for {technology}/{aligner}: {group['threads'].unique()}")
            yield TECHNOLOGY_TITLES[technology] if index == 0 else "", aligner, group


def input_data_table(runs: pd.DataFrame) -> pd.DataFrame:
    inputs = runs.drop_duplicates(["sample", "read_technology"])
    return pd.DataFrame({
        "sample": inputs["sample"],
        "technology": inputs["read_technology"].map(TECHNOLOGY_TITLES),
        "reads_m": (inputs["input_reads"] / 1e6).map("{:.2f}".format),
        "bases_gb": (inputs["input_bases"] / 1e9).map("{:.2f}".format),
        "coverage_x": (inputs["input_bases"] / GENOME_SIZE).map("{:.1f}".format),
        "mean_read_length_kb": (inputs["input_bases"] / inputs["input_reads"] / 1e3).map("{:.1f}".format),
    })


def thread_hours_table(runs: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame([{
        "technology": technology,
        "aligner": aligner,
        "threads": int(group["threads"].iloc[0]),
        "runtime_hours": mean_sd(group["runtime_hours"], 2),
        "thread_hours": mean_sd(group["thread_hours"], 1),
    } for technology, aligner, group in groups(runs)])


def alignment_summary_table(runs: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame([{
        "technology": technology or "",
        "aligner": aligner,
        "threads": int(group["threads"].iloc[0]),
        **{metric: mean_sd(group[metric], decimals) for metric, decimals in SUMMARY_METRICS.items()},
    } for technology, aligner, group in groups(runs)])


def correlation_table(runs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for technology in TECHNOLOGY_ORDER:
        group = runs.loc[runs["read_technology"] == technology]
        x, y = group["cigar_yield_percent"], group["error_percent"]
        pearson, spearman = stats.pearsonr(x, y), stats.spearmanr(x, y)
        regression = stats.linregress(x, y)
        shapiro = stats.shapiro(y - (regression.intercept + regression.slope * x))
        rows.append({
            "technology": technology,
            "n": len(group),
            "pearson_r": pearson.statistic,
            "pearson_p": pearson.pvalue,
            "spearman_rho": spearman.statistic,
            "spearman_p": spearman.pvalue,
            "shapiro_residual_W": shapiro.statistic,
            "shapiro_residual_p": shapiro.pvalue,
            "r_squared": regression.rvalue ** 2,
            "regression_slope": regression.slope,
            "pearson_spearman_difference": abs(pearson.statistic - spearman.statistic),
        })
    return pd.DataFrame(rows)


def main() -> int:
    runs = load_runs()
    outputs = {
        "report_per_run_values_30x.tsv": runs,
        "report_input_sequencing_data_30x.tsv": input_data_table(runs),
        "report_thread_hours_30x.tsv": thread_hours_table(runs),
        "report_alignment_summary_30x.tsv": alignment_summary_table(runs),
        "report_cigar_yield_error_correlation_diagnostics_30x.tsv": correlation_table(runs),
    }
    for name, table in outputs.items():
        table.to_csv(OUTPUT_DIR / name, sep="\t", index=False)
        if name != "report_per_run_values_30x.tsv":
            print(f"\n== {name}\n{table.to_string(index=False)}")
    print(f"\nWrote {len(outputs)} tables to {OUTPUT_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
