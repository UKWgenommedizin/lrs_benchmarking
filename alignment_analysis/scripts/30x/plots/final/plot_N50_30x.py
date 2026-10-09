#!/usr/bin/env python3
"""Plot the percentage of Q20 and Q30 reads per GIAB sample and technology.

Two bars per sample/technology: the share of reads whose mean base
quality is >= Q20 and >= Q30.

Metric choice
-------------
A read counts as Q20 (Q30) when the arithmetic mean of its Phred base
qualities is >= 20 (>= 30), as computed by `scripts/fastq_quality.awk`.
Q20/Q30 percentages = Q20/Q30 reads / total reads.

Input subset
------------
samtools stats of the 30x alignments does not report per-read base
quality, and the whole-genome 30x FASTQs are not stored locally. The
values therefore come from `tables/fastq_summary.tsv`, which was
computed on the 1,000-read test subsets (`HG00x.{ont,pb}.1k.fastq.gz`)
drawn from the same GIAB datasets. The script checks that every
sample/technology is present exactly once and that the stored
percentages match the stored read counts.
"""

from __future__ import annotations

from pathlib import Path
import sys

PROJECT = next(
    parent for parent in Path(__file__).resolve().parents
    if (parent / "alignment_analysis").is_dir()
)
sys.path.insert(0, str(PROJECT / "alignment_analysis" / "scripts"))

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib.patches import Patch

from utils.plot_style import (
    SAMPLE_COLORS,
    SAMPLE_ORDER,
    TECHNOLOGY_ORDER,
    TECHNOLOGY_TITLES,
    apply_style,
    save_figure,
)

TABLE_DIR = PROJECT / "alignment_analysis" / "tables"
FIGURE_DIR = PROJECT / "alignment_analysis" / "figures" / "30x" / "final"

INPUT_TSV = TABLE_DIR / "fastq_summary.tsv"
OUTPUT_SUMMARY = TABLE_DIR / "30x" / "derived" / "plot_data" / "q20_q30_reads_30x.tsv"
OUTPUT_PNG = FIGURE_DIR / "16_q20_q30_reads_30x.png"
OUTPUT_PDF = OUTPUT_PNG.with_suffix(".pdf")

TECHNOLOGY_CODES = {"ont": "ONT", "pb": "PacBio"}
QUALITY_LEVELS = ["Q20", "Q30"]

BAR_WIDTH = 0.38

# Sized to match the other report Figure 1 input panels.
PANEL_WIDTH_IN = 3.5
PANEL_HEIGHT_IN = 3.1


def darken_color(color, amount=0.35):
    rgb = np.array(mcolors.to_rgb(color))
    black = np.array([0.0, 0.0, 0.0])
    return tuple(rgb + (black - rgb) * amount)


def lighten_color(color, amount=0.55):
    rgb = np.array(mcolors.to_rgb(color))
    white = np.array([1.0, 1.0, 1.0])
    return tuple(rgb + (white - rgb) * amount)


def load_and_verify_quality() -> pd.DataFrame:
    if not INPUT_TSV.exists():
        raise FileNotFoundError(f"Input TSV not found:\n{INPUT_TSV}")

    data = pd.read_csv(INPUT_TSV, sep="\t")

    required_columns = {
        "file", "sample", "technology", "read_count",
        "Q20_reads", "Q20_percent", "Q30_reads", "Q30_percent",
    }
    missing_columns = required_columns.difference(data.columns)
    if missing_columns:
        raise ValueError(f"Missing required columns: {sorted(missing_columns)}")

    data["read_technology"] = data["technology"].map(TECHNOLOGY_CODES)
    plot_data = data.loc[
        data["sample"].isin(SAMPLE_ORDER) & data["read_technology"].isin(TECHNOLOGY_ORDER)
    ].copy()

    for column in ["read_count", "Q20_reads", "Q20_percent", "Q30_reads", "Q30_percent"]:
        plot_data[column] = pd.to_numeric(plot_data[column], errors="raise")

    duplicated = plot_data.duplicated(subset=["sample", "read_technology"], keep=False)
    if duplicated.any():
        raise ValueError("Duplicated sample/technology rows found:\n" + plot_data.loc[duplicated].to_string(index=False))

    expected_rows = len(SAMPLE_ORDER) * len(TECHNOLOGY_ORDER)
    if len(plot_data) != expected_rows:
        raise ValueError(f"Expected {expected_rows} sample x technology rows, found {len(plot_data)}.")

    # -----------------------------------------------------------
    # Verify the stored percentages against the stored read counts
    # and that Q30 reads are a subset of Q20 reads.
    # -----------------------------------------------------------

    print()
    print("Q20/Q30 verification (per sample/technology):")

    for level in QUALITY_LEVELS:
        recomputed = 100 * plot_data[f"{level}_reads"] / plot_data["read_count"]
        mismatch = (recomputed - plot_data[f"{level}_percent"]).abs() > 0.01
        if mismatch.any():
            raise ValueError(
                f"{level}_percent does not match {level}_reads / read_count:\n"
                + plot_data.loc[mismatch].to_string(index=False)
            )
        plot_data[f"{level}_percent"] = recomputed

    if (plot_data["Q30_reads"] > plot_data["Q20_reads"]).any():
        raise ValueError("Q30 reads exceed Q20 reads for at least one sample/technology.")

    for _, row in plot_data.iterrows():
        print(
            f"  {row['sample']} {row['read_technology']:6s} {row['file']:24s} "
            f"reads={row['read_count']:>6d}  Q20={row['Q20_percent']:6.2f}%  Q30={row['Q30_percent']:6.2f}%"
        )

    print()
    print("Q20/Q30 values are consistent for every sample/technology group.")

    plot_data["sample"] = pd.Categorical(plot_data["sample"], SAMPLE_ORDER, ordered=True)
    plot_data["read_technology"] = pd.Categorical(plot_data["read_technology"], TECHNOLOGY_ORDER, ordered=True)

    return plot_data.sort_values(["read_technology", "sample"])[
        ["sample", "read_technology", "file", "read_count", "Q20_reads", "Q20_percent", "Q30_reads", "Q30_percent"]
    ]


def main() -> int:
    plot_data = load_and_verify_quality()

    OUTPUT_SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    plot_data.to_csv(OUTPUT_SUMMARY, sep="\t", index=False)

    print()
    print("Plot data:")
    print(plot_data.to_string(index=False, float_format=lambda value: f"{value:.2f}"))

    apply_style()

    x_positions = np.arange(len(SAMPLE_ORDER))
    offsets = {"Q20": -BAR_WIDTH / 2, "Q30": BAR_WIDTH / 2}
    figure, axes = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(PANEL_WIDTH_IN, PANEL_HEIGHT_IN))

    # All values lie between ~90 and 100 %, so the axis is truncated
    # (as in panels b and c) to keep the differences visible.
    y_min_display = 80.0
    y_max_display = 104.0
    value_font_size = 6.0

    figure.patch.set_facecolor("white")
    figure.subplots_adjust(left=0.17, right=0.98, bottom=0.12, top=0.84, wspace=0.12)

    for panel_index, technology in enumerate(TECHNOLOGY_ORDER):
        axis = axes[panel_index]
        technology_data = (
            plot_data[plot_data["read_technology"] == technology].set_index("sample").reindex(SAMPLE_ORDER)
        )

        for level in QUALITY_LEVELS:
            values = technology_data[f"{level}_percent"].to_numpy(dtype=float)
            for position, sample, value in zip(x_positions, SAMPLE_ORDER, values):
                sample_color = SAMPLE_COLORS[sample]
                axis.bar(
                    position + offsets[level], value, width=BAR_WIDTH,
                    color=sample_color if level == "Q20" else lighten_color(sample_color),
                    edgecolor="white", linewidth=0.6, zorder=3,
                )
                axis.text(
                    position + offsets[level], value + (y_max_display - y_min_display) * 0.012, f"{value:.1f}",
                    ha="center", va="bottom", rotation=90, fontsize=value_font_size,
                    color=darken_color(sample_color, amount=0.15), zorder=7,
                )

        axis.set_title(TECHNOLOGY_TITLES[technology], fontsize=10, pad=5)
        axis.set_xticks(x_positions, SAMPLE_ORDER, fontsize=8)
        axis.set_xlim(-0.6, len(SAMPLE_ORDER) - 0.4)
        axis.set_ylim(y_min_display, y_max_display)
        axis.set_yticks([80, 85, 90, 95, 100])
        axis.tick_params(axis="y", labelsize=8)
        for spine_name, spine in axis.spines.items():
            spine.set_visible(spine_name in ("left", "bottom"))
        axis.spines["left"].set_linewidth(0.7)
        axis.spines["bottom"].set_linewidth(0.7)
        axis.set_facecolor("white")

    axes[0].set_ylabel("Reads (%)", fontsize=9)

    # Neutral legend: saturated = Q20, light = Q30 (colour encodes sample).
    legend_handles = [
        Patch(facecolor="#555555", edgecolor="white", label="≥ Q20"),
        Patch(facecolor=lighten_color("#555555"), edgecolor="white", label="≥ Q30"),
    ]
    figure.legend(
        handles=legend_handles, loc="upper center", bbox_to_anchor=(0.575, 0.985),
        ncol=2, frameon=False, fontsize=8, handlelength=1.2, columnspacing=1.2,
    )

    save_figure(figure, OUTPUT_PNG, OUTPUT_PDF)

    print()
    print(f"Input: {INPUT_TSV}")
    print(f"Plot data: {OUTPUT_SUMMARY}")
    print(f"PNG: {OUTPUT_PNG}")
    print(f"PDF: {OUTPUT_PDF}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
