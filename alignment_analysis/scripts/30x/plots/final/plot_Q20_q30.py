#!/usr/bin/env python3
"""Plot the percentage of Q20 and Q30 bases per GIAB sample and technology.

Two bars per sample/technology: the share of all 30x input bases whose
Phred base quality is >= Q20 and >= Q30.

Metric choice
-------------
Per-base, not per-read: the values come from the per-cycle base-quality
(FFQ) histogram of the full samtools stats of each 30x CRAM that retains
every input read, which covers every input base but does not keep reads
together. A per-read mean-quality Q20/Q30 would need the 30x FASTQs,
which are only on the server.

Input
-----
`tables/30x/derived/read_qc_30x.tsv`, built by
`scripts/30x/pipeline/build_read_qc_30x.py` (which checks that the
histogram base total equals the verified raw input bases and that two
independent CRAMs of the same input give identical values). This script
re-checks completeness against the raw input tables and that Q30 <= Q20.
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
    rotated_xticks,
    save_figure,
)

TABLE_DIR = PROJECT / "alignment_analysis" / "tables"
FIGURE_DIR = PROJECT / "alignment_analysis" / "figures" / "30x" / "final"

INPUT_TSV = TABLE_DIR / "30x" / "derived" / "read_qc_30x.tsv"
RAW_BASES_TSV = TABLE_DIR / "30x" / "derived" / "plot_data" / "raw_base_counts_30x.tsv"
OUTPUT_SUMMARY = TABLE_DIR / "30x" / "derived" / "plot_data" / "q20_q30_bases_30x.tsv"
OUTPUT_PNG = FIGURE_DIR / "16_q20_q30_reads_30x.png"
OUTPUT_PDF = OUTPUT_PNG.with_suffix(".pdf")

QUALITY_LEVELS = ["Q20", "Q30"]

BAR_WIDTH = 0.38

# Smaller canvas + larger fonts so the text stays legible when the
# panel is scaled to Figure 1 (e/f) in the report.
PANEL_WIDTH_IN = 3.0
PANEL_HEIGHT_IN = 2.7


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
        raise FileNotFoundError(
            f"Input TSV not found:\n{INPUT_TSV}\nRun scripts/30x/pipeline/build_read_qc_30x.py first."
        )

    data = pd.read_csv(INPUT_TSV, sep="\t")
    plot_data = data.loc[
        data["sample"].isin(SAMPLE_ORDER) & data["read_technology"].isin(TECHNOLOGY_ORDER)
    ].copy()

    duplicated = plot_data.duplicated(subset=["sample", "read_technology"], keep=False)
    if duplicated.any():
        raise ValueError("Duplicated sample/technology rows found:\n" + plot_data.loc[duplicated].to_string(index=False))

    expected_rows = len(SAMPLE_ORDER) * len(TECHNOLOGY_ORDER)
    if len(plot_data) != expected_rows:
        raise ValueError(f"Expected {expected_rows} sample x technology rows, found {len(plot_data)}.")

    raw_bases = pd.read_csv(RAW_BASES_TSV, sep="\t")
    plot_data = plot_data.merge(raw_bases[["sample", "read_technology", "raw_input_bases"]],
                                on=["sample", "read_technology"], how="left")
    if (plot_data["total_bases"] != plot_data["raw_input_bases"]).any():
        raise ValueError("read_qc_30x.tsv does not cover all 30x input bases.")
    if (plot_data["q30_bases_percent"] > plot_data["q20_bases_percent"]).any():
        raise ValueError("Q30 bases exceed Q20 bases for at least one sample/technology.")

    plot_data = plot_data.rename(columns={"q20_bases_percent": "Q20_percent", "q30_bases_percent": "Q30_percent"})

    print()
    print("Q20/Q30 bases from all 30x input bases (per sample/technology):")
    for _, row in plot_data.iterrows():
        print(
            f"  {row['sample']} {row['read_technology']:6s} bases={row['total_bases']:>12d}  "
            f"Q20={row['Q20_percent']:6.2f}%  Q30={row['Q30_percent']:6.2f}%  ({row['stats_file']})"
        )

    plot_data["sample"] = pd.Categorical(plot_data["sample"], SAMPLE_ORDER, ordered=True)
    plot_data["read_technology"] = pd.Categorical(plot_data["read_technology"], TECHNOLOGY_ORDER, ordered=True)

    return plot_data.sort_values(["read_technology", "sample"])[
        ["sample", "read_technology", "total_bases", "Q20_percent", "Q30_percent", "mean_base_quality", "stats_file"]
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

    # All values lie between ~86 and 100 %, so the axis is truncated
    # (as in panels b and c) to keep the differences visible.
    y_min_display = 80.0
    y_max_display = 107.0
    value_font_size = 8.0

    figure.patch.set_facecolor("white")
    figure.subplots_adjust(left=0.19, right=0.98, bottom=0.22, top=0.80, wspace=0.12)

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

        axis.set_title(TECHNOLOGY_TITLES[technology], fontsize=12, pad=5)
        rotated_xticks(axis, x_positions, SAMPLE_ORDER, fontsize=10)
        axis.set_xlim(-0.6, len(SAMPLE_ORDER) - 0.4)
        axis.set_ylim(y_min_display, y_max_display)
        axis.set_yticks([80, 85, 90, 95, 100])
        axis.tick_params(axis="y", labelsize=10)
        for spine_name, spine in axis.spines.items():
            spine.set_visible(spine_name in ("left", "bottom"))
        axis.spines["left"].set_linewidth(0.7)
        axis.spines["bottom"].set_linewidth(0.7)
        axis.set_facecolor("white")

    axes[0].set_ylabel("Bases (%)", fontsize=11)

    # Neutral legend: saturated = Q20, light = Q30 (colour encodes sample).
    legend_handles = [
        Patch(facecolor="#555555", edgecolor="white", label="≥ Q20"),
        Patch(facecolor=lighten_color("#555555"), edgecolor="white", label="≥ Q30"),
    ]
    figure.legend(
        handles=legend_handles, loc="upper center", bbox_to_anchor=(0.575, 0.985),
        ncol=2, frameon=False, fontsize=10, handlelength=1.2, columnspacing=1.2,
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
