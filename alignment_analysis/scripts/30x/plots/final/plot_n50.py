#!/usr/bin/env python3
"""Plot the read-length N50 (kb) per GIAB sample and technology.

One bar per sample/technology.

Metric choice
-------------
N50 = the length L such that reads of length >= L contain at least half
of all sequenced bases. Unlike the mean, it weights reads by the bases
they contribute, so it reflects where most of the sequence sits.

Input
-----
All 30x input reads: `tables/30x/derived/read_qc_30x.tsv`, built by
`scripts/30x/pipeline/build_read_qc_30x.py` from the read-length (RL)
histogram of the full samtools stats of each 30x CRAM that retains every
input read. The script re-checks that read count and total bases equal the
verified raw input of each dataset before plotting.
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
RAW_READS_TSV = TABLE_DIR / "30x" / "derived" / "plot_data" / "raw_read_counts_30x.tsv"
RAW_BASES_TSV = TABLE_DIR / "30x" / "derived" / "plot_data" / "raw_base_counts_30x.tsv"
OUTPUT_SUMMARY = TABLE_DIR / "30x" / "derived" / "plot_data" / "read_n50_30x.tsv"
OUTPUT_PNG = FIGURE_DIR / "17_read_n50_30x.png"
OUTPUT_PDF = OUTPUT_PNG.with_suffix(".pdf")

BAR_WIDTH = 1.0

# Smaller canvas + larger fonts so the text stays legible when the
# panel is scaled to Figure 1 (e/f) in the report.
PANEL_WIDTH_IN = 3.0
PANEL_HEIGHT_IN = 2.7


def darken_color(color, amount=0.35):
    rgb = np.array(mcolors.to_rgb(color))
    black = np.array([0.0, 0.0, 0.0])
    return tuple(rgb + (black - rgb) * amount)


def load_and_verify_n50() -> pd.DataFrame:
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

    raw = pd.read_csv(RAW_READS_TSV, sep="\t").merge(
        pd.read_csv(RAW_BASES_TSV, sep="\t"), on=["sample", "read_technology"]
    )
    plot_data = plot_data.merge(raw[["sample", "read_technology", "raw_input_reads", "raw_input_bases"]],
                                on=["sample", "read_technology"], how="left")
    mismatch = (plot_data["reads"] != plot_data["raw_input_reads"]) | (plot_data["total_bases"] != plot_data["raw_input_bases"])
    if mismatch.any():
        raise ValueError("read_qc_30x.tsv does not cover the full 30x input:\n" + plot_data.loc[mismatch].to_string(index=False))
    if ((plot_data["read_n50"] <= 0) | (plot_data["read_n50"] > plot_data["max_read_length"])).any():
        raise ValueError("Read N50 outside (0, max read length] for at least one row.")

    print()
    print("Read N50 from all 30x input reads (per sample/technology):")
    for _, row in plot_data.iterrows():
        print(
            f"  {row['sample']} {row['read_technology']:6s} reads={row['reads']:>8d}  "
            f"bases={row['total_bases']:>12d}  N50={row['read_n50']:>7d}  ({row['stats_file']})"
        )

    plot_data = plot_data.rename(columns={"reads": "read_count", "read_n50": "N50", "read_n50_kb": "N50_kb"})
    plot_data["sample"] = pd.Categorical(plot_data["sample"], SAMPLE_ORDER, ordered=True)
    plot_data["read_technology"] = pd.Categorical(plot_data["read_technology"], TECHNOLOGY_ORDER, ordered=True)

    return plot_data.sort_values(["read_technology", "sample"])[
        ["sample", "read_technology", "read_count", "total_bases", "N50", "N50_kb", "stats_file"]
    ]


def main() -> int:
    plot_data = load_and_verify_n50()

    OUTPUT_SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    plot_data.to_csv(OUTPUT_SUMMARY, sep="\t", index=False)

    print()
    print("Plot data:")
    print(plot_data.to_string(index=False, float_format=lambda value: f"{value:.3f}"))

    apply_style()

    x_positions = np.arange(len(SAMPLE_ORDER))
    figure, axes = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(PANEL_WIDTH_IN, PANEL_HEIGHT_IN))

    # Bars encode length, so the axis starts at zero.
    y_max_display = 40.0
    value_font_size = 9.0

    figure.patch.set_facecolor("white")
    figure.subplots_adjust(left=0.19, right=0.98, bottom=0.22, top=0.88, wspace=0.12)

    for panel_index, technology in enumerate(TECHNOLOGY_ORDER):
        axis = axes[panel_index]
        technology_data = (
            plot_data[plot_data["read_technology"] == technology].set_index("sample").reindex(SAMPLE_ORDER)
        )
        values = technology_data["N50_kb"].to_numpy(dtype=float)

        for position, sample, value in zip(x_positions, SAMPLE_ORDER, values):
            axis.bar(
                position, value, width=BAR_WIDTH,
                color=SAMPLE_COLORS[sample], edgecolor="white", linewidth=0.6, zorder=3,
            )
            axis.text(
                position, value + y_max_display * 0.012, f"{value:.1f}",
                ha="center", va="bottom", fontsize=value_font_size,
                color=darken_color(SAMPLE_COLORS[sample], amount=0.15), zorder=7,
            )

        axis.set_title(TECHNOLOGY_TITLES[technology], fontsize=12, pad=5)
        rotated_xticks(axis, x_positions, SAMPLE_ORDER, fontsize=10)
        axis.set_xlim(-0.5, len(SAMPLE_ORDER) - 0.5)
        axis.set_ylim(0.0, y_max_display)
        axis.tick_params(axis="y", labelsize=10)
        for spine_name, spine in axis.spines.items():
            spine.set_visible(spine_name in ("left", "bottom"))
        axis.spines["left"].set_linewidth(0.7)
        axis.spines["bottom"].set_linewidth(0.7)
        axis.set_facecolor("white")

    axes[0].set_ylabel("Read N50 (kb)", fontsize=11)

    save_figure(figure, OUTPUT_PNG, OUTPUT_PDF)

    print()
    print(f"Input: {INPUT_TSV}")
    print(f"Plot data: {OUTPUT_SUMMARY}")
    print(f"PNG: {OUTPUT_PNG}")
    print(f"PDF: {OUTPUT_PDF}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
