#!/usr/bin/env python3
"""Plot CIGAR-aligned, mismatch and CIGAR-unaligned bases (%) as a stacked bar.

Three-segment stack (aligned non-mismatch / mismatch / CIGAR-unaligned)
per sample per aligner, from a truncated 88-100% baseline.

Every segment is a share of the verified raw input bases of that
sample/technology (the total_length that all aligners retaining unmapped
reads agree on), the same denominator as the mapped/unmapped panels and
the input-normalized CIGAR yield. VACmap drops unmapped reads from its
output, so dividing by its own total_length would overstate its yield;
the bases of those dropped reads are counted as CIGAR-unaligned here.
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
    ALIGNER_ORDER,
    HALF_WIDTH_IN,
    PANEL_LABEL_SIZE_PT,
    PANEL_TICK_SIZE,
    PANEL_TITLE_SIZE,
    PANEL_VALUE_SIZE,
    SAMPLE_COLORS,
    SAMPLE_ORDER,
    TECHNOLOGY_ORDER,
    TECHNOLOGY_TITLES,
    rotated_xticks,
    apply_style,
    clean_spines,
    panel_legend,
    sample_legend_handles,
    save_figure,
)

TABLE_DIR = PROJECT / "alignment_analysis" / "tables" / "30x"
FIGURE_DIR = PROJECT / "alignment_analysis" / "figures" / "30x" / "final"

INPUT_TSV = TABLE_DIR / "final" / "alignment_benchmark_30x.tsv"
OUTPUT_SUMMARY = TABLE_DIR / "derived" / "plot_data" / "cigar_mapped_mismatch_unaligned_percent.tsv"
OUTPUT_PNG = FIGURE_DIR / "06_cigar_composition_30x.png"
OUTPUT_PDF = OUTPUT_PNG.with_suffix(".pdf")

SAMPLE_OFFSETS = {"HG002": -0.27, "HG003": 0.00, "HG004": 0.27}
BAR_WIDTH = 0.26


def lighten_color(color, amount=0.5):
    rgb = np.array(mcolors.to_rgb(color))
    white = np.array([1.0, 1.0, 1.0])
    return tuple(rgb + (white - rgb) * amount)


def load_data() -> pd.DataFrame:
    if not INPUT_TSV.exists():
        raise FileNotFoundError(f"Input TSV not found:\n{INPUT_TSV}")

    data = pd.read_csv(INPUT_TSV, sep="\t")
    data["aligner"] = data["aligner"].replace({"VACMap": "VACmap"})

    required_columns = {
        "sample", "read_technology", "aligner", "total_length", "bases_mapped_cigar", "mismatches", "reads_unmapped",
    }
    missing_columns = required_columns.difference(data.columns)
    if missing_columns:
        raise ValueError(f"Missing required columns: {sorted(missing_columns)}")

    for column in ["total_length", "bases_mapped_cigar", "mismatches", "reads_unmapped"]:
        data[column] = pd.to_numeric(data[column], errors="raise")

    if (data["total_length"] <= 0).any():
        raise ValueError("At least one total_length value is <= 0.")

    data = data.loc[
        data["sample"].isin(SAMPLE_ORDER)
        & data["read_technology"].isin(TECHNOLOGY_ORDER)
        & data["aligner"].isin(ALIGNER_ORDER)
    ].copy()

    # Verified raw input per sample/technology: the total_length that all
    # aligners retaining unmapped reads (reads_unmapped > 0) agree on.
    retaining = data.loc[data["reads_unmapped"] > 0]
    raw_input = retaining.groupby(["sample", "read_technology"])["total_length"].agg(["nunique", "first"])
    if (raw_input["nunique"] != 1).any():
        raise ValueError(
            "Aligners retaining unmapped reads disagree on total_length:\n"
            + raw_input.loc[raw_input["nunique"] != 1].to_string()
        )
    data = data.merge(
        raw_input["first"].rename("raw_input_bases").reset_index(), on=["sample", "read_technology"], how="left",
    )
    if data["raw_input_bases"].isna().any():
        raise ValueError("No aligner retains unmapped reads for at least one sample/technology.")
    if (data["bases_mapped_cigar"] > data["raw_input_bases"]).any():
        raise ValueError("bases_mapped_cigar exceeds the verified raw input for at least one row.")

    data["cigar_mapped_percent"] = data["bases_mapped_cigar"] / data["raw_input_bases"] * 100.0
    data["cigar_unaligned_percent"] = 100.0 - data["cigar_mapped_percent"]
    data["mismatch_total_percent"] = data["mismatches"] / data["raw_input_bases"] * 100.0
    data["aligned_non_mismatch_percent"] = (
        (data["bases_mapped_cigar"] - data["mismatches"]) / data["raw_input_bases"] * 100.0
    )

    plot_data = data.loc[
        data["sample"].isin(SAMPLE_ORDER)
        & data["read_technology"].isin(TECHNOLOGY_ORDER)
        & data["aligner"].isin(ALIGNER_ORDER),
        [
            "sample", "read_technology", "aligner", "total_length", "raw_input_bases", "bases_mapped_cigar",
            "mismatches", "cigar_mapped_percent", "cigar_unaligned_percent",
            "mismatch_total_percent", "aligned_non_mismatch_percent",
        ],
    ].copy()

    duplicated = plot_data.duplicated(subset=["sample", "read_technology", "aligner"], keep=False)
    if duplicated.any():
        raise ValueError("Duplicated sample/technology/aligner rows found:\n" + plot_data.loc[duplicated].to_string(index=False))

    expected_rows = len(SAMPLE_ORDER) * len(TECHNOLOGY_ORDER) * len(ALIGNER_ORDER)
    if len(plot_data) != expected_rows:
        raise ValueError(f"Expected {expected_rows} rows, found {len(plot_data)}.")

    plot_data["sample"] = pd.Categorical(plot_data["sample"], SAMPLE_ORDER, ordered=True)
    plot_data["read_technology"] = pd.Categorical(plot_data["read_technology"], TECHNOLOGY_ORDER, ordered=True)
    plot_data["aligner"] = pd.Categorical(plot_data["aligner"], ALIGNER_ORDER, ordered=True)
    return plot_data.sort_values(["read_technology", "aligner", "sample"])


def main() -> int:
    plot_data = load_data()

    OUTPUT_SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    plot_data.to_csv(OUTPUT_SUMMARY, sep="\t", index=False)

    apply_style()

    figure, axes = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=True, figsize=(HALF_WIDTH_IN, 6.3))
    x_positions = np.arange(len(ALIGNER_ORDER))

    for panel_index, technology in enumerate(TECHNOLOGY_ORDER):
        axis = axes[panel_index]
        technology_data = plot_data[plot_data["read_technology"] == technology]

        for sample in SAMPLE_ORDER:
            sample_data = (
                technology_data[technology_data["sample"] == sample].set_index("aligner").reindex(ALIGNER_ORDER)
            )
            aligned_values = sample_data["aligned_non_mismatch_percent"].to_numpy(dtype=float)
            mismatch_values = sample_data["mismatch_total_percent"].to_numpy(dtype=float)
            unaligned_values = sample_data["cigar_unaligned_percent"].to_numpy(dtype=float)
            cigar_mapped_values = sample_data["cigar_mapped_percent"].to_numpy(dtype=float)
            positions = x_positions + SAMPLE_OFFSETS[sample]

            base_color = SAMPLE_COLORS[sample]
            mismatch_color = lighten_color(base_color, amount=0.35)
            unaligned_color = lighten_color(base_color, amount=0.72)

            axis.bar(positions, aligned_values, width=BAR_WIDTH, color=base_color, edgecolor="white", linewidth=0.3, zorder=3)
            mismatch_bars = axis.bar(
                positions, mismatch_values, width=BAR_WIDTH, bottom=aligned_values,
                color=mismatch_color, edgecolor="white", linewidth=0.3, zorder=3,
            )
            axis.bar(
                positions, unaligned_values, width=BAR_WIDTH, bottom=aligned_values + mismatch_values,
                color=unaligned_color, edgecolor="white", linewidth=0.3, zorder=3,
            )

            for mismatch_bar, mismatch_value in zip(mismatch_bars, mismatch_values):
                if np.isnan(mismatch_value):
                    continue
                x_center = mismatch_bar.get_x() + mismatch_bar.get_width() / 2
                y_center = mismatch_bar.get_y() + mismatch_bar.get_height() / 2
                if mismatch_value >= 0.5:
                    axis.text(
                        x_center, y_center, f"{mismatch_value:.2f}",
                        ha="center", va="center", fontsize=PANEL_VALUE_SIZE, color="black", zorder=7,
                    )

            for position, mapped_value in zip(positions, cigar_mapped_values):
                if np.isnan(mapped_value):
                    continue
                # Default rotation_mode centers the rotated label's box over
                # its own bar instead of letting it lean onto the next one.
                axis.text(
                    position, 100.2, f"{mapped_value:.2f}",
                    ha="center", va="bottom", fontsize=PANEL_VALUE_SIZE, color=base_color,
                    rotation=45, clip_on=False, zorder=8,
                )

        axis.set_title(TECHNOLOGY_TITLES[technology], pad=20, fontsize=PANEL_TITLE_SIZE)
        rotated_xticks(axis, x_positions, ALIGNER_ORDER, fontsize=PANEL_TICK_SIZE, tick_length=2)
        axis.set_xlim(-0.5, len(ALIGNER_ORDER) - 0.5)
        axis.set_ylim(88, 100)
        axis.set_yticks([88, 90, 92, 94, 96, 98, 100])
        axis.tick_params(axis="y", labelsize=PANEL_TICK_SIZE, length=2, pad=1.5)
        clean_spines(axis)
        axis.set_facecolor("white")

    figure.supylabel("CIGAR-aligned and unaligned bases (%)", fontsize=PANEL_LABEL_SIZE_PT, x=0.01)

    panel_legend(figure, sample_legend_handles(), y=0.975)
    component_legend_handles = [
        Patch(facecolor="#555555", edgecolor="none", label="Aligned, non-mismatch"),
        Patch(facecolor="#999999", edgecolor="none", label="Mismatch"),
        Patch(facecolor="#DDDDDD", edgecolor="none", label="CIGAR-unaligned"),
    ]
    panel_legend(figure, component_legend_handles, y=0.953)

    figure.patch.set_facecolor("white")
    figure.subplots_adjust(left=0.12, right=0.99, bottom=0.067, top=0.885, hspace=0.28)

    save_figure(figure, OUTPUT_PNG, OUTPUT_PDF)

    print(f"Input: {INPUT_TSV}")
    print(f"Plot data: {OUTPUT_SUMMARY}")
    print(f"PNG: {OUTPUT_PNG}")
    print(f"PDF: {OUTPUT_PDF}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
