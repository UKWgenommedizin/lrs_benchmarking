#!/usr/bin/env python3
"""Draw report Figure 2 (30x) as a single compound figure.

Horizontal bars with one shared aligner axis per sequencing technology:
rows are technologies (ONT, PacBio HiFi), columns are metrics, so each
aligner reads straight across all four bar panels. The CIGAR-composition
column is widest so its aligned / mismatch / unaligned segments stay
legible, and the CIGAR-aligned yield vs mismatch error scatter sits
underneath it as panel e, closing the figure on the relationship the bar
panels build up to.

Panels:
    a  alignment error rate (%)
    b  mapped / unmapped reads (%), truncated 88-100% baseline
    c  mapped / unmapped bases (%), truncated 88-100% baseline
    d  CIGAR composition: aligned non-mismatch / mismatch / unaligned (%)
    e  input-normalized CIGAR-aligned yield vs mismatch error

Data loading and validation are reused from the standalone scripts, so
the plotted values are identical to those figures. Panel letters are drawn
here; the report only adds the surrounding box and caption.
"""

from __future__ import annotations

from pathlib import Path
import sys

PROJECT = next(
    parent for parent in Path(__file__).resolve().parents
    if (parent / "alignment_analysis").is_dir()
)
sys.path.insert(0, str(PROJECT / "alignment_analysis" / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.patches import Patch

import plot_alignment_error_rate_30x as error_rate
import plot_cigar_composition_30x as cigar_composition
import plot_cigar_yield_vs_error_30x as cigar_scatter
import plot_mapped_unmapped_bases_30x as mapped_bases
import plot_mapped_unmapped_reads_30x as mapped_reads
from utils.plot_style import (
    ALIGNER_MARKERS,
    ALIGNER_ORDER,
    PANEL_LABEL_SIZE,
    PANEL_LABEL_SIZE_PT,
    PANEL_TICK_SIZE,
    PANEL_TITLE_SIZE,
    PANEL_VALUE_SIZE,
    SAMPLE_COLORS,
    SAMPLE_ORDER,
    TECHNOLOGY_ORDER,
    TECHNOLOGY_TITLES,
    aligner_marker_handles,
    apply_style,
    clean_spines,
    panel_legend,
    sample_marker_handles,
    save_figure,
    subtle_grid,
)

FIGURE_DIR = PROJECT / "alignment_analysis" / "figures" / "30x" / "final"
OUTPUT_PNG = FIGURE_DIR / "fig2_alignment_performance_30x.png"
OUTPUT_PDF = OUTPUT_PNG.with_suffix(".pdf")

# Legend text at 7 pt, matching Figures 3 and 4 (Nature: 5-7 pt figure text).
FIGURE2_LEGEND_SIZE = 7.0

FIGURE_WIDTH_IN = 6.3
FIGURE_HEIGHT_IN = 6.7

# HG002 on top within each aligner group (the y axis is inverted).
SAMPLE_OFFSETS = {"HG002": -0.27, "HG003": 0.0, "HG004": 0.27}
BAR_HEIGHT = 0.26

MISMATCH_LIGHTEN = 0.35
SECONDARY_LIGHTEN = 0.72

# column key -> (panel letter, x label, x limits, x ticks)
COLUMNS = {
    "error": ("a", "Alignment error rate (%)", None, [0, 1, 2, 3]),
    "reads": ("b", "Mapped reads (%)", (88, 100), [88, 92, 96, 100]),
    "bases": ("c", "Mapped bases (%)", (88, 100), [88, 92, 96, 100]),
    "cigar": ("d", "CIGAR-aligned bases (%)", (88, 100), [88, 92, 96, 100]),
}
COLUMN_WIDTHS = [1.0, 1.0, 1.0, 1.9]


def lighten_color(color, amount: float) -> tuple[float, float, float]:
    rgb = np.array(mcolors.to_rgb(color))
    return tuple(rgb + (np.ones(3) - rgb) * amount)


def per_sample(subset: pd.DataFrame, sample: str, column: str) -> np.ndarray:
    return (
        subset[subset["sample"] == sample]
        .set_index("aligner")[column]
        .reindex(ALIGNER_ORDER)
        .to_numpy(dtype=float)
    )


def axis_width_points(figure: plt.Figure, axis: plt.Axes) -> float:
    return axis.get_position().width * figure.get_figwidth() * 72.0


def label_fits(figure: plt.Figure, axis: plt.Axes, span: float, text: str, fontsize: float) -> bool:
    low, high = axis.get_xlim()
    segment_points = span / (high - low) * axis_width_points(figure, axis)
    return segment_points >= len(text) * fontsize * 0.52 + 1.0


# ------------------------------------------------------------
# a-d: horizontal bar panels
# ------------------------------------------------------------


def draw_error(axis: plt.Axes, subset: pd.DataFrame, x_max: float) -> None:
    y_positions = np.arange(len(ALIGNER_ORDER))
    for sample in SAMPLE_ORDER:
        axis.barh(
            y_positions + SAMPLE_OFFSETS[sample], per_sample(subset, sample, error_rate.METRIC),
            height=BAR_HEIGHT, color=SAMPLE_COLORS[sample], edgecolor="white", linewidth=0.3, zorder=3,
        )
    axis.set_xlim(0, x_max)


def draw_mapped_unmapped(axis: plt.Axes, subset: pd.DataFrame, mapped_column: str, unmapped_column: str) -> None:
    y_positions = np.arange(len(ALIGNER_ORDER))
    for sample in SAMPLE_ORDER:
        positions = y_positions + SAMPLE_OFFSETS[sample]
        mapped = per_sample(subset, sample, mapped_column)
        unmapped = per_sample(subset, sample, unmapped_column)
        color = SAMPLE_COLORS[sample]
        axis.barh(positions, mapped, height=BAR_HEIGHT, color=color, edgecolor="white", linewidth=0.3, zorder=3)
        axis.barh(
            positions, unmapped, height=BAR_HEIGHT, left=mapped,
            color=lighten_color(color, SECONDARY_LIGHTEN), edgecolor="white", linewidth=0.3, zorder=3,
        )


def draw_cigar(figure: plt.Figure, axis: plt.Axes, subset: pd.DataFrame) -> None:
    y_positions = np.arange(len(ALIGNER_ORDER))
    for sample in SAMPLE_ORDER:
        positions = y_positions + SAMPLE_OFFSETS[sample]
        aligned = per_sample(subset, sample, "aligned_non_mismatch_percent")
        mismatch = per_sample(subset, sample, "mismatch_total_percent")
        unaligned = per_sample(subset, sample, "cigar_unaligned_percent")
        cigar_mapped = per_sample(subset, sample, "cigar_mapped_percent")
        color = SAMPLE_COLORS[sample]

        axis.barh(positions, aligned, height=BAR_HEIGHT, color=color, edgecolor="white", linewidth=0.3, zorder=3)
        axis.barh(
            positions, mismatch, height=BAR_HEIGHT, left=aligned,
            color=lighten_color(color, MISMATCH_LIGHTEN), edgecolor="white", linewidth=0.3, zorder=3,
        )
        axis.barh(
            positions, unaligned, height=BAR_HEIGHT, left=aligned + mismatch,
            color=lighten_color(color, SECONDARY_LIGHTEN), edgecolor="white", linewidth=0.3, zorder=3,
        )

        for position, left, mismatch_value, mapped_value in zip(positions, aligned, mismatch, cigar_mapped):
            if np.isfinite(mismatch_value):
                text = f"{mismatch_value:.2f}"
                if label_fits(figure, axis, mismatch_value, text, PANEL_VALUE_SIZE - 0.5):
                    axis.text(
                        left + mismatch_value / 2, position, text,
                        ha="center", va="center", fontsize=PANEL_VALUE_SIZE - 0.5, color="black", zorder=7,
                    )
                else:
                    # Segment too thin for its value: same black label on a
                    # white tag over the dark segment, with a leader line
                    # pointing at the mismatch segment it belongs to.
                    axis.annotate(
                        text, xy=(left + mismatch_value / 2, position), xytext=(left - 1.0, position),
                        ha="right", va="center", fontsize=PANEL_VALUE_SIZE - 0.5, color="black", zorder=7,
                        bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="none"),
                        arrowprops=dict(arrowstyle="-", color="black", linewidth=0.5, shrinkA=0, shrinkB=0),
                    )
            if np.isfinite(mapped_value):
                axis.text(
                    100.25, position, f"{mapped_value:.2f}",
                    ha="left", va="center", fontsize=PANEL_VALUE_SIZE - 0.5, color=color,
                    clip_on=False, zorder=8,
                )


def style_bar_axis(axis: plt.Axes, column: str, show_aligners: bool, show_xlabel: bool) -> None:
    _, x_label, x_limits, x_ticks = COLUMNS[column]
    if x_limits is not None:
        axis.set_xlim(*x_limits)
    if x_ticks is not None:
        axis.set_xticks(x_ticks)
    axis.set_ylim(len(ALIGNER_ORDER) - 0.5, -0.5)
    axis.set_yticks(range(len(ALIGNER_ORDER)), ALIGNER_ORDER)
    axis.tick_params(axis="both", labelsize=PANEL_TICK_SIZE, length=2, pad=1.5)
    if not show_aligners:
        # The y axis is shared along the row: hide the repeated names only.
        axis.tick_params(axis="y", length=0, labelleft=False)
    if show_xlabel:
        axis.set_xlabel(x_label, fontsize=PANEL_LABEL_SIZE_PT, labelpad=2)
    clean_spines(axis)
    axis.set_facecolor("white")


def draw_bar_block(figure: plt.Figure, grid) -> list[list[plt.Axes]]:
    error_data = error_rate.load_plot_data()
    reads_data = mapped_reads.load_data()
    bases_data = mapped_bases.load_data()
    cigar_data = cigar_composition.load_data()
    error_max = float(error_data[error_rate.METRIC].max()) * 1.08

    axes_grid = []
    for row_index, technology in enumerate(TECHNOLOGY_ORDER):
        row_axes = []
        for column_index, column in enumerate(COLUMNS):
            share_x = axes_grid[0][column_index] if row_index else None
            share_y = row_axes[0] if column_index else None
            axis = figure.add_subplot(grid[row_index, column_index], sharex=share_x, sharey=share_y)
            row_axes.append(axis)
        axes_grid.append(row_axes)

        is_bottom = row_index == len(TECHNOLOGY_ORDER) - 1
        for column_index, column in enumerate(COLUMNS):
            style_bar_axis(row_axes[column_index], column, show_aligners=column_index == 0, show_xlabel=is_bottom)
            # sharex hides nothing by itself; keep tick labels on both rows
            # because the truncated baselines need to be readable in each.
            row_axes[column_index].tick_params(axis="x", labelbottom=True)

        draw_error(row_axes[0], error_data[error_data["read_technology"] == technology], error_max)
        draw_mapped_unmapped(
            row_axes[1], reads_data[reads_data["read_technology"] == technology],
            "mapped_reads_percent", "unmapped_reads_percent",
        )
        draw_mapped_unmapped(
            row_axes[2], bases_data[bases_data["read_technology"] == technology],
            "mapped_bases_percent", "unmapped_bases_percent",
        )
        draw_cigar(figure, row_axes[3], cigar_data[cigar_data["read_technology"] == technology])

        # Technology label on the left, outside the aligner names.
        first = row_axes[0].get_position()
        figure.text(
            0.012, first.y0 + first.height / 2, TECHNOLOGY_TITLES[technology],
            rotation=90, ha="center", va="center", fontsize=PANEL_TITLE_SIZE + 0.5, fontweight="bold",
        )

    return axes_grid


# ------------------------------------------------------------
# e: CIGAR-aligned yield vs mismatch error
# ------------------------------------------------------------


def draw_scatter(figure: plt.Figure, grid) -> list[plt.Axes]:
    plot_data = cigar_scatter.load_data()
    axes = [figure.add_subplot(grid[0, 0])]
    axes.append(figure.add_subplot(grid[0, 1], sharey=axes[0]))

    for axis, technology in zip(axes, TECHNOLOGY_ORDER):
        technology_data = plot_data[plot_data["read_technology"] == technology]
        for _, row in technology_data.iterrows():
            axis.scatter(
                row["input_normalized_cigar_yield_percent"], row["error_percent"],
                marker=ALIGNER_MARKERS[str(row["aligner"])], color=SAMPLE_COLORS[str(row["sample"])],
                edgecolor="black", linewidth=0.4, s=16, zorder=3,
            )

        axis.text(
            0.03, 0.96,
            f"$n$ = {len(technology_data)}",
            transform=axis.transAxes, ha="left", va="top", fontsize=PANEL_TICK_SIZE,
        )
        axis.set_title(TECHNOLOGY_TITLES[technology], pad=2, fontsize=PANEL_TITLE_SIZE)
        axis.set_xlabel("Input-normalized CIGAR-aligned yield (%)", fontsize=PANEL_LABEL_SIZE_PT, labelpad=2)
        axis.tick_params(axis="both", labelsize=PANEL_TICK_SIZE, length=2, pad=1.5)
        subtle_grid(axis, "both")
        clean_spines(axis)
        axis.set_facecolor("white")

    axes[1].tick_params(axis="y", labelleft=False)
    axes[0].set_ylabel("Mismatch error (%)", fontsize=PANEL_LABEL_SIZE_PT, labelpad=2)
    return axes


# ------------------------------------------------------------
# assembly
# ------------------------------------------------------------


def panel_letter_at(figure: plt.Figure, axis: plt.Axes, letter: str, x: float, dy: float = 0.012) -> None:
    figure.text(
        x, axis.get_position().y1 + dy, letter,
        fontsize=PANEL_LABEL_SIZE + 1, fontweight="bold", ha="left", va="bottom",
    )


def main() -> int:
    apply_style()
    figure = plt.figure(figsize=(FIGURE_WIDTH_IN, FIGURE_HEIGHT_IN))
    figure.patch.set_facecolor("white")

    outer = GridSpec(
        2, 1, figure=figure, height_ratios=[2.3, 0.95], hspace=0.30,
        left=0.125, right=0.945, top=0.905, bottom=0.055,
    )
    bar_grid = GridSpecFromSubplotSpec(
        2, 4, subplot_spec=outer[0], width_ratios=COLUMN_WIDTHS, hspace=0.20, wspace=0.14,
    )
    scatter_grid = GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[1], wspace=0.06)

    bar_axes = draw_bar_block(figure, bar_grid)
    scatter_axes = draw_scatter(figure, scatter_grid)

    for column_index, (letter, *_rest) in enumerate(COLUMNS.values()):
        axis = bar_axes[0][column_index]
        x = 0.02 if column_index == 0 else axis.get_position().x0 - 0.012
        panel_letter_at(figure, axis, letter, x)
    panel_letter_at(figure, scatter_axes[0], "e", 0.02, dy=0.03)

    # Legends: sample colour + segment shade for a-d, sample colour +
    # aligner shape for e.
    grey = "#555555"
    shade_handles = [
        Patch(facecolor=grey, edgecolor="none", label="Mapped (b, c) / aligned, non-mismatch (d)"),
        Patch(facecolor=lighten_color(grey, MISMATCH_LIGHTEN), edgecolor="none", label="Mismatch (d)"),
        Patch(facecolor=lighten_color(grey, SECONDARY_LIGHTEN), edgecolor="none",
              label="Unmapped (b, c) / CIGAR-unaligned (d)"),
    ]
    sample_handles = [
        Patch(facecolor=SAMPLE_COLORS[sample], edgecolor="none", label=sample) for sample in SAMPLE_ORDER
    ]
    panel_legend(figure, sample_handles, y=0.968, fontsize=FIGURE2_LEGEND_SIZE)
    panel_legend(figure, shade_handles, y=0.945, fontsize=FIGURE2_LEGEND_SIZE)

    scatter_top = scatter_axes[0].get_position().y1
    marker_handles = sample_marker_handles() + aligner_marker_handles()
    for handle in marker_handles:
        handle.set_markersize(4)
    panel_legend(figure, marker_handles, y=scatter_top + 0.022, fontsize=FIGURE2_LEGEND_SIZE)

    save_figure(figure, OUTPUT_PNG, OUTPUT_PDF)
    print(f"PNG: {OUTPUT_PNG}")
    print(f"PDF: {OUTPUT_PDF}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
