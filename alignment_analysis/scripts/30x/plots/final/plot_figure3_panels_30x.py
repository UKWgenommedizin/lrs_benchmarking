#!/usr/bin/env python3
"""Draw the compound-figure panels for report Figure 3 (30x).

Figure 3 combines what were separate report figures (CIGAR yield vs
mismatch error, MAPQ 0 reads, memory, wall-clock runtime) plus a new
thread-hours panel. Each panel is drawn at its printed size in the same
compact style as the Figure 2 panels (ONT | PacBio HiFi side by side,
no internal panel letters); the letters and the surrounding box are added
in LaTeX. The bar panels share one sample legend, drawn once as its own
strip (fig3_legend_30x) and placed above the panels in LaTeX, so the bar
panels carry no legend of their own.

Data loading and validation are reused from the standalone scripts, so
the plotted values are identical to those figures.

Panels:
    a  input-normalized CIGAR-aligned yield vs mismatch error (full width)
    b  MAPQ 0 reads (%)
    c  configured RAM limit and measured peak RSS (GB)
    d  wall-clock runtime (h), with thread count per aligner
    e  thread-hours (allocated threads x wall-clock hours)
    f  CIGAR insertion events per 100 kb mapped
    g  CIGAR deletion events per 100 kb mapped
    legend  shared sample + RAM-limit legend for panels b-e
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
from matplotlib.lines import Line2D
from scipy import stats

import plot_cigar_yield_vs_error_30x as cigar_scatter
import plot_indels_30x as indels
import plot_memory_30x as memory
import plot_mq0_reads_30x as mq0
import plot_runtime_30x as runtime
from utils.plot_style import (
    ALIGNER_MARKERS,
    ALIGNER_ORDER,
    HALF_WIDTH_IN,
    NEUTRAL_GRAY,
    PANEL_LABEL_SIZE_PT,
    PANEL_TICK_SIZE,
    PANEL_TITLE_SIZE,
    SAMPLE_COLORS,
    SAMPLE_ORDER,
    TECHNOLOGY_ORDER,
    TECHNOLOGY_TITLES,
    aligner_marker_handles,
    apply_style,
    clean_spines,
    panel_legend,
    rotated_xticks,
    sample_legend_handles,
    sample_marker_handles,
    save_figure,
    subtle_grid,
)

FIGURE_DIR = PROJECT / "alignment_analysis" / "figures" / "30x" / "final"
PLOT_DATA_DIR = PROJECT / "alignment_analysis" / "tables" / "30x" / "derived" / "plot_data"

OUTPUT_THREAD_HOURS = PLOT_DATA_DIR / "thread_hours_30x_plotting_values.tsv"

# Full-width panel spans the report's text block inside the Figure 3 box.
FULL_PANEL_WIDTH_IN = 6.3
BAR_PANEL_HEIGHT_IN = 1.9

# Bar panels print at ~0.47 of the text block, so they use larger fonts than
# the shared PANEL_* sizes (which Figure 2 also uses) to stay readable.
BAR_TICK_SIZE = 7.0
BAR_LABEL_SIZE = 7.5
BAR_TITLE_SIZE = 7.5
BAR_VALUE_SIZE = 6.0
BAR_LEGEND_SIZE = 7.0

SAMPLE_OFFSETS = {"HG002": -0.27, "HG003": 0.0, "HG004": 0.27}
BAR_WIDTH = 0.26
LIMIT_COLOR = "#222222"


def output_paths(stem: str) -> tuple[Path, Path]:
    png = FIGURE_DIR / f"{stem}.png"
    return png, png.with_suffix(".pdf")


def grouped_bars(axis: plt.Axes, subset: pd.DataFrame, value_column: str, aligners=ALIGNER_ORDER) -> None:
    for aligner_index, aligner in enumerate(aligners):
        values = (
            subset[subset["aligner"] == aligner]
            .set_index("sample")[value_column]
            .reindex(SAMPLE_ORDER)
        )
        for sample in SAMPLE_ORDER:
            value = values.loc[sample]
            if not np.isfinite(value):
                continue
            axis.bar(
                aligner_index + SAMPLE_OFFSETS[sample], value,
                width=BAR_WIDTH, color=SAMPLE_COLORS[sample],
                edgecolor="white", linewidth=0.3, zorder=3,
            )


def style_bar_axis(axis: plt.Axes, technology: str, y_max: float) -> None:
    axis.set_title(TECHNOLOGY_TITLES[technology], pad=2, fontsize=BAR_TITLE_SIZE)
    axis.set_ylim(0, y_max)
    axis.set_xlim(-0.5, len(ALIGNER_ORDER) - 0.5)
    rotated_xticks(axis, range(len(ALIGNER_ORDER)), ALIGNER_ORDER, fontsize=BAR_TICK_SIZE, tick_length=2)
    axis.tick_params(axis="y", labelsize=BAR_TICK_SIZE, length=2, pad=1.5)
    clean_spines(axis)
    axis.set_facecolor("white")


def finish_bar_figure(figure: plt.Figure, stem: str) -> None:
    figure.patch.set_facecolor("white")
    figure.subplots_adjust(left=0.13, right=0.99, bottom=0.27, top=0.90, wspace=0.08)
    save_figure(figure, *output_paths(stem))


def ram_limit_handle() -> Line2D:
    return Line2D(
        [0], [0], color=LIMIT_COLOR, linestyle=(0, (2.5, 1.5)), linewidth=0.8, label="RAM limit (d)",
    )


def plot_shared_legend() -> None:
    """One legend strip for all bar panels, placed above them in LaTeX."""
    figure = plt.figure(figsize=(FULL_PANEL_WIDTH_IN, 0.25))
    figure.legend(
        handles=sample_legend_handles() + [ram_limit_handle()],
        frameon=False, ncols=4, loc="center", fontsize=BAR_LEGEND_SIZE,
        handlelength=1.2, handleheight=0.9, handletextpad=0.4, columnspacing=1.4,
        borderaxespad=0, borderpad=0,
    )
    figure.patch.set_facecolor("white")
    save_figure(figure, *output_paths("fig3_legend_30x"))


# ------------------------------------------------------------
# a: CIGAR-aligned yield vs mismatch error
# ------------------------------------------------------------


def plot_cigar_scatter() -> None:
    plot_data = cigar_scatter.load_data()
    figure, axes = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(FULL_PANEL_WIDTH_IN, 2.35))

    for axis, technology in zip(axes, TECHNOLOGY_ORDER):
        technology_data = plot_data[plot_data["read_technology"] == technology]
        for _, row in technology_data.iterrows():
            axis.scatter(
                row["input_normalized_cigar_yield_percent"], row["error_percent"],
                marker=ALIGNER_MARKERS[str(row["aligner"])], color=SAMPLE_COLORS[str(row["sample"])],
                edgecolor="black", linewidth=0.4, s=16, zorder=3,
            )

        spearman = stats.spearmanr(
            technology_data["input_normalized_cigar_yield_percent"], technology_data["error_percent"]
        )
        axis.text(
            0.03, 0.96,
            f"Spearman $\\rho$ = {spearman.statistic:.2f}, $P$ = {spearman.pvalue:.3f}\n$n$ = {len(technology_data)}",
            transform=axis.transAxes, ha="left", va="top", fontsize=PANEL_TICK_SIZE,
        )
        axis.set_title(TECHNOLOGY_TITLES[technology], pad=2, fontsize=PANEL_TITLE_SIZE)
        axis.tick_params(axis="both", labelsize=PANEL_TICK_SIZE, length=2, pad=1.5)
        subtle_grid(axis, "both")
        clean_spines(axis)
        axis.set_facecolor("white")

    axes[0].set_ylabel("Mismatch error (%)", fontsize=PANEL_LABEL_SIZE_PT, labelpad=2)
    figure.supxlabel("Input-normalized CIGAR-aligned yield (%)", fontsize=PANEL_LABEL_SIZE_PT, y=0.04)

    handles = sample_marker_handles() + aligner_marker_handles()
    for handle in handles:
        handle.set_markersize(4)
    panel_legend(figure, handles, y=0.92)

    figure.patch.set_facecolor("white")
    figure.subplots_adjust(left=0.07, right=0.99, bottom=0.16, top=0.83, wspace=0.08)
    save_figure(figure, *output_paths("fig3a_cigar_yield_vs_error_30x"))


# ------------------------------------------------------------
# b: MAPQ 0 reads
# ------------------------------------------------------------


def plot_mq0() -> None:
    plot_data = mq0.load_data()
    y_max = float(plot_data["reads_mq0_percent"].max()) * 1.08
    figure, axes = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(HALF_WIDTH_IN, BAR_PANEL_HEIGHT_IN))

    for axis, technology in zip(axes, TECHNOLOGY_ORDER):
        subset = plot_data[plot_data["read_technology"] == technology]
        grouped_bars(axis, subset, "reads_mq0_percent")
        style_bar_axis(axis, technology, y_max)

    axes[0].set_ylabel("MQ0 reads (%)", fontsize=BAR_LABEL_SIZE, labelpad=2)
    finish_bar_figure(figure, "fig3b_mq0_reads_30x")


# ------------------------------------------------------------
# c: configured RAM limit + measured peak RSS
# ------------------------------------------------------------


def plot_memory() -> None:
    data = memory.load_data()
    memory.prepare_plot_data(data)  # validates that each aligner has one configured limit
    limits = data.groupby("aligner")["command_ram_limit_gb"].first().reindex(ALIGNER_ORDER)

    y_max = float(limits.max()) * 1.12
    figure, axes = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(HALF_WIDTH_IN, BAR_PANEL_HEIGHT_IN))

    for axis, technology in zip(axes, TECHNOLOGY_ORDER):
        subset = data[data["read_technology"] == technology]
        grouped_bars(axis, subset, "peak_ram_gb")

        for aligner_index, aligner in enumerate(ALIGNER_ORDER):
            limit = limits.loc[aligner]
            axis.hlines(
                limit, aligner_index - 0.42, aligner_index + 0.42,
                colors=LIMIT_COLOR, linestyles=(0, (2.5, 1.5)), linewidth=0.8, zorder=4,
            )
            axis.text(
                aligner_index, limit + y_max * 0.015, f"{limit:.0f}",
                ha="center", va="bottom", fontsize=BAR_VALUE_SIZE, color=LIMIT_COLOR,
            )
            if subset.loc[subset["aligner"] == aligner, "peak_ram_gb"].isna().all():
                axis.text(
                    aligner_index, y_max * 0.02, "n.m.",
                    ha="center", va="bottom", fontsize=BAR_VALUE_SIZE, color=NEUTRAL_GRAY,
                )

        style_bar_axis(axis, technology, y_max)

    axes[0].set_ylabel("Memory (GB)", fontsize=BAR_LABEL_SIZE, labelpad=2)
    finish_bar_figure(figure, "fig3c_memory_30x")


# ------------------------------------------------------------
# d/e: wall-clock runtime and thread-hours
# ------------------------------------------------------------


def annotate_threads(axis: plt.Axes, subset: pd.DataFrame, value_column: str, y_max: float) -> None:
    for aligner_index, aligner in enumerate(ALIGNER_ORDER):
        aligner_data = subset[subset["aligner"] == aligner]
        threads = "/".join(str(t) for t in sorted(aligner_data["threads"].unique()))
        axis.text(
            aligner_index, aligner_data[value_column].max() + y_max * 0.02, f"{threads} thr.",
            ha="center", va="bottom", fontsize=BAR_VALUE_SIZE, color="#555555",
        )


def plot_runtime_and_thread_hours() -> None:
    data = runtime.load_data()
    data["thread_hours"] = data["threads"] * data["wallclock_runtime_hours"]

    export = data[["sample", "read_technology", "aligner", "threads", "wallclock_runtime_hours", "thread_hours"]].copy()
    export["read_technology"] = pd.Categorical(export["read_technology"], TECHNOLOGY_ORDER, ordered=True)
    export["aligner"] = pd.Categorical(export["aligner"], ALIGNER_ORDER, ordered=True)
    export = export.sort_values(["read_technology", "aligner", "sample"])
    OUTPUT_THREAD_HOURS.parent.mkdir(parents=True, exist_ok=True)
    export.to_csv(OUTPUT_THREAD_HOURS, sep="\t", index=False)

    panels = [
        ("wallclock_runtime_hours", "Wall-clock runtime (h)", "fig3d_runtime_30x", True),
        ("thread_hours", "Thread-hours", "fig3e_thread_hours_30x", False),
    ]
    for value_column, label, stem, show_threads in panels:
        y_max = float(data[value_column].max()) * 1.15
        figure, axes = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(HALF_WIDTH_IN, BAR_PANEL_HEIGHT_IN))

        for axis, technology in zip(axes, TECHNOLOGY_ORDER):
            subset = data[data["read_technology"] == technology]
            grouped_bars(axis, subset, value_column)
            if show_threads:
                annotate_threads(axis, subset, value_column, y_max)
            style_bar_axis(axis, technology, y_max)

        axes[0].set_ylabel(label, fontsize=BAR_LABEL_SIZE, labelpad=2)
        finish_bar_figure(figure, stem)


# ------------------------------------------------------------
# f/g: CIGAR insertion and deletion events
# ------------------------------------------------------------


def plot_indel_events() -> None:
    data = indels.load_data()
    panels = [
        ("insertion_events_per_100kb", "Insertions per 100 kb", "fig3f_insertions_30x"),
        ("deletion_events_per_100kb", "Deletions per 100 kb", "fig3g_deletions_30x"),
    ]
    # One y range for both panels so insertion and deletion rates compare directly.
    y_max = float(data[[column for column, _, _ in panels]].max().max()) * 1.08
    for value_column, label, stem in panels:
        figure, axes = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(HALF_WIDTH_IN, BAR_PANEL_HEIGHT_IN))

        for axis, technology in zip(axes, TECHNOLOGY_ORDER):
            subset = data[data["read_technology"] == technology]
            grouped_bars(axis, subset, value_column)
            for aligner_index, aligner in enumerate(ALIGNER_ORDER):
                if subset.loc[subset["aligner"] == aligner, value_column].isna().any():
                    axis.text(
                        aligner_index, y_max * 0.02, "n/a",
                        ha="center", va="bottom", fontsize=BAR_VALUE_SIZE, color=NEUTRAL_GRAY,
                    )
            style_bar_axis(axis, technology, y_max)

        axes[0].set_ylabel(label, fontsize=BAR_LABEL_SIZE, labelpad=2)
        finish_bar_figure(figure, stem)


def main() -> int:
    apply_style()
    plot_cigar_scatter()
    plot_mq0()
    plot_memory()
    plot_runtime_and_thread_hours()
    plot_indel_events()
    plot_shared_legend()

    for stem in [
        "fig3a_cigar_yield_vs_error_30x", "fig3b_mq0_reads_30x", "fig3c_memory_30x",
        "fig3d_runtime_30x", "fig3e_thread_hours_30x", "fig3f_insertions_30x",
        "fig3g_deletions_30x", "fig3_legend_30x",
    ]:
        print(f"PDF: {output_paths(stem)[1]}")
    print(f"Plot data: {OUTPUT_THREAD_HOURS}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
