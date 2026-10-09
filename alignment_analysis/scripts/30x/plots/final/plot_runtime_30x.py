#!/usr/bin/env python3
"""Plot observed wall-clock runtime (hours) for the 30x benchmark aligners.

Grouped bar chart (one bar per HG002/HG003/HG004 sample within each
aligner) from a true zero baseline, styled like the other half-width
30x panels (e.g. 01_alignment_error_rate_30x).
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
    apply_style,
    clean_spines,
    panel_legend,
    rotated_xticks,
    sample_legend_handles,
    save_figure,
)

TABLE_DIR = PROJECT / "alignment_analysis" / "tables" / "30x"
FIGURE_DIR = PROJECT / "alignment_analysis" / "figures" / "30x" / "final"

INPUT_TSV = TABLE_DIR / "final" / "alignment_benchmark_30x.tsv"
OUTPUT_DATA = TABLE_DIR / "derived" / "plot_data" / "wallclock_runtime_30x_plotting_values.tsv"
OUTPUT_PNG = FIGURE_DIR / "05_runtime_30x.png"
OUTPUT_PDF = OUTPUT_PNG.with_suffix(".pdf")

BAR_WIDTH = 0.26
SAMPLE_OFFSETS = {"HG002": -0.27, "HG003": 0.0, "HG004": 0.27}


def load_data() -> pd.DataFrame:
    if not INPUT_TSV.is_file():
        raise FileNotFoundError(f"Benchmark table not found: {INPUT_TSV}")

    data = pd.read_csv(INPUT_TSV, sep="\t")
    required = {
        "sample", "read_technology", "aligner", "configuration",
        "wallclock_runtime_hours", "threads", "wallclock_runtime_source_log",
    }
    missing = required - set(data.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    data["aligner"] = data["aligner"].replace({"VACMap": "VACmap", "VG": "VG Giraffe", "vg": "VG Giraffe"})
    data["read_technology"] = data["read_technology"].replace(
        {"PacBio HiFi": "PacBio", "PB": "PacBio", "pb": "PacBio", "ont": "ONT"}
    )
    data = data[
        data["sample"].isin(SAMPLE_ORDER)
        & data["read_technology"].isin(TECHNOLOGY_ORDER)
        & data["aligner"].isin(ALIGNER_ORDER)
    ].copy()

    data["wallclock_runtime_hours"] = pd.to_numeric(data["wallclock_runtime_hours"], errors="raise")
    data["threads"] = pd.to_numeric(data["threads"], errors="raise").astype(int)
    if data["wallclock_runtime_hours"].isna().any():
        rows = data.loc[data["wallclock_runtime_hours"].isna(), ["sample", "read_technology", "aligner"]]
        raise ValueError("Missing runtime values:\n" + rows.to_string(index=False))
    if (data["wallclock_runtime_hours"] < 0).any():
        raise ValueError("Wall-clock runtime cannot be negative")

    keys = ["sample", "read_technology", "aligner"]
    duplicates = data.duplicated(keys, keep=False)
    if duplicates.any():
        raise ValueError("Duplicate sample/technology/aligner rows:\n" + data.loc[duplicates, keys].to_string(index=False))

    expected = pd.MultiIndex.from_product([SAMPLE_ORDER, TECHNOLOGY_ORDER, ALIGNER_ORDER], names=keys)
    observed = pd.MultiIndex.from_frame(data[keys])
    missing_combinations = expected.difference(observed)
    if len(missing_combinations):
        raise ValueError(f"Missing benchmark combinations: {list(missing_combinations)}")

    if "runtime_configuration" in data.columns:
        runtime_configuration = data["runtime_configuration"].astype("string")
        data["runtime_matches_benchmark_configuration"] = runtime_configuration.eq(data["configuration"].astype("string"))
    else:
        data["runtime_configuration"] = pd.NA
        data["runtime_matches_benchmark_configuration"] = pd.NA

    return data


def add_panel(axis: plt.Axes, subset: pd.DataFrame, y_max: float) -> None:
    x_positions = np.arange(len(ALIGNER_ORDER), dtype=float)

    for sample in SAMPLE_ORDER:
        sample_data = (
            subset[subset["sample"] == sample].set_index("aligner").reindex(ALIGNER_ORDER)
        )
        values = sample_data["wallclock_runtime_hours"].to_numpy(dtype=float)
        valid = np.isfinite(values)
        axis.bar(
            x_positions[valid] + SAMPLE_OFFSETS[sample], values[valid],
            width=BAR_WIDTH,
            color=SAMPLE_COLORS[sample], edgecolor="white", linewidth=0.3,
            zorder=3,
        )

    # runtimes are only comparable alongside the CPU allocation, so note it above each group
    for x, aligner in zip(x_positions, ALIGNER_ORDER):
        aligner_data = subset[subset["aligner"] == aligner]
        threads = "/".join(str(t) for t in sorted(aligner_data["threads"].unique()))
        axis.text(
            x, aligner_data["wallclock_runtime_hours"].max() + y_max * 0.02, f"{threads} thr.",
            ha="center", va="bottom", fontsize=PANEL_VALUE_SIZE, color="#555555",
        )

    axis.set_ylim(0, y_max)
    axis.set_xlim(-0.5, len(ALIGNER_ORDER) - 0.5)
    rotated_xticks(axis, x_positions, ALIGNER_ORDER, fontsize=PANEL_TICK_SIZE, tick_length=2)
    axis.tick_params(axis="y", labelsize=PANEL_TICK_SIZE, length=2, pad=1.5)
    clean_spines(axis)
    axis.set_facecolor("white")


def main() -> int:
    data = load_data()
    output_columns = [
        "sample", "read_technology", "aligner", "configuration",
        "runtime_configuration", "runtime_matches_benchmark_configuration",
        "wallclock_runtime_hours", "wallclock_runtime_seconds",
        "wallclock_runtime_hms", "wallclock_runtime_source_log", "threads",
    ]
    output_columns = [column for column in output_columns if column in data.columns]
    plot_data = data[output_columns].copy()
    plot_data["read_technology"] = pd.Categorical(plot_data["read_technology"], TECHNOLOGY_ORDER, ordered=True)
    plot_data["aligner"] = pd.Categorical(plot_data["aligner"], ALIGNER_ORDER, ordered=True)
    plot_data["sample"] = pd.Categorical(plot_data["sample"], SAMPLE_ORDER, ordered=True)
    plot_data = plot_data.sort_values(["read_technology", "aligner", "sample"])

    OUTPUT_DATA.parent.mkdir(parents=True, exist_ok=True)
    plot_data.to_csv(OUTPUT_DATA, sep="\t", index=False)

    apply_style()

    figure, axes = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(HALF_WIDTH_IN, 1.75))
    y_max = float(plot_data["wallclock_runtime_hours"].max()) * 1.15

    for axis, technology in zip(axes, TECHNOLOGY_ORDER):
        subset = plot_data[plot_data["read_technology"] == technology]
        add_panel(axis, subset, y_max)
        axis.set_title(TECHNOLOGY_TITLES[technology], pad=2, fontsize=PANEL_TITLE_SIZE)

    axes[0].set_ylabel("Wall-clock runtime (h)", fontsize=PANEL_LABEL_SIZE_PT, labelpad=2)

    panel_legend(figure, sample_legend_handles(), y=0.93)

    figure.patch.set_facecolor("white")
    figure.subplots_adjust(left=0.12, right=0.99, bottom=0.25, top=0.82, wspace=0.08)

    save_figure(figure, OUTPUT_PNG, OUTPUT_PDF)

    print(f"Input: {INPUT_TSV}")
    print(f"Plot data: {OUTPUT_DATA}")
    print(f"PNG: {OUTPUT_PNG}")
    print(f"PDF: {OUTPUT_PDF}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
