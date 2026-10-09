#!/usr/bin/env python3
"""Version B (individual points) of the main assembler-benchmark candidate figure.

Companion to assembler_benchmark_bars.py (Version A, grouped bars): both
read the same table (utils/benchmark_data.py), plot the same five metrics on
identical y-axis scales, and use the same x-axis geometry, so the two are a
fair side-by-side comparison rather than incidentally different. See
assembler_benchmark_bars.py's docstring for why there is no runtime panel.

Each condition shows its three individual GIAB samples as dodged points
(never averaged away) plus a short black horizontal mean marker. No
confidence interval is drawn: n = 3 specific GIAB benchmark genomes, not a
random sample, so an interval would misrepresent what these three points are.
"""

from __future__ import annotations

from pathlib import Path
import sys

PROJECT = next(
    parent for parent in Path(__file__).resolve().parents
    if (parent / "assembly_analysis").is_dir()
)
sys.path.insert(0, str(PROJECT / "assembly_analysis" / "scripts"))

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ModuleNotFoundError as error:
    raise SystemExit(
        f"A required plotting package is missing: {error.name}\n"
        f"Python interpreter: {sys.executable}\n"
    ) from error

from utils.benchmark_data import (
    METRICS,
    build_plot_data,
    compute_y_limits,
    verify_plot_data_matches_source,
    write_plot_data,
)
from utils.plot_style import (
    CONFIGURATION_ORDER,
    CONFIGURATION_POSITIONS,
    FULL_WIDTH_IN,
    PANEL_LABEL_SIZE_PT,
    PANEL_TICK_SIZE,
    SAMPLE_COLORS,
    SAMPLE_ORDER,
    apply_style,
    clean_spines,
    draw_configuration_axis,
    mean_marker_handle,
    panel_legend,
    panel_letter,
    sample_marker_handles,
    save_figure,
)

OUTPUT_PNG = PROJECT / "assembly_analysis" / "figures" / "30x" / "final" / "assembler_benchmark_points.png"
OUTPUT_PDF = OUTPUT_PNG.with_suffix(".pdf")

SAMPLE_OFFSETS = {"HG002": -0.14, "HG003": 0.0, "HG004": 0.14}
MEAN_HALF_WIDTH = 0.22
POINT_SIZE = 18

X_MIN = -0.5
X_MAX = max(CONFIGURATION_POSITIONS.values()) + 0.5


def dodged_points(axis: plt.Axes, data, metric: str) -> None:
    for assembler, technology in CONFIGURATION_ORDER:
        x = CONFIGURATION_POSITIONS[(assembler, technology)]
        subset = data[(data["assembler"] == assembler) & (data["technology"] == technology)]
        values = subset.set_index("sample")[metric].reindex(SAMPLE_ORDER)

        for sample in SAMPLE_ORDER:
            axis.scatter(
                x + SAMPLE_OFFSETS[sample], values.loc[sample],
                s=POINT_SIZE, color=SAMPLE_COLORS[sample],
                edgecolor="black", linewidth=0.4, zorder=3,
            )

        mean_value = values.mean()
        axis.plot(
            [x - MEAN_HALF_WIDTH, x + MEAN_HALF_WIDTH], [mean_value, mean_value],
            color="black", linewidth=1.4, zorder=4, solid_capstyle="butt",
        )


def style_axis(axis: plt.Axes, y_label: str, y_max: float) -> None:
    axis.set_ylim(0, y_max)
    axis.set_xlim(X_MIN, X_MAX)
    axis.tick_params(axis="y", labelsize=PANEL_TICK_SIZE, length=2, pad=1.5)
    axis.set_ylabel(y_label, fontsize=PANEL_LABEL_SIZE_PT, labelpad=2)
    clean_spines(axis)
    axis.set_facecolor("white")
    draw_configuration_axis(axis)


def main() -> int:
    data = build_plot_data()
    verify_plot_data_matches_source(data)
    write_plot_data(data)
    y_limits = compute_y_limits(data)

    apply_style()
    figure, axes = plt.subplots(nrows=2, ncols=3, figsize=(FULL_WIDTH_IN, 6.2))
    axes.flat[5].set_visible(False)

    for axis, (metric, label, letter) in zip(axes.flat, METRICS):
        dodged_points(axis, data, metric)
        style_axis(axis, label, y_limits[metric])
        panel_letter(axis, letter)

    panel_legend(
        figure, sample_marker_handles() + [mean_marker_handle()], y=0.965, ncols=4,
    )

    figure.patch.set_facecolor("white")
    figure.subplots_adjust(left=0.06, right=0.99, bottom=0.08, top=0.92, hspace=0.55, wspace=0.30)

    save_figure(figure, OUTPUT_PNG, OUTPUT_PDF)

    print(f"Rows: {len(data)}")
    print(f"PNG: {OUTPUT_PNG}")
    print(f"PDF: {OUTPUT_PDF}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
