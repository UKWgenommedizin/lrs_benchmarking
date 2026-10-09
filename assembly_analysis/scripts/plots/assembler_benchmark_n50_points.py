#!/usr/bin/env python3
"""Candidate A2: assembler_benchmark_bars.py with panel (a) as log-scale points.

Companion to assembler_benchmark_n50_logbars.py (candidate A1, same log
scale, grouped bars instead of points for panel a) -- see that script's
docstring for why N50 needs a log scale at all. Panels b-f are identical to
assembler_benchmark_bars.py (reused directly). Panel (a) shows the three
individual HG002/HG003/HG004 observations per condition as dodged points
(reusing assembler_benchmark_points.py's dodged_points, which already draws
the mean marker and omits any confidence interval), read off the same log
scale and the same y-limits as candidate A1 for a fair comparison.

Generate both A1 and A2 and compare; this script does not pick a winner.
"""

from __future__ import annotations

from pathlib import Path
import sys

PROJECT = next(
    parent for parent in Path(__file__).resolve().parents
    if (parent / "assembly_analysis").is_dir()
)
sys.path.insert(0, str(PROJECT / "assembly_analysis" / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ModuleNotFoundError as error:
    raise SystemExit(
        f"A required plotting package is missing: {error.name}\n"
        f"Python interpreter: {sys.executable}\n"
    ) from error

import assembler_benchmark_bars as bars
import assembler_benchmark_points as points
from assembler_benchmark_n50_logbars import N50_LOG_YMAX, N50_LOG_YMIN
from utils.benchmark_data import FINAL_METRICS, compute_y_limits, print_source_audit
from utils.plot_style import (
    FULL_WIDTH_IN,
    apply_style,
    clean_spines,
    draw_configuration_axis,
    panel_legend,
    panel_letter,
    sample_marker_handles,
    mean_marker_handle,
    save_figure,
    style_log_y_axis,
    subtle_grid,
)

OUTPUT_PNG = (
    PROJECT / "assembly_analysis" / "figures" / "30x" / "final" / "assembler_benchmark_n50_points.png"
)
OUTPUT_PDF = OUTPUT_PNG.with_suffix(".pdf")

N50_METRIC, N50_LABEL, N50_LETTER = FINAL_METRICS[0]
PANELS_B_TO_F = FINAL_METRICS[1:]


def main() -> int:
    data = bars.build_plot_data()
    bars.verify_plot_data_matches_source(data)
    print_source_audit(data, FINAL_METRICS)
    bars.write_plot_data(data)
    y_limits = compute_y_limits(data, PANELS_B_TO_F)

    apply_style()
    figure, axes = plt.subplots(nrows=2, ncols=3, figsize=(FULL_WIDTH_IN, 6.2))

    axis_a = axes.flat[0]
    points.dodged_points(axis_a, data, N50_METRIC)
    style_log_y_axis(axis_a, "N50 (Mb, log scale)", N50_LOG_YMIN, N50_LOG_YMAX)
    axis_a.set_xlim(bars.X_MIN, bars.X_MAX)
    axis_a.set_facecolor("white")
    clean_spines(axis_a)
    subtle_grid(axis_a, "y")
    draw_configuration_axis(axis_a)
    panel_letter(axis_a, N50_LETTER)

    for axis, (metric, label, letter) in zip(axes.flat[1:], PANELS_B_TO_F):
        bars.grouped_bars(axis, data, metric)
        bars.style_axis(axis, label, y_limits[metric], metric)
        panel_letter(axis, letter)

    panel_legend(
        figure, sample_marker_handles() + [mean_marker_handle()], y=0.965, ncols=4,
    )

    figure.patch.set_facecolor("white")
    figure.subplots_adjust(left=0.06, right=0.99, bottom=0.08, top=0.92, hspace=0.55, wspace=0.30)

    save_figure(figure, OUTPUT_PNG, OUTPUT_PDF)

    print(f"\nPNG: {OUTPUT_PNG}")
    print(f"PDF: {OUTPUT_PDF}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
