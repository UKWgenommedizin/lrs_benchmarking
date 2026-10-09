#!/usr/bin/env python3
"""Candidate A1: assembler_benchmark_bars.py with panel (a) on a log y-axis.

Panels b-f are identical to assembler_benchmark_bars.py (reused directly,
not re-implemented, so they cannot drift from the adopted figure). Only
panel (a) changes: N50 spans ~0.14-33.6 Mb across configurations (~2.4
orders of magnitude), so a linear zero-baseline axis compresses everything
but Flye-ONT to a near-invisible sliver (see assembler_benchmark_bars.py's
docstring). This version keeps grouped bars -- individual HG002/HG003/HG004
bars, not means -- but reads them off a log scale instead.

Companion: assembler_benchmark_n50_points.py (candidate A2, same log scale,
individual points instead of bars for panel a). Generate both and compare
before choosing; this script does not pick a winner.
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
from utils.benchmark_data import FINAL_METRICS, compute_y_limits, print_source_audit
from utils.plot_style import (
    FULL_WIDTH_IN,
    apply_style,
    clean_spines,
    draw_configuration_axis,
    panel_legend,
    panel_letter,
    sample_legend_handles,
    save_figure,
    style_log_y_axis,
    subtle_grid,
)

OUTPUT_PNG = (
    PROJECT / "assembly_analysis" / "figures" / "30x" / "final" / "assembler_benchmark_n50_logbars.png"
)
OUTPUT_PDF = OUTPUT_PNG.with_suffix(".pdf")

N50_METRIC, N50_LABEL, N50_LETTER = FINAL_METRICS[0]
PANELS_B_TO_F = FINAL_METRICS[1:]

# Round decades bracketing the actual N50 range (0.1377-33.6 Mb), not a
# hand-tuned fit -- see the module docstring.
N50_LOG_YMIN = 0.1
N50_LOG_YMAX = 50.0


def main() -> int:
    data = bars.build_plot_data()
    bars.verify_plot_data_matches_source(data)
    print_source_audit(data, FINAL_METRICS)
    bars.write_plot_data(data)
    y_limits = compute_y_limits(data, PANELS_B_TO_F)

    apply_style()
    figure, axes = plt.subplots(nrows=2, ncols=3, figsize=(FULL_WIDTH_IN, 6.2))

    axis_a = axes.flat[0]
    bars.grouped_bars(axis_a, data, N50_METRIC)
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

    panel_legend(figure, sample_legend_handles(), y=0.965)

    figure.patch.set_facecolor("white")
    figure.subplots_adjust(left=0.06, right=0.99, bottom=0.08, top=0.92, hspace=0.55, wspace=0.30)

    save_figure(figure, OUTPUT_PNG, OUTPUT_PDF)

    print(f"\nPNG: {OUTPUT_PNG}")
    print(f"PDF: {OUTPUT_PDF}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
