#!/usr/bin/env python3
"""Shared single-panel renderer backing the six standalone QUAST panel
scripts (assembler_quast_nga50.py, _genome_fraction.py,
_duplication_ratio.py, _misassemblies.py, _mismatches.py, _indels.py).

Each standalone script is its own committed .py file with its own final
PNG/PDF -- so the six panels can be recombined and captioned outside this
repo -- but all six draw through this one render_panel() so the bar style,
sample colors and QUAST cross-check are identical across all six and
cannot drift apart: all import assembler_benchmark_bars.py's grouped_bars
/ style_axis directly (the same grouped-bar drawing already used for
assembler_benchmark_bars.py's own FASTA-stats panels and reused again by
assembler_quast_quality_panel.py), and all load values through
utils.benchmark_data.load_quast_values(), which re-parses each row's own
raw QUAST report.tsv and compares it to assembly_benchmark_30x.tsv. Each
bar is one individual GIAB sample observation (HG002/HG003/HG004), never
an average -- there is no mean marker to draw, unlike the earlier
point-based design.

Panel letters are intentionally NOT drawn here. This follows the
convention already established by
alignment_analysis/scripts/30x/plots/final/plot_figure3_panels_30x.py:
standalone panels are undecorated, and letters / grid layout are added
when the panels are composed into a figure elsewhere (e.g. in LaTeX).
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
    from matplotlib.text import Text
except ModuleNotFoundError as error:
    raise SystemExit(
        f"A required plotting package is missing: {error.name}\n"
        f"Python interpreter: {sys.executable}\n"
    ) from error

import assembler_benchmark_bars as bars_base
from utils.benchmark_data import QUAST_FINAL_METRICS, load_quast_values, print_source_audit
from utils.plot_style import (
    HALF_WIDTH_IN,
    apply_style,
    panel_legend,
    sample_legend_handles,
    save_figure,
)

PANEL_FIGSIZE = (HALF_WIDTH_IN, 3.1)

# These panels are placed in the report at 0.32\linewidth (three across in
# Figure 4), not the ~0.49\linewidth (roughly native HALF_WIDTH_IN size,
# ~1:1 scale) used for the other report figures' panels -- so at the
# PANEL_* sizes from utils.plot_style, Figure 4's text renders visibly
# smaller on the page than Figures 2/3's, even though the underlying point
# sizes are identical. TEXT_SCALE pre-enlarges every text element here so
# that after LaTeX's extra ~0.63x shrink (0.32 / 0.49), the effective
# on-page size matches the rest of the report again.
TEXT_SCALE = 1.6

# plotted column -> the assembly_benchmark_30x.tsv column it is loaded from.
# Only nga50_mb differs from its source column: it is a unit conversion
# (bp / 1e6) of the table's nga50 column for axis readability, not a rescale.
RAW_COLUMN_FOR = {
    "nga50_mb": "nga50",
    "genome_fraction_pct": "genome_fraction_pct",
    "duplication_ratio": "duplication_ratio",
    "misassemblies": "misassemblies",
    "mismatches_per_100kbp": "mismatches_per_100kbp",
    "indels_per_100kbp": "indels_per_100kbp",
}

# Reuse the compound figure's own labels so the two can't drift apart.
LABEL_FOR = {column: label for column, label, _ in QUAST_FINAL_METRICS}


def render_panel(plot_column: str, output_png: Path):
    """Load, cross-check, print the source audit for, and save one panel.

    Returns the loaded DataFrame so the caller's main() can report row
    counts alongside the PNG/PDF paths.
    """
    raw_column = RAW_COLUMN_FOR[plot_column]
    label = LABEL_FOR[plot_column]

    data = load_quast_values([raw_column])  # raises on any table-vs-report.tsv mismatch
    if plot_column != raw_column:
        data[plot_column] = data[raw_column] / 1_000_000

    print_source_audit(data, [(plot_column, label, "")])

    y_max = float(data[plot_column].max()) * 1.12

    apply_style()
    figure, axis = plt.subplots(figsize=PANEL_FIGSIZE)
    bars_base.grouped_bars(axis, data, plot_column)
    bars_base.style_axis(axis, label, y_max, plot_column)
    axis.grid(False)  # drop style_axis's default gray horizontal gridlines

    panel_legend(
        figure, sample_legend_handles(), y=0.90, ncols=3,
    )
    figure.patch.set_facecolor("white")
    figure.subplots_adjust(left=0.22, right=0.97, bottom=0.30, top=0.80)

    for text in figure.findobj(match=lambda artist: isinstance(artist, Text)):
        text.set_fontsize(text.get_fontsize() * TEXT_SCALE)

    output_pdf = output_png.with_suffix(".pdf")
    save_figure(figure, output_png, output_pdf)
    return data
