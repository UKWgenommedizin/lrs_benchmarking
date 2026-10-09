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
The same convention covers the sample legend: the panels carry none of
their own, and render_shared_legend() draws it once as its own strip
(assembler_quast_legend.py) that LaTeX places above the six panels, as
fig3_legend_30x does for Figure 3. The outermost SINGLE-TECHNOLOGY /
HYBRID x-axis tier is drawn only on the bottom-row panels (d-f), so the
grouping is labelled once per column of the 2x3 figure rather than six
times (show_strategy=True in those three scripts).
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
    sample_legend_handles,
    save_figure,
)

PANEL_FIGSIZE = (HALF_WIDTH_IN, 3.1)

# Shared legend strip: same width and font size as Figure 3's
# fig3_legend_30x, so both legends read identically on the page (LaTeX
# includes both at natural size).
LEGEND_FIGSIZE = (6.3, 0.25)
LEGEND_FONT_SIZE = 7.0

# These panels are placed in the report at 0.32\linewidth (three across in
# Figure 4), not the ~0.49\linewidth (roughly native HALF_WIDTH_IN size,
# ~1:1 scale) used for the other report figures' panels -- so at the
# PANEL_* sizes from utils.plot_style, Figure 4's text renders visibly
# smaller on the page than Figures 2/3's, even though the underlying point
# sizes are identical. TEXT_SCALE pre-enlarges every text element here so
# that after LaTeX's extra ~0.63x shrink (0.32 / 0.49), the effective
# on-page size matches the rest of the report again.
TEXT_SCALE = 1.75

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


def apply_broken_y_axis(axis, lower_max: float, upper_min: float, upper_max: float,
                        lower_ticks: list[float], upper_ticks: list[float]) -> None:
    """Show [0, lower_max] and [upper_min, upper_max] on one axis with a break.

    A piecewise-linear y scale gives the lower range LOWER_FRACTION of the
    axis height, a thin gap stands for the omitted (lower_max, upper_min)
    range, and the upper range takes the rest. Bars stay true-zero-based
    and are drawn exactly as before; a white band with break marks over the
    gap shows that bars crossing it are interrupted. Using one axis (not two
    stacked ones) keeps the shared 3-tier configuration labels positioned
    exactly as in the other five panels.
    """
    import numpy as np
    from matplotlib.lines import Line2D
    from matplotlib.patches import Rectangle

    lower_fraction, gap_fraction = 0.66, 0.08
    upper_fraction = 1.0 - lower_fraction - gap_fraction
    gap_bottom, gap_top = lower_fraction, lower_fraction + gap_fraction

    def forward(y):
        y = np.asarray(y, dtype=float)
        return np.where(
            y <= lower_max, y / lower_max * lower_fraction,
            np.where(
                y < upper_min,
                gap_bottom + (y - lower_max) / (upper_min - lower_max) * gap_fraction,
                gap_top + (y - upper_min) / (upper_max - upper_min) * upper_fraction,
            ),
        )

    def inverse(f):
        f = np.asarray(f, dtype=float)
        return np.where(
            f <= gap_bottom, f / lower_fraction * lower_max,
            np.where(
                f < gap_top,
                lower_max + (f - gap_bottom) / gap_fraction * (upper_min - lower_max),
                upper_min + (f - gap_top) / upper_fraction * (upper_max - upper_min),
            ),
        )

    axis.set_yscale("function", functions=(forward, inverse))
    axis.set_ylim(0, upper_max)
    axis.set_yticks(lower_ticks + upper_ticks)
    axis.yaxis.set_major_formatter(
        __import__("matplotlib.ticker", fromlist=["FuncFormatter"]).FuncFormatter(
            lambda value, _: f"{value:g}"
        )
    )

    # White band over the gap (also interrupts the bars that cross it),
    # extended left over the spine so the spine is visibly cut.
    axis.add_patch(Rectangle(
        (-0.03, gap_bottom), 1.03, gap_fraction, transform=axis.transAxes,
        facecolor="white", edgecolor="none", zorder=5, clip_on=False,
    ))
    for y in (gap_bottom, gap_top):
        axis.add_line(Line2D(
            [-0.025, 0.025], [y - 0.012, y + 0.012], transform=axis.transAxes,
            color="black", linewidth=0.8, zorder=6, clip_on=False,
        ))


def render_panel(plot_column: str, output_png: Path, y_break: dict | None = None,
                 show_strategy: bool = False):
    """Load, cross-check, print the source audit for, and save one panel.

    y_break: optional keyword arguments for apply_broken_y_axis(), for a
    metric whose values span orders of magnitude (NGA50).
    show_strategy: draw the SINGLE-TECHNOLOGY / HYBRID tier (bottom-row panels).

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
    bars_base.style_axis(axis, label, y_max, plot_column, show_strategy=show_strategy)
    axis.grid(False)  # drop style_axis's default gray horizontal gridlines
    if y_break is not None:
        values = data[plot_column]
        hidden = values[(values > y_break["lower_max"]) & (values < y_break["upper_min"])]
        if len(hidden) or values.max() > y_break["upper_max"]:
            raise ValueError(f"Broken axis would hide values: {sorted(hidden)} / max {values.max()}")
        apply_broken_y_axis(axis, **y_break)

    figure.patch.set_facecolor("white")
    figure.subplots_adjust(left=0.22, right=0.97, bottom=0.24, top=0.97)

    for text in figure.findobj(match=lambda artist: isinstance(artist, Text)):
        text.set_fontsize(text.get_fontsize() * TEXT_SCALE)

    output_pdf = output_png.with_suffix(".pdf")
    save_figure(figure, output_png, output_pdf)
    return data


def render_shared_legend(output_png: Path) -> None:
    """One sample-color legend strip for all six panels, placed above them in LaTeX."""
    apply_style()
    figure = plt.figure(figsize=LEGEND_FIGSIZE)
    figure.legend(
        handles=sample_legend_handles(),
        frameon=False, ncols=3, loc="center", fontsize=LEGEND_FONT_SIZE,
        handlelength=1.2, handleheight=0.9, handletextpad=0.4, columnspacing=1.4,
        borderaxespad=0, borderpad=0,
    )
    figure.patch.set_facecolor("white")
    save_figure(figure, output_png, output_png.with_suffix(".pdf"))
