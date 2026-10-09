#!/usr/bin/env python3
"""Draw the compound-figure panels for report Figure 5 (30x).

Figure 5 places the alignment- and assembly-based strategies side by side as
categorical dot plots: one mark per HG002/HG003/HG004 sample (color = sample,
shape = sequencing input) plus a subdued mean +/- SD of the three samples
beside them. Panels are drawn at their printed size with no internal letters;
the letters, row titles and column headers are added in LaTeX, and the legend
is drawn once as its own strip (fig5_legend_30x) placed above the panels.

Values come from plot_cross_strategy_comparison_30x.build_panel_data(), so
they are identical to the cross-strategy tables (panel_data.tsv / summary.tsv);
nothing is recalculated here.

Metric definitions (verified against the analysis code and samtools/QUAST runs):

  a  NM edit distance per 100 kbp: samtools stats `mismatches` (sum of the
     aligner's NM tags = substituted + inserted + deleted bases) divided by
     `bases mapped (cigar)` (M/I/=/X read bases, soft clips excluded) x 1e5.
     NOT a substitution-only measure.
  b  QUAST `# mismatches per 100 kbp`: substitutions only, per 100 kbp of
     aligned contig bases.
  c  CIGAR indel EVENTS per 100 kbp: sum of samtools stats ID records
     (number of insertion + deletion operations per length, 1-300 bp;
     samtools skips longer indels) divided by `bases mapped (cigar)` x 1e5.
  d  QUAST `# indels per 100 kbp`: indel events per 100 kbp of aligned contig
     bases.
  e  CIGAR-aligned yield: `bases mapped (cigar)` / largest samtools
     `total_length` (primary read bases) among the four aligners for that
     sample and technology x 100.
  f  QUAST genome fraction: % of GRCh38 bases covered by aligned contigs.

samtools stats: no -F/-f/-q filters, so secondary alignments are skipped
(samtools default), supplementary alignments and MAPQ 0 are included,
unmapped reads are excluded. QUAST 5.3.0: --large --min-contig 500 (min
alignment 500 bp, min identity 95 %, ambiguity "one").
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
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator, StrMethodFormatter

import plot_cross_strategy_comparison_30x as cross
from utils.plot_style import (
    ALIGNER_ORDER,
    SAMPLE_COLORS,
    SAMPLE_ORDER,
    apply_style,
    clean_spines,
    save_figure,
    subtle_grid,
)

OUT_DIR = PROJECT / "figures" / "30x" / "cross_strategy" / "panels"

# ---------------------------------------------------------------- layout
# Reads panels get more width than contig panels (8 vs 5 categories) so the
# spacing per category is about equal; LaTeX includes them at 0.56 and
# 0.42 \linewidth (6.5 in text block), i.e. at ~1:1 scale.
READ_PANEL_SIZE = (3.65, 2.05)
CONTIG_PANEL_SIZE = (2.75, 2.05)
READ_ADJUST = dict(left=0.155, right=0.995, bottom=0.20, top=0.86, wspace=0.07)
CONTIG_ADJUST = dict(left=0.205, right=0.99, bottom=0.20, top=0.86, wspace=0.10)

# ---------------------------------------------------------------- typography
TICK_SIZE = 7.5
LABEL_SIZE = 8.0
TITLE_SIZE = 8.0
LEGEND_SIZE = 8.0

# ---------------------------------------------------------------- categories
READ_COLUMNS = [("ONT", "ONT"), ("HiFi", "PacBio HiFi")]
CONTIG_COLUMNS = [
    ("ONT", "ONT", ["Flye", "GoldRush"]),
    ("HiFi", "PacBio HiFi", ["Flye", "GoldRush"]),
    ("ONT+HiFi", "ONT +\nHiFi", ["Verkko"]),
]
TICK_NAMES = {
    "minimap2": "mini-\nmap2",
    "VACmap": "VAC-\nmap",
    "VG Giraffe": "VG\nGiraffe",
    "GoldRush": "Gold-\nRush",
    "Verkko": "Verkko†",
}

# Shape = sequencing input; sizes balance the visual area of each shape.
INPUT_MARKER = {"ONT": "o", "HiFi": "s", "ONT+HiFi": "D"}
INPUT_LABEL = {"ONT": "ONT", "HiFi": "PacBio HiFi", "ONT+HiFi": "ONT + HiFi"}
MARKER_AREA = {"o": 15, "s": 12.5, "D": 11}

# ---------------------------------------------------------------- dots + summary
# Offsets in inches on the printed panel, so a triplet and its mean +/- SD
# look identical in every column; the group is centred on the tick.
SAMPLE_X_IN = dict(zip(SAMPLE_ORDER, [-0.085, -0.03, 0.025]))
SUMMARY_X_IN = 0.095
SUMMARY_COLOR = "#6E6E6E"
SUMMARY_STYLE = dict(
    fmt="_", color=SUMMARY_COLOR, markersize=6, markeredgewidth=1.2,
    ecolor=SUMMARY_COLOR, elinewidth=0.7, capsize=1.8, capthick=0.7,
)
MIN_CAPPED_SD_PT = 3.0

# Verkko: very light hatching marks the diploid (both-haplotype) assembly.
VERKKO_HATCH_COLOR = "#EDEDED"
VERKKO_HATCH = "////"

# (tidy panel, column, y label, output stem, thousands separator)
PANELS = [
    ("a", "left", "NM edit-distance bases\nper 100 kbp mapped", "fig5a_read_edit_distance_30x", True),
    ("a", "right", "QUAST mismatches\nper 100 kbp aligned", "fig5b_contig_mismatches_30x", False),
    ("b", "left", "CIGAR indel events\nper 100 kbp mapped", "fig5c_read_indels_30x", False),
    ("b", "right", "QUAST indel events\nper 100 kbp aligned", "fig5d_contig_indels_30x", False),
    ("c", "left", "CIGAR-aligned yield\n(% of input bases)", "fig5e_read_yield_30x", False),
    ("c", "right", "Genome fraction\n(% of GRCh38)", "fig5f_contig_genome_fraction_30x", False),
]


def output_paths(stem: str) -> tuple[Path, Path]:
    png = OUT_DIR / f"{stem}.png"
    return png, png.with_suffix(".pdf")


def summaries(data: pd.DataFrame) -> pd.DataFrame:
    return data.groupby(["method", "technology"])["value"].agg(["mean", "std", "count"])


def y_top(data: pd.DataFrame, is_percent: bool) -> tuple[float, list[float]]:
    """Linear axis from zero to just above the highest point or error bar."""
    stats = summaries(data)
    if (stats["count"] != 3).any():
        raise ValueError("Every configuration must have exactly three samples")
    highest = max(float(data["value"].max()), float((stats["mean"] + stats["std"]).max()))
    if is_percent:
        return 104.0, [0, 20, 40, 60, 80, 100]
    ticks = MaxNLocator(nbins=4, steps=[1, 2, 2.5, 5, 10]).tick_values(0, highest * 1.03)
    ticks = [t for t in ticks if t >= 0]
    while ticks[-1] < highest * 1.03:
        ticks.append(ticks[-1] + (ticks[1] - ticks[0]))
    return float(ticks[-1]), ticks


def dot_plot(axis: plt.Axes, subset: pd.DataFrame, groups: list[str], technology: str) -> None:
    figure = axis.get_figure()
    position = axis.get_position()
    units_per_inch = len(groups) / (position.width * figure.get_figwidth())
    y_low, y_high = axis.get_ylim()
    points_per_unit = position.height * figure.get_figheight() * 72 / (y_high - y_low)
    marker = INPUT_MARKER[technology]

    for index, group in enumerate(groups):
        values = subset[subset["method"] == group].set_index("sample")["value"].reindex(SAMPLE_ORDER)
        for sample in SAMPLE_ORDER:
            axis.scatter(
                index + SAMPLE_X_IN[sample] * units_per_inch, values.loc[sample],
                marker=marker, s=MARKER_AREA[marker], color=SAMPLE_COLORS[sample],
                edgecolor="black", linewidth=0.4, zorder=4,
            )
        mean, sd = float(values.mean()), float(values.std(ddof=1))
        style = dict(SUMMARY_STYLE)
        if sd * points_per_unit < MIN_CAPPED_SD_PT:
            # Caps on a whisker this short merge into a blob behind the mean tick.
            style.update(capsize=0)
        axis.errorbar(index + SUMMARY_X_IN * units_per_inch, mean, yerr=sd, zorder=3, **style)


def style_axis(axis: plt.Axes, title: str, groups: list[str], top: float, ticks: list[float]) -> None:
    axis.set_title(title, pad=3, fontsize=TITLE_SIZE)
    axis.set_ylim(0, top)
    axis.set_yticks(ticks)
    axis.set_xlim(-0.5, len(groups) - 0.5)
    axis.set_xticks(range(len(groups)), [TICK_NAMES.get(g, g) for g in groups])
    axis.tick_params(axis="x", labelsize=TICK_SIZE, length=2, pad=2)
    axis.tick_params(axis="y", labelsize=TICK_SIZE, length=2, pad=1.5)
    for label in axis.get_xticklabels():
        label.set_linespacing(0.95)
    subtle_grid(axis, "y")
    clean_spines(axis)
    axis.set_facecolor("white")


def shade_verkko(axis: plt.Axes) -> None:
    with plt.rc_context({"hatch.linewidth": 0.5, "hatch.color": VERKKO_HATCH_COLOR}):
        axis.axvspan(-0.5, 0.5, facecolor="none", edgecolor=VERKKO_HATCH_COLOR,
                     hatch=VERKKO_HATCH, linewidth=0, zorder=0.5)


def plot_panel(data: pd.DataFrame, column: str, y_label: str, stem: str, thousands: bool) -> None:
    top, ticks = y_top(data, "%" in y_label)
    if column == "left":
        size, adjust = READ_PANEL_SIZE, READ_ADJUST
        columns = [(tech, title, ALIGNER_ORDER) for tech, title in READ_COLUMNS]
    else:
        size, adjust = CONTIG_PANEL_SIZE, CONTIG_ADJUST
        columns = CONTIG_COLUMNS

    figure, axes = plt.subplots(
        nrows=1, ncols=len(columns), sharey=True, figsize=size,
        gridspec_kw={"width_ratios": [len(groups) for _, _, groups in columns]},
    )
    figure.subplots_adjust(**adjust)

    for axis, (technology, title, groups) in zip(axes, columns):
        style_axis(axis, title, groups, top, ticks)
        if technology == "ONT+HiFi":
            shade_verkko(axis)
        subset = data[data["technology"] == technology]
        dot_plot(axis, subset, groups, technology)

    axes[0].set_ylabel(y_label, fontsize=LABEL_SIZE, labelpad=3, linespacing=1.1)
    if thousands:
        axes[0].yaxis.set_major_formatter(StrMethodFormatter("{x:,.0f}"))
    figure.patch.set_facecolor("white")
    save_figure(figure, *output_paths(stem))


def plot_shared_legend() -> None:
    figure = plt.figure(figsize=(6.5, 0.28))
    hidden = figure.add_axes([0, 0, 0.01, 0.01])
    hidden.set_visible(False)
    mean_sd = hidden.errorbar([0], [0], yerr=[1], label="Mean ± SD (n = 3 samples)", **SUMMARY_STYLE)

    samples = [Patch(facecolor=SAMPLE_COLORS[s], edgecolor="black", linewidth=0.4, label=s)
               for s in SAMPLE_ORDER]
    inputs = [
        Line2D([0], [0], marker=INPUT_MARKER[t], linestyle="none", markerfacecolor="white",
               markeredgecolor="black", markeredgewidth=0.6,
               markersize=MARKER_AREA[INPUT_MARKER[t]] ** 0.5 * 1.1, label=INPUT_LABEL[t])
        for t in INPUT_MARKER
    ]
    figure.legend(
        handles=samples + inputs + [mean_sd],
        frameon=False, ncols=7, loc="center", fontsize=LEGEND_SIZE,
        handlelength=1.0, handleheight=0.8, handletextpad=0.35, columnspacing=1.1,
        borderaxespad=0, borderpad=0,
    )
    figure.patch.set_facecolor("white")
    save_figure(figure, *output_paths("fig5_legend_30x"))


def main() -> int:
    panel_data = cross.build_panel_data(cross.load_alignment(), cross.load_assembly())
    if panel_data["value"].isna().any():
        raise ValueError("Missing values in the cross-strategy panel data")

    apply_style()
    for key, column, y_label, stem, thousands in PANELS:
        data = panel_data[(panel_data["panel"] == key) & (panel_data["column"] == column)]
        plot_panel(data, column, y_label, stem, thousands)
        print(f"PDF: {output_paths(stem)[1]}  ({len(data)} values, y max {y_top(data, '%' in y_label)[0]:g})")
    plot_shared_legend()
    print(f"PDF: {output_paths('fig5_legend_30x')[1]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
