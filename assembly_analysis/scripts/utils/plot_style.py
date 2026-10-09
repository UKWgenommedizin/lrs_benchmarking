"""Shared publication-style helpers for the 30x final assembly figures.

Scoped to assembly_analysis so it can evolve independently of
alignment_analysis/scripts/utils/plot_style.py (see assembly_analysis/README.md,
"It does not run the assemblers themselves" / the assessment-vs-analysis
split). SAMPLE_ORDER and SAMPLE_COLORS are intentionally identical values to
that module's -- same GIAB samples, same meaning, same report -- so a reader
sees one consistent color mapping across every figure in the thesis.
"""

from __future__ import annotations

import matplotlib.font_manager as font_manager
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.patches import Patch
from pathlib import Path

# ============================================================
# SAMPLE / ASSEMBLER IDENTITY (fixed across every figure)
# ============================================================

SAMPLE_ORDER = ["HG002", "HG003", "HG004"]

SAMPLE_COLORS = {
    "HG002": "#0072B2",
    "HG003": "#D55E00",
    "HG004": "#009E73",
}

ASSEMBLER_ORDER = ["flye", "goldrush", "verkko"]
ASSEMBLER_TITLES = {"flye": "Flye", "goldrush": "GoldRush", "verkko": "Verkko"}

# One assembly configuration per (assembler, technology) combination that
# actually exists in assembly_benchmark_30x.tsv. Verkko is hybrid (single
# assembly per sample, not split by technology) -- see
# assembly_analysis/README.md, "Scientific comparison rules" #2/#3.
CONFIGURATION_ORDER = [
    ("flye", "ont"),
    ("flye", "pb"),
    ("goldrush", "ont"),
    ("goldrush", "pb"),
    ("verkko", "hybrid"),
]

CONFIGURATION_LABELS = {
    ("flye", "ont"): "Flye\nONT",
    ("flye", "pb"): "Flye\nHiFi",
    ("goldrush", "ont"): "Goldrush\nONT",
    ("goldrush", "pb"): "Goldrush\nHiFi",
    ("verkko", "hybrid"): "Verkko\nHybrid",
}

# Short, single-line technology label for the innermost x-axis tier of the
# 3-tier benchmark figures (assembler_benchmark_bars.py / _points.py):
# tier 1 = SINGLE-TECHNOLOGY / HYBRID, tier 2 = assembler name, tier 3 = this.
TECHNOLOGY_TIER_LABELS = {
    ("flye", "ont"): "ONT",
    ("flye", "pb"): "HiFi",
    ("goldrush", "ont"): "ONT",
    ("goldrush", "pb"): "HiFi",
    ("verkko", "hybrid"): "ONT + HiFi",
}

# x position per configuration: evenly spaced (1 unit apart) within the
# single-technology block, then an extra gap before Verkko so the hybrid
# column is visibly set apart by spacing alone, not a background wash or a
# bracket glyph (see benchmark_data.py / assembler_benchmark_bars.py).
VERKKO_GAP = 1.6
CONFIGURATION_POSITIONS = {
    ("flye", "ont"): 0.0,
    ("flye", "pb"): 1.0,
    ("goldrush", "ont"): 2.0,
    ("goldrush", "pb"): 3.0,
    ("verkko", "hybrid"): 3.0 + VERKKO_GAP,
}
VERKKO_SEPARATOR_X = 3.0 + VERKKO_GAP / 2

# Tier-2 (assembler) and tier-1 (input-strategy) spans, as (center_x, x_left,
# x_right) over the configurations each groups.
ASSEMBLER_SPANS = {
    "flye": (0.5, -0.5, 1.5),
    "goldrush": (2.5, 1.5, 3.5),
    "verkko": (3.0 + VERKKO_GAP, 3.0 + VERKKO_GAP - 0.5, 3.0 + VERKKO_GAP + 0.5),
}
STRATEGY_SPANS = {
    "single-technology": (1.5, -0.5, 3.5),
    "hybrid": (3.0 + VERKKO_GAP, 3.0 + VERKKO_GAP - 0.5, 3.0 + VERKKO_GAP + 0.5),
}

NEUTRAL_GRAY = "#8C8C8C"
GRID_COLOR = "#E5E5E5"

# A light neutral wash behind a subset of x-axis categories that differ in
# kind from the rest (e.g. Verkko's combined-input, haplotype-resolved
# assembly vs. the single-technology Flye/Goldrush runs). Deliberately not a
# bracket-with-end-ticks: that glyph is the established convention for a
# statistical-significance comparison between bars in Nature/Cell/Science
# figures, and reusing it for a grouping label risks being misread as one.
# Kept here for scripts that still use it; assembler_benchmark_bars.py /
# _points.py use spacing + a thin separator line instead (per request).
HYBRID_SPAN_COLOR = NEUTRAL_GRAY
HYBRID_SPAN_ALPHA = 0.15


# ============================================================
# FONT DETECTION
#
# Arial -> Helvetica -> Liberation Sans -> DejaVu Sans. Resolved once,
# programmatically, so an unavailable family in the chain never reaches
# matplotlib's findfont and never emits a warning.
# ============================================================

_FONT_PREFERENCE = ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"]


def detect_font_family() -> str:
    available = {font.name for font in font_manager.fontManager.ttflist}

    higher_priority = set(_FONT_PREFERENCE) - {"DejaVu Sans"}
    if not available & higher_priority:
        for font_path in font_manager.findSystemFonts():
            lower_path = font_path.lower()
            if "italic" in lower_path or "oblique" in lower_path:
                continue
            try:
                name = font_manager.get_font(font_path).family_name
            except Exception:
                continue
            if name in _FONT_PREFERENCE:
                font_manager.fontManager.addfont(font_path)
                available.add(name)

    for candidate in _FONT_PREFERENCE:
        if candidate in available:
            return candidate
    return "DejaVu Sans"


FONT_FAMILY = detect_font_family()


# ============================================================
# SIZE TARGETS (matches alignment_analysis's final-figure sizing)
# ============================================================

FULL_WIDTH_IN = 7.16
HALF_WIDTH_IN = 3.2

BASE_FONT_SIZE = 7.0
AXIS_LABEL_SIZE = 7.5
TITLE_SIZE = 8.0
TICK_LABEL_SIZE = 6.5
LEGEND_SIZE = 9.0

PANEL_TICK_SIZE = 6.0
PANEL_LABEL_SIZE_PT = 6.5
PANEL_TITLE_SIZE = 6.5
PANEL_LEGEND_SIZE = 6.0


def apply_style() -> None:
    """Set the shared rcParams. Call once near the top of each script."""
    plt.rcParams.update(
        {
            "font.family": FONT_FAMILY,
            "font.size": BASE_FONT_SIZE,
            "axes.labelsize": AXIS_LABEL_SIZE,
            "axes.titlesize": TITLE_SIZE,
            "xtick.labelsize": TICK_LABEL_SIZE,
            "ytick.labelsize": TICK_LABEL_SIZE,
            "legend.fontsize": LEGEND_SIZE,
            "axes.linewidth": 0.7,
            "xtick.major.width": 0.7,
            "ytick.major.width": 0.7,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "axes.unicode_minus": True,
        }
    )


# ============================================================
# AXIS / SPINE / GRID
# ============================================================


def clean_spines(axis: plt.Axes, keep=("left", "bottom")) -> None:
    for spine_name, spine in axis.spines.items():
        spine.set_visible(spine_name in keep)
    for spine_name in keep:
        axis.spines[spine_name].set_color("black")
        axis.spines[spine_name].set_linewidth(0.7)
        axis.spines[spine_name].set_zorder(5)


def subtle_grid(axis: plt.Axes, axis_direction: str = "y") -> None:
    axis.grid(
        axis=axis_direction,
        color=GRID_COLOR,
        linewidth=0.45,
        linestyle="-",
        zorder=0,
    )
    axis.set_axisbelow(True)


def rotated_xticks(
    axis: plt.Axes,
    positions,
    labels,
    fontsize: float | None = None,
    rotation: float = 0,
    tick_length: float | None = None,
) -> None:
    """Tick labels for the x-axis. Default rotation=0 fits the two-line
    "Assembler\\nTechnology" configuration labels without needing to slant
    them; pass a nonzero rotation for single-line labels that need it.
    """
    axis.set_xticks(positions, labels)
    for label in axis.get_xticklabels():
        if rotation:
            label.set_rotation(rotation)
            label.set_rotation_mode("anchor")
            label.set_horizontalalignment("right")
        else:
            label.set_horizontalalignment("center")
        if fontsize is not None:
            label.set_fontsize(fontsize)
    axis.tick_params(axis="x", pad=3)
    if tick_length is not None:
        axis.tick_params(axis="x", length=tick_length)


def panel_legend(figure: plt.Figure, handles, y: float, ncols: int | None = None) -> None:
    """Compact, frameless figure-level legend row centered at height y."""
    figure.legend(
        handles=handles,
        frameon=False,
        ncols=ncols or len(handles),
        loc="lower center",
        bbox_to_anchor=(0.5, y),
        fontsize=PANEL_LEGEND_SIZE,
        handlelength=0.9,
        handleheight=0.9,
        handletextpad=0.35,
        columnspacing=0.9,
        borderaxespad=0,
        borderpad=0,
    )


def sample_legend_handles(colors: dict[str, str] | None = None) -> list[Patch]:
    colors = colors or SAMPLE_COLORS
    return [
        Patch(facecolor=colors[sample], edgecolor="none", label=sample)
        for sample in SAMPLE_ORDER
    ]


def hybrid_span_handle(label: str = "Combined ONT + HiFi input") -> Patch:
    return Patch(
        facecolor=HYBRID_SPAN_COLOR, alpha=HYBRID_SPAN_ALPHA, edgecolor="none", label=label
    )


def sample_marker_handles(colors: dict[str, str] | None = None):
    from matplotlib.lines import Line2D

    colors = colors or SAMPLE_COLORS
    return [
        Line2D(
            [0], [0], marker="o", linestyle="none",
            markerfacecolor=colors[sample], markeredgecolor="black", markeredgewidth=0.4,
            markersize=5, label=sample,
        )
        for sample in SAMPLE_ORDER
    ]


def mean_marker_handle(label: str = "Mean of 3 samples"):
    from matplotlib.lines import Line2D

    return Line2D([0], [0], color="black", linewidth=1.4, label=label)


def panel_letter(
    axis: plt.Axes, letter: str, x: float = -0.14, y: float = 1.08, fontsize: float | None = None,
) -> None:
    axis.text(
        x, y, letter, transform=axis.transAxes,
        fontsize=fontsize or (PANEL_TITLE_SIZE + 1.5), fontweight="bold",
        ha="left", va="bottom",
    )


def vertical_separator(axis: plt.Axes, x: float, y_bottom: float = -0.62, y_top: float = 1.0) -> None:
    """A single thin rule marking a category boundary -- not a bracket.

    y in axes-fraction so it clears whatever multi-tier x-axis text sits
    below the plot box, regardless of each panel's own y-data range.
    """
    from matplotlib.transforms import blended_transform_factory

    transform = blended_transform_factory(axis.transData, axis.transAxes)
    axis.plot(
        [x, x], [y_bottom, y_top], transform=transform,
        color=GRID_COLOR, linewidth=0.8, linestyle="-", clip_on=False, zorder=0,
    )


def style_log_y_axis(axis: plt.Axes, y_label: str, y_min: float, y_max: float) -> None:
    """Log-scale y-axis for the N50 panel-a candidates (N50 spans ~2.4
    orders of magnitude across configurations -- see
    assembler_benchmark_n50_logbars.py / _points.py). A continuous log
    transform, not a broken axis: major ticks at powers of ten, unlabeled
    minor ticks at the 2x/5x subdivisions within each decade, both from
    matplotlib's own LogLocator rather than a hand-picked tick set.
    """
    axis.set_yscale("log")
    axis.set_ylim(y_min, y_max)
    axis.yaxis.set_major_locator(mticker.LogLocator(base=10))
    axis.yaxis.set_minor_locator(mticker.LogLocator(base=10, subs=(2, 5)))
    axis.yaxis.set_major_formatter(mticker.FuncFormatter(lambda value, _: f"{value:g}"))
    axis.yaxis.set_minor_formatter(mticker.NullFormatter())
    axis.tick_params(axis="y", which="major", labelsize=PANEL_TICK_SIZE, length=2, pad=1.5)
    axis.tick_params(axis="y", which="minor", length=1.2)
    axis.set_ylabel(y_label, fontsize=PANEL_LABEL_SIZE_PT, labelpad=2)


def draw_configuration_axis(axis: plt.Axes, show_strategy: bool = True) -> None:
    """The shared 3-tier x-axis for assembler_benchmark_bars.py / _points.py.

    Tier 1 (outermost): SINGLE-TECHNOLOGY / HYBRID input strategy.
    Tier 2: assembler name (Verkko carries the dagger here, once).
    Tier 3 (the actual matplotlib ticks): ONT / HiFi / ONT + HiFi.

    Verkko is set apart by x-spacing (VERKKO_GAP) and one thin separator
    line, not a background wash or a bracket-with-end-ticks -- called out
    verbatim in both scripts so a change here can't make them drift apart.

    show_strategy=False drops tier 1 (used by the QUAST single panels, where
    the Verkko separator and the report caption already carry that split).
    """
    from matplotlib.transforms import blended_transform_factory

    positions = [CONFIGURATION_POSITIONS[config] for config in CONFIGURATION_ORDER]
    labels = [TECHNOLOGY_TIER_LABELS[config] for config in CONFIGURATION_ORDER]
    rotated_xticks(axis, positions, labels, fontsize=PANEL_TICK_SIZE, tick_length=2)

    transform = blended_transform_factory(axis.transData, axis.transAxes)
    tier2_y = -0.16
    tier1_y = -0.27

    for assembler, (center, _, _) in ASSEMBLER_SPANS.items():
        label = ASSEMBLER_TITLES[assembler] + ("†" if assembler == "verkko" else "")
        axis.text(
            center, tier2_y, label, transform=transform, ha="center", va="top",
            fontsize=PANEL_TICK_SIZE, color="black", clip_on=False,
        )

    if not show_strategy:
        vertical_separator(axis, VERKKO_SEPARATOR_X, y_bottom=tier2_y - 0.08, y_top=1.0)
        return

    strategy_titles = {"single-technology": "SINGLE-TECHNOLOGY", "hybrid": "HYBRID"}
    for strategy, (center, _, _) in STRATEGY_SPANS.items():
        axis.text(
            center, tier1_y, strategy_titles[strategy], transform=transform, ha="center", va="top",
            fontsize=PANEL_TICK_SIZE - 0.5, color=NEUTRAL_GRAY, clip_on=False,
        )

    vertical_separator(axis, VERKKO_SEPARATOR_X, y_bottom=tier1_y - 0.08, y_top=1.0)


# ============================================================
# SAVE
# ============================================================


def save_figure(figure: plt.Figure, output_png: Path, output_pdf: Path) -> None:
    output_png.parent.mkdir(parents=True, exist_ok=True)
    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_png, dpi=300, bbox_inches="tight", facecolor="white")
    figure.savefig(output_pdf, bbox_inches="tight", facecolor="white")
    plt.close(figure)
