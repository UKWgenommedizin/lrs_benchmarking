#!/usr/bin/env python3
"""Final main assembler-benchmark figure -- grouped bars, 6 panels.

Adopted after an explicit A/B round against assembler_benchmark_points.py
(individual points): each bar is an actual HG002/HG003/HG004 observation,
never an average, which makes bars defensible on their own and consistent
with the rest of the report, so bars were kept rather than mixing geometries
across panels.

Panels (a-f): N50 (Mb), contig count, genome fraction (%), misassemblies,
mismatches / 100 kbp, indels / 100 kbp. NGA50 and runtime are both excluded,
not silently substituted:
- NGA50 (used in the earlier comparison round) spans ~65 kb-27 Mb across
  configurations -- on a shared, honest, zero-baseline linear scale most
  bars render as near-invisible slivers next to Flye-ONT, so it does not
  read honestly as a bar panel.
- Runtime is not measured at all for the actual 30x run -- see
  utils/benchmark_data.py's module docstring for what was checked.
N50 and contig count come from the same assembly_benchmark_30x.tsv as the
reference-based metrics, but are QUAST's FASTA-only "general" stats
(computed before QUAST ever aligns to the reference), not reference-based
values -- also documented in utils/benchmark_data.py.

Each bar is one individual GIAB sample (HG002/HG003/HG004), never a mean.
The Verkko caveat is not drawn inside the plot area (only "Verkko†" on
the tier-2 label); the caption text belongs in the report:

    † Verkko was assembled using combined 30× ONT and 30× PacBio HiFi
    input and is shown separately from the single-technology Flye and
    GoldRush assemblies. Reference-based statistics should therefore be
    interpreted in the context of its hybrid, haplotype-resolved assembly
    strategy.
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
    import matplotlib.ticker as mticker
except ModuleNotFoundError as error:
    raise SystemExit(
        f"A required plotting package is missing: {error.name}\n"
        f"Python interpreter: {sys.executable}\n"
    ) from error

from utils.benchmark_data import (
    FINAL_METRICS,
    build_plot_data,
    compute_y_limits,
    print_source_audit,
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
    panel_legend,
    panel_letter,
    sample_legend_handles,
    save_figure,
    subtle_grid,
)

OUTPUT_PNG = PROJECT / "assembly_analysis" / "figures" / "30x" / "final" / "assembler_benchmark_bars.png"
OUTPUT_PDF = OUTPUT_PNG.with_suffix(".pdf")

SAMPLE_OFFSETS = {"HG002": -0.27, "HG003": 0.0, "HG004": 0.27}
BAR_WIDTH = 0.26

X_MIN = -0.5
X_MAX = max(CONFIGURATION_POSITIONS.values()) + 0.5

# Panels whose counts run into the thousands get a comma thousands
# separator on the y-axis (10,000 / 20,000 / ...) instead of bare digits.
THOUSANDS_SEPARATOR_METRICS = {"contigs", "misassemblies"}


def grouped_bars(axis: plt.Axes, data, metric: str) -> None:
    for assembler, technology in CONFIGURATION_ORDER:
        x = CONFIGURATION_POSITIONS[(assembler, technology)]
        subset = data[(data["assembler"] == assembler) & (data["technology"] == technology)]
        values = subset.set_index("sample")[metric].reindex(SAMPLE_ORDER)

        for sample in SAMPLE_ORDER:
            axis.bar(
                x + SAMPLE_OFFSETS[sample], values.loc[sample], width=BAR_WIDTH,
                color=SAMPLE_COLORS[sample], edgecolor="white", linewidth=0.3, zorder=3,
            )


def style_axis(axis: plt.Axes, y_label: str, y_max: float, metric: str,
               show_strategy: bool = True) -> None:
    axis.set_ylim(0, y_max)
    axis.set_xlim(X_MIN, X_MAX)
    axis.tick_params(axis="y", labelsize=PANEL_TICK_SIZE, length=2, pad=1.5)
    axis.set_ylabel(y_label, fontsize=PANEL_LABEL_SIZE_PT, labelpad=2)
    if metric in THOUSANDS_SEPARATOR_METRICS:
        axis.yaxis.set_major_formatter(mticker.StrMethodFormatter("{x:,.0f}"))
    clean_spines(axis)
    subtle_grid(axis, "y")
    axis.set_facecolor("white")
    draw_configuration_axis(axis, show_strategy=show_strategy)


def main() -> int:
    data = build_plot_data()
    verify_plot_data_matches_source(data)
    print_source_audit(data, FINAL_METRICS)
    write_plot_data(data)
    y_limits = compute_y_limits(data, FINAL_METRICS)

    apply_style()
    figure, axes = plt.subplots(nrows=2, ncols=3, figsize=(FULL_WIDTH_IN, 6.2))

    for axis, (metric, label, letter) in zip(axes.flat, FINAL_METRICS):
        grouped_bars(axis, data, metric)
        style_axis(axis, label, y_limits[metric], metric)
        panel_letter(axis, letter)

    panel_legend(figure, sample_legend_handles(), y=0.965)

    figure.patch.set_facecolor("white")
    figure.subplots_adjust(left=0.06, right=0.99, bottom=0.08, top=0.92, hspace=0.55, wspace=0.30)

    save_figure(figure, OUTPUT_PNG, OUTPUT_PDF)

    print(f"\nRows: {len(data)}")
    print(f"PNG: {OUTPUT_PNG}")
    print(f"PDF: {OUTPUT_PDF}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
