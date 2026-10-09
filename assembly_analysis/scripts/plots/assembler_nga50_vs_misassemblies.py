#!/usr/bin/env python3
"""Main Figure A -- contiguity vs structural correctness (QUAST only).

One point per HG002/HG003/HG004 assembly: x = QUAST misassemblies,
y = QUAST NGA50 (Mb, log scale; the values span ~65 kb to ~27 Mb).
Color = sample (fixed thesis colors); marker shape = assembler; filled =
ONT input, open = HiFi input, diamond = Verkko (ONT + HiFi hybrid). Each
configuration's three samples cluster, so every cluster is direct-labeled
instead of relying on a marker legend.

Values come from assembly_benchmark_30x.tsv, cross-checked against each raw
QUAST report.tsv (utils.benchmark_data.load_quast_values).

Caption notes (belong in the report, not in the figure):

    † Verkko was assembled using combined 30× ONT and 30× PacBio HiFi input
    and produces haplotype-resolved assemblies. It is therefore shown
    separately from the single-technology Flye and GoldRush runs.
    Reference-based statistics should be interpreted in the context of this
    different assembly representation.

    NGA50 is computed against the haploid GRCh38 length. For assemblies
    with duplication ratio > 1 (Verkko ~2.05, Flye HiFi ~1.33) redundant
    aligned blocks count toward that length, and misassemblies are absolute
    counts over more sequence, so both axes favour / penalise those
    assemblies relative to a collapsed haploid assembly.
"""

from __future__ import annotations

from pathlib import Path
import sys

PROJECT = next(
    parent for parent in Path(__file__).resolve().parents
    if (parent / "assembly_analysis").is_dir()
)
sys.path.insert(0, str(PROJECT / "assembly_analysis" / "scripts"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D

from utils.benchmark_data import load_quast_values, print_source_audit
from utils.plot_style import (
    CONFIGURATION_ORDER,
    NEUTRAL_GRAY,
    PANEL_LABEL_SIZE_PT,
    PANEL_LEGEND_SIZE,
    PANEL_TICK_SIZE,
    SAMPLE_COLORS,
    SAMPLE_ORDER,
    apply_style,
    clean_spines,
    save_figure,
    style_log_y_axis,
)

OUTPUT_PNG = PROJECT / "assembly_analysis" / "figures" / "30x" / "final" / "assembler_nga50_vs_misassemblies.png"
OUTPUT_PDF = OUTPUT_PNG.with_suffix(".pdf")

METRICS = [
    ("nga50", "NGA50 (bp)", ""),
    ("misassemblies", "Misassemblies", ""),
]

MARKERS = {"flye": "o", "goldrush": "s", "verkko": "D"}
FILLED = {"ont": True, "pb": False, "hybrid": True}
CLUSTER_LABELS = {
    ("flye", "ont"): "Flye ONT",
    ("flye", "pb"): "Flye HiFi",
    ("goldrush", "ont"): "GoldRush ONT",
    ("goldrush", "pb"): "GoldRush HiFi",
    ("verkko", "hybrid"): "Verkko† ONT + HiFi",
}
# Label placement relative to the cluster (offset in points, alignment):
# right of the cluster by default; Verkko sits at the right edge, so above.
LABEL_PLACEMENT = {
    ("verkko", "hybrid"): ((0, 9), "center", "bottom"),
}
DEFAULT_PLACEMENT = ((9, 0), "left", "center")

Y_MIN_MB, Y_MAX_MB = 0.05, 50


def main() -> int:
    data = load_quast_values([column for column, _, _ in METRICS])
    print_source_audit(data, METRICS)
    data["nga50_mb"] = data["nga50"] / 1_000_000

    apply_style()
    figure, axis = plt.subplots(figsize=(4.4, 3.4))

    for assembler, technology in CONFIGURATION_ORDER:
        subset = data[(data["assembler"] == assembler) & (data["technology"] == technology)]
        for sample in SAMPLE_ORDER:
            row = subset[subset["sample"] == sample].iloc[0]
            color = SAMPLE_COLORS[sample]
            axis.scatter(
                row["misassemblies"], row["nga50_mb"],
                marker=MARKERS[assembler], s=34,
                facecolors=color if FILLED[technology] else "white",
                edgecolors=color if not FILLED[technology] else "white",
                linewidths=1.2 if not FILLED[technology] else 0.6,
                zorder=3,
            )

        # geometric-mean y (log axis), rightmost x of the three samples
        y_center = (subset["nga50_mb"].prod()) ** (1 / len(subset))
        offset, ha, va = LABEL_PLACEMENT.get((assembler, technology), DEFAULT_PLACEMENT)
        if (assembler, technology) in LABEL_PLACEMENT:
            anchor = (subset["misassemblies"].mean(), subset["nga50_mb"].max())
        else:
            anchor = (subset["misassemblies"].max(), y_center)
        axis.annotate(
            CLUSTER_LABELS[(assembler, technology)], anchor, xytext=offset, textcoords="offset points",
            ha=ha, va=va, fontsize=PANEL_TICK_SIZE,
            color=NEUTRAL_GRAY if assembler == "verkko" else "black",
        )

    style_log_y_axis(axis, "NGA50 (Mb)", Y_MIN_MB, Y_MAX_MB)
    axis.set_xlim(0, 21_000)
    axis.xaxis.set_major_formatter(mticker.StrMethodFormatter("{x:,.0f}"))
    axis.tick_params(axis="x", labelsize=PANEL_TICK_SIZE, length=2, pad=1.5)
    axis.set_xlabel("Misassemblies", fontsize=PANEL_LABEL_SIZE_PT, labelpad=2)
    clean_spines(axis)
    axis.set_facecolor("white")

    handles = [
        Line2D([0], [0], marker="o", linestyle="none", markersize=5,
               markerfacecolor=SAMPLE_COLORS[sample], markeredgecolor="white", label=sample)
        for sample in SAMPLE_ORDER
    ] + [
        Line2D([0], [0], marker="o", linestyle="none", markersize=5,
               markerfacecolor=NEUTRAL_GRAY, markeredgecolor="white", label="ONT"),
        Line2D([0], [0], marker="o", linestyle="none", markersize=5, markeredgewidth=1.1,
               markerfacecolor="white", markeredgecolor=NEUTRAL_GRAY, label="HiFi"),
    ]
    figure.legend(
        handles=handles, frameon=False, ncols=len(handles), loc="lower center",
        bbox_to_anchor=(0.5, 0.95), fontsize=PANEL_LEGEND_SIZE,
        handletextpad=0.2, columnspacing=0.9, borderaxespad=0,
    )

    figure.patch.set_facecolor("white")
    figure.subplots_adjust(left=0.13, right=0.97, bottom=0.13, top=0.92)
    save_figure(figure, OUTPUT_PNG, OUTPUT_PDF)

    print(f"\nPNG: {OUTPUT_PNG}\nPDF: {OUTPUT_PDF}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
