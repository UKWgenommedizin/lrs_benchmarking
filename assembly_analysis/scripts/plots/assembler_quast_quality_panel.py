#!/usr/bin/env python3
"""Main Figure B -- QUAST quality metrics, grouped bars, 2x2.

Panels: (a) genome fraction (%), (b) duplication ratio, (c) mismatches /
100 kbp, (d) indels / 100 kbp. QUAST is the sole source; every bar is one
HG002/HG003/HG004 observation, never a mean. Values are read from
assembly_benchmark_30x.tsv and cross-checked against each raw QUAST
report.tsv (utils.benchmark_data.load_quast_values).

The duplication-ratio panel carries a dashed line at 1.0 (no redundant
sequence relative to the haploid reference). Verkko's ~2.05 is expected for a
FASTA holding both haplotypes; it is not an error.

Caption notes (belong in the report, not in the figure):

    † Verkko was assembled using combined 30× ONT and 30× PacBio HiFi input
    and produces haplotype-resolved assemblies. It is therefore shown
    separately from the single-technology Flye and GoldRush runs.
    Reference-based statistics should be interpreted in the context of this
    different assembly representation.

    Mismatches and indels are relative to GRCh38 and include true sample
    variation, not only assembly errors.
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

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from assembler_benchmark_bars import grouped_bars, style_axis
from utils.benchmark_data import load_quast_values, print_source_audit
from utils.plot_style import (
    FULL_WIDTH_IN,
    NEUTRAL_GRAY,
    apply_style,
    panel_legend,
    panel_letter,
    sample_legend_handles,
    save_figure,
)

OUTPUT_PNG = PROJECT / "assembly_analysis" / "figures" / "30x" / "final" / "assembler_quast_quality_panel.png"
OUTPUT_PDF = OUTPUT_PNG.with_suffix(".pdf")

METRICS = [
    ("genome_fraction_pct", "Genome fraction (%)", "a"),
    ("duplication_ratio", "Duplication ratio", "b"),
    ("mismatches_per_100kbp", "Mismatches / 100 kbp", "c"),
    ("indels_per_100kbp", "Indels / 100 kbp", "d"),
]

Y_MAX = {
    "genome_fraction_pct": 100.0,
    "duplication_ratio": 2.4,
}


def main() -> int:
    data = load_quast_values([column for column, _, _ in METRICS])
    print_source_audit(data, METRICS)

    apply_style()
    figure, axes = plt.subplots(nrows=2, ncols=2, figsize=(FULL_WIDTH_IN * 0.72, 5.4))

    for axis, (metric, label, letter) in zip(axes.flat, METRICS):
        grouped_bars(axis, data, metric)
        y_max = Y_MAX.get(metric, float(data[metric].max()) * 1.12)
        style_axis(axis, label, y_max, metric)
        panel_letter(axis, letter)
        if metric == "duplication_ratio":
            axis.axhline(1.0, color=NEUTRAL_GRAY, linewidth=0.7, linestyle=(0, (3, 2)), zorder=4)

    panel_legend(figure, sample_legend_handles(), y=0.965)
    figure.patch.set_facecolor("white")
    figure.subplots_adjust(left=0.10, right=0.99, bottom=0.10, top=0.91, hspace=0.62, wspace=0.34)
    save_figure(figure, OUTPUT_PNG, OUTPUT_PDF)

    print(f"\nPNG: {OUTPUT_PNG}\nPDF: {OUTPUT_PDF}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
