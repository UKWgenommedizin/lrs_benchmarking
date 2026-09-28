#!/usr/bin/env python3
"""Final QUAST-metrics figure -- built on assembler_benchmark_points.py.

Not a redesign: reuses assembler_benchmark_points.py's point-drawing
(dodged_points) and axis styling (style_axis) directly, so the point style,
sample colors, mean markers, SINGLE-TECHNOLOGY/HYBRID separation, legend,
typography and layout are identical to that script, not reimplemented.

Panels (a-f): NGA50 (Mb), genome fraction (%), duplication ratio,
misassemblies, mismatches / 100 kbp, indels / 100 kbp -- QUAST
reference-based metrics plus duplication ratio (raw N50 / contig count, the
FASTA-only stats used in assembler_benchmark_bars.py, are intentionally not
here; this is the QUAST-metrics figure, not a repeat of that one).

Duplication ratio matters on its own: it is what explains why Flye-HiFi's
N50/contig-count look so different from Flye-ONT's (unpurged haplotig
duplication -- duplication ratio ~1.33 with 0 N-gaps, vs Flye-ONT's ~1.01)
and why Verkko's ~2.0-2.1 is a different, *intentional* thing (true diploid
representation). It was already a QUAST report.tsv field and already a
column in assembly_benchmark_30x.tsv -- extract_quast_metrics.py and
build_assembly_benchmark_table.py already captured it -- but no figure had
plotted it until now.

Data loading uses utils/benchmark_data.py's load_quast_values(), which is
stricter than the equality check in verify_plot_data_matches_source(): for
every row it re-parses that row's own raw QUAST report.tsv (named in
source_file) and compares the field named in QUAST_REPORT_FIELDS against
the table value directly, not table-vs-table. Prints the source audit
(utils/benchmark_data.py's print_source_audit) before saving.
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

import assembler_benchmark_points as points_base
from utils.benchmark_data import QUAST_FINAL_METRICS, load_quast_values, print_source_audit
from utils.plot_style import (
    FULL_WIDTH_IN,
    apply_style,
    mean_marker_handle,
    panel_legend,
    panel_letter,
    sample_marker_handles,
    save_figure,
)

OUTPUT_PNG = PROJECT / "assembly_analysis" / "figures" / "30x" / "final" / "assembler_quast_final_points.png"
OUTPUT_PDF = OUTPUT_PNG.with_suffix(".pdf")
PLOT_DATA_TSV = (
    PROJECT / "assembly_analysis" / "tables" / "30x" / "derived" / "plot_data"
    / "assembler_quast_final_points_data.tsv"
)

RAW_COLUMNS = [
    "nga50", "genome_fraction_pct", "duplication_ratio",
    "misassemblies", "mismatches_per_100kbp", "indels_per_100kbp",
]


def main() -> int:
    data = load_quast_values(RAW_COLUMNS)  # raises on any table-vs-report.tsv mismatch
    data["nga50_mb"] = data["nga50"] / 1_000_000  # unit conversion for the axis, not a rescale

    print_source_audit(data, QUAST_FINAL_METRICS)

    PLOT_DATA_TSV.parent.mkdir(parents=True, exist_ok=True)
    columns = ["sample", "assembler", "technology", *RAW_COLUMNS, "nga50_mb", "source_file"]
    data[columns].to_csv(PLOT_DATA_TSV, sep="\t", index=False)

    y_limits = {
        column: float(data[column].max()) * 1.12 for column, _, _ in QUAST_FINAL_METRICS
    }

    apply_style()
    figure, axes = plt.subplots(nrows=2, ncols=3, figsize=(FULL_WIDTH_IN, 6.2))

    for axis, (metric, label, letter) in zip(axes.flat, QUAST_FINAL_METRICS):
        points_base.dodged_points(axis, data, metric)
        points_base.style_axis(axis, label, y_limits[metric])
        panel_letter(axis, letter)

    panel_legend(
        figure, sample_marker_handles() + [mean_marker_handle()], y=0.965, ncols=4,
    )

    figure.patch.set_facecolor("white")
    figure.subplots_adjust(left=0.06, right=0.99, bottom=0.08, top=0.92, hspace=0.55, wspace=0.30)

    save_figure(figure, OUTPUT_PNG, OUTPUT_PDF)

    print(f"\nRows: {len(data)}")
    print(f"Plot data: {PLOT_DATA_TSV}")
    print(f"PNG: {OUTPUT_PNG}")
    print(f"PDF: {OUTPUT_PDF}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
