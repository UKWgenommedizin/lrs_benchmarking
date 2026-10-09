#!/usr/bin/env python3
"""Plot the QUAST reference-agreement panel for Flye, Goldrush and Verkko.

Composite figure, 2x2 grid, one subplot per reference-agreement metric from
the benchmark workflow diagram:

    genome fraction (%), misassemblies, mismatches / 100 kbp, indels / 100 kbp

Each subplot is a grouped bar chart: one bar per HG002/HG003/HG004 sample,
grouped by assembly configuration (Flye-ONT, Flye-HiFi, Goldrush-ONT,
Goldrush-HiFi, Verkko-Hybrid) on the x-axis. Verkko is haplotype-resolved
and hybrid (combined ONT + PacBio HiFi input, no per-technology split) --
see assembly_analysis/README.md, "Scientific comparison rules".
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
    import pandas as pd
except ModuleNotFoundError as error:
    raise SystemExit(
        f"A required plotting package is missing: {error.name}\n"
        f"Python interpreter: {sys.executable}\n"
    ) from error

from utils.plot_style import (
    CONFIGURATION_LABELS,
    CONFIGURATION_ORDER,
    FULL_WIDTH_IN,
    HYBRID_SPAN_ALPHA,
    HYBRID_SPAN_COLOR,
    PANEL_LABEL_SIZE_PT,
    PANEL_TICK_SIZE,
    PANEL_TITLE_SIZE,
    SAMPLE_COLORS,
    SAMPLE_ORDER,
    apply_style,
    clean_spines,
    hybrid_span_handle,
    panel_legend,
    rotated_xticks,
    sample_legend_handles,
    save_figure,
)

INPUT_FILE = (
    PROJECT / "assembly_analysis" / "tables" / "30x" / "final"
    / "assembly_benchmark_30x.tsv"
)
OUTPUT_PNG = (
    PROJECT / "assembly_analysis" / "figures" / "30x" / "final"
    / "01_reference_agreement_panel_30x.png"
)
OUTPUT_PDF = OUTPUT_PNG.with_suffix(".pdf")

METRICS = [
    ("genome_fraction_pct", "Genome fraction (%)"),
    ("misassemblies", "Misassemblies"),
    ("mismatches_per_100kbp", "Mismatches / 100 kbp"),
    ("indels_per_100kbp", "Indels / 100 kbp"),
]

SAMPLE_OFFSETS = {"HG002": -0.27, "HG003": 0.0, "HG004": 0.27}
BAR_WIDTH = 0.26

# Verkko combines both input technologies into one assembly per sample (~2x
# the raw coverage of a single Flye/Goldrush ONT- or HiFi-only run) and is
# haplotype-resolved (diploid), aligned by QUAST against one haploid
# reference -- see assembly_analysis/README.md, "Scientific comparison
# rules" #1/#2/#4. Both of those are properties of the *comparison setup*
# that apply to every metric in this figure, not just some of them. An
# earlier version of this script marked only misassemblies and
# mismatches/100 kbp as affected (reasoning that heterozygous SNPs are more
# common than heterozygous indels, so a diploid-vs-haploid-reference
# comparison should inflate mismatches more) -- that is a plausible
# hypothesis from population genetics, not something measured on this
# dataset (it would need Verkko run single-technology, or QUAST run per
# haplotype, to actually isolate). Rather than imply a false precision with
# a per-panel marker, the caveat is applied uniformly: one dagger next to
# the legend swatch, one footnote, covering all four panels equally. This
# matches how phased/diploid-vs-collapsed assembly comparisons are usually
# handled in the literature -- either the reference-based comparison is
# avoided altogether (per-haplotype alignment, or reference-free k-mer QV),
# or, when it is reported, one blanket methodological caveat is stated for
# all reference-based metrics rather than a magnitude ranking per metric.
#
# No dagger/symbol here: the shaded column already marks every affected bar,
# in every panel, uniformly -- a symbol earns its place when it flags a
# *subset*; repeating it next to something already visually distinct in 100%
# of cases is a second notation system saying the same thing as the first.
# The footnote refers to "the shaded columns" directly instead.
FOOTNOTE = (
    "The shaded columns above are Verkko's hybrid assembly (30x ONT + 30x PacBio HiFi combined) and\n"
    "haplotype-resolved, but QUAST here aligns it to a single haploid reference: all four Verkko values\n"
    "should be read with this in mind, not as directly input- or ploidy-comparable to the\n"
    "single-technology Flye/Goldrush runs shown alongside them."
)


def load_plot_data() -> pd.DataFrame:
    if not INPUT_FILE.exists():
        raise FileNotFoundError(
            f"The input table was not found:\n{INPUT_FILE}\n"
            "Run assembly_analysis/scripts/metrics/extract_quast_metrics.py "
            "and build_assembly_benchmark_table.py first."
        )

    data = pd.read_csv(INPUT_FILE, sep="\t")
    metric_columns = [column for column, _ in METRICS]
    required_columns = {"assembler", "sample", "technology", *metric_columns}
    missing_columns = required_columns.difference(data.columns)
    if missing_columns:
        raise ValueError(f"The input table is missing these columns: {sorted(missing_columns)}")

    for column in metric_columns:
        data[column] = pd.to_numeric(data[column], errors="raise")

    duplicated = data.duplicated(subset=["assembler", "sample", "technology"], keep=False)
    if duplicated.any():
        raise ValueError(
            "Duplicated assembler/sample/technology observations found:\n"
            + data.loc[duplicated, ["assembler", "sample", "technology"]].to_string(index=False)
        )

    expected = pd.MultiIndex.from_tuples(
        [(assembler, technology, sample) for assembler, technology in CONFIGURATION_ORDER for sample in SAMPLE_ORDER],
        names=["assembler", "technology", "sample"],
    )
    observed = pd.MultiIndex.from_frame(data[["assembler", "technology", "sample"]])
    missing = expected.difference(observed)
    if len(missing):
        raise ValueError(f"Missing assembler/technology/sample combinations: {list(missing)}")

    return data


def grouped_bars(axis: plt.Axes, data: pd.DataFrame, metric: str) -> None:
    for config_index, (assembler, technology) in enumerate(CONFIGURATION_ORDER):
        subset = data[(data["assembler"] == assembler) & (data["technology"] == technology)]
        values = subset.set_index("sample")[metric].reindex(SAMPLE_ORDER)

        for sample in SAMPLE_ORDER:
            axis.bar(
                config_index + SAMPLE_OFFSETS[sample],
                values.loc[sample],
                width=BAR_WIDTH,
                color=SAMPLE_COLORS[sample],
                edgecolor="white",
                linewidth=0.3,
                zorder=3,
            )


VERKKO_CONFIG_INDEX = CONFIGURATION_ORDER.index(("verkko", "hybrid"))


def shade_hybrid_column(axis: plt.Axes) -> None:
    """Wash the Verkko column with a light neutral background.

    Distinguishes it from the single-technology Flye/Goldrush columns without
    a bracket-and-line glyph, which in a Nature/Cell/Science-style figure
    reads as a statistical-significance comparison, not a category grouping.
    zorder sits below the bars (3) and above the plain white axes background.
    """
    axis.axvspan(
        VERKKO_CONFIG_INDEX - 0.5, VERKKO_CONFIG_INDEX + 0.5,
        color=HYBRID_SPAN_COLOR, alpha=HYBRID_SPAN_ALPHA, zorder=0.5, linewidth=0,
    )


def style_axis(axis: plt.Axes, label: str, y_max: float) -> None:
    axis.set_ylim(0, y_max)
    axis.set_xlim(-0.5, len(CONFIGURATION_ORDER) - 0.5)
    positions = range(len(CONFIGURATION_ORDER))
    labels = [CONFIGURATION_LABELS[config] for config in CONFIGURATION_ORDER]
    rotated_xticks(axis, positions, labels, fontsize=PANEL_TICK_SIZE, tick_length=2)
    axis.tick_params(axis="y", labelsize=PANEL_TICK_SIZE, length=2, pad=1.5)
    axis.set_ylabel(label, fontsize=PANEL_LABEL_SIZE_PT, labelpad=2)
    clean_spines(axis)
    axis.set_facecolor("white")


def main() -> int:
    data = load_plot_data()
    apply_style()

    figure, axes = plt.subplots(nrows=2, ncols=2, figsize=(FULL_WIDTH_IN, 5.15))

    for axis, (metric, label) in zip(axes.flat, METRICS):
        shade_hybrid_column(axis)
        grouped_bars(axis, data, metric)
        y_max = float(data[metric].max()) * 1.12
        style_axis(axis, label, y_max)
        axis.set_title(label, pad=2, fontsize=PANEL_TITLE_SIZE)

    panel_legend(
        figure,
        sample_legend_handles() + [hybrid_span_handle()],
        y=0.958, ncols=4,
    )

    figure.text(
        0.5, 0.012, FOOTNOTE,
        ha="center", va="bottom", multialignment="center",
        fontsize=PANEL_TICK_SIZE - 0.5, color="#444444", style="italic",
    )

    figure.patch.set_facecolor("white")
    figure.subplots_adjust(left=0.07, right=0.99, bottom=0.16, top=0.905, hspace=0.6, wspace=0.22)

    save_figure(figure, OUTPUT_PNG, OUTPUT_PDF)

    print(f"Input: {INPUT_FILE}")
    print(f"PNG: {OUTPUT_PNG}")
    print(f"PDF: {OUTPUT_PDF}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
