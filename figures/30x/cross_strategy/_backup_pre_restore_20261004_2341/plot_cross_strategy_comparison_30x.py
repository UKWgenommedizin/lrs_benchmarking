#!/usr/bin/env python3
"""Figure 5: cross-strategy comparison of alignment- and assembly-based analysis.

3 rows x 2 columns (30x, GIAB HG002/HG003/HG004, ONT and PacBio HiFi):

    row a  Sequence disagreement metrics   reads: NM edit distance | contigs: QUAST mismatches
    row b  Indel disagreement              reads: CIGAR I+D events | contigs: QUAST indel events
    row c  Reference representation        reads: CIGAR-aligned yield | contigs: genome fraction

Left column = alignment-based, right column = assembly-based. Each subpanel has
its own y-axis because no row pairs two numerically identical quantities.

Every value is read from the curated benchmark tables; nothing is estimated,
interpolated or imputed, and no statistical test is run between strategies.

Column mapping (all verified on load, see load_alignment()):

  a, reads    samtools stats `mismatches` / `bases_mapped_cigar` x 1e5
              (final/alignment_benchmark_30x.tsv). samtools takes `mismatches`
              from the NM tag, i.e. edit distance = substitutions + inserted +
              deleted bases. No substitution-only column exists: NM minus the
              1-300 bp indel bases (`inserted_bases`, `deleted_bases`) is not
              one either, because CIGAR I/D operations >300 bp are outside the
              samtools ID histogram and stay in the remainder.
  a, contigs  QUAST `mismatches_per_100kbp` (substitutions only, per 100 kbp
              of aligned contig bases) (assembly_benchmark_30x.tsv).
  b, reads    (`insertion_events` + `deletion_events`) / `bases_mapped_cigar`
              x 1e5 (derived/alignment_benchmark_30x_indel_recovery.tsv).
              CIGAR I/D OPERATIONS (events), 1-300 bp, not indel bases.
  b, contigs  QUAST `indels_per_100kbp` (events; consecutive indel positions
              are one indel).
  c, reads    `input_normalized_cigar_yield_percent`
              (derived/plot_data/input_normalized_cigar_yield_percent.tsv):
              CIGAR-mapped read bases / input sequencing bases (largest samtools
              `total_length` among the four aligners per sample and technology).
  c, contigs  QUAST `genome_fraction_pct` (% of GRCh38 covered).

samtools stats was run with no -F/-f/-q/-d options: secondary alignments are
skipped, primary and supplementary records counted, MAPQ 0 included.

Outputs (figures/30x/cross_strategy/):
  cross_strategy_comparison_30x.{pdf,svg,png}       main figure (PNG 600 dpi)
  panels/cross_strategy_comparison_30x_panel_{A,B,C}.{pdf,svg,png}
  cross_strategy_comparison_30x_panel_data.tsv      tidy plotted values
  cross_strategy_comparison_30x_summary.tsv         mean / SD / n per configuration
  cross_strategy_comparison_30x_warnings.txt        interpretation warnings
"""

from __future__ import annotations

import os
from pathlib import Path
import sys

MPL_CACHE = Path("/tmp/matplotlib-lrs-benchmarking")
MPL_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_CACHE))

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.legend_handler import HandlerTuple
    from matplotlib.lines import Line2D
    from matplotlib.ticker import FuncFormatter, MultipleLocator
    import numpy as np
    import pandas as pd
except ModuleNotFoundError as error:
    raise SystemExit(
        f"Missing plotting package: {error.name}\nPython: {sys.executable}"
    ) from error


def find_project_root() -> Path:
    """Find the repository without depending on a username or script depth."""
    for parent in Path(__file__).resolve().parents:
        if (parent / "alignment_analysis").is_dir() and (parent / "assembly_analysis").is_dir():
            return parent
    raise RuntimeError("Could not locate the lrs_benchmarking project root")


PROJECT = find_project_root()
sys.path.insert(0, str(PROJECT / "alignment_analysis" / "scripts"))
from utils.plot_style import (  # noqa: E402
    FULL_WIDTH_IN,
    SAMPLE_COLORS,
    SAMPLE_ORDER,
    apply_style,
)

ALN_TABLES = PROJECT / "alignment_analysis" / "tables" / "30x"
ALN_FINAL_TSV = ALN_TABLES / "final" / "alignment_benchmark_30x.tsv"
ALN_INDEL_TSV = ALN_TABLES / "derived" / "alignment_benchmark_30x_indel_recovery.tsv"
ALN_YIELD_TSV = ALN_TABLES / "derived" / "plot_data" / "input_normalized_cigar_yield_percent.tsv"
ASM_TSV = PROJECT / "assembly_analysis" / "tables" / "30x" / "final" / "assembly_benchmark_30x.tsv"

OUT_DIR = PROJECT / "figures" / "30x" / "cross_strategy"
PANEL_DIR = OUT_DIR / "panels"
STEM = "cross_strategy_comparison_30x"
TIDY_TSV = OUT_DIR / f"{STEM}_panel_data.tsv"
SUMMARY_TSV = OUT_DIR / f"{STEM}_summary.tsv"
WARNINGS_TXT = OUT_DIR / f"{STEM}_warnings.txt"

ALIGNERS = ["minimap2", "pbmm2", "VACmap", "VG Giraffe"]
ALN_TECHNOLOGIES = ["ONT", "PacBio"]
# Grouped by sequencing input (ONT, HiFi, then hybrid), as for the aligners.
ASM_CONFIGS = [
    ("flye", "ont"),
    ("goldrush", "ont"),
    ("flye", "pb"),
    ("goldrush", "pb"),
    ("verkko", "hybrid"),
]
ASM_NAMES = {"flye": "Flye", "goldrush": "GoldRush", "verkko": "Verkko"}
ASM_REPRESENTATION = {"flye": "collapsed haploid", "goldrush": "collapsed haploid",
                      "verkko": "diploid (both haplotypes)"}
TECH_LABEL = {"ONT": "ONT", "PacBio": "HiFi", "ont": "ONT", "pb": "HiFi", "hybrid": "ONT+HiFi"}
TICK_NAMES = {"VG Giraffe": "VG\nGiraffe"}  # horizontal tick labels, two lines

STRATEGY_ALN = "Alignment-based"
STRATEGY_ASM = "Assembly-based"

# ---------------------------------------------------------------- typography
# Sizes are on the 7.16 in (journal double-column) canvas. At the thesis
# include width (\linewidth = 6.5 in) they print at x0.91, so ticks stay >5 pt.
HEADER_SIZE = 9.5
LETTER_SIZE = 10.5
ROW_TITLE_SIZE = 8.0
SUBTITLE_SIZE = 7.0
AXIS_LABEL_SIZE = 7.0
TICK_SIZE = 6.3
GROUP_SIZE = 6.5
NOTE_SIZE = 5.8
LEGEND_SIZE = 6.8
INK = "#1A1A1A"
INK_SECONDARY = "#4D4D4D"
INK_MUTED = "#7A7A7A"
AXIS_LW = 0.6

# ---------------------------------------------------------------- marks
# Sample = color (shared Okabe-Ito mapping); sequencing input = marker shape.
TECH_MARKER = {"ONT": "o", "HiFi": "s", "ONT+HiFi": "D"}
TECH_MARKER_AREA = {"o": 30, "s": 26, "D": 22}  # equal visual weight per shape
SAMPLE_OFFSETS = dict(zip(SAMPLE_ORDER, [-0.27, -0.11, 0.05]))
SUMMARY_OFFSET = 0.25   # mean +/- SD drawn beside, not over, the samples
MEAN_HALF_WIDTH = 0.085
MEAN_LW = 1.5
SD_LW = 0.7
SD_COLOR = "#8A8A8A"
GROUP_GAP = 0.5         # extra x space between sequencing-input groups

# Verkko: very light band with sparse hatching, so the points stay dominant.
VERKKO_CONFIG = "Verkko (ONT+HiFi)"
VERKKO_FACE = "#F6F6F6"
VERKKO_HATCH_EDGE = "#E2E2E2"
VERKKO_HATCH = "//"
# Wording note: the recorded Verkko command used --hifi/--nano only (no trio or
# Hi-C), so the assembly holds both haplotypes but they are not phased into two
# haplotype sets. Edit here if you prefer "haplotype-resolved".
VERKKO_NOTE = "diploid,\nboth haplotypes"

WARNINGS: list[str] = []


def warn(message: str) -> None:
    WARNINGS.append(message)
    print(f"WARNING: {message}", file=sys.stderr)


# ============================================================
# LOADING
# ============================================================


def read_tsv(path: Path, required: set[str]) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Input table not found: {path}")
    data = pd.read_csv(path, sep="\t")
    missing = required - set(data.columns)
    if missing:
        raise ValueError(f"{path.name}: missing required columns {sorted(missing)}")
    return data


def normalise_aligner_keys(data: pd.DataFrame) -> pd.DataFrame:
    data = data.copy()
    data["aligner"] = data["aligner"].replace({"VACMap": "VACmap", "VG": "VG Giraffe"})
    data["read_technology"] = data["read_technology"].replace(
        {"PacBio HiFi": "PacBio", "PB": "PacBio"}
    )
    keys = ["sample", "read_technology", "aligner"]
    data = data[
        data["sample"].isin(SAMPLE_ORDER)
        & data["read_technology"].isin(ALN_TECHNOLOGIES)
        & data["aligner"].isin(ALIGNERS)
    ]
    duplicated = data.duplicated(keys, keep=False)
    if duplicated.any():
        raise ValueError("Duplicate aligner rows:\n" + data.loc[duplicated, keys].to_string())
    expected = pd.MultiIndex.from_product([SAMPLE_ORDER, ALN_TECHNOLOGIES, ALIGNERS], names=keys)
    # Reindex so an absent row becomes NaN (drawn as "n/a"), never a value.
    return data.set_index(keys).reindex(expected).reset_index()


def load_alignment() -> pd.DataFrame:
    final = normalise_aligner_keys(read_tsv(
        ALN_FINAL_TSV,
        {"sample", "read_technology", "aligner", "mismatches", "bases_mapped_cigar", "error_rate"},
    ))
    indel = normalise_aligner_keys(read_tsv(
        ALN_INDEL_TSV,
        {"sample", "read_technology", "aligner", "insertion_events", "deletion_events",
         "insertion_events_per_100kb", "deletion_events_per_100kb", "bases_mapped_cigar"},
    ))
    yield_ = normalise_aligner_keys(read_tsv(
        ALN_YIELD_TSV,
        {"sample", "read_technology", "aligner", "common_input_bases", "bases_mapped_cigar",
         "input_normalized_cigar_yield_percent"},
    ))

    keys = ["sample", "read_technology", "aligner"]
    data = final[keys + ["mismatches", "bases_mapped_cigar", "error_rate"]]
    data = data.merge(
        indel[keys + ["insertion_events", "deletion_events", "insertion_events_per_100kb",
                      "deletion_events_per_100kb", "bases_mapped_cigar"]]
        .rename(columns={"bases_mapped_cigar": "bases_mapped_cigar_indel_table"}),
        on=keys, validate="one_to_one",
    )
    data = data.merge(
        yield_[keys + ["common_input_bases", "input_normalized_cigar_yield_percent",
                       "bases_mapped_cigar"]]
        .rename(columns={"bases_mapped_cigar": "bases_mapped_cigar_yield_table"}),
        on=keys, validate="one_to_one",
    )

    # The three source tables must share the denominator, and every stored
    # rate must reproduce from its raw counts.
    for other in ["bases_mapped_cigar_indel_table", "bases_mapped_cigar_yield_table"]:
        both = data[["bases_mapped_cigar", other]].dropna()
        if not np.array_equal(both.iloc[:, 0].to_numpy(), both.iloc[:, 1].to_numpy()):
            raise ValueError(f"bases_mapped_cigar disagrees between source tables ({other})")
    rate = data["mismatches"] / data["bases_mapped_cigar"]
    if not np.allclose(rate, data["error_rate"], rtol=1e-3, equal_nan=True):
        raise ValueError("mismatches / bases_mapped_cigar does not reproduce samtools error_rate")
    for event in ["insertion", "deletion"]:
        recomputed = data[f"{event}_events"] / data["bases_mapped_cigar"] * 1e5
        if not np.allclose(recomputed, data[f"{event}_events_per_100kb"], rtol=1e-6,
                           equal_nan=True):
            raise ValueError(f"{event}_events_per_100kb not reproducible from raw counts")
    recomputed_yield = data["bases_mapped_cigar"] / data["common_input_bases"] * 100
    if not np.allclose(recomputed_yield, data["input_normalized_cigar_yield_percent"],
                       rtol=1e-9, equal_nan=True):
        raise ValueError("input_normalized_cigar_yield_percent not reproducible")

    data["strategy"] = STRATEGY_ALN
    data["method"] = data["aligner"]
    data["technology"] = data["read_technology"].map(TECH_LABEL)
    data["configuration"] = data["method"] + " (" + data["technology"] + ")"
    data["representation"] = "read alignment"
    return data


def load_assembly() -> pd.DataFrame:
    data = read_tsv(
        ASM_TSV,
        {"assembler", "sample", "technology", "genome_fraction_pct",
         "mismatches_per_100kbp", "indels_per_100kbp", "duplication_ratio"},
    )
    keys = ["assembler", "sample", "technology"]
    duplicated = data.duplicated(keys, keep=False)
    if duplicated.any():
        raise ValueError("Duplicate assembly rows:\n" + data.loc[duplicated, keys].to_string())
    expected = pd.MultiIndex.from_tuples(
        [(asm, sample, tech) for asm, tech in ASM_CONFIGS for sample in SAMPLE_ORDER], names=keys
    )
    data = data.set_index(keys).reindex(expected).reset_index()
    data["strategy"] = STRATEGY_ASM
    data["method"] = data["assembler"].map(ASM_NAMES)
    data["technology"] = data["technology"].map(TECH_LABEL)
    data["configuration"] = data["method"] + " (" + data["technology"] + ")"
    data["representation"] = data["assembler"].map(ASM_REPRESENTATION)
    return data


# ============================================================
# TIDY PANEL DATA
# ============================================================

TIDY_COLUMNS = [
    "panel", "column", "strategy", "method", "technology", "configuration", "sample",
    "representation", "metric", "y_label", "value", "unit", "count_type", "numerator",
    "denominator", "source_file", "source_columns",
]

# Each row of the figure: the metric shown in each column, its axis label and
# its exact definition (the definitions also go into the tidy table).
PANELS = {
    "A": dict(
        title="Sequence disagreement metrics",
        aln=dict(
            value="nm_per_100kbp", metric="NM edit distance per 100 kbp",
            subtitle="Reads · NM edit distance (substitutions + indel bases)",
            y_label="Edit distance /\n100 kbp aligned bases", unit="per 100 kbp",
            count_type="bases (substituted + inserted + deleted)",
            numerator="samtools stats mismatches = summed NM tag",
            denominator="CIGAR-mapped read bases (bases_mapped_cigar: M/=/X/I)",
            source=ALN_FINAL_TSV, columns="mismatches; bases_mapped_cigar"),
        asm=dict(
            value="mismatches_per_100kbp", metric="QUAST mismatches per 100 kbp",
            subtitle="Contigs · substitutions only",
            y_label="QUAST mismatches /\n100 kbp aligned contig bases", unit="per 100 kbp",
            count_type="substituted bases",
            numerator="QUAST # mismatches (substitutions only)",
            denominator="QUAST total aligned contig length",
            source=ASM_TSV, columns="mismatches_per_100kbp"),
        ylim=None,
    ),
    "B": dict(
        title="Indel disagreement",
        aln=dict(
            value="indel_events_per_100kbp", metric="CIGAR indel events per 100 kbp",
            subtitle="Reads", y_label="Indel events / 100 kbp", unit="per 100 kbp",
            count_type="events (CIGAR I/D operations, 1-300 bp)",
            numerator="insertion_events + deletion_events (samtools stats ID records)",
            denominator="CIGAR-mapped read bases (bases_mapped_cigar: M/=/X/I)",
            source=ALN_INDEL_TSV, columns="insertion_events; deletion_events; bases_mapped_cigar"),
        asm=dict(
            value="indels_per_100kbp", metric="QUAST indel events per 100 kbp",
            subtitle="Contigs", y_label="QUAST indels / 100 kbp", unit="per 100 kbp",
            count_type="events (consecutive indel positions counted once)",
            numerator="QUAST # indels",
            denominator="QUAST total aligned contig length",
            source=ASM_TSV, columns="indels_per_100kbp"),
        ylim=None,
    ),
    "C": dict(
        title="Reference representation",
        aln=dict(
            value="input_normalized_cigar_yield_percent", metric="CIGAR-aligned yield",
            subtitle="Reads", y_label="CIGAR-aligned yield (%)", unit="%",
            count_type="bases",
            numerator="CIGAR-mapped read bases (bases_mapped_cigar: M/=/X/I)",
            denominator="input sequencing bases (common_input_bases: largest samtools "
                        "total_length among the four aligners)",
            source=ALN_YIELD_TSV, columns="input_normalized_cigar_yield_percent"),
        asm=dict(
            value="genome_fraction_pct", metric="Genome fraction",
            subtitle="Contigs", y_label="Genome fraction (%)", unit="%",
            count_type="reference bases",
            numerator="GRCh38 bases covered by aligned contigs",
            denominator="GRCh38 reference length",
            source=ASM_TSV, columns="genome_fraction_pct"),
        ylim=(0, 100),
    ),
}


def tidy(frame: pd.DataFrame, panel: str, column: str, spec: dict) -> pd.DataFrame:
    out = frame[["strategy", "method", "technology", "configuration", "sample",
                 "representation"]].copy()
    out["value"] = frame[spec["value"]].to_numpy(dtype=float)
    out["panel"] = panel.lower()
    out["column"] = column
    for key in ["metric", "unit", "count_type", "numerator", "denominator"]:
        out[key] = spec[key]
    out["y_label"] = spec["y_label"].replace("\n", " ")
    out["source_file"] = str(spec["source"].relative_to(PROJECT))
    out["source_columns"] = spec["columns"]
    return out[TIDY_COLUMNS]


def build_panel_data(aln: pd.DataFrame, asm: pd.DataFrame) -> pd.DataFrame:
    aln = aln.copy()
    aln["nm_per_100kbp"] = aln["mismatches"] / aln["bases_mapped_cigar"] * 1e5
    aln["indel_events_per_100kbp"] = (
        (aln["insertion_events"] + aln["deletion_events"]) / aln["bases_mapped_cigar"] * 1e5
    )
    parts = []
    for key, spec in PANELS.items():
        parts.append(tidy(aln, key, "left", spec["aln"]))
        parts.append(tidy(asm, key, "right", spec["asm"]))
    panel_data = pd.concat(parts, ignore_index=True)
    for _, row in panel_data[panel_data["value"].isna()].iterrows():
        warn(f"Panel {row['panel']}: no value for {row['configuration']} {row['sample']} "
             "(left blank, marked n/a)")
    return panel_data


def summarise(panel_data: pd.DataFrame) -> pd.DataFrame:
    """Descriptive mean / SD (ddof=1) / n across the three samples, as drawn."""
    keys = ["panel", "strategy", "configuration", "metric", "unit"]
    summary = (panel_data.groupby(keys, sort=False)["value"]
               .agg(mean="mean", sd=lambda v: v.std(ddof=1), n="count",
                    min="min", max="max").reset_index())
    return summary


# ============================================================
# PLOTTING
# ============================================================


def split_config(config: str) -> tuple[str, str]:
    """'VG Giraffe (ONT)' -> ('VG Giraffe', 'ONT')."""
    method, tech = config.rsplit(" (", 1)
    return method, tech.rstrip(")")


def group_positions(configs: list[str]) -> tuple[np.ndarray, list[tuple[str, float, float]]]:
    """x positions with a gap between input groups, plus (input, first x, last x)."""
    positions, groups, x = [], [], 0.0
    for i, config in enumerate(configs):
        tech = split_config(config)[1]
        if i and tech != groups[-1][0]:
            x += GROUP_GAP
        if not groups or tech != groups[-1][0]:
            groups.append([tech, x, x])
        groups[-1][2] = x
        positions.append(x)
        x += 1.0
    return np.array(positions), [tuple(g) for g in groups]


def aln_configs() -> list[str]:
    return [f"{a} ({TECH_LABEL[t]})" for t in ALN_TECHNOLOGIES for a in ALIGNERS]


def asm_configs() -> list[str]:
    return [f"{ASM_NAMES[a]} ({TECH_LABEL[t]})" for a, t in ASM_CONFIGS]


X_PAD = 0.55  # data units of x padding outside the first and last slot


def x_span(configs: list[str]) -> float:
    positions, _ = group_positions(configs)
    return positions[-1] - positions[0] + 2 * X_PAD


def draw_strip(axis: plt.Axes, data: pd.DataFrame, configs: list[str]) -> None:
    """Per-sample points, with the mean +/- SD summary drawn beside them."""
    x_positions, axis.input_groups = group_positions(configs)
    for x, config in zip(x_positions, configs):
        subset = data[data["configuration"] == config]
        if config == VERKKO_CONFIG:
            axis.axvspan(x - 0.47, x + 0.47, facecolor=VERKKO_FACE,
                         edgecolor=VERKKO_HATCH_EDGE, hatch=VERKKO_HATCH, linewidth=0,
                         zorder=0)
            axis.verkko_x = x
        values = subset["value"].to_numpy(dtype=float)
        finite = values[np.isfinite(values)]
        if len(finite) >= 2:
            mean, sd = finite.mean(), finite.std(ddof=1)
            xs = x + SUMMARY_OFFSET
            axis.vlines(xs, mean - sd, mean + sd, color=SD_COLOR, linewidth=SD_LW,
                        zorder=2, capstyle="butt")
            axis.hlines(mean, xs - MEAN_HALF_WIDTH, xs + MEAN_HALF_WIDTH, color=INK,
                        linewidth=MEAN_LW, zorder=2.5, capstyle="butt")
        for _, row in subset.iterrows():
            xp = x + SAMPLE_OFFSETS[row["sample"]]
            marker = TECH_MARKER[row["technology"]]
            if np.isfinite(row["value"]):
                axis.scatter(xp, row["value"], s=TECH_MARKER_AREA[marker], marker=marker,
                             color=SAMPLE_COLORS[row["sample"]], edgecolor="white",
                             linewidth=0.5, zorder=3, clip_on=False)
            else:
                axis.text(xp, 0, "n/a", rotation=90, fontsize=5, color=INK_MUTED,
                          ha="center", va="bottom", zorder=3)

    axis.set_xticks(x_positions)
    axis.set_xticklabels([TICK_NAMES.get(split_config(c)[0], split_config(c)[0])
                          for c in configs], fontsize=TICK_SIZE, color=INK,
                         linespacing=0.95, va="top")
    axis.tick_params(axis="x", length=0, pad=3)
    axis.tick_params(axis="y", labelsize=TICK_SIZE, length=2.5, width=AXIS_LW, pad=2,
                     colors=INK)
    axis.set_xlim(x_positions[0] - X_PAD, x_positions[-1] + X_PAD)
    for side in ("top", "right"):
        axis.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        axis.spines[side].set_linewidth(AXIS_LW)
        axis.spines[side].set_color(INK)


def nice_top(value: float) -> tuple[float, float]:
    """(upper limit, tick step) with ~4 intervals of a 1/2/2.5/5 x 10^k step."""
    raw_step = value / 4
    magnitude = 10 ** np.floor(np.log10(raw_step))
    step = next(m * magnitude for m in (1, 2, 2.5, 5, 10) if m * magnitude >= raw_step)
    return float(np.ceil(value / step) * step), float(step)


def data_top(data: pd.DataFrame) -> float:
    """Highest point or mean + SD whisker."""
    top = np.nanmax(data["value"].to_numpy(dtype=float))
    for _, group in data.groupby("configuration"):
        finite = group["value"].dropna().to_numpy(dtype=float)
        if len(finite) >= 2:
            top = max(top, finite.mean() + finite.std(ddof=1))
    return float(top)


def thousands(value: float, _pos) -> str:
    return f"{value:,.0f}"


def set_y_axis(axis: plt.Axes, data: pd.DataFrame, ylim) -> None:
    """Linear axis from zero; fixed 0-100 for percentages."""
    if ylim is not None:
        axis.set_ylim(*ylim)
        axis.yaxis.set_major_locator(MultipleLocator(25))
    else:
        top, step = nice_top(data_top(data) * 1.04)
        axis.set_ylim(0, top)
        axis.yaxis.set_major_locator(MultipleLocator(step))
    axis.yaxis.set_major_formatter(FuncFormatter(thousands))


def add_group_brackets(figure: plt.Figure) -> None:
    """Thin sequencing-input bracket + label under each axis's tool names.

    Needs the rendered tick-label extents, so it runs just before saving.
    """
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    grouped = [axis for axis in figure.axes if getattr(axis, "input_groups", None)]
    # Axes in the same row share one bracket height (lowest tick label in the row).
    row_bottom: dict[float, float] = {}
    for axis in grouped:
        row = round(axis.get_position().y0, 3)
        boxes = [label.get_window_extent(renderer) for label in axis.get_xticklabels()
                 if label.get_text()]
        row_bottom[row] = min(row_bottom.get(row, np.inf), min(box.y0 for box in boxes))
    tick_px = 2.0 * figure.dpi / 72  # bracket end ticks, 2 pt
    for axis in grouped:
        bottom_px = row_bottom[round(axis.get_position().y0, 3)]
        to_axes = axis.transAxes.inverted()
        y_line = to_axes.transform((0, bottom_px - 3.0 * figure.dpi / 72))[1]
        y_tick = to_axes.transform((0, bottom_px - 3.0 * figure.dpi / 72 + tick_px))[1]
        transform = axis.get_xaxis_transform()  # x in data, y in axes fraction
        for tech, x0, x1 in axis.input_groups:
            xa, xb = x0 - 0.36, x1 + 0.36
            axis.add_line(Line2D([xa, xa, xb, xb], [y_tick, y_line, y_line, y_tick],
                                 transform=transform, color=INK_MUTED, linewidth=0.6,
                                 solid_joinstyle="miter", clip_on=False))
            axis.annotate(tech, xy=((x0 + x1) / 2, y_line), xycoords=transform,
                          xytext=(0, -1.5), textcoords="offset points", ha="center",
                          va="top", fontsize=GROUP_SIZE, color=INK_SECONDARY,
                          annotation_clip=False)


def legend_handles() -> tuple[list, list]:
    """Two legend groups: sample (color) and technology / summary (shape, mean +/- SD)."""
    samples = [Line2D([], [], linestyle="", marker="o", markersize=5.2,
                      markerfacecolor=SAMPLE_COLORS[s], markeredgecolor="white",
                      markeredgewidth=0.5, label=s) for s in SAMPLE_ORDER]
    techs = [Line2D([], [], linestyle="", marker=TECH_MARKER[t],
                    markersize={"o": 5.2, "s": 4.8, "D": 4.4}[TECH_MARKER[t]],
                    markerfacecolor=INK_MUTED, markeredgecolor="white",
                    markeredgewidth=0.5, label=t) for t in ["ONT", "HiFi", "ONT+HiFi"]]
    sd = Line2D([], [], linestyle="", marker="|", markersize=8, markeredgewidth=SD_LW,
                markeredgecolor=SD_COLOR)
    mean = Line2D([], [], linestyle="", marker="_", markersize=6, markeredgewidth=MEAN_LW,
                  markeredgecolor=INK)
    return samples, techs + [(sd, mean)]


def add_legends(figure: plt.Figure, y_top: float) -> None:
    samples, techs = legend_handles()
    common = dict(frameon=False, fontsize=LEGEND_SIZE, title_fontsize=LEGEND_SIZE,
                  handletextpad=0.3, columnspacing=1.1, borderaxespad=0, borderpad=0,
                  handlelength=1.2, alignment="left", loc="upper left")
    first = figure.legend(handles=samples, ncols=3, title="Sample",
                          bbox_to_anchor=(0.17, y_top), **common)
    second = figure.legend(handles=techs, labels=["ONT", "HiFi", "ONT+HiFi",
                                                  "Mean ± SD (n = 3)"],
                           ncols=4, title="Technology / summary",
                           handler_map={tuple: HandlerTuple(ndivide=1)},
                           bbox_to_anchor=(0.47, y_top), **common)
    for legend in (first, second):
        legend.get_title().set_fontweight("bold")
        legend.get_title().set_color(INK)


def plot_rows(panel_data: pd.DataFrame, keys: list[str], stem: Path) -> list[Path]:
    """Lay out len(keys) rows x 2 columns on an explicit inch grid and save."""
    width = FULL_WIDTH_IN
    left_margin, right_margin, col_gap = 0.74, 0.06, 0.80
    header_h, title_h, subtitle_h = 0.34, 0.20, 0.23
    axes_h, below_h, row_gap = 1.30, 0.42, 0.02
    legend_h = 0.42
    row_h = title_h + subtitle_h + axes_h + below_h
    height = header_h + len(keys) * row_h + (len(keys) - 1) * row_gap + legend_h
    figure = plt.figure(figsize=(width, height))

    # Column widths proportional to x span, so one tool slot has the same
    # physical width in both columns.
    span_l, span_r = x_span(aln_configs()), x_span(asm_configs())
    plot_w = width - left_margin - right_margin - col_gap
    w_l, w_r = plot_w * span_l / (span_l + span_r), plot_w * span_r / (span_l + span_r)
    x_l, x_r = left_margin, left_margin + w_l + col_gap

    def fig_x(inches: float) -> float:
        return inches / width

    def fig_y(inches_from_top: float) -> float:
        return 1 - inches_from_top / height

    # Strategy headers: centred over each column, short rule just beneath.
    for x0, w, label in [(x_l, w_l, STRATEGY_ALN), (x_r, w_r, STRATEGY_ASM)]:
        figure.text(fig_x(x0 + w / 2), fig_y(0.20), label, ha="center", va="baseline",
                    fontsize=HEADER_SIZE, fontweight="bold", color=INK)
        figure.add_artist(Line2D([fig_x(x0), fig_x(x0 + w)], [fig_y(0.27)] * 2,
                                 transform=figure.transFigure, color=INK, linewidth=0.7))

    top = header_h
    for key in keys:
        spec = PANELS[key]
        # Panel letter and row title share one baseline, flush left.
        baseline = fig_y(top + title_h - 0.03)
        figure.text(fig_x(0.04), baseline, key.lower(), ha="left", va="baseline",
                    fontsize=LETTER_SIZE, fontweight="bold", color=INK)
        figure.text(fig_x(0.27), baseline, spec["title"], ha="left", va="baseline",
                    fontsize=ROW_TITLE_SIZE, fontweight="bold", color=INK)
        axes_top = top + title_h + subtitle_h
        bottom = fig_y(axes_top + axes_h)
        for side, x0, w, configs in [("aln", x_l, w_l, aln_configs()),
                                     ("asm", x_r, w_r, asm_configs())]:
            side_spec = spec[side]
            axis = figure.add_axes((fig_x(x0), bottom, fig_x(w), axes_h / height))
            strategy = STRATEGY_ALN if side == "aln" else STRATEGY_ASM
            data = panel_data[(panel_data["panel"] == key.lower())
                              & (panel_data["strategy"] == strategy)]
            draw_strip(axis, data, configs)
            set_y_axis(axis, data, spec["ylim"])
            axis.set_ylabel(side_spec["y_label"], fontsize=AXIS_LABEL_SIZE, color=INK,
                            labelpad=3, linespacing=1.05)
            # Raised clear of the top y tick label.
            axis.annotate(side_spec["subtitle"], xy=(0, 1), xycoords="axes fraction",
                          xytext=(0, 6), textcoords="offset points", ha="left",
                          va="bottom", fontsize=SUBTITLE_SIZE, color=INK_SECONDARY)
            if hasattr(axis, "verkko_x"):
                # Right-aligned to the band's right edge so it stays on the canvas.
                axis.annotate(VERKKO_NOTE, xy=(axis.verkko_x + 0.47, 1), xycoords=(
                    axis.transData, axis.transAxes), xytext=(0, 2),
                    textcoords="offset points", ha="right", va="bottom",
                    fontsize=NOTE_SIZE, fontstyle="italic", color=INK_MUTED,
                    linespacing=0.95, multialignment="right", annotation_clip=False)
        top = axes_top + axes_h + below_h + row_gap

    for axis in figure.axes:
        # Same y-label x in every row; anchor = right edge of the rotated text.
        axis.yaxis.set_label_coords(-0.36 / (axis.get_position().width * width), 0.5)

    add_legends(figure, fig_y(height - legend_h + 0.10))
    paths = save(figure, stem)
    plt.close(figure)
    return paths


def save(figure: plt.Figure, stem: Path) -> list[Path]:
    """Exact canvas size (no tight bbox), so the PDF is the journal width."""
    stem.parent.mkdir(parents=True, exist_ok=True)
    add_group_brackets(figure)
    paths = [stem.with_suffix(s) for s in (".pdf", ".svg", ".png")]
    for path in paths:
        figure.savefig(path, dpi=600, facecolor="white")
    return paths


# ============================================================
# INTERPRETATION WARNINGS (definitions verified against the source tables)
# ============================================================


def record_definition_warnings(aln: pd.DataFrame, asm: pd.DataFrame) -> None:
    residual = ((aln["mismatches"] - aln["inserted_bases"] - aln["deleted_bases"])
                / aln["bases_mapped_cigar"] * 1e5) if "inserted_bases" in aln else None
    warn(
        "Panel a: no substitution-only read metric exists in the source tables. samtools "
        "stats 'mismatches' is the summed NM tag (substitutions + inserted + deleted bases), "
        "so the left subpanel is NM edit distance per 100 kbp of CIGAR-mapped read bases, "
        "labelled 'Edit distance'. QUAST '# mismatches per 100 kbp' counts substitutions only "
        "per aligned contig base. The two are related but not the same quantity; separate "
        "y-axes, no merging, no test. NM minus the 1-300 bp indel bases is NOT a substitution "
        "count either, because CIGAR I/D operations >300 bp are outside the samtools ID "
        "histogram and stay in the remainder"
        + (f" (remainder: {residual.min():.0f}-{residual.max():.0f} per 100 kbp)"
           if residual is not None else "") + "."
    )
    warn(
        "Panel b: both columns are indel EVENT counts (not indel bases) per 100 kbp of aligned "
        "query bases. Reads: CIGAR I + D operations of 1-300 bp from samtools stats ID records "
        "(longer operations are not counted). Contigs: QUAST # indels, counted only inside "
        "alignments; indels >200 bp (QUAST default) become local misassemblies. Raw reads "
        "(each locus seen ~30x, per-read error) and consensus contigs are different query "
        "units, so the rates are shown on separate y-axes and not tested against each other."
    )
    warn(
        "Panel c: CIGAR-aligned yield (denominator: input sequencing bases) and genome "
        "fraction (denominator: GRCh38 length) answer different questions and have different "
        "denominators. They are drawn on the same 0-100% scale for readability only; they are "
        "not 'completeness' on both sides and are not compared statistically."
    )
    dup = asm.loc[asm["assembler"] == "verkko", "duplication_ratio"]
    warn(
        "Verkko: QUAST evaluated the combined assembly.fasta, which holds sequence from both "
        f"haplotypes (QUAST duplication ratio {dup.min():.2f}-{dup.max():.2f}), whereas Flye "
        "and GoldRush are collapsed haploid. The recorded command "
        "(assemblers/whole_genome_asm/hybrid.assembly.verkko.smk) uses only --hifi and --nano, "
        "with no trio or Hi-C phasing input: haplotypes are separated locally within graph "
        "bubbles but not phased into two chromosome-scale haplotype sets. The figure therefore "
        "says 'diploid, both haplotypes' rather than 'haplotype-resolved'. Its genome fraction "
        "can be raised by sequence from both haplotypes and is not equivalent to that of a "
        "collapsed assembly. Verkko input depth is 'unknown' in the source table."
    )
    warn(
        "Not compared by design: mapped reads % or mapped bases % vs genome fraction %; MQ0 "
        "vs misassemblies; secondary alignments vs duplication ratio; read N50 vs NGA50; "
        "CIGAR yield vs NGA50; NM edit distance vs QUAST mismatches."
    )
    warn(
        "Statistics: HG002/HG003/HG004 are three biological samples (n = 3 per "
        "configuration), not technical replicates; mean +/- SD (ddof = 1) is descriptive "
        "only and no tests are run."
    )


def main() -> int:
    apply_style()
    plt.rcParams["hatch.linewidth"] = 0.4
    aln = load_alignment()
    aln = aln.merge(
        normalise_aligner_keys(read_tsv(ALN_INDEL_TSV, {"inserted_bases", "deleted_bases"}))
        [["sample", "read_technology", "aligner", "inserted_bases", "deleted_bases"]],
        on=["sample", "read_technology", "aligner"], validate="one_to_one",
    )
    asm = load_assembly()
    record_definition_warnings(aln, asm)

    panel_data = build_panel_data(aln, asm)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    panel_data.to_csv(TIDY_TSV, sep="\t", index=False, float_format="%.6g")
    summarise(panel_data).to_csv(SUMMARY_TSV, sep="\t", index=False, float_format="%.6g")

    outputs = plot_rows(panel_data, list(PANELS), OUT_DIR / STEM)
    for key in PANELS:
        outputs += plot_rows(panel_data, [key], PANEL_DIR / f"{STEM}_panel_{key}")

    WARNINGS_TXT.write_text(
        "".join(f"{i}. {message}\n" for i, message in enumerate(WARNINGS, 1)),
        encoding="utf-8",
    )
    print(f"Panel data: {TIDY_TSV.relative_to(PROJECT)}")
    print(f"Summary:    {SUMMARY_TSV.relative_to(PROJECT)}")
    print(f"Warnings:   {WARNINGS_TXT.relative_to(PROJECT)}")
    for path in outputs:
        print(f"Figure:     {path.relative_to(PROJECT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
