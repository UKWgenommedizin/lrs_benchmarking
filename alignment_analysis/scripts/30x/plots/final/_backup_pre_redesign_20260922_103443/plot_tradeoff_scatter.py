#!/usr/bin/env python3
"""Within-strategy trade-off scatters (GoldRush-style 2x2 layout), 30x.

Top row = alignment-based, bottom row = assembly-based. Each scatter only
uses axes that are defined identically for every point in it, so no
aligner is ever placed on an assembler axis or vice versa.

  a  Aligners:   CIGAR-aligned yield (% of input bases) vs error rate (%)
                 (input_normalized_cigar_yield_percent, error_percent)
  b  Aligners:   peak RAM (GB) vs allocated thread-hours
                 (peak_ram_gb; wallclock_runtime_hours x threads).
                 Only configurations with MEASURED peak RAM are plotted.
  c  Assemblers: QUAST NGA50 (Mb, log) vs QUAST misassemblies
                 (nga50, misassemblies)
  d  Assemblers: peak RAM vs thread-hours, from
                 assembly_analysis/tables/assembler_performance_30x.tsv when
                 complete; otherwise the panel states that it is not measured.

Encoding (same as the thesis' other scatters): color = GIAB sample,
marker shape = tool, filled = ONT, open = PacBio HiFi, diamond = Verkko
(ONT+HiFi). Nothing is estimated or imputed.
"""

from __future__ import annotations

import os
from pathlib import Path
import sys

MPL_CACHE = Path("/tmp/matplotlib-lrs-benchmarking")
MPL_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_CACHE))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

PROJECT = next(
    parent for parent in Path(__file__).resolve().parents
    if (parent / "alignment_analysis").is_dir() and (parent / "assembly_analysis").is_dir()
)
sys.path.insert(0, str(PROJECT / "alignment_analysis" / "scripts"))
from utils.plot_style import (  # noqa: E402
    FULL_WIDTH_IN,
    SAMPLE_COLORS,
    SAMPLE_ORDER,
    apply_style,
    clean_spines,
    panel_letter,
    subtle_grid,
)

ALN_FINAL_TSV = PROJECT / "alignment_analysis" / "tables" / "30x" / "final" / "alignment_benchmark_30x.tsv"
ALN_YIELD_TSV = (PROJECT / "alignment_analysis" / "tables" / "30x" / "derived" / "plot_data"
                 / "input_normalized_cigar_yield_percent.tsv")
ASM_TSV = PROJECT / "assembly_analysis" / "tables" / "30x" / "final" / "assembly_benchmark_30x.tsv"
ASM_PERF_TSV = PROJECT / "assembly_analysis" / "tables" / "assembler_performance_30x.tsv"

OUT_DIR = PROJECT / "figures" / "30x" / "cross_strategy"
STEM = "tradeoff_scatter_30x"

ALIGNERS = ["minimap2", "pbmm2", "VACmap", "VG Giraffe"]
TECHNOLOGIES = ["ONT", "PacBio"]
ALIGNER_MARKERS = {"minimap2": "o", "pbmm2": "s", "VACmap": "^", "VG Giraffe": "v"}

ASM_CONFIGS = [("flye", "ont"), ("flye", "pb"), ("goldrush", "ont"), ("goldrush", "pb"),
               ("verkko", "hybrid")]
ASM_NAMES = {"flye": "Flye", "goldrush": "GoldRush", "verkko": "Verkko"}
ASM_MARKERS = {"flye": "o", "goldrush": "s", "verkko": "D"}
TECH_SHORT = {"ONT": "ONT", "PacBio": "HiFi", "ont": "ONT", "pb": "HiFi", "hybrid": "ONT+HiFi"}

LABEL_SIZE = 6.0


def draw_point(axis, x, y, marker, sample, filled):
    color = SAMPLE_COLORS[sample]
    axis.scatter(
        x, y, marker=marker, s=30, zorder=3, clip_on=False,
        facecolors=color if filled else "white",
        edgecolors="white" if filled else color,
        linewidths=0.5 if filled else 1.1,
    )


def label_cluster(axis, subset, xcol, ycol, text, offset, log_y=False):
    x = subset[xcol].mean()
    y = np.exp(np.log(subset[ycol]).mean()) if log_y else subset[ycol].mean()
    axis.annotate(text, (x, y), xytext=offset, textcoords="offset points",
                  fontsize=LABEL_SIZE, color="#333333", ha="center", va="center", zorder=4)


def label_at(axis, subset, xcol, ycol, text, position):
    """Label a cluster at a fixed data position; draw a leader line if asked."""
    x_text, y_text, ha, leader = position
    target = (subset[xcol].mean(), subset[ycol].mean())
    axis.annotate(
        text, target, xytext=(x_text, y_text), textcoords="data", fontsize=LABEL_SIZE,
        color="#333333", ha=ha, va="center", zorder=4,
        arrowprops=dict(arrowstyle="-", color="#9A9A9A", linewidth=0.5,
                        shrinkA=1, shrinkB=4) if leader else None,
    )


# ============================================================
# DATA
# ============================================================


def load_alignment() -> pd.DataFrame:
    final = pd.read_csv(ALN_FINAL_TSV, sep="\t")
    final["aligner"] = final["aligner"].replace({"VACMap": "VACmap"})
    yield_ = pd.read_csv(ALN_YIELD_TSV, sep="\t")
    keys = ["sample", "read_technology", "aligner"]
    data = final[keys + ["error_percent", "threads", "wallclock_runtime_hours", "peak_ram_gb",
                         "bases_mapped_cigar"]].merge(
        yield_[keys + ["input_normalized_cigar_yield_percent", "bases_mapped_cigar"]],
        on=keys, suffixes=("", "_yield"), validate="one_to_one")
    if not (data["bases_mapped_cigar"] == data["bases_mapped_cigar_yield"]).all():
        raise ValueError("bases_mapped_cigar differs between the alignment source tables")
    data = data[data["sample"].isin(SAMPLE_ORDER) & data["aligner"].isin(ALIGNERS)
                & data["read_technology"].isin(TECHNOLOGIES)]
    if len(data) != len(SAMPLE_ORDER) * len(ALIGNERS) * len(TECHNOLOGIES):
        raise ValueError(f"Expected 24 alignment rows, found {len(data)}")
    data["thread_hours"] = data["wallclock_runtime_hours"] * data["threads"]
    return data.drop(columns="bases_mapped_cigar_yield")


def load_assembly() -> pd.DataFrame:
    data = pd.read_csv(ASM_TSV, sep="\t")
    keep = pd.MultiIndex.from_tuples(ASM_CONFIGS)
    data = data[pd.MultiIndex.from_frame(data[["assembler", "technology"]]).isin(keep)
                & data["sample"].isin(SAMPLE_ORDER)].copy()
    if len(data) != len(ASM_CONFIGS) * len(SAMPLE_ORDER):
        raise ValueError(f"Expected 15 assembly rows, found {len(data)}")
    if data[["nga50", "misassemblies"]].isna().any().any():
        raise ValueError("Missing NGA50/misassemblies in assembly table")
    data["nga50_mb"] = data["nga50"] / 1e6
    return data


def load_assembly_performance() -> pd.DataFrame | None:
    if not ASM_PERF_TSV.is_file():
        return None
    perf = pd.read_csv(ASM_PERF_TSV, sep="\t")
    if len(perf) < len(ASM_CONFIGS) * len(SAMPLE_ORDER) or \
            perf[["threads", "wall_clock_hours", "peak_rss_gb"]].isna().any().any():
        return None
    perf["thread_hours"] = perf["wall_clock_hours"] * perf["threads"]
    return perf


# ============================================================
# PANELS
# ============================================================

# (x, y, horizontal alignment, leader line) in data coordinates.
ALN_LABELS_A = {
    ("VACmap", "PacBio"): (0.68, 98.75, "center", False),
    ("VG Giraffe", "PacBio"): (0.95, 98.2, "center", False),
    ("pbmm2", "PacBio"): (1.55, 99.15, "left", True),
    ("minimap2", "PacBio"): (1.75, 99.75, "left", True),
    ("VACmap", "ONT"): (1.59, 97.55, "center", False),
    ("VG Giraffe", "ONT"): (2.45, 92.75, "left", False),
    ("pbmm2", "ONT"): (2.64, 98.15, "center", False),
    ("minimap2", "ONT"): (3.05, 99.88, "center", False),
}


def panel_a(axis, aln):
    for (aligner, tech), subset in aln.groupby(["aligner", "read_technology"]):
        for _, row in subset.iterrows():
            draw_point(axis, row["error_percent"], row["input_normalized_cigar_yield_percent"],
                       ALIGNER_MARKERS[aligner], row["sample"], tech == "ONT")
        label_at(axis, subset, "error_percent", "input_normalized_cigar_yield_percent",
                 f"{aligner} {TECH_SHORT[tech]}", ALN_LABELS_A[(aligner, tech)])
    axis.set_xlim(0, None)
    axis.set_ylim(90, 100)
    axis.set_xlabel("Error rate (%)")
    axis.set_ylabel("CIGAR-aligned yield (%)\n(of input sequencing bases)")
    axis.set_title("Alignment yield vs accuracy", fontsize=7, pad=4)
    axis.text(0.98, 0.04, "y-axis starts at 90%", transform=axis.transAxes, ha="right",
              va="bottom", fontsize=5.5, color="#777777", style="italic")


def panel_b(axis, aln):
    measured = aln[aln["peak_ram_gb"].notna()]
    unmeasured = sorted(aln.loc[aln["peak_ram_gb"].isna(), "aligner"].unique())
    positions = {("pbmm2", "ONT"): (118, 68, "left", True),
                 ("pbmm2", "PacBio"): (80, 55, "right", True),
                 ("minimap2", "PacBio"): (99, 30, "center", False),
                 ("minimap2", "ONT"): (158, 30, "center", False)}
    for (aligner, tech), subset in measured.groupby(["aligner", "read_technology"]):
        for _, row in subset.iterrows():
            draw_point(axis, row["thread_hours"], row["peak_ram_gb"],
                       ALIGNER_MARKERS[aligner], row["sample"], tech == "ONT")
        label_at(axis, subset, "thread_hours", "peak_ram_gb",
                 f"{aligner} {TECH_SHORT[tech]}", positions[(aligner, tech)])
    axis.set_xlim(0, 180)
    axis.set_ylim(0, 75)
    axis.set_xlabel("Allocated thread-hours (wall-clock h × threads)")
    axis.set_ylabel("Peak RAM (GB)")
    axis.set_title("Alignment resource use", fontsize=7, pad=4)
    if unmeasured:
        axis.text(0.98, 0.04, "Peak RAM not measured: " + ", ".join(unmeasured),
                  transform=axis.transAxes, ha="right", va="bottom", fontsize=5.5,
                  color="#777777", style="italic")


ASM_LABEL_OFFSETS = {
    ("flye", "ont"): (0, -11), ("flye", "pb"): (0, 10), ("goldrush", "ont"): (0, 10),
    ("goldrush", "pb"): (0, -10), ("verkko", "hybrid"): (0, -11),
}


def panel_c(axis, asm):
    for (assembler, tech) in ASM_CONFIGS:
        subset = asm[(asm["assembler"] == assembler) & (asm["technology"] == tech)]
        for _, row in subset.iterrows():
            draw_point(axis, row["misassemblies"], row["nga50_mb"], ASM_MARKERS[assembler],
                       row["sample"], tech != "pb")
        name = f"{ASM_NAMES[assembler]} {TECH_SHORT[tech]}" + ("$^\\dagger$" if assembler == "verkko" else "")
        label_cluster(axis, subset, "misassemblies", "nga50_mb", name,
                      ASM_LABEL_OFFSETS[(assembler, tech)], log_y=True)
    axis.set_yscale("log")
    axis.set_ylim(0.03, 60)
    axis.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:g}"))
    axis.set_xlim(0, 20000)
    axis.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    axis.set_xlabel("Misassemblies (QUAST, extensive)")
    axis.set_ylabel("NGA50 (Mb, log scale)")
    axis.set_title("Assembly contiguity vs structural errors", fontsize=7, pad=4)


def panel_d(axis, perf):
    axis.set_title("Assembly resource use", fontsize=7, pad=4)
    axis.set_xlabel("Allocated thread-hours (wall-clock h × threads)")
    axis.set_ylabel("Peak RAM (GB)")
    if perf is None:
        axis.set_xticks([])
        axis.set_yticks([])
        axis.text(0.5, 0.5, "Not yet measured\n(assembler benchmark runs pending)",
                  transform=axis.transAxes, ha="center", va="center", fontsize=6.5,
                  color="#777777", style="italic")
        return
    for _, row in perf.iterrows():
        draw_point(axis, row["thread_hours"], row["peak_rss_gb"], ASM_MARKERS[row["assembler"]],
                   row["sample"], row["technology"] != "pb")
    axis.set_xlim(0, None)
    axis.set_ylim(0, None)


def _marker(label, marker, filled=True):
    return Line2D([], [], linestyle="", marker=marker, markersize=4.5,
                  markerfacecolor="#7F7F7F" if filled else "white",
                  markeredgecolor="white" if filled else "#7F7F7F",
                  markeredgewidth=0.5 if filled else 1.1, label=label)


def _heading(label):
    return Line2D([], [], linestyle="", marker="", label=label)


def legend_rows():
    """Three legend rows so shared shapes (circle, square) are never ambiguous."""
    samples = [_heading("$\\bf{Sample}$")] + [
        Line2D([], [], linestyle="", marker="o", markersize=4.5, markerfacecolor=SAMPLE_COLORS[s],
               markeredgecolor="white", label=s) for s in SAMPLE_ORDER
    ] + [_heading("   $\\bf{Input}$"), _marker("ONT (filled)", "o"),
         _marker("HiFi (open)", "o", filled=False)]
    aligners = [_heading("$\\bf{Aligners\\ (a,\\ b)}$")] + [
        _marker(a, m) for a, m in ALIGNER_MARKERS.items()]
    assemblers = [_heading("$\\bf{Assemblers\\ (c,\\ d)}$")] + [
        _marker(ASM_NAMES[a], ASM_MARKERS[a]) for a in ["flye", "goldrush", "verkko"]]
    return [samples, aligners, assemblers]


def tidy(aln, asm, perf) -> pd.DataFrame:
    rows = []
    for _, r in aln.iterrows():
        cfg = f"{r['aligner']} ({TECH_SHORT[r['read_technology']]})"
        rows += [("a", "Alignment-based", cfg, r["sample"], "error_percent", r["error_percent"]),
                 ("a", "Alignment-based", cfg, r["sample"], "input_normalized_cigar_yield_percent",
                  r["input_normalized_cigar_yield_percent"])]
        if pd.notna(r["peak_ram_gb"]):
            rows += [("b", "Alignment-based", cfg, r["sample"], "thread_hours", r["thread_hours"]),
                     ("b", "Alignment-based", cfg, r["sample"], "peak_ram_gb", r["peak_ram_gb"])]
    for _, r in asm.iterrows():
        cfg = f"{ASM_NAMES[r['assembler']]} ({TECH_SHORT[r['technology']]})"
        rows += [("c", "Assembly-based", cfg, r["sample"], "misassemblies", r["misassemblies"]),
                 ("c", "Assembly-based", cfg, r["sample"], "nga50_mb", r["nga50_mb"])]
    if perf is not None:
        for _, r in perf.iterrows():
            cfg = f"{ASM_NAMES[r['assembler']]} ({TECH_SHORT[r['technology']]})"
            rows += [("d", "Assembly-based", cfg, r["sample"], "thread_hours", r["thread_hours"]),
                     ("d", "Assembly-based", cfg, r["sample"], "peak_rss_gb", r["peak_rss_gb"])]
    return pd.DataFrame(rows, columns=["panel", "strategy", "configuration", "sample", "metric", "value"])


def main() -> int:
    apply_style()
    aln, asm, perf = load_alignment(), load_assembly(), load_assembly_performance()

    figure, axes = plt.subplots(2, 2, figsize=(FULL_WIDTH_IN, 6.0))
    figure.subplots_adjust(left=0.10, right=0.98, top=0.93, bottom=0.20, hspace=0.55, wspace=0.32)
    panel_a(axes[0, 0], aln)
    panel_b(axes[0, 1], aln)
    panel_c(axes[1, 0], asm)
    panel_d(axes[1, 1], perf)
    for axis, letter in zip(axes.flat, "abcd"):
        clean_spines(axis)
        subtle_grid(axis, "both")
        panel_letter(axis, letter, x=-0.22, y=1.03)
    for row, label in [(0, "Alignment-based"), (1, "Assembly-based")]:
        box = axes[row, 0].get_position()
        figure.text(0.01, (box.y0 + box.y1) / 2, label, rotation=90, ha="left", va="center",
                    fontsize=8, fontweight="bold")
    for row_index, handles in enumerate(legend_rows()):
        figure.legend(handles=handles, loc="lower center",
                      bbox_to_anchor=(0.5, 0.075 - 0.032 * row_index), ncols=len(handles),
                      frameon=False, fontsize=6.2, handletextpad=0.3, columnspacing=1.1,
                      handlelength=1.0)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    tidy(aln, asm, perf).to_csv(OUT_DIR / f"{STEM}_panel_data.tsv", sep="\t", index=False,
                                float_format="%.6g")
    for suffix in (".pdf", ".svg", ".png"):
        figure.savefig(OUT_DIR / f"{STEM}{suffix}", dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(figure)
    print(f"Figure: {(OUT_DIR / STEM).relative_to(PROJECT)}.{{pdf,svg,png}}")
    print("Panel d:", "plotted" if perf is not None else "not measured (placeholder)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
