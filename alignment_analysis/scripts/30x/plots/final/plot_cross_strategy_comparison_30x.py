#!/usr/bin/env python3
"""Cross-strategy comparison of alignment- and assembly-based long-read analysis.

Figure: "Cross-strategy comparison of alignment- and assembly-based long-read
genome analysis" (30x, GIAB HG002/HG003/HG004, ONT and PacBio HiFi).

Every value is read from the curated benchmark tables; nothing is estimated,
interpolated or imputed. Where an aligner and an assembler metric do not share
a definition, they are drawn in adjacent, separately-labelled subpanels and
never merged onto one axis under one name.

Column mapping (verified against the source tables, see the warnings file
for why each pair is or is not mergeable):

Panel A  substitutions
  alignment : (samtools stats `mismatches` - `inserted_bases` - `deleted_bases`)
              / `bases_mapped_cigar` x 1e5. samtools takes `mismatches` from
              the NM tag (edit distance = substitutions + inserted + deleted
              bases); the indel bases (length x count summed over the samtools
              stats ID records, derived/alignment_benchmark_30x_indel_recovery.tsv)
              are subtracted to leave substitutions. ID records stop at 300 bp,
              so bases of longer indels stay in the result: it is an upper
              bound on the substitution rate.
  assembly  : QUAST `mismatches_per_100kbp` (substitutions only, per 100 kbp
              of `Total aligned length`) (assembly_benchmark_30x.tsv).
  -> both substitution rates, but raw reads vs consensus contigs:
     method-specific subpanels, independent y-axes.

Panel B  indel disagreement
  alignment : (`insertion_events` + `deletion_events`) / `bases_mapped_cigar`
              x 1e5 (derived/alignment_benchmark_30x_indel_recovery.tsv).
              Numerator = number of CIGAR I/D OPERATIONS (events, not bases):
              quality_check_aligners_indels.py sums the per-length counts of
              the samtools stats "ID" records, which cover 1-300 bp; samtools
              skips longer indels (stats.c: `if (ncig <= nindels)`).
  assembly  : QUAST `indels_per_100kbp` (indel events, not bases, per
              100 kbp of aligned contig bases).
  -> both are event counts per 100 kbp of aligned query bases, but the
     query is raw reads vs consensus contigs: separate subpanels with
     independent y-axes, no merged series and no statistical test.

Panel C  alignment yield and assembly completeness
  alignment : `input_normalized_cigar_yield_percent`
              (derived/plot_data/input_normalized_cigar_yield_percent.tsv):
              bases_mapped_cigar / largest samtools `total_length` (primary
              read bases) among the four aligners for that sample+technology.
  assembly  : QUAST `genome_fraction_pct` (denominator = GRCh38).
  -> distinct metrics and denominators, distinct subpanels.

Alignment filters (samtools stats, verified in the .smk rules and the stats
headers): no -F/-f/-q options, so secondary alignments are skipped (samtools
default), supplementary alignments and MAPQ 0 records are counted, unmapped
reads are excluded. `bases_mapped_cigar` = M/I/=/X read bases (soft clips
excluded). QUAST 5.3.0 --large --min-contig 500: min alignment 500 bp, min
identity 95 %, ambiguity "one" (quast.log).

Verkko (2.3.2, --hifi and --nano only, no trio or Hi-C data) was evaluated as
its single combined assembly.fasta (6.08-6.14 Gb, QUAST duplication ratio
2.04-2.09 vs haploid GRCh38; no per-haplotype FASTAs were produced). Phasing
status is not stated: it cannot be read from the evaluated FASTA alone.

Panel D  computational cost (allocated thread-hours)
  alignment : `wallclock_runtime_hours` x `threads` (final table). Measured
              CPU time exists only for minimap2/pbmm2, so CPU-hours are not
              used for anyone.
  assembly  : `wall_clock_hours` x `threads` from
              assembly_analysis/tables/assembler_performance_30x.tsv (written
              by parse_assembler_time_v.py). If that table is absent or
              incomplete, Panel D is NOT drawn and a warning is emitted.
  Supplementary RAM panel only if peak RAM is complete for every tool.
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
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
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
from matplotlib import transforms  # noqa: E402
from matplotlib.ticker import MaxNLocator, StrMethodFormatter  # noqa: E402
from utils.plot_style import (  # noqa: E402
    SAMPLE_COLORS,
    SAMPLE_ORDER,
    apply_style,
    clean_spines,
    panel_letter,
    rotated_xticks,
    subtle_grid,
)

ALN_TABLES = PROJECT / "alignment_analysis" / "tables" / "30x"
ALN_FINAL_TSV = ALN_TABLES / "final" / "alignment_benchmark_30x.tsv"
ALN_INDEL_TSV = ALN_TABLES / "derived" / "alignment_benchmark_30x_indel_recovery.tsv"
ALN_YIELD_TSV = ALN_TABLES / "derived" / "plot_data" / "input_normalized_cigar_yield_percent.tsv"
ASM_TSV = PROJECT / "assembly_analysis" / "tables" / "30x" / "final" / "assembly_benchmark_30x.tsv"
ASM_PERF_TSV = PROJECT / "assembly_analysis" / "tables" / "assembler_performance_30x.tsv"

OUT_DIR = PROJECT / "figures" / "30x" / "cross_strategy"
PANEL_DIR = OUT_DIR / "panels"
STEM = "cross_strategy_comparison_30x"
TIDY_TSV = OUT_DIR / f"{STEM}_panel_data.tsv"
WARNINGS_TXT = OUT_DIR / f"{STEM}_warnings.txt"

ALIGNERS = ["minimap2", "pbmm2", "VACmap", "VG Giraffe"]
ALN_TECHNOLOGIES = ["ONT", "PacBio"]
# Grouped by sequencing input (ONT, HiFi, then hybrid), as for the aligners,
# so one input bracket under the axis labels each group.
ASM_CONFIGS = [
    ("flye", "ont"),
    ("goldrush", "ont"),
    ("flye", "pb"),
    ("goldrush", "pb"),
    ("verkko", "hybrid"),
]
ASM_NAMES = {"flye": "Flye", "goldrush": "GoldRush", "verkko": "Verkko"}
TECH_LABEL = {"ONT": "ONT", "PacBio": "HiFi", "ont": "ONT", "pb": "HiFi", "hybrid": "ONT+HiFi"}

# Technology is encoded by marker shape; sample by color (shared palette).
TECH_MARKER = {"ONT": "o", "HiFi": "s", "ONT+HiFi": "D"}
TECH_MARKER_AREA = {"o": 18, "s": 15, "D": 13}  # equal visual weight per shape
TECH_BRACKET_LABEL = {"ONT": "ONT", "HiFi": "PacBio HiFi", "ONT+HiFi": "ONT +\nHiFi"}
SAMPLE_OFFSETS = dict(zip(SAMPLE_ORDER, [-0.2, 0.0, 0.2]))
TICK_NAMES = {"Verkko": "Verkko\u2020"}

# Drawn at the report's 6.5 in text width and included at \linewidth (1:1),
# so these are the printed sizes.
FIGURE_WIDTH_IN = 6.5
TICK_SIZE = 8.0
AXIS_LABEL_SIZE = 7.5
SUBTITLE_SIZE = 7.5
ROW_TITLE_SIZE = 9.0
HEADER_SIZE = 9.5
LEGEND_SIZE = 7.5
BRACKET_SIZE = 7.5
TICK_ROTATION = 35
BRACKET_DROP_PT = 34  # below the axis line, clears the rotated names
GROUP_GAP = 0.6  # extra x-axis space between sequencing-input groups

# Mean +/- SD: a subdued descriptive overlay behind the observations.
SUMMARY_COLOR = "#9A9A9A"

# Verkko: very light hatching, so the observations stay dominant.
VERKKO_FACE = "#FBFBFB"
VERKKO_HATCH_EDGE = "#E6E6E6"
VERKKO_LEGEND = "Verkko diploid assembly (hatched)"

STRATEGY_ALN = "Alignment-based"
STRATEGY_ASM = "Assembly-based"

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
    # Reindex so an absent row becomes NaN (shown as "n/a"), never a value.
    return data.set_index(keys).reindex(expected).reset_index()


def load_alignment() -> pd.DataFrame:
    final = normalise_aligner_keys(read_tsv(
        ALN_FINAL_TSV,
        {"sample", "read_technology", "aligner", "mismatches", "bases_mapped_cigar",
         "error_rate", "threads", "wallclock_runtime_hours", "cpu_time_seconds", "peak_ram_gb"},
    ))
    indel = normalise_aligner_keys(read_tsv(
        ALN_INDEL_TSV,
        {"sample", "read_technology", "aligner", "insertion_events", "deletion_events",
         "insertion_events_per_100kb", "deletion_events_per_100kb", "inserted_bases",
         "deleted_bases", "mismatches", "bases_mapped_cigar"},
    ))
    yield_ = normalise_aligner_keys(read_tsv(
        ALN_YIELD_TSV,
        {"sample", "read_technology", "aligner", "common_input_bases", "bases_mapped_cigar",
         "input_normalized_cigar_yield_percent"},
    ))

    keys = ["sample", "read_technology", "aligner"]
    data = final[keys + ["mismatches", "bases_mapped_cigar", "error_rate", "threads",
                         "wallclock_runtime_hours", "cpu_time_seconds", "peak_ram_gb"]]
    data = data.merge(
        indel[keys + ["insertion_events", "deletion_events", "insertion_events_per_100kb",
                      "deletion_events_per_100kb", "inserted_bases", "deleted_bases",
                      "mismatches", "bases_mapped_cigar"]]
        .rename(columns={"bases_mapped_cigar": "bases_mapped_cigar_indel_table",
                         "mismatches": "mismatches_indel_table"}),
        on=keys, validate="one_to_one",
    )
    data = data.merge(
        yield_[keys + ["common_input_bases", "input_normalized_cigar_yield_percent",
                       "bases_mapped_cigar"]]
        .rename(columns={"bases_mapped_cigar": "bases_mapped_cigar_yield_table"}),
        on=keys, validate="one_to_one",
    )

    # Consistency checks between the three source tables: same denominator,
    # and the stored per-100kb rates must reproduce from the raw counts.
    for column, other in [("bases_mapped_cigar", "bases_mapped_cigar_indel_table"),
                          ("bases_mapped_cigar", "bases_mapped_cigar_yield_table"),
                          ("mismatches", "mismatches_indel_table")]:
        both = data[[column, other]].dropna()
        if not np.array_equal(both.iloc[:, 0].to_numpy(), both.iloc[:, 1].to_numpy()):
            raise ValueError(f"{column} disagrees between source tables ({other})")
    substitutions = data["mismatches"] - data["inserted_bases"] - data["deleted_bases"]
    if (substitutions < 0).any():
        raise ValueError("NM minus indel bases is negative: indel bases exceed NM")
    rate = data["mismatches"] / data["bases_mapped_cigar"]
    if not np.allclose(rate, data["error_rate"], rtol=1e-3, equal_nan=True):
        raise ValueError("mismatches / bases_mapped_cigar does not reproduce samtools error_rate")
    for event in ["insertion", "deletion"]:
        recomputed = data[f"{event}_events"] / data["bases_mapped_cigar"] * 1e5
        if not np.allclose(recomputed, data[f"{event}_events_per_100kb"], rtol=1e-6, equal_nan=True):
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
        {"assembler", "sample", "technology", "representation", "genome_fraction_pct",
         "mismatches_per_100kbp", "indels_per_100kbp", "source_file"},
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
    return data


def load_assembly_performance() -> pd.DataFrame | None:
    """Return per-sample assembler thread-hours/RAM, or None if not complete."""
    if not ASM_PERF_TSV.is_file():
        warn(
            "Panel D NOT generated: assembler computational-resource table is absent "
            f"({ASM_PERF_TSV.relative_to(PROJECT)}). Run the *.benchmark.smk assembler "
            "workflows (parse_assembler_time_v.py) to produce it. No assembler runtime "
            "is estimated or substituted."
        )
        return None
    perf = read_tsv(ASM_PERF_TSV, {"assembler", "sample", "technology", "threads",
                                   "wall_clock_hours", "peak_rss_gb"})
    keys = ["assembler", "sample", "technology"]
    if perf.duplicated(keys).any():
        raise ValueError("Duplicate rows in assembler performance table")
    expected = pd.MultiIndex.from_tuples(
        [(asm, sample, tech) for asm, tech in ASM_CONFIGS for sample in SAMPLE_ORDER], names=keys
    )
    perf = perf.set_index(keys).reindex(expected).reset_index()
    missing = perf[perf[["threads", "wall_clock_hours"]].isna().any(axis=1)]
    if len(missing):
        warn(
            "Panel D NOT generated: assembler performance table is incomplete for "
            + ", ".join(missing[keys].astype(str).agg(".".join, axis=1))
        )
        return None
    perf["thread_hours"] = perf["wall_clock_hours"] * perf["threads"]
    perf["strategy"] = STRATEGY_ASM
    perf["method"] = perf["assembler"].map(ASM_NAMES)
    perf["technology"] = perf["technology"].map(TECH_LABEL)
    perf["configuration"] = perf["method"] + " (" + perf["technology"] + ")"
    return perf


# ============================================================
# TIDY PANEL DATA
# ============================================================

TIDY_COLUMNS = [
    "panel", "subpanel", "strategy", "method", "technology", "configuration", "sample",
    "representation", "metric", "value", "unit", "numerator", "denominator",
    "source_file", "source_columns",
]


def tidy(frame: pd.DataFrame, panel: str, subpanel: str, value_column: str, metric: str,
         unit: str, numerator: str, denominator: str, source: Path,
         source_columns: str) -> pd.DataFrame:
    out = frame[["strategy", "method", "technology", "configuration", "sample",
                 "representation"]].copy()
    out["value"] = frame[value_column].to_numpy(dtype=float)
    out["panel"] = panel
    out["subpanel"] = subpanel
    out["metric"] = metric
    out["unit"] = unit
    out["numerator"] = numerator
    out["denominator"] = denominator
    out["source_file"] = str(source.relative_to(PROJECT))
    out["source_columns"] = source_columns
    return out[TIDY_COLUMNS]


def build_panel_data(aln: pd.DataFrame, asm: pd.DataFrame,
                     perf: pd.DataFrame | None) -> pd.DataFrame:
    aln = aln.copy()
    aln["substitutions_per_100kbp"] = (
        (aln["mismatches"] - aln["inserted_bases"] - aln["deleted_bases"])
        / aln["bases_mapped_cigar"] * 1e5
    )
    aln["indel_events_per_100kbp"] = (
        (aln["insertion_events"] + aln["deletion_events"]) / aln["bases_mapped_cigar"] * 1e5
    )
    aln["thread_hours"] = aln["wallclock_runtime_hours"] * aln["threads"]

    asm_source = ASM_TSV
    parts = [
        tidy(aln, "A", "A1", "substitutions_per_100kbp",
             "Substitutions per 100 kbp (NM minus 1-300 bp indel bases; upper bound)",
             "per 100 kbp",
             "samtools stats mismatches (NM) - inserted_bases - deleted_bases (ID records, 1-300 bp)",
             "CIGAR-mapped read bases (bases_mapped_cigar)", ALN_INDEL_TSV,
             "mismatches; inserted_bases; deleted_bases; bases_mapped_cigar"),
        tidy(asm, "A", "A2", "mismatches_per_100kbp", "QUAST mismatches per 100 kbp",
             "per 100 kbp", "QUAST # mismatches (substitutions only)",
             "QUAST total aligned contig length", asm_source, "mismatches_per_100kbp"),
        tidy(aln, "B", "B1", "indel_events_per_100kbp", "Indel events per 100 kbp",
             "per 100 kbp", "CIGAR I + D operations (samtools stats ID records)",
             "CIGAR-mapped read bases (bases_mapped_cigar)", ALN_INDEL_TSV,
             "insertion_events; deletion_events; bases_mapped_cigar"),
        tidy(asm, "B", "B2", "indels_per_100kbp", "Indel events per 100 kbp", "per 100 kbp",
             "QUAST # indels (events; consecutive indel bases counted once)",
             "QUAST total aligned contig length", asm_source, "indels_per_100kbp"),
        tidy(aln, "C", "C1", "input_normalized_cigar_yield_percent", "CIGAR-aligned yield",
             "%", "CIGAR-mapped read bases (bases_mapped_cigar)",
             "input sequencing bases (common_input_bases)", ALN_YIELD_TSV,
             "input_normalized_cigar_yield_percent"),
        tidy(asm, "C", "C2", "genome_fraction_pct", "Genome fraction", "%",
             "GRCh38 bases covered by aligned contigs", "GRCh38 reference length",
             asm_source, "genome_fraction_pct"),
    ]
    if perf is not None:
        parts += [
            tidy(aln, "D", "D1", "thread_hours", "Allocated thread-hours", "thread-h",
                 "wall-clock hours x allocated threads", "per sample", ALN_FINAL_TSV,
                 "wallclock_runtime_hours; threads"),
            tidy(perf.assign(representation=perf["assembler"].map(
                     {"flye": "collapsed", "goldrush": "collapsed",
                      "verkko": "haplotype_resolved"})),
                 "D", "D2", "thread_hours", "Allocated thread-hours", "thread-h",
                 "wall-clock hours x allocated threads", "per sample", ASM_PERF_TSV,
                 "wall_clock_hours; threads"),
        ]
    panel_data = pd.concat(parts, ignore_index=True)

    missing = panel_data[panel_data["value"].isna()]
    for _, row in missing.iterrows():
        warn(f"Panel {row['subpanel']}: no value for {row['configuration']} {row['sample']} "
             "(left blank, marked n/a)")
    return panel_data


# ============================================================
# PLOTTING
# ============================================================


def split_config(config: str) -> tuple[str, str]:
    method, tech = config.rsplit(" (", 1)
    return method, tech.rstrip(")")


def draw_strip(axis: plt.Axes, data: pd.DataFrame, configs: list[str],
               shade_configs: tuple[str, ...] = ()) -> None:
    """Raw per-sample points over a subdued mean +/- SD per configuration."""
    x_positions = grouped_positions(configs)
    for x, config in zip(x_positions, configs):
        subset = data[data["configuration"] == config]
        if config in shade_configs:
            with plt.rc_context({"hatch.linewidth": 0.4}):
                axis.axvspan(x - 0.45, x + 0.45, facecolor=VERKKO_FACE,
                             edgecolor=VERKKO_HATCH_EDGE, hatch="////", linewidth=0, zorder=0)
        values = subset["value"].to_numpy(dtype=float)
        finite = values[np.isfinite(values)]
        if len(finite) >= 2:
            mean, sd = finite.mean(), finite.std(ddof=1)
            axis.errorbar(x, mean, yerr=sd, fmt="none", ecolor=SUMMARY_COLOR,
                          elinewidth=0.8, capsize=2.0, capthick=0.8, zorder=2)
            axis.hlines(mean, x - 0.28, x + 0.28, color=SUMMARY_COLOR, linewidth=1.0, zorder=2)
        for _, row in subset.iterrows():
            offset = SAMPLE_OFFSETS[row["sample"]]
            if np.isfinite(row["value"]):
                marker = TECH_MARKER[row["technology"]]
                axis.scatter(x + offset, row["value"], s=TECH_MARKER_AREA[marker], marker=marker,
                             color=SAMPLE_COLORS[row["sample"]], edgecolor="white",
                             linewidth=0.4, zorder=3, clip_on=False)
            else:
                axis.text(x + offset, 0, "n/a", rotation=90, fontsize=5, color="#999999",
                          ha="center", va="bottom", zorder=3)
    labels = [TICK_NAMES.get(split_config(c)[0], split_config(c)[0]) for c in configs]
    rotated_xticks(axis, x_positions, labels, fontsize=TICK_SIZE, rotation=TICK_ROTATION,
                   tick_length=2)
    axis.tick_params(axis="y", labelsize=TICK_SIZE, length=2, pad=1.5)
    axis.set_xlim(-0.6, x_positions[-1] + 0.6)
    clean_spines(axis)
    subtle_grid(axis)
    add_input_brackets(axis, configs, x_positions)


def grouped_positions(configs: list[str]) -> np.ndarray:
    """Unit spacing within a sequencing input, plus GROUP_GAP between inputs."""
    techs = [split_config(c)[1] for c in configs]
    gaps = [0.0] + [GROUP_GAP if techs[i] != techs[i - 1] else 0.0 for i in range(1, len(techs))]
    return np.arange(len(configs), dtype=float) + np.cumsum(gaps)


def add_input_brackets(axis: plt.Axes, configs: list[str], x_positions: np.ndarray) -> None:
    """One bracket per sequencing input below the method names."""
    figure = axis.get_figure()
    below = transforms.ScaledTranslation(0, -BRACKET_DROP_PT / 72, figure.dpi_scale_trans)
    transform = transforms.blended_transform_factory(axis.transData, axis.transAxes) + below
    techs = [split_config(c)[1] for c in configs]
    start = 0
    for index in range(1, len(techs) + 1):
        if index < len(techs) and techs[index] == techs[start]:
            continue
        x0, x1 = x_positions[start] - 0.38, x_positions[index - 1] + 0.38
        axis.plot([x0, x0, x1, x1], [0.03, 0, 0, 0.03], transform=transform,
                  color="#555555", linewidth=0.6, clip_on=False)
        axis.text((x0 + x1) / 2, -0.03, TECH_BRACKET_LABEL[techs[start]], transform=transform,
                  ha="center", va="top", fontsize=BRACKET_SIZE, color="#333333",
                  linespacing=0.95)
        start = index


def aln_configs() -> list[str]:
    return [f"{a} ({TECH_LABEL[t]})" for t in ALN_TECHNOLOGIES for a in ALIGNERS]


def asm_configs() -> list[str]:
    return [f"{ASM_NAMES[a]} ({TECH_LABEL[t]})" for a, t in ASM_CONFIGS]


VERKKO_CONFIG = f"{ASM_NAMES['verkko']} ({TECH_LABEL['hybrid']})"

PANEL_SPECS = {
    "A": dict(
        title="Substitutions",
        left_ylabel="Substituted bases\nper 100 kbp mapped",
        left_title="Reads: substitutions\nper 100 kbp CIGAR-mapped read bases",
        right_ylabel="QUAST mismatches\nper 100 kbp aligned",
        right_title="Contigs: substitutions only\nper 100 kbp aligned contig bases",
        percent=False,
    ),
    "B": dict(
        title="Indel disagreement (events, not bases)",
        left_ylabel="CIGAR indel events\nper 100 kbp mapped",
        left_title="Reads: CIGAR I + D operations, 1\u2013300 bp\nper 100 kbp CIGAR-mapped read bases",
        right_ylabel="QUAST indel events\nper 100 kbp aligned",
        right_title="Contigs: QUAST indel events\nper 100 kbp aligned contig bases",
        percent=False,
    ),
    "C": dict(
        title="Alignment yield and assembly completeness",
        left_ylabel="CIGAR-aligned yield\n(% of input bases)",
        left_title="Reads: CIGAR-mapped bases\n% of input sequencing bases",
        right_ylabel="Genome fraction\n(% of GRCh38)",
        right_title="Contigs: QUAST genome fraction\n% of GRCh38 covered",
        percent=True,
    ),
    "D": dict(
        title="Computational cost",
        left_ylabel="Allocated thread-hours",
        left_title="Alignment (wall-clock h x threads)",
        right_ylabel="Allocated thread-hours",
        right_title="Assembly (wall-clock h x threads)",
        percent=False,
    ),
}


def set_y_axis(axis: plt.Axes, data: pd.DataFrame, percent: bool) -> None:
    """Independent linear axis from zero to just above the highest point or error bar."""
    values = data["value"].to_numpy(dtype=float)
    stats = data.groupby("configuration")["value"].agg(["mean", "std"])
    highest = float(np.nanmax(np.r_[values, (stats["mean"] + stats["std"]).to_numpy()]))
    if percent:
        axis.set_ylim(0, 103)
        axis.set_yticks([0, 20, 40, 60, 80, 100])
        return
    ticks = [t for t in MaxNLocator(nbins=4, steps=[1, 2, 2.5, 5, 10]).tick_values(0, highest)
             if t >= 0]
    step = ticks[1] - ticks[0]
    while ticks[-1] < highest * 1.02:
        ticks.append(ticks[-1] + step)
    axis.set_ylim(0, ticks[-1])
    axis.set_yticks(ticks)
    if ticks[-1] >= 1000:
        axis.yaxis.set_major_formatter(StrMethodFormatter("{x:,.0f}"))


def draw_panel_pair(left: plt.Axes, right: plt.Axes, panel_data: pd.DataFrame, key: str,
                    letter: bool = True) -> None:
    spec = PANEL_SPECS[key]
    left_data = panel_data[(panel_data["panel"] == key) & (panel_data["strategy"] == STRATEGY_ALN)]
    right_data = panel_data[(panel_data["panel"] == key) & (panel_data["strategy"] == STRATEGY_ASM)]
    draw_strip(left, left_data, aln_configs())
    draw_strip(right, right_data, asm_configs(), shade_configs=(VERKKO_CONFIG,))

    for axis, data, label in [(left, left_data, spec["left_ylabel"]),
                              (right, right_data, spec["right_ylabel"])]:
        set_y_axis(axis, data, spec["percent"])
        axis.set_ylabel(label, fontsize=AXIS_LABEL_SIZE, labelpad=3)
    left.set_title(spec["left_title"], fontsize=SUBTITLE_SIZE, pad=4, linespacing=1.15)
    right.set_title(spec["right_title"], fontsize=SUBTITLE_SIZE, pad=4, linespacing=1.15)
    if letter:
        # Fixed offsets in points above the axes (clear of the two-line
        # subtitles) and left of the y-label, whatever the row height.
        figure = left.get_figure()
        for dx, text, size in [(-46, key.lower(), ROW_TITLE_SIZE + 1),
                               (-34, spec["title"], ROW_TITLE_SIZE)]:
            shift = transforms.ScaledTranslation(dx / 72, 27 / 72, figure.dpi_scale_trans)
            left.text(0, 1, text, transform=left.transAxes + shift, fontsize=size,
                      fontweight="bold", ha="left", va="bottom")


def legend_handles() -> list:
    handles = [Line2D([], [], linestyle="", marker="o", markersize=5,
                      markerfacecolor=SAMPLE_COLORS[s], markeredgecolor="white", label=s)
               for s in SAMPLE_ORDER]
    handles += [Line2D([], [], linestyle="", marker=TECH_MARKER[t],
                       markersize=TECH_MARKER_AREA[TECH_MARKER[t]] ** 0.5 * 1.15,
                       markerfacecolor="#7F7F7F", markeredgecolor="white",
                       label=TECH_BRACKET_LABEL[t].replace("\n", " "))
                for t in ["ONT", "HiFi", "ONT+HiFi"]]
    handles += [
        Line2D([], [], color=SUMMARY_COLOR, linewidth=1.0, marker="|", markersize=6,
               markeredgewidth=0.8, label="Mean \u00b1 SD (n = 3 samples)"),
        Patch(facecolor=VERKKO_FACE, edgecolor="#BDBDBD", hatch="////", linewidth=0.4,
              label=VERKKO_LEGEND),
    ]
    # figure.legend fills column-first: interleave so row 1 = samples + mean,
    # row 2 = sequencing inputs + Verkko note.
    samples, inputs, (mean_sd, verkko) = handles[:3], handles[3:6], handles[6:]
    ordered = []
    for top, bottom in zip(samples + [mean_sd], inputs + [verkko]):
        ordered += [top, bottom]
    return ordered


def add_strategy_headers(figure: plt.Figure, left: plt.Axes, right: plt.Axes) -> None:
    for axis, label in [(left, "Read aligners"), (right, "Genome assemblers")]:
        box = axis.get_position()
        figure.text((box.x0 + box.x1) / 2, HEADER_Y + 0.004, label, ha="center", va="bottom",
                    fontsize=HEADER_SIZE, fontweight="bold")
        figure.add_artist(Line2D([box.x0, box.x1], [HEADER_Y, HEADER_Y],
                                 transform=figure.transFigure, color="black", linewidth=0.8))


def save(figure: plt.Figure, stem: Path) -> list[Path]:
    stem.parent.mkdir(parents=True, exist_ok=True)
    paths = [stem.with_suffix(s) for s in (".pdf", ".svg", ".png")]
    for path in paths:
        figure.savefig(path, dpi=600, bbox_inches="tight", facecolor="white")
    return paths


HEADER_Y = 0.962


def plot_main(panel_data: pd.DataFrame, panels: list[str]) -> list[Path]:
    n_rows = len(panels)
    figure = plt.figure(figsize=(FIGURE_WIDTH_IN, 1.94 * n_rows + 0.72))
    grid = figure.add_gridspec(n_rows, 2, width_ratios=[8, 5], wspace=0.30, hspace=1.70,
                               left=0.112, right=0.995, top=0.865, bottom=0.17)
    first_pair = None
    for row, key in enumerate(panels):
        left = figure.add_subplot(grid[row, 0])
        right = figure.add_subplot(grid[row, 1])
        draw_panel_pair(left, right, panel_data, key)
        first_pair = first_pair or (left, right)
    add_strategy_headers(figure, *first_pair)
    # No in-figure title: the title belongs in the manuscript caption.
    figure.legend(handles=legend_handles(), loc="lower center", bbox_to_anchor=(0.5, 0.0),
                  ncols=4, frameon=False, fontsize=LEGEND_SIZE, handletextpad=0.4,
                  columnspacing=1.4)
    paths = save(figure, OUT_DIR / STEM)
    plt.close(figure)
    return paths


def plot_single_panels(panel_data: pd.DataFrame, panels: list[str]) -> list[Path]:
    paths = []
    for key in panels:
        figure = plt.figure(figsize=(FIGURE_WIDTH_IN, 3.0))
        grid = figure.add_gridspec(1, 2, width_ratios=[8, 5], wspace=0.28,
                                   left=0.09, right=0.99, top=0.78, bottom=0.36)
        left = figure.add_subplot(grid[0, 0])
        right = figure.add_subplot(grid[0, 1])
        draw_panel_pair(left, right, panel_data, key, letter=False)
        figure.suptitle(f"{key.lower()}  {PANEL_SPECS[key]['title']}", fontsize=8.5,
                        fontweight="bold", x=0.02, ha="left", y=0.99)
        for axis, label in [(left, STRATEGY_ALN), (right, STRATEGY_ASM)]:
            axis.text(0.5, 1.30, label, transform=axis.transAxes, ha="center",
                      fontsize=7.5, fontweight="bold")
        figure.legend(handles=legend_handles(), loc="lower center",
                      bbox_to_anchor=(0.5, 0.0), ncols=4, frameon=False, fontsize=6)
        paths += save(figure, PANEL_DIR / f"{STEM}_panel_{key}")
        plt.close(figure)
    return paths


def plot_ram(aln: pd.DataFrame, perf: pd.DataFrame | None) -> list[Path]:
    """Supplementary RAM panel, only if peak RAM is complete for every tool."""
    aln_missing = aln[aln["peak_ram_gb"].isna()]
    if len(aln_missing):
        tools = sorted(aln_missing["configuration"].unique())
        warn("Supplementary RAM panel NOT generated: aligner peak_ram_gb missing for "
             + ", ".join(tools) + ". No RAM value is substituted (e.g. the "
             "command_ram_limit_gb allocation is a limit, not a measurement).")
        return []
    if perf is None or perf["peak_rss_gb"].isna().any():
        warn("Supplementary RAM panel NOT generated: assembler peak RAM unavailable.")
        return []
    data = pd.concat([
        aln.assign(value=aln["peak_ram_gb"])[["strategy", "configuration", "sample",
                                              "technology", "value"]],
        perf.assign(value=perf["peak_rss_gb"])[["strategy", "configuration", "sample",
                                                "technology", "value"]],
    ])
    figure = plt.figure(figsize=(FIGURE_WIDTH_IN, 3.0))
    grid = figure.add_gridspec(1, 2, width_ratios=[8, 5], wspace=0.28,
                               left=0.09, right=0.99, top=0.85, bottom=0.36)
    left, right = figure.add_subplot(grid[0, 0]), figure.add_subplot(grid[0, 1])
    draw_strip(left, data[data["strategy"] == STRATEGY_ALN], aln_configs())
    draw_strip(right, data[data["strategy"] == STRATEGY_ASM], asm_configs(),
               shade_configs=(VERKKO_CONFIG,))
    top = np.nanmax(data["value"]) * 1.12
    for axis in (left, right):
        axis.set_ylim(0, top)
    left.set_ylabel("Peak RAM (GB)")
    figure.legend(handles=legend_handles(), loc="lower center", bbox_to_anchor=(0.5, 0.0),
                  ncols=4, frameon=False, fontsize=6)
    paths = save(figure, PANEL_DIR / f"{STEM}_supplementary_peak_ram")
    plt.close(figure)
    return paths


# ============================================================
# STATIC WARNINGS (definitions verified against the source tables)
# ============================================================


def record_definition_warnings(aln: pd.DataFrame) -> None:
    warn(
        "Panel A: both sides are substitution rates. Aligner value = samtools stats "
        "'mismatches' (summed NM tag = substitutions + inserted + deleted bases) minus the "
        "inserted and deleted bases of the samtools stats ID records, per CIGAR-mapped read "
        "base. ID records cover indels of 1-300 bp only, so bases of longer indels remain "
        "in the aligner value, which is therefore an upper bound. QUAST '# mismatches per "
        "100 kbp' counts substitutions only per aligned contig base. Raw reads (with "
        "sequencing errors) vs consensus contigs: separate subpanels, independent y-axes; "
        "not merged, not tested."
    )
    warn(
        "Panel B: both sides are indel EVENT counts (not indel bases) per 100 kbp of aligned "
        "query bases; reads count CIGAR I/D operations of 1-300 bp only (samtools stats ID "
        "records; longer indels are skipped by samtools). Independent y-axes. They remain "
        "separate subpanels because the "
        "query differs (individual raw reads, each locus counted ~30x, vs consensus contigs) "
        "and QUAST only counts indels inside alignments >=95% identity (indels >200 bp become "
        "local misassemblies). No statistical comparison is made across strategies."
    )
    warn(
        "Panel C: CIGAR-aligned yield (denominator: input sequencing bases) and genome "
        "fraction (denominator: GRCh38) are different metrics with different denominators "
        "(row title: 'Alignment yield and assembly completeness'). They are not combined into "
        "one score and not compared statistically."
    )
    cpu_missing = sorted(aln.loc[aln["cpu_time_seconds"].isna(), "configuration"].unique())
    if cpu_missing:
        warn(
            "Panel D: measured CPU time is missing for " + ", ".join(cpu_missing)
            + "; all tools are therefore reported as allocated thread-hours "
            "(wall-clock hours x allocated threads), not CPU-hours."
        )
    warn(
        "Verkko (--hifi/--nano only, no trio or Hi-C) was evaluated as its single combined "
        "assembly.fasta (6.08-6.14 Gb, duplication ratio 2.04-2.09); no per-haplotype FASTA "
        "exists and phasing status is not inferred. Flye and GoldRush are collapsed-haploid. Duplication ratio and absolute misassembly counts are "
        "deliberately not plotted; Verkko's genome fraction and per-100 kbp rates are shown on "
        "a hatched background and should not be read as equivalent to collapsed assemblies. "
        "Verkko input depth is recorded as 'unknown' in the source table."
    )
    warn(
        "Optional SV figure NOT generated: only alignment-based Truvari results exist "
        "(truvari/*/summary.json, HG002 only); no assembly-based SV calls benchmarked against "
        "the GIAB SV truth set were found."
    )
    warn(
        "Optional cross-strategy heatmap / PCA NOT generated: none of the candidate shared "
        "metrics is harmonized across all 39 configuration-sample rows (mismatch/indel "
        "definitions differ, assembler cost and SV metrics are absent, peak RAM incomplete). "
        "No z-scoring or imputation was performed."
    )
    warn(
        "Not compared by design: mapped reads % or mapped bases % vs genome fraction %; MQ0 vs "
        "misassemblies; secondary alignments vs duplication ratio; read N50 vs NGA50; CIGAR "
        "yield vs NGA50; unmapped reads vs misassemblies."
    )
    warn(
        "Statistics: HG002/HG003/HG004 are three biological samples (n = 3 per configuration), "
        "not technical replicates; mean +/- SD is descriptive only and no tests are run."
    )


def main() -> int:
    apply_style()
    aln = load_alignment()
    asm = load_assembly()
    perf = load_assembly_performance()
    record_definition_warnings(aln)

    panel_data = build_panel_data(aln, asm, perf)
    panels = ["A", "B", "C"] + (["D"] if perf is not None else [])

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    panel_data.to_csv(TIDY_TSV, sep="\t", index=False, float_format="%.6g")

    outputs = plot_main(panel_data, panels)
    outputs += plot_single_panels(panel_data, panels)
    outputs += plot_ram(aln, perf)

    WARNINGS_TXT.write_text(
        "".join(f"{i}. {message}\n" for i, message in enumerate(WARNINGS, 1)),
        encoding="utf-8",
    )
    print(f"Panels drawn: {', '.join(panels)}")
    print(f"Panel data: {TIDY_TSV.relative_to(PROJECT)}")
    print(f"Warnings:   {WARNINGS_TXT.relative_to(PROJECT)}")
    for path in outputs:
        print(f"Figure:     {path.relative_to(PROJECT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
