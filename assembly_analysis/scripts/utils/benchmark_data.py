"""Shared data loading for the assembler-benchmark figures.

assembler_benchmark_bars.py (final, 6-panel), assembler_benchmark_points.py
(the earlier 5-panel NGA50-based comparison round that led to the bars-only
decision), and any future variant must plot exactly the same underlying
values -- this module is the single place that data is loaded, validated,
and the plotting subset is derived, so scripts cannot drift apart from each
other or from the source table.

Runtime is intentionally absent: none of the five production assembler
workflows under assemblers/whole_genome_asm/ (ont.assembly.flye2.smk,
pb.assembly.flye2.smk, ont.assembly.goldrush.smk, pb.assembly.goldrush.smk,
hybrid.assembly.verkko.smk) capture wall-clock time, CPU, or RAM for the
actual 30x whole-genome run -- there is no `benchmark:` directive in any of
them. The only timing-adjacent files that exist are from an unrelated chr21
smoke test, not this benchmark, so no runtime panel is built rather than
substituting unmeasured data.

N50 and contig count (FINAL_METRICS panels a/b) come from the same
assembly_benchmark_30x.tsv as the reference-based metrics, but that is not
a shortcut: QUAST computes `# contigs` and `N50` directly from the assembly
FASTA, before it ever aligns anything to the reference -- see QUAST's own
report.tsv field order (these "general" stats are printed above "Genome
fraction (%)" and the other reference-based fields) and
docs/f2/archive/README_assemblers_legacy.md, which documents the same
grouping. There is no separate non-QUAST FASTA-stats tool anywhere in this
repo (`calculate_quality_metrics.py` in assembly_analysis/README.md's
"Intended structure" is a documented-but-not-yet-built placeholder, not an
existing source) -- so QUAST's `contigs`/`n50` columns *are* the
reference-independent, FASTA-derived values, not a stand-in for them.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from utils.plot_style import CONFIGURATION_ORDER, SAMPLE_ORDER


def find_repo_root(start: Path) -> Path:
    start = start.resolve()
    for candidate in [start, *start.parents]:
        if (candidate / "CONSTITUTION.md").is_file():
            return candidate
    raise RuntimeError("Could not locate lrs_benchmarking repository root")


PROJECT_ROOT = find_repo_root(Path(__file__))

SOURCE_TABLE = (
    PROJECT_ROOT / "assembly_analysis" / "tables" / "30x" / "final" / "assembly_benchmark_30x.tsv"
)
PLOT_DATA_TABLE = (
    PROJECT_ROOT / "assembly_analysis" / "tables" / "30x" / "derived" / "plot_data"
    / "assembler_benchmark_plot_data.tsv"
)

# (plotting column, panel label, panel letter) for the earlier 5-panel
# NGA50-based comparison round (assembler_benchmark_points.py still uses
# this). nga50_mb is a unit conversion (bp / 1e6) of the source table's
# nga50 column for axis readability, not a rescale of the underlying value.
METRICS = [
    ("nga50_mb", "NGA50 (Mb)", "a"),
    ("genome_fraction_pct", "Genome fraction (%)", "b"),
    ("misassemblies", "Misassemblies", "c"),
    ("mismatches_per_100kbp", "Mismatches / 100 kbp", "d"),
    ("indels_per_100kbp", "Indels / 100 kbp", "e"),
]

# The final, adopted 6-panel metric set (assembler_benchmark_bars.py):
# N50 and contig count replace NGA50 (not measurable assembler-wide in a
# way that reads honestly at a shared linear scale -- see the earlier
# comparison round) and runtime stays absent (not measured at all).
FINAL_METRICS = [
    ("n50_mb", "N50 (Mb)", "a"),
    ("contigs", "Contig count", "b"),
    ("genome_fraction_pct", "Genome fraction (%)", "c"),
    ("misassemblies", "Misassemblies", "d"),
    ("mismatches_per_100kbp", "Mismatches / 100 kbp", "e"),
    ("indels_per_100kbp", "Indels / 100 kbp", "f"),
]

RAW_NUMERIC_COLUMNS = [
    "nga50", "n50", "contigs", "genome_fraction_pct", "misassemblies",
    "mismatches_per_100kbp", "indels_per_100kbp",
]

PLOT_DATA_COLUMNS = [
    "sample", "assembler", "technology", "assembly_strategy",
    "n50", "n50_mb", "contigs", "nga50", "nga50_mb", "genome_fraction_pct",
    "misassemblies", "mismatches_per_100kbp", "indels_per_100kbp", "source_file",
]


def _validate_complete(plot_data: pd.DataFrame) -> None:
    expected = pd.MultiIndex.from_tuples(
        [
            (assembler, technology, sample)
            for assembler, technology in CONFIGURATION_ORDER
            for sample in SAMPLE_ORDER
        ],
        names=["assembler", "technology", "sample"],
    )
    observed = pd.MultiIndex.from_frame(plot_data[["assembler", "technology", "sample"]])
    missing = expected.difference(observed)
    if len(missing):
        raise ValueError(f"Missing assembler/technology/sample combinations: {list(missing)}")

    duplicated = plot_data.duplicated(subset=["assembler", "technology", "sample"], keep=False)
    if duplicated.any():
        raise ValueError(
            "Duplicated assembler/technology/sample observations found:\n"
            + plot_data.loc[duplicated, ["assembler", "technology", "sample"]].to_string(index=False)
        )


def build_plot_data() -> pd.DataFrame:
    if not SOURCE_TABLE.exists():
        raise FileNotFoundError(
            f"The source table was not found:\n{SOURCE_TABLE}\n"
            "Run assembly_analysis/scripts/metrics/build_assembly_benchmark_table.py first."
        )

    data = pd.read_csv(SOURCE_TABLE, sep="\t")
    for column in RAW_NUMERIC_COLUMNS:
        data[column] = pd.to_numeric(data[column], errors="raise")

    data["nga50_mb"] = data["nga50"] / 1_000_000
    data["n50_mb"] = data["n50"] / 1_000_000
    data["assembly_strategy"] = data["technology"].map(
        lambda technology: "hybrid" if technology == "hybrid" else "single-technology"
    )

    plot_data = data[PLOT_DATA_COLUMNS].copy()
    _validate_complete(plot_data)
    plot_data = plot_data.sort_values(["assembler", "technology", "sample"]).reset_index(drop=True)
    return plot_data


def verify_plot_data_matches_source(plot_data: pd.DataFrame) -> None:
    """Re-read the source table independently and assert every plotted
    number is identical to it -- not just "derived from it once".
    """
    fresh_source = pd.read_csv(SOURCE_TABLE, sep="\t")
    for column in RAW_NUMERIC_COLUMNS:
        fresh_source[column] = pd.to_numeric(fresh_source[column], errors="raise")

    merged = plot_data.merge(
        fresh_source[["assembler", "sample", "technology", *RAW_NUMERIC_COLUMNS]],
        on=["assembler", "sample", "technology"],
        suffixes=("", "_source_check"),
    )
    for column in RAW_NUMERIC_COLUMNS:
        if not (merged[column] == merged[f"{column}_source_check"]).all():
            raise ValueError(
                f"assembler_benchmark_plot_data.tsv value for '{column}' does not match "
                f"{SOURCE_TABLE.name} -- data integrity check failed."
            )


def compute_y_limits(plot_data: pd.DataFrame, metrics=METRICS) -> dict[str, float]:
    """One y_max per metric (per-metric true-zero baseline)."""
    return {column: float(plot_data[column].max()) * 1.12 for column, _, _ in metrics}


def write_plot_data(plot_data: pd.DataFrame, path: Path = PLOT_DATA_TABLE) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    plot_data[PLOT_DATA_COLUMNS].to_csv(path, sep="\t", index=False)


# The QUAST-only 6-panel metric set (assembler_quast_final_points.py):
# NGA50 back in (as the points version's log-friendly-enough scale showed it
# can read honestly with individual points + a shared y-axis, unlike the
# all-bars figure) plus duplication ratio, which explains why Flye-HiFi's
# N50/contig-count look so different from Flye-ONT's (see
# assembler_benchmark_bars.py's docstring) and was already sitting in
# assembly_benchmark_30x.tsv but never plotted. N50 and contig count are
# left out here on purpose -- this figure is QUAST reference-based +
# duplication ratio only, not a repeat of the FASTA-stats panels.
QUAST_FINAL_METRICS = [
    ("nga50_mb", "NGA50 (Mb)", "a"),
    ("genome_fraction_pct", "Genome fraction (%)", "b"),
    ("duplication_ratio", "Duplication ratio", "c"),
    ("misassemblies", "Misassemblies", "d"),
    ("mismatches_per_100kbp", "Mismatches / 100 kbp", "e"),
    ("indels_per_100kbp", "Indels / 100 kbp", "f"),
]

# assembly_benchmark_30x.tsv column -> the QUAST report.tsv field it was taken
# from, for the QUAST-only figures (assembler_quast_final_points.py).
QUAST_REPORT_FIELDS = {
    "nga50": "NGA50",
    "misassemblies": "# misassemblies",
    "genome_fraction_pct": "Genome fraction (%)",
    "duplication_ratio": "Duplication ratio",
    "mismatches_per_100kbp": "# mismatches per 100 kbp",
    "indels_per_100kbp": "# indels per 100 kbp",
}


def load_quast_values(columns: list[str]) -> pd.DataFrame:
    """One row per assembler/technology/sample with the requested columns
    from assembly_benchmark_30x.tsv, each cross-checked against the raw QUAST
    report.tsv named in that row's source_file (not just against the table).
    """
    data = pd.read_csv(SOURCE_TABLE, sep="\t")
    for column in columns:
        data[column] = pd.to_numeric(data[column], errors="raise")

    for _, row in data.iterrows():
        report_path = PROJECT_ROOT / row["source_file"]
        report = dict(
            line.rstrip("\n").split("\t", 1)
            for line in report_path.read_text(encoding="utf-8").splitlines()
            if "\t" in line
        )
        for column in columns:
            raw = float(report[QUAST_REPORT_FIELDS[column]])
            if raw != row[column]:
                raise ValueError(
                    f"{row['assembler']} {row['technology']} {row['sample']} {column}: "
                    f"table={row[column]} but {row['source_file']}={raw}"
                )

    plot_data = data[["sample", "assembler", "technology", *columns, "source_file"]].copy()
    _validate_complete(plot_data)
    return plot_data.sort_values(["assembler", "technology", "sample"]).reset_index(drop=True)


def print_source_audit(plot_data: pd.DataFrame, metrics) -> None:
    """Print the exact value used for every sample x assembler x technology
    x metric, with its source file, before a figure is saved -- so the
    printed values can be checked against assembly_benchmark_30x.tsv /
    the raw QUAST report.tsv by eye, not just by the automated equality
    check in verify_plot_data_matches_source().
    """
    print(f"{'assembler':<9} {'technology':<11} {'sample':<7} {'metric':<20} {'value':>14}  source_file")
    for _, row in plot_data.iterrows():
        for column, label, _ in metrics:
            print(
                f"{row['assembler']:<9} {row['technology']:<11} {row['sample']:<7} "
                f"{label:<20} {row[column]:>14}  {row['source_file']}"
            )
