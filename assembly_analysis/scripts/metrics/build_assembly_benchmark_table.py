#!/usr/bin/env python3
"""Build the canonical, curated 30x assembly benchmark table from QUAST metrics.

Input -> assembly_analysis/tables/assembler_metrics.tsv
         (the long/tidy table written by extract_quast_metrics.py: one row
         per assembler/dataset/metric)

Output -> assembly_analysis/tables/30x/final/assembly_benchmark_30x.tsv
          (one row per assembler/sample/technology, wide format), mirroring
          the aligner benchmark table at
          alignment_analysis/tables/30x/final/alignment_benchmark_30x.tsv

The column set follows the QUAST categories from the benchmark workflow
diagram (reference agreement: genome fraction, misassemblies, mismatches and
indels per 100 kbp; contiguity: NG50/NGA50, contig count, total length) plus
the extra context columns assembly_analysis/README.md's "Scientific
comparison rules" call for (technology, depth, representation, provenance).

Usage:
    python assembly_analysis/scripts/metrics/build_assembly_benchmark_table.py
    python assembly_analysis/scripts/metrics/build_assembly_benchmark_table.py \\
        --input /path/to/assembler_metrics.tsv \\
        --output /path/to/assembly_benchmark_30x.tsv
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def find_repo_root(start: Path) -> Path:
    start = start.resolve()
    for candidate in [start, *start.parents]:
        if (candidate / "CONSTITUTION.md").is_file():
            return candidate
    raise RuntimeError("Could not locate lrs_benchmarking repository root")


PROJECT_ROOT = find_repo_root(Path(__file__))

DEFAULT_INPUT = PROJECT_ROOT / "assembly_analysis" / "tables" / "assembler_metrics.tsv"
DEFAULT_OUTPUT = (
    PROJECT_ROOT / "assembly_analysis" / "tables" / "30x" / "final" / "assembly_benchmark_30x.tsv"
)

# QUAST report.tsv metric name -> curated output column name. Extend this
# when a new metric earns a place in the canonical table; every other metric
# stays available, per-row, in assembler_metrics.tsv.
METRIC_COLUMNS = {
    "# contigs": "contigs",
    "Total length": "total_length_bp",
    "Largest contig": "largest_contig_bp",
    "Genome fraction (%)": "genome_fraction_pct",
    "Duplication ratio": "duplication_ratio",
    "# misassemblies": "misassemblies",
    "# mismatches per 100 kbp": "mismatches_per_100kbp",
    "# indels per 100 kbp": "indels_per_100kbp",
    "N50": "n50",
    "NG50": "ng50",
    "NA50": "na50",
    "NGA50": "nga50",
    "L50": "l50",
    "LG50": "lg50",
    "# unaligned contigs": "unaligned_contigs",
    "Unaligned length": "unaligned_length_bp",
}

CONTEXT_COLUMNS = ["assembler", "sample", "technology", "depth", "representation"]
FIELDNAMES = CONTEXT_COLUMNS + list(METRIC_COLUMNS.values()) + ["source_file"]


def read_long_table(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def pivot(long_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    groups: dict[tuple[str, str], dict[str, str]] = {}
    order: list[tuple[str, str]] = []

    for row in long_rows:
        key = (row["assembler"], row["dataset"])
        if key not in groups:
            groups[key] = {column: "" for column in FIELDNAMES}
            for column in CONTEXT_COLUMNS:
                groups[key][column] = row[column]
            groups[key]["source_file"] = row["source_file"]
            order.append(key)

        column = METRIC_COLUMNS.get(row["metric"])
        if column is not None:
            groups[key][column] = row["value"]

    return [groups[key] for key in order]


def write_wide_table(rows: list[dict[str, str]], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    rows = sorted(rows, key=lambda row: (row["assembler"], row["sample"], row["technology"]))
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help=f"Long/tidy QUAST metrics table (default: {DEFAULT_INPUT})",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Curated wide output TSV path (default: {DEFAULT_OUTPUT})",
    )
    args = parser.parse_args()

    if not args.input.is_file():
        raise SystemExit(
            f"{args.input} not found. Run "
            "assembly_analysis/scripts/metrics/extract_quast_metrics.py first."
        )

    long_rows = read_long_table(args.input)
    wide_rows = pivot(long_rows)
    write_wide_table(wide_rows, args.output)

    print(f"Built {len(wide_rows)} assembler/dataset row(s) from {len(long_rows)} metric rows.")
    try:
        print(f"Wrote: {args.output.relative_to(PROJECT_ROOT)}")
    except ValueError:
        print(f"Wrote: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
