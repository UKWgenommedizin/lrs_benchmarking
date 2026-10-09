#!/usr/bin/env python3
"""Extract QUAST final-report metrics into a .tsv table.

Reads every ``report.tsv`` produced by the QUAST assessment workflow
Input-> (``assemblers/whole_genome_asm/assessment/assembly_quality_quast.smk``)

Extract-> assembly_quality/quast/{assembler}/{dataset}/report.tsv

Combines all into a common .tsv table

    assembly_analysis/tables/assembler_metrics.tsv

``assembly_quality/`` is gitignored raw-tool-output staging, not tracked
history -- it is populated locally (or synced from the server layout
documented in assemblers/whole_genome_asm/assessment/README.md)

Usage:
    python assembly_analysis/scripts/metrics/extract_quast_metrics.py
    python assembly_analysis/scripts/metrics/extract_quast_metrics.py \\
        --quast-root /path/to/assembly_quality/quast \\
        --output /path/to/assembler_metrics.tsv
"""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

#Find in the repository root
def find_repo_root(start: Path) -> Path:
    start = start.resolve()
    for candidate in [start, *start.parents]:
        if (candidate / "CONSTITUTION.md").is_file():
            return candidate
    raise RuntimeError("Could not locate lrs_benchmarking repository root")


PROJECT_ROOT = find_repo_root(Path(__file__))

#Output directory
DEFAULT_QUAST_ROOT = PROJECT_ROOT / "assembly_quality" / "quast"
DEFAULT_OUTPUT = PROJECT_ROOT / "assembly_analysis" / "tables" / "assembler_metrics.tsv"

#Readme contain the interpretation of duplication ratio and genome fraction correctly assembly_analysis/README.md, in this line can be also add other assemblers when needed. 
#In case new assemblers are extracted but not named here, deafult name "unknown" will be used
ASSEMBLER_REPRESENTATION = {
    "flye": "collapsed",
    "goldrush": "collapsed",
    "verkko": "haplotype_resolved",}

# Matches dataset directory names such as "HG002.ont.30x"
DATASET_RE = re.compile(
    r"^(?P<sample>[A-Za-z0-9]+)(?:\.(?P<technology>ont|pb))?(?:\.(?P<depth>[0-9]+x))?$")

#Input data set to be extracted
FIELDNAMES = [
    "assembler",
    "dataset",
    "sample",
    "technology",
    "depth",
    "representation",
    "metric",
    "value",
    "source_file",]

#Categorize dataset directory based on the sample, technology and depth
def parse_dataset(dataset: str) -> tuple[str, str, str]:

    match = DATASET_RE.match(dataset)
    if not match:
        return dataset, "unknown", "unknown"

    sample = match.group("sample")
    technology = match.group("technology") or "hybrid"
    depth = match.group("depth") or "unknown"
    return sample, technology, depth

#Extract the values of the parameters from the report.tsv QUAST file
def read_report(path: Path) -> list[tuple[str, str]]:
    """Return (metric, value) pairs from a QUAST report.tsv, in file order."""
    rows = []
    with path.open(encoding="utf-8") as handle:
        for line in csv.reader(handle, delimiter="\t"):
            if len(line) < 2:
                continue
            metric, value = line[0], line[1]
            if metric == "Assembly":
                continue
            rows.append((metric, value))
    return rows


def discover_reports(quast_root: Path) -> list[Path]:
    return sorted(quast_root.glob("*/*/report.tsv"))

#Make a list using the datasets input labels
def build_table(quast_root: Path) -> list[dict[str, str]]:
    records = []
    for report_path in discover_reports(quast_root):
        dataset_dir = report_path.parent
        assembler = dataset_dir.parent.name
        dataset = dataset_dir.name
        sample, technology, depth = parse_dataset(dataset)
        representation = ASSEMBLER_REPRESENTATION.get(assembler, "unknown")
        source_file = report_path.relative_to(PROJECT_ROOT).as_posix()

        for metric, value in read_report(report_path):
            records.append(
                {
                    "assembler": assembler,
                    "dataset": dataset,
                    "sample": sample,
                    "technology": technology,
                    "depth": depth,
                    "representation": representation,
                    "metric": metric,
                    "value": value,
                    "source_file": source_file,
                }
            )
    return records

#Write a table using the list and dictionary names extractd before
def write_table(records: list[dict[str, str]], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, delimiter="\t")
        writer.writeheader()
        writer.writerows(records)

#Make all the quast input report.tsv an integer, in case not .tsv is found not find files should be printed
def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--quast-root",
        type=Path,
        default=DEFAULT_QUAST_ROOT,
        help="Directory containing {assembler}/{dataset}/report.tsv "
        f"(default: {DEFAULT_QUAST_ROOT})",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Output tidy TSV path (default: {DEFAULT_OUTPUT})",
    )
    args = parser.parse_args()

    reports = discover_reports(args.quast_root)
    if not reports:
        raise SystemExit(
            f"No report.tsv files found under {args.quast_root}. "
            "Has the QUAST assessment workflow "
            "(assemblers/whole_genome_asm/assessment/assembly_quality_quast.smk) "
            "been run, and its output placed/synced under assembly_quality/quast/?"
        )
    
    records = build_table(args.quast_root)
    write_table(records, args.output)

    print(f"Parsed {len(reports)} report.tsv file(s) into {len(records)} metric rows.")
    try:
        print(f"Wrote: {args.output.relative_to(PROJECT_ROOT)}")
    except ValueError:
        print(f"Wrote: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
