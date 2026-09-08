#!/usr/bin/env python3

"""
Build one combined alignment summary table from samtools stats reports.

Currently supported alignment configurations:

1. ONT reads aligned with minimap2 map-ont
2. ONT reads aligned with pbmm2 SUBREAD
3. ONT reads aligned with VACmap
4. ONT reads aligned with VG Giraffe
5. PacBio reads aligned with minimap2 map-hifi
6. PacBio reads aligned with pbmm2 CCS
7. PacBio reads aligned with VACmap
8. PacBio reads aligned with VG Giraffe

Input files are searched in:

    alignment_analysis/tables/
    cram/

The combined table is written to:

    alignment_analysis/tables/alignment_summary.tsv
"""

import argparse
import csv
import re
import sys
from pathlib import Path


#Project paths
PROJECT = Path(__file__).resolve().parents[2]
TABLES = PROJECT / "alignment_analysis" / "tables"
CRAM = PROJECT / "cram"
OUTPUT = TABLES / "alignment_summary.tsv"


#Supported alignment method tags
METHODS = [
    "mm2-ont",
    "mm2-pb",
    "pbmm2-ccs",
    "pbmm2-subread",
    "vacmap-ont",
    "vacmap-pb",
    "vg-ont",
    "vg-pb",]


#Read the SN section of a samtools stats file
def parse_sn_file(path: Path) -> dict[str, str]:
    """Read the SN summary fields from one samtools stats report."""

    values: dict[str, str] = {}

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("SN\t"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 3:
                continue

            name = fields[1].strip().rstrip(":")
            value = fields[2].strip()

            values[name] = value

    return values


#Convert one samtools value to a number
def get_number(
    statistics: dict[str, str],
    name: str,
    *,
    integer: bool = False,
    default=None,):
    """
    Return one numeric SN value.

    Samtools values sometimes contain explanatory text after the number.
    Only the first field is converted.
    """

    raw_value = statistics.get(name)

    if raw_value is None:
        return default

    number_text = raw_value.split()[0].replace(",", "")

    try:
        number = float(number_text)
    except ValueError:
        return default

    return int(number) if integer else number


#Determine the sample, technology and alignment method from the filename
def describe_file(path: Path) -> tuple[str, str, str]:
    """
    Return:

        sample
        technology
        method tag
    """

    name = path.name.lower()

    sample_match = re.search(r"(hg00[234])", name)

    if sample_match is None:
        raise ValueError(f"Sample could not be detected from filename: {name}")

    sample = sample_match.group(1).upper()

    if ".ont." in name:
        technology = "ont"
    elif ".pb." in name:
        technology = "pb"
    else:
        raise ValueError(f"Technology could not be detected from filename: {name}")

    method_tag = None

    for method in METHODS:
        if method in name:
            method_tag = method
            break

    if method_tag is None:
        if technology == "ont":
            method_tag = "mm2-ont"
        else:
            method_tag = "mm2-pb"

    return sample, technology, method_tag


#Determine aligner metadata
def describe_alignment(
    technology: str,
    method_tag: str,) -> tuple[str, str, str]:
    """
    Return:

        aligner
        preset
        configuration
    """

    if method_tag == "mm2-ont":
        return "minimap2", "map-ont", "mm2-ont"

    if method_tag == "mm2-pb":
        return "minimap2", "map-hifi", "mm2-pb"

    if method_tag == "pbmm2-ccs":
        return "pbmm2", "CCS/HIFI", "pbmm2-pb"

    if method_tag == "pbmm2-subread":
        return "pbmm2", "SUBREAD", "pbmm2-ont"

    if method_tag == "vacmap-ont":
        return "VACmap", "vacmap-ont", "vacmap-ont"

    if method_tag == "vacmap-pb":
        return "VACmap", "vacmap-pb", "vacmap-pb"

    if method_tag == "vg-ont":
        return "VG Giraffe", "vg-ont", "vg-ont"

    if method_tag == "vg-pb":
        return "VG Giraffe", "vg-pb", "vg-pb"

    raise ValueError(f"Unsupported alignment method: {method_tag}")


#Select one statistics file for every sample, technology and alignment method
def choose_files() -> dict[tuple[str, str, str], Path]:
    """
    Find all supported samtools stats reports.

    The selection key is:

        sample, technology, alignment method

    The .1k statistics file is preferred when multiple files exist.
    """

    selected: dict[tuple[str, str, str], Path] = {}

    search_directories = [
        TABLES,
        CRAM,]

    for directory in search_directories:
        if not directory.is_dir():
            continue

        for path in sorted(directory.rglob("*")):
            if not path.is_file():
                continue

            if not path.name.endswith(
                (
                    ".samtools_stats.txt",
                    ".cram.stats",
                    ".stats.txt",
                    ".stats",)):
                continue

            try:
                sample, technology, method_tag = describe_file(path)
            except ValueError:
                continue

            dataset = (sample, technology, method_tag)

            if dataset not in selected:
                selected[dataset] = path
                continue

            current_name = selected[dataset].name
            new_name = path.name

            if ".1k." in new_name and ".1k." not in current_name:
                selected[dataset] = path

    return selected


#Build one table row
def build_row(
    path: Path,
    sample: str,
    technology: str,
    method_tag: str,
) -> dict[str, object]:
    """Create one combined summary row."""

    stats = parse_sn_file(path)

    raw_total = get_number(
        stats,
        "raw total sequences",
        integer=True,)

    reads_mapped = get_number(
        stats,
        "reads mapped",
        integer=True,)

    if raw_total is None:
        raise ValueError("missing 'raw total sequences'")

    if reads_mapped is None:
        raise ValueError("missing 'reads mapped'")

    reads_unmapped = get_number(
        stats,
        "reads unmapped",
        integer=True,
        default=raw_total - reads_mapped,)

    total_length = get_number(
        stats,
        "total length",
        integer=True,
        default=0,)

    bases_mapped = get_number(
        stats,
        "bases mapped",
        integer=True,
        default=0,)

    bases_mapped_cigar = get_number(
        stats,
        "bases mapped (cigar)",
        integer=True,
        default=bases_mapped,)

    mismatches = get_number(
        stats,
        "mismatches",
        integer=True,
        default=0,)

    error_rate = get_number(
        stats,
        "error rate",
        default=0.0,)

    average_length = get_number(
        stats,
        "average length",
        default=0.0,)

    maximum_length = get_number(
        stats,
        "maximum length",
        integer=True,
        default=0,)

    non_primary_alignments = get_number(
        stats,
        "non-primary alignments",
        integer=True,
        default=0,)

    mapped_reads_percent = (
        reads_mapped / raw_total * 100
        if raw_total > 0
        else 0.0)

    mapped_bases_percent = (
        bases_mapped / total_length * 100
        if total_length > 0
        else 0.0)

    error_percent = error_rate * 100

    aligner, preset, configuration = describe_alignment(
        technology,
        method_tag,)

    technology_label = {
        "ont": "ONT",
        "pb": "PacBio",}[technology]

    try:
        statistics_file = str(path.relative_to(PROJECT))
    except ValueError:
        statistics_file = str(path)

    return {
        "sample": sample,
        "read_technology": technology_label,
        "aligner": aligner,
        "preset": preset,
        "configuration": configuration,
        "statistics_file": statistics_file,
        "raw_total_sequences": raw_total,
        "reads_mapped": reads_mapped,
        "reads_unmapped": reads_unmapped,
        "mapped_reads_percent": f"{mapped_reads_percent:.4f}",
        "total_length": total_length,
        "bases_mapped": bases_mapped,
        "bases_mapped_cigar": bases_mapped_cigar,
        "mapped_bases_percent": f"{mapped_bases_percent:.4f}",
        "mismatches": mismatches,
        "error_rate": f"{error_rate:.8f}",
        "error_percent": f"{error_percent:.4f}",
        "insertions": None,
        "deletions": None,
        "average_length": f"{average_length:.2f}",
        "maximum_length": maximum_length,
        "non_primary_alignments": non_primary_alignments,}


#Sort rows consistently
def sorting_key(row: dict[str, object]) -> tuple[int, int, int]:
    """Define the presentation order in the output table."""

    sample_order = {
        "HG002": 0,
        "HG003": 1,
        "HG004": 2,}

    technology_order = {
        "ONT": 0,
        "PacBio": 1,}

    configuration_order = {
        "mm2-ont": 0,
        "mm2-pb": 1,
        "pbmm2-ont": 2,
        "pbmm2-pb": 3,
        "vacmap-ont": 4,
        "vacmap-pb": 5,
        "vg-ont": 6,
        "vg-pb": 7,}

    return (
        sample_order.get(str(row["sample"]), 99),
        technology_order.get(str(row["read_technology"]), 99),
        configuration_order.get(str(row["configuration"]), 99),)


#Main workflow
def main() -> int:
    """Build and write the combined TSV table."""

    global PROJECT, TABLES, CRAM, OUTPUT

    parser = argparse.ArgumentParser(
        description="Combine alignment statistics for minimap2, pbmm2, VACMap and VG.")
    parser.add_argument(
        "--project",
        type=Path,
        default=PROJECT,
        help="Repository root containing alignment_analysis/ and cram/.")
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Output TSV path. Defaults to alignment_analysis/tables/alignment_summary.tsv.")
    args = parser.parse_args()

    PROJECT = args.project.resolve()
    TABLES = PROJECT / "alignment_analysis" / "tables"
    CRAM = PROJECT / "cram"
    OUTPUT = (args.out if args.out is not None else
        TABLES / "alignment_summary.tsv").resolve()

    if not TABLES.is_dir():
        print(
            f"ERROR: tables directory not found: {TABLES}",
            file=sys.stderr,
        )
        return 1

    selected = choose_files()

    if not selected:
        print(
            "ERROR: no samtools statistics files were found.",
            file=sys.stderr,
        )
        return 1

    rows: list[dict[str, object]] = []

    for dataset, path in selected.items():
        sample, technology, method_tag = dataset

        try:
            row = build_row(
                path,
                sample,
                technology,
                method_tag,)
        except ValueError as error:
            print(
                f"WARNING: {path.name}: {error}",
                file=sys.stderr,)
            continue

        rows.append(row)

    rows.sort(key=sorting_key)

    fieldnames = [
        "sample",
        "read_technology",
        "aligner",
        "preset",
        "configuration",
        "statistics_file",
        "raw_total_sequences",
        "reads_mapped",
        "reads_unmapped",
        "mapped_reads_percent",
        "total_length",
        "bases_mapped",
        "bases_mapped_cigar",
        "mapped_bases_percent",
        "mismatches",
        "error_rate",
        "error_percent",
        "insertions",
        "deletions",
        "average_length",
        "maximum_length",
        "non_primary_alignments",]

    TABLES.mkdir(parents=True, exist_ok=True)

    with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="raise",)

        writer.writeheader()
        writer.writerows(rows)

    print(f"Alignment summary written to: {OUTPUT}")
    print(f"Data rows written: {len(rows)}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
