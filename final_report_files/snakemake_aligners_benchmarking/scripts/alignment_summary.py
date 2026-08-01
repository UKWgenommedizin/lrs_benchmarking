#!/usr/bin/env python3

import argparse
from pathlib import Path
import math
import pandas as pd


def parse_number(value):
    value = str(value).strip().replace(",", "")
    try:
        if "." in value or "e" in value.lower():
            return float(value)
        return int(value)
    except ValueError:
        return math.nan


def parse_samtools_stats(path):
    data = {}

    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("SN\t"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue

            key = parts[1].strip().rstrip(":").lower()
            value = parse_number(parts[2])
            data[key] = value

    name = Path(path).name.replace(".stats.txt", "")
    parts = name.split(".", 2)

    if len(parts) == 3:
        sample, technology, aligner = parts
    else:
        sample, technology, aligner = "unknown", "unknown", name

    raw_total_sequences = data.get("raw total sequences", math.nan)
    reads_mapped = data.get("reads mapped", math.nan)
    reads_unmapped = data.get("reads unmapped", math.nan)
    total_length = data.get("total length", math.nan)
    bases_mapped = data.get("bases mapped", math.nan)
    bases_mapped_cigar = data.get("bases mapped (cigar)", math.nan)
    mismatches = data.get("mismatches", math.nan)
    error_rate = data.get("error rate", math.nan)
    insertions = data.get("insertions", math.nan)
    deletions = data.get("deletions", math.nan)
    average_length = data.get("average length", math.nan)
    maximum_length = data.get("maximum length", math.nan)
    non_primary_alignments = data.get("non-primary alignments", math.nan)

    mapped_reads_percent = (
        reads_mapped / raw_total_sequences * 100
        if raw_total_sequences and not math.isnan(raw_total_sequences)
        else math.nan
    )

    mapped_bases_percent = (
        bases_mapped_cigar / total_length * 100
        if total_length and not math.isnan(total_length)
        else math.nan
    )

    error_percent = error_rate * 100 if not math.isnan(error_rate) else math.nan

    return {
        "sample": sample,
        "read_technology": technology,
        "aligner": aligner,
        "statistics_file": str(path),
        "raw_total_sequences": raw_total_sequences,
        "reads_mapped": reads_mapped,
        "reads_unmapped": reads_unmapped,
        "mapped_reads_percent": mapped_reads_percent,
        "total_length": total_length,
        "bases_mapped": bases_mapped,
        "bases_mapped_cigar": bases_mapped_cigar,
        "mapped_bases_percent": mapped_bases_percent,
        "mismatches": mismatches,
        "error_rate": error_rate,
        "error_percent": error_percent,
        "insertions": insertions,
        "deletions": deletions,
        "average_length": average_length,
        "maximum_length": maximum_length,
        "non_primary_alignments": non_primary_alignments,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stats", nargs="+", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    rows = [parse_samtools_stats(path) for path in args.stats]
    table = pd.DataFrame(rows)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
