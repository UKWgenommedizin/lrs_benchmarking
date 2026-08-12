#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd


def add(rows, check, status, detail):
    rows.append({"check": check, "status": status, "detail": detail})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t")
    rows = []

    add(rows, "rows_exist", "PASS" if len(df) > 0 else "FAIL", f"rows={len(df)}")

    required = ["sample", "read_technology", "aligner", "statistics_file"]
    missing_required = [column for column in required if column not in df.columns]
    add(
        rows,
        "required_columns",
        "PASS" if not missing_required else "FAIL",
        "missing=" + ",".join(missing_required) if missing_required else "all required columns present",
    )

    if not missing_required:
        duplicates = df.duplicated(subset=["sample", "read_technology", "aligner"]).sum()
        add(rows, "duplicate_conditions", "PASS" if duplicates == 0 else "FAIL", f"duplicates={duplicates}")

    missing_values = int(df.isna().sum().sum())
    add(rows, "missing_values", "PASS" if missing_values == 0 else "WARN", f"missing_values={missing_values}")

    numeric_columns = [
        "raw_total_sequences",
        "reads_mapped",
        "reads_unmapped",
        "total_length",
        "bases_mapped",
        "bases_mapped_cigar",
        "error_rate",
        "error_percent",
    ]

    existing_numeric = [column for column in numeric_columns if column in df.columns]
    negative_count = 0

    for column in existing_numeric:
        values = pd.to_numeric(df[column], errors="coerce")
        negative_count += int((values < 0).sum())

    add(rows, "negative_numeric_values", "PASS" if negative_count == 0 else "FAIL", f"negative_values={negative_count}")

    report = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
