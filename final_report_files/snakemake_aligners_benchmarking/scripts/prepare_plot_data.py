#!/usr/bin/env python3

import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--quality-report", required=True)
    parser.add_argument("--cleaned-data", required=True)
    parser.add_argument("--pivoted-data", required=True)
    parser.add_argument("--arrays", required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t")

    numeric_columns = [
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
        "non_primary_alignments",
    ]

    for column in numeric_columns:
        if column in df.columns:
            df[column] = pd.to_numeric(df[column], errors="coerce")

    Path(args.cleaned_data).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.cleaned_data, sep="\t", index=False)

    pivot = df.copy()
    pivot.to_csv(args.pivoted_data, sep="\t", index=False)

    np.savez(
        args.arrays,
        error_percent=df["error_percent"].to_numpy() if "error_percent" in df.columns else np.array([]),
        bases_mapped=df["bases_mapped"].to_numpy() if "bases_mapped" in df.columns else np.array([]),
        bases_mapped_cigar=df["bases_mapped_cigar"].to_numpy() if "bases_mapped_cigar" in df.columns else np.array([]),
    )


if __name__ == "__main__":
    main()
