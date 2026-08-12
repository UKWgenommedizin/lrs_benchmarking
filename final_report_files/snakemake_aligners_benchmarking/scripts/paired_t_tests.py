#!/usr/bin/env python3

import argparse
from pathlib import Path
import itertools
import pandas as pd
from scipy.stats import ttest_rel


def safe_ttest(a, b):
    if len(a) < 2 or len(b) < 2:
        return float("nan"), float("nan")
    result = ttest_rel(a, b, nan_policy="omit")
    return result.statistic, result.pvalue


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t")

    metrics = [
        "error_percent",
        "mapped_reads_percent",
        "mapped_bases_percent",
        "bases_mapped",
        "bases_mapped_cigar",
    ]

    rows = []

    for technology in sorted(df["read_technology"].dropna().unique()):
        tech_df = df[df["read_technology"] == technology]

        aligners = sorted(tech_df["aligner"].dropna().unique())

        for metric in metrics:
            if metric not in tech_df.columns:
                continue

            metric_df = tech_df[["sample", "aligner", metric]].copy()
            metric_df[metric] = pd.to_numeric(metric_df[metric], errors="coerce")

            wide = metric_df.pivot_table(index="sample", columns="aligner", values=metric, aggfunc="mean")

            for aligner_a, aligner_b in itertools.combinations(aligners, 2):
                if aligner_a not in wide.columns or aligner_b not in wide.columns:
                    continue

                paired = wide[[aligner_a, aligner_b]].dropna()
                stat, pvalue = safe_ttest(paired[aligner_a], paired[aligner_b])

                rows.append(
                    {
                        "read_technology": technology,
                        "metric": metric,
                        "aligner_a": aligner_a,
                        "aligner_b": aligner_b,
                        "n_pairs": len(paired),
                        "t_statistic": stat,
                        "p_value": pvalue,
                    }
                )

    out = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
