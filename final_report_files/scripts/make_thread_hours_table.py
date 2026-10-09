#!/usr/bin/env python3
"""Build the thread-hours table (allocated threads x wall-clock hours) shown after Figure 8.

Input:  alignment_analysis/tables/30x/derived/plot_data/wallclock_runtime_30x_plotting_values.tsv
        (written by alignment_analysis/scripts/30x/plots/final/plot_runtime_30x.py)
Output: f2-thesis-report/tables/thread_hours_30x.tsv, mean $\\pm$ SD across HG002/HG003/HG004,
        formatted like tables/alignment_summary_table_30x.tsv.
"""

from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
INPUT_TSV = (
    REPO / "alignment_analysis" / "tables" / "30x" / "derived" / "plot_data"
    / "wallclock_runtime_30x_plotting_values.tsv"
)
OUTPUT_TSV = REPO.parent / "f2-thesis-report" / "tables" / "thread_hours_30x.tsv"

TECHNOLOGY_ORDER = ["ONT", "PacBio"]
ALIGNER_ORDER = ["minimap2", "pbmm2", "VACmap", "VG Giraffe"]
TECHNOLOGY_LABELS = {"ONT": "ONT", "PacBio": "PacBio HiFi"}


def mean_sd(values: pd.Series, precision: int) -> str:
    return f"{values.mean():.{precision}f} $\\pm$ {values.std(ddof=1):.{precision}f}"


def main() -> int:
    data = pd.read_csv(INPUT_TSV, sep="\t")
    data["thread_hours"] = data["wallclock_runtime_hours"] * data["threads"]

    rows = []
    for technology in TECHNOLOGY_ORDER:
        for index, aligner in enumerate(ALIGNER_ORDER):
            runs = data[(data["read_technology"] == technology) & (data["aligner"] == aligner)]
            if len(runs) != 3:
                raise ValueError(f"Expected 3 samples for {technology}/{aligner}, found {len(runs)}")
            threads = runs["threads"].unique()
            if len(threads) != 1:
                raise ValueError(f"Mixed thread counts for {technology}/{aligner}: {threads}")
            rows.append({
                "technology": TECHNOLOGY_LABELS[technology] if index == 0 else "",
                "aligner": aligner,
                "threads": int(threads[0]),
                "runtime_hours": mean_sd(runs["wallclock_runtime_hours"], 2),
                "thread_hours": mean_sd(runs["thread_hours"], 1),
            })

    OUTPUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(OUTPUT_TSV, sep="\t", index=False)
    print(f"Wrote {OUTPUT_TSV}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
