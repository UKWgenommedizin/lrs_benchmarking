#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


ALIGNER_ORDER = ["minimap2", "pbmm2", "vacmap", "vg_giraffe"]
TECH_ORDER = ["ont", "pb"]


def clean_label(value):
    mapping = {
        "ont": "ONT",
        "pb": "PacBio HiFi",
        "minimap2": "minimap2",
        "pbmm2": "pbmm2",
        "vacmap": "VACmap",
        "vg_giraffe": "vg Giraffe",
    }
    return mapping.get(value, value)


def plot_metric(df, metric, ylabel, title, png, pdf):
    df = df.copy()

    if metric not in df.columns:
        raise ValueError(f"Missing metric column: {metric}")

    df[metric] = pd.to_numeric(df[metric], errors="coerce")
    df = df.dropna(subset=[metric])

    if df.empty:
        raise ValueError(f"No valid numeric data for metric: {metric}")

    aligners = [a for a in ALIGNER_ORDER if a in set(df["aligner"])]
    technologies = [t for t in TECH_ORDER if t in set(df["read_technology"])]

    fig, ax = plt.subplots(figsize=(10, 6))

    x_positions = list(range(len(aligners)))
    offsets = {}

    if len(technologies) == 1:
        offsets[technologies[0]] = 0
    else:
        offsets = {technologies[0]: -0.15, technologies[1]: 0.15}

    for tech in technologies:
        tech_df = df[df["read_technology"] == tech]

        means = []
        xs = []

        for i, aligner in enumerate(aligners):
            values = tech_df[tech_df["aligner"] == aligner][metric]
            means.append(values.mean())
            xs.append(i + offsets.get(tech, 0))

            sample_values = tech_df[tech_df["aligner"] == aligner]
            ax.scatter(
                [i + offsets.get(tech, 0)] * len(sample_values),
                sample_values[metric],
                alpha=0.75,
                s=45,
            )

        ax.plot(xs, means, marker="o", linewidth=2, label=clean_label(tech))

    ax.set_xticks(x_positions)
    ax.set_xticklabels([clean_label(a) for a in aligners], rotation=25, ha="right")
    ax.set_ylabel(ylabel)
    ax.set_xlabel("Alignment configuration")
    ax.set_title(title)
    ax.legend(title="Read technology")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()

    Path(png).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cleaned-data", required=True)
    parser.add_argument("--pivoted-data", required=True)
    parser.add_argument("--arrays", required=True)
    parser.add_argument("--tests", required=True)
    parser.add_argument("--error-png", required=True)
    parser.add_argument("--error-pdf", required=True)
    parser.add_argument("--mapped-png", required=True)
    parser.add_argument("--mapped-pdf", required=True)
    parser.add_argument("--cigar-png", required=True)
    parser.add_argument("--cigar-pdf", required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.cleaned_data, sep="\t")

    plot_metric(
        df=df,
        metric="error_percent",
        ylabel="Alignment error rate (%)",
        title="Alignment error rate across aligners",
        png=args.error_png,
        pdf=args.error_pdf,
    )

    plot_metric(
        df=df,
        metric="bases_mapped",
        ylabel="Mapped bases",
        title="Mapped bases across aligners",
        png=args.mapped_png,
        pdf=args.mapped_pdf,
    )

    plot_metric(
        df=df,
        metric="bases_mapped_cigar",
        ylabel="Mapped bases from CIGAR",
        title="CIGAR mapped bases across aligners",
        png=args.cigar_png,
        pdf=args.cigar_pdf,
    )


if __name__ == "__main__":
    main()
