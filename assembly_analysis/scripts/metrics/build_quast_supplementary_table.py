#!/usr/bin/env python3
"""Supplementary table: the full QUAST statistics for every 30x assembly.

Reads the long-format assembler_metrics.tsv (written by extract_quast_metrics.py
from assembly_quality/quast/{assembler}/{dataset}/report.tsv) and pivots the
QUAST fields below into one row per sample x assembler x input. Values are
copied verbatim; QUAST's "-" (metric undefined, e.g. NGA90 when the assembly
aligns to < 90% of the reference) becomes NA. Nothing is recomputed.

    assembly_analysis/tables/30x/final/assembly_quast_supplementary_30x.tsv

Usage:
    python assembly_analysis/scripts/metrics/extract_quast_metrics.py
    python assembly_analysis/scripts/metrics/build_quast_supplementary_table.py
"""
from __future__ import annotations

import csv
from pathlib import Path


def find_repo_root(start: Path) -> Path:
    start = start.resolve()
    for candidate in [start, *start.parents]:
        if (candidate / "CONSTITUTION.md").is_file():
            return candidate
    raise RuntimeError("Could not locate lrs_benchmarking repository root")


PROJECT_ROOT = find_repo_root(Path(__file__))
SOURCE = PROJECT_ROOT / "assembly_analysis" / "tables" / "assembler_metrics.tsv"
OUTPUT = PROJECT_ROOT / "assembly_analysis" / "tables" / "30x" / "final" / "assembly_quast_supplementary_30x.tsv"

# (output column, QUAST report.tsv field)
QUAST_FIELDS = [
    ("total_length_bp", "Total length"),
    ("sequence_count", "# contigs"),
    ("largest_contig_bp", "Largest contig"),
    ("largest_alignment_bp", "Largest alignment"),
    ("n50_bp", "N50"),
    ("ng50_bp", "NG50"),
    ("nga50_bp", "NGA50"),
    ("l50", "L50"),
    ("lga50", "LGA50"),
    ("ns_per_100kbp", "# N's per 100 kbp"),
    ("genome_fraction_pct", "Genome fraction (%)"),
    ("duplication_ratio", "Duplication ratio"),
    ("misassemblies", "# misassemblies"),
    ("local_misassemblies", "# local misassemblies"),
    ("scaffold_gap_ext_mis", "# scaffold gap ext. mis."),
    ("scaffold_gap_loc_mis", "# scaffold gap loc. mis."),
    ("mismatches_per_100kbp", "# mismatches per 100 kbp"),
    ("indels_per_100kbp", "# indels per 100 kbp"),
    ("unaligned_length_bp", "Unaligned length"),
]

ASSEMBLER_NAMES = {"flye": "Flye", "goldrush": "GoldRush", "verkko": "Verkko"}
INPUT_NAMES = {"ont": "ONT", "pb": "HiFi", "hybrid": "ONT+HiFi"}
REPRESENTATION = {
    "flye": "single-technology assembly",
    "goldrush": "single-technology assembly, ntLink-scaffolded (N50/NG50/count are scaffold-level)",
    "verkko": "hybrid assembly, both haplotypes in one FASTA (diploid)",
}
ROW_ORDER = [("flye", "ont"), ("flye", "pb"), ("goldrush", "ont"), ("goldrush", "pb"), ("verkko", "hybrid")]


def main() -> int:
    records: dict[tuple[str, str, str], dict[str, str]] = {}
    with SOURCE.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row["assembler"], row["technology"], row["sample"])
            entry = records.setdefault(key, {"source_file": row["source_file"]})
            entry[row["metric"]] = row["value"]

    fields = ["sample", "assembler", "input", "representation", *(name for name, _ in QUAST_FIELDS), "source_file"]
    rows = []
    for assembler, technology in ROW_ORDER:
        for sample in ("HG002", "HG003", "HG004"):
            metrics = records.get((assembler, technology, sample))
            if metrics is None:
                raise SystemExit(f"Missing QUAST report for {assembler} {technology} {sample}")
            row = {
                "sample": sample,
                "assembler": ASSEMBLER_NAMES[assembler],
                "input": INPUT_NAMES[technology],
                "representation": REPRESENTATION[assembler],
                "source_file": metrics["source_file"],
            }
            for column, quast_field in QUAST_FIELDS:
                value = metrics.get(quast_field)
                if value is None:
                    raise SystemExit(f"'{quast_field}' absent from {metrics['source_file']}")
                row[column] = "NA" if value.strip() == "-" else value
            rows.append(row)

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} rows: {OUTPUT.relative_to(PROJECT_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
