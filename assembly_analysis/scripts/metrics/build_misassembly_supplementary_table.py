#!/usr/bin/env python3
"""Supplementary table: extensive, local and combined QUAST misassemblies
for every 30x assembly.

The main QUAST figure (assembler_quast_misassemblies.py) plots QUAST's
"# misassemblies" field, which counts extensive misassemblies only
(relocations + translocations + inversions). Several published benchmarks
(e.g. Chen2021) instead report extensive + local misassemblies, so this
table lists both fields and their per-assembly sum to make those numbers
comparable. Values are read directly from each
assembly_quality/quast/{assembler}/{dataset}/report.tsv; the combined total
is computed per assembly (extensive + local only -- no mismatches, indels or
scaffold-gap categories). The extensive count is cross-checked against the
"misassemblies" column of assembly_benchmark_30x.tsv so the table cannot
drift from the main figure.

    assembly_analysis/tables/30x/final/assembly_misassemblies_supplementary_30x.tsv

Usage:
    python assembly_analysis/scripts/metrics/build_misassembly_supplementary_table.py
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
QUAST_ROOT = PROJECT_ROOT / "assembly_quality" / "quast"
BENCHMARK_TABLE = PROJECT_ROOT / "assembly_analysis" / "tables" / "30x" / "final" / "assembly_benchmark_30x.tsv"
OUTPUT = PROJECT_ROOT / "assembly_analysis" / "tables" / "30x" / "final" / "assembly_misassemblies_supplementary_30x.tsv"

# (assembler dir, technology, QUAST dataset dir template, display assembler, display input)
CONFIGURATIONS = [
    ("flye", "ont", "{sample}.ont.30x", "Flye", "ONT"),
    ("flye", "pb", "{sample}.pb.30x", "Flye", "HiFi"),
    ("goldrush", "ont", "{sample}.ont.30x", "GoldRush", "ONT"),
    ("goldrush", "pb", "{sample}.pb.30x", "GoldRush", "HiFi"),
    ("verkko", "hybrid", "{sample}", "Verkko", "ONT+HiFi"),
]
SAMPLES = ("HG002", "HG003", "HG004")


def read_report(path: Path) -> dict[str, str]:
    return dict(
        line.rstrip("\n").split("\t", 1)
        for line in path.read_text(encoding="utf-8").splitlines()
        if "\t" in line
    )


def load_benchmark_extensive() -> dict[tuple[str, str, str], int]:
    with BENCHMARK_TABLE.open(encoding="utf-8") as handle:
        return {
            (row["assembler"], row["technology"], row["sample"]): int(row["misassemblies"])
            for row in csv.DictReader(handle, delimiter="\t")
        }


def main() -> int:
    benchmark = load_benchmark_extensive()
    fields = ["sample", "assembler", "input", "extensive_misassemblies",
              "local_misassemblies", "combined_misassemblies", "source_file"]
    rows = []
    for assembler, technology, dataset, assembler_name, input_name in CONFIGURATIONS:
        for sample in SAMPLES:
            report_path = QUAST_ROOT / assembler / dataset.format(sample=sample) / "report.tsv"
            report = read_report(report_path)
            extensive = int(report["# misassemblies"])
            local = int(report["# local misassemblies"])

            expected = benchmark[(assembler, technology, sample)]
            if extensive != expected:
                raise SystemExit(
                    f"{assembler} {technology} {sample}: report.tsv # misassemblies={extensive} "
                    f"but assembly_benchmark_30x.tsv misassemblies={expected}"
                )

            rows.append({
                "sample": sample,
                "assembler": assembler_name,
                "input": input_name,
                "extensive_misassemblies": extensive,
                "local_misassemblies": local,
                "combined_misassemblies": extensive + local,
                "source_file": str(report_path.relative_to(PROJECT_ROOT)),
            })

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"{'sample':<6} {'assembler':<9} {'input':<9} {'extensive':>9} {'local':>7} {'combined':>9}")
    for row in rows:
        print(f"{row['sample']:<6} {row['assembler']:<9} {row['input']:<9} "
              f"{row['extensive_misassemblies']:>9} {row['local_misassemblies']:>7} "
              f"{row['combined_misassemblies']:>9}")
    print(f"Wrote {len(rows)} rows: {OUTPUT.relative_to(PROJECT_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
