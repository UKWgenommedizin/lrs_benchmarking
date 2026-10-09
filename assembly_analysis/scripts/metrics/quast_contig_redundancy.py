#!/usr/bin/env python3
"""Measure redundant (haplotig-like) sequence in each assembly from QUAST alignments.

Input -> {quast_root}/{assembler}/{dataset}/contigs_reports/all_alignments_assembly.tsv
         (QUAST's per-contig alignments to the reference; only each contig's
         chosen placement, Best_group == True, is used)

Output -> assembly_analysis/tables/30x/derived/comparability/quast_contig_redundancy_30x.tsv

For every assembly:

- reference depth: how many *different* contigs cover each aligned reference
  base. A collapsed haploid assembly covers almost every base once; a
  diploid assembly with both haplotypes covers it about twice.
- contained contigs: contigs with >= 95% of their aligned bases lying where
  a *longer* contig also aligns -- a second copy of a region the assembly
  already has (the typical haplotig pattern).
- contig N50 / NG50 with and without the contained contigs, to show what
  removing that redundancy (e.g. with purge_dups) would and would not change.
  Removing contigs never joins the remaining ones, so NG50 can only stay the
  same or drop; N50 rises because the total length shrinks.

Usage:
    python assembly_analysis/scripts/metrics/quast_contig_redundancy.py \\
        --quast-root "/mnt/d/Master Germany/F2 internship/quast/quast"
"""
from __future__ import annotations

import argparse
import csv
import re
import statistics
from collections import defaultdict
from pathlib import Path


def find_repo_root(start: Path) -> Path:
    start = start.resolve()
    for candidate in [start, *start.parents]:
        if (candidate / "CONSTITUTION.md").is_file():
            return candidate
    raise RuntimeError("Could not locate lrs_benchmarking repository root")


PROJECT_ROOT = find_repo_root(Path(__file__))
DEFAULT_QUAST_ROOT = PROJECT_ROOT / "assembly_quality" / "quast"
DEFAULT_OUTPUT = (
    PROJECT_ROOT / "assembly_analysis" / "tables" / "30x" / "derived" / "comparability"
    / "quast_contig_redundancy_30x.tsv"
)
REFERENCE_LENGTH = 3_100_885_153  # GRCh38 GIAB v3 analysis set, as reported by QUAST
CONTAINED_FRACTION = 0.95
DATASET_RE = re.compile(r"^(?P<sample>HG\d+)(?:\.(?P<technology>ont|pb))?\.?(?P<depth>\d+x)?$")


def read_alignments(path: Path):
    contig_lengths: dict[str, int] = {}
    intervals: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        next(handle)
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[0] == "CONTIG":
                contig_lengths[fields[1]] = int(fields[2])
            elif len(fields) >= 9 and fields[8] == "True":
                start, end = sorted((int(fields[0]), int(fields[1])))
                intervals[fields[4]].append((start, end + 1, fields[5]))
    return contig_lengths, intervals


def sweep(contig_lengths, intervals):
    """Reference bases per depth, and per-contig aligned / shadowed bases."""
    depth_bases: dict[int, int] = defaultdict(int)
    aligned: dict[str, int] = defaultdict(int)
    shadowed: dict[str, int] = defaultdict(int)
    for chrom_intervals in intervals.values():
        events = []
        for start, end, contig in chrom_intervals:
            events.append((start, 1, contig))
            events.append((end, -1, contig))
        events.sort(key=lambda event: (event[0], event[1]))
        active: dict[str, int] = defaultdict(int)
        previous = None
        for position, change, contig in events:
            if previous is not None and position > previous and active:
                span = position - previous
                contigs = list(active)
                depth_bases[len(contigs)] += span
                # The top-ranked contig (longest; name breaks ties) owns the
                # segment, every other contig there is shadowed by it.
                owner = max(contigs, key=lambda c: (contig_lengths[c], c))
                for c in contigs:
                    aligned[c] += span
                    if c != owner:
                        shadowed[c] += span
            active[contig] += change
            if active[contig] == 0:
                del active[contig]
            previous = position
    return depth_bases, aligned, shadowed


def nx(lengths: list[int], target: float) -> int:
    running = 0
    for length in sorted(lengths, reverse=True):
        running += length
        if running >= target:
            return length
    return 0


def summarise(assembler: str, dataset: str, path: Path) -> dict:
    contig_lengths, intervals = read_alignments(path)
    depth_bases, aligned, shadowed = sweep(contig_lengths, intervals)

    covered = sum(depth_bases.values())
    contained = [
        c for c in contig_lengths
        if aligned.get(c, 0) > 0 and shadowed.get(c, 0) >= CONTAINED_FRACTION * aligned[c]
    ]
    contained_set = set(contained)
    all_lengths = list(contig_lengths.values())
    kept_lengths = [length for c, length in contig_lengths.items() if c not in contained_set]
    total = sum(all_lengths)
    contained_length = sum(contig_lengths[c] for c in contained)

    match = DATASET_RE.match(dataset)
    return {
        "assembler": assembler,
        "sample": match.group("sample"),
        "technology": match.group("technology") or "hybrid",
        "contigs": len(all_lengths),
        "assembly_length_bp": total,
        "reference_bases_covered": covered,
        "depth1_pct": round(100 * depth_bases.get(1, 0) / covered, 2),
        "depth2_pct": round(100 * depth_bases.get(2, 0) / covered, 2),
        "depth3plus_pct": round(100 * sum(v for k, v in depth_bases.items() if k >= 3) / covered, 2),
        "contained_contigs": len(contained),
        "contained_length_bp": contained_length,
        "contained_length_pct_of_assembly": round(100 * contained_length / total, 2),
        "contained_median_length_bp": int(statistics.median(contig_lengths[c] for c in contained)) if contained else 0,
        "n50_all": nx(all_lengths, total / 2),
        "n50_without_contained": nx(kept_lengths, sum(kept_lengths) / 2),
        "ng50_all": nx(all_lengths, REFERENCE_LENGTH / 2),
        "ng50_without_contained": nx(kept_lengths, REFERENCE_LENGTH / 2),
        "length_without_contained_bp": sum(kept_lengths),
        "source_file": path.as_posix(),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--quast-root", type=Path, default=DEFAULT_QUAST_ROOT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    paths = sorted(args.quast_root.glob("*/*/contigs_reports/all_alignments_assembly.tsv"))
    paths = [p for p in paths if DATASET_RE.match(p.parent.parent.name) and ".1k" not in p.parent.parent.name]
    if not paths:
        raise SystemExit(f"No contigs_reports/all_alignments_assembly.tsv under {args.quast_root}")

    rows = []
    for path in paths:
        dataset_dir = path.parent.parent
        rows.append(summarise(dataset_dir.parent.name, dataset_dir.name, path))
        print(f"done {dataset_dir.parent.name}/{dataset_dir.name}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
