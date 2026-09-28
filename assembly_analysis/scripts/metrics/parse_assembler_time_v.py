#!/usr/bin/env python3
"""Parse ``/usr/bin/time -v`` output file(s) into one row of a performance
table (wall-clock time, peak RAM) for the assembler-benchmark re-run.

Used by the five *.benchmark.smk rule variants
(assemblers/whole_genome_asm/{ont,pb}.assembly.{flye2,goldrush}.benchmark.smk,
hybrid.assembly.verkko.benchmark.smk). Snakemake's own ``benchmark:``
directive is deliberately not used for this: it samples the host-side
process tree, but every assembler here runs inside ``docker run`` --
Docker's real memory usage lives in a separate cgroup that a host-side
sampler cannot see, so ``benchmark:`` would report a near-zero, meaningless
number. Wrapping the actual tool binary with ``/usr/bin/time -v`` *inside*
the container avoids that; this is the same mechanism GoldRush's own
Makefile already uses internally (``track_time=1``), reused here uniformly
for Flye and Verkko too.

One time file per call for Flye/Verkko (the whole tool invocation is
wrapped once). Multiple time files per call for GoldRush, one per internal
Makefile stage (silver_path, golden_path, goldpolish, tigmint, 5x ntLink
rounds, ...) -- GoldRush has no single end-to-end number, so this script
sums each stage's elapsed wall-clock time and takes the max of each
stage's peak RSS, which is the standard way to report a multi-stage
pipeline's total wall-clock time and peak (not summed) memory footprint.

The elapsed-time line's format, confirmed from a real GoldRush time file
in this repository, is NOT the standard GNU coreutils colon format
(``h:mm:ss`` / ``m:ss``); it is unit-letter-suffixed and space-separated,
e.g. ``0m 30.22s`` or ``7m 57.79s``. The regex below accepts that format
and, defensively, the standard colon format too, in case a different
container's ``time`` build is used for Flye/Verkko.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


def find_repo_root(start: Path) -> Path:
    start = start.resolve()
    for candidate in [start, *start.parents]:
        if (candidate / "CONSTITUTION.md").is_file():
            return candidate
    raise RuntimeError("Could not locate lrs_benchmarking repository root")


PROJECT_ROOT = find_repo_root(Path(__file__))
DEFAULT_OUTPUT = (
    PROJECT_ROOT / "assembly_analysis" / "tables" / "assembler_performance_30x.tsv"
)

ELAPSED_LABEL = re.compile(
    r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*(.+)$", re.MULTILINE
)
MAX_RSS_LABEL = re.compile(r"Maximum resident set size \(kbytes\):\s*(\d+)", re.MULTILINE)

# Unit-letter-suffixed format actually observed in this repo's time files,
# e.g. "0m 30.22s", "7m 57.79s", or (untested but consistent) "1h 5m 12.34s".
_LETTER_FORMAT = re.compile(
    r"^\s*(?:(?P<hours>\d+)h\s*)?(?:(?P<minutes>\d+)m\s*)?(?P<seconds>[\d.]+)s\s*$"
)
# Standard GNU coreutils colon format, e.g. "1:15:33" or "0:30.22", kept as
# a fallback in case a different `time` build is used for Flye/Verkko.
_COLON_FORMAT = re.compile(
    r"^\s*(?:(?P<hours>\d+):)?(?P<minutes>\d+):(?P<seconds>[\d.]+)\s*$"
)


def parse_elapsed_seconds(raw: str) -> float:
    raw = raw.strip()
    match = _LETTER_FORMAT.match(raw) or _COLON_FORMAT.match(raw)
    if not match:
        raise ValueError(f"Could not parse elapsed-time value: {raw!r}")
    parts = match.groupdict()
    hours = float(parts["hours"]) if parts["hours"] else 0.0
    minutes = float(parts["minutes"]) if parts["minutes"] else 0.0
    seconds = float(parts["seconds"])
    return hours * 3600.0 + minutes * 60.0 + seconds


def parse_time_file(path: Path) -> tuple[float, float]:
    """Return (elapsed_seconds, peak_rss_kbytes) for one /usr/bin/time -v file."""
    text = path.read_text(encoding="utf-8")

    elapsed_match = ELAPSED_LABEL.search(text)
    rss_match = MAX_RSS_LABEL.search(text)
    if not elapsed_match:
        raise ValueError(f"{path}: no 'Elapsed (wall clock) time' line found")
    if not rss_match:
        raise ValueError(f"{path}: no 'Maximum resident set size' line found")

    elapsed_seconds = parse_elapsed_seconds(elapsed_match.group(1))
    peak_rss_kbytes = float(rss_match.group(1))
    return elapsed_seconds, peak_rss_kbytes


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assembler", required=True, choices=["flye", "goldrush", "verkko"])
    parser.add_argument("--sample", required=True, help="e.g. HG002")
    parser.add_argument("--technology", required=True, choices=["ont", "pb", "hybrid"])
    parser.add_argument("--threads", required=True, type=int)
    parser.add_argument(
        "--time-file", required=True, nargs="+", type=Path,
        help="One /usr/bin/time -v output file (Flye/Verkko), or several -- "
        "one per internal stage (GoldRush).",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    missing = [str(f) for f in args.time_file if not f.is_file()]
    if missing:
        raise FileNotFoundError(f"time file(s) not found: {missing}")

    elapsed_total = 0.0
    peak_rss_max = 0.0
    for time_file in args.time_file:
        elapsed_seconds, peak_rss_kbytes = parse_time_file(time_file)
        elapsed_total += elapsed_seconds
        peak_rss_max = max(peak_rss_max, peak_rss_kbytes)
        print(
            f"  {time_file.name}: elapsed={elapsed_seconds:.2f}s "
            f"peak_rss={peak_rss_kbytes / 1_048_576:.3f} GB"
        )

    row = {
        "assembler": args.assembler,
        "sample": args.sample,
        "technology": args.technology,
        "threads": args.threads,
        "wall_clock_seconds": round(elapsed_total, 2),
        "wall_clock_hours": round(elapsed_total / 3600.0, 4),
        "peak_rss_gb": round(peak_rss_max / 1_048_576, 4),
        "n_stages_summed": len(args.time_file),
        "source_time_files": ";".join(f.name for f in args.time_file),
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    file_exists = args.output.is_file()
    with args.output.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t")
        if not file_exists:
            writer.writeheader()
        writer.writerow(row)

    print(f"\n{args.assembler} {args.sample} {args.technology}:")
    print(f"  wall_clock: {row['wall_clock_seconds']} s ({row['wall_clock_hours']} h)")
    print(f"  peak_rss:   {row['peak_rss_gb']} GB")
    print(f"  appended to: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
