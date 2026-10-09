#!/usr/bin/env python3
"""Build 30x input read QC (read N50, base-quality Q20/Q30) from samtools stats.

The whole-genome 30x FASTQs are only on the server, but the full samtools
stats of the 30x CRAMs hold two complete input histograms:

- RL  (read length -> read count): gives read count, total bases, read N50
  and maximum length exactly.
- FFQ (cycle x Phred quality -> base count): gives the share of all input
  bases with quality >= Q20 / >= Q30 and the mean base quality.

These describe the 30x input only when the CRAM kept every input read, so
for each sample/technology the source is a CRAM from an aligner that
retains unmapped reads (minimap2 for PacBio, VG Giraffe for ONT; VACmap
drops unmapped reads and is never used). A stats file is accepted only if

- its RL read count and RL bases equal the SN raw total sequences / total
  length and the verified raw input reads/bases of that dataset, and
- its FFQ base total equals the same raw input bases.

When a second, independent CRAM of the same input is available (VG for
PacBio, a second copy of the VG stats for ONT), its values are recomputed
and must be identical, otherwise the script stops.

samtools stats records quality per base, not per read, so Q20/Q30 here
are percentages of *bases*; a per-read mean-quality Q20/Q30 needs the
FASTQs on the server.

Output -> alignment_analysis/tables/30x/derived/read_qc_30x.tsv

Usage:
    python alignment_analysis/scripts/30x/pipeline/build_read_qc_30x.py
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

PROJECT = next(parent for parent in Path(__file__).resolve().parents if (parent / "alignment_analysis").is_dir())

STATS_DIR = PROJECT / "alignment_analysis" / "statistics_cram_files"
VG_VACMAP_STATS_DIR = PROJECT / "samtools_stats_30x_Christian" / "vacmap_vg_stats"
PLOT_DATA_DIR = PROJECT / "alignment_analysis" / "tables" / "30x" / "derived" / "plot_data"
RAW_READS_TSV = PLOT_DATA_DIR / "raw_read_counts_30x.tsv"
RAW_BASES_TSV = PLOT_DATA_DIR / "raw_base_counts_30x.tsv"
OUTPUT_TSV = PROJECT / "alignment_analysis" / "tables" / "30x" / "derived" / "read_qc_30x.tsv"

SAMPLES = ["HG002", "HG003", "HG004"]

# read_technology -> (file-name code, candidate stats files in priority order)
SOURCES = {
    "ONT": [
        VG_VACMAP_STATS_DIR / "{sample}.ont.30x.hg38.vg-ont.cram.stats",
        STATS_DIR / "{sample}.ont.30x.hg38.vg-ont.cram.stats",
    ],
    "PacBio": [
        STATS_DIR / "{sample}.pb.30x.hg38.mm2-pb.cram.stats",
        STATS_DIR / "{sample}.pb.30x.hg38.vg-pb.cram.stats",
    ],
}


def parse_stats(path: Path) -> dict:
    """Stream one samtools stats file, keeping only SN, RL and FFQ totals."""
    sn: dict[str, float] = {}
    read_lengths: dict[int, int] = {}
    quality_bases: dict[int, int] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("SN\t"):
                fields = line.rstrip("\n").split("\t")
                sn[fields[1].rstrip(":")] = float(fields[2])
            elif line.startswith("RL\t"):
                _, length, count = line.rstrip("\n").split("\t")[:3]
                read_lengths[int(length)] = read_lengths.get(int(length), 0) + int(count)
            elif line.startswith("FFQ\t"):
                counts = line.rstrip("\n").split("\t")[2:]
                for quality, count in enumerate(counts):
                    if count and count != "0":
                        quality_bases[quality] = quality_bases.get(quality, 0) + int(count)
    return {"sn": sn, "rl": read_lengths, "ffq": quality_bases}


def summarise(parsed: dict) -> dict:
    rl, ffq = parsed["rl"], parsed["ffq"]
    reads = sum(rl.values())
    bases = sum(length * count for length, count in rl.items())

    cumulative = 0
    read_n50 = None
    for length in sorted(rl, reverse=True):
        cumulative += length * rl[length]
        if cumulative * 2 >= bases:
            read_n50 = length
            break

    quality_total = sum(ffq.values())
    return {
        "reads": reads,
        "total_bases": bases,
        "read_n50": read_n50,
        "max_read_length": max(rl) if rl else None,
        "quality_bases": quality_total,
        "q20_bases_percent": 100.0 * sum(c for q, c in ffq.items() if q >= 20) / quality_total if quality_total else None,
        "q30_bases_percent": 100.0 * sum(c for q, c in ffq.items() if q >= 30) / quality_total if quality_total else None,
        "mean_base_quality": sum(q * c for q, c in ffq.items()) / quality_total if quality_total else None,
    }


def complete(parsed: dict, summary: dict, raw_reads: int, raw_bases: int) -> tuple[bool, str]:
    sn = parsed["sn"]
    checks = {
        "RL reads == SN raw total sequences": summary["reads"] == sn.get("raw total sequences"),
        "RL bases == SN total length": summary["total_bases"] == sn.get("total length"),
        "RL reads == raw input reads": summary["reads"] == raw_reads,
        "RL bases == raw input bases": summary["total_bases"] == raw_bases,
        "FFQ bases == raw input bases": summary["quality_bases"] == raw_bases,
    }
    failed = [name for name, ok in checks.items() if not ok]
    return (not failed, "; ".join(failed) if failed else "complete")


def main() -> int:
    raw_reads = pd.read_csv(RAW_READS_TSV, sep="\t").set_index(["sample", "read_technology"])["raw_input_reads"]
    raw_bases = pd.read_csv(RAW_BASES_TSV, sep="\t").set_index(["sample", "read_technology"])["raw_input_bases"]

    rows = []
    for technology, templates in SOURCES.items():
        for sample in SAMPLES:
            expected_reads = int(raw_reads[(sample, technology)])
            expected_bases = int(raw_bases[(sample, technology)])
            accepted = []
            for template in templates:
                path = Path(str(template).format(sample=sample))
                if not path.is_file():
                    print(f"  {sample} {technology}: missing {path.relative_to(PROJECT)}")
                    continue
                parsed = parse_stats(path)
                summary = summarise(parsed)
                ok, detail = complete(parsed, summary, expected_reads, expected_bases)
                print(f"  {sample} {technology}: {path.relative_to(PROJECT)} -> {detail}")
                if ok:
                    accepted.append((path, summary))

            if not accepted:
                raise SystemExit(f"No complete samtools stats file for {sample} {technology}.")

            source_path, summary = accepted[0]
            crosscheck = "no second complete file"
            if len(accepted) > 1:
                other_path, other = accepted[1]
                keys = ["reads", "total_bases", "read_n50", "max_read_length",
                        "q20_bases_percent", "q30_bases_percent", "mean_base_quality"]
                differing = [key for key in keys if other[key] != summary[key]]
                if differing:
                    raise SystemExit(
                        f"{sample} {technology}: {source_path.name} and {other_path.name} disagree on {differing}"
                    )
                crosscheck = f"identical to {other_path.relative_to(PROJECT)}"

            rows.append({
                "sample": sample,
                "read_technology": technology,
                "reads": summary["reads"],
                "total_bases": summary["total_bases"],
                "mean_read_length": summary["total_bases"] / summary["reads"],
                "read_n50": summary["read_n50"],
                "read_n50_kb": summary["read_n50"] / 1000,
                "max_read_length": summary["max_read_length"],
                "q20_bases_percent": summary["q20_bases_percent"],
                "q30_bases_percent": summary["q30_bases_percent"],
                "mean_base_quality": summary["mean_base_quality"],
                "stats_file": source_path.relative_to(PROJECT).as_posix(),
                "crosscheck": crosscheck,
            })

    table = pd.DataFrame(rows)
    OUTPUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(OUTPUT_TSV, sep="\t", index=False)
    print()
    print(table.drop(columns=["stats_file", "crosscheck"]).to_string(index=False))
    print(f"\nWrote: {OUTPUT_TSV.relative_to(PROJECT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
