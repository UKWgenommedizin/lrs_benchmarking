#!/usr/bin/env python3
"""FASTA-level comparability statistics for the 30x assembler benchmark.

Run on the server where the assembly FASTAs live (the same paths QUAST read,
see assemblers/whole_genome_asm/assessment/assembly_quality_quast.smk).
Nothing is overwritten: original assemblies are only read.

For every sample x assembler x technology it computes, from the FASTA itself:

  level = scaffold        the FASTA exactly as QUAST evaluated it
  level = contig_gapsplit sequences split at every run of N/n, empty pieces
                          dropped (GoldRush: also written as a derived FASTA)

and, where the files exist:

  level = haplotype1 / haplotype2   Verkko assembly.haplotype{1,2}.fasta
  level = header_prefix:<prefix>    Verkko contigs grouped by the name prefix
                                    Verkko itself assigned (haplotype1-,
                                    haplotype2-, unassigned-, ...)

Length statistics are reported twice: over all sequences, and over sequences
>= 500 bp, which matches QUAST's --min-contig 500 in the benchmark run.

For Flye it additionally summarises assembly_info.txt (alt_group column and
coverage) to test for retained alternative/haplotig sequence.

Usage (server, from the repository root):
    python assembly_analysis/scripts/metrics/fasta_comparability_stats.py
    python assembly_analysis/scripts/metrics/fasta_comparability_stats.py \\
        --gapsplit-dir assembly_quality/derived/comparability/goldrush_gapsplit
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

# Same roots as assembly_quality_quast.smk -- the FASTAs QUAST actually read.
YU_ROOT = Path("/data/genmedbfx/yu_j/lrs_benchmarking")
STOIBER_ROOT = Path("/home/stoiber_l/smbshare/lrs_benchmarking")
ASSEMBLY_ROOTS = {
    "flye": YU_ROOT / "assemblies" / "flye",
    "goldrush": STOIBER_ROOT / "goldrush",
    "verkko": STOIBER_ROOT / "verkko",
}

SAMPLES = ["HG002", "HG003", "HG004"]
MIN_CONTIG = 500  # QUAST --min-contig in the benchmark run

DEFAULT_OUT_DIR = PROJECT_ROOT / "assembly_analysis" / "tables" / "30x" / "derived" / "comparability"
DEFAULT_GAPSPLIT_DIR = PROJECT_ROOT / "assembly_quality" / "derived" / "comparability" / "goldrush_gapsplit"

NON_N = re.compile(r"[^Nn]+")

STAT_FIELDS = [
    "sample", "assembler", "technology", "level", "source_path",
    "seq_count_all", "seq_count_ge500", "total_length_all", "total_length_ge500",
    "non_n_length", "n_count", "n_per_100kb", "n_gap_runs",
    "n50_ge500", "n90_ge500", "max_length", "n50_all", "n90_all",
]


def datasets():
    """(sample, assembler, technology, fasta_path) for the 15 benchmark runs."""
    for sample in SAMPLES:
        for assembler in ("flye", "goldrush"):
            for tech, token in (("ONT", "ont"), ("HiFi", "pb")):
                yield sample, assembler, tech, ASSEMBLY_ROOTS[assembler] / f"{sample}.{token}.30x" / "assembly.fasta"
        yield sample, "verkko", "ONT+HiFi", ASSEMBLY_ROOTS["verkko"] / sample / "assembly.fasta"


def read_fasta(path: Path):
    name, chunks = None, []
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].strip(), []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def nx(lengths: list[int], fraction: float) -> int | str:
    if not lengths:
        return "NA"
    ordered = sorted(lengths, reverse=True)
    target = sum(ordered) * fraction
    running = 0
    for length in ordered:
        running += length
        if running >= target:
            return length
    return ordered[-1]


def summarise(lengths: list[int], non_n: int, n_count: int, gap_runs: int) -> dict:
    kept = [length for length in lengths if length >= MIN_CONTIG]
    total_all = sum(lengths)
    return {
        "seq_count_all": len(lengths),
        "seq_count_ge500": len(kept),
        "total_length_all": total_all,
        "total_length_ge500": sum(kept),
        "non_n_length": non_n,
        "n_count": n_count,
        "n_per_100kb": f"{n_count / total_all * 1e5:.2f}" if total_all else "NA",
        "n_gap_runs": gap_runs,
        "n50_ge500": nx(kept, 0.5),
        "n90_ge500": nx(kept, 0.9),
        "max_length": max(lengths) if lengths else "NA",
        "n50_all": nx(lengths, 0.5),
        "n90_all": nx(lengths, 0.9),
    }


def wrap(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def fasta_levels(path: Path, gapsplit_out: Path | None, verkko_prefixes: bool):
    """Scaffold-level and gap-split statistics in one streaming pass."""
    scaffold_lengths, piece_lengths = [], []
    n_total, gap_runs_total = 0, 0
    prefixes: dict[str, list[int]] = {}
    out = None
    if gapsplit_out is not None:
        gapsplit_out.parent.mkdir(parents=True, exist_ok=True)
        out = gapsplit_out.open("w")
    try:
        for name, seq in read_fasta(path):
            scaffold_lengths.append(len(seq))
            n_count = seq.count("N") + seq.count("n")
            n_total += n_count
            pieces = [(m.start(), m.end()) for m in NON_N.finditer(seq)]
            # gap runs = interior N runs between non-N pieces (+ terminal runs, if any)
            if n_count:
                gap_runs_total += len(re.findall(r"[Nn]+", seq))
            seq_id = name.split()[0]
            for index, (start, end) in enumerate(pieces, 1):
                piece_lengths.append(end - start)
                if out is not None:
                    out.write(f">{seq_id}_piece{index} {seq_id}:{start + 1}-{end}\n{wrap(seq[start:end])}\n")
            if verkko_prefixes:
                prefix = seq_id.split("-", 1)[0] if "-" in seq_id else seq_id
                prefixes.setdefault(prefix, []).append(len(seq))
    finally:
        if out is not None:
            out.close()

    non_n = sum(piece_lengths)
    yield "scaffold", str(path), summarise(scaffold_lengths, non_n, n_total, gap_runs_total)
    yield "contig_gapsplit", str(gapsplit_out or path), summarise(piece_lengths, non_n, 0, 0)
    for prefix, lengths in sorted(prefixes.items()):
        yield f"header_prefix:{prefix}", str(path), summarise(lengths, sum(lengths), 0, 0)


def flye_info_summary(sample: str, tech: str, fasta: Path) -> dict:
    info = fasta.parent / "assembly_info.txt"
    row = {"sample": sample, "technology": tech, "assembly_info": str(info)}
    if not info.is_file():
        row["status"] = "assembly_info.txt not found"
        return row
    records = []
    with info.open() as handle:
        header = handle.readline().lstrip("#").split("\t")
        header = [column.strip() for column in header]
        for line in handle:
            fields = dict(zip(header, line.rstrip("\n").split("\t")))
            records.append(fields)
    lengths = [int(r["length"]) for r in records]
    covs = [float(r["cov."]) for r in records]
    # length-weighted median coverage
    pairs = sorted(zip(covs, lengths))
    half, running, weighted_median = sum(lengths) / 2, 0, None
    for cov, length in pairs:
        running += length
        if running >= half:
            weighted_median = cov
            break
    alt = [r for r in records if r.get("alt_group", "*") not in ("*", "")]
    low_cov = [l for c, l in zip(covs, lengths) if weighted_median and c < 0.75 * weighted_median]
    row.update({
        "status": "ok",
        "contigs": len(records),
        "total_length": sum(lengths),
        "alt_group_column_present": "alt_group" in header,
        "alt_group_contigs": len(alt),
        "alt_group_length": sum(int(r["length"]) for r in alt),
        "length_weighted_median_cov": weighted_median,
        "contigs_cov_lt_0.75x_median": len(low_cov),
        "length_cov_lt_0.75x_median": sum(low_cov),
        "circular_contigs": sum(r.get("circ.") == "Y" for r in records),
    })
    log = fasta.parent / "flye.log"
    if log.is_file():
        text = log.read_text(errors="replace")
        row["log_keep_haplotypes"] = "keep-haplotypes" in text or "keep_haplotypes" in text
        cmd = [l for l in text.splitlines() if "flye" in l and "--" in l and ("hifi" in l or "nano" in l)]
        row["log_command_line"] = cmd[0].strip() if cmd else "NA"
    return row


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--gapsplit-dir", type=Path, default=DEFAULT_GAPSPLIT_DIR)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    stats_rows, flye_rows, missing = [], [], []
    for sample, assembler, tech, fasta in datasets():
        if not fasta.is_file():
            missing.append(str(fasta))
            print(f"MISSING {fasta}")
            continue
        print(f"reading {assembler} {tech} {sample}: {fasta}")
        gapsplit = None
        if assembler == "goldrush":
            gapsplit = args.gapsplit_dir / f"{sample}.{'ont' if tech == 'ONT' else 'pb'}.30x" / "assembly.gapsplit.fasta"
        for level, source, stats in fasta_levels(fasta, gapsplit, verkko_prefixes=(assembler == "verkko")):
            stats_rows.append({"sample": sample, "assembler": assembler, "technology": tech,
                               "level": level, "source_path": source, **stats})

        if assembler == "verkko":
            work = fasta.parent / "work"
            for hap in ("haplotype1", "haplotype2"):
                hap_fasta = work / f"assembly.{hap}.fasta"
                if hap_fasta.is_file():
                    for level, source, stats in fasta_levels(hap_fasta, None, verkko_prefixes=False):
                        if level == "scaffold":
                            stats_rows.append({"sample": sample, "assembler": assembler, "technology": tech,
                                               "level": hap, "source_path": source, **stats})
                else:
                    print(f"  no {hap_fasta} (expected without Hi-C/trio/Pore-C phasing input)")
        if assembler == "flye":
            flye_rows.append(flye_info_summary(sample, tech, fasta))

    stats_path = args.out_dir / "fasta_comparability_stats.tsv"
    with stats_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=STAT_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(stats_rows)
    print(f"wrote {stats_path}")

    if flye_rows:
        fields = sorted({key for row in flye_rows for key in row}, key=lambda k: (k not in ("sample", "technology"), k))
        flye_path = args.out_dir / "flye_assembly_info_summary.tsv"
        with flye_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", restval="NA")
            writer.writeheader()
            writer.writerows(flye_rows)
        print(f"wrote {flye_path}")

    if missing:
        print(f"{len(missing)} FASTA(s) missing -- their rows are absent, not zero.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
