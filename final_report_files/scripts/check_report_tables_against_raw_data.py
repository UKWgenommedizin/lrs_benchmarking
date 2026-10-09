#!/usr/bin/env python3
"""Recompute every data table of the F2 thesis report from raw files and compare cell by cell.

Tables checked (TSVs read by ``f2-thesis-report/final_f2_report.tex``):
  Table 1  tab:input_data                          tables/input_sequencing_data_30x.tsv
  Table 2  tab:alignment_summary                   tables/alignment_summary_table_30x.tsv
  Table S1 tab:supp_cigar_correlation_diagnostics  tables/cigar_yield_error_correlation_diagnostics.tsv

Raw sources (nothing is read from the intermediate summary tables):
  samtools stats SN     samtools_stats_30x_Christian/*.cram.stats.SN.txt   (minimap2, pbmm2)
                        samtools_stats_30x_Christian/vacmap_vg_stats/*.cram.stats (VACmap, VG Giraffe)
  minimap2 runtime/RAM  run_metrics/mm2.run_metrics.tsv                    ("Real time", "Peak RSS")
  pbmm2 runtime/RAM     alignment_analysis/logs/Logs_files_pbmm2/*.log     ("Run Time", "Peak RSS", threads)
  VACmap/VG runtime     alignment_analysis/logs/logs_VG_VACmap/*.map_sort.log (START/END timestamps, -t)

Runtime/RAM for a run is taken from the log of the same preset. When that log does not exist the
log of the other preset of the same aligner is used (as the analysis pipeline did) and the run is
flagged SUBSTITUTED_PRESET; every Table 2 cell that depends on such a run gets verdict
MATCH_WITH_SUBSTITUTED_PRESET instead of MATCH.

Outputs in final_report_files/tables/report_verification/:
  per_run_raw_metrics.tsv          24 runs: raw counts, derived metrics, runtime/RAM and its source log
  table1_input_sequencing_data_recomputed.tsv
  table2_alignment_summary_recomputed.tsv
  tableS1_cigar_yield_error_correlation_recomputed.tsv
  report_tables_check.tsv          one row per report cell: report value, recomputed value, verdict
Verdicts: MATCH; MATCH_WITH_SUBSTITUTED_PRESET; ROUNDING_EDGE (unrounded values differ by less than
half a printed unit but round differently); UNVERIFIABLE_NOT_IN_LOGS; MISMATCH.
Exit status 1 if any cell is MISMATCH.
"""

import argparse
import csv
import re
import statistics
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path

from scipy import stats

REPO = Path(__file__).resolve().parents[2]
REPORT_DIR = REPO.parent / "f2-thesis-report"
STATS_DIR = REPO / "samtools_stats_30x_Christian"
MM2_METRICS = REPO / "run_metrics" / "mm2.run_metrics.tsv"
PBMM2_LOGS = REPO / "alignment_analysis" / "logs" / "Logs_files_pbmm2"
VG_VACMAP_LOGS = REPO / "alignment_analysis" / "logs" / "logs_VG_VACmap"
DEFAULT_OUT = REPO / "final_report_files" / "tables" / "report_verification"

GENOME_SIZE = 3.1e9
SAMPLES = ["HG002", "HG003", "HG004"]
TECHS = {"ont": "ONT", "pb": "PacBio HiFi"}
# aligner -> {tech: configuration used in the report}
ALIGNERS = {
    "minimap2": {"ont": "mm2-ont", "pb": "mm2-pb"},
    "pbmm2": {"ont": "pbmm2-ont", "pb": "pbmm2-pb"},
    "VACmap": {"ont": "vacmap-ont", "pb": "vacmap-pb"},
    "VG Giraffe": {"ont": "vg-ont", "pb": "vg-pb"},
}
OTHER_PRESET = {"mm2-ont": "mm2-pb", "mm2-pb": "mm2-ont", "pbmm2-ont": "pbmm2-pb", "pbmm2-pb": "pbmm2-ont"}

# decimals printed in the PDF
T1_PRECISION = {"reads_m": 2, "bases_gb": 2, "coverage_x": 1, "mean_read_length_kb": 1}
T2_PRECISION = {"error_percent": 2, "mapped_reads_percent": 2, "mapped_bases_percent": 2, "cigar_yield_percent": 2,
                "mq0_reads_percent": 2, "runtime_hours": 2, "peak_ram_gb": 1}
S1_PRECISION = {"n": 0, "pearson_r": 3, "pearson_p": 3, "spearman_rho": 3, "spearman_p": 4, "r_squared": 3,
                "shapiro_residual_W": 3, "shapiro_residual_p": 4}
T2_RUNTIME_COLUMNS = {"runtime_hours", "peak_ram_gb", "threads"}


# ---------------------------------------------------------------- raw readers

def read_sn(path):
    sn = {}
    for line in path.read_text().splitlines():
        if line.startswith("SN\t"):
            key, value = line.split("\t")[1:3]
            sn[key.rstrip(":")] = float(value)
    return sn


def sn_path(sample, tech, config):
    if config.startswith(("vacmap", "vg")):
        return STATS_DIR / "vacmap_vg_stats" / f"{sample}.{tech}.30x.hg38.{config}.cram.stats"
    return STATS_DIR / f"{sample}_{tech}_30x.hg38.{config}.cram.stats.SN.txt"


def read_mm2_metrics():
    runs = {}
    for line in MM2_METRICS.read_text().splitlines():
        sample, tech, _cov, config, text = line.split("\t")
        real = float(re.search(r"Real time: ([\d.]+) sec", text).group(1))
        rss = float(re.search(r"Peak RSS: ([\d.]+) GB", text).group(1))
        runs[(sample, tech, config)] = (real, rss, None, f"{MM2_METRICS.relative_to(REPO)} [{sample} {tech} {config}]")
    return runs


def parse_duration(text):
    units = {"d": 86400, "h": 3600, "m": 60, "s": 1}
    return sum(int(n) * units[u] for n, u in re.findall(r"(\d+)([dhms])", text))


def read_pbmm2_log(path):
    text = path.read_text()
    runtime = parse_duration(re.search(r"Run Time: ([^\n]+)", text).group(1))
    rss = float(re.search(r"Peak RSS: ([\d.]+) GB", text).group(1))
    threads = int(re.search(r"Using (\d+) threads", text).group(1))
    return runtime, rss, threads, str(path.relative_to(REPO))


def read_map_sort_log(path):
    text = path.read_text()
    stamp = lambda tag: datetime.fromisoformat(re.search(rf"^\[([^\]]+)\] {tag} ", text, re.M).group(1))
    threads = re.search(r"(?:^|\s)(?:-t|--threads)\s+(\d+)", text)
    return ((stamp("END") - stamp("START")).total_seconds(), None,
            int(threads.group(1)) if threads else None, str(path.relative_to(REPO)))


def runtime_for(sample, tech, config, mm2):
    """Return (seconds, peak_gb, threads, source, preset_used)."""
    if config.startswith(("vacmap", "vg")):
        return (*read_map_sort_log(VG_VACMAP_LOGS / f"{sample}.{tech}.30x.hg38.{config}.map_sort.log"), config)
    for preset in (config, OTHER_PRESET[config]):
        if config.startswith("mm2") and (sample, tech, preset) in mm2:
            return (*mm2[(sample, tech, preset)], preset)
        log = PBMM2_LOGS / f"{sample}.{tech}.30x.hg38.{preset}.cram.log"
        if config.startswith("pbmm2") and log.exists():
            return (*read_pbmm2_log(log), preset)
    return None, None, None, "no log found", None


# ---------------------------------------------------------------- recomputation

def collect_runs():
    mm2 = read_mm2_metrics()
    runs = []
    for sample in SAMPLES:
        for tech in TECHS:
            for aligner, configs in ALIGNERS.items():
                config = configs[tech]
                sn = read_sn(sn_path(sample, tech, config))
                seconds, rss, threads, source, preset_used = runtime_for(sample, tech, config, mm2)
                runs.append({
                    "sample": sample, "tech": tech, "aligner": aligner, "configuration": config,
                    "raw_total_sequences": int(sn["raw total sequences"]), "total_length": int(sn["total length"]),
                    "reads_mapped": int(sn["reads mapped"]), "bases_mapped": int(sn["bases mapped"]),
                    "bases_mapped_cigar": int(sn["bases mapped (cigar)"]), "reads_mq0": int(sn["reads MQ0"]),
                    "error_percent": 100 * sn["error rate"],
                    "runtime_seconds": seconds, "peak_ram_gb": rss, "threads_logged": threads,
                    "runtime_source": source, "runtime_preset_used": preset_used,
                    "runtime_flag": "OK" if preset_used == config else "SUBSTITUTED_PRESET",
                    "sn_file": str(sn_path(sample, tech, config).relative_to(REPO)),
                })
    # verified input = totals shared by the runs that keep unmapped reads (all but VACmap)
    for sample in SAMPLES:
        for tech in TECHS:
            group = [r for r in runs if (r["sample"], r["tech"]) == (sample, tech)]
            keep = Counter((r["raw_total_sequences"], r["total_length"]) for r in group if r["aligner"] != "VACmap")
            (reads, bases), n = keep.most_common(1)[0]
            if n != 3:
                raise SystemExit(f"{sample} {tech}: aligners keeping unmapped reads disagree on input: {keep}")
            for r in group:
                r["input_reads"], r["input_bases"] = reads, bases
                r["mapped_reads_percent"] = 100 * r["reads_mapped"] / reads
                r["mapped_bases_percent"] = 100 * r["bases_mapped"] / bases
                r["cigar_yield_percent"] = 100 * r["bases_mapped_cigar"] / bases
                r["mq0_reads_percent"] = 100 * r["reads_mq0"] / r["reads_mapped"]
                r["runtime_hours"] = r["runtime_seconds"] / 3600 if r["runtime_seconds"] else None
    return runs


def table1(runs):
    rows = []
    for tech in TECHS:
        for sample in SAMPLES:
            r = next(r for r in runs if (r["sample"], r["tech"]) == (sample, tech))
            rows.append({"sample": sample, "technology": TECHS[tech],
                         "reads_m": r["input_reads"] / 1e6, "bases_gb": r["input_bases"] / 1e9,
                         "coverage_x": r["input_bases"] / GENOME_SIZE,
                         "mean_read_length_kb": r["input_bases"] / r["input_reads"] / 1e3})
    return rows


def table2(runs):
    rows = []
    for tech in TECHS:
        for aligner in ALIGNERS:
            group = [r for r in runs if r["tech"] == tech and r["aligner"] == aligner]
            threads = sorted({r["threads_logged"] for r in group if r["threads_logged"]})
            row = {"technology": TECHS[tech], "aligner": aligner,
                   "threads": "/".join(map(str, threads)) or "not in logs",
                   "substituted_runs": ",".join(f"{r['sample']}:{r['runtime_preset_used']}" for r in group
                                                if r["runtime_flag"] != "OK") or "-"}
            for metric in T2_PRECISION:
                values = [r[metric] for r in group if r[metric] is not None]
                row[metric] = (statistics.mean(values), statistics.stdev(values)) if len(values) == 3 else None
            rows.append(row)
    return rows


def table_s1(runs):
    rows = []
    for tech in TECHS:
        group = [r for r in runs if r["tech"] == tech]
        x = [r["cigar_yield_percent"] for r in group]
        y = [r["error_percent"] for r in group]
        reg = stats.linregress(x, y)
        residuals = [yi - (reg.intercept + reg.slope * xi) for xi, yi in zip(x, y)]
        pearson, spearman, shapiro = stats.pearsonr(x, y), stats.spearmanr(x, y), stats.shapiro(residuals)
        rows.append({"technology": "PacBio" if tech == "pb" else "ONT", "n": len(x),
                     "pearson_r": pearson.statistic, "pearson_p": pearson.pvalue,
                     "spearman_rho": spearman.statistic, "spearman_p": spearman.pvalue,
                     "r_squared": reg.rvalue ** 2,
                     "shapiro_residual_W": shapiro.statistic, "shapiro_residual_p": shapiro.pvalue})
    return rows


# ---------------------------------------------------------------- comparison

def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, rows, header=None):
    header = header or list(rows[0])
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def fmt(value, decimals):
    return f"{value:.{decimals}f}"


def fmt_pm(pair, decimals):
    return "n.m." if pair is None else f"{pair[0]:.{decimals}f} ± {pair[1]:.{decimals}f}"


def compare(report_dir):
    runs = collect_runs()
    t1, t2, s1 = table1(runs), table2(runs), table_s1(runs)
    checks = []

    def add(table, key, column, reported, recomputed, verdict=None, exact=None):
        verdict = verdict or ("MATCH" if reported == recomputed else "MISMATCH")
        # same value, printed differently only because it sits on a rounding boundary
        if verdict == "MISMATCH" and exact and abs(exact[0] - exact[1]) < 0.5 * 10 ** -S1_PRECISION.get(column, 2):
            verdict = f"ROUNDING_EDGE (report {exact[0]:.6f}, raw {exact[1]:.6f})"
        checks.append({"table": table, "row": key, "column": column, "report_value": reported,
                       "recomputed_value": recomputed, "verdict": verdict})

    for rep in read_tsv(report_dir / "tables" / "input_sequencing_data_30x.tsv"):
        new = next(r for r in t1 if (r["sample"], r["technology"]) == (rep["sample"], rep["technology"]))
        for col, d in T1_PRECISION.items():
            add("Table 1", f"{rep['sample']} {rep['technology']}", col, fmt(float(rep[col]), d), fmt(new[col], d),
                exact=(float(rep[col]), new[col]))

    tech = None
    for rep in read_tsv(report_dir / "tables" / "alignment_summary_table_30x.tsv"):
        tech = rep["technology"] or tech  # the report leaves repeated technology cells empty
        new = next(r for r in t2 if (r["technology"], r["aligner"]) == (tech, rep["aligner"]))
        key = f"{tech} {rep['aligner']}"
        substituted = new["substituted_runs"] != "-"
        logged = new["threads"]
        add("Table 2", key, "threads", rep["threads"], logged,
            "UNVERIFIABLE_NOT_IN_LOGS" if logged == "not in logs" else None)
        for col, d in T2_PRECISION.items():
            reported = rep[col].replace("$\\pm$", "±")
            recomputed = fmt_pm(new[col], d)
            verdict = None
            if reported == recomputed and substituted and col in T2_RUNTIME_COLUMNS:
                verdict = f"MATCH_WITH_SUBSTITUTED_PRESET ({new['substituted_runs']})"
            add("Table 2", key, col, reported, recomputed, verdict)

    for rep in read_tsv(report_dir / "tables" / "cigar_yield_error_correlation_diagnostics.tsv"):
        new = next(r for r in s1 if r["technology"] == rep["technology"])
        for col, d in S1_PRECISION.items():
            add("Table S1", rep["technology"], col, fmt(float(rep[col]), d), fmt(new[col], d),
                exact=(float(rep[col]), new[col]))

    return runs, t1, t2, s1, checks


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--report-dir", type=Path, default=REPORT_DIR, help="thesis report repo (default: %(default)s)")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    runs, t1, t2, s1, checks = compare(args.report_dir)

    write_tsv(args.out / "per_run_raw_metrics.tsv", runs)
    write_tsv(args.out / "table1_input_sequencing_data_recomputed.tsv",
              [{k: (fmt(v, T1_PRECISION[k]) if k in T1_PRECISION else v) for k, v in r.items()} for r in t1])
    write_tsv(args.out / "table2_alignment_summary_recomputed.tsv",
              [{k: (fmt_pm(v, T2_PRECISION[k]) if k in T2_PRECISION else v) for k, v in r.items()} for r in t2])
    write_tsv(args.out / "tableS1_cigar_yield_error_correlation_recomputed.tsv", s1)
    write_tsv(args.out / "report_tables_check.tsv", checks)

    print(f"Report : {args.report_dir}")
    for table in ("Table 1", "Table 2", "Table S1"):
        verdicts = Counter(c["verdict"].split(" ")[0] for c in checks if c["table"] == table)
        print(f"{table:8s}: " + ", ".join(f"{v} {n}" for v, n in sorted(verdicts.items())))
    for c in checks:
        if c["verdict"] != "MATCH":
            print(f"  {c['table']:8s} {c['row']:24s} {c['column']:22s} report={c['report_value']:18s} "
                  f"raw={c['recomputed_value']:18s} {c['verdict']}")
    print("Runtime/RAM taken from another preset's log:")
    for r in runs:
        if r["runtime_flag"] != "OK":
            print(f"  {r['sample']} {TECHS[r['tech']]} {r['configuration']}: used {r['runtime_preset_used']} "
                  f"({r['runtime_source']})")
    print(f"Outputs: {args.out}")
    return 1 if any(c["verdict"] == "MISMATCH" for c in checks) else 0


if __name__ == "__main__":
    sys.exit(main())
