# 30x alignment tables

This directory is organized by the role of each table. Start with the file in
`final/` for routine analysis.

## Canonical table

- `final/alignment_benchmark_30x.tsv` — primary 30x benchmarking table. It has
  24 rows and combines alignment statistics with the best available runtime,
  CPU, RAM, and indel-event measurements. In this table, `insertions` and
  `deletions` mean the number of insertion and deletion events calculated from
  samtools `ID` records; they are not inserted/deleted base totals.

View its performance columns in a terminal:

```bash
cut -f1-5,30-34,40-41 final/alignment_benchmark_30x.tsv \
  | column -t -s $'\t' \
  | less -S
```

View the indel event columns:

```bash
cut -f1-5,18-19 final/alignment_benchmark_30x.tsv \
  | column -t -s $'\t' \
  | less -S
```

## Directory guide

- `final/` — approved tables intended for reporting and routine analysis.
- `source/` — source and active intermediate tables consumed by scripts.
- `derived/` — tables generated for plots, diagnostics, and statistical tests.
- `archive/` — superseded tables and backups retained for provenance.

## Important supporting tables

- `source/alignment_summary_30x_Samtools_Christian.tsv` — core samtools
  alignment statistics used by many plotting scripts.
- `source/report_*_30x.tsv` — the tables printed in the F2 thesis report
  (input data, thread-hours, alignment summary, CIGAR-yield/error correlation
  diagnostics) plus `report_per_run_values_30x.tsv` with the unrounded per-run
  values behind them. Rebuild with
  `python3 alignment_analysis/scripts/30x/pipeline/build_report_tables_30x.py`.
- `source/mm2.run_metrics.tsv` — raw minimap2 runtime records.
- `source/alignment_summary_30x_combined_metrics_without_ram.tsv` — input to
  the runtime/RAM merge script.
- `source/alignment_summary_30x_combined_trusted_metrics.tsv` — trusted metric
  subset used by the MAPQ analysis.
- `derived/alignment_metrics_30x.tsv` — wider alignment-quality table,
  including indel and coverage fields where available.
- `derived/alignment_benchmark_30x_indel_recovery.tsv` — detailed recovery
  table containing event counts, affected-base counts, and normalized rates.
- `derived/alignment_benchmark_30x_indel_recovery_report.tsv` — provenance and
  recovery status for each benchmark row.

## Notes

- Files in `archive/` are not current analysis inputs.
- A filename containing `.backup` or `.bak` is historical and should not be
  used for new analysis.
- Measured peak RAM is not available for every aligner in the canonical table;
  `command_ram_limit_gb` is a configured limit, not measured usage.
- Runtime values are populated for all 24 combinations. The
  `runtime_configuration` and `wallclock_runtime_source_log` columns preserve
  the provenance of each value and are exported with the runtime plot data.
- Indel event counts are available for all 24 configurations. The full
  ONT minimap2/pbmm2 stats files are kept outside the repository and are
  passed as a second `--stats-dir` to `quality_check_aligners_indels.py`.
  The HG004 ONT VG stats file in `statistics_cram_files/` was replaced by the
  complete copy from `samtools_stats_30x_Christian/vacmap_vg_stats/`; the
  earlier copy had been truncated at 216 MiB, before its `ID` records.

After regenerating the canonical runtime table, restore its recovered indel
event columns with:

```bash
python alignment_analysis/scripts/30x/pipeline/merge_indels_into_benchmark.py
```
