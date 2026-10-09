# Final report workspace

This directory contains the F2 report and workflow material used to support it.

## Report

| File | Role |
|---|---|
| [`final_f2_report.tex`](final_f2_report.tex) | Canonical LaTeX source |
| [`commands_latex_files.tex`](commands_latex_files.tex) | Supporting LaTeX command/reference material |
| [`final_f2_report.pdf`](final_f2_report.pdf) | Compiled report |

LaTeX auxiliary files (`.aux`, `.fls`, `.fdb_latexmk`, `.out`, `.toc`) are generated and ignored by Git.

## Table provenance

The thesis report (`../../f2-thesis-report/final_f2_report.tex`) reads its tables from TSV files. How each one was produced:

| Report table | TSV in `f2-thesis-report/tables/` | Produced by (in this repo) | Pipeline output |
|---|---|---|---|
| Table S1 `tab:input_data` | `input_sequencing_data_30x.tsv` | no generator script; added by hand in f2-thesis-report commit `38180fb` | — |
| Table 1 `tab:alignment_summary` | `alignment_summary_table_30x.tsv` | `alignment_analysis/scripts/30x/pipeline/build_summary_table_30x.py` (input: `alignment_analysis/tables/30x/final/alignment_benchmark_30x.tsv`) | `alignment_analysis/tables/30x/derived/alignment_summary_table_30x.tsv` |
| Table S2 `tab:thread_hours` | `thread_hours_30x.tsv` | `final_report_files/scripts/make_thread_hours_table.py` | — |
| not in the report (TSV kept; previously "Table S1") `tab:supp_cigar_correlation_diagnostics` | `cigar_yield_error_correlation_diagnostics.tsv` | `alignment_analysis/scripts/30x/plots/exploratory/normalize_cigar_yield_error.py` (input: `alignment_analysis/tables/30x/source/alignment_summary_30x_Samtools_Christian.tsv`) | `alignment_analysis/tables/30x/derived/plot_data/cigar_yield_error_correlation_diagnostics.tsv` |

Figure 1e (Q20/Q30) and 1f (read N50) use every 30x input read: `alignment_analysis/scripts/30x/pipeline/build_read_qc_30x.py` reads the RL and FFQ histograms of the full samtools stats of the 30x CRAMs that keep all input reads (minimap2 for PacBio, VG Giraffe for ONT), checks that they total exactly the raw input reads/bases, cross-checks a second independent CRAM, and writes `alignment_analysis/tables/30x/derived/read_qc_30x.tsv`. Q20/Q30 are per-base percentages; per-read mean quality would need the 30x FASTQs on the server. The earlier 1,000-read-subset values (`alignment_analysis/tables/fastq_summary.tsv`) are no longer used by the report.

### Independent check against raw data

[`scripts/check_report_tables_against_raw_data.py`](scripts/check_report_tables_against_raw_data.py) recomputes all three tables directly from raw files (samtools stats SN files in `../samtools_stats_30x_Christian/`, `../run_metrics/mm2.run_metrics.tsv`, pbmm2 logs and VG/VACmap `map_sort` logs in `../alignment_analysis/logs/`) without using any intermediate summary table, and compares every cell at the precision printed in the PDF.

```bash
python3 final_report_files/scripts/check_report_tables_against_raw_data.py   # exit 1 on any MISMATCH
```

Outputs in [`tables/report_verification/`](tables/report_verification/):

| File | Content |
|---|---|
| `report_tables_check.tsv` | One row per report cell: report value, recomputed value, verdict |
| `per_run_raw_metrics.tsv` | The 24 runs: raw counts, derived metrics, runtime/RAM and the log it came from |
| `table1_input_sequencing_data_recomputed.tsv` | Table 1 from raw data |
| `table2_alignment_summary_recomputed.tsv` | Table 2 from raw data |
| `tableS1_cigar_yield_error_correlation_recomputed.tsv` | Table S1 from raw data (unrounded) |

Status (2026-09-24): no MISMATCH. Points to disclose or fix in the report:

- **Table 2 runtime/RAM mix presets** for 3 runs because the matching log does not exist: PacBio minimap2 HG003 and HG004 use the `mm2-ont` run (ONT preset on HiFi reads), PacBio pbmm2 HG002 uses the `pbmm2-ont` log. All other Table 2 metrics come from the correct preset.
- **Table 2 threads**: minimap2 "64" is the Snakemake rule total (48 mapping + 16 sorting threads); minimap2 and VG thread counts are not recorded in the raw logs (VG: `threads: 16` in `*.read_mapping.vg.smk`).
- **Table 2 runtime definitions differ**: minimap2 = minimap2 real time, pbmm2 = pbmm2 "Run Time" (minute resolution), VACmap/VG = START-END of the whole map+sort job.
- **Table S1 PacBio R²**: raw data give 0.16451 (prints 0.165); the report prints 0.164 because the pipeline used error rates rounded to 4 decimals.

[`scripts/Raw_data_metrics.py`](scripts/Raw_data_metrics.py) is the earlier generic checker that searched the repository for the rows of the old template Table 1 in `final_f2_report.tex` (outputs `tables/report_table1_extracted.tsv`, `tables/report_tables_provenance.tsv`) and showed it was placeholder data.

## Workflows

- [`snakemake_aligners_benchmarking/`](snakemake_aligners_benchmarking/README.md): integrated alignment workflow, active configuration, sample sheet, environments, and validation.
- `snakemake_assemblers_benchmarking/`: historical/working assembly workflow tree. The maintained assembly entry point is documented under [`../assemblers/`](../assemblers/README.md).
- Historical Snakemake training projects have been moved to [`../archive/tutorials/`](../archive/tutorials/) so they cannot be mistaken for production workflows.

The report currently has no active `\includegraphics` statements; figure inclusion should use paths documented in the relevant analysis module.
