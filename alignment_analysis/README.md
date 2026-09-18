# Long-Read Alignment Benchmarking

This module contains the analysis layer for the long-read alignment benchmark.

The active mapper workflows remain at the repository root in accordance with the project constitution. This directory contains metric extraction, result tables, figures, and analysis documentation.

## Repository explorer

<!-- AUTO_REPOSITORY_TREE_START -->
Generated from Git-tracked files. Expand only the directory you need. On GitHub, press **`t`** for fast filename search.

- [`COMMAND_LOG.md`](COMMAND_LOG.md)
- [`README.md`](README.md)
- [`README_assemblers.md`](README_assemblers.md)

<details>
<summary><b>figures/</b> — 44 files</summary>

- [`.gitignore`](figures/.gitignore)
- [`absolute_mapped_bases_by_sample.pdf`](figures/absolute_mapped_bases_by_sample.pdf)
- [`absolute_mapped_bases_by_sample.png`](figures/absolute_mapped_bases_by_sample.png)
- [`alignment_error_rate.pdf`](figures/alignment_error_rate.pdf)
- [`alignment_error_rate.png`](figures/alignment_error_rate.png)
- [`alignment_error_rate_by_technology.pdf`](figures/alignment_error_rate_by_technology.pdf)
- [`alignment_error_rate_by_technology.png`](figures/alignment_error_rate_by_technology.png)
- [`alignment_error_rate_comparison.pdf`](figures/alignment_error_rate_comparison.pdf)
- [`alignment_error_rate_comparison.png`](figures/alignment_error_rate_comparison.png)
- [`alignment_error_rate_four_configurations.pdf`](figures/alignment_error_rate_four_configurations.pdf)
- [`alignment_error_rate_four_configurations.png`](figures/alignment_error_rate_four_configurations.png)
- [`alignment_error_rate_per_file.pdf`](figures/alignment_error_rate_per_file.pdf)
- [`alignment_error_rate_per_file.png`](figures/alignment_error_rate_per_file.png)
- [`alignment_sequences_per_file.pdf`](figures/alignment_sequences_per_file.pdf)
- [`alignment_sequences_per_file.png`](figures/alignment_sequences_per_file.png)
- [`bases_mapped_cigar_by_technology.pdf`](figures/bases_mapped_cigar_by_technology.pdf)
- [`bases_mapped_cigar_by_technology.png`](figures/bases_mapped_cigar_by_technology.png)
- [`cigar_mapped_percent_by_technology.pdf`](figures/cigar_mapped_percent_by_technology.pdf)
- [`cigar_mapped_percent_by_technology.png`](figures/cigar_mapped_percent_by_technology.png)
- [`cigar_mapped_percent_diagnostic.pdf`](figures/cigar_mapped_percent_diagnostic.pdf)
- [`cigar_mapped_percent_diagnostic.png`](figures/cigar_mapped_percent_diagnostic.png)
- [`final_alignment_error_rate.pdf`](figures/final_alignment_error_rate.pdf)
- [`final_alignment_error_rate.png`](figures/final_alignment_error_rate.png)
- [`general_technology_comparison.pdf`](figures/general_technology_comparison.pdf)
- [`general_technology_comparison.png`](figures/general_technology_comparison.png)
- [`mapped_bases_by_configuration.pdf`](figures/mapped_bases_by_configuration.pdf)
- [`mapped_bases_by_configuration.png`](figures/mapped_bases_by_configuration.png)
- [`mapped_bases_cigar_reference_style.pdf`](figures/mapped_bases_cigar_reference_style.pdf)
- [`mapped_bases_cigar_reference_style.png`](figures/mapped_bases_cigar_reference_style.png)
- [`mapped_bases_cigar_reference_style_values.tsv`](figures/mapped_bases_cigar_reference_style_values.tsv)
- [`mapped_bases_example_style_true_values.pdf`](figures/mapped_bases_example_style_true_values.pdf)
- [`mapped_bases_example_style_true_values.png`](figures/mapped_bases_example_style_true_values.png)
- [`mapped_bases_percent_by_sample.pdf`](figures/mapped_bases_percent_by_sample.pdf)
- [`mapped_bases_percent_by_sample.png`](figures/mapped_bases_percent_by_sample.png)
- [`mapped_bases_reference_style_free_y.pdf`](figures/mapped_bases_reference_style_free_y.pdf)
- [`mapped_bases_reference_style_free_y.png`](figures/mapped_bases_reference_style_free_y.png)
- [`mapped_bases_reference_style_shared_y.pdf`](figures/mapped_bases_reference_style_shared_y.pdf)
- [`mapped_bases_reference_style_shared_y.png`](figures/mapped_bases_reference_style_shared_y.png)
- [`mapped_bases_reference_true_values.pdf`](figures/mapped_bases_reference_true_values.pdf)
- [`mapped_bases_reference_true_values.png`](figures/mapped_bases_reference_true_values.png)
- [`mean_read_length.png`](figures/mean_read_length.png)
- [`median_read_length.png`](figures/median_read_length.png)
- [`n50_comparison.png`](figures/n50_comparison.png)
- [`q20_q30_comparison.png`](figures/q20_q30_comparison.png)

</details>

<details>
<summary><b>scripts/</b> — 37 files</summary>

- [`01_count_fastq_reads.sh`](scripts/01_count_fastq_reads.sh)
- [`README.md`](scripts/README.md)
- [`bar_aligment_error_rate_per_file.py`](scripts/bar_aligment_error_rate_per_file.py)
- [`build_alignment_summary.py`](scripts/build_alignment_summary.py)
- [`fastq_length_summary.sh`](scripts/fastq_length_summary.sh)
- [`fastq_quality.awk`](scripts/fastq_quality.awk)
- [`plot_absolute_mapped_bases.py`](scripts/plot_absolute_mapped_bases.py)
- [`plot_alignment_error_rate.py`](scripts/plot_alignment_error_rate.py)
- [`plot_bases_mapped_cigar_by_technology.py`](scripts/plot_bases_mapped_cigar_by_technology.py)
- [`plot_cigar_mapped_by_technology.py`](scripts/plot_cigar_mapped_by_technology.py)
- [`plot_cigar_mapped_percent_check.py`](scripts/plot_cigar_mapped_percent_check.py)
- [`plot_final_alignment_benchmark.py`](scripts/plot_final_alignment_benchmark.py)
- [`plot_general_technology_comparison.py`](scripts/plot_general_technology_comparison.py)
- [`plot_mapped_bases_by_configuration.py`](scripts/plot_mapped_bases_by_configuration.py)
- [`plot_mapped_bases_cigar_final.py`](scripts/plot_mapped_bases_cigar_final.py)
- [`plot_mapped_bases_example_style_true_values.py`](scripts/plot_mapped_bases_example_style_true_values.py)
- [`plot_mapped_bases_percent.py`](scripts/plot_mapped_bases_percent.py)
- [`plot_mapped_bases_reference_style.py`](scripts/plot_mapped_bases_reference_style.py)
- [`plot_mapped_bases_reference_true_values.py`](scripts/plot_mapped_bases_reference_true_values.py)
- [`plot_mean_length.R`](scripts/plot_mean_length.R)
- [`plot_median_length.R`](scripts/plot_median_length.R)
- [`plot_n50.R`](scripts/plot_n50.R)
- [`plot_quality.R`](scripts/plot_quality.R)
- [`run_HG002_ont_test.sh`](scripts/run_HG002_ont_test.sh)
- [`run_HG002_pbmm2_ccs.sh`](scripts/run_HG002_pbmm2_ccs.sh)
- [`run_minimap2_cross_preset.sh`](scripts/run_minimap2_cross_preset.sh)
- [`run_one_alignment.sh`](scripts/run_one_alignment.sh)
- [`run_pbmm2_ccs_indexed.sh`](scripts/run_pbmm2_ccs_indexed.sh)
- [`run_pbmm2_hifi_indexed.sh`](scripts/run_pbmm2_hifi_indexed.sh)
- [`run_pbmm2_subread_indexed.sh`](scripts/run_pbmm2_subread_indexed.sh)
- [`setup_samples_try.sh`](scripts/setup_samples_try.sh)

<details open>
<summary><b>30x/</b> — 2 files</summary>

- [`quality_check_aligners.py`](scripts/30x/quality_check_aligners.py)
- [`quality_check_aligners_indels.py`](scripts/30x/quality_check_aligners_indels.py)

</details>

<details open>
<summary><b>legacy/</b> — 1 file</summary>

- [`.gitkeep`](scripts/legacy/.gitkeep)

</details>

<details open>
<summary><b>plots/</b> — 1 file</summary>

- [`.gitkeep`](scripts/plots/.gitkeep)

</details>

<details open>
<summary><b>utils/</b> — 1 file</summary>

- [`.gitkeep`](scripts/utils/.gitkeep)

</details>

<details open>
<summary><b>validation/</b> — 1 file</summary>

- [`.gitkeep`](scripts/validation/.gitkeep)

</details>

</details>

<details>
<summary><b>tables/</b> — 71 files</summary>

- [`HG002.ont.1k.flagstat.txt`](tables/HG002.ont.1k.flagstat.txt)
- [`HG002.ont.1k.idxstats.tsv`](tables/HG002.ont.1k.idxstats.tsv)
- [`HG002.ont.1k.mm2-pb.flagstat.txt`](tables/HG002.ont.1k.mm2-pb.flagstat.txt)
- [`HG002.ont.1k.mm2-pb.idxstats.tsv`](tables/HG002.ont.1k.mm2-pb.idxstats.tsv)
- [`HG002.ont.1k.mm2-pb.samtools_stats.txt`](tables/HG002.ont.1k.mm2-pb.samtools_stats.txt)
- [`HG002.ont.1k.pbmm2-ccs.flagstat.txt`](tables/HG002.ont.1k.pbmm2-ccs.flagstat.txt)
- [`HG002.ont.1k.pbmm2-ccs.idxstats.tsv`](tables/HG002.ont.1k.pbmm2-ccs.idxstats.tsv)
- [`HG002.ont.1k.pbmm2-subread.flagstat.txt`](tables/HG002.ont.1k.pbmm2-subread.flagstat.txt)
- [`HG002.ont.1k.pbmm2-subread.idxstats.tsv`](tables/HG002.ont.1k.pbmm2-subread.idxstats.tsv)
- [`HG002.ont.1k.samtools_stats.txt`](tables/HG002.ont.1k.samtools_stats.txt)
- [`HG002.ont.quality.tsv`](tables/HG002.ont.quality.tsv)
- [`HG002.ont.read_lengths.txt`](tables/HG002.ont.read_lengths.txt)
- [`HG002.pb.1k.flagstat.txt`](tables/HG002.pb.1k.flagstat.txt)
- [`HG002.pb.1k.idxstats.tsv`](tables/HG002.pb.1k.idxstats.tsv)
- [`HG002.pb.1k.mm2-ont.flagstat.txt`](tables/HG002.pb.1k.mm2-ont.flagstat.txt)
- [`HG002.pb.1k.mm2-ont.idxstats.tsv`](tables/HG002.pb.1k.mm2-ont.idxstats.tsv)
- [`HG002.pb.1k.mm2-ont.samtools_stats.txt`](tables/HG002.pb.1k.mm2-ont.samtools_stats.txt)
- [`HG002.pb.1k.pbmm2-ccs.flagstat.txt`](tables/HG002.pb.1k.pbmm2-ccs.flagstat.txt)
- [`HG002.pb.1k.pbmm2-ccs.idxstats.tsv`](tables/HG002.pb.1k.pbmm2-ccs.idxstats.tsv)
- [`HG002.pb.1k.pbmm2-ccs.samtools_stats.txt`](tables/HG002.pb.1k.pbmm2-ccs.samtools_stats.txt)
- [`HG002.pb.1k.pbmm2-subread.flagstat.txt`](tables/HG002.pb.1k.pbmm2-subread.flagstat.txt)
- [`HG002.pb.1k.pbmm2-subread.idxstats.tsv`](tables/HG002.pb.1k.pbmm2-subread.idxstats.tsv)
- [`HG002.pb.1k.samtools_stats.txt`](tables/HG002.pb.1k.samtools_stats.txt)
- [`HG003.ont.1k.flagstat.txt`](tables/HG003.ont.1k.flagstat.txt)
- [`HG003.ont.1k.idxstats.tsv`](tables/HG003.ont.1k.idxstats.tsv)
- [`HG003.ont.1k.mm2-pb.flagstat.txt`](tables/HG003.ont.1k.mm2-pb.flagstat.txt)
- [`HG003.ont.1k.mm2-pb.idxstats.tsv`](tables/HG003.ont.1k.mm2-pb.idxstats.tsv)
- [`HG003.ont.1k.mm2-pb.samtools_stats.txt`](tables/HG003.ont.1k.mm2-pb.samtools_stats.txt)
- [`HG003.ont.1k.pbmm2-ccs.flagstat.txt`](tables/HG003.ont.1k.pbmm2-ccs.flagstat.txt)
- [`HG003.ont.1k.pbmm2-ccs.idxstats.tsv`](tables/HG003.ont.1k.pbmm2-ccs.idxstats.tsv)
- [`HG003.ont.1k.pbmm2-subread.flagstat.txt`](tables/HG003.ont.1k.pbmm2-subread.flagstat.txt)
- [`HG003.ont.1k.pbmm2-subread.idxstats.tsv`](tables/HG003.ont.1k.pbmm2-subread.idxstats.tsv)
- [`HG003.ont.1k.samtools_stats.txt`](tables/HG003.ont.1k.samtools_stats.txt)
- [`HG003.pb.1k.flagstat.txt`](tables/HG003.pb.1k.flagstat.txt)
- [`HG003.pb.1k.idxstats.tsv`](tables/HG003.pb.1k.idxstats.tsv)
- [`HG003.pb.1k.mm2-ont.flagstat.txt`](tables/HG003.pb.1k.mm2-ont.flagstat.txt)
- [`HG003.pb.1k.mm2-ont.idxstats.tsv`](tables/HG003.pb.1k.mm2-ont.idxstats.tsv)
- [`HG003.pb.1k.mm2-ont.samtools_stats.txt`](tables/HG003.pb.1k.mm2-ont.samtools_stats.txt)
- [`HG003.pb.1k.pbmm2-ccs.flagstat.txt`](tables/HG003.pb.1k.pbmm2-ccs.flagstat.txt)
- [`HG003.pb.1k.pbmm2-ccs.idxstats.tsv`](tables/HG003.pb.1k.pbmm2-ccs.idxstats.tsv)
- [`HG003.pb.1k.pbmm2-subread.flagstat.txt`](tables/HG003.pb.1k.pbmm2-subread.flagstat.txt)
- [`HG003.pb.1k.pbmm2-subread.idxstats.tsv`](tables/HG003.pb.1k.pbmm2-subread.idxstats.tsv)
- [`HG003.pb.1k.samtools_stats.txt`](tables/HG003.pb.1k.samtools_stats.txt)
- [`HG004.ont.1k.flagstat.txt`](tables/HG004.ont.1k.flagstat.txt)
- [`HG004.ont.1k.idxstats.tsv`](tables/HG004.ont.1k.idxstats.tsv)
- [`HG004.ont.1k.mm2-pb.flagstat.txt`](tables/HG004.ont.1k.mm2-pb.flagstat.txt)
- [`HG004.ont.1k.mm2-pb.idxstats.tsv`](tables/HG004.ont.1k.mm2-pb.idxstats.tsv)
- [`HG004.ont.1k.mm2-pb.samtools_stats.txt`](tables/HG004.ont.1k.mm2-pb.samtools_stats.txt)
- [`HG004.ont.1k.pbmm2-ccs.flagstat.txt`](tables/HG004.ont.1k.pbmm2-ccs.flagstat.txt)
- [`HG004.ont.1k.pbmm2-ccs.idxstats.tsv`](tables/HG004.ont.1k.pbmm2-ccs.idxstats.tsv)
- [`HG004.ont.1k.pbmm2-subread.flagstat.txt`](tables/HG004.ont.1k.pbmm2-subread.flagstat.txt)
- [`HG004.ont.1k.pbmm2-subread.idxstats.tsv`](tables/HG004.ont.1k.pbmm2-subread.idxstats.tsv)
- [`HG004.ont.1k.samtools_stats.txt`](tables/HG004.ont.1k.samtools_stats.txt)
- [`HG004.pb.1k.flagstat.txt`](tables/HG004.pb.1k.flagstat.txt)
- [`HG004.pb.1k.idxstats.tsv`](tables/HG004.pb.1k.idxstats.tsv)
- [`HG004.pb.1k.mm2-ont.flagstat.txt`](tables/HG004.pb.1k.mm2-ont.flagstat.txt)
- [`HG004.pb.1k.mm2-ont.idxstats.tsv`](tables/HG004.pb.1k.mm2-ont.idxstats.tsv)
- [`HG004.pb.1k.mm2-ont.samtools_stats.txt`](tables/HG004.pb.1k.mm2-ont.samtools_stats.txt)
- [`HG004.pb.1k.pbmm2-ccs.flagstat.txt`](tables/HG004.pb.1k.pbmm2-ccs.flagstat.txt)
- [`HG004.pb.1k.pbmm2-ccs.idxstats.tsv`](tables/HG004.pb.1k.pbmm2-ccs.idxstats.tsv)
- [`HG004.pb.1k.pbmm2-subread.flagstat.txt`](tables/HG004.pb.1k.pbmm2-subread.flagstat.txt)
- [`HG004.pb.1k.pbmm2-subread.idxstats.tsv`](tables/HG004.pb.1k.pbmm2-subread.idxstats.tsv)
- [`HG004.pb.1k.samtools_stats.txt`](tables/HG004.pb.1k.samtools_stats.txt)
- [`alignment_summary.tsv`](tables/alignment_summary.tsv)
- [`alignment_summary_30x.tsv`](tables/alignment_summary_30x.tsv)
- [`fastq_length_summary.tsv`](tables/fastq_length_summary.tsv)
- [`fastq_quality_all.tsv`](tables/fastq_quality_all.tsv)
- [`fastq_read_counts.tsv`](tables/fastq_read_counts.tsv)
- [`fastq_sha256.txt`](tables/fastq_sha256.txt)
- [`fastq_summary.tsv`](tables/fastq_summary.tsv)
- [`software_versions.txt`](tables/software_versions.txt)

</details>

<!-- AUTO_REPOSITORY_TREE_END -->

## Experimental design

The production comparison covers:

- samples: `HG002`, `HG003`, `HG004`
- technologies: ONT and PacBio HiFi
- aligners: minimap2, pbmm2, VACMap and VG Giraffe

For the complete 30x design:

```text
3 samples x 2 technologies x 4 aligners = 24 alignments
```

## Active mapper tags

| Aligner | ONT | PacBio HiFi |
|---|---|---|
| minimap2 | `mm2-ont` | `mm2-pb` |
| pbmm2 | `pbmm2-ont` | `pbmm2-pb` |
| VACMap | `vacmap-ont` | `vacmap-pb` |
| VG Giraffe | `vg-ont` | `vg-pb` |

Legacy pbmm2 aliases may still be recognized by the metric parser, but new results should use the canonical tags above.

## Mapper workflows

The canonical mapper workflows are intentionally kept at repository root:

```text
ont.read_mapping.minimap2.smk
pb.read_mapping.minimap2.smk
ont.read_mapping.pbmm2.smk
pb.read_mapping.pbmm2.smk
ont.read_mapping.vacmap.smk
pb.read_mapping.vacmap.smk
ont.read_mapping.vg.smk
pb.read_mapping.vg.smk
```

Do not move these files without first updating the constitution, shared-header references, documentation, and all downstream path assumptions.

## Canonical metric extraction

The production 30x alignment summary is generated by:

```text
alignment_analysis/scripts/30x/quality_check_aligners.py
```

From the repository root:

```bash
python3 alignment_analysis/scripts/30x/quality_check_aligners.py --project "$PWD"
```

When the matching reference FASTA and CRAMs are available, provide the exact reference used for mapping so CRAM-derived metrics can be calculated consistently.

The main generated table is:

```text
alignment_analysis/tables/alignment_summary.tsv
```

For final 30x analyses, a coverage-filtered table may be written as:

```text
alignment_analysis/tables/alignment_summary_30x.tsv
```

Missing metrics must remain `NA`; they must never be replaced with invented zeros.

## Metric groups

The analysis can include, where data are available and scientifically comparable:

- mapped reads and unmapped reads
- mapped-read percentage
- mapped bases and CIGAR-mapped bases
- MQ0 reads
- MAPQ summaries
- secondary and supplementary alignments
- mismatches / error rate
- insertion and deletion metrics
- clipping
- coverage and breadth
- runtime, memory and threads
- tool / Docker provenance

Not every metric is necessarily available for every historical run. Completeness should therefore be checked before plotting a metric across aligners.

## Figures

Generated alignment figures belong in:

```text
alignment_analysis/figures/
```

Plotting scripts currently under `alignment_analysis/scripts/` can be migrated to `alignment_analysis/scripts/plots/` only through the guarded migration procedure documented in `docs/f2/PATH_STABILITY.md`.

The reorganization kit keeps compatibility entry points when scripts are moved so older commands do not immediately break.

## Scientific cautions

### VACMap denominator

If a VACMap stats file contains only mapped reads, a reported `100%` mapped-read percentage is not comparable with aligners whose stats include all input reads. Confirm that the denominator represents the original FASTQ read set before interpreting a 100% value as mapping efficiency.

### MAPQ

MAPQ is useful descriptively, but aligners can calibrate mapping quality differently. Do not rank mapper accuracy solely by mean or median MAPQ.

### Error-rate terminology

If mismatch/error fields are derived from `samtools stats`, document them as such. NM-derived mismatch counts should not automatically be described as pure substitution counts.

## Directory structure

```text
alignment_analysis/
├── README.md
├── COMMAND_LOG.md
├── scripts/
│   ├── README.md
│   ├── 30x/
│   │   └── quality_check_aligners.py
│   ├── plots/          # organized plotting scripts
│   ├── validation/     # optional validation runners
│   ├── utils/          # optional helpers
│   └── legacy/         # deprecated but retained scripts
├── tables/
└── figures/
```

The canonical 30x extractor stays in `scripts/30x/`; it is not moved by the reorganization.
