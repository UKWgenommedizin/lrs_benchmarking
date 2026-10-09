# Assembly Analysis

This module contains downstream metric aggregation, statistical analysis and figure generation for the whole-genome assembly benchmark.

It does **not** run the assemblers themselves. Production assembly workflows are under `assemblers/whole_genome_asm/`, and external assessment tools are run under `assemblers/whole_genome_asm/assessment/`.

## Repository explorer

<!-- AUTO_REPOSITORY_TREE_START -->
Generated from version-control candidates. Directories remain compact; use this module's curated sections for canonical files.

- [`README.md`](README.md)

- **figures/** — 31 files
- **scripts/** — 28 files; [`guide`](scripts/README.md)
- **tables/** — 11 files

<!-- AUTO_REPOSITORY_TREE_END -->

## Intended structure

```text
assembly_analysis/
├── README.md
├── scripts/
│   ├── README.md
│   ├── metrics/
│   │   ├── calculate_quality_metrics.py
│   │   ├── extract_quast_metrics.py
│   │   ├── extract_busco_metrics.py
│   │   └── extract_merqury_metrics.py
│   └── plots/
├── tables/
└── figures/
```

Only files that actually exist should be documented as implemented. The names above define the preferred location when those scripts are present.

## Raw assembly metric extraction

FASTA-derived metrics can include:

- total assembly length
- number of contigs / sequences
- largest and smallest sequence
- mean and median sequence length
- N50 / L50
- N90 / L90
- nucleotide counts

These descriptive metrics should be kept distinct from reference-aware metrics and k-mer/gene completeness metrics.

## External assessment metrics

Expected sources include:

```text
assembly_quality/quast/
assembly_quality/busco/
assembly_quality/merqury/
```

The final table should retain provenance so that each value can be traced back to a specific source file, assembler, sample, technology and assembly representation.

## Canonical summary table

Recommended final location:

```text
assembly_analysis/tables/assembler_metrics.tsv
```

A long/tidy or well-documented wide table is preferable to multiple disconnected hand-edited spreadsheets.

### Implemented: `extract_quast_metrics.py`

```text
assembly_analysis/scripts/metrics/extract_quast_metrics.py
```

Reads every `assembly_quality/quast/{assembler}/{dataset}/report.tsv` and
writes the long-format table above (columns: `assembler`, `dataset`,
`sample`, `technology`, `depth`, `representation`, `metric`, `value`,
`source_file`).

```bash
python assembly_analysis/scripts/metrics/extract_quast_metrics.py
```

`assembly_quality/` is gitignored raw-tool-output staging, not tracked
history -- populate it locally, or from the server layout documented in
`assemblers/whole_genome_asm/assessment/README.md`, before running this
script. Only the resulting table under `assembly_analysis/tables/` is
committed. If a new assembler directory appears under
`assembly_quality/quast/`, the script picks it up automatically; add its
representation (`collapsed` / `haplotype_resolved`) to
`ASSEMBLER_REPRESENTATION` in the script once known, otherwise it is tagged
`unknown`.

### Implemented: `build_assembly_benchmark_table.py`

```text
assembly_analysis/scripts/metrics/build_assembly_benchmark_table.py
```

Pivots `assembler_metrics.tsv` into the curated, report-ready benchmark
table, mirroring `alignment_analysis/tables/30x/final/alignment_benchmark_30x.tsv`:

```text
assembly_analysis/tables/30x/final/assembly_benchmark_30x.tsv
```

One row per assembler/sample/technology. Columns: the context fields
(`assembler`, `sample`, `technology`, `depth`, `representation`) plus the
QUAST reference-agreement and contiguity metrics from the benchmark
workflow diagram (`contigs`, `total_length_bp`, `largest_contig_bp`,
`genome_fraction_pct`, `duplication_ratio`, `misassemblies`,
`mismatches_per_100kbp`, `indels_per_100kbp`, `n50`, `ng50`, `na50`,
`nga50`, `l50`, `lg50`, `unaligned_contigs`, `unaligned_length_bp`) and a
`source_file` provenance column. Metrics not in this curated set stay
available per-row in `assembler_metrics.tsv`.

```bash
python assembly_analysis/scripts/metrics/extract_quast_metrics.py
python assembly_analysis/scripts/metrics/build_assembly_benchmark_table.py
```

Verkko's dataset directories carry no depth token in their name (e.g.
`HG002`, not `HG002.ont.30x`), so `depth` is reported as `unknown` for its
rows rather than assumed.

## Figures

Generated assembly figures belong in:

```text
assembly_analysis/figures/30x/final/
```

mirroring `alignment_analysis/figures/30x/final/`. Plotting scripts belong in:

```text
assembly_analysis/scripts/plots/
```

and share style/color helpers from `assembly_analysis/scripts/utils/plot_style.py`
(a scope-local counterpart to `alignment_analysis/scripts/utils/plot_style.py`,
reusing the same `SAMPLE_ORDER`/`SAMPLE_COLORS` values so sample colors stay
consistent across every figure in the thesis).

Recommended figure families include:

- contiguity: NG50/N50, contig count, assembly span
- reference agreement: NGA50 / misassembly-related metrics where valid
- sequence accuracy: QV / reference discordance metrics
- completeness: BUSCO and k-mer completeness
- computational resources: runtime and peak RAM
- ntLink before/after comparisons

### Implemented: `reference_agreement_panel.py`

```text
assembly_analysis/scripts/plots/reference_agreement_panel.py
```

Reads `assembly_analysis/tables/30x/final/assembly_benchmark_30x.tsv` and
writes a 2x2 composite figure (genome fraction, misassemblies, mismatches
per 100 kbp, indels per 100 kbp), one grouped bar chart per metric with bars
grouped by sample within each of the five assembly configurations
(Flye-ONT, Flye-HiFi, Goldrush-ONT, Goldrush-HiFi, Verkko-Hybrid):

```text
assembly_analysis/figures/30x/final/01_reference_agreement_panel_30x.png
assembly_analysis/figures/30x/final/01_reference_agreement_panel_30x.pdf
```

```bash
python assembly_analysis/scripts/plots/reference_agreement_panel.py
```

### Implemented: `assembler_benchmark_bars.py` (final main figure)

```text
assembly_analysis/scripts/plots/assembler_benchmark_bars.py
```

The adopted main assembler-benchmark figure, chosen after an explicit A/B
round against `assembler_benchmark_points.py` (below): grouped bars were
kept because each bar is an actual HG002/HG003/HG004 observation, never a
mean, which is defensible on its own and consistent with the rest of the
report -- not mixed with points across panels, to avoid a second visual
grammar the reader has to learn. Panels, all from
`assembly_benchmark_30x.tsv` via the shared loader
(`assembly_analysis/scripts/utils/benchmark_data.py`):

- **a** N50 (Mb) · **b** contig count · **c** genome fraction (%) ·
  **d** misassemblies · **e** mismatches / 100 kbp · **f** indels / 100 kbp.
- NGA50 (used in the earlier comparison round, see below) and wall-clock
  runtime are both excluded, not silently substituted -- NGA50 spans
  ~65 kb-27 Mb across configurations and does not read honestly on a shared
  zero-baseline linear scale (most bars render as near-invisible slivers
  next to Flye-ONT); no measured runtime source exists for the actual 30x
  run. **N50 has the same dynamic-range issue as NGA50 did** -- panel (a)
  is still dominated by Flye-ONT, with Flye-HiFi/Goldrush/Verkko compressed
  near zero. This was surfaced, not hidden, when the figure was reviewed;
  no axis change has been made without being asked for.
- N50 and contig count are QUAST's FASTA-only "general" stats, computed
  before QUAST ever aligns to the reference -- not a stand-in for a
  separate FASTA-derived table. There is no other such table anywhere in
  this repo; see `utils/benchmark_data.py`'s module docstring for how this
  was confirmed (checked against `docs/f2/archive/README_assemblers_legacy.md`
  and a repo-wide search for any independent FASTA-statistics tool).
- Verkko is set apart by x-spacing and a thin separator line (not a
  background wash or a bracket), with a 3-tier x-axis
  (SINGLE-TECHNOLOGY/HYBRID -> assembler -> ONT/HiFi/ONT+HiFi) and a dagger
  on the "Verkko" tier-2 label. No explanatory paragraph is drawn inside
  the plot area; the caption text belongs in the report/figure legend:

  > † Verkko was assembled using combined 30× ONT and 30× PacBio HiFi input
  > and is shown separately from the single-technology Flye and GoldRush
  > assemblies. Reference-based statistics should therefore be interpreted
  > in the context of its hybrid, haplotype-resolved assembly strategy.

```bash
python assembly_analysis/scripts/plots/assembler_benchmark_bars.py
```

Prints the exact value used for every sample × assembler × technology ×
metric (with its source file) before saving, and separately asserts every
plotted number matches `assembly_benchmark_30x.tsv` exactly
(`verify_plot_data_matches_source`, re-reads the source table independently
rather than trusting the value because it was derived from it once).

Outputs:

```text
assembly_analysis/figures/30x/final/assembler_benchmark_bars.png
assembly_analysis/figures/30x/final/assembler_benchmark_bars.pdf
assembly_analysis/tables/30x/derived/plot_data/assembler_benchmark_plot_data.tsv
```

Publication-refinement pass applied: "GoldRush" capitalized consistently
(was "Goldrush"); subtle horizontal gridlines added to every panel; contig
count and misassemblies (the two panels with counts in the thousands) get
comma thousands separators on the y-axis (10,000 / 20,000 / ...). "HiFi" is
kept as the tier-3 axis label rather than spelled out as "PacBio HiFi" --
that longer label collided with the adjacent column at this spacing when
tried; "HiFi" is unambiguous in this context, and "PacBio HiFi" is used in
full where there's room (the caption note above).

Assembly-FASTA provenance was independently re-verified for this pass, not
just asserted: `quast.log` for all 15 QUAST runs was checked (not just
`report.tsv`), confirming every run's input was literally
`.../{assembler}/{dataset}/assembly.fasta` in that assembler's dedicated
production output directory -- e.g. Verkko's was
`.../verkko/HG002/assembly.fasta`, not a haplotype-split, unitig, or graph
file. Same conclusion as `utils/benchmark_data.py`'s docstring, now checked
against the actual per-run logs rather than the general QUAST output
convention alone.

### QUAST-only figure set (supersedes `assembler_benchmark_bars.py` as the main figures)

QUAST is the only evaluation source. Every mark is one HG002/HG003/HG004 observation, and each plotted value is cross-checked against its raw `report.tsv` (`load_quast_values` in `utils/benchmark_data.py`). Contig N50 and contig count were dropped from the main figures. They mix contig, scaffold (GoldRush) and pooled-diploid (Verkko) values; see `tables/30x/derived/comparability/`.

| Output | Script | Content |
|---|---|---|
| `figures/30x/final/assembler_nga50_vs_misassemblies.{png,pdf}` | `scripts/plots/assembler_nga50_vs_misassemblies.py` | Main Figure A: NGA50 (log) vs misassemblies scatter |
| `figures/30x/final/assembler_quast_quality_panel.{png,pdf}` | `scripts/plots/assembler_quast_quality_panel.py` | Main Figure B: genome fraction, duplication ratio, mismatches / 100 kbp and indels / 100 kbp |
| `tables/30x/final/assembly_quast_supplementary_30x.tsv` | `scripts/metrics/build_quast_supplementary_table.py` | Full QUAST statistics, values copied verbatim (`-` becomes NA) |

Caption caveats are in each script's docstring. NGA50 uses the haploid GRCh38 length, and misassemblies are absolute counts. Both therefore behave differently for Verkko (duplication about 2.05) and Flye HiFi (about 1.33).

**NGAx curves (supplementary):** not yet generated. Only `report.tsv` was synced locally for the 30× runs, and it holds just NGA50/NGA90, not the curve. Copy each run's `aligned_stats/` directory (QUAST's `NGAx_plot.pdf`) from the server's `assembly_quality/quast/{assembler}/{dataset}/`.

### `assembler_benchmark_n50_logbars.py` / `assembler_benchmark_n50_points.py` (N50 panel-a candidates)

```text
assembly_analysis/scripts/plots/assembler_benchmark_n50_logbars.py
assembly_analysis/scripts/plots/assembler_benchmark_n50_points.py
```

N50 has the same dynamic-range problem NGA50 did (~0.14-33.6 Mb across
configurations, ~2.4 orders of magnitude): on the bars script's linear
zero-baseline panel (a), Flye-ONT dwarfs the other four configurations.
These two scripts hold panels b-f fixed -- imported directly from
`assembler_benchmark_bars.py`/`assembler_benchmark_points.py`, not
re-implemented, so they cannot silently diverge from the adopted figure --
and vary only panel (a):

- **`_logbars.py`**: grouped bars (still individual HG002/HG003/HG004
  bars, never a mean) on a log y-axis, `style_log_y_axis()` in
  `utils/plot_style.py` (major ticks at each power of ten with clean "0.1"
  /"1"/"10" labels, unlabeled minor ticks at the 2x/5x subdivisions -- a
  continuous transform, not a broken axis).
- **`_points.py`**: the same log y-axis, but panel (a) reuses
  `assembler_benchmark_points.py`'s `dodged_points()` -- three individual
  points plus a thin black mean marker, no confidence interval.

Both use the same `N50_LOG_YMIN`/`N50_LOG_YMAX` bounds (0.1-50 Mb, round
decades bracketing the real range) so they are a fair comparison of bars
vs. points specifically, not incidentally different scales.

```bash
python assembly_analysis/scripts/plots/assembler_benchmark_n50_logbars.py
python assembly_analysis/scripts/plots/assembler_benchmark_n50_points.py
```

Outputs `assembler_benchmark_n50_logbars.png`/`.pdf` and
`assembler_benchmark_n50_points.png`/`.pdf` in the same figures directory.
Not adopted as the main figure by default -- generated for side-by-side
comparison against the linear-scale panel (a) in `assembler_benchmark_bars.png`.

### `assembler_benchmark_points.py` (earlier comparison round, kept for reference)

```text
assembly_analysis/scripts/plots/assembler_benchmark_points.py
```

The individual-points alternative from the A/B round above: same 5 metrics
(NGA50 instead of N50/contig count -- it was built before that swap),
identical y-axis scales and x-axis geometry to the bars version at the
time, so the comparison was fair. Not updated to the final 6-metric set --
kept as-is as the record of what was compared, not as a maintained second
version of the adopted figure.

```bash
python assembly_analysis/scripts/plots/assembler_benchmark_points.py
```

Outputs `assembler_benchmark_points.png`/`.pdf` in the same directory. Also
writes `assembler_benchmark_plot_data.tsv` -- the same full schema as the
bars script (the loader is shared and always includes every column, not
just the ones a given script plots), so running either script after the
other does not change that file's contents.

## Scientific comparison rules

1. Compare the same biological sample and input regime whenever possible.
2. Keep ONT, PacBio HiFi and hybrid inputs explicit.
3. Do not treat ntLink as an independent assembler.
4. Record whether an assembly is collapsed or haplotype-resolved.
5. Never replace unavailable metrics with zero.
6. Preserve exact source-file and tool provenance.
