# lrs_benchmarking

Reproducible benchmarking workflows and analyses for long-read sequencing data.

This repository contains shared project workflows as well as the F2 benchmarking work focused on three scientific stages:

1. **Long-read alignment benchmarking**
2. **Whole-genome assembly benchmarking and assembly-quality assessment**
3. **Variant-calling benchmarking** (current / next project phase)

The repository is governed by [`CONSTITUTION.md`](CONSTITUTION.md). In particular, Snakemake workflows are executed from the **repository root**, use pinned Docker images, and communicate through files on disk.

## Repository explorer

<!-- AUTO_REPOSITORY_TREE_START -->
Generated from Git-tracked files. Expand only the area you need. On GitHub, press **`t`** to search tracked files by name.

<details>
<summary><b>Top-level files</b></summary>

- [`.gitignore`](.gitignore)
- [`CONSTITUTION.md`](CONSTITUTION.md)
- [`README.md`](README.md)
- [`Worflow and code implemented`](Worflow and code implemented)
- [`caller_tools_inventory.csv`](caller_tools_inventory.csv)
- [`cuteSV.hg38.smk`](cuteSV.hg38.smk)
- [`header.smk`](header.smk)
- [`header_mapper.smk`](header_mapper.smk)
- [`ilmn.snv_calling.clair3.smk`](ilmn.snv_calling.clair3.smk)
- [`lrs_benchmarking_tools_selection.R.Rmd`](lrs_benchmarking_tools_selection.R.Rmd)
- [`lrs_benchmarking_tools_selection.xlsx`](lrs_benchmarking_tools_selection.xlsx)
- [`ont.read_mapping.minimap2.smk`](ont.read_mapping.minimap2.smk)
- [`ont.read_mapping.pbmm2.smk`](ont.read_mapping.pbmm2.smk)
- [`ont.read_mapping.vacmap.smk`](ont.read_mapping.vacmap.smk)
- [`ont.read_mapping.vg.smk`](ont.read_mapping.vg.smk)
- [`ont.snv_calling.clair3.smk`](ont.snv_calling.clair3.smk)
- [`ont.snv_calling.deepvariant.smk`](ont.snv_calling.deepvariant.smk)
- [`ont.snv_filtering.v3.smk`](ont.snv_filtering.v3.smk)
- [`pb.read_mapping.minimap2.smk`](pb.read_mapping.minimap2.smk)
- [`pb.read_mapping.pbmm2.smk`](pb.read_mapping.pbmm2.smk)
- [`pb.read_mapping.vacmap.smk`](pb.read_mapping.vacmap.smk)
- [`pb.read_mapping.vg.smk`](pb.read_mapping.vg.smk)
- [`pb.snv_calling.clair3.smk`](pb.snv_calling.clair3.smk)
- [`pb.snv_calling.deepvariant.smk`](pb.snv_calling.deepvariant.smk)
- [`pbsv.hg38.smk`](pbsv.hg38.smk)
- [`run_happy.smk`](run_happy.smk)
- [`sawfish.hg38.smk`](sawfish.hg38.smk)
- [`sniffles2.hg38.smk`](sniffles2.hg38.smk)
- [`sv_f1_barplots.png`](sv_f1_barplots.png)
- [`truvari_anno.smk`](truvari_anno.smk)
- [`wgs.snv_calling.deepvariant.smk`](wgs.snv_calling.deepvariant.smk)

</details>

<details>
<summary><b>.githooks/</b> — 1 tracked file</summary>

- [`pre-commit`](.githooks/pre-commit)

</details>

<details>
<summary><b>.github/</b> — 1 tracked file</summary>

- `workflows/` — 1 file

</details>

<details>
<summary><b>SV aligners call/</b> — 11 tracked files</summary>

- `configs/` — 1 file
- `docs/` — 4 files
- `scripts/` — 1 file
- `workflow/` — 4 files
- [`CONSTITUTION_NICOLAS.md`](SV aligners call/CONSTITUTION_NICOLAS.md)

</details>

<details>
<summary><b>Workflow_and_code_implemented/</b> — 6 tracked files</summary>

- Documentation: [`README.md`](Workflow_and_code_implemented/README.md)
- [`ANALYSIS_INTERPRETATION.md`](Workflow_and_code_implemented/ANALYSIS_INTERPRETATION.md)
- [`FILE_INPUT_OUTPUT_MAP.md`](Workflow_and_code_implemented/FILE_INPUT_OUTPUT_MAP.md)
- [`PROJECT_WORKFLOW.md`](Workflow_and_code_implemented/PROJECT_WORKFLOW.md)
- [`REPRODUCTION_CHECKLIST.md`](Workflow_and_code_implemented/REPRODUCTION_CHECKLIST.md)
- [`REUSABLE_CODE_PATTERNS.md`](Workflow_and_code_implemented/REUSABLE_CODE_PATTERNS.md)

</details>

<details>
<summary><b>alignment_analysis/</b> — 154 tracked files</summary>

- Documentation: [`README.md`](alignment_analysis/README.md)
- `figures/` — 44 files
- `scripts/` — 36 files
- `tables/` — 71 files
- [`COMMAND_LOG.md`](alignment_analysis/COMMAND_LOG.md)
- [`README_assemblers.md`](alignment_analysis/README_assemblers.md)

</details>

<details>
<summary><b>archive/</b> — 17 tracked files</summary>

- `legacy/` — 13 files
- [`ont.read_mapping.ntlink.smk`](archive/ont.read_mapping.ntlink.smk)
- [`ont.read_mapping.quicked.smk`](archive/ont.read_mapping.quicked.smk)
- [`pb.read_mapping.ntlink.smk`](archive/pb.read_mapping.ntlink.smk)
- [`pb.read_mapping.quicked.smk`](archive/pb.read_mapping.quicked.smk)

</details>

<details>
<summary><b>assemblers/</b> — 41 tracked files</summary>

- Documentation: [`README.md`](assemblers/README.md)
- `benchmark_chr21_real/` — 14 files
- `containers/` — 7 files
- `envs/` — 1 file
- `scripts/` — 6 files
- `whole_genome_asm/` — 11 files
- [`samples.tsv`](assemblers/samples.tsv)

</details>

<details>
<summary><b>assembly_analysis/</b> — 6 tracked files</summary>

- Documentation: [`README.md`](assembly_analysis/README.md)
- `figures/` — 1 file
- `scripts/` — 3 files
- `tables/` — 1 file

</details>

<details>
<summary><b>docs/</b> — 5 tracked files</summary>

- `f2/` — 3 files
- [`mapper_cmds.sh`](docs/mapper_cmds.sh)
- [`retag_docker_images.sh`](docs/retag_docker_images.sh)

</details>

<details>
<summary><b>figures/</b> — 3 tracked files</summary>

- [`repo1.2_SV_FP.Rmd`](figures/repo1.2_SV_FP.Rmd)
- [`repo1_SV_type_differences.Rmd`](figures/repo1_SV_type_differences.Rmd)
- [`repo1_SV_type_differences.html`](figures/repo1_SV_type_differences.html)

</details>

<details>
<summary><b>final_report_files/</b> — 24 tracked files</summary>

- `snakemake_aligners_benchmarking/` — 24 files

</details>

<details>
<summary><b>happy_results/</b> — 146 tracked files</summary>

- `146 direct files` — use GitHub's **`t`** file finder or the module README to locate a specific file

</details>

<details>
<summary><b>mapper_legacy/</b> — 4 tracked files</summary>

- [`ont.read_mapping.graphaligner.smk`](mapper_legacy/ont.read_mapping.graphaligner.smk)
- [`ont.read_mapping.parahat.smk`](mapper_legacy/ont.read_mapping.parahat.smk)
- [`pb.read_mapping.graphaligner.smk`](mapper_legacy/pb.read_mapping.graphaligner.smk)
- [`pb.read_mapping.parahat.smk`](mapper_legacy/pb.read_mapping.parahat.smk)

</details>

<details>
<summary><b>run_metrics/</b> — 2 tracked files</summary>

- [`create_run_metrics.sh`](run_metrics/create_run_metrics.sh)
- [`mm2.run_metrics.tsv`](run_metrics/mm2.run_metrics.tsv)

</details>

<details>
<summary><b>samtools_stats_30x_Christian/</b> — 25 tracked files</summary>

- `25 direct files` — use GitHub's **`t`** file finder or the module README to locate a specific file

</details>

<details>
<summary><b>sawfish/</b> — 5 tracked files</summary>

- `HG002.ont.30x.hg38.pbmm2-pb.discover_dir/` — 1 file
- `HG002.pb.1k.hg38.pbmm2-ont.discover_dir/` — 1 file
- `HG002.pb.1k.hg38.pbmm2-pb.discover_dir/` — 1 file
- `HG002.pb.30x.hg38.pbmm2-ont.discover_dir/` — 1 file
- `HG002.pb.30x.hg38.pbmm2-pb.discover_dir/` — 1 file

</details>

<details>
<summary><b>scripts/</b> — 2 tracked files</summary>

- [`install_git_hooks.sh`](scripts/install_git_hooks.sh)
- [`update_repository_tree.py`](scripts/update_repository_tree.py)

</details>

<details>
<summary><b>truvari/</b> — 1091 tracked files</summary>

- `HG002.ont.1k.hg38.mm2-ont.cuteSV/` — 13 files
- `HG002.ont.1k.hg38.mm2-ont.sniffles2/` — 13 files
- `HG002.ont.1k.hg38.mm2-pb.cuteSV/` — 13 files
- `HG002.ont.1k.hg38.mm2-pb.sniffles2/` — 13 files
- `HG002.ont.1k.hg38.pbmm2-ont.cuteSV/` — 13 files
- `HG002.ont.1k.hg38.pbmm2-ont.pbsv/` — 13 files
- `HG002.ont.1k.hg38.pbmm2-ont.sniffles2/` — 13 files
- `HG002.ont.1k.hg38.pbmm2-pb.cuteSV/` — 13 files
- `HG002.ont.1k.hg38.pbmm2-pb.pbsv/` — 13 files
- `HG002.ont.1k.hg38.pbmm2-pb.sniffles2/` — 13 files
- `HG002.ont.30x.hg38.mm2-ont.cuteSV/` — 33 files
- `HG002.ont.30x.hg38.mm2-ont.sniffles2/` — 33 files
- `HG002.ont.30x.hg38.mm2-pb.cuteSV/` — 33 files
- `HG002.ont.30x.hg38.mm2-pb.sniffles2/` — 33 files
- `HG002.ont.30x.hg38.pbmm2-ont.cuteSV/` — 33 files
- `HG002.ont.30x.hg38.pbmm2-ont.pbsv/` — 33 files
- `HG002.ont.30x.hg38.pbmm2-ont.sniffles2/` — 33 files
- `HG002.ont.30x.hg38.pbmm2-pb.cuteSV/` — 33 files
- `HG002.ont.30x.hg38.pbmm2-pb.pbsv/` — 33 files
- `HG002.ont.30x.hg38.pbmm2-pb.sawfish/` — 33 files
- `HG002.ont.30x.hg38.pbmm2-pb.sniffles2/` — 33 files
- `HG002.pb.1k.hg38.mm2-ont.cuteSV/` — 13 files
- `HG002.pb.1k.hg38.mm2-ont.sniffles2/` — 13 files
- `HG002.pb.1k.hg38.mm2-pb.cuteSV/` — 13 files
- `HG002.pb.1k.hg38.mm2-pb.sniffles2/` — 13 files
- `HG002.pb.1k.hg38.pbmm2-ont.cuteSV/` — 13 files
- `HG002.pb.1k.hg38.pbmm2-ont.pbsv/` — 13 files
- `HG002.pb.1k.hg38.pbmm2-ont.sawfish/` — 13 files
- `HG002.pb.1k.hg38.pbmm2-ont.sniffles2/` — 13 files
- `HG002.pb.1k.hg38.pbmm2-pb.cuteSV/` — 13 files
- `HG002.pb.1k.hg38.pbmm2-pb.pbsv/` — 13 files
- `HG002.pb.1k.hg38.pbmm2-pb.sawfish/` — 13 files
- `HG002.pb.1k.hg38.pbmm2-pb.sniffles2/` — 13 files
- `HG002.pb.30x.hg38.mm2-ont.cuteSV/` — 33 files
- `HG002.pb.30x.hg38.mm2-ont.sniffles2/` — 33 files
- `HG002.pb.30x.hg38.mm2-pb.cuteSV/` — 33 files
- `HG002.pb.30x.hg38.mm2-pb.sniffles2/` — 33 files
- `HG002.pb.30x.hg38.pbmm2-ont.cuteSV/` — 33 files
- `HG002.pb.30x.hg38.pbmm2-ont.pbsv/` — 33 files
- `HG002.pb.30x.hg38.pbmm2-ont.sawfish/` — 33 files
- `HG002.pb.30x.hg38.pbmm2-ont.sniffles2/` — 33 files
- `HG002.pb.30x.hg38.pbmm2-pb.cuteSV/` — 33 files
- `HG002.pb.30x.hg38.pbmm2-pb.pbsv/` — 33 files
- `HG002.pb.30x.hg38.pbmm2-pb.sawfish/` — 33 files
- `HG002.pb.30x.hg38.pbmm2-pb.sniffles2/` — 33 files
- `46 direct files` — use GitHub's **`t`** file finder or the module README to locate a specific file

</details>

<details>
<summary><b>variant_calling_analysis/</b> — 5 tracked files</summary>

- Documentation: [`README.md`](variant_calling_analysis/README.md)
- `figures/` — 1 file
- `scripts/` — 2 files
- `tables/` — 1 file

</details>

<details>
<summary><b>vcf_called/</b> — 25 tracked files</summary>

- `snv_indel/` — 24 files
- [`.DS_Store`](vcf_called/.DS_Store)

</details>

<!-- AUTO_REPOSITORY_TREE_END -->

## Benchmark overview

```text
FASTQ
  |
  v
LONG-READ ALIGNMENT
  minimap2 | pbmm2 | VACMap | VG Giraffe
  |
  +--> CRAM + alignment QC / metrics
  |        |
  |        v
  |     alignment_analysis/
  |
  +-------------------------------+
  |                               |
  v                               v
WHOLE-GENOME ASSEMBLY          VARIANT CALLING
  Flye                         benchmarking module
  GoldRush                     (developed separately)
  Verkko
  |
  +--> ntLink scaffolding
  |
  v
ASSEMBLY QUALITY ASSESSMENT
  QUAST | BUSCO | Merqury
  |
  v
assembly_analysis/
```

## F2 benchmark design

The production benchmark uses GIAB samples **HG002, HG003 and HG004** with:

- Oxford Nanopore Technologies (ONT)
- PacBio HiFi
- production 30x datasets where available

### Alignment benchmark

| Aligner | ONT tag | PacBio HiFi tag |
|---|---|---|
| minimap2 | `mm2-ont` | `mm2-pb` |
| pbmm2 | `pbmm2-ont` | `pbmm2-pb` |
| VACMap | `vacmap-ont` | `vacmap-pb` |
| VG Giraffe | `vg-ont` | `vg-pb` |

Production alignment analysis is documented in [`alignment_analysis/README.md`](alignment_analysis/README.md).

### Assembly benchmark

| Tool | Role | Input |
|---|---|---|
| Flye | long-read assembler | ONT or PacBio HiFi |
| GoldRush | long-read assembler | ONT or PacBio HiFi |
| Verkko | hybrid assembler | ONT + PacBio HiFi |
| ntLink | post-assembly scaffolder | draft assembly + long reads |

Whole-genome assembler execution is documented in [`assemblers/whole_genome_asm/README.md`](assemblers/whole_genome_asm/README.md).

Assembly-quality assessment is documented in [`assemblers/whole_genome_asm/assessment/README.md`](assemblers/whole_genome_asm/assessment/README.md).

Downstream metric aggregation and figure generation belong in [`assembly_analysis/`](assembly_analysis/README.md).

### Variant-calling benchmark

The repository already contains variant-calling workflows contributed across the wider project. Those existing workflows are **not moved by the F2 repository reorganization**.

The F2 analysis layer for future/current variant-calling comparisons is documented in [`variant_calling_analysis/README.md`](variant_calling_analysis/README.md).

## Repository layout

```text
lrs_benchmarking/
├── README.md
├── CONSTITUTION.md
├── header_mapper.smk
├── header_assembler.smk
│
├── ont.read_mapping.*.smk          # active mapper workflows: keep at root
├── pb.read_mapping.*.smk           # active mapper workflows: keep at root
│
├── alignment_analysis/
│   ├── README.md
│   ├── scripts/
│   │   └── 30x/quality_check_aligners.py
│   ├── tables/
│   └── figures/
│
├── assemblers/
│   ├── README.md
│   └── whole_genome_asm/
│       ├── README.md
│       ├── *.assembly.*.smk
│       └── assessment/
│           ├── README.md
│           ├── assembly_quality_quast.smk
│           ├── assembly_quality_busco.smk       # when present
│           └── assembly_quality_merqury.smk     # when present
│
├── assembly_analysis/
│   ├── README.md
│   ├── scripts/
│   │   ├── metrics/
│   │   └── plots/
│   ├── tables/
│   └── figures/
│
├── variant_calling_analysis/
│   ├── README.md
│   ├── scripts/
│   │   ├── metrics/
│   │   └── plots/
│   ├── tables/
│   └── figures/
│
├── docs/
│   └── f2/
│       ├── REPOSITORY_LAYOUT.md
│       └── PATH_STABILITY.md
│
└── existing shared project files and workflows
```

## Execution rule

Run Snakemake from the repository root unless a workflow explicitly documents otherwise:

```bash
cd /path/to/lrs_benchmarking
snakemake --snakefile path/to/workflow.smk --dry-run --printshellcmds
```

Do not relocate active mapper workflows or the whole-genome assembler workflows merely for aesthetics. Some workflows use shared headers and repository-root-relative paths, so moving them without updating and validating every dependency can break execution.

## Analysis-output convention

Each benchmark domain keeps its analysis products together:

```text
alignment_analysis/
    scripts/
    tables/
    figures/

assembly_analysis/
    scripts/
    tables/
    figures/

variant_calling_analysis/
    scripts/
    tables/
    figures/
```

Workflow-generated biological data remain in their established production locations such as `cram/`, `assemblies/`, and `assembly_quality/`.

## Reproducibility

Before committing workflow changes:

```bash
git diff --check
snakemake --snakefile path/to/workflow.smk --dry-run --printshellcmds
```

Before moving any existing file:

```bash
git grep -n "old/path/or/filename"
git mv old/path new/path
```

Then update every reference, validate again, and commit the move separately from unrelated scientific changes.
