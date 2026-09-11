# lrs_benchmarking

Reproducible benchmarking workflows and analyses for long-read sequencing data.

This repository contains shared project workflows as well as the F2 benchmarking work focused on three scientific stages:

1. **Long-read alignment benchmarking**
2. **Whole-genome assembly benchmarking and assembly-quality assessment**
3. **Variant-calling benchmarking** (current / next project phase)

The repository is governed by [`CONSTITUTION.md`](CONSTITUTION.md). In particular, Snakemake workflows are executed from the **repository root**, use pinned Docker images, and communicate through files on disk.

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
