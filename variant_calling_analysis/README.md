# Variant-Calling Benchmark Analysis

This directory is reserved for the variant-calling benchmarking component of the F2 project.

The wider repository already contains variant-calling workflows created or modified by multiple contributors. The repository reorganization **does not move, rename or claim ownership of those existing workflows**.

This module is the location for the F2 comparison layer as it is developed:

```text
variant_calling_analysis/
├── README.md
├── scripts/
│   ├── metrics/
│   └── plots/
├── tables/
└── figures/
```

## Repository explorer

<!-- AUTO_REPOSITORY_TREE_START -->
Generated from Git-tracked files. Expand only the directory you need. On GitHub, press **`t`** for fast filename search.

- [`README.md`](README.md)

<details>
<summary><b>figures/</b> — 1 file</summary>

- [`.gitkeep`](figures/.gitkeep)

</details>

<details>
<summary><b>scripts/</b> — 2 files</summary>


<details open>
<summary><b>metrics/</b> — 1 file</summary>

- [`.gitkeep`](scripts/metrics/.gitkeep)

</details>

<details open>
<summary><b>plots/</b> — 1 file</summary>

- [`.gitkeep`](scripts/plots/.gitkeep)

</details>

</details>

<details>
<summary><b>tables/</b> — 1 file</summary>

- [`.gitkeep`](tables/.gitkeep)

</details>

<!-- AUTO_REPOSITORY_TREE_END -->

## Intended responsibilities

- collect caller outputs for a defined benchmark design
- extract comparable metrics
- preserve caller / aligner / sample / technology provenance
- benchmark against appropriate truth sets and confident regions
- create standardized summary tables
- create figures from the canonical tables

## Existing caller workflows

Existing root-level or project-level caller `.smk` files remain in their current locations unless a dedicated, reviewed refactor is performed later.

This prevents the F2 organizational cleanup from breaking workflows that belong to the wider repository.

## Future README expansion

When the final caller set is fixed, document:

- callers included in the comparison
- exact versions and Docker images
- required aligned inputs
- truth sets / confident regions
- evaluation tools and metric definitions
- canonical output table
- figure scripts
- known caller-specific caveats
