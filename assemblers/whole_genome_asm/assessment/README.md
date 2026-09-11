# Assembly Quality Assessment

This directory contains workflow-level quality assessment for completed whole-genome assemblies.

The assessment is intentionally separated from `assembly_analysis/`:

- `assemblers/whole_genome_asm/assessment/` **runs assessment tools**
- `assembly_analysis/` **aggregates metrics, performs downstream analysis and creates figures**

## Repository explorer

<!-- AUTO_REPOSITORY_TREE_START -->
Generated from Git-tracked files. Expand only the directory you need. On GitHub, press **`t`** for fast filename search.

- [`README.md`](README.md)
- [`assembly_quality_busco.smk`](assembly_quality_busco.smk)
- [`assembly_quality_quast.smk`](assembly_quality_quast.smk)

<!-- AUTO_REPOSITORY_TREE_END -->

## Input contract

Assessment workflows discover completed assemblies using the standardized layout:

```text
assemblies/{assembler}/{dataset}/assembly.fasta
```

The assessed methods can include:

- Flye
- GoldRush
- Verkko
- ntLink outputs, interpreted as scaffolder results rather than an independent assembler

## Complementary assessment tools

| Tool | Main scientific role |
|---|---|
| QUAST / QUAST-LG | contiguity and reference-based structural agreement |
| BUSCO | conserved-gene completeness |
| Merqury | k-mer completeness and consensus quality (QV) |

No single metric should be used as the sole assembly-quality criterion.

## QUAST

Workflow:

```text
assembly_quality_quast.smk
```

Expected output pattern:

```text
assembly_quality/quast/{assembler}/{dataset}/report.tsv
assembly_quality/quast/{assembler}/{dataset}/quast.log
```

The current workflow uses QUAST 5.3.0 with a pinned Biocontainers image.

The reference path can be supplied through Snakemake config. Relative paths should be resolved against the repository root; absolute server paths are also supported by the current path-resolution pattern.

Dry-run example:

```bash
snakemake --snakefile assemblers/whole_genome_asm/assessment/assembly_quality_quast.smk --cores 16 --resources mem_gb=128 --config reference=reference/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta --dry-run --printshellcmds
```

## BUSCO

Workflow when present:

```text
assembly_quality_busco.smk
```

For human whole-genome assemblies, the F2 workflow is designed around a primate BUSCO lineage. The lineage directory must be supplied explicitly and should contain a valid `dataset.cfg`.

Example dry-run pattern:

```bash
snakemake --snakefile assemblers/whole_genome_asm/assessment/assembly_quality_busco.smk --cores 16 --resources mem_gb=64 --config busco_lineage=/path/to/primates_odb12.2 --dry-run --printshellcmds
```

Do not copy the example lineage path literally; use the real extracted lineage directory on the execution server.

## Merqury

Workflow when present:

```text
assembly_quality_merqury.smk
```

Merqury should use the intended trusted k-mer source for the scientific comparison. The exact input data, k-mer database construction, container image, and output contract should be documented in the workflow and here once finalized.

Do not invent missing Merqury provenance or substitute a different k-mer source silently.

## Interpretation

### Contiguity

N50 is useful descriptively but is not sufficient to establish assembly correctness. Reference-aware metrics such as NG50 / NGA50, structural discrepancies, genome fraction and duplication should be considered where appropriate.

### Completeness

BUSCO completeness and k-mer completeness measure different properties and should be reported separately.

### Verkko representation

A haplotype-resolved / diploid Verkko output can have a different total representation from a collapsed assembly. Duplicated BUSCOs, total assembly size and related metrics must therefore be interpreted in the context of assembly representation rather than treated automatically as errors.

## Downstream aggregation

Final cross-tool tables and plots belong in:

```text
assembly_analysis/
```

This keeps raw tool outputs reproducible while allowing downstream analysis scripts to evolve independently.
