# Long-Read Assembly Benchmark

This directory contains assembly workflows, local / reduced validation infrastructure, and production whole-genome assembly workflows.

## Benchmark scope

Production assembly benchmarking uses GIAB samples `HG002`, `HG003`, and `HG004` with ONT and PacBio HiFi data.

| Method | Role | Input |
|---|---|---|
| Flye | long-read assembler | ONT or PacBio HiFi |
| GoldRush | long-read assembler | ONT or PacBio HiFi |
| Verkko | hybrid assembler | ONT + PacBio HiFi |
| ntLink | post-assembly scaffolder | draft assembly + long reads |

`ntLink` is evaluated as a scaffolding / post-assembly step, not as an independent fourth assembler.

## Two workflow layers

### Production whole-genome workflows

Canonical production workflows are under:

```text
assemblers/whole_genome_asm/
```

See [`whole_genome_asm/README.md`](whole_genome_asm/README.md) for exact execution instructions.

### Modular / reduced validation workflow

The existing `assemblers/Snakefile`, `config/`, `rules/`, `samples.tsv`, and related scripts support prepared-Chr21 and local validation workflows.

They are retained because they serve a different validation purpose from the standalone production whole-genome workflows.

## Assembly quality assessment

QUAST, BUSCO and Merqury assessment workflows belong under:

```text
assemblers/whole_genome_asm/assessment/
```

See [`whole_genome_asm/assessment/README.md`](whole_genome_asm/assessment/README.md).

## Downstream analysis

Metric aggregation, quality-summary tables and assembler figures belong in:

```text
assembly_analysis/
```

This separation keeps **workflow execution** distinct from **result analysis and visualization**.

## Production outputs

The standardized assembly target is:

```text
assemblies/{assembler}/{dataset}/assembly.fasta
```

Verkko may use sample-level output naming because it consumes both technologies for the same sample.

## Important path rule

Run production whole-genome workflows from the repository root. Several standalone workflows include the shared root-level `header_assembler.smk` using their current relative location. Moving these `.smk` files would require a coordinated code change and dry-run validation.

For this reason, the repository reorganization deliberately leaves all production assembler `.smk` files in place.
