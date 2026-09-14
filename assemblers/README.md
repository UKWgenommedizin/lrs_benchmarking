# Long-Read Assembly Benchmark

This directory contains assembly workflows, local / reduced validation infrastructure, and production whole-genome assembly workflows.

## Repository explorer

<!-- AUTO_REPOSITORY_TREE_START -->
Generated from Git-tracked files. Expand only the directory you need. On GitHub, press **`t`** for fast filename search.

- [`README.md`](README.md)
- [`samples.tsv`](samples.tsv)

<details open>
<summary><b>benchmark_chr21_real/</b> — 14 files</summary>

- [`BENCHMARK_CONTRACT.md`](benchmark_chr21_real/BENCHMARK_CONTRACT.md)
- [`README.md`](benchmark_chr21_real/README.md)
- [`Snakefile.assemblers`](benchmark_chr21_real/Snakefile.assemblers)
- [`Snakefile.inputs`](benchmark_chr21_real/Snakefile.inputs)

<details>
<summary><b>config/</b> — 5 files</summary>

- [`config.assemblers.yaml`](benchmark_chr21_real/config/config.assemblers.yaml)
- [`containers.yaml`](benchmark_chr21_real/config/containers.yaml)
- [`dockerhub_images.lock.tsv`](benchmark_chr21_real/config/dockerhub_images.lock.tsv)
- [`samples.assemblers.tsv`](benchmark_chr21_real/config/samples.assemblers.tsv)
- [`sources.resolved.tsv`](benchmark_chr21_real/config/sources.resolved.tsv)

</details>

<details>
<summary><b>scripts/</b> — 4 files</summary>

- [`filter_sam_start_window.py`](benchmark_chr21_real/scripts/filter_sam_start_window.py)
- [`normalize_chr21_bam_to_coverage.py`](benchmark_chr21_real/scripts/normalize_chr21_bam_to_coverage.py)
- [`resolve_chr21_sources.py`](benchmark_chr21_real/scripts/resolve_chr21_sources.py)
- [`validate_extracted_chr21_bam.py`](benchmark_chr21_real/scripts/validate_extracted_chr21_bam.py)

</details>

<details open>
<summary><b>user_settings/</b> — 1 file</summary>

- [`input_data_path.example.yaml`](benchmark_chr21_real/user_settings/input_data_path.example.yaml)

</details>

</details>

<details>
<summary><b>containers/</b> — 7 files</summary>


<details open>
<summary><b>flye2/</b> — 1 file</summary>

- [`Dockerfile`](containers/flye2/Dockerfile)

</details>

<details open>
<summary><b>goldrush/</b> — 2 files</summary>

- [`Dockerfile`](containers/goldrush/Dockerfile)
- [`build_goldrush_container.py`](containers/goldrush/build_goldrush_container.py)

</details>

<details open>
<summary><b>ntlink/</b> — 2 files</summary>

- [`Dockerfile`](containers/ntlink/Dockerfile)
- [`build_ntlink_container.py`](containers/ntlink/build_ntlink_container.py)

</details>

<details open>
<summary><b>verkko2/</b> — 2 files</summary>

- [`Dockerfile`](containers/verkko2/Dockerfile)
- [`build_verkko_container.py`](containers/verkko2/build_verkko_container.py)

</details>

</details>

<details open>
<summary><b>envs/</b> — 1 file</summary>

- [`flye.yaml`](envs/flye.yaml)

</details>

<details>
<summary><b>scripts/</b> — 6 files</summary>

- [`check_assembler_container_definitions.sh`](scripts/check_assembler_container_definitions.sh)
- [`workflow_status.py`](scripts/workflow_status.py)

<details open>
<summary><b>chr21/</b> — 4 files</summary>

- [`filter_sam_start_window.py`](scripts/chr21/filter_sam_start_window.py)
- [`normalize_chr21_bam_to_coverage.py`](scripts/chr21/normalize_chr21_bam_to_coverage.py)
- [`select_one_primary_per_qname.py`](scripts/chr21/select_one_primary_per_qname.py)
- [`validate_extracted_chr21_bam.py`](scripts/chr21/validate_extracted_chr21_bam.py)

</details>

</details>

<details open>
<summary><b>whole_genome_asm/</b> — 11 files</summary>

- [`README.md`](whole_genome_asm/README.md)
- [`hybrid.assembly.verkko.smk`](whole_genome_asm/hybrid.assembly.verkko.smk)
- [`ont.assembly.flye2.smk`](whole_genome_asm/ont.assembly.flye2.smk)
- [`ont.assembly.goldrush.smk`](whole_genome_asm/ont.assembly.goldrush.smk)
- [`ont.scaffolding.ntlink.smk`](whole_genome_asm/ont.scaffolding.ntlink.smk)
- [`pb.assembly.flye2.smk`](whole_genome_asm/pb.assembly.flye2.smk)
- [`pb.assembly.goldrush.smk`](whole_genome_asm/pb.assembly.goldrush.smk)
- [`pb.scaffolding.ntlink.smk`](whole_genome_asm/pb.scaffolding.ntlink.smk)

<details>
<summary><b>assessment/</b> — 3 files</summary>

- [`README.md`](whole_genome_asm/assessment/README.md)
- [`assembly_quality_busco.smk`](whole_genome_asm/assessment/assembly_quality_busco.smk)
- [`assembly_quality_quast.smk`](whole_genome_asm/assessment/assembly_quality_quast.smk)

</details>

</details>

<!-- AUTO_REPOSITORY_TREE_END -->

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
