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

The workflow uses QUAST 5.3.0 with the pinned Biocontainers image:

```text
quay.io/biocontainers/quast:5.3.0--py313pl5321h5ca1c30_2
```

Production assessment is restricted to `.30x` datasets. Reduced or test datasets containing `.1k`, `.chr21`, `.localtest`, or `smoke` are excluded.

### Reference genome

The workflow uses the GRCh38 reference:

```text
GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta
```

The reference is an external project resource and is not assumed to be stored inside the `lrs_benchmarking` repository.

The project constitution identifies the shared reference location as:

```text
~/smb/Analyses/Reference_sequence/hg38_KGGM/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta
```

For execution, provide the resolved absolute path with:

```bash
--config reference=/absolute/path/to/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta
```

The workflow verifies that the configured FASTA exists, derives its parent directory and mounts that directory read-only into the QUAST container.

When Snakemake is executed through a container, use the resolved absolute server path rather than a literal `~` path.

### Server execution

On the UKW server, Snakemake is provided through the Docker image:

```text
snakemake-orchestrator
```

The current image provides Snakemake 9.19.0, verified with:

```bash
docker run --rm snakemake-orchestrator snakemake --version
```

which returns:

```text
9.19.0
```

The assessment Snakefiles are executed by this Snakemake environment and launch their respective pinned QUAST and BUSCO Docker containers.

The exact production `docker run` command and server-side bind mounts are environment-specific and are therefore not hard-coded here.

The repository root remains the canonical Snakemake working directory. Repository and external-resource paths must be visible consistently to the Snakemake execution environment and the Docker host.

### Dry run

Once the repository and external reference are available in the Snakemake execution environment:

```bash
snakemake --snakefile assemblers/whole_genome_asm/assessment/assembly_quality_quast.smk --cores 16 --resources mem_gb=128 --config reference=/absolute/path/to/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta --dry-run --printshellcmds
```

A missing or inaccessible reference causes the workflow to stop before QUAST execution.

Docker exit code `125` together with missing repository or reference files generally indicates a bind-mount or execution-path problem rather than a QUAST installation problem.

### Production run

After the dry run completes successfully, execute the workflow from the configured Snakemake environment by removing `--dry-run`:

```bash
snakemake --snakefile assemblers/whole_genome_asm/assessment/assembly_quality_quast.smk --cores 16 --resources mem_gb=128 --config reference=/absolute/path/to/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta --printshellcmds
```

On the UKW server, this command is executed within the configured `snakemake-orchestrator` environment. The outer Docker invocation and server-side bind mounts are infrastructure-specific and are therefore not hard-coded here.


## BUSCO

Workflow:

```text
assembly_quality_busco.smk
```

Expected output pattern:

```text
assembly_quality/busco/{assembler}/{dataset}/busco/run_{lineage}/short_summary.json
assembly_quality/busco/{assembler}/{dataset}/busco.log
```

The workflow uses BUSCO 6.1.0 with the pinned Biocontainers image:

```text
quay.io/biocontainers/busco:6.1.0--pyhdfd78af_2
```

Production assessment is restricted to `.30x` datasets. Reduced or test datasets containing `.1k`, `.chr21`, `.localtest`, or `smoke` are excluded.

For human whole-genome assemblies, the F2 workflow is designed around a primate BUSCO lineage.

The lineage directory must be supplied explicitly:

```bash
--config busco_lineage=/absolute/path/to/primates_odb12.2
```

The workflow verifies that the lineage directory exists and contains a valid `dataset.cfg`.

The lineage directory is mounted read-only into the BUSCO container. BUSCO is executed in offline genome mode using Miniprot.

### Dry run

Once the lineage directory is available in the Snakemake execution environment:

```bash
snakemake --snakefile assemblers/whole_genome_asm/assessment/assembly_quality_busco.smk --cores 16 --resources mem_gb=64 --config busco_lineage=/absolute/path/to/primates_odb12.2 --dry-run --printshellcmds
```

The exact BUSCO lineage path is server-specific and is not hard-coded in the workflow or this README.

### Production run

After the dry run completes successfully, execute the workflow from the configured Snakemake environment by removing `--dry-run`:

```bash
snakemake --snakefile assemblers/whole_genome_asm/assessment/assembly_quality_busco.smk --cores 16 --resources mem_gb=64 --config busco_lineage=/absolute/path/to/primates_odb12.2 --printshellcmds
```

On the UKW server, this command is executed within the configured `snakemake-orchestrator` environment. The exact BUSCO lineage path and outer Docker invocation are server-specific and are not hard-coded in the workflow or this README.

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