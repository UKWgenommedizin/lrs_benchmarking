# Constitution of the lrs_benchmarking Project

## Preamble

This constitution establishes the governing principles, architectural rules, naming standards, quality requirements, and decision-making priorities for the lrs_benchmarking repository. All workflows, specifications, and code contributions must conform to these rules. No implementation may violate the constitution, and any proposed change to the constitution itself requires explicit review and consensus.

---

## Article I — Repository Structure

**I.1.** The repository root is the single working directory (`CWD`) for all Snakemake executions. No workflow may reference paths outside `CWD` except:
- Reference databases and truth sets (via configurable paths)
- Docker registries (via configurable image tags)

**I.2.** Mapping workflow files follow the pattern `<platform>.read_mapping.<tool>.smk`, where:
- `<platform>` is one of: `ont`, `pb`
- `<tool>` is the mapper name (e.g., `minimap2`, `pbmm2`, `graphaligner`)

**I.3.** SNV/indel calling workflow files follow the pattern `<platform>.snv_calling.<tool>.smk`, where:
- `<platform>` is one of: `ont`, `pb`, `ilmn`, `wgs`
- `<tool>` is the caller name (e.g., `clair3`, `deepvariant`)

**I.4.** SV calling workflow files follow the pattern `<tool>.<reference>.smk` (e.g., `sniffles2.hg38.smk`, `cuteSV.hg38.smk`).

**I.5.** Benchmarking workflow files follow the pattern `<process>.smk` (e.g., `run_happy.smk`).

**I.6.** Non-mapping workflow files must include `header.smk` as their first executable line. Mapping workflow files (`<platform>.read_mapping.<tool>.smk`) must include `header_mapper.smk` as their first executable line. `header_mapper.smk` provides `CWD` only, with no VarCAD dependency.

**I.7.** The constitution resides in `CONSTITUTION.md` at the repository root. It is immutable except by amendment per Article X.

**I.8.** Documentation lives in `docs/`. Analysis notebooks and figures live in `figures/`.

---

## Article II — Execution Environment

**II.1.** All tools must be invoked via Docker containers. No host-installed binaries or conda environments may be called directly. This includes mappers, samtools, and all utility tools.

**II.2.** Every Docker image must be pinned to a specific tag (never `latest`).

**II.3.** Execution uses `docker run` directly — not `srun`, not Slurm, not any HPC job scheduler.

**II.4.** Resource constraints (CPU, memory) are specified via Docker flags `--cpus` and `-m`.

**II.5.** Each rule must set `--rm`, `--workdir /tmp`, and bind-mount `{CWD}:{CWD}` at minimum. Reference databases are mounted read-only.

**II.6.** The user ID and group ID must be propagated via `-u $UID:$(id -g)` to ensure file ownership matches the host user.

---

## Article III — Input, Output, and Reference Data Locations

**III.1.** All data resides on the group SMB share, accessible only from the remote workstation **`genmedbfx`** via the `~/smb/` mount point. The paths below are valid on that machine only. The three distinct path categories are:

### Reference database — `~/smb/Analyses/Reference_sequence/hg38_KGGM/`
Contains the reference genome:

| Resource | Path |
|---|---|
| Reference FASTA (hg38) | `~/smb/Analyses/Reference_sequence/hg38_KGGM/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta` |

Other reference-associated files (alternative FASTA, SDF template, stratification BEDs) were used in the original `run_happy.smk` workflow at paths under `/mnt/storage/db/` — these will need to be located or re-configured for the new benchmarking runs.

### Raw sequence data — `~/smb/Bioinformatics/lrs_benchmarking/fastq/`
Contains FASTQ files for the benchmarking samples. The expected input pattern in workflows is:
```
fastq/<dataset>.fastq.gz
```
where `<dataset>` follows the pattern `<sample>.<platform>.<coverage>`, e.g. `HG002.ont.30x`.

Previously used samples: HG002, HG003, HG004 (the GIAB Ashkenazi trio). The "1k" variant (e.g., `HG002.ont.1k`) refers to downsampled FASTQs for testing.

### Tool-specific data within CWD
Certain truth sets and reference annotations live directly in the working directory:

| File | Purpose |
|---|---|
| `happy_data/NA24385_GRCh38_1_22_v4.2.1_benchmark.vcf.gz` | GIAB small variant truth set (HG002) |
| `happy_data/NA24385_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed` | High-confidence regions for small variants |
| `GRCh38_HG2-T2TQ100-V1.1_stvar.svtype.vcf.gz` | T2T-Q100 SV truth set |
| `GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.only_autosomes.bed` | SV benchmark regions |
| `sniffles2/human_GRCh38_no_alt_analysis_set.trf.bed` | Tandem repeat annotations for Sniffles2 |

**III.2.** Output directory structure is fixed:

| Content | Directory |
|---|---|
| Aligned CRAMs + QC | `cram/` |
| Temporary mapping files | `cram/tmp/` |
| SNV/indel VCFs | `vcf_called/snv_indel/` |
| SNV/indel intermediates | `vcf_called/snv_indel/tmp/` |
| SV VCFs (Sniffles2) | `sniffles2/` |
| SV VCFs (cuteSV) | `cuteSV/` |
| SV VCFs (pbsv) | `pbsv/` |
| SV VCFs (Sawfish) | `sawfish/` |
| Truvari results | `truvari/` |
| hap.py results | `happy_results/` |
| Run metrics | `run_metrics/` |

---

## Article IV — File Naming Convention

**IV.1.** All intermediate and output files must follow the pattern:
```
<dataset>.<reference>.<mapper_tag>.<optional_suffix>
```

**IV.2.** The `<dataset>` field encodes sample, platform, and coverage as dot-separated tokens:
```
<sample>.<platform>.<coverage>
```
Example: `HG002.ont.30x`

**IV.3.** The `<reference>` field is always `hg38`.

**IV.4.** The `<mapper_tag>` identifies the mapper and target platform as `<tool>-<platform>`.

Registered mapper tags (immutable once assigned):

| Mapper | ONT tag | PacBio tag |
|---|---|---|
| minimap2 | `mm2-ont` | `mm2-pb` |
| pbmm2 | `pbmm2-ont` | `pbmm2-pb` |
| GraphAligner | `ga-ont` | `ga-pb` |
| ParaHAT | `parahat-ont` | `parahat-pb` |
| VACmap | `vacmap-ont` | `vacmap-pb` |
| VG | `vg-ont` | `vg-pb` |

**IV.5.** SNV/indel caller tags follow the pattern `<tool>-<platform>` (e.g., `clair3-ont`, `dv-ont-woPG`). The caller tag is appended to the mapper-tagged filename:
```
<dataset>.<reference>.<mapper_tag>.<caller_tag>.vcf.gz
```

**IV.6.** SV caller outputs follow the pattern:
```
<sv_dir>/<dataset>.<reference>.<mapper_tag>.<sv_tool>.vcf.gz
```

**IV.7.** Benchmarking outputs follow the pattern:
```
happy_results/<dataset>.<reference>.<mapper_tag>.<caller_tag>.{summary,extended,roc}.csv
```

---

## Article V — Data Flow and Workflow Independence

**V.1.** Workflows are independent, composable units. No workflow may `include` another workflow file (except `header.smk`).

**V.2.** Workflows communicate only through files on disk using the naming conventions in Article IV. No workflow may pass data to another via environment variables, shared memory, or inter-process communication.

**V.3.** A workflow discovers its inputs by globbing the filesystem. It must not require a manifest or input list to be provided externally.

**V.4.** The execution order is:
1. Read mapping → produces CRAMs in `cram/`
2. SNV/indel calling → consumes CRAMs, produces VCFs in `vcf_called/snv_indel/`
3. SV calling → consumes CRAMs, produces VCFs + Truvari results
4. Benchmarking (hap.py) → consumes VCFs, produces `happy_results/`

Each stage is optional. Any workflow may be run independently if its inputs already exist.

**V.5.** No workflow may modify or delete files produced by another workflow.

**V.6.** Because a new mapper produces CRAMs with a novel `<mapper_tag>`, it is automatically discovered by all downstream SNV callers, SV callers, and benchmarking workflows — no downstream code changes are required.

---

## Article VI — Quality and Validation

**VI.1.** Every output-producing rule must include a post-execution validation check that causes the rule to fail (non-zero exit code) if the output is empty, truncated, or otherwise invalid.

**VI.2.** Mandatory validation checks:
- **CRAM**: `[[ $(du -b {output.cram} | cut -f 1) -le 64 ]] && exit 101`
- **VCF**: `[[ $(bcftools view -H {output.vcf_gz} | wc -l) -le 1 ]] && exit 101`
- **Stats**: `[[ $(du -b {output.stats} | cut -f 1) -lt 5000 ]] && exit 101`

**VI.3.** Every rule must produce a log file capturing stdout and stderr.

**VI.4.** Every rule must log start time, end time, and container hostname.

**VI.5.** Snakemake workflow files must be syntactically valid and pass `snakemake -n` (dry-run) before commit.

---

## Article VII — Containerization Requirements

**VII.1.** Each mapper's Docker image must contain:
- The mapper binary
- Samtools (for sorting, indexing, and QC)
- bcftools (for VCF validation)
- All shared libraries required by the above

**VII.2.** If samtools/bcftools cannot be bundled in the mapper image, a separate utility image must be provided and used in a multi-step or piped workflow.

**VII.3.** Docker images are referenced by full registry path and pinned tag. The default internal registry is `storage-node:5000`; public images are pulled from Docker Hub.

**VII.4.** The Docker `ENTRYPOINT` may be overridden via `--entrypoint` in the `docker run` command when the container provides multiple tools.

---

## Article VIII — Adding a New Mapper

**VIII.1.** A new mapper requires:
1. A new `.smk` file following the naming convention in Article I.2
2. A unique `<mapper_tag>` registered in Article IV.4
3. A Docker image meeting Article VII requirements
4. Successful test run on the `1k` downsample dataset before production runs

**VIII.2.** The new mapper workflow must produce exactly the same set of output files as the existing mapping workflows:
- `cram/<dataset>.<reference>.<mapper_tag>.cram`
- `cram/<dataset>.<reference>.<mapper_tag>.cram.crai`
- `cram/<dataset>.<reference>.<mapper_tag>.cram.idxstats`
- `cram/<dataset>.<reference>.<mapper_tag>.cram.stats`

**VIII.3.** The new mapper workflow must not require changes to any downstream workflow. If downstream changes are needed, the mapper is not yet compliant.

**VIII.4.** The new mapper must accept FASTQ input (`fastq/<dataset>.fastq.gz`) and produce sorted CRAM output with the read group header `@RG\tID:{dataset}\tSM:{dataset}`.

**VIII.5.** The mapper binary inside the container is invoked via `docker run ... <image> <mapper_command>` with input/output via pipes and bind mounts. No `srun` or Slurm is used.

---

## Article IX — Git Workflow

**IX.1.** All development occurs on feature branches branched from `main`.

**IX.2.** Branch naming: `feature/<description>` or `fix/<description>`.

**IX.3.** Before merging, a feature branch must:
- Pass `snakemake -n` (dry-run) on all new/modified workflows
- Produce output files that integrate with existing downstream workflows without modification
- Include or update relevant documentation
- Have been tested on the `1k` downsample dataset

**IX.4.** Commit messages must be descriptive and reference the feature being implemented.

---

## Article X — Binding Nature and Amendment

**X.1.** This constitution is binding on the AI agent (Cline). It exists to constrain the agent's behavior, prevent hallucinated decisions, and enforce adherence to documented rules and patterns.

**X.2.** The AI agent may not deviate from or override any article in this constitution under any circumstances, unless explicitly instructed to do so by the human user in a direct message.

**X.3.** Human users may amend or override the constitution at any time by direct instruction. No pull request or formal process is required for human-initiated changes.

**X.4.** The constitution is the immutable ground truth for all agent-driven work. All specifications, implementations, and reviews produced by the agent are measured against it.

---

# Nicolás Amendment — SV Caller and Docker Validation Rules

## Purpose

This amendment extends the existing project constitution for the current SV caller and variant-caller workflow validation task. The goal is to validate workflow structure, Docker/container references, and Snakemake dry-run behavior without executing full real-data variant calling unless explicitly approved by the supervisor.

## Article XI — Scope of the Variant Caller Validation Task

1. The validation task focuses on the variant-calling workflows already present in the repository.
2. The caller groups to validate are Clair3, DeepVariant, Sniffles2, cuteSV, pbsv, and Sawfish.
3. Clair3 and DeepVariant are SNV/indel callers.
4. Sniffles2, cuteSV, pbsv, and Sawfish are structural variant callers.
5. The technologies considered are ONT and PacBio HiFi/PacBio-style data, depending on the workflow.
6. Producing final VCF files is outside the scope of this validation stage unless explicitly requested.

## Article XII — Docker and Container Validation Rules

1. The correct term for this task is Docker image or container validation, not Dockerfile testing, unless an actual Dockerfile is being built with docker build.
2. Static Docker validation is mandatory. This means inspecting workflow files for Docker image references, container commands, mounted paths, image tags, and private registry dependencies.
3. Runtime Docker validation is optional in local WSL if Docker is not available. Missing Docker inside WSL must be documented as an environment limitation, not treated as a workflow failure.
4. Public Docker images may be checked locally when Docker is available.
5. Private institutional images under storage-node:5000 must not be marked as failed only because they are inaccessible from the local laptop. They may require the institute workstation, network, or registry access.
6. The Docker image check script must not execute full variant calling. It may only check Docker availability, list image references, optionally pull public images, and optionally run lightweight help or version checks.

Approved Docker validation files:

- SV aligners call/docs/docker_references_static_scan.txt
- SV aligners call/docs/docker_status_note.md
- SV aligners call/scripts/check_variant_caller_images.sh

Public image examples detected in the project:

- hkubal/clair3
- google/deepvariant
- google/deepvariant GPU image
- nfcore/bcl2fastq
- ensemblorg/ensembl-vep

Private registry examples detected in the project:

- storage-node:5000/own/genetic_data_analysis
- storage-node:5000/own/varcad
- storage-node:5000/illumina/hap.py

## Article XIII — Dry-Run-Only Execution Policy

1. For the current task, Snakemake dry-run validation is considered adequate unless the supervisor explicitly requests real execution.
2. The default validation command must include -n -p.
3. Commands without -n must not be executed on real data unless explicitly approved.
4. Dry-run validation is acceptable because it verifies Snakefile syntax, rule connectivity, input and output resolution, Docker command construction, resource declarations, expected workflow DAG, and file naming consistency.
5. Dry-run logs should be saved under SV aligners call/logs/.
6. DAG or rulegraph outputs should be saved under SV aligners call/docs/ or a dedicated DAG folder if created later.

## Article XIV — SV Caller-Specific Rules

1. SV callers consume mapped CRAM or BAM files produced by the read-mapping workflows.
2. SV callers must not remap reads internally unless that behavior is explicitly part of the caller.
3. SV caller outputs must follow the existing project naming logic.
4. The registered SV caller directories are sniffles2, cuteSV, pbsv, and sawfish.
5. SV caller dry-runs must confirm that each caller can resolve the expected mapped inputs from the cram, bam, or equivalent alignment output directory.
6. SV benchmarking with Truvari is downstream of SV calling and must not be mixed with caller execution rules unless the workflow explicitly defines that stage.
7. Caller-specific assumptions must be documented. pbsv and Sawfish are primarily associated with PacBio HiFi/PacBio-style data. Sniffles2 and cuteSV are general long-read SV callers and may be evaluated with ONT and PacBio-derived alignments.

## Article XV — Required Deliverables for This Validation Task

Required files:

- SV aligners call/CONSTITUTION_NICOLAS.md
- SV aligners call/docs/docker_references_static_scan.txt
- SV aligners call/docs/docker_status_note.md
- SV aligners call/scripts/check_variant_caller_images.sh

Recommended additional files:

- SV aligners call/docs/variant_caller_inventory.tsv
- SV aligners call/logs/caller_dryrun_logs
- SV aligners call/docs/caller_DAG_or_rulegraph_outputs

Large outputs, real VCFs, BAMs, CRAMs, indexes, and intermediate variant-calling results must not be committed.

## Article XVI — Git Rules for This Task

1. All work for this task must remain on feature/variant-caller-dryrun-validation.
2. Files must be added intentionally.
3. Do not use git add .
4. Use explicit adds such as git add "SV aligners call".
5. The feature branch may be pushed without merging into main.

## Article XVII — Final Decision Rule

1. If a workflow can be dry-run and its Docker/container references are documented, the validation objective for this stage is satisfied.
2. Real-data execution is out of scope unless requested by the supervisor.
3. Missing local Docker access, private registry restrictions, or unavailable full datasets must be reported transparently as environment limitations, not hidden or bypassed.
