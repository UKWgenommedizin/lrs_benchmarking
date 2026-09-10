# Long-Read Assembly Benchmarking

This README documents the whole-genome long-read assembly benchmarking workflow
for HG002, HG003 and HG004 using Oxford Nanopore (ONT) and PacBio HiFi data.

It covers:

- the existing Flye, GoldRush, Verkko and ntLink workflows;
- the standardized assembly output layout;
- the new QUAST-LG assembly-quality assessment workflow;
- commands for running the existing assembler workflows;
- commands for running and checking QUAST-LG;
- recommended quality metrics and interpretation notes.

This README is intentionally separate from `alignment_analysis/README.md`, which
documents the long-read **alignment** benchmark for minimap2, pbmm2, VACMap and
VG Giraffe.

---

## 1. Current clean QUAST pull request

The current clean QUAST pull request introduces or updates only:

```text
alignment_analysis/README_assemblers.md
assemblers/whole_genome_asm/assessment/assembly_quality_quast.smk
```

The pull request does **not** modify the Flye, GoldRush, Verkko, ntLink or
aligner workflows. Those files are documented below because they produce the
assemblies that QUAST evaluates.

---

## 2. Benchmark overview

The assembly benchmark uses:

```text
Samples:
HG002
HG003
HG004

Read technologies:
ONT
PacBio HiFi
```

Assembly/scaffolding methods:

| Method | Input | Role |
|---|---|---|
| Flye | ONT or PacBio HiFi | Independent long-read assembly |
| GoldRush | ONT or PacBio HiFi | Independent long-read assembly |
| Verkko | ONT + PacBio HiFi | Hybrid long-read assembly |
| ntLink | GoldRush draft + long reads | Scaffolding / post-assembly refinement |
| QUAST-LG | Completed assembly + reference | Assembly-quality assessment |

`ntLink` is treated as a scaffolder rather than as an independent assembler.

---

## 3. Relevant repository files

```text
lrs_benchmarking/
├── alignment_analysis/
│   ├── README.md
│   └── README_assemblers.md
│
├── assemblers/
│   ├── README.md
│   ├── Snakefile
│   ├── config/
│   │   ├── base.yaml
│   │   ├── local.yaml
│   │   └── server.yaml
│   │
│   └── whole_genome_asm/
│       ├── ont.assembly.flye2.smk
│       ├── pb.assembly.flye2.smk
│       ├── ont.assembly.Goldrush.smk
│       ├── pb.assembly.Goldrush.smk
│       ├── ont.assembly.ntlink.smk
│       ├── pb.assembly.ntlink.smk
│       ├── hybrid.assembly.verkko.smk
│       └── assessment/
│           └── assembly_quality_quast.smk
│
├── fastq/
├── assemblies/
└── assembly_quality/
```

The large FASTQ files, assembly outputs and QUAST results are normally generated
on the server and should not be pushed to GitHub unless they are intentionally
small test data or summary outputs.

---

## 4. Current tool versions and Docker images

The current whole-genome workflows use the following versions/images:

| Tool | Version | Docker image |
|---|---:|---|
| Flye | 2.9.6 | `nicolasardila1/lrs-flye2:2.9.6` |
| GoldRush | 1.2.2-ntlinkfix | `nicolasardila1/lrs-goldrush:1.2.2-ntlinkfix` |
| ntLink | 1.3.11 | `nicolasardila1/lrs-ntlink:1.3.11` |
| Verkko | 2.3.2 | `nicolasardila1/lrs-verkko2:2.3.2` |
| QUAST-LG | 5.3.0 | `quay.io/biocontainers/quast:5.3.0--py313pl5321h5ca1c30_2` |

Keep versions and container images pinned when producing benchmark results.

---

## 5. Production input naming convention

Whole-genome production FASTQs use:

```text
{sample}.{technology}.30x.fastq.gz
```

Examples:

```text
HG002.ont.30x.fastq.gz
HG002.pb.30x.fastq.gz
HG003.ont.30x.fastq.gz
HG003.pb.30x.fastq.gz
HG004.ont.30x.fastq.gz
HG004.pb.30x.fastq.gz
```

Use:

```text
ont = Oxford Nanopore
pb  = PacBio HiFi
```

Do not mix production 30x results with reduced test datasets such as:

```text
*.1k.*
*.chr21.*
*SMOKE*
*.localtest.*
```

---

## 6. Standardized assembly outputs

The standalone whole-genome workflows write assemblies under:

```text
assemblies/{assembler}/{dataset}/assembly.fasta
```

Examples:

```text
assemblies/flye/HG002.ont.30x/assembly.fasta
assemblies/flye/HG002.pb.30x/assembly.fasta

assemblies/goldrush/HG002.ont.30x/assembly.fasta
assemblies/goldrush/HG002.pb.30x/assembly.fasta

assemblies/ntlink/HG002.ont.30x/assembly.fasta
assemblies/ntlink/HG002.pb.30x/assembly.fasta
```

Verkko uses paired ONT + PacBio data and produces one assembly per sample:

```text
assemblies/verkko/HG002/assembly.fasta
assemblies/verkko/HG003/assembly.fasta
assemblies/verkko/HG004/assembly.fasta
```

The new QUAST workflow discovers completed assemblies from this standardized
layout.

---

# 7. Running the existing whole-genome assembler workflows

Run the standalone production workflows from the repository root:

```bash
cd /path/to/lrs_benchmarking
```

On the server this should be replaced by the actual repository path.

Before running any production workflow:

```bash
command -v snakemake
command -v docker
docker info
```

A dry run is strongly recommended first.

---

## 7.1 Flye — ONT

Workflow:

```text
assemblers/whole_genome_asm/ont.assembly.flye2.smk
```

Dry run:

```bash
snakemake \
  --snakefile assemblers/whole_genome_asm/ont.assembly.flye2.smk \
  --dry-run \
  --printshellcmds
```

Run:

```bash
snakemake \
  --snakefile assemblers/whole_genome_asm/ont.assembly.flye2.smk \
  --cores 32 \
  --resources mem_gb=240 \
  --rerun-incomplete \
  --printshellcmds \
  --show-failed-logs
```

The ONT Flye workflow uses:

```text
Flye 2.9.6
--nano-hq
32 threads
up to 240 GB RAM for whole-genome production data
```

---

## 7.2 Flye — PacBio HiFi

Workflow:

```text
assemblers/whole_genome_asm/pb.assembly.flye2.smk
```

Dry run:

```bash
snakemake \
  --snakefile assemblers/whole_genome_asm/pb.assembly.flye2.smk \
  --dry-run \
  --printshellcmds
```

Run:

```bash
snakemake \
  --snakefile assemblers/whole_genome_asm/pb.assembly.flye2.smk \
  --cores 32 \
  --resources mem_gb=240 \
  --rerun-incomplete \
  --printshellcmds \
  --show-failed-logs
```

The PacBio Flye workflow uses:

```text
Flye 2.9.6
--pacbio-hifi
32 threads
up to 240 GB RAM for whole-genome production data
```

---

## 7.3 GoldRush — ONT

Workflow:

```text
assemblers/whole_genome_asm/ont.assembly.Goldrush.smk
```

Dry run:

```bash
snakemake \
  --snakefile assemblers/whole_genome_asm/ont.assembly.Goldrush.smk \
  --dry-run \
  --printshellcmds
```

Run:

```bash
snakemake \
  --snakefile assemblers/whole_genome_asm/ont.assembly.Goldrush.smk \
  --cores 32 \
  --resources mem_mb=64000 \
  --rerun-incomplete \
  --printshellcmds \
  --show-failed-logs
```

Current GoldRush configuration:

```text
GoldRush image: nicolasardila1/lrs-goldrush:1.2.2-ntlinkfix
Threads: 32
Memory: 64 GB
Genome size: 3e9
ONT m parameter: 5000
```

GoldRush requires an uncompressed FASTQ internally. The workflow prepares this
from the gzipped production input.

The workflow also removes previous internal
`goldrush_intermediate_files/` state before a new run so stale Make/ntLink
checkpoints are not reused.

---

## 7.4 GoldRush — PacBio HiFi

Workflow:

```text
assemblers/whole_genome_asm/pb.assembly.Goldrush.smk
```

Dry run:

```bash
snakemake \
  --snakefile assemblers/whole_genome_asm/pb.assembly.Goldrush.smk \
  --dry-run \
  --printshellcmds
```

Run:

```bash
snakemake \
  --snakefile assemblers/whole_genome_asm/pb.assembly.Goldrush.smk \
  --cores 32 \
  --resources mem_mb=64000 \
  --rerun-incomplete \
  --printshellcmds \
  --show-failed-logs
```

Current PacBio GoldRush configuration:

```text
GoldRush image: nicolasardila1/lrs-goldrush:1.2.2-ntlinkfix
Threads: 32
Memory: 64 GB
Genome size: 3e9
PacBio m parameter: 10000
```

---

## 7.5 ntLink — ONT

Workflow:

```text
assemblers/whole_genome_asm/ont.assembly.ntlink.smk
```

ntLink uses the corresponding GoldRush assembly as its draft.

Run:

```bash
snakemake \
  --snakefile assemblers/whole_genome_asm/ont.assembly.ntlink.smk \
  --cores 16 \
  --resources mem_mb=32000 \
  --rerun-incomplete \
  --printshellcmds \
  --show-failed-logs
```

Example inputs:

```text
assemblies/goldrush/HG002.ont.30x/assembly.fasta
fastq/HG002.ont.30x.fastq.gz
```

Output:

```text
assemblies/ntlink/HG002.ont.30x/assembly.fasta
```

Current configuration:

```text
ntLink 1.3.11
Threads: 16
Memory: 32 GB
k=32
w=100
z=1000
```

---

## 7.6 ntLink — PacBio HiFi

Workflow:

```text
assemblers/whole_genome_asm/pb.assembly.ntlink.smk
```

Run:

```bash
snakemake \
  --snakefile assemblers/whole_genome_asm/pb.assembly.ntlink.smk \
  --cores 16 \
  --resources mem_mb=32000 \
  --rerun-incomplete \
  --printshellcmds \
  --show-failed-logs
```

Example inputs:

```text
assemblies/goldrush/HG002.pb.30x/assembly.fasta
fastq/HG002.pb.30x.fastq.gz
```

Output:

```text
assemblies/ntlink/HG002.pb.30x/assembly.fasta
```

---

## 7.7 Verkko — hybrid ONT + PacBio HiFi

Workflow:

```text
assemblers/whole_genome_asm/hybrid.assembly.verkko.smk
```

Verkko requires both technologies for the same sample.

Required example inputs:

```text
fastq/HG002.ont.30x.fastq.gz
fastq/HG002.pb.30x.fastq.gz
```

Dry run:

```bash
snakemake \
  --snakefile assemblers/whole_genome_asm/hybrid.assembly.verkko.smk \
  --configfile assemblers/config/server.yaml \
  --dry-run \
  --printshellcmds
```

Run:

```bash
snakemake \
  --snakefile assemblers/whole_genome_asm/hybrid.assembly.verkko.smk \
  --configfile assemblers/config/server.yaml \
  --cores 32 \
  --resources mem_mb=200000 \
  --rerun-incomplete \
  --printshellcmds \
  --show-failed-logs
```

Current configuration:

```text
Verkko 2.3.2
Threads: 32
Snakemake memory resource: 200 GB
Verkko --local-memory: 64 GB
```

Output example:

```text
assemblies/verkko/HG002/assembly.fasta
```

---

# 8. New QUAST-LG assembly-quality workflow

New workflow:

```text
assemblers/whole_genome_asm/assessment/assembly_quality_quast.smk
```

This workflow evaluates **already completed assemblies**. It does not run or
modify Flye, GoldRush, Verkko or ntLink.

It automatically searches for:

```text
assemblies/{assembler}/{dataset}/assembly.fasta
```

and accepts the following assembly output directories:

```text
flye
goldrush
verkko
ntlink
```

---

## 8.1 QUAST configuration

Current configuration:

```text
QUAST version: 5.3.0

Docker image:
quay.io/biocontainers/quast:5.3.0--py313pl5321h5ca1c30_2

Mode:
--large

Threads:
16

Memory:
128 GB

Temporary Docker filesystem:
50 GB

Minimum contig length:
500 bp
```

The workflow performs reference-based assessment against the project GRCh38
reference.

It contains a project default reference path, but the recommended execution is
to provide the exact reference explicitly through:

```text
--config reference=/path/to/reference.fasta
```

This avoids accidentally using a different reference.

---

## 8.2 Check completed assemblies before QUAST

From the repository root:

```bash
find assemblies \
  -mindepth 3 \
  -maxdepth 3 \
  -type f \
  -name assembly.fasta \
  -size +0c \
  -print | sort
```

A local reduced test assembly may also be used for workflow validation if it is
already present in the standardized layout, for example:

```text
assemblies/flye/HG002.ont.1k/assembly.fasta
```

Do not compare the resulting 1k test metrics scientifically with production 30x
whole-genome assemblies.

---

## 8.3 QUAST dry run

```bash
snakemake \
  --snakefile assemblers/whole_genome_asm/assessment/assembly_quality_quast.smk \
  --cores 16 \
  --resources mem_gb=128 \
  --config reference=/path/to/GRCh38_reference.fasta \
  --dry-run \
  --printshellcmds
```

The dry run should list one QUAST target for each completed assembly discovered
under `assemblies/`.

---

## 8.4 Run QUAST for all discovered assemblies

```bash
snakemake \
  --snakefile assemblers/whole_genome_asm/assessment/assembly_quality_quast.smk \
  --cores 16 \
  --resources mem_gb=128 \
  --config reference=/path/to/GRCh38_reference.fasta \
  --rerun-incomplete \
  --printshellcmds \
  --show-failed-logs
```

Replace:

```text
/path/to/GRCh38_reference.fasta
```

with the exact GRCh38 FASTA used for the benchmark.

---

## 8.5 Run QUAST for one assembly

Example for Flye HG002 ONT:

```bash
snakemake \
  --snakefile assemblers/whole_genome_asm/assessment/assembly_quality_quast.smk \
  "assembly_quality/quast/flye/HG002.ont.30x/report.tsv" \
  --cores 16 \
  --resources mem_gb=128 \
  --config reference=/path/to/GRCh38_reference.fasta \
  --rerun-incomplete \
  --printshellcmds \
  --show-failed-logs
```

Change the target path to run QUAST for another assembler or dataset.

---

## 8.6 QUAST outputs

Each completed assembly produces:

```text
assembly_quality/quast/{assembler}/{dataset}/report.tsv
assembly_quality/quast/{assembler}/{dataset}/quast.log
```

Examples:

```text
assembly_quality/quast/flye/HG002.ont.30x/report.tsv
assembly_quality/quast/goldrush/HG002.ont.30x/report.tsv
assembly_quality/quast/ntlink/HG002.ont.30x/report.tsv
assembly_quality/quast/verkko/HG002/report.tsv
```

The workflow fails if `report.tsv` is missing, empty, or does not contain N50.

---

## 8.7 Check QUAST results

List all completed reports:

```bash
find assembly_quality/quast \
  -type f \
  -name report.tsv \
  -size +0c \
  -print | sort
```

Open one report:

```bash
column -t -s $'	' \
  assembly_quality/quast/flye/HG002.ont.30x/report.tsv | less -S
```

Open its log:

```bash
less assembly_quality/quast/flye/HG002.ont.30x/quast.log
```

Show the main comparison metrics:

```bash
grep -E \
'^(# contigs|Largest contig|Total length|N50|N90|L50|L90|NG50|NGA50|Genome fraction|Duplication ratio|# misassemblies|# mismatches per 100 kbp|# indels per 100 kbp)' \
assembly_quality/quast/flye/HG002.ont.30x/report.tsv
```

Use the generated `report.tsv` as the source of truth.

---

# 9. Main assembly-quality metrics

For the final assembler comparison, focus on:

| Category | Main metrics |
|---|---|
| Contiguity | `# contigs`, largest contig, total length, N50/L50, N90/L90, NG50/NGA50 |
| Reference agreement | genome fraction, misassemblies, duplication ratio |
| Base-level accuracy | mismatches per 100 kbp, indels per 100 kbp |

Do not rank assemblies using N50 alone.

Additional analyses such as BUSCO, k-mer completeness, Merqury QV, runtime and
peak RAM should be generated separately when available.

---

# 10. Recommended order

```text
30x FASTQ
   ↓
Flye / GoldRush / Verkko
   ↓
ntLink where applicable
   ↓
assembly.fasta
   ↓
QUAST-LG
   ↓
report.tsv
   ↓
figures / statistical comparison
```

In practice:

```text
1. Confirm the production assemblies exist.
2. Dry-run the QUAST workflow.
3. Run QUAST-LG with the project GRCh38 reference.
4. Check report.tsv and quast.log.
5. Extract the selected metrics for comparison.
```

---

# 11. Quick validation commands

Check completed assemblies:

```bash
find assemblies \
  -type f \
  -name assembly.fasta \
  -size +0c \
  -print | sort
```

Check completed QUAST reports:

```bash
find assembly_quality/quast \
  -type f \
  -name report.tsv \
  -size +0c \
  -print | sort
```

---

# 12. Reproducibility

- Run the workflows from the repository root.
- Keep assembler and QUAST Docker images pinned.
- Use the same GRCh38 reference for every QUAST comparison.
- Do not mix reduced test datasets with production 30x results.
- Keep the raw `report.tsv` files and logs.
- Report missing values as missing, not zero.
- Treat ntLink as a scaffolding/post-assembly step.

---

# 13. Files directly changed by the current QUAST pull request

```text
alignment_analysis/README_assemblers.md
assemblers/whole_genome_asm/assessment/assembly_quality_quast.smk
```

The other assembler workflow files described above are existing project files
and are not modified by this clean QUAST pull request.
