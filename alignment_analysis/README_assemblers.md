# Long-Read Assembly Benchmarking

This directory contains the Snakemake workflows used for whole-genome long-read
assembly benchmarking in the `lrs_benchmarking` project.

The production benchmark uses HG002, HG003 and HG004 with Oxford Nanopore (ONT)
and PacBio HiFi data. Flye, GoldRush and Verkko are used for assembly, ntLink
for post-assembly scaffolding, and QUAST-LG for assembly-quality assessment.

For read-mapping benchmarking with minimap2, pbmm2, VACMap and VG Giraffe, see
`alignment_analysis/README.md`.

## Repository structure

```text
lrs_benchmarking/
├── assemblers/
│   ├── README.md
│   ├── Snakefile
│   ├── samples.tsv
│   ├── config/
│   │   ├── base.yaml
│   │   ├── local.yaml
│   │   └── server.yaml
│   ├── rules/
│   │   ├── common.smk
│   │   ├── validation.smk
│   │   ├── chr21.smk
│   │   ├── flye.smk
│   │   ├── goldrush.smk
│   │   ├── ntlink.smk
│   │   ├── verkko.smk
│   │   └── assessment.smk
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
├── fastq/
├── assemblies/
└── assembly_quality/
```

The workflows under `assemblers/whole_genome_asm/` are the production
whole-genome workflows. The modular `assemblers/Snakefile` workflow is retained
for prepared-Chr21/local validation.

## Benchmark design

| Method | Input | Role |
|---|---|---|
| Flye | ONT or PacBio HiFi | Long-read assembly |
| GoldRush | ONT or PacBio HiFi | Long-read assembly |
| Verkko | ONT + PacBio HiFi | Hybrid long-read assembly |
| ntLink | GoldRush assembly + long reads | Scaffolding |
| QUAST-LG | Completed assembly + GRCh38 | Assembly-quality assessment |

`ntLink` is not treated as an independent assembler.

Production FASTQs use:

```text
HG002.ont.30x.fastq.gz
HG002.pb.30x.fastq.gz
HG003.ont.30x.fastq.gz
HG003.pb.30x.fastq.gz
HG004.ont.30x.fastq.gz
HG004.pb.30x.fastq.gz
```

Reduced `.1k`, `.chr21`, smoke-test and `.localtest` datasets are for validation
only and should not be mixed with the production 30x comparison.

## Tool versions

| Tool | Version | Docker image |
|---|---:|---|
| Flye | 2.9.6 | `nicolasardila1/lrs-flye2:2.9.6` |
| GoldRush | 1.2.2-ntlinkfix | `nicolasardila1/lrs-goldrush:1.2.2-ntlinkfix` |
| ntLink | 1.3.11 | `nicolasardila1/lrs-ntlink:1.3.11` |
| Verkko | 2.3.2 | `nicolasardila1/lrs-verkko2:2.3.2` |
| QUAST-LG | 5.3.0 | `quay.io/biocontainers/quast:5.3.0--py313pl5321h5ca1c30_2` |

## Production output layout

```text
assemblies/{assembler}/{dataset}/assembly.fasta
```

Examples:

```text
assemblies/flye/HG002.ont.30x/assembly.fasta
assemblies/goldrush/HG002.ont.30x/assembly.fasta
assemblies/ntlink/HG002.ont.30x/assembly.fasta
assemblies/verkko/HG002/assembly.fasta
```

## Running the whole-genome workflows

Run from the repository root:

```bash
cd /path/to/lrs_benchmarking
```

Check the environment:

```bash
command -v snakemake
command -v docker
docker info
```

### Flye

ONT:

```bash
snakemake   --snakefile assemblers/whole_genome_asm/ont.assembly.flye2.smk   --cores 32   --resources mem_gb=240   --rerun-incomplete   --printshellcmds   --show-failed-logs
```

PacBio HiFi:

```bash
snakemake   --snakefile assemblers/whole_genome_asm/pb.assembly.flye2.smk   --cores 32   --resources mem_gb=240   --rerun-incomplete   --printshellcmds   --show-failed-logs
```

### GoldRush

ONT:

```bash
snakemake   --snakefile assemblers/whole_genome_asm/ont.assembly.Goldrush.smk   --cores 32   --resources mem_mb=64000   --rerun-incomplete   --printshellcmds   --show-failed-logs
```

PacBio HiFi:

```bash
snakemake   --snakefile assemblers/whole_genome_asm/pb.assembly.Goldrush.smk   --cores 32   --resources mem_mb=64000   --rerun-incomplete   --printshellcmds   --show-failed-logs
```

### ntLink

ONT:

```bash
snakemake   --snakefile assemblers/whole_genome_asm/ont.assembly.ntlink.smk   --cores 16   --resources mem_mb=32000   --rerun-incomplete   --printshellcmds   --show-failed-logs
```

PacBio HiFi:

```bash
snakemake   --snakefile assemblers/whole_genome_asm/pb.assembly.ntlink.smk   --cores 16   --resources mem_mb=32000   --rerun-incomplete   --printshellcmds   --show-failed-logs
```

### Verkko

```bash
snakemake   --snakefile assemblers/whole_genome_asm/hybrid.assembly.verkko.smk   --configfile assemblers/config/server.yaml   --cores 32   --resources mem_mb=200000   --rerun-incomplete   --printshellcmds   --show-failed-logs
```

## QUAST-LG assembly-quality assessment

Workflow:

```text
assemblers/whole_genome_asm/assessment/assembly_quality_quast.smk
```

It discovers completed assemblies under:

```text
assemblies/{assembler}/{dataset}/assembly.fasta
```

and writes:

```text
assembly_quality/quast/{assembler}/{dataset}/report.tsv
assembly_quality/quast/{assembler}/{dataset}/quast.log
```

Current settings:

```text
QUAST-LG: 5.3.0
Threads: 16
Memory: 128 GB
Mode: --large
Minimum contig length: 500 bp
Docker tmpfs: 50 GB
```

Check available assemblies:

```bash
find assemblies -type f -name assembly.fasta -size +0c -print | sort
```

Dry run:

```bash
snakemake   --snakefile assemblers/whole_genome_asm/assessment/assembly_quality_quast.smk   --cores 16   --resources mem_gb=128   --config reference=/path/to/GRCh38_reference.fasta   --dry-run   --printshellcmds
```

Run all discovered assemblies:

```bash
snakemake   --snakefile assemblers/whole_genome_asm/assessment/assembly_quality_quast.smk   --cores 16   --resources mem_gb=128   --config reference=/path/to/GRCh38_reference.fasta   --rerun-incomplete   --printshellcmds   --show-failed-logs
```

Run one assembly, for example Flye HG002 ONT:

```bash
snakemake   --snakefile assemblers/whole_genome_asm/assessment/assembly_quality_quast.smk   "assembly_quality/quast/flye/HG002.ont.30x/report.tsv"   --cores 16   --resources mem_gb=128   --config reference=/path/to/GRCh38_reference.fasta   --rerun-incomplete   --printshellcmds
```

Check completed reports:

```bash
find assembly_quality/quast -type f -name report.tsv -size +0c -print | sort
```

Open one report:

```bash
column -t -s $'	'   assembly_quality/quast/flye/HG002.ont.30x/report.tsv | less -S
```

Open its log:

```bash
less assembly_quality/quast/flye/HG002.ont.30x/quast.log
```

The workflow fails if `report.tsv` is missing, empty or does not contain N50.

## Main assembly-quality metrics

| Category | Metrics |
|---|---|
| Contiguity | number of contigs, largest contig, total length, N50/L50, N90/L90, NG50/NGA50 |
| Reference agreement | genome fraction, misassemblies, duplication ratio |
| Base-level accuracy | mismatches per 100 kbp, indels per 100 kbp |

N50 should not be interpreted alone.

## Recommended analysis order

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
figures and statistical comparison
```

## Prepared-Chr21 / local validation

The modular workflow under:

```text
assemblers/Snakefile
assemblers/config/
assemblers/rules/
```

is retained for prepared-Chr21 benchmarking and local validation.

Example dry run:

```bash
cd assemblers

snakemake   --snakefile Snakefile   --configfiles config/local.yaml   --cores 4   --resources mem_mb=12000   --dry-run   --printshellcmds
```

## Reproducibility notes

- Run production whole-genome workflows from the repository root.
- Keep Docker images and tool versions pinned.
- Use the same GRCh38 reference for every QUAST comparison.
- Do not mix reduced validation datasets with production 30x results.
- Keep raw QUAST reports and logs.
- Treat missing values as missing rather than zero.
- Treat ntLink as a scaffolding/post-assembly step.
