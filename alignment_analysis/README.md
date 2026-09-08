# Long-Read Assembly Benchmarking

This directory documents the whole-genome long-read assembly benchmark. The
benchmark compares Flye, GoldRush, Verkko and ntLink using ONT and PacBio HiFi
data from HG002, HG003 and HG004.

This README concerns assembly results only. It does not describe read-mapping
workflows or alignment statistics.

## Assembly methods

| Assembler | Input | Description |
|---|---|---|
| Flye | ONT or PacBio HiFi | Independent long-read assembly for each technology. |
| GoldRush | ONT or PacBio HiFi | Independent long-read assembly for each technology. |
| ntLink | ONT or PacBio HiFi plus GoldRush draft | Scaffolding/refinement of the corresponding GoldRush assembly. |
| Verkko | ONT and PacBio HiFi | Hybrid assembly using both technologies for the same sample. |

## Repository structure

```text
lrs_benchmarking_wgs_clean/
├── assemblers/
│   ├── whole_genome_asm/
│   │   ├── ont.assembly.flye2.smk
│   │   ├── pb.assembly.flye2.smk
│   │   ├── ont.assembly.Goldrush.smk
│   │   ├── pb.assembly.Goldrush.smk
│   │   ├── ont.assembly.ntlink.smk
│   │   ├── pb.assembly.ntlink.smk
│   │   └── hybrid.assembly.verkko.smk
│   ├── config/
│   └── README.md
├── fastq/
│   ├── HG002.ont.30x.fastq.gz
│   ├── HG002.pb.30x.fastq.gz
│   ├── HG003.ont.30x.fastq.gz
│   ├── HG003.pb.30x.fastq.gz
│   ├── HG004.ont.30x.fastq.gz
│   └── HG004.pb.30x.fastq.gz
├── assemblies/
│   ├── flye/
│   ├── goldrush/
│   ├── ntlink/
│   └── verkko/
├── results_assemblers/
│   ├── quality_metrics/
│   ├── run_metrics/
│   └── scripts/
└── figures/
```

The FASTQ files and assembly outputs are normally kept on the server and are
not pushed to GitHub. The repository should contain workflows, scripts,
configuration files, documentation and appropriately sized summary tables.

## Input naming convention

Use the following pattern for whole-genome input files:

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

Use `ont` for Oxford Nanopore and `pb` for PacBio HiFi. Verkko requires both
technologies for the same sample. Files containing `.1k`, `.chr21`, `smoke` or
`localtest` are reduced test data and must not be confused with production
whole-genome 30x inputs.

## Assembly outputs

The expected output pattern is:

```text
assemblies/{assembler}/{sample}.{technology}.30x/assembly.fasta
```

For Verkko, the output is sample-level because both technologies are used:

```text
assemblies/verkko/{sample}.30x/assembly.fasta
```

ntLink uses the matching GoldRush assembly as its draft and should be run only
after the corresponding GoldRush assembly is available.

## Assembly-quality metrics

The combined table should contain one row per assembly and the following
parameters:

| Parameter | Meaning |
|---|---|
| `sample` | HG002, HG003 or HG004. |
| `technology` | ONT, PacBio HiFi or hybrid. |
| `assembler` | Flye, GoldRush, ntLink or Verkko. |
| `assembly_file` | Path to the assembly FASTA. |
| `number_of_sequences` | Number of contigs or scaffolds. |
| `total_length` | Total assembly length in base pairs. |
| `largest_sequence` | Length of the largest contig or scaffold. |
| `smallest_sequence` | Length of the smallest contig or scaffold. |
| `mean_length` | Mean sequence length. |
| `median_length` | Median sequence length. |
| `N50` | Length at which 50% of the assembly is contained in sequences of that length or longer. |
| `L50` | Number of sequences required to reach 50% of the assembly length. |
| `N90` | Length at which 90% of the assembly is contained in sequences of that length or longer. |
| `L90` | Number of sequences required to reach 90% of the assembly length. |
| `A`, `C`, `G`, `T`, `N` | Base composition counts. |
| `runtime_seconds` | Wall-clock assembly time, when available. |
| `peak_ram_mb` | Maximum memory used, when available. |
| `threads` | Number of threads used. |
| `software_version` | Version of the assembler. |

N50 and L50 should always be interpreted together with total length and the
number of sequences. A high N50 alone does not prove that an assembly is more
complete or more accurate.

## Running the assembly workflows

Run production workflows from the repository root:

```bash
cd /path/to/lrs_benchmarking_wgs_clean
```

Perform a dry run before submitting a production job:

```bash
snakemake \
    --snakefile assemblers/whole_genome_asm/ont.assembly.flye2.smk \
    --dry-run \
    --printshellcmds
```

The production workflows are run according to the project execution plan. The
workflow commands are kept separate from the quality-metric extraction command
because assembly generation and assembly assessment are different stages.

## Extracting metrics and creating one combined table

After the assemblies have been generated, run the quality-metric extraction
from the repository root. The extractor should scan the `assemblies/` directory
recursively, identify Flye, GoldRush, ntLink and Verkko FASTA files, and write
one table containing all detected assemblies.

```bash
cd /path/to/lrs_benchmarking_wgs_clean

python3 results_assemblers/scripts/build_assembly_summary.py \
    --assemblies assemblies \
    --out results_assemblers/quality_metrics/assembly_quality_summary.tsv
```

The final table is:

```text
results_assemblers/quality_metrics/assembly_quality_summary.tsv
```

If the current extractor is named
`results_assemblers/quality_metrics/calculate_quality_metrics.py` and accepts
only one FASTA at a time, it must first be extended with a directory-scanning
or table-aggregation mode. The desired behavior is one command that scans all
assembly FASTA files and writes one TSV table, rather than manually creating a
separate table for every assembler.

## Checking the combined table

```bash
test -s results_assemblers/quality_metrics/assembly_quality_summary.tsv
head -n 2 results_assemblers/quality_metrics/assembly_quality_summary.tsv
column -t -s $'\t' \
    results_assemblers/quality_metrics/assembly_quality_summary.tsv | less -S
```

Check the detected assemblies with:

```bash
python3 - <<'PY'
import pandas as pd

path = "results_assemblers/quality_metrics/assembly_quality_summary.tsv"
table = pd.read_csv(path, sep="\t")
print(table.shape)
print(table[["sample", "technology", "assembler"]].to_string(index=False))
PY
```

The final TSV is the source of truth for figures and summary tables. Missing
metrics should be recorded as `NA`, not as zero. Runtime and RAM should only be
reported when they were measured by the workflow or scheduler.


