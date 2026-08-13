# Chr21 Long-Read Assemblers

Snakemake workflow for chromosome 21 benchmarking with:

- Flye
- GoldRush
- ntLink
- Verkko

## Supervisor quick start

The workflow is located in:

    assemblers/

For normal server use, the supervisor mainly edits:

    config/server.yaml
    samples.tsv

The rule files under rules/ normally do not need to be changed just to use
different server paths.

## Server input paths

Edit:

    config/server.yaml

Change:

    input_root: "/PATH/TO/WGS_READS"
    reference: "/PATH/TO/GRCh38.fa"

input_root is the directory containing the original whole-genome FASTQs.

reference is the GRCh38 FASTA used for chromosome 21 mapping.

## Samples

samples.tsv contains metadata only:

    sample  technology
    HG002   ont
    HG002   pb
    HG003   ont
    HG003   pb
    HG004   ont
    HG004   pb

Do not put FASTQ paths or server directories in samples.tsv.

Technology values:

    ont = Oxford Nanopore
    pb  = PacBio HiFi

Verkko requires both ONT and PB for each biological sample.

## Default input filenames

The default server filename pattern is:

    {sample}.{technology}.30x.fastq.gz

Examples:

    HG002.ont.30x.fastq.gz
    HG002.pb.30x.fastq.gz
    HG003.ont.30x.fastq.gz
    HG003.pb.30x.fastq.gz
    HG004.ont.30x.fastq.gz
    HG004.pb.30x.fastq.gz

If server filenames differ, override whole_genome_pattern in
config/server.yaml.

Example:

    input_root: "/data/reads"
    whole_genome_pattern: "{sample}/{sample}.{technology}.fastq.gz"

There is no need to modify the Snakemake rules just to change input paths.

## Server preprocessing

Server mode automatically performs:

    whole-genome FASTQ
        ->
    minimap2 mapping
        ->
    primary Chr21 + MAPQ filtering
        ->
    one primary alignment per read
        ->
    Chr21 BAM validation
        ->
    deterministic normalization to 30x
        ->
    results/chr21/{sample}.{technology}.fastq.gz
        ->
    assemblers

Chr21 preprocessing is performed once per sample/technology pair.

## Final outputs

Flye:

    results/flye/{sample}.{technology}.fasta

GoldRush:

    results/goldrush/{sample}.{technology}.fasta

ntLink:

    results/ntlink/{sample}.{technology}.fasta

Verkko:

    results/verkko/{sample}.fasta

ntLink uses the matching GoldRush assembly as its draft.

Verkko combines PacBio HiFi and ONT data from the same sample.

## Docker images

    Flye      nicolasardila1/lrs-flye2:2.9.6
    GoldRush  nicolasardila1/lrs-goldrush:1.2.2
    ntLink    nicolasardila1/lrs-ntlink:1.3.11
    Verkko    nicolasardila1/lrs-verkko2:2.3.2

## Server prerequisites

Check:

    command -v git
    command -v python3
    command -v snakemake
    command -v docker
    command -v minimap2
    command -v samtools
    docker info

## Server dry-run

From the assemblers directory:

    snakemake \
      --configfiles config/server.yaml \
      --cores 128 \
      --resources mem_mb=240000 \
      --dry-run \
      --printshellcmds

With the current samples, the expected server DAG contains 41 jobs.

Do not start the real workflow if the dry-run shows unexpected paths.

## First server test

Test one Flye path first:

    snakemake \
      results/flye/HG002.ont.fasta \
      --configfiles config/server.yaml \
      --cores 32 \
      --resources mem_mb=64000 \
      --rerun-incomplete \
      --printshellcmds \
      --show-failed-logs

Then test Verkko:

    snakemake \
      results/verkko/HG002.fasta \
      --configfiles config/server.yaml \
      --cores 32 \
      --resources mem_mb=72000 \
      --rerun-incomplete \
      --printshellcmds \
      --show-failed-logs

## Full workflow

After targeted tests pass:

    snakemake \
      --configfiles config/server.yaml \
      --cores 128 \
      --resources mem_mb=240000 \
      --rerun-incomplete \
      --printshellcmds \
      --show-failed-logs

## Current status

    READY FOR SERVER TESTING = YES
    FULL SERVER VALIDATION   = PENDING

The final end-to-end validation must be performed on the server with the real
Chr21 datasets.
