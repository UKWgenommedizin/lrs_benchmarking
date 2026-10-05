# Assembler RAM / Time Measurement

Companion to [`../README.md`](../README.md). This folder only measures
**wall-clock time and peak RAM** for the whole-genome assemblers. It never
produces assemblies used in the report.

```text
ont.assembly.flye2.ram_time.smk
pb.assembly.flye2.ram_time.smk
ont.assembly.goldrush.ram_time.smk
pb.assembly.goldrush.ram_time.smk
hybrid.assembly.verkko.ram_time.smk
```

## Why this folder exists

The production workflows in the parent directory (`ont.assembly.flye2.smk`,
`pb.assembly.flye2.smk`, `ont.assembly.goldrush.smk`,
`pb.assembly.goldrush.smk`, `hybrid.assembly.verkko.smk`) never recorded
time or RAM for the 30x whole-genome runs, and no full-run log survives to
recover them (see the module docstring of
`assembly_analysis/scripts/utils/benchmark_data.py`).

Each `*.ram_time.smk` re-runs its production counterpart with the **same
Docker image, version, memory ceiling and assembler parameters**. The only
difference is the thread count: **64 threads** by default (`ram_time_threads`)
instead of the production 32, so the numbers describe a 64-thread run.

## How time and RAM are measured

Snakemake's `benchmark:` directive is not used: it samples the host process
tree and cannot see memory used inside a Docker container, so it would
report a near-zero RAM value. Instead:

- **Flye and Verkko**: the assembler is wrapped in `/usr/bin/time -v`
  *inside* the container (`--entrypoint /usr/bin/time`). Before the
  multi-hour assembly starts, the run step checks that the image contains
  `/usr/bin/time` and `python3` and fails within seconds otherwise.
- **GoldRush**: already times every internal stage with `time -v` when
  `track_time=1` is set (as in production). The record step **sums** the
  elapsed time of all stages and takes the **maximum** peak RSS across them.
  The run step checks that the image contains `gzip` and `python3` up front,
  and fails if the final `*ntLink-5rounds.polished.fa` is missing, so a
  partial run is never reported as a full one.

## Rules: run -> record

| Rule | What it does | Kept? |
|---|---|---|
| `<tool>_run` | Runs the assembler inside a scratch directory declared as `temp(directory(...))` | No |
| `<tool>_record` | Parses the time file(s) into `{dataset}.ram_time.tsv` with `assembly_analysis/scripts/metrics/parse_assembler_time_v.py`, run with `python3` inside the assembler's own image (CONSTITUTION II.1) | **Yes** |

- The scratch directory (assembly, intermediates, time files, GoldRush's
  decompressed FASTQ) is deleted by Snakemake as soon as `<tool>_record`
  succeeds, and also when the run fails. Do not pass `--notemp`.
- If `<tool>_record` fails, the scratch directory is kept so the time files
  can be re-parsed without re-running the assembly.
- `<tool>_run` refuses to start if Snakemake grants fewer threads than
  `ram_time_threads`, so always pass `--cores 64` (or more).
- Nothing is ever written to `assemblies/`.

## Configuration

No paths are hard-coded. Everything is resolved relative to the directory
you launch Snakemake from (the repository root). Optional `--config` keys:

| Key | Default | Purpose |
|---|---|---|
| `ram_time_threads` | `64` | Threads given to the assembler and to `docker --cpus` |
| `ram_time_dir` | `assemblers_ram_time` | Where `.ram_time.tsv` files and rule logs are written |
| `ram_time_scratch_dir` | `assemblers_ram_time/scratch` | Temporary assembly workspace, deleted after each run |

Both directories must stay **inside** the repository root
(CONSTITUTION I.1); the workflow refuses to start otherwise. No
`--configfile` is needed.

Example with a different scratch directory:

```bash
snakemake --snakefile assemblers/whole_genome_asm/ram_time/ont.assembly.flye2.ram_time.smk --cores 64 --resources mem_gb=240 --config ram_time_scratch_dir=assemblers_ram_time/scratch_run2
```

## Inputs and outputs

Inputs are the production FASTQs, discovered with the same filters as the
production workflows (`.1k`, `.chr21.`, `localtest` and `smoke` datasets are
skipped):

```text
fastq/{sample}.{ont,pb}.30x.fastq.gz
```

Verkko needs both the ONT and the PacBio HiFi 30x FASTQ for every sample and
stops with an error if one is missing.

Outputs, one small TSV per dataset plus its logs:

```text
assemblers_ram_time/{flye,goldrush,verkko}/{dataset}.ram_time.tsv
assemblers_ram_time/{flye,goldrush,verkko}/{dataset}.{run,record}.log
```

Columns: `assembler`, `sample`, `technology`, `threads`,
`wall_clock_seconds`, `wall_clock_hours`, `peak_rss_gb`,
`n_stages_summed` (1 for Flye/Verkko, >1 for GoldRush), `source_time_files`.

## Where to run

Launch from the **repository root** on the host that has the `fastq/`
inputs and the Docker daemon used for the production runs. The repository
path differs between machines; it has been run from, for example:

```text
/data/genmedbfx/yu_j/lrs_benchmarking      # Flye + reference + QUAST
/home/stoiber_l/smbshare/lrs_benchmarking  # GoldRush / Verkko outputs
```

Check before running:

```bash
cd /path/to/lrs_benchmarking    # replace with the real path on this host
ls -lh fastq/*.30x.fastq.gz      # inputs are visible
docker info >/dev/null && echo "Docker daemon reachable"
```

### Verify which FASTQs each workflow will use

The FASTQs must be in `fastq/` directly under the repository root (the
directory you launch from). Flye and GoldRush take **every**
`fastq/*.fastq.gz` whose name contains `.ont.` (or `.pb.`), not only `30x`,
so any extra file there (e.g. `HG003.ont.60x.fastq.gz`) becomes an extra
multi-hour job. List each input with the real file it points to (symlinks
are followed) and its size:

```bash
for f in fastq/*.fastq.gz; do case "$f" in *.1k*|*.chr21.*|*localtest*|*smoke*|*SMOKE*) continue;; esac; printf '%-35s %6s  %s\n' "$f" "$(du -hL "$f" | cut -f1)" "$(readlink -f "$f")"; done
```

A whole-genome 30x FASTQ is tens of GB. A small file, or a target path
containing `chr21`, `1k` or `localtest`, means the symlink points to a test
subset, and the run would measure that subset instead of the whole genome.

Then do a dry run of each workflow. At startup, each workflow prints the `.ram_time.tsv` targets it found, one per
input dataset, followed by the job count:

```bash
for smk in assemblers/whole_genome_asm/ram_time/*.ram_time.smk; do echo "=== $smk"; snakemake --snakefile "$smk" --cores 64 --dry-run --quiet rules 2>&1 | grep -E 'ram_time.tsv|_run|_record|Error'; done
```

Expected output, one target per sample (shown for HG002):

```text
=== .../hybrid.assembly.verkko.ram_time.smk
['.../assemblers_ram_time/verkko/HG002.ram_time.tsv']
verkko_record        1
verkko_run           1
=== .../ont.assembly.flye2.ram_time.smk
['.../assemblers_ram_time/flye/HG002.ont.30x.ram_time.tsv']
...
```

If a list is empty (`[]`), the workflow found no input: you are not in the
repository root, or `fastq/` is missing or named differently. If it lists
datasets you did not intend to run, move those FASTQs out of `fastq/`
before the real run.

## Commands

Always do a dry-run first (add `--dry-run`), then the real run.

### Flye 2.9.6 (`nicolasardila1/lrs-flye2:2.9.6`)

```bash
# ONT
snakemake --snakefile assemblers/whole_genome_asm/ram_time/ont.assembly.flye2.ram_time.smk --cores 64 --resources mem_gb=240 --rerun-incomplete --printshellcmds --show-failed-logs

# PacBio HiFi
snakemake --snakefile assemblers/whole_genome_asm/ram_time/pb.assembly.flye2.ram_time.smk --cores 64 --resources mem_gb=240 --rerun-incomplete --printshellcmds --show-failed-logs
```

### GoldRush 1.2.2-ntlinkfix (`nicolasardila1/lrs-goldrush:1.2.2-ntlinkfix`)

```bash
# ONT
snakemake --snakefile assemblers/whole_genome_asm/ram_time/ont.assembly.goldrush.ram_time.smk --cores 64 --resources mem_mb=64000 --rerun-incomplete --printshellcmds --show-failed-logs

# PacBio HiFi
snakemake --snakefile assemblers/whole_genome_asm/ram_time/pb.assembly.goldrush.ram_time.smk --cores 64 --resources mem_mb=64000 --rerun-incomplete --printshellcmds --show-failed-logs
```

### Verkko 2.3.2 (`nicolasardila1/lrs-verkko2:2.3.2`)

```bash
snakemake --snakefile assemblers/whole_genome_asm/ram_time/hybrid.assembly.verkko.ram_time.smk --cores 64 --resources mem_mb=200000 --rerun-incomplete --printshellcmds --show-failed-logs
```

## Combining the results

After all runs finish (15 `.ram_time.tsv` files across the three assemblers,
all samples and technologies), merge them into one table:

```bash
{
    head -n 1 "$(ls assemblers_ram_time/*/*.ram_time.tsv | head -n 1)"
    tail -n +2 -q assemblers_ram_time/*/*.ram_time.tsv
} > assembly_analysis/tables/assembler_performance_30x.tsv
```

Thread-hours are `threads x wall_clock_hours`, both already in the table.

## Validation before commit

From the repository root, for every modified workflow:

```bash
git diff --check
snakemake --snakefile assemblers/whole_genome_asm/ram_time/<workflow>.ram_time.smk --cores 64 --dry-run --printshellcmds
```
