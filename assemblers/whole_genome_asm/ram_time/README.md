# Assembler RAM / Time Measurement

Companion to [`../README.md`](../README.md). This folder only measures the
**computational cost** of the whole-genome assemblers: wall-clock time, peak
RAM, CPU time and peak scratch disk use. It never produces assemblies used in
the report.

```text
ont.assembly.flye2.ram_time.smk
pb.assembly.flye2.ram_time.smk
ont.assembly.goldrush.ram_time.smk
pb.assembly.goldrush.ram_time.smk
hybrid.assembly.verkko.ram_time.smk
sample_container_resources.sh   # background memory/disk sampler used by all five
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

`time -v` gives wall-clock time, CPU time, exit status and the peak RAM of the
**largest single process**. That last value is wrong for a tool that runs
many processes at once, so every run also starts
`sample_container_resources.sh` in the background:

- Every **30 s** it reads the memory of the **whole container** (all
  processes together, page cache excluded) with `docker stats`.
- Every **5 min** it measures the size of the scratch directory with `du`,
  run *inside* the assembler container through `docker exec` (CONSTITUTION
  II.1), plus once more after the assembler finishes.
- It waits for the named container (`ramtime-<tool>-<dataset>`) to start and
  stops when it exits; the run rule also stops it if anything fails.

Before the assembly starts, the run step also checks that `du` is in the
image and that `docker stats` works on the host, and fails within seconds
otherwise. After the run it fails if the sampler recorded no samples.

The same checks and the sampler apply to all three tools. Only what is
timed differs (one `time -v` file for Flye/Verkko, one per stage for
GoldRush).

## Rules: run -> record

| Rule | What it does | Kept? |
|---|---|---|
| `<tool>_run` | Runs the assembler inside a scratch directory declared as `temp(directory(...))` | No |
| `<tool>_record` | Parses the time file(s) and the resource samples into `{dataset}.ram_time.tsv` with `assembly_analysis/scripts/metrics/parse_assembler_time_v.py`, run with `python3` inside the assembler's own image (CONSTITUTION II.1), and copies the raw files to `{dataset}.raw/` | **Yes** |

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

Outputs, one small TSV per dataset, the raw files it was computed from, and
the logs:

```text
assemblers_ram_time/{flye,goldrush,verkko}/{dataset}.ram_time.tsv
assemblers_ram_time/{flye,goldrush,verkko}/{dataset}.raw/   # time -v file(s) + resource_samples.tsv
assemblers_ram_time/{flye,goldrush,verkko}/{dataset}.{run,record}.log
```

`{dataset}.raw/` is a few KB. It holds the files the table is computed from,
so every number in `{dataset}.ram_time.tsv` can be checked later by hand.

Columns: `assembler`, `sample`, `technology`, `threads`,
`wall_clock_seconds`, `wall_clock_hours`, `peak_rss_gb`, `cpu_hours`,
`cpu_efficiency`, `exit_status`, `container_peak_mem_gb`,
`peak_scratch_disk_gb`, `n_mem_samples`, `n_stages_summed` (1 for
Flye/Verkko, >1 for GoldRush), `source_time_files`.

- `peak_rss_gb`: peak RAM of the largest single process (`time -v`), max
  over GoldRush stages.

- `cpu_hours`: user + system CPU time, summed over GoldRush stages. Unlike
  wall-clock time it barely depends on the thread count, so it stays
  comparable with the 32-thread production runs.
- `cpu_efficiency`: CPU time / (wall-clock time x `threads`). 1.0 means all
  threads were busy for the whole run; long single-threaded stages lower it.
- `exit_status`: always 0. If any time file reports a non-zero exit status,
  the record step fails and no row is written, so a failed run is never
  reported.
- `container_peak_mem_gb`: highest sampled memory of the whole container.
- `peak_scratch_disk_gb`: largest sampled size of the scratch directory
  (intermediates, and for GoldRush the decompressed FASTQ).
- `n_mem_samples`: number of 30 s memory samples behind
  `container_peak_mem_gb` (about 120 per hour of run time).

Which peak RAM to report:

| Tool | Use | Why |
|---|---|---|
| Flye | `peak_rss_gb` | Essentially one multi-threaded process; `time -v` sees every byte and misses no spike |
| GoldRush | `peak_rss_gb` | Stages run one after another; the max over stages is the peak |
| Verkko | `container_peak_mem_gb` | Runs many jobs in parallel; `peak_rss_gb` sees only the largest one and undercounts |

For Flye and GoldRush the two values should be close; a large gap is worth
a look in `{dataset}.raw/resource_samples.tsv`. Sampled values can miss
spikes shorter than 30 s (memory) or 5 min (disk), so they are slightly
below the true peak.

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
docker stats --no-stream >/dev/null && echo "docker stats works (needed by the sampler)"
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

### If a run is interrupted

After an interruption (e.g. Ctrl-C), check that its container is gone
before re-running; a leftover one makes the next run fail with "name is
already in use":

```bash
docker ps -a --filter name=ramtime- --format '{{.Names}} {{.Status}}'
docker rm -f ramtime-<tool>-<dataset>    # only if listed and not wanted
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
