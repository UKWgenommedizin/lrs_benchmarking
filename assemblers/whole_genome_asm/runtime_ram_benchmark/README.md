# Assembler Performance Benchmark (Wall-Clock Time / Peak RAM)

Companion to [`../README.md`](../README.md). This folder holds only the
**measurement-only** workflows; they never produce assemblies used in the
report. The production workflows in the parent directory (`ont.assembly.flye2.smk`, `pb.assembly.flye2.smk`,
`ont.assembly.goldrush.smk`, `pb.assembly.goldrush.smk`,
`hybrid.assembly.verkko.smk`) were never instrumented to record wall-clock
time or peak RAM for the actual 30x whole-genome runs -- there is no
`benchmark:` directive in any of them, and no full-run log survives on
disk to recover the numbers after the fact (see
`assembly_analysis/scripts/utils/benchmark_data.py`'s module docstring for
the fuller investigation).

This is why these five files exist:

```text
ont.assembly.flye2.benchmark.smk
pb.assembly.flye2.benchmark.smk
ont.assembly.goldrush.benchmark.smk
pb.assembly.goldrush.benchmark.smk
hybrid.assembly.verkko.benchmark.smk
```

Each is a **runtime / peak-RAM measurement-only variant** of its
production counterpart: same Docker image/version, same memory ceiling,
same `goldrush run` / `flye` / `verkko` parameters. They run at **64
threads** (`benchmark_threads`, configurable), not the production 32, so
these numbers describe a 64-thread run. The only result kept is the
time/RAM table; the regenerated assembly is always deleted. See
each file's own header comment for the reasoning specific to it.

## Why not just re-run the production workflows with `benchmark:` added?

Snakemake's `benchmark:` directive samples the *host-side* process tree.
Every assembler here runs inside `docker run`, and Docker's real memory
usage lives in a separate cgroup that a host-side sampler cannot see --
`benchmark:` would report a near-zero, meaningless RAM number for a
containerized process. Instead:

- **Flye and Verkko**: the actual tool binary is wrapped in
  `/usr/bin/time -v` *inside* the container (`--entrypoint /usr/bin/time`).
  Each rule preflights that `/usr/bin/time` exists in the image before
  starting the multi-hour assembly, and fails fast with a clear error if
  not.
- **GoldRush**: already self-instruments every internal Makefile stage
  (silver_path, golden_path, goldpolish, tigmint, 5x ntLink rounds, ...)
  via `command time -v -o <stage>.time <command>` whenever `track_time=1`
  is set -- which the production rule already passes. The benchmark
  variant collects every `*.time` file after a successful run and
  aggregates them: **sums** each stage's elapsed wall-clock time (there is
  no single end-to-end number, only per-stage ones) and takes the **max**
  of each stage's peak RSS (memory is not cumulative across sequential
  stages the way time is).

## Rules: run -> record (scratch deleted automatically)

Each workflow has two rules per dataset:

| Rule | What it does | Kept? |
|---|---|---|
| `<tool>_run` | Runs the assembler under `/usr/bin/time -v` (GoldRush: `track_time=1`) inside a scratch directory declared as `temp(directory(...))` | No |
| `<tool>_record` | Parses the time file(s) into `{dataset}.perf.tsv` (parser runs with `python3` inside the assembler's own image, CONSTITUTION II.1) | **Yes** |

Because the scratch directory is `temp()`, Snakemake deletes it
(assembly, intermediates, raw time files, GoldRush's decompressed FASTQ)
as soon as `<tool>_record` succeeds, and also removes it if the run
fails. A later `snakemake` call does not regenerate it once `.perf.tsv`
exists. If `<tool>_record` fails, the scratch directory is kept so the
time files can be re-parsed without re-running the assembly. Do not pass
`--notemp`, which disables this deletion. QUAST metrics for these
assemblies already exist, so nothing in the scratch tree is needed after
the measurement. **This never touches `assemblies/`.**

`<tool>_run` refuses to start if Snakemake granted fewer threads than
`benchmark_threads` (e.g. `--cores 64` with the default of 64), so a
measurement is never silently taken at a different thread count.

## Configuration (no hard-coded paths)

All locations and the thread count come from `--config` (or a
`--configfile`); relative paths are resolved against the repository root.
Both directories must stay **inside** the repository root
(CONSTITUTION I.1) -- the workflow refuses to start otherwise:

| Key | Default | Purpose |
|---|---|---|
| `benchmark_threads` | `64` | Threads given to the assembler and to `docker --cpus` |
| `benchmark_dir` | `assemblers_benchmark` | Where `.perf.tsv` and the per-rule logs are written |
| `benchmark_scratch_dir` | `assemblers_benchmark/scratch` | Temporary assembly workspace, deleted after each run |

Example:

```bash
snakemake --snakefile assemblers/whole_genome_asm/runtime_ram_benchmark/ont.assembly.flye2.benchmark.smk \
    --cores 64 --resources mem_gb=240 \
    --config benchmark_scratch_dir=assemblers_benchmark/scratch_run2
```

## Input / output convention

Inputs are the same production FASTQs, discovered the same way as the
production workflows:

```text
fastq/
```

Outputs are one small TSV per dataset (plus its rule logs), never the
assembly itself:

```text
{benchmark_dir}/{assembler}/{dataset}.perf.tsv
{benchmark_dir}/{assembler}/{dataset}.{run,record}.log
```

Columns: `assembler`, `sample`, `technology`, `threads`,
`wall_clock_seconds`, `wall_clock_hours`, `peak_rss_gb`,
`n_stages_summed` (1 for Flye/Verkko, >1 for GoldRush),
`source_time_files`.

Example:

```text
assemblers_benchmark/flye/HG002.ont.30x.perf.tsv
assemblers_benchmark/goldrush/HG002.pb.30x.perf.tsv
assemblers_benchmark/verkko/HG002.perf.tsv
```

`assemblers_benchmark/` is kept separate from `assemblies/` and never
mixed into it. Only the `.perf.tsv` files carry results.

## Run from the repository root

Like the production workflows, these files use repository-root-relative
paths (`CWD = os.getcwd()`), so
they must be launched with the repository root as the working directory.

**Adjust the path below to wherever this repository actually lives on
the execution host** -- per `../assessment/README.md`, that is not
necessarily the same path as on this machine. The production workflows
have been run from at least these locations:

```text
/data/genmedbfx/yu_j/lrs_benchmarking      # Flye + reference + QUAST
/home/stoiber_l/smbshare/lrs_benchmarking  # GoldRush / Verkko outputs
```

Run Snakemake from the host that owns the Docker daemon actually used
for the production runs -- not from inside another container, and not
from a machine that only has this repository checked out without the
`fastq/` inputs or Docker access. Verify before running:

```bash
cd /path/to/lrs_benchmarking   # <- replace with the real path on this host
pwd                             # confirm it matches
ls -lh fastq/*.30x.fastq.gz     # confirm inputs are visible
docker info >/dev/null && echo "Docker daemon reachable"
```

## Flye (ONT / PacBio HiFi)

Same version as the production workflow:

```text
Flye 2.9.6
Docker: nicolasardila1/lrs-flye2:2.9.6
```

Dry-run ONT:

```bash
snakemake --snakefile assemblers/whole_genome_asm/runtime_ram_benchmark/ont.assembly.flye2.benchmark.smk --cores 64 --resources mem_gb=240 --dry-run --printshellcmds
```

Run ONT:

```bash
snakemake --snakefile assemblers/whole_genome_asm/runtime_ram_benchmark/ont.assembly.flye2.benchmark.smk --cores 64 --resources mem_gb=240 --rerun-incomplete --printshellcmds --show-failed-logs
```

Dry-run / run PacBio HiFi by replacing the Snakefile with:

```text
assemblers/whole_genome_asm/runtime_ram_benchmark/pb.assembly.flye2.benchmark.smk
```

## GoldRush (ONT / PacBio HiFi)

Same version as the production workflow:

```text
GoldRush 1.2.2-ntlinkfix
Docker: nicolasardila1/lrs-goldrush:1.2.2-ntlinkfix
```

Dry-run ONT:

```bash
snakemake --snakefile assemblers/whole_genome_asm/runtime_ram_benchmark/ont.assembly.goldrush.benchmark.smk --cores 64 --resources mem_mb=64000 --dry-run --printshellcmds
```

Run ONT:

```bash
snakemake --snakefile assemblers/whole_genome_asm/runtime_ram_benchmark/ont.assembly.goldrush.benchmark.smk --cores 64 --resources mem_mb=64000 --rerun-incomplete --printshellcmds --show-failed-logs
```

Dry-run / run PacBio HiFi by replacing the Snakefile with:

```text
assemblers/whole_genome_asm/runtime_ram_benchmark/pb.assembly.goldrush.benchmark.smk
```

If a run fails partway (as one already has -- see the ntLink pairing
error discussed for `HG002.ont.30x` while designing this benchmark), the
rule will refuse to trust the partial `*.time` files as a total and exit
with an error before writing a `.perf.tsv` row, rather than silently
reporting an incomplete run's time as if it were the full pipeline's.

## Verkko

Same version as the production workflow:

```text
Verkko 2.3.2
Docker: nicolasardila1/lrs-verkko2:2.3.2
```

Dry-run:

```bash
snakemake --snakefile assemblers/whole_genome_asm/runtime_ram_benchmark/hybrid.assembly.verkko.benchmark.smk --configfile assemblers/config/server.yaml --cores 64 --resources mem_mb=200000 --dry-run --printshellcmds
```

Run:

```bash
snakemake --snakefile assemblers/whole_genome_asm/runtime_ram_benchmark/hybrid.assembly.verkko.benchmark.smk --configfile assemblers/config/server.yaml --cores 64 --resources mem_mb=200000 --rerun-incomplete --printshellcmds --show-failed-logs
```

Verkko requires matching ONT and PacBio HiFi 30x FASTQs for each sample,
same as the production workflow.

## Combining the per-dataset tables

Each run produces one `.perf.tsv` per dataset (15 total across all three
assemblers, all samples and technologies). To combine them into one
table once all runs are complete:

```bash
{
    head -n 1 "$(ls assemblers_benchmark/*/*.perf.tsv | head -n 1)"
    tail -n +2 -q assemblers_benchmark/*/*.perf.tsv
} > assembly_analysis/tables/30x/final/assembler_performance_30x.tsv
```

Thread-hours (threads x wall_clock_hours) is not computed by
`parse_assembler_time_v.py` -- both columns it needs
(`threads`, `wall_clock_hours`) are already in the combined table, so it
is a one-line calculation from there rather than something this
benchmark tooling needs to produce.

## Validation before commit

For every modified workflow, from the repository root:

```bash
git diff --check
snakemake --snakefile assemblers/whole_genome_asm/runtime_ram_benchmark/<workflow>.benchmark.smk --dry-run --printshellcmds
```

These workflows resolve every path from the repository root
(`CWD = os.getcwd()`) and have no relative `include:`, so they live in this
separate `runtime_ram_benchmark/` folder to keep them apart from the
production assembly workflows. Always launch them from the repository root.
