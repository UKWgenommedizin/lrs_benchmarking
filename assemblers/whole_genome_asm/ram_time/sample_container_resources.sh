#!/usr/bin/env bash
# ************************************************************************************************
# Background resource sampler for the *.ram_time.smk workflows.
#
# Usage: sample_container_resources.sh CONTAINER_NAME SCRATCH_DIR OUT_TSV [INTERVAL_S] [DISK_EVERY_N]
#
# Every INTERVAL_S seconds (default 30) records the memory used by the whole
# container (`docker stats`, all processes together, page cache excluded), and
# every DISK_EVERY_N samples (default 10, i.e. 5 min) the size of SCRATCH_DIR,
# measured with `du` *inside* the running container via `docker exec`
# (CONSTITUTION II.1: no host-installed tools besides docker).
#
# Why: /usr/bin/time -v only reports the largest single process. That is
# fine for Flye and the sequential GoldRush stages, but Verkko runs many jobs
# in parallel inside one container, so only the container total is its real
# peak RAM.
#
# Output TSV: epoch_seconds, container_mem_bytes, scratch_disk_bytes (NA when
# not measured in that sample). Waits for the container to start and exits
# by itself once it has stopped; the calling rule also kills it on exit.
# ************************************************************************************************

set -uo pipefail

NAME="$1"
SCRATCH="$2"
OUT="$3"
INTERVAL="${4:-30}"
DISK_EVERY="${5:-10}"

printf 'epoch_seconds\tcontainer_mem_bytes\tscratch_disk_bytes\n' > "$OUT"

# Convert a docker size string ("512KiB", "1.5GiB", "800MB", "0B") to bytes.
to_bytes() {
    awk -v v="$1" 'BEGIN {
        if (match(v, /^[0-9.]+/) == 0) { print "NA"; exit }
        n = substr(v, 1, RLENGTH); u = substr(v, RLENGTH + 1)
        f["B"] = 1
        f["KiB"] = 1024; f["MiB"] = 1024^2; f["GiB"] = 1024^3; f["TiB"] = 1024^4
        f["kB"] = 1e3;   f["KB"] = 1e3;     f["MB"] = 1e6;      f["GB"] = 1e9; f["TB"] = 1e12
        if (!(u in f)) { print "NA"; exit }
        printf "%.0f\n", n * f[u]
    }'
}

measure_disk() {
    local disk
    disk=$(docker exec "$NAME" du -sb "$SCRATCH" 2>/dev/null | cut -f1)
    [[ "$disk" =~ ^[0-9]+$ ]] || disk="NA"
    printf '%s\tNA\t%s\n' "$(date +%s)" "$disk" >> "$OUT"
}

trap 'exit 0' TERM INT

seen=0
misses=0
i=0
disk_pid=""
while true; do
    raw=$(docker stats --no-stream --format '{{.MemUsage}}' "$NAME" 2>/dev/null)
    if [[ -n "$raw" && "$raw" != --* ]]; then
        seen=1
        misses=0
        printf '%s\t%s\tNA\n' "$(date +%s)" "$(to_bytes "${raw%% / *}")" >> "$OUT"
        # du runs in the background so a slow scan of a large scratch
        # directory never pauses the memory sampling; a new scan starts
        # only after the previous one has finished.
        if (( i % DISK_EVERY == 0 )) && { [[ -z "$disk_pid" ]] || ! kill -0 "$disk_pid" 2>/dev/null; }; then
            measure_disk &
            disk_pid=$!
        fi
        i=$(( i + 1 ))
    elif (( seen )); then
        # Three misses in a row (about 1.5 min at 30 s): the container has
        # stopped. A single failed call can be a Docker hiccup, not the end.
        misses=$(( misses + 1 ))
        (( misses >= 3 )) && exit 0
    fi
    # Background sleep + wait so a TERM from the rule is handled immediately.
    sleep "$INTERVAL" &
    wait $!
done
