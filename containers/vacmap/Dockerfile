# =============================================================================
# VACmap 1.2.0 + RAM/TIME MEASUREMENT TOOLS
# =============================================================================
#
# Purpose:
#
# The ram_time workflows (alignment_analysis/run_metrics/ram_time/) need the
# same tools inside the VACmap image that the VG Giraffe image provides:
#
#   /usr/bin/time   GNU time, `time -v` is the container entrypoint per stage
#   python3         runs parse_assembler_time_v.py (stdlib only) in vacmap_record
#   du -sb          GNU coreutils du, used for the scratch-disk peak and by
#                   sample_container_resources.sh
#
# Solution:
#
#   Start from the exact production image (schimar/lrs-vacmap:v1.2.0), so the
#   VACmap and samtools binaries are identical to the production mapping
#   runs, and only install whichever of the three tools is missing.
# =============================================================================

ARG BASE_IMAGE=schimar/lrs-vacmap:v1.2.0

FROM ${BASE_IMAGE}

LABEL org.opencontainers.image.title="VACmap aligner (ram_time)"
LABEL org.opencontainers.image.description="VACmap 1.2.0 production image plus GNU time, python3 and GNU du for RAM/time measurement"
LABEL org.opencontainers.image.version="1.2.0"
LABEL org.opencontainers.image.vendor="Nicolas Ardila F2 benchmarking workflow"
LABEL org.opencontainers.image.base.name="schimar/lrs-vacmap:v1.2.0"


# -----------------------------------------------------------------------------
# Install only the missing tools (base distro is detected, not assumed)
# -----------------------------------------------------------------------------

USER root
RUN set -eux; \
    pkgs=""; \
    [ -x /usr/bin/time ] || pkgs="$pkgs time"; \
    command -v python3 >/dev/null 2>&1 || pkgs="$pkgs python3"; \
    du -sb /tmp >/dev/null 2>&1 || pkgs="$pkgs coreutils"; \
    if [ -z "$pkgs" ]; then \
        echo "All ram_time tools already present"; \
    elif command -v apt-get >/dev/null 2>&1; then \
        apt-get update; \
        DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends $pkgs; \
        rm -rf /var/lib/apt/lists/*; \
    elif command -v apk >/dev/null 2>&1; then \
        apk add --no-cache $pkgs; \
    elif command -v dnf >/dev/null 2>&1; then \
        dnf install -y $pkgs && dnf clean all; \
    elif command -v yum >/dev/null 2>&1; then \
        yum install -y $pkgs && yum clean all; \
    else \
        echo "ERROR: no supported package manager to install:$pkgs"; exit 1; \
    fi


# -----------------------------------------------------------------------------
# Build-time checks
# -----------------------------------------------------------------------------
#
# 1. The exact check the ram_time run step does before mapping.
# 2. The aligner and samtools from the production image are still there.
# -----------------------------------------------------------------------------

RUN set -eux; \
    command -v /usr/bin/time && command -v python3 && du -sb /tmp >/dev/null; \
    /usr/bin/time -v true; \
    python3 --version; \
    command -v vacmap; \
    samtools --version | head -n 1
