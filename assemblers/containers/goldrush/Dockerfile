# =============================================================================
# GOLDRUSH 1.2.2 + GNU GZIP COMPATIBILITY FIX
# =============================================================================
#
# Problem:
#
# GoldRush bundles ntLink, which invokes:
#
#     gzip -f -cd <reads>
#
# Some read files supplied internally to ntLink are already uncompressed
# FASTQ files.
#
# The upstream BioContainers image provides BusyBox gzip. BusyBox returns
# "invalid magic" for plain FASTQ input instead of passing it through.
#
# GNU gzip provides the passthrough behavior expected by this pipeline.
#
# Solution:
#
#   1. Obtain GNU gzip from Debian 12 in a temporary build stage.
#   2. Copy only the GNU gzip executable into /usr/local/bin.
#   3. Keep the original BusyBox installation untouched.
#
# /usr/local/bin appears before /usr/bin in PATH, therefore the GNU gzip
# executable is selected by GoldRush/ntLink.
# =============================================================================


# -----------------------------------------------------------------------------
# Global build arguments
# -----------------------------------------------------------------------------
#
# IMPORTANT:
# An ARG used by a FROM instruction must be declared before the first FROM.
# -----------------------------------------------------------------------------

ARG BASE_IMAGE=quay.io/biocontainers/goldrush:1.2.2--py312h248b54c_1


# -----------------------------------------------------------------------------
# Stage 1 — obtain GNU gzip
# -----------------------------------------------------------------------------

FROM debian:12-slim AS gzip_source

RUN apt-get update \
    && apt-get install -y --no-install-recommends gzip \
    && rm -rf /var/lib/apt/lists/*


# -----------------------------------------------------------------------------
# Stage 2 — original GoldRush BioContainer
# -----------------------------------------------------------------------------

FROM ${BASE_IMAGE}

LABEL org.opencontainers.image.title="GoldRush assembler"
LABEL org.opencontainers.image.description="GoldRush 1.2.2 with GNU gzip compatibility fix"
LABEL org.opencontainers.image.version="1.2.2-gzipfix1"
LABEL org.opencontainers.image.vendor="Nicolas Ardila F2 benchmarking workflow"
LABEL org.opencontainers.image.base.name="quay.io/biocontainers/goldrush:1.2.2--py312h248b54c_1"


# -----------------------------------------------------------------------------
# GNU gzip compatibility fix
# -----------------------------------------------------------------------------

COPY --from=gzip_source /usr/bin/gzip /usr/local/bin/gzip


# -----------------------------------------------------------------------------
# Build-time regression tests
# -----------------------------------------------------------------------------
#
# These checks ensure that:
#
#   1. GoldRush is still installed.
#   2. /usr/local/bin/gzip is selected instead of BusyBox gzip.
#   3. GNU gzip can pass an uncompressed FASTQ through `gzip -f -cd`.
#
# The third test reproduces the exact operation that failed on the server.
# -----------------------------------------------------------------------------

RUN set -eux; \
    command -v goldrush; \
    test "$(command -v gzip)" = "/usr/local/bin/gzip"; \
    gzip --version; \
    printf '@test_read\nACGT\n+\nIIII\n' > /tmp/plain.fastq; \
    gzip -f -cd /tmp/plain.fastq > /tmp/passthrough.fastq; \
    cmp /tmp/plain.fastq /tmp/passthrough.fastq; \
    rm -f /tmp/plain.fastq /tmp/passthrough.fastq


WORKDIR /work

ENTRYPOINT []

CMD ["/bin/bash"]
