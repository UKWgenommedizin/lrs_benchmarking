# =============================================================================
# GOLDRUSH 1.2.2 + ntLink PLAIN-FASTQ COMPATIBILITY FIX
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
LABEL org.opencontainers.image.description="GoldRush 1.2.2 with ntLink plain-FASTQ compatibility fix"
LABEL org.opencontainers.image.version="1.2.2-ntlinkfix"
LABEL org.opencontainers.image.vendor="Nicolas Ardila F2 benchmarking workflow"
LABEL org.opencontainers.image.base.name="quay.io/biocontainers/goldrush:1.2.2--py312h248b54c_1"


# -----------------------------------------------------------------------------
# GNU gzip compatibility fix
# -----------------------------------------------------------------------------

COPY --from=gzip_source /usr/bin/gzip /usr/local/bin/gzip

# The workflow deliberately gives ntLink an uncompressed FASTQ. Patch the
# bundled ntLink command so it reads that file directly and ignores empty
# checkpoint mappings.
USER root
RUN set -eux; \
    python3 - <<'PY'
from pathlib import Path

path = Path("/usr/local/bin/share/ntlink-1.3.11-1/ntLink")
text = path.read_text()
old = "$(gzip) -cd $(reads)"
count = text.count(old)
assert count > 0, "Expected ntLink gzip command was not found"
path.write_text(text.replace(old, "cat $(reads)"))
print("Patched ntLink plain-FASTQ handling:", count, "occurrence(s)")
PY

RUN set -eux; \
    python3 - <<'PY'
from pathlib import Path

path = Path("/usr/local/bin/share/ntlink-1.3.11-1/bin/ntlink_pair.py")
text = path.read_text()
old = """        for i, m in enumerate(mappings):
            contig_runs.append(m.contig_id)
            accepted_anchor_contigs[m.contig_id] = ContigRun(m.contig_id, ntlink_utils.parse_minimizers(m.list_hits))
"""
new = """        for i, m in enumerate(mappings):
            if not m.list_hits.strip():
                continue
            hits = ntlink_utils.parse_minimizers(m.list_hits)
            if not hits:
                continue
            contig_runs.append(m.contig_id)
            accepted_anchor_contigs[m.contig_id] = ContigRun(m.contig_id, hits)
"""
assert text.count(old) == 1, "Expected ntLink checkpoint parser block was not found"
text = text.replace(old, new)
old = """        length_read = max(read_mapping_positions)
        self.tally_pairs_from_mappings(accepted_anchor_contigs, contig_runs, length_read, pairs)
"""
new = """        if not read_mapping_positions:
            return
        length_read = max(read_mapping_positions)
        self.tally_pairs_from_mappings(accepted_anchor_contigs, contig_runs, length_read, pairs)
"""
assert text.count(old) == 1, "Expected ntLink empty-mapping guard location was not found"
path.write_text(text.replace(old, new))
print("Patched ntLink empty checkpoint mapping handling")
PY


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
    grep -qF 'cat $(reads)' /usr/local/bin/share/ntlink-1.3.11-1/ntLink; \
    ! grep -qF '$(gzip) -cd $(reads)' /usr/local/bin/share/ntlink-1.3.11-1/ntLink; \
    grep -qF 'if not m.list_hits.strip():' /usr/local/bin/share/ntlink-1.3.11-1/bin/ntlink_pair.py; \
    grep -qF 'if not hits:' /usr/local/bin/share/ntlink-1.3.11-1/bin/ntlink_pair.py; \
    grep -qF 'if not read_mapping_positions:' /usr/local/bin/share/ntlink-1.3.11-1/bin/ntlink_pair.py; \
    printf '@test_read\nACGT\n+\nIIII\n' > /tmp/plain.fastq; \
    gzip -f -cd /tmp/plain.fastq > /tmp/passthrough.fastq; \
    cmp /tmp/plain.fastq /tmp/passthrough.fastq; \
    rm -f /tmp/plain.fastq /tmp/passthrough.fastq


WORKDIR /work

ENTRYPOINT []

CMD ["/bin/bash"]
