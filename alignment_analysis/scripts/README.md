# Alignment Analysis Scripts

This directory contains scripts supporting the long-read alignment benchmark.

## Canonical production script

```text
30x/quality_check_aligners.py
```

This is the canonical metric-extraction script for the 30x alignment benchmark.

## Recommended organization

```text
scripts/
├── 30x/           # canonical production extraction
├── plots/         # figure-generation scripts
├── validation/    # targeted test / validation runners
├── utils/         # FASTQ and small helper utilities
└── legacy/        # retained historical scripts, not canonical
```

The migration helper supplied with the repository-organization kit can move selected analysis scripts while leaving compatibility entry points at their old locations.

Do not move `30x/quality_check_aligners.py` unless its project-root detection and every documented command are updated and validated.
