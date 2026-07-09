##
# header_mapper.smk
# Minimal header for read mapping workflows.
# Provides CWD only — no VarCAD dependency.
# Constitution: Article I.6

import os
from snakemake.utils import min_version

CWD = os.getcwd()
print("Current working directory: " + CWD)

# Optional dataset filter for testing (e.g., --config dataset_filter=1k)
try:
    DATASET_FILTER = config["dataset_filter"]
except (KeyError, NameError):
    DATASET_FILTER = None
