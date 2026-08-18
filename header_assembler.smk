##
# header_assembler.smk
# Minimal shared header for WGS assembly workflows.
##

import os

############################
# Current working directory

CWD = os.getcwd()
print("Current working directory: " + CWD)


############################
# Optional dataset filter

try:
    DATASET_FILTER = config["dataset_filter"]
except (KeyError, NameError):
    DATASET_FILTER = None
