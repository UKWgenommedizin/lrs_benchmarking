##
# header_mapper.smk
# Minimal header for read mapping workflows.
# Provides CWD only — no VarCAD dependency.
# Constitution: Article I.6

import os

CWD = os.getcwd()
print("Current working directory: " + CWD)