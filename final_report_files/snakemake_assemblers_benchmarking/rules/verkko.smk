#####################################################################################
# VERKKO HYBRID ASSEMBLY RULES
# Verkko combines the ONT and PacBio HiFi datasets belonging to the same sample
# to perform hybrid diploid genome assembly.
#####################################################################################


# -----------------------------------------------------------------------------
# Future Rule 5: Verkko hybrid assembly
# -----------------------------------------------------------------------------

# Planned inputs:
#   - Input validation marker
#   - PacBio HiFi FASTQ belonging to one sample
#   - ONT FASTQ belonging to the same sample
#
# Planned outputs:
#   - Final Verkko assembly graph
#   - Primary or haplotype assembly FASTA files
#   - Assembly logs and intermediate workflow files
#
# This rule is intentionally inactive until the installation, supported ONT
# input type, output structure, and required computing resources are validated.
