##
# ont.read_mapping.ntlink.smk
# mapper_tag: ntlink-ont
#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# BLOCKER — HUMAN REVIEW REQUIRED BEFORE THIS FILE CAN BE IMPLEMENTED
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
# ntLink (github.com/bcgsc/ntLink) is a genome assembly *scaffolder*, not a
# read-to-reference mapper.
#
# ntLink maps long reads to *draft assembly contigs* (in order to scaffold
# them), not to a fixed linear reference genome like hg38.  It operates via
# a Makefile interface (ntLink scaffold target=<assembly> reads=<reads>),
# outputs a scaffolded FASTA and a PAF-like mapping TSV, and has no concept
# of producing a CRAM, BAM, or SAM against a reference genome.
#
# It therefore cannot fulfill the Art. VIII.2 output requirements:
#   cram/<dataset>.hg38.ntlink-ont.cram
#   cram/<dataset>.hg38.ntlink-ont.cram.crai
#   etc.
#
# Possible resolutions — please choose one and update this file:
#
#   (A) The intended tool is a different aligner with a similar name (e.g.
#       ntJoin, ntHits, or another BCGSC tool that does reference mapping).
#       → Provide the correct tool name and GitHub URL.
#
#   (B) ntLink is in your benchmarking project in a different capacity
#       (e.g. assembly quality evaluation, not read mapping). In that case
#       the mapper_tag "ntlink-ont" should be removed from Art. IV.4 and
#       this smk file should not exist.
#
#   (C) You want to benchmark ntLink's *pair*-mode PAF mappings indirectly
#       by piping them through a conversion step. If so, the workflow
#       approach changes substantially from the standard CRAM pipeline.
#       → Describe the intended evaluation approach.
#
# This file is intentionally a stub so that snakemake -n will fail loudly.
#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

include: "header.smk"

raise WorkflowError(
    "ont.read_mapping.ntlink.smk: ntLink is a genome assembly scaffolder, "
    "not a read-to-reference mapper. This workflow cannot be executed until the "
    "correct tool or use-case is identified. See the BLOCKER comment at the top "
    "of this file for resolution options."
)
