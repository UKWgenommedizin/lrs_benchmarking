# Docker image validation status

A static scan of the Snakemake workflow files was performed to identify Docker image usage in the variant-calling and benchmarking workflows.

The scan detected public Docker images such as Clair3, DeepVariant, nf-core bcl2fastq, and Ensembl VEP, as well as private institutional images hosted under `storage-node:5000`.

Local Docker execution was not completed because the `docker` command was not available inside the WSL environment at the time of testing.

This does not block Snakemake dry-run validation, because dry-run mode validates workflow structure, rules, inputs, outputs, parameters, and shell commands without executing Docker containers or processing real data.

Docker pull/run validation should be repeated later in an environment where Docker is available and where the private `storage-node:5000` registry can be accessed.
