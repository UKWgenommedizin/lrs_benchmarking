#####################################################################################
# INPUT AND ENVIRONMENT VALIDATION
# This module checks the sample table, sequencing files, and Conda environment
# definitions before any assembler or scaffolding rule is executed.
#####################################################################################


# -----------------------------------------------------------------------------
# Rule 1: Validate workflow inputs
# -----------------------------------------------------------------------------

# Check the required sample table, FASTQ files, and environment definitions to
# prevent the workflow from failing later because of an incomplete setup
rule validate_inputs:
    output:
        VALIDATION_OK
    run:
        Path(output[0]).parent.mkdir(
            parents=True,
            exist_ok=True
        )

        errors = []

        # Verify that samples.tsv exists
        if not Path(SAMPLE_SHEET).exists():
            errors.append(
                f"Missing sample sheet: {SAMPLE_SHEET}"
            )

        # Verify that every FASTQ declared in samples.tsv exists
        for row in SAMPLE_ROWS:
            fastq = Path(row["fastq"])

            if not fastq.exists():
                errors.append(
                    f"Missing FASTQ file: {fastq}"
                )

        # Build the list of required environments from active_assemblers
        required_environments = []

        if "flye" in ACTIVE_ASSEMBLERS:
            required_environments.append(
                ENV_FLYE
            )

        if "goldrush" in ACTIVE_ASSEMBLERS:
            required_environments.append(
                ENV_GOLDRUSH
            )

        if "ntlink" in ACTIVE_ASSEMBLERS:
            required_environments.append(
                ENV_NTLINK
            )

        if "verkko" in ACTIVE_ASSEMBLERS:
            required_environments.append(
                ENV_VERKKO
            )

        # Verify that every required environment YAML file exists
        for environment_file in required_environments:
            if not Path(environment_file).exists():
                errors.append(
                    f"Missing Conda environment file: "
                    f"{environment_file}"
                )

        # Stop the workflow and display all detected problems together
        if errors:
            message = "\n".join(errors)

            raise RuntimeError(
                "\nASSEMBLER WORKFLOW VALIDATION FAILED\n"
                "Fix these problems before running the workflow:\n\n"
                f"{message}\n"
            )

        # Create the marker file after all validations pass
        Path(output[0]).write_text(
            "Assembler input validation passed.\n"
        )
