#####################################################################################
# INPUT AND ENVIRONMENT VALIDATION
#
# Validates:
#   - portable sample metadata
#   - automatically resolved source FASTQs
#   - server reference when required
#   - required Conda environment definitions
#####################################################################################


rule validate_inputs:
    output:
        VALIDATION_OK

    run:
        Path(output[0]).parent.mkdir(
            parents=True,
            exist_ok=True
        )

        errors = []

        # ---------------------------------------------------------------------
        # Sample sheet
        # ---------------------------------------------------------------------

        sample_sheet_path = Path(
            SAMPLE_SHEET
        )

        if not sample_sheet_path.is_absolute():
            sample_sheet_path = (
                PROJECT_DIR
                / sample_sheet_path
            )

        if not sample_sheet_path.exists():
            errors.append(
                f"Missing sample sheet: "
                f"{sample_sheet_path}"
            )


        # ---------------------------------------------------------------------
        # Input root
        # ---------------------------------------------------------------------

        if not INPUT_ROOT.exists():
            errors.append(
                f"Missing input_root directory: "
                f"{INPUT_ROOT}"
            )


        # ---------------------------------------------------------------------
        # Automatically resolved SOURCE FASTQs
        # ---------------------------------------------------------------------

        for row in SAMPLE_ROWS:
            sample = row["sample"]
            technology = row["technology"]

            fastq = Path(
                source_fastq_for_values(
                    sample,
                    technology
                )
            )

            if not fastq.exists():
                errors.append(
                    "Missing source FASTQ for "
                    f"{sample} {technology}: "
                    f"{fastq}"
                )

            elif not fastq.is_file():
                errors.append(
                    "Input is not a regular FASTQ file for "
                    f"{sample} {technology}: "
                    f"{fastq}"
                )


        # ---------------------------------------------------------------------
        # Reference genome required only for WGS/server mode
        # ---------------------------------------------------------------------

        if INPUT_MODE == "whole_genome":
            raw_reference = config.get(
                "reference"
            )

            if not raw_reference:
                errors.append(
                    "Server/whole_genome mode requires "
                    "'reference' in config/server.yaml"
                )

            else:
                reference_path = Path(
                    raw_reference
                ).expanduser()

                if not reference_path.is_absolute():
                    reference_path = (
                        PROJECT_DIR
                        / reference_path
                    )

                if not reference_path.exists():
                    errors.append(
                        "Missing reference FASTA: "
                        f"{reference_path}"
                    )


        # ---------------------------------------------------------------------
        # Required assembler environments
        # ---------------------------------------------------------------------

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

        for environment_file in required_environments:
            if not Path(
                environment_file
            ).exists():
                errors.append(
                    "Missing Conda environment file: "
                    f"{environment_file}"
                )


        # ---------------------------------------------------------------------
        # Report all errors together
        # ---------------------------------------------------------------------

        if errors:
            message = "\n".join(
                f"- {error}"
                for error in errors
            )

            raise RuntimeError(
                "\n"
                "ASSEMBLER WORKFLOW VALIDATION FAILED\n"
                "\n"
                f"Input mode: {INPUT_MODE}\n"
                f"Input root: {INPUT_ROOT}\n"
                "\n"
                "Fix these problems before running "
                "the workflow:\n\n"
                f"{message}\n"
            )


        # ---------------------------------------------------------------------
        # Validation marker
        # ---------------------------------------------------------------------

        Path(output[0]).write_text(
            "Assembler input validation passed.\n"
            f"input_mode={INPUT_MODE}\n"
            f"input_root={INPUT_ROOT}\n"
            f"samples={len(SAMPLES)}\n"
            f"datasets={len(SAMPLE_TECH_PAIRS)}\n"
        )
