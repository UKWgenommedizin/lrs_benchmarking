#####################################################################################
# AUTOMATIC CHROMOSOME-21 PREPROCESSING
#
# SERVER MODE ONLY
#
# whole-genome 30x FASTQ
#       ↓
# minimap2
#       ↓
# chr21 + primary + MAPQ filtering
#       ↓
# unique primary read selection
#       ↓
# filtered Chr21 BAM validation
#       ↓
# deterministic whole-read normalization to 30x
#       ↓
# results/chr21/{sample}.{technology}.fastq.gz
#
# Local/preextracted mode does not execute any rules in this module.
#####################################################################################


if INPUT_MODE == "whole_genome":


    # -------------------------------------------------------------------------
    # Map WGS reads and retain validated primary Chr21 alignments
    # -------------------------------------------------------------------------

    rule prepare_chr21_bam:
        input:
            validation=VALIDATION_OK,
            fastq=source_fastq_for,
            reference=lambda wildcards: str(REFERENCE)

        output:
            bam=temp(
                "results/work/chr21/"
                "{sample}.{technology}.bam"
            ),
            bai=temp(
                "results/work/chr21/"
                "{sample}.{technology}.bam.bai"
            ),
            mapped=(
                f"{PROGRESS_DIR}/chr21/"
                "{sample}.{technology}.mapped.done"
            )

        params:
            preset=lambda wildcards: (
                "map-ont"
                if wildcards.technology == "ont"
                else "map-hifi"
            ),
            chromosome=CHROMOSOME,
            chromosome_length=CHROMOSOME_LENGTH,
            minimum_mapq=MIN_MAPQ,
            excluded_flags=2308,
            map_threads=MAP_THREADS,
            sort_threads=SORT_THREADS,
            running=(
                f"{PROGRESS_DIR}/chr21/"
                "{sample}.{technology}.mapping.running"
            ),
            name_tmp=(
                "results/work/chr21/tmp/"
                "{sample}.{technology}.name"
            ),
            coordinate_tmp=(
                "results/work/chr21/tmp/"
                "{sample}.{technology}.coordinate"
            )

        threads:
            MAPPING_THREADS

        resources:
            mem_mb=MAPPING_MEM_MB

        log:
            (
                f"{LOG_DIR}/chr21/"
                "{sample}.{technology}.mapping.log"
            )

        shell:
            r"""
            set -euo pipefail

            mkdir -p \
                "$(dirname '{output.bam}')" \
                "$(dirname '{output.mapped}')" \
                "$(dirname '{log}')" \
                "results/work/chr21/tmp"

            rm -f \
                "{output.bam}" \
                "{output.bai}" \
                "{output.mapped}" \
                "{params.running}" \
                "{output.bam}.partial" \
                "{output.bam}.partial.bai"

            touch "{params.running}"

            trap \
                'rc=$?; rm -f "{params.running}"; exit $rc' \
                EXIT

            echo "============================================================" \
                > "{log}"
            echo "CHR21 MAPPING" \
                >> "{log}"
            echo "sample={wildcards.sample}" \
                >> "{log}"
            echo "technology={wildcards.technology}" \
                >> "{log}"
            echo "preset={params.preset}" \
                >> "{log}"
            echo "chromosome={params.chromosome}" \
                >> "{log}"
            echo "minimum_mapq={params.minimum_mapq}" \
                >> "{log}"
            echo "input={input.fastq}" \
                >> "{log}"
            echo "reference={input.reference}" \
                >> "{log}"
            echo "============================================================" \
                >> "{log}"

            minimap2 \
                -ax "{params.preset}" \
                -t "{params.map_threads}" \
                "{input.reference}" \
                "{input.fastq}" \
                2>> "{log}" \
            | python3 \
                scripts/chr21/filter_sam_start_window.py \
                "{params.chromosome}" \
                1 \
                "{params.chromosome_length}" \
                "{params.excluded_flags}" \
                "{params.minimum_mapq}" \
                2>> "{log}" \
            | samtools view \
                -b \
                - \
                2>> "{log}" \
            | samtools sort \
                -n \
                -@ "{params.sort_threads}" \
                -T "{params.name_tmp}" \
                -O SAM \
                - \
                2>> "{log}" \
            | python3 \
                scripts/chr21/select_one_primary_per_qname.py \
                2>> "{log}" \
            | samtools sort \
                -@ "{params.sort_threads}" \
                -T "{params.coordinate_tmp}" \
                -o "{output.bam}.partial" \
                - \
                2>> "{log}"

            samtools quickcheck \
                "{output.bam}.partial"

            samtools index \
                -@ "{params.sort_threads}" \
                "{output.bam}.partial"

            samtools idxstats \
                "{output.bam}.partial" \
                >> "{log}"

            mv \
                "{output.bam}.partial" \
                "{output.bam}"

            mv \
                "{output.bam}.partial.bai" \
                "{output.bai}"

            printf \
                "MAPPING COMPLETE\nsample=%s\ntechnology=%s\nSTATUS=PASS\n" \
                "{wildcards.sample}" \
                "{wildcards.technology}" \
                > "{output.mapped}"

            rm -f "{params.running}"

            trap - EXIT
            """


    # -------------------------------------------------------------------------
    # Validate the extracted Chr21 BAM scientifically
    # -------------------------------------------------------------------------

    rule validate_chr21_bam:
        input:
            bam=(
                "results/work/chr21/"
                "{sample}.{technology}.bam"
            ),
            bai=(
                "results/work/chr21/"
                "{sample}.{technology}.bam.bai"
            ),
            mapped=(
                f"{PROGRESS_DIR}/chr21/"
                "{sample}.{technology}.mapped.done"
            )

        output:
            summary=(
                "results/validation/chr21/"
                "{sample}.{technology}.tsv"
            ),
            filtered=(
                f"{PROGRESS_DIR}/chr21/"
                "{sample}.{technology}.filtered.done"
            )

        params:
            chromosome=CHROMOSOME,
            chromosome_length=CHROMOSOME_LENGTH,
            minimum_mapq=MIN_MAPQ,
            target_coverage=TARGET_COVERAGE

        log:
            (
                f"{LOG_DIR}/chr21/"
                "{sample}.{technology}.validation.log"
            )

        shell:
            r"""
            set -euo pipefail

            mkdir -p \
                "$(dirname '{output.summary}')" \
                "$(dirname '{output.filtered}')" \
                "$(dirname '{log}')"

            rm -f \
                "{output.summary}" \
                "{output.filtered}"

            python3 \
                scripts/chr21/validate_extracted_chr21_bam.py \
                --bam "{input.bam}" \
                --summary "{output.summary}" \
                --checkpoint "{output.filtered}" \
                --sample "{wildcards.sample}" \
                --technology "{wildcards.technology}" \
                --chromosome "{params.chromosome}" \
                --chromosome-length "{params.chromosome_length}" \
                --minimum-mapq "{params.minimum_mapq}" \
                --target-coverage "{params.target_coverage}" \
                > "{log}" \
                2>&1
            """


    # -------------------------------------------------------------------------
    # Deterministically normalize filtered Chr21 reads to 30x
    # -------------------------------------------------------------------------

    rule normalize_chr21_to_30x:
        input:
            bam=(
                "results/work/chr21/"
                "{sample}.{technology}.bam"
            ),
            filtered=(
                f"{PROGRESS_DIR}/chr21/"
                "{sample}.{technology}.filtered.done"
            )

        output:
            fastq=(
                f"{CHR21_OUTPUT_DIR}/"
                "{sample}.{technology}.fastq.gz"
            ),
            names=(
                "results/provenance/chr21/"
                "{sample}.{technology}.selected_names.txt"
            ),
            summary=(
                "results/validation/chr21/"
                "{sample}.{technology}.normalized.tsv"
            ),
            normalized=(
                f"{PROGRESS_DIR}/chr21/"
                "{sample}.{technology}.normalized.done"
            )

        params:
            chromosome=CHROMOSOME,
            chromosome_length=CHROMOSOME_LENGTH,
            target_coverage=TARGET_COVERAGE,
            seed=lambda wildcards: seed_for_values(
                wildcards.sample,
                wildcards.technology
            )

        log:
            (
                f"{LOG_DIR}/chr21/"
                "{sample}.{technology}.normalization.log"
            )

        shell:
            r"""
            set -euo pipefail

            mkdir -p \
                "$(dirname '{output.fastq}')" \
                "$(dirname '{output.names}')" \
                "$(dirname '{output.summary}')" \
                "$(dirname '{output.normalized}')" \
                "$(dirname '{log}')"

            # Makes the rule safe to rerun if a previous attempt produced only
            # part of its multi-file output set.
            rm -f \
                "{output.fastq}" \
                "{output.names}" \
                "{output.summary}" \
                "{output.normalized}" \
                "{output.fastq}.partial" \
                "{output.names}.partial"

            python3 \
                scripts/chr21/normalize_chr21_bam_to_coverage.py \
                --input-bam "{input.bam}" \
                --output-fastq "{output.fastq}" \
                --selected-names "{output.names}" \
                --summary "{output.summary}" \
                --checkpoint "{output.normalized}" \
                --log "{log}" \
                --sample "{wildcards.sample}" \
                --technology "{wildcards.technology}" \
                --chromosome "{params.chromosome}" \
                --chromosome-length "{params.chromosome_length}" \
                --target-coverage "{params.target_coverage}" \
                --seed "{params.seed}" \
                --selection-policy primary_MAPQ20_source
            """
