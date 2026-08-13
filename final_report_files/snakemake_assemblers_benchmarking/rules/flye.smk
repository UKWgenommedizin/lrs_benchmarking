#####################################################################################
# FLYE ASSEMBLY RULES
# Flye performs de novo assembly independently for ONT and PacBio HiFi reads.
#####################################################################################

rule flye_assembly:
    input:
        validation=VALIDATION_OK,
        fastq=fastq_for

    output:
        assembly=(
            "results/assemblies/"
            "{sample}.{technology}.flye/"
            "assembly.fasta"
        ),
        done=(
            f"{PROGRESS_DIR}/flye/"
            "{sample}.{technology}.done"
        )

    conda:
        ENV_FLYE

    params:
        outdir=(
            "results/assemblies/"
            "{sample}.{technology}.flye"
        ),
        read_option=flye_read_option,
        genome_size=GENOME_SIZE,
        running=lambda wc: (
            f"{PROGRESS_DIR}/flye/"
            f"{wc.sample}.{wc.technology}.running"
        ),
        failed=lambda wc: (
            f"{PROGRESS_DIR}/flye/"
            f"{wc.sample}.{wc.technology}.failed"
        )

    threads:
        FLYE_THREADS

    resources:
        mem_mb=FLYE_MEM_MB

    log:
        (
            f"{LOG_DIR}/flye/"
            "{sample}.{technology}.log"
        )

    shell:
        r"""
        set -euo pipefail

        mkdir -p \
            "{params.outdir}" \
            "$(dirname "{output.done}")" \
            "$(dirname "{log}")"

        rm -f \
            "{output.done}" \
            "{params.running}" \
            "{params.failed}"

        printf \
            "STATUS=RUNNING\nsample=%s\ntechnology=%s\nstarted=%s\n" \
            "{wildcards.sample}" \
            "{wildcards.technology}" \
            "$(date -Is)" \
            > "{params.running}"

        trap 'rc=$?;
              rm -f "{params.running}";
              if [[ $rc -ne 0 ]]; then
                  printf "STATUS=FAILED\nsample=%s\ntechnology=%s\nfailed=%s\n" \
                      "{wildcards.sample}" \
                      "{wildcards.technology}" \
                      "$(date -Is)" \
                      > "{params.failed}";
              fi;
              exit $rc' EXIT

        exec > "{log}" 2>&1

        flye \
            {params.read_option} "{input.fastq}" \
            --genome-size "{params.genome_size}" \
            --threads {threads} \
            --out-dir "{params.outdir}"

        test -s "{output.assembly}"

        printf \
            "STATUS=PASS\nsample=%s\ntechnology=%s\ncompleted=%s\n" \
            "{wildcards.sample}" \
            "{wildcards.technology}" \
            "$(date -Is)" \
            > "{output.done}"

        rm -f \
            "{params.running}" \
            "{params.failed}"

        trap - EXIT
        """
