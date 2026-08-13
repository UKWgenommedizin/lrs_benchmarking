#####################################################################################
# NTLINK SCAFFOLDING RULES
#####################################################################################

rule ntlink_scaffold:
    input:
        validation=VALIDATION_OK,
        draft="results/goldrush/{sample}.{technology}.fasta",
        reads=lambda wildcards: fastq_for(wildcards)

    output:
        assembly="results/ntlink/{sample}.{technology}.fasta",
        done=PROGRESS_DIR + "/ntlink/{sample}.{technology}.done"

    threads:
        NTLINK_THREADS

    resources:
        mem_mb=NTLINK_MEM_MB

    conda:
        ENV_NTLINK

    params:
        workdir="results/work/ntlink/{sample}.{technology}",
        running=lambda wc: (
            f"{PROGRESS_DIR}/ntlink/"
            f"{wc.sample}.{wc.technology}.running"
        ),
        failed=lambda wc: (
            f"{PROGRESS_DIR}/ntlink/"
            f"{wc.sample}.{wc.technology}.failed"
        )

    log:
        LOG_DIR + "/ntlink/{sample}.{technology}.log"

    shell:
        r"""
        set -euo pipefail

        mkdir -p \
            "{params.workdir}" \
            "$(dirname "{output.assembly}")" \
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

        cp "{input.draft}" \
           "{params.workdir}/draft.fasta"

        gzip -cd "{input.reads}" \
            > "{params.workdir}/reads.fastq"

        (
            cd "{params.workdir}"

            ntLink scaffold \
                target=draft.fasta \
                reads='reads.fastq' \
                t={threads}
        )

        NTLINK_ASSEMBLY=(
            "{params.workdir}/"
            "draft.fasta.k32.w100.z1000.ntLink.scaffolds.fa"
        )

        NTLINK_ASSEMBLY="${{NTLINK_ASSEMBLY[0]}}"

        if [[ ! -s "$NTLINK_ASSEMBLY" ]]; then
            echo "ERROR: ntLink expected scaffold assembly was not produced:"
            echo "$NTLINK_ASSEMBLY"
            exit 1
        fi

        cp "$NTLINK_ASSEMBLY" \
           "{output.assembly}"

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
