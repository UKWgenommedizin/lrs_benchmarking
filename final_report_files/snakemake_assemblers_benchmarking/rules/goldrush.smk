###############################################################################
# GOLDRUSH ASSEMBLY RULES
###############################################################################

rule goldrush_assembly:
    input:
        validation=VALIDATION_OK,
        reads=lambda wildcards: fastq_for(wildcards)

    output:
        assembly=(
            "results/smoke_test/goldrush/"
            "{sample}.{technology}/assembly.fasta"
        ),
        done=(
            f"{PROGRESS_DIR}/goldrush/"
            "{sample}.{technology}.done"
        )

    params:
        workdir=lambda wildcards: (
            f"results/smoke_test/goldrush/"
            f"{wildcards.sample}.{wildcards.technology}/work"
        ),
        prefix=lambda wildcards: (
            f"{wildcards.sample}_{wildcards.technology}_goldrush"
        ),
        reads_basename=lambda wildcards: (
            f"{wildcards.sample}_{wildcards.technology}"
        ),
        min_read_length=lambda wildcards: (
            5000 if wildcards.technology == "ont" else 10000
        ),
        running=lambda wc: (
            f"{PROGRESS_DIR}/goldrush/"
            f"{wc.sample}.{wc.technology}.running"
        ),
        failed=lambda wc: (
            f"{PROGRESS_DIR}/goldrush/"
            f"{wc.sample}.{wc.technology}.failed"
        )

    threads:
        GOLDRUSH_THREADS

    resources:
        mem_mb=GOLDRUSH_MEM_MB

    conda:
        ENV_GOLDRUSH

    log:
        (
            f"{LOG_DIR}/goldrush/"
            "{sample}.{technology}.log"
        )

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

        gzip -cd "{input.reads}" \
            > "{params.workdir}/{params.reads_basename}.fastq"

        cd "{params.workdir}"

        goldrush run \
            reads="{params.reads_basename}" \
            G={config[genome_size]} \
            t={threads} \
            m={params.min_read_length} \
            P=0 \
            p="{params.prefix}" \
            track_time=1

        FINAL_ASSEMBLY=$(
            find goldrush_intermediate_files \
                -type f \
                -name "*ntLink-5rounds.polished.fa" \
                -print \
                | sort \
                | tail -n 1
        )

        if [[ -z "$FINAL_ASSEMBLY" || ! -s "$FINAL_ASSEMBLY" ]]; then
            echo "ERROR: GoldRush final polished assembly was not produced."
            exit 1
        fi

        cp "$FINAL_ASSEMBLY" \
           "$OLDPWD/{output.assembly}"

        test -s "$OLDPWD/{output.assembly}"

        printf \
            "STATUS=PASS\nsample=%s\ntechnology=%s\ncompleted=%s\n" \
            "{wildcards.sample}" \
            "{wildcards.technology}" \
            "$(date -Is)" \
            > "$OLDPWD/{output.done}"

        rm -f \
            "$OLDPWD/{params.running}" \
            "$OLDPWD/{params.failed}"

        trap - EXIT
        """
