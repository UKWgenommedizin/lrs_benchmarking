#####################################################################################
# FLYE ASSEMBLY RULES
# Flye performs de novo assembly independently for ONT and PacBio HiFi reads.
#####################################################################################


# -----------------------------------------------------------------------------
# Rule 2: Flye de novo assembly
# -----------------------------------------------------------------------------

# Use the technology-specific Flye input option and generate one assembly
# directory and one assembly FASTA for every sample and technology combination
rule flye_assembly:
    input:
        validation=VALIDATION_OK,
        fastq=fastq_for
    output:
        assembly=(
            "results/assemblies/"
            "{sample}.{technology}.flye/"
            "assembly.fasta"
        )
    conda:
        ENV_FLYE
    params:
        outdir=(
            "results/assemblies/"
            "{sample}.{technology}.flye"
        ),
        read_option=flye_read_option,
        genome_size=GENOME_SIZE
    threads:
        DEFAULT_THREADS
    shell:
        r"""
        mkdir -p {params.outdir}

        flye \
            {params.read_option} {input.fastq} \
            --genome-size {params.genome_size} \
            --threads {threads} \
            --out-dir {params.outdir}
        """
