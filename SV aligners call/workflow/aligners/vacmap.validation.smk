REFERENCE = config.get(
    "reference",
    "reference/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta"
)

ONT_FASTQ = config.get(
    "ont_fastq",
    "samples_try/HG002.ont.1k.fastq.gz"
)

PB_FASTQ = config.get(
    "pb_fastq",
    "samples_try/HG002.pb.1k.fastq.gz"
)

OUTPUT_DIR = config.get(
    "output_dir",
    "SV aligners call/results/alignments"
)

VACMAP_BIN = config.get(
    "vacmap_bin",
    "conda run -n vacmap_env vacmap"
)


rule all:
    input:
        f"{OUTPUT_DIR}/HG002.ont.vacmap.validation.sorted.bam",
        f"{OUTPUT_DIR}/HG002.ont.vacmap.validation.sorted.bam.bai",
        f"{OUTPUT_DIR}/HG002.pb.vacmap.validation.sorted.bam",
        f"{OUTPUT_DIR}/HG002.pb.vacmap.validation.sorted.bam.bai"


rule vacmap_align_ont:
    input:
        ref=REFERENCE,
        fastq=ONT_FASTQ
    output:
        bam=f"{OUTPUT_DIR}/HG002.ont.vacmap.validation.sorted.bam",
        bai=f"{OUTPUT_DIR}/HG002.ont.vacmap.validation.sorted.bam.bai"
    threads: 4
    params:
        mode="H",
        platform="ONT",
        sample="HG002",
        vacmap_bin=VACMAP_BIN
    shell:
        r"""
        mkdir -p "{OUTPUT_DIR}"

        {params.vacmap_bin} \
          -ref "{input.ref}" \
          -read "{input.fastq}" \
          -mode {params.mode} \
          -t {threads} \
          --rg-id HG002_ONT \
          --rg-sm {params.sample} \
          --rg-pl {params.platform} \
          --force \
          -o "{output.bam}"

        samtools index "{output.bam}"
        """


rule vacmap_align_pb:
    input:
        ref=REFERENCE,
        fastq=PB_FASTQ
    output:
        bam=f"{OUTPUT_DIR}/HG002.pb.vacmap.validation.sorted.bam",
        bai=f"{OUTPUT_DIR}/HG002.pb.vacmap.validation.sorted.bam.bai"
    threads: 4
    params:
        mode="L",
        platform="PACBIO",
        sample="HG002",
        vacmap_bin=VACMAP_BIN
    shell:
        r"""
        mkdir -p "{OUTPUT_DIR}"

        {params.vacmap_bin} \
          -ref "{input.ref}" \
          -read "{input.fastq}" \
          -mode {params.mode} \
          -t {threads} \
          --rg-id HG002_PB \
          --rg-sm {params.sample} \
          --rg-pl {params.platform} \
          --force \
          -o "{output.bam}"

        samtools index "{output.bam}"
        """
