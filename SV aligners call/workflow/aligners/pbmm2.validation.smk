REFERENCE = config.get("reference", "reference/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta")
ONT_FASTQ = config.get("ont_fastq", "samples_try/HG002.ont.1k.fastq.gz")
PB_FASTQ = config.get("pb_fastq", "samples_try/HG002.pb.1k.fastq.gz")
OUTPUT_DIR = config.get("output_dir", "SV aligners call/results/alignments")

rule all:
    input:
        f"{OUTPUT_DIR}/HG002.ont.pbmm2.validation.sorted.bam",
        f"{OUTPUT_DIR}/HG002.ont.pbmm2.validation.sorted.bam.bai",
        f"{OUTPUT_DIR}/HG002.pb.pbmm2.validation.sorted.bam",
        f"{OUTPUT_DIR}/HG002.pb.pbmm2.validation.sorted.bam.bai"

rule pbmm2_align_ont:
    input:
        ref=REFERENCE,
        fastq=ONT_FASTQ
    output:
        bam=f"{OUTPUT_DIR}/HG002.ont.pbmm2.validation.sorted.bam",
        bai=f"{OUTPUT_DIR}/HG002.ont.pbmm2.validation.sorted.bam.bai"
    threads: 4
    params:
        preset="SUBREAD"
    shell:
        r"""
        mkdir -p "{OUTPUT_DIR}"
        pbmm2 align --preset {params.preset} --sort -j {threads} \
          "{input.ref}" "{input.fastq}" "{output.bam}"
        samtools index "{output.bam}"
        """

rule pbmm2_align_pb:
    input:
        ref=REFERENCE,
        fastq=PB_FASTQ
    output:
        bam=f"{OUTPUT_DIR}/HG002.pb.pbmm2.validation.sorted.bam",
        bai=f"{OUTPUT_DIR}/HG002.pb.pbmm2.validation.sorted.bam.bai"
    threads: 4
    params:
        preset="HIFI"
    shell:
        r"""
        mkdir -p "{OUTPUT_DIR}"
        pbmm2 align --preset {params.preset} --sort -j {threads} \
          "{input.ref}" "{input.fastq}" "{output.bam}"
        samtools index "{output.bam}"
        """
