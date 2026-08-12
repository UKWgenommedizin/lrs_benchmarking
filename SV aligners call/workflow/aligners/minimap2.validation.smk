REFERENCE = config.get("reference", "reference/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta")
ONT_FASTQ = config.get("ont_fastq", "samples_try/HG002.ont.1k.fastq.gz")
PB_FASTQ = config.get("pb_fastq", "samples_try/HG002.pb.1k.fastq.gz")
OUTPUT_DIR = config.get("output_dir", "SV aligners call/results/alignments")

rule all:
    input:
        f"{OUTPUT_DIR}/HG002.ont.minimap2.validation.sorted.bam",
        f"{OUTPUT_DIR}/HG002.ont.minimap2.validation.sorted.bam.bai",
        f"{OUTPUT_DIR}/HG002.pb.minimap2.validation.sorted.bam",
        f"{OUTPUT_DIR}/HG002.pb.minimap2.validation.sorted.bam.bai"

rule minimap2_map_ont:
    input:
        ref=REFERENCE,
        fastq=ONT_FASTQ
    output:
        bam=f"{OUTPUT_DIR}/HG002.ont.minimap2.validation.sorted.bam",
        bai=f"{OUTPUT_DIR}/HG002.ont.minimap2.validation.sorted.bam.bai"
    threads: 4
    params:
        preset="map-ont"
    shell:
        r"""
        mkdir -p "{OUTPUT_DIR}"
        minimap2 -ax {params.preset} -t {threads} "{input.ref}" "{input.fastq}" \
          | samtools sort -@ {threads} -o "{output.bam}"
        samtools index "{output.bam}"
        """

rule minimap2_map_pb:
    input:
        ref=REFERENCE,
        fastq=PB_FASTQ
    output:
        bam=f"{OUTPUT_DIR}/HG002.pb.minimap2.validation.sorted.bam",
        bai=f"{OUTPUT_DIR}/HG002.pb.minimap2.validation.sorted.bam.bai"
    threads: 4
    params:
        preset="map-hifi"
    shell:
        r"""
        mkdir -p "{OUTPUT_DIR}"
        minimap2 -ax {params.preset} -t {threads} "{input.ref}" "{input.fastq}" \
          | samtools sort -@ {threads} -o "{output.bam}"
        samtools index "{output.bam}"
        """
