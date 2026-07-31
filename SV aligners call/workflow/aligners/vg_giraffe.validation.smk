ONT_FASTQ = config.get("ont_fastq", "samples_try/HG002.ont.1k.fastq.gz")
PB_FASTQ = config.get("pb_fastq", "samples_try/HG002.pb.1k.fastq.gz")
OUTPUT_DIR = config.get("output_dir", "SV aligners call/results/alignments")

VG_GBZ = config.get("vg_gbz", "final_report_files/snakemake_aligners_benchmarking/data/graph/graph.gbz")
VG_MIN = config.get("vg_min", "final_report_files/snakemake_aligners_benchmarking/data/graph/graph.min")
VG_DIST = config.get("vg_dist", "final_report_files/snakemake_aligners_benchmarking/data/graph/graph.dist")
VG_ZIPCODES = config.get("vg_zipcodes", "final_report_files/snakemake_aligners_benchmarking/data/graph/graph.zipcodes")

def zipcodes_arg(wildcards):
    if VG_ZIPCODES:
        return f"-z {VG_ZIPCODES}"
    return ""

rule all:
    input:
        f"{OUTPUT_DIR}/HG002.ont.vg_giraffe.validation.sorted.bam",
        f"{OUTPUT_DIR}/HG002.ont.vg_giraffe.validation.sorted.bam.bai",
        f"{OUTPUT_DIR}/HG002.pb.vg_giraffe.validation.sorted.bam",
        f"{OUTPUT_DIR}/HG002.pb.vg_giraffe.validation.sorted.bam.bai"

rule vg_giraffe_align_ont:
    input:
        fastq=ONT_FASTQ,
        gbz=VG_GBZ,
        minimizer=VG_MIN,
        dist=VG_DIST
    output:
        bam=f"{OUTPUT_DIR}/HG002.ont.vg_giraffe.validation.sorted.bam",
        bai=f"{OUTPUT_DIR}/HG002.ont.vg_giraffe.validation.sorted.bam.bai"
    threads: 4
    params:
        preset="r10",
        zipcodes=zipcodes_arg
    shell:
        r"""
        mkdir -p "{OUTPUT_DIR}"
        vg giraffe -Z "{input.gbz}" -m "{input.minimizer}" -d "{input.dist}" \
          {params.zipcodes} -f "{input.fastq}" -t {threads} -b {params.preset} -o BAM \
          | samtools sort -@ {threads} -o "{output.bam}"
        samtools index "{output.bam}"
        """

rule vg_giraffe_align_pb:
    input:
        fastq=PB_FASTQ,
        gbz=VG_GBZ,
        minimizer=VG_MIN,
        dist=VG_DIST
    output:
        bam=f"{OUTPUT_DIR}/HG002.pb.vg_giraffe.validation.sorted.bam",
        bai=f"{OUTPUT_DIR}/HG002.pb.vg_giraffe.validation.sorted.bam.bai"
    threads: 4
    params:
        preset="hifi",
        zipcodes=zipcodes_arg
    shell:
        r"""
        mkdir -p "{OUTPUT_DIR}"
        vg giraffe -Z "{input.gbz}" -m "{input.minimizer}" -d "{input.dist}" \
          {params.zipcodes} -f "{input.fastq}" -t {threads} -b {params.preset} -o BAM \
          | samtools sort -@ {threads} -o "{output.bam}"
        samtools index "{output.bam}"
        """
