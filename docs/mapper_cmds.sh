
# ============================================================
# lrs_benchmarking — Mapper Docker Images & CLI References
# ============================================================
# Docker images are from schimar/* on Docker Hub.
# Run snakemake -s <workflow.smk> to execute.
# ============================================================

# ---------- minimap2 (ONT: mm2-ont, PB: mm2-pb) ----------
docker run --rm -it -v $PWD:/data --entrypoint minimap2 schimar/lrs-minimap2-ntlink:v2.28-1.3.9 --help
# ONT: snakemake -s ont.read_mapping.minimap2.smk
# PB:  snakemake -s pb.read_mapping.minimap2.smk

# ---------- pbmm2 (ONT: pbmm2-ont, PB: pbmm2-pb) ----------
docker run --rm -it -v $PWD:/data --entrypoint pbmm2 schimar/lrs-pbmm2:v1.17.0 --help
# ONT: snakemake -s ont.read_mapping.pbmm2.smk
# PB:  snakemake -s pb.read_mapping.pbmm2.smk

# ---------- VACmap (ONT: vacmap-ont, PB: vacmap-pb) ----------
docker run --rm -it -v $PWD:/data schimar/lrs-vacmap:v1.2.0 --help
# ONT: snakemake -s ont.read_mapping.vacmap.smk
# PB:  snakemake -s pb.read_mapping.vacmap.smk

# ---------- VG (ONT: vg-ont, PB: vg-pb) ----------
docker run --rm -it -v $PWD:/data schimar/lrs-vg:1.73.0 help
docker run --rm -it -v $PWD:/data schimar/lrs-vg:1.73.0 map --help
# ONT: snakemake -s ont.read_mapping.vg.smk
# PB:  snakemake -s pb.read_mapping.vg.smk
# Note: VG requires pre-built graph indices. Run `snakemake -s ont.read_mapping.vg.smk build_vg_index` first.

# ============================================================
# [LEGACY] — Deprecated mappers, not maintained
# ============================================================

# ---------- GraphAligner (ONT: ga-ont, PB: ga-pb) ----------
# SIGSEGV on full-genome hg38 — shelved. Workflows in mapper_legacy/.
docker run --rm -it -v $PWD:/data schimar/lrs-graphaligner:v1.0.20 --help
# ONT: snakemake -s mapper_legacy/ont.read_mapping.graphaligner.smk
# PB:  snakemake -s mapper_legacy/pb.read_mapping.graphaligner.smk

# ---------- ParaHAT (ONT: parahat-ont, PB: parahat-pb) ----------
# SIGSEGV on full-genome hg38 — shelved. Workflows in mapper_legacy/.
docker run --rm -it -v $PWD:/data --entrypoint ParaHAT-indexer schimar/lrs-parahat:v1.0.0-cuda <ref.fa> <index_dir>
docker run --rm -it -v $PWD:/data schimar/lrs-parahat:v1.0.0-cuda -n 4 ParaHAT-aligner -t 8 <index_dir> <reads.fastq> <ref.fa>
# ONT: snakemake -s mapper_legacy/ont.read_mapping.parahat.smk
# PB:  snakemake -s mapper_legacy/pb.read_mapping.parahat.smk
