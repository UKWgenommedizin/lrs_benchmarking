# ============================================================
# lrs_benchmarking — Mapper Docker Images & CLI References
# ============================================================
# All mappers used in the new benchmarking workflows.
# Docker images are from schimar/* on Docker Hub.
# Run snakemake -s <workflow.smk> to execute.
# ============================================================

# ---------- minimap2 (ONT: mm2-ont, PB: mm2-pb) ----------
docker run --rm -it -v $PWD:/data --entrypoint minimap2 schimar/lrs-minimap2-ntlink:latest --help
# version 2.28-r1209
# ONT: snakemake -s ont.read_mapping.minimap2.smk
# PB:  snakemake -s pb.read_mapping.minimap2.smk

# ---------- ntLink (ONT: ntlink-ont, PB: ntlink-pb) ----------
docker run --rm -it -v $PWD:/data --entrypoint ntLink schimar/lrs-minimap2-ntlink:latest --help
# v1.3.9
# ONT: snakemake -s ont.read_mapping.ntlink.smk
# PB:  snakemake -s pb.read_mapping.ntlink.smk

# ---------- ParaHAT (ONT: parahat-ont, PB: parahat-pb) ----------
# Index first:
docker run --rm -it -v $PWD:/data --entrypoint ParaHAT-indexer schimar/lrs-parahat:latest <ref.fa> <index_dir>
# v0.1.1
# Align:
docker run --rm -it -v $PWD:/data schimar/lrs-parahat:latest -n 4 ParaHAT-aligner -t 8 <index_dir> <reads.fastq> <ref.fa>
# 0.1.1 
# ONT: snakemake -s ont.read_mapping.parahat.smk
# PB:  snakemake -s pb.read_mapping.parahat.smk

# ---------- QuickEd (ONT: quicked-ont, PB: quicked-pb) ----------
docker run --rm -it -v $PWD:/data schimar/lrs-quicked:latest --help
docker run --rm -it -v $PWD:/data --entrypoint generate_dataset schimar/lrs-quicked:latest --help
# ONT: snakemake -s ont.read_mapping.quicked.smk
# PB:  snakemake -s pb.read_mapping.quicked.smk

# ---------- VACmap (ONT: vacmap-ont, PB: vacmap-pb) ----------
docker run --rm -it -v $PWD:/data schimar/lrs-vacmap:latest --help
# ONT: snakemake -s ont.read_mapping.vacmap.smk
# PB:  snakemake -s pb.read_mapping.vacmap.smk

# ---------- GraphAligner (ONT: ga-ont, PB: ga-pb) ----------
docker run --rm -it -v $PWD:/data schimar/lrs-graphaligner:latest --help
# ONT: snakemake -s ont.read_mapping.graphaligner.smk
# PB:  snakemake -s pb.read_mapping.graphaligner.smk

# ---------- VG (ONT: vg-ont, PB: vg-pb) ----------
docker run --rm -it -v $PWD:/data schimar/lrs-vg:latest help
docker run --rm -it -v $PWD:/data schimar/lrs-vg:latest map --help
# ONT: snakemake -s ont.read_mapping.vg.smk
# PB:  snakemake -s pb.read_mapping.vg.smk
# Note: VG requires pre-built graph indices. Run `snakemake -s ont.read_mapping.vg.smk build_vg_index` first.