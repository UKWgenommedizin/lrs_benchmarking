# ************************************************************************************************
# * Snakefile for SNP calling (test version)
# ************************************************************************************************

import os
import math


#Define variables
GDA_VERSION = "3.0"
DEEPVARIANT_VERSION = "1.9.0"
DB_DIR = "/mnt/storage/db"
FASTA = "references/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta"


#Set FASTA_PATH
FASTA_PATH = DB_DIR + "/" + FASTA


#Get current working directory
CWD = os.getcwd()


# Create wildcards
DATASETS,MAPPERS = glob_wildcards(CWD + "/cram/{dataset, [A-Za-z0-9\\_\\-\\.]+}.hg38.{mapper, [A-Za-z0-9\\-]+}.cram")


#Logging
print("Version of Docker image 'genetic_data_analysis': " + GDA_VERSION)
print("Version of Docker image 'deepvariant': " + DEEPVARIANT_VERSION)
print("Database directory: " + DB_DIR)
print("Path to reference genome: " + FASTA_PATH)

print("Current workding directory: " + CWD)
print("Datasets: " + ','.join(DATASETS))
print("Mappers: " + ','.join(MAPPERS))


#Decide on GPU or CPU only workdlow
if "gpus" in workflow.global_resources:
	ruleorder: deepvariant_gpu_postprocess_variants > deepvariant_cpu_only
	print("Multi-step GPU mode")

else:
	ruleorder: deepvariant_cpu_only > deepvariant_gpu_postprocess_variants
	print("Single-step CPU mode")


#Define parallelization
MAKE_EXAMPLES_SCATTER_INDICES = [ "%03d" % element for element in list(range(0, 64)) ]
CALL_VARIANTS_SCATTER_INDICES = [ "%03d" % element for element in list(range(0, 16)) ]

#print(len(MAKE_EXAMPLES_SCATTER_INDICES))
#print(*MAKE_EXAMPLES_SCATTER_INDICES)

#print(len(CALL_VARIANTS_SCATTER_INDICES))
#print(*CALL_VARIANTS_SCATTER_INDICES)


# *** Define Output
OUTPUT = [];

OUTPUT = OUTPUT + expand(CWD + "/vcf_called/snv_indel/{dataset}.hg38.{mapper}.dv-ont-woPG.vcf.gz", zip, dataset=DATASETS, mapper=MAPPERS)

OUTPUT = OUTPUT + expand(CWD + "/vcf_called/snv_indel/{dataset}.hg38.{mapper}.dv-ont-woPG.visual_report.html", zip, dataset=DATASETS, mapper=MAPPERS)

OUTPUT = OUTPUT + expand(CWD + "/vcf_called/snv_indel/{dataset}.hg38.{mapper}.dv-ont-woPG.vcf.gz.stats", zip, dataset=DATASETS, mapper=MAPPERS)


# ************************************************************************************************


rule all:
    input: OUTPUT

rule test:
    shell:  print(OUTPUT) # print(DATASETS),


# ************************************************************************************************	
# Rules
# ************************************************************************************************		

rule deepvariant_gpu_make_examples:
    input:   cram="{cwd}/cram/{dataset}.hg38.{mapper}.cram", crai="{cwd}/cram/{dataset}.hg38.{mapper}.cram.crai", fasta=FASTA_PATH
    output:  examples=temp(expand("{{cwd}}/vcf_called/snv_indel/tmp/{{dataset}}.hg38.{{mapper}}.dv-ont-woPG/make_examples.tfrecord-00{scatter_index}-of-00064.gz", scatter_index=MAKE_EXAMPLES_SCATTER_INDICES)), gvcf=temp(expand("{{cwd}}/vcf_called/snv_indel/tmp/{{dataset}}.hg38.{{mapper}}.dv-ont-woPG/gvcf.tfrecord-00{scatter_index}-of-00064.gz", scatter_index=MAKE_EXAMPLES_SCATTER_INDICES))
    params:  examples="{cwd}/vcf_called/snv_indel/tmp/{dataset}.hg38.{mapper}.dv-ont-woPG/make_examples.tfrecord@64.gz", gvcf="{cwd}/vcf_called/snv_indel/tmp/{dataset}.hg38.{mapper}.dv-ont-woPG/gvcf.tfrecord@64.gz", log="{cwd}/vcf_called/snv_indel/{dataset}.hg38.{mapper}.dv-ont-woPG.make_examples.log"
    message: "executing {rule} with output {output} and input {input}"
    priority: 6
    threads: len(MAKE_EXAMPLES_SCATTER_INDICES)
    resources:
	    mem_gb=160
    shell:   "umask 0027; \
		mkdir -p $(dirname {params.examples}); \
		srun -p all -c {threads} --mem={resources.mem_gb}GB \
		docker run --cpus {threads} -m {resources.mem_gb}g -u $UID:1002 --workdir /tmp --rm -v {CWD}:{CWD} -v {DB_DIR}:{DB_DIR}:ro google/deepvariant:{DEEPVARIANT_VERSION} /bin/bash -c \" \
			printf 'Container ID:\\t'; hostname; \
			printf 'Start time:\\t'; date; \
			umask 0027; \
			( time seq 0 63 | \
			parallel \
				-q \
				--halt 2 \
				--line-buffer /opt/deepvariant/bin/make_examples \
					--mode calling \
					--ref {input.fasta} \
					--reads {input.cram} \
					--examples {params.examples} \
					--gvcf {params.gvcf} \
					--checkpoint '/opt/models/ont_r104' \
					--alt_aligned_pileup 'diff_channels' \
					--max_reads_per_partition '600' \
					--min_mapping_quality '5' \
					--parse_sam_aux_fields \
					--partition_size '25000' \
					--phase_reads \
					--pileup_image_width '99' \
					--norealign_reads \
					--sort_by_haplotypes \
					--track_ref_reads \
					--trim_reads_for_pileup \
					--vsc_min_fraction_indels '0.12' \
					--vsc_min_fraction_snps '0.08' \
					--task {{}} ); \
			printf 'End time:\\t'; date; \" \
		&> {params.log};"

#--normalize_reads \
#/opt/deepvariant/bin/make_examples --mode calling --ref "/mnt/storage/db/references/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta" --reads "/mnt/storage/genetic_data/WGS_LR/Tests/lrs_benchmarking/cram/HG002_ont_30x_1k.hg38.mm2-pb.cram" --examples "/mnt/storage/genetic_data/WGS_LR/Tests/lrs_benchmarking/vcf_called/deepvariant/tmp/HG002_ont_30x_1k.hg38.mm2-pb/intermediate_results_dir/make_examples.tfrecord@64.gz" --checkpoint "/opt/models/ont_r104" --alt_aligned_pileup "diff_channels" --call_small_model_examples --gvcf "/mnt/storage/genetic_data/WGS_LR/Tests/lrs_benchmarking/vcf_called/deepvariant/tmp/HG002_ont_30x_1k.hg38.mm2-pb/intermediate_results_dir/gvcf.tfrecord@64.gz" --max_reads_per_partition "600" --min_mapping_quality "5" --parse_sam_aux_fields --partition_size "25000" --phase_reads --pileup_image_width "99" --norealign_reads --regions "chr1:0-10000000" --small_model_indel_gq_threshold "17" --small_model_snp_gq_threshold "9" --small_model_vaf_context_window_size "51" --sort_by_haplotypes --track_ref_reads --trained_small_model_path "/opt/smallmodels/ont_r104" --trim_reads_for_pileup --vsc_min_fraction_indels "0.12" --vsc_min_fraction_snps "0.08"


rule deepvariant_gpu_call_variants:
    input:   examples=expand("{{cwd}}/vcf_called/snv_indel/tmp/{{dataset}}.hg38.{{mapper}}.dv-ont-woPG/make_examples.tfrecord-00{scatter_index}-of-00064.gz", scatter_index=MAKE_EXAMPLES_SCATTER_INDICES)
    output:  call_variants_output=temp(expand("{{cwd}}/vcf_called/snv_indel/tmp/{{dataset}}.hg38.{{mapper}}.dv-ont-woPG/call_variants-00{scatter_index}-of-00016.tfrecord.gz", scatter_index=CALL_VARIANTS_SCATTER_INDICES))
    params:  examples="{cwd}/vcf_called/snv_indel/tmp/{dataset}.hg38.{mapper}.dv-ont-woPG/make_examples.tfrecord@64.gz", call_variants_output="{cwd}/vcf_called/snv_indel/tmp/{dataset}.hg38.{mapper}.dv-ont-woPG/call_variants.tfrecord.gz", log="{cwd}/vcf_called/snv_indel/{dataset}.hg38.{mapper}.dv-ont-woPG.call_variants.log"
    message: "executing {rule} with output {output} and input {input}"
    priority: 10
    threads: 32
    resources:
	    gpus=1,
	    mem_gb=60
    shell:   "umask 0027; \
		mkdir -p $(dirname {params.call_variants_output}); \
		srun -p gpu-partition -c {threads} --mem={resources.mem_gb}GB --gpus 1 --priority='TOP' /bin/bash -c \" \
		docker run --cpus {threads} -m {resources.mem_gb}g --runtime=nvidia -e NVIDIA_VISIBLE_DEVICES=\$CUDA_VISIBLE_DEVICES -u $UID:1002 --workdir /tmp --rm -v {CWD}:{CWD} google/deepvariant:{DEEPVARIANT_VERSION}-gpu /bin/bash -c \\\" \
			printf 'Container ID:\\t'; hostname; \
			printf 'Start time:\\t'; date; \
			umask 0027; \
			( time \
			/opt/deepvariant/bin/call_variants \
				--examples {params.examples} \
				--outfile {params.call_variants_output} \
				--checkpoint '/opt/models/ont_r104' ); \
			printf 'End time:\\t'; date; \\\" \" \
		&> {params.log};"

#/opt/deepvariant/bin/call_variants --outfile "/mnt/storage/genetic_data/WGS_LR/Tests/lrs_benchmarking/vcf_called/deepvariant/tmp/HG002_ont_30x_1k.hg38.mm2-pb/intermediate_results_dir/call_variants_output.tfrecord.gz" --examples "/mnt/storage/genetic_data/WGS_LR/Tests/lrs_benchmarking/vcf_called/deepvariant/tmp/HG002_ont_30x_1k.hg38.mm2-pb/intermediate_results_dir/make_examples.tfrecord@64.gz" --checkpoint "/opt/models/ont_r104"


rule deepvariant_gpu_postprocess_variants:
    input:   call_variants_output=expand("{{cwd}}/vcf_called/snv_indel/tmp/{{dataset}}.hg38.{{mapper}}.dv-ont-woPG/call_variants-00{scatter_index}-of-00016.tfrecord.gz", scatter_index=CALL_VARIANTS_SCATTER_INDICES), gvcf=expand("{{cwd}}/vcf_called/snv_indel/tmp/{{dataset}}.hg38.{{mapper}}.dv-ont-woPG/gvcf.tfrecord-00{scatter_index}-of-00064.gz", scatter_index=MAKE_EXAMPLES_SCATTER_INDICES), fasta=FASTA_PATH
    output:  vcf_gz="{cwd}/vcf_called/snv_indel/{dataset}.hg38.{mapper}.dv-ont-woPG.vcf.gz"
    params:  call_variants_output="{cwd}/vcf_called/snv_indel/tmp/{dataset}.hg38.{mapper}.dv-ont-woPG/call_variants.tfrecord.gz", gvcf="{cwd}/vcf_called/snv_indel/tmp/{dataset}.hg38.{mapper}.dv-ont-woPG/gvcf.tfrecord@64.gz", log="{cwd}/vcf_called/snv_indel/{dataset}.hg38.{mapper}.dv-ont-woPG.postprocess_variants.log"
    message: "executing {rule} with output {output} and input {input}"
    priority: 8
    threads: 16
    resources:
	    mem_gb=160
    shell:   "umask 0027; \
		mkdir -p $(dirname {output.vcf_gz}); \
                srun -p all -c {threads} --mem={resources.mem_gb}GB \
                docker run --cpus {threads} -m {resources.mem_gb}g -u $UID:1002 --workdir /tmp --rm -v {CWD}:{CWD} -v {DB_DIR}:{DB_DIR}:ro google/deepvariant:{DEEPVARIANT_VERSION} /bin/bash -c \" \
			printf 'Container ID:\\t'; hostname; \
			printf 'Start time:\\t'; date; \
			umask 0027; \
			time \
			/opt/deepvariant/bin/postprocess_variants \
				--ref {input.fasta} \
				--infile {params.call_variants_output} \
				--outfile {output.vcf_gz} \
				--cpus '{threads}'; \
			[[ \$(bcftools view -H {output.vcf_gz} | wc -l) -le 1 ]] && exit 101 || echo 'File size: OK'; \
			printf 'End time:\\t'; date; \" \
		&> {params.log};"

#/opt/deepvariant/bin/postprocess_variants --ref "/mnt/storage/db/references/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta" --infile "/mnt/storage/genetic_data/WGS_LR/Tests/lrs_benchmarking/vcf_called/deepvariant/tmp/HG002_ont_30x_1k.hg38.mm2-pb/intermediate_results_dir/call_variants_output.tfrecord.gz" --outfile "/mnt/storage/genetic_data/WGS_LR/Tests/lrs_benchmarking/vcf_called/deepvariant/HG002_ont_30x_1k.hg38.mm2-pb.vcf.gz" --cpus "64" --small_model_cvo_records "/mnt/storage/genetic_data/WGS_LR/Tests/lrs_benchmarking/vcf_called/deepvariant/tmp/HG002_ont_30x_1k.hg38.mm2-pb/intermediate_results_dir/make_examples_call_variant_outputs.tfrecord@64.gz" --gvcf_outfile "/mnt/storage/genetic_data/WGS_LR/Tests/lrs_benchmarking/gvcf_called/deepvariant/HG002_ont_30x_1k.hg38.mm2-pb.g.vcf.gz" --nonvariant_site_tfrecord_path "/mnt/storage/genetic_data/WGS_LR/Tests/lrs_benchmarking/vcf_called/deepvariant/tmp/HG002_ont_30x_1k.hg38.mm2-pb/intermediate_results_dir/gvcf.tfrecord@64.gz"


#[[ \$(du -b {output.gvcf_gz} | cut -f 1) -le 16384 ]] && exit 99 || echo 'File size: OK'; \
#[[ \$(du -b {output.vcf_gz} | cut -f 1) -le 16384 ]] && exit 99 || echo 'File size: OK'; \


rule deepvariant_cpu_only:
    input:   cram="{cwd}/cram/{dataset}.hg38.{mapper}.cram", crai="{cwd}/cram/{dataset}.hg38.{mapper}.cram.crai", fasta=FASTA_PATH
    output:  vcf_gz="{cwd}/vcf_called/snv_indel/{dataset}.hg38.{mapper}.dv-ont-woPG.vcf.gz"
    params:  intermediate_results_dir="{cwd}/vcf_called/snv_indel/tmp/{dataset}.hg38.{mapper}.dv-ont-woPG/intermediate_results_dir", logging_dir="{cwd}/vcf_called/snv_indel/tmp/{dataset}.hg38.{mapper}.dv-ont-woPG/logging_dir"
    message: "executing {rule} with output {output} and input {input}"
    threads: 64
    resources:
	    mem_gb=160
    shell:   "umask 0027; \
		mkdir -p {params.intermediate_results_dir}; \
		mkdir -p {params.logging_dir}; \
                srun -p all -c {threads} --mem={resources.mem_gb}GB \
		docker run --cpus {threads} -m {resources.mem_gb}g -u $UID:1002 --workdir /tmp --rm -v {CWD}:{CWD} -v {DB_DIR}:{DB_DIR}:ro google/deepvariant:{DEEPVARIANT_VERSION} /bin/bash -c \" \
			printf 'Container ID:\\t'; hostname; \
			printf 'Start time:\\t'; date; \
			umask 0027; \
			/opt/deepvariant/bin/run_deepvariant \
				--num_shards {threads} \
				--model_type ONT_R104 \
				--reads {input.cram} \
				--ref {input.fasta} \
				--intermediate_results_dir {params.intermediate_results_dir} \
                                --logging_dir {params.logging_dir} \
				--novcf_stats_report \
				--output_vcf {output.vcf_gz}; \
			[[ \$(bcftools view -H {output.vcf_gz} | wc -l) -le 1 ]] && exit 101 || echo 'File size: OK'; \
			printf 'End time:\\t'; date; \" \
		&> {output.vcf_gz}.log;"

#--make_examples_extra_args 'normalize_reads=True' \
#--dry_run
#--only_keep_pass \
#--make_examples_extra_args '' \
#--call_variants_extra_args '' \
#--postprocess_variants_extra_args '' \


rule deepvariant_visual_report:
    input:   vcf_gz="{cwd}/vcf_called/snv_indel/{dataset}.hg38.{mapper}.dv-ont-woPG.vcf.gz"
    output:  visual_report="{cwd}/vcf_called/snv_indel/{dataset}.hg38.{mapper}.dv-ont-woPG.visual_report.html"
    params:  prefix="{cwd}/vcf_called/snv_indel/{dataset}.hg38.{mapper}.dv-ont-woPG", log="{cwd}/vcf_called/snv_indel/{dataset}.hg38.{mapper}.dv-ont-woPG.visual_report.log"
    message: "executing {rule} with output {output} and input {input}"
    priority: 2
    threads: 2
    resources:
	    mem_gb=10
    shell:   "umask 0027; \
		mkdir -p $(dirname {output.visual_report}); \
		srun -p all -c {threads} --mem={resources.mem_gb}GB \
		docker run --cpus {threads} -m {resources.mem_gb}g -u $UID:1002 --workdir /tmp --rm -v {CWD}:{CWD} google/deepvariant:{DEEPVARIANT_VERSION} /bin/bash -c \" \
			printf 'Container ID:\\t'; hostname; \
			printf 'Start time:\\t'; date; \
			umask 0027; \
			time \
			/opt/deepvariant/bin/vcf_stats_report \
				--input_vcf {input.vcf_gz} \
				--outfile_base {params.prefix}; \
			printf 'End time:\\t'; date; \" \
		&> {params.log};"


rule bcftools_stats:
    input:   vcf_gz="{cwd}/vcf_called/snv_indel/{dataset}.hg38.{mapper}.dv-ont-woPG.vcf.gz", fasta=FASTA_PATH
    output:  stats="{cwd}/vcf_called/snv_indel/{dataset}.hg38.{mapper}.dv-ont-woPG.vcf.gz.stats"
    message: "executing {rule} with {output} and input {input}"
    threads: 2
    resources:
	    mem_gb=4
    shell:   "umask 0027; \
		mkdir -p $(dirname {output.stats}); \
                srun -p all -c {threads} --mem={resources.mem_gb}GB \
                docker run --cpus {threads} -m {resources.mem_gb}g -u $UID:1002 --workdir /tmp --rm -v {CWD}:{CWD} -v {DB_DIR}:{DB_DIR}:ro storage-node:5000/own/genetic_data_analysis:{GDA_VERSION} /bin/bash -c \" \
			printf 'Container ID:\\t'; hostname; \
			printf 'Start time:\\t'; date; \
			umask 0027; \
			bcftools stats \
				--threads {threads} \
				-F {input.fasta} \
				-s - \
				{input.vcf_gz} \
				> {output.stats}; \
			printf 'End time:\\t'; date; \" \
                &> {output.stats}.log;"
			

# ************************************************************************************************
