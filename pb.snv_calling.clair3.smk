# ************************************************************************************************
# * Snakefile for SNP calling (test version)
# ************************************************************************************************

import os
import math


#Define variables
GDA_VERSION = "3.0"
CLAIR3_VERSION = "v1.1.2"
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
print("Version of Docker image 'clair3': " + CLAIR3_VERSION)
print("Database directory: " + DB_DIR)
print("Path to reference genome: " + FASTA_PATH)

print("Current workding directory: " + CWD)
print("Datasets: " + ','.join(DATASETS))
print("Mappers: " + ','.join(MAPPERS))


# *** Define Output
OUTPUT = [];

OUTPUT = OUTPUT + expand(CWD + "/vcf_called/snv_indel/tmp/{dataset}.hg38.{mapper}.clair3-revio/merge_output.vcf.gz", zip, dataset=DATASETS, mapper=MAPPERS)

OUTPUT = OUTPUT + expand(CWD + "/vcf_called/snv_indel/{dataset}.hg38.{mapper}.clair3-revio.vcf.gz", zip, dataset=DATASETS, mapper=MAPPERS)

OUTPUT = OUTPUT + expand(CWD + "/vcf_called/snv_indel/{dataset}.hg38.{mapper}.clair3-revio.vcf.gz.stats", zip, dataset=DATASETS, mapper=MAPPERS)


# ************************************************************************************************


rule all:
    input: OUTPUT

rule test:
    shell:  print(OUTPUT) # print(DATASETS),


# ************************************************************************************************	
# Rules
# ************************************************************************************************		


rule run_clair3_revio:
    input:   cram="{cwd}/cram/{dataset}.hg38.{mapper}.cram", crai="{cwd}/cram/{dataset}.hg38.{mapper}.cram.crai", fasta=FASTA_PATH
    output:  vcf_gz="{cwd}/vcf_called/snv_indel/tmp/{dataset}.hg38.{mapper}.clair3-revio/merge_output.vcf.gz"
    params:  output_dir="{cwd}/vcf_called/snv_indel/tmp/{dataset}.hg38.{mapper}.clair3-revio"
    message: "executing {rule} with output {output} and input {input}"
    log:     "{cwd}/vcf_called/snv_indel/{dataset}.hg38.{mapper}.clair3-revio.log"
    threads: 64
    resources:
	    mem_gb=160
    shell:   "umask 0027; \
		mkdir -p $(dirname {output.vcf_gz}); \
                srun -p all -c {threads} --mem={resources.mem_gb}GB \
		docker run --cpus {threads} -m {resources.mem_gb}g -u $UID:1002 --workdir /tmp --rm -v {CWD}:{CWD} -v {DB_DIR}:{DB_DIR}:ro hkubal/clair3:{CLAIR3_VERSION} /bin/bash -c \" \
			printf 'Container ID:\\t'; hostname; \
			printf 'Start time:\\t'; date; \
			umask 0027; \
			/opt/bin/run_clair3.sh \
				--bam_fn={input.cram} \
				--ref_fn={input.fasta} \
				--threads={threads} \
				--platform=hifi \
				--model_path=/opt/models/hifi_revio \
				--output={params.output_dir}; \
			[[ \$(zgrep -v '^#' {output.vcf_gz} | wc -l) -le 1 ]] && exit 101 || echo 'File size: OK'; \
			printf 'End time:\\t'; date; \" \
		&> {log};"


rule copy_clair3_vcf_gz:
    input:   vcf_gz="{cwd}/vcf_called/snv_indel/tmp/{dataset}.hg38.{mapper}.clair3-revio/merge_output.vcf.gz"
    output:  vcf_gz="{cwd}/vcf_called/snv_indel/{dataset}.hg38.{mapper}.clair3-revio.vcf.gz"
    message: "executing {rule} with {output} and input {input}"
    log:     "{cwd}/vcf_called/snv_indel/{dataset}.hg38.{mapper}.clair3-revio.vcf.gz.log"
    threads: 2
    resources:
	    mem_gb=4
    shell:   "umask 0027; \
		mkdir -p $(dirname {output.vcf_gz}); \
		srun -p all -c {threads} --mem={resources.mem_gb}GB \
		docker run --cpus {threads} -m {resources.mem_gb}g -u $UID:1002 --workdir /tmp --rm -v {CWD}:{CWD} -v {DB_DIR}:{DB_DIR}:ro storage-node:5000/own/genetic_data_analysis:{GDA_VERSION} /bin/bash -c \" \
			printf 'Container ID:\\t'; hostname; \
			printf 'Start time:\\t'; date; \
			umask 0027; \
			cp {input.vcf_gz} {output.vcf_gz}; \
			[[ \$(bcftools view -H {output.vcf_gz} | wc -l) -le 1 ]] && exit 101 || echo 'File size: OK'; \
			tabix {output.vcf_gz}; \
			printf 'End time:\\t'; date; \" \
		&> {log};"


rule bcftools_stats:
    input:   vcf_gz="{cwd}/vcf_called/snv_indel/{dataset}.hg38.{mapper}.clair3-revio.vcf.gz", fasta=FASTA_PATH
    output:  stats="{cwd}/vcf_called/snv_indel/{dataset}.hg38.{mapper}.clair3-revio.vcf.gz.stats"
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
