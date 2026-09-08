# ************************************************************************************************
# * Snakefile for SV calling with cuteSV
# ************************************************************************************************

import os
import math


#Define variables
DB_DIR = "/mnt/storage/db"
FASTA = "references/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta"


#Set FASTA_PATH
FASTA_PATH = DB_DIR + "/" + FASTA


#Get current working directory
CWD = os.getcwd()


# Create wildcards
DATASETS,MAPPERS = glob_wildcards(CWD + "/cram/{dataset,[A-Za-z0-9\-\._]+}.hg38.{mapper,[A-Za-z0-9\-]+}.cram")


#Logging
print("Database directory: " + DB_DIR)
print("Path to reference genome: " + FASTA_PATH)

print("Current workding directory: " + CWD)
print("Datasets: " + ','.join(DATASETS))
print("Mappers: " + ','.join(MAPPERS))



# *** Define Output
OUTPUT = [];

OUTPUT = OUTPUT + expand(CWD + "/cuteSV/{dataset}.hg38.{mapper}.cuteSV.vcf.gz", zip, dataset=DATASETS, mapper=MAPPERS)

OUTPUT = OUTPUT + expand(CWD + "/truvari/{dataset}.hg38.{mapper}.cuteSV/summary.json", zip, dataset=DATASETS, mapper=MAPPERS)


# ************************************************************************************************


rule all:
    input: OUTPUT

rule test:
    shell:  print(OUTPUT) # print(DATASETS),


# ************************************************************************************************	
# Rules
# ************************************************************************************************		

rule cuteSV_call:
    input:   cram="{cwd}/cram/{dataset}.hg38.{mapper}.cram", crai="{cwd}/cram/{dataset}.hg38.{mapper}.cram.crai", fasta=FASTA_PATH
    output:  vcf="{cwd}/cuteSV/{dataset}.hg38.{mapper}.cuteSV.vcf"
    params:  work_dir="{cwd}/cuteSV/{dataset}.hg38.{mapper}.work"
    log:     "{cwd}/cuteSV/{dataset}.hg38.{mapper}.cuteSV.vcf.log"
    message: "executing {rule} with output {output} and input {input}"
    threads: 16
    resources:
	    mem_gb=32
    shell:   "umask 0027; \
		mkdir -p {params.work_dir}; \
		srun -p all -c {threads} --mem={resources.mem_gb}GB /bin/bash -c \" \
			printf 'Container ID:\\t'; hostname; \
			printf 'Start time:\\t'; date; \
			umask 0027; \
			/mnt/storage/groups/genetics/VarCAD-dev/external/Miniforge3/condabin/mamba run --live-stream -n cutesv cuteSV \
				--threads {threads} \
				--sample {wildcards.dataset} \
				--genotype \
				{input.cram} \
				{input.fasta} \
				{output.vcf} \
				{params.work_dir}; \
			printf 'End time:\\t'; date; \" \
		&> {log};"


rule bgzip_index_vcf:
    input:   vcf="{cwd}/cuteSV/{dataset}.hg38.{mapper}.cuteSV.vcf"
    output:  vcf_gz="{cwd}/cuteSV/{dataset}.hg38.{mapper}.cuteSV.vcf.gz"
    log:     "{cwd}/cuteSV/{dataset}.hg38.{mapper}.cuteSV.vcf.gz.log"
    message: "executing {rule} with output {output} and input {input}"
    threads: 4
    resources:
	    mem_gb=8
    shell:   "umask 0027; \
		srun -p all -c {threads} --mem={resources.mem_gb}GB /bin/bash -c \" \
			printf 'Container ID:\\t'; hostname; \
			printf 'Start time:\\t'; date; \
			umask 0027; \
			bgzip -c -@ {threads} {input.vcf} > {output.vcf_gz}; \
			tabix -p vcf {output.vcf_gz}; \
			printf 'End time:\\t'; date; \" \
		&> {log};"


rule truvari_bench_refine:
    input:   comp_vcf_gz="{cwd}/cuteSV/{dataset}.hg38.{mapper}.cuteSV.vcf.gz", base_vcf_gz="{cwd}/GRCh38_HG2-T2TQ100-V1.1_stvar.svtype.vcf.gz", bed="{cwd}/GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.only_autosomes.bed", fasta=FASTA_PATH
    output:  summary="{cwd}/truvari/{dataset}.hg38.{mapper}.cuteSV/summary.json"
    params:  output_dir="{cwd}/truvari/{dataset}.hg38.{mapper}.cuteSV"
    log:     "{cwd}/truvari/{dataset}.hg38.{mapper}.cuteSV.log"
    message: "executing {rule} with output {output} and input {input}"
    threads: 4
    resources:
	    mem_gb=8
    shell:   "umask 0027; \
		mkdir -p $(dirname {params.output_dir}); \
		rm -rf {params.output_dir}; \
		srun -p all -c {threads} --mem={resources.mem_gb}GB /bin/bash -c \" \
			printf 'Container ID:\\t'; hostname; \
			printf 'Start time:\\t'; date; \
			umask 0027; \
			/mnt/storage/groups/genetics/VarCAD-dev/external/Miniforge3/condabin/mamba run --live-stream -n truvari truvari bench \
				--reference {input.fasta} \
				--includebed {input.bed} \
				--base {input.base_vcf_gz} \
				--comp {input.comp_vcf_gz} \
				--output {params.output_dir} \
				--passonly \
				--pick ac \
				--dup-to-ins \
				--refine; \
			printf 'End time:\\t'; date; \" \
		&> {log};"


# ************************************************************************************************
