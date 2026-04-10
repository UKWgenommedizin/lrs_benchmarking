# ************************************************************************************************
# * Snakefile for SNP calling (test version)
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


# Create wildcards: mapper should be either pbmm2-pb or pbmm2-ont
DATASETS,MAPPERS = glob_wildcards(CWD + "/cram/{dataset,[A-Za-z0-9\-\._]+}.hg38.{mapper,pbmm2-pb|pbmm2-ont}.cram")


#Logging
print("Database directory: " + DB_DIR)
print("Path to reference genome: " + FASTA_PATH)

print("Current workding directory: " + CWD)
print("Datasets: " + ','.join(DATASETS))
print("Mappers: " + ','.join(MAPPERS))



# *** Define Output
OUTPUT = [];

OUTPUT = OUTPUT + expand(CWD + "/pbsv/{dataset}.hg38.{mapper}.svsig.gz", zip, dataset=DATASETS, mapper=MAPPERS)

OUTPUT = OUTPUT + expand(CWD + "/pbsv/{dataset}.hg38.{mapper}.pbsv.vcf.gz", zip, dataset=DATASETS, mapper=MAPPERS)

OUTPUT = OUTPUT + expand(CWD + "/truvari/{dataset}.hg38.{mapper}.pbsv/summary.json", zip, dataset=DATASETS, mapper=MAPPERS)


# ************************************************************************************************


rule all:
    input: OUTPUT

rule test:
    shell:  print(OUTPUT) # print(DATASETS),


# ************************************************************************************************	
# Rules
# ************************************************************************************************		

rule pbsv_discover:
    input:   cram="{cwd}/cram/{dataset}.hg38.{mapper}.cram", crai="{cwd}/cram/{dataset}.hg38.{mapper}.cram.crai", fasta=FASTA_PATH
    output:  svsig="{cwd}/pbsv/{dataset}.hg38.{mapper}.svsig.gz"
    params:  tandem_repeats="{cwd}/pbsv/human_GRCh38_no_alt_analysis_set.trf.bed",
             tmp_bam="{cwd}/pbsv/{dataset}.hg38.{mapper}.tmp.bam"
    log:     "{cwd}/pbsv/{dataset}.hg38.{mapper}.svsig.gz.log"
    message: "executing {rule} with output {output} and input {input}"
    threads: 8
    resources:
	    mem_gb=32
    shell:   "umask 0027; \
		mkdir -p $(dirname {output.svsig}); \
		srun -p all -c {threads} --mem={resources.mem_gb}GB /bin/bash -c \" \
			printf 'Container ID:\\t'; hostname; \
			printf 'Start time:\\t'; date; \
			umask 0027; \
			samtools view -b -@ {threads} -T {input.fasta} -o {params.tmp_bam} {input.cram}; \
			samtools index {params.tmp_bam}; \
			/mnt/storage/groups/genetics/VarCAD-dev/external/Miniforge3/condabin/mamba run --live-stream -n pbsv pbsv discover \
				--tandem-repeats {params.tandem_repeats} \
				{params.tmp_bam} \
				{output.svsig}; \
			rm -f {params.tmp_bam} {params.tmp_bam}.bai; \
			printf 'End time:\\t'; date; \" \
		&> {log};"


rule pbsv_call:
    input:   svsig="{cwd}/pbsv/{dataset}.hg38.{mapper}.svsig.gz", fasta=FASTA_PATH
    output:  vcf="{cwd}/pbsv/{dataset}.hg38.{mapper}.pbsv.vcf"
    params:  
    log:     "{cwd}/pbsv/{dataset}.hg38.{mapper}.pbsv.vcf.log"
    message: "executing {rule} with output {output} and input {input}"
    threads: 8
    resources:
            mem_gb=32
    shell:   "umask 0027; \
		mkdir -p $(dirname {output.vcf}); \
		srun -p all -c {threads} --mem={resources.mem_gb}GB /bin/bash -c \" \
			printf 'Container ID:\\t'; hostname; \
			printf 'Start time:\\t'; date; \
			umask 0027; \
			/mnt/storage/groups/genetics/VarCAD-dev/external/Miniforge3/condabin/mamba run --live-stream -n pbsv pbsv call \
				--num-threads {threads} \
				{input.fasta} \
				{input.svsig} \
				{output.vcf}; \
			printf 'End time:\\t'; date; \" \
		&> {log};"


rule bgzip_index_vcf:
    input:   vcf="{cwd}/pbsv/{dataset}.hg38.{mapper}.pbsv.vcf"
    output:  vcf_gz="{cwd}/pbsv/{dataset}.hg38.{mapper}.pbsv.vcf.gz"
    params:  
    log:     "{cwd}/pbsv/{dataset}.hg38.{mapper}.pbsv.vcf.gz.log"
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
    input:   comp_vcf_gz="{cwd}/pbsv/{dataset}.hg38.{mapper}.pbsv.vcf.gz", base_vcf_gz="{cwd}/GRCh38_HG2-T2TQ100-V1.1_stvar.svtype.vcf.gz", bed="{cwd}/GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.only_autosomes.bed", fasta=FASTA_PATH
    output:  summary="{cwd}/truvari/{dataset}.hg38.{mapper}.pbsv/summary.json"
    params:  output_dir="{cwd}/truvari/{dataset}.hg38.{mapper}.pbsv"
    log:     "{cwd}/truvari/{dataset}.hg38.{mapper}.pbsv.log"
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


rule bcftools_stats:
    input:   vcf_gz="{cwd}/vcf_called/deepvariant/{dataset}.hg38.vcf.gz", fasta=FASTA_PATH
    output:  stats="{cwd}/vcf_called/deepvariant/{dataset}.hg38.vcf.gz.stats"
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
