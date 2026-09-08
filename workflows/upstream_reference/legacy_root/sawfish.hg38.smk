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


# Create wildcards
DATASETS,MAPPERS = glob_wildcards(CWD + "/cram/{dataset,[A-Za-z0-9._-]+}.hg38.{mapper,[A-Za-z0-9-]+}.cram")


#Logging
print("Database directory: " + DB_DIR)
print("Path to reference genome: " + FASTA_PATH)

print("Current workding directory: " + CWD)
print("Datasets: " + ','.join(DATASETS))
print("Mappers: " + ','.join(MAPPERS))



# *** Define Output
OUTPUT = [];

OUTPUT = OUTPUT + expand(CWD + "/sawfish/{dataset}.hg38.{mapper}.discover_dir/candidate.sv.bcf", zip, dataset=DATASETS, mapper=MAPPERS)

OUTPUT = OUTPUT + expand(CWD + "/sawfish/{dataset}.hg38.{mapper}.joint_call_dir/genotyped.sv.vcf.gz", zip, dataset=DATASETS, mapper=MAPPERS)

OUTPUT = OUTPUT + expand(CWD + "/sawfish/{dataset}.hg38.{mapper}.sawfish.vcf.gz", zip, dataset=DATASETS, mapper=MAPPERS)

OUTPUT = OUTPUT + expand(CWD + "/truvari/{dataset}.hg38.{mapper}.sawfish/summary.json", zip, dataset=DATASETS, mapper=MAPPERS)


# ************************************************************************************************


rule all:
    input: OUTPUT

rule test:
    shell:  print(OUTPUT) # print(DATASETS),


# ************************************************************************************************	
# Rules
# ************************************************************************************************		

rule sawfish_discover:
    input:   cram="{cwd}/cram/{dataset}.hg38.{mapper}.cram", crai="{cwd}/cram/{dataset}.hg38.{mapper}.cram.crai", fasta=FASTA_PATH, expected_cn="{cwd}/sawfish/expected_cn.hg38.XY.bed", cnv_excluded_regions="{cwd}/sawfish/annotation_and_common_cnv.hg38.bed.gz"
    output:  bcf="{cwd}/sawfish/{dataset}.hg38.{mapper}.discover_dir/candidate.sv.bcf"
    params:  discover_dir="{cwd}/sawfish/{dataset}.hg38.{mapper}.discover_dir"
    log:     "{cwd}/sawfish/{dataset}.hg38.{mapper}.discover_dir.log"
    message: "executing {rule} with output {output} and input {input}"
    threads: 16
    resources:
	    mem_gb=160
    shell:   "umask 0027; \
		mkdir -p {params.discover_dir}; \
		rm -r {params.discover_dir}; \
		srun -p all -c {threads} --mem={resources.mem_gb}GB /bin/bash -c \" \
			printf 'Container ID:\\t'; hostname; \
			printf 'Start time:\\t'; date; \
			umask 0027; \
			/mnt/storage/groups/genetics/VarCAD-dev/bin/mamba run --live-stream -n sawfish sawfish discover \
				--threads {threads} \
				--ref {input.fasta} \
				--bam {input.cram} \
				--expected-cn {input.expected_cn} \
				--cnv-excluded-regions {input.cnv_excluded_regions} \
				--output-dir {params.discover_dir}; \
			printf 'End time:\\t'; date; \" \
		&> {log};"


rule sawfish_joint_call:
    input:   bcf="{cwd}/sawfish/{dataset}.hg38.{mapper}.discover_dir/candidate.sv.bcf"
    output:  vcf_gz="{cwd}/sawfish/{dataset}.hg38.{mapper}.joint_call_dir/genotyped.sv.vcf.gz"
    params:  discover_dir="{cwd}/sawfish/{dataset}.hg38.{mapper}.discover_dir", joint_call_dir="{cwd}/sawfish/{dataset}.hg38.{mapper}.joint_call_dir"
    log:     "{cwd}/sawfish/{dataset}.hg38.{mapper}.joint_call_dir.log"
    message: "executing {rule} with output {output} and input {input}"
    threads: 16
    resources:
            mem_gb=160
    shell:   "umask 0027; \
		mkdir -p {params.joint_call_dir}; \
		rm -r {params.joint_call_dir}; \
		srun -p all -c {threads} --mem={resources.mem_gb}GB /bin/bash -c \" \
			printf 'Container ID:\\t'; hostname; \
			printf 'Start time:\\t'; date; \
			umask 0027; \
			/mnt/storage/groups/genetics/VarCAD-dev/bin/mamba run --live-stream -n sawfish sawfish joint-call \
				--threads {threads} \
				--sample {params.discover_dir} \
				--output-dir {params.joint_call_dir}; \
			printf 'End time:\\t'; date; \" \
		&> {log};"


rule copy_vcf_gz:
    input:   vcf_gz="{cwd}/sawfish/{dataset}.hg38.{mapper}.joint_call_dir/genotyped.sv.vcf.gz"
    output:  vcf_gz="{cwd}/sawfish/{dataset}.hg38.{mapper}.sawfish.vcf.gz"
    params:  
    log:     "{cwd}/sawfish/{dataset}.hg38.{mapper}.sawfish.vcf.gz.log"
    message: "executing {rule} with output {output} and input {input}"
    threads: 2
    resources:
	    mem_gb=4
    shell:   "umask 0027; \
		srun -p all -c {threads} --mem={resources.mem_gb}GB /bin/bash -c \" \
			printf 'Container ID:\\t'; hostname; \
			printf 'Start time:\\t'; date; \
			umask 0027; \
			cp {input.vcf_gz} {output.vcf_gz}; \
			cp {input.vcf_gz}.tbi {output.vcf_gz}.tbi; \
			printf 'End time:\\t'; date; \" \
		&> {log};"


rule truvari_bench_refine:
    input:   comp_vcf_gz="{cwd}/sawfish/{dataset}.hg38.{mapper}.sawfish.vcf.gz", base_vcf_gz="{cwd}/GRCh38_HG2-T2TQ100-V1.1_stvar.svtype.vcf.gz", bed="{cwd}/GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.only_autosomes.bed", fasta=FASTA_PATH
    output:  summary="{cwd}/truvari/{dataset}.hg38.{mapper}.sawfish/summary.json"
    params:  output_dir="{cwd}/truvari/{dataset}.hg38.{mapper}.sawfish"
    log:     "{cwd}/truvari/{dataset}.hg38.{mapper}.sawfish.log"
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
