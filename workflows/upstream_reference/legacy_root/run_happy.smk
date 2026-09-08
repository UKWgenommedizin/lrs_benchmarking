# ************************************************************************************************
# 
# ************************************************************************************************


import os
import math


DB_DIR = "/mnt/storage/db"
FASTA = "references/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta"
HAPPY_VERSION = "0.3.15"


#Set FASTA_PATH
FASTA_PATH = DB_DIR + "/" + FASTA


#Get current working directory
CWD = os.getcwd()


# Create wildcards
DATASETS,MAPPERS,CALLERS = glob_wildcards(CWD + "/vcf_called/snv_indel/{dataset, [A-Za-z0-9\\_\\-\\.]+}.hg38.{mapper, [A-Za-z0-9\\-]+}.{caller, [A-Za-z0-9\\-]+}.vcf.gz")


print(*DATASETS)
print(*MAPPERS)
print(*CALLERS)



# *** Define Output
OUTPUT = [];

OUTPUT = OUTPUT + expand(CWD + "/happy_results/{dataset}.hg38.{mapper}.{caller}.happy.sh", zip, dataset=DATASETS, mapper=MAPPERS, caller=CALLERS)

OUTPUT = OUTPUT + expand(CWD + "/happy_results/{dataset}.hg38.{mapper}.{caller}.roc.all.csv.gz", zip, dataset=DATASETS, mapper=MAPPERS, caller=CALLERS)

OUTPUT.append(CWD + "/happy_results/reppy.sh")

OUTPUT.append(CWD + "/happy_results/reppy.html")


# ************************************************************************************************

rule all:
    input: OUTPUT

rule test:
    shell:  print(OUTPUT) # print(DATASETS),

# ************************************************************************************************	
# Rules
# ************************************************************************************************		


#Create hap.py script
rule create_happy_scripts:
    input:   vcf_gz="{cwd}/vcf_called/snv_indel/{dataset}.hg38.{mapper}.{caller}.vcf.gz", truth_set="{cwd}/happy_data/NA24385_GRCh38_1_22_v4.2.1_benchmark.vcf.gz", highconf_bed="{cwd}/happy_data/NA24385_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed", stratification="/mnt/storage/db/GIAB/hg38/genome-stratifications/v3.1/v3.1-GRCh38-v4.2.1-stratifications.tsv", fasta="/mnt/storage/db/deepvariant/GRCh38_hs38d1.fa", sdf="/mnt/storage/db/deepvariant/GRCh38_hs38d1.fa.sdf"
    output:  script="{cwd}/happy_results/{dataset}.hg38.{mapper}.{caller}.happy.sh"
    params:  basename="{cwd}/happy_results/{dataset}.hg38.{mapper}.{caller}"
    message: "executing {rule} following output {output} and input {input}"
    threads: 8
    shell:   "mkdir -p $(dirname {output.script});\
		touch {output.script}; \
                rm {output.script}; \
		printf \" \
		docker run --cpus {threads} -m 128g -u $UID:$GROUPS --rm -v {CWD}:{CWD} -v {DB_DIR}:{DB_DIR}:ro storage-node:5000/illumina/hap.py:{HAPPY_VERSION} \
			/opt/hap.py/bin/hap.py \
                        {input.truth_set} \
                        {input.vcf_gz} \
                        -r {input.fasta} \
			-f {input.highconf_bed} \
                        --stratification {input.stratification} \
                        --engine=vcfeval \
                        --engine-vcfeval-template {input.sdf} \
			--pass-only \
			--roc QUAL \
			--threads {threads} \
			-o {params.basename} \" \
		> {output.script};"

#-l chr20 \


#Run hap.py script
rule run_happy_script:
    input:   script="{cwd}/happy_results/{dataset}.hg38.{mapper}.{caller}.happy.sh", vcf_gz="{cwd}/vcf_called/snv_indel/{dataset}.hg38.{mapper}.{caller}.vcf.gz", truth_set="{cwd}/happy_data/NA24385_GRCh38_1_22_v4.2.1_benchmark.vcf.gz", highconf_bed="{cwd}/happy_data/NA24385_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed", stratification="/mnt/storage/db/GIAB/hg38/genome-stratifications/v3.1/v3.1-GRCh38-v4.2.1-stratifications.tsv", fasta="/mnt/storage/db/deepvariant/GRCh38_hs38d1.fa", sdf="/mnt/storage/db/deepvariant/GRCh38_hs38d1.fa.sdf"
    output:  roc_all="{cwd}/happy_results/{dataset}.hg38.{mapper}.{caller}.roc.all.csv.gz"
    params:  
    log:     "{cwd}/happy_results/{dataset}.hg38.{mapper}.{caller}.happy.log"
    message: "executing {rule} following output {output} and input {input}"
    threads: 8
    shell:   "mkdir -p $(dirname {output.roc_all});\
		touch {output.roc_all}; \
		rm {output.roc_all}; \
                srun -c {threads} --mem=128GB /bin/bash -c \" \
			printf 'Container ID:\\t'; hostname; \
			printf 'Start time:\\t'; date; \
			umask 0027; \
			/bin/bash {input.script}; \
			printf 'End time:\\t'; date; \" \
		&> {log};"


#Create rep.py script
rule create_reppy_script:
    input:   roc_all=expand("{{cwd}}/happy_results/{dataset}.hg38.{mapper}.{caller}.roc.all.csv.gz", zip, dataset=DATASETS, mapper=MAPPERS, caller=CALLERS)
    output:  script="{cwd}/happy_results/reppy.sh"
    params:  html="{cwd}/happy_results/reppy.html", input=list(set(expand("{dataset}.{mapper}.{caller}:{{cwd}}/happy_results/{dataset}.hg38.{mapper}.{caller}.roc.all.csv.gz", zip, dataset=DATASETS, mapper=MAPPERS, caller=CALLERS)))
    message: "executing {rule} following output {output} and input {input}"
    threads: 1
    shell:   "touch {output.script}; \
		rm {output.script}; \
		umask 0027; \
		printf '/opt/benchmarking-tools-master/reporting/basic/bin/rep.py \
			{params.input} \
			--roc-max-datapoints 1000 --roc-resolution 0.00001 --min-recall 0.5 --min-precision 0.5 \
			-o {params.html};' \
		> {output.script};"


#sed -i 's&_ont_30x&-ont-30x&g' {output.script}; \
#sed -i 's&/HG002-ont-30x&/HG002_ont_30x&g' {output.script}; \
#sed -i 's&_pb_30x&-pb-30x&g' {output.script}; \
#sed -i 's&/HG002-pb-30x&/HG002_pb_30x&g' {output.script};"


#Run rep.py script
rule run_reppy_script:
    input:   script="{cwd}/happy_results/reppy.sh"
    output:  html="{cwd}/happy_results/reppy.html"
    params:  
    log:     "{cwd}/happy_results/reppy.log"
    message: "executing {rule} following output {output} and input {input}"
    threads: 1
    shell:   "mkdir -p $(dirname {output.html}); \
		    touch {output.html}; \
		    rm {output.html}; \
		    srun -c {threads} --mem=32GB \
		    docker run --cpus {threads} -m 32g -u $UID:$GROUPS --rm -v {CWD}:{CWD} -v {DB_DIR}:{DB_DIR}:ro storage-node:5000/illumina/hap.py:{HAPPY_VERSION} /bin/bash -c \" \
		    	printf 'Container ID:\\t'; hostname; \
			printf 'Start time:\\t'; date; \
			umask 0027; \
			/bin/bash {input.script}; \
			printf 'End time:\\t'; date; \" \
		&> {log};"


# ************************************************************************************************
