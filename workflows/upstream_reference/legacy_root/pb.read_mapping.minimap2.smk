# ************************************************************************************************
# 
# ************************************************************************************************

#################
# Include header

include: "header.smk"


#VarCAD_db processing version
VARCAD_DB_VERSION="1.0.0"
print("VarCAD_db procssing version: " + VARCAD_DB_VERSION)


# Create wildcards
DATASETS_FASTQ, = glob_wildcards(CWD + "/fastq/{dataset, [A-Za-z0-9\-\_\.]+}.fastq.gz")
DATASETS_BAM, = glob_wildcards(CWD + "/bam_unmapped/{dataset, [A-Za-z0-9\-\_\.]+}.bam")

DATASETS = DATASETS_FASTQ + DATASETS_BAM


# *** Define Output
OUTPUT = []

OUTPUT = OUTPUT + expand(CWD + "/cram/{dataset}." + VARCAD_GENOME_BUILD + ".mm2-pb.cram", zip, dataset=DATASETS)

#OUTPUT = OUTPUT + expand(CWD + "/cram/{dataset}." + VARCAD_GENOME_BUILD + ".mm2-pb.cram.md5", zip, dataset=DATASETS)

OUTPUT = OUTPUT + expand(CWD + "/cram/{dataset}." + VARCAD_GENOME_BUILD + ".mm2-pb.cram.idxstats", zip, dataset=DATASETS)

#OUTPUT = OUTPUT + expand(CWD + "/cram/{dataset}." + VARCAD_GENOME_BUILD + ".mm2-pb.cram.flagstat", zip, dataset=DATASETS)

OUTPUT = OUTPUT + expand(CWD + "/cram/{dataset}." + VARCAD_GENOME_BUILD + ".mm2-pb.cram.stats", zip, dataset=DATASETS)


# ************************************************************************************************


rule all:
    input: OUTPUT

rule test:
    shell:  print(OUTPUT) # print(DATASETS),


# ************************************************************************************************	
# Rules
# ************************************************************************************************		


ruleorder: minimap2_pb_unmapped_bam > minimap2_pb_fastq


rule minimap2_pb_fastq:
	input:	fastq="{cwd}/fastq/{dataset}.fastq.gz",
		fasta=VARCAD_DB_PATH + "/{genome_build}/Reference_sequence/" + VARCAD_DB_VERSION + "/{genome_build}.fasta.gz"
	output:	cram="{cwd}/cram/{dataset}.{genome_build}.mm2-pb.cram"
	params:	prefix="{cwd}/cram/tmp/{dataset}.{genome_build}.mm2-pb"
	log:	"{cwd}/cram/{dataset}.{genome_build}.mm2-pb.cram.log"
	message:"executing {rule} with output {output} and input {input}"
	threads:64
	resources:
		mem_gb=128
	shell:	"umask 0027; \
		mkdir -p $(dirname {output.cram})/tmp; \
		srun -p all -c {threads} --mem={resources.mem_gb}GB /bin/bash -c \" \
			printf 'Container ID:\\t'; hostname; \
			printf 'Start time:\\t'; date; \
			umask 0027; \
			/mnt/storage/groups/heinz/tools/minimap2-2.28_x64-linux/minimap2 \
				-t {threads} \
				-R '@RG\\tID:{wildcards.dataset}\\tSM:{wildcards.dataset}' \
				-a \
				-y \
				-x map-hifi \
				-L \
				--cs \
				--MD \
				{input.fasta} \
				{input.fastq} | \
			grep -v '^@PG' | \
			\\$VARCAD_PATH/bin/samtools sort \
				-@ {threads} \
				--reference {input.fasta} \
				-m 1G \
				-T {params.prefix} \
				--no-PG \
				-O cram \
				-o {output.cram} \
				-; \
			[[ \$(du -b {output.cram} | cut -f 1) -le 64 ]] && exit 101 || echo 'File size: OK'; \
			\\$VARCAD_PATH/bin/samtools index -@ {threads} {output.cram}; \
			printf 'End time:\\t'; date; \" \
		&> {log};"


rule minimap2_pb_unmapped_bam:
	input:	unmapped_bam="{cwd}/bam_unmapped/{dataset}.bam",
		fasta=VARCAD_DB_PATH + "/{genome_build}/Reference_sequence/" + VARCAD_DB_VERSION + "/{genome_build}.fasta.gz"
	output:	cram="{cwd}/cram/{dataset}.{genome_build}.mm2-pb.cram"
	params:	prefix="{cwd}/cram/tmp/{dataset}.{genome_build}.mm2-pb"
	log:	"{cwd}/cram/{dataset}.{genome_build}.mm2-pb.cram.log"
	message:"executing {rule} with output {output} and input {input}"
	threads:64
	resources:
		mem_gb=128
	shell:	"umask 0027; \
		mkdir -p $(dirname {output.cram})/tmp; \
		srun -p all -c {threads} --mem={resources.mem_gb}GB /bin/bash -c \" \
			printf 'Container ID:\\t'; hostname; \
			printf 'Start time:\\t'; date; \
			umask 0027; \
			\\$VARCAD_PATH/bin/samtools fastq \
				-T ac,ec,ma,ML,MM,np,qe,qs,rq,sn,we,ws,zm \
				{input.unmapped_bam} | \
			/mnt/storage/groups/heinz/tools/minimap2-2.28_x64-linux/minimap2 \
				-t {threads} \
				-R '@RG\\tID:{wildcards.dataset}\\tSM:{wildcards.dataset}' \
				-a \
				-y \
				-x map-hifi \
				-L \
				--cs \
				--MD \
				{input.fasta} \
				- | \
			grep -v '^@PG' | \
			\\$VARCAD_PATH/bin/samtools sort \
				-@ {threads} \
				-m 1G \
				-T {params.prefix} \
				--no-PG \
				--reference {input.fasta} \
				-O cram \
				-o {output.cram} \
				-; \
			[[ \$(du -b {output.cram} | cut -f 1) -le 64 ]] && exit 101 || echo 'File size: OK'; \
			\\$VARCAD_PATH/bin/samtools index -@ {threads} {output.cram}; \
			printf 'End time:\\t'; date; \" \
		&> {log};"


rule md5sum_cram:
	input:	cram="{cwd}/cram/{dataset}.{genome_build}.mm2-pb.cram"
	output:	md5="{cwd}/cram/{dataset}.{genome_build}.mm2-pb.cram.md5"
	log:	"{cwd}/cram/{dataset}.{genome_build}.mm2-pb.cram.md5.log"
	message:"executing {rule} with output {output} and input {input}"
	threads:1
	resources:
		mem_gb=2
	shell:	"umask 0027; \
		mkdir -p $(dirname {output}); \
		srun -p all -c {threads} --mem={resources.mem_gb}GB \
		docker run --cpus {threads} -m {resources.mem_gb}g -u $UID:1002 --workdir /tmp --rm -v {CWD}:{CWD} -v {VARCAD_DB_PATH}:{VARCAD_DB_PATH}:ro {DOCKER_VARCAD} /bin/bash -c \" \
			printf 'Container ID:\\t'; hostname; \
			printf 'Start time:\\t'; date; \
			umask 0027; \
			md5sum {input.cram} | awk '{{print \\$1}}' > {output.md5}; \
			printf 'End time:\\t'; date; \" \
		&> {log};"


rule samtools_idxstats:
	input:	cram="{cwd}/cram/{dataset}.{genome_build}.mm2-pb.cram"
	output:	idxstats="{cwd}/cram/{dataset}.{genome_build}.mm2-pb.cram.idxstats"
	log:	"{cwd}/cram/{dataset}.{genome_build}.mm2-pb.cram.idxstats.log"
	message:"executing {rule} with output {output} and input {input}"
	threads:2
	resources:
		mem_gb=4
	shell:	"umask 0027; \
		srun -p all -c {threads} --mem={resources.mem_gb}GB \
		docker run --cpus {threads} -m {resources.mem_gb}g -u $UID:1002 --workdir /tmp --rm -v {CWD}:{CWD} -v {VARCAD_DB_PATH}:{VARCAD_DB_PATH}:ro {DOCKER_VARCAD} /bin/bash -c \" \
			printf 'Container ID:\\t'; hostname; \
			printf 'Start time:\\t'; date; \
			umask 0027; \
			\\$VARCAD_PATH/bin/samtools idxstats \
				-@ {threads} \
				{input.cram} \
				> {output.idxstats}; \
			[[ \$(du -b {output.idxstats} | cut -f 1) -le 0 ]] && exit 101 || echo 'File size: OK'; \
			printf 'End time:\\t'; date; \" \
		&> {log};"


rule samtools_flagstat:
	input:	cram="{cwd}/cram/{dataset}.{genome_build}.mm2-pb.cram"
	output:	flagstat="{cwd}/cram/{dataset}.{genome_build}.mm2-pb.cram.flagstat"
	log:	"{cwd}/cram/{dataset}.{genome_build}.mm2-pb.cram.flagstat.log"
	message:"executing {rule} with output {output} and input {input}"
	threads:2
	resources:
		mem_gb=4
	shell:	"umask 0027; \
		srun -p all -c {threads} --mem={resources.mem_gb}GB \
		docker run --cpus {threads} -m {resources.mem_gb}g -u $UID:1002 --workdir /tmp --rm -v {CWD}:{CWD} -v {VARCAD_DB_PATH}:{VARCAD_DB_PATH}:ro {DOCKER_VARCAD} /bin/bash -c \" \
			printf 'Container ID:\\t'; hostname; \
			printf 'Start time:\\t'; date; \
			umask 0027; \
			\\$VARCAD_PATH/bin/samtools idxstats \
				-@ {threads} \
				{input.cram} \
				> {output.flagstat}; \
			[[ \$(du -b {output.flagstat} | cut -f 1) -le 0 ]] && exit 101 || echo 'File size: OK'; \
			printf 'End time:\\t'; date; \" \
		&> {log};"


rule samtools_stats:
	input:	cram="{cwd}/cram/{dataset}.{genome_build}.mm2-pb.cram",
		fasta=VARCAD_DB_PATH + "/{genome_build}/Reference_sequence/" + VARCAD_DB_VERSION + "/{genome_build}.fasta.gz"
	output:	stats="{cwd}/cram/{dataset}.{genome_build}.mm2-pb.cram.stats"
	log:	"{cwd}/cram/{dataset}.{genome_build}.mm2-pb.cram.stats.log"
	message:"executing {rule} with output {output} and input {input}"
	threads:16
	resources:
		mem_gb=32
	shell:	"umask 0027; \
		srun -p all -c {threads} --mem={resources.mem_gb}GB /bin/bash -c \" \
			printf 'Container ID:\\t'; hostname; \
			printf 'Start time:\\t'; date; \
			umask 0027; \
			\\$VARCAD_PATH/bin/samtools stats \
				-@ {threads} \
				--reference {input.fasta} \
				--remove-overlaps \
				{input.cram} \
				> {output.stats}; \
			[[ \$(du -b {output.stats} | cut -f 1) -lt 5000 ]] && exit 101 || echo 'File size: OK'; \
			printf 'End time:\\t'; date; \" \
		&> {log};"


# ************************************************************************************************
