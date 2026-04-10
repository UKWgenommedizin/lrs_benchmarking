##

import os
import math


###################
# VarCAD variables

VARCAD_PATH = os.path.expandvars("$VARCAD_PATH")
print("VarCAD directory: " + VARCAD_PATH)


configfile: VARCAD_PATH + "/config/VarCAD.yaml"


VARCAD_VERSION = config['Version']
print("VarCAD version: " + VARCAD_VERSION)

VARCAD_DOCKER = config['Docker']
print("VarCAD docker: " + str(VARCAD_DOCKER))

VARCAD_DB_PATH = os.path.expandvars(config['VarCAD_db'])
print("VarCAD_db directory: " + VARCAD_DB_PATH)

VARCAD_GENOME_BUILD = config['Genome_build']
print("VarCAD genome build: " + VARCAD_GENOME_BUILD)


############################
# Current working directory

CWD = os.getcwd()
print("Current workding directory: " + CWD)


######################
# Illumina RunFolders

ILLUMINA_RUN_FOLDERS_PATH = "/mnt/storage/genetic_data/WGS/RunFolders"
print("Path to the Illumina run folders: " + ILLUMINA_RUN_FOLDERS_PATH)


################
# Docker images

DOCKER_VARCAD = "storage-node:5000/own/varcad:" + VARCAD_VERSION
print("Docker image for VarCAD: " + DOCKER_VARCAD)

DOCKER_BCL2FASTQ2 = "nfcore/bcl2fastq:2.20.0.422"
print("Docker image for bcl2fastq2: " + DOCKER_BCL2FASTQ2)

DOCKER_VEP = "ensemblorg/ensembl-vep:release_112.0"
print("Docker image for VEP: " + DOCKER_VEP)


#########################
# Unproccessed databases 


#Reference sequence
DOWNLOADS = {}

#hg38
if VARCAD_GENOME_BUILD == "hg38":

	DOWNLOADS[VARCAD_GENOME_BUILD + '/Reference_sequence/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz'] = "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz"

	VEP_DB_VERSION = "112"
	print("VEP data version: " + VEP_DB_VERSION)


#hs1
if VARCAD_GENOME_BUILD == "hs1":

	DOWNLOADS[VARCAD_GENOME_BUILD + '/Reference_sequence/hs1.fa.gz'] = "https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz"


print("Downloads:")
print("File path : URL")
for k,v in DOWNLOADS.items():
	print(k + " : " + v)


