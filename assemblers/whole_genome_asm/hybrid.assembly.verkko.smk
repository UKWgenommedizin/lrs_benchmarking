# ************************************************************************************************
# Verkko hybrid whole-genome assembly
# ************************************************************************************************

#################
# Include shared assembler header
include: "../../header_assembler.smk"

#################
# Verkko version
VERKKO_VERSION = "2.3.2"
DOCKER_VERKKO = "nicolasardila1/lrs-verkko2:" + VERKKO_VERSION

print("Verkko version: " + VERKKO_VERSION)

#####################
# Discover samples and create wildcards

# Discover whole-genome samples with 30x ONT and PacBio HiFi reads.
# Verkko is run only for samples for which both technologies are available.

ONT_SAMPLES, = glob_wildcards(CWD + r"/fastq/{sample,[A-Za-z0-9_-]+}.ont.30x.fastq.gz")

PB_SAMPLES, = glob_wildcards(CWD + r"/fastq/{sample,[A-Za-z0-9_-]+}.pb.30x.fastq.gz")

ONT_SAMPLE_SET = set(ONT_SAMPLES)
PB_SAMPLE_SET = set(PB_SAMPLES)

SAMPLES = sorted(ONT_SAMPLE_SET & PB_SAMPLE_SET)

######################
#Input samples and unpaired control checkpoint

MISSING_PB = sorted(ONT_SAMPLE_SET - PB_SAMPLE_SET)
MISSING_ONT = sorted(PB_SAMPLE_SET - ONT_SAMPLE_SET)

if MISSING_PB:
    raise ValueError("Missing PacBio HiFi input for samples: " + ", ".join(MISSING_PB))

if MISSING_ONT:
    raise ValueError("Missing ONT input for samples: " + ", ".join(MISSING_ONT))

if not SAMPLES:
    raise ValueError("No paired Verkko WGS inputs were found. " "Expected files such as "
        f"{CWD}/fastq/HG002.pb.30x.fastq.gz and "
        f"{CWD}/fastq/HG002.ont.30x.fastq.gz")

##############
# Targets

OUTPUT = []

OUTPUT += expand(CWD + "/assemblies/verkko/{sample}/assembly.fasta", sample=SAMPLES)

rule all:
    input:
        OUTPUT

print("Discover samples and create wildcards")
print(OUTPUT)


################
# Prevent local test sample from being requested explicitly

wildcard_constraints:
    sample = r"[A-Za-z0-9_-]+"


################
# Verkko resource requirements

# The server has approximately 128 CPU threads and 256 GB RAM total.
# The existing server configuration assigns 32 threads and 72 GB RAM
# to the Verkko workflow.
# Verkko local jobs are limited to 64 GB memory.
def get_verkko_memory(wildcards):
    return 72000

################
# Rules

rule verkko_assemble:
    input:
        hifi = CWD + "/fastq/{sample}.pb.30x.fastq.gz",
        ont  = CWD + "/fastq/{sample}.ont.30x.fastq.gz"

    output:
        fasta = CWD + "/assemblies/verkko/{sample}/assembly.fasta"

    params:
        outdir = CWD + "/assemblies/verkko/{sample}/work",
        local_memory_gb = 64

    log:
        CWD + "/assemblies/verkko/{sample}/verkko.log"

    message:
        "executing {rule} with output {output} and input {input}"

    # The existing server configuration assigns 32 threads to Verkko.
    threads: 32

    resources:
        mem_mb = get_verkko_memory

    shell:
        """
        mkdir -p "$(dirname "{log}")"

        (
            set -euo pipefail

            echo "[$(date -Is)] START verkko_assemble {wildcards.sample}"
            echo "Sample: {wildcards.sample}"
            echo "Read technologies: PacBio HiFi + ONT"
            echo "HiFi input: {input.hifi}"
            echo "ONT input: {input.ont}"
            echo "Threads: {threads}"
            echo "Memory: {resources.mem_mb} MB"
            echo "Container hostname: verkko-{wildcards.sample}"

            mkdir -p "{params.outdir}"

            HIFI_REAL="$(readlink -f "{input.hifi}")"
            ONT_REAL="$(readlink -f "{input.ont}")"

            echo "Resolved HiFi input: $HIFI_REAL"
            echo "Resolved ONT input: $ONT_REAL"

            [[ -r "$HIFI_REAL" ]] || {{
                echo "ERROR: Resolved HiFi input is not readable: $HIFI_REAL"
                exit 101
            }}

            [[ -r "$ONT_REAL" ]] || {{
                echo "ERROR: Resolved ONT input is not readable: $ONT_REAL"
                exit 101
            }}

            docker run --rm \
                --hostname verkko-{wildcards.sample} \
                --cpus {threads} \
                -m {resources.mem_mb}m \
                --user "$(id -u):$(id -g)" \
                -e HOME=/tmp \
                -e TMPDIR=/tmp \
                --workdir {CWD} \
                -v {CWD}:{CWD} \
                -v "$HIFI_REAL:/tmp/verkko_hifi.fastq.gz:ro" \
                -v "$ONT_REAL:/tmp/verkko_ont.fastq.gz:ro" \
                {DOCKER_VERKKO} \
                verkko \
                    -d "{params.outdir}" \
                    --hifi "/tmp/verkko_hifi.fastq.gz" \
                    --nano "/tmp/verkko_ont.fastq.gz" \
                    --local \
                    --local-cpus {threads} \
                    --local-memory {params.local_memory_gb}

            [[ -s "{params.outdir}/assembly.fasta" ]] || {{
                echo "ERROR: Verkko assembly.fasta is missing or empty"
                echo "Expected: {params.outdir}/assembly.fasta"
                exit 101;
            }}

            cp "{params.outdir}/assembly.fasta" "{output.fasta}"

            [[ -s "{output.fasta}" ]] || {{
                echo "ERROR: Verkko final assembly.fasta is missing or empty"
                exit 101;
            }}

            grep -q '^>' "{output.fasta}" || {{
                echo "ERROR: Verkko output does not appear to be a valid FASTA"
                exit 101;
            }}

            echo "[$(date -Is)] END verkko_assemble {wildcards.sample}"

        ) > "{log}" 2>&1
        """