#!/usr/bin/env bash

# ============================================================
# Turrialba Fresh Cheese Metagenomics
# Example local configuration
# ============================================================
#
# Copy this file before running scripts that source the
# repository configuration:
#
#   cp config/config.example.sh config/config.sh
#
# Then edit paths according to the local HPC environment.
#
# config/config.sh is intentionally excluded from Git.

# ------------------------------------------------------------
# Repository location
# ------------------------------------------------------------

CONFIG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_HOME="${PROJECT_HOME:-$(cd "${CONFIG_DIR}/.." && pwd)}"

# ------------------------------------------------------------
# Input and working-data locations
# ------------------------------------------------------------

# Directory containing the original paired-end FASTQ files.
RAW_DIR="${RAW_DIR:-/path/to/raw_fastq}"

# Main high-capacity working directory for intermediate and
# final metagenomic results.
SCRATCH_DIR="${SCRATCH_DIR:-/path/to/Shotgun_MAGs_Turrialba}"

# Alias used by later workflow modules.
PROJECT_DATA_ROOT="${PROJECT_DATA_ROOT:-${SCRATCH_DIR}}"

# ------------------------------------------------------------
# Metadata and logs
# ------------------------------------------------------------

SAMPLES="${SAMPLES:-${PROJECT_HOME}/metadata/samples.tsv}"
LOG_DIR="${LOG_DIR:-${PROJECT_HOME}/logs}"

# ------------------------------------------------------------
# Container infrastructure
# ------------------------------------------------------------
#
# The default below corresponds to the UCR HPC installation
# used for the thesis. Change it on other systems.

CTN="${CTN:-/opt/ohpc/pub/containers/BIO}"

FASTQC="${FASTQC:-${CTN}/fastqc-0.12.1.sif}"
FASTP="${FASTP:-${CTN}/fastp-0.23.3.sif}"
MULTIQC="${MULTIQC:-${CTN}/multiqc-1.25.1.sif}"

METAWRAP="${METAWRAP:-${CTN}/metawrap-1.2.sif}"
MEGAHIT="${MEGAHIT:-${CTN}/megahit-1.2.9.sif}"
SPADES="${SPADES:-${CTN}/spades-4.0.0.sif}"
METABAT2="${METABAT2:-${CTN}/metabat2-2.17.sif}"

CHECKM2="${CHECKM2:-${SCRATCH_DIR}/containers/checkm2-1.1.0-r1.sif}"
CHECKM2_DB="${CHECKM2_DB:-${SCRATCH_DIR}/databases/checkm2_1.1.0/CheckM2_database/uniref100.KO.1.dmnd}"

GUNC="${GUNC:-${CTN}/gunc-1.0.6.sif}"
GUNC_DB="${GUNC_DB:-${SCRATCH_DIR}/databases/gunc/gunc_db_progenomes2.1.dmnd}"

MAXBIN2="${MAXBIN2:-${SCRATCH_DIR}/containers/binners/maxbin2-2.2.7.sif}"
CONCOCT="${CONCOCT:-${SCRATCH_DIR}/containers/binners/concoct-1.1.0.sif}"
DASTOOL="${DASTOOL:-${SCRATCH_DIR}/containers/refinement/das_tool-1.1.7-r44.sif}"

MDMCLEANER="${MDMCLEANER:-${SCRATCH_DIR}/containers/mdmcleaner-0.8.7.sif}"
MDMCLEANER_DB="${MDMCLEANER_DB:-${SCRATCH_DIR}/databases/mdmcleaner}"
MDMCLEANER_CONFIG="${MDMCLEANER_CONFIG:-${PROJECT_HOME}/config/mdmcleaner/mdmcleaner.config}"

GTDBTK="${GTDBTK:-${CTN}/gtdbtk-2.4.0.sif}"
DREP="${DREP:-${CTN}/drep-3.5.0.sif}"

BOWTIE2="${BOWTIE2:-${CTN}/bowtie2-2.5.4.sif}"
SAMTOOLS="${SAMTOOLS:-${CTN}/samtools-1.21.sif}"

ABRICATE="${ABRICATE:-${CTN}/abricate-1.0.1.sif}"
EGGNOG="${EGGNOG:-${CTN}/eggnog-mapper-2.1.12.sif}"

# ------------------------------------------------------------
# Sequencing parameters
# ------------------------------------------------------------

READ_LENGTH="${READ_LENGTH:-151}"
FASTQ_EXT="${FASTQ_EXT:-fastq}"

# ------------------------------------------------------------
# Resource limits used during the thesis
# ------------------------------------------------------------

MAX_CPUS="${MAX_CPUS:-64}"
MAX_WALLTIME="${MAX_WALLTIME:-24:00:00}"

# Note:
# SLURM account, partition, QOS, and e-mail settings are not
# defined here because #SBATCH directives are interpreted by
# Slurm before this shell configuration is sourced.
