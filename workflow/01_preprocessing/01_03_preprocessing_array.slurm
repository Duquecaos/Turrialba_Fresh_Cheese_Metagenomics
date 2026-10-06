#!/bin/bash

#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --array=1-18%4
#SBATCH --job-name=shotgun_preproc
#SBATCH --output=%x-%A_%a.out
#SBATCH --error=%x-%A_%a.err

set -euo pipefail

# ============================================================
# CONFIGURACIÓN GENERAL
# ============================================================

# Locate repository configuration portably.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
PROJECT_CONFIG="${PROJECT_CONFIG:-${REPO_ROOT}/config/config.sh}"

if [[ ! -r "$PROJECT_CONFIG" ]]; then
    echo "ERROR: configuration file not found: $PROJECT_CONFIG" >&2
    echo "Create config/config.sh from config/config.example.sh before running this script." >&2
    exit 2
fi

source "$PROJECT_CONFIG"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

# ============================================================
# OBTENER MUESTRA DESDE metadata/samples.tsv
# SLURM_ARRAY_TASK_ID 1 corresponde a la primera muestra
# ============================================================

LINE=$((SLURM_ARRAY_TASK_ID + 1))

SAMPLE=$(awk -F '\t' -v line="${LINE}" \
    'NR==line {print $1}' "${SAMPLES}")

R1=$(awk -F '\t' -v line="${LINE}" \
    'NR==line {print $2}' "${SAMPLES}")

R2=$(awk -F '\t' -v line="${LINE}" \
    'NR==line {print $3}' "${SAMPLES}")

if [[ -z "${SAMPLE}" ]]; then
    echo "ERROR: no se pudo obtener la muestra para array ID ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

# ============================================================
# DIRECTORIOS
# ============================================================

RAW_QC_OUT="${SCRATCH_DIR}/01_fastqc_raw/${SAMPLE}"
FASTP_OUT="${SCRATCH_DIR}/02_fastp/${SAMPLE}"
CLEAN_QC_OUT="${SCRATCH_DIR}/03_fastqc_clean/${SAMPLE}"

mkdir -p \
    "${RAW_QC_OUT}" \
    "${FASTP_OUT}" \
    "${CLEAN_QC_OUT}"

CLEAN_R1="${FASTP_OUT}/${SAMPLE}_1.clean.fastq.gz"
CLEAN_R2="${FASTP_OUT}/${SAMPLE}_2.clean.fastq.gz"

FASTP_HTML="${FASTP_OUT}/${SAMPLE}_fastp.html"
FASTP_JSON="${FASTP_OUT}/${SAMPLE}_fastp.json"

# ============================================================
# INFORMACIÓN
# ============================================================

echo "============================================================"
echo "SHOTGUN PREPROCESSING"
echo "============================================================"
echo "Sample:       ${SAMPLE}"
echo "Array ID:     ${SLURM_ARRAY_TASK_ID}"
echo "Job ID:       ${SLURM_JOB_ID}"
echo "Nodo:         $(hostname)"
echo "Inicio:       $(date)"
echo "CPUs:         ${SLURM_CPUS_PER_TASK}"
echo
echo "R1:           ${R1}"
echo "R2:           ${R2}"
echo "============================================================"
echo

# ============================================================
# VALIDAR INPUT
# ============================================================

if [[ ! -f "${R1}" ]]; then
    echo "ERROR: no existe R1: ${R1}" >&2
    exit 1
fi

if [[ ! -f "${R2}" ]]; then
    echo "ERROR: no existe R2: ${R2}" >&2
    exit 1
fi

echo "===== INPUT ====="
ls -lh "${R1}" "${R2}"
echo

# ============================================================
# 1. FASTQC RAW
# ============================================================

RAW_QC_R1="${RAW_QC_OUT}/${SAMPLE}_1_fastqc.zip"
RAW_QC_R2="${RAW_QC_OUT}/${SAMPLE}_2_fastqc.zip"

if [[ -s "${RAW_QC_R1}" && -s "${RAW_QC_R2}" ]]; then

    echo "===== FASTQC RAW ====="
    echo "Resultados existentes. Se omite este paso."
    echo

else

    echo "===== RUNNING FASTQC RAW ====="

    srun --exact \
        -n 1 \
        -c "${SLURM_CPUS_PER_TASK}" \
        apptainer exec \
        --bind "${RAW_DIR}:${RAW_DIR}" \
        --bind "${SCRATCH_DIR}:${SCRATCH_DIR}" \
        "${FASTQC}" \
        fastqc \
        --threads "${SLURM_CPUS_PER_TASK}" \
        --outdir "${RAW_QC_OUT}" \
        "${R1}" \
        "${R2}"

    echo

fi

# ============================================================
# 2. FASTP
# ============================================================

if [[ -s "${CLEAN_R1}" && -s "${CLEAN_R2}" && \
      -s "${FASTP_JSON}" ]]; then

    echo "===== FASTP ====="
    echo "Resultados existentes. Se omite este paso."
    echo

else

    echo "===== RUNNING FASTP ====="

    srun --exact \
        -n 1 \
        -c "${SLURM_CPUS_PER_TASK}" \
        apptainer exec \
        --bind "${RAW_DIR}:${RAW_DIR}" \
        --bind "${SCRATCH_DIR}:${SCRATCH_DIR}" \
        "${FASTP}" \
        fastp \
        --in1 "${R1}" \
        --in2 "${R2}" \
        --out1 "${CLEAN_R1}" \
        --out2 "${CLEAN_R2}" \
        --thread "${SLURM_CPUS_PER_TASK}" \
        --detect_adapter_for_pe \
        --trim_poly_g \
        --poly_g_min_len 10 \
        --qualified_quality_phred 20 \
        --unqualified_percent_limit 40 \
        --n_base_limit 5 \
        --length_required 50 \
        --html "${FASTP_HTML}" \
        --json "${FASTP_JSON}" \
        --report_title "fastp ${SAMPLE}"

    echo

fi

# ============================================================
# VALIDAR FASTP
# ============================================================

if [[ ! -s "${CLEAN_R1}" ]]; then
    echo "ERROR: fastp no produjo ${CLEAN_R1}" >&2
    exit 1
fi

if [[ ! -s "${CLEAN_R2}" ]]; then
    echo "ERROR: fastp no produjo ${CLEAN_R2}" >&2
    exit 1
fi

echo "===== CLEAN READS ====="
ls -lh "${CLEAN_R1}" "${CLEAN_R2}"
echo

# ============================================================
# 3. FASTQC CLEAN
# ============================================================

CLEAN_QC_R1="${CLEAN_QC_OUT}/${SAMPLE}_1.clean_fastqc.zip"
CLEAN_QC_R2="${CLEAN_QC_OUT}/${SAMPLE}_2.clean_fastqc.zip"

if [[ -s "${CLEAN_QC_R1}" && -s "${CLEAN_QC_R2}" ]]; then

    echo "===== FASTQC CLEAN ====="
    echo "Resultados existentes. Se omite este paso."
    echo

else

    echo "===== RUNNING FASTQC CLEAN ====="

    srun --exact \
        -n 1 \
        -c "${SLURM_CPUS_PER_TASK}" \
        apptainer exec \
        --bind "${SCRATCH_DIR}:${SCRATCH_DIR}" \
        "${FASTQC}" \
        fastqc \
        --threads "${SLURM_CPUS_PER_TASK}" \
        --outdir "${CLEAN_QC_OUT}" \
        "${CLEAN_R1}" \
        "${CLEAN_R2}"

    echo

fi

# ============================================================
# RESUMEN
# ============================================================

echo "===== FINAL OUTPUT ====="

echo
echo "--- RAW FASTQC ---"
ls -lh "${RAW_QC_OUT}"

echo
echo "--- FASTP ---"
ls -lh "${FASTP_OUT}"

echo
echo "--- CLEAN FASTQC ---"
ls -lh "${CLEAN_QC_OUT}"

echo

echo "============================================================"
echo "PREPROCESSING COMPLETADO"
echo "Sample: ${SAMPLE}"
echo "Fin: $(date)"
echo "============================================================"
