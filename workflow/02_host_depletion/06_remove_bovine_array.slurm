#!/bin/bash

#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=08:00:00

#SBATCH --array=1-18%2
#SBATCH --job-name=hostrm

#SBATCH --output=%x-%A_%a.out
#SBATCH --error=%x-%A_%a.err

#SBATCH --mail-type=FAIL

set -euo pipefail

# ============================================================
# CONFIGURACION
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

LINE=$((SLURM_ARRAY_TASK_ID + 1))

SAMPLE=$(awk -F '\t' -v line="${LINE}" \
    'NR==line {print $1}' "${SAMPLES}")

if [[ -z "${SAMPLE}" ]]; then
    echo "ERROR: no se pudo obtener muestra para array ID ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

INPUT="${SCRATCH_DIR}/02_fastp/${SAMPLE}"

R1="${INPUT}/${SAMPLE}_1.clean.fastq.gz"
R2="${INPUT}/${SAMPLE}_2.clean.fastq.gz"

INDEX_PREFIX="${SCRATCH_DIR}/reference/Bos_taurus/bowtie2_index/ARS-UCD2.0"

OUT="${SCRATCH_DIR}/04_host_removed/${SAMPLE}"
TMP="${SCRATCH_DIR}/tmp/${SAMPLE}_host"

mkdir -p "${OUT}"
mkdir -p "${TMP}"

HOSTFREE_R1="${OUT}/${SAMPLE}_1.hostfree.fastq.gz"
HOSTFREE_R2="${OUT}/${SAMPLE}_2.hostfree.fastq.gz"

BT2_LOG="${OUT}/${SAMPLE}_bowtie2_bovine.log"
SUMMARY="${OUT}/${SAMPLE}_host_removal_summary.tsv"

# ============================================================
# INFO
# ============================================================

echo "============================================================"
echo "BOS TAURUS HOST REMOVAL"
echo "============================================================"
echo "Sample:       ${SAMPLE}"
echo "Array ID:     ${SLURM_ARRAY_TASK_ID}"
echo "Job ID:       ${SLURM_JOB_ID}"
echo "Nodo:         $(hostname)"
echo "Inicio:       $(date)"
echo "CPUs:         ${SLURM_CPUS_PER_TASK}"
echo "============================================================"
echo

# ============================================================
# VALIDACIONES
# ============================================================

if [[ ! -s "${R1}" ]]; then
    echo "ERROR: no existe ${R1}" >&2
    exit 1
fi

if [[ ! -s "${R2}" ]]; then
    echo "ERROR: no existe ${R2}" >&2
    exit 1
fi

if [[ ! -s "${INDEX_PREFIX}.1.bt2" ]]; then
    echo "ERROR: índice Bowtie2 no encontrado" >&2
    exit 1
fi

# ============================================================
# EVITAR REPETIR RESULTADOS COMPLETOS
# ============================================================

if [[ -s "${HOSTFREE_R1}" && \
      -s "${HOSTFREE_R2}" && \
      -s "${BT2_LOG}" ]]; then

    echo "Resultados host-free existentes para ${SAMPLE}."
    echo "Se omite el mapeo."

else

    rm -f "${HOSTFREE_R1}" "${HOSTFREE_R2}"

    echo "===== RUNNING HOST REMOVAL ====="

    srun --exact \
        -n 1 \
        -c "${SLURM_CPUS_PER_TASK}" \
        bash -c "

    apptainer exec \
        --bind '${SCRATCH_DIR}:${SCRATCH_DIR}' \
        '${BOWTIE2}' \
        bowtie2 \
        --very-sensitive \
        --threads 8 \
        -x '${INDEX_PREFIX}' \
        -1 '${R1}' \
        -2 '${R2}' \
        2> '${BT2_LOG}' \
    |
    apptainer exec \
        --bind '${SCRATCH_DIR}:${SCRATCH_DIR}' \
        '${SAMTOOLS}' \
        samtools view \
        -@ 2 \
        -b \
        -f 12 \
        -F 256 \
        - \
    |
    apptainer exec \
        --bind '${SCRATCH_DIR}:${SCRATCH_DIR}' \
        '${SAMTOOLS}' \
        samtools sort \
        -n \
        -@ 2 \
        -m 2G \
        -T '${TMP}/${SAMPLE}' \
        - \
    |
    apptainer exec \
        --bind '${SCRATCH_DIR}:${SCRATCH_DIR}' \
        '${SAMTOOLS}' \
        samtools fastq \
        -@ 2 \
        -1 '${HOSTFREE_R1}' \
        -2 '${HOSTFREE_R2}' \
        -0 /dev/null \
        -s /dev/null \
        -n \
        -

    "

fi

# ============================================================
# VALIDAR OUTPUT
# ============================================================

if [[ ! -s "${HOSTFREE_R1}" || ! -s "${HOSTFREE_R2}" ]]; then
    echo "ERROR: FASTQ host-free incompleto para ${SAMPLE}" >&2
    exit 1
fi

# ============================================================
# OBTENER ESTADISTICAS DEL LOG DE BOWTIE2
# ============================================================

CLEAN_PAIRS=$(awk '/reads; of these:/ {print $1; exit}' "${BT2_LOG}")

ALIGN_RATE=$(awk '/overall alignment rate/ {gsub("%","",$1); print $1}' "${BT2_LOG}")

# Contar pares host-free:
# número de líneas FASTQ / 4
# zcat se realiza solo sobre R1 porque R1 y R2 deben ser equivalentes.

HOSTFREE_PAIRS=$(gzip -cd "${HOSTFREE_R1}" | awk 'END {print NR/4}')

HOSTFREE_PCT=$(awk \
    -v clean="${CLEAN_PAIRS}" \
    -v hostfree="${HOSTFREE_PAIRS}" \
    'BEGIN {printf "%.4f", (hostfree/clean)*100}')

REMOVED_PCT=$(awk \
    -v keep="${HOSTFREE_PCT}" \
    'BEGIN {printf "%.4f", 100-keep}')

# ============================================================
# SUMMARY
# ============================================================

printf "sample\tclean_pairs\thostfree_pairs\thostfree_pct\tremoved_pct\tbowtie2_alignment_rate\n" \
    > "${SUMMARY}"

printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${SAMPLE}" \
    "${CLEAN_PAIRS}" \
    "${HOSTFREE_PAIRS}" \
    "${HOSTFREE_PCT}" \
    "${REMOVED_PCT}" \
    "${ALIGN_RATE}" \
    >> "${SUMMARY}"

echo
echo "===== SUMMARY ====="
column -t "${SUMMARY}"
echo

echo "===== HOST-FREE FILES ====="
ls -lh "${HOSTFREE_R1}" "${HOSTFREE_R2}"
echo

# ============================================================
# LIMPIEZA
# ============================================================

rm -rf "${TMP}"

echo "============================================================"
echo "HOST REMOVAL COMPLETADO"
echo "Sample: ${SAMPLE}"
echo "Fin: $(date)"
echo "============================================================"
