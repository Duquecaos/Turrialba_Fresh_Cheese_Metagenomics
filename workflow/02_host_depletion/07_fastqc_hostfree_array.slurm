#!/bin/bash

#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00

#SBATCH --array=1-18%4
#SBATCH --job-name=qc_hostfree

#SBATCH --output=%x-%A_%a.out
#SBATCH --error=%x-%A_%a.err

#SBATCH --mail-type=FAIL

set -euo pipefail

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
    echo "ERROR: no se pudo obtener muestra" >&2
    exit 1
fi

INPUT="${SCRATCH_DIR}/04_host_removed/${SAMPLE}"

R1="${INPUT}/${SAMPLE}_1.hostfree.fastq.gz"
R2="${INPUT}/${SAMPLE}_2.hostfree.fastq.gz"

OUT="${INPUT}/fastqc"

mkdir -p "${OUT}"

ZIP1="${OUT}/${SAMPLE}_1.hostfree_fastqc.zip"
ZIP2="${OUT}/${SAMPLE}_2.hostfree_fastqc.zip"

echo "============================================================"
echo "FASTQC HOST-FREE"
echo "Sample:   ${SAMPLE}"
echo "Job ID:   ${SLURM_JOB_ID}"
echo "Nodo:     $(hostname)"
echo "Inicio:   $(date)"
echo "============================================================"
echo

if [[ ! -s "${R1}" || ! -s "${R2}" ]]; then
    echo "ERROR: host-free FASTQ faltante para ${SAMPLE}" >&2
    exit 1
fi

if [[ -s "${ZIP1}" && -s "${ZIP2}" ]]; then

    echo "FastQC host-free ya existe para ${SAMPLE}."
    echo "Se omite."

else

    echo "===== RUNNING FASTQC ====="

    srun --exact \
        -n 1 \
        -c "${SLURM_CPUS_PER_TASK}" \
        apptainer exec \
        --bind "${SCRATCH_DIR}:${SCRATCH_DIR}" \
        "${FASTQC}" \
        fastqc \
        --threads "${SLURM_CPUS_PER_TASK}" \
        --outdir "${OUT}" \
        "${R1}" \
        "${R2}"

fi

echo
echo "===== OUTPUT ====="
ls -lh "${OUT}"

echo
echo "============================================================"
echo "FASTQC HOST-FREE COMPLETADO"
echo "Sample: ${SAMPLE}"
echo "Fin: $(date)"
echo "============================================================"
