#!/bin/bash

#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --job-name=multiqc_preproc
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

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

OUT="${SCRATCH_DIR}/QC_multiqc"

mkdir -p "${OUT}"

echo "=========================================="
echo "MULTIQC - SHOTGUN PREPROCESSING"
echo "=========================================="
echo "Inicio: $(date)"
echo "Nodo:   $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo
echo "Output: ${OUT}"
echo "=========================================="
echo

echo "===== INPUT COUNTS ====="

echo -n "FastQC RAW zip: "
find "${SCRATCH_DIR}/01_fastqc_raw" \
    -name "*_fastqc.zip" | wc -l

echo -n "fastp JSON: "
find "${SCRATCH_DIR}/02_fastp" \
    -name "*_fastp.json" | wc -l

echo -n "FastQC CLEAN zip: "
find "${SCRATCH_DIR}/03_fastqc_clean" \
    -name "*_fastqc.zip" | wc -l

echo

echo "===== RUNNING MULTIQC ====="

srun --exact \
    -n 1 \
    -c "${SLURM_CPUS_PER_TASK}" \
    apptainer exec \
    --bind "${SCRATCH_DIR}:${SCRATCH_DIR}" \
    "${MULTIQC}" \
    multiqc \
    "${SCRATCH_DIR}/01_fastqc_raw" \
    "${SCRATCH_DIR}/02_fastp" \
    "${SCRATCH_DIR}/03_fastqc_clean" \
    --outdir "${OUT}" \
    --filename "multiqc_preprocessing.html" \
    --title "Shotgun Turrialba - Preprocessing QC" \
    --force

echo

echo "===== MULTIQC OUTPUT ====="

ls -lh "${OUT}"

echo

echo "=========================================="
echo "MULTIQC COMPLETADO"
echo "Fin: $(date)"
echo "=========================================="
