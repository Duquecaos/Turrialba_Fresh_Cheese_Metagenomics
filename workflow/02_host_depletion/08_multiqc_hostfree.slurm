#!/bin/bash

#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00

#SBATCH --job-name=multiqc_hostfree

#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

#SBATCH --mail-type=END,FAIL

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

INPUT="${SCRATCH_DIR}/04_host_removed"
OUT="${SCRATCH_DIR}/QC_multiqc_hostfree"

mkdir -p "${OUT}"

echo "============================================================"
echo "MULTIQC HOST-FREE"
echo "Inicio: $(date)"
echo "Nodo:   $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "============================================================"
echo

echo -n "FastQC host-free encontrados: "

find "${INPUT}" \
    -name "*.hostfree_fastqc.zip" | wc -l

echo
echo "===== RUNNING MULTIQC ====="

srun --exact \
    -n 1 \
    -c "${SLURM_CPUS_PER_TASK}" \
    apptainer exec \
    --bind "${SCRATCH_DIR}:${SCRATCH_DIR}" \
    "${MULTIQC}" \
    multiqc \
    "${INPUT}" \
    --outdir "${OUT}" \
    --filename "multiqc_hostfree.html" \
    --title "Shotgun Turrialba - Host-free Reads" \
    --force

echo
echo "===== OUTPUT ====="

ls -lh "${OUT}"

echo
echo "============================================================"
echo "MULTIQC HOST-FREE COMPLETADO"
echo "Fin: $(date)"
echo "============================================================"
