#!/bin/bash

#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00

#SBATCH --array=1-6%2
#SBATCH --job-name=bt2_coasm_index

#SBATCH --output=%x-%A_%a.out
#SBATCH --error=%x-%A_%a.err

#SBATCH --mail-type=END,FAIL

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

COASSEMBLY_NAMES=(L1 L2 L3 M1 M2 M3)

ARRAY_INDEX=$((SLURM_ARRAY_TASK_ID - 1))

GROUP_NAME="${COASSEMBLY_NAMES[${ARRAY_INDEX}]}"

ASSEMBLY_DIR="${SCRATCH_DIR}/05_assembly/${GROUP_NAME}"

CONTIGS="${ASSEMBLY_DIR}/final.contigs.fa"

INDEX_DIR="${ASSEMBLY_DIR}/bowtie2_index"

INDEX_PREFIX="${INDEX_DIR}/${GROUP_NAME}"

mkdir -p "${INDEX_DIR}"

# ============================================================
# INFORMACION
# ============================================================

echo "============================================================"
echo "BOWTIE2 INDEX - COASSEMBLY"
echo "============================================================"
echo "Coassembly:   ${GROUP_NAME}"
echo "Array ID:     ${SLURM_ARRAY_TASK_ID}"
echo "Job ID:       ${SLURM_JOB_ID}"
echo "Nodo:         $(hostname)"
echo "Inicio:       $(date)"
echo "CPUs:         ${SLURM_CPUS_PER_TASK}"
echo
echo "Contigs:      ${CONTIGS}"
echo "Index:        ${INDEX_PREFIX}"
echo "============================================================"
echo

# ============================================================
# VALIDAR CONTIGS
# ============================================================

if [[ ! -s "${CONTIGS}" ]]; then
    echo "ERROR: no existe ${CONTIGS}" >&2
    exit 1
fi

echo "===== INPUT ====="
ls -lh "${CONTIGS}"
echo

# ============================================================
# SI EL INDICE YA EXISTE, NO REPETIR
# ============================================================

NFILES=$(find "${INDEX_DIR}" \
    -maxdepth 1 \
    \( -name "${GROUP_NAME}*.bt2" -o -name "${GROUP_NAME}*.bt2l" \) \
    2>/dev/null | wc -l)

if [[ "${NFILES}" -ge 6 ]]; then

    echo "Indice Bowtie2 ya existente."
    echo "Se omite bowtie2-build."

else

    echo "===== BUILDING INDEX ====="

    srun --exact \
        -n 1 \
        -c "${SLURM_CPUS_PER_TASK}" \
        apptainer exec \
        --bind "${SCRATCH_DIR}:${SCRATCH_DIR}" \
        "${BOWTIE2}" \
        bowtie2-build \
        --threads "${SLURM_CPUS_PER_TASK}" \
        "${CONTIGS}" \
        "${INDEX_PREFIX}"

fi

echo

# ============================================================
# VALIDAR
# ============================================================

echo "===== INDEX FILES ====="

ls -lh "${INDEX_DIR}"

echo

NFILES=$(find "${INDEX_DIR}" \
    -maxdepth 1 \
    \( -name "${GROUP_NAME}*.bt2" -o -name "${GROUP_NAME}*.bt2l" \) \
    | wc -l)

echo "Index files detected: ${NFILES}"

if [[ "${NFILES}" -lt 6 ]]; then
    echo "ERROR: índice incompleto para ${GROUP_NAME}" >&2
    exit 1
fi

echo
echo "============================================================"
echo "BOWTIE2 INDEX COMPLETADO"
echo "Coassembly: ${GROUP_NAME}"
echo "Fin: $(date)"
echo "============================================================"
