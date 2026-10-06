#!/bin/bash

#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=02:00:00

#SBATCH --job-name=86c_eggnog_merge
#SBATCH --output=86c_eggnog_merge_%j.out
#SBATCH --error=86c_eggnog_merge_%j.err

#SBATCH --mail-type=END,FAIL

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

echo "============================================================"
echo "PASO 86c - CONSOLIDACION EGGNOG"
echo "Inicio: $(date)"
echo "Host: $(hostname)"
echo "Job: ${SLURM_JOB_ID}"
echo "============================================================"

srun --exact \
    -n 1 \
    -c 1 \
    python3 \
    "${SCRIPT_DIR}/86c_consolidate_incomplete_eggnog.py"

echo
echo "============================================================"
echo "PASO 86c FINALIZÓ CORRECTAMENTE."
echo "ES SEGURO SALIR."
echo "Fin: $(date)"
echo "============================================================"
