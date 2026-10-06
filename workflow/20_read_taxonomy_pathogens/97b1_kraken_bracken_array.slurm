#!/bin/bash

#SBATCH --job-name=97b1_k2_all
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=160G
#SBATCH --time=06:00:00

#SBATCH --array=0-17%2

#SBATCH --mail-type=END,FAIL

#SBATCH --output=97b1_k2_all_%A_%a.out
#SBATCH --error=97b1_k2_all_%A_%a.err

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

export LC_ALL=C
export LANG=C
export LANGUAGE=C

export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"

WORKER="${SCRIPT_DIR}/97b1_kraken_bracken_worker.sh"

srun \
    --exact \
    -n 1 \
    -c "$SLURM_CPUS_PER_TASK" \
    "$WORKER"

echo
echo "============================================================"
echo "97B1 ARRAY TASK ${SLURM_ARRAY_TASK_ID} FINALIZO"
echo "ES SEGURO SALIR"
echo "============================================================"

