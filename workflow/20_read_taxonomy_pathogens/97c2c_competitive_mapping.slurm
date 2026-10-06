#!/bin/bash

#SBATCH --job-name=97c2c_compmap
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=12G
#SBATCH --time=03:00:00
#SBATCH --array=0-10%3

#SBATCH --mail-type=END,FAIL

#SBATCH --output=97c2c_compmap_%A_%a.out
#SBATCH --error=97c2c_compmap_%A_%a.err

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

export LC_ALL=C
export LANG=C
export LANGUAGE=C
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"

WORKER="${SCRIPT_DIR}/97c2c_competitive_mapping_worker.sh"

srun \
    --exact \
    -n 1 \
    -c "$SLURM_CPUS_PER_TASK" \
    "$WORKER" \
    "$SLURM_ARRAY_TASK_ID"

echo
echo "============================================================"
echo "97C2C ARRAY TASK ${SLURM_ARRAY_TASK_ID} FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"
