#!/bin/bash

#SBATCH --job-name=98d1b_bact_map
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --array=0-17%8

#SBATCH --mail-type=END,FAIL

#SBATCH --output=98d1b_bact_map_%A_%a.out
#SBATCH --error=98d1b_bact_map_%A_%a.err

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"

srun \
    --exact \
    -n 1 \
    -c "$SLURM_CPUS_PER_TASK" \
    "${SCRIPT_DIR}/98d1b_map_one_sample.sh"

echo
echo "============================================================"
echo "98D1B ARRAY TASK ${SLURM_ARRAY_TASK_ID} FINALIZO"
echo "ES SEGURO SALIR"
echo "============================================================"
