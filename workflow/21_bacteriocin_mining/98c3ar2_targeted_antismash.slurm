#!/bin/bash

#SBATCH --job-name=98c3ar2_antismash
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=06:00:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=98c3ar2_antismash_%j.out
#SBATCH --error=98c3ar2_antismash_%j.err

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"

srun \
    --exact \
    -n 1 \
    -c "$SLURM_CPUS_PER_TASK" \
    "${SCRIPT_DIR}/98c3a_targeted_antismash.sh"

echo
echo "============================================================"
echo "98C3AR2 FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"
