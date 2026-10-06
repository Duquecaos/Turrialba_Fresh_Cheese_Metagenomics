#!/bin/bash

#SBATCH --job-name=98d1c_bact_matrix
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:20:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=98d1c_bact_matrix_%j.out
#SBATCH --error=98d1c_bact_matrix_%j.err

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

srun \
    --exact \
    -n 1 \
    -c "$SLURM_CPUS_PER_TASK" \
    "${SCRIPT_DIR}/98d1c_build_final_matrices.sh"

echo
echo "============================================================"
echo "98D1C FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"
