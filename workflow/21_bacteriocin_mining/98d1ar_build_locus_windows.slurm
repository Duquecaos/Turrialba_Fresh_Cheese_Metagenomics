#!/bin/bash

#SBATCH --job-name=98d1ar_bact_ref
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=98d1ar_bact_ref_%j.out
#SBATCH --error=98d1ar_bact_ref_%j.err

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"

srun \
    --exact \
    -n 1 \
    -c "$SLURM_CPUS_PER_TASK" \
    "${SCRIPT_DIR}/98d1ar_build_locus_windows.sh"

echo
echo "============================================================"
echo "98D1AR FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"
