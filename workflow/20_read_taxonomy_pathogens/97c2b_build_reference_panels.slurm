#!/bin/bash

#SBATCH --job-name=97c2b_ref_panels
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=97c2b_ref_panels_%j.out
#SBATCH --error=97c2b_ref_panels_%j.err

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

export LC_ALL=C
export LANG=C
export LANGUAGE=C
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"

SCRIPT="${SCRIPT_DIR}/97c2b_build_reference_panels.sh"

srun \
    --exact \
    -n 1 \
    -c "$SLURM_CPUS_PER_TASK" \
    "$SCRIPT"

echo
echo "============================================================"
echo "97C2B FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"
