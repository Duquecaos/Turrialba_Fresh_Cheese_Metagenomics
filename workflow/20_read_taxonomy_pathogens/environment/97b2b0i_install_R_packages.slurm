#!/bin/bash

#SBATCH --job-name=97b2b0i_R_pkgs
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=97b2b0i_R_pkgs_%j.out
#SBATCH --error=97b2b0i_R_pkgs_%j.err

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

export LC_ALL=C
export LANG=C
export LANGUAGE=C
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"

SCRIPT="${SCRIPT_DIR}/97b2b0i_install_R_packages.sh"

srun \
    --exact \
    -n 1 \
    -c "$SLURM_CPUS_PER_TASK" \
    "$SCRIPT"

echo
echo "============================================================"
echo "97B2B0I FINALIZO CORRECTAMENTE"
echo "R PACKAGES READY"
echo "ES SEGURO SALIR"
echo "============================================================"

