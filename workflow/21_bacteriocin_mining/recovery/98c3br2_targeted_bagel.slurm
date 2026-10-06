#!/bin/bash

#SBATCH --job-name=98c3br2_bagel
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=98c3br2_bagel_%j.out
#SBATCH --error=98c3br2_bagel_%j.err

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

srun \
    --exact \
    -n 1 \
    -c "$SLURM_CPUS_PER_TASK" \
    "${SCRIPT_DIR}/98c3br2_targeted_bagel.sh"

echo
echo "============================================================"
echo "98C3BR2 FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"
