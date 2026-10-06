#!/bin/bash

#SBATCH --job-name=98c1r_bact_fix
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=98c1r_bact_fix_%j.out
#SBATCH --error=98c1r_bact_fix_%j.err

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

srun \
    --exact \
    -n 1 \
    -c "$SLURM_CPUS_PER_TASK" \
    "${SCRIPT_DIR}/98c1r_reintegrate_corrected.sh"

echo
echo "============================================================"
echo "98C1R FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"
