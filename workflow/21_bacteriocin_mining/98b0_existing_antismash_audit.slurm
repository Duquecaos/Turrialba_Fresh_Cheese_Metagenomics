#!/bin/bash

#SBATCH --job-name=98b0_as_existing
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=98b0_as_existing_%j.out
#SBATCH --error=98b0_as_existing_%j.err

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

SCRIPT="${SCRIPT_DIR}/98b0_existing_antismash_audit.sh"

srun \
    --exact \
    -n 1 \
    -c "$SLURM_CPUS_PER_TASK" \
    "$SCRIPT"

echo
echo "============================================================"
echo "98B0 FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"
