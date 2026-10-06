#!/bin/bash

#SBATCH --job-name=98b2_exact_bact
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=01:00:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=98b2_exact_bact_%j.out
#SBATCH --error=98b2_exact_bact_%j.err

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

srun \
    --exact \
    -n 1 \
    -c "$SLURM_CPUS_PER_TASK" \
    "${SCRIPT_DIR}/98b2_exact_candidate_review.sh"

echo
echo "============================================================"
echo "98B2 FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"
