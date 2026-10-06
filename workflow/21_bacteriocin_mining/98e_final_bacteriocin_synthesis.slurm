#!/bin/bash

#SBATCH --job-name=98e_bact_final
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:20:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=98e_bact_final_%j.out
#SBATCH --error=98e_bact_final_%j.err

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

srun \
    --exact \
    -n 1 \
    -c "$SLURM_CPUS_PER_TASK" \
    "${SCRIPT_DIR}/98e_final_bacteriocin_synthesis.sh"

echo
echo "============================================================"
echo "98E FINALIZO CORRECTAMENTE"
echo "RAMA COMPUTACIONAL DE BACTERIOCINAS/RiPP CERRADA"
echo "ES SEGURO SALIR"
echo "============================================================"
