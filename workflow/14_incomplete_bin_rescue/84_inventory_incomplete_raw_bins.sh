#!/bin/bash
#SBATCH --job-name=84_bin_inventory
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=84_bin_inventory_%j.out
#SBATCH --error=84_bin_inventory_%j.err
#SBATCH --mail-type=END,FAIL

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

echo "============================================================"
echo "PASO 84 - INVENTARIO DE BINS INCOMPLETOS"
echo "Inicio: $(date)"
echo "Host: $(hostname)"
echo "Job: ${SLURM_JOB_ID}"
echo "============================================================"

srun --exact \
  -n 1 \
  -c "$SLURM_CPUS_PER_TASK" \
  python3 \
  "${SCRIPT_DIR}/84_inventory_incomplete_raw_bins.py"

echo
echo "============================================================"
echo "PASO 84 FINALIZÓ CORRECTAMENTE."
echo "ES SEGURO SALIR."
echo "Fin: $(date)"
echo "============================================================"
