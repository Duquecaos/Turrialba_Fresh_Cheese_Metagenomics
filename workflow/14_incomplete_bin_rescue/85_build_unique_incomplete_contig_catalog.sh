#!/bin/bash
#SBATCH --job-name=85_contigs
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=85_contigs_%j.out
#SBATCH --error=85_contigs_%j.err
#SBATCH --mail-type=END,FAIL

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

echo "============================================================"
echo "PASO 85 - CATALOGO UNICO DE CONTIGS"
echo "Inicio: $(date)"
echo "Host: $(hostname)"
echo "Job: ${SLURM_JOB_ID}"
echo "============================================================"

srun --exact \
  -n 1 \
  -c "$SLURM_CPUS_PER_TASK" \
  python3 \
  "${SCRIPT_DIR}/85_build_unique_incomplete_contig_catalog.py"

echo
echo "============================================================"
echo "PASO 85 FINALIZÓ CORRECTAMENTE."
echo "ES SEGURO SALIR."
echo "Fin: $(date)"
echo "============================================================"
