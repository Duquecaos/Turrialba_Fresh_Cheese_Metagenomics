#!/bin/bash

#SBATCH --job-name=97c2d_final
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:20:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=97c2d_final_%j.out
#SBATCH --error=97c2d_final_%j.err

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

SCRIPT="${SCRIPT_DIR}/97c2d_finalize_pathogen_validation.py"

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"
OUT="${ROOT}/98_read_taxonomy/97C2D_final_pathogen_validation"

rm -rf "$OUT"

python3 "$SCRIPT"

[[ -s "${OUT}/97C2D_final_validation.tsv" ]] || exit 20

printf \
"COMPLETED\t%s\tjob=%s\n" \
"$(date -Iseconds)" \
"$SLURM_JOB_ID" \
> "${OUT}/97C2D_COMPLETE.ok"

echo
echo "============================================================"
echo "97C2D FINAL VALIDATION"
echo "============================================================"

column -t -s $'\t' \
"${OUT}/97C2D_final_validation.tsv"

echo
echo "============================================================"
echo "97C2D FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"
