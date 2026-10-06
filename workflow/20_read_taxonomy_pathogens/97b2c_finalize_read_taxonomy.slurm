#!/bin/bash

#SBATCH --job-name=97b2c_tax_final
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:20:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=97b2c_tax_final_%j.out
#SBATCH --error=97b2c_tax_final_%j.err

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

export LC_ALL=C
export LANG=C
export LANGUAGE=C

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

SCRIPT="${SCRIPT_DIR}/97b2c_finalize_read_taxonomy.py"

OUT="${ROOT}/98_read_taxonomy/97B2C_final_taxonomy"

rm -rf "$OUT"

python3 \
    "$SCRIPT" \
    "$ROOT"

[[ -s "${OUT}/97B2C_global_summary.tsv" ]] || {
    echo "ERROR: summary missing" >&2
    exit 20
}

[[ -s "${OUT}/97B2C_thesis_results_snapshot.txt" ]] || {
    echo "ERROR: thesis snapshot missing" >&2
    exit 21
}

printf \
"COMPLETED\t%s\tjob=%s\n" \
"$(date -Iseconds)" \
"$SLURM_JOB_ID" \
> "${OUT}/97B2C_COMPLETE.ok"

echo
echo "============================================================"
echo "97B2C GLOBAL SUMMARY"
echo "============================================================"

column -t -s $'\t' \
"${OUT}/97B2C_global_summary.tsv"

echo
echo "============================================================"
echo "97B2C TAXONOMIA SHOTGUN CERRADA"
echo "ES SEGURO SALIR"
echo "============================================================"
