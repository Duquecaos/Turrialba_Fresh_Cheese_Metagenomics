#!/bin/bash

#SBATCH --job-name=97b2a_tax_matrix
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=97b2a_tax_matrix_%j.out
#SBATCH --error=97b2a_tax_matrix_%j.err

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

export LC_ALL=C
export LANG=C
export LANGUAGE=C

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

SCRIPT="${SCRIPT_DIR}/97b2a_taxonomy_matrix_qc.py"

OUT="${ROOT}/98_read_taxonomy/97B2A_matrix_QC"

rm -rf "$OUT"

python3 \
    "$SCRIPT" \
    "$ROOT"


[[ -s "${OUT}/97B2A_global_summary.tsv" ]] || {
    echo "ERROR: global summary missing" >&2
    exit 2
}


[[ -s "${OUT}/97B2A_genus_relative_cellular_microbiota.tsv" ]] || {
    echo "ERROR: genus matrix missing" >&2
    exit 3
}


printf \
"COMPLETED\t%s\tjob=%s\n" \
"$(date -Iseconds)" \
"${SLURM_JOB_ID}" \
> "${OUT}/97B2A_COMPLETE.ok"


echo
echo "============================================================"
echo "97B2A GLOBAL SUMMARY"
echo "============================================================"

column -t -s $'\t' \
    "${OUT}/97B2A_global_summary.tsv"

echo
echo "============================================================"
echo "97B2A FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"

