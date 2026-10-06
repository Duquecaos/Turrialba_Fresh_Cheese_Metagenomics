#!/bin/bash

#SBATCH --job-name=95c2_16S_QC
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=95c2_16S_QC_%j.out
#SBATCH --error=95c2_16S_QC_%j.err

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

PYTHON_SCRIPT="${SCRIPT_DIR}/95c2_bacterial_16S_QC.py"

OUT="${ROOT}/97_16S_shotgun/95C2_bacterial_QC"

rm -rf "$OUT"

echo "============================================================"
echo "95C2 - BACTERIAL 16S TAXONOMIC QC"
echo "Inicio: $(date)"
echo "Job: ${SLURM_JOB_ID}"
echo "============================================================"

python3 \
    "$PYTHON_SCRIPT" \
    "$ROOT"

[[ -s "${OUT}/95C2_global_summary.tsv" ]] || {
    echo "ERROR: global summary faltante" >&2
    exit 2
}

[[ -s "${OUT}/95C2_exact_representatives_QC.tsv" ]] || {
    echo "ERROR: QC table faltante" >&2
    exit 3
}

N=$(
    tail -n +2 \
    "${OUT}/95C2_exact_representatives_QC.tsv" \
    | wc -l
)

if [[ "$N" -ne 148 ]]
then
    echo "ERROR: representatives=$N expected=148" >&2
    exit 4
fi

printf \
"COMPLETED\t%s\tjob=%s\n" \
"$(date -Iseconds)" \
"${SLURM_JOB_ID}" \
> "${OUT}/95C2_COMPLETE.ok"

echo
echo "============================================================"
echo "GLOBAL SUMMARY"
echo "============================================================"

column -t -s $'\t' \
    "${OUT}/95C2_global_summary.tsv"

echo
echo "============================================================"
echo "95C2 FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"

