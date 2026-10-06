#!/bin/bash

#SBATCH --job-name=95c2b_16S_final
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=95c2b_16S_final_%j.out
#SBATCH --error=95c2b_16S_final_%j.err

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

PYTHON_SCRIPT="${SCRIPT_DIR}/95c2b_final_bacterial_filter.py"

OUT="${ROOT}/97_16S_shotgun/95C2B_final_bacterial_catalog"

rm -rf "$OUT"

echo "============================================================"
echo "95C2B - FINAL BACTERIAL 16S FILTER"
echo "Inicio: $(date)"
echo "Job: ${SLURM_JOB_ID}"
echo "============================================================"

python3 \
    "$PYTHON_SCRIPT" \
    "$ROOT"

[[ -s "${OUT}/95C2B_global_summary.tsv" ]] || {
    echo "ERROR: global summary faltante" >&2
    exit 2
}

[[ -s "${OUT}/95C2B_final_bacterial_quantification_targets.fasta" ]] || {
    echo "ERROR: final target FASTA faltante/vacio" >&2
    exit 3
}

CONTRADICTIONS=$(
    tail -n +2 \
    "${OUT}/95C2B_remaining_filter_contradictions.tsv" \
    | wc -l
)

if [[ "$CONTRADICTIONS" -ne 0 ]]
then
    echo "ERROR: filter contradictions=$CONTRADICTIONS" >&2
    exit 4
fi

printf \
"COMPLETED\t%s\tjob=%s\n" \
"$(date -Iseconds)" \
"${SLURM_JOB_ID}" \
> "${OUT}/95C2B_COMPLETE.ok"

echo
echo "============================================================"
echo "GLOBAL SUMMARY"
echo "============================================================"

column -t -s $'\t' \
    "${OUT}/95C2B_global_summary.tsv"

echo
echo "============================================================"
echo "95C2B FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"

