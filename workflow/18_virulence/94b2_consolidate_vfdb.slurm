#!/bin/bash

#SBATCH --job-name=94b2_VFDB_cons
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=94b2_VFDB_cons_%j.out
#SBATCH --error=94b2_VFDB_cons_%j.err

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

PYTHON_SCRIPT="${SCRIPT_DIR}/94b2_consolidate_vfdb.py"

OUT="${ROOT}/96_virulence/94B2_VFDB_consolidated"

rm -rf "$OUT"

echo "============================================================"
echo "94B2 - VFDB CONSOLIDATION"
echo "Inicio: $(date)"
echo "Job: ${SLURM_JOB_ID}"
echo "============================================================"

python3 \
    "$PYTHON_SCRIPT" \
    "$ROOT"

[[ -s "${OUT}/94B2_global_summary.tsv" ]] || {
    echo "ERROR: global summary faltante" >&2
    exit 2
}

[[ -s "${OUT}/94B2_community_high90_90.tsv" ]] || {
    echo "ERROR: high90/90 table faltante" >&2
    exit 3
}

N80=$(
    tail -n +2 \
    "${OUT}/94B2_community_all80_80.tsv" \
    | wc -l
)

if [[ "$N80" -ne 1075 ]]
then
    echo "ERROR: 80/80 rows=$N80 expected=1075" >&2
    exit 4
fi

printf \
"COMPLETED\t%s\tjob=%s\n" \
"$(date -Iseconds)" \
"${SLURM_JOB_ID}" \
> "${OUT}/94B2_COMPLETE.ok"

echo
echo "============================================================"
echo "GLOBAL SUMMARY"
echo "============================================================"

column -t -s $'\t' \
    "${OUT}/94B2_global_summary.tsv"

echo
echo "============================================================"
echo "FINAL18 SUMMARY"
echo "============================================================"

column -t -s $'\t' \
    "${OUT}/94B2_final18_summary.tsv"

echo
echo "============================================================"
echo "94B2 FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"

