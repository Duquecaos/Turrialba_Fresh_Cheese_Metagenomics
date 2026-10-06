#!/bin/bash

#SBATCH --job-name=94c_VFDB_mobile
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=94c_VFDB_mobile_%j.out
#SBATCH --error=94c_VFDB_mobile_%j.err

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

PYTHON_SCRIPT="${SCRIPT_DIR}/94c_integrate_vfdb_mobile_context.py"

OUT="${ROOT}/96_virulence/94C_VFDB_mobile_context"

rm -rf "$OUT"

echo "============================================================"
echo "94C - VFDB + MOBILE CONTEXT"
echo "Inicio: $(date)"
echo "Job: ${SLURM_JOB_ID}"
echo "============================================================"

python3 \
    "$PYTHON_SCRIPT" \
    "$ROOT"

[[ -s "${OUT}/94C_global_summary.tsv" ]] || {
    echo "ERROR: global summary faltante" >&2
    exit 2
}

[[ -s "${OUT}/94C_VFDB_high90_90_mobile_integrated.tsv" ]] || {
    echo "ERROR: master table faltante" >&2
    exit 3
}

N=$(
    tail -n +2 \
    "${OUT}/94C_VFDB_high90_90_mobile_integrated.tsv" \
    | wc -l
)

if [[ "$N" -ne 51 ]]
then
    echo "ERROR: master rows=$N expected=51" >&2
    exit 4
fi

printf \
"COMPLETED\t%s\tjob=%s\n" \
"$(date -Iseconds)" \
"${SLURM_JOB_ID}" \
> "${OUT}/94C_COMPLETE.ok"

echo
echo "============================================================"
echo "GLOBAL SUMMARY"
echo "============================================================"

column -t -s $'\t' \
    "${OUT}/94C_global_summary.tsv"

echo
echo "============================================================"
echo "94C FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"

