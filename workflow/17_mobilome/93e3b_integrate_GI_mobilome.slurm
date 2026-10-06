#!/bin/bash

#SBATCH --job-name=93e3b_GI_mobile
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=93e3b_GI_mobile_%j.out
#SBATCH --error=93e3b_GI_mobile_%j.err

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

PYTHON_SCRIPT="${SCRIPT_DIR}/93e3b_integrate_GI_mobilome.py"

OUT="${ROOT}/95_mobilome/93E3B_GI_mobile_integration"

rm -rf "$OUT"

echo "============================================================"
echo "93E3B - GI + MOBILOME INTEGRATION"
echo "Inicio: $(date)"
echo "Job: ${SLURM_JOB_ID}"
echo "============================================================"

[[ -s "$PYTHON_SCRIPT" ]] || {
    echo "ERROR: falta $PYTHON_SCRIPT" >&2
    exit 2
}

python3 "$PYTHON_SCRIPT" "$ROOT"


[[ -s "${OUT}/93E3B_global_summary.tsv" ]] || {
    echo "ERROR: global summary faltante" >&2
    exit 3
}

[[ -s "${OUT}/93E3B_GI_mobile_integrated.tsv" ]] || {
    echo "ERROR: master table faltante" >&2
    exit 4
}


N=$(
    tail -n +2 \
    "${OUT}/93E3B_GI_mobile_integrated.tsv" \
    | wc -l
)

if [[ "$N" -ne 1251 ]]
then
    echo "ERROR: master rows=$N, expected=1251" >&2
    exit 5
fi


printf \
"COMPLETED\t%s\tjob=%s\n" \
"$(date -Iseconds)" \
"${SLURM_JOB_ID}" \
> "${OUT}/93E3B_COMPLETE.ok"


echo
echo "============================================================"
echo "GLOBAL SUMMARY"
echo "============================================================"

column -t -s $'\t' \
"${OUT}/93E3B_global_summary.tsv"


echo
echo "============================================================"
echo "EVIDENCE CLASSES"
echo "============================================================"

column -t -s $'\t' \
"${OUT}/93E3B_evidence_class_summary.tsv"


echo
echo "============================================================"
echo "93E3B FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"

