#!/bin/bash

#SBATCH --job-name=95d2b_16S_final
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=95d2b_16S_final_%j.out
#SBATCH --error=95d2b_16S_final_%j.err

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

SCRIPT="${SCRIPT_DIR}/95d2b_16S_final_support.py"

OUT="${ROOT}/97_16S_shotgun/95D2B_final_support"

rm -rf "$OUT"

echo "============================================================"
echo "95D2B - FINAL 16S MAPPING SUPPORT"
echo "Inicio: $(date)"
echo "Job: ${SLURM_JOB_ID}"
echo "============================================================"

python3 \
    "$SCRIPT" \
    "$ROOT"

[[ -s "${OUT}/95D2B_global_summary.tsv" ]] || {
    echo "ERROR: global summary missing" >&2
    exit 2
}

[[ -s "${OUT}/95D2B_locus_sample_final_support.tsv" ]] || {
    echo "ERROR: classification table missing" >&2
    exit 3
}

N=$(
    tail -n +2 \
    "${OUT}/95D2B_locus_sample_final_support.tsv" \
    | wc -l
)

if [[ "$N" -ne 363 ]]
then
    echo "ERROR: classifications=$N expected=363" >&2
    exit 4
fi

printf \
"COMPLETED\t%s\tjob=%s\n" \
"$(date -Iseconds)" \
"${SLURM_JOB_ID}" \
> "${OUT}/95D2B_COMPLETE.ok"

echo
echo "============================================================"
echo "GLOBAL SUMMARY"
echo "============================================================"

column -t -s $'\t' \
    "${OUT}/95D2B_global_summary.tsv"

echo
echo "============================================================"
echo "95D2B FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"

