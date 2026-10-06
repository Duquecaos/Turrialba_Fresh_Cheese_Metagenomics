#!/bin/bash

#SBATCH --job-name=95b2_16S_cons
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=95b2_16S_cons_%j.out
#SBATCH --error=95b2_16S_cons_%j.err

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

PYTHON_SCRIPT="${SCRIPT_DIR}/95b2_consolidate_16S.py"

OUT="${ROOT}/97_16S_shotgun/95B2_consolidated"

rm -rf "$OUT"

echo "============================================================"
echo "95B2 - CONSOLIDACION 16S SHOTGUN"
echo "Inicio: $(date)"
echo "Job: ${SLURM_JOB_ID}"
echo "============================================================"

python3 \
    "$PYTHON_SCRIPT" \
    "$ROOT"

[[ -s "${OUT}/95B2_global_summary.tsv" ]] || {
    echo "ERROR: global summary faltante" >&2
    exit 2
}

[[ -s "${OUT}/95B2_coassembly_16S_canonical.tsv" ]] || {
    echo "ERROR: catalogo canónico faltante" >&2
    exit 3
}

N=$(
    tail -n +2 \
    "${OUT}/95B2_coassembly_16S_canonical.tsv" \
    | wc -l
)

if [[ "$N" -ne 158 ]]
then
    echo "ERROR: primary 16S rows=$N expected=158" >&2
    exit 4
fi

printf \
"COMPLETED\t%s\tjob=%s\n" \
"$(date -Iseconds)" \
"${SLURM_JOB_ID}" \
> "${OUT}/95B2_COMPLETE.ok"

echo
echo "============================================================"
echo "GLOBAL SUMMARY"
echo "============================================================"

column -t -s $'\t' \
"${OUT}/95B2_global_summary.tsv"

echo
echo "============================================================"
echo "95B2 FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"

