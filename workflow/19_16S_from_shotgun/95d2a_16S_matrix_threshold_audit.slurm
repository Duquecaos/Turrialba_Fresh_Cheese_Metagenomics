#!/bin/bash

#SBATCH --job-name=95d2a_16S_audit
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=95d2a_16S_audit_%j.out
#SBATCH --error=95d2a_16S_audit_%j.err

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

SCRIPT="${SCRIPT_DIR}/95d2a_16S_matrix_threshold_audit.py"

OUT="${ROOT}/97_16S_shotgun/95D2A_matrix_threshold_audit"

rm -rf "$OUT"

echo "============================================================"
echo "95D2A - 16S MATRIX + THRESHOLD AUDIT"
echo "Inicio: $(date)"
echo "Job: ${SLURM_JOB_ID}"
echo "============================================================"

python3 \
    "$SCRIPT" \
    "$ROOT"

[[ -s "${OUT}/95D2A_global_summary.tsv" ]] || {
    echo "ERROR: global summary missing" >&2
    exit 2
}

[[ -s "${OUT}/95D2A_threshold_sensitivity_by_sample.tsv" ]] || {
    echo "ERROR: threshold audit missing" >&2
    exit 3
}

[[ -s "${OUT}/95D2A_cluster_matrix_normalized_mean_depth_all.tsv" ]] || {
    echo "ERROR: cluster matrix missing" >&2
    exit 4
}

printf \
"COMPLETED\t%s\tjob=%s\n" \
"$(date -Iseconds)" \
"${SLURM_JOB_ID}" \
> "${OUT}/95D2A_COMPLETE.ok"

echo
echo "============================================================"
echo "GLOBAL SUMMARY"
echo "============================================================"

column -t -s $'\t' \
    "${OUT}/95D2A_global_summary.tsv"

echo
echo "============================================================"
echo "95D2A FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"

