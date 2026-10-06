#!/bin/bash

#SBATCH --job-name=97b2br_stats_fix
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=97b2br_stats_fix_%j.out
#SBATCH --error=97b2br_stats_fix_%j.err

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

export LC_ALL=C
export LANG=C
export LANGUAGE=C

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

source \
"${ROOT}/98_read_taxonomy/R_environment.sh"

SCRIPT="${SCRIPT_DIR}/97b2br_statistics_correction.R"

OUT="${ROOT}/98_read_taxonomy/97B2BR_statistics_correction"

rm -rf "$OUT"

Rscript \
    "$SCRIPT" \
    "$ROOT"

[[ -s "${OUT}/97B2BR_global_summary.tsv" ]] || {
    echo "ERROR: summary missing" >&2
    exit 20
}

printf \
"COMPLETED\t%s\tjob=%s\n" \
"$(date -Iseconds)" \
"$SLURM_JOB_ID" \
> "${OUT}/97B2BR_COMPLETE.ok"

echo
echo "============================================================"
echo "97B2BR GLOBAL SUMMARY"
echo "============================================================"

column -t -s $'\t' \
    "${OUT}/97B2BR_global_summary.tsv"

echo
echo "============================================================"
echo "97B2BR FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"
