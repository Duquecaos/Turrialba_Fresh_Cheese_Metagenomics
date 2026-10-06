#!/bin/bash

#SBATCH --job-name=97c1_path_check
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:20:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=97c1_path_check_%j.out
#SBATCH --error=97c1_path_check_%j.err

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

export LC_ALL=C
export LANG=C
export LANGUAGE=C

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

SCRIPT="${SCRIPT_DIR}/97c1_direct_kraken_pathogen_evidence.py"

OUT="${ROOT}/98_read_taxonomy/97C1_direct_pathogen_evidence"

rm -rf "$OUT"

python3 \
    "$SCRIPT" \
    "$ROOT"


[[ -s "${OUT}/97C1_target_summary.tsv" ]] || {
    echo "ERROR: target summary missing" >&2
    exit 20
}


[[ -s "${OUT}/97C1_all_target_sample_evidence.tsv" ]] || {
    echo "ERROR: evidence table missing" >&2
    exit 21
}


printf \
"COMPLETED\t%s\tjob=%s\n" \
"$(date -Iseconds)" \
"$SLURM_JOB_ID" \
> "${OUT}/97C1_COMPLETE.ok"


echo
echo "============================================================"
echo "97C1 TARGET SUMMARY"
echo "============================================================"

column -t -s $'\t' \
    "${OUT}/97C1_target_summary.tsv"

echo
echo "============================================================"
echo "97C1 FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"
