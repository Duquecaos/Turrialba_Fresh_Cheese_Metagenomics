#!/bin/bash

#SBATCH --job-name=95e_16S_synthesis
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=95e_16S_synthesis_%j.out
#SBATCH --error=95e_16S_synthesis_%j.err

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

SCRIPT="${SCRIPT_DIR}/95e_16S_thesis_synthesis.py"

OUT="${ROOT}/97_16S_shotgun/95E_thesis_synthesis"

rm -rf "$OUT"

echo "============================================================"
echo "95E - FINAL SHOTGUN 16S THESIS SYNTHESIS"
echo "Inicio: $(date)"
echo "Job: ${SLURM_JOB_ID}"
echo "============================================================"

python3 \
    "$SCRIPT" \
    "$ROOT"


[[ -s "${OUT}/95E_global_summary.tsv" ]] || {
    echo "ERROR: global summary missing" >&2
    exit 2
}


[[ -s "${OUT}/95E_master_16S_locus_sample_table.tsv" ]] || {
    echo "ERROR: master table missing" >&2
    exit 3
}


[[ -s "${OUT}/figures/95E_Fig1_mapping_support_by_sample.svg" ]] || {
    echo "ERROR: figure 1 missing" >&2
    exit 4
}


[[ -s "${OUT}/figures/95E_Fig2_top15_taxa_signal_heatmap.svg" ]] || {
    echo "ERROR: figure 2 missing" >&2
    exit 5
}


MASTER_ROWS=$(
    tail -n +2 \
        "${OUT}/95E_master_16S_locus_sample_table.tsv" \
    | wc -l
)

if [[ "$MASTER_ROWS" -ne 363 ]]
then
    echo "ERROR: master rows=$MASTER_ROWS expected=363" >&2
    exit 6
fi


MAG_ROWS=$(
    tail -n +2 \
        "${OUT}/95E_table_MAG_16S_summary.tsv" \
    | wc -l
)

if [[ "$MAG_ROWS" -ne 3 ]]
then
    echo "ERROR: MAG summary rows=$MAG_ROWS expected=3" >&2
    exit 7
fi


printf \
"COMPLETED\t%s\tjob=%s\n" \
"$(date -Iseconds)" \
"${SLURM_JOB_ID}" \
> "${OUT}/95E_COMPLETE.ok"


echo
echo "============================================================"
echo "GLOBAL SUMMARY"
echo "============================================================"

column -t -s $'\t' \
    "${OUT}/95E_global_summary.tsv"

echo
echo "============================================================"
echo "95E FINALIZO CORRECTAMENTE"
echo "OBJETIVO SHOTGUN 16S CERRADO"
echo "ES SEGURO SALIR"
echo "============================================================"

