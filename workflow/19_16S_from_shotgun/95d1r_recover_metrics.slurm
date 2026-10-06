#!/bin/bash

#SBATCH --job-name=95d1r_16S_recover
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=95d1r_16S_recover_%j.out
#SBATCH --error=95d1r_16S_recover_%j.err

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

D0="${ROOT}/97_16S_shotgun/95D0_mapping_audit"
D1="${ROOT}/97_16S_shotgun/95D1_locus_quantification"

SCRIPT="${SCRIPT_DIR}/95d1r_recover_metrics.py"

echo "============================================================"
echo "95D1R - RECOVER METRICS FROM EXISTING DEPTH FILES"
echo "Inicio: $(date)"
echo "Job: ${SLURM_JOB_ID}"
echo "============================================================"

# ============================================================
# Validate all 36 existing depth files before parsing
# ============================================================

MISSING=0
EMPTY=0
SAMPLES=0

tail -n +2 \
    "${D0}/95D0_bam_manifest.tsv" \
| while IFS=$'\t' read -r \
    GROUP \
    SAMPLE \
    BAM \
    SIZE \
    INDEX \
    HEADER \
    IDX
do

    SAMPLES=$((SAMPLES + 1))

    for SUFFIX in \
        "depth.all.tsv" \
        "depth.mapq10.tsv"
    do

        F="${D1}/results/${SAMPLE}/${SAMPLE}.${SUFFIX}"

        if [[ ! -e "$F" ]]; then
            echo "MISSING: $F" >&2
            MISSING=$((MISSING + 1))

        elif [[ ! -s "$F" ]]; then
            echo "EMPTY: $F" >&2
            EMPTY=$((EMPTY + 1))
        fi

    done

done

# More reliable independent count outside pipeline subshell:
NDEPTH=$(
    find "${D1}/results" \
        -mindepth 2 \
        -maxdepth 2 \
        -type f \
        \( \
            -name '*.depth.all.tsv' \
            -o \
            -name '*.depth.mapq10.tsv' \
        \) \
    | wc -l
)

echo "depth_files_found=${NDEPTH}"

if [[ "$NDEPTH" -ne 36 ]]
then
    echo "ERROR: expected 36 depth files" >&2
    exit 2
fi


NONEMPTY=$(
    find "${D1}/results" \
        -mindepth 2 \
        -maxdepth 2 \
        -type f \
        \( \
            -name '*.depth.all.tsv' \
            -o \
            -name '*.depth.mapq10.tsv' \
        \) \
        -size +0c \
    | wc -l
)

echo "nonempty_depth_files=${NONEMPTY}"

if [[ "$NONEMPTY" -ne 36 ]]
then
    echo "ERROR: expected 36 nonempty depth files" >&2
    exit 3
fi


python3 \
    "$SCRIPT" \
    "$ROOT"


[[ -s "${D1}/95D1R_global_summary.tsv" ]] || {
    echo "ERROR: recovery summary missing" >&2
    exit 4
}


NMARK=$(
    find "${D1}/results" \
        -mindepth 2 \
        -maxdepth 2 \
        -type f \
        -name '95D1R_RECOVERED.ok' \
    | wc -l
)

if [[ "$NMARK" -ne 18 ]]
then
    echo "ERROR: recovered markers=$NMARK expected=18" >&2
    exit 5
fi


printf \
"COMPLETED\t%s\tjob=%s\n" \
"$(date -Iseconds)" \
"${SLURM_JOB_ID}" \
> "${D1}/95D1R_COMPLETE.ok"


echo
echo "============================================================"
echo "GLOBAL RECOVERY SUMMARY"
echo "============================================================"

column -t -s $'\t' \
    "${D1}/95D1R_global_summary.tsv"

echo
echo "============================================================"
echo "95D1R FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"

