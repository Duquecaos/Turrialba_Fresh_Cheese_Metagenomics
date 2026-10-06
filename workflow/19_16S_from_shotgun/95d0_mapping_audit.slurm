#!/bin/bash

#SBATCH --job-name=95d0_16S_mapaudit
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=95d0_16S_mapaudit_%j.out
#SBATCH --error=95d0_16S_mapaudit_%j.err

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

APPTAINER="${APPTAINER:-/opt/ohpc/pub/libs/apptainer/1.4.1/bin/apptainer}"
SAMTOOLS="${SAMTOOLS:-/opt/ohpc/pub/containers/BIO/samtools-1.21.sif}"

PYTHON_SCRIPT="${SCRIPT_DIR}/95d0_mapping_audit.py"

OUT="${ROOT}/97_16S_shotgun/95D0_mapping_audit"

rm -rf "$OUT"

mkdir -p \
    "${OUT}/headers" \
    "${OUT}/idxstats"

echo "============================================================"
echo "95D0 - EXISTING BAM MAPPING AUDIT"
echo "Inicio: $(date)"
echo "Job: ${SLURM_JOB_ID}"
echo "============================================================"

# ------------------------------------------------------------
# Headers + idxstats from all 18 BAMs
# ------------------------------------------------------------

COUNT=0

for BAM in \
    "${ROOT}"/06_mapping/*/*/*.sorted.bam
do

    [[ -s "$BAM" ]] || continue

    SAMPLE="$(
        basename "$BAM" .sorted.bam
    )"

    echo "Auditing ${SAMPLE}"

    "$APPTAINER" exec \
        --bind "${ROOT}:${ROOT}" \
        "$SAMTOOLS" \
        samtools view \
            -H \
            "$BAM" \
        > "${OUT}/headers/${SAMPLE}.header.sam"

    "$APPTAINER" exec \
        --bind "${ROOT}:${ROOT}" \
        "$SAMTOOLS" \
        samtools idxstats \
            "$BAM" \
        > "${OUT}/idxstats/${SAMPLE}.idxstats.tsv"

    COUNT=$((COUNT + 1))

done

if [[ "$COUNT" -ne 18 ]]
then
    echo "ERROR: audited BAMs=$COUNT expected=18" >&2
    exit 2
fi

python3 \
    "$PYTHON_SCRIPT" \
    "$ROOT"

[[ -s "${OUT}/95D0_global_summary.tsv" ]] || {
    echo "ERROR: global summary faltante" >&2
    exit 3
}

MISSING=$(
    awk -F'\t' '
        NR>1 && $12 != 0 {n++}
        END {print n+0}
    ' "${OUT}/95D0_mapping_strategy.tsv"
)

if [[ "$MISSING" -ne 0 ]]
then
    echo "ERROR: BAMs con target contigs ausentes=$MISSING" >&2
    exit 4
fi

printf \
"COMPLETED\t%s\tjob=%s\n" \
"$(date -Iseconds)" \
"${SLURM_JOB_ID}" \
> "${OUT}/95D0_COMPLETE.ok"

echo
echo "============================================================"
echo "GLOBAL SUMMARY"
echo "============================================================"

column -t -s $'\t' \
    "${OUT}/95D0_global_summary.tsv"

echo
echo "============================================================"
echo "95D0 FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"

