#!/bin/bash

#SBATCH --job-name=95d1_16S_quant
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --array=0-17%6

#SBATCH --mail-type=END,FAIL

#SBATCH --output=95d1_16S_quant_%A_%a.out
#SBATCH --error=95d1_16S_quant_%A_%a.err

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

APPTAINER="${APPTAINER:-/opt/ohpc/pub/libs/apptainer/1.4.1/bin/apptainer}"
SAMTOOLS="${SAMTOOLS:-/opt/ohpc/pub/containers/BIO/samtools-1.21.sif}"

D0="${ROOT}/97_16S_shotgun/95D0_mapping_audit"
BASE="${ROOT}/97_16S_shotgun/95D1_locus_quantification"

MANIFEST="${D0}/95D0_bam_manifest.tsv"

# Historical workflow note:
# The original 95D1 parser file was empty. The array job therefore
# generated the samtools depth files but failed before metric
# validation. Metrics are reconstructed from those native depth files
# by the canonical recovery step 95D1R.

ROW=$(
    sed -n \
        "$((SLURM_ARRAY_TASK_ID + 2))p" \
        "$MANIFEST"
)

# Defensive CRLF sanitation
ROW="${ROW//$'\r'/}"

[[ -n "$ROW" ]] || {
    echo "ERROR: manifest row missing" >&2
    exit 2
}

IFS=$'\t' read -r \
    GROUP \
    SAMPLE \
    BAM \
    BAM_SIZE \
    INDEX_EXISTS \
    HEADER_FILE \
    IDXSTATS_FILE \
    <<< "$ROW"


[[ -s "$BAM" ]] || {
    echo "ERROR: BAM missing: $BAM" >&2
    exit 3
}


BED="${D0}/beds/${GROUP}.16S_targets.bed"

[[ -s "$BED" ]] || {
    echo "ERROR: BED missing: $BED" >&2
    exit 4
}


OUTDIR="${BASE}/results/${SAMPLE}"

mkdir -p "$OUTDIR"


DEPTH_ALL="${OUTDIR}/${SAMPLE}.depth.all.tsv"
DEPTH_Q10="${OUTDIR}/${SAMPLE}.depth.mapq10.tsv"


echo "============================================================"
echo "95D1 - 16S LOCUS QUANTIFICATION"
echo "Sample: ${SAMPLE}"
echo "Coassembly: ${GROUP}"
echo "BAM: ${BAM}"
echo "BED: ${BED}"
echo "============================================================"


# ------------------------------------------------------------
# All primary reported alignments
# ------------------------------------------------------------

"$APPTAINER" exec \
    --bind "${ROOT}:${ROOT}" \
    "$SAMTOOLS" \
    samtools depth \
        -q 0 \
        -Q 0 \
        -b "$BED" \
        "$BAM" \
    > "$DEPTH_ALL"


# ------------------------------------------------------------
# Sensitivity analysis: MAPQ >=10
# ------------------------------------------------------------

"$APPTAINER" exec \
    --bind "${ROOT}:${ROOT}" \
    "$SAMTOOLS" \
    samtools depth \
        -q 0 \
        -Q 10 \
        -b "$BED" \
        "$BAM" \
    > "$DEPTH_Q10"


# ------------------------------------------------------------
# Validate native depth outputs
#
# Historical 95D1 generated these files successfully, then failed
# because the intended parser was an empty file. The downstream
# canonical step 95D1R reconstructs all locus metrics from these
# existing depth files without rerunning samtools depth.
# ------------------------------------------------------------

[[ -f "$DEPTH_ALL" ]] || {
    echo "ERROR: depth ALL output missing" >&2
    exit 5
}

[[ -f "$DEPTH_Q10" ]] || {
    echo "ERROR: depth MAPQ10 output missing" >&2
    exit 6
}


DEPTH_ROWS_ALL=$(
    wc -l < "$DEPTH_ALL"
)

DEPTH_ROWS_Q10=$(
    wc -l < "$DEPTH_Q10"
)


printf \
"COMPLETED_DEPTH_GENERATION\t%s\tjob=%s_%s\tsample=%s\tdepth_rows_all=%s\tdepth_rows_q10=%s\n" \
"$(date -Iseconds)" \
"${SLURM_ARRAY_JOB_ID}" \
"${SLURM_ARRAY_TASK_ID}" \
"${SAMPLE}" \
"${DEPTH_ROWS_ALL}" \
"${DEPTH_ROWS_Q10}" \
> "${OUTDIR}/95D1_DEPTH_COMPLETE.ok"


echo
echo "============================================================"
echo "DEPTH GENERATION SUMMARY"
echo "============================================================"

printf "sample\t%s\n" "$SAMPLE"
printf "depth_all\t%s\n" "$DEPTH_ALL"
printf "depth_all_rows\t%s\n" "$DEPTH_ROWS_ALL"
printf "depth_mapq10\t%s\n" "$DEPTH_Q10"
printf "depth_mapq10_rows\t%s\n" "$DEPTH_ROWS_Q10"

echo
echo "Locus metrics are reconstructed downstream by 95D1R."

echo
echo "============================================================"
echo "95D1 ${SAMPLE} DEPTH GENERATION FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"

