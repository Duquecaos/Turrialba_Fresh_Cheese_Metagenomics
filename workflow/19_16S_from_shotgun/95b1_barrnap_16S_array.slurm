#!/bin/bash

#SBATCH --job-name=95b1_barrnap16S
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --array=0-23%8

#SBATCH --mail-type=END,FAIL

#SBATCH --output=95b1_barrnap16S_%A_%a.out
#SBATCH --error=95b1_barrnap16S_%A_%a.err

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

APPTAINER="${APPTAINER:-/opt/ohpc/pub/libs/apptainer/1.4.1/bin/apptainer}"
BARRNAP="${BARRNAP:-/opt/ohpc/pub/containers/BIO/barrnap-0.9.sif}"

BASE="${ROOT}/97_16S_shotgun/95B1_barrnap"

MANIFEST="${BASE}/95B1_input_manifest.tsv"

PARSER="${SCRIPT_DIR}/95b1_parse_barrnap_16s.py"

ROW=$(
    sed -n \
    "$((SLURM_ARRAY_TASK_ID + 2))p" \
    "$MANIFEST"
)

# Protección frente a manifiestos CRLF generados por csv.DictWriter.
ROW="${ROW//$'\\r'/}"

if [[ -z "$ROW" ]]; then
    echo "ERROR: no manifest row for ${SLURM_ARRAY_TASK_ID}" >&2
    exit 2
fi

IFS=$'\t' read -r \
    ARRAY_INDEX \
    SOURCE_TYPE \
    SOURCE_ID \
    COASSEMBLY \
    FASTA \
    <<< "$ROW"

TASK_ID="${SOURCE_TYPE}__${SOURCE_ID}"

TASKOUT="${BASE}/results/${TASK_ID}"

mkdir -p "$TASKOUT"

GFF="${TASKOUT}/${TASK_ID}.barrnap.gff"
BERR="${TASKOUT}/${TASK_ID}.barrnap.stderr"

echo "============================================================"
echo "95B1 - BARRNAP 16S"
echo "Array index: ${ARRAY_INDEX}"
echo "Source type: ${SOURCE_TYPE}"
echo "Source ID: ${SOURCE_ID}"
echo "Coassembly: ${COASSEMBLY}"
echo "FASTA: ${FASTA}"
echo "Job: ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "============================================================"

[[ -s "$FASTA" ]] || {
    echo "ERROR: FASTA missing $FASTA" >&2
    exit 3
}

[[ -s "$PARSER" ]] || {
    echo "ERROR: parser missing $PARSER" >&2
    exit 4
}

set +e

"$APPTAINER" exec \
    --bind "${ROOT}:${ROOT}" \
    "$BARRNAP" \
    barrnap \
        --kingdom bac \
        --threads 1 \
        "$FASTA" \
    > "$GFF" \
    2> "$BERR"

RC=$?

set -e

if [[ "$RC" -ne 0 ]]; then
    echo "ERROR: Barrnap exit code=${RC}" >&2
    cat "$BERR" >&2
    exit "$RC"
fi

python3 \
    "$PARSER" \
    "$SOURCE_TYPE" \
    "$SOURCE_ID" \
    "$COASSEMBLY" \
    "$FASTA" \
    "$GFF" \
    "$TASKOUT"

TABLE="${TASKOUT}/${TASK_ID}.16S.tsv"
SEQFA="${TASKOUT}/${TASK_ID}.16S.fasta"
SUMMARY="${TASKOUT}/${TASK_ID}.summary.tsv"

[[ -e "$TABLE" ]] || {
    echo "ERROR: missing 16S table" >&2
    exit 5
}

[[ -e "$SEQFA" ]] || {
    echo "ERROR: missing 16S FASTA" >&2
    exit 6
}

[[ -s "$SUMMARY" ]] || {
    echo "ERROR: missing summary" >&2
    exit 7
}

{
    echo -e "field\tvalue"
    echo -e "source_type\t${SOURCE_TYPE}"
    echo -e "source_id\t${SOURCE_ID}"
    echo -e "coassembly\t${COASSEMBLY}"
    echo -e "input_fasta\t${FASTA}"
    echo -e "barrnap_kingdom\tbac"
    echo -e "barrnap_exit_code\t${RC}"
    echo -e "analysis_target\t16S_rRNA"
    echo -e "length_ge1300_interpretation\tnear_full_length_candidate_not_proof_of_completeness"
    echo -e "edge_touching_interpretation\tpossible_boundary_truncation"
    echo -e "negative_MAG_16S_equals_taxon_absence\tNO"
} > "${TASKOUT}/${TASK_ID}.metadata.tsv"

sha256sum \
    "$GFF" \
    "$TABLE" \
    "$SEQFA" \
> "${TASKOUT}/${TASK_ID}.sha256.tsv"

printf \
"COMPLETED\t%s\tjob=%s_%s\tsource=%s:%s\n" \
"$(date -Iseconds)" \
"${SLURM_ARRAY_JOB_ID}" \
"${SLURM_ARRAY_TASK_ID}" \
"${SOURCE_TYPE}" \
"${SOURCE_ID}" \
> "${TASKOUT}/95B1_COMPLETE.ok"

echo
echo "============================================================"
echo "SUMMARY"
echo "============================================================"

column -t -s $'\t' \
"$SUMMARY"

echo
echo "============================================================"
echo "95B1 ${TASK_ID} FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"

