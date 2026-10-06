#!/bin/bash

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

export LC_ALL=C
export LANG=C
export LANGUAGE=C


ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

BASE="${ROOT}/98_bacteriocin_mining/98D1_competitive_mapping_v2"

REFDIR="${BASE}/reference"

D0="${ROOT}/98_bacteriocin_mining/98D0_deduplication"

C="${ROOT}/98_bacteriocin_mining/98C_FINAL_synthesis"

OUT="${BASE}/98D1C_final_matrices"

PY="${SCRIPT_DIR}/98d1c_build_final_matrices.py"


rm -rf "$OUT"

mkdir -p "$OUT"


echo "============================================================"
echo "98D1C - FINAL 18x30 BACTERIOCIN/RiPP MATRICES"
echo "NO READ MAPPING"
echo "START=$(date -Iseconds)"
echo "============================================================"


# ============================================================
# Inputs
# ============================================================

for F in \
    "${REFDIR}/98D1AR_competitive_locus_windows_30_metadata.tsv" \
    "${REFDIR}/98D1AR_sample_manifest.tsv" \
    "${D0}/98D0_context_redundancy_groups.tsv" \
    "${C}/98C_FINAL_primary_contig_synthesis.tsv" \
    "$PY"
do

    [[ -s "$F" ]] || {
        echo "ERROR: missing $F" >&2
        exit 20
    }

done


NCOMPLETE=$(
    find "${BASE}/samples" \
        -mindepth 2 \
        -maxdepth 2 \
        -type f \
        -name '98D1B_COMPLETE.ok' \
    | wc -l
)


NMETRICS=$(
    find "${BASE}/samples" \
        -mindepth 2 \
        -maxdepth 2 \
        -type f \
        -name '*.98D1B_locus_metrics.tsv' \
    | wc -l
)


echo "COMPLETE_SAMPLES=${NCOMPLETE}/18"
echo "METRIC_FILES=${NMETRICS}/18"


[[ "$NCOMPLETE" -eq 18 ]] || {
    echo "ERROR: incomplete mapping samples" >&2
    exit 21
}


[[ "$NMETRICS" -eq 18 ]] || {
    echo "ERROR: incomplete metric files" >&2
    exit 22
}


# ============================================================
# Aggregate
# ============================================================

python3 \
    "$PY" \
    "$BASE" \
    "${REFDIR}/98D1AR_competitive_locus_windows_30_metadata.tsv" \
    "${REFDIR}/98D1AR_sample_manifest.tsv" \
    "${D0}/98D0_context_redundancy_groups.tsv" \
    "${C}/98C_FINAL_primary_contig_synthesis.tsv" \
    "$OUT"


# ============================================================
# Validate main products
# ============================================================

OUTPUTS=(
    "${OUT}/98D1C_long_18x30.tsv"
    "${OUT}/98D1C_matrix_detection_class.tsv"
    "${OUT}/98D1C_matrix_robust_binary.tsv"
    "${OUT}/98D1C_matrix_breadth_ALL.tsv"
    "${OUT}/98D1C_matrix_breadth_MAPQ10.tsv"
    "${OUT}/98D1C_matrix_mean_depth_ALL.tsv"
    "${OUT}/98D1C_matrix_mean_depth_MAPQ10.tsv"
    "${OUT}/98D1C_matrix_mean_depth_MAPQ10_per_million_input_pairs.tsv"
    "${OUT}/98D1C_reference_summary_30.tsv"
    "${OUT}/98D1C_existing_ATTRLOC_summary.tsv"
    "${OUT}/98D1C_additional_specialized_summary.tsv"
    "${OUT}/98D1C_additional_exploratory_summary.tsv"
    "${OUT}/98D1C_sample_summary_18.tsv"
    "${OUT}/98D1C_ATTRLOC002_vs_ATTRLOC006.tsv"
    "${OUT}/98D1C_global_summary.tsv"
    "${OUT}/98D1C_methodological_scope.tsv"
)


for F in "${OUTPUTS[@]}"
do

    [[ -s "$F" ]] || {
        echo "ERROR: missing output $F" >&2
        exit 30
    }

done


LONG_ROWS=$(
    awk '
        NR > 1 {n++}
        END {print n+0}
    ' "${OUT}/98D1C_long_18x30.tsv"
)


[[ "$LONG_ROWS" -eq 540 ]] || {
    echo "ERROR: expected 540 long rows; found $LONG_ROWS" >&2
    exit 31
}


printf \
"COMPLETED\t%s\n" \
"$(date -Iseconds)" \
> "${OUT}/98D1C_COMPLETE.ok"


echo
echo "============================================================"
echo "98D1C GLOBAL SUMMARY"
echo "============================================================"

column -t -s $'\t' \
"${OUT}/98D1C_global_summary.tsv"


echo
echo "============================================================"
echo "98D1C FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"
