#!/bin/bash

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

export LC_ALL=C
export LANG=C
export LANGUAGE=C


ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

D1="${ROOT}/98_bacteriocin_mining/98D1_competitive_mapping_v2/98D1C_final_matrices"

C="${ROOT}/98_bacteriocin_mining/98C_FINAL_synthesis"

ATTR="${ROOT}/49_bacteriocin_locus_matrices/locus_summary_by_ATTRLOC.tsv"

OUT="${ROOT}/98_bacteriocin_mining/98E_FINAL_thesis_synthesis"

PY="${SCRIPT_DIR}/98e_final_bacteriocin_synthesis.py"


rm -rf "$OUT"

mkdir -p "$OUT"


echo "============================================================"
echo "98E - FINAL BACTERIOCIN/RiPP THESIS SYNTHESIS"
echo "NO PREDICTION"
echo "NO READ MAPPING"
echo "START=$(date -Iseconds)"
echo "============================================================"


# ============================================================
# Inputs
# ============================================================

INPUTS=(
    "${D1}/98D1C_reference_summary_30.tsv"
    "${D1}/98D1C_global_summary.tsv"
    "${C}/98C_FINAL_global_summary.tsv"
    "${C}/98C_FINAL_known_candidate_recurrence.tsv"
    "$ATTR"
    "$PY"
)


for F in "${INPUTS[@]}"
do

    [[ -s "$F" ]] || {
        echo "ERROR: missing input $F" >&2
        exit 20
    }

done


# ============================================================
# Integration
# ============================================================

python3 \
    "$PY" \
    "${D1}/98D1C_reference_summary_30.tsv" \
    "${D1}/98D1C_global_summary.tsv" \
    "${C}/98C_FINAL_global_summary.tsv" \
    "${C}/98C_FINAL_known_candidate_recurrence.tsv" \
    "$ATTR" \
    "$OUT"


# ============================================================
# Validate outputs
# ============================================================

OUTPUTS=(
    "${OUT}/98E_FINAL_context_catalog_30.tsv"
    "${OUT}/98E_FINAL_known_candidate_recurrences.tsv"
    "${OUT}/98E_FINAL_existing_ATTRLOC_7.tsv"
    "${OUT}/98E_FINAL_additional_specialized_20.tsv"
    "${OUT}/98E_FINAL_exploratory_3.tsv"
    "${OUT}/98E_FINAL_additional_BAGEL_candidate_level.tsv"
    "${OUT}/98E_FINAL_specialized_without_robust_mapping.tsv"
    "${OUT}/98E_FINAL_key_findings.tsv"
    "${OUT}/98E_FINAL_reporting_guardrails.tsv"
    "${OUT}/98E_FINAL_thesis_ready_summary.txt"
)


for F in "${OUTPUTS[@]}"
do

    [[ -s "$F" ]] || {
        echo "ERROR: missing output $F" >&2
        exit 30
    }

done


NCONTEXT=$(
    awk '
        NR > 1 {n++}
        END {print n+0}
    ' "${OUT}/98E_FINAL_context_catalog_30.tsv"
)


[[ "$NCONTEXT" -eq 30 ]] || {
    echo "ERROR: final context catalog=$NCONTEXT" >&2
    exit 31
}


printf \
"COMPLETED\t%s\n" \
"$(date -Iseconds)" \
> "${OUT}/98E_FINAL_COMPLETE.ok"


echo
echo "============================================================"
echo "98E KEY FINDINGS"
echo "============================================================"

column -t -s $'\t' \
"${OUT}/98E_FINAL_key_findings.tsv"


echo
echo "============================================================"
echo "THESIS-READY SUMMARY"
echo "============================================================"

cat \
"${OUT}/98E_FINAL_thesis_ready_summary.txt"


echo
echo "============================================================"
echo "98E FINALIZO CORRECTAMENTE"
echo "RAMA COMPUTACIONAL DE BACTERIOCINAS/RiPP CERRADA"
echo "ES SEGURO SALIR"
echo "============================================================"
