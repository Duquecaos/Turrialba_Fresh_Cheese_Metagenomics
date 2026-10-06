#!/bin/bash

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

export LC_ALL=C
export LANG=C
export LANGUAGE=C

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

OUT="${ROOT}/98_bacteriocin_mining/98C3C_BAGEL_candidate_level"

BAGEL_SUMMARY="${ROOT}/98_bacteriocin_mining/98C3BR2R_BAGEL_salvage/98C3BR2R_target_BAGEL_summary.tsv"

PROTEINS="${ROOT}/98_bacteriocin_mining/98C1R_positive_screen_corrected/98C1_positive_proteins.faa"

SESSION="${ROOT}/98_bacteriocin_mining/98C3BR2_targeted_BAGEL/session"

PY="${SCRIPT_DIR}/98c3c_bagel_candidate_level.py"


rm -rf "$OUT"

mkdir -p "$OUT"


echo "============================================================"
echo "98C3C - BAGEL CANDIDATE-LEVEL REVIEW"
echo "NO BAGEL REEXECUTION"
echo "START=$(date -Iseconds)"
echo "============================================================"


for F in \
    "$BAGEL_SUMMARY" \
    "$PROTEINS" \
    "$PY"
do

    [[ -s "$F" ]] || {
        echo "ERROR: missing $F" >&2
        exit 20
    }

done


[[ -d "$SESSION" ]] || {
    echo "ERROR: BAGEL session missing: $SESSION" >&2
    exit 21
}


python3 \
    "$PY" \
    "$BAGEL_SUMMARY" \
    "$PROTEINS" \
    "$SESSION" \
    "$OUT"


for F in \
    "${OUT}/98C3C_target_candidate_review.tsv" \
    "${OUT}/98C3C_BAGEL_ORF_sequence_matches.tsv" \
    "${OUT}/98C3C_candidate_evidence_lines.tsv" \
    "${OUT}/98C3C_artifact_inventory.tsv" \
    "${OUT}/98C3C_global_summary.tsv"
do

    [[ -s "$F" ]] || {
        echo "ERROR: output missing or empty: $F" >&2
        exit 30
    }

done


cat > "${OUT}/98C3C_methodological_scope.tsv" <<'EOF'
field	value
analysis	candidate_level_review_of_three_validated_BAGEL_AOI_positive_incomplete_contexts
BAGEL_reexecuted	NO
candidate_sequence_source	98C1R_positive_proteins
BAGEL_input_queryfolder_excluded_from_sequence_evidence	YES
sequence_mapping	exact_or_full_target_containment_in_BAGEL_predicted_ORF
strict_candidate_evidence	BAGEL_bacteriocin_HMM_or_bacteriocin_BLAST_or_predict_artifact
GeneTable_only_candidate_support	NO
AOI_table_only_candidate_support	NO
BAGEL_subcomponents_counted_as_independent_method_families	NO
BAGEL_method_family_count	ONE
locus_level_AOI_support	functional_confirmation_NO
candidate_level_BAGEL_support	functional_confirmation_NO
expression	NOT_ASSESSED
antimicrobial_activity	NOT_ASSESSED
amplicon_data_used	NO
EOF


printf \
"COMPLETED\t%s\n" \
"$(date -Iseconds)" \
> "${OUT}/98C3C_COMPLETE.ok"


echo
echo "============================================================"
echo "98C3C GLOBAL SUMMARY"
echo "============================================================"

column -t -s $'\t' \
"${OUT}/98C3C_global_summary.tsv"


echo
echo "============================================================"
echo "98C3C FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"
