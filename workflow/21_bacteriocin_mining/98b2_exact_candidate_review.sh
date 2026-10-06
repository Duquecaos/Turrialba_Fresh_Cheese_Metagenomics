#!/bin/bash

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

export LC_ALL=C
export LANG=C
export LANGUAGE=C

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

OUT="${ROOT}/98_bacteriocin_mining/98B2_exact_candidate_review"

ATTRLOC="${ROOT}/49_bacteriocin_locus_matrices/locus_summary_by_ATTRLOC.tsv"

PARSER="${SCRIPT_DIR}/98b2_exact_candidate_review.py"


ANTI=$(
    find \
    "${ROOT}/98_bacteriocin_mining" \
    -type f \
    -name '98B0_ATTRLOC_antismash_comparison.tsv' \
    -print \
    -quit
)


mkdir -p "$OUT"


echo "============================================================"
echo "98B2 - EXACT CANDIDATE-LEVEL REVIEW"
echo "START=$(date -Iseconds)"
echo "============================================================"


if [[ ! -s "$ATTRLOC" ]]
then
    echo "ERROR: ATTRLOC table missing" >&2
    exit 20
fi


if [[ -z "${ANTI:-}" || ! -s "$ANTI" ]]
then
    echo "ERROR: 98B0 antiSMASH comparison not found" >&2
    exit 21
fi


if [[ ! -s "$PARSER" ]]
then
    echo "ERROR: parser missing" >&2
    exit 22
fi


echo "ATTRLOC=$ATTRLOC"
echo "ANTISMASH=$ANTI"
echo


python3 \
    "$PARSER" \
    "$ROOT" \
    "$ATTRLOC" \
    "$ANTI" \
    "$OUT"


for F in \
    "${OUT}/98B2_candidate_exact_evidence.tsv" \
    "${OUT}/98B2_locus_exact_summary.tsv" \
    "${OUT}/98B2_global_summary.tsv"
do

    if [[ ! -s "$F" ]]
    then
        echo "ERROR: missing output $F" >&2
        exit 30
    fi

done


cat > "${OUT}/98B2_methodological_scope.tsv" <<'EOF'
field	value
analysis	exact_candidate_level_review_of_existing_bacteriocin_RiPP_evidence
new_prediction	NO
antiSMASH_role	same_contig_RiPP_context_plus_candidate_translation_crosscheck
BAGEL_role	restricted_to_run_of_the_actual_MAG_and_candidate_contig
Comparippson_role	candidate_identifier_or_protein_sequence_crosscheck
Master38	role_as_candidate_sequence_registry_not_independent_evidence
sequence_exact_match	same_amino_acid_sequence
sequence_contained_match	one_protein_sequence_contains_the_other_with_minimum_length_20aa
BAGEL_cross_MAG_hits	counted_as_support	NO
generic_ORF_number_hits_counted_without_contig	counted_as_support	NO
same_contig_antismash_region	equivalent_to_exact_candidate_detection	NO
candidate_multimethod_support	equivalent_to_functional_bacteriocin_validation	NO
expression	NOT_ASSESSED
antimicrobial_activity	NOT_ASSESSED
amplicon_data_used	NO
EOF


printf \
"COMPLETED\t%s\n" \
"$(date -Iseconds)" \
> "${OUT}/98B2_COMPLETE.ok"


echo
echo "============================================================"
echo "GLOBAL SUMMARY"
echo "============================================================"

column -t -s $'\t' \
"${OUT}/98B2_global_summary.tsv"


echo
echo "============================================================"
echo "98B2 FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"
