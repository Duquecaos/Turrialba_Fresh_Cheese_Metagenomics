#!/bin/bash

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

export LC_ALL=C
export LANG=C
export LANGUAGE=C

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

OUT="${ROOT}/98_bacteriocin_mining/98B1_existing_evidence_integration"

ATTRLOC="${ROOT}/49_bacteriocin_locus_matrices/locus_summary_by_ATTRLOC.tsv"

PARSER="${SCRIPT_DIR}/98b1_integrate_existing_bacteriocin_evidence.py"

mkdir -p "$OUT"


echo "============================================================"
echo "98B1 - EXISTING BACTERIOCIN EVIDENCE INTEGRATION"
echo "START=$(date -Iseconds)"
echo "============================================================"


for F in \
    "$ATTRLOC" \
    "$PARSER"
do

    if [[ ! -s "$F" ]]
    then
        echo "ERROR: missing $F" >&2
        exit 20
    fi

done


python3 \
    "$PARSER" \
    "$ROOT" \
    "$ATTRLOC" \
    "$OUT"


for F in \
    "${OUT}/98B1_source_inventory.tsv" \
    "${OUT}/98B1_ATTRLOC_evidence_matrix.tsv" \
    "${OUT}/98B1_match_details.tsv" \
    "${OUT}/98B1_global_summary.tsv"
do

    if [[ ! -s "$F" ]]
    then
        echo "ERROR: missing output $F" >&2
        exit 30
    fi

done


cat > "${OUT}/98B1_methodological_scope.tsv" <<'EOF'
field	value
analysis	integration_of_existing_bacteriocin_RiPP_computational_evidence
new_prediction	NO
antismash	existing_7.1.0_results_only
BAGEL	existing_outputs_only
Comparippson	existing_outputs_only
source_hit	contig_or_candidate_string_detected_in_relevant_existing_output
source_hit_not_equivalent_to	functional_validation
independent_method_families	antiSMASH_BAGEL_Comparippson
BAGEL35_and_BAGEL36_counted_as_one_family	YES
diagnostic_and_master_tables_not_counted_as_independent_methods	YES
same_contig_support_not_equivalent_to	exact_candidate_gene_overlap
expression	NOT_ASSESSED
activity	NOT_ASSESSED
amplicon_data_used	NO
EOF


printf \
"COMPLETED\t%s\n" \
"$(date -Iseconds)" \
> "${OUT}/98B1_COMPLETE.ok"


echo
echo "============================================================"
echo "GLOBAL SUMMARY"
echo "============================================================"

column -t -s $'\t' \
"${OUT}/98B1_global_summary.tsv"


echo
echo "============================================================"
echo "98B1 FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"
