#!/bin/bash

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

IN="${ROOT}/98_bacteriocin_mining/98B2_exact_candidate_review/98B2_candidate_exact_evidence.tsv"

OUT="${ROOT}/98_bacteriocin_mining/98B2R_strict_comparippson"

PY="${SCRIPT_DIR}/98b2r_strict_comparippson.py"

rm -rf "$OUT"
mkdir -p "$OUT"

python3 \
    "$PY" \
    "$ROOT" \
    "$IN" \
    "$OUT"


for F in \
    "${OUT}/98B2R_candidate_strict_evidence.tsv" \
    "${OUT}/98B2R_comparippson_result_lines.tsv" \
    "${OUT}/98B2R_locus_strict_summary.tsv" \
    "${OUT}/98B2R_global_summary.tsv"
do

    [[ -s "$F" ]] || {
        echo "ERROR: missing $F" >&2
        exit 20
    }

done


cat > "${OUT}/98B2R_methodological_scope.tsv" <<'EOF'
field	value
analysis	strict_correction_of_Comparippson_evidence
Comparippson_input_FASTA_counted_as_evidence	NO
Comparippson_candidate_registry_counted_as_evidence	NO
Comparippson_result_tables_used	all_hits_best_hits_top10
Comparippson_result_support	exact_candidate_identifier_present_in_real_result_table
generic_ORF_only_matching	NO
antiSMASH_candidate_support	exact_or_contained_translation_in_RiPP_region
BAGEL_candidate_support	actual_same_MAG_run_plus_candidate_level_match
specialized_method_families	antiSMASH_BAGEL_Comparippson
expression	NOT_ASSESSED
antimicrobial_activity	NOT_ASSESSED
amplicon_data_used	NO
EOF


printf \
"COMPLETED\t%s\n" \
"$(date -Iseconds)" \
> "${OUT}/98B2R_COMPLETE.ok"


echo
echo "============================================================"
column -t -s $'\t' \
"${OUT}/98B2R_global_summary.tsv"
echo "============================================================"
echo "98B2R FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
