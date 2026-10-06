#!/bin/bash

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

export LC_ALL=C
export LANG=C
export LANGUAGE=C

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

OUT="${ROOT}/98_bacteriocin_mining/98C_FINAL_synthesis"

B2="${ROOT}/98_bacteriocin_mining/98B2R_strict_comparippson"

C2="${ROOT}/98_bacteriocin_mining/98C2_positive_context"

A3="${ROOT}/98_bacteriocin_mining/98C3A_targeted_antismash/summary"

BR="${ROOT}/98_bacteriocin_mining/98C3BR2R_BAGEL_salvage"

C3="${ROOT}/98_bacteriocin_mining/98C3C_BAGEL_candidate_level"

PY="${SCRIPT_DIR}/98c_final_synthesis.py"


rm -rf "$OUT"

mkdir -p "$OUT"


INPUTS=(
    "${B2}/98B2R_global_summary.tsv"
    "${B2}/98B2R_locus_strict_summary.tsv"
    "${C2}/98C2_positive_context_summary.tsv"
    "${A3}/98C3A_primary_contig_antismash_summary.tsv"
    "${A3}/98C3A_antismash_regions.tsv"
    "${BR}/98C3BR2R_target_BAGEL_summary.tsv"
    "${C3}/98C3C_target_candidate_review.tsv"
    "${C3}/98C3C_candidate_evidence_lines.tsv"
)


for F in "${INPUTS[@]}"
do

    [[ -s "$F" ]] || {
        echo "ERROR: missing input $F" >&2
        exit 20
    }

done


python3 \
    "$PY" \
    "${B2}/98B2R_global_summary.tsv" \
    "${B2}/98B2R_locus_strict_summary.tsv" \
    "${C2}/98C2_positive_context_summary.tsv" \
    "${A3}/98C3A_primary_contig_antismash_summary.tsv" \
    "${A3}/98C3A_antismash_regions.tsv" \
    "${BR}/98C3BR2R_target_BAGEL_summary.tsv" \
    "${C3}/98C3C_target_candidate_review.tsv" \
    "${C3}/98C3C_candidate_evidence_lines.tsv" \
    "$OUT"


for F in \
    "${OUT}/98C_FINAL_primary_protein_synthesis.tsv" \
    "${OUT}/98C_FINAL_primary_contig_synthesis.tsv" \
    "${OUT}/98C_FINAL_known_candidate_recurrence.tsv" \
    "${OUT}/98C_FINAL_BAGEL_candidate_level.tsv" \
    "${OUT}/98C_FINAL_specialized_followup_negative_contigs.tsv" \
    "${OUT}/98C_FINAL_global_summary.tsv" \
    "${OUT}/98C_FINAL_methodological_scope.tsv"
do

    [[ -s "$F" ]] || {
        echo "ERROR: missing output $F" >&2
        exit 30
    }

done


printf \
"COMPLETED\t%s\n" \
"$(date -Iseconds)" \
> "${OUT}/98C_FINAL_COMPLETE.ok"


echo
echo "============================================================"
echo "98C FINAL GLOBAL SUMMARY"
echo "============================================================"

column -t -s $'\t' \
"${OUT}/98C_FINAL_global_summary.tsv"


echo
echo "============================================================"
echo "98C FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"
