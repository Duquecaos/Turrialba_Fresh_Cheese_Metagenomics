#!/bin/bash

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

export LC_ALL=C
export LANG=C
export LANGUAGE=C

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

OLD="${ROOT}/98_bacteriocin_mining/98C1_positive_screen"

OUT="${ROOT}/98_bacteriocin_mining/98C1R_positive_screen_corrected"

CATALOG="${ROOT}/70_incomplete_gene_catalog/incomplete_lt90_new.prodigal.faa"

EGGNOG="${ROOT}/72_incomplete_eggnog_consolidated/incomplete_lt90_eggnog_annotations_698054.tsv"

MASTER="${ROOT}/38_bacteriocin_master/structural_candidates_protein.faa"

PY="${SCRIPT_DIR}/98c1_integrate_positive_screen.py"


rm -rf "$OUT"
mkdir -p "$OUT"


for F in \
    "$CATALOG" \
    "$EGGNOG" \
    "$MASTER" \
    "$PY" \
    "${OLD}/98C1_GA_hmmscan.tbl" \
    "${OLD}/98C1_noGA_hmmscan.tbl" \
    "${OLD}/98C1_blastp_short_raw.tsv"
do

    [[ -s "$F" ]] || {
        echo "ERROR: missing $F" >&2
        exit 20
    }

done


python3 \
    "$PY" \
    "$CATALOG" \
    "$EGGNOG" \
    "${OLD}/98C1_GA_hmmscan.tbl" \
    "${OLD}/98C1_noGA_hmmscan.tbl" \
    "${OLD}/98C1_blastp_short_raw.tsv" \
    "$MASTER" \
    "$OUT"


cat > "${OUT}/98C1R_methodological_scope.tsv" <<'EOF'
field	value
analysis	corrected_reintegration_of_98C1_existing_raw_results
HMM_reexecuted	NO
BLAST_reexecuted	NO
qcov_definition	abs_qend_minus_qstart_plus_1_divided_by_query_length
qcov_capped_at_100_percent	YES
original_alignment_length_used_as_qcov	NO
new_sequence_search_performed	NO
reads_remapped	NO
amplicon_data_used	NO
next_step	98C2_positive_contig_and_bin_context
EOF


printf \
"COMPLETED\t%s\n" \
"$(date -Iseconds)" \
> "${OUT}/98C1R_COMPLETE.ok"


echo "============================================================"
echo "98C1R GLOBAL SUMMARY"
echo "============================================================"

column -t -s $'\t' \
"${OUT}/98C1_global_summary.tsv"

echo
echo "98C1R FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
