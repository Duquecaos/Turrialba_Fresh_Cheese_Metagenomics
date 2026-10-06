#!/bin/bash

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

export LC_ALL=C
export LANG=C
export LANGUAGE=C

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

OUT="${ROOT}/98_bacteriocin_mining/98C2_positive_context"

POSITIVE="${ROOT}/98_bacteriocin_mining/98C1R_positive_screen_corrected/98C1_positive_proteins.tsv"

GENES="${ROOT}/70_incomplete_gene_catalog/gene_to_contig_context.tsv"

MEMBERSHIP="${ROOT}/69_incomplete_unique_contigs/bin_to_unique_contig_membership.tsv"

CONTIGCAT="${ROOT}/69_incomplete_unique_contigs/unique_contig_catalog_lt90.tsv"

EGGNOG="${ROOT}/72_incomplete_eggnog_consolidated/incomplete_lt90_eggnog_annotations_698054.tsv"

PY="${SCRIPT_DIR}/98c2_positive_context.py"


rm -rf "$OUT"

mkdir -p "$OUT"


for F in \
    "$POSITIVE" \
    "$GENES" \
    "$MEMBERSHIP" \
    "$CONTIGCAT" \
    "$EGGNOG" \
    "$PY"
do

    if [[ ! -s "$F" ]]
    then
        echo "ERROR: missing $F" >&2
        exit 20
    fi

done


python3 \
    "$PY" \
    "$POSITIVE" \
    "$GENES" \
    "$MEMBERSHIP" \
    "$CONTIGCAT" \
    "$EGGNOG" \
    "$OUT"


for F in \
    "${OUT}/98C2_positive_context_summary.tsv" \
    "${OUT}/98C2_neighbor_genes.tsv" \
    "${OUT}/98C2_positive_bin_memberships.tsv" \
    "${OUT}/98C2_primary_target_contigs.tsv" \
    "${OUT}/98C2_exploratory_noGA_contigs.tsv" \
    "${OUT}/98C2_multihit_primary_contigs.tsv" \
    "${OUT}/98C2_global_summary.tsv"
do

    [[ -s "$F" ]] || {
        echo "ERROR: missing output $F" >&2
        exit 30
    }

done


cat > "${OUT}/98C2_methodological_scope.tsv" <<'EOF'
field	value
analysis	genomic_and_bin_context_reconstruction_for_98C1_positive_proteins
neighbor_window	plus_minus_5_predicted_genes_on_same_unique_contig
primary_followup	homology_or_curated_GA_HMM_positive
supplementary_noGA_only	exploratory_not_primary_for_specialized_prediction
represented_in_final18	separates_redundant_context_from_additional_incomplete_context
additional_to_final18	equivalent_to_novel_bacteriocin	NO
eggNOG_context_categories	auxiliary_descriptive_flags_only
generic_transporter_alone	sufficient_for_bacteriocin_call	NO
contaminated_bin_context	weakens_bin_or_taxon_attribution_but_not_sequence_level_positive_evidence
negative_in_incomplete_bin	equivalent_to_biological_absence	NO
expression	NOT_ASSESSED
antimicrobial_activity	NOT_ASSESSED
reads_remapped	NO
amplicon_data_used	NO
EOF


printf \
"COMPLETED\t%s\n" \
"$(date -Iseconds)" \
> "${OUT}/98C2_COMPLETE.ok"


echo
echo "============================================================"
echo "98C2 GLOBAL SUMMARY"
echo "============================================================"

column -t -s $'\t' \
"${OUT}/98C2_global_summary.tsv"


echo
echo "============================================================"
echo "98C2 FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"
