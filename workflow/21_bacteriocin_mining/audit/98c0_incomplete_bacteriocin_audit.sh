#!/bin/bash

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

export LC_ALL=C
export LANG=C
export LANGUAGE=C

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

OUT="${ROOT}/98_bacteriocin_mining/98C0_incomplete_audit"

CATALOG="${ROOT}/70_incomplete_gene_catalog/incomplete_lt90_new.prodigal.faa"

EGGNOG="${ROOT}/72_incomplete_eggnog_consolidated/incomplete_lt90_eggnog_annotations_698054.tsv"

MASTER="${ROOT}/38_bacteriocin_master/structural_candidates_protein.faa"

PY="${SCRIPT_DIR}/98c0_audit_incomplete_catalog.py"

HMMROOT="${ROOT}/databases/antismash_7.1.0_package/antismash/detection/hmm_detection/data"

APPTAINER="${APPTAINER:-/opt/ohpc/pub/libs/apptainer/1.4.1/bin/apptainer}"

HMMER="${HMMER:-/opt/ohpc/pub/containers/BIO/hmmer-3.4.sif}"

DIAMOND="${DIAMOND:-/opt/ohpc/pub/containers/BIO/diamond-2.1.9.sif}"


rm -rf "$OUT"

mkdir -p "$OUT"


echo "============================================================"
echo "98C0 - INCOMPLETE BACTERIOCIN CATALOG AUDIT"
echo "START=$(date -Iseconds)"
echo "============================================================"


# ============================================================
# Required inputs
# ============================================================

for F in \
    "$CATALOG" \
    "$EGGNOG" \
    "$MASTER" \
    "$PY"
do

    if [[ ! -s "$F" ]]
    then

        echo "ERROR: missing $F" >&2
        exit 20

    fi

done


# ============================================================
# Protein catalog + eggNOG audit
# ============================================================

python3 \
    "$PY" \
    "$CATALOG" \
    "$EGGNOG" \
    "$MASTER" \
    "$OUT"


# ============================================================
# HMM inventory
# ============================================================

HMM_TABLE="${OUT}/98C0_bacteriocin_HMM_audit.tsv"

printf \
"hmm\tpresent\tbytes\tNAME\tLENG\thas_GA\thas_TC\thas_NC\n" \
> "$HMM_TABLE"


HMMS=(
    "BacteriocIIc_cy.hmm"
    "Bacteriocin_II.hmm"
    "Bacteriocin_IIc.hmm"
    "Bacteriocin_IId.hmm"
    "Bacteriocin_IIi.hmm"
    "Lactococcin.hmm"
    "Lactococcin_972.hmm"
    "LcnG-beta.hmm"
    "glycocin.hmm"
    "lasso.hmm"
    "thio_amide.hmm"
    "thiostrepton.hmm"
)


for H in "${HMMS[@]}"
do

    F="${HMMROOT}/${H}"

    if [[ -s "$F" ]]
    then

        BYTES=$(
            stat -c %s "$F"
        )

        NAME=$(
            awk '
                $1 == "NAME" {
                    print $2
                    exit
                }
            ' "$F"
        )

        LENG=$(
            awk '
                $1 == "LENG" {
                    print $2
                    exit
                }
            ' "$F"
        )

        if grep -q '^GA[[:space:]]' "$F"
        then
            GA=YES
        else
            GA=NO
        fi

        if grep -q '^TC[[:space:]]' "$F"
        then
            TC=YES
        else
            TC=NO
        fi

        if grep -q '^NC[[:space:]]' "$F"
        then
            NC=YES
        else
            NC=NO
        fi

        printf \
        "%s\tYES\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$H" \
        "$BYTES" \
        "${NAME:-NA}" \
        "${LENG:-NA}" \
        "$GA" \
        "$TC" \
        "$NC" \
        >> "$HMM_TABLE"

    else

        printf \
        "%s\tNO\t0\tNA\tNA\tNA\tNA\tNA\n" \
        "$H" \
        >> "$HMM_TABLE"

    fi

done


# ============================================================
# Tool versions
# ============================================================

TOOL_TABLE="${OUT}/98C0_tool_versions.tsv"

printf \
"tool\tcontainer\tversion\n" \
> "$TOOL_TABLE"


HMM_VERSION=$(
    "$APPTAINER" exec \
    "$HMMER" \
    hmmsearch -h \
    2>&1 \
    | head -n 2 \
    | tail -n 1 \
    || true
)


DIAMOND_VERSION=$(
    "$APPTAINER" exec \
    "$DIAMOND" \
    diamond version \
    2>&1 \
    | head -n 1 \
    || true
)


printf \
"HMMER\t%s\t%s\n" \
"$HMMER" \
"${HMM_VERSION:-UNRESOLVED}" \
>> "$TOOL_TABLE"


printf \
"DIAMOND\t%s\t%s\n" \
"$DIAMOND" \
"${DIAMOND_VERSION:-UNRESOLVED}" \
>> "$TOOL_TABLE"


# ============================================================
# Candidate mapping / metadata files near incomplete workflow
# ============================================================

MAPFILES="${OUT}/98C0_candidate_mapping_files.txt"

find "$ROOT" \
    -maxdepth 2 \
    -type f \
    \( \
        -name '*.tsv' \
        -o \
        -name '*.txt' \
        -o \
        -name '*.csv' \
    \) \
    -print \
| grep -E \
'/(69|70|71|72|73|74|75)_[^/]+/' \
| grep -Ei \
'(map|manifest|metadata|source|catalog|gene|contig|bin|partial|summary|QC)' \
| sort \
> "$MAPFILES" \
|| true


# ============================================================
# Keyword audit in consolidated eggNOG
# Auxiliary evidence only.
# ============================================================

KEYWORDS="${OUT}/98C0_eggnog_bacteriocin_keyword_hits.tsv"

{
    grep -Ei \
    'bacterioc|lactococc|pediocin|nisin|lanthipept|lantibiotic|RiPP|thiopept|lassopept|microviridin|ranthipept|sactipept|glycocin' \
    "$EGGNOG" \
    || true
} > "$KEYWORDS"


N_KEYWORD=$(
    wc -l < "$KEYWORDS"
)


# ============================================================
# Operational summary
# ============================================================

CATALOG_N=$(
    awk -F $'\t' '
        $1 == "catalog_sequences" {
            print $2
        }
    ' "${OUT}/98C0_global_summary.tsv"
)


MATCH_N=$(
    awk -F $'\t' '
        $1 == "eggnog_unique_catalog_ID_matches" {
            print $2
        }
    ' "${OUT}/98C0_global_summary.tsv"
)


HMM_PRESENT=$(
    awk -F $'\t' '
        NR > 1 && $2 == "YES" {
            n++
        }
        END {
            print n+0
        }
    ' "$HMM_TABLE"
)


{
    printf "metric\tvalue\n"
    printf "catalog_sequences\t%s\n" "${CATALOG_N:-NA}"
    printf "eggnog_catalog_ID_matches\t%s\n" "${MATCH_N:-NA}"
    printf "bacteriocin_HMMs_expected\t12\n"
    printf "bacteriocin_HMMs_present\t%s\n" "$HMM_PRESENT"
    printf "eggnog_keyword_hits\t%s\n" "$N_KEYWORD"
    printf "HMM_search_performed\tNO\n"
    printf "DIAMOND_search_performed\tNO\n"
    printf "new_bacteriocin_prediction_performed\tNO\n"
    printf "reads_remapped\tNO\n"
    printf "amplicon_data_used\tNO\n"
    printf "next_step\t98C1_positive_screen_HMM_plus_similarity\n"

} > "${OUT}/98C0_operational_summary.tsv"


# ============================================================
# Scope
# ============================================================

cat > "${OUT}/98C0_methodological_scope.tsv" <<'EOF'
field	value
analysis	audit_before_positive_evidence_search_in_incomplete_bins
catalog_scope	proteins_predicted_from_lt90_percent_complete_bins
incomplete_bin_negative_result_interpretation	NOT_BIOLOGICAL_ABSENCE
HMM_profiles	role_as_profile_screen_not_full_antiSMASH_prediction
DIAMOND	role_as_homology_screen_against_existing_candidate_proteins
eggNOG	role_as_auxiliary_annotation_evidence
generic_transporter_annotation_sufficient_for_bacteriocin_call	NO
new_prediction	NO
expression	NOT_ASSESSED
antimicrobial_activity	NOT_ASSESSED
amplicon_data_used	NO
EOF


# ============================================================
# Validate
# ============================================================

for F in \
    "${OUT}/98C0_global_summary.tsv" \
    "${OUT}/98C0_catalog_length_distribution.tsv" \
    "${OUT}/98C0_catalog_header_examples.txt" \
    "${OUT}/98C0_header_pattern_summary.tsv" \
    "${OUT}/98C0_eggnog_schema.txt" \
    "${OUT}/98C0_bacteriocin_HMM_audit.tsv" \
    "${OUT}/98C0_tool_versions.tsv" \
    "${OUT}/98C0_operational_summary.tsv" \
    "${OUT}/98C0_methodological_scope.tsv"
do

    if [[ ! -s "$F" ]]
    then

        echo "ERROR: missing output $F" >&2
        exit 30

    fi

done


printf \
"COMPLETED\t%s\n" \
"$(date -Iseconds)" \
> "${OUT}/98C0_COMPLETE.ok"


echo
echo "============================================================"
echo "OPERATIONAL SUMMARY"
echo "============================================================"

column -t -s $'\t' \
"${OUT}/98C0_operational_summary.tsv"


echo
echo "============================================================"
echo "98C0 FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"
